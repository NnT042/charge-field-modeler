#[compute]
#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 0) buffer PhotonBuffer {
    vec4 data[];
};

layout(std140, set = 0, binding = 1) uniform SimParams {
    vec4 sim_center_radius;
    float dt;
    float time_scale;
    float photon_ratio;
    float direction_strength;
    vec4 field_direction;
    uint photon_count;
    uint frame_number;
    uint segment_count;
    uint direction_mask;
    vec4 spawn_params;
    vec4 spin_axis_emit;   // .xyz = spin axis, .w = emit_enabled
    vec4 focus_pos_r;      // .xyz = focus pos, .w = collision_radius
    vec4 focus_vel_rate;   // .xyz = focus vel, .w = emit_rate
};

uint hash(uint x) {
    x ^= x >> 16u;
    x *= 0x45d9f3bu;
    x ^= x >> 16u;
    x *= 0x45d9f3bu;
    x ^= x >> 16u;
    return x;
}

float rand01(inout uint seed) {
    seed = hash(seed);
    return float(seed) / 4294967295.0;
}

vec3 random_direction(inout uint seed) {
    float theta = rand01(seed) * 6.2831853;
    float z = 2.0 * rand01(seed) - 1.0;
    float r = sqrt(1.0 - z * z);
    return vec3(r * cos(theta), r * sin(theta), z);
}

void main() {
    uint idx = gl_GlobalInvocationID.x;
    if (idx >= photon_count) return;

    uint base = idx * 2u;
    vec4 pos_e = data[base];
    vec4 vel_f = data[base + 1u];

    vec3 pos = pos_e.xyz;
    vec3 vel = vel_f.xyz;

    pos += vel * dt * time_scale;

    vec3 offset = pos - sim_center_radius.xyz;
    float dist = length(offset);

    if (dist > sim_center_radius.w) {
        uint seed = hash(idx ^ (frame_number * 1999u));
        vec3 rand_vel = random_direction(seed);

        // Pole emission: spawn charge at ghost sphere poles, traveling inward
        bool emitted = false;
        if (focus_vel_rate.w > 0.0 && rand01(seed) < focus_vel_rate.w) {
            vec3 axis = spin_axis_emit.xyz;
            float axis_len = length(axis);
            if (axis_len > 1e-6) {
                axis = axis / axis_len;
                float ghost_r = spawn_params.x;

                // Pick north or south pole
                float pole_sign = (rand01(seed) > 0.5) ? 1.0 : -1.0;

                // Cone spread around pole (~30° half-angle = polar cap)
                vec3 perp1, perp2;
                if (abs(axis.y) < 0.9) {
                    perp1 = normalize(cross(axis, vec3(0.0, 1.0, 0.0)));
                } else {
                    perp1 = normalize(cross(axis, vec3(1.0, 0.0, 0.0)));
                }
                perp2 = cross(axis, perp1);

                float spread = sqrt(rand01(seed)) * 0.52;  // ~30° cone, sqrt for uniform area
                float phi = rand01(seed) * 6.2831853;
                float sa = sin(spread);
                vec3 spawn_dir = normalize(
                    axis * pole_sign * cos(spread)
                    + (perp1 * cos(phi) + perp2 * sin(phi)) * sa
                );

                pos = focus_pos_r.xyz + spawn_dir * ghost_r;
                vel = normalize(-axis * pole_sign + random_direction(seed) * 0.08);
                emitted = true;
            }
        }

        if (!emitted) {
            vec3 chosen_dir = field_direction.xyz;
            bool has_mask = (direction_mask != 0u);

            if (has_mask) {
                int count = 0;
                for (uint b = 0u; b < 6u; b++) {
                    if ((direction_mask & (1u << b)) != 0u) count++;
                }
                int pick = int(rand01(seed) * float(count));
                pick = min(pick, count - 1);

                int seen = 0;
                uint chosen = 0u;
                for (uint b = 0u; b < 6u; b++) {
                    if ((direction_mask & (1u << b)) != 0u) {
                        if (seen == pick) { chosen = b; break; }
                        seen++;
                    }
                }

                float sign_val = ((chosen & 1u) == 0u) ? 1.0 : -1.0;
                uint axis = chosen >> 1u;
                chosen_dir = vec3(0.0);
                if (axis == 0u) chosen_dir.x = sign_val;
                else if (axis == 1u) chosen_dir.y = sign_val;
                else chosen_dir.z = sign_val;

                vec3 perp1, perp2;
                if (axis == 0u)      { perp1 = vec3(0,1,0); perp2 = vec3(0,0,1); }
                else if (axis == 1u) { perp1 = vec3(1,0,0); perp2 = vec3(0,0,1); }
                else                 { perp1 = vec3(1,0,0); perp2 = vec3(0,1,0); }

                float disc_r = spawn_params.x * 0.75 * sqrt(rand01(seed));
                float angle = rand01(seed) * 6.2831853;
                vec3 disc_offset = (cos(angle) * perp1 + sin(angle) * perp2) * disc_r;
                float R = sim_center_radius.w;
                float axis_dist = sqrt(max(R * R - disc_r * disc_r, 1.0));
                pos = sim_center_radius.xyz + chosen_dir * axis_dist + disc_offset;
            } else {
                vec3 spawn_dir = random_direction(seed);
                pos = sim_center_radius.xyz + spawn_dir * sim_center_radius.w;
            }

            vec3 inward = normalize(sim_center_radius.xyz - pos);
            vec3 travel_dir = has_mask ? -chosen_dir : chosen_dir;
            vec3 base_dir = mix(inward, travel_dir, direction_strength);
            vel = normalize(mix(base_dir, rand_vel, field_direction.w));

            if (spawn_params.y > 0.0) {
                pos += vel * rand01(seed) * 2.0 * sim_center_radius.w;
            }
        }

        uint flags = floatBitsToUint(vel_f.w);
        flags = (flags & ~5u) | (rand01(seed) > photon_ratio ? 1u : 0u);
        flags |= 2u;
        vel_f.w = uintBitsToFloat(flags);
    }

    data[base] = vec4(pos, pos_e.w);
    data[base + 1u] = vec4(vel, vel_f.w);
}
