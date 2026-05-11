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
    vec4 spawn_params; // .x = effective_radius
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

        vec3 chosen_dir = field_direction.xyz;
        bool has_mask = (direction_mask != 0u);

        if (has_mask) {
            // Pick one of the active cardinal directions at random
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

            // +X=0, -X=1, +Y=2, -Y=3, +Z=4, -Z=5
            float sign_val = ((chosen & 1u) == 0u) ? 1.0 : -1.0;
            uint axis = chosen >> 1u;
            chosen_dir = vec3(0.0);
            if (axis == 0u) chosen_dir.x = sign_val;
            else if (axis == 1u) chosen_dir.y = sign_val;
            else chosen_dir.z = sign_val;

            // Spawn on a disc perpendicular to chosen_dir at sim boundary
            vec3 perp1, perp2;
            if (axis == 0u)      { perp1 = vec3(0,1,0); perp2 = vec3(0,0,1); }
            else if (axis == 1u) { perp1 = vec3(1,0,0); perp2 = vec3(0,0,1); }
            else                 { perp1 = vec3(1,0,0); perp2 = vec3(0,1,0); }

            float disc_r = spawn_params.x * 1.1 * sqrt(rand01(seed));
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

        uint flags = floatBitsToUint(vel_f.w);
        flags = (flags & ~5u) | (rand01(seed) > photon_ratio ? 1u : 0u);
        flags |= 2u;
        vel_f.w = uintBitsToFloat(flags);
    }

    data[base] = vec4(pos, pos_e.w);
    data[base + 1u] = vec4(vel, vel_f.w);
}
