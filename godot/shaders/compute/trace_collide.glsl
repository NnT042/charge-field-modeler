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

layout(std140, set = 0, binding = 2) uniform FocusState {
    vec4 base_solid;
};

layout(std430, set = 0, binding = 3) buffer ImpulseBuffer {
    vec4 impulses[];
};

layout(std430, set = 0, binding = 4) readonly buffer TraceSegments {
    vec4 segments[];
};

void main() {
    uint idx = gl_GlobalInvocationID.x;
    if (idx >= photon_count) return;

    uint base = idx * 2u;
    vec3 pos = data[base].xyz;
    vec3 vel = data[base + 1u].xyz;
    uint flags = floatBitsToUint(data[base + 1u].w);
    float energy = data[base].w;
    vec4 impulse = vec4(0.0);
    bool hit = false;

    bool zone_mode = spin_axis_emit.w > 0.5;

    if (segment_count > 0u) {
        uint last_hit_seg = 0xFFFFFFFFu;

        for (uint bounce = 0u; bounce < 3u; bounce++) {
            float best_dist = 1e30;
            vec3 best_normal = vec3(0.0);
            vec3 best_surf_vel = vec3(0.0);
            float best_radius = 1.0;
            vec3 best_closest = vec3(0.0);
            uint best_seg = 0xFFFFFFFFu;

            for (uint s = 0u; s < segment_count; s++) {
                if (s == last_hit_seg) continue;

                uint si = s * 4u;
                vec3 a_pos = segments[si].xyz;
                float cap_r = segments[si].w;
                vec3 b_pos = segments[si + 1u].xyz;
                vec3 a_vel = segments[si + 2u].xyz;
                vec3 b_vel = segments[si + 3u].xyz;

                vec3 ab = b_pos - a_pos;
                float ab_len2 = dot(ab, ab);
                float t = 0.0;
                if (ab_len2 > 1e-12) {
                    t = clamp(dot(pos - a_pos, ab) / ab_len2, 0.0, 1.0);
                }
                vec3 closest = a_pos + t * ab;
                vec3 diff = pos - closest;
                float dist = length(diff);

                if (dist < cap_r && dist < best_dist) {
                    best_dist = dist;
                    best_closest = closest;
                    best_radius = cap_r;
                    best_normal = (dist > 1e-6) ? diff / dist : vec3(0.0, 1.0, 0.0);
                    best_surf_vel = mix(a_vel, b_vel, t);
                    best_seg = s;
                }
            }

            if (best_dist >= best_radius) break;

            vec3 v_rel = vel - best_surf_vel;
            float vn = dot(v_rel, best_normal);
            if (vn >= 0.0) break;

            last_hit_seg = best_seg;

            if (zone_mode) {
                vec3 axis = spin_axis_emit.xyz;
                vec3 to_hit = pos - focus_pos_r.xyz;
                float th_len = length(to_hit);
                vec3 hit_dir = (th_len > 1e-6) ? to_hit / th_len : best_normal;
                float cos_axis = abs(dot(hit_dir, axis));
                float equatorial = 1.0 - cos_axis * cos_axis;

                if (equatorial < 0.25) {
                    // Polar: channel through without reflecting
                    pos = best_closest + best_normal * (best_radius + 0.01);
                    vec3 to_center = focus_pos_r.xyz - pos;
                    float tc_len = length(to_center);
                    if (tc_len > 1e-6) {
                        vel = normalize(vel + (to_center / tc_len) * 0.3);
                    }
                } else {
                    // Equatorial: tangent deflection
                    v_rel = v_rel - 2.0 * vn * best_normal;
                    vel = v_rel + best_surf_vel;

                    float surf_speed = length(best_surf_vel);
                    if (surf_speed > 1e-6) {
                        vec3 tangent = best_surf_vel / surf_speed;
                        float chi_sign = ((flags & 1u) == 0u) ? 1.0 : -1.0;
                        vel += chi_sign * tangent * equatorial * 0.4;
                    }

                    vec3 radial_out = pos - focus_pos_r.xyz;
                    float ro_len = length(radial_out);
                    if (ro_len > 1e-6) {
                        vel += (radial_out / ro_len) * equatorial * 0.25;
                    }

                    vel = normalize(vel);
                    pos = best_closest + best_normal * (best_radius + 0.01);
                    impulse += vec4(-2.0 * vn * best_normal * energy, 1.0);
                    flags |= 4u;
                    hit = true;
                }

            } else {
                // Legacy specular reflection
                v_rel = v_rel - 2.0 * vn * best_normal;
                vel = v_rel + best_surf_vel;

                float chi_str = spawn_params.z;
                if (chi_str > 0.0) {
                    float surf_speed = length(best_surf_vel);
                    if (surf_speed > 1e-6) {
                        vec3 tangent = best_surf_vel / surf_speed;
                        float chi_sign = ((flags & 1u) == 0u) ? 1.0 : -1.0;
                        vel += chi_sign * tangent * chi_str;
                    }
                }

                vel = normalize(vel);
                pos = best_closest + best_normal * (best_radius + 0.01);
                impulse += vec4(-2.0 * vn * best_normal * energy, 1.0);
                flags |= 4u;
                hit = true;
            }
        }
    }

    // Always check base_solid sphere (not just when segments absent)
    if (!hit) {
        vec3 center = base_solid.xyz;
        float radius = base_solid.w;
        vec3 to_photon = pos - center;
        float dist = length(to_photon);

        if (dist < radius) {
            vec3 normal = normalize(to_photon);
            float vn = dot(vel, normal);
            if (vn < 0.0) {
                vel = vel - 2.0 * vn * normal;
                pos = center + normal * (radius + 0.01);
                impulse = vec4(-2.0 * vn * normal * energy, 1.0);
                flags |= 4u;
            }
        }
    }

    data[base] = vec4(pos, energy);
    data[base + 1u] = vec4(vel, uintBitsToFloat(flags));
    impulses[idx] = impulse;
}
