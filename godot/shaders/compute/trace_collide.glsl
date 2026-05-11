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

    if (segment_count > 0u) {
        float best_dist = 1e30;
        vec3 best_normal = vec3(0.0);
        vec3 best_surf_vel = vec3(0.0);
        float best_radius = 1.0;
        vec3 best_closest = vec3(0.0);

        for (uint s = 0u; s < segment_count; s++) {
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
            }
        }

        if (best_dist < best_radius) {
            vec3 v_rel = vel - best_surf_vel;
            float vn = dot(v_rel, best_normal);
            if (vn < 0.0) {
                v_rel = v_rel - 2.0 * vn * best_normal;
                vel = v_rel + best_surf_vel;
                vel = normalize(vel);
                pos = best_closest + best_normal * (best_radius + 0.01);
                impulse = vec4(-2.0 * vn * best_normal * energy, 1.0);
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
