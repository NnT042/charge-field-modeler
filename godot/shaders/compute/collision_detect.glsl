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
    vec4 spin_axis_emit;
    vec4 focus_pos_r;
    vec4 focus_vel_rate;
};

layout(std140, set = 0, binding = 2) uniform FocusState {
    vec4 base_solid;
};

layout(std430, set = 0, binding = 3) buffer ImpulseBuffer {
    vec4 impulses[];
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

    data[base] = vec4(pos, energy);
    data[base + 1u] = vec4(vel, uintBitsToFloat(flags));
    impulses[idx] = impulse;
}
