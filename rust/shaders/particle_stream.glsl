// Bench kernel: particles as 2x vec4 (pos.xyz+w, vel.xyz+w), advance and
// reflect off an axis-aligned box. Mirrors the PhotonBuffer layout used by
// godot/shaders/compute/*.glsl.
#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 0) buffer Particles {
    vec4 data[];
};

layout(std140, set = 0, binding = 1) uniform Params {
    float dt;
    float box_min;
    float box_max;
    uint n;
};

void main() {
    // wgpu caps a single dispatch dimension at 65535 workgroups, so large N
    // (tens of millions of particles) dispatch as a 2-D workgroup grid;
    // recover the linear thread index from gl_NumWorkGroups rather than
    // assuming gl_GlobalInvocationID.x alone covers it.
    uint linear_wg = gl_WorkGroupID.x + gl_WorkGroupID.y * gl_NumWorkGroups.x;
    uint idx = linear_wg * 256u + gl_LocalInvocationID.x;
    if (idx >= n) return;

    uint base = idx * 2u;
    vec3 pos = data[base].xyz;
    vec3 vel = data[base + 1u].xyz;

    pos += vel * dt;

    for (int i = 0; i < 3; i++) {
        if (pos[i] < box_min) {
            pos[i] = box_min + (box_min - pos[i]);
            vel[i] = -vel[i];
        } else if (pos[i] > box_max) {
            pos[i] = box_max - (pos[i] - box_max);
            vel[i] = -vel[i];
        }
    }

    data[base] = vec4(pos, data[base].w);
    data[base + 1u] = vec4(vel, data[base + 1u].w);
}
