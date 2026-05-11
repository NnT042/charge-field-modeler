#[compute]
#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 0) buffer ImpulseBuffer {
    vec4 impulses[];
};

layout(std430, set = 0, binding = 1) buffer ReductionOutput {
    vec4 partial_sums[];
};

layout(std140, set = 0, binding = 2) uniform ReductionParams {
    uint element_count;
    uint _pad0;
    uint _pad1;
    uint _pad2;
};

shared vec4 sdata[256];

void main() {
    uint tid = gl_LocalInvocationID.x;
    uint gid = gl_GlobalInvocationID.x;

    sdata[tid] = (gid < element_count) ? impulses[gid] : vec4(0.0);
    barrier();

    for (uint s = 128u; s > 0u; s >>= 1u) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        barrier();
    }

    if (tid == 0u) {
        partial_sums[gl_WorkGroupID.x] = sdata[0];
    }
}
