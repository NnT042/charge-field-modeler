#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 0) buffer XBuffer {
    float x[];
};

layout(std430, set = 0, binding = 1) buffer YBuffer {
    float y[];
};

layout(std140, set = 0, binding = 2) uniform Params {
    float a;
    uint n;
    uint _pad0;
    uint _pad1;
};

void main() {
    uint gid = gl_GlobalInvocationID.x;
    if (gid >= n) return;
    y[gid] = a * x[gid] + y[gid];
}
