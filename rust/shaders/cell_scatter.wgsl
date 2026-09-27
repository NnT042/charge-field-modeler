// DSMC binning pass 3: scatter particle indices into sorted order using a
// per-cell atomic cursor initialised to the exclusive-scan offsets.
// WGSL for the same reason as cell_count.wgsl (naga's GLSL frontend has no
// atomic support).

struct Params {
    n: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> cell_of: array<u32>;
@group(0) @binding(1) var<storage, read_write> cursor: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> sorted: array<u32>;
@group(0) @binding(3) var<uniform> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let idx = linear_wg * 256u + lid.x;
    if (idx >= params.n) {
        return;
    }

    let cell = cell_of[idx];
    let slot = atomicAdd(&cursor[cell], 1u);
    sorted[slot] = idx;
}
