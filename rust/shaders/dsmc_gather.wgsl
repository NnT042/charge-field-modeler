// DSMC step 4: reorder photons into cell-sorted order using CellSorter's
// `sorted` permutation (dead photons land in the trailing "dead" bucket, so
// after this pass the live photons are contiguous in [0, n_live)).

struct GpuPhoton {
    a: vec4<f32>,
    b: vec4<f32>,
    c: vec4<u32>,
};

struct GatherParams {
    n: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> sorted: array<u32>;
@group(0) @binding(1) var<storage, read> src: array<GpuPhoton>;
@group(0) @binding(2) var<storage, read_write> dst: array<GpuPhoton>;
@group(0) @binding(3) var<uniform> params: GatherParams;

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let k = linear_wg * 256u + lid.x;
    if (k >= params.n) {
        return;
    }
    dst[k] = src[sorted[k]];
}
