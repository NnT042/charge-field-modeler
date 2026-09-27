// Second half of the multi-level scan: adds each block's (already scanned)
// offset back onto every element of that block. No barrier needed here
// (every invocation only touches its own element), but kept as WGSL to
// match scan_block.wgsl's buffer layout.

struct Params {
    n: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read_write> data: array<u32>;
@group(0) @binding(1) var<storage, read> block_sums: array<u32>;
@group(0) @binding(2) var<uniform> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    // Same 2-D-dispatch-safe linear index as scan_block.wgsl.
    let linear_wg = wid.x + wid.y * nwg.x;
    let gid = linear_wg * 256u + lid.x;
    if (gid >= params.n) {
        return;
    }
    let block_id = gid / 256u;
    data[gid] = data[gid] + block_sums[block_id];
}
