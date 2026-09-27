// DSMC binning pass 1: compute each particle's flat cell index in a
// grid_dim^3 grid of unit cells and atomically tally per-cell counts.
//
// This is WGSL rather than GLSL because naga's GLSL frontend (as of naga
// 30.0.1 / wgpu 30.0.1) does not implement atomic built-ins at all
// (atomicAdd et al. come back "Unknown function"), while its WGSL frontend
// supports them fully.

struct Params {
    n: u32,
    grid_dim: u32,
    _pad0: u32,
    _pad1: u32,
};

@group(0) @binding(0) var<storage, read> positions: array<vec4<f32>>;
@group(0) @binding(1) var<storage, read_write> counts: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> cell_of: array<u32>;
@group(0) @binding(3) var<uniform> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    // wgpu caps a single dispatch dimension at 65535 workgroups, so large N
    // dispatches as a 2-D workgroup grid; recover the linear thread index
    // from num_workgroups rather than assuming a 1-D grid.
    let linear_wg = wid.x + wid.y * nwg.x;
    let idx = linear_wg * 256u + lid.x;
    if (idx >= params.n) {
        return;
    }

    let p = positions[idx].xyz;
    let gd = i32(params.grid_dim);
    let cx = clamp(i32(floor(p.x)), 0, gd - 1);
    let cy = clamp(i32(floor(p.y)), 0, gd - 1);
    let cz = clamp(i32(floor(p.z)), 0, gd - 1);
    let cell = u32(cx) + u32(cy) * params.grid_dim + u32(cz) * params.grid_dim * params.grid_dim;

    cell_of[idx] = cell;
    atomicAdd(&counts[cell], 1u);
}
