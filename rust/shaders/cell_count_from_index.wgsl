// DSMC binning pass 1, precomputed-index variant: tally per-cell counts from
// a `cell_of` buffer some other kernel already computed (e.g. dsmc_move.wgsl,
// which needs the body/lane/boundary logic to decide a photon's cell anyway).
// Companion to cell_count.wgsl's position-based pass; CellSorter picks
// whichever the caller asked for.

struct Params {
    n: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read> cell_of_in: array<u32>;
@group(0) @binding(1) var<storage, read_write> counts: array<atomic<u32>>;
@group(0) @binding(2) var<storage, read_write> cell_of_out: array<u32>;
@group(0) @binding(3) var<uniform> params: Params;

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let idx = linear_wg * 256u + lid.x;
    if (idx >= params.n) {
        return;
    }
    let cell = cell_of_in[idx];
    cell_of_out[idx] = cell;
    atomicAdd(&counts[cell], 1u);
}
