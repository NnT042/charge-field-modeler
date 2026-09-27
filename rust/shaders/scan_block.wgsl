// One level of a multi-level exclusive prefix sum (Hillis-Steele scan per
// 256-element block). Writes the exclusive-scanned values back in place and
// emits each block's total sum to `block_sums` for the caller to scan
// recursively and add back with add_block_sums.wgsl.
//
// This is WGSL rather than GLSL: naga's GLSL frontend always lowers the
// `barrier()` builtin to `Barrier::all()`, and naga's SPIR-V backend then
// picks *Subgroup* as the execution scope whenever SUB_GROUP is one of the
// set flags (see naga's `write_control_barrier`) — so a GLSL `barrier()`
// only synchronizes within a subgroup (32 threads on this NVIDIA hardware),
// not the full 256-thread workgroup this scan needs. WGSL's
// `workgroupBarrier()` lowers to a real WORK_GROUP-only barrier and has no
// such downgrade.

struct Params {
    n: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
};

@group(0) @binding(0) var<storage, read_write> data: array<u32>;
@group(0) @binding(1) var<storage, read_write> block_sums: array<u32>;
@group(0) @binding(2) var<uniform> params: Params;

var<workgroup> temp: array<u32, 256>;

@compute @workgroup_size(256)
fn main(@builtin(local_invocation_id) lid3: vec3<u32>, @builtin(workgroup_id) wid3: vec3<u32>, @builtin(num_workgroups) nwg3: vec3<u32>) {
    let tid = lid3.x;
    // wgpu caps a single dispatch dimension at 65535 workgroups; recover the
    // linear block index from a possibly-2-D workgroup grid rather than
    // relying on global_invocation_id.x, which only reflects workgroup_id.x.
    let block_id = wid3.x + wid3.y * nwg3.x;
    let gid = block_id * 256u + tid;

    let val = select(0u, data[gid], gid < params.n);
    temp[tid] = val;
    workgroupBarrier();

    var offset = 1u;
    loop {
        if (offset >= 256u) {
            break;
        }
        let src = i32(tid) - i32(offset);
        let t = select(0u, temp[u32(src)], src >= 0);
        workgroupBarrier();
        temp[tid] = temp[tid] + t;
        workgroupBarrier();
        offset = offset << 1u;
    }

    let inclusive = temp[tid];
    let exclusive = inclusive - val;

    if (gid < params.n) {
        data[gid] = exclusive;
    }
    if (tid == 255u) {
        block_sums[block_id] = inclusive;
    }
}
