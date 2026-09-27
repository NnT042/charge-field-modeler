// DSMC step 5: one thread per grid cell. NTC (no-time-counter) candidate
// selection + pairwise exchange, per doc §3.1/§7. dsmc_common.wgsl is
// concatenated ahead of this source by engine.rs.
//
// Determinism: cells hold "tens" of photons (doc §5), so each thread copies
// its cell's slot list into a small local array and insertion-sorts it by
// photon id before pairing -- this is what makes results independent of the
// GPU's atomic-scatter order in CellSorter (doc §7 step 5).
//
// Simplification vs the full doc: cells with more than MAX_CELL live
// photons only pair among the first MAX_CELL (by id); this is a documented
// B1 limitation, not expected to bite at gate densities.

const MAX_CELL: u32 = 64u;

struct GpuPhoton {
    a: vec4<f32>,
    b: vec4<f32>,
    c: vec4<u32>,
};

struct GpuBody {
    center: vec4<f32>,
    axis: vec4<f32>,
    emit: vec4<f32>,
    lane: vec4<f32>,
};

struct Params {
    n_bodies: u32, n_live: u32, capacity: u32, step: u32,
    dt: f32, sigma: f32, n_ambient: f32, photon_fraction: f32,
    pole_intake_bias: f32, cancel_redirect: f32, spin_transfer: f32, b_max: f32,
    boundary_mode: u32, boundary_extent: f32, rate_kernel: u32, mesh_mode: u32,
    bounce_cap: u32, d5_far_zone_gate: u32, grid_dim_x: u32, grid_dim_y: u32,
    grid_dim_z: u32, cell_h: f32, grid_min_x: f32, grid_min_y: f32,
    grid_min_z: f32, seed_hi: u32, seed_lo: u32, num_cells: u32,
    domain_center_x: f32, domain_center_y: f32, domain_center_z: f32, _pad2: f32,
};

const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const FLAG_COLLIDED_SINCE_CROSS: u32 = 1u << 13u;
const TALLY_SCALE: f32 = 4096.0;

// NOTE: `photons` here is the POST-GATHER buffer -- dsmc_gather.wgsl already
// wrote `dst[k] = src[sorted[k]]`, so position `off+k` (from `offsets`/
// `counts`, which describe that same sorted layout) already IS the photon;
// there is no `sorted` permutation left to apply here. (An earlier version
// of this kernel re-indexed through `sorted` a second time, which silently
// scrambled pairings and showed up as non-reproducible collision counts
// across identical-seed runs -- see g1_determinism.)
@group(0) @binding(0) var<storage, read_write> photons: array<GpuPhoton>;
@group(0) @binding(1) var<storage, read> offsets: array<u32>;
@group(0) @binding(2) var<storage, read> counts: array<u32>;
@group(0) @binding(3) var<storage, read> bodies: array<GpuBody>;
@group(0) @binding(4) var<storage, read> nearest_body: array<u32>;
@group(0) @binding(5) var<storage, read> within_bmax: array<u32>;
@group(0) @binding(6) var<storage, read> free_volume: array<f32>;
@group(0) @binding(7) var<storage, read_write> cell_remainder: array<f32>;
@group(0) @binding(8) var<storage, read_write> tallies: array<atomic<i32>>;
@group(0) @binding(9) var<uniform> params: Params;

fn chi_of(flags: u32) -> f32 {
    return select(-1.0, 1.0, (flags & FLAG_SPECIES_PLUS) != 0u);
}

@compute @workgroup_size(64)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let cell = linear_wg * 64u + lid.x;
    // The grid's own cells only -- the trailing "dead" bucket (index
    // grid_num_cells) is never dispatched over (see engine.rs).
    if (cell >= params.num_cells) {
        return;
    }
    let count = counts[cell];
    if (count < 2u || params.sigma <= 0.0) {
        return;
    }
    let off = offsets[cell];
    let capn = min(count, MAX_CELL);

    var idxs: array<u32, MAX_CELL>;
    for (var k = 0u; k < capn; k = k + 1u) {
        idxs[k] = off + k;
    }
    // Insertion sort by photon id (capn is small: tens of entries).
    for (var a = 1u; a < capn; a = a + 1u) {
        let key = idxs[a];
        let key_id = bitcast<u32>(photons[key].a.w);
        var b = a;
        loop {
            if (b == 0u) {
                break;
            }
            let prev = idxs[b - 1u];
            if (bitcast<u32>(photons[prev].a.w) <= key_id) {
                break;
            }
            idxs[b] = idxs[b - 1u];
            b = b - 1u;
        }
        idxs[b] = key;
    }

    var axis = vec3<f32>(0.0, 0.0, 1.0);
    var body_center = vec3<f32>(0.0, 0.0, 0.0);
    if (params.n_bodies > 0u) {
        let nb = bodies[nearest_body[cell]];
        axis = nb.axis.xyz;
        body_center = nb.center.xyz;
    }
    let within = within_bmax[cell] != 0u;
    let vcell = max(free_volume[cell], 1e-6);
    let g_max: f32 = select(1.0, 2.0, params.rate_kernel == 1u);

    let nf = f32(count);
    var m_f = 0.5 * nf * (nf - 1.0) * params.sigma * g_max * params.dt / vcell + cell_remainder[cell];
    var m_i = i32(floor(m_f));
    m_f = m_f - f32(m_i);
    cell_remainder[cell] = m_f;
    let cap_m = i32(capn) * i32(capn);
    if (m_i > cap_m) {
        m_i = cap_m;
    }
    if (m_i <= 0) {
        return;
    }
    let m = u32(m_i);

    // Hash-mixed so distinct (cell, step) pairs never share an RNG stream
    // (a linear combination aliases: 7919*dc == -104729*ds mod 2^32).
    let key_base = hash_u32(cell ^ hash_u32(params.step ^ 0x9E3779B9u));
    for (var c = 0u; c < m; c = c + 1u) {
        let draw0 = c * 8u;
        let ui = rand01(params.seed_hi, params.seed_lo, key_base, draw0);
        let uj0 = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 1u);
        let i = min(u32(ui * f32(capn)), capn - 1u);
        var j = min(u32(uj0 * f32(capn)), capn - 1u);
        if (j == i) {
            j = (j + 1u) % capn;
        }

        let slot_a = idxs[i];
        let slot_b = idxs[j];
        var pa = photons[slot_a];
        var pb = photons[slot_b];
        let da = pa.b.xyz;
        let db = pb.b.xyz;
        let g = length(da - db);
        let accept_p: f32 = select(1.0, g / g_max, params.rate_kernel == 1u);
        if (rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 2u) >= accept_p) {
            continue;
        }

        let a_emitted = (pa.c[0] & FLAG_KIND_EMITTED) != 0u;
        let b_emitted = (pb.c[0] & FLAG_KIND_EMITTED) != 0u;
        if (params.d5_far_zone_gate != 0u && !a_emitted && !b_emitted && !within) {
            continue;
        }

        let rel = da - db;
        let rel_len = length(rel);
        if (rel_len < 1e-6) {
            continue;
        }
        let rel_hat = rel / rel_len;
        let nhat = sample_contact_normal(rel_hat, params.seed_hi, params.seed_lo, key_base, draw0 + 3u);

        let new_da_raw = da - dot(da, nhat) * nhat + dot(db, nhat) * nhat;
        let new_db_raw = db - dot(db, nhat) * nhat + dot(da, nhat) * nhat;
        let len_a = max(length(new_da_raw), 1e-6);
        let len_b = max(length(new_db_raw), 1e-6);
        var new_da = new_da_raw / len_a;
        var new_db = new_db_raw / len_b;

        let drift = length(new_da_raw - new_da) + length(new_db_raw - new_db);
        atomicAdd(&tallies[67], i32(drift * TALLY_SCALE)); // TALLY_MOM_DRIFT_SUM
        atomicAdd(&tallies[62], 1); // TALLY_MOM_DRIFT_COUNT
        atomicAdd(&tallies[61], 1); // TALLY_COLLISIONS

        let chi_a = chi_of(pa.c[0]);
        let chi_b = chi_of(pb.c[0]);
        var mesh = chi_a * chi_b;
        var dir_sign = 1.0;
        if (params.mesh_mode == 1u) {
            dir_sign = sign(dot(da, db));
            if (dir_sign == 0.0) {
                dir_sign = 1.0;
            }
            mesh = mesh * dir_sign;
        }
        let spin_delta_a = params.spin_transfer * chi_b * dir_sign;
        let spin_delta_b = params.spin_transfer * chi_a * dir_sign;

        if (params.cancel_redirect > 0.0 && params.n_bodies > 0u) {
            let pos_a = pa.a.xyz - body_center;
            let pos_b = pb.a.xyz - body_center;
            if (mesh < 0.0) {
                let pole_a = sign(dot(pos_a, axis)) * axis;
                let blended_a = mix(new_da, pole_a, params.cancel_redirect);
                if (dot(blended_a, blended_a) > 1e-10) {
                    new_da = normalize(blended_a);
                }
                let pole_b = sign(dot(pos_b, axis)) * axis;
                let blended_b = mix(new_db, pole_b, params.cancel_redirect);
                if (dot(blended_b, blended_b) > 1e-10) {
                    new_db = normalize(blended_b);
                }
            } else if (mesh > 0.0) {
                let flat_a = new_da - dot(new_da, axis) * axis;
                if (dot(flat_a, flat_a) > 1e-9) {
                    let blended_a = mix(new_da, normalize(flat_a), params.cancel_redirect);
                    if (dot(blended_a, blended_a) > 1e-10) {
                        new_da = normalize(blended_a);
                    }
                }
                let flat_b = new_db - dot(new_db, axis) * axis;
                if (dot(flat_b, flat_b) > 1e-9) {
                    let blended_b = mix(new_db, normalize(flat_b), params.cancel_redirect);
                    if (dot(blended_b, blended_b) > 1e-10) {
                        new_db = normalize(blended_b);
                    }
                }
            }
        }

        let bounces_a = pa.c[1];
        let bounces_b = pb.c[1];
        let cap = params.bounce_cap;
        if (cap == 0u || bounces_a < cap) {
            pa.b = vec4<f32>(new_da, clamp(pa.b.w + spin_delta_a, -2.0, 2.0));
            pa.c[1] = bounces_a + 1u;
            pa.c[0] = pa.c[0] | FLAG_COLLIDED_SINCE_CROSS;
        }
        if (cap == 0u || bounces_b < cap) {
            pb.b = vec4<f32>(new_db, clamp(pb.b.w + spin_delta_b, -2.0, 2.0));
            pb.c[1] = bounces_b + 1u;
            pb.c[0] = pb.c[0] | FLAG_COLLIDED_SINCE_CROSS;
        }
        photons[slot_a] = pa;
        photons[slot_b] = pb;
    }
}
