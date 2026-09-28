// DSMC step 5: one WORKGROUP per cell (B2 rewrite; dsmc_common.wgsl is
// concatenated ahead of this source by engine.rs).
//
// B1 limitation this replaces: the old kernel paired only the first
// MAX_CELL=64 photons (by id) while computing the NTC candidate count from
// the FULL cell population, so any cell over 64 was biased -- exactly where
// the densification physics (doc §1) lives. Fix (doc §11b, work item 1):
// load up to MAX_CELL id-ordered (id, slot) pairs into workgroup shared
// memory, bitonic-sort them in parallel, then have thread 0 run the
// sequential NTC candidate loop over the full sorted list. Cells whose live
// count exceeds MAX_CELL still only pair among the first MAX_CELL (a correct
// fallback is prohibitively slow per doc's own analysis), but now increment
// TALLY_CELL_OVERFLOW so the host can print a loud warning every step instead
// of silently mismodeling dense cells.
//
// Determinism: unchanged in spirit -- pairing depends only on the photon ids
// present in a cell, never on CellSorter's atomic-scatter order, because the
// shared-memory sort re-establishes id order before any pairing happens.

const MAX_CELL: u32 = 2048u;
const SORT_LOCAL: u32 = 256u;
const SORT_STRIDE: u32 = MAX_CELL / SORT_LOCAL; // 8

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
    domain_center_x: f32, domain_center_y: f32, domain_center_z: f32, exchange_rule: u32,
    density_accumulate: u32, _pad9a: u32, _pad9b: u32, _pad9c: u32,
};

const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const FLAG_COLLIDED_SINCE_CROSS: u32 = 1u << 13u;
const TALLY_SCALE: f32 = 4096.0;
// Tally indices (must stay in sync with gpu_types.rs).
const TALLY_COLLISIONS_IDX: u32 = 61u;
const TALLY_MOM_DRIFT_COUNT_IDX: u32 = 62u;
const TALLY_MOM_DRIFT_SUM_IDX: u32 = 67u;
const TALLY_CELL_OVERFLOW_IDX: u32 = 68u;
const TALLY_EXCHANGE_DRIFT_SUM_IDX: u32 = 69u;
const TALLY_REDIRECT_DRIFT_SUM_IDX: u32 = 70u;
const TALLY_LOW_RANK64_IDX: u32 = 71u;
// The pre-rewrite kernel's fixed candidate cap -- used only to measure how
// concentrated collisions are among the lowest-id photons per cell (see
// TALLY_LOW_RANK64_IDX doc in gpu_types.rs), not as a limit here.
const RANK64: u32 = 64u;

// NOTE: `photons` here is the POST-GATHER buffer -- dsmc_gather.wgsl already
// wrote `dst[k] = src[sorted[k]]`, so position `off+k` (from `offsets`/
// `counts`, which describe that same sorted layout) already IS the photon;
// there is no `sorted` permutation left to apply here.
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

var<workgroup> sh_ids: array<u32, MAX_CELL>;
var<workgroup> sh_slots: array<u32, MAX_CELL>;

fn chi_of(flags: u32) -> f32 {
    return select(-1.0, 1.0, (flags & FLAG_SPECIES_PLUS) != 0u);
}

@compute @workgroup_size(SORT_LOCAL)
fn main(
    @builtin(workgroup_id) wid: vec3<u32>,
    @builtin(local_invocation_id) lid: vec3<u32>,
    @builtin(num_workgroups) nwg: vec3<u32>,
) {
    let cell = wid.x + wid.y * nwg.x;
    if (cell >= params.num_cells) {
        return;
    }
    let count = counts[cell];
    // sigma<=0 and count<2 are workgroup-uniform (both depend only on `cell`
    // / global params), so this early-out ahead of any barrier is legal.
    if (count < 2u || params.sigma <= 0.0) {
        return;
    }
    let off = offsets[cell];
    let capn = min(count, MAX_CELL);
    // Sort only the next power of two >= capn, not all MAX_CELL slots: cells
    // typically hold a handful of photons, and sorting 2048 padded entries in
    // every one of ~900k cells made collide ~100x slower than move. `len` is
    // workgroup-uniform (depends only on `cell`), so the barrier loops below
    // stay uniform. Bitonic order over a padded prefix is identical, so
    // results are bit-for-bit unchanged.
    var len = 2u;
    loop {
        if (len >= capn) { break; }
        len = len << 1u;
    }
    let strides = (len + SORT_LOCAL - 1u) / SORT_LOCAL;
    if (count > MAX_CELL && lid.x == 0u) {
        atomicAdd(&tallies[TALLY_CELL_OVERFLOW_IDX], 1);
    }

    // --- Load (id, slot) pairs into shared memory, pad the tail with the
    // max sentinel so the bitonic network (which needs a power-of-two
    // length) sorts real entries to the front. ---
    for (var m = 0u; m < strides; m = m + 1u) {
        let i = lid.x + m * SORT_LOCAL;
        if (i >= len) { continue; }
        if (i < capn) {
            sh_ids[i] = bitcast<u32>(photons[off + i].a.w);
            sh_slots[i] = off + i;
        } else {
            sh_ids[i] = 0xFFFFFFFFu;
            sh_slots[i] = 0xFFFFFFFFu;
        }
    }
    workgroupBarrier();

    // --- Parallel bitonic sort by id over MAX_CELL elements. Shared memory
    // is workgroup-visible, so a thread may touch indices it doesn't "own"
    // (i, i^j) as long as every (i, ixj) pair with i<ixj is only handled
    // once -- guaranteed by the `ixj > i` guard below, independent of which
    // thread's loop iteration reaches it. ---
    var k = 2u;
    loop {
        if (k > len) { break; }
        var j = k >> 1u;
        loop {
            if (j == 0u) { break; }
            for (var m = 0u; m < strides; m = m + 1u) {
                let i = lid.x + m * SORT_LOCAL;
                let ixj = i ^ j;
                if (i < len && ixj > i) {
                    let ascending = (i & k) == 0u;
                    let vi = sh_ids[i];
                    let vx = sh_ids[ixj];
                    if ((vi > vx) == ascending) {
                        sh_ids[i] = vx;
                        sh_ids[ixj] = vi;
                        let si = sh_slots[i];
                        sh_slots[i] = sh_slots[ixj];
                        sh_slots[ixj] = si;
                    }
                }
            }
            workgroupBarrier();
            j = j >> 1u;
        }
        k = k << 1u;
    }
    // sh_ids/sh_slots[0..capn) are now id-ordered; [capn..MAX_CELL) hold the
    // sentinel and are never touched below.

    if (lid.x != 0u) {
        return;
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

    // Hash-mixed so distinct (cell, step) pairs never share an RNG stream.
    let key_base = hash_u32(cell ^ hash_u32(params.step ^ 0x9E3779B9u));
    for (var c = 0u; c < m; c = c + 1u) {
        let draw0 = c * 8u;
        let ui = rand01(params.seed_hi, params.seed_lo, key_base, draw0);
        // j drawn uniformly from the n-1 slots != i (doc's fix -- the old
        // "bump j to j+1 when j==i" double-weighted i+1).
        let uj = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 1u);
        let i = min(u32(ui * f32(capn)), capn - 1u);
        let n_minus_1 = max(capn, 2u) - 1u;
        let j_off = min(u32(uj * f32(n_minus_1)), n_minus_1 - 1u);
        let j = (i + 1u + j_off) % capn;

        let slot_a = sh_slots[i];
        let slot_b = sh_slots[j];
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

        // --- D6 exchange: SwapRenorm (legacy, |v|=c only) or RotateRelative
        // (conserves |v|=c AND d1+d2 exactly; doc §4 D6). ---
        var new_da: vec3<f32>;
        var new_db: vec3<f32>;
        var renorm_drift = 0.0;
        if (params.exchange_rule == 1u) {
            // RotateRelative. V = mean, g = relative; |V|>0 case rotates g
            // about V-hat (always perpendicular for equal-speed vectors);
            // the exact head-on case (|V|~=0) rotates g about a random axis
            // perpendicular to g instead.
            let v_vec = (da + db) * 0.5;
            let g_vec = da - db;
            let v_len = length(v_vec);
            // Glancing-biased scattering angle: same u-skew as
            // sample_contact_normal (arp.pdf), mapped through the
            // hard-sphere impact-parameter relation b=sin(chi/2)... doc's
            // exact mapping: b = sqrt(u), chi = pi - 2*asin(b).
            let u0 = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 3u);
            let one_minus_u0 = 1.0 - u0;
            let u = 1.0 - one_minus_u0 * one_minus_u0 * one_minus_u0;
            let b_impact = sqrt(clamp(u, 0.0, 1.0));
            let chi_ang = PI - 2.0 * asin(clamp(b_impact, 0.0, 1.0));
            let s = select(-1.0, 1.0, rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 4u) < 0.5);

            var g_rot: vec3<f32>;
            if (v_len > 1e-6) {
                let vhat = v_vec / v_len;
                let theta = s * chi_ang;
                g_rot = g_vec * cos(theta) + cross(vhat, g_vec) * sin(theta);
            } else {
                let ghat = g_vec / max(length(g_vec), 1e-6);
                let helper = perp_helper(ghat);
                let t = normalize(cross(ghat, helper));
                let bperp = cross(ghat, t);
                let psi = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 5u) * TAU;
                let raxis = t * cos(psi) + bperp * sin(psi);
                g_rot = g_vec * cos(chi_ang) + cross(raxis, g_vec) * sin(chi_ang);
            }
            let raw_da = v_vec + g_rot * 0.5;
            let raw_db = v_vec - g_rot * 0.5;
            new_da = normalize(raw_da); // defensive re-normalise: f32 error only
            new_db = normalize(raw_db);
            renorm_drift = length(raw_da - new_da) + length(raw_db - new_db);
        } else {
            // SwapRenorm (legacy exchange rule, Phase-1-like): sample a
            // grazing-biased contact normal, swap its component between the
            // pair, renormalise both to c. Does not conserve d1+d2 exactly
            // (measured ~0.30 photon-momenta/collision in G1) -- reported
            // for comparison, not used for physics runs.
            let rel_hat = rel / rel_len;
            let nhat = sample_contact_normal(rel_hat, params.seed_hi, params.seed_lo, key_base, draw0 + 3u);
            let raw_da = da - dot(da, nhat) * nhat + dot(db, nhat) * nhat;
            let raw_db = db - dot(db, nhat) * nhat + dot(da, nhat) * nhat;
            let len_a = max(length(raw_da), 1e-6);
            let len_b = max(length(raw_db), 1e-6);
            new_da = raw_da / len_a;
            new_db = raw_db / len_b;
            renorm_drift = length(raw_da - new_da) + length(raw_db - new_db);
        }

        let exchange_drift = length((new_da + new_db) - (da + db));
        atomicAdd(&tallies[TALLY_MOM_DRIFT_SUM_IDX], i32(renorm_drift * TALLY_SCALE));
        atomicAdd(&tallies[TALLY_EXCHANGE_DRIFT_SUM_IDX], i32(exchange_drift * TALLY_SCALE));
        atomicAdd(&tallies[TALLY_MOM_DRIFT_COUNT_IDX], 1); // TALLY_MOM_DRIFT_COUNT (shared denominator)
        atomicAdd(&tallies[TALLY_COLLISIONS_IDX], 1);
        if (max(i, j) < RANK64) {
            atomicAdd(&tallies[TALLY_LOW_RANK64_IDX], 1);
        }

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

        let pre_redirect_sum = new_da + new_db;
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
        let redirect_drift = length((new_da + new_db) - pre_redirect_sum);
        atomicAdd(&tallies[TALLY_REDIRECT_DRIFT_SUM_IDX], i32(redirect_drift * TALLY_SCALE));

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
