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
    d: vec4<f32>,
    e: vec4<f32>,
    f: vec4<f32>,
    g: vec4<f32>,
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
    density_accumulate: u32, emit_stiffness_level: u32, stiffness_shed: u32, sigma_size_scaling: u32,
    debug_force_equal_mass_exchange: u32, g_max_catch: f32, rough_beta: f32, rough_beta_twist: f32,
    rough_kappa: f32, twist_loss_to_normal: u32, speed_mode: u32, spin_edge0: f32,
    glancing_bias: u32, stacking: u32, stack_max_level: u32, unstack_edge: f32,
    body_model: u32, tau_body: f32, grain_edge: f32, body_twist_loss: u32,
    body_beta: f32, body_beta_twist: f32, _pad14a: u32, maze_mode: u32,
};

// Doc §15.47: bit 13 of maze_mode = local (nearest-neighbour) pairing. The
// cell's candidate count stays NTC's, but each candidate pairs a photon with
// a near neighbour in a Morton order of 4x4x4 sub-cells instead of anyone in
// the cell, so a dense directed flow (disc, vortex wall, beam tube) shadows
// what is behind it at sub-cell scale (user, 2026-11-17).
fn local_pairing() -> bool {
    return ((params.maze_mode >> 13u) & 1u) != 0u;
}

fn spread2(v: u32) -> u32 {
    // 2-bit value -> bits 0 and 3 (for a 3-axis Morton interleave).
    return (v & 1u) | ((v & 2u) << 2u);
}

fn subcell_key(cell: u32, x: vec3<f32>) -> u32 {
    let dx = params.grid_dim_x;
    let dy = params.grid_dim_y;
    let cx = cell % dx;
    let cy = (cell / dx) % dy;
    let cz = cell / (dx * dy);
    let lo = vec3<f32>(params.grid_min_x, params.grid_min_y, params.grid_min_z) + vec3<f32>(f32(cx), f32(cy), f32(cz)) * params.cell_h;
    let s = clamp(floor((x - lo) / params.cell_h * 4.0), vec3<f32>(0.0), vec3<f32>(3.0));
    return spread2(u32(s.x)) | (spread2(u32(s.y)) << 1u) | (spread2(u32(s.z)) << 2u);
}

const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const FLAG_COLLIDED_SINCE_CROSS: u32 = 1u << 13u;
// Doc §15: a photon-photon collision happened during the current body
// passage (cleared by dsmc_move on entering a body envelope).
const FLAG_PP_PASSAGE: u32 = 1u << 14u;
const TALLY_SCALE: f32 = 4096.0;
// Tally indices (must stay in sync with gpu_types.rs).
const TALLY_COLLISIONS_IDX: u32 = 61u;
const TALLY_MOM_DRIFT_COUNT_IDX: u32 = 62u;
const TALLY_MOM_DRIFT_SUM_IDX: u32 = 67u;
const TALLY_CELL_OVERFLOW_IDX: u32 = 68u;
// Doc §15.41: collisions by the two photons' emitter tags (0 never
// captured, 1 + last emitting body), [min tag*5 + max tag], len 25; a
// collision between two different bodies' emission marks both photons.
const TALLY_XCOLL_IDX: u32 = 2964u;
const FLAG_XHIT: u32 = 1u << 21u;
const TALLY_PAR_CELLS_IDX: u32 = 2315u; // doc §15.13: cells run in parallel rounds
const TALLY_EXCHANGE_DRIFT_SUM_IDX: u32 = 69u;
const TALLY_REDIRECT_DRIFT_SUM_IDX: u32 = 70u;
const TALLY_LOW_RANK64_IDX: u32 = 71u;
// E1 (doc §11d) -- must stay in sync with gpu_types.rs.
const TALLY_DEFLECT_SUM_IDX: u32 = 194u;
const TALLY_DEFLECT_SUMSQ_IDX: u32 = 202u;
const TALLY_DEFLECT_COUNT_IDX: u32 = 210u;
// Doc §13 -- must stay in sync with gpu_types.rs.
const TALLY_ACCEPT_OVERFLOW_IDX: u32 = 218u;
const TALLY_PIN_DEFICIT_COUNT_IDX: u32 = 219u;
const TALLY_PIN_DEFICIT_SUM_IDX: u32 = 220u;
const TALLY_EDGE_OVER_C_IDX: u32 = 221u;
const TALLY_TWIST_TO_NORMAL_IDX: u32 = 228u;
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
// Doc §15.13: this cell's NTC candidate count, broadcast to the workgroup.
var<workgroup> sh_m: u32;
// Cells with more candidates than this run them in parallel rounds of
// disjoint pairs; below it, thread 0 runs the original serial loop
// (bit-identical to before).
const PAR_MIN_CANDIDATES: u32 = 8192u;

fn chi_of(flags: u32) -> f32 {
    return select(-1.0, 1.0, (flags & FLAG_SPECIES_PLUS) != 0u);
}



// Doc §13: one rough-sphere contact. Mirrors crate::rough_sphere::collide_full
// (plus contact_normal) line for line; rs_gpu_matches_cpu_reference compares
// the two, so keep them in sync. m = R = 2^k, I = kappa*m*R^2.
fn rough_collide(pa: ptr<function, GpuPhoton>, pb: ptr<function, GpuPhoton>, key_base: u32, draw0: u32) {
    // Stacking: unit linear mass, contact radius of the outer level. Without
    // stacking: E1 stiffness, m = R = 2^k.
    let ra = exp2(f32((*pa).c[3]));
    let rb = exp2(f32((*pb).c[3]));
    let ma = select(ra, 1.0, params.stacking != 0u);
    let mb = select(rb, 1.0, params.stacking != 0u);
    let kappa = params.rough_kappa;
    let ia = kappa * ma * ra * ra;
    let ib = kappa * mb * rb * rb;
    var va = (*pa).b.xyz * (*pa).d.w;
    var vb = (*pb).b.xyz * (*pb).d.w;
    var wa = (*pa).d.xyz;
    var wb = (*pb).d.xyz;
    let wa_pre = wa;
    let wb_pre = wb;
    let e_rot_in = 0.5 * ia * dot(wa, wa) + 0.5 * ib * dot(wb, wb);
    let u = va - vb;
    let ulen = length(u);
    if (ulen < 1e-6) {
        return;
    }
    // Contact normal: geometric hard-sphere impact parameter, b/d = sqrt(rho)
    // uniform over the collision disk, azimuth psi around u.
    let uh = u / ulen;
    let helper = perp_helper(uh);
    let t = normalize(cross(uh, helper));
    let bp = cross(uh, t);
    var rho = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 3u);
    if (params.glancing_bias != 0u) {
        // Same skew as RotateRelative / sample_contact_normal (arp.pdf).
        let om = 1.0 - rho;
        rho = 1.0 - om * om * om;
    }
    let psi = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 4u) * TAU;
    let e_perp = t * cos(psi) + bp * sin(psi);
    let n = normalize(uh * sqrt(max(1.0 - rho, 0.0)) + e_perp * sqrt(rho));

    let g = u + ra * cross(wa, n) + rb * cross(wb, n);
    let gn = dot(g, n);
    if (gn <= 0.0) {
        return;
    }
    let g_t = g - gn * n;
    let mu = ma * mb / (ma + mb);
    let j_n = -mu * 2.0 * gn * n; // e = 1
    let k_t = (1.0 / ma) * (1.0 + 1.0 / kappa) + (1.0 / mb) * (1.0 + 1.0 / kappa);
    let j_t = -(1.0 + params.rough_beta) * g_t / k_t;
    let j = j_n + j_t;
    va = va + j / ma;
    vb = vb - j / mb;
    let nj = cross(n, j);
    // Spin kicks tracked separately so the stacking axis (the kick
    // direction) never comes from subtracting two nearly equal f32 spins.
    var kick_a = ra * nj / ia;
    var kick_b = rb * nj / ib;
    wa = wa + kick_a;
    wb = wb + kick_b;

    // Torsional (drilling) grip about n.
    let omega_n = dot(wa - wb, n);
    let tq = -(1.0 + params.rough_beta_twist) * omega_n / (1.0 / ia + 1.0 / ib);
    let e_rot_before = 0.5 * ia * dot(wa, wa) + 0.5 * ib * dot(wb, wb);
    wa = wa + tq * n / ia;
    wb = wb - tq * n / ib;
    kick_a = kick_a + tq * n / ia;
    kick_b = kick_b - tq * n / ib;
    let lost = e_rot_before - (0.5 * ia * dot(wa, wa) + 0.5 * ib * dot(wb, wb));
    if (params.twist_loss_to_normal != 0u && lost > 0.0) {
        let w_n = dot(va - vb, n);
        let jx = mu * (-w_n - sqrt(w_n * w_n + 2.0 * lost / mu));
        va = va + jx * n / ma;
        vb = vb - jx * n / mb;
        atomicAdd(&tallies[TALLY_TWIST_TO_NORMAL_IDX], i32(lost * TALLY_SCALE));
    }

    // D7 speed handling (PinnedC): both speeds back to c, both spins scaled
    // by one factor so the pair's spin energy is what it came in with.
    // BodyFree (§15.2): each photon back to its own incoming speed instead.
    // Mirrors rough_sphere::pin_pair_to.
    if (params.speed_mode != 0u) {
        if (params.speed_mode == 1u) {
            va = normalize(va);
            vb = normalize(vb);
        } else {
            va = normalize(va) * (*pa).d.w;
            vb = normalize(vb) * (*pb).d.w;
        }
        let e_rot = 0.5 * ia * dot(wa, wa) + 0.5 * ib * dot(wb, wb);
        if (e_rot > 1e-12) {
            let sc = sqrt(e_rot_in / e_rot);
            wa = wa * sc;
            wb = wb * sc;
            // w' - w_pre = sc*(w_pre + kick) - w_pre = (sc - 1)*w_pre + sc*kick
            kick_a = (sc - 1.0) * wa_pre + sc * kick_a;
            kick_b = (sc - 1.0) * wb_pre + sc * kick_b;
        } else {
            let fl = length(nj);
            if (fl < 1e-9) {
                if (e_rot_in > 1e-12) {
                    atomicAdd(&tallies[TALLY_PIN_DEFICIT_COUNT_IDX], 1);
                    atomicAdd(&tallies[TALLY_PIN_DEFICIT_SUM_IDX], i32(e_rot_in * TALLY_SCALE));
                }
            } else {
                let ax = nj / fl;
                wa = ax * sqrt(e_rot_in / ia);
                wb = ax * sqrt(e_rot_in / ib);
                kick_a = wa - wa_pre;
                kick_b = wb - wb_pre;
            }
        }
    }
    let sa = length(va);
    let sb = length(vb);
    let da = va / max(sa, 1e-12);
    let db = vb / max(sb, 1e-12);
    // Edge > c BEFORE stacking resolves it (the stacking trigger rate).
    if (length(wa) * ra > 1.0 + 1e-4) {
        atomicAdd(&tallies[TALLY_EDGE_OVER_C_IDX], 1);
    }
    if (length(wb) * rb > 1.0 + 1e-4) {
        atomicAdd(&tallies[TALLY_EDGE_OVER_C_IDX], 1);
    }
    if (params.stacking != 0u) {
        wa = restack(pa, wa, kick_a);
        wb = restack(pb, wb, kick_b);
    }
    // Scalar spin slot = helicity of the axial level, +-1 at birth.
    (*pa).b = vec4<f32>(da, axial_helicity(pa, wa, da));
    (*pb).b = vec4<f32>(db, axial_helicity(pb, wb, db));
    (*pa).d = vec4<f32>(wa, sa);
    (*pb).d = vec4<f32>(wb, sb);
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
            let pid = bitcast<u32>(photons[off + i].a.w);
            sh_ids[i] = pid;
            if (local_pairing()) {
                // Sort by sub-cell first (Morton), then id; the id is masked
                // below 0x3FFFFFF so no real key equals the sentinel.
                sh_ids[i] = (subcell_key(cell, photons[off + i].a.xyz) << 26u) | (pid & 0x3FFFFFEu);
            }
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

    var axis = vec3<f32>(0.0, 0.0, 1.0);
    var body_center = vec3<f32>(0.0, 0.0, 0.0);
    if (params.n_bodies > 0u) {
        let nb = bodies[nearest_body[cell]];
        axis = nb.axis.xyz;
        body_center = nb.center.xyz;
    }
    let within = within_bmax[cell] != 0u;
    let vcell = max(free_volume[cell], 1e-6);
    let g_max: f32 = select(1.0, params.g_max_catch, params.rate_kernel == 1u);

    // E1 size-driven cross-section (doc §11d, arp.pdf): sigma_pair =
    // sigma*((m1+m2)/2)^2, photon radius proportional to stiffness mass m.
    // NTC must draw candidates using the WORST-CASE (largest) pair
    // cross-section so the accept/reject ratio below stays <= 1; k can only
    // ever be as large as emit_stiffness_level (Transfer shedding moves a
    // level DOWN from its current holder, never invents a higher one), so
    // that config value is a safe, cheap ceiling -- see DsmcConfig doc.
    // Off -> size_max = 1 -> sigma_eff = sigma (bit-identical to HEAD).
    let k_max_ntc = min(params.emit_stiffness_level, 7u);
    let size_max: f32 = select(1.0, exp2(f32(k_max_ntc)), params.sigma_size_scaling != 0u);
    let sigma_eff = params.sigma * size_max * size_max;

    if (lid.x == 0u) {
        let nf = f32(count);
        var m_f = 0.5 * nf * (nf - 1.0) * sigma_eff * g_max * params.dt / vcell + cell_remainder[cell];
        var m_i = i32(floor(m_f));
        m_f = m_f - f32(m_i);
        cell_remainder[cell] = m_f;
        let cap_m = i32(capn) * i32(capn);
        if (m_i > cap_m) {
            m_i = cap_m;
        }
        sh_m = u32(max(m_i, 0));
    }
    let m = workgroupUniformLoad(&sh_m);
    if (m == 0u) {
        return;
    }

    // Hash-mixed so distinct (cell, step) pairs never share an RNG stream.
    let key_base = hash_u32(cell ^ hash_u32(params.step ^ 0x9E3779B9u));

    // Doc §15.13: dense cells (near a strong emitter) would serialize the
    // whole step on one thread. Run their candidates in rounds instead: each
    // round pairs position i with i ^ mask (a random perfect matching of the
    // power-of-two padded list, so every unordered pair is equally likely
    // across rounds), each valid pair is a candidate with probability `sel`,
    // and the rounds and `sel` are set so the expected candidate count per
    // pair equals the serial NTC's m / (n(n-1)/2).
    if (m > PAR_MIN_CANDIDATES) {
        let pairs_total = f32(capn) * f32(capn - 1u) * 0.5;
        let lm1 = f32(len - 1u);
        var rounds = u32(ceil(f32(m) * lm1 / pairs_total));
        var sel = f32(m) * lm1 / (f32(rounds) * pairs_total);
        // §15.47 local: masks 1..7 pair within Morton blocks of 8; nearly
        // every i < capn has a valid partner, ~capn/2 pairs per round.
        let lp = local_pairing();
        if (lp) {
            let per_round = max(f32(capn / 2u), 1.0);
            rounds = u32(ceil(f32(m) / per_round));
            sel = f32(m) / (f32(rounds) * per_round);
        }
        for (var r = 0u; r < rounds; r = r + 1u) {
            let rkey = hash_u32(key_base ^ hash_u32(r * 0x85EBCA6Bu ^ 0x51ED27u));
            var mask = 1u + (hash_u32(rkey) % (len - 1u));
            if (lp) {
                mask = 1u + (hash_u32(rkey) % min(len - 1u, 7u));
            }
            for (var q = 0u; q < strides; q = q + 1u) {
                let i = lid.x + q * SORT_LOCAL;
                let j = i ^ mask;
                if (i < j && j < capn) {
                    let pkey = hash_u32(rkey ^ (i * 0x9E3779B9u));
                    if (rand01(params.seed_hi, params.seed_lo, pkey, 0u) < sel) {
                        try_pair(i, j, pkey, 8u, g_max, size_max, within, axis, body_center);
                    }
                }
            }
            storageBarrier();
            workgroupBarrier();
        }
        if (lid.x == 0u) {
            atomicAdd(&tallies[TALLY_PAR_CELLS_IDX], 1);
        }
        return;
    }

    if (lid.x != 0u) {
        return;
    }
    for (var c = 0u; c < m; c = c + 1u) {
        let draw0 = c * 8u;
        let ui = rand01(params.seed_hi, params.seed_lo, key_base, draw0);
        // j drawn uniformly from the n-1 slots != i (doc's fix -- the old
        // "bump j to j+1 when j==i" double-weighted i+1).
        let uj = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 1u);
        let i = min(u32(ui * f32(capn)), capn - 1u);
        let n_minus_1 = max(capn, 2u) - 1u;
        let j_off = min(u32(uj * f32(n_minus_1)), n_minus_1 - 1u);
        var j = (i + 1u + j_off) % capn;
        if (local_pairing()) {
            // §15.47: a neighbour 1-3 places away in the Morton order, on a
            // random side (the other side if that one runs off the list).
            let o = 1u + min(u32(uj * 6.0), 5u) % 3u;
            let up = u32(uj * 6.0) >= 3u;
            if (up && i + o < capn) {
                j = i + o;
            } else if (i >= o) {
                j = i - o;
            } else if (i + o < capn) {
                j = i + o;
            } else {
                j = (i + 1u) % capn;
            }
        }
        try_pair(i, j, key_base, draw0, g_max, size_max, within, axis, body_center);
    }
}

// One NTC candidate pair (sorted positions i != j in this cell's shared
// list): acceptance, then the collision. Shared by the serial loop and the
// parallel rounds of dense cells (doc §15.13).
// Doc §15.48 beam-fate tally (layout in dsmc_move.wgsl / gpu_types).
const FLAG_BEAM: u32 = 1u << 22u;
const FLAG_BEAM_COLL: u32 = 1u << 23u;
const TALLY_BEAM: u32 = 2989u;

// First photon-photon collision of a neutron's beam photon since its
// release: region (0 inside the emitter's envelope, 1 inside another
// body's, 2 outside every envelope) x partner (0 ambient, 1 the emitter's
// own emission, 2 another body's emission); for region 1 also the
// photon's distance from that body's axis.
fn beam_first(p: ptr<function, GpuPhoton>, partner_tag: u32) {
    if (((*p).c[0] & FLAG_BEAM) == 0u || ((*p).c[0] & FLAG_BEAM_COLL) != 0u) {
        return;
    }
    (*p).c[0] = (*p).c[0] | FLAG_BEAM_COLL;
    let e = ((*p).c[0] >> 17u) & 7u;
    let x = (*p).a.xyz;
    var region = 2u;
    var rho = 0.0;
    for (var b = 0u; b < params.n_bodies; b = b + 1u) {
        let d = x - bodies[b].center.xyz;
        if (length(d) < bodies[b].center.w) {
            if (b + 1u == e) {
                if (region == 2u) { region = 0u; }
            } else {
                region = 1u;
                let ax = bodies[b].axis.xyz;
                rho = length(d - ax * dot(d, ax));
            }
        }
    }
    var pc = 2u;
    if (partner_tag == 0u) { pc = 0u; } else if (partner_tag == e) { pc = 1u; }
    atomicAdd(&tallies[TALLY_BEAM + 4u + region * 3u + pc], 1);
    if (region == 1u) {
        atomicAdd(&tallies[TALLY_BEAM + 13u + min(u32(rho / 0.05), 9u)], 1);
    }
}

// Doc §15.49 near-exit counter: every collision of a beam photon (flag set
// from its release until its next capture) while it is within 1 of its
// emitter's envelope and within 45 deg of the emitter's axis. Bins of 0.1 in
// the distance outside the envelope x partner (0 ambient, 1 a fresh beam
// photon of the same emitter that has not yet collided, 2 any other photon
// of the same emitter, 3 another body's photon). Layout TALLY_BEAM+27+bin*4+k.
fn beam_near(p: ptr<function, GpuPhoton>, q_c0: u32) {
    if (((*p).c[0] & FLAG_BEAM) == 0u) {
        return;
    }
    let e = ((*p).c[0] >> 17u) & 7u;
    if (e == 0u || e > params.n_bodies) {
        return;
    }
    let d = (*p).a.xyz - bodies[e - 1u].center.xyz;
    let ax = bodies[e - 1u].axis.xyz;
    let z = abs(dot(d, ax));
    let rho = length(d - ax * dot(d, ax));
    let s = length(d) - bodies[e - 1u].center.w;
    if (s < 0.0 || s >= 1.0 || z < rho) {
        return;
    }
    let qe = (q_c0 >> 17u) & 7u;
    var k = 3u;
    if (qe == 0u) {
        k = 0u;
    } else if (qe == e) {
        k = select(2u, 1u, (q_c0 & FLAG_BEAM) != 0u && (q_c0 & FLAG_BEAM_COLL) == 0u);
    }
    atomicAdd(&tallies[TALLY_BEAM + 27u + min(u32(s / 0.1), 9u) * 4u + k], 1);
}

// Doc §15.52 pipe flags (dsmc_move.wgsl): PIPED = turned down a body's hole
// and not yet collided; LANE = piped in protected-lane mode (no collisions
// until it leaves the envelope).
const FLAG_PIPED: u32 = 1u << 24u;
const FLAG_LANE: u32 = 1u << 25u;
const TALLY_FAR: u32 = 3056u;
const TALLY_CHAN: u32 = 3080u;

// Doc §15.53 channel spin rule (halbach.pdf; test knob). Bit 14 of maze_mode
// on, eta = bits 16..23 / 255. Called for a pair holding a piped photon,
// instead of the scattering exchange. Nothing is randomised:
//  - co-moving (d_a.d_b >= 0): side by side, nothing happens;
//  - head-on, lab axial spins parallel (helicities opposite: photon meets
//    antiphoton, "opposite twice, and so the same"): they pass through and
//    spin each other up, eta of the relative linear energy going into spin;
//  - head-on, lab axial spins antiparallel (same species): they catch and
//    force one another apart, i.e. swap velocities (bounce back), and eta of
//    their outer spin energy goes to linear.
// The relative velocity alone is rescaled, so momentum is exact; the energy
// moved is exact (the restack keeps a level's energy). A bounced photon
// reverses with its spin unchanged, so its helicity flips (reflected photons
// are antiphotons).
fn chan_collide(pa: ptr<function, GpuPhoton>, pb: ptr<function, GpuPhoton>) {
    var bi = ((*pa).c[2] & 7u);
    if (bi == 0u) {
        bi = ((*pb).c[2] & 7u);
    }
    let slot = TALLY_CHAN + min(max(bi, 1u) - 1u, 3u) * 4u;
    let da = (*pa).b.xyz;
    let db = (*pb).b.xyz;
    if (dot(da, db) >= 0.0) {
        atomicAdd(&tallies[slot + 2u], 1);
        return;
    }
    let eta = f32((params.maze_mode >> 16u) & 255u) / 255.0;
    let kappa = params.rough_kappa;
    let ra = exp2(f32((*pa).c[3]));
    let rb = exp2(f32((*pb).c[3]));
    let ia = kappa * select(ra, 1.0, params.stacking != 0u) * ra * ra;
    let ib = kappa * select(rb, 1.0, params.stacking != 0u) * rb * rb;
    let va = da * (*pa).d.w;
    let vb = db * (*pb).d.w;
    let vcm = (va + vb) * 0.5;
    var u = va - vb;
    let wa0 = (*pa).d.xyz;
    let wb0 = (*pb).d.xyz;
    var wa = wa0;
    var wb = wb0;
    let ea = 0.5 * ia * dot(wa0, wa0);
    let eb = 0.5 * ib * dot(wb0, wb0);
    // Doc §15.55 edge mode (bits 24..25): 1 = no bounce, a spin-down pair
    // passes and its spin energy goes to linear along its own directions;
    // 2 = also reverse which pairing spins up (gear reading).
    let edge = (params.maze_mode >> 24u) & 3u;
    let up = ((*pa).b.w * (*pb).b.w < 0.0) != (edge == 2u);
    // Doc §15.54 quantum mode (bit 15; photon5.pdf "one particle loses a
    // spin and one gains", limat.pdf): a spin-up gives each photon one
    // saturated level, (1/2) kappa edge0^2, paid from the relative linear
    // energy if it can be afforded (else a plain pass); a spin-down strips
    // each photon's whole outer level into linear energy.
    let quantum = ((params.maze_mode >> 15u) & 1u) != 0u;
    let quant = 0.5 * kappa * params.spin_edge0 * params.spin_edge0;
    var de = 0.0;
    if (up) {
        // Pass through; relative linear energy (1/4)|u|^2 (unit masses)
        // gives eta of itself to the two spins, half each.
        let erel = 0.25 * dot(u, u);
        if (quantum) {
            de = select(0.0, 2.0 * quant, erel >= 2.0 * quant);
        } else {
            de = eta * erel;
        }
        u = u * sqrt(max(1.0 - de / max(erel, 1e-12), 0.0));
        var axa = select(sign((*pa).b.w) * da, wa0 / max(length(wa0), 1e-12), length(wa0) > 1e-9);
        var axb = select(sign((*pb).b.w) * db, wb0 / max(length(wb0), 1e-12), length(wb0) > 1e-9);
        wa = axa * sqrt(2.0 * (ea + 0.5 * de) / ia);
        wb = axb * sqrt(2.0 * (eb + 0.5 * de) / ib);
        atomicAdd(&tallies[slot], 1);
    } else {
        // Bounce (velocity swap), and eta of the outer spin energy to linear.
        // Quantum: only a stacked level (c[3] > 0) can be lost; a base
        // photon keeps its own axial spin.
        let keepa = select(1.0 - eta, select(1.0, 0.0, (*pa).c[3] > 0u), quantum);
        let keepb = select(1.0 - eta, select(1.0, 0.0, (*pb).c[3] > 0u), quantum);
        de = -(1.0 - keepa) * ea - (1.0 - keepb) * eb;
        wa = wa0 * sqrt(keepa);
        wb = wb0 * sqrt(keepb);
        let ul = length(u);
        u = select(-u, u, edge != 0u) * sqrt(max(ul * ul - 4.0 * de, 0.0)) / max(ul, 1e-12);
        atomicAdd(&tallies[slot + 1u], 1);
    }
    atomicAdd(&tallies[slot + 3u], i32(de * TALLY_SCALE));
    if (params.stacking != 0u) {
        wa = restack(pa, wa, wa - wa0);
        wb = restack(pb, wb, wb - wb0);
    }
    let nva = vcm + 0.5 * u;
    let nvb = vcm - 0.5 * u;
    let sa = length(nva);
    let sb = length(nvb);
    let nda = nva / max(sa, 1e-12);
    let ndb = nvb / max(sb, 1e-12);
    (*pa).b = vec4<f32>(nda, axial_helicity(pa, wa, nda));
    (*pb).b = vec4<f32>(ndb, axial_helicity(pb, wb, ndb));
    (*pa).d = vec4<f32>(wa, sa);
    (*pb).d = vec4<f32>(wb, sb);
}

// Doc §15.52: collisions within 1 outside body 0's envelope and within 45
// deg of its axis, by pole, s bin (<0.5, <1) and partner pair.
fn far_tally(x: vec3<f32>, ta: u32, tb: u32) {
    if (params.n_bodies == 0u) {
        return;
    }
    let d = x - bodies[0].center.xyz;
    let ax = bodies[0].axis.xyz;
    let z = dot(d, ax);
    let rho = length(d - ax * z);
    let s = length(d) - bodies[0].center.w;
    if (s < 0.0 || s >= 1.0 || abs(z) < rho) {
        return;
    }
    let ca = select(select(2u, 1u, ta == 1u), 0u, ta == 0u);
    let cb = select(select(2u, 1u, tb == 1u), 0u, tb == 0u);
    let lo = min(ca, cb);
    let hi = max(ca, cb);
    let k = select(select(5u, 3u + hi - 1u, lo == 1u), hi, lo == 0u);
    let pole = select(1u, 0u, z >= 0.0);
    atomicAdd(&tallies[TALLY_FAR + pole * 12u + select(1u, 0u, s < 0.5) * 6u + k], 1);
}

fn try_pair(i: u32, j: u32, key_base: u32, draw0: u32, g_max: f32, size_max: f32, within: bool, axis: vec3<f32>, body_center: vec3<f32>) {
    let slot_a = sh_slots[i];
    let slot_b = sh_slots[j];
    var pa = photons[slot_a];
    var pb = photons[slot_b];
    // Doc §15.11: photons held in a body's tank take no collisions.
    if (((pa.c[0] | pb.c[0]) & (1u << 20u)) != 0u) {
        return;
    }
    // Doc §15.52 protected lane: a lane photon takes no collisions.
    if (((pa.c[0] | pb.c[0]) & FLAG_LANE) != 0u) {
        return;
    }
    let da = pa.b.xyz;
    let db = pb.b.xyz;
    // E1 (doc §11d): pre-collision stiffness level/mass (m = 2^k). Read
    // here, before either the exchange or shedding below can change them,
    // since both the exchange math and NTC's size-scaling need the
    // PRE-collision values.
    let k_a = pa.c[3];
    let k_b = pb.c[3];
    let m_a = exp2(f32(k_a));
    let m_b = exp2(f32(k_b));

    // Relative VELOCITY (speed 1 unless RoughSphere+Free, where it
    // bit-matches the old length(da - db)).
    let g = length(da * pa.d.w - db * pb.d.w);
    var accept_p: f32 = select(1.0, g / g_max, params.rate_kernel == 1u);
    if (accept_p > 1.0) {
        atomicAdd(&tallies[TALLY_ACCEPT_OVERFLOW_IDX], 1);
    }
    if (params.sigma_size_scaling != 0u) {
        // sigma_pair / sigma_max = ((m1+m2)/2)^2 / size_max^2; sigma
        // itself cancels (doc §11d).
        let size_pair = (m_a + m_b) * 0.5;
        accept_p = accept_p * (size_pair * size_pair) / (size_max * size_max);
    }
    if (rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 2u) >= accept_p) {
        return;
    }

    let a_emitted = (pa.c[0] & FLAG_KIND_EMITTED) != 0u;
    let b_emitted = (pb.c[0] & FLAG_KIND_EMITTED) != 0u;
    if (params.d5_far_zone_gate != 0u && !a_emitted && !b_emitted && !within) {
        return;
    }

    // Doc §15.53: in channel spin mode a pair holding a piped photon spin-
    // exchanges instead of scattering (and is not counted as a collision).
    if (((params.maze_mode >> 14u) & 1u) != 0u && ((pa.c[0] | pb.c[0]) & FLAG_PIPED) != 0u) {
        chan_collide(&pa, &pb);
        photons[slot_a] = pa;
        photons[slot_b] = pb;
        return;
    }

    {
        let tag_a = min((pa.c[0] >> 17u) & 7u, 4u);
        let tag_b = min((pb.c[0] >> 17u) & 7u, 4u);
        atomicAdd(&tallies[TALLY_XCOLL_IDX + min(tag_a, tag_b) * 5u + max(tag_a, tag_b)], 1);
        if (tag_a != 0u && tag_b != 0u && tag_a != tag_b) {
            pa.c[0] = pa.c[0] | FLAG_XHIT;
            pb.c[0] = pb.c[0] | FLAG_XHIT;
        }
        let ca0 = pa.c[0];
        let cb0 = pb.c[0];
        beam_near(&pa, cb0);
        beam_near(&pb, ca0);
        beam_first(&pa, (pb.c[0] >> 17u) & 7u);
        beam_first(&pb, (pa.c[0] >> 17u) & 7u);
        far_tally(pa.a.xyz, (pa.c[0] >> 17u) & 7u, (pb.c[0] >> 17u) & 7u);
        // Doc §15.52: a collision ends a photon's guided run down a pipe.
        pa.c[0] = pa.c[0] & ~FLAG_PIPED;
        pb.c[0] = pb.c[0] & ~FLAG_PIPED;
    }

    if (params.exchange_rule == 2u) {
        // Doc §13 RoughSphere: contact mechanics replace the exchange,
        // the gear redirect, spin_transfer and shedding.
        rough_collide(&pa, &pb, key_base, draw0);
        atomicAdd(&tallies[TALLY_COLLISIONS_IDX], 1);
        if (max(i, j) < RANK64) {
            atomicAdd(&tallies[TALLY_LOW_RANK64_IDX], 1);
        }
        let bounces_a = pa.c[1];
        let bounces_b = pb.c[1];
        pa.c[1] = bounces_a + 1u;
        pb.c[1] = bounces_b + 1u;
        pa.c[0] = pa.c[0] | FLAG_COLLIDED_SINCE_CROSS | FLAG_PP_PASSAGE;
        pb.c[0] = pb.c[0] | FLAG_COLLIDED_SINCE_CROSS | FLAG_PP_PASSAGE;
        photons[slot_a] = pa;
        photons[slot_b] = pb;
        return;
    }

    let rel = da - db;
    let rel_len = length(rel);
    if (rel_len < 1e-6) {
        return;
    }

    // --- D6 exchange: SwapRenorm (legacy, |v|=c only) or RotateRelative
    // (conserves |v|=c AND d1+d2 exactly; doc §4 D6). ---
    var new_da: vec3<f32>;
    var new_db: vec3<f32>;
    var renorm_drift = 0.0;
    if (params.exchange_rule == 1u) {
        // RotateRelative, generalised for E1 stiffness (doc §4 D6, §11d).
        // Glancing-biased scattering angle: same u-skew as
        // sample_contact_normal (arp.pdf), mapped through the
        // hard-sphere impact-parameter relation b=sin(chi/2)... doc's
        // exact mapping: b = sqrt(u), chi = pi - 2*asin(b). Shared by
        // both the equal-mass and generalised branches below.
        let u0 = rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 3u);
        let one_minus_u0 = 1.0 - u0;
        let u = 1.0 - one_minus_u0 * one_minus_u0 * one_minus_u0;
        let b_impact = sqrt(clamp(u, 0.0, 1.0));
        let chi_ang = PI - 2.0 * asin(clamp(b_impact, 0.0, 1.0));
        let s = select(-1.0, 1.0, rand01(params.seed_hi, params.seed_lo, key_base, draw0 + 4u) < 0.5);

        // Gate-only negative control (doc §11d gate 1): force the
        // equal-mass branch even when k_a != k_b, to show the momentum
        // check below is sensitive to using the wrong (mass-blind) rule.
        let force_equal = params.debug_force_equal_mass_exchange != 0u;
        if (k_a == k_b || force_equal) {
            // Equal mass (includes the k=0 HEAD default: this branch is
            // BIT-IDENTICAL to the pre-E1 code). V = mean, g = relative;
            // |V|>0 case rotates g about V-hat (always perpendicular for
            // equal-speed vectors); the exact head-on case (|V|~=0)
            // rotates g about a random axis perpendicular to g instead.
            let v_vec = (da + db) * 0.5;
            let g_vec = da - db;
            let v_len = length(v_vec);

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
            // Unequal mass (E1, doc §11d): P = m1*d1 + m2*d2 is
            // conserved. d1 lies on the cone d1.Phat = const, so rotating
            // it about Phat by the SAME glancing-biased angle keeps
            // |d1'|=1 and leaves |P - m1*d1'| invariant (= m2 exactly,
            // mod f32), so d2' = (P - m1*d1')/m2 already has |d2'| ~= 1
            // before the defensive renormalise below (doc's "assert
            // pre-renorm |d2'| within 1e-4 of 1" -- checked via the
            // existing renorm_drift/TALLY_MOM_DRIFT_SUM tally, which is
            // exactly |raw_db|-1 here since raw_da is already unit-length
            // to f32 rounding). |P|~=0 is only possible when m1=m2
            // (excluded by this branch), so no separate near-zero case.
            let p_vec = m_a * da + m_b * db;
            let p_len = max(length(p_vec), 1e-6);
            let phat = p_vec / p_len;
            let theta = s * chi_ang;
            let cos_t = cos(theta);
            let sin_t = sin(theta);
            // Full Rodrigues rotation (da is NOT generally perpendicular
            // to phat, unlike the equal-mass g_vec/vhat pair above).
            let raw_da = da * cos_t + cross(phat, da) * sin_t + phat * dot(phat, da) * (1.0 - cos_t);
            new_da = normalize(raw_da);
            let raw_db = (p_vec - m_a * new_da) / m_b;
            new_db = normalize(raw_db);
            renorm_drift = length(raw_da - new_da) + length(raw_db - new_db);
        }
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

    // E1 (doc §11d): per-photon deflection angle from the exchange step
    // alone (before the gear redirect below can rotate it further),
    // bucketed by that photon's OWN pre-collision stiffness level --
    // e1_heavy_deflects_less compares these buckets directly. Tallied
    // for both exchange rules (SwapRenorm collisions land in bucket 0
    // since unshed/ambient photons stay at k=0).
    let defl_a = acos(clamp(dot(da, new_da), -1.0, 1.0));
    let defl_b = acos(clamp(dot(db, new_db), -1.0, 1.0));
    let bucket_a = min(k_a, 7u);
    let bucket_b = min(k_b, 7u);
    atomicAdd(&tallies[TALLY_DEFLECT_SUM_IDX + bucket_a], i32(defl_a * TALLY_SCALE));
    atomicAdd(&tallies[TALLY_DEFLECT_SUMSQ_IDX + bucket_a], i32(defl_a * defl_a * TALLY_SCALE));
    atomicAdd(&tallies[TALLY_DEFLECT_COUNT_IDX + bucket_a], 1);
    atomicAdd(&tallies[TALLY_DEFLECT_SUM_IDX + bucket_b], i32(defl_b * TALLY_SCALE));
    atomicAdd(&tallies[TALLY_DEFLECT_SUMSQ_IDX + bucket_b], i32(defl_b * defl_b * TALLY_SCALE));
    atomicAdd(&tallies[TALLY_DEFLECT_COUNT_IDX + bucket_b], 1);

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
        pa.c[0] = pa.c[0] | FLAG_COLLIDED_SINCE_CROSS | FLAG_PP_PASSAGE;
    }
    if (cap == 0u || bounces_b < cap) {
        pb.b = vec4<f32>(new_db, clamp(pb.b.w + spin_delta_b, -2.0, 2.0));
        pb.c[1] = bounces_b + 1u;
        pb.c[0] = pb.c[0] | FLAG_COLLIDED_SINCE_CROSS | FLAG_PP_PASSAGE;
    }

    // E1 shedding (doc §11d, Transfer mode): move one stacked level from
    // the higher-k photon to the lower-k one, AFTER the exchange (which
    // used the pre-collision k_a/k_b/m_a/m_b above) -- conserves the
    // total stacked-level count across the box. Equal k: no change.
    // Independent of the bounce cap (k is a stiffness bookkeeping field,
    // not a ricochet-limited direction update).
    if (params.stiffness_shed == 1u) {
        if (k_a > k_b) {
            pa.c[3] = k_a - 1u;
            pb.c[3] = k_b + 1u;
        } else if (k_b > k_a) {
            pb.c[3] = k_b - 1u;
            pa.c[3] = k_a + 1u;
        }
    }

    photons[slot_a] = pa;
    photons[slot_b] = pb;
}
