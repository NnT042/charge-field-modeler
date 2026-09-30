// DSMC step 2: advance every live photon by c*dt, test the nearest body for
// absorption/lane-crossing, apply the domain boundary, and emit this
// photon's cell index (dead -> bucket `num_cells`) for CellSorter.
// dsmc_common.wgsl is concatenated ahead of this source by engine.rs.

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
    body_beta: f32, body_beta_twist: f32, _pad14a: u32, _pad14b: u32,
};

const FLAG_ALIVE: u32 = 1u;
const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const SOURCE_SHIFT: u32 = 3u;
const SOURCE_MASK: u32 = 3u << SOURCE_SHIFT;
const CROSSED_SHIFT: u32 = 5u;
const ARRIVED_SHIFT: u32 = 9u;
const FLAG_COLLIDED_SINCE_CROSS: u32 = 1u << 13u;
const HIST_BINS: u32 = 36u;
const HIST_MAX_BODIES: u32 = 4u;
const HIST_ARRIVAL_BLOCK: u32 = HIST_MAX_BODIES * HIST_BINS;
const HIST_EMISSION_BLOCK_MV: u32 = 0u;
const TALLY_SCALE: f32 = 4096.0;
// B2b item 3 (doc §8) tally offsets -- must stay in sync with gpu_types.rs.
const TALLY_ABSORBED_MOM: u32 = 72u;
const TALLY_ESCAPED_BAND: u32 = 84u;
const TALLY_ESCAPED_BAND_SPIN: u32 = 114u;
const TALLY_ESCAPED_CROSSED_BODY: u32 = 144u;
const TALLY_ESCAPED_CROSSED_BODY_SPIN: u32 = 152u;
const TALLY_SURVIVAL: u32 = 160u;
// E1 (doc §11d) -- must stay in sync with gpu_types.rs.
const TALLY_ESCAPED_EMITTED_K_SUM: u32 = 192u;
const TALLY_ESCAPED_EMITTED_K_COUNT: u32 = 193u;
// Doc §13 -- must stay in sync with gpu_types.rs.
const TALLY_ESC_ALL: u32 = 222u;
const TALLY_ESC_SPEED_SUM: u32 = 223u;
const TALLY_ESC_SPEED_SUMSQ: u32 = 224u;
const TALLY_ESC_EDGE_SUM: u32 = 225u;
const TALLY_ESC_EDGE_SUMSQ: u32 = 226u;
const TALLY_ESC_EDGE_OVER_C: u32 = 227u;
const TALLY_TWIST_TO_NORMAL_IDX: u32 = 228u;
const TALLY_ESC_LEVEL: u32 = 235u; // [k], len 4
// Doc §14 shells -- must stay in sync with gpu_types.rs.
const TALLY_SH_ENTRIES: u32 = 240u;      // [body]
const TALLY_SH_EXITS: u32 = 244u;        // [body], exits after >= 1 contact
const TALLY_SH_THROUGH: u32 = 248u;      // [body], pole-to-pole passages
const TALLY_SH_CONTACTS: u32 = 252u;     // [body], contact count
const TALLY_SH_CLEAN: u32 = 256u;        // [body], passages with no contact
const TALLY_SH_EXIT_BAND: u32 = 260u;    // [body*3 + band]
const TALLY_SH_EXIT_BAND_HEL: u32 = 272u;// [body*3 + band], helicity sum (fixed)
const TALLY_SH_DP: u32 = 284u;           // [body*3 + axis], momentum to body (fixed)
const TALLY_SH_TORQUE: u32 = 296u;       // [body*3 + axis], (x - c) x dp (fixed)
const TALLY_SH_WORK: u32 = 308u;         // [body], energy given to photons (fixed)
// Passage state packed in c[2] (birth_step is unused): bits 0-2 = inside
// body index + 1 (0 = outside), bits 3-8 = entry latitude bin, bits 9-12 =
// entry angle-to-pole bin (doc §15), bits 13-31 = contacts this passage.
const TALLY_BALL_OVERLAP: u32 = 312u;
const TALLY_PASS: u32 = 313u;         // [(abin*4 + outcome)*2 + spec], body 0, len 72
const TALLY_PASS_THROUGH: u32 = 385u; // [abin*2 + spec], body 0, len 18
const TALLY_PASS_HEL: u32 = 403u;     // [outcome*2 + spec], body 0, len 8
const FLAG_PP_PASSAGE: u32 = 1u << 14u;
const FLAG_HIT_CLEAN: u32 = 1u << 15u;    // first ball hit came before any pp collision
const FLAG_HIT_AFTER_PP: u32 = 1u << 16u; // first ball hit came after a pp collision

@group(0) @binding(0) var<storage, read_write> photons: array<GpuPhoton>;
@group(0) @binding(1) var<storage, read> bodies: array<GpuBody>;
@group(0) @binding(2) var<storage, read_write> cell_of: array<u32>;
@group(0) @binding(3) var<storage, read_write> tallies: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> hist: array<atomic<u32>>;
@group(0) @binding(5) var<storage, read_write> density: array<atomic<u32>>;
// Doc §14 shells: per body 8 vec4 -- [2k] = (origin.xyz, radius),
// [2k+1] = (omega.xyz, 0) for levels k = 0..3 (9..12), world space,
// uploaded every step from baryon_shells::ShellBody::frames.
@group(0) @binding(6) var<storage, read> shells: array<vec4<f32>>;
// Doc §15 ball: per body 9 keyframes x 3 vec4 -- (centre.xyz, radius),
// (omega.xyz, 0), (true centre velocity.xyz, 0) at sub-step s*dt/8.
@group(0) @binding(7) var<storage, read> ballkeys: array<vec4<f32>>;
@group(0) @binding(8) var<uniform> params: Params;

fn hist_bin(lat_deg: f32) -> u32 {
    let raw = i32(floor((lat_deg + 90.0) / 5.0));
    return u32(clamp(raw, 0, i32(HIST_BINS) - 1));
}

fn species_idx(flags: u32) -> u32 {
    return select(0u, 1u, (flags & FLAG_SPECIES_PLUS) != 0u);
}

// Doc §8: 0 = polar (>60 deg), 1 = mid (30-60 deg), 2 = equatorial (<30 deg).
fn escape_band(abs_lat_deg: f32) -> u32 {
    if (abs_lat_deg > 60.0) {
        return 0u;
    } else if (abs_lat_deg >= 30.0) {
        return 1u;
    }
    return 2u;
}


fn shell_origin(b: u32, k: u32) -> vec3<f32> { return shells[b * 8u + 2u * k].xyz; }
fn shell_radius(b: u32, k: u32) -> f32 { return shells[b * 8u + 2u * k].w; }
fn shell_omega(b: u32, k: u32) -> vec3<f32> { return shells[b * 8u + 2u * k + 1u].xyz; }

// Innermost shell of body b containing x, or 4 when outside the body.
fn shell_region(b: u32, x: vec3<f32>) -> u32 {
    for (var k = 0u; k < 4u; k = k + 1u) {
        let d = x - shell_origin(b, k);
        let r = shell_radius(b, k);
        if (dot(d, d) <= r * r) {
            return k;
        }
    }
    return 4u;
}

// Matter velocity in shell k: sum_{j>=k} w_j x (x - o_j)
// (baryon_shells::matter_velocity).
fn shell_velocity(b: u32, k: u32, x: vec3<f32>) -> vec3<f32> {
    var v = vec3<f32>(0.0);
    for (var j = k; j < 4u; j = j + 1u) {
        v = v + cross(shell_omega(b, j), x - shell_origin(b, j));
    }
    return v;
}

fn body_of(x: vec3<f32>) -> u32 {
    for (var i = 0u; i < params.n_bodies; i = i + 1u) {
        if (params.body_model == 2u) {
            // Ball mode: the envelope (level-12 ghost) only does bookkeeping.
            let d = x - bodies[i].center.xyz;
            if (dot(d, d) <= bodies[i].center.w * bodies[i].center.w) {
                return i;
            }
        } else if (shell_region(i, x) < 4u) {
            return i;
        }
    }
    return 4u;
}

fn bk(b: u32, s: u32, j: u32) -> vec4<f32> { return ballkeys[(b * 9u + s) * 3u + j]; }

// Mirrors rough_sphere::moving_sphere_hit.
fn moving_sphere_hit(x0: vec3<f32>, v: vec3<f32>, c0: vec3<f32>, vc: vec3<f32>, r: f32, t_max: f32) -> f32 {
    let d = x0 - c0;
    let u = v - vc;
    let a = dot(u, u);
    let bq = dot(d, u);
    let c = dot(d, d) - r * r;
    if (c <= 0.0 || a < 1e-18 || bq >= 0.0) {
        return -1.0;
    }
    let disc = bq * bq - a * c;
    if (disc < 0.0) {
        return -1.0;
    }
    let t = (-bq - sqrt(disc)) / a;
    if (t >= 0.0 && t <= t_max) {
        return t;
    }
    return -1.0;
}

fn body_lat_deg(b: u32, x: vec3<f32>) -> f32 {
    let body = bodies[b];
    let d = x - body.center.xyz;
    let colat = acos(clamp(dot(d, body.axis.xyz) / max(length(d), 1e-6), -1.0, 1.0));
    return 90.0 - degrees(colat);
}

// Doc §14/§14.8: one photon vs body-grain contact. Mirrors
// rough_sphere::wall_contact_spinning line for line (an infinitely massive
// grain moving at gv and spinning at gw, unit radius), with the body contact
// law (body_beta / body_beta_twist / body_twist_loss), then the §13.7
// restack. Stacking mass model: m = 1 when stacking, else m = R = 2^k.
// Returns the momentum given to the body.
fn wall_contact(p: ptr<function, GpuPhoton>, gv: vec3<f32>, gw: vec3<f32>, n: vec3<f32>) -> vec3<f32> {
    return wall_contact_general(p, gv, gw, gw, n);
}

// rough_sphere::wall_contact_general: the partner surface slides at
// gv + sw x n; the twist grip acts on (w - tw).n.
fn wall_contact_general(p: ptr<function, GpuPhoton>, gv: vec3<f32>, sw: vec3<f32>, tw: vec3<f32>, n: vec3<f32>) -> vec3<f32> {
    let r = exp2(f32((*p).c[3]));
    let m = select(r, 1.0, params.stacking != 0u);
    let kappa = params.rough_kappa;
    let inertia = kappa * m * r * r;
    var v = (*p).b.xyz * (*p).d.w;
    var w = (*p).d.xyz;
    let v_in = v;
    let w_in = w;
    let e_rot_in = 0.5 * inertia * dot(w, w);
    let g = (v - gv) + r * cross(w, n) + cross(sw, n);
    let gn = dot(g, n);
    if (gn <= 0.0) {
        return vec3<f32>(0.0);
    }
    let g_t = g - gn * n;
    let j_n = -m * 2.0 * gn * n; // e = 1
    let k_t = (1.0 / m) * (1.0 + 1.0 / kappa);
    let j_t = -(1.0 + params.body_beta) * g_t / k_t;
    let j = j_n + j_t;
    v = v + j / m;
    let nj = cross(n, j);
    var kick = r * nj / inertia;
    w = w + kick;
    let omega_n = dot(w - tw, n);
    let e_before = 0.5 * inertia * dot(w, w);
    let dtw = -(1.0 + params.body_beta_twist) * omega_n * n;
    w = w + dtw;
    kick = kick + dtw;
    let lost = e_before - 0.5 * inertia * dot(w, w);
    if (params.body_twist_loss != 0u && lost > 0.0) {
        let w_n = dot(v - gv, n);
        let jx = m * (-w_n - sqrt(w_n * w_n + 2.0 * lost / m));
        v = v + jx * n / m;
        atomicAdd(&tallies[TALLY_TWIST_TO_NORMAL_IDX], i32(lost * TALLY_SCALE));
    }
    if (params.speed_mode == 1u) {
        v = normalize(v);
        let e_rot = 0.5 * inertia * dot(w, w);
        if (e_rot > 1e-12) {
            let sc = sqrt(e_rot_in / e_rot);
            kick = (sc - 1.0) * w_in + sc * kick;
            w = w * sc;
        } else {
            let fl = length(nj);
            if (fl > 1e-9 && e_rot_in > 1e-12) {
                w = nj / fl * sqrt(2.0 * e_rot_in / inertia);
                kick = w - w_in;
            }
        }
    }
    let s = length(v);
    let d = v / max(s, 1e-12);
    if (params.stacking != 0u) {
        w = restack(p, w, kick);
    }
    (*p).b = vec4<f32>(d, axial_helicity(p, w, d));
    (*p).d = vec4<f32>(w, s);
    return -(v - v_in) * m;
}

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let idx = linear_wg * 256u + lid.x;
    if (idx >= params.n_live) {
        return;
    }

    var ph = photons[idx];
    if ((ph.c[0] & FLAG_ALIVE) == 0u) {
        cell_of[idx] = params.num_cells;
        return;
    }

    // Doc §14 shells: body-grain contact at the start of the step, before
    // the photon moves (so this step's path uses the post-contact velocity).
    let shells_mode = params.body_model == 1u;
    if (shells_mode) {
        let x0 = ph.a.xyz;
        let bi = body_of(x0);
        if (bi < 4u) {
            let k = shell_region(bi, x0);
            let gv = shell_velocity(bi, k, x0);
            let u = ph.b.xyz * ph.d.w - gv;
            let ulen = length(u);
            let pid = bitcast<u32>(ph.a.w);
            let key = hash_u32(pid ^ hash_u32(params.step ^ 0xB0D1E5u));
            let p_hit = 1.0 - exp(-params.tau_body * ulen * params.dt);
            if (ulen > 1e-6 && rand01(params.seed_hi, params.seed_lo, key, 0u) < p_hit) {
                var rho = rand01(params.seed_hi, params.seed_lo, key, 1u);
                if (params.glancing_bias != 0u) {
                    let om = 1.0 - rho;
                    rho = 1.0 - om * om * om;
                }
                let psi = rand01(params.seed_hi, params.seed_lo, key, 2u) * TAU;
                let uh = u / ulen;
                let helper = perp_helper(uh);
                let t = normalize(cross(uh, helper));
                let bp = cross(uh, t);
                let e_perp = t * cos(psi) + bp * sin(psi);
                let n = normalize(uh * sqrt(max(1.0 - rho, 0.0)) + e_perp * sqrt(rho));
                let v_before = ph.b.xyz * ph.d.w;
                let e_before = 0.5 * dot(v_before, v_before);
                // §14.8: the grain spins with its shell (direction of the
                // shell matter's total angular velocity) at edge speed
                // grain_edge on unit radius. grain_edge = 0: non-spinning.
                var omega_k = vec3<f32>(0.0);
                for (var j = k; j < 4u; j = j + 1u) {
                    omega_k = omega_k + shell_omega(bi, j);
                }
                var gw = vec3<f32>(0.0);
                if (length(omega_k) > 1e-9) {
                    gw = normalize(omega_k) * params.grain_edge;
                }
                let dp = wall_contact(&ph, gv, gw, n);
                if (dot(dp, dp) > 0.0) {
                    let bc = bodies[bi].center.xyz;
                    let tq = cross(x0 - bc, dp);
                    for (var a = 0u; a < 3u; a = a + 1u) {
                        atomicAdd(&tallies[TALLY_SH_DP + bi * 3u + a], i32(dp[a] * TALLY_SCALE));
                        atomicAdd(&tallies[TALLY_SH_TORQUE + bi * 3u + a], i32(tq[a] * TALLY_SCALE));
                    }
                    let v_after = ph.b.xyz * ph.d.w;
                    atomicAdd(&tallies[TALLY_SH_WORK + bi], i32((0.5 * dot(v_after, v_after) - e_before) * TALLY_SCALE));
                    atomicAdd(&tallies[TALLY_SH_CONTACTS + bi], 1);
                    // Count the contact against this passage (bits 9-31).
                    ph.c[2] = ph.c[2] + (1u << 13u);
                }
            }
        }
    }

    var pos = ph.a.xyz;
    let dir = ph.b.xyz;
    let id = bitcast<u32>(ph.a.w);
    let kind_emitted = (ph.c[0] & FLAG_KIND_EMITTED) != 0u;
    let spec = species_idx(ph.c[0]);
    let source_body = (ph.c[0] & SOURCE_MASK) >> SOURCE_SHIFT;

    // Path length this step: |v|*dt (|v| = 1 except RoughSphere+Free).
    let speed = ph.d.w;
    var remaining = params.dt * speed;

    // Doc §15 ball: exact hit test against each body's moving base ball over
    // 8 linear sub-intervals, one hit per step. On a hit the photon moves to
    // the contact point, takes the contact, and flies the rest of the step on
    // its new velocity.
    var ball_moved = false;
    if (params.body_model == 2u) {
        let v0 = ph.b.xyz * ph.d.w;
        let h = params.dt / 8.0;
        var best_t = 1e30;
        var best_b = 4u;
        var best_s = 0u;
        for (var b = 0u; b < params.n_bodies; b = b + 1u) {
            let d0 = pos - bk(b, 0u, 0u).xyz;
            let rb0 = bk(b, 0u, 0u).w;
            if (dot(d0, d0) < rb0 * rb0) {
                atomicAdd(&tallies[TALLY_BALL_OVERLAP], 1);
            }
            for (var si = 0u; si < 8u; si = si + 1u) {
                let t0 = f32(si) * h;
                let c0 = bk(b, si, 0u).xyz;
                let c1 = bk(b, si + 1u, 0u).xyz;
                let vc = (c1 - c0) / h;
                let t = moving_sphere_hit(pos + v0 * t0, v0, c0, vc, bk(b, si, 0u).w, h);
                if (t >= 0.0 && t0 + t < best_t) {
                    best_t = t0 + t;
                    best_b = b;
                    best_s = si;
                }
                if (best_b == b) {
                    break; // later sub-intervals of this ball are later in time
                }
            }
        }
        if (best_b < 4u) {
            let b = best_b;
            let si = best_s;
            let f = (best_t - f32(si) * h) / h;
            let xh = pos + v0 * best_t;
            let ch = mix(bk(b, si, 0u).xyz, bk(b, si + 1u, 0u).xyz, f);
            let rb = bk(b, si, 0u).w;
            let om = mix(bk(b, si, 1u).xyz, bk(b, si + 1u, 1u).xyz, f);
            let vcen = mix(bk(b, si, 2u).xyz, bk(b, si + 1u, 2u).xyz, f);
            let n = normalize(ch - xh);
            var tw = vec3<f32>(0.0);
            if (length(om) > 1e-9) {
                tw = normalize(om) * params.grain_edge;
            }
            let pp_before = (ph.c[0] & FLAG_PP_PASSAGE) != 0u;
            let dp = wall_contact_general(&ph, vcen, rb * om, tw, n);
            if (dot(dp, dp) > 0.0) {
                let bc = bodies[b].center.xyz;
                let tq = cross(xh - bc, dp);
                for (var a = 0u; a < 3u; a = a + 1u) {
                    atomicAdd(&tallies[TALLY_SH_DP + b * 3u + a], i32(dp[a] * TALLY_SCALE));
                    atomicAdd(&tallies[TALLY_SH_TORQUE + b * 3u + a], i32(tq[a] * TALLY_SCALE));
                }
                atomicAdd(&tallies[TALLY_SH_CONTACTS + b], 1);
                if ((ph.c[2] >> 13u) == 0u && (ph.c[0] & (FLAG_HIT_CLEAN | FLAG_HIT_AFTER_PP)) == 0u) {
                    ph.c[0] = ph.c[0] | select(FLAG_HIT_CLEAN, FLAG_HIT_AFTER_PP, pp_before);
                }
                ph.c[2] = ph.c[2] + (1u << 13u);
            }
            let v1 = ph.b.xyz * ph.d.w;
            pos = xh + v1 * (params.dt - best_t);
            ball_moved = true;
        }
    }

    var hit_body: i32 = -1;
    var hit_t = remaining;
    // Absorber bodies only; shells are transparent regions.
    let n_absorbers = select(params.n_bodies, 0u, shells_mode || params.body_model == 2u);
    for (var i = 0u; i < n_absorbers; i = i + 1u) {
        let b = bodies[i];
        let t = ray_sphere_hit(pos, dir, b.center.xyz, b.center.w, remaining);
        if (t >= 0.0 && t < hit_t) {
            hit_t = t;
            hit_body = i32(i);
        }
    }

    var dead = false;
    var crossed_any = false;

    if (hit_body >= 0) {
        let i = u32(hit_body);
        let body = bodies[i];
        let hit_pos = pos + dir * hit_t;
        let r_in = body.center.w;

        // Counter-stream survival (doc §8, "measured not integrated"): BEFORE
        // this event changes any flags, snapshot which lanes this photon has
        // crossed and whether it has collided since its last crossing. A
        // photon that crossed body k's lane earlier and now ARRIVES at a
        // DIFFERENT body i (whether it is absorbed here or crosses again) is
        // one sample of "did it survive the gap uninterrupted".
        let crossed_mask_before = (ph.c[0] >> CROSSED_SHIFT) & 0xFu;
        let collided_before = select(0u, 1u, (ph.c[0] & FLAG_COLLIDED_SINCE_CROSS) != 0u);
        for (var k = 0u; k < 4u; k = k + 1u) {
            if (((crossed_mask_before >> k) & 1u) != 0u && k != i) {
                atomicAdd(&tallies[TALLY_SURVIVAL + (k * 4u + i) * 2u + collided_before], 1);
            }
        }

        // Arrival latitude in the hit body's own frame (ambient only, per doc
        // §9 G0). Tallied once per body (first contact only), matching the
        // CPU reference's `probe_arrived` guard in `march_stack`.
        let arrived_bit = 1u << (ARRIVED_SHIFT + i);
        if (!kind_emitted && (ph.c[0] & arrived_bit) == 0u) {
            ph.c[0] = ph.c[0] | arrived_bit;
            let d0 = hit_pos - body.center.xyz;
            let colat0 = acos(clamp(dot(d0, body.axis.xyz) / r_in, -1.0, 1.0));
            let lat0 = 90.0 - degrees(colat0);
            atomicAdd(&hist[HIST_ARRIVAL_BLOCK + i * HIST_BINS + hist_bin(lat0)], 1u);
        }

        let through_lane_rho = body.lane.x;
        let cos_max_angle = body.lane.y;
        var crossed = false;
        // Doc §2/§7: "a photon may cross body i's lane unless it was emitted
        // by body i" -- applies to every photon, not just ambient (an
        // emitted photon from a NEIGHBOUR must be able to cross this body's
        // lane, which is the whole through-exhaust mechanism a p-n-p stack
        // measures). Only self-emission is excluded.
        let self_lane = kind_emitted && source_body == i;
        if (!self_lane && through_lane_rho > 0.0) {
            let d = hit_pos - body.center.xyz;
            let axial = dot(dir, body.axis.xyz);
            let along = dot(d, body.axis.xyz);
            let rho = length(d - body.axis.xyz * along);
            let incoming = dot(dir, d / r_in);
            if (rho < through_lane_rho * r_in && abs(axial) > cos_max_angle && incoming < 0.0) {
                let t2 = -2.0 * dot(d, dir);
                if (t2 > 0.0) {
                    let exitv = d + dir * t2;
                    let exit_len = max(length(exitv), 1e-6);
                    pos = body.center.xyz + exitv * (r_in / exit_len) * 1.001;
                    crossed = true;
                    crossed_any = true;
                    ph.c[0] = ph.c[0] | (1u << (CROSSED_SHIFT + i));
                    // Doc §8: crossing a lane CLEARS "collided since last
                    // crossing" -- the counter-stream survival window resets.
                    ph.c[0] = ph.c[0] & ~FLAG_COLLIDED_SINCE_CROSS;
                    atomicAdd(&tallies[32 + i * 2 + spec], 1); // TALLY_CROSSINGS
                }
            }
        }

        if (!crossed) {
            // Absorbed.
            dead = true;
            let already_crossed = (ph.c[0] >> CROSSED_SHIFT) != 0u;
            let kind_i = select(0u, 1u, kind_emitted);
            atomicAdd(&tallies[i * 4u + kind_i * 2u + spec], 1); // TALLY_ABSORBED
            if (already_crossed) {
                atomicAdd(&tallies[16 + i * 4u + kind_i * 2u + spec], 1); // TALLY_ABSORBED_CROSSED
            }
            // Absorbed momentum (doc §8): Sigma dir.xyz per body -- the
            // open-field absorber force observable.
            atomicAdd(&tallies[TALLY_ABSORBED_MOM + i * 3u + 0u], i32(dir.x * TALLY_SCALE));
            atomicAdd(&tallies[TALLY_ABSORBED_MOM + i * 3u + 1u], i32(dir.y * TALLY_SCALE));
            atomicAdd(&tallies[TALLY_ABSORBED_MOM + i * 3u + 2u], i32(dir.z * TALLY_SCALE));
        } else {
            // Continue the remaining path length from the far side (dt is
            // small relative to body spacing, so a single hit test per step
            // is an adequate approximation -- see engine.rs doc comment).
            pos = pos + dir * (remaining - hit_t);
        }
    } else if (!ball_moved) {
        pos = pos + dir * remaining;
    }

    // Doc §14/§15 passage bookkeeping after the move.
    if ((shells_mode || params.body_model == 2u) && !dead) {
        let prev_inside = ph.c[2] & 7u;
        let now_body = body_of(pos);
        let now_inside = select(0u, now_body + 1u, now_body < 4u);
        if (prev_inside != now_inside) {
            if (prev_inside != 0u) {
                let b = prev_inside - 1u;
                let contacts = ph.c[2] >> 13u;
                if (b == 0u) {
                    // Doc §15 passage outcome by entry angle to the pole.
                    let abin = (ph.c[2] >> 9u) & 15u;
                    var outcome = 0u;
                    if ((ph.c[0] & FLAG_HIT_CLEAN) != 0u) {
                        outcome = 2u;
                    } else if ((ph.c[0] & FLAG_HIT_AFTER_PP) != 0u) {
                        outcome = 3u;
                    } else if ((ph.c[0] & FLAG_PP_PASSAGE) != 0u) {
                        outcome = 1u;
                    }
                    atomicAdd(&tallies[TALLY_PASS + (abin * 4u + outcome) * 2u + spec], 1);
                    atomicAdd(&tallies[TALLY_PASS_HEL + outcome * 2u + spec], i32(ph.b.w * TALLY_SCALE));
                    let lat_in0 = f32((ph.c[2] >> 3u) & 63u) * 5.0 - 87.5;
                    let lat_out0 = body_lat_deg(b, pos);
                    if ((lat_in0 > 60.0 && lat_out0 < -60.0) || (lat_in0 < -60.0 && lat_out0 > 60.0)) {
                        atomicAdd(&tallies[TALLY_PASS_THROUGH + abin * 2u + spec], 1);
                    }
                }
                if (contacts == 0u) {
                    atomicAdd(&tallies[TALLY_SH_CLEAN + b], 1);
                } else {
                    atomicAdd(&tallies[TALLY_SH_EXITS + b], 1);
                    let lat_out = body_lat_deg(b, pos);
                    let band = escape_band(abs(lat_out));
                    atomicAdd(&tallies[TALLY_SH_EXIT_BAND + b * 3u + band], 1);
                    atomicAdd(&tallies[TALLY_SH_EXIT_BAND_HEL + b * 3u + band], i32(ph.b.w * TALLY_SCALE));
                    atomicAdd(&hist[HIST_EMISSION_BLOCK_MV + b * HIST_BINS + hist_bin(lat_out)], 1u);
                    let lat_in = f32((ph.c[2] >> 3u) & 63u) * 5.0 - 87.5;
                    if ((lat_in > 60.0 && lat_out < -60.0) || (lat_in < -60.0 && lat_out > 60.0)) {
                        atomicAdd(&tallies[TALLY_SH_THROUGH + b], 1);
                    }
                }
            }
            var packed = 0u;
            if (now_inside != 0u) {
                let b = now_inside - 1u;
                let lat_in = body_lat_deg(b, pos);
                atomicAdd(&tallies[TALLY_SH_ENTRIES + b], 1);
                atomicAdd(&hist[HIST_ARRIVAL_BLOCK + b * HIST_BINS + hist_bin(lat_in)], 1u);
                let ax = abs(dot(ph.b.xyz, bodies[b].axis.xyz));
                let abin = min(u32(degrees(acos(clamp(ax, 0.0, 1.0))) / 10.0), 8u);
                packed = now_inside | (hist_bin(lat_in) << 3u) | (abin << 9u);
            }
            ph.c[2] = packed;
            ph.c[0] = ph.c[0] & ~(FLAG_PP_PASSAGE | FLAG_HIT_CLEAN | FLAG_HIT_AFTER_PP);
        }
    }

    if (!dead) {
        let center = vec3<f32>(params.domain_center_x, params.domain_center_y, params.domain_center_z);
        if (params.boundary_mode == 1u) {
            // PeriodicBox: wrap.
            let half = params.boundary_extent;
            let full = 2.0 * half;
            var rel = pos - center;
            rel = rel - full * floor((rel + half) / full);
            pos = center + rel;
        } else {
            if (length(pos - center) > params.boundary_extent) {
                dead = true;
                let kind_i = select(0u, 1u, kind_emitted);
                let crossed_i = select(0u, 1u, crossed_any || (ph.c[0] >> CROSSED_SHIFT) != 0u);
                atomicAdd(&tallies[40 + source_body * 4u + crossed_i * 2u + spec], 1); // TALLY_ESCAPED (legacy shape, kept)

                // Doc §8: escape topology by source-body-slot x species x
                // |latitude| band. Emitted -> that body's own frame; ambient
                // (no source body) -> the stack frame (axis Z through the
                // domain centroid, per doc "ambient source = none -> use
                // band in stack frame").
                var band_axis = vec3<f32>(0.0, 0.0, 1.0);
                var band_center = center;
                var slot = 4u; // "none"
                if (kind_emitted) {
                    let sb = bodies[source_body];
                    band_axis = sb.axis.xyz;
                    band_center = sb.center.xyz;
                    slot = source_body;
                }
                let rel_e = pos - band_center;
                let rlen = max(length(rel_e), 1e-6);
                let colat_e = acos(clamp(dot(rel_e, band_axis) / rlen, -1.0, 1.0));
                let lat_e = 90.0 - degrees(colat_e);
                let band = escape_band(abs(lat_e));
                let eb_idx = slot * 2u * 3u + spec * 3u + band;
                atomicAdd(&tallies[TALLY_ESCAPED_BAND + eb_idx], 1);
                atomicAdd(&tallies[TALLY_ESCAPED_BAND_SPIN + eb_idx], i32(ph.b.w * TALLY_SCALE));

                // Of escaped photons, per-crossed-body-i breakdown (doc §8:
                // "keep crossed body i for i=0..3"). NOT mutually exclusive:
                // a photon that crossed more than one lane increments more
                // than one counter.
                let cmask = (ph.c[0] >> CROSSED_SHIFT) & 0xFu;
                for (var k = 0u; k < 4u; k = k + 1u) {
                    if (((cmask >> k) & 1u) != 0u) {
                        atomicAdd(&tallies[TALLY_ESCAPED_CROSSED_BODY + k * 2u + spec], 1);
                        atomicAdd(&tallies[TALLY_ESCAPED_CROSSED_BODY_SPIN + k * 2u + spec], i32(ph.b.w * TALLY_SCALE));
                    }
                }

                // Doc §13 (D7): speed and spin-edge-speed spread of
                // everything that escapes.
                let edge = length(ph.d.xyz) * exp2(f32(ph.c[3]));
                atomicAdd(&tallies[TALLY_ESC_ALL], 1);
                atomicAdd(&tallies[TALLY_ESC_SPEED_SUM], i32(speed * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_ESC_SPEED_SUMSQ], i32(speed * speed * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_ESC_EDGE_SUM], i32(edge * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_ESC_EDGE_SUMSQ], i32(edge * edge * TALLY_SCALE));
                if (edge > 1.0 + 1e-4) {
                    atomicAdd(&tallies[TALLY_ESC_EDGE_OVER_C], 1);
                }
                atomicAdd(&tallies[TALLY_ESC_LEVEL + min(ph.c[3], 3u)], 1);

                // E1 (doc §11d): mean stiffness level of escaping EMITTED
                // photons (plain integer sum -- k in 0..7, no scale needed).
                if (kind_emitted) {
                    atomicAdd(&tallies[TALLY_ESCAPED_EMITTED_K_SUM], i32(ph.c[3]));
                    atomicAdd(&tallies[TALLY_ESCAPED_EMITTED_K_COUNT], 1);
                }
            }
        }
    }

    if (dead) {
        ph.c[0] = ph.c[0] & ~FLAG_ALIVE;
        photons[idx] = ph;
        cell_of[idx] = params.num_cells;
        return;
    }

    ph.a = vec4<f32>(pos, ph.a.w);
    photons[idx] = ph;

    // Isotropy tallies (G1): fixed-point sum of dir.xyz and dz^2 over every
    // live photon, sampled once per step at move time (post-collide from the
    // previous step).
    atomicAdd(&tallies[63], i32(dir.x * TALLY_SCALE));
    atomicAdd(&tallies[64], i32(dir.y * TALLY_SCALE));
    atomicAdd(&tallies[65], i32(dir.z * TALLY_SCALE));
    atomicAdd(&tallies[66], i32(dir.z * dir.z * TALLY_SCALE));

    let gd = vec3<i32>(i32(params.grid_dim_x), i32(params.grid_dim_y), i32(params.grid_dim_z));
    let rel = (pos - vec3<f32>(params.grid_min_x, params.grid_min_y, params.grid_min_z)) / params.cell_h;
    let cx = clamp(i32(floor(rel.x)), 0, gd.x - 1);
    let cy = clamp(i32(floor(rel.y)), 0, gd.y - 1);
    let cz = clamp(i32(floor(rel.z)), 0, gd.z - 1);
    let cell_idx = u32(cx) + u32(cy) * params.grid_dim_x + u32(cz) * params.grid_dim_x * params.grid_dim_y;
    cell_of[idx] = cell_idx;

    // Density map (doc §8): time-averaged counts by kind x species, u32
    // per-cell accumulators added each step over the averaging window ONLY
    // (density_accumulate is 0 during warm-up). Never reset mid-window --
    // the host clears `density` once when the window starts and reads it
    // back once at the end.
    if (params.density_accumulate != 0u) {
        let kind_i2 = select(0u, 1u, kind_emitted);
        atomicAdd(&density[cell_idx * 4u + kind_i2 * 2u + spec], 1u);
    }
}
