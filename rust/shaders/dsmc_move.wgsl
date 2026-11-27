// DSMC step 2: advance every live photon by c*dt, test the nearest body for
// absorption/lane-crossing, apply the domain boundary, and emit this
// photon's cell index (dead -> bucket `num_cells + 1`, held photons ->
// parking bucket `num_cells`, which collide never visits) for CellSorter.
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
    body_beta: f32, body_beta_twist: f32, grain_mode: u32, maze_mode: u32,
    mix_c0: f32, mix_c1: f32, mix_c2: f32, mix_flags: u32,
    mix_mass_scale: f32, wind_x: f32, wind_y: f32, wind_z: f32,
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
const BALL_MAX_HITS: u32 = 6u;
const TALLY_PASS: u32 = 313u;         // [(abin*4 + outcome)*2 + spec], body 0, len 72
const TALLY_PASS_THROUGH: u32 = 385u; // [abin*2 + spec], body 0, len 18
const TALLY_PASS_HEL: u32 = 403u;
const TALLY_HIT_LAT: u32 = 412u;      // [bin], body 0, 1-degree bins, len 180
const TALLY_EXIT_LAT: u32 = 592u;     // [spec*180 + bin], body 0, len 360
const TALLY_MAZE: u32 = 1456u;     // [outcome], body 0, len 5 (4 = patch fallbacks)
const TALLY_MAZE_LAT: u32 = 1461u; // [bin], body 0, 1-degree bins, len 180
const TALLY_MAZE_AZ: u32 = 1641u;  // [bin], body 0, 10-degree bins, len 36
const TALLY_FLUX: u32 = 952u;         // [((k*18 + lat)*2 + spec)*2 + dir], body 0, len 504     // [outcome*2 + spec], body 0, len 8
const FLAG_PP_PASSAGE: u32 = 1u << 14u;
const FLAG_HIT_CLEAN: u32 = 1u << 15u;    // first ball hit came before any pp collision
const FLAG_HIT_AFTER_PP: u32 = 1u << 16u; // first ball hit came after a pp collision
// Doc §15.10: bits 17-19 = 1 + the body whose maze last emitted the photon
// (0 = never captured).
const FLAG_EMIT_SHIFT: u32 = 17u;
const FLAG_EMIT_MASK: u32 = 7u << 17u;
const TALLY_FEED: u32 = 1677u; // [capturing body*5 + emitter tag], len 20
const TALLY_LINE: u32 = 1697u; // [(plane*3 + tag class)*2 + dir], len 246
const LINE_U0: f32 = 2.0;
const LINE_DU: f32 = 0.25;
const LINE_N: i32 = 41;
// Doc §15.11 intake and holding tank.
const FLAG_HELD: u32 = 1u << 20u;  // held in a body's tank; c[1] = release step, c[2] & 7 = body + 1
const TALLY_HOLD: u32 = 1943u;     // [b*2 + (0 in, 1 out)], len 8
const TALLY_CAPT_LAT: u32 = 1951u; // [(tag != 0)*180 + bin], body 0 envelope captures only, len 360
const TALLY_ENV_CAPT: u32 = 2311u; // [b], len 4
const TALLY_HOLE_PASS: u32 = 2316u; // [b], len 4 (doc §15.24)
// Doc §15.40: photon spin angular momentum exchanged with body b. In:
// [(b*4 + cat)*3 + axis], cat = (entry north of the equator)*2 + (species
// plus), spin of photons taken in or deflected; out: [b*3 + axis], spin of
// the photons sent back out. Scaled by SPIN_SCALE.
const TALLY_SPIN_IN: u32 = 2320u;
const TALLY_SPIN_OUT: u32 = 2368u;
const SPIN_SCALE: f32 = 4194304.0;
// Photon radius in world units: level-1 over level-12 amplitude (1/512) of
// a body of radius 1 (baryon_shells::photon_radius).
const PHOTON_R_WORLD: f32 = 0.001953125;
// Doc §15.41: where each body's push comes from. Momentum and torque (about
// the body centre) brought in by the photons it takes in or deflects:
// [(b*20 + cat)*6 + (0..2 momentum, 3..5 torque)], cat = (emitter tag
// 0 never captured / 1 + last emitting body)*4 + (entry north of the
// equator)*2 + (hit a photon of another body's emission since it was last
// emitted); counts [b*20 + cat]; recoil of what it sends out [b*6 + ...].
const TALLY_SRC_IN: u32 = 2380u;
const TALLY_SRC_N: u32 = 2860u;
const TALLY_SRC_OUT: u32 = 2940u;
// Set by dsmc_collide when two photons last emitted by different bodies
// collide; cleared when a body emits the photon.
const FLAG_XHIT: u32 = 1u << 21u;
// Doc §15.48 beam-fate tally: bit 22 = released trapped (polar beam) by a
// neutron; bit 23 = has had a photon-photon collision since that release.
const FLAG_BEAM: u32 = 1u << 22u;
const FLAG_BEAM_COLL: u32 = 1u << 23u;
// Doc §15.52 pipe flags (also in dsmc_collide.wgsl).
const FLAG_PIPED: u32 = 1u << 24u;
const FLAG_LANE: u32 = 1u << 25u;
// [0..4) beam releases per body; [4..13) first collision region x partner
// (dsmc_collide); [13..23) rho of first collisions inside another body;
// [23..27) beam captured: other body clean / other collided / own clean /
// own collided. Must match gpu_types::TALLY_BEAM.
const TALLY_BEAM: u32 = 2989u;

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

const BALLKEY_STRIDE: u32 = 51u;
fn bk(b: u32, s: u32, j: u32) -> vec4<f32> { return ballkeys[b * BALLKEY_STRIDE + s * 4u + j]; }
// Doc §15.45 trace copies (gpu_types::BALLKEY_COPY_*): copy 0 is the real
// ball; copies 1.. are the same body at fixed time offsets along its trace.
const BALL_COPY_MAX: u32 = 64u;
fn n_copies(b: u32) -> u32 { return max(u32(ballkeys[4u * BALLKEY_STRIDE][b]), 1u); }
fn bkc(b: u32, m: u32, s: u32, j: u32) -> vec4<f32> {
    if (m == 0u) {
        return bk(b, s, j);
    }
    return ballkeys[4u * BALLKEY_STRIDE + 1u + (b * (BALL_COPY_MAX - 1u) + m - 1u) * 36u + s * 4u + j];
}
// Level-9 phase at the step start of copy m (for the swept maze frame).
fn copy_th9(b: u32, m: u32) -> f32 {
    if (m == 0u) {
        return ballkeys[b * BALLKEY_STRIDE + 36u].z;
    }
    return bkc(b, m, 0u, 1u).w;
}

fn quat_rot(q: vec4<f32>, v: vec3<f32>) -> vec3<f32> {
    let u = q.xyz;
    let s = q.w;
    return 2.0 * dot(u, v) * u + (s * s - dot(u, u)) * v + 2.0 * s * cross(u, v);
}

fn axis_rot(a: vec3<f32>, th: f32, v: vec3<f32>) -> vec3<f32> {
    let c = cos(th);
    let s = sin(th);
    return v * c + cross(a, v) * s + a * dot(a, v) * (1.0 - c);
}

// Doc §15.8: summed angular velocity (photon spin units) of the unresolved
// inner levels of ball b at random phases. Mirrors baryon_shells::inner_twist.
fn inner_twist(b: u32, q: vec4<f32>, key: u32) -> vec3<f32> {
    let m = u32(ballkeys[b * BALLKEY_STRIDE + 36u].x);
    var th: array<f32, 8>;
    for (var j = 0u; j < m; j = j + 1u) {
        th[j] = TAU * rand01(params.seed_hi, params.seed_lo, key, 40u + j);
    }
    var sum = vec3<f32>(0.0);
    for (var j = 0u; j < m; j = j + 1u) {
        let lv = ballkeys[b * BALLKEY_STRIDE + 37u + j];
        var w = lv.xyz * lv.w;
        for (var k = j + 1u; k < m; k = k + 1u) {
            let lk = ballkeys[b * BALLKEY_STRIDE + 37u + k];
            w = axis_rot(lk.xyz, th[k], w);
        }
        sum = sum + quat_rot(q, w);
    }
    return sum;
}

// Doc §15.9 maze emission law. Mirrors maze.rs (rot90, turn, spin_path, arc,
// in_patch, emit) and the maze branch of shell_gates::cpu_ball_step.
const MAZE_TRIES: u32 = 64u;
const PATCH_MARGIN: f32 = 1e-4;

fn maze_axis(i: u32) -> vec3<f32> {
    if (i == 1u) {
        return vec3<f32>(1.0, 0.0, 0.0);
    }
    if (i == 2u) {
        return vec3<f32>(0.0, 1.0, 0.0);
    }
    return vec3<f32>(0.0, 0.0, 1.0);
}

fn rot90(a: vec3<f32>, sg: f32, v: vec3<f32>) -> vec3<f32> {
    return a * dot(a, v) + sg * cross(a, v);
}

fn maze_turn(sg: vec4<f32>, v: vec3<f32>) -> vec3<f32> {
    var w = v;
    for (var i = 0u; i < 4u; i = i + 1u) {
        w = rot90(maze_axis(i), sg[i], w);
    }
    return w;
}

fn maze_spin(sg: vec4<f32>, s0: vec3<f32>) -> vec3<f32> {
    var s = s0;
    for (var i = 0u; i < 4u; i = i + 1u) {
        let a = maze_axis(i);
        s = rot90(a, sg[i], s);
        let n = dot(s, a);
        s = s + (sg[i] * abs(n) - n) * a;
    }
    return s;
}

fn arc_rot(a: vec3<f32>, b: vec3<f32>, v: vec3<f32>) -> vec3<f32> {
    let c = dot(a, b);
    if (c < -0.9999) {
        let t = select(vec3<f32>(0.0, 1.0, 0.0), vec3<f32>(1.0, 0.0, 0.0), abs(a.x) < 0.9);
        let ax = normalize(cross(a, t));
        return 2.0 * dot(ax, v) * ax - v;
    }
    let k = cross(a, b);
    return v * c + cross(k, v) + k * (dot(k, v) / (1.0 + c));
}

// Mathis's table over all 16 (a, x, y, z) states: escape iff x == z,
// upright iff y == z.
fn in_patch(d0: vec3<f32>, h: f32) -> bool {
    let s0 = d0 * h;
    for (var i = 0u; i < 16u; i = i + 1u) {
        let sg = vec4<f32>(
            select(-1.0, 1.0, (i & 1u) != 0u),
            select(-1.0, 1.0, (i & 2u) != 0u),
            select(-1.0, 1.0, (i & 4u) != 0u),
            select(-1.0, 1.0, (i & 8u) != 0u)
        );
        let e = dot(maze_turn(sg, d0), d0);
        let u = dot(maze_spin(sg, s0), s0);
        if (abs(e) < PATCH_MARGIN || abs(u) < PATCH_MARGIN || (e >= 0.0) != (sg.y == sg.w) || (u >= 0.0) != (sg.z == sg.w)) {
            return false;
        }
    }
    return true;
}

// World -> maze frame (body orient, then the swept level-9 phase th).
fn maze_local(oq: vec4<f32>, th: f32, v: vec3<f32>) -> vec3<f32> {
    return axis_rot(vec3<f32>(0.0, 0.0, 1.0), -th, quat_rot(vec4<f32>(-oq.xyz, oq.w), v));
}

fn maze_world(oq: vec4<f32>, th: f32, v: vec3<f32>) -> vec3<f32> {
    return quat_rot(oq, axis_rot(vec3<f32>(0.0, 0.0, 1.0), th, v));
}

// Doc §15.9: capture at the ball and re-emission from the core through the
// maze. The photon's helicity picks the patch hemisphere, a start d0 is drawn
// in the patch, its spin is carried rigidly from its arrival direction to
// d0, and it runs the canonical maze with the body's signs. Escapes leave
// along the maze exit, trapped photons out of the opposite pole. Speed is
// kept relative to the ball surface; inner stacked levels turn rigidly.
// Returns the emission point.
fn maze_contact(p: ptr<function, GpuPhoton>, b: u32, xh: vec3<f32>, ch: vec3<f32>, rb: f32, om: vec3<f32>, vcen: vec3<f32>, t_step: f32, key: u32, th0: f32) -> vec3<f32> {
    let base = b * BALLKEY_STRIDE;
    let bmeta = ballkeys[base + 36u];
    let oq = ballkeys[base + 45u];
    let sg = ballkeys[base + 46u];
    var th = 0.0;
    if ((params.maze_mode & 255u) == 2u) {
        th = th0 + bmeta.w * t_step;
    }
    // Doc §15.65/§15.66: bit 26 = averaged maze (default). Each capture turns the
    // maze by a random angle about the body axis, so its phase against the
    // trace is not locked (the maze is averaged around the axis).
    if (((params.maze_mode >> 26u) & 1u) != 0u) {
        th = th + TAU * rand01(params.seed_hi, params.seed_lo, key, 250u);
    }
    let dir = (*p).b.xyz;
    let v_in = dir * (*p).d.w;
    let dl = maze_local(oq, th, dir);
    let wl = maze_local(oq, th, (*p).d.xyz);
    let h = select(1.0, -1.0, dot(wl, dl) < 0.0);
    var d0 = normalize(vec3<f32>(1.0, 0.0, 1.0)) * -h;
    var fallback = true;
    // §15.11 control (mode 3): isotropic re-emission, no maze. Treated as a
    // maze whose start is the exit (d = d0, spin carried rigidly).
    let iso = (params.maze_mode & 255u) == 3u;
    if (iso) {
        let zz = 2.0 * rand01(params.seed_hi, params.seed_lo, key, 200u) - 1.0;
        let pz = TAU * rand01(params.seed_hi, params.seed_lo, key, 201u);
        let rz = sqrt(max(1.0 - zz * zz, 0.0));
        d0 = vec3<f32>(rz * cos(pz), rz * sin(pz), zz);
        fallback = false;
    }
    for (var t = 0u; t < select(MAZE_TRIES, 0u, iso); t = t + 1u) {
        let u = rand01(params.seed_hi, params.seed_lo, key, 50u + 2u * t);
        let ph = TAU * rand01(params.seed_hi, params.seed_lo, key, 51u + 2u * t);
        let r = sqrt(max(1.0 - u * u, 0.0));
        let c = vec3<f32>(r * cos(ph), r * sin(ph), -h * u);
        if (in_patch(c, h)) {
            d0 = c;
            fallback = false;
            break;
        }
    }
    let s0 = arc_rot(dl, d0, wl);
    var d = maze_turn(sg, d0);
    var s = maze_spin(sg, s0);
    if (iso) {
        d = d0;
        s = s0;
    }
    let escaped = dot(d, d0) >= 0.0;
    let upright = dot(s, s0) >= 0.0;
    // Doc §15.66 test knob: bit 27 = rigid spin. The maze picks the exit,
    // but the spin is carried rigidly from d0 to the exit (no imprint).
    if (!iso && ((params.maze_mode >> 27u) & 1u) != 0u) {
        s = arc_rot(d0, select(vec3<f32>(0.0, 0.0, h), d, escaped), s0);
    }
    // Doc §15.48: beam-fate tally. A beam photon reaching a ball is
    // captured here; then the flag is reset for this release.
    if (((*p).c[0] & FLAG_BEAM) != 0u) {
        let src = ((*p).c[0] >> 17u) & 7u;
        let own = src == b + 1u;
        let coll = ((*p).c[0] & FLAG_BEAM_COLL) != 0u;
        atomicAdd(&tallies[TALLY_BEAM + 23u + select(0u, 2u, own) + select(0u, 1u, coll)], 1);
    }
    let beam_now = !escaped && bodies[b].lane.z == 1.0;
    (*p).c[0] = ((*p).c[0] & ~(FLAG_BEAM | FLAG_BEAM_COLL)) | select(0u, FLAG_BEAM, beam_now);
    if (beam_now) {
        atomicAdd(&tallies[TALLY_BEAM + min(b, 3u)], 1);
    }
    var out_l = d;
    if (!escaped) {
        out_l = vec3<f32>(0.0, 0.0, h);
    }
    var out_w = maze_world(oq, th, out_l);
    // Doc §15.66/§15.67: bit 28 = escapes keep heading (default). An escaped photon
    // leaves in the hemisphere it was travelling toward (mirrored across the
    // equator if not), so photons in at the south leave north of the equator
    // (rhrule.pdf, magmom.pdf: emitted at 30 N / 30 S by species). Its spin is
    // mirrored the same way (helicity kept).
    var mirror_ax = vec3<f32>(0.0);
    if (escaped && ((params.maze_mode >> 28u) & 1u) != 0u) {
        let axw = bodies[b].axis.xyz;
        let zo = dot(out_w, axw);
        if (zo * dot(v_in, axw) < 0.0) {
            out_w = out_w - 2.0 * zo * axw;
            mirror_ax = axw;
        }
    }
    // Doc §15.32: bit 12 = through release. A trapped photon leaves by the
    // pole it was heading toward on entry (the neutron as a pipe: in one
    // side, out the other, in its own direction of travel).
    if (!escaped && ((params.maze_mode >> 12u) & 1u) != 0u) {
        let axw = bodies[b].axis.xyz;
        out_w = axw * select(-1.0, 1.0, dot(v_in, axw) >= 0.0);
    }
    let k = (*p).c[3];
    if (k > 0u) {
        var q = arc_rot(dl, d0, maze_local(oq, th, (*p).e.xyz));
        if (!iso) { q = maze_turn(sg, q); }
        (*p).e = vec4<f32>(maze_world(oq, th, q), (*p).e.w);
    }
    if (k > 1u) {
        var q = arc_rot(dl, d0, maze_local(oq, th, (*p).f.xyz));
        if (!iso) { q = maze_turn(sg, q); }
        (*p).f = vec4<f32>(maze_world(oq, th, q), (*p).f.w);
    }
    if (k > 2u) {
        var q = arc_rot(dl, d0, maze_local(oq, th, (*p).g.xyz));
        if (!iso) { q = maze_turn(sg, q); }
        (*p).g = vec4<f32>(maze_world(oq, th, q), (*p).g.w);
    }
    // Doc §15.14: bit 8 of maze_mode = MazeSpeed::EjectC, the baryon ejects
    // at c relative to its own surface (charge3.html); otherwise the photon
    // keeps its incoming speed relative to the surface (MazeSpeed::Relative).
    var s_rel = length(v_in - (vcen + cross(om, xh - ch)));
    if (((params.maze_mode >> 8u) & 1u) != 0u) {
        s_rel = 1.0;
    }
    var xe = ch + out_w * rb * (1.0 + 1e-3);
    var v_out = vcen + cross(om, xe - ch) + out_w * s_rel;
    // Doc §15.29: bit 10 = MazeSpeed::EjectSurface, at c from the envelope
    // (center.w = envelope radius) along the exit direction.
    if (((params.maze_mode >> 10u) & 1u) != 0u) {
        xe = bodies[b].center.xyz + out_w * bodies[b].center.w * (1.0 + 1e-3);
        v_out = out_w;
        if (!escaped) {
            // A trapped photon leaves along the pole. Spread its exit point
            // uniformly over the polar opening (radius 1/3 of the envelope,
            // the measured trace hole) instead of the single pole point,
            // which would make an infinitely thin beam.
            let ax = out_w;
            var e1 = normalize(cross(ax, vec3<f32>(1.0, 0.0, 0.0)));
            if (abs(ax.x) > 0.9) {
                e1 = normalize(cross(ax, vec3<f32>(0.0, 1.0, 0.0)));
            }
            let e2 = cross(ax, e1);
            // Doc §15.46: opening radius per body (lane.w; 0 = legacy 1/3).
            let beam = select(1.0 / 3.0, bodies[b].lane.w, bodies[b].lane.w > 0.0);
            let rr = beam * sqrt(rand01(params.seed_hi, params.seed_lo, key, 300u));
            let pa = TAU * rand01(params.seed_hi, params.seed_lo, key, 301u);
            let off = e1 * (rr * cos(pa)) + e2 * (rr * sin(pa));
            let R = bodies[b].center.w;
            let nrm = off + ax * sqrt(max(1.0 - rr * rr, 0.0));
            xe = bodies[b].center.xyz + nrm * R * (1.0 + 1e-3);
            // Doc §15.30: bit 11 = leave along the surface normal there (a
            // cone of half-angle asin(1/3) = 19.5 deg) instead of parallel
            // to the axis.
            if (((params.maze_mode >> 11u) & 1u) != 0u) {
                v_out = nrm;
            }
        }
    }
    // Doc §15.27: bit 9 = MazeSpeed::EjectRest, at c in the body frame.
    if (((params.maze_mode >> 9u) & 1u) != 0u) {
        v_out = out_w;
    }
    if (params.speed_mode == 1u) {
        v_out = normalize(v_out);
    }
    var w_out = maze_world(oq, th, s);
    w_out = w_out - 2.0 * dot(w_out, mirror_ax) * mirror_ax;
    let sp = length(v_out);
    let dn = v_out / max(sp, 1e-12);
    (*p).b = vec4<f32>(dn, axial_helicity(p, w_out, dn));
    (*p).d = vec4<f32>(w_out, sp);
    if (b == 0u) {
        atomicAdd(&tallies[TALLY_MAZE + select(2u, 0u, escaped) + select(1u, 0u, upright)], 1);
        if (fallback) {
            atomicAdd(&tallies[TALLY_MAZE + 4u], 1);
        }
        let lat = degrees(asin(clamp(dot(out_w, bodies[0].axis.xyz), -1.0, 1.0)));
        atomicAdd(&tallies[TALLY_MAZE_LAT + fine_lat_bin(lat)], 1);
        if (escaped) {
            let lo = quat_rot(vec4<f32>(-oq.xyz, oq.w), out_w);
            let az = degrees(atan2(lo.y, lo.x)) + 360.0;
            atomicAdd(&tallies[TALLY_MAZE_AZ + min(u32(floor(az - 360.0 * floor(az / 360.0)) / 10.0), 35u)], 1);
        }
    }
    return xe;
}

// Doc §15.11: body b takes the photon in at xh (ball contact point, or the
// core surface on its arrival side for an envelope capture). Tallies who
// emitted it last and where it was taken in; then either parks it in the
// tank (FLAG_HELD, released at step c[1]; returns the end-of-step core
// centre) or re-emits it at once through maze_contact (returns the emission
// point). `held` reports which.
fn capture(p: ptr<function, GpuPhoton>, b: u32, xh: vec3<f32>, ch: vec3<f32>, rb: f32, om: vec3<f32>, vcen: vec3<f32>, t_step: f32, key: u32, held: ptr<function, bool>, th0: f32) -> vec3<f32> {
    (*p).c[0] = (*p).c[0] & ~(FLAG_PIPED | FLAG_LANE);
    let tag = ((*p).c[0] >> FLAG_EMIT_SHIFT) & 7u;
    atomicAdd(&tallies[TALLY_FEED + b * 5u + tag], 1);
    let hold_steps = u32(ballkeys[b * BALLKEY_STRIDE + 47u].y);
    if (hold_steps > 0u) {
        (*p).c[0] = (*p).c[0] | FLAG_HELD;
        (*p).c[1] = params.step + hold_steps;
        (*p).c[2] = ((*p).c[2] & ~7u) | (b + 1u);
        atomicAdd(&tallies[TALLY_HOLD + b * 2u], 1);
        *held = true;
        return bk(b, 8u, 0u).xyz;
    }
    *held = false;
    let xe = maze_contact(p, b, xh, ch, rb, om, vcen, t_step, key, th0);
    (*p).c[0] = ((*p).c[0] & ~(FLAG_EMIT_MASK | FLAG_XHIT)) | ((b + 1u) << FLAG_EMIT_SHIFT);
    return xe;
}

fn cell_index(pos: vec3<f32>) -> u32 {
    let gd = vec3<i32>(i32(params.grid_dim_x), i32(params.grid_dim_y), i32(params.grid_dim_z));
    let rel = (pos - vec3<f32>(params.grid_min_x, params.grid_min_y, params.grid_min_z)) / params.cell_h;
    let cx = clamp(i32(floor(rel.x)), 0, gd.x - 1);
    let cy = clamp(i32(floor(rel.y)), 0, gd.y - 1);
    let cz = clamp(i32(floor(rel.z)), 0, gd.z - 1);
    return u32(cx) + u32(cy) * params.grid_dim_x + u32(cz) * params.grid_dim_x * params.grid_dim_y;
}

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

fn fine_lat_bin(lat_deg: f32) -> u32 {
    return u32(clamp(i32(floor(lat_deg + 90.0)), 0, 179));
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
    let m = photon_mass((*p).c[0], (*p).c[3]);
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
    // Body contacts are pinned only in PinnedC; BodyFree leaves them free.
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

// Doc §15.40: a photon's own spin angular momentum in world units:
// kappa m (r R_ph)^2 x angular velocity, with the spin vector w in photon
// spin units (edge speed r|w| in c), so L = kappa m r^2 R_ph w.
fn photon_spin_am(p: ptr<function, GpuPhoton>) -> vec3<f32> {
    let r = exp2(f32((*p).c[3]));
    let m = photon_mass((*p).c[0], (*p).c[3]);
    return params.rough_kappa * m * r * r * PHOTON_R_WORLD * (*p).d.xyz;
}

fn tally_spin(base: u32, l: vec3<f32>) {
    for (var a = 0u; a < 3u; a = a + 1u) {
        atomicAdd(&tallies[base + a], i32(l[a] * SPIN_SCALE));
    }
}

fn spin_cat(b: u32, x: vec3<f32>, p: ptr<function, GpuPhoton>) -> u32 {
    let north = dot(x - bodies[b].center.xyz, bodies[b].axis.xyz) >= 0.0;
    return select(0u, 2u, north) + select(0u, 1u, ((*p).c[0] & FLAG_SPECIES_PLUS) != 0u);
}

// Doc §15.41 source category of a photon reaching body b at x (read before
// the body re-tags it).
fn src_cat(b: u32, x: vec3<f32>, p: ptr<function, GpuPhoton>) -> u32 {
    let north = dot(x - bodies[b].center.xyz, bodies[b].axis.xyz) >= 0.0;
    let tag = min(((*p).c[0] >> FLAG_EMIT_SHIFT) & 7u, 4u);
    return tag * 4u + select(0u, 2u, north) + select(0u, 1u, ((*p).c[0] & FLAG_XHIT) != 0u);
}

fn tally_src_in(b: u32, cat: u32, mom: vec3<f32>, tq: vec3<f32>) {
    let base = TALLY_SRC_IN + (b * 20u + cat) * 6u;
    for (var a = 0u; a < 3u; a = a + 1u) {
        atomicAdd(&tallies[base + a], i32(mom[a] * TALLY_SCALE));
        atomicAdd(&tallies[base + 3u + a], i32(tq[a] * TALLY_SCALE));
    }
    atomicAdd(&tallies[TALLY_SRC_N + b * 20u + cat], 1);
}

fn tally_src_out(b: u32, mom: vec3<f32>, tq: vec3<f32>) {
    for (var a = 0u; a < 3u; a = a + 1u) {
        atomicAdd(&tallies[TALLY_SRC_OUT + b * 6u + a], i32(mom[a] * TALLY_SCALE));
        atomicAdd(&tallies[TALLY_SRC_OUT + b * 6u + 3u + a], i32(tq[a] * TALLY_SCALE));
    }
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
        cell_of[idx] = params.num_cells + 1u; // dead bucket (num_cells = held-photon parking, §15.11)
        return;
    }
    if ((idx & 63u) == 0u && (params.mix_flags != 0u || params.mix_c0 < 1.0)) {
        atomicAdd(&tallies[TALLY_LIVE_LEVEL + min(ph.c[3], 3u)], 1);
        atomicAdd(&tallies[TALLY_LIVE_BIRTH + ((ph.c[0] >> FLAG_BIRTH_SHIFT) & 3u)], 1);
    }

    // Doc §15.11 holding tank: a held photon rides the core until its
    // release step, then leaves through the maze at the step start and
    // moves on normally.
    if ((ph.c[0] & FLAG_HELD) != 0u) {
        let hb = (ph.c[2] & 7u) - 1u;
        if (params.step >= ph.c[1]) {
            let ch0 = bk(hb, 0u, 0u).xyz;
            let rb0 = bk(hb, 0u, 0u).w;
            let rkey = hash_u32(bitcast<u32>(ph.a.w) ^ hash_u32(params.step ^ 0x4E1Eu));
            ph.c[0] = ph.c[0] & ~FLAG_HELD;
            ph.c[1] = 0u;
            let xe = maze_contact(&ph, hb, ch0 - ph.b.xyz * rb0, ch0, rb0, bk(hb, 0u, 1u).xyz, bk(hb, 0u, 2u).xyz, 0.0, rkey, copy_th9(hb, 0u));
            ph.c[0] = (ph.c[0] & ~(FLAG_EMIT_MASK | FLAG_XHIT)) | ((hb + 1u) << FLAG_EMIT_SHIFT);
            ph.a = vec4<f32>(xe, ph.a.w);
            atomicAdd(&tallies[TALLY_HOLD + hb * 2u + 1u], 1);
        } else {
            let xc = bk(hb, 8u, 0u).xyz;
            ph.a = vec4<f32>(xc, ph.a.w);
            photons[idx] = ph;
            cell_of[idx] = params.num_cells;
            return;
        }
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

    // Doc §15.52 swirl: a piped photon that has not collided keeps to its
    // helix down the hole. Each step its velocity turns about the body axis
    // by (azimuthal speed x dt / helix radius); the walls take the change.
    if (params.body_model == 2u && (ph.c[0] & FLAG_PIPED) != 0u && (ph.c[2] & 7u) != 0u) {
        let bg = (ph.c[2] & 7u) - 1u;
        let hk = ballkeys[bg * BALLKEY_STRIDE + 50u];
        if (ballkeys[bg * BALLKEY_STRIDE + 49u].w != 0.0 && hk.x > 0.0) {
            // Doc §15.54 radius-held helix: the photon's own stacked-spin
            // trace. The velocity turns about the axis by |v_perp| dt / r_h
            // each step (r_h fixed), so the straight steps are chords of a
            // circle of radius r_h about a guide line parallel to the axis:
            // no radial drift. The body takes each turn (the fan reading,
            // rhrule.pdf), as at the mouth, so momentum stays consistent.
            let ax = bodies[bg].axis.xyz;
            let v_old = ph.b.xyz;
            let vp = length(v_old - dot(v_old, ax) * ax);
            let phi = sign(ballkeys[bg * BALLKEY_STRIDE + 49u].w) * vp * ph.d.w * params.dt / (hk.x * bodies[bg].center.w);
            let v_new = normalize(v_old * cos(phi) + cross(ax, v_old) * sin(phi) + ax * dot(ax, v_old) * (1.0 - cos(phi)));
            let m = photon_mass(ph.c[0], ph.c[3]);
            let dpg = -(v_new - v_old) * ph.d.w * m;
            let tqg = cross(ph.a.xyz - bodies[bg].center.xyz, dpg);
            for (var a = 0u; a < 3u; a = a + 1u) {
                atomicAdd(&tallies[TALLY_SH_DP + bg * 3u + a], i32(dpg[a] * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_SH_TORQUE + bg * 3u + a], i32(tqg[a] * TALLY_SCALE));
            }
            ph.b = vec4<f32>(v_new, ph.b.w);
        } else if (ballkeys[bg * BALLKEY_STRIDE + 49u].w != 0.0) {
            let ax = bodies[bg].axis.xyz;
            let rel = ph.a.xyz - bodies[bg].center.xyz;
            let perp = rel - dot(rel, ax) * ax;
            let pl = length(perp);
            if (pl > 1e-4) {
                let v_old = ph.b.xyz;
                let vaz = dot(v_old, cross(ax, perp / pl));
                let phi = vaz * ph.d.w * params.dt / pl;
                let v_new = v_old * cos(phi) + cross(ax, v_old) * sin(phi) + ax * dot(ax, v_old) * (1.0 - cos(phi));
                let m = photon_mass(ph.c[0], ph.c[3]);
                let dpg = -(v_new - v_old) * ph.d.w * m;
                let tqg = cross(rel, dpg);
                for (var a = 0u; a < 3u; a = a + 1u) {
                    atomicAdd(&tallies[TALLY_SH_DP + bg * 3u + a], i32(dpg[a] * TALLY_SCALE));
                    atomicAdd(&tallies[TALLY_SH_TORQUE + bg * 3u + a], i32(tqg[a] * TALLY_SCALE));
                }
                ph.b = vec4<f32>(v_new, ph.b.w);
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

    // Doc §15 ball: exact hit tests against each body's moving base ball over
    // 8 linear sub-intervals. Up to BALL_MAX_HITS contacts per step: after a
    // contact the photon continues from the contact point on its new velocity
    // and is tested again. A ball moving faster than c can catch the photon
    // again (snowplow). A photon found touching or inside the ball while the
    // ball still closes on it is pushed back to the surface and takes the
    // contact there (tallied as TALLY_BALL_OVERLAP). Mirrors
    // shell_gates::cpu_ball_step.
    var ball_moved = false;
    if (params.body_model == 2u) {
        let h = params.dt / 8.0;
        var x = pos;
        var tc = 0.0;
        for (var iter = 0u; iter < BALL_MAX_HITS; iter = iter + 1u) {
            let v = ph.b.xyz * ph.d.w;
            var best_t = 1e30;
            var best_b = 4u;
            var best_s = 0u;
            var best_m = 0u;
            var touching = false;
            // Doc §15.45: copy 0 is the real ball; copies 1.. are tested only
            // when the photon is near the body (the trace stays inside the
            // envelope).
            let reach = params.dt * ph.d.w + 0.25;
            for (var b = 0u; b < params.n_bodies; b = b + 1u) {
                let ncb = n_copies(b);
                let near_d = x - bodies[b].center.xyz;
                let near_r = bodies[b].center.w + reach;
                let near = dot(near_d, near_d) <= near_r * near_r;
                for (var m = 0u; m < ncb; m = m + 1u) {
                    if (m > 0u && !near) {
                        break;
                    }
                    let si0 = min(u32(tc / h), 7u);
                    let f0 = (tc - f32(si0) * h) / h;
                    let c_now = mix(bkc(b, m, si0, 0u).xyz, bkc(b, m, si0 + 1u, 0u).xyz, f0);
                    let rb = bkc(b, m, si0, 0u).w;
                    let vc_now = (bkc(b, m, si0 + 1u, 0u).xyz - bkc(b, m, si0, 0u).xyz) / h;
                    let d = x - c_now;
                    if (dot(d, d) <= rb * rb * (1.0 + 2e-4) && dot(d, v - vc_now) < 0.0) {
                        best_t = tc;
                        best_b = b;
                        best_s = si0;
                        best_m = m;
                        touching = true;
                        break;
                    }
                    for (var si = si0; si < 8u; si = si + 1u) {
                        let ts = max(f32(si) * h, tc);
                        let c0 = mix(bkc(b, m, si, 0u).xyz, bkc(b, m, si + 1u, 0u).xyz, (ts - f32(si) * h) / h);
                        let vc = (bkc(b, m, si + 1u, 0u).xyz - bkc(b, m, si, 0u).xyz) / h;
                        let t = moving_sphere_hit(x + v * (ts - tc), v, c0, vc, rb, f32(si + 1u) * h - ts);
                        var found = false;
                        if (t >= 0.0 && ts + t < best_t) {
                            best_t = ts + t;
                            best_b = b;
                            best_s = si;
                            best_m = m;
                            found = true;
                        }
                        if (found || ts >= best_t) {
                            break;
                        }
                    }
                }
                if (touching) {
                    break;
                }
            }
            if (best_b == 4u) {
                break;
            }
            let b = best_b;
            let si = best_s;
            let bm = best_m;
            let f = (best_t - f32(si) * h) / h;
            var xh = x + v * (best_t - tc);
            let ch = mix(bkc(b, bm, si, 0u).xyz, bkc(b, bm, si + 1u, 0u).xyz, f);
            let rb = bkc(b, bm, si, 0u).w;
            if (touching) {
                let dd = xh - ch;
                if (dot(dd, dd) > 1e-20) {
                    xh = ch + normalize(dd) * rb;
                }
                atomicAdd(&tallies[TALLY_BALL_OVERLAP], 1);
            }
            let om = mix(bkc(b, bm, si, 1u).xyz, bkc(b, bm, si + 1u, 1u).xyz, f);
            let vcen = mix(bkc(b, bm, si, 2u).xyz, bkc(b, bm, si + 1u, 2u).xyz, f);
            let pp_before = (ph.c[0] & FLAG_PP_PASSAGE) != 0u;
            var dp = vec3<f32>(0.0);
            var touched = false;
            var x_next = xh;
            var now_held = false;
            // Doc §15.20: a maze capture's torque about the body centre uses
            // the photon's true lines of action, absorbed at the contact xh
            // and re-emitted from x_next (the ball surface opposite).
            var tq_maze = vec3<f32>(0.0);
            var has_tq_maze = false;
            if ((params.maze_mode & 255u) != 0u && ballkeys[b * BALLKEY_STRIDE + 47u].x > 0.5) {
                // §15.9/§15.11: capture, then the tank or the maze.
                let v_before = ph.b.xyz * ph.d.w;
                let mkey = hash_u32(bitcast<u32>(ph.a.w) ^ hash_u32(params.step ^ 0x3A2Eu ^ (iter * 0x9E3779B9u)));
                tally_spin(TALLY_SPIN_IN + (b * 4u + spin_cat(b, xh, &ph)) * 3u, photon_spin_am(&ph));
                let scat = src_cat(b, xh, &ph);
                x_next = capture(&ph, b, xh, ch, rb, om, vcen, best_t, mkey, &now_held, copy_th9(b, bm));
                if (!now_held) {
                    tally_spin(TALLY_SPIN_OUT + b * 3u, photon_spin_am(&ph));
                }
                let m = photon_mass(ph.c[0], ph.c[3]);
                let bcm = bodies[b].center.xyz;
                tally_src_in(b, scat, v_before * m, cross(xh - bcm, v_before * m));
                if (now_held) {
                    dp = v_before * m;
                    tq_maze = cross(xh - bcm, v_before * m);
                } else {
                    dp = -(ph.b.xyz * ph.d.w - v_before) * m;
                    tq_maze = cross(xh - bcm, v_before * m) - cross(x_next - bcm, ph.b.xyz * ph.d.w * m);
                    tally_src_out(b, -ph.b.xyz * ph.d.w * m, -cross(x_next - bcm, ph.b.xyz * ph.d.w * m));
                }
                has_tq_maze = true;
                touched = true;
            } else {
                let n = normalize(ch - xh);
                var tw = vec3<f32>(0.0);
                if (params.grain_mode == 0u) {
                    if (length(om) > 1e-9) {
                        tw = normalize(om) * params.grain_edge;
                    }
                } else {
                    // §15.8: the ball's rigid spin in photon spin units, plus
                    // (mode 2) the unresolved inner levels at random phases.
                    tw = om * ballkeys[b * BALLKEY_STRIDE + 36u].y;
                    if (params.grain_mode == 2u) {
                        let q = normalize(mix(bkc(b, bm, si, 3u), bkc(b, bm, si + 1u, 3u), f));
                        let key = hash_u32(bitcast<u32>(ph.a.w) ^ hash_u32(params.step ^ 0x7A157u ^ (iter * 0x9E3779B9u)));
                        tw = tw + inner_twist(b, q, key);
                    }
                }
                dp = wall_contact_general(&ph, vcen, rb * om, tw, n);
                touched = dot(dp, dp) > 0.0;
            }
            if (touched) {
                let bc = bodies[b].center.xyz;
                let tq = select(cross(xh - bc, dp), tq_maze, has_tq_maze);
                for (var a = 0u; a < 3u; a = a + 1u) {
                    atomicAdd(&tallies[TALLY_SH_DP + b * 3u + a], i32(dp[a] * TALLY_SCALE));
                    atomicAdd(&tallies[TALLY_SH_TORQUE + b * 3u + a], i32(tq[a] * TALLY_SCALE));
                }
                atomicAdd(&tallies[TALLY_SH_CONTACTS + b], 1);
                if (b == 0u) {
                    atomicAdd(&tallies[TALLY_HIT_LAT + fine_lat_bin(body_lat_deg(0u, xh))], 1);
                }
                if ((ph.c[2] >> 13u) == 0u && (ph.c[0] & (FLAG_HIT_CLEAN | FLAG_HIT_AFTER_PP)) == 0u) {
                    ph.c[0] = ph.c[0] | select(FLAG_HIT_CLEAN, FLAG_HIT_AFTER_PP, pp_before);
                }
                ph.c[2] = ph.c[2] + (1u << 13u);
            }
            x = x_next;
            tc = best_t;
            ball_moved = true;
            if (now_held) {
                tc = params.dt;
                break;
            }
        }
        if (ball_moved) {
            pos = x + (ph.b.xyz * ph.d.w) * (params.dt - tc);
        }
        if ((ph.c[0] & FLAG_HELD) != 0u) {
            ph.a = vec4<f32>(pos, ph.a.w);
            photons[idx] = ph;
            cell_of[idx] = params.num_cells;
            return;
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

    // Doc §15.6 charge-pause map: count this step's crossings of the body-0
    // shells r = 0.5 .. 5 (gpu_types FLUX_RADII), by the latitude at the end
    // of the step (a step is 0.05 long, small next to the shell spacing).
    if (params.body_model == 2u && params.n_bodies > 0u) {
        let c0 = bodies[0].center.xyz;
        let s0 = length(ph.a.xyz - c0);
        let s1 = length(pos - c0);
        let lb = u32(clamp(i32(floor((body_lat_deg(0u, pos) + 90.0) / 10.0)), 0, 17));
        let radii = array<f32, 7>(0.5, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0);
        for (var k = 0u; k < 7u; k = k + 1u) {
            let r = radii[k];
            if ((s0 < r) != (s1 < r)) {
                let dir = select(1u, 0u, s1 > s0);
                atomicAdd(&tallies[TALLY_FLUX + ((k * 18u + lb) * 2u + spec) * 2u + dir], 1);
            }
        }
    }

    // Doc §15.10 two-body line flux: crossings of the planes perpendicular to
    // the body-0 -> body-1 axis at u = -LINE_U0 + j*LINE_DU (u from body 0's
    // centre), inside a cylinder of radius 0.5 about the axis, by the
    // photon's emitter tag (0 never captured, 1 body 0, 2 body 1) and
    // direction (0 toward body 1).
    if (params.body_model == 2u && params.n_bodies > 1u) {
        let c0 = bodies[0].center.xyz;
        let e = normalize(bodies[1].center.xyz - c0);
        let a0 = ph.a.xyz - c0;
        let a1 = pos - c0;
        let u0 = dot(a0, e);
        let u1 = dot(a1, e);
        if (u0 != u1) {
            let jlo = max(i32(ceil((min(u0, u1) + LINE_U0) / LINE_DU)), 0);
            let jhi = min(i32(floor((max(u0, u1) + LINE_U0) / LINE_DU)), LINE_N - 1);
            let tc = min((ph.c[0] >> FLAG_EMIT_SHIFT) & 7u, 2u);
            let dir = select(1u, 0u, u1 > u0);
            for (var j = jlo; j <= jhi; j = j + 1) {
                let up = f32(j) * LINE_DU - LINE_U0;
                let xc = mix(a0, a1, (up - u0) / (u1 - u0));
                let perp = xc - up * e;
                if (dot(perp, perp) <= 0.25) {
                    atomicAdd(&tallies[TALLY_LINE + (u32(j) * 3u + tc) * 2u + dir], 1);
                }
            }
        }
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
                    if (b == 0u) {
                        atomicAdd(&tallies[TALLY_EXIT_LAT + spec * 180u + fine_lat_bin(lat_out)], 1);
                    }
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
            ph.c[0] = ph.c[0] & ~(FLAG_PP_PASSAGE | FLAG_HIT_CLEAN | FLAG_HIT_AFTER_PP | FLAG_PIPED | FLAG_LANE);
            // Doc §15.51 pipe test: a photon entering the polar mouth heading
            // inward within the pipe angle of the axis is turned onto the axis
            // (down the trace's clear hole). The body takes the sideways
            // momentum of the turn, attributed to the photon's source.
            if (params.body_model == 2u && now_inside != 0u && prev_inside == 0u) {
                let b = now_inside - 1u;
                let pk = ballkeys[b * BALLKEY_STRIDE + 49u];
                if (pk.x > 0.5) {
                    let ce = bodies[b].center.xyz;
                    let rr = bodies[b].center.w;
                    let ax = bodies[b].axis.xyz;
                    let rel = pos - ce;
                    let za = dot(rel, ax);
                    let perp = rel - za * ax;
                    let ca = dot(ph.b.xyz, ax);
                    if (ca * za < 0.0 && abs(ca) >= pk.y && dot(perp, perp) <= pk.z * pk.z * rr * rr) {
                        let m = photon_mass(ph.c[0], ph.c[3]);
                        var nd = ax * sign(ca);
                        // Doc §15.52 swirl: a helix about the axis, turning
                        // with the body (pk.w = signed helix angle).
                        let pl = length(perp);
                        if (pk.w != 0.0 && pl > 1e-4) {
                            nd = nd * cos(pk.w) + cross(ax, perp / pl) * sin(pk.w);
                        }
                        // Doc §15.70 funnel test: aim at the axis point
                        // [50].z radii past the centre (the vortex cone).
                        let fk = ballkeys[b * BALLKEY_STRIDE + 50u];
                        if (pk.w == 0.0 && fk.y > 0.5) {
                            let dv = ce + ax * sign(ca) * fk.z * rr - pos;
                            if (dot(dv, dv) > 1e-8) { nd = normalize(dv); }
                        }
                        let dpp = (ph.b.xyz - nd) * ph.d.w * m;
                        let tqp = cross(rel, dpp);
                        tally_src_in(b, src_cat(b, pos, &ph), dpp, tqp);
                        for (var a = 0u; a < 3u; a = a + 1u) {
                            atomicAdd(&tallies[TALLY_SH_DP + b * 3u + a], i32(dpp[a] * TALLY_SCALE));
                            atomicAdd(&tallies[TALLY_SH_TORQUE + b * 3u + a], i32(tqp[a] * TALLY_SCALE));
                        }
                        ph.b = vec4<f32>(nd, axial_helicity(&ph, ph.d.xyz, nd));
                        ph.c[0] = ph.c[0] | FLAG_PIPED | select(0u, FLAG_LANE, pk.x > 1.5);
                        atomicAdd(&tallies[TALLY_HOLE_PASS + b], 1);
                    }
                }
            }
            // Doc §15.11: envelope intake. A photon entering body b's
            // envelope is taken in (envelope mode: with probability p
            // anywhere; polar mode: within `aperture` of the axis) and
            // funnelled to the core.
            if (params.body_model == 2u && (params.maze_mode & 255u) != 0u && now_inside != 0u && prev_inside == 0u) {
                let b = now_inside - 1u;
                let ik = ballkeys[b * BALLKEY_STRIDE + 47u];
                var take = false;
                let ekey = hash_u32(bitcast<u32>(ph.a.w) ^ hash_u32(params.step ^ 0xE17Eu));
                let lk = ballkeys[b * BALLKEY_STRIDE + 48u];
                if (lk.y > 0.5) {
                    // Doc §15.32: a recycling electron. Every photon that
                    // reaches it is taken in and re-emitted at once from its
                    // surface, at c, into its equatorial disc (half-width
                    // lk.z radians about the spin axis), random azimuth.
                    let v_before = ph.b.xyz * ph.d.w;
                    let pos_in = pos;
                    let ce = bodies[b].center.xyz;
                    let ax = bodies[b].axis.xyz;
                    var e1 = normalize(cross(ax, vec3<f32>(1.0, 0.0, 0.0)));
                    if (abs(ax.x) > 0.9) {
                        e1 = normalize(cross(ax, vec3<f32>(0.0, 1.0, 0.0)));
                    }
                    let e2 = cross(ax, e1);
                    let phe = TAU * rand01(params.seed_hi, params.seed_lo, ekey, 400u);
                    let late = (2.0 * rand01(params.seed_hi, params.seed_lo, ekey, 401u) - 1.0) * lk.z;
                    let dir = (e1 * cos(phe) + e2 * sin(phe)) * cos(late) + ax * sin(late);
                    pos = ce + dir * bodies[b].center.w * (1.0 + 1e-3);
                    ph.b = vec4<f32>(dir, axial_helicity(&ph, ph.d.xyz, dir));
                    ph.d = vec4<f32>(ph.d.xyz, 1.0);
                    ph.c[0] = (ph.c[0] & ~(FLAG_EMIT_MASK | FLAG_XHIT)) | ((b + 1u) << FLAG_EMIT_SHIFT);
                    let m = photon_mass(ph.c[0], ph.c[3]);
                    let v_after = dir;
                    let dp = -(v_after - v_before) * m;
                    let tqe = cross(pos_in - ce, v_before * m) - cross(pos - ce, v_after * m);
                    for (var a = 0u; a < 3u; a = a + 1u) {
                        atomicAdd(&tallies[TALLY_SH_DP + b * 3u + a], i32(dp[a] * TALLY_SCALE));
                        atomicAdd(&tallies[TALLY_SH_TORQUE + b * 3u + a], i32(tqe[a] * TALLY_SCALE));
                    }
                    atomicAdd(&tallies[TALLY_ENV_CAPT + b], 1);
                    atomicAdd(&tallies[TALLY_SH_CONTACTS + b], 1);
                } else if (ik.x > 0.5 && ik.w > 0.5 && ik.w < 1.5) {
                    take = rand01(params.seed_hi, params.seed_lo, ekey, 0u) < ik.z;
                    // Doc §15.24: the through-charge hole passes polar
                    // entries straight through.
                    let hole = ballkeys[b * BALLKEY_STRIDE + 48u].x;
                    if (hole > 0.0) {
                        let relh = pos - bodies[b].center.xyz;
                        let axh = bodies[b].axis.xyz;
                        let perph = relh - dot(relh, axh) * axh;
                        if (dot(perph, perph) <= hole * hole) {
                            take = false;
                            atomicAdd(&tallies[TALLY_HOLE_PASS + b], 1);
                        }
                    }
                } else if (ik.x > 0.5 && ik.w > 1.5 && ik.w < 2.5) {
                    let rel = pos - bodies[b].center.xyz;
                    let ax = bodies[b].axis.xyz;
                    let perp = rel - dot(rel, ax) * ax;
                    take = dot(perp, perp) <= ik.z * ik.z;
                } else if (ik.w > 2.5) {
                    // Doc §15.35 PolarDeflect: in through the polar opening;
                    // anywhere else the photon is batted off the outer spin,
                    // a sphere of the envelope radius turning about the axis
                    // with signed equatorial surface speed lk.w (rough
                    // contact law, as for the core ball), and stays outside.
                    let ce = bodies[b].center.xyz;
                    let rr = bodies[b].center.w;
                    let rel = pos - ce;
                    let ax = bodies[b].axis.xyz;
                    let perp = rel - dot(rel, ax) * ax;
                    if (ik.x > 0.5 && dot(perp, perp) <= ik.z * ik.z * rr * rr) {
                        take = true;
                    } else {
                        let n = -normalize(rel);
                        let om_env = ax * (lk.w / rr);
                        let tw = om_env * ballkeys[b * BALLKEY_STRIDE + 36u].y;
                        tally_spin(TALLY_SPIN_IN + (b * 4u + spin_cat(b, pos, &ph)) * 3u, photon_spin_am(&ph));
                        let scat_d = src_cat(b, pos, &ph);
                        let dpd = wall_contact_general(&ph, vec3<f32>(0.0), om_env * rr, tw, n);
                        tally_spin(TALLY_SPIN_OUT + b * 3u, photon_spin_am(&ph));
                        pos = ce + normalize(rel) * rr * (1.0 + 1e-3);
                        ph.c[2] = 0u;
                        let tqd = cross(pos - ce, dpd);
                        tally_src_in(b, scat_d, dpd, tqd);
                        for (var a = 0u; a < 3u; a = a + 1u) {
                            atomicAdd(&tallies[TALLY_SH_DP + b * 3u + a], i32(dpd[a] * TALLY_SCALE));
                            atomicAdd(&tallies[TALLY_SH_TORQUE + b * 3u + a], i32(tqd[a] * TALLY_SCALE));
                        }
                        atomicAdd(&tallies[TALLY_HOLE_PASS + b], 1);
                    }
                }
                if (take) {
                    let lat_in = body_lat_deg(b, pos);
                    tally_spin(TALLY_SPIN_IN + (b * 4u + spin_cat(b, pos, &ph)) * 3u, photon_spin_am(&ph));
                    let ch8 = bk(b, 8u, 0u).xyz;
                    let rb8 = bk(b, 8u, 0u).w;
                    var eheld = false;
                    let v_before = ph.b.xyz * ph.d.w;
                    let pos_in = pos;
                    // Tally the capture at the envelope entry latitude.
                    let xh8 = ch8 - ph.b.xyz * rb8;
                    let tag = (ph.c[0] >> FLAG_EMIT_SHIFT) & 7u;
                    let scat = src_cat(b, pos, &ph);
                    if (b == 0u) {
                        atomicAdd(&tallies[TALLY_CAPT_LAT + select(0u, 180u, tag != 0u) + fine_lat_bin(lat_in)], 1);
                        atomicAdd(&tallies[TALLY_HIT_LAT + fine_lat_bin(lat_in)], 1);
                    }
                    atomicAdd(&tallies[TALLY_FEED + b * 5u + tag], 1);
                    atomicAdd(&tallies[TALLY_ENV_CAPT + b], 1);
                    atomicAdd(&tallies[TALLY_SH_CONTACTS + b], 1);
                    let hold_steps = u32(ik.y);
                    if (hold_steps > 0u) {
                        ph.c[0] = ph.c[0] | FLAG_HELD;
                        ph.c[1] = params.step + hold_steps;
                        atomicAdd(&tallies[TALLY_HOLD + b * 2u], 1);
                        pos = ch8;
                    } else {
                        pos = maze_contact(&ph, b, xh8, ch8, rb8, bk(b, 8u, 1u).xyz, bk(b, 8u, 2u).xyz, params.dt, ekey, copy_th9(b, 0u));
                        tally_spin(TALLY_SPIN_OUT + b * 3u, photon_spin_am(&ph));
                        ph.c[0] = (ph.c[0] & ~(FLAG_EMIT_MASK | FLAG_XHIT)) | ((b + 1u) << FLAG_EMIT_SHIFT);
                    }
                    let m = photon_mass(ph.c[0], ph.c[3]);
                    let v_after = select(ph.b.xyz * ph.d.w, vec3<f32>(0.0), (ph.c[0] & FLAG_HELD) != 0u);
                    let dp = -(v_after - v_before) * m;
                    // Doc §15.20: torque about the body centre, absorbed on
                    // the photon's incoming line (through pos_in) and
                    // re-emitted from the core surface (pos).
                    let bce = bodies[b].center.xyz;
                    let tqe = cross(pos_in - bce, v_before * m) - cross(pos - bce, v_after * m);
                    tally_src_in(b, scat, v_before * m, cross(pos_in - bce, v_before * m));
                    tally_src_out(b, -v_after * m, -cross(pos - bce, v_after * m));
                    for (var a = 0u; a < 3u; a = a + 1u) {
                        atomicAdd(&tallies[TALLY_SH_DP + b * 3u + a], i32(dp[a] * TALLY_SCALE));
                        atomicAdd(&tallies[TALLY_SH_TORQUE + b * 3u + a], i32(tqe[a] * TALLY_SCALE));
                    }
                    ph.c[2] = ph.c[2] + (1u << 13u);
                }
            }
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
        cell_of[idx] = params.num_cells + 1u; // dead bucket (num_cells = held-photon parking, §15.11)
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

    let cell_idx = cell_index(pos);
    if ((ph.c[0] & FLAG_HELD) != 0u) {
        // §15.11: captured into the tank at the envelope this step.
        cell_of[idx] = params.num_cells;
        return;
    }
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
