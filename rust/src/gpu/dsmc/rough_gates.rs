//! Doc §13 gates for `ExchangeRule::RoughSphere` (GPU side). The CPU half
//! (per-collision conservation, the RS3 sign table) lives in
//! `crate::rough_sphere`.
//!
//! - `rs_gpu_matches_cpu_reference`: many isolated two-photon cells, one step
//!   each. Every accepted GPU collision must equal `rough_sphere::collide_full`
//!   (+ `pin_speed`) fed the same contact normal. The kernel's counter-based
//!   RNG is replicated here bit for bit to get that normal. Negative control:
//!   the same comparison against a wrong law (beta = 0.5) must fail.
//! - `rs1_gpu_box_conservation`: periodic box, mixed levels. Free speed
//!   conserves energy and momentum (f64 readback); beta = 0 must lose energy;
//!   PinnedC keeps every |v| = c and conserves energy but not momentum.
//! - `rs2_gpu_equipartition`: a perfectly rough gas relaxes to
//!   <E_rot> = <E_lin> with a Maxwellian speed distribution
//!   (<v^4>/<v^2>^2 = 5/3). Smooth spheres (beta = -1) must stay at the
//!   birth ratio.

use glam::DVec3;

use super::config::{Boundary, DsmcConfig, EmissionMode, ExchangeRule, MeshMode, RateKernel, RoughParams, SpeedMode, StiffnessShed};
use super::engine::{DsmcEngine, Totals};
use super::gpu_types::{GpuPhoton, FLAG_ALIVE, FLAG_SPECIES_PLUS};
use crate::gpu::context::GpuContext;
use crate::rough_sphere::{collide_full, contact_normal, pin_pair, Ball, ContactLaw};

fn hash_u32(x: u32) -> u32 {
    let mut v = x;
    v ^= v >> 16;
    v = v.wrapping_mul(0x7feb_352d);
    v ^= v >> 15;
    v = v.wrapping_mul(0x846c_a68b);
    v ^= v >> 16;
    v
}

/// Bit-exact replica of dsmc_common.wgsl `rand01`.
fn rand01(seed: u64, id: u32, draw: u32) -> f32 {
    let (seed_hi, seed_lo) = ((seed >> 32) as u32, seed as u32);
    let mut h = hash_u32(seed_hi ^ 0x9e37_79b9);
    h = hash_u32(h ^ seed_lo);
    h = hash_u32(h ^ id);
    h = hash_u32(h ^ draw);
    (h as f32) * (1.0f32 / 4_294_967_296.0)
}

struct Lcg(u64);
impl Lcg {
    fn f(&mut self) -> f64 {
        self.0 = self.0.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64
    }
    fn unit(&mut self) -> DVec3 {
        let z = 2.0 * self.f() - 1.0;
        let p = std::f64::consts::TAU * self.f();
        let r = (1.0 - z * z).max(0.0).sqrt();
        DVec3::new(r * p.cos(), r * p.sin(), z)
    }
}

fn box_config(half_extent: f64, h: f64, dt: f64, sigma: f64, capacity: u32, seed: u64, rough: RoughParams) -> DsmcConfig {
    DsmcConfig {
        boundary: Boundary::PeriodicBox { half_extent },
        bodies: vec![],
        h,
        dt,
        capacity,
        seed,
        n_ambient: 0.0,
        sigma,
        photon_fraction: 0.5,
        pole_intake_bias: 0.0,
        cancel_redirect: 0.0,
        spin_transfer: 0.0,
        swing_boost: 0.0,
        rate_kernel: RateKernel::Catch,
        mesh_mode: MeshMode::ChiralityOnly,
        exchange_rule: ExchangeRule::RoughSphere,
        emission_mode: EmissionMode::Forced,
        bounce_cap: 0,
        d5_far_zone_gate: false,
        b_max: 1.0,
        emit_stiffness_level: 0,
        stiffness_shed: StiffnessShed::None,
        sigma_size_scaling: false,
        debug_force_equal_mass_exchange: false,
        rough,
    }
}

fn make_photon(id: u32, pos: DVec3, v: DVec3, w: DVec3, k: u32, chi: f64) -> GpuPhoton {
    let s = v.length();
    let d = v / s;
    let r = (2.0f64).powi(k as i32);
    let mut flags = FLAG_ALIVE;
    if chi > 0.0 {
        flags |= FLAG_SPECIES_PLUS;
    }
    GpuPhoton {
        a: [pos.x as f32, pos.y as f32, pos.z as f32, f32::from_bits(id)],
        b: [d.x as f32, d.y as f32, d.z as f32, (w.dot(d) * r) as f32],
        c: [flags, 0, 0, k],
        d: [w.x as f32, w.y as f32, w.z as f32, s as f32],
    }
}

fn to_ball(p: &GpuPhoton) -> Ball {
    let v = p.velocity();
    let w = p.omega();
    Ball::level(p.c[3], DVec3::new(v[0], v[1], v[2]), DVec3::new(w[0], w[1], w[2]))
}

fn law_of(r: &RoughParams) -> ContactLaw {
    ContactLaw { e: 1.0, beta: r.beta, beta_twist: r.beta_twist, twist_loss_to_normal: r.twist_loss_to_normal, kappa: r.kappa }
}

struct MatchResult {
    accepted: usize,
    checked: usize,
    max_err: f64,
    max_err_wrong_law: f64,
    unaccepted_changed: usize,
}

fn run_pair_match(rough: RoughParams, seed: u64) -> MatchResult {
    const HALF: f64 = 5.0;
    const H: f64 = 1.0;
    const DT: f64 = 0.01;
    // One NTC candidate per two-photon cell: m_f = sigma*g_max*dt/V = 1.5.
    let sigma = 1.5 / (rough.g_max * DT);
    let config = box_config(HALF, H, DT, sigma, 4096, seed, rough);
    let ctx = GpuContext::new().expect("RS gates need a GPU adapter");
    let mut engine = DsmcEngine::new(&ctx, config);
    let (dim, min) = (engine.grid.dim, engine.grid.min);
    let free_speed = rough.speed_mode == SpeedMode::Free;

    let mut rng = Lcg(seed ^ 0xABCD);
    let mut photons = Vec::new();
    let mut cells = Vec::new();
    let mut id = 0u32;
    for cz in (0..dim[2]).step_by(2) {
        for cy in (0..dim[1]).step_by(2) {
            for cx in 0..dim[0] {
                let c = engine.grid.cell_center(cx, cy, cz);
                let cell = cx + cy * dim[0] + cz * dim[0] * dim[1];
                for off in [-0.2, 0.2] {
                    let k = (rng.f() * 3.0) as u32;
                    let speed = if free_speed { 0.5 + rng.f() } else { 1.0 };
                    let v = rng.unit() * speed;
                    let w = rng.unit() * (1.5 * rng.f()) / (2.0f64).powi(k as i32);
                    let chi = if rng.f() < 0.5 { 1.0 } else { -1.0 };
                    photons.push(make_photon(id, c + DVec3::splat(off), v, w, k, chi));
                    id += 1;
                }
                cells.push(cell);
            }
        }
    }
    let _ = min;
    let before = photons.clone();
    engine.write_photons_for_test(&ctx, &photons);
    let step = engine.step;
    let mut totals = Totals::default();
    engine.step_collide(&ctx, &mut totals);
    let after = engine.read_photons(&ctx);
    let by_id = |id: u32| after.iter().find(|p| p.id() == id).copied().expect("photon lost");

    let law = law_of(&rough);
    let wrong = ContactLaw { beta: 0.5, ..law };
    let mut res = MatchResult { accepted: 0, checked: 0, max_err: 0.0, max_err_wrong_law: 0.0, unaccepted_changed: 0 };
    for (pair, &cell) in cells.iter().enumerate() {
        let (p0, p1) = (before[2 * pair], before[2 * pair + 1]);
        let key_base = hash_u32(cell ^ hash_u32(step ^ 0x9E37_79B9));
        let ui = rand01(seed, key_base, 0);
        let i = ((ui * 2.0) as u32).min(1);
        let (pa, pb) = if i == 0 { (p0, p1) } else { (p1, p0) };
        let (qa, qb) = (by_id(pa.id()), by_id(pb.id()));
        // Acceptance, in the kernel's f32 arithmetic.
        let va32 = [pa.b[0] * pa.d[3], pa.b[1] * pa.d[3], pa.b[2] * pa.d[3]];
        let vb32 = [pb.b[0] * pb.d[3], pb.b[1] * pb.d[3], pb.b[2] * pb.d[3]];
        let g32 = ((va32[0] - vb32[0]).powi(2) + (va32[1] - vb32[1]).powi(2) + (va32[2] - vb32[2]).powi(2)).sqrt();
        let p_acc = g32 / rough.g_max as f32;
        let u_acc = rand01(seed, key_base, 2);
        if (u_acc - p_acc).abs() < 1e-4 {
            continue; // too close to call across f32 rounding
        }
        if u_acc >= p_acc {
            let unchanged = (0..4).all(|c| (qa.d[c] - pa.d[c]).abs() < 1e-7 && (qb.d[c] - pb.d[c]).abs() < 1e-7);
            if !unchanged {
                res.unaccepted_changed += 1;
            }
            continue;
        }
        res.accepted += 1;
        let rho = rand01(seed, key_base, 3) as f64;
        let psi = rand01(seed, key_base, 4) as f64 * std::f64::consts::TAU;
        let (ba0, bb0) = (to_ball(&pa), to_ball(&pb));
        let u = ba0.v - bb0.v;
        let n = contact_normal(u, rho, psi);
        let slip_n = crate::rough_sphere::contact_slip(&ba0, &bb0, n).dot(n);
        if slip_n.abs() < 1e-4 {
            continue; // grazing: separating-or-not is decided by f32 rounding
        }
        let expect = |l: &ContactLaw| {
            let (mut a, mut b) = (ba0, bb0);
            let e_rot_in = a.e_rot(l.kappa) + b.e_rot(l.kappa);
            let out = collide_full(&mut a, &mut b, n, l);
            if !free_speed && slip_n > 0.0 {
                pin_pair(&mut a, &mut b, l.kappa, e_rot_in, n.cross(out.j));
            }
            (a, b)
        };
        let err_vs = |(a, b): (Ball, Ball)| {
            let (ga, gb) = (to_ball(&qa), to_ball(&qb));
            let e = |x: DVec3, y: DVec3| (x - y).length();
            e(a.v, ga.v).max(e(b.v, gb.v)).max(e(a.w * a.r, ga.w * ga.r)).max(e(b.w * b.r, gb.w * gb.r))
        };
        res.checked += 1;
        res.max_err = res.max_err.max(err_vs(expect(&law)));
        res.max_err_wrong_law = res.max_err_wrong_law.max(err_vs(expect(&wrong)));
    }
    res
}

#[test]
fn rs_gpu_matches_cpu_reference() {
    let edge = RoughParams::default();
    let twist = RoughParams { beta: -1.0, beta_twist: 0.0, twist_loss_to_normal: true, ..RoughParams::default() };
    let both = RoughParams { beta: 1.0, beta_twist: 0.0, twist_loss_to_normal: true, ..RoughParams::default() };
    let cases = [
        ("edge, PinnedC", edge),
        ("edge, Free", RoughParams { speed_mode: SpeedMode::Free, g_max: 4.0, ..edge }),
        ("twist->normal, PinnedC", twist),
        ("twist->normal, Free", RoughParams { speed_mode: SpeedMode::Free, g_max: 4.0, ..twist }),
        ("edge+twist, Free", RoughParams { speed_mode: SpeedMode::Free, g_max: 4.0, ..both }),
    ];
    for (i, (name, rough)) in cases.into_iter().enumerate() {
        let r = run_pair_match(rough, 0x5EED_0000 + i as u64);
        println!(
            "RS gpu-vs-cpu [{name}]: accepted {} checked {} max err {:.2e} (wrong law: {:.2e}), unaccepted-but-changed {}",
            r.accepted, r.checked, r.max_err, r.max_err_wrong_law, r.unaccepted_changed
        );
        assert!(r.checked > 50, "{name}: too few collisions checked ({})", r.checked);
        assert!(r.max_err < 2e-3, "{name}: GPU disagrees with the CPU reference by {}", r.max_err);
        assert_eq!(r.unaccepted_changed, 0, "{name}: a rejected candidate changed state");
        // The wrong-law control only differs when the edge grip is engaged.
        if rough.beta != -1.0 {
            assert!(r.max_err_wrong_law > 0.05, "{name}: negative control did not fail ({})", r.max_err_wrong_law);
        }
    }
}

struct BoxStats {
    e_lin: f64,
    e_rot: f64,
    p: DVec3,
    p_scale: f64,
    v2: f64,
    v4: f64,
    max_speed_err: f64,
    n: usize,
}

fn box_stats(ph: &[GpuPhoton], kappa: f64) -> BoxStats {
    let mut s = BoxStats { e_lin: 0.0, e_rot: 0.0, p: DVec3::ZERO, p_scale: 0.0, v2: 0.0, v4: 0.0, max_speed_err: 0.0, n: 0 };
    for q in ph.iter().filter(|q| q.is_alive()) {
        let (el, er) = q.energies(kappa);
        let m = q.stiffness_mass();
        let v = q.velocity();
        let v = DVec3::new(v[0], v[1], v[2]);
        s.e_lin += el;
        s.e_rot += er;
        s.p += m * v;
        s.p_scale += m * v.length();
        s.v2 += v.length_squared();
        s.v4 += v.length_squared().powi(2);
        s.max_speed_err = s.max_speed_err.max((q.speed() - 1.0).abs());
        s.n += 1;
    }
    s
}

fn run_box(rough: RoughParams, n: u32, levels: u32, steps: u32, seed: u64) -> (BoxStats, BoxStats, Totals) {
    const HALF: f64 = 5.0;
    let sigma = 0.01;
    let config = box_config(HALF, 0.5, 0.05, sigma, n + 16, seed, rough);
    let ctx = GpuContext::new().expect("RS gates need a GPU adapter");
    let mut engine = DsmcEngine::new(&ctx, config);
    let mut rng = Lcg(seed);
    let photons: Vec<GpuPhoton> = (0..n)
        .map(|id| {
            let pos = DVec3::new(rng.f(), rng.f(), rng.f()) * (2.0 * HALF) - DVec3::splat(HALF);
            let k = (rng.f() * levels as f64) as u32;
            let d = rng.unit();
            let chi = if rng.f() < 0.5 { 1.0 } else { -1.0 };
            let w = d * (chi * rough.spin_edge0 / (2.0f64).powi(k as i32));
            make_photon(id, pos, d, w, k, chi)
        })
        .collect();
    engine.write_photons_for_test(&ctx, &photons);
    let s0 = box_stats(&engine.read_photons(&ctx), rough.kappa);
    let mut totals = Totals::default();
    for _ in 0..steps {
        engine.step_collide(&ctx, &mut totals);
    }
    let s1 = box_stats(&engine.read_photons(&ctx), rough.kappa);
    (s0, s1, totals)
}

#[test]
fn rs1_gpu_box_conservation() {
    // Mixed levels: light photons equilibrate faster (v^2 ~ 1/m), so the
    // relative-speed bound needs headroom well above 2.
    let free = RoughParams { speed_mode: SpeedMode::Free, g_max: 8.0, ..RoughParams::default() };
    let cases = [
        ("edge beta=1, Free", free, true),
        ("twist lock->normal + edge, Free", RoughParams { beta_twist: 0.0, twist_loss_to_normal: true, ..free }, true),
        ("NEG CONTROL edge beta=0, Free", RoughParams { beta: 0.0, ..free }, false),
    ];
    for (name, rough, conserving) in cases {
        let (s0, s1, t) = run_box(rough, 200_000, 3, 150, 0x0B0C_5001);
        let e0 = s0.e_lin + s0.e_rot;
        let de = (s1.e_lin + s1.e_rot - e0) / e0;
        let dp = (s1.p - s0.p).length() / s0.p_scale;
        let per_photon = 2.0 * t.collisions as f64 / s0.n as f64;
        println!(
            "RS1 gpu [{name}]: {:.1} collisions/photon, dE/E={de:+.2e}, |dP|/sum m|v|={dp:.2e}, accept overflow {}",
            per_photon, t.accept_overflow
        );
        assert!(per_photon > 5.0, "{name}: too few collisions ({per_photon})");
        assert_eq!(t.accept_overflow, 0, "{name}: g_max too small");
        assert!(dp < 1e-4, "{name}: momentum drift {dp}");
        if conserving {
            assert!(de.abs() < 1e-4, "{name}: energy drift {de}");
        } else {
            assert!(de < -0.05, "{name}: negative control did not lose energy ({de})");
        }
    }

    // PinnedC: every speed stays c, energy is conserved exactly (the pair's
    // spin energy is restored), and momentum is NOT (the price of pinning).
    let pinned = RoughParams::default();
    let (s0, s1, t) = run_box(pinned, 200_000, 3, 150, 0x0B0C_5002);
    let e0 = s0.e_lin + s0.e_rot;
    let de = (s1.e_lin + s1.e_rot - e0) / e0;
    let dp = (s1.p - s0.p).length() / s0.p_scale;
    println!(
        "RS1 gpu [edge, PinnedC]: dE/E={de:+.2e}, fallback failures {}, max||v|-1|={:.1e}, |dP|/sum m|v|={dp:.2e}, edge>c updates {}",
        t.pin_deficit_count, s1.max_speed_err, t.edge_over_c
    );
    assert!(s1.max_speed_err < 1e-5);
    assert_eq!(t.pin_deficit_count, 0);
    assert!(de.abs() < 1e-4, "PinnedC energy drift {de}");
}

#[test]
fn rs2_gpu_equipartition() {
    let free = RoughParams { speed_mode: SpeedMode::Free, g_max: 6.0, ..RoughParams::default() };
    let (s0, s1, t) = run_box(free, 200_000, 1, 300, 0x0B0C_5003);
    let ratio0 = s0.e_rot / s0.e_lin;
    let ratio = s1.e_rot / s1.e_lin;
    let kurt = (s1.v4 / s1.n as f64) / (s1.v2 / s1.n as f64).powi(2);
    println!(
        "RS2 gpu rough: {:.1} collisions/photon, E_rot/E_lin {ratio0:.3} -> {ratio:.4} (expect 1), <v^4>/<v^2>^2 = {kurt:.4} (Maxwell 5/3 = 1.6667), overflow {}",
        2.0 * t.collisions as f64 / s0.n as f64,
        t.accept_overflow
    );
    assert_eq!(t.accept_overflow, 0);
    assert!((ratio - 1.0).abs() < 0.02, "rough gas did not equipartition: {ratio}");
    assert!((kurt - 5.0 / 3.0).abs() < 0.03, "speed distribution not Maxwellian: {kurt}");

    let smooth = RoughParams { beta: -1.0, ..free };
    let (s0, s1, _) = run_box(smooth, 200_000, 1, 300, 0x0B0C_5004);
    let r_smooth = s1.e_rot / s1.e_lin;
    println!("RS2 gpu smooth (control): E_rot/E_lin {:.3} -> {r_smooth:.4} (must not move)", s0.e_rot / s0.e_lin);
    assert!((r_smooth - s0.e_rot / s0.e_lin).abs() < 1e-3, "smooth spheres moved spin energy");
}
