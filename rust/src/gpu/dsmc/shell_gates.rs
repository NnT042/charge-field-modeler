//! Doc §14 gates for permeable shell bodies (`BodyModel::Shells`).
//!
//! - `shell_gpu_matches_cpu_contact`: photons placed inside a body, one step
//!   at huge tau_body. Every GPU body contact must equal
//!   `rough_sphere::wall_contact` (+ restack) fed the same grain velocity
//!   (`baryon_shells::matter_velocity`) and the same contact normal. The
//!   kernel RNG is replicated bit for bit. Negative control: a CPU reference
//!   with the wrong grain velocity (outer levels dropped) must disagree.
//! - `shell_bg0_transparent`: tau_body = 0 means every passage is clean, with
//!   no contacts and no momentum delivered. tau_body > 0 must register
//!   contacts (control).
//! - `shell_bg1_static_body_conserves_photon_energy`: a non-rotating body
//!   (all level signs 0) is an elastic scatterer, so photon energy (linear +
//!   spin, f64 readback) is conserved, and the momentum the photons lose
//!   equals the tallied momentum delivered to the body. A rotating body must
//!   do work (control).
//! - `shell_bg2_bg3_fan`: a uniform fan (only level 12 spinning) takes
//!   angular momentum *out of* the body in its rotation sense. The torque on
//!   the body opposes the spin. Mirroring the fan (sign flip plus ambient
//!   chirality swap) flips the torque and helicity and leaves the exit bands
//!   unchanged (BG3).

use glam::{DQuat, DVec3};

use super::config::{
    BodyConfig, BodyModel, Boundary, DsmcConfig, EmissionMode, ExchangeRule, GpuClass, GpuOrient, MeshMode, RateKernel, RoughParams,
    SpeedMode, StiffnessShed,
};
use super::engine::{DsmcEngine, Totals};
use super::gpu_types::{GpuPhoton, FLAG_ALIVE, FLAG_SPECIES_PLUS};
use crate::baryon_shells::{matter_velocity, region_of, PROTON_STATES};
use crate::gpu::context::GpuContext;
use crate::rough_sphere::{contact_normal, level_radius, wall_contact, Ball, ContactLaw, StackRule, Stacked};

fn hash_u32(x: u32) -> u32 {
    let mut v = x;
    v ^= v >> 16;
    v = v.wrapping_mul(0x7feb_352d);
    v ^= v >> 15;
    v = v.wrapping_mul(0x846c_a68b);
    v ^= v >> 16;
    v
}

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

fn twist_rough(stacking: bool, mode: SpeedMode) -> RoughParams {
    RoughParams {
        beta: -1.0,
        beta_twist: 0.0,
        twist_loss_to_normal: true,
        speed_mode: mode,
        g_max: if mode == SpeedMode::Free { 6.0 } else { 2.0 },
        stacking,
        ..RoughParams::default()
    }
}

#[allow(clippy::too_many_arguments)]
fn shell_config(
    boundary: Boundary,
    signs: [f64; 4],
    tau_body: f64,
    rough: RoughParams,
    sigma: f64,
    n_a: f64,
    photon_fraction: f64,
    capacity: u32,
    seed: u64,
) -> DsmcConfig {
    DsmcConfig {
        boundary,
        bodies: vec![BodyConfig { center: DVec3::ZERO, class: GpuClass::Proton, orient: GpuOrient::PoleOn, through_lane_rho: 0.0 }],
        h: 0.25,
        dt: 0.02,
        capacity,
        seed,
        n_ambient: n_a,
        sigma,
        photon_fraction,
        pole_intake_bias: 0.0,
        cancel_redirect: 0.0,
        spin_transfer: 0.0,
        swing_boost: 0.0,
        rate_kernel: RateKernel::Catch,
        mesh_mode: MeshMode::DirRel,
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
        body_model: BodyModel::Shells { tau_body },
        shell_signs: vec![signs],
    }
}

fn make_photon(id: u32, pos: DVec3, v: DVec3, w: DVec3, chi: f64) -> GpuPhoton {
    let s = v.length();
    let d = v / s;
    let mut flags = FLAG_ALIVE;
    if chi > 0.0 {
        flags |= FLAG_SPECIES_PLUS;
    }
    GpuPhoton {
        a: [pos.x as f32, pos.y as f32, pos.z as f32, f32::from_bits(id)],
        b: [d.x as f32, d.y as f32, d.z as f32, w.dot(d) as f32],
        c: [flags, 0, 0, 0],
        d: [w.x as f32, w.y as f32, w.z as f32, s as f32],
        e: [0.0; 4],
        f: [0.0; 4],
        g: [0.0; 4],
    }
}

fn to_stacked(p: &GpuPhoton) -> Stacked {
    let v = p.velocity();
    let mut st = Stacked::new(DVec3::new(v[0], v[1], v[2]), DVec3::ZERO);
    st.k = p.c[3];
    for lvl in 0..=st.k.min(3) {
        let w = p.level_omega(lvl);
        st.w[lvl as usize] = DVec3::new(w[0], w[1], w[2]);
    }
    st
}

fn run_contact_match(rough: RoughParams, seed: u64) -> (usize, f64, f64, f64) {
    let n: u32 = 4000;
    let tau = 1.0e5;
    let config = shell_config(Boundary::PeriodicBox { half_extent: 2.0 }, PROTON_STATES[0], tau, rough, 0.0, 0.0, 0.5, n + 16, seed);
    let ctx = GpuContext::new().expect("shell gates need a GPU adapter");
    let mut engine = DsmcEngine::new(&ctx, config);
    let body = engine.shell_bodies[0].clone();
    let frames = body.frames();
    let mut rng = Lcg(seed ^ 0x51E1);
    let free = rough.speed_mode == SpeedMode::Free;
    let photons: Vec<GpuPhoton> = (0..n)
        .map(|id| {
            let pos = loop {
                let x = DVec3::new(2.0 * rng.f() - 1.0, 2.0 * rng.f() - 1.0, 2.0 * rng.f() - 1.0) * 0.98;
                if x.length() < 0.98 {
                    break x;
                }
            };
            let speed = if free { 0.5 + rng.f() } else { 1.0 };
            let v = rng.unit() * speed;
            let w = rng.unit() * rng.f();
            make_photon(id, pos, v, w, if rng.f() < 0.5 { 1.0 } else { -1.0 })
        })
        .collect();
    engine.write_photons_for_test(&ctx, &photons);
    let step = engine.step;
    let mut totals = Totals::default();
    engine.step_collide(&ctx, &mut totals);
    let after = engine.read_photons(&ctx);
    let law = ContactLaw { e: 1.0, beta: rough.beta, beta_twist: rough.beta_twist, twist_loss_to_normal: rough.twist_loss_to_normal, kappa: rough.kappa };
    let rule = StackRule { max_level: rough.stack_max_level, stack_tol: 1e-4, unstack_edge: rough.unstack_edge, kappa: rough.kappa };
    let (mut checked, mut max_err, mut max_err_wrong) = (0usize, 0.0f64, 0.0f64);
    let mut dp_cpu = DVec3::ZERO;
    for p in &photons {
        let id = p.id();
        let q = after.iter().find(|q| q.id() == id).expect("photon lost");
        let x = DVec3::new(p.a[0] as f64, p.a[1] as f64, p.a[2] as f64);
        let Some(k) = region_of(&frames, x) else { continue };
        let gv = matter_velocity(&frames, k, x);
        let st0 = to_stacked(p);
        let u = st0.v - gv;
        let key = hash_u32(id ^ hash_u32(step ^ 0xB0D1E5));
        if u.length() < 1e-3 || rand01(0, 0, 0).is_nan() {
            continue;
        }
        let _ = rand01(engine.config.seed, key, 0); // p_hit ~ 1 at this tau
        let mut rho = rand01(engine.config.seed, key, 1);
        if rough.glancing_bias {
            let om = 1.0 - rho;
            rho = 1.0 - om * om * om;
        }
        let psi = rand01(engine.config.seed, key, 2) as f64 * std::f64::consts::TAU;
        let nrm = contact_normal(u, rho as f64, psi);
        let b0 = st0.outer_ball();
        let slip_n = ((b0.v - gv) + b0.r * b0.w.cross(nrm)).dot(nrm);
        if slip_n.abs() < 1e-4 {
            continue;
        }
        let expect = |gv_used: DVec3| {
            let mut st = st0;
            let mut b: Ball = st.outer_ball();
            let out = wall_contact(&mut b, gv_used, nrm, &law, !free);
            st.set_outer(&b);
            if out.touched && rough.stacking {
                st.restack(out.kick, &rule);
            }
            (st, out.dp_body)
        };
        let gpu = to_stacked(q);
        let err = |st: &Stacked| {
            if st.k != gpu.k {
                return f64::INFINITY;
            }
            let mut e = (st.v - gpu.v).length();
            for lvl in 0..=st.k as usize {
                e = e.max((st.w[lvl] - gpu.w[lvl]).length() * level_radius(lvl as u32));
            }
            e
        };
        let (good, dp) = expect(gv);
        dp_cpu += dp;
        let wrong_gv = frames[k].omega.cross(x - frames[k].origin);
        let (bad, _) = expect(wrong_gv);
        checked += 1;
        max_err = max_err.max(err(&good));
        if k < 3 {
            max_err_wrong = max_err_wrong.max(err(&bad));
        }
    }
    // Tally consistency over ALL photons (including the few borderline
    // ones the per-photon comparison skips): the delivered-momentum tally
    // must equal minus the photons' own momentum change (unit mass here).
    let _ = dp_cpu;
    let mut dp_state = DVec3::ZERO;
    for p in &photons {
        let q = after.iter().find(|q| q.id() == p.id()).expect("photon lost");
        let (v0, v1) = (p.velocity(), q.velocity());
        dp_state -= DVec3::new(v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]);
    }
    let dp_gpu = DVec3::new(totals.sh_dp[0][0], totals.sh_dp[0][1], totals.sh_dp[0][2]);
    let dp_rel = (dp_gpu - dp_state).length() / dp_state.length().max(1e-9);
    (checked, max_err, max_err_wrong, dp_rel)
}

#[test]
fn shell_gpu_matches_cpu_contact() {
    for (name, rough) in [
        ("twist, PinnedC", twist_rough(false, SpeedMode::PinnedC)),
        ("twist, PinnedC, stacking", twist_rough(true, SpeedMode::PinnedC)),
        ("edge, Free", RoughParams { speed_mode: SpeedMode::Free, g_max: 6.0, ..RoughParams::default() }),
        ("edge+twist, PinnedC, stacking, glancing", RoughParams { beta: 1.0, glancing_bias: true, ..twist_rough(true, SpeedMode::PinnedC) }),
    ] {
        let (checked, err, wrong, dp_rel) = run_contact_match(rough, 0x5E11_0001);
        println!("shell contact gpu-vs-cpu [{name}]: checked {checked}, max err {err:.2e} (wrong grain velocity: {wrong:.2e}), delivered-momentum rel diff {dp_rel:.1e}");
        assert!(checked > 1000);
        assert!(err < 2e-3, "{name}: GPU body contact disagrees with CPU: {err}");
        assert!(wrong > 0.05, "{name}: negative control did not fail");
        assert!(dp_rel < 1e-2, "{name}: delivered momentum tally disagrees: {dp_rel}");
    }
}

fn run_flux(config: DsmcConfig, steps: u32) -> (Totals, Vec<GpuPhoton>) {
    let ctx = GpuContext::new().expect("shell gates need a GPU adapter");
    let mut engine = DsmcEngine::new(&ctx, config);
    let mut t = Totals::default();
    for _ in 0..steps {
        t.add(&engine.step_flux(&ctx));
    }
    let ph = engine.read_photons(&ctx);
    (t, ph)
}

#[test]
fn shell_bg0_transparent() {
    let rough = twist_rough(true, SpeedMode::PinnedC);
    let mk = |tau| shell_config(Boundary::SphereEscape { r_d: 3.0 }, PROTON_STATES[0], tau, rough, 0.0, 2000.0, 2.0 / 3.0, 3_000_000, 0xB60);
    let (t0, _) = run_flux(mk(0.0), 400);
    let (t1, _) = run_flux(mk(5.0), 400);
    println!(
        "BG0: tau=0 entries {} clean {} contacts {} exits {} |dp| {:.3}; tau=5 contacts {} exits {}",
        t0.sh_entries[0],
        t0.sh_clean[0],
        t0.sh_contacts[0],
        t0.sh_exits[0],
        (t0.sh_dp[0][0].powi(2) + t0.sh_dp[0][1].powi(2) + t0.sh_dp[0][2].powi(2)).sqrt(),
        t1.sh_contacts[0],
        t1.sh_exits[0]
    );
    assert!(t0.sh_entries[0] > 10_000 && t0.sh_clean[0] > 10_000);
    assert_eq!(t0.sh_contacts[0], 0);
    assert_eq!(t0.sh_exits[0], 0);
    assert!(t1.sh_contacts[0] > 10_000 && t1.sh_exits[0] > 5_000, "control: a permeable body with tau > 0 must scatter");
}

fn photon_energy(ph: &[GpuPhoton], kappa: f64) -> (f64, DVec3) {
    let mut e = 0.0;
    let mut p = DVec3::ZERO;
    for q in ph.iter().filter(|q| q.is_alive()) {
        let (l, r) = q.energies_stacked(kappa);
        e += l + r;
        let v = q.velocity();
        p += DVec3::new(v[0], v[1], v[2]);
    }
    (e, p)
}

#[test]
fn shell_bg1_static_body_conserves_photon_energy() {
    let rough = RoughParams { beta: 1.0, speed_mode: SpeedMode::Free, g_max: 6.0, stacking: true, ..RoughParams::default() };
    for (name, signs, conserving) in [("static", [0.0; 4], true), ("CONTROL rotating", PROTON_STATES[0], false)] {
        let n: u32 = 200_000;
        let config = shell_config(Boundary::PeriodicBox { half_extent: 2.0 }, signs, 3.0, rough, 0.0, 0.0, 0.5, n + 16, 0xB61);
        let ctx = GpuContext::new().expect("shell gates need a GPU adapter");
        let mut engine = DsmcEngine::new(&ctx, config);
        let mut rng = Lcg(0xB61);
        let photons: Vec<GpuPhoton> = (0..n)
            .map(|id| {
                let pos = DVec3::new(rng.f(), rng.f(), rng.f()) * 4.0 - DVec3::splat(2.0);
                let d = rng.unit();
                let chi = if rng.f() < 0.5 { 1.0 } else { -1.0 };
                make_photon(id, pos, d, d * chi, chi)
            })
            .collect();
        engine.write_photons_for_test(&ctx, &photons);
        let (e0, p0) = photon_energy(&engine.read_photons(&ctx), rough.kappa);
        let mut t = Totals::default();
        for _ in 0..200 {
            engine.step_collide(&ctx, &mut t);
        }
        let (e1, p1) = photon_energy(&engine.read_photons(&ctx), rough.kappa);
        let de = (e1 - e0) / e0;
        let dp_body = DVec3::new(t.sh_dp[0][0], t.sh_dp[0][1], t.sh_dp[0][2]);
        let balance = ((p1 - p0) + dp_body).length() / dp_body.length().max(1.0);
        println!(
            "BG1 [{name}]: contacts {} dE/E {de:+.2e}, |dP_photons + dP_body|/|dP_body| {balance:.2e} (|dP_body| {:.1})",
            t.sh_contacts[0],
            dp_body.length()
        );
        assert!(t.sh_contacts[0] > 100_000);
        assert!(balance < 1e-2, "{name}: momentum given to the body does not match what the photons lost");
        if conserving {
            assert!(de.abs() < 1e-4, "static body changed photon energy: {de}");
        } else {
            assert!(de.abs() > 1e-2, "control: a spinning body must do work on the photons ({de})");
        }
    }
}

#[test]
fn shell_bg2_bg3_fan() {
    let rough = twist_rough(true, SpeedMode::PinnedC);
    let rough_edge = RoughParams { beta: 1.0, ..rough };
    let mut rows = Vec::new();
    for (name, r) in [("twist", rough), ("edge+twist", rough_edge)] {
        for (sign, pf) in [(1.0, 2.0 / 3.0), (-1.0, 1.0 / 3.0)] {
            let config =
                shell_config(Boundary::SphereEscape { r_d: 3.0 }, [0.0, 0.0, 0.0, sign], 4.0, r, 0.0, 2000.0, pf, 3_000_000, 0xB62);
            let (t, _) = run_flux(config, 800);
            let exits: i64 = t.sh_exits[0].max(1);
            let bands: Vec<f64> = (0..3).map(|b| t.sh_exit_band[0][b] as f64 / exits as f64).collect();
            let hel: f64 = (0..3).map(|b| t.sh_exit_band_hel[0][b]).sum::<f64>() / exits as f64;
            let tz = t.sh_torque[0][2] / exits as f64;
            println!(
                "BG2/3 fan [{name}] sign {sign:+} pf {pf:.2}: exits {} bands p/m/e {:.3}/{:.3}/{:.3} hel {hel:+.4} torque_z/exit {tz:+.5}",
                t.sh_exits[0], bands[0], bands[1], bands[2]
            );
            rows.push((name, sign, bands, hel, tz, exits));
        }
    }
    for pair in rows.chunks(2) {
        let (a, b) = (&pair[0], &pair[1]);
        // BG2: the torque on the body opposes its spin (photons carry
        // angular momentum away in the rotation sense).
        assert!(a.4 < 0.0 && b.4 > 0.0, "{}: torque does not oppose the fan's spin ({} / {})", a.0, a.4, b.4);
        // BG3: mirror flips torque and helicity, keeps bands.
        let se = |x: f64, n: i64| (x * (1.0 - x) / n as f64).sqrt();
        for band in 0..3 {
            let z = (a.2[band] - b.2[band]).abs() / (se(a.2[band], a.5).powi(2) + se(b.2[band], b.5).powi(2)).sqrt();
            assert!(z < 4.0, "{}: mirrored fan changed exit band {band} by {z:.1} sigma", a.0);
        }
        assert!((a.4 + b.4).abs() < 0.1 * a.4.abs(), "{}: torque not antisymmetric", a.0);
        assert!(a.3 * b.3 <= 0.0 || (a.3 + b.3).abs() < 0.02, "{}: helicity not antisymmetric", a.0);
    }
}

#[allow(dead_code)]
fn _unused(_: DQuat) {}

/// Doc §14 first physics run: a single mechanical baryon (all four proton
/// and four neutron sign states) in the 2:1 ambient field, sweeping tau_body.
/// Twist grip, PinnedC, stacking on, photon-photon lambda = 3. The questions:
/// does exit topology split proton vs neutron, does anything pass
/// pole-to-pole, and which bands do contacted photons leave through? An
/// isotropic flux crossing a sphere exits 0.134 / 0.366 / 0.500
/// (polar / mid / equatorial by area).
/// Env: CFM_B_TAUS (default "1,4,16"), CFM_B_BIAS (default "0,0.5"),
/// CFM_B_WARMUP / CFM_B_AVG (400 / 1600).
#[test]
#[ignore]
fn report_b_shell_single_body() {
    use crate::baryon_shells::NEUTRON_STATES;
    use std::time::Instant;
    let env_list = |k: &str, d: &str| -> Vec<f64> {
        std::env::var(k).unwrap_or(d.to_string()).split(',').filter_map(|x| x.trim().parse().ok()).collect()
    };
    let taus = env_list("CFM_B_TAUS", "1,4,16");
    let biases = env_list("CFM_B_BIAS", "0,0.5");
    let warmup: u32 = std::env::var("CFM_B_WARMUP").ok().and_then(|v| v.parse().ok()).unwrap_or(400);
    let avg: u32 = std::env::var("CFM_B_AVG").ok().and_then(|v| v.parse().ok()).unwrap_or(1600);
    const N_A: f64 = 400.0;
    const LAMBDA: f64 = 3.0;
    let rough = twist_rough(true, SpeedMode::PinnedC);
    println!("=== B: single mechanical baryon (shells), twist/PinnedC/stacking, lambda={LAMBDA}, n_a={N_A}, R_D=6, warmup {warmup}, avg {avg} ===");
    println!("isotropic crossing baseline bands: polar 0.134 mid 0.366 equat 0.500");
    println!(
        "{:<9} {:<10} {:>5} {:>4} | {:>7} {:>6} {:>6} | {:>6} {:>6} {:>6} | {:>7} {:>7} {:>7} | {:>7} {:>8} {:>8} {:>6}",
        "class", "state", "tau", "bias", "entry/s", "cont%", "c/pass", "polar", "mid", "equat", "hel_p", "hel_m", "hel_e", "thru%", "Tz/exit", "W/exit", "steady"
    );
    let ctx = GpuContext::new().expect("GPU");
    for &bias in &biases {
        for &tau in &taus {
            for (class, states) in [("proton", PROTON_STATES), ("neutron", NEUTRON_STATES)] {
                for signs in states {
                    let started = Instant::now();
                    let mut config =
                        shell_config(Boundary::SphereEscape { r_d: 6.0 }, signs, tau, rough, 1.0 / (LAMBDA * N_A), N_A, 2.0 / 3.0, 6_000_000, 0xB0B0);
                    config.dt = 0.05;
                    config.pole_intake_bias = bias;
                    let mut engine = DsmcEngine::new(&ctx, config);
                    let rep = engine.run_steady(&ctx, warmup, avg, 200);
                    let t = &rep.totals;
                    let exits = t.sh_exits[0].max(1) as f64;
                    let band = |b: usize| t.sh_exit_band[0][b] as f64 / exits;
                    let hel = |b: usize| t.sh_exit_band_hel[0][b] / (t.sh_exit_band[0][b].max(1) as f64);
                    let st: String = signs.iter().zip(['a', 'x', 'y', 'z']).map(|(s, c)| format!("{}{}", if *s > 0.0 { '+' } else { '-' }, c)).collect();
                    println!(
                        "{:<9} {:<10} {:>5} {:>4} | {:>7.0} {:>6.1} {:>6.2} | {:>6.3} {:>6.3} {:>6.3} | {:>+7.3} {:>+7.3} {:>+7.3} | {:>7.3} {:>+8.4} {:>+8.4} {:>6}  ({:.0}s)",
                        class,
                        st,
                        tau,
                        bias,
                        t.sh_entries[0] as f64 / (avg as f64 * 0.05),
                        100.0 * t.sh_exits[0] as f64 / t.sh_entries[0].max(1) as f64,
                        t.sh_contacts[0] as f64 / exits,
                        band(0),
                        band(1),
                        band(2),
                        hel(0),
                        hel(1),
                        hel(2),
                        100.0 * t.sh_through[0] as f64 / exits,
                        t.sh_torque[0][2] / exits,
                        t.sh_work[0] / exits,
                        rep.steady_ok,
                        started.elapsed().as_secs_f64()
                    );
                }
            }
        }
    }
}
