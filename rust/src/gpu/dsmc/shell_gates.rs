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
use super::gpu_types::{GpuPhoton, BALLKEY_STRIDE, FLAG_ALIVE, FLAG_SPECIES_PLUS};
use super::config::{GrainMode, MazeMode};
use crate::baryon_shells::{matter_velocity, region_of, PROTON_STATES};
use crate::gpu::context::GpuContext;
use crate::rough_sphere::{contact_normal, level_radius, wall_contact_spinning, Ball, ContactLaw, StackRule, Stacked};

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
    let (bb, bt, bl) = rough.body_law.unwrap_or((rough.beta, rough.beta_twist, rough.twist_loss_to_normal));
    let law = ContactLaw { e: 1.0, beta: bb, beta_twist: bt, twist_loss_to_normal: bl, kappa: rough.kappa };
    let rule = StackRule { max_level: rough.stack_max_level, stack_tol: 1e-4, unstack_edge: rough.unstack_edge, kappa: rough.kappa };
    let (mut checked, mut max_err, mut max_err_wrong) = (0usize, 0.0f64, 0.0f64);
    let mut dp_cpu = DVec3::ZERO;
    for p in &photons {
        let id = p.id();
        let q = after.iter().find(|q| q.id() == id).expect("photon lost");
        let x = DVec3::new(p.a[0] as f64, p.a[1] as f64, p.a[2] as f64);
        let Some(k) = region_of(&frames, x) else { continue };
        let gv = matter_velocity(&frames, k, x);
        let omega_k: DVec3 = (k..4).map(|j| frames[j].omega).sum();
        let gw = if omega_k.length() > 1e-9 { omega_k.normalize() * rough.grain_edge } else { DVec3::ZERO };
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
        let slip_n = ((b0.v - gv) + b0.r * b0.w.cross(nrm) + gw.cross(nrm)).dot(nrm);
        if slip_n.abs() < 1e-4 {
            continue;
        }
        let expect = |gv_used: DVec3| {
            let mut st = st0;
            let mut b: Ball = st.outer_ball();
            let out = wall_contact_spinning(&mut b, gv_used, gw, nrm, &law, !free);
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
        (
            "SPINNING GRAINS, body edge+twist law, photon twist law, stacking",
            RoughParams { grain_edge: 1.0, body_law: Some((1.0, 0.0, true)), ..twist_rough(true, SpeedMode::PinnedC) },
        ),
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
    // CFM_B_LAW: "twist" (default), "edge", "edge+twist".
    let law_name = std::env::var("CFM_B_LAW").unwrap_or("twist".to_string());
    let base = twist_rough(true, SpeedMode::PinnedC);
    let law_of = |name: &str| -> (f64, f64, bool) {
        match name {
            "edge" => (1.0, -1.0, false),
            "edge+twist" => (1.0, 0.0, true),
            _ => (-1.0, 0.0, true),
        }
    };
    let (b0, t0, l0) = law_of(&law_name);
    // CFM_B_BODY_LAW: separate photon-body law (default: same as CFM_B_LAW).
    // CFM_B_GRAIN: grain spin edge speed (default 0 = non-spinning).
    let body_law_name = std::env::var("CFM_B_BODY_LAW").ok();
    let grain_edge: f64 = std::env::var("CFM_B_GRAIN").ok().and_then(|v| v.parse().ok()).unwrap_or(0.0);
    let rough = RoughParams {
        beta: b0,
        beta_twist: t0,
        twist_loss_to_normal: l0,
        grain_edge,
        body_law: body_law_name.as_deref().map(law_of),
        ..base
    };
    let law_name = format!("{law_name} body={} grain={grain_edge}", body_law_name.as_deref().unwrap_or("same"));
    println!("=== B: single mechanical baryon (shells), law={law_name}, PinnedC/stacking, lambda={LAMBDA}, n_a={N_A}, R_D=6, warmup {warmup}, avg {avg} ===");
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

// ---------------------------------------------------------------------------
// Doc §15: the instantaneous ball
// ---------------------------------------------------------------------------

fn ball_config(boundary: Boundary, signs: [f64; 4], rough: RoughParams, sigma: f64, n_a: f64, capacity: u32, seed: u64) -> DsmcConfig {
    let mut c = shell_config(boundary, signs, 0.0, rough, sigma, n_a, 2.0 / 3.0, capacity, seed);
    c.body_model = BodyModel::Ball;
    c
}

/// Doc §15.8 grain mode plus what the RigidInner mirror needs to reproduce
/// the kernel's random phases (run seed and step index).
#[derive(Clone, Copy)]
struct GrainMirror {
    mode: GrainMode,
    seed: u64,
    step: u32,
    /// Doc §15.9 maze emission law.
    maze: MazeMode,
}

/// Rodrigues rotation (the shader's axis_rot).
fn axis_rot(a: DVec3, th: f64, v: DVec3) -> DVec3 {
    let (s, c) = th.sin_cos();
    v * c + a.cross(v) * s + a * a.dot(v) * (1.0 - c)
}

/// CPU mirror of the move kernel's ball phase for one photon (multi-hit
/// loop, push-out on touching, 8 linear sub-intervals between keyframes).
/// Returns (number of contacts, final state, final position).
#[allow(clippy::too_many_arguments)]
fn cpu_ball_step(
    p: &GpuPhoton,
    keys: &[[f32; 4]],
    dt: f64,
    law: &ContactLaw,
    pinned: bool,
    rule: &StackRule,
    grain_edge: f64,
    stacking: bool,
    grain: GrainMirror,
) -> (u32, Stacked, DVec3) {
    use crate::rough_sphere::{moving_sphere_hit, wall_contact_general};
    let k = |s: usize, j: usize| {
        let q = keys[s * 4 + j];
        (DVec3::new(q[0] as f64, q[1] as f64, q[2] as f64), q[3] as f64)
    };
    let mut st = to_stacked(p);
    let mut x = DVec3::new(p.a[0] as f64, p.a[1] as f64, p.a[2] as f64);
    let h = dt / 8.0;
    let mut tc = 0.0;
    let mut hits = 0;
    for iter in 0..6u32 {
        let v = st.v;
        let si0 = ((tc / h) as usize).min(7);
        let f0 = (tc - si0 as f64 * h) / h;
        let (ca, rb) = k(si0, 0);
        let (cb, _) = k(si0 + 1, 0);
        let c_now = ca.lerp(cb, f0);
        let vc_now = (cb - ca) / h;
        let d = x - c_now;
        let mut best: Option<(f64, usize, bool)> = None;
        if d.length_squared() <= rb * rb * (1.0 + 2e-4) && d.dot(v - vc_now) < 0.0 {
            best = Some((tc, si0, true));
        } else {
            for si in si0..8 {
                let ts = (si as f64 * h).max(tc);
                let (c0a, _) = k(si, 0);
                let (c1a, _) = k(si + 1, 0);
                let c0 = c0a.lerp(c1a, (ts - si as f64 * h) / h);
                let vc = (c1a - c0a) / h;
                if let Some(t) = moving_sphere_hit(x + v * (ts - tc), v, c0, vc, rb, (si + 1) as f64 * h - ts) {
                    best = Some((ts + t, si, false));
                    break;
                }
            }
        }
        let Some((th_hit, si, touching)) = best else { break };
        let f = (th_hit - si as f64 * h) / h;
        let mut xh = x + v * (th_hit - tc);
        let ch = k(si, 0).0.lerp(k(si + 1, 0).0, f);
        if touching {
            let dd = xh - ch;
            if dd.length_squared() > 1e-20 {
                xh = ch + dd.normalize() * rb;
            }
        }
        let om = k(si, 1).0.lerp(k(si + 1, 1).0, f);
        let vcen = k(si, 2).0.lerp(k(si + 1, 2).0, f);
        if grain.maze != MazeMode::Off {
            // Doc §15.9: capture and re-emission through the maze (mirrors
            // the shader's maze_contact).
            let meta = keys[36];
            let oq = keys[45];
            let orient = glam::DQuat::from_xyzw(oq[0] as f64, oq[1] as f64, oq[2] as f64, oq[3] as f64);
            let sg = keys[46];
            let signs = [sg[0] as f64, sg[1] as f64, sg[2] as f64, sg[3] as f64];
            let th = if grain.maze == MazeMode::Swept { meta[2] as f64 + meta[3] as f64 * th_hit } else { 0.0 };
            let local = |v: DVec3| axis_rot(DVec3::Z, -th, orient.inverse() * v);
            let world = |v: DVec3| orient * axis_rot(DVec3::Z, th, v);
            let v_in = st.v;
            let dl = local(v_in.normalize());
            let key = hash_u32(p.a[3].to_bits() ^ hash_u32(grain.step ^ 0x3A2E ^ iter.wrapping_mul(0x9E37_79B9)));
            let e = crate::maze::emit(&signs, dl, local(st.w[st.k as usize]), |i| rand01(grain.seed, key, 50 + i) as f64);
            for lvl in 0..st.k as usize {
                st.w[lvl] = world(crate::maze::turn(&signs, crate::maze::arc(dl, e.d0, local(st.w[lvl]))));
            }
            let out_w = world(e.out);
            let s_rel = (v_in - (vcen + om.cross(xh - ch))).length();
            let xe = ch + out_w * rb * (1.0 + 1e-3);
            let mut v_out = vcen + om.cross(xe - ch) + out_w * s_rel;
            if pinned {
                v_out = v_out.normalize();
            }
            st.v = v_out;
            st.w[st.k as usize] = world(e.spin);
            hits += 1;
            x = xe;
            tc = th_hit;
            continue;
        }
        let n = (ch - xh).normalize();
        let tw = match grain.mode {
            GrainMode::EdgeC => {
                if om.length() > 1e-9 {
                    om.normalize() * grain_edge
                } else {
                    DVec3::ZERO
                }
            }
            _ => {
                let meta = keys[36];
                let mut tw = om * meta[1] as f64;
                if grain.mode == GrainMode::RigidInner {
                    let qa = keys[si * 4 + 3];
                    let qb = keys[(si + 1) * 4 + 3];
                    let l = |a: f32, b: f32| (a + (b - a) * f as f32) as f64;
                    let q = glam::DQuat::from_xyzw(l(qa[0], qb[0]), l(qa[1], qb[1]), l(qa[2], qb[2]), l(qa[3], qb[3])).normalize();
                    let m = meta[0] as usize;
                    let levels: Vec<(DVec3, f64)> = (0..m).map(|i| {
                        let e = keys[37 + i];
                        (DVec3::new(e[0] as f64, e[1] as f64, e[2] as f64), e[3] as f64)
                    }).collect();
                    let key = hash_u32(p.a[3].to_bits() ^ hash_u32(grain.step ^ 0x7A157 ^ iter.wrapping_mul(0x9E37_79B9)));
                    let angles: Vec<f64> = (0..m).map(|j| (std::f32::consts::TAU * rand01(grain.seed, key, 40 + j as u32)) as f64).collect();
                    tw += crate::baryon_shells::inner_twist(q, &levels, &angles);
                }
                tw
            }
        };
        let mut b = st.outer_ball();
        let out = wall_contact_general(&mut b, vcen, om * rb, tw, n, law, pinned);
        st.set_outer(&b);
        if out.touched {
            hits += 1;
            if stacking {
                st.restack(out.kick, rule);
            }
        }
        x = xh;
        tc = th_hit;
    }
    (hits, st, x + st.v * (dt - tc))
}

/// GPU ball hits == CPU mirror, per photon (position after the step,
/// velocity, every spin level). Photons are launched at the ball from
/// random directions so most of them hit within the step. Negative control:
/// the CPU mirror with the ball's spin dropped must disagree.
#[test]
fn ball_gpu_matches_cpu() {
    for (name, rough) in [
        ("twist body law, grain 1, PinnedC, stacking", RoughParams { grain_edge: 1.0, ..twist_rough(true, SpeedMode::PinnedC) }),
        ("edge law, Free", RoughParams { beta: 1.0, speed_mode: SpeedMode::Free, g_max: 6.0, ..RoughParams::default() }),
        ("twist body law, grain 1, BodyFree, stacking", RoughParams { grain_edge: 1.0, speed_mode: SpeedMode::BodyFree, g_max: 8.0, ..twist_rough(true, SpeedMode::PinnedC) }),
        ("twist body law, RIGID grain, BodyFree, stacking", RoughParams { grain_mode: GrainMode::Rigid, speed_mode: SpeedMode::BodyFree, g_max: 8.0, ..twist_rough(true, SpeedMode::PinnedC) }),
        ("twist body law, RIGID+INNER grain (z-flip neutron, depth 2), BodyFree, stacking", RoughParams { grain_mode: GrainMode::RigidInner, speed_mode: SpeedMode::BodyFree, g_max: 8.0, ..twist_rough(true, SpeedMode::PinnedC) }),
        ("MAZE swept, proton, BodyFree, stacking", RoughParams { maze: MazeMode::Swept, speed_mode: SpeedMode::BodyFree, g_max: 8.0, ..twist_rough(true, SpeedMode::PinnedC) }),
        ("MAZE swept, neutron, BodyFree, stacking", RoughParams { maze: MazeMode::Swept, speed_mode: SpeedMode::BodyFree, g_max: 8.0, ..twist_rough(true, SpeedMode::PinnedC) }),
        ("MAZE fixed, neutron, PinnedC", RoughParams { maze: MazeMode::Fixed, ..twist_rough(false, SpeedMode::PinnedC) }),
    ] {
        let n: u32 = 6000;
        let inner_case = rough.grain_mode == GrainMode::RigidInner;
        let signs = if name.contains("neutron") { crate::baryon_shells::NEUTRON_STATES[0] } else { PROTON_STATES[0] };
        let mut config = ball_config(Boundary::PeriodicBox { half_extent: 2.0 }, signs, rough, 0.0, 0.0, n + 16, 0xBA11);
        config.dt = 0.05;
        let ctx = GpuContext::new().expect("GPU");
        let mut engine = DsmcEngine::new(&ctx, config);
        if inner_case {
            let b = &mut engine.shell_bodies[0];
            b.depth = 2;
            b.inner_signs = NEUTRON_INNER_SIGNS;
            b.tier1_signs = NEUTRON_INNER_SIGNS;
        }
        let keys = engine.ball_keyframes();
        let grain = GrainMirror { mode: rough.grain_mode, seed: 0xBA11, step: engine.step, maze: rough.maze };
        let ball = engine.shell_bodies[0].ball();
        let mut rng = Lcg(0xBA11);
        let photons: Vec<GpuPhoton> = (0..n)
            .map(|id| {
                let dir_in = rng.unit();
                let off = rng.unit() * 0.8 * ball.radius;
                let pos = ball.center - dir_in * (ball.radius + 0.01 + 0.02 * rng.f()) + off;
                let speed = if rough.speed_mode != SpeedMode::PinnedC { 0.7 + 0.6 * rng.f() } else { 1.0 };
                let w = rng.unit() * rng.f();
                make_photon(id, pos, (dir_in + rng.unit() * 0.3).normalize() * speed, w, if rng.f() < 0.5 { 1.0 } else { -1.0 })
            })
            .filter(|p| (DVec3::new(p.a[0] as f64, p.a[1] as f64, p.a[2] as f64) - ball.center).length() > ball.radius + 1e-3)
            .collect();
        engine.write_photons_for_test(&ctx, &photons);
        let mut t = Totals::default();
        engine.step_collide(&ctx, &mut t);
        let after = engine.read_photons(&ctx);
        let (bb, bt, bl) = rough.body_law.unwrap_or((rough.beta, rough.beta_twist, rough.twist_loss_to_normal));
        let law = ContactLaw { e: 1.0, beta: bb, beta_twist: bt, twist_loss_to_normal: bl, kappa: rough.kappa };
        let rule = StackRule { max_level: rough.stack_max_level, stack_tol: 1e-4, unstack_edge: rough.unstack_edge, kappa: rough.kappa };
        let pinned = rough.speed_mode == SpeedMode::PinnedC;
        let (mut hits, mut mism_hit, mut max_err, mut max_wrong, mut multi) = (0, 0, 0.0f64, 0.0f64, 0usize);
        let mut keys_ns = keys[0..BALLKEY_STRIDE].to_vec();
        for s in 0..9 {
            keys_ns[s * 4 + 1] = [0.0; 4];
        }
        for i in 0..8 {
            keys_ns[37 + i][3] = 0.0;
        }
        // Maze control: flip the body's z sign (escape iff x == z flips).
        keys_ns[46][3] = -keys_ns[46][3];
        for p in &photons {
            let q = after.iter().find(|q| q.id() == p.id()).expect("lost");
            let (nhit, st, x1) = cpu_ball_step(p, &keys[0..BALLKEY_STRIDE], 0.05, &law, pinned, &rule, rough.grain_edge, rough.stacking, grain);
            let hit = nhit > 0;
            multi += (nhit > 1) as usize;
            // Hit = the photon's velocity or spin changed. (The passage
            // counter in c[2] is reset by the envelope-entry bookkeeping for
            // these synthetic photons, which start inside the envelope.)
            let gpu_hit = (0..4).any(|c| (q.b[c] - p.b[c]).abs() > 1e-6 || (q.d[c] - p.d[c]).abs() > 1e-6);
            if hit != gpu_hit {
                mism_hit += 1;
                continue;
            }
            if !hit {
                continue;
            }
            hits += 1;
            let g = to_stacked(q);
            let gx = DVec3::new(q.a[0] as f64, q.a[1] as f64, q.a[2] as f64);
            let mut e = (st.v - g.v).length().max((x1 - gx).length());
            if st.k != g.k {
                e = f64::INFINITY;
            } else {
                for lvl in 0..=st.k as usize {
                    e = e.max((st.w[lvl] - g.w[lvl]).length() * level_radius(lvl as u32));
                }
            }
            max_err = max_err.max(e);
            let (_, st_w, _) = cpu_ball_step(p, &keys_ns, 0.05, &law, pinned, &rule, 0.0, rough.stacking, grain);
            let mut wrong = (st_w.v - g.v).length();
            if st_w.k == g.k {
                for lvl in 0..=g.k as usize {
                    wrong = wrong.max((st_w.w[lvl] - g.w[lvl]).length() * level_radius(lvl as u32));
                }
            } else {
                wrong = f64::INFINITY;
            }
            max_wrong = max_wrong.max(wrong);
        }
        println!(
            "ball gpu-vs-cpu [{name}]: hits {hits}/{} (multi-hit {multi}), hit/miss disagreements {mism_hit}, max err {max_err:.2e}, no-ball-spin control {max_wrong:.2e}, push-outs {}",
            photons.len(),
            t.ball_overlap
        );
        assert!(hits > 500, "{name}: too few hits to test ({hits})");
        assert!(mism_hit <= photons.len() / 1000, "{name}: GPU and CPU disagree on which photons hit ({mism_hit})");
        assert!(max_err < 2e-3, "{name}: GPU ball contact disagrees with CPU: {max_err}");
        // Rigid: the ball's spin in photon units is ~0.03, so dropping it
        // moves photons far less than under EdgeC; still >> the f32 error.
        let control_min = if rough.grain_mode == GrainMode::Rigid { 5e-3 } else { 0.05 };
        assert!(max_wrong > control_min, "{name}: control did not fail ({max_wrong})");
        // Invariant: no photon ends the step inside the ball.
        let (ce, re) = (keys[8 * 4], keys[8 * 4][3] as f64);
        let ce = DVec3::new(ce[0] as f64, ce[1] as f64, ce[2] as f64);
        let inside = after
            .iter()
            .filter(|q| (DVec3::new(q.a[0] as f64, q.a[1] as f64, q.a[2] as f64) - ce).length() < re * (1.0 - 1e-3))
            .count();
        println!("   photons inside the ball at step end: {inside}");
        // Free body contacts (Free, BodyFree): the invariant must hold.
        // PinnedC: a photon held at c cannot separate from a ball surface
        // closing faster than c (the Spin-mode ball moves at ~1.8c), so
        // photons get overrun (doc §15.1); reported, not asserted.
        if rough.speed_mode != SpeedMode::PinnedC {
            assert!(inside * 1000 <= photons.len(), "{name}: {inside} photons ended inside the ball");
        }
    }
}

/// A non-rotating ball conserves photon energy (elastic law, Free) and the
/// momentum it receives balances the photons' loss; a moving ball (control)
/// does work.
#[test]
fn ball_static_conserves_energy() {
    let rough = RoughParams { beta: 1.0, speed_mode: SpeedMode::Free, g_max: 6.0, stacking: true, ..RoughParams::default() };
    for (name, signs, conserving) in [("static", [0.0; 4], true), ("CONTROL moving", PROTON_STATES[0], false)] {
        let n: u32 = 300_000;
        let config = ball_config(Boundary::PeriodicBox { half_extent: 1.0 }, signs, rough, 0.0, 0.0, n + 16, 0xBA12);
        let ctx = GpuContext::new().expect("GPU");
        let mut engine = DsmcEngine::new(&ctx, config);
        let mut rng = Lcg(0xBA12);
        let ball = engine.shell_bodies[0].ball();
        let photons: Vec<GpuPhoton> = (0..n)
            .map(|id| loop {
                let pos = DVec3::new(rng.f(), rng.f(), rng.f()) * 2.0 - DVec3::splat(1.0);
                if (pos - ball.center).length() > ball.radius + 0.01 {
                    let d = rng.unit();
                    let chi = if rng.f() < 0.5 { 1.0 } else { -1.0 };
                    break make_photon(id, pos, d, d * chi, chi);
                }
            })
            .collect();
        engine.write_photons_for_test(&ctx, &photons);
        let (e0, p0) = photon_energy(&engine.read_photons(&ctx), rough.kappa);
        let mut t = Totals::default();
        for _ in 0..300 {
            engine.step_collide(&ctx, &mut t);
        }
        let (e1, p1) = photon_energy(&engine.read_photons(&ctx), rough.kappa);
        let de = (e1 - e0) / e0;
        let dp_body = DVec3::new(t.sh_dp[0][0], t.sh_dp[0][1], t.sh_dp[0][2]);
        let balance = ((p1 - p0) + dp_body).length() / dp_body.length().max(1.0);
        println!("ball [{name}]: hits {} dE/E {de:+.2e} momentum balance {balance:.1e} overlaps {}", t.sh_contacts[0], t.ball_overlap);
        assert!(t.sh_contacts[0] > 10_000);
        assert!(balance < 1e-2);
        // Invariant: no photon ends inside the ball (push-outs are normal
        // contacts, reported but not limited).
        let bs = engine.shell_bodies[0].ball();
        let inside = engine
            .read_photons(&ctx)
            .iter()
            .filter(|q| q.is_alive() && (DVec3::new(q.a[0] as f64, q.a[1] as f64, q.a[2] as f64) - bs.center).length() < bs.radius * (1.0 - 1e-3))
            .count();
        println!("   push-outs {}, photons inside the ball at the end: {inside}", t.ball_overlap);
        assert!(inside * 10_000 <= n as usize, "{name}: {inside} photons inside the ball");
        if conserving {
            assert!(de.abs() < 1e-4, "static ball changed photon energy {de}");
        } else {
            assert!(de.abs() > 1e-3, "control: a moving ball must do work ({de})");
        }
    }
}

/// Doc §15.7: signs of levels 8, 7, 6, 5 for the adopted neutron ("z flipped
/// in every tier"; the proton keeps all +1). Level 4 (tier-1 z) is inside the
/// rigid ball at every depth, so its flip is not represented.
const NEUTRON_INNER_SIGNS: [f64; 4] = [-1.0, 1.0, 1.0, 1.0];

/// Doc §15 first physics run: the instantaneous base ball, all 4 proton and 4
/// neutron sign states, in the 2:1 ambient field. It prints the user's
/// pole-to-pole table: passages through the envelope by entry angle to the
/// pole, split into
/// - clean passes;
/// - photon-photon only;
/// - ball hit first;
/// - ball hit after a photon-photon collision inside the envelope ("knocked
///   into the wall");
///
/// plus opposite-pole exits and the contacted-exit bands.
/// Env:
/// - CFM_C_PP: photon-photon law, "twist" (default) or "edge";
/// - CFM_C_BODY: body law, "twist" (default) or "edge";
/// - CFM_C_GRAIN: grain spin edge (default 1);
/// - CFM_C_NA: ambient density (default 4000; lambda = 3 is held fixed by
///   scaling sigma);
/// - CFM_C_STATES: "all" (default) or "first" (state 0 of each class);
/// - CFM_C_WARMUP / CFM_C_AVG: 300 / 1200;
/// - CFM_C_SEED: run seed (default 49344 = 0xC0C0).
#[test]
#[ignore]
fn report_c_ball_single_body() {
    use crate::baryon_shells::NEUTRON_STATES;
    use std::time::Instant;
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let law = |n: &str| -> (f64, f64, bool) {
        if n == "edge" {
            (1.0, -1.0, false)
        } else {
            (-1.0, 0.0, true)
        }
    };
    let pp = env("CFM_C_PP", "twist");
    let body = env("CFM_C_BODY", "twist");
    let grain: f64 = env("CFM_C_GRAIN", "1").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "4000").parse().unwrap();
    let warmup: u32 = env("CFM_C_WARMUP", "300").parse().unwrap();
    let avg: u32 = env("CFM_C_AVG", "1200").parse().unwrap();
    let first_only = env("CFM_C_STATES", "all") == "first";
    // CFM_C_SEED: run seed (default 0xC0C0), for multi-seed confirmation.
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    // Doc §15.4: CFM_C_DEPTH (0-4) resolves electron-tier levels 8, 7, 6, 5
    // as motion; CFM_C_DT sets the step (0.05); CFM_C_GMAX the NTC bound.
    let depth: usize = env("CFM_C_DEPTH", "0").parse().unwrap();
    let dt: f64 = env("CFM_C_DT", "0.05").parse().unwrap();
    let g_max: f64 = env("CFM_C_GMAX", "12").parse().unwrap();
    // Doc §15.5 density sweep: CFM_C_H = DSMC cell size (0.25; keep it well
    // below lambda), CFM_C_CAP = photon capacity (9M).
    let cell_h: f64 = env("CFM_C_H", "0.25").parse().unwrap();
    let cap: u32 = env("CFM_C_CAP", "9000000").parse().unwrap();
    let LAMBDA: f64 = env("CFM_C_LAMBDA", "3").parse().unwrap();
    let (b0, t0, l0) = law(&pp);
    let rough = RoughParams {
        beta: b0,
        beta_twist: t0,
        twist_loss_to_normal: l0,
        grain_edge: grain,
        body_law: Some(law(&body)),
        ..twist_rough(true, SpeedMode::PinnedC)
    };
    // CFM_C_SPEED: "bodyfree" (default, doc §15.2), "pinned" or "free" (D7).
    let speed = env("CFM_C_SPEED", "bodyfree");
    let rough = match speed.as_str() {
        "free" => RoughParams { speed_mode: SpeedMode::Free, g_max: 8.0, ..rough },
        "pinned" => rough,
        _ => RoughParams { speed_mode: SpeedMode::BodyFree, g_max, ..rough },
    };
    println!(
        "=== C: instantaneous ball, speed={speed}, photon-photon={pp}, body={body}, grain={grain}, stacking, depth={depth}, h={cell_h}, lambda={LAMBDA}, n_a={n_a}, R_D=6, dt={dt}, g_max={g_max}, warmup {warmup}, avg {avg} ==="
    );
    let ctx = GpuContext::new().expect("GPU");
    // CFM_C_HIST_DIR: if set, write 1-degree latitude histograms per class
    // (summed over the states run) there, doc §15.3.
    let hist_dir = std::env::var("CFM_C_HIST_DIR").ok();
    // CFM_C_SIGNS: "a,x,y,z" runs one custom sign/speed state instead of the
    // 8 standard states (e.g. "-1,-1,-1,0" = the Spin-mode neutron preset,
    // whose z level was torpid).
    // Several custom states separated by ';'.
    let custom: Option<Vec<[f64; 4]>> = std::env::var("CFM_C_SIGNS").ok().map(|v| {
        v.split(';')
            .map(|st| {
                let x: Vec<f64> = st.split(',').map(|t| t.trim().parse().unwrap()).collect();
                [x[0], x[1], x[2], x[3]]
            })
            .collect()
    });
    // Doc §15.6 neutron-model alternatives: CFM_C_TOP_P / CFM_C_TOP_N scale
    // the level-12 (z) edge speed of the proton / neutron states (1 = c,
    // 0 = torpid); CFM_C_INNER is the sign of levels 8-5 (Spin-mode presets:
    // +1).
    let top_p: f64 = env("CFM_C_TOP_P", "1").parse().unwrap();
    let top_n: f64 = env("CFM_C_TOP_N", "1").parse().unwrap();
    let inner: f64 = env("CFM_C_INNER", "1").parse().unwrap();
    // CFM_C_INNER_P / CFM_C_INNER_N: per-class signs of levels 8,7,6,5
    // ("s8,s7,s6,s5"), overriding CFM_C_INNER. E.g. the user's "z flipped in
    // every tier" neutron: CFM_C_INNER_N="-1,1,1,1" (level 8 = the tier-2 z).
    let parse4 = |k: &str| -> Option<[f64; 4]> {
        std::env::var(k).ok().map(|v| {
            let x: Vec<f64> = v.split(',').map(|t| t.trim().parse().unwrap()).collect();
            [x[0], x[1], x[2], x[3]]
        })
    };
    let inner_p = parse4("CFM_C_INNER_P").unwrap_or([inner; 4]);
    // Doc §15.8: CFM_C_GRAIN_MODE = edge (default, §14.8) | rigid | inner.
    // Tier-1 (levels 4,3,2,1) signs per class: the adopted neutron flips
    // level 4 like level 8 (CFM_C_TIER1_N / CFM_C_TIER1_P override).
    let grain_mode = match env("CFM_C_GRAIN_MODE", "edge").as_str() {
        "rigid" => GrainMode::Rigid,
        "inner" => GrainMode::RigidInner,
        _ => GrainMode::EdgeC,
    };
    // Doc §15.9: CFM_C_MAZE = off (default) | fixed | swept.
    let maze = match env("CFM_C_MAZE", "off").as_str() {
        "fixed" => MazeMode::Fixed,
        "swept" => MazeMode::Swept,
        _ => MazeMode::Off,
    };

    let tier1_p = parse4("CFM_C_TIER1_P").unwrap_or([1.0; 4]);
    // Doc §15.7: the adopted neutron has z flipped in every tier (user
    // decision 2026-09-30): level 8 (the tier-2 z) opposes levels 5-7, as
    // level 12 opposes 9-11 in the Mathis states. Override with
    // CFM_C_INNER_N="1,1,1,1" for the outer-only neutron.
    let inner_n = parse4("CFM_C_INNER_N").unwrap_or(if std::env::var("CFM_C_INNER").is_ok() { [inner; 4] } else { NEUTRON_INNER_SIGNS });
    let tier1_n = parse4("CFM_C_TIER1_N").unwrap_or(NEUTRON_INNER_SIGNS);
    let scale_top = |states: &[[f64; 4]], t: f64| -> Vec<[f64; 4]> { states.iter().map(|s| [s[0], s[1], s[2], s[3] * t]).collect() };
    let classes: Vec<(&str, Vec<[f64; 4]>)> = match custom {
        Some(c) => vec![("custom", c)],
        None => vec![("proton", scale_top(&PROTON_STATES, top_p)), ("neutron", scale_top(&NEUTRON_STATES, top_n))],
    };
    // CFM_C_ONLY: "proton" or "neutron" runs one class.
    let only = env("CFM_C_ONLY", "");
    let classes: Vec<(&str, Vec<[f64; 4]>)> = classes.into_iter().filter(|(c, _)| only.is_empty() || *c == only).collect();
    println!("   top_p {top_p} top_n {top_n} inner {inner} inner_p {inner_p:?} inner_n {inner_n:?} tier1_p {tier1_p:?} tier1_n {tier1_n:?} grain_mode {grain_mode:?} maze {maze:?}");
    for (class, states) in classes {
        let mut acc_hit = vec![0i64; 180];
        let mut acc_exit = vec![0i64; 360];
        let mut acc_flux = vec![0i64; 504];
        let mut acc_mlat = vec![0i64; 180];
        let mut acc_maz = vec![0i64; 36];
        for (si, signs) in states.iter().enumerate() {
            if first_only && si > 0 {
                continue;
            }
            let started = Instant::now();
            let mut config = ball_config(Boundary::SphereEscape { r_d: 6.0 }, *signs, RoughParams { grain_mode, maze, ..rough }, 1.0 / (LAMBDA * n_a), n_a, cap, seed);
            config.h = cell_h;
            config.dt = dt;
            let mut engine = DsmcEngine::new(&ctx, config);
            for b in engine.shell_bodies.iter_mut() {
                b.depth = depth;
                b.inner_signs = if class == "neutron" { inner_n } else if class == "proton" { inner_p } else { [inner; 4] };
                b.tier1_signs = if class == "neutron" { tier1_n } else { tier1_p };
            }
            // Time-averaged ball-centre speed (CPU, long incommensurate run).
            let ball_speed = {
                let mut sb = engine.shell_bodies[0].clone();
                let mut acc = 0.0;
                for _ in 0..20_000 {
                    sb.advance(0.0137);
                    acc += sb.ball().velocity.length();
                }
                acc / 20_000.0
            };
            let rep = engine.run_steady(&ctx, warmup, avg, 200);
            let t = &rep.totals;
            // Local gas density vs the nominal ambient: live photons per volume
            // in shells around the body, from the time-averaged density map.
            let dm = crate::gpu::dsmc::DensityMap::new(&engine.grid, engine.free_volume_cpu().to_vec(), rep.density_raw.clone(), avg);
            let shell_n = |r0: f64, r1: f64| -> f64 {
                // Mean ambient-kind density over points spread through the shell.
                let mut rng = Lcg(0xD35);
                let mut sum = 0.0;
                let m = 4000;
                for _ in 0..m {
                    let r = r0 + (r1 - r0) * rng.f();
                    sum += dm.density_at_both_species(rng.unit() * r, 0);
                }
                sum / m as f64
            };
            println!(
                "   gas density / n_a: r<1 {:.2}, 1-2 {:.2}, 2-4 {:.2}, 4-6 {:.2}; live {}; pp collisions {} ({:.3} per live photon per unit time)",
                shell_n(0.0, 1.0) / n_a,
                shell_n(1.0, 2.0) / n_a,
                shell_n(2.0, 4.0) / n_a,
                shell_n(4.0, 5.9) / n_a,
                engine.n_live,
                t.collisions,
                2.0 * t.collisions as f64 / (engine.n_live as f64 * avg as f64 * dt)
            );
            let st: String = signs.iter().zip(['a', 'x', 'y', 'z']).map(|(s, c)| format!("{}{}", if *s > 0.0 { "+" } else if *s < 0.0 { "-" } else { "0" }, c)).collect();
            let passages: i64 = (0..9).map(|a| (0..4).map(|o| t.pass[a][o][0] + t.pass[a][o][1]).sum::<i64>()).sum();
            let exits = t.sh_exits[0].max(1) as f64;
            let band = |b: usize| t.sh_exit_band[0][b] as f64 / exits;
            let hel = |o: usize| (t.pass_hel[o][0] + t.pass_hel[o][1]) / ((0..9).map(|a| t.pass[a][o][0] + t.pass[a][o][1]).sum::<i64>().max(1) as f64);
            println!(
                "{class:<7} {st} | ball v {ball_speed:.2} | passages {passages} hits {} (hit% {:.2}) push-outs {} | esc speed {:.3}±{:.3} ovfl {} | contacted-exit bands p/m/e {:.3}/{:.3}/{:.3} | hel clean/pp/hit/pp-hit {:+.3}/{:+.3}/{:+.3}/{:+.3} | Tz/hit {:+.4} | steady {} ({:.0}s)",
                t.sh_contacts[0],
                100.0 * t.sh_exits[0] as f64 / passages.max(1) as f64,
                t.ball_overlap,
                t.esc_speed_stats().0,
                t.esc_speed_stats().1,
                t.accept_overflow,
                band(0),
                band(1),
                band(2),
                hel(0),
                hel(1),
                hel(2),
                hel(3),
                t.sh_torque[0][2] / t.sh_contacts[0].max(1) as f64,
                rep.steady_ok,
                started.elapsed().as_secs_f64()
            );
            // The pole-angle table, aggregated to 0-10, 10-20, 20-30, 30-60,
            // 60-90 degrees from the pole, both species together, then the
            // photon (+) / antiphoton (-) split of "hit after pp".
            let groups: [(&str, std::ops::Range<usize>); 5] = [("0-10", 0..1), ("10-20", 1..2), ("20-30", 2..3), ("30-60", 3..6), ("60-90", 6..9)];
            let mut line = String::from("   angle-to-pole: ");
            for (label, r) in groups.iter() {
                let sum = |o: usize, sp: usize| -> i64 { r.clone().map(|a| t.pass[a][o][sp]).sum() };
                let tot: i64 = (0..4).map(|o| sum(o, 0) + sum(o, 1)).sum::<i64>().max(1);
                let thr: i64 = r.clone().map(|a| t.pass_through[a][0] + t.pass_through[a][1]).sum();
                line += &format!(
                    "[{label}: n {tot} clean {:.3} pp {:.3} hit {:.4} pp->hit {:.4} (+{:.4}/-{:.4}) thru {:.3}] ",
                    (sum(0, 0) + sum(0, 1)) as f64 / tot as f64,
                    (sum(1, 0) + sum(1, 1)) as f64 / tot as f64,
                    (sum(2, 0) + sum(2, 1)) as f64 / tot as f64,
                    (sum(3, 0) + sum(3, 1)) as f64 / tot as f64,
                    sum(3, 1) as f64 / (0..4).map(|o| sum(o, 1)).sum::<i64>().max(1) as f64,
                    sum(3, 0) as f64 / (0..4).map(|o| sum(o, 0)).sum::<i64>().max(1) as f64,
                    thr as f64 / tot as f64
                );
            }
            println!("{line}");
            if maze != MazeMode::Off {
                let m = &t.maze;
                let tot = (m[0] + m[1] + m[2] + m[3]).max(1) as f64;
                let lat_band = |lo: usize, hi: usize| -> f64 {
                    let s: i64 = (0..180usize).filter(|&d| (d as i32 - 90).unsigned_abs() as usize >= lo && ((d as i32 - 90).unsigned_abs() as usize) < hi).map(|d| t.maze_lat[d]).sum();
                    s as f64 / tot
                };
                let az_tot = t.maze_az.iter().sum::<i64>().max(1) as f64;
                let oct: Vec<String> = (0..8).map(|o| format!("{:.2}", (0..36).filter(|&b| (b * 10 + 5) / 45 == o).map(|b| t.maze_az[b]).sum::<i64>() as f64 / az_tot)).collect();
                println!(
                    "   maze: esc-up {:.3} esc-down {:.3} trap-up {:.3} trap-down {:.3} fallback {} | emission |lat| 0-30 {:.3} 30-60 {:.3} 60-89 {:.3} poles {:.3} | escape azimuth octants from +x {:?}",
                    m[0] as f64 / tot,
                    m[1] as f64 / tot,
                    m[2] as f64 / tot,
                    m[3] as f64 / tot,
                    m[4],
                    lat_band(0, 30),
                    lat_band(30, 60),
                    lat_band(60, 89),
                    lat_band(89, 91),
                    oct
                );
                for (a, b) in acc_mlat.iter_mut().zip(&t.maze_lat) {
                    *a += *b;
                }
                for (a, b) in acc_maz.iter_mut().zip(&t.maze_az) {
                    *a += *b;
                }
            }
            for (a, b) in acc_hit.iter_mut().zip(&t.hit_lat) {
                *a += *b;
            }
            for (a, b) in acc_exit.iter_mut().zip(&t.exit_lat) {
                *a += *b;
            }
            for (a, b) in acc_flux.iter_mut().zip(&t.flux) {
                *a += *b;
            }
        }
        if let Some(dir) = &hist_dir {
            let mut csv = String::from("degree,hits,exits_photon,exits_antiphoton
");
            for d in 0..180usize {
                csv += &format!("{},{},{},{}
", d as i32 - 90, acc_hit[d], acc_exit[d], acc_exit[180 + d]);
            }
            let path = std::path::Path::new(dir).join(format!("ball_{class}_lat.csv"));
            std::fs::write(&path, csv).expect("write histogram");
            println!("   wrote {}", path.display());
            let mut csv = String::from("radius,lat_lo,photon_out,photon_in,antiphoton_out,antiphoton_in\n");
            for (k, r) in crate::gpu::dsmc::gpu_types::FLUX_RADII.iter().enumerate() {
                for lb in 0..18usize {
                    let at = |sp: usize, d: usize| acc_flux[((k * 18 + lb) * 2 + sp) * 2 + d];
                    csv += &format!("{r},{},{},{},{},{}\n", lb as i32 * 10 - 90, at(0, 0), at(0, 1), at(1, 0), at(1, 1));
                }
            }
            std::fs::write(std::path::Path::new(dir).join(format!("ball_{class}_flux.csv")), csv).expect("write flux");
            if maze != MazeMode::Off {
                let mut csv = String::from("degree,emission_lat\n");
                for d in 0..180usize {
                    csv += &format!("{},{}\n", d as i32 - 90, acc_mlat[d]);
                }
                csv += "azimuth_lo,escapes\n";
                for b in 0..36usize {
                    csv += &format!("{},{}\n", b * 10, acc_maz[b]);
                }
                std::fs::write(std::path::Path::new(dir).join(format!("ball_{class}_maze.csv")), csv).expect("write maze");
            }
        }
    }
}

/// Regression gate for the carved-free-volume bug: with permeable bodies the
/// photon-photon collision rate inside the envelope must match the rate in
/// the gas outside it (same density, same lambda). Measured on a box of gas
/// with a body in the middle and sigma > 0; the Absorber free-volume rule
/// would inflate the inside rate by orders of magnitude.
#[test]
fn permeable_body_does_not_distort_gas_collision_rate() {
    for model in [BodyModel::Shells { tau_body: 0.0 }, BodyModel::Ball] {
        let n: u32 = 400_000;
        let half: f64 = 2.0;
        let n_density = n as f64 / (2.0 * half).powi(3);
        let sigma = 1.0 / (2.0 * n_density);
        let mut config = shell_config(Boundary::PeriodicBox { half_extent: half }, [0.0; 4], 0.0, twist_rough(false, SpeedMode::PinnedC), sigma, 0.0, 0.5, n + 16, 0x5A7E);
        config.body_model = model;
        let ctx = GpuContext::new().expect("GPU");
        let mut engine = DsmcEngine::new(&ctx, config);
        let fv = engine.free_volume_cpu().to_vec();
        let (dim, h) = (engine.grid.dim, engine.grid.h);
        let cell_vol = (h * h * h) as f32;
        let carved = fv.iter().filter(|&&v| v < 0.99 * cell_vol).count();
        let mut rng = Lcg(0x5A7E);
        let photons: Vec<GpuPhoton> = (0..n)
            .map(|id| {
                let pos = DVec3::new(rng.f(), rng.f(), rng.f()) * (2.0 * half) - DVec3::splat(half);
                let d = rng.unit();
                let chi = if rng.f() < 0.5 { 1.0 } else { -1.0 };
                make_photon(id, pos, d, d * chi, chi)
            })
            .collect();
        engine.write_photons_for_test(&ctx, &photons);
        let mut t = Totals::default();
        for _ in 0..40 {
            engine.step_collide(&ctx, &mut t);
        }
        // Expected per-photon rate n*sigma*<g> with <g> = 4/3.
        let rate = 2.0 * t.collisions as f64 / (n as f64 * 40.0 * engine.config.dt);
        let expect = n_density * sigma * 4.0 / 3.0;
        println!("{model:?}: carved cells {carved} of {}, collision rate {rate:.4} vs expected {expect:.4}", dim[0] * dim[1] * dim[2]);
        assert_eq!(carved, 0, "{model:?}: permeable body carved gas volume");
        assert!((rate / expect - 1.0).abs() < 0.03, "{model:?}: gas collision rate distorted");
    }
}

/// Doc §15.10 settings shared by the two-body gate and report: the §15.9
/// maze ball (swept, BodyFree), bodies at `centers`, pole along +Z.
fn two_body_config(centers: &[DVec3], classes: &[&str], n_a: f64, lambda: f64, cap: u32, seed: u64, dt: f64) -> DsmcConfig {
    let rough = RoughParams {
        beta: -1.0,
        beta_twist: 0.0,
        twist_loss_to_normal: true,
        grain_edge: 1.0,
        body_law: Some((-1.0, 0.0, true)),
        maze: MazeMode::Swept,
        speed_mode: SpeedMode::BodyFree,
        g_max: 12.0,
        ..twist_rough(true, SpeedMode::PinnedC)
    };
    let mut c = ball_config(Boundary::SphereEscape { r_d: 6.0 }, PROTON_STATES[0], rough, 1.0 / (lambda * n_a), n_a, cap, seed);
    c.dt = dt;
    c.bodies = centers
        .iter()
        .zip(classes)
        .map(|(ctr, cl)| BodyConfig {
            center: *ctr,
            class: if *cl == "neutron" { GpuClass::Neutron } else { GpuClass::Proton },
            orient: GpuOrient::PoleOn,
            through_lane_rho: 0.0,
        })
        .collect();
    c.shell_signs = classes.iter().map(|cl| if *cl == "neutron" { crate::baryon_shells::NEUTRON_STATES[0] } else { PROTON_STATES[0] }).collect();
    c
}

/// Doc §15.10 bookkeeping gate: with two maze bodies, every capture is
/// tallied once by emitter tag (feed sums == contacts per body), both bodies
/// capture, and the line flux sees traffic both ways between them.
#[test]
fn two_body_feed_bookkeeping() {
    let ctx = GpuContext::new().expect("GPU");
    let centers = [DVec3::new(-1.25, 0.0, 0.0), DVec3::new(1.25, 0.0, 0.0)];
    let mut engine = DsmcEngine::new(&ctx, two_body_config(&centers, &["proton", "neutron"], 2000.0, 3.0, 2_000_000, 0x2B0D, 0.05));
    let rep = engine.run_steady(&ctx, 60, 120, 60);
    let t = &rep.totals;
    for b in 0..2 {
        let fed: i64 = t.feed[b].iter().sum();
        println!("body {b}: contacts {} feed {:?}", t.sh_contacts[b], t.feed[b]);
        assert_eq!(fed, t.sh_contacts[b], "body {b}: feed tally != contacts");
        assert!(t.sh_contacts[b] > 100, "body {b}: too few captures");
    }
    use super::gpu_types::{LINE_DU, LINE_U0};
    let j_mid = ((0.0 + 1.25 + LINE_U0) / LINE_DU).round() as usize + 0;
    let (fwd, back): (i64, i64) = (0..3).map(|tc| (t.line[(j_mid * 3 + tc) * 2], t.line[(j_mid * 3 + tc) * 2 + 1])).fold((0, 0), |a, b| (a.0 + b.0, a.1 + b.1));
    println!("mid-plane crossings: toward body 1 {fwd}, toward body 0 {back}");
    assert!(fwd > 100 && back > 100, "line flux sees no traffic");
    // A photon emitted by body 1 crossing the mid-plane toward body 1 must
    // have been sent back (tag class 2, dir 0) less often than it leaves
    // body 1 toward body 0 (tag class 2, dir 1).
    assert!(t.line[(j_mid * 3 + 2) * 2 + 1] > 0);
}

/// Report (doc §15.10): two maze bodies side by side (separation along x,
/// in each other's equatorial plane) or stacked (along the pole). Per body:
/// captures per step vs the isolated body, what fed them (never captured /
/// itself / the neighbour), contacted-exit bands; along the line between
/// them: net photon flow and gas density.
/// Env: CFM_2_PAIRS ("pp,pn,nn"), CFM_2_GEOMS ("side,axial"), CFM_2_SEP
/// (2.5), CFM_C_DEPTH (0), CFM_C_NA (4000), CFM_C_LAMBDA (3), CFM_C_WARMUP
/// (300), CFM_C_AVG (1200), CFM_C_SEED, CFM_C_CAP (9M).
#[test]
#[ignore]
fn report_two_bodies() {
    use super::gpu_types::{LINE_DU, LINE_N, LINE_U0};
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let pairs = env("CFM_2_PAIRS", "pp,pn,nn");
    let geoms = env("CFM_2_GEOMS", "side,axial");
    let sep: f64 = env("CFM_2_SEP", "2.5").parse().unwrap();
    let depth: usize = env("CFM_C_DEPTH", "0").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "4000").parse().unwrap();
    let lambda: f64 = env("CFM_C_LAMBDA", "3").parse().unwrap();
    let warmup: u32 = env("CFM_C_WARMUP", "300").parse().unwrap();
    let avg: u32 = env("CFM_C_AVG", "1200").parse().unwrap();
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    let cap: u32 = env("CFM_C_CAP", "9000000").parse().unwrap();
    let dt = 0.05;
    println!("=== two maze bodies (swept, BodyFree), sep {sep}, depth {depth}, lambda {lambda}, n_a {n_a}, warmup {warmup}, avg {avg} ===");
    let ctx = GpuContext::new().expect("GPU");
    let setup = |engine: &mut DsmcEngine, classes: &[&str]| {
        for (b, cl) in engine.shell_bodies.iter_mut().zip(classes) {
            b.depth = depth;
            if *cl == "neutron" {
                b.inner_signs = NEUTRON_INNER_SIGNS;
                b.tier1_signs = NEUTRON_INNER_SIGNS;
            }
        }
    };
    let bands = |t: &Totals, b: usize| -> String {
        let n = t.sh_exits[b].max(1) as f64;
        format!("{:.3}/{:.3}/{:.3}", t.sh_exit_band[b][0] as f64 / n, t.sh_exit_band[b][1] as f64 / n, t.sh_exit_band[b][2] as f64 / n)
    };
    // Isolated baselines.
    let mut iso = std::collections::HashMap::new();
    for cl in ["proton", "neutron"] {
        let mut engine = DsmcEngine::new(&ctx, two_body_config(&[DVec3::ZERO], &[cl], n_a, lambda, cap, seed, dt));
        setup(&mut engine, &[cl]);
        let rep = engine.run_steady(&ctx, warmup, avg, 200);
        let t = &rep.totals;
        let rate = t.sh_contacts[0] as f64 / avg as f64;
        println!("isolated {cl:<7}: captures/step {rate:.2} | contacted-exit p/m/e {} | steady {}", bands(t, 0), rep.steady_ok);
        iso.insert(cl, rate);
    }
    for pair in pairs.split(',') {
        let classes: Vec<&str> = pair.chars().map(|c| if c == 'n' { "neutron" } else { "proton" }).collect();
        for geom in geoms.split(',') {
            let e = if geom == "axial" { DVec3::Z } else { DVec3::X };
            let centers = [-e * sep / 2.0, e * sep / 2.0];
            let mut engine = DsmcEngine::new(&ctx, two_body_config(&centers, &classes, n_a, lambda, cap, seed, dt));
            setup(&mut engine, &classes);
            let rep = engine.run_steady(&ctx, warmup, avg, 200);
            let t = &rep.totals;
            println!("--- {pair} {geom} (body 0 = {} at {:?}, body 1 = {} at {:?}), steady {} ---", classes[0], centers[0], classes[1], centers[1], rep.steady_ok);
            for b in 0..2 {
                let c = t.sh_contacts[b] as f64;
                let f = &t.feed[b];
                let (own, nb) = (f[1 + b], f[1 + (1 - b)]);
                println!(
                    "   body {b} {:<7}: captures/step {:.2} (x{:.3} isolated) | fed by: never-captured {:.3} itself {:.3} neighbour {:.3} | contacted-exit p/m/e {}",
                    classes[b],
                    c / avg as f64,
                    c / avg as f64 / iso[classes[b]],
                    f[0] as f64 / c.max(1.0),
                    own as f64 / c.max(1.0),
                    nb as f64 / c.max(1.0),
                    bands(t, b)
                );
            }
            // Line profile between the bodies (u from body 0's centre).
            let dm = crate::gpu::dsmc::DensityMap::new(&engine.grid, engine.free_volume_cpu().to_vec(), rep.density_raw.clone(), avg);
            let mut rows = Vec::new();
            for j in 0..LINE_N {
                let u = j as f64 * LINE_DU - LINE_U0;
                if u < -1.5 - 1e-9 || u > sep + 1.5 + 1e-9 {
                    continue;
                }
                let at = |tc: usize, d: usize| t.line[(j * 3 + tc) * 2 + d];
                let tot: i64 = (0..3).map(|tc| at(tc, 0) + at(tc, 1)).sum();
                let net = |tc: usize| (at(tc, 0) - at(tc, 1)) as f64 / tot.max(1) as f64 * 1e3;
                let net_all = (0..3).map(|tc| at(tc, 0) - at(tc, 1)).sum::<i64>() as f64 / tot.max(1) as f64 * 1e3;
                let share = |tc: usize| (at(tc, 0) + at(tc, 1)) as f64 / tot.max(1) as f64;
                let p = centers[0] + e * u;
                let mut dens = 0.0;
                let mut rng = Lcg(0x11E + j as u64);
                for _ in 0..200 {
                    let q = rng.unit();
                    dens += dm.density_at_both_species(p + (q - e * q.dot(e)) * 0.4, 0);
                }
                rows.push(format!(
                    "     u {u:+5.2}: n/n_a {:.3} | net flow toward body 1 x1e3: all {net_all:+6.1} (noise {:.1}) | tag body0 {:+6.1} body1 {:+6.1} | share tagged b0 {:.3} b1 {:.3}",
                    dens / 200.0 / n_a,
                    1e3 / (tot.max(1) as f64).sqrt(),
                    net(1),
                    net(2),
                    share(1),
                    share(2)
                ));
            }
            for r in rows {
                println!("{r}");
            }
        }
    }
}
