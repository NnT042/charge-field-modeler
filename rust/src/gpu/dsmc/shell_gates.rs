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
    // Doc §15.11: CFM_2_INTAKE = ball (default) | env<p> | polar<a>.
    let intake_s = env("CFM_2_INTAKE", "ball");
    let intake = if let Some(p) = intake_s.strip_prefix("env") {
        super::config::Intake::Envelope { p: p.parse().unwrap() }
    } else if let Some(a) = intake_s.strip_prefix("polar") {
        super::config::Intake::Polar { aperture: a.parse().unwrap() }
    } else {
        super::config::Intake::Ball
    };
    let two_body_config = |c: &[DVec3], cl: &[&str], n_a: f64, lambda: f64, cap: u32, seed: u64, dt: f64| -> DsmcConfig {
        let mut cfg = two_body_config(c, cl, n_a, lambda, cap, seed, dt);
        cfg.rough.intake = intake;
        cfg
    };
    println!("=== two maze bodies (swept, BodyFree), intake {intake:?}, sep {sep}, depth {depth}, lambda {lambda}, n_a {n_a}, warmup {warmup}, avg {avg} ===");
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

/// Doc §15.11 bookkeeping gate for envelope intake and the holding tank:
/// with p = 1 every photon entering the envelope is taken in (no clean
/// passages); held photons are neither lost nor collided (tank in - out ==
/// the host count == photons flagged held); a cap holds the tank near it;
/// the isotropic control emits isotropically.
#[test]
fn intake_tank_bookkeeping() {
    use super::config::Intake;
    let ctx = GpuContext::new().expect("GPU");
    for (name, maze, cap) in [("swept, tank 0.5, no cap", MazeMode::Swept, 0u32), ("swept, tank 0.5, cap 40", MazeMode::Swept, 40), ("isotropic, tank 0.5", MazeMode::Isotropic, 0)] {
        let mut config = two_body_config(&[DVec3::ZERO], &["proton"], 2000.0, 3.0, 2_000_000, 0x1A7E, 0.05);
        config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, hold_time: 0.5, hold_cap: cap, maze, ..config.rough };
        let mut engine = DsmcEngine::new(&ctx, config);
        let mut t = Totals::default();
        let mut max_held = 0i64;
        for _ in 0..240 {
            let st = engine.step_flux(&ctx);
            max_held = max_held.max(engine.held[0]);
            t.add(&st);
        }
        let flagged = engine.read_photons(&ctx).iter().filter(|q| q.c[0] & FLAG_ALIVE != 0 && q.c[0] & (1 << 20) != 0).count() as i64;
        let (tin, tout) = (t.hold[0][0], t.hold[0][1]);
        let em = &t.maze_lat;
        let tot: i64 = em.iter().sum();
        let eq: i64 = (60..120).map(|d| em[d]).sum();
        println!(
            "[{name}] env captures {} ball contacts {} clean passages {} | tank in {tin} out {tout} host held {} flagged {flagged} max held {max_held} | emission |lat|<30 {:.3}",
            t.env_capt[0],
            t.sh_contacts[0] - t.env_capt[0],
            t.sh_clean[0],
            engine.held[0],
            eq as f64 / tot.max(1) as f64
        );
        assert!(t.env_capt[0] > 1000, "{name}: too few captures");
        if cap == 0 {
            assert_eq!(t.sh_clean[0], 0, "{name}: a photon crossed the envelope without being taken in");
        }
        assert_eq!(tin - tout, engine.held[0], "{name}: tank tally != host count");
        assert_eq!(flagged, engine.held[0], "{name}: held photons lost or leaked");
        // The cap is soft: the host closes the intake one step after the
        // count passes it, so a step's intake (a few hundred here) can land
        // on top.
        if cap > 0 {
            let per_step = 250;
            assert!(max_held <= cap as i64 + 3 * per_step, "{name}: cap not respected ({max_held})");
            assert!(t.sh_clean[0] > 0, "{name}: a full tank should turn photons away");
        }
        let f = eq as f64 / tot.max(1) as f64;
        if maze == MazeMode::Isotropic {
            assert!((f - 0.5).abs() < 0.04, "{name}: isotropic emission not isotropic ({f})");
        } else {
            assert!(f > 0.99, "{name}: maze emission left the disc ({f})");
        }
    }
}

/// Doc §15.14 gate: under MazeSpeed::EjectC every maze emission leaves at
/// c relative to the ball surface, so no photon can outrun the fastest
/// surface point by more than c (the ball's centre speed is ~2c, its spin
/// adds at most ~|w| r_b); under Relative, repeated recaptures in a dense
/// field random-walk speeds upward (Fermi-style) past that bound.
#[test]
fn maze_eject_c_bounds_speed() {
    use super::config::{Intake, MazeSpeed};
    let ctx = GpuContext::new().expect("GPU");
    let mut max_by = Vec::new();
    for speed in [MazeSpeed::EjectC, MazeSpeed::Relative] {
        let mut config = two_body_config(&[DVec3::ZERO], &["proton"], 2000.0, 0.1, 2_000_000, 0xE7EC, 0.02);
        config.boundary = Boundary::SphereEscape { r_d: 2.0 };
        config.h = 0.1;
        config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, maze_speed: speed, ..config.rough };
        let mut engine = DsmcEngine::new(&ctx, config);
        for _ in 0..1500 {
            engine.step_flux(&ctx);
        }
        let mx = engine.read_photons(&ctx).iter().filter(|q| q.c[0] & FLAG_ALIVE != 0).map(|q| q.d[3]).fold(0.0f32, f32::max);
        println!("{speed:?}: max photon speed {mx:.2}");
        max_by.push(mx);
    }
    assert!(max_by[0] < 4.0, "EjectC: a photon outran surface + c ({})", max_by[0]);
    assert!(max_by[1] > max_by[0], "Relative should heat past EjectC ({} vs {})", max_by[1], max_by[0]);
}

/// Doc §15.24 gate: the through-charge hole. With hole 0 nothing passes;
/// with hole 1/3 the passes are the share of envelope entries that land on
/// the two polar caps inside the cylinder (1 - sqrt(8/9) = 0.057 of the
/// area for isotropic arrivals; the body's own emission skews it a little).
#[test]
fn through_hole_passes_polar_entries() {
    use super::config::{Intake, MazeSpeed};
    let ctx = GpuContext::new().expect("GPU");
    let mut share = Vec::new();
    for hole in [0.0, 1.0 / 3.0] {
        let mut config = two_body_config(&[DVec3::ZERO], &["proton"], 2000.0, 1.0, 2_000_000, 0x401E, 0.05);
        config.boundary = Boundary::SphereEscape { r_d: 2.5 };
        config.h = 0.25;
        config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, maze_speed: MazeSpeed::EjectC, through_hole: hole, ..config.rough };
        let mut engine = DsmcEngine::new(&ctx, config);
        for _ in 0..200 {
            engine.step_flux(&ctx);
        }
        let mut t = Totals::default();
        for _ in 0..400 {
            t.add(&engine.step_flux(&ctx));
        }
        let s = t.hole_pass[0] as f64 / (t.hole_pass[0] + t.env_capt[0]).max(1) as f64;
        println!("hole {hole:.3}: passes {}, envelope captures {}, pass share {s:.4}", t.hole_pass[0], t.env_capt[0]);
        share.push(s);
    }
    assert_eq!(share[0], 0.0, "no hole, no passes");
    assert!(share[1] > 0.03 && share[1] < 0.10, "hole 1/3 pass share {} (expected ~0.057)", share[1]);
}

/// Doc §15.17 gate: `set_body_centers` moves a body without moving the
/// domain. After a 1.0 shift the ball keyframes sit at the new centre, the
/// body keeps capturing at its old rate (uniform ambient), and the live
/// photon count (set by the fixed escape sphere) is unchanged.
#[test]
fn moving_body_follows() {
    use super::config::Intake;
    let ctx = GpuContext::new().expect("GPU");
    let mut config = two_body_config(&[DVec3::ZERO], &["proton"], 2000.0, 1.0, 2_000_000, 0x40BE, 0.05);
    config.boundary = Boundary::SphereEscape { r_d: 4.0 };
    config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, ..config.rough };
    let mut engine = DsmcEngine::new(&ctx, config);
    for _ in 0..150 {
        engine.step_flux(&ctx);
    }
    let rate = |engine: &mut DsmcEngine| -> (f64, u32) {
        let mut t = Totals::default();
        for _ in 0..200 {
            t.add(&engine.step_flux(&ctx));
        }
        (t.env_capt[0] as f64 / 200.0, engine.n_live)
    };
    let (r0, live0) = rate(&mut engine);
    let to = DVec3::new(1.0, 0.0, 0.0);
    engine.set_body_centers(&ctx, &[to]);
    for _ in 0..100 {
        engine.step_flux(&ctx);
    }
    let keys = engine.ball_keyframes();
    let ball = DVec3::new(keys[0][0] as f64, keys[0][1] as f64, keys[0][2] as f64);
    let (r1, live1) = rate(&mut engine);
    println!("captures/step {r0:.1} -> {r1:.1}, live {live0} -> {live1}, ball centre {ball:?}");
    // The ball orbits inside the envelope (levels 9-12), up to ~0.9 from the
    // body centre; before the move it was within that of the origin.
    assert!((ball - to).length() < 0.95, "ball keyframe did not follow ({ball:?})");
    assert_eq!(engine.shell_bodies[0].center, to);
    assert!((r1 / r0 - 1.0).abs() < 0.1, "capture rate changed after the move ({r0} -> {r1})");
    assert!((live1 as f64 / live0 as f64 - 1.0).abs() < 0.02, "domain moved ({live0} -> {live1})");
}

/// Doc §15.20 gate: `set_body_poses` turns a body. A proton turned so its
/// spin axis points along +x emits its maze disc about the new equator
/// (|d.x| small), not the old one; the unturned control does the reverse.
#[test]
fn turned_body_emits_in_its_own_frame() {
    use super::config::Intake;
    let ctx = GpuContext::new().expect("GPU");
    let mut res = Vec::new();
    for turned in [false, true] {
        let mut config = two_body_config(&[DVec3::ZERO], &["proton"], 2000.0, 3.0, 2_000_000, 0x7E2D, 0.05);
        config.boundary = Boundary::SphereEscape { r_d: 4.0 };
        config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, ..config.rough };
        let mut engine = DsmcEngine::new(&ctx, config);
        if turned {
            let q = DQuat::from_rotation_arc(DVec3::Z, DVec3::X);
            engine.set_body_poses(&ctx, &[DVec3::ZERO], Some(&[q]));
        }
        for _ in 0..120 {
            engine.step_flux(&ctx);
        }
        // Photons last emitted by body 0, still near it.
        let ph = engine.read_photons(&ctx);
        let near: Vec<DVec3> = ph
            .iter()
            .filter(|q| q.c[0] & FLAG_ALIVE != 0 && (q.c[0] >> 17) & 7 == 1)
            .filter(|q| DVec3::new(q.a[0] as f64, q.a[1] as f64, q.a[2] as f64).length() < 1.5)
            .map(|q| DVec3::new(q.b[0] as f64, q.b[1] as f64, q.b[2] as f64))
            .collect();
        let n = near.len().max(1) as f64;
        let eq_z = near.iter().filter(|d| d.z.abs() < 0.5).count() as f64 / n;
        let eq_x = near.iter().filter(|d| d.x.abs() < 0.5).count() as f64 / n;
        println!("turned {turned}: {} recent emissions, within 30 deg of the z-equator {eq_z:.3}, of the x-equator {eq_x:.3}", near.len());
        assert!(near.len() > 200, "too few emitted photons to judge");
        res.push((eq_z, eq_x));
    }
    assert!(res[0].0 > res[0].1 + 0.15, "unturned body: disc not about the z-equator ({:?})", res[0]);
    assert!(res[1].1 > res[1].0 + 0.15, "turned body: disc did not follow the axis ({:?})", res[1]);
}

/// Report (doc §15.11): one body with envelope / polar / ball intake, tank,
/// maze or isotropic emission. The option-2 test: with intake equally likely
/// everywhere on the envelope, does the latitude of the never-captured
/// photons' intake focus toward the poles when the body emits a disc,
/// compared with the isotropic-emission control and the no-field
/// expectation? Env: CFM_I_RUNS = ';'-separated "intake:emission:hold",
/// intake = ball | env<p> | polar<a>, emission = swept | iso; CFM_I_CLASS
/// (proton), CFM_C_LAMBDA (3), CFM_C_NA (4000), CFM_C_DEPTH (0),
/// CFM_C_WARMUP (300), CFM_C_AVG (1200), CFM_C_SEED, CFM_C_CAP.
#[test]
#[ignore]
fn report_intake() {
    use super::config::Intake;
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let runs = env("CFM_I_RUNS", "env1:swept:0;env1:iso:0");
    let class = env("CFM_I_CLASS", "proton");
    let lambda: f64 = env("CFM_C_LAMBDA", "3").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "4000").parse().unwrap();
    let depth: usize = env("CFM_C_DEPTH", "0").parse().unwrap();
    let warmup: u32 = env("CFM_C_WARMUP", "300").parse().unwrap();
    let avg: u32 = env("CFM_C_AVG", "1200").parse().unwrap();
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    let cap: u32 = env("CFM_C_CAP", "9000000").parse().unwrap();
    // CFM_C_DT (0.05; keep below lambda), CFM_C_RD (escape radius, 6),
    // CFM_C_BLOCK (steps per printed block of the recaptured share, 0 = off).
    let dt: f64 = env("CFM_C_DT", "0.05").parse().unwrap();
    let r_d: f64 = env("CFM_C_RD", "6").parse().unwrap();
    let block: u32 = env("CFM_C_BLOCK", "0").parse().unwrap();
    println!("=== intake report: {class}, dt {dt}, r_d {r_d}, lambda {lambda}, n_a {n_a}, depth {depth}, warmup {warmup}, avg {avg} ===");
    println!("   (no-field expectation for intake latitude on a sphere, isotropic arrivals: |lat| 0-30 0.500, 30-60 0.366, 60-90 0.134)");
    let ctx = GpuContext::new().expect("GPU");
    for run in runs.split(';') {
        let f: Vec<&str> = run.split(':').collect();
        let intake = if let Some(p) = f[0].strip_prefix("env") {
            Intake::Envelope { p: p.parse().unwrap() }
        } else if let Some(a) = f[0].strip_prefix("polar") {
            Intake::Polar { aperture: a.parse().unwrap() }
        } else {
            Intake::Ball
        };
        let maze = if f[1] == "iso" { MazeMode::Isotropic } else { MazeMode::Swept };
        let hold: f64 = f.get(2).map(|h| h.parse().unwrap()).unwrap_or(0.0);
        let mut config = two_body_config(&[DVec3::ZERO], &[class.as_str()], n_a, lambda, cap, seed, dt);
        // CFM_I_SPEED: ejectc (default) | rel (doc §15.14).
        let maze_speed = if env("CFM_I_SPEED", "ejectc") == "ejectc" { super::config::MazeSpeed::EjectC } else { super::config::MazeSpeed::Relative };
        config.rough = RoughParams { intake, maze, maze_speed, hold_time: hold, ..config.rough };
        // CFM_C_H: DSMC cell size (0.25); keep it at or below lambda.
        config.h = env("CFM_C_H", "0.25").parse().unwrap();
        config.boundary = Boundary::SphereEscape { r_d };
        // CFM_I_PHI: photon (vs antiphoton) fraction of the ambient (2/3).
        if let Ok(phi) = std::env::var("CFM_I_PHI") {
            config.photon_fraction = phi.parse().unwrap();
        }
        let mut engine = DsmcEngine::new(&ctx, config);
        for b in engine.shell_bodies.iter_mut() {
            b.depth = depth;
            if class == "neutron" {
                b.inner_signs = NEUTRON_INNER_SIGNS;
                b.tier1_signs = NEUTRON_INNER_SIGNS;
            }
        }
        let t0 = std::time::Instant::now();
        for _ in 0..warmup {
            engine.step_flux(&ctx);
        }
        let mut t = Totals::default();
        let mut bt = Totals::default();
        let mut trace = String::new();
        for s in 0..avg {
            let st = engine.step_flux(&ctx);
            t.add(&st);
            bt.add(&st);
            if block > 0 && (s + 1) % block == 0 {
                let c = bt.feed[0].iter().sum::<i64>().max(1) as f64;
                trace += &format!(" {:.3}", 1.0 - bt.feed[0][0] as f64 / c);
                bt = Totals::default();
            }
        }
        if block > 0 {
            println!("   recaptured share by block of {block} steps:{trace}");
        }
        struct Rep { steady_ok: &'static str }
        let rep = Rep { steady_ok: "n/a" };
        println!("   wall {:.0} s, live {}", t0.elapsed().as_secs_f64(), engine.n_live);
        {
            // Speed census of the live free photons at the end (doc §15.14).
            let ph = engine.read_photons(&ctx);
            let mut sp: Vec<f64> = ph.iter().filter(|q| q.c[0] & FLAG_ALIVE != 0 && q.c[0] & (1 << 20) == 0).map(|q| q.d[3] as f64).collect();
            sp.sort_by(|a, b| a.partial_cmp(b).unwrap());
            let n = sp.len().max(1) as f64;
            let near: Vec<f64> = ph.iter().filter(|q| q.c[0] & FLAG_ALIVE != 0 && q.c[0] & (1 << 20) == 0 && DVec3::new(q.a[0] as f64, q.a[1] as f64, q.a[2] as f64).length() < 2.0).map(|q| q.d[3] as f64).collect();
            println!(
                "   speeds: mean {:.3} median {:.3} p99 {:.2} max {:.1} | >2c {:.4} >4c {:.5} | mean within r<2: {:.3} ({} photons)",
                sp.iter().sum::<f64>() / n,
                sp[sp.len() / 2],
                sp[(sp.len() * 99) / 100],
                sp.last().copied().unwrap_or(0.0),
                sp.iter().filter(|&&s| s > 2.0).count() as f64 / n,
                sp.iter().filter(|&&s| s > 4.0).count() as f64 / n,
                near.iter().sum::<f64>() / near.len().max(1) as f64,
                near.len()
            );
        }
        let t = &t;
        let caps = t.feed[0].iter().sum::<i64>();
        let bands = |lo: usize| -> [f64; 3] {
            let h = &t.capt_lat[lo..lo + 180];
            let tot = h.iter().sum::<i64>().max(1) as f64;
            let b = |a: i32, z: i32| (0..180usize).filter(|&d| { let l = (d as i32 - 90).abs(); l >= a && l < z }).map(|d| h[d]).sum::<i64>() as f64 / tot;
            [b(0, 30), b(30, 60), b(60, 91)]
        };
        let fresh = bands(0);
        let re = bands(180);
        let exits = t.sh_exits[0].max(1) as f64;
        println!(
            "--- {run} (h {}): captures/step {:.1} (envelope {:.1}, tank releases {:.1}, fresh {:.3}, recaptured {:.3}) | held {} | steady {} ---",
            engine.config.h,
            caps as f64 / avg as f64,
            t.env_capt[0] as f64 / avg as f64,
            t.hold[0][1] as f64 / avg as f64,
            t.feed[0][0] as f64 / caps.max(1) as f64,
            1.0 - t.feed[0][0] as f64 / caps.max(1) as f64,
            engine.held[0],
            rep.steady_ok
        );
        println!(
            "   intake latitude, fresh: 0-30 {:.3} 30-60 {:.3} 60-90 {:.3} (polar/expect {:.3}, equatorial/expect {:.3}) | recaptured: {:.3}/{:.3}/{:.3}",
            fresh[0], fresh[1], fresh[2], fresh[2] / 0.134, fresh[0] / 0.5, re[0], re[1], re[2]
        );
        println!(
            "   contacted-exit bands p/m/e {:.3}/{:.3}/{:.3} | pp collisions/step {:.0} | maze esc-up/esc-down/trap-up/trap-down {:?}",
            t.sh_exit_band[0][0] as f64 / exits,
            t.sh_exit_band[0][1] as f64 / exits,
            t.sh_exit_band[0][2] as f64 / exits,
            t.collisions as f64 / avg as f64,
            &t.maze[0..4]
        );
        // Net radial flow by band at r = 1.5, 2, 3 (x1e3, pol / mid / eq).
        let mut line = String::from("   net radial flow x1e3 (pol/mid/eq):");
        for (k, r) in crate::gpu::dsmc::gpu_types::FLUX_RADII.iter().enumerate() {
            if ![1.5, 2.0, 3.0, 5.0].contains(r) {
                continue;
            }
            line += &format!(" r={r}:");
            for (lo, hi) in [(60, 90), (30, 60), (0, 30)] {
                let (mut o, mut i) = (0i64, 0i64);
                for lb in 0..18usize {
                    let a = (lb as i32 * 10 - 90 + 5).abs();
                    if a >= lo && a < hi {
                        for sp in 0..2 {
                            o += t.flux[((k * 18 + lb) * 2 + sp) * 2];
                            i += t.flux[((k * 18 + lb) * 2 + sp) * 2 + 1];
                        }
                    }
                }
                line += &format!(" {:+.1}", 1e3 * (o - i) as f64 / (o + i).max(1) as f64);
            }
        }
        println!("{line}");
        // Doc §15.14: the circulation in absolute terms, net crossings per
        // capture (polar in, equatorial out) at r = 1.5 and 2.
        let mut circ = String::from("   net flow per capture (pol in / eq out):");
        for (k, r) in crate::gpu::dsmc::gpu_types::FLUX_RADII.iter().enumerate() {
            if ![1.5, 2.0].contains(r) {
                continue;
            }
            let band = |lo: i32, hi: i32| -> i64 {
                let mut net = 0i64;
                for lb in 0..18usize {
                    let a = (lb as i32 * 10 - 90 + 5).abs();
                    if a >= lo && a < hi {
                        for sp in 0..2 {
                            net += t.flux[((k * 18 + lb) * 2 + sp) * 2] - t.flux[((k * 18 + lb) * 2 + sp) * 2 + 1];
                        }
                    }
                }
                net
            };
            circ += &format!(" r={r}: {:.3} / {:.3}", -band(60, 90) as f64 / caps.max(1) as f64, band(0, 30) as f64 / caps.max(1) as f64);
        }
        println!("{circ}");
        // Doc §15.15: north / south split (signed latitude > 60 / < -60):
        // where the body emits, where fresh intake arrives, and the net
        // flow per capture at r = 1.5 (+ = inward).
        let em = &t.maze_lat;
        let em_tot = em.iter().sum::<i64>().max(1) as f64;
        let em_n = (150..180).map(|d| em[d]).sum::<i64>() as f64 / em_tot;
        let em_s = (0..30).map(|d| em[d]).sum::<i64>() as f64 / em_tot;
        let fr = &t.capt_lat[0..180];
        let fr_tot = fr.iter().sum::<i64>().max(1) as f64;
        let fr_n = (150..180).map(|d| fr[d]).sum::<i64>() as f64 / fr_tot;
        let fr_s = (0..30).map(|d| fr[d]).sum::<i64>() as f64 / fr_tot;
        let k15 = crate::gpu::dsmc::gpu_types::FLUX_RADII.iter().position(|&r| r == 1.5).unwrap();
        let net_in = |lbs: std::ops::Range<usize>| -> f64 {
            let mut net = 0i64;
            for lb in lbs {
                for sp in 0..2 {
                    net += t.flux[((k15 * 18 + lb) * 2 + sp) * 2 + 1] - t.flux[((k15 * 18 + lb) * 2 + sp) * 2];
                }
            }
            net as f64 / caps.max(1) as f64
        };
        println!(
            "   north/south: emission {:.3} / {:.3} | fresh intake {:.3} / {:.3} (no field 0.067 each) | net inflow per capture at r=1.5 {:+.3} / {:+.3}",
            em_n, em_s, fr_n, fr_s, net_in(15..18), net_in(0..3)
        );
    }
}

/// Report (doc §15.12): a short stack of up to 4 maze bodies in a row, along
/// the poles (axial) or in each other's equatorial plane (side). The loop
/// test of §15.11b per body: envelope captures with the maze disc against
/// the isotropic-emission control, under polar-only intake, so a closed
/// emission -> steering -> intake loop shows as more captures with the
/// maze. Also what fed each body (never captured / itself / others).
/// Env: CFM_S_CLASSES ("pnpn"), CFM_S_GEOM (axial | side), CFM_S_SEP (2.2),
/// CFM_S_INTAKE (polar0.6), CFM_S_EMIT ("swept,iso"), CFM_S_HOLD (0), plus
/// CFM_C_LAMBDA / NA / DEPTH / WARMUP / AVG / SEED / CAP.
#[test]
#[ignore]
fn report_stack() {
    use super::config::Intake;
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let classes_s = env("CFM_S_CLASSES", "pnpn");
    let classes: Vec<&str> = classes_s.chars().map(|c| if c == 'n' { "neutron" } else { "proton" }).collect();
    let geom = env("CFM_S_GEOM", "axial");
    let sep: f64 = env("CFM_S_SEP", "2.2").parse().unwrap();
    let intake_s = env("CFM_S_INTAKE", "polar0.6");
    let intake = if let Some(p) = intake_s.strip_prefix("env") {
        Intake::Envelope { p: p.parse().unwrap() }
    } else if let Some(a) = intake_s.strip_prefix("polar") {
        Intake::Polar { aperture: a.parse().unwrap() }
    } else {
        Intake::Ball
    };
    let hold: f64 = env("CFM_S_HOLD", "0").parse().unwrap();
    let lambda: f64 = env("CFM_C_LAMBDA", "0.3").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "4000").parse().unwrap();
    let depth: usize = env("CFM_C_DEPTH", "0").parse().unwrap();
    let warmup: u32 = env("CFM_C_WARMUP", "300").parse().unwrap();
    let avg: u32 = env("CFM_C_AVG", "1200").parse().unwrap();
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    let cap: u32 = env("CFM_C_CAP", "9000000").parse().unwrap();
    let e = if geom == "side" { DVec3::X } else { DVec3::Z };
    let n = classes.len();
    let centers: Vec<DVec3> = (0..n).map(|i| e * (i as f64 - (n as f64 - 1.0) / 2.0) * sep).collect();
    println!("=== stack {classes_s} {geom}, sep {sep}, intake {intake:?}, hold {hold}, lambda {lambda}, n_a {n_a}, depth {depth}, warmup {warmup}, avg {avg} ===");
    let ctx = GpuContext::new().expect("GPU");
    let mut rates: Vec<(String, Vec<f64>)> = Vec::new();
    for emit in env("CFM_S_EMIT", "swept,iso").split(',') {
        let maze = if emit == "iso" { MazeMode::Isotropic } else { MazeMode::Swept };
        let mut config = two_body_config(&centers, &classes, n_a, lambda, cap, seed, 0.05);
        config.rough = RoughParams { intake, maze, hold_time: hold, ..config.rough };
        let mut engine = DsmcEngine::new(&ctx, config);
        for (b, cl) in engine.shell_bodies.iter_mut().zip(&classes) {
            b.depth = depth;
            if *cl == "neutron" {
                b.inner_signs = NEUTRON_INNER_SIGNS;
                b.tier1_signs = NEUTRON_INNER_SIGNS;
            }
        }
        let rep = engine.run_steady(&ctx, warmup, avg, 200);
        let t = &rep.totals;
        println!("--- emission {emit}, steady {} ---", rep.steady_ok);
        let mut r = Vec::new();
        for b in 0..n {
            let caps: i64 = t.feed[b].iter().sum();
            let own = t.feed[b][1 + b];
            let others: i64 = (0..n).filter(|&o| o != b).map(|o| t.feed[b][1 + o]).sum();
            let ex = t.sh_exits[b].max(1) as f64;
            println!(
                "   body {b} {:<7} at {:+.2}: envelope captures/step {:.1} (all {:.1}) | fed by never-captured {:.3} itself {:.3} others {:.3} | contacted-exit p/m/e {:.3}/{:.3}/{:.3}",
                classes[b],
                centers[b].dot(e),
                t.env_capt[b] as f64 / avg as f64,
                caps as f64 / avg as f64,
                t.feed[b][0] as f64 / caps.max(1) as f64,
                own as f64 / caps.max(1) as f64,
                others as f64 / caps.max(1) as f64,
                t.sh_exit_band[b][0] as f64 / ex,
                t.sh_exit_band[b][1] as f64 / ex,
                t.sh_exit_band[b][2] as f64 / ex
            );
            r.push(t.env_capt[b] as f64 / avg as f64);
        }
        rates.push((emit.to_string(), r));
    }
    if rates.len() == 2 {
        let line: Vec<String> = (0..n).map(|b| format!("{} {:+.1}%", classes[b], 100.0 * (rates[0].1[b] / rates[1].1[b] - 1.0))).collect();
        println!("   loop gain per body ({} vs {} envelope captures): {}", rates[0].0, rates[1].0, line.join(", "));
    }
}

/// Doc §15.13 gate: in cells crowded past PAR_MIN_CANDIDATES the collide
/// kernel runs NTC candidates in parallel rounds of disjoint pairs. The
/// collision rate must still match n*sigma*<g> (<g> = 4/3 for isotropic
/// unit speeds), and the dense path must actually be taken.
#[test]
fn dense_cells_parallel_collision_rate() {
    let half: f64 = 2.0;
    let per_cell: u32 = 1200;
    let n: u32 = per_cell * 64; // h = 1: 64 unit cells
    let n_density = n as f64 / (2.0 * half).powi(3);
    let sigma = 0.5; // ~14k candidates per cell per step at dt 0.02
    let mut config = shell_config(Boundary::PeriodicBox { half_extent: half }, [0.0; 4], 0.0, twist_rough(false, SpeedMode::PinnedC), sigma, 0.0, 0.5, n + 16, 0xDE5E);
    config.body_model = BodyModel::Ball;
    config.bodies.clear();
    config.shell_signs.clear();
    config.h = 1.0;
    let ctx = GpuContext::new().expect("GPU");
    let mut engine = DsmcEngine::new(&ctx, config);
    let mut rng = Lcg(0xDE5E);
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
    let steps = 20;
    for _ in 0..steps {
        engine.step_collide(&ctx, &mut t);
    }
    let rate = 2.0 * t.collisions as f64 / (n as f64 * steps as f64 * engine.config.dt);
    let expect = n_density * sigma * 4.0 / 3.0;
    println!("dense cells: parallel cell-steps {} of {}, collision rate {rate:.2} vs expected {expect:.2}", t.par_cells, 64 * steps, );
    assert!(t.par_cells as u32 > 60 * steps, "dense path not taken ({})", t.par_cells);
    assert!((rate / expect - 1.0).abs() < 0.03, "dense-cell collision rate distorted");
}

/// Report (doc §15.15): the force between two maze bodies vs separation.
/// Envelope intake p = 1, EjectC. The force on each body is the momentum
/// delivered by its contacts and captures (absorption + emission recoil,
/// pause.html), per unit time, projected on the line from body 0 to body 1.
/// "toward" > 0 means the body is pushed toward the other one. Isolated
/// baselines first (a lone body should feel ~0). Env: CFM_F_PAIR (pn),
/// CFM_F_GEOM (axial | side), CFM_F_SEPS ("2,2.5,3,4"), CFM_C_LAMBDA (0.05),
/// CFM_C_NA (32000), CFM_C_H (0.05), CFM_C_DT (lambda/5), CFM_I_PHI (2/3),
/// CFM_F_MARGIN (domain radius beyond the bodies, 2), CFM_C_WARMUP (1500),
/// CFM_C_AVG (3000), CFM_C_SEED.
#[test]
#[ignore]
fn report_pair_force() {
    use super::config::{Intake, MazeSpeed};
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let pair = env("CFM_F_PAIR", "pn");
    let geom = env("CFM_F_GEOM", "axial");
    let seps: Vec<f64> = env("CFM_F_SEPS", "2,2.5,3,4").split(',').filter(|x| !x.is_empty()).map(|x| x.parse().unwrap()).collect();
    let lambda: f64 = env("CFM_C_LAMBDA", "0.05").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "32000").parse().unwrap();
    let h: f64 = env("CFM_C_H", "0.05").parse().unwrap();
    let dt: f64 = env("CFM_C_DT", &format!("{}", lambda / 5.0)).parse().unwrap();
    let phi: f64 = env("CFM_I_PHI", "0.6666666667").parse().unwrap();
    let margin: f64 = env("CFM_F_MARGIN", "2").parse().unwrap();
    // CFM_F_RD: fixed escape radius for every run (0 = sep/2 + 1 + margin),
    // so the isolated self-force can be subtracted at the same domain.
    let rd_fixed: f64 = env("CFM_F_RD", "0").parse().unwrap();
    let rd = |half: f64| if rd_fixed > 0.0 { rd_fixed } else { half + 1.0 + margin };
    let warmup: u32 = env("CFM_C_WARMUP", "1500").parse().unwrap();
    let avg: u32 = env("CFM_C_AVG", "3000").parse().unwrap();
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    let classes: Vec<&str> = pair.chars().map(|c| if c == 'n' { "neutron" } else { "proton" }).collect();
    let e = if geom == "axial" { DVec3::Z } else { DVec3::X };
    println!("=== pair force: {pair} {geom}, lambda {lambda}, n_a {n_a}, h {h}, dt {dt}, phi {phi:.3}, warmup {warmup}, avg {avg} ===");
    let ctx = GpuContext::new().expect("GPU");
    let run = |centers: &[DVec3], cls: &[&str], r_d: f64| -> (Totals, f64) {
        // Buffer: 1.6x the ambient count, under the 2 GiB binding limit.
        let cap = ((n_a * 4.0 / 3.0 * std::f64::consts::PI * r_d.powi(3) * 1.6) as u32).min(19_000_000);
        let mut config = two_body_config(centers, cls, n_a, lambda, cap, seed, dt);
        config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, maze_speed: MazeSpeed::EjectC, ..config.rough };
        config.h = h;
        config.photon_fraction = phi;
        config.boundary = Boundary::SphereEscape { r_d };
        let mut engine = DsmcEngine::new(&ctx, config);
        for (b, cl) in engine.shell_bodies.iter_mut().zip(cls) {
            b.depth = 0;
            if *cl == "neutron" {
                b.inner_signs = NEUTRON_INNER_SIGNS;
                b.tier1_signs = NEUTRON_INNER_SIGNS;
            }
        }
        // CFM_F_BLOCK > 0: print every body's force (per unit time, along
        // the pair axis e) per block of steps, warm-up included, to see
        // when it settles.
        let block: u32 = std::env::var("CFM_F_BLOCK").ok().and_then(|b| b.parse().ok()).unwrap_or(0);
        let mut bt = Totals::default();
        let mut trace = String::new();
        let mut t = Totals::default();
        for s in 0..warmup + avg {
            let st = engine.step_flux(&ctx);
            if s >= warmup {
                t.add(&st);
            }
            if block > 0 {
                bt.add(&st);
                if (s + 1) % block == 0 {
                    let tb = block as f64 * dt;
                    trace += &format!(" [t {:.0}:", (s + 1) as f64 * dt);
                    for b in 0..centers.len() {
                        trace += &format!(" {:+.0}", DVec3::new(bt.sh_dp[b][0], bt.sh_dp[b][1], bt.sh_dp[b][2]).dot(e) / tb);
                    }
                    trace += "]";
                    bt = Totals::default();
                }
            }
        }
        if block > 0 {
            println!("   force along e by block:{trace}");
        }
        (t, avg as f64 * dt)
    };
    let fmt = |t: &Totals, b: usize, time: f64| -> (DVec3, f64) {
        let f = DVec3::new(t.sh_dp[b][0], t.sh_dp[b][1], t.sh_dp[b][2]) / time;
        let caps = t.feed[b].iter().sum::<i64>() as f64 / time;
        (f, caps)
    };
    let mut uniq = classes.clone();
    uniq.dedup();
    for cl in uniq {
        let (t, time) = run(&[DVec3::ZERO], &[cl], rd(0.0));
        let (f, caps) = fmt(&t, 0, time);
        println!("isolated {cl:<7}: captures/time {caps:.0} | force/time ({:+.1}, {:+.1}, {:+.1}) = ({:+.4}, {:+.4}, {:+.4}) per capture", f.x, f.y, f.z, f.x / caps, f.y / caps, f.z / caps);
    }
    for &sep in &seps {
        let centers = [-e * sep / 2.0, e * sep / 2.0];
        let (t, time) = run(&centers, &classes, rd(sep / 2.0));
        // CFM_F_SOLO=1: also run each body alone at its own pair position in
        // the same domain and report the interaction = pair - solo (removes
        // the self-force and its domain-position dependence).
        let solo: Vec<Option<(DVec3, f64)>> = if env("CFM_F_SOLO", "0") == "1" {
            (0..2).map(|b| { let (ts, tm) = run(&centers[b..b + 1], &classes[b..b + 1], rd(sep / 2.0)); Some(fmt(&ts, 0, tm)) }).collect()
        } else {
            vec![None, None]
        };
        let mut line = format!("sep {sep:.2}:");
        for b in 0..2 {
            let (f, caps) = fmt(&t, b, time);
            let toward = if b == 0 { f.dot(e) } else { -f.dot(e) };
            let fed_nb = t.feed[b][1 + (1 - b)] as f64 / t.feed[b].iter().sum::<i64>().max(1) as f64;
            line += &format!(
                " | body {b} {}: captures/time {caps:.0}, toward other {toward:+.1} ({:+.4} per capture), sideways {:.1}, fed by other {fed_nb:.3}",
                classes[b],
                toward / caps,
                (f - e * f.dot(e)).length()
            );
            if let Some((fs, cs)) = solo[b] {
                let ts = if b == 0 { fs.dot(e) } else { -fs.dot(e) };
                line += &format!(" [solo here: toward {ts:+.1}, captures {cs:.0}; INTERACTION toward other {:+.1} ({:+.4} per solo capture)]", toward - ts, (toward - ts) / cs);
            }
        }
        println!("{line}");
    }
}

/// Report (doc §15.17): a free deuteron. A proton and a neutron start pole
/// to pole and move under the photon traffic alone. Self-force rule (user
/// decision (b), 2026-10-01): each body's own isolated self-force, measured
/// first in the same field and domain, is subtracted, so only the
/// interaction moves the bodies. Per step, body i gets impulse
/// J = dp_i - F_self,i dt and v += J / M - damp v dt, x += v dt (translation
/// only, orientation fixed; axial motion only unless CFM_D_AXIAL=0). Env:
/// CFM_F_PAIR (pn), CFM_F_GEOM (axial | side), CFM_D_SEP0 (3), CFM_D_MASS
/// (1e6, momentum-per-velocity in photon units: the calibration knob),
/// CFM_D_DAMP (0), CFM_D_FIXT (250 time units with the pair held),
/// CFM_D_CALWARMT (250) + CFM_D_CALT (150) for each self-force,
/// CFM_D_STEPS (5000), CFM_D_LOG (100), CFM_D_CONTACT (2: hard contact
/// at envelope touching; 0 = none), CFM_D_LAYOUT / CFM_D_ORIENTS /
/// CFM_D_SEPS (stacks, §15.18),
/// CFM_C_LAMBDA (0.1), CFM_C_NA (8000), CFM_C_H (0.1), CFM_C_DT (lambda/5),
/// CFM_I_PHI (2/3), CFM_F_RD (4.5), CFM_C_SEED.
#[test]
#[ignore]
fn report_deuteron() {
    use super::config::{Intake, MazeSpeed};
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let pair = env("CFM_F_PAIR", "pn");
    let geom = env("CFM_F_GEOM", "axial");
    let sep0: f64 = env("CFM_D_SEP0", "3").parse().unwrap();
    let mass: f64 = env("CFM_D_MASS", "1e6").parse().unwrap();
    let damp: f64 = env("CFM_D_DAMP", "0").parse().unwrap();
    // Warm-ups and the calibration are in time units (doc §15.17: the
    // neutron's self-force relaxes over ~150-200 time units at lambda 0.1
    // and 0.3, so anything shorter measures a transient).
    let lambda_pre: f64 = env("CFM_C_LAMBDA", "0.1").parse().unwrap();
    let dt_pre: f64 = env("CFM_C_DT", &format!("{}", lambda_pre / 5.0)).parse().unwrap();
    let steps_of = |k: &str, d: &str| -> u32 { (env(k, d).parse::<f64>().unwrap() / dt_pre).round() as u32 };
    let fix = steps_of("CFM_D_FIXT", "250");
    let cal = steps_of("CFM_D_CALT", "150");
    let cal_warm = steps_of("CFM_D_CALWARMT", "250");
    let steps: u32 = env("CFM_D_STEPS", "5000").parse().unwrap();
    let log: u32 = env("CFM_D_LOG", "100").parse().unwrap();
    let axial_only = env("CFM_D_AXIAL", "1") == "1";
    // CFM_D_FIXCOM=1 (default): only the separation evolves. The pair's
    // centre stays at the domain centre (the common-mode push of the open,
    // non-Newtonian forces is reported, not integrated), which keeps both
    // bodies where the centre-calibrated self-force applies.
    let fix_com = env("CFM_D_FIXCOM", "1") == "1";
    // CFM_D_CONTACT (doc §15.18): bodies are hard spheres of this centre
    // distance (2 = envelopes touching; 0 = off, they pass through each
    // other). Overlap is resolved by an elastic, equal-mass contact.
    let contact: f64 = env("CFM_D_CONTACT", "2").parse().unwrap();
    let lambda: f64 = env("CFM_C_LAMBDA", "0.1").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "8000").parse().unwrap();
    let h: f64 = env("CFM_C_H", "0.1").parse().unwrap();
    let dt: f64 = env("CFM_C_DT", &format!("{}", lambda / 5.0)).parse().unwrap();
    let phi: f64 = env("CFM_I_PHI", "0.6666666667").parse().unwrap();
    let r_d: f64 = env("CFM_F_RD", "4.5").parse().unwrap();
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    // Doc §15.23: 'e' in CFM_F_PAIR is an electron (lepton ball, radius
    // CFM_D_RE, mass CFM_D_ME, no intake; spin axis from CFM_D_AXES).
    let classes: Vec<&str> = pair.chars().map(|c| if c == 'n' { "neutron" } else if c == 'e' { "electron" } else { "proton" }).collect();
    let r_e: f64 = env("CFM_D_RE", "0.2").parse().unwrap();
    let m_e: f64 = env("CFM_D_ME", "3800").parse().unwrap();
    let e = if geom == "axial" { DVec3::Z } else { DVec3::X };
    println!("=== deuteron: {pair} {geom}, sep0 {sep0}, mass {mass:e}, damp {damp}, fixed centre {fix_com}, lambda {lambda}, n_a {n_a}, h {h}, dt {dt}, phi {phi:.3}, r_d {r_d} ===");
    let ctx = GpuContext::new().expect("GPU");
    // CFM_D_INTAKE: env<p> (default env1) | polar<aperture> (doc §15.23).
    let intake_s = env("CFM_D_INTAKE", "env1");
    let d_intake = if let Some(a) = intake_s.strip_prefix("polar") { Intake::Polar { aperture: a.parse().unwrap() } } else { Intake::Envelope { p: intake_s.trim_start_matches("env").parse().unwrap() } };
    // CFM_D_HOLE (doc §15.24): through-charge hole radius, 0 = off.
    let d_hole: f64 = env("CFM_D_HOLE", "0").parse().unwrap();
    // CFM_D_SPEED: surface (default, doc §15.29) | rest (§15.27, recaptured by the
    // inner ball) | ejectc (§15.14).
    let d_speed = match env("CFM_D_SPEED", "surface").as_str() { "ejectc" => MazeSpeed::EjectC, "rest" => MazeSpeed::EjectRest, _ => MazeSpeed::EjectSurface };
    println!("intake {d_intake:?}, through hole {d_hole}, maze speed {d_speed:?}");
    let make = |centers: &[DVec3], cls: &[&str], ors: &[GpuOrient]| -> DsmcEngine {
        let cap = ((n_a * 4.0 / 3.0 * std::f64::consts::PI * r_d.powi(3) * 1.6) as u32).min(19_000_000);
        let mut config = two_body_config(centers, cls, n_a, lambda, cap, seed, dt);
        for ((bc, o), cl) in config.bodies.iter_mut().zip(ors).zip(cls) {
            bc.orient = *o;
            if *cl == "electron" {
                bc.class = GpuClass::Electron;
            }
        }
        // CFM_D_POLAR=normal (doc §15.30): trapped releases leave along the
        // surface normal over the polar opening.
        config.rough = RoughParams { intake: d_intake, maze_speed: d_speed, through_hole: d_hole, polar_release_normal: env("CFM_D_POLAR", "axial") == "normal", ..config.rough };
        config.h = h;
        config.photon_fraction = phi;
        config.boundary = Boundary::SphereEscape { r_d };
        let mut engine = DsmcEngine::new(&ctx, config);
        for (b, cl) in engine.shell_bodies.iter_mut().zip(cls) {
            b.depth = 0;
            if *cl == "neutron" {
                b.inner_signs = NEUTRON_INNER_SIGNS;
                b.tier1_signs = NEUTRON_INNER_SIGNS;
            }
        }
        for (b, cl) in cls.iter().enumerate() {
            if *cl == "electron" {
                apply_level_speeds(&mut engine);
        engine.make_lepton(&ctx, b, r_e);
            }
        }
        apply_level_speeds(&mut engine);
        engine
    };
    // 1. Self-force per class, body alone at the domain centre.
    let mut f_self = std::collections::HashMap::new();
    // Doc §15.20: the isolated self-torque, in the body frame (calibrated
    // pole up, frame = identity), subtracted like the self-force.
    let mut t_self: std::collections::HashMap<&str, DVec3> = std::collections::HashMap::new();
    let mut uniq = classes.clone();
    uniq.sort();
    uniq.dedup();
    // CFM_D_FSELF="neutron:x,y,z;proton:x,y,z" reuses settled self-forces
    // measured earlier in the same field and domain (skips calibration).
    for item in env("CFM_D_FSELF", "").split(';').filter(|q| !q.is_empty()) {
        let (cl, f) = item.split_once(':').unwrap();
        let c: Vec<f64> = f.split(',').map(|q| q.parse().unwrap()).collect();
        let key = *uniq.iter().find(|u| **u == cl).expect("CFM_D_FSELF class not in the stack");
        println!("self-force {key:<7}: ({:+.1}, {:+.1}, {:+.1}) per unit time (given)", c[0], c[1], c[2]);
        f_self.insert(key, DVec3::new(c[0], c[1], c[2]));
        if c.len() >= 6 {
            println!("self-torque {key:<7}: ({:+.1}, {:+.1}, {:+.1}) per unit time (given)", c[3], c[4], c[5]);
            t_self.insert(key, DVec3::new(c[3], c[4], c[5]));
        }
    }
    let inertia: f64 = env("CFM_D_INERTIA", "0").parse().unwrap();
    for cl in uniq.clone() {
        let want_torque = inertia > 0.0 || env("CFM_D_FREEZE", "0") == "1";
        if f_self.contains_key(cl) && (!want_torque || t_self.contains_key(cl)) {
            continue;
        }
        let mut engine = make(&[DVec3::ZERO], &[cl], &[GpuOrient::PoleOn]);
        for _ in 0..cal_warm {
            engine.step_flux(&ctx);
        }
        let mut t = Totals::default();
        for _ in 0..cal {
            t.add(&engine.step_flux(&ctx));
        }
        let f = DVec3::new(t.sh_dp[0][0], t.sh_dp[0][1], t.sh_dp[0][2]) / (cal as f64 * dt);
        let tq = DVec3::new(t.sh_torque[0][0], t.sh_torque[0][1], t.sh_torque[0][2]) / (cal as f64 * dt);
        println!("self-force {cl:<7}: ({:+.1}, {:+.1}, {:+.1}) per unit time", f.x, f.y, f.z);
        println!("self-torque {cl:<7}: ({:+.1}, {:+.1}, {:+.1}) per unit time", tq.x, tq.y, tq.z);
        f_self.insert(cl, f);
        t_self.insert(cl, tq);
    }
    // 2. The stack, held while the field settles, then free. Bodies sit on
    // the axis e in CFM_F_PAIR order with gaps CFM_D_SEPS (default all
    // CFM_D_SEP0), centred on the domain.
    // CFM_D_ORIENTS: one letter per body, u = pole up (+z, default),
    // d = pole down (-z). CFM_D_LAYOUT=alpha (doc §15.18, deut.pdf's He
    // sandwich): CFM_F_PAIR must be "pnnp"; protons on the axis at z = +-a
    // meeting pole to pole (top one upside down), neutrons side by side at
    // x = -+b, antiparallel (each sits under/over its own proton, as in the
    // free deuteron). CFM_D_A, CFM_D_B (2.16 each: p-n distance 3.05).
    let n = classes.len();
    assert!(n >= 2 && n <= 4, "2-4 bodies");
    let layout = env("CFM_D_LAYOUT", "chain");
    let parse_ors = |d: &str| -> Vec<GpuOrient> { env("CFM_D_ORIENTS", d).chars().map(|c| if c == 'd' { GpuOrient::PoleDown } else { GpuOrient::PoleOn }).collect() };
    // CFM_D_LAYOUT=custom (doc §15.21): CFM_D_POS "x,y,z;x,y,z;..." and
    // CFM_D_AXES "x,y,z;..." (spin axes, default +z), any orientation.
    let parse_v3s = |k: &str| -> Vec<DVec3> {
        env(k, "").split(';').filter(|q| !q.is_empty()).map(|q| { let c: Vec<f64> = q.split(',').map(|v| v.trim().parse().unwrap()).collect(); DVec3::new(c[0], c[1], c[2]) }).collect()
    };
    let custom_axes: Option<Vec<DVec3>> = if layout == "custom" {
        let ax = parse_v3s("CFM_D_AXES");
        Some(if ax.is_empty() { vec![DVec3::Z; n] } else { ax.iter().map(|a| a.normalize()).collect() })
    } else {
        None
    };
    let (mut x, ors): (Vec<DVec3>, Vec<GpuOrient>) = if layout == "custom" {
        let pos = parse_v3s("CFM_D_POS");
        assert_eq!(pos.len(), n, "CFM_D_POS needs one position per body");
        (pos, vec![GpuOrient::PoleOn; n])
    } else if layout == "sandwich" {
        // Mathis's He sandwich as drawn (user, 2026-10-01): the neutrons
        // stand side by side between the proton discs, antiparallel; both
        // protons turn the same way about the nuclear axis. Compact: every
        // neighbour touching (n-n = p-n = 2 at CFM_D_B 1, CFM_D_A sqrt 3).
        // Order pnnp; CFM_D_ORIENTS default "uudu".
        assert_eq!(pair, "pnnp", "sandwich layout needs CFM_F_PAIR=pnnp");
        let a: f64 = env("CFM_D_A", "1.7320508").parse().unwrap();
        let b: f64 = env("CFM_D_B", "1").parse().unwrap();
        (vec![DVec3::new(0.0, 0.0, -a), DVec3::new(b, 0.0, 0.0), DVec3::new(-b, 0.0, 0.0), DVec3::new(0.0, 0.0, a)], parse_ors("uudu"))
    } else if layout == "alpha" {
        assert_eq!(pair, "pnnp", "alpha layout needs CFM_F_PAIR=pnnp");
        let a: f64 = env("CFM_D_A", "2.16").parse().unwrap();
        let b: f64 = env("CFM_D_B", "2.16").parse().unwrap();
        (
            vec![DVec3::new(0.0, 0.0, -a), DVec3::new(b, 0.0, 0.0), DVec3::new(-b, 0.0, 0.0), DVec3::new(0.0, 0.0, a)],
            vec![GpuOrient::PoleOn, GpuOrient::PoleOn, GpuOrient::PoleDown, GpuOrient::PoleDown],
        )
    } else {
        let gaps: Vec<f64> = env("CFM_D_SEPS", "").split(',').filter(|g| !g.is_empty()).map(|g| g.parse().unwrap()).collect();
        let gaps: Vec<f64> = if gaps.is_empty() { vec![sep0; n - 1] } else { gaps };
        assert_eq!(gaps.len(), n - 1);
        let mut u = vec![0.0f64; n];
        for k in 1..n {
            u[k] = u[k - 1] + gaps[k - 1];
        }
        let mid = (u[0] + u[n - 1]) / 2.0;
        let os: Vec<GpuOrient> = parse_ors(&"u".repeat(n));
        (u.iter().map(|&q| e * (q - mid)).collect(), os)
    };
    let free3d = layout != "chain";
    let axial_only = axial_only && !free3d;
    // Each body's self-force is its class's (calibrated pole up), turned
    // by the same rotation that orients the body.
    // Doc §15.20 rotation: body frames q (spin axis = q Z). With
    // CFM_D_INERTIA = I > 0 the axis turns under the measured torque minus
    // the frame-turned self-torque. Angular momentum L = L_s a + I w_perp:
    // CFM_D_LSPIN = L_s (0 = plain rigid turning; > 0 a gyroscope that
    // precesses). The axial part is dropped (inner spin rates are driven).
    // CFM_D_RDAMP damps w_perp.
    let l_spin: f64 = env("CFM_D_LSPIN", "0").parse().unwrap();
    let rdamp: f64 = env("CFM_D_RDAMP", "0").parse().unwrap();
    let mut q: Vec<glam::DQuat> = (0..n)
        .map(|b| glam::DQuat::from_rotation_arc(DVec3::Z, custom_axes.as_ref().map(|a| a[b]).unwrap_or(ors[b].spin_axis())))
        .collect();
    // CFM_D_FREEZE=1 (doc §15.21): hold every body (no translation or
    // rotation) and report the mean interaction force and torque on each,
    // the stability screen. With hard contact, a layout whose interaction
    // forces all point into the cluster (compression) holds; any body pulled
    // outward means it comes apart.
    let freeze = env("CFM_D_FREEZE", "0") == "1";
    let mut l_ang: Vec<DVec3> = (0..n).map(|b| q[b] * DVec3::Z * l_spin).collect();
    println!("rotation: inertia {inertia:e}, L_spin {l_spin:e}, damp {rdamp}");
    println!("layout {layout}: orients {ors:?}, start {:?}", x.iter().map(|p| (p.x, p.y, p.z)).collect::<Vec<_>>());
    let mut v = vec![DVec3::ZERO; n];
    let is_e: Vec<bool> = classes.iter().map(|c| *c == "electron").collect();
    let m_b: Vec<f64> = is_e.iter().map(|&e| if e { m_e } else { mass }).collect();
    let rad: Vec<f64> = is_e.iter().map(|&e| if e { r_e } else { contact / 2.0 }).collect();
    let m_tot: f64 = m_b.iter().sum();
    let mut engine = make(&x, &classes, &ors);
    engine.set_body_poses(&ctx, &x, Some(&q));
    for _ in 0..fix {
        engine.step_flux(&ctx);
    }
    let mut win = vec![DVec3::ZERO; n];
    let mut twin = vec![DVec3::ZERO; n];
    let mut fsum = vec![DVec3::ZERO; n];
    let mut tsum = vec![DVec3::ZERO; n];
    let mut steps_done = 0u32;
    let mut caps = vec![0i64; n];
    let mut n_contacts = 0u64;
    let pairs_ij: Vec<(usize, usize)> = (0..n).flat_map(|i| (i + 1..n).map(move |k| (i, k))).collect();
    let n_g = if free3d { pairs_ij.len() } else { n - 1 };
    let mut gap_hist: Vec<Vec<f64>> = vec![Vec::new(); n_g];
    for s in 1..=steps {
        let t = engine.step_flux(&ctx);
        steps_done = s;
        for b in 0..n {
            let dp = DVec3::new(t.sh_dp[b][0], t.sh_dp[b][1], t.sh_dp[b][2]);
            let mut j = dp - (q[b] * f_self[classes[b]]) * dt;
            if axial_only {
                j = e * j.dot(e);
            }
            win[b] += j;
            caps[b] += t.feed[b].iter().sum::<i64>();
            if !freeze {
                let vb = v[b];
                v[b] = vb + j / m_b[b] - vb * damp * dt;
                x[b] += v[b] * dt;
            }
            fsum[b] += j;
        }
        if contact > 0.0 {
            let mut contacts_now = 0;
            for i in 0..n {
                for k in i + 1..n {
                    let d = x[k] - x[i];
                    let r = d.length();
                    let c_ik = rad[i] + rad[k];
                    if r < c_ik && r > 1e-9 {
                        // Elastic, mass-weighted (equal masses: the old
                        // half-and-half push and normal-velocity swap).
                        let nrm = d / r;
                        let (wi, wk) = (m_b[k] / (m_b[i] + m_b[k]), m_b[i] / (m_b[i] + m_b[k]));
                        x[i] -= nrm * (c_ik - r) * wi;
                        x[k] += nrm * (c_ik - r) * wk;
                        let vn = (v[k] - v[i]).dot(nrm);
                        if vn < 0.0 {
                            v[i] += nrm * vn * 2.0 * wi;
                            v[k] -= nrm * vn * 2.0 * wk;
                        }
                        contacts_now += 1;
                    }
                }
            }
            n_contacts += contacts_now;
        }
        if fix_com {
            let com = (0..n).map(|b| x[b] * m_b[b]).sum::<DVec3>() / m_tot;
            let vcom = (0..n).map(|b| v[b] * m_b[b]).sum::<DVec3>() / m_tot;
            for b in 0..n {
                x[b] -= com;
                v[b] -= vcom;
            }
        }
        if freeze && t_self.len() == f_self.len() {
            for b in 0..n {
                let a = q[b] * DVec3::Z;
                let tq = DVec3::new(t.sh_torque[b][0], t.sh_torque[b][1], t.sh_torque[b][2]);
                let jr = tq - (q[b] * t_self[classes[b]]) * dt;
                tsum[b] += jr - a * jr.dot(a);
            }
        }
        if inertia > 0.0 && !freeze {
            for b in 0..n {
                if is_e[b] {
                    continue;
                }
                let a = q[b] * DVec3::Z;
                let tq = DVec3::new(t.sh_torque[b][0], t.sh_torque[b][1], t.sh_torque[b][2]);
                let jr = tq - (q[b] * t_self[classes[b]]) * dt;
                let jr = jr - a * jr.dot(a);
                twin[b] += jr;
                l_ang[b] += jr;
                // Damping acts on the tumble part of L only.
                let lperp = l_ang[b] - a * l_ang[b].dot(a);
                l_ang[b] -= lperp * rdamp * dt;
                let w = (l_ang[b] - a * l_ang[b].dot(a)) / inertia;
                let ang = w.length() * dt;
                if ang > 0.0 {
                    q[b] = (glam::DQuat::from_axis_angle(w / w.length(), ang) * q[b]).normalize();
                }
                // L itself is unchanged by the turn (only torque changes it);
                // the driven inner spins hold its component along the new
                // axis at L_s. Keeping L (not rebuilding it from w) is what
                // couples spin and tilt: a gyroscope precesses.
                let a2 = q[b] * DVec3::Z;
                l_ang[b] = l_ang[b] - a2 * l_ang[b].dot(a2) + a2 * l_spin;
            }
        }
        engine.set_body_poses(&ctx, &x, Some(&q));
        let g: Vec<f64> = if free3d {
            pairs_ij.iter().map(|&(i, k)| (x[k] - x[i]).length()).collect()
        } else {
            (0..n - 1).map(|k| (x[k + 1] - x[k]).dot(e)).collect()
        };
        for k in 0..g.len() {
            gap_hist[k].push(g[k]);
        }
        if s % log == 0 {
            let tw = log as f64 * dt;
            let mut line = format!("t {:7.2}: gaps", s as f64 * dt);
            for k in 0..g.len() {
                line += &format!(" {:.4}", g[k]);
            }
            if free3d {
                line += " | interaction (x,z):";
                for b in 0..n {
                    line += &format!(" {} ({:+.0},{:+.0})", &classes[b][..1], win[b].x / tw, win[b].z / tw);
                }
                line += " | max lateral y";
                line += &format!(" {:.3}", x.iter().map(|p| p.y.abs()).fold(0.0, f64::max));
            } else {
                line += " | gap rates";
                for k in 0..n - 1 {
                    line += &format!(" {:+.5}", (v[k + 1] - v[k]).dot(e));
                }
                line += " | interaction along +e:";
                for b in 0..n {
                    line += &format!(" {} {:+.0}", &classes[b][..1], win[b].dot(e) / tw);
                }
            }
            if inertia > 0.0 {
                line += " | axes";
                for b in 0..n {
                    let a = q[b] * DVec3::Z;
                    line += &format!(" ({:+.2},{:+.2},{:+.2})", a.x, a.y, a.z);
                }
                line += " | net torque";
                for b in 0..n {
                    line += &format!(" {:.0}", twin[b].length() / tw);
                }
                twin = vec![DVec3::ZERO; n];
            }
            // Each electron: distance to and latitude on its nearest baryon
            // (latitude in that baryon's frame, + = toward its pole).
            for b in (0..n).filter(|&b| is_e[b]) {
                if let Some(k) = (0..n).filter(|&k| !is_e[k]).min_by(|&i, &j| (x[i] - x[b]).length().total_cmp(&(x[j] - x[b]).length())) {
                    let d = x[b] - x[k];
                    let lat = (d.normalize().dot(q[k] * DVec3::Z)).clamp(-1.0, 1.0).asin().to_degrees();
                    line += &format!(" | e{b} on {}{k}: r {:.2} lat {lat:+.0}", &classes[k][..1], d.length());
                }
            }
            line += &format!(" | contact steps {n_contacts}");
            line += " | captures/time";
            for b in 0..n {
                line += &format!(" {:.0}", caps[b] as f64 / tw);
            }
            println!("{line}");
            win = vec![DVec3::ZERO; n];
            caps = vec![0; n];
        }
        if (contact == 0.0 && g.iter().any(|&q| q < 0.5)) || x.iter().any(|p| p.length() > r_d - 1.2) {
            println!("stopped at t {:.2}: gaps {:?}, centres {:?}", s as f64 * dt, g, x.iter().map(|p| p.dot(e)).collect::<Vec<_>>());
            break;
        }
    }
    for k in 0..n_g {
        let (i, j) = if free3d { pairs_ij[k] } else { (k, k + 1) };
        let h = &gap_hist[k];
        let (mn, mx) = h.iter().fold((f64::MAX, f64::MIN), |(a, b), &q| (a.min(q), b.max(q)));
        let tail = &h[h.len() * 3 / 4..];
        println!("distance {i}-{j} ({}-{}): range {mn:.3} .. {mx:.3}; last quarter mean {:.3}", classes[i], classes[j], tail.iter().sum::<f64>() / tail.len().max(1) as f64);
    }
    println!("drift of the stack centre {:+.3}", (x.iter().sum::<DVec3>() / n as f64).dot(e));
    // Doc §15.21 screen summary: mean interaction force per body over the
    // run, its component toward the cluster centre (+ = compressive), and
    // toward each touching neighbour; mean tilting torque.
    let time = steps_done.max(1) as f64 * dt;
    let cen = x.iter().sum::<DVec3>() / n as f64;
    for b in 0..n {
        let f = fsum[b] / time;
        let inward = (cen - x[b]).try_normalize().map(|u| f.dot(u)).unwrap_or(0.0);
        let mut nb = String::new();
        for k in 0..n {
            if k != b && (x[k] - x[b]).length() < 3.6 {
                nb += &format!(" {}{}:{:+.0}", &classes[k][..1], k, f.dot((x[k] - x[b]).normalize()));
            }
        }
        println!(
            "screen body {b} {:<7} at ({:+.2},{:+.2},{:+.2}) axis ({:+.2},{:+.2},{:+.2}): force ({:+.0},{:+.0},{:+.0}) |{:.0}| inward {inward:+.0} | toward neighbours{nb} | tilt torque {:.0}",
            classes[b],
            x[b].x, x[b].y, x[b].z,
            (q[b] * DVec3::Z).x, (q[b] * DVec3::Z).y, (q[b] * DVec3::Z).z,
            f.x, f.y, f.z, f.length(),
            (tsum[b] / time).length()
        );
    }
}

/// Report (doc §15.22): the hydrogen atom. A proton (maze, envelope intake,
/// held at the origin, pole +z) and a free electron: a lepton ball of radius
/// CFM_H_RE that photons bounce off and that takes nothing in, mass
/// CFM_H_ME, moved only by its contact impulses. It cannot enter the
/// proton's envelope (hard contact at 1 + r_e). diatom.pdf: the electron
/// should be driven off near the equator by emitted charge and drawn to a
/// pole by incoming charge, then circle that pole. Env: CFM_H_START
/// ("x,y,z", 2.5,0,0), CFM_H_RE (0.2), CFM_H_ME (3800 = 7e6/1836),
/// CFM_H_SPIN (+1/-1, electron spin sense about +z), CFM_C_LAMBDA (0.3),
/// CFM_C_NA (8000), CFM_C_H (0.1), CFM_C_DT (lambda/5), CFM_I_PHI (2/3),
/// CFM_F_RD (6), CFM_H_FIXT (250 time units held), CFM_H_STEPS (5000),
/// CFM_H_LOG (250), CFM_C_SEED.
#[test]
#[ignore]
fn report_hydrogen() {
    use super::config::{Intake, MazeSpeed};
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let v3 = |q: &str| -> DVec3 { let c: Vec<f64> = q.split(',').map(|v| v.trim().parse().unwrap()).collect(); DVec3::new(c[0], c[1], c[2]) };
    let start = v3(&env("CFM_H_START", "2.5,0,0"));
    let r_e: f64 = env("CFM_H_RE", "0.2").parse().unwrap();
    let m_e: f64 = env("CFM_H_ME", "3800").parse().unwrap();
    let spin: f64 = env("CFM_H_SPIN", "1").parse().unwrap();
    let lambda: f64 = env("CFM_C_LAMBDA", "0.3").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "8000").parse().unwrap();
    let h: f64 = env("CFM_C_H", "0.1").parse().unwrap();
    let dt: f64 = env("CFM_C_DT", &format!("{}", lambda / 5.0)).parse().unwrap();
    let phi: f64 = env("CFM_I_PHI", "0.6666666667").parse().unwrap();
    let r_d: f64 = env("CFM_F_RD", "6").parse().unwrap();
    let fix = (env("CFM_H_FIXT", "250").parse::<f64>().unwrap() / dt).round() as u32;
    let steps: u32 = env("CFM_H_STEPS", "5000").parse().unwrap();
    let log: u32 = env("CFM_H_LOG", "250").parse().unwrap();
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    println!("=== hydrogen: electron r {r_e}, m {m_e}, spin {spin:+}, start {start:?}, lambda {lambda}, n_a {n_a}, h {h}, dt {dt}, phi {phi:.3}, r_d {r_d} ===");
    let ctx = GpuContext::new().expect("GPU");
    let cap = ((n_a * 4.0 / 3.0 * std::f64::consts::PI * r_d.powi(3) * 1.6) as u32).min(19_000_000);
    let mut config = two_body_config(&[DVec3::ZERO, start], &["proton", "proton"], n_a, lambda, cap, seed, dt);
    config.bodies[1].class = GpuClass::Electron;
    // CFM_H_INTAKE: env<p> (default env1) | polar<aperture> (doc §15.23:
    // a polar drain, photons taken in only within `aperture` of the axis).
    let intake_s = env("CFM_H_INTAKE", "env1");
    let intake = if let Some(a) = intake_s.strip_prefix("polar") { Intake::Polar { aperture: a.parse().unwrap() } } else { Intake::Envelope { p: intake_s.trim_start_matches("env").parse().unwrap() } };
    let hole: f64 = env("CFM_H_HOLE", "0").parse().unwrap();
    println!("intake {intake:?}, through hole {hole}");
    // CFM_H_SPEED: surface (default, doc §15.29) | rest (§15.27, recaptured by the
    // inner ball) | ejectc (§15.14).
    let mspeed = match env("CFM_H_SPEED", "surface").as_str() { "ejectc" => MazeSpeed::EjectC, "rest" => MazeSpeed::EjectRest, _ => MazeSpeed::EjectSurface };
    println!("maze speed {mspeed:?}");
    config.rough = RoughParams { intake, maze_speed: mspeed, through_hole: hole, ..config.rough };
    config.h = h;
    config.photon_fraction = phi;
    config.boundary = Boundary::SphereEscape { r_d };
    let mut engine = DsmcEngine::new(&ctx, config);
    engine.shell_bodies[0].depth = 0;
    let q_e = if spin < 0.0 { DQuat::from_rotation_arc(DVec3::Z, -DVec3::Z) } else { DQuat::IDENTITY };
    engine.make_lepton(&ctx, 1, r_e);
    apply_level_speeds(&mut engine);
    engine.set_body_poses(&ctx, &[DVec3::ZERO, start], Some(&[DQuat::IDENTITY, q_e]));
    for _ in 0..fix {
        engine.step_flux(&ctx);
    }
    let mut x = start;
    let mut v = DVec3::ZERO;
    let mut win = DVec3::ZERO;
    let mut hits = 0i64;
    let mut bands = [0u64; 6]; // |lat| 0-30 / 30-60 / 60-90, then south / north caps (>60), then contact steps
    let mut r_hist = Vec::new();
    let r_min = 1.0 + r_e;
    for s in 1..=steps {
        let t = engine.step_flux(&ctx);
        let dp = DVec3::new(t.sh_dp[1][0], t.sh_dp[1][1], t.sh_dp[1][2]);
        hits += t.sh_contacts[1];
        win += dp;
        v += dp / m_e;
        x += v * dt;
        // Hard contact with the proton's envelope.
        let r = x.length();
        if r < r_min {
            let n = x / r.max(1e-9);
            x = n * r_min;
            let vn = v.dot(n);
            if vn < 0.0 {
                v -= n * vn * 2.0;
            }
            bands[5] += 1;
        }
        if x.length() > r_d - 1.0 {
            println!("electron escaped the domain at t {:.1}: {x:?}", s as f64 * dt);
            break;
        }
        engine.set_body_poses(&ctx, &[DVec3::ZERO, x], Some(&[DQuat::IDENTITY, q_e]));
        let r = x.length();
        let lat = (x.z / r).clamp(-1.0, 1.0).asin().to_degrees();
        bands[if lat.abs() < 30.0 { 0 } else if lat.abs() < 60.0 { 1 } else { 2 }] += 1;
        if lat < -60.0 {
            bands[3] += 1;
        } else if lat > 60.0 {
            bands[4] += 1;
        }
        r_hist.push(r);
        if s % log == 0 {
            let tw = log as f64 * dt;
            let f = win / tw;
            let rh = x / r;
            let f_r = f.dot(rh);
            // Poleward component: along -d(lat) direction toward the nearer pole.
            let e_lat = (DVec3::Z - rh * rh.z).try_normalize().unwrap_or(DVec3::X);
            let f_pole = f.dot(e_lat) * lat.signum();
            // Swirl about the proton axis (+ = the proton's own turning sense).
            let e_az = DVec3::Z.cross(rh).try_normalize().unwrap_or(DVec3::Y);
            println!(
                "t {:7.1}: r {r:.3} lat {lat:+6.1} speed {:.4} (radial {:+.4} swirl {:+.4}) | force/time radial {f_r:+.1} poleward {f_pole:+.1} swirl {:+.1} | photon hits/time {:.0} | proton contacts {}",
                s as f64 * dt,
                v.length(),
                v.dot(rh),
                v.dot(e_az),
                f.dot(e_az),
                hits as f64 / tw,
                bands[5]
            );
            win = DVec3::ZERO;
            hits = 0;
        }
    }
    let n = r_hist.len().max(1) as f64;
    println!(
        "time share |lat| 0-30 {:.3} 30-60 {:.3} 60-90 {:.3} (uniform on a sphere 0.500/0.366/0.134) | south cap {:.3} north cap {:.3} | mean r {:.3}, last half {:.3}",
        bands[0] as f64 / n,
        bands[1] as f64 / n,
        bands[2] as f64 / n,
        bands[3] as f64 / n,
        bands[4] as f64 / n,
        r_hist.iter().sum::<f64>() / n,
        r_hist[r_hist.len() / 2..].iter().sum::<f64>() / (r_hist.len() - r_hist.len() / 2).max(1) as f64
    );
}


/// Doc §15.28: CFM_LEVEL_SPEEDS = "s9,s10,s11,s12", the edge speed of each
/// baryon level as a fraction of c (default all 1). Applied to every
/// non-lepton shell body.
fn apply_level_speeds(engine: &mut DsmcEngine) {
    let v: Vec<f64> = std::env::var("CFM_LEVEL_SPEEDS").unwrap_or("1,1,1,1".into()).split(',').map(|q| q.trim().parse().unwrap()).collect();
    assert_eq!(v.len(), 4, "CFM_LEVEL_SPEEDS needs 4 values (levels 9-12)");
    for sb in engine.shell_bodies.iter_mut().filter(|b| !b.lepton) {
        sb.speeds = [v[0], v[1], v[2], v[3]];
    }
    if v.iter().any(|&q| q != 1.0) {
        println!("level speeds (9-12, fraction of c): {v:?}");
    }
}

/// Report (doc §15.26): the meridional flow map. Bodies are held on the z
/// axis through the domain centre; after the warm-up, photon snapshots are
/// binned by (R, z), R = distance from the axis (an average over azimuth,
/// so a pair must sit on the axis). Per cell: count, mean velocity (v_R,
/// v_z, swirl v_phi = + along the right-handed turn about +z) and counts by
/// emitter tag (0 = never emitted by a body, b + 1 = last emitted by body b).
/// Written as JSON for the flow-map page. The user's question: does the
/// disc outflow come back to the poles as a cone vortex, and how much of a
/// partner stream reaches the core?
/// Env: CFM_M_LAYOUT (p | n | pp_par | pp_nn | pp_ss), CFM_M_SEP (3, centre
/// distance of a pair), CFM_C_LAMBDA (0.1), CFM_C_NA (8000), CFM_C_H (0.1),
/// CFM_F_RD (6), CFM_M_HOLE (0), CFM_M_FIXT (250 time units), CFM_M_SNAPS
/// (80), CFM_M_EVERY (10 steps), CFM_M_OUT (JSON path), CFM_C_SEED.
#[test]
#[ignore]
fn report_flowmap() {
    use super::config::{Intake, MazeSpeed};
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let layout = env("CFM_M_LAYOUT", "p");
    let sep: f64 = env("CFM_M_SEP", "3").parse().unwrap();
    let lambda: f64 = env("CFM_C_LAMBDA", "0.1").parse().unwrap();
    let n_a: f64 = env("CFM_C_NA", "8000").parse().unwrap();
    let h: f64 = env("CFM_C_H", "0.1").parse().unwrap();
    let dt: f64 = env("CFM_C_DT", &format!("{}", lambda / 5.0)).parse().unwrap();
    let phi: f64 = env("CFM_I_PHI", "0.6666666667").parse().unwrap();
    let r_d: f64 = env("CFM_F_RD", "6").parse().unwrap();
    let hole: f64 = env("CFM_M_HOLE", "0").parse().unwrap();
    let fix = (env("CFM_M_FIXT", "250").parse::<f64>().unwrap() / dt).round() as u32;
    let snaps: u32 = env("CFM_M_SNAPS", "80").parse().unwrap();
    let every: u32 = env("CFM_M_EVERY", "10").parse().unwrap();
    let out = env("CFM_M_OUT", "flowmap.json");
    let seed: u64 = env("CFM_C_SEED", "49344").parse().unwrap();
    let up = DQuat::IDENTITY;
    let down = DQuat::from_rotation_arc(DVec3::Z, -DVec3::Z);
    let s2 = sep / 2.0;
    let (centers, classes, orients): (Vec<DVec3>, Vec<&str>, Vec<DQuat>) = match layout.as_str() {
        "p" => (vec![DVec3::ZERO], vec!["proton"], vec![up]),
        "n" => (vec![DVec3::ZERO], vec!["neutron"], vec![up]),
        // Parallel: both poles up (the lower north faces the upper south).
        "pp_par" => (vec![DVec3::new(0.0, 0.0, -s2), DVec3::new(0.0, 0.0, s2)], vec!["proton", "proton"], vec![up, up]),
        // North faces north.
        "pp_nn" => (vec![DVec3::new(0.0, 0.0, -s2), DVec3::new(0.0, 0.0, s2)], vec!["proton", "proton"], vec![up, down]),
        // South faces south.
        "pp_ss" => (vec![DVec3::new(0.0, 0.0, -s2), DVec3::new(0.0, 0.0, s2)], vec!["proton", "proton"], vec![down, up]),
        _ => panic!("CFM_M_LAYOUT: p | n | pp_par | pp_nn | pp_ss"),
    };
    println!("=== flow map: {layout}, sep {sep}, lambda {lambda}, n_a {n_a}, h {h}, dt {dt}, phi {phi:.3}, r_d {r_d}, hole {hole}, {snaps} snapshots every {every} steps ===");
    let ctx = GpuContext::new().expect("GPU");
    let cap = ((n_a * 4.0 / 3.0 * std::f64::consts::PI * r_d.powi(3) * 1.6) as u32).min(19_000_000);
    let mut config = two_body_config(&centers, &classes, n_a, lambda, cap, seed, dt);
    // CFM_M_SPEED: surface (default, doc §15.29) | rest (§15.27, recaptured by the
    // inner ball) | ejectc (§15.14).
    let mspeed = match env("CFM_M_SPEED", "surface").as_str() { "ejectc" => MazeSpeed::EjectC, "rest" => MazeSpeed::EjectRest, _ => MazeSpeed::EjectSurface };
    println!("maze speed {mspeed:?}");
    config.rough = RoughParams { intake: Intake::Envelope { p: 1.0 }, maze_speed: mspeed, through_hole: hole, polar_release_normal: env("CFM_M_POLAR", "axial") == "normal", ..config.rough };
    config.h = h;
    config.photon_fraction = phi;
    config.boundary = Boundary::SphereEscape { r_d };
    let mut engine = DsmcEngine::new(&ctx, config);
    for (b, cl) in engine.shell_bodies.iter_mut().zip(&classes) {
        b.depth = 0;
        if *cl == "neutron" {
            b.inner_signs = NEUTRON_INNER_SIGNS;
            b.tier1_signs = NEUTRON_INNER_SIGNS;
        }
    }
    apply_level_speeds(&mut engine);
    engine.set_body_poses(&ctx, &centers, Some(&orients));
    for _ in 0..fix {
        engine.step_flux(&ctx);
    }
    // Grid: R in [0, r_max), z in [-r_max, r_max), square cells DR.
    const DR: f64 = 0.1;
    let r_max = (r_d - 1.0).min(5.0);
    let nr = (r_max / DR).round() as usize;
    let nz = (2.0 * r_max / DR).round() as usize;
    let cells = nr * nz;
    let mut cnt = vec![0f64; cells];
    let mut s_vr = vec![0f64; cells];
    let mut s_vz = vec![0f64; cells];
    let mut s_vp = vec![0f64; cells];
    let mut tags = vec![[0f64; 3]; cells]; // ambient, body 0, body 1
    let mut plus = vec![0f64; cells];
    let mut t = Totals::default();
    // Body 0 ball kinematics, sampled every step: mean ball velocity along
    // the axis by the ball's own (distance from axis, height) in the body
    // frame. Emission inherits the ball velocity (EjectC: v = v_ball + c).
    let mut bk = vec![[0f64; 3]; 4 * 4]; // [r band 0-0.25-0.5-0.75-1][z band] -> (n, sum vz, sum |v|)
    for k in 0..snaps {
        for _ in 0..every {
            t.add(&engine.step_flux(&ctx));
            let sb = &engine.shell_bodies[0];
            let bs = sb.ball();
            let rel = sb.orient.inverse() * (bs.center - sb.center);
            let vl = sb.orient.inverse() * bs.velocity;
            let ir = ((rel.x * rel.x + rel.y * rel.y).sqrt() / 0.25).min(3.0) as usize;
            let iz = (((rel.z + 1.0) / 0.5).max(0.0) as usize).min(3);
            let c = &mut bk[ir * 4 + iz];
            c[0] += 1.0;
            c[1] += vl.z;
            c[2] += vl.length();
        }
        for p in engine.read_photons(&ctx) {
            // Skip dead photons and photons held in a tank (flag bit 20).
            if p.c[0] & FLAG_ALIVE == 0 || p.c[0] & (1u32 << 20) != 0 {
                continue;
            }
            let (x, y, z) = (p.a[0] as f64, p.a[1] as f64, p.a[2] as f64);
            let r = (x * x + y * y).sqrt();
            if r >= r_max || z.abs() >= r_max {
                continue;
            }
            let i = ((r / DR) as usize) * nz + ((z + r_max) / DR) as usize;
            let sp = p.d[3] as f64;
            let (vx, vy, vz) = (p.b[0] as f64 * sp, p.b[1] as f64 * sp, p.b[2] as f64 * sp);
            let (vr, vp) = if r > 1e-9 { ((x * vx + y * vy) / r, (x * vy - y * vx) / r) } else { (0.0, 0.0) };
            cnt[i] += 1.0;
            s_vr[i] += vr;
            s_vz[i] += vz;
            s_vp[i] += vp;
            let tag = ((p.c[0] >> 17) & 7) as usize;
            tags[i][tag.min(2)] += 1.0;
            if p.c[0] & FLAG_SPECIES_PLUS != 0 {
                plus[i] += 1.0;
            }
        }
        if k % 20 == 19 {
            println!("snapshot {}/{snaps}", k + 1);
        }
    }
    let steps = (snaps * every) as f64;
    let caps: Vec<String> = (0..centers.len()).map(|b| format!("{:.0}", t.env_capt[b] as f64 / (steps * dt))).collect();
    println!("envelope captures per unit time: {}", caps.join(", "));
    // Body 0 maze outcomes [escaped upright, escaped flipped, trapped
    // upright, trapped flipped, patch fallbacks] and the emission direction
    // latitude at the ball (before the ball's own velocity is added).
    println!("body 0 maze outcomes {:?}", t.maze);
    for ir in 0..4 {
        let row: Vec<String> = (0..4).map(|iz| { let c = bk[ir * 4 + iz]; if c[0] > 0.0 { format!("z{:+.2}..: share {:.3} vz {:+.2} |v| {:.2}", iz as f64 * 0.5 - 1.0, c[0] / (snaps * every) as f64, c[1] / c[0], c[2] / c[0]) } else { "-".into() } }).collect();
        println!("ball at axis dist {:.2}-{:.2}: {}", ir as f64 * 0.25, ir as f64 * 0.25 + 0.25, row.join(" | "));
    }
    let ml = &t.maze_lat;
    if ml.len() == 180 {
        let band = |lo: i32, hi: i32| -> i64 { (0..180).filter(|&k| { let l = k as i32 - 90; l >= lo && l < hi }).map(|k| ml[k]).sum() };
        let tot = ml.iter().sum::<i64>().max(1) as f64;
        println!(
            "body 0 emission latitude at the ball: N 60-90 {:.3}, N 30-60 {:.3}, eq +-30 {:.3}, S 30-60 {:.3}, S 60-90 {:.3}",
            band(60, 91) as f64 / tot, band(30, 60) as f64 / tot, band(-30, 30) as f64 / tot, band(-60, -30) as f64 / tot, band(-91, -60) as f64 / tot
        );
    }
    // Text digest: density and mean flow on a few probe points (2x2 cells).
    let vol = |ir: usize| 2.0 * std::f64::consts::PI * (ir as f64 + 0.5) * DR * DR * DR;
    let probe = |rr: f64, zz: f64| -> String {
        let ir = (rr / DR) as usize;
        let iz = ((zz + r_max) / DR) as usize;
        let (mut n, mut a, mut bz, mut c, mut v) = (0.0, 0.0, 0.0, 0.0, 0.0);
        for dr in 0..2 {
            for dz in 0..2 {
                let jr = (ir + dr).min(nr - 1);
                let i = jr * nz + (iz + dz).min(nz - 1);
                n += cnt[i];
                a += s_vr[i];
                bz += s_vz[i];
                c += s_vp[i];
                v += vol(jr);
            }
        }
        let nn = n.max(1.0);
        format!("({rr:.1},{zz:+.1}): dens {:.3} vR {:+.3} vz {:+.3} vphi {:+.3}", n / (snaps as f64 * v * n_a), a / nn, bz / nn, c / nn)
    };
    for (rr, zz) in [(0.0, 1.2), (0.0, 2.0), (0.0, 3.0), (0.0, -1.2), (0.0, -2.0), (0.6, 1.2), (1.0, 1.0), (1.3, 0.0), (2.0, 0.0), (3.0, 0.0), (1.5, 1.5), (2.5, 2.5)] {
        println!("probe {}", probe(rr, zz));
    }
    let mut js = String::new();
    js += &format!(
        "{{\"layout\":\"{layout}\",\"sep\":{sep},\"lambda\":{lambda},\"n_a\":{n_a},\"phi\":{phi},\"r_d\":{r_d},\"hole\":{hole},\"snaps\":{snaps},\"dr\":{DR},\"r_max\":{r_max},\"nr\":{nr},\"nz\":{nz},"
    );
    js += &format!(
        "\"bodies\":[{}],",
        (0..centers.len()).map(|b| format!("{{\"class\":\"{}\",\"z\":{},\"axis_z\":{}}}", classes[b], centers[b].z, (orients[b] * DVec3::Z).z)).collect::<Vec<_>>().join(",")
    );
    let col = |f: &dyn Fn(usize) -> f64| -> String { (0..cells).map(|i| format!("{:.4}", f(i))).collect::<Vec<_>>().join(",") };
    js += &format!("\"dens\":[{}],", col(&|i| cnt[i] / (snaps as f64 * vol(i / nz) * n_a)));
    js += &format!("\"count\":[{}],", col(&|i| cnt[i]));
    js += &format!("\"vr\":[{}],", col(&|i| s_vr[i] / cnt[i].max(1.0)));
    js += &format!("\"vz\":[{}],", col(&|i| s_vz[i] / cnt[i].max(1.0)));
    js += &format!("\"vphi\":[{}],", col(&|i| s_vp[i] / cnt[i].max(1.0)));
    js += &format!("\"tag_amb\":[{}],", col(&|i| tags[i][0] / cnt[i].max(1.0)));
    js += &format!("\"tag_b0\":[{}],", col(&|i| tags[i][1] / cnt[i].max(1.0)));
    js += &format!("\"tag_b1\":[{}],", col(&|i| tags[i][2] / cnt[i].max(1.0)));
    js += &format!("\"plus\":[{}]}}", col(&|i| plus[i] / cnt[i].max(1.0)));
    std::fs::write(&out, js).unwrap();
    println!("wrote {out} ({nr} x {nz} cells, index = iR * nz + iz)");
}

/// Report (doc §15.28): the core ball's path in the body frame for sub-c
/// outer levels (CPU only). The user's Spin-mode observation: with the
/// outermost level at c the core runs a wide loop that avoids whole
/// regions; at 0.5c or 0.05c the loop precesses round the axis like a fan
/// blade. Prints, per level-12 speed: where the ball spends its time
/// (distance from axis x height), its mean axial velocity there, and the
/// azimuth drift of the loop. Env: CFM_BL_SPEEDS (";"-separated sets of
/// "s9,s10,s11,s12" fractions), CFM_BL_T (time span, 400), CFM_BL_DT (0.002).
#[test]
#[ignore]
fn report_ball_loop() {
    use crate::baryon_shells::{ShellBody, PROTON_STATES};
    let env = |k: &str, d: &str| std::env::var(k).unwrap_or(d.to_string());
    let span: f64 = env("CFM_BL_T", "400").parse().unwrap();
    let dt: f64 = env("CFM_BL_DT", "0.002").parse().unwrap();
    for sp in env("CFM_BL_SPEEDS", "1,1,1,1;1,1,1,0.5;1,1,1,0.05").split(';') {
        let v: Vec<f64> = sp.split(',').map(|q| q.parse().unwrap()).collect();
        let mut b = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 1.0, PROTON_STATES[0]);
        b.speeds = [v[0], v[1], v[2], v[3]];
        // [r band 0.125 wide, 8 bands][z band 0.25 wide, 8 bands]: n, sum vz
        let mut grid = vec![[0f64; 2]; 64];
        let mut az_prev: Option<f64> = None;
        let mut az_unwrap = 0.0;
        let mut lz = 0.0; // mean angular momentum of the ball path about z
        let mut vz_r: [f64; 2] = [0.0; 2]; // correlation sum(vz * (r - mean r)) pieces
        let mut sum_r = 0.0;
        let mut n = 0.0;
        let steps = (span / dt) as usize;
        for _ in 0..steps {
            b.advance(dt);
            let bs = b.ball();
            let p = bs.center;
            let v = bs.velocity;
            let r = (p.x * p.x + p.y * p.y).sqrt();
            let ir = ((r / 0.125) as usize).min(7);
            let iz = (((p.z + 1.0) / 0.25).max(0.0) as usize).min(7);
            grid[ir * 8 + iz][0] += 1.0;
            grid[ir * 8 + iz][1] += v.z;
            let az = p.y.atan2(p.x);
            if let Some(a0) = az_prev {
                let mut d = az - a0;
                if d > std::f64::consts::PI {
                    d -= std::f64::consts::TAU;
                }
                if d < -std::f64::consts::PI {
                    d += std::f64::consts::TAU;
                }
                az_unwrap += d;
            }
            az_prev = Some(az);
            lz += p.x * v.y - p.y * v.x;
            vz_r[0] += v.z * r;
            vz_r[1] += v.z;
            sum_r += r;
            n += 1.0;
        }
        let cov = vz_r[0] / n - (vz_r[1] / n) * (sum_r / n);
        println!("=== level 9-12 speeds {v:?} (fraction of c): mean axial velocity {:+.3}, cov(v_z, distance from axis) {:+.4}, mean L_z of path {:+.3}, net azimuth turns {:+.1} in {span} time units", vz_r[1] / n, cov, lz / n, az_unwrap / std::f64::consts::TAU);
        println!("  time share by distance from axis (rows, 0.125 bands) x height (cols -1..1, 0.25 bands); mean v_z in brackets");
        for ir in 0..8 {
            let row: Vec<String> = (0..8)
                .map(|iz| {
                    let c = grid[ir * 8 + iz];
                    if c[0] > 0.0 { format!("{:5.3}({:+.2})", c[0] / n, c[1] / c[0]) } else { "     -      ".to_string() }
                })
                .collect();
            println!("  r {:.3}-{:.3}: {}", ir as f64 * 0.125, ir as f64 * 0.125 + 0.125, row.join(" "));
        }
    }
}
