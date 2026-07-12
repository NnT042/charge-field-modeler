//! Scenario harness for headless physics iteration.
//!
//! Drives `AtomCore` through canonical setups (hydrogen capture, H₂ bonding,
//! alpha/helium presets) and reduces trajectories to summary metrics that
//! tests can assert on and `atom_lab` (bin) can print. This is the fast
//! feedback loop: `cargo test --release --manifest-path rust/Cargo.toml`.

use crate::atom_core::{load_histogram_csv, config_dir, AtomCore};
use glam::DVec3;

/// Default radius below which the orbiter counts as captured (sim units;
/// contact for p+e is 1.3, hydrogen spawn is ~2.06).
pub const CAPTURE_RADIUS: f64 = 2.0;
/// Radius beyond which a previously-captured orbiter counts as escaped.
pub const ESCAPE_RADIUS: f64 = 12.0;
/// |v_rad| above this, inside contact+0.2, counts as a wall bounce.
pub const BOUNCE_SPEED: f64 = 0.5;

#[derive(Debug, Clone, Default)]
pub struct ScenarioMetrics {
    /// First sim time with r < CAPTURE_RADIUS, if any.
    pub capture_time: Option<f64>,
    /// Orbit stats over the post-settle window.
    pub mean_orbit_r: f64,
    pub orbit_r_stddev: f64,
    pub v_tan_mean: f64,
    pub v_rad_mean: f64,
    pub theta_pole_mean_deg: f64,
    /// Extremes over the whole run.
    pub min_r: f64,
    pub max_r: f64,
    pub final_r: f64,
    /// r exceeded ESCAPE_RADIUS after having been captured.
    pub escaped: bool,
    /// Any POST-SETTLE sample with |v_rad| > BOUNCE_SPEED within contact+0.2
    /// (persistent wall chaos; the initial touchdown is expected to bounce).
    pub bounced: bool,
}

/// A core with the three standard profiles registered from the shipped
/// histogram CSVs, matching what atom_mode.gd loads at runtime.
pub fn standard_core() -> AtomCore {
    let dir = config_dir();
    let mut core = AtomCore::new();
    let proton = load_histogram_csv(&dir.join("histogram_proton.csv"));
    core.register_profile("proton", 1.0, 1.0, &proton);
    let neutron = load_histogram_csv(&dir.join("histogram_neutron.csv"));
    core.register_profile("neutron", 1.0, 1.0, &neutron);
    let electron = load_histogram_csv(&dir.join("histogram_electron.csv"));
    core.register_profile("electron", 1.0 / 1836.0, 0.3, &electron);
    core.running = true;
    core
}

/// Spawn a hydrogen atom: proton at `base` with pole along `axis`, electron
/// released near the chosen pole with a small tangential velocity.
/// `spin_sign` flips the proton's axial spin (H₂ combo matrix needs both).
/// Returns (proton_id, electron_id).
pub fn spawn_hydrogen(
    core: &mut AtomCore,
    base: DVec3,
    axis: DVec3,
    electron_on_positive_pole: bool,
    spin_sign: f64,
) -> (usize, usize) {
    let axis = axis.normalize();
    let p_id = core.profile_id_by_name("proton").expect("proton profile");
    let e_id = core.profile_id_by_name("electron").expect("electron profile");

    let p_spin = crate::atom_core::default_spin_rate("proton") * spin_sign;
    let proton = core
        .spawn_particle_ex(p_id, base, DVec3::ZERO, axis, p_spin)
        .expect("spawn proton");

    // Mirror atom_mode.gd's _spawn_hydrogen(): 0.5 off-axis, 2.0 along the
    // pole, tangential kick of 0.5 perpendicular to both.
    let pole_dir = if electron_on_positive_pole { axis } else { -axis };
    let (right, forward) = crate::atom_core::build_frame(pole_dir);
    let e_pos = base + right * 0.5 + pole_dir * 2.0;
    let e_vel = forward * 0.5;
    let e_spin = crate::atom_core::default_spin_rate("electron") * spin_sign;
    let electron = core
        .spawn_particle_ex(e_id, e_pos, e_vel, pole_dir, e_spin)
        .expect("spawn electron");

    (proton, electron)
}

/// Run `n_steps`, sampling the (center, orbiter) pair every `sample_every`
/// steps, and reduce to metrics. Orbit statistics (means/stddev) are taken
/// over the samples after `settle_frac` of the run so transients don't
/// pollute them; extremes and capture/escape/bounce cover the whole run.
pub fn run_pair(
    core: &mut AtomCore,
    center: usize,
    orbiter: usize,
    n_steps: usize,
    sample_every: usize,
    settle_frac: f64,
) -> ScenarioMetrics {
    let sample_every = sample_every.max(1);
    let contact = core.contact_distance(center, orbiter);
    let n_samples = n_steps / sample_every;
    let settle_start = ((n_samples as f64) * settle_frac.clamp(0.0, 0.95)) as usize;

    let mut m = ScenarioMetrics {
        min_r: f64::MAX,
        ..Default::default()
    };
    let mut window: Vec<(f64, f64, f64, f64)> = Vec::new(); // r, v_tan, v_rad, theta

    for s in 0..n_samples {
        core.step_n(sample_every);

        let d = core.particles[orbiter].position - core.particles[center].position;
        let r = d.length();
        let d_hat = d / r.max(1e-12);
        let v_rel = core.particles[orbiter].velocity - core.particles[center].velocity;
        let v_rad = v_rel.dot(d_hat);
        let v_tan = (v_rel - d_hat * v_rad).length();
        let cos_theta = core.particles[center].pole_axis().dot(d_hat);
        let theta_deg = cos_theta.abs().clamp(0.0, 1.0).acos().to_degrees();

        m.min_r = m.min_r.min(r);
        m.max_r = m.max_r.max(r);
        m.final_r = r;

        if m.capture_time.is_none() && r < CAPTURE_RADIUS {
            m.capture_time = Some(core.time);
        }
        if m.capture_time.is_some() && r > ESCAPE_RADIUS {
            m.escaped = true;
        }
        if s >= settle_start {
            if r < contact + 0.2 && v_rad.abs() > BOUNCE_SPEED {
                m.bounced = true;
            }
            window.push((r, v_tan, v_rad, theta_deg));
        }
    }

    if !window.is_empty() {
        let n = window.len() as f64;
        m.mean_orbit_r = window.iter().map(|w| w.0).sum::<f64>() / n;
        m.v_tan_mean = window.iter().map(|w| w.1).sum::<f64>() / n;
        m.v_rad_mean = window.iter().map(|w| w.2).sum::<f64>() / n;
        m.theta_pole_mean_deg = window.iter().map(|w| w.3).sum::<f64>() / n;
        let var = window
            .iter()
            .map(|w| (w.0 - m.mean_orbit_r).powi(2))
            .sum::<f64>()
            / n;
        m.orbit_r_stddev = var.sqrt();
    }
    if m.min_r == f64::MAX {
        m.min_r = 0.0;
    }
    m
}

/// Free-space radial equilibrium radius at polar angle `theta` (radians)
/// for an orbiter around a center, poles parallel:
///   attraction (G_q + I_q·A_c²(θ))/r²  =  repulsion C_q·E_c(θ)·R_o(θ)/r⁴
///   ⇒  r_eq = sqrt( C_q·E_c(θ)·R_o(θ) / (G_q + I_q·A_c²(θ)) )
/// Where r_eq < contact, the orbiter rests ON the boundary (diatom.pdf:
/// captured electrons sit at the nuclear boundary); where r_eq > contact,
/// the 1/r⁴ emission wall stands it off.
pub fn equilibrium_radius(
    core: &AtomCore,
    center_profile: usize,
    orbiter_profile: usize,
    theta: f64,
) -> f64 {
    let cq = &core.couplings;
    let cos_t = theta.cos();
    let c = &core.profiles[center_profile];
    let o = &core.profiles[orbiter_profile];
    let e_c = c.emission.sample(cos_t);
    let a_c = c.absorption.sample(cos_t);
    let r_o = o.absorption.sample(cos_t);
    let denom = cq.g_q + cq.intake * a_c * a_c;
    if denom <= 0.0 {
        return f64::INFINITY;
    }
    (cq.c_q * e_c * r_o / denom).max(0.0).sqrt()
}

/// Smallest polar angle (degrees) where emission reaches `threshold` of peak —
/// the edge of the polar tunnel.
pub fn tunnel_edge_angle_deg(table: &crate::atom_core::EmissionTable, threshold: f64) -> f64 {
    for tenth_deg in 0..=900 {
        let theta = (tenth_deg as f64 / 10.0).to_radians();
        if table.sample(theta.cos()) >= threshold {
            return theta.to_degrees();
        }
    }
    90.0
}

// ── H₂ bonding scenarios ──────────────────────────────────────────────────

/// Wall-riding orbit geometry the M1 physics settles into (see
/// hydrogen_orbit_stable_long_run): r = contact (1.3), θ ≈ 11° from the
/// pole, v_tan = COROT_V_MAX circling the pole axis.
const RIDE_R: f64 = 1.3;
const RIDE_THETA_DEG: f64 = 11.0;

/// Spawn a PRE-FORMED hydrogen atom: proton plus electron already in the
/// wall-riding polar orbit, so bonding runs don't wait out the capture
/// transient. `electron_dir` is the unit direction from the proton to the
/// occupied pole. Returns (proton_id, electron_id).
pub fn spawn_formed_hydrogen(
    core: &mut AtomCore,
    p_pos: DVec3,
    proton_pole: DVec3,
    electron_dir: DVec3,
    spin_sign: f64,
) -> (usize, usize) {
    let p_id = core.profile_id_by_name("proton").expect("proton profile");
    let e_id = core.profile_id_by_name("electron").expect("electron profile");
    let proton_pole = proton_pole.normalize();
    let electron_dir = electron_dir.normalize();

    let p_spin = crate::atom_core::default_spin_rate("proton") * spin_sign;
    let proton = core
        .spawn_particle_ex(p_id, p_pos, DVec3::ZERO, proton_pole, p_spin)
        .expect("spawn proton");

    // Electron on the wall at θ from the occupied pole, circling the axis
    // at the corotation speed in the direction the proton's spin drags it.
    let theta = RIDE_THETA_DEG.to_radians();
    let (lat_dir, _) = crate::atom_core::build_frame(electron_dir);
    let e_offset = electron_dir * (RIDE_R * theta.cos()) + lat_dir * (RIDE_R * theta.sin());
    let e_pos = p_pos + e_offset;
    // v_corot direction = (spin vector) × (offset from axis)
    let spin_vec = proton_pole * spin_sign;
    let lateral = lat_dir * (RIDE_R * theta.sin());
    let v_dir = spin_vec.cross(lateral).normalize();
    let e_vel = v_dir * crate::atom_core::COROT_V_MAX;
    let e_spin = crate::atom_core::default_spin_rate("electron") * spin_sign;
    let electron = core
        .spawn_particle_ex(e_id, e_pos, e_vel, electron_dir, e_spin)
        .expect("spawn electron");

    (proton, electron)
}

/// Spawn two hydrogen atoms stacked on the Y axis (the bond in Mathis is
/// POLAR — diatom.pdf: atoms align their poles; the bond forms when the
/// facing poles are bare and the electrons sit on the OUTER poles).
/// Returns [proton_a, electron_a, proton_b, electron_b].
/// Combo encoding for the 8-way matrix: spin pair (a,b) ∈ {+,−}² ×
/// electrons outside/between.
pub fn spawn_h2(
    core: &mut AtomCore,
    electrons_outside: bool,
    spin_a: f64,
    spin_b: f64,
    separation: f64,
) -> [usize; 4] {
    let a_pos = DVec3::new(0.0, -separation / 2.0, 0.0);
    let b_pos = DVec3::new(0.0, separation / 2.0, 0.0);
    // Both protons poles +Y (parallel; anti-parallel cases are covered by
    // the spin signs since the histograms are bilaterally symmetric).
    let (a_e_dir, b_e_dir) = if electrons_outside {
        (-DVec3::Y, DVec3::Y) // away from the partner atom
    } else {
        (DVec3::Y, -DVec3::Y) // sandwiched between the protons
    };
    let (pa, ea) = spawn_formed_hydrogen(core, a_pos, DVec3::Y, a_e_dir, spin_a);
    let (pb, eb) = spawn_formed_hydrogen(core, b_pos, DVec3::Y, b_e_dir, spin_b);
    [pa, ea, pb, eb]
}

#[derive(Debug, Clone, Default)]
pub struct BondMetrics {
    pub initial_d: f64,
    /// Proton–proton distance stats over the post-settle window.
    pub mean_d: f64,
    pub d_stddev: f64,
    pub min_d: f64,
    pub max_d: f64,
    pub final_d: f64,
    /// Both electrons stayed within 2.0 of their own proton.
    pub electrons_retained: bool,
}

impl BondMetrics {
    /// Bonded ⇔ settled at molecular standoff, vibrating (M5 plan: the bond
    /// SHOULD vibrate — rigid means pressure too high) but bounded, never
    /// collapsing to nuclear contact (min_d > 2.2 ⇒ no ambient fusion),
    /// and still holding at the end.
    pub fn bonded(&self) -> bool {
        self.mean_d >= 2.5
            && self.mean_d <= 6.0
            && self.d_stddev / self.mean_d < 0.35
            && self.min_d > 2.2
            && self.final_d < 7.0
    }

    /// Repelled ⇔ driven out and parked/oscillating clearly beyond bond
    /// range (the stoppered-channel standoff), or monotonically leaving.
    pub fn repelled(&self) -> bool {
        self.mean_d > 7.0
            || (self.final_d > 1.5 * self.initial_d && self.final_d >= self.max_d - 1e-9)
    }
}

/// Run and reduce an H₂ (or any two-atom) run to bond metrics.
pub fn run_bond(
    core: &mut AtomCore,
    ids: [usize; 4],
    n_steps: usize,
    sample_every: usize,
    settle_frac: f64,
) -> BondMetrics {
    let [pa, ea, pb, eb] = ids;
    let sample_every = sample_every.max(1);
    let n_samples = n_steps / sample_every;
    let settle_start = ((n_samples as f64) * settle_frac.clamp(0.0, 0.95)) as usize;

    let mut m = BondMetrics {
        initial_d: core.pair_distance(pa, pb),
        min_d: f64::MAX,
        electrons_retained: true,
        ..Default::default()
    };
    let mut window: Vec<f64> = Vec::new();

    for s in 0..n_samples {
        core.step_n(sample_every);
        let d = core.pair_distance(pa, pb);
        m.min_d = m.min_d.min(d);
        m.max_d = m.max_d.max(d);
        m.final_d = d;
        if core.pair_distance(pa, ea) > 2.0 || core.pair_distance(pb, eb) > 2.0 {
            m.electrons_retained = false;
        }
        if s >= settle_start {
            window.push(d);
        }
    }
    if !window.is_empty() {
        let n = window.len() as f64;
        m.mean_d = window.iter().sum::<f64>() / n;
        let var = window.iter().map(|d| (d - m.mean_d).powi(2)).sum::<f64>() / n;
        m.d_stddev = var.sqrt();
    }
    if m.min_d == f64::MAX {
        m.min_d = 0.0;
    }
    m
}

/// One sampled trajectory frame (all particles).
pub struct TrajectoryRow {
    pub t: f64,
    pub positions: Vec<DVec3>,
    pub velocities: Vec<DVec3>,
}

/// Run and record full trajectories for offline inspection (atom_lab CSV).
pub fn run_recording(
    core: &mut AtomCore,
    n_steps: usize,
    sample_every: usize,
) -> Vec<TrajectoryRow> {
    let sample_every = sample_every.max(1);
    let n_samples = n_steps / sample_every;
    let mut rows = Vec::with_capacity(n_samples);
    for _ in 0..n_samples {
        core.step_n(sample_every);
        rows.push(TrajectoryRow {
            t: core.time,
            positions: core.particles.iter().map(|p| p.position).collect(),
            velocities: core.particles.iter().map(|p| p.velocity).collect(),
        });
    }
    rows
}

// ── Scenario tests ────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// M0 baseline canary: with the current force model the electron must
    /// fall into the polar intake (capture) and stay bound. This pinned the
    /// user-observed behavior of session 27 BEFORE the M1 physics changes.
    #[test]
    fn hydrogen_capture() {
        let mut core = standard_core();
        let (p, e) = spawn_hydrogen(&mut core, DVec3::ZERO, DVec3::Y, true, 1.0);
        let m = run_pair(&mut core, p, e, 100_000, 50, 0.4);
        assert!(
            m.capture_time.is_some(),
            "electron never captured: {m:?}"
        );
        assert!(!m.escaped, "electron escaped after capture: {m:?}");
    }

    /// M1 gate: the captured electron must ORBIT the pole (circle the drain),
    /// not sit dead on the contact wall. Radius at the nuclear boundary
    /// (contact = 1.3), sustained tangential motion at the corotation speed,
    /// inside the polar tunnel, no escape, no persistent bouncing.
    #[test]
    fn hydrogen_orbit_stable_long_run() {
        let mut core = standard_core();
        let (p, e) = spawn_hydrogen(&mut core, DVec3::ZERO, DVec3::Y, true, 1.0);
        let m = run_pair(&mut core, p, e, 500_000, 100, 0.4);
        assert!(m.capture_time.is_some(), "no capture: {m:?}");
        assert!(!m.escaped, "escaped: {m:?}");
        assert!(!m.bounced, "persistent contact-wall bouncing: {m:?}");
        assert!(
            m.mean_orbit_r >= 1.1 && m.mean_orbit_r <= 1.8,
            "orbit radius off nuclear boundary: {m:?}"
        );
        assert!(
            m.orbit_r_stddev / m.mean_orbit_r < 0.15,
            "orbit not stable: {m:?}"
        );
        assert!(m.v_tan_mean > 0.05, "no sustained tangential motion: {m:?}");
        // Corotation drag should hold v_tan at the (saturated) vortex speed.
        assert!(
            (m.v_tan_mean - crate::atom_core::COROT_V_MAX).abs() < 0.05,
            "v_tan should sit at corotation ceiling: {m:?}"
        );
        // And the orbit must live inside the polar tunnel, not the equator.
        let edge = tunnel_edge_angle_deg(&core.profiles[0].emission, 0.25);
        assert!(
            m.theta_pole_mean_deg < edge,
            "orbit drifted out of the polar tunnel (θ={} edge={edge}): {m:?}",
            m.theta_pole_mean_deg
        );
    }

    /// M1 gate: the locked constants must be consistent with the radial
    /// equilibrium equation and the wall-riding orbit it predicts:
    /// r_eq(θ) < contact inside the tunnel (electron rests ON the nuclear
    /// boundary — diatom.pdf), r_eq(θ) > contact at mid-latitudes (the
    /// 1/r⁴ emission wall confines it to the pole channel), and the
    /// simulated steady state lands within 15% of the prediction.
    #[test]
    fn derived_constants_equilibrium() {
        let mut core = standard_core();
        let p_prof = core.profile_id_by_name("proton").unwrap();
        let e_prof = core.profile_id_by_name("electron").unwrap();

        // The loaded proton histogram puts E ≥ 0.25·peak at ~46° — wider than
        // the ~25-28° "tunnel half-angle" quoted in PHYSICS_REFERENCE.md
        // (different threshold convention); the orbit must simply live well
        // inside it, which hydrogen_orbit_stable_long_run asserts.
        let edge = tunnel_edge_angle_deg(&core.profiles[p_prof].emission, 0.25);
        assert!(
            (15.0..=55.0).contains(&edge),
            "polar tunnel edge out of plausible band: got {edge}"
        );

        let contact = 1.3; // r_p(1.0) + r_e(0.3)

        // Find the wall-riding crossover θ*: the largest angle from the pole
        // where the free equilibrium still sits inside the contact wall.
        let mut crossover_deg = 0.0;
        for tenth in 0..=900 {
            let theta_deg = tenth as f64 / 10.0;
            let r_eq = equilibrium_radius(&core, p_prof, e_prof, theta_deg.to_radians());
            if r_eq < contact {
                crossover_deg = theta_deg;
            } else if crossover_deg > 0.0 {
                break;
            }
        }
        // The wall-riding channel must exist and be a polar feature.
        assert!(
            crossover_deg > 5.0 && crossover_deg < edge,
            "wall-riding crossover θ*={crossover_deg}° should be inside the tunnel (edge={edge}°)"
        );
        // Mid-latitudes: the emission wall stands the electron off.
        for theta_deg in [55.0, 65.0, 75.0] {
            let r_eq = equilibrium_radius(&core, p_prof, e_prof, f64::to_radians(theta_deg));
            assert!(
                r_eq > contact,
                "emission wall at θ={theta_deg}° should stand off beyond contact: r_eq={r_eq}"
            );
        }

        // The simulation must agree: steady wall-riding orbit at the contact
        // boundary, inside the crossover channel.
        let (p, e) = spawn_hydrogen(&mut core, DVec3::ZERO, DVec3::Y, true, 1.0);
        let m = run_pair(&mut core, p, e, 200_000, 100, 0.5);
        let predicted = contact; // wall-riding prediction
        assert!(
            (m.mean_orbit_r - predicted).abs() / predicted < 0.15,
            "sim orbit r={} vs predicted {predicted}: {m:?}",
            m.mean_orbit_r
        );
        assert!(
            m.theta_pole_mean_deg < crossover_deg,
            "orbit θ={}° should sit inside the wall-riding channel θ*={crossover_deg}°: {m:?}",
            m.theta_pole_mean_deg
        );
    }

    /// M2 gate: of the 8 spin/pole combinations (4 spin pairings × electron
    /// placement), exactly the 4 with electrons on the OUTER poles bond and
    /// the 4 with electrons between the protons repel. The discriminant
    /// emerges from the stream-cushion + stoppered-vortex forces — there is
    /// no coded bonding rule anywhere.
    #[test]
    fn h2_bond_matrix() {
        let mut bonds = 0;
        let mut repels = 0;
        for combo in 0..8usize {
            let mut core = standard_core();
            let outside = combo < 4;
            let spin_a = if combo & 0b01 == 0 { 1.0 } else { -1.0 };
            let spin_b = if combo & 0b10 == 0 { 1.0 } else { -1.0 };
            let ids = spawn_h2(&mut core, outside, spin_a, spin_b, 6.0);
            let m = run_bond(&mut core, ids, 250_000, 200, 0.4);
            assert!(
                m.electrons_retained,
                "combo {combo}: an atom lost its electron: {m:?}"
            );
            if outside {
                assert!(
                    m.bonded() && !m.repelled(),
                    "combo {combo} (electrons outside) should bond: {m:?}"
                );
                bonds += 1;
            } else {
                assert!(
                    m.repelled() && !m.bonded(),
                    "combo {combo} (electrons between) should repel: {m:?}"
                );
                repels += 1;
            }
        }
        assert_eq!((bonds, repels), (4, 4));
    }

    /// M3 gate: an isolated alpha stays rigid (constituent geometry exact),
    /// conserves COM momentum under a kick (no external forces), and never
    /// goes non-finite.
    #[test]
    fn alpha_holds_and_conserves() {
        let mut core = standard_core();
        let gid = core
            .spawn_preset("alpha", DVec3::ZERO, DVec3::new(0.3, 0.1, 0.0), DVec3::Y)
            .expect("alpha preset");
        let members = core.groups[gid].members.clone();
        let initial_dists: Vec<f64> = pair_dists(&core, &members);
        let p0 = momentum(&core);

        core.step_n(100_000);

        for p in &core.particles {
            assert!(
                p.position.is_finite() && p.velocity.is_finite(),
                "non-finite state"
            );
        }
        let final_dists = pair_dists(&core, &members);
        for (a, b) in initial_dists.iter().zip(&final_dists) {
            assert!(
                (a - b).abs() < 1e-9,
                "rigidity violated: pair distance {a} -> {b}"
            );
        }
        let p1 = momentum(&core);
        assert!(
            (p0 - p1).length() < 1e-9,
            "momentum not conserved: {p0:?} -> {p1:?}"
        );
    }

    /// Two free protons, poles parallel, side by side: the disc-collision
    /// cushion must stand them off — parked and parallel, NOT gravitating
    /// into a mutual orbit ("deranged ornaments", session 28 user report)
    /// and NOT fusing. Ambient-pressure fusion is a star's job.
    #[test]
    fn two_protons_stand_off() {
        let mut core = standard_core();
        let pid = core.profile_id_by_name("proton").unwrap();
        core.spawn_particle(pid, DVec3::new(-3.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
        core.spawn_particle(pid, DVec3::new(3.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
        let m = run_pair(&mut core, 0, 1, 200_000, 100, 0.4);
        assert!(
            m.mean_orbit_r > 4.0 && m.mean_orbit_r < 6.5,
            "protons should park at the disc-collision standoff: {m:?}"
        );
        assert!(
            m.v_tan_mean < 0.02,
            "protons should not circle each other: {m:?}"
        );
        assert!(
            (m.theta_pole_mean_deg - 90.0).abs() < 5.0,
            "poles should stay parallel (equator-facing): {m:?}"
        );
        assert!(m.min_r > 2.2, "protons must not fuse at ambient: {m:?}");
    }

    /// Stretch presets (C/N/O/Ne/Ar) smoke test: spawn, run, stay rigid
    /// and finite. Rigidity invariant SPLIT (session-31 addendum A7):
    /// with independent per-alpha roll, cross-alpha member distances
    /// legitimately change under `advance_display` (post azimuths
    /// decorrelate as each alpha rolls at its own random phase) — so
    /// PHYSICS rigidity (`step_n`, no phases advance in `step`) is
    /// checked all-pairs over a probe, and DISPLAY rigidity
    /// (`advance_display`) is checked only WITHIN one alpha (which stays
    /// one rigid piece).
    #[test]
    fn heavier_presets_smoke() {
        for name in ["carbon", "tri_alpha", "nitrogen", "oxygen", "neon", "argon"] {
            let mut core = standard_core();
            let gid = core
                .spawn_preset(name, DVec3::ZERO, DVec3::new(0.1, 0.0, 0.05), DVec3::Y)
                .unwrap_or_else(|| panic!("{name} preset missing"));
            let members = core.groups[gid].members.clone();
            let expected = match name {
                // Session-32 carbon: 2 core alphas (8) + 2 plug pairs (4).
                "carbon" => 12,
                "tri_alpha" => 12, // bare 3-stack (harness structure)
                "nitrogen" => 14,
                "oxygen" => 16,
                "neon" => 20,   // center alpha + 4 carousel alphas
                "argon" => 36,  // 9 alphas: axial line of 5 + 4 carousel
                _ => unreachable!(),
            };
            assert_eq!(members.len(), expected, "{name} constituent count");

            // Polar plug geometry (session-29), rewritten in AlphaUnit
            // terms (session-31 addendum A7) and FUSED-pair terms
            // (session 34, see plug_pair): LONE plugs (nitrogen) are
            // single-member alphas — protons edge-on (rest_axis·Y ≈ 0,
            // disc feeds the hole), neutrons pole-on-axis; PAIRED plugs
            // (carbon, oxygen) are one 2-member fused alpha per pole —
            // proton pole toward its partner, neutron pole on the stack
            // axis, side by side at the same height. All ride the
            // carousel.
            let member_particle =
                |core: &AtomCore, gid: usize, ai: usize, k: usize| -> usize {
                    core.groups[gid].members[core.groups[gid].alphas[ai].members[k]]
                };
            let profile_name_of = |core: &AtomCore, pid: usize| -> String {
                core.profiles[core.particles[pid].profile_id].name.clone()
            };
            let lone_plugs: Vec<usize> = (0..core.groups[gid].alphas.len())
                .filter(|&ai| core.groups[gid].alphas[ai].members.len() == 1)
                .collect();
            let pair_plugs: Vec<usize> = (0..core.groups[gid].alphas.len())
                .filter(|&ai| core.groups[gid].alphas[ai].members.len() == 2)
                .collect();
            for &ai in &lone_plugs {
                let a = &core.groups[gid].alphas[ai];
                assert!(a.orbits_core, "plug alpha should ride the carousel");
                let m = member_particle(&core, gid, ai, 0);
                match profile_name_of(&core, m).as_str() {
                    "proton" => assert!(
                        a.rest_axis.dot(DVec3::Y).abs() < 1e-9,
                        "plug proton must be edge-on to the stack"
                    ),
                    "neutron" => assert!(
                        a.rest_axis.dot(DVec3::Y).abs() > 1.0 - 1e-9,
                        "plug neutron must keep its pole on the stack axis"
                    ),
                    other => panic!("unexpected plug profile {other}"),
                }
            }
            for &ai in &pair_plugs {
                assert!(
                    core.groups[gid].alphas[ai].orbits_core,
                    "plug pair should ride the carousel"
                );
                let (m0, m1) = (
                    member_particle(&core, gid, ai, 0),
                    member_particle(&core, gid, ai, 1),
                );
                let (mp, mn) = if profile_name_of(&core, m0) == "proton" {
                    (m0, m1)
                } else {
                    (m1, m0)
                };
                assert_eq!(profile_name_of(&core, mp), "proton");
                assert_eq!(profile_name_of(&core, mn), "neutron");
                let pp = core.particles[mp].position;
                let pn = core.particles[mn].position;
                assert!(
                    (pp.y - pn.y).abs() < 1e-9,
                    "pair must sit at the same height (side by side)"
                );
                assert!(
                    (pp - pn).length() > 0.5,
                    "pair members must sit beside each other, not overlap"
                );
                let to_partner = (pn - pp).normalize();
                assert!(
                    core.particles[mp].pole_axis().dot(to_partner) > 0.99,
                    "plug proton pole should point at its fused neutron"
                );
                assert!(
                    core.particles[mn].pole_axis().dot(DVec3::Y).abs() > 1.0 - 1e-9,
                    "plug neutron pole must stay on the stack axis"
                );
            }
            match name {
                "nitrogen" => {
                    assert_eq!(lone_plugs.len(), 2, "nitrogen has 2 lone plugs");
                    assert!(pair_plugs.is_empty());
                }
                "carbon" => {
                    assert_eq!(pair_plugs.len(), 2, "carbon has 2 fused plug pairs");
                    assert!(lone_plugs.is_empty());
                }
                "oxygen" => {
                    assert_eq!(pair_plugs.len(), 2, "oxygen has 2 fused plug pairs");
                    assert!(lone_plugs.is_empty());
                }
                _ => {}
            }

            // Physics rigidity: all-pairs distances over a probe of the
            // first 8 members unchanged under step_n (phases don't
            // advance in step).
            let probe = &members[..8.min(members.len())];
            let before_step = pair_dists(&core, probe);
            // Big carousel presets (Ne/Ar) pay O(n²) per step — a shorter
            // run keeps the ~10 s test gate while still proving rigidity.
            let steps = if expected > 16 { 10_000 } else { 50_000 };
            core.step_n(steps);
            let after_step = pair_dists(&core, probe);
            for (a, b) in before_step.iter().zip(&after_step) {
                assert!((a - b).abs() < 1e-9, "{name} physics rigidity violated");
            }

            // Display rigidity: WITHIN one alpha, all pair distances
            // unchanged under advance_display (an alpha is rigid).
            let alpha0_ids: Vec<usize> = core.groups[gid].alphas[0]
                .members
                .iter()
                .map(|&k| core.groups[gid].members[k])
                .collect();
            let before_disp = pair_dists(&core, &alpha0_ids);
            core.advance_display(0.8);
            let after_disp = pair_dists(&core, &alpha0_ids);
            for (a, b) in before_disp.iter().zip(&after_disp) {
                assert!(
                    (a - b).abs() < 1e-9,
                    "{name} display rigidity violated within alpha 0"
                );
            }

            for p in &core.particles {
                assert!(p.position.is_finite(), "{name}: non-finite state");
            }
        }
    }

    /// Bond detection matches the h2 matrix: electrons-outside reads as a
    /// molecular bond, electron-between is a stoppered channel (no bond),
    /// and the equator-facing two-proton standoff is repulsion, not a bond.
    #[test]
    fn molecular_bond_detection_matches_h2_matrix() {
        let mut core = standard_core();
        let ids = spawn_h2(&mut core, true, 1.0, 1.0, 4.7);
        let bonds = core.molecular_bonds();
        assert!(
            bonds.contains(&(ids[0], ids[2])) || bonds.contains(&(ids[2], ids[0])),
            "electrons-outside H2 should read as bonded: {bonds:?}"
        );

        let mut core = standard_core();
        spawn_h2(&mut core, false, 1.0, 1.0, 4.7);
        let bonds = core.molecular_bonds();
        assert!(
            bonds.is_empty(),
            "electron-between = stoppered channel, no bond: {bonds:?}"
        );

        let mut core = standard_core();
        let p_id = core.profile_id_by_name("proton").unwrap();
        core.spawn_particle(p_id, DVec3::new(-2.95, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
        core.spawn_particle(p_id, DVec3::new(2.95, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
        assert!(
            core.molecular_bonds().is_empty(),
            "side-by-side protons at standoff are repelling, not bonded"
        );
    }

    fn pair_dists(core: &AtomCore, ids: &[usize]) -> Vec<f64> {
        let mut out = Vec::new();
        for (n, &a) in ids.iter().enumerate() {
            for &b in &ids[n + 1..] {
                out.push(core.pair_distance(a, b));
            }
        }
        out
    }

    fn momentum(core: &AtomCore) -> DVec3 {
        core.particles
            .iter()
            .map(|p| p.velocity * core.profiles[p.profile_id].mass)
            .sum()
    }

    /// M3 gate: a free electron released over the alpha's top proton pole is
    /// captured into the same wall-riding orbit as hydrogen — constituent-
    /// level nearfield working through the rigid group.
    #[test]
    fn alpha_captures_electron() {
        let mut core = standard_core();
        let gid = core
            .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        // Members: [proton_bottom, proton_top, neutron, neutron]
        let top_proton = core.groups[gid].members[1];
        let e_id = core.profile_id_by_name("electron").unwrap();
        let e = core
            .spawn_particle(
                e_id,
                DVec3::new(0.5, 3.3, 0.0), // 2.0 above the top proton at (0,1.3,0)
                DVec3::new(0.0, 0.0, 0.4),
                DVec3::Y,
            )
            .unwrap();
        let m = run_pair(&mut core, top_proton, e, 300_000, 100, 0.5);
        assert!(m.capture_time.is_some(), "no capture on alpha pole: {m:?}");
        assert!(!m.escaped, "escaped: {m:?}");
        assert!(
            m.mean_orbit_r >= 1.1 && m.mean_orbit_r <= 2.0,
            "electron should ride near the top proton's boundary: {m:?}"
        );
        assert!(
            m.orbit_r_stddev / m.mean_orbit_r < 0.2,
            "orbit not stable: {m:?}"
        );
    }

    /// M3 gate: helium — alpha plus two electrons riding the two outer
    /// proton poles — holds over a long run: both electrons stay bound to
    /// their poles, the nucleus stays put (no self-propulsion).
    #[test]
    fn helium_stable() {
        let mut core = standard_core();
        let gid = core
            .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        let bottom_proton = core.groups[gid].members[0];
        let top_proton = core.groups[gid].members[1];
        let e_id = core.profile_id_by_name("electron").unwrap();
        // Pre-placed wall-riding electrons on the two outer poles, circling
        // opposite ways to match the shared stack spin direction.
        let theta = RIDE_THETA_DEG.to_radians();
        let lat = RIDE_R * theta.sin();
        let ax = RIDE_R * theta.cos();
        let e_top = core
            .spawn_particle(
                e_id,
                DVec3::new(lat, 1.3 + ax, 0.0),
                DVec3::new(0.0, 0.0, -crate::atom_core::COROT_V_MAX),
                DVec3::Y,
            )
            .unwrap();
        let e_bot = core
            .spawn_particle(
                e_id,
                DVec3::new(lat, -1.3 - ax, 0.0),
                DVec3::new(0.0, 0.0, -crate::atom_core::COROT_V_MAX),
                -DVec3::Y,
            )
            .unwrap();

        core.step_n(300_000);

        let d_top = core.pair_distance(top_proton, e_top);
        let d_bot = core.pair_distance(bottom_proton, e_bot);
        assert!(
            d_top < 2.0 && d_bot < 2.0,
            "helium electrons should stay bound: d_top={d_top} d_bot={d_bot}"
        );
        let com_drift = core.groups[gid].com.length();
        assert!(
            com_drift < 1.0,
            "helium nucleus should not self-propel: drift={com_drift}"
        );
        for p in &core.particles {
            assert!(p.position.is_finite(), "non-finite state");
        }
    }

    // ── Part 2: sim-driven nuclei validation harness (session-31 addendum
    // A11) ──────────────────────────────────────────────────────────────

    /// A11's uniform-boost sweep (7 boosts, 1x-12x) FAILED monotonically:
    /// adjacent alphas' facing protons rest only 1.15 apart
    /// (`ALPHA_PITCH − NUCLEON_PITCH = 3.75 − 2.6`) — deep inside the
    /// H₂-scale stream-cushion standoff (~3.5 for bare facing poles, see
    /// `Couplings::default` doc), and because a UNIFORM
    /// `intra_nucleus_boost` scales every attractive AND repulsive term by
    /// the same factor, it cannot move that force-balance point — drift
    /// *worsened monotonically* with boost (1233% at ×1 → 2056% at ×12).
    ///
    /// A13 (session-31 round 3) redesign: replace the uniform boost with
    /// **channeling attenuation** of `c_q`/`stream` (bb2.pdf: "attraction
    /// must always be explained as loss of repulsion") plus a
    /// **nuclear_ambient** shadow term (nuclear.pdf: "the charge field is
    /// both the initial pressure and the subsequent glue") — see
    /// `Couplings::channeling`/`Couplings::nuclear_ambient` docs and
    /// `compute_forces`.
    ///
    /// Session-31 round 4 (binding v2 hardening): the round-3 winner
    /// (channeling=0.9, ambient=2.0) held a QUIET nucleus but the user
    /// found two failure modes it missed — a spontaneous late "Jenga"
    /// collapse, and trivial destruction by a stray particle spawned
    /// nearby. Root cause: `channeling_factor`'s falloff hit a hard cliff
    /// at `r = 2·NUCLEON_PITCH` with no recapture basin beyond it (see
    /// `CHANNEL_TAIL`'s doc). `CHANNEL_TAIL` extends that taper, which
    /// changes the force balance enough to warrant a fresh grid: channeling
    /// {0.85, 0.9, 0.95} — 1.0 is EXCLUDED because with the longer tail it
    /// collapses catastrophically once ambient ≥ ~5 (full repulsion
    /// cancellation lets the never-attenuated disc-aware contact spring
    /// fire, drift in the millions of percent) — × nuclear_ambient
    /// {2, 5, 10, 20}, still 20k steps (~10s release). A SECOND pass
    /// dimension is added: of the quiet-sweep passers, which also survive
    /// `run_flyby_scenario` (an external close pass + a parked neighbor)?
    /// This runs the SAME full 15k+15k duration `nucleus_survives_flyby`
    /// uses, not a shorter proxy — a short 8k+8k window was tried first and
    /// found to be a false-positive trap: several combos "looked" robust
    /// at 8k+8k but the parked neutron's pull is a SUSTAINED perturbation,
    /// and those combos only started destabilizing between step ~16k and
    /// ~27k, well past an 8k+8k window (see this test's and `nucleus_
    /// survives_flyby`'s doc for the concrete numbers). The chosen defaults
    /// must clear BOTH bars — picked by the largest margin (lowest
    /// worst-case max-drift across quiet AND flyby), ties broken toward
    /// mid-range ambient, away from the collapse-adjacent high end the
    /// quiet table's NONFIN/eject rows expose. A final 250k-step QUIET run
    /// at the chosen defaults confirms (past the ~160k-step undamped-breathing horizon)
    /// the "Jenga" horizon (the user watched for minutes before the late
    /// collapse; the 20k-step sweep alone wouldn't necessarily have caught
    /// a failure that only shows up that late).
    #[test]
    fn alpha_stays_bound() {
        const STEPS: usize = 20_000;
        const SAMPLE_EVERY: usize = 50;
        let channelings = [0.85, 0.9, 0.95];
        let ambients = [2.0, 5.0, 10.0, 20.0];

        println!("\n-- quiet sweep (RigidAlpha carbon, {STEPS} steps) --");
        println!(
            "{:>5} {:>6} {:>10} {:>9} {:>9} {:>10}",
            "chan", "amb", "max_drift", "min_d", "max_d", "class"
        );

        // (channeling, ambient, max_rel_drift)
        let mut passing: Vec<(f64, f64, f64)> = Vec::new();
        let mut all_results: Vec<(f64, f64, f64)> = Vec::new();
        let mut best: Option<(f64, f64, f64)> = None;

        for &channeling in &channelings {
            for &ambient in &ambients {
                let mut core = standard_core();
                core.couplings.channeling = channeling;
                core.couplings.nuclear_ambient = ambient;
                let gid = core
                    .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("carbon preset");
                core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
                let n_alphas = core.groups[gid].alphas.len();
                assert_eq!(n_alphas, 3, "carbon should have 3 alphas");

                let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                let pairs: Vec<(usize, usize)> = (0..n_alphas)
                    .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                    .collect();
                let d0: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                let mut min_d = d0.clone();
                let mut max_d = d0.clone();

                let n_samples = STEPS / SAMPLE_EVERY;
                for _ in 0..n_samples {
                    core.step_n(SAMPLE_EVERY);
                    for (k, &(a, b)) in pairs.iter().enumerate() {
                        let d = (com(&core, a) - com(&core, b)).length();
                        min_d[k] = min_d[k].min(d);
                        max_d[k] = max_d[k].max(d);
                    }
                }

                let rel_drifts: Vec<f64> = (0..pairs.len())
                    .map(|k| {
                        let up = (max_d[k] - d0[k]) / d0[k];
                        let down = (min_d[k] - d0[k]) / d0[k];
                        if up.abs() >= down.abs() { up } else { down }
                    })
                    .collect();
                let max_rel_drift = rel_drifts.iter().fold(0.0f64, |acc, &d| acc.max(d.abs()));
                let any_eject = rel_drifts.iter().any(|&d| d > 0.30);
                let any_collapse = rel_drifts.iter().any(|&d| d < -0.30);
                let finite = core.particles.iter().all(|p| p.position.is_finite());

                let class = if !finite {
                    "NONFIN"
                } else if any_eject && any_collapse {
                    "shear"
                } else if any_eject {
                    "eject"
                } else if any_collapse {
                    "collapse"
                } else {
                    "stable"
                };

                println!(
                    "{channeling:>5.2} {ambient:>6.1} {:>9.1}% {:>9.3} {:>9.3} {:>10}",
                    max_rel_drift * 100.0,
                    min_d.iter().cloned().fold(f64::MAX, f64::min),
                    max_d.iter().cloned().fold(0.0f64, f64::max),
                    class,
                );

                assert!(
                    finite,
                    "channeling={channeling} ambient={ambient}: non-finite state"
                );

                all_results.push((channeling, ambient, max_rel_drift));
                if max_rel_drift <= 0.30 {
                    passing.push((channeling, ambient, max_rel_drift));
                }
                if best.map_or(true, |(_, _, bd)| max_rel_drift < bd) {
                    best = Some((channeling, ambient, max_rel_drift));
                }
            }
        }

        let (bc, ba, bd) = best.expect("sweep always has a best");

        if passing.is_empty() {
            // HONESTY CLAUSE (A13): nothing binds — print a failure-step
            // force breakdown for the best combo (nearest member pair
            // between the two adjacent alphas) showing WHICH force is
            // unbalanced, then fail. Do NOT loosen the ±30% to pass.
            println!(
                "\nno combo binds within +/-30% (best: channeling={bc} ambient={ba} drift={:.1}%)",
                bd * 100.0
            );
            let mut core = standard_core();
            core.couplings.channeling = bc;
            core.couplings.nuclear_ambient = ba;
            let gid = core
                .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let d0_01 = (com(&core, 0) - com(&core, 1)).length();

            let members_of =
                |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].members.clone();
            let (m0, m1) = (members_of(&core, 0), members_of(&core, 1));
            let closest_pair = |core: &AtomCore| -> (usize, usize) {
                let mut best_pair = (
                    core.groups[gid].members[m0[0]],
                    core.groups[gid].members[m1[0]],
                );
                let mut best_r = f64::MAX;
                for &ki in &m0 {
                    for &kj in &m1 {
                        let pi = core.groups[gid].members[ki];
                        let pj = core.groups[gid].members[kj];
                        let r = core.pair_distance(pi, pj);
                        if r < best_r {
                            best_r = r;
                            best_pair = (pi, pj);
                        }
                    }
                }
                best_pair
            };

            println!(
                "\n{:>6} {:>8} {:>8} {:>10} {:>10} {:>10} {:>10} {:>10} {:>10}",
                "step", "d01", "drift%", "r", "channel", "f_grav", "f_charge", "f_intake",
                "f_stream"
            );
            const FSAMPLE: usize = 200;
            for s in 0..(STEPS / FSAMPLE) {
                core.step_n(FSAMPLE);
                let d01 = (com(&core, 0) - com(&core, 1)).length();
                let drift = (d01 - d0_01) / d0_01;
                let (pi, pj) = closest_pair(&core);
                let fb = core.pair_force_breakdown(pi, pj);
                println!(
                    "{:>6} {:>8.3} {:>7.1}% {:>10.4} {:>10.4} {:>10.3} {:>10.3} {:>10.3} {:>10.3}",
                    s * FSAMPLE,
                    d01,
                    drift * 100.0,
                    fb[0],
                    fb[1],
                    fb[2],
                    fb[3],
                    fb[5],
                    fb[6],
                );
                if drift.abs() > 0.30 {
                    println!(
                        "  ^ drift exceeded +/-30% at step {} (contact term f_contact={:.3})",
                        s * FSAMPLE,
                        fb[7]
                    );
                    break;
                }
            }
            panic!(
                "no (channeling, nuclear_ambient) combo binds carbon within +/-30%; \
                 best combo channeling={bc} ambient={ba} drift={:.1}% — see tables above",
                bd * 100.0
            );
        }

        println!(
            "\n{} combo(s) bind quietly within +/-30%; best: channeling={bc} ambient={ba} drift={:.1}%",
            passing.len(),
            bd * 100.0
        );

        // Second pass dimension (work item 3): of the quiet-sweep passers,
        // which also survive `run_flyby_scenario` — an external close pass
        // + a parked neighbor?
        //
        // DEVIATION from the spec's "short version: 8k+8k steps": measured
        // empirically, the short window is NOT a reliable proxy. The
        // parked "swat" neutron's pull is a SUSTAINED perturbation, not a
        // transient one — several combos (e.g. channeling=0.85,
        // ambient=2.0) looked perfectly robust at 8k+8k (1.8% drift) but
        // then slowly destabilized and blew past the recapture basin
        // between step ~16k and ~27k of the full run (see
        // `nucleus_survives_flyby`'s FAILED run against that combo: 63.5%
        // final drift). A short-window "robust" verdict that the full
        // scenario immediately contradicts is worse than useless — it
        // launders a bad default past this test only to have `nucleus_
        // survives_flyby` fail it a moment later. So this uses the SAME
        // full 15k+15k duration as `nucleus_survives_flyby` itself (it's
        // still cheap: ~20s for the whole passing-combo set, well inside
        // budget) — the two tests must agree on what "robust" means.
        const FLYBY_STEPS: usize = 15_000;
        const FLYBY_SAMPLE: usize = 1_000;
        println!(
            "\n-- flyby/swat robustness ({FLYBY_STEPS}+{FLYBY_STEPS} steps, quiet-sweep passers only) --"
        );
        println!(
            "{:>5} {:>6} {:>10} {:>10} {:>8} {:>10}",
            "chan", "amb", "quiet_dr", "flyby_dr", "finite", "class"
        );
        // (channeling, ambient, quiet_drift, flyby_drift)
        let mut robust: Vec<(f64, f64, f64, f64)> = Vec::new();
        for &(channeling, ambient, quiet_drift) in &passing {
            let (flyby_drift, finite) =
                run_flyby_scenario(channeling, ambient, FLYBY_STEPS, FLYBY_STEPS, FLYBY_SAMPLE);
            let class = if !finite {
                "NONFIN"
            } else if flyby_drift > FLYBY_DRIFT_LIMIT {
                "eject"
            } else {
                "robust"
            };
            println!(
                "{channeling:>5.2} {ambient:>6.1} {:>9.1}% {:>9.1}% {finite:>8} {:>10}",
                quiet_drift * 100.0,
                flyby_drift * 100.0,
                class,
            );
            if finite && flyby_drift <= FLYBY_DRIFT_LIMIT {
                robust.push((channeling, ambient, quiet_drift, flyby_drift));
            }
        }

        // LONG-HORIZON VETO (session-32 Phase B lesson): combos that won
        // this 20k-sweep + 15k-flyby chooser but FAILED the 8-seed ×
        // 1M-step gate (`report_long_horizon_drift`) — the short-horizon
        // metrics have a blind spot past the 250k confirm, and defaults
        // must never ship on them alone. Each entry carries its gate
        // evidence; re-run the gate to challenge an entry.
        //
        // SESSION-34 RE-DERIVATION: the flow-network NO-STARVE fix
        // (charge_flow.rs pass 2) changed every live flow amplitude and
        // therefore the tension/align forces, so ALL pre-v2 gate
        // evidence was retired and the five contenders were re-gated
        // under the new network. Results (8 seeds × 1M each):
        //   (0.85,  5.0) 8/8 all-cold [3.34/6.68, KE 5e-4] → UN-VETOED,
        //                and as best surviving margin (13.7%/12.5%) it
        //                is the shipped default again (as in session 32).
        //   (0.85, 10.0) 8/8 all-cold [3.24/6.47] — clean runner-up,
        //                not vetoed (worse margins, never wins).
        //
        // Standing vetoes (network-v2 evidence unless noted):
        // - (0.9, 5.0): the naive post-fix chooser winner — seeds k=0/k=4
        //   collapse at 800k/350k.
        // - (0.9, 2.0): seeds collapse 2/8 (also failed 4/8 in the
        //   session-32 starved-network era — consistently bad).
        // - (0.85, 2.0): no collapse by 1M, but 3/8 seeds heat
        //   MONOTONICALLY to KE ≈ 128 (vs 4e-4 cold) with growing tilt —
        //   the pre-collapse whirl signature; vetoed for failing to
        //   reach a stable attractor. Low ambient keeps losing the long
        //   game in every era.
        // - (0.95, 2.0): STARVED-NETWORK evidence only (4/8 collapse,
        //   session 33) — not re-gated under v2 because it cannot win
        //   this chooser anyway (15.0% quiet loses to 13.7%); re-gate it
        //   before ever un-vetoing.
        const LONG_HORIZON_VETO: &[(f64, f64)] =
            &[(0.9, 5.0), (0.9, 2.0), (0.85, 2.0), (0.95, 2.0)];

        // Choose defaults = the combo with the largest margin against BOTH
        // failure modes (lowest worst-case max-drift across quiet AND
        // flyby), ties broken toward mid-range ambient — away from the
        // collapse-adjacent high end the quiet table's NONFIN/eject rows
        // expose.
        let mid_ambient = ambients.iter().sum::<f64>() / ambients.len() as f64;
        let mut chosen: Option<(f64, f64, f64, f64)> = None; // (chan, amb, quiet_drift, flyby_drift)
        for &(c, a, qd, fd) in &robust {
            if LONG_HORIZON_VETO
                .iter()
                .any(|&(vc, va)| (c - vc).abs() < 1e-9 && (a - va).abs() < 1e-9)
            {
                println!(
                    "  (skipping channeling={c} ambient={a}: long-horizon veto — \
                     failed the 8-seed 1M gate, see LONG_HORIZON_VETO)"
                );
                continue;
            }
            let margin = qd.max(fd);
            let take = match chosen {
                None => true,
                Some((_, ca, cqd, cfd)) => {
                    let cur_margin = cqd.max(cfd);
                    if (margin - cur_margin).abs() < 1e-9 {
                        (a - mid_ambient).abs() < (ca - mid_ambient).abs()
                    } else {
                        margin < cur_margin
                    }
                }
            };
            if take {
                chosen = Some((c, a, qd, fd));
            }
        }

        let (def_chan, def_amb) = match chosen {
            Some((c, a, qd, fd)) => {
                println!(
                    "\nchosen defaults: channeling={c} nuclear_ambient={a} \
                     (quiet drift {:.1}%, flyby drift {:.1}%) — {} of {} quiet-passers also robust",
                    qd * 100.0,
                    fd * 100.0,
                    robust.len(),
                    passing.len()
                );
                (c, a)
            }
            None => {
                // HONESTY CLAUSE (work item 3 / CONSTRAINTS): quiet binding
                // works but NOTHING survives the flyby/swat perturbation.
                // Report it plainly and fall back to the best quiet combo
                // for `Couplings::default` (it still stops the "Jenga"
                // quiet collapse); `nucleus_survives_flyby` documents the
                // unresolved flyby fragility and must be `#[ignore]`d with
                // that analysis if it fails at these defaults — see its
                // doc comment.
                println!(
                    "\nno quiet-passing combo survives the flyby/swat perturbation; \
                     falling back to the best quiet combo channeling={bc} ambient={ba} \
                     ({:.1}% quiet drift) for Couplings::default — see \
                     nucleus_survives_flyby for the documented flyby-robustness gap",
                    bd * 100.0
                );
                (bc, ba)
            }
        };

        // Assert the shipped Couplings::default MATCHES this sweep's
        // chosen winner — any future default change must re-earn this by
        // re-running (and, if the winner moves, updating) the sweep.
        let def = crate::atom_core::Couplings::default();
        assert!(
            (def.channeling - def_chan).abs() < 1e-9 && (def.nuclear_ambient - def_amb).abs() < 1e-9,
            "Couplings::default (channeling={}, nuclear_ambient={}) does not match this sweep's \
             chosen winner (channeling={def_chan}, nuclear_ambient={def_amb}) — update \
             Couplings::default to match, or re-derive the winner if the sweep changed",
            def.channeling,
            def.nuclear_ambient
        );
        let def_result = all_results
            .iter()
            .find(|&&(c, a, _)| {
                (c - def.channeling).abs() < 1e-12 && (a - def.nuclear_ambient).abs() < 1e-12
            })
            .unwrap_or_else(|| {
                panic!(
                    "sweep grid must include the Couplings defaults \
                     (channeling={}, nuclear_ambient={})",
                    def.channeling, def.nuclear_ambient
                )
            });
        assert!(
            def_result.2 <= 0.30,
            "Couplings::default (channeling={}, nuclear_ambient={}) no longer binds carbon quietly: \
             max drift {:.1}% > 30%",
            def.channeling,
            def.nuclear_ambient,
            def_result.2 * 100.0
        );

        // Final confirmation (work item 3): a 250k-step QUIET run at the
        // chosen defaults — the "Jenga" horizon the user watched for
        // minutes before the spontaneous late collapse. No external
        // perturbation; this purely re-checks binding over a much longer
        // horizon than the 20k-step sweep, since the user's report was a
        // LATE failure a short sweep wouldn't necessarily catch.
        // 250k: past the ~160k-step pumping horizon where the UNDAMPED
        // breathing mode ejected an alpha (session-31 round 5 — the 100k
        // confirm sat comfortably below the horizon and false-passed;
        // ALPHA_TRANS_RELAX is the fix, this is the regression net).
        const CONFIRM_STEPS: usize = 250_000;
        const CONFIRM_SAMPLE: usize = 200;
        let mut core = standard_core();
        core.couplings.channeling = def.channeling;
        core.couplings.nuclear_ambient = def.nuclear_ambient;
        let gid = core
            .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
        let n_alphas = core.groups[gid].alphas.len();
        let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
        let pairs: Vec<(usize, usize)> = (0..n_alphas)
            .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
            .collect();
        let d0: Vec<f64> = pairs
            .iter()
            .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
            .collect();
        let mut confirm_max_drift = 0.0f64;
        let mut confirm_finite = true;
        for _ in 0..(CONFIRM_STEPS / CONFIRM_SAMPLE) {
            core.step_n(CONFIRM_SAMPLE);
            for (k, &(a, b)) in pairs.iter().enumerate() {
                let d = (com(&core, a) - com(&core, b)).length();
                confirm_max_drift = confirm_max_drift.max(((d - d0[k]) / d0[k]).abs());
            }
            confirm_finite &= core.particles.iter().all(|p| p.position.is_finite());
        }
        println!(
            "\n250k-step quiet confirmation at defaults (channeling={} ambient={}): \
             max_drift={:.1}% finite={}",
            def.channeling,
            def.nuclear_ambient,
            confirm_max_drift * 100.0,
            confirm_finite
        );
        assert!(
            confirm_finite,
            "non-finite state over the 250k-step ('Jenga horizon') confirmation run"
        );
        assert!(
            confirm_max_drift <= 0.30,
            "carbon nucleus did not survive the 250k-step ('Jenga horizon') confirmation run \
             at chosen defaults: max drift {:.1}% > 30%",
            confirm_max_drift * 100.0
        );
    }

    // ── Multi-seed transient stability (session-31 round 5, work item 2)
    // ────────────────────────────────────────────────────────────────

    /// Burn `k` throwaway "alpha" preset spawns (each consumes one
    /// `roll_phase` draw + one `carousel_phase` draw from `core.rng` —
    /// see `spawn_preset`) so the REAL spawn that follows lands on a
    /// different random roll/carousel phase set. `clear_particles` does
    /// NOT reset the rng, so this is the cheapest way to walk through
    /// distinct seed states without threading a seed parameter through
    /// `AtomCore::new`.
    fn burn_seed_offset(core: &mut AtomCore, k: usize) {
        for _ in 0..k {
            core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("throwaway alpha for seed offset");
            core.clear_particles();
        }
    }

    /// Re-drive one failing seed and print: the drift timeline, each
    /// alpha's axis tilt (`orientation * Y` dotted with world `Y` — 1.0 =
    /// still pole-aligned with the stack axis, drifting below that means
    /// the torque-driven tilt mechanism, not just translation, is in
    /// play), the channeling factor for the worst-drifting pair's closest
    /// member pair, and the full `pair_force_breakdown` at the sample the
    /// drift limit is crossed (A13 honesty-clause style diagnostics, same
    /// pattern as `print_flyby_diagnostics`).
    fn print_transient_failure_diagnostics(
        k: usize,
        steps: usize,
        sample_every: usize,
        drift_limit: f64,
    ) {
        let mut core = standard_core();
        burn_seed_offset(&mut core, k);
        let gid = core
            .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);

        let n_alphas = core.groups[gid].alphas.len();
        let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
        let pairs: Vec<(usize, usize)> = (0..n_alphas)
            .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
            .collect();
        let d0: Vec<f64> = pairs
            .iter()
            .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
            .collect();
        let members_of =
            |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].members.clone();
        let closest_pair = |core: &AtomCore, ai: usize, aj: usize| -> (usize, usize) {
            let (mi, mj) = (members_of(core, ai), members_of(core, aj));
            let mut best_pair = (
                core.groups[gid].members[mi[0]],
                core.groups[gid].members[mj[0]],
            );
            let mut best_r = f64::MAX;
            for &ki in &mi {
                for &kj in &mj {
                    let pi = core.groups[gid].members[ki];
                    let pj = core.groups[gid].members[kj];
                    let r = core.pair_distance(pi, pj);
                    if r < best_r {
                        best_r = r;
                        best_pair = (pi, pj);
                    }
                }
            }
            best_pair
        };

        println!(
            "\n-- transient failure diagnostics (seed k={k}) --\n\
             {:>8} {:>8} {:>28} {:>10}",
            "step", "drift%", "axis.Y per alpha", "channel"
        );
        for s in 0..(steps / sample_every) {
            core.step_n(sample_every);
            let mut worst_idx = 0usize;
            let mut worst_drift = 0.0f64;
            for (idx, &(a, b)) in pairs.iter().enumerate() {
                let d = (com(&core, a) - com(&core, b)).length();
                let drift = ((d - d0[idx]) / d0[idx]).abs();
                if drift > worst_drift {
                    worst_drift = drift;
                    worst_idx = idx;
                }
            }
            let (wa, wb) = pairs[worst_idx];
            let (pi, pj) = closest_pair(&core, wa, wb);
            let fb = core.pair_force_breakdown(pi, pj);
            let tilts: Vec<f64> = (0..n_alphas)
                .map(|ai| (core.groups[gid].alphas[ai].orientation * DVec3::Y).dot(DVec3::Y))
                .collect();
            println!(
                "{:>8} {:>7.1}% {:>28} {:>10.4}",
                (s + 1) * sample_every,
                worst_drift * 100.0,
                format!("{tilts:.3?}"),
                fb[1],
            );
            let finite = core.particles.iter().all(|p| p.position.is_finite());
            if !finite || worst_drift > drift_limit {
                println!(
                    "  ^ seed k={k} crossed +/-{:.0}% at step {} (worst pair {wa}-{wb}, \
                     r={:.4} channel={:.4} f_grav={:.3} f_charge={:.3} f_ambient={:.3} \
                     f_intake={:.3} f_stream={:.3} f_contact={:.3}, finite={finite})",
                    drift_limit * 100.0,
                    (s + 1) * sample_every,
                    fb[0], fb[1], fb[2], fb[3], fb[4], fb[5], fb[6], fb[7],
                );
                break;
            }
        }
    }

    /// Root-cause regression test (session-31 round 5 — "which alpha
    /// detaches, or whether any does, depends on the random initial roll
    /// phases"): the ~160k-step "Jenga" transient walks the equilibrium
    /// spacing 3.75 → 3.29 → ~3.71 as the SEEDED alpha spin decays
    /// (`ALPHA_SPIN_RELAX`, τ=2s); which pair crosses the channeling basin
    /// during that migration is a matter of seed-phase luck, so a single
    /// fixed-seed run (as `alpha_stays_bound`'s 250k confirm uses) cannot
    /// certify the fix. `alpha_kinematic_state`'s rest-entry seeding fix
    /// (dropping the display roll/carousel rate from the seeded
    /// velocity/angular_velocity — see its doc comment) removes the
    /// fictional injected spin that drove the migration; this drives 8
    /// distinct roll/carousel phase draws (`burn_seed_offset`, k=0..8)
    /// through 150k steps each (the dissolve probe's migration window
    /// closes by ~100k, so 150k covers the whole transient with margin
    /// well short of the 250k-step "Jenga" horizon, keeping 8 seeds
    /// affordable) and requires every inter-alpha spacing to stay within
    /// +/-35% for EVERY seed — not just the one `alpha_stays_bound`
    /// happens to draw.
    #[test]
    fn alpha_transient_survives_all_seeds() {
        const STEPS: usize = 150_000;
        const SAMPLE_EVERY: usize = 2_000;
        const DRIFT_LIMIT: f64 = 0.35;

        println!("\n-- multi-seed transient stability (RigidAlpha carbon, {STEPS} steps/seed) --");
        println!("{:>4} {:>10} {:>8}", "seed", "max_drift", "finite");

        let mut failures: Vec<(usize, f64)> = Vec::new();
        for k in 0..8usize {
            let mut core = standard_core();
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            assert_eq!(n_alphas, 3, "carbon should have 3 alphas");

            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();

            let mut max_drift = 0.0f64;
            let mut finite = true;
            'outer: for _ in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                for (idx, &(a, b)) in pairs.iter().enumerate() {
                    let d = (com(&core, a) - com(&core, b)).length();
                    max_drift = max_drift.max(((d - d0[idx]) / d0[idx]).abs());
                }
                finite = core.particles.iter().all(|p| p.position.is_finite());
                if !finite || max_drift > DRIFT_LIMIT {
                    break 'outer;
                }
            }

            println!("{k:>4} {:>9.1}% {finite:>8}", max_drift * 100.0);
            if !finite || max_drift > DRIFT_LIMIT {
                failures.push((k, max_drift));
            }
        }

        if !failures.is_empty() {
            println!(
                "\n{} of 8 seed(s) failed: {:?}",
                failures.len(),
                failures
                    .iter()
                    .map(|&(k, d)| format!("k={k} drift={:.1}%", d * 100.0))
                    .collect::<Vec<_>>()
            );
            // HONESTY CLAUSE (A13): print full diagnostics for every
            // failing seed rather than just the first, since different
            // seeds may fail via different pairs/mechanisms.
            for &(k, _) in &failures {
                print_transient_failure_diagnostics(k, STEPS, SAMPLE_EVERY, DRIFT_LIMIT);
            }
        }

        assert!(
            failures.is_empty(),
            "alpha transient did not survive all 8 seed phases within +/-{:.0}%: {:?} — \
             see per-seed diagnostics above",
            DRIFT_LIMIT * 100.0,
            failures
                .iter()
                .map(|&(k, d)| format!("k={k} drift={:.1}%", d * 100.0))
                .collect::<Vec<_>>()
        );
    }

    /// Session-32 long-horizon drift report. The user's live carbon
    /// collapse at ~300-400 sim-time units = 600k-800k steps (dt=0.0005)
    /// sits far beyond every green harness (150k transient seeds, 250k
    /// quiet confirm) — the round-6 fixes stretched the horizon ~4x but
    /// something still drifts unboundedly. This drives 8 seed phases
    /// through 1M steps each and logs the full state signature at every
    /// sample (spacings, kinetic energy, alpha spin/relative-velocity,
    /// worst axis tilt, channeling factor of the tightest member pair) so
    /// the slow mechanism shows itself BEFORE the spacing blows up.
    /// Report only — no assertions. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_long_horizon_drift --nocapture`
    /// Optional env overrides for the ignored long-run reports, so a
    /// CANDIDATE default can be driven through the 1M gate without
    /// editing `Couplings::default` before it has earned the change
    /// (session-32 LONG_HORIZON_VETO discipline). Recognized:
    /// CFM_CHAN, CFM_AMB, CFM_TENSION, CFM_ALIGN, CFM_SUCTION,
    /// CFM_EMITSCALE, CFM_GAP.
    fn apply_env_overrides(core: &mut AtomCore) {
        let get = |k: &str| std::env::var(k).ok().and_then(|s| s.parse::<f64>().ok());
        if let Some(v) = get("CFM_CHAN") {
            core.couplings.channeling = v;
        }
        if let Some(v) = get("CFM_AMB") {
            core.couplings.nuclear_ambient = v;
        }
        if let Some(v) = get("CFM_TENSION") {
            core.flow_tension = v;
        }
        if let Some(v) = get("CFM_ALIGN") {
            core.flow_align = v;
        }
        if let Some(v) = get("CFM_SUCTION") {
            core.flow_suction = v;
        }
        if let Some(v) = get("CFM_EMITSCALE") {
            core.flow_emit_scale = v;
        }
        if let Some(v) = get("CFM_GAP") {
            core.plug_pair_gap = v;
        }
        if let Some(v) = get("CFM_FLOWEVERY") {
            core.flow_solve_every = (v as usize).max(1);
        }
        if let Some(v) = get("CFM_PLUGLOCK") {
            core.plug_orient_lock = v > 0.5;
        }
        if let Some(v) = get("CFM_AMBALIGN") {
            core.ambient_align = v;
        }
        if let Some(v) = get("CFM_BONDSUCTION") {
            core.bond_suction = v;
        }
        if let Some(v) = get("CFM_BONDTENSION") {
            core.bond_tension = v;
        }
        // Directional ambient flow along −Y (charge travels downward), so a
        // nucleus's +Y pole faces upstream and becomes the intake socket.
        // Magnitude = bias strength (see solve_charge_flow). 0 = isotropic.
        if let Some(v) = get("CFM_AMBFLOW") {
            core.ambient_charge_dir = DVec3::new(0.0, -v, 0.0);
        }
    }

    #[test]
    #[ignore]
    fn report_long_horizon_drift() {
        const STEPS: usize = 1_000_000;
        const SAMPLE_EVERY: usize = 25_000;

        for k in 0..8usize {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            if k == 0 {
                println!(
                    "couplings: channeling={} nuclear_ambient={} flow_tension={} flow_align={}",
                    core.couplings.channeling,
                    core.couplings.nuclear_ambient,
                    core.flow_tension,
                    core.flow_align,
                );
            }
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();

            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();
            let members_of =
                |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].members.clone();
            let closest_pair = |core: &AtomCore, ai: usize, aj: usize| -> (usize, usize) {
                let (mi, mj) = (members_of(core, ai), members_of(core, aj));
                let mut best_pair = (
                    core.groups[gid].members[mi[0]],
                    core.groups[gid].members[mj[0]],
                );
                let mut best_r = f64::MAX;
                for &ki in &mi {
                    for &kj in &mj {
                        let pi = core.groups[gid].members[ki];
                        let pj = core.groups[gid].members[kj];
                        let r = core.pair_distance(pi, pj);
                        if r < best_r {
                            best_r = r;
                            best_pair = (pi, pj);
                        }
                    }
                }
                best_pair
            };

            println!(
                "\n== seed k={k} (d0={d0:.4?}) ==\n\
                 {:>8} {:>7} {:>22} {:>10} {:>8} {:>8} {:>7} {:>8}",
                "step", "worst%", "spacings", "KE", "mean|w|", "mean|dv|", "tilt", "channel"
            );
            for s in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                let mut worst_idx = 0usize;
                let mut worst_drift = 0.0f64;
                let spacings: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                for (idx, &d) in spacings.iter().enumerate() {
                    let drift = ((d - d0[idx]) / d0[idx]).abs();
                    if drift > worst_drift {
                        worst_drift = drift;
                        worst_idx = idx;
                    }
                }
                let v_mean = (0..n_alphas)
                    .map(|ai| core.groups[gid].alphas[ai].velocity)
                    .fold(DVec3::ZERO, |acc, v| acc + v)
                    / n_alphas as f64;
                let mean_w = (0..n_alphas)
                    .map(|ai| core.groups[gid].alphas[ai].angular_velocity.length())
                    .sum::<f64>()
                    / n_alphas as f64;
                let mean_dv = (0..n_alphas)
                    .map(|ai| (core.groups[gid].alphas[ai].velocity - v_mean).length())
                    .sum::<f64>()
                    / n_alphas as f64;
                let min_tilt = (0..n_alphas)
                    .map(|ai| (core.groups[gid].alphas[ai].orientation * DVec3::Y).dot(DVec3::Y))
                    .fold(f64::MAX, f64::min);
                let (wa, wb) = pairs[worst_idx];
                let (pi, pj) = closest_pair(&core, wa, wb);
                let fb = core.pair_force_breakdown(pi, pj);
                println!(
                    "{:>8} {:>6.1}% {:>22} {:>10.4e} {:>8.5} {:>8.5} {:>7.4} {:>8.4}",
                    (s + 1) * SAMPLE_EVERY,
                    worst_drift * 100.0,
                    format!("{spacings:.3?}"),
                    core.total_kinetic_energy(),
                    mean_w,
                    mean_dv,
                    min_tilt,
                    fb[1],
                );
                let finite = core.particles.iter().all(|p| p.position.is_finite());
                if !finite || worst_drift > 1.0 {
                    println!(
                        "  ^ seed k={k} COLLAPSED at step {} (worst pair {wa}-{wb}, \
                         r={:.4} channel={:.4} f_grav={:.3} f_charge={:.3} \
                         f_ambient={:.3} f_intake={:.3} f_stream={:.3} f_contact={:.3}, \
                         finite={finite})",
                        (s + 1) * SAMPLE_EVERY,
                        fb[0], fb[1], fb[2], fb[3], fb[4], fb[5], fb[6], fb[7],
                    );
                    break;
                }
            }
        }
    }

    /// Session-32 onset zoom. `report_long_horizon_drift` showed every
    /// seed detonates at 525k-700k steps with an ABRUPT signature: quiet
    /// at ~1e-4 relative velocity for hundreds of k-steps, then KE x10 in
    /// one 25k sample, then a ~75k-step death spiral. Channeling stays at
    /// its maximum until AFTER the KE spike, so channeling erosion is the
    /// consequence, not the trigger. This re-drives seed 0 (deterministic),
    /// fast-forwards through the quiet phase, then samples every 500 steps
    /// through the detonation window logging: per-alpha roll azimuth (and
    /// relative azimuth between adjacent alphas — the posts sweep past
    /// each other as rolls decorrelate), axial vs tumble spin split,
    /// closest member pair between adjacent alphas with its profile kinds
    /// and force breakdown. Report only. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_onset_zoom --nocapture`
    #[test]
    #[ignore]
    fn report_onset_zoom() {
        const QUIET_STEPS: usize = 540_000;
        const ZOOM_SAMPLES: usize = 220;
        const SAMPLE_EVERY: usize = 500;

        let mut core = standard_core();
        burn_seed_offset(&mut core, 0);
        let gid = core
            .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
        let n_alphas = core.groups[gid].alphas.len();

        // Roll azimuth: where the alpha's body X axis points in the world
        // XZ plane (valid while the alpha axis stays ~Y, which holds until
        // the detonation is already underway).
        let azimuth = |core: &AtomCore, ai: usize| -> f64 {
            let x = core.groups[gid].alphas[ai].orientation * DVec3::X;
            x.z.atan2(x.x).to_degrees()
        };
        let members_of =
            |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].members.clone();
        // Closest member pair between adjacent alphas ai/aj: particle ids,
        // distance, and profile ids (0=proton 1=neutron).
        let closest = |core: &AtomCore, ai: usize, aj: usize| -> (usize, usize, f64) {
            let (mi, mj) = (members_of(core, ai), members_of(core, aj));
            let mut best = (0usize, 0usize, f64::MAX);
            for &ki in &mi {
                for &kj in &mj {
                    let pi = core.groups[gid].members[ki];
                    let pj = core.groups[gid].members[kj];
                    let r = core.pair_distance(pi, pj);
                    if r < best.2 {
                        best = (pi, pj, r);
                    }
                }
            }
            best
        };

        core.step_n(QUIET_STEPS);
        println!(
            "\n== onset zoom, seed k=0, from step {QUIET_STEPS} ==\n\
             {:>7} {:>10} {:>8} {:>8} {:>7} {:>7} {:>7} {:>17} {:>17}",
            "step", "KE", "rel01", "rel12", "w_ax", "w_tum", "tilt", "close01", "close12"
        );
        for s in 0..ZOOM_SAMPLES {
            core.step_n(SAMPLE_EVERY);
            let step = QUIET_STEPS + (s + 1) * SAMPLE_EVERY;
            let az: Vec<f64> = (0..n_alphas).map(|ai| azimuth(&core, ai)).collect();
            let rel01 = (az[1] - az[0]).rem_euclid(360.0);
            let rel12 = (az[2] - az[1]).rem_euclid(360.0);
            let (mut w_ax, mut w_tum) = (0.0f64, 0.0f64);
            let mut min_tilt = f64::MAX;
            for ai in 0..n_alphas {
                let a = &core.groups[gid].alphas[ai];
                let axis = a.orientation * DVec3::Y;
                let ax = a.angular_velocity.dot(axis);
                w_ax += ax.abs();
                w_tum += (a.angular_velocity - axis * ax).length();
                min_tilt = min_tilt.min(axis.dot(DVec3::Y));
            }
            w_ax /= n_alphas as f64;
            w_tum /= n_alphas as f64;
            let fmt_close = |core: &AtomCore, ai: usize, aj: usize| -> String {
                let (pi, pj, r) = closest(core, ai, aj);
                let (qi, qj) = (
                    core.particles[pi].profile_id,
                    core.particles[pj].profile_id,
                );
                let kind = |q: usize| if q == 1 { "n" } else { "p" };
                let fb = core.pair_force_breakdown(pi, pj);
                // net radial push-pull: charge+stream out, ambient+intake+grav in
                let net = fb[3] + fb[6] - fb[2] - fb[4] - fb[5];
                format!("{}{}{:5.2} net{:+.3}", kind(qi), kind(qj), r, net)
            };
            println!(
                "{:>7} {:>10.4e} {:>7.1} {:>7.1} {:>7.4} {:>7.4} {:>7.4} {:>17} {:>17}",
                step,
                core.total_kinetic_energy(),
                rel01,
                rel12,
                w_ax,
                w_tum,
                min_tilt,
                fmt_close(&core, 0, 1),
                fmt_close(&core, 1, 2),
            );
        }
    }

    /// Session-32 buckling-mode probe. `report_onset_zoom` showed the
    /// RigidAlpha carbon stack sits at an unstable equilibrium: tumble
    /// angular velocity grows EXPONENTIALLY from the numerical noise floor
    /// (e-fold ~9k steps) at frozen geometry — a compressed-column
    /// buckling saddle that velocity damping can slow but never stabilize.
    /// This applies a controlled 1e-4 perturbation to the MIDDLE alpha
    /// (pure tilt about X / pure lateral shear along X / control) and logs
    /// the growth curve plus the mode shape (per-alpha axis tilt and com
    /// lateral offset), so the fix can target the actual unstable
    /// direction. Report only. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_buckling_mode --nocapture`
    #[test]
    #[ignore]
    fn report_buckling_mode() {
        const SETTLE: usize = 20_000;
        const STEPS: usize = 160_000;
        const SAMPLE_EVERY: usize = 4_000;
        const EPS: f64 = 1e-4;

        for scenario in ["control", "tilt-mid", "shear-mid"] {
            let mut core = standard_core();
            let gid = core
                .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            core.step_n(SETTLE);

            match scenario {
                "tilt-mid" => {
                    let a = &mut core.groups[gid].alphas[1];
                    a.orientation = (glam::DQuat::from_axis_angle(DVec3::X, EPS)
                        * a.orientation)
                        .normalize();
                }
                "shear-mid" => {
                    core.groups[gid].alphas[1].com += DVec3::new(EPS, 0.0, 0.0);
                }
                _ => {}
            }

            let com0: Vec<DVec3> = (0..n_alphas)
                .map(|ai| core.groups[gid].alphas[ai].com)
                .collect();
            println!(
                "\n== {scenario} (eps={EPS}) ==\n\
                 {:>7} {:>10} {:>9} {:>9} {:>26} {:>26}",
                "step", "KE", "w_tum", "dv", "axis.x per alpha", "com.x offset per alpha"
            );
            for s in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                let v_mean = (0..n_alphas)
                    .map(|ai| core.groups[gid].alphas[ai].velocity)
                    .fold(DVec3::ZERO, |acc, v| acc + v)
                    / n_alphas as f64;
                let mut w_tum = 0.0f64;
                let mut dv = 0.0f64;
                let mut ax_x = Vec::new();
                let mut off_x = Vec::new();
                for ai in 0..n_alphas {
                    let a = &core.groups[gid].alphas[ai];
                    let axis = a.orientation * DVec3::Y;
                    let w_axial = axis * a.angular_velocity.dot(axis);
                    w_tum += (a.angular_velocity - w_axial).length();
                    dv += (a.velocity - v_mean).length();
                    ax_x.push(axis.x);
                    off_x.push(a.com.x - com0[ai].x);
                }
                w_tum /= n_alphas as f64;
                dv /= n_alphas as f64;
                println!(
                    "{:>7} {:>10.4e} {:>9.2e} {:>9.2e} {:>26} {:>26}",
                    SETTLE + (s + 1) * SAMPLE_EVERY,
                    core.total_kinetic_energy(),
                    w_tum,
                    dv,
                    format!("{ax_x:+.5?}"),
                    format!("{off_x:+.5?}"),
                );
            }
        }
    }

    /// Session-32 stability-vs-time probe. `report_buckling_mode` showed
    /// tilt/shear kicks DECAY at t~180k, yet `report_onset_zoom` showed
    /// noise-floor tumble growing exponentially from t~598k — the
    /// equilibrium starts stable and slowly drifts across a bifurcation.
    /// Two creeping candidates from the drift report: alpha roll spin
    /// (0.0387->0.0406, rotor whirl threshold) and stack compression
    /// (spacings tighten ~0.5%, Euler buckling load). This kicks the
    /// middle alpha (1e-4 tilt) at a series of checkpoints along the seed-0
    /// trajectory and logs whether each kick decays or grows, plus the
    /// slow-state (spin, spacings, net compression on the closest pair) at
    /// each checkpoint. Report only. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_stability_vs_time --nocapture`
    #[test]
    #[ignore]
    fn report_stability_vs_time() {
        const CHECKPOINTS: [usize; 5] = [100_000, 300_000, 450_000, 530_000, 570_000];
        const OBSERVE: usize = 36_000;
        const SAMPLE_EVERY: usize = 3_000;
        const EPS: f64 = 1e-4;

        let mut core = standard_core();
        let gid = core
            .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
        let n_alphas = core.groups[gid].alphas.len();

        let mut now = 0usize;
        for &cp in CHECKPOINTS.iter() {
            core.step_n(cp - now);
            now = cp;

            // Slow-state snapshot at the checkpoint.
            let w_ax_mean = (0..n_alphas)
                .map(|ai| {
                    let a = &core.groups[gid].alphas[ai];
                    a.angular_velocity.dot(a.orientation * DVec3::Y).abs()
                })
                .sum::<f64>()
                / n_alphas as f64;
            let d01 = (core.groups[gid].alphas[0].com - core.groups[gid].alphas[1].com)
                .length();
            let d12 = (core.groups[gid].alphas[1].com - core.groups[gid].alphas[2].com)
                .length();
            println!(
                "\n== checkpoint {cp}: w_ax={w_ax_mean:.5} d01={d01:.4} d12={d12:.4} \
                 KE={:.4e} ==",
                core.total_kinetic_energy()
            );

            // Kick and observe. NOTE: the kick itself perturbs the
            // trajectory that later checkpoints ride on — acceptable, the
            // kicks are tiny and decayed kicks leave ~nothing behind.
            {
                let a = &mut core.groups[gid].alphas[1];
                a.orientation = (glam::DQuat::from_axis_angle(DVec3::X, EPS)
                    * a.orientation)
                    .normalize();
            }
            println!("{:>9} {:>9} {:>9}", "step", "w_tum", "dv");
            for s in 0..(OBSERVE / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                now += SAMPLE_EVERY;
                let v_mean = (0..n_alphas)
                    .map(|ai| core.groups[gid].alphas[ai].velocity)
                    .fold(DVec3::ZERO, |acc, v| acc + v)
                    / n_alphas as f64;
                let mut w_tum = 0.0f64;
                let mut dv = 0.0f64;
                for ai in 0..n_alphas {
                    let a = &core.groups[gid].alphas[ai];
                    let axis = a.orientation * DVec3::Y;
                    let w_axial = axis * a.angular_velocity.dot(axis);
                    w_tum += (a.angular_velocity - w_axial).length();
                    dv += (a.velocity - v_mean).length();
                }
                println!(
                    "{:>9} {:>9.2e} {:>9.2e}",
                    cp + (s + 1) * SAMPLE_EVERY,
                    w_tum / n_alphas as f64,
                    dv / n_alphas as f64,
                );
            }
        }
    }

    /// Session-32 neutron-post anatomy experiment (start-here item 1a,
    /// user design question): compare `PostAnatomy::Axial` (session-31
    /// modeling choice) vs `PostAnatomy::Radial` (deut.pdf charge-channel/
    /// self-balancing-regulator reading) across the full stability
    /// battery. Report only — adoption of Radial requires re-earning the
    /// sweep tables and updating the default. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_post_anatomy --nocapture`
    #[test]
    #[ignore]
    fn report_post_anatomy() {
        use crate::atom_core::PostAnatomy;

        for anatomy in [PostAnatomy::Axial, PostAnatomy::Radial] {
            println!("\n===== {anatomy:?} =====");

            // (a) FreeNucleon lone alpha: intra-alpha cohesion with no
            // rigid constraint (nucleon_balance's scenario at the
            // concluded boost=1.0).
            {
                let mut core = standard_core();
                apply_env_overrides(&mut core);
                core.post_anatomy = anatomy;
                let gid = core
                    .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("alpha preset");
                core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::FreeNucleon);
                let members = core.groups[gid].members.clone();
                let d0 = pair_dists(&core, &members);
                core.step_n(20_000);
                let d1 = pair_dists(&core, &members);
                let max_rel = d0
                    .iter()
                    .zip(&d1)
                    .map(|(a, b)| (a - b).abs() / a.max(1e-9))
                    .fold(0.0f64, f64::max);
                let finite = core.particles.iter().all(|p| p.position.is_finite());
                println!(
                    "FreeNucleon lone alpha, 20k steps: max_pair_drift={:.1}% \
                     KE={:.4e} finite={finite}",
                    max_rel * 100.0,
                    core.total_kinetic_energy(),
                );
            }

            // (b) RigidAlpha carbon quiet run, 100k steps.
            {
                let mut core = standard_core();
                apply_env_overrides(&mut core);
                core.post_anatomy = anatomy;
                let gid = core
                    .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("carbon preset");
                core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
                let n_alphas = core.groups[gid].alphas.len();
                let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                let pairs: Vec<(usize, usize)> = (0..n_alphas)
                    .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                    .collect();
                let d0: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                let mut max_drift = 0.0f64;
                for _ in 0..50 {
                    core.step_n(2_000);
                    for (idx, &(a, b)) in pairs.iter().enumerate() {
                        let d = (com(&core, a) - com(&core, b)).length();
                        max_drift = max_drift.max(((d - d0[idx]) / d0[idx]).abs());
                    }
                }
                let w_ax = (0..n_alphas)
                    .map(|ai| {
                        let a = &core.groups[gid].alphas[ai];
                        a.angular_velocity.dot(a.orientation * DVec3::Y).abs()
                    })
                    .sum::<f64>()
                    / n_alphas as f64;
                let finite = core.particles.iter().all(|p| p.position.is_finite());
                println!(
                    "RigidAlpha carbon quiet, 100k steps: max_drift={:.1}% \
                     KE={:.4e} roll={w_ax:.4} finite={finite}",
                    max_drift * 100.0,
                    core.total_kinetic_energy(),
                );
            }

            // (c) Flyby + swat at the shipped couplings defaults.
            {
                let def = crate::atom_core::Couplings::default();
                let (drift, finite) = run_flyby_scenario_anatomy(
                    def.channeling,
                    def.nuclear_ambient,
                    15_000,
                    15_000,
                    500,
                    anatomy,
                );
                println!(
                    "Flyby/swat at defaults: max_drift={:.1}% finite={finite} \
                     (limit {:.0}%)",
                    drift * 100.0,
                    FLYBY_DRIFT_LIMIT * 100.0
                );
            }

            // (d) Seed-phase robustness: 4 seeds x 200k quiet steps.
            {
                print!("Transient seeds (200k steps): ");
                for k in 0..4usize {
                    let mut core = standard_core();
                    apply_env_overrides(&mut core);
                    core.post_anatomy = anatomy;
                    burn_seed_offset(&mut core, k);
                    let gid = core
                        .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                        .expect("carbon preset");
                    core.set_nucleus_dynamics(
                        crate::atom_core::NucleusDynamics::RigidAlpha,
                    );
                    let n_alphas = core.groups[gid].alphas.len();
                    let com =
                        |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                    let pairs: Vec<(usize, usize)> = (0..n_alphas)
                        .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                        .collect();
                    let d0: Vec<f64> = pairs
                        .iter()
                        .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                        .collect();
                    let mut max_drift = 0.0f64;
                    for _ in 0..40 {
                        core.step_n(5_000);
                        for (idx, &(a, b)) in pairs.iter().enumerate() {
                            let d = (com(&core, a) - com(&core, b)).length();
                            max_drift =
                                max_drift.max(((d - d0[idx]) / d0[idx]).abs());
                        }
                    }
                    print!("k={k}:{:.1}% ", max_drift * 100.0);
                }
                println!();
            }
        }
    }

    /// Session-32 plugged-carbon stability report: the corrected carbon
    /// (2 core alphas + a proton/neutron plug pair each pole) in
    /// RigidAlpha across 4 seed phases. Plugs are single-member alpha
    /// units — independent force-held bodies — so this measures both
    /// core-stack binding AND plug retention (a brand-new question; the
    /// old 3-stack had no plugs). Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_carbon_stability --nocapture`
    #[test]
    #[ignore]
    fn report_carbon_stability() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        for k in 0..4usize {
            let mut core = standard_core();
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();

            let mut worst = 0.0f64;
            let mut worst_pair = (0usize, 0usize);
            for _ in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                for (idx, &(a, b)) in pairs.iter().enumerate() {
                    let d = (com(&core, a) - com(&core, b)).length();
                    let drift = ((d - d0[idx]) / d0[idx]).abs();
                    if drift > worst {
                        worst = drift;
                        worst_pair = (a, b);
                    }
                }
            }
            let finite = core.particles.iter().all(|p| p.position.is_finite());
            // Alpha-unit sizes tell which pair kind drifted (4 = core
            // alpha, 1 = plug).
            let kind = |ai: usize| {
                if core.groups[gid].alphas[ai].members.len() == 4 {
                    "core"
                } else {
                    "plug"
                }
            };
            println!(
                "seed k={k}: worst_drift={:.1}% over {STEPS} steps \
                 (pair {}-{} = {}-{}) finite={finite}",
                worst * 100.0,
                worst_pair.0,
                worst_pair.1,
                kind(worst_pair.0),
                kind(worst_pair.1),
            );
        }
    }

    /// Session-32 plug-retention diagnostic: WHERE does the plugged
    /// carbon go, and WHAT force ejects it? Per-alpha COM trajectories
    /// (cylindrical about the group COM) plus rest-pose pair force
    /// breakdowns for every plug↔socket-proton, plug↔plug-partner and
    /// core↔core pair — the same anatomy view that found the post
    /// interpenetration bug in Phase B. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_retention --nocapture`
    #[test]
    #[ignore]
    fn report_plug_retention() {
        // ── Rest-pose force anatomy ──
        let mut core = standard_core();
        apply_env_overrides(&mut core);
        println!(
            "knobs: gap={} suction={} emit_scale={} tension={} align={}",
            core.plug_pair_gap,
            core.flow_suction,
            core.flow_emit_scale,
            core.flow_tension,
            core.flow_align,
        );
        let gid = core
            .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
        core.running = true;
        core.step_n(200);

        let alphas = &core.groups[gid].alphas;
        let label = |core: &AtomCore, ai: usize| -> String {
            let a = &core.groups[gid].alphas[ai];
            if a.members.len() == 4 {
                format!("core{ai}")
            } else {
                let side = if a.com.y < 0.0 { "-" } else { "+" };
                format!("pair{side}")
            }
        };
        // Pairs of interest (session-34 fused pairs): each plug-pair
        // MEMBER vs its nearest core proton (the socket), plus the
        // intra-pair row — which compute_forces now SKIPS as pre-fused
        // in RigidAlpha; the breakdown still evaluates it, showing what
        // fusion is suppressing.
        let mut pairs: Vec<(usize, usize, String)> = Vec::new();
        let core_protons: Vec<usize> = alphas
            .iter()
            .filter(|a| a.members.len() == 4)
            .flat_map(|a| a.members.clone())
            .filter(|&m| core.profiles[core.particles[m].profile_id].name == "proton")
            .collect();
        let plug_alphas: Vec<usize> = (0..alphas.len())
            .filter(|&ai| alphas[ai].members.len() == 2)
            .collect();
        for &ai in &plug_alphas {
            let members = core.groups[gid].alphas[ai].members.clone();
            for &m in &members {
                let initial = core.profiles[core.particles[m].profile_id].name
                    [..1]
                    .to_uppercase();
                let socket = core_protons
                    .iter()
                    .copied()
                    .min_by(|&a, &b| {
                        core.pair_distance(m, a)
                            .partial_cmp(&core.pair_distance(m, b))
                            .unwrap()
                    })
                    .unwrap();
                pairs.push((
                    m,
                    socket,
                    format!("{}.{initial}↔socket", label(&core, ai)),
                ));
            }
            pairs.push((
                members[0],
                members[1],
                format!("{} intra(SKIPPED)", label(&core, ai)),
            ));
        }

        println!(
            "\n-- carbon rest-pose pair anatomy (RigidAlpha, t≈200 steps) --\n\
             {:>18} {:>6} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9}",
            "pair", "r", "channel", "grav", "charge", "ambient", "intake", "stream", "contact", "tension"
        );
        for (i, j, name) in &pairs {
            let fb = core.pair_force_breakdown(*i, *j);
            println!(
                "{name:>18} {:>6.3} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>+9.4}",
                fb[0], fb[1], fb[2], fb[3], fb[4], fb[5], fb[6], fb[7], fb[8],
            );
        }

        // ── Trajectories: seeds 0 and 2 (plug-plug and core-core worst
        // cases in report_carbon_stability) ──
        for k in [0usize, 2] {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            core.running = true;
            let n_alphas = core.groups[gid].alphas.len();
            let labels: Vec<String> =
                (0..n_alphas).map(|ai| label(&core, ai)).collect();
            println!(
                "\n== seed k={k} per-alpha (lat, y) about group COM ==\n{:>6} {:>10} {}",
                "step",
                "KE",
                labels
                    .iter()
                    .map(|l| format!("{l:>14}"))
                    .collect::<String>(),
            );
            for s in 0..15 {
                core.step_n(20_000);
                let com: DVec3 = (0..n_alphas)
                    .map(|ai| core.groups[gid].alphas[ai].com)
                    .sum::<DVec3>()
                    / n_alphas as f64;
                let row: String = (0..n_alphas)
                    .map(|ai| {
                        let d = core.groups[gid].alphas[ai].com - com;
                        format!(
                            "{:>14}",
                            format!("({:.2},{:+.2})", (d.x * d.x + d.z * d.z).sqrt(), d.y)
                        )
                    })
                    .collect();
                println!(
                    "{:>6} {:>10.4} {row}",
                    (s + 1) * 20_000,
                    core.total_kinetic_energy(),
                );
            }
        }
    }

    /// Session-32 post-sign-fix calibration: flow_tension (align =
    /// tension/2, the earned ratio) swept against plugged-carbon
    /// retention (4 seeds × 300k) AND tri_alpha quiet drift — the fixed
    /// sign doubles tension on symmetric links instead of cancelling it,
    /// so the old tension=1.0 default must be re-earned, not assumed. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_tension_sweep --nocapture`
    #[test]
    #[ignore]
    fn report_plug_tension_sweep() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        let drift_run = |preset: &str, tension: f64, k: usize| -> (f64, f64, bool) {
            let mut core = standard_core();
            apply_env_overrides(&mut core); // honors CFM_GAP for the
                                            // post-pump-kill retest
            core.flow_tension = tension;
            core.flow_align = tension * 0.5;
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();
            let mut worst = 0.0f64;
            core.running = true;
            for _ in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                for (idx, &(a, b)) in pairs.iter().enumerate() {
                    let d = (com(&core, a) - com(&core, b)).length();
                    worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                }
            }
            let finite = core.particles.iter().all(|p| p.position.is_finite());
            (worst, core.total_kinetic_energy(), finite)
        };

        println!(
            "\n{:>8} | {:>40} | {:>40}",
            "tension", "carbon worst drift (k=0..3)", "tri_alpha worst drift (k=0..3)"
        );
        for &tension in &[0.5, 1.0, 2.0, 4.0, 8.0] {
            let fmt = |preset: &str| -> String {
                (0..4)
                    .map(|k| {
                        let (w, ke, finite) = drift_run(preset, tension, k);
                        format!(
                            "{:>6.0}%{}(KE {ke:.0})",
                            w * 100.0,
                            if finite { "" } else { "!" }
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            println!("{tension:>8.1} | {:>40} | {:>40}", fmt("carbon"), fmt("tri_alpha"));
        }
    }

    /// Session-32 energy-pump ablation: the plugged carbon gains KE
    /// unboundedly (up to ~4000 over 300k) at EVERY flow_tension while
    /// tri_alpha stays cold at every one — some force term injects
    /// energy only in plug configurations. Toggle the suspects one at a
    /// time and watch end-state KE + worst drift. Known-suspect notes:
    /// vortex has no reaction force (Opus audit); the flow solver is 8
    /// steps stale, making tension slightly non-conservative; align
    /// torques act on near-zero-inertia 1-nucleon plug alphas. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_energy_ablation --nocapture`
    #[test]
    #[ignore]
    fn report_plug_energy_ablation() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        type Ablate = fn(&mut AtomCore);
        let cases: &[(&str, Ablate)] = &[
            ("baseline", |_c| {}),
            ("align=0", |c| c.flow_align = 0.0),
            ("tension=0 align=0", |c| {
                c.flow_tension = 0.0;
                c.flow_align = 0.0;
            }),
            ("vortex=0", |c| c.couplings.vortex = 0.0),
            ("corot=0", |c| c.couplings.corot = 0.0),
            ("torque=0", |c| c.couplings.torque = 0.0),
            ("vortex=0 align=0", |c| {
                c.couplings.vortex = 0.0;
                c.flow_align = 0.0;
            }),
            ("vortex=0 corot=0", |c| {
                c.couplings.vortex = 0.0;
                c.couplings.corot = 0.0;
            }),
            // Session-34: flow-solver staleness as an injector — the
            // Phase B forces read amplitudes up to 8 steps old, and the
            // no-starve network roughly doubled those amplitudes.
            ("solve_every=1", |c| c.flow_solve_every = 1),
            ("solve_every=1 vortex=0", |c| {
                c.flow_solve_every = 1;
                c.couplings.vortex = 0.0;
            }),
            // Session-34 round 2: the first matrix left KE elevated in
            // EVERY row — the never-ablated velocity-dependent suspect
            // is the doppler factor (asymmetric clamp 0.2..5.0
            // rectifies pair oscillation into net heating; the plug
            // pair at r≈0.9 is the perfect rectifier). intake included
            // for completeness.
            ("drag=0", |c| c.couplings.drag = 0.0),
            ("drag=0 vortex=0 corot=0", |c| {
                c.couplings.drag = 0.0;
                c.couplings.vortex = 0.0;
                c.couplings.corot = 0.0;
            }),
            ("intake=0", |c| c.couplings.intake = 0.0),
            ("drag=0 torque=0", |c| {
                c.couplings.drag = 0.0;
                c.couplings.torque = 0.0;
            }),
        ];

        println!(
            "\n{:>20} | {:>26} | {:>26}",
            "ablation", "k=0 drift / KE", "k=1 drift / KE"
        );
        for (name, ablate) in cases {
            let mut cols = Vec::new();
            for k in 0..2usize {
                let mut core = standard_core();
                apply_env_overrides(&mut core); // CFM_GAP etc.
                ablate(&mut core);
                burn_seed_offset(&mut core, k);
                let gid = core
                    .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("carbon preset");
                core.set_nucleus_dynamics(
                    crate::atom_core::NucleusDynamics::RigidAlpha,
                );
                let n_alphas = core.groups[gid].alphas.len();
                let com =
                    |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                let pairs: Vec<(usize, usize)> = (0..n_alphas)
                    .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                    .collect();
                let d0: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                let mut worst = 0.0f64;
                core.running = true;
                for _ in 0..(STEPS / SAMPLE_EVERY) {
                    core.step_n(SAMPLE_EVERY);
                    for (idx, &(a, b)) in pairs.iter().enumerate() {
                        let d = (com(&core, a) - com(&core, b)).length();
                        worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                    }
                }
                cols.push(format!(
                    "{:>7.0}% / {:>10.1}",
                    worst * 100.0,
                    core.total_kinetic_energy()
                ));
            }
            println!("{name:>20} | {:>26} | {:>26}", cols[0], cols[1]);
        }
    }

    /// Session-33 torque-war diagnostic: the two orientation authorities
    /// on each carbon plug — the "equator toward charge" gear-mesh torque
    /// and the flow-align torque — summed over all same-group partners,
    /// tracked with the plug's pole orientation and its alpha's angular
    /// velocity. The energy-pump ablation showed retention is best with
    /// the charge torque OFF and much worse with flow-align off; this
    /// report shows whether they actually fight (opposed directions,
    /// comparable magnitudes) and at what orientation the war breaks out.
    /// Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_torque_war --nocapture`
    #[test]
    #[ignore]
    fn report_plug_torque_war() {
        let mut core = standard_core();
        let gid = core
            .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
        core.running = true;

        // Session-34 fused pairs: the plug is one 2-member alpha; track
        // its proton member (members[0]) as the orientation probe.
        let plug_alphas: Vec<usize> = (0..core.groups[gid].alphas.len())
            .filter(|&ai| core.groups[gid].alphas[ai].members.len() == 2)
            .collect();
        let group_members: Vec<usize> = core.groups[gid].members.clone();
        // Rest pole per plug, captured at spawn (world frame).
        let rest_pole: Vec<DVec3> = plug_alphas
            .iter()
            .map(|&ai| {
                let m = core.groups[gid].alphas[ai].members[0];
                core.particles[m].pole_axis()
            })
            .collect();
        let label = |core: &AtomCore, ai: usize| -> String {
            let m = core.groups[gid].alphas[ai].members[0];
            let name = &core.profiles[core.particles[m].profile_id].name;
            let side = if core.groups[gid].alphas[ai].com.y < 0.0 { "-" } else { "+" };
            format!("{}{side}", &name[..1].to_uppercase())
        };
        let labels: Vec<String> =
            plug_alphas.iter().map(|&ai| label(&core, ai)).collect();

        println!(
            "\nper plug: tilt°(pole vs rest) |tau_charge| |tau_flow| cos(charge,flow) |w|"
        );
        println!(
            "{:>7} {}",
            "step",
            labels
                .iter()
                .map(|l| format!("{l:>38}"))
                .collect::<String>()
        );
        for s in 0..15 {
            core.step_n(if s == 0 { 200 } else { 20_000 });
            let row: String = plug_alphas
                .iter()
                .enumerate()
                .map(|(pi, &ai)| {
                    let m = core.groups[gid].alphas[ai].members[0];
                    let pole = core.particles[m].pole_axis();
                    let tilt = pole.dot(rest_pole[pi]).clamp(-1.0, 1.0).acos().to_degrees();
                    let mut tau_c = DVec3::ZERO;
                    let mut tau_f = DVec3::ZERO;
                    for &other in &group_members {
                        if other == m {
                            continue;
                        }
                        let tb = core.pair_torque_breakdown(m, other);
                        tau_c += tb[0];
                        tau_f += tb[2];
                    }
                    let cosang = if tau_c.length() > 1e-12 && tau_f.length() > 1e-12 {
                        tau_c.normalize().dot(tau_f.normalize())
                    } else {
                        0.0
                    };
                    let w = core.groups[gid].alphas[ai].angular_velocity.length();
                    format!(
                        "{:>38}",
                        format!(
                            "{tilt:>5.1}° {:.4} {:.4} {cosang:+.2} {w:.3}",
                            tau_c.length(),
                            tau_f.length(),
                        )
                    )
                })
                .collect();
            println!(
                "{:>7} {row}",
                if s == 0 { 200 } else { s * 20_000 },
            );
        }
    }

    /// Session-33 gyroscopic-stiffness sweep (wig.pdf mechanism, see
    /// `AtomCore::gyro_spin`): plugged-carbon retention and tri_alpha
    /// quiet drift vs the spin angular momentum knob. gyro=0 is the
    /// classical baseline (must reproduce report_carbon_stability).
    /// Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_gyro_sweep --nocapture`
    #[test]
    #[ignore]
    fn report_gyro_sweep() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        let drift_run = |preset: &str, gyro: f64, k: usize| -> (f64, f64, bool) {
            let mut core = standard_core();
            core.gyro_spin = gyro;
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();
            let mut worst = 0.0f64;
            core.running = true;
            for _ in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                for (idx, &(a, b)) in pairs.iter().enumerate() {
                    let d = (com(&core, a) - com(&core, b)).length();
                    worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                }
            }
            let finite = core.particles.iter().all(|p| p.position.is_finite());
            (worst, core.total_kinetic_energy(), finite)
        };

        println!(
            "\n{:>8} | {:>44} | {:>44}",
            "gyro", "carbon worst drift (k=0..3)", "tri_alpha worst drift (k=0..3)"
        );
        for &gyro in &[0.0, 0.5, 2.0, 8.0, 32.0, 128.0] {
            let fmt = |preset: &str| -> String {
                (0..4)
                    .map(|k| {
                        let (w, ke, finite) = drift_run(preset, gyro, k);
                        format!(
                            "{:>5.0}%{}(KE {ke:.0})",
                            w * 100.0,
                            if finite { "" } else { "!" }
                        )
                    })
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            println!(
                "{gyro:>8.1} | {:>44} | {:>44}",
                fmt("carbon"),
                fmt("tri_alpha")
            );
        }
    }

    /// Session-33 factored plug-experiment matrix: the three candidate
    /// mechanisms (plug-pair gap, close-range ambient saturation,
    /// align-to-stream torque target) swept ORTHOGONALLY against
    /// plugged-carbon retention, 4 seeds each — after serially stacking
    /// them one-at-a-time made carbon worse than the committed baseline
    /// (the whirl-saga lesson: attribute before adopting). gyro_spin is
    /// excluded (its sweep already rejected it). Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_matrix --nocapture`
    #[test]
    #[ignore]
    fn report_plug_matrix() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        println!(
            "\n{:>5} {:>6} {:>7} | {:>52}",
            "gap", "satur", "stream", "carbon worst drift / end KE (k=0..3)"
        );
        for &gap in &[0.8, 0.35] {
            for &sat in &[0.0, crate::atom_core::NUCLEON_PITCH] {
                for &stream in &[false, true] {
                    let cols: Vec<String> = (0..4)
                        .map(|k| {
                            let mut core = standard_core();
                            core.plug_pair_gap = gap;
                            core.ambient_sat_r = sat;
                            core.align_to_stream = stream;
                            burn_seed_offset(&mut core, k);
                            let gid = core
                                .spawn_preset(
                                    "carbon",
                                    DVec3::ZERO,
                                    DVec3::ZERO,
                                    DVec3::Y,
                                )
                                .expect("carbon preset");
                            core.set_nucleus_dynamics(
                                crate::atom_core::NucleusDynamics::RigidAlpha,
                            );
                            let n_alphas = core.groups[gid].alphas.len();
                            let com = |core: &AtomCore, ai: usize| {
                                core.groups[gid].alphas[ai].com
                            };
                            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                                .flat_map(|a| {
                                    ((a + 1)..n_alphas).map(move |b| (a, b))
                                })
                                .collect();
                            let d0: Vec<f64> = pairs
                                .iter()
                                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                                .collect();
                            let mut worst = 0.0f64;
                            core.running = true;
                            for _ in 0..(STEPS / SAMPLE_EVERY) {
                                core.step_n(SAMPLE_EVERY);
                                for (idx, &(a, b)) in pairs.iter().enumerate() {
                                    let d =
                                        (com(&core, a) - com(&core, b)).length();
                                    worst =
                                        worst.max(((d - d0[idx]) / d0[idx]).abs());
                                }
                            }
                            format!(
                                "{:>5.0}%/{:<6.0}",
                                worst * 100.0,
                                core.total_kinetic_energy()
                            )
                        })
                        .collect();
                    println!(
                        "{gap:>5.2} {sat:>6.2} {:>7} | {:>52}",
                        stream,
                        cols.join(" ")
                    );
                }
            }
        }
    }

    /// Session-33 follow-up to report_plug_matrix: its best cell
    /// (plug_pair_gap=0.35, no ambient saturation, line torque target)
    /// kills the pair collapse-bounce KE pump (end KE 5-24 vs 17-578)
    /// but the cooled plugs still slide off the sockets (137-466%) —
    /// the socket tension (~0.16) can't beat the plug proton's ~0.85
    /// outward charge push. Sweep Phase C1 ambient confinement on top:
    /// the old "not confinement-fixable" verdict (best ~109% @ 0.5)
    /// predates both the tension sign fix and the pump kill. tri_alpha
    /// swept alongside — confinement also compresses the bare stack and
    /// must not wreck it. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_confine --nocapture`
    #[test]
    #[ignore]
    fn report_plug_confine() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        let run = |preset: &str, confine: f64, k: usize| -> (f64, f64) {
            let mut core = standard_core();
            core.plug_pair_gap = 0.35;
            core.ambient_confine = confine;
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();
            let mut worst = 0.0f64;
            core.running = true;
            for _ in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                for (idx, &(a, b)) in pairs.iter().enumerate() {
                    let d = (com(&core, a) - com(&core, b)).length();
                    worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                }
            }
            (worst, core.total_kinetic_energy())
        };

        println!(
            "\n{:>8} | {:>52} | {:>52}",
            "confine",
            "carbon (gap 0.35) worst drift / KE (k=0..3)",
            "tri_alpha worst drift / KE (k=0..3)"
        );
        for &confine in &[0.0, 0.25, 0.5, 1.0, 2.0, 4.0] {
            let fmt = |preset: &str| -> String {
                (0..4)
                    .map(|k| {
                        let (w, ke) = run(preset, confine, k);
                        format!("{:>5.0}%/{:<6.0}", w * 100.0, ke)
                    })
                    .collect::<Vec<_>>()
                    .join(" ")
            };
            println!(
                "{confine:>8.2} | {:>52} | {:>52}",
                fmt("carbon"),
                fmt("tri_alpha")
            );
        }
    }

    /// Session-34 Phase B.2: throughput-scaled suction sweep — the
    /// plug-retention missing piece (phos.pdf: a FED socket has
    /// suction; the static intake term measures ~0.016 at the carbon
    /// socket regardless of live flow). flow_suction scales the intake
    /// pull + channeling force by each body's live FlowState::mult;
    /// flow_emit_scale does the same to the charge push ("a fed funnel
    /// pushes harder" — expected to oppose retention, swept alongside
    /// so the matrix attributes both). Runs at the report_plug_matrix
    /// best cell (plug_pair_gap=0.35 — pump dead, plugs slide off
    /// quietly at 137-466%). tri_alpha swept as the do-no-harm control.
    /// Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_suction_sweep --nocapture`
    #[test]
    #[ignore]
    fn report_plug_suction_sweep() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;

        let run = |preset: &str, suction: f64, emit: f64, k: usize| -> (f64, f64) {
            let mut core = standard_core();
            core.plug_pair_gap = 0.35;
            core.flow_suction = suction;
            core.flow_emit_scale = emit;
            burn_seed_offset(&mut core, k);
            let gid = core
                .spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            let n_alphas = core.groups[gid].alphas.len();
            let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
            let pairs: Vec<(usize, usize)> = (0..n_alphas)
                .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                .collect();
            let d0: Vec<f64> = pairs
                .iter()
                .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                .collect();
            let mut worst = 0.0f64;
            core.running = true;
            for _ in 0..(STEPS / SAMPLE_EVERY) {
                core.step_n(SAMPLE_EVERY);
                for (idx, &(a, b)) in pairs.iter().enumerate() {
                    let d = (com(&core, a) - com(&core, b)).length();
                    worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                }
            }
            (worst, core.total_kinetic_energy())
        };

        println!(
            "\n{:>8} {:>5} | {:>52} | {:>52}",
            "suction",
            "emit",
            "carbon (gap 0.35) worst drift / KE (k=0..3)",
            "tri_alpha worst drift / KE (k=0..3)"
        );
        for &suction in &[0.0, 1.0, 2.0, 4.0, 8.0] {
            for &emit in &[0.0, 1.0] {
                if suction == 0.0 && emit != 0.0 {
                    continue;
                }
                let fmt = |preset: &str| -> String {
                    (0..4)
                        .map(|k| {
                            let (w, ke) = run(preset, suction, emit, k);
                            format!("{:>5.0}%/{:<6.0}", w * 100.0, ke)
                        })
                        .collect::<Vec<_>>()
                        .join(" ")
                };
                println!(
                    "{suction:>8.1} {emit:>5.1} | {:>52} | {:>52}",
                    fmt("carbon"),
                    fmt("tri_alpha")
                );
            }
        }
    }

    /// Session-34 static well scan: does the carbon plug pair sit in a
    /// POTENTIAL WELL at all? The full energy-pump ablation matrix
    /// (report_plug_energy_ablation) showed no single-term injector and
    /// no retaining cell — if the static force landscape is monotone
    /// outward, retention is impossible regardless of dynamics quality
    /// and the campaign must move to missing physics instead of pump
    /// hunting. Displaces the +y plug pair rigidly (axially +
    /// laterally), converges the flow network at each geometry, and
    /// prints the net force component along the displacement (negative
    /// = restoring). Swept over (suction, tension) knob combos. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_plug_well --nocapture`
    #[test]
    #[ignore]
    fn report_plug_well() {
        let deltas: &[f64] = &[
            -2.4, -2.0, -1.6, -1.2, -0.8, -0.4, -0.2, 0.0, 0.15, 0.3, 0.5,
            0.8, 1.2, 1.6, 2.4, 3.2,
        ];
        println!(
            "\n{:>8} {:>8} {:>5} | {}",
            "suction",
            "tension",
            "axis",
            deltas
                .iter()
                .map(|d| format!("{d:>8.2}"))
                .collect::<String>()
        );
        for &(suction, tension) in &[
            (0.0, 1.0),
            (4.0, 1.0),
            (8.0, 1.0),
            (0.0, 4.0),
            (0.0, 8.0),
            (4.0, 4.0),
            (8.0, 8.0),
        ] {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            core.plug_pair_gap = 0.35;
            core.flow_suction = suction;
            core.flow_tension = tension;
            core.flow_align = tension * 0.5;
            let gid = core
                .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
            // Converge the flow network at the rest geometry without
            // moving anything.
            for _ in 0..40 {
                core.solve_charge_flow();
            }
            let pos0: Vec<DVec3> =
                core.particles.iter().map(|p| p.position).collect();
            // The +y plug pair: members of 1-nucleon alphas with com.y > 0.
            // Session-34 fused pair: the +y plug is ONE 2-member alpha.
            let plug_members: Vec<usize> = (0..core.groups[gid].alphas.len())
                .filter(|&ai| {
                    core.groups[gid].alphas[ai].members.len() == 2
                        && core.groups[gid].alphas[ai].com.y > 0.0
                })
                .flat_map(|ai| core.groups[gid].alphas[ai].members.clone())
                .collect();
            assert_eq!(plug_members.len(), 2, "expected p+n plug pair on +y");

            // Per-member force vectors at the rest pose (δ=0) — the
            // lateral scan found a huge displacement-flat force that
            // must be identified, not guessed at.
            if suction == 0.0 && tension == 1.0 {
                for (i, p) in core.particles.iter_mut().enumerate() {
                    p.position = pos0[i];
                    p.velocity = DVec3::ZERO;
                }
                for _ in 0..40 {
                    core.solve_charge_flow();
                }
                core.compute_forces();
                for &m in &plug_members {
                    let p = &core.particles[m];
                    println!(
                        "  rest member {m} ({}) pos=({:+.2},{:+.2},{:+.2}) pole=({:+.2},{:+.2},{:+.2}) F=({:+.3},{:+.3},{:+.3})",
                        core.profiles[p.profile_id].name,
                        p.position.x, p.position.y, p.position.z,
                        p.pole_axis().x, p.pole_axis().y, p.pole_axis().z,
                        p.force_accum.x, p.force_accum.y, p.force_accum.z,
                    );
                }
            }

            for (axis, dir) in [("axial", DVec3::Y), ("later", DVec3::X)] {
                let row: String = deltas
                    .iter()
                    .map(|&delta| {
                        for (i, p) in core.particles.iter_mut().enumerate() {
                            p.position = pos0[i];
                            p.velocity = DVec3::ZERO;
                        }
                        for &m in &plug_members {
                            core.particles[m].position += dir * delta;
                        }
                        for _ in 0..40 {
                            core.solve_charge_flow();
                        }
                        core.compute_forces();
                        let f: f64 = plug_members
                            .iter()
                            .map(|&m| core.particles[m].force_accum.dot(dir))
                            .sum();
                        format!("{f:>8.3}")
                    })
                    .collect();
                println!("{suction:>8.1} {tension:>8.1} {axis:>5} | {row}");
            }
        }
    }

    /// Session-35 momentum + energy flux audit. The papers say charge is an
    /// OPEN field — "the charge field, by itself, doesn't conserve energy"
    /// (cc.pdf) — and its momentum "is always radially out from the center"
    /// (pause.html). Our sim is CLOSED: it has no channel for charge to carry
    /// momentum or energy out, so (a) the proton→neutron emission asymmetry
    /// and the reaction-less vortex leave a NET force on an isolated nucleus
    /// (the in-app carbon core drifts off-camera), and (b) the only
    /// velocity-dependent term is a Doppler MULTIPLIER on the conservative
    /// charge force that can INJECT energy, with no −k·v sink to remove it.
    ///
    /// This measures both, so the fix is chosen from numbers:
    ///   PART A — with the uniform ambient globals zeroed (true isolation),
    ///   Σ force_accum over the whole nucleus at the rest pose. For an
    ///   isolated body this MUST be ≈0. Corot and flow-tension are applied as
    ///   equal-and-opposite pairs, so ablating them must NOT change the sum
    ///   (a built-in consistency check); only the un-paired vortex and the
    ///   central charge/intake asymmetry can move it. Reports the residual
    ///   vector after each ablation → attribution.
    ///   PART B — a normal dynamics run (default app conditions: ambient on,
    ///   plug_orient_lock off), sampling total KE, nucleus COM drift, and COM
    ///   speed. Shows the pump climbing with no counterpart.
    /// Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_momentum_energy_audit --nocapture`
    #[test]
    #[ignore]
    fn report_momentum_energy_audit() {
        // Whole-nucleus aggregates from the group's flat member list.
        fn nucleus_members(core: &AtomCore, gid: usize) -> Vec<usize> {
            core.groups[gid].members.clone()
        }
        fn net_force(core: &AtomCore, members: &[usize]) -> DVec3 {
            members.iter().map(|&m| core.particles[m].force_accum).sum()
        }
        fn com(core: &AtomCore, members: &[usize]) -> DVec3 {
            let mut num = DVec3::ZERO;
            let mut den = 0.0;
            for &m in members {
                let mass = core.profiles[core.particles[m].profile_id].mass;
                num += core.particles[m].position * mass;
                den += mass;
            }
            num / den
        }
        fn com_vel(core: &AtomCore, members: &[usize]) -> DVec3 {
            let mut num = DVec3::ZERO;
            let mut den = 0.0;
            for &m in members {
                let mass = core.profiles[core.particles[m].profile_id].mass;
                num += core.particles[m].velocity * mass;
                den += mass;
            }
            num / den
        }

        // ---- PART A: static momentum-leak attribution (isolated) ----
        // Zero the uniform ambient globals so Σforce is purely the internal
        // pairwise budget — a uniform per-mass field is a real external
        // force, not a leak, and would swamp the residual we care about.
        println!("\n=== PART A: net force on isolated carbon at rest pose ===");
        println!("(ambient globals zeroed; an isolated body MUST sum to ~0)");
        type Ablate = fn(&mut AtomCore);
        let ablations: &[(&str, Ablate)] = &[
            ("full internal      ", |_c| {}),
            ("vortex=0           ", |c| c.couplings.vortex = 0.0),
            ("vortex=0 corot=0   ", |c| {
                c.couplings.vortex = 0.0;
                c.couplings.corot = 0.0;
            }),
            ("+tension=0 (=>resid)", |c| {
                c.couplings.vortex = 0.0;
                c.couplings.corot = 0.0;
                c.flow_tension = 0.0;
                c.flow_align = 0.0;
            }),
        ];
        // `raw leak` = group_self_force (the un-paired intra-nucleus net
        // force, BEFORE the session-35 mass-weighted cancellation).
        // `net ΣF` = Σ force_accum over the members AFTER cancellation — must
        // be ≈0 for the isolated body (verifies the fix).
        println!(
            "{:>21} | {:>10} | {:>10} | {:>28}",
            "ablation", "raw leak", "net ΣF", "raw leak vector (x,y,z)"
        );
        for (label, ablate) in ablations {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            core.plug_pair_gap = 0.35;
            let gid = core
                .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("carbon preset");
            core.set_nucleus_dynamics(
                crate::atom_core::NucleusDynamics::RigidAlpha,
            );
            core.ambient_gravity = DVec3::ZERO;
            core.ambient_charge = DVec3::ZERO;
            ablate(&mut core);
            for _ in 0..40 {
                core.solve_charge_flow();
            }
            core.compute_forces();
            let members = nucleus_members(&core, gid);
            let raw = core.group_self_force[gid];
            let net = net_force(&core, &members);
            println!(
                "{label} | {:>10.4} | {:>10.4} | ({:+9.4},{:+9.4},{:+9.4})",
                raw.length(),
                net.length(),
                raw.x,
                raw.y,
                raw.z
            );
        }

        // ---- PART B: dynamics KE + COM-drift trace (default app cond.) ----
        println!("\n=== PART B: KE + COM drift over a normal run (defaults) ===");
        println!(
            "{:>10} | {:>12} | {:>12} | {:>12}",
            "step", "total_KE", "COM_drift", "COM_speed"
        );
        const STEPS: usize = 120_000;
        const SAMPLE_EVERY: usize = 10_000;
        let mut core = standard_core();
        apply_env_overrides(&mut core);
        core.plug_pair_gap = 0.35;
        let gid = core
            .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);
        let members = nucleus_members(&core, gid);
        let com0 = com(&core, &members);
        for s in 0..=STEPS {
            if s % SAMPLE_EVERY == 0 {
                let drift = (com(&core, &members) - com0).length();
                let speed = com_vel(&core, &members).length();
                println!(
                    "{s:>10} | {:>12.3} | {:>12.4} | {:>12.6}",
                    core.total_kinetic_energy(),
                    drift,
                    speed
                );
            }
            core.step();
        }
    }

    /// Session-35 "upright + stable" probe. User feedback on the in-app
    /// carbon: the nucleus is "hooked to a mad gyroscope and refuses to
    /// align with the ambient field". The papers say a nucleus turns to
    /// "stand straight up, to align to E" (dielec.pdf) in a DIRECTIONAL
    /// ambient field (pasta.pdf: Earth's field "is pointing straight up
    /// everywhere"). This measures whether `apply_ambient_alignment` gives
    /// the tumbling nucleus that orientation reference, and how RigidAlpha
    /// (free alphas, emergent carousel) compares with RigidLock (rigid
    /// fused body — uf4.pdf: "alphas can't be broken and rearranged").
    ///
    /// Ambient field tilted 45° off the spawn axis so alignment has real
    /// work to do. Columns: KE (churn), upright = |n̂·â| (1 = principal
    /// axis parallel to the field), spread = max |alpha_com − group_com|
    /// (shape: a rigid nucleus holds it, a churning one grows it). Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_upright --nocapture`
    #[test]
    #[ignore]
    fn report_upright() {
        use crate::atom_core::NucleusDynamics;
        // Principal (long) axis of a member cloud about its mean position.
        fn principal_axis(core: &AtomCore, members: &[usize], seed: DVec3) -> DVec3 {
            let mut c = DVec3::ZERO;
            for &m in members {
                c += core.particles[m].position;
            }
            c /= members.len() as f64;
            let offs: Vec<DVec3> =
                members.iter().map(|&m| core.particles[m].position - c).collect();
            let mut n = seed;
            for _ in 0..12 {
                let mut v = DVec3::ZERO;
                for o in &offs {
                    v += *o * o.dot(n);
                }
                let vl = v.length();
                if vl < 1e-12 {
                    break;
                }
                n = v / vl;
            }
            n
        }

        let amb_dir = DVec3::new(0.0, 1.0, 1.0).normalize(); // 45° off +Y spawn
        for (mode_name, mode) in
            [("RigidAlpha", NucleusDynamics::RigidAlpha), ("RigidLock", NucleusDynamics::RigidLock)]
        {
            for &align in &[0.0, 20.0, 80.0] {
                let mut core = standard_core();
                apply_env_overrides(&mut core);
                core.plug_pair_gap = 0.35;
                core.ambient_charge_dir = amb_dir;
                core.ambient_align = align;
                let gid = core
                    .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("carbon preset");
                core.set_nucleus_dynamics(mode);
                let members = core.groups[gid].members.clone();
                core.running = true;
                println!(
                    "\n== {mode_name}  ambient_align={align}  (field 45° off spawn axis) ==\n\
                     {:>7} | {:>9} | {:>7} | {:>7}",
                    "step", "KE", "upright", "spread"
                );
                let mut n_hat = principal_axis(&core, &members, amb_dir);
                for s in 0..=200_000 {
                    if s % 20_000 == 0 {
                        n_hat = principal_axis(&core, &members, n_hat);
                        let upright = n_hat.dot(amb_dir).abs();
                        let mut c = DVec3::ZERO;
                        for &m in &members {
                            c += core.particles[m].position;
                        }
                        c /= members.len() as f64;
                        let spread = members
                            .iter()
                            .map(|&m| (core.particles[m].position - c).length())
                            .fold(0.0f64, f64::max);
                        println!(
                            "{s:>7} | {:>9.3} | {:>7.3} | {:>7.3}",
                            core.total_kinetic_energy(),
                            upright,
                            spread
                        );
                    }
                    core.step();
                }
            }
        }
    }

    /// Session-35 element FORCE MAP — the "forces coming in and out of a
    /// particular element" (user's Phase-2 aim: lay oxygen in open space,
    /// have it form O₂ on its own). Probe the force a lone test proton (the
    /// universal plug) feels at a grid of positions around a RigidLock,
    /// N-S-aligned nucleus. The radial force sign IS the map: NEGATIVE =
    /// pulled IN (a polar intake socket — where a partner atom plugs in and
    /// bonds); POSITIVE = pushed OUT (equatorial emission disc — the
    /// standoff/repel). The probe faces the socket (pole toward center),
    /// like an approaching plug. RigidLock so the source is a clean, stable
    /// rigid field; the probe↔nucleus pairs are INTER-group, so un-skipped
    /// and fully real — this is exactly the chemistry force, not a faked
    /// motion. Ambient globals zeroed to isolate the element's own field.
    /// Run: `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_element_field_map --nocapture`
    #[test]
    #[ignore]
    fn report_element_field_map() {
        use crate::atom_core::NucleusDynamics;
        let thetas_deg: &[f64] = &[0.0, 15.0, 30.0, 45.0, 60.0, 90.0];
        let radii: &[f64] = &[3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0];
        for element in ["carbon", "oxygen"] {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            core.plug_pair_gap = 0.35;
            core.ambient_gravity = DVec3::ZERO;
            core.ambient_charge = DVec3::ZERO;
            core.spawn_preset(element, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(NucleusDynamics::RigidLock);
            for _ in 0..40 {
                core.solve_charge_flow();
            }
            let p_id = core.profile_id_by_name("proton").unwrap();
            let probe = core
                .spawn_particle_ex(
                    p_id,
                    DVec3::new(0.0, 6.0, 0.0),
                    DVec3::ZERO,
                    -DVec3::Y,
                    0.0,
                )
                .expect("probe");
            println!(
                "\n== {element}: radial force on a test proton  (− = intake/bond site, + = emit/repel) ==\n\
                 {:>11} | {}",
                "θ from pole",
                radii
                    .iter()
                    .map(|r| format!("{:>8}", format!("r={r:.1}")))
                    .collect::<String>()
            );
            for &th_deg in thetas_deg {
                let th = th_deg.to_radians();
                // +Y = pole axis, +X = equator.
                let dir = DVec3::new(th.sin(), th.cos(), 0.0);
                let row: String = radii
                    .iter()
                    .map(|&r| {
                        let pos = dir * r;
                        // Distance to nearest nucleus body: inside contact
                        // range the hard-sphere spring dominates and the
                        // number is meaningless — mark it buried instead.
                        let min_d = core
                            .particles
                            .iter()
                            .take(core.particles.len() - 1) // exclude the probe
                            .map(|p| (p.position - pos).length())
                            .fold(f64::INFINITY, f64::min);
                        if min_d < 2.4 {
                            return format!("{:>8}", "·buried");
                        }
                        core.particles[probe].position = pos;
                        core.particles[probe].orientation =
                            glam::DQuat::from_rotation_arc(DVec3::Y, -dir);
                        core.particles[probe].velocity = DVec3::ZERO;
                        core.compute_forces();
                        let f_r = core.particles[probe].force_accum.dot(dir);
                        format!("{f_r:>8.2}")
                    })
                    .collect();
                println!("{th_deg:>10.0}° | {row}");
            }
        }
    }

    /// Session-35 flow readout: per-nucleon throughput of a lone nucleus.
    /// The "make them suck" bonding vortex scales the intake pull by each
    /// emitter's `FlowState::mult` (output ÷ free-field baseline). If the
    /// pole plug that channels the whole stack's through-charge reads
    /// mult ≈ 1, it has NO amplified vortex to reach a partner — the gap the
    /// user flagged ("a proton on the pole of a larger stack ought to feel a
    /// vortex bigger than its own spin"). Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_nucleus_flow --nocapture`
    #[test]
    #[ignore]
    fn report_nucleus_flow() {
        use crate::atom_core::NucleusDynamics;
        for element in ["carbon", "oxygen"] {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            core.plug_pair_gap = 0.35;
            let gid = core
                .spawn_preset(element, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
            for _ in 0..60 {
                core.solve_charge_flow();
            }
            println!(
                "\n== {element} per-nucleon flow (RigidAlpha) ==\n\
                 {:>4} {:>8} {:>7} {:>6} {:>7} {:>7} {:>7}",
                "id", "name", "y", "mult", "intake", "out_pol", "lateral"
            );
            let members = core.groups[gid].members.clone();
            for m in members {
                let p = &core.particles[m];
                let f = &p.flow;
                println!(
                    "{m:>4} {:>8} {:>7.2} {:>6.2} {:>7.2} {:>7.2} {:>7.2}",
                    core.profiles[p.profile_id].name,
                    p.position.y,
                    f.mult,
                    f.intake[0] + f.intake[1],
                    f.through(),
                    f.lateral,
                );
            }
        }
    }

    /// Session-35 OH bond probe (user: "see how OH fares and work up from
    /// there"). A hydrogen atom approaches the +Y pole of an oxygen nucleus
    /// with a gentle inward velocity. Tests the user's DYNAMIC bond
    /// mechanism: the field draws H into O's polar VORTEX (far stronger than
    /// H's own spin, reaching past the plug proton's disc), H and O spin into
    /// sync, and they hold as a unit — LOOSER than fusion (a breakable
    /// directional link), which a static bare-proton probe (see
    /// report_element_field_map) cannot show because the sync force is
    /// velocity-dependent (`corot`). Columns: O–H distance, lateral offset
    /// from the pole axis, azimuth about the axis (advancing ⇒ captured,
    /// orbiting the vortex), H axial spin. On-/off-axis start × suction
    /// off/on (throughput may amplify the socket — the plug proton channels
    /// the whole stack's intake). Bond ⇔ d settles bounded, not repelled
    /// (d grows) nor fused (d→contact ~2). Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_oh_bond --nocapture`
    #[test]
    #[ignore]
    fn report_oh_bond() {
        use crate::atom_core::NucleusDynamics;
        // Spawn OUTSIDE the repulsive pole wall. report_oh_wellscan showed
        // a giant charge/stream wall for d<~9.75 (dPole<3.4); the only
        // attractive region — the through-charge tension shelf — sits at
        // d≈10-13 (dPole≈4-6). The old test started H at d=8 (INSIDE the
        // wall), so it was blasted out before it could fall into the well.
        let approaches = [
            ("on-axis ", DVec3::new(0.0, 12.5, 0.0)),
            ("off-axis", DVec3::new(3.0, 12.0, 0.0)),
        ];
        // Pure tension sweep. bond_suction is INERT at the pole (terminal
        // nucleon mult=1 → 1+gain·(mult−1)=1; wellscan: suck8/20/50
        // identical) and intake is 100-1000× weaker than the wall anyway.
        // Tension is the only lever that makes a well, so sweep it alone.
        // Gain is O(1) (flow_tension_pair already carries O(1) couplings);
        // watch for the stiff-torque blowup at the top of the range.
        let combos = [
            ("base ", 0.0f64, 0.0f64),
            ("t2   ", 0.0, 2.0),
            ("t4   ", 0.0, 4.0),
            ("t6   ", 0.0, 6.0),
        ];
        for (cname, bond_suction, bond_tension) in combos {
            for (name, start) in approaches {
                let mut core = standard_core();
                apply_env_overrides(&mut core);
                core.plug_pair_gap = 0.35;
                core.bond_suction = bond_suction;
                core.bond_tension = bond_tension;
                // Default directional flow on for the diagnostic (charge
                // down −Y → +Y pole is the intake/in-line socket) unless an
                // env CFM_AMBFLOW already set it.
                if core.ambient_charge_dir == DVec3::ZERO {
                    core.ambient_charge_dir = DVec3::new(0.0, -2.0, 0.0);
                }
                let gid = core
                    .spawn_preset("oxygen", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("oxygen");
                core.set_nucleus_dynamics(NucleusDynamics::RigidLock);
                // H approaches pole-first (proton pole toward O, electron on
                // the far side), with a gentle inward velocity.
                let toward = (DVec3::ZERO - start).normalize();
                let (hp, he) =
                    spawn_formed_hydrogen(&mut core, start, toward, -toward, 1.0);
                // Very gentle nudge — the shelf is shallow (~1-2 units of
                // force); a fast approach just overshoots into the wall and
                // reflects back out. Let the well do the work.
                let v_in = toward * 0.01;
                core.particles[hp].velocity += v_in;
                core.particles[he].velocity += v_in;
                core.running = true;
                let ocom = core.groups[gid].com;
                println!(
                    "\n== OH  {cname}  {name}  suction={bond_suction} tension={bond_tension} ==\n\
                     {:>7} | {:>7} | {:>6} | {:>7} | {:>7} | {:>6}",
                    "step", "O-H d", "lat", "azim°", "H_spin", "e_d"
                );
                let mut d_min = f64::INFINITY;
                let mut d_final = 0.0f64;
                for s in 0..=300_000 {
                    if s % 25_000 == 0 {
                        let rel = core.particles[hp].position - ocom;
                        let d = rel.length();
                        let lat = (rel.x * rel.x + rel.z * rel.z).sqrt();
                        let az = rel.z.atan2(rel.x).to_degrees();
                        let spin = core.particles[hp]
                            .angular_velocity
                            .dot(core.particles[hp].pole_axis());
                        let e_d = (core.particles[he].position
                            - core.particles[hp].position)
                            .length();
                        println!(
                            "{s:>7} | {d:>7.2} | {lat:>6.2} | {az:>7.1} | {spin:>7.2} | {e_d:>6.2}"
                        );
                    }
                    core.step();
                    let d = (core.particles[hp].position - ocom).length();
                    d_min = d_min.min(d);
                    d_final = d;
                }
                // Verdict: H held on the tension shelf (well minimum d≈11,
                // dPole≈4.6) vs drifted off. Started at d≈12.5, so "held"
                // means it stayed on the shelf (roughly 9-15) instead of
                // escaping. Wall reflection would send an overshoot to large
                // d, so d_final near the well = a real (weak/long) bond.
                let verdict = if d_final < 16.0 {
                    "HELD"
                } else {
                    "free"
                };
                println!(
                    "  -> {verdict}  d_final={d_final:.2}  d_min={d_min:.2}"
                );
            }
        }
    }

    /// Session-35 cont-2: STATIC axial well scan for the in-line OH bond.
    /// The dynamic report_oh_bond showed H stalling at d≈8 and drifting
    /// off — no attractive well. report_element_field_map hid the seat
    /// (it marks everything within 2.4 of a nucleon "·buried"). This walks
    /// a lone test proton (pole axial, absorption facing O) straight down
    /// the +Y intake pole INTO that region and prints the NET radial force
    /// + a coarse attribution, so we can see whether the intake vortex
    /// makes a well at the backed-out seat and what (if anything)
    /// bond_suction adds. Negative Fy = inward (toward O = a bond pull).
    /// A stable well = Fy < 0 for d beyond the seat, crossing to Fy > 0
    /// (contact/cushion) as d shrinks. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_oh_wellscan --nocapture`
    #[test]
    #[ignore]
    fn report_oh_wellscan() {
        use crate::atom_core::NucleusDynamics;
        // d = distance from O COM (origin) up the +Y axis. O's outer pole
        // nucleon sits ~6.5 out; scan from just above it out to clear space.
        let ds: Vec<f64> = (0..=26).map(|k| 6.75 + 0.25 * k as f64).collect();
        for (sname, suction, tension) in [
            ("base      ", 0.0f64, 0.0f64),
            ("suck8     ", 8.0, 0.0),
            ("suck20    ", 20.0, 0.0),
            ("suck50    ", 50.0, 0.0),
            ("suck20+t2 ", 20.0, 2.0),
        ] {
            let mut core = standard_core();
            apply_env_overrides(&mut core);
            core.plug_pair_gap = 0.35;
            core.ambient_gravity = DVec3::ZERO;
            core.ambient_charge = DVec3::ZERO;
            core.bond_suction = suction;
            core.bond_tension = tension;
            if core.ambient_charge_dir == DVec3::ZERO {
                core.ambient_charge_dir = DVec3::new(0.0, -2.0, 0.0);
            }
            core.spawn_preset("oxygen", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("oxygen");
            core.set_nucleus_dynamics(NucleusDynamics::RigidLock);
            for _ in 0..40 {
                core.solve_charge_flow();
            }
            let p_id = core.profile_id_by_name("proton").unwrap();
            let probe = core
                .spawn_particle_ex(
                    p_id,
                    DVec3::new(0.0, ds[0], 0.0),
                    DVec3::ZERO,
                    DVec3::Y,
                    0.0,
                )
                .expect("probe");
            let n_o = probe; // O nucleons are ids 0..probe
            println!(
                "\n== OH wellscan  {sname}  suction={suction} tension={tension} ==\n\
                 {:>6} | {:>6} | {:>9} | {:>8} | {:>8} | {:>8}",
                "d", "dPole", "Fy_net", "intake", "repel", "tens"
            );
            for &d in &ds {
                let pos = DVec3::new(0.0, d, 0.0);
                core.particles[probe].position = pos;
                core.particles[probe].orientation = glam::DQuat::IDENTITY;
                core.particles[probe].velocity = DVec3::ZERO;
                for _ in 0..40 {
                    core.solve_charge_flow();
                }
                core.compute_forces();
                let fy = core.particles[probe].force_accum.y;
                // nearest O nucleon + coarse term attribution.
                let mut d_pole = f64::INFINITY;
                let (mut intake, mut repel, mut tens) = (0.0f64, 0.0f64, 0.0f64);
                for o in 0..n_o {
                    let dd = (core.particles[o].position - pos).length();
                    d_pole = d_pole.min(dd);
                    let b = core.pair_force_breakdown(o, probe);
                    intake += b[5];
                    repel += b[3] + b[6] + b[7];
                    tens += b[8];
                }
                println!(
                    "{d:>6.2} | {d_pole:>6.2} | {fy:>9.3} | {intake:>8.3} | {repel:>8.3} | {tens:>8.3}"
                );
            }
        }
    }

    /// Session-33 negative control (haf.pdf free prediction): a bare
    /// FOUR-alpha stack "can't hold together" — external side-charge
    /// overwhelms the weak bare-stack channel. The current force model
    /// has no side-charge term, so today this is a MEASUREMENT of the
    /// gap, not a passing test: if quad_alpha holds as comfortably as
    /// tri_alpha, the model is missing the mechanism that sets the
    /// stack-height limit (expected to arrive with Phase C ambient
    /// surface effects). When it dissolves for the RIGHT reason, this
    /// becomes an asserting test. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_quad_alpha --nocapture`
    #[test]
    #[ignore]
    fn report_quad_alpha() {
        const STEPS: usize = 300_000;
        const SAMPLE_EVERY: usize = 10_000;
        for preset in ["tri_alpha", "quad_alpha"] {
            for k in 0..4usize {
                let mut core = standard_core();
                burn_seed_offset(&mut core, k);
                let gid = core
                    .spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("preset");
                core.set_nucleus_dynamics(
                    crate::atom_core::NucleusDynamics::RigidAlpha,
                );
                let n_alphas = core.groups[gid].alphas.len();
                let com =
                    |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                let pairs: Vec<(usize, usize)> = (0..n_alphas)
                    .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                    .collect();
                let d0: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                let mut worst = 0.0f64;
                core.running = true;
                for _ in 0..(STEPS / SAMPLE_EVERY) {
                    core.step_n(SAMPLE_EVERY);
                    for (idx, &(a, b)) in pairs.iter().enumerate() {
                        let d = (com(&core, a) - com(&core, b)).length();
                        worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                    }
                }
                println!(
                    "{preset:>10} k={k}: worst_drift={:>6.1}% KE={:.4} over {STEPS} steps",
                    worst * 100.0,
                    core.total_kinetic_energy(),
                );
            }
        }
    }

    // ── Flyby / swat robustness (session-31 round 4, work item 2/3)
    // ────────────────────────────────────────────────────────────────

    /// Every inter-alpha center spacing must stay within this fraction of
    /// its initial value for a RigidAlpha carbon nucleus to "survive" a
    /// flyby + parked-neighbor scenario. Looser than the quiet-sweep 30%
    /// (`alpha_stays_bound`) because an external close pass legitimately
    /// perturbs the stack more than internal jitter alone.
    const FLYBY_DRIFT_LIMIT: f64 = 0.35;

    /// Spawn a RigidAlpha carbon nucleus at the given (channeling,
    /// nuclear_ambient) couplings, then Phase A fly a free proton close
    /// past the stack and Phase B park a free neutron AT REST right next to
    /// it (the user's "swat" report — trivial destruction by a stray
    /// spawned particle). Returns the worst (max abs) inter-alpha
    /// center-spacing drift and whether the whole run stayed finite.
    /// `standard_core()` seeds a fixed RNG, so this is fully deterministic
    /// — a failing run can be re-driven with the same arguments to
    /// reconstruct exact diagnostics (see `print_flyby_diagnostics`).
    fn run_flyby_scenario(
        channeling: f64,
        ambient: f64,
        phase_a_steps: usize,
        phase_b_steps: usize,
        sample_every: usize,
    ) -> (f64, bool) {
        run_flyby_scenario_anatomy(
            channeling,
            ambient,
            phase_a_steps,
            phase_b_steps,
            sample_every,
            crate::atom_core::PostAnatomy::default(),
        )
    }

    /// `run_flyby_scenario` with an explicit internal-post [`PostAnatomy`]
    /// (session-32 anatomy experiment — `report_post_anatomy`).
    fn run_flyby_scenario_anatomy(
        channeling: f64,
        ambient: f64,
        phase_a_steps: usize,
        phase_b_steps: usize,
        sample_every: usize,
        anatomy: crate::atom_core::PostAnatomy,
    ) -> (f64, bool) {
        let mut core = standard_core();
        core.post_anatomy = anatomy;
        core.couplings.channeling = channeling;
        core.couplings.nuclear_ambient = ambient;
        let gid = core
            .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);

        let n_alphas = core.groups[gid].alphas.len();
        let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
        let pairs: Vec<(usize, usize)> = (0..n_alphas)
            .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
            .collect();
        let d0: Vec<f64> = pairs
            .iter()
            .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
            .collect();

        let p_id = core.profile_id_by_name("proton").unwrap();
        let n_id = core.profile_id_by_name("neutron").unwrap();

        // Phase A: a free proton passes near the stack (spawned at (8,6,0),
        // inbound velocity (-0.4,-0.3,0) — closes to a near miss then
        // continues on, mirroring the user's "spawn a stray particle
        // nearby" report). Simplest handling: let it fly through/off
        // rather than removing it (removal isn't needed for the drift
        // assertion and avoids remove_particle's swap_remove reindexing).
        core.spawn_particle(
            p_id,
            DVec3::new(8.0, 6.0, 0.0),
            DVec3::new(-0.4, -0.3, 0.0),
            DVec3::Y,
        )
        .expect("spawn flyby proton");

        let mut max_drift = 0.0f64;
        let mut finite = true;

        let n_samples_a = phase_a_steps / sample_every;
        for _ in 0..n_samples_a {
            core.step_n(sample_every);
            for (idx, &(a, b)) in pairs.iter().enumerate() {
                let d = (com(&core, a) - com(&core, b)).length();
                let drift = ((d - d0[idx]) / d0[idx]).abs();
                max_drift = max_drift.max(drift);
            }
            finite &= core
                .particles
                .iter()
                .all(|p| p.position.is_finite() && p.velocity.is_finite());
        }

        // Phase B: a free neutron parks AT REST right next to the nucleus
        // — the user's "swat" (trivial destruction by a nearby stray
        // particle with no velocity at all, the minimal perturbation).
        core.spawn_particle(n_id, DVec3::new(4.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y)
            .expect("spawn parked neutron");

        let n_samples_b = phase_b_steps / sample_every;
        for _ in 0..n_samples_b {
            core.step_n(sample_every);
            for (idx, &(a, b)) in pairs.iter().enumerate() {
                let d = (com(&core, a) - com(&core, b)).length();
                let drift = ((d - d0[idx]) / d0[idx]).abs();
                max_drift = max_drift.max(drift);
            }
            finite &= core
                .particles
                .iter()
                .all(|p| p.position.is_finite() && p.velocity.is_finite());
        }

        (max_drift, finite)
    }

    /// Deterministic re-run of `run_flyby_scenario` that prints the full
    /// per-sample drift table AND `pair_force_breakdown` for the worst
    /// (closest) member pair between the two most-drifted alphas, breaking
    /// as soon as the drift exceeds `FLYBY_DRIFT_LIMIT` (A13 honesty
    /// clause style — see `alpha_stays_bound`'s failure path).
    fn print_flyby_diagnostics(
        channeling: f64,
        ambient: f64,
        phase_a_steps: usize,
        phase_b_steps: usize,
        sample_every: usize,
    ) {
        let mut core = standard_core();
        core.couplings.channeling = channeling;
        core.couplings.nuclear_ambient = ambient;
        let gid = core
            .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);

        let n_alphas = core.groups[gid].alphas.len();
        let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
        let pairs: Vec<(usize, usize)> = (0..n_alphas)
            .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
            .collect();
        let d0: Vec<f64> = pairs
            .iter()
            .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
            .collect();
        let members_of =
            |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].members.clone();
        let closest_pair = |core: &AtomCore, ai: usize, aj: usize| -> (usize, usize) {
            let (mi, mj) = (members_of(core, ai), members_of(core, aj));
            let mut best_pair = (
                core.groups[gid].members[mi[0]],
                core.groups[gid].members[mj[0]],
            );
            let mut best_r = f64::MAX;
            for &ki in &mi {
                for &kj in &mj {
                    let pi = core.groups[gid].members[ki];
                    let pj = core.groups[gid].members[kj];
                    let r = core.pair_distance(pi, pj);
                    if r < best_r {
                        best_r = r;
                        best_pair = (pi, pj);
                    }
                }
            }
            best_pair
        };

        let p_id = core.profile_id_by_name("proton").unwrap();
        let n_id = core.profile_id_by_name("neutron").unwrap();
        core.spawn_particle(
            p_id,
            DVec3::new(8.0, 6.0, 0.0),
            DVec3::new(-0.4, -0.3, 0.0),
            DVec3::Y,
        )
        .expect("spawn flyby proton");

        println!(
            "\n{:>6} {:>7} {:>9} {:>10} {:>10} {:>10} {:>10} {:>10}",
            "step", "phase", "drift%", "r", "channel", "f_grav", "f_charge", "f_stream"
        );

        let mut step_total = 0usize;
        let mut phase = "A(flyby)";
        let mut phase_b_started = false;
        let total_samples = (phase_a_steps + phase_b_steps) / sample_every;
        let switch_at_sample = phase_a_steps / sample_every;

        for s in 0..total_samples {
            if s == switch_at_sample && !phase_b_started {
                core.spawn_particle(n_id, DVec3::new(4.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y)
                    .expect("spawn parked neutron");
                phase = "B(swat)";
                phase_b_started = true;
            }
            core.step_n(sample_every);
            step_total += sample_every;

            let mut worst_idx = 0usize;
            let mut worst_drift = 0.0f64;
            for (idx, &(a, b)) in pairs.iter().enumerate() {
                let d = (com(&core, a) - com(&core, b)).length();
                let drift = ((d - d0[idx]) / d0[idx]).abs();
                if drift > worst_drift {
                    worst_drift = drift;
                    worst_idx = idx;
                }
            }
            let (wa, wb) = pairs[worst_idx];
            let (pi, pj) = closest_pair(&core, wa, wb);
            let fb = core.pair_force_breakdown(pi, pj);
            println!(
                "{:>6} {:>7} {:>8.1}% {:>10.4} {:>10.4} {:>10.3} {:>10.3} {:>10.3}",
                step_total, phase, worst_drift * 100.0, fb[0], fb[1], fb[2], fb[3], fb[6],
            );

            let finite = core
                .particles
                .iter()
                .all(|p| p.position.is_finite() && p.velocity.is_finite());
            if !finite {
                println!("  ^ non-finite state at step {step_total} (contact term f_contact={:.3})", fb[7]);
                break;
            }
            if worst_drift > FLYBY_DRIFT_LIMIT {
                println!(
                    "  ^ drift exceeded +/-{:.0}% at step {step_total} (contact term f_contact={:.3})",
                    FLYBY_DRIFT_LIMIT * 100.0,
                    fb[7]
                );
                break;
            }
        }
    }

    /// HARD gate (session-31 round 4, work item 2): a RigidAlpha carbon
    /// nucleus must survive an external close pass AND a parked neighbor,
    /// not just internal jitter — the user's two failure reports were (a)
    /// a spontaneous late "Jenga" collapse over a long quiet run
    /// (`alpha_stays_bound` covers the quiet case) and (b) trivial
    /// destruction by spawning a stray particle nearby. This drives (b)
    /// directly: Phase A flies a free proton close past the stack, Phase B
    /// parks a free neutron at rest right next to it. Root cause (session-
    /// 31 round 4 diagnosis): the OLD channeling falloff hit a hard cliff
    /// at `r = 2·NUCLEON_PITCH` — any excursion past it regained full
    /// molecular repulsion with nothing pulling it back, so a flyby's
    /// gravitational/charge nudge (or a parked neighbor's steady pull) that
    /// displaced an alpha past the cliff was a one-way ejection. The
    /// `CHANNEL_TAIL` extension (`channeling_factor`) gives displaced
    /// alphas a restoring basin instead. Uses `Couplings::default()` so
    /// this test automatically re-validates whichever combo `alpha_stays_
    /// bound`'s sweep selects.
    #[test]
    fn nucleus_survives_flyby() {
        const PHASE_A_STEPS: usize = 15_000;
        const PHASE_B_STEPS: usize = 15_000;
        const SAMPLE_EVERY: usize = 1_000;

        let def = crate::atom_core::Couplings::default();
        let (max_drift, finite) = run_flyby_scenario(
            def.channeling,
            def.nuclear_ambient,
            PHASE_A_STEPS,
            PHASE_B_STEPS,
            SAMPLE_EVERY,
        );

        if !finite || max_drift > FLYBY_DRIFT_LIMIT {
            println!(
                "\nnucleus_survives_flyby FAILED (channeling={} ambient={}): \
                 max_drift={:.1}% finite={}",
                def.channeling,
                def.nuclear_ambient,
                max_drift * 100.0,
                finite
            );
            print_flyby_diagnostics(
                def.channeling,
                def.nuclear_ambient,
                PHASE_A_STEPS,
                PHASE_B_STEPS,
                SAMPLE_EVERY,
            );
        }
        assert!(finite, "non-finite state during flyby+swat scenario");
        assert!(
            max_drift <= FLYBY_DRIFT_LIMIT,
            "inter-alpha spacing drifted {:.1}% (> {:.0}%) during flyby+swat: \
             a displaced alpha was not recaptured — see diagnostics above",
            max_drift * 100.0,
            FLYBY_DRIFT_LIMIT * 100.0
        );
    }

    /// HARD gate: RigidAlpha neon (the "six-sided closure", nuclear.pdf —
    /// unreactive because the axial hole is surrounded by four carousel
    /// charge maxima) must stay inert once forces (not a kinematic
    /// constraint) hold it together: a probe proton dropped toward a pole
    /// from 12 r out must never be captured — never come within 2.0 of ANY
    /// member particle over the run.
    #[test]
    fn neon_is_inert() {
        let mut core = standard_core();
        let gid = core
            .spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);

        let p_id = core.profile_id_by_name("proton").unwrap();
        // Probe approaches the pole (the axial hole) from 12 r out with a
        // modest inbound velocity — enough to close the distance and
        // actually probe the closure within a test-sized step budget.
        let probe = core
            .spawn_particle(
                p_id,
                DVec3::new(0.0, 12.0, 0.0),
                DVec3::new(0.0, -0.5, 0.0),
                -DVec3::Y,
            )
            .expect("probe proton");

        let members = core.groups[gid].members.clone();
        let mut min_dist = f64::MAX;
        for _ in 0..600 {
            core.step_n(50);
            for &m in &members {
                min_dist = min_dist.min(core.pair_distance(probe, m));
            }
        }
        for p in &core.particles {
            assert!(p.position.is_finite(), "non-finite state");
        }
        assert!(
            min_dist > 2.0,
            "probe proton should never be captured by neon's six-sided closure: min_dist={min_dist}"
        );
    }

    /// REPORT scenario (session-31 addendum A11, NOT asserted beyond
    /// finiteness): does a RigidAlpha neon's carousel spin up on its own
    /// from forces alone (alphas spawn with their kinematic roll_rate as a
    /// real initial angular_velocity, §2.2 — this checks whether the FORCE
    /// model, not the initial condition, sustains/organizes a net ring
    /// rotation)? Prints net carousel ω, inter-alpha spacing drift, and
    /// total KE over the run — the deliverable is the settling behavior,
    /// not a pass/fail. Run with:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored carousel_self_organizes --nocapture`
    ///
    /// OBSERVED with the A13 binding-v2 defaults (channeling=0.9,
    /// nuclear_ambient=2.0, session-31 round 3): it does NOT self-organize.
    /// Net carousel ω decays monotonically 0.335 → 0 by ~step 25k (the
    /// initial kinematic roll is damped away, never force-sustained) and
    /// carousel spacing drifts to +705% over 40k steps (the ring migrates
    /// outward; KE decays smoothly 21.5 → 8.8, no explosion). Mechanism:
    /// binding v2 is a funnel-mouth effect — the channeling smoothstep dies
    /// at 2·NUCLEON_PITCH = 5.2, but adjacent carousel alphas sit
    /// CAROUSEL_R·√2 ≈ 6.36 apart and the core sits 4.5 from each carousel
    /// alpha center with near-perpendicular pole geometry (cos²θ ≈ 0), so
    /// neither channeling nor the weak 1/r² nuclear_ambient reaches the
    /// carousel level. Stays #[ignore]d report-style per A13 ("unless it
    /// genuinely stabilizes AND spins"). Carousel-level binding would need
    /// the edge-to-hole funnel coupling the CAROUSEL_R doc already flags as
    /// a future refinement.
    #[test]
    #[ignore]
    fn carousel_self_organizes() {
        let mut core = standard_core();
        let gid = core
            .spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon preset");
        core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::RigidAlpha);

        let axis = core.groups[gid].orientation * DVec3::Y;
        let carousel_alphas: Vec<usize> = core.groups[gid]
            .alphas
            .iter()
            .enumerate()
            .filter(|(_, a)| a.orbits_core)
            .map(|(i, _)| i)
            .collect();
        let initial_spacing: Vec<f64> = carousel_alphas
            .windows(2)
            .map(|w| {
                (core.groups[gid].alphas[w[0]].com - core.groups[gid].alphas[w[1]].com).length()
            })
            .collect();

        const STEPS: usize = 40_000;
        const SAMPLE_EVERY: usize = 200;
        println!("\n{:>8} {:>12} {:>14} {:>10}", "step", "net_omega", "spacing_drift", "total_KE");
        for s in 0..(STEPS / SAMPLE_EVERY) {
            core.step_n(SAMPLE_EVERY);
            let net_omega: f64 = carousel_alphas
                .iter()
                .map(|&ai| core.groups[gid].alphas[ai].angular_velocity.dot(axis))
                .sum::<f64>()
                / carousel_alphas.len().max(1) as f64;
            let spacing_drift: f64 = carousel_alphas
                .windows(2)
                .zip(&initial_spacing)
                .map(|(w, &d0)| {
                    let d = (core.groups[gid].alphas[w[0]].com
                        - core.groups[gid].alphas[w[1]].com)
                        .length();
                    (d - d0).abs() / d0
                })
                .fold(0.0f64, f64::max);
            let ke = core.total_kinetic_energy();
            println!(
                "{:>8} {:>12.4} {:>13.1}% {:>10.4}",
                s * SAMPLE_EVERY,
                net_omega,
                spacing_drift * 100.0,
                ke
            );
            assert!(net_omega.is_finite() && ke.is_finite(), "non-finite settling metric");
        }
    }

    /// REPORT scenario (FreeNucleon; session-31 addendum A11, NOT asserted
    /// beyond finiteness): with EVERY nucleon free — no rigid constraint at
    /// all — does a lone alpha hold together under the boosted force model,
    /// and at what `intra_nucleus_boost`? The deliverable is the tuned-
    /// constant table for a future FreeNucleon promotion, not pass/fail.
    /// Run with:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored nucleon_balance --nocapture`
    #[test]
    #[ignore]
    fn nucleon_balance() {
        const STEPS: usize = 20_000;
        println!(
            "\n{:>6} {:>16} {:>10} {:>8}",
            "boost", "max_pair_drift", "KE", "finite"
        );
        for boost in [1.0, 2.0, 3.0, 4.0, 6.0, 9.0, 12.0] {
            let mut core = standard_core();
            core.couplings.intra_nucleus_boost = boost;
            let gid = core
                .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("alpha preset");
            core.set_nucleus_dynamics(crate::atom_core::NucleusDynamics::FreeNucleon);
            let members = core.groups[gid].members.clone();
            let d0 = pair_dists(&core, &members);

            core.step_n(STEPS);

            let d1 = pair_dists(&core, &members);
            let max_rel = d0
                .iter()
                .zip(&d1)
                .map(|(a, b)| (a - b).abs() / a.max(1e-9))
                .fold(0.0f64, f64::max);
            let ke = core.total_kinetic_energy();
            let finite = core.particles.iter().all(|p| p.position.is_finite());
            println!(
                "{boost:>6.1} {:>15.1}% {ke:>10.4} {finite:>8}",
                max_rel * 100.0
            );
        }
    }
}
