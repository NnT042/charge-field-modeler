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
        for name in ["carbon", "nitrogen", "oxygen", "neon", "argon"] {
            let mut core = standard_core();
            let gid = core
                .spawn_preset(name, DVec3::ZERO, DVec3::new(0.1, 0.0, 0.05), DVec3::Y)
                .unwrap_or_else(|| panic!("{name} preset missing"));
            let members = core.groups[gid].members.clone();
            let expected = match name {
                "carbon" => 12,
                "nitrogen" => 14,
                "oxygen" => 16,
                "neon" => 20,   // center alpha + 4 carousel alphas
                "argon" => 36,  // 9 alphas: axial line of 5 + 4 carousel
                _ => unreachable!(),
            };
            assert_eq!(members.len(), expected, "{name} constituent count");

            // Polar plug geometry (session-29), rewritten in AlphaUnit
            // terms (session-31 addendum A7): plug alphas are single-
            // member, orbits_core == true; plug PROTONS sit edge-on to
            // the stack (rest_axis·Y ≈ 0, disc feeds the hole), plug
            // NEUTRONS keep their pole on the stack axis (rest_axis ≈
            // +Y, graphene.pdf: neutrons channel pole-to-pole); oxygen's
            // plugs come as side-by-side proton+neutron pairs.
            let profile_name_of = |core: &AtomCore, gid: usize, ai: usize| -> String {
                let g = &core.groups[gid];
                let k = g.alphas[ai].members[0];
                core.profiles[core.particles[g.members[k]].profile_id]
                    .name
                    .clone()
            };
            let plug_alphas: Vec<usize> = core.groups[gid]
                .alphas
                .iter()
                .enumerate()
                .filter(|(_, a)| a.members.len() == 1)
                .map(|(i, _)| i)
                .collect();
            for &ai in &plug_alphas {
                let a = &core.groups[gid].alphas[ai];
                assert!(a.orbits_core, "plug alpha should ride the carousel");
                match profile_name_of(&core, gid, ai).as_str() {
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
            match name {
                "nitrogen" => assert_eq!(plug_alphas.len(), 2, "nitrogen has 2 plug alphas"),
                "oxygen" => {
                    assert_eq!(plug_alphas.len(), 4, "oxygen has 4 plug alphas");
                    let mut by_y: Vec<(f64, usize)> = plug_alphas
                        .iter()
                        .map(|&ai| (core.groups[gid].alphas[ai].rest_center.y, ai))
                        .collect();
                    by_y.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
                    for pair in by_y.chunks(2) {
                        let (ya, ai_a) = pair[0];
                        let (yb, ai_b) = pair[1];
                        assert!(
                            (ya - yb).abs() < 1e-9,
                            "pair must sit at the same height (side by side)"
                        );
                        let (proton_ai, neutron_ai) =
                            if profile_name_of(&core, gid, ai_a) == "proton" {
                                (ai_a, ai_b)
                            } else {
                                (ai_b, ai_a)
                            };
                        let g = &core.groups[gid];
                        let ca = g.alphas[proton_ai].rest_center;
                        let cb = g.alphas[neutron_ai].rest_center;
                        assert!(
                            (ca - cb).length() > 0.5,
                            "pair members must sit beside each other, not overlap"
                        );
                        let to_partner = (cb - ca).normalize();
                        assert!(
                            g.alphas[proton_ai].rest_axis.dot(to_partner) > 0.99,
                            "plug proton rest_axis should point at its paired neutron"
                        );
                    }
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
    /// `compute_forces`. This test sweeps both dimensions on the same
    /// carbon nucleus (3 axially-stacked alphas, all `orbits_core == false`
    /// — no carousel confound) in RigidAlpha. Pass = SOME combo keeps every
    /// inter-alpha center-pair distance within ±30% of its initial value
    /// over 20k steps (the full sweep runs in ~12s release — inside the
    /// ~30s `cargo test` gate now that this is a hard, non-ignored test).
    #[test]
    fn alpha_stays_bound() {
        const STEPS: usize = 20_000;
        const SAMPLE_EVERY: usize = 50;
        let channelings = [0.7, 0.9, 1.0];
        let ambients = [0.0, 2.0, 5.0, 10.0, 20.0];

        println!(
            "\n{:>5} {:>6} {:>10} {:>9} {:>9} {:>10}",
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
                    .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
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
                    "{channeling:>5.1} {ambient:>6.1} {:>9.1}% {:>9.3} {:>9.3} {:>10}",
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
                .spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
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

        // Session-31 round 3 RESULT: 9 of 15 combos bind. Winner
        // channeling=0.9, nuclear_ambient=2.0 at 2.5% max drift — promoted
        // to `Couplings::default` (see its doc for the full table's edge
        // analysis: ambient=0 ejects on residual repulsion; channeling=1.0
        // with ambient>=5 cancels repulsion completely, collapses into the
        // disc-aware contact wall, and the never-attenuated contact spring
        // ejects at 1e4-1e6% drift). Assert the shipped DEFAULTS are a
        // passing combo, so any future default change must re-earn this.
        println!(
            "\n{} combo(s) bind within +/-30%; best: channeling={bc} ambient={ba} drift={:.1}%",
            passing.len(),
            bd * 100.0
        );
        let def = crate::atom_core::Couplings::default();
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
            "Couplings::default (channeling={}, nuclear_ambient={}) no longer binds carbon: \
             max drift {:.1}% > 30%",
            def.channeling,
            def.nuclear_ambient,
            def_result.2 * 100.0
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
