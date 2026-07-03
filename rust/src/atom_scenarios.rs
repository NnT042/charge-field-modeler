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

    /// M3 gates: alpha rigidity + conservation, electron capture on alpha,
    /// helium stability.
    #[test]
    #[ignore = "M3: needs rigid groups + presets"]
    fn alpha_holds_and_conserves() {
        unimplemented!("filled in at M3");
    }

    #[test]
    #[ignore = "M3: needs rigid groups + presets"]
    fn alpha_captures_electron() {
        unimplemented!("filled in at M3");
    }

    #[test]
    #[ignore = "M3: needs rigid groups + presets"]
    fn helium_stable() {
        unimplemented!("filled in at M3");
    }
}
