//! Pure-Rust physics core for atom mode — no Godot dependencies.
//!
//! Everything simulation-related lives here so it can run headless:
//! the Godot node `AtomSim` (atom_sim.rs) is a thin wrapper, and the
//! scenario harness (atom_scenarios.rs) + `atom_lab` bin drive this
//! directly for fast physics iteration without opening the editor.

use glam::{DQuat, DVec3};
use std::f64::consts::{FRAC_PI_2, TAU};
use std::path::{Path, PathBuf};

// ── Simple PRNG (xorshift64) ─────────────────────────────────────────────

pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed.max(1),
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

// ── Emission Table ────────────────────────────────────────────────────────
// 91 bins: index 0 = pole (θ=0°), index 90 = equator (θ=90°).
// Bilateral symmetry assumed — both poles identical.

#[derive(Clone)]
pub struct EmissionTable {
    pub bins: Vec<f32>,
}

impl EmissionTable {
    pub fn sample(&self, cos_theta: f64) -> f64 {
        let theta = cos_theta.abs().acos(); // bilateral fold: 0..π/2
        let t = theta / FRAC_PI_2;
        let idx = ((t * (self.bins.len() - 1) as f64) as usize).min(self.bins.len() - 1);
        self.bins[idx] as f64
    }

    pub fn from_csv_values(csv: &[f32]) -> Self {
        assert!(csv.len() >= 180, "Need at least 180 CSV values (-90° to +89°)");
        let n = csv.len();
        let mut bins = Vec::with_capacity(91);
        for theta in 0..=90u32 {
            let north_idx = n.saturating_sub(1).min((180u32.saturating_sub(theta)) as usize);
            let south_idx = (theta as usize).min(n - 1);
            let north = csv[north_idx];
            let south = csv[south_idx];
            bins.push((north + south) / 2.0);
        }
        let max = bins.iter().copied().fold(0.0f32, f32::max);
        if max > 0.0 {
            for v in &mut bins {
                *v /= max;
            }
        }
        Self { bins }
    }

    pub fn complement(&self) -> Self {
        Self {
            bins: self.bins.iter().map(|v| 1.0 - v).collect(),
        }
    }
}

// ── Emission CDF (for weighted sampling) ─────────────────────────────────

#[derive(Clone)]
pub struct EmissionCdf {
    cdf: Vec<f64>,
}

impl EmissionCdf {
    pub fn from_bins(bins: &[f32]) -> Self {
        let n = bins.len();
        let mut cdf = Vec::with_capacity(n);
        let mut sum = 0.0f64;
        for (i, &val) in bins.iter().enumerate() {
            let theta = FRAC_PI_2 * i as f64 / (n - 1).max(1) as f64;
            sum += val as f64 * theta.sin().max(0.001);
            cdf.push(sum);
        }
        if sum > 0.0 {
            for v in &mut cdf {
                *v /= sum;
            }
        }
        Self { cdf }
    }

    pub fn sample(&self, u: f64) -> f64 {
        let u = u.clamp(0.0, 0.9999);
        let idx = self.cdf.partition_point(|&v| v < u).min(self.cdf.len() - 1);
        FRAC_PI_2 * idx as f64 / (self.cdf.len() - 1).max(1) as f64
    }
}

// ── Particle Profile ──────────────────────────────────────────────────────

#[derive(Clone)]
pub struct ParticleProfile {
    pub name: String,
    pub mass: f64,
    pub radius: f64,
    pub emission: EmissionTable,
    pub absorption: EmissionTable,
    pub emission_cdf: EmissionCdf,
}

// ── Sim Particle ──────────────────────────────────────────────────────────

pub struct SimParticle {
    pub profile_id: usize,
    pub position: DVec3,
    pub velocity: DVec3,
    pub orientation: DQuat,
    pub angular_velocity: DVec3,
    pub force_accum: DVec3,
    pub torque_accum: DVec3,
}

impl SimParticle {
    pub fn pole_axis(&self) -> DVec3 {
        self.orientation * DVec3::Y
    }
}

// ── VFX Particle (charge emission sprinkler) ─────────────────────────────

pub struct VfxParticle {
    position: DVec3,
    velocity: DVec3,
    age: f64,
    color: (f32, f32, f32),
}

pub const SOFTENING: f64 = 0.05;
pub const DAMPING: f64 = 0.9999;
pub const ANGULAR_DAMPING: f64 = 0.998;
pub const MIN_RENDER_RADIUS: f32 = 0.15;
pub const ENVELOPE_MIN_R: f64 = 0.25;
pub const CONTACT_STIFFNESS: f64 = 100.0;

// ── Force couplings ──────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Couplings {
    /// Gravity coupling (1/r², isotropic, attractive).
    pub g_q: f64,
    /// Charge coupling (1/r⁴, emission × absorption, repulsive).
    pub c_q: f64,
    /// Ambient isotropic charge pressure (pushes into charge shadows).
    pub ambient_pressure: f64,
    /// Torque coupling (sin 2θ equator-alignment).
    pub torque: f64,
    /// Equatorial vortex drag (tangential, ∝ E(θ), dead at poles).
    pub vortex: f64,
    /// Doppler velocity correction on charge.
    pub drag: f64,
    /// Polar intake (radial 1/r² sink + channeling toward pole axis, ∝ A²).
    pub intake: f64,
}

impl Default for Couplings {
    fn default() -> Self {
        Self {
            g_q: 1.0,
            c_q: 500.0,
            ambient_pressure: 0.0,
            torque: 0.1,
            vortex: 0.1,
            drag: 0.3,
            intake: 0.5,
        }
    }
}

// ── AtomCore: the whole simulation, engine-free ──────────────────────────

pub struct AtomCore {
    pub profiles: Vec<ParticleProfile>,
    pub particles: Vec<SimParticle>,
    pub couplings: Couplings,
    pub ambient_gravity: DVec3,
    pub ambient_charge: DVec3,
    pub running: bool,
    pub dt: f64,
    pub time: f64,
    rng: Rng,
    vfx_particles: Vec<VfxParticle>,
    pub vfx_enabled: bool,
}

impl Default for AtomCore {
    fn default() -> Self {
        Self::new()
    }
}

impl AtomCore {
    pub fn new() -> Self {
        Self {
            profiles: Vec::new(),
            particles: Vec::new(),
            couplings: Couplings::default(),
            ambient_gravity: DVec3::ZERO,
            ambient_charge: DVec3::ZERO,
            running: false,
            dt: 0.0005,
            time: 0.0,
            rng: Rng::new(0xDEAD_BEEF_CAFE),
            vfx_particles: Vec::new(),
            vfx_enabled: false,
        }
    }

    // ── Profile management ────────────────────────────────────────────

    pub fn register_profile(&mut self, name: &str, mass: f64, radius: f64, csv: &[f32]) -> usize {
        let emission = EmissionTable::from_csv_values(csv);
        let absorption = emission.complement();
        let emission_cdf = EmissionCdf::from_bins(&emission.bins);
        let id = self.profiles.len();
        self.profiles.push(ParticleProfile {
            name: name.to_string(),
            mass,
            radius,
            emission,
            absorption,
            emission_cdf,
        });
        id
    }

    pub fn profile_id_by_name(&self, name: &str) -> Option<usize> {
        self.profiles.iter().position(|p| p.name == name)
    }

    // ── Particle management ───────────────────────────────────────────

    pub fn spawn_particle(
        &mut self,
        profile_id: usize,
        pos: DVec3,
        vel: DVec3,
        pole_dir: DVec3,
    ) -> Option<usize> {
        if profile_id >= self.profiles.len() {
            return None;
        }
        let spin_rate = default_spin_rate(&self.profiles[profile_id].name);
        self.spawn_particle_ex(profile_id, pos, vel, pole_dir, spin_rate)
    }

    /// Spawn with an explicit signed axial spin rate (rad/s about the pole).
    /// Needed by the H₂ 8-combination matrix where spin sign matters.
    pub fn spawn_particle_ex(
        &mut self,
        profile_id: usize,
        pos: DVec3,
        vel: DVec3,
        pole_dir: DVec3,
        spin_rate: f64,
    ) -> Option<usize> {
        if profile_id >= self.profiles.len() {
            return None;
        }
        let orientation = orientation_from_pole(pole_dir);
        let pole = orientation * DVec3::Y;

        let id = self.particles.len();
        self.particles.push(SimParticle {
            profile_id,
            position: pos,
            velocity: vel,
            orientation,
            angular_velocity: pole * spin_rate,
            force_accum: DVec3::ZERO,
            torque_accum: DVec3::ZERO,
        });
        Some(id)
    }

    pub fn remove_particle(&mut self, id: usize) {
        if id < self.particles.len() {
            self.particles.swap_remove(id);
        }
    }

    pub fn clear_particles(&mut self) {
        self.particles.clear();
        self.vfx_particles.clear();
    }

    pub fn set_particle_pole(&mut self, id: usize, target: DVec3) {
        if let Some(p) = self.particles.get_mut(id) {
            if target.length_squared() > 1e-9 {
                p.orientation = orientation_from_pole(target);
            }
        }
    }

    // ── Calibration & diagnostics ─────────────────────────────────────

    /// Indices of (heaviest, lightest) particles by mass, or None if fewer
    /// than two distinct particles. Used to pick the orbit center vs orbiter.
    pub fn center_and_orbiter(&self) -> Option<(usize, usize)> {
        if self.particles.len() < 2 {
            return None;
        }
        let mut center = 0usize;
        let mut orbiter = 0usize;
        let mut m_max = f64::MIN;
        let mut m_min = f64::MAX;
        for (i, p) in self.particles.iter().enumerate() {
            let m = self.profiles[p.profile_id].mass;
            if m > m_max {
                m_max = m;
                center = i;
            }
            if m < m_min {
                m_min = m;
                orbiter = i;
            }
        }
        if center == orbiter {
            return None;
        }
        Some((center, orbiter))
    }

    /// Auto-calibrate gravity coupling for a circular polar orbit.
    /// Sets G_q so gravity alone supplies the centripetal acceleration:
    /// a_c = v_tan²/r = G_q·m_center/r² ⇒ G_q = v_tan²·r / m_center.
    /// Returns the new G_q (0.0 if it couldn't calibrate).
    pub fn auto_calibrate_polar(&mut self) -> f64 {
        let Some((center, orbiter)) = self.center_and_orbiter() else {
            return 0.0;
        };
        let m_center = self.profiles[self.particles[center].profile_id].mass;
        let d = self.particles[orbiter].position - self.particles[center].position;
        let r = d.length();
        if r < 1e-9 || m_center < 1e-12 {
            return 0.0;
        }
        let d_hat = d / r;
        let v_rel = self.particles[orbiter].velocity - self.particles[center].velocity;
        let v_tan = (v_rel - d_hat * v_rel.dot(d_hat)).length();
        let g = v_tan * v_tan * r / m_center;
        self.couplings.g_q = g;
        g
    }

    /// Orbit diagnostics: [r, v_tan, v_rad, suggested_G_q, emission_on_axis, theta_deg].
    /// theta_deg = angle from center's nearest pole to orbiter, in degrees.
    pub fn orbit_info(&self) -> [f64; 6] {
        let Some((center, orbiter)) = self.center_and_orbiter() else {
            return [0.0; 6];
        };
        let cp = &self.particles[center];
        let op = &self.particles[orbiter];
        let m_center = self.profiles[cp.profile_id].mass;
        let d = op.position - cp.position;
        let r = d.length().max(1e-9);
        let d_hat = d / r;
        let v_rel = op.velocity - cp.velocity;
        let v_rad = v_rel.dot(d_hat);
        let v_tan = (v_rel - d_hat * v_rad).length();
        let suggested_g = v_tan * v_tan * r / m_center.max(1e-12);
        let cos_theta = cp.pole_axis().dot(d_hat);
        let emission = self.profiles[cp.profile_id].emission.sample(cos_theta);
        let theta_deg = cos_theta.abs().acos().to_degrees();
        [r, v_tan, v_rad, suggested_g, emission, theta_deg]
    }

    /// Force breakdown on the orbiter from the center particle.
    /// [F_gravity, F_charge, F_intake_radial, F_channel, F_vortex_eq,
    ///  A_center, A²_center, speed].
    /// Signed scalars: positive = repulsive (away), negative = attractive.
    /// Channel is always ≥ 0 (magnitude toward the axis).
    pub fn force_breakdown(&self) -> [f64; 8] {
        let Some((ci, oi)) = self.center_and_orbiter() else {
            return [0.0; 8];
        };
        let cp = &self.particles[ci];
        let op = &self.particles[oi];
        let c_prof = &self.profiles[cp.profile_id];
        let o_prof = &self.profiles[op.profile_id];

        let d_vec = op.position - cp.position;
        let r = d_vec.length().max(SOFTENING);
        let r2 = r * r;
        let r4 = r2 * r2;
        let d_hat = d_vec / r;

        let pole = cp.pole_axis();
        let cos_theta = pole.dot(d_hat);
        let emission_c = c_prof.emission.sample(cos_theta);
        let absorption_c = c_prof.absorption.sample(cos_theta);
        let pole_o = op.pole_axis();
        let absorption_o = o_prof.absorption.sample(pole_o.dot(-d_hat));
        let ai2 = absorption_c * absorption_c;

        let cq = &self.couplings;
        let f_grav = -(cq.g_q * c_prof.mass * o_prof.mass / r2);
        let v_rel = op.velocity - cp.velocity;
        let v_rad = v_rel.dot(d_hat);
        let doppler = (1.0 - cq.drag * v_rad).clamp(0.2, 5.0);
        let f_charge = cq.c_q * c_prof.mass * o_prof.mass * emission_c * absorption_o / r4 * doppler;
        let f_intake = -(cq.intake * c_prof.mass * o_prof.mass * ai2 / r2);

        let along = d_vec.dot(pole);
        let lateral = d_vec - pole * along;
        let lat_dist = lateral.length();
        let f_channel = if lat_dist > 1e-9 {
            cq.intake * c_prof.mass * o_prof.mass * ai2 * lat_dist / (r * r2)
        } else {
            0.0
        };

        let spin = cp.angular_velocity.dot(pole);
        let sin_vec = pole.cross(d_hat);
        let sin_theta = sin_vec.length();
        let f_vortex = if sin_theta > 1e-9 {
            cq.vortex * c_prof.mass * emission_c * absorption_o * spin * sin_theta / r4
        } else {
            0.0
        };

        let speed = op.velocity.length();

        [
            f_grav, f_charge, f_intake, f_channel, f_vortex,
            absorption_c, ai2, speed,
        ]
    }

    pub fn total_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for p in &self.particles {
            let m = self.profiles[p.profile_id].mass;
            ke += 0.5 * m * p.velocity.length_squared();
        }
        ke
    }

    pub fn pair_distance(&self, a: usize, b: usize) -> f64 {
        match (self.particles.get(a), self.particles.get(b)) {
            (Some(a), Some(b)) => (a.position - b.position).length(),
            _ => -1.0,
        }
    }

    /// Contact distance (sum of render-clamped radii) for a particle pair.
    pub fn contact_distance(&self, a: usize, b: usize) -> f64 {
        let ra = self
            .particles
            .get(a)
            .map(|p| self.profiles[p.profile_id].radius.max(MIN_RENDER_RADIUS as f64))
            .unwrap_or(0.0);
        let rb = self
            .particles
            .get(b)
            .map(|p| self.profiles[p.profile_id].radius.max(MIN_RENDER_RADIUS as f64))
            .unwrap_or(0.0);
        ra + rb
    }

    // ── Simulation step (velocity Verlet) ─────────────────────────────

    pub fn step(&mut self) {
        if !self.running || self.particles.is_empty() {
            return;
        }
        let dt = self.dt;
        let n = self.particles.len();

        // 1. Forces at current state
        self.compute_forces();

        // 2. Half-step velocities + full-step positions
        for i in 0..n {
            let p = &mut self.particles[i];
            let inv_m = 1.0 / self.profiles[p.profile_id].mass;
            let a = p.force_accum * inv_m;
            let alpha = p.torque_accum;

            p.velocity += a * (dt * 0.5);
            p.angular_velocity += alpha * (dt * 0.5);
            p.position += p.velocity * dt;

            let w = p.angular_velocity;
            let w_len = w.length();
            if w_len > 1e-12 {
                let rot = DQuat::from_axis_angle(w / w_len, w_len * dt);
                p.orientation = (rot * p.orientation).normalize();
            }
        }

        // 3. Forces at new positions
        self.compute_forces();

        // 4. Complete velocity step + damping
        for i in 0..n {
            let p = &mut self.particles[i];
            let inv_m = 1.0 / self.profiles[p.profile_id].mass;
            let a = p.force_accum * inv_m;
            let alpha = p.torque_accum;

            p.velocity += a * (dt * 0.5);
            p.angular_velocity += alpha * (dt * 0.5);
            p.velocity *= DAMPING;
            // Preserve axial spin (intrinsic), only damp tumble/precession
            let pole = p.pole_axis();
            let w_axial = pole * p.angular_velocity.dot(pole);
            let w_tumble = p.angular_velocity - w_axial;
            p.angular_velocity = w_axial + w_tumble * ANGULAR_DAMPING;
        }

        self.time += dt;
    }

    pub fn step_n(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    // ── Force computation ─────────────────────────────────────────────

    fn compute_forces(&mut self) {
        let n = self.particles.len();
        for p in &mut self.particles {
            p.force_accum = DVec3::ZERO;
            p.torque_accum = DVec3::ZERO;
        }

        let cq = self.couplings;

        // Pairwise forces
        for i in 0..n {
            for j in (i + 1)..n {
                let d_vec = self.particles[j].position - self.particles[i].position;
                let r2 = d_vec.length_squared();
                let r = r2.sqrt().max(SOFTENING);
                let r2s = r * r;
                let r4s = r2s * r2s;
                let d_hat = d_vec / r;

                let pi_prof = &self.profiles[self.particles[i].profile_id];
                let pj_prof = &self.profiles[self.particles[j].profile_id];

                let pole_i = self.particles[i].pole_axis();
                let pole_j = self.particles[j].pole_axis();

                // θ_A: angle from A's pole to direction toward B
                let cos_theta_i = pole_i.dot(d_hat);
                // θ_B: angle from B's pole to direction toward A
                let cos_theta_j = pole_j.dot(-d_hat);

                let emission_i = pi_prof.emission.sample(cos_theta_i);
                let emission_j = pj_prof.emission.sample(cos_theta_j);
                let absorption_j = pj_prof.absorption.sample(cos_theta_j);
                let absorption_i = pi_prof.absorption.sample(cos_theta_i);

                // Doppler correction: approaching particles sweep through more
                // charge photons per unit time, retreating ones sweep fewer.
                // This asymmetry extracts kinetic energy on each close pass.
                let v_rel = self.particles[j].velocity - self.particles[i].velocity;
                let v_radial = v_rel.dot(d_hat);
                let doppler = (1.0 - cq.drag * v_radial).clamp(0.2, 5.0);

                // Gravity: G_q * m_i * m_j / r²  (attractive)
                let f_grav = cq.g_q * pi_prof.mass * pj_prof.mass / r2s;

                // Charge: C_q * m_emitter * m_receiver * E(θ) * R(θ) / r⁴ × doppler
                // Receiver mass models cross-section: electron intercepts 1/1836
                // the photons a proton would (Mathis's Dalton).
                let mass_prod = pi_prof.mass * pj_prof.mass;
                let f_charge_on_j =
                    cq.c_q * mass_prod * emission_i * absorption_j / r4s * doppler;
                let f_charge_on_i =
                    cq.c_q * mass_prod * emission_j * absorption_i / r4s * doppler;

                // Ambient pressure: pushes particles into charge shadows
                let shadow = (1.0 - emission_i) * (1.0 - emission_j);
                let f_ambient = cq.ambient_pressure * shadow / r2s;

                // Polar intake: the proton recycles charge through its poles,
                // creating a focused inward flow.  A(θ)² sharpens the profile
                // to a narrow polar cone matching Mathis's vortex description.
                let ai2 = absorption_i * absorption_i;
                let aj2 = absorption_j * absorption_j;

                // Radial intake: pulls toward the emitter, 1/r² (flow sink).
                let f_intake_on_j = cq.intake * pi_prof.mass * pj_prof.mass * ai2 / r2s;
                let f_intake_on_i = cq.intake * pj_prof.mass * pi_prof.mass * aj2 / r2s;

                // d_hat points from i to j.
                // Gravity pulls j toward i: along -d_hat
                // Charge pushes j away from i: along +d_hat
                // Ambient pulls j toward i: along -d_hat
                // Intake pulls j toward i: along -d_hat
                let net_on_j = (f_charge_on_j - f_grav - f_ambient - f_intake_on_j) * d_hat;
                let net_on_i = (f_grav + f_ambient + f_intake_on_i - f_charge_on_i) * d_hat;

                self.particles[j].force_accum += net_on_j;
                self.particles[i].force_accum += net_on_i;

                // Channeling: the converging intake flow funnels particles
                // toward the emitter's pole axis.  Force ∝ A²·sin(θ)/r²,
                // directed toward the nearest point on the axis.
                if cq.intake.abs() > 1e-12 {
                    // Channel j toward i's pole axis
                    let along_i = d_vec.dot(pole_i);
                    let lateral_i = d_vec - pole_i * along_i;
                    let lat_dist_i = lateral_i.length();
                    if lat_dist_i > 1e-9 {
                        let toward_axis = -lateral_i / lat_dist_i;
                        let f_chan = cq.intake * pi_prof.mass
                            * pj_prof.mass * ai2 * lat_dist_i / (r * r2s);
                        self.particles[j].force_accum += toward_axis * f_chan;
                    }
                    // Channel i toward j's pole axis
                    let d_vec_ji = -d_vec;
                    let along_j = d_vec_ji.dot(pole_j);
                    let lateral_j = d_vec_ji - pole_j * along_j;
                    let lat_dist_j = lateral_j.length();
                    if lat_dist_j > 1e-9 {
                        let toward_axis = -lateral_j / lat_dist_j;
                        let f_chan = cq.intake * pj_prof.mass
                            * pi_prof.mass * aj2 * lat_dist_j / (r * r2s);
                        self.particles[i].force_accum += toward_axis * f_chan;
                    }

                    // NOTE: intake vortex (tangential drag from spinning intake)
                    // removed — ω≈19 multiplier injects too much energy, causing
                    // slingshot escape.  The electron spawns with tangential
                    // velocity; radial intake + channeling should produce orbit
                    // without tangential drive.  Revisit with a self-limiting
                    // drag-toward-corotation model if orbits decay.
                }

                // Vortex force: spinning emission carries tangential momentum.
                // Photons leave the surface at the spin velocity, creating
                // a "sprinkler" drag that captures nearby particles into orbit.
                if cq.vortex.abs() > 1e-12 {
                    let spin_i = self.particles[i].angular_velocity.dot(pole_i);
                    let sin_vec_i = pole_i.cross(d_hat);
                    let sin_theta_i = sin_vec_i.length();
                    if sin_theta_i > 1e-9 {
                        let tangent_i = sin_vec_i / sin_theta_i;
                        let f_vort = cq.vortex * mass_prod
                            * emission_i * absorption_j
                            * spin_i * sin_theta_i / r4s;
                        self.particles[j].force_accum += tangent_i * f_vort;
                    }

                    let spin_j = self.particles[j].angular_velocity.dot(pole_j);
                    let sin_vec_j = pole_j.cross(-d_hat);
                    let sin_theta_j = sin_vec_j.length();
                    if sin_theta_j > 1e-9 {
                        let tangent_j = sin_vec_j / sin_theta_j;
                        let f_vort = cq.vortex * mass_prod
                            * emission_j * absorption_i
                            * spin_j * sin_theta_j / r4s;
                        self.particles[i].force_accum += tangent_j * f_vort;
                    }
                }

                // Contact repulsion: hard-sphere boundary at sum of radii
                let ri = pi_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let rj = pj_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let r_contact = ri + rj;
                if r < r_contact {
                    let overlap = r_contact - r;
                    let f_contact = CONTACT_STIFFNESS * overlap;
                    self.particles[j].force_accum += d_hat * f_contact;
                    self.particles[i].force_accum -= d_hat * f_contact;
                }

                // Torque: aligns equator toward charge source.
                // Uses -sin(2θ) restoring form: stable with equator facing charge.
                // torque = -2 * dot(pole, d_hat) * cross(pole, d_hat) * |f| * coupling
                let torque_j = pole_j.cross(-d_hat)
                    * (-2.0 * pole_j.dot(-d_hat) * f_charge_on_j * cq.torque);
                let torque_i = pole_i.cross(d_hat)
                    * (-2.0 * pole_i.dot(d_hat) * f_charge_on_i * cq.torque);

                self.particles[j].torque_accum += torque_j;
                self.particles[i].torque_accum += torque_i;
            }
        }

        // Ambient environmental forces
        for p in &mut self.particles {
            let m = self.profiles[p.profile_id].mass;
            p.force_accum += self.ambient_gravity * m;
            p.force_accum += self.ambient_charge * m;
        }
    }

    // ── Rendering buffers (pure data, converted by the Godot wrapper) ──

    /// Meridional cross-section ring for a particle's emission or absorption
    /// profile. World-space vertices for a closed LINE_STRIP.
    /// `plane` 0 or 1 selects two orthogonal meridional planes.
    /// `power` raises the profile value (e.g. 2.0 for A²).
    pub fn profile_ring(
        &self,
        particle_id: usize,
        use_absorption: bool,
        base_radius: f32,
        scale: f32,
        num_points: usize,
        plane: usize,
        power: f32,
    ) -> Vec<DVec3> {
        if particle_id >= self.particles.len() {
            return Vec::new();
        }
        let p = &self.particles[particle_id];
        let prof = &self.profiles[p.profile_id];
        let table = if use_absorption { &prof.absorption } else { &prof.emission };
        let center = p.position;
        let pole = p.pole_axis();
        let (right, forward) = build_frame(pole);
        let perp = if plane == 0 { right } else { forward };
        let n = num_points.max(12);
        let mut out = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let alpha = (i as f64 / n as f64) * TAU;
            let dir = pole * alpha.cos() + perp * alpha.sin();
            let cos_theta = alpha.cos();
            let val = (table.sample(cos_theta) as f32).powf(power);
            let r = base_radius + val * scale;
            out.push(center + dir * r as f64);
        }
        out
    }

    /// 16 floats per particle: 12 (3x4 transform) + 4 (RGBA color).
    /// Row-major with interleaved origin (Godot MultiMesh format).
    pub fn build_multimesh_buffer(&self) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::with_capacity(self.particles.len() * 16);
        for p in &self.particles {
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            push_transform_color(&mut buf, p, scale, profile_color(&profile.name), 1.0);
        }
        buf
    }

    /// Transform+color buffer for particles of one profile type only.
    pub fn build_multimesh_buffer_for_profile(&self, profile_id: usize) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::new();
        for p in &self.particles {
            if p.profile_id != profile_id {
                continue;
            }
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            push_transform_color(&mut buf, p, scale, (1.0, 1.0, 1.0), 1.0);
        }
        buf
    }

    /// Line buffer for pole axis indicators.
    /// 2 vertices per particle (center and pole tip), 6 floats each (pos xyz + color rgb).
    pub fn build_pole_indicator_buffer(&self) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::with_capacity(self.particles.len() * 12);
        for p in &self.particles {
            let profile = &self.profiles[p.profile_id];
            let r = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            let pos = p.position;
            let pole = p.pole_axis();
            let tip = pos + pole * (r as f64 * 1.5);

            buf.extend_from_slice(&[
                pos.x as f32, pos.y as f32, pos.z as f32,
                1.0, 0.9, 0.2,
            ]);
            buf.extend_from_slice(&[
                tip.x as f32, tip.y as f32, tip.z as f32,
                1.0, 0.4, 0.1,
            ]);
        }
        buf
    }

    /// Surface-of-revolution mesh from a profile's emission table.
    /// Packed: [vert_count, idx_count, verts...(x,y,z,nx,ny,nz,r,g,b,a), indices...]
    pub fn build_profile_mesh(
        &self,
        profile_id: usize,
        lon_segments: usize,
        lat_segments: usize,
    ) -> Vec<f32> {
        let prof = match self.profiles.get(profile_id) {
            Some(p) => p,
            None => return Vec::new(),
        };

        let lon = lon_segments.max(8);
        let lat = lat_segments.max(4);
        let rings = lat * 2 + 1;
        let verts_per_ring = lon + 1;
        let total_verts = rings * verts_per_ring;
        let quad_count = (rings - 1) * lon;
        let total_indices = quad_count * 6;

        let floats_per_vert = 10; // pos(3) + normal(3) + color(4)
        let mut buf: Vec<f32> =
            Vec::with_capacity(2 + total_verts * floats_per_vert + total_indices);

        buf.push(total_verts as f32);
        buf.push(total_indices as f32);

        let (type_r, type_g, type_b) = profile_color(&prof.name);

        for lat_idx in 0..rings {
            let theta = std::f64::consts::PI * lat_idx as f64 / (rings - 1).max(1) as f64;
            let cos_theta = theta.cos();
            let sin_theta = theta.sin();

            let emission = prof.emission.sample(cos_theta);
            // Lathe approach: fixed height (y = cos θ), emission modulates width only
            let r_cross = sin_theta * (ENVELOPE_MIN_R + emission * (1.0 - ENVELOPE_MIN_R));
            let y_pos = cos_theta;

            let intensity = (0.3f32 + emission as f32 * 0.7).min(1.0);
            let alpha = 0.25 + emission as f32 * 0.25;

            for lon_idx in 0..=lon {
                let phi = TAU * lon_idx as f64 / lon as f64;
                let x = (r_cross * phi.cos()) as f32;
                let y = y_pos as f32;
                let z = (r_cross * phi.sin()) as f32;

                buf.extend_from_slice(&[x, y, z]);

                let len = (x * x + y * y + z * z).sqrt();
                if len > 1e-6 {
                    buf.extend_from_slice(&[x / len, y / len, z / len]);
                } else {
                    buf.extend_from_slice(&[0.0, if lat_idx == 0 { 1.0 } else { -1.0 }, 0.0]);
                }

                buf.extend_from_slice(&[
                    type_r * intensity,
                    type_g * intensity,
                    type_b * intensity,
                    alpha,
                ]);
            }
        }

        for lat_idx in 0..(rings - 1) {
            for lon_idx in 0..lon {
                let tl = (lat_idx * verts_per_ring + lon_idx) as f32;
                let tr = tl + 1.0;
                let bl = ((lat_idx + 1) * verts_per_ring + lon_idx) as f32;
                let br = bl + 1.0;
                buf.extend_from_slice(&[tl, bl, tr, tr, bl, br]);
            }
        }

        buf
    }

    pub fn count_particles_with_profile(&self, profile_id: usize) -> usize {
        self.particles
            .iter()
            .filter(|p| p.profile_id == profile_id)
            .count()
    }

    // ── VFX: charge emission sprinkler ───────────────────────────────────

    /// Advance VFX particles and return the MultiMesh buffer
    /// (12 transform + 4 color floats per particle).
    pub fn advance_vfx(
        &mut self,
        delta: f64,
        emit_per_particle: usize,
        speed: f64,
        lifetime: f64,
    ) -> Vec<f32> {
        if !self.vfx_enabled || self.particles.is_empty() {
            self.vfx_particles.clear();
            return Vec::new();
        }

        let max_pool = 2048usize;

        // Age and move existing
        for vp in &mut self.vfx_particles {
            vp.age += delta;
            vp.position += vp.velocity * delta;
        }
        self.vfx_particles.retain(|vp| vp.age < lifetime);

        // Collect per-particle data to avoid borrow conflicts
        let emit_info: Vec<(DVec3, DVec3, DVec3, DVec3, f64, (f32, f32, f32), usize)> = self
            .particles
            .iter()
            .map(|p| {
                let prof = &self.profiles[p.profile_id];
                let pole = p.pole_axis();
                let (right, forward) = build_frame(pole);
                let radius = prof.radius.max(MIN_RENDER_RADIUS as f64);
                let color = profile_color(&prof.name);
                (p.position, pole, right, forward, radius, color, p.profile_id)
            })
            .collect();

        for (pos, pole, right, forward, radius, color, pid) in &emit_info {
            for _ in 0..emit_per_particle {
                if self.vfx_particles.len() >= max_pool {
                    break;
                }
                let theta = self.profiles[*pid].emission_cdf.sample(self.rng.next_f64());
                let north = self.rng.next_f64() > 0.5;
                let phi = self.rng.next_f64() * TAU;

                let cos_t = theta.cos();
                let sin_t = theta.sin();
                let local_y = if north { cos_t } else { -cos_t };
                let dir =
                    *right * (sin_t * phi.cos()) + *pole * local_y + *forward * (sin_t * phi.sin());
                let spawn_pos = *pos + dir * *radius;

                self.vfx_particles.push(VfxParticle {
                    position: spawn_pos,
                    velocity: dir * speed,
                    age: 0.0,
                    color: *color,
                });
            }
        }

        // Build MultiMesh buffer: 12 (transform) + 4 (color) per particle
        let n = self.vfx_particles.len();
        let mut buf = Vec::with_capacity(n * 16);
        let scale = 0.03f32;

        for vp in &self.vfx_particles {
            let fade = ((1.0 - vp.age / lifetime) as f32).max(0.0);
            buf.extend_from_slice(&[
                scale, 0.0, 0.0, vp.position.x as f32,
                0.0, scale, 0.0, vp.position.y as f32,
                0.0, 0.0, scale, vp.position.z as f32,
                vp.color.0, vp.color.1, vp.color.2, fade * 0.8,
            ]);
        }

        buf
    }

    pub fn vfx_count(&self) -> usize {
        self.vfx_particles.len()
    }

    pub fn set_vfx_enabled(&mut self, enabled: bool) {
        self.vfx_enabled = enabled;
        if !enabled {
            self.vfx_particles.clear();
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────

/// Quaternion rotating +Y onto `target` (falls back to +Y for zero vectors).
pub fn orientation_from_pole(target: DVec3) -> DQuat {
    let target = if target.length_squared() > 1e-9 {
        target.normalize()
    } else {
        DVec3::Y
    };
    if (target - DVec3::Y).length_squared() < 1e-9 {
        DQuat::IDENTITY
    } else if (target + DVec3::Y).length_squared() < 1e-9 {
        DQuat::from_rotation_z(std::f64::consts::PI)
    } else {
        DQuat::from_rotation_arc(DVec3::Y, target)
    }
}

pub fn build_frame(pole: DVec3) -> (DVec3, DVec3) {
    let right = if pole.dot(DVec3::X).abs() < 0.9 {
        pole.cross(DVec3::X).normalize()
    } else {
        pole.cross(DVec3::Z).normalize()
    };
    let forward = right.cross(pole);
    (right, forward)
}

pub fn default_spin_rate(name: &str) -> f64 {
    match name {
        "proton" => TAU * 3.0,
        "neutron" => TAU * 3.0,
        "electron" => TAU * 5.0,
        _ => TAU,
    }
}

pub fn profile_color(name: &str) -> (f32, f32, f32) {
    match name {
        "proton" => (0.92, 0.30, 0.20),
        "neutron" => (0.35, 0.50, 0.92),
        "electron" => (0.20, 0.90, 0.35),
        _ => (0.7, 0.7, 0.7),
    }
}

fn push_transform_color(
    buf: &mut Vec<f32>,
    p: &SimParticle,
    scale: f32,
    color: (f32, f32, f32),
    alpha: f32,
) {
    let pos = p.position;
    let x = p.orientation * DVec3::X;
    let y = p.orientation * DVec3::Y;
    let z = p.orientation * DVec3::Z;

    buf.extend_from_slice(&[
        x.x as f32 * scale, x.y as f32 * scale, x.z as f32 * scale, pos.x as f32,
        y.x as f32 * scale, y.y as f32 * scale, y.z as f32 * scale, pos.y as f32,
        z.x as f32 * scale, z.y as f32 * scale, z.z as f32 * scale, pos.z as f32,
    ]);
    buf.extend_from_slice(&[color.0, color.1, color.2, alpha]);
}

// ── Headless CSV loading (tests + atom_lab bin) ──────────────────────────

/// Load a `degree,percentage` histogram CSV, returning the value column.
/// Panics with a clear message on malformed files (test/dev tool, not UI path).
pub fn load_histogram_csv(path: &Path) -> Vec<f32> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", path.display(), e));
    let mut values = Vec::with_capacity(181);
    for line in text.lines().skip(1) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let value = line
            .split(',')
            .nth(1)
            .unwrap_or_else(|| panic!("malformed CSV line in {}: {}", path.display(), line));
        values.push(value.trim().parse::<f32>().unwrap_or_else(|e| {
            panic!("bad value in {}: {} ({})", path.display(), line, e)
        }));
    }
    values
}

/// The godot/config directory, resolved relative to the Cargo manifest so
/// tests and bins work regardless of the working directory.
pub fn config_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../godot/config")
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_flat_csv() -> Vec<f32> {
        vec![1.0; 181]
    }

    fn make_proton_like_csv() -> Vec<f32> {
        let mut csv = vec![0.0; 181];
        // Equator (index 90) = bright, poles (0 and 180) = dim
        for i in 0..181 {
            let deg = i as f64 - 90.0;
            let frac = 1.0 - (deg.abs() / 90.0);
            csv[i] = frac as f32;
        }
        csv
    }

    #[test]
    fn emission_table_bilateral_symmetry() {
        let csv = make_proton_like_csv();
        let table = EmissionTable::from_csv_values(&csv);
        let pole_val = table.sample(1.0);
        let anti_pole = table.sample(-1.0);
        assert!(
            (pole_val - anti_pole).abs() < 1e-6,
            "Bilateral: pole={} vs anti-pole={}",
            pole_val,
            anti_pole
        );
    }

    #[test]
    fn emission_table_equator_vs_pole() {
        let csv = make_proton_like_csv();
        let table = EmissionTable::from_csv_values(&csv);
        let equator = table.sample(0.0);
        let pole = table.sample(1.0);
        assert!(
            equator > pole,
            "Equator should be brighter: eq={} pole={}",
            equator,
            pole
        );
    }

    #[test]
    fn complement_inverse() {
        let csv = make_proton_like_csv();
        let emission = EmissionTable::from_csv_values(&csv);
        let absorption = emission.complement();
        let e = emission.sample(0.5);
        let a = absorption.sample(0.5);
        assert!(
            (e + a - 1.0).abs() < 1e-6,
            "emission + absorption should be 1.0: {} + {} = {}",
            e,
            a,
            e + a
        );
    }

    #[test]
    fn flat_emission_gives_uniform_one() {
        let table = EmissionTable::from_csv_values(&make_flat_csv());
        for i in 0..=10 {
            let cos_t = (i as f64 - 5.0) / 5.0;
            let val = table.sample(cos_t);
            assert!(
                (val - 1.0).abs() < 1e-5,
                "Flat CSV should give 1.0 everywhere, got {} at cos_theta={}",
                val,
                cos_t
            );
        }
    }

    #[test]
    fn histogram_csvs_load_and_register() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        for name in ["proton", "neutron", "electron"] {
            let csv = load_histogram_csv(&dir.join(format!("histogram_{name}.csv")));
            assert!(csv.len() >= 180, "{name} CSV too short: {}", csv.len());
            core.register_profile(name, 1.0, 1.0, &csv);
        }
        assert_eq!(core.profiles.len(), 3);
        // Proton emission must peak at the equator, vanish near the pole
        let proton = &core.profiles[0];
        assert!(proton.emission.sample(0.0) > 0.5);
        assert!(proton.emission.sample(1.0) < 0.1);
    }
}
