//! Field Calibration Mode — the freed-DOF rigid body (CM-0).
//!
//! See `docs/FIELD_CALIBRATION_MODE.md`. Spin mode froze the particle in place;
//! here a baked-loop hitbox ([`crate::hitbox::BakedLoop`]) becomes a rigid body
//! that can translate, tumble, and change its outer-spin rate under collisions.
//!
//! The single-photon impulse splits three ways — the whole point of the mode:
//!   1. **Linear** — `impulse / mass` → drift.
//!   2. **Transverse torque** — the part of `r × J` perpendicular to the pole
//!      tumbles the whole body (pole wander).
//!   3. **Pole-axis torque** — the part along the pole pumps `outer_spin` (the
//!      swing rate of level `L+1`), routed through a lumped `spin_coupling`
//!      coefficient that the mode MEASURES (decision 1), never hardcodes.
//!
//! Moments of inertia are computed from the baked geometry, so the enormous
//! amplitude of a baryon loop makes each hit's Δ(outer_spin) tiny *for free* —
//! the "trillions of collisions to move it" physics is geometric, not a magic
//! constant. `outer_spin` saturates at ±1 (= c); reaching it flags a
//! transmutation (the spin stack would unlock the next level — e.g. a neutron's
//! L12 hitting c unlocks L13 and it becomes a proton).

use glam::{DQuat, DVec3};

use crate::hitbox::{bake_loop, recommended_samples, BakedLoop};
use crate::types::level_amplitude;

pub struct CalibrationParticle {
    pub hitbox: BakedLoop,
    /// Center of mass, natural units.
    pub position: DVec3,
    /// Body pose. The world pole is `orientation * hitbox.swing_axis`.
    pub orientation: DQuat,
    /// Accumulated linear drift velocity.
    pub lin_velocity: DVec3,
    /// Whole-body tumble, world frame, rad/time.
    pub ang_velocity: DVec3,
    /// Swing-level angular velocity in [-1, 1]; ±1 = c = saturation.
    pub outer_spin: f64,
    /// Accumulated swing angle (rad) fed to `BakedLoop::swept_point`.
    pub swing_phase: f64,
    /// Particle mass (natural units). Scales both linear and rotational inertia.
    pub mass: f64,
    /// Lumped coefficient routing pole-axis torque into `outer_spin`. MEASURED
    /// by the calibration mode (decision 1); 1.0 is a scaffold placeholder.
    pub spin_coupling: f64,
    /// Set once `|outer_spin|` reaches c — the next spin level would unlock.
    pub transmuted: bool,
    /// orbit_radius of the swing level, used to convert `outer_spin` (fraction
    /// of c) into an angular rate for `swing_phase`.
    swing_orbit_radius: f64,
    /// Moment of inertia about the pole axis (drives `outer_spin` response).
    i_spin: f64,
    /// Moment of inertia about a transverse COM axis (drives tumble response).
    i_transverse: f64,
}

impl CalibrationParticle {
    pub fn new(hitbox: BakedLoop, mass: f64) -> Self {
        let mass = mass.max(1e-12);
        // Precession swings rotate in place (orbit_radius = amplitude); orbital
        // swings ride a string (orbit_radius = amplitude/2), per SpinLevel.
        let swing_orbit_radius = if hitbox.swing_is_precession {
            level_amplitude(hitbox.swing_level)
        } else {
            level_amplitude(hitbox.swing_level) / 2.0
        }
        .max(1e-9);

        let (i_spin, i_transverse) = Self::compute_inertia(&hitbox, mass);

        Self {
            hitbox,
            position: DVec3::ZERO,
            orientation: DQuat::IDENTITY,
            lin_velocity: DVec3::ZERO,
            ang_velocity: DVec3::ZERO,
            outer_spin: 0.0,
            swing_phase: 0.0,
            mass,
            spin_coupling: 1.0,
            transmuted: false,
            swing_orbit_radius,
            i_spin,
            i_transverse,
        }
    }

    /// Convenience: a proton (L12 loop, L13 precession swing).
    pub fn proton(mass: f64) -> Self {
        Self::new(bake_loop(12, recommended_samples(12)), mass)
    }

    /// Convenience: a neutron (L11 loop, non-relativistic L12 orbital swing).
    pub fn neutron(mass: f64) -> Self {
        Self::new(bake_loop(11, recommended_samples(11)), mass)
    }

    /// Convenience: an electron (L8 loop, L9 precession swing).
    pub fn electron(mass: f64) -> Self {
        Self::new(bake_loop(8, recommended_samples(8)), mass)
    }

    /// Inertia about the pole axis and about a transverse COM axis, computed
    /// from the swept mass distribution (loop points + kite offset). The
    /// transverse value is an axisymmetric estimate (`½·perp² + along²` per unit
    /// mass) — adequate for CM-0; refine to a full tensor in CM-3 if needed.
    fn compute_inertia(h: &BakedLoop, mass: f64) -> (f64, f64) {
        let n = h.points.len().max(1);
        let m_per = mass / n as f64;
        let axis = h.swing_axis.normalize_or_zero();
        let mut i_spin = 0.0;
        let mut i_trans = 0.0;
        for &p in &h.points {
            let q = p + h.swing_offset; // position relative to COM in the pole frame
            let along = q.dot(axis);
            let perp2 = (q - axis * along).length_squared();
            i_spin += m_per * perp2;
            i_trans += m_per * (0.5 * perp2 + along * along);
        }
        (i_spin.max(1e-12), i_trans.max(1e-12))
    }

    /// World-space pole (spin) axis.
    pub fn pole_axis(&self) -> DVec3 {
        (self.orientation * self.hitbox.swing_axis).normalize_or_zero()
    }

    /// The live collision geometry: baked loop, swung to `swing_phase`, then
    /// posed by the body's orientation and position.
    pub fn world_points(&self) -> Vec<DVec3> {
        self.hitbox
            .swept_points(self.swing_phase)
            .into_iter()
            .map(|p| self.position + self.orientation * p)
            .collect()
    }

    /// Apply a single photon impulse `J` delivered at world point `contact`.
    /// Splits into linear drift, transverse tumble, and pole-axis spin pump.
    pub fn apply_impulse(&mut self, contact: DVec3, impulse: DVec3) {
        self.lin_velocity += impulse / self.mass;

        let r = contact - self.position;
        let tau = r.cross(impulse);
        let pole = self.pole_axis();
        let tau_spin = tau.dot(pole);
        let tau_perp = tau - pole * tau_spin;

        self.outer_spin += self.spin_coupling * tau_spin / self.i_spin;
        self.ang_velocity += tau_perp / self.i_transverse;

        if self.outer_spin.abs() >= 1.0 {
            self.outer_spin = self.outer_spin.clamp(-1.0, 1.0);
            self.transmuted = true;
        }
    }

    /// World velocity of a point rigidly attached to the body: drift + tumble +
    /// the outer swing. `time_scale` converts `outer_spin` (fraction of c) into a
    /// swing angular rate, matching `integrate`.
    pub fn contact_velocity(&self, point: DVec3, time_scale: f64) -> DVec3 {
        let r = point - self.position;
        let omega_swing = self.pole_axis() * (self.outer_spin * time_scale / self.swing_orbit_radius);
        self.lin_velocity + (self.ang_velocity + omega_swing).cross(r)
    }

    /// Unit direction the outer-spin surface is moving at `point` (independent of
    /// swing rate magnitude). ZERO if the particle isn't spinning or the point
    /// lies on the pole axis — then the gear catch is neutral.
    fn swing_tangent(&self, point: DVec3) -> DVec3 {
        if self.outer_spin.abs() < 1e-12 {
            return DVec3::ZERO;
        }
        let axis = self.pole_axis();
        let r = point - self.position;
        let rp = r - axis * r.dot(axis);
        (axis.cross(rp) * self.outer_spin.signum()).normalize_or_zero()
    }

    /// Apply one field photon under the gear rule (model B, see
    /// docs/FIELD_CALIBRATION_MODE.md). Momentum transfers linearly + as tumble
    /// (Newtonian); the pole-axis spin pump is `spin_gain · f · chirality / i_spin`
    /// where the catch factor `f = (1 − dir·t̂)/2` is 1 for a photon opposing the
    /// surface (gears catch), 0 for one co-moving with it (slips), and `chirality`
    /// (+1 photon / −1 antiphoton) sets augment vs cancel. Inertia resists, so a
    /// heavy baryon barely moves per hit.
    ///
    /// `dir` is the photon travel direction; `chirality` ∈ {+1, −1}.
    pub fn apply_photon(
        &mut self,
        contact: DVec3,
        dir: DVec3,
        chirality: f64,
        momentum: f64,
        spin_gain: f64,
    ) {
        let dir = dir.normalize_or_zero();
        let j = dir * momentum;

        // Newtonian: linear drift + transverse tumble (chirality-independent).
        self.lin_velocity += j / self.mass;
        let r = contact - self.position;
        let tau = r.cross(j);
        let pole = self.pole_axis();
        let tau_perp = tau - pole * tau.dot(pole);
        self.ang_velocity += tau_perp / self.i_transverse;

        // Gear rule: pole-axis spin pump, catch strength × chirality, inertia-resisted.
        let t = self.swing_tangent(contact);
        let f = if t.length_squared() > 1e-18 {
            ((1.0 - dir.dot(t)) * 0.5).clamp(0.0, 1.0)
        } else {
            0.5
        };
        self.outer_spin += spin_gain * f * chirality / self.i_spin;
        if self.outer_spin.abs() >= 1.0 {
            self.outer_spin = self.outer_spin.clamp(-1.0, 1.0);
            self.transmuted = true;
        }
    }

    /// DIAGNOSTIC (measurement, not a rule): the presented-orientation
    /// distribution *as the field actually samples it*. Each (loop point, swing
    /// phase) is weighted by its gear catch factor `f` against the photon
    /// direction(s) at the current `outer_spin`, then binned by cos(presented
    /// spin axis · pole). `field_dir = Some(d)` is a directional field; `None`
    /// averages over `dir_samples` isotropic directions. Averages over
    /// `swing_samples` swing phases. Returns `bins` fractions summing to 1.
    ///
    /// This answers "does the hit-weighting break the loop's orientation
    /// symmetry as the disc spins up?" without touching any collision rule.
    pub fn orientation_participation(
        &self,
        field_dir: Option<DVec3>,
        bins: usize,
        swing_samples: usize,
        dir_samples: usize,
    ) -> Vec<f64> {
        let bins = bins.max(1);
        let mut h = vec![0.0; bins];
        let n = self.hitbox.points.len();
        if n == 0 {
            return h;
        }
        let pole = self.pole_axis();
        let swing_samples = swing_samples.max(1);

        let dirs: Vec<DVec3> = match field_dir {
            Some(d) => vec![d.normalize_or_zero()],
            None => {
                let m = dir_samples.max(1);
                (0..m)
                    .map(|k| {
                        // Deterministic Fibonacci-sphere sampling.
                        let z = 1.0 - 2.0 * (k as f64 + 0.5) / m as f64;
                        let r = (1.0 - z * z).max(0.0).sqrt();
                        let phi = k as f64 * 2.399_963_229_728_653;
                        DVec3::new(r * phi.cos(), r * phi.sin(), z)
                    })
                    .collect()
            }
        };

        // Physical time_scale: surface reaches c at outer_spin = 1.
        let omega = pole * (self.outer_spin / self.swing_orbit_radius);
        let inv_dirs = 1.0 / dirs.len() as f64;

        for si in 0..swing_samples {
            let phi = si as f64 / swing_samples as f64 * std::f64::consts::TAU;
            let srot = DQuat::from_axis_angle(self.hitbox.swing_axis, phi);
            for i in 0..n {
                let local = srot * (self.hitbox.points[i] + self.hitbox.swing_offset);
                let world_p = self.position + self.orientation * local;
                let vhat = omega.cross(world_p - self.position).normalize_or_zero();
                let a = (self.orientation * (srot * self.hitbox.orientations[i])) * DVec3::Y;
                let c = a.dot(pole).clamp(-1.0, 1.0);
                let idx = ((((c + 1.0) * 0.5) * bins as f64).floor() as usize).min(bins - 1);
                let mut w = 0.0;
                for d in &dirs {
                    w += if vhat.length_squared() > 1e-18 {
                        ((1.0 - d.dot(vhat)) * 0.5).clamp(0.0, 1.0)
                    } else {
                        0.5
                    };
                }
                h[idx] += w * inv_dirs;
            }
        }

        let sum: f64 = h.iter().sum();
        if sum > 0.0 {
            for x in &mut h {
                *x /= sum;
            }
        }
        h
    }

    /// DIAGNOSTIC (measurement, not a rule): net pole-axis tangential torque per
    /// photon from isotropic bombardment, catch-weighted by relative speed
    /// (gears catch ∝ |v_photon − v_surface|). This is the VELOCITY channel: as
    /// the swing surface approaches c, head-on photons (opposing the surface)
    /// have relative speed ~c+s and co-moving ones ~c−s, so the catch is
    /// asymmetric in a way that does NOT swing-average out. A negative return
    /// growing with `outer_spin` is a spin-DOWN drag — a candidate self-limiter
    /// that owes nothing to chirality. No collision rule is changed here.
    pub fn swing_drag_torque(&self, dir_samples: usize, swing_samples: usize) -> f64 {
        let n = self.hitbox.points.len();
        if n == 0 {
            return 0.0;
        }
        let pole = self.pole_axis();
        let omega = pole * (self.outer_spin / self.swing_orbit_radius);
        let swing_samples = swing_samples.max(1);
        let m = dir_samples.max(1);
        let dirs: Vec<DVec3> = (0..m)
            .map(|k| {
                let z = 1.0 - 2.0 * (k as f64 + 0.5) / m as f64;
                let r = (1.0 - z * z).max(0.0).sqrt();
                let phi = k as f64 * 2.399_963_229_728_653;
                DVec3::new(r * phi.cos(), r * phi.sin(), z)
            })
            .collect();

        let mut torque = 0.0;
        let mut count = 0.0;
        for si in 0..swing_samples {
            let phi = si as f64 / swing_samples as f64 * std::f64::consts::TAU;
            let srot = DQuat::from_axis_angle(self.hitbox.swing_axis, phi);
            for i in 0..n {
                let world_p =
                    self.position + self.orientation * (srot * (self.hitbox.points[i] + self.hitbox.swing_offset));
                let r = world_p - self.position;
                let v = omega.cross(r); // surface velocity
                let vt = v.length();
                if vt < 1e-12 {
                    // still surface: isotropic sum of tangential momentum is zero.
                    count += m as f64;
                    continue;
                }
                let that = v / vt; // surface tangential direction
                let lever = (r - pole * r.dot(pole)).length();
                for d in &dirs {
                    let catch = (*d - v).length(); // relative speed (c=1)
                    torque += catch * d.dot(that) * lever; // +along surface = spin-up
                    count += 1.0;
                }
            }
        }
        if count > 0.0 {
            torque / count
        } else {
            0.0
        }
    }

    /// DIAGNOSTIC (measurement, not a rule): net whole-body torque from a
    /// DIRECTIONAL field along `field_dir`, catch-weighted by relative speed and
    /// swing-averaged. The component along the pole is the spin drag; the
    /// component perpendicular to the pole is the ALIGNMENT torque — does the
    /// field turn the pole? Tilt the particle (set `orientation`) and read the
    /// component about the tilt axis to see if there's a restoring alignment.
    pub fn directional_torque(&self, field_dir: DVec3, swing_samples: usize) -> DVec3 {
        let n = self.hitbox.points.len();
        if n == 0 {
            return DVec3::ZERO;
        }
        let d = field_dir.normalize_or_zero();
        let pole = self.pole_axis();
        let omega = pole * (self.outer_spin / self.swing_orbit_radius);
        let swing_samples = swing_samples.max(1);
        let mut tau = DVec3::ZERO;
        let mut count = 0.0;
        for si in 0..swing_samples {
            let phi = si as f64 / swing_samples as f64 * std::f64::consts::TAU;
            let srot = DQuat::from_axis_angle(self.hitbox.swing_axis, phi);
            for i in 0..n {
                let world_p =
                    self.position + self.orientation * (srot * (self.hitbox.points[i] + self.hitbox.swing_offset));
                let r = world_p - self.position;
                let catch = (d - omega.cross(r)).length(); // relative speed (c=1)
                tau += r.cross(d * catch);
                count += 1.0;
            }
        }
        if count > 0.0 {
            tau / count
        } else {
            DVec3::ZERO
        }
    }

    /// Integrate one step. `time_scale` matches the spin engine's convention
    /// (natural radians per real second at v = c for orbit_radius 1).
    pub fn integrate(&mut self, dt: f64, time_scale: f64) {
        self.position += self.lin_velocity * dt;

        if self.ang_velocity.length_squared() > 1e-30 {
            let dq = DQuat::from_scaled_axis(self.ang_velocity * dt);
            self.orientation = (dq * self.orientation).normalize();
        }

        let omega_swing = self.outer_spin * time_scale / self.swing_orbit_radius;
        self.swing_phase += omega_swing * dt;
    }

    /// orbit_radius of the swing level (converts `outer_spin` to a swing rate).
    pub fn swing_orbit_radius(&self) -> f64 {
        self.swing_orbit_radius
    }

    /// Moment of inertia about the pole axis (geometry-derived).
    pub fn i_spin(&self) -> f64 {
        self.i_spin
    }
}

/// xorshift64* → uniform f64 in [0, 1). Deterministic; no std rng dependency.
pub fn xorshift64(state: &mut u64) -> f64 {
    let mut x = *state;
    x ^= x >> 12;
    x ^= x << 25;
    x ^= x >> 27;
    *state = x;
    ((x.wrapping_mul(0x2545F4914F6CDD1D) >> 11) as f64) / ((1u64 << 53) as f64)
}

/// Uniform point on the unit sphere.
pub fn rand_unit_vec(state: &mut u64) -> DVec3 {
    let z = xorshift64(state) * 2.0 - 1.0;
    let phi = xorshift64(state) * std::f64::consts::TAU;
    let r = (1.0 - z * z).max(0.0).sqrt();
    DVec3::new(r * phi.cos(), r * phi.sin(), z)
}

/// Temperature/density presets for the ambient field (CM-1). "heat is photon
/// density" (heat.html), so temperature scales the number density (→ flux).
/// Room is anchored to the Mathis average charge-field density (photon3.pdf,
/// see `units::AMBIENT_PHOTON_NUMBER_DENSITY_PER_M3`); `density_multiple` is
/// relative to that baseline. Composition follows PHYSICS_REFERENCE §5: Earth
/// presets 2/3 photon, the void 50/50.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum FieldPreset {
    ColdVoid,
    Room293K,
    Hot1000K,
    SolarCore,
}

impl FieldPreset {
    pub fn label(self) -> &'static str {
        match self {
            FieldPreset::ColdVoid => "cold_void",
            FieldPreset::Room293K => "room_293K",
            FieldPreset::Hot1000K => "hot_1000K",
            FieldPreset::SolarCore => "solar_core",
        }
    }

    /// Number density relative to the Mathis ambient baseline (Room = 1×). Solar
    /// is display-compressed; the true Sun/Earth ratio is 84,986× (pause.html,
    /// `units::SUN_EARTH_CHARGE_DENSITY_RATIO`).
    pub fn density_multiple(self) -> f64 {
        match self {
            FieldPreset::ColdVoid => 0.03,
            FieldPreset::Room293K => 1.0,
            FieldPreset::Hot1000K => 10.0,
            FieldPreset::SolarCore => 100.0,
        }
    }

    /// Photon fraction (χ = +1); the rest are antiphotons.
    pub fn photon_fraction(self) -> f64 {
        match self {
            FieldPreset::ColdVoid => 0.5,
            _ => 2.0 / 3.0,
        }
    }

    pub fn cycle(self) -> Self {
        match self {
            FieldPreset::ColdVoid => FieldPreset::Room293K,
            FieldPreset::Room293K => FieldPreset::Hot1000K,
            FieldPreset::Hot1000K => FieldPreset::SolarCore,
            FieldPreset::SolarCore => FieldPreset::ColdVoid,
        }
    }
}

/// A parameterized ambient charge field that pelts a particle with photons.
///
/// This is the CM-1/CM-2 driver: `flux` scales with temperature/density,
/// `photon_fraction` is the photon:antiphoton mix (Earth ≈ 2/3 photon per
/// PHYSICS_REFERENCE §5; 0.5 = balanced void), and `direction_bias` makes the
/// field directional (None = isotropic). Momentum and spin_gain are scaffold
/// magnitudes awaiting the sourced photon momentum + measured coupling.
pub struct AmbientField {
    /// Expected photon contacts per unit time.
    pub flux: f64,
    /// Fraction of photons (chirality +1); the rest are antiphotons (−1).
    pub photon_fraction: f64,
    /// None = isotropic; Some(dir) = all photons travel along `dir`.
    pub direction_bias: Option<DVec3>,
    /// Per-photon momentum (natural units, placeholder).
    pub momentum: f64,
    /// Lumped spin-transfer coefficient (the thing CM-2 measures).
    pub spin_gain: f64,
    /// Swing-rate time_scale, matching `CalibrationParticle::integrate`.
    pub time_scale: f64,
}

impl AmbientField {
    /// Room-temperature Earth-like default: 2/3 photon, isotropic. Magnitudes
    /// are placeholders (calibrated in CM-1).
    pub fn room_default(time_scale: f64) -> Self {
        Self {
            flux: 200.0,
            photon_fraction: 2.0 / 3.0,
            direction_bias: None,
            momentum: 0.02,
            spin_gain: 0.5,
            time_scale,
        }
    }

    /// Deliver one time step's worth of photon contacts to `p`.
    pub fn tick(&self, p: &mut CalibrationParticle, dt: f64, rng: &mut u64) {
        let expected = (self.flux * dt).max(0.0);
        let mut n = expected.floor() as u64;
        if xorshift64(rng) < expected - expected.floor() {
            n += 1;
        }
        if n == 0 {
            return;
        }
        let pts = p.world_points();
        if pts.is_empty() {
            return;
        }
        for _ in 0..n {
            let idx = (xorshift64(rng) * pts.len() as f64) as usize % pts.len();
            let contact = pts[idx];
            let dir = match self.direction_bias {
                Some(b) => b.normalize_or_zero(),
                None => rand_unit_vec(rng),
            };
            let chirality = if xorshift64(rng) < self.photon_fraction { 1.0 } else { -1.0 };
            p.apply_photon(contact, dir, chirality, self.momentum, self.spin_gain);
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    /// A small Z-pole particle for clean impulse-decomposition tests: loop L3,
    /// swing L4 (Z orbital) → pole = Z under identity orientation.
    fn z_pole_particle() -> CalibrationParticle {
        let p = CalibrationParticle::new(bake_loop(3, 16), 1.0);
        assert!((p.pole_axis() - DVec3::Z).length() < 1e-9, "expected Z pole");
        p
    }

    #[test]
    fn free_drift_advances_position() {
        let mut p = z_pole_particle();
        p.lin_velocity = DVec3::X * 2.0;
        p.integrate(0.5, 1.0);
        assert!((p.position - DVec3::X).length() < 1e-9, "got {:?}", p.position);
    }

    #[test]
    fn swing_phase_advances_with_outer_spin() {
        let mut p = z_pole_particle();
        p.outer_spin = 0.5;
        let ts = 1.0;
        p.integrate(2.0, ts);
        let expected = 0.5 * ts / p.swing_orbit_radius * 2.0;
        assert!((p.swing_phase - expected).abs() < 1e-12, "{} vs {}", p.swing_phase, expected);
    }

    #[test]
    fn radial_impulse_is_pure_linear() {
        let mut p = z_pole_particle();
        // r parallel to J → zero torque.
        let contact = p.position + DVec3::X;
        p.apply_impulse(contact, DVec3::X * 3.0);
        assert!((p.lin_velocity - DVec3::X * 3.0).length() < 1e-9);
        assert!(p.ang_velocity.length() < 1e-12, "no torque expected");
        assert!(p.outer_spin.abs() < 1e-12, "no spin pump expected");
    }

    #[test]
    fn pole_axis_torque_pumps_spin_only() {
        let mut p = z_pole_particle();
        // r = X, J = Y → r×J = Z (= pole): pure spin torque.
        let contact = p.position + DVec3::X;
        p.apply_impulse(contact, DVec3::Y * 0.001);
        assert!(p.outer_spin > 0.0, "spin should increase, got {}", p.outer_spin);
        assert!(p.ang_velocity.length() < 1e-12, "no tumble expected, got {:?}", p.ang_velocity);
    }

    #[test]
    fn transverse_torque_tumbles_only() {
        let mut p = z_pole_particle();
        // r = X, J = Z → r×J = -Y (⊥ pole): pure tumble.
        let contact = p.position + DVec3::X;
        p.apply_impulse(contact, DVec3::Z * 0.5);
        assert!(p.ang_velocity.length() > 1e-9, "tumble expected");
        assert!(p.outer_spin.abs() < 1e-12, "no spin pump, got {}", p.outer_spin);
    }

    /// The marquee physics check: the same hit with the same lever pumps a light
    /// loop's spin far more than a heavy baryon loop's — geometric inertia, no
    /// tuned constant. (electron L8 vs proton L12, both Z-pole precession swings.)
    #[test]
    fn heavy_loop_resists_spin_pump() {
        let mut light = CalibrationParticle::electron(1.0);
        let mut heavy = CalibrationParticle::proton(1.0);
        assert!((light.pole_axis() - DVec3::Z).length() < 1e-9);
        assert!((heavy.pole_axis() - DVec3::Z).length() < 1e-9);

        let contact_l = light.position + DVec3::X;
        let contact_h = heavy.position + DVec3::X;
        light.apply_impulse(contact_l, DVec3::Y * 0.001);
        heavy.apply_impulse(contact_h, DVec3::Y * 0.001);

        assert!(heavy.i_spin() > light.i_spin() * 10.0,
            "proton pole inertia should dwarf electron's: {} vs {}", heavy.i_spin(), light.i_spin());
        assert!(light.outer_spin > heavy.outer_spin * 10.0,
            "same hit should pump the light loop far more: {} vs {}",
            light.outer_spin, heavy.outer_spin);
    }

    #[test]
    fn reaching_c_flags_transmutation() {
        let mut p = z_pole_particle();
        // Hammer the pole axis until outer_spin saturates.
        for _ in 0..100_000 {
            let contact = p.position + DVec3::X;
            p.apply_impulse(contact, DVec3::Y * 1.0);
            if p.transmuted {
                break;
            }
        }
        assert!(p.transmuted, "should saturate to c and flag transmutation");
        assert!(p.outer_spin.abs() <= 1.0 + 1e-12, "outer_spin clamped at c");
    }

    #[test]
    fn swing_over_full_turn_returns_geometry() {
        // Sanity: swinging a precession loop by TAU returns the same world points.
        let mut p = CalibrationParticle::proton(1.0);
        let before = p.world_points();
        p.swing_phase = TAU;
        let after = p.world_points();
        let max_drift = before
            .iter()
            .zip(after.iter())
            .map(|(a, b)| (*a - *b).length())
            .fold(0.0f64, f64::max);
        assert!(max_drift < 1e-6, "TAU swing should be identity, max drift {}", max_drift);
    }

    // --- gear-rule collisions (model B) ---

    /// A photon opposing the surface motion (gears catch, f≈1) pumps far more
    /// spin than one co-moving with it (slips, f≈0).
    #[test]
    fn opposing_photon_catches_comoving_slips() {
        // Z-pole particle spinning +; at contact +X the surface moves +Y.
        let mut opp = z_pole_particle();
        let mut com = z_pole_particle();
        opp.outer_spin = 0.5;
        com.outer_spin = 0.5;
        let contact_o = opp.position + DVec3::X;
        let contact_c = com.position + DVec3::X;
        // opposing surface (+Y) → dir -Y; co-moving → dir +Y.
        opp.apply_photon(contact_o, -DVec3::Y, 1.0, 0.02, 0.5);
        com.apply_photon(contact_c, DVec3::Y, 1.0, 0.02, 0.5);
        let d_opp = opp.outer_spin - 0.5;
        let d_com = com.outer_spin - 0.5;
        assert!(d_opp > 0.0, "opposing hit should pump spin, got {}", d_opp);
        assert!(d_opp > d_com.abs() * 100.0,
            "opposing (f≈1) should dwarf co-moving (f≈0): {} vs {}", d_opp, d_com);
    }

    /// Chirality sets the sign of the spin pump (augment vs cancel).
    #[test]
    fn chirality_flips_spin_pump_sign() {
        let mut pho = z_pole_particle();
        let mut anti = z_pole_particle();
        pho.outer_spin = 0.5;
        anti.outer_spin = 0.5;
        let c = DVec3::X;
        pho.apply_photon(pho.position + c, -DVec3::Y, 1.0, 0.02, 0.5);
        anti.apply_photon(anti.position + c, -DVec3::Y, -1.0, 0.02, 0.5);
        assert!((pho.outer_spin - 0.5) > 0.0);
        assert!((anti.outer_spin - 0.5) < 0.0);
    }

    /// An imbalanced field (photon-rich) spins the particle up from rest, while
    /// a balanced field leaves it near zero — the augment/cancel statistic.
    #[test]
    fn imbalanced_field_spins_up_balanced_does_not() {
        let make = || CalibrationParticle::new(bake_loop(3, 32), 1.0);
        let ts = 1.0;

        let mut balanced_p = make();
        let bal = AmbientField {
            flux: 500.0, photon_fraction: 0.5, direction_bias: None,
            momentum: 0.0, spin_gain: 0.5, time_scale: ts,
        };
        let mut rng_b = 0x1234_5678_9abc_def0u64;
        for _ in 0..400 {
            bal.tick(&mut balanced_p, 0.01, &mut rng_b);
            balanced_p.integrate(0.01, ts);
        }

        let mut imbalanced_p = make();
        let imb = AmbientField { photon_fraction: 0.9, ..bal_like(ts) };
        let mut rng_i = 0x1234_5678_9abc_def0u64;
        for _ in 0..400 {
            imb.tick(&mut imbalanced_p, 0.01, &mut rng_i);
            imbalanced_p.integrate(0.01, ts);
        }

        assert!(imbalanced_p.outer_spin.abs() > balanced_p.outer_spin.abs() * 5.0,
            "imbalanced should spin up far more than balanced: {} vs {}",
            imbalanced_p.outer_spin, balanced_p.outer_spin);
        assert!(imbalanced_p.outer_spin > 0.0, "photon-rich field spins up (+)");
    }

    /// A directional field pushes net linear drift along its travel direction.
    #[test]
    fn directional_field_drives_drift() {
        let mut p = CalibrationParticle::new(bake_loop(3, 32), 1.0);
        let field = AmbientField {
            flux: 500.0, photon_fraction: 0.5, direction_bias: Some(DVec3::X),
            momentum: 0.02, spin_gain: 0.0, time_scale: 1.0,
        };
        let mut rng = 0xABCD_1234_5678_9012u64;
        for _ in 0..200 {
            field.tick(&mut p, 0.01, &mut rng);
        }
        assert!(p.lin_velocity.x > 0.0, "field along +X should drift +X, got {:?}", p.lin_velocity);
        assert!(p.lin_velocity.x > p.lin_velocity.y.abs() * 5.0, "drift should be mostly along X");
    }

    #[test]
    fn field_presets_scale_flux_and_composition() {
        // Density rises monotonically with temperature.
        assert!(FieldPreset::ColdVoid.density_multiple() < FieldPreset::Room293K.density_multiple());
        assert!(FieldPreset::Room293K.density_multiple() < FieldPreset::Hot1000K.density_multiple());
        assert!(FieldPreset::Hot1000K.density_multiple() < FieldPreset::SolarCore.density_multiple());
        // Void is balanced; Earth presets are photon-rich (spin-up bias).
        assert_eq!(FieldPreset::ColdVoid.photon_fraction(), 0.5);
        assert!(FieldPreset::Room293K.photon_fraction() > 0.5);
        // Cycle visits all four and returns.
        let mut p = FieldPreset::ColdVoid;
        for _ in 0..4 {
            p = p.cycle();
        }
        assert_eq!(p, FieldPreset::ColdVoid);
    }

    /// At rest (outer_spin = 0) the surface is still, so the catch factor is
    /// uniform and hit-weighting reduces to the bare orientation histogram
    /// (the swing rotation is about the pole, preserving cos-to-pole). This is
    /// the control: any asymmetry at outer_spin > 0 is a genuine hit-weighting
    /// effect, not an artifact.
    #[test]
    fn participation_at_rest_matches_bare_histogram() {
        let p = CalibrationParticle::new(bake_loop(3, 128), 1.0);
        assert_eq!(p.outer_spin, 0.0);
        let bare = p.hitbox.orientation_histogram(16);
        let part = p.orientation_participation(Some(DVec3::X), 16, 8, 1);
        // Equal up to bin-boundary FP jitter (a stray sample near a cutoff);
        // any real hit-weighting skew would be an order of magnitude larger.
        for i in 0..16 {
            assert!((bare[i] - part[i]).abs() < 5e-3,
                "at rest, participation should match bare histogram at bin {}: {} vs {}",
                i, bare[i], part[i]);
        }
    }

    /// Participation is always a valid distribution, spinning or not.
    #[test]
    fn participation_is_a_distribution_when_spinning() {
        let mut p = CalibrationParticle::new(bake_loop(3, 128), 1.0);
        p.outer_spin = 0.9;
        let h = p.orientation_participation(Some(DVec3::X), 16, 8, 1);
        let sum: f64 = h.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "must sum to 1, got {}", sum);
        assert!(h.iter().all(|&v| v >= 0.0));
    }

    #[test]
    fn swing_drag_is_zero_at_rest() {
        let p = CalibrationParticle::new(bake_loop(3, 128), 1.0);
        assert_eq!(p.outer_spin, 0.0);
        assert!(p.swing_drag_torque(48, 8).abs() < 1e-9, "no drag at rest");
    }

    #[test]
    fn directional_torque_has_no_spin_drag_at_rest() {
        // At rest the catch is uniform, so the pole (spin) component vanishes.
        let p = CalibrationParticle::new(bake_loop(3, 128), 1.0);
        let tau = p.directional_torque(DVec3::X, 16);
        assert!(tau.dot(p.pole_axis()).abs() < 1e-9,
            "no spin drag at rest, got pole torque {}", tau.dot(p.pole_axis()));
    }

    fn bal_like(ts: f64) -> AmbientField {
        AmbientField {
            flux: 500.0, photon_fraction: 0.5, direction_bias: None,
            momentum: 0.0, spin_gain: 0.5, time_scale: ts,
        }
    }
}
