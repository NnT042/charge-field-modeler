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
use crate::units;

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

        // Model A: velocity-channel drag — the pole component of the momentum
        // torque, which model B discarded, restored with the relative-speed catch
        // weight |c·dir − v_surface| (c = 1 natural units). Head-on photons catch
        // harder than co-moving ones can chase, so an isotropic field nets a
        // spin-DOWN torque growing with outer_spin (see swing_drag_torque). At
        // rest the isotropic mean is zero; per-photon scatter remains as a
        // thermal floor.
        let omega = pole * (self.outer_spin / self.swing_orbit_radius);
        let v_surf = omega.cross(r);
        let catch = (dir - v_surf).length();
        let tau_pole = r.cross(dir * (momentum * catch)).dot(pole);
        self.outer_spin += self.spin_coupling * tau_pole / self.i_spin;

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

    // --- Alignment-channel diagnostics (measurement only) ---
    //
    // Geometry conventions shared by everything below:
    //   - Field photons travel along unit `d` (the field DIRECTION of travel).
    //     "Upstream" = −d (where the photons are coming FROM).
    //   - "North" aperture sits at `position + h·pole`, "south" at
    //     `position − h·pole`, where `h` = the loop's pole extent (see
    //     `pole_extent`).
    //   - Chirality gate (venus2.pdf): photons (χ=+1) are captured only at the
    //     SOUTH aperture; antiphotons (χ=−1) only at the NORTH aperture.
    //   - Lambert cosine capture weight: the south mouth's outward normal is
    //     `−pole`, and it faces the incoming stream when `(−pole)·(−d) > 0`
    //     i.e. `d·pole > 0`, so `w_s = max(0, d·pole)`; symmetrically
    //     `w_n = max(0, −d·pole)`.

    /// Half-extent of the baked loop along the pole axis, in the body
    /// (pre-orientation) frame: `h = max_i |local_point_i · local_pole|` with
    /// `local_pole = hitbox.swing_axis` (unnormalized-safe) and
    /// `local_point_i = points[i] + swing_offset` — derived the same way
    /// `pole_axis()` is, just without applying `self.orientation`. The swing
    /// rotation (about `swing_axis`) preserves the along-axis component, so no
    /// swing average is needed here. Returns `(h, used_fallback)`; if the loop
    /// is flat/equatorial (`h < 1e-6`, e.g. no measurable pole extent) falls
    /// back to `swing_orbit_radius()` and flags it.
    fn pole_extent(&self) -> (f64, bool) {
        let local_pole = self.hitbox.swing_axis.normalize_or_zero();
        let h = self
            .hitbox
            .points
            .iter()
            .map(|&p| (p + self.hitbox.swing_offset).dot(local_pole).abs())
            .fold(0.0_f64, f64::max);
        if h < 1e-6 {
            (self.swing_orbit_radius, true)
        } else {
            (h, false)
        }
    }

    /// DIAGNOSTIC (measurement, not a rule): net whole-body torque about
    /// `self.position` from a directional field, per unit photon momentum —
    /// like `directional_torque`, but each contribution is weighted by an
    /// upstream-exposure ("Lambert shadow") factor `expose = max(0, −d·r̂)`
    /// with `r̂ = (world_p − position).normalize_or_zero()`: points on the
    /// downstream side of the body (facing away from the stream) are treated
    /// as self-occluded and contribute nothing, points squarely upstream
    /// contribute fully. This tests whether pure geometric shadowing — with NO
    /// chirality gate and no aperture structure — produces a net alignment
    /// torque toward pole-parallel-to-field on its own.
    pub fn surface_shadow_torque(&self, field_dir: DVec3, swing_samples: usize) -> DVec3 {
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
                let world_p = self.position
                    + self.orientation * (srot * (self.hitbox.points[i] + self.hitbox.swing_offset));
                let r = world_p - self.position;
                let rhat = r.normalize_or_zero();
                let expose = (-d.dot(rhat)).max(0.0);
                let catch = (d - omega.cross(r)).length(); // relative speed (c=1)
                tau += r.cross(d * (catch * expose));
                count += 1.0;
            }
        }
        if count > 0.0 {
            tau / count
        } else {
            DVec3::ZERO
        }
    }

    /// DIAGNOSTIC (measurement, not a rule) — Scheme A: absorb at the aperture
    /// mouth. Chirality-gated, Lambert-weighted: photons (χ=+1, fraction
    /// `photon_fraction`) capture only at the south mouth `r_s = −h·pole`;
    /// antiphotons (fraction `1 − photon_fraction`) only at the north mouth
    /// `r_n = +h·pole`. `τ = photon_fraction·w_s·(r_s×d) + (1−photon_fraction)
    /// ·w_n·(r_n×d)`, per unit incoming photon momentum. No swing average —
    /// the apertures ride the pole and don't move under the swing rotation.
    pub fn intake_mouth_torque(&self, field_dir: DVec3, photon_fraction: f64) -> DVec3 {
        let d = field_dir.normalize_or_zero();
        let pole = self.pole_axis();
        let (h, _fallback) = self.pole_extent();
        let r_s = -pole * h;
        let r_n = pole * h;
        let w_s = d.dot(pole).max(0.0);
        let w_n = (-d.dot(pole)).max(0.0);
        photon_fraction * w_s * r_s.cross(d) + (1.0 - photon_fraction) * w_n * r_n.cross(d)
    }

    /// DIAGNOSTIC (measurement, not a rule) — Scheme B: captured charge
    /// channels from the mouth toward the core and re-emits as a symmetric
    /// equatorial ring (ring recoil cancels, so it deposits no torque of its
    /// own). Approximation: model the momentum deposit as applied at the
    /// midpoint of the mouth→center segment rather than at the mouth itself —
    /// same formula as `intake_mouth_torque` but with `r_s/2`, `r_n/2`.
    pub fn intake_channel_torque(&self, field_dir: DVec3, photon_fraction: f64) -> DVec3 {
        let d = field_dir.normalize_or_zero();
        let pole = self.pole_axis();
        let (h, _fallback) = self.pole_extent();
        let r_s = -pole * h * 0.5;
        let r_n = pole * h * 0.5;
        let w_s = d.dot(pole).max(0.0);
        let w_n = (-d.dot(pole)).max(0.0);
        photon_fraction * w_s * r_s.cross(d) + (1.0 - photon_fraction) * w_n * r_n.cross(d)
    }

    /// DIAGNOSTIC (measurement, not a rule) — Scheme C: pole-to-pole
    /// through-charge (venus2.pdf). The photon enters the mouth traveling
    /// along `d` and exits at the OPPOSITE pole redirected along that pole's
    /// outward direction. Honest two-point bookkeeping: the entry deposits
    /// `+d` at the entry mouth (`r_entry`), and launching the photon back out
    /// deposits the recoil `−exit_dir` at the exit tip (`r_exit`). For a south
    /// entry (photons): `r_entry = r_s`, `r_exit = r_n`, `exit_dir = +pole`.
    /// For a north entry (antiphotons): `r_entry = r_n`, `r_exit = r_s`,
    /// `exit_dir = −pole`. In both cases `r_exit` is itself along the pole
    /// axis, so `r_exit × (−exit_dir) = 0` exactly — the exit recoil produces
    /// NO torque (the exit tip sits ON the axis it's recoiling along). The
    /// term is kept explicit in code for bookkeeping honesty even though it
    /// algebraically vanishes. Net per-photon:
    /// `τ = photon_fraction·w_s·[r_s×d + r_n×(−pole)]
    ///    + (1−photon_fraction)·w_n·[r_n×d + r_s×pole]`.
    pub fn intake_through_torque(&self, field_dir: DVec3, photon_fraction: f64) -> DVec3 {
        let d = field_dir.normalize_or_zero();
        let pole = self.pole_axis();
        let (h, _fallback) = self.pole_extent();
        let r_s = -pole * h;
        let r_n = pole * h;
        let w_s = d.dot(pole).max(0.0);
        let w_n = (-d.dot(pole)).max(0.0);

        // South entry (photons): enters at r_s along d, exits at r_n along
        // +pole. Exit recoil = r_n × (−pole); r_n ∥ pole so this is exactly
        // zero.
        let south_entry = r_s.cross(d);
        let south_exit = r_n.cross(-pole);
        let south = south_entry + south_exit;

        // North entry (antiphotons): enters at r_n along d, exits at r_s
        // along −pole. Exit recoil = r_s × pole; r_s ∥ pole so this is also
        // exactly zero.
        let north_entry = r_n.cross(d);
        let north_exit = r_s.cross(pole);
        let north = north_entry + north_exit;

        photon_fraction * w_s * south + (1.0 - photon_fraction) * w_n * north
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

/// SI anchoring for an AmbientField (CM-1 second half). The sim can't deliver
/// ~1e11 contacts/s, so one sim photon stands in for `aggregation` (K) real
/// photons: per-photon magnitudes (momentum, spin_gain) scale UP by K, the
/// contact rate scales DOWN by K. Pump and drag both scale linearly in their
/// per-photon magnitude, so the equilibrium outer_spin is K-invariant (only
/// noise granularity grows with K — report_si_calibration checks this).
/// `seconds_per_time_unit` converts sim time to SI: flux_sim contacts per
/// natural time unit represent flux_si·density contacts per second.
pub struct SiCalibration {
    pub flux_si_hz: f64,
    pub aggregation: f64,
    pub momentum_natural: f64,
    pub seconds_per_time_unit: f64,
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

    /// Build a field whose magnitudes derive from the sourced SI anchors.
    /// `spin_gain_per_photon` is the one remaining free coupling (CM-2
    /// measures it); `sim_momentum` picks the working macro-photon impulse
    /// (K = sim_momentum / true momentum); `sim_flux` picks the working
    /// contact rate per natural time unit.
    pub fn from_si(
        preset: FieldPreset,
        particle_mass_kg: f64,
        spin_gain_per_photon: f64,
        sim_momentum: f64,
        sim_flux: f64,
        time_scale: f64,
    ) -> (Self, SiCalibration) {
        let momentum_natural = units::photon_momentum_natural(particle_mass_kg);
        let aggregation = sim_momentum / momentum_natural;
        let flux_si_hz = units::recycle_flux_hz(particle_mass_kg) * preset.density_multiple();
        let seconds_per_time_unit = sim_flux * aggregation / flux_si_hz;

        let field = Self {
            flux: sim_flux,
            photon_fraction: preset.photon_fraction(),
            direction_bias: None,
            momentum: sim_momentum,
            spin_gain: spin_gain_per_photon * aggregation,
            time_scale,
        };
        let cal = SiCalibration {
            flux_si_hz,
            aggregation,
            momentum_natural,
            seconds_per_time_unit,
        };
        (field, cal)
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
    ///
    /// NOTE (model A landed): `momentum=0.0` here, changed from the original
    /// `0.02`. The pump term this test isolates doesn't depend on `momentum` at
    /// all, but the new velocity-channel drag term (added in `apply_photon`
    /// alongside the pump, unconditionally) does — and for this exact geometry
    /// (opposing photon, outer_spin=0.5) the drag is a same-order-of-magnitude
    /// NEGATIVE contribution that ate most of the old 100x headroom (measured
    /// ratio dropped to ~27x with momentum=0.02). That's expected: it's the same
    /// self-limiter `report_spin_equilibrium` measures. Zeroing momentum isolates
    /// the chirality gear-catch mechanism this test is actually about, leaving
    /// the drag/pump interaction to the model-A-specific tests.
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
        opp.apply_photon(contact_o, -DVec3::Y, 1.0, 0.0, 0.5);
        com.apply_photon(contact_c, DVec3::Y, 1.0, 0.0, 0.5);
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

    /// The SI momentum-current identity `from_si` is built to preserve exactly:
    /// sim momentum-per-tick (`field.momentum · field.flux`) must equal the true
    /// per-photon momentum times the true SI contact rate times the sim-time→SI
    /// conversion, since `seconds_per_time_unit = sim_flux·K / flux_si_hz` and
    /// `K = sim_momentum / momentum_natural` cancel algebraically.
    #[test]
    fn si_momentum_current_preserved() {
        let (field, cal) = AmbientField::from_si(
            FieldPreset::Room293K,
            units::PROTON_MASS_KG,
            1.0,
            0.02,
            200.0,
            1.0,
        );
        let lhs = field.momentum * field.flux;
        let rhs = cal.momentum_natural * cal.flux_si_hz * cal.seconds_per_time_unit;
        let rel = (lhs - rhs).abs() / lhs.abs().max(1e-300);
        assert!(rel < 1e-9, "momentum current mismatch: lhs={:e} rhs={:e} rel={:e}", lhs, rhs, rel);
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

    // --- Alignment-channel diagnostics: fast unit tests ---

    /// A field exactly along (or exactly against) the untilted pole gives both
    /// aperture position vectors (`r_s`, `r_n`) parallel to `d`, so the lever
    /// `r × d` is zero regardless of aperture weights — the mouth-scheme
    /// torque must vanish at exact alignment.
    #[test]
    fn intake_mouth_zero_at_exact_alignment() {
        let p = CalibrationParticle::new(bake_loop(12, 256), 1.0); // untilted, pole = +Z
        let tau = p.intake_mouth_torque(DVec3::Z, 2.0 / 3.0);
        assert!(tau.length() < 1e-9, "pole parallel to field should give zero lever, got {:?}", tau);
        let tau2 = p.intake_mouth_torque(-DVec3::Z, 2.0 / 3.0);
        assert!(tau2.length() < 1e-9, "pole antiparallel to field should also give zero lever, got {:?}", tau2);
    }

    /// Chirality gate: at a tilt where the NORTH aperture is geometrically
    /// favored (`w_n > 0`, `w_s = 0`), a pure-photon field (`photon_fraction =
    /// 1.0`) must still return exactly zero — photons (χ=+1) are never
    /// captured at the north mouth, so its geometric weight is discarded
    /// entirely, not merely down-weighted. Flipping to pure-antiphoton
    /// (`photon_fraction = 0.0`) activates the north term; hand-computed
    /// against the documented formula.
    #[test]
    fn intake_mouth_chirality_gate_isolates_active_term() {
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        p.orientation = DQuat::from_axis_angle(DVec3::X, 90f64.to_radians()); // pole -> -Y
        let pole = p.pole_axis();
        // d chosen so the north aperture is geometrically favored: w_n > 0, w_s == 0.
        let d = (DVec3::Y + DVec3::Z).normalize();
        let w_s = d.dot(pole).max(0.0);
        let w_n = (-d.dot(pole)).max(0.0);
        assert_eq!(w_s, 0.0, "expected south weight clipped to zero for this d");
        assert!(w_n > 0.0, "expected north weight to be geometrically favored");

        // Pure photon field: south is zero (w_s=0) AND the chirality gate
        // discards the north term's coefficient -> net torque must be exactly zero.
        let pure_photon = p.intake_mouth_torque(d, 1.0);
        assert!(pure_photon.length() < 1e-9,
            "pure-photon field should ignore the geometrically-favored north aperture: got {:?}", pure_photon);

        // Pure antiphoton field activates the north-only term; hand-compute it.
        let (h, _) = p.pole_extent();
        let r_n = pole * h;
        let expected_north_only = w_n * r_n.cross(d);
        let pure_anti = p.intake_mouth_torque(d, 0.0);
        assert!((pure_anti - expected_north_only).length() < 1e-9,
            "expected north-only term {:?}, got {:?}", expected_north_only, pure_anti);
    }

    /// Regression check on `surface_shadow_torque`'s aggregate loop: with
    /// `swing_samples = 1` (no averaging, `srot` = identity) the method must
    /// equal a direct hand-rolled recomputation of its documented formula
    /// (`expose = max(0, -d·r̂)`, contribution `r × (d · catch · expose)`).
    /// Guards against sign/typo regressions without re-deriving the physics.
    #[test]
    fn surface_shadow_matches_hand_rolled_formula() {
        let mut p = CalibrationParticle::new(bake_loop(3, 32), 1.0);
        p.orientation = DQuat::from_axis_angle(DVec3::X, 0.7);
        p.outer_spin = 0.4;
        let d = DVec3::new(0.3, -0.6, 0.74).normalize();

        let pole = p.pole_axis();
        let omega = pole * (p.outer_spin / p.swing_orbit_radius());
        let n = p.hitbox.points.len();
        let mut expected = DVec3::ZERO;
        for i in 0..n {
            let world_p = p.position + p.orientation * (p.hitbox.points[i] + p.hitbox.swing_offset);
            let r = world_p - p.position;
            let rhat = r.normalize_or_zero();
            let expose = (-d.dot(rhat)).max(0.0);
            let catch = (d - omega.cross(r)).length();
            expected += r.cross(d * (catch * expose));
        }
        expected /= n as f64;

        let actual = p.surface_shadow_torque(d, 1);
        assert!((actual - expected).length() < 1e-9, "expected {:?}, got {:?}", expected, actual);
    }

    fn bal_like(ts: f64) -> AmbientField {
        AmbientField {
            flux: 500.0, photon_fraction: 0.5, direction_bias: None,
            momentum: 0.0, spin_gain: 0.5, time_scale: ts,
        }
    }

    // --- model A: velocity-channel drag ---

    /// Deterministic Fibonacci-sphere directions, same pattern as
    /// `swing_drag_torque`'s isotropic sampling.
    fn fib_sphere(m: usize) -> Vec<DVec3> {
        (0..m)
            .map(|k| {
                let z = 1.0 - 2.0 * (k as f64 + 0.5) / m as f64;
                let r = (1.0 - z * z).max(0.0).sqrt();
                let phi = k as f64 * 2.399_963_229_728_653;
                DVec3::new(r * phi.cos(), r * phi.sin(), z)
            })
            .collect()
    }

    /// At rest, an isotropic bombardment's pole-axis drag torque must average to
    /// zero (per-photon it's nonzero — that's the thermal floor — but the mean
    /// over an isotropic direction set at a fixed contact point vanishes). Each
    /// photon is applied to a *fresh* rest-state particle so later hits don't
    /// contaminate the "at rest" measurement (matches how `swing_drag_torque`
    /// itself is a stateless snapshot at a given `outer_spin`).
    #[test]
    fn drag_has_zero_mean_at_rest() {
        let template = CalibrationParticle::new(bake_loop(3, 128), 1.0);
        let pts = template.world_points();
        let stride = (pts.len() / 12).max(1);
        let contacts: Vec<DVec3> = pts.iter().step_by(stride).copied().collect();
        let dirs = fib_sphere(64);
        let momentum = 0.02;

        let mut sum = 0.0;
        let mut sum_sq = 0.0;
        let mut n = 0u64;
        for &contact in &contacts {
            for &dir in &dirs {
                let mut p = CalibrationParticle::new(bake_loop(3, 128), 1.0);
                assert_eq!(p.outer_spin, 0.0);
                p.apply_photon(contact, dir, 1.0, momentum, 0.0);
                let d = p.outer_spin; // delta from the rest state
                sum += d;
                sum_sq += d * d;
                n += 1;
            }
        }
        let mean = sum / n as f64;
        let rms = (sum_sq / n as f64).sqrt();
        assert!(rms > 1e-12, "expected nonzero per-photon thermal scatter, got rms {}", rms);
        assert!(mean.abs() < 0.05 * rms,
            "isotropic drag should have ~zero mean at rest: mean={} rms={} (n={})", mean, rms, n);
    }

    /// A balanced (photon_fraction 0.5, spin_gain 0.0) isotropic field drags a
    /// spun-up particle back down toward zero — the drag term alone, no pump.
    #[test]
    fn balanced_field_drags_spin_down() {
        // Smaller bake than the full proton (still L12/Z-pole) to keep the test fast.
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        assert!((p.pole_axis() - DVec3::Z).length() < 1e-9);
        p.outer_spin = 0.8;
        let field = AmbientField {
            flux: 200.0,
            photon_fraction: 0.5,
            direction_bias: None,
            momentum: 0.02,
            spin_gain: 0.0,
            time_scale: 1.0,
        };
        let dt = 0.02;
        let mut rng = 0x9E37_79B9_7F4A_7C15u64;
        let mut settled = false;
        for _ in 0..20_000 {
            field.tick(&mut p, dt, &mut rng);
            p.integrate(dt, field.time_scale);
            if p.outer_spin < 0.4 {
                settled = true;
                break;
            }
        }
        assert!(settled, "drag should pull outer_spin below 0.4, stuck at {}", p.outer_spin);
        assert!(!p.transmuted, "pure drag should never transmute, got outer_spin {}", p.outer_spin);
    }

    /// Mirror of the above starting from negative spin: drag restores toward
    /// zero from either side, no sign-dependent gate.
    #[test]
    fn drag_restores_negative_spin_toward_zero() {
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        p.outer_spin = -0.8;
        let field = AmbientField {
            flux: 200.0,
            photon_fraction: 0.5,
            direction_bias: None,
            momentum: 0.02,
            spin_gain: 0.0,
            time_scale: 1.0,
        };
        let dt = 0.02;
        let mut rng = 0x1D87_2B23_45FF_0011u64;
        let mut settled = false;
        for _ in 0..20_000 {
            field.tick(&mut p, dt, &mut rng);
            p.integrate(dt, field.time_scale);
            if p.outer_spin > -0.4 {
                settled = true;
                break;
            }
        }
        assert!(settled, "drag should pull outer_spin above -0.4, stuck at {}", p.outer_spin);
        assert!(!p.transmuted, "pure drag should never transmute, got outer_spin {}", p.outer_spin);
    }

    /// Evaluate `swing_drag_torque` at a trial `outer_spin` without permanently
    /// disturbing the particle (restores afterward). Returns torque per photon
    /// WITHOUT the momentum factor (see the diagnostic's doc comment).
    fn drag_at_spin(p: &mut CalibrationParticle, s: f64) -> f64 {
        let saved = p.outer_spin;
        p.outer_spin = s;
        let t = p.swing_drag_torque(48, 16);
        p.outer_spin = saved;
        t
    }

    /// Bisect for the equilibrium `outer_spin` where the chirality pump balances
    /// the velocity-channel drag: `pump + spin_coupling·momentum·drag_at_spin(s)
    /// == 0`. `pump` is signed (positive for a photon-rich field); the drag term
    /// is negative-growing for positive spin (see `swing_drag_torque`'s doc), so
    /// for a photon-rich field the root sits at positive `s`.
    fn bisect_equilibrium(p: &mut CalibrationParticle, pump: f64, spin_coupling: f64, momentum: f64) -> f64 {
        let g = |p: &mut CalibrationParticle, s: f64| -> f64 { pump + spin_coupling * momentum * drag_at_spin(p, s) };
        let sign = pump.signum();
        let (mut lo, mut hi) = if sign >= 0.0 { (0.0, 0.995) } else { (-0.995, 0.0) };
        let mut glo = g(p, lo);
        let ghi = g(p, hi);
        if glo == 0.0 {
            return lo;
        }
        if glo.signum() == ghi.signum() {
            // No sign change in range (drag never catches the pump) — report the
            // boundary closest to balance so the caller can see it's saturating.
            return if ghi.abs() < glo.abs() { hi } else { lo };
        }
        for _ in 0..60 {
            let mid = 0.5 * (lo + hi);
            let gmid = g(p, mid);
            if gmid.signum() == glo.signum() {
                lo = mid;
                glo = gmid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }

    /// REPORT (ignored by default): validates the model-A scaling law end to
    /// end. Predicted equilibrium from the augment/cancel pump vs. the measured
    /// `swing_drag_torque` drag (via bisection) against the actual settled
    /// `outer_spin` of a ticked `AmbientField` simulation, swept over field
    /// imbalance and density. Prints tables; no hard physics asserts (see
    /// `report_*` convention elsewhere in the repo) beyond basic sanity.
    #[test]
    #[ignore]
    fn report_spin_equilibrium() {
        // Proton-shaped bake (L12/Z-pole), reduced samples for tractable runtime;
        // geometry (pole axis, swing kind) is identical to the full bake.
        let make_particle = || CalibrationParticle::new(bake_loop(12, 256), 1.0);

        let run_sweep = |label: &str, rows: Vec<(f64, f64, f64)>| {
            // rows: (photon_fraction, flux, spin_gain)
            println!(
                "\n== {label} ==\n{:>10} {:>8} {:>10} {:>12} {:>14} {:>12}",
                "p_photon", "flux", "spin_gain", "predicted_s*", "sim_settled_s", "settle_step"
            );
            for (photon_fraction, flux, spin_gain) in rows {
                let momentum = 0.02;
                let mut p = make_particle();
                let coupling = p.spin_coupling;
                let pump = 0.5 * spin_gain * (2.0 * photon_fraction - 1.0);
                let predicted = bisect_equilibrium(&mut p, pump, coupling, momentum);

                let mut sim = make_particle();
                let field = AmbientField {
                    flux,
                    photon_fraction,
                    direction_bias: None,
                    momentum,
                    spin_gain,
                    time_scale: 1.0,
                };
                let dt = (2.0 / flux).clamp(0.005, 0.05); // aim ~1-4 photons/tick
                let steps = 150_000usize;
                let mut rng = 0xC0FF_EE12_3456_789Au64;
                let mut trace = Vec::with_capacity(steps);
                for _ in 0..steps {
                    field.tick(&mut sim, dt, &mut rng);
                    sim.integrate(dt, field.time_scale);
                    trace.push(sim.outer_spin);
                }
                let tail_start = steps - steps / 5; // last 20%
                let settled: f64 = trace[tail_start..].iter().sum::<f64>() / (steps - tail_start) as f64;
                let target = 0.9 * settled;
                let settle_step = trace
                    .iter()
                    .position(|&v| {
                        if settled >= 0.0 {
                            v >= target
                        } else {
                            v <= target
                        }
                    })
                    .unwrap_or(steps);

                println!(
                    "{:>10.3} {:>8.1} {:>10.3} {:>12.4} {:>14.4} {:>12}",
                    photon_fraction, flux, spin_gain, predicted, settled, settle_step
                );
            }
        };

        // Imbalance sweep at flux 200.
        run_sweep(
            "imbalance sweep (flux=200)",
            vec![
                (0.5, 200.0, 2.0),
                (0.583, 200.0, 2.0),
                (0.667, 200.0, 2.0),
                (0.75, 200.0, 2.0),
                (1.0, 200.0, 2.0),
            ],
        );

        // Density sweep at photon_fraction 0.667.
        run_sweep(
            "density sweep (photon_fraction=0.667)",
            vec![(0.667, 50.0, 2.0), (0.667, 200.0, 2.0), (0.667, 800.0, 2.0)],
        );
    }

    // --- CM-1 second half: SI flux calibration ---

    /// Ticks `field` against a fresh particle from `make` for 150k steps and
    /// returns `(settled_s*, settle_step, dt)`, same convention as
    /// `report_spin_equilibrium`'s inner sweep.
    fn run_settle(
        make: &dyn Fn() -> CalibrationParticle,
        field: &AmbientField,
        seed: u64,
        steps: usize,
    ) -> (f64, usize, f64, bool) {
        let dt = (2.0 / field.flux).clamp(0.005, 0.05);
        let mut rng = seed;
        let mut sim = make();
        let mut trace = Vec::with_capacity(steps);
        for _ in 0..steps {
            field.tick(&mut sim, dt, &mut rng);
            sim.integrate(dt, field.time_scale);
            trace.push(sim.outer_spin);
        }
        let tail_start = steps - steps / 5; // last 20%
        let settled: f64 = trace[tail_start..].iter().sum::<f64>() / (steps - tail_start) as f64;
        let target = 0.9 * settled;
        let settle_step = trace
            .iter()
            .position(|&v| if settled >= 0.0 { v >= target } else { v <= target })
            .unwrap_or(steps);
        (settled, settle_step, dt, sim.transmuted)
    }

    /// REPORT (ignored by default): the natural↔SI flux calibration end to end.
    /// Table 1 prints the sourced anchors per particle (mass, true recycle flux,
    /// true per-photon momentum, the macro-photon aggregation K, and the sim↔SI
    /// time conversion) for the working sim knobs (momentum=0.02, flux=200).
    /// Table 2 runs the same settle loop as `report_spin_equilibrium` through
    /// `from_si`, converting settle time to SI seconds. Table 3 checks that the
    /// settled equilibrium is K-invariant — quartering `sim_momentum` (and
    /// quadrupling `sim_flux` to hold the represented physical density fixed)
    /// must reproduce the same settled `outer_spin` within noise. Run B gets
    /// 2× the steps so both runs receive the same physical photon dose (see
    /// the inline comment at the run_settle calls).
    #[test]
    #[ignore]
    fn report_si_calibration() {
        // --- Table 1: sourced anchors per particle ---
        println!(
            "\n== Table 1: SI anchors (Room, sim_momentum=0.02, flux_sim=200) ==\n{:>10} {:>12} {:>14} {:>16} {:>12} {:>18}",
            "particle", "mass_kg", "flux_si_Hz", "momentum_nat", "K", "sec_per_tu"
        );
        let particles: [(&str, f64); 3] = [
            ("proton", units::PROTON_MASS_KG),
            ("neutron", units::NEUTRON_MASS_KG),
            ("electron", units::ELECTRON_MASS_KG),
        ];
        for (label, mass_kg) in particles {
            let momentum_natural = units::photon_momentum_natural(mass_kg);
            let flux_si = units::recycle_flux_hz(mass_kg) * FieldPreset::Room293K.density_multiple();
            let k = 0.02 / momentum_natural;
            let sec_per_tu = 200.0 * k / flux_si;
            println!(
                "{:>10} {:>12.4e} {:>14.4e} {:>16.4e} {:>12.4e} {:>18.4e}",
                label, mass_kg, flux_si, momentum_natural, k, sec_per_tu
            );
        }

        // --- Table 2: Room-preset settling per particle ---
        struct ParticleSpec {
            label: &'static str,
            mass_kg: f64,
            make: Box<dyn Fn() -> CalibrationParticle>,
        }
        let specs: Vec<ParticleSpec> = vec![
            ParticleSpec {
                label: "proton",
                mass_kg: units::PROTON_MASS_KG,
                make: Box::new(|| CalibrationParticle::new(bake_loop(12, 256), 1.0)),
            },
            ParticleSpec {
                label: "neutron",
                mass_kg: units::NEUTRON_MASS_KG,
                make: Box::new(|| CalibrationParticle::new(bake_loop(11, 256), 1.0)),
            },
            ParticleSpec {
                label: "electron",
                mass_kg: units::ELECTRON_MASS_KG,
                make: Box::new(|| CalibrationParticle::new(bake_loop(8, 256), 1.0)),
            },
        ];

        println!(
            "\n== Table 2: Room-preset settling ==\n{:>10} {:>12} {:>12} {:>14} {:>16} {:>18}",
            "particle", "settled_s*", "settle_step", "settle_t_nat", "settle_t_SI_s", "per_photon_gain"
        );
        for spec in &specs {
            let momentum_natural = units::photon_momentum_natural(spec.mass_kg);
            let k = 0.02 / momentum_natural;
            let spin_gain_per_photon = 2.0 / k;
            let (field, cal) = AmbientField::from_si(
                FieldPreset::Room293K,
                spec.mass_kg,
                spin_gain_per_photon,
                0.02,
                200.0,
                1.0,
            );
            let (settled, settle_step, dt, transmuted) =
                run_settle(&*spec.make, &field, 0xC0FF_EE12_3456_789Au64, 150_000);
            let settle_t_nat = settle_step as f64 * dt;
            let settle_t_si = settle_t_nat * cal.seconds_per_time_unit;
            println!(
                "{:>10} {:>12.4} {:>12} {:>14.4} {:>16.4e} {:>18.4e}{}",
                spec.label,
                settled,
                settle_step,
                settle_t_nat,
                settle_t_si,
                spin_gain_per_photon,
                if transmuted { "  (TRANSMUTED — hit c)" } else { "" }
            );
        }

        // --- Table 3: K-invariance (proton, Room) ---
        println!("\n== Table 3: K-invariance (proton, Room) ==");
        let mass_kg = units::PROTON_MASS_KG;
        let momentum_natural = units::photon_momentum_natural(mass_kg);
        let k_a = 0.02 / momentum_natural;
        let spin_gain_per_photon = 2.0 / k_a;

        let (field_a, cal_a) =
            AmbientField::from_si(FieldPreset::Room293K, mass_kg, spin_gain_per_photon, 0.02, 200.0, 1.0);
        let (field_b, cal_b) =
            AmbientField::from_si(FieldPreset::Room293K, mass_kg, spin_gain_per_photon, 0.005, 800.0, 1.0);

        let make_proton: Box<dyn Fn() -> CalibrationParticle> =
            Box::new(|| CalibrationParticle::new(bake_loop(12, 256), 1.0));
        // Equal PHYSICAL photon dose, not equal step count: per step, run A
        // delivers dt·flux·K = 0.01·200·K = 2K real photons while run B
        // delivers 0.005·800·(K/4) = 1K — half the dose. With equal steps run
        // B is still converging when the tail window opens, which shows up as
        // a spurious ~15% "K-dependence" (seen 2026-07-16). Doubling B's steps
        // equalizes the represented physical exposure.
        let (settled_a, step_a, _dt_a, _tr_a) =
            run_settle(&*make_proton, &field_a, 0xC0FF_EE12_3456_789Au64, 150_000);
        let (settled_b, step_b, _dt_b, _tr_b) =
            run_settle(&*make_proton, &field_b, 0xC0FF_EE12_3456_789Au64, 300_000);

        let rel_diff = (settled_a - settled_b).abs() / settled_a.abs().max(settled_b.abs()).max(1e-12);
        println!(
            "{:>8} {:>10} {:>10} {:>10} {:>12}",
            "run", "K", "flux_sim", "settled_s*", "settle_step"
        );
        println!("{:>8} {:>10.4e} {:>10.1} {:>10.4} {:>12}", "A", cal_a.aggregation, field_a.flux, settled_a, step_a);
        println!("{:>8} {:>10.4e} {:>10.1} {:>10.4} {:>12}", "B", cal_b.aggregation, field_b.flux, settled_b, step_b);
        println!("relative difference = {:.4}", rel_diff);

        assert!(
            rel_diff < 0.05,
            "K-invariance violated: settled_a={} settled_b={} rel_diff={}",
            settled_a, settled_b, rel_diff
        );
    }

    // --- Alignment-channel report ---

    /// Locate zero-crossing equilibria in a `(theta, value)` series and
    /// classify each as stable/unstable. Because the intake schemes are exact
    /// (not swing-averaged noise), some crossings sit exactly ON a sample
    /// point (e.g. theta=90 where both aperture weights are analytically
    /// zero) rather than between two nonzero-opposite-sign samples, so both
    /// cases are handled: an exact-zero grid point is an equilibrium if its
    /// nearest nonzero neighbors on either side disagree in sign; otherwise a
    /// linear-interpolated crossing between two nonzero opposite-sign samples.
    /// A crossing from `+` (below theta*) to `-` (above theta*) is STABLE
    /// under `dtheta/dt = tau_x` (positive tau_x increases theta, so it pushes
    /// theta UP toward theta* from below and DOWN toward theta* from above);
    /// the reverse crossing is UNSTABLE (pushes away from theta* both sides).
    fn find_equilibria(thetas: &[f64], values: &[f64]) -> Vec<(f64, bool)> {
        let is_zero = |v: f64| v.abs() < 1e-9;
        let n = values.len();
        let mut out = Vec::new();
        for i in 0..n {
            if is_zero(values[i]) {
                let before = (0..i).rev().map(|j| values[j]).find(|&v| !is_zero(v));
                let after = (i + 1..n).map(|j| values[j]).find(|&v| !is_zero(v));
                if let (Some(b), Some(a)) = (before, after) {
                    if b.signum() != a.signum() {
                        out.push((thetas[i], b > 0.0 && a < 0.0));
                    }
                }
            } else if i + 1 < n && !is_zero(values[i + 1]) && values[i].signum() != values[i + 1].signum() {
                let (t0, t1, v0, v1) = (thetas[i], thetas[i + 1], values[i], values[i + 1]);
                let theta_star = t0 + (t1 - t0) * v0 / (v0 - v1);
                out.push((theta_star, v0 > 0.0 && v1 < 0.0));
            }
        }
        // De-dup near-identical crossings (e.g. an exact-zero grid point that
        // also triggers the interpolated branch on an adjacent iteration).
        out.dedup_by(|a, b| (a.0 - b.0).abs() < 1e-6);
        out
    }

    /// Translate an equilibrium theta (degrees) to a plain-English polarity
    /// reading, per the tilt convention documented in `report_alignment_channels`.
    fn polarity_of(theta_deg: f64) -> &'static str {
        if theta_deg.abs() < 1.0 {
            "north mouth upstream"
        } else if (theta_deg - 180.0).abs() < 1.0 {
            "south mouth upstream"
        } else if (theta_deg - 90.0).abs() < 1.0 {
            "axis perpendicular to field"
        } else {
            "intermediate tilt"
        }
    }

    /// Print one verdict line for a channel's (theta, tau_x) series.
    fn print_verdict(name: &str, thetas: &[f64], values: &[f64]) {
        let is_zero = |v: f64| v.abs() < 1e-9;
        if values.iter().all(|&v| is_zero(v)) {
            println!("{:<16} no alignment torque (zero everywhere)", name);
            return;
        }
        // Domain ENDPOINTS (theta = 0 and 180) are equilibria whenever tau_x
        // vanishes there (it must, for a pole||field configuration), but the
        // interior-crossing scan can't see them — no sample exists past the
        // edge. Classify them from the adjacent interior sign: negative tau_x
        // just above 0 pushes theta back DOWN to 0 (attractor); positive tau_x
        // just below 180 pushes theta UP to 180 (attractor). This is where the
        // real story lives for chirality-blind channels (2026-07-16 sweep:
        // surface_shadow's attractors at 0/180 were invisible to the interior
        // scan, leaving only the 90-degree repeller in the printout).
        if let Some(&first) = values.iter().find(|v| !is_zero(**v)) {
            if is_zero(values[0]) {
                let stable = first < 0.0;
                println!(
                    "{:<16} theta_eq=   0.00 deg  {:<9} -> {}",
                    name,
                    if stable { "STABLE" } else { "UNSTABLE" },
                    polarity_of(0.0)
                );
            }
        }
        let eqs = find_equilibria(thetas, values);
        for (theta_star, stable) in eqs {
            println!(
                "{:<16} theta_eq={:>7.2} deg  {:<9} -> {}",
                name,
                theta_star,
                if stable { "STABLE" } else { "UNSTABLE" },
                polarity_of(theta_star)
            );
        }
        if let Some(&last) = values.iter().rev().find(|v| !is_zero(**v)) {
            if is_zero(values[values.len() - 1]) {
                let stable = last > 0.0;
                println!(
                    "{:<16} theta_eq= 180.00 deg  {:<9} -> {}",
                    name,
                    if stable { "STABLE" } else { "UNSTABLE" },
                    polarity_of(180.0)
                );
            }
        }
    }

    /// REPORT (ignored by default): measures four candidate alignment-torque
    /// channels — none of which touch any collision rule — across a tilt sweep
    /// against a FIXED field direction, to see which (if any) produce a
    /// restoring torque toward pole-parallel-to-field, and toward which
    /// polarity. `directional_torque` (existing diagnostic) already found ZERO
    /// alignment torque from the pure velocity-catch mechanism at every tilt;
    /// this report checks four alternatives: (1) upstream self-shadowing on
    /// the bare surface, and Scheme A/B/C intake-aperture models gated by
    /// chirality (venus2.pdf).
    #[test]
    #[ignore]
    fn report_alignment_channels() {
        let d = DVec3::NEG_Z; // field DIRECTION of travel: photons move -Z (stream source sits at +Z, above).

        println!("\n=== Alignment channel report ===");
        println!("Field direction d = {:?} (photons travel toward -Z; the stream source is above, at +Z).", d);
        println!("Tilt: p.orientation = DQuat::from_axis_angle(+X, theta). Untilted pole = +Z.");
        println!("  theta=0    -> north pole (+Z) faces the stream source (mouth-upstream for ANTIPHOTON capture at N).");
        println!("  theta=180  -> south pole faces the stream source (mouth-upstream for PHOTON capture at S).");
        println!("Sign convention: the tilt is a rotation about +X, so ang_velocity_x = tau_x / I_transverse (I>0)");
        println!("  and dtheta/dt = ang_velocity_x for this geometry. => POSITIVE tau_x INCREASES theta.");
        println!("  A STABLE equilibrium at theta* therefore needs tau_x > 0 for theta < theta* (pushes UP toward");
        println!("  theta*) and tau_x < 0 for theta > theta* (pushes DOWN toward theta*) -- a +-to-- crossing.");

        let make = |theta_deg: f64| -> CalibrationParticle {
            let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
            p.orientation = DQuat::from_axis_angle(DVec3::X, theta_deg.to_radians());
            p
        };

        let probe = CalibrationParticle::proton(1.0);
        let (h, fell_back) = probe.pole_extent();
        println!(
            "\npole extent h = {:.6} natural units (fallback to swing_orbit_radius triggered: {})",
            h, fell_back
        );

        let thetas: Vec<f64> = (0..=180).step_by(15).map(|x| x as f64).collect();

        // --- Table 1: surface_shadow (3 outer_spin values) + intake_* @ Earth (2/3) ---
        println!("\n== Table 1: surface_shadow @ outer_spin={{0.0,0.6,0.9}}; intake_* @ photon_fraction=2/3 (Earth) ==");
        println!(
            "{:>6} {:>14} {:>14} {:>14} {:>14} {:>14} {:>14}",
            "theta", "shadow_s0.0", "shadow_s0.6", "shadow_s0.9", "mouth_2/3", "channel_2/3", "through_2/3"
        );
        let mut shadow_s00 = Vec::new();
        let mut shadow_s06 = Vec::new();
        let mut shadow_s09 = Vec::new();
        let mut mouth_23 = Vec::new();
        let mut channel_23 = Vec::new();
        let mut through_23 = Vec::new();
        for &td in &thetas {
            let mut p = make(td);
            p.outer_spin = 0.0;
            let s00 = p.surface_shadow_torque(d, 32).x;
            p.outer_spin = 0.6;
            let s06 = p.surface_shadow_torque(d, 32).x;
            p.outer_spin = 0.9;
            let s09 = p.surface_shadow_torque(d, 32).x;
            let m = p.intake_mouth_torque(d, 2.0 / 3.0).x;
            let c = p.intake_channel_torque(d, 2.0 / 3.0).x;
            let t = p.intake_through_torque(d, 2.0 / 3.0).x;
            println!(
                "{:>6.0} {:>14.6e} {:>14.6e} {:>14.6e} {:>14.6e} {:>14.6e} {:>14.6e}",
                td, s00, s06, s09, m, c, t
            );
            shadow_s00.push(s00);
            shadow_s06.push(s06);
            shadow_s09.push(s09);
            mouth_23.push(m);
            channel_23.push(c);
            through_23.push(t);
        }

        // --- Table 2: intake_* @ photon_fraction 0.5 and 1.0 ---
        println!("\n== Table 2: intake_* @ photon_fraction=0.5 (balanced) and 1.0 (pure photon) ==");
        println!(
            "{:>6} {:>14} {:>14} {:>14} {:>14} {:>14} {:>14}",
            "theta", "mouth_0.5", "channel_0.5", "through_0.5", "mouth_1.0", "channel_1.0", "through_1.0"
        );
        let mut mouth_50 = Vec::new();
        let mut channel_50 = Vec::new();
        let mut through_50 = Vec::new();
        let mut mouth_10 = Vec::new();
        let mut channel_10 = Vec::new();
        let mut through_10 = Vec::new();
        for &td in &thetas {
            let p = make(td);
            let m5 = p.intake_mouth_torque(d, 0.5).x;
            let c5 = p.intake_channel_torque(d, 0.5).x;
            let t5 = p.intake_through_torque(d, 0.5).x;
            let m1 = p.intake_mouth_torque(d, 1.0).x;
            let c1 = p.intake_channel_torque(d, 1.0).x;
            let t1 = p.intake_through_torque(d, 1.0).x;
            println!(
                "{:>6.0} {:>14.6e} {:>14.6e} {:>14.6e} {:>14.6e} {:>14.6e} {:>14.6e}",
                td, m5, c5, t5, m1, c1, t1
            );
            mouth_50.push(m5);
            channel_50.push(c5);
            through_50.push(t5);
            mouth_10.push(m1);
            channel_10.push(c1);
            through_10.push(t1);
        }

        // --- Table 3: verdicts ---
        println!("\n== Table 3: verdicts (representative columns: surface_shadow@outer_spin=0.9, intake_*@photon_fraction=2/3) ==");
        print_verdict("surface_shadow", &thetas, &shadow_s09);
        print_verdict("intake_mouth", &thetas, &mouth_23);
        print_verdict("intake_channel", &thetas, &channel_23);
        print_verdict("intake_through", &thetas, &through_23);

        // --- Sanity asserts ---
        let i0 = 0;
        let ilast = thetas.len() - 1;
        assert!(mouth_23[i0].abs() < 1e-12, "intake_mouth should be zero at theta=0, got {}", mouth_23[i0]);
        assert!(mouth_23[ilast].abs() < 1e-12, "intake_mouth should be zero at theta=180, got {}", mouth_23[ilast]);
        assert!(channel_23[i0].abs() < 1e-12, "intake_channel should be zero at theta=0, got {}", channel_23[i0]);
        assert!(channel_23[ilast].abs() < 1e-12, "intake_channel should be zero at theta=180, got {}", channel_23[ilast]);
        assert!(through_23[i0].abs() < 1e-12, "intake_through should be zero at theta=0, got {}", through_23[i0]);
        assert!(through_23[ilast].abs() < 1e-12, "intake_through should be zero at theta=180, got {}", through_23[ilast]);

        // surface_shadow at theta=0, outer_spin=0: small relative to its theta=90 magnitude.
        let idx90 = thetas.iter().position(|&t| t == 90.0).unwrap();
        let bound = (shadow_s00[idx90].abs() * 1e-3).max(1e-9);
        assert!(
            shadow_s00[i0].abs() < bound,
            "surface_shadow at theta=0,outer_spin=0 should be small vs theta=90: {} vs bound {} (theta90 val {})",
            shadow_s00[i0], bound, shadow_s00[idx90]
        );
    }
}
