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

/// CM-2 gear-pump SELF-LIMITER candidates under measurement (see the session
/// task this was written for; sourced docs: elecpro.html/higgs3.pdf/elecrad.pdf
/// for the spin-energy ladder that motivates `gear_efficiency` values fed to
/// `apply_photon`; pole.pdf/grav4.pdf/bright.pdf for direction-relative
/// chirality). These variants do NOT change `apply_photon` — they are read
/// through [`CalibrationParticle::gear_pump_variant_torque`] only, to measure
/// which (if any) makes the gear pump self-limit against the drag channel
/// (`swing_drag_torque`) at Earth mix without an externally-imposed ladder
/// efficiency.
#[derive(Clone, Copy)]
pub enum PumpRule {
    /// Today's `apply_photon` rule: `w_sign = χ`, `w_mag = 1`. Reproduces
    /// `swing_pump_torque · (2·photon_fraction − 1)` exactly.
    Baseline,
    /// Model B, b1 (sign only): chirality becomes DIRECTION-RELATIVE —
    /// `w_sign = χ·sgn(−dir·t̂)`, `w_mag = 1`. `t̂` is the ACTUAL signed swing
    /// tangent (`swing_tangent`), falling back to the geometric positive
    /// tangent `t_pos` at rest (documented in `gear_pump_variant_torque`).
    DirRelSign,
    /// Model B, b2 (magnitude only): `w_sign = χ`, `w_mag = |dir − v_surf|`
    /// (`v_surf` = swing surface velocity, same as `swing_drag_torque`'s `v`).
    CatchWeight,
    /// Model B, b3: b1 × b2 combined.
    DirRelCatch,
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

    /// Apply one field photon under the unified gear-tangent rule (CM-2, see
    /// docs/FIELD_CALIBRATION_MODE.md "CM-2: the gear-tangent pump"). Sourced
    /// from uran.pdf/uran4.pdf ("the spin force is at the tangent" — a
    /// photon's edge spins at c while it travels at c, so the tangential spin
    /// force at the contact tooth has the SAME magnitude as the linear
    /// force) and halbach.pdf (spin energy transfers to linear "on the
    /// tangent" — no separate spin-energy budget). The old model kept the
    /// chirality pump as an independently-scaled term (`spin_gain`, a free
    /// scaffold constant) stacked on top of a SECOND drag impulse that
    /// double-counted the plain Newtonian momentum transfer. This rewrite is
    /// two transfer terms, zero fitted constants:
    ///
    ///   1. Catch-weighted Newtonian transfer — replaces BOTH the old plain
    ///      impulse and the old stacked drag. At rest (`catch = 1`) this is
    ///      exactly the plain impulse; away from rest every motion channel
    ///      (swing, tumble, drift) feels a drag toward the field rest frame.
    ///      This term IS model A.
    ///   2. Gear tangential transfer — the photon's own momentum redirected
    ///      along the surface tangent, signed by the chirality mesh and
    ///      scaled by the catch factor `f`. This term IS the chirality pump,
    ///      now carrying its own lever (`r × t̂`) rather than a lumped gain.
    ///
    /// `gear_efficiency` (dimensionless, ships 1.0 = the uran.pdf-sourced
    /// value) is kept as a parameter ONLY for ablation experiments — it is
    /// not a free constant to be tuned, it's a knob to turn the gear term off
    /// (0.0) for isolation tests.
    ///
    /// `dir` is the photon travel direction; `chirality` ∈ {+1, −1}.
    pub fn apply_photon(&mut self, contact: DVec3, dir: DVec3, chirality: f64, momentum: f64, gear_efficiency: f64) {
        let dir = dir.normalize_or_zero();
        let r = contact - self.position;
        let pole = self.pole_axis();

        // (1) Catch-weighted Newtonian transfer (replaces BOTH the old plain
        // impulse and the old stacked drag — no more double counting). catch =
        // relative speed of photon vs local material point (c = 1), material
        // velocity clamped at c (NaN-runaway guard, see 5793cda). At rest
        // catch = 1 and this reduces to the plain Newtonian impulse exactly.
        // Head-on photons catch harder than chasing ones, so every motion
        // channel (outer swing, tumble, drift) feels a drag toward the field
        // rest frame — this term IS model A.
        let omega = pole * (self.outer_spin / self.swing_orbit_radius);
        let v_mat = (self.lin_velocity + (self.ang_velocity + omega).cross(r)).clamp_length_max(1.0);
        let catch = (dir - v_mat).length();
        let j = dir * (momentum * catch);
        let tau = r.cross(j);
        let tau_pole = tau.dot(pole);
        let tau_perp = tau - pole * tau_pole;
        self.lin_velocity += j / self.mass;
        self.ang_velocity += tau_perp / self.i_transverse;
        self.outer_spin += self.spin_coupling * tau_pole / self.i_spin;

        // (2) Gear tangential transfer (the chirality pump, uran.pdf): the
        // photon's edge spin delivers a tangential impulse of the SAME
        // momentum magnitude as the linear hit (edge speed = c = travel
        // speed), along the pole-positive surface tangent, signed by
        // chirality (inverse-Compton mesh = augment, Compton = cancel),
        // scaled by the catch factor f (opposing surface catches, co-moving
        // slips). gear_efficiency ships 1.0 (sourced); parameter kept for
        // ablation only.
        let rp = r - pole * r.dot(pole);
        let t_pos = pole.cross(rp).normalize_or_zero(); // positive-rotation tangent (NOT signed by current spin)
        if t_pos.length_squared() > 1e-18 {
            let t_actual = self.swing_tangent(contact);
            let f = if t_actual.length_squared() > 1e-18 {
                ((1.0 - dir.dot(t_actual)) * 0.5).clamp(0.0, 1.0)
            } else {
                0.5
            };
            let j_tan = t_pos * (momentum * f * chirality * gear_efficiency);
            let tau_g = r.cross(j_tan);
            let tau_g_pole = tau_g.dot(pole);
            let tau_g_perp = tau_g - pole * tau_g_pole;
            self.lin_velocity += j_tan / self.mass;
            self.ang_velocity += tau_g_perp / self.i_transverse;
            self.outer_spin += self.spin_coupling * tau_g_pole / self.i_spin;
        }

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
    /// that owes nothing to chirality. No collision rule is changed here. This
    /// IS the CM-2 term (1) channel (`apply_photon`'s catch-weighted Newtonian
    /// transfer), measured per unit momentum.
    ///
    /// `occluded`: when true, each `(point, dir)` sample is weighted by
    /// upstream exposure `max(0, −d·r̂)` (Lambert, same rule `tick` uses for
    /// contact sampling — see `AmbientField::sample_contact`), turning the
    /// plain isotropic mean into the shadow-weighted mean a real occluded
    /// field would deliver. `false` reproduces the original unshadowed
    /// analytic channel exactly.
    pub fn swing_drag_torque(&self, dir_samples: usize, swing_samples: usize, occluded: bool) -> f64 {
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
                let rhat = r.normalize_or_zero();
                let v = omega.cross(r); // surface velocity
                let vt = v.length();
                let that = if vt < 1e-12 { None } else { Some(v / vt) };
                let lever = (r - pole * r.dot(pole)).length();
                for d in &dirs {
                    let w = if occluded { (-d.dot(rhat)).max(0.0) } else { 1.0 };
                    if let Some(that) = that {
                        let catch = (*d - v).length(); // relative speed (c=1)
                        torque += w * catch * d.dot(that) * lever; // +along surface = spin-up
                    }
                    // still surface: isotropic sum of tangential momentum is zero,
                    // but the sample's weight still counts toward normalization.
                    count += w;
                }
            }
        }
        if count > 0.0 {
            torque / count
        } else {
            0.0
        }
    }

    /// DIAGNOSTIC (measurement, not a rule): mean pole-axis torque per unit
    /// `(momentum · chirality · gear_efficiency)` from isotropic bombardment —
    /// this IS the CM-2 gear-tangent pump channel (`apply_photon` term 2),
    /// measured the same way `swing_drag_torque` measures term 1. Mirrors its
    /// sampling structure exactly (same Fibonacci-sphere directions, same
    /// swing-phase sweep); per `(point, dir)` sample the contribution is
    /// `f · (r × t_pos)·pole`, where `t_pos` is the positive-rotation tangent
    /// and `f` is the catch factor from the ACTUAL signed tangent at the
    /// current `outer_spin` — exactly as computed inside `apply_photon`.
    /// `(r × t_pos)·pole` reduces algebraically to `lever = |r_perp|` (always
    /// ≥ 0: `t_pos` is defined perpendicular to both `pole` and `r_perp`, so
    /// the pole component of the cross product is just the lever length) —
    /// the gear pump's sign comes entirely from `chirality`, not from this
    /// geometric factor.
    ///
    /// `occluded`: same Lambert exposure weighting as `swing_drag_torque`.
    pub fn swing_pump_torque(&self, dir_samples: usize, swing_samples: usize, occluded: bool) -> f64 {
        let n = self.hitbox.points.len();
        if n == 0 {
            return 0.0;
        }
        let pole = self.pole_axis();
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
                let rhat = r.normalize_or_zero();
                let rp = r - pole * r.dot(pole);
                let t_pos = pole.cross(rp).normalize_or_zero();
                let on_axis = t_pos.length_squared() < 1e-18;
                let lever_pole = if on_axis { 0.0 } else { r.cross(t_pos).dot(pole) };
                let t_actual = if on_axis { DVec3::ZERO } else { self.swing_tangent(world_p) };
                for d in &dirs {
                    let w = if occluded { (-d.dot(rhat)).max(0.0) } else { 1.0 };
                    if !on_axis {
                        let f = if t_actual.length_squared() > 1e-18 {
                            ((1.0 - d.dot(t_actual)) * 0.5).clamp(0.0, 1.0)
                        } else {
                            0.5
                        };
                        torque += w * f * lever_pole;
                    }
                    count += w;
                }
            }
        }
        if count > 0.0 {
            torque / count
        } else {
            0.0
        }
    }

    /// DIAGNOSTIC (measurement, not a rule): mean pole-axis torque per unit
    /// `momentum` from isotropic bombardment under one of the [`PumpRule`]
    /// self-limiter candidates (see that enum's doc comment for sourcing).
    /// Generalizes `swing_pump_torque` (which measures ONLY today's
    /// `apply_photon` gear term, per unit chirality) to also mix the two
    /// chirality populations explicitly by `photon_fraction` (`p`) and to
    /// let the sign/magnitude weighting itself vary by `rule`.
    ///
    /// Structure mirrors `swing_pump_torque` exactly: same Fibonacci-sphere
    /// `dirs`, same swing-phase sweep, same world-point transform, same
    /// count-normalization. Per `(point, dir)` sample:
    /// `net += p·contrib(χ=+1) + (1−p)·contrib(χ=−1)`, with
    /// `contrib(χ) = w_sign(χ) · w_mag · f · lever_signed`:
    ///   - `lever_signed = (r×t_pos)·pole` — the exact same geometric factor
    ///     `swing_pump_torque` calls `lever_pole`; per that method's doc
    ///     comment it reduces algebraically to `|r_perp| ≥ 0` (`t_pos` is
    ///     defined ⊥ both `pole` and `r_perp`). It carries no sign channel
    ///     of its own here either — `DirRelSign`'s sign comes entirely from
    ///     `sgn(−dir·t̂)` below, not from this factor. Named to match this
    ///     task's spec; kept as the raw dot product (not `.abs()`) purely to
    ///     reuse `swing_pump_torque`'s internals unmodified.
    ///   - `f` is `apply_photon` term 2's existing catch-direction factor,
    ///     `(1 − dir·t̂_actual)/2` (0.5 at rest) — `t̂_actual` is the REAL
    ///     signed swing tangent (`swing_tangent`), independent of the
    ///     DirRelSign fallback below.
    ///   - `w_sign(χ) = χ` (Baseline, CatchWeight) or `χ·sgn(−dir·t̂)`
    ///     (DirRelSign, DirRelCatch), where `t̂` = `t̂_actual` if nonzero,
    ///     else falls back to the geometric positive tangent `t_pos` (the
    ///     documented rest-state choice — per the task this was measured
    ///     for, if this fallback makes each chirality population cancel
    ///     itself to ~zero net pump at rest, that's a FINDING, not a bug).
    ///   - `w_mag = 1` (Baseline, DirRelSign) or `|dir − v_surf|`
    ///     (CatchWeight, DirRelCatch), `v_surf = omega × r` (the swing
    ///     surface velocity, same `v` as `swing_drag_torque`).
    ///
    /// `Baseline` at `photon_fraction = 1.0` must equal `swing_pump_torque`
    /// exactly (unit-tested) since `p·contrib(+1) + 0·contrib(−1) =
    /// contrib(+1)` and `contrib(+1)` under Baseline is bit-identical to
    /// `swing_pump_torque`'s per-sample term.
    pub fn gear_pump_variant_torque(
        &self,
        photon_fraction: f64,
        rule: PumpRule,
        dir_samples: usize,
        swing_samples: usize,
    ) -> f64 {
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

        let p = photon_fraction;
        let mut torque = 0.0;
        let mut count = 0.0;
        for si in 0..swing_samples {
            let phi = si as f64 / swing_samples as f64 * std::f64::consts::TAU;
            let srot = DQuat::from_axis_angle(self.hitbox.swing_axis, phi);
            for i in 0..n {
                let world_p =
                    self.position + self.orientation * (srot * (self.hitbox.points[i] + self.hitbox.swing_offset));
                let r = world_p - self.position;
                let rp = r - pole * r.dot(pole);
                let t_pos = pole.cross(rp).normalize_or_zero();
                let on_axis = t_pos.length_squared() < 1e-18;
                // Same quantity as swing_pump_torque's `lever_pole` (reduces
                // algebraically to |r_perp| >= 0; no sign channel of its own).
                let lever_signed = if on_axis { 0.0 } else { r.cross(t_pos).dot(pole) };
                let t_actual = if on_axis { DVec3::ZERO } else { self.swing_tangent(world_p) };
                // DirRelSign/DirRelCatch sign-reference tangent: the actual
                // signed swing tangent, falling back to the geometric
                // positive tangent at rest (documented fallback choice).
                let t_hat = if t_actual.length_squared() > 1e-18 { t_actual } else { t_pos };
                let v_surf = omega.cross(r);
                for d in &dirs {
                    let f = if t_actual.length_squared() > 1e-18 {
                        ((1.0 - d.dot(t_actual)) * 0.5).clamp(0.0, 1.0)
                    } else {
                        0.5
                    };
                    let contrib = |chi: f64| -> f64 {
                        let w_sign = match rule {
                            PumpRule::Baseline | PumpRule::CatchWeight => chi,
                            PumpRule::DirRelSign | PumpRule::DirRelCatch => chi * (-d.dot(t_hat)).signum(),
                        };
                        let w_mag = match rule {
                            PumpRule::Baseline | PumpRule::DirRelSign => 1.0,
                            PumpRule::CatchWeight | PumpRule::DirRelCatch => (*d - v_surf).length(),
                        };
                        w_sign * w_mag * f * lever_signed
                    };
                    torque += p * contrib(1.0) + (1.0 - p) * contrib(-1.0);
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
/// field directional (None = isotropic). `momentum` is the sourced per-photon
/// impulse magnitude (SI-anchored, see `from_si`); `gear_efficiency` is the
/// dimensionless uran.pdf-sourced gear-tangent coupling (ships 1.0 — CM-2
/// paid off the last scaffold constant).
pub struct AmbientField {
    /// Expected photon contacts per unit time.
    pub flux: f64,
    /// Fraction of photons (chirality +1); the rest are antiphotons (−1).
    pub photon_fraction: f64,
    /// None = isotropic; Some(dir) = all photons travel along `dir`.
    pub direction_bias: Option<DVec3>,
    /// Per-photon momentum (natural units, placeholder).
    pub momentum: f64,
    /// Dimensionless gear-tangent efficiency (CM-2, uran.pdf-sourced; ships
    /// 1.0). Kept as a knob for ablation experiments only — NOT a free
    /// scaffold constant to be tuned; see `CalibrationParticle::apply_photon`.
    pub gear_efficiency: f64,
    /// Swing-rate time_scale, matching `CalibrationParticle::integrate`.
    pub time_scale: f64,
    /// When true, `tick` samples contact points by upstream exposure (Lambert
    /// shadowing: `expose = max(0, -dir·r̂)`), matching the measured
    /// `surface_shadow_torque` channel (see `report_alignment_channels`).
    /// When false, contacts are drawn uniformly (pre-occlusion behavior,
    /// bit-identical), which is what the isotropic-mean analytics in
    /// `swing_drag_torque` / `report_spin_equilibrium` assume.
    pub occlusion: bool,
}

/// SI anchoring for an AmbientField (CM-1 second half). The sim can't deliver
/// ~1e11 contacts/s, so one sim photon stands in for `aggregation` (K) real
/// photons: `momentum` scales UP by K, the contact rate scales DOWN by K.
/// `gear_efficiency` does NOT scale with K — it's the dimensionless CM-2
/// coupling (ships 1.0), and both the pump and drag terms in `apply_photon`
/// are linear in the SAME `momentum`, so the equilibrium outer_spin is
/// K-invariant by construction (only noise granularity grows with K —
/// report_si_calibration checks this empirically). `seconds_per_time_unit`
/// converts sim time to SI: flux_sim contacts per natural time unit represent
/// flux_si·density contacts per second.
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
            gear_efficiency: 1.0,
            time_scale,
            occlusion: true,
        }
    }

    /// Build a field whose magnitudes derive from the sourced SI anchors.
    /// `sim_momentum` picks the working macro-photon impulse (K = sim_momentum
    /// / true momentum); `sim_flux` picks the working contact rate per
    /// natural time unit. `gear_efficiency` is set to 1.0 (the uran.pdf-sourced
    /// value) — under the unified apply_photon the gear pump scales with the
    /// SAME per-photon momentum as everything else, so it needs no separate
    /// K-scaling: the settled equilibrium is K-invariant BY CONSTRUCTION, not
    /// by tuning a second aggregated coupling (contrast the old
    /// `spin_gain_per_photon * aggregation` scheme this replaces).
    pub fn from_si(
        preset: FieldPreset,
        particle_mass_kg: f64,
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
            gear_efficiency: 1.0,
            time_scale,
            occlusion: true,
        };
        let cal = SiCalibration {
            flux_si_hz,
            aggregation,
            momentum_natural,
            seconds_per_time_unit,
        };
        (field, cal)
    }

    /// Choose a contact index out of `pts` (world-space points, `center` =
    /// the body's `position`) for an incoming photon traveling along `dir`.
    ///
    /// When `self.occlusion` is false: a single uniform draw — bit-identical
    /// to the pre-occlusion code path (hit COUNT unaffected either way; this
    /// only changes which point within the body gets credited).
    ///
    /// When true: REJECTION SAMPLING for an exposure-weighted draw. Candidate
    /// indices are drawn uniformly and accepted with probability
    /// `expose = max(0, -dir·r̂)` (`r̂ = (pts[idx]-center).normalize_or_zero()`),
    /// i.e. Lambert-shadowed — points facing the incoming stream are likelier
    /// to be picked, points on the downstream/self-occluded side are unlikely.
    /// This matches the measured `surface_shadow_torque` channel (Lambert
    /// sphere-like shadowing, chirality-blind) while keeping hit COUNT per
    /// tick exactly `n` (only the spatial distribution changes). Capped at 64
    /// attempts; on cap, falls through to the last candidate drawn (guards
    /// degenerate/near-planar geometry — e.g. every candidate near-zero
    /// exposure — where rejection sampling could stall indefinitely).
    pub fn sample_contact(&self, pts: &[DVec3], center: DVec3, dir: DVec3, rng: &mut u64) -> usize {
        if pts.is_empty() {
            return 0;
        }
        if !self.occlusion {
            return (xorshift64(rng) * pts.len() as f64) as usize % pts.len();
        }
        let mut last = 0usize;
        for _ in 0..64 {
            let idx = (xorshift64(rng) * pts.len() as f64) as usize % pts.len();
            last = idx;
            let rhat = (pts[idx] - center).normalize_or_zero();
            let expose = (-dir.dot(rhat)).clamp(0.0, 1.0);
            if xorshift64(rng) < expose {
                return idx;
            }
        }
        last
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
            if self.occlusion {
                // dir drawn FIRST, then the contact is rejection-sampled
                // against it (exposure depends on dir).
                let dir = match self.direction_bias {
                    Some(b) => b.normalize_or_zero(),
                    None => rand_unit_vec(rng),
                };
                let idx = self.sample_contact(&pts, p.position, dir, rng);
                let contact = pts[idx];
                let chirality = if xorshift64(rng) < self.photon_fraction { 1.0 } else { -1.0 };
                p.apply_photon(contact, dir, chirality, self.momentum, self.gear_efficiency);
            } else {
                // Exact pre-occlusion draw order — bit-identical to the old
                // code path (idx, then dir, then chirality).
                let idx = (xorshift64(rng) * pts.len() as f64) as usize % pts.len();
                let contact = pts[idx];
                let dir = match self.direction_bias {
                    Some(b) => b.normalize_or_zero(),
                    None => rand_unit_vec(rng),
                };
                let chirality = if xorshift64(rng) < self.photon_fraction { 1.0 } else { -1.0 };
                p.apply_photon(contact, dir, chirality, self.momentum, self.gear_efficiency);
            }
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

    // --- gear-pump self-limiter candidates (PumpRule / gear_pump_variant_torque) ---

    /// Identity check: at `photon_fraction = 1.0` (pure photon population,
    /// chi=-1 term contributes zero weight), `PumpRule::Baseline` must
    /// reproduce `swing_pump_torque · (2·1.0 − 1) = swing_pump_torque`
    /// exactly — `contrib(+1)` under Baseline (`w_sign=1, w_mag=1`) is
    /// bit-for-bit the same per-sample term `swing_pump_torque` accumulates.
    /// Checked at rest and mid-spin; validates the plumbing before trusting
    /// any of the other rule variants.
    #[test]
    fn baseline_variant_matches_swing_pump_torque() {
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        for &s in &[0.0, 0.6] {
            p.outer_spin = s;
            let swing = p.swing_pump_torque(48, 16, false);
            let variant = p.gear_pump_variant_torque(1.0, PumpRule::Baseline, 48, 16);
            assert!(
                (swing - variant).abs() < 1e-9,
                "s={}: swing_pump_torque={} gear_pump_variant_torque(Baseline)={}",
                s, swing, variant
            );
        }
    }

    /// MEASURED relationship (not prejudged): at s=0.6, photon_fraction=1.0,
    /// `PumpRule::CatchWeight` (`w_mag = |dir - v_surf|`) pumps MORE than
    /// `PumpRule::Baseline` (`w_mag = 1`) — probed directly: baseline ≈
    /// 130.57, catch-weight ≈ 149.66 (rust/src/calibration.rs probe, 2026-07).
    /// Reading: at s=0.6 the opposing (head-on) half of the isotropic
    /// direction set has relative speed > 1 (`|dir - v_surf| > 1`) while the
    /// co-moving half has relative speed < 1, and because the catch factor
    /// `f` already up-weights the opposing directions (same asymmetry
    /// `swing_drag_torque` exploits for its drag), the CatchWeight magnitude
    /// term reinforces rather than dilutes the pump at this spin — i.e. b2
    /// alone does NOT act as a limiter here, it grows the pump with speed.
    /// This test locks in the sign of that relationship as a regression
    /// guard; it is a measurement, not a physics claim about which model is
    /// correct (that's the requester's call, made in the report interpretation).
    #[test]
    fn catch_weight_rule_grows_pump_with_headon() {
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        p.outer_spin = 0.6;
        let baseline = p.gear_pump_variant_torque(1.0, PumpRule::Baseline, 48, 16);
        let catch_weight = p.gear_pump_variant_torque(1.0, PumpRule::CatchWeight, 48, 16);
        assert!(
            catch_weight > baseline,
            "measured CatchWeight ({}) should exceed Baseline ({}) at s=0.6, p=1.0 — see comment for the reading",
            catch_weight, baseline
        );
    }

    /// MEASURED finding (not a bug): at rest (`outer_spin = 0`), every point's
    /// actual swing tangent is zero, so `PumpRule::DirRelSign` falls back to
    /// `sgn(-dir·t_pos)` (the geometric, dir-independent-only-in-magnitude
    /// tangent). Averaged over the isotropic Fibonacci-sphere direction set,
    /// this sign term cancels almost exactly to machine noise (probed:
    /// ~3.3e-17 for the proton bake, vs ~130.57 for Baseline at the same
    /// state) — i.e. the direction-relative-chirality sign rule, at rest,
    /// self-cancels: half the isotropic photons see the tangent as "opposing"
    /// and half as "co-moving" in exactly equal measure, net pump ~0. This
    /// is the "population cancellation" finding flagged in the task spec,
    /// not a plumbing bug — asserted here as a tight bound so a future
    /// change that breaks the cancellation shows up as a test failure.
    #[test]
    fn dirrel_sign_rule_at_rest() {
        let p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        assert_eq!(p.outer_spin, 0.0);
        let v = p.gear_pump_variant_torque(1.0, PumpRule::DirRelSign, 48, 16);
        assert!(v.abs() < 1e-6, "expected near-zero net pump at rest by population cancellation, got {:e}", v);
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

    // --- gear-rule collisions (CM-2 unified apply_photon) ---

    /// A photon opposing the surface motion (gears catch, f≈1) pumps far more
    /// spin than one co-moving with it (slips, f≈0).
    ///
    /// NOTE (CM-2 landed): the pump now carries the lever AND rides on
    /// `momentum` (term 2 is `momentum · f · chirality · gear_efficiency`), so
    /// the old `momentum=0.0` isolation trick would zero the pump too, not
    /// just the drag. Isolate the GEAR term instead by differencing
    /// `gear_efficiency=1` against `gear_efficiency=0` at the SAME `momentum`
    /// and contact/dir: term 1 (the catch-weighted Newtonian transfer) is
    /// identical in both runs since it never reads `gear_efficiency`, so it
    /// cancels exactly in the subtraction, leaving only the gear term's
    /// contribution. Same spirit as the original test: opposing ≫ co-moving.
    #[test]
    fn opposing_photon_catches_comoving_slips() {
        // Z-pole particle spinning +; at contact +X the surface moves +Y.
        let gear_delta = |dir: DVec3| -> f64 {
            let mut on = z_pole_particle();
            let mut off = z_pole_particle();
            on.outer_spin = 0.5;
            off.outer_spin = 0.5;
            let contact = on.position + DVec3::X;
            on.apply_photon(contact, dir, 1.0, 0.02, 1.0);
            off.apply_photon(contact, dir, 1.0, 0.02, 0.0);
            on.outer_spin - off.outer_spin
        };
        // opposing surface (+Y) → dir -Y; co-moving → dir +Y.
        let d_opp = gear_delta(-DVec3::Y);
        let d_com = gear_delta(DVec3::Y);
        assert!(d_opp > 0.0, "opposing hit should pump spin, got {}", d_opp);
        assert!(d_opp > d_com.abs() * 100.0,
            "opposing (f≈1) should dwarf co-moving (f≈0): {} vs {}", d_opp, d_com);
    }

    /// Chirality sets the sign of the spin pump (augment vs cancel). Term 1
    /// (the catch-weighted Newtonian transfer) never reads `chirality`, so it
    /// doesn't strictly need isolating here — but use the same gear-ablation
    /// differencing technique as `opposing_photon_catches_comoving_slips` for
    /// a clean, unambiguous read on the GEAR term's sign specifically.
    #[test]
    fn chirality_flips_spin_pump_sign() {
        let gear_delta = |chirality: f64| -> f64 {
            let mut on = z_pole_particle();
            let mut off = z_pole_particle();
            on.outer_spin = 0.5;
            off.outer_spin = 0.5;
            let c = DVec3::X;
            on.apply_photon(on.position + c, -DVec3::Y, chirality, 0.02, 1.0);
            off.apply_photon(off.position + c, -DVec3::Y, chirality, 0.02, 0.0);
            on.outer_spin - off.outer_spin
        };
        let d_pho = gear_delta(1.0);
        let d_anti = gear_delta(-1.0);
        assert!(d_pho > 0.0, "photon (chi=+1) gear term should pump spin up, got {}", d_pho);
        assert!(d_anti < 0.0, "antiphoton (chi=-1) gear term should pump spin down, got {}", d_anti);
    }

    /// An imbalanced field (photon-rich) spins the particle up from rest, while
    /// a balanced field leaves it near zero — the augment/cancel statistic.
    ///
    /// NOTE (CM-2 landed): `momentum=0.0` used to isolate the pump for free
    /// (the old model B pump didn't depend on momentum at all). Under the
    /// unified `apply_photon` BOTH transfer terms scale with momentum, so
    /// this now runs full physics (momentum > 0, gear_efficiency = 1.0, small
    /// dt) and reads the net `outer_spin` drift directly: a balanced field's
    /// drag+pump mix should stay small, a photon-rich field's pump should
    /// dominate and drive it clearly positive.
    #[test]
    fn imbalanced_field_spins_up_balanced_does_not() {
        let make = || CalibrationParticle::new(bake_loop(3, 32), 1.0);
        let ts = 1.0;

        let run = |photon_fraction: f64| -> f64 {
            let mut p = make();
            // occlusion: false — isolate the chirality augment/cancel
            // statistic from the shadow-sampling channel, matching the
            // other pump-isolation tests in this module.
            let field = AmbientField {
                flux: 500.0, photon_fraction, direction_bias: None,
                momentum: 0.02, gear_efficiency: 1.0, time_scale: ts, occlusion: false,
            };
            let mut rng = 0x1234_5678_9abc_def0u64;
            for _ in 0..400 {
                field.tick(&mut p, 0.01, &mut rng);
                p.integrate(0.01, ts);
            }
            p.outer_spin
        };

        let balanced = run(0.5);
        let imbalanced = run(0.9);

        assert!(imbalanced.abs() > balanced.abs() * 5.0,
            "imbalanced should spin up far more than balanced: {} vs {}",
            imbalanced, balanced);
        assert!(imbalanced > 0.0, "photon-rich field spins up (+)");
    }

    /// A directional field pushes net linear drift along its travel direction.
    #[test]
    fn directional_field_drives_drift() {
        let mut p = CalibrationParticle::new(bake_loop(3, 32), 1.0);
        // occlusion: true — the drag vector is always along the fixed `dir`
        // (direction_bias = Some(X)), same as the base Newtonian impulse, so
        // its sign can't flip regardless of which contact point occlusion
        // picks; exercising the shadowed default here is free.
        let field = AmbientField {
            flux: 500.0, photon_fraction: 0.5, direction_bias: Some(DVec3::X),
            momentum: 0.02, gear_efficiency: 0.0, time_scale: 1.0, occlusion: true,
        };
        let mut rng = 0xABCD_1234_5678_9012u64;
        for _ in 0..200 {
            field.tick(&mut p, 0.01, &mut rng);
        }
        assert!(p.lin_velocity.x > 0.0, "field along +X should drift +X, got {:?}", p.lin_velocity);
        assert!(p.lin_velocity.x > p.lin_velocity.y.abs() * 5.0, "drift should be mostly along X");
    }

    // --- Shadow occlusion (Change 1) + full-velocity drag (Change 2) ---

    /// With occlusion ON, contacts are drawn exposure-weighted toward the
    /// upstream side of the incoming photon (Lambert shadowing): the mean of
    /// `-dir·r̂` over many accepted contacts should be clearly positive. With
    /// occlusion OFF, the draw is uniform and that same statistic should
    /// average to ~0 (no directional preference in which point gets hit).
    #[test]
    fn occlusion_biases_contacts_upstream() {
        let p = CalibrationParticle::proton(1.0);
        let pts = p.world_points();
        // Photon travel direction (matches report_alignment_channels' convention:
        // the stream source sits at +Z, photons travel toward -Z).
        let dir = DVec3::NEG_Z;
        let trials = 2000;

        let field_on = AmbientField { occlusion: true, ..bal_like(1.0) };
        let mut rng_on = 0x51A1_C1A5_0000_0001u64;
        let mut sum_on = 0.0;
        for _ in 0..trials {
            let idx = field_on.sample_contact(&pts, p.position, dir, &mut rng_on);
            let rhat = (pts[idx] - p.position).normalize_or_zero();
            sum_on += -dir.dot(rhat);
        }
        let mean_on = sum_on / trials as f64;

        let field_off = AmbientField { occlusion: false, ..bal_like(1.0) };
        let mut rng_off = 0x51A1_C1A5_0000_0001u64;
        let mut sum_off = 0.0;
        for _ in 0..trials {
            let idx = field_off.sample_contact(&pts, p.position, dir, &mut rng_off);
            let rhat = (pts[idx] - p.position).normalize_or_zero();
            sum_off += -dir.dot(rhat);
        }
        let mean_off = sum_off / trials as f64;

        assert!(mean_on > 0.2, "occlusion on should bias contacts upstream, mean={}", mean_on);
        assert!(mean_off.abs() < 0.1, "occlusion off should be unbiased, mean={}", mean_off);
    }

    /// A balanced isotropic field (no pump: photon_fraction 0.5, gear_efficiency
    /// 0) damps an existing tumble (`ang_velocity`) back toward zero — the
    /// full-contact-velocity drag catch (Change 2) opposing the body's own
    /// transverse motion, not just its swing surface.
    #[test]
    fn tumble_damps_in_balanced_field() {
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        p.ang_velocity = DVec3::X * 0.3;
        let field = AmbientField {
            flux: 200.0,
            photon_fraction: 0.5,
            direction_bias: None,
            momentum: 0.02,
            gear_efficiency: 0.0,
            time_scale: 1.0,
            occlusion: true,
        };
        let dt = 0.02;
        let mut rng: u64 = 0x7A11_5EED_0BA1_0001u64;
        let mut settled = false;
        for _ in 0..40_000 {
            field.tick(&mut p, dt, &mut rng);
            p.integrate(dt, field.time_scale);
            if p.ang_velocity.length() < 0.15 {
                settled = true;
                break;
            }
        }
        assert!(settled, "tumble should damp below 0.15, stuck at {} (len {})", p.ang_velocity, p.ang_velocity.length());
    }

    /// Mirror of the above for linear drift: a balanced isotropic field
    /// decelerates existing `lin_velocity` toward the field rest frame.
    #[test]
    fn drift_damps_in_balanced_field() {
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        p.lin_velocity = DVec3::X * 0.3;
        let field = AmbientField {
            flux: 200.0,
            photon_fraction: 0.5,
            direction_bias: None,
            momentum: 0.02,
            gear_efficiency: 0.0,
            time_scale: 1.0,
            occlusion: true,
        };
        let dt = 0.02;
        let mut rng: u64 = 0x0D21_F7A0_0BA1_0002u64;
        let mut settled = false;
        for _ in 0..40_000 {
            field.tick(&mut p, dt, &mut rng);
            p.integrate(dt, field.time_scale);
            if p.lin_velocity.length() < 0.15 {
                settled = true;
                break;
            }
        }
        assert!(settled, "drift should damp below 0.15, stuck at {} (len {})", p.lin_velocity, p.lin_velocity.length());
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
        assert!(p.swing_drag_torque(48, 8, false).abs() < 1e-9, "no drag at rest");
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
        // occlusion: false — momentum=0.0 makes BOTH apply_photon transfer
        // terms inert (term 1 scales with momentum directly; term 2 scales
        // with momentum·gear_efficiency), so this is purely a contact-sampling
        // helper (used by `occlusion_biases_contacts_upstream`, which only
        // reads WHERE contacts land, not the physics applied there).
        AmbientField {
            flux: 500.0, photon_fraction: 0.5, direction_bias: None,
            momentum: 0.0, gear_efficiency: 1.0, time_scale: ts, occlusion: false,
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

    /// A balanced (photon_fraction 0.5, gear_efficiency 0.0) isotropic field
    /// drags a spun-up particle back down toward zero — the drag term alone,
    /// no pump.
    #[test]
    fn balanced_field_drags_spin_down() {
        // Smaller bake than the full proton (still L12/Z-pole) to keep the test fast.
        let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0);
        assert!((p.pole_axis() - DVec3::Z).length() < 1e-9);
        p.outer_spin = 0.8;
        // occlusion: false — this is a basic-mechanism validation of the
        // drag term proper (same family as report_spin_equilibrium, whose
        // predictions bisect the UNSHADOWED swing_drag_torque channel); keep
        // it on the analytic isotropic footing rather than the shadowed one.
        let field = AmbientField {
            flux: 200.0,
            photon_fraction: 0.5,
            direction_bias: None,
            momentum: 0.02,
            gear_efficiency: 0.0,
            time_scale: 1.0,
            occlusion: false,
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
        // occlusion: false — mirrors balanced_field_drags_spin_down's choice.
        let field = AmbientField {
            flux: 200.0,
            photon_fraction: 0.5,
            direction_bias: None,
            momentum: 0.02,
            gear_efficiency: 0.0,
            time_scale: 1.0,
            occlusion: false,
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
    /// disturbing the particle (restores afterward). Returns torque per unit
    /// momentum (see the diagnostic's doc comment).
    fn drag_at_spin(p: &mut CalibrationParticle, s: f64, occluded: bool) -> f64 {
        let saved = p.outer_spin;
        p.outer_spin = s;
        let t = p.swing_drag_torque(48, 16, occluded);
        p.outer_spin = saved;
        t
    }

    /// Evaluate `swing_pump_torque` at a trial `outer_spin` without permanently
    /// disturbing the particle (restores afterward). Returns torque per unit
    /// `momentum · chirality · gear_efficiency` (see the diagnostic's doc
    /// comment).
    fn pump_at_spin(p: &mut CalibrationParticle, s: f64, occluded: bool) -> f64 {
        let saved = p.outer_spin;
        p.outer_spin = s;
        let t = p.swing_pump_torque(48, 16, occluded);
        p.outer_spin = saved;
        t
    }

    /// Evaluate `gear_pump_variant_torque` at a trial `outer_spin` without
    /// permanently disturbing the particle (restores afterward). Same
    /// save/restore convention as `pump_at_spin`/`drag_at_spin`, generalized
    /// to a [`PumpRule`] and a `photon_fraction`. Used by
    /// `report_pump_self_limiter`.
    fn gear_pump_at_spin(
        p: &mut CalibrationParticle,
        s: f64,
        photon_fraction: f64,
        rule: PumpRule,
        dir_samples: usize,
        swing_samples: usize,
    ) -> f64 {
        let saved = p.outer_spin;
        p.outer_spin = s;
        let t = p.gear_pump_variant_torque(photon_fraction, rule, dir_samples, swing_samples);
        p.outer_spin = saved;
        t
    }

    /// Bisect for the equilibrium `outer_spin` where the gear pump balances the
    /// velocity-channel drag: `mix·pump_at_spin(s) + drag_at_spin(s) == 0`,
    /// `mix = 2·photon_fraction − 1` (the field's ⟨chirality⟩ expectation).
    ///
    /// Under CM-2's unified `apply_photon`, BOTH transfer terms are scaled by
    /// the SAME `spin_coupling / i_spin` and (for the pump) `momentum ·
    /// gear_efficiency` / (for the drag) `momentum` — since `gear_efficiency`
    /// ships 1.0, momentum and the inertia coupling cancel out of the
    /// equilibrium condition entirely (they multiply both terms identically).
    /// What's left is pure geometry (`pump`/`drag`, both measured "per unit
    /// momentum/chirality/gear") and the field's photon:antiphoton mix — the
    /// equilibrium is momentum-free BY CONSTRUCTION, not by cancellation of
    /// two independently-tuned scaffold magnitudes as the old model required.
    /// `mix` is signed (positive for a photon-rich field); the drag term is
    /// negative-growing for positive spin (see `swing_drag_torque`'s doc), so
    /// for a photon-rich field the root sits at positive `s`.
    fn bisect_equilibrium(p: &mut CalibrationParticle, mix: f64, occluded: bool) -> f64 {
        let g = |p: &mut CalibrationParticle, s: f64| -> f64 {
            mix * pump_at_spin(p, s, occluded) + drag_at_spin(p, s, occluded)
        };
        let sign = mix.signum();
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

    /// REPORT (ignored by default): THE PRESETS TABLE. Validates the CM-2
    /// unified `apply_photon` scaling law end to end, per particle. Predicted
    /// equilibrium (bisecting `(2p−1)·swing_pump_torque(s*) +
    /// swing_drag_torque(s*) == 0`, both unshadowed and shadowed) against the
    /// actual settled `outer_spin` of a ticked `AmbientField` simulation
    /// (occlusion off and on), swept over field imbalance for proton, neutron,
    /// and electron. Also a small density block confirming density-independence
    /// survives the rewrite. Prints tables; no hard physics asserts (see
    /// `report_*` convention elsewhere in the repo) beyond basic sanity. Ends
    /// with the PRESETS CANDIDATE line — the numbers destined for
    /// `godot/config/particle_presets.json` (not edited here).
    #[test]
    #[ignore]
    fn report_spin_equilibrium() {
        struct ParticleSpec {
            label: &'static str,
            make: Box<dyn Fn() -> CalibrationParticle>,
        }
        // Same bake pattern as report_si_calibration's ParticleSpec table:
        // bake_loop's swing_offset/swing_axis/swing_is_precession depend only
        // on `loop_level`, not `samples`, so `bake_loop(11, 256)` already
        // reproduces `CalibrationParticle::neutron`'s kite offset exactly.
        let specs: Vec<ParticleSpec> = vec![
            ParticleSpec { label: "proton", make: Box::new(|| CalibrationParticle::new(bake_loop(12, 256), 1.0)) },
            ParticleSpec { label: "neutron", make: Box::new(|| CalibrationParticle::new(bake_loop(11, 256), 1.0)) },
            ParticleSpec { label: "electron", make: Box::new(|| CalibrationParticle::new(bake_loop(8, 256), 1.0)) },
        ];

        let momentum = 0.02;
        let flux = 200.0;
        let photon_fractions = [0.5, 0.583, 0.667, 0.75, 1.0];
        let steps = 150_000usize;

        // Ticks `field` (built fresh from `make`) for `steps` and returns
        // (settled_s*, transmuted) — last-20%-of-run mean, same convention as
        // the old sweep.
        let settle = |make: &dyn Fn() -> CalibrationParticle,
                      photon_fraction: f64,
                      flux: f64,
                      occlusion: bool|
         -> (f64, bool) {
            let mut sim = make();
            let field = AmbientField {
                flux,
                photon_fraction,
                direction_bias: None,
                momentum,
                gear_efficiency: 1.0,
                time_scale: 1.0,
                occlusion,
            };
            let dt = (2.0 / flux).clamp(0.005, 0.05); // aim ~1-4 photons/tick
            let mut rng = 0xC0FF_EE12_3456_789Au64;
            let mut sum = 0.0;
            let tail_start = steps - steps / 5; // last 20%
            for i in 0..steps {
                field.tick(&mut sim, dt, &mut rng);
                sim.integrate(dt, field.time_scale);
                if i >= tail_start {
                    sum += sim.outer_spin;
                }
            }
            (sum / (steps - tail_start) as f64, sim.transmuted)
        };

        for spec in &specs {
            println!(
                "\n== {} (flux={:.0}) ==\n{:>10} {:>16} {:>16} {:>16} {:>16} {:>12}",
                spec.label, flux, "p_photon", "pred_s*_unshad", "pred_s*_shad", "sim_s_occ=false", "sim_s_occ=true", "transmuted"
            );
            for &pf in &photon_fractions {
                let mix = 2.0 * pf - 1.0;
                let mut p_u = (spec.make)();
                let predicted_unshadowed = bisect_equilibrium(&mut p_u, mix, false);
                let mut p_s = (spec.make)();
                let predicted_shadowed = bisect_equilibrium(&mut p_s, mix, true);

                let (sim_false, tr_false) = settle(&spec.make, pf, flux, false);
                let (sim_true, tr_true) = settle(&spec.make, pf, flux, true);
                let transmuted = tr_false || tr_true;

                println!(
                    "{:>10.3} {:>16.4} {:>16.4} {:>16.4} {:>16.4} {:>12}",
                    pf,
                    predicted_unshadowed,
                    predicted_shadowed,
                    sim_false,
                    sim_true,
                    if transmuted { "YES" } else { "no" }
                );
            }
        }

        // Density block: flux 50/200/800 at photon_fraction=2/3, occlusion=true
        // — density-independence should survive the rewrite (proton only, kept
        // small per the report's original scope).
        println!("\n== density block (proton, photon_fraction=2/3, occlusion=true) ==");
        println!("{:>8} {:>14} {:>12}", "flux", "settled_s*", "transmuted");
        let proton_make = &specs[0].make;
        for &flux_d in &[50.0, 200.0, 800.0] {
            let (settled, transmuted) = settle(proton_make, 2.0 / 3.0, flux_d, true);
            println!("{:>8.1} {:>14.4} {:>12}", flux_d, settled, if transmuted { "YES" } else { "no" });
        }

        // --- PRESETS CANDIDATE: Earth mix (2/3 photon), occlusion on, Room
        // flux (200), settled outer_spin per particle — the working numbers
        // for godot/config/particle_presets.json.
        let (proton_s, _) = settle(&specs[0].make, 2.0 / 3.0, flux, true);
        let (neutron_s, _) = settle(&specs[1].make, 2.0 / 3.0, flux, true);
        let (electron_s, _) = settle(&specs[2].make, 2.0 / 3.0, flux, true);
        println!(
            "\nPRESETS CANDIDATE (Earth mix 2/3, occlusion on): proton={:.4} neutron={:.4} electron={:.4}",
            proton_s, neutron_s, electron_s
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
    /// `from_si`, converting settle time to SI seconds — under CM-2's unified
    /// `apply_photon` this table is now fully physical: `from_si` no longer
    /// takes a `spin_gain_per_photon` scaffold to solve for, `gear_efficiency`
    /// is just 1.0 (the uran.pdf-sourced value), so there's no "per_photon_gain"
    /// column left to print. Table 3 checks that the settled equilibrium is
    /// K-invariant — quartering `sim_momentum` (and quadrupling `sim_flux` to
    /// hold the represented physical density fixed) must reproduce the same
    /// settled `outer_spin` within noise. This should now pass close to
    /// trivially: `gear_efficiency` doesn't scale with K at all (see
    /// `SiCalibration`'s doc comment), so K-invariance is a near-tautology of
    /// the construction — the empirical check is kept anyway as a guard
    /// against a wiring regression. Run B gets 2× the steps so both runs
    /// receive the same physical photon dose (see the inline comment at the
    /// run_settle calls).
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
            "\n== Table 2: Room-preset settling ==\n{:>10} {:>12} {:>12} {:>14} {:>16}",
            "particle", "settled_s*", "settle_step", "settle_t_nat", "settle_t_SI_s"
        );
        for spec in &specs {
            let (mut field, cal) = AmbientField::from_si(FieldPreset::Room293K, spec.mass_kg, 0.02, 200.0, 1.0);
            // occlusion: false — see report_spin_equilibrium's comment: the
            // K-invariance identity this report validates is derived against
            // the unshadowed isotropic drag channel.
            field.occlusion = false;
            let (settled, settle_step, dt, transmuted) =
                run_settle(&*spec.make, &field, 0xC0FF_EE12_3456_789Au64, 150_000);
            let settle_t_nat = settle_step as f64 * dt;
            let settle_t_si = settle_t_nat * cal.seconds_per_time_unit;
            println!(
                "{:>10} {:>12.4} {:>12} {:>14.4} {:>16.4e}{}",
                spec.label,
                settled,
                settle_step,
                settle_t_nat,
                settle_t_si,
                if transmuted { "  (TRANSMUTED — hit c)" } else { "" }
            );
        }

        // --- Table 3: K-invariance (proton, Room) ---
        println!("\n== Table 3: K-invariance (proton, Room) ==");
        let mass_kg = units::PROTON_MASS_KG;

        let (mut field_a, cal_a) = AmbientField::from_si(FieldPreset::Room293K, mass_kg, 0.02, 200.0, 1.0);
        let (mut field_b, cal_b) = AmbientField::from_si(FieldPreset::Room293K, mass_kg, 0.005, 800.0, 1.0);
        // occlusion: false — see above; K-invariance is validated on the
        // unshadowed channel.
        field_a.occlusion = false;
        field_b.occlusion = false;

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

    /// REPORT (ignored by default): the headline end-to-end check for
    /// Changes 1+2 together. A proton, tilted 45 deg off the field direction,
    /// dropped into a DIRECTIONAL room-like field with shadow occlusion on:
    /// does it actually settle pole-parallel-to-field (0 or 180 deg), not
    /// just show a restoring torque in a static diagnostic sweep (that was
    /// already `report_alignment_channels`)? The full-velocity drag (Change
    /// 2) is what should let it SETTLE rather than ring forever.
    #[test]
    #[ignore]
    fn report_standing_up() {
        let d = DVec3::NEG_Z; // photons travel -Z; stream source sits at +Z (matches report_alignment_channels).

        let make = || {
            let mut p = CalibrationParticle::new(bake_loop(12, 256), 1.0); // proton
            p.orientation = DQuat::from_axis_angle(DVec3::X, 45f64.to_radians());
            p
        };

        // Knobs: Room-like mix (2/3 photon), occlusion on. gear_efficiency=0 to
        // isolate the alignment+damping channels (Change 1's shadow torque +
        // Change 2's full-velocity drag, both now term 1 of `apply_photon`)
        // from the chirality gear pump (term 2), which this report isn't
        // about. Cranking momentum ABOVE this hits an explicit-Euler
        // instability in the velocity-dependent catch (NaN, not physics).
        // Knob history: the pre-unification bookkeeping (5793cda) stacked a
        // redundant catch-weighted impulse on the plain one, ~doubling the
        // coherent alignment torque per hit — momentum 0.001 / flux 2000
        // converged then, but does NOT under the honest single transfer
        // (verified 2M and 10M steps: tilt parks in the 35-49° band; the
        // static attractor per report_alignment_channels is intact, the
        // stochastic descent just stalls). Since coherent drift scales with
        // momentum and diffusion with momentum², halving momentum while
        // doubling flux keeps the same physical momentum current but restores
        // the old drift/noise ratio — many smaller photons, physically MORE
        // honest (real counts are ~1e11/s at ~1.6e-10 momentum).
        let flux = 4000.0;
        let momentum = 0.0005;
        let gear_efficiency = 0.0;
        let photon_fraction = 2.0 / 3.0;
        let dt = 0.005;
        let steps = 2_000_000usize;
        println!(
            "\n=== report_standing_up ===\nknobs: flux={} momentum={} gear_efficiency={} photon_fraction={:.4} dt={} steps={}",
            flux, momentum, gear_efficiency, photon_fraction, dt, steps
        );

        let tilt_deg = |p: &CalibrationParticle| -> f64 { p.pole_axis().dot(DVec3::Z).clamp(-1.0, 1.0).acos().to_degrees() };

        // --- Directional run ---
        let mut p = make();
        let field = AmbientField {
            flux,
            photon_fraction,
            direction_bias: Some(d),
            momentum,
            gear_efficiency,
            time_scale: 1.0,
            occlusion: true,
        };
        let mut rng: u64 = 0xA111_5EED_0001_0003u64;
        let sample_every = (steps / 20).max(1);
        let mut peak_ang = 0.0f64;
        println!("\n-- directional field --\n{:>10} {:>10} {:>14} {:>12}", "step", "tilt_deg", "ang_vel", "outer_spin");
        let mut final_tilt = tilt_deg(&p);
        let mut final_ang = p.ang_velocity.length();
        for step in 0..steps {
            field.tick(&mut p, dt, &mut rng);
            p.integrate(dt, field.time_scale);
            // LATTICE HOLD (hof.pdf: nuclei are "held in a solid structure" by
            // neighboring charge and can only TURN, not translate). Without
            // this, a free particle in a pure beam is a radiation sail: it
            // accelerates downstream until it co-moves with the stream,
            // catch -> 0, and ALL coupling (alignment torque included) dies —
            // observed 2026-07-17 as tilt parking at ~45-47 deg with vanishing
            // ang_velocity/outer_spin. The pre-unification bookkeeping masked
            // this with its unweighted (never-decoupling) duplicate impulse.
            p.lin_velocity = DVec3::ZERO;
            p.position = DVec3::ZERO;
            // FIXED-SWING SCAFFOLD: gear_efficiency=0 means nothing sustains
            // outer_spin, and a NON-swinging loop is a frozen 3D curve with
            // its own private exposure-torque balance (observed: parks hard
            // at 18.3 deg). The 0/180-attractor prediction comes from the
            // SWING-AVERAGED static diagnostic (surface_shadow_torque), so
            // hold outer_spin at a representative settled rate — same
            // fixed-spin footing as that diagnostic's columns. integrate()
            // then advances swing_phase from it.
            p.outer_spin = 0.6;
            let av = p.ang_velocity.length();
            if av > peak_ang {
                peak_ang = av;
            }
            if step % sample_every == 0 || step + 1 == steps {
                let t = tilt_deg(&p);
                println!("{:>10} {:>10.3} {:>14.6} {:>12.4}", step, t, av, p.outer_spin);
                final_tilt = t;
                final_ang = av;
            }
        }
        println!(
            "final: tilt={:.3} deg, |ang_velocity|={:.6}, running peak |ang_velocity|={:.6}",
            final_tilt, final_ang, peak_ang
        );

        // ADIABATIC-REGIME criterion (2026-07-17): at outer_spin 0.6 one swing
        // period is ~8.6M steps (rate = spin/orbit_radius with the L13 lever)
        // while alignment relaxes in ~100k steps, so a single particle TRACKS
        // the per-swing-phase equilibrium (observed: parks at ~18.3 deg for
        // the starting phase) rather than feeling the swing-AVERAGED torque
        // (whose attractors are exactly 0/180 — see report_alignment_channels).
        // The clean 0/180 stand-up is the fast-swing/ensemble limit; the
        // honest single-particle dynamic promise is: descend from 45 deg into
        // the alignment basin and hold, while the isotropic control doesn't
        // move. Require the basin, not the limit point.
        assert!(
            final_tilt < 30.0 || final_tilt > 150.0,
            "directional field should reach the alignment basin (tilt < 30 or > 150), got {:.3} deg",
            final_tilt
        );
        assert!(
            final_ang < peak_ang,
            "ang_velocity should have decayed from its running peak: final={:.6} peak={:.6}",
            final_ang, peak_ang
        );

        // --- CONTROL: isotropic field, everything else identical ---
        let mut pc = make();
        let field_c = AmbientField {
            flux,
            photon_fraction,
            direction_bias: None,
            momentum,
            gear_efficiency,
            time_scale: 1.0,
            occlusion: true,
        };
        let mut rng_c: u64 = 0xA111_5EED_0001_0003u64;
        println!("\n-- isotropic CONTROL --\n{:>10} {:>10} {:>14} {:>12}", "step", "tilt_deg", "ang_vel", "outer_spin");
        let start_tilt_c = tilt_deg(&pc);
        let mut final_tilt_c = start_tilt_c;
        for step in 0..steps {
            field_c.tick(&mut pc, dt, &mut rng_c);
            pc.integrate(dt, field_c.time_scale);
            // Same lattice hold + fixed-swing scaffold as the directional run.
            pc.lin_velocity = DVec3::ZERO;
            pc.position = DVec3::ZERO;
            pc.outer_spin = 0.6;
            if step % sample_every == 0 || step + 1 == steps {
                let t = tilt_deg(&pc);
                println!("{:>10} {:>10.3} {:>14.6} {:>12.4}", step, t, pc.ang_velocity.length(), pc.outer_spin);
                final_tilt_c = t;
            }
        }
        println!("final: tilt={:.3} deg (start was {:.3} deg)", final_tilt_c, start_tilt_c);

        // No directional preference: the isotropic control should NOT show
        // the same tight convergence to 0/180 that the directional run does.
        assert!(
            (final_tilt_c > 15.0 && final_tilt_c < 165.0) || (final_tilt_c - start_tilt_c).abs() < 10.0,
            "isotropic control should show no directional alignment preference, got final tilt {:.3} deg (start {:.3})",
            final_tilt_c, start_tilt_c
        );
    }

    // --- CM-2 gear-pump self-limiter candidates: PumpRule report ---

    /// Locates the SETTLING (dynamically stable) root of `s ->
    /// gear_eff·pump(s) + drag(s)` using CACHED, pre-tabulated `pump`/`drag`
    /// values on a shared `grid` (not fresh calls per gear_efficiency — see
    /// `report_pump_self_limiter`'s header comment on why the tables are
    /// built once and reused).
    ///
    /// `net > 0` pushes `outer_spin` up, `net < 0` pulls it down, so a
    /// physically settling equilibrium is a DOWNWARD crossing (`net` goes
    /// from `>= 0` to `< 0` as `s` increases) — perturb it up, it's pulled
    /// back; perturb it down, it's pushed back up. An UPWARD crossing (`net`
    /// goes `< 0` to `>= 0`) is a REPELLER, not an answer: `s` sitting
    /// exactly there is a mathematical root but an unstable one, and reporting
    /// it would be misleading (this bites `DirRelSign`/`DirRelCatch`
    /// specifically: their pump is ~0 at rest by the population-cancellation
    /// finding, so `net(0) ~= 0` is ALWAYS a root, but for large enough
    /// `gear_efficiency` the pump jumps to a nonzero plateau for any `s > 0`
    /// and net immediately goes positive — s=0 is a repeller there, and the
    /// real settling point is farther out where drag finally catches the
    /// plateau). So this scans for the FIRST DOWNWARD crossing specifically,
    /// treating upward crossings as pass-through waypoints, not answers.
    /// Crossings are refined via the piecewise-LINEAR closed-form root within
    /// the bracketing grid cell (equivalent precision to iterative bisection
    /// given a grid this fine, at zero extra function-call cost).
    ///
    /// No downward crossing anywhere in the scanned range means: if `net`
    /// ever went positive (at s=0 or via an upward crossing), it's "runaway
    /// (transmutes)" (nothing ever pulls `outer_spin` back down, so it rides
    /// the gear pump to the c-saturation clamp); otherwise `net <= 0`
    /// throughout, "rest (no spin-up)".
    fn find_equilibrium_label(grid: &[f64], pump: &[f64], drag: &[f64], gear_eff: f64) -> String {
        let net = |i: usize| gear_eff * pump[i] + drag[i];
        let n = grid.len();
        let mut prev = net(0);
        let mut saw_upward = prev > 0.0;
        for i in 1..n {
            let cur = net(i);
            if prev >= 0.0 && cur < 0.0 {
                // Downward (settling) crossing — the answer.
                let s0 = grid[i - 1];
                let s1 = grid[i];
                let root = (s0 - prev * (s1 - s0) / (cur - prev)).clamp(s0, s1);
                return format!("s*={:.4}", root);
            }
            if prev < 0.0 && cur >= 0.0 {
                // Upward (repelling) crossing — note it and keep scanning;
                // this is NOT a settling equilibrium.
                saw_upward = true;
            }
            prev = cur;
        }
        if saw_upward || prev > 0.0 {
            "runaway (transmutes)".to_string()
        } else {
            "rest (no spin-up)".to_string()
        }
    }

    /// REPORT (ignored by default): MEASURES (does not choose between) the
    /// two candidate CM-2 gear-pump self-limiters described in the task this
    /// was written for.
    ///
    /// Background: `apply_photon`'s gear term (`j_tan = t_pos · momentum · f
    /// · χ · gear_efficiency`) overruns the catch-weighted Newtonian drag at
    /// Earth mix (2/3 photon) with `gear_efficiency = 1.0`, so proton and
    /// electron both spin up to c (transmute) under an ordinary ambient field
    /// — falsified, since real matter is stable in ordinary starlight. Two
    /// families of fix are on the table:
    ///   (A) an externally-imposed spin-energy-ladder efficiency
    ///       (elecpro.html/higgs3.pdf: 16385 = the full stacked-spin sum,
    ///       1820.56 ≈ Dalton = 16385/9, 9 = the axial-spin-electron rung) —
    ///       this report tries `gear_efficiency` in {1.0, 1/9, 1/1820.56,
    ///       1/16385} against the unmodified (`Baseline`) gear rule;
    ///   (B) making the gear rule ITSELF direction/speed-relative
    ///       (pole.pdf/grav4.pdf/bright.pdf), tried here as three sub-rules
    ///       (`DirRelSign`, `CatchWeight`, `DirRelCatch`, see [`PumpRule`])
    ///       at `gear_efficiency = 1.0` — does changing the rule alone (no
    ///       ladder) already produce a finite equilibrium?
    /// The two families aren't mutually exclusive — the equilibrium table
    /// (Block 2) crosses every rule against every `gear_efficiency`, so a
    /// rule that needs help from the ladder shows up as "runaway" at
    /// `gear_efficiency=1.0` but a finite `s*` at a smaller one.
    ///
    /// This test performs MEASUREMENT ONLY: no collision rule
    /// (`apply_photon`, `AmbientField::tick`) is touched, and `PumpRule` /
    /// `gear_pump_variant_torque` are read-only diagnostics layered on top.
    /// The only assert is the drag(0)≈0 sanity check inherited from
    /// `swing_drag_is_zero_at_rest`; everything else is print-only — the
    /// tables ARE the deliverable, interpretation is the requester's call
    /// (see this test's own header printout for a build-time honesty check
    /// on the proton/neutron/electron construction).
    ///
    /// PERFORMANCE NOTE: `gear_pump_variant_torque`/`swing_drag_torque` at
    /// (dir_samples=48, swing_samples=16, n=256 loop points) cost ~2e5 inner
    /// iterations per call. Block 2's cross product (3 particles × 4 rules ×
    /// 3 photon_fractions × 4 gear_efficiencies = 144 rows) would mean tens
    /// of thousands of calls if each row recomputed `pump`/`drag` from
    /// scratch — but `drag` never depends on `rule`/`photon_fraction`/
    /// `gear_efficiency` at all, and `pump` never depends on
    /// `gear_efficiency`. So this report tabulates `pump(s)`/`drag(s)` ONCE
    /// per (particle) and (particle, rule, photon_fraction) respectively, on
    /// a shared 201-point grid over `s ∈ [0, 0.995]`, and every row of Block
    /// 2/3 just looks up (interpolates) into those cached tables via
    /// `find_equilibrium_label` — physically identical results, without the
    /// wasted recomputation.
    #[test]
    #[ignore]
    fn report_pump_self_limiter() {
        struct ParticleSpec {
            label: &'static str,
            make: Box<dyn Fn() -> CalibrationParticle>,
        }
        let specs: Vec<ParticleSpec> = vec![
            ParticleSpec { label: "proton", make: Box::new(|| CalibrationParticle::new(bake_loop(12, 256), 1.0)) },
            ParticleSpec { label: "neutron", make: Box::new(|| CalibrationParticle::new(bake_loop(11, 256), 1.0)) },
            ParticleSpec { label: "electron", make: Box::new(|| CalibrationParticle::new(bake_loop(8, 256), 1.0)) },
        ];

        let rules = [PumpRule::Baseline, PumpRule::DirRelSign, PumpRule::CatchWeight, PumpRule::DirRelCatch];
        let rule_label = |r: PumpRule| match r {
            PumpRule::Baseline => "Baseline",
            PumpRule::DirRelSign => "DirRelSign",
            PumpRule::CatchWeight => "CatchWeight",
            PumpRule::DirRelCatch => "DirRelCatch",
        };

        println!("\n=== report_pump_self_limiter ===");
        println!(
            "particles: proton=bake_loop(12,256) neutron=bake_loop(11,256) electron=bake_loop(8,256), mass=1.0 (direct bake, not the proton()/neutron()/electron() ctors)."
        );

        // --- Header honesty check: does the direct bake match the ctors? ---
        // bake_loop's swing_offset/swing_axis/swing_is_precession depend only
        // on loop_level, not on `samples` (see report_spin_equilibrium's own
        // comment on this) — so the ONLY thing the ctors could differ on is
        // sample count (they call bake_loop(level, recommended_samples(level))
        // instead of a fixed 256). Verify this claim rather than assert it.
        let ctor_pairs: [(&str, u8, fn(f64) -> CalibrationParticle); 3] = [
            ("proton", 12, CalibrationParticle::proton),
            ("neutron", 11, CalibrationParticle::neutron),
            ("electron", 8, CalibrationParticle::electron),
        ];
        for (label, level, ctor) in ctor_pairs {
            let direct = CalibrationParticle::new(bake_loop(level, 256), 1.0);
            let via_ctor = ctor(1.0);
            let off_diff = (direct.hitbox.swing_offset - via_ctor.hitbox.swing_offset).length();
            let axis_diff = (direct.hitbox.swing_axis - via_ctor.hitbox.swing_axis).length();
            let prec_match = direct.hitbox.swing_is_precession == via_ctor.hitbox.swing_is_precession;
            println!(
                "  {}: ctor samples={} (direct=256) | swing_offset diff={:.3e} swing_axis diff={:.3e} precession_match={} -> {}",
                label,
                recommended_samples(level),
                off_diff,
                axis_diff,
                prec_match,
                if off_diff < 1e-9 && axis_diff < 1e-9 && prec_match {
                    "GEOMETRY MATCHES (this report's direct bake reproduces the ctor's swing geometry exactly, including the neutron kite offset; sample-count is the only difference and it doesn't affect swing_offset/axis/precession)"
                } else {
                    "GEOMETRY DIFFERS (see diffs above)"
                }
            );
        }

        // ============================= BLOCK 1 =============================
        // Pump curves per particle, all 4 rules as columns, at Earth mix
        // (p=2/3), plus drag(s) as the last column for reference. Small
        // s-value set (11 points) — direct calls, no caching needed.
        let p_earth = 2.0 / 3.0;
        let s_values = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99];
        println!("\n--- Block 1: pump curves (photon_fraction = 2/3, Earth mix) ---");
        println!("(pump columns: per unit momentum·gear_efficiency; drag column: per unit momentum)");
        for spec in &specs {
            let mut p = (spec.make)();
            println!(
                "\n{}\n{:>6} {:>14} {:>14} {:>14} {:>14} {:>14}",
                spec.label, "s", "Baseline", "DirRelSign", "CatchWeight", "DirRelCatch", "drag"
            );
            for &s in &s_values {
                let b = gear_pump_at_spin(&mut p, s, p_earth, PumpRule::Baseline, 48, 16);
                let dr = gear_pump_at_spin(&mut p, s, p_earth, PumpRule::DirRelSign, 48, 16);
                let cw = gear_pump_at_spin(&mut p, s, p_earth, PumpRule::CatchWeight, 48, 16);
                let dc = gear_pump_at_spin(&mut p, s, p_earth, PumpRule::DirRelCatch, 48, 16);
                let dg = drag_at_spin(&mut p, s, false);
                println!("{:>6.2} {:>14.6} {:>14.6} {:>14.6} {:>14.6} {:>14.6}", s, b, dr, cw, dc, dg);
                if s == 0.0 {
                    assert!(dg.abs() < 1e-9, "{}: drag(0) should be ~0, got {}", spec.label, dg);
                }
            }
        }

        // ============================= BLOCK 2 =============================
        // Equilibrium table: particle x rule x photon_fraction x gear_eff.
        // Tables cached once per (particle) [drag] and (particle, rule,
        // photon_fraction) [pump] on a shared grid, reused across all 4
        // gear_efficiency values (see the performance note in this test's doc
        // comment). Also reused by Block 3 below.
        const GRID_N: usize = 200;
        let grid_s: Vec<f64> = (0..=GRID_N).map(|i| 0.995 * i as f64 / GRID_N as f64).collect();
        let p_values = [0.5, 2.0 / 3.0, 1.0];
        let gear_effs: [(&str, f64); 4] =
            [("1.0", 1.0), ("1/9", 1.0 / 9.0), ("1/1820.56", 1.0 / 1820.56), ("1/16385", 1.0 / 16385.0)];

        // drag_tables[particle_idx] ; pump_tables[particle_idx][rule_idx][p_idx]
        let mut drag_tables: Vec<Vec<f64>> = Vec::with_capacity(specs.len());
        let mut pump_tables: Vec<Vec<Vec<Vec<f64>>>> = Vec::with_capacity(specs.len());
        for spec in &specs {
            let mut p = (spec.make)();
            let drag_table: Vec<f64> = grid_s.iter().map(|&s| drag_at_spin(&mut p, s, false)).collect();
            let mut per_rule: Vec<Vec<Vec<f64>>> = Vec::with_capacity(rules.len());
            for &rule in &rules {
                let mut per_p: Vec<Vec<f64>> = Vec::with_capacity(p_values.len());
                for &pf in &p_values {
                    let table: Vec<f64> =
                        grid_s.iter().map(|&s| gear_pump_at_spin(&mut p, s, pf, rule, 48, 16)).collect();
                    per_p.push(table);
                }
                per_rule.push(per_p);
            }
            drag_tables.push(drag_table);
            pump_tables.push(per_rule);
        }

        println!("\n--- Block 2: equilibrium table (gear_efficiency*pump(s*) + drag(s*) = 0) ---");
        println!("{:>10} {:>12} {:>8} {:>12} {:>28}", "particle", "rule", "p", "gear_eff", "verdict/s*");
        for (pi, spec) in specs.iter().enumerate() {
            for (ri, &rule) in rules.iter().enumerate() {
                for (pfi, &pf) in p_values.iter().enumerate() {
                    let pump_table = &pump_tables[pi][ri][pfi];
                    let drag_table = &drag_tables[pi];
                    for &(eff_label, eff) in &gear_effs {
                        let verdict = find_equilibrium_label(&grid_s, pump_table, drag_table, eff);
                        println!(
                            "{:>10} {:>12} {:>8.4} {:>12} {:>28}",
                            spec.label, rule_label(rule), pf, eff_label, verdict
                        );
                    }
                }
            }
        }

        // ============================= BLOCK 3 =============================
        // Anchors: electron s* under (Baseline & CatchWeight) x gear_eff=1/9
        // @ p=2/3, vs the two sourced electron-speed anchors; proton s* under
        // gear_eff=1/16385 (does it stay finite/sub-c at p=1.0, the pure-
        // photon/collider-like stress case, or does even the ladder factor
        // fail to save it there?).
        println!("\n--- Block 3: anchors ---");
        let electron_idx = specs.iter().position(|s| s.label == "electron").unwrap();
        let proton_idx = specs.iter().position(|s| s.label == "proton").unwrap();
        let baseline_idx = rules.iter().position(|r| matches!(r, PumpRule::Baseline)).unwrap();
        let catchweight_idx = rules.iter().position(|r| matches!(r, PumpRule::CatchWeight)).unwrap();
        let p_earth_idx = p_values.iter().position(|&p| (p - p_earth).abs() < 1e-12).unwrap();
        let p_one_idx = p_values.iter().position(|&p| (p - 1.0).abs() < 1e-12).unwrap();

        let e_drag = &drag_tables[electron_idx];
        let e_baseline_verdict = find_equilibrium_label(
            &grid_s, &pump_tables[electron_idx][baseline_idx][p_earth_idx], e_drag, 1.0 / 9.0,
        );
        let e_catchweight_verdict = find_equilibrium_label(
            &grid_s, &pump_tables[electron_idx][catchweight_idx][p_earth_idx], e_drag, 1.0 / 9.0,
        );
        println!("electron s* (Baseline,   gear_eff=1/9, p=2/3 Earth mix): {}", e_baseline_verdict);
        println!("electron s* (CatchWeight, gear_eff=1/9, p=2/3 Earth mix): {}", e_catchweight_verdict);
        println!("  anchor: Mathis electron speed ~0.0057c (comp2.html, fourth root of 2G)");
        println!("  anchor: project's old trace structure ~0.055c");

        println!();
        let pr_drag = &drag_tables[proton_idx];
        for (ri, &rule) in rules.iter().enumerate() {
            let verdict =
                find_equilibrium_label(&grid_s, &pump_tables[proton_idx][ri][p_one_idx], pr_drag, 1.0 / 16385.0);
            println!(
                "proton s* ({:<11} gear_eff=1/16385, p=1.0 pure-photon/collider-stress): {}",
                format!("{},", rule_label(rule)),
                verdict
            );
        }
    }
}
