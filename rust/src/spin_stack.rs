//! Spin stack engine: stores the active SpinLevels for the focus particle and
//! composes their rotations into a world position + orientation each tick.
//!
//! Per PROJECT_DESIGN §"Spin Stack Engine" and PHYSICS_REFERENCE §2:
//!
//! - Each level's orbit plane is intrinsic to its spin type and FIXED in the
//!   lab frame. Inner levels move the orbit's center (via additive sum) but
//!   do not re-orient its plane — this is the "outside the gyroscopic
//!   influence" property: each spin's amplitude doubles precisely so its
//!   plane is decoupled from the inner level.
//! - Position offset for level k is therefore `lab_rot_k * orbital_start_k *
//!   amp_k`, summed across levels. Only the moving center accumulates.
//! - Orientation is the composition of all levels' rotations — that's how
//!   the visible sphere's pose evolves.
//! - The visible sphere mesh always represents the base photon (radius 1);
//!   higher tiers manifest as the path it traces.

use glam::{DQuat, DVec3};

use crate::types::{level_amplitude, level_spec, SpinType, Tier};

pub struct SpinLevel {
    pub level: u8,
    pub spin_type: SpinType,
    pub tier: Tier,
    pub amplitude: f64,
    /// Angular velocity in natural units. Range [-1.0, 1.0]; sign carries
    /// chirality, magnitude 1.0 = c (saturation, unlocks the next level).
    pub angular_velocity: f64,
    /// Accumulated rotation angle in radians (unbounded — wraps naturally
    /// through quaternion construction).
    pub current_angle: f64,
}

impl SpinLevel {
    pub fn new(level: u8) -> Self {
        let (spin_type, tier) = level_spec(level);
        Self {
            level,
            spin_type,
            tier,
            amplitude: level_amplitude(level),
            angular_velocity: 0.0,
            current_angle: 0.0,
        }
    }

    /// Whether this level has saturated to ±c, allowing the next to activate.
    pub fn is_at_c(&self) -> bool {
        self.angular_velocity.abs() >= 1.0 - 1e-9
    }
}

pub struct SpinStack {
    levels: Vec<SpinLevel>,
}

impl SpinStack {
    /// Level 1 (axial, photon tier) is always present.
    pub fn new() -> Self {
        Self {
            levels: vec![SpinLevel::new(1)],
        }
    }

    pub fn level_count(&self) -> usize {
        self.levels.len()
    }

    pub fn get(&self, level_1based: u8) -> Option<&SpinLevel> {
        if level_1based == 0 {
            None
        } else {
            self.levels.get((level_1based - 1) as usize)
        }
    }

    /// Set the level's angular velocity (clamped to [-1, 1]). Returns false
    /// if the level isn't currently active.
    pub fn set_velocity(&mut self, level_1based: u8, v: f64) -> bool {
        if level_1based == 0 || level_1based as usize > self.levels.len() {
            return false;
        }
        self.levels[(level_1based - 1) as usize].angular_velocity = v.clamp(-1.0, 1.0);
        true
    }

    /// Push the next level onto the stack if the topmost is saturated and
    /// we're not already at level 12. Returns the new level index, or None.
    /// The new level inherits the top's velocity so the 2:1 frequency ratio
    /// between adjacent levels is established immediately.
    pub fn activate_next(&mut self) -> Option<u8> {
        let top = self.levels.last()?;
        if !top.is_at_c() || top.level >= 12 {
            return None;
        }
        let inherited_v = top.angular_velocity;
        let next = top.level + 1;
        let mut new_level = SpinLevel::new(next);
        new_level.angular_velocity = inherited_v;
        self.levels.push(new_level);
        Some(next)
    }

    pub fn can_activate_next(&self) -> bool {
        self.levels
            .last()
            .map_or(false, |l| l.is_at_c() && l.level < 12)
    }

    /// World-position of the outermost level's orbit center in natural units.
    /// This is the sum of orbital contributions from all levels except the
    /// outermost, and represents where the ghost sphere should be centered.
    /// The particle always lies on the ghost sphere's equatorial circle at
    /// distance `effective_radius()` from this center.
    pub fn outer_orbit_center(&self) -> DVec3 {
        if self.levels.len() <= 1 {
            return DVec3::ZERO;
        }
        let mut center = DVec3::ZERO;
        for l in &self.levels[..self.levels.len() - 1] {
            let lab_rot = DQuat::from_axis_angle(l.spin_type.rotation_axis(), l.current_angle);
            let start = l.spin_type.orbital_start();
            if start.length_squared() > 0.0 {
                center += lab_rot * (start * l.amplitude);
            }
        }
        center
    }

    /// Reset to a single idle axial level — equivalent to a fresh stack.
    pub fn reset(&mut self) {
        self.levels.truncate(1);
        if let Some(l) = self.levels.first_mut() {
            l.angular_velocity = 0.0;
            l.current_angle = 0.0;
        }
    }

    /// Advance every active level's angle by `velocity * dt * time_scale / amplitude`.
    /// `time_scale` is in (natural radians) per (real second) at v=c for the
    /// axial level (amplitude = 1). Each orbital level's angular velocity is
    /// ω = v/r = angular_velocity / amplitude, so outer levels rotate slower —
    /// this is what produces spirograph patterns instead of ellipses.
    pub fn advance(&mut self, dt: f64, time_scale: f64) {
        let scaled = dt * time_scale;
        for l in &mut self.levels {
            l.current_angle += l.angular_velocity * scaled / l.amplitude;
        }
    }

    /// Compose the stack into (world position of base particle, orientation).
    pub fn compose(&self) -> (DVec3, DQuat) {
        let mut position = DVec3::ZERO;
        let mut orientation = DQuat::IDENTITY;

        for l in &self.levels {
            let lab_rot = DQuat::from_axis_angle(l.spin_type.rotation_axis(), l.current_angle);

            // Position: each level's orbit plane is fixed in the lab frame,
            // so the offset uses this level's rotation alone. Inner levels
            // only contribute the moving center (via the running sum).
            let start = l.spin_type.orbital_start();
            if start.length_squared() > 0.0 {
                position += lab_rot * (start * l.amplitude);
            }

            // Orientation: pre-multiply so inner spins (added first) apply
            // first to a body vector and outer spins ride on top. Iterating
            // [axial, X, Y, Z] this builds q = R_z * R_y * R_x * R_a, so the
            // visible body's pole only depends on the outer spins above it —
            // axial alone never moves the pole, X tilts it in YZ, etc.
            orientation = lab_rot * orientation;
        }

        (position, orientation)
    }

    /// Effective radius = amplitude of the outermost active level.
    /// At rest (level 1 only) this is 1.0 (the base photon radius).
    pub fn effective_radius(&self) -> f64 {
        self.levels.last().map_or(1.0, |l| l.amplitude)
    }

    /// Spin signature like "+a-x+y+z". Levels with ω=0 are skipped (they
    /// exist in the stack but contribute no chirality yet).
    pub fn signature(&self) -> String {
        let mut s = String::new();
        for l in &self.levels {
            let sign = if l.angular_velocity > 0.0 {
                '+'
            } else if l.angular_velocity < 0.0 {
                '-'
            } else {
                continue;
            };
            s.push(sign);
            s.push(l.spin_type.signature_letter());
        }
        if s.is_empty() {
            "(no spin)".to_string()
        } else {
            s
        }
    }

    /// Classification based on the highest *saturated* level. A particle is
    /// only "at" a tier when that tier's spins are all at ±c.
    pub fn classification(&self) -> &'static str {
        let max_saturated = self
            .levels
            .iter()
            .filter(|l| l.is_at_c())
            .map(|l| l.level)
            .max()
            .unwrap_or(0);
        match max_saturated {
            0 => "B-photon (at rest)",
            1..=3 => "B-photon (sub-photon)",
            4 => "B-photon (full photon)",
            5..=7 => "Sub-electron",
            8 => "Electron",
            9..=11 => "Sub-baryon",
            12 => "Baryon",
            _ => "(unknown)",
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    const EPS: f64 = 1e-9;

    fn approx_zero(v: f64) -> bool {
        v.abs() < EPS
    }

    fn approx_eq(a: f64, b: f64) -> bool {
        (a - b).abs() < EPS
    }

    /// Build a stack with `level` rungs activated. Saturates each rung at +1
    /// in turn so `activate_next` succeeds.
    fn stack_with_levels(level: u8) -> SpinStack {
        let mut s = SpinStack::new();
        while (s.level_count() as u8) < level {
            let top = s.level_count() as u8;
            s.set_velocity(top, 1.0);
            s.activate_next();
        }
        s
    }

    #[test]
    fn axial_only_does_not_translate() {
        let mut s = SpinStack::new();
        s.set_velocity(1, 1.0);
        s.advance(1.0, TAU);
        let (pos, _) = s.compose();
        assert!(approx_zero(pos.x) && approx_zero(pos.y) && approx_zero(pos.z),
            "axial spin must not translate the particle, got {:?}", pos);
    }

    #[test]
    fn x_spin_alone_orbits_in_yz_plane() {
        for k in 0..16u32 {
            let mut s = stack_with_levels(2);
            s.set_velocity(1, 0.0);
            s.set_velocity(2, 1.0);
            s.advance(0.13 * k as f64, TAU);
            let (pos, _) = s.compose();
            assert!(approx_zero(pos.x), "k={}: x={} (expected 0)", k, pos.x);
            assert!(approx_eq(pos.length(), 2.0),
                "k={}: |pos|={} (expected 2)", k, pos.length());
        }
    }

    /// Regression for the "balloon yo-yo" bug: an active axial spin must
    /// NOT re-orient the X-spin's orbit plane. Pre-fix, position composition
    /// used the accumulated quaternion, so axial Y-rotation precessed the
    /// YZ-plane orbit around Y, producing 3D motion. Lab-frame composition
    /// keeps the orbit plane intrinsic to the spin type.
    #[test]
    fn x_spin_with_active_axial_still_orbits_in_yz_plane() {
        for k in 0..16u32 {
            let mut s = stack_with_levels(2);
            s.set_velocity(1, 1.0);
            s.set_velocity(2, 1.0);
            s.advance(0.13 * k as f64, TAU);
            let (pos, _) = s.compose();
            assert!(approx_zero(pos.x), "k={}: x={} (expected 0)", k, pos.x);
            assert!(approx_eq(pos.length(), 2.0),
                "k={}: |pos|={} (expected 2)", k, pos.length());
        }
    }

    /// Regression for the same-frequency ellipse bug: when X (amp=2) and Y
    /// (amp=4) are both at v=c but advancing at the same angular rate, their
    /// sum is an ellipse. With ω = v/r, after dt=1/TAU seconds at time_scale=TAU,
    /// X advances π radians and Y advances π/2 radians.
    /// X orbital_start=Z: R_X(π)*Z*2 = (0,0,-2); R_Y(π/2)*Z*4 = (4,0,0) → total (4,0,-2).
    #[test]
    fn orbital_angular_rate_scales_with_inverse_amplitude() {
        let mut s = stack_with_levels(3);
        s.set_velocity(1, 0.0);
        s.set_velocity(2, 1.0);
        s.set_velocity(3, 1.0);
        // dt=1, time_scale=TAU → θ_X = TAU/2 = π, θ_Y = TAU/4 = π/2
        s.advance(1.0, TAU);
        let (pos, _) = s.compose();
        assert!(approx_eq(pos.x, 4.0),  "pos.x={} (expected 4)",  pos.x);
        assert!(approx_zero(pos.y),     "pos.y={} (expected 0)",  pos.y);
        assert!(approx_eq(pos.z, -2.0), "pos.z={} (expected -2)", pos.z);
    }

    /// Verify that X+Y at saturation produces a spirograph rather than an ellipse.
    /// With θ_X = 2*θ_Y (ω ratio 2:1), dist² = 20 + 32*sin(θ_Y)*cos²(θ_Y),
    /// which ranges from ≈7.68 (dist≈2.77) to ≈32.32 (dist≈5.69).
    /// The same-frequency bug produces dist² = 20 + 8*sin(2θ), ranging only
    /// from ≈3.46 to ≈5.29 — a distinct, tighter band.
    #[test]
    fn x_y_combined_path_is_spirograph_not_ellipse() {
        let mut s = stack_with_levels(3);
        s.set_velocity(1, 0.0);
        s.set_velocity(2, 1.0);
        s.set_velocity(3, 1.0);
        let (mut min_d, mut max_d) = (f64::INFINITY, f64::NEG_INFINITY);
        for _ in 0..400 {
            s.advance(0.02, TAU);
            let (pos, _) = s.compose();
            let d = pos.length();
            if d < min_d { min_d = d; }
            if d > max_d { max_d = d; }
        }
        // Fixed (ω_X=2*ω_Y): max≈5.69, min≈2.77.
        // Buggy (same rate):  max≈5.29, min≈3.46 — both bounds would fail.
        assert!(max_d > 5.5, "max distance {} should approach 5.69", max_d);
        assert!(min_d < 3.2, "min distance {} should approach 2.77", min_d);
    }

    /// Regression for the orientation-wobble bug observed 2026-04-27:
    /// with both axial and X-spin running, the visible body's Y pole was
    /// precessing around the lab Y axis instead of staying in the YZ plane.
    /// Root cause was orientation composition order: `orientation * lab_rot`
    /// applied X first then axial-around-lab-Y to a body vector, dragging
    /// the X-tilted pole around lab Y. Inner spins must apply first to the
    /// body so outer spins ride on top — pre-multiply (`lab_rot * orientation`).
    #[test]
    fn body_pole_stays_in_yz_plane_under_axial_plus_x() {
        for k in 0..32u32 {
            let mut s = stack_with_levels(2);
            s.set_velocity(1, 1.0);
            s.set_velocity(2, 1.0);
            s.advance(0.07 * k as f64, TAU);
            let (_, orient) = s.compose();
            let pole = orient * DVec3::Y;
            assert!(approx_zero(pole.x),
                "k={}: pole.x={} (expected 0; pole={:?})", k, pole.x, pole);
            assert!(approx_eq(pole.length(), 1.0),
                "k={}: pole not unit length: |pole|={}", k, pole.length());
        }
    }
}
