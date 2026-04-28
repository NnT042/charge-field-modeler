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
    pub fn activate_next(&mut self) -> Option<u8> {
        let top = self.levels.last()?;
        if !top.is_at_c() || top.level >= 12 {
            return None;
        }
        let next = top.level + 1;
        self.levels.push(SpinLevel::new(next));
        Some(next)
    }

    pub fn can_activate_next(&self) -> bool {
        self.levels
            .last()
            .map_or(false, |l| l.is_at_c() && l.level < 12)
    }

    /// Advance every active level's angle by `velocity * dt * time_scale`.
    /// `time_scale` is in (natural radians) per (real second) at v=c.
    pub fn advance(&mut self, dt: f64, time_scale: f64) {
        let scaled = dt * time_scale;
        for l in &mut self.levels {
            l.current_angle += l.angular_velocity * scaled;
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
