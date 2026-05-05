//! Spin stack engine: stores the active SpinLevels for the focus particle and
//! composes their rotations into a world position + orientation each tick.
//!
//! Per PROJECT_DESIGN §"Spin Stack Engine" and PHYSICS_REFERENCE §2:
//!
//! - Position is computed hierarchically: each orbital level shifts the
//!   inner structure outward by `orbital_start * orbit_radius`, then rotates
//!   the whole thing. This means the Y level carries the X ghost (and the
//!   particle inside it) around as a unit — like a ball rolling inside a
//!   larger ball. Orbit radius = amplitude/2 = previous amplitude, so the
//!   inner ghost's far edge just touches the outer ghost's surface.
//! - Orientation is the composition of all levels' rotations — that's how
//!   the visible sphere's pose evolves.
//! - The visible sphere mesh always represents the base photon (radius 1);
//!   higher tiers manifest as the path it traces.

use glam::{DQuat, DVec3};

use crate::types::{level_amplitude, level_spec, SpinType, Tier};

pub struct SpinLevel {
    pub level: u8,
    pub spin_type: SpinType,
    #[allow(dead_code)]
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

    /// True if this level precesses the inner tier around a new axis
    /// without adding orbital offset. The structure stays centered and
    /// rotates as a unit — spinning the fan blades, not swinging them.
    pub fn is_precession(&self) -> bool {
        self.spin_type == SpinType::Axial && self.tier != Tier::Photon
    }

    /// True if this level traces an orbital path (X/Y/Z spins).
    /// Precession levels rotate in place and don't produce ghost spheres.
    pub fn is_orbital(&self) -> bool {
        self.spin_type != SpinType::Axial
    }

    /// Radius used for angular rate scaling: ω_actual = v × time_scale / orbit_radius.
    /// Axial and precession levels use amplitude (larger structure = slower).
    /// Orbital levels use amplitude/2 (the actual orbit radius).
    pub fn orbit_radius(&self) -> f64 {
        if self.spin_type == SpinType::Axial {
            self.amplitude
        } else {
            self.amplitude / 2.0
        }
    }

    /// Rotation axis. Tier 2+ axial precesses around Z (the nozzle void axis
    /// of the tier 1 fan-blade structure). Easy to swap for experimentation.
    pub fn rotation_axis(&self) -> DVec3 {
        if self.is_precession() {
            DVec3::Z
        } else {
            self.spin_type.rotation_axis()
        }
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
    /// The new level starts at ω=0 so the user can bring it up manually.
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

    /// World-position of the outermost ghost sphere's center. In the
    /// hierarchical nesting model the outermost shell is always at the origin.
    pub fn outer_orbit_center(&self) -> DVec3 {
        DVec3::ZERO
    }

    /// Reset to a single idle axial level — equivalent to a fresh stack.
    pub fn reset(&mut self) {
        self.levels.truncate(1);
        if let Some(l) = self.levels.first_mut() {
            l.angular_velocity = 0.0;
            l.current_angle = 0.0;
        }
    }

    /// Advance every active level's angle by `velocity * dt * time_scale / orbit_radius`.
    /// `time_scale` is in (natural radians) per (real second) at v=c for the
    /// axial level (orbit_radius = 1). Each orbital level's angular velocity is
    /// ω = v/r where r = orbit_radius = amplitude/2, so outer levels rotate
    /// slower — this is what produces spirograph patterns instead of ellipses.
    pub fn advance(&mut self, dt: f64, time_scale: f64) {
        let scaled = dt * time_scale;
        for l in &mut self.levels {
            l.current_angle += l.angular_velocity * scaled / l.orbit_radius();
        }
    }

    /// Compose the stack into (world position of base particle, orientation).
    ///
    /// Three kinds of level:
    ///   - Pure axial (tier 1 A): orientation only, position unchanged.
    ///   - Precession (tier 2+ A): rotates accumulated position around a new
    ///     axis without offset — spins the fan blades in place.
    ///   - Orbital (X/Y/Z): offsets outward then rotates — standard nesting.
    pub fn compose(&self) -> (DVec3, DQuat) {
        let mut position = DVec3::ZERO;
        let mut orientation = DQuat::IDENTITY;

        for l in &self.levels {
            let rot_axis = l.rotation_axis();
            let lab_rot = DQuat::from_axis_angle(rot_axis, l.current_angle);

            if l.is_orbital() {
                let start = l.spin_type.orbital_start();
                position = lab_rot * (position + start * l.orbit_radius());
            } else if l.is_precession() {
                position = lab_rot * position;
            }

            orientation = lab_rot * l.spin_type.tumble_alignment() * orientation;
        }

        (position, orientation)
    }

    /// Ghost sphere data for each orbital (X/Y/Z) level:
    /// (center, amplitude, spin_type, tier).
    /// Precession levels (tier 2+ axial) rotate in place and have no ghost.
    /// Center is in natural units, computed by running compose on all levels
    /// ABOVE the ghost's own level (outermost ghost is always at origin).
    pub fn ghost_spheres(&self) -> Vec<(DVec3, f64, SpinType, Tier, DQuat)> {
        let mut result = Vec::new();
        for (idx, level) in self.levels.iter().enumerate() {
            if !level.is_orbital() {
                continue;
            }
            let mut center = DVec3::ZERO;
            let mut orientation = DQuat::IDENTITY;
            for outer in &self.levels[idx + 1..] {
                if outer.is_orbital() {
                    let lab_rot = DQuat::from_axis_angle(
                        outer.spin_type.rotation_axis(),
                        outer.current_angle,
                    );
                    let start = outer.spin_type.orbital_start();
                    center = lab_rot * (center + start * outer.orbit_radius());
                    orientation = lab_rot * orientation;
                } else if outer.is_precession() {
                    let lab_rot = DQuat::from_axis_angle(
                        outer.rotation_axis(),
                        outer.current_angle,
                    );
                    center = lab_rot * center;
                    orientation = lab_rot * orientation;
                }
            }
            result.push((center, level.amplitude, level.spin_type, level.tier, orientation));
        }
        result
    }

    /// The outermost active orbital (non-axial) level with non-zero spin,
    /// or None if no orbital level is spinning. Skips idle (ω≈0) levels so
    /// wavelength reflects the actual spinning structure, not inactive ones.
    pub fn outermost_orbital(&self) -> Option<&SpinLevel> {
        self.levels
            .iter()
            .rev()
            .find(|l| l.is_orbital() && l.angular_velocity.abs() > 1e-12)
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
            assert!(approx_eq(pos.length(), 1.0),
                "k={}: |pos|={} (expected 1, orbit_radius=amp/2=1)", k, pos.length());
        }
    }

    /// Regression for the "balloon yo-yo" bug: an active axial spin must
    /// NOT re-orient the X-spin's orbit plane.
    #[test]
    fn x_spin_with_active_axial_still_orbits_in_yz_plane() {
        for k in 0..16u32 {
            let mut s = stack_with_levels(2);
            s.set_velocity(1, 1.0);
            s.set_velocity(2, 1.0);
            s.advance(0.13 * k as f64, TAU);
            let (pos, _) = s.compose();
            assert!(approx_zero(pos.x), "k={}: x={} (expected 0)", k, pos.x);
            assert!(approx_eq(pos.length(), 1.0),
                "k={}: |pos|={} (expected 1, orbit_radius=amp/2=1)", k, pos.length());
        }
    }

    /// With ω = v/orbit_radius, X (orbit_radius=1) and Y (orbit_radius=2)
    /// advance at different rates. dt=0.5, time_scale=TAU:
    /// θ_X = 0.5*TAU/1 = π, θ_Y = 0.5*TAU/2 = π/2.
    /// Hierarchical: inner = R_X(π)*Z = (0,0,-1), then
    /// R_Y(π/2) * ((0,0,-1) + Z*2) = R_Y(π/2) * (0,0,1) = (1,0,0).
    #[test]
    fn orbital_angular_rate_scales_with_inverse_orbit_radius() {
        let mut s = stack_with_levels(3);
        s.set_velocity(1, 0.0);
        s.set_velocity(2, 1.0);
        s.set_velocity(3, 1.0);
        s.advance(0.5, TAU);
        let (pos, _) = s.compose();
        assert!(approx_eq(pos.x, 1.0),  "pos.x={} (expected 1)",  pos.x);
        assert!(approx_zero(pos.y),     "pos.y={} (expected 0)",  pos.y);
        assert!(approx_zero(pos.z),     "pos.z={} (expected 0)",  pos.z);
    }

    /// X+Y with hierarchical nesting: the path is a sine wave on the Y
    /// ghost's equator. dist² = 5 + 4*cos(2θ_Y), ranging from 1 to 3.
    /// (Rotation preserves length, so dist depends only on the inner offset.)
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
        assert!(max_d > 2.9, "max distance {} should approach 3.0", max_d);
        assert!(min_d < 1.1, "min distance {} should approach 1.0", min_d);
    }

    #[test]
    fn ghost_spheres_x_only_at_origin() {
        let s = stack_with_levels(2);
        let ghosts = s.ghost_spheres();
        assert_eq!(ghosts.len(), 1, "axial+X → one orbital ghost");
        let (center, amp, st, _, _) = &ghosts[0];
        assert!(approx_zero(center.length()), "X ghost center at origin, got {:?}", center);
        assert!(approx_eq(*amp, 2.0), "X amplitude should be 2, got {}", amp);
        assert_eq!(*st, SpinType::X);
    }

    #[test]
    fn ghost_spheres_x_y_positions() {
        let mut s = stack_with_levels(3);
        s.set_velocity(1, 0.0);
        s.set_velocity(2, 1.0);
        s.set_velocity(3, 1.0);
        // θ_Y = 0 → X ghost center = R_Y(0) * (ZERO + Z*2) = Z*2 = (0,0,2)
        let ghosts = s.ghost_spheres();
        assert_eq!(ghosts.len(), 2);
        let (x_center, x_amp, _, _, _) = &ghosts[0];
        let (y_center, y_amp, _, _, _) = &ghosts[1];
        assert!(approx_eq(*x_amp, 2.0));
        assert!(approx_eq(*y_amp, 4.0));
        assert!(approx_zero(y_center.length()), "Y ghost at origin");
        assert!(approx_eq(x_center.z, 2.0), "X ghost at (0,0,2), got {:?}", x_center);
        assert!(approx_zero(x_center.x));
        assert!(approx_zero(x_center.y));

        // Advance so θ_Y = π/2 → X ghost center = R_Y(π/2)*Z*2 = X*2 = (2,0,0)
        s.advance(0.5, TAU);
        let ghosts2 = s.ghost_spheres();
        let (x_center2, _, _, _, _) = &ghosts2[0];
        assert!(approx_eq(x_center2.x, 2.0),
            "after θ_Y=π/2, X ghost at (2,0,0), got {:?}", x_center2);
        assert!(approx_zero(x_center2.y));
        assert!(approx_zero(x_center2.z));
    }

    #[test]
    fn ghost_sphere_inner_edge_touches_outer() {
        let s = stack_with_levels(3);
        let ghosts = s.ghost_spheres();
        let (x_center, x_amp, _, _, _) = &ghosts[0];
        let (y_center, y_amp, _, _, _) = &ghosts[1];
        let dist = (*x_center - *y_center).length();
        assert!(approx_eq(dist + x_amp, *y_amp),
            "inner ghost far edge should touch outer: dist={} + x_amp={} vs y_amp={}",
            dist, x_amp, y_amp);
    }

    /// Regression for the orientation-wobble bug observed 2026-04-27:
    /// with both axial and X-spin running, the visible body's Y pole was
    /// precessing around the lab Y axis instead of staying in the YZ plane.
    /// Root cause was orientation composition order: `orientation * lab_rot`
    /// applied X first then axial-around-lab-Y to a body vector, dragging
    /// the X-tilted pole around lab Y. Inner spins must apply first to the
    /// body so outer spins ride on top — pre-multiply (`lab_rot * orientation`).
    #[test]
    fn outermost_orbital_returns_correct_level() {
        let mut s = stack_with_levels(3);
        s.set_velocity(3, 1.0);
        let outer = s.outermost_orbital().unwrap();
        assert_eq!(outer.spin_type, SpinType::Y);
        assert_eq!(outer.level, 3);
    }

    #[test]
    fn outermost_orbital_skips_idle_levels() {
        let s = stack_with_levels(3);
        // Level 3 activated but at ω=0 — outermost spinning orbital is level 2
        let outer = s.outermost_orbital().unwrap();
        assert_eq!(outer.spin_type, SpinType::X);
        assert_eq!(outer.level, 2);
    }

    #[test]
    fn outermost_orbital_none_for_axial_only() {
        let s = SpinStack::new();
        assert!(s.outermost_orbital().is_none());
    }

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

    #[test]
    fn tier_amplitude_resets_at_boundary() {
        use crate::types::level_amplitude;
        assert_eq!(level_amplitude(4), 8.0);
        assert_eq!(level_amplitude(5), 8.0, "A2 should equal Z1 amplitude");
        assert_eq!(level_amplitude(6), 16.0);
        assert_eq!(level_amplitude(7), 32.0);
        assert_eq!(level_amplitude(8), 64.0, "Z2");
        assert_eq!(level_amplitude(9), 64.0, "A3 should equal Z2 amplitude");
        assert_eq!(level_amplitude(12), 512.0, "Z3");
    }

    #[test]
    fn tier2_axial_translates_particle() {
        let mut s = stack_with_levels(5);
        s.set_velocity(5, 1.0);
        s.advance(0.25, TAU);
        let (pos, _) = s.compose();
        assert!(pos.length() > 0.1,
            "tier 2 axial should orbit, not spin in place; got {:?}", pos);
    }

    #[test]
    fn tier2_axial_no_ghost_sphere() {
        let s = stack_with_levels(5);
        let ghosts = s.ghost_spheres();
        let has_axial_ghost = ghosts.iter().any(|(_, _, st, _, _)| {
            *st == SpinType::Axial
        });
        assert!(!has_axial_ghost,
            "precession levels rotate in place — no ghost sphere");
    }

    #[test]
    fn tier2_axial_precesses_existing_structure() {
        let mut s = stack_with_levels(5);
        for lvl in 1..=4 { s.set_velocity(lvl, 1.0); }
        s.set_velocity(5, 0.0);
        s.advance(0.1, TAU);
        let (pos_before, _) = s.compose();

        s.set_velocity(5, 1.0);
        s.advance(0.25, TAU);
        let (pos_after, _) = s.compose();

        assert!(pos_before.length() > 0.1, "tier 1 structure should have extent");
        assert!((pos_after - pos_before).length() > 0.1,
            "A2 precession should rotate the tier 1 structure to a new position");
    }
}
