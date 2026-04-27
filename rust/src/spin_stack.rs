//! Spin stack engine: stores the active SpinLevels for the focus particle and
//! composes their nested rotations into a world position + orientation each
//! tick.
//!
//! The core algorithm (per `skills/cfm-dev/SKILL.md`):
//!
//! - For each level top-down, accumulate the level's rotation into a running
//!   quaternion.
//! - For non-axial levels, the orbital offset starts perpendicular to the
//!   rotation axis (per `SpinType::orbital_start`), is scaled by the level's
//!   amplitude, then rotated by the accumulated frame and added to position.
//! - The result: position is the focus particle's center, accumulated quat is
//!   its orientation. The visible sphere mesh always represents the base
//!   photon (radius 1) — the higher tiers manifest as the path it traces.

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
        let mut accumulated = DQuat::IDENTITY;

        for l in &self.levels {
            let rot = DQuat::from_axis_angle(l.spin_type.rotation_axis(), l.current_angle);
            accumulated *= rot;

            let start = l.spin_type.orbital_start();
            if start.length_squared() > 0.0 {
                let offset = accumulated * (start * l.amplitude);
                position += offset;
            }
        }

        (position, accumulated)
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
