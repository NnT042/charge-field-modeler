use std::f64::consts::TAU;

use godot::classes::{INode3D, Node3D};
use godot::prelude::*;

use crate::spin_stack::SpinStack;
use crate::types::{level_amplitude, level_spec};

/// The focus particle. Owns a `SpinStack` and applies its composed position +
/// orientation to the underlying Node3D each physics tick.
///
/// Time-scale semantics: `time_scale` is "natural radians per real second
/// when ω = c". Default = TAU → one full rotation per real second when level
/// 1 is at ±c (matches the M1 visual feel).
#[derive(GodotClass)]
#[class(base=Node3D)]
pub struct FocusParticle {
    base: Base<Node3D>,
    stack: SpinStack,
    time_scale: f64,
}

#[godot_api]
impl INode3D for FocusParticle {
    fn init(base: Base<Node3D>) -> Self {
        Self {
            base,
            stack: SpinStack::new(),
            time_scale: TAU,
        }
    }

    fn physics_process(&mut self, delta: f64) {
        self.stack.advance(delta, self.time_scale);
        let (pos, rot) = self.stack.compose();

        let g_pos = Vector3::new(pos.x as f32, pos.y as f32, pos.z as f32);
        let g_rot = Quaternion::new(
            rot.x as f32,
            rot.y as f32,
            rot.z as f32,
            rot.w as f32,
        );

        let mut base = self.base_mut();
        base.set_position(g_pos);
        base.set_quaternion(g_rot);
    }
}

#[godot_api]
impl FocusParticle {
    #[func]
    fn level_count(&self) -> i32 {
        self.stack.level_count() as i32
    }

    #[func]
    fn set_level_velocity(&mut self, level: i32, v: f64) -> bool {
        if level <= 0 || level > 12 {
            return false;
        }
        self.stack.set_velocity(level as u8, v)
    }

    #[func]
    fn get_level_velocity(&self, level: i32) -> f64 {
        if level <= 0 {
            return 0.0;
        }
        self.stack
            .get(level as u8)
            .map_or(0.0, |l| l.angular_velocity)
    }

    #[func]
    fn can_activate_next(&self) -> bool {
        self.stack.can_activate_next()
    }

    /// Activate the next level. Returns the new level index, or 0 if not
    /// allowed (top level not at ±c, or already at level 12).
    #[func]
    fn activate_next(&mut self) -> i32 {
        self.stack.activate_next().map(|n| n as i32).unwrap_or(0)
    }

    #[func]
    fn signature(&self) -> GString {
        GString::from(&self.stack.signature())
    }

    #[func]
    fn classification(&self) -> GString {
        GString::from(self.stack.classification())
    }

    /// Effective radius in natural units (r=1 for the base photon).
    #[func]
    fn effective_radius(&self) -> f64 {
        self.stack.effective_radius()
    }

    /// `time_scale` in natural-radians per real-second at v=c.
    /// TAU (~6.283) gives one revolution per second at level-1 saturation.
    #[func]
    fn set_time_scale(&mut self, scale: f64) {
        self.time_scale = scale.max(0.0);
    }

    #[func]
    fn get_time_scale(&self) -> f64 {
        self.time_scale
    }

    /// UI helper: spin type label ("axial" / "x" / "y" / "z") for an
    /// absolute level index 1..=12. Returns "" for out-of-range.
    #[func]
    fn spin_type_label(&self, level: i32) -> GString {
        if level <= 0 || level > 12 {
            return GString::default();
        }
        let (st, _) = level_spec(level as u8);
        GString::from(st.label())
    }

    /// UI helper: tier label ("photon" / "electron" / "baryon").
    #[func]
    fn tier_label(&self, level: i32) -> GString {
        if level <= 0 || level > 12 {
            return GString::default();
        }
        let (_, t) = level_spec(level as u8);
        GString::from(t.label())
    }

    /// UI helper: orbital amplitude in natural units for an absolute level.
    #[func]
    fn level_amplitude(&self, level: i32) -> f64 {
        if level <= 0 || level > 12 {
            return 0.0;
        }
        level_amplitude(level as u8)
    }
}
