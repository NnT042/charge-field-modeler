use std::f64::consts::{PI, TAU};

use glam::DVec3;
use godot::classes::{INode3D, Node3D};
use godot::prelude::*;

use crate::path_trace::PathTrace;
use crate::spin_stack::SpinStack;
use crate::types::{level_amplitude, level_spec};
use crate::units;

/// Default trace ring-buffer capacity. At 60Hz physics tick with `min_step`
/// filtering, this gives several seconds of history at modest spin rates and
/// caps GPU work for the line-strip rebuild.
const PATH_TRACE_CAPACITY: usize = 8192;
/// Default minimum movement (natural units) between consecutive recorded
/// samples. The base photon has radius 1, so 0.005 r is fine but not absurd.
const PATH_TRACE_MIN_STEP: f64 = 0.005;
/// Natural-to-world-unit scale factor. The photon mesh has Godot radius 0.5,
/// representing a photon of natural radius 1 — so 1 natural unit = 0.5 world
/// units. All positions returned to GDScript are multiplied by this factor.
const WORLD_SCALE: f32 = 0.5;

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
    path: PathTrace,
    surface_path: PathTrace,
    linear_enabled: bool,
    linear_speed: f64,
    linear_offset: DVec3,
    linear_dir_preset: u8,
    crest_markers: Vec<DVec3>,
    trough_markers: Vec<DVec3>,
    prev_outer_half_period: i64,
}

const LINEAR_DIR_COUNT: u8 = 4;

fn direction_vec(preset: u8) -> DVec3 {
    match preset {
        0 => DVec3::X,
        1 => DVec3::Y,
        2 => DVec3::Z,
        3 => DVec3::new(1.0, 1.0, 1.0).normalize(),
        _ => DVec3::X,
    }
}

fn direction_label(preset: u8) -> &'static str {
    match preset {
        0 => "+X",
        1 => "+Y",
        2 => "+Z",
        3 => "XYZ",
        _ => "+X",
    }
}

#[godot_api]
impl INode3D for FocusParticle {
    fn init(base: Base<Node3D>) -> Self {
        Self {
            base,
            stack: SpinStack::new(),
            time_scale: TAU,
            path: PathTrace::new(PATH_TRACE_CAPACITY, PATH_TRACE_MIN_STEP),
            surface_path: PathTrace::new(PATH_TRACE_CAPACITY, PATH_TRACE_MIN_STEP),
            linear_enabled: false,
            linear_speed: 0.0,
            linear_offset: DVec3::ZERO,
            linear_dir_preset: 0,
            crest_markers: Vec::new(),
            trough_markers: Vec::new(),
            prev_outer_half_period: i64::MIN,
        }
    }

    fn physics_process(&mut self, delta: f64) {
        self.stack.advance(delta, self.time_scale);
        let (local_pos, rot) = self.stack.compose();

        if self.linear_enabled {
            self.linear_offset +=
                direction_vec(self.linear_dir_preset) * self.linear_speed * delta * self.time_scale;
        }

        let world_pos = local_pos + self.linear_offset;
        self.path.record(world_pos, 1.0);

        let pole = rot * DVec3::Y;
        self.surface_path.record(world_pos + pole, 1.0);

        if self.linear_enabled {
            if let Some(outer) = self.stack.outermost_orbital() {
                if outer.angular_velocity.abs() > 1e-12 {
                    let hp = (outer.current_angle / PI).floor() as i64;
                    if self.prev_outer_half_period != i64::MIN
                        && hp != self.prev_outer_half_period
                    {
                        if hp & 1 == 0 {
                            if self.crest_markers.len() >= 256 {
                                self.crest_markers.remove(0);
                            }
                            self.crest_markers.push(world_pos);
                        } else {
                            if self.trough_markers.len() >= 256 {
                                self.trough_markers.remove(0);
                            }
                            self.trough_markers.push(world_pos);
                        }
                    }
                    self.prev_outer_half_period = hp;
                }
            }
        }

        let g_pos = Vector3::new(
            world_pos.x as f32 * WORLD_SCALE,
            world_pos.y as f32 * WORLD_SCALE,
            world_pos.z as f32 * WORLD_SCALE,
        );
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

    /// Path-trace samples in chronological order (oldest -> newest), in world
    /// space (natural units × WORLD_SCALE). Returned as `PackedVector3Array`
    /// so the GDScript renderer can feed it straight to `ImmediateMesh`.
    #[func]
    fn get_path_points(&self) -> PackedVector3Array {
        let pts = self.path.points_chronological();
        let mut arr = PackedVector3Array::new();
        arr.resize(pts.len());
        for (i, p) in pts.iter().enumerate() {
            arr[i] = Vector3::new(
                p.x as f32 * WORLD_SCALE,
                p.y as f32 * WORLD_SCALE,
                p.z as f32 * WORLD_SCALE,
            );
        }
        arr
    }

    /// Per-sample particle radius in world space, same order as `get_path_points`.
    #[func]
    fn get_path_radii(&self) -> PackedFloat32Array {
        let samples = self.path.samples_chronological();
        let mut arr = PackedFloat32Array::new();
        arr.resize(samples.len());
        for (i, s) in samples.iter().enumerate() {
            arr[i] = s.radius as f32 * WORLD_SCALE;
        }
        arr
    }

    #[func]
    fn clear_path(&mut self) {
        self.clear_all_traces();
    }

    #[func]
    fn set_path_enabled(&mut self, e: bool) {
        self.path.set_enabled(e);
        self.surface_path.set_enabled(e);
        if !e {
            self.clear_all_traces();
        }
    }

    #[func]
    fn is_path_enabled(&self) -> bool {
        self.path.enabled()
    }

    #[func]
    fn get_surface_path_points(&self) -> PackedVector3Array {
        let pts = self.surface_path.points_chronological();
        let mut arr = PackedVector3Array::new();
        arr.resize(pts.len());
        for (i, p) in pts.iter().enumerate() {
            arr[i] = Vector3::new(
                p.x as f32 * WORLD_SCALE,
                p.y as f32 * WORLD_SCALE,
                p.z as f32 * WORLD_SCALE,
            );
        }
        arr
    }

    /// World-space center of the outermost level's orbit.
    #[func]
    fn outer_orbit_center(&self) -> Vector3 {
        let c = self.stack.outer_orbit_center() + self.linear_offset;
        Vector3::new(
            c.x as f32 * WORLD_SCALE,
            c.y as f32 * WORLD_SCALE,
            c.z as f32 * WORLD_SCALE,
        )
    }

    #[func]
    fn ghost_sphere_count(&self) -> i32 {
        self.stack.ghost_spheres().len() as i32
    }

    #[func]
    fn ghost_sphere_center(&self, index: i32) -> Vector3 {
        let ghosts = self.stack.ghost_spheres();
        if let Some((c, _, _)) = ghosts.get(index as usize) {
            let world_c = *c + self.linear_offset;
            Vector3::new(
                world_c.x as f32 * WORLD_SCALE,
                world_c.y as f32 * WORLD_SCALE,
                world_c.z as f32 * WORLD_SCALE,
            )
        } else {
            Vector3::ZERO
        }
    }

    #[func]
    fn ghost_sphere_radius(&self, index: i32) -> f32 {
        let ghosts = self.stack.ghost_spheres();
        ghosts
            .get(index as usize)
            .map_or(0.0, |(_, amp, _)| *amp as f32 * WORLD_SCALE)
    }

    #[func]
    fn ghost_sphere_spin_label(&self, index: i32) -> GString {
        let ghosts = self.stack.ghost_spheres();
        ghosts
            .get(index as usize)
            .map_or(GString::default(), |(_, _, st)| GString::from(st.label()))
    }

    #[func]
    fn set_linear_enabled(&mut self, enabled: bool) {
        self.linear_enabled = enabled;
        self.linear_offset = DVec3::ZERO;
        self.clear_all_traces();
    }

    #[func]
    fn is_linear_enabled(&self) -> bool {
        self.linear_enabled
    }

    #[func]
    fn set_linear_speed(&mut self, speed: f64) {
        self.linear_speed = speed.clamp(0.0, 1.0);
    }

    #[func]
    fn get_linear_speed(&self) -> f64 {
        self.linear_speed
    }

    #[func]
    fn set_linear_direction_preset(&mut self, preset: i32) {
        let p = (preset.max(0) as u8).min(LINEAR_DIR_COUNT - 1);
        if p != self.linear_dir_preset {
            self.linear_dir_preset = p;
            self.linear_offset = DVec3::ZERO;
            self.clear_all_traces();
        }
    }

    #[func]
    fn get_linear_direction_preset(&self) -> i32 {
        self.linear_dir_preset as i32
    }

    #[func]
    fn get_linear_direction_label(&self) -> GString {
        GString::from(direction_label(self.linear_dir_preset))
    }

    #[func]
    fn get_linear_direction_count(&self) -> i32 {
        LINEAR_DIR_COUNT as i32
    }

    #[func]
    fn get_linear_offset(&self) -> Vector3 {
        Vector3::new(
            self.linear_offset.x as f32 * WORLD_SCALE,
            self.linear_offset.y as f32 * WORLD_SCALE,
            self.linear_offset.z as f32 * WORLD_SCALE,
        )
    }

    /// Wavelength in natural units: the spatial period of the outermost
    /// orbital wave during linear motion. This is the peak-to-peak distance
    /// visible in the path trace. Returns 0 if no orbital spin or no
    /// linear motion.
    #[func]
    fn get_wavelength(&self) -> f64 {
        if !self.linear_enabled || self.linear_speed < 1e-12 {
            return 0.0;
        }
        match self.stack.outermost_orbital() {
            Some(outer) if outer.angular_velocity.abs() > 1e-12 => {
                self.linear_speed * TAU * outer.orbit_radius() / outer.angular_velocity.abs()
            }
            _ => 0.0,
        }
    }

    /// Characteristic wavelength in natural units at speed c. The spatial
    /// period of the outermost orbital wave. Doesn't require linear
    /// motion — this is the intrinsic wavelength of the spin structure.
    #[func]
    fn get_characteristic_wavelength(&self) -> f64 {
        match self.stack.outermost_orbital() {
            Some(outer) if outer.angular_velocity.abs() > 1e-12 => {
                TAU * outer.orbit_radius() / outer.angular_velocity.abs()
            }
            _ => 0.0,
        }
    }

    /// Wavelength in nanometers. Uses the linear-motion wavelength when
    /// moving, otherwise falls back to the characteristic (speed=c) wavelength.
    #[func]
    fn get_wavelength_nm(&self) -> f64 {
        let wl_nat = self.get_wavelength();
        let wl = if wl_nat > 0.0 {
            wl_nat
        } else {
            self.get_characteristic_wavelength()
        };
        units::wavelength_nm(wl)
    }

    /// EM band label with visible-color detail (e.g. "visible (green)", "near-IR").
    #[func]
    fn get_em_band(&self) -> GString {
        let nm = self.get_wavelength_nm();
        GString::from(units::em_band_detail(nm))
    }

    /// Formatted SI wavelength string with auto-scaled units (e.g. "502.0 nm", "2.00 μm").
    #[func]
    fn get_wavelength_si(&self) -> GString {
        GString::from(units::format_wavelength_si(self.get_wavelength_nm()).as_str())
    }

    /// Approximate spectral color for the current wavelength. Returns a Godot
    /// Color — spectral RGB in visible range, neutral gray outside.
    #[func]
    fn get_wavelength_color(&self) -> Color {
        let (r, g, b) = units::wavelength_to_rgb(self.get_wavelength_nm());
        Color::from_rgb(r, g, b)
    }

    fn clear_all_traces(&mut self) {
        self.path.clear();
        self.surface_path.clear();
        self.crest_markers.clear();
        self.trough_markers.clear();
        self.prev_outer_half_period = i64::MIN;
    }

    #[func]
    fn get_crest_markers(&self) -> PackedVector3Array {
        let mut arr = PackedVector3Array::new();
        arr.resize(self.crest_markers.len());
        for (i, p) in self.crest_markers.iter().enumerate() {
            arr[i] = Vector3::new(
                p.x as f32 * WORLD_SCALE,
                p.y as f32 * WORLD_SCALE,
                p.z as f32 * WORLD_SCALE,
            );
        }
        arr
    }

    #[func]
    fn get_trough_markers(&self) -> PackedVector3Array {
        let mut arr = PackedVector3Array::new();
        arr.resize(self.trough_markers.len());
        for (i, p) in self.trough_markers.iter().enumerate() {
            arr[i] = Vector3::new(
                p.x as f32 * WORLD_SCALE,
                p.y as f32 * WORLD_SCALE,
                p.z as f32 * WORLD_SCALE,
            );
        }
        arr
    }

    /// Reset the spin stack to a single idle axial level and clear the path.
    #[func]
    fn reset_stack(&mut self) {
        self.stack.reset();
        self.clear_all_traces();
        self.linear_enabled = false;
        self.linear_speed = 0.0;
        self.linear_offset = DVec3::ZERO;
        self.linear_dir_preset = 0;
    }
}
