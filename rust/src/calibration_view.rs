//! Godot-facing eyeball node for Field Calibration Mode (CM-0 dev scene).
//!
//! Wraps a [`CalibrationParticle`] as a Node3D so a scene can watch the baked
//! loop swing, tumble, and drift as a free rigid body. This is a DISPLAY node —
//! it fakes idle motion and a display-only swing pace so there is something to
//! see before the real collision engine exists. No physics is measured here.
//!
//! Coordinates are normalized so every particle frames the same regardless of
//! its true amplitude (a proton loop spans ~512 natural units, an electron ~64).

use std::f64::consts::TAU;

use glam::DVec3;
use godot::classes::{INode3D, Node3D};
use godot::prelude::*;

use crate::calibration::{AmbientField, CalibrationParticle, FieldPreset};

/// Natural→world unit scale, matching FocusParticle (photon mesh radius 0.5).
const WORLD_SCALE: f64 = 0.5;
/// The baked loop's max radius is normalized to this many natural units so all
/// particles frame similarly. Also the reported `effective_radius` for the
/// camera rig's auto-fit.
const DISPLAY_EXTENT_NATURAL: f64 = 8.0;
/// Full swing revolutions per second at `outer_spin = c`, for display pacing
/// only. The true swing rate is far slower (huge orbit radius); real
/// calibration will set `time_scale` honestly.
const DISPLAY_SWING_REV_PER_SEC: f64 = 0.25;
/// Reference photon flux at Room293K (1× density). Higher/lower presets scale
/// this by `FieldPreset::density_multiple`. A display-pacing magnitude, not yet
/// the SI-calibrated flux (that awaits the natural↔SI density conversion).
const REF_FLUX: f64 = 150.0;

#[derive(GodotClass)]
#[class(base=Node3D)]
pub struct CalibrationView {
    base: Base<Node3D>,
    particle: CalibrationParticle,
    kind: i32,
    /// natural → normalized-natural (before WORLD_SCALE).
    display_scale: f64,
    /// max swept radius of the baked loop, natural units (for scale-aware nudges).
    baked_max_r: f64,
    time_scale: f64,
    paused: bool,
    rng: u64,
    field_on: bool,
    field_balanced: bool,
    field_directional: bool,
    preset: FieldPreset,
    field: AmbientField,
    field_rng: u64,
}

impl CalibrationView {
    fn build(kind: i32) -> (CalibrationParticle, f64, f64, f64) {
        let mut p = match kind {
            1 => CalibrationParticle::neutron(1.0),
            2 => CalibrationParticle::electron(1.0),
            _ => CalibrationParticle::proton(1.0),
        };
        // Idle display motion: a swinging disc + a gentle tumble so the pole wanders.
        p.outer_spin = 0.5;
        p.ang_velocity = DVec3::new(1.0, 0.3, 0.0).normalize() * 0.15;

        let baked_max_r = p
            .hitbox
            .points
            .iter()
            .map(|q| (*q + p.hitbox.swing_offset).length())
            .fold(0.0_f64, f64::max)
            .max(1e-6);
        let display_scale = DISPLAY_EXTENT_NATURAL / baked_max_r;
        // omega = outer_spin·time_scale/swing_orbit_radius; want omega(c)=REV·TAU.
        let time_scale = DISPLAY_SWING_REV_PER_SEC * TAU * p.swing_orbit_radius();
        (p, display_scale, baked_max_r, time_scale)
    }

    /// Display-tuned ambient field: spin_gain scaled by pole inertia so every
    /// particle spins up on a watchable timescale (physically it's trillions of
    /// hits; here we accelerate it for the eyeball). Magnitudes are placeholders.
    fn make_field(
        p: &CalibrationParticle,
        preset: FieldPreset,
        balanced: bool,
        directional: bool,
        ts: f64,
        max_r: f64,
    ) -> AmbientField {
        AmbientField {
            flux: REF_FLUX * preset.density_multiple(),
            photon_fraction: if balanced { 0.5 } else { preset.photon_fraction() },
            direction_bias: if directional { Some(DVec3::X) } else { None },
            momentum: max_r * 0.0005,
            spin_gain: p.i_spin() * 0.005,
            time_scale: ts,
        }
    }

    fn rebuild(&mut self, kind: i32) {
        let (p, ds, mr, ts) = Self::build(kind);
        self.field =
            Self::make_field(&p, self.preset, self.field_balanced, self.field_directional, ts, mr);
        self.particle = p;
        self.kind = kind;
        self.display_scale = ds;
        self.baked_max_r = mr;
        self.time_scale = ts;
    }

    /// xorshift64* — deterministic, no std rng dependency.
    fn next_rand(&mut self) -> f64 {
        let mut x = self.rng;
        x ^= x >> 12;
        x ^= x << 25;
        x ^= x >> 27;
        self.rng = x;
        ((x.wrapping_mul(0x2545F4914F6CDD1D) >> 11) as f64) / ((1u64 << 53) as f64)
    }

    fn rand_unit(&mut self) -> DVec3 {
        let z = self.next_rand() * 2.0 - 1.0;
        let phi = self.next_rand() * TAU;
        let r = (1.0 - z * z).max(0.0).sqrt();
        DVec3::new(r * phi.cos(), r * phi.sin(), z)
    }

    fn to_world(&self, p: DVec3) -> Vector3 {
        let s = self.display_scale * WORLD_SCALE;
        Vector3::new((p.x * s) as f32, (p.y * s) as f32, (p.z * s) as f32)
    }
}

#[godot_api]
impl INode3D for CalibrationView {
    fn init(base: Base<Node3D>) -> Self {
        let (particle, display_scale, baked_max_r, time_scale) = Self::build(0);
        let preset = FieldPreset::Room293K;
        let field = Self::make_field(&particle, preset, false, false, time_scale, baked_max_r);
        Self {
            base,
            particle,
            kind: 0,
            display_scale,
            baked_max_r,
            time_scale,
            paused: false,
            // Fixed non-zero seeds (Math.random-free; determinism is fine here).
            rng: 0x9E3779B97F4A7C15,
            field_on: false,
            field_balanced: false,
            field_directional: false,
            preset,
            field,
            field_rng: 0xD1B54A32D192ED03,
        }
    }

    fn physics_process(&mut self, delta: f64) {
        if !self.paused {
            if self.field_on {
                self.field.tick(&mut self.particle, delta, &mut self.field_rng);
            }
            self.particle.integrate(delta, self.time_scale);
        }
    }
}

#[godot_api]
impl CalibrationView {
    /// Baked loop, live-swung and posed, in world units. Fed to a LINE_STRIP.
    #[func]
    fn get_loop_points(&self) -> PackedVector3Array {
        let pts = self.particle.world_points();
        let mut arr = PackedVector3Array::new();
        arr.resize(pts.len());
        for (i, p) in pts.iter().enumerate() {
            arr[i] = self.to_world(*p);
        }
        arr
    }

    /// A short segment along the world pole axis from the COM, for orientation.
    #[func]
    fn get_pole_segment(&self) -> PackedVector3Array {
        let c = self.particle.position;
        let tip = c + self.particle.pole_axis() * DISPLAY_EXTENT_NATURAL;
        let mut a = PackedVector3Array::new();
        a.resize(2);
        a[0] = self.to_world(c);
        a[1] = self.to_world(tip);
        a
    }

    #[func]
    fn set_kind(&mut self, kind: i32) {
        self.rebuild(kind.clamp(0, 2));
    }

    #[func]
    fn reset(&mut self) {
        let k = self.kind;
        self.rebuild(k);
    }

    #[func]
    fn bump_outer_spin(&mut self, d: f64) {
        self.particle.outer_spin = (self.particle.outer_spin + d).clamp(-1.0, 1.0);
        if self.particle.outer_spin.abs() >= 1.0 {
            self.particle.transmuted = true;
        }
    }

    #[func]
    fn bump_tumble(&mut self, d: f64) {
        // Nudge the tumble rate about the current axis (or seed one if idle).
        let axis = if self.particle.ang_velocity.length_squared() > 1e-12 {
            self.particle.ang_velocity.normalize()
        } else {
            DVec3::new(1.0, 0.3, 0.0).normalize()
        };
        let new_rate = (self.particle.ang_velocity.length() + d).max(0.0);
        self.particle.ang_velocity = axis * new_rate;
    }

    /// Apply one random photon-like impulse, scaled to the particle's size so
    /// drift/tumble/spin changes are visible across proton/neutron/electron.
    #[func]
    fn nudge(&mut self) {
        let pts = self.particle.world_points();
        if pts.is_empty() {
            return;
        }
        let idx = (self.next_rand() * pts.len() as f64) as usize % pts.len();
        let contact = pts[idx];
        let dir = self.rand_unit();
        let mag = self.baked_max_r * 0.05 * self.particle.mass;
        self.particle.apply_impulse(contact, dir * mag);
    }

    #[func]
    fn toggle_pause(&mut self) {
        self.paused = !self.paused;
    }

    #[func]
    fn toggle_field(&mut self) {
        self.field_on = !self.field_on;
    }

    #[func]
    fn toggle_balance(&mut self) {
        self.field_balanced = !self.field_balanced;
        self.field.photon_fraction = if self.field_balanced { 0.5 } else { 2.0 / 3.0 };
    }

    #[func]
    fn toggle_directional(&mut self) {
        self.field_directional = !self.field_directional;
        self.field.direction_bias = if self.field_directional {
            Some(DVec3::X)
        } else {
            None
        };
    }

    #[func]
    fn cycle_preset(&mut self) {
        self.preset = self.preset.cycle();
        self.field.flux = REF_FLUX * self.preset.density_multiple();
        self.field.photon_fraction = if self.field_balanced {
            0.5
        } else {
            self.preset.photon_fraction()
        };
    }

    #[func]
    fn get_status_text(&self) -> GString {
        let name = match self.kind {
            1 => "Neutron  (L11 loop -> L12 orbital swing: kite rake, traps)",
            2 => "Electron (L8 loop -> L9 precession swing)",
            _ => "Proton   (L12 loop -> L13 precession swing: sweeps disc, emits)",
        };
        let drift = self.particle.lin_velocity.length();
        let tumble = self.particle.ang_velocity.length();
        let phi = self.particle.swing_phase / std::f64::consts::PI;
        let pole = self.particle.pole_axis();
        let field_state = if self.field_on {
            format!(
                "ON  mix {:.2} ({})  {}",
                self.field.photon_fraction,
                if self.field_balanced { "balanced" } else { "photon-rich -> spin-up" },
                if self.field_directional { "directional +X" } else { "isotropic" }
            )
        } else {
            "OFF".to_string()
        };
        let preset_line = format!(
            "temp: {} ({:.2}x density)",
            self.preset.label(),
            self.preset.density_multiple()
        );
        let text = format!(
            "{name}\n\
             outer_spin: {:+.3} c    swing phi: {:.2} pi    {}\n\
             tumble: {:.3} rad/s    drift |v|: {:.3}    pole: ({:+.2},{:+.2},{:+.2})\n\
             field: {}    transmuted: {}\n\
             {}\n\
             [Z/X/C] type   [Up/Dn] spin   [Lt/Rt] tumble   [Space] nudge\n\
             [G] field {}   [T] temp   [B] balance   [V] directional   [P] {}   [R] reset   [F] fit",
            self.particle.outer_spin,
            phi,
            if self.paused { "(PAUSED)" } else { "" },
            tumble,
            drift,
            pole.x,
            pole.y,
            pole.z,
            field_state,
            if self.particle.transmuted { "YES (L-unlock)" } else { "no" },
            preset_line,
            if self.field_on { "off" } else { "on" },
            if self.paused { "resume" } else { "pause" },
        );
        GString::from(text.as_str())
    }

    // --- OrbitCameraRig compatibility (it calls these on its focus node) ---

    #[func]
    fn is_linear_enabled(&self) -> bool {
        false
    }

    #[func]
    fn get_linear_offset(&self) -> Vector3 {
        Vector3::ZERO
    }

    #[func]
    fn effective_radius(&self) -> f64 {
        DISPLAY_EXTENT_NATURAL
    }
}
