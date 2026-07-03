//! Godot-facing wrapper around the pure-Rust physics core.
//!
//! All simulation logic lives in `atom_core.rs` (headless-testable).
//! This node only converts between Godot types and core types — keep it
//! logic-free so tests and the editor always see the same physics.

use crate::atom_core::AtomCore;
use glam::DVec3;
use godot::classes::INode;
use godot::prelude::*;

fn dv(v: Vector3) -> DVec3 {
    DVec3::new(v.x as f64, v.y as f64, v.z as f64)
}

fn gv(v: DVec3) -> Vector3 {
    Vector3::new(v.x as f32, v.y as f32, v.z as f32)
}

fn packed(buf: &[f32]) -> PackedFloat32Array {
    PackedFloat32Array::from(buf)
}

#[derive(GodotClass)]
#[class(base=Node)]
pub struct AtomSim {
    base: Base<Node>,
    core: AtomCore,
}

#[godot_api]
impl INode for AtomSim {
    fn init(base: Base<Node>) -> Self {
        Self {
            base,
            core: AtomCore::new(),
        }
    }
}

#[godot_api]
impl AtomSim {
    // ── Profile management ────────────────────────────────────────────

    #[func]
    fn register_profile(
        &mut self,
        name: GString,
        mass: f64,
        radius: f64,
        csv_data: PackedFloat32Array,
    ) -> i32 {
        let data: Vec<f32> = (0..csv_data.len()).map(|i| csv_data[i]).collect();
        self.core.register_profile(&name.to_string(), mass, radius, &data) as i32
    }

    #[func]
    fn get_profile_count(&self) -> i32 {
        self.core.profiles.len() as i32
    }

    #[func]
    fn get_profile_name(&self, id: i32) -> GString {
        self.core
            .profiles
            .get(id as usize)
            .map(|p| GString::from(&p.name))
            .unwrap_or_default()
    }

    #[func]
    fn get_profile_mass(&self, profile_id: i32) -> f64 {
        self.core
            .profiles
            .get(profile_id as usize)
            .map(|p| p.mass)
            .unwrap_or(0.0)
    }

    // ── Particle management ───────────────────────────────────────────

    #[func]
    fn spawn_particle(
        &mut self,
        profile_id: i32,
        pos: Vector3,
        vel: Vector3,
        pole_dir: Vector3,
    ) -> i32 {
        self.core
            .spawn_particle(profile_id as usize, dv(pos), dv(vel), dv(pole_dir))
            .map(|id| id as i32)
            .unwrap_or(-1)
    }

    #[func]
    fn remove_particle(&mut self, id: i32) {
        self.core.remove_particle(id as usize);
    }

    #[func]
    fn clear_particles(&mut self) {
        self.core.clear_particles();
    }

    #[func]
    fn get_particle_count(&self) -> i32 {
        self.core.particles.len() as i32
    }

    // ── Simulation control ────────────────────────────────────────────

    #[func]
    fn set_running(&mut self, running: bool) {
        self.core.running = running;
    }

    #[func]
    fn is_running(&self) -> bool {
        self.core.running
    }

    #[func]
    fn set_timestep(&mut self, dt: f64) {
        self.core.dt = dt.max(1e-7);
    }

    #[func]
    fn get_timestep(&self) -> f64 {
        self.core.dt
    }

    #[func]
    fn set_gravity_coupling(&mut self, g: f64) {
        self.core.couplings.g_q = g;
    }

    #[func]
    fn get_gravity_coupling(&self) -> f64 {
        self.core.couplings.g_q
    }

    #[func]
    fn set_charge_coupling(&mut self, c: f64) {
        self.core.couplings.c_q = c;
    }

    #[func]
    fn get_charge_coupling(&self) -> f64 {
        self.core.couplings.c_q
    }

    #[func]
    fn set_ambient_gravity(&mut self, g: Vector3) {
        self.core.ambient_gravity = dv(g);
    }

    #[func]
    fn set_ambient_charge(&mut self, c: Vector3) {
        self.core.ambient_charge = dv(c);
    }

    #[func]
    fn set_ambient_pressure(&mut self, p: f64) {
        self.core.couplings.ambient_pressure = p;
    }

    #[func]
    fn get_ambient_pressure(&self) -> f64 {
        self.core.couplings.ambient_pressure
    }

    #[func]
    fn set_torque_coupling(&mut self, t: f64) {
        self.core.couplings.torque = t;
    }

    #[func]
    fn get_torque_coupling(&self) -> f64 {
        self.core.couplings.torque
    }

    #[func]
    fn set_vortex_coupling(&mut self, v: f64) {
        self.core.couplings.vortex = v;
    }

    #[func]
    fn get_vortex_coupling(&self) -> f64 {
        self.core.couplings.vortex
    }

    #[func]
    fn set_drag_coupling(&mut self, d: f64) {
        self.core.couplings.drag = d;
    }

    #[func]
    fn get_drag_coupling(&self) -> f64 {
        self.core.couplings.drag
    }

    #[func]
    fn set_intake_coupling(&mut self, i: f64) {
        self.core.couplings.intake = i;
    }

    #[func]
    fn get_intake_coupling(&self) -> f64 {
        self.core.couplings.intake
    }

    #[func]
    fn set_corot_coupling(&mut self, c: f64) {
        self.core.couplings.corot = c;
    }

    #[func]
    fn get_corot_coupling(&self) -> f64 {
        self.core.couplings.corot
    }

    /// Auto-calibrate gravity coupling for a circular polar orbit.
    /// See AtomCore::auto_calibrate_polar for the derivation.
    #[func]
    fn auto_calibrate_polar(&mut self) -> f64 {
        self.core.auto_calibrate_polar()
    }

    /// Orbit diagnostics for the HUD:
    /// [r, v_tan, v_rad, suggested_G_q, emission_on_axis, theta_deg].
    #[func]
    fn get_orbit_info(&self) -> PackedFloat32Array {
        let info = self.core.orbit_info();
        let out: Vec<f32> = info.iter().map(|&v| v as f32).collect();
        packed(&out)
    }

    /// Force breakdown on the orbiter from the center particle.
    /// [F_gravity, F_charge, F_intake_radial, F_channel, F_vortex_eq,
    ///  A_center, A²_center, speed]. Positive = repulsive, negative = attractive.
    #[func]
    fn get_force_breakdown(&self) -> PackedFloat32Array {
        let fb = self.core.force_breakdown();
        let out: Vec<f32> = fb.iter().map(|&v| v as f32).collect();
        packed(&out)
    }

    #[func]
    fn get_time(&self) -> f64 {
        self.core.time
    }

    // ── Per-particle accessors ────────────────────────────────────────

    #[func]
    fn get_particle_position(&self, id: i32) -> Vector3 {
        self.core
            .particles
            .get(id as usize)
            .map(|p| gv(p.position))
            .unwrap_or(Vector3::ZERO)
    }

    #[func]
    fn get_particle_velocity(&self, id: i32) -> Vector3 {
        self.core
            .particles
            .get(id as usize)
            .map(|p| gv(p.velocity))
            .unwrap_or(Vector3::ZERO)
    }

    #[func]
    fn get_particle_force(&self, id: i32) -> Vector3 {
        self.core
            .particles
            .get(id as usize)
            .map(|p| gv(p.force_accum))
            .unwrap_or(Vector3::ZERO)
    }

    #[func]
    fn get_particle_pole(&self, id: i32) -> Vector3 {
        self.core
            .particles
            .get(id as usize)
            .map(|p| gv(p.pole_axis()))
            .unwrap_or(Vector3::UP)
    }

    #[func]
    fn get_particle_speed(&self, id: i32) -> f64 {
        self.core
            .particles
            .get(id as usize)
            .map(|p| p.velocity.length())
            .unwrap_or(0.0)
    }

    #[func]
    fn set_particle_position(&mut self, id: i32, pos: Vector3) {
        if let Some(p) = self.core.particles.get_mut(id as usize) {
            p.position = dv(pos);
        }
    }

    #[func]
    fn set_particle_velocity(&mut self, id: i32, vel: Vector3) {
        if let Some(p) = self.core.particles.get_mut(id as usize) {
            p.velocity = dv(vel);
        }
    }

    #[func]
    fn set_particle_pole(&mut self, id: i32, pole: Vector3) {
        self.core.set_particle_pole(id as usize, dv(pole));
    }

    // ── Bulk readback for rendering ───────────────────────────────────

    #[func]
    fn get_positions(&self) -> PackedFloat32Array {
        let mut arr = PackedFloat32Array::new();
        arr.resize(self.core.particles.len() * 3);
        for (i, p) in self.core.particles.iter().enumerate() {
            arr[i * 3] = p.position.x as f32;
            arr[i * 3 + 1] = p.position.y as f32;
            arr[i * 3 + 2] = p.position.z as f32;
        }
        arr
    }

    #[func]
    fn get_pole_axes(&self) -> PackedFloat32Array {
        let mut arr = PackedFloat32Array::new();
        arr.resize(self.core.particles.len() * 3);
        for (i, p) in self.core.particles.iter().enumerate() {
            let pole = p.pole_axis();
            arr[i * 3] = pole.x as f32;
            arr[i * 3 + 1] = pole.y as f32;
            arr[i * 3 + 2] = pole.z as f32;
        }
        arr
    }

    #[func]
    fn get_profile_ids(&self) -> PackedInt32Array {
        let mut arr = PackedInt32Array::new();
        arr.resize(self.core.particles.len());
        for (i, p) in self.core.particles.iter().enumerate() {
            arr[i] = p.profile_id as i32;
        }
        arr
    }

    /// Meridional cross-section ring for a particle's emission/absorption
    /// profile. See AtomCore::profile_ring.
    #[func]
    fn get_profile_ring(
        &self,
        particle_id: i32,
        use_absorption: bool,
        base_radius: f32,
        scale: f32,
        num_points: i32,
        plane: i32,
        power: f32,
    ) -> PackedVector3Array {
        let ring = self.core.profile_ring(
            particle_id.max(0) as usize,
            use_absorption,
            base_radius,
            scale,
            num_points.max(0) as usize,
            plane.max(0) as usize,
            power,
        );
        let mut arr = PackedVector3Array::new();
        arr.resize(ring.len());
        for (i, v) in ring.iter().enumerate() {
            arr[i] = gv(*v);
        }
        arr
    }

    /// 16 floats per particle: 12 (3x4 transform) + 4 (RGBA color).
    /// Row-major with interleaved origin (Godot MultiMesh format).
    #[func]
    fn build_multimesh_buffer(&self) -> PackedFloat32Array {
        packed(&self.core.build_multimesh_buffer())
    }

    /// Build transform+color buffer for particles of one profile type only.
    #[func]
    fn build_multimesh_buffer_for_profile(&self, profile_id: i32) -> PackedFloat32Array {
        packed(&self.core.build_multimesh_buffer_for_profile(profile_id.max(0) as usize))
    }

    /// Line buffer for pole axis indicators (2 verts per particle, pos+rgb).
    #[func]
    fn build_pole_indicator_buffer(&self) -> PackedFloat32Array {
        packed(&self.core.build_pole_indicator_buffer())
    }

    /// Surface-of-revolution mesh from a profile's emission table.
    /// [vert_count, idx_count, verts...(x,y,z,nx,ny,nz,r,g,b,a), indices...]
    #[func]
    fn build_profile_mesh(
        &self,
        profile_id: i32,
        lon_segments: i32,
        lat_segments: i32,
    ) -> PackedFloat32Array {
        packed(&self.core.build_profile_mesh(
            profile_id.max(0) as usize,
            lon_segments.max(0) as usize,
            lat_segments.max(0) as usize,
        ))
    }

    #[func]
    fn count_particles_with_profile(&self, profile_id: i32) -> i32 {
        self.core.count_particles_with_profile(profile_id.max(0) as usize) as i32
    }

    // ── Simulation step ───────────────────────────────────────────────

    #[func]
    fn step(&mut self) {
        self.core.step();
    }

    #[func]
    fn step_n(&mut self, steps: i32) {
        self.core.step_n(steps.max(0) as usize);
    }

    // ── Diagnostics ───────────────────────────────────────────────────

    #[func]
    fn get_total_kinetic_energy(&self) -> f64 {
        self.core.total_kinetic_energy()
    }

    #[func]
    fn get_pair_distance(&self, a: i32, b: i32) -> f64 {
        self.core.pair_distance(a.max(0) as usize, b.max(0) as usize)
    }

    // ── VFX: charge emission sprinkler ───────────────────────────────────

    #[func]
    fn advance_vfx(
        &mut self,
        delta: f64,
        emit_per_particle: i32,
        speed: f64,
        lifetime: f64,
    ) -> PackedFloat32Array {
        packed(&self.core.advance_vfx(
            delta,
            emit_per_particle.max(0) as usize,
            speed,
            lifetime,
        ))
    }

    #[func]
    fn get_vfx_count(&self) -> i32 {
        self.core.vfx_count() as i32
    }

    #[func]
    fn set_vfx_enabled(&mut self, enabled: bool) {
        self.core.set_vfx_enabled(enabled);
    }

    #[func]
    fn is_vfx_enabled(&self) -> bool {
        self.core.vfx_enabled
    }
}
