use glam::{DVec3, DQuat};
use godot::classes::INode;
use godot::prelude::*;
use std::f64::consts::{FRAC_PI_2, TAU};

// ── Simple PRNG (xorshift64) ─────────────────────────────────────────────

struct Rng {
    state: u64,
}

impl Rng {
    fn new(seed: u64) -> Self {
        Self {
            state: seed.max(1),
        }
    }

    fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

// ── Emission Table ────────────────────────────────────────────────────────
// 91 bins: index 0 = pole (θ=0°), index 90 = equator (θ=90°).
// Bilateral symmetry assumed — both poles identical.

#[derive(Clone)]
struct EmissionTable {
    bins: Vec<f32>,
}

impl EmissionTable {
    fn sample(&self, cos_theta: f64) -> f64 {
        let theta = cos_theta.abs().acos(); // bilateral fold: 0..π/2
        let t = theta / FRAC_PI_2;
        let idx = ((t * (self.bins.len() - 1) as f64) as usize).min(self.bins.len() - 1);
        self.bins[idx] as f64
    }

    fn from_csv_values(csv: &[f32]) -> Self {
        assert!(csv.len() >= 180, "Need at least 180 CSV values (-90° to +89°)");
        let n = csv.len();
        let mut bins = Vec::with_capacity(91);
        for theta in 0..=90u32 {
            let north_idx = n.saturating_sub(1).min((180u32.saturating_sub(theta)) as usize);
            let south_idx = (theta as usize).min(n - 1);
            let north = csv[north_idx];
            let south = csv[south_idx];
            bins.push((north + south) / 2.0);
        }
        let max = bins.iter().copied().fold(0.0f32, f32::max);
        if max > 0.0 {
            for v in &mut bins {
                *v /= max;
            }
        }
        Self { bins }
    }

    fn complement(&self) -> Self {
        Self {
            bins: self.bins.iter().map(|v| 1.0 - v).collect(),
        }
    }
}

// ── Emission CDF (for weighted sampling) ─────────────────────────────────

#[derive(Clone)]
struct EmissionCdf {
    cdf: Vec<f64>,
}

impl EmissionCdf {
    fn from_bins(bins: &[f32]) -> Self {
        let n = bins.len();
        let mut cdf = Vec::with_capacity(n);
        let mut sum = 0.0f64;
        for (i, &val) in bins.iter().enumerate() {
            let theta = FRAC_PI_2 * i as f64 / (n - 1).max(1) as f64;
            sum += val as f64 * theta.sin().max(0.001);
            cdf.push(sum);
        }
        if sum > 0.0 {
            for v in &mut cdf {
                *v /= sum;
            }
        }
        Self { cdf }
    }

    fn sample(&self, u: f64) -> f64 {
        let u = u.clamp(0.0, 0.9999);
        let idx = self.cdf.partition_point(|&v| v < u).min(self.cdf.len() - 1);
        FRAC_PI_2 * idx as f64 / (self.cdf.len() - 1).max(1) as f64
    }
}

// ── Particle Profile ──────────────────────────────────────────────────────

#[derive(Clone)]
struct ParticleProfile {
    name: String,
    mass: f64,
    radius: f64,
    emission: EmissionTable,
    absorption: EmissionTable,
    emission_cdf: EmissionCdf,
}

// ── Sim Particle ──────────────────────────────────────────────────────────

struct SimParticle {
    profile_id: usize,
    position: DVec3,
    velocity: DVec3,
    orientation: DQuat,
    angular_velocity: DVec3,
    force_accum: DVec3,
    torque_accum: DVec3,
}

impl SimParticle {
    fn pole_axis(&self) -> DVec3 {
        self.orientation * DVec3::Y
    }
}

// ── AtomSim GodotClass ───────────────────────────────────────────────────

// ── VFX Particle (charge emission sprinkler) ─────────────────────────────

struct VfxParticle {
    position: DVec3,
    velocity: DVec3,
    age: f64,
    color: (f32, f32, f32),
}

const SOFTENING: f64 = 0.05;
const DAMPING: f64 = 0.9999;
const ANGULAR_DAMPING: f64 = 0.998;
const MIN_RENDER_RADIUS: f32 = 0.15;
const ENVELOPE_MIN_R: f64 = 0.25;
const CONTACT_STIFFNESS: f64 = 100.0;

#[derive(GodotClass)]
#[class(base=Node)]
pub struct AtomSim {
    base: Base<Node>,
    profiles: Vec<ParticleProfile>,
    particles: Vec<SimParticle>,
    g_q: f64,
    c_q: f64,
    ambient_gravity: DVec3,
    ambient_charge: DVec3,
    ambient_pressure: f64,
    torque_coupling: f64,
    vortex_coupling: f64,
    drag_coupling: f64,
    running: bool,
    dt: f64,
    time: f64,
    rng: Rng,
    vfx_particles: Vec<VfxParticle>,
    vfx_enabled: bool,
}

#[godot_api]
impl INode for AtomSim {
    fn init(base: Base<Node>) -> Self {
        Self {
            base,
            profiles: Vec::new(),
            particles: Vec::new(),
            g_q: 1.0,
            c_q: 4.0,
            ambient_gravity: DVec3::ZERO,
            ambient_charge: DVec3::ZERO,
            ambient_pressure: 0.0,
            torque_coupling: 0.1,
            vortex_coupling: 0.1,
            drag_coupling: 0.3,
            running: false,
            dt: 0.0005,
            time: 0.0,
            rng: Rng::new(0xDEAD_BEEF_CAFE),
            vfx_particles: Vec::new(),
            vfx_enabled: true,
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
        let emission = EmissionTable::from_csv_values(&data);
        let absorption = emission.complement();
        let emission_cdf = EmissionCdf::from_bins(&emission.bins);
        let id = self.profiles.len();
        self.profiles.push(ParticleProfile {
            name: name.to_string(),
            mass,
            radius,
            emission,
            absorption,
            emission_cdf,
        });
        id as i32
    }

    #[func]
    fn get_profile_count(&self) -> i32 {
        self.profiles.len() as i32
    }

    #[func]
    fn get_profile_name(&self, id: i32) -> GString {
        self.profiles
            .get(id as usize)
            .map(|p| GString::from(&p.name))
            .unwrap_or_default()
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
        let pid = profile_id as usize;
        if pid >= self.profiles.len() {
            return -1;
        }
        let target = DVec3::new(pole_dir.x as f64, pole_dir.y as f64, pole_dir.z as f64);
        let target = if target.length_squared() > 1e-9 {
            target.normalize()
        } else {
            DVec3::Y
        };
        let orientation = if (target - DVec3::Y).length_squared() < 1e-9 {
            DQuat::IDENTITY
        } else if (target + DVec3::Y).length_squared() < 1e-9 {
            DQuat::from_rotation_z(std::f64::consts::PI)
        } else {
            DQuat::from_rotation_arc(DVec3::Y, target)
        };

        let pole = orientation * DVec3::Y;
        let spin_rate = default_spin_rate(&self.profiles[pid].name);

        let id = self.particles.len();
        self.particles.push(SimParticle {
            profile_id: pid,
            position: DVec3::new(pos.x as f64, pos.y as f64, pos.z as f64),
            velocity: DVec3::new(vel.x as f64, vel.y as f64, vel.z as f64),
            orientation,
            angular_velocity: pole * spin_rate,
            force_accum: DVec3::ZERO,
            torque_accum: DVec3::ZERO,
        });
        id as i32
    }

    #[func]
    fn remove_particle(&mut self, id: i32) {
        let idx = id as usize;
        if idx < self.particles.len() {
            self.particles.swap_remove(idx);
        }
    }

    #[func]
    fn clear_particles(&mut self) {
        self.particles.clear();
        self.vfx_particles.clear();
    }

    #[func]
    fn get_particle_count(&self) -> i32 {
        self.particles.len() as i32
    }

    // ── Simulation control ────────────────────────────────────────────

    #[func]
    fn set_running(&mut self, running: bool) {
        self.running = running;
    }

    #[func]
    fn is_running(&self) -> bool {
        self.running
    }

    #[func]
    fn set_timestep(&mut self, dt: f64) {
        self.dt = dt.max(1e-7);
    }

    #[func]
    fn get_timestep(&self) -> f64 {
        self.dt
    }

    #[func]
    fn set_gravity_coupling(&mut self, g: f64) {
        self.g_q = g;
    }

    #[func]
    fn get_gravity_coupling(&self) -> f64 {
        self.g_q
    }

    #[func]
    fn set_charge_coupling(&mut self, c: f64) {
        self.c_q = c;
    }

    #[func]
    fn get_charge_coupling(&self) -> f64 {
        self.c_q
    }

    #[func]
    fn set_ambient_gravity(&mut self, g: Vector3) {
        self.ambient_gravity = DVec3::new(g.x as f64, g.y as f64, g.z as f64);
    }

    #[func]
    fn set_ambient_charge(&mut self, c: Vector3) {
        self.ambient_charge = DVec3::new(c.x as f64, c.y as f64, c.z as f64);
    }

    #[func]
    fn set_ambient_pressure(&mut self, p: f64) {
        self.ambient_pressure = p;
    }

    #[func]
    fn get_ambient_pressure(&self) -> f64 {
        self.ambient_pressure
    }

    #[func]
    fn set_torque_coupling(&mut self, t: f64) {
        self.torque_coupling = t;
    }

    #[func]
    fn get_torque_coupling(&self) -> f64 {
        self.torque_coupling
    }

    #[func]
    fn set_vortex_coupling(&mut self, v: f64) {
        self.vortex_coupling = v;
    }

    #[func]
    fn get_vortex_coupling(&self) -> f64 {
        self.vortex_coupling
    }

    #[func]
    fn set_drag_coupling(&mut self, d: f64) {
        self.drag_coupling = d;
    }

    #[func]
    fn get_drag_coupling(&self) -> f64 {
        self.drag_coupling
    }

    /// Auto-calibrate gravity coupling for a circular polar orbit.
    ///
    /// Finds the heaviest particle (the "center") and the lightest (the
    /// "orbiter"), measures their separation r and the orbiter's tangential
    /// speed v_tan, then sets G_q so gravity alone supplies the centripetal
    /// acceleration: a_c = v_tan²/r = G_q·m_center/r²  ⇒  G_q = v_tan²·r / m_center.
    ///
    /// Near the pole the charge force ≈ 0 (emission vanishes on-axis), so the
    /// orbiter is held almost purely by gravity — this lands it in the right
    /// dartboard. Returns the new G_q (0.0 if it couldn't calibrate).
    #[func]
    fn auto_calibrate_polar(&mut self) -> f64 {
        let Some((center, orbiter)) = self.center_and_orbiter() else {
            return 0.0;
        };
        let m_center = self.profiles[self.particles[center].profile_id].mass;
        let d = self.particles[orbiter].position - self.particles[center].position;
        let r = d.length();
        if r < 1e-9 || m_center < 1e-12 {
            return 0.0;
        }
        let d_hat = d / r;
        let v_rel = self.particles[orbiter].velocity - self.particles[center].velocity;
        let v_tan = (v_rel - d_hat * v_rel.dot(d_hat)).length();
        let g = v_tan * v_tan * r / m_center;
        self.g_q = g;
        g
    }

    /// Orbit diagnostics for the HUD:
    /// [r, v_tan, v_rad, suggested_G_q, emission_on_axis].
    /// emission_on_axis is the center's charge emission sampled toward the
    /// orbiter — ≈0 means the orbiter sits in the field-free polar channel.
    #[func]
    fn get_orbit_info(&self) -> PackedFloat32Array {
        let Some((center, orbiter)) = self.center_and_orbiter() else {
            return PackedFloat32Array::from(&[0.0f32; 5][..]);
        };
        let cp = &self.particles[center];
        let op = &self.particles[orbiter];
        let m_center = self.profiles[cp.profile_id].mass;
        let d = op.position - cp.position;
        let r = d.length().max(1e-9);
        let d_hat = d / r;
        let v_rel = op.velocity - cp.velocity;
        let v_rad = v_rel.dot(d_hat);
        let v_tan = (v_rel - d_hat * v_rad).length();
        let suggested_g = v_tan * v_tan * r / m_center.max(1e-12);
        let cos_theta = cp.pole_axis().dot(d_hat);
        let emission = self.profiles[cp.profile_id].emission.sample(cos_theta);
        PackedFloat32Array::from(
            &[
                r as f32,
                v_tan as f32,
                v_rad as f32,
                suggested_g as f32,
                emission as f32,
            ][..],
        )
    }

    #[func]
    fn get_time(&self) -> f64 {
        self.time
    }

    // ── Per-particle accessors ────────────────────────────────────────

    #[func]
    fn get_particle_position(&self, id: i32) -> Vector3 {
        self.particles
            .get(id as usize)
            .map(|p| Vector3::new(p.position.x as f32, p.position.y as f32, p.position.z as f32))
            .unwrap_or(Vector3::ZERO)
    }

    #[func]
    fn get_particle_velocity(&self, id: i32) -> Vector3 {
        self.particles
            .get(id as usize)
            .map(|p| Vector3::new(p.velocity.x as f32, p.velocity.y as f32, p.velocity.z as f32))
            .unwrap_or(Vector3::ZERO)
    }

    #[func]
    fn get_particle_force(&self, id: i32) -> Vector3 {
        self.particles
            .get(id as usize)
            .map(|p| {
                Vector3::new(
                    p.force_accum.x as f32,
                    p.force_accum.y as f32,
                    p.force_accum.z as f32,
                )
            })
            .unwrap_or(Vector3::ZERO)
    }

    #[func]
    fn get_particle_pole(&self, id: i32) -> Vector3 {
        self.particles
            .get(id as usize)
            .map(|p| {
                let pole = p.pole_axis();
                Vector3::new(pole.x as f32, pole.y as f32, pole.z as f32)
            })
            .unwrap_or(Vector3::UP)
    }

    #[func]
    fn get_particle_speed(&self, id: i32) -> f64 {
        self.particles
            .get(id as usize)
            .map(|p| p.velocity.length())
            .unwrap_or(0.0)
    }

    #[func]
    fn set_particle_position(&mut self, id: i32, pos: Vector3) {
        if let Some(p) = self.particles.get_mut(id as usize) {
            p.position = DVec3::new(pos.x as f64, pos.y as f64, pos.z as f64);
        }
    }

    #[func]
    fn set_particle_velocity(&mut self, id: i32, vel: Vector3) {
        if let Some(p) = self.particles.get_mut(id as usize) {
            p.velocity = DVec3::new(vel.x as f64, vel.y as f64, vel.z as f64);
        }
    }

    #[func]
    fn set_particle_pole(&mut self, id: i32, pole: Vector3) {
        if let Some(p) = self.particles.get_mut(id as usize) {
            let target = DVec3::new(pole.x as f64, pole.y as f64, pole.z as f64);
            if target.length_squared() > 1e-9 {
                let target = target.normalize();
                if (target - DVec3::Y).length_squared() < 1e-9 {
                    p.orientation = DQuat::IDENTITY;
                } else if (target + DVec3::Y).length_squared() < 1e-9 {
                    p.orientation = DQuat::from_rotation_z(std::f64::consts::PI);
                } else {
                    p.orientation = DQuat::from_rotation_arc(DVec3::Y, target);
                }
            }
        }
    }

    // ── Bulk readback for rendering ───────────────────────────────────

    #[func]
    fn get_positions(&self) -> PackedFloat32Array {
        let mut arr = PackedFloat32Array::new();
        arr.resize(self.particles.len() * 3);
        for (i, p) in self.particles.iter().enumerate() {
            arr[i * 3] = p.position.x as f32;
            arr[i * 3 + 1] = p.position.y as f32;
            arr[i * 3 + 2] = p.position.z as f32;
        }
        arr
    }

    #[func]
    fn get_pole_axes(&self) -> PackedFloat32Array {
        let mut arr = PackedFloat32Array::new();
        arr.resize(self.particles.len() * 3);
        for (i, p) in self.particles.iter().enumerate() {
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
        arr.resize(self.particles.len());
        for (i, p) in self.particles.iter().enumerate() {
            arr[i] = p.profile_id as i32;
        }
        arr
    }

    /// 16 floats per particle: 12 (3x4 transform) + 4 (RGBA color).
    /// Row-major with interleaved origin (Godot MultiMesh format).
    #[func]
    fn build_multimesh_buffer(&self) -> PackedFloat32Array {
        let floats_per = 16usize;
        let mut buf: Vec<f32> = Vec::with_capacity(self.particles.len() * floats_per);

        for p in &self.particles {
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            let pos = p.position;

            let x = p.orientation * DVec3::X;
            let y = p.orientation * DVec3::Y;
            let z = p.orientation * DVec3::Z;

            buf.extend_from_slice(&[
                x.x as f32 * scale, x.y as f32 * scale, x.z as f32 * scale, pos.x as f32,
                y.x as f32 * scale, y.y as f32 * scale, y.z as f32 * scale, pos.y as f32,
                z.x as f32 * scale, z.y as f32 * scale, z.z as f32 * scale, pos.z as f32,
            ]);

            let (r, g, b) = profile_color(&profile.name);
            buf.extend_from_slice(&[r, g, b, 1.0]);
        }

        PackedFloat32Array::from(buf.as_slice())
    }

    /// Build a line buffer for pole axis indicators.
    /// 2 vertices per particle (center and pole tip), 6 floats each (pos xyz + color rgb).
    #[func]
    fn build_pole_indicator_buffer(&self) -> PackedFloat32Array {
        let mut buf: Vec<f32> = Vec::with_capacity(self.particles.len() * 12);
        for p in &self.particles {
            let profile = &self.profiles[p.profile_id];
            let r = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            let pos = p.position;
            let pole = p.pole_axis();
            let tip = pos + pole * (r as f64 * 1.5);

            buf.extend_from_slice(&[
                pos.x as f32, pos.y as f32, pos.z as f32,
                1.0, 0.9, 0.2,
            ]);
            buf.extend_from_slice(&[
                tip.x as f32, tip.y as f32, tip.z as f32,
                1.0, 0.4, 0.1,
            ]);
        }
        PackedFloat32Array::from(buf.as_slice())
    }

    // ── Simulation step ───────────────────────────────────────────────

    #[func]
    fn step(&mut self) {
        if !self.running || self.particles.is_empty() {
            return;
        }
        let dt = self.dt;
        let n = self.particles.len();

        // 1. Forces at current state
        self.compute_forces();

        // 2. Half-step velocities + full-step positions
        for i in 0..n {
            let p = &mut self.particles[i];
            let inv_m = 1.0 / self.profiles[p.profile_id].mass;
            let a = p.force_accum * inv_m;
            let alpha = p.torque_accum;

            p.velocity += a * (dt * 0.5);
            p.angular_velocity += alpha * (dt * 0.5);
            p.position += p.velocity * dt;

            let w = p.angular_velocity;
            let w_len = w.length();
            if w_len > 1e-12 {
                let rot = DQuat::from_axis_angle(w / w_len, w_len * dt);
                p.orientation = (rot * p.orientation).normalize();
            }
        }

        // 3. Forces at new positions
        self.compute_forces();

        // 4. Complete velocity step + damping
        for i in 0..n {
            let p = &mut self.particles[i];
            let inv_m = 1.0 / self.profiles[p.profile_id].mass;
            let a = p.force_accum * inv_m;
            let alpha = p.torque_accum;

            p.velocity += a * (dt * 0.5);
            p.angular_velocity += alpha * (dt * 0.5);
            p.velocity *= DAMPING;
            // Preserve axial spin (intrinsic), only damp tumble/precession
            let pole = p.pole_axis();
            let w_axial = pole * p.angular_velocity.dot(pole);
            let w_tumble = p.angular_velocity - w_axial;
            p.angular_velocity = w_axial + w_tumble * ANGULAR_DAMPING;
        }

        self.time += dt;
    }

    #[func]
    fn step_n(&mut self, steps: i32) {
        for _ in 0..steps.max(0) {
            self.step();
        }
    }

    // ── Diagnostics ───────────────────────────────────────────────────

    #[func]
    fn get_total_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for p in &self.particles {
            let m = self.profiles[p.profile_id].mass;
            ke += 0.5 * m * p.velocity.length_squared();
        }
        ke
    }

    #[func]
    fn get_pair_distance(&self, a: i32, b: i32) -> f64 {
        let pa = self.particles.get(a as usize);
        let pb = self.particles.get(b as usize);
        match (pa, pb) {
            (Some(a), Some(b)) => (a.position - b.position).length(),
            _ => -1.0,
        }
    }

    // ── Profile mesh generation ─────────────────────────────────────────

    /// Build a surface-of-revolution mesh from a profile's emission table.
    /// Returns a packed buffer: [vert_count, idx_count, verts...(x,y,z,nx,ny,nz,r,g,b,a), indices...]
    #[func]
    fn build_profile_mesh(
        &self,
        profile_id: i32,
        lon_segments: i32,
        lat_segments: i32,
    ) -> PackedFloat32Array {
        let pid = profile_id as usize;
        let prof = match self.profiles.get(pid) {
            Some(p) => p,
            None => return PackedFloat32Array::new(),
        };

        let lon = lon_segments.max(8) as usize;
        let lat = lat_segments.max(4) as usize;
        let rings = lat * 2 + 1;
        let verts_per_ring = lon + 1;
        let total_verts = rings * verts_per_ring;
        let quad_count = (rings - 1) * lon;
        let total_indices = quad_count * 6;

        let floats_per_vert = 10; // pos(3) + normal(3) + color(4)
        let mut buf: Vec<f32> =
            Vec::with_capacity(2 + total_verts * floats_per_vert + total_indices);

        buf.push(total_verts as f32);
        buf.push(total_indices as f32);

        let (type_r, type_g, type_b) = profile_color(&prof.name);

        for lat_idx in 0..rings {
            let theta =
                std::f64::consts::PI * lat_idx as f64 / (rings - 1).max(1) as f64;
            let cos_theta = theta.cos();
            let sin_theta = theta.sin();

            let emission = prof.emission.sample(cos_theta);
            // Lathe approach: fixed height (y = cos θ), emission modulates width only
            let r_cross = sin_theta * (ENVELOPE_MIN_R + emission * (1.0 - ENVELOPE_MIN_R));
            let y_pos = cos_theta;

            let intensity = (0.3f32 + emission as f32 * 0.7).min(1.0);
            let alpha = 0.25 + emission as f32 * 0.25;

            for lon_idx in 0..=lon {
                let phi = 2.0 * std::f64::consts::PI * lon_idx as f64 / lon as f64;
                let x = (r_cross * phi.cos()) as f32;
                let y = y_pos as f32;
                let z = (r_cross * phi.sin()) as f32;

                buf.extend_from_slice(&[x, y, z]);

                let len = (x * x + y * y + z * z).sqrt();
                if len > 1e-6 {
                    buf.extend_from_slice(&[x / len, y / len, z / len]);
                } else {
                    buf.extend_from_slice(&[0.0, if lat_idx == 0 { 1.0 } else { -1.0 }, 0.0]);
                }

                buf.extend_from_slice(&[
                    type_r * intensity,
                    type_g * intensity,
                    type_b * intensity,
                    alpha,
                ]);
            }
        }

        for lat_idx in 0..(rings - 1) {
            for lon_idx in 0..lon {
                let tl = (lat_idx * verts_per_ring + lon_idx) as f32;
                let tr = tl + 1.0;
                let bl = ((lat_idx + 1) * verts_per_ring + lon_idx) as f32;
                let br = bl + 1.0;
                buf.extend_from_slice(&[tl, bl, tr, tr, bl, br]);
            }
        }

        PackedFloat32Array::from(buf.as_slice())
    }

    #[func]
    fn count_particles_with_profile(&self, profile_id: i32) -> i32 {
        let pid = profile_id as usize;
        self.particles
            .iter()
            .filter(|p| p.profile_id == pid)
            .count() as i32
    }

    /// Build transform+color buffer for particles of one profile type only.
    #[func]
    fn build_multimesh_buffer_for_profile(&self, profile_id: i32) -> PackedFloat32Array {
        let pid = profile_id as usize;
        let mut buf: Vec<f32> = Vec::new();

        for p in &self.particles {
            if p.profile_id != pid {
                continue;
            }
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            let pos = p.position;

            let bx = p.orientation * DVec3::X;
            let by = p.orientation * DVec3::Y;
            let bz = p.orientation * DVec3::Z;

            buf.extend_from_slice(&[
                bx.x as f32 * scale, bx.y as f32 * scale, bx.z as f32 * scale, pos.x as f32,
                by.x as f32 * scale, by.y as f32 * scale, by.z as f32 * scale, pos.y as f32,
                bz.x as f32 * scale, bz.y as f32 * scale, bz.z as f32 * scale, pos.z as f32,
            ]);

            buf.extend_from_slice(&[1.0, 1.0, 1.0, 1.0]);
        }

        PackedFloat32Array::from(buf.as_slice())
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
        if !self.vfx_enabled || self.particles.is_empty() {
            self.vfx_particles.clear();
            return PackedFloat32Array::new();
        }

        let max_pool = 2048usize;

        // Age and move existing
        for vp in &mut self.vfx_particles {
            vp.age += delta;
            vp.position += vp.velocity * delta;
        }
        self.vfx_particles.retain(|vp| vp.age < lifetime);

        // Collect per-particle data to avoid borrow conflicts
        let emit_info: Vec<(DVec3, DVec3, DVec3, DVec3, f64, (f32, f32, f32), usize)> = self
            .particles
            .iter()
            .map(|p| {
                let prof = &self.profiles[p.profile_id];
                let pole = p.pole_axis();
                let (right, forward) = build_frame(pole);
                let radius = prof.radius.max(MIN_RENDER_RADIUS as f64);
                let color = profile_color(&prof.name);
                (p.position, pole, right, forward, radius, color, p.profile_id)
            })
            .collect();

        let epp = emit_per_particle.max(0) as usize;
        for (pos, pole, right, forward, radius, color, pid) in &emit_info {
            for _ in 0..epp {
                if self.vfx_particles.len() >= max_pool {
                    break;
                }
                let theta = self.profiles[*pid].emission_cdf.sample(self.rng.next_f64());
                let north = self.rng.next_f64() > 0.5;
                let phi = self.rng.next_f64() * TAU;

                let cos_t = theta.cos();
                let sin_t = theta.sin();
                let local_y = if north { cos_t } else { -cos_t };
                let dir =
                    *right * (sin_t * phi.cos()) + *pole * local_y + *forward * (sin_t * phi.sin());
                let spawn_pos = *pos + dir * *radius;

                self.vfx_particles.push(VfxParticle {
                    position: spawn_pos,
                    velocity: dir * speed,
                    age: 0.0,
                    color: *color,
                });
            }
        }

        // Build MultiMesh buffer: 12 (transform) + 4 (color) per particle
        let n = self.vfx_particles.len();
        let mut buf = Vec::with_capacity(n * 16);
        let scale = 0.03f32;

        for vp in &self.vfx_particles {
            let fade = ((1.0 - vp.age / lifetime) as f32).max(0.0);
            buf.extend_from_slice(&[
                scale, 0.0, 0.0, vp.position.x as f32,
                0.0, scale, 0.0, vp.position.y as f32,
                0.0, 0.0, scale, vp.position.z as f32,
                vp.color.0, vp.color.1, vp.color.2, fade * 0.8,
            ]);
        }

        PackedFloat32Array::from(buf.as_slice())
    }

    #[func]
    fn get_vfx_count(&self) -> i32 {
        self.vfx_particles.len() as i32
    }

    #[func]
    fn set_vfx_enabled(&mut self, enabled: bool) {
        self.vfx_enabled = enabled;
        if !enabled {
            self.vfx_particles.clear();
        }
    }

    #[func]
    fn is_vfx_enabled(&self) -> bool {
        self.vfx_enabled
    }
}

// ── Force computation ─────────────────────────────────────────────────────

impl AtomSim {
    /// Indices of (heaviest, lightest) particles by mass, or None if fewer
    /// than two distinct particles. Used to pick the orbit center vs orbiter.
    fn center_and_orbiter(&self) -> Option<(usize, usize)> {
        if self.particles.len() < 2 {
            return None;
        }
        let mut center = 0usize;
        let mut orbiter = 0usize;
        let mut m_max = f64::MIN;
        let mut m_min = f64::MAX;
        for (i, p) in self.particles.iter().enumerate() {
            let m = self.profiles[p.profile_id].mass;
            if m > m_max {
                m_max = m;
                center = i;
            }
            if m < m_min {
                m_min = m;
                orbiter = i;
            }
        }
        if center == orbiter {
            return None;
        }
        Some((center, orbiter))
    }

    fn compute_forces(&mut self) {
        let n = self.particles.len();
        for p in &mut self.particles {
            p.force_accum = DVec3::ZERO;
            p.torque_accum = DVec3::ZERO;
        }

        // Pairwise forces
        for i in 0..n {
            for j in (i + 1)..n {
                let d_vec = self.particles[j].position - self.particles[i].position;
                let r2 = d_vec.length_squared();
                let r = r2.sqrt().max(SOFTENING);
                let r2s = r * r;
                let r4s = r2s * r2s;
                let d_hat = d_vec / r;

                let pi_prof = &self.profiles[self.particles[i].profile_id];
                let pj_prof = &self.profiles[self.particles[j].profile_id];

                let pole_i = self.particles[i].pole_axis();
                let pole_j = self.particles[j].pole_axis();

                // θ_A: angle from A's pole to direction toward B
                let cos_theta_i = pole_i.dot(d_hat);
                // θ_B: angle from B's pole to direction toward A
                let cos_theta_j = pole_j.dot(-d_hat);

                let emission_i = pi_prof.emission.sample(cos_theta_i);
                let emission_j = pj_prof.emission.sample(cos_theta_j);
                let absorption_j = pj_prof.absorption.sample(cos_theta_j);
                let absorption_i = pi_prof.absorption.sample(cos_theta_i);

                // Doppler correction: approaching particles sweep through more
                // charge photons per unit time, retreating ones sweep fewer.
                // This asymmetry extracts kinetic energy on each close pass.
                let v_rel = self.particles[j].velocity - self.particles[i].velocity;
                let v_radial = v_rel.dot(d_hat);
                let doppler = (1.0 - self.drag_coupling * v_radial).clamp(0.2, 5.0);

                // Gravity: G_q * m_i * m_j / r²  (attractive)
                let f_grav = self.g_q * pi_prof.mass * pj_prof.mass / r2s;

                // Charge on j from i: C_q * m_i * E_i(θ_i) * R_j(θ_j) / r⁴ × doppler
                let f_charge_on_j =
                    self.c_q * pi_prof.mass * emission_i * absorption_j / r4s * doppler;
                // Charge on i from j: C_q * m_j * E_j(θ_j) * R_i(θ_i) / r⁴ × doppler
                let f_charge_on_i =
                    self.c_q * pj_prof.mass * emission_j * absorption_i / r4s * doppler;

                // Ambient pressure: pushes particles into charge shadows
                let shadow = (1.0 - emission_i) * (1.0 - emission_j);
                let f_ambient = self.ambient_pressure * shadow / r2s;

                // d_hat points from i to j.
                // Gravity pulls j toward i: along -d_hat
                // Charge pushes j away from i: along +d_hat
                // Ambient pulls j toward i: along -d_hat
                let net_on_j =
                    (f_charge_on_j - f_grav - f_ambient) * d_hat;
                let net_on_i =
                    (f_grav + f_ambient - f_charge_on_i) * d_hat;

                self.particles[j].force_accum += net_on_j;
                self.particles[i].force_accum += net_on_i;

                // Vortex force: spinning emission carries tangential momentum.
                // Photons leave the surface at the spin velocity, creating
                // a "sprinkler" drag that captures nearby particles into orbit.
                if self.vortex_coupling.abs() > 1e-12 {
                    let spin_i = self.particles[i].angular_velocity.dot(pole_i);
                    let sin_vec_i = pole_i.cross(d_hat);
                    let sin_theta_i = sin_vec_i.length();
                    if sin_theta_i > 1e-9 {
                        let tangent_i = sin_vec_i / sin_theta_i;
                        let f_vort = self.vortex_coupling * pi_prof.mass
                            * emission_i * absorption_j
                            * spin_i * sin_theta_i / r4s;
                        self.particles[j].force_accum += tangent_i * f_vort;
                    }

                    let spin_j = self.particles[j].angular_velocity.dot(pole_j);
                    let sin_vec_j = pole_j.cross(-d_hat);
                    let sin_theta_j = sin_vec_j.length();
                    if sin_theta_j > 1e-9 {
                        let tangent_j = sin_vec_j / sin_theta_j;
                        let f_vort = self.vortex_coupling * pj_prof.mass
                            * emission_j * absorption_i
                            * spin_j * sin_theta_j / r4s;
                        self.particles[i].force_accum += tangent_j * f_vort;
                    }
                }

                // Contact repulsion: hard-sphere boundary at sum of radii
                let ri = pi_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let rj = pj_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let r_contact = ri + rj;
                if r < r_contact {
                    let overlap = r_contact - r;
                    let f_contact = CONTACT_STIFFNESS * overlap;
                    self.particles[j].force_accum += d_hat * f_contact;
                    self.particles[i].force_accum -= d_hat * f_contact;
                }

                // Torque: aligns equator toward charge source.
                // Uses -sin(2θ) restoring form: stable with equator facing charge.
                // torque = -2 * dot(pole, d_hat) * cross(pole, d_hat) * |f| * coupling
                let torque_j = pole_j.cross(-d_hat)
                    * (-2.0 * pole_j.dot(-d_hat) * f_charge_on_j * self.torque_coupling);
                let torque_i = pole_i.cross(d_hat)
                    * (-2.0 * pole_i.dot(d_hat) * f_charge_on_i * self.torque_coupling);

                self.particles[j].torque_accum += torque_j;
                self.particles[i].torque_accum += torque_i;
            }
        }

        // Ambient environmental forces
        for p in &mut self.particles {
            let m = self.profiles[p.profile_id].mass;
            p.force_accum += self.ambient_gravity * m;
            p.force_accum += self.ambient_charge * m;
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────

fn build_frame(pole: DVec3) -> (DVec3, DVec3) {
    let right = if pole.dot(DVec3::X).abs() < 0.9 {
        pole.cross(DVec3::X).normalize()
    } else {
        pole.cross(DVec3::Z).normalize()
    };
    let forward = right.cross(pole);
    (right, forward)
}

fn default_spin_rate(name: &str) -> f64 {
    match name {
        "proton" => TAU * 3.0,
        "neutron" => TAU * 3.0,
        "electron" => TAU * 5.0,
        _ => TAU,
    }
}

fn profile_color(name: &str) -> (f32, f32, f32) {
    match name {
        "proton" => (0.92, 0.30, 0.20),
        "neutron" => (0.35, 0.50, 0.92),
        "electron" => (0.20, 0.90, 0.35),
        _ => (0.7, 0.7, 0.7),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn make_flat_csv() -> Vec<f32> {
        vec![1.0; 181]
    }

    fn make_proton_like_csv() -> Vec<f32> {
        let mut csv = vec![0.0; 181];
        // Equator (index 90) = bright, poles (0 and 180) = dim
        for i in 0..181 {
            let deg = i as f64 - 90.0;
            let frac = 1.0 - (deg.abs() / 90.0);
            csv[i] = frac as f32;
        }
        csv
    }

    #[test]
    fn emission_table_bilateral_symmetry() {
        let csv = make_proton_like_csv();
        let table = EmissionTable::from_csv_values(&csv);
        let pole_val = table.sample(1.0);
        let anti_pole = table.sample(-1.0);
        assert!(
            (pole_val - anti_pole).abs() < 1e-6,
            "Bilateral: pole={} vs anti-pole={}",
            pole_val,
            anti_pole
        );
    }

    #[test]
    fn emission_table_equator_vs_pole() {
        let csv = make_proton_like_csv();
        let table = EmissionTable::from_csv_values(&csv);
        let equator = table.sample(0.0);
        let pole = table.sample(1.0);
        assert!(
            equator > pole,
            "Equator should be brighter: eq={} pole={}",
            equator,
            pole
        );
    }

    #[test]
    fn complement_inverse() {
        let csv = make_proton_like_csv();
        let emission = EmissionTable::from_csv_values(&csv);
        let absorption = emission.complement();
        let e = emission.sample(0.5);
        let a = absorption.sample(0.5);
        assert!(
            (e + a - 1.0).abs() < 1e-6,
            "emission + absorption should be 1.0: {} + {} = {}",
            e,
            a,
            e + a
        );
    }

    #[test]
    fn flat_emission_gives_uniform_one() {
        let table = EmissionTable::from_csv_values(&make_flat_csv());
        for i in 0..=10 {
            let cos_t = (i as f64 - 5.0) / 5.0;
            let val = table.sample(cos_t);
            assert!(
                (val - 1.0).abs() < 1e-5,
                "Flat CSV should give 1.0 everywhere, got {} at cos_theta={}",
                val,
                cos_t
            );
        }
    }
}
