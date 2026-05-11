use godot::classes::INode;
use godot::prelude::*;

const WORLD_SCALE: f32 = 0.5;

#[repr(C)]
#[derive(Clone, Copy)]
pub struct FieldPhoton {
    pub position: [f32; 3],
    pub energy: f32,
    pub velocity: [f32; 3],
    pub flags: u32,
}

const FLAG_CHIRALITY: u32 = 1;
const FLAG_ACTIVE: u32 = 1 << 1;
#[allow(dead_code)]
const FLAG_BOUNCED: u32 = 1 << 2;

#[derive(Clone, Copy)]
struct RicochetPhoton {
    position: [f32; 3],
    velocity: [f32; 3],
    age: f32,
    exited_ghost: bool,
}

#[derive(GodotClass)]
#[class(base=Node)]
pub struct FieldSim {
    base: Base<Node>,
    photons: Vec<FieldPhoton>,
    sim_radius: f32,
    photon_ratio: f32,
    field_direction: [f32; 3],
    direction_strength: f32,
    direction_scatter: f32,
    direction_mask: u32,
    effective_radius: f32,
    running: bool,
    last_collision_count: i32,
    ricochets: Vec<RicochetPhoton>,
    exit_holes: Vec<[f32; 3]>,
    targeted_rate: f32,
    targeted_accum: f32,
    targeted_impulse: [f32; 3],
    targeted_hits: i32,
    stream_jitter: f32,
    rng_state: u32,
}

#[godot_api]
impl INode for FieldSim {
    fn init(base: Base<Node>) -> Self {
        Self {
            base,
            photons: Vec::new(),
            sim_radius: 10.0,
            photon_ratio: 0.67,
            field_direction: [0.0, -1.0, 0.0],
            direction_strength: 0.0,
            direction_scatter: 0.3,
            direction_mask: 0,
            effective_radius: 1.0,
            running: false,
            last_collision_count: 0,
            ricochets: Vec::new(),
            exit_holes: Vec::new(),
            targeted_rate: 0.0,
            targeted_accum: 0.0,
            targeted_impulse: [0.0; 3],
            targeted_hits: 0,
            stream_jitter: 0.0,
            rng_state: 0xBADC0FFE,
        }
    }
}

fn xorshift32(state: &mut u32) -> u32 {
    let mut x = *state;
    x ^= x << 13;
    x ^= x >> 17;
    x ^= x << 5;
    *state = x;
    x
}

fn rand_f32(state: &mut u32) -> f32 {
    (xorshift32(state) & 0x7F_FFFF) as f32 / 0x7F_FFFF as f32
}

fn rand_in_sphere(state: &mut u32, radius: f32) -> [f32; 3] {
    let r = radius * rand_f32(state).cbrt();
    let theta = std::f32::consts::TAU * rand_f32(state);
    let cos_phi = 2.0 * rand_f32(state) - 1.0;
    let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
    [
        r * sin_phi * theta.cos(),
        r * sin_phi * theta.sin(),
        r * cos_phi,
    ]
}

fn rand_direction(state: &mut u32) -> [f32; 3] {
    let theta = std::f32::consts::TAU * rand_f32(state);
    let cos_phi = 2.0 * rand_f32(state) - 1.0;
    let sin_phi = (1.0 - cos_phi * cos_phi).sqrt();
    [sin_phi * theta.cos(), sin_phi * theta.sin(), cos_phi]
}


#[godot_api]
impl FieldSim {
    #[func]
    fn init_field(&mut self, count: i32) {
        let n = count.max(0) as usize;
        self.photons.clear();
        self.photons.reserve(n);

        let mut rng: u32 = 0xDEAD_BEEF;
        let r = self.sim_radius;
        let ratio = self.photon_ratio;
        let force_respawn = self.direction_mask != 0;

        for _ in 0..n {
            let position = if force_respawn {
                [r * 2.0, 0.0, 0.0]
            } else {
                rand_in_sphere(&mut rng, r)
            };
            let velocity = rand_direction(&mut rng);
            let chirality = if rand_f32(&mut rng) > ratio {
                FLAG_CHIRALITY
            } else {
                0
            };
            self.photons.push(FieldPhoton {
                position,
                energy: 1.0,
                velocity,
                flags: FLAG_ACTIVE | chirality,
            });
        }
    }

    #[func]
    fn get_init_buffer(&self) -> PackedByteArray {
        let byte_len = self.photons.len() * std::mem::size_of::<FieldPhoton>();
        let src = unsafe {
            std::slice::from_raw_parts(self.photons.as_ptr() as *const u8, byte_len)
        };
        PackedByteArray::from(src)
    }

    #[func]
    fn parse_readback(&mut self, data: PackedByteArray) {
        let photon_size = std::mem::size_of::<FieldPhoton>();
        let count = data.len() / photon_size;
        self.photons.resize(
            count,
            FieldPhoton {
                position: [0.0; 3],
                energy: 1.0,
                velocity: [0.0, 0.0, 1.0],
                flags: 0,
            },
        );

        let src = data.as_slice();
        let dst = unsafe {
            std::slice::from_raw_parts_mut(
                self.photons.as_mut_ptr() as *mut u8,
                count * photon_size,
            )
        };
        dst.copy_from_slice(&src[..count * photon_size]);

        // Exit hole detection moved to detect_exit_holes_from_impulses
    }

    #[func]
    fn build_multimesh_buffer(&self) -> PackedFloat32Array {
        let floats_per: usize = 16;
        let total = self.photons.len() + self.ricochets.len();
        let mut vec: Vec<f32> = Vec::with_capacity(total * floats_per);

        for p in self.photons.iter() {
            let px = p.position[0] * WORLD_SCALE;
            let py = p.position[1] * WORLD_SCALE;
            let pz = p.position[2] * WORLD_SCALE;
            vec.extend_from_slice(&[
                1.0, 0.0, 0.0, px,
                0.0, 1.0, 0.0, py,
                0.0, 0.0, 1.0, pz,
            ]);
            if p.flags & FLAG_CHIRALITY == 0 {
                vec.extend_from_slice(&[0.3, 0.5, 1.0, 0.8]);
            } else {
                vec.extend_from_slice(&[1.0, 0.3, 0.3, 0.8]);
            }
        }

        for r in self.ricochets.iter() {
            let px = r.position[0] * WORLD_SCALE;
            let py = r.position[1] * WORLD_SCALE;
            let pz = r.position[2] * WORLD_SCALE;
            vec.extend_from_slice(&[
                1.0, 0.0, 0.0, px,
                0.0, 1.0, 0.0, py,
                0.0, 0.0, 1.0, pz,
                0.2, 1.0, 0.3, 1.0,
            ]);
        }

        PackedFloat32Array::from(vec.as_slice())
    }

    #[func]
    fn sum_partial_impulses(&mut self, data: PackedByteArray) -> Vector3 {
        let float_size = std::mem::size_of::<f32>();
        let vec4_size = float_size * 4;
        let count = data.len() / vec4_size;
        let src = data.as_slice();

        let mut sum = [0.0f32; 3];
        let mut hit_count = 0.0f32;
        for i in 0..count {
            let offset = i * vec4_size;
            let x = f32::from_le_bytes(src[offset..offset + 4].try_into().unwrap());
            let y = f32::from_le_bytes(src[offset + 4..offset + 8].try_into().unwrap());
            let z = f32::from_le_bytes(src[offset + 8..offset + 12].try_into().unwrap());
            let w = f32::from_le_bytes(src[offset + 12..offset + 16].try_into().unwrap());
            sum[0] += x;
            sum[1] += y;
            sum[2] += z;
            hit_count += w;
        }

        self.last_collision_count = hit_count as i32;
        Vector3::new(sum[0], sum[1], sum[2])
    }

    #[func]
    fn get_last_collision_count(&self) -> i32 {
        self.last_collision_count
    }

    #[func]
    fn get_photon_count(&self) -> i32 {
        self.photons.len() as i32
    }

    #[func]
    fn set_sim_radius(&mut self, radius: f64) {
        self.sim_radius = radius.max(1.0) as f32;
    }

    #[func]
    fn get_sim_radius(&self) -> f64 {
        self.sim_radius as f64
    }

    #[func]
    fn adapt_sim_radius(&mut self, effective_radius: f64) {
        self.sim_radius = (effective_radius * 3.0).max(10.0) as f32;
        self.effective_radius = effective_radius.max(1.0) as f32;
    }

    #[func]
    fn set_photon_ratio(&mut self, ratio: f64) {
        self.photon_ratio = ratio.clamp(0.0, 1.0) as f32;
    }

    #[func]
    fn get_photon_ratio(&self) -> f64 {
        self.photon_ratio as f64
    }

    #[func]
    fn set_field_direction(&mut self, dir: Vector3) {
        let len = (dir.x * dir.x + dir.y * dir.y + dir.z * dir.z).sqrt();
        if len > 1e-6 {
            self.field_direction = [dir.x / len, dir.y / len, dir.z / len];
        }
    }

    #[func]
    fn set_direction_strength(&mut self, strength: f64) {
        self.direction_strength = strength.clamp(0.0, 1.0) as f32;
    }

    #[func]
    fn get_direction_strength(&self) -> f64 {
        self.direction_strength as f64
    }

    #[func]
    fn set_direction_scatter(&mut self, scatter: f64) {
        self.direction_scatter = scatter.clamp(0.0, 1.0) as f32;
    }

    #[func]
    fn get_direction_scatter(&self) -> f64 {
        self.direction_scatter as f64
    }

    #[func]
    fn set_direction_mask(&mut self, mask: i32) {
        self.direction_mask = mask as u32 & 0x3F;
    }

    #[func]
    fn get_direction_mask(&self) -> i32 {
        self.direction_mask as i32
    }

    #[func]
    fn set_stream_jitter(&mut self, v: f64) {
        self.stream_jitter = v.clamp(0.0, 1.0) as f32;
    }

    #[func]
    fn clear_field(&mut self) {
        self.photons.clear();
        self.ricochets.clear();
        self.exit_holes.clear();
        self.targeted_accum = 0.0;
        self.targeted_impulse = [0.0; 3];
        self.targeted_hits = 0;
    }

    #[func]
    fn set_running(&mut self, running: bool) {
        self.running = running;
    }

    #[func]
    fn is_running(&self) -> bool {
        self.running
    }

    #[func]
    fn set_targeted_rate(&mut self, rate: f64) {
        self.targeted_rate = rate.max(0.0) as f32;
    }

    #[func]
    fn get_targeted_rate(&self) -> f64 {
        self.targeted_rate as f64
    }

    #[func]
    fn tick_targeted(&mut self, center: Vector3, radius: f32, dt: f32, _time_scale: f32, eff_radius: f32) -> Vector3 {
        let cull_radius = eff_radius * 1.5;
        let cull_r2 = cull_radius * cull_radius;
        let max_age = cull_radius * 2.0;
        let ricochet_speed = cull_radius * 0.5;
        let step = dt * ricochet_speed;
        let eff_r2 = eff_radius * eff_radius;

        let cx = center.x;
        let cy = center.y;
        let cz = center.z;

        let mut new_exits: Vec<[f32; 3]> = Vec::new();

        self.ricochets.retain_mut(|r| {
            r.position[0] += r.velocity[0] * step;
            r.position[1] += r.velocity[1] * step;
            r.position[2] += r.velocity[2] * step;
            r.age += step;
            if r.age > max_age {
                return false;
            }

            // Exit hole: check if ricochet crossed the outermost ghost sphere (at origin)
            if !r.exited_ghost {
                let d2 = r.position[0] * r.position[0]
                    + r.position[1] * r.position[1]
                    + r.position[2] * r.position[2];
                if d2 > eff_r2 {
                    r.exited_ghost = true;
                    let inv = 1.0 / d2.sqrt();
                    new_exits.push([
                        r.position[0] * inv,
                        r.position[1] * inv,
                        r.position[2] * inv,
                    ]);
                }
            }

            // Cull relative to focus particle center
            let dx = r.position[0] - cx;
            let dy = r.position[1] - cy;
            let dz = r.position[2] - cz;
            dx * dx + dy * dy + dz * dz < cull_r2
        });

        self.exit_holes.extend_from_slice(&new_exits);

        // Spawn new collisions based on rate
        self.targeted_accum += self.targeted_rate * dt;
        let to_spawn = self.targeted_accum as i32;
        self.targeted_accum -= to_spawn as f32;

        let cx = center.x;
        let cy = center.y;
        let cz = center.z;
        let mut frame_impulse = [0.0f32; 3];

        for _ in 0..to_spawn {
            // Random point on sphere surface
            let normal = rand_direction(&mut self.rng_state);
            let pos = [
                cx + normal[0] * (radius + 0.01),
                cy + normal[1] * (radius + 0.01),
                cz + normal[2] * (radius + 0.01),
            ];

            // Random inbound direction (hemisphere facing inward)
            let mut inbound = rand_direction(&mut self.rng_state);
            let dot = inbound[0] * normal[0] + inbound[1] * normal[1] + inbound[2] * normal[2];
            if dot > 0.0 {
                inbound[0] = -inbound[0];
                inbound[1] = -inbound[1];
                inbound[2] = -inbound[2];
            }

            // Specular reflection: v_out = v_in - 2*dot(v_in, n)*n
            let vn = inbound[0] * normal[0] + inbound[1] * normal[1] + inbound[2] * normal[2];
            let velocity = [
                inbound[0] - 2.0 * vn * normal[0],
                inbound[1] - 2.0 * vn * normal[1],
                inbound[2] - 2.0 * vn * normal[2],
            ];

            // Impulse = -2 * vn * normal (energy = 1)
            frame_impulse[0] += -2.0 * vn * normal[0];
            frame_impulse[1] += -2.0 * vn * normal[1];
            frame_impulse[2] += -2.0 * vn * normal[2];

            self.ricochets.push(RicochetPhoton {
                position: pos,
                velocity,
                age: 0.0,
                exited_ghost: false,
            });
            self.targeted_hits += 1;
        }

        self.targeted_impulse[0] += frame_impulse[0];
        self.targeted_impulse[1] += frame_impulse[1];
        self.targeted_impulse[2] += frame_impulse[2];

        Vector3::new(frame_impulse[0], frame_impulse[1], frame_impulse[2])
    }

    #[func]
    fn get_targeted_hits(&self) -> i32 {
        self.targeted_hits
    }

    #[func]
    fn get_targeted_impulse(&self) -> Vector3 {
        Vector3::new(
            self.targeted_impulse[0],
            self.targeted_impulse[1],
            self.targeted_impulse[2],
        )
    }

    #[func]
    fn get_ricochet_count(&self) -> i32 {
        self.ricochets.len() as i32
    }

    #[func]
    fn get_total_display_count(&self) -> i32 {
        (self.photons.len() + self.ricochets.len()) as i32
    }

    #[func]
    fn get_exit_hole_count(&self) -> i32 {
        self.exit_holes.len() as i32
    }

    #[func]
    fn build_exit_holes_buffer(&self, eff_radius: f32) -> PackedFloat32Array {
        let floats_per = 16usize;
        let mut vec: Vec<f32> = Vec::with_capacity(self.exit_holes.len() * floats_per);
        let r = (eff_radius + 0.02) * WORLD_SCALE;
        let scale = WORLD_SCALE;

        for dir in &self.exit_holes {
            let px = dir[0] * r;
            let py = dir[1] * r;
            let pz = dir[2] * r;
            vec.extend_from_slice(&[
                scale, 0.0, 0.0, px,
                0.0, scale, 0.0, py,
                0.0, 0.0, scale, pz,
                1.0, 0.5, 0.1, 1.0,
            ]);
        }

        PackedFloat32Array::from(vec.as_slice())
    }

    #[func]
    fn detect_exit_holes_from_impulses(&mut self, impulse_data: PackedByteArray) {
        let src = impulse_data.as_slice();
        let vec4_size = 16usize;
        for i in 0..self.photons.len() {
            let offset = i * vec4_size;
            if offset + vec4_size > src.len() {
                break;
            }
            let w = f32::from_le_bytes(src[offset + 12..offset + 16].try_into().unwrap());
            if w > 0.5 {
                let p = &self.photons[i];
                let d2 = p.position[0] * p.position[0]
                    + p.position[1] * p.position[1]
                    + p.position[2] * p.position[2];
                if d2 > 0.01 {
                    let inv = 1.0 / d2.sqrt();
                    self.exit_holes.push([
                        p.position[0] * inv,
                        p.position[1] * inv,
                        p.position[2] * inv,
                    ]);
                }
            }
        }
    }

    #[func]
    fn clear_exit_holes(&mut self) {
        self.exit_holes.clear();
    }

    #[func]
    fn reset_targeted(&mut self) {
        self.ricochets.clear();
        self.exit_holes.clear();
        self.targeted_accum = 0.0;
        self.targeted_impulse = [0.0; 3];
        self.targeted_hits = 0;
    }

    #[func]
    fn pack_sim_params(
        &self,
        center: Vector3,
        dt: f32,
        time_scale: f32,
        frame: i32,
        segment_count: i32,
    ) -> PackedFloat32Array {
        let mut buf = PackedFloat32Array::new();
        buf.resize(20);

        buf[0] = center.x;
        buf[1] = center.y;
        buf[2] = center.z;
        buf[3] = self.sim_radius;

        buf[4] = dt;
        buf[5] = time_scale;
        buf[6] = self.photon_ratio;
        buf[7] = self.direction_strength;

        buf[8] = self.field_direction[0];
        buf[9] = self.field_direction[1];
        buf[10] = self.field_direction[2];
        buf[11] = self.direction_scatter;

        buf[12] = f32::from_bits(self.photons.len() as u32);
        buf[13] = f32::from_bits(frame as u32);
        buf[14] = f32::from_bits(segment_count as u32);
        buf[15] = f32::from_bits(self.direction_mask);

        buf[16] = self.effective_radius;
        buf[17] = self.stream_jitter;
        buf[18] = 0.0;
        buf[19] = 0.0;

        buf
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn photon_struct_is_32_bytes() {
        assert_eq!(std::mem::size_of::<FieldPhoton>(), 32);
    }

    #[test]
    fn photon_struct_is_aligned() {
        assert_eq!(std::mem::align_of::<FieldPhoton>(), 4);
    }

    #[test]
    fn xorshift_never_zero() {
        let mut state: u32 = 1;
        for _ in 0..1000 {
            let v = xorshift32(&mut state);
            assert_ne!(v, 0);
        }
    }

    #[test]
    fn rand_f32_in_range() {
        let mut state: u32 = 0xCAFE;
        for _ in 0..1000 {
            let v = rand_f32(&mut state);
            assert!((0.0..=1.0).contains(&v));
        }
    }

    #[test]
    fn rand_in_sphere_within_radius() {
        let mut state: u32 = 0xBEEF;
        let radius = 5.0f32;
        for _ in 0..1000 {
            let p = rand_in_sphere(&mut state, radius);
            let dist = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
            assert!(dist <= radius + 1e-5);
        }
    }

    #[test]
    fn rand_direction_is_unit() {
        let mut state: u32 = 0xFACE;
        for _ in 0..1000 {
            let d = rand_direction(&mut state);
            let len = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            assert!((len - 1.0).abs() < 1e-5);
        }
    }

}
