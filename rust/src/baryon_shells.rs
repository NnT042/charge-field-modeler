//! Mechanical baryon, layered-shell model (docs/GPU_DSMC_DESIGN.md §14).
//!
//! The baryon-tier levels 9-12 (axial precession, x, y, z) become four nested
//! permeable shells that photons pass through and collide with. The
//! kinematics are exactly `spin_stack.rs`'s, normalised so level 12's
//! amplitude (512) is the body radius 1:
//! - Level j maps inner-frame points by T_j(p) = R_j(theta_j) (p + s_j).
//!   R_j is a rotation about its local axis. s_j is the orbital start offset
//!   times the orbit radius (zero for the level-9 precession).
//! - Outer levels carry inner ones, so a point fixed in level k's inner frame
//!   sits at T_12 o ... o T_k (p) in body coordinates.
//! - Level k's shell region is its Spin-mode ghost sphere: radius amplitude_k,
//!   centred at o_k = T_{>k}(0). Level 9's region is the core of radius
//!   amplitude_9, centred at its rotation centre.
//! - Matter carried by levels >= k moves at V_k(x) = sum_{j>=k} w_j x (x - o_j),
//!   where w_j = R_{>j} a_j * dtheta_j/dt. That is exact for nested rigid
//!   rotations; `kinematics_match_finite_differences` checks it against the
//!   composed motion.
//!
//! Body coordinates have the pole along +Z. `orient` maps body to world.

use glam::{DQuat, DVec3};

use crate::types::{level_amplitude, SpinType};

/// Baryon-tier levels, innermost first.
pub const SHELL_LEVELS: [u8; 4] = [9, 10, 11, 12];

fn level_spin_type(i: usize) -> SpinType {
    [SpinType::Axial, SpinType::X, SpinType::Y, SpinType::Z][i]
}

/// Amplitude normalised to level 12 = 1: 1/8, 1/4, 1/2, 1.
pub fn shell_amplitude(i: usize) -> f64 {
    level_amplitude(SHELL_LEVELS[i]) / level_amplitude(12)
}

/// spin_stack's `orbit_radius`: amplitude for the axial precession, amplitude/2
/// for orbital levels.
pub fn shell_orbit_radius(i: usize) -> f64 {
    if i == 0 {
        shell_amplitude(0)
    } else {
        shell_amplitude(i) / 2.0
    }
}

/// Local rotation axis (spin_stack `rotation_axis`: tier-2+ axial precesses
/// about Z).
fn local_axis(i: usize) -> DVec3 {
    if i == 0 {
        DVec3::Z
    } else {
        level_spin_type(i).rotation_axis()
    }
}

/// Local offset added before rotating (spin_stack `compose`).
fn local_offset(i: usize) -> DVec3 {
    if i == 0 {
        DVec3::ZERO
    } else {
        level_spin_type(i).orbital_start() * shell_orbit_radius(i)
    }
}

/// Sign patterns (a, x, y, z) for levels 9-12 (PHYSICS_REFERENCE §3).
pub const PROTON_STATES: [[f64; 4]; 4] = [[1.0, 1.0, 1.0, 1.0], [-1.0, -1.0, -1.0, -1.0], [-1.0, 1.0, 1.0, 1.0], [1.0, -1.0, -1.0, -1.0]];
pub const NEUTRON_STATES: [[f64; 4]; 4] = [[-1.0, -1.0, -1.0, 1.0], [1.0, 1.0, 1.0, -1.0], [-1.0, 1.0, 1.0, -1.0], [1.0, -1.0, -1.0, 1.0]];

/// One level's world-space state at an instant.
#[derive(Clone, Copy, Debug)]
pub struct LevelFrame {
    /// Rotation centre = shell centre (world).
    pub origin: DVec3,
    /// Shell radius (world units; body radius = `ShellBody::radius`).
    pub radius: f64,
    /// World angular velocity of this level's own rotation.
    pub omega: DVec3,
}

#[derive(Clone, Debug)]
pub struct ShellBody {
    pub center: DVec3,
    pub orient: DQuat,
    /// Body radius in world units (the level-12 amplitude).
    pub radius: f64,
    /// Rotation sense per level (a, x, y, z), each at edge speed c.
    pub signs: [f64; 4],
    /// Current rotation angle per level.
    pub angles: [f64; 4],
}

impl ShellBody {
    pub fn new(center: DVec3, orient: DQuat, radius: f64, signs: [f64; 4]) -> Self {
        ShellBody { center, orient, radius, signs, angles: [0.0; 4] }
    }

    /// Angular rate of level i at edge speed c (spin_stack `advance` with
    /// v = +-1): sign / orbit radius, in world units.
    pub fn rate(&self, i: usize) -> f64 {
        self.signs[i] / (shell_orbit_radius(i) * self.radius)
    }

    pub fn advance(&mut self, dt: f64) {
        for i in 0..4 {
            self.angles[i] += self.rate(i) * dt;
        }
    }

    fn rot(&self, i: usize) -> DQuat {
        DQuat::from_axis_angle(local_axis(i), self.angles[i])
    }

    /// T_{>i} as (rotation, translation) in body coordinates (unit radius):
    /// maps level-i output-frame points to the body frame.
    fn outer_transform(&self, i: usize) -> (DQuat, DVec3) {
        let mut r = DQuat::IDENTITY;
        let mut t = DVec3::ZERO;
        // Compose inner-to-outer: x -> R_j (x + s_j), applied for j = i+1..3.
        for j in (i + 1)..4 {
            let rj = self.rot(j);
            // New map: x -> R_j((r x + t) + s_j) = (R_j r) x + R_j (t + s_j)
            t = rj * (t + local_offset(j));
            r = rj * r;
        }
        (r, t)
    }

    /// World-space frames of all four levels at the current angles.
    pub fn frames(&self) -> [LevelFrame; 4] {
        let mut out = [LevelFrame { origin: DVec3::ZERO, radius: 0.0, omega: DVec3::ZERO }; 4];
        for i in 0..4 {
            let (r, t) = self.outer_transform(i);
            let axis_body = r * local_axis(i);
            out[i] = LevelFrame {
                origin: self.center + self.orient * (t * self.radius),
                radius: shell_amplitude(i) * self.radius,
                omega: self.orient * axis_body * self.rate(i),
            };
        }
        out
    }

    /// World position of a point `p` fixed in level i's input frame (body
    /// units), carried by levels i..=3.
    pub fn carry(&self, i: usize, p: DVec3) -> DVec3 {
        let ri = self.rot(i);
        let q = ri * (p + local_offset(i));
        let (r, t) = self.outer_transform(i);
        self.center + self.orient * ((r * q + t) * self.radius)
    }
}

/// Innermost shell containing `x`, or None outside the body.
pub fn region_of(frames: &[LevelFrame; 4], x: DVec3) -> Option<usize> {
    (0..4).find(|&k| (x - frames[k].origin).length_squared() <= frames[k].radius * frames[k].radius)
}

/// Velocity of body matter in shell k at x: sum_{j>=k} w_j x (x - o_j).
pub fn matter_velocity(frames: &[LevelFrame; 4], k: usize, x: DVec3) -> DVec3 {
    (k..4).map(|j| frames[j].omega.cross(x - frames[j].origin)).sum()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// V_k(x) must equal the time derivative of the composed Spin-mode motion
    /// of a point fixed in level k's frame, for every level, sign pattern and
    /// a spread of angles. The negative control drops the outer levels'
    /// contribution and must fail.
    #[test]
    fn kinematics_match_finite_differences() {
        let orient = DQuat::from_rotation_arc(DVec3::Z, DVec3::new(1.0, 0.3, -0.2).normalize());
        let mut worst: f64 = 0.0;
        let mut worst_wrong: f64 = 0.0;
        for signs in PROTON_STATES.iter().chain(NEUTRON_STATES.iter()) {
            let mut body = ShellBody::new(DVec3::new(0.3, -1.0, 2.0), orient, 1.7, *signs);
            for step in 0..40 {
                body.angles = [0.37 * step as f64, 1.1 + 0.21 * step as f64, 0.5 + 0.13 * step as f64, 2.0 - 0.07 * step as f64];
                let frames = body.frames();
                for k in 0..4 {
                    let p = DVec3::new(0.03 * (k as f64 + 1.0), -0.02, 0.05);
                    let x = body.carry(k, p);
                    let h = 1e-6;
                    let mut plus = body.clone();
                    let mut minus = body.clone();
                    for i in 0..4 {
                        plus.angles[i] += body.rate(i) * h;
                        minus.angles[i] -= body.rate(i) * h;
                    }
                    let fd = (plus.carry(k, p) - minus.carry(k, p)) / (2.0 * h);
                    let v = matter_velocity(&frames, k, x);
                    worst = worst.max((fd - v).length() / fd.length().max(1e-9));
                    let wrong = frames[k].omega.cross(x - frames[k].origin);
                    if k < 3 {
                        worst_wrong = worst_wrong.max((fd - wrong).length() / fd.length().max(1e-9));
                    }
                }
            }
        }
        println!("shell kinematics: worst rel err {worst:.2e}, own-level-only control {worst_wrong:.2e}");
        assert!(worst < 1e-6, "matter velocity disagrees with the composed motion: {worst}");
        assert!(worst_wrong > 0.1, "negative control should fail");
    }

    /// Geometry: shells nest off-centre, each inner ghost touching the next
    /// from inside (spin_stack: orbit radius = amplitude/2), and level 12's
    /// shell is the whole body.
    #[test]
    fn shells_nest_like_ghost_spheres() {
        let mut body = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 1.0, PROTON_STATES[0]);
        for s in 0..50 {
            body.advance(0.037 * s as f64);
            let f = body.frames();
            assert!((f[3].origin).length() < 1e-12 && (f[3].radius - 1.0).abs() < 1e-12);
            for k in 1..3 {
                // Inner ghost k touches ghost k+1 from inside.
                let d = (f[k].origin - f[k + 1].origin).length();
                assert!((d + f[k].radius - f[k + 1].radius).abs() < 1e-12, "level {k} ghost not internally tangent");
            }
            // The core (level 9) sits at level 10's orbit point, inside ghost 10.
            let d = (f[0].origin - f[1].origin).length();
            assert!((d - shell_orbit_radius(1)).abs() < 1e-12);
        }
    }

    /// Every level's own rotation runs its shell's edge at c.
    #[test]
    fn each_level_edge_at_c() {
        let body = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 2.5, NEUTRON_STATES[2]);
        let f = body.frames();
        for i in 0..4 {
            let edge = f[i].omega.length() * shell_orbit_radius(i) * body.radius;
            assert!((edge - 1.0).abs() < 1e-12, "level {i} edge {edge}");
        }
    }
}
