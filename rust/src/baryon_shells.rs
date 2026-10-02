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
    /// Doc §15.4: how many electron-tier levels (8, 7, 6, 5 in that order)
    /// are resolved as motion instead of being folded into the base ball.
    /// 0 = the §15 ball (levels 1-8, radius 1/8); 1 = levels 1-7 (1/16)
    /// carried by level 8; 2 = levels 1-6 (1/32); 3 = levels 1-5 (1/64);
    /// 4 = levels 1-4 (1/64) spun in place by the level-5 precession.
    /// Only `ball()` uses it; the shells model ignores it.
    pub depth: usize,
    /// Rotation sense of levels 8, 7, 6, 5 (Spin-mode presets: all +1).
    pub inner_signs: [f64; 4],
    /// Current rotation angle of levels 8, 7, 6, 5.
    pub inner_angles: [f64; 4],
    /// Doc §15.8: rotation sense of the photon-tier levels 4, 3, 2, 1 (used
    /// only by the unresolved-inner-levels surface twist; +1 in every
    /// Spin-mode preset; the adopted neutron flips level 4).
    pub tier1_signs: [f64; 4],
    /// Doc §15.22: a lepton (electron) body. Mathis's electron is a photon
    /// with only the level-9 axial spin: one ball of radius `radius` at
    /// `center`, spinning about the axis at edge speed c, with no orbiting
    /// levels. It takes nothing in (too big to be recycled, diatom.pdf).
    pub lepton: bool,
}

impl ShellBody {
    pub fn new(center: DVec3, orient: DQuat, radius: f64, signs: [f64; 4]) -> Self {
        ShellBody { center, orient, radius, signs, angles: [0.0; 4], depth: 0, inner_signs: [1.0; 4], inner_angles: [0.0; 4], tier1_signs: [1.0; 4], lepton: false }
    }

    /// Angular rate of electron-tier level `8 - m` at edge speed c.
    pub fn inner_rate(&self, m: usize) -> f64 {
        self.inner_signs[m] / (inner_geom(m).2 * self.radius)
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
        for m in 0..self.depth {
            self.inner_angles[m] += self.inner_rate(m) * dt;
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

/// Doc §15: the instantaneous base particle. Levels 1-8 collapse to a ball of
/// radius amplitude_8 / amplitude_12 = 1/8 (the level-9 core), carried by
/// levels 9-12.
#[derive(Clone, Copy, Debug)]
pub struct BallState {
    pub center: DVec3,
    pub radius: f64,
    /// Centre velocity: sum_{j>=10} w_j x (o_9 - o_j). Level 9 spins the
    /// ball in place, so it adds nothing here.
    pub velocity: DVec3,
    /// Ball spin: sum_{j>=9} w_j.
    pub omega: DVec3,
    /// Doc §15.8: orientation of the ball's inner frame (the composed
    /// rotation of every resolved level, then `orient`): maps directions
    /// fixed in the unresolved inner stack to world.
    pub frame: DQuat,
}

/// (local axis, local offset, orbit radius) of electron-tier level 8 - m, in
/// body units (level-12 amplitude = 1), from Spin mode's own `SpinLevel`
/// rules: level 8 = Z orbit, 7 = Y, 6 = X, 5 = axial precession about Z.
fn inner_geom(m: usize) -> (DVec3, DVec3, f64) {
    let l = crate::spin_stack::SpinLevel::new(8 - m as u8);
    let scale = level_amplitude(12);
    let orbit = l.orbit_radius() / scale;
    let offset = if l.is_orbital() { l.spin_type.orbital_start() * orbit } else { DVec3::ZERO };
    (l.rotation_axis(), offset, orbit)
}

/// Doc §15.8: summed world angular velocity (photon spin units) of the
/// unresolved inner levels at the given phases. `levels` is innermost first
/// (index 0 = level 1); level j's axis is carried by the rotations of the
/// unresolved levels above it and then by `frame`:
/// w_j = frame * R_m(th_m) ... R_{j+1}(th_{j+1}) a_j * rate_j.
pub fn inner_twist(frame: DQuat, levels: &[(DVec3, f64)], angles: &[f64]) -> DVec3 {
    let m = levels.len();
    let mut sum = DVec3::ZERO;
    for j in 0..m {
        let mut w = levels[j].0 * levels[j].1;
        for k in (j + 1)..m {
            w = DQuat::from_axis_angle(levels[k].0, angles[k]) * w;
        }
        sum += frame * w;
    }
    sum
}

/// Base-ball radius at resolution `depth` (body units): the amplitude of the
/// outermost unresolved level, 8 - depth.
pub fn ball_radius_at_depth(depth: usize) -> f64 {
    level_amplitude(8 - depth as u8) / level_amplitude(12)
}

impl ShellBody {
    /// The instantaneous base ball at `self.depth`. The resolved chain is
    /// levels (9 - depth)..=12, innermost first; each maps x -> R_j(x + s_j)
    /// about a fixed local axis (spin_stack `compose`). The ball centre is the
    /// chain applied to the origin; every chain level carries the ball, so
    /// V = sum_j w_j x (c - o_j) and Omega = sum_j w_j.
    pub fn ball(&self) -> BallState {
        if self.lepton {
            let axis = self.orient * DVec3::Z;
            let sign = if self.signs[0] < 0.0 { -1.0 } else { 1.0 };
            return BallState { center: self.center, radius: self.radius, velocity: DVec3::ZERO, omega: axis * (sign / self.radius), frame: self.orient };
        }
        // (rotation, offset, axis, rate) per chain level, innermost first.
        let mut chain: Vec<(DQuat, DVec3, DVec3, f64)> = Vec::with_capacity(8);
        for m in (0..self.depth).rev() {
            let (axis, offset, _) = inner_geom(m);
            chain.push((DQuat::from_axis_angle(axis, self.inner_angles[m]), offset, axis, self.inner_rate(m)));
        }
        for i in 0..4 {
            chain.push((self.rot(i), local_offset(i), local_axis(i), self.rate(i)));
        }
        let n = chain.len();
        let mut center_b = DVec3::ZERO;
        let mut rot_total = DQuat::IDENTITY;
        for (r, s, _, _) in &chain {
            center_b = *r * (center_b + *s);
            rot_total = *r * rot_total;
        }
        let center = self.center + self.orient * (center_b * self.radius);
        let mut velocity = DVec3::ZERO;
        let mut omega = DVec3::ZERO;
        for j in 0..n {
            // T_{>j}: levels above j, as (rotation, translation).
            let (mut r, mut t) = (DQuat::IDENTITY, DVec3::ZERO);
            for (rk, sk, _, _) in &chain[j + 1..] {
                t = *rk * (t + *sk);
                r = *rk * r;
            }
            let o_j = self.center + self.orient * (t * self.radius);
            let w_j = self.orient * (r * chain[j].2) * chain[j].3;
            velocity += w_j.cross(center - o_j);
            omega += w_j;
        }
        BallState { center, radius: ball_radius_at_depth(self.depth) * self.radius, velocity, omega, frame: self.orient * rot_total }
    }

    /// Doc §15.8: photon radius in world units (level-1 amplitude over the
    /// level-12 amplitude, times the body radius). Angular velocities times
    /// this are in the photon's own spin units (edge speed in c for a
    /// level-1 photon), which is what the twist grip compares.
    pub fn photon_radius(&self) -> f64 {
        self.radius * level_amplitude(1) / level_amplitude(12)
    }

    /// Doc §15.8: the unresolved inner levels 1..=(8 - depth), innermost
    /// first, as (local rotation axis, signed angular rate in photon spin
    /// units), from Spin mode's `SpinLevel` rules (rate = sign / orbit
    /// radius, at edge speed c).
    pub fn unresolved_levels(&self) -> Vec<(DVec3, f64)> {
        let m = 8 - self.depth;
        (1..=m as u8)
            .map(|lv| {
                let l = crate::spin_stack::SpinLevel::new(lv);
                let orbit = l.orbit_radius() / level_amplitude(12);
                let sign = if lv <= 4 { self.tier1_signs[(4 - lv) as usize] } else { self.inner_signs[(8 - lv) as usize] };
                (l.rotation_axis(), sign * (level_amplitude(1) / level_amplitude(12)) / orbit)
            })
            .collect()
    }

    /// World position of a point `p` (body units) fixed in the base ball's
    /// frame at `self.depth` (for kinematics tests).
    pub fn carry_ball(&self, p: DVec3) -> DVec3 {
        let mut x = p;
        for m in (0..self.depth).rev() {
            let (axis, offset, _) = inner_geom(m);
            x = DQuat::from_axis_angle(axis, self.inner_angles[m]) * (x + offset);
        }
        for i in 0..4 {
            x = self.rot(i) * (x + local_offset(i));
        }
        self.center + self.orient * (x * self.radius)
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
    /// The ball's centre velocity is the time derivative of its centre, and
    /// its surface point velocity V_c + Omega x (x - c) matches the composed
    /// motion of a point fixed on the ball.
    #[test]
    fn ball_kinematics_match_finite_differences() {
        let mut worst: f64 = 0.0;
        for signs in PROTON_STATES.iter().chain(NEUTRON_STATES.iter()) {
            let mut body = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 1.0, *signs);
            for step in 0..30 {
                body.angles = [0.3 * step as f64, 0.7 + 0.11 * step as f64, 0.2 * step as f64, 1.0 - 0.05 * step as f64];
                let b = body.ball();
                let h = 1e-6;
                let mut plus = body.clone();
                let mut minus = body.clone();
                for i in 0..4 {
                    plus.angles[i] += body.rate(i) * h;
                    minus.angles[i] -= body.rate(i) * h;
                }
                let fd_c = (plus.ball().center - minus.ball().center) / (2.0 * h);
                worst = worst.max((fd_c - b.velocity).length() / fd_c.length().max(1e-9));
                // A point on the ball surface, fixed in the level-9 input frame.
                let p = DVec3::new(0.08, -0.05, 0.07).normalize() * shell_amplitude(0);
                let x = body.carry(0, p);
                let fd_x = (plus.carry(0, p) - minus.carry(0, p)) / (2.0 * h);
                let v = b.velocity + b.omega.cross(x - b.center);
                worst = worst.max((fd_x - v).length() / fd_x.length().max(1e-9));
            }
        }
        println!("ball kinematics worst rel err {worst:.2e}");
        assert!(worst < 1e-6);
    }

    /// Doc §15.4: the deeper balls. Depth 0 reproduces the §15 ball exactly;
    /// every depth's centre velocity and surface velocity match finite
    /// differences of the composed motion; each depth's ball stays inside the
    /// previous depth's ball; and Spin mode's own level rules agree with this
    /// module for levels 9-12.
    #[test]
    fn deeper_ball_kinematics() {
        let mut worst: f64 = 0.0;
        for signs in PROTON_STATES.iter().chain(NEUTRON_STATES.iter()) {
            for depth in 0..=4 {
                let mut body = ShellBody::new(DVec3::new(0.2, 0.1, -0.3), DQuat::from_rotation_arc(DVec3::Z, DVec3::new(0.3, 1.0, 0.2).normalize()), 1.3, *signs);
                body.depth = depth;
                for step in 0..25 {
                    body.angles = [0.3 * step as f64, 0.7 + 0.11 * step as f64, 0.2 * step as f64, 1.0 - 0.05 * step as f64];
                    body.inner_angles = [0.9 * step as f64, 0.4 + 1.3 * step as f64, 2.0 - 0.7 * step as f64, 0.17 * step as f64];
                    let b = body.ball();
                    if depth == 0 {
                        let f = body.frames();
                        assert!((b.center - f[0].origin).length() < 1e-12 && (b.radius - f[0].radius).abs() < 1e-12);
                        let v0 = matter_velocity(&f, 0, f[0].origin);
                        assert!((b.velocity - v0).length() < 1e-12);
                    } else {
                        let mut up = body.clone();
                        up.depth = depth - 1;
                        let bp = up.ball();
                        assert!((b.center - bp.center).length() + b.radius <= bp.radius * (1.0 + 1e-9) + 1e-12, "depth {depth} ball leaves its parent");
                    }
                    let h = 1e-6;
                    let (mut plus, mut minus) = (body.clone(), body.clone());
                    plus.advance(h);
                    minus.advance(-h);
                    let fd_c = (plus.ball().center - minus.ball().center) / (2.0 * h);
                    worst = worst.max((fd_c - b.velocity).length() / fd_c.length().max(1e-9));
                    let p = DVec3::new(0.08, -0.05, 0.07).normalize() * ball_radius_at_depth(depth);
                    let x = body.carry_ball(p);
                    assert!(((x - b.center).length() - b.radius).abs() < 1e-9);
                    let fd_x = (plus.carry_ball(p) - minus.carry_ball(p)) / (2.0 * h);
                    worst = worst.max((fd_x - (b.velocity + b.omega.cross(x - b.center))).length() / fd_x.length().max(1e-9));
                }
            }
        }
        for i in 0..4 {
            let l = crate::spin_stack::SpinLevel::new(SHELL_LEVELS[i]);
            assert!((l.rotation_axis() - local_axis(i)).length() < 1e-12);
            assert!((l.orbit_radius() / level_amplitude(12) - shell_orbit_radius(i)).abs() < 1e-12);
        }
        println!("deeper ball kinematics worst rel err {worst:.2e}");
        assert!(worst < 1e-6);
    }

    /// Doc §15.4 report: per depth, the ball's centre speed (mean / max), its
    /// fastest surface point, the fraction of time it blocks the pole axis,
    /// its centre-latitude spread, and the keyframe error: how far the true
    /// centre strays from the straight line between two keyframes dt/8 apart
    /// (dt = 0.05), relative to the ball radius.
    #[test]
    #[ignore]
    fn report_deeper_ball_map() {
        for (name, signs) in [("proton +a+x+y+z", PROTON_STATES[0]), ("neutron -a-x-y+z", NEUTRON_STATES[0])] {
            for depth in 0..=4 {
                let mut body = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 1.0, signs);
                body.depth = depth;
                let (n, dt) = (400_000, 0.00731);
                let h = 0.05 / 8.0;
                let (mut vs, mut vmax, mut smax, mut blocked, mut kerr) = (0.0, 0.0f64, 0.0f64, 0u64, 0.0f64);
                let mut lat_hist = [0u64; 9];
                for _ in 0..n {
                    body.advance(dt);
                    let b = body.ball();
                    let v = b.velocity.length();
                    vs += v;
                    vmax = vmax.max(v);
                    smax = smax.max(v + b.omega.length() * b.radius);
                    let c = b.center;
                    if (c.x * c.x + c.y * c.y).sqrt() < b.radius {
                        blocked += 1;
                    }
                    let lat = (c.z / c.length().max(1e-12)).clamp(-1.0, 1.0).asin().to_degrees().abs();
                    lat_hist[((lat / 10.0) as usize).min(8)] += 1;
                    let mut mid = body.clone();
                    mid.advance(h / 2.0);
                    let mut end = body.clone();
                    end.advance(h);
                    let lin = (b.center + end.ball().center) / 2.0;
                    kerr = kerr.max((mid.ball().center - lin).length() / b.radius);
                }
                println!(
                    "{name} depth {depth}: ball r {:.4} | centre speed mean {:.2} max {:.2} c | fastest surface {:.2} c | axis blocked {:.4} | keyframe err max {:.3} r | |lat| 10-deg bins {:?}",
                    ball_radius_at_depth(depth),
                    vs / n as f64,
                    vmax,
                    smax,
                    blocked as f64 / n as f64,
                    kerr,
                    lat_hist.iter().map(|&x| format!("{:.2}", x as f64 / n as f64)).collect::<Vec<_>>()
                );
            }
        }
    }

    /// Doc §15.8: the unresolved levels carry Spin mode's rates in photon
    /// units (levels 1 and 2 spin at exactly edge c for a level-1 photon,
    /// level 4 at 1/4, level 8 at 1/32), follow the signs, and the ball
    /// frame composes the resolved rotations (checked by carrying a fixed
    /// direction through `carry_ball`).
    #[test]
    fn unresolved_levels_rates_and_frame() {
        let mut body = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 1.0, NEUTRON_STATES[0]);
        body.tier1_signs = [-1.0, 1.0, 1.0, 1.0];
        body.inner_signs = [-1.0, 1.0, 1.0, 1.0];
        let lv = body.unresolved_levels();
        assert_eq!(lv.len(), 8);
        let rates: Vec<f64> = lv.iter().map(|l| l.1).collect();
        let want = [1.0, 1.0, 0.5, -0.25, 0.125, 0.125, 0.0625, -0.03125];
        for (a, b) in rates.iter().zip(want) {
            assert!((a - b).abs() < 1e-12, "rates {rates:?}");
        }
        for depth in 0..=4 {
            body.depth = depth;
            body.angles = [0.3, 1.1, -0.4, 2.0];
            body.inner_angles = [0.7, -1.3, 0.2, 0.9];
            assert_eq!(body.unresolved_levels().len(), 8 - depth);
            let b = body.ball();
            let d = DVec3::new(0.3, -0.5, 0.8).normalize() * 1e-3;
            let moved = body.carry_ball(d) - body.carry_ball(DVec3::ZERO);
            assert!((moved - b.frame * d).length() < 1e-12, "frame at depth {depth}");
        }
        let f = DQuat::from_rotation_y(0.4);
        let one = [(DVec3::Z, 2.0)];
        assert!((inner_twist(f, &one, &[1.0]) - f * DVec3::Z * 2.0).length() < 1e-12);
    }

    /// Report: where the ball's centre goes (the "donut"). A long run with
    /// incommensurate phases, binned by the centre's polar angle and
    /// cylindrical radius. Prints the fraction of time the ball blocks the
    /// polar axis (its surface within rho < r_ball of the axis).
    #[test]
    #[ignore]
    fn report_ball_swept_map() {
        for (name, signs) in [("proton +a+x+y+z", PROTON_STATES[0]), ("neutron -a-x-y+z", NEUTRON_STATES[0])] {
            let mut body = ShellBody::new(DVec3::ZERO, DQuat::IDENTITY, 1.0, signs);
            let (n, dt) = (2_000_000, 0.0137);
            let mut lat_hist = [0u64; 18];
            let mut rho_hist = [0u64; 10];
            let mut axis_blocked = 0u64;
            let mut speed_sum = 0.0;
            for _ in 0..n {
                body.advance(dt);
                let b = body.ball();
                let c = b.center;
                let r = c.length().max(1e-12);
                let lat = (c.z / r).clamp(-1.0, 1.0).asin().to_degrees();
                lat_hist[(((lat + 90.0) / 10.0) as usize).min(17)] += 1;
                let rho = (c.x * c.x + c.y * c.y).sqrt();
                rho_hist[((rho / 0.1) as usize).min(9)] += 1;
                if rho < b.radius {
                    axis_blocked += 1;
                }
                speed_sum += b.velocity.length();
            }
            println!("=== ball swept map: {name} ({n} samples) ===");
            println!("  centre |r| is always {:.4} +- (varies); mean centre speed {:.3} c; axis blocked {:.4} of the time", body.ball().center.length(), speed_sum / n as f64, axis_blocked as f64 / n as f64);
            println!("  centre latitude (10 deg bins, -90..90): {:?}", lat_hist.iter().map(|&h| format!("{:.3}", h as f64 / n as f64)).collect::<Vec<_>>());
            println!("  centre cylindrical radius (0.1 bins): {:?}", rho_hist.iter().map(|&h| format!("{:.3}", h as f64 / n as f64)).collect::<Vec<_>>());
        }
    }

}
