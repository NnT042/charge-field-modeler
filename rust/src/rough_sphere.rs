//! Rough-sphere photon contact mechanics: the CPU f64 reference
//! (docs/GPU_DSMC_DESIGN.md §13).
//!
//! This is the textbook perfectly-rough-sphere collision (Bryan 1894; Pidduck
//! 1922; Chapman & Cowling ch. 11), generalised to unequal masses and radii
//! with normal restitution `e` and tangential restitution `beta`:
//! - `beta = -1`: frictionless (smooth) spheres. Spins never change.
//! - `beta = +1`: perfectly rough spheres, "gears that don't slip". The
//!   contact slip reverses.
//!
//! Energy is conserved only when both |e| and |beta| equal 1. Linear momentum
//! and total angular momentum (orbital + spin about any fixed point) are
//! conserved for every e and beta.
//!
//! The GPU kernel (`dsmc_collide.wgsl`, exchange rule `RoughSphere`) mirrors
//! [`collide`] line for line. The gates compare the two directly
//! (`rs_gpu_matches_cpu_reference`), so keep them in sync.
//!
//! Units: speeds in c. Photon radius R_k = 2^k and mass m_k = 2^k at stacked
//! level k (bottom.pdf). With k = 0 at R = 1, the edge speed is |w|*R, in c.

use glam::DVec3;

/// Moment-of-inertia factor for a uniform solid sphere, I = kappa*m*R^2.
pub const KAPPA_SOLID: f64 = 0.4;

#[derive(Clone, Copy, Debug)]
pub struct Ball {
    pub v: DVec3,
    pub w: DVec3,
    pub m: f64,
    pub r: f64,
}

impl Ball {
    /// A photon at stacked level `k`, velocity `v` and spin vector `w`.
    pub fn level(k: u32, v: DVec3, w: DVec3) -> Self {
        let s = (2.0f64).powi(k as i32);
        Ball { v, w, m: s, r: s }
    }

    pub fn inertia(&self, kappa: f64) -> f64 {
        kappa * self.m * self.r * self.r
    }

    pub fn e_lin(&self) -> f64 {
        0.5 * self.m * self.v.length_squared()
    }

    pub fn e_rot(&self, kappa: f64) -> f64 {
        0.5 * self.inertia(kappa) * self.w.length_squared()
    }

    /// Surface (edge) speed of the spin: |w| * R, in c.
    pub fn edge_speed(&self) -> f64 {
        self.w.length() * self.r
    }
}

/// Relative velocity of the two contact points: surface of `a` minus surface
/// of `b`. `n` is the unit normal from a's centre to b's centre.
pub fn contact_slip(a: &Ball, b: &Ball, n: DVec3) -> DVec3 {
    (a.v - b.v) + a.r * a.w.cross(n) + b.r * b.w.cross(n)
}

/// Applies one collision in place and returns the impulse J on `a` (`b`
/// receives -J). If the pair is separating along `n` (g.n <= 0), nothing
/// changes and zero is returned.
pub fn collide(a: &mut Ball, b: &mut Ball, n: DVec3, e: f64, beta: f64, kappa: f64) -> DVec3 {
    let g = contact_slip(a, b, n);
    let gn = g.dot(n);
    if gn <= 0.0 {
        return DVec3::ZERO;
    }
    let g_t = g - gn * n;
    let mu = a.m * b.m / (a.m + b.m);
    let j_n = -mu * (1.0 + e) * gn * n;
    // A tangential impulse J_t changes the slip by
    // J_t * (1/ma + 1/mb + Ra^2/Ia + Rb^2/Ib).
    let k_t = (1.0 / a.m) * (1.0 + 1.0 / kappa) + (1.0 / b.m) * (1.0 + 1.0 / kappa);
    let j_t = -(1.0 + beta) * g_t / k_t;
    let j = j_n + j_t;
    a.v += j / a.m;
    b.v -= j / b.m;
    // Torques about each centre: contact point at +Ra*n on a, -Rb*n on b,
    // with force J on a and -J on b, so both torques are +R*(n x J).
    let nj = n.cross(j);
    a.w += a.r * nj / a.inertia(kappa);
    b.w += b.r * nj / b.inertia(kappa);
    j
}

/// Full contact law: sliding grip (`beta`, the classical rough sphere) plus
/// torsional ("drilling") grip about the contact normal (`beta_twist`, the
/// same -1..+1 convention applied to the relative spin about n). Torsional
/// grip is a pure couple: it exchanges spin between the two photons and never
/// touches linear motion. When `twist_loss_to_normal` is set, whatever energy
/// the twist grip removes (only when |beta_twist| < 1) is paid back into the
/// normal separation (halbach.pdf: "repulsion implies a real loss of local
/// spin... energy has transferred from spin energy to linear"). That keeps
/// energy and both momenta exact.
#[derive(Clone, Copy, Debug)]
pub struct ContactLaw {
    pub e: f64,
    pub beta: f64,
    pub beta_twist: f64,
    pub twist_loss_to_normal: bool,
    pub kappa: f64,
}

impl ContactLaw {
    pub fn rough(beta: f64) -> Self {
        ContactLaw { e: 1.0, beta, beta_twist: -1.0, twist_loss_to_normal: false, kappa: KAPPA_SOLID }
    }
}

/// Result of [`collide_full`].
#[derive(Clone, Copy, Debug, Default)]
pub struct ContactOutcome {
    /// Normal + sliding impulse on `a` (before the twist step).
    pub j: DVec3,
    /// Energy the twist step moved into normal separation.
    pub twist_to_normal: f64,
}

/// [`collide`] followed by the torsional step.
pub fn collide_full(a: &mut Ball, b: &mut Ball, n: DVec3, law: &ContactLaw) -> ContactOutcome {
    if contact_slip(a, b, n).dot(n) <= 0.0 {
        return ContactOutcome::default();
    }
    let j = collide(a, b, n, law.e, law.beta, law.kappa);
    let ia = a.inertia(law.kappa);
    let ib = b.inertia(law.kappa);
    let omega_n = (a.w - b.w).dot(n);
    let t = -(1.0 + law.beta_twist) * omega_n / (1.0 / ia + 1.0 / ib);
    let e_before = a.e_rot(law.kappa) + b.e_rot(law.kappa);
    a.w += t * n / ia;
    b.w -= t * n / ib;
    let lost = e_before - (a.e_rot(law.kappa) + b.e_rot(law.kappa));
    if law.twist_loss_to_normal && lost > 0.0 {
        let mu = a.m * b.m / (a.m + b.m);
        let w_n = (a.v - b.v).dot(n); // < 0 after the bounce: separating
        let jx = mu * (-w_n - (w_n * w_n + 2.0 * lost / mu).sqrt());
        a.v += jx * n / a.m;
        b.v -= jx * n / b.m;
        return ContactOutcome { j, twist_to_normal: lost };
    }
    ContactOutcome { j, twist_to_normal: 0.0 }
}

/// D7 option (b), PinnedC (doc §13.3; neutron.pdf "collisions affect only
/// the spins, not c"). Call after [`collide_full`] on a pair that entered at
/// |v| = 1: both velocities are projected back to |v| = 1 and both spins are
/// scaled by one common factor so the pair's rotational energy equals
/// `e_rot_before`. Energy is then exact (both linear energies are unchanged),
/// momentum is not, and the mechanics still decide the directions and how
/// spin is shared. If both spins came out exactly zero, they take
/// `fallback_axis` with half the energy each. Returns true if that fallback
/// was impossible (a zero axis).
pub fn pin_pair(a: &mut Ball, b: &mut Ball, kappa: f64, e_rot_before: f64, fallback_axis: DVec3) -> bool {
    a.v = a.v.normalize();
    b.v = b.v.normalize();
    let e_rot = a.e_rot(kappa) + b.e_rot(kappa);
    if e_rot > 1e-12 {
        let s = (e_rot_before / e_rot).sqrt();
        a.w *= s;
        b.w *= s;
        return false;
    }
    let fl = fallback_axis.length();
    if fl < 1e-9 {
        return e_rot_before > 1e-12;
    }
    let ax = fallback_axis / fl;
    a.w = ax * (e_rot_before / a.inertia(kappa)).sqrt();
    b.w = ax * (e_rot_before / b.inertia(kappa)).sqrt();
    false
}

/// Contact normal for a pair approaching with relative velocity `u`, with
/// impact parameter b/d = sqrt(rho) (uniform over the collision disk) and
/// azimuth `psi` around u. This is the geometric hard-sphere distribution,
/// with no glancing skew imposed.
pub fn contact_normal(u: DVec3, rho: f64, psi: f64) -> DVec3 {
    let uh = u.normalize();
    let helper = if uh.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let t = uh.cross(helper).normalize();
    let bp = uh.cross(t);
    let e_perp = t * psi.cos() + bp * psi.sin();
    (uh * (1.0 - rho).max(0.0).sqrt() + e_perp * rho.sqrt()).normalize()
}

/// Result of [`wall_contact`].
#[derive(Clone, Copy, Debug, Default)]
pub struct WallOutcome {
    /// Whether the contact happened (the pair was approaching along n).
    pub touched: bool,
    /// Linear momentum given to the body (= minus the photon's change).
    pub dp_body: DVec3,
    /// The photon's spin change (the stacking kick).
    pub kick: DVec3,
}

/// Doc §14: a photon meets an infinitely massive, non-spinning grain of body
/// matter moving at `grain_v`. This is the same contact law as
/// [`collide_full`] in the M -> infinity limit: mu = m, the grain adds no
/// slip term, and the twist grip acts against a non-rotating face. With
/// `pinned`, the photon's speed returns to 1 and its own spin energy is
/// restored (D7 PinnedC; the grain has no spin to share). `n` points from the
/// photon's centre toward the grain.
pub fn wall_contact(b: &mut Ball, grain_v: DVec3, n: DVec3, law: &ContactLaw, pinned: bool) -> WallOutcome {
    wall_contact_spinning(b, grain_v, DVec3::ZERO, n, law, pinned)
}

/// [`wall_contact`] against a grain that spins with angular velocity
/// `grain_w` at unit radius (doc §14.8: body matter carries its shell's spin,
/// so its surface moves at `grain_w x n`). The sliding slip gains the
/// grain's surface term. The twist grip acts on the RELATIVE spin about n:
/// (w - grain_w).n -> -beta_twist * (w - grain_w).n. A chirality-matched
/// photon therefore passes and a mismatched one catches. The grain is
/// infinitely heavy, so only the photon changes.
pub fn wall_contact_spinning(b: &mut Ball, grain_v: DVec3, grain_w: DVec3, n: DVec3, law: &ContactLaw, pinned: bool) -> WallOutcome {
    let g = (b.v - grain_v) + b.r * b.w.cross(n) + grain_w.cross(n);
    let gn = g.dot(n);
    if gn <= 0.0 {
        return WallOutcome::default();
    }
    let m = b.m;
    let inertia = b.inertia(law.kappa);
    let e_rot_in = b.e_rot(law.kappa);
    let w_in = b.w;
    let v_in = b.v;
    let g_t = g - gn * n;
    let j_n = -m * (1.0 + law.e) * gn * n;
    let k_t = (1.0 / m) * (1.0 + 1.0 / law.kappa);
    let j_t = -(1.0 + law.beta) * g_t / k_t;
    let j = j_n + j_t;
    b.v += j / m;
    let nj = n.cross(j);
    b.w += b.r * nj / inertia;
    // Twist grip on the relative spin about n.
    let omega_n = (b.w - grain_w).dot(n);
    let e_before = b.e_rot(law.kappa);
    b.w -= (1.0 + law.beta_twist) * omega_n * n;
    let lost = e_before - b.e_rot(law.kappa);
    if law.twist_loss_to_normal && lost > 0.0 {
        let w_n = (b.v - grain_v).dot(n); // < 0: separating
        let jx = m * (-w_n - (w_n * w_n + 2.0 * lost / m).sqrt());
        b.v += jx * n / m;
    }
    let mut kick = b.w - w_in;
    if pinned {
        b.v = b.v.normalize();
        let e_rot = b.e_rot(law.kappa);
        if e_rot > 1e-12 {
            let sc = (e_rot_in / e_rot).sqrt();
            // w' - w_in = (sc - 1) w_in + sc * (w - w_in), without cancellation.
            kick = (sc - 1.0) * w_in + sc * kick;
            b.w *= sc;
        } else if nj.length() > 1e-9 && e_rot_in > 1e-12 {
            b.w = nj.normalize() * (2.0 * e_rot_in / inertia).sqrt();
            kick = b.w - w_in;
        }
    }
    WallOutcome { touched: true, dp_body: -(b.v - v_in) * m, kick }
}

// ---------------------------------------------------------------------
// Spin stacking (doc §13.7, user rule 2026-09-29)
// ---------------------------------------------------------------------

/// A photon with stacked spin levels: level 0 = axial (radius 1), then x, y,
/// z at radius 2, 4, 8 (PHYSICS_REFERENCE §2). Only the OUTER level `k`
/// collides (neut2.pdf: "collisions affect only the outer spins"). Linear
/// mass is 1 at every level. The outer level's contact radius is 2^k and its
/// inertia is kappa * 1 * (2^k)^2. Inner levels are locked at edge speed c,
/// each keeping its own axis.
#[derive(Clone, Copy, Debug)]
pub struct Stacked {
    pub v: DVec3,
    /// Spin vector per level; entries above `k` are zero.
    pub w: [DVec3; 4],
    pub k: u32,
}

/// Stack/unstack settings.
#[derive(Clone, Copy, Debug)]
pub struct StackRule {
    /// Highest level index allowed (3 = the photon tier's z-spin).
    pub max_level: u32,
    /// Stack only when the edge exceeds c by more than this.
    pub stack_tol: f64,
    /// Unstack when a stacked outer level's edge falls below this (in c).
    pub unstack_edge: f64,
    pub kappa: f64,
}

impl Default for StackRule {
    fn default() -> Self {
        StackRule { max_level: 3, stack_tol: 1e-4, unstack_edge: 0.01, kappa: KAPPA_SOLID }
    }
}

/// What one [`Stacked::restack`] call did.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct RestackEvents {
    pub stacked: u32,
    pub unstacked: u32,
    /// The outer edge is still above c at `max_level` (nowhere to stack).
    pub capped: bool,
}

pub fn level_radius(k: u32) -> f64 {
    (2.0f64).powi(k as i32)
}

impl Stacked {
    pub fn new(v: DVec3, w0: DVec3) -> Self {
        Stacked { v, w: [w0, DVec3::ZERO, DVec3::ZERO, DVec3::ZERO], k: 0 }
    }

    /// The colliding body: unit mass, outer-level radius and spin.
    pub fn outer_ball(&self) -> Ball {
        Ball { v: self.v, w: self.w[self.k as usize], m: 1.0, r: level_radius(self.k) }
    }

    pub fn set_outer(&mut self, b: &Ball) {
        self.v = b.v;
        self.w[self.k as usize] = b.w;
    }

    /// Total energy: linear + every level's spin.
    pub fn energy(&self, kappa: f64) -> f64 {
        let mut e = 0.5 * self.v.length_squared();
        for i in 0..=self.k as usize {
            let r = level_radius(i as u32);
            e += 0.5 * kappa * r * r * self.w[i].length_squared();
        }
        e
    }

    /// Applies the user's rule after a collision. `dw` is the change the
    /// collision made to the outer spin (the kick that overfilled it).
    /// - Stack: while the outer edge exceeds c (and k < max), lock the outer
    ///   level at edge c, keeping its axis, and move the excess energy into a
    ///   new level k+1 at twice the radius. The new axis is the kick
    ///   direction `dw` (super.html: the new spin is centred on the collision
    ///   point), falling back to the old axis. The new spin runs at whatever
    ///   fraction of c the excess buys.
    /// - Unstack (only when nothing stacked this call): a stacked outer level
    ///   spun down below `unstack_edge` is removed. Its energy goes back into
    ///   the level below, which becomes the outer level again (limat.pdf
    ///   "stripped down").
    ///
    /// Energy is exact in both directions.
    pub fn restack(&mut self, dw: DVec3, rule: &StackRule) -> RestackEvents {
        let kappa = rule.kappa;
        let mut ev = RestackEvents::default();
        loop {
            let k = self.k as usize;
            let r = level_radius(self.k);
            let w = self.w[k];
            let edge = w.length() * r;
            if edge <= 1.0 + rule.stack_tol {
                break;
            }
            if self.k >= rule.max_level {
                ev.capped = true;
                break;
            }
            // E_cap = 0.5 * I * (1/r)^2 = 0.5 * kappa (unit mass).
            let excess = 0.5 * kappa * r * r * w.length_squared() - 0.5 * kappa;
            let axis_old = w / w.length();
            self.w[k] = axis_old / r;
            let axis_new = if dw.length() > 1e-9 { dw.normalize() } else { axis_old };
            let r_new = 2.0 * r;
            self.k += 1;
            self.w[self.k as usize] = axis_new * (2.0 * excess.max(0.0) / (kappa * r_new * r_new)).sqrt();
            ev.stacked += 1;
        }
        if ev.stacked == 0 {
            while self.k > 0 {
                let k = self.k as usize;
                let r = level_radius(self.k);
                if self.w[k].length() * r >= rule.unstack_edge {
                    break;
                }
                let e_out = 0.5 * kappa * r * r * self.w[k].length_squared();
                self.w[k] = DVec3::ZERO;
                self.k -= 1;
                let kb = self.k as usize;
                let rb = level_radius(self.k);
                let e_b = 0.5 * kappa * rb * rb * self.w[kb].length_squared() + e_out;
                let len = self.w[kb].length();
                if len > 0.0 {
                    self.w[kb] = self.w[kb] / len * (2.0 * e_b / (kappa * rb * rb)).sqrt();
                }
                ev.unstacked += 1;
            }
        }
        ev
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rng(state: &mut u64) -> f64 {
        *state ^= *state << 13;
        *state ^= *state >> 7;
        *state ^= *state << 17;
        (*state >> 11) as f64 / (1u64 << 53) as f64
    }

    fn rand_unit(s: &mut u64) -> DVec3 {
        let z = 2.0 * rng(s) - 1.0;
        let p = std::f64::consts::TAU * rng(s);
        let r = (1.0 - z * z).max(0.0).sqrt();
        DVec3::new(r * p.cos(), r * p.sin(), z)
    }

    /// RS1 (CPU half): per-collision conservation, over random states,
    /// unequal masses and levels. Linear and angular momentum are always
    /// conserved; energy only at beta = +-1. beta = 0 is the negative
    /// control, and must lose energy.
    #[test]
    fn rs1_per_collision_conservation() {
        let mut s = 0x1234_5678_9abc_def1u64;
        let kappa = KAPPA_SOLID;
        let mut max_dp: f64 = 0.0;
        let mut max_dl: f64 = 0.0;
        let mut max_de_elastic: f64 = 0.0;
        let mut min_de_beta0: f64 = 0.0;
        for i in 0..20_000 {
            let ka = (i % 3) as u32;
            let kb = ((i / 3) % 3) as u32;
            let mut a = Ball::level(ka, rand_unit(&mut s) * (0.5 + rng(&mut s)), rand_unit(&mut s) * 2.0 * rng(&mut s));
            let mut b = Ball::level(kb, rand_unit(&mut s) * (0.5 + rng(&mut s)), rand_unit(&mut s) * 2.0 * rng(&mut s));
            let u = a.v - b.v;
            if u.length() < 1e-3 {
                continue;
            }
            let n = contact_normal(u, rng(&mut s), std::f64::consts::TAU * rng(&mut s));
            let beta = [1.0, -1.0, 0.0][i % 3];
            // Place a at the origin and b at contact, so the orbital angular
            // momentum is computed from real contact geometry.
            let xa = DVec3::ZERO;
            let xb = (a.r + b.r) * n;
            let p0 = a.m * a.v + b.m * b.v;
            let l0 = xa.cross(a.m * a.v) + xb.cross(b.m * b.v) + a.inertia(kappa) * a.w + b.inertia(kappa) * b.w;
            let e0 = a.e_lin() + b.e_lin() + a.e_rot(kappa) + b.e_rot(kappa);
            let law = ContactLaw { e: 1.0, beta, beta_twist: [1.0, 0.0, 0.0][(i / 9) % 3], twist_loss_to_normal: (i / 27) % 2 == 0, kappa };
            let dissipative = beta == 0.0 || (law.beta_twist == 0.0 && !law.twist_loss_to_normal);
            collide_full(&mut a, &mut b, n, &law);
            let p1 = a.m * a.v + b.m * b.v;
            let l1 = xa.cross(a.m * a.v) + xb.cross(b.m * b.v) + a.inertia(kappa) * a.w + b.inertia(kappa) * b.w;
            let e1 = a.e_lin() + b.e_lin() + a.e_rot(kappa) + b.e_rot(kappa);
            max_dp = max_dp.max((p1 - p0).length());
            max_dl = max_dl.max((l1 - l0).length());
            if !dissipative {
                max_de_elastic = max_de_elastic.max(((e1 - e0) / e0).abs());
            } else if beta == 0.0 {
                min_de_beta0 = min_de_beta0.min(e1 - e0);
            }
        }
        println!("RS1 cpu: max|dP|={max_dp:.2e} max|dL|={max_dl:.2e} max|dE/E|(beta=+-1)={max_de_elastic:.2e} most-negative dE(beta=0)={min_de_beta0:.3e}");
        assert!(max_dp < 1e-12, "linear momentum not conserved: {max_dp}");
        assert!(max_dl < 1e-12, "angular momentum not conserved: {max_dl}");
        assert!(max_de_elastic < 1e-12, "energy not conserved at beta=+-1: {max_de_elastic}");
        assert!(min_de_beta0 < -1e-3, "negative control: beta=0 must dissipate, got {min_de_beta0}");
    }

    fn stacked_state(s: &mut u64, k: u32, edge_scale: f64) -> Stacked {
        let mut p = Stacked::new(rand_unit(s), rand_unit(s));
        for i in 0..=k as usize {
            p.w[i] = rand_unit(s) / level_radius(i as u32);
        }
        p.k = k;
        p.w[k as usize] *= edge_scale;
        p
    }

    /// Stack/unstack gates: energy exact through every event, multi-level
    /// stacks, the top cap, inner levels locked at edge c with their own
    /// axes, and a stack -> unstack round trip back to the axial level.
    #[test]
    fn stacking_conserves_energy_and_round_trips() {
        let rule = StackRule::default();
        let kappa = rule.kappa;
        let mut s = 0x57AC_u64;
        let (mut n_stack, mut n_multi, mut n_cap, mut n_unstack) = (0, 0, 0, 0);
        let mut max_de: f64 = 0.0;
        for i in 0..20_000 {
            let k = (i % 4) as u32;
            let scale = [0.005, 0.5, 1.0, 1.7, 3.0, 9.0][i % 6];
            let mut p = stacked_state(&mut s, k, scale);
            let e0 = p.energy(kappa);
            let ev = p.restack(rand_unit(&mut s), &rule);
            max_de = max_de.max(((p.energy(kappa) - e0) / e0).abs());
            n_stack += (ev.stacked > 0) as u32;
            n_multi += (ev.stacked > 1) as u32;
            n_cap += ev.capped as u32;
            n_unstack += ev.unstacked;
            for lvl in 0..p.k as usize {
                let edge = p.w[lvl].length() * level_radius(lvl as u32);
                assert!((edge - 1.0).abs() < 1e-9, "inner level {lvl} not locked at c: {edge}");
            }
            if !ev.capped {
                assert!(p.w[p.k as usize].length() * level_radius(p.k) <= 1.0 + rule.stack_tol + 1e-12);
            }
            for lvl in (p.k as usize + 1)..4 {
                assert_eq!(p.w[lvl], DVec3::ZERO);
            }
        }
        println!("stacking: stack {n_stack} (multi {n_multi}) capped {n_cap} unstack {n_unstack} max|dE/E| {max_de:.1e}");
        assert!(max_de < 1e-12);
        assert!(n_stack > 1000 && n_multi > 100 && n_cap > 100 && n_unstack > 100);

        let mut p = Stacked::new(DVec3::X, DVec3::X * 1.3); // new edge = sqrt(1.3^2 - 1) = 0.83: one stack
        let e0 = p.energy(kappa);
        let ev = p.restack(DVec3::Y, &rule);
        assert_eq!((ev.stacked, p.k), (1, 1));
        assert!((p.w[0] - DVec3::X).length() < 1e-12, "axial level keeps its axis at edge c");
        assert!(p.w[1].normalize().dot(DVec3::Y) > 0.999999, "new level takes the kick axis");
        assert!((p.energy(kappa) - e0).abs() < 1e-12);
        p.w[1] *= 0.001;
        let ev = p.restack(DVec3::ZERO, &rule);
        assert_eq!((ev.unstacked, p.k), (1, 0));
        assert!(p.w[0].normalize().dot(DVec3::X) > 0.999999);
        assert!(p.w[0].length() >= 1.0);
    }

    /// Wall contact: a non-moving grain with a perfectly elastic law conserves
    /// photon energy (beta = +-1, Free); a moving grain does work on it; the
    /// pinned version keeps |v| = 1 and the photon's spin energy.
    #[test]
    fn wall_contact_energy() {
        let mut s = 0x3A11u64;
        let law = ContactLaw::rough(1.0);
        let twist = ContactLaw { beta_twist: 0.0, twist_loss_to_normal: true, ..ContactLaw::rough(-1.0) };
        let mut max_de: f64 = 0.0;
        let mut max_moving: f64 = 0.0;
        for _ in 0..5000 {
            for l in [&law, &twist] {
                let mut b = Ball::level(0, rand_unit(&mut s), rand_unit(&mut s) * rng(&mut s));
                let e0 = b.e_lin() + b.e_rot(KAPPA_SOLID);
                let n = contact_normal(b.v, rng(&mut s), 6.28 * rng(&mut s));
                wall_contact(&mut b, DVec3::ZERO, n, l, false);
                max_de = max_de.max((b.e_lin() + b.e_rot(KAPPA_SOLID) - e0).abs() / e0);
                let mut c = Ball::level(0, rand_unit(&mut s), rand_unit(&mut s) * rng(&mut s));
                let e0 = c.e_lin() + c.e_rot(KAPPA_SOLID);
                let gv = rand_unit(&mut s) * 0.8;
                let n = contact_normal(c.v - gv, rng(&mut s), 6.28 * rng(&mut s));
                wall_contact(&mut c, gv, n, l, false);
                max_moving = max_moving.max((c.e_lin() + c.e_rot(KAPPA_SOLID) - e0).abs() / e0);
                let mut d = Ball::level(0, rand_unit(&mut s), rand_unit(&mut s) * rng(&mut s));
                let er0 = d.e_rot(KAPPA_SOLID);
                let n = contact_normal(d.v - gv, rng(&mut s), 6.28 * rng(&mut s));
                if wall_contact(&mut d, gv, n, l, true).touched {
                    assert!((d.v.length() - 1.0).abs() < 1e-12);
                    assert!((d.e_rot(KAPPA_SOLID) - er0).abs() < 1e-9 * (1.0 + er0));
                }
            }
        }
        println!("wall contact: static grain max|dE/E| {max_de:.1e}; moving grain (control) {max_moving:.2}");
        assert!(max_de < 1e-12);
        assert!(max_moving > 0.05, "a moving grain must do work on the photon");
    }

    /// Spinning grain, twist grip: a photon whose spin about n matches the
    /// grain's passes untouched; the opposite sense is locked to the grain's
    /// (reversed).
    #[test]
    fn spinning_grain_twist_matches_and_catches() {
        let twist = ContactLaw { beta_twist: 0.0, twist_loss_to_normal: true, ..ContactLaw::rough(-1.0) };
        let n = DVec3::X;
        let v = DVec3::new(1.0, 0.0, 0.0);
        let mut same = Ball::level(0, v, DVec3::X);
        let out = wall_contact_spinning(&mut same, DVec3::ZERO, DVec3::X, n, &twist, false);
        assert!(out.touched && (same.w - DVec3::X).length() < 1e-12, "matched spin must mesh");
        let mut opp = Ball::level(0, v, -DVec3::X);
        wall_contact_spinning(&mut opp, DVec3::ZERO, DVec3::X, n, &twist, false);
        // Lock against an infinitely heavy grain takes the RELATIVE spin to
        // zero: the mismatched photon is dragged all the way to the grain's
        // sense (its handedness about n is reversed).
        assert!((opp.w.dot(n) - 1.0).abs() < 1e-12, "mismatched spin must be locked to the grain's");
    }

    /// Smooth spheres (beta = -1) must leave both spins untouched: the
    /// tangential grip is the only spin<->linear coupling.
    #[test]
    fn rs_smooth_spheres_never_touch_spin() {
        let mut s = 99u64;
        for _ in 0..1000 {
            let mut a = Ball::level(0, rand_unit(&mut s), rand_unit(&mut s));
            let mut b = Ball::level(1, rand_unit(&mut s), rand_unit(&mut s));
            let (wa, wb) = (a.w, b.w);
            let u = a.v - b.v;
            if u.length() < 1e-3 {
                continue;
            }
            let n = contact_normal(u, rng(&mut s), 6.28 * rng(&mut s));
            collide(&mut a, &mut b, n, 1.0, -1.0, KAPPA_SOLID);
            assert!((a.w - wa).length() < 1e-14 && (b.w - wb).length() < 1e-14);
        }
    }

    /// Averages over the geometric impact-parameter distribution for two
    /// photons at speed c, edge speed c, with directions at angle
    /// `theta_deg` and spins set by `spin_model`.
    #[derive(Clone, Copy, Debug, Default)]
    pub struct PairStats {
        pub slip_t: f64,
        pub d_e_rot: f64,
        pub d_e_lin: f64,
        pub d_u: f64,
        pub edge_over_c: f64,
    }

    /// Axis of a photon's colliding (outer) spin.
    #[derive(Clone, Copy, Debug, PartialEq)]
    pub enum SpinModel {
        /// Axial spin about the direction of travel. chi = +1 is clockwise
        /// viewed from behind, i.e. w along +d (right-hand rule, quark.html).
        Axial,
        /// End-over-end spin with its axis perpendicular to travel at a random
        /// azimuth. It carries no helicity; this is the control.
        EndOverEnd,
    }

    pub fn pair_stats(theta_deg: f64, chi_a: f64, chi_b: f64, model: SpinModel, law: &ContactLaw, samples: usize) -> PairStats {
        let mut s = 0xDEAD_BEEF_u64 ^ (theta_deg.to_bits());
        let th = theta_deg.to_radians();
        let da = DVec3::X;
        let db = DVec3::new(th.cos(), th.sin(), 0.0);
        let mut acc = PairStats::default();
        let mut n_ok = 0.0;
        for _ in 0..samples {
            let spin_axis = |d: DVec3, chi: f64, s: &mut u64| -> DVec3 {
                match model {
                    SpinModel::Axial => chi * d,
                    SpinModel::EndOverEnd => {
                        let helper = if d.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                        let t = d.cross(helper).normalize();
                        let bp = d.cross(t);
                        let p = std::f64::consts::TAU * rng(s);
                        chi * (t * p.cos() + bp * p.sin())
                    }
                }
            };
            let wa = spin_axis(da, chi_a, &mut s);
            let wb = spin_axis(db, chi_b, &mut s);
            let mut a = Ball::level(0, da, wa);
            let mut b = Ball::level(0, db, wb);
            let u = a.v - b.v;
            let n = contact_normal(u, rng(&mut s), std::f64::consts::TAU * rng(&mut s));
            let g = contact_slip(&a, &b, n);
            let g_t = g - g.dot(n) * n;
            let (er0, el0, u0) = (a.e_rot(KAPPA_SOLID) + b.e_rot(KAPPA_SOLID), a.e_lin() + b.e_lin(), u.length());
            collide_full(&mut a, &mut b, n, law);
            acc.slip_t += g_t.length();
            acc.d_e_rot += a.e_rot(KAPPA_SOLID) + b.e_rot(KAPPA_SOLID) - er0;
            acc.d_e_lin += a.e_lin() + b.e_lin() - el0;
            acc.d_u += (a.v - b.v).length() - u0;
            acc.edge_over_c += ((a.edge_speed() > 1.0 + 1e-9) as u32 + (b.edge_speed() > 1.0 + 1e-9) as u32) as f64 * 0.5;
            n_ok += 1.0;
        }
        PairStats {
            slip_t: acc.slip_t / n_ok,
            d_e_rot: acc.d_e_rot / n_ok,
            d_e_lin: acc.d_e_lin / n_ok,
            d_u: acc.d_u / n_ok,
            edge_over_c: acc.edge_over_c / n_ok,
        }
    }

    /// RS3 gate. The sign pattern of grav4.pdf/halbach.pdf ("DirRel": the
    /// pair catches, i.e. spins down, when chi1*chi2*sign(d1.d2) < 0) is
    /// checked against each contact channel, with catch = -dE_rot averaged
    /// over impact parameters. With axial spin:
    /// - the twist (drilling) channel must reproduce DirRel;
    /// - the edge (sliding) channel, the classical rough sphere, must give
    ///   the exact inverse. Externally meshing gears counter-rotate, and two
    ///   same-helicity photons meeting head-on do counter-rotate in the lab.
    /// End-over-end spin carries no helicity, so it must show no chirality
    /// dependence at all (control).
    #[test]
    fn rs3_sign_patterns() {
        let twist = ContactLaw { beta_twist: 0.0, twist_loss_to_normal: true, ..ContactLaw::rough(-1.0) };
        let edge = ContactLaw::rough(1.0);
        for theta in [30.0, 60.0, 120.0, 150.0, 180.0] {
            let head_on = theta > 90.0;
            let catch_diff = |law: &ContactLaw, model| {
                let same = pair_stats(theta, 1.0, 1.0, model, law, 100_000);
                let opp = pair_stats(theta, 1.0, -1.0, model, law, 100_000);
                (-same.d_e_rot) - (-opp.d_e_rot)
            };
            let tw = catch_diff(&twist, SpinModel::Axial);
            let ed = catch_diff(&edge, SpinModel::Axial);
            let eoe = catch_diff(&edge, SpinModel::EndOverEnd);
            println!("theta={theta}: catch(same)-catch(opp): twist {tw:+.4}  edge {ed:+.4}  end-over-end {eoe:+.4}");
            // DirRel: same-chirality catches harder than opposite iff head-on.
            assert_eq!(tw > 0.0, head_on, "twist channel should follow DirRel at {theta} deg");
            assert!(tw.abs() > 0.01);
            assert_eq!(ed > 0.0, !head_on, "edge channel should be the inverse of DirRel at {theta} deg");
            assert!(ed.abs() > 0.01);
            assert!(eoe.abs() < 0.005, "end-over-end spin must be chirality-blind, got {eoe}");
        }
    }

    /// RS3 report: the two-photon sign table (same/opposite chirality x
    /// meeting angle) for each spin-axis model. Prints only; the gate that
    /// asserts on it is `rs3_axial_sign_pattern`.
    #[test]
    #[ignore]
    fn report_rs3_sign_table() {
        for model in [SpinModel::Axial, SpinModel::EndOverEnd] {
            let laws = [
                ("edge grip only (classical rough sphere, beta=+1)", ContactLaw::rough(1.0)),
                ("twist grip only (beta=-1, beta_twist=+1)", ContactLaw { beta_twist: 1.0, ..ContactLaw::rough(-1.0) }),
                ("twist lock, loss->normal (beta=-1, beta_twist=0)", ContactLaw { beta_twist: 0.0, twist_loss_to_normal: true, ..ContactLaw::rough(-1.0) }),
                ("edge + twist lock->normal (beta=+1, beta_twist=0)", ContactLaw { beta_twist: 0.0, twist_loss_to_normal: true, ..ContactLaw::rough(1.0) }),
            ];
            for (name, law) in laws {
                println!("\n=== RS3 {model:?} spin, {name} (speed c, edge speed c, k=0, 200k impact samples) ===");
                println!("{:>6} {:>5} | {:>9} {:>10} {:>10} {:>10} {:>8}", "angle", "chi", "slip|g_t|", "dE_rot", "dE_lin", "d|u|", "edge>c");
                println!("(d|u| > 0 = pushed apart harder = repulsion; dE_rot < 0 = spun down)");
                for theta in [30.0, 60.0, 90.0, 120.0, 150.0, 180.0] {
                    for (label, cb) in [("same", 1.0), ("opp", -1.0)] {
                        let st = pair_stats(theta, 1.0, cb, model, &law, 200_000);
                        println!(
                            "{:>6.0} {:>5} | {:>9.4} {:>+10.4} {:>+10.4} {:>+10.4} {:>8.3}",
                            theta, label, st.slip_t, st.d_e_rot, st.d_e_lin, st.d_u, st.edge_over_c
                        );
                    }
                }
            }
        }
    }
}
