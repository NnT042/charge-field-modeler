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

/// [`collide`] followed by the torsional step. Returns the energy the twist
/// step moved into the normal separation (0 unless `twist_loss_to_normal`).
pub fn collide_full(a: &mut Ball, b: &mut Ball, n: DVec3, law: &ContactLaw) -> f64 {
    if contact_slip(a, b, n).dot(n) <= 0.0 {
        return 0.0;
    }
    collide(a, b, n, law.e, law.beta, law.kappa);
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
        return lost;
    }
    0.0
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
