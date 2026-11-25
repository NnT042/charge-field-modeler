//! The quark.html emission maze on the CPU (docs/MATHIS_SPIN_PATHS.md §A).
//!
//! Mathis (overhaul of QCD, 2008): emission leaves the surface of the core and
//! passes outward through the four stacked spins a, x, y, z. At each level it
//! turns right (+) or left (-). The sign combination decides whether it
//! escapes and whether it comes back "upside-down" (reversed spin). He shows
//! only a one-plane example and states the result of "the full analysis" as a
//! 16-state table (`MATHIS_TABLE`).
//!
//! The table depends only on the relative signs of x, y, z (a never matters):
//! escape iff x == z, upright iff y == z. This module searches a family of 3-D
//! turn rules for ones that reproduce it: turn angle, the turn axis of each
//! level, lab-fixed vs carried axes, level order, initial direction and spin,
//! and the escape / upright criteria.

use glam::{DQuat, DVec3};

/// (a, x, y, z) signs, class, escapes, returns upright.
pub const MATHIS_TABLE: [([f64; 4], &str, bool, bool); 16] = [
    ([1.0, 1.0, 1.0, 1.0], "proton", true, true),
    ([-1.0, -1.0, -1.0, -1.0], "proton", true, true),
    ([-1.0, 1.0, 1.0, 1.0], "proton", true, true),
    ([1.0, -1.0, -1.0, -1.0], "proton", true, true),
    ([-1.0, -1.0, -1.0, 1.0], "neutron", false, false),
    ([1.0, 1.0, 1.0, -1.0], "neutron", false, false),
    ([-1.0, 1.0, 1.0, -1.0], "neutron", false, false),
    ([1.0, -1.0, -1.0, 1.0], "neutron", false, false),
    ([1.0, -1.0, 1.0, -1.0], "anti-proton", true, false),
    ([-1.0, 1.0, -1.0, 1.0], "anti-proton", true, false),
    ([1.0, 1.0, -1.0, 1.0], "anti-proton", true, false),
    ([-1.0, -1.0, 1.0, -1.0], "anti-proton", true, false),
    ([1.0, -1.0, 1.0, 1.0], "anti-neutron", false, true),
    ([1.0, 1.0, -1.0, -1.0], "anti-neutron", false, true),
    ([-1.0, -1.0, 1.0, 1.0], "anti-neutron", false, true),
    ([-1.0, 1.0, -1.0, -1.0], "anti-neutron", false, true),
];

/// Classify an (escape, upright) outcome with Mathis's names.
pub fn class_of(escape: bool, upright: bool) -> &'static str {
    match (escape, upright) {
        (true, true) => "proton",
        (false, false) => "neutron",
        (true, false) => "anti-proton",
        (false, true) => "anti-neutron",
    }
}

/// One candidate turn rule.
#[derive(Clone, Copy, Debug)]
pub struct MazeRule {
    /// Turn per level (radians); the level's sign sets the sense.
    pub angle: f64,
    /// Turn axis of the a, x, y, z levels (lab frame).
    pub axes: [DVec3; 4],
    /// true: each level's axis is carried along by the turns already taken
    /// (body-frame composition); false: axes stay fixed in the lab.
    pub carried: bool,
    /// true: the emission meets z, y, x, a (outermost first) instead of
    /// a, x, y, z.
    pub outer_first: bool,
    /// Initial emission direction and photon spin.
    pub d0: DVec3,
    pub s0: DVec3,
    /// What each level does to the photon's spin (the direction always turns
    /// rigidly).
    pub spin: SpinRule,
}

/// Per-level spin interaction. `Rigid`: the spin turns with the direction
/// (helicity can never change). The others are contact effects at the
/// level, acting on the spin component n = s.axis along the level's axis
/// (unit photon spin, level spin = its sign at the same edge speed):
/// `Set`: n -> sign*|n| (the level imprints its sense);
/// `Flip`: n -> -n when n opposes the sign (a mismatched spin is reversed);
/// `Rough`: n -> 2*sign - n, the §13 twist grip at full grip against a
/// level spinning at sign (the relative spin reverses).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpinRule {
    Rigid,
    Set,
    Flip,
    Rough,
}

/// Escape criteria: how the final direction d must relate to d0.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Escape {
    /// d.d0 > 0.5: still heading out.
    Outward,
    /// d.d0 > -0.5: not turned back.
    NotBack,
    /// d == d0: exactly the original line.
    Same,
    /// d.d0 >= 0: not reversed (sideways counts as escape).
    NotReversed,
}

/// Upright criteria, fixed in advance: the final spin s against the initial
/// spin s0, the particle pole (the a-level axis) or the emission direction d0.
/// Upright is always the positive side (no polarity fitting).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Upright {
    /// s.s0 > 0.5.
    Aligned,
    /// s.s0 > -0.5.
    NotFlipped,
    /// s.s0 > 0.
    NotReversed,
    /// s.pole > 0 (pole = the a-level axis).
    PoleUp,
    /// s.d0 > 0.
    OutUp,
}

/// A result this close to a criterion threshold is decided by rounding, not
/// geometry; the rule is then rejected.
const AMBIGUOUS: f64 = 1e-6;

/// Run the maze for one sign state: final (direction, spin).
pub fn run(rule: &MazeRule, signs: &[f64]) -> (DVec3, DVec3) {
    let axes: Vec<DVec3> = (0..signs.len()).map(|i| rule.axes[i % 4]).collect();
    run_levels(rule, &axes, signs)
}

/// The maze through an explicit list of level axes (e.g. all 12 Spin-mode
/// levels), innermost first; `rule.axes` is ignored.
pub fn run_levels(rule: &MazeRule, level_axes: &[DVec3], signs: &[f64]) -> (DVec3, DVec3) {
    let n = signs.len();
    let order: Vec<usize> = if rule.outer_first { (0..n).rev().collect() } else { (0..n).collect() };
    let mut q = DQuat::IDENTITY;
    let mut s = rule.s0;
    for &i in &order {
        let base = level_axes[i];
        let axis = (if rule.carried { q * base } else { base }).normalize();
        let r = DQuat::from_axis_angle(axis, signs[i] * rule.angle);
        q = r * q;
        s = r * s;
        let n = s.dot(axis);
        let n_new = match rule.spin {
            SpinRule::Rigid => n,
            SpinRule::Set => signs[i] * n.abs(),
            SpinRule::Flip => {
                if n * signs[i] < 0.0 {
                    -n
                } else {
                    n
                }
            }
            SpinRule::Rough => 2.0 * signs[i] - n,
        };
        s += (n_new - n) * axis;
    }
    (q * rule.d0, s)
}

/// (escapes, upright), or None if any criterion sits on its threshold.
pub fn judge(rule: &MazeRule, signs: &[f64], e: Escape, u: Upright) -> Option<(bool, bool)> {
    let (d, s) = run(rule, signs);
    let dd = d.dot(rule.d0);
    let cmp = |v: f64, thr: f64| -> Option<bool> {
        if (v - thr).abs() < AMBIGUOUS {
            None
        } else {
            Some(v > thr)
        }
    };
    let esc = match e {
        Escape::Outward => cmp(dd, 0.5)?,
        Escape::NotBack => cmp(dd, -0.5)?,
        Escape::Same => (d - rule.d0).length() < AMBIGUOUS,
        Escape::NotReversed => cmp(dd, 0.0)?,
    };
    let up = match u {
        Upright::Aligned => cmp(s.dot(rule.s0), 0.5)?,
        Upright::NotFlipped => cmp(s.dot(rule.s0), -0.5)?,
        Upright::NotReversed => cmp(s.dot(rule.s0), 0.0)?,
        Upright::PoleUp => cmp(s.dot(rule.axes[0]), 0.0)?,
        Upright::OutUp => cmp(s.dot(rule.d0), 0.0)?,
    };
    Some((esc, up))
}

/// Classify a run through explicit levels with the given criteria.
pub fn judge_levels(rule: &MazeRule, level_axes: &[DVec3], signs: &[f64], e: Escape, u: Upright) -> Option<(bool, bool)> {
    let (d, s) = run_levels(rule, level_axes, signs);
    let dd = d.dot(rule.d0);
    let cmp = |v: f64, thr: f64| -> Option<bool> {
        if (v - thr).abs() < AMBIGUOUS {
            None
        } else {
            Some(v > thr)
        }
    };
    let esc = match e {
        Escape::Outward => cmp(dd, 0.5)?,
        Escape::NotBack => cmp(dd, -0.5)?,
        Escape::Same => (d - rule.d0).length() < AMBIGUOUS,
        Escape::NotReversed => cmp(dd, 0.0)?,
    };
    let up = match u {
        Upright::Aligned => cmp(s.dot(rule.s0), 0.5)?,
        Upright::NotFlipped => cmp(s.dot(rule.s0), -0.5)?,
        Upright::NotReversed => cmp(s.dot(rule.s0), 0.0)?,
        Upright::PoleUp => cmp(s.dot(rule.axes[0]), 0.0)?,
        Upright::OutUp => cmp(s.dot(rule.d0), 0.0)?,
    };
    Some((esc, up))
}

/// How many of the 16 table states a rule gets right (both escape and
/// upright), plus how many get the escape column right.
pub fn score(rule: &MazeRule, e: Escape, u: Upright) -> (usize, usize) {
    let mut both = 0;
    let mut esc_ok = 0;
    for (signs, _, esc, up) in MATHIS_TABLE.iter() {
        let Some((ge, gu)) = judge(rule, signs, e, u) else {
            return (0, 0);
        };
        if ge == *esc {
            esc_ok += 1;
            if gu == *up {
                both += 1;
            }
        }
    }
    (both, esc_ok)
}

/// The search family.
pub fn candidates() -> Vec<(MazeRule, Escape, Upright)> {
    use std::f64::consts::PI;
    let unit = [DVec3::X, DVec3::Y, DVec3::Z];
    let d0s = [
        DVec3::X,
        DVec3::Y,
        DVec3::Z,
        DVec3::new(1.0, 1.0, 0.0).normalize(),
        DVec3::new(1.0, 0.0, 1.0).normalize(),
        DVec3::new(0.0, 1.0, 1.0).normalize(),
        DVec3::new(1.0, 1.0, 1.0).normalize(),
    ];
    let mut out = Vec::new();
    for angle in [PI / 4.0, PI / 3.0, PI / 2.0, 2.0 * PI / 3.0, PI] {
        for ia in 0..3 {
            for ix in 0..3 {
                for iy in 0..3 {
                    for iz in 0..3 {
                        let axes = [unit[ia], unit[ix], unit[iy], unit[iz]];
                        for carried in [false, true] {
                            for outer_first in [false, true] {
                                for d0 in d0s {
                                    for s0 in [DVec3::X, DVec3::Y, DVec3::Z, d0] {
                                        for spin in [SpinRule::Rigid, SpinRule::Set, SpinRule::Flip, SpinRule::Rough] {
                                        let rule = MazeRule { angle, axes, carried, outer_first, d0, s0, spin };
                                        for e in [Escape::Outward, Escape::NotBack, Escape::Same, Escape::NotReversed] {
                                            for u in [Upright::Aligned, Upright::NotFlipped, Upright::NotReversed, Upright::PoleUp, Upright::OutUp] {
                                                out.push((rule, e, u));
                                            }
                                        }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    out
}

fn axis_name(v: DVec3) -> &'static str {
    if v.x.abs() > 0.9 {
        "X"
    } else if v.y.abs() > 0.9 {
        "Y"
    } else {
        "Z"
    }
}

pub fn describe(rule: &MazeRule, e: Escape, u: Upright) -> String {
    format!(
        "spin {:?}, angle {:.0} deg, axes a/x/y/z = {}/{}/{}/{}, {}, {}, d0 ({:.2},{:.2},{:.2}), s0 ({:.2},{:.2},{:.2}), escape {:?}, upright {:?}",
        rule.spin,
        rule.angle.to_degrees(),
        axis_name(rule.axes[0]),
        axis_name(rule.axes[1]),
        axis_name(rule.axes[2]),
        axis_name(rule.axes[3]),
        if rule.carried { "carried axes" } else { "lab axes" },
        if rule.outer_first { "z..a order" } else { "a..z order" },
        rule.d0.x,
        rule.d0.y,
        rule.d0.z,
        rule.s0.x,
        rule.s0.y,
        rule.s0.z,
        e,
        u
    )
}

// ---- Doc §15.9: the maze as the ball's emission law ----

/// The canonical rule's level axes (a, x, y, z) in the maze frame: Spin
/// mode's rotation axes, lab-fixed.
pub const MAZE_AXES: [DVec3; 4] = [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z];

/// A start whose table outcome is decided closer than this to a threshold is
/// rejected (f32 on the GPU).
pub const PATCH_MARGIN: f64 = 1e-4;

/// Rejection-sampling tries for the start direction before the fallback.
pub const PATCH_TRIES: u32 = 64;

/// Exact 90-degree turn of `v` about the unit axis `axis`, sense `sign`.
pub fn rot90(axis: DVec3, sign: f64, v: DVec3) -> DVec3 {
    axis * axis.dot(v) + sign * axis.cross(v)
}

/// The rigid part of the canonical maze (the direction's path): the four
/// turns a, x, y, z in order.
pub fn turn(signs: &[f64; 4], v: DVec3) -> DVec3 {
    let mut v = v;
    for i in 0..4 {
        v = rot90(MAZE_AXES[i], signs[i], v);
    }
    v
}

/// The canonical maze's spin path (Set contact spin at each level).
pub fn spin_path(signs: &[f64; 4], s0: DVec3) -> DVec3 {
    let mut s = s0;
    for i in 0..4 {
        let a = MAZE_AXES[i];
        s = rot90(a, signs[i], s);
        let n = s.dot(a);
        s += (signs[i] * n.abs() - n) * a;
    }
    s
}

/// Rotation taking unit `a` to unit `b` (shortest arc), applied to `v`.
/// Antiparallel: a half turn about an axis perpendicular to `a`.
pub fn arc(a: DVec3, b: DVec3, v: DVec3) -> DVec3 {
    let c = a.dot(b);
    if c < -0.9999 {
        let t = if a.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
        let ax = a.cross(t).normalize();
        return 2.0 * ax.dot(v) * ax - v;
    }
    let k = a.cross(b);
    v * c + k.cross(v) + k * (k.dot(v) / (1.0 + c))
}

/// True if a start `d0` with helical spin `h*d0` reproduces the whole of
/// Mathis's table under the canonical maze (with `PATCH_MARGIN`).
pub fn in_patch(d0: DVec3, h: f64) -> bool {
    let s0 = d0 * h;
    for (signs, _, esc, up) in MATHIS_TABLE.iter() {
        let e = turn(signs, d0).dot(d0);
        let u = spin_path(signs, s0).dot(s0);
        if e.abs() < PATCH_MARGIN || u.abs() < PATCH_MARGIN || (e >= 0.0) != *esc || (u >= 0.0) != *up {
            return false;
        }
    }
    true
}

/// Fallback start when rejection sampling finds none: the patch centre used
/// by the regression test, mirrored to the south for h = +1.
pub fn patch_fallback(h: f64) -> DVec3 {
    DVec3::new(1.0, 0.0, 1.0).normalize() * -h
}

/// One maze emission in the maze frame.
#[derive(Clone, Copy, Debug)]
pub struct Emission {
    /// Start direction drawn in the patch.
    pub d0: DVec3,
    /// Maze exit direction (before the trapped-photon release).
    pub d: DVec3,
    /// Outgoing direction: `d` if it escaped, else the opposite pole.
    pub out: DVec3,
    /// Outgoing spin (the maze's).
    pub spin: DVec3,
    pub escaped: bool,
    pub upright: bool,
    pub fallback: bool,
}

/// The ball's emission law (doc §15.9), in the maze frame. A captured
/// photon arriving along `dl` with outer spin `wl` is re-emitted from the
/// core: its helicity picks the patch hemisphere (h < 0: north, h >= 0:
/// south), a start d0 is drawn in the patch by rejection (`draw(i)` gives
/// uniform [0, 1) numbers, indices 0..2*PATCH_TRIES), its spin is carried
/// from dl to d0 rigidly, and it runs the canonical maze. Escaped photons
/// leave along the maze exit; trapped ones funnel to the axis and leave by
/// the opposite pole (neutron.pdf), with the maze's spin.
pub fn emit(signs: &[f64; 4], dl: DVec3, wl: DVec3, mut draw: impl FnMut(u32) -> f64) -> Emission {
    let h = if wl.dot(dl) < 0.0 { -1.0 } else { 1.0 };
    let mut d0 = patch_fallback(h);
    let mut fallback = true;
    for t in 0..PATCH_TRIES {
        let u = draw(2 * t);
        let ph = std::f64::consts::TAU * draw(2 * t + 1);
        let r = (1.0 - u * u).max(0.0).sqrt();
        let c = DVec3::new(r * ph.cos(), r * ph.sin(), -h * u);
        if in_patch(c, h) {
            d0 = c;
            fallback = false;
            break;
        }
    }
    let s0 = arc(dl, d0, wl);
    let d = turn(signs, d0);
    let spin = spin_path(signs, s0);
    let escaped = d.dot(d0) >= 0.0;
    let upright = spin.dot(s0) >= 0.0;
    let out = if escaped { d } else { DVec3::new(0.0, 0.0, h) };
    Emission { d0, d, out, spin, escaped, upright, fallback }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The table is exactly "escape iff x == z, upright iff y == z", with a
    /// irrelevant, and every one of the 16 sign states appears once.
    #[test]
    fn table_structure() {
        let mut seen = std::collections::HashSet::new();
        for (s, _, esc, up) in MATHIS_TABLE.iter() {
            assert_eq!(*esc, s[1] == s[3]);
            assert_eq!(*up, s[2] == s[3]);
            seen.insert(s.iter().map(|v| (*v > 0.0) as u8).collect::<Vec<_>>());
        }
        assert_eq!(seen.len(), 16);
    }

    /// The canonical maze rule found by the search reproduces Mathis's
    /// 16-state table with Spin mode's own axes (a and z about Z, x about X,
    /// y about Y): 90-degree turns about lab-fixed axes in the order a, x, y,
    /// z, emission starting at 45 degrees between the pole and the x axis,
    /// helical photon spin (s0 = -d0), escape = not turned back, upright =
    /// spin not reversed. Both contact spin rules (Set, Flip) work. The
    /// same rule with rigid spin (no contact) must fail: rigid turning
    /// preserves helicity and cannot tell neutron from anti-neutron.
    #[test]
    fn spin_mode_maze_reproduces_mathis_table() {
        let d0 = DVec3::new(1.0, 0.0, 1.0).normalize();
        let base = MazeRule {
            angle: std::f64::consts::FRAC_PI_2,
            axes: [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z],
            carried: false,
            outer_first: false,
            d0,
            s0: -d0,
            spin: SpinRule::Set,
        };
        for spin in [SpinRule::Set, SpinRule::Flip] {
            let (both, _) = score(&MazeRule { spin, ..base }, Escape::NotReversed, Upright::NotReversed);
            assert_eq!(both, 16, "{spin:?}");
        }
        let (rigid, _) = score(&MazeRule { spin: SpinRule::Rigid, ..base }, Escape::NotReversed, Upright::NotReversed);
        assert!(rigid < 16, "rigid turning must not reproduce the table");
    }

    /// Report: for rules that get the escape column right, what does the
    /// final spin look like per state, against several references? If no
    /// reference separates upright from upside-down, rigid turning cannot
    /// produce the upright column.
    #[test]
    #[ignore]
    fn report_maze_upright_diagnostic() {
        let cands = candidates();
        let mut n = 0;
        let mut best_by_ref = [0usize; 3];
        let names = ["s.s0", "s.d0", "s.pole"];
        let mut shown = 0;
        for (rule, e, u) in &cands {
            if *u != Upright::Aligned {
                continue; // spin criteria re-derived below
            }
            let (_, esc) = score(rule, *e, *u);
            if esc < 16 {
                continue;
            }
            n += 1;
            // For each reference r, the upright column is separable if the
            // sign of s.r splits the table exactly (either polarity).
            for (ri, best) in best_by_ref.iter_mut().enumerate() {
                let mut agree = 0;
                for (signs, _, _, up) in MATHIS_TABLE.iter() {
                    let (d, s) = run(rule, signs);
                    let _ = d;
                    let r = [rule.s0, rule.d0, rule.axes[0]][ri];
                    if (s.dot(r) > 1e-9) == *up {
                        agree += 1;
                    }
                }
                *best = (*best).max(agree);
            }
            if shown < 3 {
                shown += 1;
                println!("-- {}", describe(rule, *e, *u));
                for (signs, class, _, up) in MATHIS_TABLE.iter() {
                    let (d, s) = run(rule, signs);
                    println!("   {class:<12} {:?} up={up} d=({:+.2},{:+.2},{:+.2}) s=({:+.2},{:+.2},{:+.2})", signs, d.x, d.y, d.z, s.x, s.y, s.z);
                }
            }
        }
        println!("escape-perfect rules: {n}");
        for (i, b) in best_by_ref.iter().enumerate() {
            println!("  best upright split using sign of {}: {b}/16", names[i]);
        }
    }

    /// Report: Spin mode's own level axes (baryon tier: a = precession about
    /// Z, x about X, y about Y, z about Z; spin_stack `rotation_axis`), with
    /// the emission start direction and spin free over all 26 cube
    /// directions. Which rules reproduce the table?
    #[test]
    #[ignore]
    fn report_maze_spin_mode_axes() {
        use std::f64::consts::PI;
        let mut dirs = Vec::new();
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    if (i, j, k) != (0, 0, 0) {
                        dirs.push(DVec3::new(i as f64, j as f64, k as f64).normalize());
                    }
                }
            }
        }
        let axes = [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z];
        let mut pass = Vec::new();
        let mut tried = 0;
        for angle in [PI / 4.0, PI / 3.0, PI / 2.0, 2.0 * PI / 3.0, PI] {
            for carried in [false, true] {
                for outer_first in [false, true] {
                    for spin in [SpinRule::Rigid, SpinRule::Set, SpinRule::Flip, SpinRule::Rough] {
                        for &d0 in &dirs {
                            for s0 in dirs.iter().copied() {
                                let rule = MazeRule { angle, axes, carried, outer_first, d0, s0, spin };
                                for e in [Escape::Outward, Escape::NotBack, Escape::Same, Escape::NotReversed] {
                                    for u in [Upright::Aligned, Upright::NotFlipped, Upright::NotReversed, Upright::PoleUp, Upright::OutUp] {
                                        tried += 1;
                                        if score(&rule, e, u).0 == 16 {
                                            pass.push((rule, e, u));
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
        println!("=== Spin-mode axes Z/X/Y/Z: {tried} rules, {} reproduce the table ===", pass.len());
        for (r, e, u) in pass.iter().take(40) {
            println!("  {}", describe(r, *e, *u));
        }
    }

    /// Report: the rules that reproduce the table with Spin mode's axes, run
    /// through all 12 levels (three tiers, each level on its Spin-mode axis).
    /// Classes tested: the 4 proton states and 4 Mathis neutron states with
    /// inner tiers all +1, and the adopted neutron (z flipped in every tier:
    /// levels 4 and 8 = -1 as well).
    #[test]
    #[ignore]
    fn report_maze_twelve_levels() {
        use std::f64::consts::PI;
        let level_axes: Vec<DVec3> = (1..=12u8).map(|l| crate::spin_stack::SpinLevel::new(l).rotation_axis()).collect();
        println!("level axes 1..12: {:?}", level_axes.iter().map(|a| if a.x.abs() > 0.9 { 'X' } else if a.y.abs() > 0.9 { 'Y' } else { 'Z' }).collect::<String>());
        let mut dirs = Vec::new();
        for i in -1..=1 {
            for j in -1..=1 {
                for k in -1..=1 {
                    if (i, j, k) != (0, 0, 0) {
                        dirs.push(DVec3::new(i as f64, j as f64, k as f64).normalize());
                    }
                }
            }
        }
        let axes = [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z];
        let mut pass = Vec::new();
        for spin in [SpinRule::Set, SpinRule::Flip] {
            for &d0 in &dirs {
                for s0 in dirs.iter().copied() {
                    let rule = MazeRule { angle: PI / 2.0, axes, carried: false, outer_first: false, d0, s0, spin };
                    for u in [Upright::Aligned, Upright::NotFlipped, Upright::NotReversed, Upright::PoleUp, Upright::OutUp] {
                        if score(&rule, Escape::NotReversed, u).0 == 16 {
                            pass.push((rule, u));
                        }
                    }
                }
            }
        }
        println!("=== {} passing 4-level rules, run through 12 levels ===", pass.len());
        let tiered = |outer: [f64; 4], inner: [f64; 4]| -> Vec<f64> {
            let mut v = inner.to_vec();
            v.extend_from_slice(&inner);
            v.extend_from_slice(&outer);
            v
        };
        let groups: [(&str, Vec<Vec<f64>>); 3] = [
            ("proton (inner +1)", crate::baryon_shells::PROTON_STATES.iter().map(|s| tiered(*s, [1.0; 4])).collect()),
            ("Mathis neutron (inner +1)", crate::baryon_shells::NEUTRON_STATES.iter().map(|s| tiered(*s, [1.0; 4])).collect()),
            ("z-flip-every-tier neutron", crate::baryon_shells::NEUTRON_STATES.iter().map(|s| tiered(*s, [1.0, 1.0, 1.0, -1.0])).collect()),
        ];
        for (rule, u) in &pass {
            let mut line = format!("  {:?} s0=({:+.1},{:+.1},{:+.1}) d0=({:+.1},{:+.1},{:+.1}) {:?}:", rule.spin, rule.s0.x, rule.s0.y, rule.s0.z, rule.d0.x, rule.d0.y, rule.d0.z, u);
            for (name, states) in &groups {
                let classes: Vec<&str> = states
                    .iter()
                    .map(|sg| match judge_levels(rule, &level_axes, sg, Escape::NotReversed, *u) {
                        Some((e, up)) => class_of(e, up),
                        None => "ambiguous",
                    })
                    .collect();
                line += &format!(" | {name}: {}", classes.join(","));
            }
            println!("{line}");
        }
    }

    /// Report: the canonical rule applied to photons arriving from every
    /// direction (isotropic d_in, helical spin +-d_in), per sign state: the
    /// fractions (escape upright / escape upside-down / trapped upright /
    /// trapped upside-down), escape = not turned back relative to d_in,
    /// upright = spin not reversed relative to s_in. Does each class still
    /// sort photons the way the table says, on average?
    #[test]
    #[ignore]
    fn report_maze_isotropic() {
        let base = MazeRule {
            angle: std::f64::consts::FRAC_PI_2,
            axes: [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z],
            carried: false,
            outer_first: false,
            d0: DVec3::X,
            s0: DVec3::X,
            spin: SpinRule::Set,
        };
        let n = 20_000;
        let mut seed = 0x1234_5678_9abc_def0u64;
        let mut f = move || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut dirs = Vec::with_capacity(n);
        for _ in 0..n {
            let z = 2.0 * f() - 1.0;
            let ph = std::f64::consts::TAU * f();
            let r = (1.0 - z * z).sqrt();
            dirs.push((DVec3::new(r * ph.cos(), r * ph.sin(), z), if f() < 0.5 { 1.0 } else { -1.0 }));
        }
        for spin in [SpinRule::Set, SpinRule::Flip] {
            println!("=== {spin:?} ===");
            for (signs, class, _, _) in MATHIS_TABLE.iter() {
                let mut c = [0usize; 4];
                for (d, h) in &dirs {
                    let rule = MazeRule { d0: *d, s0: *d * *h, spin, ..base };
                    let (df, sf) = run(&rule, signs);
                    let esc = df.dot(*d) >= 0.0;
                    let up = sf.dot(rule.s0) >= 0.0;
                    c[(!esc as usize) * 2 + (!up as usize)] += 1;
                }
                let fr = |i: usize| c[i] as f64 / n as f64;
                println!("  {class:<12} {:?}: esc-up {:.3} esc-down {:.3} trap-up {:.3} trap-down {:.3}", signs, fr(0), fr(1), fr(2), fr(3));
            }
        }
    }

    /// Report: map of start directions d0 (latitude from the equator of the
    /// pole Z, azimuth from the x axis) for which the canonical rule
    /// reproduces the table, with the spin helical (s0 = -d0 and +d0).
    #[test]
    #[ignore]
    fn report_maze_start_map() {
        for spin in [SpinRule::Set, SpinRule::Flip] {
            for h in [-1.0, 1.0] {
                println!("=== {spin:?}, s0 = {h:+} d0: states right out of 16, rows = latitude -80..80 by 10, cols = azimuth 0..350 by 10 ===");
                for li in -8..=8 {
                    let lat = (li as f64 * 10.0).to_radians();
                    let mut row = format!("{:+4} ", li * 10);
                    for ai in 0..36 {
                        let az = (ai as f64 * 10.0).to_radians();
                        let d0 = DVec3::new(lat.cos() * az.cos(), lat.cos() * az.sin(), lat.sin());
                        let rule = MazeRule {
                            angle: std::f64::consts::FRAC_PI_2,
                            axes: [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z],
                            carried: false,
                            outer_first: false,
                            d0,
                            s0: d0 * h,
                            spin,
                        };
                        let (b, _) = score(&rule, Escape::NotReversed, Upright::NotReversed);
                        row += if b == 16 { "#" } else if b >= 14 { "+" } else if b >= 12 { "." } else { " " };
                    }
                    println!("{row}");
                }
            }
        }
    }

    /// Report: for starts inside the patch that reproduces the table
    /// (s0 = -d0, northern wedge around +x), the final direction of the
    /// emission per class: latitude histogram (|lat| bands 0-30 / 30-60 /
    /// 60-90), mean latitude, and the fraction ending in each hemisphere.
    #[test]
    #[ignore]
    fn report_maze_exit_directions() {
        let mut starts = Vec::new();
        for li in 1..=8 {
            for ai in -18..18 {
                let lat = (li as f64 * 10.0).to_radians();
                let az = (ai as f64 * 10.0).to_radians();
                let d0 = DVec3::new(lat.cos() * az.cos(), lat.cos() * az.sin(), lat.sin());
                let rule = MazeRule {
                    angle: std::f64::consts::FRAC_PI_2,
                    axes: [DVec3::Z, DVec3::X, DVec3::Y, DVec3::Z],
                    carried: false,
                    outer_first: false,
                    d0,
                    s0: -d0,
                    spin: SpinRule::Set,
                };
                if score(&rule, Escape::NotReversed, Upright::NotReversed).0 == 16 {
                    starts.push(rule);
                }
            }
        }
        println!("=== {} start directions in the patch ===", starts.len());
        for class in ["proton", "neutron", "anti-proton", "anti-neutron"] {
            let mut bands = [0usize; 3];
            let (mut lat_sum, mut n, mut north) = (0.0, 0usize, 0usize);
            let mut az_hist = [0usize; 8];
            for rule in &starts {
                for (signs, c, _, _) in MATHIS_TABLE.iter() {
                    if *c != class {
                        continue;
                    }
                    let (d, _) = run(rule, signs);
                    let lat = d.z.clamp(-1.0, 1.0).asin().to_degrees();
                    bands[((lat.abs() / 30.0) as usize).min(2)] += 1;
                    lat_sum += lat;
                    north += (lat > 0.0) as usize;
                    let az = d.y.atan2(d.x).to_degrees().rem_euclid(360.0);
                    az_hist[((az / 45.0) as usize).min(7)] += 1;
                    n += 1;
                }
            }
            let f = |x: usize| x as f64 / n as f64;
            println!(
                "  {class:<12} |lat| 0-30 {:.2} 30-60 {:.2} 60-90 {:.2} | mean lat {:+.1} | north {:.2} | azimuth octants from +x {:?}",
                f(bands[0]),
                f(bands[1]),
                f(bands[2]),
                lat_sum / n as f64,
                f(north),
                az_hist.iter().map(|&x| format!("{:.2}", f(x))).collect::<Vec<_>>()
            );
        }
    }

    /// Doc §15.9: the exact-turn emission helpers agree with `run()`, the
    /// patch test agrees with `score()`, both fallback centres are in the
    /// patch, and emission from the patch reproduces the table for every
    /// state (helical arrivals), with every escape within 30 degrees of the
    /// equator and every trapped photon sent out of the opposite pole.
    #[test]
    fn emission_law_matches_canonical_maze() {
        let base = MazeRule {
            angle: std::f64::consts::FRAC_PI_2,
            axes: MAZE_AXES,
            carried: false,
            outer_first: false,
            d0: DVec3::X,
            s0: DVec3::X,
            spin: SpinRule::Set,
        };
        let mut seed = 0x51u64;
        let mut f = move || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        let mut agree = 0;
        let mut inside = 0;
        for _ in 0..4000 {
            let z = 2.0 * f() - 1.0;
            let ph = std::f64::consts::TAU * f();
            let r = (1.0 - z * z).sqrt();
            let d0 = DVec3::new(r * ph.cos(), r * ph.sin(), z);
            let h = if f() < 0.5 { -1.0 } else { 1.0 };
            let rule = MazeRule { d0, s0: d0 * h, ..base };
            for (signs, _, _, _) in MATHIS_TABLE.iter() {
                let (d, s) = run(&rule, signs);
                assert!((d - turn(signs, d0)).length() < 1e-12);
                assert!((s - spin_path(signs, d0 * h)).length() < 1e-12);
            }
            let perfect = score(&rule, Escape::NotReversed, Upright::NotReversed).0 == 16;
            if perfect == in_patch(d0, h) {
                agree += 1;
            }
            inside += in_patch(d0, h) as usize;
        }
        assert!(agree >= 3995, "patch test vs score: {agree}/4000");
        assert!(inside > 40, "patch too small to sample: {inside}/4000");
        for h in [-1.0, 1.0] {
            assert!(in_patch(patch_fallback(h), h), "fallback h {h}");
        }
        for (signs, _, esc, up) in MATHIS_TABLE.iter() {
            for h in [-1.0, 1.0] {
                for _ in 0..50 {
                    let e = emit(signs, DVec3::X, DVec3::X * h, |_| f());
                    assert_eq!(e.escaped, *esc);
                    assert_eq!(e.upright, *up);
                    if e.escaped {
                        assert!(e.out.z.abs() <= 0.5 + 1e-9, "escape latitude {}", e.out.z.asin().to_degrees());
                    } else {
                        assert!((e.out - DVec3::new(0.0, 0.0, h)).length() < 1e-12);
                    }
                }
            }
        }
    }

    /// Report (doc §15.9): patch acceptance per hemisphere and the tries
    /// needed, from uniform hemisphere draws.
    #[test]
    #[ignore]
    fn report_patch_sampling() {
        let mut seed = 0x77u64;
        let mut f = move || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
            (seed >> 11) as f64 / (1u64 << 53) as f64
        };
        for h in [-1.0, 1.0] {
            let n = 200_000;
            let mut acc = 0;
            let mut sum = DVec3::ZERO;
            for _ in 0..n {
                let u = f();
                let ph = std::f64::consts::TAU * f();
                let r = (1.0 - u * u).sqrt();
                let c = DVec3::new(r * ph.cos(), r * ph.sin(), -h * u);
                if in_patch(c, h) {
                    acc += 1;
                    sum += c;
                }
            }
            let p = acc as f64 / n as f64;
            println!(
                "h {h:+}: acceptance {:.4} | P(fallback after {PATCH_TRIES} tries) {:.2e} | mean start dir {:?}",
                p,
                (1.0 - p).powi(PATCH_TRIES as i32),
                sum / acc as f64
            );
        }
    }

    /// Report: which 3-D turn rules reproduce Mathis's table.
    #[test]
    #[ignore]
    fn report_maze_search() {
        let cands = candidates();
        let mut perfect = Vec::new();
        let mut esc_perfect = 0;
        let mut hist = [0usize; 17];
        for (rule, e, u) in &cands {
            let (both, esc) = score(rule, *e, *u);
            hist[both] += 1;
            if esc == 16 {
                esc_perfect += 1;
            }
            if both == 16 {
                perfect.push((*rule, *e, *u));
            }
        }
        println!("=== maze search: {} candidate rules ===", cands.len());
        println!("score histogram (states right out of 16): {:?}", hist);
        println!("rules with the escape column right: {esc_perfect}");
        println!("rules reproducing the whole table: {}", perfect.len());
        // Group the perfect rules by their structural choices.
        let mut groups: std::collections::BTreeMap<String, usize> = std::collections::BTreeMap::new();
        for (r, e, u) in &perfect {
            let key = format!(
                "spin {:?}, angle {:.0}, axes {}{}{}{}, {}, {}, {:?}/{:?}",
                r.spin,
                r.angle.to_degrees(),
                axis_name(r.axes[0]),
                axis_name(r.axes[1]),
                axis_name(r.axes[2]),
                axis_name(r.axes[3]),
                if r.carried { "carried" } else { "lab" },
                if r.outer_first { "z-first" } else { "a-first" },
                e,
                u
            );
            *groups.entry(key).or_default() += 1;
        }
        for (k, n) in &groups {
            println!("  {k}: {n} (d0/s0 variants)");
        }
        for (r, e, u) in perfect.iter().take(12) {
            println!("  e.g. {}", describe(r, *e, *u));
        }
    }
}

#[cfg(test)]
mod tc66_probe {
    use super::*;
    /// Doc §15.66 probe: mean outgoing direction (maze frame) per arrival
    /// helicity for each baryon state, random arrival directions.
    #[test]
    #[ignore]
    fn report_maze_axial_bias() {
        let mut s: u64 = 0x9E3779B97F4A7C15;
        let mut rnd = move || { s ^= s << 13; s ^= s >> 7; s ^= s << 17; (s >> 11) as f64 / (1u64 << 53) as f64 };
        for (signs, name, _, _) in MATHIS_TABLE.iter().take(8) {
            let sg = [signs[0], signs[1], signs[2], signs[3]];
            for h in [-1.0f64, 1.0] {
                let n = 20000;
                let (mut m, mut esc, mut zesc, mut ztr) = (DVec3::ZERO, 0, 0.0, 0.0);
                let mut d0z = 0.0;
                for _ in 0..n {
                    let z = 2.0 * rnd() - 1.0; let p = std::f64::consts::TAU * rnd();
                    let r = (1.0f64 - z * z).max(0.0).sqrt();
                    let dl = DVec3::new(r * p.cos(), r * p.sin(), z);
                    let wl = dl * h;
                    let mut draws: Vec<f64> = (0..2 * PATCH_TRIES).map(|_| rnd()).collect();
                    let e = emit(&sg, dl, wl, |i| std::mem::take(&mut draws[i as usize]));
                    m += e.out; d0z += e.spin.dot(e.out);
                    if e.escaped { esc += 1; zesc += e.spin.z; } else { ztr += e.spin.dot(e.out); }
                }
                let nf = n as f64;
                println!("{name:8} {:?} h {h:+} mean out ({:+.3},{:+.3},{:+.3}) esc {:.2} out spin z {:+.3} out helicity(trapped)  {:+.3} helicity {:+.3}",
                    sg, m.x / nf, m.y / nf, m.z / nf, esc as f64 / nf, zesc / nf, ztr / nf, d0z / nf);
            }
        }
    }
}
