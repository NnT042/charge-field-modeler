//! Baked-loop hitbox for Field Calibration Mode (CM-0).
//!
//! See `docs/FIELD_CALIBRATION_MODE.md` §"The baked-loop hitbox". A particle's
//! collision shape is NOT re-composed from 12 spin levels each tick. Instead it
//! is **one baked loop + one live swing**, exploiting that `SpinStack::compose`
//! applies levels inner→outer so the outermost level is always last:
//!
//! ```text
//! final_pos = R_swing(angle) · ( baked_loop_point + swing_offset )
//! ```
//!
//! - **Baked loop:** freeze levels `1..=L` at ±c and sweep the base photon
//!   through exactly one period of level `L`. Amplitudes climb in exact powers
//!   of 2, so every inner level completes an integer number of turns within that
//!   period — the loop closes exactly and the bake is lossless.
//! - **Live swing:** the next level up (`L+1`) rotates the baked loop about the
//!   pole (Z) at a collision-driven rate (`outer_spin`). Precession swings
//!   (proton L13, electron L9) rotate in place (offset 0); orbital swings
//!   (neutron L12) offset the loop out on a string first — the "kite".
//!
//! This module is pure kinematics reused from the spin ladder; it introduces no
//! free physical constants. The collision/torque coefficients that ride on top
//! are measured by the calibration mode, not baked in here.

use std::f64::consts::TAU;

use glam::{DQuat, DVec3};

use crate::spin_stack::SpinStack;
use crate::types::{level_amplitude, level_spec, SpinType};

/// Target sample density around the fastest inner loop when auto-sizing a bake.
const POINTS_PER_INNER_REV: usize = 16;
/// Hard cap on baked polyline length (matches the GPU collision buffer budget).
const MAX_BAKE_SAMPLES: usize = 4096;

/// A particle's collision shape, baked once from its saturated spin loop.
pub struct BakedLoop {
    /// Base-photon positions over exactly one period of `loop_level`, in the
    /// swing-axis (pole) frame, natural units. Closed curve: the implicit wrap
    /// from the last point back to the first closes it (no duplicate endpoint).
    pub points: Vec<DVec3>,
    /// Base-photon orientation at each point (parallel to `points`). The base
    /// particle tumbles many times per loop, so this is where the "meets the
    /// field upright / sideways / upside-down" information lives — the raw
    /// material for the direction-relative-chirality sign rule (model A).
    pub orientations: Vec<DQuat>,
    /// Outermost level swept to form the loop (levels `1..=loop_level` at ±c).
    pub loop_level: u8,
    /// The level that swings this loop at runtime (`loop_level + 1`).
    pub swing_level: u8,
    /// Swing rotation axis (the pole). Z for our proton/neutron/electron.
    pub swing_axis: DVec3,
    /// Constant offset applied BEFORE the swing rotation. Non-zero only when the
    /// swing level is orbital (the neutron's L12 kite-on-a-string); ZERO when it
    /// is an axial precession (proton L13, electron L9).
    pub swing_offset: DVec3,
    /// True if the swing is an in-place precession (offset 0), false if orbital.
    pub swing_is_precession: bool,
}

impl BakedLoop {
    /// World-frame position of a baked loop point after swinging by `angle`
    /// (radians) about the pole. Precession → pure rotation; orbital → offset
    /// out on the string, then rotate.
    pub fn swept_point(&self, loop_point: DVec3, angle: f64) -> DVec3 {
        let rot = DQuat::from_axis_angle(self.swing_axis, angle);
        rot * (loop_point + self.swing_offset)
    }

    /// The whole baked loop swung to `angle`, as a fresh polyline. Convenience
    /// for feeding a collision pass; hot paths should call `swept_point` per
    /// point to avoid the allocation.
    pub fn swept_points(&self, angle: f64) -> Vec<DVec3> {
        let rot = DQuat::from_axis_angle(self.swing_axis, angle);
        self.points
            .iter()
            .map(|&p| rot * (p + self.swing_offset))
            .collect()
    }

    /// The base particle's spin axis at baked point `i` (loop/pole frame): the
    /// level-1 axial spin (local Y) carried by the composed orientation. This is
    /// the "face" the particle presents to the field at that point.
    pub fn spin_axis(&self, i: usize) -> DVec3 {
        self.orientations[i] * DVec3::Y
    }

    /// Distribution of the presented spin axis over the loop: `bins` buckets of
    /// cos(angle between spin axis and the pole), index 0 = cos −1 ("upside-
    /// down"), last = cos +1 ("upright"), middle = "sideways". Fractions of the
    /// loop, summing to 1. This is the zeroth-order "probability of meeting the
    /// field in each orientation" the user asked about — uniform over the route;
    /// model A will re-weight it by actual hit likelihood.
    pub fn orientation_histogram(&self, bins: usize) -> Vec<f64> {
        let bins = bins.max(1);
        let mut h = vec![0.0; bins];
        let axis = self.swing_axis.normalize_or_zero();
        let n = self.orientations.len();
        if n == 0 {
            return h;
        }
        for q in &self.orientations {
            let c = (*q * DVec3::Y).dot(axis).clamp(-1.0, 1.0);
            let idx = (((c + 1.0) * 0.5) * bins as f64).floor() as usize;
            h[idx.min(bins - 1)] += 1.0;
        }
        let inv = 1.0 / n as f64;
        for v in &mut h {
            *v *= inv;
        }
        h
    }
}

/// Recommended sample count to resolve the fan-blade fine structure of a loop
/// baked at `loop_level`. The loop winds the fastest translating orbital level
/// `orbit_radius(loop_level)` times, so we scale points with that ratio and cap
/// at the GPU buffer budget.
pub fn recommended_samples(loop_level: u8) -> usize {
    // orbit_radius of the outermost baked level = amplitude/2 for an orbital,
    // which is also the winding count of the unit-radius fastest inner orbital.
    let winds = (level_amplitude(loop_level) / 2.0).max(1.0) as usize;
    (winds * POINTS_PER_INNER_REV).clamp(8, MAX_BAKE_SAMPLES)
}

/// Geometry of the swing level (`loop_level + 1`): axis, pre-rotation offset,
/// and whether it is an in-place precession.
fn swing_geometry(loop_level: u8) -> (u8, DVec3, DVec3, bool) {
    let swing_level = loop_level + 1;
    let (swing_type, _) = level_spec(swing_level);
    // Axial levels above the photon tier precess about the pole (Z); orbital
    // levels rotate about their own axis and carry an outward offset.
    let is_precession = swing_type == SpinType::Axial && swing_level > 4;
    if is_precession {
        (swing_level, DVec3::Z, DVec3::ZERO, true)
    } else {
        let axis = swing_type.rotation_axis();
        let offset = swing_type.orbital_start() * (level_amplitude(swing_level) / 2.0);
        (swing_level, axis, offset, false)
    }
}

/// Bake the loop for a particle whose outermost saturated level is `loop_level`.
///
/// Levels `1..=loop_level` are driven to +c; the base photon is swept through
/// exactly one full period of `loop_level` and its `samples` positions recorded
/// at angles `k·(2π/samples)`, `k = 0..samples`. `samples` is clamped to at
/// least 8; pass `recommended_samples(loop_level)` for real particles.
pub fn bake_loop(loop_level: u8, samples: usize) -> BakedLoop {
    bake_loop_signed(loop_level, samples, &[])
}

/// Like `bake_loop` but with per-level chirality signs.
///
/// `signs[i]` gives the chirality for level `i+1` (1-indexed in the stack).
/// +1.0 = CW, −1.0 = CCW. Levels beyond `signs.len()` default to +1.0.
/// The magnitude is clamped to exactly 1.0 (saturated at c); only the sign
/// is used. This lets us sweep all 16 baryon chirality configurations from
/// PHYSICS_REFERENCE §3.
pub fn bake_loop_signed(loop_level: u8, samples: usize, signs: &[f64]) -> BakedLoop {
    let saturated: Vec<f64> = signs.iter().map(|&s| {
        let s = s.signum();
        if s == 0.0 { 1.0 } else { s }
    }).collect();
    bake_loop_velocities(loop_level, samples, &saturated)
}

/// Like `bake_loop_signed` but the values are raw angular velocities, not just
/// signs. `velocities[i]` is passed directly to `set_velocity` for level `i+1`
/// (clamped to [-1, 1] by the stack). Levels beyond `velocities.len()` default
/// to +1.0. Use this for sub-c spin rates like the neutron hypothesis
/// (|L12| = 0.05).
pub fn bake_loop_velocities(loop_level: u8, samples: usize, velocities: &[f64]) -> BakedLoop {
    let samples = samples.max(8);
    let loop_level = loop_level.clamp(1, 15);

    let mut stack = SpinStack::new();
    while (stack.level_count() as u8) < loop_level {
        let top = stack.level_count() as u8;
        stack.set_velocity(top, 1.0);
        stack.activate_next();
    }
    for lvl in 1..=loop_level {
        let v = velocities
            .get((lvl - 1) as usize)
            .copied()
            .unwrap_or(1.0);
        let v = if v == 0.0 { 1.0 } else { v };
        stack.set_velocity(lvl, v);
    }

    let r_outer = stack
        .get(loop_level)
        .map(|l| l.orbit_radius())
        .unwrap_or(1.0);
    let dt = r_outer * TAU / samples as f64;

    let mut points = Vec::with_capacity(samples);
    let mut orientations = Vec::with_capacity(samples);
    for _ in 0..samples {
        let (pos, orient) = stack.compose();
        points.push(pos);
        orientations.push(orient);
        stack.advance(dt, 1.0);
    }

    let (swing_level, swing_axis, swing_offset, swing_is_precession) = swing_geometry(loop_level);

    BakedLoop {
        points,
        orientations,
        loop_level,
        swing_level,
        swing_axis,
        swing_offset,
        swing_is_precession,
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const EPS: f64 = 1e-9;

    /// L=2 bake (axial + X, both at c). Level 1 axial doesn't translate, so the
    /// loop is purely the X orbital: a unit circle in the plane x = 0.
    #[test]
    fn bake_l2_is_unit_circle_in_x0_plane() {
        let baked = bake_loop(2, 8);
        assert_eq!(baked.points.len(), 8);
        for p in &baked.points {
            assert!(p.x.abs() < EPS, "point off x=0 plane: {:?}", p);
            assert!((p.length() - 1.0).abs() < EPS, "point not on unit circle: {:?}", p);
        }
        // X orbital starts along +Z (orbital_start(X) = Z).
        let p0 = baked.points[0];
        assert!((p0 - DVec3::Z).length() < EPS, "loop should start at +Z, got {:?}", p0);
    }

    /// The 4-sample L=2 loop hits the four cardinal points of the x=0 circle,
    /// confirming exact angular spacing (θ = 0, π/2, π, 3π/2 about X).
    #[test]
    fn bake_l2_spans_cardinal_points() {
        let baked = bake_loop(2, 8);
        let has = |t: DVec3| baked.points.iter().any(|p| (*p - t).length() < 1e-6);
        assert!(has(DVec3::new(0.0, 0.0, 1.0)));
        assert!(has(DVec3::new(0.0, -1.0, 0.0)));
        assert!(has(DVec3::new(0.0, 0.0, -1.0)));
        assert!(has(DVec3::new(0.0, 1.0, 0.0)));
    }

    /// Power-of-2 exactness: sweeping one extra step past the last sample must
    /// return to the first sample (the loop closes with no drift).
    #[test]
    fn baked_loop_closes_exactly() {
        // Re-run the bake manually with one extra step and compare to points[0].
        let baked = bake_loop(3, 8);
        // Rebuild the same stack and advance a full period, then one more sample
        // worth, landing back on the start.
        let mut stack = SpinStack::new();
        while (stack.level_count() as u8) < 3 {
            let top = stack.level_count() as u8;
            stack.set_velocity(top, 1.0);
            stack.activate_next();
        }
        for lvl in 1..=3 {
            stack.set_velocity(lvl, 1.0);
        }
        let r_outer = stack.get(3).unwrap().orbit_radius();
        let dt = r_outer * TAU / 8.0;
        let start = stack.compose().0;
        for _ in 0..8 {
            stack.advance(dt, 1.0);
        }
        let back = stack.compose().0;
        assert!((start - back).length() < 1e-9, "loop did not close: {:?} vs {:?}", start, back);
        assert!((baked.points[0] - start).length() < EPS);
    }

    /// Proton swing (L13) is an axial precession: swing about Z, no offset.
    #[test]
    fn proton_swing_is_precession_about_z() {
        let baked = bake_loop(12, 64);
        assert_eq!(baked.swing_level, 13);
        assert!(baked.swing_is_precession);
        assert!((baked.swing_axis - DVec3::Z).length() < EPS);
        assert!(baked.swing_offset.length() < EPS, "precession must have no offset");
    }

    /// Neutron swing (L12) is a Z-orbital: kite-on-a-string, offset X·256.
    #[test]
    fn neutron_swing_is_orbital_kite() {
        let baked = bake_loop(11, 64);
        assert_eq!(baked.swing_level, 12);
        assert!(!baked.swing_is_precession);
        assert!((baked.swing_axis - DVec3::Z).length() < EPS);
        // orbital_start(Z) = X, orbit_radius(L12) = amplitude(12)/2 = 512/2 = 256.
        let expected = DVec3::X * 256.0;
        assert!((baked.swing_offset - expected).length() < EPS,
            "neutron offset should be X·256, got {:?}", baked.swing_offset);
    }

    /// Electron swing (L9) is an axial precession: swing about Z, no offset.
    #[test]
    fn electron_swing_is_precession() {
        let baked = bake_loop(8, 64);
        assert_eq!(baked.swing_level, 9);
        assert!(baked.swing_is_precession);
        assert!(baked.swing_offset.length() < EPS);
    }

    /// A precession swing preserves each point's distance from the pole axis
    /// (pure rotation about Z), while an orbital swing translates the center.
    #[test]
    fn precession_swing_preserves_radius() {
        let baked = bake_loop(12, 64);
        let p = baked.points[3];
        for &angle in &[0.0, 0.5, 1.3, 2.7, TAU - 0.1] {
            let swung = baked.swept_point(p, angle);
            assert!((swung.length() - p.length()).abs() < 1e-6,
                "precession changed radius at angle {}", angle);
        }
    }

    #[test]
    fn orbital_swing_moves_center() {
        let baked = bake_loop(11, 64);
        // The loop's centroid should sit ~256 out and orbit the pole as it swings.
        let p = baked.points[0];
        let at0 = baked.swept_point(p, 0.0);
        let at_half = baked.swept_point(p, std::f64::consts::PI);
        // Two opposite swing phases straddle the pole: their midpoint is near Z axis.
        let mid = (at0 + at_half) * 0.5;
        assert!(mid.x.abs() < 1e-6 && mid.y.abs() < 1e-6,
            "kite midpoint should lie on the pole axis, got {:?}", mid);
        // And the point is genuinely flung out by the string.
        assert!(at0.length() > 100.0, "kite point should be far out, got {}", at0.length());
    }

    #[test]
    fn recommended_samples_scales_and_caps() {
        // L2: winds = amp(2)/2 = 1 → 16 points.
        assert_eq!(recommended_samples(2), 16);
        // L12: winds = amp(12)/2 = 256 → 256·16 = 4096, at the cap.
        assert_eq!(recommended_samples(12), MAX_BAKE_SAMPLES);
        assert!(recommended_samples(11) <= MAX_BAKE_SAMPLES);
    }

    #[test]
    fn bake_clamps_tiny_sample_counts() {
        let baked = bake_loop(2, 1);
        assert!(baked.points.len() >= 8);
    }

    #[test]
    fn orientations_captured_parallel_to_points() {
        let baked = bake_loop(3, 32);
        assert_eq!(baked.orientations.len(), baked.points.len());
        for i in 0..baked.points.len() {
            assert!((baked.spin_axis(i).length() - 1.0).abs() < 1e-9,
                "presented spin axis should be unit length at {}", i);
        }
    }

    #[test]
    fn orientation_histogram_is_a_distribution() {
        let baked = bake_loop(3, 256);
        let h = baked.orientation_histogram(16);
        assert_eq!(h.len(), 16);
        let sum: f64 = h.iter().sum();
        assert!((sum - 1.0).abs() < 1e-9, "histogram must sum to 1, got {}", sum);
        // The base particle tumbles over the loop, so orientation is spread
        // across several buckets — not parked in one.
        let occupied = h.iter().filter(|&&v| v > 1e-6).count();
        assert!(occupied >= 3, "expected a spread of orientations, got {} buckets", occupied);
    }
}
