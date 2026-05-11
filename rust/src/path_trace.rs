//! Path trace ring buffer: records the focus particle's worldline as a
//! bounded sequence of sampled positions for visualization.
//!
//! Sampling is distance-thresholded — a new point is recorded only when the
//! particle has moved at least `min_step` from the last recorded point. This
//! keeps the line strip's vertex density roughly uniform across time scales
//! (slow motion produces sparse history, fast motion fills the buffer).
//!
//! When `capacity` points have been recorded and a new sample arrives, the
//! oldest is overwritten in classic ring-buffer fashion. `points()` returns
//! the live samples in chronological order (oldest -> newest).

use glam::DVec3;

#[derive(Clone, Copy)]
pub struct TraceSample {
    pub pos: DVec3,
    pub radius: f64,
    pub vel: DVec3,
}

pub struct PathTrace {
    buf: Vec<TraceSample>,
    /// Index of the next slot to write. Equal to `len` until the buffer fills,
    /// then wraps modulo `capacity`.
    head: usize,
    /// Number of valid samples in `buf` (saturates at `capacity`).
    len: usize,
    capacity: usize,
    min_step: f64,
    enabled: bool,
}

impl PathTrace {
    pub fn new(capacity: usize, min_step: f64) -> Self {
        Self {
            buf: Vec::with_capacity(capacity),
            head: 0,
            len: 0,
            capacity: capacity.max(1),
            min_step: min_step.max(0.0),
            enabled: true,
        }
    }

    pub fn enabled(&self) -> bool {
        self.enabled
    }

    pub fn set_enabled(&mut self, e: bool) {
        self.enabled = e;
    }

    pub fn set_min_step(&mut self, step: f64) {
        self.min_step = step.max(0.0);
    }

    #[cfg(test)]
    pub fn len(&self) -> usize {
        self.len
    }

    pub fn clear(&mut self) {
        self.buf.clear();
        self.head = 0;
        self.len = 0;
    }

    /// Record a position if tracing is enabled and the threshold is met.
    /// The very first sample is always recorded; subsequent samples must be
    /// at least `min_step` from the previous recorded point.
    pub fn record(&mut self, pos: DVec3, radius: f64, vel: DVec3) {
        if !self.enabled {
            return;
        }
        if let Some(last) = self.last() {
            if (pos - last.pos).length_squared() < self.min_step * self.min_step {
                return;
            }
        }
        let sample = TraceSample { pos, radius, vel };
        if self.buf.len() < self.capacity {
            self.buf.push(sample);
        } else {
            self.buf[self.head] = sample;
        }
        self.head = (self.head + 1) % self.capacity;
        if self.len < self.capacity {
            self.len += 1;
        }
    }

    /// Most recently recorded sample, if any.
    pub fn last(&self) -> Option<TraceSample> {
        if self.len == 0 {
            return None;
        }
        let idx = (self.head + self.capacity - 1) % self.capacity;
        Some(self.buf[idx])
    }

    /// Live samples in chronological order: oldest first, newest last.
    pub fn samples_chronological(&self) -> Vec<TraceSample> {
        if self.len == 0 {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(self.len);
        if self.len < self.capacity {
            out.extend_from_slice(&self.buf);
        } else {
            out.extend_from_slice(&self.buf[self.head..]);
            out.extend_from_slice(&self.buf[..self.head]);
        }
        out
    }

    /// Positions only (backward-compatible accessor).
    pub fn points_chronological(&self) -> Vec<DVec3> {
        self.samples_chronological().iter().map(|s| s.pos).collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const R: f64 = 1.0;

    const V: DVec3 = DVec3::ZERO;

    #[test]
    fn first_point_always_recorded() {
        let mut t = PathTrace::new(8, 0.5);
        t.record(DVec3::new(0.0, 0.0, 0.0), R, V);
        assert_eq!(t.len(), 1);
    }

    #[test]
    fn distance_threshold_skips_close_samples() {
        let mut t = PathTrace::new(8, 0.5);
        t.record(DVec3::ZERO, R, V);
        t.record(DVec3::new(0.1, 0.0, 0.0), R, V);
        t.record(DVec3::new(0.4, 0.0, 0.0), R, V);
        assert_eq!(t.len(), 1);
        t.record(DVec3::new(0.6, 0.0, 0.0), R, V);
        assert_eq!(t.len(), 2);
    }

    #[test]
    fn ring_buffer_wraps_and_preserves_chronological_order() {
        let mut t = PathTrace::new(4, 0.0);
        for i in 0..6 {
            t.record(DVec3::new(i as f64, 0.0, 0.0), R, V);
        }
        assert_eq!(t.len(), 4);
        let pts = t.points_chronological();
        assert_eq!(pts.len(), 4);
        assert_eq!(pts[0].x, 2.0);
        assert_eq!(pts[1].x, 3.0);
        assert_eq!(pts[2].x, 4.0);
        assert_eq!(pts[3].x, 5.0);
    }

    #[test]
    fn clear_resets_state() {
        let mut t = PathTrace::new(4, 0.0);
        t.record(DVec3::ZERO, R, V);
        t.record(DVec3::X, R, V);
        t.clear();
        assert_eq!(t.len(), 0);
        assert!(t.last().is_none());
        t.record(DVec3::Y, R, V);
        assert_eq!(t.points_chronological(), vec![DVec3::Y]);
    }

    #[test]
    fn disabled_drops_samples() {
        let mut t = PathTrace::new(4, 0.0);
        t.set_enabled(false);
        t.record(DVec3::ZERO, R, V);
        t.record(DVec3::X, R, V);
        assert_eq!(t.len(), 0);
        t.set_enabled(true);
        t.record(DVec3::Y, R, V);
        assert_eq!(t.len(), 1);
    }

    #[test]
    fn samples_store_radius() {
        let mut t = PathTrace::new(4, 0.0);
        t.record(DVec3::ZERO, 1.0, V);
        t.record(DVec3::X, 2.0, V);
        let samples = t.samples_chronological();
        assert_eq!(samples.len(), 2);
        assert_eq!(samples[0].radius, 1.0);
        assert_eq!(samples[1].radius, 2.0);
    }
}
