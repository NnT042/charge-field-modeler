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

pub struct PathTrace {
    buf: Vec<DVec3>,
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
    pub fn record(&mut self, pos: DVec3) {
        if !self.enabled {
            return;
        }
        if let Some(last) = self.last() {
            if (pos - last).length_squared() < self.min_step * self.min_step {
                return;
            }
        }
        if self.buf.len() < self.capacity {
            self.buf.push(pos);
        } else {
            self.buf[self.head] = pos;
        }
        self.head = (self.head + 1) % self.capacity;
        if self.len < self.capacity {
            self.len += 1;
        }
    }

    /// Most recently recorded sample, if any.
    pub fn last(&self) -> Option<DVec3> {
        if self.len == 0 {
            return None;
        }
        // `head` points one past the newest, modulo capacity.
        let idx = (self.head + self.capacity - 1) % self.capacity;
        Some(self.buf[idx])
    }

    /// Live samples in chronological order: oldest first, newest last.
    pub fn points_chronological(&self) -> Vec<DVec3> {
        if self.len == 0 {
            return Vec::new();
        }
        let mut out = Vec::with_capacity(self.len);
        if self.len < self.capacity {
            // Buffer hasn't wrapped yet — `buf` is already chronological.
            out.extend_from_slice(&self.buf);
        } else {
            // After wrap, `head` points at the oldest sample.
            out.extend_from_slice(&self.buf[self.head..]);
            out.extend_from_slice(&self.buf[..self.head]);
        }
        out
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn first_point_always_recorded() {
        let mut t = PathTrace::new(8, 0.5);
        t.record(DVec3::new(0.0, 0.0, 0.0));
        assert_eq!(t.len(), 1);
    }

    #[test]
    fn distance_threshold_skips_close_samples() {
        let mut t = PathTrace::new(8, 0.5);
        t.record(DVec3::ZERO);
        t.record(DVec3::new(0.1, 0.0, 0.0)); // below threshold, skipped
        t.record(DVec3::new(0.4, 0.0, 0.0)); // still below threshold from origin
        assert_eq!(t.len(), 1);
        t.record(DVec3::new(0.6, 0.0, 0.0)); // 0.6 from origin, above threshold
        assert_eq!(t.len(), 2);
    }

    #[test]
    fn ring_buffer_wraps_and_preserves_chronological_order() {
        let mut t = PathTrace::new(4, 0.0);
        for i in 0..6 {
            t.record(DVec3::new(i as f64, 0.0, 0.0));
        }
        assert_eq!(t.len(), 4);
        let pts = t.points_chronological();
        // Oldest two (0.0, 1.0) should have been overwritten; we keep 2..=5.
        assert_eq!(pts.len(), 4);
        assert_eq!(pts[0].x, 2.0);
        assert_eq!(pts[1].x, 3.0);
        assert_eq!(pts[2].x, 4.0);
        assert_eq!(pts[3].x, 5.0);
    }

    #[test]
    fn clear_resets_state() {
        let mut t = PathTrace::new(4, 0.0);
        t.record(DVec3::ZERO);
        t.record(DVec3::X);
        t.clear();
        assert_eq!(t.len(), 0);
        assert!(t.last().is_none());
        // After clear, the next record should land at index 0 chronologically.
        t.record(DVec3::Y);
        assert_eq!(t.points_chronological(), vec![DVec3::Y]);
    }

    #[test]
    fn disabled_drops_samples() {
        let mut t = PathTrace::new(4, 0.0);
        t.set_enabled(false);
        t.record(DVec3::ZERO);
        t.record(DVec3::X);
        assert_eq!(t.len(), 0);
        t.set_enabled(true);
        t.record(DVec3::Y);
        assert_eq!(t.len(), 1);
    }
}
