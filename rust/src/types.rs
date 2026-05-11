//! Shared types for the spin stack: spin axis kinds, tier hierarchy, and the
//! geometric mappings (level → spin type, level → tier, level → amplitude).
//!
//! Per Mathis, levels 1..=16 cycle Axial / X / Y / Z four times across four
//! tiers. Levels 1-8 are all photon spins (charge photon + high photon).
//! Levels 9-12 span electron → meson → baryon-at-rest.
//! Levels 13-16 are suprabaryon: 13 is the proton engine starter,
//! 14-16 are unstable uberon territory (D meson and above).

use glam::{DQuat, DVec3};

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum SpinType {
    /// Rotation about the particle's own axis. No orbital offset — the frame
    /// rotates in place. By convention we use the local Y axis.
    Axial,
    X,
    Y,
    Z,
}

#[derive(Copy, Clone, Debug, PartialEq, Eq)]
pub enum Tier {
    Photon,
    Electron,
    Baryon,
    SuperBaryon,
}

impl SpinType {
    /// Rotation axis in the level's local (already-composed) frame.
    pub fn rotation_axis(self) -> DVec3 {
        match self {
            SpinType::Axial => DVec3::Y,
            SpinType::X => DVec3::X,
            SpinType::Y => DVec3::Y,
            SpinType::Z => DVec3::Z,
        }
    }

    /// Initial orbital offset direction (perpendicular to rotation_axis) in
    /// the level's local frame. The center traces a circle in the plane
    /// perpendicular to the rotation axis. Axial has no offset.
    ///
    /// X starts along Z (depth/blue axis): at angle=0 the particle is at
    /// (0, 0, amplitude) — on the ghost sphere's equatorial ring, along the
    /// depth axis, matching Wheeler's side-view convention.
    pub fn orbital_start(self) -> DVec3 {
        match self {
            SpinType::Axial => DVec3::ZERO,
            SpinType::X => DVec3::Z,
            SpinType::Y => DVec3::Z,
            SpinType::Z => DVec3::X,
        }
    }

    /// Alignment quaternion for end-over-end tumble: rotates the body's Y
    /// pole to the orbital_start direction so the pole points radially
    /// outward from the orbit center. Identity for Axial (no orbital offset).
    pub fn tumble_alignment(self) -> DQuat {
        let start = self.orbital_start();
        if start.length_squared() < 1e-9 {
            DQuat::IDENTITY
        } else {
            DQuat::from_rotation_arc(DVec3::Y, start)
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            SpinType::Axial => "axial",
            SpinType::X => "x",
            SpinType::Y => "y",
            SpinType::Z => "z",
        }
    }

    /// Single-letter tag used in spin signature strings (e.g., "+a-x+y+z").
    pub fn signature_letter(self) -> char {
        match self {
            SpinType::Axial => 'a',
            SpinType::X => 'x',
            SpinType::Y => 'y',
            SpinType::Z => 'z',
        }
    }
}

impl Tier {
    #[allow(dead_code)]
    pub fn label(self) -> &'static str {
        level_tier_label(match self {
            Tier::Photon => 1,
            Tier::Electron => 5,
            Tier::Baryon => 10,
            Tier::SuperBaryon => 13,
        })
    }
}

/// For an absolute level index 1..=16, return the (spin type, tier) pair.
/// Levels outside this range are clamped at the super-baryon tier.
pub fn level_spec(level: u8) -> (SpinType, Tier) {
    let tier = match level {
        1..=4 => Tier::Photon,
        5..=8 => Tier::Electron,
        9..=12 => Tier::Baryon,
        _ => Tier::SuperBaryon,
    };
    let spin_type = match (level.saturating_sub(1)) % 4 {
        0 => SpinType::Axial,
        1 => SpinType::X,
        2 => SpinType::Y,
        _ => SpinType::Z,
    };
    (spin_type, tier)
}

/// Level-specific tier label that accounts for the electron/meson split
/// within tier 3. Level 9 (axial) is still an electron; 10-12 are meson-class.
pub fn level_tier_label(level: u8) -> &'static str {
    match level {
        1..=4 => "charge photon",
        5..=8 => "high photon",
        9 => "electron",
        10..=12 => "meson",
        13..=16 => "baryon",
        _ => "(unknown)",
    }
}

/// Geometric amplitude for a level (natural units, base photon r = 1).
/// Within each tier, amplitudes double: A=base, X=2×base, Y=4×base, Z=8×base.
/// Each tier's base equals the previous tier's Z amplitude, so the new axial
/// wraps the previous tier's outer extent without a gap.
///
///   Tier 1 (charge photon):    1,   2,   4,    8
///   Tier 2 (high photon):      8,  16,  32,   64
///   Tier 3 (meson):           64, 128, 256,  512
///   Tier 4 (baryon):         512,1024,2048, 4096
pub fn level_amplitude(level: u8) -> f64 {
    if level == 0 || level > 16 {
        return 0.0;
    }
    let tier = ((level - 1) / 4) as u32;
    let pos = ((level - 1) % 4) as u32;
    let tier_base = 8u64.pow(tier) as f64;
    tier_base * (1u64 << pos) as f64
}
