//! Shared types for the spin stack: spin axis kinds, tier hierarchy, and the
//! geometric mappings (level → spin type, level → tier, level → amplitude).
//!
//! Per Mathis, levels 1..=12 cycle Axial / X / Y / Z four times across three
//! tiers. Levels 1-8 are all photon spins (charge photon + high photon).
//! Levels 9-12 span electron → meson → baryon-at-rest. A hypothetical level 13
//! (baryon axial) would be the proton/neutron engine starter.

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
    pub fn label(self) -> &'static str {
        match self {
            Tier::Photon => "charge photon",
            Tier::Electron => "high photon",
            Tier::Baryon => "meson",
        }
    }
}

/// For an absolute level index 1..=12, return the (spin type, tier) pair.
/// Levels outside this range are clamped at the baryon tier.
pub fn level_spec(level: u8) -> (SpinType, Tier) {
    let tier = match level {
        1..=4 => Tier::Photon,
        5..=8 => Tier::Electron,
        _ => Tier::Baryon,
    };
    let spin_type = match (level.saturating_sub(1)) % 4 {
        0 => SpinType::Axial,
        1 => SpinType::X,
        2 => SpinType::Y,
        _ => SpinType::Z,
    };
    (spin_type, tier)
}

/// Geometric amplitude for a level (natural units, base photon r = 1).
/// Within each tier, amplitudes double: A=base, X=2×base, Y=4×base, Z=8×base.
/// Each tier's base equals the previous tier's Z amplitude, so the new axial
/// wraps the previous tier's outer extent without a gap.
///
///   Tier 1 (charge photon):  1,  2,  4,   8
///   Tier 2 (high photon):    8, 16, 32,  64
///   Tier 3 (meson):         64,128,256, 512
pub fn level_amplitude(level: u8) -> f64 {
    if level == 0 || level > 12 {
        return 0.0;
    }
    let tier = ((level - 1) / 4) as u32;
    let pos = ((level - 1) % 4) as u32;
    let tier_base = 8u64.pow(tier) as f64;
    tier_base * (1u64 << pos) as f64
}
