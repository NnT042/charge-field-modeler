//! Shared types for the spin stack: spin axis kinds, tier hierarchy, and the
//! geometric mappings (level → spin type, level → tier, level → amplitude).
//!
//! Per `docs/PHYSICS_REFERENCE.md`, levels 1..=12 cycle Axial / X / Y / Z four
//! times across three tiers (photon / electron / baryon), with orbital radius
//! doubling at each step (1, 2, 4, 8, 16, ..., 2048).

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
            Tier::Photon => "photon",
            Tier::Electron => "electron",
            Tier::Baryon => "baryon",
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

/// Geometric orbital radius for a level (natural units, base photon r = 1).
/// Doubles each level: 1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024, 2048.
pub fn level_amplitude(level: u8) -> f64 {
    if level == 0 || level > 12 {
        return 0.0;
    }
    (1u64 << (level - 1)) as f64
}
