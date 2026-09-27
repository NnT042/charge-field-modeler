//! Config surface for the GPU DSMC photon engine (docs/GPU_DSMC_DESIGN.md §4, §11).
//!
//! Every `[DECISION]` in the design doc is a field here rather than a code
//! branch, so a physics run picks it at construction time. B1 gates only
//! exercise a subset (see the doc's §9 gate table); the rest exist so B2/B3
//! do not need a config rewrite.

use glam::DVec3;

/// D1 — collision-rate kernel (doc §4). `Constant`: rate ∝ σ·n·c, blind to
/// relative motion (Phase-1's rule). `Catch`: rate ∝ σ·n·g with
/// g = |d1 - d2| ∈ [0,2], the standard kinetic rate.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum RateKernel {
    Constant,
    Catch,
}

/// D2 — gear mesh sign (doc §4). `ChiralityOnly`: mesh = χ1·χ2 (Phase-1,
/// direction-blind). `DirRel`: mesh = χ1·χ2·sign(d1·d2), matching the
/// `DirRelSign` rule `apply_photon` already uses for photon-particle hits.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MeshMode {
    ChiralityOnly,
    DirRel,
}

/// Per-body emission mode (doc §3.3). `Forced`: emit at the Phase-1 analogue
/// rate (lines with b <= B_MAX). `Pump`: emit at own absorption rate,
/// smoothed as an EMA over `PUMP_EMA_STEPS` steps.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EmissionMode {
    Forced,
    Pump,
}

pub const PUMP_EMA_STEPS: f64 = 50.0;

/// Nucleon seating in a stack (recycling.rs `StackOrient`, duplicated here so
/// the GPU module has no dependency on the CPU Phase-1 module).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GpuOrient {
    PoleOn,
    EdgeOn,
}

impl GpuOrient {
    pub fn spin_axis(self) -> DVec3 {
        match self {
            GpuOrient::PoleOn => DVec3::Z,
            GpuOrient::EdgeOn => DVec3::X,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GpuClass {
    Proton,
    Neutron,
    Electron,
}

impl GpuClass {
    pub fn chirality_sign(self) -> f64 {
        match self {
            GpuClass::Proton => 1.0,
            GpuClass::Neutron => -1.0,
            GpuClass::Electron => 1.0,
        }
    }

    /// (peak0_deg, peak1_deg_or_none, sigma_deg) — same numbers as
    /// recycling.rs's EXIT_PEAK_LAT_DEG / ELECTRON_EXIT_PEAK_LAT_DEG.
    pub fn exit_profile(self) -> (f64, Option<f64>, f64) {
        match self {
            GpuClass::Proton | GpuClass::Neutron => (7.0, Some(58.0), 5.0),
            GpuClass::Electron => (3.0, None, 2.5),
        }
    }

    pub fn id(self) -> u32 {
        match self {
            GpuClass::Proton => 0,
            GpuClass::Neutron => 1,
            GpuClass::Electron => 2,
        }
    }
}

/// One static body in the stack (doc §5).
#[derive(Clone, Copy, Debug)]
pub struct BodyConfig {
    pub center: DVec3,
    pub class: GpuClass,
    pub orient: GpuOrient,
    /// Through-lane radius (fraction of R_IN). 0 = lane shut for this body.
    pub through_lane_rho: f64,
}

/// The domain boundary. `SphereEscape` is the doc's normal geometry
/// (escape at R_D, centred on the body centroid). `PeriodicBox` is a G1-only
/// convenience: a cubic box with wraparound, so an ambient-only isotropic
/// population never needs injection/escape bookkeeping to hold `n_a` fixed.
#[derive(Clone, Copy, Debug)]
pub enum Boundary {
    SphereEscape { r_d: f64 },
    PeriodicBox { half_extent: f64 },
}

#[derive(Clone, Debug)]
pub struct DsmcConfig {
    pub boundary: Boundary,
    pub bodies: Vec<BodyConfig>,

    /// Cell size (world units). Doc default 0.25 for physics runs; gates may
    /// use a coarser grid for speed.
    pub h: f64,
    /// c·dt per step (c = 1 in these units).
    pub dt: f64,
    /// Max resident photon count (buffer sizing).
    pub capacity: u32,
    pub seed: u64,

    // --- physical knobs ---
    /// Ambient number density (photons per unit volume) used for injection
    /// flux and the analytic G1/G0 comparisons.
    pub n_ambient: f64,
    /// Cross-section σ (1/λ · n_a already folded out — this is the bare σ
    /// used as rate = σ·n·g). σ = 0 disables all collisions (G0).
    pub sigma: f64,
    pub photon_fraction: f64,
    pub pole_intake_bias: f64,
    pub cancel_redirect: f64,
    pub spin_transfer: f64,
    pub swing_boost: f64,

    pub rate_kernel: RateKernel,
    pub mesh_mode: MeshMode,
    pub emission_mode: EmissionMode,
    /// D4: per-photon bounce cap; 0 disables the cap.
    pub bounce_cap: u32,
    /// D5: far-zone (r > B_MAX) ambient-ambient collisions off when true.
    pub d5_far_zone_gate: bool,
    pub b_max: f64,
}

impl DsmcConfig {
    pub fn cos_max_angle() -> f64 {
        // cos(20 deg), same constant as recycling.rs::THROUGH_COS_MAX_ANGLE.
        0.939_692_620_785_908_4
    }
}
