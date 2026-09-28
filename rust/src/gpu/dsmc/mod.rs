//! GPU DSMC photon engine, phase B1 (docs/GPU_DSMC_DESIGN.md). Engine core +
//! validation gates G0 (transport) and G1 (kinetics); see the design doc for
//! the physics this ports from `recycling.rs`.

pub mod config;
pub mod density;
pub mod engine;
pub mod gpu_types;
pub mod precompute;
pub mod sigma_calib;

pub use config::{BodyConfig, Boundary, DsmcConfig, EmissionMode, ExchangeRule, GpuClass, GpuOrient, MeshMode, RateKernel, StiffnessShed};
pub use density::DensityMap;
pub use engine::{DsmcEngine, SteadyBlock, SteadyReport, Totals};
pub use sigma_calib::{calibrate_sigma_for_tau_in, measure_taus, SigmaCalibration};
