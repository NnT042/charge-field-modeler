//! GPU DSMC photon engine, phase B1 (docs/GPU_DSMC_DESIGN.md). Engine core +
//! validation gates G0 (transport) and G1 (kinetics); see the design doc for
//! the physics this ports from `recycling.rs`.

pub mod config;
pub mod density;
pub mod engine;
pub mod gpu_types;
pub mod precompute;
pub mod sigma_calib;
#[cfg(test)]
mod rough_gates;
#[cfg(test)]
mod shell_gates;

pub use config::{BodyConfig, BodyModel, Boundary, DsmcConfig, EmissionMode, ExchangeRule, GpuClass, GrainMode, GpuOrient, MeshMode, RateKernel, RoughParams, SpeedMode, StiffnessShed};
pub use density::DensityMap;
pub use engine::{DsmcEngine, SteadyBlock, SteadyReport, Totals};
pub use sigma_calib::{calibrate_sigma_for_tau_in, measure_taus, SigmaCalibration};
