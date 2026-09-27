//! GPU DSMC photon engine, phase B1 (docs/GPU_DSMC_DESIGN.md). Engine core +
//! validation gates G0 (transport) and G1 (kinetics); see the design doc for
//! the physics this ports from `recycling.rs`.

pub mod config;
pub mod engine;
pub mod gpu_types;
pub mod precompute;

pub use config::{BodyConfig, Boundary, DsmcConfig, EmissionMode, GpuClass, GpuOrient, MeshMode, RateKernel};
pub use engine::{DsmcEngine, Totals};
