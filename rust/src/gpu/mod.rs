//! Headless GPU compute harness (feature `gpu`, off by default).
//!
//! This exists so Rust tests and reports can drive GPU compute shaders
//! *without* Godot: Godot's `RenderingDevice` API does not exist when the
//! engine is not running, so it cannot be used from `cargo test`. This module
//! wraps `wgpu` (Vulkan preferred, headless — no surface/window) instead, and
//! is infrastructure for a future GPU DSMC (Direct Simulation Monte Carlo)
//! photon-collision engine.
//!
//! Nothing here is linked into the default build or the Godot `.dll`: the
//! `wgpu`/`pollster`/`bytemuck` dependencies are all `optional = true` and
//! gated behind the `gpu` Cargo feature.
//!
//! Run the tests:
//! ```text
//! cargo test --release --manifest-path rust/Cargo.toml --features gpu gpu::
//! ```
//! Run the bench (ignored by default, prints numbers):
//! ```text
//! cargo test --release --manifest-path rust/Cargo.toml --features gpu gpu::tests::bench_particle_stream -- --ignored --nocapture
//! ```

pub mod context;
pub mod kernel;
pub mod sort;

pub use context::GpuContext;
pub use kernel::Kernel;
pub use sort::CellSorter;

#[cfg(test)]
mod tests;
