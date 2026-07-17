use godot::prelude::*;

pub mod atom_core;
pub mod atom_scenarios;
pub mod calibration;
mod calibration_view;
pub mod charge_flow;
mod atom_sim;
mod field_sim;
mod focus_particle;
pub mod hitbox;
mod path_trace;
pub mod recycling;
mod spin_stack;
mod types;
mod units;

struct ChargeFieldModelerExtension;

#[gdextension]
unsafe impl ExtensionLibrary for ChargeFieldModelerExtension {}
