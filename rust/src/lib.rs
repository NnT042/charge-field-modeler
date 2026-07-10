use godot::prelude::*;

pub mod atom_core;
pub mod atom_scenarios;
pub mod charge_flow;
mod atom_sim;
mod field_sim;
mod focus_particle;
mod path_trace;
mod spin_stack;
mod types;
mod units;

struct ChargeFieldModelerExtension;

#[gdextension]
unsafe impl ExtensionLibrary for ChargeFieldModelerExtension {}
