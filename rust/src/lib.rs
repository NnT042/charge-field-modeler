use godot::prelude::*;

mod focus_particle;
mod path_trace;
mod spin_stack;
mod types;
mod units;

struct ChargeFieldModelerExtension;

#[gdextension]
unsafe impl ExtensionLibrary for ChargeFieldModelerExtension {}
