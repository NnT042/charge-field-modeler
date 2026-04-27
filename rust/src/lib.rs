use godot::prelude::*;

mod focus_particle;
mod spin_stack;
mod types;

struct ChargeFieldModelerExtension;

#[gdextension]
unsafe impl ExtensionLibrary for ChargeFieldModelerExtension {}
