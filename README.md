# Charge Field Modeler

A 3D interactive physics visualization tool for exploring Miles Mathis's Charge Field theory through mechanically accurate simulation of stacked photon spins.

## What Is This?

The Charge Field Modeler (CFM) lets you build particles from the ground up — starting with a single spinning photon and stacking spin levels one at a time until you've constructed electrons, protons, and neutrons. Then it lets you throw them into a simulated charge field and watch what happens.

This is not a simplified educational animation. Every rotation is a real quaternion composition, every collision is tracked and summed, and the emergent behaviors (charge emission patterns, molecular bonding, field-driven motion) arise naturally from the simulation rules rather than being pre-programmed.

### Features (Phase 1)

- **Stacked Spin Visualization:** Watch a sphere acquire up to 12 levels of nested rotation, from axial spin through end-over-end tumbling in three orthogonal planes, across three tiers (photon → electron → baryon)
- **Path Tracing:** See the actual path the base particle follows through space, from simple circles to complex spirograph-like curves that blur into solid shapes at full speed
- **Charge Field Simulation:** Surround your particle with hundreds of thousands to millions of charge photons and observe collision patterns, charge channeling, and impulse transfer
- **Particle Classification:** The simulator identifies what you've built (photon, electron, proton, neutron) based on its spin state signature
- **Field Presets:** Simulate conditions at Earth's surface, in solar wind, in vacuum, or custom environments
- **SI Unit Display:** All measurements shown in both natural units and standard SI

## Tech Stack

- **Godot 4** — Rendering, UI, camera, scene management
- **Rust** (via GDExtension) — Spin stack computation, quaternion math, collision logic
- **Vulkan compute shaders** — GPU-parallel field photon simulation
- **Target hardware floor:** AMD RX 550 4GB (scales up to high-end GPUs for larger field simulations)

## Building

### Prerequisites

- Godot 4.x (latest stable)
- Rust toolchain (stable, latest)
- Vulkan SDK
- `gdext` crate (godot-rust GDExtension bindings)

### Steps

```bash
# Clone the repository
git clone https://github.com/[username]/charge-field-modeler.git
cd charge-field-modeler

# Build the Rust GDExtension
cd rust
cargo build --release

# Open the Godot project
# (Godot will find the GDExtension library automatically via .gdextension config)
cd ../godot
godot --editor
```

## Project Documentation

- [Project Design Document](docs/PROJECT_DESIGN.md) — Full technical specification
- [Physics Reference](docs/PHYSICS_REFERENCE.md) — Charge field theory as it applies to this simulation
- [Development Skill](skills/cfm-dev/SKILL.md) — Conventions and patterns for AI-assisted development

## Theory Background

This project visualizes the physics described in Miles Mathis's charge field theory papers, available at [milesmathis.com](http://milesmathis.com). The key ideas:

- All particles are real spinning spheres with calculable radii (no point particles)
- Particles can carry up to four independent, nested spins (axial, x, y, z), each at double the amplitude of the previous
- Stacking three tiers of four spins produces the hierarchy: photon → electron → baryon
- The charge field is composed of real photons (B-photons) traveling at c, mediating all electromagnetic interactions through direct bombardment
- The difference between protons and neutrons is which spin chirality combinations allow or block charge emission
- There is no strong force, no weak force, no quarks, no gluons — only stacked spins and charge recycling

## Status

🚧 **Early development** — Currently building the foundational spin stack engine and renderer.

## License

[TBD]

## Acknowledgments

Physics framework by Miles Mathis. All errors in implementation are ours.
