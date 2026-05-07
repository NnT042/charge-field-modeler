# Charge Field Modeler

A 3D interactive physics visualization tool built on Miles Mathis's Charge Field theory. Build particles from scratch by stacking mechanical spins, then watch them interact in a real-time charge field.

## What Is This?

The Charge Field Modeler (CFM) is a hands-on simulator for a specific alternative physics framework: the idea that every particle in the universe — from photons to protons — is built from the same building block, just spinning in more and more complex ways.

You start with a single spinning sphere (a photon) and add layers of rotation one at a time. First it spins on its axis. Then it tumbles end-over-end. Then it tumbles on a second axis, then a third. Four spins make a "full" photon. Stack another four spins on top and you get an electron. Four more and you get a proton or neutron.

The simulator computes all of this mechanically — real quaternion math, real nested orbits — so the particle's path through space is an emergent consequence of its spin structure, not an animation. At high spin speeds the path trace fills in to approximate the larger sphere that an external observer would "see," demonstrating how a tiny spinning object can present as a much larger effective particle.

### Why Build This?

Standard particle physics treats electrons, quarks, and force carriers as fundamentally different objects with properties assigned by fiat. In Mathis's framework, those properties emerge from geometry: a proton is literally an electron with extra spins, and the 1820:1 mass ratio between them falls out of the energy formula for nested rotation. Whether you find that compelling or not, it's a framework that makes specific, testable geometric predictions — and this simulator is a tool for exploring them visually and mechanically.

The long-term goal is a full charge field simulation: millions of photons bombarding your constructed particle, producing emergent forces, charge channeling, and field-driven motion without hard-coding any of it.

### Features (Current Build)

- **16-Level Spin Stacking:** Build particles from a single axial photon up through electrons, baryons, and beyond — 4 tiers of 4 spins each
- **Path Tracing:** Multiple visualization modes (line, tube, surface, snake) show the actual path the base particle traces through space
- **Ghost Spheres:** Transparent overlays show the effective size of the particle at each spin level
- **Particle Classification:** The simulator identifies what you've built (B-photon, charge photon, electron, meson, proton, neutron) based on its spin state
- **Chirality Detection:** At the baryon tier, the +/- pattern of spins determines whether you've built a proton, neutron, anti-proton, or anti-neutron
- **SI Unit Display:** Wavelength, EM band classification, and spectral color computed from the spin geometry
- **Linear Motion:** Toggle linear velocity to see how the spinning structure produces wave-like motion and measurable wavelength

## Charge Field Theory in Brief

If you've wandered in here and have no idea what "charge field theory" means, here's the short version.

### The Building Block

In mainstream physics, the photon is a massless, point-like "quantum of energy." In Mathis's framework, the photon is a real, physical sphere with a measurable radius (roughly 6.67 x 10^-26 m, called the B-photon). It spins on its axis at the speed of light and travels at the speed of light. That's it — that's the only fundamental particle.

### Stacked Spins

A photon can spin in up to four independent ways:

1. **Axial (a):** Spinning on its own axis, like a top
2. **X-spin:** Tumbling end-over-end in the YZ plane
3. **Y-spin:** Tumbling end-over-end in the XZ plane
4. **Z-spin:** Tumbling end-over-end in the XY plane

Each spin level orbits at double the amplitude of the one inside it, and each must reach the speed of light before the next one can activate. A photon with all four spins at c is a "full" or "charge" photon — what we usually call visible light.

### The Particle Hierarchy

Stack three tiers of four spins and you get the familiar zoo:

| Tier | Spins | What You Get |
|------|-------|--------------|
| 1 (levels 1-4) | a, x, y, z | Photon (from infrared to gamma, depending on spin speed) |
| 2 (levels 5-8) | a2, x2, y2, z2 | Electron |
| 3 (levels 9-12) | a3, x3, y3, z3 | Proton or Neutron (determined by chirality pattern) |
| 4 (levels 13-16) | a4, x4, y4, z4 | Suprabaryon particles (D mesons, uberons) |

The electron is a photon with extra spins. The proton is an electron with extra spins. The mass ratio between proton and electron (1820:1) comes from the energy formula for nested rotation, not from any free parameter.

### The Charge Field

Electromagnetism, gravity, and the nuclear force are all mediated by the same mechanism: a field of real photons (B-photons) streaming through all matter at c. These photons collide with particles, transfer momentum, and are recycled through the spin structure. Whether charge is channeled through or deflected by a particle depends on the spin geometry — which is why protons and neutrons behave differently despite having the same structural tier.

There are no virtual particles, no exchange bosons, no gluons, and no quarks in this framework. All forces reduce to mechanical collisions between spinning spheres.

### Further Reading

For the full treatment, Miles Mathis's papers are available at [milesmathis.com](http://milesmathis.com). Some entry points relevant to what this simulator models:

- [Unifying the Electron and Proton](http://milesmathis.com/elecpro.html) — How spin stacking produces the 1820:1 mass ratio
- [What is the Charge Field?](http://milesmathis.com/charge.html) — Overview of the photon-based field theory
- [What is Pi?](http://milesmathis.com/pi2.html) — The relationship between circular motion and linear traversal

## Tech Stack

- **Godot 4** — Rendering, UI, camera, scene management
- **Rust** (via GDExtension) — Spin stack computation, quaternion math, collision logic
- **Vulkan compute shaders** — GPU-parallel field photon simulation
- **Target hardware floor:** AMD RX 550 4GB (scales up to high-end GPUs)

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
cd ../godot
godot --editor
```

## Project Documentation

- [Project Design Document](docs/PROJECT_DESIGN.md) — Full technical specification
- [Physics Reference](docs/PHYSICS_REFERENCE.md) — Charge field theory as it applies to this simulation
- [Development Skill](skills/cfm-dev/SKILL.md) — Conventions and patterns for AI-assisted development

## Status

Active development — spin stack engine and renderer complete through 16 levels (4 tiers). GPU-accelerated charge field simulation is the next major milestone.

## License

MIT

## Acknowledgments

Physics framework by Miles Mathis. All errors in implementation are ours.
