---
name: cfm-dev
description: Development skill for the Charge Field Modeler (CFM), a Godot 4 + Rust GDExtension 3D physics visualization tool. Use this skill whenever working on the CFM project — including Rust spin stack math, Godot scene/UI work, Vulkan compute shaders for field simulation, quaternion operations, collision detection, or any task related to visualizing Miles Mathis's charge field theory with stacked spins. Also trigger when the user mentions photon spins, spin stacking, charge field simulation, B-photons, proton/neutron spin states, or the Faraday disc model of charge recycling. This skill covers the full stack: Rust GDExtension (gdext crate), Godot 4 GDScript, GLSL compute shaders, and project architecture.
---

# Charge Field Modeler — Development Skill

## Quick Orientation

You are working on the **Charge Field Modeler (CFM)**, a 3D interactive physics simulation built with:
- **Godot 4.x** — Scene management, UI, rendering, camera
- **Rust GDExtension** (gdext crate) — Spin stack math, quaternion composition, collision logic
- **Vulkan compute shaders** (GLSL) — GPU-parallel field photon simulation

**Before writing any code**, read these project documents in order:
1. `docs/PROJECT_DESIGN.md` — Full technical specification (architecture, data structures, unit system, milestones)
2. `docs/PHYSICS_REFERENCE.md` — The physics model (spin rules, particle types, collision behavior, predictions to test)
3. This file — Development conventions and patterns

The project knowledge base contains the source papers by Miles Mathis. Search it when you need to verify a physics claim or look up a specific derivation.

---

## Project Structure

```
charge-field-modeler/
├── godot/                  # Godot 4 project
│   ├── scenes/             # .tscn scene files
│   ├── scripts/            # GDScript (.gd)
│   ├── shaders/compute/    # Vulkan GLSL compute shaders
│   ├── shaders/visual/     # Godot .gdshader files
│   └── resources/presets/  # Field simulation presets (.tres)
├── rust/                   # Rust GDExtension crate
│   ├── Cargo.toml
│   └── src/                # Rust source files
├── docs/                   # Design docs and physics reference
└── skills/cfm-dev/         # This skill
```

---

## Development Conventions

### Rust (Simulation Core)

- Use the `gdext` crate (godot-rust) for GDExtension bindings
- All physics math in **f64** internally; convert to f32 only at the Godot boundary
- Quaternion library: use `glam` crate (f64 via `DMat4`, `DQuat`, `DVec3`)
- Unit system: **natural units** internally (r=1, c=1). SI conversion happens in a dedicated `units.rs` module, never inline in simulation code
- Error handling: `Result` types for initialization; simulation hot paths should not allocate or panic
- Naming: `snake_case` for functions/variables, `PascalCase` for types, match the data structure names in PROJECT_DESIGN.md exactly

### Godot (Presentation Layer)

- GDScript for UI logic, scene wiring, camera control
- Call into Rust via GDExtension class methods — never duplicate simulation logic in GDScript
- UI panels as separate scenes composed into the main scene
- Use Godot signals for UI → Rust communication (e.g., slider changed → call Rust method)
- Shaders go in `shaders/` — compute shaders are raw GLSL, visual shaders use Godot's .gdshader format

### Compute Shaders (GPU Field Sim)

- GLSL 450 targeting Vulkan
- Dispatched through Godot's `RenderingDevice` API (accessed from Rust or GDScript)
- Work group size: 256 (good default for GCN architecture / AMD RX 550)
- Buffer layout must match the `FieldPhoton` struct (64 bytes per photon, tightly packed)
- Always include bounds checking — thread ID may exceed particle count in the last work group

### Git Conventions

- Branch per milestone: `m1-spinning-sphere`, `m2-spin-stacking`, etc.
- Commit messages: imperative mood, reference milestone (`M1: Add axial spin slider`)
- Tag releases at milestone completion: `v0.1-m1`, `v0.2-m2`, etc.

---

## Key Implementation Patterns

### Quaternion Spin Composition

The core algorithm that computes the focus particle's position from its spin stack:

```
// Pseudocode — implement in Rust
fn compute_position(spin_stack: &[SpinLevel]) -> DVec3 {
    let mut position = DVec3::ZERO;
    let mut accumulated_rotation = DQuat::IDENTITY;
    
    for level in spin_stack {
        // Update this level's angle based on angular velocity
        // level.current_angle += level.angular_velocity * dt;
        
        // Create rotation quaternion for this level
        let axis = level.spin_type.axis(); // X, Y, or Z unit vector
        let rotation = DQuat::from_axis_angle(axis, level.current_angle);
        
        // Compose with all previous rotations
        accumulated_rotation = accumulated_rotation * rotation;
        
        // The offset vector points along the spin plane's normal, length = amplitude
        // For axial spin, no positional offset (rotation in place)
        // For X/Y/Z spins, offset is in the spin plane
        if level.spin_type != SpinType::Axial {
            let offset = accumulated_rotation * DVec3::new(level.amplitude, 0.0, 0.0);
            position += offset;
        }
    }
    
    position
}
```

**Critical note:** The axial spin rotates the particle in place — it doesn't change position. The X, Y, Z spins create orbital motion at the given amplitude. Each orbital motion is in the reference frame defined by all previous rotations. This is what makes the path complex: the Y-spin orbits in a plane that is itself being rotated by the X-spin, which is itself being rotated by the axial spin.

### Tier Boundary Detection

For collision detection, each tier (Photon, Electron, Baryon) has an effective radius equal to the amplitude of its outermost spin level:
- Photon tier: 8r (z-spin of levels 1-4)
- Electron tier: 128r (z-spin of levels 5-8)
- Baryon tier: 2048r (z-spin of levels 9-12)

A field photon at position P collides with the focus particle at position C when:
```
distance(P, C) < tier_effective_radius
```

Check from outermost tier inward. The first tier boundary crossed determines the interaction type.

### GPU Buffer Management

Field photon data lives in GPU storage buffers (SSBO). The pattern:

1. **Initialization:** Allocate buffer, fill with randomized photon data (positions, velocities, spins per preset)
2. **Per frame:** Dispatch compute shader, which reads and writes the same buffer in-place
3. **Collision readback:** A separate small buffer holds collision records. After compute dispatch, read this back to CPU for impulse summing. Keep this buffer small — only store aggregate data, not per-photon collision history.

Minimize CPU↔GPU data transfer. The photon positions never need to come back to CPU (they're rendered directly from the GPU buffer as point sprites). Only the collision aggregates need readback.

---

## Common Tasks

### "Add a new spin level to the UI"
1. In Godot: duplicate the spin slider scene, wire its signal to the Rust GDExtension
2. In Rust: extend `FocusParticle::add_spin_level()` to push a new `SpinLevel` onto the stack
3. Verify the activation rule: check that the previous level's angular_velocity == ±1.0

### "Change the field preset"
1. Edit or create a `.tres` resource in `godot/resources/presets/`
2. The preset stores: photon_ratio (float, 0-1 where 1.0 = all CW), net_direction (Vector3), energy_mean (float), energy_variance (float), density_default (int)
3. In Rust: `FieldSim::apply_preset()` reads these values and reinitializes the GPU buffer

### "Add a new particle type classification"
1. In Rust `types.rs`: add the variant to `ParticleType` enum
2. In `focus_particle.rs`: update `classify()` to check the spin state signature against the known states in PHYSICS_REFERENCE.md
3. In Godot: update the readout panel to display the new classification

### "Optimize performance"
- First: profile with Godot's built-in profiler and `cargo flamegraph`
- Spin computation: should be <1μs per frame for 12 levels. If not, check for unnecessary allocations
- GPU compute: increase work group size if occupancy is low (check with `radeon_top` on Linux or Radeon GPU Profiler on Windows)
- Rendering: field photon rendering should be the bottleneck. Reduce point sprite size or cull distant photons before adding complexity elsewhere
- Memory: field photon buffer is the largest allocation. If VRAM is tight, reduce particle count before reducing precision

---

## Testing Approach

### Unit Tests (Rust)
- Quaternion composition: verify that 4 levels of spin at known angles produce the expected position
- Spin activation rules: verify that level N+1 cannot activate until level N is at c
- Unit conversion: verify natural → SI → natural round-trips within floating point tolerance
- Collision detection: verify tier boundary checks with known positions

### Integration Tests
- Spin stack + rendering: visually verify that the path trace matches the expected spirograph patterns
- Field sim + collision: verify that collision records are generated when field photons cross tier boundaries
- Impulse summing: verify that a symmetric field produces near-zero net impulse

### Physics Validation (Manual / Semi-Automated)
- Proton spin state: set up a 12-level stack with proton chirality, run field sim, measure emission angle distribution. Target: peak near 30° from equator.
- Neutron spin state: same setup with neutron chirality, verify no coherent equatorial emission.
- H2 formation: place two proton-configured particles in a field, observe if they orient equator-to-equator.

---

## Pitfalls and Gotchas

1. **Don't confuse spin amplitude with effective radius.** The amplitude is the orbital radius of that spin level. The effective radius is the amplitude of the *outermost* active level. A particle with only levels 1-2 active has effective radius 2, not 1+2=3.

2. **Quaternion multiplication order matters.** In glam, `a * b` applies rotation `b` first, then `a`. For our composition, we want inner spins applied first, so accumulate as `accumulated = accumulated * new_level_rotation`.

3. **The axial spin has no positional offset.** It rotates the particle in place. Only X, Y, Z spins move the particle's center to a new position. The axial spin still matters because it rotates the *reference frame* for all subsequent spins.

4. **GPU buffer alignment.** The FieldPhoton struct must be exactly 64 bytes with proper alignment for Vulkan SSBO. Pad explicitly — don't rely on Rust struct layout matching GLSL layout.

5. **The RX 550 has limited VRAM.** Always check available VRAM before allocating field photon buffers. Display a warning in the UI if allocation would exceed 75% of available VRAM. Default to 500K photons, not millions.

6. **Spin chirality signs are relative to the particle, not the world.** +a means clockwise when viewed from the particle's north pole. When composing rotations, the "north pole" rotates with each level, so the axis of each level is in the rotated frame of all previous levels.

7. **The 30° angle is a prediction, not an input.** Never hard-code emission angles. The simulation should produce them from the spin geometry and collision mechanics. If it doesn't, the model needs adjustment, not a hack.
