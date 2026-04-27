# Charge Field Modeler — Phase Roadmap

## Phase 1: Core Simulation (Current)

Build the foundational visualization and simulation engine. By the end of Phase 1, a user can construct a baryon from scratch, run a field simulation, and observe emergent charge behavior.

### Milestones

**M1 — Spinning Sphere** (Foundation)
- Godot project with Rust GDExtension compiling and loading
- Single PBR-shaded sphere in a 3D scene
- Axial spin slider with real-time animation
- Orbit camera with zoom/pan/rotate
- Basic readout panel (spin state, angular velocity)

**M2 — Spin Stacking: Photon Tier** (Levels 1-4)
- Quaternion composition engine for 4 spin levels
- Spin activation rules enforced in UI
- Path trace visualization (ring buffer + line strip)
- Ghost sphere overlay (transparent effective radius)
- Time scale slider (logarithmic, slow-motion to full-speed)
- Readout: effective radius, wavelength, particle classification
- Linear motion toggle with wavelength visualization

**M3 — Full Spin Stack** (Levels 5-12)
- Extend engine to 12 spin levels
- Tier transitions with auto-zoom
- Performance optimization for deep nesting
- GPU-assisted path computation if needed at high level counts
- Spin state signature display (+a+x-y+z... notation)
- Particle type auto-classification (photon/electron/proton/neutron)

**M4 — Field Simulation**
- Vulkan compute shader pipeline (update, collide, sum)
- GPU buffer management for field photons
- Point sprite rendering with color-coded chirality
- Collision detection against focus particle tier boundaries
- Impulse tracking with rolling average display
- Field presets (Earth, Solar Wind, Vacuum, Custom)
- Density slider with VRAM monitoring
- Collision statistics panel

**M5 — Multi-Particle Interactions**
- Multiple simultaneous focus particles with full spin stacks
- Inter-particle field interactions (each particle's emission affects others)
- Test scenarios:
  - Electron near proton (should be captured in orbital, too large to channel through)
  - Two protons (should orient into H2 configuration)
  - Two protons + two neutrons (should form alpha particle / He nucleus)
- Variable field conditions for each test

**M6 — Verification and Polish**
- Measure proton emission angle distribution (target: ~30° from equator)
- Verify neutron charge blocking
- Verify electron behavior in charge currents
- Calibrate unit system against known constants
- UI polish, tooltips, keyboard shortcuts
- Data export (CSV/JSON) for external analysis
- User documentation and tutorial overlay

---

## Phase 2: Nuclear Construction (Future)

Build atomic nuclei by assembling protons and neutrons according to the disc-stacking model from the nuclear papers. Verify stability predictions.

Possible features:
- Drag-and-drop nucleus builder
- Stability analysis (does the configuration hold together in the field sim?)
- Element identification from nuclear configuration
- Periodic table integration — select an element, see its predicted structure
- Charge channeling visualization through nuclear interior

---

## Phase 3: Macro-Scale Charge Field (Future)

Scale the simulation up to model charge field behavior at macro scales — planetary charge recycling, the 30° latitude charge maximum, stellar charge processing.

Possible features:
- Abstracted particle representations for large-scale fields
- Gravitational effects from charge field density gradients
- Planetary charge recycling model
- Solar charge processing model

---

## Phase 4: Experimental Sandbox (Future)

Free-form experimentation environment for testing charge field predictions.

Possible features:
- Particle accelerator mode (collide particles at controlled energies)
- Detector simulation (measure decay products)
- Conductivity simulation (charge flow through lattice structures)
- Custom experiment scripting

---

## Development Approach

This project has no deadline. Quality and accuracy take precedence over speed. Each milestone should be thoroughly tested before moving to the next. The field simulation predictions (30° emission, H2 orientation, alpha particle stability) are the ultimate validation criteria — if they don't emerge, the model needs adjustment, not shortcuts.

AI-assisted development (Claude Code, Cowork) is the primary development method. The skill file (`skills/cfm-dev/SKILL.md`) and project documents provide enough context for any session to pick up where the previous one left off.
