# Charge Field Modeler — Project Design Document

## Overview
The Charge Field Modeler (CFM) is a 3D interactive physics visualization tool built on Miles Mathis's Charge Field theory. It renders real-time mechanical simulations of photon spin stacking, from a single base photon up through electrons and baryons, and provides a field simulation environment for testing charge interactions.

This is not a toy or a simplified educational animation. The goal is mechanical accuracy: every rotation is a real nested quaternion composition, every collision is tracked and summed, and emergent behaviors (proton charge emission patterns, neutron neutrality, molecular bonding) should arise naturally from the simulation rather than being hard-coded.

---

## Architecture
**Engine:** Godot 4.x (latest stable)

**Simulation Core:** Rust via GDExtension (gdext crate)

**GPU Compute:** Vulkan compute shaders (GLSL) dispatched through Godot's RenderingDevice API

**Target Hardware Floor:** AMD RX 550 4GB (GCN architecture, Vulkan 1.1)

**Target OS:** Windows 10/11 primary, Linux secondary

### Why This Stack
- **Godot** provides the UI layer (sliders, panels, readouts), orbit camera, scene management, and Vulkan rendering pipeline without us having to build any of that from scratch.
- **Rust GDExtension** provides safe, fast, zero-overhead computation for the spin stack math, quaternion operations, and CPU-side simulation logic.
- **Vulkan compute shaders** handle the embarrassingly-parallel field simulation (millions of independent photon updates per frame) on the GPU.

### Repository Structure

```plaintext
charge-field-modeler/
├── godot/                      # Godot 4 project root
│   ├── project.godot
│   ├── scenes/
│   │   ├── main.tscn           # Main application scene
│   │   ├── ui/                 # UI scenes (panels, sliders, readouts)
│   │   ├── particles/          # Focus particle + field rendering scenes
│   │   └── camera/             # Orbit camera rig
│   ├── scripts/                # GDScript for UI logic and scene wiring
│   ├── shaders/
│   │   ├── compute/            # Vulkan compute shaders (GLSL)
│   │   │   ├── field_update.glsl
│   │   │   ├── collision_detect.glsl
│   │   │   └── impulse_sum.glsl
│   │   └── visual/             # Rendering shaders
│   │       ├── ghost_sphere.gdshader
│   │       ├── path_trace.gdshader
│   │       └── field_photon.gdshader
│   └── resources/              # Materials, textures, presets
│       └── presets/
│           ├── earth_ground.tres
│           ├── solar_wind.tres
│           └── vacuum.tres
├── rust/                       # Rust GDExtension crate
│   ├── Cargo.toml
│   └── src/
│       ├── lib.rs              # GDExtension entry point
│       ├── spin_stack.rs       # Spin level management and quaternion math
│       ├── focus_particle.rs   # Focus particle state and position computation
│       ├── field_sim.rs        # Field simulation dispatch and readback
│       ├── collision.rs        # Collision detection and impulse tracking
│       ├── units.rs            # Unit system (natural <-> SI conversion)
│       └── types.rs            # Shared data types
├── docs/
│   ├── PROJECT_DESIGN.md       # This file
│   ├── PHYSICS_REFERENCE.md    # Charge field theory reference (spin rules, etc.)
│   └── PHASE_ROADMAP.md        # Phase-by-phase development plan
├── skills/
│   └── cfm-dev/
│       └── SKILL.md            # Skill file for Claude Code/Cowork sessions
└── README.md

```

---

## Unit System

### Internal (Natural Units)
All simulation math uses natural units:

- **Base photon radius:** r = 1
- **Speed of light:** c = 1
- **Time:** 1 tick = time for light to cross 1r
Spin level amplitudes (radius of the effective sphere at each level):

| Level | Spin Type | Amplitude | Tier |
| --- | --- | --- | --- |
| 1 | Axial (a) | 1 | Photon |
| 2 | X-spin | 2 | Photon |
| 3 | Y-spin | 4 | Photon |
| 4 | Z-spin | 8 | Photon |
| 5 | Axial (a) | 8 | Electron |
| 6 | X-spin | 16 | Electron |
| 7 | Y-spin | 32 | Electron |
| 8 | Z-spin | 64 | Electron |
| 9 | Axial (a) | 64 | Baryon |
| 10 | X-spin | 128 | Baryon |
| 11 | Y-spin | 256 | Baryon |
| 12 | Z-spin | 512 | Baryon |
| 13 | Axial (a) | 512 | SuperBaryon |
| 14 | X-spin | 1024 | SuperBaryon |
| 15 | Y-spin | 2048 | SuperBaryon |
| 16 | Z-spin | 4096 | SuperBaryon |

**Note:** Each tier's axial amplitude equals the previous tier's Z amplitude, so nested tiers connect without gaps. The geometric ratio (512/64 = 8) is far from the known proton/electron mass ratio (~1820). From Mathis (elecpro.html): the 1820 ratio comes from *energy* scaling, not radius. A particle with spins a+x+y+z has energy 16385 vs 9 for a-spin only (ratio 1820.56). Geometric amplitudes define the spin *path*; the energy formula governs mass/energy readouts.

### Display (SI Units)
For UI readouts, convert using the B-photon radius as anchor:

- B-photon radius ≈ 6.67 × 10⁻²⁶ m (≈ G × proton radius)
- Proton radius ≈ 8.7 × 10⁻¹⁶ m (from Mathis)
- Electron radius ≈ proton radius / 1820
Display values: effective radius (m), wavelength (m), energy (eV or J), velocity (m/s or fraction of c).

---

## Core Data Structures

### SpinLevel

```plaintext
struct SpinLevel {
    level: u8,              // 1-16
    spin_type: SpinType,    // Axial, X, Y, Z
    tier: Tier,             // Photon, Electron, Baryon, SuperBaryon
    chirality: Chirality,   // CW (+) or CCW (-)
    angular_velocity: f64,  // current spin rate, range [0.0, 1.0] where 1.0 = c
    current_angle: f64,     // current rotation angle in radians
    amplitude: f64,         // orbital radius = 2^(level-1) in natural units
    quaternion: Quat,       // current rotation state
}

```

### FocusParticle

```plaintext
struct FocusParticle {
    spin_stack: Vec<SpinLevel>,     // active spin levels (1 to 16)
    position: Vec3,                 // world position
    velocity: Vec3,                 // linear velocity (for wave visualization)
    net_impulse: Vec3,              // accumulated impulse from field collisions
    impulse_history: RingBuffer,    // rolling window of impulse sums
    effective_radius: f64,          // current outermost spin amplitude
    particle_type: ParticleType,    // Photon, Electron, Proton, Neutron, etc.
}

```

### FieldPhoton (GPU-side, packed for compute shader)

```plaintext
// 64 bytes per photon, tightly packed for GPU buffer
struct FieldPhoton {
    position: [f32; 3],        // 12 bytes
    velocity: [f32; 3],        // 12 bytes
    outer_spin_axis: [f32; 3], // 9 bytes (normalized)
    outer_spin_speed: f32,     // 4 bytes
    chirality: u32,            // 4 bytes (0 = CW/photon, 1 = CCW/antiphoton)
    energy: f32,               // 4 bytes
    spin_level_count: u32,     // 4 bytes (how many levels this photon has)
    flags: u32,                // 4 bytes (collision state, active, etc.)
    _padding: [f32; 3],        // 12 bytes (alignment to 64)
}

```

### CollisionRecord

```plaintext
struct CollisionRecord {
    impact_point: Vec3,        // where on the focus particle's structure
    impact_tier: Tier,         // which tier boundary was hit
    incoming_velocity: Vec3,   // photon velocity at moment of impact
    incoming_energy: f32,      // photon energy
    incoming_chirality: Chirality,
    outcome: CollisionOutcome, // Deflected, Channeled, PassedThrough
    impulse_transferred: Vec3, // momentum change applied to focus particle
}

```

---

## Spin Stack Engine

### Quaternion Composition
Each frame, the focus particle's position is computed by summing every active spin level's lab-frame orbital offset:

1. Start at origin.
2. For each active spin level (1 to N):a. Update the level's `current_angle` by `angular_velocity dt`*.**b. Build this level's lab-frame quaternion from *`current_angle`* and its spin axis.**c. Add *`lab_quat orbital_start amplitude`* to the running position. (Axial levels have a zero *`orbital_start`* and contribute no offset.)*
3. *The orientation output is the composition of every level's lab-frame quaternion (separate from the position sum).*
*Each level's orbit plane is intrinsic to its spin type and ****fixed in the lab frame**** — inner levels move the orbit's center (via the running sum) but do NOT re-orient its plane. This is the "outside the gyroscopic influence" property from PHYSICS_REFERENCE §2: amplitudes double precisely so each level is decoupled from the inner one. So: axial = rotation in place, X = end-over-end in YZ plane, Y = end-over-end in XZ plane, Z = end-over-end in XY plane — regardless of what inner levels are doing.*

### *Spin Activation Rules*
- *Level 1 (axial) is always available*
- *Level N+1 becomes available only when Level N reaches angular_velocity = ±1.0 (i.e., ±c)*
- *Once activated, a level cannot be deactivated (spin energy is conserved)*
- *The UI enforces this: next level's slider appears (or un-grays) only when current level is maxed*

### *Path Tracing*
*The path trace records the focus particle's world position into a fixed-capacity ring buffer (default 2048 samples) each physics tick. Sampling is distance-thresholded — a new sample is recorded only when the particle has moved at least *`min_step`* (default 0.005 r) from the last one — so vertex density stays roughly uniform across time scales. When the buffer fills, the oldest sample is overwritten. The Rust side (*`path_trace.rs`*) owns the buffer; the Godot side (*`path_trace_line.gd`*) reads *`get_path_points()`* each frame and rebuilds an unshaded line strip with per-vertex alpha fade (oldest = transparent, newest = opaque). At slow playback the trace is a visible spirograph curve; at high speeds the strip fills in densely enough to approximate the swept volume — visually demonstrating why a spinning particle "looks like" a larger sphere. Hotkeys: *`T`* toggles, *`C`* clears.*

*The ghost sphere overlay renders a transparent sphere at the current effective radius (= amplitude of the outermost active spin level), giving the user the "from the outside" view simultaneously with the detailed internal motion.*

---

## *Field Simulation*

### *GPU Compute Pipeline*
*The field simulation runs as a Vulkan compute shader dispatched through Godot's RenderingDevice API. Each frame:*

1. ***Update pass**** (*`field_update.glsl`*): Each thread updates one photon's position based on its velocity. Boundary conditions (simulation volume) wrap or reflect photons at the edges.*
2. ***Collision pass**** (*`collision_detect.glsl`*): Each thread checks its photon against the focus particle's effective radii (one per tier). If a photon penetrates a tier boundary:*
- *Compare outer spin chirality and axis alignment*
- *Determine outcome (deflect, channel through, pass through)*
- *Write collision data to a collision buffer*
- *Update photon velocity/position based on outcome*
3. ***Impulse sum pass**** (*`impulse_sum.glsl`*): Parallel reduction over the collision buffer to sum all impulse vectors for this frame. Result is read back to CPU and applied to the focus particle's net_impulse.*

### *Scalable Particle Count*
*The field photon count is configurable at runtime:*

- ***Default:**** 500,000 (comfortable on RX 550 4GB)*
- ***Medium:**** 2,000,000 (fills ~128MB VRAM, good for RX 550 with headroom)*
- ***High:**** 5,000,000 (~320MB VRAM, mid-range cards)*
- ***Extreme:**** 50,000,000 (~3.2GB VRAM, high-end cards only)*
*The compute shader work group size is fixed (e.g., 256 threads); the dispatch count scales with particle count. The UI should display current VRAM usage and warn before allocating more than 75% of available VRAM.*

### *Field Presets*
***Earth Ground Level:***

- *Photon/antiphoton ratio: ~2:1 (67% CW, 33% CCW)*
- *Net field direction: configurable (default: +X in simulation coordinates)*
- *Energy distribution: thermal (configurable temperature parameter)*
- *Density: user-adjustable*
***Solar Wind:***

- *Higher energy, more directional*
- *Stronger photon bias*
***Vacuum (Cosmic Background):***

- *Equal photon/antiphoton ratio (50/50)*
- *Isotropic (no net direction)*
- *Low energy*
***Custom:***

- *User sets all parameters manually*

---

## *Collision Model*

### *Tier-Based Detection*
*Collision detection is hierarchical:*

1. *Check against outermost tier boundary (effective radius of highest active spin level)*
2. *If photon penetrates, check against next inner tier*
3. *Continue until photon is deflected or passes through all tiers*

### *Interaction Rules*
*At each tier boundary, the outcome depends on:*

- ***Size comparison:**** A field photon with fewer spin levels (smaller effective radius) can potentially pass through a larger structure's outer spin boundary. A photon of equal or greater size cannot.*
- ***Spin alignment:**** Matching chirality (both CW or both CCW) = repulsive deflection. Opposite chirality = partial cancellation of angular momentum, allowing deeper penetration.*
- ***Angle of incidence:**** Photons approaching along the spin axis (poles) interact differently than those approaching equatorially, due to the spin geometry.*

### *Impulse Tracking*
*Every collision transfers momentum. The impulse is:*

- *Direction: based on deflection angle*
- *Magnitude: proportional to the field photon's energy and the degree of interaction*
*Impulses are summed over a configurable window (default: 60 frames) and applied as a net force to the focus particle. This produces emergent motion:*

- *Symmetric field → no net motion (particle at rest)*
- *Asymmetric field → net drift in direction of lower field density*
- *Strong directional field → sustained acceleration (charge-driven motion)*
*A running display shows: total impulse magnitude, net impulse direction, collision rate (hits/second), and collision breakdown by tier and outcome type.*

---

## *Rendering*

### *Focus Particle*
- *Single sphere mesh with PBR material*
- *Shaded to clearly show surface curvature and spin direction*
- *Visible at all zoom levels (scales to minimum screen size at extreme zoom-out)*

### *Ghost Sphere*
- *Transparent sphere at current effective radius*
- *Additive blend, subtle grid or wireframe overlay*
- *Updates in real-time as spin levels are added*
- *Optional: show ghost spheres for each tier simultaneously (nested transparent shells)*

### *Path Trace*
- *Ring buffer of world positions rendered as a line strip or trail mesh*
- *Color encodes time (recent = bright, old = fading)*
- *At high speed, density of points creates the appearance of a solid volume*
- *Configurable trail length (number of frames retained)*

### *Field Photons*
- *Rendered as point sprites (GPU instanced billboards)*
- *Soft circular texture, not full sphere geometry*
- *Color-coded: blue = photon (CW), red = antiphoton (CCW)*
- *Brightness encodes energy level*
- *Size on screen is fixed (does not scale with distance — they're meant to convey density, not individual geometry)*

### *Camera*
- *Orbit camera with smooth zoom*
- *Auto-zoom when spin levels are added (pull back to keep effective radius visible)*
- *Manual override for zoom, pan, rotate*
- *"Follow particle" mode for when field sim produces emergent motion*
- *Quick-zoom presets: "see base particle," "see current tier," "see full structure"*

---

## *UI Layout*

### *Control Panel (Left Side)*
- ***Spin Level Controls:**** One slider per active level, showing chirality (±) and angular velocity. "Add Level" button appears when current highest level is at c. Each slider labeled with level number, spin type, and tier.*
- ***Speed Control:**** Global time scale slider (logarithmic). Range from "watch individual tumble" to "full speed" where paths blur into solid shapes.*
- ***View Toggles:**** Ghost sphere on/off, path trace on/off, nested tier shells on/off, field photon visibility on/off.*

### *Readout Panel (Right Side)*
- ***Current State:**** Effective radius (natural + SI), effective wavelength (natural + SI), particle classification (photon/electron/proton/neutron/unknown), spin state signature (e.g., +a+x+y-z...).*
- ***Field Sim Stats:**** Collision rate, net impulse vector, impulse magnitude, collision breakdown pie chart.*
- ***Performance:**** FPS, VRAM usage, active field photon count.*

### *Field Sim Controls (Bottom)*
- ***Preset selector:**** Earth, Solar Wind, Vacuum, Custom*
- ***Density slider:**** Photon count (with VRAM warning)*
- ***Ratio slider:**** Photon/antiphoton balance*
- ***Direction controls:**** Net field direction*
- ***Start/Stop/Reset buttons***

### *Linear Motion Toggle (Top)*
- *Switch between "at rest" view (particle spins in place, Einstein's thought experiment) and "in motion" view (particle translates linearly, showing wave behavior)*
- *When in motion: velocity slider and direction control*
- *Wavelength readout updates in real-time based on outermost spin + linear velocity*

---

## *Phase 1 Milestones*

### *M1: Spinning Sphere*
- [x] *Godot project scaffold with Rust GDExtension building and loading*
- [x] *Single sphere rendered with PBR shading*
- [x] *Axial spin slider (±c range) with real-time rotation*
- [x] *Basic orbit camera*

### *M2: Spin Stacking (Levels 1-4, Photon Tier)*
- [x] *Quaternion-based spin composition for 4 levels*
- [x] *Level activation rules (previous must be at c)*
- [x] *Path trace visualization*
- [x] *Ghost sphere overlay*
- [x] *Speed control (time scale slider)*
- [x] *Readout panel: radius, wavelength, classification (wavelength deferred to SI units module)*

### *M3: Full Spin Stack + Particle Zoo*
- [x] *Extend spin stack to 12 levels (3 tiers)*
- [x] *Tier-tabbed slider panel with auto-activate*
- [x] *SI units module (wavelength, EM band classification)*
- [x] *EM spectrum readout with visible-spectrum trace tint*
- [x] *Auto-zoom on level addition*
- [x] *Particle classification (B-photon → charge photon → high photon → electron → meson → muon → baryon)*

### *M3.5: Baryon Extension (Levels 13-16)*
- [x] *Extend spin stack to 16 levels (4 tiers, SuperBaryon tier)*
- [x] *Baryon tab with dark red/brown ghost sphere palette*
- [x] *Chirality-based classification (proton/neutron/anti-proton/anti-neutron from ± sign pattern)*
- [x] *D meson and uberon labels for suprabaryon levels*
- [x] *Performance: adaptive trace min_step scaling with effective radius*
- [x] *Performance: ArrayMesh bulk upload replaces per-vertex ImmediateMesh calls*
- [x] *Spin annul button (0) — truncate stack from any level*
- [x] *Snake trace mode — quarter-rotation fading tube tail*
- [ ] *GPU-accelerated position computation (for future field sim at scale)*

### *M4: Field Simulation*
- [ ] *GPU compute pipeline for field photon updates*
- [ ] *Collision detection against focus particle*
- [ ] *Impulse tracking and net force computation*
- [ ] *Field presets (Earth, Solar Wind, Vacuum, Custom)*
- [ ] *Point sprite rendering for field photons*
- [ ] *Collision statistics display*

### *M5: Multi-Particle Interactions*
- [ ] *Support multiple focus particles (each with full spin stack)*
- [ ] *Inter-particle field interactions*
- [ ] *Test scenarios: electron near proton, two protons, proton + neutron pairs*
- [ ] *Alpha particle (helium nucleus) formation test*
- [ ] *H2 molecule orientation test*

### *M6: Polish and Verification*
- [ ] *Compare proton emission pattern against 30° prediction*
- [ ] *Verify neutron charge blocking behavior*
- [ ] *Calibrate unit system against known physical constants*
- [ ] *UI polish, tooltips, help overlay*
- [ ] *Export simulation data (CSV or JSON) for external analysis*
- [ ] *README and user documentation*

---

## *Open Questions*
1. ***Radius vs. mass/energy scaling: RESOLVED.**** The geometric doubling (1, 2, 4, 8...) correctly describes the spin *path amplitudes* — the actual orbital radii. The 1820 ratio is an *energy* ratio, not a radius ratio. From "Unifying the Electron and Proton" (Mathis, 2008): a particle with axial spin only has energy 9 (1 rest + 8 spin). Adding x, y, z spins multiplies energy via [1 + (8 × 16)/2], [1 + (8 × 16 × 32)/2²], [1 + (8 × 16 × 32 × 64)/2⁴] = 9, 65, 1025, 16385. The ratio 16385/9 = 1820.56, matching the nucleon/electron mass ratio. The electron at rest is a proton stripped of its outer three spins (x, y, z at the baryon tier). **Implementation:** geometric amplitudes for path computation, energy formula for mass/energy readouts. The `SpinLevel` struct should carry both `amplitude` (geometric, for position) and `energy_contribution` (from this formula, for display).
2. **Collision penetration depth:** When a field photon enters the focus particle's outer tier, how deep does it go? Is this deterministic (based on size and spin comparison) or probabilistic? The papers describe charge channeling through nuclei as a definite mechanical process, suggesting deterministic.
3. **The 30° emission angle:** Is this an input to the model or an emergent output? The magnetic moment paper describes it as arising from the spin geometry. Our simulation should reproduce it without hard-coding it. If it doesn't, either our spin model or our collision model needs adjustment.
4. **Multi-particle compute budget:** When simulating multiple baryons, each with its own field interactions, the compute cost scales quadratically. Spatial partitioning (octree or spatial hash) will be needed. Design for this from the start even though it's an M5 concern.
5. **Photon spin level for field particles:** Field photons in the sim — how many spin levels do they carry? A 1-level photon is just a spinning sphere. A 4-level photon has full z-spin wavelength. The outer spin is what matters for most interactions, but inner spins affect penetration behavior. Default to 4-level (full photon tier) with outer-spin-only collision math, and optionally allow full spin resolution for detailed studies.