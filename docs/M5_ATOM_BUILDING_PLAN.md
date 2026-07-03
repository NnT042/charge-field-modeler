# M5: Atom-Building Mode — Design Document

## Overview

M5 introduces a new simulation mode: **atom-building**. Instead of visualizing a single particle's internal spin stack (Phase 1), we simulate multiple particles interacting via their **charge emission profiles**. Each particle is abstracted to a point mass with an angular emission histogram — the data exported from Phase 1. Forces between particles emerge from the overlap of these profiles with the unified field (gravity + charge). No bonding rules are prescribed; H₂ formation, electron capture, and nuclear structure should arise from the physics.

Target: start with free hydrogen (proton + electron), demonstrate spontaneous H₂ bond formation, then scale to heavier elements via nuclear presets.

---

## Physics Foundation (from Mathis research)

### Two-Field Unified Model

Every interaction is a vector sum of two fields:

| Field | Direction | Cause | Radial falloff |
|-------|-----------|-------|----------------|
| Solo gravity | Inward (apparent attraction) | Matter expansion | 1/r² |
| Charge (E/M) | Outward (always repulsive) | Photon bombardment | 1/r⁴ |

Gravity dominates at distance (slower falloff). Charge dominates up close (stronger per unit area). The crossover creates stable structures: orbits, bonds, nuclear arrangements.

At the proton-electron scale (from Mathis, `quantumg.html`):
- Gravity toward electron: 9.81 m/s²
- Charge repulsion at Bohr distance: 0.00064 m/s²
- Net: 9.80 m/s² apparent attraction

### Electron Position

Electrons do **not** orbit the nucleus at the Bohr radius. Per Mathis (`diatom.pdf`, `fine4.pdf`):

- The Bohr radius (corrected: 9 × 10⁻⁹ m, 177× the standard value) is the **capture radius** — the limit of the proton's charge vortex.
- Electrons are channeled in the charge field like photons, but are too large (too wide a spin radius) to be recycled through the proton.
- They orbit the **pole of a specific proton**, right at the nuclear boundary — like a ping pong ball circling a drain.
- When captured at one pole, the opposite pole's vortex mostly closes. One electron per proton, one pole per electron.
- Four possible states per hydrogen atom: electron at (top/bottom) × (left-spin/right-spin).

**Why the electron doesn't fall through the pole channel:**
Phase 1 histogram data confirms hypothesis (2) from session 22 — the electron's equatorial emission disc (~76° wide at 25% intensity) is significantly wider than the proton's polar intake tunnel (~56° cone). The electron's charge pressure footprint, not its physical body, determines the contact boundary. All three particle types share nearly identical polar tunnel widths (~25-28° from pole), since the innermost spin levels define the axial null zone. The particles differ in what they emit *outside* that tunnel.

### Charge Channeling: Protons vs Neutrons

| Property | Proton | Neutron |
|----------|--------|---------|
| Charge path | Pole → equator (2D, planar) | Pole → pole (1D, linear) |
| Analogy | Spinning fan / lawn sprinkler | Lightning rod |
| Relative channeling strength | 1.0 | 0.685 (from magnetic moment) |
| Orientation in nucleus | Equator faces inward | Pole faces inward |
| Through-charge | Minimal | Primary mode |

Protons create equatorial emission (the magnetic vector). Neutrons create through-charge (the electrical vector). This distinction drives nuclear architecture.

### H₂ Bond Mechanics

Per Mathis (`diatom.pdf`), the diatomic hydrogen bond is **polar, not equatorial**:

1. Each hydrogen atom has an electron at one pole, creating charge asymmetry.
2. The pole without the electron has a charge minimum (low recycling).
3. When two atoms align with electrons on the **outside** (away from each other), a charge minimum forms between them.
4. The ambient charge field pushes them together from above and below — a push, not a pull.
5. If electrons are **between** the protons, competing vortices drive atoms apart.

Bonding rule: 4 of 8 possible spin/pole combinations create bonds. The discriminant is electron position (outside vs between). This should emerge from the force model without being prescribed.

### Nuclear Diagrams (for presets)

Mathis models nuclei as plug-and-socket arrangements of disks:
- **Alphas** (He-4): two protons hole-to-hole + two neutrons. Fundamental nuclear building block.
- **Carousel level**: equatorial ring of protons flung outward by angular momentum.
- **Axial level**: pole-stacked protons/neutrons channeling through-charge.
- Bonds are protons plugging into charge holes. Structure defines bonding, not electrons.

Heavier elements are built from alphas, with protons and neutrons plugging into pole and carousel positions. These structures are **fused** — forced together under extreme pressure (stellar furnace), then stabilized by charge channeling. They will be represented as rigid presets, not dynamically simulated at the nuclear scale.

---

## Measured Emission Profiles (Phase 1 Data)

All profiles collected at 4M+ hits, 1° bin resolution, from Phase 1 field simulation (session 22, 2026-05-24). Source files in `reference/`.

### Proton

Source: `reference/histogram_proton_2026-05-24T22-39-24.csv`

- **Shape**: flying saucer — strong equatorial disc with bimodal shoulder peaks
- **Absolute peak**: +1° at 2.36%
- **Central dip**: 0° at 1.98% (equatorial plane is transition zone, not maximum)
- **Bimodal bump**: +29° to +34° at 0.98-1.09% (secondary peak, consistent with earlier 7°/58° bimodal finding from proton pole frame)
- **Valley**: ±20-25° drops to 0.60-0.74%
- **Polar minimum**: < 0.05% beyond ±62° latitude
- **Polar tunnel**: ±28° half-angle from pole axis (emission < 0.05)
- **Chirality asymmetry**: sharp cliff at +3°→+4° (1.68→1.37), more gradual on negative side

### Neutron

Source: `reference/histogram_neutron_2026-05-24T22-44-03.csv`

- **Shape**: squashed tin can — broad flat emission with sharp edge cutoffs
- **Peak**: ±57-58° at 1.21-1.23% (near the cutoff edge, NOT at equator)
- **Equatorial value**: 0° at 0.53% — the *minimum* of the active emission band
- **Plateau**: ±55° range, steady 0.53-1.18% (gradually slopes up from equator to edges)
- **Edge cliff**: ±60°, drops from ~1.0 to 0.35 in one bin
- **Polar leakage**: ~0.002% at ±80° — 4× the proton's polar emission (through-charge signature)
- **Polar tunnel**: ±27° half-angle from pole axis

The neutron profile slopes upward from equator to ±57°, opposite of the proton. Equatorial emission minimum + enhanced polar leakage = pole-to-pole channeling confirmed.

### Electron

Source: `reference/histogram_electron_(at_rest)_2026-05-24T22-30-33.csv`

- **Shape**: fat frisbee with 30° echo
- **Absolute peak**: -3° at 2.24%
- **Central dip**: 0° at 1.71% (same pattern as proton)
- **Core disc**: ±5° at >2.0 (90% of peak)
- **Strong disc**: ±17° at >1.0 (45% of peak)
- **30° bump**: -29° at 1.11% — scale-model echo of proton's bimodal structure (same stacked-spin geometry, fewer levels, much less charge recycling)
- **±9° minimum**: sharp asymmetric drop, +8° cliff from 1.41→0.97 (chirality signature)
- **Polar tunnel**: ±25° half-angle from pole axis
- **Disc extent at 25% of peak**: ±42° from equator (76° full width)

### Profile Assumptions and Caveats

These histograms have parameters baked in that are subject to future tuning:
- Outer spin angular velocity (the dominant free parameter)
- Collision model details from Phase 1 (tier boundaries, deflection rules)
- Chirality preset for the spin stack

The qualitative shapes (proton disc + bimodal bump, neutron tin can, electron razor disc) should be robust to parameter changes. Exact peak positions and widths may shift.

---

## Data Structures

### ParticleProfile

Defines a particle type. Loaded from JSON, referenced by ID at runtime.

```rust
struct ParticleProfile {
    name: String,                         // "proton", "neutron", "electron", "alpha", ...
    profile_type: ProfileType,            // Free or NuclearPreset
    mass: f64,                            // natural units (proton = 1.0)
    radius: f64,                          // natural units (proton = 1.0)
    emission: EmissionTable,              // intensity vs polar angle
    absorption: EmissionTable,            // how it receives charge (derived or hand-tuned)
    total_emission: f64,                  // integral over sphere (for normalization)
    constituents: Option<Vec<Constituent>>, // for nuclear presets only
}

enum ProfileType {
    Free,            // single baryon/electron, fully dynamic
    NuclearPreset,   // locked multi-baryon structure
}
```

### EmissionTable

The histogram, query-optimized. Axially symmetric → 1D table indexed by polar angle θ.

```rust
struct EmissionTable {
    bins: Vec<f32>,       // intensity values, normalized so max = 1.0
    bin_count: usize,     // 180 = 1° resolution (validated as sufficient by Phase 1 data)
    symmetry: Symmetry,
}

enum Symmetry {
    Bilateral,  // north ≈ south — table covers 0°..90°, mirrored (proton, electron, neutron)
    Polar,      // north ≠ south — table covers 0°..180° (some nuclear presets)
}

impl EmissionTable {
    fn sample(&self, cos_theta: f32) -> f32 {
        let theta = cos_theta.acos();
        let t = theta / PI;
        let idx = (t * self.bin_count as f32) as usize;
        self.bins[idx.min(self.bin_count - 1)]
    }
}
```

### Constituent (for nuclear presets)

```rust
struct Constituent {
    profile_id: usize,        // index into profile registry
    local_position: Vec3,     // offset from preset center
    local_orientation: Quat,  // rotation relative to preset frame
}
```

Composite emission profile is precomputed offline by superimposing all constituent emissions in their locked orientations.

### SimParticle (runtime instance)

```rust
struct SimParticle {
    profile_id: usize,        // which ParticleProfile
    position: Vec3,
    velocity: Vec3,
    orientation: Quat,         // pole axis direction
    angular_velocity: Vec3,    // spin rate and axis
    mass: f64,                 // copied from profile (or overridden)
}
```

### Profile Data Format (on disk)

JSON, consistent with existing `particle_presets.json`:

```json
{
    "name": "proton",
    "type": "free",
    "mass": 1.0,
    "radius": 1.0,
    "emission": {
        "symmetry": "bilateral",
        "bins": [0.05, 0.08, 0.12, "...90 values, pole to equator..."]
    },
    "absorption": {
        "symmetry": "bilateral",
        "bins": ["...derived or hand-tuned..."]
    }
}
```

Nuclear presets add `constituents`:

```json
{
    "name": "alpha",
    "type": "nuclear_preset",
    "mass": 4.0,
    "radius": 2.0,
    "emission": { "...composite, precomputed..." },
    "constituents": [
        { "profile": "proton", "position": [0, 0, 0.5], "orientation": [0, 0, 0, 1] },
        { "profile": "proton", "position": [0, 0, -0.5], "orientation": [1, 0, 0, 0] },
        { "profile": "neutron", "position": [0.5, 0, 0], "orientation": [0, 0, 0, 1] },
        { "profile": "neutron", "position": [-0.5, 0, 0], "orientation": [0, 0, 0, 1] }
    ]
}
```

---

## Force Model

### Pairwise Force: A acting on B

```
Given:
  d     = unit vector from A to B
  r     = distance between A and B
  θ_A   = angle from A's pole axis to d        (what part of A faces B)
  θ_B   = angle from B's pole axis to -d       (what part of B receives)

Force on B from A:
  F_gravity =  G_q * m_A * m_B / r²   * d      (toward A — attractive)
  F_charge  = -C_q * m_A * E_A(θ_A) * R_B(θ_B) / r⁴  * d  (away from A — repulsive)
  F_net     =  F_gravity + F_charge
```

Where:
- `E_A(θ_A)` = A's emission intensity at angle θ_A (table lookup)
- `R_B(θ_B)` = B's absorption/reception intensity at angle θ_B (table lookup)
- `G_q`, `C_q` = coupling constants (calibrated so that known behaviors emerge)
- θ is computed from `dot(pole_axis, d)` — two dot products per pair

### Absorption Profile

Charge enters at poles, exits at equator. First approximation:

```
absorption(θ) = 1.0 - normalized_emission(θ)
```

Where emission is high (equator), absorption is low. Where emission is low (poles), absorption is high. This gives:
- Proton aimed pole-first at another proton's equator: max bombardment × min resistance = max repulsion
- Two protons pole-to-pole: min bombardment × max absorption = min repulsion → gravity wins → alpha formation

**Exception: neutron.** Its pole-to-pole channeling means incoming polar charge passes *through* rather than being absorbed. The neutron needs either a hand-tuned absorption table or a `channeling_factor` that bleeds incoming pole charge to the opposite pole.

### Torque

Charge hitting B off-center of its mass creates torque, aligning particles:

```
torque_on_B = offset_B × F_charge
```

Where `offset_B` is the displacement from B's center to the effective point of charge pressure, derived from the asymmetry of B's absorption profile at the current angle. Equatorial bombardment twists the particle until its pole faces the charge source (minimizing force → stable orientation).

### Ambient Field

Constant background, representing Earth-surface conditions:

- **Gravity**: uniform downward, 9.81 m/s² (scaled to sim units)
- **Charge pressure (up)**: ~0.00955 m/s² (Mathis, `atmo.html`) — this levitates the atmosphere
- **Isotropic charge background**: uniform pressure from all directions (the cosmic charge field). This is what creates bonds — it fills the space around particles and pushes them together wherever a charge shadow (minimum) exists between them.

The ambient isotropic pressure is the critical ingredient for H₂ bonding. Without it, two hydrogen atoms with electrons on the outside would have a charge minimum between them but nothing to push them together. The ambient field provides that push.

### Integration: Velocity Verlet

Symplectic integrator — conserves energy over long runs, no drift, dead simple:

```
v(t + dt/2) = v(t) + a(t) * dt/2           // half-step velocity
x(t + dt)   = x(t) + v(t + dt/2) * dt      // full-step position
a(t + dt)   = F(x(t + dt)) / m             // recompute forces at new position
v(t + dt)   = v(t + dt/2) + a(t + dt) * dt/2  // complete velocity step
```

Same structure for angular dynamics (orientation as quaternion, angular velocity as vec3, torque as angular force).

**Why Verlet over RK4:** RK4 is 4th-order accurate but not symplectic — it slowly gains or loses energy over long runs. Verlet is only 2nd-order but exactly conserves the Hamiltonian, which matters more for watching bonds form and hold over thousands of ticks. We're not doing high-eccentricity orbits where RK4's accuracy advantage pays off.

**Timestep:** adaptive — shrink when particles are close (strong forces, fast dynamics), grow when they're far apart (weak forces, slow drift). Floor at `dt_min` to prevent tunneling through each other.

---

## Calibration

### Coupling Constants

Two free parameters: `G_q` (gravity coupling) and `C_q` (charge coupling). These set the force scale.

**Known constraints:**
- At r = 100,000 × proton radius (Bohr distance): charge is negligible vs gravity. Net force ≈ gravity only.
- At r = 1 × proton radius (contact): charge dominates enormously.
- The crossover distance (where F_gravity = F_charge for a proton-proton pair) should be somewhere in the range of 10-1000× proton radius.

**Calibration procedure:**
1. Place two protons at varying distances, measure equilibrium separation.
2. Adjust `G_q / C_q` ratio until equilibrium matches the expected scale.
3. Add an electron to one proton's pole — it should settle at approximately proton-surface distance, not fly away or collapse.
4. Bring two hydrogen atoms together — they should form H₂ when electrons are on the outside and repel when electrons are between.

### Ambient Field Strength

The isotropic background charge pressure determines how strongly bonds are pushed together. Too weak: particles find charge minima but drift. Too strong: everything collapses. Calibrate by:
- H₂ bond: the bond should hold but vibrate (thermal motion). If it's rigid, pressure is too high. If it drifts apart, too low.
- Electron capture: an electron approaching a proton's pole should spiral in and settle, not overshoot and escape.

---

## Architecture

### Mode Switch

M5 is a new mode within the existing program, not a separate binary. The main scene gets a mode toggle:
- **Phase 1 mode**: single focus particle with spin stack + field sim (existing)
- **Phase 2 mode**: atom-building — multi-particle, profile-based forces, no internal spin rendering

Mode switch clears the 3D viewport and loads the appropriate scene tree. Rust backend detects mode and runs the appropriate simulation loop.

### Simulation Loop (Phase 2)

```
each tick:
    for each pair (A, B):
        compute F_net(A→B) and F_net(B→A)
        compute torque on A from B, and on B from A
    for each particle:
        add ambient field forces
        integrate (Verlet)
        update orientation from angular velocity
    render
```

For N particles, pairwise force is O(N²). Acceptable for small counts (< 100 particles). For larger scenes, spatial partitioning (octree or grid) with cutoff distance.

### Rendering (Phase 2)

Each particle is rendered as a **mesh generated from its emission profile**:
- Revolve the 1D emission histogram around the pole axis → 3D surface of revolution
- Proton: flying saucer shape
- Neutron: squashed tin can
- Electron: thin disc/frisbee
- Alpha preset: composite shape from constituent superposition

Color-code by type. Opacity or glow intensity encodes charge emission strength. Orientation matches the particle's quaternion.

Electrons attached to proton poles rendered at contact distance, visually small, orbiting the pole.

---

## Particle Types at Launch

### Free Particles

| Type | Mass | Radius | Profile source | Notes |
|------|------|--------|---------------|-------|
| Proton | 1.0 | 1.0 | histogram_proton CSV | Flying saucer |
| Neutron | 1.0 | 1.0 | histogram_neutron CSV | Tin can, 0.685 channeling factor |
| Electron | 1/1836 | 1/1836 | histogram_electron CSV | Frisbee, nearly massless |

### Nuclear Presets (future, post-H₂)

| Type | Constituents | Notes |
|------|-------------|-------|
| Alpha (He-4) | 2p + 2n, hole-to-hole | Fundamental nuclear building block |
| Carbon-12 | 3 alphas + 4 carousel protons | First carousel-level element |
| Nitrogen-14 | Carbon + south proton + north neutron | Ammonia target |
| Oxygen-16 | 4 alphas + pole protons | Water target |

Presets are rigid: internal structure is locked, only the composite emission profile interacts with the field. Fusion (building presets dynamically) is a stretch goal — "star furnace mode."

---

## Implementation Phases

### Phase 2a: Foundation

- [x] Profile data format and loader (CSV → EmissionTable, bilateral symmetry)
- [x] Export Phase 1 histograms to profile data (proton, neutron, electron CSVs in godot/config/)
- [x] SimParticle struct with position, velocity, orientation, angular_velocity
- [x] Pairwise force computation (gravity 1/r² + charge 1/r⁴ with profile lookup)
- [x] Velocity Verlet integrator (linear + angular, with softening and damping)
- [x] Ambient field (uniform gravity + isotropic charge pressure + shadow-based pairwise attraction)
- [x] Basic 3D rendering (MultiMesh spheres + pole axis indicators, colored by type)
- [x] Torque model (sin(2θ) restoring torque, aligns equator toward charge source)
- [x] Mode switching (Tab key toggles Phase 1 ↔ Phase 2)
- [x] On-screen diagnostics (particle positions, velocities, distances, KE, coupling constants)

### Phase 2b: Hydrogen — COMPLETE (session 28)

- [x] Spawn proton + electron, verify electron capture at pole
      (`hydrogen_capture` test; wall-riding orbit at r=1.3, θ≈11°,
      v_tan = corotation speed — diatom.pdf's "circling the drain" AT the
      nuclear boundary)
- [x] Spawn two hydrogen atoms, verify H₂ bond formation
      (polar bond, vibrating around d≈4.7 — `h2_bond_matrix`)
- [x] Test all 8 spin/pole combinations — confirm 4 bond, 4 repel
      (emerges from stream-cushion + stoppered-vortex forces, no coded rule)
- [x] Calibrate G_q / C_q ratio for stable bonds
      (constants LOCKED with Mathis-sourced derivations in
      `Couplings::default()`, rust/src/atom_core.rs; pinned by
      `derived_constants_equilibrium`)
- [x] Calibrate ambient field for realistic bond vibration
      (DEVIATION: the bond standoff comes from the stream-collision cushion
      vs gravity+intake, not from P_amb — the pairwise shadow term stays
      available at 0 for environment effects)
- [x] Profile-based mesh rendering (surface of revolution from histogram)

### Phase 2c: Presets and Heavier Elements — COMPLETE except molecules (session 28)

- [x] Nuclear preset data format and loader
      (DEVIATION: Rust const data in `atom_core::preset_constituents()`,
      not JSON — presets must be spawnable headlessly by tests, and GDScript
      only ever spawns by name)
- [x] Composite profile baking — REJECTED in favor of rigid groups:
      constituents stay real particles (per-constituent nearfield lets an
      electron capture at one specific proton's pole) integrated as one
      rigid body. A baked profile has no constituent-level nearfield.
- [x] Alpha particle preset (2p short stack + 2n posts, oxygen.pdf geometry;
      `alpha_holds_and_conserves`, `alpha_captures_electron`)
- [x] Helium atom: alpha + 2 captured electrons (`helium_stable`)
- [x] Carbon, Nitrogen, Oxygen presets (3-alpha stack + polar plugs;
      `heavier_presets_smoke`)
- [ ] Simple molecule tests (H₂O, NH₃)

### Phase 2d: Tools and Polish

Session 28 delivered the polished HUD (scenes/ui/atom_hud.tscn: toolbar,
scenario buttons, readout panel, debug-hidden tuning sliders, mode-switch
buttons in both modes) and the headless test loop (docs/M5_TEST_LOOP.md).
Remaining 2d items:

- [ ] Particle spawning UI (click to place, set type, set orientation)
- [ ] Real-time force vector visualization (arrows showing F_gravity, F_charge, F_net)
- [ ] Energy readout (kinetic + potential, verify conservation)
- [ ] Bond detection and display (highlight when particles are in bound state)
- [ ] Export simulation state (positions, velocities, forces as CSV/JSON)
- [ ] Star furnace mode (stretch goal: crank ambient pressure, watch fusion)

---

## Open Questions

1. **Neutron absorption model:** Simple complement (1 - emission) doesn't capture pole-to-pole through-channeling. May need a hand-tuned absorption table or a `channeling_factor` that redirects incoming polar charge to the exit pole, reducing absorption and creating the apparent neutrality.

2. **Electron spin dynamics:** Does the electron's orientation (pole axis direction) matter for force calculation at this abstraction level? In Phase 1 it was a full spin stack. In Phase 2, do we track its orientation or treat it as a freely-spinning disc whose time-averaged profile is symmetric?

3. **Scale rendering:** Electrons are 1/1836 the mass and ~1/1836 the radius of protons. At proton-scale zoom they're invisible. At electron-scale zoom the proton fills the screen. Need a rendering strategy — possibly exaggerated electron size with a different visual treatment (glow, halo) rather than physical scale.

4. **Timestep tuning:** The force ratio between contact distance (charge-dominated) and Bohr distance (gravity-dominated) spans ~10²¹. Adaptive timestep needs careful floor/ceiling to avoid both tunneling and stalling.

5. **Isotropic pressure implementation:** The ambient charge field pushes from all directions. Computationally, this is equivalent to a restoring force toward charge minima — but implementing it as an actual omnidirectional field vs. an effective potential changes the torque behavior. Need to decide which approach gives more physical results.

---

## References (Mathis papers, via search_mathis)

- `diatom.pdf` — H₂ bond mechanics, electron position, 4/8 bonding rule
- `fine4.pdf` — Bohr radius as capture limit, electron at nuclear boundary
- `quantumg.html` — Quantum gravity, unified field force at proton-electron scale
- `magneton.html` — Corrected Bohr radius (177×), electron radius corrections
- `elec3.html` — Electron radius = 1/c²
- `neutron.pdf` — Neutron charge channeling, pole-to-pole throughput
- `stack.html` — Nuclear construction without strong force, alpha formation
- `per4.pdf` — Period 4 nuclear diagrams, Iron/Copper carousel structure
- `ethane.pdf` — Carbon nuclear diagram, molecular plug-and-socket bonds
- `phos.pdf` — Neutron as linear channeler vs proton as planar channeler
- `gauss.pdf` / `gauss2.pdf` — Unified field equations, gravity-charge equivalence
- `disp2.pdf` — Charge field as foundational, E = C/g
- `atmo.html` — Earth charge field strength (0.00955 m/s²)
- `water2.pdf` — Water/ice bond mechanics, charge channeling through Oxygen
- `ammon.pdf` — Ammonia bond angles from nuclear structure + neutron position
