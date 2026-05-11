# M4: Field Simulation — Implementation Plan

## Overview

The charge field is a cloud of B-photons streaming through the simulation volume at c, colliding with the actual base solid underlying the focus particle (the fundamental spinning sphere whose trajectory the ghost spheres merely illustrate), transferring momentum, and being recycled when they leave the volume. Ghost spheres are visualization-only — they show the current orbital trajectory of each spin level but are not physical surfaces. The focus particle should drift, spin, and react to asymmetric bombardment — all emergent, nothing hard-coded.

Target: 500K photons at 60 FPS on an RX 550 4GB. Scale to 5M+ on better hardware.

---

## Architecture & Data Flow

```
                         ┌──────────────────────────────┐
   CPU (Rust)            │        GPU (Vulkan 1.1)       │
   ───────────           │                                │
   compose()             │  ┌──────────────────────┐      │
   base_solid_pos() ─────┼─>│ Uniform: FocusState   │      │
                         │  └──────────┬───────────┘      │
                         │             │                   │
                         │  ┌──────────▼───────────┐      │
                         │  │ Compute: field_update │      │
                         │  │ (advance + recycle)   │      │
                         │  └──────────┬───────────┘      │
                         │             │                   │
                         │  ┌──────────▼───────────┐      │
                         │  │ Compute: collide      │      │
                         │  │ (detect + deflect)    │      │
                         │  └──────────┬───────────┘      │
                         │             │                   │
                         │  ┌──────────▼───────────┐      │
                         │  │ Compute: impulse_sum  │      │
                         │  │ (parallel reduction)  │      │
                         │  └──────────┬───────────┘      │
                         │             │                   │
   apply_impulse() <─────┼─────────────┘ (4-float readback)│
                         │                                │
                         │  ┌──────────────────────┐      │
                         │  │ Render: point sprites │      │
   MultiMesh.set_buffer()│  │ (from CPU readback)   │      │
   ──────────────────────┼─>│                       │      │
                         │  └──────────────────────┘      │
                         └──────────────────────────────┘
```

**Key decision:** Use a local RenderingDevice for compute, CPU readback for rendering via MultiMesh. The readback is 500K × 12 bytes = 6MB ≈ 0.05ms on the RX 550's 112 GB/s bus. This avoids the complexity of sharing SSBOs between compute and render contexts. If readback becomes a bottleneck at higher photon counts, we can switch to `RenderingServer.get_rendering_device()` and share the buffer directly.

---

## Photon Data Structure

### GPU-side (SSBO, std430)

```glsl
// 32 bytes per photon, naturally aligned for vec4 access
struct Photon {
    vec4 pos_energy;   // xyz = position (natural units), w = energy
    vec4 vel_flags;    // xyz = velocity (unit direction × speed), w = packed flags
};
// flags (as uint via floatBitsToUint):
//   bit 0:     chirality (0 = CW/photon, 1 = CCW/antiphoton)
//   bit 1:     active (1 = alive, 0 = dead/recycling)
//   bits 2-3:  collision state (0 = none, 1 = deflected, 2 = channeled)
//   bits 4-7:  spin_level_count (0-4 for field photons)
//   bits 8-31: reserved
```

### Rust-side (CPU, for initialization and readback)

```rust
#[repr(C)]
pub struct FieldPhoton {
    pub position: [f32; 3],
    pub energy: f32,
    pub velocity: [f32; 3],
    pub flags: u32,
}
// 32 bytes, matches GPU layout exactly
```

### Memory budget

| Count | Buffer Size | RX 550 VRAM % |
|-------|-------------|---------------|
| 100K  | 3.2 MB      | 0.08%         |
| 500K  | 16 MB       | 0.4%          |
| 2M    | 64 MB       | 1.6%          |
| 5M    | 160 MB      | 4.0%          |
| 50M   | 1.6 GB      | 40%           |

Default: **500K**. Safe headroom at all tiers. The impulse and reduction buffers add negligible overhead (< 1MB total).

---

## Simulation Volume

A sphere centered on the focus particle's world position. Radius adapts to the particle's effective radius:

```
sim_radius = max(effective_radius * 3.0, 10.0)
```

The 3× factor gives photons room to approach from all angles. The floor of 10.0 ensures a reasonable volume at level 1 (effective_radius = 1).

### Boundary recycling

When a photon exits the sphere (distance from center > sim_radius):

1. Place it at a random point on the sphere surface
2. Set velocity to c (1.0) pointing inward:
   - **Isotropic field:** random inward direction (cosine-weighted toward center)
   - **Biased field:** inward direction skewed toward `field_direction`
3. Reset collision state flags
4. Preserve chirality (or re-roll based on `photon_ratio`)

This keeps the photon count constant with zero allocation. The photon "teleports" from the exit point to a new entry point. At steady state the volume maintains roughly uniform density.

### Initial distribution

On sim start, place all photons randomly within the volume:
- Position: uniform random within the sphere
- Velocity: random direction × c
- Chirality: CW with probability `photon_ratio`, CCW otherwise
- Energy: drawn from the preset's energy distribution

---

## Compute Shader Pipeline

Three dispatches per physics frame, in order. Workgroup size: 256 (standard for GCN).

### Pass 1: field_update.glsl

Advance positions, apply boundary conditions, recycle exited photons.

```glsl
#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 0) buffer PhotonBuffer {
    vec4 data[];  // 2 vec4s per photon
};

layout(std140, set = 0, binding = 1) uniform SimParams {
    vec4 sim_center_radius;   // xyz = center, w = radius
    float dt;
    float time_scale;
    float photon_ratio;       // fraction CW (0.0-1.0)
    float direction_strength; // 0 = isotropic, 1 = fully directional
    vec4 field_direction;     // xyz = direction, w = spawn_speed (should be 1.0 = c)
    uint photon_count;
    uint frame_number;        // for RNG seed
    uint _pad0;
    uint _pad1;
};

uint hash(uint x) {
    x ^= x >> 16; x *= 0x45d9f3bu;
    x ^= x >> 16; x *= 0x45d9f3bu;
    x ^= x >> 16; return x;
}

float rand01(inout uint seed) {
    seed = hash(seed);
    return float(seed) / 4294967295.0;
}

vec3 random_direction(inout uint seed) {
    float theta = rand01(seed) * 6.2831853;
    float z = 2.0 * rand01(seed) - 1.0;
    float r = sqrt(1.0 - z * z);
    return vec3(r * cos(theta), r * sin(theta), z);
}

void main() {
    uint idx = gl_GlobalInvocationID.x;
    if (idx >= photon_count) return;

    uint base = idx * 2;
    vec4 pos_e = data[base];
    vec4 vel_f = data[base + 1];

    vec3 pos = pos_e.xyz;
    vec3 vel = vel_f.xyz;

    // Advance position
    pos += vel * dt * time_scale;

    // Boundary check: recycle if outside sim sphere
    vec3 to_center = pos - sim_center_radius.xyz;
    float dist = length(to_center);

    if (dist > sim_center_radius.w) {
        // Respawn on sphere surface
        uint seed = hash(idx ^ (frame_number * 1999u));
        vec3 spawn_dir = random_direction(seed);
        pos = sim_center_radius.xyz + spawn_dir * sim_center_radius.w;

        // Velocity: inward with optional bias
        vec3 inward = normalize(sim_center_radius.xyz - pos);
        vec3 random_vel = random_direction(seed);
        // Blend between random-inward and field direction
        vec3 base_dir = mix(inward, field_direction.xyz, direction_strength);
        vel = normalize(mix(base_dir, random_vel, 0.3)) * field_direction.w;

        // Re-roll chirality
        uint flags = floatBitsToUint(vel_f.w);
        flags = (flags & ~1u) | (rand01(seed) > photon_ratio ? 1u : 0u);
        flags |= 2u; // set active
        vel_f.w = uintBitsToFloat(flags);
    }

    data[base] = vec4(pos, pos_e.w);
    data[base + 1] = vec4(vel, vel_f.w);
}
```

### Pass 2: collision_detect.glsl

Test each photon against the base solid — the actual physical sphere underlying the focus particle. Ghost spheres are visualization-only and play no role in collisions. The base solid's position is its composed world-space location (the result of all stacked spin transforms applied to the origin), and its radius is always 1.0 in natural units (the fundamental building block).

**IMPORTANT:** The collision target is NOT the ghost sphere hierarchy. Ghost spheres merely illustrate orbital trajectories — they are not physical surfaces. The base solid is the real particle that B-photons collide with. Its position changes every physics frame as the spin stack composes.

```glsl
#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 0) buffer PhotonBuffer {
    vec4 data[];
};

layout(std140, set = 0, binding = 1) uniform SimParams {
    // ... same as above ...
};

layout(std140, set = 0, binding = 2) uniform FocusState {
    vec4 base_solid;   // xyz = world-space position of the base solid, w = radius (1.0)
    // Future: may add velocity for relativistic corrections
};

layout(std430, set = 0, binding = 3) buffer ImpulseBuffer {
    vec4 impulses[];   // per-photon impulse (xyz = direction, w = magnitude)
};

void main() {
    uint idx = gl_GlobalInvocationID.x;
    if (idx >= photon_count) return;

    uint base = idx * 2;
    vec3 pos = data[base].xyz;
    vec3 vel = data[base + 1].xyz;
    float energy = data[base].w;
    vec4 impulse = vec4(0.0);

    // Test against the base solid (the real physical sphere)
    vec3 center = base_solid.xyz;
    float radius = base_solid.w;
    vec3 to_photon = pos - center;
    float dist = length(to_photon);

    if (dist < radius) {
        // Collision: reflect velocity about sphere normal
        vec3 normal = normalize(to_photon);
        float vn = dot(vel, normal);

        if (vn < 0.0) {
            // Photon is moving inward — deflect it
            vel = vel - 2.0 * vn * normal;

            // Push photon just outside to prevent re-collision
            pos = center + normal * (radius + 0.01);

            // Impulse transferred = change in momentum (2 × normal component)
            impulse = vec4(-2.0 * vn * normal * energy, abs(2.0 * vn * energy));
        }
    }

    // Write back
    data[base] = vec4(pos, energy);
    data[base + 1] = vec4(vel, data[base + 1].w);
    impulses[idx] = impulse;
}
```

### Pass 3: impulse_sum.glsl

Parallel reduction of per-photon impulses into a single net impulse vector. Two-stage:
1. Workgroup-level reduction via shared memory
2. Cross-workgroup reduction via a second dispatch (or atomics)

```glsl
#version 450

layout(local_size_x = 256) in;

layout(std430, set = 0, binding = 3) buffer ImpulseBuffer {
    vec4 impulses[];
};

layout(std430, set = 0, binding = 4) buffer ReductionOutput {
    vec4 partial_sums[];  // one per workgroup
};

shared vec4 sdata[256];

uniform uint element_count;

void main() {
    uint tid = gl_LocalInvocationID.x;
    uint gid = gl_GlobalInvocationID.x;

    sdata[tid] = (gid < element_count) ? impulses[gid] : vec4(0.0);
    barrier();

    for (uint s = 128u; s > 0u; s >>= 1u) {
        if (tid < s) sdata[tid] += sdata[tid + s];
        barrier();
    }

    if (tid == 0) {
        partial_sums[gl_WorkGroupID.x] = sdata[0];
    }
}
```

The CPU reads back `partial_sums` (at most ceil(500K/256) = 1954 vec4s = 31KB) and sums them. This avoids a second GPU dispatch and the 31KB readback is negligible.

---

## Rendering

### Point Sprites via MultiMesh

Use a MultiMesh with a tiny billboard quad (or IcoSphere at very low poly). Each photon gets a position and color.

```
MultiMeshInstance3D
  └─ MultiMesh
       ├─ mesh: QuadMesh (0.02 × 0.02) or SphereMesh (r=0.01, 4 segments)
       ├─ instance_count: 500000
       ├─ use_colors: true
       └─ transform_format: TRANSFORM_3D
```

**Update loop (per frame):**
1. GPU compute finishes → `rd.buffer_get_data(photon_buffer)` → PackedByteArray (16MB)
2. Rust function `build_multimesh_buffer(raw_data, photon_count)`:
   - Parses the raw photon data (32 bytes/photon)
   - For each photon, writes 12 floats (identity basis + position) + 1 float (packed color)
   - Returns PackedFloat32Array
3. GDScript: `multimesh.set_buffer(packed_array)` — single bulk upload

**Color encoding:**
- CW photon (chirality 0): blue `Color(0.3, 0.5, 1.0, 0.8)`
- CCW antiphoton (chirality 1): red `Color(1.0, 0.3, 0.3, 0.8)`
- Collision flash: white for 1 frame on collision

**Material:**
```
StandardMaterial3D:
  shading_mode: UNSHADED
  vertex_color_use_as_albedo: true
  transparency: ALPHA
  billboard_mode: BILLBOARD_ENABLED
  disable_receive_shadows: true
```

**Brightness/size:**
Fixed screen size (billboard). Doesn't scale with distance — photons convey density, not individual geometry. At 500K photons in a sphere of radius ~30 world units, the visual density is roughly 1 photon per 0.2 cubic world units. Visible but not overwhelming.

---

## CPU-Side Management (Rust)

### field_sim.rs

```rust
pub struct FieldSim {
    photons: Vec<FieldPhoton>,     // CPU mirror (for init + readback parsing)
    count: usize,
    sim_radius: f64,
    photon_ratio: f64,             // fraction CW
    field_direction: DVec3,
    direction_strength: f64,
    net_impulse: DVec3,            // latest frame's net impulse
    collision_count: u32,          // latest frame's collision count
    running: bool,
}
```

**Exposed to GDScript via FocusParticle or a separate FieldSim node:**
- `init_field(count: i32)` — allocate and randomize photons
- `get_init_buffer() -> PackedByteArray` — raw bytes for GPU upload
- `parse_readback(data: PackedByteArray)` — extract positions for MultiMesh
- `build_multimesh_buffer() -> PackedFloat32Array` — ready for MultiMesh.set_buffer()
- `sum_partial_impulses(data: PackedByteArray) -> Vector3` — CPU-side final reduction
- `get_net_impulse() -> Vector3` — latest net impulse for readout
- `get_collision_count() -> i32` — latest collision count
- `pack_base_solid() -> PackedFloat32Array` — pack base solid position + radius for GPU uniform (single vec4)
- Field parameter setters: ratio, direction, strength, sim_radius

### Impulse application

Each frame, after readback:
1. `sum_partial_impulses()` computes net impulse from the reduction buffer
2. GDScript passes net impulse to `FocusParticle.apply_field_impulse(impulse)`
3. FocusParticle adds impulse to its velocity/position (if field-driven motion is enabled)

The impulse is averaged over a configurable window (default 60 frames) for stability.

---

## GDScript Integration

### field_sim.gd (new Node, child of Main)

Responsibilities:
- RenderingDevice lifecycle (init, dispatch, cleanup)
- Buffer creation and management
- Per-frame dispatch loop
- Readback and handoff to Rust
- MultiMesh update

```
@export var focus_particle_path: NodePath
@export var default_photon_count: int = 500000
@export var auto_start: bool = false

var _rd: RenderingDevice
var _update_shader: RID
var _collide_shader: RID
var _reduce_shader: RID
var _update_pipeline: RID
var _collide_pipeline: RID
var _reduce_pipeline: RID
var _photon_buffer: RID
var _params_buffer: RID
var _focus_state_buffer: RID
var _impulse_buffer: RID
var _reduction_buffer: RID
var _update_uniform_set: RID
var _focus_uniform_set: RID
var _reduce_uniform_set: RID
var _running: bool = false
```

### field_renderer.gd (new Node, child of Main)

Responsibilities:
- MultiMeshInstance3D management
- Reading photon positions from Rust
- Color assignment by chirality
- Visibility toggle

---

## UI Controls

### Field Panel (Bottom of screen, per PROJECT_DESIGN)

```
┌─────────────────────────────────────────────────────┐
│ Field: [Start] [Stop] [Reset]   Preset: [Earth ▼]  │
│                                                      │
│ Density: ────●──── 500K    VRAM: 16MB / 4096MB      │
│ CW/CCW:  ────●──── 67%     Collisions: 1,204/s     │
│ Direction: ──●──── 0.3     Net impulse: 0.042       │
│ Dir axis: [+Y ▼]           Impulse dir: ↗ (+X,+Y)  │
└─────────────────────────────────────────────────────┘
```

### Field Presets

| Preset      | Ratio | Dir Strength | Dir Axis | Energy Dist |
|-------------|-------|-------------|----------|-------------|
| Earth       | 0.67  | 0.3         | -Y       | thermal     |
| Solar Wind  | 0.75  | 0.7         | +X       | high-energy |
| Vacuum      | 0.50  | 0.0         | —        | uniform-low |
| Custom      | user  | user        | user     | user        |

**Earth:** Slight downward bias (gravity-like), 2:1 CW/CCW ratio.
**Solar Wind:** Stronger directional bias, higher energy, more CW-dominated.
**Vacuum:** Perfectly isotropic, equal CW/CCW, low energy (cosmic background).

---

## Performance Budget (RX 550)

**Hardware:** 8 CUs × 64 SPs = 512 shader processors. 112 GB/s memory bandwidth. 1.2 TFLOPS.

### Per-frame at 500K photons

| Stage              | Data moved   | Estimated time |
|--------------------|-------------|----------------|
| field_update       | 32MB R + 32MB W | 0.57ms     |
| collision_detect   | 32MB R + 32MB W + 8MB W (impulse) | 0.64ms |
| impulse_sum        | 8MB R + 31KB W | 0.07ms      |
| CPU readback (pos) | 6MB          | 0.05ms         |
| MultiMesh upload   | 26MB         | 0.23ms         |
| **Total GPU compute** | | **~1.3ms** |
| **Total CPU overhead** | | **~0.3ms** |

Well within the 16.6ms frame budget. At 2M photons, multiply by 4: ~5.2ms GPU, ~1.2ms CPU. Still comfortable.

### VRAM at 500K

| Buffer           | Size   |
|------------------|--------|
| Photon SSBO      | 16 MB  |
| Impulse SSBO     | 8 MB   |
| Reduction SSBO   | 31 KB  |
| Params uniform   | 128 B  |
| FocusState uniform | 16 B |
| MultiMesh buffer | 26 MB  |
| **Total**        | **~50 MB** |

4% of the 4GB budget. Plenty of room.

---

## Implementation Phases

### Phase 1: Visible Field (Sessions 1-2)

**Goal:** Dots on screen, moving through a sphere, recycling at the boundary.

1. `field_sim.rs`: FieldPhoton struct, random init, `get_init_buffer()`, `parse_readback()`, `build_multimesh_buffer()`
2. `field_update.glsl`: Position advance + boundary recycling (no collisions yet)
3. `field_sim.gd`: RenderingDevice setup, buffer creation, single compute dispatch, readback loop
4. `field_renderer.gd`: MultiMesh setup, per-frame buffer update
5. UI: Start/Stop/Reset buttons, photon count slider, field visibility toggle (hotkey F)
6. Verify: dots stream through the volume at c, recycle at boundary, hold 60 FPS

**Checkpoint:** A cloud of blue/red dots flowing through the space around the spinning focus particle.

### Phase 2: Collisions (Sessions 3-4)

**Goal:** Photons bounce off the base solid, net impulse displayed.

1. `collision_detect.glsl`: Base solid collision test, velocity reflection, impulse write
2. `impulse_sum.glsl`: Parallel reduction + CPU final sum
3. Upload base solid position + radius each frame (pack from Rust via composed world position)
4. `FocusParticle.apply_field_impulse()`: accumulate impulse, apply as velocity change
5. Visual: collision flash (white for 1 frame), or just let the deflection be visible
6. UI: collision rate counter, net impulse magnitude + direction arrow
7. Test: symmetric field → no net motion. Asymmetric → drift.

**Checkpoint:** Visible deflections at the base solid, particle drifts in asymmetric fields.

### Phase 3: Field Presets & Chirality (Sessions 5-6)

**Goal:** Full field parameter control, chirality-dependent interactions.

1. Implement preset system (Earth, Solar Wind, Vacuum, Custom)
2. CW/CCW ratio slider with visual feedback (blue/red balance changes)
3. Direction bias with axis selector
4. Chirality-aware collision: same-chirality = hard deflection, opposite = partial penetration
5. Collision statistics: rate, angle distribution, impulse histogram
6. VRAM usage display (actual vs. allocated)
7. Impulse history graph (rolling 60-frame window)

**Checkpoint:** Earth preset shows slight downward drift. Solar wind pushes the particle. Vacuum is at rest.

### Phase 4: Verification & Polish (Session 7)

**Goal:** Sanity-check against Mathis predictions.

1. Does the proton emit charge at ~30° from the poles? (Measure collision angle distribution)
2. Does the neutron block more charge than the proton? (Compare collision rates by chirality pattern)
3. Does an asymmetric field produce linear drift proportional to the asymmetry?
4. Performance at 2M and 5M photons — does the pipeline scale?
5. UI polish: tooltips, keyboard shortcuts, field stats panel

---

## Open Questions for Implementation

1. **Compute dispatch timing:** Should the compute dispatch happen in `_physics_process` (synchronized with the spin stack advance) or `_process` (every render frame)? Physics sync is more correct but might miss frames. Recommend: dispatch in `_physics_process`, readback in `_process` (1-frame latency is imperceptible).

2. **Base solid position packing:** The collision shader needs the base solid's world-space position, updated every frame. The base solid position is the composed result of all spin transforms (i.e., where the actual particle currently sits in space). We need a `pack_base_solid() -> PackedFloat32Array` accessor on FocusParticle that returns [cx,cy,cz,r] — a single vec4 (position + radius). Ghost spheres are visualization-only and are NOT uploaded to the collision shader.

3. **MultiMesh vs. custom SSBO rendering:** MultiMesh with CPU readback works for 500K but may bottleneck at 5M+ (26MB upload per frame × 5 = 130MB). At that scale, switch to a custom vertex shader that reads directly from the compute SSBO. This requires using the global RenderingDevice instead of a local one. Design the abstraction layer so the switch is painless.

4. **Photon spin levels:** PROJECT_DESIGN specifies field photons can have 1-4 spin levels. For M4 MVP, all field photons are single-level (radius 1, speed c). Adding multi-level field photons is an M5/M6 refinement that affects penetration depth.

5. **Collision response sophistication:** MVP uses simple specular reflection off the base solid. The full Mathis model has charge channeling (photon enters one pole, exits the other), partial absorption, and spin-dependent deflection angles. The collision target is always the base solid — sophistication means varying the deflection *response*, not the collision *target*. Build collision_detect.glsl with a `collision_mode` uniform so we can swap between simple and sophisticated models without rewriting the shader.

6. **Impulse window averaging:** Raw per-frame impulse is noisy. Averaging over 60 frames (1 second at 60 FPS) smooths the signal. The averaging window should be configurable. Implement as a ring buffer in Rust (like PathTrace but for impulse vectors).

7. **Sim radius adaptation:** When the user adds spin levels, the effective radius jumps (e.g., 8 → 64 when activating tier 2). The sim radius should grow accordingly, but existing photons are now inside the "new" particle. Options: (a) recycle all photons on level change, (b) let them collide and deflect naturally. Option (b) is more physical and produces a visible "shockwave" when a new tier activates — could be a cool effect.

---

## File Inventory

```
NEW FILES:
  rust/src/field_sim.rs                    — FieldPhoton, buffer packing, impulse math
  godot/scripts/field_sim.gd               — GPU compute dispatch, RenderingDevice management
  godot/scripts/field_renderer.gd          — MultiMesh point sprite rendering
  godot/shaders/compute/field_update.glsl  — Position advance + boundary recycling
  godot/shaders/compute/collision_detect.glsl — Ghost sphere collision + deflection
  godot/shaders/compute/impulse_sum.glsl   — Parallel impulse reduction

MODIFIED FILES:
  rust/src/lib.rs                          — Register FieldSim class
  rust/src/focus_particle.rs               — Add pack_base_solid(), apply_field_impulse()
  godot/scripts/hud.gd                     — Field panel UI, preset selector, stats display
  godot/scenes/main.tscn                   — Add FieldSim + FieldRenderer nodes
```
