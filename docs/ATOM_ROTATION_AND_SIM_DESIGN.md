# Atom Mode — Alpha Sub-Unit Refactor & Sim-Driven Nuclei

**Status:** DESIGN — ready for implementation.
**Author:** Opus (session 30, 2026-07-07), for Fable to review + implement.
**Scope:** (Part 1) fix the carousel rotation bug by giving each alpha its
own rotational frame; (Part 2) make nuclear presets force-driven, staged
from rigid-alpha to full-nucleon.

This doc is self-contained. Read `docs/PROJECT_DESIGN.md` and
`docs/PHYSICS_REFERENCE.md` for background, and the memory notes
`project-start-here` / `project-session29-learnings` for the visual-rotation
history — but the diagnosis below was re-verified fresh against the current
build and supersedes earlier guesses.

---

## 0. Verified diagnosis (do not re-litigate)

The rotation bug was reproduced empirically against the loaded build. Method
(re-runnable): spawn `neon`, log a carousel alpha's members, call
`advance_display(π/2 / CAROUSEL_VIS_RATE)`, log again. Result over a 90°
carousel orbit:

| Quantity | before | after 90° | correct? |
|---|---|---|---|
| Proton **pole** (donut-hole axis) | (0.81, 0, 0.59) | (0.59, 0, −0.81) | ✅ tracks the orbit |
| Disc **face pattern** ("seam" = local X) | (0, −1, 0) | (0, −1, 0) | ❌ world-locked, never spins |
| Neutron post **vertical offset** in the piece | +0.70 / −0.70 | +0.70 / −0.70 | ❌ pinned top/bottom, never revolves |

### Root cause

There is exactly **one** rotation channel for group members: `car_rot =
DQuat::from_rotation_y(carousel_phase)`, a spin about the **nucleus** Y axis
(`sync_group_members`, atom_core.rs:983, 989, 998). It is forced to do two
physically distinct jobs at once:

1. **Orbit** each non-core alpha around the nucleus axis (the carousel), and
2. **Roll** each alpha's neutron posts around **that alpha's own** axis.

For the **core** alpha these two axes coincide (the alpha's stack axis *is*
the nucleus Y), so `car_rot` does both correctly — this is why the core
"looks right." For any **carousel / sideways / cap** alpha, its own axis is
radial or edge-on, **not** Y, so `car_rot` orbits the piece fine (poles do
track) but cannot roll its posts about the alpha's own axis — they stay
welded to the top/bottom of the piece. Separately, the current WIP froze
`display_spin_phase` for carousel members (atom_core.rs:1040) to kill an
earlier "twirl," which is why their disc faces are world-locked.

**The load-bearing missing thing:** a second rotational DOF — *each alpha
spins about its own axis, independent of the carousel orbit*. Every prior
rewrite re-plumbed the single `car_rot` channel, which is why fixing one
symptom always reintroduced another (twirl ↔ frozen ↔ wobble). One phase
variable cannot cover two independent rotations.

---

## Part 1 — Alpha sub-unit refactor (the presentable fix)

Chosen approach (user, session 30): **refactor to alpha sub-units now** — do
not band-aid the flat member list. The alpha-frame structure introduced here
is also the substrate Part 2 needs, so there is no throwaway work.

### 1.1 The idea in one paragraph

Model a nucleus as a set of **alpha units**. Each alpha rotates **rigidly
about its own axis** by a `roll_phase` (this both spins the two proton discs
in place and revolves the two neutron posts around the axis — a rigid spin of
an off-axis disc does both). Non-core alphas **additionally** have their
center + axis revolve around the nucleus axis by the shared `carousel_phase`.
Two independent phases, advanced in wall clock. The core-vs-carousel
asymmetry disappears by construction: every alpha rolls about its own axis;
only the orbiting ones also carousel.

### 1.2 Data model

Introduce an `AlphaUnit`. Convention: **the alpha's own stack axis is its
local +Y** (same convention as a free particle's pole).

```rust
pub struct AlphaUnit {
    /// Nucleon particle ids in this alpha (2 protons + 2 posts for a full
    /// alpha; 1 for a lone plug; etc.).
    pub members: Vec<usize>,
    /// Each member's rest position in the ALPHA's own frame (+Y = axis):
    ///   proton  → (0, ±NUCLEON_PITCH/2, 0)
    ///   post    → (POST_R, 0, 0) and (−POST_R, 0, 0)  [revolve via roll]
    pub member_local_pos: Vec<DVec3>,
    /// Each member's pole in the alpha frame. Discs whose axis is the alpha
    /// axis use +Y; a plug proton pointing at its partner uses its own dir.
    pub member_local_pole: Vec<DVec3>,
    /// Rest axis of this alpha in the NUCLEUS frame (Y core, radial u
    /// carousel, X sideways/plug-edge-on).
    pub rest_axis: DVec3,
    /// Rest center of this alpha in the NUCLEUS frame (0 core, u·CAROUSEL_R
    /// carousel, (0,y,0) axial cap/connector, (0,y,±gap) paired plug).
    pub rest_center: DVec3,
    /// Does this alpha's center+axis orbit the nucleus axis (carousel)?
    pub orbits_core: bool,
    /// Kinematic roll about the alpha's own axis (rad, wall clock).
    pub roll_phase: f64,
    pub roll_rate: f64,

    // ── Part 2 promotes the alpha to a real body (unused in RigidLock) ──
    pub com: DVec3,
    pub velocity: DVec3,
    pub orientation: DQuat,      // maps alpha +Y → world axis
    pub angular_velocity: DVec3,
    pub mass: f64,
    pub inertia: f64,
}
```

`RigidGroup` keeps its nucleus-level fields (`com`, `velocity`,
`orientation`, `angular_velocity`, `mass`, `inertia`, `carousel_phase`,
`carousel_rate`, `skin_reach`, `disc_exits`) and **replaces** the flat
`members / local_offsets / local_orients / member_spin / member_spin_phase /
carousel` arrays with `pub alphas: Vec<AlphaUnit>`. Keep a flat
`members()` accessor (flatten `alphas[*].members`) so existing call sites
(rendering buffer, occlusion, flow VFX) keep working.

### 1.3 Kinematics (exact — this is the whole fix)

For a member `m` in alpha `a`, in **RigidLock** (kinematic) mode:

```
R_a   = orientation_from_pole(a.rest_axis)      // alpha +Y → rest_axis
Roll  = DQuat::from_rotation_y(a.roll_phase)     // about alpha-local +Y
Car   = if a.orbits_core { DQuat::from_rotation_y(carousel_phase) }
        else { DQuat::IDENTITY }

// position in the nucleus frame:
inner       = R_a * (Roll * member_local_pos[m])     // roll revolves posts
nucleus_pos = Car * (a.rest_center + inner)           // carousel orbits piece

// orientation in the nucleus frame:
q_nucleus   = Car * R_a * Roll * orientation_from_pole(member_local_pole[m])

// nucleus → world:
p.position    = com + nucleus_orientation * nucleus_pos
p.orientation = (nucleus_orientation * q_nucleus).normalize()
```

Key properties (all verify directly with a kinematic test):
- **Protons** sit on the alpha axis → `Roll` leaves them in place but spins
  their faces. ✅ visible disc spin.
- **Posts** sit off-axis → `Roll` revolves them around the alpha axis,
  staying 180° apart. ✅ fixes the pinned neutrons on *every* alpha.
- **Non-core alphas** get `Car` on center+axis → the piece orbits the core
  and its poles track radially. ✅ (already worked; still works.)
- **Core alpha** has `orbits_core = false`, so `Car = I`; its posts revolve
  purely from its own `roll_phase`. ✅ decoupled from the carousel.

`render_orientation` (atom_core.rs:3320) for group members becomes simply
`p.orientation` — the roll already lives in `p.orientation`. **Do NOT** add a
separate `rotY(display_spin_phase)` on top, and **do NOT** reintroduce
swing–twist extraction (it is singular near 180° configs — session-29). The
per-member `display_spin_phase` field and its freeze logic
(atom_core.rs:1040, 866) are **deleted**; disc spin now comes from
`roll_phase`.

### 1.4 advance_display

```rust
pub fn advance_display(&mut self, delta: f64) {
    // free particles: unchanged (parallel-transported display_orientation)
    for gi in 0..self.groups.len() {
        let g = &mut self.groups[gi];
        g.carousel_phase = (g.carousel_phase + g.carousel_rate * delta) % TAU;
        for a in &mut g.alphas {
            a.roll_phase = (a.roll_phase + a.roll_rate * delta) % TAU;
        }
        self.sync_group_members(gi);
    }
}
```

Rates: `roll_rate = DISPLAY_SPIN_RATE` (0.9), `carousel_rate =
CAROUSEL_VIS_RATE`. **Revert the WIP** `CAROUSEL_VIS_RATE = DISPLAY_SPIN_RATE`
back to an independent value (was 0.35) — with roll and orbit now separate
DOFs there is no reason to lock them equal, and equal rates make the whole
thing look like one solid gear. Randomize **initial `roll_phase` per alpha**
and `carousel_phase` per nucleus so pieces don't turn in lockstep
("synchronized ballet"). All members *within* one alpha share that alpha's
roll (they are one rigid alpha — coherence is correct here).

### 1.5 Rebuilding the presets in alpha-unit terms

Rewrite `preset_constituents` (atom_core.rs:378) as a builder that returns
`Vec<AlphaUnit>` rather than a flat `Vec<Constituent>`. Each helper below
returns one `AlphaUnit` (or a small list). Geometry constants keep their
current values (NUCLEON_PITCH 2.6, ALPHA_PITCH 3.75, PLUG_PAIR_GAP 0.8;
introduce `POST_R = 0.7`). **Reset `CAROUSEL_R` from the WIP 8.0 back toward
~4.5** — the 8.0 was a workaround for skin-circle overlap, not a geometry
decision; solve overlap on the skin side (§1.7), not by flinging the alphas
out where they look detached.

| Alpha helper | rest_axis | rest_center | orbits_core | notes |
|---|---|---|---|---|
| core `alpha(y)` | +Y | (0, y, 0) | false | protons ±pitch/2 on axis, posts ±POST_R |
| `carousel_alpha(φ)` | u=(cosφ,0,sinφ) | u·CAROUSEL_R | **true** | radial stack, posts revolve about u |
| `sideways_alpha(y)` | +X | (0, y, 0) | **true** | edge-on connector; axis sweeps with carousel |
| `cap alpha(y)` | +Y | (0, y, 0) | false | argon caps, parallel to core |
| `plug_proton(y,z)` | edge-on / toward partner | (0, y, z) | **true** | 1 member; paired plug revolves about pole hole |
| `plug_neutron(y,z)` | +Y (pole-on axis) | (0, y, z) | **true** | 1 member; graphene.pdf |

The physics rationale for each piece is unchanged from the current comments
(oxygen.pdf posts, nuclear.pdf carousel, graphene.pdf plug neutron, etc.) —
carry those comments over. Preset compositions (`carbon` = 3 core alphas,
`neon` = core + 4 carousel, `argon` = caps + connectors + center + 4
carousel, `nitrogen`/`oxygen` plugs) are identical; only the grouping
changes.

### 1.6 What else touches the flat member arrays (update these)

- `sync_group_members` (:962) — rewritten per §1.3.
- `compute_group_reach` / `compute_disc_exits` / `compute_skin_segments`
  (:892–896, :907) — currently classify pieces by lateral offset (`>1.5 →
  carousel`). Re-key them off `AlphaUnit` membership directly (an alpha *is*
  a skin segment) — cleaner and removes the magic 1.5.
- Flow VFX path builders (`flow_group_through` / `flow_group_disc`) and
  `build_group_carousel_overlay` — consume `alphas` instead of the lateral
  test. `get_group_carousel_phase` (used by atom_mode.gd:565 overlay
  rotation) stays.
- `spawn_preset` (:823) — build `alphas`, then `sync_group_members`.
- Inertia/mass aggregation — sum over `alphas[*].members`.

### 1.7 Skin overlap (the real reason CAROUSEL_R got bumped)

Neighboring carousel alphas' skin bands overlapped at CAROUSEL_R 4.5, so the
WIP pushed them to 8.0. Instead: clip each alpha's skin band to an azimuthal
wedge (±45° for 4-fold carousel) or shrink the band radius to
`min(field_reach, 0.9·half_spacing)` so adjacent bands kiss but don't
interpenetrate. Keep CAROUSEL_R at the physically-motivated ~4.5
(nearest carousel proton pole just outside the center disc edge, edge-to-hole
— four.pdf).

### 1.8 Acceptance tests (headless, assert don't eyeball)

Add to `mod tests`:
- `alpha_roll_revolves_posts`: spawn any preset; advance_display a quarter
  roll period; a core alpha's post offset (in the nucleus XZ plane) rotates
  ~90°, protons stay on axis.
- `carousel_orbit_tracks_pole`: neon carousel proton pole rotates with the
  orbit AND its posts revolve about the *alpha* axis (the exact case that
  fails today — see §0 table).
- `core_and_carousel_symmetric`: a core alpha and a carousel alpha, given
  equal roll rates, both revolve their posts about their own axes (no
  special-casing).
- `no_lockstep`: two alphas in one nucleus have different roll phases after
  spawn.
- Keep both headless smokes green (scene + `dev_smoke_atom.gd`).

---

## Part 2 — Sim-driven nuclei (staged)

Chosen approach (user, session 30): **a mixture.** Alphas stay stuck
together for now but **push on each other like the H₂ free protons, only
stronger**; add an **alternate mode** that frees the individual nucleons for
experimentation until a stable balance of settings is found. So: a
**granularity switch**, three levels, sharing the AlphaUnit structure from
Part 1.

```rust
pub enum NucleusDynamics {
    RigidLock,    // current + Part 1 kinematics. Presentable default.
    RigidAlpha,   // alphas are rigid bodies; forces act BETWEEN alphas.
    FreeNucleon,  // every nucleon free; forces alone hold alphas together.
}
```

### 2.1 The one rule that changes: the intra-group skip

Today `compute_forces` skips **all** same-group pairs (atom_core.rs:1528):

```rust
if a[i].group.is_some() && a[i].group == a[j].group { continue; }
```

Replace with a mode-aware predicate:

| Mode | Skip pair (i,j) when… | Forces that act |
|---|---|---|
| RigidLock | same group (any) | none internal (current) |
| RigidAlpha | same **alpha** | inter-alpha (× `INTRA_NUCLEUS_BOOST`) + external |
| FreeNucleon | never | all pairs (× boost) + external |

`INTRA_NUCLEUS_BOOST` (start ~2–4×, tunable) is the "stronger than H₂"
coupling that binds alphas tighter than a loose molecule. It multiplies the
same pairwise terms already used for free protons (charge E(θ)R(θ)/r⁴·doppler,
intake, stream cushion, corotation — the LOCKED force table in
`project-m5-plan`). No new force law; the nuclear binding is the *same*
physics at a stronger coupling, which is the whole point — if it can't hold a
nucleus together with the molecular terms scaled up, the force model is wrong
and we want to know.

### 2.2 Integration per mode

- **RigidLock:** unchanged (§1). Kinematic phases drive motion.
- **RigidAlpha:** each `AlphaUnit` integrates as a rigid body (port the
  existing group `kick`/`drift` aggregation at atom_core.rs:1297–1358 down
  one level — aggregate its members' forces into an alpha COM force + torque,
  integrate `com`/`orientation`, then place members rigidly from the alpha
  frame). The nucleus `RigidGroup` is **no longer a rigid lock** in this
  mode — it is a *container* of alphas coupled only by forces. `carousel_phase`
  / `roll_phase` become **initial conditions**; real ω takes over. Keep
  `GROUP_SPIN_RELAX`-style damping available per-alpha but start it low — you
  *want* the carousel to be able to turn here.
- **FreeNucleon:** no rigid constraint; every nucleon integrates freely with
  the normal `kick`/`drift`. Intra-alpha cohesion must come from the force
  model; if it won't hold, add an explicit short-range attractive bond
  (a tunable spring toward the alpha rest geometry) as a scaffold while
  tuning, then try to remove it.

### 2.3 Validation harness (do this before wiring visuals)

Build all of it in `atom_scenarios.rs` first (fast `cargo test --release`
loop), mirroring how the M5 force model was locked:
- **`alpha_stays_bound`** (RigidAlpha): spawn one alpha's worth of alphas at
  rest; over N steps, inter-alpha spacing stays within ±X% (doesn't fly
  apart or collapse). Sweep `INTRA_NUCLEUS_BOOST`.
- **`carousel_self_organizes`** (RigidAlpha): spawn neon with the carousel
  alphas at rest (no scripted phase); the equatorial charge output should
  spin the ring up to a steady carousel rate. Assert net carousel ω settles
  nonzero and bounded.
- **`neon_is_inert`**: a probe proton approaching neon's pole is not
  captured (the six-sided closure — nuclear.pdf).
- **`nucleon_balance`** (FreeNucleon): report the settling — is there a
  stable configuration, and at what settings? This is the sandbox the user
  asked for; the deliverable is the tuned constant table, not a pass/fail.
- Only after a mode passes its harness: expose a HUD toggle
  (RigidLock ⇄ RigidAlpha ⇄ FreeNucleon) so the user can watch it.

### 2.4 Phase-2 endgame (why this matters)

Once RigidAlpha holds nuclei together *by force* and reproduces the carousel
emergently, the same machinery lets you place two nuclei + electrons in a
box, integrate, and watch them bond/repel — i.e. the periodic-table /
chemistry-experiment goal. FreeNucleon is the research mode for confirming
the alpha itself is a force equilibrium, not an imposed constraint.

---

## Constants & sync points

- New: `POST_R = 0.7`, `INTRA_NUCLEUS_BOOST` (~2–4, tune), `NucleusDynamics`.
- Revert WIP: `CAROUSEL_R` 8.0 → ~4.5 (fix skin overlap in the skin builder,
  §1.7); `CAROUSEL_VIS_RATE` back to an independent ~0.35 (not
  `= DISPLAY_SPIN_RATE`).
- Delete: per-member `display_spin_phase` + its freeze branch
  (atom_core.rs:866, 1040); `member_spin` / `member_spin_phase` arrays.
- GDScript stays in sync: `advance_display(_delta)` call (atom_mode.gd:401),
  `get_group_carousel_phase` overlay (atom_mode.gd:565), any preset spawn
  args. The h2 scenario RIDE_AXIAL/RIDE_LATERAL/COROT_V literals
  (project-start-here gotcha) are unaffected.

## Build / test loop (unchanged, and one trap)

- `cargo test --release --manifest-path rust/Cargo.toml` (~10s gate; never
  `cd` into rust/; always **release** build — the `.gdextension` points at
  `target/release/`).
- Then two headless smokes: the scene, and
  `-s res://scripts/dev_smoke_atom.gd`.
- **Build trap:** a running Godot locks `charge_field_modeler.dll`; the
  linker then silently writes `~charge_field_modeler.dll` and Godot keeps
  running the OLD binary (this is exactly what happened this session — a
  14:19 rebuild couldn't replace the 17:14 DLL). **Close Godot before
  `cargo build --release`**, or confirm the DLL mtime advanced.

## Suggested implementation order

1. `AlphaUnit` struct + `preset_constituents → Vec<AlphaUnit>` builders (§1.2, §1.5).
2. `sync_group_members` + `advance_display` + `render_orientation` (§1.3–1.4).
3. Update skin/flow/overlay call sites to consume `alphas` (§1.6).
4. Part-1 acceptance tests (§1.8) — must pass, plus both smokes. **Ship the
   presentable fix here.**
5. `NucleusDynamics` enum + mode-aware skip predicate (§2.1).
6. RigidAlpha integration + `atom_scenarios` harness (§2.2–2.3).
7. FreeNucleon mode + tuning sandbox (§2.2–2.3).
8. HUD mode toggle + skin overlap cleanup + CAROUSEL_R restore (§1.7).

---

## Session-31 implementation addendum (Fable — BINDING, refines the above)

Reviewed against the live code (78 tests green at baseline, WIP diff in
place). These decisions refine §1.2/§1.6/§1.8 where the doc and the code
disagree or where the doc would break a working invariant. Where this
addendum conflicts with the sections above, the addendum wins.

### A1. Keep the flat rest arrays as DERIVED data (refines §1.2, §1.6)

`RigidGroup` KEEPS `members: Vec<usize>`, `local_offsets: Vec<DVec3>`,
`local_orients: Vec<DQuat>`, `member_spin: Vec<f64>` — but they become
**rest-pose data derived from `alphas` at spawn** (roll = 0, carousel = 0
canonical pose; `local_offsets[k] = a.rest_center + R_a * member_local_pos`,
`local_orients[k] = R_a * orientation_from_pole(member_local_pole)`).
They are never mutated after spawn. Rationale: `segment_sources`,
`march_reach`, `compute_disc_exits`, `build_group_carousel_overlay`, and the
flow-VFX emitter loop all consume rest geometry (azimuth-averaged or drawn in
a phase-rotated child node) — keeping the rest arrays means ZERO changes to
those consumers' math for Part 1.

DELETE: `member_spin_phase`, `carousel: Vec<bool>` (RigidGroup),
`Constituent` + the old builders, `SimParticle::display_spin_phase`, and the
WIP freeze branch in `advance_display`. `AlphaUnit::members` holds **k
indices into the flat arrays** (not particle ids); particle id = 
`group.members[k]`.

### A2. Alpha builders — exact rest frames (refines §1.5)

All alphas use the local convention: axis = +Y, protons at
`(0, ∓NUCLEON_PITCH/2, 0)` pole `+Y`, posts at `(±POST_R, 0, 0)` pole `+Y`
(a post's pole IS the alpha axis, matching today's geometry in every
builder). `POST_R = 0.7` (new const, replaces the four literal 0.7s).

| Builder | rest_axis | rest_center | orbits_core |
|---|---|---|---|
| `core_alpha(y)` (core, caps) | `+Y` | `(0,y,0)` | false |
| `carousel_alpha(φ)` | `u=(cosφ,0,sinφ)` | `u·CAROUSEL_R` | true |
| `sideways_alpha(y)` | `+X` | `(0,y,0)` | true |
| `plug_proton(y,z)` | `(0,0,−sgn z)` if z≠0 else `+X` | `(0,y,z)` | true |
| `plug_neutron(y,z)` | `+Y` | `(0,y,z)` | true |

Plugs are single-member alphas: `member_local_pos = [(0,0,0)]`,
`member_local_pole = [+Y]` — the old world-frame pole becomes the
rest_axis, so `Roll` spins the plug disc about its own pole and `Car`
rides it around the socket exactly as the old `carousel: true` flag did.
Carry over every physics comment (oxygen.pdf posts, graphene.pdf plug
neutron, nuclear.pdf carousel, etc.) onto the new builders.

### A3. Member velocities in sync_group_members (completes §1.3)

Keep the rigid field `v = g.velocity + ω_group × world_off` and the member
axial spin `ω_member = ω_group + pole·member_spin[k]` as today, then add the
kinematic display terms (parity with the old car_omega handling):

```
axis_w   = nucleus_orientation * (Car * a.rest_axis)        // alpha axis, world
center_w = com + nucleus_orientation * (Car * a.rest_center) // alpha center, world
if a.orbits_core:
    car_w = (nucleus_orientation * DVec3::Y) * g.carousel_rate
    v += car_w × world_off;      ω_member += car_w
v += (axis_w * a.roll_rate) × (p.position − center_w)
ω_member += axis_w * a.roll_rate
```

### A4. Skin segmentation re-key (refines §1.6 — carbon must stay ONE tube)

§1.6's "an alpha *is* a skin segment" is WRONG for the axial stack: the
`skin_segments_match_structure` test (and the user-approved v3 skin) needs
carbon to read as ONE tube, not three. Correct re-key:

- An alpha belongs to the **carousel level** iff its `rest_center` lateral
  distance `> 1.5` (compares 0.8 (plug pair) vs 4.5 (CAROUSEL_R) — no longer
  a per-member magic number). All carousel-level alphas' members form the
  single carousel segment, as today.
- All other alphas' members go through the existing y-sorted `SEG_GAP`
  clustering UNCHANGED (that keeps carbon whole and splits argon's
  center/connector/cap pieces at the 1.75 gaps).

### A5. Disc exits re-key (refines §1.6 — kills the pairing hack)

Replace the sorted-consecutive-proton-pairing in `compute_disc_exits` with:
one disc center per alpha with `|rest_axis·Y| > 0.9` AND ≥ 2 proton members,
at `y = rest_center.y`. Same output as today for every preset (carbon 3,
neon 1, argon 3), no fragile pairing. The θ-crossing solve against
`skin_reach` stays as is.

### A6. Overlay re-key + §1.7 overlap fix

`build_group_carousel_overlay`: iterate carousel-level alphas directly (one
circle per alpha, unit = its members) instead of azimuth-clustering the
member offsets. Clamp each circle radius:
`r_disp = min(marched r_disp, 0.45 · CAROUSEL_R · √2)` (adjacent carousel
centers sit `CAROUSEL_R·√2` apart on a 4-fold ring — circles kiss, never
interpenetrate, §1.7). Restore the carousel piece's wide-disc skin band in
`build_group_skin_mesh` (delete the WIP `if *carousel { continue; }`) — the
band is axisymmetric and does not overlap per-unit; only the overlay circles
needed the clamp.

### A7. Rigidity invariant SPLIT (refines §1.8 — heavier_presets_smoke)

With independent per-alpha roll, cross-alpha member distances legitimately
change under `advance_display` (post azimuths decorrelate). The old blanket
assertion would fail by design. New invariant in `heavier_presets_smoke`:

- **Physics rigidity:** all-pairs distances over the first 8 members
  unchanged under `step_n(N)` (phases don't advance in `step`).
- **Display rigidity:** WITHIN one alpha, all pair distances unchanged under
  `advance_display(0.8)` (an alpha is rigid; check the first full alpha).

The plug-geometry assertions in that test rewrite in `AlphaUnit` terms:
nitrogen/oxygen plug alphas are single-member, `orbits_core == true`, plug
proton `rest_axis·Y ≈ 0`, plug neutron `rest_axis ≈ +Y`; oxygen pairs sit at
equal `rest_center.y`, `|Δrest_center| > 0.5`, proton rest_axis pointing at
the paired neutron.

### A8. Test bookkeeping

- REWRITE `display_phases_are_wall_clock`: `step()` advances no roll/carousel
  phase; `advance_display(dt)` advances `carousel_phase` by
  `CAROUSEL_VIS_RATE·dt` and every `roll_phase` by `DISPLAY_SPIN_RATE·dt`.
- DELETE the WIP carousel-freeze assertions and
  `carousel_circles_tangent_at_new_radius`; ADD an overlay-clamp test
  (every overlay circle radius ≤ 0.45·CAROUSEL_R·√2 + ε at R = 4.5).
- ADD §1.8's four tests. For `carousel_orbit_tracks_pole` use the §0 method:
  neon, quarter carousel orbit via advance_display with roll_rate zeroed on
  the probed alpha (set it directly in the test) so orbit and roll effects
  are separable; then a second advance with roll to see the posts revolve
  about the alpha axis.
- `heavier_presets_smoke` per A7. Everything else must pass UNCHANGED —
  buffer formats, GDScript API (`get_group_carousel_phase`,
  `advance_display`, overlay packing) are frozen for Part 1.

### A9. Part 2 — where the mode lives (refines §2.1)

`NucleusDynamics` is a GLOBAL setting on `AtomCore`
(`pub dynamics: NucleusDynamics`, default `RigidLock`), exposed through
`AtomSim` as `set_nucleus_dynamics(mode: i32)` / `get_nucleus_dynamics()
-> i32` (0/1/2). Per-group mixing is not needed for the sandbox.
`SimParticle` gains `pub alpha: Option<usize>` (index into its group's
`alphas`), set at spawn — the skip predicate needs it:

```
RigidLock:   skip iff same group
RigidAlpha:  skip iff same group AND same alpha index
FreeNucleon: never skip
```

Non-skipped same-group pairs get every pairwise force term multiplied by
`INTRA_NUCLEUS_BOOST` (new pub const, start 3.0, sweep in the harness).

### A10. Part 2 — mode transitions (refines §2.2)

`set_nucleus_dynamics` seeds state on entry:
- **→ RigidAlpha / FreeNucleon:** capture each alpha's CURRENT kinematic
  POSE as its body state: `orientation = g.orientation·Car·R_a·Roll`,
  `com = g.com + g.orientation·(Car·rest_center)`. Velocities seed from
  the group's REAL rigid field only (`velocity` = group v + group ω ×
  offset, `angular_velocity` = group ω) — do NOT add the display roll or
  carousel rates: they are render-channel readability rates, and seeding
  them as physical momentum injects fictional energy whose ~2 s decay
  drags the force equilibrium through a phase-dependent migration that
  can shed an alpha (the session-31 round-5 dissolution). FreeNucleon
  additionally leaves each member particle with its current
  position/velocity/ω (already correct from the last sync) and simply
  stops rigid syncing.
- **→ RigidLock:** snap back to the kinematic pose (resume phase-driven
  sync from the stored `roll_phase`/`carousel_phase`). Positions jump;
  that's fine for a sandbox toggle — document it in the HUD tooltip.
- In RigidAlpha, `kick`/`drift` integrate ALPHAS (aggregate member force →
  alpha COM force + torque about the alpha com; scalar inertia
  `Σ m(|local|² + 0.4 r²)` about the alpha center) and members are placed
  rigidly from the alpha body frame each drift. The nucleus-level
  `RigidGroup` com/orientation stop integrating (container only). Keep a
  weak per-alpha damping `ALPHA_SPIN_RELAX = 0.5` (1/s) — well below the
  nucleus lock's 20.0, the carousel must stay able to turn (§2.2).
- In FreeNucleon, members integrate as free particles (normal kick/drift
  paths); `sync_group_members` and the group branch of `advance_display`
  are gated OFF for non-RigidLock modes. Composite skins/flow keep
  anchoring to the (now frozen) group frame — acceptable for the sandbox,
  noted as a known visual limitation.

### A11. Part 2 — harness honesty (refines §2.3)

The harness must not pretend the physics passes before it's tuned:
- `alpha_stays_bound` (RigidAlpha, one nucleus's alphas): HARD test — but
  it sweeps `INTRA_NUCLEUS_BOOST` over {1, 2, 3, 4, 6} internally and
  asserts SOME value keeps max inter-alpha drift within ±30% over the run;
  it prints the per-boost table so the winner is visible in test output.
- `neon_is_inert` (RigidAlpha): HARD test — a probe proton approaching
  neon's pole from 12 r does not come within capture range (< 2.0) in N
  steps.
- `carousel_self_organizes`, `nucleon_balance`: `#[ignore]`d report
  scenarios (run with `cargo test --release -- --ignored <name>
  --nocapture`); they print settling metrics (net carousel ω, spacing
  drift, KE) instead of asserting — the deliverable is the tuned-constant
  table, and hard-asserting emergence before tuning would just be a red
  suite. Promote to hard tests once the user signs off on constants.

### A13. Binding v2 — channeling + ambient, NOT a boost (session 31 round 3)

The A11 harness proved a uniform multiplier cannot bind alphas (it scales
attraction and repulsion equally; alphas ejected at boost 1–12,
monotonically worse). The Mathis-faithful redesign, agreed with the user:

**Sources.** bb2.pdf: "attraction must always be explained as loss of
repulsion" (spin cancellations lower the between-field's repulsive
energy). strong.html: "charge is channeled through the nucleus by baryon
spin, and so does not cause a repulsion between protons... There is no
charge field within the nucleus." nuclear.pdf: "the charge field is both
the initial pressure and the subsequent glue."

**Three mechanisms, all same-group-gated (molecular force table frozen):**

1. **Channeling attenuation.** For non-skipped same-group pairs
   (RigidAlpha/FreeNucleon), plugged neighbors route charge through each
   other instead of colliding:
   `C = channeling · max(cos²θ_i, cos²θ_j) · c_dist(r)` where
   `cosθ = pole·d̂` (a partner presenting its HOLE channels; two facing
   equators — the molecular repel config — get C≈0), and `c_dist` is a
   smoothstep from 1 at r ≤ NUCLEON_PITCH (2.6, inside the funnel mouth)
   to 0 at r ≥ 2·NUCLEON_PITCH. Apply `c_q × (1−C)` and `stream × (1−C)`.
   g_q, intake, torque, vortex, corot stay full — intake IS the
   channeled flow (diamag.pdf: polar protons are "fans, pulling charge
   in"), and it is what holds the edge-to-hole carousel plugs.
2. **Nuclear ambient pressure.** Same-group pairs use a new coupling
   `nuclear_ambient` in place of the (default-0) molecular
   `ambient_pressure` in the EXISTING shadow term
   `P·(1−E_i)(1−E_j)·(1−occ)/r²` — the ambient field pushing fused
   bodies into each other's charge shadows: the "subsequent glue".
3. **Disc-aware contact.** Facing protons of adjacent alphas rest at
   1.15 while sphere-contact fires at r_i+r_j = 2.0 — an unquenchable
   ~85-unit repulsion that no channeling can touch (contact is never
   attenuated, correctly). But the bodies ARE discs (the entire model is
   planar emission); a pole-on approach nests into the funnel and only
   the thin disc waist can collide: effective contact radius per side
   `r_eff = r · (POLE_HALF_THICKNESS + (1−POLE_HALF_THICKNESS)·sinθ)`
   with `POLE_HALF_THICKNESS = 0.35`, θ = angle of d̂ from that body's
   pole. Equator-on contact unchanged (1.0); pole-on stacks contact at
   0.35·r. GLOBAL (not group-gated) — a disc is a disc — but if any
   molecular scenario test moves, fall back to same-group-gated and
   report which test forced it.

**New Couplings fields:** `channeling: f64` (default 1.0 = full effect
strength, sweep dimension) and `nuclear_ambient: f64` (sweep; expect
O(1–20)). `intra_nucleus_boost` default DROPS to 1.0 — the uniform-boost
experiment is concluded; field kept for sweeps.

**Harness (rework A11's):** `alpha_stays_bound` sweeps
channeling {0.7, 0.9, 1.0} × nuclear_ambient {0, 2, 5, 10, 20} on
RigidAlpha carbon, ≥20k steps, pass = some combo holds every inter-alpha
spacing within ±30%; print the full drift table. If a combo passes: set
its values as the Couplings defaults, un-ignore the test, and run
`carousel_self_organizes` with the winners (report net carousel ω).
If nothing passes: report honestly — best combo, failure mode (eject vs
collapse vs lateral shear), and which force is unbalanced at the failure.
`neon_is_inert` and the whole molecular suite must stay green untouched.

### A12. Part 2 — HUD toggle

`atom_hud.tscn` gains a 3-state control (OptionButton or 3 buttons):
RigidLock / RigidAlpha / FreeNucleon → `atom_sim.set_nucleus_dynamics(i)`.
GDScript only calls the setter — no physics logic in GDScript. Remember
the Variant gotcha: use explicit types for values returned from
GDExtension calls, not `:=`.
