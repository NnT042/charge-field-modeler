# Field Calibration Mode — Design Plan

**Status:** Draft plan (2026-07-14). Not yet implemented. Author session: fresh start after
the session-35 bonding campaign.

## Why this mode exists

The Atom mode (M5) tried to abstract nuclear structure into emission-profile forces before we
had measured how a single particle actually responds to the field. That skipped steps, and the
result was a long campaign of tuning knobs against explosions we couldn't ground.

Field Calibration Mode is a **wind tunnel for one particle at a time**. It is a developer-only
instrument (never shipped to end users) whose entire job is to produce **numbers and formulae**
that the revived Atom mode can plug in with confidence:

- How far does a proton drift, and how fast does it turn, per unit time, under a given field?
- Does the equatorial vortex (pole intake → 30° emission → ricochet → re-intake) form on its
  own from the collision rules, or must it be imposed?
- What field density corresponds to "room temperature," so every later experiment runs at a
  standard, physically meaningful heat?
- **What is the *correct* outer-spin rate for a given environment?** Spin mode currently uses a
  sub-c outer rate chosen by eye ("this trace looks cool," see [[feedback_subc_top_level]]). If
  the outer spin is a live variable driven by collisions, its equilibrium value under a
  standard field is a *measurable* number — a chance to replace the eyeballed rate with a
  derived one.
- How does a neutron behave over its lifetime, and can we reproduce the ~15-minute decay as an
  emergent collision statistic?

It sits **beside** the spin-stacking Atom pipeline, not on top of it. It borrows the particle's
*shape* (the trace hitbox) from spin mode, but gives that shape mechanical freedoms the spin
mode deliberately froze.

## Relationship to existing code

| Reuse from | What we take | What changes |
|---|---|---|
| Spin mode (`spin_stack.rs`, `path_trace.rs`) | The swept-volume trace → disc-and-hole hitbox per particle type ([[project_hitbox_abstraction]]) | The particle is now a rigid body that can translate **and** re-orient; outer rotation becomes a tracked state, not a fixed display spin |
| M4 field sim (`field_sim.rs`, `collision_detect.glsl`, `impulse_sum.glsl`) | GPU photon buffer, collision-against-hitbox pass, impulse reduction | Photons gain **per-photon spin state** and **N-bounce persistence** (they don't die on first contact); chirality-dependent collision rules replace specular reflection |
| M5 Atom mode | Nothing structural — this mode *feeds* M5, it isn't built on it | — |
| `units.rs` | Natural↔SI conversion | Add a density↔temperature conversion (new) |

Deliberately **not** reused: the emission-profile force model, the flow network, the bond
knobs, RigidLock/RigidAlpha. Those are the abstractions we're trying to *derive*, so they can't
be inputs here.

## The particle model in calibration mode

A calibration particle is a rigid body with these state variables:

```
position      : Vec3        // free to move (collisions shift it)
orientation   : Quat        // free to turn (collisions torque it)
lin_velocity  : Vec3        // accumulated linear drift
ang_velocity  : Vec3        // tumbling of the whole body (pole wander)
outer_spin    : f64         // rate of the OUTERMOST spin level, in units of c
hitbox        : HitboxKind  // disc+hole (proton), thin disc (electron), polar-only (neutron)
emitter       : EmissionProfile  // where this body throws charge out (from spin-mode histograms)
```

Two mechanical freedoms that spin mode held back, now released:

1. **Translation & re-orientation via collisions.** Every photon impact transfers momentum
   (linear) and, because impacts are off-center tangential edge-hits, torque (angular). This is
   what lets us measure "how much a proton is nudged and flipped per second."

2. **Outer-rotation rate as a live variable (`outer_spin`).** This is the open design fork you
   flagged in note 2. Two candidate models — the box is built to test both:
   - **(A) Variable rate:** collisions can spin the outer level up or down (`outer_spin`
     changes). Physically motivated by Mathis's inverse-Compton/Compton rules (`neut2.pdf`,
     `halbach.pdf`): same-chirality head-on hits pump spin energy in, opposite-chirality hits
     bleed it out. This is required for the neutron spin-up experiment (note 7) to even be
     expressible.
   - **(B) Fixed rate + linear stacking:** the outer spin is pinned at its characteristic value
     and *all* momentum change shows up as linear/tumbling motion on top of a rigid rotor.
     Simpler, and may be what the data actually supports.

   **Decision (locked):** build (A). It strictly generalizes (B) — with the spin-coupling
   coefficient at 0 it *is* (B) — so we measure the coefficient rather than assume it. Two
   payoffs beyond the (A)/(B) question:
   - **Runaway/collapse as a tuning signal.** If `outer_spin` runs away to infinity or decays to
     zero under a physically reasonable field, that's a red flag that a collision rule or a
     coefficient is wrong — a built-in sanity check on the whole rig, not just a result.
   - **Deriving the "true" outer rate.** Where `outer_spin` *settles* under the standard
     `room_293K` field is a candidate for the correct sub-c outer-spin rate — the number spin
     mode currently eyeballs. If it settles cleanly and reproducibly, we feed it back and retire
     the "looks cool" value.

## The baked-loop hitbox (locked 2026-07-15)

The hitbox is not re-composed from 12 spin levels each tick. Instead each particle is
**one baked loop + one live swing**, exploiting that `compose()` applies levels inner→outer,
so the outermost level is always the last operation:

```
final_pos = R_swing(angle) · ( baked_loop_point + swing_offset )
```

- **Baked loop.** Freeze levels `1..=L` at ±c (that saturation *is* the particle's identity)
  and sweep the base photon through exactly one period of level `L`. Because amplitudes climb
  in exact powers of 2 (`types.rs::level_amplitude`), every inner level completes an integer
  (power-of-2) number of turns within that period, so **one period is a perfectly closed loop —
  the bake is lossless, not an approximation.** Store it once as a polyline in the pole frame.
  (Sampling caveat: the loop winds the fastest translating level `orbit_radius(L)` times, so a
  high-`L` bake needs thousands of points to resolve the fan-blade fine structure — that inner
  asymmetry is the torque handle, so it must be resolved, not smoothed.)
- **Live swing.** The next level up (`L+1`) rotates the baked loop about the pole (Z) at a
  rate that collisions drive up and down. This is the *only* per-tick geometry update besides
  the rigid-body pose. `outer_spin` is that swing rate.

The three particles, unified:

| Particle | Baked loop (`L`) | Live swing (`L+1`) | Swing kind | Behavior |
|---|---|---|---|---|
| Electron | L8 loop  | L9  | precession (offset 0)      | rides the current |
| Proton   | L12 loop | L13 | precession (offset 0)      | L13 sweeps loop → **equatorial disc, emits** |
| Neutron  | L11 loop | L12 | **orbital** (offset X·256) | non-rel L12 swings the L11 kite → **polar rake, traps charge** |

**The proton/neutron split is the swing *kind*, and it's load-bearing.** The proton's swing
(L13, A4 axial) is a *precession* — spin the loop in place about the pole → charge flung to an
equatorial disc. The neutron's swing (L12, Z-orbital) is an *orbital* — offset the loop out on
a `r₁₂ = 256` string and swing it like a kite → charge raked up and down through the poles and
recycled, never reaching an emitting equator. This reproduces `neutron.pdf` verbatim: *"no
normal emission at the equator (like a proton)... charge gets funneled back to the core... and
escapes back out one of the poles."* Transforms:

- Proton:  `pos = R_z(θ₁₃) · loop`                (precession, no offset)
- Neutron: `pos = R_z(θ₁₂) · (loop + X·256)`      (orbital kite-on-a-string)

### Neutron decay falls out for free — L12 first-passage to c

The canon neutron is a saturated **L11 loop swung by a *non-relativistic* L12, with no L13**
(charge trapped, no equatorial emission). Neutron → proton is therefore just **collisions
pumping that non-relativistic L12 up to c.** The instant L12 saturates, the spin stack's own
rule (`is_at_c()` → `activate_next()`, `spin_stack.rs`) **unlocks L13** — and L13 turning on is
exactly what sweeps the loop into the emitting disc. The particle is now a proton. "Emergent
L13" (decision 2a) is not bolted on; it is the automatic consequence of L12 crossing the c
ceiling. This is the same c-ceiling **transmutation event** flagged for the runaway problem —
for the neutron, *crossing it IS the decay.*

Lifetime, made simulable: **~15 min = the mean first-passage time for L12 to random-walk up to
the c barrier** under room-temp flux and the augment/cancel torque law (decision 1). L12 sits at
a stable non-relativistic equilibrium (augment ≈ cancel); rare fluctuations drive it to c. Where
that equilibrium sits relative to the barrier sets 15 min, and it must scale with density in a
specific way — the curve-fit that decides hypothesis (a) vs Mathis's single-positron-hit (b).

**Open, to measure not assume:** is the flip *sharp* (kinematics says the swept envelope is
rate-invariant, so pumping L12 barely changes shape until L13 snaps on at c) or *gradual* (the
earlier histogram-walk, which was the *post*-c L13-ramp regime)? The box resolves it, and that
answer *is* decision 2's open sub-question "is adding L13 geometrically legal?" sharpened to
"does L12 actually reach c under realistic flux, and does L13 auto-activate when it does?"

### Why runaway is bounded, and the energy scale

Inertia sets the *timescale*, not the equilibrium: a proton's moment of inertia about the pole
is `I ≈ m·R²` with `m_p = r_p²` and `R` = amplitude 512, so each charge-photon hit changes
`outer_spin` by a negligible amount — but tiny steps under a *net* torque still run away. The
actual regulator is the **uranus.pdf augment/cancel rule** (decision 1): same-chirality
head-to-head pumps spin up, opposite-chirality side-to-side bleeds it down. Isotropic field →
pump ≈ bleed → finite equilibrium (the measurable "true outer rate"); directional field → drift
toward sync. The torque law is self-limiting:

```
dω/dt = (1/I) · Σ τ_z(hit)   →   fluid limit:   dω/dt = k · Φ · f(ω ; field)
```
with `Φ` = field flux (density × relative speed, fed by the temperature knob), `f` a restoring
function with a fixed point at ω_sync, and `k` = the spin-coupling coefficient the mode *measures*
(never hardcoded). The c ceiling is not an overflow guard — it is saturation → next-level unlock
(= transmutation, flagged).

**Sourced energy scale (per CLAUDE.md constants rule).** The Dalton `D = 1821` is the
proton/electron mass ratio and equals `r_p/r_e`, with `m_p = r_p²` (elec3.html, proton3.html).
The stacked-spin energy ladder is `1, 9, 65, 1025, 16385` for no-spin→axial→x→y→z, and
`16385/9 = 1820.56 ≈ D` (elecpro.html); the proton-to-charge-photon factor is larger still —
many more rungs of the same doubling ladder. The user's "1821³" is an order-of-magnitude framing
of that hugeness; the exact torque coefficient `k` and `I` are measured outputs, not inputs.

## The ambient field model

### Density = temperature = heat (note 5)

Grounded in `heat.html` (*"heat is photon density"*) and `photon3.pdf` / "Redefining the Photon"
(average charge-field density derivable from c: `D = ∛(2r)·n/V`).

Plan:
- Expose **one primary knob: temperature**, backed internally by photon number density (photons
  per simulation volume).
- Calibrate the mapping so a labeled **"Earth room temperature ≈ 293 K"** preset corresponds to
  the Mathis-derived ambient density. *Implementation note:* pull the actual density figure from
  "Redefining the Photon" via `search_mathis` at build time — do not guess the constant (per
  CLAUDE.md: constants tied to charge theory need a source).
- Presets ladder: `cold_void` (near-zero), `room_293K` (default for all standard runs),
  `hot_1000K`, `solar_core` (stress test). Every calibration experiment reports which preset it
  ran at, so results are comparable.
- Composition stays from `PHYSICS_REFERENCE §5`: 2/3 photon / 1/3 antiphoton for Earth presets,
  50/50 for void.

### Per-photon spin momentum (note 6)

Each field photon carries its own spin state (axis + rate + chirality), not just a chirality
bit. This is the new ingredient that makes flows able to **spin each other up or down** and to
**spin up the focus particle**:

- When two flows meet, apply the `uranus.pdf` rule: opposing spins **head-to-head** → augment
  (spin-up); opposing spins **side-to-side** (same travel direction) → cancel (spin-down).
- Photon spin can transfer to the focus particle's `outer_spin` and `ang_velocity` on contact.
  The target observable: does a free particle placed in a directional field **slowly sync its
  rotation to the ambient and settle pole-down?** (note 6). That's a direct, watchable output.

## The collision & ricochet engine (note 4 — the hard part)

This is the heart of the mode and the biggest departure from M4. Today's field photons die on
first contact. Here they **persist for N ricochets** and only leave the accounting when they
exit the "immediate area" and rejoin the anonymous field.

The vortex hypothesis, made testable, is a pipeline of three collision events:

1. **Engine ejection.** The focus particle throws charge out along its emission profile — the
   equatorial disc and the pole→30° channels. (Emission directions come from the spin-mode
   histograms; we are *not* re-deriving them here, we're tracking what happens *after* ejection,
   which nothing has modeled yet.)
2. **First ricochet.** An ejected photon collides with an ambient field photon moving up/down
   the axis. The head-to-head vs. glancing geometry (rule above) sends both into complementary
   directions. **Measure:** what fraction of ricocheted charge lands inside the ~60°-wide dead
   zone above a pole intake?
3. **Second ricochet / capture.** A photon in that dead zone meets an axis-aligned photon
   descending the pole. Outcome: either it's turned *down the pole* (vortex closes — charge
   recycled) or flung to the cone wall (lost). **Measure:** the down-pole capture fraction.

If steps 2–3 produce a self-sustaining recycle without us imposing it, the Faraday-disc vortex
is **emergent** and we have the mechanism the whole project has been assuming. If it doesn't
close on its own, we learn exactly which collision rule or geometry is missing — which is just
as valuable.

**Two implementation strategies — build the analytic one first, validate against brute force
when the new GPU lands:**

- **(i) Probabilistic ricochet (default, card-independent).** For each tracked dot, generate a
  *probability zone* of where it most likely hits next (a function of local field density) and
  at what angle, sample a hit from it, and recurse N times. Cheap, runs on the current RX 550,
  and gets the charge-accounting logic correct without needing to render a million live balls.
- **(ii) Brute-force ricochet (ground truth, needs the new card).** Literally simulate ~10⁶ pool
  balls with full pairwise collision physics and watch where they go. Slower but assumption-free.

The plan: implement (i) first to nail the accounting and telemetry, then once the new GPU is in,
run (ii) on a few scenarios and **validate the probability model against the direct simulation**.
If they agree, we trust (i) for fast iteration and keep (ii) for spot-checks; if they diverge,
(ii) tells us what the analytic zone got wrong.

Engineering shape (shared):
- Photons get a `bounce_count` and a `home_cell` (spatial-hash cell of origin). They're
  simulated with full collision physics while within R of the focus particle; past that R and
  past `bounce_count == N`, they're retired back to the ambient pool.
- A spatial hash (already flagged as needed in PROJECT_DESIGN open-question 4) is a prerequisite
  for (ii) and useful for (i)'s density lookups — build it, not an afterthought.
- Start with small N (2–3) and a modest particle count to get the accounting right before
  scaling.

## Measurement & telemetry (note 3 — the actual point of the mode)

The mode is worthless without rigorous, time-averaged readouts. Every experiment logs, over a
configurable time window:

- **Drift:** mean & variance of Δposition per unit time (the "how much does it shift" number).
- **Turn:** mean & variance of Δorientation per unit time (pole wander rate); and time-to-align
  when a directional field is on.
- **Spin response:** `outer_spin` trajectory vs. time; the spin-coupling coefficient from the
  (A)/(B) test above.
- **Vortex accounting:** ejected-charge fate breakdown — % re-captured down-pole, % to cone
  wall, % escaped; dead-zone occupancy fraction.
- **Collision census:** hits/sec by region (disc, pole, dead zone), by outcome (augment/cancel/
  deflect), by chirality match.

All exportable to CSV/JSON (M6 already wants this) so the numbers can be curve-fit offline into
the force laws Atom mode will use. **The deliverable of this mode is a table of fitted constants,
not a pretty simulation.**

## Experiment ladder (notes 7 & 8)

Run in order; each depends on the previous being calibrated.

1. **Proton, isotropic room-temp field.** Baseline drift/turn. Should be ≈0 net drift by
   symmetry — a correctness check on the whole rig. Confirm the 30° emission and equatorial disc
   are consistent with spin-mode predictions ([[project_30deg_validation.md]]).
2. **Proton, directional field.** Measure drift-down-gradient and alignment rate. This is where
   the vortex recycle (note 4) is tested.
3. **Electron.** Too large to pass proton channels but pushed by the field
   (`PHYSICS_REFERENCE §7.3`). Measure how it rides the current. Confirm it can't self-emit a
   coherent disc.
4. **Neutron — the full campaign (note 7).** Polar-only emission; weak north pole.
   - Measure lifetime-relevant statistics: mean time between "right kind" of pole hits at a
     given field density.
   - **Three decay hypotheses, tested head-to-head — the divergence from Mathis is
     *histogram-driven*, not a whim.** Our own spin-path histogram work shows that **adding
     level-13 (uberon-tier axial) motion to a Mathis-configured neutron shifts its emission
     histogram markedly toward a proton disc.** That is the empirical reason to keep an
     alternative to Mathis's single-hit model on the table. The three:
     - *(a) Emergent spin-up (your primary model):* does sustained collision slowly drive the
       stack toward starting a level-13 axial rotation? **Open sub-question: is adding level 13
       even geometrically legal with the spins in a Mathis-approved neutron configuration?** The
       box is where we find out — build a neutron, try to energize level 13, and watch whether
       the histogram walks toward the proton disc or destabilizes.
     - *(b) Mathis's hit model (`quark.html`, `electron.pdf`):* a single glancing positron/lepton
       hit on the weak north pole reverses the z-spin (level-12 chirality), flipping
       neutron→proton in one event.
     - *(c) Fallback / cross-check:* we already have a **positron emission histogram** — spawn a
       positron body and literally fling it at the neutron's north pole, and see if the collision
       reproduces (b) mechanically rather than by fiat.
   - Whichever reproduces ~15 min at room-temp density is the mechanism. Sweep density and check
     the lifetime scales the way the winning model predicts. **This experiment is the single
     strongest reason to build the mode** — nothing else we have can express it.
5. **Hydrogen (proton + electron).** Electron capture into the polar eddy; measure standoff
   distance and orientation purely from the field, no coded bond.
6. **H₂.** Two protons — do they orient equators-facing on their own (`PHYSICS_REFERENCE §7.4`)?
   Measure bond distance and angle.
7. **Deuterium** (p + n bound) and **Helium/alpha** (2p + 2n). Read off distances and pole
   orientations from watching them settle. If the field numbers from experiments 1–6 are right,
   these should assemble without per-configuration tuning — the test that we're finally ready to
   go back to Atom mode.

## Decisions (locked 2026-07-14)

1. **Outer-rotation model → variable-rate (A).** Measure the spin-coupling coefficient rather
   than assume it. Treat runaway/collapse of `outer_spin` as a tuning signal, and treat its
   settled value under `room_293K` as a candidate for the correct outer-spin rate (see model
   section above).
2. **Neutron decay → test three hypotheses, don't bake in an answer.** Emergent level-13 spin-up
   (primary, motivated by the histogram walking toward a proton disc), Mathis's single positron
   hit reversing the z-spin, and a positron-histogram collision cross-check. The open geometric
   question — whether level 13 can even be added to a Mathis-legal neutron — is itself an
   experiment output.
3. **Ricochet → start shallow (N=2–3), analytic first.** Build the probabilistic-zone model on
   the current card; validate against brute-force simulation once the new GPU arrives.
4. **Location → new mode in the same app, dev-gated behind a flag,** with a **quick toggle/hotkey**
   so it's fast to re-enter during iteration (this is a heavily-iterated instrument, not a
   one-shot view).

## Build order (when coding starts)

- **CM-0:** Mode scaffold + dev flag; load a trace hitbox as a free rigid body; single-photon
  collision → linear + angular impulse. (Proves the freed-DOF particle works.)
  - *Started 2026-07-15:* `hitbox.rs` — `bake_loop(L)` sweeps levels `1..=L` (saturated) into a
    closed polyline; `BakedLoop::swept_point()` applies the live swing (`R_z(θ)·(loop+offset)`),
    per the "baked-loop hitbox" section above. Proton `L=12`/swing L13, neutron `L=11`/swing L12,
    electron `L=8`/swing L9.
- **CM-1:** Density=temperature calibration + presets; pull the Mathis density constant.
- **CM-2:** Per-photon spin state + augment/cancel collision rule; the "sync to ambient" readout.
- **CM-3:** N-bounce ricochet engine + spatial hash + vortex accounting telemetry.
- **CM-4:** Full time-averaged measurement panel + CSV/JSON export.
- **CM-5:** Experiment ladder 1→7, each as a saved scenario with logged outputs.

Each CM step ends with a committed diagnostic and a numbers table, same discipline as the
session-3x campaign — but this time the numbers feed forward into Atom mode instead of chasing
an explosion.

## Model A: velocity-channel drag (implemented)

**Measurement history.** The orientation-participation diagnostic showed the hit-weighting
channel (which orientation the loop presents to the field) is dead as a self-limiter — it barely
moves the histogram even at high `outer_spin`. `swing_drag_torque` found the real limiter
instead: model B's gear-catch pump (`apply_photon`) computed `tau_perp` and *discarded* the
pole-axis component of the plain momentum torque (`r × j`). Restoring it, weighted by the
photon's speed **relative to the moving swing surface** (`|c·dir − v_surface|`, c = 1), gives a
torque that grows with `outer_spin` and is signed to oppose it — a genuine drag, not another pump.

**The rule.** In `apply_photon`, after the chirality pump: `tau_pole = (r × (dir · momentum ·
catch)) · pole`, `outer_spin += spin_coupling · tau_pole / i_spin`. `catch = |c·dir − v_surface|`
where `v_surface = ω × r` is the swing surface's velocity at the contact point. Head-on photons
have a larger relative speed than co-moving ones can chase, so an isotropic field's mean torque
is nonzero and spin-down — unlike the chirality pump, this channel owes nothing to chirality.

**The scaling law.** Equilibrium is where the augment/cancel pump balances the drag:
`0.5·spin_gain·(2p−1) + spin_coupling·momentum·swing_drag_torque(s*) = 0` (p = photon_fraction).
`report_spin_equilibrium` bisects that relation per field setting and compares it against a
ticked `AmbientField` simulation's actual settled `outer_spin`. Across an imbalance sweep
(p = 0.5…1.0, flux 200) and a density sweep (flux 50/200/800, p = 0.667), predicted and simulated
equilibria agree to within a few percent, and the density sweep confirms the settled value is
flux-independent while `settle_step` (time-to-settle) shrinks with flux — exactly what a
density-independent equilibrium with density-dependent approach rate predicts.

**Thermal floor.** At rest, or for any single photon, the drag term is generally nonzero — only
the *isotropic mean* vanishes (`drag_has_zero_mean_at_rest`). Per-photon scatter around zero is
expected noise, not a bug; a balanced field's simulated settled spin jitters near zero for the
same reason (see the p=0.5 row above).
