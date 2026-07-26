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

## CM-1 second half: SI flux calibration (implemented)

**Sourced anchors (`units.rs`).** `e = 1.602e-19 C`, and `1 C = 2e-7 kg/s` (SI Ampere
definition), so the elementary charge recycles mass at `3.204e-26 kg/s` (xpart.pdf, fine3.pdf
"Fine Structure Constant") — the proton recycles ~19.2x its own mass in charge photons every
second. Dividing by the charge-photon mass `2.75e-37 kg` (photon3.pdf `D/n`) gives the proton's
true recycle flux, ~1.16e11 photons/s; that photon mass cross-checks against photon.html's
`m_p/1821³ ≈ 2.77e-37` to <1%.

**Macro-photon aggregation.** The sim can't tick 1e11 contacts/s, so one sim photon stands in
for `K` real photons: `AmbientField::from_si` scales the per-photon magnitudes (momentum,
spin_gain) UP by `K = sim_momentum / true_momentum_natural` while the contact rate stays at the
chosen sim `flux`. Because the chirality pump and the velocity-channel drag are both linear in
per-photon magnitude, the settled equilibrium `outer_spin` is K-invariant — only the noise floor
coarsens as K grows (`report_si_calibration` Table 3 checks this directly, <15% tolerance).

**Time mapping.** `seconds_per_time_unit = sim_flux · K / flux_si_hz` converts natural
simulation time into SI seconds by matching the real momentum-current the sim run represents.

**`report_si_calibration`** (ignored, `--nocapture` to read) prints: Table 1, the sourced
anchors and derived `K`/time-conversion per particle at Room density; Table 2, the settle
trajectory (same 150k-step loop as `report_spin_equilibrium`) run through `from_si`, with
settle time converted to SI seconds; Table 3, the K-invariance check for the proton.

## Alignment channels measured (2026-07-16) — shadow aligns, intake goes broadside

`directional_torque` already established the clean negative: the velocity-catch channel
produces ZERO alignment torque at every tilt. Four candidate channels were then measured
(diagnostics only, no rule change): `surface_shadow_torque` (downstream self-occlusion of
the surface channel, chirality-blind) and three chirality-gated pole-intake bookkeepings
(`intake_mouth_torque` absorb-at-mouth, `intake_channel_torque` channel-to-core +
equatorial ring emission, `intake_through_torque` pole-to-pole through charge — the last
reduces exactly to the first because the exit recoil is parallel to the pole). Sweep:
`report_alignment_channels`, proton bake, tilt 0-180 deg vs a fixed stream.

**Results.** SHADOW is the standing-up channel: stable attractors at 0 AND 180 deg (pole
parallel to field, either polarity), repeller at 90 — the dielec.pdf "stand straight up"
torque, present at all outer_spin values tested. The INTAKE bookkeepings are all
weathervane-mechanical: an aperture that absorbs momentum is pushed downstream, so every
variant is STABLE BROADSIDE (90 deg) — no absorption bookkeeping can turn a mouth INTO
the stream (hof.pdf's feeding polarity). In a photon-rich field the intake channel also
destabilizes the south-mouth-upstream pole twice as strongly as the north (2:1 mix).

**Open question (polarity):** shadow alignment is polarity-degenerate; Mathis's
south-mouth-into-the-stream selection needs a mechanism beyond passive absorption —
likely a suction/pressure-deficit term (capture REMOVES ambient pushback at the mouth,
flipping the sign of the mouth torque). Model and measure before wiring any intake rule.

## Shadow occlusion wired + full-velocity drag (2026-07-16)

The two measured-but-inert channels above (`surface_shadow_torque`'s alignment restoring
force, `swing_drag_torque`'s velocity-dependent drag) are now live collision rules, not just
diagnostics.

**Rule 1 — exposure-weighted contact sampling.** `AmbientField` gained `pub occlusion: bool`.
When true, `tick` no longer draws a contact point uniformly from `p.world_points()`; it draws
the photon direction `dir` first, then REJECTION-SAMPLES the contact so a candidate point is
accepted with probability `expose = max(0, -dir·r̂)` (`r̂` from the body's center to the
candidate). Points facing the incoming stream are likelier to be hit, points on the
self-occluded downstream side are unlikely — the same Lambert-shadow law
`surface_shadow_torque` already measured. Hit COUNT per tick is unchanged (still exactly `n`);
only which point gets credited shifts. Capped at 64 draws with a last-candidate fallback for
degenerate geometry. `occlusion = false` reproduces the old uniform draw bit-for-bit (used by
`room_default`/`from_si` as the new default `true`, but pinned back to `false` in
`report_spin_equilibrium`/`report_si_calibration` — see below).

**Rule 2 — full-contact-velocity drag catch.** `apply_photon`'s model-A drag catch,
`|c·dir − v_surface|`, used to compute `v_surface` from the swing rotation alone
(`omega × r`), so a tumbling or drifting body felt no drag on those channels — an aligning
particle driven by Rule 1's restoring torque would have oscillated forever, an undamped
pendulum. `v_surface` is now the FULL local material velocity (`contact_velocity`: drift +
tumble + swing), and the catch-weighted "drag impulse" is split the same three ways the plain
Newtonian impulse already is: pole component → `outer_spin` (unchanged from before), transverse
component → `ang_velocity` (pendulum settling), full vector → `lin_velocity` (drift decelerated
toward the field rest frame). Isotropic mean of a catch-weighted impulse is anti-parallel to
whatever velocity feeds it (`E_dir[dir·|dir−v|] ≈ −v/3` for small `v`), so this is genuine
damping on all three channels, not another pump. **Numerical note:** unlike `outer_spin`
(clamped to ±1 = c), `ang_velocity`/`lin_velocity` had no speed ceiling, and `catch` GROWS with
`|v_surface|` once it's no longer small relative to c — an explicit-Euler feedback loop that
diverges to NaN over a long run if nothing caps it. `v_surface` is now clamped to length ≤ 1
before computing `catch` (the same physical ceiling `outer_spin` already respects).

**Tests.** `occlusion_biases_contacts_upstream` confirms the exposure bias directly via a new
`AmbientField::sample_contact` helper (mean of `-dir·r̂` over accepted contacts: clearly
positive with occlusion on, ~0 with it off). `tumble_damps_in_balanced_field` and
`drift_damps_in_balanced_field` confirm Rule 2 actually damps existing `ang_velocity`/
`lin_velocity` in a balanced (no-pump) isotropic field. `report_standing_up` (ignored) is the
headline: a proton tilted 45° in a directional Room-mix field with occlusion on settles into a
fluctuating near-zero-tilt band (this is a stochastic system — think Kramers double-well
hopping between the 0° and 180° attractors at high noise, not a clean deterministic settle;
tuning `momentum` down and `flux` up improves the drift/diffusion ratio, since diffusion grows
with `momentum²` but coherent drift only with `momentum`) while the isotropic control shows no
directional preference at all.

**Shadow-induced equilibrium shift.** `report_spin_equilibrium` now runs its Earth-mix
flux-200 settle twice — occlusion off vs on — everything else identical. Unshadowed settled
`outer_spin* = 0.3633`; shadowed `= 0.4963`; **relative shift ≈ 37%**. This is a genuinely
large shift (occlusion changes which points a random-direction photon can hit, which changes
the chirality-pump statistics even under an isotropic-direction draw), so `report_spin_
equilibrium` and `report_si_calibration` pin `occlusion = false` to keep validating against the
analytic `swing_drag_torque` prediction (which models the unshadowed channel); the shadowed
number is now on record here rather than silently baked into those reports.

## CM-2: the gear-tangent pump — zero fitted constants (2026-07-16)

`apply_photon` carried two scaffold artifacts: `spin_gain` (a free coupling with no source)
and a drag term that double-counted the plain Newtonian impulse (added once unweighted, then
again catch-weighted). uran.pdf/uran4.pdf source the fix directly — "the spin force is at the
tangent": a photon's edge spins at c while it travels at c, so the tangential spin force at the
contact tooth has the SAME magnitude as the linear force. halbach.pdf adds that spin energy
transfers to linear "on the tangent," with no separate spin-energy budget. So the chirality pump
isn't an independent coupling; it's the photon's own momentum redirected along the surface
tangent, signed by the chirality mesh. `apply_photon` is now exactly two terms: (1) a
catch-weighted Newtonian transfer (linear + transverse + pole, replacing both the old plain
impulse and the old stacked drag — this term IS model A, and reduces to the plain impulse
exactly at rest where `catch = 1`), and (2) a gear tangential transfer of the SAME per-hit
momentum along the pole-positive tangent, signed by chirality and scaled by the catch factor
`f` — this term IS the chirality pump, now carrying its own geometric lever (`r × t̂`) instead of
a lumped gain. `gear_efficiency` (dimensionless, ships 1.0 — the uran.pdf-sourced value) replaces
`spin_gain`; it survives only as an ablation knob, not a magnitude to fit.

**Equilibrium is momentum-free by construction.** Both terms in `apply_photon` are scaled by the
SAME `spin_coupling / i_spin`, and the pump additionally by `momentum · gear_efficiency` where
the drag has just `momentum` — since `gear_efficiency = 1.0`, momentum and the inertia coupling
cancel out of the equilibrium condition entirely. The bisected equation is now
`(2p−1)·swing_pump_torque(s*) + swing_drag_torque(s*) = 0` — pure geometry and photon:antiphoton
mix, nothing else. `swing_pump_torque` is the new analytic counterpart to `swing_drag_torque`,
sampled the same way (isotropic Fibonacci-sphere directions, swing-phase sweep, optional Lambert
exposure weighting via `occluded`); both `swing_drag_torque` and `swing_pump_torque` gained that
`occluded` parameter so `report_spin_equilibrium` can bisect both the unshadowed AND shadowed
channels per particle.

**PRESETS CANDIDATE (Earth mix 2/3 photon, occlusion on, flux 200)**: `proton = 0.9996`,
`neutron = 0.1625`, `electron = 0.9973` — proton and electron settle to nearly the SAME value, a
strong confirmation that their normalized drag curves are identical as this refactor predicts;
the neutron's non-axisymmetric kite-offset swing sits in a completely different regime. Both
proton and electron TRANSMUTE at Earth mix under these settled numbers — worth flagging rather
than treating as expected, since a tied (momentum-free) pump was not obviously guaranteed to
still saturate at c for both. Density-independence survives intact (proton settled `0.9996` at
flux 50/200/800, bit-identical to 4 decimals). K-invariance is now close to a tautology
(`gear_efficiency` doesn't scale with K at all) and the empirical check passes at 0.5% relative
difference (tolerance 5%). The neutron's predicted-vs-simulated agreement is notably worse than
proton/electron's (~3-4% off) — off by ~80% at several photon fractions — because its kite-offset
orbital swing breaks the isotropic swing-average the bisection assumes far more than the clean
axisymmetric precession swings do; flagged for follow-up, not patched here.

`photon.html`'s size-differential edge-hit rule (bigger loops catch harder) is realized
geometrically now, not probabilistically: `f` (catch factor), the lever `r × t̂`, and occlusion's
exposure weighting are all measured from the baked geometry directly — nothing hardcodes which
particle "should" catch more.

## CM-2 self-limiter wired: DirRelSign × ladder (2026-07-16, session 36)

The gear-tangent pump above was momentum-free but NOT self-limiting: at `gear_efficiency = 1.0`
it overruns the catch-weighted Newtonian drag at Earth mix and both proton and electron transmute
under an ordinary ambient field (flagged, above) — falsified, since real matter is stable in
ordinary starlight. A measurement campaign (commits 663d8d8 "Pump self-limiter measured: ladder
efficiency x direction-relative sign wins" and 6301db2 "Critical gear efficiency measured: baryon
rest is a TRUE zero") tried two independent families of fix — (A) an externally-imposed
spin-energy-ladder efficiency, (B) making the gear rule itself direction-relative — crossed
against each other, and the user signed off on the combined rule. It is now wired into
`apply_photon` directly (no longer measurement-only).

**The rule.** `apply_photon` term 2's chirality sign is no longer `chirality` alone; it's the
DirRelSign EFFECTIVE chirality, `chi_eff = chirality * (-dir·t_hat).signum()`, where `t_hat` is
the actual signed swing tangent at the contact if the body is spinning, else the geometric
positive tangent `t_pos` at rest (pole.pdf/bright.pdf: "a photon going in reverse is automatically
an antiphoton"). This exactly mirrors the measurement diagnostic `gear_pump_variant_torque`'s
`PumpRule::DirRelSign` arm. Two properties fall out of it, both measured not assumed: the pump is
EXACTLY zero at rest in an isotropic field (population cancellation — half the isotropic photons
see the fallback tangent as opposing, half as co-moving, in equal measure, to machine epsilon),
and it plateaus once spinning rather than growing without bound — a structural self-limiter that
owes nothing to an externally-tuned magnitude.

The magnitude channel is a new per-particle field, `CalibrationParticle::gear_efficiency`
(dimensionless `E_photon / E_platform`, elecpro.html/higgs3.pdf's spin-energy ladder
`1, 9, 65, 1025, 16385` for no-spin → axial → x → y → z). Under the **locked-platform reading**
— the baked levels are locked in and deliver the full force of THEIR ladder rung; only the live
top swing is the free variable — the sourced constants are:

| Particle | `GEAR_LADDER_*` | Ladder rung |
|---|---|---|
| Electron | `1/9`     | axial term |
| Proton   | `1/16385` | full stacked sum |
| Neutron  | `1/1025`  | one named rung below the proton's (candidate) |

This is a property of the TARGET particle's own locked platform, not the field, so it lives on
`CalibrationParticle` (installed by the `proton`/`neutron`/`electron` constructors), not on
`AmbientField` — `AmbientField::gear_efficiency` is gone.

**Baryon rest is a TRUE zero, not a small number.** Because both the DirRelSign pump and the
drag vanish at rest, whether a particle spins up at all is a threshold contest:
`eff_crit := min over s>0 of -drag(s)/pump(s)`. Below `eff_crit` the net torque is negative at
EVERY `s > 0`, so the deterministic equilibrium is exactly zero. Measured (6301db2), Earth mix
(`p=2/3`) | pure photon (`p=1.0`): proton `eff_crit = 1/908 | 1/2723`, neutron `eff_crit =
1/454 | 1/1361` (roughly 2x easier to pump than the proton), electron `eff_crit = 1/908 | 1/2723`
(identical to the proton to 4 digits — the proton==electron normalized-curve finding again). Both
`1/16385` (proton) and `1/1025` (neutron) sit far below their respective `eff_crit` at every
mix — the proton rests at TRUE zero at every ambient mix up to and including pure photon; the
neutron rests at Earth mix but crosses into a tiny live equilibrium (`s* ≈ 0.0007`, needs a
fine near-zero grid to resolve — a uniform 200-point scan over `[0, 0.995]` is too coarse) only
at `photon_fraction = 1.0`, where `1/1025 > eff_crit(1/1361)`.

**Electron anchor holds.** At Earth mix, occlusion off: `s* ≈ 0.0504` (predicted, bisected against
the analytic `gear_pump_variant_torque`/`swing_drag_torque` channels) vs `s* ≈ 0.0498` (simulated
`AmbientField` settle) — within 1.2% of each other, and both within ~6% of the `~0.055c` outer
structure measured in the M2-era traces, with zero fitted constants (Mathis's `~0.0057c` is the
electron's LINEAR speed — a different observable, not the right anchor for the swing rate).

**Verification.** `report_spin_equilibrium`'s bisection helper (`predict_equilibrium_dirrelsign`)
now scans for the first DOWNWARD (stable) crossing of `gear_efficiency · gear_pump_at_spin(s, p,
DirRelSign) + drag_at_spin(s)` rather than plain-bisecting `mix · Baseline_pump(s) + drag(s)`,
since population cancellation makes `net(0) ≈ 0` ALWAYS true under DirRelSign — an upward
crossing at `s = 0` is a repeller, not the answer (same convention `find_equilibrium_label`
already used for `report_pump_self_limiter`/`report_critical_gear_efficiency`). Where the
prediction is a true zero, the report asserts the simulation's settled `|outer_spin|` against a
small noise floor rather than a relative error against exactly 0.0. K-invariance and
density-independence checks (`report_si_calibration` Table 3, `report_spin_equilibrium`'s density
block) now key off the ELECTRON rather than the proton — the proton's settled value is ~0 at
every K/density under its wired ladder rung, which would make either check trivially "pass"
without exercising the momentum-scaling identity at all; the electron's finite live equilibrium
is the particle these checks can actually discriminate on (and needed averaging over ~16
independent seeds to separate real K-invariance, confirmed at 2.2% relative difference, from the
noisier single-seed estimate a smaller-amplitude live equilibrium carries).

**Open question (unchanged from 663d8d8).** Which named rung governs the neutron's live L12
swing — the full `1/16385` or the rung below its bake, `1/1025` — decides the neutron decay
first-passage rate (β-decay as a ~1e14-contact rarity crossing L12 to c). `1/1025` is the
locked-platform CANDIDATE adopted here (one rung below the proton's full stack), not a settled
answer; discriminable later by a decay first-passage test comparing predicted mean lifetimes
under each candidate against the measured ~15 minutes.

## The axial-cone intake branch (2026-07-19, commit 7516470)

**The knife-edge.** `apply_photon` term 2's sign, `chi_eff = chirality * (-dir·t_hat).signum()`,
reads `t_hat` (the swing tangent) against `dir` (the photon direction). At exact pole-axis
incidence — `dir` parallel or antiparallel to the pole — `t_hat` lies exactly in the plane
perpendicular to the pole, so `dir·t_hat == +/-0.0` and `signum` evaluates IEEE signed-zero
bookkeeping instead of geometry. `report_vortex_feedback_probes` (B1-B3) confirmed the pathology
directly: the sign at exact axial incidence is floating-point-bit-pattern-dependent, theta=180
does not mirror to the opposite attractor, and theta=1/5 degrees do not connect continuously to
the axial value.

**The loop-resolved measurement.** `report_axial_channel_geometry` (commit 7516470) resolved what
the sign SHOULD be by replaying the same collision against the baked loop's own tangent field
(the base photon's actual closed path at c=1) instead of the swing tangent — a vector generically
NOT confined to the pole-perpendicular plane, so the dot product against an axial `dir` is
well-defined and non-degenerate. The result: axial incidence delivers a chirality-signed pole
torque, per entry lane (north = photon travels toward `-pole`, i.e. enters at the north pole;
south = travels toward `+pole`) and per particle, and the response is measured FLAT out to ~15
degrees off-axis before the swing rule's own geometry takes over.

**Lane constants.** Each pair is the section-[2] mean pole torque per caught photon at `s=0,
chi=+1`, normalized so the stronger lane is exactly `±1.0`:

| Particle | raw T (north / south) | `(AXIAL_GEAR_*)` north, south |
|---|---|---|
| Proton   | 22.264 / 35.002 | (0.636, 1.0) |
| Electron |  1.082 /  4.799 | (0.226, 1.0) |
| Neutron  |  5.048 / -3.798 | (1.0, -0.752) |

The neutron's sign flip is real, not a normalization artifact — the kite torques OPPOSITELY per
entry lane. Sourced from halbach.pdf (photons moving along the pole carry spin in the right
plane; edge hits land at the tangent) and pole.pdf (the polar intake is a charge engine).

**The cone.** `apply_photon` blends smoothstep-wise from the axial lane constant inside a 10
degree half-angle cone around the pole to the untouched oblique `sign_oblique` swing rule outside
a 20 degree half-angle cone, using `mu = dir·pole` to pick both the blend weight and the
north/south lane. Beyond 20 degrees the blend weight is exactly `0.0`, so the oblique path is
bit-identical to the pre-cone rule — this branch adds axial-incidence physics without perturbing
anything the swing rule already covered. The axial branch's magnitude still rides the same
`gear_efficiency` ladder rung as the oblique branch (see the table above the self-limiter
section); only the sign/lane-ratio channel is new. `report_vortex_feedback`'s THETA=0 CAVEAT is
resolved by this change; theta=0 rows produced before it remain untrustworthy history.

---

## The electrical / magnetic split, and the electron class (2026-07-24)

Two CM-3 additions in `rust/src/recycling.rs`. Both are measurement-side; the gear
geometry signed off in commit `ebf637e` is unchanged, and the committed
`histogram_{proton,neutron}_recycling.csv` are byte-identical after the refactor
(the export reproduces eqEsc 0.585 / 0.308 and peaks 5 deg / -46 deg exactly).

### ParticleClass replaces BaryonClass

`BaryonClass {Proton, Neutron}` became `ParticleClass {Proton, Neutron, Electron}`,
and the exit-latitude profile became a per-class property (`exit_lat_profile`,
`exit_peaks_deg`, `exit_sigma_deg`) threaded through the TAU calibration
(`tau_in_integral`, `tau_self_emitted_integral`) and `CollisionField`.

The electron is NOT a small baryon and the profile difference is structural, not
fitted. It is level 9 - baryon-tier AXIAL spin only, no x/y/z (PHYSICS_REFERENCE
section 3b; photon.html "the electron is spinning axially; the proton is spinning
axially plus x, y, and z"). One spin traces ONE exit family, so the bimodal
7/58 disc is unavailable to it by construction. Input profile: single Gaussian at
+-3 deg, sigma 2.5. It also has no level-12 z-spin at all, so its chirality tag comes
from the axial spin sense - an orientation, not a class invariant.

Bit-identity for baryons is deliberate and load-bearing: the peak-selector RNG draw
is consumed for every class (for n=2 peaks the index reproduces the old `u < 0.5`
branch), and `exit_lat_profile` accumulates one term at a time and divides rather
than multiplying by 0.25, so the two-peak sum matches the old expression to the bit.
An ULP shift there would propagate through sigma_e and could flip a collision coin.

### Photon.spin - the magnetic degree of freedom

`Photon.chirality` was a frozen tag that `scatter_off` read but never wrote, so the
old `chi_escaped_sum` was identically `class_sign x escaped` and could not show spin
cancellation at all. That is why the magnetic-neutrality item stalled: the hook was
measuring nothing.

`Photon.spin` is now mutable, initialized to the chirality and updated per collision
by ONE rule - the photon picks up the partner's spin, `spin += spin_transfer *
chi_partner`, clamped to +-2. Same sign stacks, opposite cancels; no case split.
Sourced from magmom.pdf ("those spins can stack or cancel") and voyag.pdf ("the
spins will offset as a sum, and the field will be magnetically flat", with "a
leftover or resultant field spin" when the populations do not balance). No new RNG
draws are consumed, so gear-on geometry is unchanged and gears-off stays a strict
no-op (`spin_transfer` 0 in `GearConfig::off()`).

The split this enables is the one Mathis insists on: magnetism is SPIN and is
geometry-blind, electricity is the linear motion c and is a genuine VECTOR sum, and
"photon collisions affect only the spins, not c" (neutron.pdf).

### report_magnetic_neutrality (ignored, ~87s, 16 cells)

class x photon_fraction {0.5, 2/3} x spin_transfer {0, 0.1, 0.25, 0.5}, at
cancel_redirect 0.5.

CONFIRMED, and knob-INDEPENDENT:
- The electrical channel is provably collision-independent. Incoherent z/xy ratio is
  0.517-0.518 for the proton and 1.032-1.042 for the neutron across ALL spin_transfer
  values - exactly as "collisions affect only the spins, not c" requires. Proton flux
  sits in the equatorial plane, neutron flux is axial, factor ~2.0 apart.
- Both channels NULL-CONTROL at the balanced field: at photon_fraction 0.5 the proton
  and neutron are indistinguishable in both observables (z/xy 0.741 both; magnetic
  ratio 1.000 at every spin_transfer). Matches neutron.pdf's "balanced field ->
  indistinguishable".
- spin_transfer = 0 gives M/esc = exactly +-1.0000, verifying the channel is wired and
  nothing else moves spin.
- The proton/neutron magnetic asymmetry exists ONLY in an imbalanced field. Mechanism:
  the proton's monochiral chi=+1 exhaust STACKS with the photon-majority ambient
  (M_p/esc rises 1.00 -> 1.59) while the neutron's chi=-1 exhaust CANCELS against it
  (|M_n/esc| falls toward 0.88).

NOT CONFIRMED - the magnitude is a FIT, not a derivation. |M_p/M_n| at Earth 2/3 grows
monotonically with spin_transfer (1.000 / 1.218 / 1.438 / 1.799) and passes through the
experimental |mu_p/mu_n| = 1.4600 at spin_transfer ~= 0.25. Since 1.46 sits mid-sweep and
nothing independently pins spin_transfer, only the ORDERING and the FIELD-DEPENDENCE are
results. magmom.pdf reaches 1.458 by a different route (an angle differential plus the
proton's faster spin), which this model does not implement.

REFUTED for the current architecture - balanced-field neutrality. neutron.pdf demands
"near-total spin cancellation ... the neutron is then neutral regarding magnetism" in a
balanced field, with charge still recycled. Measured: |M_n/esc| stays at 1.00-1.31 at
photon_fraction 0.5 for every spin_transfer, with escaped ~38.4k. Diagnosis: at balance the
ambient mean chirality is 0, so the meshing drift vanishes by construction and the
neutron's own monochiral exhaust survives untouched. Genuine neutrality needs AMBIENT
charge to enter one pole and leave the other RETAINING ITS SPECIES - venus2.pdf's
"through charge" - which the pump cannot express: it absorbs ambient and re-emits
monochiral charge from the fixed disc. This is the SAME architectural limit that blocked
the polar-exit profile until the through-channel moved the redirect to the march level.
The fix is a real pass-through limb, not a knob.

### report_electron_and_swing_control (ignored, ~83s, 6 cells)

class x swing_boost {0.05, 0.0} at Earth 2/3, cancel_redirect 0.5.

- SWING CONTROL: INSENSITIVE for all three classes (mean |delta| per bin 0.010-0.020,
  delta-eqEsc <= 0.002). `DEFAULT_SWING_BOOST` = 0.05 is the ELECTRON's measured settled
  outer_spin while CM-2 measured baryons at a TRUE ZERO at Earth mix, so the exported
  baryon CSVs had been generated with a swing the baryons should not have. It turns out
  not to be load-bearing - the proton disc and neutron butterfly survive swing = 0
  unchanged, so the exported CSVs STAND. Worth keeping on record as a resolved doubt.
- ELECTRON VALIDATED against its own trace profile. The input was a bare +-3 deg core with
  NO mid-latitude plateau; the engine REGENERATED the trace's shoulder - engine 0.383
  vs trace 0.359 over 45-70 deg from the pole (~7%), with overall mean |delta| = 0.057 across
  the folded profile, much closer than either baryon (proton 0.276, neutron 0.196). So
  the electron trace's mid-latitude tail IS collective scattering, and the collective and
  single-particle derivations AGREE for the electron. Exit topology: eqEsc 0.687, midEsc
  0.267, polEsc 0.047, peak at the equator - a clean single disc, as one axial spin should
  give.
- Residual: the engine puts ~3% of electron emission inside 25 deg of the pole where the trace
  is exactly 0.000. Same filled-polar-cone signature the baryon profiles show; it matters
  because Atom mode's binding story is a wall-riding electron orbit in an emission-free
  polar tunnel.

`histogram_electron_recycling.csv` is now exported alongside the two baryon files, so all
three Atom-mode profiles can come from one engine. The trace CSVs remain untouched and
`atom_scenarios::standard_core` still loads them - nothing is swapped.

---

## The through-charge limb (2026-07-24, later)

Built to answer the one prediction `report_magnetic_neutrality` found REFUTED:
neutron.pdf's balanced-field magnetic neutrality. That refutation's diagnosis was
architectural - every exit the pump offers is monochiral disc re-emission, so at
balance there is nothing mixed for the counter-stream to cancel against.

### The limb

In `march`, at the `r < R_IN` branch: an AMBIENT photon reaching the surface inside
the polar lane, travelling near-axially and inward, is NOT absorbed. It crosses to
the far side of the body along its own straight line (exit point from
t = -2(pos.dir), so no teleport) and continues, RETAINING ITS CHIRALITY AND SPIN.
That retention is the entire point. venus2.pdf: through charge "goes straight through
the body from pole to pole, avoiding lateral recycling". neutron.pdf: photons in the
south, antiphotons in the north, "then they go out the other pole".

Selection is geometric, per salt.pdf ("the charge that enters nearest the center of
the hole or pole [is] most likely to pass through. Charge that comes in nearer the
edges, or that enters on an angle, will be forced by centrifugal forces into the
equatorial whirlpool"):
- `THROUGH_LANE_RHO` = 1/3 of body radius. NOT tuned - it is the MEASURED hole from
  `report_axial_channel_geometry` (proton hole 127.5 against a disc reaching 382.5,
  ~11% of disc area). On a unit sphere that subtends asin(1/3) = 19.5 deg, landing on
  the same 10-20 deg axial cone already wired into `apply_photon` - two independent
  routes to one aperture.
- `THROUGH_COS_MAX_ANGLE` = cos(20 deg) axial-incidence gate, matching apply_photon's
  COS_EDGE.
- CLASS-INDEPENDENT on purpose (the same donut/hole map was measured for proton and
  neutron), so any p/n asymmetry has to emerge.

Consumes no RNG draws, so `through_lane_rho` = 0 is a strict no-op.

### report_through_charge (ignored, ~260s)

FIRST RESULT WAS A STARVATION ARTIFACT, and the report now says so in-source. At the
full ambient ensemble the limb fired 0-1 times per 40k launches. That is geometry, not
physics: the lane is b < 1/3 against an ensemble spread over b <= B_MAX = 6, so
P(b < lane) = 3.09e-3, and independently P(|dir.z| > cos20) = 6.03e-2, giving a joint
1.86e-4 - about 7 expected per 40k. Any verdict read off those rows would have
repeated the v1 lesson (never trust a transport null until the report proves its own
optical depth).

THE INTERPRETABLE RATE, conditioned on ARRIVAL rather than on the b<=6 launch
ensemble: of charge that actually reaches the body, 0.67% is on a through-trajectory
(11.1% by impact parameter, times the 6.0% axial gate).

LANE-ONLY PROBE (the conditional physics at real statistics). Launches only on
through-trajectories - b area-weighted inside the lane, direction axial with tilt
< 20 deg, entry pole set by species per neutron.pdf - marching through a converged
gas. ~500 of 20k survive the collisions and cross. Through-exhaust spin per escaped
photon:

| class | balanced | earth-2/3 |
|---|---|---|
| proton  | +0.271 | +0.451 |
| neutron | -0.149 | +0.078 |

against the DISC limb's monochiral +-1.0 to 1.6. So the through limb does supply a
LOW-SPIN output channel, and the neutron's gearing perturbs it least - 1.8x closer to
zero than the proton at balance, 5.8x at Earth, emergent from a class-independent
aperture.

CAVEAT KEPT IN THE OUTPUT: at a balanced field the probe population is launched 50/50
by species, so its net spin is ~0 BEFORE reaching the body. A near-zero exhaust there
does NOT show that something was cancelled - it shows how much the particle's gearing
BREAKS an already-balanced mixture. Only the comparative claim is defensible.

AXIAL CHANNELLING SWEEP (N raised to 250k - at 40k the expected count is ~7, far too
few for a trend). venus2.pdf's "the nucleus actually channels charge along the
axis... it is pushed there by charge streams" = `pole_intake_bias`:

| bias | proton (% of arrivals) | neutron (% of arrivals) |
|---|---|---|
| 0.5 | 0.049% | 0.168% |
| 2.0 | 0.494% | 0.885% |
| 8.0 | 4.632% | 4.732% |

Monotonic, roughly quadratic in the bias, ~100x over the swept range - and it
SATURATES near ~4.7% for both classes, so channelling opens the lane but does not
itself create the p/n asymmetry.

### VERDICT: the limb is built and correct, but it cannot be the neutrality mechanism
### for an ISOLATED particle

Total magnetic output barely moves with the limb on (M_tot/esc 1.3095 -> 1.3114), because
the disc limb still carries ~95-99% of the output even under strong channelling. So
balanced-field neutrality is still NOT delivered, and the reason is now quantified rather
than architectural.

This is consistent with venus2.pdf rather than in conflict with it: through charge is
"normally a minor complication" and dominates specifically at the NUCLEAR level, where
"much more charge passes through" because "the distance from pole to pole is so much
shorter" AND the nucleus channels charge axially. In a nucleus the nucleons are stacked
pole-to-pole and the axial lane is fed by the NEIGHBOURS' exhaust, not by isotropic
ambient - a multi-body configuration. Testing whether through charge delivers neutrality
therefore needs an axial STACK, which is Atom mode, not the single-particle cell.

---

## CM-4: the axial stack (2026-07-24, later still)

The single-particle line ended by pointing at a multi-body test, so this is it: a
second body on the pole axis, feeding the through-lane with its own exhaust
instead of relying on isotropic ambient. `report_axial_stack` (~70s, 5 stack
cells + 2 baselines).

### The seating is not a free parameter

graphene.pdf fixes it: "the proton is plugged in with its equator pointing down.
But the neutron is plugged in with its pole pointing down. This is because
protons channel charge pole to equator, while neutrons channel pole to pole."
Atom mode already asserts exactly this for plugs (`atom_scenarios`: plug protons
edge-on `rest_axis·Y ≈ 0` "disc feeds the hole", plug neutrons `|rest_axis·Y| ≈ 1`
"pole on the stack axis"), so CM-4 and Atom mode agree on nuclear geometry by
construction rather than by luck. `stack_seating_matches_atom_mode_plug_convention`
locks the two together.

Consequence worth stating plainly: a nuclear stack aims a proton's EQUATORIAL
DISC - the limb carrying 95-99% of its output - straight into the neutron's
polar lane.

### Architecture

`march_stack` is a separate path from `march`, so every single-body bit-identity
anchor is untouched. Three differences forced by having more than one centre:
absorption and the lane gate are tested per body in that body's OWN frame; the
chirality gears are evaluated with the photon rotated into the field body's local
frame and the new direction rotated back (skipping this would give an edge-on
proton a disc lying in the stack's xy plane instead of its own); escape is
measured from the stack centroid.

Phase-1 simplification, and it matters quantitatively below: only the feeder has
a converged gas. The probe is a passive absorber, so its own field never scatters
arriving charge back out of its lane, which makes every on-lane rate here an
UPPER bound. Body-level shadowing IS included, since real photons are marched
through the real geometry.

### CONFIRMED: a neighbour turns through-charge into a percent-level limb

At `sep = NUCLEON_PITCH = 2.6`, nuclear seating, per arrival at the probe:

| channel | arrivals | crossings | on-lane %/arrival |
|---|---|---|---|
| feeder's emitted exhaust | 1970 | 390 | 19.797% |
| ambient (same cell) | 3181 | 26 | 0.817% |

**24x**, and that is the only defensible enhancement figure: one field
configuration, one geometry, two launch populations. The tempting 204x against
the solo-probe baseline (0.097%) is NOT defensible - in the solo cell the probe
owns the gas and scatters its own arrivals out of the lane, in the stack cells it
owns none, and that difference alone accounts for 0.097% -> 0.817%.

In absolute terms the probe intercepts 4.92% of the feeder's total output and
passes 0.97% of the feeder's WHOLE exhaust through its lane, against ~0.001% of
launched ambient in the isolated cell. That is the number that decides whether
the limb can ever offset the disc.

Separation sweep, and it converges to a computable ceiling:

| sep | feeder reach | feeder %/arrival |
|---|---|---|
| 2.2 | 6.20% | 22.048% |
| 2.6 | 4.92% | 19.797% |
| 5.2 | 1.42% | 10.229% |

For a collimated beam the lane is a pure area gate, `(lane/R_IN)^2` = 11.1% of
the cross-section, so 11.1% is the far-field ceiling - which sep 5.2 sits on. A
NEAR feeder BEATS it by illuminating the polar cap preferentially (1/r^2 alone
weights the near pole 2.25x over the rim of the illuminated cap at sep 2.6). So
venus2.pdf's "the distance from pole to pole is so much shorter" clause is here
as a measured geometric effect, not an assertion.

### CONFIRMED, and emergent: the seating asymmetry

The lane gate is class-INDEPENDENT - nothing was told to prefer neutrons - yet:

- edge-on PROTON probe: **0 crossings from 1882 arrivals** (0.000%), against the
  pole-on neutron's 19.797%. graphene.pdf's "protons channel pole to equator,
  neutrons pole to pole" falls out of the seating.
- pole-on feeder control (what a stack test would assume WITHOUT graphene.pdf):
  reaches the probe with only 0.39% of its output vs the edge-on 4.92%, **13x**
  less, because a proton's poles are nearly dark. The paper's seating is what
  couples the stack at all.

`stack_lane_gate_is_evaluated_in_each_bodys_own_frame` pins this in a unit test:
one axial photon, pole-on body passes it, edge-on body absorbs it.

### REFUTED: neutrality via a p -> n feed. It goes BACKWARDS

The well-fed lane is monochiral (species mix 0.031) because a proton's exhaust is
all one species by construction, and its exhaust spin per escaped photon is
**1.437** - ABOVE the isolated limb's 0.27-0.45, because a crossing photon keeps
its spin and then stacks more of the same sign in the feeder's monochiral gas.

So the neighbour feed fixes the RATE problem and makes the SPECIES problem worse.
Two mechanisms have now failed for the same reason, which is the useful part: the
lane needs OPPOSITE-species feeds, and a single proton neighbour cannot supply
them.

### Where that points

ammon.pdf on a real nucleus: "charge moving pole-to-pole through the alphas,
south to north and north to south. The south to north channel is stronger, but
both exist, giving us both charge and anticharge." In this model the neutron's
own emission IS the anticharge (`chirality_sign` -1) while the protons' is
charge. So the candidate cancellation is between the neutron's through-charge (+)
and its own emission (-) - a WHOLE-PARTICLE sum, not a within-limb one, which is
also what voyag.pdf's "the spins will offset as a sum" actually says.

That is a 3-body p-n-p cell and a different observable (total magnetic output of
the middle body), not a bigger version of this measurement.

### Count discipline

Isotropic launches put only `(R_IN/B_MAX)^2` = 2.8% of photons on the body and
~0.1% of those on the lane, so 40k launches yield ~7 crossings. The isotropic
channels therefore run at 600k (baseline) and 400k (in-cell ambient), and the
report prints the expected count and the Poisson error on every small-count row -
the baseline's 3 and 4 crossings are "order 0.1%", nothing finer. This is the
same trap the through-charge report fell into on its first run.

One correction to earlier notes: the 0.67% isolated figure quoted in the
through-charge section is the COLLISION-FREE GEOMETRIC BOUND. The measured
isolated rate at these gear settings is 0.049-0.168%, and part [1] here
independently gives ~0.1%. The ~7x gap is scattering loss.

Suite 157/0/48 (3 new unit tests: seating matches Atom mode, lane gate is
per-body-frame, stack-of-one matches single-body geometry and two bodies shadow).

## CM-4b: the p-n-p whole-particle sum (2026-07-26)

`report_pnp_stack` (ignored, ~22s, 5 configs x 2 converged gases). The question
CM-4 could not ask: does the middle neutron's TOTAL magnetic output - its own
emission (chirality -1, the anticharge in this model) plus the through-exhaust
its proton neighbours feed it (+1) - sum toward zero, per ammon.pdf ("both
charge and anticharge") and voyag.pdf ("the spins will offset as a sum")?

### Results

- **CONFIRMED, the sign**: through-exhaust escapes at spin +1.49 to +1.52 per
  photon against the neutron's own -1.07. Species opposition is real and has
  the right sign to cancel. The whole-particle sum moves TOWARD zero when fed
  (p-n shift +0.029 on a -1.072 base).
- **REFUTED, the magnitude at bare-pair rates**: the shift is ~2%. Linear
  extrapolation says zeroing the sum needs f* ~ 379 proton-feeder equivalents
  = ~190x the p-n-p feed. Neutrality is NOT delivered by a bare pitch-2.6
  triplet at Phase-1 rates. The 190x is the number the nuclear papers' dense
  alpha-stack field must supply - the Phase-2 question, stated quantitatively.
- **EMERGENT - charge hand-off**: the far feeder sits on the exit axis and
  absorbs ~79% of the middle body's through-exhaust (p-n-p: 1379 crossings,
  283 stack-escapes vs p-n's 713/704). Charge crossing the middle body is
  mostly handed to the next body in line - the nuclear channel picture, found
  in the data, and the reason the escaped-sum dose-response is non-monotonic
  in feeder count.
- **Momentum (bodies static; the force a dynamic run would feel)**: one-sided
  p-n pushes the middle body off the feeder at +0.88 z-momentum per absorbed
  photon (3096 absorbed per 80k feeder budget); symmetric p-n-p cancels to
  +0.0003/absorbed while each side still bears the load. Absorber-force-only
  per cc.pdf (see project_momentum_open_field). New `StackOutcome.absorbed_by`
  carries the identity; unit test `stack_absorber_identity_for_momentum_tally`.
- **Counter-stream optical depth (what Phase 1 does not march)**: a lane
  photon crossing the north gap toward the far feeder runs tau = 0.755 through
  that feeder's disc gas - only 47% would survive uncollided. South gap: 7%
  loss (far feeder is distant), middle body's own gas: ~2% per gap. So the
  counter-stream collision channel is a factor-2 effect on inter-body hand-off
  and feeds the equatorial disc (gears funnel augments equatorial) - the
  alpha-densification mechanism, quantified but not yet marched.

### Caveats

- The |shift| column carries ~0.01 of non-lane variation (lane-shut null moved
  +0.012 with zero crossings: per-config seeds + shadowing). The clean signal
  is the thru column (704 escapes at +1.49); totals are illustrative.
- The n-n-n same-species control STARVED rather than failed: pole-on neutron
  feeders barely couple (reach 288 vs 3764 per feeder), so the species
  conclusion rests on the thru-sign opposition, which needs no control.
- Equal per-body emission budgets assumed when combining populations (19x
  mass/sec scales with mass, p ~ n).

Suite 158/0/49.
