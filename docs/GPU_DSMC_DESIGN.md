# GPU DSMC Photon Engine — Design (Phase 2, step 3)

Status: design, 2026-09-27. Builds on the wgpu harness (`rust/src/gpu/`, commit 9ae82b3).
Owner decisions still open are marked **[DECISION]**; the engine implements every
option so the choice is a config flag, not a rewrite.

---

## 0. Glossary

**Simulation method**
- **DSMC:** Direct Simulation Monte Carlo. Keep many real particles in flight at
  once; each time step, move them, sort them into grid cells, and collide randomly
  chosen pairs *within the same cell* at the rate kinetic theory gives. The standard
  method for rarefied gases (Bird 1994).
- **Monte Carlo:** answering a question by sampling many random cases and counting
  outcomes.
- **Phase 1:** the existing CPU model in `recycling.rs`. It marches *one* photon at a
  time through a precomputed, frozen density field, then re-tallies the field and
  repeats until it stops changing.
- **Phantom partner:** in Phase 1, when a photon "collides", its partner's direction is
  invented from the local field's average. The partner doesn't exist, so it never
  recoils.
- **Cell:** one cube of the collision grid (0.25 body radii on a side). Only photons in
  the same cell can collide.
- **NTC:** No-Time-Counter, Bird's recipe for how many candidate pairs to test per
  cell per step so the collision rate comes out right.
- **Cross-section σ:** the effective target area a photon presents; collision rate =
  σ × density × closing speed. *Not the same σ as "3σ" below.*
- **Mean free path λ:** average distance a photon travels between collisions,
  = 1/(σ·n). The engine's main physical knob.
- **TAU (τ), optical depth:** the expected number of collisions along a given path.
  τ = 1 means a photon crossing that path collides once on average; the survival
  fraction is e^(−τ).
- **Reciprocity:** one collision involves two photons, so the rate at which A-photons
  hit B-photons must equal the rate at which B-photons get hit by A-photons.
- **Renormalise:** rescale a direction vector back to length 1, i.e. speed c.
- **Steady state:** the flow has settled, so inflow equals outflow and averages stop
  drifting.
- **Kernel (GPU):** one small program the GPU runs on thousands of items in parallel.
  "Rate kernel" (D1) instead means the formula for the collision rate.
- **Seed / salt:** the number a random-number stream starts from. Same seed = identical
  run. `CFM_SEED_SALT` gives a different, equally valid set of random draws.

**Validation**
- **Gate:** a pass/fail check the engine must pass before its physics numbers are
  trusted. It is named for the gate you pass through, not a logic gate.
  - **G0 transport:** collisions off; photons must fly, hit, cross and escape
    exactly as the CPU Phase-1 code predicts.
  - **G1 kinetics:** in a plain box of photons, the measured collision rate must match
    textbook kinetic theory.
  - **G2 single body, G3 p-n-p stack:** report-only comparisons with Phase 1, which
    are *expected* to differ once the partner is real.
- **z / "3σ" (statistics):** how many standard deviations apart two measurements are.
  Noise alone rarely exceeds 3–4; a real disagreement grows with more samples.
- **Family-wise threshold:** when you make ~100 comparisons, a few will hit 3σ by
  luck, so the gate uses 4σ.

**Mathis mechanics (as modelled)**
- **Chirality χ:** a photon's spin sense, +1 (photon) or −1 (antiphoton). A frozen
  *tag*: it says which species the photon is.
- **Spin (carried):** the photon's magnetic degree of freedom. It starts equal to χ and
  is stacked or cancelled by collisions. The magnetic observables sum this, not χ.
- **Species:** photon (χ = +1) vs antiphoton (χ = −1). Near Earth the ambient field is
  2/3 photon (`photon_fraction`).
- **Ambient vs emitted:** ambient = the background charge field flowing in from
  outside; emitted = charge a body puts out (its recycled exhaust).
- **Gears / gearing:** Mathis's picture of two spinning photons edge-hitting like cog
  wheels. Whether the cogs *mesh* decides the outcome:
  - **Augment** (they mesh): spins stack; the charge escapes normally, and is nudged
    toward the body's equatorial disc.
  - **Cancel** (they catch): spins cancel; the charge is blocked, and is nudged toward
    the nearer pole (the through-channel).
- **Mesh sign:** the number that decides augment (+) or cancel (−). Phase 1 uses
  χ₁·χ₂; D2 is about whether travel direction also enters.
- **Redirect (`cancel_redirect`):** how strongly a cancel or augment nudges the
  photon's direction toward the pole or the disc (0 = none, 1 = fully).
- **Spin transfer (`spin_transfer`):** how much of the partner's spin a photon picks up
  per collision.
- **Glancing hit:** photons hit edge to edge, almost never head-on (arp.pdf). So
  collisions mostly change direction slightly and change spin, rather than bouncing
  photons back.
- **Lane / through-lane / through-charge:** the narrow polar hole (radius 1/3 of the
  body) where charge arriving nearly along the spin axis passes straight through the
  body pole-to-pole, keeping its species and spin, instead of being absorbed.
- **Disc:** a body's equatorial emission, where ~95–99% of its output leaves.
- **Seating:** how a body sits in a stack. **PoleOn** = spin axis along the stack (the
  lane faces its neighbours; neutrons). **EdgeOn** = spin axis sideways (the disc
  fires along the stack; protons).
- **Feeder / probe / middle:** in the p-n-p stack, the two outer protons are *feeders*
  whose disc exhaust feeds the *middle* neutron. *Probe* = whichever body is being
  measured.
- **Hand-off:** charge that crosses the middle body's lane and is then absorbed by the
  far feeder, i.e. passed along the stack. It measured ~78% in Phase 1.
- **f\*:** how many times more feed the middle neutron needs for its escaping spin sum
  to reach zero (neutral). Phase 1: ~170–180× the p-n-p feed.
- **Densification:** gas piling up in and around a stack so the feed gets stronger.
  This is the mechanism that would have to supply f\*.
- **Pump mode vs forced mode:** pump = each body emits exactly what it absorbs, so it
  is an engine and intake drives output. Forced = fixed emission rate (a Phase-1
  comparison mode).
- **Intake tilt (`pole_intake_bias`):** ambient photons are steered along the pole
  axis by species (photons enter south, antiphotons north; venus2.pdf).
- **Near zone / B_MAX = 6:** the region within 6 body radii, where Phase 1 does its
  detailed bookkeeping.
- **R_D = 12:** radius of the DSMC domain. Photons crossing it outward have escaped.
- **Ricochet cap:** Phase 1's maximum collisions per photon (30).
- **Swing boost:** a small sideways kick given to emitted photons by the body's outer
  spin (0.05 c).

## 1. Why, and what Phase 1 cannot do

Phase 1 (`recycling.rs`) marches one test photon at a time through a *frozen*
collision field (an analytic density × a tallied correction), iterating the field
to self-consistency. Three limits came out of the July campaign:

1. **The partner is a phantom.** `scatter_off` samples a partner direction from the
   cell's mean direction + coherence. Only the marching photon changes; the partner is
   never updated, so no collision pushes back on the population it came from.
2. **One body owns the gas.** In `march_stack`, only `field_body` has a gas; the probe
   is a passive absorber, and the far feeder's gas is ignored. The counter-stream optical
   depth (τ ≈ 0.76 north gap) was *integrated* but never *marched* (CM-4b part [5]).
3. **Densities can't respond.** Nuclear densification — the 170–180× feed f* that
   neutrality needs — has to *emerge* from the gas piling up, and a frozen field can't
   pile up.

**DSMC** (Direct Simulation Monte Carlo, Bird 1994) fixes all three at once. It keeps
millions of photons in flight simultaneously, sorts them into grid cells each step, and
collides randomly-chosen *real pairs* within each cell at the kinetic rate. Every body
emits and absorbs, every photon is someone's partner, and densities are whatever
the flow makes them.

## 2. What carries over unchanged (the sourced rules)

All of these are ported from `recycling.rs` as-is, same constants, same sources:

| Rule | Phase-1 source | Notes |
|---|---|---|
| Emission sampler | `spawn_emitted` | Bimodal 7°/58° exit latitudes (electron single-peak), same-hemisphere sign (deviation #1), swing boost × class sign, r = 1.001, body frame → stack via `StackBody::to_stack` |
| Ambient species + intake tilt | `spawn_ambient` | χ = +1 at `photon_fraction`; direction tilted by `pole_intake_bias · χ · ẑ` (venus2/neutron) |
| Through-lane gate | `march_stack` | ρ < `THROUGH_LANE_RHO`·r, \|dir·axis\| > cos 20°, heading inward, in *that body's* frame; straight-line crossing, species and spin retained |
| Absorption | `march`/`march_stack` | any surface hit that fails the lane gate |
| Contact normal | `sample_contact_normal` | glancing-biased (arp.pdf), u = 1−(1−u₀)³ |
| Gear redirect | `scatter_off` | cancel → lerp toward nearer pole by `cancel_redirect`; augment → lerp toward the equatorial plane; in a body's local frame |
| Spin transfer | `scatter_off` | Δspin = `spin_transfer` × partner chirality, clamp ±2 |
| Speed | — | \|dir\| = c = 1 always; collisions change direction and spin only ("photon collisions affect only the spins, not c", neutron.pdf) |

The lane rule generalises from the two Phase-1 versions: a photon may cross body *i*'s
lane **unless it was emitted by body *i***. That reproduces single-body (only ambient
crosses) and stack (foreign exhaust crosses the probe) exactly.

## 3. What changes because the partner is real

### 3.1 Both photons update
One pair collision updates both members. Exchange: sample the contact normal n̂ about
the relative direction, swap normal components
(`d1' = d1 − (d1·n̂)n̂ + (d2·n̂)n̂`, and symmetrically for d2), renormalise both to
c, then apply gear redirect and spin transfer to each, using its partner's chirality.
Renormalising keeps |v| = c (a sourced constraint) at the cost of exact linear-momentum
conservation per collision. The drift is measured and reported (§9, G1).

### 3.2 One cross-section, and why Phase 1 can't be reproduced exactly
Phase 1 calibrates **two** couplings from one TAU: σ_e (ambient photons scattering off
emitted density) and σ_a (emitted scattering off ambient). In a pair engine an e–a
collision is *one* event that hits both photons, so the two rates are tied by
reciprocity (rate_e·n_a = rate_a·n_e). Phase 1's scheme is non-reciprocal: no
pair-collision engine can match it at finite TAU. This is a known, expected difference,
not a bug to tune away.

The engine's physical knob is the **far-field ambient mean free path λ** (body radii):
σ = 1/(λ · n_a,sim). Optional per-pair-type multipliers (ee, ea, aa) exist for
diagnostics only.

### 3.3 Steady state instead of iteration
The simulation time-steps a continuous flow. Ambient photons enter through the domain
boundary; bodies emit; photons leave by absorption or escape. After a warm-up of about
2 domain-crossing times, tallies are time-averaged. Self-consistency is automatic.

- **Ambient injection:** a uniform isotropic field of density n_a has an inward flux of
  n_a·c/4 per unit area through any sphere. Each step, inject
  `n_a · c · π R_D² · Δt` photons (fractional remainder carried) at uniform points on
  the boundary sphere, with cosine-weighted inward directions (then the species tilt).
  With a nonzero tilt the field is deliberately non-isotropic, as in Phase 1.
- **Emission, forced mode (Phase-1 analogue):** each body emits at the rate its near
  zone receives ambient in Phase 1, where equal photon counts meant emitted rate =
  the rate of lines with b ≤ B_MAX = `n_a · c · π B_MAX²`.
- **Emission, pump mode (v4 analogue):** each body's emission rate = its own
  absorption rate, smoothed as an exponential moving average over 50 steps (intake
  determines output; photon number conserved by the surface).

## 4. **[DECISION]** Physics choices the engine exposes

The first two are the big ones. **D2 in particular could change the neutron result.**

**D1 — Collision-rate kernel.** `Constant`: rate ∝ σ n c, blind to relative motion
(Phase 1). `Catch`: rate ∝ σ n g, with g = |d₁ − d₂| ∈ [0, 2] — the standard kinetic
rate, where co-moving photons almost never meet and head-on ones meet most. `Catch` is
the same "catch" factor `apply_photon` already uses for photon-on-particle hits
(calibration.rs, CM-2). Mathis (bright.pdf): collisions concentrate where exiting
streams meet incoming ones. *Recommend `Catch` for physics runs, `Constant` for
the gates.*

**D2 — Gear mesh sign.** Phase 1: mesh = χ₁·χ₂ (direction-blind). Mathis says the
relative direction flips it:
- grav4.pdf: *"If a photon and antiphoton are moving in the same linear direction and
  edge-hit, they spin one another down. If they meet head-to-head, they spin one
  another up."*
- halbach.pdf: *"If opposing spins (photon and antiphoton) meet from opposite
  directions, they are opposite twice, and so are the same… If same spins (photon and
  photon) meet from opposite directions, they are opposite once, so they spin one
  another down."*
- bright.pdf: photons moving the opposite way "will look like right spinners."

`DirRel`: mesh = χ₁·χ₂·sign(d₁·d₂), and the spin picked up = k·χ₂·sign(d₁·d₂). That
is the photon-photon twin of `DirRelSign`, which `apply_photon` *already* uses for
photon-particle hits (commits 663d8d8/6301db2). So today the codebase is inconsistent:
direction-aware for particle hits, direction-blind for photon-photon.

Why it matters: a neutron's exhaust (χ = −1) leaving a pole *head-on* into incoming
ambient photons (χ = +1) is a **cancel** under Phase 1 but an **augment** under
`DirRel`. The polar funnel and the magnetic-cancellation results were all measured
under the direction-blind rule. *Recommend: gates in `ChiralityOnly`, then run every
physics config under both and report side by side, before changing any Phase-1 code.*

**D3 — Whose frame the gears use.** Gear redirect ("nearer pole", "equatorial plane")
needs a spin axis. Phase 1 used the field body's. With every body owning gas, the
engine uses **the nearest body's frame**, precomputed per cell (the body whose vortex
dominates locally). Alternative: density-weighted blend. *Recommend nearest.*

**D4 — Ricochet cap.** Phase 1 caps collisions per photon (`n_ricochet` = 30). The
engine keeps a per-photon bounce counter; `cap = 0` disables it. A capped photon still
acts as a partner (its partner updates) but doesn't change itself. *Gates: 30;
physics: off (report how many photons would have hit it).*

**D5 — Far-zone ambient-ambient gate.** Phase 1 switches off aa collisions beyond
B_MAX (a uniform isotropic gas is invariant under self-scattering). The engine keeps it
as a flag, off by default: with real pairs, the far field *is* the physics and costs
nothing extra. *Gates: on; physics: off.*

**D6 — Exchange rule (added after B1).** The §3.1 swap-then-renormalise keeps
|v| = c, but G1 measured it creating or destroying **~0.30 photon-momenta per
collision** (mean |Σ drift| over both photons). That averages out in an isotropic box,
but not in directed streams, where it could fake the momentum-push observable.

A rule exists that conserves both speed c and momentum exactly. For two equal-speed
photons, the relative vector g = d₁ − d₂ is always ⟂ to the mean V = (d₁ + d₂)/2 (since
|d₁| = |d₂|). Rotating g about V by a scattering angle χ gives d₁,₂' = V ± g'/2, with
|d₁,₂'| = 1 exactly and d₁' + d₂' unchanged. The glancing bias (arp.pdf) becomes a
small-χ-biased angle distribution. `SwapRenorm` (Phase-1-like) vs `RotateRelative`.
*Recommend `RotateRelative` for physics runs: sourced constraints (speed c) and
conservation both hold, and no parameter is added.*

## 5. Geometry and discretisation

- **Bodies:** static unit spheres (R_IN = 1) at fixed centres with a seating
  (`PoleOn`/`EdgeOn`), class and chirality — `StackBody` as-is. Up to 4 bodies.
- **Domain:** sphere of radius R_D = 12 about the stack centroid. Phase 1 used 30, but
  its collisions are concentrated within ~6 radii, and a smaller sphere puts photons
  where the physics is. Escape = leaving R_D. This is a known difference, reported
  alongside results.
- **Cells:** uniform Cartesian grid, h = 0.25, over the domain's bounding cube
  (96³ = 884,736 cells). Per-cell data is precomputed once on CPU:
  - free volume V_c = h³ minus body overlap, by Monte Carlo sampling (cells fully inside
    a body are never used);
  - index of the nearest body (D3);
  - whether the cell is inside B_MAX of any body (D5);
  - the lane-gap region tag (§8).
- **Time step:** c·Δt = 0.05 (Phase 1's near-surface ds; ≤ h/5). Warm-up ≈ 2 crossing
  times ≈ 1000 steps; averaging ≥ 2000 steps.
- **Precision:** f32 absolute positions. Error ~1e-6 at |x| ≤ 12 is far below h; no
  cell-relative encoding needed at this scale.

## 6. Photon storage (48 bytes, 3 × vec4)

| Word | Contents |
|---|---|
| `a` (vec4) | pos.xyz, `id` (u32 bitcast): a unique id, used for deterministic pairing |
| `b` (vec4) | dir.xyz, spin (f32) |
| `c` (vec4<u32>) | `flags`, `bounces`, `birth_step`, spare |

`flags` bits: 0 alive · 1 kind (0 ambient, 1 emitted) · 2 species (χ: 0 → −1, 1 → +1)
· 3–4 source body · 5–8 crossed-lane mask (bit *i* = crossed body *i*) · 9–12
arrived mask · 13 collided-since-last-crossing · 14–31 spare.

A single binding caps at 2 GB, so the capacity is ≤ 44M photons. The prototype targets
10–20M. Two-buffer (SoA) splitting is deferred until needed.

## 7. One step (all GPU-resident, one submit)

1. **Inject + emit.** New photons are written at slots `[n_live, n_live + k)`: ambient
   at the boundary, emitted at each body's surface. They get unique ids from a global
   counter.
2. **Move** (one thread per slot). Advance by c·Δt. Segment–sphere test against each
   body, taking the first hit along the segment:
   - if the photon passes the lane gate, it crosses and carries on with the remaining
     step length;
   - otherwise it is absorbed and tallied, and becomes dead.
   Leaving R_D = escape: tallied, dead. Then compute `cell_of` (dead → bucket
   `num_cells`).
3. **Sort.** `CellSorter` over `num_cells + 1` buckets. It is generalised to take a
   precomputed `cell_of` from the move kernel instead of computing it from position.
   Dead photons land in the last bucket, so after gather the live photons are
   contiguous in `[0, n_live)` and no free list is needed. `n_live` = offset of the
   dead bucket.
4. **Gather.** Reorder photons into sorted order (ping-pong buffers). This keeps cells
   contiguous and memory-coherent for collide and for next step's sort.
5. **Collide** (one thread per cell).
   - Load the cell's photons; order them by `id` (small in-register or shared-memory
     sort; cells hold ~tens), so pairing doesn't depend on atomic scatter order and runs
     are bit-reproducible.
   - NTC candidate count: `M = ½·N(N−1)·σ·g_max·Δt / V_c + remainder_c`, carrying the
     fractional part in a per-cell buffer. g_max = 2 for `Catch`; for `Constant`, use
     g ≡ 1 and accept every candidate.
   - For each candidate: pick i ≠ j via the RNG; accept with probability g/g_max; apply
     §3.1 plus the D2/D3/D4/D5 rules; set the collided bit; tally.
6. **Tally readback.** Small per-step tally buffers are read back each step and
   accumulated in f64 on the CPU (§8). Density maps accumulate on GPU and are read at
   the end.

**RNG:** counter-based PCG hash of (run seed, step, photon id or cell index, draw
counter). Nothing depends on thread scheduling, so a run is deterministic for a given
seed, and `CFM_SEED_SALT` works as it does for Phase 1.

## 8. Observables (tallies)

These are per-step integer counts (u32 atomics) plus fixed-point i32 sums (scale 2¹²,
reset every step, so no overflow).

- **Per body:** absorbed count by (kind, species, source body, crossed-mask ≠ 0);
  absorbed momentum Σdir.xyz (open-field absorber force); emission count; lane
  crossings by species.
- **Escapes:** count and spin sum by source body × crossed-body mask × species × exit
  latitude band (in the source body's frame).
- **The CM-4b headline set** (for the p-n-p config), computed from the above:
  - own M/esc of the middle body;
  - through-exhaust M/esc;
  - crossings and stack escapes → the hand-off fraction;
  - momentum per absorbed;
  - the lane-shut null.
- **Counter-stream survival, measured not integrated:** of the photons that crossed
  the middle body's lane and then traversed the north gap region (tagged cells between
  the middle and far bodies), the fraction that reach the far body's surface with the
  collided bit still clear. Compare with Phase 1's integrated exp(−0.76) = 0.47.
- **Density maps:** time-averaged photon counts per cell by kind and species → n_e,
  n_a and the densification profile. That is the Phase-2 quantity: **does the stack
  pile up enough gas to close the 170–180× f\* gap?**
- **Conservation checks:** injected + emitted = absorbed + escaped + Δ(n_live), every
  step (exact integer identity); per-collision momentum drift from renormalisation.

## 9. Validation gates

The engine is not trusted for physics until G0 and G1 pass.

| Gate | Configuration | Pass criterion |
|---|---|---|
| **G0 transport** | σ = 0 (no collisions): single body; then the p-n-p stack | Must match a CPU Phase-1 march at sigma = 0 on the same geometry (emission exit-latitude histogram; ambient arrival-latitude histogram; per-body reach, lane-crossing and hand-off fractions) within 3σ statistical, every bin. Integer conservation identity exact. |
| **G1 kinetics** | Ambient-only box, no bodies, isotropic, D5 off | Measured collisions per photon per unit path = σ·n (`Constant`) and σ·n·⟨g⟩ with ⟨g⟩ = 4/3 (`Catch`, isotropic unit vectors), within 2%. Direction distribution stays isotropic. Report momentum drift per collision. |
| **G2 single body** | One proton, then one neutron, pump mode, `Constant` + `ChiralityOnly` + cap 30 + D5 on; λ chosen so TAU_in matches Phase-1 TAU = 3 | **Report-only.** Exit-band topology (proton disc vs neutron butterfly), escape fraction and spin/escaped vs the Phase-1 cell, with differences attributed (reciprocity §3.2, real partners §3.1). |
| **G3 p-n-p** | The CM-4b configurations (n alone, p-n, p-n-p, n-n-n, lane shut), same settings as G2 | **Report-only.** The headline table vs the e0dcae7 re-baseline ranges, including measured vs integrated counter-stream survival. This is where Phase 2's "probe owns gas, far feeder's gas live" first shows up. |

After G3: physics runs with D1 = `Catch`, D2 in both variants, cap off, D5 off.

## 10. Performance budget (RTX 3060, from the 9ae82b3 bench)

Per step at 20M photons (48 B):
- move: ~6 ms (bandwidth-bound; 48 B read + write);
- sort: ~25 ms (random-input worst case measured 45 ms per 20M; coherent input
  after gather should be faster);
- gather: ~6 ms;
- collide: ~10–20 ms (estimate).

That is about 50 ms per step, so 3000 steps ≈ 2.5 min per configuration. The five-config
G3 table takes about 13 min per seed. Multi-seed ranges are affordable overnight or
by halving N.

## 11. Build phases

- **B1 — engine core + G0 + G1.** `rust/src/gpu/dsmc/` (config, host orchestration,
  tallies), WGSL kernels (`dsmc_inject_emit`, `dsmc_move`, `dsmc_gather`,
  `dsmc_collide`, shared RNG/photon-layout snippets concatenated host-side, since WGSL
  has no `#include`), CellSorter generalised to take precomputed `cell_of`. Gate tests
  live under `recycling.rs`'s test module behind `#[cfg(all(test, feature = "gpu"))]`
  so they can call the private Phase-1 samplers and `march_stack` for the CPU reference.
- **B2 — observables + G2 + G3 reports** (ignored tests, CM-4b-style printout).
- **B3 — physics runs:** the D1/D2 grid on the p-n-p stack plus density maps. That
  answers whether densification closes f\*.

## 11b. B1 findings (2026-09-27)

- **Gates pass:** G0 single body, G0 p-n-p (isotropic and bias 0.5), G1 (`Constant`
  0.96%, `Catch` 0.37% rate error), and determinism.
- **Two bugs caught by G0:**
  1. The boundary injector tilted cosine-weighted boundary directions, which is only
     correct for an untilted field. It now uses the Phase-1 construction (direction
     first, then a uniform line over the ⟂ disk).
  2. The first G0 stack test let feeder *exhaust* inflate the middle body's
     lane-crossing tally against an ambient-only CPU reference. This produced a
     spurious 40σ "disagreement" that was nearly waved through as a geometry effect.
     Emission and ambient now run separately.
- **Finding:** at σ = 0 no ambient photon crosses the p-n-p middle lane in either engine.
  The feeders shadow it completely. Through-charge in the stack is entirely
  feeder-fed or collision-fed.
- **Gate statistics:** a per-bin 3σ rule over ~100 comparisons false-alarms on noise.
  The gate is now 4σ family-wise, with emission at 1M photons per side (a 4.2σ bin at
  20k fell to 1.7σ at 1M).
- **Lane rule:** a photon crosses body *i*'s lane unless it was emitted by body *i*.
  `march_stack` itself has no exclusion; the difference only matters for a body's own
  exhaust re-entering its own lane after scattering.
- **Collide RNG key** was a linear combination of (cell, step) that aliases; it is now
  hash-mixed.
- **B2 must fix first — per-cell cap.** `dsmc_collide` pairs only the first 64 photons
  (by id) while computing the collision count from the full population. In any cell
  over 64, those 64 absorb all of the cell's collisions. That is biased exactly where
  the densification physics lives.
  - Fix: one workgroup per cell, a shared-memory bitonic sort by id (≤ 1024), then the
    sequential candidate loop over the full sorted list.
  - Cells over 1024: a loud tally (refine h), never a silent cap.
  - Also D6 `RotateRelative`, pump mode exercised, the G2/G3 reports and density maps.

## 11c. B2a done (2026-09-27)

- **Collide rewrite:** one workgroup per cell, a shared-memory bitonic sort of up to
  2048 photons by id, and a candidate loop over the full cell. The pair draw is fixed
  (j = i+1 was double-weighted). Cells over 2048 print a warning every step, never
  silently.
- **G1-dense** (~300 photons per cell) passes. An aggregate rate check can't see the old
  bug (the old kernel counted candidates from the full cell, so its rate was right).
  The gate instead measures the fraction of collisions where both partners have an
  id-rank below 64: 0.045 measured vs (64/300)² ≈ 0.045 expected, where the old kernel
  gives exactly 1.
- **D6 `RotateRelative`:** the total box momentum is read back in f64 before and
  after 200 steps. Over ~530k collisions it changes by 1.0e-4 (f32 rounding) vs
  2.8e2 for `SwapRenorm`, and max ||d|−1| < 1e-5. The fixed-point drift tally has a
  resolution of only 1/4096 per collision, so it is report-only, never a gate.
- **Still to do (B2b/B2c):** pump-mode test, the full §8 observables (density maps,
  counter-stream survival), σ calibration, and the G2/G3 reports with the first B3
  table.

## 11d. B2b result + experiment queue (2026-09-27, dddcc12)

**Single-body λ sweep** (2.9M photons, 1 seed, physics defaults, `report_dsmc_steady_single_body`):
- At σ = 0, proton and neutron are identical (same emission geometry), with spin ±1.
- **Phase-1 gear rule (`ChiralityOnly`):** proton stays a disc (equatorial 0.50 → 0.64
  at λ = 1); neutron goes polar (0.81 at λ = 1) and its escaping spin swings −1 →
  +0.70 through zero. This reproduces the Phase-1 picture with real partners.
- **`DirRel`:** the proton goes polar too (0.52 polar, 0.11 equatorial at λ = 1); the
  neutron goes polar (0.49) but its escaping spin barely cancels (−0.89…−1.04).
  Both observational signatures are lost, which points to a missing ingredient
  rather than a finished verdict.
- TAU_in stays 0.006–0.1 at every λ: an isolated body cannot densify its own gas.

**E1 — emitted-photon stiffness (user hypothesis, 2026-09-27).** The disc is
squashed because emitted photons deflect as easily as ambient ones. But photons
leaving a 12-spin baryon should carry spin energy from the emission spin-up and push
through several field hits.
- **Sources:**
  - pause.html: the proton "resists the spins of the photons better because it is
    larger";
  - bottom.pdf: a stacked spin "will double its mass equivalence in any
    collision", so stiffness is quantised, m = 2^k;
  - arp.pdf caveat: larger photons "suffer more collisions" but each changes them
    less, "these two factors offset" (for energy).
- **Model:** a per-photon stiffness m in a generalised RotateRelative (unequal-mass,
  momentum-conserving, both photons still at speed c, so the stiff photon deflects
  less). Emitted photons start at m = 2^k and shed it through collisions.
  - Optional: σ ∝ size² to test the arp.pdf offset.
- **Run:** the λ sweep under DirRel, k ∈ {0,1,2,3}, with and without σ scaling.
- **Questions:** does the proton disc survive? Does the neutron's magnetic
  cancellation come back? Does the offset erase the effect?

**E1 result (824fb46): stiffness is not the missing ingredient.** Under DirRel it
only partly resists the polar squash (proton λ = 1 polar 0.52 → 0.41 at k0 = 3).
Even permanent, unshed stiffness leaves the equatorial share at 0.20 (vs 0.50
collision-free). arp.pdf size scaling sheds stiffness almost at once and makes it
worse. Neutron spin is unchanged.

**E2 — redirect ablation (the decisive one).** With `cancel_redirect = 0`, the exit
topology is **identical** in every case: ~0.14 polar / 0.37 mid / 0.49 equatorial.
This holds for proton and neutron, both D2 modes, intake tilt on or off, and
λ = 3 or 1. It is exactly the collision-free emission geometry.
- Real (glancing, momentum-conserving) collisions do not reshape where charge
  leaves.
- **Every topology effect so far comes from the gear-redirect device:** the proton disc
  sharpening, the neutron polar funnel, and the DirRel polar squash.
- That device borrows neutron.pdf's *internal* funnel and applies it to every
  photon-photon hit in the gas, keyed to the nearest body's axis. So Phase 1's
  "neutron butterfly emergent from chirality" is emergent from chirality *through an
  imposed redirect*.
- The magnetic channel does respond to collisions without the redirect:
  - ChiralityOnly neutron spin/esc −0.38 (λ = 3) → +0.97 (λ = 1);
  - DirRel neutron stays −1.0…−1.13, since head-on meets with incoming photons spin
    it *up* (grav4.pdf);
  - DirRel proton drops to 0.03 at λ = 1.
- **Open question for the user:** the redirect is the whole topology engine and needs
  a sourced replacement. The candidate is PHYSICS_REFERENCE §5's edge-hit rule:
  "same-spin deflects ALONG the spin direction; opposite-spin deflects AGAINST it",
  i.e. a sideways deflection set by each photon's spin *axis*, which the model does
  not yet carry.

**Then B2c:** the G2/G3 tables and the p-n-p stack, run under the Phase-1-equivalent
settings (known reference) and with the redirect off (no imposed funnel), in both
D2 modes.

**B2c result (G2/G3, λ = 3, 1 seed; every body emits and owns its gas).**

| p-n-p headline | S1 compat | S2 P1-gear | S3 DirRel | S4 DirRel, no redirect | Phase-1 (e0dcae7) |
|---|---|---|---|---|---|
| through M/esc | 1.34 | 1.48 | 0.91 | 0.69 | 1.46–1.51 |
| hand-off to far feeder | 77.2% | 71.5% | 74.1% | 72.3% | 77–78.5% |
| north-gap survival | 0.875 | 0.881 | 0.907 | 0.729 | ≈0.47 (integrated) |
| densification (emitted gas, 1 < r < 1.5, vs n alone) | 1.56× | 2.22× | 1.63× | 1.37× | needs ~170–180× (f\*) |
| total-sum shift, p-n-p vs n alone | 0.032 | 0.079 | 0.010 | 0.001 | — |

- **Hand-off is robust.** ~72–77% in every setting, and S1 lands on Phase 1. The far
  feeder capturing through-exhaust is geometric, not a gear effect.
- **Through-exhaust spin needs the old gear rule.** It stays high only under the
  Phase-1 gear rule; DirRel and the no-redirect setting cut it roughly in half.
- **Counter-stream survival** is ~0.88 here vs Phase-1's integrated 0.47. The
  difference is the reciprocity regime (§3.2): a pump-mode exhaust is optically thin,
  so Phase-1's TAU = 3 emitted gas does not exist here.
- **Densification is real but small.** The stack piles emitted gas up only 1.4–2.5×
  around the middle neutron, nowhere near the ~175× neutrality needs. The whole-particle
  spin sum barely moves (|shift| ≤ 0.08). In this model, stacking three nucleons does
  not deliver neutrality in any setting.
- **Caveats:**
  - `mom/abs(mid)` (0.08–0.19) averages over *all* photons the middle body absorbs
    (ambient, own and feeder). Phase-1's +0.88 was feeder photons only, so the two are
    NOT comparable until absorbed momentum is split by source.
  - p-n `surv(S)` = 0.000 under S2 is almost certainly 0/0 (no south-bound crossers).
  - The S1 lane-shut row failed its steady check (z_absorbed −4.1).
  - Single seed.


## 13. Next: photons with full mechanics (proposal, 2026-09-28)

User directive: photons get whatever state the mechanics need (spin axis, stacked
spins, precise energies); optimise later. E2 showed why: with photons reduced to
direction + scalar spin, the only thing shaping exit topology was the imposed redirect.

### 13.1 Photon state
- Position, velocity (see D7 for the speed), and **mass/size level k** (m = 2^k,
  radius ∝ 2^k; bottom.pdf).
- **Spin vector ω**: its axis is the spin axis, its sense is the chirality, and its
  edge speed is |ω|·r ≤ c.
- **Stacked levels:** when an edge speed would exceed c, the spin stacks a new level at
  2× radius instead (quant.html/half.pdf quantisation, PHYSICS_REFERENCE §2). The
  photon tier is levels 1–4 (axial, x, y, z).
- **Energy bookkeeping per level** (the §3b quantum-spin-energy formula), so
  spin ↔ linear transfer is accounted in real units.
- Roughly 96–128 B per photon. That is fine at a few million photons.

### 13.2 Contact mechanics: perfectly rough spheres
The textbook model for colliding *spinning* balls (Bryan 1894; Pidduck 1922; Chapman &
Cowling ch. 11).
- At contact, the relative surface velocity is v_rel = (v₁ − v₂) + ω₁×(r₁n̂) + ω₂×(r₂n̂).
  The impulse has a normal part (the bounce) and a tangential part (the surfaces
  grip).
- A tangential restitution β sets the grip: −1 = frictionless, +1 = perfectly rough,
  "gears that don't slip".
- Linear momentum, angular momentum and energy are conserved together, and spin ↔
  linear energy exchange happens on the tangent. This is halbach.pdf's "it happens on
  the tangent… energy has transferred from spin energy to linear".

**Why this is promising:**
- **D2 becomes derived instead of a rule.** Whether two surfaces *mesh* (tangential
  surface speeds match, little exchange) or *catch* (they oppose, hard grip) depends
  jointly on both spin senses and on relative direction. That is Mathis's "opposite
  twice = same", falling out of the geometry.
- **The sideways deflection comes for free.** "Same-spin deflects ALONG the spin
  direction; opposite-spin deflects AGAINST it" (PHYSICS_REFERENCE §5) is the
  tangential impulse; nothing needs to be imposed.
- **The redirect device can be retired**, and E2's question gets a real answer: does
  topology emerge from the mechanics?

### 13.3 **[DECISION] D7 — the speed of light in collisions**
Rough-sphere spin ↔ linear transfer changes a photon's *speed*. The sources pull
both ways:
- neutron.pdf ("collisions affect only the spins, not c") and arp.pdf ("spin changes,
  not velocity changes");
- photon.html ("c is the average speed, not the maximum… newly emitted photons may be
  going slightly over c") and halbach.pdf (repulsion is an *increase in linear
  energy*).

The options:
- **(a) Free speed:** full rough-sphere; speeds scatter about c.
- **(b) Speed pinned at c:** the linear part of the exchange is projected back to
  |v| = c, and the energy difference goes into spin (energy conserved, momentum
  approximately conserved).
- **(c) Implement both** and measure the speed spread (a) produces. If it stays
  within a few % of c, arp.pdf's picture holds either way.

*Recommend (c).*

### 13.4 Gates before physics
- **RS1:** a rough-sphere box conserves energy, linear and angular momentum.
- **RS2:** the textbook result that a perfectly rough gas equipartitions translational
  and rotational energy is reproduced (a sharp, well-known check).
- **RS3:** a two-photon collision table (same/opposite spin × co-moving/head-on)
  reproduces the sign pattern in grav4.pdf/halbach.pdf. That verifies D2 emerges and
  isn't assumed.
- **Then:** rerun E2 (single body, no redirect) and the G3 stack with the
  mechanics doing the work.

### 13.5 Built (2026-09-29): user chose D7 (c), both speed modes

**Sources pinned:**
- photon2.html: "let the tangential velocity of the spin of the photon equal c", so
  photons are born at edge speed c (`RoughParams::spin_edge0 = 1`), helical, with
  ω along +d for chirality +1 (clockwise from behind, quark.html right-hand rule).
- super.html (2013 update): a photon "stacks a second spin on top of the first because
  it cannot spin any faster on the first axis. It has reached a spin velocity of c",
  and the new spin is centred on the collision point. That is the stacking rule
  (§13.6).

**Contact law: two channels, not one.** A spinning-ball contact can grip in two ways:
- **Edge (sliding) grip `beta`:** the classical rough sphere of §13.2.
- **Twist (drilling) grip `beta_twist`:** the contact faces twist against each other
  about the line of centres. This is a pure couple, so spin is exchanged but linear
  motion is untouched. `twist_loss_to_normal` pays any energy it removes into normal
  separation (halbach.pdf: "repulsion implies a real loss of local spin… transferred
  from spin energy to linear").

**RS3 result (CPU, `rough_sphere::rs3_sign_patterns`).** The §13.2 claim that D2 falls
out of the textbook model is **wrong as stated**:
- **Twist channel:** reproduces grav4/halbach DirRel exactly. Same spins head-on catch,
  spin down and push apart; opposite spins head-on pass clean; the pattern flips for
  co-moving pairs.
- **Edge channel:** gives the **exact inverse**. Two same-helicity photons meeting
  head-on counter-rotate in the lab frame, and counter-rotating gears mesh.
- **Control:** end-over-end spin (no helicity) is chirality-blind.
- With both channels on, the edge channel dominates (about 2× larger).

So which way the gear sign goes depends on which grip photons have: an empirical
question. E3 runs all of them.

**D7 (b) had to be reformulated.** Pinning each photon's speed and booking the change
into its *own* spin ran out of spin energy:
- a κ = 0.4 sphere at edge c holds 0.2 spin vs 0.5 linear;
- 664k deficits in a 200k-photon box, i.e. energy created.

`PinnedC` is now pair-level: after the contact, both speeds go back to c and both spins
are scaled by one factor so the pair's spin energy equals what it came in with. That
matches neutron.pdf's "collisions affect only the spins, not c".
- Energy is exact by construction.
- Momentum is not: 0.9% of Σm|v| drift over ~19 collisions/photon.

**Gates (all pass, each with a control that fails):**

| Gate | Result | Control |
|---|---|---|
| RS1 CPU per collision | P, L, E exact to 1e-14 | partial grip β = 0 loses energy |
| GPU vs CPU, 5 law × D7 combos | ≤ 4e-6 per collision | wrong law (β = 0.5) differs by 0.7–5.5 |
| RS1 GPU box, Free | dE/E ~1e-6, dP ~1e-9 | β = 0 loses 82% |
| RS1 GPU box, PinnedC | dE/E 1.6e-7, all \|v\| = c to 2e-7 | — |
| RS2 equipartition | E_rot/E_lin 0.400 → 0.9996, ⟨v⁴⟩/⟨v²⟩² = 1.669 (Maxwell 5/3) | smooth spheres stay at 0.400 |

**Edge speeds exceed c constantly under edge grip.** At birth, E_rot/E_lin = κ = 0.4,
below equipartition, so edge grip spins photons *up*, and a large fraction of
post-collision edges exceed c. That is exactly super.html's stacking trigger, so
stacking is not optional under edge grip (§13.6).

## 12. Out of scope for this step

Moving or rotating bodies (atom-mode coupling), the Godot live view, >44M photons,
spatially varying `photon_fraction`, electron bodies in stacks.
