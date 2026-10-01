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

### 13.6 E3: E2 rerun with the mechanics doing the work (2026-09-29)

**Setup:** `report_e3_rough_single_body`. Single body, λ = 3, pump emission, no
redirect, no spin_transfer. Runs cover 4 contact laws × both D7 modes × both
impact distributions × proton/neutron, at 2400 warmup + 2000 averaging steps, 1 seed.
- "glanc" is the RotateRelative/arp.pdf grazing skew, i.e. E2's transport.
- "geom" is the true geometric hard-sphere impact distribution: near-isotropic
  scattering, diffusive at λ = 1.

| law (glancing, PinnedC) | polar | mid | equat | proton hel/esc | neutron hel/esc | esc edge>c |
|---|---|---|---|---|---|---|
| smooth (spin inert, null) | 0.138 | 0.366 | 0.497 | +0.31 | −0.31 | 0 |
| edge grip | 0.191 | 0.400 | 0.408 | +0.02 | −0.00 | 0.40 |
| twist grip → normal | 0.137 | 0.366 | 0.497 | **+0.40** | **−0.16** | 0.43 |
| edge + twist | 0.152 | 0.378 | 0.469 | +0.03 | +0.01 | 0.40 |

- **The null reproduces.** Smooth spheres give E2's collision-free geometry exactly.
  Helicity decays to ±0.31 even with inert spin, because directions scatter while ω
  doesn't (helicity = ω·d).
- **Topology is chirality-blind in every law.** Proton and neutron rows agree to
  within 0.003 in all 32 runs.
  - Edge grip squashes both toward the pole by 5 points, via spin → linear
    tangential kicks.
  - Twist grip touches no topology at all.
  - With emission geometry identical, **the disc vs polar-funnel distinction does not
    emerge from photon-photon mechanics at a single body.** It has to come from the
    emitter: the body's own spin geometry and channelling, i.e. where Mathis puts it.
- **The magnetic channel is where the mechanics act.**
  - Twist grip preserves the majority-chirality exhaust (proton +0.40 vs +0.31
    smooth) and spins the minority down (neutron −0.16 vs −0.31). That is the 2:1
    ambient field despinning the opposite species (antimatter-instability
    reference), and it emerges from DirRel-signed twist contacts.
  - Edge grip scrambles helicity to zero for both.
- **D7.** Free-mode speeds drift down under edge grip: escaping 0.79 ± 0.27 c
  (glancing), 0.86 ± 0.28 c (geometric). Spin soaks up linear energy on the way to
  equipartition, and Free runs never pass the steady check (the gas keeps cooling).
  - Photons do not stay "within a few % of c" under a free rough-sphere law. That
    favours (b) PinnedC, or a law that doesn't pump linear into spin.
  - Twist-only Free stays at 1.05 ± 0.34 c.
- **Edge speeds above c:** 40–67% of escaping photons under edge grip; 43% under
  pinned twist, because the restore step re-inflates the pair's spin. super.html's
  stacking trigger fires constantly, so §13.7 is needed before these numbers are
  final.
- The geometric impact distribution gives the same picture with weaker effects
  (escaping helicity is 0.20 even for smooth spheres).

**G3 p-n-p stack, rough-sphere settings** (`report_g3_pnp_stack`,
`CFM_G23_SETTINGS=S5,S6`, glancing, PinnedC, λ = 3, 1 seed):

| p-n-p headline | S5 twist | S6 edge+twist | B2c S1–S4 | Phase-1 |
|---|---|---|---|---|
| through M/esc | 0.34 | 0.02 | 0.69–1.48 | 1.46–1.51 |
| hand-off | 71.4% | 67.5% | 72–77% | 77–78.5% |
| north / south survival | 0.75 / 0.77 | 0.78 / 0.80 | ~0.73–0.91 | ≈0.47 |
| densification (1 < r < 1.5, vs n alone) | 1.37× | 1.41× | 1.4–2.2× | needs ~175× |
| total-sum shift, p-n-p vs n alone | 0.007 | 0.002 | ≤ 0.08 | — |

- Hand-off and densification are unchanged by the photon mechanics. They are
  geometric/transport properties of an optically thin pump-mode stack. **Real contact
  mechanics do not close the ~175× densification gap.**
- The through-lane feed is almost pure χ+: 618 : 11 in the p-n-p. That is the feeders'
  proton exhaust.
  - Under twist grip it keeps a third of its helicity (0.34).
  - Under edge grip it is scrambled (0.02).
- A lone neutron's own exhaust is despun by the field: −0.16 under S5, ≈0 under S6.

### 13.7 **[DECISION] S1 — spin stacking accounting (proposal)**
super.html: at edge speed c, extra spin energy stacks a new level centred on the
collision point. bottom.pdf: each stacked level doubles collision mass (m, R → 2m, 2R).
The proposed event, applied after the contact whenever the outer edge exceeds c:
- The old outer level locks at edge c into a per-photon inner-energy store. Inner
  levels don't collide (neut2.pdf "collisions affect only the outer spins").
- k → k + 1.
- The new outer spin's axis is the contact torque axis n×J.

The open question is where the cost comes from. At speed c, doubling m doubles the
linear energy ½mc². Options:
- **(i) Energy-exact:**
  - The lock energy ½I/R² and the linear cost ½m_k c² are paid from the spin excess.
    A k = 0 photon needs an edge of ≥ 1.87c to afford a stack.
  - Short of that, the edge stays above c (pinned mode) or the excess goes to speed
    (free mode).
  - Momentum is not conserved at the stacking event.
- **(ii) Cap-only:**
  - The edge is clamped at c and the excess banked as inner energy; no mass or size
    change.
  - It's simple, but ignores bottom.pdf.
- **(iii) Stack without a linear cost:** the doubled mass is a collision property
  only, and energy is tallied separately. This breaks the conservation gates.

Unstacking (limat.pdf, collisions "strip" high photons down) is the reverse event and
needs a threshold. *Recommend (i) with unstacking at outer spin → 0.*

**Decided 2026-09-29 (user): none of (i)–(iii). There is no doubling cost; energy is
conserved and levels are fractional.** Built as `RoughParams::stacking`
(`rough_sphere::Stacked::restack` is the CPU reference):
- **Levels:** 0–3 (axial, x, y, z), radius 2^k, each with its own spin vector. Inner
  levels live in the photon's e/f/g words, and the record is now 112 B.
- **Collision:** only the outer level collides (neut2.pdf). Linear mass is 1 at every
  level, the contact radius is 2^k, and I = κ·(2^k)².
- **Stack:** when a collision pushes the outer edge above c, that level locks at edge c
  (keeping its axis). The excess energy starts level k+1, with its axis along the
  collision's spin kick (super.html) and its edge speed at whatever fraction of c the
  excess buys. This can repeat, up to the z-level (3); beyond that the edge stays
  above c and is tallied as `capped`.
- **Unstack:** a stacked outer level spun below 0.01c is removed. Its energy returns
  to the level below, which becomes the outer level again (limat.pdf).
- **Property:** at edge c every level holds the same energy (½κ), because mass is 1
  and I ∝ R². So a stack takes new-edge² = old-edge² − 1. Mathis's per-level energy
  table (9, 65, 1025, 16385) is *not* built in; whether anything like it emerges is a
  measurement.
- **Consequence:** a stacked photon's outer axis is generally sideways to travel, so
  its collisions carry no handedness (panel A in the grip discussion). Handedness
  survives in the locked axial level, which `b.w` now reports.
- **Grip default:** twist grip for the axial level (user, 2026-09-29, after
  distinguishing lab-frame from travel-relative "same spin").
- **Gates:**
  - CPU: energy exact through stacking, multi-level stacks, the top cap and unstacking
    (1e-15); round trip back to the axial level.
  - GPU vs CPU with stacking: ≤ 1.4e-6 on every level. The stack axis is built from
    the impulse terms, not from the f32 difference of spins before and after, which
    had tilted it by 1e-3.
  - Stacking box: dE/E ≤ 5e-7, inner levels stay at edge c.
  - Early sign: Free + edge grip drives the box to the top level (57% at L3, 207k
    capped events) as spin fills toward equipartition. PinnedC edge sits mostly at
    L0–L1.

**E3 with stacking** (glancing, λ = 3, 1 seed, `CFM_E3_STACKING=1`):

| law, D7 | polar/mid/equat (p) | same for n | hel p / n | escape levels L0/L1/L2/L3 | stacks/photon |
|---|---|---|---|---|---|
| smooth (null) | 0.138/0.366/0.497 | same | +0.31 / −0.31 | 1/0/0/0 | 0 |
| twist, PinnedC | 0.137/0.364/0.499 | 0.138/0.365/0.497 | +0.32 / −0.25 | .48/.47/.05/.01 | 0.59 |
| twist, Free | 0.138/0.365/0.498 | 0.138/0.364/0.498 | +0.24 / −0.15 | .75/.25/0/0 | 0.25 |
| edge, PinnedC | 0.279/0.400/0.321 | 0.278/0.401/0.321 | +0.04 / −0.04 | .54/.38/.07/.01 | 0.56 |
| edge+twist, PinnedC | 0.215/0.403/0.382 | 0.217/0.400/0.383 | +0.03 / −0.03 | .57/.33/.08/.03 | 0.57 |
| edge or edge+twist, Free | **runaway** | | | | |

- **Stacking is common:** about half of escaping photons carry a second level under
  PinnedC. Almost none reach L3 on a single body at λ = 3, and unstacking is rare
  (1e-4 per photon).
- **Topology is still chirality-blind** (p = n within 0.003) under every law, with or
  without stacking. Stacking strengthens the edge-grip polar squash (0.19 → 0.28,
  equatorial 0.41 → 0.32), but for both species alike. That confirms §13.6: the
  disc/funnel split must come from the emitter.
- **Twist grip keeps its magnetic asymmetry** (+0.32 / −0.25 pinned). It is weaker than
  without stacking (+0.40 / −0.16), because stacked photons collide on sideways outer
  axes that carry no handedness.
- **Free speed + edge grip is a runaway.** The stacked spin reservoir drains linear
  energy, photons slow, and the population outgrows 6M within 30–40 s. Together with
  the no-stacking Free result (0.79c), **D7 favours PinnedC** whenever edge grip is
  on. Twist-only Free stays near c (1.03 ± 0.34).

## 14. Next: a mechanical baryon on the new engine (proposal, 2026-09-29)

**Why.** Across §13 the photon gas turned out to set the *magnetic* channel (helicity),
but never the *topology*. Proton and neutron exits are identical in every contact law
while the emitter is a fixed-profile absorb/re-emit sphere. Densification also needs a
real trap. Both point at the emitter. The current bodies are placeholders: they absorb
whatever touches them and re-emit from a hard-coded latitude histogram (7°/58°, the
Spin-mode trace).

**What Mathis says the body does:**
- **pole.pdf:** "any spinning sphere will have more angular momentum at the equator
  than the poles. Given an ambient field of smaller particles (the photons) and
  permeability of the larger body by the smaller, this will create a net vector out at
  the equator and a net vector in at the poles… the photons will feel a 'negative
  pressure' at the poles and will be sucked in." Also, the greater field "will sort
  the smaller particles, sending the photons to the south pole… and antiphotons to the
  north pole".
- **rhrule.pdf:** "protons act like fans (the photons don't just channel through the
  nucleus, they channel through the individual nucleons as well)". Pole-to-equator is
  the primary channel; pole-to-pole through-charge is the neutron's primary scheme.
- **galcharge.pdf / core2.pdf:** incoming photons "gain energy by stacking on more
  outer spins" as they come in at the poles. Through-charge meeting the opposite
  stream is spun *up* (opposite twice). §13.7 stacking and the twist grip already
  carry both.
- **graphene.pdf:** the nucleus is "spun up by the ambient charge field", and emission
  leaves "like a spinning lawn sprinkler".

So the body is **permeable, spinning matter that photons pass through and collide
with.** Emission is scattered ambient, not absorbed-and-reborn charge, and the
equatorial preference should come from the radial growth of surface speed (ω×r).
Nothing imposes it.

### 14.1 [DECISION] B1 — how much of the 12-level stack is modelled
- **(a) Uniform fan.** One permeable sphere rotating rigidly about its axis at edge c.
  Photons inside collide with "body matter" at velocity Ω×r. This is the minimal
  mechanism of pole.pdf. It has no chirality structure, so proton and neutron are
  identical by construction. It is a control.
- **(b) Layered shells (recommended).** The baryon-tier levels 9–12 (axial, x, y, z)
  become four nested permeable shells at radii ⅛, ¼, ½, 1. Each moves with that level's
  *instantaneous* angular velocity, taken from `spin_stack.rs` (the same sourced
  kinematics that drives Spin mode: hierarchy, orbit radius = amplitude/2, precession).
  Each has its own chirality sign from the proton/neutron patterns in
  PHYSICS_REFERENCE §3 (each class has four listed states, e.g. proton +a+x+y+z, neutron −a−x−y+z; all four get run). A photon inside shell k
  collides with matter moving at that level's local velocity and carrying that level's
  spin. That is where proton and neutron can finally differ mechanically: the
  neutron's reversed outer z-spin meets photons the inner levels have already
  channelled.
- **(c) Full trace.** The Spin-mode swept path (spirograph capsules) at full
  resolution, with photons hitting capsule segments. This is the most literal version,
  and the old 7°/58° histogram came from it. It is heavy and time-dependent, and the
  trace at levels 9–12 is enormous relative to a photon. It is best kept as a later
  validation of (b), not the workhorse.

### 14.2 [DECISION] B2 — what a photon-body contact is
Use the same rough-sphere contact (twist grip default, PinnedC). The partner is an
infinitely massive element of the shell, moving at the shell's local velocity and
spinning with the shell's level spin. There is no new physics rule: the contact law is
the one gated in §13.
- **Alternative:** a smooth "fan-blade" deflection (tangential velocity kick only). It
  is simpler but repeats the redirect mistake of imposing the outcome.

### 14.3 [DECISION] B3 — how permeable the body is
This means how many body contacts a photon suffers crossing a shell: an optical depth
τ_body per shell. Mathis gives no number, so it is a **swept knob**:
- τ_body → 0: transparent, no channelling.
- τ_body → ∞: an opaque ball that reflects.

The prediction to test is that a window of τ_body produces pole-in / equator-out with
proton ≠ neutron. We report how wide that window is, never a best fit.
- Possible later anchor: the charge-recycling rate (19× mass/s,
  reference_charge_recycling_rate) as a throughput constraint.

### 14.4 [DECISION] B4 — does the body feel the photons (back-reaction)
- **(a) Driven (recommended first):** level spins are held fixed (steady state, "spun
  up by the ambient field"). The net torque and force the photons exert are tallied
  per level. A consistent model should show the torque roughly balancing, which is a
  check rather than an input.
- **(b) Free:** levels gain and lose angular momentum from hits, and the body can
  drift. This is the free-drift focus particle from the backlog, and the natural step
  after (a).

### 14.5 What goes away, what stays
- **Retired:** pump mode, the absorb/re-emit emission histogram, the through-lane gate
  (`lane.x`, cos 20°). Through-charge must now *emerge* as photons crossing pole to
  pole, and the neutron's "primary scheme" is a prediction.
- **Kept:** the ambient injector, the DSMC photon gas (§13 mechanics + stacking), and
  every observable (exit bands, helicity, survival, density maps, absorbed momentum
  → now "momentum delivered to the body").
- **Gates first:**
  - **BG0:** a transparent body (τ = 0) leaves the ambient field untouched (flux and
    isotropy match G0).
  - **BG1:** in a non-rotating opaque body, spin and linear momentum given to the body
    balance the photons' change (conservation across the body boundary).
  - **BG2:** under the uniform fan (a), photon angular momentum about the axis is
    delivered in the rotation sense, and the equatorial exit excess grows with τ
    (pole.pdf's mechanism, measured).
  - **BG3:** mirror symmetry. Proton vs anti-proton (all signs flipped) in a mirrored
    field gives mirrored results.
- **First physics run:** single proton vs neutron under (b). Does topology split? Does
  a neutron pass pole-to-pole? What exit latitudes appear (compare 7°/58° and the
  30°-cone claim), and does the stack densify?

### 14.6 Decided (user, 2026-09-29) and the build plan
- **B1 = (b) layered shells** to start. The user added a planned detour: rebuild
  Spin-mode tracing (option c) on the new mechanics and rerun its histograms, for
  comparison.
- **B2 =** the same rough-sphere law. **B3 =** sweep τ_body. **B4 =** driven body
  first.

**Kinematics, reused exactly from `spin_stack.rs`:**
- Levels 9–12 have amplitudes 64 / 128 / 256 / 512 (`types::level_amplitude`),
  normalised so level 12 = the body radius 1. That gives ⅛, ¼, ½, 1.
- Level 9 is a *precession* (rotates in place about local Z, orbit radius = its
  amplitude). Levels 10–12 are orbital (X, Y, Z spins, orbit radius = amplitude/2).
  At c, the angular rates are 8, 8, 4 and 2 per unit time.
- **Shell regions are the Spin-mode ghost spheres.** Level k's sphere (radius
  amplitude_k) is centred where the levels outside k put its frame origin, so shells
  nest off-centre, each touching the next from inside. A photon belongs to the
  innermost shell containing it.
- **Matter velocity in shell k:** V(x) = Σ_{j≥k} ω_j × (x − o_j). Here ω_j is level
  j's world angular velocity (its local axis carried by the levels outside it), and
  o_j is its rotation centre (its frame origin). This is exact for nested rigid
  rotations, and a CPU gate checks it against finite differences of the composed
  motion.
- **Chirality:** each level's sign ±1 sets its rotation sense. The four proton and
  four neutron patterns (PHYSICS_REFERENCE §3) all get run.

**Body contact:**
- A photon inside the body meets body "grains" at rate τ_body·|v − V(x)|.
- The contact normal is sampled like photon-photon contacts.
- The grain is infinitely massive, moves at V(x), and has no spin of its own. Its
  spin at photon scale is negligible, because Ω·R_photon ≪ c.
- The contact law is the same `ContactLaw` (§13), with D7 PinnedC or Free.
- Consequence: level chirality enters **only through the direction of the matter
  flow V**. The photon's own spin couples to it through the grip.

**Engine changes:**
- `BodyModel::{Absorber (legacy), Shells{tau_body}}`.
- In Shells mode bodies are transparent regions with no absorption and no emission.
  All charge is scattered ambient.
- New per-body tallies:
  - entries;
  - exits after ≥1 contact;
  - entry and exit latitude histograms;
  - contacts per passage;
  - pole-to-pole passages;
  - helicity by exit band;
  - momentum and angular momentum delivered to the body (driven-body torque check).

**Gates:**
- the kinematics finite-difference gate;
- GPU contact vs CPU mirror;
- BG0 transparent;
- BG1 non-rotating body conserves photon energy, and delivered momentum balances;
- BG2 uniform fan gives angular momentum in the rotation sense;
- BG3 mirror symmetry.

### 14.7 First results (2026-09-29, fd00989 + report_b_shell_single_body)

**Setup:** single body, all 4 proton + 4 neutron sign states, λ = 3, PinnedC,
stacking on, 1600 averaging steps, 1 seed. The isotropic-crossing baseline is
polar 0.134 / mid 0.366 / equatorial 0.500.

Contacted-exit bands:

| body-contact law | τ_body | polar | mid | equat | p vs n |
|---|---|---|---|---|---|
| twist | 1 | 0.108 | 0.348 | 0.543 | identical |
| twist | 4 | 0.101 | 0.355 | 0.543 | identical |
| twist | 64 | 0.119 | 0.381 | 0.500 | identical |
| **edge** | **1** | **0.042** | **0.169** | **0.789** | identical |
| edge | 4 | 0.074 | 0.321 | 0.605 | identical |
| edge + twist | 1 | 0.044 | 0.179 | 0.777 | identical |

- **The disc emerges mechanically under edge grip.** With nothing imposed (no
  emission histogram, no redirect), photons scattered by the spinning shells leave
  79% equatorial and 4% polar at τ_body = 1, vs 50% / 13% for an isotropic crossing.
  That is pole.pdf's sprinkler, measured.
  - It weakens as the body gets opaque (τ = 4: 60%) and is gone by τ = 64.
  - Twist grip against a non-spinning grain gives only a weak lean (54%).
- **Proton and neutron are identical in every run.** The torque's sign follows the
  outer z-level; nothing else depends on the sign pattern. The species-split intake
  (bias 0.5) changes nothing either.
- Exiting helicity is ≈ 0 and pole-to-pole passage is ≤ 4% for both.
- **Why they are identical:** in this model body handedness enters only as the
  direction of the matter flow V(x). The outer z-level dominates V because the flow
  speed grows with radius. Reversing a flow field gives the same scattering statistics
  unless something in the photons prefers one sense. The photon's own helicity could
  do that, but a few contacts scramble it to ≈ 0 before it matters. The body and the
  photon never exchange *handedness* information in a way that survives.
- **Open (user):**
  - (i) Should photon-body contacts use edge grip (the sprinkler) while photon-photon
    contacts keep twist grip (the DirRel table)?
  - (ii) Where should body handedness enter mechanically? Candidates:
    - grains that carry the level's spin (twist grip against a spinning face, so
      chirality-matched photons pass and mismatched ones catch);
    - level-resolved contact rates;
    - the neutron's "returns with reversed spin" as a real mechanism.

### 14.8 Spinning grains + separate body law (2026-09-29, 43dba30)
- **Built:** `RoughParams::grain_edge`. Body grains spin with their shell (axis =
  the shell matter's total angular velocity, edge speed c). The twist grip acts on
  the *relative* spin, so a matched photon meshes, and a mismatched one is locked to
  the grain's sense, i.e. reversed (the neutron.pdf "returns with reversed spin"
  mechanism, per contact). `RoughParams::body_law` separates the photon-body law
  from the photon-photon one. GPU == CPU to 2e-6.
- **Result: proton == neutron still, in all 48 runs** (3 body laws × τ 1, 4 × 8
  states), with grain spin on.
- **Correction to §14.7:** the emergent disc comes from the **photon-photon edge
  grip**, not from the body-contact law (τ = 1, grain 0, contacted-exit bands):

| photon-photon | photon-body | polar | mid | equat |
|---|---|---|---|---|
| twist | edge | 0.101 | 0.342 | 0.557 |
| edge | twist | 0.041 | 0.149 | 0.810 |
| edge | edge | 0.043 | 0.171 | 0.786 |

  Under edge grip the gas carries stacked, spinning photons, and those interact with
  the moving shells very differently from twist-grip photons. The mechanism is not yet
  diagnosed; next diagnostics are the level distribution at entry and contacts by
  level. Spinning grains with an edge body law *remove* most of the disc (0.53).
- **Why p and n stay identical, structurally:**
  1. **Volume.** The inner levels' shells are small: level 11's region is ⅛ of the
     body volume, level 10's 1/64, level 9's 1/512. So the outer z-level sets the
     matter flow almost everywhere, and proton +a+x+y+z and neutron −a−x−y+z share
     the outer z-sense.
  2. **Symmetry.** In an isotropic ambient field, a body whose every level is
     reversed (proton state 1 vs state 0) is statistically a rotated copy. Photon
     helicity (a pseudoscalar) can only tell them apart if the body has a
     pseudoscalar property coupling to it. Relative level senses (ω_i · ω_j) are
     scalars, so p vs n is a *scalar* structural difference: inner levels co- vs
     counter-rotating with the outer one. That must show up through the inner
     shells' own flow, and in a volume-filling model hardly any photons sample it.
- **Implication for the next step (user, 2026-09-29: "stop and think about using
  what we find to build a new spin trace"):** the p/n difference lives in the inner
  levels, which a volume-averaged shell model under-weights. The Spin-mode trace
  (option c) resolves exactly where the inner levels' matter is at each instant, and
  photons that reach it (through the polar channel) are the ones that can see the
  difference. That makes the trace rebuild the natural next test, not just a
  cross-check.

## 15. New spin trace: the instantaneous ball (2026-09-29)

User decisions:
1. **Targets:** try both and compare. Instantaneous ball vs the §14 smeared shells.
2. **Levels:** start with levels 1–8 as one spherical base particle (the electron-tier
   ball, as in the earlier papers), traced through levels 9–12. Go deeper if results
   say a factor is missing. The user notes the tiers look fractal.
3. **Pole-to-pole:**
   - A perfectly axial photon flies through untouched, but that is rare.
   - Watch whether the entry angle to the pole is steep enough to brush the "donut"
     wall.
   - Watch how often a shallow (potential through-exiting) photon has a photon-photon
     hit inside the envelope that pushes it into the wall it was angled away from.
     That could be tied to the neutron's charge suppression.
4. **Time scale:** hits only on the instantaneous location of the real body, never a
   swept volume. Misses must stay in the simulation, because they hit each other.

**Model:**
- **The ball.** Levels 1–8 collapse to a ball of radius amplitude₈/amplitude₁₂ = ⅛.
  That is exactly §14's core shell: centre o₉ (where levels 10–12 carry the level-9
  frame origin).
  - Centre velocity: V_c = Σ_{j≥10} ω_j × (o₉ − o_j).
  - Spin: Ω_b = Σ_{j≥9} ω_j. Level 9's own precession spins the ball in place.
- **All other matter is gone.** Levels 9–12 are *motion*, not material. The envelope
  (level-12 ghost, radius 1) is used only for passage bookkeeping.
- **Continuous collision.** The CPU uploads 9 ball keyframes per step (centre, V_c,
  Ω_b at 8 sub-intervals). The move kernel solves photon-segment vs moving-sphere in
  each sub-interval, allowing one hit per step, then continues the rest of the step
  on the new velocity.
- **Contact.** The geometric normal is the photon→centre direction at the hit, with
  no sampling.
  - The sliding slip uses the true surface velocity V_c + Ω_b × (x − c).
  - The twist grip uses the grain-spin knob (0, or along Ω_b at edge c): the ball's
    photon-scale surface spin (levels 1–8), which is not resolved.
  - Same rough-sphere law, D7 PinnedC, stacking.
- **Every photon stays in the DSMC gas.** Photon-photon collisions inside the envelope
  are live and flagged per passage.
- **Statistics.** Hits scale with n_a. Raising n_a while scaling σ down keeps λ (the
  physics) fixed, and a λ-held-fixed comparison checks that.

**Passage tallies (body 0):** entry-angle-to-pole bin (10° steps, folded to 0–90°)
× outcome × species:
- outcomes: 0 clean; 1 photon-photon only; 2 ball hit first; 3 ball hit after a
  photon-photon hit;
- plus the opposite-pole exit count per angle bin and species, and exit helicity.

**Gates:**
- GPU vs CPU continuous hit + contact (replicated keyframes);
- a static ball conserves photon energy and delivered momentum balances;
- a ball with radius → 0 is never hit;
- the swept-centre map of the ball (CPU) shows the donut's polar hole, for comparison
  with Spin mode's disc-and-hole.

### 15.1 Bug found, results rerun, first proton/neutron split (2026-09-29, 4f10674)

**Bug: permeable bodies carved gas volume.** The per-cell precompute still treated
every body as a solid sphere (the Absorber rule). Cells inside a permeable envelope
got near-zero free volume, so their photon-photon candidate count was capped only by
the per-cell limit.
- The in-envelope photon-photon rate came out 2.8× the expected rate at λ = 3, and
  ~40× at λ = 10⁶.
- It was caught because the §15 passage table showed 92% of passages with a
  photon-photon collision, vs ~45% expected.
- Fixed: permeable bodies use the body-free volume. Regression gate
  `permeable_body_does_not_distort_gas_collision_rate`: 0 cells carved, rate 0.657 vs
  0.667 expected.

**§14 shell results rerun (supersede §14.7/§14.8 numbers).** τ_body = 1, grain 0,
PinnedC, stacking; contacted-exit polar/mid/equatorial:

| photon-photon | photon-body | bands | pole-to-pole |
|---|---|---|---|
| twist | twist | 0.084 / 0.333 / 0.583 | 17–19% |
| twist | edge | 0.065 / 0.315 / 0.620 | 7–8% |
| edge | twist | 0.083 / 0.334 / 0.583 | 17–19% |
| edge | edge | 0.064 / 0.315 / 0.621 | 6–10% |

- **The 79–81% "disc" was an artifact of the bug.** So was the attribution to
  photon-photon edge grip. The real shell sprinkler is a modest equatorial lean
  (58–62% vs 50% isotropic), slightly stronger with edge body contacts.
- Pole-to-pole passage is real (7–19% of contacted passages at τ = 1).
- **Proton == neutron still holds in every shell run.**

**§15 instantaneous ball (report_c, all 8 states, λ = 3, n_a = 4000, 1200 steps,
1 seed; free speed; twist grip; grain 1; stacking).** Gas density near the body
is 1.1× n_a, and ~46% of passages have a photon-photon hit, as expected.

| | contacted-exit bands (p/m/e) | helicity of ball-hit exits | torque / hit |
|---|---|---|---|
| proton, all 4 states | 0.105–0.107 / 0.339–0.345 / **0.548–0.556** | **+0.23 to +0.24** | ±0.65–0.67 |
| neutron, all 4 states | 0.109–0.116 / 0.352–0.360 / **0.526–0.532** | **+0.11 to +0.13** | ±0.70–0.71 |

- **This is the first proton/neutron split from mechanics alone,** consistent across
  all four states of each class. There is no imposed profile or redirect, and the
  only input difference is the level-sign pattern.
  - The neutron exits ~2.5 points *less* equatorially (~8σ).
  - The neutron carries about *half* the helicity (0.12 vs 0.23).
  - In words: the neutron sprinkles less disc and less magnetism.
- **The pole-angle table is identical for p and n:**
  - Through-passage is 8.7% for entries within 10° of the pole, 5% at 10–20°,
    1.7% at 20–30°, and ≈0 beyond. That is the built-in polar hole of the swept
    donut.
  - "Photon-photon then ball hit" (the user's knocked-into-the-wall case) is
    1.1–1.3% of passages at every angle.
  - Photon vs antiphoton show no split. The neutron's suppression is **not** yet
    visible in the pole channel; it shows up in *how the ball treats the photons it
    hits*.
- **Why the ball separates p and n when the shells could not:** the ball's real
  instantaneous velocity and spin history differs between sign patterns (mean centre
  speed 1.78c proton vs 1.67c neutron). Time-averaging into shells erased exactly
  that. This vindicates the user's instantaneous-hit decision.
- **PinnedC + ball:** the ball moves at ~1.8c, so a photon held at c cannot get
  away. 88% of hits are push-outs (snowplow), and some photons are overrun. The split
  persists (equatorial 0.57 vs 0.53, helicity 0.09–0.13 vs 0.04–0.07), but this
  combination is inconsistent. **Open decision:** free speed for body contacts
  (photon.html: "newly emitted photons may be going slightly over c"), or reconsider
  the superluminal compounding of Spin-mode level speeds.
- **Caveats:** one seed, one λ, one τ-free geometry, and levels 1–8 as a rigid ball.
  Magnitudes are small (the band split is 2.5 points); the helicity split is large.

### 15.2 BodyFree speed mode: the ball may throw photons off c (2026-09-29)

**User decision:** option 1 from §15.1. Body contacts are free, so a baryon surface
can throw a photon faster (or slower) than c. The sources:
- charge3.html: "The speed of light is the speed at which it is ejected by the
  baryon"; the photons emitted by the electron "are going somewhat faster than c.
  The speed of light may be an average."
- photon.html: newly emitted photons may go slightly over c.

The backup, if the speed question gets out of hand: revisit how Spin mode compounds
level speeds so the ball stays under c.

**Rule (`SpeedMode::BodyFree`, GPU speed_mode 2):**
- **Body contacts:** full rough-sphere mechanics, with no pin.
- **Photon-photon collisions:** still pinned, but each photon returns to its *own
  incoming speed*, not to c (`rough_sphere::pin_pair_to`). Each linear energy is
  unchanged, and the pair's spin energy is restored with one common factor as
  before. Energy stays exact, and with every photon at c this is PinnedC exactly.
- Only a body contact can change a photon's speed. An off-c speed is carried until
  the next body contact or escape.
- g_max is 12 (8 left a handful of overflow candidates per run, out of 47M
  collisions).

**Gates:**
- `rs_gpu_matches_cpu_reference` has a BodyFree + stacking case with speeds of
  0.5–1.5: error 1.0e-6.
- `ball_gpu_matches_cpu` has a BodyFree case: 0 hit/miss disagreements, error
  8.9e-6, 0 photons left inside the ball.

**Results** (`report_c_ball_single_body`, CFM_C_SPEED=bodyfree, twist/twist, grain
1, stacking, λ = 3, n_a = 4000, warmup 400, avg 1200; all 8 states on seed 0xC0C0,
plus seeds 1–3 on the first state of each class):

| | push-outs / run | escape speed | contacted-exit equatorial | helicity of ball-hit exits | torque / hit |
|---|---|---|---|---|---|
| proton (4 states + 3 seeds) | 8–30 | 1.001 ± 0.058 | **0.547–0.557** | **+0.136 to +0.147** | ±0.66–0.68 |
| neutron (4 states + 3 seeds) | 3–20 | 1.001 ± **0.063–0.064** | **0.528–0.535** | **+0.066 to +0.079** | ±0.71–0.73 |

- **The snowplow is gone:** push-outs drop from ~175k per run (PinnedC) to ≤30, and
  no photons are overrun.
- **c emerges as an average:** photons leave at 1.001c with a ~6% spread, never
  imposed. The neutron's spread is consistently wider (0.063–0.064 vs 0.058–0.059,
  in every state and seed). That is a small third p/n tell.
- **The p/n split holds under all three speed rules** (Free, PinnedC, BodyFree) and
  across seeds. The neutron is ~2 points less equatorial and carries ~half the
  ball-hit helicity (0.07 vs 0.14; the absolute helicity is lower than Free because
  photon-photon collisions are pinned).
- **The pole table is unchanged and identical for p and n:**
  - through-passage 9.4% (0–10°), 5.7% (10–20°), 1.7% (20–30°), 0.3% (30–60°);
  - pp→hit 0.8–1.1% at every angle;
  - no photon/antiphoton split.
- **BodyFree is the new default** for the ball report.

### 15.3 Trace-style histograms: new ball vs old Spin-mode trace (2026-09-29)

**New tallies:**
- `TALLY_HIT_LAT`: 180 one-degree bins, body 0, the latitude of each ball-contact
  point about the body centre and pole. This is what Spin mode's exit histogram
  binned: the position of each hit photon relative to the focus centre
  (field_sim.rs).
- `TALLY_EXIT_LAT`: 360 bins, [species][bin], the latitude of each contacted exit on
  the envelope.
- `report_c_ball_single_body` writes per-class CSVs (summed over states) with
  `CFM_C_HIST_DIR`. `CFM_C_SIGNS="a,x,y,z"` runs one custom state; 0 means a torpid
  level.

**Old Spin-mode presets** (godot/config/particle_presets.json):
- Levels 1–8 are all +1 for every particle, so going deeper can use +1 without a
  new assumption.
- The old neutron preset (−a−x−y) ran its z level at **0** (torpid), while the
  proton had z at c. The new neutron states have z at ±c.

**Results** (BodyFree, all 4 states per class, avg 2400; custom torpid-z run avg
4800). Folded |latitude|:

| | eq <30° | mid 30–60° | polar >60° | mean \|lat\| | peaks |
|---|---|---|---|---|---|
| OLD trace proton | 0.714 | 0.285 | 0.002 | 20.1° | 30° |
| **NEW ball proton, contact points** | **0.740** | **0.260** | **0.000** | **19.5°** | **28°** |
| OLD trace neutron (z torpid) | 0.406 | 0.576 | 0.019 | 34.1° | 10, 24, 55° |
| NEW ball, same torpid-z preset | 0.660 | 0.340 | 0.000 | 22.6° | 28° |
| NEW ball neutron (z at c), contact points | 0.628 | 0.372 | 0.000 | 24.0° | 2, 11, 30° |
| NEW ball proton, contacted exits | 0.551 | 0.341 | 0.107 | 30.1° | ~4–10° |
| NEW ball neutron, contacted exits | 0.529 | 0.360 | 0.112 | 31.2° | ~4–17° |
| OLD recycling proton / neutron | 0.584 / 0.308 | 0.345 / 0.482 | 0.071 / 0.210 | 28.1 / 41.5° | 5 / 44° |

- **The new ball reproduces the old proton trace almost exactly:** eq 0.74 vs 0.71,
  mean latitude 19.5° vs 20.1°, peak 28° vs 30°, no polar hits. The new engine's
  kinematics and the old trace agree where they should.
- **The old neutron trace is NOT reproduced**, even with the identical torpid-z
  preset (eq 0.66 vs 0.41). The cause is the hit geometry:
  - The old trace hit photons against a "snake": capsules along the particle's
    recent path, held static over one full turn of the second-outermost orbital
    (`snake_trace_length`). That weights hits by *where the path is*.
  - The instantaneous ball weights hits by *where the ball spends its time*. For the
    proton's path the two agree; for the torpid-z neutron path they do not.
  - The old neutron's mid-latitude butterfly is therefore largely a product of the
    swept-volume method the user chose to retire (§15).
- **New p/n split in the contact geometry itself:** the neutron's ball spends more
  time at mid-latitudes (contact eq 0.63 vs 0.74, mean 24.0° vs 19.5°).
  - Photon-photon scattering on the way out washes most of this out: the
    contacted-exit split is only 0.529 vs 0.551.
  - No polar exit excess for the neutron (0.112 vs 0.107), unlike the old recycling
    profile (0.21). Photon vs antiphoton exits are identical within noise.
- Overlay figure: scratch ball_vs_trace.svg (not committed).

### 15.4 Going deeper: resolving levels 8, 7, 6, 5 (2026-09-29)

**Model.** `ShellBody::depth` (0–4) resolves electron-tier levels 8, 7, 6, 5 as
motion, using Spin mode's own `SpinLevel` rules (level 8 = Z orbit, 7 = Y, 6 = X,
5 = axial precession about Z). The base ball is the outermost *unresolved* level's
sphere.

| depth | ball | radius | centre speed, p mean/max | n mean/max | fastest surface |
|---|---|---|---|---|---|
| 0 | levels 1–8 (spun by 9) | 1/8 | 1.78 / 2.30c | 1.67 / 2.52c | 3.8–4.3c |
| 1 | levels 1–7 (carried by 8) | 1/16 | 2.21 / 3.83c | 1.62 / 3.03c | 3.9–5.4c |
| 2 | levels 1–6 | 1/32 | 2.53 / 4.73c | 1.90 / 3.21c | 4.3–6.2c |
| 4 | levels 1–4 (spun by 5) | 1/64 | 2.72 / 4.86c | 2.19 / 4.01c | 5.5–6.1c |

- **Tests:** `deeper_ball_kinematics` checks depth 0 against the §15 ball exactly,
  checks centre and surface velocity against finite differences at every depth (4e-9),
  and checks that every ball nests inside its parent.
- **Reports:** `report_deeper_ball_map`; `CFM_C_DEPTH`, `CFM_C_DT` and `CFM_C_GMAX`
  in the ball report.
- **Keyframe interpolation error:** ≤2% of the ball radius at depths ≤2, 6% at depth
  4, so depth 4 runs at dt = 0.025.
- **Inner levels use +1 sign** (every Spin-mode preset).
- **The polar axis is never blocked at any depth.**
- **Geometric p/n difference appears with depth:** the proton ball reaches the 40–50°
  latitude band 12–13% of the time, the neutron 5–6% (identical 2% at depth 0).
- **Speed compounds with depth,** faster for the proton. This is the "speed out of
  hand" risk; photons still escape at 1.000c on average at every depth.

**Results** (BodyFree, first state of each class, λ = 3; depth 1/2/4 averaged over
4800/9600/38400 steps; hits are per 2400 steps at dt 0.05, or the same physical time
at depth 4):

| depth | hits/2400 steps p / n | n/p | contacted-exit equatorial p / n | ball-hit helicity p / n | contact eq p / n | contact mean \|lat\| p / n |
|---|---|---|---|---|---|---|
| 0 | 44.1k / 44.2k | 1.00 | 0.547–0.557 / 0.528–0.535 | +0.14 / +0.07 | 0.740 / 0.628 | 19.5° / 24.0° |
| 1 | 13.7k / 11.1k | **0.81** | 0.537 / 0.523 | +0.062 / +0.070 | 0.718 / 0.647 | 19.7° / 22.5° |
| 2 | 3.9k / 3.1k | **0.81** | 0.542 / 0.534 | +0.148 / +0.155 | 0.723 / 0.660 | 18.8° / 21.4° |
| 4 | 1.01k / 0.86k | **0.85** | 0.548 / 0.541 | +0.099 / +0.115 | 0.721 / 0.688 | 19.1° / 20.9° |

- **The depth-0 magnetic (helicity) split does not survive resolution.** From depth 1
  on, p and n helicity agree within noise (±0.006–0.011). The resolved inner levels
  (+1 for both classes, and the fastest: rates 16–64 vs 2–8) dominate the ball's
  surface twist, which carries the helicity. The depth-0 split came from folding the
  level-9 precession into the ball's own spin.
- **The contacted-exit equatorial split shrinks with depth:** 2 points at depth 0,
  1.4 (3σ) at depth 1, then 0.7–0.8 (under 2σ). Photon-photon scattering on the way
  out washes out most of the geometry.
- **The proton's contact profile is depth-invariant:** eq 0.72–0.74, mean latitude
  19–20° at every depth. That matches the old Spin-mode trace, and the user's
  observation that the tiers are fractals of each other.
- **The neutron's contact geometry stays more mid-latitude** (0.65–0.69 vs 0.72),
  converging slowly toward the proton with depth.
- **New robust p/n difference: once the inner levels move, the neutron is hit
  15–19% less often.** At depth 0 the rates are equal. Its slower ball sweeps less
  volume, so a neutron engages (recycles) less of the ambient field than a proton,
  in every depth tested. A candidate mechanical piece of neutron neutrality: less
  intake rather than a different exit profile.
- **Pole table unchanged at every depth:** through-passage 10% / 5.8% / 1.8%. With
  hits 30× rarer at depth 4 and the same through fractions, the "through" channel is
  a property of the photon gas and the envelope geometry, not of the ball.
- **Cost:** hits scale slower than the ball's cross-section (a faster ball sweeps
  more). Depth 4 needs ~40× the run time of depth 0 for similar statistics.

### 15.5 Density sweep: does photon-photon crowding drive the neutron difference? (2026-09-30)

**Setup:** depth 1 (the depth where the intake difference shows), BodyFree, first
state of each class, λ = 3 → 1 → 0.3 (photon-photon mean free path, body radius 1).
- λ = 0.3 needs a finer grid (h = 0.1 instead of 0.25) and n_a = 15000 (13.4M live
  photons; the 2 GiB binding limit caps capacity near 19M). Warmup 2400, avg 4800.
- Collision rates came out at 0.434 / 1.30 / 4.36 per photon per unit time, as
  expected for each λ.

| λ | through-pass (0–10° / 10–20°) | clean passages | hits per passage (hit + pp→hit), p / n | share of hits after a pp collision, p | n/p hit rate | contacted-exit eq, p / n | ball-hit helicity, p / n |
|---|---|---|---|---|---|---|---|
| 3 | 10% / 5.8% | 57% | 0.0089 / 0.0073 | 30% | 0.81 | 0.537 / 0.523 | +0.06 / +0.07 |
| 1 | 3.4% / 2.3% | 21% | 0.0090 / 0.0074 | 64% | 0.82 | 0.540 / 0.529 | +0.05 / +0.05 |
| 0.3 | 0.4% / 0.4% | 2.6% | 0.0086 / 0.0070 | 95% | 0.805 | 0.548 / 0.537 | +0.04 / +0.04 |

- **Density changes which hits come after a photon-photon collision, not how often
  the ball is hit.**
  - At λ = 0.3, 95% of ball hits follow a photon-photon collision inside the
    envelope: the user's "knocked into the wall" channel dominates.
  - But hits per passage are flat (0.0086–0.0090), and hits per photon density are
    flat too.
- **The neutron's lower intake is density-independent** (n/p 0.81, 0.82, 0.805). It
  is geometric (a slower ball sweeps less space), not caused by internal
  photon-photon hits. The internal-hit hypothesis is not supported, as the user
  suspected.
- **The contacted-exit equatorial split is about 1 point at every λ.** At λ = 0.3
  (100k hits) that is ~5σ: small but real, and also density-independent.
- **The pole-to-pole channel closes with density:** 10% → 3.4% → 0.4% for near-pole
  entries, identical for p and n. It is a gas-transparency effect.
- **Density scrambles helicity.** The photon-photon passage helicity falls from
  +0.066 to +0.002, and the ball-hit helicity falls with it; p ≈ n throughout.
- **Speed tail:** g_max-overflow pairs grow with density (proton 685 of 7.0e9
  collisions at λ = 0.3; neutron 3). Negligible bias, but the runaway lives near the
  faster proton ball. Escape speed is unchanged (1.000 ± 0.050 p, ± 0.037 n).

### 15.6 Neutron-model alternatives, all states; charge-pause map (2026-09-30)

**Question (user):** which neutron are we using? The new engine uses Mathis's QCD
sign states with all four baryon levels at c. Old Spin mode ran the neutron with
z = 0 (torpid); the old recycling engine used z = ±0.05c. None was fully
satisfying, so each variant was run through the same test.

**Report additions:**
- `CFM_C_TOP_P` / `CFM_C_TOP_N`: level-12 edge speed per class.
- `CFM_C_INNER`: sign of levels 8–5.
- `CFM_C_INNER_P` / `CFM_C_INNER_N`: per-class "s8,s7,s6,s5".
- `CFM_C_ONLY`: run one class. `CFM_C_SIGNS` accepts ';'-separated states.
- Ball-centre mean speed per state.
- `TALLY_FLUX` (952, len 504) and `FLUX_RADII`: shell crossings by radius ×
  latitude × species × direction; the CSV is `ball_<class>_flux.csv`.

**All at λ = 3, BodyFree, all 4 states per class.** Depth 1 at avg 2400, depth 2 at
avg 4800.

| model | depth | ball speed by state | hits n/p (state-summed) | contacted-exit eq, p / n (state mean) | ball-hit helicity by state |
|---|---|---|---|---|---|
| N1: all levels at c, Mathis signs | 1 | p 2.21/1.82/1.93/2.22, n 1.62/2.21/1.89/2.22 | **0.99** | 0.550 / 0.530 | p +.06/+.15/+.14/+.06, n +.07/+.04/+.06/+.03 |
| N2: neutron z torpid (0) | 1 | n 1.43/1.95/1.65/1.95 | 0.90 | – / 0.515 | n +.14/+.07/+.11/+.06 |
| N3: neutron z 0.05c | 1 | n 1.42/1.95/1.65/1.95 | 0.88 | – / 0.517 | n ≈ N2 |
| N4: BOTH tops at 0.05c | 1 | identical sets | 1.00 | 0.515 / 0.517 | identical |
| N5: inner levels −1 | 1 | p/n speeds permuted vs N1 | 0.99 | 0.547 / 0.530 | permuted |
| N6: user's "z flipped in every tier" neutron (level 8 = −1) | 1 | = N5 neutron | – | – / 0.530 | n +.02/+.07/+.03/+.05 |
| N6 | 2 | p 2.53/1.99/2.17/2.41, n 2.47/1.97/2.47/2.11 | 0.99 | 0.546 / 0.541 | p +.14/+.05/+.14/+.16, **n −.10/−.08/−.09/+.01** |
| control 1: PROTON with the tier-2 z flip | 2 | 2.07/2.45/2.52/2.13 | – | 0.545 / – | **+.03/−.08/−.08/−.07** |
| control 2: neutron, no flip | 2 | 1.90/2.44/2.14/2.46 | – | – / 0.535 | +.16/+.14/+.01/+.14 |

**The §15.4 "neutron hit 15–19% less" is retracted: it was a single-state
artifact.**
- Hits track the ball's mean speed, which depends on each state's outer signs
  *relative to the inner +1 levels*, not on the class.
- Summed over its 4 states, each class spans the same speeds (1.6–2.2c), and
  n/p = 0.99.
- Flipping the inner sign (N5) permutes which states are fast.

**What survives every state and both inner signs: the neutron's contacted exits are
~2 points less equatorial** (0.530 vs 0.547–0.550 at depth 1; 0.5–1 point at
depth 2).

**Torpid or sub-c top (N2/N3):**
- It lowers the disc for whichever particle has it (0.515).
- If both classes have it (N4), p and n are identical, because Mathis's p and n
  states differ only in z relative to a, x, y.
- The old "neutron = no outer spin" difference was therefore a
  spun-top-vs-torpid-top comparison, not a sign-pattern effect.
- magmom.pdf ("the proton must be spinning faster than the neutron") would support
  a slower neutron top, but nothing yet fixes how much slower.

**User's N6 ("z flipped in every tier"):**
- At depth 2 (levels 8 and 7 moving), the neutron's ball-hit photons come out with
  *reversed* helicity against the net-positive ambient field (−0.08 to −0.10 in 3 of
  4 states, errors ±0.012), while the proton's are +0.05 to +0.16.
- Controls:
  - A proton given the same tier-2 z flip also reverses (3 of 4 states).
  - A neutron without it does not.
- So the reversal comes from the inner tier's z-vs-rest arrangement. Under the
  user's definition, where neutron-ness propagates to every tier, it becomes a
  clean p/n magnetic split with the sign neutron.pdf needs ("returns with reversed
  spin").
- Level 4 (tier-1 z) sits inside the rigid ball at every depth we can run, so its
  flip is untested.
- **Caveat:** the output helicity is set by the twist grip against the ball's
  surface spin, whose direction comes from the §14.8 grain-spin rule (grains carry
  the ball's net Ω at edge c). That rule is an open assumption.

**Charge-pause map (λ = 3, depth 1): null.**
- Net radial flux per latitude band is ≤1e-3 of crossings and within 1–2σ at every
  radius 0.5–5.
- There is no polar-in / equatorial-out cell. The only consistent signal is a
  uniform slight inflow near the outer boundary, from how the ambient gas is
  injected there.
- The ball touches <1% of the photons that pass through the envelope, too little to
  drive a circulation.
- A charge-pause needs the recycling (emission) side, or a much stronger
  interaction.

### 15.7 Adopted neutron: z flipped in every tier; confirmation (2026-09-30)

**User decision:** adopt the "z flipped in every tier" neutron. In the ball report
the neutron's levels 8, 7, 6, 5 default to (−1, +1, +1, +1) (`NEUTRON_INNER_SIGNS`)
on top of the Mathis outer states. The proton keeps all +1. The outer-only neutron
is `CFM_C_INNER_N="1,1,1,1"`. Level 4 is not represented (inside the ball).

**Ball-hit helicity by state** (λ = 3, BodyFree, all 4 states per class; ±0.012 at
depth 2, ±0.018 at depth 3):

| run | proton | neutron | class means p / n |
|---|---|---|---|
| depth 2, seed 0xC0C0 (§15.6) | +.14 / +.05 / +.14 / +.16 | −.10 / −.08 / −.09 / +.01 | +0.12 / −0.07 |
| depth 2, seed 1 | +.13 / +.05 / +.14 / +.17 | −.09 / −.08 / −.11 / .00 | +0.12 / −0.07 |
| depth 2, seed 2 | +.16 / +.07 / +.13 / +.15 | −.09 / −.10 / −.11 / +.03 | +0.13 / −0.07 |
| depth 3 (dt 0.025, 16000 steps) | +.17 / +.10 / +.15 / +.17 | +.02 / +.08 / +.05 / +.12 | +0.15 / +0.07 |

- **Depth 2 is reproducible across 3 seeds.** The adopted neutron reverses the
  helicity of the photons it touches in 3 of 4 states (−0.08 to −0.11); the
  +a−x−y+z state stays at ~0 in every seed.
- **At depth 3 the reversal weakens to a halving.** Resolving level 6 (X spin,
  rate 64, +1 for both classes) makes the fastest resolved level common to p and n.
  The neutron's ball-hit helicity is still about half the proton's in every state,
  so the p/n split survives, but not as a sign reversal.
- **The depth trend tracks the §14.8 grain-spin rule.** The surface twist a photon
  meets is the direction of the ball's summed Ω, which the fastest resolved levels
  dominate. The magnitude, and even the sign, of the magnetic split therefore
  depends on how the ball's surface twist is modelled. **Deciding that rule
  mechanically is now the key open question** (see the options in the session
  notes).
- **The contacted-exit equatorial split fades with depth:** ~1 point at depth 2
  (0.542 vs 0.533–0.535, both seeds), ~0 at depth 3 (0.545 vs 0.543).
- Hits, ball speeds and torque per hit: no class-level difference beyond the state
  spread.

### 15.8 The surface-spin rule: rigid and rigid + inner levels (2026-09-30)

**Question:** §15.7 showed the magnetic split depends on the §14.8 grain rule (twist
reference = direction of the ball's summed Ω at edge c). The user asked for option
2 (rigid) to get a fair go, and then option 3 (inner levels) to see if it helps.

**Build:**
- `RoughParams::grain_mode: GrainMode {EdgeC, Rigid, RigidInner}` (GpuParams
  `grain_mode`, formerly `_pad14a`).
- **Rigid:** the twist reference is the ball's actual angular velocity in photon
  spin units, Ω × photon radius, with photon radius = level-1 / level-12 amplitude
  = 1/512 of the body radius. No assumed edge speed.
- **RigidInner:** Rigid plus the unresolved levels 1..(8 − depth), using the
  Spin-mode `SpinLevel` rates and the class signs:
  - rates in photon units: 1, 1, ½, ¼ (levels 1–4), ⅛, ⅛, 1/16, 1/32 (levels 5–8);
  - each level's axis is carried by the unresolved levels above it and the ball
    frame, each at a random phase per contact (hash RNG; their periods are far
    below the keyframe spacing);
  - `ShellBody::{tier1_signs, unresolved_levels, photon_radius}`,
    `BallState::frame`, `baryon_shells::inner_twist`.
- **Keyframes** are now 4 vec4 per keyframe (the frame quaternion is added) plus
  the inner-level table; `BALLKEY_STRIDE` = 45.
- **The adopted neutron flips level 4 as well** (`CFM_C_TIER1_N`, default
  (−1, +1, +1, +1)).
- **Gates:**
  - `ball_gpu_matches_cpu` adds Rigid (error 8.9e-6) and RigidInner (z-flip neutron,
    depth 2, random phases mirrored exactly: error 4.8e-5).
  - The control compares spins too; under Rigid, dropping the ball's spin changes
    photons by only 2.5e-2.
  - `unresolved_levels_rates_and_frame` checks the rates and the frame.
- **Report:** `CFM_C_GRAIN_MODE` = edge | rigid | inner.

**Results** (λ = 3, BodyFree, all 4 states per class; ball-hit helicity, state
means; errors ±0.006 at depths 0–1, ±0.006 at depth 2, ±0.009 at depth 3):

| grain rule | depth | proton | neutron | contacted-exit eq, p / n |
|---|---|---|---|---|
| EdgeC (§15.7) | 2 | +0.12 | −0.07 | 0.542 / 0.534 |
| EdgeC (§15.7) | 3 | +0.15 | +0.07 | 0.545 / 0.543 |
| Rigid | 1 | +0.018 | +0.021 | 0.548 / 0.534 |
| Rigid | 2 | +0.026 | +0.016 | 0.548 / 0.533 |
| RigidInner | 0 | +0.017 | +0.017 | 0.550 / 0.528 |
| RigidInner | 1 | +0.017 | +0.014 | 0.549 / 0.530 |
| RigidInner | 2 | +0.022 | +0.015 | 0.543 / 0.539 |
| RigidInner | 3 | +0.035 | +0.016 | 0.546 / 0.550 |

- **Under both mechanical rules the ball erases handedness for both classes.**
  Ball-hit helicity falls to +0.01–0.04 (from +0.137 on arrival), with no p/n
  difference beyond ~1.5σ at any depth.
  - **Rigid:** the ball's rigid spin is only ~0.03–0.14 in photon units, so the
    twist grip despins photons.
  - **RigidInner:** the inner levels give a twist of order c (levels 1–2 spin at
    exactly edge c in photon units), but their fast, mutually rotating axes make it
    a random direction at each contact.
  - The coherent part (level 4 and level 5, both along the frame's Z) is about 0.1–
    0.4 and flips for the neutron. Projected onto contact normals all over the
    ball, it averages out.
- **The EdgeC split (§15.6–15.7) came from the assumed coherent twist at edge c.**
  No mechanical version of the rule reproduces it.
- **The equatorial exit split survives the rule change at shallow depth:** ~1.4–2.2
  points at depths 0–2, gone at depth 3, as under EdgeC. It is geometric, not a
  twist effect.
- **Reading:** a photon that only bounces off the base ball cannot carry the
  neutron's magnetic signature out. In neutron.pdf the photons "go in, rattle around
  in the neutron, taking the four given spin paths" before leaving. The magnetic
  channel needs photons to travel *through* the spin structure (channelling /
  recycling), which the rigid-ball abstraction does not allow. That points back to
  the emitter / channelling design question.
- The report default stays EdgeC for continuity; Rigid / RigidInner are the
  mechanical references.

### 15.9 The maze as the ball's emission law (2026-09-30)

**Question:** docs/MATHIS_SPIN_PATHS.md showed on the CPU that Mathis's quark.html
maze reproduces his 16-state table as *core emission*. With the fixed Spin-mode
frame the emission is a one-sided lobe. The user chose (1) a **swept** frame, turned
about the pole by the core's level-9 axial phase (the sprinkler), and (2) the
neutron.pdf release: trapped photons funnel to the axis and leave by the
*opposite* pole with the maze's spin.

**Build:**
- `RoughParams::maze: MazeMode {Off, Fixed, Swept}` (GpuParams `maze_mode`,
  formerly `_pad14b`). This mode needs the Ball body model.
- **Every ball contact is a capture.** The photon is re-emitted by
  `maze::emit` (shader `maze_contact`):
  - **Hemisphere:** the photon's helicity in the maze frame picks the start
    patch: negative → north, around +x; positive → south, around −x.
  - **Start d0:** drawn by rejection from the patch, i.e. the starts that
    reproduce all 16 table rows with helical spin, with margin 1e-4.
    Acceptance is about 12% per hemisphere, so 64 tries almost always succeed;
    the fallback is the patch centre (about 4e-4 of captures).
  - **Spin transport:** the photon's spin is carried rigidly from its arrival
    direction to d0.
  - **The maze:** the canonical maze runs with the body's (a, x, y, z) signs:
    90° turns about Z, X, Y, Z, with Set contact spin.
  - **Outcome:** escape = d·d0 ≥ 0, and escapes leave along the maze exit.
    Trapped photons leave along (0, 0, h), the pole opposite the start
    hemisphere.
  - **Speed:** kept relative to the ball surface, so the emission rides the
    moving core. Under PinnedC the speed is c instead.
  - **Stacked levels:** inner stacked spins turn rigidly with the photon.
    |spin| is kept at every level.
- **Maze frame:** `orient`, times a rotation about the pole by θ₉ + ω₉·t
  (Swept) or by 0 (Fixed).
  - `BALLKEY_STRIDE` = 47. The meta vec4's z and w hold θ₉ at the step start
    (mod 2π) and ω₉; [45] = `orient`; [46] = signs.
- **Tallies:** `TALLY_MAZE` (outcome × upright, plus fallbacks),
  `TALLY_MAZE_LAT` (emission latitude), `TALLY_MAZE_AZ` (escape azimuth in the
  body frame). `TALLY_LEN` = 1677.
- **Gates:**
  - `maze::emission_law_matches_canonical_maze`: the exact 90° helpers equal
    `run()`; the patch test agrees with `score()`; the fallbacks lie in the
    patch; every state gives the table outcome; escapes are within 30° of the
    equator; trapped photons leave by the opposite pole.
  - `ball_gpu_matches_cpu` adds three maze cases (swept p, swept n, fixed n
    under PinnedC). Max error is 1.2–1.9e-5 with 0 hit/miss disagreements. The
    control (body z sign flipped) is off by 2–7.
- **Report:** `CFM_C_MAZE` = off | fixed | swept. It writes
  `ball_<class>_maze.csv`.

**Results** (λ = 3, BodyFree, default seed, all 4 states per class, adopted
neutron; contacted-exit bands at the R_D = 6 envelope, after the gas):

| run | hit % | proton p / m / e | neutron p / m / e | hit helicity p / n |
|---|---|---|---|---|
| no maze, Rigid (§15.8), depth 2 | — | — / — / 0.548 | — / — / 0.533 | +0.026 / +0.016 |
| maze **swept**, depth 2 | 0.20–0.26 | 0.053 / 0.276 / **0.671** | **0.214** / 0.488 / 0.298 | +0.092 / −0.001 |
| maze fixed, depth 2 | 0.20–0.26 | 0.063 / 0.296 / 0.642 | 0.211 / 0.491 / 0.299 | +0.078 / +0.008 |
| maze **swept**, depth 0 | 2.9–3.0 | 0.045 / 0.239 / **0.717** | **0.251** / 0.571 / 0.179 | +0.095 / +0.002 |

- **At the ball:**
  - Every proton state emits 100% within 30° of the equator. Every neutron
    state emits 100% from the poles.
  - **Swept vs fixed azimuth:** the swept frame fills the ring evenly (escape
    octants 0.09–0.16). The fixed frame gives **two opposite lobes**, not one:
    both helicities arrive, and each picks its own hemisphere and side. Which
    quadrants they fill depends on the state.
  - By the r = 6 envelope the latitude bands barely differ between swept and
    fixed, because the ball's orbital motion and the gas already smear
    azimuth.
- **The proton/neutron electrical split is now large:**
  - Contacted proton exits are 0.67–0.72 equatorial; neutron exits are
    0.18–0.30 equatorial and 0.21–0.25 polar.
  - Before, both classes were ~0.54 equatorial.
  - What the construction fixes: proton = disc, neutron = poles, from
    Mathis's table plus neutron.pdf's release rule, with nothing fitted.
  - What the run adds: the split survives transport through the gas. The
    neutron's polar jet widens a lot (100% → 21–25% polar at r = 6); the
    proton disc widens less.
- **The helicity split is back, mechanically.**
  - Photons that hit carry +0.09 out for protons and about 0 for neutrons.
  - Neutron states pair by the z sign: z = +1 gives +0.02 to +0.04, and z = −1
    gives −0.02 to −0.04. The trapped release goes along the pole chosen by
    the arrival helicity, and the last Set contact imprints the z sign along
    that pole.
- **Charge-pause flux map** (net (out − in)/(out + in) × 10³ by radius and
  band; the noise at r = 1 is ~1.6 at the poles and ~0.6 at the equator):
  - **Depth 2:** null, as before. The ball, radius 1/32, captures 0.25% of
    passages. What remains is the uniform +1 to +1.7 drift toward the
    boundary seen in every run.
  - **Depth 0** (ball radius 1/8, 3% of passages):

    | class | band | r = 1 | r = 2 | r = 4 |
    |---|---|---|---|---|
    | proton | poles | net intake −7.5 / −8.4 (S / N) | −4.6 / −1.7 | ≈ baseline |
    | proton | mid | −5.2 / −2.8 | | |
    | proton | equator | net outflow +5.4 | +2.5 | ≈ baseline |
    | neutron | poles | net outflow +11.8 / +25.0 (S / N) | | north pole still +4.3 |
    | neutron | mid | +2.8 / +18.3 | | |
    | neutron | equator | net intake −12.5 | | |

    - **The proton shows Mathis's recycling pattern:** intake at the poles,
      emission at the equator. The intake side was not imposed. It is where
      the capture profile exceeds the emission profile.
    - **The neutron is reversed:** intake at the equator, emission at the
      poles.
    - The neutron's north/south asymmetry comes from the ambient gas's net
      arrival helicity (+0.137), which picks the release pole. It is a
      property of the gas, not of the neutron.
  - **No clean in/out reversal (charge-pause) at a finite radius:** the
    signals decay into the boundary drift by r ≈ 4. The emitter handles only
    3% of the gas. A real pause probably needs the full recycling rate
    (Mathis: ~19× the mass per second) or a second body feeding the first.
- **Open questions:**
  - Neutron intake: in Mathis the neutron's through-charge enters at a pole.
    Here capture happens wherever the ball is hit (mostly equatorial). That
    makes the neutron's equator a net sink, which may be wrong.
  - Emission strength (the intake cross-section) is set by the ball radius,
    i.e. by the depth.
  - Upright fractions for protons are 0.51–0.77 by state, because arriving
    spins aren't helical. That is secondary.

### 15.10 Two maze bodies side by side (2026-09-30)

**Question:** do two §15.9 maze bodies feed each other, and does a charge-pause
form between them?

**Build:**
- **Emitter tag:** photon flag bits 17–19 hold 1 + the body whose maze last
  emitted the photon (0 = never captured).
- **`TALLY_FEED`** [capturing body × tag] counts every capture by who emitted
  the photon last.
- **`TALLY_LINE`** counts crossings of planes perpendicular to the body-0 →
  body-1 axis, every 0.25 from u = −2, inside radius 0.5 of the axis. It
  splits them by tag class (never captured / body 0 / body 1) and direction.
  `TALLY_LEN` = 1943.
- **Gate:** `two_body_feed_bookkeeping`. The feed sums equal the contacts per
  body, both bodies capture, and traffic crosses the mid-plane both ways.
- **Report:** `report_two_bodies` runs the isolated baselines, then pairs
  pp / pn / nn, side by side (separation along x, in each other's equatorial
  plane) or axial (along the pole). Env: `CFM_2_PAIRS`, `CFM_2_GEOMS`,
  `CFM_2_SEP`, plus the usual `CFM_C_*`.
- **Setup:** swept maze, BodyFree, depth 0, λ = 3, first sign state per class,
  adopted neutron.

**Results** (share of each body's captures that its neighbour emitted last;
separation 2.5 / 2.05, body envelopes of radius 1):

| pair | side by side | axial |
|---|---|---|
| p – p | 0.002 & 0.002 / 0.002 & 0.002 | 0.001 & 0.000 / 0.000 & 0.000 |
| p – n (p, n) | 0.000, 0.002 / 0.000, 0.003 | **0.004**, 0.000 / **0.005**, 0.000 |
| n – n (lower/left, upper/right) | 0.000 & 0.000 / 0.001 & 0.000 | **0.004** & 0.001 / **0.006** & 0.002 |

- **The feeding sorts by geometry the way Mathis's nucleus needs:**
  - A proton's equatorial ring feeds a side neighbour (p–p side; the neutron
    in p–n side) and never an axial one.
  - A neutron's polar jets feed an axial neighbour and never a side one. In
    p–n axial, the neutron's jet pours into the proton's pole, which is the
    proton's net-intake region (§15.9).
  - In p–n side, net flow along the line runs from the proton into the
    neutron's equator (the neutron's net intake region), +16 to +29 × 10⁻³.
  - Nothing here is fitted; it follows from the §15.9 emission law.
- **But the coupling is about 100× too weak to matter:**
  - The neighbour supplies ≤ 0.6% of captures.
  - Capture rates stay at 0.98–1.015 × isolated.
  - Gas density along the line stays at 0.95–0.99 n_a, with no pile-up in the
    gap.
  - The tag-decomposed line flux shows each body's own emission pattern
    unchanged by its neighbour: the streams pass through each other at λ = 3.
  - There is no charge-pause: the net flow between the bodies is just the sum
    of the two emission patterns.
- **Why:** each ball captures ~3% of the photons crossing its envelope. Only
  ~4–5% of the gas near a body has ever been re-emitted (the "share tagged"
  column). In Mathis's picture the near field of a nucleon is mostly its own
  recycled charge, so the intake has to capture most of what passes through,
  not ~3%.
- **Unexplained:** in n–n axial the lower body is fed ~3× more than the upper
  one (about 90 vs 30 captures), even though the release-pole bias
  (+0.137 arrival helicity → north) should favour the upper body. The counts
  are small.
- **Next lever:** the intake cross-section / emitter strength (open since
  §15.9). Until a body recycles most of the charge crossing it, neighbours
  cannot feed each other enough for a charge-pause or binding.

### 15.11 Intake scaffolding and the first emergent recycling loop (2026-09-30)

**User direction:** the end goal is MATHIS_SPIN_PATHS §D option 2, intake that
*emerges* from the field. Mathis: the equatorial emission drives incoming charge
to the poles, with feedback up to a stable maximum. Options 1 (a capture rule) and
3 (a holding tank) are scaffolding. They are allowed as steps because simulating
10¹¹ real photons is not. First test: can a body emit hard enough to start the
loop?

**Build:**
- `RoughParams::intake: Intake {Ball, Envelope{p}, Polar{aperture}}`.
  - A photon entering a maze body's envelope (radius 1) is taken in with
    probability p anywhere (Envelope), or when it enters within `aperture` of
    the axis (Polar).
  - It is funnelled to the core and leaves through the maze. Ball contacts
    still capture as in §15.9.
- **Holding tank:** `hold_time`, `hold_time_neutron`, `hold_cap`.
  - Captured photons are flagged held (flag bit 20; c[1] = release step) and
    ride the core, taking no collisions (collide skips pairs with a held
    photon).
  - At the release step they leave through the maze.
  - The engine tracks the count per body (`DsmcEngine::held`) and closes the
    intake while it is at the cap. The cap is soft: it lags one step.
- `MazeMode::Isotropic`: control that re-emits in a random direction, spin
  carried rigidly.
- Keyframe [47] per body = (intake open, hold steps, intake parameter,
  intake mode); `BALLKEY_STRIDE` = 48.
- Tallies:
  - `TALLY_HOLD` [b × in/out];
  - `TALLY_CAPT_LAT` (body-0 intake latitude, by never-captured / recaptured);
  - `TALLY_ENV_CAPT`;
  - `TALLY_LEN` = 2315.
- **Gate `intake_tank_bookkeeping`:**
  - p = 1 leaves no clean passages;
  - tank in − out = host count = photons flagged held;
  - a cap turns photons away;
  - the isotropic control emits 0.502 within 30° of the equator (0.5 expected);
  - the maze emits 1.000 there.
  - Default behaviour is unchanged: all earlier gates are bit-identical.
- **Reports:** `report_intake` (`CFM_I_RUNS` = "intake:emission:hold;...",
  `CFM_I_CLASS`) and `report_two_bodies` + `CFM_2_INTAKE`.

**Results** (depth 0, BodyFree, n_a = 4000, first sign state):
- **Emitter strength:** with envelope p = 1 a proton recycles 609 photons per
  step at λ = 3. Ball-only recycles 18.3 and polar aperture 1/3 recycles 54.
- **The option-2 test.** With intake equally likely everywhere on the
  envelope, where does the *never-captured* intake arrive? No field gives
  |lat| 0–30 / 30–60 / 60–90 = 0.500 / 0.366 / 0.134. Ratios to that
  expectation:

  | λ | proton maze polar / equatorial | isotropic control | maze ÷ control, polar | recaptured share |
  |---|---|---|---|---|
  | 3 | 1.100 / 0.949 | 1.043 / 0.980 | **1.05** | 3% |
  | 1 | 1.213 / 0.891 | 1.057 / 0.973 | **1.15** | 9% |
  | 0.3 | 1.456 / 0.783 | 1.092 / 0.955 | **1.33** | 27% |
  | 0.3, neutron | 0.705 / 1.161 | 1.072 / 0.972 | **0.66** | 22% |

  - **The proton's disc pushes its own intake toward the poles**, and the push
    grows with gas density: it is collisional, as in Mathis's "drives off
    ambient photons that are incoming, by direct bombardment" (season.pdf).
    Nothing in the rule prefers the poles.
  - **The neutron's polar jets do the reverse:** its intake moves to the
    equator. The mechanism is general (intake moves away from wherever a body
    emits). With our neutron release rule, that contradicts Mathis's
    pole-to-pole neutron, which takes in at a pole. This is a tension to
    resolve, not a confirmation.
  - The ~4–9% polar excess in the isotropic control is an unexplained
    baseline, common to both runs, so compare maze against control rather
    than against 0.134.
  - Emission shields intake: captures fall with density (609 → 458 per step
    from λ = 3 to 0.3). The isotropic control falls the same way, so this is
    not maze-specific.
- **Two bodies at λ = 1, envelope p = 1, separation 2.5.** Share of captures
  fed by the neighbour:

  | pair | side by side | axial |
  |---|---|---|
  | p – p | 0.096 / 0.098 | 0.042 / 0.050 |
  | p – n (p, n) | 0.051, 0.097 | **0.116**, 0.050 |
  | n – n | 0.050 / 0.050 | **0.111 / 0.169** |

  - This is up 20–30× from §15.10 (≤ 0.006).
  - The geometry sorting survives: the disc feeds a side neighbour, and the
    jets feed an axial one.
  - In p–n axial the neutron's jet crosses the gap into the proton with net
    flow −0.4 to −0.6 (neutron-tagged), a one-way axial handoff. In
    symmetric pairs the two emission streams meet head-on in the gap (e.g.
    p–p side at the midpoint: +0.35 / −0.33, net ≈ 0).
  - Capture totals stay 0.93–1.05 × isolated: fed charge replaces fresh
    charge rather than adding to it.
- **Status of the loop:**
  - Emission now shapes intake, which is the first half of Mathis's loop,
    and emerges.
  - The intake *amount* is still set by the scaffold (p = 1 at the envelope).
- **Next steps toward option 2:**
  - Lower p (or tie capture to the local incoming flux) and see whether the
    focusing still steers what gets in.
  - Test whether a self-limiting maximum appears (Mathis's "stabilizes").
  - Resolve the neutron intake tension.
  - The tank is built but untested for physics.

#### 15.11b Lower p, and does the loop close? (2026-09-30)

`TALLY_CAPT_LAT` now counts envelope captures only (ball captures sit at the ball
radius; they stay in `TALLY_HIT_LAT`).

**Lower p (proton, λ = 0.3, maze vs isotropic control, same p):**

| p | captures/step | fresh polar intake ratio | equatorial ratio |
|---|---|---|---|
| 1 | 458 | 1.33 | 0.82 |
| 0.3 | 171 | 1.09 | 0.94 |
| 0.1 | 74 | 1.045 | 0.98 |

- The field still steers intake at every p.
- The steering is a **linear response** to emission strength: about 0.06% more
  polar intake per captured-and-re-emitted photon per step. There is no
  threshold and no runaway.
- With uniform capture probability, steering changes *where* charge arrives
  but not *how much* is taken in, so it cannot feed back.

**Loop closure with polar-only intake.** Intake only near the poles, so
poleward steering should raise intake, emission, and then steering:

| aperture | λ | envelope captures/step, maze vs isotropic | change |
|---|---|---|---|
| 1/3 | 0.3 | 34.5 vs 33.9 | +1.8% |
| 0.6 | 0.3 | 114.0 vs 108.9 | +4.7% |
| 1/3 | 1 | 36.0 vs 35.4 | +1.7% |
| 0.6 | 1 | 125.0 vs 121.3 | +3.1% |

- **The loop closes:** the disc raises its own polar intake. Statistical noise
  is ~0.3%.
- **Its gain is small:** roughly 0.04 extra intake per extra emitted photon,
  far below the ~1 needed for Mathis's "accelerates … up to some calculable
  maximum".
- The gain rises with aperture and with density. Raising it means more
  emitted-vs-incoming collisions just outside the envelope: a denser near
  field, stronger emission per body, or a longer interaction zone.

#### 15.11c The holding tank does not raise emission per body (2026-09-30)

**Fix:** held photons all sat in the core cell (~30k at a hold of one level-12
cycle). The collide kernel drew O(count²) candidate pairs there and rejected
every one, so the step ran ~10× slower and the cell overflowed MAX_CELL.
- Held photons now go to a **parking bucket** (index `num_cells`). Dead
  photons go to `num_cells + 1`.
- Compaction keeps the parking bucket, and collide never visits it.
- The sorter has `num_cells + 2` buckets.
- Suite: 210 passed / 0 failed.

**Test** (proton, λ = 0.3, hold = π, i.e. one level-12 spin cycle, and
hold = 10):

| run | captures/step | releases/step | held | fresh polar ratio (no field = 1) |
|---|---|---|---|---|
| env p = 1, no tank (§15.11) | 458 | — | 0 | 1.456 |
| env p = 1, hold π | 448 | 444 | 28k | 1.492 |
| env p = 1, hold 10 | 439 | 556* | 90k | 1.500 |

\* Still draining warm-up captures. The window starts at step 300, and a hold
of 10 is 200 steps.

**Polar-only intake (aperture 0.6), maze vs isotropic envelope captures:**
- no tank: +4.7%;
- hold π: +5.2% (114.5 vs 108.8);
- hold 10: +5.6% (114.2 vs 108.1).

**Reading:**
- **In steady state a tank re-emits exactly what it takes in**, so it cannot
  raise emission per body. Measured: releases = captures (444 vs 448).
- It removes rate × τ photons from the gas (28k–90k here).
- It nudges the steering up by ~2–3% (polar 1.456 → 1.49) and the loop gain
  by under a point. That is consistent with fewer gas photons near the core,
  not with stronger emission.
- **Emission per body is capped by what arrives.** Raising it needs more
  arriving charge: a denser ambient, or a neighbour's emission. Holding it
  longer doesn't help.
- The tank's real use is the capacity / held-charge targets (EMC saturation,
  15 m_e held). It is not a loop amplifier.

### 15.12 Short stacks: does lining up nucleons raise the loop gain? (2026-09-30)

**Setup:**
- `report_stack`: 4 maze bodies in a row, separation 2.2 (envelopes 0.2
  apart).
- Axial (along the poles) or side (equatorial plane), λ = 0.3, polar-only
  intake (aperture 0.6), depth 0, no tank.
- Loop gain per body = envelope captures with the maze disc vs the
  isotropic-emission control (single proton: +4.7%, §15.11b).

| stack | bodies (in row order) | gain per body |
|---|---|---|
| p-n-p-n axial | p (end), n, **p (between two n)**, n (end) | +4.5%, −11.6%, **+15.2%**, −9.3% |
| p-p-p-p axial | | +0.4%, −3.0%, −3.3%, +1.2% |
| p-n-p-n side | | +5.8%, −7.1%, +6.1%, −7.5% |
| p-p-p-p side | | +4.5%, +3.6%, +3.3%, +4.7% |

- **Alternating p-n axially is the one arrangement that amplifies the
  loop.** The proton between two neutrons triples its single-body gain
  (+15%), and 20% of its intake comes from neighbours. The neutron jets
  pour into its polar intake: Mathis's nuclear axis with neutrons feeding
  protons.
- **Protons pole-to-pole shade each other** (−3%). Side-by-side rows behave
  like isolated bodies.
- **Neutrons lose intake (−7% to −12%) under polar-only intake**, because
  their own polar jets block their polar intake. This is the same tension as
  §15.11: our neutron releases charge where Mathis's neutron takes it in.
- Even the best case stays far below a self-sustaining gain (~1).

### 15.13 Denser local field: the recycling loop starts to self-confine (2026-09-30)

**User's lead:** the local charge field is probably far denser than we have
run.

**Mathis on density:**
- An average charge-field density of 1.54e-29 kg/m³, about 5.6e7 photons/m³
  (photon3.pdf). This is a universal average.
- Heat = photon density (heat.html), with no room-temperature figure.
- At proton scale the field is "quite full" (fine3.pdf).
- Converting to a proton-scale mean free path hinges on the photon's
  effective size, which his figures leave open by ~17 orders of magnitude:
  the 1e-24 m radius vs the c²-scaled field size. So λ stays a free
  parameter.
- The usable target is that a proton's near field is dominated by its own
  emission. The measured proxy is the **recaptured share** r of intake: each
  photon recycles 1/(1 − r) times on average before escaping.

**Run:** single proton, envelope intake p = 1, maze disc vs isotropic control.
Finer cells (h = 0.1) and n_a = 8000 (`CFM_C_H` added to `report_intake`).

| λ | captures/step, disc / control | recaptured share, disc / control | recycles per photon, disc | fresh polar intake, disc ÷ control |
|---|---|---|---|---|
| 0.3 | 951 / 877 | 0.298 / 0.266 | 1.42 | 1.28 |
| 0.15 | **957 / 776** | **0.533 / 0.440** | **2.14** | 1.34 |

- **The disc now holds its charge better than isotropic emission.**
  - At λ = 0.15 the disc body recycles each photon 2.1 times on average,
    against 1.8 for the control.
  - It takes in 23% more in total, because its core ball catches back about
    200 of its own emitted photons per step (21 for the control).
  - The concentrated disc scatters more of its outflow back onto the core.
  - This is the first sign of the self-confinement Mathis describes: trapped
    photons held in by the boundary (fine3.pdf).
- The recaptured share rose 0.03 → 0.09 → 0.27/0.30 → 0.53 as λ went 3 → 1
  → 0.3 → 0.15. A near field dominated by emission (r → 1) looks reachable at
  higher density.
- Polar steering holds at about 1.3× the control. The h = 0.1 run gave 1.28
  vs 1.33 at h = 0.25, so it is resolution-stable.
- **λ = 0.1 did not finish.** At h = 0.1 the cells at the core overflow
  MAX_CELL already at λ = 0.15 (warned once), and at 0.1 the collide pass
  stalled for over an hour. It was stopped.
- **Going denser needs a better grid first:** finer cells near the body, or
  a cap on per-cell candidate counts.

## 12. Out of scope for this step

Moving or rotating bodies (atom-mode coupling), the Godot live view, >44M photons,
spatially varying `photon_fraction`, electron bodies in stacks.
