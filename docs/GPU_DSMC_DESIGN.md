# GPU DSMC Photon Engine — Design (Phase 2, step 3)

Status: design, 2026-09-27. Builds on the wgpu harness (`rust/src/gpu/`, commit 9ae82b3).
Owner decisions still open are marked **[DECISION]**; the engine implements every
option so the choice is a config flag, not a rewrite.

---

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

## 12. Out of scope for this step

Moving or rotating bodies (atom-mode coupling), the Godot live view, >44M photons,
spatially varying `photon_fraction`, electron bodies in stacks.
