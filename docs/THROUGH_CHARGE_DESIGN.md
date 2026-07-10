# Through-Charge Flow Model — Design (Session 32)

**Status: Phase A implemented (flow network + diagnostics, no force
changes). Phases B/C are design-locked but not built.**

The user's top priority after the whirl-instability fix: the force model
treats emission/absorption as static geometric profiles, so every effect
that depends on charge being *conducted through* a body and re-emitted
elsewhere is invisible to it. That gap is what made the radial-post
anatomy experiment a false negative (the deut.pdf regulator is a flow
effect), it is why FreeNucleon alphas dissolve (posts are held by the
charge flowing through them, not by profile geometry), and it blocks the
carousel funnel drive and stray-flier immunity.

## Paper basis (all via search_mathis, session 32)

1. **The circuit.** Charge enters at the poles and exits at the equator;
   the pole→pole component is "through charge" (graphene.pdf,
   nuclear.pdf Faraday-disk analogy). Both directions flow (photons
   S→N, antiphotons N→S); on Earth the split is 2/3 photons, so S→N
   dominates (graphene.pdf, ammon.pdf).
2. **Protons vs neutrons.** "Protons channel charge pole to equator,
   while neutrons channel pole to pole" — protons are 2D "fans",
   neutrons 1D "lightning rods" (phos.pdf, graphene.pdf, neutron.pdf).
   Neutrons channel at **67%** of a proton (haf.pdf; matches the
   magnetic-moment ratio 1.913/2.793).
3. **Posts divert.** In-alpha neutrons sit with poles along the alpha
   axis and "divert part of the charge that the protons would otherwise
   have channeled equatorially ... back to the axial channel"
   (atmo2.pdf). This is the quantitative role of the posts, and it
   supports the Axial anatomy (the session-32 A/B verdict).
4. **Binding is the stream.** "The vertical charge streams coming up and
   down the poles normally tend to tie them [stacked alphas] together
   pretty tightly" (graphene.pdf). Breaking alignment breaks the
   channel; the channel resists (strong.html, bb2.pdf — already the
   basis of the channeling factor).
5. **Conservation and stress.** A chain channels no more than its
   narrowest cross-section; excess must be recycled out laterally where
   parallel alphas allow it (salt.pdf). A nucleus channeling more than
   its configuration supports is stressed and reconfigures/bonds to
   disperse the excess (graphene.pdf).
6. **Intruder repulsion.** "The strong charge stream exiting the
   carousel level drives off all intruders, including ambient photons
   and free electrons" (auger.pdf) — the stray-immunity mechanism.
7. **Capacities.** An alpha pulls in ~2 proton-units of charge
   (graphene.pdf); polar neutrons raise tritium's through-charge to
   1.37 (deut.pdf); bare stacks channel weakly without polar plug
   "fans", and ≥4 stacked alphas fail because side-charge overwhelms
   the channel (haf.pdf) — a free prediction to test.

## Model

### State (per particle, recomputed every step)

Each particle carries a `FlowState`:

- `intake[2]` — charge/sec entering at pole− (south, index 0) and pole+
  (north, index 1).
- `through` — pole→pole component conducted along the axis.
- `lateral` — pole→equator component emitted at the disc.
- `mult` — output multiplier = total output / free-field baseline; what
  Phase B scales forces with.

Unsigned scalar flow in v1 (no photon/antiphoton split yet; the split
only sets *which* pole dominates, not the binding topology). Free-field
baseline: a lone proton recycles 1.0 unit/sec (0.5 per pole), neutron
0.67, electron 1/1836.

### Ports and routing

- Proton: intake at both poles; output = equatorial disc (`lateral`)
  plus a pass-through fraction along the axis (`through`). Routing
  fraction is configuration-driven, not fixed: output that a plugged
  axial neighbor can accept leaves axially, the remainder exits at the
  disc (salt.pdf: "depending on the circumstances, it can do either").
- Neutron: intake at both poles, output at the opposite pole only,
  capacity 0.67. In-alpha posts additionally *divert*: they capture a
  fraction of their protons' disc output (geometric coupling of the
  disc plane to the post pole) and feed it back to the axial channel
  (atmo2.pdf).
- Electron: negligible conduit (1/1836); riders damp a pole's flow
  (existing `ELECTRON_FLOW_DAMP` VFX constant becomes physical here).

### Coupling (who feeds whom)

Directed port-to-port coupling `g(i.out → j.in) ∈ [0,1]`: the fraction
of i's output at a port captured by j's intake port. Product of:

- distance falloff: smoothstep 1 → 0 over `NUCLEON_PITCH →
  CHANNEL_TAIL·NUCLEON_PITCH` (reuses the channeling geometry — the
  funnel mouth IS the capture cross-section);
- alignment: `cos²` of the emitting port's direction vs the line to the
  receiving port, times `cos²` of the receiving pole vs that line
  (pole-on capture = 1, side-on = 0);
- competition: a port's captures are normalized so no more than 100% of
  an output is consumed, nearest-strongest first.

### Solve

One relaxation sweep per physics step, seeded with last step's values
(Gauss–Seidel; flows settle over tens of steps ≈ instantly at dt=5e-4,
and lag is physical — charge takes time to re-route):

```
intake_j = ambient_exposure_j + Σ_i g(i→j) · output_i(prev)
through/lateral routing per particle kind and plug context
output_j = min(intake_j, capacity_j)        // conservation
stress_j = intake_j − capacity_j (≥ 0)      // graphene.pdf overload
mult_j   = output_j / baseline_j
```

Ambient exposure per pole is reduced by occlusion (existing machinery)
and by being plugged into a neighbor (a plugged pole receives the
neighbor's stream *instead of* ambient, not in addition).

### Phase B — forces (design-locked, not yet built)

1. **Stream tension**: each directed link carries momentum flux
   ∝ `g·output`; interrupting it costs pressure → pair attraction along
   the link, plus **alignment stiffness**: force/torque down the
   gradient of `g` (the cos² terms), pulling plugged ports back into
   line. This is the restoring stiffness the whirl saga showed the
   stack lacks, and the deut.pdf post regulator: a disc dip raises
   `g(disc→post)`, the diverted flow's reaction nudges the disc back.
2. **Throughput-scaled emission**: existing charge push/stream terms
   scale with the emitter's live `mult` (a fed funnel pushes harder
   than a starved one).
3. **Output push on third parties** (Phase C): a body standing in a
   port's output stream feels `∝ output · profile(θ) / r²` — auger.pdf
   intruder repulsion = stray-flier immunity. Carousel funnel drive =
   the tangential component of the core disc's output captured by a
   carousel alpha's axial hole.

Re-earn discipline: every phase re-runs the full battery (quiet sweep,
flyby, 8-seed transient, 1M-step long-horizon, molecular locked tests)
before its defaults ship.

### Phase A deliverables (this session)

- `FlowState` on `SimParticle` + `solve_charge_flow` sweep in `step`.
- No force changes. Diagnostics only:
  - `report_flow_network` (ignored test): print the flow table for
    alpha / carbon / neon / argon — verify the carbon stack carries a
    coherent axial through-stream, posts divert disc output, carousel
    alphas capture the core's lateral output, totals conserve.
  - `AtomSim::get_flow_summary(id)` for HUD/inspection later.
- Sanity tests: conservation (Σ output ≤ Σ intake + ambient), free
  proton mult = 1, plugged stack through-flow > lone alpha's.
