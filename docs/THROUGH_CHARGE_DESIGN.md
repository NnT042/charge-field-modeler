# Through-Charge Flow Model — Design (Session 32)

**Status: Phases A and B implemented. Phase C is design-locked but not
built.**

## Session-32 addendum: Phase B tension SIGN BUG (found and fixed)

The original `flow_tension_pair` signed the two stream directions
asymmetrically: links where the pair's `j` end was the EMITTER pushed
the pair APART instead of pulling it together. Consequences, all
verified by `report_plug_retention` before/after:

- Symmetric links (core alpha ↔ core alpha, both feeding each other)
  partially SELF-CANCELLED — the earlier tension sweeps looked monotone
  but measured roughly half the intended glue, and the "tension supplies
  part of the glue" Phase-B outcome understated the fixed-sign strength.
- Asymmetric links were live wrong: the plugged-carbon socket stream
  BLEW THE PLUG NEUTRON OFF (`plugN↔socket` signed tension +0.12
  repulsive at rest; −0.27 attractive after the fix). This was the
  dominant plug-retention killer — electron.pdf's "charge wind" glue
  was acting as a leaf-blower.
- Post-fix: plugs hug the stack axis in RigidAlpha trajectories; carbon
  worst-drift halved (218–409% → 109–127% on 3 of 4 seeds).

Re-earn outcome (fixed sign): three successive chooser winners FAILED
the 8-seed × 1M gate — (0.85, 2.0) 2/8-fail (325k/575k), the shipped
(0.85, 5.0) 1/8-fail (700k; its pre-fix 8/8 was earned against the
buggy force), and (0.95, 2.0) 4/8-fail (250k–600k). All three joined
`LONG_HORIZON_VETO`; the failure pattern is ambient-driven (every
low-ambient combo wins short-horizon and loses the long game). **New
defaults: channeling=0.85, nuclear_ambient=10.0** — 16.5% quiet /
15.9% flyby in the chooser, and 8/8 × 1M with every seed converging to
the same cold compressed attractor (spacings 3.24/6.48, KE ≈ 5.6e-4,
perfect alignment).
`flow_tension` stays 1.0 / `flow_align` 0.5: tension=2.0 was also gated
and failed (7/8, k=5 at 375k), and the post-fix sweep
(`report_plug_tension_sweep`) shows tri_alpha stable 10–13% at ALL
tensions 0.5–8 while plugged carbon fails at ALL of them (96–825%, KE
pumping to ~4000) — plug retention is NOT a tension-magnitude problem.
The remaining blocker is an energy pump specific to plug configurations
(`report_plug_energy_ablation`): tri_alpha ends 300k runs at KE ≈ 0
while carbon ends hot in EVERY single-term ablation; best retention at
`couplings.torque = 0` (61%/125% vs baseline 116%/109%) and much worse
at `flow_align = 0` (452–482%) → leading hypothesis: the charge
"equator toward charge" torque and the flow-align torque fight over the
plug proton's edge-on orientation, and the loser precesses forever,
pumping energy through the 8-step-stale flow solve. That torque
contract is the next work item, ahead of Phase C2.

Phase B outcome (session 32): stream tension + channel-alignment
stiffness live on same-group pairs (`AtomCore::flow_tension = 1.0`,
`flow_align = 0.5`; `flow_tension_pair` in charge_flow.rs, hooked into
`compute_forces`). The sweep winner moved to channeling=0.9,
nuclear_ambient=2.0 with 7.8% quiet / 7.8% flyby — the tension supplies
part of the glue the higher ambient had been compensating for, and the
balanced margins beat the pre-tension 12.3%/11.9%. `pair_force_breakdown`
gained a 9th element (signed tension along d_hat).

Phase B findings on the FreeNucleon lone alpha (the stretch goal — NOT
yet achieved, by design honesty):

- The dominant "dissolution" was a GEOMETRY bug, not physics: the two
  rest-pose posts (±0.7, disc-shaped contact radius 1.0) were
  interpenetrating, and FreeNucleon's first step fired them out on a
  58-unit compressed contact spring. Fixed via the neutron-as-rod
  same-group contact rule (phos.pdf "neutrons are 1D"), keeping POST_R
  at the earned 0.7 (widening it to 1.05 wrecked the RigidAlpha flyby
  envelope — 1/11 robust — and was reverted).
- With contact fixed: posts are still net-repelled (the neutron's REAL
  measured equatorial emission leak, ~0.53 vs 1.23 polar peak, drives
  un-channeled post-post charge repulsion 1.47 at rest), and once posts
  leave, the protons collapse into contact — oxygen.pdf's "neutrons
  keep the protons apart" is load-bearing.
- Tension response is monotone and real (free-alpha KE 5.9 → 0.9 across
  tension 0 → 16) but link tension alone cannot confine: at 16 the
  posts nearly stall while the over-pulled PROTONS eject axially.
- Conclusion: the missing piece is AMBIENT CONFINEMENT — the external
  field pressing inward on the composite (haf.pdf side-charge, salt.pdf
  pressure), which is the SAME mechanism as Phase C's stray-immunity
  bulk occlusion. FreeNucleon cohesion therefore lands with Phase C,
  not with more link tension.

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
