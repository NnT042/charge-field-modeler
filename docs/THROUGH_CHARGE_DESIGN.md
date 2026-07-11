# Through-Charge Flow Model — Design (Session 32)

**Status: Phases A and B implemented (B.2 as off-default knobs).
Phase C is design-locked but not built.**

## Session-34 addendum: Phase B.2, the starved network, and the no-starve fix

Phase B.2 (throughput-scaled forces) landed as two factored runtime
knobs, both shipping 0.0: `flow_suction` scales the intake pull +
channeling force by the acting body's live `FlowState::mult`
(`1 + gain·(mult−1)`, clamped ≥ 0 — phos.pdf: an unfed socket "doesn't
have much pull, or suction"; diatom.pdf: a stoppered pole's vortex
"mostly dries up"), and `flow_emit_scale` applies the same form to the
charge push ("a fed funnel pushes harder"). A free particle has
mult = 1 by construction, so both scalings are exactly transparent
outside a live network. Env overrides for gate runs: CFM_SUCTION,
CFM_EMITSCALE, CFM_GAP, CFM_FLOWEVERY.

**The v1 suction sweep was a clean negative** (5 gains × 2 emit modes
× 4 seeds × 300k, `report_plug_suction_sweep`): carbon retention did
not improve at ANY gain, and the anatomy showed why — the network said
the plug proton was STARVED (mult < 0.75), so at gain ≥ 4 the scaling
zeroed the plug pair's own vortex glue (0.61 of the pull holding the
plug neutron). Root cause in pass 2 of `solve_charge_flow`: ambient
was displaced by the full GEOMETRIC capture demand while the pole was
credited only the competition-shared, cos²-gated delivered stream —
every plugged pole lost more than it gained. Carbon ran at mult
0.53–0.80 everywhere (total intake 7.3 vs 10.0 for the same particles
free): bound structures channeled LESS than free ones, inverting
graphene.pdf's picture of polar plugs as "fans, increasing the charge
streams coming in".

**No-starve rule (the fix):** ambient displacement is capped at
`delivered / ambient` — deut.pdf's leaky hose ("like a hose that
hasn't been screwed in all the way"): an incomplete plug lets ambient
keep seeping in. A pole's intake can never drop below its free-field
ambient; plugging can feed, never starve. Carbon now runs mult ≥ 1.0
everywhere (fed sockets 1.13, stack-end intake protons 1.33, plug
neutrons at capacity with stress 0.44), and the stream tension —
already throughput-scaled by construction — roughly DOUBLED at the
plug sockets for free (plugP↔socket −0.16 → −0.30).

**Re-earn:** the network change tripped `alpha_stays_bound` by design
(all pre-v2 long-horizon evidence was earned against starved flows).
All five contender combos were re-gated 8 seeds × 1M: **(0.85, 5.0) —
the original session-32 default — passes 8/8 all-cold and has the best
surviving short-horizon margins (13.7%/12.5%), so `nuclear_ambient`
returned 10.0 → 5.0.** (0.85, 10.0) is a clean runner-up (8/8
all-cold); (0.9, 5.0) and (0.9, 2.0) collapse 2/8; (0.85, 2.0) stays
bounded but 3/8 seeds heat monotonically to KE ≈ 128 (the pre-collapse
whirl signature). The coherent story: ambient 10.0 had been
compensating for exactly the glue the starved network under-delivered.

**The injector hunt** (`report_plug_energy_ablation`, full matrix
including drag/intake/solve-cadence cases): NO single-term injector —
doppler's asymmetric clamp (0.2–5.0) injects the most heat but is also
the system's main damping (drag=0 collapses KE 382→37 while drift
explodes to 2445–3212%). The pump framing was retired in favor of a
STATIC question: does the plug configuration have a potential well at
all? (`report_plug_well`, displacement scan of the pair with the flow
network converged at each geometry.)

**The hidden 150-unit blast (the real session-33/34 retention
killer):** the pair spawned as two 1-nucleon alphas with the proton's
pole aimed AT its partner. The neutron's real equatorial emission leak
(~0.5) fired ~150 units of charge push down the proton's pole-on
absorption maximum at contact range — invisible in every anatomy
table, which prints only the force on the second body of each pair. No
azimuth escapes it: at hole-sharing range (0.7–1.6) some equator
always blasts some partner at 1/r⁴. Resolution per the papers:
**FUSION** — a p+n in contact sharing a pole IS deuterium (deut.pdf),
and fused units "can't be broken and rearranged" (uf4.pdf), the same
pre-fusion rule that exempts every core alpha's members. `plug_pair()`
now spawns carbon/oxygen plugs as ONE 2-member alpha each (same world
anatomy, honest mass-2 inertia); nitrogen's lone plugs are unchanged.

**The seat:** the old funnel-mouth spawn (2.6 beyond the end proton,
inherited from stacked-alpha pitch) sits 1.2 units up a monotone-
inward axial slope — the rigid pair FELL into the core and scattered
it (the fall was always there; the intra-pair contact bounce used to
disguise it). The extended well scan puts the true equilibrium ≈ 1.37
beyond the end proton → `PLUG_SEAT = 1.4` (ammon.pdf "that seventh
proton is in tight"; the graphene.pdf diagrams nest the plugs in the
hole). At the seat the static landscape is stable in BOTH axes (axial
stiffness ≈ 89/unit; lateral walls ±88 within 0.5 units).

**Open (end of session 34):** translation is statically solved;
dynamics still melt. The pair's free rigid-body rotation re-aims the
±90-unit close-range force structure on every tip (KE pumps to
400–1400; high tension binds the nucleus into a hot disordered droplet
instead of ejecting). Probe knob `plug_orient_lock` (ships false,
CFM_PLUGLOCK) slaves ≤2-member alpha orientation in RigidAlpha per
graphene.pdf ("the proton and neutron ... stayed in line"). The
deeper gap: the model has no energy sink beyond doppler/contact
damping — transients cannot radiate away. A paper-based dissipation
channel is the standing design question for the next session.

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
## Session-33 addendum 2: the plug-retention campaign (pm)

The post-sign-fix energy pump was root-caused NOT to a torque war but
to the plug PAIR: its rest distance (2×PLUG_PAIR_GAP = 1.6) is not a
force equilibrium — net inward ≈ 1.4 (ambient shadow, gravity, intake;
charge repulsion is exactly zero pole-on) with nothing opposing before
the contact wall at r ≈ 0.7. The pair is dropped onto the stiffness-100
contact spring: an undamped collapse-bounce oscillator. Spawning at
contact (`plug_pair_gap = 0.35`) kills the pump (end KE 5–24 vs
17–578).

Mechanisms tried and REJECTED with data (all remain as off-by-default
runtime knobs; reports named):
- tension magnitude (`report_plug_tension_sweep`), gyroscopic spin
  stiffness `gyro_spin` (`report_gyro_sweep` — diverges at small S,
  hot at large S), align-to-stream torque target `align_to_stream`
  (capture-gate version rerouted the network; torque-only version
  interacts badly with the tight gap), close-range ambient saturation
  `ambient_sat_r`, and Phase C1 `ambient_confine` on the pump-killed
  cell (`report_plug_confine` — monotonically worse; the pre-fix "not
  confinement-fixable" verdict stands).
- Factored attribution matrix: `report_plug_matrix` (the serial
  stacking of these "obviously right" fixes was worse than baseline —
  interactions dominate).

With the pump dead the plugs slide off QUIETLY: the socket's static
intake force measures a useless ~0.016 regardless of how much flow the
socket actually channels, while the plug proton feels ~0.85 outward
charge push. Per phos.pdf an unfed socket "doesn't have much pull, or
suction" — and ours is permanently unfed as far as the FORCE model is
concerned. The missing piece is exactly Phase B item 2, designed above
and never built: **throughput-scaled emission/intake** (scale the
static force profiles by the live `FlowState::mult`). That is the next
work item; then re-run the matrix best cell and the full re-earn chain,
and reconsider `plug_pair_gap = 0.35` as default.

Also landed: `quad_alpha` preset + `report_quad_alpha` (haf.pdf free
prediction — a bare 4-stack "can't hold together"; currently documents
that the model does NOT yet dissolve it, pending real side-charge), and
the post-sign-fix radial-post retest (`report_post_anatomy`): Radial
now BEATS Axial in every rigid harness (9.9%/10.5%/10.6% vs 15.9%) but
still cannot self-hold in FreeNucleon (1023%) — the anatomy default
stays Axial until the suction mechanism exists.

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
