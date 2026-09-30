# Mathis on the spin paths: how charge moves through a baryon

Research notes for the channelling / recycling design (GPU_DSMC_DESIGN §15.8 ended
with: a photon that only bounces off the base ball cannot carry the neutron's
magnetic signature; the photons have to travel *through* the spin structure).
Everything below is paraphrased from the Mathis Library mirror, with source links.
Direct quotes are kept to a few words.

Mathis gives **two descriptions** of what happens inside a baryon. They were
written years apart and he never explicitly joins them. A third, photon-scale
picture fills in the capture step.

---

## A. The emission maze (quark.html, 2008: "overhaul of QCD")

Source: https://milesmathis.com/quark.html

- **Setup.** Assume every baryon emits equally from the surface of its core.
  That emission must travel outward through the four stacked spins (a, x, y, z).
  Some sign combinations let it through and some block it.
- **Method.** He works in one plane. The core emits outward like a pinwheel or
  sprinkler (one dt of emission is an expanding circle).
  - At each spin level the emission must "turn right or left", with right = "+".
  - He considers the levels one at a time, with magnitudes 2, 4, 8. Neither the
    end-over-end x nor the y level can block anything by itself, because the
    particle can't outrun its own emission at c. Blocking first becomes possible
    at the z level (magnitude 8).
  - Worked example: three right turns bring the emission back to where it
    started, so it is trapped.
  - He calls the one-plane version "a slight oversimplification" and gives only
    the result of "the full analysis". **The full 3-D analysis is not in the
    paper.**
- **Result (the 16 baryon states):**

  | class | states | what happens to the emission |
  |---|---|---|
  | proton | +a+x+y+z, −a−x−y−z, −a+x+y+z, +a−x−y−z | escapes with the same direction and chirality it was emitted with |
  | neutron | −a−x−y+z, +a+x+y−z, −a+x+y−z, +a−x−y+z | cannot escape, *and* returns "upside-down" (reversed spin), so its energy cancels |
  | anti-proton | +a−x+y−z, −a+x−y+z, +a+x−y+z, −a−x+y−z | escapes in the same direction but upside-down |
  | anti-neutron | +a−x+y+z, +a+x−y−z, −a−x+y+z, −a+x−y−z | trapped, but not reversed, so it does not fully cancel (mass difference) |

- **Consequences he draws:**
  - The neutron's emission "cannot get beyond its z-spin".
  - The neutron is heavier because it keeps the emission's energy; the baryon
    emission's mass equivalence is ~2.3e-30 kg.
  - "A neutron is always a proton with a reversed outer spin." Since the field
    photons' spins interact directly with the baryon's *outer* spins, the
    neutron's magnetic moment has the opposite sign.
  - From the μp/μn difference, the proton supplies about **1/3 of the charge
    field at its own surface**; 2/3 is ambient.
  - Beta decay is a positron hit that flips the neutron's z-spin (weak2.html
    repeats this).

**Note for our model.** In this paper the neutron differs *only in the outer z*.
The adopted "z flipped in every tier" neutron (§15.7) is our extension of the rule
down the tiers, not something Mathis states.

---

## B. The recycling channel (later papers: pole, aurora, magmom, neutron, rhrule, venus2, graphene, corio)

**The engine.** A spinning sphere that is permeable to smaller particles has more
angular momentum at the equator than at the poles. That creates a net outward
vector at the equator and a net inward one at the poles, a "charge engine" with
photons as the working gas. The larger body sorts the field: photons go in the
south pole, antiphotons in the north.
- https://milesmathis.com/pole.pdf

**Where the push comes from (aurora.pdf).**
- The outward (centrifugal) push comes from the whole spin *axis*, not the
  centre, so incoming charge starts being pushed back out as soon as it enters a
  pole. That is why charge is emitted at all latitudes, not only the equator.
- Intake is focused at the poles; outflow is spread over the larger equatorial
  area.
- **Spin speed sets the shape.** A fast spinner compresses outgoing charge into a
  disk "as we see when modeling the proton": all emission below 30° and no
  crossover boundary (charge-pause). Slow spinners (planets) spread emission and
  push the crossover toward the poles.
- https://milesmathis.com/aurora.pdf

**The ±30° exits (magmom.pdf).**
- Photon charge exits at +30°, antiphoton charge at −30°; averaged, they sit at
  the equator.
- Proton and neutron emissions differ by 90° (equator vs pole). He treats that as
  half the photon/antiphoton spin difference, which gives the ~1.5 moment ratio
  (refined to 1.458).
- Also: the proton must spin a little faster than the neutron (it doesn't hold
  charge as long), so it recycles a bit more charge per second.
- https://milesmathis.com/magmom.pdf

**The neutron (neutron.pdf).**
- It pulls charge in at the poles like the proton, but its spin stack blocks the
  normal equatorial exit. Charge is "funneled back" to the core / central spin
  axis and can only leave by a pole: **in one pole, out the opposite one.**
- It is not a straight-through tube. The photons "go in, rattle around"
  following the four spin paths, then exit in an orderly way.
- As they leave, they collide with incoming photons and the spins cancel.
  - In a balanced field the cancellation is near-total, so the neutron is
    non-magnetic, but it still has a strong polar (sub-electrical) field, since
    collisions affect only spin, not c.
  - In Earth's 2:1 field the cancellation is partial.
- In the nucleus, neutrons sit with their axes at 90° to the charge channel and
  act as "stoppers". Protons are "spin fans" that channel *in the plane of their
  spin*.
- https://milesmathis.com/neutron.pdf
- graphene.pdf states it compactly: protons channel pole to equator, neutrons
  pole to pole, and neutrons channel less, by the ratio of their magnetic
  moments. https://milesmathis.com/graphene.pdf

**Through charge (venus2, rhrule, core2, salt).**
- A minority of photons, entering very near the axis on the right trajectory,
  pass straight pole to pole. This matters more at nuclear scale, where the
  nucleus actively pushes charge along the axis.
- Entries near the hole's centre tend to pass through; off-centre or angled
  entries are forced by centrifugal effects into the equatorial "whirlpool".
- https://milesmathis.com/venus2.pdf, https://milesmathis.com/rhrule.pdf,
  https://milesmathis.com/salt.pdf

**Individual photons vs the average path (corio.html).**
- Single photons don't follow the drawn curve. Some pass straight through; some
  hit antiphotons head-on and lose spin.
- The many collisions sum into a curve from the pole to just off the equator,
  which the median photon follows. Each collision also pushes the body.
- https://milesmathis.com/corio.html

**Collisions inside the channel.**
- Recycling compresses charge, and photons collide edge to edge in the channel
  (push2.pdf).
- Opposite streams meeting head-on spin each other *up* ("opposite twice"), while
  same-direction photon/antiphoton pairs spin each other down (rhrule.pdf). This is
  the §13 RS3 / twist-grip result we already built.

**The pole "hole".**
- He deduces rather than observes that the proton is porous to photons but not
  electrons: the proton is ~6e9 times the photon's size but only 1821 times the
  electron's.
- The pole is where photons tend to enter, because angular momentum is lowest
  there. https://milesmathis.com/nuclear.pdf, https://milesmathis.com/bohr3.pdf

---

## C. The capture step at photon scale (photon.html, galcharge, bottom)

- A photon with a z-spin (orthogonal to its motion) acts like a "scoop or an
  intake valve". Smaller photons are temporarily absorbed, and once it has
  taken its fill they spill back out. That exhaust is charge. This is how a
  large photon becomes an electron.
  https://milesmathis.com/photon.html
- Photons coming into a body at the poles gain energy by stacking on more outer
  spins. https://milesmathis.com/galcharge.pdf
- The proton's z-spin is already at c. The in-and-out charge stream acts as an
  extension of the particle beyond the z-spin "edge".
  https://milesmathis.com/bottom.pdf

---

## What is and isn't specified

| question | Mathis gives | not given |
|---|---|---|
| where charge enters | the poles; photons south, antiphotons north | the entry cone width |
| where it leaves (proton) | ±30°, compressed into a disk for a fast spinner | the path from pole to exit |
| where it leaves (neutron) | the opposite pole, after following the spin paths | what "rattle around" means level by level |
| what decides p vs n | the sign combination of the four spins (maze); the neutron = the proton with z reversed | the full 3-D maze rule (only the one-plane example) |
| magnetic signature | exit spin: same (p), reversed (n); outer spins interact with field spins | how a level reverses a photon's spin mechanically |
| how much is recycled | ~1/3 of the field at the proton's surface; neutron less, by the moment ratio | a rate law |
| inner tiers (levels 1–8) | spins stack as fractal tiers (photon2.html) | whether the p/n pattern repeats in inner tiers |

The two models run in **opposite directions**:
- **A** is emission *from the core outward* through the levels.
- **B** is ambient charge *from the poles inward, then out*.

One way to reconcile them: incoming charge reaches the core region (B, "funneled
back to the core"), then leaves through the levels by the maze rule (A). Mathis
implies this (neutron.pdf points back to the maze paper for how the blocking
works) but never states it as one mechanism.

---

## Candidate mechanical representations (for discussion)

1. **Maze first, as a CPU check.**
   - Implement the level-by-level turn rule in 3-D: at each level the photon's
     direction turns ±90° about that level's axis, with the sign from the state.
   - See whether *any* reasonable 3-D version reproduces the 16-state table:
     escape vs trapped, and upright vs upside-down.
   - Cheap, falsifiable, and it would settle what "the full analysis" must have
     been before we build anything on the GPU.
   - It also tests the z-flip-every-tier neutron directly: apply the rule to the
     inner tiers and see what comes out.
2. **Carried by the levels (permeable-body transport).**
   - Inside the envelope, photons are dragged by the local matter velocity of the
     level region they are in (the §14 shell kinematics, V_k(x)), instead of
     bouncing off a ball.
   - The centrifugal push from the whole axis, and the 30° disk for a fast
     spinner, should then emerge. p vs n enters through the sign pattern of the
     velocities.
   - This is a flow model (Mathis B).
3. **Scoop and release (capture step).**
   - A level's z-spin captures photons temporarily (a dwell time), and they are
     released along the maze exit direction with the exit spin (upright or
     upside-down).
   - This is closest to "go in, rattle around, exit in an orderly manner", and
     gives the charge-pause a real source.

Option 1 is the natural first step. It is a pure-logic test of Mathis's own
classification, and its outcome decides the rule that options 2 and 3 would
carry out mechanically.

---

## Results: the maze on the CPU (2026-09-30, `rust/src/maze.rs`)

**Reading the table first.** Mathis's 16-state table depends only on the relative
signs of x, y and z; the sign of a never matters.

| class | odd level | escapes | returns upright |
|---|---|---|---|
| proton | none (x = y = z) | yes | yes |
| neutron | z | no | no |
| anti-proton | y | yes | no |
| anti-neutron | x | no | yes |

So **escape ⇔ x = z** and **upright ⇔ y = z**. The outer z-spin gates the
*direction* against x and the *spin* against y. (`table_structure` test.)

**Search.**
- The family:
  - turn angle 45–180°;
  - each level's axis ∈ {X, Y, Z};
  - lab-fixed vs carried axes; level order a→z vs z→a;
  - start direction and photon spin;
  - escape criteria vs d0, and upright criteria vs s0, the pole or d0, fixed in
    advance, with upright always the positive side.
- Rules whose results sit within 1e-6 of a threshold are rejected; an earlier
  pass was fooled by exactly that rounding.

**Rigid turning (spin carried along with the path): 0 of 907,200 reproduce the
table.**
- The escape column is easy: 1,236 rules get it right.
- The upright column is not: among escape-perfect rules the best is 8/16, i.e.
  chance.
- A rigid turn preserves helicity, so a trapped photon always comes back the
  same way. Neutron (trapped, upside-down) and anti-neutron (trapped, upright)
  cannot be told apart.

**Contact spin rules reproduce it.** Each level also acts on the photon's spin
component along its axis. Three rules were tried:
- *Set*: the level imprints its sense.
- *Flip*: a mismatched spin component is reversed.
- *Rough*: the §13 twist grip at full grip.

Of 3.6M rules, 24 score 16/16. With **Spin mode's own axes** (a and z about Z, x
about X, y about Y; spin_stack `rotation_axis`), 20 rules pass, all with:
- 90° turns about lab-fixed axes, met in the order a → x → y → z;
- emission starting at **45° between the pole and the x axis**, ±(X + Z)/√2;
- escape = "not turned back"; *Set* or *Flip*;
- several initial spins work, including a helical photon (s0 = −d0).

This is locked in by `spin_mode_maze_reproduces_mathis_table`, which checks that
the rigid version of the same rule fails.
- **Reading:** the geometry we already use contains Mathis's classification, but
  only if photons interact with each level *by contact* (spin changed relative to
  path). This matches the §15.8 finding: bouncing off a rigid body can't carry the
  p/n magnetic signature.

**Twelve levels (`report_maze_twelve_levels`).** Running the passing rules as one
long maze through all 12 Spin-mode levels breaks the classification:
- the four proton states come out as four different classes;
- the Mathis neutron and the z-flip-every-tier neutron give identical results
  (flipping levels 4 and 8 changes nothing under these rules).

So the maze is a **single-tier** construct. Emission starts at the surface of the
"particle" made by the inner tiers (our base ball) and passes the four baryon
levels, as Mathis describes it. How tiers would chain is open. One option: each
tier emits from its own core, and a tier's output becomes the next tier's input
with a fresh 45° start.

**Next (for discussion):** make the maze mechanical inside the DSMC ball model.
A photon that reaches the base ball is carried out through the four baryon
levels by this rule (turns + contact spin) instead of bouncing:
- escapes are emitted with the maze's exit direction and spin;
- trapped ones are held and released at a pole (neutron.pdf), upright or
  upside-down per the rule.

That is option 3 ("scoop and release") with the rule now fixed by Mathis's own
table.
