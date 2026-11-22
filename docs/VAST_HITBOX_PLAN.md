# Plan: hit-box photon runs on rented GPUs (2026-11-22)

Follows doc §15.59 (tc59). Status: **ready for your approval; nothing rented or run yet.**

## Where things stand

Step 0 ran on our own machines (results in `reference/tc59`, commits 1174623 to 52a5c6d).

- **The pile-up was the empty start, not the hit box.**
  - A wraparound box stays flat at 1.00 with the hit box on, and the escape ball slowly settles.
  - The pre-filled start (`CFM_P_PREFILL=1`) removes the pile-up: flat from the first snapshot.
- **Pre-filled hit-box p-n runs at h 0.25:**

| run | proton push | neutron push | proton take-in /tu | over-full cells |
|---|---|---|---|---|
| hit box only | −114 (empty start: −121) | +4 | 5,030 | none |
| hit box + mass | −5,373 | +155 | 6,150 | a few, at the ball's edge |
| mass only (no hit box) | −9,480 | +450 | 8,900 | none |

- **The lone-proton leftover force is sideways** (~134 across the pair, ~12 along it), so the pair numbers are not swamped as I first said.
- **The remaining over-full cells** (hit box + mass) are in the cell column on the pole axis at the ball's edge, 4–5 from the bodies, where the outgoing beam meets the incoming field. None are at the bodies.
- **The open question is resolution.**
  - With the hit box on, the average photon travels ~0.017 between hits (two level-3 photons ~0.005), but the cells are 0.25 across.
  - Collision partners in one cell can be ~15 collision distances apart, which blurs the push near a body.
  - We can't get cells down to 0.017: that needs ~1.8 billion photons, and the engine caps at 19 million.
  - What we can do is shrink the cells in steps and see whether the answers stop moving.

## The ladder

Keep the ball at R_D 6, so the edge stays 4–5 from the bodies, and shrink the cells. The pre-filled start is used throughout. Settings otherwise as tc59 (tc56 standard, frozen p-n at 2.05, even mix, random inner spins).

| rung | h | n_a | photons per cell | live photons | work vs rung 1 |
|---|---|---|---|---|---|
| 1 | 0.25 | 8,000 | ~125 | 7.2M | 1× |
| 2 | 0.177 | 12,000 | ~66 | 10.9M | ~2× (1.5× photons, 1.4× steps) |
| 3 | 0.125 | 16,000 | ~31 | 14.5M | ~4× (2× photons, 2× steps) |

Rung 3 is the §15.57 "hi-res" setting, already proven on the Vast 5090. The time step stays at dt 0.06 (c·dt < h holds down to h 0.125, and §15.57 ran this way). Rung 1 repeats on the 5090 what we already have, which checks that the rented machine gives our numbers.

**Runs:**

| rung | runs | count |
|---|---|---|
| 1 | p-n hit box only; p-n hit box + mass; H atom hit box + mass, seeds 49344 / 1 / 2 (does the electron really leave?) | 5 |
| 2 | the two p-n runs | 2 |
| 3 | the two p-n runs | 2 |

Each p-n run also prints the lone proton and lone neutron, and the crowded-cell locator (`CFM_P_CELLS=1`) is on in all of them.

**Pass test:**

- the proton push and take-in change by less than ~20% from one rung to the next;
- no over-full cells within 2 of a body;
- the leftover force along the pair axis stays small next to the push.

If the numbers still move at rung 3, the hit-box answer can't be trusted at any size we can run, and the next step is a design change (option B below), not more hardware.

## Cost and time (estimate)

These come from our own timings: hit-box p-n ~4.3 h on the 3060 and ~4.6 h on ISABELLE's 5060. I'm assuming the 5090 is ~3–4× faster, running 4 jobs at once ($0.52/h, §15.57).

| rung | wall clock | cost |
|---|---|---|
| 1 | ~4–5 h | ~$3 |
| 2 | ~5–6 h | ~$3 |
| 3 | ~9–10 h | ~$5 |
| total | ~one day | **~$11–15** |

The autostop script pauses the instance 30 min after the queue ends (12 h cap per start). I'll time rung 1 and re-quote before going on.

## What I need from you

1. **Approve the ladder** (~$11–15).
2. **Start the Vast instance** (61730368, or a fresh 5090) from your account; payment and account steps are yours. Once it's up I bootstrap it over ssh from the pushed `main` (the tc59 code is all pushed), run rung 1, and bring the results back. Destroying the instance stays your button.

## Option B, for later or instead

**B: hold the field's average mean free path at 0.3.** Scale the base cross-section down so the mixed field collides as often as the field we've always used. The hit box then only decides which photons hit harder and more often relative to the rest.

- It runs fine at h 0.25 on our own machines, so no rental is needed.
- It needs a new knob and your sign-off on what it means physically.
- It's the fallback if the ladder doesn't settle.

## Not in this plan

- Why the beam bunches at the ball's edge (it's 4–5 from the bodies; the locator keeps watch).
- The axis-restoring mechanism.
- More than 4 photon levels.
- Changing defaults: everything stays behind the `CFM_P_*` knobs.
