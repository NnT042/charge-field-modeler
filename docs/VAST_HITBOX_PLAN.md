# Plan: hit-box photon runs on rented GPUs (draft, 2026-11-21)

Follows doc §15.59 (tc59). Status: **draft for the user to approve; nothing has been run or rented yet.**

## Why this plan exists

With the hit box on (`CFM_P_SIZE=1`, collision size = outer spin radius), the h 0.25 runs failed:

- **Size only (ISABELLE, 4.6 h):** the p-n push was −121, but a lone proton, which should feel no net force, showed about 150. The grid overflowed 298 times.
- **Size + mass (ISABELLE):** the photons piled up until up to 823 cells held more than the 2,048-photon cell limit. The lone-body leftover force reached ~60,000. After 100 min the GPU dropped out.
- **Field box with no bodies, 24 time units:** the density in the middle was already ~11× the density at the edge, where before it was flat.

## Step 0 first (on our own machines, cheap): what is the pile-up?

**Rented time can't fix this until we know the cause.** If the pile-up is a model artifact, more resolution only fills the bigger GPU faster.

**Leading suspect: the escape boundary.** The ball is fed from its edge as if the outside field streamed freely. With the hit box on, the inside is far more collisional:

- the mean free path for two level-3 photons is λ/64 ≈ 0.005;
- the mix average is ≈ λ/18 ≈ 0.017, against 0.3 before.

So photons arrive faster than they can random-walk back out, and the inward streams meet and stall in the middle. A real field would be just as collisional outside the ball as inside it, with no inward stream.

**Tests, all at h 0.25, all on ESTELLE:**

1. **Long box run:** the no-body field box with the hit box on, run for 200 time units instead of 24, logging density by radius every 10 time units. If the core keeps climbing, it's a runaway, not a steady state.
2. **Periodic box:** the same mix in the periodic box (wraparound, no inflow or escape). If the density stays flat there, the boundary is the cause.
3. **Speed by radius:** add photon linear speed by radius to the box report. If photons slow down in the middle, part of the pile-up is a traffic jam (density ∝ 1/speed), not just trapping.

Time: under an hour each. No new physics, only reporting.

## Decision for the user (after step 0)

**How collisional should a field of bigger photons be?** Right now the knob keeps the photon *count* per volume (n_a 8000) and lets the bigger hit box raise the collision rate up to 64×. Three options:

- **A. Keep it** (bigger photons = more collisions at the same count). This is the literal reading of the hit box and of arp.pdf, which says larger photons collide more. It needs cells smaller than the mean free path (h ≲ 0.017) with ~10+ photons per cell. That is ~1.8 billion photons in today's R_D 6 ball, about 200 GB of photon records, and the engine's 2 GiB storage binding caps us at 19M photons. **Not reachable by rented hardware** without a new engine design (photon weighting or a continuum outer zone).
- **B. Keep the field's average mean free path at 0.3.** Scale the base cross-section down so the mix-averaged collision rate matches the field we've used all along; the hit box then only sets *who* hits whom harder and more often relative to the rest. This runs at h 0.25 to 0.125 on hardware we have or can rent. It is a knob, labelled as one, and needs your sign-off.
- **C. In between.** Pick a field mean free path (say 0.1) and run a resolution ladder to see whether the answers settle.

**Boundary fix, needed for A and C, and for B if step 0 still shows a pile-up.** Feed the ball from a reservoir held at the field density, or give the ball an outer buffer shell where the field relaxes before it reaches the bodies. This is a modelling choice, also yours.

## Rented runs (after the decisions above)

Hardware: one RTX 5090 on Vast (§15.57: $0.52/h, 4 jobs at once, same 2 GiB binding limit as ours, so the 19M-photon cap stands). Scripts: `scripts/vast/bootstrap.sh`, `autostop.sh` (pauses 30 min after the queue ends, 12 h cap). Destroying the instance stays your button.

**Resolution ladder.** The same runs at three cell sizes, raising n_a so each cell keeps ~125 photons, and shrinking R_D where needed to stay under 19M. The time step shrinks with h (c·dt < h), so finer rungs also need more steps:

| h | n_a | R_D | photons (approx.) | relative work |
|---|---|---|---|---|
| 0.25 | 8,000 | 6 | 11.6M | 1× (today) |
| 0.177 | 22,600 | 4.5 | 13.8M | ~1.7× (1.2× photons, 1.4× steps) |
| 0.125 | 64,000 | 3.2 | 14.0M | ~2.4× (1.2× photons, 2× steps) |

R_D 3.2 leaves only ~1.2 of field past the p-n pair, so h 0.125 is for the lone proton only.

**Runs per rung, with the hit box on and the mix chosen in the decision:**

- lone proton (check: leftover force ≲ 10% of the pair signal);
- frozen p-n at 2.05 (push, take-in rate, spin torque);
- H atom (proton plus free electron), 3 seeds, to settle whether the heavy mix really pushes the electron out (tc59 saw it on one seed).

**Pass/fail per rung:**

- no cell over the 2,048 limit;
- density flat in the field away from the bodies;
- lone-body leftover force small next to the signal;
- the answer changes by less than ~20% between rungs.

**Cost guess** (rough; scaled from q59_a's 4.6 h on the RTX 5060 and §15.57's 5090 speed):

- **h 0.25:** ~1.5 h per run, all 6 runs in one 4-job batch, about $2.
- **h 0.177:** ~2.5 h per run, about $3 per batch.
- **h 0.125:** ~3.5 h for the lone proton alone, about $2.

Ladder total: **~$10–15**, plus about half a day of wall clock. I'll measure the first rung before quoting the rest.

## Not in this plan

- The axis-restoring mechanism (separate design).
- More than 4 photon levels (the level-3 cap). It matters once the field is resolved; noted for later.
- Changing defaults: everything above stays behind the `CFM_P_*` knobs until you sign off.
