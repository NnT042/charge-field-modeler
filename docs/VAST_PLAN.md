# Rented-GPU plan (Vast.ai) — draft 2026-11-18

Nothing is rented yet. This is the list of what a big GPU buys us and what each run would cost, so
we can pick the first experiment before setting up an account instance.

## What limits us today

Current runs (tc48): mean free path λ 0.3, cell h 0.25, step dt 0.06, domain radius R_D 6, 8000
ambient photons per unit volume, about 7M photons, 125 per cell. One frozen pair run takes about 10
minutes on the 3060 (12 GB) or the 5060 (8 GB). DSMC is limited by memory bandwidth, not
arithmetic: every step reads and writes every photon (64 bytes) several times (sort, collide,
move).

Rules from §15.14 (converged field): h about λ, dt at most λ/5. In report_deuteron λ sets the
cross-section σ only; the photon count (n_a 8000) does not change. So halving λ costs:
- cells: (h_old/h_new)^3 = 8x for the same domain (cheap: the cell arrays are small)
- photons: the same
- steps: 2x, plus a longer warm-up (about R_D²/λ)
- total: roughly 2-3x the work, and about the same memory.

Correction (2026-11-19): an earlier draft assumed the photons per cell stay fixed (16x the work).
They don't. The λ 0.15 runs (§15.50) ran locally on the 5060 in a few hours each. Rented GPUs
pay off only where the photon count itself must grow (finer h at FIXED photons per cell, the
furnace, many seeds).

## Candidate rented GPUs

| GPU | memory | bandwidth (TB/s) | vs 3060 | rough $/h (median) |
|---|---|---|---|---|
| RTX 3060 (ours) | 12 GB | 0.36 | 1x | — |
| RTX 4090 | 24 GB | 1.0 | ~2.8x | ~0.35-0.5 |
| RTX 5090 | 32 GB | 1.8 | ~5x | ~0.71-0.79 |
| H100 80GB | 80 GB | 3.3 | ~9x | ~1.73-1.93 |

Bandwidth ratio is the speed-up to expect, not the shader-core count.

## Experiments ranked by what they would teach us

### A. Resolution check of today's results — DONE locally (§15.50)
λ 0.15/0.3/0.6/1.0 and h 0.125 at λ 0.3 ran on our own GPUs. Halving h at λ 0.3 raises the push
11-28% (above the seed spread); local pairing recovers most of that. The sign and trends hold at
every λ. Not needed on a rented GPU.

### B. Fine grid near the bodies (the "high-resolution grid")
The beam is 0.222 wide and the clear hole 0.222; the cell is 0.25. At h 0.05 (λ 0.05) the beam spans
~4 cells, so a scattered photon cannot pair with a beam photon 0.2 away (the random-pairing
artefact local pairing was built to dodge). Uniform h 0.05 over R_D 6 is 125x the cells: ~900M
photons, ~60 GB of photons alone. Too big even for an H100. Two ways out:
- a smaller domain around the pair (R_D 3: 8x fewer cells, but the free run hit R_D 6 and the
  boundary feed must stay isotropic), or
- a nested grid: fine cells only within ~1.5 of each body. That's an engine change (§15.13 noted it
  as the lever before). It's the right long-term answer, and it also makes the 3060 useful at fine
  resolution.

### C. Many seeds in parallel (cheapest per answer)
A 4090/5090 runs several 7M-photon jobs side by side. Good for anything where the effect is near
the seed spread (local pairing, beam width, tilt torques).

### D. "Stellar furnace" free-for-all (for fun, and for what emerges)
User idea (2026-11-18): a dense, hot field with full spin-up/spin-down mechanics, many free
nucleons, run long, and watch whether any of them lock their axial streams together (Mathis:
fused nucleons are bound by one shared through-stream along the main axis).
- Bodies free to move AND rotate (rotation from §15.20, gyroscopic L_s), trace bodies, through
  release, emergent intake, no freezing, no layout.
- Dense field: λ 0.1 (the density at which the H atom and free H2 held, §15.25/15.32).
- Engine limits to lift first:
  - The emitter tag is 3 bits (photon flag bits 17-19), so at most 7 bodies are tagged.
  - The arrival histograms take at most 4 bodies (HIST_MAX_BODIES).
  - A "furnace" of, say, 8 protons + 8 neutrons needs 16+ tags, or untagged bodies with tallies
    off.
- Scale (corrected 2026-11-19): at the same n_a, R_D 8 is (8/6)^3 ≈ 2.4x the photons (~17M, fits
  our 12 GB card). The cost is time: dt ≤ λ/5 = 0.02 and a long free run. For "stellar" density,
  raise n_a itself (more photons per volume, not just a smaller λ). That is the knob that needs a
  big card: 10x n_a ≈ 170M photons ≈ 11 GB of photons plus buffers → 5090/H100.
- Readouts: pair distances over time, which pairs stay within 2.5 for 100+ tu, axis alignment of
  close pairs (pole-to-pole vs equator-to-pole), whether a shared axial stream forms (flux along the
  line of centres).

## Engine work before renting (all useful on our own GPUs too)
1. Near-exit counter results (§15.49) decide whether random pairing is inventing beam collisions.
   Done: mostly not (it adds ~14%); the medium is thick (§15.49-15.50).
2. Lift the 7-body tag limit (widen the tag field or drop the tag on extra bodies).
3. Headless "multi-GPU / multi-seed" runner on Linux: the exe already runs self-contained; check
   Vulkan on the rented box with `vulkaninfo` (NVIDIA_DRIVER_CAPABILITIES must include graphics).
4. Optional: the nested near-body grid (B), which is the bigger win.

## Status 2026-11-19 (night): the photon cap comes first

- **The engine holds at most ~19M photons per run on any card.** The photon buffer is one
  storage binding, and our cards report a 2 GiB binding limit (`report_deuteron` clamps the
  capacity to 19M for that reason). A rented 24-80 GB card does not raise this by itself.
  Check `max_storage_buffer_binding_size` on the rented box (the gpu tests print it); some
  NVIDIA drivers report 4 GiB, which would give ~2x.
- So experiment B (high resolution) needs engine work first: split the photon buffer across
  several bindings (raises the cap with memory), or the nested near-body grid (lowers the
  photons needed). Both are useful on our own cards too.
- What a rented card can do with no engine work: C (many seeds or settings side by side; a
  24 GB card fits several of our ~7M-photon jobs at once), and runs up to 19M photons faster.
- The overnight local queue (reference/tc56) includes a hi-res pilot at the cap: n_a 16000
  (~14.5M photons) and h 0.125, same mean free path. If the lane / edge results move a lot
  between h 0.25 and h 0.125, the buffer split is worth doing before renting.

## Morning checklist (user does the account and payment steps)
1. Vast.ai account, add a small credit (user only; Claude does not enter payment details).
2. Add an SSH public key to the account (ESTELLE's `~/.ssh/id_ed25519.pub`, or a new key).
3. Pick an offer: RTX 4090 or 5090, reliability ≥ 0.98, an Ubuntu 22.04/24.04 CUDA template.
   Set `NVIDIA_DRIVER_CAPABILITIES=all` in the template env (otherwise there is no Vulkan
   inside the container).
4. On the box: `apt install libvulkan1 vulkan-tools`, `vulkaninfo --summary` (must list the
   NVIDIA GPU, not llvmpipe), install rustup, clone the repo, build with
   `cargo build --manifest-path rust/Cargo.toml --profile fast --features gpu --tests`.
5. Run the gates once (`gpu::dsmc`, 23 pass), note the binding limit, then a queue script
   modelled on reference/tc56/q56e.sh; scp results back.
6. Destroy the instance when done (storage bills while stopped).

## Setup steps (to walk through with the user later)
Pick a 4090/5090 offer with a CUDA/Ubuntu template → add an SSH key → confirm Vulkan
(`vulkaninfo --summary`) → build the Linux exe (or build on the box) → copy a queue script like
reference/tc48/q48.sh → scp results back. Destroy the instance when done (storage still bills
while stopped).
