# M5 Test Loop — Headless Physics Iteration

How to iterate on atom-mode physics **without opening the Godot editor**.
All simulation logic lives in pure Rust (`rust/src/atom_core.rs`); the Godot
node `AtomSim` is a thin wrapper, so the headless harness exercises exactly
the physics the app runs.

## The loop

1. **Physics gates** (~10 s, run after every change):

   ```
   cargo test --release --manifest-path rust/Cargo.toml
   ```

   Scenario tests live in `rust/src/atom_scenarios.rs`:
   - `hydrogen_capture` — baseline canary: electron falls into the polar intake
   - `hydrogen_orbit_stable_long_run` — wall-riding orbit at r=1.3, θ≈11°,
     v_tan = corotation speed, 500k steps
   - `derived_constants_equilibrium` — locked constants satisfy the radial
     equilibrium equation (wall-riding inside the tunnel, standoff outside)
   - `h2_bond_matrix` — 8 spin/pole combos: 4 bond, 4 repel, no coded rule
   - `alpha_holds_and_conserves`, `alpha_captures_electron`, `helium_stable`,
     `heavier_presets_smoke` — nuclear presets

2. **Trajectory dumps** (when a test fails and you need to see *why*):

   ```
   cargo run --release --manifest-path rust/Cargo.toml --bin atom_lab -- <scenario> [steps] [sample_every] [out.csv]
   ```

   Scenarios: `hydrogen`, `protons`, `h2` (full 8-combo verdict table),
   `h2:N` (single combo 0–7), `alpha`, `helium`.
   Any `key=value` argument overrides a coupling for sweeps without
   recompiling: `g_q c_q intake vortex drag torque corot stream p_amb dt`.
   Without a CSV path it prints summary `ScenarioMetrics`; with one it dumps
   per-frame positions/velocities plus r, v_tan, v_rad, θ.

3. **Integration gate** (before committing UI/wiring changes):

   ```
   cargo build --release --manifest-path rust/Cargo.toml
   & "D:\App\Godot\Godot_v4.6.1-stable_win64.exe" --headless --path godot res://scenes/atom_mode.tscn --quit-after 300
   & "D:\App\Godot\Godot_v4.6.1-stable_win64.exe" --headless --path godot res://scenes/main.tscn --quit-after 300
   ```

   Pass = exit 0 and no `SCRIPT ERROR` lines. (`field_sim` prints one
   intended warning headless — the GPU field sim needs a real device.)

   The scene smoke never leaves the default "protons" scenario, so also run
   the scripted smoke, which drives preset/bond scenario switches
   (rigid groups → composite skins, h2 → bond bridges):

   ```
   & "D:\App\Godot\Godot_v4.6.1-stable_win64.exe" --headless --path godot -s res://scripts/dev_smoke_atom.gd
   ```

   Pass = exit 0 and `[smoke] OK`.

4. Human visual check only at milestone boundaries.

## Rules

- The locked coupling constants live in `Couplings::default()`
  (rust/src/atom_core.rs) with their Mathis-sourced derivations. **Change
  them only with a failing scenario test as justification** — the tuning
  panel (HUD → Debug toggle) is for experiments, not for shipping values.
- Always build **release**: the .gdextension points at `target/release/`.
- Never `cd` into rust/ — use `--manifest-path rust/Cargo.toml`.
