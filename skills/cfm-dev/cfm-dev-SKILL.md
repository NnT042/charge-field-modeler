---
name: cfm-dev
description: >
  Development skill for the Charge Field Modeler (Godot 4 + Rust GDExtension + Vulkan).
  Use this skill for ALL coding tasks in the charge-field-modeler project: Rust spin stack
  math, GDScript UI wiring, Vulkan compute shaders, quaternion operations, Godot scene
  setup, and physics simulation logic. Also use when the user mentions spin stacking,
  photon tiers, charge field simulation, focus particles, or field photons. This skill
  includes a local LLM delegation pattern — use it to offload boilerplate generation to
  a local Qwen model via LM Studio to conserve Claude Code tokens.
---

# Charge Field Modeler — Development Skill

## Project Context

Read these before starting any session:
- `docs/PROJECT_DESIGN.md` — architecture, data structures, unit system, milestones
- `docs/PHYSICS_REFERENCE.md` — charge field theory rules driving the simulation

The stack: Godot 4.x + Rust via gdext + Vulkan 1.1 compute shaders.
Target hardware floor: AMD RX 550 4GB.

---

## Local LLM Delegation

A Qwen 2.5 Coder 14B model runs in LM Studio on the dev machine.
Use it to save Claude Code tokens on mechanical tasks.

### How to call it

```bash
./scripts/delegate_to_local.sh "Your prompt here"
```

Or for longer prompts:
```bash
cat <<'EOF' | ./scripts/delegate_to_local.sh -
Here is a Rust struct:

pub struct SpinLevel {
    pub level: u8,
    pub spin_type: SpinType,
    // ...
}

Write an impl block that provides:
- new() constructor with defaults
- a method to compute amplitude as 2^(level-1)
- Display trait implementation
EOF
```

The script hits `localhost:1234` (LM Studio's OpenAI-compatible endpoint).
If it fails, remind the user to start the server via LM Studio's system tray icon.

### What to delegate (DO delegate these)

| Task | Why it's safe to offload |
|------|--------------------------|
| Rust impl blocks from a struct definition | Narrow scope, clear spec |
| GDScript UI boilerplate (signals, sliders, panels) | Repetitive, well-defined |
| Unit test skeletons when function signatures exist | Mechanical |
| GLSL shader stubs from a compute pipeline spec | Template-driven |
| Match arms, enum variants, From/Into impls | Pure boilerplate |
| Doc comments and inline documentation | Low-stakes, easy to review |
| Repetitive file scaffolding (scene files, resource files) | Copy-paste territory |
| Simple math utilities (lerp, clamp, remap) | One-function scope |

### What to keep in Claude Code (do NOT delegate these)

| Task | Why it needs Claude |
|------|---------------------|
| Quaternion composition logic for spin stacking | Core physics, must match theory |
| Godot ↔ Rust GDExtension interop wiring | Cross-boundary, fiddly, error-prone |
| Vulkan RenderingDevice pipeline setup | Complex API, easy to get wrong |
| Debugging build failures or runtime crashes | Needs full project context |
| Architecture decisions (what goes where) | Needs PROJECT_DESIGN.md context |
| Collision/impulse math | Physics correctness matters |
| Performance optimization passes | Needs profiling context and tradeoff reasoning |
| Anything touching the unit system or SI conversion | Must match Mathis's constants exactly |

### Delegation workflow

1. **Identify** a mechanical task during your work
2. **Write a tight prompt** — include the struct/signature/spec the local model needs
3. **Call the script** and capture output
4. **Review the output** — fix anything wrong before integrating
5. **Integrate** into the codebase

Always review. The local model is fast but not infallible. Catch type errors,
missing imports, and wrong gdext API usage before committing.

### If the server is down

If `delegate_to_local.sh` returns a connection error:
- Tell the user: "LM Studio server isn't responding. Right-click the LM Studio icon in
  your system tray and click 'Start Server', then I'll retry."
- Do NOT silently fall back to doing the work yourself — the whole point is token savings.
  Wait for confirmation, then retry.

---

## Coding Standards

### Rust (gdext)

- All public types derive `GodotClass` where they need to be exposed to Godot
- Use `#[godot_api]` for methods callable from GDScript
- Quaternion math: use `glam` crate (already a gdext dependency), NOT nalgebra
- Error handling: `Result<T, E>` internally, log errors at the Godot boundary
- All spin-related math in natural units (r=1, c=1). SI conversion only in `units.rs`

### GDScript

- Static typing everywhere: `var speed: float = 0.0`
- Signals for all UI → simulation communication
- No physics logic in GDScript — it calls into Rust

### GLSL Compute Shaders

- Target Vulkan 1.1 / SPIR-V 1.3
- Use `layout(local_size_x = 64)` as default workgroup size (tuned for GCN)
- Field photon data as SSBOs, not UBOs (too large for uniform buffers)
- Compile with: `glslc -fshader-stage=compute -o output.spv input.glsl`

---

## Milestone Reference

Current milestones (check PROJECT_DESIGN.md for latest status):

- **M1**: Spinning sphere — Godot scaffold, Rust GDExtension, PBR sphere, orbit camera
- **M2**: Spin stacking levels 1-4 — quaternion composition, path trace, ghost sphere
- **M3**: Full spin stack levels 5-12 — tier transitions, GPU-accelerated path compute
- **M4**: Field simulation — GPU compute pipeline, collision detection, impulse tracking
- **M5**: Multi-particle interactions — inter-particle fields, alpha/H2 formation tests
- **M6**: Polish and verification — compare against theory predictions

---

## Session Startup Checklist

1. Read `docs/PROJECT_DESIGN.md` for current architecture
2. Check git log for recent changes: `git log --oneline -10`
3. Verify Rust builds: `cd rust && cargo build`
4. Verify Godot project opens: check `godot/project.godot` exists
5. Test local LLM connection: `./scripts/delegate_to_local.sh "Say hello"`
