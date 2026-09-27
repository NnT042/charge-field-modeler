# Charge Field Modeler — Claude Code Context

## Read First

Before any coding session, read:
- `docs/PROJECT_DESIGN.md` — architecture, milestones, current status
- `docs/PHYSICS_REFERENCE.md` — charge field theory rules driving the simulation

For development tasks (Rust, GDScript, shaders, Godot wiring), load the `cfm-dev` skill.

---

## Looking Up Mathis's Physics

**Use the `search_mathis` MCP tool. Do not browse milesmathis.com.**

The Mathis Library is a home-lab-hosted MCP server attached as an account-level
connector (not a project `.mcp.json` entry). Its tools may be deferred — load them with
ToolSearch (`search_mathis`) before calling. It holds complete mirrors of both of
Mathis's sites and returns passages with source URLs for verification.

| Tool | Use |
|---|---|
| `search_mathis` | Science papers (milesmathis.com) — the default for physics lookups |
| `search_mathis_all` | Both sites at once |
| `search_mathis_art` | Art/history/genealogy site (mileswmathis.com) — rarely relevant here |
| `get_page` | Read a whole paper by URL once a search hit looks relevant |
| `list_pages` | Browse papers by title if you're not sure what to search for |

```
search_mathis("stacked spins angular momentum")
search_mathis("charge field photon recycling proton")
search_mathis("axial spin radius electron")
```

Never dispatch a web search or sub-agent to milesmathis.com — the mirror covers the
full site and is available instantly.

---

## Key Rules

- **Spin math stays in Rust.** No physics logic in GDScript — it calls into Rust.
- **Units:** All internal math in natural units (r=1, c=1). SI conversion only in `units.rs`.
- **Quaternions:** Use `glam`, not `nalgebra`.
- **Boilerplate:** Offload to the local Qwen model via `./scripts/delegate_to_local.sh`.
  Do not use Claude Code tokens for mechanical tasks the local model can handle.
- **Constants must match Mathis.** Any numerical constant tied to charge field theory
  needs a source — use `search_mathis` to verify before committing.

---

## Stack

Godot 4 + Rust (gdext) + Vulkan 1.1 compute shaders. Target hardware: AMD RX 550 4GB.
