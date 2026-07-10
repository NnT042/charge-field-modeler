//! Through-charge flow network — Phase A (session 32).
//!
//! Models charge CONDUCTED through bodies and re-emitted elsewhere, per
//! docs/THROUGH_CHARGE_DESIGN.md (paper basis cited there): protons
//! channel pole→equator ("2D fans"), neutrons pole→pole at 67% ("1D
//! lightning rods", haf.pdf), in-alpha posts divert disc output back to
//! the axial channel (atmo2.pdf), and stacked alphas are tied together
//! by the vertical through-streams (graphene.pdf).
//!
//! Phase A is STATE + DIAGNOSTICS only: `solve_charge_flow` runs a
//! relaxation sweep during stepping and fills each particle's
//! [`FlowState`]; no force term reads it yet. Phase B (stream tension +
//! alignment stiffness + throughput-scaled emission) builds on this
//! state and re-earns the full stability battery before shipping.

use glam::DVec3;

use crate::atom_core::{AtomCore, CHANNEL_TAIL, NUCLEON_PITCH};

/// Free-field recycling baseline of one proton (natural units of
/// charge/sec — everything else is expressed relative to this).
pub const FLOW_BASELINE_PROTON: f64 = 1.0;
/// Neutrons channel at 67% of a proton (haf.pdf; the magnetic-moment
/// ratio 1.913/2.793).
pub const FLOW_NEUTRON_FACTOR: f64 = 0.67;
/// A fed proton can channel above its free-field baseline — deut.pdf
/// puts tritium's through-charge at 1.37 and graphene.pdf gives an
/// alpha ~2 proton-units of intake — but not without limit; beyond this
/// the excess registers as `FlowState::stress` (graphene.pdf: an
/// over-channeling nucleus is stressed and must reconfigure).
pub const FLOW_PROTON_CAPACITY: f64 = 2.0;
/// Cap on the fraction of a proton's output that can leave through ONE
/// pole as through-charge when fully demanded downstream; the disc fan
/// keeps the rest (pole-to-equator is "the main and default route",
/// fourier.pdf). Two fully-plugged poles pass 2×0.35 = 0.7, leaving
/// ≥30% lateral — the proton never stops being a fan.
pub const FLOW_PROTON_POLE_PASS: f64 = 0.35;
/// Run the relaxation sweep every N physics steps. Flow re-routing is
/// not instantaneous physically, and at dt = 5e-4 a stale-by-8-steps
/// network is 4 ms of sim time — far below any mechanical timescale in
/// the sandbox.
pub const FLOW_SOLVE_EVERY: usize = 8;

/// Live charge-flow state of one particle. Pole index convention:
/// 0 = south (−pole_axis), 1 = north (+pole_axis).
#[derive(Debug, Clone, Copy)]
pub struct FlowState {
    /// Charge/sec entering at each pole (ambient + captured streams).
    pub intake: [f64; 2],
    /// Through-charge leaving each pole (pole index = EXIT side).
    pub out_pole: [f64; 2],
    /// Pole→equator output emitted at the disc (0 for neutrons).
    pub lateral: f64,
    /// Intake beyond capacity (graphene.pdf overload — future
    /// reconfiguration/instability signal). 0 in a healthy network.
    pub stress: f64,
    /// Total output relative to the free-field baseline — the Phase B
    /// force multiplier. 1.0 for a free particle by construction.
    pub mult: f64,
}

impl Default for FlowState {
    fn default() -> Self {
        Self {
            intake: [0.0; 2],
            out_pole: [0.0; 2],
            lateral: 0.0,
            stress: 0.0,
            mult: 1.0,
        }
    }
}

impl FlowState {
    pub fn through(&self) -> f64 {
        self.out_pole[0] + self.out_pole[1]
    }
    pub fn output(&self) -> f64 {
        self.through() + self.lateral
    }
}

/// Distance falloff of port-to-port capture: 1 inside the funnel mouth
/// (`NUCLEON_PITCH`), smoothstep to 0 at `CHANNEL_TAIL·NUCLEON_PITCH` —
/// the same geometry as the channeling factor (the funnel mouth IS the
/// capture cross-section).
fn capture_falloff(r: f64) -> f64 {
    let lo = NUCLEON_PITCH;
    let hi = CHANNEL_TAIL * NUCLEON_PITCH;
    let t = ((r - lo) / (hi - lo)).clamp(0.0, 1.0);
    1.0 - t * t * (3.0 - 2.0 * t)
}

/// One directed capture edge found by the pass-1 geometry scan:
/// emitter particle `i`'s port (`0`/`1` = pole exits, `2` = disc) →
/// receiver particle `j`'s intake pole `p`, with raw geometric coupling
/// `g ∈ [0,1]`.
struct CaptureEdge {
    i: usize,
    port: usize,
    j: usize,
    p: usize,
    g: f64,
}

impl AtomCore {
    /// (baseline, capacity, is_pole_to_pole_conduit) for a profile.
    /// Unknown profiles fall back to mass-scaled proton behavior (mass
    /// is 1.0 for nucleons and 1/1836 for the electron, so the electron
    /// comes out as the negligible conduit it should be).
    fn flow_params(&self, profile_id: usize) -> (f64, f64, bool) {
        let prof = &self.profiles[profile_id];
        match prof.name.as_str() {
            "proton" => (FLOW_BASELINE_PROTON, FLOW_PROTON_CAPACITY, false),
            "neutron" => (
                FLOW_BASELINE_PROTON * FLOW_NEUTRON_FACTOR,
                FLOW_BASELINE_PROTON * FLOW_NEUTRON_FACTOR,
                true,
            ),
            _ => (prof.mass, prof.mass, false),
        }
    }

    /// One Gauss–Seidel-style relaxation sweep of the flow network,
    /// reading last sweep's outputs and writing fresh
    /// intake/routing/mult into every particle's `FlowState`. Called
    /// from `step` every [`FLOW_SOLVE_EVERY`] steps; converges to the
    /// steady network over a few sweeps and tracks geometry changes
    /// smoothly after that.
    pub fn solve_charge_flow(&mut self) {
        let n = self.particles.len();
        if n == 0 {
            return;
        }
        let max_r = CHANNEL_TAIL * NUCLEON_PITCH;

        // Pass 1: geometric capture edges, from LAST sweep's outputs.
        let mut edges: Vec<CaptureEdge> = Vec::new();
        // Total geometric demand per emitter port (for output-side
        // normalization: a port's output cannot be captured more than
        // once) and per receiver pole (for ambient replacement: a
        // plugged pole receives the stream INSTEAD of ambient).
        let mut port_demand = vec![[0.0f64; 3]; n];
        let mut pole_capture = vec![[0.0f64; 2]; n];

        for j in 0..n {
            let pos_j = self.particles[j].position;
            let pole_j = self.particles[j].pole_axis();
            for i in 0..n {
                if i == j {
                    continue;
                }
                let d_vec = pos_j - self.particles[i].position;
                let r = d_vec.length();
                if r >= max_r || r < 1e-9 {
                    continue;
                }
                let d_hat = d_vec / r;
                let fall = capture_falloff(r);
                let pole_i = self.particles[i].pole_axis();

                for p in 0..2usize {
                    // Receiving pole faces outward along ±pole_j; it
                    // captures flow arriving INTO that face.
                    let face = if p == 1 { pole_j } else { -pole_j };
                    let recv = (-d_hat).dot(face).max(0.0).powi(2);
                    if recv < 1e-6 {
                        continue;
                    }
                    // Emitter pole ports: through-charge leaves along
                    // ±pole_i and only reaches j if j lies that way.
                    for port in 0..2usize {
                        let out_dir = if port == 1 { pole_i } else { -pole_i };
                        let emit = out_dir.dot(d_hat).max(0.0).powi(2);
                        let g = fall * emit * recv;
                        if g > 1e-6 {
                            edges.push(CaptureEdge { i, port, j, p, g });
                            port_demand[i][port] += g;
                            pole_capture[j][p] += g;
                        }
                    }
                    // Emitter disc port: lateral output travels
                    // radially in the disc plane (both sides) — sin²
                    // of the emitter pole vs the line.
                    let emit_disc = 1.0 - pole_i.dot(d_hat).powi(2);
                    let g = fall * emit_disc * recv;
                    if g > 1e-6 {
                        edges.push(CaptureEdge { i, port: 2, j, p, g });
                        port_demand[i][2] += g;
                        pole_capture[j][p] += g;
                    }
                }
            }
        }

        // Pass 2: intakes = ambient (reduced where a pole is plugged)
        // + normalized captured streams from last sweep's outputs.
        let prev: Vec<FlowState> = self.particles.iter().map(|p| p.flow).collect();
        let mut intake = vec![[0.0f64; 2]; n];
        for j in 0..n {
            let (baseline, _, _) = self.flow_params(self.particles[j].profile_id);
            for p in 0..2usize {
                let plugged = pole_capture[j][p].min(1.0);
                intake[j][p] = 0.5 * baseline * (1.0 - plugged);
            }
        }
        for e in &edges {
            let out = match e.port {
                2 => prev[e.i].lateral,
                q => prev[e.i].out_pole[q],
            };
            if out <= 0.0 {
                continue;
            }
            // Output-side normalization: the port's captures share its
            // actual output, never exceeding it.
            let share = e.g / port_demand[e.i][e.port].max(1.0);
            intake[e.j][e.p] += out * share;
        }

        // Pass 3: conservation, routing, multipliers.
        for j in 0..n {
            let (baseline, capacity, conduit) =
                self.flow_params(self.particles[j].profile_id);
            let total_in = intake[j][0] + intake[j][1];
            let output = total_in.min(capacity);
            let stress = (total_in - capacity).max(0.0);

            let flow = &mut self.particles[j].flow;
            flow.intake = intake[j];
            flow.stress = stress;
            flow.mult = if baseline > 1e-12 { output / baseline } else { 0.0 };

            if conduit {
                // Neutron: 1D pipe — what enters one pole exits the
                // other; no disc output (neutron.pdf).
                let split = if total_in > 1e-12 {
                    intake[j][0] / total_in
                } else {
                    0.5
                };
                flow.out_pole = [output * (1.0 - split), output * split];
                flow.lateral = 0.0;
            } else {
                // Proton (fan): each pole passes through-charge only to
                // the extent something downstream captures it, capped
                // at FLOW_PROTON_POLE_PASS; the disc keeps the rest.
                let mut through = [0.0f64; 2];
                for q in 0..2usize {
                    through[q] =
                        output * FLOW_PROTON_POLE_PASS * port_demand[j][q].min(1.0);
                }
                flow.out_pole = through;
                flow.lateral = (output - through[0] - through[1]).max(0.0);
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::atom_scenarios::standard_core;
    use crate::atom_core::NucleusDynamics;

    /// Settle the flow network (several sweeps happen over this many
    /// steps at FLOW_SOLVE_EVERY cadence).
    fn settle(core: &mut AtomCore, steps: usize) {
        core.running = true;
        core.step_n(steps);
    }

    #[test]
    fn free_proton_flow_is_baseline() {
        let mut core = standard_core();
        let pid = core.profile_id_by_name("proton").unwrap();
        core.spawn_particle(pid, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .unwrap();
        settle(&mut core, 100);
        let f = core.particles[0].flow;
        assert!(
            (f.mult - 1.0).abs() < 1e-9,
            "a free proton recycles exactly its baseline: mult={}",
            f.mult
        );
        assert!(
            f.through() < 1e-12,
            "nothing downstream — all output is lateral: through={}",
            f.through()
        );
        assert!(f.stress < 1e-12, "no overload on a free particle");
    }

    #[test]
    fn flow_network_conserves_and_feeds_the_stack() {
        let mut core = standard_core();
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        settle(&mut core, 400);

        // Conservation: every particle's output ≤ its intake, and
        // nothing is NaN/negative.
        let mut stack_through = 0.0f64;
        let mut post_intake = 0.0f64;
        for p in &core.particles {
            let f = p.flow;
            assert!(f.output() <= f.intake[0] + f.intake[1] + 1e-9);
            assert!(f.intake[0] >= 0.0 && f.intake[1] >= 0.0);
            assert!(f.mult.is_finite() && f.stress.is_finite());
            let name = core.profiles[p.profile_id].name.as_str();
            if name == "proton" {
                stack_through += f.through();
            } else if name == "neutron" {
                post_intake += f.intake[0] + f.intake[1];
            }
        }
        // The stacked protons carry a real axial through-stream
        // (graphene.pdf: the vertical streams tie the alphas together)…
        assert!(
            stack_through > 0.1,
            "carbon stack should carry through-charge: {stack_through}"
        );
        // …and the posts receive flow (atmo2.pdf diversion) despite
        // sitting broadside to the discs.
        assert!(
            post_intake > 0.1,
            "posts should divert/capture flow: {post_intake}"
        );
    }

    #[test]
    fn plugged_stack_outchannels_lone_alpha() {
        let mut core = standard_core();
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        settle(&mut core, 400);
        let lone_through: f64 = core.particles.iter().map(|p| p.flow.through()).sum();

        let mut core = standard_core();
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        settle(&mut core, 400);
        let carbon_through: f64 =
            core.particles.iter().map(|p| p.flow.through()).sum();

        assert!(
            carbon_through > lone_through * 1.5,
            "a plugged stack channels more than a lone alpha \
             (haf.pdf: bare stacks channel weakly): lone={lone_through} \
             carbon={carbon_through}"
        );
    }

    /// Session-32 Phase A network dump: the flow table for the preset
    /// ladder. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_flow_network --nocapture`
    #[test]
    #[ignore]
    fn report_flow_network() {
        for preset in ["alpha", "carbon", "neon", "argon"] {
            let mut core = standard_core();
            core.spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("preset");
            core.set_nucleus_dynamics(NucleusDynamics::RigidLock);
            settle(&mut core, 800);

            println!(
                "\n== {preset} ==\n{:>3} {:>8} {:>6} {:>8} {:>8} {:>8} {:>8} {:>8} {:>7} {:>7}",
                "id", "kind", "alpha", "in_S", "in_N", "out_S", "out_N", "lateral", "mult", "stress"
            );
            let mut tot_in = 0.0f64;
            let mut tot_out = 0.0f64;
            for (id, p) in core.particles.iter().enumerate() {
                let f = p.flow;
                tot_in += f.intake[0] + f.intake[1];
                tot_out += f.output();
                println!(
                    "{id:>3} {:>8} {:>6} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>7.3} {:>7.4}",
                    core.profiles[p.profile_id].name,
                    p.alpha.map(|a| a.to_string()).unwrap_or_default(),
                    f.intake[0],
                    f.intake[1],
                    f.out_pole[0],
                    f.out_pole[1],
                    f.lateral,
                    f.mult,
                    f.stress,
                );
            }
            println!("   totals: intake={tot_in:.4} output={tot_out:.4}");
        }
    }
}
