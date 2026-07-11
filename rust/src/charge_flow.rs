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

use crate::atom_core::{AtomCore, NucleusDynamics, CHANNEL_TAIL, NUCLEON_PITCH};

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

/// Phase B default: stream-tension coupling. Each capture link carries
/// `flow` charge/sec; interrupting it costs pressure, so the pair is
/// pulled together with force `FLOW_TENSION × flow` along the link
/// (graphene.pdf: "the vertical charge streams ... tie them together
/// pretty tightly"). No 1/r² — a captured stream is a conduit, not a
/// radiant field; the capture falloff already bounds its range.
/// Runtime field `AtomCore::flow_tension` (seeded from this) so the
/// calibration sweep and the debug panel can vary it.
pub const DEFAULT_FLOW_TENSION: f64 = 1.0;
/// Phase B default: channel-alignment stiffness. The cos² capture
/// gates are differentiated into restoring torques on both ends of
/// every link — the configuration-space stiffness the session-32 whirl
/// analysis proved velocity damping cannot provide, and the deut.pdf
/// post-regulator mechanism. Runtime field `AtomCore::flow_align`.
pub const DEFAULT_FLOW_ALIGN: f64 = 0.5;

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
    /// Total geometric capture demand per output port ([south pole,
    /// north pole, disc]) from the last sweep — the share normalizer
    /// Phase B needs to evaluate a single link's flow at force time
    /// without re-scanning the whole network.
    pub demand: [f64; 3],
}

impl Default for FlowState {
    fn default() -> Self {
        Self {
            intake: [0.0; 2],
            out_pole: [0.0; 2],
            lateral: 0.0,
            stress: 0.0,
            mult: 1.0,
            demand: [0.0; 3],
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
            flow.demand = port_demand[j];

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

/// Phase C1 default: ambient-confinement pressure. 0.0 = OFF — the
/// constant ships only after `report_confinement_calibration` finds a
/// value that binds the FreeNucleon alpha AND re-earns the full
/// RigidAlpha battery (the surface push also compresses stacks).
pub const DEFAULT_AMBIENT_CONFINE: f64 = 0.0;

impl AtomCore {
    /// Phase C1 (session 32, docs/THROUGH_CHARGE_DESIGN.md): ambient
    /// confinement — the external charge field presses on the SURFACE of
    /// a composite (haf.pdf: a stack "will be getting hit by charge from
    /// the sides"; strong.html: "there is no charge field within the
    /// nucleus", so interior members feel nothing). For every member of
    /// every group (in the force-simulated modes), push INWARD along the
    /// outward direction from the group's live center, scaled by the
    /// member's EXPOSURE — how unblocked its outward view is by
    /// farther-out members. A buried member is shielded; a member
    /// drifting out of the cluster becomes exposed and gets pressed back
    /// in. External-field force: no internal reaction pair (the momentum
    /// comes from the ambient field), and a symmetric composite feels no
    /// net force.
    ///
    /// This is deliberately NOT pairwise shadow glue (`nuclear_ambient`):
    /// pairwise glue scales with every interior pair and over-compresses
    /// the stack long before it can hold a surface; surface pressure
    /// confines without touching interior balance.
    pub(crate) fn apply_ambient_confinement(&mut self) {
        let k = self.ambient_confine;
        if k.abs() < 1e-12 || self.dynamics == NucleusDynamics::RigidLock {
            return;
        }
        for gi in 0..self.groups.len() {
            let members: Vec<usize> = self.groups[gi].members.clone();
            if members.len() < 2 {
                continue;
            }
            let (com, _) = self.live_group_frame(gi);
            for &i in &members {
                let out_vec = self.particles[i].position - com;
                let s = out_vec.length();
                if s < 1e-9 {
                    continue; // dead center — fully buried
                }
                let u = out_vec / s;
                // Blocking: members beyond i along u shadow its outward
                // view; solid-angle proxy (r_j/d)² times cos² alignment
                // of the blocker with the outward direction.
                let mut blocked = 0.0f64;
                for &j in &members {
                    if j == i {
                        continue;
                    }
                    let d_vec = self.particles[j].position - self.particles[i].position;
                    let along = d_vec.dot(u);
                    if along <= 0.0 {
                        continue;
                    }
                    let d2 = d_vec.length_squared().max(1e-9);
                    let rj = self.profiles[self.particles[j].profile_id].radius;
                    let cos2 = (along * along) / d2;
                    blocked += (rj * rj / d2).min(1.0) * cos2;
                }
                let exposure = 1.0 / (1.0 + blocked);
                let mi = self.profiles[self.particles[i].profile_id].mass;
                self.particles[i].force_accum -= u * (k * exposure * mi);
            }
        }
    }

    /// Phase B: stream tension + channel-alignment stiffness for one
    /// pair, evaluated at CURRENT geometry with flow amplitudes (and
    /// port-demand normalizers) from the last solver sweep. Returns
    /// `(force_on_j, torque_on_i, torque_on_j)`; the force on i is the
    /// exact opposite (the stream pulls both ends together — no
    /// momentum injection, unlike the vortex term's known defect).
    ///
    /// Both directions (i feeds j, j feeds i) are summed. Per capture
    /// link the energy is `−k·fall(r)·emit(θ_e)·recv(θ_v)·flow`, so:
    /// - tension: `k·flow_link` along the line (falloff bounds range);
    /// - torques: analytic derivatives of the cos²/sin² gates —
    ///   `τ = −(∂E/∂c)(u×w)` for each orientation-dependent gate,
    ///   restoring the pole/disc alignment that keeps the channel open
    ///   (the deut.pdf regulator; the configuration-space stiffness the
    ///   whirl saga proved missing).
    pub(crate) fn flow_tension_pair(
        &self,
        i: usize,
        j: usize,
        d_hat: DVec3,
        r: f64,
    ) -> (DVec3, DVec3, DVec3) {
        let k_t = self.flow_tension;
        let k_a = self.flow_align;
        if k_t.abs() < 1e-12 && k_a.abs() < 1e-12 {
            return (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        }
        let fall = capture_falloff(r);
        if fall < 1e-6 {
            return (DVec3::ZERO, DVec3::ZERO, DVec3::ZERO);
        }

        let mut f_on_j = DVec3::ZERO;
        let mut tau = [DVec3::ZERO, DVec3::ZERO]; // [on i, on j]

        // (emitter, receiver, d̂ from emitter to receiver). Tension on j
        // is ALWAYS −d_hat (toward i), whichever end emits — a stream
        // pulls both of its ends together. The original Phase B code
        // signed the j-emits branch +d_hat, turning every link where j
        // is the emitter into REPULSION; that inverted force was the
        // plug-retention killer (report_plug_retention, session 32: the
        // socket's out-pole stream blew the plug neutron off instead of
        // gluing it — electron.pdf's "charge wind" acts as a glue).
        for &(e, v, dir) in &[(i, j, d_hat), (j, i, -d_hat)] {
            let pole_e = self.particles[e].pole_axis();
            let pole_v = self.particles[v].pole_axis();
            let fe = &self.particles[e].flow;
            let tau_e_idx = if e == i { 0 } else { 1 };
            let tau_v_idx = 1 - tau_e_idx;

            for p in 0..2usize {
                let face = if p == 1 { pole_v } else { -pole_v };
                let c_recv = (-dir).dot(face).max(0.0);
                if c_recv < 1e-6 {
                    continue;
                }
                let recv = c_recv * c_recv;

                // Pole (through-charge) ports of the emitter.
                for q in 0..2usize {
                    let out = fe.out_pole[q];
                    if out <= 1e-12 {
                        continue;
                    }
                    let out_dir = if q == 1 { pole_e } else { -pole_e };
                    let c_emit = out_dir.dot(dir).max(0.0);
                    if c_emit < 1e-6 {
                        continue;
                    }
                    let emit = c_emit * c_emit;
                    let share = 1.0 / fe.demand[q].max(1.0);
                    let flow = fall * emit * recv * out * share;
                    f_on_j -= d_hat * (k_t * flow);
                    // Alignment: receiver face toward the arrival line,
                    // emitter exit toward the departure line.
                    let k_link = k_a * fall * out * share;
                    tau[tau_v_idx] +=
                        (face.cross(-dir)) * (2.0 * k_link * emit * c_recv);
                    tau[tau_e_idx] +=
                        (out_dir.cross(dir)) * (2.0 * k_link * recv * c_emit);
                }

                // Disc (lateral) port of the emitter: sin² gate — the
                // restoring torque keeps the disc PLANE containing the
                // line to the receiver.
                let out = fe.lateral;
                if out > 1e-12 {
                    let c_disc = pole_e.dot(dir);
                    let emit = 1.0 - c_disc * c_disc;
                    if emit > 1e-6 {
                        let share = 1.0 / fe.demand[2].max(1.0);
                        let flow = fall * emit * recv * out * share;
                        f_on_j -= d_hat * (k_t * flow);
                        let k_link = k_a * fall * out * share;
                        tau[tau_v_idx] +=
                            (face.cross(-dir)) * (2.0 * k_link * emit * c_recv);
                        tau[tau_e_idx] -=
                            (pole_e.cross(dir)) * (2.0 * k_link * recv * c_disc);
                    }
                }
            }
        }
        (f_on_j, tau[0], tau[1])
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

    /// All pairwise distances between the given particle ids.
    fn pair_dists(core: &AtomCore, members: &[usize]) -> Vec<f64> {
        let mut d = Vec::new();
        for a in 0..members.len() {
            for b in (a + 1)..members.len() {
                d.push(core.pair_distance(members[a], members[b]));
            }
        }
        d
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

    /// Session-32 Phase B calibration: sweep flow_tension × flow_align
    /// against BOTH regimes that must hold simultaneously — the
    /// FreeNucleon lone alpha (intra-alpha cohesion, currently
    /// dissolves at ~2500% drift) and the RigidAlpha carbon quiet run
    /// (inter-alpha binding, must not over-compress or destabilize).
    /// For the FreeNucleon rows the worst post-proton pair's force
    /// breakdown is printed at the end so the ejection mechanism is
    /// visible, not just the drift. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_flow_calibration --nocapture`
    #[test]
    #[ignore]
    fn report_flow_calibration() {
        println!(
            "\n{:>8} {:>7} | {:>12} {:>10} | {:>12} {:>10}",
            "tension", "align", "free_drift", "free_KE", "rigid_drift", "rigid_KE"
        );
        for &tension in &[0.0, 1.0, 2.0, 4.0, 8.0, 16.0] {
            for &align in &[0.0, 0.5 * tension] {
                if tension == 0.0 && align != 0.0 {
                    continue;
                }
                // FreeNucleon lone alpha, 20k steps.
                let mut core = standard_core();
                core.flow_tension = tension;
                core.flow_align = align;
                let gid = core
                    .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("alpha preset");
                core.set_nucleus_dynamics(NucleusDynamics::FreeNucleon);
                let members = core.groups[gid].members.clone();
                let d0 = pair_dists(&core, &members);
                core.running = true;
                core.step_n(20_000);
                let d1 = pair_dists(&core, &members);
                let free_drift = d0
                    .iter()
                    .zip(&d1)
                    .map(|(a, b)| (a - b).abs() / a.max(1e-9))
                    .fold(0.0f64, f64::max);
                let free_ke = core.total_kinetic_energy();

                // RigidAlpha 3-stack quiet, 20k steps (the binding
                // harness structure — see `tri_alpha`'s doc).
                let mut core = standard_core();
                core.flow_tension = tension;
                core.flow_align = align;
                let gid = core
                    .spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("tri_alpha preset");
                core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
                core.running = true;
                let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                let n_alphas = core.groups[gid].alphas.len();
                let pairs: Vec<(usize, usize)> = (0..n_alphas)
                    .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                    .collect();
                let d0: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                core.step_n(20_000);
                let rigid_drift = pairs
                    .iter()
                    .enumerate()
                    .map(|(k, &(a, b))| {
                        let d = (com(&core, a) - com(&core, b)).length();
                        ((d - d0[k]) / d0[k]).abs()
                    })
                    .fold(0.0f64, f64::max);
                let rigid_ke = core.total_kinetic_energy();

                println!(
                    "{tension:>8.1} {align:>7.2} | {:>11.1}% {free_ke:>10.4} | {:>11.1}% {rigid_ke:>10.4}",
                    free_drift * 100.0,
                    rigid_drift * 100.0,
                );
            }
        }

        // Ejection anatomy at tension=0: what actually blows the free
        // alpha apart? Print every intra-alpha pair's breakdown at the
        // rest pose after a short settle.
        let mut core = standard_core();
        core.flow_tension = 0.0;
        core.flow_align = 0.0;
        let gid = core
            .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        core.set_nucleus_dynamics(NucleusDynamics::FreeNucleon);
        core.running = true;
        core.step_n(100);
        let members = core.groups[gid].members.clone();
        println!(
            "\n-- free-alpha pair anatomy at t≈0 (tension off) --\n\
             {:>5} {:>6} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>8} {:>9}",
            "pair", "r", "channel", "grav", "charge", "ambient", "intake", "stream", "contact", "tension"
        );
        for a in 0..members.len() {
            for b in (a + 1)..members.len() {
                let (i, j) = (members[a], members[b]);
                let fb = core.pair_force_breakdown(i, j);
                println!(
                    "{a}-{b:<3} {:>6.3} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} {:>+9.4}",
                    fb[0], fb[1], fb[2], fb[3], fb[4], fb[5], fb[6], fb[7], fb[8],
                );
            }
        }
    }

    /// Session-32 Phase B: watch WHERE the FreeNucleon alpha goes at a
    /// given tension instead of only the endpoint drift — per-member
    /// positions (cylindrical: radial offset from the alpha axis, y) and
    /// all pair distances, sampled through the run. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_free_alpha_trajectory --nocapture`
    #[test]
    #[ignore]
    fn report_free_alpha_trajectory() {
        for &tension in &[0.0, 8.0, 16.0] {
            let mut core = standard_core();
            core.flow_tension = tension;
            core.flow_align = tension * 0.5;
            let gid = core
                .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                .expect("alpha preset");
            core.set_nucleus_dynamics(NucleusDynamics::FreeNucleon);
            core.running = true;
            let members = core.groups[gid].members.clone();

            println!(
                "\n== tension={tension} ==\n{:>6} {:>10} {:>44} {:>36}",
                "step", "KE", "members (lat,y)", "pair r (01 02 03 12 13 23)"
            );
            for s in 0..10 {
                core.step_n(4_000);
                let com: DVec3 = members
                    .iter()
                    .map(|&m| core.particles[m].position)
                    .sum::<DVec3>()
                    / members.len() as f64;
                let cyl: Vec<String> = members
                    .iter()
                    .map(|&m| {
                        let d = core.particles[m].position - com;
                        format!("({:.2},{:+.2})", (d.x * d.x + d.z * d.z).sqrt(), d.y)
                    })
                    .collect();
                let dists: Vec<String> = pair_dists(&core, &members)
                    .iter()
                    .map(|d| format!("{d:.2}"))
                    .collect();
                println!(
                    "{:>6} {:>10.4} {:>44} {:>36}",
                    (s + 1) * 4_000,
                    core.total_kinetic_energy(),
                    cyl.join(" "),
                    dists.join(" "),
                );
            }
        }
    }

    /// Session-32 Phase C1 calibration: ambient-confinement pressure
    /// swept against the three regimes it must serve simultaneously —
    /// FreeNucleon lone alpha (cohesion: the arch needs a surface push),
    /// RigidAlpha tri_alpha (the earned binding tables must not be
    /// over-compressed), and RigidAlpha plugged carbon (plug retention,
    /// report_carbon_stability's 218-409% baseline). Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_confinement_calibration --nocapture`
    #[test]
    #[ignore]
    fn report_confinement_calibration() {
        println!(
            "\n{:>8} | {:>12} {:>10} | {:>12} | {:>12}",
            "confine", "free_drift", "free_KE", "tri_drift", "carbon_drift"
        );
        for &confine in &[0.0, 0.25, 0.5, 1.0, 2.0, 4.0] {
            // FreeNucleon lone alpha, 30k steps.
            let free = {
                let mut core = standard_core();
                core.ambient_confine = confine;
                let gid = core
                    .spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("alpha preset");
                core.set_nucleus_dynamics(NucleusDynamics::FreeNucleon);
                let members = core.groups[gid].members.clone();
                let d0 = pair_dists(&core, &members);
                core.running = true;
                core.step_n(30_000);
                let d1 = pair_dists(&core, &members);
                let drift = d0
                    .iter()
                    .zip(&d1)
                    .map(|(a, b)| (a - b).abs() / a.max(1e-9))
                    .fold(0.0f64, f64::max);
                (drift, core.total_kinetic_energy())
            };

            // RigidAlpha drift over 30k steps for a preset's alpha units.
            let rigid_drift = |preset: &str, confine: f64| -> f64 {
                let mut core = standard_core();
                core.ambient_confine = confine;
                let gid = core
                    .spawn_preset(preset, DVec3::ZERO, DVec3::ZERO, DVec3::Y)
                    .expect("preset");
                core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
                core.running = true;
                let com = |core: &AtomCore, ai: usize| core.groups[gid].alphas[ai].com;
                let n_alphas = core.groups[gid].alphas.len();
                let pairs: Vec<(usize, usize)> = (0..n_alphas)
                    .flat_map(|a| ((a + 1)..n_alphas).map(move |b| (a, b)))
                    .collect();
                let d0: Vec<f64> = pairs
                    .iter()
                    .map(|&(a, b)| (com(&core, a) - com(&core, b)).length())
                    .collect();
                let mut worst = 0.0f64;
                for _ in 0..15 {
                    core.step_n(2_000);
                    for (idx, &(a, b)) in pairs.iter().enumerate() {
                        let d = (com(&core, a) - com(&core, b)).length();
                        worst = worst.max(((d - d0[idx]) / d0[idx]).abs());
                    }
                }
                worst
            };
            let tri = rigid_drift("tri_alpha", confine);
            let carbon = rigid_drift("carbon", confine);

            println!(
                "{confine:>8.2} | {:>11.1}% {:>10.4} | {:>11.1}% | {:>11.1}%",
                free.0 * 100.0,
                free.1,
                tri * 100.0,
                carbon * 100.0,
            );
        }
    }

    /// Session-32 Phase A network dump: the flow table for the preset
    /// ladder. Run:
    /// `cargo test --release --manifest-path rust/Cargo.toml -- --ignored
    ///  report_flow_network --nocapture`
    #[test]
    #[ignore]
    fn report_flow_network() {
        for preset in ["alpha", "carbon", "tri_alpha", "neon", "argon"] {
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
