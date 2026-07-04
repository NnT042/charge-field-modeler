//! Pure-Rust physics core for atom mode — no Godot dependencies.
//!
//! Everything simulation-related lives here so it can run headless:
//! the Godot node `AtomSim` (atom_sim.rs) is a thin wrapper, and the
//! scenario harness (atom_scenarios.rs) + `atom_lab` bin drive this
//! directly for fast physics iteration without opening the editor.

use glam::{DMat3, DQuat, DVec3};
use std::f64::consts::{FRAC_PI_2, TAU};
use std::path::{Path, PathBuf};

// ── Simple PRNG (xorshift64) ─────────────────────────────────────────────

pub struct Rng {
    state: u64,
}

impl Rng {
    pub fn new(seed: u64) -> Self {
        Self {
            state: seed.max(1),
        }
    }

    pub fn next_u64(&mut self) -> u64 {
        let mut x = self.state;
        x ^= x << 13;
        x ^= x >> 7;
        x ^= x << 17;
        self.state = x;
        x
    }

    pub fn next_f64(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / ((1u64 << 53) as f64)
    }
}

// ── Emission Table ────────────────────────────────────────────────────────
// 91 bins: index 0 = pole (θ=0°), index 90 = equator (θ=90°).
// Bilateral symmetry assumed — both poles identical.

#[derive(Clone)]
pub struct EmissionTable {
    pub bins: Vec<f32>,
}

impl EmissionTable {
    pub fn sample(&self, cos_theta: f64) -> f64 {
        self.sample_theta(cos_theta.abs().acos()) // bilateral fold: 0..π/2
    }

    fn sample_theta(&self, theta: f64) -> f64 {
        let t = (theta / FRAC_PI_2).clamp(0.0, 1.0);
        let idx = ((t * (self.bins.len() - 1) as f64) as usize).min(self.bins.len() - 1);
        self.bins[idx] as f64
    }

    /// dE/dθ (per radian) at the folded polar angle, central difference.
    pub fn d_dtheta(&self, cos_theta: f64) -> f64 {
        let theta = cos_theta.abs().acos();
        let h = 0.02;
        let t1 = (theta - h).max(0.0);
        let t2 = (theta + h).min(FRAC_PI_2);
        if t2 - t1 < 1e-9 {
            return 0.0;
        }
        (self.sample_theta(t2) - self.sample_theta(t1)) / (t2 - t1)
    }

    pub fn from_csv_values(csv: &[f32]) -> Self {
        assert!(csv.len() >= 180, "Need at least 180 CSV values (-90° to +89°)");
        let n = csv.len();
        let mut bins = Vec::with_capacity(91);
        for theta in 0..=90u32 {
            let north_idx = n.saturating_sub(1).min((180u32.saturating_sub(theta)) as usize);
            let south_idx = (theta as usize).min(n - 1);
            let north = csv[north_idx];
            let south = csv[south_idx];
            bins.push((north + south) / 2.0);
        }
        let max = bins.iter().copied().fold(0.0f32, f32::max);
        if max > 0.0 {
            for v in &mut bins {
                *v /= max;
            }
        }
        Self { bins }
    }

    pub fn complement(&self) -> Self {
        Self {
            bins: self.bins.iter().map(|v| 1.0 - v).collect(),
        }
    }
}

// ── Emission CDF (for weighted sampling) ─────────────────────────────────

#[derive(Clone)]
pub struct EmissionCdf {
    cdf: Vec<f64>,
}

impl EmissionCdf {
    pub fn from_bins(bins: &[f32]) -> Self {
        let n = bins.len();
        let mut cdf = Vec::with_capacity(n);
        let mut sum = 0.0f64;
        for (i, &val) in bins.iter().enumerate() {
            let theta = FRAC_PI_2 * i as f64 / (n - 1).max(1) as f64;
            sum += val as f64 * theta.sin().max(0.001);
            cdf.push(sum);
        }
        if sum > 0.0 {
            for v in &mut cdf {
                *v /= sum;
            }
        }
        Self { cdf }
    }

    pub fn sample(&self, u: f64) -> f64 {
        let u = u.clamp(0.0, 0.9999);
        let idx = self.cdf.partition_point(|&v| v < u).min(self.cdf.len() - 1);
        FRAC_PI_2 * idx as f64 / (self.cdf.len() - 1).max(1) as f64
    }
}

// ── Particle Profile ──────────────────────────────────────────────────────

#[derive(Clone)]
pub struct ParticleProfile {
    pub name: String,
    pub mass: f64,
    pub radius: f64,
    pub emission: EmissionTable,
    pub absorption: EmissionTable,
    pub emission_cdf: EmissionCdf,
    /// A²-weighted CDF — samples the polar intake cone for the vortex cloud.
    pub intake_cdf: EmissionCdf,
}

// ── Sim Particle ──────────────────────────────────────────────────────────

pub struct SimParticle {
    pub profile_id: usize,
    pub position: DVec3,
    pub velocity: DVec3,
    pub orientation: DQuat,
    pub angular_velocity: DVec3,
    pub force_accum: DVec3,
    pub torque_accum: DVec3,
    /// Index into AtomCore::groups when this particle is a rigid-composite
    /// constituent (nuclear preset); None for free particles.
    pub group: Option<usize>,
}

impl SimParticle {
    pub fn pole_axis(&self) -> DVec3 {
        self.orientation * DVec3::Y
    }
}

// ── Nuclear presets (rigid composites) ───────────────────────────────────

#[derive(Clone)]
pub struct Constituent {
    pub profile_name: &'static str,
    pub local_pos: DVec3,
    pub local_pole: DVec3,
    pub spin_sign: f64,
    /// Carousel member: rides around the group's stack axis (the alpha's
    /// neutron posts are not immovable — they roll around the midsection,
    /// riding the disc outputs between the two protons, staying 180° apart).
    pub carousel: bool,
}

/// One alpha block centered at `y` on the stack axis: two protons in a
/// short stack (CDs stacked hole-to-hole, spinning the SAME direction so
/// charge channels through pole-to-pole as a dipole —
/// milesmathis.com/oxygen.pdf), with the two neutrons as posts between the
/// disks keeping the protons from turning (oxygen.pdf). Posts sit off-axis
/// so the axial charge channel stays open (nuclear.pdf: the hole in the CD
/// is the recycling channel). Exact post positions are a modeling choice
/// within those constraints.
fn alpha_block(y: f64) -> Vec<Constituent> {
    vec![
        Constituent {
            profile_name: "proton",
            local_pos: DVec3::new(0.0, y - 0.9, 0.0),
            local_pole: DVec3::Y,
            spin_sign: 1.0,
            carousel: false,
        },
        Constituent {
            profile_name: "proton",
            local_pos: DVec3::new(0.0, y + 0.9, 0.0),
            local_pole: DVec3::Y,
            spin_sign: 1.0,
            carousel: false,
        },
        Constituent {
            profile_name: "neutron",
            local_pos: DVec3::new(-0.5, y, 0.0),
            local_pole: DVec3::Y,
            spin_sign: 1.0,
            carousel: true,
        },
        Constituent {
            profile_name: "neutron",
            local_pos: DVec3::new(0.5, y, 0.0),
            local_pole: DVec3::Y,
            spin_sign: 1.0,
            carousel: true,
        },
    ]
}

/// Polar plug proton: plugged into a stack pole with its pole PERPENDICULAR
/// to the stack axis, so its equatorial disc output feeds the stack's open
/// polar channel (phos.pdf plug-and-socket; ammon.pdf: N = C-stack + proton
/// in the south pole, neutron in the north; O = protons in both poles).
/// `z` offsets the plug off-axis so a proton+neutron pair can share the
/// polar hole side by side (0 for a lone plug); a paired proton points its
/// pole AT its partner across the hole. Plugs ride the carousel like the
/// alpha posts — they roll around the polar socket rather than sitting
/// welded (session-29 user note).
fn plug_proton(y: f64, z: f64) -> Constituent {
    let local_pole = if z.abs() > 1e-9 {
        DVec3::new(0.0, 0.0, -z.signum()) // toward the paired plug
    } else {
        DVec3::X
    };
    Constituent {
        profile_name: "proton",
        local_pos: DVec3::new(0.0, y, z),
        local_pole,
        spin_sign: 1.0,
        carousel: true,
    }
}

/// Polar plug neutron: pole ON the stack axis — "the neutron is plugged in
/// with its pole pointing down … protons channel charge pole to equator,
/// while neutrons channel pole to pole" (graphene.pdf), so its axial pole
/// pulls charge into the stack's hole (atmo2.pdf: paired neutrons are
/// "pulling charge into the axial holes"). `z` offsets it off-axis to sit
/// side by side with a paired plug proton (0 for a lone plug). Rides the
/// carousel like the alpha posts.
fn plug_neutron(y: f64, z: f64) -> Constituent {
    Constituent {
        profile_name: "neutron",
        local_pos: DVec3::new(0.0, y, z),
        local_pole: DVec3::Y,
        spin_sign: 1.0,
        carousel: true,
    }
}

/// Half-gap between the members of a proton+neutron pair sharing a polar
/// hole ("two baryons in the hole fill the hole much better" — atmo2.pdf).
/// Same nestling scale as the alpha's neutron posts (±0.5 off-axis).
const PLUG_PAIR_GAP: f64 = 0.55;

/// Alpha stack pitch: adjacent alpha centers along the axis. Tight enough
/// that the disks read as plugged (nuclear.pdf), loose enough to see the
/// blocks.
const ALPHA_PITCH: f64 = 2.6;

pub fn preset_constituents(name: &str) -> Option<Vec<Constituent>> {
    match name {
        "alpha" => Some(alpha_block(0.0)),
        // Carbon: three alphas stacked (nuclear.pdf: "Carbon blocks — three
        // alphas stacked"; the single-stack limit that makes C the basis of
        // life).
        "carbon" => Some(
            [
                alpha_block(-ALPHA_PITCH),
                alpha_block(0.0),
                alpha_block(ALPHA_PITCH),
            ]
            .concat(),
        ),
        // Nitrogen: carbon stack + 7th proton plugged in the south pole
        // (edge-on, disc feeding the hole) and the balancing neutron in the
        // north (ammon.pdf), pole-down per graphene.pdf.
        "nitrogen" => {
            let mut c = preset_constituents("carbon")?;
            c.push(plug_proton(-(ALPHA_PITCH + 1.8), 0.0));
            c.push(plug_neutron(ALPHA_PITCH + 1.8, 0.0));
            Some(c)
        }
        // Oxygen: carbon stack + BOTH poles capped by a proton+neutron PAIR
        // (oxygen.pdf: the 7th and 8th protons go on the ends because four
        // alphas can't stack; atmo2.pdf: their neutrons are paired with
        // them in the hole). Each pair sits side by side — proton edge-on
        // (disc feeds the hole), neutron pole-down (channels axially).
        "oxygen" => {
            let mut c = preset_constituents("carbon")?;
            let y = ALPHA_PITCH + 1.8;
            c.push(plug_proton(-y, -PLUG_PAIR_GAP));
            c.push(plug_neutron(-y, PLUG_PAIR_GAP));
            c.push(plug_proton(y, -PLUG_PAIR_GAP));
            c.push(plug_neutron(y, PLUG_PAIR_GAP));
            Some(c)
        }
        _ => None,
    }
}

/// Preset names for UI listings.
pub fn preset_names() -> &'static [&'static str] {
    &["alpha", "carbon", "nitrogen", "oxygen"]
}

/// A rigidly-locked composite (nucleus). Constituents remain real particles
/// (forces sampled per-constituent — an electron can capture at one specific
/// proton's pole), but they integrate as one rigid body. Intra-group pair
/// forces are skipped: the nucleus is pre-fused, its internal balance is
/// not simulated (uf4.pdf: "the alphas can't be broken and rearranged").
pub struct RigidGroup {
    pub members: Vec<usize>,
    pub local_offsets: Vec<DVec3>,
    pub local_orients: Vec<DQuat>,
    /// Intrinsic axial spin rate of each member (rad/s about its own pole).
    pub member_spin: Vec<f64>,
    /// Accumulated axial spin phase per member — members visibly rotate
    /// about their own poles even though the group frame is rigid.
    pub member_spin_phase: Vec<f64>,
    /// Members that ride the carousel around the group's stack axis
    /// (alpha neutron posts, polar plugs). Kinematic: preserves ALL pair
    /// distances — every off-axis member rides the carousel at the same
    /// phase, everything else is on-axis. Riders' orientations rotate with
    /// the ride (they roll around the axis, not slide).
    pub carousel: Vec<bool>,
    pub carousel_rate: f64,
    pub carousel_phase: f64,
    pub com: DVec3,
    pub velocity: DVec3,
    pub orientation: DQuat,
    pub angular_velocity: DVec3,
    pub mass: f64,
    /// Scalar inertia approximation: Σ m(|offset|² + 0.4 r²).
    pub inertia: f64,
    /// Azimuth-averaged skin reach vs polar angle (SKIN_REACH_RINGS samples
    /// over θ ∈ [0, π]), computed once at spawn — the group is rigid.
    /// Feeds the composite skin mesh, the emission-ring overlay, and the
    /// flow-parcel exit paths (parcels die AT this boundary).
    pub skin_reach: Vec<f64>,
    /// (ring radius, y) per alpha block: where each alpha's disc meets the
    /// reach surface. Feeds the ring overlay and the flow disc exits.
    pub disc_exits: Vec<(f64, f64)>,
}

/// Samples in a group's stored `skin_reach` table.
pub const SKIN_REACH_RINGS: usize = 49;

/// Linear interpolation into a reach-vs-θ table (θ ∈ [0, π]).
pub fn reach_at_theta(table: &[f64], theta: f64) -> f64 {
    if table.is_empty() {
        return 0.0;
    }
    if table.len() == 1 {
        return table[0];
    }
    let t = (theta / std::f64::consts::PI).clamp(0.0, 1.0) * (table.len() - 1) as f64;
    let i = (t as usize).min(table.len() - 2);
    let frac = t - i as f64;
    table[i] * (1.0 - frac) + table[i + 1] * frac
}

// ── VFX Particle (charge emission sprinkler) ─────────────────────────────

pub struct VfxParticle {
    position: DVec3,
    velocity: DVec3,
    age: f64,
    lifetime: f64,
    color: (f32, f32, f32),
    base_scale: f32,
    /// Scale growth per second — smoke expands as it drifts.
    grow: f32,
    kind: VfxKind,
}

/// Motion model for a cloud particle.
enum VfxKind {
    /// Free ballistic puff (bond stream bridges).
    Smoke,
    /// Charge parcel on a complete recycling path: polar funnel in,
    /// through the body/stack, out where the emission maps send it (disc
    /// spray, far pole, alpha ring, or captured by a cap electron).
    /// `pts` is an arc-length-uniform polyline in the anchor's local
    /// frame (pole = +Y, fast axial spin excluded); the whole path swirls
    /// about the pole at `swirl` rad/s so intake AND exhaust visibly
    /// corotate. Color lerps from the intake tint to `color_end` as the
    /// parcel recycles through.
    Flow {
        anchor: VfxAnchor,
        swirl: f64,
        color_end: (f32, f32, f32),
        pts: Vec<DVec3>,
    },
}

/// What a flow parcel's path is expressed relative to.
enum VfxAnchor {
    /// A free particle: position + pole frame (build_frame — the visual
    /// swirl is applied separately at a readable rate, not the physical
    /// TAU·3 axial spin).
    Particle(usize),
    /// A rigid group: com + orientation. One frame for the whole nucleus,
    /// so axial paths span the full stack.
    Group(usize),
}

pub const SOFTENING: f64 = 0.05;
pub const ANGULAR_DAMPING: f64 = 0.998;
/// Charge-field lock on nuclei (1/s): the ambient field that fused a
/// nucleus also LOCKS its orientation (nuclear.pdf: "The charge field then
/// locks them into these configurations") — group tumble relaxes as
/// exp(−RELAX·t) instead of letting every passing particle pump angular
/// momentum into the nucleus ("wild turning", session-29 user report).
/// Torques still act, so slow molecular alignment remains possible; the
/// steady state is ω ≈ τ/(I·RELAX).
pub const GROUP_SPIN_RELAX: f64 = 20.0;
pub const MIN_RENDER_RADIUS: f32 = 0.15;
pub const ENVELOPE_MIN_R: f64 = 0.25;
pub const CONTACT_STIFFNESS: f64 = 100.0;

/// Corotation speed ceiling (fraction of c=1). The polar vortex spins with
/// the emitter, but matter riding it stays well sub-c — the render-scaled
/// axial spin rate (TAU·3 ≈ 19 rad/s) would otherwise demand superluminal
/// corotation at the contact wall, which is exactly the energy injection
/// that made the session-27 intake vortex slingshot. Sub-c cap per the
/// spin-model rule that top-level motion stays below c (see
/// docs/PHYSICS_REFERENCE.md; observed orbital speeds ~0.05c).
pub const COROT_V_MAX: f64 = 0.25;

// ── Force couplings ──────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Couplings {
    /// Gravity coupling (1/r², isotropic, attractive).
    pub g_q: f64,
    /// Charge coupling (1/r⁴, emission × absorption, repulsive).
    pub c_q: f64,
    /// Ambient isotropic charge pressure (pushes into charge shadows).
    pub ambient_pressure: f64,
    /// Torque coupling (sin 2θ equator-alignment).
    pub torque: f64,
    /// Equatorial vortex drag (tangential, ∝ E(θ), dead at poles).
    pub vortex: f64,
    /// Doppler velocity correction on charge.
    pub drag: f64,
    /// Polar intake (radial 1/r² sink + channeling toward pole axis, ∝ A²).
    pub intake: f64,
    /// Corotation drag: pulls tangential velocity toward the local speed of
    /// the emitter's polar charge vortex (self-limiting tangential drive).
    pub corot: f64,
    /// Polar stream-collision cushion: repulsion between two facing intake
    /// streams (∝ A_i²·A_j²·(m_i·m_j)²/r⁴). Sets the molecular standoff.
    pub stream: f64,
}

impl Default for Couplings {
    /// LOCKED defaults (M5 lockdown). These are not sliders any more — the
    /// relationships between them are pinned by the scenario tests
    /// (`derived_constants_equilibrium`, `hydrogen_orbit_stable_long_run`,
    /// `h2_bond_matrix`). Change only with a failing test as justification.
    ///
    /// Derivation chain:
    /// - `g_q = 1.0` — definition of the natural force unit (all other
    ///   couplings are expressed relative to gravity at r=1).
    /// - `c_q = 500` — sets the equatorial emission wall. With the loaded
    ///   histograms this puts the free equilibrium radius
    ///   r_eq = √(C_q·E(θ)·R/(G_q+I_q·A²)) INSIDE the contact wall for polar
    ///   angles (⇒ captured electrons rest ON the nuclear boundary,
    ///   milesmathis.com/diatom.pdf & neon.pdf: electrons sit at the nuclear
    ///   boundary, circling the pole) and far OUTSIDE it near the equator
    ///   (⇒ equatorial repulsion; bond standoff scale for M2).
    /// - `intake = 0.5` — polar recycling inflow at half gravity strength at
    ///   A²=1; the capture funnel. Charge recycles in at the poles and out
    ///   at the equator (quantumg.html; the vortex capture picture in
    ///   diatom.pdf where the corrected Bohr radius 9e-9 m is the CAPTURE
    ///   limit of this vortex, not an orbit radius — fine4.pdf).
    /// - `corot = 0.5` — relaxation rate of the corotation drag. Steady state
    ///   is set by COROT_V_MAX, not by this value (self-limiting drag), so
    ///   any O(0.1–1) value lands the same orbit; 0.5 converges within the
    ///   capture transient.
    /// - `drag = 0.3` (doppler) and `vortex = 0.1` (equatorial sprinkler) —
    ///   retained from M4 field-sim calibration; both are dissipative or
    ///   pole-dead and do not move the polar orbit (verified by sweep,
    ///   session 28).
    /// - `torque = 0.1` — equator-toward-charge alignment (gear-mesh torque,
    ///   Compton/inverse-Compton collision mechanics).
    /// - `stream = 24.5` — head-on collision of two facing polar intake
    ///   streams (pole-to-pole charge "meeting head-to-head",
    ///   milesmathis.com/fourier.pdf & jup3.pdf). 1/r⁴ (product of two 1/r²
    ///   stream densities), ∝ (m_i·m_j)² since recycling throughput scales
    ///   with mass. Sets the H₂ bond standoff: d = √(stream·A⁴/(G_q+2·I_q))
    ///   ≈ 3.5 for bare facing poles — molecular, well outside nuclear
    ///   contact (2.0). Also why bare protons refuse to fuse at ambient
    ///   pressure (nuclear.pdf: alphas need stars).
    /// - `ambient_pressure = 0.0` — the pairwise shadow term stays available
    ///   for environment effects, but the H₂ bond emerges from
    ///   gravity+intake attraction vs the stream cushion without it.
    fn default() -> Self {
        Self {
            g_q: 1.0,
            c_q: 500.0,
            ambient_pressure: 0.0,
            torque: 0.1,
            vortex: 0.1,
            drag: 0.3,
            intake: 0.5,
            corot: 0.5,
            stream: 24.5,
        }
    }
}

// ── AtomCore: the whole simulation, engine-free ──────────────────────────

pub struct AtomCore {
    pub profiles: Vec<ParticleProfile>,
    pub particles: Vec<SimParticle>,
    pub groups: Vec<RigidGroup>,
    pub couplings: Couplings,
    pub ambient_gravity: DVec3,
    pub ambient_charge: DVec3,
    pub running: bool,
    pub dt: f64,
    pub time: f64,
    rng: Rng,
    vfx_particles: Vec<VfxParticle>,
    pub vfx_enabled: bool,
    /// Wall-clock accumulator for cloud animation (independent of sim time
    /// so the visual rotation rate doesn't change with substeps).
    vfx_time: f64,
}

impl Default for AtomCore {
    fn default() -> Self {
        Self::new()
    }
}

impl AtomCore {
    pub fn new() -> Self {
        Self {
            profiles: Vec::new(),
            particles: Vec::new(),
            groups: Vec::new(),
            couplings: Couplings::default(),
            ambient_gravity: DVec3::ZERO,
            ambient_charge: DVec3::ZERO,
            running: false,
            dt: 0.0005,
            time: 0.0,
            rng: Rng::new(0xDEAD_BEEF_CAFE),
            vfx_particles: Vec::new(),
            vfx_enabled: true,
            vfx_time: 0.0,
        }
    }

    // ── Profile management ────────────────────────────────────────────

    pub fn register_profile(&mut self, name: &str, mass: f64, radius: f64, csv: &[f32]) -> usize {
        let emission = EmissionTable::from_csv_values(csv);
        let absorption = emission.complement();
        let emission_cdf = EmissionCdf::from_bins(&emission.bins);
        let a2_bins: Vec<f32> = absorption.bins.iter().map(|v| v * v).collect();
        let intake_cdf = EmissionCdf::from_bins(&a2_bins);
        let id = self.profiles.len();
        self.profiles.push(ParticleProfile {
            name: name.to_string(),
            mass,
            radius,
            emission,
            absorption,
            emission_cdf,
            intake_cdf,
        });
        id
    }

    pub fn profile_id_by_name(&self, name: &str) -> Option<usize> {
        self.profiles.iter().position(|p| p.name == name)
    }

    // ── Particle management ───────────────────────────────────────────

    pub fn spawn_particle(
        &mut self,
        profile_id: usize,
        pos: DVec3,
        vel: DVec3,
        pole_dir: DVec3,
    ) -> Option<usize> {
        if profile_id >= self.profiles.len() {
            return None;
        }
        let spin_rate = default_spin_rate(&self.profiles[profile_id].name);
        self.spawn_particle_ex(profile_id, pos, vel, pole_dir, spin_rate)
    }

    /// Spawn with an explicit signed axial spin rate (rad/s about the pole).
    /// Needed by the H₂ 8-combination matrix where spin sign matters.
    pub fn spawn_particle_ex(
        &mut self,
        profile_id: usize,
        pos: DVec3,
        vel: DVec3,
        pole_dir: DVec3,
        spin_rate: f64,
    ) -> Option<usize> {
        if profile_id >= self.profiles.len() {
            return None;
        }
        let orientation = orientation_from_pole(pole_dir);
        let pole = orientation * DVec3::Y;

        let id = self.particles.len();
        self.particles.push(SimParticle {
            profile_id,
            position: pos,
            velocity: vel,
            orientation,
            angular_velocity: pole * spin_rate,
            force_accum: DVec3::ZERO,
            torque_accum: DVec3::ZERO,
            group: None,
        });
        Some(id)
    }

    /// Spawn a rigid nuclear preset at `pos` with its axis along `axis`.
    /// Returns the group id, or None for an unknown preset name.
    pub fn spawn_preset(
        &mut self,
        name: &str,
        pos: DVec3,
        vel: DVec3,
        axis: DVec3,
    ) -> Option<usize> {
        let constituents = preset_constituents(name)?;
        let orientation = orientation_from_pole(axis);
        let gid = self.groups.len();

        let mut members = Vec::new();
        let mut local_offsets = Vec::new();
        let mut local_orients = Vec::new();
        let mut member_spin = Vec::new();
        let mut carousel = Vec::new();
        let mut mass = 0.0;
        let mut inertia = 0.0;

        for c in &constituents {
            let profile_id = self.profile_id_by_name(c.profile_name)?;
            let prof_mass = self.profiles[profile_id].mass;
            let prof_radius = self.profiles[profile_id].radius;
            let spin = default_spin_rate(c.profile_name) * c.spin_sign;
            let world_pos = pos + orientation * c.local_pos;
            let world_pole = orientation * c.local_pole;
            let id = self.spawn_particle_ex(profile_id, world_pos, vel, world_pole, spin)?;
            self.particles[id].group = Some(gid);
            members.push(id);
            local_offsets.push(c.local_pos);
            local_orients.push(orientation_from_pole(c.local_pole));
            member_spin.push(spin);
            carousel.push(c.carousel);
            mass += prof_mass;
            inertia += prof_mass * (c.local_pos.length_squared() + 0.4 * prof_radius * prof_radius);
        }

        let n_members = members.len();
        self.groups.push(RigidGroup {
            members,
            local_offsets,
            local_orients,
            member_spin,
            member_spin_phase: vec![0.0; n_members],
            carousel,
            // Posts ride the disc outputs — they roll in time with the
            // proton spin (oxygen.pdf: the neutrons keep the protons from
            // turning by rolling with them, cohering the two fields).
            carousel_rate: default_spin_rate("proton"),
            carousel_phase: 0.0,
            com: pos,
            velocity: vel,
            orientation,
            angular_velocity: DVec3::ZERO,
            mass,
            inertia: inertia.max(1e-9),
            skin_reach: Vec::new(),
            disc_exits: Vec::new(),
        });
        self.sync_group_members(gid);
        let reach = self.compute_group_reach(gid, SKIN_REACH_RINGS);
        self.groups[gid].skin_reach = reach;
        let exits = self.compute_disc_exits(gid);
        self.groups[gid].disc_exits = exits;
        Some(gid)
    }

    /// Reposition a group's members from the group frame and give them the
    /// rigid-body velocity field (v_com + ω×r) plus their intrinsic axial
    /// spin — velocity-dependent forces (doppler, corotation) on
    /// constituents need correct member velocities.
    fn sync_group_members(&mut self, gid: usize) {
        struct MemberSync {
            id: usize,
            offset: DVec3,
            orient: DQuat,
            spin: f64,
            spin_phase: f64,
            carousel: bool,
        }
        let (data, com, vel, orientation, omega, car_rot, car_omega) = {
            let g = &self.groups[gid];
            let data: Vec<MemberSync> = (0..g.members.len())
                .map(|k| MemberSync {
                    id: g.members[k],
                    offset: g.local_offsets[k],
                    orient: g.local_orients[k],
                    spin: g.member_spin[k],
                    spin_phase: g.member_spin_phase[k],
                    carousel: g.carousel[k],
                })
                .collect();
            let car_rot = DQuat::from_rotation_y(g.carousel_phase);
            // Carousel angular velocity in world space (about the stack axis)
            let car_omega = (g.orientation * DVec3::Y) * g.carousel_rate;
            (data, g.com, g.velocity, g.orientation, g.angular_velocity, car_rot, car_omega)
        };
        for ms in data {
            let local_off = if ms.carousel { car_rot * ms.offset } else { ms.offset };
            let world_off = orientation * local_off;
            let p = &mut self.particles[ms.id];
            p.position = com + world_off;
            // Rigid frame × (carousel ride) × constituent frame ×
            // accumulated axial spin — members visibly rotate about their
            // own poles, and carousel riders ROLL around the stack axis
            // (their orientation tracks the ride, so an edge-on plug keeps
            // facing its partner all the way around).
            let member_orient = if ms.carousel { car_rot * ms.orient } else { ms.orient };
            p.orientation =
                (orientation * member_orient * DQuat::from_rotation_y(ms.spin_phase)).normalize();
            p.velocity = vel + omega.cross(world_off);
            if ms.carousel {
                p.velocity += car_omega.cross(world_off);
            }
            let pole = p.pole_axis();
            p.angular_velocity = omega + pole * ms.spin;
            if ms.carousel {
                p.angular_velocity += car_omega;
            }
        }
    }

    pub fn remove_particle(&mut self, id: usize) {
        if id >= self.particles.len() {
            return;
        }
        // Constituents of a rigid preset can't be removed individually.
        if self.particles[id].group.is_some() {
            return;
        }
        let last = self.particles.len() - 1;
        self.particles.swap_remove(id);
        // The particle formerly at `last` now lives at `id` — fix group refs.
        if id != last {
            if let Some(g) = self.particles[id].group {
                for m in &mut self.groups[g].members {
                    if *m == last {
                        *m = id;
                    }
                }
            }
        }
    }

    pub fn clear_particles(&mut self) {
        self.particles.clear();
        self.groups.clear();
        self.vfx_particles.clear();
    }

    pub fn set_particle_pole(&mut self, id: usize, target: DVec3) {
        if let Some(p) = self.particles.get_mut(id) {
            if target.length_squared() > 1e-9 {
                p.orientation = orientation_from_pole(target);
            }
        }
    }

    // ── Calibration & diagnostics ─────────────────────────────────────

    /// Indices of (heaviest, lightest) particles by mass, or None if fewer
    /// than two distinct particles. Used to pick the orbit center vs orbiter.
    pub fn center_and_orbiter(&self) -> Option<(usize, usize)> {
        if self.particles.len() < 2 {
            return None;
        }
        let mut center = 0usize;
        let mut orbiter = 0usize;
        let mut m_max = f64::MIN;
        let mut m_min = f64::MAX;
        for (i, p) in self.particles.iter().enumerate() {
            let m = self.profiles[p.profile_id].mass;
            if m > m_max {
                m_max = m;
                center = i;
            }
            if m < m_min {
                m_min = m;
                orbiter = i;
            }
        }
        if center == orbiter {
            return None;
        }
        Some((center, orbiter))
    }

    /// Auto-calibrate gravity coupling for a circular polar orbit.
    /// Sets G_q so gravity alone supplies the centripetal acceleration:
    /// a_c = v_tan²/r = G_q·m_center/r² ⇒ G_q = v_tan²·r / m_center.
    /// Returns the new G_q (0.0 if it couldn't calibrate).
    pub fn auto_calibrate_polar(&mut self) -> f64 {
        let Some((center, orbiter)) = self.center_and_orbiter() else {
            return 0.0;
        };
        let m_center = self.profiles[self.particles[center].profile_id].mass;
        let d = self.particles[orbiter].position - self.particles[center].position;
        let r = d.length();
        if r < 1e-9 || m_center < 1e-12 {
            return 0.0;
        }
        let d_hat = d / r;
        let v_rel = self.particles[orbiter].velocity - self.particles[center].velocity;
        let v_tan = (v_rel - d_hat * v_rel.dot(d_hat)).length();
        let g = v_tan * v_tan * r / m_center;
        self.couplings.g_q = g;
        g
    }

    /// Orbit diagnostics: [r, v_tan, v_rad, suggested_G_q, emission_on_axis, theta_deg].
    /// theta_deg = angle from center's nearest pole to orbiter, in degrees.
    pub fn orbit_info(&self) -> [f64; 6] {
        let Some((center, orbiter)) = self.center_and_orbiter() else {
            return [0.0; 6];
        };
        let cp = &self.particles[center];
        let op = &self.particles[orbiter];
        let m_center = self.profiles[cp.profile_id].mass;
        let d = op.position - cp.position;
        let r = d.length().max(1e-9);
        let d_hat = d / r;
        let v_rel = op.velocity - cp.velocity;
        let v_rad = v_rel.dot(d_hat);
        let v_tan = (v_rel - d_hat * v_rad).length();
        let suggested_g = v_tan * v_tan * r / m_center.max(1e-12);
        let cos_theta = cp.pole_axis().dot(d_hat);
        let emission = self.profiles[cp.profile_id].emission.sample(cos_theta);
        let theta_deg = cos_theta.abs().acos().to_degrees();
        [r, v_tan, v_rad, suggested_g, emission, theta_deg]
    }

    /// Force breakdown on the orbiter from the center particle.
    /// [F_gravity, F_charge, F_intake_radial, F_channel, F_vortex_eq,
    ///  A_center, A²_center, speed].
    /// Signed scalars: positive = repulsive (away), negative = attractive.
    /// Channel is always ≥ 0 (magnitude toward the axis).
    pub fn force_breakdown(&self) -> [f64; 8] {
        let Some((ci, oi)) = self.center_and_orbiter() else {
            return [0.0; 8];
        };
        let cp = &self.particles[ci];
        let op = &self.particles[oi];
        let c_prof = &self.profiles[cp.profile_id];
        let o_prof = &self.profiles[op.profile_id];

        let d_vec = op.position - cp.position;
        let r = d_vec.length().max(SOFTENING);
        let r2 = r * r;
        let r4 = r2 * r2;
        let d_hat = d_vec / r;

        let pole = cp.pole_axis();
        let cos_theta = pole.dot(d_hat);
        let emission_c = c_prof.emission.sample(cos_theta);
        let absorption_c = c_prof.absorption.sample(cos_theta);
        let pole_o = op.pole_axis();
        let absorption_o = o_prof.absorption.sample(pole_o.dot(-d_hat));
        let ai2 = absorption_c * absorption_c;

        let cq = &self.couplings;
        let f_grav = -(cq.g_q * c_prof.mass * o_prof.mass / r2);
        let v_rel = op.velocity - cp.velocity;
        let v_rad = v_rel.dot(d_hat);
        let doppler = (1.0 - cq.drag * v_rad).clamp(0.2, 5.0);
        let f_charge = cq.c_q * c_prof.mass * o_prof.mass * emission_c * absorption_o / r4 * doppler;
        let f_intake = -(cq.intake * c_prof.mass * o_prof.mass * ai2 / r2);

        let along = d_vec.dot(pole);
        let lateral = d_vec - pole * along;
        let lat_dist = lateral.length();
        let f_channel = if lat_dist > 1e-9 {
            cq.intake * c_prof.mass * o_prof.mass * ai2 * lat_dist / (r * r2)
        } else {
            0.0
        };

        let spin = cp.angular_velocity.dot(pole);
        let sin_vec = pole.cross(d_hat);
        let sin_theta = sin_vec.length();
        let f_vortex = if sin_theta > 1e-9 {
            cq.vortex * c_prof.mass * emission_c * absorption_o * spin * sin_theta / r4
        } else {
            0.0
        };

        let speed = op.velocity.length();

        [
            f_grav, f_charge, f_intake, f_channel, f_vortex,
            absorption_c, ai2, speed,
        ]
    }

    pub fn total_kinetic_energy(&self) -> f64 {
        let mut ke = 0.0;
        for p in &self.particles {
            let m = self.profiles[p.profile_id].mass;
            ke += 0.5 * m * p.velocity.length_squared();
        }
        ke
    }

    pub fn pair_distance(&self, a: usize, b: usize) -> f64 {
        match (self.particles.get(a), self.particles.get(b)) {
            (Some(a), Some(b)) => (a.position - b.position).length(),
            _ => -1.0,
        }
    }

    /// Contact distance (sum of render-clamped radii) for a particle pair.
    pub fn contact_distance(&self, a: usize, b: usize) -> f64 {
        let ra = self
            .particles
            .get(a)
            .map(|p| self.profiles[p.profile_id].radius.max(MIN_RENDER_RADIUS as f64))
            .unwrap_or(0.0);
        let rb = self
            .particles
            .get(b)
            .map(|p| self.profiles[p.profile_id].radius.max(MIN_RENDER_RADIUS as f64))
            .unwrap_or(0.0);
        ra + rb
    }

    // ── Simulation step (velocity Verlet) ─────────────────────────────

    pub fn step(&mut self) {
        if !self.running || self.particles.is_empty() {
            return;
        }
        let dt = self.dt;

        // 1. Forces at current state
        self.compute_forces();

        // 2. Half-step velocities + full-step positions
        self.kick(dt * 0.5);
        self.drift(dt);

        // 3. Forces at new positions
        self.compute_forces();

        // 4. Complete velocity step.
        // No blanket velocity damping: dissipation comes only from physical
        // channels — doppler drag (radial) and corotation drag (tangential).
        // A flat 0.9999/step multiplier was what killed tangential orbit
        // velocity through session 27.
        self.kick(dt * 0.5);
        for i in 0..self.particles.len() {
            if self.particles[i].group.is_some() {
                continue;
            }
            let p = &mut self.particles[i];
            // Preserve axial spin (intrinsic), only damp tumble/precession
            let pole = p.pole_axis();
            let w_axial = pole * p.angular_velocity.dot(pole);
            let w_tumble = p.angular_velocity - w_axial;
            p.angular_velocity = w_axial + w_tumble * ANGULAR_DAMPING;
        }

        self.time += dt;
    }

    /// Velocity (and angular velocity) update for free particles and groups.
    fn kick(&mut self, dt: f64) {
        for i in 0..self.particles.len() {
            if self.particles[i].group.is_some() {
                continue;
            }
            let p = &mut self.particles[i];
            let inv_m = 1.0 / self.profiles[p.profile_id].mass;
            p.velocity += p.force_accum * inv_m * dt;
            p.angular_velocity += p.torque_accum * dt;
        }
        // Groups: aggregate member forces into a COM force + torque.
        for gi in 0..self.groups.len() {
            let (f, tau) = {
                let g = &self.groups[gi];
                let mut f = DVec3::ZERO;
                let mut tau = DVec3::ZERO;
                for &m in &g.members {
                    let p = &self.particles[m];
                    f += p.force_accum;
                    tau += (p.position - g.com).cross(p.force_accum) + p.torque_accum;
                }
                (f, tau)
            };
            let g = &mut self.groups[gi];
            g.velocity += f / g.mass * dt;
            g.angular_velocity += tau / g.inertia * dt;
            // Charge-field lock — see GROUP_SPIN_RELAX.
            g.angular_velocity *= (-GROUP_SPIN_RELAX * dt).exp();
        }
    }

    /// Position (and orientation) update for free particles and groups.
    fn drift(&mut self, dt: f64) {
        for i in 0..self.particles.len() {
            if self.particles[i].group.is_some() {
                continue;
            }
            let p = &mut self.particles[i];
            p.position += p.velocity * dt;
            let w = p.angular_velocity;
            let w_len = w.length();
            if w_len > 1e-12 {
                let rot = DQuat::from_axis_angle(w / w_len, w_len * dt);
                p.orientation = (rot * p.orientation).normalize();
            }
        }
        for gi in 0..self.groups.len() {
            {
                let g = &mut self.groups[gi];
                g.com += g.velocity * dt;
                let w = g.angular_velocity;
                let w_len = w.length();
                if w_len > 1e-12 {
                    let rot = DQuat::from_axis_angle(w / w_len, w_len * dt);
                    g.orientation = (rot * g.orientation).normalize();
                }
                // Kinematic phases: carousel ride + member axial spins.
                g.carousel_phase = (g.carousel_phase + g.carousel_rate * dt) % TAU;
                for k in 0..g.member_spin_phase.len() {
                    g.member_spin_phase[k] =
                        (g.member_spin_phase[k] + g.member_spin[k] * dt) % TAU;
                }
            }
            self.sync_group_members(gi);
        }
    }

    pub fn step_n(&mut self, steps: usize) {
        for _ in 0..steps {
            self.step();
        }
    }

    // ── Force computation ─────────────────────────────────────────────

    /// Per-pair channel occlusion in [0,1]: 1 when a third particle sits in
    /// the polar channel between the pair (within lateral 0.3 of the line,
    /// full fade by 0.8). Models the STOPPERED vortex of diatom.pdf — a
    /// captured electron riding a pole blocks that pole's intake stream.
    fn compute_occlusion(&self) -> Vec<f64> {
        let n = self.particles.len();
        let mut occ = vec![0.0f64; n * n];
        if n < 3 {
            return occ;
        }
        for i in 0..n {
            for j in (i + 1)..n {
                let a = self.particles[i].position;
                let b = self.particles[j].position;
                let ab = b - a;
                let len2 = ab.length_squared();
                if len2 < 1e-12 {
                    continue;
                }
                let mut worst = 0.0f64;
                for k in 0..n {
                    if k == i || k == j {
                        continue;
                    }
                    let t = (self.particles[k].position - a).dot(ab) / len2;
                    if t <= 0.05 || t >= 0.95 {
                        continue;
                    }
                    let closest = a + ab * t;
                    let lat = (self.particles[k].position - closest).length();
                    let block = ((0.8 - lat) / 0.5).clamp(0.0, 1.0);
                    worst = worst.max(block);
                }
                occ[i * n + j] = worst;
                occ[j * n + i] = worst;
            }
        }
        occ
    }

    /// Per-particle, per-pole intake occlusion in [0,1] for the cloud VFX:
    /// `[north, south]` where north = +pole_axis. Same lateral falloff as
    /// the pairwise `compute_occlusion` (full block within 0.3 of the
    /// channel line, faded out by 0.8): anything sitting over a pole within
    /// ~3 units stoppers that pole's visible intake tornado — a wall-riding
    /// electron, a sideways baryon, or the next nucleon in a rigid stack
    /// (interior stack poles read as plugged; only the open ends drain).
    pub fn pole_occlusion(&self) -> Vec<[f64; 2]> {
        self.pole_occlusion_blockers()
            .into_iter()
            .map(|p| [p[0].0, p[1].0])
            .collect()
    }

    /// Per-pole occlusion plus WHO stoppers it, when the dominant blocker
    /// is an ELECTRON: a cap rider doesn't deaden the intake, it CAPTURES
    /// and disperses it — the flow VFX shows that instead of nothing.
    pub fn pole_occlusion_blockers(&self) -> Vec<[(f64, Option<usize>); 2]> {
        let n = self.particles.len();
        let mut occ = vec![[(0.0f64, None); 2]; n];
        for i in 0..n {
            let pole = self.particles[i].pole_axis();
            let a = self.particles[i].position;
            for (side, dir) in [pole, -pole].into_iter().enumerate() {
                let mut worst = 0.0f64;
                let mut who: Option<usize> = None;
                for (k, other) in self.particles.iter().enumerate() {
                    if k == i {
                        continue;
                    }
                    let rel = other.position - a;
                    let along = rel.dot(dir);
                    if along <= 0.05 || along >= 3.0 {
                        continue;
                    }
                    let lat = (rel - dir * along).length();
                    let block = ((0.8 - lat) / 0.5).clamp(0.0, 1.0);
                    if block > worst {
                        worst = block;
                        who = (self.profiles[other.profile_id].mass < 0.5).then_some(k);
                    }
                }
                occ[i][side] = (worst, who);
            }
        }
        occ
    }

    /// Molecular bond detection: pairs of BARYONS standing at molecular
    /// range with both poles facing along the pair axis — the
    /// stream-cushion geometry H₂ settles into (diatom.pdf;
    /// h2_bond_matrix). Exclusions carry the physics:
    /// - same-group pairs are FUSED, not bonded (a molecule is "backed out
    ///   a distance" — ethane.pdf/ammon.pdf; fusion is a forcefully filled
    ///   plug — deut.pdf);
    /// - a stoppered channel is no bond (electron-between drives atoms
    ///   apart), so occluded pairs are excluded — which also makes
    ///   electron-capped helium poles correctly inert;
    /// - electrons don't count (their streams are 1/1836 — no cushion).
    pub fn molecular_bonds(&self) -> Vec<(usize, usize)> {
        const D_MIN: f64 = 2.5;
        const D_MAX: f64 = 8.0;
        /// |cos| of pole vs pair axis — poles must face each other.
        const ALIGN: f64 = 0.7;
        let n = self.particles.len();
        let mut out = Vec::new();
        if n < 2 {
            return out;
        }
        let occ = self.compute_occlusion();
        for i in 0..n {
            if self.profiles[self.particles[i].profile_id].mass < 0.5 {
                continue;
            }
            for j in (i + 1)..n {
                if self.profiles[self.particles[j].profile_id].mass < 0.5 {
                    continue;
                }
                if self.particles[i].group.is_some()
                    && self.particles[i].group == self.particles[j].group
                {
                    continue;
                }
                let d_vec = self.particles[j].position - self.particles[i].position;
                let d = d_vec.length();
                if !(D_MIN..=D_MAX).contains(&d) {
                    continue;
                }
                let d_hat = d_vec / d;
                if self.particles[i].pole_axis().dot(d_hat).abs() < ALIGN
                    || self.particles[j].pole_axis().dot(d_hat).abs() < ALIGN
                {
                    continue;
                }
                if occ[i * n + j] > 0.5 {
                    continue;
                }
                out.push((i, j));
            }
        }
        out
    }

    fn compute_forces(&mut self) {
        let n = self.particles.len();
        for p in &mut self.particles {
            p.force_accum = DVec3::ZERO;
            p.torque_accum = DVec3::ZERO;
        }

        let cq = self.couplings;
        let occlusion = self.compute_occlusion();

        // Pairwise forces
        for i in 0..n {
            for j in (i + 1)..n {
                // Rigid-composite constituents don't interact internally:
                // the nucleus is pre-fused (uf4.pdf — alphas can't be broken
                // and rearranged); its internal balance isn't simulated.
                if self.particles[i].group.is_some()
                    && self.particles[i].group == self.particles[j].group
                {
                    continue;
                }
                let d_vec = self.particles[j].position - self.particles[i].position;
                let r2 = d_vec.length_squared();
                let r = r2.sqrt().max(SOFTENING);
                let r2s = r * r;
                let r4s = r2s * r2s;
                let d_hat = d_vec / r;

                let pi_prof = &self.profiles[self.particles[i].profile_id];
                let pj_prof = &self.profiles[self.particles[j].profile_id];

                let pole_i = self.particles[i].pole_axis();
                let pole_j = self.particles[j].pole_axis();

                // θ_A: angle from A's pole to direction toward B
                let cos_theta_i = pole_i.dot(d_hat);
                // θ_B: angle from B's pole to direction toward A
                let cos_theta_j = pole_j.dot(-d_hat);

                let emission_i = pi_prof.emission.sample(cos_theta_i);
                let emission_j = pj_prof.emission.sample(cos_theta_j);
                let absorption_j = pj_prof.absorption.sample(cos_theta_j);
                let absorption_i = pi_prof.absorption.sample(cos_theta_i);

                // Doppler correction: approaching particles sweep through more
                // charge photons per unit time, retreating ones sweep fewer.
                // This asymmetry extracts kinetic energy on each close pass.
                let v_rel = self.particles[j].velocity - self.particles[i].velocity;
                let v_radial = v_rel.dot(d_hat);
                let doppler = (1.0 - cq.drag * v_radial).clamp(0.2, 5.0);

                // Gravity: G_q * m_i * m_j / r²  (attractive)
                let f_grav = cq.g_q * pi_prof.mass * pj_prof.mass / r2s;

                // Charge: C_q * m_emitter * m_receiver * E(θ) * R(θ) / r⁴ × doppler
                // Receiver mass models cross-section: electron intercepts 1/1836
                // the photons a proton would (Mathis's Dalton).
                let mass_prod = pi_prof.mass * pj_prof.mass;
                let f_charge_on_j =
                    cq.c_q * mass_prod * emission_i * absorption_j / r4s * doppler;
                let f_charge_on_i =
                    cq.c_q * mass_prod * emission_j * absorption_i / r4s * doppler;

                // Meridional emission-pressure gradient (transverse confinement).
                // The radial charge push above is the r-component of a photon
                // flux whose density varies with latitude; the θ-component
                // pushes DOWN the gradient — into the polar channel where
                // E(θ) ≈ 0. Two jobs: (a) keeps a polar orbiter from drifting
                // into the equatorial 1/r⁴ wall (the session-26 "transverse
                // confinement" gap), and (b) makes the charge force the
                // gradient of Φ ∝ E(θ)·R/3r³ — conservative, so close passes
                // no longer pump energy into the orbiter.
                // θ̂ = (d̂·cosθ − pole_nearest)/sinθ, singular only on-axis
                // where dE/dθ = 0 anyway.
                {
                    // i's emission gradient acting on j
                    let pole_eff_i = pole_i * cos_theta_i.signum();
                    let sin_i = (1.0 - cos_theta_i * cos_theta_i).max(0.0).sqrt();
                    if sin_i > 1e-6 {
                        let theta_hat = (d_hat * cos_theta_i.abs() - pole_eff_i) / sin_i;
                        let de = pi_prof.emission.d_dtheta(cos_theta_i);
                        let f_conf = -cq.c_q * mass_prod * absorption_j * de / (3.0 * r4s);
                        self.particles[j].force_accum += theta_hat * f_conf;
                        self.particles[i].force_accum -= theta_hat * f_conf;
                    }
                    // j's emission gradient acting on i
                    let pole_eff_j = pole_j * cos_theta_j.signum();
                    let sin_j = (1.0 - cos_theta_j * cos_theta_j).max(0.0).sqrt();
                    if sin_j > 1e-6 {
                        let theta_hat = (-d_hat * cos_theta_j.abs() - pole_eff_j) / sin_j;
                        let de = pj_prof.emission.d_dtheta(cos_theta_j);
                        let f_conf = -cq.c_q * mass_prod * absorption_i * de / (3.0 * r4s);
                        self.particles[i].force_accum += theta_hat * f_conf;
                        self.particles[j].force_accum -= theta_hat * f_conf;
                    }
                }

                // Channel occlusion for this pair (stoppered vortex).
                let occ = occlusion[i * n + j];

                // Ambient pressure: pushes particles into charge shadows.
                // A stoppered channel has no shadow minimum (diatom.pdf).
                let shadow = (1.0 - emission_i) * (1.0 - emission_j);
                let f_ambient = cq.ambient_pressure * shadow * (1.0 - occ) / r2s;

                // Polar intake: the proton recycles charge through its poles,
                // creating a focused inward flow.  A(θ)² sharpens the profile
                // to a narrow polar cone matching Mathis's vortex description.
                // A particle stoppering the channel blocks the stream.
                let ai2 = absorption_i * absorption_i;
                let aj2 = absorption_j * absorption_j;

                // Radial intake: pulls toward the emitter, 1/r² (flow sink).
                let f_intake_on_j =
                    cq.intake * pi_prof.mass * pj_prof.mass * ai2 * (1.0 - occ) / r2s;
                let f_intake_on_i =
                    cq.intake * pj_prof.mass * pi_prof.mass * aj2 * (1.0 - occ) / r2s;

                // Stream-collision cushion: opposing charge streams meet
                // head-on between the pair. Two symmetric channels:
                // - A_i²·A_j² — facing polar INTAKE streams (pole-to-pole
                //   charge meeting head-to-head — fourier.pdf, jup3.pdf);
                //   sets the H₂ bond standoff.
                // - E_i²·E_j² — facing EMISSION discs colliding edge-to-edge;
                //   this is why two bare protons repel to a standoff instead
                //   of gravitating into mutual orbit at ambient pressure.
                // (Pole-facing-disc mixes are NOT collisions — the disc blows
                // into the intake and feeds it, handled by intake/channeling.)
                // Pressure ∝ product of the two stream densities ⇒ 1/r⁴, and
                // ∝ (m_i·m_j)² since recycling throughput scales with mass
                // (the electron's 1/1836 stream is no cushion at all). A
                // stoppering particle back-scatters both streams — the
                // blocked channel pushes HARDER (rotor back-pressure):
                // ×(1+2·occ). This is what stands bonded atoms off at
                // molecular distance and drives electron-between atoms apart
                // (4-bond/4-repel matrix).
                let ei2 = emission_i * emission_i;
                let ej2 = emission_j * emission_j;
                let f_stream = cq.stream * (mass_prod * mass_prod)
                    * (ai2 * aj2 + ei2 * ej2)
                    * (1.0 + 2.0 * occ)
                    / r4s;

                // d_hat points from i to j.
                // Gravity pulls j toward i: along -d_hat
                // Charge pushes j away from i: along +d_hat
                // Ambient pulls j toward i: along -d_hat
                // Intake pulls j toward i: along -d_hat
                // Stream cushion pushes apart: along +d_hat for j
                let net_on_j =
                    (f_charge_on_j + f_stream - f_grav - f_ambient - f_intake_on_j) * d_hat;
                let net_on_i =
                    (f_grav + f_ambient + f_intake_on_i - f_charge_on_i - f_stream) * d_hat;

                self.particles[j].force_accum += net_on_j;
                self.particles[i].force_accum += net_on_i;

                // Channeling: the converging intake flow funnels particles
                // toward the emitter's pole axis.  Force ∝ A²·sin(θ)/r²,
                // directed toward the nearest point on the axis.
                if cq.intake.abs() > 1e-12 {
                    // Channel j toward i's pole axis
                    let along_i = d_vec.dot(pole_i);
                    let lateral_i = d_vec - pole_i * along_i;
                    let lat_dist_i = lateral_i.length();
                    if lat_dist_i > 1e-9 {
                        let toward_axis = -lateral_i / lat_dist_i;
                        let f_chan = cq.intake * pi_prof.mass
                            * pj_prof.mass * ai2 * lat_dist_i / (r * r2s);
                        self.particles[j].force_accum += toward_axis * f_chan;
                    }
                    // Channel i toward j's pole axis
                    let d_vec_ji = -d_vec;
                    let along_j = d_vec_ji.dot(pole_j);
                    let lateral_j = d_vec_ji - pole_j * along_j;
                    let lat_dist_j = lateral_j.length();
                    if lat_dist_j > 1e-9 {
                        let toward_axis = -lateral_j / lat_dist_j;
                        let f_chan = cq.intake * pj_prof.mass
                            * pi_prof.mass * aj2 * lat_dist_j / (r * r2s);
                        self.particles[i].force_accum += toward_axis * f_chan;
                    }

                }

                // Corotation drag: the polar charge vortex corotates with the
                // emitter's axial spin (diatom.pdf — captured electrons "circle
                // the drain or the pole"; the vortex is what holds them there).
                // Field angular velocity decays as (r_emitter/r)²; the drag
                // pulls the orbiter's tangential velocity toward the local
                // field velocity, capped at COROT_V_MAX (sub-c). Self-limiting:
                // the force vanishes as v_tan → v_corot, so unlike the removed
                // session-27 intake vortex it cannot accelerate past corotation
                // (no slingshot). Gated by A²(θ) to the polar funnel.
                if cq.corot.abs() > 1e-12 {
                    let v_tan_vec = v_rel - d_hat * v_radial;

                    // i's vortex acting on j
                    let spin_i = self.particles[i].angular_velocity.dot(pole_i);
                    let omega_i = spin_i * (pi_prof.radius / r).powi(2);
                    let mut v_corot = pole_i.cross(d_vec) * omega_i;
                    let vc_len = v_corot.length();
                    if vc_len > COROT_V_MAX {
                        v_corot *= COROT_V_MAX / vc_len;
                    }
                    let f_corot =
                        (v_corot - v_tan_vec) * (cq.corot * mass_prod * ai2 / r2s);
                    self.particles[j].force_accum += f_corot;
                    self.particles[i].force_accum -= f_corot;

                    // j's vortex acting on i (relative velocities negate)
                    let spin_j = self.particles[j].angular_velocity.dot(pole_j);
                    let omega_j = spin_j * (pj_prof.radius / r).powi(2);
                    let mut v_corot_j = pole_j.cross(-d_vec) * omega_j;
                    let vcj_len = v_corot_j.length();
                    if vcj_len > COROT_V_MAX {
                        v_corot_j *= COROT_V_MAX / vcj_len;
                    }
                    let f_corot_i =
                        (v_corot_j + v_tan_vec) * (cq.corot * mass_prod * aj2 / r2s);
                    self.particles[i].force_accum += f_corot_i;
                    self.particles[j].force_accum -= f_corot_i;
                }

                // Vortex force: spinning emission carries tangential momentum.
                // Photons leave the surface at the spin velocity, creating
                // a "sprinkler" drag that captures nearby particles into orbit.
                if cq.vortex.abs() > 1e-12 {
                    let spin_i = self.particles[i].angular_velocity.dot(pole_i);
                    let sin_vec_i = pole_i.cross(d_hat);
                    let sin_theta_i = sin_vec_i.length();
                    if sin_theta_i > 1e-9 {
                        let tangent_i = sin_vec_i / sin_theta_i;
                        let f_vort = cq.vortex * mass_prod
                            * emission_i * absorption_j
                            * spin_i * sin_theta_i / r4s;
                        self.particles[j].force_accum += tangent_i * f_vort;
                    }

                    let spin_j = self.particles[j].angular_velocity.dot(pole_j);
                    let sin_vec_j = pole_j.cross(-d_hat);
                    let sin_theta_j = sin_vec_j.length();
                    if sin_theta_j > 1e-9 {
                        let tangent_j = sin_vec_j / sin_theta_j;
                        let f_vort = cq.vortex * mass_prod
                            * emission_j * absorption_i
                            * spin_j * sin_theta_j / r4s;
                        self.particles[i].force_accum += tangent_j * f_vort;
                    }
                }

                // Contact repulsion: hard-sphere boundary at sum of radii,
                // with near-critical normal damping while approaching —
                // an undamped spring against the 1/1836-mass electron
                // (dt·ω ≈ 0.2) scatters it chaotically ("bunny hopping").
                let ri = pi_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let rj = pj_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let r_contact = ri + rj;
                if r < r_contact {
                    let overlap = r_contact - r;
                    let mut f_contact = CONTACT_STIFFNESS * overlap;
                    if v_radial < 0.0 {
                        let m_red = mass_prod / (pi_prof.mass + pj_prof.mass);
                        let c_n = 2.0 * (CONTACT_STIFFNESS * m_red).sqrt();
                        f_contact -= c_n * v_radial; // approaching ⇒ extra repulsion
                    }
                    self.particles[j].force_accum += d_hat * f_contact;
                    self.particles[i].force_accum -= d_hat * f_contact;
                }

                // Torque: aligns equator toward charge source.
                // Uses -sin(2θ) restoring form: stable with equator facing charge.
                // torque = -2 * dot(pole, d_hat) * cross(pole, d_hat) * |f| * coupling
                let torque_j = pole_j.cross(-d_hat)
                    * (-2.0 * pole_j.dot(-d_hat) * f_charge_on_j * cq.torque);
                let torque_i = pole_i.cross(d_hat)
                    * (-2.0 * pole_i.dot(d_hat) * f_charge_on_i * cq.torque);

                self.particles[j].torque_accum += torque_j;
                self.particles[i].torque_accum += torque_i;
            }
        }

        // Ambient environmental forces
        for p in &mut self.particles {
            let m = self.profiles[p.profile_id].mass;
            p.force_accum += self.ambient_gravity * m;
            p.force_accum += self.ambient_charge * m;
        }
    }

    // ── Rendering buffers (pure data, converted by the Godot wrapper) ──

    /// Meridional cross-section ring for a particle's emission or absorption
    /// profile. World-space vertices for a closed LINE_STRIP.
    /// `plane` 0 or 1 selects two orthogonal meridional planes.
    /// `power` raises the profile value (e.g. 2.0 for A²).
    pub fn profile_ring(
        &self,
        particle_id: usize,
        use_absorption: bool,
        base_radius: f32,
        scale: f32,
        num_points: usize,
        plane: usize,
        power: f32,
    ) -> Vec<DVec3> {
        if particle_id >= self.particles.len() {
            return Vec::new();
        }
        let p = &self.particles[particle_id];
        let prof = &self.profiles[p.profile_id];
        let table = if use_absorption { &prof.absorption } else { &prof.emission };
        let center = p.position;
        let pole = p.pole_axis();
        let (right, forward) = build_frame(pole);
        let perp = if plane == 0 { right } else { forward };
        let n = num_points.max(12);
        let mut out = Vec::with_capacity(n + 1);
        for i in 0..=n {
            let alpha = (i as f64 / n as f64) * TAU;
            let dir = pole * alpha.cos() + perp * alpha.sin();
            let cos_theta = alpha.cos();
            let val = (table.sample(cos_theta) as f32).powf(power);
            let r = base_radius + val * scale;
            out.push(center + dir * r as f64);
        }
        out
    }

    /// 16 floats per particle: 12 (3x4 transform) + 4 (RGBA color).
    /// Row-major with interleaved origin (Godot MultiMesh format).
    pub fn build_multimesh_buffer(&self) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::with_capacity(self.particles.len() * 16);
        for p in &self.particles {
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            push_transform_color(&mut buf, p, scale, profile_color(&profile.name), 1.0);
        }
        buf
    }

    /// Transform+color buffer for particles of one profile type only.
    pub fn build_multimesh_buffer_for_profile(&self, profile_id: usize) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::new();
        for p in &self.particles {
            if p.profile_id != profile_id {
                continue;
            }
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            push_transform_color(&mut buf, p, scale, (1.0, 1.0, 1.0), 1.0);
        }
        buf
    }

    /// Instance buffer for the field-extent SKINS: FREE particles only.
    /// Rigid-group constituents are fused — pushed inside each other's
    /// free-field reach and locked there (nuclear.pdf: "great forces
    /// pushing baryons into configurations they couldn't otherwise
    /// achieve"), recycling as one unit. Their individual free-field skins
    /// don't exist any more, for the same reason intra-group pair forces
    /// are skipped in compute_forces.
    pub fn build_skin_multimesh_buffer_for_profile(&self, profile_id: usize) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::new();
        for p in &self.particles {
            if p.profile_id != profile_id || p.group.is_some() {
                continue;
            }
            let profile = &self.profiles[p.profile_id];
            let scale = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            push_transform_color(&mut buf, p, scale, (1.0, 1.0, 1.0), 1.0);
        }
        buf
    }

    /// Line buffer for pole axis indicators.
    /// 2 vertices per particle (center and pole tip), 6 floats each (pos xyz + color rgb).
    pub fn build_pole_indicator_buffer(&self) -> Vec<f32> {
        let mut buf: Vec<f32> = Vec::with_capacity(self.particles.len() * 12);
        for p in &self.particles {
            let profile = &self.profiles[p.profile_id];
            let r = (profile.radius as f32).max(MIN_RENDER_RADIUS);
            let pos = p.position;
            let pole = p.pole_axis();
            let tip = pos + pole * (r as f64 * 1.5);

            buf.extend_from_slice(&[
                pos.x as f32, pos.y as f32, pos.z as f32,
                1.0, 0.9, 0.2,
            ]);
            buf.extend_from_slice(&[
                tip.x as f32, tip.y as f32, tip.z as f32,
                1.0, 0.4, 0.1,
            ]);
        }
        buf
    }

    /// Surface-of-revolution mesh from a profile's emission table.
    /// Packed: [vert_count, idx_count, verts...(x,y,z,nx,ny,nz,r,g,b,a), indices...]
    pub fn build_profile_mesh(
        &self,
        profile_id: usize,
        lon_segments: usize,
        lat_segments: usize,
    ) -> Vec<f32> {
        let prof = match self.profiles.get(profile_id) {
            Some(p) => p,
            None => return Vec::new(),
        };

        let lon = lon_segments.max(8);
        let lat = lat_segments.max(4);
        let rings = lat * 2 + 1;
        let verts_per_ring = lon + 1;
        let total_verts = rings * verts_per_ring;
        let quad_count = (rings - 1) * lon;
        let total_indices = quad_count * 6;

        let floats_per_vert = 10; // pos(3) + normal(3) + color(4)
        let mut buf: Vec<f32> =
            Vec::with_capacity(2 + total_verts * floats_per_vert + total_indices);

        buf.push(total_verts as f32);
        buf.push(total_indices as f32);

        let (type_r, type_g, type_b) = profile_color(&prof.name);

        for lat_idx in 0..rings {
            let theta = std::f64::consts::PI * lat_idx as f64 / (rings - 1).max(1) as f64;
            let cos_theta = theta.cos();
            let sin_theta = theta.sin();

            let emission = prof.emission.sample(cos_theta);
            // Lathe approach: fixed height (y = cos θ), emission modulates width only
            let r_cross = sin_theta * (ENVELOPE_MIN_R + emission * (1.0 - ENVELOPE_MIN_R));
            let y_pos = cos_theta;

            let intensity = (0.3f32 + emission as f32 * 0.7).min(1.0);
            let alpha = 0.25 + emission as f32 * 0.25;

            for lon_idx in 0..=lon {
                let phi = TAU * lon_idx as f64 / lon as f64;
                let x = (r_cross * phi.cos()) as f32;
                let y = y_pos as f32;
                let z = (r_cross * phi.sin()) as f32;

                buf.extend_from_slice(&[x, y, z]);

                let len = (x * x + y * y + z * z).sqrt();
                if len > 1e-6 {
                    buf.extend_from_slice(&[x / len, y / len, z / len]);
                } else {
                    buf.extend_from_slice(&[0.0, if lat_idx == 0 { 1.0 } else { -1.0 }, 0.0]);
                }

                buf.extend_from_slice(&[
                    type_r * intensity,
                    type_g * intensity,
                    type_b * intensity,
                    alpha,
                ]);
            }
        }

        for lat_idx in 0..(rings - 1) {
            for lon_idx in 0..lon {
                let tl = (lat_idx * verts_per_ring + lon_idx) as f32;
                let tr = tl + 1.0;
                let bl = ((lat_idx + 1) * verts_per_ring + lon_idx) as f32;
                let br = bl + 1.0;
                buf.extend_from_slice(&[tl, bl, tr, tr, bl, br]);
            }
        }

        buf
    }

    /// Field-extent skin: a translucent surface of revolution whose RADIUS
    /// at each latitude is the REACH of this profile's charge push — the
    /// distance where the emission force on a unit test absorber decays to
    /// the natural force unit (gravity of two unit masses at r=1):
    ///     C_q·m·E(θ)/r⁴ = 1  ⇒  r(θ) = (C_q·m·E(θ))^¼.
    /// Opacity carries strength (alpha ∝ E(θ)/E_max). Together these
    /// restore the wireframe rings' information — how strong, how far — as
    /// a soft envelope: a proton reads as a wide equatorial ledge
    /// (~4.7 natural units at C_q=500) pinching down to the body at the
    /// polar holes, where the skin goes transparent and the intake tornado
    /// shows through.
    /// Mesh-local units are PARTICLE RADII (the instance transform scales
    /// by profile radius, same as the body mesh). Same packing as
    /// `build_profile_mesh`:
    /// [vert_count, idx_count, verts(x,y,z,nx,ny,nz,r,g,b,a)…, indices…].
    pub fn build_field_skin_mesh(
        &self,
        profile_id: usize,
        lon_segments: usize,
        lat_segments: usize,
    ) -> Vec<f32> {
        let prof = match self.profiles.get(profile_id) {
            Some(p) => p,
            None => return Vec::new(),
        };

        let lon = lon_segments.max(8);
        let lat = lat_segments.max(4);
        let rings = lat * 2 + 1;
        let verts_per_ring = lon + 1;
        let total_verts = rings * verts_per_ring;
        let quad_count = (rings - 1) * lon;
        let total_indices = quad_count * 6;

        let mut buf: Vec<f32> = Vec::with_capacity(2 + total_verts * 10 + total_indices);
        buf.push(total_verts as f32);
        buf.push(total_indices as f32);

        let (type_r, type_g, type_b) = profile_color(&prof.name);
        let render_r = prof.radius.max(MIN_RENDER_RADIUS as f64);
        let e_max = prof
            .emission
            .bins
            .iter()
            .fold(0.0f64, |m, &v| m.max(v as f64))
            .max(1e-6);
        let cq_m = self.couplings.c_q * prof.mass;

        for lat_idx in 0..rings {
            let theta = std::f64::consts::PI * lat_idx as f64 / (rings - 1).max(1) as f64;
            let cos_theta = theta.cos();
            let sin_theta = theta.sin();

            let emission = prof.emission.sample(cos_theta);
            // Reach in natural units, clamped to hug the body where the
            // emission vanishes (the polar holes), then converted to
            // mesh-local units.
            let reach = (cq_m * emission).max(0.0).powf(0.25).max(render_r * 1.06);
            let rho = reach / render_r;

            let strength = (emission / e_max).clamp(0.0, 1.0) as f32;
            let alpha = 0.03 + 0.20 * strength;

            for lon_idx in 0..=lon {
                let phi = TAU * lon_idx as f64 / lon as f64;
                let nx = sin_theta * phi.cos();
                let ny = cos_theta;
                let nz = sin_theta * phi.sin();

                buf.extend_from_slice(&[
                    (rho * nx) as f32,
                    (rho * ny) as f32,
                    (rho * nz) as f32,
                ]);
                buf.extend_from_slice(&[nx as f32, ny as f32, nz as f32]);
                buf.extend_from_slice(&[
                    0.4 + 0.6 * type_r,
                    0.4 + 0.6 * type_g,
                    0.4 + 0.6 * type_b,
                    alpha,
                ]);
            }
        }

        for lat_idx in 0..(rings - 1) {
            for lon_idx in 0..lon {
                let tl = (lat_idx * verts_per_ring + lon_idx) as f32;
                let tr = tl + 1.0;
                let bl = ((lat_idx + 1) * verts_per_ring + lon_idx) as f32;
                let br = bl + 1.0;
                buf.extend_from_slice(&[tl, bl, tr, tr, bl, br]);
            }
        }

        buf
    }

    /// March the group's summed charge push outward to find its reach
    /// surface (the composite-skin boundary): azimuth-averaged reach per
    /// polar angle, group-local frame, never dipping inside the member
    /// bodies' silhouette. Called once at spawn (the group is rigid) and
    /// stored as `RigidGroup::skin_reach`.
    pub fn compute_group_reach(&self, group_idx: usize, rings: usize) -> Vec<f64> {
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        struct Src {
            pos: DVec3,
            pole: DVec3,
            mass: f64,
            radius: f64,
            profile: usize,
        }
        let srcs: Vec<Src> = g
            .members
            .iter()
            .enumerate()
            .map(|(k, &pid)| {
                let profile = self.particles[pid].profile_id;
                Src {
                    pos: g.local_offsets[k],
                    pole: g.local_orients[k] * DVec3::Y,
                    mass: self.profiles[profile].mass,
                    radius: self.profiles[profile].radius.max(MIN_RENDER_RADIUS as f64),
                    profile,
                }
            })
            .collect();

        let cq = self.couplings.c_q;
        let push = |x: DVec3| -> f64 {
            let mut f = DVec3::ZERO;
            for s in &srcs {
                let d_vec = x - s.pos;
                let r = d_vec.length().max(SOFTENING);
                let d_hat = d_vec / r;
                let e = self.profiles[s.profile].emission.sample(s.pole.dot(d_hat));
                f += d_hat * (cq * s.mass * e / (r * r * r * r));
            }
            f.length()
        };

        const F_REF: f64 = 1.0;
        const DR: f64 = 0.08;
        const R_MAX: f64 = 40.0;
        const AZ_SAMPLES: usize = 8;

        let rings = rings.max(2);
        let mut reach = Vec::with_capacity(rings);
        for lat_idx in 0..rings {
            let theta = std::f64::consts::PI * lat_idx as f64 / (rings - 1) as f64;
            let (sin_t, cos_t) = (theta.sin(), theta.cos());
            let mut sum = 0.0;
            for az in 0..AZ_SAMPLES {
                let phi = TAU * az as f64 / AZ_SAMPLES as f64;
                let dir = DVec3::new(sin_t * phi.cos(), cos_t, sin_t * phi.sin());
                let hull = srcs
                    .iter()
                    .map(|s| s.pos.dot(dir) + s.radius)
                    .fold(0.3f64, f64::max);
                let mut r = hull;
                while r < R_MAX && push(dir * r) >= F_REF {
                    r += DR;
                }
                sum += r;
            }
            reach.push(sum / AZ_SAMPLES as f64);
        }
        reach
    }

    /// Composite nucleus skin: the reach envelope of a WHOLE rigid group —
    /// the surface where the group's SUMMED charge push on a unit absorber
    /// falls to the natural force unit (same F=1 reference as the free
    /// skin), from the group's stored `skin_reach` table. Fused
    /// constituents draw no individual skins (they recycle as one unit);
    /// this is the field the locked configuration projects instead.
    /// GROUP-LOCAL units — render it riding the group transform, unscaled.
    /// Gold, to read as "one fused unit" against the per-type skins of
    /// free particles. Same packing as `build_profile_mesh`.
    pub fn build_group_skin_mesh(
        &self,
        group_idx: usize,
        lon_segments: usize,
        lat_segments: usize,
    ) -> Vec<f32> {
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        let lon = lon_segments.max(8);
        let lat = lat_segments.max(4);
        let rings = lat * 2 + 1;
        let verts_per_ring = lon + 1;
        let total_verts = rings * verts_per_ring;
        let total_indices = (rings - 1) * lon * 6;

        let table = if g.skin_reach.len() >= 2 {
            g.skin_reach.clone()
        } else {
            self.compute_group_reach(group_idx, SKIN_REACH_RINGS)
        };
        let reach: Vec<f64> = (0..rings)
            .map(|lat_idx| {
                let theta =
                    std::f64::consts::PI * lat_idx as f64 / (rings - 1).max(1) as f64;
                reach_at_theta(&table, theta)
            })
            .collect();
        let reach_max = reach.iter().fold(1e-6f64, |m, &v| v.max(m));

        let mut buf: Vec<f32> = Vec::with_capacity(2 + total_verts * 10 + total_indices);
        buf.push(total_verts as f32);
        buf.push(total_indices as f32);

        const GOLD: (f32, f32, f32) = (0.95, 0.80, 0.50);
        for (lat_idx, &rho) in reach.iter().enumerate() {
            let theta = std::f64::consts::PI * lat_idx as f64 / (rings - 1).max(1) as f64;
            let (sin_t, cos_t) = (theta.sin(), theta.cos());
            // reach ∝ F^¼ ⇒ reach⁴ recovers the field-strength analog the
            // free skin maps to opacity.
            let strength = ((rho / reach_max).powi(4)).clamp(0.0, 1.0) as f32;
            let alpha = 0.03 + 0.20 * strength;
            for lon_idx in 0..=lon {
                let phi = TAU * lon_idx as f64 / lon as f64;
                let nx = sin_t * phi.cos();
                let ny = cos_t;
                let nz = sin_t * phi.sin();
                buf.extend_from_slice(&[
                    (rho * nx) as f32,
                    (rho * ny) as f32,
                    (rho * nz) as f32,
                ]);
                buf.extend_from_slice(&[nx as f32, ny as f32, nz as f32]);
                buf.extend_from_slice(&[GOLD.0, GOLD.1, GOLD.2, alpha]);
            }
        }

        for lat_idx in 0..(rings - 1) {
            for lon_idx in 0..lon {
                let tl = (lat_idx * verts_per_ring + lon_idx) as f32;
                let tr = tl + 1.0;
                let bl = ((lat_idx + 1) * verts_per_ring + lon_idx) as f32;
                let br = bl + 1.0;
                buf.extend_from_slice(&[tl, bl, tr, tr, bl, br]);
            }
        }

        buf
    }

    /// Where each alpha block's disc meets the composite reach surface:
    /// (ring radius, y) per alpha, group-local. Computed once at spawn and
    /// stored as `RigidGroup::disc_exits` — feeds the ring overlay and the
    /// flow-VFX disc exit paths.
    pub fn compute_disc_exits(&self, group_idx: usize) -> Vec<(f64, f64)> {
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        if g.skin_reach.len() < 2 {
            return Vec::new();
        }
        // Alpha disc latitudes: cluster the AXIAL protons' heights (each
        // alpha = two stacked protons ±0.9 around its center; plug protons
        // are edge-on and don't count).
        let mut ys: Vec<f64> = Vec::new();
        for (k, &pid) in g.members.iter().enumerate() {
            let prof = &self.profiles[self.particles[pid].profile_id];
            let axial = (g.local_orients[k] * DVec3::Y).dot(DVec3::Y).abs() > 0.9;
            if prof.name == "proton" && axial {
                ys.push(g.local_offsets[k].y);
            }
        }
        if ys.is_empty() {
            return Vec::new();
        }
        ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
        // Axial protons come as one stacked pair per alpha block, and
        // blocks don't interleave — consecutive sorted heights pair up.
        // (A gap-based clustering fails here: adjacent-block protons are
        // CLOSER (0.8) than an alpha's own pair (1.8) at ALPHA_PITCH 2.6.)
        let mut centers: Vec<f64> = Vec::new();
        let mut i = 0;
        while i < ys.len() {
            if i + 1 < ys.len() {
                centers.push((ys[i] + ys[i + 1]) / 2.0);
                i += 2;
            } else {
                centers.push(ys[i]);
                i += 1;
            }
        }

        // For each alpha center find where the skin surface crosses that
        // height: y(θ) = reach(θ)·cosθ runs from +reach to −reach over
        // θ ∈ [0, π] — take the first crossing.
        let table = &g.skin_reach;
        let n = table.len();
        let mut rings: Vec<(f64, f64)> = Vec::new(); // (ring radius, y)
        for &yc in &centers {
            for i in 0..(n - 1) {
                let th0 = std::f64::consts::PI * i as f64 / (n - 1) as f64;
                let th1 = std::f64::consts::PI * (i + 1) as f64 / (n - 1) as f64;
                let y0 = table[i] * th0.cos();
                let y1 = table[i + 1] * th1.cos();
                if (y0 - yc) * (y1 - yc) <= 0.0 && (y0 - y1).abs() > 1e-12 {
                    let f = ((y0 - yc) / (y0 - y1)).clamp(0.0, 1.0);
                    let th = th0 + (th1 - th0) * f;
                    let r = reach_at_theta(table, th);
                    rings.push((r * th.sin(), r * th.cos()));
                    break;
                }
            }
        }
        rings
    }

    /// Max-emission ring overlay for the composite skin: one circle per
    /// ALPHA BLOCK, drawn ON the reach surface at that alpha's disc
    /// latitude — "this is where this alpha's disc pushes farthest", the
    /// quantitative line the old wireframe rings carried. Group-local
    /// units, same frame as the skin mesh.
    /// Packed: [ring_count, pts_per_ring, ring_count·pts_per_ring × (x,y,z)].
    pub fn build_group_emission_rings(&self, group_idx: usize) -> Vec<f32> {
        const PTS: usize = 65;
        let rings = match self.groups.get(group_idx) {
            Some(g) if !g.disc_exits.is_empty() => g.disc_exits.clone(),
            Some(_) => self.compute_disc_exits(group_idx),
            None => return Vec::new(),
        };
        if rings.is_empty() {
            return Vec::new();
        }
        let mut buf = Vec::with_capacity(2 + rings.len() * PTS * 3);
        buf.push(rings.len() as f32);
        buf.push(PTS as f32);
        for (rr, y) in rings {
            for p in 0..PTS {
                let phi = TAU * p as f64 / (PTS - 1) as f64;
                buf.extend_from_slice(&[
                    (rr * phi.cos()) as f32,
                    y as f32,
                    (rr * phi.sin()) as f32,
                ]);
            }
        }
        buf
    }

    pub fn count_particles_with_profile(&self, profile_id: usize) -> usize {
        self.particles
            .iter()
            .filter(|p| p.profile_id == profile_id)
            .count()
    }

    // ── Charge cloud VFX ─────────────────────────────────────────────────
    // Charge-recycling FLOW visualization: every dot is a photon parcel on
    // a COMPLETE path — pulled down the polar funnel, through the body (or
    // the whole nuclear stack), and back out where the emission maps send
    // it: the equatorial disc for protons, the far pole for neutrons and
    // through-charge, an alpha's max-emission ring for nuclei, or captured
    // and dispersed by a cap electron riding the pole. Paths are kinematic
    // illustrations of the force model (probability maps + composite
    // reach), not collision-checked trajectories — spin mode does that.

    /// Advance the cloud particles by wall-clock `delta` seconds and return
    /// the MultiMesh buffer (12 transform + 4 color floats per particle).
    pub fn advance_clouds(&mut self, delta: f64) -> Vec<f32> {
        const MAX_POOL: usize = 12_000;
        /// Parcel speed along its path (natural units / wall second).
        const FLOW_SPEED: f64 = 2.0;
        /// Recycling loops per second through a free baryon (split across
        /// the two poles). Electrons recycle a sliver of this.
        const FLOW_FREE_PER_SEC: f64 = 260.0;
        /// Loops per second per alpha engine in a nucleus — more engines
        /// in the stack, denser flow.
        const FLOW_GROUP_PER_ENGINE: f64 = 130.0;
        /// Destination weight for passing straight through the stack vs
        /// peeling off at one alpha disc (weight 1 each).
        const THROUGH_WEIGHT: f64 = 0.7;
        /// Visible swirl of the whole path about the pole (rad/s wall
        /// clock) — far below the physical spin rate so CW/CCW reads.
        const FLOW_SWIRL: f64 = 1.6;
        const INTAKE_COLOR: (f32, f32, f32) = (0.30, 0.65, 1.0);
        const THROUGH_COLOR: (f32, f32, f32) = (0.75, 0.90, 1.00);
        const ELECTRON_COLOR: (f32, f32, f32) = (0.20, 0.90, 0.35);
        const PROTON_COLOR: (f32, f32, f32) = (0.92, 0.30, 0.20);

        if !self.vfx_enabled || self.particles.is_empty() {
            self.vfx_particles.clear();
            return Vec::new();
        }

        self.vfx_time += delta;

        // Per-frame anchor frames. Particle frames come from the pole only
        // (build_frame) — the physical TAU·3 axial spin is excluded and the
        // readable FLOW_SWIRL is applied in the path eval instead.
        struct Emitter {
            pos: DVec3,
            rot: DQuat,
            radius: f64,
            color: (f32, f32, f32),
            pid: usize,
            spin_sign: f64,
            group: Option<usize>,
        }
        let emitters: Vec<Emitter> = self
            .particles
            .iter()
            .map(|p| {
                let prof = &self.profiles[p.profile_id];
                let pole = p.pole_axis();
                let (right, forward) = build_frame(pole);
                Emitter {
                    pos: p.position,
                    rot: DQuat::from_mat3(&DMat3::from_cols(right, pole, forward)),
                    radius: prof.radius.max(MIN_RENDER_RADIUS as f64),
                    color: profile_color(&prof.name),
                    pid: p.profile_id,
                    spin_sign: p.angular_velocity.dot(pole).signum(),
                    group: p.group,
                }
            })
            .collect();
        let group_frames: Vec<(DVec3, DQuat)> = self
            .groups
            .iter()
            .map(|g| (g.com, g.orientation))
            .collect();

        // Age + kinematic path evaluation.
        for vp in &mut self.vfx_particles {
            vp.age += delta;
            match vp.kind {
                VfxKind::Smoke => vp.position += vp.velocity * delta,
                VfxKind::Flow {
                    ref anchor,
                    swirl,
                    ref pts,
                    ..
                } => {
                    let frame = match *anchor {
                        VfxAnchor::Particle(i) => emitters.get(i).map(|e| (e.pos, e.rot)),
                        VfxAnchor::Group(gi) => group_frames.get(gi).copied(),
                    };
                    let Some((origin, rot)) = frame else {
                        vp.age = vp.lifetime + 1.0; // anchor gone — cull
                        continue;
                    };
                    let u = (vp.age / vp.lifetime).clamp(0.0, 1.0);
                    let t = u * (pts.len() - 1) as f64;
                    let i = (t as usize).min(pts.len().saturating_sub(2));
                    let f = t - i as f64;
                    let local = pts[i].lerp(pts[i + 1], f);
                    // Whole-path swirl about the pole.
                    let a = swirl * vp.age;
                    let (s, c) = a.sin_cos();
                    let sw = DVec3::new(
                        local.x * c + local.z * s,
                        local.y,
                        -local.x * s + local.z * c,
                    );
                    vp.position = origin + rot * sw;
                }
            }
        }
        self.vfx_particles.retain(|vp| vp.age < vp.lifetime);

        // Fractional spawn counts via random rounding.
        let spawn_count = |rng: &mut Rng, rate: f64| -> usize {
            let x = rate * delta;
            let base = x.floor() as usize;
            base + usize::from(rng.next_f64() < x.fract())
        };

        let occ = self.pole_occlusion_blockers();

        // FREE particles: complete recycling loops — in at an open pole,
        // through the body, out on the far hemisphere at the angle the
        // emission CDF picks (protons → the disc, neutrons → the far pole,
        // electrons → their tiny flat disc).
        for (ei, em) in emitters.iter().enumerate() {
            if em.group.is_some() {
                continue; // fused members flow as their GROUP, below
            }
            let mass = self.profiles[em.pid].mass;
            let rate_scale = if mass < 0.5 { 0.2 } else { 1.0 };
            for (side, &(occ_v, blocker)) in occ[ei].iter().enumerate() {
                let entry = if side == 0 { 1.0 } else { -1.0 };
                let n = spawn_count(
                    &mut self.rng,
                    FLOW_FREE_PER_SEC * 0.5 * (1.0 - occ_v) * rate_scale,
                );
                for _ in 0..n {
                    if self.vfx_particles.len() >= MAX_POOL {
                        break;
                    }
                    let theta_in =
                        self.profiles[em.pid].intake_cdf.sample(self.rng.next_f64());
                    let theta_out =
                        self.profiles[em.pid].emission_cdf.sample(self.rng.next_f64());
                    let reach = (self.couplings.c_q
                        * mass
                        * self.profiles[em.pid].emission.sample(theta_out.cos()))
                    .max(0.0)
                    .powf(0.25)
                    .max(em.radius * 1.5);
                    let (pts, len) = flow_free_loop(
                        &mut self.rng,
                        em.radius,
                        entry,
                        theta_in,
                        theta_out,
                        reach,
                        em.spin_sign,
                    );
                    self.vfx_particles.push(VfxParticle {
                        position: em.pos + em.rot * pts[0],
                        velocity: DVec3::ZERO,
                        age: 0.0,
                        lifetime: len / FLOW_SPEED * (0.85 + 0.3 * self.rng.next_f64()),
                        color: INTAKE_COLOR,
                        base_scale: (0.016 * em.radius.max(0.35)) as f32,
                        grow: 0.0,
                        kind: VfxKind::Flow {
                            anchor: VfxAnchor::Particle(ei),
                            swirl: em.spin_sign * FLOW_SWIRL,
                            color_end: em.color,
                            pts,
                        },
                    });
                }
                // Electron-capped pole: the inflow doesn't die, it gets
                // CAPTURED at the rider and dispersed off its disc.
                if blocker.is_some() && occ_v > 0.5 {
                    let n = spawn_count(&mut self.rng, FLOW_FREE_PER_SEC * 0.3 * rate_scale);
                    for _ in 0..n {
                        if self.vfx_particles.len() >= MAX_POOL {
                            break;
                        }
                        let (pts, len) =
                            flow_capture(&mut self.rng, entry * em.radius, entry, em.radius);
                        self.vfx_particles.push(VfxParticle {
                            position: em.pos + em.rot * pts[0],
                            velocity: DVec3::ZERO,
                            age: 0.0,
                            lifetime: len / FLOW_SPEED * (0.85 + 0.3 * self.rng.next_f64()),
                            color: INTAKE_COLOR,
                            base_scale: (0.016 * em.radius.max(0.35)) as f32,
                            grow: 0.0,
                            kind: VfxKind::Flow {
                                anchor: VfxAnchor::Particle(ei),
                                swirl: em.spin_sign * FLOW_SWIRL,
                                color_end: ELECTRON_COLOR,
                                pts,
                            },
                        });
                    }
                }
            }
        }

        // NUCLEI: one engine per alpha block. Parcels enter an OPEN stack
        // end, ride the axial channel, and either peel off at an alpha's
        // max-emission ring (weight 1 per alpha) or pass all the way
        // through and out the far pole (THROUGH_WEIGHT) — more engines,
        // proportionally denser flow.
        for gi in 0..self.groups.len() {
            let (com, g_rot) = group_frames[gi];
            let axis = g_rot * DVec3::Y;
            let (mut y_top, mut y_bot) = (f64::MIN, f64::MAX);
            for (k, &m) in self.groups[gi].members.iter().enumerate() {
                let r = self.profiles[self.particles[m].profile_id].radius;
                let y = self.groups[gi].local_offsets[k].y;
                y_top = y_top.max(y + r);
                y_bot = y_bot.min(y - r);
            }
            let exits = self.groups[gi].disc_exits.clone();
            let engines = exits.len().max(1);
            let reach_n = reach_at_theta(&self.groups[gi].skin_reach, 0.0);
            let reach_s = reach_at_theta(&self.groups[gi].skin_reach, std::f64::consts::PI);
            for side in 0..2 {
                let (entry, tip_y, far_tip, reach_out) = if side == 0 {
                    (1.0, y_top, y_bot, reach_s)
                } else {
                    (-1.0, y_bot, y_top, reach_n)
                };
                // End occlusion: a particle sitting over this stack end
                // stoppers it; an electron rider captures instead.
                let tip_world = com + axis * tip_y;
                let dir_out = axis * entry;
                let mut occ_v = 0.0f64;
                let mut blocker: Option<usize> = None;
                for (k, p) in self.particles.iter().enumerate() {
                    if p.group == Some(gi) {
                        continue;
                    }
                    let rel = p.position - tip_world;
                    let along = rel.dot(dir_out);
                    if along <= 0.05 || along >= 3.0 {
                        continue;
                    }
                    let lat = (rel - dir_out * along).length();
                    let block = ((0.8 - lat) / 0.5).clamp(0.0, 1.0);
                    if block > occ_v {
                        occ_v = block;
                        blocker = (self.profiles[p.profile_id].mass < 0.5).then_some(k);
                    }
                }
                let end_rate = FLOW_GROUP_PER_ENGINE * engines as f64 * 0.5;
                let n = spawn_count(&mut self.rng, end_rate * (1.0 - occ_v));
                for _ in 0..n {
                    if self.vfx_particles.len() >= MAX_POOL {
                        break;
                    }
                    let w_total = exits.len() as f64 + THROUGH_WEIGHT;
                    let pick = self.rng.next_f64() * w_total;
                    let (pts, len, color_end) = if pick >= exits.len() as f64 {
                        let (pts, len) = flow_group_through(
                            &mut self.rng, entry, tip_y, far_tip, reach_out,
                        );
                        (pts, len, THROUGH_COLOR)
                    } else {
                        let (rr, ry) = exits[pick as usize];
                        let (pts, len) =
                            flow_group_disc(&mut self.rng, entry, tip_y, rr, ry);
                        (pts, len, PROTON_COLOR)
                    };
                    self.vfx_particles.push(VfxParticle {
                        position: com + g_rot * pts[0],
                        velocity: DVec3::ZERO,
                        age: 0.0,
                        lifetime: len / FLOW_SPEED * (0.85 + 0.3 * self.rng.next_f64()),
                        color: INTAKE_COLOR,
                        base_scale: 0.016,
                        grow: 0.0,
                        kind: VfxKind::Flow {
                            anchor: VfxAnchor::Group(gi),
                            swirl: FLOW_SWIRL,
                            color_end,
                            pts,
                        },
                    });
                }
                // Cap electron on this end: captured + dispersed inflow.
                if blocker.is_some() && occ_v > 0.5 {
                    let n = spawn_count(&mut self.rng, end_rate * 0.5);
                    for _ in 0..n {
                        if self.vfx_particles.len() >= MAX_POOL {
                            break;
                        }
                        let (pts, len) = flow_capture(&mut self.rng, tip_y, entry, 1.0);
                        self.vfx_particles.push(VfxParticle {
                            position: com + g_rot * pts[0],
                            velocity: DVec3::ZERO,
                            age: 0.0,
                            lifetime: len / FLOW_SPEED * (0.85 + 0.3 * self.rng.next_f64()),
                            color: INTAKE_COLOR,
                            base_scale: 0.016,
                            grow: 0.0,
                            kind: VfxKind::Flow {
                                anchor: VfxAnchor::Group(gi),
                                swirl: FLOW_SWIRL,
                                color_end: ELECTRON_COLOR,
                                pts,
                            },
                        });
                    }
                }
            }
        }

        // Bond stream bridges — a molecular bond IS two facing polar
        // streams meeting head-on (fourier.pdf/jup3.pdf): show charge
        // flowing from both poles into the collision plane at the
        // midpoint. This is what visually separates a BONDED neighbor
        // (backed out, streams meeting across the gap) from a FUSED one
        // (no gap, no streams — one body).
        const BRIDGE_PER_SEC: f64 = 70.0; // per bond end
        const BRIDGE_SPEED: f64 = 2.4;
        const BRIDGE_COLOR: (f32, f32, f32) = (0.80, 0.92, 1.0);
        for (i, j) in self.molecular_bonds() {
            let mid = (emitters[i].pos + emitters[j].pos) * 0.5;
            for &a in &[i, j] {
                let em = &emitters[a];
                let gap = (mid - em.pos).length();
                let travel = gap - em.radius;
                if travel <= 0.05 {
                    continue;
                }
                let axis = (mid - em.pos) / gap;
                let (r1, r2) = build_frame(axis);
                let n_spawn = spawn_count(&mut self.rng, BRIDGE_PER_SEC);
                for _ in 0..n_spawn {
                    if self.vfx_particles.len() >= MAX_POOL {
                        break;
                    }
                    let ang = self.rng.next_f64() * TAU;
                    let lat_off = self.rng.next_f64() * 0.12;
                    let start = em.pos
                        + axis * em.radius
                        + (r1 * ang.cos() + r2 * ang.sin()) * lat_off;
                    self.vfx_particles.push(VfxParticle {
                        position: start,
                        velocity: axis * BRIDGE_SPEED,
                        age: 0.0,
                        // Dies AT the collision plane — the cushion.
                        lifetime: travel / BRIDGE_SPEED,
                        color: BRIDGE_COLOR,
                        base_scale: 0.022 * em.radius.max(0.35) as f32,
                        grow: 0.5, // splashes wider as it nears the plane
                        kind: VfxKind::Smoke,
                    });
                }
            }
        }

        // Build MultiMesh buffer: 12 (transform) + 4 (color) per particle
        let n = self.vfx_particles.len();
        let mut buf = Vec::with_capacity(n * 16);

        for vp in &self.vfx_particles {
            let t = (vp.age / vp.lifetime) as f32;
            // Quick fade-in; stays bright most of the trip and drops off
            // near the end — the parcel visibly REACHES its exit before
            // dying there.
            let fade = (t / 0.12).min(1.0) * (1.0 - t * t * t).max(0.0);
            let scale = (vp.base_scale * (1.0 + vp.grow * vp.age as f32)).max(0.004);
            // Flow parcels recolor as they recycle: intake tint on the way
            // in, exit tint on the way out.
            let (cr, cg, cb, brightness) = match vp.kind {
                VfxKind::Flow { color_end, .. } => (
                    vp.color.0 + (color_end.0 - vp.color.0) * t,
                    vp.color.1 + (color_end.1 - vp.color.1) * t,
                    vp.color.2 + (color_end.2 - vp.color.2) * t,
                    0.9f32, // tiny dots need the extra alpha
                ),
                VfxKind::Smoke => (vp.color.0, vp.color.1, vp.color.2, 0.6),
            };
            buf.extend_from_slice(&[
                scale, 0.0, 0.0, vp.position.x as f32,
                0.0, scale, 0.0, vp.position.y as f32,
                0.0, 0.0, scale, vp.position.z as f32,
                cr, cg, cb, fade * brightness,
            ]);
        }

        buf
    }

    pub fn vfx_count(&self) -> usize {
        self.vfx_particles.len()
    }

    pub fn set_vfx_enabled(&mut self, enabled: bool) {
        self.vfx_enabled = enabled;
        if !enabled {
            self.vfx_particles.clear();
        }
    }
}

// ── Helpers ───────────────────────────────────────────────────────────────

/// Quaternion rotating +Y onto `target` (falls back to +Y for zero vectors).
pub fn orientation_from_pole(target: DVec3) -> DQuat {
    let target = if target.length_squared() > 1e-9 {
        target.normalize()
    } else {
        DVec3::Y
    };
    if (target - DVec3::Y).length_squared() < 1e-9 {
        DQuat::IDENTITY
    } else if (target + DVec3::Y).length_squared() < 1e-9 {
        DQuat::from_rotation_z(std::f64::consts::PI)
    } else {
        DQuat::from_rotation_arc(DVec3::Y, target)
    }
}

pub fn build_frame(pole: DVec3) -> (DVec3, DVec3) {
    let right = if pole.dot(DVec3::X).abs() < 0.9 {
        pole.cross(DVec3::X).normalize()
    } else {
        pole.cross(DVec3::Z).normalize()
    };
    let forward = right.cross(pole);
    (right, forward)
}

// ── Charge-flow path builders (VFX) ──────────────────────────────────────
// All paths are polylines in an anchor's LOCAL frame (pole = +Y),
// resampled to uniform arc length so parcels travel at constant speed.

/// Resample a polyline to `n` points spaced uniformly by arc length.
/// Returns (points, total length).
fn resample_polyline_uniform(raw: &[DVec3], n: usize) -> (Vec<DVec3>, f64) {
    let n = n.max(2);
    if raw.len() < 2 {
        let p = raw.first().copied().unwrap_or(DVec3::ZERO);
        return (vec![p; n], 0.0);
    }
    let mut cum = Vec::with_capacity(raw.len());
    let mut total = 0.0;
    cum.push(0.0);
    for w in raw.windows(2) {
        total += (w[1] - w[0]).length();
        cum.push(total);
    }
    if total < 1e-9 {
        return (vec![raw[0]; n], 0.0);
    }
    let mut out = Vec::with_capacity(n);
    let mut seg = 0usize;
    for k in 0..n {
        let d = total * k as f64 / (n - 1) as f64;
        while seg + 2 < cum.len() && cum[seg + 1] < d {
            seg += 1;
        }
        let span = (cum[seg + 1] - cum[seg]).max(1e-12);
        let f = ((d - cum[seg]) / span).clamp(0.0, 1.0);
        out.push(raw[seg].lerp(raw[seg + 1], f));
    }
    (out, total)
}

/// Complete free-particle recycling loop: polar funnel in at `entry`
/// (+1 north / −1 south), through the body, out on the FAR hemisphere at
/// `theta_out` from the exit pole (sampled from the emission CDF — protons
/// leave via the disc, neutrons via the far pole), ending at `reach`.
fn flow_free_loop(
    rng: &mut Rng,
    radius: f64,
    entry: f64,
    theta_in: f64,
    theta_out: f64,
    reach: f64,
    spin: f64,
) -> (Vec<DVec3>, f64) {
    let mut raw: Vec<DVec3> = Vec::with_capacity(12);
    let phi0 = rng.next_f64() * TAU;
    let top = 2.6 * radius;
    let mouth = (top * theta_in.sin()).max(0.1 * radius);
    for k in 0..4 {
        let f = k as f64 / 3.0;
        let h = top - (top - radius) * f;
        let taper = ((h / radius - 1.0) / 1.6).clamp(0.0, 1.0);
        let rho = mouth * taper.powf(0.65) + 0.05 * radius;
        let phi = phi0 + spin * 2.2 * f;
        raw.push(DVec3::new(rho * phi.cos(), entry * h, rho * phi.sin()));
    }
    raw.push(DVec3::new(0.0, entry * 0.45 * radius, 0.0));
    raw.push(DVec3::new(0.0, -entry * 0.55 * radius, 0.0));
    let phi_e = phi0 + spin * 2.9;
    let (st, ct) = (theta_out.sin(), theta_out.cos());
    for j in 0..4 {
        let f = j as f64 / 3.0;
        let rr = radius * 0.95 + (reach - radius * 0.95) * f;
        let phi = phi_e + spin * 1.2 * f;
        raw.push(DVec3::new(
            st * phi.cos() * rr,
            -entry * ct * rr,
            st * phi.sin() * rr,
        ));
    }
    resample_polyline_uniform(&raw, 16)
}

/// Inflow toward an electron-capped pole (`tip_y` = signed local y of the
/// pole tip): down the funnel, intercepted at the rider's altitude, then
/// dispersed sideways off its little disc.
fn flow_capture(rng: &mut Rng, tip_y: f64, entry: f64, radius: f64) -> (Vec<DVec3>, f64) {
    let mut raw: Vec<DVec3> = Vec::with_capacity(7);
    let phi0 = rng.next_f64() * TAU;
    let mouth = radius * (0.5 + 0.8 * rng.next_f64());
    for k in 0..3 {
        let f = k as f64 / 2.0;
        let h = tip_y + entry * radius * (1.9 - 1.5 * f);
        let rho = mouth * (1.0 - f) + 0.10 * radius;
        let phi = phi0 + 1.8 * f;
        raw.push(DVec3::new(rho * phi.cos(), h, rho * phi.sin()));
    }
    let phi_e = phi0 + 2.4;
    for j in 0..3 {
        let f = j as f64 / 2.0;
        let rho = radius * (0.25 + 1.3 * f);
        let h = tip_y + entry * radius * (0.35 - 0.15 * f);
        let phi = phi_e + 1.4 * f;
        raw.push(DVec3::new(rho * phi.cos(), h, rho * phi.sin()));
    }
    resample_polyline_uniform(&raw, 12)
}

/// Through-charge in a nucleus: open end (`tip_y`) → the whole axial
/// channel → out the far pole to the composite reach — the stack's
/// pole-to-pole channel.
fn flow_group_through(
    rng: &mut Rng,
    entry: f64,
    tip_y: f64,
    far_tip_y: f64,
    reach_out: f64,
) -> (Vec<DVec3>, f64) {
    let mut raw: Vec<DVec3> = Vec::with_capacity(8);
    let phi0 = rng.next_f64() * TAU;
    let mouth = 0.6 + 0.9 * rng.next_f64();
    for k in 0..3 {
        let f = k as f64 / 2.0;
        let h = tip_y + entry * (1.9 - 1.7 * f);
        let rho = mouth * (1.0 - f) + 0.08;
        let phi = phi0 + 2.0 * f;
        raw.push(DVec3::new(rho * phi.cos(), h, rho * phi.sin()));
    }
    raw.push(DVec3::new(0.0, (tip_y + far_tip_y) * 0.5, 0.0));
    raw.push(DVec3::new(0.0, far_tip_y, 0.0));
    let out_r = reach_out.max(far_tip_y.abs() + 1.2);
    let phi_e = phi0 + 0.8;
    raw.push(DVec3::new(
        0.12 * phi_e.cos(),
        far_tip_y - entry * 0.6,
        0.12 * phi_e.sin(),
    ));
    raw.push(DVec3::new(
        0.18 * (phi_e + 0.7).cos(),
        -entry * out_r,
        0.18 * (phi_e + 0.7).sin(),
    ));
    resample_polyline_uniform(&raw, 16)
}

/// Disc exit in a nucleus: open end → axial channel down to one alpha's
/// latitude → flung out to that alpha's max-emission ring on the skin.
fn flow_group_disc(
    rng: &mut Rng,
    entry: f64,
    tip_y: f64,
    ring_r: f64,
    ring_y: f64,
) -> (Vec<DVec3>, f64) {
    let mut raw: Vec<DVec3> = Vec::with_capacity(9);
    let phi0 = rng.next_f64() * TAU;
    let mouth = 0.6 + 0.9 * rng.next_f64();
    for k in 0..3 {
        let f = k as f64 / 2.0;
        let h = tip_y + entry * (1.9 - 1.7 * f);
        let rho = mouth * (1.0 - f) + 0.08;
        let phi = phi0 + 2.0 * f;
        raw.push(DVec3::new(rho * phi.cos(), h, rho * phi.sin()));
    }
    raw.push(DVec3::new(0.0, (tip_y + ring_y) * 0.5, 0.0));
    raw.push(DVec3::new(0.0, ring_y, 0.0));
    let yj = ring_y + (rng.next_f64() - 0.5) * 0.3;
    let phi_e = phi0 + 2.6;
    for j in 0..3 {
        let f = (j as f64 + 1.0) / 3.0;
        let rho = 0.4 + (ring_r - 0.4).max(0.3) * f;
        let phi = phi_e + 1.3 * f;
        raw.push(DVec3::new(
            rho * phi.cos(),
            ring_y + (yj - ring_y) * f,
            rho * phi.sin(),
        ));
    }
    resample_polyline_uniform(&raw, 16)
}

pub fn default_spin_rate(name: &str) -> f64 {
    match name {
        "proton" => TAU * 3.0,
        "neutron" => TAU * 3.0,
        "electron" => TAU * 5.0,
        _ => TAU,
    }
}

pub fn profile_color(name: &str) -> (f32, f32, f32) {
    match name {
        "proton" => (0.92, 0.30, 0.20),
        "neutron" => (0.35, 0.50, 0.92),
        "electron" => (0.20, 0.90, 0.35),
        _ => (0.7, 0.7, 0.7),
    }
}

fn push_transform_color(
    buf: &mut Vec<f32>,
    p: &SimParticle,
    scale: f32,
    color: (f32, f32, f32),
    alpha: f32,
) {
    let pos = p.position;
    let x = p.orientation * DVec3::X;
    let y = p.orientation * DVec3::Y;
    let z = p.orientation * DVec3::Z;

    buf.extend_from_slice(&[
        x.x as f32 * scale, x.y as f32 * scale, x.z as f32 * scale, pos.x as f32,
        y.x as f32 * scale, y.y as f32 * scale, y.z as f32 * scale, pos.y as f32,
        z.x as f32 * scale, z.y as f32 * scale, z.z as f32 * scale, pos.z as f32,
    ]);
    buf.extend_from_slice(&[color.0, color.1, color.2, alpha]);
}

// ── Trace-model reorientation ─────────────────────────────────────────────

/// Reorient a saved spin-mode path trace so the pattern's symmetry axis —
/// the axial hole the body spins around, which in spin mode is the
/// outermost orbital level's LAB axis (X, Y or Z depending on the stack) —
/// lands on +Y, the pole axis atom-mode bodies spin about
/// (`member_spin_phase` / axial spin rotate meshes about local +Y).
///
/// The axis is recovered from the data itself, so traces saved before this
/// fix reorient too: for a surface of revolution the covariance eigenvalue
/// along the symmetry axis is distinct from the two (equal) transverse
/// ones — take the eigenvector whose eigenvalue is farthest from the other
/// two. Works for both oblate (disc) and prolate (spindle) patterns.
pub fn reorient_trace_to_pole(points: &[DVec3]) -> Vec<DVec3> {
    if points.len() < 8 {
        return points.to_vec();
    }
    let n = points.len() as f64;
    let centroid: DVec3 = points.iter().copied().sum::<DVec3>() / n;
    let mut cov = [[0.0f64; 3]; 3];
    for p in points {
        let d = *p - centroid;
        let v = [d.x, d.y, d.z];
        for (r, row) in cov.iter_mut().enumerate() {
            for (c, cell) in row.iter_mut().enumerate() {
                *cell += v[r] * v[c];
            }
        }
    }
    let (eigvals, eigvecs) = jacobi_eigen_3x3(cov);
    // Symmetry axis = the eigenvalue farthest from both others.
    let mut best = 0usize;
    let mut best_gap = f64::MIN;
    for k in 0..3 {
        let gap = (0..3)
            .filter(|&j| j != k)
            .map(|j| (eigvals[k] - eigvals[j]).abs())
            .fold(f64::MAX, f64::min);
        if gap > best_gap {
            best_gap = gap;
            best = k;
        }
    }
    let mut axis = eigvecs[best];
    if axis.length_squared() < 1e-12 {
        return points.to_vec();
    }
    axis = axis.normalize();
    if axis.dot(DVec3::Y) < 0.0 {
        axis = -axis; // hemisphere choice is free (body of revolution)
    }
    let rot = DQuat::from_rotation_arc(axis, DVec3::Y);
    points.iter().map(|p| rot * *p).collect()
}

/// Eigen-decomposition of a symmetric 3×3 matrix via cyclic Jacobi
/// rotations (Numerical Recipes convention). Returns (eigenvalues,
/// eigenvectors) with `eigenvectors[k]` paired to `eigenvalues[k]`.
fn jacobi_eigen_3x3(mut a: [[f64; 3]; 3]) -> ([f64; 3], [DVec3; 3]) {
    fn mat_mul(l: &[[f64; 3]; 3], r: &[[f64; 3]; 3]) -> [[f64; 3]; 3] {
        let mut out = [[0.0; 3]; 3];
        for i in 0..3 {
            for j in 0..3 {
                for (k, lk) in l[i].iter().enumerate() {
                    out[i][j] += lk * r[k][j];
                }
            }
        }
        out
    }
    let mut v = [[0.0f64; 3]; 3];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }
    for _ in 0..64 {
        // Largest off-diagonal element.
        let (mut p, mut q, mut max) = (0usize, 1usize, 0.0f64);
        for r in 0..3 {
            for c in (r + 1)..3 {
                if a[r][c].abs() > max {
                    max = a[r][c].abs();
                    p = r;
                    q = c;
                }
            }
        }
        if max < 1e-13 {
            break;
        }
        let tau = (a[q][q] - a[p][p]) / (2.0 * a[p][q]);
        let t = tau.signum() / (tau.abs() + (1.0 + tau * tau).sqrt());
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;
        let mut j = [[0.0f64; 3]; 3];
        for (i, row) in j.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        j[p][p] = c;
        j[q][q] = c;
        j[p][q] = s;
        j[q][p] = -s;
        let jt = [
            [j[0][0], j[1][0], j[2][0]],
            [j[0][1], j[1][1], j[2][1]],
            [j[0][2], j[1][2], j[2][2]],
        ];
        a = mat_mul(&mat_mul(&jt, &a), &j);
        v = mat_mul(&v, &j);
    }
    (
        [a[0][0], a[1][1], a[2][2]],
        [
            DVec3::new(v[0][0], v[1][0], v[2][0]),
            DVec3::new(v[0][1], v[1][1], v[2][1]),
            DVec3::new(v[0][2], v[1][2], v[2][2]),
        ],
    )
}

// ── Headless CSV loading (tests + atom_lab bin) ──────────────────────────

/// Load a `degree,percentage` histogram CSV, returning the value column.
/// Panics with a clear message on malformed files (test/dev tool, not UI path).
pub fn load_histogram_csv(path: &Path) -> Vec<f32> {
    let text = std::fs::read_to_string(path)
        .unwrap_or_else(|e| panic!("cannot read {}: {}", path.display(), e));
    let mut values = Vec::with_capacity(181);
    for line in text.lines().skip(1) {
        let line = line.trim();
        if line.is_empty() {
            continue;
        }
        let value = line
            .split(',')
            .nth(1)
            .unwrap_or_else(|| panic!("malformed CSV line in {}: {}", path.display(), line));
        values.push(value.trim().parse::<f32>().unwrap_or_else(|e| {
            panic!("bad value in {}: {} ({})", path.display(), line, e)
        }));
    }
    values
}

/// The godot/config directory, resolved relative to the Cargo manifest so
/// tests and bins work regardless of the working directory.
pub fn config_dir() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../godot/config")
}

// ── Tests ─────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn make_flat_csv() -> Vec<f32> {
        vec![1.0; 181]
    }

    fn make_proton_like_csv() -> Vec<f32> {
        let mut csv = vec![0.0; 181];
        // Equator (index 90) = bright, poles (0 and 180) = dim
        for i in 0..181 {
            let deg = i as f64 - 90.0;
            let frac = 1.0 - (deg.abs() / 90.0);
            csv[i] = frac as f32;
        }
        csv
    }

    #[test]
    fn emission_table_bilateral_symmetry() {
        let csv = make_proton_like_csv();
        let table = EmissionTable::from_csv_values(&csv);
        let pole_val = table.sample(1.0);
        let anti_pole = table.sample(-1.0);
        assert!(
            (pole_val - anti_pole).abs() < 1e-6,
            "Bilateral: pole={} vs anti-pole={}",
            pole_val,
            anti_pole
        );
    }

    #[test]
    fn emission_table_equator_vs_pole() {
        let csv = make_proton_like_csv();
        let table = EmissionTable::from_csv_values(&csv);
        let equator = table.sample(0.0);
        let pole = table.sample(1.0);
        assert!(
            equator > pole,
            "Equator should be brighter: eq={} pole={}",
            equator,
            pole
        );
    }

    #[test]
    fn complement_inverse() {
        let csv = make_proton_like_csv();
        let emission = EmissionTable::from_csv_values(&csv);
        let absorption = emission.complement();
        let e = emission.sample(0.5);
        let a = absorption.sample(0.5);
        assert!(
            (e + a - 1.0).abs() < 1e-6,
            "emission + absorption should be 1.0: {} + {} = {}",
            e,
            a,
            e + a
        );
    }

    #[test]
    fn flat_emission_gives_uniform_one() {
        let table = EmissionTable::from_csv_values(&make_flat_csv());
        for i in 0..=10 {
            let cos_t = (i as f64 - 5.0) / 5.0;
            let val = table.sample(cos_t);
            assert!(
                (val - 1.0).abs() < 1e-5,
                "Flat CSV should give 1.0 everywhere, got {} at cos_theta={}",
                val,
                cos_t
            );
        }
    }

    #[test]
    fn histogram_csvs_load_and_register() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        for name in ["proton", "neutron", "electron"] {
            let csv = load_histogram_csv(&dir.join(format!("histogram_{name}.csv")));
            assert!(csv.len() >= 180, "{name} CSV too short: {}", csv.len());
            core.register_profile(name, 1.0, 1.0, &csv);
        }
        assert_eq!(core.profiles.len(), 3);
        // Proton emission must peak at the equator, vanish near the pole
        let proton = &core.profiles[0];
        assert!(proton.emission.sample(0.0) > 0.5);
        assert!(proton.emission.sample(1.0) < 0.1);
    }

    /// Oblate (disc) trace wound around the X axis: the hole must land on Y.
    #[test]
    fn trace_reorients_disc_hole_to_y() {
        let mut pts = Vec::new();
        for i in 0..400 {
            let a = i as f64 * 0.377;
            let jig = (((i * 7) % 13) as f64 / 13.0 - 0.5) * 0.3; // thin along X
            pts.push(DVec3::new(jig, a.cos(), a.sin()));
        }
        let out = reorient_trace_to_pole(&pts);
        let max_y = out.iter().map(|p| p.y.abs()).fold(0.0, f64::max);
        assert!(max_y < 0.2, "hole axis should land on Y, max|y|={max_y}");
        let max_ring = out
            .iter()
            .map(|p| (p.x * p.x + p.z * p.z).sqrt())
            .fold(0.0, f64::max);
        assert!(
            (max_ring - 1.0).abs() < 0.05,
            "ring radius should be preserved: {max_ring}"
        );
    }

    /// Prolate (spindle) trace along Z: the long axis must land on Y.
    #[test]
    fn trace_reorients_spindle_to_y() {
        let mut pts = Vec::new();
        for i in 0..400 {
            let a = i as f64 * 0.377;
            let h = ((i % 41) as f64 / 40.0 - 0.5) * 4.0; // long along Z
            pts.push(DVec3::new(0.5 * a.cos(), 0.5 * a.sin(), h));
        }
        let out = reorient_trace_to_pole(&pts);
        let max_y = out.iter().map(|p| p.y.abs()).fold(0.0, f64::max);
        assert!(max_y > 1.8, "long axis should land on Y, max|y|={max_y}");
        let max_lat = out
            .iter()
            .map(|p| (p.x * p.x + p.z * p.z).sqrt())
            .fold(0.0, f64::max);
        assert!(max_lat < 0.6, "transverse spread should stay small: {max_lat}");
    }

    /// A trace already aligned to Y must come back (near-)unchanged.
    #[test]
    fn trace_reorient_is_stable_for_canonical_input() {
        let mut pts = Vec::new();
        for i in 0..400 {
            let a = i as f64 * 0.377;
            let jig = (((i * 7) % 13) as f64 / 13.0 - 0.5) * 0.3;
            pts.push(DVec3::new(a.cos(), jig, a.sin()));
        }
        let out = reorient_trace_to_pole(&pts);
        // Data-driven axis recovery on jittered data is exact only to
        // ~milliradians — "stable" means no visible rotation, not bitwise.
        for (a, b) in pts.iter().zip(&out) {
            assert!(
                (*a - *b).length() < 0.01,
                "canonical trace moved: {a:?} -> {b:?}"
            );
        }
    }

    /// The wall-riding electron geometry (r=1.3, θ≈11°) must fully occlude
    /// the pole it rides and leave the far pole open.
    #[test]
    fn pole_occlusion_stoppers_ridden_pole_only() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let e_csv = load_histogram_csv(&dir.join("histogram_electron.csv"));
        let p_id = core.register_profile("proton", 1.0, 1.0, &p_csv);
        let e_id = core.register_profile("electron", 1.0 / 1836.0, 0.3, &e_csv);
        core.spawn_particle(p_id, DVec3::ZERO, DVec3::ZERO, DVec3::Y);
        // Rider over the +Y pole at the session-28 orbit geometry.
        core.spawn_particle(e_id, DVec3::new(0.248, 1.276, 0.0), DVec3::ZERO, DVec3::Y);
        let occ = core.pole_occlusion();
        assert!(
            occ[0][0] > 0.9,
            "north pole should be stoppered by the rider: {:?}",
            occ[0]
        );
        assert!(
            occ[0][1] < 1e-9,
            "south pole should stay open: {:?}",
            occ[0]
        );
    }

    /// Electron-capped pole: the inflow is CAPTURED at the rider — no
    /// parcel entering the capped side passes through the body; the open
    /// pole still runs full recycling loops.
    #[test]
    fn flow_captured_at_electron_capped_pole() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let e_csv = load_histogram_csv(&dir.join("histogram_electron.csv"));
        let p_id = core.register_profile("proton", 1.0, 1.0, &p_csv);
        let e_id = core.register_profile("electron", 1.0 / 1836.0, 0.3, &e_csv);
        core.spawn_particle(p_id, DVec3::ZERO, DVec3::ZERO, DVec3::Y);
        core.spawn_particle(e_id, DVec3::new(0.248, 1.276, 0.0), DVec3::ZERO, DVec3::Y);
        for _ in 0..60 {
            core.advance_clouds(1.0 / 60.0);
        }
        let (mut north_in, mut north_through, mut south_in) = (0usize, 0usize, 0usize);
        for vp in &core.vfx_particles {
            if let VfxKind::Flow {
                anchor: VfxAnchor::Particle(0),
                ref pts,
                ..
            } = vp.kind
            {
                let first = pts[0];
                let last = pts[pts.len() - 1];
                if first.y > 0.0 {
                    north_in += 1;
                    if last.y < -0.5 {
                        north_through += 1;
                    }
                } else {
                    south_in += 1;
                }
            }
        }
        assert!(north_in > 0, "capped pole should show captured inflow");
        assert_eq!(north_through, 0, "no parcel may pass the stoppered channel");
        assert!(south_in > 0, "open pole should run full recycling loops");
    }

    /// Nucleus flow: parcels enter the open stack ends; some peel off at
    /// alpha disc rings, some pass through the whole stack and out the far
    /// pole.
    #[test]
    fn group_flow_disc_and_through_exits() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        for _ in 0..90 {
            core.advance_clouds(1.0 / 60.0);
        }
        let (mut disc, mut through) = (0usize, 0usize);
        for vp in &core.vfx_particles {
            if let VfxKind::Flow {
                anchor: VfxAnchor::Group(0),
                ref pts,
                ..
            } = vp.kind
            {
                let last = pts[pts.len() - 1];
                let lat = (last.x * last.x + last.z * last.z).sqrt();
                if lat > 2.0 {
                    disc += 1;
                } else if last.y.abs() > 4.6 {
                    through += 1;
                }
            }
        }
        assert!(disc > 0, "expected disc-ring exits, got none");
        assert!(through > 0, "expected through-channel exits, got none");
    }

    /// The charge-field lock: a spun-up nucleus relaxes back to rest
    /// instead of tumbling forever.
    #[test]
    fn group_tumble_relaxes() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        core.running = true;
        core.groups[0].angular_velocity = DVec3::new(3.0, 1.0, 2.0);
        core.step_n(2000); // 1.0 sim s ⇒ e^{−20} decay
        let w = core.groups[0].angular_velocity.length();
        assert!(w < 0.05, "nucleus should relax to rest, |ω|={w}");
    }

    /// One max-emission ring per alpha block: alpha → 1, carbon → 3.
    #[test]
    fn emission_rings_one_per_alpha() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        let rings = core.build_group_emission_rings(0);
        assert!(!rings.is_empty(), "alpha should have an emission ring");
        assert_eq!(rings[0] as usize, 1, "alpha = one alpha block, one ring");

        core.clear_particles();
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        let rings = core.build_group_emission_rings(0);
        assert_eq!(rings[0] as usize, 3, "carbon = three alphas, three rings");
    }

    /// Composite alpha skin: a sane lathe whose equatorial reach exceeds a
    /// single free proton's (two fused discs push farther, ≈ ×2^¼).
    #[test]
    fn group_skin_reach_sane_for_alpha() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        let p_id = core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        let (lon, lat) = (16usize, 8usize);
        let buf = core.build_group_skin_mesh(0, lon, lat);
        let rings = lat * 2 + 1;
        let vpr = lon + 1;
        assert_eq!(buf[0] as usize, rings * vpr);
        let ring_r = |ring: usize| -> f64 {
            let o = 2 + ring * vpr * 10;
            DVec3::new(buf[o] as f64, buf[o + 1] as f64, buf[o + 2] as f64).length()
        };
        let eq = ring_r(rings / 2);
        let single =
            (core.couplings.c_q * core.profiles[p_id].emission.sample(0.0)).powf(0.25);
        assert!(
            eq > single * 1.05,
            "alpha equator reach {eq} should exceed a single proton's {single}"
        );
        assert!(eq < 12.0, "alpha equator reach implausibly large: {eq}");
        for ring in 0..rings {
            let r = ring_r(ring);
            assert!(r.is_finite() && r < 40.0, "ring {ring} reach bad: {r}");
        }
    }

    /// Fused constituents draw no free-field skin: an alpha (2p+2n) plus
    /// one free proton must produce a 3-instance body buffer but a
    /// 1-instance skin buffer for the proton profile.
    #[test]
    fn fused_constituents_have_no_skin_instances() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        let p_id = core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        core.spawn_particle(p_id, DVec3::new(6.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
        let bodies = core.build_multimesh_buffer_for_profile(p_id).len() / 16;
        let skins = core.build_skin_multimesh_buffer_for_profile(p_id).len() / 16;
        assert_eq!(bodies, 3, "2 fused + 1 free proton bodies");
        assert_eq!(skins, 1, "only the free proton gets a skin");
    }

    /// Skin radius carries the reach law r(θ) = (C_q·m·E(θ))^¼, clamped to
    /// hug the body at the polar holes.
    #[test]
    fn field_skin_reach_matches_quarter_power_law() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let p_id = core.register_profile("proton", 1.0, 1.0, &p_csv);
        let (lon, lat) = (16usize, 8usize);
        let buf = core.build_field_skin_mesh(p_id, lon, lat);
        let rings = lat * 2 + 1;
        let verts_per_ring = lon + 1;
        assert_eq!(buf[0] as usize, rings * verts_per_ring);

        let vert_pos = |ring: usize| -> DVec3 {
            let o = 2 + ring * verts_per_ring * 10;
            DVec3::new(buf[o] as f64, buf[o + 1] as f64, buf[o + 2] as f64)
        };
        // Equator (middle ring): reach = (C_q · E_eq)^¼ in particle radii.
        let e_eq = core.profiles[p_id].emission.sample(0.0);
        let expected = (core.couplings.c_q * e_eq).powf(0.25);
        let rho_eq = vert_pos(rings / 2).length();
        assert!(
            (rho_eq - expected).abs() < 0.02 * expected,
            "equator reach {rho_eq} != (C_q·E)^¼ = {expected}"
        );
        // Pole: emission ≈ 0 ⇒ skin hugs the body.
        let rho_pole = vert_pos(0).length();
        assert!(
            rho_pole < 1.5,
            "polar skin should pinch to the body, got {rho_pole}"
        );
        assert!(
            rho_eq > 3.0 * rho_pole,
            "skin must read as a wide equatorial ledge: eq={rho_eq} pole={rho_pole}"
        );
    }
}
