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
    /// Index into `group`'s `RigidGroup::alphas` when this particle is a
    /// group constituent; None for free particles. Set at spawn
    /// (`spawn_preset`) — the Part 2 mode-aware force skip predicate needs
    /// it (session-31 addendum A9: RigidAlpha skips same-ALPHA pairs, not
    /// same-group pairs).
    pub alpha: Option<usize>,
    /// Render frame for FREE particles: its Y axis smoothly tracks the
    /// physics pole (incremental parallel transport — no swing–twist
    /// extraction, which is singular near 180° configurations and made
    /// bodies flip wildly, session-29), while spinning about the pole at
    /// the readable display rate. Physics orientation is untouched.
    pub display_orientation: DQuat,
    /// Live through-charge flow state (session 32, Phase A —
    /// docs/THROUGH_CHARGE_DESIGN.md). Filled by `solve_charge_flow`
    /// during stepping; no force term reads it yet.
    pub flow: crate::charge_flow::FlowState,
}

impl SimParticle {
    pub fn pole_axis(&self) -> DVec3 {
        self.orientation * DVec3::Y
    }
}

// ── Nuclear presets (rigid composites) ───────────────────────────────────

/// One member of an [`AlphaSpec`], in the ALPHA's own rest frame
/// (+Y = the alpha's stack axis — same convention as a free particle's
/// pole). `spin_sign` is the intrinsic physical axial-spin sign.
#[derive(Clone)]
struct AlphaMemberSpec {
    profile_name: &'static str,
    local_pos: DVec3,
    local_pole: DVec3,
    spin_sign: f64,
}

/// Builder-time description of one alpha unit: its rest frame in the
/// NUCLEUS (group) frame, whether it rides the carousel, and its members
/// in the alpha's OWN frame. `preset_alphas` assembles nuclei from these;
/// `spawn_preset` turns them into runtime [`AlphaUnit`]s.
#[derive(Clone)]
struct AlphaSpec {
    rest_axis: DVec3,
    rest_center: DVec3,
    orbits_core: bool,
    members: Vec<AlphaMemberSpec>,
}

/// Radius off-axis for an alpha's neutron posts (or a sideways alpha's
/// posts) — how far the posts sit from the alpha's own stack axis, in the
/// alpha's own rest frame. Replaces the four `0.7` literals that used to
/// appear separately in every builder.
///
/// Session 32: briefly widened to 1.05 because at 0.7 the two posts sat
/// 1.4 apart with DISC-shaped contact radius 1.0 each — interpenetrating
/// at rest, so the moment FreeNucleon activated intra-alpha forces the
/// compressed contact spring (58 force units) fired both posts out
/// instantly (the real reason `nucleon_balance` read ~2500% drift at
/// every boost). But widening the posts degraded the RigidAlpha flyby
/// envelope (1/11 robust — inter-alpha geometry left its earned basin),
/// so the fix moved to the CONTACT SHAPE instead: same-group neutrons
/// collide as thin RODS (phos.pdf "neutrons are 1D... lightning rods"),
/// not discs — see the contact section in `compute_forces`. Geometry
/// stays at the earned 0.7.
const POST_R: f64 = 0.7;

/// Pole anatomy of the alpha's INTERNAL neutron posts (session-32
/// experiment, user design question from end of session 31).
///
/// - `Axial` (the session-31 modeling choice, NOT pinned by the papers):
///   post pole parallel to the alpha stack axis.
/// - `Radial`: post pole along its own lateral offset — sideways to the
///   through-stream, the pole facing the region where the two proton
///   discs' emission crosses. The user reads deut.pdf this way: posts are
///   CHARGE CHANNELS ("not only acting as posts... also charge channels";
///   neutrons 1D "lightning rods", protons 2D fans) acting as
///   self-balancing regulators — a disc dip feeds more flow into that
///   side's neutron pole, which flings it back out toward the disc as a
///   restoring nudge.
///
/// Does NOT affect the polar plug neutrons (`plug_neutron`): those are
/// pole-on-axis per graphene.pdf — a DIFFERENT, paper-pinned position.
/// Selected via `AtomCore::post_anatomy`; takes effect at `spawn_preset`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PostAnatomy {
    #[default]
    Axial,
    Radial,
}

/// The standard 4-nucleon alpha shape, in the ALPHA's OWN frame: two
/// protons stacked hole-to-hole on the axis (CDs spinning the SAME
/// direction so charge channels through pole-to-pole as a dipole —
/// milesmathis.com/oxygen.pdf), with the two neutrons as posts straddling
/// the midsection, off-axis, keeping the protons from turning (oxygen.pdf).
/// Posts sit off-axis so the axial charge channel stays open (nuclear.pdf:
/// the hole in the CD is the recycling channel) — their exact starting
/// azimuth doesn't matter (each alpha's `roll_phase` is a free, randomized
/// DOF, session-31): what matters is they're off-axis and roll with the
/// alpha. Post POLE orientation is the open [`PostAnatomy`] question.
/// Shared by every alpha builder below (core, cap, connector, carousel) —
/// only the alpha's `rest_axis`/`rest_center`/`orbits_core` differ.
fn alpha_members(anatomy: PostAnatomy) -> Vec<AlphaMemberSpec> {
    let post_pole = |lateral: DVec3| match anatomy {
        PostAnatomy::Axial => DVec3::Y,
        PostAnatomy::Radial => lateral.normalize(),
    };
    vec![
        AlphaMemberSpec {
            profile_name: "proton",
            local_pos: DVec3::new(0.0, -NUCLEON_PITCH / 2.0, 0.0),
            local_pole: DVec3::Y,
            spin_sign: 1.0,
        },
        AlphaMemberSpec {
            profile_name: "proton",
            local_pos: DVec3::new(0.0, NUCLEON_PITCH / 2.0, 0.0),
            local_pole: DVec3::Y,
            spin_sign: 1.0,
        },
        AlphaMemberSpec {
            profile_name: "neutron",
            local_pos: DVec3::new(-POST_R, 0.0, 0.0),
            local_pole: post_pole(DVec3::new(-POST_R, 0.0, 0.0)),
            spin_sign: 1.0,
        },
        AlphaMemberSpec {
            profile_name: "neutron",
            local_pos: DVec3::new(POST_R, 0.0, 0.0),
            local_pole: post_pole(DVec3::new(POST_R, 0.0, 0.0)),
            spin_sign: 1.0,
        },
    ]
}

/// One alpha block centered at `y` on the stack axis: parallel to the
/// core (used for the core alpha itself and, in argon, the top/bottom
/// caps — nuclear.pdf: "Level one is the center disk... Level 3b is the
/// caps top and bottom"). Does not ride the carousel.
fn core_alpha(y: f64, anatomy: PostAnatomy) -> AlphaSpec {
    AlphaSpec {
        rest_axis: DVec3::Y,
        rest_center: DVec3::new(0.0, y, 0.0),
        orbits_core: false,
        members: alpha_members(anatomy),
    }
}

/// One alpha mounted ON THE CAROUSEL at azimuth `phi_deg`: its stack axis
/// points RADIALLY (the center disk's equatorial output feeds this alpha's
/// axial hole — nuclear.pdf's first carousel configuration, Neon), and the
/// whole block rides the carousel around the core. Charge is flung out
/// equatorially through these (dielec.pdf/diamag.pdf).
fn carousel_alpha(phi_deg: f64, anatomy: PostAnatomy) -> AlphaSpec {
    let phi = phi_deg.to_radians();
    let u = DVec3::new(phi.cos(), 0.0, phi.sin()); // radial stack axis
    AlphaSpec {
        rest_axis: u,
        rest_center: u * CAROUSEL_R,
        orbits_core: true,
        members: alpha_members(anatomy),
    }
}

/// Connector alpha: mounted SIDEWAYS on the stack axis between the center
/// and a cap — facing like the polar plugs (edge-on, its disc feeding the
/// axial channel), riding the carousel rotation and spinning on its own
/// pole (nuclear.pdf Argon: "Level 3a is the posts up and down" — the
/// session-29 diagram reading has them turned 90° to the core).
fn sideways_alpha(y: f64, anatomy: PostAnatomy) -> AlphaSpec {
    AlphaSpec {
        rest_axis: DVec3::X,
        rest_center: DVec3::new(0.0, y, 0.0),
        orbits_core: true,
        members: alpha_members(anatomy),
    }
}

/// Polar plug proton: plugged into a stack pole with its pole PERPENDICULAR
/// to the stack axis, so its equatorial disc output feeds the stack's open
/// polar channel (phos.pdf plug-and-socket; ammon.pdf: N = C-stack + proton
/// in the south pole, neutron in the north; O = protons in both poles).
/// LONE plug only (nitrogen) — paired plugs are [`plug_pair`], one FUSED
/// alpha. Single-member alpha: its own rest axis IS the old world-frame
/// pole, so `Roll` spins the plug disc about its own pole and the carousel
/// `Car` rides it around the socket exactly as the old `carousel: true`
/// flag did (session-29 user note).
fn plug_proton(y: f64, z: f64) -> AlphaSpec {
    let rest_axis = if z.abs() > 1e-9 {
        DVec3::new(0.0, 0.0, -z.signum()) // toward the paired plug
    } else {
        DVec3::X
    };
    AlphaSpec {
        rest_axis,
        rest_center: DVec3::new(0.0, y, z),
        orbits_core: true,
        members: vec![AlphaMemberSpec {
            profile_name: "proton",
            local_pos: DVec3::ZERO,
            local_pole: DVec3::Y,
            spin_sign: 1.0,
        }],
    }
}

/// Polar plug PAIR: a proton+neutron sharing a polar hole as ONE FUSED
/// 2-member alpha (session 34). The session-32 model spawned the pair as
/// two independent 1-nucleon alphas, and the static well scan
/// (report_plug_well) found that configuration mechanically impossible:
/// at hole-sharing range (0.7–1.6) some equator always blasts some
/// partner at 1/r⁴ — with the pole-at-partner azimuth, the neutron's
/// REAL equatorial emission leak (~0.5) fired ~150 units of charge push
/// straight down the proton's pole-on absorption maximum, two orders of
/// magnitude above all available glue (invisible in the anatomy tables,
/// which only print the force on the second body of each pair). No
/// azimuth escapes this: rotating the proton tangential just redirects
/// the blast onto the neutron via the proton's full disc fan.
///
/// The papers resolve it as FUSION: a p+n in contact sharing a pole IS
/// the deuterium configuration (deut.pdf), and fused units "can't be
/// broken and rearranged" (uf4.pdf) — the same pre-fusion rule that
/// already exempts every core alpha's members from mutual forces in
/// RigidAlpha. graphene.pdf's carbon plugs "stayed in line" — they move
/// as one. Geometry is unchanged from the paired singles (proton at
/// −z with pole toward its partner, neutron at +z pole-on-axis); the
/// pair also gets honest mass-2 inertia, retiring the near-zero-inertia
/// torque overshoot of 1-nucleon plug alphas.
fn plug_pair(y: f64, gap: f64) -> AlphaSpec {
    AlphaSpec {
        rest_axis: DVec3::Y,
        rest_center: DVec3::new(0.0, y, 0.0),
        orbits_core: true,
        members: vec![
            AlphaMemberSpec {
                profile_name: "proton",
                local_pos: DVec3::new(0.0, 0.0, -gap),
                local_pole: DVec3::Z, // toward the paired neutron
                spin_sign: 1.0,
            },
            AlphaMemberSpec {
                profile_name: "neutron",
                local_pos: DVec3::new(0.0, 0.0, gap),
                local_pole: DVec3::Y, // pole-on-axis (graphene.pdf)
                spin_sign: 1.0,
            },
        ],
    }
}

/// Polar plug neutron: pole ON the stack axis — "the neutron is plugged in
/// with its pole pointing down … protons channel charge pole to equator,
/// while neutrons channel pole to pole" (graphene.pdf), so its axial pole
/// pulls charge into the stack's hole (atmo2.pdf: paired neutrons are
/// "pulling charge into the axial holes"). `z` offsets it off-axis to sit
/// side by side with a paired plug proton (0 for a lone plug). Single-
/// member alpha, rides the carousel like the alpha posts.
fn plug_neutron(y: f64, z: f64) -> AlphaSpec {
    AlphaSpec {
        rest_axis: DVec3::Y,
        rest_center: DVec3::new(0.0, y, z),
        orbits_core: true,
        members: vec![AlphaMemberSpec {
            profile_name: "neutron",
            local_pos: DVec3::ZERO,
            local_pole: DVec3::Y,
            spin_sign: 1.0,
        }],
    }
}

/// Half-gap between the members of a proton+neutron pair sharing a polar
/// hole ("two baryons in the hole fill the hole much better" — atmo2.pdf).
/// Same nestling scale as the alpha's neutron posts (±POST_R off-axis).
///
/// Session-33 finding (report_plug_retention anatomy): the pair rest
/// distance (1.6) is NOT a force equilibrium — net inward pull ≈ 1.4
/// (nuclear-ambient shadow + gravity + intake vs a small stream cushion;
/// charge repulsion is exactly zero pole-on) with nothing opposing it
/// until the disc-aware contact wall at r ≈ 0.7. The collapse-bounce on
/// the contact spring is a KE pump. Experiments narrowing the gap to
/// 0.35 (spawn at contact) and saturating the close-range ambient made
/// 300k retention WORSE, not better (see `AtomCore::plug_pair_gap` /
/// `ambient_sat_r` — runtime knobs for the factored matrix sweep,
/// report_plug_matrix), so the baseline stays until a cell of that
/// matrix earns a change.
const PLUG_PAIR_GAP: f64 = 0.8;

/// Axial seat of a FUSED plug pair beyond the end proton of its stack —
/// the MEASURED axial force equilibrium of the pair against the socket
/// (session 34, report_plug_well extended scan: net force crosses zero
/// ≈1.37 beyond the end proton, i.e. δ ≈ −1.23 from the old 2.6
/// funnel-mouth spawn; slope ≈ −89/unit, contact wall 0.3 further in).
/// The pair spawns 1.4 out — a residual ~2.5 inward settles it the last
/// ~0.03. The old funnel-mouth seat (2.6, inherited from stacked-alpha
/// pitch) left the pair 1.2 units up a monotone-inward slope: it FELL,
/// gained ~KE 10, and scattered the nucleus (the deep source of every
/// "plug retention" failure once the pair itself was fused). Paper
/// basis for the tight seat: ammon.pdf "that seventh proton is in
/// tight, so it channels with less loss"; the graphene.pdf diagrams
/// nest the plugs in the pole hole itself.
const PLUG_SEAT: f64 = 1.4;

/// Radius of the carousel level: distance from the stack axis to a
/// carousel alpha's center. Sets the nearest carousel proton pole just
/// outside the center disk's edge — plugged edge-to-hole (nuclear.pdf,
/// four.pdf: "all disks fit together edge to hole, like male and female
/// sockets"). Proportional rescale; exact funnel coupling of the carousel
/// level is a future refinement.
const CAROUSEL_R: f64 = 4.5;

/// Stacked-nucleon pitch: a fused neighbor parks at the MOUTH of the
/// source's intake funnel (2.6 r — see the funnel constants in the
/// VFX), its polar skin edge (1.06 r) hanging down inside the funnel
/// where the density still feeds momentum back (session-29 decision;
/// nuclear.pdf: fusion forces baryons inside the free-field standoff,
/// but not into contact).
pub(crate) const NUCLEON_PITCH: f64 = 2.6;

/// Channeling recapture tail multiplier (session-31 round 4 — binding v2
/// hardening). `channeling_factor`'s distance falloff used to hit exactly
/// zero at `2·NUCLEON_PITCH` (5.2) — a CLIFF, not a taper: an alpha
/// displaced (or shoved) even slightly past that radius instantly regained
/// full molecular `c_q`/`stream` repulsion with no attractive term to pull
/// it back, so every excursion past the cliff was a one-way ejection. This
/// is the root cause of both user-observed failures — the spontaneous late
/// "Jenga" collapse (one alpha random-walks past 5.2 during ordinary
/// thermal jitter and never returns) and the trivial destruction by a
/// nearby stray particle (any nudge past the cliff is terminal). Extending
/// the taper to `CHANNEL_TAIL · NUCLEON_PITCH` (7.8) gives a displaced
/// alpha a restoring basin: channeling — and therefore some of the binding
/// glue — is still partially in effect out past the old cutoff, so a
/// perturbed alpha is pulled back rather than ejected outright.
pub(crate) const CHANNEL_TAIL: f64 = 3.0;

/// Disc-aware contact fraction (A13, session-31 round 3): the bodies in
/// this model are DISCS, not spheres — the whole force model is planar
/// emission (`EmissionTable`, equator-bright/pole-dark). Facing protons of
/// adjacent alphas rest well inside the old sphere-contact wall
/// (r_i+r_j = 2.0, an unquenchable ~85-unit repulsion no channeling term
/// could touch), because a pole-on approach nests into the partner's
/// funnel and only the thin disc waist can actually collide, while an
/// equator-on approach presents the full radius. Effective per-side
/// contact radius:
/// `r_eff = r · (POLE_HALF_THICKNESS + (1 − POLE_HALF_THICKNESS)·sinθ)`,
/// θ = angle between the separation direction and that body's pole
/// (sinθ=0 pole-on ⇒ r_eff=POLE_HALF_THICKNESS·r; sinθ=1 equator-on ⇒
/// r_eff=r, unchanged from the old sphere-sum contact). 0.35 is the
/// waist-to-radius ratio the funnel-mouth geometry implies for a nested
/// pole-on pair (nuclear.pdf — fusion packs baryons past the free-field
/// standoff without true contact). SAME-GROUP-GATED in `compute_forces`
/// (not global): applying it to every pair, including molecular
/// proton–electron contact, broke the LOCKED M5 molecular suite
/// (`hydrogen_capture`, `derived_constants_equilibrium`,
/// `hydrogen_orbit_stable_long_run`, `alpha_captures_electron`) by
/// thinning the contact wall a near-pole-on electron hits, so it is
/// restricted to RigidAlpha/FreeNucleon intra-nucleus pairs where the
/// disc-nesting argument was actually made.
const POLE_HALF_THICKNESS: f64 = 0.35;

/// Alpha stack pitch: adjacent alpha centers along the axis. Keeps the
/// inter-block proton gap proportional to the old layout under the new
/// pitch.
const ALPHA_PITCH: f64 = 3.75;

/// Build a nuclear preset in alpha-unit terms: a list of [`AlphaSpec`]s
/// (rest frame + orbits_core + members), consumed by `spawn_preset`.
/// Preset compositions are unchanged from the flat-constituent era —
/// only the grouping into alpha units changed (session-31).
fn preset_alphas(
    name: &str,
    anatomy: PostAnatomy,
    plug_gap: f64,
) -> Option<Vec<AlphaSpec>> {
    // A symmetric stack of `n` core alphas about the origin, ALPHA_PITCH
    // apart — the axial spine the element presets build on.
    let stack = |n: usize| -> Vec<AlphaSpec> {
        (0..n)
            .map(|k| {
                let y = (k as f64 - (n as f64 - 1.0) / 2.0) * ALPHA_PITCH;
                core_alpha(y, anatomy)
            })
            .collect()
    };
    match name {
        "alpha" => Some(vec![core_alpha(0.0, anatomy)]),
        // Carbon (session-32 correction, user + meth.pdf/graphene.pdf):
        // TWO stacked alphas with a proton+neutron PLUG PAIR at each
        // pole — 6p 6n. graphene.pdf: "in the case of Carbon, we will
        // get a proton on each pole, as well as a neutron... the proton
        // is plugged in with its equator pointing down [edge-on, disc
        // feeds the hole], the neutron with its pole pointing down";
        // "the proton and neutron in the original configuration of
        // Carbon didn't spread out. They stayed in line." meth.pdf shows
        // stronger neighbors (Oxygen) can BREAK Carbon into other forms
        // — the shape is environmental; this is the default. The old
        // "three alphas stacked" reading (nuclear.pdf) is kept as the
        // `tri_alpha` harness structure below: per haf.pdf a bare stack
        // channels weakly without polar plug fans, which is exactly why
        // it fought the binding model for two sessions.
        "carbon" => {
            let mut a = stack(2);
            // End proton of the 2-stack sits at ALPHA_PITCH/2 + 1.3;
            // the FUSED plug pair (session 34, see plug_pair) parks
            // PLUG_SEAT beyond it — the measured force equilibrium
            // (report_plug_well), not the stacked-alpha funnel mouth.
            let y = ALPHA_PITCH / 2.0 + 1.3 + PLUG_SEAT;
            a.push(plug_pair(-y, plug_gap));
            a.push(plug_pair(y, plug_gap));
            Some(a)
        }
        // The bare three-alpha stack — NOT an element (session-32): kept
        // as the canonical stability-harness structure (the hard case
        // every binding regression showed up in) and as haf.pdf's
        // negative exemplar: a stack with no polar plug fans "will be
        // channeling weakly... relying only on ambient field potential".
        // Per haf.pdf bare 1-3 stacks ARE viable weak channelers
        // (helium is Mathis's own example) — what they must NOT drive is
        // the choice of global defaults (session-33 harness note).
        "tri_alpha" => Some(stack(3)),
        // The bare FOUR-alpha stack — haf.pdf's free negative
        // prediction: "Once you stack four alphas, the external charge
        // overwhelms the charge being channeled, and the nucleus can't
        // hold together." The model SHOULD dissolve this once
        // side-charge (Phase C ambient surface effects) is real;
        // report_quad_alpha tracks whether it does (session 33: it does
        // NOT yet — the side-charge term is what's missing).
        "quad_alpha" => Some(stack(4)),
        // Nitrogen: three-alpha core + 7th proton plugged in the south
        // pole (edge-on, disc feeding the hole) and the balancing
        // neutron in the north (ammon.pdf), pole-down per graphene.pdf.
        "nitrogen" => {
            let mut a = stack(3);
            // End proton sits at ALPHA_PITCH + 1.3; the plug parks one
            // funnel mouth (2.6) beyond it.
            a.push(plug_proton(-(ALPHA_PITCH + 3.9), 0.0));
            a.push(plug_neutron(ALPHA_PITCH + 3.9, 0.0));
            Some(a)
        }
        // Oxygen: three-alpha core + BOTH poles capped by a
        // proton+neutron PAIR (oxygen.pdf: the 7th and 8th protons go on
        // the ends because four alphas can't stack; atmo2.pdf: their
        // neutrons are paired with them in the hole). Each pair sits
        // side by side — proton edge-on (disc feeds the hole), neutron
        // pole-down (channels axially).
        "oxygen" => {
            let mut a = stack(3);
            // End proton sits at ALPHA_PITCH + 1.3; the FUSED plug pair
            // parks PLUG_SEAT beyond it (measured equilibrium — see
            // carbon above and report_plug_well).
            let y = ALPHA_PITCH + 1.3 + PLUG_SEAT;
            a.push(plug_pair(-y, plug_gap));
            a.push(plug_pair(y, plug_gap));
            Some(a)
        }
        // Neon: THE first carousel configuration (nuclear.pdf) — one
        // center alpha with four carousel alphas plugged edge-to-hole
        // around its equator. The axial charge hole top and bottom is
        // "surrounded by four charge maxima": unreactive, six-sided.
        "neon" => Some(vec![
            core_alpha(0.0, anatomy),
            carousel_alpha(0.0, anatomy),
            carousel_alpha(90.0, anatomy),
            carousel_alpha(180.0, anatomy),
            carousel_alpha(270.0, anatomy),
        ]),
        // Argon: Neon's carousel + the full axial line on the same center
        // disk — "nine disks... Level one is the center disk. Level two
        // consists of the four carousel disks. Level 3a is the posts up
        // and down. Level 3b is the caps top and bottom" (nuclear.pdf).
        // The connectors (3a) sit SIDEWAYS like the polar plugs; the caps
        // (3b) are parallel to the core.
        "argon" => Some(vec![
            core_alpha(-2.0 * ALPHA_PITCH, anatomy),  // cap (parallel to core)
            sideways_alpha(-ALPHA_PITCH, anatomy),    // connector (edge-on)
            core_alpha(0.0, anatomy),                 // center
            sideways_alpha(ALPHA_PITCH, anatomy),     // connector (edge-on)
            core_alpha(2.0 * ALPHA_PITCH, anatomy),   // cap
            carousel_alpha(0.0, anatomy),
            carousel_alpha(90.0, anatomy),
            carousel_alpha(180.0, anatomy),
            carousel_alpha(270.0, anatomy),
        ]),
        _ => None,
    }
}

/// Preset names for UI listings. `tri_alpha` is the bare 3-stack harness
/// structure (not an element) — listed so the sandbox can compare it
/// against the corrected plugged carbon live.
pub fn preset_names() -> &'static [&'static str] {
    &["alpha", "carbon", "tri_alpha", "nitrogen", "oxygen", "neon", "argon"]
}

/// One alpha unit of a nucleus: rolls rigidly about its OWN axis
/// (`roll_phase`), and — if non-core — additionally has its center+axis
/// revolve around the nucleus axis (`RigidGroup::carousel_phase`,
/// gated by `orbits_core`). Two independent kinematic DOFs replace the
/// single `car_rot` channel that used to do both jobs at once and could
/// only get one of them right per alpha (docs/ATOM_ROTATION_AND_SIM_DESIGN.md
/// §0). See `AtomCore::sync_group_members` for the exact composition.
pub struct AlphaUnit {
    /// Indices into the group's flat `members`/`local_offsets`/etc arrays
    /// (k, NOT particle ids — particle id = `group.members[members[i]]`).
    pub members: Vec<usize>,
    /// Each member's rest position in the ALPHA's own frame (+Y = axis),
    /// parallel to `members`.
    pub member_local_pos: Vec<DVec3>,
    /// Each member's pole in the alpha frame, parallel to `members`.
    pub member_local_pole: Vec<DVec3>,
    /// Rest axis of this alpha in the NUCLEUS frame (Y core/cap, radial u
    /// carousel, X sideways connector, edge-on/pole-on for plugs).
    pub rest_axis: DVec3,
    /// Rest center of this alpha in the NUCLEUS frame.
    pub rest_center: DVec3,
    /// Does this alpha's center+axis orbit the nucleus axis (carousel)?
    pub orbits_core: bool,
    /// Kinematic roll about the alpha's own axis (rad, wall clock).
    pub roll_phase: f64,
    pub roll_rate: f64,

    // ── Part 2 promotes the alpha to a real body (unused in RigidLock,
    // docs/ATOM_ROTATION_AND_SIM_DESIGN.md Part 2) — populated at spawn
    // so the fields exist and are sane, but nothing reads them yet. ──
    pub com: DVec3,
    pub velocity: DVec3,
    pub orientation: DQuat,
    pub angular_velocity: DVec3,
    /// Gyroscopic precession rate (session 33, wig.pdf): when
    /// `AtomCore::gyro_spin` > 0, the transverse component of this
    /// alpha's torque produces this orientation DRIFT instead of
    /// accumulating transverse angular momentum — recomputed fresh every
    /// `kick`, applied in `drift`, never integrated ("think of a
    /// spinning wheel, which resists being pushed sideways").
    pub precession: DVec3,
    pub mass: f64,
    pub inertia: f64,

    /// Three-ring skin radii (session-32, user-proposed): quantitative,
    /// computed ONCE at spawn from this alpha's OWN reach surface
    /// (`march_reach` over just its members, alpha frame — see
    /// `compute_alpha_ring_radii`). `ring_disc_r` = reach-surface radius
    /// at the proton disc planes (y = ±NUCLEON_PITCH/2); `ring_mid_r` =
    /// equatorial reach at the alpha midplane where the neutron posts
    /// live. 0.0 = not a standard 4-nucleon alpha, draw no rings. The
    /// ring GEOMETRY is anchored to live member particle positions in
    /// `build_group_alpha_rings`, so the readout works in every
    /// `NucleusDynamics` mode (closes the session-31 hidden-overlay gap).
    pub ring_disc_r: f64,
    pub ring_mid_r: f64,
}

/// A rigidly-locked composite (nucleus). Constituents remain real particles
/// (forces sampled per-constituent — an electron can capture at one specific
/// proton's pole), but they integrate as one rigid body. Intra-group pair
/// forces are skipped: the nucleus is pre-fused, its internal balance is
/// not simulated (uf4.pdf: "the alphas can't be broken and rearranged").
pub struct RigidGroup {
    /// Alpha sub-units — the source of truth for kinematics
    /// (`sync_group_members`). The flat arrays below are REST-POSE data
    /// derived from `alphas` once at spawn (roll = 0, carousel phase = 0)
    /// and never mutated after — they exist so `segment_sources`,
    /// `march_reach`, `compute_disc_exits`, `build_group_carousel_overlay`,
    /// and the flow-VFX emitters (all azimuth-averaged or drawn in a
    /// phase-rotated child node) don't need to change for Part 1
    /// (session-31 addendum A1).
    pub alphas: Vec<AlphaUnit>,
    pub members: Vec<usize>,
    pub local_offsets: Vec<DVec3>,
    pub local_orients: Vec<DQuat>,
    /// Intrinsic axial spin rate of each member (rad/s about its own pole).
    pub member_spin: Vec<f64>,
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
    /// Skin segments: the nucleus is fused from distinct pieces (axial
    /// stacks, sideways connectors/plugs, the carousel level), and each
    /// draws its OWN shell — a single averaged envelope smears them into
    /// an uninformative blob (session-29 user note on Neon).
    pub skin_segments: Vec<SkinSegment>,
}

/// One piece of a nucleus for skin purposes.
pub struct SkinSegment {
    /// Member indices (k into members/local_offsets/local_orients).
    pub members: Vec<usize>,
    /// Lathe center (group-local y); carousel segment sits at 0.
    pub y_center: f64,
    /// The equatorial carousel level (drawn as the wide disc).
    pub carousel: bool,
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

/// One member as a charge-field source (group-local), for skin marching.
struct FieldSrc {
    pos: DVec3,
    pole: DVec3,
    mass: f64,
    radius: f64,
    profile: usize,
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
        /// True for paths that END at the field-reach boundary (free
        /// loops, group through-paths, group disc exits) — these are
        /// eligible to leave a `SkinMote` behind when they die. False for
        /// electron-capture paths, which end at the rider, not the skin.
        dies_at_skin: bool,
    },
    /// A boundary mote: session-31 replacement for the translucent lathe
    /// "skin" surface around a composite nucleus. Flow parcels already
    /// die exactly at the field-reach boundary, so instead of a closed
    /// surface (which the user flagged as reading like a "gold sausage"
    /// blob), a fraction of those deaths leave a brief bright dot pinned
    /// AT the death point — hundreds of these per second sketch the
    /// reach surface as a living point-cloud, without ever drawing a
    /// continuous membrane. `local` is fixed in the anchor's frame (the
    /// mote rides the nucleus/particle); `swirl` is a slow residual spin
    /// about the anchor's Y axis, same formula as `Flow`'s whole-path
    /// swirl but at the mote's own (much slower) rate.
    SkinMote {
        anchor: VfxAnchor,
        local: DVec3,
        swirl: f64,
    },
}

/// What a flow parcel's path is expressed relative to.
#[derive(Clone, Copy)]
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
/// Charge-field lock on a single ALPHA in `NucleusDynamics::RigidAlpha`
/// (1/s) — the Part 2 analog of `GROUP_SPIN_RELAX`, but deliberately much
/// weaker: in this mode the nucleus `RigidGroup` is just a container, and
/// the carousel motion is supposed to emerge from real inter-alpha forces
/// (docs/ATOM_ROTATION_AND_SIM_DESIGN.md §2.2 — "you *want* the carousel
/// to be able to turn here"). This only bleeds off numerical spin-up, not
/// the emergent ω.
pub const ALPHA_SPIN_RELAX: f64 = 0.5;

/// Session-33 gyroscopic spin stiffness default — ships OFF (0.0,
/// classical torque response) until `report_gyro_sweep` earns a value
/// against plugged-carbon retention AND the tri_alpha battery. See
/// `AtomCore::gyro_spin` for the mechanism (wig.pdf).
pub const DEFAULT_GYRO_SPIN: f64 = 0.0;

/// Translational charge-field lock on RigidAlpha nuclei (1/s): damps each
/// alpha's velocity toward the nucleus's mass-weighted mean velocity —
/// relative motion only, whole-nucleus translation is untouched. The same
/// ambient lock that suppresses tumble (nuclear.pdf: "The charge field
/// then locks them into these configurations") dissipates internal
/// oscillation; without it the bound stack's breathing mode is UNDAMPED
/// and the velocity-dependent pair terms slowly pump it to ejection at a
/// ~160k-step horizon (session-31 round-5 "Jenga" reproduction —
/// dev_probe_dissolve.gd / probe_app_loop_carbon_rigidalpha).
pub const ALPHA_TRANS_RELAX: f64 = 0.5;
/// VISIBLE axial-spin rate of bodies (rad/s WALL clock — substep
/// independent), far below the physical TAU·3 rad/sim-s so the eye can
/// track it. For free particles this drives `display_orientation`; for
/// group members it is each alpha's `roll_phase` rate — a rigid roll of
/// an off-axis disc both spins the protons in place and revolves the
/// posts around the axis (docs/ATOM_ROTATION_AND_SIM_DESIGN.md §1.1).
/// Initial phases are randomized so alphas don't turn in lockstep.
pub const DISPLAY_SPIN_RATE: f64 = 0.9;
/// VISIBLE carousel ride rate (rad/s WALL clock) — independent of
/// `DISPLAY_SPIN_RATE`. Roll (an alpha spinning about its own axis) and
/// orbit (a non-core alpha's center+axis revolving around the nucleus
/// axis) are separate kinematic DOFs; locking their rates equal made the
/// whole nucleus read as one solid gear (session-31).
pub const CAROUSEL_VIS_RATE: f64 = 0.35;
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

/// Default intra-nucleus coupling multiplier (Part 2, session-31 addendum
/// A9). ORIGINALLY the sole nuclear-binding mechanism — every charge-field
/// pairwise term on a non-skipped same-group pair multiplied uniformly by
/// this ("push on each other like the H₂ free protons, only stronger").
/// The `alpha_stays_bound` sweep (A11) proved a UNIFORM boost cannot bind a
/// nucleus: it scales attraction and repulsion equally, so it only
/// rescales the existing molecular repulsive standoff instead of moving
/// it — drift got monotonically WORSE with boost (1233% at ×1 to 2056% at
/// ×12; see the historical sweep table this const's old value came from).
/// Concluded (A13, session-31 round 3): binding needs channeling
/// attenuation of `c_q`/`stream`, not a uniform multiplier. This const
/// drops to 1.0 (no-op); it stays as a `Couplings` field
/// (`intra_nucleus_boost`) so the harness can still sweep it as an
/// orthogonal knob on top of channeling.
pub const INTRA_NUCLEUS_BOOST: f64 = 1.0;

// ── Force couplings ──────────────────────────────────────────────────────

#[derive(Clone, Copy, Debug)]
pub struct Couplings {
    /// Gravity coupling (1/r², isotropic, attractive).
    pub g_q: f64,
    /// Charge coupling (1/r⁴, emission × absorption, repulsive).
    pub c_q: f64,
    /// Ambient isotropic charge pressure (pushes into charge shadows).
    /// Molecular-scale environment term; same-group pairs use
    /// `nuclear_ambient` instead (A13).
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
    /// Intra-nucleus coupling multiplier for non-skipped same-group pairs
    /// in `NucleusDynamics::RigidAlpha` / `FreeNucleon` (Part 2, session-31
    /// addendum A9). See `INTRA_NUCLEUS_BOOST` for the default/rationale;
    /// a `Couplings` field (not a bare const use site) so the
    /// `alpha_stays_bound` harness can sweep it at runtime. Session-31
    /// round 3 (A13) concluded a UNIFORM boost cannot bind a nucleus — it
    /// scales attraction and repulsion equally, so it only rescales the
    /// existing molecular repulsive standoff instead of moving it. The
    /// default drops to 1.0 (no-op); the field stays as a sweep hook for
    /// `c_q`/`stream` alongside channeling attenuation (below).
    pub intra_nucleus_boost: f64,
    /// Channeling attenuation strength for non-skipped same-group pairs
    /// (A13, session-31 round 3). Plugged neighbors route charge through
    /// each other's pole channel instead of colliding equator-to-equator:
    /// bb2.pdf — "attraction must always be explained as loss of
    /// repulsion" (spin cancellations lower the between-field's repulsive
    /// energy); strong.html — "charge is channeled through the nucleus by
    /// baryon spin, and so does not cause a repulsion between protons...
    /// There is no charge field within the nucleus." At `channeling = 1.0`
    /// a pole-on pair within the funnel mouth (r ≤ NUCLEON_PITCH) has its
    /// `c_q`/`stream` terms fully cancelled; two facing equators (the
    /// molecular repel configuration) are unaffected (C≈0). Default 0.85
    /// (round 5, re-swept after the rest-entry seeding fix — see
    /// `Couplings::default` doc) — the `alpha_stays_bound` sweep winner:
    /// full cancellation (1.0) plus any real ambient glue collapses the
    /// stack into the contact wall (see `Couplings::default` doc); the
    /// 15% residual repulsion is the cushion that keeps the bound stack
    /// off contact.
    pub channeling: f64,
    /// Nuclear ambient charge pressure (A13, session-31 round 3): the
    /// same-group replacement for `ambient_pressure` in the existing
    /// shadow term `P·(1−E_i)(1−E_j)·(1−occ)/r²` — nuclear.pdf: "the
    /// charge field is both the initial pressure and the subsequent
    /// glue," i.e. the field ambient to the fused nucleus pushes bodies
    /// into each other's charge shadows, the "glue" on top of channeling's
    /// "loss of repulsion". Default 2.0 (unchanged across round 3/4/5) —
    /// the `alpha_stays_bound` sweep winner (paired with `channeling =
    /// 0.85` as of round 5); 0 ejects, ≥5 with full channeling collapses
    /// (see `Couplings::default` doc).
    pub nuclear_ambient: f64,
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
    /// - `intra_nucleus_boost = 1.0`, `channeling = 0.85`,
    ///   `nuclear_ambient = 2.0` — nuclear-binding-only terms (A13, session-
    ///   31 round 3; RETUNED round 4; RETUNED round 5), gated to
    ///   non-skipped same-group pairs only; they do not touch the
    ///   molecular defaults above.
    ///
    ///   Round 3 shipped `channeling = 0.9`, `nuclear_ambient = 2.0` (2.5%
    ///   quiet drift, best of 9 of 15 passing combos on the OLD
    ///   channeling-falloff cliff at `r = 2·NUCLEON_PITCH`). Round 4
    ///   (session-31, binding-v2 hardening) found that cliff was the root
    ///   cause of two user-observed failures — a spontaneous late "Jenga"
    ///   collapse, and trivial destruction by a stray particle spawned
    ///   nearby — because any alpha displaced past the cliff regained full
    ///   molecular repulsion with no restoring basin (see `CHANNEL_TAIL`'s
    ///   doc). Extending the falloff to `CHANNEL_TAIL·NUCLEON_PITCH`
    ///   changes the force balance enough to warrant a fresh 3×4 sweep —
    ///   channeling {0.85, 0.9, 0.95} (1.0 EXCLUDED: with the longer tail
    ///   it collapses catastrophically once ambient ≥ ~5, drift in the
    ///   millions of percent, same disc-aware-contact-spring mechanism as
    ///   round 3) × `nuclear_ambient` {2, 5, 10, 20} — plus a SECOND pass
    ///   dimension: of the 12/12 quiet passers (all ≤30% drift over 20k
    ///   steps), which also survive an external flyby + parked-neighbor
    ///   perturbation (`nucleus_survives_flyby`'s scenario, full 15k+15k).
    ///   Round 4 picked `channeling=0.9, nuclear_ambient=5.0` (14.9% quiet /
    ///   13.5% flyby margin) over the more-tempting `channeling=0.85,
    ///   ambient=2.0` because at the time the latter read a SUSTAINED
    ///   63.5% flyby drift at full duration (only 1.8% on a short 8k+8k
    ///   proxy that later proved an unreliable window).
    ///
    ///   Round 5 (session-31 round 5, "Jenga" root-cause fix) found that
    ///   63.5%-drift result was itself an artifact: `alpha_kinematic_state`
    ///   seeded each alpha's `angular_velocity`/`velocity` with the
    ///   cosmetic DISPLAY roll/carousel rate on top of the real group
    ///   motion, injecting fictional spin energy whose slow decay
    ///   (`ALPHA_SPIN_RELAX`, τ=2s) walked the equilibrium spacing and made
    ///   the flyby result seed/timing-dependent. With that seeding fix
    ///   (real group kinematics only — see the function's doc), the round-4
    ///   3×4 sweep was re-run in full: `channeling=0.85, nuclear_ambient=2.0`
    ///   now has by far the best worst-case margin — 1.5% quiet drift, 2.7%
    ///   flyby drift — vs. the round-4 pick's 14.4%/17.2% and every other
    ///   combo in the 9/12 quiet-and-flyby-robust set (0.95 at ambient 5/10
    ///   still ejects on the flyby; 0.9/20 ejects too). `alpha_stays_bound`
    ///   asserts these defaults are that sweep's chosen winner AND holds
    ///   over a confirming 250k-step quiet run (past the ~160k-step
    ///   undamped-breathing horizon the seeding fix and `ALPHA_TRANS_RELAX`
    ///   jointly address); `nucleus_survives_flyby` re-validates the
    ///   flyby/swat case directly; `alpha_transient_survives_all_seeds`
    ///   re-validates across 8 distinct roll/carousel seed phases (the
    ///   axis the round-4 "Jenga" report actually varied along).
    ///
    ///   Session 32 (whirl-instability fix) RETUNED `nuclear_ambient`
    ///   2.0 → 5.0. Attenuating same-group `corot`/`vortex` by the
    ///   channeling factor (see the attenuation block in `compute_forces`)
    ///   killed the slow corot-driven whirl instability that dissolved
    ///   every seed at 225k-825k steps (user-observed "carbon dissolves at
    ///   ~300-400 sim time"), but also removed the incidental corot drag
    ///   that had been aiding flyby/swat recovery: at ambient=2.0 the
    ///   quiet drift stays a superb 1.5% but flyby drift degrades to
    ///   33.2%; ambient=5.0 balances both at 12.3%/11.9% — the sweep
    ///   chooser's worst-case-margin winner (11 of 12 quiet passers are
    ///   also flyby-robust post-fix).
    ///
    ///   Session 32 Phase B (through-charge flow tension,
    ///   docs/THROUGH_CHARGE_DESIGN.md): with stream tension +
    ///   channel-alignment stiffness active on same-group pairs
    ///   (flow_tension=1.0, flow_align=0.5), the 20k sweep's winner moved
    ///   to `channeling=0.9, nuclear_ambient=2.0` (7.8%/7.8%) — but that
    ///   combo FAILED the 8-seed × 1M-step long-horizon gate (4/8 seeds
    ///   collapse at 475k-825k), so it is on the sweep chooser's
    ///   LONG_HORIZON_VETO list and the defaults stayed at the then
    ///   gate-surviving `channeling=0.85, nuclear_ambient=5.0`
    ///   (12.3%/11.9% with tension). Lesson encoded in the chooser:
    ///   short-horizon margins alone must never ship a default.
    ///
    ///   Session 33 (tension SIGN FIX, THROUGH_CHARGE_DESIGN.md
    ///   addendum) RETUNED `nuclear_ambient` 5.0 → 10.0. Fixing the
    ///   inverted j-emits tension branch reshaped the long-horizon
    ///   landscape: the post-fix chooser winner (0.85, 2.0) failed the
    ///   1M gate (2/8), and the old (0.85, 5.0) default now fails it
    ///   too (1/8, k=0 collapse at 700k) — both vetoed. (0.85, 10.0)
    ///   is the next-best robust combo (16.5% quiet / 15.9% flyby) and
    ///   passes the gate 8/8, every seed converging to the same cold
    ///   compressed attractor (spacings 3.24/6.48, KE ≈ 5.6e-4,
    ///   perfect axis alignment, channel at max).
    ///
    ///   Session 34 (flow-network NO-STARVE fix, charge_flow.rs pass 2)
    ///   RETUNED `nuclear_ambient` 10.0 → 5.0 — back to the session-32
    ///   value. The starved network had been systematically
    ///   under-delivering stream glue (carbon ran at mult 0.53-0.80;
    ///   tension measured half its intended strength at the plug
    ///   sockets), and ambient 10.0 was compensating for exactly that
    ///   deficit. With flows floored at free-field ambient (deut.pdf
    ///   leaky hose; graphene.pdf "fans"), the full 5-combo 1M re-gate
    ///   found: (0.85, 5.0) 8/8 all-cold (spacings 3.34/6.68,
    ///   KE ≈ 5e-4) AND best surviving short-horizon margins
    ///   (13.7%/12.5%) → WINNER; (0.85, 10.0) also 8/8 all-cold but
    ///   worse margins (16.8%/15.5%); (0.9, 5.0) 2/8 collapse;
    ///   (0.9, 2.0) 2/8 collapse; (0.85, 2.0) bounded but 3/8 heat
    ///   monotonically to KE ≈ 128 by 1M (the pre-collapse whirl
    ///   signature) — all three vetoed with v2 evidence.
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
            intra_nucleus_boost: INTRA_NUCLEUS_BOOST,
            channeling: 0.85,
            nuclear_ambient: 5.0,
        }
    }
}

/// Nucleus dynamics mode (Part 2, session-31 addendum A9): a GLOBAL setting
/// on `AtomCore` (not per-group — per-group mixing is not needed for the
/// sandbox), exposed to Godot as `AtomSim::set_nucleus_dynamics(i32)` /
/// `get_nucleus_dynamics() -> i32` (0/1/2). Controls both the intra-nucleus
/// force skip predicate (`AtomCore::compute_forces`) and which level
/// integrates rigidly in `kick`/`drift` (docs/ATOM_ROTATION_AND_SIM_DESIGN.md
/// Part 2, §2.1–§2.2).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum NucleusDynamics {
    /// Current + Part 1 kinematics: the whole nucleus is pre-fused and
    /// phase-driven (uf4.pdf: "the alphas can't be broken and rearranged").
    /// Presentable default.
    RigidLock,
    /// Alphas are rigid bodies; forces act BETWEEN alphas within a nucleus
    /// (boosted by `Couplings::intra_nucleus_boost`). The nucleus
    /// `RigidGroup` becomes a container only — no longer integrated.
    RigidAlpha,
    /// Every nucleon is a free force participant (boosted, same as above);
    /// intra-alpha cohesion — if any — comes from the SAME force model, not
    /// a kinematic constraint. The research mode for confirming the alpha
    /// itself is a force equilibrium (§2.4).
    FreeNucleon,
}

impl Default for NucleusDynamics {
    fn default() -> Self {
        NucleusDynamics::RigidLock
    }
}

/// Default share of skin-boundary deaths that leave a mote behind — the
/// density knob for the boundary point-cloud (session-31: replaces the
/// lathe skin surface entirely). Promoted from a function-local const to a
/// runtime `AtomCore::mote_fraction` field (session-31 round 4, work item
/// 5) so the debug HUD can tune it live; this const stays as the seeded
/// default.
pub const DEFAULT_MOTE_FRACTION: f64 = 0.6;
/// Default mote lifetime (s) — long enough to read as a lingering point on
/// the boundary, short enough that the cloud stays "live". See
/// `DEFAULT_MOTE_FRACTION` for why this is now a runtime field
/// (`AtomCore::mote_lifetime`).
pub const DEFAULT_MOTE_LIFETIME: f64 = 2.8;
/// Default mote size multiplier. Motes are deliberately oversized
/// stand-ins for the hundreds of real parcels they represent — user
/// request, RX 550 budget (a few hundred big dots read better than
/// thousands of true-scale ones on this hardware). See
/// `DEFAULT_MOTE_FRACTION` for why this is now a runtime field
/// (`AtomCore::mote_scale`).
pub const DEFAULT_MOTE_SCALE: f32 = 3.0;

// ── AtomCore: the whole simulation, engine-free ──────────────────────────

pub struct AtomCore {
    pub profiles: Vec<ParticleProfile>,
    pub particles: Vec<SimParticle>,
    pub groups: Vec<RigidGroup>,
    /// Per-group scratch (indexed by group id): the net INTRA-group pairwise
    /// force accumulated during `compute_forces`. An isolated nucleus cannot
    /// self-propel — cc.pdf: the charge field is OPEN, and the momentum
    /// imbalance from the proton→neutron emission asymmetry (and any other
    /// un-paired internal term) is radiated away, isotropically on average,
    /// so it must not translate the COM. `step()` cancels it mass-weighted
    /// (translation only; internal torques/shape forces are preserved).
    /// Sized to `groups.len()` each `compute_forces`; empty otherwise.
    pub(crate) group_self_force: Vec<DVec3>,
    pub couplings: Couplings,
    /// Global nucleus dynamics mode (Part 2 sandbox toggle). See
    /// `NucleusDynamics`; default `RigidLock`.
    pub dynamics: NucleusDynamics,
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
    /// Runtime knobs for `advance_clouds`' skin-mote conversion (session-31
    /// round 4, work item 5 — debug-UI sliders). Seeded from
    /// `DEFAULT_MOTE_FRACTION`/`DEFAULT_MOTE_SCALE`/`DEFAULT_MOTE_LIFETIME`.
    pub mote_fraction: f64,
    pub mote_scale: f32,
    pub mote_lifetime: f64,
    /// Internal neutron-post pole anatomy for alpha builders (session-32
    /// experiment — see [`PostAnatomy`]). Takes effect at `spawn_preset`;
    /// changing it does not retrofit already-spawned groups.
    pub post_anatomy: PostAnatomy,
    /// Physics step counter — cadence source for the through-charge flow
    /// sweep (`charge_flow::FLOW_SOLVE_EVERY`).
    pub(crate) step_count: u64,
    /// Phase B flow couplings (session 32, docs/THROUGH_CHARGE_DESIGN.md):
    /// stream-tension and channel-alignment stiffness, applied to
    /// non-skipped same-group pairs in `compute_forces`. Runtime fields so
    /// the calibration sweep and the debug panel can vary them.
    pub flow_tension: f64,
    pub flow_align: f64,
    /// Phase C1 ambient-confinement pressure (session 32 — see
    /// `apply_ambient_confinement`). Ships 0.0 (off) until the
    /// calibration re-earns the battery at a nonzero value.
    pub ambient_confine: f64,
    /// Session-35 directional ambient charge field. The ambient bath is not
    /// isotropic — on Earth it "is pointing straight up everywhere"
    /// (pasta.pdf); nuclei "align their intake vortices to incoming charge"
    /// and "turn to stand straight up, to align to E" (dielec.pdf;
    /// raman.pdf: proper orientation is "inline with it"). `ambient_align`
    /// is the aligning-torque strength; `ambient_charge_dir` is the field
    /// AXIS (normalized on use; the pole-to-pole channel is bidirectional,
    /// so a nucleus aligns to ±dir, whichever is nearer). Both ship 0 (the
    /// isotropic, no-preferred-orientation status quo) until calibrated —
    /// see `apply_ambient_alignment`.
    pub ambient_align: f64,
    pub ambient_charge_dir: DVec3,
    /// Session-33 gyroscopic spin stiffness (wig.pdf: "The increased
    /// angular momentum acts to prevent the protons from turning...
    /// Think of a spinning wheel, which resists being pushed sideways").
    /// Each alpha unit carries spin angular momentum `gyro_spin × mass`
    /// along its own axis; in `kick` (RigidAlpha) the TRANSVERSE torque
    /// component then produces gyroscopic PRECESSION (an instantaneous
    /// orientation drift, `AlphaUnit::precession`) instead of
    /// accumulating transverse angular momentum — which removes both the
    /// near-zero-inertia overshoot of 1-nucleon plug alphas and the
    /// torque-war energy pump (report_plug_torque_war/-energy_ablation).
    /// Ships 0.0 (off — classical response) until calibrated.
    pub gyro_spin: f64,
    /// Session-33 plug-experiment knobs (report_plug_matrix): runtime
    /// fields so the factored sweep can vary them without rebuilds; ALL
    /// default to the earned baseline.
    ///
    /// Polar plug-pair half-gap at spawn (see PLUG_PAIR_GAP).
    pub plug_pair_gap: f64,
    /// Close-range saturation floor for the pairwise nuclear-ambient
    /// shadow glue: its 1/r² is evaluated at max(r, this). 0.0 = off
    /// (baseline); NUCLEON_PITCH = saturate inside the funnel mouth
    /// (far-field argument — the shadow glue overestimates once the
    /// pair excludes ambient photons from its own gap).
    pub ambient_sat_r: f64,
    /// Flow-align receiver torque target: false = line to the emitter
    /// (baseline), true = the emitter's pole-port STREAM direction
    /// (deut.pdf "align to that charge stream" — removes the ~17°
    /// frustration of the off-axis plug pair; see flow_tension_pair).
    pub align_to_stream: bool,
    /// Phase B.2 (session 34): throughput-scaled suction gain — the
    /// intake pull and channeling force from a same-group body scale by
    /// `1 + gain·(mult − 1)` of its live `FlowState::mult`
    /// (charge_flow::DEFAULT_FLOW_SUCTION; phos.pdf socket suction).
    pub flow_suction: f64,
    /// Phase B.2: same form on the charge push — a fed emitter pushes
    /// harder (charge_flow::DEFAULT_FLOW_EMIT_SCALE).
    pub flow_emit_scale: f64,
    /// Session-35 INTER-atom bonding vortex. `flow_suction` above is
    /// same-group (internal, earned at 0). This is its external twin: a
    /// nucleon channeling the whole stack's through-charge has an intake
    /// vortex far stronger than its own spin (user), reaching PAST its own
    /// disc to draw in a bonding partner. Scales the intake/channeling/
    /// corotation pull an INTER-group emitter exerts by `1 + gain·(mult−1)`
    /// of its live `FlowState::mult` — emission is NOT amplified (we make
    /// nuclei SUCK, not push). A free particle has mult=1 so H₂ (both ends
    /// free) is untouched; only multi-nucleon nuclei reach out. Ships 0
    /// until the OH/O₂ bond earns it. See `compute_forces`.
    pub bond_suction: f64,
    /// Session-35 Phase C: INTER-group through-charge tension gain. The
    /// same `flow_tension_pair` glue that binds a nucleus's own posts
    /// (atom_core.rs same-group branch) IS the covalent/H bond — poll.pdf:
    /// "water/water bonds ... are just channels in the stream. You could
    /// say that about any bond." water2.pdf: a bonding proton "isn't
    /// plugging into an alpha, it is simply aligning itself to the charge
    /// field" — which is exactly this force's channel tension + its cos²/
    /// sin² restoring torques. Same-group tension is unscaled (gain 1);
    /// this scales the INTER-group copy by `bond_tension`. It only bites
    /// where a real channel exists: the force is proportional to the
    /// emitter's `out_pole` flow, ~0 without the directional network and
    /// only significant once a stacked nucleus builds mult>1 through-charge
    /// (CFM_AMBFLOW). bond_suction draws the partner in from long range;
    /// this locks it at the seat (capture_falloff maxes r=7.8). A free
    /// particle emits almost nothing polar, so H₂ and the locked M5
    /// molecular suite are untouched. Ships 0 until OH/O₂ earns it.
    pub bond_tension: f64,
    /// Session-35 cont-2: INTER-group bond dissipation (radial). The OH
    /// tension well (above) is conservative — report_oh_bond showed H
    /// falling in, overshooting the shallow minimum into the pole wall, and
    /// reflecting back out; nothing sheds the approach energy so it never
    /// settles. Mathis's charge field DOES dissipate: charge-collision "spin
    /// damping" (neut2.pdf) and cog-opposition drag when a body moves
    /// against a charge stream (venus2.pdf: "like opposite cogs meeting …
    /// causes slowing"; matched/co-moving profiles feel NO drag, which is
    /// why this can't fight the co-rotating carousel — eccen.html: the
    /// charge wind "resists eccentricity, pushing orbits back toward
    /// equilibrium"). Modeled as a drag on the RADIAL relative velocity of
    /// an inter-group pair, gated by `capture_falloff` (same bond range as
    /// the tension) and applied equal-and-opposite (momentum conserved).
    /// Two atoms at rest in a bond have zero relative velocity → no drag,
    /// so it damps the APPROACH without decaying the bond. INTER-group only
    /// → the intra-nucleus carousel is untouched. Ships 0.
    pub bond_damp: f64,
    /// Session-35 cont-3: TANGENTIAL companion to `bond_damp`. REFUTED as a
    /// capture aid — kept as a documented dead-end (ships 0 → radial-only).
    /// The idea: radial damping alone can't capture because H escapes
    /// TANGENTIALLY, so also drag the component of inter-group relative
    /// velocity PERPENDICULAR to the line. But `v_rel` is measured against
    /// each individual (spinning, orbiting) oxygen NUCLEON, not the local
    /// charge stream — so this "drag" pulls H toward each nucleon's orbital
    /// velocity, i.e. it PUMPS H up to the nucleus's rotation and slings it
    /// off. report_oh_bond (session-35 cont-3): dt=20 ejected H to d≈3800-4300
    /// vs radial-only d≈1300; tangential-only flung it out immediately. The
    /// real lesson: the ejector is the `corot`/`vortex` sprinkler, and the
    /// bond is a spin-SYNC capture problem, not a tension-well + drag problem
    /// (the well merely feeds H into the sprinkler faster the deeper it gets).
    /// A correct dissipation must be relative to the CHARGE STREAM (venus2.pdf
    /// co-moving = no drag), not the nucleon — future work.
    pub bond_damp_tan: f64,
    /// Session-35 cont-2 diagnostic knob: scales the INTER-group tension
    /// ALIGN TORQUE only (not the tension force). Ships 1.0 (no change to
    /// the bond_tension behavior). report_oh_bond showed the inter-group
    /// align torque pumps rotational energy (non-conservative at this dt) —
    /// setting this to 0 isolates the tension FORCE well + bond_damp to
    /// test whether the torque is the sole capture blocker. Same-group
    /// torque is always unscaled.
    pub bond_align: f64,
    /// Flow-solver cadence in physics steps (session-34 knob, seeded
    /// from charge_flow::FLOW_SOLVE_EVERY). The Phase B forces read
    /// flow amplitudes up to this many steps stale; the lag makes them
    /// slightly non-conservative, and the no-starve network's larger
    /// amplitudes raised the stakes — 1 = solve every step (the
    /// lag-free reference for the energy-pump ablation).
    pub flow_solve_every: usize,
    /// Session-34 experiment (ships false): slave plug-pair (≤2-member
    /// alpha) orientation in RigidAlpha — graphene.pdf: the carbon
    /// plugs "stayed in line". See the kick() block for the rationale.
    pub plug_orient_lock: bool,
}

impl Default for AtomCore {
    fn default() -> Self {
        Self::new()
    }
}

/// Channeling attenuation factor `C` for a non-skipped same-group pair
/// (A13, session-31 round 3): `C = channeling · max(cos²θ_i, cos²θ_j) ·
/// c_dist(r)`. `cos_theta_i`/`cos_theta_j` are each `pole·d̂` toward the
/// OTHER particle (the same convention `compute_forces` already uses for
/// `cos_theta_i`/`cos_theta_j`) — a pole-on presentation (either side)
/// drives the `cos²` term to 1 (routes charge through the hole instead of
/// colliding), while two facing equators (cos ≈ 0 on both sides, the
/// molecular repel configuration) drive it to 0 regardless of distance.
/// `c_dist` is a smoothstep from 1 at `r ≤ NUCLEON_PITCH` (inside the
/// funnel mouth) to 0 at `r ≥ CHANNEL_TAIL·NUCLEON_PITCH` (7.8) — channeling
/// is a short-range, plugged-neighbor effect, not a whole-nucleus one.
///
/// Session-31 round 4: the falloff used to hit 0 at exactly `2·NUCLEON_PITCH`
/// (5.2), a hard cliff — a displaced alpha that crossed it regained full
/// repulsion instantly with nothing pulling it back, producing one-way
/// ejections (the user's "Jenga" collapse and flyby fragility; see
/// `CHANNEL_TAIL`'s doc). Stretching the taper to `CHANNEL_TAIL` gives
/// excursions past the old cutoff a restoring basin instead of a wall.
fn channeling_factor(channeling: f64, cos_theta_i: f64, cos_theta_j: f64, r: f64) -> f64 {
    let lo = NUCLEON_PITCH;
    let hi = CHANNEL_TAIL * NUCLEON_PITCH;
    let t = ((r - lo) / (hi - lo)).clamp(0.0, 1.0);
    let c_dist = 1.0 - t * t * (3.0 - 2.0 * t); // smoothstep, 1 at r=lo, 0 at r=hi
    let pole_align = cos_theta_i.powi(2).max(cos_theta_j.powi(2));
    channeling * pole_align * c_dist
}

impl AtomCore {
    pub fn new() -> Self {
        Self {
            profiles: Vec::new(),
            particles: Vec::new(),
            groups: Vec::new(),
            group_self_force: Vec::new(),
            couplings: Couplings::default(),
            ambient_align: 0.0,
            ambient_charge_dir: DVec3::ZERO,
            dynamics: NucleusDynamics::RigidLock,
            ambient_gravity: DVec3::ZERO,
            ambient_charge: DVec3::ZERO,
            running: false,
            dt: 0.0005,
            time: 0.0,
            rng: Rng::new(0xDEAD_BEEF_CAFE),
            vfx_particles: Vec::new(),
            vfx_enabled: true,
            vfx_time: 0.0,
            mote_fraction: DEFAULT_MOTE_FRACTION,
            mote_scale: DEFAULT_MOTE_SCALE,
            mote_lifetime: DEFAULT_MOTE_LIFETIME,
            post_anatomy: PostAnatomy::default(),
            step_count: 0,
            flow_tension: crate::charge_flow::DEFAULT_FLOW_TENSION,
            flow_align: crate::charge_flow::DEFAULT_FLOW_ALIGN,
            ambient_confine: crate::charge_flow::DEFAULT_AMBIENT_CONFINE,
            gyro_spin: DEFAULT_GYRO_SPIN,
            plug_pair_gap: PLUG_PAIR_GAP,
            ambient_sat_r: 0.0,
            align_to_stream: false,
            flow_suction: crate::charge_flow::DEFAULT_FLOW_SUCTION,
            flow_emit_scale: crate::charge_flow::DEFAULT_FLOW_EMIT_SCALE,
            bond_suction: 0.0,
            bond_tension: 0.0,
            bond_damp: 0.0,
            bond_damp_tan: 0.0,
            bond_align: 1.0,
            flow_solve_every: crate::charge_flow::FLOW_SOLVE_EVERY,
            plug_orient_lock: false,
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
        // Random display phase so bodies never turn in lockstep.
        let phase0 = self.rng.next_f64() * TAU;

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
            alpha: None,
            display_orientation: (orientation * DQuat::from_rotation_y(phase0)).normalize(),
            flow: crate::charge_flow::FlowState::default(),
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
        let alpha_specs = preset_alphas(name, self.post_anatomy, self.plug_pair_gap)?;
        let orientation = orientation_from_pole(axis);
        let gid = self.groups.len();

        let mut members = Vec::new();
        let mut local_offsets = Vec::new();
        let mut local_orients = Vec::new();
        let mut member_spin = Vec::new();
        let mut mass = 0.0;
        let mut inertia = 0.0;
        let mut alphas: Vec<AlphaUnit> = Vec::new();

        for (alpha_idx, spec) in alpha_specs.iter().enumerate() {
            // Rest pose (roll = 0, carousel = 0): the alpha's own frame
            // mapped straight into the nucleus frame by its rest_axis.
            let r_a = orientation_from_pole(spec.rest_axis);
            let mut a_members = Vec::new();
            let mut a_local_pos = Vec::new();
            let mut a_local_pole = Vec::new();
            let mut a_mass = 0.0;
            let mut a_inertia = 0.0;
            for ms in &spec.members {
                let profile_id = self.profile_id_by_name(ms.profile_name)?;
                let prof_mass = self.profiles[profile_id].mass;
                let prof_radius = self.profiles[profile_id].radius;
                let spin = default_spin_rate(ms.profile_name) * ms.spin_sign;
                let nucleus_local_pos = spec.rest_center + r_a * ms.local_pos;
                let nucleus_local_pole = r_a * ms.local_pole;
                let world_pos = pos + orientation * nucleus_local_pos;
                let world_pole = orientation * nucleus_local_pole;
                let id =
                    self.spawn_particle_ex(profile_id, world_pos, vel, world_pole, spin)?;
                self.particles[id].group = Some(gid);
                self.particles[id].alpha = Some(alpha_idx);
                let k = members.len();
                members.push(id);
                local_offsets.push(nucleus_local_pos);
                local_orients.push(orientation_from_pole(nucleus_local_pole));
                member_spin.push(spin);
                mass += prof_mass;
                inertia += prof_mass
                    * (nucleus_local_pos.length_squared() + 0.4 * prof_radius * prof_radius);
                a_members.push(k);
                a_local_pos.push(ms.local_pos);
                a_local_pole.push(ms.local_pole);
                a_mass += prof_mass;
                a_inertia +=
                    prof_mass * (ms.local_pos.length_squared() + 0.4 * prof_radius * prof_radius);
            }
            // Each alpha's roll starts at a random phase so alphas don't
            // turn in lockstep (the "synchronized ballet", session-29);
            // members WITHIN one alpha share its roll — that coherence is
            // correct (they're one rigid piece).
            let roll_phase = self.rng.next_f64() * TAU;
            alphas.push(AlphaUnit {
                members: a_members,
                member_local_pos: a_local_pos,
                member_local_pole: a_local_pole,
                rest_axis: spec.rest_axis,
                rest_center: spec.rest_center,
                orbits_core: spec.orbits_core,
                roll_phase,
                roll_rate: DISPLAY_SPIN_RATE,
                com: spec.rest_center,
                velocity: DVec3::ZERO,
                orientation: orientation_from_pole(spec.rest_axis),
                angular_velocity: DVec3::ZERO,
                precession: DVec3::ZERO,
                mass: a_mass,
                inertia: a_inertia.max(1e-9),
                ring_disc_r: 0.0,
                ring_mid_r: 0.0,
            });
        }

        let carousel_phase = self.rng.next_f64() * TAU;
        self.groups.push(RigidGroup {
            alphas,
            members,
            local_offsets,
            local_orients,
            member_spin,
            // Carousel orbit — DISPLAY rate, wall clock, readable, substep
            // independent (docs/ATOM_ROTATION_AND_SIM_DESIGN.md §1.4).
            carousel_rate: CAROUSEL_VIS_RATE,
            carousel_phase,
            com: pos,
            velocity: vel,
            orientation,
            angular_velocity: DVec3::ZERO,
            mass,
            inertia: inertia.max(1e-9),
            skin_reach: Vec::new(),
            disc_exits: Vec::new(),
            skin_segments: Vec::new(),
        });
        self.sync_group_members(gid);
        let reach = self.compute_group_reach(gid, SKIN_REACH_RINGS);
        self.groups[gid].skin_reach = reach;
        let exits = self.compute_disc_exits(gid);
        self.groups[gid].disc_exits = exits;
        let segments = self.compute_skin_segments(gid);
        self.groups[gid].skin_segments = segments;
        self.compute_alpha_ring_radii(gid);
        Some(gid)
    }

    /// Spawn-time pass for the three-ring alpha skin (session-32): for
    /// every standard 4-nucleon alpha (2 protons + posts), march THIS
    /// alpha's own reach surface (its members only, in the alpha's own
    /// frame — `member_local_pos`/`member_local_pole` are already
    /// alpha-frame) and store two quantitative radii on the [`AlphaUnit`]:
    /// the reach radius at the proton disc planes and the equatorial reach
    /// at the midplane. Radii are static spawn data (the reach surface is
    /// a rest-pose property, like the group `skin_reach`); the rings'
    /// live placement comes from member particle state in
    /// `build_group_alpha_rings`.
    fn compute_alpha_ring_radii(&mut self, gid: usize) {
        let n_alphas = self.groups[gid].alphas.len();
        for ai in 0..n_alphas {
            let srcs: Vec<FieldSrc> = {
                let g = &self.groups[gid];
                let a = &g.alphas[ai];
                let proton_count = a
                    .members
                    .iter()
                    .filter(|&&k| {
                        self.profiles[self.particles[g.members[k]].profile_id].name
                            == "proton"
                    })
                    .count();
                if proton_count < 2 {
                    continue;
                }
                a.members
                    .iter()
                    .enumerate()
                    .map(|(i, &k)| {
                        let profile = self.particles[g.members[k]].profile_id;
                        FieldSrc {
                            pos: a.member_local_pos[i],
                            pole: a.member_local_pole[i],
                            mass: self.profiles[profile].mass,
                            radius: self.profiles[profile]
                                .radius
                                .max(MIN_RENDER_RADIUS as f64),
                            profile,
                        }
                    })
                    .collect()
            };
            let table = self.march_reach(&srcs, DVec3::ZERO, SKIN_REACH_RINGS);
            if table.len() < 2 {
                continue;
            }
            let mid_r = reach_at_theta(&table, std::f64::consts::PI / 2.0);
            // Disc ring: where the reach surface crosses the (north)
            // proton disc plane y = +NUCLEON_PITCH/2 — same crossing
            // logic as `compute_disc_exits`, in the alpha frame.
            let yc = NUCLEON_PITCH / 2.0;
            let n = table.len();
            let mut disc_r = 0.0f64;
            for i in 0..(n - 1) {
                let th0 = std::f64::consts::PI * i as f64 / (n - 1) as f64;
                let th1 = std::f64::consts::PI * (i + 1) as f64 / (n - 1) as f64;
                let y0 = table[i] * th0.cos();
                let y1 = table[i + 1] * th1.cos();
                if (y0 - yc) * (y1 - yc) <= 0.0 && (y0 - y1).abs() > 1e-12 {
                    let f = ((y0 - yc) / (y0 - y1)).clamp(0.0, 1.0);
                    let th = th0 + (th1 - th0) * f;
                    disc_r = reach_at_theta(&table, th) * th.sin();
                    break;
                }
            }
            let a = &mut self.groups[gid].alphas[ai];
            a.ring_mid_r = mid_r;
            a.ring_disc_r = disc_r;
        }
    }

    /// Split a group into skin segments: alphas whose `rest_center` sits
    /// off the stack axis (lateral > 1.5 — compares PLUG_PAIR_GAP 0.8 vs
    /// CAROUSEL_R 4.5, so it's no longer a per-member magic number,
    /// session-31 addendum A4) all belong to ONE carousel segment; the
    /// rest cluster along the axis with a gap threshold that keeps a
    /// contiguous alpha stack together. Intra-alpha member gaps are ≤ 1.3
    /// (proton↔post), inter-piece gaps ≥ 1.75 (argon center↔connector),
    /// so 1.5 splits pieces and keeps a contiguous carbon stack whole.
    pub fn compute_skin_segments(&self, group_idx: usize) -> Vec<SkinSegment> {
        const SEG_GAP: f64 = 1.5;
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        let mut carousel: Vec<usize> = Vec::new();
        let mut axial: Vec<(f64, usize)> = Vec::new();
        for a in &g.alphas {
            let lateral = (a.rest_center.x * a.rest_center.x + a.rest_center.z * a.rest_center.z)
                .sqrt();
            if lateral > 1.5 {
                carousel.extend(a.members.iter().copied());
            } else {
                for &k in &a.members {
                    axial.push((g.local_offsets[k].y, k));
                }
            }
        }
        axial.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        let mut segments: Vec<SkinSegment> = Vec::new();
        let mut cluster: Vec<(f64, usize)> = Vec::new();
        let flush = |cluster: &mut Vec<(f64, usize)>, segments: &mut Vec<SkinSegment>| {
            if cluster.is_empty() {
                return;
            }
            let y_center =
                cluster.iter().map(|c| c.0).sum::<f64>() / cluster.len() as f64;
            segments.push(SkinSegment {
                members: cluster.iter().map(|c| c.1).collect(),
                y_center,
                carousel: false,
            });
            cluster.clear();
        };
        for (y, k) in axial {
            if let Some(&(last_y, _)) = cluster.last() {
                if y - last_y > SEG_GAP {
                    flush(&mut cluster, &mut segments);
                }
            }
            cluster.push((y, k));
        }
        flush(&mut cluster, &mut segments);
        if !carousel.is_empty() {
            segments.push(SkinSegment {
                members: carousel,
                y_center: 0.0,
                carousel: true,
            });
        }
        segments
    }

    /// Reposition a group's members from each alpha's frame and give them
    /// the rigid-body velocity field (v_com + ω×r) plus their intrinsic
    /// axial spin AND the kinematic display terms (roll, carousel) —
    /// velocity-dependent forces (doppler, corotation) on constituents
    /// need correct member velocities. Exact kinematics from
    /// docs/ATOM_ROTATION_AND_SIM_DESIGN.md §1.3 + session-31 addendum A3:
    ///
    /// ```text
    /// R_a   = orientation_from_pole(a.rest_axis)      // alpha +Y → rest_axis
    /// Roll  = DQuat::from_rotation_y(a.roll_phase)     // about alpha-local +Y
    /// Car   = if a.orbits_core { rotation_y(carousel_phase) } else { IDENTITY }
    ///
    /// inner       = R_a * (Roll * member_local_pos)     // roll revolves posts
    /// nucleus_pos = Car * (a.rest_center + inner)        // carousel orbits piece
    /// q_nucleus   = Car * R_a * Roll * orientation_from_pole(member_local_pole)
    ///
    /// p.position    = com + nucleus_orientation * nucleus_pos
    /// p.orientation = (nucleus_orientation * q_nucleus).normalize()
    /// ```
    ///
    /// Two independent rotational DOFs per alpha (roll about its OWN axis,
    /// carousel about the NUCLEUS axis) replace the single `car_rot`
    /// channel that used to try to do both jobs at once and could only get
    /// one of them right per alpha (§0 of the design doc).
    fn sync_group_members(&mut self, gid: usize) {
        struct AlphaSync {
            rest_axis: DVec3,
            rest_center: DVec3,
            orbits_core: bool,
            roll_phase: f64,
            roll_rate: f64,
            member_ids: Vec<usize>,
            local_pos: Vec<DVec3>,
            local_pole: Vec<DVec3>,
            spins: Vec<f64>,
        }
        let (alpha_data, com, vel, nucleus_orientation, omega, carousel_phase, carousel_rate) = {
            let g = &self.groups[gid];
            let alpha_data: Vec<AlphaSync> = g
                .alphas
                .iter()
                .map(|a| AlphaSync {
                    rest_axis: a.rest_axis,
                    rest_center: a.rest_center,
                    orbits_core: a.orbits_core,
                    roll_phase: a.roll_phase,
                    roll_rate: a.roll_rate,
                    member_ids: a.members.iter().map(|&k| g.members[k]).collect(),
                    local_pos: a.member_local_pos.clone(),
                    local_pole: a.member_local_pole.clone(),
                    spins: a.members.iter().map(|&k| g.member_spin[k]).collect(),
                })
                .collect();
            (
                alpha_data,
                g.com,
                g.velocity,
                g.orientation,
                g.angular_velocity,
                g.carousel_phase,
                g.carousel_rate,
            )
        };
        for a in &alpha_data {
            let r_a = orientation_from_pole(a.rest_axis);
            let roll = DQuat::from_rotation_y(a.roll_phase);
            let car = if a.orbits_core {
                DQuat::from_rotation_y(carousel_phase)
            } else {
                DQuat::IDENTITY
            };
            // Alpha axis/center in WORLD space, riding the carousel — used
            // by the roll velocity term below (addendum A3).
            let axis_w = nucleus_orientation * (car * a.rest_axis);
            let center_w = com + nucleus_orientation * (car * a.rest_center);
            let car_w = if a.orbits_core {
                (nucleus_orientation * DVec3::Y) * carousel_rate
            } else {
                DVec3::ZERO
            };
            for i in 0..a.member_ids.len() {
                let inner = r_a * (roll * a.local_pos[i]);
                let nucleus_pos = car * (a.rest_center + inner);
                let q_nucleus = car * r_a * roll * orientation_from_pole(a.local_pole[i]);
                let world_off = nucleus_orientation * nucleus_pos;

                let p = &mut self.particles[a.member_ids[i]];
                p.position = com + world_off;
                p.orientation = (nucleus_orientation * q_nucleus).normalize();

                let mut v = vel + omega.cross(world_off);
                let mut ang = omega + p.pole_axis() * a.spins[i];
                if a.orbits_core {
                    v += car_w.cross(world_off);
                    ang += car_w;
                }
                let roll_w = axis_w * a.roll_rate;
                v += roll_w.cross(p.position - center_w);
                ang += roll_w;
                p.velocity = v;
                p.angular_velocity = ang;
            }
        }
    }

    // ── Part 2: nucleus dynamics modes (docs/ATOM_ROTATION_AND_SIM_DESIGN.md
    // §2, session-31 addendum A9/A10) ───────────────────────────────────

    pub fn get_nucleus_dynamics(&self) -> NucleusDynamics {
        self.dynamics
    }

    /// Switch nucleus dynamics mode, seeding/derived state so the next
    /// `step()` is physically continuous in the chosen mode (addendum A10):
    /// - **→ RigidAlpha / FreeNucleon:** every alpha's rigid-body state
    ///   (`com`, `orientation`, `velocity`, `angular_velocity`) is captured
    ///   from its group's CURRENT phase-driven kinematic pose
    ///   (`roll_phase`/`carousel_phase`). Those phases are frozen while
    ///   away from RigidLock (`advance_display` gates the group branch off
    ///   below), so this is well-defined and idempotent regardless of the
    ///   PREVIOUS mode. Member particles already sit at that exact pose
    ///   (from the last `sync_group_members` call), so nothing else needs
    ///   seeding for FreeNucleon; RigidAlpha additionally re-places members
    ///   from the freshly seeded alpha frame (algebraically a no-op here,
    ///   kept for robustness — see `place_alpha_members`). FreeNucleon also
    ///   resets group members' `display_orientation` to their current
    ///   physics orientation: `advance_display` applies the free-particle
    ///   display treatment to them too in this mode (so they keep visibly
    ///   spinning), and a STALE display frame (last touched at spawn) could
    ///   sit near-antiparallel to the current physics pole —
    ///   `from_rotation_arc` is only stable for nearby vectors.
    /// - **→ RigidLock:** resumes phase-driven sync (`sync_group_members`)
    ///   — positions SNAP to the kinematic pose. Documented sandbox
    ///   behavior (a HUD toggle would warn the user — A12, not this
    ///   commit).
    ///
    /// A no-op if `mode` is already the current mode.
    pub fn set_nucleus_dynamics(&mut self, mode: NucleusDynamics) {
        if mode == self.dynamics {
            return;
        }
        match mode {
            NucleusDynamics::RigidAlpha | NucleusDynamics::FreeNucleon => {
                for gid in 0..self.groups.len() {
                    let n_alphas = self.groups[gid].alphas.len();
                    for ai in 0..n_alphas {
                        let (com, orientation, velocity, angular_velocity) =
                            self.alpha_kinematic_state(gid, ai);
                        let a = &mut self.groups[gid].alphas[ai];
                        a.com = com;
                        a.orientation = orientation;
                        a.velocity = velocity;
                        a.angular_velocity = angular_velocity;
                    }
                    if mode == NucleusDynamics::RigidAlpha {
                        for ai in 0..n_alphas {
                            self.place_alpha_members(gid, ai);
                        }
                    }
                }
                if mode == NucleusDynamics::FreeNucleon {
                    for i in 0..self.particles.len() {
                        if self.particles[i].group.is_some() {
                            let phase0 = self.rng.next_f64() * TAU;
                            let p = &mut self.particles[i];
                            p.display_orientation =
                                (p.orientation * DQuat::from_rotation_y(phase0)).normalize();
                        }
                    }
                }
            }
            NucleusDynamics::RigidLock => {
                for gid in 0..self.groups.len() {
                    self.sync_group_members(gid);
                }
            }
        }
        self.dynamics = mode;
    }

    /// An alpha's rigid-body state (com, orientation, velocity,
    /// angular_velocity) reconstructed from its group's CURRENT
    /// phase-driven kinematic pose — the same composition
    /// `sync_group_members` uses (§1.3 + addendum A3), evaluated at the
    /// alpha's own center instead of per-member. Used by
    /// `set_nucleus_dynamics` to seed RigidAlpha/FreeNucleon.
    ///
    /// Derivation: `sync_group_members` shows
    /// `world_pos = center_w + alpha_orientation·member_local_pos` where
    /// `center_w = com + nucleus_orientation·(Car·rest_center)` and
    /// `alpha_orientation = nucleus_orientation·Car·R_a·Roll` — i.e. the
    /// alpha's own body frame factors cleanly out of the per-member
    /// formula. Position/orientation are seeded from that FULL kinematic
    /// pose (roll phase and carousel phase both included) — this is only
    /// a snapshot of where things are, not how fast they're "really"
    /// moving.
    ///
    /// Velocity/angular-velocity are seeded from the GROUP rigid-field
    /// term ONLY (`g.velocity`, `g.angular_velocity` composed at
    /// `center_w`) — deliberately dropping `roll_w` (per-alpha display
    /// spin) and `car_w` (display carousel orbit rate). Those are wall-
    /// clock READABILITY rates (see `ALPHA_ROLL_RATE`/carousel docs, A10
    /// display-channel note above `advance_display`), not physical
    /// momenta; injecting them as seeded angular/linear velocity pumps
    /// fictional energy into the rigid-body integrator, which the
    /// spin-relax damping (`ALPHA_SPIN_RELAX`, τ=2s) then bleeds off over
    /// ~100k+ steps — a slow, seed-phase-dependent transient that walks
    /// the equilibrium spacing and can shed an alpha (session-31 round-5
    /// "Jenga" root cause: the seeded roll/carousel rate, not the genuine
    /// channeling well at rest spacing, which is stable).
    fn alpha_kinematic_state(&self, gid: usize, ai: usize) -> (DVec3, DQuat, DVec3, DVec3) {
        let g = &self.groups[gid];
        let a = &g.alphas[ai];
        let r_a = orientation_from_pole(a.rest_axis);
        let roll = DQuat::from_rotation_y(a.roll_phase);
        let car = if a.orbits_core {
            DQuat::from_rotation_y(g.carousel_phase)
        } else {
            DQuat::IDENTITY
        };
        let nucleus_orientation = g.orientation;
        let center_w = g.com + nucleus_orientation * (car * a.rest_center);
        let center_off = center_w - g.com;

        // Physical seed: group translation + group rotation carried out
        // to this alpha's center. NO roll_w, NO car_w — those are the
        // cosmetic display rates, not physics (see doc comment above).
        let v = g.velocity + g.angular_velocity.cross(center_off);
        let ang = g.angular_velocity;

        let orientation = (nucleus_orientation * car * r_a * roll).normalize();
        (center_w, orientation, v, ang)
    }

    /// Place one alpha's members rigidly from its OWN integrated body frame
    /// (`com`/`orientation`/`velocity`/`angular_velocity`) — the
    /// `NucleusDynamics::RigidAlpha` analog of `sync_group_members`, driven
    /// by physics integration instead of phase variables (addendum A10).
    fn place_alpha_members(&mut self, gid: usize, ai: usize) {
        let (member_ids, local_pos, local_pole, spins, com, orientation, velocity, angular_velocity) = {
            let g = &self.groups[gid];
            let a = &g.alphas[ai];
            (
                a.members.iter().map(|&k| g.members[k]).collect::<Vec<usize>>(),
                a.member_local_pos.clone(),
                a.member_local_pole.clone(),
                a.members.iter().map(|&k| g.member_spin[k]).collect::<Vec<f64>>(),
                a.com,
                a.orientation,
                a.velocity,
                a.angular_velocity,
            )
        };
        for i in 0..member_ids.len() {
            let world_off = orientation * local_pos[i];
            let p = &mut self.particles[member_ids[i]];
            p.position = com + world_off;
            p.orientation = (orientation * orientation_from_pole(local_pole[i])).normalize();
            let v = velocity + angular_velocity.cross(world_off);
            // Rigid alpha ω + the member's own intrinsic axial spin, same
            // pattern as sync_group_members.
            let ang = angular_velocity + p.pole_axis() * spins[i];
            p.velocity = v;
            p.angular_velocity = ang;
        }
    }

    /// Advance the VISIBLE kinematic rotations by wall-clock `delta`:
    /// each alpha's roll phase, the shared carousel orbit, and free
    /// particles' display twist. These run at readable DISPLAY rates,
    /// deliberately far below the physical spin (TAU·3 rad/sim-s) and
    /// INDEPENDENT of the substep count — substeps scale physics time, not
    /// how fast the bodies visibly turn.
    pub fn advance_display(&mut self, delta: f64) {
        // Free particles always get the parallel-transported display
        // frame; in FreeNucleon, group members ALSO get it (addendum A10)
        // — rigid kinematic sync is suspended in that mode, so without
        // this they'd read as visibly frozen despite spinning physically.
        let free_nucleon = self.dynamics == NucleusDynamics::FreeNucleon;
        for p in &mut self.particles {
            if p.group.is_none() || free_nucleon {
                // Parallel-transport the render frame: nudge its Y onto
                // the current physics pole (consecutive poles are close,
                // so from_rotation_arc is stable), then spin about the
                // pole at the display rate.
                let pole = p.pole_axis();
                let cur_y = (p.display_orientation * DVec3::Y).normalize();
                let align = DQuat::from_rotation_arc(cur_y, pole);
                let sign = p.angular_velocity.dot(pole).signum();
                let spin =
                    DQuat::from_axis_angle(pole, sign * DISPLAY_SPIN_RATE * delta);
                p.display_orientation = (spin * align * p.display_orientation).normalize();
            }
        }
        // Group phase-driven kinematics (carousel orbit + per-alpha roll)
        // only apply in RigidLock — in RigidAlpha/FreeNucleon the group is
        // either a frozen container or a rigid body integrated in
        // `drift()`, and re-syncing from phases here would fight that
        // physics integration (addendum A10).
        if self.dynamics == NucleusDynamics::RigidLock {
            for gi in 0..self.groups.len() {
                {
                    let g = &mut self.groups[gi];
                    g.carousel_phase = (g.carousel_phase + g.carousel_rate * delta) % TAU;
                    for a in &mut g.alphas {
                        a.roll_phase = (a.roll_phase + a.roll_rate * delta) % TAU;
                    }
                }
                self.sync_group_members(gi);
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
        self.time = 0.0;
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

    /// Diagnostic: pairwise force-term breakdown between particles `i` and
    /// `j`, computed with the EXACT same skip/attenuation/boost logic
    /// `compute_forces` uses (current `self.dynamics`/`self.couplings`) —
    /// for harness failure analysis (A13 honesty clause: "print
    /// force_breakdown-style numbers at the failure step for the best
    /// combo"). Read-only; not used by the sim step itself.
    /// Returns `[r, channel, f_grav, f_charge_on_j, f_ambient,
    /// f_intake_on_j, f_stream, f_contact, f_tension]` — indices 2..=7 are
    /// UNSIGNED magnitudes of the term as it appears in `net_on_j`
    /// (gravity/ambient/intake pull j toward i; charge/stream/contact push
    /// j away from i); index 8 (session-32 flow tension) is SIGNED along
    /// d_hat (negative pulls j toward i).
    pub fn pair_force_breakdown(&self, i: usize, j: usize) -> [f64; 9] {
        let gi = self.particles[i].group;
        let gj = self.particles[j].group;
        let same_group = gi.is_some() && gi == gj;
        let base_cq = self.couplings;

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
        let cos_theta_i = pole_i.dot(d_hat);
        let cos_theta_j = pole_j.dot(-d_hat);

        let channel = if same_group {
            channeling_factor(base_cq.channeling, cos_theta_i, cos_theta_j, r)
        } else {
            0.0
        };
        let atten = (1.0 - channel) * base_cq.intra_nucleus_boost;
        let cq = if same_group {
            Couplings {
                g_q: base_cq.g_q,
                c_q: base_cq.c_q * atten,
                ambient_pressure: base_cq.nuclear_ambient,
                // Torque attenuated with the push — the gear-mesh torque
                // comes from the partner's emission COLLIDING with this
                // disc; a channeled pair routes it through instead (see
                // the compute_forces copy of this block for the full why).
                torque: base_cq.torque * atten,
                vortex: base_cq.vortex,
                drag: base_cq.drag,
                intake: base_cq.intake,
                corot: base_cq.corot,
                stream: base_cq.stream * atten,
                intra_nucleus_boost: base_cq.intra_nucleus_boost,
                channeling: base_cq.channeling,
                nuclear_ambient: base_cq.nuclear_ambient,
            }
        } else {
            base_cq
        };

        let emission_i = pi_prof.emission.sample(cos_theta_i);
        let emission_j = pj_prof.emission.sample(cos_theta_j);
        let absorption_j = pj_prof.absorption.sample(cos_theta_j);
        let absorption_i = pi_prof.absorption.sample(cos_theta_i);

        let v_rel = self.particles[j].velocity - self.particles[i].velocity;
        let v_radial = v_rel.dot(d_hat);
        let doppler = (1.0 - cq.drag * v_radial).clamp(0.2, 5.0);

        // Phase B.2 throughput scaling — lockstep with compute_forces.
        let (suction_i, epush_i) = if same_group {
            let mult_i = self.particles[i].flow.mult;
            (
                (1.0 + self.flow_suction * (mult_i - 1.0)).max(0.0),
                (1.0 + self.flow_emit_scale * (mult_i - 1.0)).max(0.0),
            )
        } else {
            (1.0, 1.0)
        };

        let f_grav = cq.g_q * pi_prof.mass * pj_prof.mass / r2s;
        let mass_prod = pi_prof.mass * pj_prof.mass;
        let f_charge_on_j =
            cq.c_q * mass_prod * emission_i * absorption_j / r4s * doppler * epush_i;

        let n = self.particles.len();
        let occlusion = self.compute_occlusion();
        let occ = occlusion[i * n + j];
        let shadow = (1.0 - emission_i) * (1.0 - emission_j);
        // Optional close-range saturation — lockstep with compute_forces.
        let r_amb = r.max(self.ambient_sat_r);
        let f_ambient = cq.ambient_pressure * shadow * (1.0 - occ) / (r_amb * r_amb);

        let ai2 = absorption_i * absorption_i;
        let f_intake_on_j =
            cq.intake * pi_prof.mass * pj_prof.mass * ai2 * suction_i * (1.0 - occ) / r2s;

        let aj2 = absorption_j * absorption_j;
        let ei2 = emission_i * emission_i;
        let ej2 = emission_j * emission_j;
        let f_stream = cq.stream
            * (mass_prod * mass_prod)
            * (ai2 * aj2 + ei2 * ej2)
            * (1.0 + 2.0 * occ)
            / r4s;

        let ri = pi_prof.radius.max(MIN_RENDER_RADIUS as f64);
        let rj = pj_prof.radius.max(MIN_RENDER_RADIUS as f64);
        let r_contact = if same_group {
            // Same-group rule incl. the session-32 neutron-as-rod
            // addendum — keep in lockstep with compute_forces.
            let sin_theta_i = (1.0 - cos_theta_i * cos_theta_i).max(0.0).sqrt();
            let sin_theta_j = (1.0 - cos_theta_j * cos_theta_j).max(0.0).sqrt();
            let ri_eff = if pi_prof.name == "neutron" {
                ri * POLE_HALF_THICKNESS
            } else {
                ri * (POLE_HALF_THICKNESS + (1.0 - POLE_HALF_THICKNESS) * sin_theta_i)
            };
            let rj_eff = if pj_prof.name == "neutron" {
                rj * POLE_HALF_THICKNESS
            } else {
                rj * (POLE_HALF_THICKNESS + (1.0 - POLE_HALF_THICKNESS) * sin_theta_j)
            };
            ri_eff + rj_eff
        } else {
            ri + rj
        };
        let f_contact = if r < r_contact {
            let overlap = r_contact - r;
            let mut fc = CONTACT_STIFFNESS * overlap;
            if v_radial < 0.0 {
                let m_red = mass_prod / (pi_prof.mass + pj_prof.mass);
                let c_n = 2.0 * (CONTACT_STIFFNESS * m_red).sqrt();
                fc -= c_n * v_radial;
            }
            fc
        } else {
            0.0
        };

        // Session-32 Phase B: signed flow-tension component on j along
        // d_hat (negative = pulls j toward i), zero for non-same-group
        // pairs (the term is gated in compute_forces).
        let f_tension = if same_group {
            let (f_on_j, _, _) = self.flow_tension_pair(i, j, d_hat, r);
            f_on_j.dot(d_hat)
        } else {
            0.0
        };


        [
            r,
            channel,
            f_grav,
            f_charge_on_j,
            f_ambient,
            f_intake_on_j,
            f_stream,
            f_contact,
            f_tension,
        ]
    }

    /// Diagnostic twin of `pair_force_breakdown` for the TORQUE channel
    /// (session-33 plug torque-war investigation): the two orientation
    /// authorities acting on a same-group pair, mirrored exactly from
    /// `compute_forces`. Returns
    /// `[charge_tau_i, charge_tau_j, flow_tau_i, flow_tau_j]` — the
    /// "equator toward charge" gear-mesh torque (attenuated by the
    /// channeling factor) and the flow-align torque (from
    /// `flow_tension_pair`) on each particle. Read-only.
    pub fn pair_torque_breakdown(&self, i: usize, j: usize) -> [DVec3; 4] {
        let gi = self.particles[i].group;
        let gj = self.particles[j].group;
        let same_group = gi.is_some() && gi == gj;
        let base_cq = self.couplings;

        let d_vec = self.particles[j].position - self.particles[i].position;
        let r = d_vec.length().max(SOFTENING);
        let r4s = (r * r) * (r * r);
        let d_hat = d_vec / r;

        let pi_prof = &self.profiles[self.particles[i].profile_id];
        let pj_prof = &self.profiles[self.particles[j].profile_id];
        let pole_i = self.particles[i].pole_axis();
        let pole_j = self.particles[j].pole_axis();
        let cos_theta_i = pole_i.dot(d_hat);
        let cos_theta_j = pole_j.dot(-d_hat);

        let (torque_coupling, c_q) = if same_group {
            let channel =
                channeling_factor(base_cq.channeling, cos_theta_i, cos_theta_j, r);
            let atten = (1.0 - channel) * base_cq.intra_nucleus_boost;
            (base_cq.torque * atten, base_cq.c_q * atten)
        } else {
            (base_cq.torque, base_cq.c_q)
        };

        let emission_i = pi_prof.emission.sample(cos_theta_i);
        let emission_j = pj_prof.emission.sample(cos_theta_j);
        let absorption_j = pj_prof.absorption.sample(cos_theta_j);
        let absorption_i = pi_prof.absorption.sample(cos_theta_i);
        let v_rel = self.particles[j].velocity - self.particles[i].velocity;
        let doppler = (1.0 - base_cq.drag * v_rel.dot(d_hat)).clamp(0.2, 5.0);
        let mass_prod = pi_prof.mass * pj_prof.mass;
        // Phase B.2 throughput scaling — lockstep with compute_forces.
        let (epush_i, epush_j) = if same_group {
            (
                (1.0 + self.flow_emit_scale * (self.particles[i].flow.mult - 1.0)).max(0.0),
                (1.0 + self.flow_emit_scale * (self.particles[j].flow.mult - 1.0)).max(0.0),
            )
        } else {
            (1.0, 1.0)
        };
        let f_charge_on_j =
            c_q * mass_prod * emission_i * absorption_j / r4s * doppler * epush_i;
        let f_charge_on_i =
            c_q * mass_prod * emission_j * absorption_i / r4s * doppler * epush_j;

        let charge_tau_j = pole_j.cross(-d_hat)
            * (-2.0 * pole_j.dot(-d_hat) * f_charge_on_j * torque_coupling);
        let charge_tau_i = pole_i.cross(d_hat)
            * (-2.0 * pole_i.dot(d_hat) * f_charge_on_i * torque_coupling);

        let (flow_tau_i, flow_tau_j) = if same_group {
            let (_, ti, tj) = self.flow_tension_pair(i, j, d_hat, r);
            (ti, tj)
        } else {
            (DVec3::ZERO, DVec3::ZERO)
        };

        [charge_tau_i, charge_tau_j, flow_tau_i, flow_tau_j]
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

        // 0. Through-charge flow relaxation sweep (session 32, Phase A —
        // docs/THROUGH_CHARGE_DESIGN.md). Runs BEFORE forces so Phase B
        // force terms will read a network consistent with the geometry
        // they act on. Every FLOW_SOLVE_EVERY steps: re-routing is not
        // instantaneous physically, and the stale window (4 ms sim time)
        // is far below any mechanical timescale here.
        self.step_count = self.step_count.wrapping_add(1);
        if self.step_count % self.flow_solve_every.max(1) as u64 == 0 {
            self.solve_charge_flow();
        }

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
        // Tumble/precession damping: free particles always; in FreeNucleon,
        // group members too — they integrate exactly like free particles in
        // that mode (kick/drift below), so they need the same damping or
        // their tumble would grow unbounded (nothing else damps it once the
        // group-level GROUP_SPIN_RELAX/ALPHA_SPIN_RELAX stop applying).
        let free_nucleon = self.dynamics == NucleusDynamics::FreeNucleon;
        for i in 0..self.particles.len() {
            if self.particles[i].group.is_some() && !free_nucleon {
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

    /// Velocity (and angular velocity) update for free particles and
    /// groups/alphas per `NucleusDynamics` (Part 2, §2.2). Free (ungrouped)
    /// particles integrate the same way in every mode; only how the GROUPED
    /// levels integrate changes.
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
        match self.dynamics {
            NucleusDynamics::RigidLock => {
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
            NucleusDynamics::RigidAlpha => {
                // The nucleus RigidGroup is a frozen CONTAINER in this
                // mode; each ALPHA integrates as its own rigid body —
                // aggregate its members' forces into an alpha COM force +
                // torque about the alpha's own com (§2.2).
                for gi in 0..self.groups.len() {
                    let n_alphas = self.groups[gi].alphas.len();
                    for ai in 0..n_alphas {
                        let (f, tau, mass, inertia) = {
                            let g = &self.groups[gi];
                            let a = &g.alphas[ai];
                            let mut f = DVec3::ZERO;
                            let mut tau = DVec3::ZERO;
                            for &k in &a.members {
                                let p = &self.particles[g.members[k]];
                                f += p.force_accum;
                                tau += (p.position - a.com).cross(p.force_accum)
                                    + p.torque_accum;
                            }
                            (f, tau, a.mass, a.inertia)
                        };
                        let a = &mut self.groups[gi].alphas[ai];
                        a.velocity += f / mass * dt;
                        if self.gyro_spin > 1e-12 {
                            // Gyroscopic response (session 33, wig.pdf):
                            // the alpha carries spin angular momentum
                            // S = gyro_spin·mass along its own axis, so a
                            // TRANSVERSE torque precesses the axis
                            // (Ω = â × τ⊥ / S, from τ = Ω × Sâ) instead of
                            // accelerating a tilt; only the AXIAL (roll)
                            // component integrates classically. Precession
                            // is recomputed fresh each kick — it is a
                            // response to the PRESENT torque, not a stored
                            // momentum, which is exactly why it cannot be
                            // pumped into an oscillation.
                            let axis = a.orientation * DVec3::Y;
                            let tau_par = axis * tau.dot(axis);
                            let tau_perp = tau - tau_par;
                            let s = self.gyro_spin * mass;
                            a.angular_velocity += tau_par / inertia * dt;
                            a.precession = axis.cross(tau_perp) / s;
                        } else {
                            a.angular_velocity += tau / inertia * dt;
                            a.precession = DVec3::ZERO;
                        }
                        // Weak per-alpha damping — well below the nucleus
                        // lock's GROUP_SPIN_RELAX (20.0): the carousel must
                        // stay free to turn here (§2.2); this only bleeds
                        // off numerical spin-up, not the emergent ω.
                        a.angular_velocity *= (-ALPHA_SPIN_RELAX * dt).exp();
                        // Session-34 experiment (ships false): plug pairs
                        // "stayed in line" (graphene.pdf) — a fused pair
                        // plugged into a socket channel does not tumble
                        // independently, so slave its orientation (the
                        // static well scans are stable in BOTH axes at
                        // the measured seat; free pair rotation is the
                        // remaining chaos channel, report_plug_retention).
                        if self.plug_orient_lock && a.members.len() <= 2 {
                            a.angular_velocity = DVec3::ZERO;
                            a.precession = DVec3::ZERO;
                        }
                    }
                    // TRANSLATIONAL charge-field lock: damp each alpha's
                    // velocity toward the nucleus's mass-weighted mean —
                    // the same ambient lock that suppresses tumble
                    // (nuclear.pdf: "The charge field then locks them into
                    // these configurations") also dissipates RELATIVE
                    // nucleon motion. Without it the compressed stack's
                    // breathing mode has ZERO damping, and the
                    // velocity-dependent pair terms (doppler, corotation)
                    // slowly pump it until an alpha crosses the channeling
                    // basin and ejects (observed horizon ~160k steps —
                    // the user-reported "Jenga" collapse, session 31
                    // round 5). Damping is RELATIVE, so whole-nucleus
                    // translation is untouched.
                    let (v_mean, m_tot) = {
                        let g = &self.groups[gi];
                        let mut p = DVec3::ZERO;
                        let mut m = 0.0;
                        for a in &g.alphas {
                            p += a.velocity * a.mass;
                            m += a.mass;
                        }
                        (p / m.max(1e-12), m)
                    };
                    if m_tot > 0.0 {
                        let decay = (-ALPHA_TRANS_RELAX * dt).exp();
                        for a in &mut self.groups[gi].alphas {
                            a.velocity = v_mean + (a.velocity - v_mean) * decay;
                        }
                    }
                }
            }
            NucleusDynamics::FreeNucleon => {
                // Every nucleon (grouped or not) integrates as a free
                // particle; the group is bookkeeping only (§2.2).
                for i in 0..self.particles.len() {
                    if self.particles[i].group.is_none() {
                        continue; // already integrated above
                    }
                    let p = &mut self.particles[i];
                    let inv_m = 1.0 / self.profiles[p.profile_id].mass;
                    p.velocity += p.force_accum * inv_m * dt;
                    p.angular_velocity += p.torque_accum * dt;
                }
            }
        }
    }

    /// Position (and orientation) update for free particles and
    /// groups/alphas per `NucleusDynamics` (Part 2, §2.2).
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
        match self.dynamics {
            NucleusDynamics::RigidLock => {
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
                        // Kinematic display phases (carousel ride, member
                        // axial spins) advance in WALL CLOCK via
                        // advance_display(), not here — visible rotation
                        // must not scale with substeps.
                    }
                    self.sync_group_members(gi);
                }
            }
            NucleusDynamics::RigidAlpha => {
                // The nucleus RigidGroup itself freezes (container only,
                // §2.2); each alpha integrates its own com/orientation and
                // places its members rigidly from that frame.
                for gi in 0..self.groups.len() {
                    let n_alphas = self.groups[gi].alphas.len();
                    for ai in 0..n_alphas {
                        let a = &mut self.groups[gi].alphas[ai];
                        a.com += a.velocity * dt;
                        // Gyroscopic precession (session 33) rides on top
                        // of the integrated angular velocity — a drift,
                        // not a momentum (see kick / AlphaUnit).
                        let w = a.angular_velocity + a.precession;
                        let w_len = w.length();
                        if w_len > 1e-12 {
                            let rot = DQuat::from_axis_angle(w / w_len, w_len * dt);
                            a.orientation = (rot * a.orientation).normalize();
                        }
                    }
                    for ai in 0..n_alphas {
                        self.place_alpha_members(gi, ai);
                    }
                }
            }
            NucleusDynamics::FreeNucleon => {
                // Every nucleon integrates freely; the group (com/
                // orientation) stays frozen and sync_group_members is NOT
                // called — rigid re-syncing would fight the free
                // integration (addendum A10).
                for i in 0..self.particles.len() {
                    if self.particles[i].group.is_none() {
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
            }
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

    // pub(crate) for the session-34 static well-scan diagnostic
    // (report_plug_well) — sim stepping still owns the call sequence.
    pub(crate) fn compute_forces(&mut self) {
        let n = self.particles.len();
        for p in &mut self.particles {
            p.force_accum = DVec3::ZERO;
            p.torque_accum = DVec3::ZERO;
        }
        // Reset the per-group intra-nucleus self-force accumulator (see the
        // field doc). Each non-skipped same-group pair adds its net applied
        // force here; `step()` cancels it so the nucleus can't self-propel.
        self.group_self_force.clear();
        self.group_self_force.resize(self.groups.len(), DVec3::ZERO);

        let base_cq = self.couplings;
        let dynamics = self.dynamics;
        let occlusion = self.compute_occlusion();

        // Pairwise forces
        for i in 0..n {
            for j in (i + 1)..n {
                let gi = self.particles[i].group;
                let gj = self.particles[j].group;
                let same_group = gi.is_some() && gi == gj;
                // Mode-aware intra-nucleus skip (Part 2, session-31
                // addendum A9):
                //   RigidLock   — skip any same-group pair. The nucleus is
                //                 pre-fused (uf4.pdf — "alphas can't be
                //                 broken and rearranged"); internal balance
                //                 isn't simulated (current/Part-1 behavior,
                //                 unchanged).
                //   RigidAlpha  — skip only same-ALPHA pairs: an alpha
                //                 stays pre-fused, but inter-alpha pairs
                //                 within the nucleus now interact (boosted
                //                 below) — "push on each other like H₂,
                //                 only stronger" (nuclear.pdf carousel).
                //   FreeNucleon — never skip; every nucleon is a free force
                //                 participant.
                let skip = same_group
                    && match dynamics {
                        NucleusDynamics::RigidLock => true,
                        NucleusDynamics::RigidAlpha => {
                            self.particles[i].alpha == self.particles[j].alpha
                        }
                        NucleusDynamics::FreeNucleon => false,
                    };
                if skip {
                    continue;
                }

                // Snapshot both members' force so we can attribute THIS
                // pair's net applied force to the group (intra-group pairs
                // only — a bound nucleus must not self-propel; see
                // `group_self_force`). Paired terms (corot, tension, contact)
                // net to zero here automatically; only the un-paired central
                // charge/intake asymmetry and vortex contribute.
                let f_snap_i = self.particles[i].force_accum;
                let f_snap_j = self.particles[j].force_accum;

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

                // Binding v2 (A13, session-31 round 3): a non-skipped
                // same-group pair (RigidAlpha inter-alpha pairs; any
                // FreeNucleon pair) gets `c_q`/`stream` ATTENUATED by the
                // channeling factor C (bb2.pdf: "attraction must always be
                // explained as loss of repulsion"; strong.html: "charge is
                // channeled through the nucleus by baryon spin... There is
                // no charge field within the nucleus") and `ambient_pressure`
                // replaced by `nuclear_ambient` (nuclear.pdf: "the charge
                // field is both the initial pressure and the subsequent
                // glue"). `intra_nucleus_boost` (default 1.0, a sweep hook
                // left over from the concluded uniform-boost experiment —
                // see its doc) still multiplies `c_q`/`stream` on top.
                // Everything else (g_q, torque, vortex, intake, corot,
                // drag) is left UNMULTIPLIED — intake IS the channeled flow
                // (diamag.pdf: polar protons are "fans, pulling charge in")
                // and is what holds the edge-to-hole carousel plugs; it
                // must not also be attenuated by channeling or it would
                // undermine the very mechanism it represents. The hard
                // contact-repulsion term further below reads
                // CONTACT_STIFFNESS directly (never `cq`), so it is
                // unaffected by any of this — fusion must not become a
                // spring-launcher.
                let cq = if same_group {
                    let channel = channeling_factor(
                        base_cq.channeling,
                        cos_theta_i,
                        cos_theta_j,
                        r,
                    );
                    let atten = (1.0 - channel) * base_cq.intra_nucleus_boost;
                    Couplings {
                        g_q: base_cq.g_q,
                        c_q: base_cq.c_q * atten,
                        ambient_pressure: base_cq.nuclear_ambient,
                        // Torque is attenuated with the push: the
                        // "equator toward charge" gear-mesh torque comes
                        // from the partner's emission COLLIDING with this
                        // disc — a channeled pair routes that charge
                        // through instead, and un-attenuated it slowly
                        // TILTS stacked alphas out of their pole-aligned
                        // plug until the channel degrades and repulsion
                        // wins (session-31 round-5 second dissolution
                        // mechanism, visible once ALPHA_TRANS_RELAX
                        // removed the fast breathing pump).
                        torque: base_cq.torque * atten,
                        // Corot and vortex are attenuated with the push for
                        // the same reason torque is: both are couplings to
                        // the partner's EMITTED swirl (gear-mesh with the
                        // collision field), and a channeled pair routes that
                        // charge through instead. Un-attenuated, corot's
                        // cross-coupled force on the closest facing-proton
                        // pairs taps the (kinematically re-imposed, hence
                        // inexhaustible) member_spin reservoir as a
                        // circulatory "whirl" instability: the stack sits
                        // just below threshold, the middle alpha's ~2x roll
                        // drive sweeps the relative roll phase across the
                        // stability boundary after a seed-dependent 200-800k
                        // steps, and the lateral shear mode then grows
                        // exponentially from the noise floor until the
                        // channeling tail releases the stored compression —
                        // the session-32 "carbon dissolves at 300-400 sim
                        // time" collapse (see report_long_horizon_drift /
                        // report_onset_zoom / report_stability_vs_time).
                        // Intake stays UN-attenuated: it IS the channeled
                        // flow (see the block comment above).
                        vortex: base_cq.vortex * atten,
                        drag: base_cq.drag,
                        intake: base_cq.intake,
                        corot: base_cq.corot * atten,
                        stream: base_cq.stream * atten,
                        intra_nucleus_boost: base_cq.intra_nucleus_boost,
                        channeling: base_cq.channeling,
                        nuclear_ambient: base_cq.nuclear_ambient,
                    }
                } else {
                    base_cq
                };

                // Phase B.2 (session 34, docs/THROUGH_CHARGE_DESIGN.md):
                // throughput-scaled force profiles. The static
                // emission/absorption tables describe a FREE-FIELD
                // particle; a body channeling flow above baseline sucks
                // and pushes harder in proportion (phos.pdf: an unfed
                // socket "doesn't have much pull, or suction"). Factor
                // `1 + gain·(mult − 1)` per body: exactly 1 for a free
                // particle (mult = 1 by construction) at any gain, so
                // the locked molecular tables only move where a live
                // network exists. Same-group gated like the rest of
                // Phase B; Phase C extends output effects to third
                // parties (auger.pdf).
                let (suction_i, suction_j, epush_i, epush_j) = if same_group {
                    let mult_i = self.particles[i].flow.mult;
                    let mult_j = self.particles[j].flow.mult;
                    (
                        (1.0 + self.flow_suction * (mult_i - 1.0)).max(0.0),
                        (1.0 + self.flow_suction * (mult_j - 1.0)).max(0.0),
                        (1.0 + self.flow_emit_scale * (mult_i - 1.0)).max(0.0),
                        (1.0 + self.flow_emit_scale * (mult_j - 1.0)).max(0.0),
                    )
                } else if self.bond_suction.abs() > 1e-12 {
                    // Session-35 INTER-atom bonding vortex: a nucleon
                    // channeling the whole stack's through-charge reaches
                    // out and SUCKS a partner in — intake/channeling/corot
                    // scaled by the emitter's throughput, emission left
                    // alone (epush = 1). Free partners (mult = 1) are
                    // untouched, so H₂ and lone particles are unaffected.
                    let mult_i = self.particles[i].flow.mult;
                    let mult_j = self.particles[j].flow.mult;
                    (
                        (1.0 + self.bond_suction * (mult_i - 1.0)).max(0.0),
                        (1.0 + self.bond_suction * (mult_j - 1.0)).max(0.0),
                        1.0,
                        1.0,
                    )
                } else {
                    (1.0, 1.0, 1.0, 1.0)
                };

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
                // × epush (Phase B.2): a fed emitter pushes harder.
                let mass_prod = pi_prof.mass * pj_prof.mass;
                let f_charge_on_j =
                    cq.c_q * mass_prod * emission_i * absorption_j / r4s * doppler * epush_i;
                let f_charge_on_i =
                    cq.c_q * mass_prod * emission_j * absorption_i / r4s * doppler * epush_j;

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
                    // i's emission gradient acting on j (× epush_i with
                    // the radial push it is the θ-component of)
                    let pole_eff_i = pole_i * cos_theta_i.signum();
                    let sin_i = (1.0 - cos_theta_i * cos_theta_i).max(0.0).sqrt();
                    if sin_i > 1e-6 {
                        let theta_hat = (d_hat * cos_theta_i.abs() - pole_eff_i) / sin_i;
                        let de = pi_prof.emission.d_dtheta(cos_theta_i);
                        let f_conf =
                            -cq.c_q * mass_prod * absorption_j * de * epush_i / (3.0 * r4s);
                        self.particles[j].force_accum += theta_hat * f_conf;
                        self.particles[i].force_accum -= theta_hat * f_conf;
                    }
                    // j's emission gradient acting on i
                    let pole_eff_j = pole_j * cos_theta_j.signum();
                    let sin_j = (1.0 - cos_theta_j * cos_theta_j).max(0.0).sqrt();
                    if sin_j > 1e-6 {
                        let theta_hat = (-d_hat * cos_theta_j.abs() - pole_eff_j) / sin_j;
                        let de = pj_prof.emission.d_dtheta(cos_theta_j);
                        let f_conf =
                            -cq.c_q * mass_prod * absorption_i * de * epush_j / (3.0 * r4s);
                        self.particles[i].force_accum += theta_hat * f_conf;
                        self.particles[j].force_accum -= theta_hat * f_conf;
                    }
                }

                // Channel occlusion for this pair (stoppered vortex).
                let occ = occlusion[i * n + j];

                // Ambient pressure: pushes particles into charge shadows.
                // A stoppered channel has no shadow minimum (diatom.pdf).
                // Optional close-range saturation (session-33 experiment,
                // `ambient_sat_r` — 0.0 = off): the shadow glue is a
                // far-field approximation, and un-saturated its 1/r² digs
                // a ~7-unit well between the pole-on carbon plug partners
                // (shadow ≈ 1, zero charge repulsion) with nothing but
                // the contact spring at the bottom (report_plug_retention).
                let shadow = (1.0 - emission_i) * (1.0 - emission_j);
                let r_amb = r.max(self.ambient_sat_r);
                let f_ambient =
                    cq.ambient_pressure * shadow * (1.0 - occ) / (r_amb * r_amb);

                // Polar intake: the proton recycles charge through its poles,
                // creating a focused inward flow.  A(θ)² sharpens the profile
                // to a narrow polar cone matching Mathis's vortex description.
                // A particle stoppering the channel blocks the stream.
                let ai2 = absorption_i * absorption_i;
                let aj2 = absorption_j * absorption_j;

                // Radial intake: pulls toward the emitter, 1/r² (flow sink).
                // × suction (Phase B.2): a FED funnel sucks harder than
                // its static absorption profile — the plug-retention
                // mechanism (phos.pdf; the static term measures ~0.016
                // at a carbon socket regardless of live flow).
                let f_intake_on_j = cq.intake * pi_prof.mass * pj_prof.mass
                    * ai2 * suction_i * (1.0 - occ) / r2s;
                let f_intake_on_i = cq.intake * pj_prof.mass * pi_prof.mass
                    * aj2 * suction_j * (1.0 - occ) / r2s;

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
                        // × suction: the same live vortex as the radial
                        // intake term (Phase B.2).
                        let f_chan = cq.intake * pi_prof.mass
                            * pj_prof.mass * ai2 * suction_i * lat_dist_i / (r * r2s);
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
                            * pi_prof.mass * aj2 * suction_j * lat_dist_j / (r * r2s);
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
                    // × suction_i: the bonding vortex that captures a
                    // partner also SYNCS it (session-35 — the spin-lock that
                    // holds A to B's vortex; same throughput scaling as the
                    // radial pull). Unity internally (flow_suction=0) and for
                    // free partners (mult=1).
                    let f_corot = (v_corot - v_tan_vec)
                        * (cq.corot * mass_prod * ai2 * suction_i / r2s);
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
                    let f_corot_i = (v_corot_j + v_tan_vec)
                        * (cq.corot * mass_prod * aj2 * suction_j / r2s);
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

                // Through-charge stream tension + channel-alignment
                // stiffness (session-32 Phase B —
                // docs/THROUGH_CHARGE_DESIGN.md): a same-group pair
                // joined by a live capture link is pulled together with
                // the link's flow, and both ends feel restoring torques
                // from the cos²/sin² capture gates (graphene.pdf: "the
                // vertical charge streams ... tie them together";
                // deut.pdf posts as self-balancing charge channels).
                // Equal-and-opposite by construction. Gated to
                // same-group pairs in Phase B so the locked molecular
                // tests are untouched; Phase C extends output-stream
                // effects to third parties (auger.pdf intruder
                // repulsion).
                // Same-group tension is the intra-nucleus post glue
                // (gain 1). Phase C (session 35): the SAME force binds an
                // external partner along a shared channel — poll.pdf's "you
                // could say that about any bond" — scaled by `bond_tension`
                // (ships 0, so the locked molecular suite is untouched).
                // Driven by the emitter's `out_pole`, it is ~0 without the
                // directional network and only grabs where a real through-
                // channel (mult>1) reaches the partner.
                let tension_scale = if same_group {
                    1.0
                } else {
                    self.bond_tension
                };
                if tension_scale.abs() > 1e-12 {
                    let (f_on_j, tau_i, tau_j) =
                        self.flow_tension_pair(i, j, d_hat, r);
                    self.particles[j].force_accum += f_on_j * tension_scale;
                    self.particles[i].force_accum -= f_on_j * tension_scale;
                    // Torque scaled additionally by bond_align for INTER-group
                    // (diagnostic: 0 = force-only well, isolates the torque
                    // pump). Same-group keeps the full torque.
                    let torque_scale = if same_group {
                        tension_scale
                    } else {
                        tension_scale * self.bond_align
                    };
                    self.particles[i].torque_accum += tau_i * torque_scale;
                    self.particles[j].torque_accum += tau_j * torque_scale;
                }

                // Inter-group bond dissipation (radial): the tension well is
                // conservative, so a captured partner overshoots and
                // reflects instead of settling. This bleeds the RADIAL
                // relative velocity within bond range (venus2.pdf cog-drag /
                // neut2.pdf spin damping — see `bond_damp`). Equal-and-
                // opposite (momentum conserved); zero at relative rest so a
                // seated bond doesn't decay. INTER-group only.
                // Baryon-only: an equal-and-opposite drag force on the
                // 1/1836-mass electron is a huge acceleration (numerically
                // stiff → the H proton/electron pair blows apart). The bond
                // drag is between the nucleons; the electron just rides.
                if !same_group
                    && (self.bond_damp.abs() > 1e-12 || self.bond_damp_tan.abs() > 1e-12)
                    && pi_prof.mass > 0.5
                    && pj_prof.mass > 0.5
                {
                    let fall = crate::charge_flow::capture_falloff(r);
                    if fall > 1e-6 {
                        let v_rel = self.particles[j].velocity
                            - self.particles[i].velocity;
                        let v_radial = v_rel.dot(d_hat);
                        let v_rad_vec = d_hat * v_radial;
                        let v_tan_vec = v_rel - v_rad_vec;
                        // Radial drag along the line + tangential drag on the
                        // perpendicular component (session-35 cont-3). Both
                        // equal-and-opposite; zero at relative rest.
                        let f_damp = v_rad_vec * (-self.bond_damp * fall)
                            + v_tan_vec * (-self.bond_damp_tan * fall);
                        self.particles[j].force_accum += f_damp;
                        self.particles[i].force_accum -= f_damp;
                    }
                }

                // Contact repulsion: hard-sphere boundary at sum of radii,
                // with near-critical normal damping while approaching —
                // an undamped spring against the 1/1836-mass electron
                // (dt·ω ≈ 0.2) scatters it chaotically ("bunny hopping").
                //
                // Disc-aware effective radii (A13, session-31 round 3): the
                // bodies are planar emitters, not spheres, so a pole-on
                // approach nests into the partner's funnel and only the
                // thin disc waist can collide, while an equator-on approach
                // presents the full radius (see `POLE_HALF_THICKNESS` doc
                // for the r_eff formula). The design doc calls for this
                // GLOBALLY (every pair, "a disc is a disc regardless of
                // fusion state") but that measurably broke the molecular
                // suite: `hydrogen_capture`, `derived_constants_equilibrium`,
                // `hydrogen_orbit_stable_long_run`, and
                // `alpha_captures_electron` all failed with disc-aware
                // contact applied to proton–electron pairs (thinning the
                // contact wall for a near-pole-on electron let it punch
                // through to unstable orbits/escape) — this could not be
                // honestly attributed to improved physics, it's a
                // regression of the LOCKED M5 molecular force table. Per
                // the doc's own fallback clause, disc-aware contact is
                // SAME-GROUP-GATED only (RigidAlpha/FreeNucleon intra-
                // nucleus pairs); every other pair keeps the plain
                // sphere-sum contact that the M5 lockdown tests pin.
                // Session-32 addendum to the same-group rule: NEUTRONS
                // collide as thin RODS, not discs — "protons are more 2D
                // while neutrons are more 1D... neutrons act more like
                // little lightning rods" (phos.pdf). The disc formula
                // gave a neutron a disc's FULL radius equator-on, which
                // put the two rest-pose posts (1.4 apart, ±POST_R) in
                // deep interpenetration: invisible while same-alpha
                // pairs are skipped, but FreeNucleon's first step fired
                // them out on a 58-unit compressed contact spring (the
                // real bulk of nucleon_balance's ~2500% drift). A rod is
                // thin in EVERY presentation at these scales, so the
                // neutron contributes POLE_HALF_THICKNESS·r flat.
                // Inter-alpha same-group distances never reach contact
                // range, so the earned RigidAlpha tables are unaffected.
                let ri = pi_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let rj = pj_prof.radius.max(MIN_RENDER_RADIUS as f64);
                let r_contact = if same_group {
                    let sin_theta_i = (1.0 - cos_theta_i * cos_theta_i).max(0.0).sqrt();
                    let sin_theta_j = (1.0 - cos_theta_j * cos_theta_j).max(0.0).sqrt();
                    let ri_eff = if pi_prof.name == "neutron" {
                        ri * POLE_HALF_THICKNESS
                    } else {
                        ri * (POLE_HALF_THICKNESS + (1.0 - POLE_HALF_THICKNESS) * sin_theta_i)
                    };
                    let rj_eff = if pj_prof.name == "neutron" {
                        rj * POLE_HALF_THICKNESS
                    } else {
                        rj * (POLE_HALF_THICKNESS + (1.0 - POLE_HALF_THICKNESS) * sin_theta_j)
                    };
                    ri_eff + rj_eff
                } else {
                    ri + rj
                };
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

                // Attribute this pair's net force injection to the group.
                // Non-zero only for the un-paired internal terms; `step()`
                // removes the accumulated total mass-weighted (translation
                // only). `same_group` ⇒ both members share group `gi`.
                if same_group {
                    let dp = (self.particles[i].force_accum - f_snap_i)
                        + (self.particles[j].force_accum - f_snap_j);
                    self.group_self_force[gi.unwrap()] += dp;
                }
            }
        }

        // Ambient environmental forces
        for p in &mut self.particles {
            let m = self.profiles[p.profile_id].mass;
            p.force_accum += self.ambient_gravity * m;
            p.force_accum += self.ambient_charge * m;
        }

        // Ambient surface confinement (session-32 Phase C1 — see
        // charge_flow::apply_ambient_confinement). No-op while
        // `ambient_confine` is 0 or in RigidLock.
        self.apply_ambient_confinement();

        // Directional-ambient alignment torque (session-35 — see
        // charge_flow::apply_ambient_alignment). Rotates each nucleus's
        // principal axis toward `ambient_charge_dir`; no-op while
        // `ambient_align` is 0 or the field is isotropic. Zero net force, so
        // it does not feed the self-force accumulator below.
        self.apply_ambient_alignment();

        // Cancel each nucleus's intra-group self-propulsion (session-35; see
        // `group_self_force`). The proton→neutron emission asymmetry (audit:
        // ~1.98 units along carbon's plug axis, un-cancelled across the two
        // like-oriented plug pairs) and the reaction-less vortex leave a net
        // force on an ISOLATED nucleus — physically impossible; cc.pdf says
        // the open charge field radiates that momentum away, isotropically on
        // average, so it cannot translate the COM. Distribute −self_force
        // mass-weighted across the members: a uniform per-mass acceleration
        // has zero torque about the COM, so this removes ONLY the common-mode
        // drift and preserves every internal (inter-alpha carousel, shape)
        // force. No-op in RigidLock (same-group pairs are skipped, so
        // `self_force` is zero) and for pairs of identical particles (their
        // emission is symmetric). Ambient GLOBALS (gravity/charge) are added
        // above and are NOT touched — a real external field still moves the
        // nucleus.
        for gi in 0..self.groups.len() {
            let f = self.group_self_force[gi];
            let gmass = self.groups[gi].mass;
            if f.length_squared() < 1e-24 || gmass <= 0.0 {
                continue;
            }
            let members = self.groups[gi].members.clone();
            for m in members {
                let pm = self.profiles[self.particles[m].profile_id].mass;
                self.particles[m].force_accum -= f * (pm / gmass);
            }
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

    /// Field sources for a subset of a group's members, group-local.
    fn segment_sources(&self, group_idx: usize, members: Option<&[usize]>) -> Vec<FieldSrc> {
        let g = &self.groups[group_idx];
        let pick: Vec<usize> = match members {
            Some(m) => m.to_vec(),
            None => (0..g.members.len()).collect(),
        };
        pick.iter()
            .map(|&k| {
                let profile = self.particles[g.members[k]].profile_id;
                FieldSrc {
                    pos: g.local_offsets[k],
                    pole: g.local_orients[k] * DVec3::Y,
                    mass: self.profiles[profile].mass,
                    radius: self.profiles[profile].radius.max(MIN_RENDER_RADIUS as f64),
                    profile,
                }
            })
            .collect()
    }

    /// March a set of member fields outward from `center` to find the
    /// reach surface (azimuth-averaged reach per polar angle, never
    /// dipping inside the members' silhouette).
    fn march_reach(&self, srcs: &[FieldSrc], center: DVec3, rings: usize) -> Vec<f64> {
        let cq = self.couplings.c_q;
        let push = |x: DVec3| -> f64 {
            let mut f = DVec3::ZERO;
            for s in srcs {
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
                    .map(|s| (s.pos - center).dot(dir) + s.radius)
                    .fold(0.3f64, f64::max);
                let mut r = hull;
                while r < R_MAX && push(center + dir * r) >= F_REF {
                    r += DR;
                }
                sum += r;
            }
            reach.push(sum / AZ_SAMPLES as f64);
        }
        reach
    }

    /// Whole-group azimuth-averaged reach from the group center — feeds
    /// the flow-parcel exits and the disc-exit ring solve. (The SKIN is
    /// drawn per segment instead; see build_group_skin_mesh.)
    pub fn compute_group_reach(&self, group_idx: usize, rings: usize) -> Vec<f64> {
        if self.groups.get(group_idx).is_none() {
            return Vec::new();
        }
        let srcs = self.segment_sources(group_idx, None);
        self.march_reach(&srcs, DVec3::ZERO, rings)
    }

    /// Per-piece skin BANDS: every fused piece draws a THICK DISC clipped
    /// to its own axial extent — no funnels or outer latitudes where a
    /// neighbor continues the structure (those buried the interior and a
    /// total-field lathe blurs into a blob, since a neighboring disc's
    /// 1/r⁴ tail fills every waist). Only the OUTERMOST pieces close with
    /// their own domes; the carousel level flares into the wide
    /// equatorial disc. Each band's radius is that piece's OWN radial
    /// reach (az-averaged, F = 1), so each disc/alpha can be picked out
    /// and its effect followed (session-29 rounds 6–8). GROUP-LOCAL
    /// units, gold, same packing as `build_profile_mesh` (disjoint bands
    /// in one surface).
    pub fn build_group_skin_mesh(
        &self,
        group_idx: usize,
        lon_segments: usize,
        _lat_segments: usize,
    ) -> Vec<f32> {
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        let lon = lon_segments.max(8);
        let verts_per_ring = lon + 1;
        const N_IN: usize = 7; // interior rows per band
        const DOME_STEP: f64 = 0.45;
        const DOME_MAX: usize = 14;
        const AZ: usize = 10;

        let segs: Vec<(Vec<usize>, bool)> = if g.skin_segments.is_empty() {
            vec![((0..g.members.len()).collect(), false)]
        } else {
            g.skin_segments
                .iter()
                .map(|s| (s.members.clone(), s.carousel))
                .collect()
        };

        struct Piece {
            srcs: Vec<FieldSrc>,
            lo: f64,
            hi: f64,
            carousel: bool,
            y_center: f64,
        }
        let mut pieces: Vec<Piece> = Vec::new();
        for (members, carousel) in &segs {
            let srcs = self.segment_sources(group_idx, Some(members));
            if srcs.is_empty() {
                continue;
            }
            let (mut lo, mut hi, mut yc) = (f64::MAX, f64::MIN, 0.0);
            for s in &srcs {
                lo = lo.min(s.pos.y - s.radius - 0.3);
                hi = hi.max(s.pos.y + s.radius + 0.3);
                yc += s.pos.y;
            }
            yc /= srcs.len() as f64;
            pieces.push(Piece {
                srcs,
                lo,
                hi,
                carousel: *carousel,
                y_center: yc,
            });
        }
        if pieces.is_empty() {
            return Vec::new();
        }

        let cq = self.couplings.c_q;
        // Az-averaged radial reach of ONE piece's field at height y.
        let radial_reach = |srcs: &[FieldSrc], y: f64| -> f64 {
            let origin = DVec3::new(0.0, y, 0.0);
            let mut sum = 0.0;
            for az in 0..AZ {
                let phi = TAU * az as f64 / AZ as f64;
                let dir = DVec3::new(phi.cos(), 0.0, phi.sin());
                let hull = srcs
                    .iter()
                    .map(|s| (s.pos - origin).dot(dir) + s.radius)
                    .fold(0.25f64, f64::max);
                let mut r = hull;
                'march: while r < 40.0 {
                    let x = origin + dir * r;
                    let mut f = DVec3::ZERO;
                    for s in srcs {
                        let d_vec = x - s.pos;
                        let rr = d_vec.length().max(SOFTENING);
                        let d_hat = d_vec / rr;
                        let e =
                            self.profiles[s.profile].emission.sample(s.pole.dot(d_hat));
                        f += d_hat * (cq * s.mass * e / (rr * rr * rr * rr));
                    }
                    if f.length() < 1.0 {
                        break 'march;
                    }
                    r += 0.08;
                }
                sum += r;
            }
            sum / AZ as f64
        };

        // Which axial pieces are outermost (get domes)?
        let mut top_piece = None;
        let mut bot_piece = None;
        for (i, p) in pieces.iter().enumerate() {
            if p.carousel {
                continue;
            }
            if top_piece.map_or(true, |t: usize| p.y_center > pieces[t].y_center) {
                top_piece = Some(i);
            }
            if bot_piece.map_or(true, |b: usize| p.y_center < pieces[b].y_center) {
                bot_piece = Some(i);
            }
        }

        // Collect (y, radius) rows per band, top-to-bottom per band.
        let mut band_rows: Vec<Vec<(f64, f64)>> = Vec::new();
        for (i, p) in pieces.iter().enumerate() {
            let mut rows: Vec<(f64, f64)> = Vec::new();
            if !p.carousel && top_piece == Some(i) {
                // Dome above: rows ascend until the piece's own reach
                // collapses, then a small apex row closes the shell.
                let mut ext: Vec<(f64, f64)> = Vec::new();
                let mut y = p.hi + DOME_STEP;
                while ext.len() < DOME_MAX {
                    let r = radial_reach(&p.srcs, y);
                    if r <= 1.0 {
                        break;
                    }
                    ext.push((y, r));
                    y += DOME_STEP;
                }
                ext.push((y, 0.3)); // apex
                ext.reverse();
                rows.extend(ext);
            }
            for k in 0..N_IN {
                let y = p.hi + (p.lo - p.hi) * k as f64 / (N_IN - 1) as f64;
                rows.push((y, radial_reach(&p.srcs, y)));
            }
            if !p.carousel && bot_piece == Some(i) {
                let mut y = p.lo - DOME_STEP;
                let mut n = 0;
                while n < DOME_MAX {
                    let r = radial_reach(&p.srcs, y);
                    if r <= 1.0 {
                        break;
                    }
                    rows.push((y, r));
                    y -= DOME_STEP;
                    n += 1;
                }
                rows.push((y, 0.3)); // apex
            }
            band_rows.push(rows);
        }

        let r_max = band_rows
            .iter()
            .flatten()
            .fold(1e-6f64, |m, &(_, r)| r.max(m));

        const GOLD: (f32, f32, f32) = (0.95, 0.80, 0.50);
        let mut verts: Vec<f32> = Vec::new();
        let mut indices: Vec<f32> = Vec::new();
        let mut vert_base = 0usize;
        for rows in &band_rows {
            if rows.len() < 2 {
                continue;
            }
            for &(y, rho) in rows {
                let strength = ((rho / r_max).powi(2)).clamp(0.0, 1.0) as f32;
                let alpha = 0.03 + 0.14 * strength;
                for lon_idx in 0..=lon {
                    let phi = TAU * lon_idx as f64 / lon as f64;
                    let (nx, nz) = (phi.cos(), phi.sin());
                    verts.extend_from_slice(&[
                        (rho * nx) as f32,
                        y as f32,
                        (rho * nz) as f32,
                    ]);
                    verts.extend_from_slice(&[nx as f32, 0.0, nz as f32]);
                    verts.extend_from_slice(&[GOLD.0, GOLD.1, GOLD.2, alpha]);
                }
            }
            for row in 0..(rows.len() - 1) {
                for lon_idx in 0..lon {
                    let tl = (vert_base + row * verts_per_ring + lon_idx) as f32;
                    let tr = tl + 1.0;
                    let bl = (vert_base + (row + 1) * verts_per_ring + lon_idx) as f32;
                    let br = bl + 1.0;
                    indices.extend_from_slice(&[tl, bl, tr, tr, bl, br]);
                }
            }
            vert_base += rows.len() * verts_per_ring;
        }

        let mut buf: Vec<f32> = Vec::with_capacity(2 + verts.len() + indices.len());
        buf.push(vert_base as f32);
        buf.push(indices.len() as f32);
        buf.extend_from_slice(&verts);
        buf.extend_from_slice(&indices);
        buf
    }

    /// Carousel dispersal overlay: one circle per carousel-level ALPHA
    /// (session-31 addendum A6 — iterate alphas directly instead of
    /// azimuth-clustering member offsets), plane ⊥ the alpha's radial
    /// axis, radius = how far that alpha's field reaches perpendicular to
    /// the main disc — "4 protruding circles following the outside parts".
    /// Group-local at carousel phase 0; render in a child node rotated by
    /// get_group_carousel_phase. Clamped to `0.45·CAROUSEL_R·√2` so
    /// adjacent carousel circles (centers `CAROUSEL_R·√2` apart on a
    /// 4-fold ring) kiss but never interpenetrate (§1.7). Packed like
    /// emission rings: [ring_count, pts_per_ring, xyz…].
    pub fn build_group_carousel_overlay(&self, group_idx: usize) -> Vec<f32> {
        const PTS: usize = 49;
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        let units: Vec<&Vec<usize>> = g
            .alphas
            .iter()
            .filter(|a| {
                let lateral =
                    (a.rest_center.x * a.rest_center.x + a.rest_center.z * a.rest_center.z)
                        .sqrt();
                lateral > 1.5
            })
            .map(|a| &a.members)
            .collect();
        if units.is_empty() {
            return Vec::new();
        }
        let max_r = 0.45 * CAROUSEL_R * std::f64::consts::SQRT_2;
        let mut buf: Vec<f32> = vec![units.len() as f32, PTS as f32];
        for unit in &units {
            let center = unit
                .iter()
                .map(|&k| g.local_offsets[k])
                .sum::<DVec3>()
                / unit.len() as f64;
            let u = DVec3::new(center.x, 0.0, center.z).normalize_or_zero();
            if u.length_squared() < 0.5 {
                continue;
            }
            let w = u.cross(DVec3::Y).normalize();
            // Perpendicular dispersal reach: march this unit's field in
            // the plane ⊥ its axis, average 4 directions.
            let srcs = self.segment_sources(group_idx, Some(unit.as_slice()));
            let mut r_disp = 0.0;
            for dir in [DVec3::Y, -DVec3::Y, w, -w] {
                let hull = srcs
                    .iter()
                    .map(|s| (s.pos - center).dot(dir) + s.radius)
                    .fold(0.3f64, f64::max);
                let cq = self.couplings.c_q;
                let mut r = hull;
                while r < 40.0 {
                    let x = center + dir * r;
                    let mut f = DVec3::ZERO;
                    for s in &srcs {
                        let d_vec = x - s.pos;
                        let rr = d_vec.length().max(SOFTENING);
                        let d_hat = d_vec / rr;
                        let e =
                            self.profiles[s.profile].emission.sample(s.pole.dot(d_hat));
                        f += d_hat * (cq * s.mass * e / (rr * rr * rr * rr));
                    }
                    if f.length() < 1.0 {
                        break;
                    }
                    r += 0.1;
                }
                r_disp += r;
            }
            r_disp = (r_disp / 4.0).min(max_r);
            for p in 0..PTS {
                let a = TAU * p as f64 / (PTS - 1) as f64;
                let pt = center + (DVec3::Y * a.cos() + w * a.sin()) * r_disp;
                buf.extend_from_slice(&[pt.x as f32, pt.y as f32, pt.z as f32]);
            }
        }
        if buf.len() <= 2 {
            return Vec::new();
        }
        buf
    }

    /// Where each AXIAL alpha's disc meets the composite reach surface:
    /// (ring radius, y) per alpha, group-local. Computed once at spawn and
    /// stored as `RigidGroup::disc_exits` — feeds the ring overlay and the
    /// flow-VFX disc exit paths. Re-keyed off `AlphaUnit` directly
    /// (session-31 addendum A5): one disc center per alpha whose
    /// `rest_axis` is axial (|rest_axis·Y| > 0.9) AND has ≥ 2 proton
    /// members, at `y = rest_center.y` — same output as the old
    /// sorted-consecutive-proton-pairing for every preset, but without the
    /// fragile pairing (adjacent-block protons used to sit closer together
    /// than an alpha's own pair, so a gap-based clustering couldn't work).
    pub fn compute_disc_exits(&self, group_idx: usize) -> Vec<(f64, f64)> {
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };
        if g.skin_reach.len() < 2 {
            return Vec::new();
        }
        let mut centers: Vec<f64> = Vec::new();
        for a in &g.alphas {
            if a.rest_axis.dot(DVec3::Y).abs() <= 0.9 {
                continue;
            }
            let proton_count = a
                .members
                .iter()
                .filter(|&&k| {
                    self.profiles[self.particles[g.members[k]].profile_id].name == "proton"
                })
                .count();
            if proton_count >= 2 {
                centers.push(a.rest_center.y);
            }
        }
        if centers.is_empty() {
            return Vec::new();
        }
        centers.sort_by(|a, b| a.partial_cmp(b).unwrap());

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

    /// Three-ring alpha skin (session-32, user-proposed): per standard
    /// alpha, two DISC rings (one at each proton, in that proton's own
    /// live disc plane — tilt readout) and one CENTER ring at the alpha
    /// midplane (equatorial reach level, where the neutron posts live)
    /// with a short radial TICK at each post's actual azimuth (live roll
    /// readout). All geometry is anchored to CURRENT member particle
    /// positions/orientations, so it tracks the bodies in every
    /// `NucleusDynamics` mode — RigidLock kinematics, RigidAlpha rigid
    /// bodies, even FreeNucleon distortion (the "rings" then deform with
    /// the alpha, which is itself a readout). Radii are the spawn-time
    /// quantitative reach values (`compute_alpha_ring_radii`).
    ///
    /// WORLD-space variable-polyline buffer:
    /// `[n_polylines, (pt_count, x,y,z × pt_count)…]`.
    pub fn build_group_alpha_rings(&self, group_idx: usize) -> Vec<f32> {
        const PTS: usize = 49;
        const TICK_IN: f64 = 0.85;
        const TICK_OUT: f64 = 1.15;
        let g = match self.groups.get(group_idx) {
            Some(g) => g,
            None => return Vec::new(),
        };

        let mut buf: Vec<f32> = vec![0.0]; // n_polylines patched at the end
        let mut count = 0u32;
        let push_polyline = |buf: &mut Vec<f32>, pts: &[DVec3]| {
            buf.push(pts.len() as f32);
            for p in pts {
                buf.extend_from_slice(&[p.x as f32, p.y as f32, p.z as f32]);
            }
        };
        let circle = |center: DVec3, normal: DVec3, r: f64| -> Vec<DVec3> {
            let e1 = normal.any_orthonormal_vector();
            let e2 = normal.cross(e1);
            (0..PTS)
                .map(|p| {
                    let phi = TAU * p as f64 / (PTS - 1) as f64;
                    center + (e1 * phi.cos() + e2 * phi.sin()) * r
                })
                .collect()
        };

        for a in &g.alphas {
            if a.ring_disc_r <= 0.0 {
                continue;
            }
            let mut protons: Vec<usize> = Vec::new();
            let mut neutrons: Vec<usize> = Vec::new();
            for &k in &a.members {
                let id = g.members[k];
                match self.profiles[self.particles[id].profile_id].name.as_str() {
                    "proton" => protons.push(id),
                    "neutron" => neutrons.push(id),
                    _ => {}
                }
            }
            if protons.len() != 2 {
                continue;
            }
            let p1 = self.particles[protons[0]].position;
            let p2 = self.particles[protons[1]].position;
            let axis_vec = p2 - p1;
            if axis_vec.length_squared() < 1e-12 {
                continue;
            }
            let axis = axis_vec.normalize();
            let center = (p1 + p2) * 0.5;

            // Disc rings in each proton's own live disc plane.
            for &pid in &protons {
                let pole = self.particles[pid].pole_axis();
                push_polyline(
                    &mut buf,
                    &circle(self.particles[pid].position, pole, a.ring_disc_r),
                );
                count += 1;
            }
            // Center ring at the midplane, ⊥ the live alpha axis.
            push_polyline(&mut buf, &circle(center, axis, a.ring_mid_r));
            count += 1;
            // Post ticks: radial dashes at each neutron's actual azimuth.
            for &nid in &neutrons {
                let d = self.particles[nid].position - center;
                let lat = d - axis * d.dot(axis);
                if lat.length_squared() > 1e-12 {
                    let u = lat.normalize();
                    push_polyline(
                        &mut buf,
                        &[
                            center + u * (a.ring_mid_r * TICK_IN),
                            center + u * (a.ring_mid_r * TICK_OUT),
                        ],
                    );
                    count += 1;
                }
            }
        }
        if count == 0 {
            return Vec::new();
        }
        buf[0] = count as f32;
        buf
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

    /// Anchor frame for a group's flow VFX and group-anchored motes. In
    /// RigidLock this is the phase-driven `(g.com, g.orientation)`. In
    /// RigidAlpha/FreeNucleon the group container is FROZEN while the real
    /// bodies move (addendum A10) — the session-31 round-4 gate hid the
    /// nuclei flow engine there because it visually disconnected. Session
    /// 32 replaces the gate with a LIVE frame derived from member particle
    /// state: center = mean member position, axis = the line through the
    /// two rest-axially-extreme members, applied as a minimal rotation on
    /// top of the frozen orientation — continuous, equal to the frozen
    /// frame while undeformed, and valid in every mode (same
    /// particle-anchored principle as `build_group_alpha_rings`).
    pub(crate) fn live_group_frame(&self, gi: usize) -> (DVec3, DQuat) {
        let g = &self.groups[gi];
        if self.dynamics == NucleusDynamics::RigidLock || g.members.len() < 2 {
            return (g.com, g.orientation);
        }
        let mut com = DVec3::ZERO;
        for &m in &g.members {
            com += self.particles[m].position;
        }
        com /= g.members.len() as f64;
        let (mut k_lo, mut k_hi) = (0usize, 0usize);
        for (k, off) in g.local_offsets.iter().enumerate() {
            if off.y < g.local_offsets[k_lo].y {
                k_lo = k;
            }
            if off.y > g.local_offsets[k_hi].y {
                k_hi = k;
            }
        }
        if k_lo == k_hi {
            return (com, g.orientation);
        }
        let live_axis = self.particles[g.members[k_hi]].position
            - self.particles[g.members[k_lo]].position;
        if live_axis.length_squared() < 1e-12 {
            return (com, g.orientation);
        }
        let frozen_axis = g.orientation * DVec3::Y;
        let align = DQuat::from_rotation_arc(frozen_axis, live_axis.normalize());
        (com, (align * g.orientation).normalize())
    }

    /// Advance the cloud particles by wall-clock `delta` seconds and return
    /// the MultiMesh buffer (12 transform + 4 color floats per particle).
    pub fn advance_clouds(&mut self, delta: f64) -> Vec<f32> {
        const MAX_POOL: usize = 40_000;
        /// Parcel speed along its path (natural units / wall second).
        const FLOW_SPEED: f64 = 2.4;
        /// Recycling loops per second through a free baryon (split across
        /// the two poles). Electrons recycle a sliver of this.
        const FLOW_FREE_PER_SEC: f64 = 800.0;
        /// Loops per second per alpha engine in a nucleus — more engines
        /// in the stack, denser flow.
        const FLOW_GROUP_PER_ENGINE: f64 = 380.0;
        /// A cap electron only DAMPS the pole flow it rides (it skims what
        /// it can carry — its stream is 1/1836 of the proton's), it does
        /// not stopper it: effective flow occlusion is scaled by this.
        const ELECTRON_FLOW_DAMP: f64 = 0.35;
        /// Each plug nucleon on a stack end boosts that end's pole-flow
        /// density: its perpendicular disc (proton) or axial channel
        /// (neutron) gathers extra field into the socket (phos.pdf
        /// plug-and-socket; atmo2.pdf "pulling charge into the axial
        /// holes"). Rate multiplier = 1 + BOOST × plugs.
        const PLUG_FLOW_BOOST: f64 = 0.5;
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
        // MOTE_FRACTION/MOTE_SCALE/MOTE_LIFETIME are runtime knobs now
        // (session-31 round 4, work item 5) — `self.mote_fraction`,
        // `self.mote_scale`, `self.mote_lifetime` fields, seeded from
        // `DEFAULT_MOTE_FRACTION`/`DEFAULT_MOTE_SCALE`/`DEFAULT_MOTE_LIFETIME`
        // in `AtomCore::new`, exposed to the HUD via `AtomSim`
        // get/set_mote_fraction/scale/lifetime for a follow-up agent to wire
        // up debug sliders.
        /// The skin gold — same tint the retired lathe surface used.
        const MOTE_COLOR: (f32, f32, f32) = (0.95, 0.80, 0.50);

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
        let group_frames: Vec<(DVec3, DQuat)> = (0..self.groups.len())
            .map(|gi| self.live_group_frame(gi))
            .collect();

        // Age + kinematic path evaluation. A skin-dying Flow parcel that
        // crosses death this frame has a chance to convert IN PLACE into a
        // SkinMote instead of being culled — flow parcels already die
        // exactly at the field-reach boundary, so letting a fraction of
        // those deaths leave a bright dot sketches the boundary as a
        // living point-cloud (session-31: replaces the lathe skin, which
        // read as a "gold sausage" blob around multi-piece nuclei). Flow
        // parcels keep priority over the pool budget, so conversion is
        // skipped once the pool is nearly full.
        let skip_motes = self.vfx_particles.len() as f64 >= MAX_POOL as f64 * 0.9;
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
                VfxKind::SkinMote {
                    ref anchor,
                    local,
                    swirl,
                } => {
                    let frame = match *anchor {
                        VfxAnchor::Particle(i) => emitters.get(i).map(|e| (e.pos, e.rot)),
                        VfxAnchor::Group(gi) => group_frames.get(gi).copied(),
                    };
                    let Some((origin, rot)) = frame else {
                        vp.age = vp.lifetime + 1.0; // anchor gone — cull
                        continue;
                    };
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

            // Read-only probe for a skin-death conversion — kept as its
            // own statement so the borrow from the match above (through
            // `pts`/`anchor`) is fully released before we write vp.kind
            // back below.
            let mote: Option<(VfxAnchor, DVec3, f64)> = if let VfxKind::Flow {
                anchor,
                swirl,
                dies_at_skin,
                pts,
                ..
            } = &vp.kind
            {
                if *dies_at_skin && vp.age >= vp.lifetime && !skip_motes {
                    Some((*anchor, *pts.last().unwrap(), *swirl))
                } else {
                    None
                }
            } else {
                None
            };
            if let Some((anchor, last, swirl)) = mote {
                if self.rng.next_f64() < self.mote_fraction {
                    // Swirl-rotated final point, evaluated at u=1.0 (i.e.
                    // at the parcel's nominal lifetime, not its possibly
                    // slightly-overshot age) — the mote's fixed anchor-local
                    // position.
                    let a_end = swirl * vp.lifetime;
                    let (se, ce) = a_end.sin_cos();
                    let mote_local = DVec3::new(
                        last.x * ce + last.z * se,
                        last.y,
                        -last.x * se + last.z * ce,
                    );
                    vp.color = MOTE_COLOR;
                    vp.base_scale *= self.mote_scale;
                    vp.grow = 0.0;
                    vp.age = 0.0;
                    vp.lifetime = self.mote_lifetime * (0.8 + 0.4 * self.rng.next_f64());
                    vp.kind = VfxKind::SkinMote {
                        anchor,
                        local: mote_local,
                        swirl: swirl * 0.25,
                    };
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
                // An electron cap only skims the stream; a baryon stopper
                // really blocks it.
                let occ_flow = if blocker.is_some() {
                    occ_v * ELECTRON_FLOW_DAMP
                } else {
                    occ_v
                };
                let n = spawn_count(
                    &mut self.rng,
                    FLOW_FREE_PER_SEC * 0.5 * (1.0 - occ_flow) * rate_scale,
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
                            dies_at_skin: true, // free recycling loop — ends at the reach boundary
                        },
                    });
                }
                // Electron-capped pole: the skimmed share of the inflow is
                // CAPTURED at the rider and dispersed off its disc — the
                // rest recycles through as usual.
                if blocker.is_some() && occ_v > 0.5 {
                    let n = spawn_count(
                        &mut self.rng,
                        FLOW_FREE_PER_SEC * 0.5 * occ_v * ELECTRON_FLOW_DAMP * rate_scale,
                    );
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
                                dies_at_skin: false, // captured at the rider, not the skin
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
        //
        // Session-31 round 4 gated this engine to RigidLock because it
        // rode the group's FROZEN phase-driven frame and visually
        // disconnected in RigidAlpha/FreeNucleon (user report). Session 32
        // removed the gate: `group_frames` now comes from
        // `live_group_frame`, which tracks the actual member bodies in
        // every mode, so the flow (and the motes it seeds) stays attached
        // to the stack it annotates. Rest-frame path data (skin_reach,
        // disc_exits, local y-extents) rides that live frame — exact while
        // the stack is undeformed, and a few-percent approximation under
        // bound-stack breathing, which is well inside the VFX bar.
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
            // Plug nucleons per end: anything sitting beyond the outermost
            // AXIAL proton is a polar plug — its disc/channel gathers extra
            // field into the socket, boosting that end's intake.
            let mut y_pmax = 0.0f64;
            for (k, &m) in self.groups[gi].members.iter().enumerate() {
                let axial = (self.groups[gi].local_orients[k] * DVec3::Y)
                    .dot(DVec3::Y)
                    .abs()
                    > 0.9;
                if axial && self.profiles[self.particles[m].profile_id].name == "proton" {
                    y_pmax = y_pmax.max(self.groups[gi].local_offsets[k].y.abs());
                }
            }
            let mut plugs = [0usize; 2]; // [top, bottom]
            if y_pmax > 0.0 {
                for off in &self.groups[gi].local_offsets {
                    if off.y > y_pmax + 0.2 {
                        plugs[0] += 1;
                    } else if off.y < -(y_pmax + 0.2) {
                        plugs[1] += 1;
                    }
                }
            }
            for side in 0..2 {
                let (entry, tip_y, far_tip, reach_out, n_plugs) = if side == 0 {
                    (1.0, y_top, y_bot, reach_s, plugs[0])
                } else {
                    (-1.0, y_bot, y_top, reach_n, plugs[1])
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
                let occ_flow = if blocker.is_some() {
                    occ_v * ELECTRON_FLOW_DAMP
                } else {
                    occ_v
                };
                let boost = 1.0 + PLUG_FLOW_BOOST * n_plugs as f64;
                let end_rate = FLOW_GROUP_PER_ENGINE * engines as f64 * 0.5 * boost;
                let n = spawn_count(&mut self.rng, end_rate * (1.0 - occ_flow));
                for _ in 0..n {
                    if self.vfx_particles.len() >= MAX_POOL {
                        break;
                    }
                    // The boosted share arrives SIDEWAYS through the plug —
                    // its perpendicular disc gathers field into the socket.
                    let side_feed =
                        n_plugs > 0 && self.rng.next_f64() < (boost - 1.0) / boost;
                    let w_total = exits.len() as f64 + THROUGH_WEIGHT;
                    let pick = self.rng.next_f64() * w_total;
                    let (pts, len, color_end) = if pick >= exits.len() as f64 {
                        let (pts, len) = flow_group_through(
                            &mut self.rng, entry, tip_y, far_tip, reach_out, side_feed,
                        );
                        (pts, len, THROUGH_COLOR)
                    } else {
                        let (rr, ry) = exits[pick as usize];
                        let (pts, len) =
                            flow_group_disc(&mut self.rng, entry, tip_y, rr, ry, side_feed);
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
                            dies_at_skin: true, // through-path or disc exit — both end at the reach boundary
                        },
                    });
                }
                // Cap electron on this end: the skimmed share is captured
                // and dispersed; the rest passed through above.
                if blocker.is_some() && occ_v > 0.5 {
                    let n =
                        spawn_count(&mut self.rng, end_rate * occ_v * ELECTRON_FLOW_DAMP);
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
                                dies_at_skin: false, // captured at the rider, not the skin
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
            // dying there. SkinMotes get their own smooth in-out instead
            // (they don't "travel", they just bloom and fade in place).
            let fade = match vp.kind {
                VfxKind::SkinMote { .. } => (std::f32::consts::PI * t).sin().max(0.0),
                _ => (t / 0.12).min(1.0) * (1.0 - t * t * t).max(0.0),
            };
            let scale = (vp.base_scale * (1.0 + vp.grow * vp.age as f32)).max(0.004);
            // Flow parcels recolor as they recycle: intake tint on the way
            // in, exit tint on the way out. SkinMotes hold the skin-gold
            // tint steady — they mark a boundary point, not a transition.
            let (cr, cg, cb, brightness) = match vp.kind {
                VfxKind::Flow { color_end, .. } => (
                    vp.color.0 + (color_end.0 - vp.color.0) * t,
                    vp.color.1 + (color_end.1 - vp.color.1) * t,
                    vp.color.2 + (color_end.2 - vp.color.2) * t,
                    0.9f32, // tiny dots need the extra alpha
                ),
                VfxKind::Smoke => (vp.color.0, vp.color.1, vp.color.2, 0.6),
                VfxKind::SkinMote { .. } => (vp.color.0, vp.color.1, vp.color.2, 0.9f32),
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
    // Impact parameter: the exit angle picks how deep through the body the
    // parcel bores — polar exits pass through the middle, equatorial exits
    // ride the outer bore and deflect. Flow spans the WHOLE diameter, not
    // a few-pixel core.
    let b = 0.8 * radius * (theta_out / FRAC_PI_2).clamp(0.0, 1.0).powf(0.7);
    let phi_b = phi0 + spin * 2.2;
    for k in 0..4 {
        let f = k as f64 / 3.0;
        let h = top - (top - radius) * f;
        let taper = ((h / radius - 1.0) / 1.6).clamp(0.0, 1.0);
        let rho = b + (mouth - b).max(0.0) * taper.powf(0.65) + 0.03 * radius;
        let phi = phi0 + spin * 2.2 * f;
        raw.push(DVec3::new(rho * phi.cos(), entry * h, rho * phi.sin()));
    }
    raw.push(DVec3::new(
        b * phi_b.cos(),
        entry * 0.45 * radius,
        b * phi_b.sin(),
    ));
    raw.push(DVec3::new(
        b * phi_b.cos(),
        -entry * 0.55 * radius,
        b * phi_b.sin(),
    ));
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
    let b = rng.next_f64() * 0.45 * radius;
    for k in 0..3 {
        let f = k as f64 / 2.0;
        let h = tip_y + entry * radius * (1.9 - 1.5 * f);
        let rho = b + (mouth - b).max(0.0) * (1.0 - f) + 0.05 * radius;
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

/// Entry prefix for a nucleus flow: either the polar funnel above the
/// stack end, or — when a plug gathers the perpendicular field — a
/// SIDEWAYS approach spiraling into the plug at the end's latitude.
fn push_group_entry(
    raw: &mut Vec<DVec3>,
    rng: &mut Rng,
    entry: f64,
    tip_y: f64,
    b: f64,
    side_feed: bool,
) -> f64 {
    let phi0 = rng.next_f64() * TAU;
    if side_feed {
        for k in 0..3 {
            let f = k as f64 / 2.0;
            let rho = b + (2.9 - b) * (1.0 - f);
            let h = tip_y + entry * (0.2 - 0.5 * f);
            let phi = phi0 + 1.6 * f;
            raw.push(DVec3::new(rho * phi.cos(), h, rho * phi.sin()));
        }
    } else {
        let mouth = 0.6 + 0.9 * rng.next_f64();
        for k in 0..3 {
            let f = k as f64 / 2.0;
            let h = tip_y + entry * (1.9 - 1.7 * f);
            let rho = b + (mouth - b).max(0.0) * (1.0 - f) + 0.04;
            let phi = phi0 + 2.0 * f;
            raw.push(DVec3::new(rho * phi.cos(), h, rho * phi.sin()));
        }
    }
    phi0
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
    side_feed: bool,
) -> (Vec<DVec3>, f64) {
    let mut raw: Vec<DVec3> = Vec::with_capacity(8);
    // Bore radius inside the stack hole — the channel fills the whole
    // hole, not a hairline on the axis.
    let b = rng.next_f64().sqrt() * 0.55;
    let phi0 = push_group_entry(&mut raw, rng, entry, tip_y, b, side_feed);
    let phi_b = phi0 + 2.4;
    raw.push(DVec3::new(
        b * phi_b.cos(),
        (tip_y + far_tip_y) * 0.5,
        b * phi_b.sin(),
    ));
    raw.push(DVec3::new(b * phi_b.cos(), far_tip_y, b * phi_b.sin()));
    let out_r = reach_out.max(far_tip_y.abs() + 1.2);
    let phi_e = phi_b + 0.8;
    raw.push(DVec3::new(
        (b * 0.7 + 0.1) * phi_e.cos(),
        far_tip_y - entry * 0.6,
        (b * 0.7 + 0.1) * phi_e.sin(),
    ));
    raw.push(DVec3::new(
        (b * 0.7 + 0.15) * (phi_e + 0.7).cos(),
        -entry * out_r,
        (b * 0.7 + 0.15) * (phi_e + 0.7).sin(),
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
    side_feed: bool,
) -> (Vec<DVec3>, f64) {
    let mut raw: Vec<DVec3> = Vec::with_capacity(9);
    let b = rng.next_f64().sqrt() * 0.55;
    let phi0 = push_group_entry(&mut raw, rng, entry, tip_y, b, side_feed);
    let phi_b = phi0 + 2.4;
    raw.push(DVec3::new(
        b * phi_b.cos(),
        (tip_y + ring_y) * 0.5,
        b * phi_b.sin(),
    ));
    raw.push(DVec3::new(b * phi_b.cos(), ring_y, b * phi_b.sin()));
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

/// Render orientation. Group members: just `p.orientation` — the roll
/// (each alpha spinning about its own axis) and the carousel orbit are
/// already baked into it by `sync_group_members` (§1.3). NO swing–twist
/// extraction and NO extra rotY on top: that decomposition is singular
/// near 180° configurations and made off-axis members flip wildly as the
/// ride swept them through (session-29). Free particles: the
/// parallel-transported display frame. The force model never reads
/// rotation about the pole, so both are purely visual.
fn render_orientation(p: &SimParticle) -> DQuat {
    if p.group.is_some() {
        p.orientation
    } else {
        p.display_orientation
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
    let rq = render_orientation(p);
    let x = rq * DVec3::X;
    let y = rq * DVec3::Y;
    let z = rq * DVec3::Z;

    // Godot's MultiMesh buffer wants the MATRIX ROWS (basis.rows[i] =
    // (x_axis.i, y_axis.i, z_axis.i)) with the origin interleaved — NOT
    // the axis vectors written out as rows. Writing axes-as-rows is the
    // TRANSPOSE, i.e. every body renders with the INVERSE rotation. That
    // bug hid for months because all visible rotations were about each
    // body's own pole (a body of revolution spinning backwards looks
    // identical, and ±pole flips are invisible on symmetric traces); the
    // carousel ride was the first big rotation about a non-pole axis and
    // rendered as a synchronized counter-spin "twirl" (session-31,
    // verified against godot mesh_storage.cpp multimesh_instance_set_transform).
    buf.extend_from_slice(&[
        x.x as f32 * scale, y.x as f32 * scale, z.x as f32 * scale, pos.x as f32,
        x.y as f32 * scale, y.y as f32 * scale, z.y as f32 * scale, pos.y as f32,
        x.z as f32 * scale, y.z as f32 * scale, z.z as f32 * scale, pos.z as f32,
    ]);
    buf.extend_from_slice(&[color.0, color.1, color.2, alpha]);
}

// ── Trace-model reorientation ─────────────────────────────────────────────

/// Reorient a saved spin-mode path trace so the pattern's symmetry axis —
/// the axial hole the body spins around, which in spin mode is the
/// outermost orbital level's LAB axis (X, Y or Z depending on the stack) —
/// lands on +Y, the pole axis atom-mode bodies spin about (an alpha's
/// `roll_phase` rotates group-member meshes about local +Y).
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
    // CENTER the pattern on the spin axis: an off-center trace hula-hoops
    // around the pole when the display spin turns it — that was the
    // "twirling in place" / synchronized-ballet report (session-29).
    points.iter().map(|p| rot * (*p - centroid)).collect()
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

    /// A trace already aligned to Y (and centered) must come back
    /// (near-)unchanged.
    #[test]
    fn trace_reorient_is_stable_for_canonical_input() {
        let mut pts = Vec::new();
        for i in 0..400 {
            let a = i as f64 * 0.377;
            // Zero-mean jitter: the reorient now CENTERS the pattern, so
            // a biased jig would read as an intentional shift.
            let jig = (((i * 7) % 13) as f64 / 12.0 - 0.5) * 0.3;
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

    /// Electron-capped pole: the rider only SKIMS the stream — part of
    /// the inflow is captured and dispersed at the electron, the rest
    /// still recycles through, at a visibly damped rate vs the open pole.
    #[test]
    fn flow_damped_and_skimmed_at_electron_capped_pole() {
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
        // Captures end high at the rider (last.y > 0.5); loops entering
        // north exit on the far hemisphere (last.y < 0).
        let (mut captures, mut north_loops, mut south_in) = (0usize, 0usize, 0usize);
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
                    if last.y > 0.5 {
                        captures += 1;
                    } else if last.y < -0.1 {
                        north_loops += 1;
                    }
                } else {
                    south_in += 1;
                }
            }
        }
        assert!(captures > 0, "rider should skim and disperse part of the inflow");
        assert!(
            north_loops > 0,
            "damped flow must still pass the electron, not be destroyed"
        );
        assert!(
            north_loops < south_in,
            "capped pole should be visibly damped: north={north_loops} south={south_in}"
        );
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
                } else if last.y.abs() > 6.2 {
                    through += 1;
                }
            }
        }
        assert!(disc > 0, "expected disc-ring exits, got none");
        assert!(through > 0, "expected through-channel exits, got none");
    }

    /// Session-31: the lathe "skin" surface is retired in favor of
    /// boundary motes — skin-dying flow parcels (free loops, group
    /// through-paths, group disc exits) leave a bright skin-gold dot
    /// behind at their death point instead of always being culled. Over a
    /// few seconds of a neon nucleus's flow, the buffer should contain at
    /// least one MOTE_COLOR instance, and the pool must stay bounded.
    #[test]
    fn flow_deaths_leave_skin_motes() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon preset");
        core.set_vfx_enabled(true);

        const MOTE_COLOR: (f32, f32, f32) = (0.95, 0.80, 0.50);
        const MAX_POOL: usize = 40_000; // mirrors advance_clouds' local const

        let mut found_mote = false;
        for _ in 0..120 {
            // 120 × 0.05s = 6.0 simulated seconds.
            let buf = core.advance_clouds(0.05);
            assert!(
                core.vfx_count() < MAX_POOL,
                "vfx pool exceeded its budget: {}",
                core.vfx_count()
            );
            for chunk in buf.chunks_exact(16) {
                let (r, g, b, a) = (chunk[12], chunk[13], chunk[14], chunk[15]);
                if a > 0.0
                    && (r - MOTE_COLOR.0).abs() < 1e-3
                    && (g - MOTE_COLOR.1).abs() < 1e-3
                    && (b - MOTE_COLOR.2).abs() < 1e-3
                {
                    found_mote = true;
                }
            }
        }
        assert!(
            found_mote,
            "expected at least one visible skin-mote instance in the buffer"
        );
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
        assert_eq!(
            rings[0] as usize,
            2,
            "carbon (session-32 shape) = two core alphas, two rings (plugs ring nothing)"
        );

        core.clear_particles();
        core.spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("tri_alpha preset");
        let rings = core.build_group_emission_rings(0);
        assert_eq!(rings[0] as usize, 3, "tri_alpha = three alphas, three rings");

        // Carousel alphas are NOT axial engines: neon rings only its
        // center alpha, argon its five axial disks.
        core.clear_particles();
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon preset");
        let rings = core.build_group_emission_rings(0);
        assert_eq!(rings[0] as usize, 1, "neon = one axial alpha, one ring");

        core.clear_particles();
        core.spawn_preset("argon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("argon preset");
        let rings = core.build_group_emission_rings(0);
        assert_eq!(
            rings[0] as usize,
            3,
            "argon = center + two caps axial (connectors are sideways), three rings"
        );
    }

    /// Three-ring alpha skin (session-32): every standard 4-nucleon alpha
    /// gets spawn-time reach radii and a live world-space polyline set —
    /// 2 disc rings + 1 center ring + 2 post ticks = 5 polylines per
    /// alpha, in EVERY dynamics mode (the whole point: the overlay must
    /// not disappear in RigidAlpha the way the group-frame rings did).
    #[test]
    fn alpha_rings_live_in_all_modes() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        core.running = true;

        // Session-32 carbon: 2 ring-bearing core alphas + 4 single-member
        // plug alphas (which draw no rings — radii stay 0).
        let ringed = core.groups[0]
            .alphas
            .iter()
            .filter(|a| a.ring_disc_r > 0.0)
            .count();
        assert_eq!(ringed, 2, "carbon = two ring-bearing core alphas");
        for a in core.groups[0].alphas.iter().filter(|a| a.ring_disc_r > 0.0) {
            assert!(
                a.ring_disc_r > 1.0 && a.ring_mid_r > 1.0,
                "spawn-time ring radii should be outside the body: disc={} mid={}",
                a.ring_disc_r,
                a.ring_mid_r
            );
            // The reach surface BULGES at the proton disc planes (that's
            // where the alpha pushes farthest) and waists slightly at the
            // midplane between them — the disc rings must sit outside the
            // center ring, or the radii sampling is wrong.
            assert!(
                a.ring_disc_r > a.ring_mid_r,
                "disc ring should bulge past the midplane waist: disc={} mid={}",
                a.ring_disc_r,
                a.ring_mid_r
            );
        }

        let parse_polylines = |buf: &[f32]| -> usize {
            assert!(!buf.is_empty(), "ring buffer should not be empty");
            let n = buf[0] as usize;
            let mut o = 1usize;
            for _ in 0..n {
                let pts = buf[o] as usize;
                assert!(pts >= 2, "each polyline needs >= 2 points");
                o += 1 + pts * 3;
            }
            assert_eq!(o, buf.len(), "buffer length must match its header");
            n
        };

        for mode in [
            NucleusDynamics::RigidLock,
            NucleusDynamics::RigidAlpha,
            NucleusDynamics::FreeNucleon,
        ] {
            core.set_nucleus_dynamics(mode);
            core.step_n(200);
            let buf = core.build_group_alpha_rings(0);
            let n = parse_polylines(&buf);
            assert_eq!(
                n, 10,
                "carbon = 2 core alphas x (2 disc + 1 center + 2 ticks) polylines in {mode:?}"
            );
        }
    }

    /// Display rotations are WALL-clock: physics stepping must not move
    /// the carousel or any alpha's roll phase (that's what made visible
    /// speed scale with the substep count) — advance_display must, at the
    /// documented rates (session-31 addendum A8).
    #[test]
    fn display_phases_are_wall_clock() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        core.running = true;
        let car0 = core.groups[0].carousel_phase;
        let roll0 = core.groups[0].alphas[0].roll_phase;
        core.step_n(2000);
        assert_eq!(
            core.groups[0].carousel_phase, car0,
            "sim steps must not advance the carousel"
        );
        assert_eq!(
            core.groups[0].alphas[0].roll_phase, roll0,
            "sim steps must not advance an alpha's roll"
        );
        core.advance_display(0.5);
        assert!(
            (core.groups[0].carousel_phase - (car0 + CAROUSEL_VIS_RATE * 0.5) % TAU).abs()
                < 1e-9,
            "advance_display should ride the carousel at CAROUSEL_VIS_RATE"
        );
        assert!(
            (core.groups[0].alphas[0].roll_phase - (roll0 + DISPLAY_SPIN_RATE * 0.5) % TAU)
                .abs()
                < 1e-9,
            "advance_display should roll the alpha at DISPLAY_SPIN_RATE"
        );

        // Roll and carousel are INDEPENDENT rates (session-31: locking
        // them equal made the whole nucleus read as one solid gear).
        assert!(
            (CAROUSEL_VIS_RATE - DISPLAY_SPIN_RATE).abs() > 1e-6,
            "carousel orbit rate must be independent of the roll rate"
        );
    }

    /// Skin segmentation: carbon reads as ONE tube; neon = center +
    /// carousel; argon = center, 2 connectors, 2 caps, carousel.
    #[test]
    fn skin_segments_match_structure() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);

        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon");
        assert_eq!(
            core.groups[0].skin_segments.len(),
            1,
            "carbon = ONE contiguous tube: the fused plug pairs seat at \
             PLUG_SEAT (1.4, inside the 1.5 segment gap) — 'that seventh \
             proton is in tight' (ammon.pdf); the session-32 three-segment \
             reading encoded the old funnel-mouth standoff"
        );

        core.clear_particles();
        core.spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("tri_alpha");
        assert_eq!(
            core.groups[0].skin_segments.len(),
            1,
            "bare 3-stack = one tube"
        );

        core.clear_particles();
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon");
        let segs = &core.groups[0].skin_segments;
        assert_eq!(segs.len(), 2, "neon = center + carousel");
        assert_eq!(segs.iter().filter(|s| s.carousel).count(), 1);

        core.clear_particles();
        core.spawn_preset("argon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("argon");
        let segs = &core.groups[0].skin_segments;
        assert_eq!(
            segs.len(),
            6,
            "argon = 2 caps + 2 connectors + center + carousel"
        );
        assert_eq!(segs.iter().filter(|s| s.carousel).count(), 1);
        // Carousel overlay: four dispersal circles riding the ride.
        let overlay = core.build_group_carousel_overlay(0);
        assert!(!overlay.is_empty(), "argon carousel should have an overlay");
        assert_eq!(overlay[0] as usize, 4, "four carousel units, four circles");
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
        let buf = core.build_group_skin_mesh(0, 16, 8);
        let n_verts = buf[0] as usize;
        assert!(n_verts > 0, "alpha skin should have geometry");
        // Generic scan (format-independent): widest radial extent is the
        // equatorial reach; the extreme-|y| vertex is the dome apex.
        let (mut r_eq, mut y_extreme, mut r_at_extreme) = (0.0f64, 0.0f64, f64::MAX);
        for v in 0..n_verts {
            let o = 2 + v * 10;
            let (x, y, z) = (buf[o] as f64, buf[o + 1] as f64, buf[o + 2] as f64);
            let radial = (x * x + z * z).sqrt();
            r_eq = r_eq.max(radial);
            if y.abs() > y_extreme {
                y_extreme = y.abs();
                r_at_extreme = radial;
            }
        }
        let single =
            (core.couplings.c_q * core.profiles[p_id].emission.sample(0.0)).powf(0.25);
        assert!(
            r_eq > single * 1.05,
            "alpha equator reach {r_eq} should exceed a single proton's {single}"
        );
        assert!(r_eq < 12.0, "alpha equator reach implausibly large: {r_eq}");
        assert!(
            r_at_extreme < 1.0,
            "skin must close at the poles (apex radius {r_at_extreme} at |y|={y_extreme})"
        );
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

    /// TEMP diagnostic probe: reproduce the APP loop for the 3-alpha
    /// stack + RigidAlpha (session-31 user report: dissolves within
    /// seconds in the app while the contiguous-step harness holds 100k
    /// steps at 2.1% — the structure was then called "carbon"; session 32
    /// renamed the bare stack `tri_alpha` and gave carbon its real
    /// plugged 2+2 shape). Mimics atom_mode.gd _process: step_n(substeps)
    /// + advance_display(dt) + advance_clouds(dt) with VFX enabled, mode
    /// toggled after a few seconds of RigidLock display time.
    #[test]
    #[ignore]
    fn probe_app_loop_carbon_rigidalpha() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        let e_csv = load_histogram_csv(&dir.join("histogram_electron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.register_profile("electron", 1.0 / 1836.0, 0.3, &e_csv);
        core.spawn_preset("tri_alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("tri_alpha");
        core.set_vfx_enabled(true);
        core.running = true;

        let alpha_centers = |core: &AtomCore| -> Vec<DVec3> {
            let g = &core.groups[0];
            g.alphas
                .iter()
                .map(|a| {
                    let sum: DVec3 = a
                        .members
                        .iter()
                        .map(|&k| core.particles[g.members[k]].position)
                        .sum();
                    sum / a.members.len() as f64
                })
                .collect()
        };
        let spacings = |core: &AtomCore| -> Vec<f64> {
            let c = alpha_centers(core);
            let mut d = Vec::new();
            for i in 0..c.len() {
                for j in (i + 1)..c.len() {
                    d.push((c[i] - c[j]).length());
                }
            }
            d
        };

        const DT_FRAME: f64 = 1.0 / 60.0;
        // 3 s of RigidLock display first (user watches, then toggles).
        for _ in 0..180 {
            core.step_n(100);
            core.advance_display(DT_FRAME);
            core.advance_clouds(DT_FRAME);
        }
        core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
        let d0 = spacings(&core);
        println!("t=0s spacings: {d0:?}");
        // 30 wall-seconds of app loop = 180k steps.
        for sec in 1..=30 {
            for _ in 0..60 {
                core.step_n(100);
                core.advance_display(DT_FRAME);
                core.advance_clouds(DT_FRAME);
            }
            let d = spacings(&core);
            let drift = d
                .iter()
                .zip(&d0)
                .map(|(a, b)| ((a - b) / b).abs())
                .fold(0.0f64, f64::max);
            println!(
                "t={sec:>2}s max_drift={:>7.1}% spacings={:?}",
                drift * 100.0,
                d.iter().map(|x| (x * 100.0).round() / 100.0).collect::<Vec<_>>()
            );
            if drift > 3.0 {
                println!("DISSOLVED — reproduction successful");
                break;
            }
        }
    }

    /// TEMP session-31 diagnostic probe (run with `-- --ignored
    /// probe_carousel_twirl --nocapture`): reproduce the app frame loop
    /// on neon and measure each body's ACTUAL angular velocity,
    /// decomposed about the core axis (world Y) and about its own
    /// alpha's current axis. Expected: carousel members show
    /// ω·Y = CAROUSEL_VIS_RATE (the ride) and ω·u = DISPLAY_SPIN_RATE
    /// (the roll); anything else is the reported twirl.
    #[test]
    #[ignore]
    fn probe_carousel_twirl() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon");
        let gid = 0;
        let proton_id = core.profile_id_by_name("proton").unwrap();

        // One core proton + one proton and one post from each carousel alpha.
        let mut tracked: Vec<(String, usize, usize)> = Vec::new(); // (label, particle, alpha)
        for (ai, a) in core.groups[gid].alphas.iter().enumerate() {
            let mut got_proton = false;
            for &k in &a.members {
                let pid = core.groups[gid].members[k];
                let is_p = core.particles[pid].profile_id == proton_id;
                if is_p && !got_proton {
                    tracked.push((
                        format!("{}[{ai}].proton", if a.orbits_core { "car" } else { "core" }),
                        pid,
                        ai,
                    ));
                    got_proton = true;
                } else if !is_p && a.orbits_core && a.members.len() == 4 {
                    tracked.push((format!("car[{ai}].post"), pid, ai));
                    break;
                }
            }
        }

        // Mirror the app loop: 100 physics substeps + one wall-clock
        // display advance per 60 fps frame.
        const FRAMES: usize = 120;
        const DT: f64 = 1.0 / 60.0;
        let mut prev_q: Vec<DQuat> =
            tracked.iter().map(|t| core.particles[t.1].orientation).collect();
        let mut prev_pos: Vec<DVec3> =
            tracked.iter().map(|t| core.particles[t.1].position).collect();
        let mut acc_w: Vec<DVec3> = vec![DVec3::ZERO; tracked.len()];
        let mut acc_orbit: Vec<f64> = vec![0.0; tracked.len()];
        for _ in 0..FRAMES {
            core.step_n(100);
            core.advance_display(DT);
            for (ti, t) in tracked.iter().enumerate() {
                let q = core.particles[t.1].orientation;
                let dq = (q * prev_q[ti].conjugate()).normalize();
                let (axis, angle) = dq.to_axis_angle();
                // Wrap to the short way around.
                let angle = if angle > std::f64::consts::PI {
                    angle - TAU
                } else {
                    angle
                };
                acc_w[ti] += axis * (angle / DT);
                let p0 = prev_pos[ti];
                let p1 = core.particles[t.1].position;
                acc_orbit[ti] +=
                    (p1.z.atan2(p1.x) - p0.z.atan2(p0.x) + std::f64::consts::PI)
                        .rem_euclid(TAU) - std::f64::consts::PI;
                prev_q[ti] = q;
                prev_pos[ti] = p1;
            }
        }
        let total_t = FRAMES as f64 * DT;
        println!(
            "\nexpected: ride(Y) = {CAROUSEL_VIS_RATE}, roll(u) = {DISPLAY_SPIN_RATE}\n\
             {:<16} {:>8} {:>8} {:>8} {:>10}",
            "body", "w_Y", "w_u", "|w|", "orbit_rate"
        );
        for (ti, t) in tracked.iter().enumerate() {
            let w = acc_w[ti] / FRAMES as f64;
            let a = &core.groups[gid].alphas[t.2];
            let car = if a.orbits_core {
                DQuat::from_rotation_y(core.groups[gid].carousel_phase)
            } else {
                DQuat::IDENTITY
            };
            let u = (core.groups[gid].orientation * (car * a.rest_axis)).normalize();
            println!(
                "{:<16} {:>8.3} {:>8.3} {:>8.3} {:>10.3}",
                t.0,
                w.dot(DVec3::Y),
                w.dot(u),
                w.length(),
                acc_orbit[ti] / total_t
            );
        }
    }

    /// Overlay-clamp (session-31 addendum A6/A8, replaces the old
    /// `carousel_circles_tangent_at_new_radius` which assumed CAROUSEL_R
    /// 8.0): every carousel dispersal circle's radius stays within the
    /// §1.7 clamp `0.45·CAROUSEL_R·√2` at the restored R = 4.5, so
    /// adjacent circles (centers `CAROUSEL_R·√2` apart on the 4-fold
    /// ring) kiss but never interpenetrate.
    #[test]
    fn carousel_overlay_circles_are_clamped() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon");
        assert!((CAROUSEL_R - 4.5).abs() < 1e-9, "CAROUSEL_R should be restored to ~4.5");
        let overlay = core.build_group_carousel_overlay(0);
        assert!(!overlay.is_empty());
        let n_units = overlay[0] as usize;
        let pts = overlay[1] as usize;
        assert_eq!(n_units, 4, "four carousel units, four circles");
        let max_r = 0.45 * CAROUSEL_R * std::f64::consts::SQRT_2;
        for u in 0..n_units {
            let base = 2 + u * pts * 3;
            // Average over the UNIQUE points only (p in 0..pts-1) — point
            // pts-1 duplicates point 0 (the circle closes at a = TAU), and
            // including it twice biases a naive average off the true
            // center enough to spuriously inflate the recovered radius.
            let mut center = DVec3::ZERO;
            for p in 0..(pts - 1) {
                let o = base + p * 3;
                center += DVec3::new(
                    overlay[o] as f64,
                    overlay[o + 1] as f64,
                    overlay[o + 2] as f64,
                );
            }
            center /= (pts - 1) as f64;
            for p in 0..pts {
                let o = base + p * 3;
                let pt = DVec3::new(
                    overlay[o] as f64,
                    overlay[o + 1] as f64,
                    overlay[o + 2] as f64,
                );
                let r = (pt - center).length();
                assert!(
                    r <= max_r + 1e-6,
                    "carousel overlay circle radius {r} exceeds clamp {max_r}"
                );
            }
        }
    }

    /// A core alpha's post offset (in the nucleus XZ plane, perpendicular
    /// to the axis) rotates ~90° over a quarter roll period; its protons
    /// stay on axis (docs/ATOM_ROTATION_AND_SIM_DESIGN.md §1.8).
    #[test]
    fn alpha_roll_revolves_posts() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        let gid = 0;
        let proton_id = core.profile_id_by_name("proton").unwrap();
        let neutron_id = core.profile_id_by_name("neutron").unwrap();
        let post = *core.groups[gid]
            .members
            .iter()
            .find(|&&id| core.particles[id].profile_id == neutron_id)
            .unwrap();
        let proton = *core.groups[gid]
            .members
            .iter()
            .find(|&&id| core.particles[id].profile_id == proton_id)
            .unwrap();

        let com = core.groups[gid].com;
        let axis = core.groups[gid].orientation * DVec3::Y;
        let perp = |p: DVec3| -> DVec3 {
            let rel = p - com;
            rel - axis * rel.dot(axis)
        };
        let post_perp0 = perp(core.particles[post].position);
        let proton_perp0 = perp(core.particles[proton].position).length();

        let quarter = (TAU / 4.0) / DISPLAY_SPIN_RATE;
        core.advance_display(quarter);

        let post_perp1 = perp(core.particles[post].position);
        let proton_perp1 = perp(core.particles[proton].position).length();

        let cos_ang = post_perp0.normalize().dot(post_perp1.normalize());
        assert!(
            cos_ang.abs() < 0.15,
            "post should revolve ~90° about the alpha axis, cos={cos_ang}"
        );
        assert!(
            proton_perp0 < 1e-6 && proton_perp1 < 1e-6,
            "protons should stay on the alpha axis: {proton_perp0} -> {proton_perp1}"
        );
    }

    /// The exact §0 failure case: a carousel alpha's proton pole (the
    /// alpha's own radial axis) tracks the carousel orbit, AND — once
    /// rolling — its posts revolve about the ALPHA's own axis (not stay
    /// pinned top/bottom of the piece, as the single `car_rot` channel
    /// used to leave them).
    #[test]
    fn carousel_orbit_tracks_pole() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon preset");
        let gid = 0;
        let proton_id = core.profile_id_by_name("proton").unwrap();
        let neutron_id = core.profile_id_by_name("neutron").unwrap();
        let ai = core.groups[gid]
            .alphas
            .iter()
            .position(|a| a.orbits_core)
            .expect("neon should have a carousel alpha");
        let protons: Vec<usize> = core.groups[gid].alphas[ai]
            .members
            .iter()
            .copied()
            .filter(|&k| core.particles[core.groups[gid].members[k]].profile_id == proton_id)
            .map(|k| core.groups[gid].members[k])
            .collect();
        let post = core.groups[gid].alphas[ai]
            .members
            .iter()
            .copied()
            .find(|&k| core.particles[core.groups[gid].members[k]].profile_id == neutron_id)
            .map(|k| core.groups[gid].members[k])
            .unwrap();
        assert_eq!(protons.len(), 2, "carousel alpha should have 2 protons");

        // Phase A: orbit only (this alpha's own roll frozen) — the proton
        // pole must track the carousel orbit.
        core.groups[gid].alphas[ai].roll_rate = 0.0;
        let pole0 = core.particles[protons[0]].pole_axis();
        let quarter_orbit = (TAU / 4.0) / CAROUSEL_VIS_RATE;
        core.advance_display(quarter_orbit);
        let pole1 = core.particles[protons[0]].pole_axis();
        assert!(
            pole0.dot(pole1).abs() < 0.15,
            "carousel proton pole should track the orbit, cos={}",
            pole0.dot(pole1)
        );

        // Phase B: roll only (freeze the carousel now) — the post must
        // revolve about THIS alpha's own axis.
        core.groups[gid].carousel_rate = 0.0;
        core.groups[gid].alphas[ai].roll_rate = DISPLAY_SPIN_RATE;
        let center =
            (core.particles[protons[0]].position + core.particles[protons[1]].position) * 0.5;
        let axis = (core.particles[protons[1]].position - core.particles[protons[0]].position)
            .normalize();
        let perp = |p: DVec3| -> DVec3 {
            let rel = p - center;
            rel - axis * rel.dot(axis)
        };
        let post_perp0 = perp(core.particles[post].position);
        let quarter_roll = (TAU / 4.0) / DISPLAY_SPIN_RATE;
        core.advance_display(quarter_roll);
        let post_perp1 = perp(core.particles[post].position);
        let cos_ang = post_perp0.normalize().dot(post_perp1.normalize());
        assert!(
            cos_ang.abs() < 0.15,
            "post should revolve ~90° about the alpha's OWN axis, cos={cos_ang}"
        );
        assert!(
            (post_perp0.length() - post_perp1.length()).abs() < 0.05,
            "post's distance from its own alpha axis should be preserved (rigid roll)"
        );
    }

    /// A core alpha (orbits_core=false) and a carousel alpha
    /// (orbits_core=true), given equal roll rates, both revolve their
    /// posts about their OWN axes by the same angle — no special-casing
    /// between core and carousel (§1.3 "Key properties").
    #[test]
    fn core_and_carousel_symmetric() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("neon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("neon preset");
        let gid = 0;
        let proton_id = core.profile_id_by_name("proton").unwrap();
        let neutron_id = core.profile_id_by_name("neutron").unwrap();
        let core_ai = core.groups[gid]
            .alphas
            .iter()
            .position(|a| !a.orbits_core)
            .unwrap();
        let car_ai = core.groups[gid]
            .alphas
            .iter()
            .position(|a| a.orbits_core)
            .unwrap();
        // Freeze the carousel orbit so only roll moves anything — an
        // apples-to-apples comparison of the two alphas' own-axis
        // revolution.
        core.groups[gid].carousel_rate = 0.0;
        core.groups[gid].alphas[core_ai].roll_rate = DISPLAY_SPIN_RATE;
        core.groups[gid].alphas[car_ai].roll_rate = DISPLAY_SPIN_RATE;

        let members_of = |g: &RigidGroup, ai: usize, particles: &[SimParticle]| -> (Vec<usize>, usize) {
            let protons: Vec<usize> = g.alphas[ai]
                .members
                .iter()
                .copied()
                .filter(|&k| particles[g.members[k]].profile_id == proton_id)
                .map(|k| g.members[k])
                .collect();
            let post = g.alphas[ai]
                .members
                .iter()
                .copied()
                .find(|&k| particles[g.members[k]].profile_id == neutron_id)
                .map(|k| g.members[k])
                .unwrap();
            (protons, post)
        };
        let (core_protons, core_post) = members_of(&core.groups[gid], core_ai, &core.particles);
        let (car_protons, car_post) = members_of(&core.groups[gid], car_ai, &core.particles);

        let perp_of = |particles: &[SimParticle], protons: &[usize], post: usize| -> DVec3 {
            let center = (particles[protons[0]].position + particles[protons[1]].position) * 0.5;
            let axis =
                (particles[protons[1]].position - particles[protons[0]].position).normalize();
            let rel = particles[post].position - center;
            rel - axis * rel.dot(axis)
        };
        let core_perp0 = perp_of(&core.particles, &core_protons, core_post);
        let car_perp0 = perp_of(&core.particles, &car_protons, car_post);

        let quarter = (TAU / 4.0) / DISPLAY_SPIN_RATE;
        core.advance_display(quarter);

        let core_perp1 = perp_of(&core.particles, &core_protons, core_post);
        let car_perp1 = perp_of(&core.particles, &car_protons, car_post);

        let core_cos = core_perp0.normalize().dot(core_perp1.normalize());
        let car_cos = car_perp0.normalize().dot(car_perp1.normalize());
        assert!(core_cos.abs() < 0.15, "core alpha post should revolve ~90°, cos={core_cos}");
        assert!(car_cos.abs() < 0.15, "carousel alpha post should revolve ~90°, cos={car_cos}");
    }

    /// Two alphas in one nucleus have different roll phases after spawn —
    /// they must not turn in lockstep (the "synchronized ballet",
    /// session-29/session-31 addendum A8).
    #[test]
    fn no_lockstep() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");
        let g = &core.groups[0];
        assert!(g.alphas.len() >= 2, "need multiple alphas");
        let phases: Vec<f64> = g.alphas.iter().map(|a| a.roll_phase).collect();
        let all_equal = phases.windows(2).all(|w| (w[0] - w[1]).abs() < 1e-9);
        assert!(!all_equal, "alphas should not share one roll phase: {phases:?}");
    }

    // ── Part 2: nucleus dynamics modes (session-31 addendum A9–A11) ──────

    /// Mode-aware intra-nucleus skip predicate, probed via `force_accum`
    /// after `compute_forces` (addendum A9):
    /// - RigidLock: no same-group pair interacts at all (current/Part-1
    ///   behavior).
    /// - RigidAlpha: inter-alpha pairs interact (carbon's three alphas feel
    ///   each other), but intra-alpha pairs still don't — isolated with a
    ///   LONE alpha (no neighbors, so any nonzero force would have to come
    ///   from its own members).
    /// - FreeNucleon: never skips — even a lone alpha's own members now
    ///   interact.
    #[test]
    fn dynamics_mode_skip_predicate() {
        let dir = config_dir();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));

        let mut core = AtomCore::new();
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("carbon", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("carbon preset");

        // RigidLock (default): no intra-nucleus force at all.
        core.compute_forces();
        assert!(
            core.particles.iter().all(|p| p.force_accum.length() < 1e-9),
            "RigidLock must skip ALL same-group pairs"
        );

        // RigidAlpha: carbon's core alphas and plugs share the axis —
        // adjacent alpha units must now feel each other.
        core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
        core.compute_forces();
        assert!(
            core.particles.iter().any(|p| p.force_accum.length() > 1e-9),
            "RigidAlpha must let inter-alpha pairs interact"
        );

        // Isolate intra-alpha: a LONE alpha (no neighbors) in RigidAlpha
        // must still see zero force — its own 4 members are one alpha.
        let mut lone = AtomCore::new();
        lone.register_profile("proton", 1.0, 1.0, &p_csv);
        lone.register_profile("neutron", 1.0, 1.0, &n_csv);
        lone.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        lone.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
        lone.compute_forces();
        assert!(
            lone.particles.iter().all(|p| p.force_accum.length() < 1e-9),
            "RigidAlpha must skip intra-alpha pairs (lone alpha has no neighbors)"
        );

        // FreeNucleon: even the lone alpha's own members now interact.
        lone.set_nucleus_dynamics(NucleusDynamics::FreeNucleon);
        lone.compute_forces();
        assert!(
            lone.particles.iter().any(|p| p.force_accum.length() > 1e-9),
            "FreeNucleon must never skip"
        );
    }

    /// `channeling_factor` (session-31 round 4, recapture tail): a
    /// same-group pole-to-hole pair (both sides presenting their pole/hole
    /// to the other, i.e. cos²θ = 1 on at least one side) at the
    /// funnel-mouth distance (r = NUCLEON_PITCH = 2.6) gets
    /// C ≈ channeling·1·1; two facing equators (cos θ = 0 on both sides)
    /// get C ≈ 0 regardless of distance; the OLD cliff radius
    /// (2·NUCLEON_PITCH = 5.2) now sits mid-taper with partial channeling
    /// still active (the recapture basin — this is the whole point of the
    /// tail extension); and r ≥ CHANNEL_TAIL·NUCLEON_PITCH = 7.8 gets C = 0
    /// regardless of angle (channeling is still short-range, just less of a
    /// cliff).
    #[test]
    fn channeling_factor_pole_vs_equator_vs_distance() {
        let channeling = 1.0;

        // Pole-to-hole: both sides present their pole (cos θ = ±1).
        let c_pole = channeling_factor(channeling, 1.0, 1.0, NUCLEON_PITCH);
        assert!(
            (c_pole - channeling).abs() < 1e-9,
            "pole-to-hole pair at r=NUCLEON_PITCH should give C≈channeling: {c_pole}"
        );

        // Facing equators: both sides present their equator (cos θ = 0).
        let c_equator = channeling_factor(channeling, 0.0, 0.0, NUCLEON_PITCH);
        assert!(
            c_equator.abs() < 1e-9,
            "facing-equator pair should give C≈0 regardless of distance: {c_equator}"
        );

        // The OLD cliff radius (2*NUCLEON_PITCH = 5.2) must now sit inside
        // the recapture basin: partial channeling still active, not zero.
        let c_old_cliff = channeling_factor(channeling, 1.0, 1.0, 2.0 * NUCLEON_PITCH);
        assert!(
            c_old_cliff > 0.05 && c_old_cliff < 1.0,
            "old cliff radius should now sit mid-taper (recapture basin), not at C=0: {c_old_cliff}"
        );

        // Far apart (r >= CHANNEL_TAIL*NUCLEON_PITCH): C = 0 even pole-on.
        let c_far = channeling_factor(channeling, 1.0, 1.0, CHANNEL_TAIL * NUCLEON_PITCH);
        assert!(
            c_far.abs() < 1e-9,
            "pole-on pair at r=CHANNEL_TAIL*NUCLEON_PITCH should give C=0 (falloff floor): {c_far}"
        );
        let c_farther = channeling_factor(channeling, 1.0, 1.0, 4.0 * NUCLEON_PITCH);
        assert!(
            c_farther.abs() < 1e-9,
            "pole-on pair beyond CHANNEL_TAIL*NUCLEON_PITCH should give C=0: {c_farther}"
        );
    }

    /// Entering RigidAlpha/FreeNucleon seeds each alpha's rigid-body state
    /// (com/orientation/velocity/angular_velocity) as finite, sane values,
    /// and members are already placed at that seeded frame (addendum A10).
    #[test]
    fn mode_transition_seeds_finite_alpha_bodies() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset(
            "neon",
            DVec3::new(1.0, 2.0, -3.0),
            DVec3::new(0.1, 0.0, 0.0),
            DVec3::Y,
        )
        .expect("neon preset");

        core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
        assert_eq!(core.get_nucleus_dynamics(), NucleusDynamics::RigidAlpha);

        for a in &core.groups[0].alphas {
            assert!(a.com.is_finite(), "alpha com not finite: {:?}", a.com);
            assert!(
                a.velocity.is_finite(),
                "alpha velocity not finite: {:?}",
                a.velocity
            );
            assert!(
                a.angular_velocity.is_finite(),
                "alpha angular_velocity not finite: {:?}",
                a.angular_velocity
            );
            assert!(
                (a.orientation.length() - 1.0).abs() < 1e-6,
                "alpha orientation should be a unit quaternion: {:?}",
                a.orientation
            );
            assert!(a.mass > 0.0, "alpha mass should be positive");
            assert!(a.inertia > 0.0, "alpha inertia should be positive");
        }

        // Members already sit at the seeded alpha frame (place_alpha_members
        // ran as part of the transition) — spot check alpha 0's first member.
        let g = &core.groups[0];
        let a0 = &g.alphas[0];
        let expected_pos = a0.com + a0.orientation * a0.member_local_pos[0];
        let actual_pos = core.particles[g.members[a0.members[0]]].position;
        assert!(
            (expected_pos - actual_pos).length() < 1e-6,
            "member position should match the seeded alpha frame: {expected_pos:?} vs {actual_pos:?}"
        );
    }

    /// Round-trip RigidLock → RigidAlpha → RigidLock: the excursion lets
    /// real rigid-body integration move the alpha (a direct velocity
    /// perturbation on the ALPHA, applied AFTER entering RigidAlpha,
    /// translates it — see below for why this replaced a roll-ω probe);
    /// returning to RigidLock SNAPS back to the kinematic pose (positions
    /// jump, documented, addendum A10) and resumes the phase-driven
    /// contract — immobile under `step()`, moves only via `advance_display`
    /// (same invariant as `display_phases_are_wall_clock`).
    ///
    /// Perturbs `alphas[0].velocity` directly rather than spawning with a
    /// nonzero group velocity, for two reasons: (1) session-31 round 5
    /// found `alpha_kinematic_state` was seeding each alpha's
    /// `angular_velocity`/`velocity` with the cosmetic DISPLAY roll/
    /// carousel rate on top of real group motion — for a lone force-free
    /// alpha (zero group velocity/angular_velocity, this test's original
    /// setup) that fictional roll was the ONLY thing that moved it, which
    /// is exactly the bug this fix removes: a symmetric alpha at its own
    /// force balance with no real perturbation should NOT move under
    /// RigidAlpha, and no longer does. (2) A nonzero SPAWN velocity sets
    /// the CONTAINER-level `RigidGroup::velocity`/`com`, which `kick`/
    /// `drift`'s `RigidLock` arm integrates too (`g.com += g.velocity *
    /// dt` every step, regardless of mode) — RigidAlpha never zeroes or
    /// even touches that container-level field (only per-alpha state), so
    /// it would silently carry through the round-trip and break the
    /// "RigidLock is immobile" check below. Perturbing the alpha's own
    /// `velocity` (populated only by `alpha_kinematic_state`/`kick`'s
    /// RigidAlpha arm, never read by RigidLock) avoids that leak.
    #[test]
    fn rigid_lock_round_trip_resumes_kinematics() {
        let dir = config_dir();
        let mut core = AtomCore::new();
        let p_csv = load_histogram_csv(&dir.join("histogram_proton.csv"));
        let n_csv = load_histogram_csv(&dir.join("histogram_neutron.csv"));
        core.register_profile("proton", 1.0, 1.0, &p_csv);
        core.register_profile("neutron", 1.0, 1.0, &n_csv);
        core.spawn_preset("alpha", DVec3::ZERO, DVec3::ZERO, DVec3::Y)
            .expect("alpha preset");
        core.running = true;

        // The kinematic pose BEFORE the trip. It stays valid as the
        // "resume point" because carousel_phase/roll_phase are frozen
        // while away from RigidLock (advance_display's group branch is
        // gated to RigidLock only), and RigidAlpha never touches the
        // GROUP-level com/orientation/velocity either (container only).
        let kinematic_pose: Vec<DVec3> = core.particles.iter().map(|p| p.position).collect();

        core.set_nucleus_dynamics(NucleusDynamics::RigidAlpha);
        // Direct perturbation of the ALPHA's own velocity (not the
        // container) — a real, non-fictional nudge that should carry the
        // alpha under rigid-body integration without leaking into the
        // container-level state RigidLock reads on re-entry.
        core.groups[0].alphas[0].velocity = DVec3::new(0.05, 0.0, 0.0);
        core.step_n(500);
        core.advance_display(0.3);

        let drifted: Vec<DVec3> = core.particles.iter().map(|p| p.position).collect();
        assert!(
            kinematic_pose
                .iter()
                .zip(&drifted)
                .any(|(a, b)| (*a - *b).length() > 1e-6),
            "RigidAlpha should let the alpha actually translate under a real velocity perturbation"
        );

        core.set_nucleus_dynamics(NucleusDynamics::RigidLock);
        for (expected, p) in kinematic_pose.iter().zip(&core.particles) {
            assert!(
                (*expected - p.position).length() < 1e-6,
                "RigidLock re-entry should snap back to the frozen kinematic pose"
            );
        }

        let before: Vec<DVec3> = core.particles.iter().map(|p| p.position).collect();
        core.step_n(2000);
        for (b, p) in before.iter().zip(&core.particles) {
            assert!(
                (*b - p.position).length() < 1e-9,
                "RigidLock must be immobile under step() (phases only move via advance_display)"
            );
        }
    }
}
