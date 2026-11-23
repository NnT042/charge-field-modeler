//! Config surface for the GPU DSMC photon engine (docs/GPU_DSMC_DESIGN.md §4, §11).
//!
//! Every `[DECISION]` in the design doc is a field here rather than a code
//! branch, so a physics run picks it at construction time. B1 gates only
//! exercise a subset (see the doc's §9 gate table); the rest exist so B2/B3
//! do not need a config rewrite.

use glam::DVec3;

/// D1 — collision-rate kernel (doc §4). `Constant`: rate ∝ σ·n·c, blind to
/// relative motion (Phase-1's rule). `Catch`: rate ∝ σ·n·g with
/// g = |d1 - d2| ∈ [0,2], the standard kinetic rate.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum RateKernel {
    Constant,
    Catch,
}

/// D2 — gear mesh sign (doc §4). `ChiralityOnly`: mesh = χ1·χ2 (Phase-1,
/// direction-blind). `DirRel`: mesh = χ1·χ2·sign(d1·d2), matching the
/// `DirRelSign` rule `apply_photon` already uses for photon-particle hits.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MeshMode {
    ChiralityOnly,
    DirRel,
}

/// Per-body emission mode (doc §3.3). `Forced`: emit at the Phase-1 analogue
/// rate (lines with b <= B_MAX). `Pump`: emit at own absorption rate,
/// smoothed as an EMA over `PUMP_EMA_STEPS` steps.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum EmissionMode {
    Forced,
    Pump,
}

pub const PUMP_EMA_STEPS: f64 = 50.0;

/// D6 — exchange rule (doc §4 D6, extended post-B1). `SwapRenorm`: sample a
/// contact normal, swap its component between d1/d2, renormalise both to c
/// (Phase-1-like; measured ~0.30 photon-momenta drift per collision in G1).
/// `RotateRelative`: rotate the relative vector g = d1-d2 about the mean
/// V = (d1+d2)/2 (or, when |V| ~ 0, about a random axis perp to g) by a
/// glancing-biased scattering angle; conserves |d|=1 and momentum exactly.
/// `RoughSphere` (doc §13) replaces the direction-only exchange with real
/// contact mechanics. Every photon carries a spin vector and a speed, and the
/// pair exchanges linear and angular momentum through the contact law in
/// [`RoughParams`]. It mirrors `crate::rough_sphere::collide_full` line for
/// line. The gear redirect, `spin_transfer` and stiffness shedding do not
/// apply in this mode (`DsmcEngine::new` asserts they are off).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum ExchangeRule {
    SwapRenorm,
    RotateRelative,
    RoughSphere,
}

/// D7 (doc §13.3): what a rough-sphere collision may do to a photon's speed.
/// `Free`: full rough-sphere mechanics, so speeds scatter about c.
/// `PinnedC`: each photon's post-collision velocity is projected back to
/// |v| = c, and the linear energy that removes (or adds) is booked into that
/// photon's own spin. Energy stays exact per photon; momentum does not.
/// `BodyFree` (doc §15.2, user decision 2026-09-29): body contacts are free,
/// so a baryon surface can throw a photon faster (or slower) than c
/// (charge3.html: "The speed of light is the speed at which it is ejected by
/// the baryon"). Photon-photon collisions stay pinned, but each photon is
/// pinned to its OWN incoming speed, not to c. With every photon at c this is
/// PinnedC exactly.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum SpeedMode {
    Free,
    PinnedC,
    BodyFree,
}

/// Doc §15.8: what spin the twist grip compares a photon with at the ball's
/// surface. `EdgeC` (§14.8): the direction of the ball's summed angular
/// velocity, at edge speed `grain_edge`. `Rigid`: the ball's actual rigid
/// angular velocity, in photon spin units (times the photon radius). No
/// assumed edge speed. `RigidInner`: `Rigid` plus the unresolved inner levels
/// 1..(8 - depth) at their Spin-mode rates and class signs, each at a random
/// phase per contact (their periods are far below the keyframe spacing).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GrainMode {
    EdgeC,
    Rigid,
    RigidInner,
}

/// Doc §15.9: the ball's emission law. `Off`: photons bounce off the ball
/// (`GrainMode` contact). Otherwise every ball contact is a capture and the
/// photon is re-emitted from the core through Mathis's maze
/// (`crate::maze::emit`: patch start by helicity, canonical 4-level maze with
/// the body's a, x, y, z signs, trapped photons out of the opposite pole).
/// `Fixed`: the maze frame is the body frame (`orient`). `Swept`: the maze
/// frame turns about the pole with the core's level-9 axial phase (the
/// sprinkler).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MazeMode {
    Off,
    Fixed,
    Swept,
    /// Doc §15.11 control: re-emit in a uniformly random direction from the
    /// core (spin carried rigidly), no maze.
    Isotropic,
}

/// Doc §15.14: the speed a maze body re-emits at. `Relative`: the photon
/// keeps its incoming speed relative to the ball surface and the surface
/// velocity is added (§15.9; repeated recaptures by a ball moving at ~2c
/// random-walk the speed upward, Fermi-style). `EjectC`: the baryon ejects
/// at c relative to its own surface (charge3.html: "The speed of light is
/// the speed at which it is ejected by the baryon"), so the lab speed is
/// |V_surface + c n| whatever the photon arrived with.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum MazeSpeed {
    Relative,
    EjectC,
    /// Doc §15.27 experiment: eject at c in the body's rest frame, without
    /// the ball's centre velocity or spin (the ball's poloidal loop is then
    /// not imprinted on the emission).
    EjectRest,
    /// Doc §15.29: eject at c from the body surface (the envelope) along the
    /// maze exit direction. Same as EjectRest without the inner ball
    /// overtaking (and recapturing) its own emission.
    EjectSurface,
}

/// Doc §15.11: where a maze body takes photons in (scaffolding toward the
/// emergent intake of MATHIS_SPIN_PATHS §D, option 2). `Ball`: only photons
/// that hit the core ball (§15.9). `Envelope { p }`: additionally, a photon
/// crossing into the body envelope (radius 1) is captured with probability
/// p anywhere on it, so any polar focus of the captures must come from the
/// field, not the rule. `Polar { aperture }`: captured only when it enters
/// within `aperture` (fraction of the radius) of the spin axis.
/// `PolarDeflect { aperture }` (doc §15.35): as Polar, but a photon that
/// reaches the envelope outside the opening is batted off the outer spin
/// (a sphere turning about the axis at surface speed c x the level-12 speed
/// fraction at the equator, rough contact law), pushing the body.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Intake {
    Ball,
    Envelope { p: f64 },
    Polar { aperture: f64 },
    PolarDeflect { aperture: f64 },
}

/// Contact law for `ExchangeRule::RoughSphere` (doc §13.2;
/// `crate::rough_sphere::ContactLaw`).
#[derive(Clone, Copy, Debug)]
pub struct RoughParams {
    /// Sliding (edge) grip: -1 = smooth, +1 = perfectly rough.
    pub beta: f64,
    /// Torsional (drilling) grip about the contact normal, same convention.
    pub beta_twist: f64,
    /// Pays the energy the twist grip removes into normal separation
    /// (halbach.pdf).
    pub twist_loss_to_normal: bool,
    /// I = kappa*m*R^2 (0.4 = uniform solid sphere).
    pub kappa: f64,
    pub speed_mode: SpeedMode,
    /// Upper bound on |v1 - v2| for NTC candidate selection. It is 2 when
    /// speeds are pinned at c; Free mode needs headroom, and any pair above
    /// the bound is tallied loudly (TALLY_ACCEPT_OVERFLOW).
    pub g_max: f64,
    /// Birth edge speed |w|*R in units of c. It is 1.0 per photon2.html:
    /// "let the tangential velocity of the spin of the photon equal c".
    pub spin_edge0: f64,
    /// Impact-parameter distribution. false: geometric (b/d = sqrt(rho), rho
    /// uniform), which scatters near-isotropically. true: rho skewed toward 1
    /// (rho = 1-(1-u)^3) exactly as `RotateRelative` does (arp.pdf "glancing"),
    /// so E3 can compare with E2 on equal transport.
    pub glancing_bias: bool,
    /// Doc §13.7 (user rule 2026-09-29): spin stacking. When the outer edge
    /// would exceed c, that level locks at c and the excess starts level k+1
    /// (twice the radius, axis = the collision's spin kick). A stacked level
    /// spun below `unstack_edge` is stripped. Linear mass is 1 at every level.
    /// Off: the E1 model (m = R = 2^k, levels fixed).
    pub stacking: bool,
    /// Highest level index (3 = photon tier z-spin; storage allows <= 3).
    pub stack_max_level: u32,
    pub unstack_edge: f64,
    /// Doc §14.8: body grains spin with their shell (axis = the shell
    /// matter's angular velocity direction) at this edge speed on unit
    /// radius. 0 = non-spinning grains (§14.7 behaviour).
    pub grain_edge: f64,
    /// Doc §15.8: the surface twist a photon meets on the ball.
    pub grain_mode: GrainMode,
    /// Doc §15.9: ball contacts become maze captures (Ball body model only).
    pub maze: MazeMode,
    pub maze_speed: MazeSpeed,
    /// Doc §15.11: intake rule and holding tank. Captured photons are held
    /// at the core for `hold_time` (x `hold_time_neutron` for neutrons;
    /// 0 = re-emitted at once), excluded from collisions, then released
    /// through the maze. While a body holds `hold_cap` photons (0 = no cap)
    /// it takes no more in (they pass or bounce as without the maze).
    pub intake: Intake,
    /// Doc §15.35: neutrons use this intake instead of `intake` when set.
    pub intake_neutron: Option<Intake>,
    pub hold_time: f64,
    pub hold_time_neutron: f64,
    pub hold_cap: u32,
    /// Doc §15.24: through-charge hole. With envelope intake, a photon that
    /// enters within `through_hole` (fraction of the radius) of the spin
    /// axis is not taken in: it crosses the body on its own line, pushing
    /// nothing (venus2.pdf through charge, "water down the drain").
    /// 0 = off.
    pub through_hole: f64,
    /// Doc §15.30: with EjectSurface, a trapped (polar) release leaves along
    /// the surface normal at its exit point on the polar opening (a cone of
    /// half-angle 19.5 deg) instead of parallel to the axis.
    pub polar_release_normal: bool,
    /// Doc §15.32: a trapped release leaves by the pole the photon was
    /// heading toward on entry (instead of the maze's helicity pole).
    pub polar_release_through: bool,
    /// Doc §15.45 trace copies (user, 2026-11-17): each ball body is tested
    /// as this many copies of its own trace ball, copy m being the body at a
    /// fixed time offset along its own trace. An AVERAGING to save compute: it
    /// stands in for a denser field (opacity x copies) while every contact
    /// keeps the trace's real geometry, holes and surface velocity. 1 = the
    /// real ball only (bit-identical to before); max gpu_types::BALL_COPY_MAX.
    pub ball_copies: u32,
    /// Doc §15.46: with EjectSurface, the radius (fraction of the envelope)
    /// of the polar opening a trapped photon's exit point is spread over.
    /// 1/3 = the hole of the trace's ball CENTRES (§15.29); the hole the ball
    /// itself leaves clear is 1/3 - 1/8 = 0.222 (report_trace_clear_hole).
    pub polar_beam: f64,
    /// Doc §15.47 (user, 2026-11-17: dense directed emission is a semisolid
    /// wall that shadows what is behind it): pair photon-photon candidates
    /// with near neighbours (Morton order of 4x4x4 sub-cells) instead of at
    /// random within the cell. The NTC candidate count is unchanged. false =
    /// random in-cell pairing (bit-identical to before).
    pub local_pairing: bool,
    /// Doc §15.51 pipe TEST (user, 2026-11-19: "a stopper, not a pipe"):
    /// bit 1 = protons, bit 2 = neutrons. A photon entering such a body's
    /// envelope within `pipe_mouth` (fraction of the radius) of the axis,
    /// heading inward within acos(`pipe_cos`) of it, is turned exactly onto
    /// the axis. It then runs down the trace's clear hole and leaves by the
    /// far pole unless a photon collision knocks it into the walls. The body
    /// takes only the sideways momentum of the turn. 0 = off (bit-identical).
    pub pipe_classes: u32,
    pub pipe_cos: f64,
    pub pipe_mouth: f64,
    /// Doc §15.52: piped photons take no photon-photon collisions until they
    /// leave the envelope (a protected lane).
    pub pipe_lane: bool,
    /// Doc §15.52 (user's helix idea, 2026-11-19): piped photons follow a
    /// helix about the axis at this angle to it (radians), turning with the
    /// body's outer spin, until their first collision. 0 = straight.
    pub pipe_swirl: f64,
    /// Doc §15.53 channel spin rule (halbach.pdf, test knob): a collision
    /// involving a piped photon does not scatter. Co-moving pairs pass
    /// untouched; a head-on pair whose lab-frame axial spins are parallel
    /// (photon meets antiphoton) passes through and spins up, taking
    /// `pipe_spin_eta` of the relative linear energy into spin; antiparallel
    /// (same species) swaps velocities (bounces back) and spins down, giving
    /// `pipe_spin_eta` of its outer spin energy to linear. Momentum and energy
    /// are conserved exactly. Ignored with `pipe_lane`.
    pub pipe_spin: bool,
    pub pipe_spin_eta: f64,
    /// Doc §15.54 quantum exchange (photon5.pdf, limat.pdf): with
    /// `pipe_spin`, a spin-up gives each photon one saturated level
    /// ((1/2) kappa edge0^2) from the relative linear energy when affordable,
    /// a spin-down strips each photon's outer level to linear. Replaces eta.
    pub pipe_quantum: bool,
    /// Doc §15.55 edge-on passing (user, 2026-11-19; test knob): with
    /// `pipe_spin`, opposed flows in the hole meet side by side, so a
    /// head-on pair never bounces. 0 = off (§15.53 bounce), 1 = every
    /// head-on pair passes; the spin-up/spin-down split is unchanged, a
    /// spin-down's energy goes to linear speed along the photons' own
    /// directions. 2 = as 1 with the split reversed (gear reading: parallel
    /// lab spins grind and despin, antiparallel mesh and spin up).
    pub pipe_edge: u32,
    /// Doc §15.54 radius-held helix (user, 2026-11-19): piped photons with a
    /// swirl follow a helix of this radius (body radii) about a guide line
    /// parallel to the axis: the linear trace of a stacked-spin photon
    /// (`trace_reach_fraction`). 0 = the §15.52 axis-centred swirl.
    pub pipe_helix_r: f64,
    /// Doc §14.8: a separate law for photon-body contacts, as (beta,
    /// beta_twist, twist_loss_to_normal). None = the photon-photon law above.
    pub body_law: Option<(f64, f64, bool)>,
    /// Doc §15.59 photon-size knob (user, 2026-11-21). Test knob; the
    /// default reproduces the engine before it.
    pub photon_mix: PhotonMix,
}

/// How the inner (already full) levels of a born-larger ambient photon point.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum InnerAxes {
    /// Along travel with the outer level's sense (chi * dir).
    Match,
    /// All +dir.
    Plus,
    /// All -dir.
    Minus,
    /// Each inner level on its own isotropic random axis.
    Random,
}

/// Doc §15.59 (user 2026-11-20/21): ambient photons born at stacking level
/// k with every inner level full (edge c), outer level at edge c along
/// travel with the chirality chi (grav4.pdf: photon/antiphoton = the outer
/// spin's left/right). Mass per birth level and a hit box set by the outer
/// spin radius (arp.pdf: larger photons suffer more collisions).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PhotonMix {
    /// Share of ambient photons born at k = 0..3 (normalised on upload).
    pub weights: [f64; 4],
    pub inner: InnerAxes,
    /// Stacking mode only: linear mass = 2^(birth level) instead of 1
    /// (each spin doubles the energy, photon2/howell3). Collision stacking
    /// later keeps the birth mass.
    pub mass_doubling: bool,
    /// Stacking mode only: pair cross-section sigma * ((R1+R2)/2)^2 with
    /// R = 2^(current outer level).
    pub size_xsec: bool,
    /// NTC ceiling level for `size_xsec` (candidates drawn at
    /// sigma * 4^ceiling). Pairs above it are clamped and tallied.
    pub size_ceiling: u32,
    /// Test knob (doc §15.59, user 2026-11-23): how the mixed field is
    /// normalised against the old one-size field.
    pub norm: MixNorm,
}

/// Doc §15.59 field normalisation for the photon mix (test knob, off by
/// default).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MixNorm {
    Off,
    /// Same field mass (user): heavier photons, so fewer of them. Every
    /// photon mass divided by the mix mean <m>, and the base cross-section
    /// divided by <m> (fewer real photons per simulated one).
    Mass,
    /// Option B: field collision rate held at the old one-size value. Base
    /// cross-section divided by the pair-averaged hit-box factor
    /// <((m_a+m_b)/2)^2> (1 without the hit box); masses unchanged.
    Path,
    /// Both: masses divided by <m> and collision rate held (as Path).
    MassPath,
}

impl Default for PhotonMix {
    fn default() -> Self {
        PhotonMix { weights: [1.0, 0.0, 0.0, 0.0], inner: InnerAxes::Match, mass_doubling: false, size_xsec: false, size_ceiling: 0, norm: MixNorm::Off }
    }
}

impl PhotonMix {
    fn shares(&self) -> [f64; 4] {
        let s: f64 = self.weights.iter().sum();
        [self.weights[0] / s, self.weights[1] / s, self.weights[2] / s, self.weights[3] / s]
    }

    /// Mix mean birth mass <m> = sum w_k 2^k.
    pub fn mean_mass(&self) -> f64 {
        self.shares().iter().enumerate().map(|(k, w)| w * 2f64.powi(k as i32)).sum()
    }

    /// Pair-averaged hit-box factor <((2^a + 2^b)/2)^2> over birth levels
    /// (1 when the hit box is off). Birth levels only: stacking in flight
    /// raises it a little, so the held rate is approximate.
    pub fn size_factor(&self) -> f64 {
        if !self.size_xsec {
            return 1.0;
        }
        let w = self.shares();
        let mut f = 0.0;
        for a in 0..4 {
            for b in 0..4 {
                let r = (2f64.powi(a) + 2f64.powi(b)) * 0.5;
                f += w[a as usize] * w[b as usize] * r * r;
            }
        }
        f
    }

    /// Multiplier on every photon mass (applied only with mass doubling).
    pub fn mass_scale(&self) -> f64 {
        match self.norm {
            MixNorm::Mass | MixNorm::MassPath if self.mass_doubling => 1.0 / self.mean_mass(),
            _ => 1.0,
        }
    }

    /// Multiplier on the base cross-section sigma.
    pub fn sigma_scale(&self) -> f64 {
        match self.norm {
            MixNorm::Off => 1.0,
            MixNorm::Mass => if self.mass_doubling { 1.0 / self.mean_mass() } else { 1.0 },
            MixNorm::Path | MixNorm::MassPath => 1.0 / self.size_factor(),
        }
    }

    /// GPU packing: cumulative shares c0..c2 and flags (bit 0 mass, bits
    /// 1-2 inner, bit 3 size, bits 4-5 ceiling).
    pub fn pack(&self) -> ([f32; 3], u32) {
        let s: f64 = self.weights.iter().sum();
        assert!(s > 0.0 && self.weights.iter().all(|w| *w >= 0.0), "photon mix weights");
        let c0 = self.weights[0] / s;
        let c1 = c0 + self.weights[1] / s;
        let c2 = c1 + self.weights[2] / s;
        let inner = match self.inner {
            InnerAxes::Match => 0,
            InnerAxes::Plus => 1,
            InnerAxes::Minus => 2,
            InnerAxes::Random => 3,
        };
        ([c0 as f32, c1 as f32, c2 as f32], self.mass_doubling as u32 | (inner << 1) | ((self.size_xsec as u32) << 3) | (self.size_ceiling.min(3) << 4))
    }
}

impl Default for RoughParams {
    fn default() -> Self {
        RoughParams { beta: 1.0, beta_twist: -1.0, twist_loss_to_normal: false, kappa: 0.4, speed_mode: SpeedMode::PinnedC, g_max: 2.0, spin_edge0: 1.0, glancing_bias: false, stacking: false, stack_max_level: 3, unstack_edge: 0.01, grain_edge: 0.0, grain_mode: GrainMode::EdgeC, maze: MazeMode::Off, maze_speed: MazeSpeed::Relative, intake: Intake::Ball, intake_neutron: None, hold_time: 0.0, hold_time_neutron: 1.0, hold_cap: 0, through_hole: 0.0, polar_release_normal: false, polar_release_through: false, ball_copies: 1, polar_beam: 1.0 / 3.0, local_pairing: false, pipe_classes: 0, pipe_cos: 0.939_692_620_785_908_4, pipe_mouth: 1.0 / 3.0 - 1.0 / 8.0, pipe_lane: false, pipe_swirl: 0.0, pipe_spin: false, pipe_spin_eta: 0.0, pipe_quantum: false, pipe_edge: 0, pipe_helix_r: 0.0, body_law: None, photon_mix: PhotonMix { weights: [1.0, 0.0, 0.0, 0.0], inner: InnerAxes::Match, mass_doubling: false, size_xsec: false, size_ceiling: 0, norm: MixNorm::Off } }
    }
}

/// E1 -- per-photon stiffness shedding rule (doc §11d). `None`: a photon's
/// stiffness level k is permanent once emitted. `Transfer` (default for
/// physics runs): each accepted collision moves one level from the higher-k
/// photon to the lower-k one (k unchanged if equal) -- conserves the total
/// stacked-level count across the box, applied AFTER the exchange (which
/// uses pre-collision k/mass).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum StiffnessShed {
    None,
    Transfer,
}

/// Doc §14: what a body is. `Absorber` is the legacy absorb/re-emit sphere
/// (pump or forced emission from a fixed exit profile). `Shells` is the
/// mechanical baryon: levels 9-12 as nested permeable shells moving with the
/// Spin-mode kinematics (`crate::baryon_shells`). Photons pass through and
/// meet body grains at rate `tau_body * |v - V(x)|`. There is no absorption
/// and no emission; all charge is scattered ambient.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum BodyModel {
    Absorber,
    Shells { tau_body: f64 },
    /// Doc §15: the instantaneous base ball (levels 1-8, radius 1/8) carried
    /// by levels 9-12. It is the only body matter. Hits are exact
    /// moving-sphere contacts; the envelope does bookkeeping only.
    Ball,
}

/// Nucleon seating in a stack (recycling.rs `StackOrient`, duplicated here so
/// the GPU module has no dependency on the CPU Phase-1 module).
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GpuOrient {
    PoleOn,
    EdgeOn,
    /// Doc §15.18: pole-on but upside down (spin axis -z), e.g. the lower
    /// deuteron of Mathis's He sandwich (deut.pdf).
    PoleDown,
}

impl GpuOrient {
    pub fn spin_axis(self) -> DVec3 {
        match self {
            GpuOrient::PoleOn => DVec3::Z,
            GpuOrient::EdgeOn => DVec3::X,
            GpuOrient::PoleDown => -DVec3::Z,
        }
    }
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
pub enum GpuClass {
    Proton,
    Neutron,
    Electron,
}

impl GpuClass {
    pub fn chirality_sign(self) -> f64 {
        match self {
            GpuClass::Proton => 1.0,
            GpuClass::Neutron => -1.0,
            GpuClass::Electron => 1.0,
        }
    }

    /// (peak0_deg, peak1_deg_or_none, sigma_deg) — same numbers as
    /// recycling.rs's EXIT_PEAK_LAT_DEG / ELECTRON_EXIT_PEAK_LAT_DEG.
    pub fn exit_profile(self) -> (f64, Option<f64>, f64) {
        match self {
            GpuClass::Proton | GpuClass::Neutron => (7.0, Some(58.0), 5.0),
            GpuClass::Electron => (3.0, None, 2.5),
        }
    }

    pub fn id(self) -> u32 {
        match self {
            GpuClass::Proton => 0,
            GpuClass::Neutron => 1,
            GpuClass::Electron => 2,
        }
    }
}

/// One static body in the stack (doc §5).
#[derive(Clone, Copy, Debug)]
pub struct BodyConfig {
    pub center: DVec3,
    pub class: GpuClass,
    pub orient: GpuOrient,
    /// Through-lane radius (fraction of R_IN). 0 = lane shut for this body.
    pub through_lane_rho: f64,
}

/// The domain boundary. `SphereEscape` is the doc's normal geometry
/// (escape at R_D, centred on the body centroid). `PeriodicBox` is a G1-only
/// convenience: a cubic box with wraparound, so an ambient-only isotropic
/// population never needs injection/escape bookkeeping to hold `n_a` fixed.
#[derive(Clone, Copy, Debug)]
pub enum Boundary {
    SphereEscape { r_d: f64 },
    PeriodicBox { half_extent: f64 },
}

#[derive(Clone, Debug)]
pub struct DsmcConfig {
    pub boundary: Boundary,
    pub bodies: Vec<BodyConfig>,

    /// Cell size (world units). Doc default 0.25 for physics runs; gates may
    /// use a coarser grid for speed.
    pub h: f64,
    /// c·dt per step (c = 1 in these units).
    pub dt: f64,
    /// Max resident photon count (buffer sizing).
    pub capacity: u32,
    pub seed: u64,

    // --- physical knobs ---
    /// Ambient number density (photons per unit volume) used for injection
    /// flux and the analytic G1/G0 comparisons.
    pub n_ambient: f64,
    /// Cross-section σ (1/λ · n_a already folded out — this is the bare σ
    /// used as rate = σ·n·g). σ = 0 disables all collisions (G0).
    pub sigma: f64,
    pub photon_fraction: f64,
    pub pole_intake_bias: f64,
    pub cancel_redirect: f64,
    pub spin_transfer: f64,
    pub swing_boost: f64,

    pub rate_kernel: RateKernel,
    pub mesh_mode: MeshMode,
    pub exchange_rule: ExchangeRule,
    pub emission_mode: EmissionMode,
    /// D4: per-photon bounce cap; 0 disables the cap.
    pub bounce_cap: u32,
    /// D5: far-zone (r > B_MAX) ambient-ambient collisions off when true.
    pub d5_far_zone_gate: bool,
    pub b_max: f64,

    // --- E1: emitted-photon stiffness (doc §11d) ---
    /// k0: the stiffness level (0..7) emitted photons are born with; m = 2^k0
    /// (bottom.pdf). Ambient photons are always born at k=0. 0 disables the
    /// feature entirely -- every photon stays at m=1, so the generalised
    /// exchange always takes the equal-mass (bit-identical-to-HEAD) branch.
    /// Also used as the NTC size-scaling upper bound (see
    /// `sigma_size_scaling`): since a Transfer collision only ever moves a
    /// level DOWN from its current holder, no photon's k can exceed the
    /// highest k present at emission, so k0 is a safe, config-known ceiling.
    pub emit_stiffness_level: u32,
    pub stiffness_shed: StiffnessShed,
    /// arp.pdf size-driven cross-section: when true, sigma_pair =
    /// sigma * ((m1+m2)/2)^2 (photon "radius" proportional to m), with the
    /// NTC candidate count using sigma_max = sigma * (2^emit_stiffness_level)^2
    /// so acceptance stays a proper accept/reject ratio. When false,
    /// sigma_pair = sigma unconditionally (bit-identical to HEAD).
    pub sigma_size_scaling: bool,
    /// Test-only negative control for `e1_unequal_exchange_conserves_momentum`
    /// (doc §11d gate 1): forces the exchange to take the EQUAL-mass branch
    /// (rotate about V=(d1+d2)/2, ignoring m1/m2) even when k1 != k2 -- the
    /// "deliberately wrong variant" the gate needs to show its momentum check
    /// is sensitive. Always false outside that one gate.
    pub debug_force_equal_mass_exchange: bool,
    /// Doc §13: contact law and D7 speed mode, used only when
    /// `exchange_rule == RoughSphere`.
    pub rough: RoughParams,
    /// Doc §14 body model. Shells requires `ExchangeRule::RoughSphere`.
    pub body_model: BodyModel,
    /// Shells: per-body level signs (a, x, y, z). An empty vector uses each
    /// class's first listed state (`baryon_shells::{PROTON,NEUTRON}_STATES[0]`).
    pub shell_signs: Vec<[f64; 4]>,
}

impl DsmcConfig {
    pub fn cos_max_angle() -> f64 {
        // cos(20 deg), same constant as recycling.rs::THROUGH_COS_MAX_ANGLE.
        0.939_692_620_785_908_4
    }
}
