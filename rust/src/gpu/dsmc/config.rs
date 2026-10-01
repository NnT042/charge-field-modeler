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
}

/// Doc §15.11: where a maze body takes photons in (scaffolding toward the
/// emergent intake of MATHIS_SPIN_PATHS §D, option 2). `Ball`: only photons
/// that hit the core ball (§15.9). `Envelope { p }`: additionally, a photon
/// crossing into the body envelope (radius 1) is captured with probability
/// p anywhere on it, so any polar focus of the captures must come from the
/// field, not the rule. `Polar { aperture }`: captured only when it enters
/// within `aperture` (fraction of the radius) of the spin axis.
#[derive(Clone, Copy, PartialEq, Debug)]
pub enum Intake {
    Ball,
    Envelope { p: f64 },
    Polar { aperture: f64 },
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
    pub hold_time: f64,
    pub hold_time_neutron: f64,
    pub hold_cap: u32,
    /// Doc §14.8: a separate law for photon-body contacts, as (beta,
    /// beta_twist, twist_loss_to_normal). None = the photon-photon law above.
    pub body_law: Option<(f64, f64, bool)>,
}

impl Default for RoughParams {
    fn default() -> Self {
        RoughParams { beta: 1.0, beta_twist: -1.0, twist_loss_to_normal: false, kappa: 0.4, speed_mode: SpeedMode::PinnedC, g_max: 2.0, spin_edge0: 1.0, glancing_bias: false, stacking: false, stack_max_level: 3, unstack_edge: 0.01, grain_edge: 0.0, grain_mode: GrainMode::EdgeC, maze: MazeMode::Off, maze_speed: MazeSpeed::Relative, intake: Intake::Ball, hold_time: 0.0, hold_time_neutron: 1.0, hold_cap: 0, body_law: None }
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
}

impl GpuOrient {
    pub fn spin_axis(self) -> DVec3 {
        match self {
            GpuOrient::PoleOn => DVec3::Z,
            GpuOrient::EdgeOn => DVec3::X,
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
