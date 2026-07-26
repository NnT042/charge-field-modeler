//! CM-3: the disc-to-pole recycling loop (diagnostic) — v3.
//!
//! The claim under test (per the papers, driving everything from local atomic
//! fields to the solar system): the roughly disc-shaped emission profile
//! meets the incoming ambient field, turns it away, and channels it along the
//! polar lines to be pulled back in and recycled — and this same loop may be
//! what sustains the particle's own outer rotation (the swing/`outer_spin`
//! level).
//!
//! Sources for the mechanism (see project MEMORY.md / docs, cite in place of
//! re-deriving from milesmathis.com per project convention):
//!   - encel2.pdf: photons only collide where matter channels them — at the
//!     nuclear boundary, into the face of incoming ambient charge. Emission
//!     density ~1/r² peaks near the emitter, so that's where the collision
//!     zone concentrates.
//!   - cu.pdf: near-nucleus steering is photon-to-photon and "takes multiple
//!     hits" — hence `n_ricochet` (bounded bounce count) rather than a single
//!     deflection.
//!   - arp.pdf: photon-photon hits are glancing, almost never head-on — the
//!     collision model below biases its sampled contact normal toward
//!     grazing incidence for this reason (see `sample_contact_normal`).
//!   - aurora.pdf: polar intake = pressure lows at the poles of any spinning
//!     sphere; a crossover boundary exists where incoming charge = outgoing.
//!   - drift.pdf / venus2.pdf: loop shape is in-at-poles, out-heaviest around
//!     ~30° latitude. venus2's photon-south/antiphoton-north asymmetry is a
//!     later refinement NOT modeled here (this is a single-species pass).
//!   - Emission profile: this project's OWN measured proton exit structure —
//!     bimodal dual peaks at 7° and 58° latitude (M2-era trace measurement,
//!     see project memory `reference_bimodal_exit_structure`).
//!
//! Units: natural (particle surface radius = 1, c = 1), matching project
//! convention (all internal math in natural units; SI lives only in
//! `units.rs`, not touched by this module).
//!
//! ## v2 knob redesign (why v1's null was uninformative)
//!
//! The v1 run (commit e1da3f2) produced an honest but knob-broken null,
//! diagnosed three ways:
//!   1. Unbiased ambient launch from r=30 aimed only ~0.03% of photons
//!      anywhere near the particle — every conditional statistic sat on
//!      ~150 events.
//!   2. sigma_eff was scaled from the *tallied* iteration-0 equator line
//!      integral (~21k), which was dominated by the launch-shell density
//!      spike at r=1.001 (every emitted photon crosses a hair-thin shell of
//!      tiny cell volume) — so virtually ALL of TAU concentrated in that
//!      skin and the rest of the domain stayed transparent, which is why
//!      arrivals barely moved between TAU=0 and TAU=10.
//!   3. Emitted photons saw effectively zero ambient optical depth
//!      (99.94–100% escaped uncollided): the recycling half of the loop was
//!      never exercised.
//!
//! v2 fixes, in order:
//!   1. **Importance-sampled ambient launch**: every ambient photon travels a
//!      straight line whose impact parameter b (closest approach to origin)
//!      is area-weighted over [0, B_MAX=6]: b = B_MAX·sqrt(u). For a uniform
//!      isotropic external field the flux through an impact-parameter annulus
//!      is ∝ b·db, so area weighting IS the physically-correct relative
//!      weighting *within* the b ≤ B_MAX ensemble — all reported ambient
//!      stats are CONDITIONAL on this near-zone ensemble (lines passing
//!      within 6 radii); absolute capture fractions can be recovered
//!      analytically later. b ≤ 1 photons aim at the surface (aim latitude
//!      recorded as before); b > 1 photons are the aim=MISS population, now
//!      guaranteed interesting (all pass within 6 radii).
//!   2. **Analytic regularized collision fields**: collision densities are
//!      analytic seeds — emitted n_e(r,θ) = profile(θ)/r² with profile = the
//!      bimodal exit-latitude mixture pdf; ambient n_a = 1 (its natural
//!      far-field value in these units) — times a per-cell multiplicative
//!      correction from the tallies, capped to [0.25, 4]×. This kills the
//!      launch-shell spike (the cap bounds it at 4× analytic) while keeping
//!      the self-consistent feedback iteration meaningful. The ambient
//!      correction is only applied for r ≤ B_MAX: a point at radius r is
//!      crossed only by lines with b ≤ r, so the b ≤ B_MAX ensemble gives
//!      COMPLETE (unbiased-shape) line coverage exactly for r ≤ B_MAX; for
//!      r > B_MAX lines with b > B_MAX are missing from the tally, so the
//!      correction there is pinned to 1 (= the true far-field value).
//!   3. **Symmetric TAU loading**: TAU_in = expected collisions for a b=0
//!      equatorial ambient photon through analytic n_e over r = 6 → 1
//!      (excluding r < 1.05, so no skin-bin dominance), TAU_out = expected
//!      collisions for an equatorial emitted photon r = 1 → R_OUT through
//!      n_a. One swept knob TAU sets both: sigma_e (ambient-through-emitted)
//!      from TAU_in = TAU, sigma_a (emitted-through-ambient) from
//!      TAU_out = TAU. Both sigmas printed.
//!   4. **Knob verification**: the report tallies and prints the MEASURED
//!      mean collisions per ambient near-zone photon and per emitted photon
//!      for every sweep cell, and prints a loud WARNING when they are far
//!      from O(TAU) instead of proceeding silently.
//!
//! ## v3: same-population collisions — the photon gas gets pressure
//!
//! v2's verdict (commit bf5beb0): every PIECE of the loop is real and
//! monotonic in TAU (screening, deflection-capture intake, 1.2% emitted
//! recycling, the geometric near-surface polar inflow window) — but the cell
//! does NOT self-organize: captured flux funnels onto the emission disc,
//! equator-aim→polar redirection stays <1%, and large-r polar columns drift
//! OUTWARD with collisions on. Diagnosis: only CROSS-population binary
//! collisions exist, so the gas is ballistic and aurora.pdf's actual intake
//! mechanism ("any spinning sphere creates field potentials with lows at the
//! poles") has no mechanical carrier. Pressure requires a gas that pushes on
//! ITSELF. v3 adds the self-collision channel:
//!
//!   1. **Two channels per photon.** Each photon marches through BOTH the
//!      other population's field (cross — sigma_e/sigma_a from TAU, exactly
//!      as v2) and its OWN population's field (self — sigma_ee/sigma_aa from
//!      the new TAU_SELF knob). Per step the two rates add; on a collision
//!      the partner population is chosen proportional to its local rate.
//!      With TAU_SELF = 0 the RNG draw sequence is bit-identical to v2, so
//!      the (TAU=3, TAU_SELF=0) sweep cell is a live v2 regression anchor.
//!   2. **TAU_SELF calibration mirrors v2 fix #3**: sigma_ee such that an
//!      equatorial emitted photon accumulates TAU_SELF expected self-hits
//!      through its own analytic profile over r = 1.05 → R_OUT; sigma_aa
//!      such that a b=0 ambient photon accumulates TAU_SELF expected
//!      self-hits crossing the near zone (r = B_MAX → 1, path B_MAX − 1)
//!      through n_a = 1. Both sigmas printed.
//!   3. **Ambient-ambient far-zone gate (r ≤ B_MAX only).** A uniform
//!      isotropic gas scattering off itself is ensemble-invariant (isotropic
//!      in → isotropic out, detailed balance), so far-zone self-collisions
//!      carry no signal — but they WOULD destroy the importance-sampled aim
//!      bookkeeping and lengthen every path. Self-collisions therefore act
//!      only where the near zone's gradients make them meaningful. The
//!      emitted population needs no gate (its tally is complete everywhere
//!      and its density dies as 1/r² anyway).
//!   4. **Knob verification per channel** (v1→v2 lesson, extended): the
//!      report prints MEASURED mean cross-collisions AND self-collisions per
//!      photon for every cell, and shouts when either active channel is far
//!      off its target. Note the shared ricochet cap truncates both channels
//!      jointly when TAU + TAU_SELF exceeds it.
//!   5. **New observables for circulation**: a meridional-flux accumulator
//!      (dir·θ̂_local — the loop signature is poleward drift, antisymmetric
//!      about the equator) and a printed PRESSURE map (the local collision
//!      rate per unit length an ambient photon feels) so the disc-ridge-high
//!      / polar-low landscape and any flow down its gradient are visible
//!      directly.
//!
//! Mechanism note (why this is not just more noise): partner sampling uses
//! the cell's mean direction and coherence, and the glancing exchange
//! transfers nothing when photon and partner co-move — so coherent streams
//! do not jostle themselves (pole.pdf gearing: "co-moving photons jostle,
//! opposing gears catch"). Pressure emerges only where flows cross or pile
//! up, which is exactly the disc ridge meeting the ambient inflow.
//!
//! ## v4: the pump surface + the charge-pause bubble
//!
//! v3's verdict (commit 1b66357): pressure is a REAL carrier — the polar-low
//! landscape is measured and the gas responds down-gradient monotonically in
//! TAU_SELF — but the cell still doesn't close: the poles stay exhaust.
//! Diagnosis: the particle was an INFINITE ONE-WAY SOURCE (emission spawned
//! unconditionally, absorption passive), so nothing coupled intake to output
//! and pressure could only redirect a wind, never wrap it into a loop.
//!
//! v4 makes the particle an ENGINE (user decision, option (a)):
//!
//!   1. **Endogenous emission budget (the pump).** In pump mode the emitted
//!      population's per-iteration launch count is no longer fixed at N — it
//!      equals the PREVIOUS iteration's measured absorption (ambient
//!      arrivals + emitted re-entries). Intake determines output; photon
//!      number is conserved by the surface in steady state; the settled
//!      budget is the engine's throughput and is printed per iteration.
//!      Sources: eccen.html ("you have two winds... The Sun takes in this
//!      wind at the poles, due to spin, and emits it at the equator");
//!      pause.html ("charge goes in at the poles... all matter is an
//!      engine"). Collision-field STRENGTH stays knob-set (TAU/TAU_SELF
//!      against the analytic seeds, as v2/v3) — the pump changes the source
//!      topology, not the coupling calibration; tying coupling strength to
//!      measured throughput is a further rung, deliberately not taken yet.
//!   2. **The charge-pause bubble (user's heliopause picture, 2026-07-18).**
//!      The emitted wind meeting the ambient inflow should form a balance
//!      boundary — pause.html's magnetopause IS "the boundary between the
//!      charge field of the Earth and the charge field of the Sun"
//!      ("the magnetopause is at the charge-pause"); startrek.pdf: two
//!      photon fields meeting head to head "causes a boundary. There must be
//!      some altitude at which the local energy (density) of one field
//!      equals that of the other. We would expect to see special effects at
//!      this boundary."; eccen.html: returning charge "will be compressed
//!      and its density will be increased". New observables, all measured
//!      not assumed: (i) the CROSSOVER SURFACE r_c(θ) — outermost radius
//!      with net inward combined radial flow, per θ (prediction: bulged at
//!      the disc where emission pushes ambient back, pinched at the poles =
//!      the intake funnel); (ii) SHELL DENSITY PROFILES per latitude band
//!      (looking for the compression ridge — note the correction cap bounds
//!      any tally ridge at 4x analytic, so a stronger shell would be
//!      under-reported; documented, revisit if the ridge saturates the
//!      cap); (iii) the ambient ABSORPTION latitude split (does intake
//!      shift poleward once the pump runs — the in-at-the-poles signature).
//!
//! ## Deliberate deviations / disambiguations from the literal spec (v1,
//! still in force)
//!
//! 1. **Same-hemisphere emission direction.** The spec draws the emission
//!    direction's latitude from "a NEW draw from the same [bimodal] mixture"
//!    independently of the start latitude. Taken fully independently, the
//!    mixture's two peaks (7°/58°) are far enough apart (58° − (−58°) = 116°
//!    > 90°) that ~12.5% of emitted photons would have direction and position
//!    latitudes on OPPOSITE hemispheres — making the "radially outward" ray
//!    actually point back into the sphere from the get-go. Fix: the
//!    hemisphere sign is drawn ONCE per photon and shared by both draws;
//!    only the peak-magnitude (7° vs 58°) is redrawn independently. Max
//!    separation is then < 90°, so the exit ray is always genuinely outward.
//! 2. **Circulation-matrix accumulator.** Dotting the cell's azimuthally-
//!    folded Cartesian direction sum against a single reference r̂ cancels
//!    identically at the equator (r̂ is pure x/y there and rotates with φ).
//!    Fix: a separate `radial_flux` accumulator deposits
//!    `weight·ds·(dir·r̂_local)` with r̂_local evaluated at each photon's own
//!    position, BEFORE the azimuthal fold. `dir_accum` itself is kept as
//!    specified for the collision partner-sampling / coherence use.

#![allow(dead_code)]

use glam::{DQuat, DVec3};
use std::f64::consts::{FRAC_PI_2, PI, TAU};
use std::sync::OnceLock;

use crate::calibration::{rand_unit_vec, xorshift64};

// ---------------------------------------------------------------------
// Grid geometry
// ---------------------------------------------------------------------

/// Particle surface radius (natural units).
pub const R_IN: f64 = 1.0;
/// Outer boundary of the tallied shell.
pub const R_OUT: f64 = 30.0;
/// Log-spaced radial bins over [R_IN, R_OUT].
pub const N_R: usize = 24;
/// Uniform polar-angle (colatitude) bins over [0, π].
pub const N_THETA: usize = 18;
const GRID_CELLS: usize = N_R * N_THETA;

/// Self-consistent field-iteration passes (feedback-loop settling).
pub const ITERS: usize = 4;
/// Pump-mode passes (v4): the endogenous emission budget needs extra
/// iterations to settle (n_emit* = A/(1 − reabsorption_fraction) is reached
/// geometrically; ~3 iters after the ambient field stabilizes).
pub const ITERS_PUMP: usize = 8;
/// Photons per population per pass (report default; unit tests use less).
pub const DEFAULT_N: usize = 200_000;
/// Tangential exit-direction boost — the settled proton-scale outer_spin
/// rate (see project memory: presets settled around ±0.065, "outer velocity:
/// non-relativistic only" landed near 0.055c). Used only as this module's
/// own default; not read from `calibration::CalibrationParticle`.
pub const DEFAULT_SWING_BOOST: f64 = 0.05;
/// Ambient photon fraction near Earth (magmom.pdf: local field ~2/3 photon,
/// 1/3 antiphoton, the Sun's left-spin imbalance). This asymmetry is what
/// breaks the chirality-gear symmetry: at 0.5 (balanced) proton and neutron
/// mesh identically with the field (neutron.pdf "balanced field → near-total
/// spin cancellation, neutron indistinguishable"); only the 2/3 imbalance
/// makes the two classes diverge. Swept vs 0.5 as the control.
pub const PHOTON_FRACTION: f64 = 2.0 / 3.0;
/// Through-channel strength for the chirality gears. On an opposite-chirality
/// ("single-negative", Compton) collision the charge is blocked from
/// equatorial escape and "gets funneled back to the central spin axis... and
/// from there it can only escape back out one of the poles" (neutron.pdf;
/// PHYSICS_REFERENCE §5, halbach.pdf). We model that funnel by blending the
/// exit direction toward the nearer pole axis: dir = lerp(dir, ±ẑ, k). 0 = no
/// funnel (soft, like augment → disc escape); 1 = full axial redirect (out
/// the pole). WHICH charge cancels is set by chirality × field imbalance, so
/// the polar-exit RATE is emergent; only the pole DIRECTION is the sourced
/// mechanism. Swept in the report as the through-channel knob.
pub const DEFAULT_CANCEL_REDIRECT: f64 = 1.0;
/// Per-collision spin transfer for the magnetic channel — how much of a
/// partner's spin a photon picks up when they mesh (magmom.pdf: "those spins
/// can stack or cancel"). This is the ONE knob on the magnetic observable, so
/// the report SWEEPS it and the qualitative predictions (balanced field →
/// neutron magnetically flat; Earth 2/3 → suppressed but nonzero) are required
/// to hold across the sweep rather than at a tuned value. 0 = frozen spin =
/// the pre-existing behaviour, in which cancellation was unmeasurable.
pub const DEFAULT_SPIN_TRANSFER: f64 = 0.25;
/// Through-charge lane radius as a fraction of the body radius — the polar
/// "hole" that charge can cross without being pulled into the equatorial
/// whirlpool (salt.pdf; venus2.pdf "only charge that travels very near the pole
/// can do this, and it has to enter the body on the right trajectory, too").
///
/// NOT a tuned knob: this is the MEASURED pass-through lane from
/// `calibration::report_axial_channel_geometry` (commit 7516470) — proton hole
/// radius 127.5 against a disc reaching 382.5, i.e. rho < 1/3 of the body
/// radius, about 11% of the disc area. On a unit sphere that lane subtends a
/// polar cap of asin(1/3) ≈ 19.5°, which lands on the same 10-20° axial cone
/// already wired into `apply_photon`'s intake branch — two independent routes to
/// the same aperture. Note the SAME map was measured for proton and neutron
/// (structural identity: the L12 loop is the L11 pattern riding the level-12
/// orbit), so the lane is deliberately CLASS-INDEPENDENT here and any
/// proton/neutron asymmetry in through-charge has to emerge from the gearing and
/// the exit profiles rather than from a per-class aperture.
pub const THROUGH_LANE_RHO: f64 = 1.0 / 3.0;
/// Axial-incidence gate for through-charge: |dir·ẑ| must exceed this, i.e. the
/// photon arrives within 20° of the pole axis. Charge arriving "on an angle" is
/// spun into the equatorial whirlpool instead (salt.pdf). 20° matches
/// `apply_photon`'s COS_EDGE so the two axial branches agree on what "axial"
/// means.
pub const THROUGH_COS_MAX_ANGLE: f64 = 0.939_692_620_785_908_4; // cos(20°)
/// Maximum ambient-launch impact parameter (importance-sampling near zone).
pub const B_MAX: f64 = 6.0;
/// Safety cap on total path length before a photon is given up as lost.
const SAFETY_PATH_CAP: f64 = 20.0 * R_OUT;
const MIN_DS: f64 = 0.02;
const DS_FRACTION: f64 = 0.05;
/// Tally-vs-analytic multiplicative correction bounds (v2 fix #2).
const CORRECTION_MIN: f64 = 0.25;
const CORRECTION_MAX: f64 = 4.0;
/// TAU_in integration bounds: r = TAU_IN_R_MIN → TAU_IN_R_MAX on the b=0
/// equatorial line. The 1.05 floor excludes the launch-skin bins (v1 lesson).
const TAU_IN_R_MIN: f64 = 1.05;
const TAU_IN_R_MAX: f64 = B_MAX;

fn r_edge(i: usize) -> f64 {
    R_IN * (R_OUT / R_IN).powf(i as f64 / N_R as f64)
}

fn theta_edge(j: usize) -> f64 {
    (j as f64 / N_THETA as f64) * PI
}

fn r_bin_index(r: f64) -> usize {
    let rc = r.clamp(R_IN, R_OUT);
    let t = (rc / R_IN).ln() / (R_OUT / R_IN).ln();
    let idx = (t * N_R as f64) as isize;
    idx.clamp(0, N_R as isize - 1) as usize
}

fn theta_bin_index(theta: f64) -> usize {
    let tc = theta.clamp(0.0, PI);
    let idx = (tc / PI * N_THETA as f64) as isize;
    idx.clamp(0, N_THETA as isize - 1) as usize
}

fn grid_index(r: f64, theta: f64) -> usize {
    r_bin_index(r) * N_THETA + theta_bin_index(theta)
}

/// Cell volumes (spherical-shell segments, integrated over full azimuth),
/// computed once and cached — grid geometry is fixed for the module.
fn cell_volumes() -> &'static [f64] {
    static CELLS: OnceLock<Vec<f64>> = OnceLock::new();
    CELLS.get_or_init(|| {
        let mut v = vec![0.0f64; GRID_CELLS];
        for i in 0..N_R {
            let r1 = r_edge(i);
            let r2 = r_edge(i + 1);
            for j in 0..N_THETA {
                let th1 = theta_edge(j);
                let th2 = theta_edge(j + 1);
                // V = (2π/3)(r2³-r1³)(cosθ1 - cosθ2), θ1 < θ2 so cosθ1 > cosθ2.
                let vol = (2.0 * PI / 3.0) * (r2.powi(3) - r1.powi(3)) * (th1.cos() - th2.cos());
                v[i * N_THETA + j] = vol;
            }
        }
        v
    })
}

// ---------------------------------------------------------------------
// Analytic seed fields + TAU calibration (v2 fixes #2 and #3)
// ---------------------------------------------------------------------

/// Measured BARYON exit-structure peaks (project memory: 7°/58° dual peaks).
/// Two peaks because the baryon carries a FOUR-level spin stack (levels 9-12,
/// PHYSICS_REFERENCE §3), which traces two distinct exit families.
const EXIT_PEAK_LAT_DEG: [f64; 2] = [7.0, 58.0];
const EXIT_SIGMA_DEG: f64 = 5.0;

/// ELECTRON exit peak — a SINGLE equatorial disc. The electron at rest is
/// level 9: baryon-tier AXIAL spin only, no x/y/z (PHYSICS_REFERENCE §3b,
/// photon.html "the electron is spinning axially; the proton is spinning
/// axially plus x, y, and z"). One spin traces one exit family, so the bimodal
/// structure is unavailable to it by construction — not a fitted difference.
/// The peak location matches the measured trace histogram (peak at ±3°); the
/// core WIDTH is the input, and whether the engine regenerates the trace's
/// low mid-latitude plateau is the validation this profile is used to test.
const ELECTRON_EXIT_PEAK_LAT_DEG: [f64; 1] = [3.0];
const ELECTRON_EXIT_SIGMA_DEG: f64 = 2.5;

/// Unnormalized exit-latitude mixture pdf for `class`: equal-weight Gaussians
/// at ±each peak. Normalization is irrelevant — sigma_e is calibrated against
/// this same profile (see `tau_in_integral`), so only the SHAPE enters the
/// physics. The 1/(2·n_peaks) prefactor keeps the baryon value bit-identical
/// to the pre-generalization `0.25·(g(7)+g(−7)+g(58)+g(−58))`.
fn exit_lat_profile(class: ParticleClass, lat_deg: f64) -> f64 {
    let sigma = class.exit_sigma_deg();
    let two_s2 = 2.0 * sigma * sigma;
    let g = |mu: f64| (-((lat_deg - mu) * (lat_deg - mu)) / two_s2).exp();
    let peaks = class.exit_peaks_deg();
    // Accumulated ONE TERM AT A TIME, and divided rather than multiplied by the
    // reciprocal, so the two-peak sum reproduces the pre-generalization
    // expression `0.25·(g(7)+g(−7)+g(58)+g(−58))` bit-for-bit: the addition
    // association matches, `0.0 + x == x` exactly, and /4.0 == *0.25 exactly.
    // An ULP shift here would propagate through the sigma_e calibration and
    // could flip a collision coin, breaking the v2/v3/v4 regression anchors.
    let mut sum = 0.0;
    for &p in peaks {
        sum += g(p);
        sum += g(-p);
    }
    sum / (2.0 * peaks.len() as f64)
}

/// Analytic seed density for a population (v2 fix #2): emitted =
/// profile(θ)/r² (free-streaming bimodal source), ambient = 1 (uniform
/// far-field value in these units).
fn analytic_density(kind: PhotonKind, class: ParticleClass, r: f64, theta: f64) -> f64 {
    match kind {
        PhotonKind::Emitted => {
            let lat_deg = 90.0 - theta.to_degrees();
            exit_lat_profile(class, lat_deg) / (r * r).max(1.0)
        }
        PhotonKind::Ambient => 1.0,
    }
}

/// ∫ analytic n_e dr along the b=0 equatorial line from TAU_IN_R_MIN to
/// TAU_IN_R_MAX: profile(0°)·(1/r_min − 1/r_max). sigma_e = TAU / this.
fn tau_in_integral(class: ParticleClass) -> f64 {
    exit_lat_profile(class, 0.0) * (1.0 / TAU_IN_R_MIN - 1.0 / TAU_IN_R_MAX)
}

/// ∫ analytic n_a dr for an equatorial emitted photon r = 1 → R_OUT:
/// 1·(R_OUT − R_IN). sigma_a = TAU / this.
fn tau_out_integral() -> f64 {
    R_OUT - R_IN
}

/// ∫ analytic n_e dr for an equatorial EMITTED photon through its OWN
/// population, r = TAU_IN_R_MIN → R_OUT (skin excluded, v1 lesson):
/// profile(0°)·(1/TAU_IN_R_MIN − 1/R_OUT). sigma_ee = TAU_SELF / this.
fn tau_self_emitted_integral(class: ParticleClass) -> f64 {
    exit_lat_profile(class, 0.0) * (1.0 / TAU_IN_R_MIN - 1.0 / R_OUT)
}

/// ∫ analytic n_a dr for a b=0 AMBIENT photon through its own population
/// across the near zone where the self channel is active (see module doc,
/// v3 fix #3): 1·(B_MAX − R_IN). sigma_aa = TAU_SELF / this.
fn tau_self_ambient_integral() -> f64 {
    B_MAX - R_IN
}

// ---------------------------------------------------------------------
// Density field (track-length estimator)
// ---------------------------------------------------------------------

#[derive(Clone)]
struct DensityField {
    /// Σ weight*ds per cell → density = track_weight / cell_volume.
    track_weight: Vec<f64>,
    /// Σ weight*ds*dir per cell (raw Cartesian vector sum) — mean-direction /
    /// coherence source for collision partner sampling ONLY (see module doc
    /// deviation #2 for why this is not used for the circulation matrix).
    dir_accum: Vec<DVec3>,
    /// Σ weight*ds*(dir·r̂_local) per cell, r̂_local evaluated at each
    /// photon's own position — the circulation-matrix radial-flow source.
    radial_flux: Vec<f64>,
    /// Σ weight*ds*(dir·θ̂_local) per cell — meridional flow (v3). θ̂ points
    /// toward INCREASING colatitude (southward), so poleward drift shows as
    /// negative in the northern hemisphere and positive in the southern.
    theta_flux: Vec<f64>,
}

impl DensityField {
    fn new() -> Self {
        Self {
            track_weight: vec![0.0; GRID_CELLS],
            dir_accum: vec![DVec3::ZERO; GRID_CELLS],
            radial_flux: vec![0.0; GRID_CELLS],
            theta_flux: vec![0.0; GRID_CELLS],
        }
    }

    fn deposit(&mut self, pos: DVec3, dir: DVec3, weight: f64, ds: f64) {
        let r = pos.length().max(1e-9);
        let theta = colatitude_of(pos, r);
        let idx = grid_index(r, theta);
        let r_hat = pos / r;
        let phi = DVec3::Z.cross(r_hat);
        // θ̂ = φ̂ × r̂ (southward); on the polar axis φ̂ degenerates and the
        // meridional component is geometrically undefined — deposit 0 there.
        let theta_hat = if phi.length_squared() < 1e-18 {
            DVec3::ZERO
        } else {
            phi.normalize().cross(r_hat)
        };
        let w = weight * ds;
        self.track_weight[idx] += w;
        self.dir_accum[idx] += dir * w;
        self.radial_flux[idx] += w * dir.dot(r_hat);
        self.theta_flux[idx] += w * dir.dot(theta_hat);
    }

    fn density_by_idx(&self, idx: usize) -> f64 {
        self.track_weight[idx] / cell_volumes()[idx]
    }

    fn density_at(&self, r: f64, theta: f64) -> f64 {
        self.density_by_idx(grid_index(r, theta))
    }

    /// (mean unit direction, coherence = |Σ dir|/Σ weight ∈ [0,1]) at a cell.
    fn mean_dir_and_coherence_at(&self, r: f64, theta: f64) -> (DVec3, f64) {
        let idx = grid_index(r, theta);
        let w = self.track_weight[idx];
        if w < 1e-12 {
            return (DVec3::Z, 0.0);
        }
        let coherence = (self.dir_accum[idx].length() / w).min(1.0);
        let mean_dir = self.dir_accum[idx].normalize_or_zero();
        let mean_dir = if mean_dir.length_squared() > 1e-12 {
            mean_dir
        } else {
            DVec3::Z
        };
        (mean_dir, coherence)
    }

    fn relax(&self, measured: &DensityField, alpha: f64) -> DensityField {
        let mut out = DensityField::new();
        for idx in 0..GRID_CELLS {
            out.track_weight[idx] =
                alpha * self.track_weight[idx] + (1.0 - alpha) * measured.track_weight[idx];
            out.dir_accum[idx] =
                self.dir_accum[idx] * alpha + measured.dir_accum[idx] * (1.0 - alpha);
            out.radial_flux[idx] =
                alpha * self.radial_flux[idx] + (1.0 - alpha) * measured.radial_flux[idx];
            out.theta_flux[idx] =
                alpha * self.theta_flux[idx] + (1.0 - alpha) * measured.theta_flux[idx];
        }
        out
    }
}

// ---------------------------------------------------------------------
// Collision field: analytic seed × capped tally correction (v2 fix #2)
// ---------------------------------------------------------------------

/// The density field one population presents to the OTHER population's
/// transport: analytic seed times a per-cell multiplicative correction from
/// the (relaxed) tallies, capped to [CORRECTION_MIN, CORRECTION_MAX].
struct CollisionField {
    /// The population this field DESCRIBES (not the one marching through it).
    kind: PhotonKind,
    /// Emitting class — selects the analytic seed's exit-latitude profile
    /// (baryon bimodal disc vs electron single disc). Irrelevant for the
    /// ambient population, whose seed is uniform.
    class: ParticleClass,
    tally: DensityField,
    /// Global scale mapping tallied track-density into analytic units;
    /// 0.0 = no usable tally → pure analytic (correction ≡ 1).
    norm: f64,
}

impl CollisionField {
    fn pure_analytic(kind: PhotonKind, class: ParticleClass) -> Self {
        Self {
            kind,
            class,
            tally: DensityField::new(),
            norm: 0.0,
        }
    }

    /// Wraps a tally. The normalization region differs by population:
    /// emitted launches are tracked completely, so the whole domain
    /// participates; ambient tallies are only shape-complete for r ≤ B_MAX
    /// (see module doc, v2 fix #2), so both the norm and the correction are
    /// restricted to that near zone (correction ≡ 1 outside it).
    fn from_tally(kind: PhotonKind, class: ParticleClass, tally: DensityField) -> Self {
        let mut measured_total = 0.0;
        let mut analytic_total = 0.0;
        for i in 0..N_R {
            let r_c = (r_edge(i) * r_edge(i + 1)).sqrt();
            if kind == PhotonKind::Ambient && r_c > B_MAX {
                continue;
            }
            for j in 0..N_THETA {
                let idx = i * N_THETA + j;
                let th_c = 0.5 * (theta_edge(j) + theta_edge(j + 1));
                measured_total += tally.track_weight[idx];
                analytic_total += analytic_density(kind, class, r_c, th_c) * cell_volumes()[idx];
            }
        }
        let norm = if measured_total > 1e-12 {
            analytic_total / measured_total
        } else {
            0.0
        };
        Self {
            kind,
            class,
            tally,
            norm,
        }
    }

    fn effective_density_at(&self, r: f64, theta: f64) -> f64 {
        let a = analytic_density(self.kind, self.class, r, theta);
        if self.norm == 0.0 {
            return a;
        }
        if self.kind == PhotonKind::Ambient && r > B_MAX {
            return a; // outside the shape-complete near zone: correction = 1
        }
        let measured = self.tally.density_at(r, theta) * self.norm;
        a * (measured / a.max(1e-300)).clamp(CORRECTION_MIN, CORRECTION_MAX)
    }
}

/// Effective (regularized) density snapshot at cell centers, for the
/// convergence metric.
fn effective_density_snapshot(
    kind: PhotonKind,
    class: ParticleClass,
    tally: &DensityField,
) -> Vec<f64> {
    let field = CollisionField::from_tally(kind, class, tally.clone());
    let mut out = vec![0.0; GRID_CELLS];
    for i in 0..N_R {
        let r_c = (r_edge(i) * r_edge(i + 1)).sqrt();
        for j in 0..N_THETA {
            let th_c = 0.5 * (theta_edge(j) + theta_edge(j + 1));
            out[i * N_THETA + j] = field.effective_density_at(r_c, th_c);
        }
    }
    out
}

fn l2_rel_change(old: &[f64], new: &[f64]) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for idx in 0..old.len() {
        num += (new[idx] - old[idx]).powi(2);
        den += old[idx].powi(2);
    }
    if den > 1e-300 {
        num.sqrt() / den.sqrt()
    } else {
        0.0
    }
}

// ---------------------------------------------------------------------
// Photons
// ---------------------------------------------------------------------

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum PhotonKind {
    Emitted,
    Ambient,
}

/// Particle spin class for the recycling test.
///
/// BARYONS (proton/neutron): identical emission GEOMETRY — every baryon-sized
/// 4-level construction traces the same bimodal equatorial disc — and the ONLY
/// difference is the level-12 (z) spin sign, the survey's discriminator (commit
/// 8b0d374: sub-c L12 settles proton outer_spin positive, neutron negative).
/// That sign enters as (a) the emitted charge's chirality tag and (b) the sign
/// of the outer-spin swing. Whether it produces a distinct exit topology is the
/// emergent question the photon-gas run answers — and it does (commit 4e68dee).
///
/// ELECTRON: structurally NOT a small baryon. It is level 9 — baryon-tier
/// AXIAL spin only, no x/y/z (PHYSICS_REFERENCE §3b) — so it has no level-12
/// z-spin at all. Two consequences, both structural rather than fitted:
/// (a) its chirality tag comes from the AXIAL spin sign instead of a z-sign,
/// and (b) with one spin rather than four it traces a SINGLE-peak disc, so the
/// bimodal exit structure is unavailable to it by construction.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum ParticleClass {
    Proton,
    Neutron,
    Electron,
}

impl ParticleClass {
    /// Sign of the spin that tags this class's emitted charge, and of its
    /// outer-spin swing. For baryons this is the level-12 z-spin sign (proton
    /// +1, neutron −1). The electron HAS no z-spin — the sign here is its
    /// axial (level-9) spin sense, which is an orientation rather than a class
    /// invariant; +1 is the reference orientation.
    fn chirality_sign(self) -> f64 {
        match self {
            ParticleClass::Proton => 1.0,
            ParticleClass::Neutron => -1.0,
            ParticleClass::Electron => 1.0,
        }
    }

    /// Exit-latitude peak magnitudes (degrees from the equator).
    fn exit_peaks_deg(self) -> &'static [f64] {
        match self {
            ParticleClass::Proton | ParticleClass::Neutron => &EXIT_PEAK_LAT_DEG,
            ParticleClass::Electron => &ELECTRON_EXIT_PEAK_LAT_DEG,
        }
    }

    fn exit_sigma_deg(self) -> f64 {
        match self {
            ParticleClass::Proton | ParticleClass::Neutron => EXIT_SIGMA_DEG,
            ParticleClass::Electron => ELECTRON_EXIT_SIGMA_DEG,
        }
    }

    fn label(self) -> &'static str {
        match self {
            ParticleClass::Proton => "proton",
            ParticleClass::Neutron => "neutron",
            ParticleClass::Electron => "electron",
        }
    }
}

/// Chirality-gear collision physics (PHYSICS_REFERENCE §5). `enabled=false`
/// is a strict no-op: no chirality is drawn or read, so the RNG stream and
/// every result are bit-identical to the pre-chirality module (preserves the
/// TAU_SELF=0 v2 anchor and all existing report/test seeds). When enabled,
/// emitted charge carries χ = class.chirality_sign(), ambient partners are sampled at
/// `photon_fraction`, and an opposite-chirality (cancel) collision reverses
/// the outward radial component by `cancel_backscatter` (turn-around).
#[derive(Clone, Copy)]
struct GearConfig {
    enabled: bool,
    photon_fraction: f64,
    /// Through-channel (axial funnel) strength — see DEFAULT_CANCEL_REDIRECT.
    cancel_redirect: f64,
    /// Species-split directed intake strength (venus2.pdf / neutron.pdf):
    /// 0 = isotropic ambient (the original ensemble); as this rises, ambient
    /// photons are steered along the pole axis by species (photon → enters
    /// south, antiphoton → enters north), building the poleward intake flow
    /// the isotropic gas lacked. This is the coupling that closes the pole
    /// limb; swept in the report.
    pole_intake_bias: f64,
    class: ParticleClass,
    /// Per-collision spin transfer for the MAGNETIC observable (magmom.pdf:
    /// "if we have a collision between an emitted photon and an ambient photon
    /// above the pole of the neutron, those spins can stack or cancel").
    /// 0 = spin is a frozen tag (the pre-existing behaviour, where the escaped
    /// chirality sum could not show cancellation at all). See `Photon.spin`.
    spin_transfer: f64,
    /// Through-charge lane radius (fraction of body radius). 0 = limb OFF, which
    /// is bit-identical to the pre-through-charge module: the branch consumes no
    /// RNG draws, so the stream is untouched. Default `THROUGH_LANE_RHO` is the
    /// measured hole, not a fitted value.
    through_lane_rho: f64,
}

impl GearConfig {
    /// The regression-preserving no-op (class Proton is irrelevant while
    /// disabled — chirality_sign is +1 so even the swing sign is unchanged).
    fn off() -> Self {
        Self {
            enabled: false,
            photon_fraction: PHOTON_FRACTION,
            cancel_redirect: 0.0,
            pole_intake_bias: 0.0,
            class: ParticleClass::Proton,
            spin_transfer: 0.0,
            through_lane_rho: 0.0,
        }
    }

    /// Chirality of a partner drawn from `kind`'s population: ambient is a
    /// 2/3-photon mix (sampled), the emitted population is monochiral at the
    /// class z-sign (a particle's own emitted charge is all one chirality, so
    /// emitted↔emitted self-collisions always augment).
    fn partner_chirality(&self, kind: PhotonKind, rng: &mut u64) -> f64 {
        match kind {
            PhotonKind::Ambient => {
                if xorshift64(rng) < self.photon_fraction {
                    1.0
                } else {
                    -1.0
                }
            }
            PhotonKind::Emitted => self.class.chirality_sign(),
        }
    }
}

#[derive(Clone, Copy, Debug)]
struct Photon {
    pos: DVec3,
    dir: DVec3,
    weight: f64,
    kind: PhotonKind,
    bounces: u32,
    /// Emitted charge's spin chirality (±1 = class chirality sign); ambient
    /// photons carry 0 here (their partner chirality is sampled at collision
    /// time, not stored). Unused while gears are disabled. This is a FROZEN
    /// TAG — it identifies which population the photon belongs to and gates
    /// the mesh, and it must not be mutated (the gear geometry keys off it).
    chirality: f64,
    /// Carried SPIN — the magnetic degree of freedom, mutable by collisions.
    /// Initialized to `chirality` and then stacked/cancelled per collision at
    /// `GearConfig.spin_transfer`.
    ///
    /// This is the split magmom.pdf insists on: magnetism is the photon's SPIN
    /// and "those spins can stack or cancel" on collision, whereas electricity
    /// is the linear motion c and "photon collisions affect only the spins, not
    /// c" (neutron.pdf). So the magnetic output degrades along the path while
    /// the electrical output (the count, and the linear-momentum vector sum)
    /// does not. Summing `chirality` at escape — as the old chi_escaped_sum
    /// did — could never show cancellation, because the tag never changed.
    spin: f64,
}

enum Termination {
    Absorbed { pos: DVec3, dir: DVec3 },
    Escaped { pos: DVec3, dir: DVec3 },
    Lost,
}

#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Band {
    Polar = 0,
    Mid = 1,
    Equatorial = 2,
}

fn lat_band(lat_deg: f64) -> Band {
    let a = lat_deg.abs();
    if a > 60.0 {
        Band::Polar
    } else if a >= 30.0 {
        Band::Mid
    } else {
        Band::Equatorial
    }
}

fn colatitude_of(pos: DVec3, r: f64) -> f64 {
    (pos.z / r).clamp(-1.0, 1.0).acos()
}

fn arrival_latitude(pos: DVec3) -> f64 {
    let r = pos.length().max(1e-9);
    90.0 - colatitude_of(pos, r).to_degrees()
}

fn latlon_to_unit(lat_deg: f64, az: f64) -> DVec3 {
    let lat = lat_deg.to_radians();
    let colat = FRAC_PI_2 - lat;
    let s = colat.sin();
    let c = colat.cos();
    DVec3::new(s * az.cos(), s * az.sin(), c)
}

fn phi_hat_at(r_hat: DVec3) -> DVec3 {
    let v = DVec3::Z.cross(r_hat);
    if v.length_squared() < 1e-18 {
        DVec3::X
    } else {
        v.normalize()
    }
}

fn rand_gaussian(rng: &mut u64, mean: f64, sigma: f64) -> f64 {
    let u1 = xorshift64(rng).max(1e-12);
    let u2 = xorshift64(rng);
    let z0 = (-2.0 * u1.ln()).sqrt() * (TAU * u2).cos();
    mean + sigma * z0
}

/// Draws a peak MAGNITUDE (0..89, unsigned) from `class`'s equal-weight peak
/// mixture. Sign/hemisphere is applied by the caller — see module doc
/// deviation #1 for why sign is shared rather than redrawn per latitude.
///
/// The selector draw is consumed for EVERY class, including single-peak ones
/// where it carries no information, so that the two-peak (baryon) RNG stream
/// stays bit-identical to the pre-generalization sampler: for n=2 the index
/// `(u·2) as usize` reproduces the old `if u < 0.5` branch exactly.
fn sample_exit_magnitude_deg(rng: &mut u64, class: ParticleClass) -> f64 {
    let peaks = class.exit_peaks_deg();
    let u = xorshift64(rng);
    let idx = ((u * peaks.len() as f64) as usize).min(peaks.len() - 1);
    rand_gaussian(rng, peaks[idx], class.exit_sigma_deg()).clamp(0.0, 89.0)
}

fn spawn_emitted(rng: &mut u64, swing_boost: f64, class: ParticleClass) -> Photon {
    let sign = if xorshift64(rng) < 0.5 { 1.0 } else { -1.0 };
    let lat0 = sign * sample_exit_magnitude_deg(rng, class);
    let az0 = xorshift64(rng) * TAU;
    let r_hat0 = latlon_to_unit(lat0, az0);
    let pos = r_hat0 * 1.001;

    // NEW draw for the exit direction's own latitude (same hemisphere sign,
    // independent peak/magnitude draw — module doc deviation #1).
    let lat1 = sign * sample_exit_magnitude_deg(rng, class);
    let dir_radial = latlon_to_unit(lat1, az0);
    let phi_hat = phi_hat_at(dir_radial);
    // Swing sense follows the settled outer_spin sign, which the survey ties
    // to the z-sign (proton +, neutron −). For a proton (sign=+1) this is
    // bit-identical to the pre-chirality path (same RNG draws, same sign).
    let swing = swing_boost * class.chirality_sign();
    let dir = (dir_radial + phi_hat * swing).normalize();

    Photon {
        pos,
        dir,
        weight: 1.0,
        kind: PhotonKind::Emitted,
        bounces: 0,
        chirality: class.chirality_sign(),
        spin: class.chirality_sign(),
    }
}

/// Latitude of the point where a straight ray from `pos` in direction `dir`
/// would cross r=R_IN, or `None` if it misses the sphere entirely.
fn ray_sphere_aim_latitude(pos: DVec3, dir: DVec3) -> Option<f64> {
    let b = 2.0 * pos.dot(dir);
    let c = pos.dot(pos) - R_IN * R_IN;
    let disc = b * b - 4.0 * c;
    if disc < 0.0 {
        return None;
    }
    let sq = disc.sqrt();
    let t1 = (-b - sq) / 2.0;
    let t2 = (-b + sq) / 2.0;
    let t = if t1 > 1e-9 {
        t1
    } else if t2 > 1e-9 {
        t2
    } else {
        return None;
    };
    let hit = pos + dir * t;
    Some(arrival_latitude(hit))
}

/// v2 fix #1: importance-sampled ambient launch. The straight-line impact
/// parameter b (closest approach to the origin) is area-weighted over
/// [0, B_MAX]: b = B_MAX·sqrt(u) — proportional to the b·db flux weighting
/// of a uniform isotropic field restricted to b ≤ B_MAX, so all photons
/// carry equal weight and the ensemble is the physically-correct
/// CONDITIONAL near-zone flux. Construction per spec: random direction d
/// uniform on the sphere, random unit e ⊥ d, line start = −d·R_OUT + e·b,
/// direction d; the photon is then advanced along d to its entry into the
/// tallied shell (r = R_OUT), where the march begins.
///
/// Species-split directed intake (gears on, pole_intake_bias > 0): the
/// ambient photon is first assigned a chirality (±1 at photon_fraction), then
/// its travel direction is tilted along the pole axis by species — a photon
/// (χ=+1) travels +z so it enters the SOUTH pole; an antiphoton (χ=−1)
/// travels −z so it enters the NORTH pole (venus2.pdf / neutron.pdf). This
/// breaks the isotropic near-zone ensemble on purpose: it is the directed
/// intake that builds the poleward flow the isotropic gas never had. With
/// bias = 0 (or gears off) the launch is the original isotropic one and the
/// RNG stream for the gears-off path is unchanged (no chirality is drawn).
fn spawn_ambient(rng: &mut u64, gears: &GearConfig) -> (Photon, Option<f64>) {
    let (chirality, d) = if gears.enabled {
        let chi = if xorshift64(rng) < gears.photon_fraction {
            1.0
        } else {
            -1.0
        };
        let base = rand_unit_vec(rng);
        let d = if gears.pole_intake_bias > 0.0 {
            (base + DVec3::Z * (gears.pole_intake_bias * chi)).normalize()
        } else {
            base
        };
        (chi, d)
    } else {
        (0.0, rand_unit_vec(rng))
    };
    let helper = if d.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let t_perp = d.cross(helper).normalize();
    let b_perp = d.cross(t_perp);
    let psi = xorshift64(rng) * TAU;
    let e = t_perp * psi.cos() + b_perp * psi.sin();
    let b = B_MAX * xorshift64(rng).sqrt();

    let start = -d * R_OUT + e * b;
    // |start + t·d|² = b² + (t − R_OUT)² ⇒ shell entry at
    // t = R_OUT − sqrt(R_OUT² − b²); nudge ε further along d so the march's
    // escape check (r > R_OUT) doesn't fire on the launch point itself.
    let t_entry = R_OUT - (R_OUT * R_OUT - b * b).sqrt();
    let pos = start + d * (t_entry + 1e-9);
    let aim_lat = ray_sphere_aim_latitude(pos, d);
    (
        Photon {
            pos,
            dir: d,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality,
            // An ambient photon carries its own species spin (0 when gears are
            // off, where no chirality is drawn). This is the counter-stream
            // whose spin cancels the neutron's polar exhaust.
            spin: chirality,
        },
        aim_lat,
    )
}

/// Samples a contact normal for a glancing hard-sphere exchange: `b` is a
/// uniformly-random unit vector perpendicular to `rel_hat`; `u` is biased
/// toward 1 (grazing) rather than uniform, per arp.pdf ("photon-photon hits
/// are glancing, almost never head-on").
fn sample_contact_normal(rel_hat: DVec3, rng: &mut u64) -> DVec3 {
    let helper = if rel_hat.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let t = rel_hat.cross(helper).normalize();
    let bperp = rel_hat.cross(t);
    let psi = xorshift64(rng) * TAU;
    let b = t * psi.cos() + bperp * psi.sin();
    let u0 = xorshift64(rng);
    let u = 1.0 - (1.0 - u0).powi(3); // skew toward u=1 (grazing)
    let nhat = b * u.sqrt() + rel_hat * (1.0 - u).sqrt();
    if nhat.length_squared() > 1e-12 {
        nhat.normalize()
    } else {
        rel_hat
    }
}

/// One glancing hard-sphere exchange against a partner sampled from
/// `field`'s local mean direction + coherence at (r, theta). Returns the
/// photon's new direction (unchanged when the pair is exactly co-moving or
/// the exchange degenerates — co-moving photons don't jostle).
///
/// Chirality gears (PHYSICS_REFERENCE §5): when `gears.enabled`, a partner
/// chirality is drawn from `field`'s population and meshed with the marching
/// photon's `chirality`. Same-chirality (augment, mesh>0) is the soft glance
/// above — the photon keeps its outward progress and tends to escape.
/// Opposite-chirality (cancel, mesh<0) is a hard bounce: the outward radial
/// component of the exit direction is partially reversed (the blocked charge
/// "has to turn around"), sending the photon back into the near zone where
/// the pressure gradient — not this rule — decides where it re-emerges.
/// When disabled, no chirality is drawn (RNG stream preserved) and behavior
/// is bit-identical to the pre-chirality module.
#[allow(clippy::too_many_arguments)]
fn scatter_off(
    field: &CollisionField,
    dir: DVec3,
    pos: DVec3,
    r: f64,
    theta: f64,
    chirality: f64,
    gears: &GearConfig,
    rng: &mut u64,
) -> (DVec3, f64) {
    let (mean_dir, coherence) = field.tally.mean_dir_and_coherence_at(r, theta);
    let mut partner_dir = mean_dir * coherence + rand_unit_vec(rng);
    if partner_dir.length_squared() < 1e-12 {
        partner_dir = rand_unit_vec(rng);
    }
    partner_dir = partner_dir.normalize();

    let mut new_dir = dir;
    let rel = dir - partner_dir;
    if rel.length_squared() > 1e-12 {
        let rel_hat = rel.normalize();
        for _ in 0..4 {
            let nhat = sample_contact_normal(rel_hat, rng);
            let candidate = dir - dir.dot(nhat) * nhat + partner_dir.dot(nhat) * nhat;
            if candidate.length_squared() > 1e-8 {
                new_dir = candidate.normalize();
                break;
            }
        }
    }

    if !gears.enabled {
        return (new_dir, 0.0);
    }

    // Chirality mesh. The partner-chirality draw only happens here (gears on),
    // so the disabled path never consumes it.
    let chi_partner = gears.partner_chirality(field.kind, rng);
    let mesh = chirality * chi_partner;

    // MAGNETIC channel (magmom.pdf: "if we have a collision between an emitted
    // photon and an ambient photon above the pole of the neutron, those spins
    // can stack or cancel"). One rule covers both: the photon picks up the
    // PARTNER's spin. Same sign stacks, opposite sign cancels — no case split
    // and no sign convention to get wrong. voyag.pdf gives the ensemble
    // consequence: "the spins will offset as a sum, and the field will be
    // magnetically flat" when the populations balance, with "a leftover or
    // resultant field spin" when they don't.
    //
    // NOTE this is the LINEAR motion's counterpart, not a replacement for it:
    // "photon collisions affect only the spins, not c" (neutron.pdf), so `dir`
    // and the escaped COUNT — the electrical output — are untouched here.
    let spin_delta = gears.spin_transfer * chi_partner;
    if gears.cancel_redirect > 0.0 {
        // SYMMETRIC gear channel (PHYSICS_REFERENCE §5: "same-spin deflects
        // ALONG the spin direction; opposite-spin deflects AGAINST it").
        // Both nudges are small per-collision; the NET over many collisions
        // decides the exit, so an augment-dominated (proton) charge channels
        // to the disc while a cancel-dominated (neutron) charge channels to
        // the poles. A one-sided (cancel-only) redirect over-funnels BOTH
        // classes polar because cancels compound over the ricochet chain;
        // the augment counter-nudge is what keeps the proton a disc.
        if mesh < 0.0 {
            // Cancel = blocked from the equator → through-channel to the
            // nearer pole (neutron.pdf "funnelled back to the central spin
            // axis… escape out one of the poles").
            let pole = if pos.z >= 0.0 { DVec3::Z } else { -DVec3::Z };
            let blended = new_dir.lerp(pole, gears.cancel_redirect);
            if blended.length_squared() > 1e-12 {
                return (blended.normalize(), spin_delta);
            }
        } else if mesh > 0.0 {
            // Augment = escapes normally → channel toward the equatorial disc
            // (the Faraday-sprinkler limb): flatten the exit toward the xy
            // plane.
            let flat = DVec3::new(new_dir.x, new_dir.y, 0.0);
            if flat.length_squared() > 1e-9 {
                let blended = new_dir.lerp(flat.normalize(), gears.cancel_redirect);
                if blended.length_squared() > 1e-12 {
                    return (blended.normalize(), spin_delta);
                }
            }
        }
    }
    (new_dir, spin_delta)
}

/// Ray-marches one photon through TWO regularized collision fields — the
/// other population's (`cross`, coupling `sigma_cross`) and its own
/// population's (`selff`, coupling `sigma_self`; v3) — depositing its own
/// track-length into `own_field`, until absorbed (r<R_IN), escaped
/// (r>R_OUT), or lost (safety path-length cap). The two channels' rates add
/// per step; on a collision the partner population is chosen proportional
/// to its local rate. Returns the termination and the photon's final
/// (cross, self) collision counts; `n_ricochet` caps their SUM.
///
/// With sigma_self = 0 the RNG draw sequence is identical to v2's
/// single-channel march (no extra draws are consumed), so TAU_SELF = 0
/// sweep cells reproduce v2 bit-for-bit at equal seeds.
#[allow(clippy::too_many_arguments)]
fn march(
    mut photon: Photon,
    cross: &CollisionField,
    sigma_cross: f64,
    selff: &CollisionField,
    sigma_self: f64,
    n_ricochet: u32,
    gears: &GearConfig,
    rng: &mut u64,
    own_field: &mut DensityField,
) -> (Termination, u32, u32, f64, u32) {
    let mut path_length = 0.0;
    let mut cross_bounces = 0u32;
    let mut self_bounces = 0u32;
    let mut through_count = 0u32;
    let term = loop {
        let r = photon.pos.length();
        if r < R_IN {
            // THROUGH-CHARGE LIMB (venus2.pdf "through charge... goes straight
            // through the body from pole to pole, avoiding lateral recycling";
            // neutron.pdf "the photons will go in the south pole and the
            // antiphotons will go in the north. Then they go out the other
            // pole"). An AMBIENT photon that reaches the surface inside the
            // polar pass-through lane, travelling near-axially, is NOT absorbed:
            // it crosses and exits the far side RETAINING ITS SPECIES AND SPIN.
            //
            // That retention is the whole point. Everything the pump re-emits is
            // monochiral charge from the fixed disc, which is why balanced-field
            // magnetic neutrality was unreachable: with no species-preserving
            // exit limb there was nothing for the counter-stream to cancel
            // against. Through-charge supplies the mixed-species exhaust.
            //
            // Selection is GEOMETRIC, per salt.pdf: "the charge that enters
            // nearest the center of the hole or pole [is] most likely to pass
            // through. Charge that comes in nearer the edges, or that enters on
            // an angle, will be forced by centrifugal forces into the equatorial
            // whirlpool." So the two gates are lane radius and axial incidence,
            // and the through-RATE is emergent from the arriving flux geometry —
            // only the lane WIDTH is imposed, and that is anchored to the
            // measured hole (see THROUGH_LANE_RHO).
            if photon.kind == PhotonKind::Ambient && gears.through_lane_rho > 0.0 {
                let rho = (photon.pos.x * photon.pos.x + photon.pos.y * photon.pos.y).sqrt();
                let pos_hat = photon.pos / r.max(1e-12);
                let incoming = photon.dir.dot(pos_hat); // < 0 = heading inward
                let axial = photon.dir.z.abs();
                if rho < gears.through_lane_rho * r
                    && axial > THROUGH_COS_MAX_ANGLE
                    && incoming < 0.0
                {
                    // Straight-line crossing: from a point at radius r, the far
                    // intersection with the same sphere is at t = −2(pos·dir).
                    // No teleport — this IS "goes straight through".
                    let t = -2.0 * photon.pos.dot(photon.dir);
                    if t > 0.0 {
                        let exit = photon.pos + photon.dir * t;
                        // Nudge just outside so the r < R_IN test doesn't refire
                        // on the exit point itself.
                        photon.pos = exit * (R_IN / exit.length().max(1e-12)) * 1.001;
                        through_count += 1;
                        path_length += t;
                        // Chirality and spin are deliberately untouched. The
                        // exiting stream now meets the incoming stream at the
                        // far pole and their spins mesh through the normal
                        // collision path — neutron.pdf's "as they leave, they
                        // collide with photons coming in, and we get spin
                        // cancellations" is then emergent, not imposed.
                        continue;
                    }
                }
            }
            break Termination::Absorbed {
                pos: photon.pos,
                dir: photon.dir,
            };
        }
        if r > R_OUT {
            break Termination::Escaped {
                pos: photon.pos,
                dir: photon.dir,
            };
        }
        if path_length > SAFETY_PATH_CAP {
            break Termination::Lost;
        }

        let ds = (DS_FRACTION * r).max(MIN_DS);
        let theta = colatitude_of(photon.pos, r);
        own_field.deposit(photon.pos, photon.dir, photon.weight, ds);

        if (sigma_cross > 0.0 || sigma_self > 0.0) && photon.bounces < n_ricochet {
            let lam_cross = if sigma_cross > 0.0 {
                sigma_cross * cross.effective_density_at(r, theta)
            } else {
                0.0
            };
            // Ambient-ambient far-zone gate (module doc, v3 fix #3): a
            // uniform isotropic gas is invariant under self-scattering, so
            // outside the near zone the self channel carries no signal and
            // is switched off to preserve the aim bookkeeping.
            let lam_self = if sigma_self > 0.0
                && !(selff.kind == PhotonKind::Ambient && r > B_MAX)
            {
                sigma_self * selff.effective_density_at(r, theta)
            } else {
                0.0
            };
            let lam = lam_cross + lam_self;
            if lam > 0.0 {
                let p_collide = 1.0 - (-lam * ds).exp();
                if xorshift64(rng) < p_collide {
                    // Channel pick draws randomness ONLY when both channels
                    // compete (keeps the sigma_self=0 stream v2-identical).
                    let use_self = if lam_cross <= 0.0 {
                        true
                    } else if lam_self <= 0.0 {
                        false
                    } else {
                        xorshift64(rng) < lam_self / lam
                    };
                    let field = if use_self { selff } else { cross };
                    let (new_dir, spin_delta) = scatter_off(
                        field,
                        photon.dir,
                        photon.pos,
                        r,
                        theta,
                        photon.chirality,
                        gears,
                        rng,
                    );
                    photon.dir = new_dir;
                    // Spins stack or cancel; |spin| is capped at the fully
                    // stacked double-spin value so a long ricochet chain can't
                    // manufacture unbounded magnetism.
                    photon.spin = (photon.spin + spin_delta).clamp(-2.0, 2.0);
                    photon.bounces += 1;
                    if use_self {
                        self_bounces += 1;
                    } else {
                        cross_bounces += 1;
                    }
                }
            }
        }

        photon.pos += photon.dir * ds;
        path_length += ds;
    };
    // The photon is taken by value, so its collision-mutated SPIN has to come
    // back out here — this is the magnetic state the tally records.
    (
        term,
        cross_bounces,
        self_bounces,
        photon.spin,
        through_count,
    )
}

// ---------------------------------------------------------------------
// Tallies
// ---------------------------------------------------------------------

#[derive(Default, Clone)]
struct PassTally {
    total: f64,
    arrived: f64,
    arrived_band: [f64; 3],
    escaped: f64,
    /// Escaped-photon count by latitude band of the escape point — the exit
    /// TOPOLOGY. For emitted charge this is the emission profile the gas
    /// actually produces: equatorial-dominated = proton disc, polar-dominated
    /// = the neutron's axial channel (if it emerges).
    escaped_band: [f64; 3],
    /// Σ chirality over escaped EMITTED photons. NOTE: `chirality` is a frozen
    /// population tag, so for a monochiral emitted population this is just
    /// class_sign × escaped and carries NO cancellation information. Kept only
    /// as a bookkeeping cross-check on `escaped`; the real magnetic observable
    /// is `spin_escaped_sum` below.
    chi_escaped_sum: f64,
    /// Σ carried SPIN over escaped EMITTED photons — THE MAGNETIC OUTPUT.
    /// Magnetism is the spin sum and is geometry-blind (magmom.pdf: spins "can
    /// add in immediately" whatever direction the photon left in), so unlike
    /// the electrical vector sum below this does not care where on the surface
    /// the charge exited. Neutron neutrality shows up here as a collapse toward
    /// 0 while `escaped` stays large — charge still recycled, magnetism
    /// cancelled (neutron.pdf).
    spin_escaped_sum: f64,
    /// Σ carried spin over escaped emitted photons, split by escape band —
    /// locates WHERE the surviving magnetism leaves.
    spin_escaped_band: [f64; 3],
    /// Σ carried spin over ABSORBED (incoming) photons, by arrival band. This
    /// is the counter-stream: where it shares a channel with the escaping
    /// exhaust, the two spin sums cancel. The proton intakes at the poles and
    /// exhausts at the equator (different channels → no cancellation); the
    /// neutron does both on the pole axis (same channel → cancellation). That
    /// channel overlap is the whole proton/neutron magnetic asymmetry, and it
    /// is GEOMETRIC — no free parameter sets it.
    spin_arrived_band: [f64; 3],
    /// Σ of escaped-photon DIRECTION vectors — the electrical output, which is
    /// a genuine VECTOR sum of linear motion (magmom.pdf: the neutron "is
    /// always emitting at a right angle to the proton, so if the charge field
    /// is defined in terms of proton emission, the neutron's charge energy will
    /// seem to be zero"). Resolved as xy-plane magnitude vs |z| by the report:
    /// a proton disc sums into xy, a neutron axial exhaust sums into z and so
    /// contributes ~nothing to a proton-defined equatorial charge field.
    /// Unaffected by collisions, which "affect only the spins, not c".
    dir_escaped_sum: DVec3,
    /// Σ |xy| and Σ |z| of escaped directions — the incoherent (magnitude)
    /// counterpart, so a cancelling vector sum can be told apart from an
    /// absence of flux.
    dir_escaped_xy_abs: f64,
    dir_escaped_z_abs: f64,
    /// THROUGH-CHARGE (ambient population only): photons that crossed the body
    /// pole-to-pole through the polar lane instead of being absorbed, and their
    /// carried spin at escape. This is the SPECIES-PRESERVING output limb — the
    /// pump's disc re-emission is monochiral, so this is the only channel whose
    /// exhaust can be a mixture and therefore the only one whose spin sum can
    /// cancel. `through_crossings` counts crossings (a photon may cross more
    /// than once); `through_photons` counts distinct photons that crossed.
    through_crossings: f64,
    through_photons: f64,
    /// Σ carried spin over escaped photons that made at least one crossing, and
    /// their escape-band split. Neutron neutrality should show up as this sum
    /// collapsing toward 0 in a BALANCED field while the count stays large.
    through_spin_escaped: f64,
    through_escaped_band: [f64; 3],
    /// Fine escape-latitude histogram of emitted charge: 180 one-degree bins,
    /// index = floor(lat_deg + 90) ∈ [0,179] (0 = −90° south pole, 179 = +89°
    /// north). Empty until the first escape (sized lazily so PassTally keeps
    /// deriving Default). This is the atom-mode emission profile the recycling
    /// engine exports — the `degree,percentage` histogram Atom mode loads.
    escaped_lat_hist: Vec<f64>,
    lost: f64,
    lz_arrived: f64,
    lz_escaped: f64,
    /// Σ CROSS-channel collision events over ALL photons — knob verification.
    collisions_cross_sum: f64,
    /// Σ SELF-channel collision events over ALL photons (v3) — knob
    /// verification for TAU_SELF.
    collisions_self_sum: f64,
    /// Σ bounces (both channels) over ABSORBED photons — "are arrivals the
    /// multiply-scattered ones?"
    arrived_bounce_sum: f64,
    // ambient-only:
    aim_eq_count: f64,
    aim_eq_arrived: f64,
    aim_eq_arrived_polar: f64,
    aim_miss_count: f64,
    aim_miss_arrived: f64,
    aim_miss_arrived_band: [f64; 3],
}

impl PassTally {
    fn record_emitted(
        &mut self,
        term: &Termination,
        cross_bounces: u32,
        self_bounces: u32,
        chirality: f64,
        spin: f64,
    ) {
        let bounces = cross_bounces + self_bounces;
        self.total += 1.0;
        self.collisions_cross_sum += cross_bounces as f64;
        self.collisions_self_sum += self_bounces as f64;
        match term {
            Termination::Absorbed { pos, dir } => {
                self.arrived += 1.0;
                self.arrived_bounce_sum += bounces as f64;
                let band = lat_band(arrival_latitude(*pos));
                self.arrived_band[band as usize] += 1.0;
                self.spin_arrived_band[band as usize] += spin;
                self.lz_arrived += pos.cross(*dir).z;
            }
            Termination::Escaped { pos, dir } => {
                self.escaped += 1.0;
                let lat = arrival_latitude(*pos);
                let band = lat_band(lat);
                self.escaped_band[band as usize] += 1.0;
                self.chi_escaped_sum += chirality;
                self.spin_escaped_sum += spin;
                self.spin_escaped_band[band as usize] += spin;
                self.dir_escaped_sum += *dir;
                self.dir_escaped_xy_abs += (dir.x * dir.x + dir.y * dir.y).sqrt();
                self.dir_escaped_z_abs += dir.z.abs();
                self.lz_escaped += pos.cross(*dir).z;
                if self.escaped_lat_hist.is_empty() {
                    self.escaped_lat_hist = vec![0.0; 180];
                }
                let idx = ((lat + 90.0).floor() as isize).clamp(0, 179) as usize;
                self.escaped_lat_hist[idx] += 1.0;
            }
            Termination::Lost => self.lost += 1.0,
        }
    }

    fn record_ambient(
        &mut self,
        term: &Termination,
        cross_bounces: u32,
        self_bounces: u32,
        aim_lat: Option<f64>,
        spin: f64,
        through: u32,
    ) {
        let bounces = cross_bounces + self_bounces;
        self.total += 1.0;
        self.collisions_cross_sum += cross_bounces as f64;
        self.collisions_self_sum += self_bounces as f64;
        self.through_crossings += through as f64;
        if through > 0 {
            self.through_photons += 1.0;
        }
        let aim_band = aim_lat.map(lat_band);
        match aim_band {
            Some(Band::Equatorial) => self.aim_eq_count += 1.0,
            None => self.aim_miss_count += 1.0,
            _ => {}
        }
        match term {
            Termination::Absorbed { pos, dir } => {
                self.arrived += 1.0;
                self.arrived_bounce_sum += bounces as f64;
                let band = lat_band(arrival_latitude(*pos));
                self.arrived_band[band as usize] += 1.0;
                self.spin_arrived_band[band as usize] += spin;
                self.lz_arrived += pos.cross(*dir).z;
                match aim_band {
                    Some(Band::Equatorial) => {
                        self.aim_eq_arrived += 1.0;
                        if band == Band::Polar {
                            self.aim_eq_arrived_polar += 1.0;
                        }
                    }
                    None => {
                        self.aim_miss_arrived += 1.0;
                        self.aim_miss_arrived_band[band as usize] += 1.0;
                    }
                    _ => {}
                }
            }
            Termination::Escaped { pos, dir } => {
                self.escaped += 1.0;
                let band = lat_band(arrival_latitude(*pos));
                self.escaped_band[band as usize] += 1.0;
                self.spin_escaped_band[band as usize] += spin;
                self.lz_escaped += pos.cross(*dir).z;
                if through > 0 {
                    self.through_spin_escaped += spin;
                    self.through_escaped_band[band as usize] += 1.0;
                }
            }
            Termination::Lost => self.lost += 1.0,
        }
    }
}

#[allow(clippy::too_many_arguments)]
#[allow(clippy::too_many_arguments)]
fn run_pass(
    kind: PhotonKind,
    n: usize,
    swing_boost: f64,
    cross: &CollisionField,
    sigma_cross: f64,
    selff: &CollisionField,
    sigma_self: f64,
    n_ricochet: u32,
    gears: &GearConfig,
    rng: &mut u64,
) -> (DensityField, PassTally, f64) {
    let mut field = DensityField::new();
    let mut tally = PassTally::default();
    let mut lz_source_total = 0.0;
    for _ in 0..n {
        match kind {
            PhotonKind::Emitted => {
                let photon = spawn_emitted(rng, swing_boost, gears.class);
                lz_source_total += photon.pos.cross(photon.dir).z * photon.weight;
                let (term, cb, sb, spin, _through) = march(
                    photon, cross, sigma_cross, selff, sigma_self, n_ricochet, gears, rng,
                    &mut field,
                );
                tally.record_emitted(&term, cb, sb, photon.chirality, spin);
            }
            PhotonKind::Ambient => {
                let (photon, aim_lat) = spawn_ambient(rng, gears);
                let (term, cb, sb, spin, through) = march(
                    photon, cross, sigma_cross, selff, sigma_self, n_ricochet, gears, rng,
                    &mut field,
                );
                tally.record_ambient(&term, cb, sb, aim_lat, spin, through);
            }
        }
    }
    (field, tally, lz_source_total)
}

/// Free-streaming (sigma=0) seed passes: their tallies seed the correction
/// and mean-direction fields for iteration 1 of every sweep cell. The
/// launch-shell density spike this produces near r=1.001 is harmless in v2:
/// the correction cap bounds its effect on collisions at 4× analytic.
fn run_iteration0(
    n: usize,
    swing_boost: f64,
    gears: &GearConfig,
    rng: &mut u64,
) -> (DensityField, DensityField) {
    let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient, gears.class);
    let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted, gears.class);
    let (emitted0, _, _) = run_pass(
        PhotonKind::Emitted,
        n,
        swing_boost,
        &pure_ambient,
        0.0,
        &pure_emitted,
        0.0,
        0,
        gears,
        rng,
    );
    let (ambient0, _, _) = run_pass(
        PhotonKind::Ambient,
        n,
        swing_boost,
        &pure_emitted,
        0.0,
        &pure_ambient,
        0.0,
        0,
        gears,
        rng,
    );
    (emitted0, ambient0)
}

struct ReportCell {
    tau_target: f64,
    /// v3 self-collision optical-depth target (0 = v2 behavior).
    tau_self: f64,
    /// v4 pump mode: emitted budget = previous iteration's absorption.
    pump: bool,
    /// Per-iteration (ambient_absorbed, emitted_reabsorbed, next_budget);
    /// empty in forced-source mode.
    budget_track: Vec<(f64, f64, usize)>,
    /// The final iteration's emitted launch count (== n in forced mode).
    n_emit_final: usize,
    n_ricochet: u32,
    sigma_e: f64,
    sigma_a: f64,
    /// Emitted-through-emitted self coupling (v3).
    sigma_ee: f64,
    /// Ambient-through-ambient self coupling (v3, near-zone gated).
    sigma_aa: f64,
    /// Per-iteration L2 relative change of the (regularized) density fields,
    /// as (emitted, ambient) pairs.
    convergence: Vec<(f64, f64)>,
    emitted_tally: PassTally,
    ambient_tally: PassTally,
    lz_source_total: f64,
    final_emitted_field: DensityField,
    final_ambient_field: DensityField,
}

#[allow(clippy::too_many_arguments)]
fn run_self_consistent(
    n: usize,
    swing_boost: f64,
    tau_target: f64,
    tau_self: f64,
    n_ricochet: u32,
    iters: usize,
    pump: bool,
    gears: &GearConfig,
    emitted0: &DensityField,
    ambient0: &DensityField,
    seed: u64,
) -> ReportCell {
    // v2 fix #3: symmetric TAU loading from the ANALYTIC seeds, both knobs
    // derived from the one swept TAU.
    let (sigma_e, sigma_a) = if tau_target > 0.0 {
        (
            tau_target / tau_in_integral(gears.class),
            tau_target / tau_out_integral(),
        )
    } else {
        (0.0, 0.0)
    };
    // v3: the self channel's couplings from the one swept TAU_SELF.
    let (sigma_ee, sigma_aa) = if tau_self > 0.0 {
        (
            tau_self / tau_self_emitted_integral(gears.class),
            tau_self / tau_self_ambient_integral(),
        )
    } else {
        (0.0, 0.0)
    };

    let mut rng = seed;
    let mut tally_emitted = emitted0.clone();
    let mut tally_ambient = ambient0.clone();
    let mut convergence = Vec::with_capacity(iters);
    let mut emitted_tally = PassTally::default();
    let mut ambient_tally = PassTally::default();
    let mut lz_source_total = 0.0;
    let mut final_emitted_field = tally_emitted.clone();
    let mut final_ambient_field = tally_ambient.clone();

    // v4 pump: emitted budget starts at 0 (nothing has been absorbed yet)
    // and becomes the previous iteration's measured absorption. Forced
    // mode keeps the fixed N source (v2/v3 behavior).
    let mut n_emit = if pump { 0 } else { n };
    let mut budget_track: Vec<(f64, f64, usize)> = Vec::new();

    for _ in 0..iters {
        let field_emitted =
            CollisionField::from_tally(PhotonKind::Emitted, gears.class, tally_emitted.clone());
        let field_ambient =
            CollisionField::from_tally(PhotonKind::Ambient, gears.class, tally_ambient.clone());

        // Emitted march through the AMBIENT field with sigma_a (cross) and
        // their OWN field with sigma_ee (self); ambient march through the
        // EMITTED field with sigma_e (cross) and their own with sigma_aa.
        let (measured_emitted, e_tally, lz_src) = run_pass(
            PhotonKind::Emitted,
            n_emit,
            swing_boost,
            &field_ambient,
            sigma_a,
            &field_emitted,
            sigma_ee,
            n_ricochet,
            gears,
            &mut rng,
        );
        let (measured_ambient, a_tally, _) = run_pass(
            PhotonKind::Ambient,
            n,
            swing_boost,
            &field_emitted,
            sigma_e,
            &field_ambient,
            sigma_aa,
            n_ricochet,
            gears,
            &mut rng,
        );

        if pump {
            let next = (a_tally.arrived + e_tally.arrived).round() as usize;
            budget_track.push((a_tally.arrived, e_tally.arrived, next));
            n_emit = next;
        }

        let new_emitted = tally_emitted.relax(&measured_emitted, 0.5);
        let new_ambient = tally_ambient.relax(&measured_ambient, 0.5);

        let old_e_snap =
            effective_density_snapshot(PhotonKind::Emitted, gears.class, &tally_emitted);
        let new_e_snap = effective_density_snapshot(PhotonKind::Emitted, gears.class, &new_emitted);
        let old_a_snap =
            effective_density_snapshot(PhotonKind::Ambient, gears.class, &tally_ambient);
        let new_a_snap = effective_density_snapshot(PhotonKind::Ambient, gears.class, &new_ambient);
        convergence.push((
            l2_rel_change(&old_e_snap, &new_e_snap),
            l2_rel_change(&old_a_snap, &new_a_snap),
        ));

        tally_emitted = new_emitted;
        tally_ambient = new_ambient;
        emitted_tally = e_tally;
        ambient_tally = a_tally;
        lz_source_total = lz_src;
        final_emitted_field = measured_emitted;
        final_ambient_field = measured_ambient;
    }

    ReportCell {
        tau_target,
        tau_self,
        pump,
        budget_track,
        n_emit_final: emitted_tally.total as usize,
        n_ricochet,
        sigma_e,
        sigma_a,
        sigma_ee,
        sigma_aa,
        convergence,
        emitted_tally,
        ambient_tally,
        lz_source_total,
        final_emitted_field,
        final_ambient_field,
    }
}

// ---------------------------------------------------------------------
// CM-4: the axial stack (multi-body)
// ---------------------------------------------------------------------
//
// The single-particle cell established that the through-charge limb WORKS
// (species- and spin-preserving pole-to-pole crossing) but is only ~0.67% of
// arriving charge, saturating near 4.7% under strong axial channelling — far
// too small to deliver magnetic neutrality, because the equatorial disc limb
// still carries 95-99% of output. venus2.pdf says exactly that for an isolated
// body ("normally a minor complication") and locates the dominance at the
// NUCLEAR level, where two things change that a one-body cell cannot express:
// the pole-to-pole distance is short, and the lane is fed by a NEIGHBOUR's
// exhaust rather than by isotropic ambient. That is what this section tests.
//
// The seating is NOT free. graphene.pdf: "the proton is plugged in with its
// equator pointing down. But the neutron is plugged in with its pole pointing
// down. This is because protons channel charge pole to equator, while neutrons
// channel pole to pole." So a nuclear stack aims a proton's EQUATORIAL DISC —
// the limb that carries ~95-99% of its output — straight into the neutron's
// polar lane. Atom mode already asserts this convention for plugs
// (atom_scenarios: plug protons edge-on `rest_axis·Y ≈ 0` "disc feeds the
// hole", plug neutrons `|rest_axis·Y| ≈ 1` "pole on the stack axis"), so the
// two modes agree on geometry by construction rather than by coincidence.

/// Stacked-nucleon centre separation in body radii — `atom_core::NUCLEON_PITCH`
/// ("a fused neighbor parks at the MOUTH of the source's intake funnel"). Kept
/// as a local copy so the recycling cell stays independent of Atom mode; the
/// stack report sweeps around it.
pub const STACK_PITCH: f64 = 2.6;

/// Seating of a nucleon in an axial stack. The two cases are physically
/// distinct limbs, not a cosmetic rotation: which of the body's own limbs
/// points along the stack axis decides what the neighbour receives.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum StackOrient {
    /// Spin axis ALONG the stack axis: the body channels pole-to-pole through
    /// the stack, so its through-lane opens onto its neighbours. The neutron's
    /// nuclear seating.
    PoleOn,
    /// Spin axis PERPENDICULAR to the stack axis: the equatorial disc fires
    /// along the stack and the through-lane points sideways, away from the
    /// neighbours. The proton's nuclear seating.
    EdgeOn,
}

impl StackOrient {
    /// The body's spin axis in stack coordinates (stack axis = +z).
    fn spin_axis(self) -> DVec3 {
        match self {
            StackOrient::PoleOn => DVec3::Z,
            StackOrient::EdgeOn => DVec3::X,
        }
    }

    /// The nuclear seating for a class, per graphene.pdf. The electron is not a
    /// nucleon; it is given the pole-on convention only so the enum is total.
    fn nuclear_for(class: ParticleClass) -> Self {
        match class {
            ParticleClass::Proton => StackOrient::EdgeOn,
            ParticleClass::Neutron | ParticleClass::Electron => StackOrient::PoleOn,
        }
    }

    fn label(self) -> &'static str {
        match self {
            StackOrient::PoleOn => "pole-on",
            StackOrient::EdgeOn => "edge-on",
        }
    }
}

/// One body in an axial stack: a centre on the stack axis plus its own spin
/// orientation. Every per-body quantity (density lookup, through-lane gate,
/// gear redirect) is evaluated in this body's frame, never the stack's.
#[derive(Clone, Copy)]
struct StackBody {
    center: DVec3,
    class: ParticleClass,
    orient: StackOrient,
}

impl StackBody {
    fn spin_axis(&self) -> DVec3 {
        self.orient.spin_axis()
    }

    /// Rotation carrying BODY-LOCAL coordinates (spin axis = +z — the frame
    /// every sampler, density field and gear rule in this module is written in)
    /// into stack coordinates. For a pole-on body this is the identity.
    fn to_stack(&self) -> DQuat {
        DQuat::from_rotation_arc(DVec3::Z, self.spin_axis())
    }

    /// (radius, colatitude) of a stack-frame point in this body's own frame —
    /// the pair every `CollisionField` lookup wants.
    fn local_polar(&self, pos: DVec3) -> (f64, f64) {
        let d = pos - self.center;
        let r = d.length();
        let ct = (d.dot(self.spin_axis()) / r.max(1e-12)).clamp(-1.0, 1.0);
        (r, ct.acos())
    }

    /// Latitude of a surface point in this body's own frame (+90 = its north
    /// pole), which is what the exit-topology bands mean.
    fn local_latitude(&self, pos: DVec3) -> f64 {
        let (_, theta) = self.local_polar(pos);
        90.0 - theta.to_degrees()
    }

    /// The two quantities the through-lane gate needs, about this body's OWN
    /// spin axis: cylindrical radius of `pos`, and the axial component of
    /// `dir`. An edge-on body's lane runs along x, so axial feed along the
    /// stack misses it — that asymmetry is the point of the measurement.
    fn lane_coords(&self, pos: DVec3, dir: DVec3) -> (f64, f64) {
        let d = pos - self.center;
        let a = self.spin_axis();
        let along = d.dot(a);
        let rho = (d - a * along).length();
        (rho, dir.dot(a))
    }
}

/// What one stacked march did, from the point of view of the body being
/// measured (`probe`).
struct StackOutcome {
    term: Termination,
    spin: f64,
    /// Reached the probe's surface at all (whether or not it then crossed).
    probe_arrived: bool,
    /// Latitude of that first contact, in the PROBE's own frame.
    probe_arrival_lat: f64,
    /// Through-lane crossings of the probe.
    probe_crossed: u32,
    /// Absorbed by some OTHER body — the neighbour shadowing/reabsorption
    /// channel, which a one-body cell has no way to charge for.
    other_absorbed: bool,
    /// WHICH body absorbed it, if any. Together with the direction stored in
    /// `Termination::Absorbed` this is the linear-momentum bookkeeping: the
    /// stream's push on a neighbour is the sum of absorbed photon directions,
    /// per cc.pdf's open field (absorber force only — emitter recoil is not a
    /// channel in this model, see [[project-momentum-open-field]]).
    absorbed_by: Option<usize>,
}

/// Ray-marches one photon through a STACK of bodies. Same transport physics as
/// `march`, with three differences forced by there being more than one centre:
///
/// 1. Absorption and the through-lane gate are tested against every body, each
///    in its own frame (`StackBody::lane_coords`), so seating matters.
/// 2. The chirality gears are geometric — `scatter_off` funnels cancels to
///    "the nearer pole" and flattens augments toward "the equatorial plane",
///    both keyed to +z. So the collision is evaluated with the photon rotated
///    into the FIELD BODY's local frame and the new direction rotated back.
///    Skipping that would give an edge-on proton a disc lying in the stack's
///    xy plane instead of its own.
/// 3. Escape is measured from the stack CENTROID, so the tallied shell encloses
///    the whole stack rather than one body.
///
/// Phase 1 deliberately gives only `field_body` a converged gas: the probe is a
/// passive absorber. The probe's own gas would scatter arriving charge, some of
/// it OUT of the lane, so the on-lane rate measured here is an UPPER bound on
/// the coupled value. Transport-level shadowing by the probe's solid body IS
/// included, since real photons are marched through the real geometry. Nor does
/// this deposit track-length into a field — the stack is not iterated to
/// self-consistency here, it reads a feed rate off converged single-body gases.
#[allow(clippy::too_many_arguments)]
fn march_stack(
    mut photon: Photon,
    bodies: &[StackBody],
    field_body: usize,
    probe: usize,
    cross: &CollisionField,
    sigma_cross: f64,
    selff: &CollisionField,
    sigma_self: f64,
    n_ricochet: u32,
    gears: &GearConfig,
    rng: &mut u64,
) -> StackOutcome {
    let mut centroid = DVec3::ZERO;
    for b in bodies {
        centroid += b.center;
    }
    centroid /= bodies.len() as f64;

    let source = &bodies[field_body];
    let q_to_stack = source.to_stack();
    let q_to_local = q_to_stack.inverse();

    let mut path_length = 0.0;
    let mut probe_arrived = false;
    let mut probe_arrival_lat = f64::NAN;
    let mut probe_crossed = 0u32;
    let mut other_absorbed = false;
    let mut absorbed_by: Option<usize> = None;

    let term = loop {
        // Nearest surface decides the step, so resolution stays fine near
        // either body rather than only near the origin.
        let mut r_min = f64::INFINITY;
        let mut hit: Option<usize> = None;
        for (i, b) in bodies.iter().enumerate() {
            let r_i = (photon.pos - b.center).length();
            if r_i < r_min {
                r_min = r_i;
            }
            if r_i < R_IN && hit.is_none() {
                hit = Some(i);
            }
        }

        if let Some(i) = hit {
            let body = &bodies[i];
            let r_i = (photon.pos - body.center).length();
            if i == probe && !probe_arrived {
                probe_arrived = true;
                probe_arrival_lat = body.local_latitude(photon.pos);
            }
            // Through-charge limb, same two geometric gates as the single-body
            // version (salt.pdf: near the centre of the hole passes, edge or
            // angled entry is thrown to the equatorial whirlpool) — but in
            // THIS body's frame. Species and spin are retained; that retention
            // is what makes a mixed-species exhaust possible at all.
            if gears.through_lane_rho > 0.0 {
                let (rho, axial) = body.lane_coords(photon.pos, photon.dir);
                let d = photon.pos - body.center;
                let incoming = photon.dir.dot(d / r_i.max(1e-12));
                if rho < gears.through_lane_rho * r_i
                    && axial.abs() > THROUGH_COS_MAX_ANGLE
                    && incoming < 0.0
                {
                    let t = -2.0 * d.dot(photon.dir);
                    if t > 0.0 {
                        let exit = d + photon.dir * t;
                        photon.pos =
                            body.center + exit * (R_IN / exit.length().max(1e-12)) * 1.001;
                        if i == probe {
                            probe_crossed += 1;
                        }
                        path_length += t;
                        continue;
                    }
                }
            }
            if i != probe {
                other_absorbed = true;
            }
            absorbed_by = Some(i);
            break Termination::Absorbed {
                pos: photon.pos,
                dir: photon.dir,
            };
        }

        if (photon.pos - centroid).length() > R_OUT {
            break Termination::Escaped {
                pos: photon.pos,
                dir: photon.dir,
            };
        }
        if path_length > SAFETY_PATH_CAP {
            break Termination::Lost;
        }

        let ds = (DS_FRACTION * r_min).max(MIN_DS);
        let (r_f, theta_f) = source.local_polar(photon.pos);

        if (sigma_cross > 0.0 || sigma_self > 0.0) && photon.bounces < n_ricochet {
            let lam_cross = if sigma_cross > 0.0 {
                sigma_cross * cross.effective_density_at(r_f, theta_f)
            } else {
                0.0
            };
            let lam_self = if sigma_self > 0.0
                && !(selff.kind == PhotonKind::Ambient && r_f > B_MAX)
            {
                sigma_self * selff.effective_density_at(r_f, theta_f)
            } else {
                0.0
            };
            let lam = lam_cross + lam_self;
            if lam > 0.0 {
                let p_collide = 1.0 - (-lam * ds).exp();
                if xorshift64(rng) < p_collide {
                    let use_self = if lam_cross <= 0.0 {
                        true
                    } else if lam_self <= 0.0 {
                        false
                    } else {
                        xorshift64(rng) < lam_self / lam
                    };
                    let field = if use_self { selff } else { cross };
                    // Into the field body's frame for the gear geometry, back
                    // out again afterwards (see doc comment, point 2).
                    let pos_local = q_to_local * (photon.pos - source.center);
                    let dir_local = q_to_local * photon.dir;
                    let (new_dir_local, spin_delta) = scatter_off(
                        field,
                        dir_local,
                        pos_local,
                        r_f,
                        theta_f,
                        photon.chirality,
                        gears,
                        rng,
                    );
                    photon.dir = q_to_stack * new_dir_local;
                    photon.spin = (photon.spin + spin_delta).clamp(-2.0, 2.0);
                    photon.bounces += 1;
                }
            }
        }

        photon.pos += photon.dir * ds;
        path_length += ds;
    };

    StackOutcome {
        term,
        spin: photon.spin,
        probe_arrived,
        probe_arrival_lat,
        probe_crossed,
        other_absorbed,
        absorbed_by,
    }
}

/// Spawns emitted charge on `body` in STACK coordinates: the body-local
/// sampler is used unchanged (so the exit profile and RNG draws are the same
/// ones the single-body cell validated) and the result is rotated onto the
/// body's actual spin axis and translated to its centre.
fn spawn_emitted_on(rng: &mut u64, swing_boost: f64, body: &StackBody) -> Photon {
    let mut p = spawn_emitted(rng, swing_boost, body.class);
    let q = body.to_stack();
    p.pos = body.center + q * p.pos;
    p.dir = q * p.dir;
    p
}

/// Ambient launch for a stack: the single-body importance ensemble, recentred
/// on the stack centroid. `B_MAX = 6` still covers both bodies of a
/// `STACK_PITCH` pair (each sits 1.3 from the centroid). The species-split
/// `pole_intake_bias` inside `spawn_ambient` steers along +z, which in stack
/// coordinates IS the stack axis — the correct axis for a nuclear channel.
fn spawn_ambient_on_stack(rng: &mut u64, gears: &GearConfig, centroid: DVec3) -> Photon {
    let (mut p, _aim) = spawn_ambient(rng, gears);
    p.pos += centroid;
    p
}

// ---------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    /// (a) TAU=0 (sigma=0, no collisions) with the v2 importance-sampled
    /// launch: every b ≤ 1 ambient photon must arrive at its straight-line
    /// aim latitude; every b > 1 photon must escape; zero collisions
    /// anywhere; zero emitted re-entries.
    ///
    /// NOTE on tolerance: sampled march positions lie exactly on the launch
    /// line, whose minimum radius is b — so false hits (b > 1 arriving) are
    /// geometrically impossible. The only discretization artifact left is a
    /// near-tangent b → 1⁻ graze whose sub-surface chord (2·sqrt(1−b²)) is
    /// thinner than the ~0.05 step near r=1, so the march can hop over the
    /// crossing and report an escape. With area-weighted b, the fraction of
    /// launches in that sliver is O(0.1%), so a small mismatch tolerance
    /// covers it without hiding a real transport bug.
    #[test]
    fn tau0_transport_matches_straight_line_geometry() {
        let mut rng = 0xA5A5_1234_9876_5432u64;
        let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);
        let n = 5000;

        let mut checked = 0;
        let mut hits = 0;
        let mut mismatches = 0;
        for _ in 0..n {
            let (photon, aim_lat) = spawn_ambient(&mut rng, &GearConfig::off());
            let mut field = DensityField::new();
            let (term, cb, sb, _spin, _through) = march(
                photon,
                &pure_emitted,
                0.0,
                &pure_ambient,
                0.0,
                0,
                &GearConfig::off(),
                &mut rng,
                &mut field,
            );
            assert_eq!(cb + sb, 0, "collision recorded with sigma=0");
            match (aim_lat, &term) {
                (Some(aim), Termination::Absorbed { pos, .. }) => {
                    let arrival = arrival_latitude(*pos);
                    if (arrival - aim).abs() >= 5.0 {
                        mismatches += 1;
                    }
                    hits += 1;
                    checked += 1;
                }
                (Some(_aim), Termination::Escaped { .. }) => {
                    // Near-tangent graze hopped by the finite step — see note.
                    mismatches += 1;
                    checked += 1;
                }
                (None, Termination::Escaped { .. }) => {
                    checked += 1;
                }
                (None, Termination::Absorbed { .. }) => {
                    panic!("b>1 photon arrived under sigma=0 (line min radius is b — impossible)");
                }
                (_, Termination::Lost) => {
                    panic!("photon lost under TAU=0 (no collisions, should never hit the path cap)");
                }
            }
        }
        assert_eq!(checked, n);
        // Area-weighted b: P(b <= 1) = (1/B_MAX)^2 = 1/36 -> ~139 of 5000.
        assert!(
            hits > 50,
            "expected ~n/36 straight-aim hits from the importance-sampled launch, got {hits}"
        );
        let mismatch_frac = mismatches as f64 / n as f64;
        assert!(
            mismatch_frac < 0.01,
            "too many aim/arrival mismatches for a sigma=0 straight-line transport: {mismatches}/{n} ({:.2}%)",
            mismatch_frac * 100.0
        );

        // Emitted photons at TAU=0 must never re-enter r<1 (straight lines
        // from the surface, always genuinely outward — see module doc
        // deviation #1 for why this is guaranteed by construction).
        let mut reentries = 0;
        for _ in 0..n {
            let photon = spawn_emitted(&mut rng, DEFAULT_SWING_BOOST, ParticleClass::Proton);
            let mut field = DensityField::new();
            let (term, cb, sb, _spin, _through) = march(
                photon,
                &pure_ambient,
                0.0,
                &pure_emitted,
                0.0,
                0,
                &GearConfig::off(),
                &mut rng,
                &mut field,
            );
            assert_eq!(cb + sb, 0, "collision recorded with sigma=0");
            if matches!(term, Termination::Absorbed { .. }) {
                reentries += 1;
            }
            assert!(
                !matches!(term, Termination::Lost),
                "emitted photon lost under TAU=0"
            );
        }
        assert_eq!(reentries, 0, "TAU=0 emitted photons should never bend back into r<1");
    }

    /// (b) Free-streaming emitted density on the equator falls off ~1/r²
    /// (ratio check between two well-separated radial bins, 20% tolerance).
    ///
    /// NOTE: emission position latitude (lat0) and exit direction latitude
    /// (lat1) are independent draws from the same bimodal mixture (module
    /// doc deviation #1), so an individual photon's trajectory's theta bin
    /// drifts from ~lat0 near the surface toward ~lat1 asymptotically. A
    /// SINGLE narrow theta bin straddling the equator is (a) populated only
    /// by the Gaussian tail of the 7°-peak (mean 7°, sigma 5°), and (b) at
    /// the exact θ=π/2 boundary `theta_bin_index` resolves to one specific
    /// 10°-wide bin, asymmetric about true 0° latitude — both effects
    /// measurably distort a single-bin ratio (empirically: ~20-21x instead
    /// of the naive ~29-30x expected for narrow-bin choices, stable across
    /// N=150k..500k, i.e. a real geometric effect, not noise). Fix:
    /// aggregate over the *report's own* equatorial band definition
    /// (|lat|<30°, theta bins 6..=11) so the comparison is symmetric about
    /// the equator and has ~6x the single-bin statistics, then compare two
    /// well-separated, non-adjacent-to-source radial bins.
    #[test]
    fn free_streaming_emitted_density_falls_like_inverse_r_squared_on_equator() {
        let mut rng = 0x1357_2468_ABCD_EF01u64;
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);
        let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);
        let (field, _tally, _lz) = run_pass(
            PhotonKind::Emitted,
            300_000,
            DEFAULT_SWING_BOOST,
            &pure_ambient,
            0.0,
            &pure_emitted,
            0.0,
            0,
            &GearConfig::off(),
            &mut rng,
        );

        let equatorial_band_js: Vec<usize> = (6..=11).collect(); // |lat| < 30 deg
        let aggregate_density = |i: usize| -> f64 {
            let mut tw = 0.0;
            let mut vol = 0.0;
            for &j in &equatorial_band_js {
                let idx = i * N_THETA + j;
                tw += field.track_weight[idx];
                vol += cell_volumes()[idx];
            }
            tw / vol
        };

        let i_a = 10usize;
        let i_b = 20usize;
        let d_a = aggregate_density(i_a);
        let d_b = aggregate_density(i_b);
        assert!(
            d_a > 0.0 && d_b > 0.0,
            "insufficient equatorial-band statistics: d_a={d_a}, d_b={d_b}"
        );

        let r_a = (r_edge(i_a) * r_edge(i_a + 1)).sqrt();
        let r_b = (r_edge(i_b) * r_edge(i_b + 1)).sqrt();
        let observed_ratio = d_a / d_b;
        let expected_ratio = (r_b / r_a).powi(2);
        let rel_err = (observed_ratio - expected_ratio).abs() / expected_ratio;
        assert!(
            rel_err < 0.20,
            "observed d_a/d_b={observed_ratio:.4}, expected~{expected_ratio:.4} (1/r²), rel_err={rel_err:.3}"
        );
    }

    /// (c) swing_boost=0, TAU=0: net delivered L_z per photon ~ 0 — no
    /// systematic angular-momentum source means no systematic delivery.
    ///
    /// NOTE: per-photon |L_z| = |(pos×dir).z| ≤ the line's impact parameter,
    /// which the v2 launch bounds at B_MAX=6 — so the sampling-noise floor
    /// for the mean over N=20k photons is ~6/sqrt(20000) ≈ 0.04. The 0.2
    /// threshold gives comfortable headroom over that floor while staying
    /// far below any physically-meaningful systematic delivery.
    #[test]
    fn zero_swing_zero_tau_delivers_no_net_angular_momentum() {
        let mut rng = 0xFEED_BEEF_1234_5678u64;
        let n = 20_000usize;
        let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);
        let (_field_e, emitted_tally, _lz_src) = run_pass(
            PhotonKind::Emitted,
            n,
            0.0,
            &pure_ambient,
            0.0,
            &pure_emitted,
            0.0,
            0,
            &GearConfig::off(),
            &mut rng,
        );
        let (_field_a, ambient_tally, _) = run_pass(
            PhotonKind::Ambient,
            n,
            0.0,
            &pure_emitted,
            0.0,
            &pure_ambient,
            0.0,
            0,
            &GearConfig::off(),
            &mut rng,
        );

        let delivered = ambient_tally.lz_arrived + emitted_tally.lz_arrived;
        let carried_off = ambient_tally.lz_escaped + emitted_tally.lz_escaped;
        let net = delivered - carried_off;
        let net_per_photon = net / n as f64;
        assert!(
            net_per_photon.abs() < 0.2,
            "expected ~0 net L_z/photon at swing_boost=0, got {net_per_photon:.5}"
        );
    }

    /// (d) v3 self-collision channel: wiring, the shared ricochet cap, and
    /// the ambient-ambient far-zone gate.
    #[test]
    fn self_collision_channel_wires_caps_and_gates() {
        let mut rng = 0xDEAD_BEEF_CAFE_0001u64;
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);
        let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);

        // Emitted photons, cross channel OFF, self coupling enormous: every
        // collision must land on the self channel, capped by n_ricochet.
        // (Launch latitudes sit on the profile peaks where the analytic
        // self-density is ~0.25, so sigma_ee=1000 collides within a step.)
        let cap = 5u32;
        let mut any_self = false;
        for _ in 0..200 {
            let photon = spawn_emitted(&mut rng, DEFAULT_SWING_BOOST, ParticleClass::Proton);
            let mut field = DensityField::new();
            let (_term, cb, sb, _spin, _through) = march(
                photon,
                &pure_ambient,
                0.0,
                &pure_emitted,
                1000.0,
                cap,
                &GearConfig::off(),
                &mut rng,
                &mut field,
            );
            assert_eq!(cb, 0, "cross collision recorded with sigma_cross=0");
            assert!(sb <= cap, "self bounces {sb} exceeded the shared ricochet cap {cap}");
            if sb > 0 {
                any_self = true;
            }
        }
        assert!(
            any_self,
            "huge sigma_ee produced zero self collisions — self channel not wired"
        );

        // Ambient far-zone gate: a line entirely outside r = B_MAX must
        // never self-collide, at ANY coupling.
        let far = Photon {
            pos: DVec3::new(B_MAX + 2.0, 0.0, 0.0),
            dir: DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 0.0,
            spin: 0.0,
        };
        let mut field = DensityField::new();
        let (term, cb, sb, _spin, _through) = march(
            far,
            &pure_emitted,
            0.0,
            &pure_ambient,
            1000.0,
            100,
            &GearConfig::off(),
            &mut rng,
            &mut field,
        );
        assert_eq!(cb, 0);
        assert_eq!(
            sb, 0,
            "ambient self-collision fired outside the near zone — far-zone gate broken"
        );
        assert!(matches!(term, Termination::Escaped { .. }));

        // Inside the near zone the same coupling must fire immediately.
        let near = Photon {
            pos: DVec3::new(2.0, 0.0, 0.0),
            dir: DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 0.0,
            spin: 0.0,
        };
        let mut field = DensityField::new();
        let (_term, _cb, sb, _spin, _through) = march(
            near,
            &pure_emitted,
            0.0,
            &pure_ambient,
            1000.0,
            100,
            &GearConfig::off(),
            &mut rng,
            &mut field,
        );
        assert!(
            sb > 0,
            "ambient self-collision failed to fire inside the near zone"
        );
    }

    /// (e) v4 pump mode at TAU=0: the endogenous emission budget must settle
    /// to the pure geometric absorption rate (P(b<=1) = 1/36 of the
    /// conditional ambient ensemble), and TAU=0 recycled photons never
    /// re-enter (straight outward lines), so the budget is ambient-fed only.
    #[test]
    fn pump_budget_tracks_geometric_absorption_at_tau0() {
        let mut rng = 0xBEE5_0000_1234_ABCDu64;
        let n = 20_000usize;
        let gears = GearConfig::off();
        let (emitted0, ambient0) = run_iteration0(n, DEFAULT_SWING_BOOST, &gears, &mut rng);
        let cell = run_self_consistent(
            n,
            DEFAULT_SWING_BOOST,
            0.0,
            0.0,
            0,
            4,
            true,
            &gears,
            &emitted0,
            &ambient0,
            0x5eed_5eed_5eed_5eedu64,
        );
        assert_eq!(cell.budget_track.len(), 4);

        let (a1, e1, next1) = cell.budget_track[0];
        assert!(a1 > 0.0, "ambient absorption must be nonzero");
        assert_eq!(e1, 0.0, "iteration 1 launches zero emitted photons");
        assert_eq!(next1, a1 as usize);

        let e = &cell.emitted_tally;
        assert_eq!(e.arrived, 0.0, "TAU=0 recycled photons must never re-enter");
        assert!(cell.n_emit_final > 0, "pump budget failed to spin up");

        let a = &cell.ambient_tally;
        let frac = a.arrived / a.total;
        assert!(
            (frac - 1.0 / 36.0).abs() < 0.01,
            "geometric absorption fraction {frac:.4} far from 1/36"
        );
    }

    /// The CM-3 v3 report: sweeps TAU (cross optical depth, v2 fix #3) x
    /// TAU_SELF (same-population optical depth, v3) x N_RICOCHET,
    /// self-consistently settles both populations' regularized density
    /// fields for ITERS passes, and prints the screening, redirection,
    /// recycling, circulation, PRESSURE, and spin-sustenance evidence for
    /// the disc-to-pole loop claim — plus per-channel knob verification
    /// (measured vs target collision counts for BOTH channels, printed
    /// loudly when off). All ambient statistics are CONDITIONAL on the
    /// importance-sampled near-zone ensemble (lines with impact parameter
    /// b ≤ B_MAX=6; see module doc, v2 fix #1).
    #[test]
    #[ignore]
    fn report_recycling_loop() {
        use std::time::Instant;
        let started = Instant::now();

        let n = DEFAULT_N;
        let swing_boost = DEFAULT_SWING_BOOST;
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;

        println!("=== CM-3 Recycling Loop Report (v4: v3 pressure gas + pump surface (endogenous emission budget) + charge-pause observables) ===");
        println!(
            "N per population = {n}, ITERS = {ITERS}, grid = {N_R}x{N_THETA} (r,theta), R_OUT = {R_OUT}, swing_boost = {swing_boost}, B_MAX = {B_MAX}"
        );
        println!(
            "AMBIENT ENSEMBLE IS CONDITIONAL: all ambient photons are straight-line launches with impact parameter b <= {B_MAX} (area-weighted = correct relative flux within the near zone). Absolute capture fractions are conditional on this ensemble."
        );
        println!(
            "TAU calibration (analytic, cross): TAU_in integral (b=0 equatorial, r={TAU_IN_R_MIN}->{TAU_IN_R_MAX}, profile(0)={:.6}) = {:.6}; TAU_out integral (equatorial emitted, r=1->{R_OUT}, n_a=1) = {:.1}",
            exit_lat_profile(ParticleClass::Proton, 0.0),
            tau_in_integral(ParticleClass::Proton),
            tau_out_integral()
        );
        println!(
            "TAU_SELF calibration (analytic, self): emitted-emitted integral (equatorial, r={TAU_IN_R_MIN}->{R_OUT}) = {:.6}; ambient-ambient integral (near-zone path, B_MAX-1) = {:.1}. Ambient self channel is GATED to r <= B_MAX (far-zone self-scattering of a uniform isotropic gas is ensemble-null).",
            tau_self_emitted_integral(ParticleClass::Proton),
            tau_self_ambient_integral()
        );

        // This report keeps the original chirality-blind physics (proton
        // disc profile, no gears) as the v4 regression anchor. The
        // proton-vs-neutron chirality comparison is a separate report.
        let gears = GearConfig::off();
        let (emitted0, ambient0) = run_iteration0(n, swing_boost, &gears, &mut seed);

        // (pump, TAU, TAU_SELF, N_RICOCHET, iters) sweep — v4. The forced
        // flagship is the v3 regression anchor; every pump cell's emission
        // budget is endogenous (= previous iteration's absorption). The
        // pump null (TAU=0, TAU_SELF=0) measures the pure geometric engine;
        // the self-only pump cell isolates pressure with cross off.
        let specs: Vec<(bool, f64, f64, u32, usize)> = vec![
            (false, 3.0, 3.0, 30, ITERS),
            (true, 0.0, 0.0, 0, ITERS_PUMP),
            (true, 3.0, 0.0, 30, ITERS_PUMP),
            (true, 3.0, 1.0, 30, ITERS_PUMP),
            (true, 3.0, 3.0, 30, ITERS_PUMP),
            (true, 3.0, 10.0, 30, ITERS_PUMP),
            (true, 0.0, 3.0, 30, ITERS_PUMP),
        ];

        let mut cells: Vec<ReportCell> = Vec::new();
        for &(pump, tau, tau_self, nr, iters) in &specs {
            seed = seed.wrapping_add(1);
            let cell = run_self_consistent(
                n, swing_boost, tau, tau_self, nr, iters, pump, &gears, &emitted0, &ambient0, seed,
            );
            println!(
                "cell {} TAU={tau:.1} TAU_SELF={tau_self:.1} Nric={nr} iters={iters} -> sigma_e={:.4} sigma_a={:.6} sigma_ee={:.4} sigma_aa={:.6} n_emit_final={}",
                if pump { "PUMP" } else { "FORCED" },
                cell.sigma_e,
                cell.sigma_a,
                cell.sigma_ee,
                cell.sigma_aa,
                cell.n_emit_final
            );
            cells.push(cell);
        }

        println!();
        println!("KNOB VERIFICATION + main sweep table (per cell; cAmb/cEm = MEASURED mean CROSS collisions per ambient/emitted photon vs TAU; sAmb/sEm = MEASURED mean SELF collisions vs TAU_SELF; bArr/bReent = mean ricochets among ambient arrivals / emitted re-entries; nEmit = final-iteration emitted budget — pump cells' emitted stats ride on this smaller sample):");
        println!(
            "{:>4} {:>5} {:>5} {:>4} {:>7} | {:>7} {:>7} {:>7} {:>7} | {:>8} {:>8} {:>8} {:>8} | {:>9} {:>8} {:>8} {:>8} {:>8} | {:>8} {:>8} | {:>6} {:>7} | {:>11} {:>9} {:>5}",
            "mode", "TAU", "TAUs", "Nric", "nEmit",
            "cAmb", "cEm", "sAmb", "sEm",
            "amb_arr", "polar", "mid", "equat",
            "redirEq>P", "missCapt", "missPol", "missMid", "missEq",
            "em_reent", "em_esc",
            "bArr", "bReent",
            "netLz/ph", "fbRatio", "sign"
        );
        let mut knob_warnings: Vec<String> = Vec::new();
        for cell in &cells {
            let a = &cell.ambient_tally;
            let e = &cell.emitted_tally;

            // .max(1.0): a pump cell whose budget collapsed to zero would
            // otherwise print NaN columns; 0/1 = 0 reads correctly as
            // "no emitted photons, no collisions".
            let c_amb = a.collisions_cross_sum / a.total;
            let c_em = e.collisions_cross_sum / e.total.max(1.0);
            let s_amb = a.collisions_self_sum / a.total;
            let s_em = e.collisions_self_sum / e.total.max(1.0);

            let amb_arr_frac = a.arrived / a.total;
            let polar_frac = a.arrived_band[Band::Polar as usize] / a.total;
            let mid_frac = a.arrived_band[Band::Mid as usize] / a.total;
            let eq_frac = a.arrived_band[Band::Equatorial as usize] / a.total;

            let redirect = if a.aim_eq_arrived > 0.0 {
                a.aim_eq_arrived_polar / a.aim_eq_arrived
            } else {
                f64::NAN
            };
            let miss_capture = if a.aim_miss_count > 0.0 {
                a.aim_miss_arrived / a.aim_miss_count
            } else {
                f64::NAN
            };
            let miss_pol = if a.aim_miss_arrived > 0.0 {
                a.aim_miss_arrived_band[Band::Polar as usize] / a.aim_miss_arrived
            } else {
                f64::NAN
            };
            let miss_mid = if a.aim_miss_arrived > 0.0 {
                a.aim_miss_arrived_band[Band::Mid as usize] / a.aim_miss_arrived
            } else {
                f64::NAN
            };
            let miss_eq = if a.aim_miss_arrived > 0.0 {
                a.aim_miss_arrived_band[Band::Equatorial as usize] / a.aim_miss_arrived
            } else {
                f64::NAN
            };

            let em_reentry = e.arrived / e.total.max(1.0);
            let em_escape = e.escaped / e.total.max(1.0);

            let b_arr = if a.arrived > 0.0 {
                a.arrived_bounce_sum / a.arrived
            } else {
                f64::NAN
            };
            let b_reent = if e.arrived > 0.0 {
                e.arrived_bounce_sum / e.arrived
            } else {
                f64::NAN
            };

            let delivered = a.lz_arrived + e.lz_arrived;
            let carried_off = a.lz_escaped + e.lz_escaped;
            let net = delivered - carried_off;
            let net_per_photon = net / n as f64;
            let fb_ratio = if cell.lz_source_total.abs() > 1e-9 {
                net / cell.lz_source_total
            } else {
                f64::NAN
            };
            let sign_label = if net_per_photon.signum() == swing_boost.signum() {
                "match"
            } else {
                "FLIP"
            };
            let nric_label = if cell.tau_target == 0.0 && cell.tau_self == 0.0 {
                "n/a".to_string()
            } else {
                cell.n_ricochet.to_string()
            };

            println!(
                "{:>4} {:>5.1} {:>5.1} {:>4} {:>7} | {:>7.3} {:>7.3} {:>7.3} {:>7.3} | {:>8.4} {:>8.4} {:>8.4} {:>8.4} | {:>9.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} | {:>8.4} {:>8.4} | {:>6.2} {:>7.2} | {:>11.6} {:>9.4} {:>5}",
                if cell.pump { "pump" } else { "frc" },
                cell.tau_target, cell.tau_self, nric_label, cell.n_emit_final,
                c_amb, c_em, s_amb, s_em,
                amb_arr_frac, polar_frac, mid_frac, eq_frac,
                redirect, miss_capture, miss_pol, miss_mid, miss_eq,
                em_reentry, em_escape,
                b_arr, b_reent,
                net_per_photon, fb_ratio, sign_label
            );

            // v2 fix #4, extended per channel: each ACTIVE channel's knob
            // must be PROVEN to work. The SHARED ricochet cap truncates both
            // channels jointly, so each channel's effective target is
            // min(its TAU, N_RICOCHET) and the band stays generous — shout
            // only when a knob is broken by >5x. (When TAU + TAU_SELF
            // exceeds the cap, both channels read low together; that shows
            // up here as ratios sliding toward 0.2, by design.)
            let checks = [
                ("cross/ambient", cell.tau_target, c_amb),
                ("cross/emitted", cell.tau_target, c_em),
                ("self/ambient", cell.tau_self, s_amb),
                ("self/emitted", cell.tau_self, s_em),
            ];
            for (label, tau_knob, measured) in checks {
                if tau_knob > 0.0 {
                    let target = tau_knob.min(cell.n_ricochet as f64);
                    let ratio = measured / target;
                    if !(0.2..=5.0).contains(&ratio) {
                        knob_warnings.push(format!(
                            "WARNING: TAU={:.1} TAU_SELF={:.1} Nric={} {} measured mean collisions {:.3} vs effective target {:.1} (ratio {:.3}) — KNOB IS OFF, treat this cell's physics columns with suspicion",
                            cell.tau_target, cell.tau_self, cell.n_ricochet, label, measured, target, ratio
                        ));
                    }
                }
            }
        }
        if knob_warnings.is_empty() {
            println!("\nKNOB CHECK: all active channels' measured collision counts are within [0.2, 5.0]x of min(their TAU, N_RICOCHET) — both optical-depth knobs are working.");
        } else {
            println!("\nKNOB CHECK FAILURES ({}):", knob_warnings.len());
            for w in &knob_warnings {
                println!("{w}");
            }
        }

        let find_cell = |pump: bool, tau: f64, tau_self: f64, nric: u32| -> &ReportCell {
            cells
                .iter()
                .find(|c| {
                    c.pump == pump
                        && c.tau_target == tau
                        && c.tau_self == tau_self
                        && c.n_ricochet == nric
                })
                .expect("sweep cell present")
        };
        let highlighted = [
            (
                "FORCED FLAGSHIP (v3 anchor: TAU=3, TAU_SELF=3, Nric=30)",
                find_cell(false, 3.0, 3.0, 30),
            ),
            ("PUMP NULL (TAU=0, TAU_SELF=0)", find_cell(true, 0.0, 0.0, 0)),
            (
                "PUMP FLAGSHIP (TAU=3, TAU_SELF=3, Nric=30)",
                find_cell(true, 3.0, 3.0, 30),
            ),
            (
                "PUMP HIGH-PRESSURE (TAU=3, TAU_SELF=10, Nric=30)",
                find_cell(true, 3.0, 10.0, 30),
            ),
        ];

        for (label, cell) in highlighted {
            println!("\n--- {label}: per-iteration convergence (L2 relative change of regularized density fields, emitted / ambient) ---");
            for (k, (ce, ca)) in cell.convergence.iter().enumerate() {
                println!("  iter {}: emitted {:.6}  ambient {:.6}", k + 1, ce, ca);
            }

            let print_matrix = |value_at: &dyn Fn(usize) -> f64| {
                print!("{:>10}", "r\\theta");
                for j in (0..N_THETA).step_by(2) {
                    print!(" {:>7.0}", theta_edge(j).to_degrees());
                }
                println!();
                for i in (0..N_R).step_by(3) {
                    print!("{:>10.3}", r_edge(i));
                    for j in (0..N_THETA).step_by(2) {
                        let val = value_at(i * N_THETA + j);
                        print!(" {val:>7.3}");
                    }
                    println!();
                }
            };

            // Combined matrices show the whole gas; the AMBIENT-ONLY split
            // is the decisive intake readout — near the axis the combined
            // numbers are dominated by the emitted population's outflow, so
            // only the split can show whether ambient flows DOWN the polar
            // column (the loop's intake limb).
            let sources: [(&str, Vec<&DensityField>); 2] = [
                (
                    "combined emitted+ambient",
                    vec![&cell.final_emitted_field, &cell.final_ambient_field],
                ),
                ("AMBIENT-ONLY", vec![&cell.final_ambient_field]),
            ];
            for (src_label, fields) in &sources {
                let flow = |accessor: fn(&DensityField) -> &Vec<f64>, idx: usize| -> f64 {
                    let w: f64 = fields.iter().map(|f| f.track_weight[idx]).sum();
                    let flux: f64 = fields.iter().map(|f| accessor(f)[idx]).sum();
                    if w > 1e-9 {
                        flux / w
                    } else {
                        0.0
                    }
                };

                println!("--- {label}: circulation matrix, {src_label} (rows = radial bin edge r, every 3rd bin; cols = theta bin edge in deg, every 2nd bin) ---");
                println!("    value = mean radial-flow component of {src_label} flux in that cell, bounded [-1,1]; + = net outward, - = net inward, 0 = no net radial bias.");
                print_matrix(&|idx| flow(|f| &f.radial_flux, idx));

                println!("--- {label}: meridional matrix, {src_label} (same layout) ---");
                println!("    value = mean theta-flow component (dir . theta_hat, theta_hat points SOUTHWARD/toward increasing colatitude). POLEWARD drift = negative at theta<90 and positive at theta>90 (antisymmetric about the equator); equatorward drift = the reverse.");
                print_matrix(&|idx| flow(|f| &f.theta_flux, idx));
            }

            println!("--- {label}: pressure map (same layout) ---");
            println!("    value = local collision rate per unit path an AMBIENT photon feels: sigma_e*n_e_eff + [r<=B_MAX] sigma_aa*n_a_eff — the landscape whose gradient the gas should flow down (disc ridge high, polar axis low if the loop's intake mechanism is real).");
            let f_e = CollisionField::from_tally(
                PhotonKind::Emitted,
                ParticleClass::Proton,
                cell.final_emitted_field.clone(),
            );
            let f_a = CollisionField::from_tally(
                PhotonKind::Ambient,
                ParticleClass::Proton,
                cell.final_ambient_field.clone(),
            );
            print_matrix(&|idx| {
                let i = idx / N_THETA;
                let j = idx % N_THETA;
                let r_c = (r_edge(i) * r_edge(i + 1)).sqrt();
                let th_c = 0.5 * (theta_edge(j) + theta_edge(j + 1));
                let mut lam = cell.sigma_e * f_e.effective_density_at(r_c, th_c);
                if r_c <= B_MAX {
                    lam += cell.sigma_aa * f_a.effective_density_at(r_c, th_c);
                }
                lam
            });

            // v4: the charge-pause crossover surface (user's heliopause
            // bubble; pause.html "the magnetopause is at the charge-pause").
            // Prediction if the intake funnel is real: r_c bulges at the
            // disc (emission pushes ambient back) and pinches at the poles.
            let crossover_row = |fields: &[&DensityField]| {
                for j in 0..N_THETA {
                    let mut rc: Option<f64> = None;
                    for i in (0..N_R).rev() {
                        let idx = i * N_THETA + j;
                        let w: f64 = fields.iter().map(|f| f.track_weight[idx]).sum();
                        let fl: f64 = fields.iter().map(|f| f.radial_flux[idx]).sum();
                        if w > 1e-9 && fl / w < -0.005 {
                            rc = Some((r_edge(i) * r_edge(i + 1)).sqrt());
                            break;
                        }
                    }
                    match rc {
                        Some(r) => print!(" {r:>6.2}"),
                        None => print!(" {:>6}", "-"),
                    }
                }
                println!();
            };
            println!("--- {label}: charge-pause crossover radius r_c(theta) — outermost radius whose net radial flow is INWARD (normalized < -0.005); '-' = none in column ---");
            print!("{:>10}", "theta_c");
            for j in 0..N_THETA {
                print!(
                    " {:>6.0}",
                    (0.5 * (theta_edge(j) + theta_edge(j + 1))).to_degrees()
                );
            }
            println!();
            print!("{:>10}", "combined");
            crossover_row(&[&cell.final_emitted_field, &cell.final_ambient_field]);
            print!("{:>10}", "ambient");
            crossover_row(&[&cell.final_ambient_field]);

            // v4: shell density profile — eccen.html "charge that does make
            // it back will be compressed and its density will be increased";
            // a compression ridge shows as a non-monotonic bump vs radius.
            // NOTE: the collision-field correction cap (4x analytic) bounds
            // how strongly any tally ridge can act back on transport — if a
            // ridge saturates the cap, revisit CORRECTION_MAX.
            println!("--- {label}: shell density profile (combined emitted+ambient track density per radial bin, by latitude band) ---");
            println!("{:>8} {:>12} {:>12} {:>12}", "r", "polar>60", "mid30-60", "equat<30");
            for i in 0..N_R {
                let mut tw = [0.0f64; 3];
                let mut vol = [0.0f64; 3];
                for j in 0..N_THETA {
                    let band = if !(3..15).contains(&j) {
                        0
                    } else if !(6..12).contains(&j) {
                        1
                    } else {
                        2
                    };
                    let idx = i * N_THETA + j;
                    tw[band] += cell.final_emitted_field.track_weight[idx]
                        + cell.final_ambient_field.track_weight[idx];
                    vol[band] += cell_volumes()[idx];
                }
                println!(
                    "{:>8.3} {:>12.6} {:>12.6} {:>12.6}",
                    (r_edge(i) * r_edge(i + 1)).sqrt(),
                    tw[0] / vol[0],
                    tw[1] / vol[1],
                    tw[2] / vol[2]
                );
            }

            if cell.pump {
                println!("--- {label}: pump budget trajectory (iter: ambient_absorbed + emitted_reabsorbed -> next emitted budget) ---");
                for (k, (aab, eab, next)) in cell.budget_track.iter().enumerate() {
                    println!("  iter {}: {aab:.0} + {eab:.0} -> {next}", k + 1);
                }
            }
        }

        let elapsed = started.elapsed();
        println!(
            "\nTotal report runtime: {:.1}s ({} sweep cells x {} photons/population x {} iters x 2 populations, plus iteration-0)",
            elapsed.as_secs_f64(),
            cells.len(),
            n,
            ITERS
        );
    }

    /// CM-3 neutron test: does chirality trapping + species-split directed
    /// intake give the neutron spin class a distinct polar EXIT TOPOLOGY?
    ///
    /// Setup answers "what is the neutron's profile?" without hardcoding one.
    /// Emission GEOMETRY is the same equatorial disc for both classes (every
    /// baryon-sized construction traces a disc). The only class difference is
    /// the chirality tag χ = z-sign (proton +1, neutron −1; the survey's
    /// discriminator) and the swing sign.
    ///
    /// Two coupled mechanisms run in the full photon gas + pressure + pump:
    ///   1. Chirality gears (k = cancel_backscatter, fixed 1.0 here): a
    ///      cross-collision meshes χ_marcher · χ_partner; opposite chirality
    ///      (cancel) reverses the outward run — the charge is TURNED BACK
    ///      (magmom.pdf "has to turn around and find its way out one of the
    ///      poles"). An earlier run confirmed this traps neutron charge ~2×
    ///      the proton on Earth's 2/3 field and not at all on a balanced field
    ///      — but with isotropic intake the trapped charge only reabsorbed; it
    ///      never reached a pole (the CM-3 "poles are exhaust, not intake"
    ///      wall).
    ///   2. Species-split directed intake (swept here, pole_intake_bias): the
    ///      ambient is no longer isotropic — photons stream in the SOUTH pole,
    ///      antiphotons the NORTH (venus2.pdf / neutron.pdf), building the
    ///      poleward flow the isotropic gas lacked. This is the limb that was
    ///      missing; the sweep asks whether restoring it lets the neutron's
    ///      trapped charge exit polar.
    ///
    /// Read: on Earth's 2/3 field, a distinct neutron profile shows as the
    /// bias rises — emitted polar-escape (polEsc) rising for the neutron above
    /// the proton, and the ambient polar radial flow (poleIn) turning INWARD
    /// (< 0 = intake limb closed). The balanced field is the null control
    /// (neutron.pdf: balanced → indistinguishable). bias = 0 reproduces the
    /// earlier isotropic-intake result.
    #[test]
    #[ignore = "long-running report; run with --ignored --nocapture"]
    fn report_neutron_recycling() {
        use std::time::Instant;
        let started = Instant::now();

        let n = 80_000usize;
        let swing_boost = DEFAULT_SWING_BOOST;
        // Some species-split directed intake stays on throughout (the intake
        // limb); the through-channel (cancel_redirect) is the swept knob.
        let intake_bias = 0.5_f64;
        let mut seed = 0x1BAD_C0DE_F00D_2027u64;

        println!("=== CM-3 Neutron Recycling Report (through-channel: does chirality-gated axial funnelling give the neutron a polar exit?) ===");
        println!(
            "N per population = {n}, pump flagship (TAU=3, TAU_SELF=3, Nric=30, iters={ITERS_PUMP}), swing_boost = {swing_boost}, pole_intake_bias = {intake_bias}, grid = {N_R}x{N_THETA}, R_OUT = {R_OUT}"
        );
        println!("Emission geometry IDENTICAL for both classes (bimodal disc ±7/±58); class difference = chirality χ=z_sign (proton +1, neutron −1) + swing sign.");
        println!("Sweeping cancel_redirect (through-channel): 0 = cancel does nothing (disc-only, the old wall); >0 funnels a CANCELLED (chirality-blocked) photon toward the nearer pole (neutron.pdf axial channel).");
        println!("Which photons cancel is set by chirality × field imbalance — so the polar-exit RATE is emergent; only the pole direction is imposed. Proton (χ=+1) augments on the 2/3 field (disc); neutron (χ=−1) cancels (funnelled polar).");
        println!("EXIT TOPOLOGY (emitted escaped, by latitude band): equatorial = proton disc (electric charge); polar = neutron axial channel (sub-electrical along poles).");
        println!("poleIn / eqIn = mean AMBIENT radial flow near the surface (r<~2.9) in the polar / equatorial band: NEGATIVE = inward = intake limb working.\n");

        // Mean radial flow of a density field over a radial-bin range and a
        // theta-band predicate (track-weighted). Negative = net inward.
        let band_radial_flow =
            |field: &DensityField, i_lo: usize, i_hi: usize, polar: bool| -> f64 {
                let mut flux = 0.0;
                let mut w = 0.0;
                for i in i_lo..i_hi {
                    for j in 0..N_THETA {
                        let is_polar = !(3..15).contains(&j);
                        if is_polar != polar {
                            continue;
                        }
                        let idx = i * N_THETA + j;
                        flux += field.radial_flux[idx];
                        w += field.track_weight[idx];
                    }
                }
                if w > 1e-9 {
                    flux / w
                } else {
                    0.0
                }
            };

        struct Row {
            class: &'static str,
            phi_label: &'static str,
            redirect: f64,
            n_emit: usize,
            esc_frac: f64,
            reent_frac: f64,
            eq_esc: f64,
            mid_esc: f64,
            pol_esc: f64,
            pole_in: f64,
            eq_in: f64,
        }

        let classes = [
            ("proton", ParticleClass::Proton),
            ("neutron", ParticleClass::Neutron),
        ];
        let phis = [("earth-2/3", 2.0 / 3.0), ("balanced", 0.5)];
        let redirects = [0.0_f64, 0.5, 1.0];

        let mut rows: Vec<Row> = Vec::new();
        for (cname, class) in classes {
            // iteration-0 seed depends only on the class (swing sign);
            // free-streaming, so redirect / bias are irrelevant here.
            let seed_gears = GearConfig {
                enabled: true,
                photon_fraction: PHOTON_FRACTION,
                cancel_redirect: 0.0,
                pole_intake_bias: 0.0,
                class,
                spin_transfer: DEFAULT_SPIN_TRANSFER,
                // Through-charge limb OFF: this report's committed numbers were
                // measured without it, so keep them reproducible.
                through_lane_rho: 0.0,
            };
            seed = seed.wrapping_add(0x100);
            let (emitted0, ambient0) = run_iteration0(n, swing_boost, &seed_gears, &mut seed);

            for (phi_label, phi) in phis {
                for &redirect in &redirects {
                    let gears = GearConfig {
                        enabled: true,
                        photon_fraction: phi,
                        cancel_redirect: redirect,
                        pole_intake_bias: intake_bias,
                        class,
                        spin_transfer: DEFAULT_SPIN_TRANSFER,
                        through_lane_rho: 0.0,
                    };
                    seed = seed.wrapping_add(1);
                    let cell = run_self_consistent(
                        n, swing_boost, 3.0, 3.0, 30, ITERS_PUMP, true, &gears, &emitted0,
                        &ambient0, seed,
                    );
                    let e = &cell.emitted_tally;
                    let tot = e.total.max(1.0);
                    let esc = e.escaped.max(0.0);
                    let esc_den = esc.max(1.0);
                    rows.push(Row {
                        class: cname,
                        phi_label,
                        redirect,
                        n_emit: cell.n_emit_final,
                        esc_frac: esc / tot,
                        reent_frac: e.arrived / tot,
                        eq_esc: e.escaped_band[Band::Equatorial as usize] / esc_den,
                        mid_esc: e.escaped_band[Band::Mid as usize] / esc_den,
                        pol_esc: e.escaped_band[Band::Polar as usize] / esc_den,
                        // Near-surface = r bins 0..8 (r ≈ 1 → 2.9).
                        pole_in: band_radial_flow(&cell.final_ambient_field, 0, 8, true),
                        eq_in: band_radial_flow(&cell.final_ambient_field, 0, 8, false),
                    });
                }
            }
        }

        println!(
            "{:>8} {:>10} {:>6} {:>8} {:>7} {:>7} | {:>7} {:>7} {:>7} | {:>8} {:>8}",
            "class", "field", "redir", "nEmit", "escFr", "reentFr", "eqEsc%", "midEsc%",
            "polEsc%", "poleIn", "eqIn"
        );
        for r in &rows {
            println!(
                "{:>8} {:>10} {:>6.2} {:>8} {:>7.4} {:>7.4} | {:>7.3} {:>7.3} {:>7.3} | {:>8.4} {:>8.4}",
                r.class, r.phi_label, r.redirect, r.n_emit, r.esc_frac, r.reent_frac,
                r.eq_esc, r.mid_esc, r.pol_esc, r.pole_in, r.eq_in
            );
        }

        println!("\n--- Proton vs neutron contrast (equatorial-escape ratio n/p, polar-escape diff n−p, neutron polar exit fraction) ---");
        println!(
            "{:>10} {:>6} {:>12} {:>12} {:>14} {:>14}",
            "field", "redir", "eqEsc n/p", "polEsc n−p", "n polEsc%", "p polEsc%"
        );
        let find = |cls: &str, phi: &str, redirect: f64| -> &Row {
            rows.iter()
                .find(|r| {
                    r.class == cls && r.phi_label == phi && (r.redirect - redirect).abs() < 1e-9
                })
                .expect("row present")
        };
        for (phi_label, _) in phis {
            for &redirect in &redirects {
                let p = find("proton", phi_label, redirect);
                let nn = find("neutron", phi_label, redirect);
                let eq_ratio = if p.eq_esc.abs() > 1e-9 {
                    nn.eq_esc / p.eq_esc
                } else {
                    f64::NAN
                };
                println!(
                    "{:>10} {:>6.2} {:>12.3} {:>12.3} {:>14.3} {:>14.3}",
                    phi_label,
                    redirect,
                    eq_ratio,
                    nn.pol_esc - p.pol_esc,
                    nn.pol_esc,
                    p.pol_esc
                );
            }
        }
        println!("\nREAD: neutron polar profile emerges if, on earth-2/3 as redirect rises, eqEsc(n/p) drops below 1 (disc suppressed) AND polEsc(n−p) rises above 0 (axial gained). Balanced rows are the null control (neutron.pdf: balanced field → indistinguishable).");

        let elapsed = started.elapsed();
        println!(
            "\nTotal report runtime: {:.1}s ({} cells x {} photons/pop x {} pump iters x 2 pops, plus 2 iteration-0 seeds)",
            elapsed.as_secs_f64(),
            rows.len(),
            n,
            ITERS_PUMP
        );
    }

    /// CM-3 THROUGH-CHARGE: does a species-preserving pole-to-pole limb deliver
    /// the balanced-field magnetic neutrality that `report_magnetic_neutrality`
    /// found REFUTED without it?
    ///
    /// The refutation's diagnosis was architectural, not numerical: every exit
    /// the pump offers is monochiral disc re-emission, so at a balanced field
    /// there is nothing mixed for the counter-stream to cancel against and
    /// |M_n/esc| sits at 1.0-1.3 no matter what `spin_transfer` does. venus2.pdf
    /// names the missing limb — "through charge... goes straight through the body
    /// from pole to pole, avoiding lateral recycling" — and neutron.pdf spells
    /// out the consequence: photons in the south, antiphotons in the north, "then
    /// they go out the other pole", and at balance that yields "near-total spin
    /// cancellation" while charge is still recycled.
    ///
    /// The limb is now wired in `march` (see the r < R_IN branch). Its aperture
    /// is the MEASURED hole (`THROUGH_LANE_RHO` = 1/3 of body radius, from
    /// report_axial_channel_geometry) and it is CLASS-INDEPENDENT on purpose, so
    /// any proton/neutron asymmetry has to emerge.
    ///
    /// Predictions under test:
    /// 1. Balanced field: the through-exhaust is a species MIXTURE, so
    ///    `through_spin_escaped / through_photons` → ~0 while the count stays
    ///    large. This is the direct test of the refuted claim.
    /// 2. Neutron more affected than proton, EMERGENTLY: the neutron's disc
    ///    exhaust is suppressed (eqEsc 0.308 vs the proton's 0.585), so
    ///    through-charge is a larger share of its total output and its
    ///    spin-cancelling component should dominate the total magnetic sum.
    /// 3. Earth 2/3 must NOT go neutral — the imbalance leaves a residual, per
    ///    neutron.pdf's "neutrons will only tamp down the magnetic field of the
    ///    exiting photons, but they will not cancel it".
    /// 4. Null: lane 0 reproduces the no-limb numbers bit-for-bit (the branch
    ///    draws no randomness).
    #[test]
    #[ignore = "long-running report; run with --ignored --nocapture"]
    fn report_through_charge() {
        use std::time::Instant;
        let started = Instant::now();

        let n = 40_000usize;
        let swing_boost = DEFAULT_SWING_BOOST;
        let redirect = 0.5_f64;
        let intake_bias = 0.5_f64;
        let mut seed = 0x7480_0006_C0DE_2407u64;

        println!("=== CM-3 Through-Charge Limb (species-preserving pole-to-pole crossing) ===");
        println!(
            "N={n}/pop, iters={ITERS}, cancel_redirect={redirect}, pole_intake_bias={intake_bias}, spin_transfer={DEFAULT_SPIN_TRANSFER}"
        );
        println!(
            "Lane aperture is the MEASURED hole: rho < {THROUGH_LANE_RHO:.4} of body radius (report_axial_channel_geometry: proton hole 127.5 / disc 382.5, ~11% of disc area), incidence within 20 deg of the axis. CLASS-INDEPENDENT by design."
        );
        println!(
            "TOTAL magnetic output M_tot = emitted disc exhaust (monochiral) + through-charge exhaust (species MIXTURE). Only the second can cancel.\n"
        );

        struct Row {
            cname: &'static str,
            phi_label: &'static str,
            lane: f64,
            esc_emit: f64,
            m_emit: f64,
            thru_photons: f64,
            thru_crossings: f64,
            m_thru: f64,
            esc_ambient: f64,
            thru_band: [f64; 3],
        }

        let classes = [
            ("proton", ParticleClass::Proton),
            ("neutron", ParticleClass::Neutron),
        ];
        let phis = [("balanced", 0.5_f64), ("earth-2/3", 2.0 / 3.0)];
        let lanes = [0.0_f64, THROUGH_LANE_RHO];

        let mut rows: Vec<Row> = Vec::new();
        for (cname, class) in classes {
            for (phi_label, phi) in phis {
                for &lane in &lanes {
                    let gears = GearConfig {
                        enabled: true,
                        photon_fraction: phi,
                        cancel_redirect: redirect,
                        pole_intake_bias: intake_bias,
                        class,
                        spin_transfer: DEFAULT_SPIN_TRANSFER,
                        through_lane_rho: lane,
                    };
                    seed = seed.wrapping_add(0x100);
                    let (emitted0, ambient0) = run_iteration0(n, swing_boost, &gears, &mut seed);
                    seed = seed.wrapping_add(1);
                    let cell = run_self_consistent(
                        n, swing_boost, 3.0, 3.0, 30, ITERS, false, &gears, &emitted0, &ambient0,
                        seed,
                    );
                    let e = &cell.emitted_tally;
                    let a = &cell.ambient_tally;
                    rows.push(Row {
                        cname,
                        phi_label,
                        lane,
                        esc_emit: e.escaped,
                        m_emit: e.spin_escaped_sum,
                        thru_photons: a.through_photons,
                        thru_crossings: a.through_crossings,
                        m_thru: a.through_spin_escaped,
                        esc_ambient: a.escaped,
                        thru_band: a.through_escaped_band,
                    });
                }
            }
        }

        println!("--- [1] Through-charge throughput ---");
        println!(
            "{:<9} {:<10} {:>7} {:>11} {:>11} {:>11} {:>26}",
            "class", "field", "lane", "thru_phot", "crossings", "%of ambient", "thru escape [pol/mid/eq]"
        );
        for r in &rows {
            println!(
                "{:<9} {:<10} {:>7.3} {:>11.0} {:>11.0} {:>10.2}% {:>8.0}{:>9.0}{:>9.0}",
                r.cname,
                r.phi_label,
                r.lane,
                r.thru_photons,
                r.thru_crossings,
                100.0 * r.thru_photons / (n as f64),
                r.thru_band[0],
                r.thru_band[1],
                r.thru_band[2],
            );
        }

        println!("\n--- [2] THE TEST: is the through-exhaust magnetically cancelled? ---");
        println!(
            "M_thru/phot is the mean surviving spin per through-photon. ~0 = species mixture cancelled (neutron.pdf's 'near-total spin cancellation'); ~±1 = monochiral, uncancelled."
        );
        println!(
            "{:<9} {:<10} {:>7} {:>12} {:>12} {:>12} {:>12}",
            "class", "field", "lane", "M_emit/esc", "M_thru", "M_thru/phot", "M_tot/esc_all"
        );
        for r in &rows {
            let m_emit_per = r.m_emit / r.esc_emit.max(1.0);
            let m_thru_per = if r.thru_photons > 0.0 {
                r.m_thru / r.thru_photons
            } else {
                0.0
            };
            let m_tot = (r.m_emit + r.m_thru) / (r.esc_emit + r.thru_photons).max(1.0);
            println!(
                "{:<9} {:<10} {:>7.3} {:>12.4} {:>12.1} {:>12.4} {:>12.4}",
                r.cname, r.phi_label, r.lane, m_emit_per, r.m_thru, m_thru_per, m_tot
            );
        }

        println!("\n--- [3] VERDICT vs the previously-REFUTED balanced-field prediction ---");
        for (cname, _) in classes {
            for (phi_label, _) in phis {
                let off = rows
                    .iter()
                    .find(|r| r.cname == cname && r.phi_label == phi_label && r.lane == 0.0);
                let on = rows
                    .iter()
                    .find(|r| r.cname == cname && r.phi_label == phi_label && r.lane > 0.0);
                if let (Some(off), Some(on)) = (off, on) {
                    let tot_off = (off.m_emit + off.m_thru)
                        / (off.esc_emit + off.thru_photons).max(1.0);
                    let tot_on =
                        (on.m_emit + on.m_thru) / (on.esc_emit + on.thru_photons).max(1.0);
                    let thru_per = if on.thru_photons > 0.0 {
                        on.m_thru / on.thru_photons
                    } else {
                        f64::NAN
                    };
                    println!(
                        "  {cname:<9} {phi_label:<10} M_tot/esc {tot_off:>8.4} (no limb) -> {tot_on:>8.4} (limb on), through-exhaust spin/photon {thru_per:>8.4}  {}",
                        if phi_label == "balanced" && thru_per.abs() < 0.15 {
                            "<== through-exhaust CANCELLED at balance"
                        } else if phi_label == "balanced" {
                            "<== still NOT cancelled at balance"
                        } else {
                            ""
                        }
                    );
                }
            }
        }

        println!("\n--- [4] NULL: lane 0 must be bit-identical to the no-limb module ---");
        println!(
            "The through branch draws no randomness, so at lane 0 every lane-0 row above must have thru_phot = 0 and crossings = 0 exactly."
        );
        let mut null_ok = true;
        for r in rows.iter().filter(|r| r.lane == 0.0) {
            if r.thru_photons != 0.0 || r.thru_crossings != 0.0 {
                null_ok = false;
                println!(
                    "  VIOLATION {} {}: thru_phot={} crossings={}",
                    r.cname, r.phi_label, r.thru_photons, r.thru_crossings
                );
            }
        }
        println!(
            "  {}",
            if null_ok {
                "clean — no crossings recorded with the lane closed"
            } else {
                "BROKEN — the limb fired with the lane closed"
            }
        );

        println!("\n--- [5] WHY [1]-[3] ARE STARVED, and the geometric rate ---");
        {
            // A straight ambient line with impact parameter b enters the sphere
            // essentially anti-parallel to its own direction d, so a
            // through-trajectory needs BOTH b < lane AND d within 20° of the
            // axis, and the two are near-independent.
            let p_b = (THROUGH_LANE_RHO / B_MAX).powi(2);
            let p_dir = 1.0 - THROUGH_COS_MAX_ANGLE; // both caps
            println!(
                "  P(b < lane) = (lane/B_MAX)^2 = {p_b:.3e};  P(|dir.z| > cos20) = {p_dir:.3e};  joint ~ {:.3e}",
                p_b * p_dir
            );
            println!(
                "  => ~{:.1} of {n} launches expected on a through-trajectory. The lane-0.333 rows above found 0-1. THE NULL IN [2]/[3] IS STATISTICS, NOT PHYSICS — do not read it as a refutation.",
                p_b * p_dir * n as f64
            );
            // The launch-relative number is conditioned on the b <= B_MAX
            // importance ensemble and so is not physically interpretable on its
            // own. Conditioning on ARRIVAL is: of charge that actually reaches
            // the body, the lane holds (lane/R_IN)^2 by impact parameter and the
            // axial gate keeps P(dir) of that.
            let p_arrive = (R_IN / B_MAX).powi(2);
            println!(
                "  CONDITIONAL ON ARRIVAL (the interpretable rate): arrivals are {:.3} of launches, and of ARRIVING charge {:.2}% is on a through-trajectory — (lane/R_IN)^2 = {:.1}% by impact parameter, times the {:.1}% axial-direction gate.",
                p_arrive,
                100.0 * p_b * p_dir / p_arrive,
                100.0 * THROUGH_LANE_RHO.powi(2),
                100.0 * p_dir
            );
            println!(
                "  Physically: in an ISOTROPIC near-zone field through-charge is ~0.02% of the flux. venus2.pdf agrees for large bodies (\"normally a minor complication\") but claims it dominates at the nuclear level for two reasons we can separate: the pole-to-pole distance is short (geometry, already in the model) AND \"the nucleus actually channels charge along the axis... it is pushed there by charge streams\" (ACTIVE channelling = pole_intake_bias). [7] sweeps that."
            );
        }

        println!("\n--- [6] LANE-ONLY PROBE: the conditional physics, at real statistics ---");
        println!(
            "Launches ONLY on through-trajectories (b area-weighted inside the lane, direction axial with tilt <20°, entry pole set by species per neutron.pdf: photons enter SOUTH, antiphotons NORTH). This answers 'IF charge crosses, does its exhaust cancel?' without pretending to measure the rate."
        );
        println!(
            "{:<9} {:<10} {:>9} {:>9} {:>11} {:>13} {:>24}",
            "class", "field", "launched", "crossed", "escaped", "spin/escaped", "exit band [pol/mid/eq]"
        );
        let n_probe = 20_000usize;
        let mut probe_seed = 0x1A4E_0000_9001_2407u64;
        let mut probe_rows: Vec<(&str, &str, f64, f64)> = Vec::new();
        for (cname, class) in classes {
            for (phi_label, phi) in phis {
                let gears = GearConfig {
                    enabled: true,
                    photon_fraction: phi,
                    cancel_redirect: redirect,
                    pole_intake_bias: intake_bias,
                    class,
                    spin_transfer: DEFAULT_SPIN_TRANSFER,
                    through_lane_rho: THROUGH_LANE_RHO,
                };
                // Converge the gas first so the probe marches through the real
                // fields rather than analytic seeds.
                probe_seed = probe_seed.wrapping_add(0x100);
                let (emitted0, ambient0) =
                    run_iteration0(n, swing_boost, &gears, &mut probe_seed);
                probe_seed = probe_seed.wrapping_add(1);
                let cell = run_self_consistent(
                    n, swing_boost, 3.0, 3.0, 30, ITERS, false, &gears, &emitted0, &ambient0,
                    probe_seed,
                );
                let (sigma_e, sigma_a, sigma_ee, sigma_aa) =
                    (cell.sigma_e, cell.sigma_a, cell.sigma_ee, cell.sigma_aa);
                let f_e = CollisionField::from_tally(
                    PhotonKind::Emitted,
                    class,
                    cell.final_emitted_field.clone(),
                );
                let f_a = CollisionField::from_tally(
                    PhotonKind::Ambient,
                    class,
                    cell.final_ambient_field.clone(),
                );

                let mut crossed = 0.0;
                let mut escaped = 0.0;
                let mut spin_sum = 0.0;
                let mut band = [0.0f64; 3];
                for _ in 0..n_probe {
                    let chi = if xorshift64(&mut probe_seed) < phi { 1.0 } else { -1.0 };
                    // Species split: a photon (χ=+1) travels +z and so enters
                    // the SOUTH pole; an antiphoton travels −z, entering NORTH.
                    let along = chi;
                    // Small axial tilt, uniform inside the 20° acceptance.
                    let cos_t = THROUGH_COS_MAX_ANGLE
                        + (1.0 - THROUGH_COS_MAX_ANGLE) * xorshift64(&mut probe_seed);
                    let sin_t = (1.0 - cos_t * cos_t).max(0.0).sqrt();
                    let psi = xorshift64(&mut probe_seed) * TAU;
                    let dir = DVec3::new(sin_t * psi.cos(), sin_t * psi.sin(), along * cos_t)
                        .normalize();
                    // Area-weighted b inside the lane. The offset MUST be
                    // perpendicular to `dir`, not to z: offsetting in the xy
                    // plane while the direction carries a tilt of up to 20° over
                    // an R_OUT=30 lever arm displaces the line by ~10 units and
                    // the photon misses the body completely. Same construction
                    // as `spawn_ambient`: start at −d·R_OUT + e·b with e ⊥ d, so
                    // the closest approach to the origin is exactly b.
                    let b = THROUGH_LANE_RHO * xorshift64(&mut probe_seed).sqrt();
                    let helper = if dir.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
                    let t_perp = dir.cross(helper).normalize();
                    let b_perp = dir.cross(t_perp);
                    let phi_b = xorshift64(&mut probe_seed) * TAU;
                    let e = t_perp * phi_b.cos() + b_perp * phi_b.sin();
                    let start = -dir * R_OUT + e * b;
                    // Advance to the tallied shell so the march's escape test
                    // doesn't fire on the launch point.
                    let t_entry = R_OUT - (R_OUT * R_OUT - b * b).sqrt();
                    let pos = start + dir * (t_entry + 1e-9);
                    let photon = Photon {
                        pos,
                        dir,
                        weight: 1.0,
                        kind: PhotonKind::Ambient,
                        bounces: 0,
                        chirality: chi,
                        spin: chi,
                    };
                    let mut scratch = DensityField::new();
                    // An AMBIENT photon marches through the EMITTED field as its
                    // cross channel (sigma_e) and its own population as self
                    // (sigma_aa) — matching run_self_consistent's wiring. Using
                    // sigma_a here would be the emitted photon's coupling.
                    let (term, _cb, _sb, spin, through) = march(
                        photon, &f_e, sigma_e, &f_a, sigma_aa, 30, &gears, &mut probe_seed,
                        &mut scratch,
                    );
                    let _ = (sigma_a, sigma_ee);
                    if through > 0 {
                        crossed += 1.0;
                        if let Termination::Escaped { pos, .. } = term {
                            escaped += 1.0;
                            spin_sum += spin;
                            band[lat_band(arrival_latitude(pos)) as usize] += 1.0;
                        }
                    }
                }
                let per = if escaped > 0.0 { spin_sum / escaped } else { f64::NAN };
                println!(
                    "{:<9} {:<10} {:>9} {:>9.0} {:>11.0} {:>13.4} {:>8.0}{:>8.0}{:>8.0}",
                    cname, phi_label, n_probe, crossed, escaped, per, band[0], band[1], band[2]
                );
                probe_rows.push((cname, phi_label, per, escaped));
            }
        }

        println!("\n  Conditional verdict — READ THE CAVEAT:");
        println!(
            "    CAVEAT: at a BALANCED field the probe population is launched 50/50 by species, so its net spin is ~0 BEFORE it reaches the body. A near-zero exhaust there therefore does NOT demonstrate that something was cancelled — it shows how much the particle's gearing BREAKS an already-balanced mixture. The defensible claim is the comparative one: the through limb supplies a LOW-SPIN output channel (|spin/escaped| 0.08-0.45) where the disc limb supplies a monochiral +-1.0-1.6, and the NEUTRON's gearing perturbs it least."
        );
        for (cname, phi_label, per, esc) in &probe_rows {
            if *phi_label == "balanced" {
                println!(
                    "    {cname:<9} balanced: through-exhaust spin/escaped = {per:>8.4} over {esc:.0} escapes (launched mixture was net ~0)"
                );
            }
        }
        for (cname, phi_label, per, esc) in &probe_rows {
            if *phi_label == "earth-2/3" {
                println!(
                    "    {cname:<9} earth-2/3: through-exhaust spin/escaped = {per:>8.4} over {esc:.0} escapes  (must stay NONZERO — neutron.pdf: 'neutrons will only tamp down the magnetic field ... they will not cancel it')"
                );
            }
        }

        println!("\n--- [7] How much AXIAL CHANNELLING would through-charge need to matter? ---");
        println!(
            "Sweeping pole_intake_bias (venus2.pdf's \"pushed there by charge streams\") at the full ambient ensemble, counting crossings per {n} launches."
        );
        // N raised well above the sweep's 40k: at baseline only ~7 crossings per
        // 40k are expected, so a 0-vs-4 difference there is noise, not a trend.
        let n_bias = 250_000usize;
        println!(
            "  N={n_bias} for this part (at 40k the expected count is ~7, far too few to read a trend from)."
        );
        println!(
            "{:<9} {:>10} {:>12} {:>14} {:>16}",
            "class", "bias", "crossings", "% of launches", "% of arrivals"
        );
        let mut bias_seed = 0x8142_0000_5A5A_2407u64;
        for (cname, class) in classes {
            for &bias in &[0.5_f64, 2.0, 8.0] {
                let gears = GearConfig {
                    enabled: true,
                    photon_fraction: 2.0 / 3.0,
                    cancel_redirect: redirect,
                    pole_intake_bias: bias,
                    class,
                    spin_transfer: DEFAULT_SPIN_TRANSFER,
                    through_lane_rho: THROUGH_LANE_RHO,
                };
                bias_seed = bias_seed.wrapping_add(0x100);
                let (emitted0, ambient0) =
                    run_iteration0(n_bias, swing_boost, &gears, &mut bias_seed);
                bias_seed = bias_seed.wrapping_add(1);
                let cell = run_self_consistent(
                    n_bias, swing_boost, 3.0, 3.0, 30, ITERS, false, &gears, &emitted0, &ambient0,
                    bias_seed,
                );
                let c = cell.ambient_tally.through_crossings;
                let arrivals = cell.ambient_tally.arrived.max(1.0);
                println!(
                    "{:<9} {:>10.1} {:>12.0} {:>13.4}% {:>15.3}%",
                    cname,
                    bias,
                    c,
                    100.0 * c / n_bias as f64,
                    100.0 * c / arrivals
                );
            }
        }

        println!(
            "\nTotal report runtime: {:.1}s ({} sweep cells + {} probe cells + bias sweep)",
            started.elapsed().as_secs_f64(),
            rows.len(),
            probe_rows.len()
        );
    }

    /// CM-4: does a NEIGHBOUR's exhaust light up the through-charge lane?
    ///
    /// The single-particle cell measured through-charge at 0.67% of arriving
    /// charge (saturating ~4.7% under strong axial channelling) — too small to
    /// matter, exactly as venus2.pdf says for an isolated body. Its claim is
    /// that this inverts at the NUCLEAR level. Two things change there, and this
    /// report separates them: the pole-to-pole distance is short (a separation
    /// sweep) and the lane is fed by a neighbour's exhaust rather than by
    /// isotropic ambient (an emitted-population launch from a second body).
    ///
    /// The seating comes from graphene.pdf, not from convenience: protons sit
    /// EDGE-ON (equatorial disc firing along the stack), neutrons POLE-ON
    /// (through-lane open along the stack). That predicts an asymmetry nobody
    /// put in by hand — the neutron's lane receives the proton's disc, the
    /// proton's lane points sideways at nothing.
    #[test]
    #[ignore]
    fn report_axial_stack() {
        use std::time::Instant;
        let started = Instant::now();

        let n = 40_000usize;
        let swing_boost = DEFAULT_SWING_BOOST;
        let redirect = 0.5_f64;
        let intake_bias = 0.5_f64;

        println!("=== CM-4 Axial Stack: is the through-lane fed by a NEIGHBOUR? ===");
        println!(
            "N={n}/population, iters={ITERS} (feeder gas), cancel_redirect={redirect}, pole_intake_bias={intake_bias}, spin_transfer={DEFAULT_SPIN_TRANSFER}, lane={THROUGH_LANE_RHO:.4}"
        );
        println!(
            "SEATING (graphene.pdf, and already asserted for plugs in atom_scenarios): proton EDGE-ON — 'plugged in with its equator pointing down', disc fires along the stack; neutron POLE-ON — 'plugged in with its pole pointing down', channels pole to pole. Stack axis = +z."
        );
        println!(
            "BASELINE, stated precisely. The isolated cell's 0.67% was a COLLISION-FREE GEOMETRIC BOUND (impact parameter times axial gate); its MEASURED rate at these gear settings was 0.049% (proton) / 0.168% (neutron) of arrivals, because scattering deflects most lane-bound charge back out. Part [1] re-measures the isolated rate here; parts [3]-[4] then compare channels WITHIN one cell, which is the only apples-to-apples comparison available (see [7]).\n"
        );

        /// One measured configuration.
        struct StackRow {
            label: String,
            sep: f64,
            /// (arrivals at probe, on-lane crossings) from the FEEDER's emitted
            /// exhaust.
            emit: (f64, f64),
            /// Same, from the ambient background.
            amb: (f64, f64),
            /// Crossing photons by species, and their escaped spin sum/count.
            chi_plus: f64,
            chi_minus: f64,
            spin_sum: f64,
            spin_n: f64,
            /// Emitted photons the feeder itself reabsorbed (shadowing channel).
            emit_self_lost: f64,
        }

        /// Converges a single body's gas and hands back everything a stacked
        /// march needs to read densities off it.
        fn feeder_gas(
            n: usize,
            swing_boost: f64,
            gears: &GearConfig,
            seed: &mut u64,
        ) -> (CollisionField, CollisionField, f64, f64, f64, f64) {
            *seed = seed.wrapping_add(0x100);
            let (e0, a0) = run_iteration0(n, swing_boost, gears, seed);
            *seed = seed.wrapping_add(1);
            let cell = run_self_consistent(
                n, swing_boost, 3.0, 3.0, 30, ITERS, false, gears, &e0, &a0, *seed,
            );
            let f_e = CollisionField::from_tally(
                PhotonKind::Emitted,
                gears.class,
                cell.final_emitted_field.clone(),
            );
            let f_a = CollisionField::from_tally(
                PhotonKind::Ambient,
                gears.class,
                cell.final_ambient_field.clone(),
            );
            (
                f_e,
                f_a,
                cell.sigma_e,
                cell.sigma_a,
                cell.sigma_ee,
                cell.sigma_aa,
            )
        }

        // -------------------------------------------------------------
        // COUNT DISCIPLINE (the v1 lesson, and the trap the through-charge
        // report fell into first time). An isotropic launch puts only
        // (R_IN/B_MAX)^2 of photons on the body at all, and ~0.67% of those on
        // the lane, so 40k launches yield ~7 crossings — a number no rate can
        // be read off. The isotropic channels therefore get their own much
        // larger N, and the expected count is printed so a starved row cannot
        // be mistaken for physics. The FEEDER channel needs no such help: it is
        // aimed, so its arrival fraction is the probe's solid angle (~4%), two
        // orders up.
        let n_iso = 600_000usize;
        let p_arrive_iso = (R_IN / B_MAX).powi(2);
        println!("--- [1] BASELINE: the probe ALONE, ambient feed only ---");
        println!(
            "  N={n_iso} here (not {n}): collision-free geometry predicts arrivals at {:.3} of launches => ~{:.0}, and at the MEASURED ~0.1% lane rate ~{:.0} crossings. At {n} that would be ~{:.0} — unreadable, which is the trap this part exists to avoid.",
            p_arrive_iso,
            p_arrive_iso * n_iso as f64,
            0.001 * p_arrive_iso * n_iso as f64,
            0.001 * p_arrive_iso * n as f64,
        );
        println!(
            "  Expect the measured arrivals to come in WELL UNDER that geometric figure: it assumes straight lines, while at TAU=3 with the gears live most ambient charge is scattered off its aim before reaching the surface. The shortfall is the scattering loss, not a bug — the same loss that puts the real lane rate ~7x below the 0.67% bound."
        );
        println!(
            "{:<22} {:>10} {:>11} {:>13} {:>16}",
            "probe", "launched", "arrivals", "crossings", "% of arrivals"
        );
        let mut base_rate: Vec<(&str, f64)> = Vec::new();
        for (pname, pclass) in [("neutron", ParticleClass::Neutron), ("proton", ParticleClass::Proton)]
        {
            let orient = StackOrient::nuclear_for(pclass);
            let gears = GearConfig {
                enabled: true,
                photon_fraction: 2.0 / 3.0,
                cancel_redirect: redirect,
                pole_intake_bias: intake_bias,
                class: pclass,
                spin_transfer: DEFAULT_SPIN_TRANSFER,
                through_lane_rho: THROUGH_LANE_RHO,
            };
            let mut seed = 0x5741_0000_C0DE_2407u64 ^ ((pclass as u64) << 40);
            let (f_e, f_a, sigma_e, _sigma_a, _sigma_ee, sigma_aa) =
                feeder_gas(n, swing_boost, &gears, &mut seed);
            let bodies = [StackBody {
                center: DVec3::ZERO,
                class: pclass,
                orient,
            }];
            let mut arrivals = 0.0f64;
            let mut crossings = 0.0f64;
            for _ in 0..n_iso {
                let p = spawn_ambient_on_stack(&mut seed, &gears, DVec3::ZERO);
                let out = march_stack(
                    p, &bodies, 0, 0, &f_e, sigma_e, &f_a, sigma_aa, 30, &gears, &mut seed,
                );
                if out.probe_arrived {
                    arrivals += 1.0;
                }
                if out.probe_crossed > 0 {
                    crossings += 1.0;
                }
            }
            let rate = 100.0 * crossings / arrivals.max(1.0);
            println!(
                "{:<22} {:>10} {:>11.0} {:>13.0} {:>15.3}%",
                format!("{pname} ({})", orient.label()),
                n_iso,
                arrivals,
                crossings,
                rate
            );
            println!(
                "    {crossings:.0} crossings is a SMALL COUNT: Poisson +-{:.0} => rate {rate:.3}% +-{:.3}%. Read this as 'order 0.1%', nothing finer. Consistent with the isolated cell's measured 0.049-0.168%; the collision-free 0.67% bound is ~7x above it, which is the scattering loss.",
                crossings.sqrt(),
                rate / crossings.max(1.0).sqrt()
            );
            base_rate.push((pname, rate));
        }

        // -------------------------------------------------------------
        println!("\n--- [2] NULL CONTROL: lane closed, stack present ---");
        {
            let gears_open = GearConfig {
                enabled: true,
                photon_fraction: 2.0 / 3.0,
                cancel_redirect: redirect,
                pole_intake_bias: intake_bias,
                class: ParticleClass::Proton,
                spin_transfer: DEFAULT_SPIN_TRANSFER,
                through_lane_rho: THROUGH_LANE_RHO,
            };
            let gears_shut = GearConfig {
                through_lane_rho: 0.0,
                ..gears_open
            };
            let mut seed = 0x9C10_0000_5EED_2407u64;
            let (f_e, f_a, _se, sigma_a, sigma_ee, _sa) =
                feeder_gas(n, swing_boost, &gears_open, &mut seed);
            let bodies = [
                StackBody {
                    center: DVec3::new(0.0, 0.0, -STACK_PITCH),
                    class: ParticleClass::Proton,
                    orient: StackOrient::EdgeOn,
                },
                StackBody {
                    center: DVec3::ZERO,
                    class: ParticleClass::Neutron,
                    orient: StackOrient::PoleOn,
                },
            ];
            let mut shut_cross = 0.0;
            let mut shut_seed = seed;
            for _ in 0..n {
                let p = spawn_emitted_on(&mut shut_seed, swing_boost, &bodies[0]);
                let out = march_stack(
                    p, &bodies, 0, 1, &f_a, sigma_a, &f_e, sigma_ee, 30, &gears_shut,
                    &mut shut_seed,
                );
                shut_cross += out.probe_crossed as f64;
            }
            println!(
                "  lane=0.000 over {n} feeder-emitted launches: {shut_cross:.0} crossings — {}",
                if shut_cross == 0.0 {
                    "clean"
                } else {
                    "BROKEN: the limb fired with the lane closed"
                }
            );
        }

        // -------------------------------------------------------------
        println!("\n--- [3] THE STACK: proton feeder -> probe, separation sweep ---");
        let n_amb = 400_000usize;
        println!(
            "Two populations per cell: {n} from the FEEDER's own emitted exhaust (the neighbour channel venus2.pdf points at) and {n_amb} ambient (the isolated channel, at the larger count for the reason in [1])."
        );
        println!(
            "They are reported ONLY per-channel and per-arrival. A combined rate would weight the two by their launch counts, which are a sampling choice, not the physical flux ratio — that ratio is not something this cell knows."
        );

        let mut rows: Vec<StackRow> = Vec::new();
        let feeder_class = ParticleClass::Proton;
        let mut seed = 0x2B0D_0000_57AC_2407u64;
        let gears = GearConfig {
            enabled: true,
            photon_fraction: 2.0 / 3.0,
            cancel_redirect: redirect,
            pole_intake_bias: intake_bias,
            class: feeder_class,
            spin_transfer: DEFAULT_SPIN_TRANSFER,
            through_lane_rho: THROUGH_LANE_RHO,
        };
        // gears.class drives both the emitted chirality and the collision
        // partner sampling, so it must be the FIELD body's class.
        assert_eq!(gears.class, feeder_class);
        let (f_e, f_a, sigma_e, sigma_a, sigma_ee, sigma_aa) =
            feeder_gas(n, swing_boost, &gears, &mut seed);

        struct Config {
            label: &'static str,
            feeder_orient: StackOrient,
            probe_class: ParticleClass,
            probe_orient: StackOrient,
            sep: f64,
        }
        let mut configs: Vec<Config> = Vec::new();
        for sep in [2.2_f64, STACK_PITCH, 5.2] {
            configs.push(Config {
                label: "p(edge)->n(pole)  NUCLEAR",
                feeder_orient: StackOrient::EdgeOn,
                probe_class: ParticleClass::Neutron,
                probe_orient: StackOrient::PoleOn,
                sep,
            });
        }
        // The nuclear proton's own lane points sideways — it should receive
        // almost nothing from an axial neighbour. Emergent, not imposed.
        configs.push(Config {
            label: "p(edge)->p(edge)  NUCLEAR",
            feeder_orient: StackOrient::EdgeOn,
            probe_class: ParticleClass::Proton,
            probe_orient: StackOrient::EdgeOn,
            sep: STACK_PITCH,
        });
        // Naive control: both poles on the stack axis, which is what a stack
        // test would assume WITHOUT graphene.pdf.
        configs.push(Config {
            label: "p(pole)->n(pole)  control",
            feeder_orient: StackOrient::PoleOn,
            probe_class: ParticleClass::Neutron,
            probe_orient: StackOrient::PoleOn,
            sep: STACK_PITCH,
        });

        for cfg in &configs {
            let bodies = [
                StackBody {
                    center: DVec3::new(0.0, 0.0, -cfg.sep),
                    class: feeder_class,
                    orient: cfg.feeder_orient,
                },
                StackBody {
                    center: DVec3::ZERO,
                    class: cfg.probe_class,
                    orient: cfg.probe_orient,
                },
            ];
            let centroid = (bodies[0].center + bodies[1].center) * 0.5;

            let mut emit = (0.0, 0.0);
            let mut amb = (0.0, 0.0);
            let mut chi_plus = 0.0;
            let mut chi_minus = 0.0;
            let mut spin_sum = 0.0;
            let mut spin_n = 0.0;
            let mut emit_self_lost = 0.0;

            for _ in 0..n {
                let p = spawn_emitted_on(&mut seed, swing_boost, &bodies[0]);
                let chi = p.chirality;
                let out = march_stack(
                    p, &bodies, 0, 1, &f_a, sigma_a, &f_e, sigma_ee, 30, &gears, &mut seed,
                );
                if out.probe_arrived {
                    emit.0 += 1.0;
                }
                if out.probe_crossed > 0 {
                    emit.1 += 1.0;
                    if chi >= 0.0 {
                        chi_plus += 1.0
                    } else {
                        chi_minus += 1.0
                    }
                    if let Termination::Escaped { .. } = out.term {
                        spin_sum += out.spin;
                        spin_n += 1.0;
                    }
                } else if out.other_absorbed {
                    emit_self_lost += 1.0;
                }
            }
            for _ in 0..n_amb {
                let p = spawn_ambient_on_stack(&mut seed, &gears, centroid);
                let chi = p.chirality;
                let out = march_stack(
                    p, &bodies, 0, 1, &f_e, sigma_e, &f_a, sigma_aa, 30, &gears, &mut seed,
                );
                if out.probe_arrived {
                    amb.0 += 1.0;
                }
                if out.probe_crossed > 0 {
                    amb.1 += 1.0;
                    if chi >= 0.0 {
                        chi_plus += 1.0
                    } else {
                        chi_minus += 1.0
                    }
                    if let Termination::Escaped { .. } = out.term {
                        spin_sum += out.spin;
                        spin_n += 1.0;
                    }
                }
            }
            rows.push(StackRow {
                label: cfg.label.to_string(),
                sep: cfg.sep,
                emit,
                amb,
                chi_plus,
                chi_minus,
                spin_sum,
                spin_n,
                emit_self_lost,
            });
        }

        println!(
            "\n{:<28} {:>5} {:>19} {:>21}",
            "config", "sep", "FEEDER arr/cross", "AMBIENT arr/cross"
        );
        for r in &rows {
            println!(
                "{:<28} {:>5.1} {:>10.0}/{:<8.0} {:>11.0}/{:<9.0}",
                r.label, r.sep, r.emit.0, r.emit.1, r.amb.0, r.amb.1
            );
        }

        println!("\n--- [4] Per-channel lane rate: which feed is efficient? ---");
        println!(
            "The neighbour channel is AIMED, the ambient channel isotropic. venus2's claim predicts the FEEDER column runs far above the ambient one per arrival. 'reach' is how much of the feeder's whole output the probe intercepts at all — the coupling strength, which decides whether an efficient lane rate matters."
        );
        println!(
            "{:<28} {:>5} {:>14} {:>15} {:>15} {:>13}",
            "config", "sep", "feeder reach", "feeder %/arr", "ambient %/arr", "enhancement"
        );
        for r in &rows {
            let fe = 100.0 * r.emit.1 / r.emit.0.max(1.0);
            let am = 100.0 * r.amb.1 / r.amb.0.max(1.0);
            println!(
                "{:<28} {:>5.1} {:>13.2}% {:>14.3}% {:>14.3}% {:>12.1}x",
                r.label,
                r.sep,
                100.0 * r.emit.0 / n as f64,
                fe,
                am,
                if am > 0.0 { fe / am } else { f64::NAN }
            );
        }

        println!(
            "\n  Why the rate FALLS with separation, and what it converges to. For a collimated beam the lane is a pure area gate: (lane/R_IN)^2 = {:.1}% of the geometric cross-section, and a distant feeder is collimated, so {:.1}% is the far-field CEILING. A NEAR feeder beats it because it illuminates the polar cap preferentially — the near pole sits at (sep - R_IN) while the rim of the illuminated cap sits at sqrt(sep^2 - R_IN^2), so 1/r^2 alone weights the pole by {:.2}x at sep {STACK_PITCH}, before obliquity. That is venus2.pdf's 'the pole-to-pole distance is short' clause, as a measured geometric effect rather than an assertion.",
            100.0 * THROUGH_LANE_RHO.powi(2),
            100.0 * THROUGH_LANE_RHO.powi(2),
            (STACK_PITCH * STACK_PITCH - R_IN * R_IN) / ((STACK_PITCH - R_IN) * (STACK_PITCH - R_IN))
        );

        println!("\n--- [5] Species mix and exhaust spin of the probe's through-charge ---");
        println!(
            "Monochiral through-charge cannot cancel anything: the feeder's emitted exhaust is all one species, so a lane fed ONLY by the neighbour re-imports the same problem the disc limb had. chi+/chi- is the test."
        );
        println!(
            "{:<28} {:>5} {:>9} {:>9} {:>12} {:>14}",
            "config", "sep", "chi+", "chi-", "mix (min/tot)", "spin/escaped"
        );
        for r in &rows {
            let tot = r.chi_plus + r.chi_minus;
            let mix = if tot > 0.0 {
                r.chi_plus.min(r.chi_minus) / tot
            } else {
                f64::NAN
            };
            let per = if r.spin_n > 0.0 {
                r.spin_sum / r.spin_n
            } else {
                f64::NAN
            };
            println!(
                "{:<28} {:>5.1} {:>9.0} {:>9.0} {:>12.3} {:>14.4}",
                r.label, r.sep, r.chi_plus, r.chi_minus, mix, per
            );
        }

        println!("\n--- [6] Shadowing: what the neighbour costs the feeder ---");
        println!(
            "Emitted charge the FEEDER reabsorbed (its own surface, plus anything the probe blocked before the lane) — a channel the one-body cell cannot charge for."
        );
        for r in &rows {
            println!(
                "  {:<28} sep {:>4.1}: {:>7.0} of {n} feeder-emitted lost to a non-probe body ({:>5.2}%)",
                r.label,
                r.sep,
                r.emit_self_lost,
                100.0 * r.emit_self_lost / n as f64
            );
        }

        println!("\n--- [7] VERDICT ---");
        for (pname, rate) in &base_rate {
            println!("  isolated baseline, {pname} probe: {rate:.3}% of arrivals on-lane");
        }
        let nuclear = rows
            .iter()
            .find(|r| r.label.starts_with("p(edge)->n(pole)") && r.sep == STACK_PITCH);
        if let (Some(nr), Some((_, base))) = (nuclear, base_rate.first()) {
            let fe = 100.0 * nr.emit.1 / nr.emit.0.max(1.0);
            let am = 100.0 * nr.amb.1 / nr.amb.0.max(1.0);
            println!(
                "  WHICH RATIO IS DEFENSIBLE. Feeder channel {fe:.3}%/arrival vs the ambient channel MEASURED IN THE SAME CELL {am:.3}%/arrival = {:.0}x. That is the honest enhancement: one field configuration, one geometry, two launch populations.",
                fe / am.max(1e-9)
            );
            println!(
                "  The tempting {:.0}x against the SOLO baseline ({base:.3}%) is NOT defensible — in the solo cell the probe owns the gas, so its own field scatters arriving charge out of its lane, while in the stack cells the gas belongs to the FEEDER and the probe has none (Phase-1 simplification, see march_stack). That difference alone moves the ambient rate from {base:.3}% to {am:.3}%, so the cross-configuration ratio is mostly an artifact of who owns the gas.",
                fe / base.max(1e-9)
            );
            println!(
                "  ABSOLUTE THROUGHPUT — the number that decides whether the limb can offset the disc: the probe intercepts {:.2}% of the feeder's total output and passes {:.0} of {n} feeder-emitted photons through its lane = {:.2}% of the feeder's WHOLE exhaust. In the isolated cell the lane carried ~0.001% of launched ambient. A neighbour therefore turns through-charge from negligible into a percent-level limb.",
                100.0 * nr.emit.0 / n as f64,
                nr.emit.1,
                100.0 * nr.emit.1 / n as f64
            );
            let proton_probe = rows
                .iter()
                .find(|r| r.label.starts_with("p(edge)->p(edge)"));
            if let Some(pp) = proton_probe {
                println!(
                    "  SEATING ASYMMETRY, emergent — the lane gate is class-INDEPENDENT, so nothing here was told to prefer neutrons: the edge-on PROTON probe passes {:.3}% of its {:.0} arrivals ({:.0} crossings) against the pole-on neutron's {fe:.3}%. graphene.pdf's 'protons channel pole to equator, neutrons pole to pole' falls out of the seating instead of being assumed.",
                    100.0 * pp.emit.1 / pp.emit.0.max(1.0),
                    pp.emit.0,
                    pp.emit.1
                );
            }
            let control = rows.iter().find(|r| r.label.starts_with("p(pole)"));
            if let Some(c) = control {
                println!(
                    "  AND THE SEATING IS WHAT COUPLES THE STACK AT ALL: a pole-on feeder reaches the probe with only {:.2}% of its output vs the edge-on {:.2}% ({:.0}x), because a proton's poles are nearly dark — its output is the disc. Stacking protons pole-on, the configuration a test would assume WITHOUT graphene.pdf, gives almost no coupling.",
                    100.0 * c.emit.0 / n as f64,
                    100.0 * nr.emit.0 / n as f64,
                    nr.emit.0 / c.emit.0.max(1.0)
                );
            }
            let mix = nr.chi_plus.min(nr.chi_minus) / (nr.chi_plus + nr.chi_minus).max(1.0);
            let per = nr.spin_sum / nr.spin_n.max(1.0);
            println!(
                "  BUT NEUTRALITY GOES BACKWARDS. The well-fed lane is MONOCHIRAL (species mix {mix:.3}) because a proton's exhaust is all one species by construction, and its exhaust spin per escaped photon is {per:.3} — ABOVE the isolated limb's 0.27-0.45, since a crossing photon keeps its spin and then stacks more of the same sign in the feeder's monochiral gas. Feeding the lane from a single proton neighbour fixes the RATE problem and makes the SPECIES problem worse."
            );
            println!(
                "  So a p->n pair cannot be the neutrality mechanism either. The lane needs opposite-species feeds, which is what ammon.pdf describes for a real nucleus: 'charge moving pole-to-pole through the alphas, south to north and north to south... giving us both charge and anticharge'. In this model the neutron's OWN emission is the anticharge (chirality_sign -1) while the protons' is charge, so the candidate cancellation is between the neutron's through-charge (+) and its own emission (-) — a whole-particle sum, not a within-limb one. That is a 3-body p-n-p cell and a different observable."
            );
        }
        println!(
            "\nTotal report runtime: {:.1}s ({} stack cells + {} baselines)",
            started.elapsed().as_secs_f64(),
            rows.len(),
            base_rate.len()
        );
    }

    /// CM-4b: the 3-body p–n–p cell, and the observable CM-4 could not ask.
    ///
    /// CM-4's refutation was that any single-species feed makes the lane MORE
    /// monochiral. ammon.pdf's picture is charge running through the stack both
    /// ways, "giving us both charge and anticharge" — and in this model the
    /// neutron's OWN emission is the anticharge (chirality −1) while the proton
    /// feeders' through-charge is +1. So the candidate cancellation is a
    /// WHOLE-PARTICLE sum over everything leaving the middle body (voyag.pdf:
    /// "the spins will offset as a sum"): own emission plus through-exhaust,
    /// NOT the within-lane mix that CM-4 measured.
    ///
    /// Phase-1 limits, stated up front:
    /// - one converged single-body gas per marched population (the field
    ///   body's); the probe still owns no gas, so through rates stay UPPER
    ///   bounds. Part [5] integrates the optical depth Phase 1 skips, so the
    ///   size of that error is printed rather than guessed at.
    /// - bodies are STATIC. The stream's push on a neighbour is tallied in [4]
    ///   (absorbed photon directions, cc.pdf open-field absorber force) but
    ///   nothing moves; matching Atom mode's stream-cushion force channel is
    ///   the eventual consumer of that number.
    /// - equal per-body emission budgets are assumed when the own and feeder
    ///   populations are combined (all nucleons recycle ~19x their mass/sec,
    ///   so to first order budgets scale with mass ≈ equal for p/n).
    #[test]
    #[ignore]
    fn report_pnp_stack() {
        use std::time::Instant;
        let started = Instant::now();

        let n_gas = 40_000usize;
        let n_own = 40_000usize;
        let n_feed = 80_000usize;
        // Equal-budget weight: a feeder launch of n_feed represents the same
        // physical budget as the middle body's n_own.
        let w = n_own as f64 / n_feed as f64;
        let swing_boost = DEFAULT_SWING_BOOST;
        let redirect = 0.5_f64;
        let intake_bias = 0.5_f64;

        let gears_for = |class: ParticleClass| GearConfig {
            enabled: true,
            photon_fraction: 2.0 / 3.0,
            cancel_redirect: redirect,
            pole_intake_bias: intake_bias,
            class,
            spin_transfer: DEFAULT_SPIN_TRANSFER,
            through_lane_rho: THROUGH_LANE_RHO,
        };
        let gears_p = gears_for(ParticleClass::Proton);
        let gears_n = gears_for(ParticleClass::Neutron);

        println!("=== CM-4b p-n-p Stack: does the WHOLE-PARTICLE spin sum move toward zero? ===");
        println!(
            "gas N={n_gas} x{ITERS} iters/class, own N={n_own}, feeder N={n_feed} (budget weight {w:.2}), cancel_redirect={redirect}, pole_intake_bias={intake_bias}, spin_transfer={DEFAULT_SPIN_TRANSFER}, lane={THROUGH_LANE_RHO:.4}"
        );
        println!(
            "OBSERVABLE: M_total/escaped over everything leaving the MIDDLE body — its own emission (neutron: chirality -1, the anticharge) plus the through-exhaust fed by its neighbours (proton: +1). CM-4 measured only the lane and found it monochiral; the sum is the different question ammon.pdf/voyag.pdf actually pose."
        );
        println!(
            "CONTROLS: n-n-n feeds the lane the SAME species as the middle body's own emission, so its sum must move AWAY from zero if species opposition is what matters; lane-shut p-n-p must reproduce the isolated sum exactly.\n"
        );

        /// Converges a single body's gas (same recipe as report_axial_stack).
        fn gas_for(
            n: usize,
            swing_boost: f64,
            gears: &GearConfig,
            seed: &mut u64,
        ) -> (CollisionField, CollisionField, f64, f64, f64, f64) {
            *seed = seed.wrapping_add(0x100);
            let (e0, a0) = run_iteration0(n, swing_boost, gears, seed);
            *seed = seed.wrapping_add(1);
            let cell = run_self_consistent(
                n, swing_boost, 3.0, 3.0, 30, ITERS, false, gears, &e0, &a0, *seed,
            );
            let f_e = CollisionField::from_tally(
                PhotonKind::Emitted,
                gears.class,
                cell.final_emitted_field.clone(),
            );
            let f_a = CollisionField::from_tally(
                PhotonKind::Ambient,
                gears.class,
                cell.final_ambient_field.clone(),
            );
            (f_e, f_a, cell.sigma_e, cell.sigma_a, cell.sigma_ee, cell.sigma_aa)
        }

        let mut seed = 0x7E57_0000_9A5C_2607u64;
        let (fe_p, fa_p, _se_p, sa_p, see_p, _saa_p) =
            gas_for(n_gas, swing_boost, &gears_p, &mut seed);
        let (fe_n, fa_n, _se_n, sa_n, see_n, _saa_n) =
            gas_for(n_gas, swing_boost, &gears_n, &mut seed);

        let body = |z: f64, class: ParticleClass| StackBody {
            center: DVec3::new(0.0, 0.0, z),
            class,
            orient: StackOrient::nuclear_for(class),
        };
        let p_at = |z: f64| body(z, ParticleClass::Proton);
        let n_at = |z: f64| body(z, ParticleClass::Neutron);

        struct PnpConfig {
            label: &'static str,
            bodies: Vec<StackBody>,
            middle: usize,
            feeders: Vec<usize>,
            lane_open: bool,
        }
        let configs = [
            PnpConfig {
                label: "n alone",
                bodies: vec![n_at(0.0)],
                middle: 0,
                feeders: vec![],
                lane_open: true,
            },
            PnpConfig {
                label: "p-n (1 feeder)",
                bodies: vec![p_at(-STACK_PITCH), n_at(0.0)],
                middle: 1,
                feeders: vec![0],
                lane_open: true,
            },
            PnpConfig {
                label: "p-n-p (2 feeders)",
                bodies: vec![p_at(-STACK_PITCH), n_at(0.0), p_at(STACK_PITCH)],
                middle: 1,
                feeders: vec![0, 2],
                lane_open: true,
            },
            PnpConfig {
                label: "n-n-n SAME-SPECIES",
                bodies: vec![n_at(-STACK_PITCH), n_at(0.0), n_at(STACK_PITCH)],
                middle: 1,
                feeders: vec![0, 2],
                lane_open: true,
            },
            PnpConfig {
                label: "p-n-p LANE SHUT",
                bodies: vec![p_at(-STACK_PITCH), n_at(0.0), p_at(STACK_PITCH)],
                middle: 1,
                feeders: vec![0, 2],
                lane_open: false,
            },
        ];

        struct PnpRow {
            label: &'static str,
            /// Own-emission channel of the middle body.
            own_esc: f64,
            own_m: f64,
            own_shadowed: f64,
            /// Through channel, summed over feeders, at feeder budget n_feed each.
            reach: f64,
            crossings: f64,
            thru_esc: f64,
            thru_m: f64,
            /// Feeder charge the middle body absorbed (its intake).
            fed_absorbed: f64,
            /// Net z-momentum the middle body ABSORBED from feeder exhaust,
            /// summed over feeders (units: photon momenta, per n_feed budget each).
            mom_mid_z: f64,
            /// Same for each feeder body itself (self-reabsorption + opposite
            /// feeder's beam), keyed south (index 0) / north.
            mom_feeders_z: f64,
            /// Through species tally of crossed-and-escaped photons.
            thru_chi_plus: f64,
            thru_chi_minus: f64,
        }

        let mut rows: Vec<PnpRow> = Vec::new();
        for (ci, cfg) in configs.iter().enumerate() {
            let shut_p;
            let shut_n;
            let (gp, gn) = if cfg.lane_open {
                (&gears_p, &gears_n)
            } else {
                shut_p = GearConfig { through_lane_rho: 0.0, ..gears_p };
                shut_n = GearConfig { through_lane_rho: 0.0, ..gears_n };
                (&shut_p, &shut_n)
            };
            let mut seed = 0x2607_0000_0B0D_1E5Eu64 ^ ((ci as u64) << 32);

            // --- Own-emission pass: field body = probe = the middle body. ---
            let mut own_esc = 0.0;
            let mut own_m = 0.0;
            let mut own_shadowed = 0.0;
            for _ in 0..n_own {
                let p = spawn_emitted_on(&mut seed, swing_boost, &cfg.bodies[cfg.middle]);
                let out = march_stack(
                    p, &cfg.bodies, cfg.middle, cfg.middle, &fa_n, sa_n, &fe_n, see_n, 30,
                    gn, &mut seed,
                );
                if let Termination::Escaped { .. } = out.term {
                    own_esc += 1.0;
                    own_m += out.spin;
                }
                if out.other_absorbed {
                    own_shadowed += 1.0;
                }
            }

            // --- Feeder passes: field body = that feeder, probe = middle. ---
            let mut reach = 0.0;
            let mut crossings = 0.0;
            let mut thru_esc = 0.0;
            let mut thru_m = 0.0;
            let mut fed_absorbed = 0.0;
            let mut mom_mid_z = 0.0;
            let mut mom_feeders_z = 0.0;
            let mut thru_chi_plus = 0.0;
            let mut thru_chi_minus = 0.0;
            for &fi in &cfg.feeders {
                let fclass = cfg.bodies[fi].class;
                let (g, fa, sa, fe, see) = match fclass {
                    ParticleClass::Proton => (gp, &fa_p, sa_p, &fe_p, see_p),
                    _ => (gn, &fa_n, sa_n, &fe_n, see_n),
                };
                for _ in 0..n_feed {
                    let p = spawn_emitted_on(&mut seed, swing_boost, &cfg.bodies[fi]);
                    let chi = p.chirality;
                    let out = march_stack(
                        p, &cfg.bodies, fi, cfg.middle, fa, sa, fe, see, 30, g, &mut seed,
                    );
                    if out.probe_arrived {
                        reach += 1.0;
                    }
                    if out.probe_crossed > 0 {
                        crossings += 1.0;
                        if let Termination::Escaped { .. } = out.term {
                            thru_esc += 1.0;
                            thru_m += out.spin;
                            if chi >= 0.0 {
                                thru_chi_plus += 1.0;
                            } else {
                                thru_chi_minus += 1.0;
                            }
                        }
                    }
                    if let (Some(bi), Termination::Absorbed { dir, .. }) =
                        (out.absorbed_by, out.term)
                    {
                        if bi == cfg.middle {
                            fed_absorbed += 1.0;
                            mom_mid_z += dir.z;
                        } else {
                            mom_feeders_z += dir.z;
                        }
                    }
                }
            }

            rows.push(PnpRow {
                label: cfg.label,
                own_esc,
                own_m,
                own_shadowed,
                reach,
                crossings,
                thru_esc,
                thru_m,
                fed_absorbed,
                mom_mid_z,
                mom_feeders_z,
                thru_chi_plus,
                thru_chi_minus,
            });
        }

        // -------------------------------------------------------------
        println!("--- [1] The whole-particle sum, per config ---");
        println!(
            "own = the middle body's emission marched through the stack; thru = feeder photons that crossed the middle lane and then escaped, scaled x{w:.2} to the same budget. M values are spin sums over escaped photons."
        );
        println!(
            "{:<22} {:>9} {:>10} {:>7} {:>9} {:>10} {:>13} {:>13}",
            "config", "own esc", "own M/esc", "cross", "thru esc", "thru M/esc", "TOTAL M/esc", "|shift|"
        );
        let mut base_m: f64 = f64::NAN;
        for r in &rows {
            let own_per = r.own_m / r.own_esc.max(1.0);
            let thru_per = if r.thru_esc > 0.0 { r.thru_m / r.thru_esc } else { f64::NAN };
            let total = (r.own_m + w * r.thru_m) / (r.own_esc + w * r.thru_esc).max(1.0);
            if r.label == "n alone" {
                base_m = total;
            }
            println!(
                "{:<22} {:>9.0} {:>10.4} {:>7.0} {:>9.0} {:>10.4} {:>13.4} {:>13.4}",
                r.label,
                r.own_esc,
                own_per,
                r.crossings,
                r.thru_esc,
                thru_per,
                total,
                (total - base_m).abs()
            );
            if r.thru_esc > 0.0 && r.thru_esc < 100.0 {
                println!(
                    "    NOTE {:.0} through-escapes is a small count: Poisson +-{:.0}, so thru M/esc carries ~{:.0}% relative error.",
                    r.thru_esc,
                    r.thru_esc.sqrt(),
                    100.0 / r.thru_esc.max(1.0).sqrt()
                );
            }
        }
        println!("  Shadowing of the middle body's own emission by its neighbours (of {n_own} launched):");
        for r in &rows {
            println!(
                "    {:<22} {:>6.0} ({:>5.2}%) — feeder reach into the middle body: {:>7.0} of {} per feeder pass",
                r.label,
                r.own_shadowed,
                100.0 * r.own_shadowed / n_own as f64,
                r.reach,
                n_feed
            );
        }

        // -------------------------------------------------------------
        println!("\n--- [2] Dose-response and the feed multiple neutrality would need ---");
        let one = rows.iter().find(|r| r.label.starts_with("p-n (1")).unwrap();
        let two = rows.iter().find(|r| r.label.starts_with("p-n-p (2")).unwrap();
        let nnn = rows.iter().find(|r| r.label.starts_with("n-n-n")).unwrap();
        let shut = rows.iter().find(|r| r.label.ends_with("SHUT")).unwrap();
        {
            let m_of = |r: &PnpRow| (r.own_m + w * r.thru_m) / (r.own_esc + w * r.thru_esc).max(1.0);
            let d1 = m_of(one) - base_m;
            let d2 = m_of(two) - base_m;
            println!(
                "  shift per proton feeder: 1 feeder {d1:+.4}, 2 feeders {d2:+.4} (linear would be {:+.4}) — sign {} the neutron's own {base_m:+.4}.",
                2.0 * d1,
                if d1 * base_m < 0.0 { "OPPOSES" } else { "does NOT oppose" }
            );
            if nnn.thru_esc < 10.0 {
                println!(
                    "  same-species control: STARVED, not failed — pole-on neutron feeders barely couple (reach {:.0} vs the edge-on proton's {:.0} per {n_feed}), so the same-species lane got no flux to move the sum with. The control is one-sided; the species conclusion rests on the thru M/esc SIGN opposing the own sign, which needs no control.",
                    nnn.reach / 2.0,
                    two.reach / 2.0
                );
            } else {
                println!(
                    "  same-species control: n-n-n total {:+.4} vs isolated {base_m:+.4} — {}",
                    m_of(nnn),
                    if (m_of(nnn) - base_m) * base_m > 0.0 {
                        "moves AWAY from zero, as species opposition requires"
                    } else {
                        "FAILS to move away from zero — the shift is not a species effect"
                    }
                );
            }
            let cap1 = 1.0 - one.thru_esc / one.crossings.max(1.0);
            let cap2 = 1.0 - two.thru_esc / two.crossings.max(1.0);
            println!(
                "  HAND-OFF (found in the data, not designed for): of the middle-lane crossings, {:.1}% (p-n) vs {:.1}% (p-n-p) are absorbed before escaping the stack — the FAR feeder sits on the exit axis and eats the through-exhaust. That is nuclear charge hand-off along the stack, and it is why the ESCAPED-sum dose-response is not monotonic in feeder count: the second feeder both doubles the feed and captures most of what crosses.",
                100.0 * cap1,
                100.0 * cap2
            );
            println!(
                "  NOISE FLOOR for the |shift| column: the lane-shut null moved {:+.4} with ZERO crossings (seed + shadowing), so config-to-config totals carry non-lane variation of that order. The thru columns ({:.0} escapes at {:+.2}) are the clean signal; the totals are illustrative.",
                m_of(shut) - base_m,
                one.thru_esc,
                one.thru_m / one.thru_esc.max(1.0)
            );
            println!(
                "  lane-shut null: total {:+.4} vs isolated {base_m:+.4} (crossings {:.0}) — {}",
                m_of(shut),
                shut.crossings,
                if shut.crossings == 0.0 { "clean" } else { "BROKEN: lane fired while shut" }
            );
            // How much feed would zero the sum, holding the per-photon spins fixed?
            let thru_m_per_feeder = w * two.thru_m / 2.0;
            if thru_m_per_feeder.abs() > 0.0 {
                let f_star = -one.own_m / thru_m_per_feeder;
                println!(
                    "  EXTRAPOLATION (linear in feed, per-photon spins held): zeroing the sum needs f* = {f_star:.0} proton-feeder equivalents = {:.0}x the p-n-p feed. A bare pitch-2.6 pair passes ~1% of a feeder's exhaust through the lane; f* says what multiple of that flux the sum actually needs. THAT is the quantity the nuclear papers claim the dense alpha-stack field supplies — the number to test when the probe owns a gas and feeds are self-consistent (Phase 2), not a verdict on Phase 1.",
                    f_star / 2.0
                );
            }
        }

        // -------------------------------------------------------------
        println!("\n--- [3] Species mix of the middle body's TOTAL output ---");
        println!(
            "CM-4's lane-only mix was 0.031 (monochiral). The whole-particle mix counts the own emission too, budget-weighted; opposition needs BOTH species present in what leaves the body."
        );
        for r in &rows {
            // The middle body is a NEUTRON in every config here, so its own
            // escaped emission is all chirality -1 by construction.
            let plus = w * r.thru_chi_plus;
            let minus = r.own_esc + w * r.thru_chi_minus;
            let tot = plus + minus;
            println!(
                "  {:<22} chi+ {:>8.0}  chi- {:>8.0}  mix {:.4}",
                r.label,
                plus,
                minus,
                if tot > 0.0 { plus.min(minus) / tot } else { f64::NAN }
            );
        }

        // -------------------------------------------------------------
        println!("\n--- [4] Momentum: what the feeder streams PUSH on (bodies are static; this is the force a dynamic run would feel) ---");
        println!(
            "Sum of absorbed photon z-directions per body (cc.pdf open field: absorber force only, no emitter recoil — see project_momentum_open_field). Units: photon momenta per feeder budget of {n_feed}."
        );
        for r in &rows {
            if r.reach == 0.0 {
                continue;
            }
            println!(
                "  {:<22} middle absorbed {:>7.0} photons, net z-push {:>+9.1} ({:+.4}/absorbed); feeders' own bodies {:>+9.1}",
                r.label,
                r.fed_absorbed,
                r.mom_mid_z,
                r.mom_mid_z / r.fed_absorbed.max(1.0),
                r.mom_feeders_z
            );
        }
        println!(
            "  Read: p-n pushes the middle body OFF the feeder (one-sided cushion); p-n-p should cancel to ~0 net by symmetry while each side still bears the load. The through-crossing itself transfers nothing here — an uncollided crossing is a free pass by construction, which is exactly the salt.pdf gate."
        );

        // -------------------------------------------------------------
        println!("\n--- [5] What Phase 1 does NOT march: optical depth across the gap ---");
        println!(
            "A lane photon from the south feeder crosses two gas columns Phase 1 ignores: the middle body's OWN gas (probe owns no gas => rates are upper bounds) and the FAR feeder's gas (the counter-stream: collisions there would gear-redirect lane charge into the equatorial whirlpool -> a denser disc, the alpha-densification effect the nuclear papers describe). Integrated straight-line tau over each column, on-axis:"
        );
        {
            let mid = n_at(0.0);
            let far = p_at(STACK_PITCH);
            let seg = |a: f64, b: f64, body: &StackBody, field: &CollisionField, sigma: f64| {
                let ds = 0.005;
                let mut tau = 0.0;
                let mut z = a;
                while z < b {
                    let (r, theta) = body.local_polar(DVec3::new(0.0, 0.0, z));
                    tau += sigma * field.effective_density_at(r, theta) * ds;
                    z += ds;
                }
                tau
            };
            // South gap: feeder surface to probe surface; north gap: probe exit
            // to far-feeder surface.
            let tau_mid_s = seg(-STACK_PITCH + R_IN, -R_IN, &mid, &fe_n, see_n)
                + seg(-STACK_PITCH + R_IN, -R_IN, &mid, &fa_n, sa_n);
            let tau_far_s = seg(-STACK_PITCH + R_IN, -R_IN, &far, &fe_p, see_p);
            let tau_mid_n = seg(R_IN, STACK_PITCH - R_IN, &mid, &fe_n, see_n)
                + seg(R_IN, STACK_PITCH - R_IN, &mid, &fa_n, sa_n);
            let tau_far_n = seg(R_IN, STACK_PITCH - R_IN, &far, &fe_p, see_p);
            println!(
                "  south gap (feeder->probe): tau_middle-gas {tau_mid_s:.3} (survival {:.2}), tau_counter-stream {tau_far_s:.3} (survival {:.2})",
                (-tau_mid_s).exp(),
                (-tau_far_s).exp()
            );
            println!(
                "  north gap (probe->far feeder): tau_middle-gas {tau_mid_n:.3} (survival {:.2}), tau_counter-stream {tau_far_n:.3} (survival {:.2})",
                (-tau_mid_n).exp(),
                (-tau_far_n).exp()
            );
            println!(
                "  Read: 1-survival is the fraction of lane charge a self-consistent Phase 2 would strip from the lane per gap crossing and hand mostly to the DISC (gears funnel augments equatorial). It scales the Phase-1 through rates DOWN and the disc density UP — both worth keeping an eye on before trusting any Phase-1 magnitude to better than this factor."
            );
        }

        println!(
            "\nTotal report runtime: {:.1}s ({} configs, 2 converged gases)",
            started.elapsed().as_secs_f64(),
            rows.len()
        );
    }

    /// Momentum bookkeeping needs the absorber's identity: a blocked photon
    /// must name the blocker, an escaping or through-crossing photon must name
    /// nobody. Pure geometry (sigma = 0).
    #[test]
    fn stack_absorber_identity_for_momentum_tally() {
        let gears = GearConfig {
            enabled: true,
            photon_fraction: PHOTON_FRACTION,
            cancel_redirect: 0.0,
            pole_intake_bias: 0.0,
            class: ParticleClass::Proton,
            spin_transfer: 0.0,
            through_lane_rho: THROUGH_LANE_RHO,
        };
        let f_e = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);
        let f_a = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);
        let mut rng = 0xD00D_0000_ABCD_2607u64;

        let pair = [
            StackBody {
                center: DVec3::new(0.0, 0.0, -STACK_PITCH),
                class: ParticleClass::Proton,
                orient: StackOrient::EdgeOn,
            },
            StackBody {
                center: DVec3::ZERO,
                class: ParticleClass::Neutron,
                orient: StackOrient::PoleOn,
            },
        ];

        // Aimed at the edge-on feeder from outside: absorbed by body 0, and the
        // stored direction is the momentum the absorber received.
        let blocked = Photon {
            pos: DVec3::new(0.0, 0.0, -(R_OUT - 0.5)),
            dir: DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 1.0,
            spin: 0.0,
        };
        let out = march_stack(
            blocked, &pair, 0, 1, &f_e, 0.0, &f_a, 0.0, 0, &gears, &mut rng,
        );
        assert_eq!(out.absorbed_by, Some(0), "the blocker must be named");
        assert!(out.other_absorbed);
        match out.term {
            Termination::Absorbed { dir, .. } => {
                assert!((dir.z - 1.0).abs() < 1e-12, "absorbed momentum is the photon's dir")
            }
            _ => panic!("expected absorption"),
        }

        // Fired outward: escapes, no absorber.
        let escaping = Photon {
            pos: DVec3::new(0.0, 3.0, 0.0),
            dir: DVec3::Y,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 1.0,
            spin: 0.0,
        };
        let out = march_stack(
            escaping, &pair, 0, 1, &f_e, 0.0, &f_a, 0.0, 0, &gears, &mut rng,
        );
        assert!(matches!(out.term, Termination::Escaped { .. }));
        assert_eq!(out.absorbed_by, None);

        // Down the pole-on body's open lane: crosses and leaves — a free pass
        // must transfer nothing, so it names no absorber either.
        let lane = Photon {
            pos: DVec3::new(0.0, 0.0, 5.0),
            dir: -DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 1.0,
            spin: 0.0,
        };
        let out = march_stack(
            lane, &[pair[1]], 0, 0, &f_e, 0.0, &f_a, 0.0, 0, &gears, &mut rng,
        );
        assert_eq!(out.probe_crossed, 1, "lane must be open on-axis");
        assert!(matches!(out.term, Termination::Escaped { .. }));
        assert_eq!(out.absorbed_by, None, "a free pass transfers no momentum");
    }

    /// The stack's seating convention must be the one graphene.pdf states and
    /// the one Atom mode already asserts for plugs — protons edge-on, neutrons
    /// pole-on. If these ever diverge, the two modes are modelling different
    /// nuclei.
    #[test]
    fn stack_seating_matches_atom_mode_plug_convention() {
        let p = StackOrient::nuclear_for(ParticleClass::Proton);
        let nn = StackOrient::nuclear_for(ParticleClass::Neutron);
        assert_eq!(p, StackOrient::EdgeOn, "graphene.pdf: proton equator down");
        assert_eq!(nn, StackOrient::PoleOn, "graphene.pdf: neutron pole down");
        // Mirrors atom_scenarios' plug assertions (rest_axis·stack ≈ 0 for a
        // proton, ≈ 1 for a neutron).
        assert!(
            p.spin_axis().dot(DVec3::Z).abs() < 1e-12,
            "plug proton must be edge-on to the stack axis"
        );
        assert!(
            nn.spin_axis().dot(DVec3::Z).abs() > 1.0 - 1e-12,
            "plug neutron must keep its pole on the stack axis"
        );
    }

    /// The through-lane gate must be evaluated in each body's OWN frame. Same
    /// axial photon, same stack: a pole-on body passes it, an edge-on body
    /// (whose lane points sideways) absorbs it. This is the geometric asymmetry
    /// the nuclear seating predicts, so it has to be real in the code and not
    /// an artifact of how the report happens to be wired.
    #[test]
    fn stack_lane_gate_is_evaluated_in_each_bodys_own_frame() {
        let gears = GearConfig {
            enabled: true,
            photon_fraction: 2.0 / 3.0,
            cancel_redirect: 0.0,
            pole_intake_bias: 0.0,
            class: ParticleClass::Proton,
            // No collisions: this is a pure geometry test.
            spin_transfer: 0.0,
            through_lane_rho: THROUGH_LANE_RHO,
        };
        let f_e = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);
        let f_a = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);

        // A photon on the stack axis heading +z straight at the probe at origin.
        let launch = |orient: StackOrient, rng: &mut u64| {
            let bodies = [StackBody {
                center: DVec3::ZERO,
                class: ParticleClass::Neutron,
                orient,
            }];
            let photon = Photon {
                pos: DVec3::new(0.0, 0.0, -5.0),
                dir: DVec3::Z,
                weight: 1.0,
                kind: PhotonKind::Ambient,
                bounces: 0,
                chirality: 1.0,
                spin: 1.0,
            };
            march_stack(
                photon, &bodies, 0, 0, &f_e, 0.0, &f_a, 0.0, 0, &gears, rng,
            )
        };

        let mut rng = 0xA5A5_0000_1234_5678u64;
        let pole_on = launch(StackOrient::PoleOn, &mut rng);
        assert!(pole_on.probe_arrived, "photon should reach the surface");
        assert_eq!(
            pole_on.probe_crossed, 1,
            "a pole-on body's lane is open along the stack axis"
        );
        assert!(
            matches!(pole_on.term, Termination::Escaped { .. }),
            "having crossed, it should leave the shell"
        );
        // Species and spin survive the crossing — the retention that makes a
        // mixed exhaust possible.
        assert_eq!(pole_on.spin, 1.0, "crossing must not alter spin");

        let edge_on = launch(StackOrient::EdgeOn, &mut rng);
        assert!(edge_on.probe_arrived, "photon should reach the surface");
        assert_eq!(
            edge_on.probe_crossed, 0,
            "an edge-on body's lane points along x, so axial feed must be absorbed"
        );
        assert!(matches!(edge_on.term, Termination::Absorbed { .. }));
    }

    /// A one-body "stack" must reproduce the single-body cell's geometry: the
    /// stacked march is a generalization, not a different model. Checks the two
    /// things that could silently diverge — where the shell boundary sits, and
    /// that a second body genuinely shadows the first.
    #[test]
    fn stack_of_one_matches_single_body_geometry_and_two_bodies_shadow() {
        let gears = GearConfig::off();
        let f_e = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Proton);
        let f_a = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Proton);
        let mut rng = 0x0BAD_0000_C0FF_EE01u64;

        let solo = [StackBody {
            center: DVec3::ZERO,
            class: ParticleClass::Proton,
            orient: StackOrient::PoleOn,
        }];
        // Straight at the body from just inside the shell: must be absorbed,
        // never escape (with the lane shut this is the plain single-body rule).
        let aimed = Photon {
            pos: DVec3::new(0.0, 0.0, -(R_OUT - 0.5)),
            dir: DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 0.0,
            spin: 0.0,
        };
        let out = march_stack(
            aimed, &solo, 0, 0, &f_e, 0.0, &f_a, 0.0, 0, &gears, &mut rng,
        );
        assert!(out.probe_arrived);
        assert!(matches!(out.term, Termination::Absorbed { .. }));
        assert_eq!(out.probe_crossed, 0, "lane is shut in GearConfig::off()");
        assert!(!out.other_absorbed, "there is no other body");

        // Now park a blocker between the launch point and the probe. The same
        // photon must be stopped by the blocker and never reach the probe.
        let pair = [
            StackBody {
                center: DVec3::new(0.0, 0.0, -STACK_PITCH),
                class: ParticleClass::Proton,
                orient: StackOrient::PoleOn,
            },
            StackBody {
                center: DVec3::ZERO,
                class: ParticleClass::Neutron,
                orient: StackOrient::PoleOn,
            },
        ];
        let blocked = Photon {
            pos: DVec3::new(0.0, 0.0, -(R_OUT - 0.5)),
            dir: DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: 0.0,
            spin: 0.0,
        };
        let out2 = march_stack(
            blocked, &pair, 0, 1, &f_e, 0.0, &f_a, 0.0, 0, &gears, &mut rng,
        );
        assert!(
            !out2.probe_arrived,
            "the blocker must shadow the probe completely on-axis"
        );
        assert!(out2.other_absorbed, "the blocker should have absorbed it");
    }

    /// The through-charge limb must (a) be a strict no-op with the lane closed,
    /// (b) actually pass charge when a photon is aimed down the lane, and (c)
    /// preserve the crossing photon's SPECIES — the retention is the whole point,
    /// since monochiral disc re-emission is what made balanced-field neutrality
    /// unreachable.
    #[test]
    fn through_charge_limb_crosses_and_preserves_species() {
        let gears_on = GearConfig {
            enabled: true,
            photon_fraction: PHOTON_FRACTION,
            cancel_redirect: 0.5,
            pole_intake_bias: 0.5,
            class: ParticleClass::Neutron,
            spin_transfer: 0.0, // isolate the limb from the meshing channel
            through_lane_rho: THROUGH_LANE_RHO,
        };
        let gears_off = GearConfig {
            through_lane_rho: 0.0,
            ..gears_on
        };
        let pure_e = CollisionField::pure_analytic(PhotonKind::Emitted, ParticleClass::Neutron);
        let pure_a = CollisionField::pure_analytic(PhotonKind::Ambient, ParticleClass::Neutron);

        // Straight down the axis, well inside the lane, no collisions. NOTE the
        // start radius must be strictly INSIDE R_OUT: launching at exactly
        // r = R_OUT trips the march's escape test on the first iteration, which
        // would make the negative assertions below pass for the wrong reason.
        // With sigma = 0 the start radius is otherwise physically irrelevant.
        const START_R: f64 = 5.0;
        let launch = |chi: f64| Photon {
            pos: DVec3::new(0.05, 0.0, -START_R),
            dir: DVec3::Z,
            weight: 1.0,
            kind: PhotonKind::Ambient,
            bounces: 0,
            chirality: chi,
            spin: chi,
        };

        for chi in [1.0_f64, -1.0] {
            let mut rng = 0xC0FF_EE00_1234_0001u64;
            let mut field = DensityField::new();
            let (term, _cb, _sb, spin, through) = march(
                launch(chi), &pure_e, 0.0, &pure_a, 0.0, 0, &gears_on, &mut rng, &mut field,
            );
            assert_eq!(through, 1, "axial lane photon did not cross (chi {chi})");
            assert!(
                matches!(term, Termination::Escaped { .. }),
                "crossing photon should leave the domain, not be absorbed"
            );
            // SPECIES PRESERVED: spin_transfer is 0 here, so the exiting spin
            // must still be the entering species sign.
            assert_eq!(spin, chi, "through-charge lost its species (chi {chi})");

            // Same photon, lane closed ⇒ absorbed, no crossing.
            let mut rng2 = 0xC0FF_EE00_1234_0001u64;
            let mut field2 = DensityField::new();
            let (term_off, _, _, _, through_off) = march(
                launch(chi), &pure_e, 0.0, &pure_a, 0.0, 0, &gears_off, &mut rng2, &mut field2,
            );
            assert_eq!(through_off, 0, "limb fired with the lane closed");
            assert!(
                matches!(term_off, Termination::Absorbed { .. }),
                "with the lane closed an axial photon must be absorbed"
            );
        }

        // Off-lane (large rho) and off-axis (oblique) photons must NOT cross,
        // per salt.pdf: charge entering near the edges or on an angle is pulled
        // into the equatorial whirlpool instead.
        let mut rng = 0xC0FF_EE00_1234_0002u64;
        let mut field = DensityField::new();
        let wide = Photon {
            pos: DVec3::new(0.8, 0.0, -START_R),
            ..launch(1.0)
        };
        let (_t, _c, _s, _sp, through_wide) = march(
            wide, &pure_e, 0.0, &pure_a, 0.0, 0, &gears_on, &mut rng, &mut field,
        );
        assert_eq!(through_wide, 0, "off-lane photon crossed — lane gate broken");

        let mut rng = 0xC0FF_EE00_1234_0003u64;
        let mut field = DensityField::new();
        let oblique_dir = DVec3::new(0.6, 0.0, 0.8);
        let oblique = Photon {
            pos: -oblique_dir * START_R,
            dir: oblique_dir,
            ..launch(1.0)
        };
        let (_t, _c, _s, _sp, through_obl) = march(
            oblique, &pure_e, 0.0, &pure_a, 0.0, 0, &gears_on, &mut rng, &mut field,
        );
        assert_eq!(
            through_obl, 0,
            "oblique photon crossed — axial-incidence gate broken"
        );
    }

    /// The two-peak profile must reproduce the pre-generalization expression
    /// BIT-FOR-BIT: an ULP shift propagates through the sigma_e calibration and
    /// can flip a collision coin, which would silently break the v2/v3/v4
    /// regression anchors.
    #[test]
    fn baryon_exit_profile_is_bit_identical_to_hardcoded_bimodal() {
        let two_s2 = 2.0 * EXIT_SIGMA_DEG * EXIT_SIGMA_DEG;
        for &lat in &[0.0, 3.5, 7.0, 21.0, 45.0, 58.0, -7.0, -58.0, 89.0] {
            let g = |mu: f64| (-((lat - mu) * (lat - mu)) / two_s2).exp();
            let old = 0.25
                * (g(EXIT_PEAK_LAT_DEG[0])
                    + g(-EXIT_PEAK_LAT_DEG[0])
                    + g(EXIT_PEAK_LAT_DEG[1])
                    + g(-EXIT_PEAK_LAT_DEG[1]));
            let new = exit_lat_profile(ParticleClass::Proton, lat);
            assert_eq!(
                old.to_bits(),
                new.to_bits(),
                "profile drifted at lat {lat}: {old:.20e} vs {new:.20e}"
            );
        }
        // Neutron shares the baryon geometry exactly.
        assert_eq!(
            exit_lat_profile(ParticleClass::Proton, 30.0).to_bits(),
            exit_lat_profile(ParticleClass::Neutron, 30.0).to_bits()
        );
    }

    /// The electron traces a SINGLE exit family (level 9 = axial spin only), so
    /// its profile must have one peak per hemisphere where the baryon has two.
    #[test]
    fn electron_profile_is_single_peak_baryon_is_bimodal() {
        assert_eq!(ParticleClass::Electron.exit_peaks_deg().len(), 1);
        assert_eq!(ParticleClass::Proton.exit_peaks_deg().len(), 2);

        let count_maxima = |class: ParticleClass| {
            let v: Vec<f64> = (0..90)
                .map(|d| exit_lat_profile(class, d as f64))
                .collect();
            (1..v.len() - 1)
                .filter(|&i| v[i] > v[i - 1] && v[i] > v[i + 1])
                .count()
        };
        assert_eq!(
            count_maxima(ParticleClass::Electron),
            1,
            "electron must show exactly one northern-hemisphere peak"
        );
        assert_eq!(
            count_maxima(ParticleClass::Proton),
            2,
            "baryon must keep the bimodal 7/58 structure"
        );
    }

    /// The magnetic channel must be inert at `spin_transfer = 0` (spin stays the
    /// frozen tag) and must actually move at a positive value. This is what the
    /// old `chi_escaped_sum` could never do — it summed the immutable tag.
    #[test]
    fn spin_channel_is_inert_at_zero_transfer_and_live_above_it() {
        let run = |k: f64, class: ParticleClass| -> (f64, f64) {
            let gears = GearConfig {
                enabled: true,
                photon_fraction: PHOTON_FRACTION,
                cancel_redirect: 0.5,
                pole_intake_bias: 0.5,
                class,
                spin_transfer: k,
                through_lane_rho: 0.0,
            };
            let mut seed = 0x51E1_2C0D_4711_0001u64;
            let n = 3000;
            let (emitted0, ambient0) = run_iteration0(n, DEFAULT_SWING_BOOST, &gears, &mut seed);
            let cell = run_self_consistent(
                n, DEFAULT_SWING_BOOST, 3.0, 3.0, 30, 2, false, &gears, &emitted0, &ambient0, seed,
            );
            let e = &cell.emitted_tally;
            (e.spin_escaped_sum, e.escaped)
        };

        for class in [ParticleClass::Proton, ParticleClass::Neutron] {
            let (m0, esc0) = run(0.0, class);
            assert!(esc0 > 0.0, "no charge escaped — cell is degenerate");
            // Frozen spin ⇒ the sum is exactly class_sign × escaped.
            assert_eq!(
                m0, class.chirality_sign() * esc0,
                "spin moved at spin_transfer = 0 for {:?}",
                class
            );

            let (m1, esc1) = run(0.25, class);
            assert!(
                (m1 / esc1.max(1.0)).abs() > 1.0,
                "spin did not stack at spin_transfer = 0.25 for {:?}: M/esc = {}",
                class,
                m1 / esc1.max(1.0)
            );
        }
    }

    /// CM-3: the ELECTRON's collective profile, plus the SWING-BOOST control
    /// the already-exported baryon CSVs need.
    ///
    /// Two jobs in one sweep (class × swing_boost), because they share cells:
    ///
    /// 1. ELECTRON. The recycling engine was baryon-tuned; this is its first
    ///    run on a level-9 particle. The INPUT is only the electron's core
    ///    single-peak disc (±3°, sigma 2.5 — one axial spin traces one exit
    ///    family, PHYSICS_REFERENCE §3b), deliberately WITHOUT the low
    ///    mid-latitude plateau the trace histogram shows out to ~45°. Whether
    ///    the engine REGENERATES that plateau collectively is the validation:
    ///    if it does, the trace's tail is a collective scattering effect and
    ///    the two derivations agree; if it doesn't, the electron needs its own
    ///    emission mechanism rather than a baryon-shaped one.
    ///
    /// 2. SWING CONTROL. `DEFAULT_SWING_BOOST` = 0.05 is the ELECTRON's
    ///    measured settled outer_spin (report_spin_equilibrium: electron s*
    ///    ≈ 0.0504 at Earth mix, matching the 0.055c trace anchor), whereas
    ///    CM-2 measured BARYONS at a TRUE ZERO at Earth mix under
    ///    DirRelSign × the 1/16385 ladder. The committed
    ///    `histogram_{proton,neutron}_recycling.csv` were generated at 0.05,
    ///    i.e. with a swing the baryons should not have. Comparing swing 0.05
    ///    vs 0.0 says whether the proton disc / neutron butterfly survive —
    ///    and therefore whether those CSVs need regenerating before Atom mode
    ///    ever loads them.
    ///
    /// Each cell's escaped-latitude histogram is compared against that class's
    /// shipped TRACE csv, so the collective and single-particle derivations can
    /// be read side by side per band.
    #[test]
    #[ignore = "long-running report; run with --ignored --nocapture"]
    fn report_electron_and_swing_control() {
        use std::time::Instant;
        let started = Instant::now();

        let n = 100_000usize;
        let redirect = 0.5_f64;
        let intake_bias = 0.5_f64;
        let mut seed = 0xE1EC_7207_5D1Au64;

        println!("=== CM-3 Electron Profile + Swing-Boost Control ===");
        println!(
            "N={n}/pop, iters={ITERS}, earth-2/3, cancel_redirect={redirect}, pole_intake_bias={intake_bias}"
        );
        println!(
            "Electron input profile: single peak ±{}° sigma {} (level 9 = axial spin only). Baryon input: bimodal ±{}/±{}° sigma {}.",
            ELECTRON_EXIT_PEAK_LAT_DEG[0],
            ELECTRON_EXIT_SIGMA_DEG,
            EXIT_PEAK_LAT_DEG[0],
            EXIT_PEAK_LAT_DEG[1],
            EXIT_SIGMA_DEG
        );
        println!(
            "swing 0.05 = the electron's measured settled outer_spin; swing 0.00 = the baryons' measured TRUE ZERO at Earth mix. The already-exported baryon CSVs used 0.05.\n"
        );

        // Folds a 180-bin signed-latitude histogram to 91 pole->equator bins
        // (the same folding atom_core::EmissionTable::from_csv_values does) and
        // max-normalizes, so engine output and trace CSV are directly
        // comparable. Index 0 = pole, 90 = equator.
        let fold_norm = |h: &[f64]| -> Vec<f64> {
            let mut bins = vec![0.0f64; 91];
            for (theta, b) in bins.iter_mut().enumerate() {
                let south = h.get(theta).copied().unwrap_or(0.0);
                let north = h.get(179usize.saturating_sub(theta)).copied().unwrap_or(0.0);
                *b = 0.5 * (south + north);
            }
            let max = bins.iter().copied().fold(0.0f64, f64::max).max(1e-12);
            for b in &mut bins {
                *b /= max;
            }
            bins
        };

        let dir = crate::atom_core::config_dir();
        let classes = [
            ("proton", ParticleClass::Proton),
            ("neutron", ParticleClass::Neutron),
            ("electron", ParticleClass::Electron),
        ];
        let swings = [0.05_f64, 0.0];
        // Probe angles measured FROM THE POLE (so 90 = equator).
        let probes = [0usize, 10, 20, 32, 45, 58, 70, 90];

        struct Row {
            cname: &'static str,
            swing: f64,
            escaped: f64,
            band: [f64; 3],
            peak_from_pole: usize,
            folded: Vec<f64>,
        }
        let mut rows: Vec<Row> = Vec::new();

        for (cname, class) in classes {
            for &swing in &swings {
                let gears = GearConfig {
                    enabled: true,
                    photon_fraction: PHOTON_FRACTION,
                    cancel_redirect: redirect,
                    pole_intake_bias: intake_bias,
                    class,
                    spin_transfer: DEFAULT_SPIN_TRANSFER,
                    through_lane_rho: 0.0,
                };
                seed = seed.wrapping_add(0x100);
                let (emitted0, ambient0) = run_iteration0(n, swing, &gears, &mut seed);
                seed = seed.wrapping_add(1);
                let cell = run_self_consistent(
                    n, swing, 3.0, 3.0, 30, ITERS, false, &gears, &emitted0, &ambient0, seed,
                );
                let e = &cell.emitted_tally;
                let raw = if e.escaped_lat_hist.is_empty() {
                    vec![0.0; 180]
                } else {
                    e.escaped_lat_hist.clone()
                };
                let folded = fold_norm(&raw);
                let peak = folded
                    .iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(i, _)| i)
                    .unwrap_or(0);
                let esc = e.escaped.max(1.0);
                rows.push(Row {
                    cname,
                    swing,
                    escaped: e.escaped,
                    band: [
                        e.escaped_band[0] / esc,
                        e.escaped_band[1] / esc,
                        e.escaped_band[2] / esc,
                    ],
                    peak_from_pole: peak,
                    folded,
                });
            }
        }

        println!("--- [1] Exit topology by class and swing ---");
        println!(
            "{:<9} {:>6} {:>9} {:>8} {:>8} {:>8} {:>14}",
            "class", "swing", "escaped", "polEsc", "midEsc", "eqEsc", "peak(from pole)"
        );
        for r in &rows {
            println!(
                "{:<9} {:>6.2} {:>9.0} {:>8.3} {:>8.3} {:>8.3} {:>14}",
                r.cname, r.swing, r.escaped, r.band[0], r.band[1], r.band[2], r.peak_from_pole
            );
        }

        println!("\n--- [2] SWING CONTROL: does killing the baryon swing change the profile? ---");
        for (cname, _) in classes {
            let a = rows.iter().find(|r| r.cname == cname && r.swing == 0.05);
            let b = rows.iter().find(|r| r.cname == cname && r.swing == 0.0);
            if let (Some(a), Some(b)) = (a, b) {
                let l1: f64 = a
                    .folded
                    .iter()
                    .zip(&b.folded)
                    .map(|(x, y)| (x - y).abs())
                    .sum::<f64>()
                    / 91.0;
                let d_eq = (a.band[2] - b.band[2]).abs();
                println!(
                    "  {cname:<9} mean|Δ| per bin = {l1:.4}, ΔeqEsc = {d_eq:.4}, peak {} -> {}  => {}",
                    a.peak_from_pole,
                    b.peak_from_pole,
                    if l1 < 0.05 && d_eq < 0.03 {
                        "INSENSITIVE (exported CSVs stand)"
                    } else {
                        "MATERIAL CHANGE (regenerate before Atom mode uses them)"
                    }
                );
            }
        }

        println!("\n--- [3] Collective profile vs the shipped TRACE csv (normalized, by angle from pole) ---");
        println!("value pairs are engine / trace");
        print!("{:<9} {:>6}", "class", "swing");
        for p in probes {
            print!("{:>14}", format!("{}deg", p));
        }
        println!();
        for r in &rows {
            let trace_path = dir.join(format!("histogram_{}.csv", r.cname));
            let trace_raw = crate::atom_core::load_histogram_csv(&trace_path);
            let trace: Vec<f64> = trace_raw.iter().map(|v| *v as f64).collect();
            let trace_folded = fold_norm(&trace);
            print!("{:<9} {:>6.2}", r.cname, r.swing);
            for p in probes {
                print!(
                    "{:>14}",
                    format!("{:.2}/{:.2}", r.folded[p.min(90)], trace_folded[p.min(90)])
                );
            }
            let l1: f64 = r
                .folded
                .iter()
                .zip(&trace_folded)
                .map(|(x, y)| (x - y).abs())
                .sum::<f64>()
                / 91.0;
            println!("   mean|Δ|={l1:.3}");
        }

        println!("\n--- [4] ELECTRON verdict: did the engine regenerate the trace's plateau? ---");
        println!(
            "The electron INPUT was a bare ±3° core with no plateau. The trace CSV has a low shoulder (~0.36-0.46 of peak) from ~45-70° off the pole and a hard zero inside ~30° of the pole. If the engine output shows that shoulder, the trace tail is collective scattering and the two derivations agree."
        );
        for r in rows.iter().filter(|r| r.cname == "electron") {
            let trace_path = dir.join("histogram_electron.csv");
            let trace_raw = crate::atom_core::load_histogram_csv(&trace_path);
            let trace: Vec<f64> = trace_raw.iter().map(|v| *v as f64).collect();
            let tf = fold_norm(&trace);
            // Shoulder = mean level 45-70 from the pole; polar cone = mean 0-25.
            let mean = |v: &[f64], a: usize, b: usize| -> f64 {
                let s: f64 = v[a..=b].iter().sum();
                s / (b - a + 1) as f64
            };
            println!(
                "  swing {:.2}: engine shoulder(45-70°)={:.3} vs trace {:.3} | engine polar cone(0-25°)={:.3} vs trace {:.3}",
                r.swing,
                mean(&r.folded, 45, 70),
                mean(&tf, 45, 70),
                mean(&r.folded, 0, 25),
                mean(&tf, 0, 25),
            );
        }

        println!(
            "\nTotal report runtime: {:.1}s ({} cells)",
            started.elapsed().as_secs_f64(),
            rows.len()
        );
    }

    /// CM-3: the MAGNETIC-vs-ELECTRICAL split, and whether neutron magnetic
    /// neutrality emerges. Two observables that Mathis keeps strictly apart:
    ///
    /// - ELECTRICAL output = the summed LINEAR motion, a genuine VECTOR sum.
    ///   magmom.pdf: "The neutron is always emitting at a right angle to the
    ///   proton, so if the charge field is defined in terms of proton emission,
    ///   the neutron's charge energy will seem to be zero." Measured here as
    ///   the escaped-direction vector sum resolved into the equatorial plane
    ///   (|xy|) vs the pole axis (|z|), alongside the INCOHERENT magnitude sums
    ///   so "cancelled by symmetry" can be told from "no flux".
    ///   ZERO free parameters — pure geometry of where charge left.
    ///
    /// - MAGNETIC output = the summed SPIN, which magmom.pdf says is
    ///   geometry-BLIND ("if we look at the neutron's emitted photons as a
    ///   matter of spin rather than linear motion, they can add in
    ///   immediately"), and which collisions can change while leaving c alone
    ///   ("photon collisions affect only the spins, not c", neutron.pdf).
    ///   Measured as `spin_escaped_sum` under the meshing rule, swept over
    ///   `spin_transfer` so no conclusion rests on one knob value.
    ///
    /// Predictions under test:
    /// 1. Proton exhaust sums into xy, neutron into z (electrical asymmetry).
    /// 2. Earth 2/3 field: the proton's monochiral χ=+1 exhaust STACKS with the
    ///    photon-majority ambient while the neutron's χ=−1 exhaust CANCELS
    ///    against it, so |M_p| > |M_n|. Compared against the experimental
    ///    μ_p/μ_n = 2.7928/1.9130 = 1.4600 (magmom.pdf derives 1.458) — but see
    ///    the honesty note in the output: agreement in MAGNITUDE would ride on
    ///    `spin_transfer`, so the ordering and the field-dependence are the
    ///    real tests, not the number.
    /// 3. Balanced field: neutron.pdf demands "near-total spin cancellation"
    ///    and a magnetically neutral neutron. This is the one EXPECTED TO FAIL
    ///    with the current architecture, and the report is built to show the
    ///    failure rather than hide it — at 0.5 the ambient mean chirality is 0,
    ///    so the meshing drift vanishes and both classes should keep |M| ≈ 1
    ///    per escaped photon. Genuine neutrality needs AMBIENT charge to enter
    ///    one pole and leave the other RETAINING ITS SPECIES (venus2.pdf
    ///    "through charge"), which the pump cannot express: it absorbs ambient
    ///    and re-emits monochiral charge from the disc. Same architectural
    ///    limit that blocked the polar-exit profile until the through-channel.
    #[test]
    #[ignore = "long-running report; run with --ignored --nocapture"]
    fn report_magnetic_neutrality() {
        use std::time::Instant;
        let started = Instant::now();

        let n = 40_000usize;
        let swing_boost = DEFAULT_SWING_BOOST;
        let redirect = 0.5_f64; // the reconciled regime (commit ebf637e)
        let intake_bias = 0.5_f64;
        let mut seed = 0x3A6E_71C0_BEEF_1046u64;

        println!("=== CM-3 Magnetic Neutrality + Electrical/Magnetic Split ===");
        println!(
            "N={n}/pop, iters={ITERS}, cancel_redirect={redirect}, pole_intake_bias={intake_bias}, swing_boost={swing_boost}"
        );
        println!(
            "ELECTRICAL = vector sum of escaped directions (geometry matters). MAGNETIC = scalar spin sum (geometry-blind, magmom.pdf)."
        );
        println!(
            "spin_transfer=0 is the NULL: spin stays the frozen tag, so M/esc must be exactly the class sign (±1) — verifies the channel is wired and nothing else moved it.\n"
        );

        struct Row {
            cname: &'static str,
            phi_label: &'static str,
            k: f64,
            escaped: f64,
            m: f64,
            m_per_esc: f64,
            spin_esc_band: [f64; 3],
            spin_arr_band: [f64; 3],
            xy_vec: f64,
            z_vec: f64,
            xy_abs: f64,
            z_abs: f64,
        }

        let classes = [
            ("proton", ParticleClass::Proton),
            ("neutron", ParticleClass::Neutron),
        ];
        let phis = [("balanced", 0.5_f64), ("earth-2/3", 2.0 / 3.0)];
        let transfers = [0.0_f64, 0.1, 0.25, 0.5];

        let mut rows: Vec<Row> = Vec::new();
        for (cname, class) in classes {
            for (phi_label, phi) in phis {
                for &k in &transfers {
                    let gears = GearConfig {
                        enabled: true,
                        photon_fraction: phi,
                        cancel_redirect: redirect,
                        pole_intake_bias: intake_bias,
                        class,
                        spin_transfer: k,
                        through_lane_rho: 0.0,
                    };
                    seed = seed.wrapping_add(0x100);
                    let (emitted0, ambient0) = run_iteration0(n, swing_boost, &gears, &mut seed);
                    seed = seed.wrapping_add(1);
                    let cell = run_self_consistent(
                        n, swing_boost, 3.0, 3.0, 30, ITERS, false, &gears, &emitted0, &ambient0,
                        seed,
                    );
                    let e = &cell.emitted_tally;
                    let a = &cell.ambient_tally;
                    let esc = e.escaped.max(1.0);
                    rows.push(Row {
                        cname,
                        phi_label,
                        k,
                        escaped: e.escaped,
                        m: e.spin_escaped_sum,
                        m_per_esc: e.spin_escaped_sum / esc,
                        spin_esc_band: e.spin_escaped_band,
                        spin_arr_band: a.spin_arrived_band,
                        xy_vec: (e.dir_escaped_sum.x * e.dir_escaped_sum.x
                            + e.dir_escaped_sum.y * e.dir_escaped_sum.y)
                            .sqrt(),
                        z_vec: e.dir_escaped_sum.z,
                        xy_abs: e.dir_escaped_xy_abs,
                        z_abs: e.dir_escaped_z_abs,
                    });
                }
            }
        }

        println!("--- [1] MAGNETIC: spin sum over escaped emitted charge ---");
        println!(
            "{:<9} {:<10} {:>5} {:>9} {:>11} {:>9}  {:>26} {:>26}",
            "class", "field", "k", "escaped", "M", "M/esc", "spin_esc [pol/mid/eq]", "spin_arr [pol/mid/eq]"
        );
        for r in &rows {
            println!(
                "{:<9} {:<10} {:>5.2} {:>9.0} {:>11.1} {:>9.4}  {:>8.1}{:>9.1}{:>9.1} {:>8.1}{:>9.1}{:>9.1}",
                r.cname,
                r.phi_label,
                r.k,
                r.escaped,
                r.m,
                r.m_per_esc,
                r.spin_esc_band[0],
                r.spin_esc_band[1],
                r.spin_esc_band[2],
                r.spin_arr_band[0],
                r.spin_arr_band[1],
                r.spin_arr_band[2],
            );
        }

        println!("\n--- [2] ELECTRICAL: vector vs incoherent sums of escaped directions ---");
        println!(
            "A DISC exhaust cancels in the vector sum by azimuthal symmetry but has large |xy| incoherent; an AXIAL exhaust survives in z. The magmom claim is about which PLANE the flux occupies, so the incoherent columns are the diagnostic ones."
        );
        println!(
            "{:<9} {:<10} {:>5} {:>11} {:>11} {:>13} {:>13} {:>11}",
            "class", "field", "k", "|xy_vec|", "z_vec", "xy_abs/esc", "z_abs/esc", "z/xy"
        );
        for r in &rows {
            let esc = r.escaped.max(1.0);
            let xa = r.xy_abs / esc;
            let za = r.z_abs / esc;
            println!(
                "{:<9} {:<10} {:>5.2} {:>11.1} {:>11.1} {:>13.4} {:>13.4} {:>11.3}",
                r.cname,
                r.phi_label,
                r.k,
                r.xy_vec,
                r.z_vec,
                xa,
                za,
                za / xa.max(1e-12),
            );
        }

        println!("\n--- [3] VERDICT: proton/neutron magnetic ratio vs the 1.46 anchor ---");
        println!(
            "Experimental |mu_p/mu_n| = 2.7928/1.9130 = 1.4600; magmom.pdf derives 1.458 from an angle differential plus the proton's faster spin. Neither of those routes is what this model uses, so treat a near-1.46 value as suggestive, NOT as a derivation."
        );
        for (phi_label, _) in phis {
            for &k in &transfers {
                let find = |c: &str| {
                    rows.iter()
                        .find(|r| r.cname == c && r.phi_label == phi_label && r.k == k)
                        .map(|r| r.m_per_esc)
                        .unwrap_or(0.0)
                };
                let mp = find("proton");
                let mn = find("neutron");
                let ratio = if mn.abs() > 1e-12 {
                    (mp / mn).abs()
                } else {
                    f64::INFINITY
                };
                println!(
                    "  {phi_label:<10} k={k:<5.2} M_p/esc={mp:>8.4}  M_n/esc={mn:>8.4}  |ratio|={ratio:>8.3}{}",
                    if (ratio - 1.46).abs() < 0.15 && k > 0.0 {
                        "   <-- within 10% of 1.46"
                    } else {
                        ""
                    }
                );
            }
        }

        println!("\n--- [4] The balanced-field test neutron.pdf demands ---");
        println!(
            "neutron.pdf: balanced field => 'near-total spin cancellation ... the neutron is then neutral regarding magnetism', while charge is still recycled. PASS would be |M_n/esc| -> ~0 at photon_fraction 0.5 with escaped still large."
        );
        for &k in &transfers {
            let find = |c: &str| {
                rows.iter()
                    .find(|r| r.cname == c && r.phi_label == "balanced" && r.k == k)
                    .map(|r| (r.m_per_esc, r.escaped))
                    .unwrap_or((0.0, 0.0))
            };
            let (mn, escn) = find("neutron");
            let (mp, _) = find("proton");
            println!(
                "  k={k:<5.2} balanced: M_n/esc={mn:>8.4} (escaped {escn:>7.0})  M_p/esc={mp:>8.4}  -> {}",
                if mn.abs() < 0.15 {
                    "NEUTRAL (prediction met)"
                } else {
                    "NOT neutral — spin survives; see the through-charge note in this test's doc comment"
                }
            );
        }

        println!(
            "\nTotal report runtime: {:.1}s ({} cells)",
            started.elapsed().as_secs_f64(),
            rows.len()
        );
    }

    /// Exports atom-mode emission profiles for proton and neutron DERIVED FROM
    /// THE COLLECTIVE RECYCLING ENGINE (not the single-particle trace). Writes
    /// `godot/config/histogram_{proton,neutron}_recycling.csv` in the same
    /// `degree,percentage` format Atom mode loads (M5). These are NEW files —
    /// the trace-derived `histogram_*.csv` the atom_scenarios tests calibrate
    /// against are left untouched; the next session decides whether to swap
    /// Atom mode over.
    ///
    /// Config = the reconciled regime: symmetric gear at cancel_redirect 0.5,
    /// Earth 2/3 field, moderate species-split intake — the proton stays a
    /// clean disc, the neutron is disc-suppressed and mid/high-latitude
    /// (matches the ~57° trace butterfly while being chirality-emergent).
    /// FORCED mode (fixed large N) for smooth per-degree statistics.
    #[test]
    #[ignore = "generates atom-mode profile CSVs; run with --ignored --nocapture"]
    fn export_recycling_profiles() {
        use std::time::Instant;
        let started = Instant::now();

        let n = 200_000usize;
        let swing_boost = DEFAULT_SWING_BOOST;
        let redirect = 0.5_f64;
        let intake_bias = 0.5_f64;
        let mut seed = 0xE111_7ED0_F00D_2027u64;

        let dir = crate::atom_core::config_dir();
        println!("=== CM-3 Recycling Profile Export (atom-mode emission histograms) ===");
        println!(
            "Forced mode, N={n}, earth-2/3, cancel_redirect={redirect}, pole_intake_bias={intake_bias}, iters={ITERS}. Writing to {}",
            dir.display()
        );

        let classes = [
            ("proton", ParticleClass::Proton),
            ("neutron", ParticleClass::Neutron),
            // Added after report_electron_and_swing_control showed the engine
            // REGENERATES the electron trace's mid-latitude shoulder (0.383 vs
            // trace 0.359) from a bare ±3° core input — the two derivations
            // agree for the electron, so all three profiles can come from one
            // engine instead of the electron staying trace-only.
            ("electron", ParticleClass::Electron),
        ];

        for (cname, class) in classes {
            let gears = GearConfig {
                enabled: true,
                photon_fraction: PHOTON_FRACTION,
                cancel_redirect: redirect,
                pole_intake_bias: intake_bias,
                class,
                spin_transfer: DEFAULT_SPIN_TRANSFER,
                // Through-charge limb OFF: this report's committed numbers were
                // measured without it, so keep them reproducible.
                through_lane_rho: 0.0,
            };
            seed = seed.wrapping_add(0x100);
            let (emitted0, ambient0) = run_iteration0(n, swing_boost, &gears, &mut seed);
            seed = seed.wrapping_add(1);
            let cell = run_self_consistent(
                n, swing_boost, 3.0, 3.0, 30, ITERS, false, &gears, &emitted0, &ambient0, seed,
            );
            let e = &cell.emitted_tally;
            let raw = if e.escaped_lat_hist.is_empty() {
                vec![0.0; 180]
            } else {
                e.escaped_lat_hist.clone()
            };
            // Light 3-bin box smooth (per-degree counts are ~1k, ~3% noise).
            let mut hist = vec![0.0f64; 180];
            for i in 0..180usize {
                let lo = i.saturating_sub(1);
                let hi = (i + 1).min(179);
                hist[i] = (raw[lo] + raw[i] + raw[hi]) / (hi - lo + 1) as f64;
            }
            let max = hist.iter().copied().fold(0.0f64, f64::max).max(1e-12);

            let esc = e.escaped.max(1.0);
            println!(
                "{cname}: escaped {} | eqEsc {:.3} midEsc {:.3} polEsc {:.3} | peak at {}° ",
                e.escaped as usize,
                e.escaped_band[Band::Equatorial as usize] / esc,
                e.escaped_band[Band::Mid as usize] / esc,
                e.escaped_band[Band::Polar as usize] / esc,
                hist.iter()
                    .enumerate()
                    .max_by(|a, b| a.1.partial_cmp(b.1).unwrap())
                    .map(|(i, _)| i as i32 - 90)
                    .unwrap_or(0),
            );

            let mut csv = String::from("degree,percentage\n");
            for d in -90i32..90 {
                let v = hist[(d + 90) as usize] / max;
                csv.push_str(&format!("{d},{v:.4}\n"));
            }
            let path = dir.join(format!("histogram_{cname}_recycling.csv"));
            std::fs::write(&path, &csv)
                .unwrap_or_else(|err| panic!("write {}: {err}", path.display()));

            // Round-trip: confirm Atom mode's loader + table accept the file.
            let reload = crate::atom_core::load_histogram_csv(&path);
            let table = crate::atom_core::EmissionTable::from_csv_values(&reload);
            assert_eq!(table.bins.len(), 91, "profile must fold to 91 pole→equator bins");
            println!("  wrote + verified {} ({} rows)", path.display(), reload.len());
        }
        println!("done in {:.1}s", started.elapsed().as_secs_f64());
    }
}
