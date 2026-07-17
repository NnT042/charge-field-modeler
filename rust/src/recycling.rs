//! CM-3: the disc-to-pole recycling loop (diagnostic) — v2.
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

use glam::DVec3;
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
/// Photons per population per pass (report default; unit tests use less).
pub const DEFAULT_N: usize = 200_000;
/// Tangential exit-direction boost — the settled proton-scale outer_spin
/// rate (see project memory: presets settled around ±0.065, "outer velocity:
/// non-relativistic only" landed near 0.055c). Used only as this module's
/// own default; not read from `calibration::CalibrationParticle`.
pub const DEFAULT_SWING_BOOST: f64 = 0.05;
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

/// Measured proton exit-structure peaks (project memory: 7°/58° dual peaks).
const EXIT_PEAK_LAT_DEG: [f64; 2] = [7.0, 58.0];
const EXIT_SIGMA_DEG: f64 = 5.0;

/// Unnormalized bimodal exit-latitude mixture pdf: equal-weight Gaussians at
/// ±7° and ±58°, sigma 5°. Normalization is irrelevant — sigma_e is
/// calibrated against this same profile (see `tau_in_integral`), so only the
/// profile's SHAPE enters the physics.
fn bimodal_lat_profile(lat_deg: f64) -> f64 {
    let two_s2 = 2.0 * EXIT_SIGMA_DEG * EXIT_SIGMA_DEG;
    let g = |mu: f64| (-((lat_deg - mu) * (lat_deg - mu)) / two_s2).exp();
    0.25 * (g(EXIT_PEAK_LAT_DEG[0])
        + g(-EXIT_PEAK_LAT_DEG[0])
        + g(EXIT_PEAK_LAT_DEG[1])
        + g(-EXIT_PEAK_LAT_DEG[1]))
}

/// Analytic seed density for a population (v2 fix #2): emitted =
/// profile(θ)/r² (free-streaming bimodal source), ambient = 1 (uniform
/// far-field value in these units).
fn analytic_density(kind: PhotonKind, r: f64, theta: f64) -> f64 {
    match kind {
        PhotonKind::Emitted => {
            let lat_deg = 90.0 - theta.to_degrees();
            bimodal_lat_profile(lat_deg) / (r * r).max(1.0)
        }
        PhotonKind::Ambient => 1.0,
    }
}

/// ∫ analytic n_e dr along the b=0 equatorial line from TAU_IN_R_MIN to
/// TAU_IN_R_MAX: profile(0°)·(1/r_min − 1/r_max). sigma_e = TAU / this.
fn tau_in_integral() -> f64 {
    bimodal_lat_profile(0.0) * (1.0 / TAU_IN_R_MIN - 1.0 / TAU_IN_R_MAX)
}

/// ∫ analytic n_a dr for an equatorial emitted photon r = 1 → R_OUT:
/// 1·(R_OUT − R_IN). sigma_a = TAU / this.
fn tau_out_integral() -> f64 {
    R_OUT - R_IN
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
}

impl DensityField {
    fn new() -> Self {
        Self {
            track_weight: vec![0.0; GRID_CELLS],
            dir_accum: vec![DVec3::ZERO; GRID_CELLS],
            radial_flux: vec![0.0; GRID_CELLS],
        }
    }

    fn deposit(&mut self, pos: DVec3, dir: DVec3, weight: f64, ds: f64) {
        let r = pos.length().max(1e-9);
        let theta = colatitude_of(pos, r);
        let idx = grid_index(r, theta);
        let r_hat = pos / r;
        let w = weight * ds;
        self.track_weight[idx] += w;
        self.dir_accum[idx] += dir * w;
        self.radial_flux[idx] += w * dir.dot(r_hat);
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
    tally: DensityField,
    /// Global scale mapping tallied track-density into analytic units;
    /// 0.0 = no usable tally → pure analytic (correction ≡ 1).
    norm: f64,
}

impl CollisionField {
    fn pure_analytic(kind: PhotonKind) -> Self {
        Self {
            kind,
            tally: DensityField::new(),
            norm: 0.0,
        }
    }

    /// Wraps a tally. The normalization region differs by population:
    /// emitted launches are tracked completely, so the whole domain
    /// participates; ambient tallies are only shape-complete for r ≤ B_MAX
    /// (see module doc, v2 fix #2), so both the norm and the correction are
    /// restricted to that near zone (correction ≡ 1 outside it).
    fn from_tally(kind: PhotonKind, tally: DensityField) -> Self {
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
                analytic_total += analytic_density(kind, r_c, th_c) * cell_volumes()[idx];
            }
        }
        let norm = if measured_total > 1e-12 {
            analytic_total / measured_total
        } else {
            0.0
        };
        Self { kind, tally, norm }
    }

    fn effective_density_at(&self, r: f64, theta: f64) -> f64 {
        let a = analytic_density(self.kind, r, theta);
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
fn effective_density_snapshot(kind: PhotonKind, tally: &DensityField) -> Vec<f64> {
    let field = CollisionField::from_tally(kind, tally.clone());
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

#[derive(Clone, Copy, Debug)]
struct Photon {
    pos: DVec3,
    dir: DVec3,
    weight: f64,
    kind: PhotonKind,
    bounces: u32,
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

/// Draws a peak MAGNITUDE (0..89, unsigned) from the equal-weight two-peak
/// mixture. Sign/hemisphere is applied by the caller — see module doc
/// deviation #1 for why sign is shared rather than redrawn per latitude.
fn sample_exit_magnitude_deg(rng: &mut u64) -> f64 {
    let peak = if xorshift64(rng) < 0.5 {
        EXIT_PEAK_LAT_DEG[0]
    } else {
        EXIT_PEAK_LAT_DEG[1]
    };
    rand_gaussian(rng, peak, EXIT_SIGMA_DEG).clamp(0.0, 89.0)
}

fn spawn_emitted(rng: &mut u64, swing_boost: f64) -> Photon {
    let sign = if xorshift64(rng) < 0.5 { 1.0 } else { -1.0 };
    let lat0 = sign * sample_exit_magnitude_deg(rng);
    let az0 = xorshift64(rng) * TAU;
    let r_hat0 = latlon_to_unit(lat0, az0);
    let pos = r_hat0 * 1.001;

    // NEW draw for the exit direction's own latitude (same hemisphere sign,
    // independent peak/magnitude draw — module doc deviation #1).
    let lat1 = sign * sample_exit_magnitude_deg(rng);
    let dir_radial = latlon_to_unit(lat1, az0);
    let phi_hat = phi_hat_at(dir_radial);
    let dir = (dir_radial + phi_hat * swing_boost).normalize();

    Photon {
        pos,
        dir,
        weight: 1.0,
        kind: PhotonKind::Emitted,
        bounces: 0,
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
fn spawn_ambient(rng: &mut u64) -> (Photon, Option<f64>) {
    let d = rand_unit_vec(rng);
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

/// Ray-marches one photon through the OTHER population's regularized
/// collision field (`other`, with coupling `sigma`), depositing its own
/// track-length into `own_field`, until absorbed (r<R_IN), escaped
/// (r>R_OUT), or lost (safety path-length cap). Returns the termination and
/// the photon's final collision (bounce) count.
fn march(
    mut photon: Photon,
    other: &CollisionField,
    sigma: f64,
    n_ricochet: u32,
    rng: &mut u64,
    own_field: &mut DensityField,
) -> (Termination, u32) {
    let mut path_length = 0.0;
    let term = loop {
        let r = photon.pos.length();
        if r < R_IN {
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

        if sigma > 0.0 && photon.bounces < n_ricochet {
            let n_density = other.effective_density_at(r, theta);
            let p_collide = 1.0 - (-sigma * n_density * ds).exp();
            if xorshift64(rng) < p_collide {
                let (mean_dir, coherence) = other.tally.mean_dir_and_coherence_at(r, theta);
                let mut partner_dir = mean_dir * coherence + rand_unit_vec(rng);
                if partner_dir.length_squared() < 1e-12 {
                    partner_dir = rand_unit_vec(rng);
                }
                partner_dir = partner_dir.normalize();

                let rel = photon.dir - partner_dir;
                if rel.length_squared() > 1e-12 {
                    let rel_hat = rel.normalize();
                    let mut new_dir = None;
                    for _ in 0..4 {
                        let nhat = sample_contact_normal(rel_hat, rng);
                        let candidate = photon.dir - photon.dir.dot(nhat) * nhat
                            + partner_dir.dot(nhat) * nhat;
                        if candidate.length_squared() > 1e-8 {
                            new_dir = Some(candidate.normalize());
                            break;
                        }
                    }
                    if let Some(dir) = new_dir {
                        photon.dir = dir;
                    }
                }
                photon.bounces += 1;
            }
        }

        photon.pos += photon.dir * ds;
        path_length += ds;
    };
    (term, photon.bounces)
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
    lost: f64,
    lz_arrived: f64,
    lz_escaped: f64,
    /// Σ collision (bounce) events over ALL photons — knob verification.
    collisions_sum: f64,
    /// Σ bounces over ABSORBED photons — "are arrivals the multiply-
    /// scattered ones?"
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
    fn record_emitted(&mut self, term: &Termination, bounces: u32) {
        self.total += 1.0;
        self.collisions_sum += bounces as f64;
        match term {
            Termination::Absorbed { pos, dir } => {
                self.arrived += 1.0;
                self.arrived_bounce_sum += bounces as f64;
                let band = lat_band(arrival_latitude(*pos));
                self.arrived_band[band as usize] += 1.0;
                self.lz_arrived += pos.cross(*dir).z;
            }
            Termination::Escaped { pos, dir } => {
                self.escaped += 1.0;
                self.lz_escaped += pos.cross(*dir).z;
            }
            Termination::Lost => self.lost += 1.0,
        }
    }

    fn record_ambient(&mut self, term: &Termination, bounces: u32, aim_lat: Option<f64>) {
        self.total += 1.0;
        self.collisions_sum += bounces as f64;
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
                self.lz_escaped += pos.cross(*dir).z;
            }
            Termination::Lost => self.lost += 1.0,
        }
    }
}

fn run_pass(
    kind: PhotonKind,
    n: usize,
    swing_boost: f64,
    other: &CollisionField,
    sigma: f64,
    n_ricochet: u32,
    rng: &mut u64,
) -> (DensityField, PassTally, f64) {
    let mut field = DensityField::new();
    let mut tally = PassTally::default();
    let mut lz_source_total = 0.0;
    for _ in 0..n {
        match kind {
            PhotonKind::Emitted => {
                let photon = spawn_emitted(rng, swing_boost);
                lz_source_total += photon.pos.cross(photon.dir).z * photon.weight;
                let (term, bounces) = march(photon, other, sigma, n_ricochet, rng, &mut field);
                tally.record_emitted(&term, bounces);
            }
            PhotonKind::Ambient => {
                let (photon, aim_lat) = spawn_ambient(rng);
                let (term, bounces) = march(photon, other, sigma, n_ricochet, rng, &mut field);
                tally.record_ambient(&term, bounces, aim_lat);
            }
        }
    }
    (field, tally, lz_source_total)
}

/// Free-streaming (sigma=0) seed passes: their tallies seed the correction
/// and mean-direction fields for iteration 1 of every sweep cell. The
/// launch-shell density spike this produces near r=1.001 is harmless in v2:
/// the correction cap bounds its effect on collisions at 4× analytic.
fn run_iteration0(n: usize, swing_boost: f64, rng: &mut u64) -> (DensityField, DensityField) {
    let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient);
    let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted);
    let (emitted0, _, _) = run_pass(
        PhotonKind::Emitted,
        n,
        swing_boost,
        &pure_ambient,
        0.0,
        0,
        rng,
    );
    let (ambient0, _, _) = run_pass(
        PhotonKind::Ambient,
        n,
        swing_boost,
        &pure_emitted,
        0.0,
        0,
        rng,
    );
    (emitted0, ambient0)
}

struct ReportCell {
    tau_target: f64,
    n_ricochet: u32,
    sigma_e: f64,
    sigma_a: f64,
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
    n_ricochet: u32,
    iters: usize,
    emitted0: &DensityField,
    ambient0: &DensityField,
    seed: u64,
) -> ReportCell {
    // v2 fix #3: symmetric TAU loading from the ANALYTIC seeds, both knobs
    // derived from the one swept TAU.
    let (sigma_e, sigma_a) = if tau_target > 0.0 {
        (tau_target / tau_in_integral(), tau_target / tau_out_integral())
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

    for _ in 0..iters {
        let field_emitted = CollisionField::from_tally(PhotonKind::Emitted, tally_emitted.clone());
        let field_ambient = CollisionField::from_tally(PhotonKind::Ambient, tally_ambient.clone());

        // Emitted march through the AMBIENT field with sigma_a; ambient
        // march through the EMITTED field with sigma_e.
        let (measured_emitted, e_tally, lz_src) = run_pass(
            PhotonKind::Emitted,
            n,
            swing_boost,
            &field_ambient,
            sigma_a,
            n_ricochet,
            &mut rng,
        );
        let (measured_ambient, a_tally, _) = run_pass(
            PhotonKind::Ambient,
            n,
            swing_boost,
            &field_emitted,
            sigma_e,
            n_ricochet,
            &mut rng,
        );

        let new_emitted = tally_emitted.relax(&measured_emitted, 0.5);
        let new_ambient = tally_ambient.relax(&measured_ambient, 0.5);

        let old_e_snap = effective_density_snapshot(PhotonKind::Emitted, &tally_emitted);
        let new_e_snap = effective_density_snapshot(PhotonKind::Emitted, &new_emitted);
        let old_a_snap = effective_density_snapshot(PhotonKind::Ambient, &tally_ambient);
        let new_a_snap = effective_density_snapshot(PhotonKind::Ambient, &new_ambient);
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
        n_ricochet,
        sigma_e,
        sigma_a,
        convergence,
        emitted_tally,
        ambient_tally,
        lz_source_total,
        final_emitted_field,
        final_ambient_field,
    }
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
        let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted);
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient);
        let n = 5000;

        let mut checked = 0;
        let mut hits = 0;
        let mut mismatches = 0;
        for _ in 0..n {
            let (photon, aim_lat) = spawn_ambient(&mut rng);
            let mut field = DensityField::new();
            let (term, bounces) = march(photon, &pure_emitted, 0.0, 0, &mut rng, &mut field);
            assert_eq!(bounces, 0, "collision recorded with sigma=0");
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
            let photon = spawn_emitted(&mut rng, DEFAULT_SWING_BOOST);
            let mut field = DensityField::new();
            let (term, bounces) = march(photon, &pure_ambient, 0.0, 0, &mut rng, &mut field);
            assert_eq!(bounces, 0, "collision recorded with sigma=0");
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
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient);
        let (field, _tally, _lz) = run_pass(
            PhotonKind::Emitted,
            300_000,
            DEFAULT_SWING_BOOST,
            &pure_ambient,
            0.0,
            0,
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
        let pure_emitted = CollisionField::pure_analytic(PhotonKind::Emitted);
        let pure_ambient = CollisionField::pure_analytic(PhotonKind::Ambient);
        let (_field_e, emitted_tally, _lz_src) =
            run_pass(PhotonKind::Emitted, n, 0.0, &pure_ambient, 0.0, 0, &mut rng);
        let (_field_a, ambient_tally, _) =
            run_pass(PhotonKind::Ambient, n, 0.0, &pure_emitted, 0.0, 0, &mut rng);

        let delivered = ambient_tally.lz_arrived + emitted_tally.lz_arrived;
        let carried_off = ambient_tally.lz_escaped + emitted_tally.lz_escaped;
        let net = delivered - carried_off;
        let net_per_photon = net / n as f64;
        assert!(
            net_per_photon.abs() < 0.2,
            "expected ~0 net L_z/photon at swing_boost=0, got {net_per_photon:.5}"
        );
    }

    /// The CM-3 v2 report: sweeps TAU (symmetric optical-depth target, v2
    /// fix #3) x N_RICOCHET (max bounces), self-consistently settles both
    /// populations' regularized density fields for ITERS passes, and prints
    /// the screening, redirection, recycling, circulation, and
    /// spin-sustenance evidence for the disc-to-pole loop claim — plus the
    /// v2 knob verification (measured vs target collision counts, printed
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

        println!("=== CM-3 Recycling Loop Report (v2: importance-sampled ambient, analytic regularized fields, symmetric TAU) ===");
        println!(
            "N per population = {n}, ITERS = {ITERS}, grid = {N_R}x{N_THETA} (r,theta), R_OUT = {R_OUT}, swing_boost = {swing_boost}, B_MAX = {B_MAX}"
        );
        println!(
            "AMBIENT ENSEMBLE IS CONDITIONAL: all ambient photons are straight-line launches with impact parameter b <= {B_MAX} (area-weighted = correct relative flux within the near zone). Absolute capture fractions are conditional on this ensemble."
        );
        println!(
            "TAU calibration (analytic): TAU_in integral (b=0 equatorial, r={TAU_IN_R_MIN}->{TAU_IN_R_MAX}, profile(0)={:.6}) = {:.6}; TAU_out integral (equatorial emitted, r=1->{R_OUT}, n_a=1) = {:.1}",
            bimodal_lat_profile(0.0),
            tau_in_integral(),
            tau_out_integral()
        );

        let (emitted0, ambient0) = run_iteration0(n, swing_boost, &mut seed);

        let tau_values = [0.0, 1.0, 3.0, 10.0];
        let ricochet_values = [3u32, 10, 30];

        let mut cells: Vec<ReportCell> = Vec::new();
        for &tau in &tau_values {
            if tau == 0.0 {
                println!("TAU=0.0 -> sigma_e = 0, sigma_a = 0 (control, ricochet sweep skipped)");
                seed = seed.wrapping_add(1);
                cells.push(run_self_consistent(
                    n, swing_boost, tau, 0, ITERS, &emitted0, &ambient0, seed,
                ));
            } else {
                println!(
                    "TAU={tau:.1} -> sigma_e (ambient-through-emitted) = {:.4}, sigma_a (emitted-through-ambient) = {:.6}",
                    tau / tau_in_integral(),
                    tau / tau_out_integral()
                );
                for &nr in &ricochet_values {
                    seed = seed.wrapping_add(1);
                    cells.push(run_self_consistent(
                        n, swing_boost, tau, nr, ITERS, &emitted0, &ambient0, seed,
                    ));
                }
            }
        }

        println!();
        println!("KNOB VERIFICATION + main sweep table (per cell; cAmb/cEm = MEASURED mean collisions per ambient/emitted photon vs target TAU; bArr/bReent = mean ricochets among ambient arrivals / emitted re-entries):");
        println!(
            "{:>5} {:>4} | {:>7} {:>7} | {:>8} {:>8} {:>8} {:>8} | {:>9} {:>8} {:>8} {:>8} {:>8} | {:>8} {:>8} | {:>6} {:>7} | {:>11} {:>9} {:>5}",
            "TAU", "Nric",
            "cAmb", "cEm",
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

            let c_amb = a.collisions_sum / a.total;
            let c_em = e.collisions_sum / e.total;

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

            let em_reentry = e.arrived / e.total;
            let em_escape = e.escaped / e.total;

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
            let nric_label = if cell.tau_target == 0.0 {
                "n/a".to_string()
            } else {
                cell.n_ricochet.to_string()
            };

            println!(
                "{:>5.1} {:>4} | {:>7.3} {:>7.3} | {:>8.4} {:>8.4} {:>8.4} {:>8.4} | {:>9.4} {:>8.4} {:>8.4} {:>8.4} {:>8.4} | {:>8.4} {:>8.4} | {:>6.2} {:>7.2} | {:>11.6} {:>9.4} {:>5}",
                cell.tau_target, nric_label,
                c_amb, c_em,
                amb_arr_frac, polar_frac, mid_frac, eq_frac,
                redirect, miss_capture, miss_pol, miss_mid, miss_eq,
                em_reentry, em_escape,
                b_arr, b_reent,
                net_per_photon, fb_ratio, sign_label
            );

            // v2 fix #4: the knob must be PROVEN to work. The ricochet cap
            // truncates the achievable count, so compare against
            // min(TAU, N_RICOCHET); the ambient ensemble mean is further
            // diluted below the b=0 equatorial calibration ray by latitude
            // profile + impact-parameter spread, so use a generous O(1)
            // band and shout only when the knob is broken by >5x.
            if cell.tau_target > 0.0 {
                let target = cell.tau_target.min(cell.n_ricochet as f64);
                for (label, measured) in [("ambient", c_amb), ("emitted", c_em)] {
                    let ratio = measured / target;
                    if !(0.2..=5.0).contains(&ratio) {
                        knob_warnings.push(format!(
                            "WARNING: TAU={:.1} Nric={} {} measured mean collisions {:.3} vs effective target {:.1} (ratio {:.3}) — KNOB IS OFF, treat this cell's physics columns with suspicion",
                            cell.tau_target, cell.n_ricochet, label, measured, target, ratio
                        ));
                    }
                }
            }
        }
        if knob_warnings.is_empty() {
            println!("\nKNOB CHECK: all sweep cells' measured collision counts are within [0.2, 5.0]x of min(TAU, N_RICOCHET) — the optical-depth knob is working.");
        } else {
            println!("\nKNOB CHECK FAILURES ({}):", knob_warnings.len());
            for w in &knob_warnings {
                println!("{w}");
            }
        }

        let control = cells
            .iter()
            .find(|c| c.tau_target == 0.0)
            .expect("TAU=0 control cell present");
        let flagship = cells
            .iter()
            .find(|c| c.tau_target == 3.0 && c.n_ricochet == 10)
            .expect("TAU=3, N_RICOCHET=10 flagship cell present");

        for (label, cell) in [
            ("CONTROL (TAU=0)", control),
            ("FLAGSHIP (TAU=3, N_RICOCHET=10)", flagship),
        ] {
            println!("\n--- {label}: per-iteration convergence (L2 relative change of regularized density fields, emitted / ambient) ---");
            for (k, (ce, ca)) in cell.convergence.iter().enumerate() {
                println!("  iter {}: emitted {:.6}  ambient {:.6}", k + 1, ce, ca);
            }

            println!("--- {label}: circulation matrix (rows = radial bin edge r, every 3rd bin; cols = theta bin edge in deg, every 2nd bin) ---");
            println!("    value = mean radial-flow component of combined emitted+ambient flux in that cell, bounded [-1,1]; + = net outward, - = net inward, 0 = no net radial bias.");
            print!("{:>10}", "r\\theta");
            for j in (0..N_THETA).step_by(2) {
                print!(" {:>7.0}", theta_edge(j).to_degrees());
            }
            println!();
            for i in (0..N_R).step_by(3) {
                print!("{:>10.3}", r_edge(i));
                for j in (0..N_THETA).step_by(2) {
                    let idx = i * N_THETA + j;
                    let combined_w = cell.final_emitted_field.track_weight[idx]
                        + cell.final_ambient_field.track_weight[idx];
                    let combined_flux = cell.final_emitted_field.radial_flux[idx]
                        + cell.final_ambient_field.radial_flux[idx];
                    let val = if combined_w > 1e-9 {
                        combined_flux / combined_w
                    } else {
                        0.0
                    };
                    print!(" {val:>7.3}");
                }
                println!();
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
}
