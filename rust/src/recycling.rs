//! CM-3: the disc-to-pole recycling loop (diagnostic).
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
//!     v2 refinement NOT modeled here (this is a single-species, v1 pass).
//!   - Emission profile: this project's OWN measured proton exit structure —
//!     bimodal dual peaks at 7° and 58° latitude (M2-era trace measurement,
//!     see project memory `reference_bimodal_exit_structure`).
//!
//! Units: natural (particle surface radius = 1, c = 1), matching project
//! convention (all internal math in natural units; SI lives only in
//! `units.rs`, not touched by this module).
//!
//! ## Design overview
//!
//! Two photon populations live in a spherical shell `r ∈ [R_IN, R_OUT]`,
//! tallied on an azimuthally-symmetric `(r, θ)` grid (`N_R` log-spaced radial
//! bins × `N_THETA` uniform polar-angle bins, θ = colatitude from the +Z
//! pole). EMITTED photons start at the surface with the bimodal exit-latitude
//! profile; AMBIENT photons start on the outer sphere heading inward
//! (cosine-weighted, isotropic external field). Each is ray-marched through
//! the *other* population's densit­y field from the previous self-consistent
//! iteration, with a collision probability driven by a swept `sigma_eff`, up
//! to `n_ricochet` scattering events per photon (`FIELD ITERATION` below).
//!
//! ## Deliberate deviations / disambiguations from the literal spec
//!
//! 1. **Same-hemisphere emission direction.** The spec draws the emission
//!    direction's latitude from "a NEW draw from the same [bimodal] mixture"
//!    independently of the start latitude. Taken fully independently, the
//!    mixture's two peaks (7°/58°) are far enough apart (58° − (−58°) = 116°
//!    > 90°) that ~12.5% of emitted photons would have direction and position
//!    latitudes on OPPOSITE hemispheres — making the "radially outward" ray
//!    actually point back into the sphere from the get-go (dot(pos,dir) < 0),
//!    which would fail the "zero emitted re-entries at TAU=0" unit test for a
//!    reason that has nothing to do with the recycling physics being tested.
//!    Fix: the hemisphere sign is drawn ONCE per photon and shared by both the
//!    position-latitude draw and the direction-latitude draw; only the
//!    peak-magnitude (7° vs 58°) is redrawn independently. Max separation is
//!    then 58°−7°=51° (or up to ~89° at extreme clamped tails), always < 90°,
//!    so the exit ray is always genuinely outward. This preserves the
//!    bimodal-profile intent (either peak, independently, per draw) while
//!    keeping the construction physically sane.
//! 2. **Circulation-matrix accumulator.** The spec's item 4 formula
//!    ("sum(dir)·r̂ per cell") is read most literally as dotting the cell's
//!    raw Cartesian direction-vector sum (`dir_accum`, folded over the full
//!    2π azimuthal range per cell) against a single reference r̂. That
//!    doesn't work: r̂(φ) rotates with azimuth, so integrating dir_accum over
//!    a full azimuthal fold cancels its x/y components identically to zero
//!    at the equator (θ=π/2, where r̂ = (cosφ, sinφ, 0) — pure x/y — so ANY
//!    azimuthally-averaged Cartesian sum dotted with a fixed r̂ is zero there
//!    by construction, regardless of actual radial flow). That would make the
//!    flagship diagnostic vacuous exactly where the paper's claim (outflow at
//!    equatorial mid-radii) needs to show up. Fix: a separate `radial_flux`
//!    accumulator is deposited per-photon as `weight*ds*dir·r̂_local`, where
//!    `r̂_local` is evaluated at THAT photon's own position (i.e. before any
//!    azimuthal folding). This is the only construction that survives the
//!    fold and gives a meaningful signed radial-flow indicator. `dir_accum`
//!    itself is kept exactly as specified for the collision partner-sampling
//!    / coherence use (an explicitly-sanctioned simplification there).
//!
//! Both are documented here and repeated in the report test's header comment
//! so they're visible without spelunking.

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
/// Safety cap on total path length before a photon is given up as lost.
const SAFETY_PATH_CAP: f64 = 20.0 * R_OUT;
const MIN_DS: f64 = 0.02;
const DS_FRACTION: f64 = 0.05;

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

    /// Mean radial-flow component in a cell, bounded [-1,1]; +1 = fully
    /// outward, -1 = fully inward, 0 = no net radial bias.
    fn circulation_at_idx(&self, idx: usize) -> f64 {
        let w = self.track_weight[idx];
        if w < 1e-9 {
            0.0
        } else {
            self.radial_flux[idx] / w
        }
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

fn l2_rel_change(
    old_e: &DensityField,
    old_a: &DensityField,
    new_e: &DensityField,
    new_a: &DensityField,
) -> f64 {
    let mut num = 0.0;
    let mut den = 0.0;
    for idx in 0..GRID_CELLS {
        let o = old_e.density_by_idx(idx) + old_a.density_by_idx(idx);
        let n = new_e.density_by_idx(idx) + new_a.density_by_idx(idx);
        num += (n - o).powi(2);
        den += o.powi(2);
    }
    if den > 1e-300 {
        num.sqrt() / den.sqrt()
    } else {
        0.0
    }
}

fn equator_line_integral(field: &DensityField) -> f64 {
    let j = theta_bin_index(FRAC_PI_2);
    let mut total = 0.0;
    for i in 0..N_R {
        let idx = i * N_THETA + j;
        let dr = r_edge(i + 1) - r_edge(i);
        total += field.density_by_idx(idx) * dr;
    }
    total
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

/// Measured proton exit-structure peaks (project memory: 7°/58° dual peaks).
const EXIT_PEAK_LAT_DEG: [f64; 2] = [7.0, 58.0];
const EXIT_SIGMA_DEG: f64 = 5.0;

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

fn cosine_weighted_hemisphere(rng: &mut u64, n: DVec3) -> DVec3 {
    let xi1 = xorshift64(rng);
    let xi2 = xorshift64(rng);
    let r = xi1.sqrt();
    let phi = xi2 * TAU;
    let x = r * phi.cos();
    let y = r * phi.sin();
    let z = (1.0 - xi1).max(0.0).sqrt();
    let helper = if n.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    let t = n.cross(helper).normalize();
    let b = n.cross(t);
    (t * x + b * y + n * z).normalize()
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

fn spawn_ambient(rng: &mut u64) -> (Photon, Option<f64>) {
    let r_hat = rand_unit_vec(rng);
    let pos = r_hat * R_OUT;
    let dir = cosine_weighted_hemisphere(rng, -r_hat);
    let aim_lat = ray_sphere_aim_latitude(pos, dir);
    (
        Photon {
            pos,
            dir,
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

/// Ray-marches one photon through `n_other`'s density field, depositing its
/// own track-length into `own_field`, until absorbed (r<R_IN), escaped
/// (r>R_OUT), or lost (safety path-length cap).
fn march(
    mut photon: Photon,
    n_other: &DensityField,
    sigma_eff: f64,
    n_ricochet: u32,
    rng: &mut u64,
    own_field: &mut DensityField,
) -> Termination {
    let mut path_length = 0.0;
    loop {
        let r = photon.pos.length();
        if r < R_IN {
            return Termination::Absorbed {
                pos: photon.pos,
                dir: photon.dir,
            };
        }
        if r > R_OUT {
            return Termination::Escaped {
                pos: photon.pos,
                dir: photon.dir,
            };
        }
        if path_length > SAFETY_PATH_CAP {
            return Termination::Lost;
        }

        let ds = (DS_FRACTION * r).max(MIN_DS);
        let theta = colatitude_of(photon.pos, r);
        own_field.deposit(photon.pos, photon.dir, photon.weight, ds);

        if sigma_eff > 0.0 && photon.bounces < n_ricochet {
            let n_density = n_other.density_at(r, theta);
            let p_collide = 1.0 - (-sigma_eff * n_density * ds).exp();
            if xorshift64(rng) < p_collide {
                let (mean_dir, coherence) = n_other.mean_dir_and_coherence_at(r, theta);
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
                    if let Some(d) = new_dir {
                        photon.dir = d;
                    }
                }
                photon.bounces += 1;
            }
        }

        photon.pos += photon.dir * ds;
        path_length += ds;
    }
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
    // ambient-only:
    aim_eq_count: f64,
    aim_eq_arrived: f64,
    aim_eq_arrived_polar: f64,
    aim_miss_count: f64,
    aim_miss_arrived: f64,
    aim_miss_arrived_band: [f64; 3],
}

impl PassTally {
    fn record_emitted(&mut self, term: &Termination) {
        self.total += 1.0;
        match term {
            Termination::Absorbed { pos, dir } => {
                self.arrived += 1.0;
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

    fn record_ambient(&mut self, term: &Termination, aim_lat: Option<f64>) {
        self.total += 1.0;
        let aim_band = aim_lat.map(lat_band);
        match aim_band {
            Some(Band::Equatorial) => self.aim_eq_count += 1.0,
            None => self.aim_miss_count += 1.0,
            _ => {}
        }
        match term {
            Termination::Absorbed { pos, dir } => {
                self.arrived += 1.0;
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
    n_other: &DensityField,
    sigma_eff: f64,
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
                let term = march(photon, n_other, sigma_eff, n_ricochet, rng, &mut field);
                tally.record_emitted(&term);
            }
            PhotonKind::Ambient => {
                let (photon, aim_lat) = spawn_ambient(rng);
                let term = march(photon, n_other, sigma_eff, n_ricochet, rng, &mut field);
                tally.record_ambient(&term, aim_lat);
            }
        }
    }
    (field, tally, lz_source_total)
}

fn run_iteration0(n: usize, swing_boost: f64, rng: &mut u64) -> (DensityField, DensityField) {
    let empty = DensityField::new();
    let (emitted0, _, _) = run_pass(PhotonKind::Emitted, n, swing_boost, &empty, 0.0, 0, rng);
    let (ambient0, _, _) = run_pass(PhotonKind::Ambient, n, swing_boost, &empty, 0.0, 0, rng);
    (emitted0, ambient0)
}

struct ReportCell {
    tau_target: f64,
    n_ricochet: u32,
    sigma_eff: f64,
    convergence: Vec<f64>,
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
    let base_integral = equator_line_integral(emitted0);
    let sigma_eff = if tau_target <= 0.0 || base_integral <= 0.0 {
        0.0
    } else {
        tau_target / base_integral
    };

    let mut rng = seed;
    let mut density_emitted = emitted0.clone();
    let mut density_ambient = ambient0.clone();
    let mut convergence = Vec::with_capacity(iters);
    let mut emitted_tally = PassTally::default();
    let mut ambient_tally = PassTally::default();
    let mut lz_source_total = 0.0;
    let mut final_emitted_field = density_emitted.clone();
    let mut final_ambient_field = density_ambient.clone();

    for _ in 0..iters {
        let (measured_emitted, e_tally, lz_src) = run_pass(
            PhotonKind::Emitted,
            n,
            swing_boost,
            &density_ambient,
            sigma_eff,
            n_ricochet,
            &mut rng,
        );
        let (measured_ambient, a_tally, _) = run_pass(
            PhotonKind::Ambient,
            n,
            swing_boost,
            &density_emitted,
            sigma_eff,
            n_ricochet,
            &mut rng,
        );
        let new_emitted = density_emitted.relax(&measured_emitted, 0.5);
        let new_ambient = density_ambient.relax(&measured_ambient, 0.5);
        convergence.push(l2_rel_change(
            &density_emitted,
            &density_ambient,
            &new_emitted,
            &new_ambient,
        ));
        density_emitted = new_emitted;
        density_ambient = new_ambient;
        emitted_tally = e_tally;
        ambient_tally = a_tally;
        lz_source_total = lz_src;
        final_emitted_field = measured_emitted;
        final_ambient_field = measured_ambient;
    }

    ReportCell {
        tau_target,
        n_ricochet,
        sigma_eff,
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

    /// (a) TAU=0 (sigma_eff=0, no collisions): every ambient photon either
    /// hits the sphere along its straight aim, or escapes — aim latitude ==
    /// arrival latitude within bin resolution; zero emitted re-entries.
    ///
    /// NOTE: the production ray-march step (`ds = 0.05*r`, min 0.02 — the
    /// spec'd adaptive scheme, kept as-is here since the report needs it
    /// fast across 21M+ photon transports) is deliberately coarse far from
    /// the surface (ds~1.5 near r=R_OUT). A small fraction of ambient aims
    /// are near-tangent grazes whose entire dip below r=1 is thinner than
    /// one step — the discretized path hops clean over the crossing and
    /// comes out the other side registering as "escaped" even though the
    /// continuous-geometry aim says "hit". This is a genuine discretization
    /// floor of the fixed step size, not a transport-logic bug, so this test
    /// tolerates a small fraction of such mismatches rather than demanding
    /// bit-for-bit agreement on every grazing incidence.
    #[test]
    fn tau0_transport_matches_straight_line_geometry() {
        let mut rng = 0xA5A5_1234_9876_5432u64;
        let empty = DensityField::new();
        let n = 5000;

        let mut checked = 0;
        let mut hits = 0;
        let mut mismatches = 0;
        for _ in 0..n {
            let (photon, aim_lat) = spawn_ambient(&mut rng);
            let mut field = DensityField::new();
            let term = march(photon, &empty, 0.0, 0, &mut rng, &mut field);
            match (aim_lat, &term) {
                (Some(aim), Termination::Absorbed { pos, .. }) => {
                    let arrival = arrival_latitude(*pos);
                    // Bin resolution ~ 180/18 = 10 deg; allow a little slack
                    // for the boundary/edge cases near tangential aims.
                    if (arrival - aim).abs() >= 5.0 {
                        mismatches += 1;
                    }
                    hits += 1;
                    checked += 1;
                }
                (Some(_aim), Termination::Escaped { .. }) => {
                    // Grazing-incidence discretization floor — see note above.
                    mismatches += 1;
                    checked += 1;
                }
                (None, Termination::Escaped { .. }) => {
                    checked += 1;
                }
                (None, Termination::Absorbed { .. }) => {
                    // Same discretization floor in reverse: a near-tangent
                    // MISS whose closest approach is just barely above 1 can
                    // still register a hit if the step lands inside r<1.
                    mismatches += 1;
                    checked += 1;
                }
                (_, Termination::Lost) => {
                    panic!("photon lost under TAU=0 (no collisions, should never hit the path cap)");
                }
            }
        }
        assert_eq!(checked, n);
        assert!(hits > 0, "expected at least some ambient photons to geometrically hit the sphere");
        let mismatch_frac = mismatches as f64 / n as f64;
        assert!(
            mismatch_frac < 0.03,
            "too many aim/arrival mismatches for a sigma_eff=0 straight-line transport: {mismatches}/{n} ({:.2}%)",
            mismatch_frac * 100.0
        );

        // Emitted photons at TAU=0 must never re-enter r<1 (straight lines
        // from the surface, always genuinely outward — see module doc
        // deviation #1 for why this is guaranteed by construction).
        let mut reentries = 0;
        for _ in 0..n {
            let photon = spawn_emitted(&mut rng, DEFAULT_SWING_BOOST);
            let mut field = DensityField::new();
            let term = march(photon, &empty, 0.0, 0, &mut rng, &mut field);
            if matches!(term, Termination::Absorbed { .. }) {
                reentries += 1;
            }
            assert!(!matches!(term, Termination::Lost), "emitted photon lost under TAU=0");
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
    /// of the naive ~29-30x expected for the (8,20) and (4,16) bin choices,
    /// stable across N=150k..500k, i.e. a real geometric effect, not noise).
    /// Fix: aggregate over the *report's own* equatorial band definition
    /// (|lat|<30°, theta bins 6..=11) so the comparison is symmetric about
    /// the equator and has ~6x the single-bin statistics, then compare two
    /// well-separated, non-adjacent-to-source radial bins.
    #[test]
    fn free_streaming_emitted_density_falls_like_inverse_r_squared_on_equator() {
        let mut rng = 0x1357_2468_ABCD_EF01u64;
        let empty = DensityField::new();
        let (field, _tally, _lz) =
            run_pass(PhotonKind::Emitted, 300_000, DEFAULT_SWING_BOOST, &empty, 0.0, 0, &mut rng);

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

    /// (c) swing_boost=0, TAU=0: net delivered L_z per photon ~ 0 (below
    /// noise floor for N=20k) — no systematic angular-momentum source means
    /// no systematic delivery.
    ///
    /// NOTE: per-photon L_z = (pos×dir).z, and most of the tallied events are
    /// ambient photons escaping back out near r=R_OUT=30, where a single
    /// event's L_z is O(10) (impact-parameter scale ~ r·sinθ). With N=20,000
    /// such events the standard error of the mean is ~10/sqrt(20000)≈0.07 —
    /// an empirically-measured run at this seed landed at -0.0796, matching
    /// that estimate. 0.05 was too tight for pure sampling noise at this N;
    /// 0.2 gives ~3x headroom over the observed value while still being a
    /// small fraction of a typical single-collision impact-parameter scale.
    #[test]
    fn zero_swing_zero_tau_delivers_no_net_angular_momentum() {
        let mut rng = 0xFEED_BEEF_1234_5678u64;
        let n = 20_000usize;
        let empty = DensityField::new();
        let (_field_e, emitted_tally, _lz_src) =
            run_pass(PhotonKind::Emitted, n, 0.0, &empty, 0.0, 0, &mut rng);
        let (_field_a, ambient_tally, _) =
            run_pass(PhotonKind::Ambient, n, 0.0, &empty, 0.0, 0, &mut rng);

        let delivered = ambient_tally.lz_arrived + emitted_tally.lz_arrived;
        let carried_off = ambient_tally.lz_escaped + emitted_tally.lz_escaped;
        let net = delivered - carried_off;
        let net_per_photon = net / n as f64;
        assert!(
            net_per_photon.abs() < 0.2,
            "expected ~0 net L_z/photon at swing_boost=0, got {net_per_photon:.5}"
        );
    }

    /// The CM-3 report: sweeps TAU (collision optical depth target) x
    /// N_RICOCHET (max bounces), self-consistently settles both populations'
    /// density fields for ITERS passes, and prints the screening,
    /// redirection, recycling, circulation, and spin-sustenance evidence for
    /// the disc-to-pole loop claim. See the module doc comment for the two
    /// deliberate deviations from the literal spec (same-hemisphere emission
    /// direction; local-frame radial-flux accumulator for the circulation
    /// matrix).
    #[test]
    #[ignore]
    fn report_recycling_loop() {
        use std::time::Instant;
        let started = Instant::now();

        let n = DEFAULT_N;
        let swing_boost = DEFAULT_SWING_BOOST;
        let mut seed = 0x9E37_79B9_7F4A_7C15u64;

        println!("=== CM-3 Recycling Loop Report ===");
        println!(
            "N per population = {n}, ITERS = {ITERS}, grid = {N_R}x{N_THETA} (r,theta), R_OUT = {R_OUT}, swing_boost = {swing_boost}"
        );

        let (emitted0, ambient0) = run_iteration0(n, swing_boost, &mut seed);
        let base_integral = equator_line_integral(&emitted0);
        println!(
            "iteration-0 equator line integral (emitted density, r=1..{R_OUT}) = {base_integral:.6}"
        );

        let tau_values = [0.0, 0.3, 1.0, 3.0, 10.0];
        let ricochet_values = [1u32, 3, 10];

        let mut cells: Vec<ReportCell> = Vec::new();
        for &tau in &tau_values {
            if tau == 0.0 {
                println!("TAU=0.0 -> sigma_eff = 0.000000 (control, ricochet sweep skipped)");
                seed = seed.wrapping_add(1);
                cells.push(run_self_consistent(
                    n, swing_boost, tau, 0, ITERS, &emitted0, &ambient0, seed,
                ));
            } else {
                let sigma_preview = tau / base_integral;
                println!("TAU={tau:.1} -> sigma_eff = {sigma_preview:.6}");
                for &nr in &ricochet_values {
                    seed = seed.wrapping_add(1);
                    cells.push(run_self_consistent(
                        n, swing_boost, tau, nr, ITERS, &emitted0, &ambient0, seed,
                    ));
                }
            }
        }

        println!();
        println!(
            "{:>6} {:>4} {:>10} | {:>8} {:>8} {:>8} {:>8} | {:>10} {:>9} {:>8} {:>8} {:>8} | {:>9} {:>9} | {:>11} {:>10} {:>6}",
            "TAU", "Nric", "sigma_eff",
            "amb_arr", "polar", "mid", "equat",
            "redirEq>P", "missCapt", "missPol", "missMid", "missEq",
            "em_reent", "em_esc",
            "netLz/ph", "fbRatio", "sign"
        );
        for cell in &cells {
            let a = &cell.ambient_tally;
            let e = &cell.emitted_tally;

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
                "{:>6.1} {:>4} {:>10.6} | {:>8.4} {:>8.4} {:>8.4} {:>8.4} | {:>10.4} {:>9.4} {:>8.4} {:>8.4} {:>8.4} | {:>9.4} {:>9.4} | {:>11.6} {:>10.4} {:>6}",
                cell.tau_target, nric_label, cell.sigma_eff,
                amb_arr_frac, polar_frac, mid_frac, eq_frac,
                redirect, miss_capture, miss_pol, miss_mid, miss_eq,
                em_reentry, em_escape,
                net_per_photon, fb_ratio, sign_label
            );
        }

        let control = cells
            .iter()
            .find(|c| c.tau_target == 0.0)
            .expect("TAU=0 control cell present");
        let flagship = cells
            .iter()
            .find(|c| c.tau_target == 3.0 && c.n_ricochet == 10)
            .expect("TAU=3, N_RICOCHET=10 flagship cell present");

        for (label, cell) in [("CONTROL (TAU=0)", control), ("FLAGSHIP (TAU=3, N_RICOCHET=10)", flagship)] {
            println!("\n--- {label}: per-iteration convergence (L2 relative change of combined density field) ---");
            for (k, c) in cell.convergence.iter().enumerate() {
                println!("  iter {}: {:.6}", k + 1, c);
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
                    let combined_w =
                        cell.final_emitted_field.track_weight[idx] + cell.final_ambient_field.track_weight[idx];
                    let combined_flux =
                        cell.final_emitted_field.radial_flux[idx] + cell.final_ambient_field.radial_flux[idx];
                    let val = if combined_w > 1e-9 { combined_flux / combined_w } else { 0.0 };
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
