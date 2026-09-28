//! Density map post-processing (doc §8): wraps `DsmcEngine::read_density_raw`
//! / `SteadyReport::density_raw` (per-cell [kind x species] u32 counts,
//! accumulated over the averaging window only) into number-density queries --
//! a ray sample (used by `calibrate_sigma_for_tau_in`, doc §11c work item 4)
//! and a CSV dump of the xz-plane slice (doc §8's densification profile).

use std::io::Write as _;

use glam::DVec3;

use super::precompute::Grid;

pub const KIND_AMBIENT: usize = 0;
pub const KIND_EMITTED: usize = 1;
pub const SPECIES_MINUS: usize = 0;
pub const SPECIES_PLUS: usize = 1;

/// Time-averaged photon density over one steady-state averaging window (doc
/// §8: "time-averaged photon counts per cell by kind and species").
pub struct DensityMap {
    dim: [u32; 3],
    min: DVec3,
    h: f64,
    /// [cell*4 + kind*2 + species], raw accumulated counts.
    counts: Vec<u32>,
    /// Free volume per cell (doc §5's Monte-Carlo body-overlap correction),
    /// same units as `h^3`.
    free_volume: Vec<f32>,
    avg_steps: u32,
}

impl DensityMap {
    pub fn new(grid: &Grid, free_volume: Vec<f32>, counts: Vec<u32>, avg_steps: u32) -> Self {
        let num_cells = grid.num_cells as usize;
        assert_eq!(counts.len(), num_cells * 4, "DensityMap: counts length must be num_cells*4");
        assert_eq!(free_volume.len(), num_cells, "DensityMap: free_volume length must be num_cells");
        Self { dim: grid.dim, min: grid.min, h: grid.h, counts, free_volume, avg_steps }
    }

    /// u32 overflow headroom check (doc §8 work item 3: "u32 won't overflow
    /// for <= 10k steps at typical cell counts -- check and assert"). Never
    /// silent: panics loudly rather than returning a quietly-wrong density,
    /// matching the project's TALLY_CELL_OVERFLOW convention.
    pub fn assert_no_overflow(&self) {
        let max_count = self.counts.iter().copied().max().unwrap_or(0);
        assert!(
            max_count < u32::MAX / 4,
            "DensityMap: a cell accumulated {max_count} counts over {} steps -- within striking distance of u32::MAX ({}); refine h or shorten the averaging window",
            self.avg_steps,
            u32::MAX
        );
    }

    fn cell_indices(&self, p: DVec3) -> [i64; 3] {
        let rel = (p - self.min) / self.h;
        [rel.x.floor() as i64, rel.y.floor() as i64, rel.z.floor() as i64]
    }

    /// Nearest-cell lookup, clamped to the grid (a ray endpoint slightly
    /// outside the grid due to f32/f64 rounding still gets a sane answer).
    fn cell_index_clamped(&self, p: DVec3) -> usize {
        let [cx, cy, cz] = self.cell_indices(p);
        let cx = cx.clamp(0, self.dim[0] as i64 - 1) as u32;
        let cy = cy.clamp(0, self.dim[1] as i64 - 1) as u32;
        let cz = cz.clamp(0, self.dim[2] as i64 - 1) as u32;
        (cx + cy * self.dim[0] + cz * self.dim[0] * self.dim[1]) as usize
    }

    pub fn count_at(&self, p: DVec3, kind: usize, species: usize) -> u32 {
        let cell = self.cell_index_clamped(p);
        self.counts[cell * 4 + kind * 2 + species]
    }

    /// Number density (photons per unit volume) at `p`: count / (avg_steps *
    /// free_volume). `free_volume` already excludes body overlap (doc §5),
    /// so a cell fully inside a body (free_volume ~ 0) returns 0 rather than
    /// a divide blow-up.
    pub fn density_at(&self, p: DVec3, kind: usize, species: usize) -> f64 {
        let cell = self.cell_index_clamped(p);
        let vcell = self.free_volume[cell] as f64;
        if vcell < 1e-9 {
            return 0.0;
        }
        self.counts[cell * 4 + kind * 2 + species] as f64 / (self.avg_steps.max(1) as f64 * vcell)
    }

    /// Density summed over species (ambient+emitted convenience: pass the
    /// desired `kind`, both species).
    pub fn density_at_both_species(&self, p: DVec3, kind: usize) -> f64 {
        self.density_at(p, kind, SPECIES_MINUS) + self.density_at(p, kind, SPECIES_PLUS)
    }

    /// Samples density at `n` equally spaced points from `p0` to `p1`
    /// (inclusive), returning (distance from p0, density) pairs. Used for
    /// TAU_in/TAU_out (doc §11c work item 4: the equatorial b=0 line).
    pub fn sample_ray(&self, p0: DVec3, p1: DVec3, n: usize, kind: usize, species: usize) -> Vec<(f64, f64)> {
        assert!(n >= 2, "sample_ray: need at least 2 samples");
        let seg = p1 - p0;
        let len = seg.length();
        (0..n)
            .map(|i| {
                let t = i as f64 / (n - 1) as f64;
                let p = p0 + seg * t;
                (t * len, self.density_at(p, kind, species))
            })
            .collect()
    }

    /// Trapezoidal integral of density along the ray p0->p1 -- the "TAU = sigma
    /// * integral(n ds)" building block (doc §11c work item 4). Caller
    /// multiplies by sigma (and any rate-kernel closing-speed factor).
    pub fn integrate_ray(&self, p0: DVec3, p1: DVec3, n: usize, kind: usize, species: usize) -> f64 {
        let samples = self.sample_ray(p0, p1, n, kind, species);
        let mut integral = 0.0;
        for w in samples.windows(2) {
            let (d0, n0) = w[0];
            let (d1, n1) = w[1];
            integral += 0.5 * (n0 + n1) * (d1 - d0);
        }
        integral
    }

    /// Dumps the xz-plane slice through the cell row nearest `y0` as CSV
    /// (`x,z,count_ambient_minus,count_ambient_plus,count_emitted_minus,count_emitted_plus`),
    /// doc §8's densification profile.
    pub fn dump_xz_csv(&self, path: &str, y0: f64) -> std::io::Result<()> {
        let cy = (((y0 - self.min.y) / self.h).floor() as i64).clamp(0, self.dim[1] as i64 - 1) as u32;
        let mut f = std::fs::File::create(path)?;
        writeln!(f, "x,z,count_amb_minus,count_amb_plus,count_emit_minus,count_emit_plus")?;
        for cz in 0..self.dim[2] {
            for cx in 0..self.dim[0] {
                let cell = (cx + cy * self.dim[0] + cz * self.dim[0] * self.dim[1]) as usize;
                let x = self.min.x + (cx as f64 + 0.5) * self.h;
                let z = self.min.z + (cz as f64 + 0.5) * self.h;
                let base = cell * 4;
                writeln!(
                    f,
                    "{x},{z},{},{},{},{}",
                    self.counts[base + KIND_AMBIENT * 2 + SPECIES_MINUS],
                    self.counts[base + KIND_AMBIENT * 2 + SPECIES_PLUS],
                    self.counts[base + KIND_EMITTED * 2 + SPECIES_MINUS],
                    self.counts[base + KIND_EMITTED * 2 + SPECIES_PLUS],
                )?;
            }
        }
        Ok(())
    }

    pub fn avg_steps(&self) -> u32 {
        self.avg_steps
    }
}
