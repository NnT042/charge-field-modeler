//! Per-cell precompute (doc §5): free volume (Monte Carlo sampling around
//! body overlap), nearest body (D3 gear frame), and the within-B_MAX flag
//! (D5). Done once on CPU at engine construction; cheap relative to a run
//! (a few thousand cells at gate scale, ~885k at the doc's physics-run h).

use glam::DVec3;
use rand::{RngExt, SeedableRng};

use super::config::{BodyConfig, Boundary, DsmcConfig};

pub struct Grid {
    pub dim: [u32; 3],
    pub min: DVec3,
    pub h: f64,
    pub num_cells: u32,
}

impl Grid {
    pub fn new(config: &DsmcConfig) -> Self {
        let (center, extent) = match config.boundary {
            Boundary::SphereEscape { r_d } => (Self::centroid(config), r_d),
            Boundary::PeriodicBox { half_extent } => (Self::centroid(config), half_extent),
        };
        let side = 2.0 * extent;
        let dim_scalar = (side / config.h).ceil().max(1.0) as u32;
        let dim = [dim_scalar, dim_scalar, dim_scalar];
        let min = center - DVec3::splat(dim_scalar as f64 * config.h / 2.0);
        let num_cells = dim[0] * dim[1] * dim[2];
        Grid { dim, min, h: config.h, num_cells }
    }

    pub(crate) fn centroid(config: &DsmcConfig) -> DVec3 {
        if config.bodies.is_empty() {
            return DVec3::ZERO;
        }
        let sum: DVec3 = config.bodies.iter().map(|b| b.center).sum();
        sum / config.bodies.len() as f64
    }

    pub fn cell_center(&self, cx: u32, cy: u32, cz: u32) -> DVec3 {
        self.min
            + DVec3::new(
                (cx as f64 + 0.5) * self.h,
                (cy as f64 + 0.5) * self.h,
                (cz as f64 + 0.5) * self.h,
            )
    }
}

pub struct PerCellData {
    pub free_volume: Vec<f32>,
    pub nearest_body: Vec<u32>,
    pub within_b_max: Vec<u32>,
}

const MC_SAMPLES: usize = 32;

pub fn precompute(grid: &Grid, bodies: &[BodyConfig], b_max: f64, seed: u64) -> PerCellData {
    let n = grid.num_cells as usize;
    let mut free_volume = vec![(grid.h * grid.h * grid.h) as f32; n];
    let mut nearest_body = vec![0u32; n];
    let mut within_b_max = vec![1u32; n];

    if bodies.is_empty() {
        return PerCellData { free_volume, nearest_body, within_b_max };
    }

    let mut rng = rand::rngs::StdRng::seed_from_u64(seed);
    let mut idx = 0usize;
    for cz in 0..grid.dim[2] {
        for cy in 0..grid.dim[1] {
            for cx in 0..grid.dim[0] {
                let center = grid.cell_center(cx, cy, cz);

                // Nearest body + within-B_MAX.
                let mut best = 0usize;
                let mut best_d = f64::INFINITY;
                for (bi, b) in bodies.iter().enumerate() {
                    let d = (center - b.center).length();
                    if d < best_d {
                        best_d = d;
                        best = bi;
                    }
                }
                nearest_body[idx] = best as u32;
                within_b_max[idx] = if best_d <= b_max { 1 } else { 0 };

                // Free volume: Monte Carlo sample the cell cube, count
                // fraction outside every body's unit sphere (R_IN = 1).
                // Skip sampling (cell fully outside everything, the common
                // case at gate domain sizes) unless the cell is close enough
                // to some body to plausibly overlap it.
                if best_d < 1.0 + grid.h {
                    let half = grid.h / 2.0;
                    let mut outside = 0usize;
                    for _ in 0..MC_SAMPLES {
                        let p = center
                            + DVec3::new(
                                rng.random_range(-half..half),
                                rng.random_range(-half..half),
                                rng.random_range(-half..half),
                            );
                        let inside_any = bodies.iter().any(|b| (p - b.center).length() < 1.0);
                        if !inside_any {
                            outside += 1;
                        }
                    }
                    free_volume[idx] = (grid.h * grid.h * grid.h) as f32 * (outside as f32 / MC_SAMPLES as f32);
                }

                idx += 1;
            }
        }
    }

    PerCellData { free_volume, nearest_body, within_b_max }
}
