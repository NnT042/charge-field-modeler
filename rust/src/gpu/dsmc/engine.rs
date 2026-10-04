//! Host orchestration of the DSMC step (doc §7). B1 exposes three
//! primitives rather than one fixed step function, because the two B1 gates
//! need different subsets of the pipeline:
//!
//! - [`DsmcEngine::inject_batch`] -- a one-off inject+emit burst.
//! - [`DsmcEngine::step_transport_only`] -- move only, no sort/gather/collide
//!   (G0: sigma = 0, so collisions never fire and there is nothing to sort
//!   for -- a burst of independent photons just needs marching to
//!   termination).
//! - [`DsmcEngine::step_collide`] -- the full move -> sort -> gather ->
//!   collide pipeline (G1, and B2/B3's continuous steady-state loop).
//!
//! Both step functions read back the small tally buffer each call (doc's
//! "Tally readback" -- a few hundred bytes, not the photon buffer) and fold
//! it into host `i64` totals so a run of thousands of steps never overflows
//! the per-step fixed-point accumulators.

use bytemuck::Zeroable;

use crate::gpu::context::GpuContext;
use crate::gpu::kernel::{workgroup_grid_for_threads, Kernel};
use crate::gpu::sort::CellSorter;

use super::config::{BodyModel, Boundary, DsmcConfig, GpuClass, EmissionMode, ExchangeRule, GpuOrient, MeshMode, RateKernel, SpeedMode, StiffnessShed};
use super::gpu_types::{GpuBody, GpuInjectParams, GpuParams, GpuPhoton, HIST_LEN, TALLY_LEN};
use super::precompute::{precompute, Grid};

const COMMON_WGSL: &str = include_str!("../../../shaders/dsmc_common.wgsl");

fn load(ctx: &GpuContext, common: bool, src: &str) -> Kernel {
    let full = if common { format!("{COMMON_WGSL}\n{src}") } else { src.to_string() };
    let module = ctx.load_wgsl_compute(&full);
    Kernel::new(ctx, &module, "main")
}

/// Host-side running totals, folded from the GPU tally readback each step.
#[derive(Default, Clone, Debug)]
pub struct Totals {
    pub absorbed: [[[i64; 2]; 2]; 4],
    pub absorbed_crossed: [[[i64; 2]; 2]; 4],
    pub crossings: [[i64; 2]; 4],
    pub escaped: [[[i64; 2]; 2]; 4],
    pub emitted: [i64; 4],
    pub injected: i64,
    pub collisions: i64,
    pub mom_drift_count: i64,
    pub dir_sum: [f64; 3],
    pub dz2_sum: f64,
    pub mom_drift_sum: f64,
    pub dir_sum_samples: i64,
    pub cell_overflow: i64,
    pub exchange_drift_sum: f64,
    pub redirect_drift_sum: f64,
    pub low_rank64_collisions: i64,

    // --- B2b work item 3 (doc §8 observables) ---
    /// Sigma dir.xyz per body -- absorbed momentum (open-field absorber force).
    pub absorbed_mom: [[f64; 3]; 4],
    /// Escapes by source-body-SLOT (0..3 = that body, 4 = "none"/ambient) x
    /// species x |latitude| band (0 polar, 1 mid, 2 equatorial).
    pub escaped_band: [[[i64; 3]; 2]; 5],
    pub escaped_band_spin: [[[f64; 3]; 2]; 5],
    /// Of escaped photons, how many had crossed body i's lane at some point
    /// (NOT mutually exclusive across i, see gpu_types.rs doc).
    pub escaped_crossed_body: [[i64; 2]; 4],
    pub escaped_crossed_body_spin: [[f64; 2]; 4],
    /// Counter-stream survival: [crossed_body][arrived_body][collided 0/1].
    pub survival: [[[i64; 2]; 4]; 4],

    // --- E1: emitted-photon stiffness (doc §11d) ---
    /// Sum/count of stiffness level k over escaped EMITTED photons (report's
    /// "mean k of escaping emitted photons").
    pub escaped_emitted_k_sum: i64,
    pub escaped_emitted_k_count: i64,
    /// Per-photon deflection angle (radians), bucketed by that photon's own
    /// pre-collision stiffness level k (0..7): sum, sum-of-squares, count.
    pub deflect_sum: [f64; 8],
    pub deflect_sumsq: [f64; 8],
    pub deflect_count: [i64; 8],

    // --- Doc §13: rough-sphere photons ---
    pub accept_overflow: i64,
    pub pin_deficit_count: i64,
    pub pin_deficit_sum: f64,
    pub edge_over_c: i64,
    pub esc_all: i64,
    pub esc_speed_sum: f64,
    pub esc_speed_sumsq: f64,
    pub esc_edge_sum: f64,
    pub esc_edge_sumsq: f64,
    pub esc_edge_over_c: i64,
    pub twist_to_normal: f64,
    pub stack_events: i64,
    pub unstack_events: i64,
    pub stack_capped: i64,
    pub esc_level: [i64; 4],

    // --- Doc §14: permeable shell bodies (per body) ---
    pub sh_entries: [i64; 4],
    pub sh_exits: [i64; 4],
    pub sh_through: [i64; 4],
    pub sh_contacts: [i64; 4],
    pub sh_clean: [i64; 4],
    pub sh_exit_band: [[i64; 3]; 4],
    pub sh_exit_band_hel: [[f64; 3]; 4],
    pub sh_dp: [[f64; 3]; 4],
    pub sh_torque: [[f64; 3]; 4],
    pub sh_work: [f64; 4],

    // --- Doc §15 ball ---
    pub ball_overlap: i64,
    /// [angle bin][outcome][species]
    pub pass: [[[i64; 2]; 4]; 9],
    pub pass_through: [[i64; 2]; 9],
    pub pass_hel: [[f64; 2]; 4],
    /// Doc §15.3: body-0 contact-point latitude, 180 one-degree bins
    /// (empty until the first readback).
    pub hit_lat: Vec<i64>,
    /// Doc §15.3: body-0 contacted-exit latitude, [spec][180 bins] flattened.
    pub exit_lat: Vec<i64>,
    /// Doc §15.6: body-0 shell crossings, gpu_types::TALLY_FLUX layout.
    pub flux: Vec<i64>,
    /// Doc §15.9: maze outcomes (gpu_types::TALLY_MAZE layout), emission
    /// latitude (180 bins) and escape azimuth (36 bins), body 0.
    pub maze: [i64; 5],
    pub maze_lat: Vec<i64>,
    pub maze_az: Vec<i64>,
    /// Doc §15.10: captures [body][emitter tag] and the two-body line flux
    /// (gpu_types::TALLY_LINE layout).
    pub feed: [[i64; 5]; 4],
    pub line: Vec<i64>,
    /// Doc §15.11: tank in/out per body, capture latitude (body 0, by
    /// never-captured / recaptured), envelope captures per body.
    pub hold: [[i64; 2]; 4],
    pub capt_lat: Vec<i64>,
    pub env_capt: [i64; 4],
    /// Doc §15.24: photons passed through the through-charge hole.
    pub hole_pass: [i64; 4],
    /// Doc §15.13: cell-steps run in parallel collision rounds.
    pub par_cells: i64,
}

impl Totals {
    /// Adds every field of `other` into `self`. Used by the steady-state
    /// driver (B2b work item 1) to maintain two running accumulators (the
    /// whole averaging window, and the current reporting block) from the same
    /// per-step fold without reading the GPU tally buffer twice.
    pub fn add(&mut self, other: &Totals) {
        for body in 0..4 {
            for kind in 0..2 {
                for species in 0..2 {
                    self.absorbed[body][kind][species] += other.absorbed[body][kind][species];
                    self.absorbed_crossed[body][kind][species] += other.absorbed_crossed[body][kind][species];
                    self.escaped[body][kind][species] += other.escaped[body][kind][species];
                }
            }
            for species in 0..2 {
                self.crossings[body][species] += other.crossings[body][species];
            }
            self.emitted[body] += other.emitted[body];
        }
        self.injected += other.injected;
        self.collisions += other.collisions;
        self.mom_drift_count += other.mom_drift_count;
        for k in 0..3 {
            self.dir_sum[k] += other.dir_sum[k];
        }
        self.dz2_sum += other.dz2_sum;
        self.mom_drift_sum += other.mom_drift_sum;
        self.dir_sum_samples += other.dir_sum_samples;
        self.cell_overflow += other.cell_overflow;
        self.exchange_drift_sum += other.exchange_drift_sum;
        self.redirect_drift_sum += other.redirect_drift_sum;
        self.low_rank64_collisions += other.low_rank64_collisions;

        for body in 0..4 {
            for a in 0..3 {
                self.absorbed_mom[body][a] += other.absorbed_mom[body][a];
            }
        }
        for slot in 0..5 {
            for species in 0..2 {
                for band in 0..3 {
                    self.escaped_band[slot][species][band] += other.escaped_band[slot][species][band];
                    self.escaped_band_spin[slot][species][band] += other.escaped_band_spin[slot][species][band];
                }
            }
        }
        for body in 0..4 {
            for species in 0..2 {
                self.escaped_crossed_body[body][species] += other.escaped_crossed_body[body][species];
                self.escaped_crossed_body_spin[body][species] += other.escaped_crossed_body_spin[body][species];
            }
            for arrived in 0..4 {
                for c in 0..2 {
                    self.survival[body][arrived][c] += other.survival[body][arrived][c];
                }
            }
        }

        self.escaped_emitted_k_sum += other.escaped_emitted_k_sum;
        self.escaped_emitted_k_count += other.escaped_emitted_k_count;
        for k in 0..8 {
            self.deflect_sum[k] += other.deflect_sum[k];
            self.deflect_sumsq[k] += other.deflect_sumsq[k];
            self.deflect_count[k] += other.deflect_count[k];
        }
        self.accept_overflow += other.accept_overflow;
        self.pin_deficit_count += other.pin_deficit_count;
        self.pin_deficit_sum += other.pin_deficit_sum;
        self.edge_over_c += other.edge_over_c;
        self.esc_all += other.esc_all;
        self.esc_speed_sum += other.esc_speed_sum;
        self.esc_speed_sumsq += other.esc_speed_sumsq;
        self.esc_edge_sum += other.esc_edge_sum;
        self.esc_edge_sumsq += other.esc_edge_sumsq;
        self.esc_edge_over_c += other.esc_edge_over_c;
        self.twist_to_normal += other.twist_to_normal;
        self.stack_events += other.stack_events;
        self.unstack_events += other.unstack_events;
        self.stack_capped += other.stack_capped;
        for k in 0..4 {
            self.esc_level[k] += other.esc_level[k];
        }
        for b in 0..4 {
            self.sh_entries[b] += other.sh_entries[b];
            self.sh_exits[b] += other.sh_exits[b];
            self.sh_through[b] += other.sh_through[b];
            self.sh_contacts[b] += other.sh_contacts[b];
            self.sh_clean[b] += other.sh_clean[b];
            self.sh_work[b] += other.sh_work[b];
            for a in 0..3 {
                self.sh_exit_band[b][a] += other.sh_exit_band[b][a];
                self.sh_exit_band_hel[b][a] += other.sh_exit_band_hel[b][a];
                self.sh_dp[b][a] += other.sh_dp[b][a];
                self.sh_torque[b][a] += other.sh_torque[b][a];
            }
        }
        self.ball_overlap += other.ball_overlap;
        for ab in 0..9 {
            for sp in 0..2 {
                self.pass_through[ab][sp] += other.pass_through[ab][sp];
                for o in 0..4 {
                    self.pass[ab][o][sp] += other.pass[ab][o][sp];
                }
            }
        }
        for o in 0..4 {
            for sp in 0..2 {
                self.pass_hel[o][sp] += other.pass_hel[o][sp];
            }
        }
        add_vec(&mut self.hit_lat, &other.hit_lat);
        add_vec(&mut self.exit_lat, &other.exit_lat);
        add_vec(&mut self.flux, &other.flux);
        for i in 0..5 {
            self.maze[i] += other.maze[i];
        }
        add_vec(&mut self.maze_lat, &other.maze_lat);
        add_vec(&mut self.maze_az, &other.maze_az);
        for b in 0..4 {
            for t in 0..5 {
                self.feed[b][t] += other.feed[b][t];
            }
        }
        add_vec(&mut self.line, &other.line);
        for b in 0..4 {
            self.hold[b][0] += other.hold[b][0];
            self.hold[b][1] += other.hold[b][1];
            self.env_capt[b] += other.env_capt[b];
            self.hole_pass[b] += other.hole_pass[b];
        }
        add_vec(&mut self.capt_lat, &other.capt_lat);
        self.par_cells += other.par_cells;
    }

    /// (mean, std) of escaping photons' speed, in c (doc §13, D7).
    pub fn esc_speed_stats(&self) -> (f64, f64) {
        mean_std(self.esc_speed_sum, self.esc_speed_sumsq, self.esc_all)
    }

    /// (mean, std) of escaping photons' spin edge speed |w|R, in c.
    pub fn esc_edge_stats(&self) -> (f64, f64) {
        mean_std(self.esc_edge_sum, self.esc_edge_sumsq, self.esc_all)
    }

    /// Sum of absorbed counts (both kind, both species) for one body -- the
    /// per-step scalar the pump EMA tracks (doc §3.3).
    pub fn absorbed_total(&self, body: usize) -> i64 {
        self.absorbed[body][0][0] + self.absorbed[body][0][1] + self.absorbed[body][1][0] + self.absorbed[body][1][1]
    }

    fn fold(&mut self, raw: &[i32], sample_count: u32) {
        const SCALE: f64 = super::gpu_types::TALLY_SCALE;
        for body in 0..4 {
            for kind in 0..2 {
                for species in 0..2 {
                    self.absorbed[body][kind][species] += raw[body * 4 + kind * 2 + species] as i64;
                    self.absorbed_crossed[body][kind][species] += raw[16 + body * 4 + kind * 2 + species] as i64;
                    self.escaped[body][kind][species] += raw[40 + body * 4 + kind * 2 + species] as i64;
                }
            }
            for species in 0..2 {
                self.crossings[body][species] += raw[32 + body * 2 + species] as i64;
            }
            self.emitted[body] += raw[56 + body] as i64;
        }
        self.injected += raw[60] as i64;
        self.collisions += raw[61] as i64;
        self.mom_drift_count += raw[62] as i64;
        self.dir_sum[0] += raw[63] as f64 / SCALE;
        self.dir_sum[1] += raw[64] as f64 / SCALE;
        self.dir_sum[2] += raw[65] as f64 / SCALE;
        self.dz2_sum += raw[66] as f64 / SCALE;
        self.mom_drift_sum += raw[67] as f64 / SCALE;
        self.cell_overflow += raw[68] as i64;
        // Never silent: a cell above MAX_CELL only pairs its first MAX_CELL
        // photons (by id), which biases exactly the dense gas we care about.
        // Callers should refine `h` if this ever prints.
        if raw[68] > 0 {
            eprintln!(
                "DSMC WARNING: {} cell(s) exceeded MAX_CELL this step -- collisions there are biased; refine h",
                raw[68]
            );
        }
        self.exchange_drift_sum += raw[69] as f64 / SCALE;
        self.redirect_drift_sum += raw[70] as f64 / SCALE;
        self.low_rank64_collisions += raw[71] as i64;
        self.dir_sum_samples += sample_count as i64;

        use super::gpu_types::{
            TALLY_ABSORBED_MOM, TALLY_ESCAPED_BAND, TALLY_ESCAPED_BAND_SPIN, TALLY_ESCAPED_CROSSED_BODY, TALLY_ESCAPED_CROSSED_BODY_SPIN,
            TALLY_SURVIVAL,
        };
        for body in 0..4 {
            for a in 0..3 {
                self.absorbed_mom[body][a] += raw[TALLY_ABSORBED_MOM + body * 3 + a] as f64 / SCALE;
            }
        }
        for slot in 0..5 {
            for species in 0..2 {
                for band in 0..3 {
                    let eb_idx = slot * 2 * 3 + species * 3 + band;
                    self.escaped_band[slot][species][band] += raw[TALLY_ESCAPED_BAND + eb_idx] as i64;
                    self.escaped_band_spin[slot][species][band] += raw[TALLY_ESCAPED_BAND_SPIN + eb_idx] as f64 / SCALE;
                }
            }
        }
        for body in 0..4 {
            for species in 0..2 {
                self.escaped_crossed_body[body][species] += raw[TALLY_ESCAPED_CROSSED_BODY + body * 2 + species] as i64;
                self.escaped_crossed_body_spin[body][species] += raw[TALLY_ESCAPED_CROSSED_BODY_SPIN + body * 2 + species] as f64 / SCALE;
            }
            for arrived in 0..4 {
                for c in 0..2 {
                    self.survival[body][arrived][c] += raw[TALLY_SURVIVAL + (body * 4 + arrived) * 2 + c] as i64;
                }
            }
        }

        use super::gpu_types::{TALLY_DEFLECT_COUNT, TALLY_DEFLECT_SUM, TALLY_DEFLECT_SUMSQ, TALLY_ESCAPED_EMITTED_K_COUNT, TALLY_ESCAPED_EMITTED_K_SUM};
        self.escaped_emitted_k_sum += raw[TALLY_ESCAPED_EMITTED_K_SUM] as i64;
        self.escaped_emitted_k_count += raw[TALLY_ESCAPED_EMITTED_K_COUNT] as i64;
        for k in 0..8 {
            self.deflect_sum[k] += raw[TALLY_DEFLECT_SUM + k] as f64 / SCALE;
            self.deflect_sumsq[k] += raw[TALLY_DEFLECT_SUMSQ + k] as f64 / SCALE;
            self.deflect_count[k] += raw[TALLY_DEFLECT_COUNT + k] as i64;
        }

        use super::gpu_types::{
            TALLY_ACCEPT_OVERFLOW, TALLY_EDGE_OVER_C, TALLY_ESC_ALL, TALLY_ESC_EDGE_OVER_C, TALLY_ESC_EDGE_SUM, TALLY_ESC_EDGE_SUMSQ,
            TALLY_ESC_SPEED_SUM, TALLY_ESC_SPEED_SUMSQ, TALLY_PIN_DEFICIT_COUNT, TALLY_PIN_DEFICIT_SUM, TALLY_TWIST_TO_NORMAL,
        };
        self.accept_overflow += raw[TALLY_ACCEPT_OVERFLOW] as i64;
        if raw[TALLY_ACCEPT_OVERFLOW] > 0 {
            eprintln!(
                "DSMC WARNING: {} collision candidate(s) had |v1-v2| > g_max this step -- rate biased; raise RoughParams::g_max",
                raw[TALLY_ACCEPT_OVERFLOW]
            );
        }
        self.pin_deficit_count += raw[TALLY_PIN_DEFICIT_COUNT] as i64;
        self.pin_deficit_sum += raw[TALLY_PIN_DEFICIT_SUM] as f64 / SCALE;
        self.edge_over_c += raw[TALLY_EDGE_OVER_C] as i64;
        self.esc_all += raw[TALLY_ESC_ALL] as i64;
        self.esc_speed_sum += raw[TALLY_ESC_SPEED_SUM] as f64 / SCALE;
        self.esc_speed_sumsq += raw[TALLY_ESC_SPEED_SUMSQ] as f64 / SCALE;
        self.esc_edge_sum += raw[TALLY_ESC_EDGE_SUM] as f64 / SCALE;
        self.esc_edge_sumsq += raw[TALLY_ESC_EDGE_SUMSQ] as f64 / SCALE;
        self.esc_edge_over_c += raw[TALLY_ESC_EDGE_OVER_C] as i64;
        self.twist_to_normal += raw[TALLY_TWIST_TO_NORMAL] as f64 / SCALE;
        use super::gpu_types::{TALLY_ESC_LEVEL, TALLY_STACK_CAPPED, TALLY_STACK_EVENTS, TALLY_UNSTACK_EVENTS};
        self.stack_events += raw[TALLY_STACK_EVENTS] as i64;
        self.unstack_events += raw[TALLY_UNSTACK_EVENTS] as i64;
        self.stack_capped += raw[TALLY_STACK_CAPPED] as i64;
        for k in 0..4 {
            self.esc_level[k] += raw[TALLY_ESC_LEVEL + k] as i64;
        }
        use super::gpu_types::{
            TALLY_SH_CLEAN, TALLY_SH_CONTACTS, TALLY_SH_DP, TALLY_SH_ENTRIES, TALLY_SH_EXITS, TALLY_SH_EXIT_BAND, TALLY_SH_EXIT_BAND_HEL,
            TALLY_SH_THROUGH, TALLY_SH_TORQUE, TALLY_SH_WORK,
        };
        for b in 0..4 {
            self.sh_entries[b] += raw[TALLY_SH_ENTRIES + b] as i64;
            self.sh_exits[b] += raw[TALLY_SH_EXITS + b] as i64;
            self.sh_through[b] += raw[TALLY_SH_THROUGH + b] as i64;
            self.sh_contacts[b] += raw[TALLY_SH_CONTACTS + b] as i64;
            self.sh_clean[b] += raw[TALLY_SH_CLEAN + b] as i64;
            self.sh_work[b] += raw[TALLY_SH_WORK + b] as f64 / SCALE;
            for a in 0..3 {
                self.sh_exit_band[b][a] += raw[TALLY_SH_EXIT_BAND + b * 3 + a] as i64;
                self.sh_exit_band_hel[b][a] += raw[TALLY_SH_EXIT_BAND_HEL + b * 3 + a] as f64 / SCALE;
                self.sh_dp[b][a] += raw[TALLY_SH_DP + b * 3 + a] as f64 / SCALE;
                self.sh_torque[b][a] += raw[TALLY_SH_TORQUE + b * 3 + a] as f64 / SCALE;
            }
        }
        use super::gpu_types::{TALLY_BALL_OVERLAP, TALLY_PASS, TALLY_PASS_HEL, TALLY_PASS_THROUGH};
        self.ball_overlap += raw[TALLY_BALL_OVERLAP] as i64;
        for ab in 0..9 {
            for sp in 0..2 {
                self.pass_through[ab][sp] += raw[TALLY_PASS_THROUGH + ab * 2 + sp] as i64;
                for o in 0..4 {
                    self.pass[ab][o][sp] += raw[TALLY_PASS + (ab * 4 + o) * 2 + sp] as i64;
                }
            }
        }
        for o in 0..4 {
            for sp in 0..2 {
                self.pass_hel[o][sp] += raw[TALLY_PASS_HEL + o * 2 + sp] as f64 / SCALE;
            }
        }
        use super::gpu_types::{TALLY_EXIT_LAT, TALLY_HIT_LAT};
        let hit: Vec<i64> = raw[TALLY_HIT_LAT..TALLY_HIT_LAT + 180].iter().map(|&x| x as i64).collect();
        let exit: Vec<i64> = raw[TALLY_EXIT_LAT..TALLY_EXIT_LAT + 360].iter().map(|&x| x as i64).collect();
        add_vec(&mut self.hit_lat, &hit);
        add_vec(&mut self.exit_lat, &exit);
        use super::gpu_types::TALLY_FLUX;
        let flux: Vec<i64> = raw[TALLY_FLUX..TALLY_FLUX + 504].iter().map(|&x| x as i64).collect();
        add_vec(&mut self.flux, &flux);
        use super::gpu_types::{TALLY_MAZE, TALLY_MAZE_AZ, TALLY_MAZE_LAT};
        for i in 0..5 {
            self.maze[i] += raw[TALLY_MAZE + i] as i64;
        }
        let ml: Vec<i64> = raw[TALLY_MAZE_LAT..TALLY_MAZE_LAT + 180].iter().map(|&x| x as i64).collect();
        let ma: Vec<i64> = raw[TALLY_MAZE_AZ..TALLY_MAZE_AZ + 36].iter().map(|&x| x as i64).collect();
        add_vec(&mut self.maze_lat, &ml);
        add_vec(&mut self.maze_az, &ma);
        use super::gpu_types::{LINE_N, TALLY_FEED, TALLY_LINE};
        for b in 0..4 {
            for t in 0..5 {
                self.feed[b][t] += raw[TALLY_FEED + b * 5 + t] as i64;
            }
        }
        let ln: Vec<i64> = raw[TALLY_LINE..TALLY_LINE + LINE_N * 6].iter().map(|&x| x as i64).collect();
        add_vec(&mut self.line, &ln);
        use super::gpu_types::{TALLY_CAPT_LAT, TALLY_ENV_CAPT, TALLY_HOLD};
        for b in 0..4 {
            self.hold[b][0] += raw[TALLY_HOLD + b * 2] as i64;
            self.hold[b][1] += raw[TALLY_HOLD + b * 2 + 1] as i64;
            self.env_capt[b] += raw[TALLY_ENV_CAPT + b] as i64;
            self.hole_pass[b] += raw[super::gpu_types::TALLY_HOLE_PASS + b] as i64;
        }
        let cl: Vec<i64> = raw[TALLY_CAPT_LAT..TALLY_CAPT_LAT + 360].iter().map(|&x| x as i64).collect();
        add_vec(&mut self.capt_lat, &cl);
        self.par_cells += raw[super::gpu_types::TALLY_PAR_CELLS] as i64;
    }
}

fn add_vec(a: &mut Vec<i64>, b: &[i64]) {
    if b.is_empty() {
        return;
    }
    if a.is_empty() {
        a.resize(b.len(), 0);
    }
    for (x, y) in a.iter_mut().zip(b) {
        *x += *y;
    }
}

fn mean_std(sum: f64, sumsq: f64, n: i64) -> (f64, f64) {
    if n <= 0 {
        return (f64::NAN, f64::NAN);
    }
    let m = sum / n as f64;
    (m, (sumsq / n as f64 - m * m).max(0.0).sqrt())
}

pub struct DsmcEngine {
    pub config: DsmcConfig,
    pub grid: Grid,
    photons: [wgpu::Buffer; 2],
    cur: usize,
    bodies_buf: wgpu::Buffer,
    /// Doc §15.17: the escape sphere's centre, fixed at construction so
    /// that moving bodies (`set_body_centers`) do not drag the domain.
    domain_ctr: glam::DVec3,
    /// Body centres the nearest-body map was last computed for.
    nb_centers: Vec<glam::DVec3>,
    body_emit_counts_buf: wgpu::Buffer,
    raw_cell_of: wgpu::Buffer,
    nearest_body_buf: wgpu::Buffer,
    within_bmax_buf: wgpu::Buffer,
    free_volume_buf: wgpu::Buffer,
    cell_remainder_buf: wgpu::Buffer,
    tallies_buf: wgpu::Buffer,
    hist_buf: wgpu::Buffer,
    /// Density map (doc §8): per-cell [kind x species] u32 counters, NOT
    /// reset every step like `tallies_buf` -- cleared once by
    /// `begin_density_accum` and read back once at the end of the averaging
    /// window (see `density.rs`).
    density_buf: wgpu::Buffer,
    /// CPU copy of the precompute free-volume map, kept for `read_density`'s
    /// count -> number-density conversion (the GPU-resident copy is
    /// write-once and only used by `dsmc_collide.wgsl`'s NTC rate).
    free_volume_cpu: Vec<f32>,
    density_accum_enabled: bool,
    inject_params_buf: wgpu::Buffer,
    params_buf: wgpu::Buffer,
    gather_params_buf: wgpu::Buffer,
    sorter: CellSorter,
    /// Doc §14: per-body shell kinematics (empty unless BodyModel::Shells)
    /// and the GPU copy of their frames (8 vec4 per body).
    pub shell_bodies: Vec<crate::baryon_shells::ShellBody>,
    /// Doc §15.11: photons currently held in each body's tank.
    pub held: Vec<i64>,
    shells_buf: wgpu::Buffer,
    ballkeys_buf: wgpu::Buffer,

    inject_emit_kernel: Kernel,
    move_kernel: Kernel,
    gather_kernel: Kernel,
    collide_kernel: Kernel,

    pub n_used: u32,
    pub n_live: u32,
    pub step: u32,
    pub id_next: u32,
    ambient_remainder: f64,
    emit_remainder: Vec<f64>,
    /// Pump-mode EMA of each body's absorption rate (doc §3.3); unused (and
    /// unexercised by the B1 gates, both of which use Forced emission) when
    /// `config.emission_mode == Forced`.
    pump_rate_ema: Vec<f64>,
}

impl DsmcEngine {
    pub fn new(ctx: &GpuContext, config: DsmcConfig) -> Self {
        if config.exchange_rule == ExchangeRule::RoughSphere {
            // Doc §13: contact mechanics replace these devices outright; a
            // leftover nonzero value would otherwise be silently ignored.
            assert!(config.cancel_redirect == 0.0, "RoughSphere: cancel_redirect must be 0 (the redirect is retired)");
            assert!(config.spin_transfer == 0.0, "RoughSphere: spin_transfer must be 0 (spin is the vector w)");
            assert!(config.stiffness_shed == StiffnessShed::None, "RoughSphere: stiffness shedding has no energy accounting");
            assert!(config.bounce_cap == 0, "RoughSphere: bounce cap not supported");
            assert!(config.rate_kernel == RateKernel::Catch, "RoughSphere: collision rate must be the kinetic |v1-v2| kernel");
            assert!(config.rough.speed_mode != SpeedMode::PinnedC || config.rough.g_max == 2.0, "PinnedC: g_max is exactly 2");
            assert!(config.rough.stack_max_level <= 3, "stacking storage holds levels 0..=3");
            assert!(!config.rough.stacking || config.emit_stiffness_level == 0, "stacking: photons are born at level 0");
        }
        if config.body_model != BodyModel::Absorber {
            assert!(config.exchange_rule == ExchangeRule::RoughSphere, "Shells bodies use the rough-sphere photon mechanics");
            assert!(config.bodies.len() <= 4, "at most 4 bodies");
            assert!(config.body_model == BodyModel::Ball || config.rough.grain_mode == super::config::GrainMode::EdgeC, "grain modes Rigid/RigidInner need the ball model");
            assert!(config.body_model == BodyModel::Ball || config.rough.maze == super::config::MazeMode::Off, "the maze emission law needs the ball model");
        }
        let grid = Grid::new(&config);
        let n_bodies = config.bodies.len();
        let mut pc = precompute(&grid, &config.bodies, config.b_max, config.seed ^ 0xB0D7);
        if config.body_model != BodyModel::Absorber {
            // Permeable bodies (doc §14/§15) occupy no gas volume: the cells
            // inside them hold gas like any other. Carving the body spheres
            // out (the Absorber rule) left near-zero free volume there and
            // inflated the photon-photon rate inside every envelope, capped
            // only by the per-cell candidate limit. Found 2026-09-29 by the
            // §15 passage tallies; the §14 shell results before the fix are
            // contaminated.
            pc.free_volume = precompute(&grid, &[], config.b_max, config.seed ^ 0xB0D7).free_volume;
        }

        let gpu_bodies: Vec<GpuBody> = config
            .bodies
            .iter()
            .map(|b| {
                let axis = b.orient.spin_axis();
                let (p0, p1, sigma) = b.class.exit_profile();
                let swing = config.swing_boost * b.class.chirality_sign();
                GpuBody {
                    center: [b.center.x as f32, b.center.y as f32, b.center.z as f32, 1.0],
                    axis: [axis.x as f32, axis.y as f32, axis.z as f32, b.class.chirality_sign() as f32],
                    emit: [p0 as f32, p1.map(|v| v as f32).unwrap_or(-1.0), sigma as f32, swing as f32],
                    lane: [b.through_lane_rho as f32, DsmcConfig::cos_max_angle() as f32, b.class.id() as f32, 0.0],
                }
            })
            .collect();
        let gpu_bodies_or_dummy = if gpu_bodies.is_empty() { vec![GpuBody::zeroed()] } else { gpu_bodies };

        let bodies_buf = ctx.storage_buffer_from_slice(&gpu_bodies_or_dummy, wgpu::BufferUsages::empty());
        let body_emit_counts_buf = ctx.zeroed_storage_buffer::<u32>(n_bodies.max(1), wgpu::BufferUsages::empty());

        let cap = config.capacity as usize;
        let dead = GpuPhoton::dead();
        let photons_init = vec![dead; cap];
        let photons = [
            ctx.storage_buffer_from_slice(&photons_init, wgpu::BufferUsages::empty()),
            ctx.storage_buffer_from_slice(&photons_init, wgpu::BufferUsages::empty()),
        ];

        let raw_cell_of = ctx.zeroed_storage_buffer::<u32>(cap, wgpu::BufferUsages::empty());
        let nearest_body_buf = ctx.storage_buffer_from_slice(&pc.nearest_body, wgpu::BufferUsages::empty());
        let within_bmax_buf = ctx.storage_buffer_from_slice(&pc.within_b_max, wgpu::BufferUsages::empty());
        let free_volume_buf = ctx.storage_buffer_from_slice(&pc.free_volume, wgpu::BufferUsages::empty());
        let cell_remainder_buf = ctx.zeroed_storage_buffer::<f32>(grid.num_cells as usize, wgpu::BufferUsages::empty());

        let tallies_buf = ctx.zeroed_storage_buffer::<i32>(TALLY_LEN, wgpu::BufferUsages::empty());
        let hist_buf = ctx.zeroed_storage_buffer::<u32>(HIST_LEN, wgpu::BufferUsages::empty());
        let density_buf = ctx.zeroed_storage_buffer::<u32>(grid.num_cells as usize * 4, wgpu::BufferUsages::empty());
        let free_volume_cpu = pc.free_volume.clone();

        let inject_params_buf = ctx.storage_buffer_from_slice(
            &[GpuInjectParams { n_ambient_inject: 0, n_emit_total: 0, id_base: 0, slot_base: 0 }],
            wgpu::BufferUsages::empty(),
        );
        let params_buf = ctx.uniform_buffer(&GpuParams::zeroed());
        let gather_params_buf = ctx.uniform_buffer(&NParamsLocal { n: 0, _p: [0; 3] });

        let sorter = CellSorter::new_generic(ctx, grid.num_cells + 2, config.capacity);

        let shell_bodies: Vec<crate::baryon_shells::ShellBody> = match config.body_model {
            BodyModel::Absorber => vec![],
            BodyModel::Shells { .. } | BodyModel::Ball => config
                .bodies
                .iter()
                .enumerate()
                .map(|(i, b)| {
                    let signs = config.shell_signs.get(i).copied().unwrap_or(match b.class {
                        GpuClass::Neutron => crate::baryon_shells::NEUTRON_STATES[0],
                        _ => crate::baryon_shells::PROTON_STATES[0],
                    });
                    let orient = glam::DQuat::from_rotation_arc(glam::DVec3::Z, b.orient.spin_axis());
                    crate::baryon_shells::ShellBody::new(b.center, orient, 1.0, signs)
                })
                .collect(),
        };
        let shells_buf = ctx.zeroed_storage_buffer::<[f32; 4]>(32, wgpu::BufferUsages::empty());
        let ballkeys_buf = ctx.zeroed_storage_buffer::<[f32; 4]>(4 * super::gpu_types::BALLKEY_STRIDE, wgpu::BufferUsages::empty());

        let inject_emit_kernel = load(ctx, true, include_str!("../../../shaders/dsmc_inject_emit.wgsl"));
        let move_kernel = load(ctx, true, include_str!("../../../shaders/dsmc_move.wgsl"));
        let gather_kernel = load(ctx, false, include_str!("../../../shaders/dsmc_gather.wgsl"));
        let collide_kernel = load(ctx, true, include_str!("../../../shaders/dsmc_collide.wgsl"));

        let emit_remainder = vec![0.0; n_bodies.max(1)];
        let pump_rate_ema = vec![0.0; n_bodies.max(1)];

        let domain_ctr = Grid::centroid(&config);
        let nb_centers: Vec<glam::DVec3> = config.bodies.iter().map(|b| b.center).collect();
        Self {
            config,
            grid,
            photons,
            cur: 0,
            bodies_buf,
            domain_ctr,
            nb_centers,
            body_emit_counts_buf,
            raw_cell_of,
            nearest_body_buf,
            within_bmax_buf,
            free_volume_buf,
            cell_remainder_buf,
            tallies_buf,
            hist_buf,
            density_buf,
            free_volume_cpu,
            density_accum_enabled: false,
            inject_params_buf,
            params_buf,
            gather_params_buf,
            sorter,
            held: vec![0; n_bodies.max(1)],
            shell_bodies,
            shells_buf,
            ballkeys_buf,
            inject_emit_kernel,
            move_kernel,
            gather_kernel,
            collide_kernel,
            n_used: 0,
            n_live: 0,
            step: 0,
            id_next: 0,
            ambient_remainder: 0.0,
            emit_remainder,
            pump_rate_ema,
        }
    }

    fn cur_buf(&self) -> &wgpu::Buffer {
        &self.photons[self.cur]
    }

    fn other_buf(&self) -> &wgpu::Buffer {
        &self.photons[1 - self.cur]
    }

    fn domain_center(&self) -> glam::DVec3 {
        self.domain_ctr
    }

    /// Doc §15.17: move the bodies (translation only; orientation and level
    /// phases are untouched). Updates the body buffer and the shell
    /// kinematics the ball keyframes come from; the nearest-body map is
    /// recomputed once any body has moved more than h/4 since the last
    /// time. The domain stays where it was built. The bodies' own motion is
    /// not added to the ball velocity (it is many orders below c).
    /// Doc §15.22: turn body `b` into a lepton (electron) of ball radius `r`
    /// and re-upload the bodies (its envelope becomes its ball).
    pub fn make_lepton(&mut self, ctx: &GpuContext, b: usize, r: f64) {
        let sb = &mut self.shell_bodies[b];
        sb.lepton = true;
        sb.radius = r;
        sb.signs = [1.0, 0.0, 0.0, 0.0];
        let centers: Vec<glam::DVec3> = self.config.bodies.iter().map(|c| c.center).collect();
        self.set_body_poses(ctx, &centers, None);
    }

    pub fn set_body_centers(&mut self, ctx: &GpuContext, centers: &[glam::DVec3]) {
        self.set_body_poses(ctx, centers, None);
    }

    /// Doc §15.20: move and turn the bodies. `orients` (one per body) sets
    /// each body's frame (spin axis = frame * Z); the inner levels keep
    /// their own phases inside it. Needs shell bodies (Ball/Shells model).
    pub fn set_body_poses(&mut self, ctx: &GpuContext, centers: &[glam::DVec3], orients: Option<&[glam::DQuat]>) {
        assert_eq!(centers.len(), self.config.bodies.len());
        for (i, c) in centers.iter().enumerate() {
            self.config.bodies[i].center = *c;
            if let Some(sb) = self.shell_bodies.get_mut(i) {
                sb.center = *c;
                if let Some(o) = orients {
                    sb.orient = o[i].normalize();
                }
            }
        }
        let gpu_bodies: Vec<GpuBody> = self
            .config
            .bodies
            .iter()
            .enumerate()
            .map(|(i, b)| {
                let axis = self.shell_bodies.get(i).map(|sb| sb.orient * glam::DVec3::Z).unwrap_or(b.orient.spin_axis());
                let (p0, p1, sigma) = b.class.exit_profile();
                let swing = self.config.swing_boost * b.class.chirality_sign();
                // Doc §15.22: a lepton's envelope is its own ball.
                let env_r = self.shell_bodies.get(i).filter(|sb| sb.lepton).map(|sb| sb.radius).unwrap_or(1.0);
                GpuBody {
                    center: [b.center.x as f32, b.center.y as f32, b.center.z as f32, env_r as f32],
                    axis: [axis.x as f32, axis.y as f32, axis.z as f32, b.class.chirality_sign() as f32],
                    emit: [p0 as f32, p1.map(|v| v as f32).unwrap_or(-1.0), sigma as f32, swing as f32],
                    lane: [b.through_lane_rho as f32, DsmcConfig::cos_max_angle() as f32, b.class.id() as f32, 0.0],
                }
            })
            .collect();
        ctx.write_buffer(&self.bodies_buf, &gpu_bodies);
        let moved = centers.iter().zip(&self.nb_centers).any(|(a, b)| (*a - *b).length() > self.grid.h / 4.0);
        if moved {
            let g = &self.grid;
            let mut nb = vec![0u32; g.num_cells as usize];
            for cz in 0..g.dim[2] {
                for cy in 0..g.dim[1] {
                    for cx in 0..g.dim[0] {
                        let p = g.cell_center(cx, cy, cz);
                        let (mut best, mut bd) = (0usize, f64::MAX);
                        for (i, c) in centers.iter().enumerate() {
                            let d = (p - *c).length();
                            if d < bd {
                                bd = d;
                                best = i;
                            }
                        }
                        nb[(cx + cy * g.dim[0] + cz * g.dim[0] * g.dim[1]) as usize] = best as u32;
                    }
                }
            }
            ctx.write_buffer(&self.nearest_body_buf, &nb);
            self.nb_centers = centers.to_vec();
        }
    }

    fn write_params(&self, ctx: &GpuContext) {
        let c = &self.config;
        let (boundary_mode, boundary_extent) = match c.boundary {
            Boundary::SphereEscape { r_d } => (0u32, r_d as f32),
            Boundary::PeriodicBox { half_extent } => (1u32, half_extent as f32),
        };
        let center = self.domain_center();
        let p = GpuParams {
            n_bodies: c.bodies.len() as u32,
            n_live: self.n_used,
            capacity: c.capacity,
            step: self.step,
            dt: c.dt as f32,
            sigma: c.sigma as f32,
            n_ambient: c.n_ambient as f32,
            photon_fraction: c.photon_fraction as f32,
            pole_intake_bias: c.pole_intake_bias as f32,
            cancel_redirect: c.cancel_redirect as f32,
            spin_transfer: c.spin_transfer as f32,
            b_max: c.b_max as f32,
            boundary_mode,
            boundary_extent,
            rate_kernel: if c.rate_kernel == RateKernel::Catch { 1 } else { 0 },
            mesh_mode: if c.mesh_mode == MeshMode::DirRel { 1 } else { 0 },
            bounce_cap: c.bounce_cap,
            d5_far_zone_gate: c.d5_far_zone_gate as u32,
            grid_dim_x: self.grid.dim[0],
            grid_dim_y: self.grid.dim[1],
            grid_dim_z: self.grid.dim[2],
            cell_h: self.grid.h as f32,
            grid_min_x: self.grid.min.x as f32,
            grid_min_y: self.grid.min.y as f32,
            grid_min_z: self.grid.min.z as f32,
            seed_hi: (c.seed >> 32) as u32,
            seed_lo: c.seed as u32,
            num_cells: self.grid.num_cells,
            domain_center_x: center.x as f32,
            domain_center_y: center.y as f32,
            domain_center_z: center.z as f32,
            exchange_rule: match c.exchange_rule {
                ExchangeRule::SwapRenorm => 0,
                ExchangeRule::RotateRelative => 1,
                ExchangeRule::RoughSphere => 2,
            },
            density_accumulate: self.density_accum_enabled as u32,
            emit_stiffness_level: c.emit_stiffness_level,
            stiffness_shed: if c.stiffness_shed == StiffnessShed::Transfer { 1 } else { 0 },
            sigma_size_scaling: c.sigma_size_scaling as u32,
            debug_force_equal_mass_exchange: c.debug_force_equal_mass_exchange as u32,
            g_max_catch: if c.exchange_rule == ExchangeRule::RoughSphere { c.rough.g_max as f32 } else { 2.0 },
            rough_beta: c.rough.beta as f32,
            rough_beta_twist: c.rough.beta_twist as f32,
            rough_kappa: c.rough.kappa as f32,
            twist_loss_to_normal: c.rough.twist_loss_to_normal as u32,
            speed_mode: match c.rough.speed_mode {
                SpeedMode::Free => 0,
                SpeedMode::PinnedC => 1,
                SpeedMode::BodyFree => 2,
            },
            spin_edge0: c.rough.spin_edge0 as f32,
            glancing_bias: c.rough.glancing_bias as u32,
            stacking: c.rough.stacking as u32,
            stack_max_level: c.rough.stack_max_level,
            unstack_edge: c.rough.unstack_edge as f32,
            body_model: match c.body_model {
                BodyModel::Absorber => 0,
                BodyModel::Shells { .. } => 1,
                BodyModel::Ball => 2,
            },
            tau_body: match c.body_model {
                BodyModel::Shells { tau_body } => tau_body as f32,
                _ => 0.0,
            },
            grain_edge: c.rough.grain_edge as f32,
            body_twist_loss: c.rough.body_law.map(|l| l.2).unwrap_or(c.rough.twist_loss_to_normal) as u32,
            body_beta: c.rough.body_law.map(|l| l.0).unwrap_or(c.rough.beta) as f32,
            body_beta_twist: c.rough.body_law.map(|l| l.1).unwrap_or(c.rough.beta_twist) as f32,
            grain_mode: match c.rough.grain_mode {
                super::config::GrainMode::EdgeC => 0,
                super::config::GrainMode::Rigid => 1,
                super::config::GrainMode::RigidInner => 2,
            },
            maze_mode: match c.rough.maze {
                super::config::MazeMode::Off => 0,
                super::config::MazeMode::Fixed => 1,
                super::config::MazeMode::Swept => 2,
                super::config::MazeMode::Isotropic => 3,
            } | match c.rough.maze_speed {
                super::config::MazeSpeed::EjectC => 256,
                super::config::MazeSpeed::EjectRest => 256 | 512,
                super::config::MazeSpeed::EjectSurface => 256 | 1024,
                super::config::MazeSpeed::Relative => 0,
            },
        };
        ctx.write_buffer(&self.params_buf, &[p]);
        self.upload_shells(ctx);
    }

    /// Doc §14: uploads every shell body's current level frames (world
    /// space) to the move kernel.
    fn upload_shells(&self, ctx: &GpuContext) {
        if self.shell_bodies.is_empty() {
            return;
        }
        let mut data = vec![[0.0f32; 4]; 32];
        for (b, body) in self.shell_bodies.iter().enumerate() {
            for (k, f) in body.frames().iter().enumerate() {
                data[b * 8 + 2 * k] = [f.origin.x as f32, f.origin.y as f32, f.origin.z as f32, f.radius as f32];
                data[b * 8 + 2 * k + 1] = [f.omega.x as f32, f.omega.y as f32, f.omega.z as f32, 0.0];
            }
        }
        ctx.write_buffer(&self.shells_buf, &data);
        if self.config.body_model == BodyModel::Ball {
            ctx.write_buffer(&self.ballkeys_buf, &self.ball_keyframes());
        }
    }

    /// Doc §15/§15.8: 9 keyframes per body over this step (s * dt/8) plus the
    /// unresolved inner levels; layout in gpu_types::BALLKEY_STRIDE.
    pub fn ball_keyframes(&self) -> Vec<[f32; 4]> {
        use super::gpu_types::BALLKEY_STRIDE;
        let mut data = vec![[0.0f32; 4]; 4 * BALLKEY_STRIDE];
        for (b, body) in self.shell_bodies.iter().enumerate() {
            let mut sb = body.clone();
            let mut prev_q: Option<glam::DQuat> = None;
            for s in 0..9 {
                let st = sb.ball();
                let base = b * BALLKEY_STRIDE + s * 4;
                data[base] = [st.center.x as f32, st.center.y as f32, st.center.z as f32, st.radius as f32];
                data[base + 1] = [st.omega.x as f32, st.omega.y as f32, st.omega.z as f32, 0.0];
                data[base + 2] = [st.velocity.x as f32, st.velocity.y as f32, st.velocity.z as f32, 0.0];
                // Keep consecutive quaternions in the same hemisphere so the
                // shader's normalised lerp takes the short way.
                let mut q = st.frame;
                if let Some(pq) = prev_q {
                    if pq.dot(q) < 0.0 {
                        q = -q;
                    }
                }
                prev_q = Some(q);
                data[base + 3] = [q.x as f32, q.y as f32, q.z as f32, q.w as f32];
                sb.advance(self.config.dt / 8.0);
            }
            let lv = body.unresolved_levels();
            let th9 = body.angles[0].rem_euclid(std::f64::consts::TAU);
            data[b * BALLKEY_STRIDE + 36] = [lv.len() as f32, body.photon_radius() as f32, th9 as f32, body.rate(0) as f32];
            let o = body.orient;
            data[b * BALLKEY_STRIDE + 45] = [o.x as f32, o.y as f32, o.z as f32, o.w as f32];
            data[b * BALLKEY_STRIDE + 46] = [body.signs[0] as f32, body.signs[1] as f32, body.signs[2] as f32, body.signs[3] as f32];
            let r = &self.config.rough;
            let neutron = self.config.bodies.get(b).map(|c| c.class == super::config::GpuClass::Neutron).unwrap_or(false);
            let hold = r.hold_time * if neutron { r.hold_time_neutron } else { 1.0 };
            // Doc §15.22: lepton bodies never take photons in.
            let open = !body.lepton && (r.hold_cap == 0 || self.held.get(b).copied().unwrap_or(0) < r.hold_cap as i64);
            let (mode, par) = match r.intake {
                super::config::Intake::Ball => (0.0, 0.0),
                super::config::Intake::Envelope { p } => (1.0, p),
                super::config::Intake::Polar { aperture } => (2.0, aperture),
            };
            data[b * BALLKEY_STRIDE + 47] = [if open { 1.0 } else { 0.0 }, (hold / self.config.dt).round() as f32, par as f32, mode];
            data[b * BALLKEY_STRIDE + 48] = [r.through_hole as f32, 0.0, 0.0, 0.0];
            for (i, (axis, rate)) in lv.iter().enumerate() {
                data[b * BALLKEY_STRIDE + 37 + i] = [axis.x as f32, axis.y as f32, axis.z as f32, *rate as f32];
            }
        }
        data
    }

    fn advance_shells(&mut self) {
        let dt = self.config.dt;
        for body in &mut self.shell_bodies {
            body.advance(dt);
        }
    }

    /// Injects `n_ambient` ambient photons and `emit_counts[i]` emitted
    /// photons for body `i`, appended right after the current live range.
    /// One dispatch; the caller decides how often to call this (a single
    /// burst for a transport-only gate, or every step for a steady-state
    /// loop).
    pub fn inject_batch(&mut self, ctx: &GpuContext, n_ambient: u32, emit_counts: &[u32]) {
        let n_emit_total: u32 = emit_counts.iter().sum();
        let total = n_ambient + n_emit_total;
        if total == 0 {
            return;
        }
        assert!(
            self.n_used + total <= self.config.capacity,
            "DSMC capacity exceeded: {} + {} > {}",
            self.n_used,
            total,
            self.config.capacity
        );

        ctx.write_buffer(&self.body_emit_counts_buf, emit_counts);
        ctx.write_buffer(
            &self.inject_params_buf,
            &[GpuInjectParams { n_ambient_inject: n_ambient, n_emit_total, id_base: self.id_next, slot_base: self.n_used }],
        );
        self.write_params(ctx);

        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        let grid = workgroup_grid_for_threads(total, 256);
        self.inject_emit_kernel.record(
            ctx,
            &mut encoder,
            &[self.cur_buf(), &self.bodies_buf, &self.body_emit_counts_buf, &self.tallies_buf, &self.hist_buf, &self.inject_params_buf],
            Some(&self.params_buf),
            grid,
        );
        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll failed");

        self.n_used += total;
        self.n_live += total;
        self.id_next += total;
    }

    /// Move-only step (no sort/gather/collide): valid whenever sigma == 0,
    /// since independent photon trajectories need no binning at all (G0).
    pub fn step_transport_only(&mut self, ctx: &GpuContext, totals: &mut Totals) {
        self.write_params(ctx);
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        let grid = workgroup_grid_for_threads(self.n_used, 256);
        self.move_kernel.record(
            ctx,
            &mut encoder,
            &[self.cur_buf(), &self.bodies_buf, &self.raw_cell_of, &self.tallies_buf, &self.hist_buf, &self.density_buf, &self.shells_buf, &self.ballkeys_buf],
            Some(&self.params_buf),
            grid,
        );
        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll failed");
        self.readback_tallies(ctx, totals, self.n_used);
        self.step += 1;
        self.advance_shells();
    }

    /// Full move -> sort -> gather -> collide step (G1, and the future
    /// steady-state loop). Compacts `n_used` down to the post-gather live
    /// count so dead slots are reclaimed for the next injection.
    pub fn step_collide(&mut self, ctx: &GpuContext, totals: &mut Totals) {
        self.write_params(ctx);
        let n = self.n_used;

        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        let grid_n = workgroup_grid_for_threads(n, 256);
        self.move_kernel.record(
            ctx,
            &mut encoder,
            &[self.cur_buf(), &self.bodies_buf, &self.raw_cell_of, &self.tallies_buf, &self.hist_buf, &self.density_buf, &self.shells_buf, &self.ballkeys_buf],
            Some(&self.params_buf),
            grid_n,
        );

        // Sort + gather ALWAYS run: besides binning for collisions they are
        // what compacts dead photons out of the slot range. Skipping them at
        // sigma = 0 leaked every dead slot (a flux-driven sigma = 0 run filled
        // any capacity in ~1000 steps). Only the collide pass is optional.
        {
            self.sorter.record_from_cell_of(ctx, &mut encoder, &self.raw_cell_of, n);
            ctx.write_buffer(&self.gather_params_buf, &[NParamsLocal { n, _p: [0; 3] }]);
            self.gather_kernel.record(
                ctx,
                &mut encoder,
                &[&self.sorter.sorted, self.cur_buf(), self.other_buf()],
                Some(&self.gather_params_buf),
                grid_n,
            );
        }
        if self.config.sigma > 0.0 {
            // One workgroup PER CELL (not 64 cells/workgroup as before): the
            // collide kernel now sorts and pairs the whole cell in workgroup
            // shared memory, so `local_size=1` here just turns
            // `workgroup_grid_for_threads` into "however many 2-D workgroups
            // it takes to cover num_cells", matching dsmc_collide.wgsl's
            // `cell = wid.x + wid.y * nwg.x` indexing.
            let grid_cells = workgroup_grid_for_threads(self.grid.num_cells, 1);
            self.collide_kernel.record(
                ctx,
                &mut encoder,
                &[
                    self.other_buf(),
                    &self.sorter.offsets,
                    &self.sorter.counts,
                    &self.bodies_buf,
                    &self.nearest_body_buf,
                    &self.within_bmax_buf,
                    &self.free_volume_buf,
                    &self.cell_remainder_buf,
                    &self.tallies_buf,
                ],
                Some(&self.params_buf),
                grid_cells,
            );
        }

        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll failed");

        {
            self.cur = 1 - self.cur;
            // n_live_after = offsets[grid.num_cells + 1] (the dead bucket's
            // exclusive-scan offset; bucket num_cells parks held photons, doc §15.11).
            let offsets: Vec<u32> = ctx.read_buffer(&self.sorter.offsets, self.grid.num_cells as usize + 2);
            let n_live_after = offsets[self.grid.num_cells as usize + 1];
            self.n_live = n_live_after;
            self.n_used = n_live_after;
        }

        self.readback_tallies(ctx, totals, n);
        self.step += 1;
        self.advance_shells();
    }

    fn readback_tallies(&mut self, ctx: &GpuContext, totals: &mut Totals, sample_count: u32) {
        let raw: Vec<i32> = ctx.read_buffer(&self.tallies_buf, TALLY_LEN);
        // Doc §15.11: photons held in each body's tank (for the intake cap).
        for b in 0..self.held.len() {
            use super::gpu_types::TALLY_HOLD;
            self.held[b] += raw[TALLY_HOLD + b * 2] as i64 - raw[TALLY_HOLD + b * 2 + 1] as i64;
        }

        // Pump mode (doc §3.3): each body's next emission rate is an EMA of
        // its OWN absorption over PUMP_EMA_STEPS steps, bootstrapped from
        // zero (pump_rate_ema starts at vec![0.0; ..] in `new`). Folded from
        // the raw per-step tally directly (not from `totals`, which the
        // caller may be accumulating cumulatively across many steps) so this
        // works identically whether the caller passes a fresh `Totals` per
        // step or one shared accumulator.
        if self.config.emission_mode == EmissionMode::Pump {
            let alpha = 1.0 / super::config::PUMP_EMA_STEPS;
            for body in 0..self.config.bodies.len().min(4) {
                let absorbed_this_step: i64 = (0..2)
                    .flat_map(|k| (0..2).map(move |s| (k, s)))
                    .map(|(k, s)| raw[body * 4 + k * 2 + s] as i64)
                    .sum();
                let ema = &mut self.pump_rate_ema[body];
                *ema += alpha * (absorbed_this_step as f64 - *ema);
            }
        }

        totals.fold(&raw, sample_count);
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.clear_buffer(&self.tallies_buf, 0, None);
        ctx.queue.submit(Some(encoder.finish()));
    }

    pub fn read_hist(&self, ctx: &GpuContext) -> Vec<u32> {
        ctx.read_buffer(&self.hist_buf, HIST_LEN)
    }

    /// Starts the density-map averaging window (doc §8): clears the per-cell
    /// [kind x species] u32 accumulator and tells `dsmc_move.wgsl` to start
    /// adding to it every step. Call once, right after warm-up, before the
    /// averaging loop -- never during warm-up (a still-transient population
    /// would bias the map).
    pub fn begin_density_accum(&mut self, ctx: &GpuContext) {
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.clear_buffer(&self.density_buf, 0, None);
        ctx.queue.submit(Some(encoder.finish()));
        self.density_accum_enabled = true;
    }

    pub fn end_density_accum(&mut self) {
        self.density_accum_enabled = false;
    }

    /// Reads back the raw per-cell [kind x species] u32 accumulator (len
    /// `grid.num_cells * 4`, index `cell*4 + kind*2 + species`). Wrap in
    /// `density::DensityMap::new` for number-density queries, ray sampling
    /// and the xz-slice CSV dump.
    pub fn read_density_raw(&self, ctx: &GpuContext) -> Vec<u32> {
        ctx.read_buffer(&self.density_buf, self.grid.num_cells as usize * 4)
    }

    pub fn free_volume_cpu(&self) -> &[f32] {
        &self.free_volume_cpu
    }

    /// Reads back the first `n_used` photon records from the live buffer.
    /// Test/report use only (e.g. G1-dense's per-photon collision-fairness
    /// check) -- not on any hot path.
    pub fn read_photons(&self, ctx: &GpuContext) -> Vec<GpuPhoton> {
        ctx.read_buffer(self.cur_buf(), self.n_used as usize)
    }

    /// Test/gate-only: overwrites every LIVE photon's stiffness level (E1,
    /// doc §11d `c[3]`) with `f(index_among_live)`, leaving everything else
    /// (position, direction, flags) untouched. Used to build synthetic mixed-
    /// stiffness populations (e.g. "half k=0, half k=3") that no ordinary
    /// emission path produces, so the E1 gates can test the exchange/shedding
    /// math directly instead of waiting on a steady-state run to populate it.
    pub fn set_stiffness_levels_for_test(&mut self, ctx: &GpuContext, f: impl Fn(usize) -> u32) {
        let mut ph = self.read_photons(ctx);
        for (i, p) in ph.iter_mut().enumerate() {
            if p.c[0] & super::gpu_types::FLAG_ALIVE != 0 {
                p.c[3] = f(i);
            }
        }
        ctx.write_buffer(self.cur_buf(), &ph);
    }

    /// Test/gate-only: replaces the whole photon population with `ph` (ids,
    /// flags and all), as if it had been injected. Used by the doc §13 gates
    /// to set up exact two-photon collisions and synthetic boxes.
    pub fn write_photons_for_test(&mut self, ctx: &GpuContext, ph: &[GpuPhoton]) {
        assert!(ph.len() as u32 <= self.config.capacity);
        ctx.write_buffer(self.cur_buf(), ph);
        self.n_used = ph.len() as u32;
        self.n_live = self.n_used;
        self.id_next = ph.iter().map(|p| p.id() + 1).max().unwrap_or(0).max(self.id_next);
    }

    /// Doc §3.3 flux formulas, exposed so gate/report code doesn't duplicate
    /// them: ambient boundary flux and the Forced-mode per-body emission
    /// rate (rate of lines with b <= B_MAX).
    pub fn ambient_flux_per_step(&self) -> f64 {
        let r_d = match self.config.boundary {
            Boundary::SphereEscape { r_d } => r_d,
            Boundary::PeriodicBox { .. } => return 0.0,
        };
        self.config.n_ambient * 1.0 * std::f64::consts::PI * r_d * r_d * self.config.dt
    }

    pub fn forced_emission_rate_per_step(&self) -> f64 {
        self.config.n_ambient * std::f64::consts::PI * self.config.b_max * self.config.b_max * self.config.dt
    }

    /// Advances the host-side fractional injection counters by one step and
    /// returns integer counts to inject this step. Forced mode uses the
    /// fixed Phase-1-analogue rate; Pump mode uses `pump_rate_ema[i]`, which
    /// `readback_tallies` updates from each step's own absorption tally (so
    /// this call must follow the PREVIOUS step's `step_collide`/
    /// `step_transport_only` for the EMA to reflect it).
    pub fn next_injection_counts(&mut self) -> (u32, Vec<u32>) {
        self.ambient_remainder += self.ambient_flux_per_step();
        let n_ambient = self.ambient_remainder.floor().max(0.0) as u32;
        self.ambient_remainder -= n_ambient as f64;

        let mut counts = Vec::with_capacity(self.config.bodies.len());
        for (i, _b) in self.config.bodies.iter().enumerate() {
            let rate = match (self.config.body_model, self.config.emission_mode) {
                // Doc §14: shell bodies emit nothing; their "emission" is
                // scattered ambient.
                (BodyModel::Shells { .. } | BodyModel::Ball, _) => 0.0,
                (_, EmissionMode::Forced) => self.forced_emission_rate_per_step(),
                (_, EmissionMode::Pump) => self.pump_rate_ema.get(i).copied().unwrap_or(0.0),
            };
            self.emit_remainder[i] += rate;
            let n = self.emit_remainder[i].floor().max(0.0) as u32;
            self.emit_remainder[i] -= n as f64;
            counts.push(n);
        }
        (n_ambient, counts)
    }

    /// One flux-driven step: advance the host injection counters, inject,
    /// then run the full collide pipeline. Used by both `run_steady` and any
    /// caller that wants the doc §3.3 continuous-flow loop without the
    /// block/steady bookkeeping. Returns this step's OWN tallies (not
    /// cumulative), so conservation and pump-EMA checks can be made per step.
    pub fn step_flux(&mut self, ctx: &GpuContext) -> Totals {
        let (n_ambient, emit_counts) = self.next_injection_counts();
        self.inject_batch(ctx, n_ambient, &emit_counts);
        let mut step_totals = Totals::default();
        self.step_collide(ctx, &mut step_totals);
        step_totals
    }

    /// One reporting block of the steady-state time series (doc B2b work
    /// item 1): live count and per-body absorption/emission rate at the END
    /// of the block, averaged over the block's steps, plus the collision
    /// rate. Rates are per unit time (steps * dt), so blocks of different
    /// length are still comparable.
    pub fn run_steady(&mut self, ctx: &GpuContext, warmup_steps: u32, avg_steps: u32, block: u32) -> SteadyReport {
        assert!(block > 0, "run_steady: block must be > 0");
        let n_bodies = self.config.bodies.len();

        for _ in 0..warmup_steps {
            self.step_flux(ctx);
        }

        // Density map (doc §8) accumulates over the averaging window only --
        // never warm-up, which is still transient by construction.
        self.begin_density_accum(ctx);

        let mut totals = Totals::default();
        let mut blocks: Vec<SteadyBlock> = Vec::new();
        let mut block_totals = Totals::default();
        let mut steps_in_block: u32 = 0;

        let mut steps_done = 0u32;
        while steps_done < avg_steps {
            let step_totals = self.step_flux(ctx);
            totals.add(&step_totals);
            block_totals.add(&step_totals);
            steps_in_block += 1;
            steps_done += 1;

            if steps_in_block == block || steps_done == avg_steps {
                let dt_block = steps_in_block as f64 * self.config.dt;
                let absorbed_rate: Vec<f64> = (0..n_bodies).map(|b| block_totals.absorbed_total(b) as f64 / dt_block).collect();
                let emitted_rate: Vec<f64> = (0..n_bodies).map(|b| block_totals.emitted[b] as f64 / dt_block).collect();
                blocks.push(SteadyBlock {
                    step: self.step,
                    live: self.n_live,
                    absorbed_total: (0..n_bodies).map(|b| block_totals.absorbed_total(b)).sum(),
                    absorbed_rate,
                    emitted_rate,
                    collisions_rate: block_totals.collisions as f64 / dt_block,
                });
                block_totals = Totals::default();
                steps_in_block = 0;
            }
        }

        let density_raw = self.read_density_raw(ctx);
        self.end_density_accum();

        let (steady_ok, steady_detail) = steady_state_check(&blocks);

        SteadyReport { totals, blocks, warmup_steps, avg_steps, block, steady_ok, steady_detail, density_raw }
    }
}

/// One block of the steady-state time series (doc §1: "records a block time
/// series ... so steady state can be judged").
#[derive(Clone, Debug)]
pub struct SteadyBlock {
    pub step: u32,
    pub live: u32,
    pub absorbed_total: i64,
    pub absorbed_rate: Vec<f64>,
    pub emitted_rate: Vec<f64>,
    pub collisions_rate: f64,
}

/// Result of `DsmcEngine::run_steady` (doc B2b work item 1).
#[derive(Clone, Debug)]
pub struct SteadyReport {
    /// Accumulated over the averaging window only (warm-up discarded).
    pub totals: Totals,
    pub blocks: Vec<SteadyBlock>,
    pub warmup_steps: u32,
    pub avg_steps: u32,
    pub block: u32,
    pub steady_ok: bool,
    pub steady_detail: String,
    /// Raw per-cell [kind x species] density accumulator (doc §8), covering
    /// exactly the averaging window. Wrap with `density::DensityMap::new`
    /// (grid geometry + `DsmcEngine::free_volume_cpu()` + `avg_steps`) for
    /// number-density queries.
    pub density_raw: Vec<u32>,
}

/// Doc's steady-state check: "last-half vs first-half block means of live
/// count and absorption within 3σ of block-to-block scatter" -- a two-sample
/// z-test using the SAMPLE standard deviation of the per-block values as the
/// scatter estimate (not a Poisson assumption: a block sums many correlated
/// per-photon events, so the block-to-block spread is the honest error bar).
/// Needs >= 4 blocks (2 per half) to form a scatter estimate at all; fewer is
/// reported as "insufficient blocks" and never silently passed.
fn steady_state_check(blocks: &[SteadyBlock]) -> (bool, String) {
    if blocks.len() < 4 {
        return (false, format!("insufficient blocks ({}) for a first-half/second-half steady check (need >= 4)", blocks.len()));
    }
    let half = blocks.len() / 2;
    let live: Vec<f64> = blocks.iter().map(|b| b.live as f64).collect();
    let absorbed: Vec<f64> = blocks.iter().map(|b| b.absorbed_total as f64).collect();

    fn half_z(values: &[f64], half: usize) -> f64 {
        let (a, b) = values.split_at(half);
        let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len() as f64;
        let var = |v: &[f64], m: f64| v.iter().map(|x| (x - m).powi(2)).sum::<f64>() / (v.len().max(2) - 1) as f64;
        let ma = mean(a);
        let mb = mean(b);
        let pooled_var = (var(a, ma) * (a.len() - 1) as f64 + var(b, mb) * (b.len() - 1) as f64) / ((a.len() + b.len() - 2).max(1) as f64);
        let se = (pooled_var * (1.0 / a.len() as f64 + 1.0 / b.len() as f64)).sqrt();
        if se > 1e-12 {
            (mb - ma) / se
        } else {
            0.0
        }
    }

    let z_live = half_z(&live, half);
    let z_absorbed = half_z(&absorbed, half);
    let ok = z_live.abs() < 3.0 && z_absorbed.abs() < 3.0;
    let detail = format!(
        "{} blocks (half={}): z_live={z_live:.2} z_absorbed={z_absorbed:.2} (limit 3.0, first-half vs second-half block means over block-to-block scatter)",
        blocks.len(),
        half
    );
    (ok, detail)
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct NParamsLocal {
    n: u32,
    _p: [u32; 3],
}

#[allow(dead_code)]
fn _assert_orient_used(_: GpuOrient) {}
