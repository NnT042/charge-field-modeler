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

use super::config::{Boundary, DsmcConfig, EmissionMode, GpuOrient, MeshMode, RateKernel};
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
}

impl Totals {
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
        self.dir_sum_samples += sample_count as i64;
    }
}

pub struct DsmcEngine {
    pub config: DsmcConfig,
    pub grid: Grid,
    photons: [wgpu::Buffer; 2],
    cur: usize,
    bodies_buf: wgpu::Buffer,
    body_emit_counts_buf: wgpu::Buffer,
    raw_cell_of: wgpu::Buffer,
    nearest_body_buf: wgpu::Buffer,
    within_bmax_buf: wgpu::Buffer,
    free_volume_buf: wgpu::Buffer,
    cell_remainder_buf: wgpu::Buffer,
    tallies_buf: wgpu::Buffer,
    hist_buf: wgpu::Buffer,
    inject_params_buf: wgpu::Buffer,
    params_buf: wgpu::Buffer,
    gather_params_buf: wgpu::Buffer,
    sorter: CellSorter,

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
        let grid = Grid::new(&config);
        let n_bodies = config.bodies.len();
        let pc = precompute(&grid, &config.bodies, config.b_max, config.seed ^ 0xB0D7);

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

        let inject_params_buf = ctx.storage_buffer_from_slice(
            &[GpuInjectParams { n_ambient_inject: 0, n_emit_total: 0, id_base: 0, slot_base: 0 }],
            wgpu::BufferUsages::empty(),
        );
        let params_buf = ctx.uniform_buffer(&GpuParams::zeroed());
        let gather_params_buf = ctx.uniform_buffer(&NParamsLocal { n: 0, _p: [0; 3] });

        let sorter = CellSorter::new_generic(ctx, grid.num_cells + 1, config.capacity);

        let inject_emit_kernel = load(ctx, true, include_str!("../../../shaders/dsmc_inject_emit.wgsl"));
        let move_kernel = load(ctx, true, include_str!("../../../shaders/dsmc_move.wgsl"));
        let gather_kernel = load(ctx, false, include_str!("../../../shaders/dsmc_gather.wgsl"));
        let collide_kernel = load(ctx, true, include_str!("../../../shaders/dsmc_collide.wgsl"));

        let emit_remainder = vec![0.0; n_bodies.max(1)];
        let pump_rate_ema = vec![0.0; n_bodies.max(1)];

        Self {
            config,
            grid,
            photons,
            cur: 0,
            bodies_buf,
            body_emit_counts_buf,
            raw_cell_of,
            nearest_body_buf,
            within_bmax_buf,
            free_volume_buf,
            cell_remainder_buf,
            tallies_buf,
            hist_buf,
            inject_params_buf,
            params_buf,
            gather_params_buf,
            sorter,
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
        Grid::centroid(&self.config)
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
            _pad2: 0.0,
        };
        ctx.write_buffer(&self.params_buf, &[p]);
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
            &[self.cur_buf(), &self.bodies_buf, &self.raw_cell_of, &self.tallies_buf, &self.hist_buf],
            Some(&self.params_buf),
            grid,
        );
        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::PollType::wait_indefinitely()).expect("poll failed");
        self.readback_tallies(ctx, totals, self.n_used);
        self.step += 1;
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
            &[self.cur_buf(), &self.bodies_buf, &self.raw_cell_of, &self.tallies_buf, &self.hist_buf],
            Some(&self.params_buf),
            grid_n,
        );

        if self.config.sigma > 0.0 {
            self.sorter.record_from_cell_of(ctx, &mut encoder, &self.raw_cell_of, n);
            ctx.write_buffer(&self.gather_params_buf, &[NParamsLocal { n, _p: [0; 3] }]);
            self.gather_kernel.record(
                ctx,
                &mut encoder,
                &[&self.sorter.sorted, self.cur_buf(), self.other_buf()],
                Some(&self.gather_params_buf),
                grid_n,
            );
            let grid_cells = workgroup_grid_for_threads(self.grid.num_cells, 64);
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

        if self.config.sigma > 0.0 {
            self.cur = 1 - self.cur;
            // n_live_after = offsets[grid.num_cells] (the dead bucket's exclusive-scan offset).
            let offsets: Vec<u32> = ctx.read_buffer(&self.sorter.offsets, self.grid.num_cells as usize + 1);
            let n_live_after = offsets[self.grid.num_cells as usize];
            self.n_live = n_live_after;
            self.n_used = n_live_after;
        }

        self.readback_tallies(ctx, totals, n);
        self.step += 1;
    }

    fn readback_tallies(&self, ctx: &GpuContext, totals: &mut Totals, sample_count: u32) {
        let raw: Vec<i32> = ctx.read_buffer(&self.tallies_buf, TALLY_LEN);
        totals.fold(&raw, sample_count);
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.clear_buffer(&self.tallies_buf, 0, None);
        ctx.queue.submit(Some(encoder.finish()));
    }

    pub fn read_hist(&self, ctx: &GpuContext) -> Vec<u32> {
        ctx.read_buffer(&self.hist_buf, HIST_LEN)
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
    /// returns integer counts to inject this step (Forced emission mode;
    /// Pump mode is a documented no-op placeholder for B2/B3 -- neither B1
    /// gate uses it).
    pub fn next_injection_counts(&mut self) -> (u32, Vec<u32>) {
        self.ambient_remainder += self.ambient_flux_per_step();
        let n_ambient = self.ambient_remainder.floor().max(0.0) as u32;
        self.ambient_remainder -= n_ambient as f64;

        let mut counts = Vec::with_capacity(self.config.bodies.len());
        for (i, _b) in self.config.bodies.iter().enumerate() {
            let rate = match self.config.emission_mode {
                EmissionMode::Forced => self.forced_emission_rate_per_step(),
                EmissionMode::Pump => self.pump_rate_ema.get(i).copied().unwrap_or(0.0),
            };
            self.emit_remainder[i] += rate;
            let n = self.emit_remainder[i].floor().max(0.0) as u32;
            self.emit_remainder[i] -= n as f64;
            counts.push(n);
        }
        (n_ambient, counts)
    }
}

#[repr(C)]
#[derive(Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
struct NParamsLocal {
    n: u32,
    _p: [u32; 3],
}

#[allow(dead_code)]
fn _assert_orient_used(_: GpuOrient) {}
