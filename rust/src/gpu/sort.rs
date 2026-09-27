//! GPU-resident counting sort of particles by grid cell — the DSMC binning
//! primitive. After `record`, `offsets[c]..offsets[c]+counts[c]` indexes the
//! slice of `sorted` holding cell `c`'s particle indices.
//!
//! All buffers are allocated once for a maximum particle count; `record`
//! encodes every pass (clear, count, scan, scatter) into a caller-supplied
//! encoder, so a whole sort is one submit with no host round trips. Order
//! within a cell is not deterministic (atomic cursors) — DSMC pair selection
//! must not depend on it.

use bytemuck::{Pod, Zeroable};

use crate::gpu::context::GpuContext;
use crate::gpu::kernel::{workgroup_grid_for_threads, Kernel};

const BLOCK_SIZE: u32 = 256;

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct CountParams {
    n: u32,
    grid_dim: u32,
    _pad0: u32,
    _pad1: u32,
}

#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable)]
struct NParams {
    n: u32,
    _pad0: u32,
    _pad1: u32,
    _pad2: u32,
}

/// One level of the multi-level exclusive scan: scans `n` elements of the
/// level's input into per-block sums held in `sums`.
struct ScanLevel {
    n: u32,
    sums: wgpu::Buffer,
    params: wgpu::Buffer,
}

pub struct CellSorter {
    cell_count: Kernel,
    cell_count_from_index: Kernel,
    scan_block: Kernel,
    add_block_sums: Kernel,
    cell_scatter: Kernel,
    pub grid_dim: u32,
    pub num_cells: u32,
    pub n_max: u32,
    pub counts: wgpu::Buffer,
    pub offsets: wgpu::Buffer,
    cursor: wgpu::Buffer,
    pub cell_of: wgpu::Buffer,
    pub sorted: wgpu::Buffer,
    count_params: wgpu::Buffer,
    count_from_index_params: wgpu::Buffer,
    scatter_params: wgpu::Buffer,
    levels: Vec<ScanLevel>,
}

impl CellSorter {
    pub fn new(ctx: &GpuContext, grid_dim: u32, n_max: u32) -> Self {
        Self::new_impl(ctx, grid_dim, grid_dim.pow(3), n_max)
    }

    /// A sorter over an arbitrary cell count with no position-based path --
    /// for callers (like the DSMC engine) whose own move kernel computes
    /// `cell_of` (bodies, lane gates, boundary) and just needs binning +
    /// counting-sort from that precomputed buffer via
    /// [`Self::record_from_cell_of`]. `num_cells` here is whatever the caller
    /// wants, e.g. a spatial grid's cell count *plus one* dead-photon bucket.
    pub fn new_generic(ctx: &GpuContext, num_cells: u32, n_max: u32) -> Self {
        Self::new_impl(ctx, 0, num_cells, n_max)
    }

    fn new_impl(ctx: &GpuContext, grid_dim: u32, num_cells: u32, n_max: u32) -> Self {
        let load = |src: &str| {
            let module = ctx.load_wgsl_compute(src);
            Kernel::new(ctx, &module, "main")
        };
        // WGSL, not GLSL: naga's GLSL frontend has no atomics, and its GLSL
        // `barrier()` lowers to a subgroup-scoped barrier (see scan_block.wgsl).
        let cell_count = load(include_str!("../../shaders/cell_count.wgsl"));
        let cell_count_from_index = load(include_str!("../../shaders/cell_count_from_index.wgsl"));
        let scan_block = load(include_str!("../../shaders/scan_block.wgsl"));
        let add_block_sums = load(include_str!("../../shaders/add_block_sums.wgsl"));
        let cell_scatter = load(include_str!("../../shaders/cell_scatter.wgsl"));

        let buf_u32 = |count: u32| ctx.zeroed_storage_buffer::<u32>(count.max(1) as usize, wgpu::BufferUsages::empty());

        // Scan levels: level 0 scans `offsets` (num_cells); each further level
        // scans the previous level's block sums until one block remains.
        let mut levels = Vec::new();
        let mut n = num_cells;
        loop {
            let num_blocks = n.div_ceil(BLOCK_SIZE);
            levels.push(ScanLevel {
                n,
                sums: buf_u32(num_blocks),
                params: ctx.uniform_buffer(&NParams { n, _pad0: 0, _pad1: 0, _pad2: 0 }),
            });
            if num_blocks <= 1 {
                break;
            }
            n = num_blocks;
        }

        Self {
            cell_count,
            cell_count_from_index,
            scan_block,
            add_block_sums,
            cell_scatter,
            grid_dim,
            num_cells,
            n_max,
            counts: buf_u32(num_cells),
            offsets: buf_u32(num_cells),
            cursor: buf_u32(num_cells),
            cell_of: buf_u32(n_max),
            sorted: buf_u32(n_max),
            count_params: ctx.uniform_buffer(&CountParams { n: n_max, grid_dim, _pad0: 0, _pad1: 0 }),
            count_from_index_params: ctx.uniform_buffer(&NParams { n: n_max, _pad0: 0, _pad1: 0, _pad2: 0 }),
            scatter_params: ctx.uniform_buffer(&NParams { n: n_max, _pad0: 0, _pad1: 0, _pad2: 0 }),
            levels,
        }
    }

    /// Encode a full sort of the first `n` particles of `positions`
    /// (vec4<f32> per particle, xyz in cell units) into `encoder`.
    pub fn record(&self, ctx: &GpuContext, encoder: &mut wgpu::CommandEncoder, positions: &wgpu::Buffer, n: u32) {
        assert!(n <= self.n_max, "CellSorter sized for {} particles, got {n}", self.n_max);
        ctx.write_buffer(&self.count_params, &[CountParams { n, grid_dim: self.grid_dim, _pad0: 0, _pad1: 0 }]);
        ctx.write_buffer(&self.scatter_params, &[NParams { n, _pad0: 0, _pad1: 0, _pad2: 0 }]);

        let cells_bytes = self.num_cells as u64 * 4;
        let grid_n = workgroup_grid_for_threads(n, BLOCK_SIZE);

        encoder.clear_buffer(&self.counts, 0, None);
        self.cell_count
            .record(ctx, encoder, &[positions, &self.counts, &self.cell_of], Some(&self.count_params), grid_n);

        self.scan_and_scatter(ctx, encoder, n, cells_bytes, grid_n);
    }

    /// Encode a full sort of the first `n` photons using an already-computed
    /// `cell_of_src` buffer (e.g. the DSMC move kernel's output, which folds
    /// in body/lane/boundary logic no generic position-based pass could know
    /// about). `cell_of_src[i] == self.num_cells - 1`'s bucket is the
    /// convention callers use for "dead" -- this sorter does not interpret
    /// it, it just bins whatever indices it is given.
    pub fn record_from_cell_of(
        &self,
        ctx: &GpuContext,
        encoder: &mut wgpu::CommandEncoder,
        cell_of_src: &wgpu::Buffer,
        n: u32,
    ) {
        assert!(n <= self.n_max, "CellSorter sized for {} particles, got {n}", self.n_max);
        ctx.write_buffer(&self.count_from_index_params, &[NParams { n, _pad0: 0, _pad1: 0, _pad2: 0 }]);
        ctx.write_buffer(&self.scatter_params, &[NParams { n, _pad0: 0, _pad1: 0, _pad2: 0 }]);

        let cells_bytes = self.num_cells as u64 * 4;
        let grid_n = workgroup_grid_for_threads(n, BLOCK_SIZE);

        encoder.clear_buffer(&self.counts, 0, None);
        self.cell_count_from_index.record(
            ctx,
            encoder,
            &[cell_of_src, &self.counts, &self.cell_of],
            Some(&self.count_from_index_params),
            grid_n,
        );

        self.scan_and_scatter(ctx, encoder, n, cells_bytes, grid_n);
    }

    fn scan_and_scatter(
        &self,
        ctx: &GpuContext,
        encoder: &mut wgpu::CommandEncoder,
        n: u32,
        cells_bytes: u64,
        grid_n: (u32, u32, u32),
    ) {
        let _ = n;
        encoder.copy_buffer_to_buffer(&self.counts, 0, &self.offsets, 0, cells_bytes);
        // Down-sweep: scan every level into its block sums...
        for (i, level) in self.levels.iter().enumerate() {
            let input = if i == 0 { &self.offsets } else { &self.levels[i - 1].sums };
            let grid = workgroup_grid_for_threads(level.n, BLOCK_SIZE);
            self.scan_block.record(ctx, encoder, &[input, &level.sums], Some(&level.params), grid);
        }
        // ...then add the scanned block offsets back, deepest level first.
        for i in (0..self.levels.len().saturating_sub(1)).rev() {
            let level = &self.levels[i];
            let input = if i == 0 { &self.offsets } else { &self.levels[i - 1].sums };
            let grid = workgroup_grid_for_threads(level.n, BLOCK_SIZE);
            self.add_block_sums.record(ctx, encoder, &[input, &level.sums], Some(&level.params), grid);
        }

        encoder.copy_buffer_to_buffer(&self.offsets, 0, &self.cursor, 0, cells_bytes);
        self.cell_scatter
            .record(ctx, encoder, &[&self.cell_of, &self.cursor, &self.sorted], Some(&self.scatter_params), grid_n);
    }
}
