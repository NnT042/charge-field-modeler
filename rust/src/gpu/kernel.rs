use std::time::{Duration, Instant};

use crate::gpu::context::GpuContext;

/// wgpu (and the underlying APIs) cap each dispatch dimension at 65535
/// workgroups. Kernels here that may need more than that for large N read
/// back the actual dispatch shape via the `num_workgroups` builtin
/// (`gl_NumWorkGroups` in GLSL) rather than assuming a 1-D grid, so pass the
/// result of this helper straight to `dispatch`/`dispatch_timed`.
pub fn workgroup_grid_for_threads(total_threads: u32, local_size: u32) -> (u32, u32, u32) {
    const MAX_DIM: u32 = 65535;
    let total_workgroups = total_threads.div_ceil(local_size).max(1);
    if total_workgroups <= MAX_DIM {
        (total_workgroups, 1, 1)
    } else {
        let x = MAX_DIM;
        let y = total_workgroups.div_ceil(x);
        (x, y, 1)
    }
}

/// A compute pipeline plus a timestamp query pair for GPU-side timing.
///
/// Binding convention: `dispatch`/`dispatch_timed` bind `buffers[0..]` to
/// bindings `0..buffers.len()` and, if `uniform` is given, it goes at binding
/// `buffers.len()`. All bindings live in bind group 0 (matching the
/// `set = 0` used throughout `rust/shaders/*.glsl` and
/// `godot/shaders/compute/*.glsl`). This means every kernel's storage buffers
/// must be declared in that same order, with the uniform (if any) declared
/// last.
pub struct Kernel {
    pipeline: wgpu::ComputePipeline,
    query_set: Option<wgpu::QuerySet>,
}

impl Kernel {
    pub fn new(ctx: &GpuContext, shader: &wgpu::ShaderModule, entry_point: &str) -> Self {
        let pipeline = ctx.device.create_compute_pipeline(&wgpu::ComputePipelineDescriptor {
            label: None,
            layout: None,
            module: shader,
            entry_point: Some(entry_point),
            compilation_options: Default::default(),
            cache: None,
        });

        let query_set = if ctx.supports_timestamps {
            Some(ctx.device.create_query_set(&wgpu::QuerySetDescriptor {
                label: Some("kernel-timestamps"),
                ty: wgpu::QueryType::Timestamp,
                count: 2,
            }))
        } else {
            None
        };

        Self { pipeline, query_set }
    }

    fn bind_group(&self, ctx: &GpuContext, buffers: &[&wgpu::Buffer], uniform: Option<&wgpu::Buffer>) -> wgpu::BindGroup {
        let layout = self.pipeline.get_bind_group_layout(0);
        let mut entries: Vec<wgpu::BindGroupEntry> = buffers
            .iter()
            .enumerate()
            .map(|(i, b)| wgpu::BindGroupEntry {
                binding: i as u32,
                resource: b.as_entire_binding(),
            })
            .collect();
        if let Some(u) = uniform {
            entries.push(wgpu::BindGroupEntry {
                binding: buffers.len() as u32,
                resource: u.as_entire_binding(),
            });
        }
        ctx.device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: None,
            layout: &layout,
            entries: &entries,
        })
    }

    /// Record a dispatch into an existing encoder without submitting, so a
    /// multi-pass pipeline (e.g. `CellSorter`) goes to the GPU as ONE submit
    /// with no host round trips. Same binding convention as `dispatch`.
    pub fn record(
        &self,
        ctx: &GpuContext,
        encoder: &mut wgpu::CommandEncoder,
        buffers: &[&wgpu::Buffer],
        uniform: Option<&wgpu::Buffer>,
        workgroups: (u32, u32, u32),
    ) {
        let bind_group = self.bind_group(ctx, buffers, uniform);
        let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
            label: None,
            timestamp_writes: None,
        });
        pass.set_pipeline(&self.pipeline);
        pass.set_bind_group(0, &bind_group, &[]);
        pass.dispatch_workgroups(workgroups.0, workgroups.1, workgroups.2);
    }

    /// Dispatch and wait for completion. Returns wall-clock time around
    /// submit+poll (not pure GPU time; see `dispatch_timed`).
    pub fn dispatch(
        &self,
        ctx: &GpuContext,
        buffers: &[&wgpu::Buffer],
        uniform: Option<&wgpu::Buffer>,
        workgroups: (u32, u32, u32),
    ) -> Duration {
        let bind_group = self.bind_group(ctx, buffers, uniform);
        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: None,
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(workgroups.0, workgroups.1, workgroups.2);
        }

        let start = Instant::now();
        ctx.queue.submit(Some(encoder.finish()));
        ctx.device.poll(wgpu::PollType::wait_indefinitely()).expect("device poll failed");
        start.elapsed()
    }

    /// Dispatch and measure GPU execution time via timestamp queries when the
    /// adapter supports them; falls back to wall-clock around submit+poll
    /// otherwise (still returned as `Duration`, just less precise).
    pub fn dispatch_timed(
        &self,
        ctx: &GpuContext,
        buffers: &[&wgpu::Buffer],
        uniform: Option<&wgpu::Buffer>,
        workgroups: (u32, u32, u32),
    ) -> Duration {
        let Some(query_set) = &self.query_set else {
            return self.dispatch(ctx, buffers, uniform, workgroups);
        };

        let bind_group = self.bind_group(ctx, buffers, uniform);
        let resolve_buffer = ctx.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("timestamp-resolve"),
            size: 2 * 8,
            usage: wgpu::BufferUsages::QUERY_RESOLVE | wgpu::BufferUsages::COPY_SRC,
            mapped_at_creation: false,
        });

        let mut encoder = ctx.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        {
            let mut pass = encoder.begin_compute_pass(&wgpu::ComputePassDescriptor {
                label: None,
                timestamp_writes: Some(wgpu::ComputePassTimestampWrites {
                    query_set,
                    beginning_of_pass_write_index: Some(0),
                    end_of_pass_write_index: Some(1),
                }),
            });
            pass.set_pipeline(&self.pipeline);
            pass.set_bind_group(0, &bind_group, &[]);
            pass.dispatch_workgroups(workgroups.0, workgroups.1, workgroups.2);
        }
        encoder.resolve_query_set(query_set, 0..2, &resolve_buffer, 0);
        ctx.queue.submit(Some(encoder.finish()));

        let timestamps: Vec<u64> = ctx.read_buffer(&resolve_buffer, 2);

        let period_ns = ctx.queue.get_timestamp_period() as f64;
        let delta_ticks = timestamps[1].saturating_sub(timestamps[0]);
        Duration::from_nanos((delta_ticks as f64 * period_ns) as u64)
    }
}
