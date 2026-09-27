use std::borrow::Cow;

use bytemuck::Pod;
use wgpu::util::DeviceExt;

/// A headless wgpu device/queue pair, usable from `cargo test` without Godot's
/// `RenderingDevice` (which does not exist in headless mode).
pub struct GpuContext {
    pub instance: wgpu::Instance,
    pub adapter: wgpu::Adapter,
    pub device: wgpu::Device,
    pub queue: wgpu::Queue,
    pub adapter_name: String,
    pub backend: String,
    pub limits: wgpu::Limits,
    pub supports_timestamps: bool,
}

impl GpuContext {
    /// Create a headless GPU context. Prefers Vulkan, falls back to any other
    /// backend wgpu can find. Returns `None` (never panics) if no adapter is
    /// available at all, so tests/benches can skip gracefully on CI or
    /// GPU-less machines.
    pub fn new() -> Option<Self> {
        // Try Vulkan first (the target backend for this project), then fall
        // back to whatever else wgpu can find headlessly.
        Self::new_with_backends(wgpu::Backends::VULKAN)
            .or_else(|| Self::new_with_backends(wgpu::Backends::PRIMARY))
    }

    fn new_with_backends(backends: wgpu::Backends) -> Option<Self> {
        let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
            backends,
            flags: wgpu::InstanceFlags::default(),
            memory_budget_thresholds: wgpu::MemoryBudgetThresholds::default(),
            backend_options: wgpu::BackendOptions::default(),
            display: None,
        });

        let adapter = pollster::block_on(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::HighPerformance,
            compatible_surface: None,
            force_fallback_adapter: false,
            apply_limit_buckets: false,
        }))
        .ok()?;

        let adapter_info = adapter.get_info();
        let adapter_features = adapter.features();
        let adapter_limits = adapter.limits();

        let mut required_features = wgpu::Features::empty();
        let supports_timestamps = adapter_features.contains(wgpu::Features::TIMESTAMP_QUERY);
        if supports_timestamps {
            required_features |= wgpu::Features::TIMESTAMP_QUERY;
        }
        // naga's GLSL frontend lowers the plain `barrier()` builtin to a
        // Barrier flag set that includes SUB_GROUP, so any GLSL compute
        // shader using barrier() + shared memory (workgroup reductions,
        // scans) needs this feature enabled even though it never uses actual
        // subgroup intrinsics. Without it, shader module creation fails
        // validation ("missing capability ... SUBGROUP | SUBGROUP_BARRIER").
        if adapter_features.contains(wgpu::Features::SUBGROUP) {
            required_features |= wgpu::Features::SUBGROUP;
        }
        if adapter_features.contains(wgpu::Features::SUBGROUP_BARRIER) {
            required_features |= wgpu::Features::SUBGROUP_BARRIER;
        }

        // Request the adapter's own limits so large storage buffers (DSMC
        // particle arrays) are not clamped to wgpu's conservative defaults.
        let required_limits = adapter_limits.clone();

        let (device, queue) = pollster::block_on(adapter.request_device(&wgpu::DeviceDescriptor {
            label: Some("cfm-headless-gpu"),
            required_features,
            required_limits: required_limits.clone(),
            experimental_features: wgpu::ExperimentalFeatures::default(),
            memory_hints: wgpu::MemoryHints::Performance,
            trace: wgpu::Trace::Off,
        }))
        .ok()?;

        Some(Self {
            instance,
            adapter,
            device,
            queue,
            adapter_name: adapter_info.name.clone(),
            backend: format!("{:?}", adapter_info.backend),
            limits: required_limits,
            supports_timestamps,
        })
    }

    /// Create a storage buffer initialised from `data`, usable as both a
    /// shader storage binding and a copy source/destination.
    pub fn storage_buffer_from_slice<T: Pod>(&self, data: &[T], extra_usage: wgpu::BufferUsages) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::cast_slice(data),
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST
                | extra_usage,
        })
    }

    /// Create a zero-initialised storage buffer sized for `count` elements of `T`.
    pub fn zeroed_storage_buffer<T: Pod>(&self, count: usize, extra_usage: wgpu::BufferUsages) -> wgpu::Buffer {
        let size = (count * std::mem::size_of::<T>()) as u64;
        self.device.create_buffer(&wgpu::BufferDescriptor {
            label: None,
            size,
            usage: wgpu::BufferUsages::STORAGE
                | wgpu::BufferUsages::COPY_SRC
                | wgpu::BufferUsages::COPY_DST
                | extra_usage,
            mapped_at_creation: false,
        })
    }

    /// Create a uniform buffer from a `Pod` params struct.
    pub fn uniform_buffer<T: Pod>(&self, data: &T) -> wgpu::Buffer {
        self.device.create_buffer_init(&wgpu::util::BufferInitDescriptor {
            label: None,
            contents: bytemuck::bytes_of(data),
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
        })
    }

    /// Overwrite a buffer's contents (must have COPY_DST usage).
    pub fn write_buffer<T: Pod>(&self, buffer: &wgpu::Buffer, data: &[T]) {
        self.queue.write_buffer(buffer, 0, bytemuck::cast_slice(data));
    }

    /// Read back `count` elements of `T` from a buffer (must have COPY_SRC usage).
    pub fn read_buffer<T: Pod + Default + Clone>(&self, buffer: &wgpu::Buffer, count: usize) -> Vec<T> {
        let size = (count * std::mem::size_of::<T>()) as u64;
        let staging = self.device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("readback-staging"),
            size,
            usage: wgpu::BufferUsages::MAP_READ | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });

        let mut encoder = self.device.create_command_encoder(&wgpu::CommandEncoderDescriptor { label: None });
        encoder.copy_buffer_to_buffer(buffer, 0, &staging, 0, size);
        self.queue.submit(Some(encoder.finish()));

        let slice = staging.slice(..);
        let (tx, rx) = std::sync::mpsc::channel();
        slice.map_async(wgpu::MapMode::Read, move |res| {
            let _ = tx.send(res);
        });
        self.device.poll(wgpu::PollType::wait_indefinitely()).expect("device poll failed");
        rx.recv().expect("map_async channel closed").expect("buffer map failed");

        let data = slice.get_mapped_range().expect("get_mapped_range failed");
        let result = bytemuck::cast_slice::<u8, T>(&data).to_vec();
        drop(data);
        staging.unmap();
        result
    }

    /// Compile a GLSL compute shader. Strips a leading Godot `#[compute]`
    /// pragma line if present (godot/shaders/compute/*.glsl all start with
    /// one; naga's GLSL frontend does not understand it).
    pub fn load_glsl_compute(&self, src: &str) -> wgpu::ShaderModule {
        let cleaned = strip_compute_pragma(src);
        self.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: None,
            source: wgpu::ShaderSource::Glsl {
                shader: Cow::Owned(cleaned),
                stage: wgpu::naga::ShaderStage::Compute,
                defines: &[],
            },
        })
    }

    /// Compile a WGSL compute shader.
    pub fn load_wgsl_compute(&self, src: &str) -> wgpu::ShaderModule {
        self.device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: None,
            source: wgpu::ShaderSource::Wgsl(Cow::Borrowed(src)),
        })
    }
}

/// Strip a leading `#[compute]` line (Godot's RenderingDevice shader pragma)
/// so the source is plain `#version 450` GLSL that naga's frontend accepts.
pub fn strip_compute_pragma(src: &str) -> String {
    let trimmed = src.trim_start();
    if let Some(rest) = trimmed.strip_prefix("#[compute]") {
        rest.trim_start_matches(['\r', '\n']).to_string()
    } else {
        src.to_string()
    }
}
