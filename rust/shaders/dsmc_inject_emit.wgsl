// DSMC step 1: inject ambient photons (boundary flux, or full-box fill for
// PeriodicBox) and emit each body's charge, into slots
// [slot_base, slot_base + n_ambient_inject + n_emit_total).
// dsmc_common.wgsl is concatenated ahead of this source by engine.rs.
//
// Tally/histogram layout matches rust/src/gpu/dsmc/gpu_types.rs exactly
// (TALLY_* / HIST_* constants) -- keep the two in sync by hand.

struct GpuPhoton {
    a: vec4<f32>,
    b: vec4<f32>,
    c: vec4<u32>,
    d: vec4<f32>,
    e: vec4<f32>,
    f: vec4<f32>,
    g: vec4<f32>,
};

struct GpuBody {
    center: vec4<f32>,
    axis: vec4<f32>,
    emit: vec4<f32>,
    lane: vec4<f32>,
};

struct InjectParams {
    n_ambient_inject: u32,
    n_emit_total: u32,
    id_base: u32,
    slot_base: u32,
};

struct Params {
    n_bodies: u32, n_live: u32, capacity: u32, step: u32,
    dt: f32, sigma: f32, n_ambient: f32, photon_fraction: f32,
    pole_intake_bias: f32, cancel_redirect: f32, spin_transfer: f32, b_max: f32,
    boundary_mode: u32, boundary_extent: f32, rate_kernel: u32, mesh_mode: u32,
    bounce_cap: u32, d5_far_zone_gate: u32, grid_dim_x: u32, grid_dim_y: u32,
    grid_dim_z: u32, cell_h: f32, grid_min_x: f32, grid_min_y: f32,
    grid_min_z: f32, seed_hi: u32, seed_lo: u32, num_cells: u32,
    domain_center_x: f32, domain_center_y: f32, domain_center_z: f32, exchange_rule: u32,
    density_accumulate: u32, emit_stiffness_level: u32, stiffness_shed: u32, sigma_size_scaling: u32,
    debug_force_equal_mass_exchange: u32, g_max_catch: f32, rough_beta: f32, rough_beta_twist: f32,
    rough_kappa: f32, twist_loss_to_normal: u32, speed_mode: u32, spin_edge0: f32,
    glancing_bias: u32, stacking: u32, stack_max_level: u32, unstack_edge: f32,
    body_model: u32, tau_body: f32, grain_edge: f32, body_twist_loss: u32,
    body_beta: f32, body_beta_twist: f32, _pad14a: u32, _pad14b: u32,
};

const FLAG_ALIVE: u32 = 1u;
const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const SOURCE_SHIFT: u32 = 3u;

const HIST_BINS: u32 = 36u;
const HIST_MAX_BODIES: u32 = 4u;
const HIST_EMISSION_BLOCK: u32 = 0u;

@group(0) @binding(0) var<storage, read_write> photons: array<GpuPhoton>;
@group(0) @binding(1) var<storage, read> bodies: array<GpuBody>;
@group(0) @binding(2) var<storage, read> body_emit_counts: array<u32>;
@group(0) @binding(3) var<storage, read_write> tallies: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> hist: array<atomic<u32>>;
@group(0) @binding(5) var<storage, read> inject_params: InjectParams;
@group(0) @binding(6) var<uniform> params: Params;

fn hist_bin(lat_deg: f32) -> u32 {
    let raw = i32(floor((lat_deg + 90.0) / 5.0));
    return u32(clamp(raw, 0, i32(HIST_BINS) - 1));
}

fn latlon_to_unit(t: vec3<f32>, bb: vec3<f32>, axis: vec3<f32>, lat_deg: f32, az: f32) -> vec3<f32> {
    let lat = radians(lat_deg);
    let colat = 1.5707963 - lat;
    let s = sin(colat);
    let cz = cos(colat);
    return t * (s * cos(az)) + bb * (s * sin(az)) + axis * cz;
}

fn sample_exit_deg(body: GpuBody, id: u32, base_draw: u32, peak_select_draw: u32) -> f32 {
    let peak0 = body.emit.x;
    let peak1 = body.emit.y;
    var peak = peak0;
    if (peak1 >= 0.0) {
        let u = rand01(params.seed_hi, params.seed_lo, id, peak_select_draw);
        if (u >= 0.5) {
            peak = peak1;
        }
    }
    let mag = rand_gaussian(params.seed_hi, params.seed_lo, id, base_draw, peak, body.emit.z);
    return clamp(mag, 0.0, 89.0);
}

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let idx = linear_wg * 256u + lid.x;
    let total = inject_params.n_ambient_inject + inject_params.n_emit_total;
    if (idx >= total) {
        return;
    }

    let slot = inject_params.slot_base + idx;
    if (slot >= params.capacity) {
        return;
    }
    let id = inject_params.id_base + idx;

    var out: GpuPhoton;

    if (idx < inject_params.n_ambient_inject) {
        // --- ambient ---
        let chi_u = rand01(params.seed_hi, params.seed_lo, id, 0u);
        let chi: f32 = select(-1.0, 1.0, chi_u < params.photon_fraction);
        let center = vec3<f32>(params.domain_center_x, params.domain_center_y, params.domain_center_z);

        var pos: vec3<f32>;
        var dir: vec3<f32>;
        if (params.boundary_mode == 1u) {
            // PeriodicBox: uniform fill of the whole box, isotropic direction
            // (used once at t=0 -- there is no ongoing boundary flux to inject
            // in a periodic domain, the population is exactly conserved).
            let ux = rand01(params.seed_hi, params.seed_lo, id, 1u) * 2.0 - 1.0;
            let uy = rand01(params.seed_hi, params.seed_lo, id, 2u) * 2.0 - 1.0;
            let uz = rand01(params.seed_hi, params.seed_lo, id, 3u) * 2.0 - 1.0;
            pos = center + vec3<f32>(ux, uy, uz) * params.boundary_extent;
            dir = rand_unit_vec(params.seed_hi, params.seed_lo, id, 4u);
        } else {
            // Boundary inflow of a HOMOGENEOUS field whose direction
            // distribution may be tilted (species intake bias). Direction is
            // drawn FIRST from the field's own distribution (spawn_ambient's
            // rule), then the line is placed uniformly over the disk of radius
            // R_D perpendicular to it and started at its sphere entry -- the
            // Phase-1 construction with B_MAX -> R_D. Drawing a cosine-weighted
            // inward direction per boundary point and tilting it afterwards is
            // only equivalent for an untilted field (it skewed the axial lane
            // flux in the tilted p-n-p G0 by ~40 sigma). Rate stays
            // n_a c pi R_D^2 for any direction distribution.
            var d = rand_unit_vec(params.seed_hi, params.seed_lo, id, 1u);
            if (params.pole_intake_bias > 0.0) {
                d = normalize(d + vec3<f32>(0.0, 0.0, params.pole_intake_bias * chi));
            }
            let t_perp = normalize(cross(d, perp_helper(d)));
            let b_perp = cross(d, t_perp);
            let psi = rand01(params.seed_hi, params.seed_lo, id, 3u) * TAU;
            let e = t_perp * cos(psi) + b_perp * sin(psi);
            let big_r = params.boundary_extent;
            let b = big_r * sqrt(rand01(params.seed_hi, params.seed_lo, id, 4u));
            let t_entry = big_r - sqrt(max(big_r * big_r - b * b, 0.0));
            // Nudge inside so the move kernel's escape test doesn't fire at once.
            pos = center - d * big_r + e * b + d * (t_entry + 1e-4);
            dir = d;
        }

        out.a = vec4<f32>(pos, bitcast<f32>(id));
        out.b = vec4<f32>(dir, chi);
        var flags = FLAG_ALIVE;
        if (chi > 0.0) {
            flags = flags | FLAG_SPECIES_PLUS;
        }
        // E1 (doc §11d): ambient photons are always born at stiffness k=0.
        out.c = vec4<u32>(flags, 0u, 0u, 0u);
        // Doc §13: helical spin along travel at edge speed spin_edge0
        // (photon2.html), sense = chirality (clockwise from behind = +dir,
        // quark.html right-hand rule). R = 1 at k = 0. Speed c.
        out.d = vec4<f32>(dir * (chi * params.spin_edge0), 1.0);
        out.e = vec4<f32>(0.0);
        out.f = vec4<f32>(0.0);
        out.g = vec4<f32>(0.0);
        atomicAdd(&tallies[60], 1); // TALLY_INJECTED
    } else {
        // --- emitted: find owning body by walking per-body counts ---
        var remaining = idx - inject_params.n_ambient_inject;
        var body_idx = 0u;
        for (var i = 0u; i < params.n_bodies; i = i + 1u) {
            let c = body_emit_counts[i];
            if (remaining < c) {
                body_idx = i;
                break;
            }
            remaining = remaining - c;
        }
        let body = bodies[body_idx];
        let axis = body.axis.xyz;
        let chirality = body.axis.w;
        let helper = perp_helper(axis);
        let t = normalize(cross(axis, helper));
        let bb = cross(axis, t);

        let sign_u = rand01(params.seed_hi, params.seed_lo, id, 0u);
        let sgn: f32 = select(-1.0, 1.0, sign_u < 0.5);

        let mag0 = sample_exit_deg(body, id, 1u, 2u);
        let lat0 = sgn * mag0;
        let az0 = rand01(params.seed_hi, params.seed_lo, id, 3u) * TAU;
        let r_hat0 = latlon_to_unit(t, bb, axis, lat0, az0);

        let mag1 = sample_exit_deg(body, id, 4u, 5u);
        let lat1 = sgn * mag1;
        let dir_radial = latlon_to_unit(t, bb, axis, lat1, az0);
        let cax = cross(axis, dir_radial);
        let phi_hat_len2 = dot(cax, cax);
        var phi_hat = vec3<f32>(1.0, 0.0, 0.0);
        if (phi_hat_len2 > 1e-10) {
            phi_hat = cax / sqrt(phi_hat_len2);
        }
        let swing = body.emit.w;
        let dir = normalize(dir_radial + phi_hat * swing);
        let pos = body.center.xyz + r_hat0 * 1.001;

        out.a = vec4<f32>(pos, bitcast<f32>(id));
        out.b = vec4<f32>(dir, chirality);
        var flags = FLAG_ALIVE | FLAG_KIND_EMITTED | (body_idx << SOURCE_SHIFT);
        if (chirality > 0.0) {
            flags = flags | FLAG_SPECIES_PLUS;
        }
        // E1 (doc §11d): emitted photons are born at stiffness level k0
        // (params.emit_stiffness_level, default 0 -- HEAD behaviour).
        out.c = vec4<u32>(flags, 0u, 0u, params.emit_stiffness_level);
        // Doc §13: same helical birth spin, radius R = 2^k0.
        let r_emit = exp2(f32(params.emit_stiffness_level));
        out.d = vec4<f32>(dir * (chirality * params.spin_edge0 / r_emit), 1.0);
        out.e = vec4<f32>(0.0);
        out.f = vec4<f32>(0.0);
        out.g = vec4<f32>(0.0);

        atomicAdd(&tallies[56u + body_idx], 1); // TALLY_EMITTED
        let bin = hist_bin(lat0);
        atomicAdd(&hist[HIST_EMISSION_BLOCK + body_idx * HIST_BINS + bin], 1u);
    }

    photons[slot] = out;
}
