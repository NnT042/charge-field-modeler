// DSMC step 2: advance every live photon by c*dt, test the nearest body for
// absorption/lane-crossing, apply the domain boundary, and emit this
// photon's cell index (dead -> bucket `num_cells`) for CellSorter.
// dsmc_common.wgsl is concatenated ahead of this source by engine.rs.

struct GpuPhoton {
    a: vec4<f32>,
    b: vec4<f32>,
    c: vec4<u32>,
    d: vec4<f32>,
};

struct GpuBody {
    center: vec4<f32>,
    axis: vec4<f32>,
    emit: vec4<f32>,
    lane: vec4<f32>,
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
    glancing_bias: u32, _pad12a: u32, _pad12b: u32, _pad12c: u32,
};

const FLAG_ALIVE: u32 = 1u;
const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const SOURCE_SHIFT: u32 = 3u;
const SOURCE_MASK: u32 = 3u << SOURCE_SHIFT;
const CROSSED_SHIFT: u32 = 5u;
const ARRIVED_SHIFT: u32 = 9u;
const FLAG_COLLIDED_SINCE_CROSS: u32 = 1u << 13u;
const HIST_BINS: u32 = 36u;
const HIST_MAX_BODIES: u32 = 4u;
const HIST_ARRIVAL_BLOCK: u32 = HIST_MAX_BODIES * HIST_BINS;
const TALLY_SCALE: f32 = 4096.0;
// B2b item 3 (doc §8) tally offsets -- must stay in sync with gpu_types.rs.
const TALLY_ABSORBED_MOM: u32 = 72u;
const TALLY_ESCAPED_BAND: u32 = 84u;
const TALLY_ESCAPED_BAND_SPIN: u32 = 114u;
const TALLY_ESCAPED_CROSSED_BODY: u32 = 144u;
const TALLY_ESCAPED_CROSSED_BODY_SPIN: u32 = 152u;
const TALLY_SURVIVAL: u32 = 160u;
// E1 (doc §11d) -- must stay in sync with gpu_types.rs.
const TALLY_ESCAPED_EMITTED_K_SUM: u32 = 192u;
const TALLY_ESCAPED_EMITTED_K_COUNT: u32 = 193u;
// Doc §13 -- must stay in sync with gpu_types.rs.
const TALLY_ESC_ALL: u32 = 222u;
const TALLY_ESC_SPEED_SUM: u32 = 223u;
const TALLY_ESC_SPEED_SUMSQ: u32 = 224u;
const TALLY_ESC_EDGE_SUM: u32 = 225u;
const TALLY_ESC_EDGE_SUMSQ: u32 = 226u;
const TALLY_ESC_EDGE_OVER_C: u32 = 227u;

@group(0) @binding(0) var<storage, read_write> photons: array<GpuPhoton>;
@group(0) @binding(1) var<storage, read> bodies: array<GpuBody>;
@group(0) @binding(2) var<storage, read_write> cell_of: array<u32>;
@group(0) @binding(3) var<storage, read_write> tallies: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> hist: array<atomic<u32>>;
@group(0) @binding(5) var<storage, read_write> density: array<atomic<u32>>;
@group(0) @binding(6) var<uniform> params: Params;

fn hist_bin(lat_deg: f32) -> u32 {
    let raw = i32(floor((lat_deg + 90.0) / 5.0));
    return u32(clamp(raw, 0, i32(HIST_BINS) - 1));
}

fn species_idx(flags: u32) -> u32 {
    return select(0u, 1u, (flags & FLAG_SPECIES_PLUS) != 0u);
}

// Doc §8: 0 = polar (>60 deg), 1 = mid (30-60 deg), 2 = equatorial (<30 deg).
fn escape_band(abs_lat_deg: f32) -> u32 {
    if (abs_lat_deg > 60.0) {
        return 0u;
    } else if (abs_lat_deg >= 30.0) {
        return 1u;
    }
    return 2u;
}

@compute @workgroup_size(256)
fn main(@builtin(workgroup_id) wid: vec3<u32>, @builtin(local_invocation_id) lid: vec3<u32>, @builtin(num_workgroups) nwg: vec3<u32>) {
    let linear_wg = wid.x + wid.y * nwg.x;
    let idx = linear_wg * 256u + lid.x;
    if (idx >= params.n_live) {
        return;
    }

    var ph = photons[idx];
    if ((ph.c[0] & FLAG_ALIVE) == 0u) {
        cell_of[idx] = params.num_cells;
        return;
    }

    var pos = ph.a.xyz;
    let dir = ph.b.xyz;
    let id = bitcast<u32>(ph.a.w);
    let kind_emitted = (ph.c[0] & FLAG_KIND_EMITTED) != 0u;
    let spec = species_idx(ph.c[0]);
    let source_body = (ph.c[0] & SOURCE_MASK) >> SOURCE_SHIFT;

    // Path length this step: |v|*dt (|v| = 1 except RoughSphere+Free).
    let speed = ph.d.w;
    var remaining = params.dt * speed;
    var hit_body: i32 = -1;
    var hit_t = remaining;
    for (var i = 0u; i < params.n_bodies; i = i + 1u) {
        let b = bodies[i];
        let t = ray_sphere_hit(pos, dir, b.center.xyz, b.center.w, remaining);
        if (t >= 0.0 && t < hit_t) {
            hit_t = t;
            hit_body = i32(i);
        }
    }

    var dead = false;
    var crossed_any = false;

    if (hit_body >= 0) {
        let i = u32(hit_body);
        let body = bodies[i];
        let hit_pos = pos + dir * hit_t;
        let r_in = body.center.w;

        // Counter-stream survival (doc §8, "measured not integrated"): BEFORE
        // this event changes any flags, snapshot which lanes this photon has
        // crossed and whether it has collided since its last crossing. A
        // photon that crossed body k's lane earlier and now ARRIVES at a
        // DIFFERENT body i (whether it is absorbed here or crosses again) is
        // one sample of "did it survive the gap uninterrupted".
        let crossed_mask_before = (ph.c[0] >> CROSSED_SHIFT) & 0xFu;
        let collided_before = select(0u, 1u, (ph.c[0] & FLAG_COLLIDED_SINCE_CROSS) != 0u);
        for (var k = 0u; k < 4u; k = k + 1u) {
            if (((crossed_mask_before >> k) & 1u) != 0u && k != i) {
                atomicAdd(&tallies[TALLY_SURVIVAL + (k * 4u + i) * 2u + collided_before], 1);
            }
        }

        // Arrival latitude in the hit body's own frame (ambient only, per doc
        // §9 G0). Tallied once per body (first contact only), matching the
        // CPU reference's `probe_arrived` guard in `march_stack`.
        let arrived_bit = 1u << (ARRIVED_SHIFT + i);
        if (!kind_emitted && (ph.c[0] & arrived_bit) == 0u) {
            ph.c[0] = ph.c[0] | arrived_bit;
            let d0 = hit_pos - body.center.xyz;
            let colat0 = acos(clamp(dot(d0, body.axis.xyz) / r_in, -1.0, 1.0));
            let lat0 = 90.0 - degrees(colat0);
            atomicAdd(&hist[HIST_ARRIVAL_BLOCK + i * HIST_BINS + hist_bin(lat0)], 1u);
        }

        let through_lane_rho = body.lane.x;
        let cos_max_angle = body.lane.y;
        var crossed = false;
        // Doc §2/§7: "a photon may cross body i's lane unless it was emitted
        // by body i" -- applies to every photon, not just ambient (an
        // emitted photon from a NEIGHBOUR must be able to cross this body's
        // lane, which is the whole through-exhaust mechanism a p-n-p stack
        // measures). Only self-emission is excluded.
        let self_lane = kind_emitted && source_body == i;
        if (!self_lane && through_lane_rho > 0.0) {
            let d = hit_pos - body.center.xyz;
            let axial = dot(dir, body.axis.xyz);
            let along = dot(d, body.axis.xyz);
            let rho = length(d - body.axis.xyz * along);
            let incoming = dot(dir, d / r_in);
            if (rho < through_lane_rho * r_in && abs(axial) > cos_max_angle && incoming < 0.0) {
                let t2 = -2.0 * dot(d, dir);
                if (t2 > 0.0) {
                    let exitv = d + dir * t2;
                    let exit_len = max(length(exitv), 1e-6);
                    pos = body.center.xyz + exitv * (r_in / exit_len) * 1.001;
                    crossed = true;
                    crossed_any = true;
                    ph.c[0] = ph.c[0] | (1u << (CROSSED_SHIFT + i));
                    // Doc §8: crossing a lane CLEARS "collided since last
                    // crossing" -- the counter-stream survival window resets.
                    ph.c[0] = ph.c[0] & ~FLAG_COLLIDED_SINCE_CROSS;
                    atomicAdd(&tallies[32 + i * 2 + spec], 1); // TALLY_CROSSINGS
                }
            }
        }

        if (!crossed) {
            // Absorbed.
            dead = true;
            let already_crossed = (ph.c[0] >> CROSSED_SHIFT) != 0u;
            let kind_i = select(0u, 1u, kind_emitted);
            atomicAdd(&tallies[i * 4u + kind_i * 2u + spec], 1); // TALLY_ABSORBED
            if (already_crossed) {
                atomicAdd(&tallies[16 + i * 4u + kind_i * 2u + spec], 1); // TALLY_ABSORBED_CROSSED
            }
            // Absorbed momentum (doc §8): Sigma dir.xyz per body -- the
            // open-field absorber force observable.
            atomicAdd(&tallies[TALLY_ABSORBED_MOM + i * 3u + 0u], i32(dir.x * TALLY_SCALE));
            atomicAdd(&tallies[TALLY_ABSORBED_MOM + i * 3u + 1u], i32(dir.y * TALLY_SCALE));
            atomicAdd(&tallies[TALLY_ABSORBED_MOM + i * 3u + 2u], i32(dir.z * TALLY_SCALE));
        } else {
            // Continue the remaining path length from the far side (dt is
            // small relative to body spacing, so a single hit test per step
            // is an adequate approximation -- see engine.rs doc comment).
            pos = pos + dir * (remaining - hit_t);
        }
    } else {
        pos = pos + dir * remaining;
    }

    if (!dead) {
        let center = vec3<f32>(params.domain_center_x, params.domain_center_y, params.domain_center_z);
        if (params.boundary_mode == 1u) {
            // PeriodicBox: wrap.
            let half = params.boundary_extent;
            let full = 2.0 * half;
            var rel = pos - center;
            rel = rel - full * floor((rel + half) / full);
            pos = center + rel;
        } else {
            if (length(pos - center) > params.boundary_extent) {
                dead = true;
                let kind_i = select(0u, 1u, kind_emitted);
                let crossed_i = select(0u, 1u, crossed_any || (ph.c[0] >> CROSSED_SHIFT) != 0u);
                atomicAdd(&tallies[40 + source_body * 4u + crossed_i * 2u + spec], 1); // TALLY_ESCAPED (legacy shape, kept)

                // Doc §8: escape topology by source-body-slot x species x
                // |latitude| band. Emitted -> that body's own frame; ambient
                // (no source body) -> the stack frame (axis Z through the
                // domain centroid, per doc "ambient source = none -> use
                // band in stack frame").
                var band_axis = vec3<f32>(0.0, 0.0, 1.0);
                var band_center = center;
                var slot = 4u; // "none"
                if (kind_emitted) {
                    let sb = bodies[source_body];
                    band_axis = sb.axis.xyz;
                    band_center = sb.center.xyz;
                    slot = source_body;
                }
                let rel_e = pos - band_center;
                let rlen = max(length(rel_e), 1e-6);
                let colat_e = acos(clamp(dot(rel_e, band_axis) / rlen, -1.0, 1.0));
                let lat_e = 90.0 - degrees(colat_e);
                let band = escape_band(abs(lat_e));
                let eb_idx = slot * 2u * 3u + spec * 3u + band;
                atomicAdd(&tallies[TALLY_ESCAPED_BAND + eb_idx], 1);
                atomicAdd(&tallies[TALLY_ESCAPED_BAND_SPIN + eb_idx], i32(ph.b.w * TALLY_SCALE));

                // Of escaped photons, per-crossed-body-i breakdown (doc §8:
                // "keep crossed body i for i=0..3"). NOT mutually exclusive:
                // a photon that crossed more than one lane increments more
                // than one counter.
                let cmask = (ph.c[0] >> CROSSED_SHIFT) & 0xFu;
                for (var k = 0u; k < 4u; k = k + 1u) {
                    if (((cmask >> k) & 1u) != 0u) {
                        atomicAdd(&tallies[TALLY_ESCAPED_CROSSED_BODY + k * 2u + spec], 1);
                        atomicAdd(&tallies[TALLY_ESCAPED_CROSSED_BODY_SPIN + k * 2u + spec], i32(ph.b.w * TALLY_SCALE));
                    }
                }

                // Doc §13 (D7): speed and spin-edge-speed spread of
                // everything that escapes.
                let edge = length(ph.d.xyz) * exp2(f32(ph.c[3]));
                atomicAdd(&tallies[TALLY_ESC_ALL], 1);
                atomicAdd(&tallies[TALLY_ESC_SPEED_SUM], i32(speed * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_ESC_SPEED_SUMSQ], i32(speed * speed * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_ESC_EDGE_SUM], i32(edge * TALLY_SCALE));
                atomicAdd(&tallies[TALLY_ESC_EDGE_SUMSQ], i32(edge * edge * TALLY_SCALE));
                if (edge > 1.0 + 1e-4) {
                    atomicAdd(&tallies[TALLY_ESC_EDGE_OVER_C], 1);
                }

                // E1 (doc §11d): mean stiffness level of escaping EMITTED
                // photons (plain integer sum -- k in 0..7, no scale needed).
                if (kind_emitted) {
                    atomicAdd(&tallies[TALLY_ESCAPED_EMITTED_K_SUM], i32(ph.c[3]));
                    atomicAdd(&tallies[TALLY_ESCAPED_EMITTED_K_COUNT], 1);
                }
            }
        }
    }

    if (dead) {
        ph.c[0] = ph.c[0] & ~FLAG_ALIVE;
        photons[idx] = ph;
        cell_of[idx] = params.num_cells;
        return;
    }

    ph.a = vec4<f32>(pos, ph.a.w);
    photons[idx] = ph;

    // Isotropy tallies (G1): fixed-point sum of dir.xyz and dz^2 over every
    // live photon, sampled once per step at move time (post-collide from the
    // previous step).
    atomicAdd(&tallies[63], i32(dir.x * TALLY_SCALE));
    atomicAdd(&tallies[64], i32(dir.y * TALLY_SCALE));
    atomicAdd(&tallies[65], i32(dir.z * TALLY_SCALE));
    atomicAdd(&tallies[66], i32(dir.z * dir.z * TALLY_SCALE));

    let gd = vec3<i32>(i32(params.grid_dim_x), i32(params.grid_dim_y), i32(params.grid_dim_z));
    let rel = (pos - vec3<f32>(params.grid_min_x, params.grid_min_y, params.grid_min_z)) / params.cell_h;
    let cx = clamp(i32(floor(rel.x)), 0, gd.x - 1);
    let cy = clamp(i32(floor(rel.y)), 0, gd.y - 1);
    let cz = clamp(i32(floor(rel.z)), 0, gd.z - 1);
    let cell_idx = u32(cx) + u32(cy) * params.grid_dim_x + u32(cz) * params.grid_dim_x * params.grid_dim_y;
    cell_of[idx] = cell_idx;

    // Density map (doc §8): time-averaged counts by kind x species, u32
    // per-cell accumulators added each step over the averaging window ONLY
    // (density_accumulate is 0 during warm-up). Never reset mid-window --
    // the host clears `density` once when the window starts and reads it
    // back once at the end.
    if (params.density_accumulate != 0u) {
        let kind_i2 = select(0u, 1u, kind_emitted);
        atomicAdd(&density[cell_idx * 4u + kind_i2 * 2u + spec], 1u);
    }
}
