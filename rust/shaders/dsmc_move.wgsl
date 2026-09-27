// DSMC step 2: advance every live photon by c*dt, test the nearest body for
// absorption/lane-crossing, apply the domain boundary, and emit this
// photon's cell index (dead -> bucket `num_cells`) for CellSorter.
// dsmc_common.wgsl is concatenated ahead of this source by engine.rs.

struct GpuPhoton {
    a: vec4<f32>,
    b: vec4<f32>,
    c: vec4<u32>,
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
    domain_center_x: f32, domain_center_y: f32, domain_center_z: f32, _pad2: f32,
};

const FLAG_ALIVE: u32 = 1u;
const FLAG_KIND_EMITTED: u32 = 2u;
const FLAG_SPECIES_PLUS: u32 = 4u;
const SOURCE_SHIFT: u32 = 3u;
const SOURCE_MASK: u32 = 3u << SOURCE_SHIFT;
const CROSSED_SHIFT: u32 = 5u;
const ARRIVED_SHIFT: u32 = 9u;
const HIST_BINS: u32 = 36u;
const HIST_MAX_BODIES: u32 = 4u;
const HIST_ARRIVAL_BLOCK: u32 = HIST_MAX_BODIES * HIST_BINS;
const TALLY_SCALE: f32 = 4096.0;

@group(0) @binding(0) var<storage, read_write> photons: array<GpuPhoton>;
@group(0) @binding(1) var<storage, read> bodies: array<GpuBody>;
@group(0) @binding(2) var<storage, read_write> cell_of: array<u32>;
@group(0) @binding(3) var<storage, read_write> tallies: array<atomic<i32>>;
@group(0) @binding(4) var<storage, read_write> hist: array<atomic<u32>>;
@group(0) @binding(5) var<uniform> params: Params;

fn hist_bin(lat_deg: f32) -> u32 {
    let raw = i32(floor((lat_deg + 90.0) / 5.0));
    return u32(clamp(raw, 0, i32(HIST_BINS) - 1));
}

fn species_idx(flags: u32) -> u32 {
    return select(0u, 1u, (flags & FLAG_SPECIES_PLUS) != 0u);
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

    var remaining = params.dt;
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
                atomicAdd(&tallies[40 + source_body * 4u + crossed_i * 2u + spec], 1); // TALLY_ESCAPED
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
    cell_of[idx] = u32(cx) + u32(cy) * params.grid_dim_x + u32(cz) * params.grid_dim_x * params.grid_dim_y;
}
