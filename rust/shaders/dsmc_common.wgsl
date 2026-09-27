// Shared WGSL for the DSMC engine kernels. Concatenated host-side ahead of
// each kernel's own source (naga/WGSL has no #include) -- see
// rust/src/gpu/dsmc/engine.rs `load_kernel`.
//
// RNG: counter-based hash (murmur3-style finalizer), NOT the CPU's
// xorshift64* -- the gates only require statistical agreement (3-sigma/2%),
// not a bit-identical stream, so a GPU-native hash keeps every thread's draw
// independent of scheduling order while staying deterministic in (seed,
// step, id, draw).

fn hash_u32(x: u32) -> u32 {
    var v = x;
    v = v ^ (v >> 16u);
    v = v * 0x7feb352du;
    v = v ^ (v >> 15u);
    v = v * 0x846ca68bu;
    v = v ^ (v >> 16u);
    return v;
}

fn rand01(seed_hi: u32, seed_lo: u32, id: u32, draw: u32) -> f32 {
    var h = hash_u32(seed_hi ^ 0x9e3779b9u);
    h = hash_u32(h ^ seed_lo);
    h = hash_u32(h ^ id);
    h = hash_u32(h ^ draw);
    return f32(h) * (1.0 / 4294967296.0);
}

const PI: f32 = 3.14159265358979;
const TAU: f32 = 6.28318530717959;

// Box-Muller, consumes draws [draw, draw+1].
fn rand_gaussian(seed_hi: u32, seed_lo: u32, id: u32, draw: u32, mean: f32, sigma: f32) -> f32 {
    let u1 = max(rand01(seed_hi, seed_lo, id, draw), 1e-7);
    let u2 = rand01(seed_hi, seed_lo, id, draw + 1u);
    let z0 = sqrt(-2.0 * log(u1)) * cos(TAU * u2);
    return mean + sigma * z0;
}

fn rand_unit_vec(seed_hi: u32, seed_lo: u32, id: u32, draw: u32) -> vec3<f32> {
    let z = rand01(seed_hi, seed_lo, id, draw) * 2.0 - 1.0;
    let phi = rand01(seed_hi, seed_lo, id, draw + 1u) * TAU;
    let r = sqrt(max(1.0 - z * z, 0.0));
    return vec3<f32>(r * cos(phi), r * sin(phi), z);
}

// Orthonormal helper vector perpendicular to `a` (a must be unit).
fn perp_helper(a: vec3<f32>) -> vec3<f32> {
    if (abs(a.x) < 0.9) {
        return vec3<f32>(1.0, 0.0, 0.0);
    }
    return vec3<f32>(0.0, 1.0, 0.0);
}

// Grazing-biased contact normal about rel_hat (recycling.rs sample_contact_normal):
// u skewed toward 1 (u = 1-(1-u0)^3) so hits are almost-tangent (arp.pdf).
fn sample_contact_normal(rel_hat: vec3<f32>, seed_hi: u32, seed_lo: u32, id: u32, draw: u32) -> vec3<f32> {
    let helper = perp_helper(rel_hat);
    let t = normalize(cross(rel_hat, helper));
    let bperp = cross(rel_hat, t);
    let psi = rand01(seed_hi, seed_lo, id, draw) * TAU;
    let b = t * cos(psi) + bperp * sin(psi);
    let u0 = rand01(seed_hi, seed_lo, id, draw + 1u);
    let one_minus_u0 = 1.0 - u0;
    let u = 1.0 - one_minus_u0 * one_minus_u0 * one_minus_u0;
    let nhat = b * sqrt(u) + rel_hat * sqrt(1.0 - u);
    let len2 = dot(nhat, nhat);
    if (len2 > 1e-10) {
        return nhat / sqrt(len2);
    }
    return rel_hat;
}

// Cosine-weighted direction about inward normal `n` (Malley's method), for
// the sphere-boundary ambient injection flux.
fn rand_cosine_hemisphere(n: vec3<f32>, seed_hi: u32, seed_lo: u32, id: u32, draw: u32) -> vec3<f32> {
    let u1 = rand01(seed_hi, seed_lo, id, draw);
    let u2 = rand01(seed_hi, seed_lo, id, draw + 1u);
    let r = sqrt(u1);
    let phi = TAU * u2;
    let lx = r * cos(phi);
    let ly = r * sin(phi);
    let lz = sqrt(max(1.0 - u1, 0.0));
    let helper = perp_helper(n);
    let t = normalize(cross(n, helper));
    let bp = cross(n, t);
    return normalize(t * lx + bp * ly + n * lz);
}

// Ray (p + s*d, s in [0, max_s]) vs unit sphere centred at c radius r_in.
// Returns the smallest positive hit distance <= max_s, or -1.0 if none.
fn ray_sphere_hit(p: vec3<f32>, d: vec3<f32>, c: vec3<f32>, r_in: f32, max_s: f32) -> f32 {
    let oc = p - c;
    let b = 2.0 * dot(oc, d);
    let cc = dot(oc, oc) - r_in * r_in;
    let disc = b * b - 4.0 * cc;
    if (disc < 0.0) {
        return -1.0;
    }
    let sq = sqrt(disc);
    let t1 = (-b - sq) * 0.5;
    let t2 = (-b + sq) * 0.5;
    var t = -1.0;
    if (t1 > 1e-6 && t1 <= max_s) {
        t = t1;
    } else if (t2 > 1e-6 && t2 <= max_s) {
        t = t2;
    }
    return t;
}
