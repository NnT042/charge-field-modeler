//! Pod/Zeroable GPU-layout structs shared by host and WGSL (doc §6).

use bytemuck::{Pod, Zeroable};

/// 48-byte photon record, 3 x vec4 (doc §6).
///
/// `a` = pos.xyz, bitcast(id) as f32 in `.w`.
/// `b` = dir.xyz, spin in `.w`.
/// `c` = [flags, bounces, birth_step, spare] (u32 x4).
///
/// flags bits: 0 alive, 1 kind (0 ambient/1 emitted), 2 species (0 -> -1, 1 -> +1),
/// 3-4 source body index, 5-8 crossed-lane mask, 9-12 arrived mask,
/// 13 collided-since-last-crossing, 14-31 spare.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug, Default)]
pub struct GpuPhoton {
    pub a: [f32; 4],
    pub b: [f32; 4],
    pub c: [u32; 4],
}

pub const FLAG_ALIVE: u32 = 1 << 0;
pub const FLAG_KIND_EMITTED: u32 = 1 << 1;
pub const FLAG_SPECIES_PLUS: u32 = 1 << 2;
pub const SOURCE_SHIFT: u32 = 3;
pub const SOURCE_MASK: u32 = 0b11 << SOURCE_SHIFT; // 2 bits -> up to 4 bodies
pub const CROSSED_SHIFT: u32 = 5;
pub const ARRIVED_SHIFT: u32 = 9;
pub const FLAG_COLLIDED_SINCE_CROSS: u32 = 1 << 13;

impl GpuPhoton {
    pub fn dead() -> Self {
        GpuPhoton {
            a: [0.0; 4],
            b: [0.0; 4],
            c: [0; 4],
        }
    }

    pub fn is_alive(&self) -> bool {
        self.c[0] & FLAG_ALIVE != 0
    }

    pub fn id(&self) -> u32 {
        self.a[3].to_bits()
    }

    pub fn pos(&self) -> [f32; 3] {
        [self.a[0], self.a[1], self.a[2]]
    }

    pub fn dir(&self) -> [f32; 3] {
        [self.b[0], self.b[1], self.b[2]]
    }

    pub fn spin(&self) -> f32 {
        self.b[3]
    }

    pub fn species(&self) -> f32 {
        if self.c[0] & FLAG_SPECIES_PLUS != 0 {
            1.0
        } else {
            -1.0
        }
    }

    pub fn is_emitted(&self) -> bool {
        self.c[0] & FLAG_KIND_EMITTED != 0
    }

    pub fn source_body(&self) -> u32 {
        (self.c[0] & SOURCE_MASK) >> SOURCE_SHIFT
    }

    pub fn crossed_mask(&self) -> u32 {
        (self.c[0] >> CROSSED_SHIFT) & 0b1111
    }

    pub fn arrived_mask(&self) -> u32 {
        (self.c[0] >> ARRIVED_SHIFT) & 0b1111
    }
}

/// One static body, GPU-resident (4 x vec4 = 64 bytes).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
pub struct GpuBody {
    /// xyz = centre, w = R_IN (always 1.0 in these units).
    pub center: [f32; 4],
    /// xyz = spin axis (unit), w = chirality sign (+-1).
    pub axis: [f32; 4],
    /// x = peak0_deg, y = peak1_deg (-1 = none), z = sigma_deg, w = swing (signed).
    pub emit: [f32; 4],
    /// x = through_lane_rho, y = cos_max_angle, z = class_id, w = spare.
    pub lane: [f32; 4],
}

/// Global per-step uniform parameters.
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
pub struct GpuParams {
    pub n_bodies: u32,
    pub n_live: u32,
    pub capacity: u32,
    pub step: u32,

    pub dt: f32,
    pub sigma: f32,
    pub n_ambient: f32,
    pub photon_fraction: f32,

    pub pole_intake_bias: f32,
    pub cancel_redirect: f32,
    pub spin_transfer: f32,
    pub b_max: f32,

    /// 0 = SphereEscape, 1 = PeriodicBox.
    pub boundary_mode: u32,
    pub boundary_extent: f32,
    /// 0 = Constant, 1 = Catch.
    pub rate_kernel: u32,
    /// 0 = ChiralityOnly, 1 = DirRel.
    pub mesh_mode: u32,

    pub bounce_cap: u32,
    pub d5_far_zone_gate: u32,
    pub grid_dim_x: u32,
    pub grid_dim_y: u32,

    pub grid_dim_z: u32,
    pub cell_h: f32,
    pub grid_min_x: f32,
    pub grid_min_y: f32,

    pub grid_min_z: f32,
    pub seed_hi: u32,
    pub seed_lo: u32,
    pub num_cells: u32,

    pub domain_center_x: f32,
    pub domain_center_y: f32,
    pub domain_center_z: f32,
    /// D6: 0 = SwapRenorm, 1 = RotateRelative.
    pub exchange_rule: u32,

    /// Density map (doc §8): 1 while the caller wants `dsmc_move.wgsl` to
    /// accumulate this step's live population into `density_buf` (the
    /// steady-state averaging window only -- never during warm-up).
    pub density_accumulate: u32,
    pub _pad9a: u32,
    pub _pad9b: u32,
    pub _pad9c: u32,
}

/// Injection/emission counts + RNG draw base for one step, computed host-side
/// (fractional remainders carried on the host so counts are exact integers).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
pub struct GpuInjectParams {
    pub n_ambient_inject: u32,
    pub n_emit_total: u32,
    pub id_base: u32,
    pub slot_base: u32,
}

/// Per-body emission counts for this step (parallel array to `GpuBody`).
#[repr(C)]
#[derive(Clone, Copy, Pod, Zeroable, Debug)]
pub struct GpuBodyEmit {
    pub count: u32,
    pub slot_base: u32,
    pub _pad0: u32,
    pub _pad1: u32,
}

/// Fixed-point tally accumulator, scale 2^12 (doc §8), reset every step.
pub const TALLY_SCALE: f64 = 4096.0;

/// Per-step atomic tallies: one flat `array<atomic<i32>>` in WGSL (u32 counts
/// fit comfortably in i32 at gate scales; fixed-point sums use `TALLY_SCALE`).
/// Offsets are shared by hand between here and `dsmc_move.wgsl` /
/// `dsmc_collide.wgsl` -- keep the three in sync.
pub const TALLY_ABSORBED: usize = 0; // [body*4 + kind*2 + species], len 16
pub const TALLY_ABSORBED_CROSSED: usize = 16; // same indexing, len 16
pub const TALLY_CROSSINGS: usize = 32; // [body*2 + species], len 8
pub const TALLY_ESCAPED: usize = 40; // [body*4 + crossed_any*2 + species], len 16
pub const TALLY_EMITTED: usize = 56; // [body], len 4
pub const TALLY_INJECTED: usize = 60; // len 1
pub const TALLY_COLLISIONS: usize = 61; // len 1
pub const TALLY_MOM_DRIFT_COUNT: usize = 62; // len 1
pub const TALLY_DIR_SUM: usize = 63; // [x,y,z] fixed-point, len 3
pub const TALLY_DZ2_SUM: usize = 66; // fixed-point, len 1
pub const TALLY_MOM_DRIFT_SUM: usize = 67; // fixed-point, len 1 (SwapRenorm-only renorm-correction magnitude, legacy)
pub const TALLY_CELL_OVERFLOW: usize = 68; // count of cells with live count > dsmc_collide::MAX_CELL this dispatch
pub const TALLY_EXCHANGE_DRIFT_SUM: usize = 69; // fixed-point: |sum((d1'+d2')-(d1+d2))| from the exchange step alone (D6)
pub const TALLY_REDIRECT_DRIFT_SUM: usize = 70; // fixed-point: same, from the gear-redirect step alone
/// Count of accepted collisions where BOTH partners' id-rank within their
/// (sorted) cell are < 64 -- the old kernel's dead-fixed MAX_CELL. In a fair
/// draw over a cell of size n>64 this should be ~(64/n)^2 (a few percent at
/// n~300); a kernel wrongly capped at 64 candidates gives ~1.0 regardless of
/// n, by construction -- a single-step-sufficient discriminator (unlike a
/// multi-step "fraction of photons ever touched", which a moving population
/// washes out within ~20 steps as cell membership reshuffles).
pub const TALLY_LOW_RANK64_COLLISIONS: usize = 71;

// --- B2b work item 3 (doc §8 observables) ---

/// Absorbed momentum: Sigma dir.xyz per body, fixed-point (open-field
/// absorber force). [body*3 + axis], len 12.
pub const TALLY_ABSORBED_MOM: usize = 72;
/// Number of source-body "slots" for the escape tallies below: bodies 0..3
/// plus slot 4 = "none" (ambient has no source body -- doc §8: "ambient
/// source = none -> use band in stack frame").
pub const ESCAPE_SRC_SLOTS: usize = 5;
/// Latitude bands for the exit/escape topology (doc §8): 0 = polar (>60 deg),
/// 1 = mid (30-60 deg), 2 = equatorial (<30 deg), of |latitude| in the
/// relevant frame (source body's own frame if emitted, stack frame -- axis Z
/// through the domain centroid -- if ambient).
pub const ESCAPE_BANDS: usize = 3;
/// Escapes by source-body-slot x species x exit-latitude band.
/// [slot*2*3 + species*3 + band], len 5*2*3 = 30.
pub const TALLY_ESCAPED_BAND: usize = 84;
/// Spin sum (fixed-point) over the same index space, len 30.
pub const TALLY_ESCAPED_BAND_SPIN: usize = 114;
/// Of escaped photons, how many had crossed body i's lane at some point
/// during their flight (a photon may set more than one bit, so these are NOT
/// mutually exclusive with each other or with TALLY_ESCAPED_BAND -- they
/// answer "crossed body i" independently for i=0..3, per doc §8's explicit
/// "keep crossed body i for i=0..3"). [crossed_body*2 + species], len 8.
pub const TALLY_ESCAPED_CROSSED_BODY: usize = 144;
/// Spin sum (fixed-point) over the same index space, len 8.
pub const TALLY_ESCAPED_CROSSED_BODY_SPIN: usize = 152;
/// Counter-stream survival (doc §8, "measured not integrated"): whenever a
/// photon that crossed body `i`'s lane at some earlier point arrives at
/// (hits, whether absorbed or crossing again) a DIFFERENT body `j`, tally
/// whether FLAG_COLLIDED_SINCE_CROSS is still clear (survived uninterrupted)
/// or set (scattered somewhere in between). [i*4*2 + j*2 + collided(0/1)],
/// len 4*4*2 = 32. Diagonal i==j entries are never written (a photon cannot
/// "arrive at" the body whose lane it just crossed without leaving first).
pub const TALLY_SURVIVAL: usize = 160;
pub const TALLY_LEN: usize = 192;

pub const HIST_BINS: usize = 36; // 5-degree bins over [-90, 90]
pub const HIST_MAX_BODIES: usize = 4;
pub const HIST_EMISSION_BLOCK: usize = 0;
pub const HIST_ARRIVAL_BLOCK: usize = HIST_MAX_BODIES * HIST_BINS;
pub const HIST_LEN: usize = HIST_MAX_BODIES * HIST_BINS * 2;

pub fn hist_bin(lat_deg: f64) -> usize {
    (((lat_deg + 90.0) / 5.0).floor() as isize).clamp(0, HIST_BINS as isize - 1) as usize
}

/// Escape/exit latitude band (doc §8): 0 = polar (>60 deg), 1 = mid
/// (30-60 deg), 2 = equatorial (<30 deg), of the ABSOLUTE latitude.
pub fn escape_band(abs_lat_deg: f64) -> usize {
    if abs_lat_deg > 60.0 {
        0
    } else if abs_lat_deg >= 30.0 {
        1
    } else {
        2
    }
}
