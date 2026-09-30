//! Pod/Zeroable GPU-layout structs shared by host and WGSL (doc §6).

use bytemuck::{Pod, Zeroable};

/// 112-byte photon record, 7 x vec4 (doc §6, grown for §13).
///
/// `a` = pos.xyz, bitcast(id) as f32 in `.w`.
/// `b` = dir.xyz (unit), scalar spin in `.w`. Under `RoughSphere` the scalar
///       is the helicity (w.dir)*R/spin_edge0, which is +-1 at birth.
/// `c` = [flags, bounces, passage state (doc §14 shells; 0 otherwise),
///       stiffness/stack level k] (u32 x4).
/// `d` = spin vector w.xyz (edge speed |w|*R_k in c, with R_k = 2^k), and
///       speed |v| in c in `.w` (1.0 except under `RoughSphere`+`Free`).
/// `e`, `f`, `g` = spin vectors of inner stacked levels 0, 1, 2 (doc §13.7;
///       each locked at edge c, zero above the outer level k). Under stacking
///       the scalar in `b.w` is the helicity of level 0 (the axial spin).
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
    pub d: [f32; 4],
    pub e: [f32; 4],
    pub f: [f32; 4],
    pub g: [f32; 4],
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
            d: [0.0; 4],
            e: [0.0; 4],
            f: [0.0; 4],
            g: [0.0; 4],
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

    /// E1 (doc §11d): per-photon stiffness level, stored whole in the spare
    /// word `c[3]` (never touched by any flag bit -- see the `flags` doc
    /// above, which only covers `c[0]`).
    pub fn stiffness_level(&self) -> u32 {
        self.c[3]
    }

    /// m = 2^k (bottom.pdf: a stacked spin doubles collision mass).
    pub fn stiffness_mass(&self) -> f64 {
        (2.0f64).powi(self.c[3] as i32)
    }

    /// Doc §13: spin vector w (radius R_k = 2^k = m).
    pub fn omega(&self) -> [f64; 3] {
        [self.d[0] as f64, self.d[1] as f64, self.d[2] as f64]
    }

    pub fn speed(&self) -> f64 {
        self.d[3] as f64
    }

    /// Velocity dir*|v| in c.
    pub fn velocity(&self) -> [f64; 3] {
        let s = self.speed();
        [self.b[0] as f64 * s, self.b[1] as f64 * s, self.b[2] as f64 * s]
    }

    /// Edge speed |w|*R_k in c.
    pub fn edge_speed(&self) -> f64 {
        let w = self.omega();
        (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]).sqrt() * self.stiffness_mass()
    }

    /// Spin vector of level `lvl` (0..=k): inner levels from e/f/g, the outer
    /// level from d.
    pub fn level_omega(&self, lvl: u32) -> [f64; 3] {
        let v = if lvl == self.c[3] {
            self.d
        } else {
            match lvl {
                0 => self.e,
                1 => self.f,
                _ => self.g,
            }
        };
        [v[0] as f64, v[1] as f64, v[2] as f64]
    }

    /// Stacking mass model (doc §13.7): (linear, rotational summed over all
    /// levels), with unit linear mass and I = kappa*R_level^2.
    pub fn energies_stacked(&self, kappa: f64) -> (f64, f64) {
        let mut rot = 0.0;
        for lvl in 0..=self.c[3].min(3) {
            let w = self.level_omega(lvl);
            let r = (2.0f64).powi(lvl as i32);
            rot += 0.5 * kappa * r * r * (w[0] * w[0] + w[1] * w[1] + w[2] * w[2]);
        }
        (0.5 * self.speed() * self.speed(), rot)
    }

    /// (linear, rotational) kinetic energy with m = R = 2^k, I = kappa*m*R^2.
    pub fn energies(&self, kappa: f64) -> (f64, f64) {
        let m = self.stiffness_mass();
        let w = self.omega();
        let w2 = w[0] * w[0] + w[1] * w[1] + w[2] * w[2];
        (0.5 * m * self.speed() * self.speed(), 0.5 * kappa * m * m * m * w2)
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
    /// E1 (doc §11d): k0, the stiffness level emitted photons are born with,
    /// and (doubling as) the NTC size-scaling ceiling (see `DsmcConfig`).
    pub emit_stiffness_level: u32,
    /// E1: 0 = StiffnessShed::None, 1 = Transfer.
    pub stiffness_shed: u32,
    /// E1: 0/1 -- arp.pdf size-driven cross-section on/off.
    pub sigma_size_scaling: u32,

    /// E1 gate-only negative control: 0/1 -- force the equal-mass exchange
    /// branch regardless of k1/k2 (see `DsmcConfig::debug_force_equal_mass_exchange`).
    pub debug_force_equal_mass_exchange: u32,
    /// NTC bound on |v1-v2| for the Catch kernel (2.0 unless RoughSphere+Free).
    pub g_max_catch: f32,
    /// Doc §13 contact law (`RoughParams`).
    pub rough_beta: f32,
    pub rough_beta_twist: f32,

    pub rough_kappa: f32,
    pub twist_loss_to_normal: u32,
    /// D7: 0 = Free, 1 = PinnedC.
    pub speed_mode: u32,
    pub spin_edge0: f32,

    /// 0/1: skew rough-sphere impact parameters toward grazing (arp.pdf).
    pub glancing_bias: u32,
    /// Doc §13.7: 0/1 stacking on (unit linear mass), highest level, unstack
    /// threshold (edge, in c).
    pub stacking: u32,
    pub stack_max_level: u32,
    pub unstack_edge: f32,

    /// Doc §14: 0 = Absorber bodies (legacy), 1 = permeable Shells.
    pub body_model: u32,
    /// Shells: body-grain contact rate per unit relative path length.
    pub tau_body: f32,
    /// §14.8: body-grain spin edge speed (0 = non-spinning grains).
    pub grain_edge: f32,
    /// Body-contact law (defaults to the photon-photon law).
    pub body_twist_loss: u32,

    pub body_beta: f32,
    pub body_beta_twist: f32,
    pub _pad14a: u32,
    pub _pad14b: u32,
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

// --- E1: emitted-photon stiffness (doc §11d) ---

/// Of escaped photons that were EMITTED (not ambient), sum of their
/// stiffness level k (plain integer sum, k in 0..7 -- no fixed-point scale
/// needed) -- the report's "mean k of escaping emitted photons".
pub const TALLY_ESCAPED_EMITTED_K_SUM: usize = 192;
pub const TALLY_ESCAPED_EMITTED_K_COUNT: usize = 193;
/// Per-photon deflection angle (radians, fixed-point TALLY_SCALE) from the
/// exchange step alone, bucketed by that photon's OWN pre-collision
/// stiffness level k (0..7) -- e1_heavy_deflects_less compares bucket 0 vs
/// bucket 3 (or any k1 vs k2). Every collision contributes two samples (one
/// per photon), tallied regardless of exchange rule (SwapRenorm collisions
/// land in bucket 0 since ambient/unshed photons stay at k=0).
pub const TALLY_DEFLECT_SUM: usize = 194; // len 8
pub const TALLY_DEFLECT_SUMSQ: usize = 202; // len 8
pub const TALLY_DEFLECT_COUNT: usize = 210; // len 8 (plain count)

// --- Doc §13: rough-sphere photons ---

/// Accepted-candidate draws whose |v1-v2| exceeded `g_max_catch` (accept
/// probability clamped at 1, which biases the rate). This must stay at zero,
/// and the host prints a loud warning when it doesn't.
pub const TALLY_ACCEPT_OVERFLOW: usize = 218;
/// PinnedC: photons whose post-collision linear energy change could not be
/// booked into their spin (rotational energy would go negative). Count, and
/// the unbooked energy (fixed-point).
pub const TALLY_PIN_DEFICIT_COUNT: usize = 219;
pub const TALLY_PIN_DEFICIT_SUM: usize = 220;
/// Post-collision photon updates with edge speed |w|R > c.
pub const TALLY_EDGE_OVER_C: usize = 221;
/// Every escape (all kinds/species): count, speed sum and sum of squares,
/// edge-speed sum and sum of squares (fixed-point), and edge > c count.
pub const TALLY_ESC_ALL: usize = 222;
pub const TALLY_ESC_SPEED_SUM: usize = 223;
pub const TALLY_ESC_SPEED_SUMSQ: usize = 224;
pub const TALLY_ESC_EDGE_SUM: usize = 225;
pub const TALLY_ESC_EDGE_SUMSQ: usize = 226;
pub const TALLY_ESC_EDGE_OVER_C: usize = 227;
/// Energy moved from twist grip into normal separation (fixed-point).
pub const TALLY_TWIST_TO_NORMAL: usize = 228;

/// Doc §13.7 stacking: stack events (levels added), unstack events, and
/// photons left above edge c at the top level.
pub const TALLY_STACK_EVENTS: usize = 232;
pub const TALLY_UNSTACK_EVENTS: usize = 233;
pub const TALLY_STACK_CAPPED: usize = 234;
/// Escapes by outer level k (0..3), all kinds.
pub const TALLY_ESC_LEVEL: usize = 235;

// --- Doc §14: permeable shell bodies (per body, 4 bodies max) ---
pub const TALLY_SH_ENTRIES: usize = 240;
/// Passages that ended after >= 1 body contact (the body's "emission").
pub const TALLY_SH_EXITS: usize = 244;
/// Passages entering above +60 deg and leaving below -60 deg (or vice versa).
pub const TALLY_SH_THROUGH: usize = 248;
pub const TALLY_SH_CONTACTS: usize = 252;
/// Passages with no contact (straight through, untouched).
pub const TALLY_SH_CLEAN: usize = 256;
/// Contacted exits by body x |latitude| band (0 polar, 1 mid, 2 equatorial).
pub const TALLY_SH_EXIT_BAND: usize = 260;
pub const TALLY_SH_EXIT_BAND_HEL: usize = 272;
/// Momentum delivered to the body [body*3 + axis] (fixed-point).
pub const TALLY_SH_DP: usize = 284;
/// Torque from delivered linear momentum, (x - center) x dp (fixed-point).
pub const TALLY_SH_TORQUE: usize = 296;
/// Linear energy given to photons by the driven body (fixed-point).
pub const TALLY_SH_WORK: usize = 308;

// --- Doc §15: instantaneous ball ---
/// Push-outs: a photon found touching or inside a ball that is still closing
/// on it is moved back to the surface and takes the contact there. Happens
/// when a ball faster than c re-catches a photon; reported every run.
pub const TALLY_BALL_OVERLAP: usize = 312;
/// Body-0 passages by entry angle to the pole (10 deg bins, 0..8) x outcome
/// (0 clean, 1 pp-only, 2 ball hit first, 3 ball hit after pp) x species:
/// [(abin*4 + outcome)*2 + spec], len 72.
pub const TALLY_PASS: usize = 313;
/// Body-0 opposite-pole exits by entry angle bin x species, len 18.
pub const TALLY_PASS_THROUGH: usize = 385;
/// Body-0 exit helicity sum by outcome x species (fixed-point), len 8.
pub const TALLY_PASS_HEL: usize = 403;

pub const TALLY_LEN: usize = 412;

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
