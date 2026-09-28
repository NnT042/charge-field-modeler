//! sigma calibration (doc §11c work item 4): pick the far-field cross-section
//! sigma so a single body's own measured TAU_in (ambient photons scattering
//! off ITS emitted density, near zone) matches a target -- the DSMC analogue
//! of Phase-1's `tau_in_integral` calibration in `recycling.rs`.

use glam::DVec3;

use crate::gpu::context::GpuContext;

use super::config::DsmcConfig;
use super::density::{DensityMap, KIND_AMBIENT, KIND_EMITTED, SPECIES_MINUS, SPECIES_PLUS};
use super::engine::DsmcEngine;

pub struct SigmaCalibration {
    pub sigma: f64,
    /// sigma * integral(n_emitted ds) along the equatorial b=0 line,
    /// r = R_MAX -> R_MIN (Phase-1's TAU_IN_R_MIN/MAX, 1.05 and B_MAX).
    /// Reported for the `Constant`-equivalent kernel (no closing-speed
    /// factor) per doc §11c work item 4's explicit "simpler and acceptable"
    /// fallback -- NOT rescaled for `Catch`'s mean closing speed.
    pub tau_in_measured: f64,
    /// sigma * integral(n_ambient ds), equatorial, r = 1 -> R_D (an emitted
    /// photon's own optical depth outbound through the ambient gas).
    pub tau_out_measured: f64,
    /// 1 / (sigma * n_ambient): the far-field ambient mean free path.
    pub lambda_far: f64,
    pub passes: u32,
    pub converged: bool,
}

fn perp_axis(axis: DVec3) -> DVec3 {
    let helper = if axis.x.abs() < 0.9 { DVec3::X } else { DVec3::Y };
    (helper - axis * helper.dot(axis)).normalize()
}

fn integrate_both_species(map: &DensityMap, p0: DVec3, p1: DVec3, n: usize, kind: usize) -> f64 {
    map.integrate_ray(p0, p1, n, kind, SPECIES_MINUS) + map.integrate_ray(p0, p1, n, kind, SPECIES_PLUS)
}

/// Doc's algorithm: run single-body steady state at a guess sigma, measure
/// TAU_in from the emitted-density map, rescale, repeat <= 3 passes until
/// within 5%. `warmup_steps`/`avg_steps` are the same knobs `run_steady`
/// takes -- keep them small enough that 3 calibration passes stay inside the
/// report's few-minutes budget (doc §11c work item 5).
pub fn calibrate_sigma_for_tau_in(
    ctx: &GpuContext,
    base_config: &DsmcConfig,
    tau_in_target: f64,
    warmup_steps: u32,
    avg_steps: u32,
) -> SigmaCalibration {
    assert_eq!(base_config.bodies.len(), 1, "sigma calibration is single-body only (doc §11c work item 4)");
    let body = base_config.bodies[0];
    let r_d = match base_config.boundary {
        super::config::Boundary::SphereEscape { r_d } => r_d,
        super::config::Boundary::PeriodicBox { .. } => panic!("sigma calibration needs a SphereEscape domain"),
    };

    const R_MIN: f64 = 1.05; // Phase-1 TAU_IN_R_MIN
    let r_max = 6.0_f64.min(r_d * 0.9); // Phase-1 TAU_IN_R_MAX = B_MAX, capped under R_D

    let axis = body.orient.spin_axis();
    let equatorial_dir = perp_axis(axis);
    let center = body.center;

    // Crude order-of-magnitude first guess: treat the near-zone emitted
    // density as if it were of order n_ambient over the R_MIN..R_MAX span
    // (refined by the rescaling loop below regardless of how rough this is).
    let mut sigma = (tau_in_target / (base_config.n_ambient.max(1e-6) * (r_max - R_MIN))).max(1e-6);

    let mut tau_in_measured = 0.0;
    let mut tau_out_measured = 0.0;
    let mut passes = 0u32;
    let mut converged = false;

    for pass in 1..=3u32 {
        passes = pass;
        let mut cfg = base_config.clone();
        cfg.sigma = sigma;
        let mut engine = DsmcEngine::new(ctx, cfg);
        let report = engine.run_steady(ctx, warmup_steps, avg_steps, avg_steps.max(1));
        let map = DensityMap::new(&engine.grid, engine.free_volume_cpu().to_vec(), report.density_raw.clone(), report.avg_steps);
        map.assert_no_overflow();

        let p_max = center + equatorial_dir * r_max;
        let p_min = center + equatorial_dir * R_MIN;
        let n_e_integral = integrate_both_species(&map, p_max, p_min, 64, KIND_EMITTED);
        tau_in_measured = sigma * n_e_integral;

        let p_surface = center + equatorial_dir * 1.0;
        let p_edge = center + equatorial_dir * r_d;
        let n_a_integral = integrate_both_species(&map, p_surface, p_edge, 64, KIND_AMBIENT);
        tau_out_measured = sigma * n_a_integral;

        if tau_in_measured < 1e-12 {
            // No emitted density reached the sample line at all (sigma too
            // small / run too short to build up a profile) -- widen sigma a
            // lot and try again rather than dividing by ~0.
            sigma *= 10.0;
            continue;
        }
        let ratio = tau_in_measured / tau_in_target;
        if (ratio - 1.0).abs() < 0.05 {
            converged = true;
            break;
        }
        sigma /= ratio.max(1e-6);
    }

    let lambda_far = 1.0 / (sigma * base_config.n_ambient.max(1e-12));
    SigmaCalibration { sigma, tau_in_measured, tau_out_measured, lambda_far, passes, converged }
}

/// TAU_in / TAU_out as OUTCOMES of a run at a given sigma (the B2b lesson:
/// in a pump-mode, reciprocal engine an isolated body's exhaust is thinner
/// than the ambient field, so Phase-1's "TAU_in = TAU_out = 3" cannot be a
/// calibration target -- forcing TAU_in = 3 drags TAU_out past 100). Both
/// are averaged over `n_az` equatorial azimuths to beat single-ray noise.
/// Same definitions as `calibrate_sigma_for_tau_in` (Constant-equivalent,
/// no closing-speed factor): TAU_in = sigma * integral(n_emitted) over
/// r = 6 -> 1.05; TAU_out = sigma * integral(n_ambient) over r = 1 -> R_D.
pub fn measure_taus(map: &DensityMap, config: &DsmcConfig, sigma: f64, n_az: usize) -> (f64, f64) {
    assert_eq!(config.bodies.len(), 1, "measure_taus is single-body only");
    let body = config.bodies[0];
    let r_d = match config.boundary {
        super::config::Boundary::SphereEscape { r_d } => r_d,
        super::config::Boundary::PeriodicBox { .. } => panic!("measure_taus needs a SphereEscape domain"),
    };
    let axis = body.orient.spin_axis();
    let e1 = perp_axis(axis);
    let e2 = axis.cross(e1);
    let (mut tin, mut tout) = (0.0, 0.0);
    for k in 0..n_az {
        let phi = std::f64::consts::TAU * (k as f64 + 0.5) / n_az as f64;
        let dir = e1 * phi.cos() + e2 * phi.sin();
        let c = body.center;
        tin += integrate_both_species(map, c + dir * 6.0_f64.min(r_d * 0.9), c + dir * 1.05, 64, KIND_EMITTED);
        tout += integrate_both_species(map, c + dir * 1.0, c + dir * r_d, 96, KIND_AMBIENT);
    }
    (sigma * tin / n_az as f64, sigma * tout / n_az as f64)
}
