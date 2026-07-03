//! Headless scenario runner / trajectory dumper for atom-mode physics.
//!
//! Never calls Godot APIs — links the crate only for atom_core/atom_scenarios.
//!
//! Usage:
//!   atom_lab <scenario> [steps] [sample_every] [out.csv]
//!
//! Scenarios: hydrogen | protons  (h2:<0-7>, alpha, helium arrive with M2/M3)
//! Without out.csv, prints summary metrics only.
//!
//!   cargo run --release --manifest-path rust/Cargo.toml --bin atom_lab -- hydrogen 200000 100 out.csv

use charge_field_modeler::atom_core::AtomCore;
use charge_field_modeler::atom_scenarios::{
    run_pair, run_recording, spawn_hydrogen, standard_core,
};
use glam::DVec3;
use std::io::Write;

fn main() {
    let all_args: Vec<String> = std::env::args().skip(1).collect();
    // Any arg of the form key=value is a coupling override (for sweeps
    // without recompiling); the rest are positional.
    let (overrides, args): (Vec<&String>, Vec<&String>) =
        all_args.iter().partition(|a| a.contains('='));

    let scenario = args.first().map(|s| s.as_str()).unwrap_or("hydrogen");
    let steps: usize = args
        .get(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(200_000);
    let sample_every: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(100);
    let out_csv = args.get(3).cloned();

    let mut core = standard_core();
    for ov in overrides {
        let (key, value) = ov.split_once('=').unwrap();
        let value: f64 = value.parse().unwrap_or_else(|e| {
            eprintln!("bad override '{ov}': {e}");
            std::process::exit(2);
        });
        let c = &mut core.couplings;
        match key {
            "g_q" => c.g_q = value,
            "c_q" => c.c_q = value,
            "intake" => c.intake = value,
            "vortex" => c.vortex = value,
            "drag" => c.drag = value,
            "torque" => c.torque = value,
            "corot" => c.corot = value,
            "p_amb" => c.ambient_pressure = value,
            "dt" => core.dt = value,
            other => {
                eprintln!("unknown coupling '{other}'");
                std::process::exit(2);
            }
        }
    }
    let pair: Option<(usize, usize)> = match scenario {
        "hydrogen" => {
            let (p, e) = spawn_hydrogen(&mut core, DVec3::ZERO, DVec3::Y, true, 1.0);
            Some((p, e))
        }
        "protons" => {
            let pid = core.profile_id_by_name("proton").unwrap();
            core.spawn_particle(pid, DVec3::new(-3.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
            core.spawn_particle(pid, DVec3::new(3.0, 0.0, 0.0), DVec3::ZERO, DVec3::Y);
            None
        }
        other => {
            eprintln!("unknown scenario '{other}' (try: hydrogen, protons)");
            std::process::exit(2);
        }
    };

    println!("scenario={scenario} steps={steps} sample_every={sample_every}");
    println!(
        "couplings: G_q={} C_q={} I_q={} V_q={} D_q={} T_q={} corot={} P_amb={}",
        core.couplings.g_q,
        core.couplings.c_q,
        core.couplings.intake,
        core.couplings.vortex,
        core.couplings.drag,
        core.couplings.torque,
        core.couplings.corot,
        core.couplings.ambient_pressure,
    );

    match out_csv {
        Some(path) => {
            dump_csv(&mut core, steps, sample_every, pair, &path);
            println!("wrote {path}");
        }
        None => {
            if let Some((center, orbiter)) = pair {
                let m = run_pair(&mut core, center, orbiter, steps, sample_every, 0.4);
                println!("{m:#?}");
            } else {
                core.step_n(steps);
                println!("t={:.4}  KE={:.6e}", core.time, core.total_kinetic_energy());
                for (i, p) in core.particles.iter().enumerate() {
                    println!(
                        "  #{i} {} pos=({:.3},{:.3},{:.3}) |v|={:.4}",
                        core.profiles[p.profile_id].name,
                        p.position.x,
                        p.position.y,
                        p.position.z,
                        p.velocity.length()
                    );
                }
            }
        }
    }
}

fn dump_csv(
    core: &mut AtomCore,
    steps: usize,
    sample_every: usize,
    pair: Option<(usize, usize)>,
    path: &str,
) {
    let n_particles = core.particles.len();
    let rows = run_recording(core, steps, sample_every);

    let file = std::fs::File::create(path).expect("create output CSV");
    let mut w = std::io::BufWriter::new(file);

    let mut header = String::from("t");
    for i in 0..n_particles {
        header.push_str(&format!(",p{i}x,p{i}y,p{i}z,v{i}x,v{i}y,v{i}z"));
    }
    if pair.is_some() {
        header.push_str(",r,v_tan,v_rad,theta_deg");
    }
    writeln!(w, "{header}").unwrap();

    for row in &rows {
        let mut line = format!("{:.6}", row.t);
        for i in 0..n_particles {
            let p = row.positions[i];
            let v = row.velocities[i];
            line.push_str(&format!(
                ",{:.6},{:.6},{:.6},{:.6},{:.6},{:.6}",
                p.x, p.y, p.z, v.x, v.y, v.z
            ));
        }
        if let Some((center, orbiter)) = pair {
            let d = row.positions[orbiter] - row.positions[center];
            let r = d.length();
            let d_hat = d / r.max(1e-12);
            let v_rel = row.velocities[orbiter] - row.velocities[center];
            let v_rad = v_rel.dot(d_hat);
            let v_tan = (v_rel - d_hat * v_rad).length();
            // Pole axis isn't recorded per-frame; approximate θ against +Y
            // (scenario spawns keep the center pole on +Y).
            let theta = (d_hat.y.abs()).clamp(0.0, 1.0).acos().to_degrees();
            line.push_str(&format!(
                ",{:.6},{:.6},{:.6},{:.3}",
                r, v_tan, v_rad, theta
            ));
        }
        writeln!(w, "{line}").unwrap();
    }
}
