//! Dev diagnostic (measurement only): does a directional field turn the pole?
//! Tilt a proton's pole by angle θ from a field along +Z and measure the
//! alignment torque (about the tilt axis X) and the spin-drag torque (about the
//! pole), swing-averaged and catch-weighted. Occlusion is NOT modeled, so any
//! alignment here comes purely from the relative-velocity (catch) asymmetry.
//!
//!   cargo run --release --bin dump_alignment -- <out.svg>

use charge_field_modeler::calibration::CalibrationParticle;
use charge_field_modeler::hitbox::bake_loop;
use glam::{DQuat, DVec3};

const W: f64 = 760.0;
const H: f64 = 460.0;
const M: f64 = 60.0;

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "alignment.svg".into());
    let mut p = CalibrationParticle::new(bake_loop(12, 768), 1.0);
    let field = DVec3::Z;
    let tilt_axis = DVec3::X;

    let thetas_deg: Vec<f64> = (0..=12).map(|k| k as f64 * 15.0).collect();
    let spins = [0.3_f64, 0.6, 0.9];
    let colors = ["#f2dc8a", "#7cf0e6", "#c9a6ff"];

    // curves[si] = Vec<(theta_deg, tau_align)>
    let mut curves: Vec<Vec<(f64, f64)>> = vec![vec![]; spins.len()];

    eprintln!("proton: alignment torque about X (turns pole toward/away from +Z field)");
    eprint!("theta_deg ");
    for &t in &thetas_deg {
        eprint!("|{:>5.0} ", t);
    }
    eprintln!();
    for (si, &os) in spins.iter().enumerate() {
        p.outer_spin = os;
        eprint!("spin {:.1}  ", os);
        for &td in &thetas_deg {
            let th = td.to_radians();
            // Tilt the pole by th about X (default pole = +Z).
            p.orientation = DQuat::from_axis_angle(tilt_axis, th);
            let tau = p.directional_torque(field, 32);
            let t_align = tau.dot(tilt_axis);
            curves[si].push((td, t_align));
            eprint!("|{:>5.1} ", t_align);
        }
        eprintln!();
    }
    // reset pose
    p.orientation = DQuat::IDENTITY;

    // Also report the pole/spin-drag component at 90° tilt vs a few spins.
    eprintln!("--- spin-drag (pole component) at theta=90 ---");
    for &os in &[0.0_f64, 0.3, 0.6, 0.9] {
        p.outer_spin = os;
        p.orientation = DQuat::from_axis_angle(tilt_axis, std::f64::consts::FRAC_PI_2);
        let tau = p.directional_torque(field, 32);
        eprintln!("  spin {:.1}: pole_torque = {:+.3}", os, tau.dot(p.pole_axis()));
    }

    // Plot tau_align vs theta.
    let maxabs = curves
        .iter()
        .flat_map(|c| c.iter().map(|(_, v)| v.abs()))
        .fold(1e-9, f64::max);
    let x = |t: f64| M + (t / 180.0) * (W - 2.0 * M);
    let y = |v: f64| H * 0.5 - (v / maxabs) * (H * 0.5 - M);

    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" viewBox=\"0 0 {:.0} {:.0}\" style=\"background:#080a10\">\n<style>text{{font-family:monospace;fill:#aeb8c2}}</style>\n",
        W, H, W, H
    ));
    // zero line + axes
    svg.push_str(&format!("<line x1=\"{:.0}\" y1=\"{:.0}\" x2=\"{:.0}\" y2=\"{:.0}\" stroke=\"#33404f\"/>\n", x(0.0), y(0.0), x(180.0), y(0.0)));
    svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"12\">pole tilt θ from field (deg) →</text>\n", x(70.0), H - 20.0));
    svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"12\">align torque +</text>\n", 8.0, M - 20.0));
    for td in [0.0, 90.0, 180.0] {
        svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"10\">{:.0}</text>\n", x(td) - 6.0, y(0.0) + 16.0, td));
    }
    for (si, c) in curves.iter().enumerate() {
        let pts: String = c.iter().map(|(t, v)| format!("{:.1},{:.1} ", x(*t), y(*v))).collect();
        svg.push_str(&format!("<polyline points=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"2\"/>\n", pts.trim(), colors[si]));
        svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"12\" fill=\"{}\">outer_spin {:.1}</text>\n", W - M - 120.0, M + si as f64 * 18.0, colors[si], spins[si]));
    }
    svg.push_str("</svg>\n");
    std::fs::write(&out, svg).expect("write svg");
    eprintln!("wrote {}", out);
}
