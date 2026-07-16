//! Dev diagnostic (measurement only): the velocity-channel spin-down drag vs
//! outer_spin, for proton / neutron / electron, over a fine grid up to the c
//! ceiling. The drag is chirality/composition-independent (pure momentum), so
//! this curve — together with the analytic chirality pump (flat in outer_spin,
//! ∝ field imbalance) — fixes the equilibrium: pump = |drag| ⇒ a finite rate.
//!
//!   cargo run --release --bin dump_drag_curves -- <out.svg>
//!
//! Prints a raw table (torque units) and a per-particle-normalized shape, then
//! renders the normalized shapes so proton/neutron/electron overlay.

use charge_field_modeler::calibration::CalibrationParticle;
use charge_field_modeler::hitbox::bake_loop;

const W: f64 = 720.0;
const H: f64 = 460.0;
const M: f64 = 60.0;

struct Series {
    name: &'static str,
    loop_level: u8,
    color: &'static str,
    drag: Vec<f64>,
}

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "drag_curves.svg".into());
    let spins: Vec<f64> = (0..=10).map(|k| k as f64 * 0.1).map(|s| s.min(0.99)).collect();

    let mut series = [
        Series { name: "proton  (L12)", loop_level: 12, color: "#f2dc8a", drag: vec![] },
        Series { name: "neutron (L11)", loop_level: 11, color: "#7cf0e6", drag: vec![] },
        Series { name: "electron (L8)", loop_level: 8, color: "#c9a6ff", drag: vec![] },
    ];

    for s in series.iter_mut() {
        let mut p = CalibrationParticle::new(bake_loop(s.loop_level, 768), 1.0);
        for &os in &spins {
            p.outer_spin = os;
            s.drag.push(p.swing_drag_torque(48, 16));
        }
    }

    // Raw + normalized tables to stderr (cache-able).
    eprint!("outer_spin ");
    for &os in &spins {
        eprint!("| {:>5.2} ", os);
    }
    eprintln!();
    for s in &series {
        eprint!("{:>13} raw ", s.name);
        for d in &s.drag {
            eprint!("|{:>6.1} ", d);
        }
        eprintln!();
        let denom = s.drag.last().cloned().unwrap_or(-1.0).abs().max(1e-9);
        eprint!("{:>13} nrm ", "");
        for d in &s.drag {
            eprint!("|{:>6.3} ", d / denom);
        }
        eprintln!();
    }

    // Normalized-shape line plot.
    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" \
         viewBox=\"0 0 {:.0} {:.0}\" style=\"background:#080a10\">\n<style>text{{font-family:monospace;fill:#aeb8c2}}</style>\n",
        W, H, W, H
    ));
    // axes (y goes 0 at top to -1 at bottom, i.e. spin-down downward)
    let x = |t: f64| M + t * (W - 2.0 * M);
    let y = |v: f64| M + (-v) * (H - 2.0 * M); // v in [-1, 0]
    svg.push_str(&format!("<line x1=\"{:.0}\" y1=\"{:.0}\" x2=\"{:.0}\" y2=\"{:.0}\" stroke=\"#33404f\"/>\n", x(0.0), y(0.0), x(1.0), y(0.0)));
    svg.push_str(&format!("<line x1=\"{:.0}\" y1=\"{:.0}\" x2=\"{:.0}\" y2=\"{:.0}\" stroke=\"#33404f\"/>\n", x(0.0), y(0.0), x(0.0), y(-1.0)));
    svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"12\">outer_spin →</text>\n", x(0.72), y(0.0) + 22.0));
    svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"12\">spin-down drag ↓ (normalized)</text>\n", M - 40.0, M - 18.0));

    for (si, s) in series.iter().enumerate() {
        let denom = s.drag.last().cloned().unwrap_or(-1.0).abs().max(1e-9);
        let pts: String = spins
            .iter()
            .zip(&s.drag)
            .map(|(&t, &d)| format!("{:.1},{:.1} ", x(t), y(d / denom)))
            .collect();
        svg.push_str(&format!(
            "<polyline points=\"{}\" fill=\"none\" stroke=\"{}\" stroke-width=\"2\"/>\n",
            pts.trim(), s.color
        ));
        svg.push_str(&format!(
            "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"13\" fill=\"{}\">{}</text>\n",
            W - M - 150.0, M + 20.0 + si as f64 * 20.0, s.color, s.name
        ));
    }
    svg.push_str("</svg>\n");
    std::fs::write(&out, svg).expect("write svg");
    eprintln!("wrote {}", out);
}
