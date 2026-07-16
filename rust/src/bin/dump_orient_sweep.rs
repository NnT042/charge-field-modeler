//! Dev diagnostic (measurement, no rule change): sweep a proton's outer_spin and
//! measure the hit-weighted presented-orientation distribution for three field
//! geometries — perpendicular, axial, isotropic. Answers "does the hit-weighting
//! break the loop's orientation symmetry as the disc spins up, and where?"
//!
//!   cargo run --release --bin dump_orient_sweep -- <out.svg>
//!
//! `mean_cos` = Σ h·cos(axis-to-pole): 0 = symmetric (upright ≡ upside-down),
//! nonzero = a net orientation bias the collisions actually sample.

use charge_field_modeler::calibration::CalibrationParticle;
use charge_field_modeler::hitbox::bake_loop;
use glam::DVec3;

const BINS: usize = 30;
const CELL_W: f64 = 300.0;
const CELL_H: f64 = 180.0;
const PAD: f64 = 34.0;
const ROW_GAP: f64 = 30.0;

fn mean_cos(h: &[f64]) -> f64 {
    let bins = h.len();
    h.iter()
        .enumerate()
        .map(|(i, &v)| {
            let c = -1.0 + (i as f64 + 0.5) / bins as f64 * 2.0;
            v * c
        })
        .sum()
}

fn bars_svg(h: &[f64], ox: f64, oy: f64, color: &str, title: &str) -> String {
    let maxh = h.iter().cloned().fold(0.0_f64, f64::max).max(1e-9);
    let inner_w = CELL_W - 2.0 * PAD;
    let inner_h = CELL_H - 2.0 * PAD;
    let bw = inner_w / h.len() as f64;
    let base_y = oy + CELL_H - PAD;
    let mut s = String::new();
    s.push_str(&format!(
        "<rect x=\"{:.0}\" y=\"{:.0}\" width=\"{:.0}\" height=\"{:.0}\" fill=\"none\" stroke=\"#1b2230\"/>\n",
        ox, oy, CELL_W, CELL_H
    ));
    s.push_str(&format!(
        "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"11\" fill=\"#aeb8c2\">{}</text>\n",
        ox + 6.0, oy + 13.0, title
    ));
    for (i, &v) in h.iter().enumerate() {
        let bh = (v / maxh) * inner_h;
        let x = ox + PAD + i as f64 * bw;
        s.push_str(&format!(
            "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"{}\" opacity=\"0.85\"/>\n",
            x, base_y - bh, bw * 0.85, bh, color
        ));
    }
    // center line (sideways) + baseline
    let cx = ox + PAD + inner_w * 0.5;
    s.push_str(&format!(
        "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#33404f\" stroke-dasharray=\"3 3\"/>\n",
        cx, oy + PAD, cx, base_y
    ));
    s
}

fn main() {
    let out = std::env::args().nth(1).unwrap_or_else(|| "orient_sweep.svg".into());

    // Proton loop (downsampled for the diagnostic).
    let mut p = CalibrationParticle::new(bake_loop(12, 1024), 1.0);

    let spins = [0.0, 0.25, 0.5, 0.75, 0.95];
    // Field geometries: perpendicular (+X, ⊥ pole Z), axial (+Z, ∥ pole), isotropic.
    let fields: [(&str, Option<DVec3>, &str); 3] = [
        ("perp field (+X)", Some(DVec3::X), "#f2dc8a"),
        ("axial field (+Z)", Some(DVec3::Z), "#8ad0f2"),
        ("isotropic", None, "#c9a6ff"),
    ];

    let width = PAD + 3.0 * CELL_W + PAD;
    let height = PAD + spins.len() as f64 * (CELL_H + ROW_GAP) + PAD;
    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" \
         viewBox=\"0 0 {:.0} {:.0}\" style=\"background:#080a10\">\n",
        width, height, width, height
    ));
    svg.push_str("<style>text{font-family:monospace}</style>\n");

    eprintln!("proton | orientation mean_cos (0 = symmetric) | swing_drag (< 0 = spin-DOWN)");
    eprintln!("outer_spin |   perp(+X) |  axial(+Z) | isotropic  |  swing_drag");
    for (ri, &os) in spins.iter().enumerate() {
        p.outer_spin = os;
        let drag = p.swing_drag_torque(64, 24, false);
        let row_y = PAD + ri as f64 * (CELL_H + ROW_GAP);
        svg.push_str(&format!(
            "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"13\" fill=\"#dfe6ec\">\
             outer_spin = {:.2} c    swing_drag = {:+.4}  (&lt;0 = spin-down)</text>\n",
            PAD, row_y - 6.0, os, drag
        ));
        let mut mc = [0.0; 3];
        for (ci, (label, dir, color)) in fields.iter().enumerate() {
            let h = match dir {
                Some(d) => p.orientation_participation(Some(*d), BINS, 48, 1),
                None => p.orientation_participation(None, BINS, 24, 48),
            };
            mc[ci] = mean_cos(&h);
            let ox = PAD + ci as f64 * CELL_W;
            let title = format!("{}   mean_cos={:+.4}", label, mc[ci]);
            svg.push_str(&bars_svg(&h, ox, row_y, color, &title));
        }
        eprintln!(
            "   {:.2}     | {:+.5} | {:+.5} | {:+.5} | {:+.5}",
            os, mc[0], mc[1], mc[2], drag
        );
    }

    // Bottom axis legend on the last row.
    let last_y = PAD + (spins.len() - 1) as f64 * (CELL_H + ROW_GAP) + CELL_H - PAD + 14.0;
    for ci in 0..3 {
        let ox = PAD + ci as f64 * CELL_W;
        svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"9\" fill=\"#8a97a5\">upside-down</text>\n", ox + PAD, last_y));
        svg.push_str(&format!("<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"9\" fill=\"#8a97a5\" text-anchor=\"end\">upright</text>\n", ox + CELL_W - PAD, last_y));
    }

    svg.push_str("</svg>\n");
    std::fs::write(&out, svg).expect("write svg");
    eprintln!("wrote {}", out);
}
