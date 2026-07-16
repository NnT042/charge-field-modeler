//! Dev utility: render the three baked calibration loops to an SVG so the hitbox
//! geometry can be eyeballed without launching Godot. Not part of the app.
//!
//!   cargo run --release --bin dump_loops -- <out.svg>
//!
//! Layout: 3 rows (proton / neutron / electron) × 2 cols (pole-view, side-view).
//! Each cell draws the recorded loop bright, plus the loop at several swing
//! phases faint — the swept disc (proton/electron) or kite crescent (neutron).

use charge_field_modeler::hitbox::{bake_loop, recommended_samples, BakedLoop};
use glam::DVec3;

const CELL: f64 = 420.0;
const PAD: f64 = 30.0;
const PHASES: usize = 16;

struct Row {
    name: &'static str,
    loop_level: u8,
    color: &'static str,
}

fn project(p: DVec3, view: usize) -> (f64, f64) {
    // view 0 = pole-view (look down Z): (x, y); view 1 = side-view: (x, z).
    match view {
        0 => (p.x, p.y),
        _ => (p.x, p.z),
    }
}

/// Build an SVG polyline `points=""` string, fitting `pts` into a CELL box at
/// origin (ox, oy) using a shared scale/center so all views of one particle
/// share proportions.
fn polyline(pts: &[DVec3], view: usize, ox: f64, oy: f64, span: f64, cx: f64, cy: f64) -> String {
    let s = (CELL - 2.0 * PAD) / span.max(1e-9);
    let half = CELL / 2.0;
    let mut out = String::with_capacity(pts.len() * 12);
    for p in pts {
        let (u, v) = project(*p, view);
        let x = ox + half + (u - cx) * s;
        let y = oy + half - (v - cy) * s; // flip y for SVG
        out.push_str(&format!("{:.1},{:.1} ", x, y));
    }
    out
}

/// Extent (span, center-u, center-v) across BOTH views so the particle keeps
/// aspect. We size on the pole-view radius (the disc face).
fn extent(baked: &BakedLoop) -> (f64, f64, f64, f64, f64) {
    let mut mn = [f64::INFINITY; 3];
    let mut mx = [f64::NEG_INFINITY; 3];
    // Sample the swept family to capture the full occupied region.
    for k in 0..PHASES {
        let ang = k as f64 / PHASES as f64 * std::f64::consts::TAU;
        for p in baked.swept_points(ang) {
            for a in 0..3 {
                let c = p[a];
                if c < mn[a] { mn[a] = c; }
                if c > mx[a] { mx[a] = c; }
            }
        }
    }
    let span = (0..3).map(|a| mx[a] - mn[a]).fold(0.0, f64::max);
    let cx = 0.5 * (mn[0] + mx[0]);
    let cy_pole = 0.5 * (mn[1] + mx[1]);
    let cy_side = 0.5 * (mn[2] + mx[2]);
    (span, cx, cy_pole, cy_side, span)
}

fn subsample(pts: Vec<DVec3>, target: usize) -> Vec<DVec3> {
    if pts.len() <= target || target == 0 {
        return pts;
    }
    let step = pts.len() / target;
    pts.into_iter().step_by(step.max(1)).collect()
}

fn main() {
    let out_path = std::env::args().nth(1).unwrap_or_else(|| "loops.svg".into());

    let rows = [
        Row { name: "Proton  (L12 loop -> L13 precession: sweeps disc)", loop_level: 12, color: "#f2dc8a" },
        Row { name: "Neutron (L11 loop -> L12 orbital kite: rakes poles)", loop_level: 11, color: "#7cf0e6" },
        Row { name: "Electron (L8 loop -> L9 precession)", loop_level: 8, color: "#c9a6ff" },
    ];

    let width = PAD + 3.0 * CELL + PAD;
    let height = PAD + 3.0 * (CELL + 40.0) + PAD;
    let mut svg = String::new();
    svg.push_str(&format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{:.0}\" height=\"{:.0}\" \
         viewBox=\"0 0 {:.0} {:.0}\" style=\"background:#080a10\">\n",
        width, height, width, height
    ));
    svg.push_str("<style>text{font-family:monospace;fill:#aeb8c2}</style>\n");
    let col_titles = ["pole view (down the spin axis)", "side view (edge-on)"];

    for (ri, row) in rows.iter().enumerate() {
        let samples = recommended_samples(row.loop_level).min(3200);
        let baked = bake_loop(row.loop_level, samples);
        let (span, cx, cy_pole, cy_side, _) = extent(&baked);

        let row_y = PAD + ri as f64 * (CELL + 40.0);
        svg.push_str(&format!(
            "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"16\" fill=\"{}\">{}</text>\n",
            PAD, row_y + 18.0, row.color, row.name
        ));

        for view in 0..2usize {
            let ox = PAD + view as f64 * CELL;
            let oy = row_y + 24.0;
            let cyc = if view == 0 { cy_pole } else { cy_side };
            // cell frame + title
            svg.push_str(&format!(
                "<rect x=\"{:.0}\" y=\"{:.0}\" width=\"{:.0}\" height=\"{:.0}\" \
                 fill=\"none\" stroke=\"#1b2230\"/>\n",
                ox, oy, CELL, CELL
            ));
            svg.push_str(&format!(
                "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"11\">{}</text>\n",
                ox + 6.0, oy + 14.0, col_titles[view]
            ));

            // faint swept family
            for k in 0..PHASES {
                let ang = k as f64 / PHASES as f64 * std::f64::consts::TAU;
                let fam = subsample(baked.swept_points(ang), 400);
                let pl = polyline(&fam, view, ox, oy, span, cx, cyc);
                svg.push_str(&format!(
                    "<polyline points=\"{}\" fill=\"none\" stroke=\"{}\" \
                     stroke-width=\"0.5\" opacity=\"0.12\"/>\n",
                    pl.trim(), row.color
                ));
            }
            // bright recorded loop (phase 0)
            let base = baked.swept_points(0.0);
            let pl = polyline(&base, view, ox, oy, span, cx, cyc);
            svg.push_str(&format!(
                "<polyline points=\"{}\" fill=\"none\" stroke=\"{}\" \
                 stroke-width=\"0.9\" opacity=\"0.95\"/>\n",
                pl.trim(), row.color
            ));
        }

        // Third column: orientation histogram — the probability the base
        // particle meets the field upright / sideways / upside-down over the loop.
        {
            let ox = PAD + 2.0 * CELL;
            let oy = row_y + 24.0;
            svg.push_str(&format!(
                "<rect x=\"{:.0}\" y=\"{:.0}\" width=\"{:.0}\" height=\"{:.0}\" \
                 fill=\"none\" stroke=\"#1b2230\"/>\n",
                ox, oy, CELL, CELL
            ));
            svg.push_str(&format!(
                "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"11\">orientation vs pole (loop probability)</text>\n",
                ox + 6.0, oy + 14.0
            ));
            let hist = baked.orientation_histogram(30);
            let maxh = hist.iter().cloned().fold(0.0_f64, f64::max).max(1e-9);
            let inner = CELL - 2.0 * PAD;
            let bw = inner / hist.len() as f64;
            let base_y = oy + CELL - PAD;
            for (bi, &v) in hist.iter().enumerate() {
                let bh = (v / maxh) * inner;
                let x = ox + PAD + bi as f64 * bw;
                svg.push_str(&format!(
                    "<rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" \
                     fill=\"{}\" opacity=\"0.85\"/>\n",
                    x, base_y - bh, bw * 0.85, bh, row.color
                ));
            }
            svg.push_str(&format!(
                "<line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#33404f\"/>\n",
                ox + PAD, base_y, ox + PAD + inner, base_y
            ));
            svg.push_str(&format!(
                "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"10\">upside-down</text>\n",
                ox + PAD, base_y + 15.0
            ));
            svg.push_str(&format!(
                "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"10\" text-anchor=\"middle\">sideways</text>\n",
                ox + CELL / 2.0, base_y + 15.0
            ));
            svg.push_str(&format!(
                "<text x=\"{:.0}\" y=\"{:.0}\" font-size=\"10\" text-anchor=\"end\">upright</text>\n",
                ox + PAD + inner, base_y + 15.0
            ));
        }
    }

    svg.push_str("</svg>\n");
    std::fs::write(&out_path, svg).expect("write svg");
    eprintln!("wrote {}", out_path);
}
