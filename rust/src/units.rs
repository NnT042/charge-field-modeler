//! SI unit conversion and EM spectrum classification.
//!
//! All simulation math runs in natural units (r=1, c=1). This module converts
//! to SI for display and classifies wavelengths into EM bands.
//!
//! Calibration anchor: a full charge photon (z-spin outermost, level 4, ω=c,
//! speed=1) has natural wavelength TAU × orbit_radius = 8π ≈ 25.13 → maps to
//! ~2000 nm (near-IR). Per Mathis ("Bbody"): charge photon peaks at ~2 μm.
//! Scale factor = 2000 / (8π) = 250/π ≈ 79.58.

use std::f64::consts::PI;

const NM_PER_NATURAL_WL: f64 = 250.0 / PI;

pub fn wavelength_nm(wavelength_natural: f64) -> f64 {
    wavelength_natural * NM_PER_NATURAL_WL
}

#[cfg(test)]
pub fn em_band(wl_nm: f64) -> &'static str {
    if wl_nm <= 0.0 {
        return "—";
    }
    match wl_nm {
        w if w < 0.01 => "gamma",
        w if w < 10.0 => "hard X-ray",
        w if w < 100.0 => "soft X-ray",
        w if w < 380.0 => "UV",
        w if w < 780.0 => "visible",
        w if w < 1_000_000.0 => "IR",
        w if w < 1_000_000_000.0 => "microwave",
        _ => "radio",
    }
}

pub fn em_band_detail(wl_nm: f64) -> &'static str {
    if wl_nm <= 0.0 {
        return "—";
    }
    match wl_nm {
        w if w < 0.01 => "gamma",
        w if w < 10.0 => "hard X-ray",
        w if w < 100.0 => "soft X-ray",
        w if w < 200.0 => "vacuum UV",
        w if w < 280.0 => "UV-C",
        w if w < 315.0 => "UV-B",
        w if w < 380.0 => "UV-A",
        w if w < 450.0 => "visible (violet)",
        w if w < 495.0 => "visible (blue)",
        w if w < 570.0 => "visible (green)",
        w if w < 590.0 => "visible (yellow)",
        w if w < 620.0 => "visible (orange)",
        w if w < 780.0 => "visible (red)",
        w if w < 1_400.0 => "near-IR",
        w if w < 3_000.0 => "short-wave IR",
        w if w < 8_000.0 => "mid-IR",
        w if w < 15_000.0 => "thermal IR",
        w if w < 50_000.0 => "far-IR",
        w if w < 300_000.0 => "sub-mm (electron scale)",
        w if w < 1_000_000.0 => "microwave (baryon scale)",
        w if w < 1_000_000_000.0 => "microwave",
        _ => "radio",
    }
}

pub fn format_wavelength_si(wl_nm: f64) -> String {
    if wl_nm <= 0.0 {
        return "—".to_string();
    }
    if wl_nm < 1.0 {
        format!("{:.1} pm", wl_nm * 1000.0)
    } else if wl_nm < 1_000.0 {
        format!("{:.1} nm", wl_nm)
    } else if wl_nm < 1_000_000.0 {
        format!("{:.2} μm", wl_nm / 1_000.0)
    } else if wl_nm < 1_000_000_000.0 {
        format!("{:.2} mm", wl_nm / 1_000_000.0)
    } else {
        format!("{:.3} m", wl_nm / 1_000_000_000.0)
    }
}

/// Approximate spectral color for a given wavelength in nm.
/// Returns (r, g, b) in [0, 1]. Outside the visible range (380–780 nm)
/// returns neutral gray (0.5, 0.5, 0.5).
///
/// Based on Dan Bruton's spectral approximation with intensity tapering
/// at the edges of human vision.
pub fn wavelength_to_rgb(wl_nm: f64) -> (f32, f32, f32) {
    if wl_nm < 380.0 || wl_nm > 780.0 {
        return (0.5, 0.5, 0.5);
    }

    let (r, g, b) = if wl_nm < 440.0 {
        (-(wl_nm - 440.0) / (440.0 - 380.0), 0.0, 1.0)
    } else if wl_nm < 490.0 {
        (0.0, (wl_nm - 440.0) / (490.0 - 440.0), 1.0)
    } else if wl_nm < 510.0 {
        (0.0, 1.0, -(wl_nm - 510.0) / (510.0 - 490.0))
    } else if wl_nm < 580.0 {
        ((wl_nm - 510.0) / (580.0 - 510.0), 1.0, 0.0)
    } else if wl_nm < 645.0 {
        (1.0, -(wl_nm - 645.0) / (645.0 - 580.0), 0.0)
    } else {
        (1.0, 0.0, 0.0)
    };

    let intensity = if wl_nm < 420.0 {
        0.3 + 0.7 * (wl_nm - 380.0) / (420.0 - 380.0)
    } else if wl_nm > 700.0 {
        0.3 + 0.7 * (780.0 - wl_nm) / (780.0 - 700.0)
    } else {
        1.0
    };

    (
        (r * intensity) as f32,
        (g * intensity) as f32,
        (b * intensity) as f32,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::TAU;

    #[test]
    fn calibration_z1_at_c() {
        // z-spin orbit_radius = 4, at c: wl_nat = TAU * 4 = 8π
        let wl_nat = TAU * 4.0;
        let nm = wavelength_nm(wl_nat);
        assert!((nm - 2000.0).abs() < 0.01, "z1 at c should be 2000 nm, got {}", nm);
    }

    #[test]
    fn calibration_x1_at_c() {
        // x-spin orbit_radius = 1, at c: wl_nat = TAU * 1 = 2π
        let wl_nat = TAU * 1.0;
        let nm = wavelength_nm(wl_nat);
        assert!((nm - 500.0).abs() < 0.01, "x1 at c should be 500 nm, got {}", nm);
    }

    #[test]
    fn calibration_x1_at_07c() {
        // x-spin orbit_radius=1, speed=1, omega=0.7: wl_nat = TAU * 1 / 0.7
        let wl_nat = TAU / 0.7;
        let nm = wavelength_nm(wl_nat);
        assert!((nm - 714.3).abs() < 1.0, "x1 at 0.7c should be ~714 nm, got {}", nm);
    }

    #[test]
    fn calibration_y1_at_c() {
        // y-spin orbit_radius = 2, at c: wl_nat = TAU * 2 = 4π
        let wl_nat = TAU * 2.0;
        let nm = wavelength_nm(wl_nat);
        assert!((nm - 1000.0).abs() < 0.01, "y1 at c should be 1000 nm, got {}", nm);
    }

    #[test]
    fn em_band_boundaries() {
        assert_eq!(em_band(0.005), "gamma");
        assert_eq!(em_band(5.0), "hard X-ray");
        assert_eq!(em_band(50.0), "soft X-ray");
        assert_eq!(em_band(300.0), "UV");
        assert_eq!(em_band(550.0), "visible");
        assert_eq!(em_band(1000.0), "IR");
        assert_eq!(em_band(5_000_000.0), "microwave");
        assert_eq!(em_band(2_000_000_000.0), "radio");
    }

    #[test]
    fn em_band_zero_and_negative() {
        assert_eq!(em_band(0.0), "—");
        assert_eq!(em_band(-1.0), "—");
    }

    #[test]
    fn em_band_detail_visible_colors() {
        assert_eq!(em_band_detail(420.0), "visible (violet)");
        assert_eq!(em_band_detail(470.0), "visible (blue)");
        assert_eq!(em_band_detail(530.0), "visible (green)");
        assert_eq!(em_band_detail(585.0), "visible (yellow)");
        assert_eq!(em_band_detail(600.0), "visible (orange)");
        assert_eq!(em_band_detail(700.0), "visible (red)");
    }

    #[test]
    fn em_band_detail_ir_subdivisions() {
        assert_eq!(em_band_detail(900.0), "near-IR");
        assert_eq!(em_band_detail(2000.0), "short-wave IR");
        assert_eq!(em_band_detail(5000.0), "mid-IR");
        assert_eq!(em_band_detail(10000.0), "thermal IR");
        assert_eq!(em_band_detail(100_000.0), "sub-mm (electron scale)");
    }

    #[test]
    fn format_wavelength_ranges() {
        assert_eq!(format_wavelength_si(0.5), "500.0 pm");
        assert_eq!(format_wavelength_si(500.0), "500.0 nm");
        assert_eq!(format_wavelength_si(2000.0), "2.00 μm");
        assert_eq!(format_wavelength_si(5_000_000.0), "5.00 mm");
        assert_eq!(format_wavelength_si(2_000_000_000.0), "2.000 m");
    }

    #[test]
    fn format_wavelength_zero() {
        assert_eq!(format_wavelength_si(0.0), "—");
        assert_eq!(format_wavelength_si(-1.0), "—");
    }

    #[test]
    fn rgb_green_at_530nm() {
        let (r, g, b) = wavelength_to_rgb(530.0);
        assert!(g > r && g > b, "530 nm should be green-dominant: r={} g={} b={}", r, g, b);
    }

    #[test]
    fn rgb_red_at_700nm() {
        let (r, g, b) = wavelength_to_rgb(650.0);
        assert!(r > g && r > b, "650 nm should be red-dominant: r={} g={} b={}", r, g, b);
    }

    #[test]
    fn rgb_blue_at_450nm() {
        let (r, g, b) = wavelength_to_rgb(450.0);
        assert!(b > g, "450 nm should be blue-dominant: r={} g={} b={}", r, g, b);
    }

    #[test]
    fn rgb_gray_outside_visible() {
        assert_eq!(wavelength_to_rgb(300.0), (0.5, 0.5, 0.5));
        assert_eq!(wavelength_to_rgb(800.0), (0.5, 0.5, 0.5));
    }

    #[test]
    fn rgb_values_in_range() {
        for nm in (380..=780).step_by(5) {
            let (r, g, b) = wavelength_to_rgb(nm as f64);
            assert!(r >= 0.0 && r <= 1.0, "r out of range at {} nm: {}", nm, r);
            assert!(g >= 0.0 && g <= 1.0, "g out of range at {} nm: {}", nm, g);
            assert!(b >= 0.0 && b <= 1.0, "b out of range at {} nm: {}", nm, b);
        }
    }
}
