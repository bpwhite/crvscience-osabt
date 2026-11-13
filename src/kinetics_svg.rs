// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #13 (`osabt kinetics`) — Multi-panel vector `.svg` visualizers for:
//! 1. **Global 2D Reversible Inhibition (`[S]` vs. `v0` + Diagnostic Eadie–Hofstee `v0` vs. `v0/[S]` Inset)**
//!    highlighting why nonlinear regression + visual Eadie–Hofstee / Lineweaver–Burk inspection work together.
//! 2. **Morrison (1969) Tight-Binding Fractional Velocity Titration** (`v_i / v_0` vs. `[I]_T`)
//!    with active-site stoichiometric tangent `[E]_T` callout.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::kinetics_core::{
    eval_inhibition_velocity, eval_morrison_velocity, GlobalInhibitionFit, KineticsPoint,
    MorrisonTightBindingFit,
};

/// Renders a two-panel `.svg` figure pairing the direct nonlinear Michaelis–Menten global fit (`left`)
/// with a diagnostic **Eadie–Hofstee (`v0` vs. `v0/[S]`)** panel (`right`) for mechanism verification.
pub fn render_global_inhibition_svg(
    points: &[KineticsPoint],
    best_fit: &GlobalInhibitionFit,
    title: &str,
) -> String {
    let width = 940.0_f64;
    let height = 450.0_f64;
    let p1_left = 68.0_f64;
    let p1_w = 390.0_f64;
    let p2_left = 535.0_f64;
    let p2_w = 355.0_f64;
    let top = 62.0_f64;
    let ph = 320.0_f64;

    let s_max = points
        .iter()
        .map(|p| p.substrate_um)
        .fold(100.0_f64, f64::max)
        * 1.06;
    let v_max = (best_fit.vmax * 1.08).max(10.0);

    let mut by_inh: BTreeMap<i64, Vec<KineticsPoint>> = BTreeMap::new();
    for &pt in points {
        let key = (pt.inhibitor_nm * 1000.0).round() as i64;
        by_inh.entry(key).or_default().push(pt);
    }

    let mut eh_x_max = 0.5_f64;
    for pt in points {
        if pt.substrate_um > 0.0 {
            eh_x_max = eh_x_max.max((pt.velocity / pt.substrate_um) * 1.15);
        }
    }

    let map_x1 = |s: f64| p1_left + (s / s_max).clamp(0.0, 1.0) * p1_w;
    let map_y = |v: f64| top + ph - (v / v_max).clamp(0.0, 1.0) * ph;
    let map_x2 = |vs: f64| p2_left + (vs / eh_x_max).clamp(0.0, 1.0) * p2_w;

    const COLORS: &[&str] = &["#2563eb", "#16a34a", "#d97706", "#dc2626", "#9333ea"];
    let mut svg = String::with_capacity(16_384);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#ffffff\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{p1_left:.0}\" y=\"26\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — {} (AICc w = {:.1}%)</text>",
        xml_escape(title),
        xml_escape(best_fit.mechanism.label()),
        best_fit.akaike_weight * 100.0
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{p1_left:.0}\" y=\"45\" font-family=\"monospace\" font-size=\"11.5\" fill=\"#475569\">Vmax = {:.2} | Km = {:.2} µM | Ki = {:.2} nM | R² = {:.4}</text>",
        best_fit.vmax, best_fit.km_um, best_fit.ki_nm, best_fit.r_squared
    );

    // Panel frames
    let _ = writeln!(
        svg,
        "  <rect x=\"{p1_left:.1}\" y=\"{top:.1}\" width=\"{p1_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p2_left:.1}\" y=\"{top:.1}\" width=\"{p2_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    let alpha = best_fit.alpha_factor.unwrap_or(1.0);
    for (idx, (&inh_key, grp_pts)) in by_inh.iter().enumerate() {
        let inh_nm = (inh_key as f64) / 1000.0;
        let col = COLORS[idx % COLORS.len()];

        // Direct Michaelis–Menten curve
        let mut poly1 = String::new();
        let mut poly2 = String::new();
        for step in 1..=90 {
            let s_val = (step as f64 / 90.0) * s_max;
            let v_pred = eval_inhibition_velocity(
                best_fit.mechanism,
                s_val,
                inh_nm,
                best_fit.vmax,
                best_fit.km_um,
                best_fit.ki_nm,
                alpha,
            );
            let _ = write!(poly1, "{:.2},{:.2} ", map_x1(s_val), map_y(v_pred));
            let _ = write!(poly2, "{:.2},{:.2} ", map_x2(v_pred / s_val), map_y(v_pred));
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.3\"/>",
            poly1.trim()
        );
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"1.8\" stroke-dasharray=\"5,3\"/>",
            poly2.trim()
        );

        for pt in grp_pts {
            let _ = writeln!(
                svg,
                "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"4.0\" fill=\"{col}\"/>",
                map_x1(pt.substrate_um),
                map_y(pt.velocity)
            );
            if pt.substrate_um > 0.0 {
                let _ = writeln!(
                    svg,
                    "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"3.8\" fill=\"{col}\"/>",
                    map_x2(pt.velocity / pt.substrate_um),
                    map_y(pt.velocity)
                );
            }
        }
    }

    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Substrate [S] (µM) — Direct Nonlinear Global Fit</text>",
        p1_left + 0.5 * p1_w,
        height - 18.0
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Eadie–Hofstee Diagnostic: v0 / [S] vs. v0</text>",
        p2_left + 0.5 * p2_w,
        height - 18.0
    );

    svg.push_str("</svg>\n");
    svg
}

/// Renders an `.svg` Morrison (1969) tight-binding titration curve (`[I]_T` vs. `v_i`) highlighting
/// the stoichiometric active-site inflection `[E]_T` and intrinsic `Ki`.
pub fn render_morrison_svg(
    observations: &[(f64, f64)],
    fit: &MorrisonTightBindingFit,
    title: &str,
) -> String {
    let width = 820.0_f64;
    let height = 430.0_f64;
    let left = 74.0_f64;
    let top = 58.0_f64;
    let pw = 520.0_f64;
    let ph = 310.0_f64;

    let i_max = observations
        .iter()
        .map(|p| p.0)
        .fold(fit.active_enzyme_et_nm * 2.0, f64::max)
        * 1.06;
    let v_max = (fit.v0_unhibited * 1.08).max(1.0);

    let map_x = |i_nm: f64| left + (i_nm / i_max).clamp(0.0, 1.0) * pw;
    let map_y = |v: f64| top + ph - (v / v_max).clamp(0.0, 1.0) * ph;

    let mut svg = String::with_capacity(8192);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#ffffff\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Morrison (1969) Tight-Binding Quadratic Fit</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    // Stoichiometric active-site titration tangent `(0, v0) -> ([E]_T, 0)`
    let x_et = map_x(fit.active_enzyme_et_nm);
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{:.1}\" x2=\"{x_et:.1}\" y2=\"{:.1}\" stroke=\"#94a3b8\" stroke-width=\"1.5\" stroke-dasharray=\"6,4\"/>",
        map_y(fit.v0_unhibited),
        map_y(0.0)
    );

    let mut poly = String::new();
    for s in 0..=120 {
        let it = (s as f64 / 120.0) * i_max;
        let vp = eval_morrison_velocity(it, fit.v0_unhibited, fit.active_enzyme_et_nm, fit.ki_app_nm);
        let _ = write!(poly, "{:.2},{:.2} ", map_x(it), map_y(vp));
    }
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#dc2626\" stroke-width=\"2.6\"/>",
        poly.trim()
    );

    for &(it, vy) in observations {
        let _ = writeln!(
            svg,
            "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"4.2\" fill=\"#1e293b\"/>",
            map_x(it),
            map_y(vy)
        );
    }

    let lx = left + pw + 18.0;
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11.5\" font-weight=\"700\" fill=\"#0f172a\">[E]_T  = {:.2} nM</text>",
        top + 28.0,
        fit.active_enzyme_et_nm
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11.5\" font-weight=\"700\" fill=\"#dc2626\">Ki_app = {:.3} nM</text>",
        top + 50.0,
        fit.ki_app_nm
    );
    if let Some(ki_int) = fit.intrinsic_ki_nm {
        let _ = writeln!(
            svg,
            "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11.5\" fill=\"#2563eb\">Ki     = {:.3} nM</text>",
            top + 72.0,
            ki_int
        );
    }
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11.5\" fill=\"#475569\">R²     = {:.4}</text>",
        top + 94.0,
        fit.r_squared
    );

    svg.push_str("</svg>\n");
    svg
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}
