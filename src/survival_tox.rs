// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #19 (`osabt survival`) — Acute ecotoxicology Probit / Logit `LC50` binomial MLE
//! and two-panel `.svg` Kaplan–Meier + Dose-Response plotter (`Finney 1971`, `OECD 203`).

use std::fmt::Write as _;

use crate::survival_core::{KaplanMeierCurve, LogRankTestResult};
use crate::wildlife_occupancy::{diag_hessian_se, nelder_mead_minimize};

/// Single concentration-mortality bioassay group.
#[derive(Debug, Clone)]
pub struct EcotoxBioassayPoint {
    pub concentration: f64,
    pub exposed_n: usize,
    pub dead_n: usize,
}

/// Link function for binomial dose-response `LCx` estimation.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum DoseResponseLink {
    Probit,
    Logit,
}

impl DoseResponseLink {
    pub fn label(self) -> &'static str {
        match self {
            Self::Probit => "Finney Probit",
            Self::Logit => "Binomial Logit",
        }
    }
}

/// Fitted `LC10`, `LC50`, and `LC90` acute toxicity summary with 95% confidence bounds.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EcotoxLc50Report {
    pub link: DoseResponseLink,
    pub intercept_alpha: f64,
    pub slope_beta: f64,
    pub lc10: f64,
    pub lc50: f64,
    pub lc50_ci95: (f64, f64),
    pub lc90: f64,
    pub log_likelihood: f64,
    pub pearson_chi2: f64,
}

/// Parses an acute bioassay CSV/TSV (`concentration,exposed_n,dead_n`).
pub fn parse_ecotox_csv(raw_csv: &str) -> Result<Vec<EcotoxBioassayPoint>, String> {
    let mut pts = Vec::new();
    let mut hdr = false;
    for (idx, line) in raw_csv.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();
        if !hdr {
            hdr = true;
            let c0 = cols[0].to_ascii_lowercase();
            if c0.contains("conc") || c0.contains("dose") {
                continue;
            }
        }
        if cols.len() < 3 {
            return Err(format!(
                "Ecotox CSV line {}: expected 3 columns (concentration,exposed_n,dead_n).",
                idx + 1
            ));
        }
        let conc = cols[0]
            .parse::<f64>()
            .map_err(|_| format!("Ecotox line {}: invalid concentration '{}'.", idx + 1, cols[0]))?;
        let exposed = cols[1]
            .parse::<usize>()
            .map_err(|_| format!("Ecotox line {}: invalid exposed_n '{}'.", idx + 1, cols[1]))?;
        let dead = cols[2]
            .parse::<usize>()
            .map_err(|_| format!("Ecotox line {}: invalid dead_n '{}'.", idx + 1, cols[2]))?;
        if dead > exposed || exposed == 0 {
            return Err(format!(
                "Ecotox line {}: dead_n ({dead}) must be <= exposed_n ({exposed}) > 0.",
                idx + 1
            ));
        }
        if conc > 0.0 {
            pts.push(EcotoxBioassayPoint {
                concentration: conc,
                exposed_n: exposed,
                dead_n: dead,
            });
        }
    }
    if pts.len() < 3 {
        return Err(String::from(
            "Ecotoxicology LC50 fitting requires at least 3 positive concentrations.",
        ));
    }
    pts.sort_by(|a, b| a.concentration.partial_cmp(&b.concentration).unwrap());
    Ok(pts)
}

/// Fits a Finney (1971) Probit or Binomial Logit concentration-mortality model (`P = F(alpha + beta * log10(C))`).
pub fn fit_ecotox_lc50(
    points: &[EcotoxBioassayPoint],
    link: DoseResponseLink,
) -> Result<EcotoxLc50Report, String> {
    if points.len() < 3 {
        return Err(String::from("At least 3 concentrations are required for LC50 fitting."));
    }
    let log_c: Vec<f64> = points.iter().map(|p| p.concentration.log10()).collect();
    let mean_lc = log_c.iter().sum::<f64>() / (log_c.len() as f64);

    let eval_prob = |z: f64| -> f64 {
        match link {
            DoseResponseLink::Probit => probit_norm_cdf(z),
            DoseResponseLink::Logit => 1.0 / (1.0 + (-z).exp()),
        }
    };

    let neg_ll = |th: &[f64]| -> f64 {
        if th[0].abs() > 25.0 || th[1].abs() > 25.0 {
            return 1e12 + th[0] * th[0] + th[1] * th[1];
        }
        let mut ll = 0.0_f64;
        for (i, pt) in points.iter().enumerate() {
            let p = eval_prob(th[0] + th[1] * log_c[i]).clamp(1e-12, 1.0 - 1e-12);
            let k = pt.dead_n as f64;
            let n = pt.exposed_n as f64;
            ll += k * p.ln() + (n - k) * (1.0 - p).ln();
        }
        -ll
    };

    let (best_th, best_nll) = nelder_mead_minimize(&neg_ll, &[-2.0 * mean_lc, 2.0], 1000, 1e-12);
    let alpha = best_th[0];
    let beta = best_th[1].max(1e-6);

    let se_ab = diag_hessian_se(&neg_ll, &[alpha, beta]);
    let x50 = -alpha / beta;
    let lc50 = 10_f64.powf(x50);
    let se_x50 = ((se_ab[0] / beta).powi(2) + (alpha * se_ab[1] / (beta * beta)).powi(2))
        .sqrt()
        .min(0.45);
    let lc50_ci95 = (
        10_f64.powf(x50 - 1.96 * se_x50),
        10_f64.powf(x50 + 1.96 * se_x50),
    );

    let (z10, z90) = match link {
        DoseResponseLink::Probit => (-1.281_551_6, 1.281_551_6),
        DoseResponseLink::Logit => (-(9.0_f64.ln()), 9.0_f64.ln()),
    };
    let lc10 = 10_f64.powf((z10 - alpha) / beta);
    let lc90 = 10_f64.powf((z90 - alpha) / beta);

    let mut chi2 = 0.0_f64;
    for (i, pt) in points.iter().enumerate() {
        let p = eval_prob(alpha + beta * log_c[i]).clamp(1e-9, 1.0 - 1e-9);
        let nf = pt.exposed_n as f64;
        let exp_d = nf * p;
        let diff = (pt.dead_n as f64) - exp_d;
        chi2 += (diff * diff) / (nf * p * (1.0 - p));
    }

    Ok(EcotoxLc50Report {
        link,
        intercept_alpha: alpha,
        slope_beta: beta,
        lc10,
        lc50,
        lc50_ci95,
        lc90,
        log_likelihood: -best_nll,
        pearson_chi2: chi2,
    })
}

/// Standard normal CDF `Phi(z)` via Abramowitz–Stegun 7.1.26 (`error < 1.5e-7`).
#[inline]
fn probit_norm_cdf(z: f64) -> f64 {
    let x = z / std::f64::consts::SQRT_2;
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let ax = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * ax);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-ax * ax).exp();
    0.5 * (1.0 + sign * y)
}

/// Renders a two-panel `.svg` figure pairing Kaplan–Meier step survival curves (`left`)
/// with the acute ecotoxicology Probit/Logit `LC50` dose-response curve (`right`).
pub fn render_survival_and_lc50_svg(
    curves: &[KaplanMeierCurve],
    log_rank: Option<&LogRankTestResult>,
    points: &[EcotoxBioassayPoint],
    lc50_rep: &EcotoxLc50Report,
    title: &str,
) -> String {
    let width = 940.0_f64;
    let height = 430.0_f64;
    let p1_x = 68.0_f64;
    let p1_w = 390.0_f64;
    let p2_x = 538.0_f64;
    let p2_w = 352.0_f64;
    let top = 58.0_f64;
    let ph = 305.0_f64;

    const COLORS: &[&str] = &["#2563eb", "#dc2626", "#16a34a", "#9333ea"];
    let mut svg = String::with_capacity(14_336);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#ffffff\"/>"
    );
    let lr_txt = log_rank
        .map(|l| format!(" (Log-Rank p = {:.4})", l.p_value))
        .unwrap_or_default();
    let _ = writeln!(
        svg,
        "  <text x=\"{p1_x:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Kaplan–Meier Survival{lr_txt} &amp; Ecotox LC50</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p1_x:.1}\" y=\"{top:.1}\" width=\"{p1_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p2_x:.1}\" y=\"{top:.1}\" width=\"{p2_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    let max_t = curves
        .iter()
        .filter_map(|c| c.steps.last().map(|s| s.time))
        .fold(1.0_f64, f64::max);
    let map_x1 = |t: f64| p1_x + (t / max_t).clamp(0.0, 1.0) * p1_w;
    let map_y = |prob: f64| top + ph - prob.clamp(0.0, 1.0) * ph;

    for (idx, c) in curves.iter().enumerate() {
        let col = COLORS[idx % COLORS.len()];
        let mut poly = String::new();
        let mut prev_s = 1.0_f64;
        for st in &c.steps {
            let xp = map_x1(st.time);
            let _ = write!(poly, "{xp:.1},{:.1} ", map_y(prev_s));
            let _ = write!(poly, "{xp:.1},{:.1} ", map_y(st.survival_s));
            if st.n_censored > 0 {
                let yc = map_y(st.survival_s);
                let _ = writeln!(
                    svg,
                    "  <line x1=\"{xp:.1}\" y1=\"{:.1}\" x2=\"{xp:.1}\" y2=\"{:.1}\" stroke=\"{col}\" stroke-width=\"1.6\"/>",
                    yc - 4.0,
                    yc + 4.0
                );
            }
            prev_s = st.survival_s;
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.4\"/>",
            poly.trim()
        );
    }

    let min_lc = points.first().map(|p| p.concentration.log10()).unwrap_or(0.0) - 0.25;
    let max_lc = points.last().map(|p| p.concentration.log10()).unwrap_or(2.0) + 0.25;
    let span_lc = (max_lc - min_lc).max(1e-6);
    let map_x2 = |lc: f64| p2_x + ((lc - min_lc) / span_lc).clamp(0.0, 1.0) * p2_w;

    let mut curve_pts = String::new();
    for k in 0..=48 {
        let lc = min_lc + (k as f64 / 48.0) * span_lc;
        let z = lc50_rep.intercept_alpha + lc50_rep.slope_beta * lc;
        let prob = match lc50_rep.link {
            DoseResponseLink::Probit => probit_norm_cdf(z),
            DoseResponseLink::Logit => 1.0 / (1.0 + (-z).exp()),
        };
        let _ = write!(curve_pts, "{:.1},{:.1} ", map_x2(lc), map_y(prob));
    }
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.4\"/>",
        curve_pts.trim()
    );
    for pt in points {
        let xp = map_x2(pt.concentration.log10());
        let yp = map_y((pt.dead_n as f64) / (pt.exposed_n as f64));
        let _ = writeln!(svg, "  <circle cx=\"{xp:.1}\" cy=\"{yp:.1}\" r=\"4.5\" fill=\"#dc2626\"/>");
    }

    let x50 = map_x2(lc50_rep.lc50.log10());
    let y50 = map_y(0.5);
    let _ = writeln!(
        svg,
        "  <line x1=\"{x50:.1}\" y1=\"{y50:.1}\" x2=\"{x50:.1}\" y2=\"{:.1}\" stroke=\"#0f172a\" stroke-dasharray=\"4,4\" stroke-width=\"1.4\"/>",
        top + ph
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Time (Days) — Survival S(t)</text>",
        p1_x + 0.5 * p1_w,
        height - 16.0
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">{} LC50 = {:.2} ug/L</text>",
        p2_x + 0.5 * p2_w,
        height - 16.0,
        lc50_rep.link.label(),
        lc50_rep.lc50
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

/// Built-in 6-concentration 96-h Copper Sulfate (`CuSO4`, `ug/L`) acute *Daphnia magna* bioassay (`n=30/conc`).
pub fn demo_ecotox_lc50_csv() -> &'static str {
    "\
# Acute 96-h CuSO4 Daphnia magna Bioassay (ug/L, n=30 per concentration)
concentration_ug_L,exposed_n,dead_n
5.0,30,1
10.0,30,4
20.0,30,11
40.0,30,20
80.0,30,27
160.0,30,29
"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_finney_probit_and_logit_lc50_cleanly() {
        let pts = parse_ecotox_csv(demo_ecotox_lc50_csv()).unwrap();
        let prob = fit_ecotox_lc50(&pts, DoseResponseLink::Probit).unwrap();
        let logit = fit_ecotox_lc50(&pts, DoseResponseLink::Logit).unwrap();
        assert!(prob.lc50 > 20.0 && prob.lc50 < 35.0);
        assert!(logit.lc50 > 20.0 && logit.lc50 < 35.0);
        assert!(prob.lc10 < prob.lc50 && prob.lc50 < prob.lc90);
    }
}
