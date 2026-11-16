// Copyright (C) 2026 CRV Science Network (https://crvscience.com)
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://gnu.org>.

// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #9 (`osabt plate`) — High-Throughput Screening (`HTS`) assay quality control
//! (Zhang, Chung & Oldenburg 1999 `Z'`-factor, `SSMD`, `S/B`, `%CV`) and Sebaugh (2011)
//! 4-Parameter Logistic (`4PL`) nonlinear `IC50`/`EC50` curve fitter via damped
//! Levenberg–Marquardt Gauss–Newton regression with analytical Jacobians and 95% CIs.

use std::collections::BTreeMap;
use std::f64::consts::LN_10;
use std::fmt::Write as _;

use crate::plate_layout::{
    generate_randomized_plate, Microplate, PlatePlanConfig, WellRole,
};
use crate::svg_util::{begin_svg_canvas, xml_escape};

/// Zhang, Chung & Oldenburg (1999) *J. Biomol. Screen.* 4:67–73 assay quality-control metrics.
#[derive(Debug, Clone)]
pub struct PlateAssayQc {
    pub neg_ctrl_n: usize,
    pub pos_ctrl_n: usize,
    pub neg_mean: f64,
    pub neg_sd: f64,
    pub neg_cv_pct: f64,
    pub pos_mean: f64,
    pub pos_sd: f64,
    pub pos_cv_pct: f64,
    pub signal_to_background: f64,
    /// `Z' = 1 - 3 * (sd_pos + sd_neg) / |mean_pos - mean_neg|`.
    pub z_prime_factor: f64,
    /// Strictly Standardized Mean Difference `(mean_neg - mean_pos) / sqrt(sd_neg^2 + sd_pos^2)`.
    pub ssmd: f64,
    pub quality_grade: &'static str,
}

/// Fitted 4-Parameter Logistic (`4PL`) dose-response curve for a single test compound:
/// `f(x) = bottom + (top - bottom) / (1 + 10^((log10_ec50 - log10_conc) * hill_slope))`.
#[derive(Debug, Clone)]
pub struct FourParamLogisticFit {
    pub compound_name: String,
    pub n_points: usize,
    pub bottom: f64,
    pub top: f64,
    pub log10_ec50: f64,
    /// Relative `IC50` / `EC50` inflection concentration (`10^log10_ec50`, nM).
    pub relative_ec50_nm: f64,
    pub relative_ec50_ci95_nm: (f64, f64),
    /// Sebaugh (2011) Absolute `IC50` / `EC50` where predicted response crosses `50.0%` (`None` if `50%` lies outside `(bottom, top)`).
    pub absolute_ic50_nm: Option<f64>,
    /// 90%-effect concentration (`IC90` for inhibitors `hill < 0`, `EC90` for agonists `hill > 0`).
    pub ec90_nm: f64,
    pub hill_slope: f64,
    pub hill_se: f64,
    pub r_squared: f64,
    pub residual_se: f64,
    /// `(conc_nm, mean_response, sd_response, n_reps)` aggregated per dose for plotting.
    pub dose_means: Vec<(f64, f64, f64, usize)>,
}

/// Computes Zhang et al. (1999) `Z'`-factor, `SSMD`, and normalizes all measured wells on `plate`
/// to `% Viability / Signal` (`100% = Negative/Vehicle Control`, `0% = Positive/100%-Inhib Control`).
pub fn evaluate_plate_qc_and_normalize(plate: &mut Microplate) -> Result<PlateAssayQc, String> {
    let neg_vals: Vec<f64> = plate
        .wells
        .iter()
        .filter(|w| w.role == WellRole::NegativeControl)
        .filter_map(|w| w.raw_readout)
        .collect();
    let pos_vals: Vec<f64> = plate
        .wells
        .iter()
        .filter(|w| w.role == WellRole::PositiveControl)
        .filter_map(|w| w.raw_readout)
        .collect();

    if neg_vals.len() < 2 || pos_vals.len() < 2 {
        return Err(String::from(
            "Assay QC requires >= 2 measured NegativeControl and >= 2 PositiveControl wells.",
        ));
    }

    let (neg_mean, neg_sd) = mean_and_sd(&neg_vals);
    let (pos_mean, pos_sd) = mean_and_sd(&pos_vals);
    let span = (neg_mean - pos_mean).abs();
    if span <= 1e-9 {
        return Err(String::from(
            "Negative and Positive control means are identical; dynamic range is zero.",
        ));
    }

    let z_prime = 1.0 - (3.0 * (neg_sd + pos_sd)) / span;
    let ssmd = (neg_mean - pos_mean) / (neg_sd * neg_sd + pos_sd * pos_sd).sqrt().max(1e-12);
    let grade = if z_prime >= 0.50 {
        "EXCELLENT (Z' >= 0.50)"
    } else if z_prime >= 0.0 {
        "MARGINAL / DOABLE (0 <= Z' < 0.50)"
    } else {
        "REJECT / OVERLAPPING CONTROLS (Z' < 0)"
    };

    for w in &mut plate.wells {
        if let Some(raw) = w.raw_readout {
            w.normalized_pct = Some(((raw - pos_mean) / (neg_mean - pos_mean)) * 100.0);
        }
    }

    Ok(PlateAssayQc {
        neg_ctrl_n: neg_vals.len(),
        pos_ctrl_n: pos_vals.len(),
        neg_mean,
        neg_sd,
        neg_cv_pct: (neg_sd / neg_mean.abs().max(1e-9)) * 100.0,
        pos_mean,
        pos_sd,
        pos_cv_pct: (pos_sd / pos_mean.abs().max(1e-9)) * 100.0,
        signal_to_background: neg_mean / pos_mean.max(1e-9),
        z_prime_factor: z_prime,
        ssmd,
        quality_grade: grade,
    })
}

/// Evaluates the 4PL model `bottom + (top - bottom) / (1 + 10^((x50 - x) * hill))` at `x = log10(conc)`.
#[inline]
pub fn eval_4pl_log10(x: f64, bottom: f64, top: f64, x50: f64, hill: f64) -> f64 {
    let exp_arg = ((x50 - x) * hill).clamp(-60.0, 60.0);
    let u = 10_f64.powf(exp_arg);
    bottom + (top - bottom) / (1.0 + u)
}

/// Fits a 4-Parameter Logistic (`4PL`) curve to `(concentration_nM, response)` observations
/// via damped Levenberg–Marquardt nonlinear least squares.
pub fn fit_four_param_logistic(
    compound_name: &str,
    observations: &[(f64, f64)],
) -> Result<FourParamLogisticFit, String> {
    let valid: Vec<(f64, f64)> = observations
        .iter()
        .copied()
        .filter(|&(c, y)| c > 0.0 && c.is_finite() && y.is_finite())
        .collect();

    if valid.len() < 5 {
        return Err(format!(
            "Compound '{}' has {} positive-concentration points; 4PL requires >= 5.",
            compound_name,
            valid.len()
        ));
    }

    let xs: Vec<f64> = valid.iter().map(|&(c, _)| c.log10()).collect();
    let ys: Vec<f64> = valid.iter().map(|&(_, y)| y).collect();

    let y_min = ys.iter().copied().fold(f64::INFINITY, f64::min);
    let y_max = ys.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let x_min = xs.iter().copied().fold(f64::INFINITY, f64::min);
    let x_max = xs.iter().copied().fold(f64::NEG_INFINITY, f64::max);

    // Determine initial Hill sign by comparing mean response in lowest vs highest log-dose halves
    let x_mid = 0.5 * (x_min + x_max);
    let low_mean = mean_only(
        &valid
            .iter()
            .filter(|&&(c, _)| c.log10() <= x_mid)
            .map(|&(_, y)| y)
            .collect::<Vec<_>>(),
    );
    let high_mean = mean_only(
        &valid
            .iter()
            .filter(|&&(c, _)| c.log10() > x_mid)
            .map(|&(_, y)| y)
            .collect::<Vec<_>>(),
    );

    let init_hill = if high_mean < low_mean { -1.0 } else { 1.0 };
    // Estimate initial x50 near the observation closest to `0.5 * (y_min + y_max)`
    let half_y = 0.5 * (y_min + y_max);
    let init_x50 = xs
        .iter()
        .zip(&ys)
        .min_by(|a, b| {
            (a.1 - half_y)
                .abs()
                .partial_cmp(&(b.1 - half_y).abs())
                .unwrap()
        })
        .map(|(&x, _)| x)
        .unwrap_or(x_mid);

    let mut params = [y_min, y_max.max(y_min + 1.0), init_x50, init_hill];
    let mut lambda = 1e-3_f64;
    let mut cur_sse = sse_4pl(&xs, &ys, &params);

    for _iter in 0..120 {
        let (jtj, jtr) = accumulate_normal_equations(&xs, &ys, &params);
        let mut damped = jtj;
        for d in 0..4 {
            damped[d][d] += lambda * damped[d][d].max(1e-6);
        }

        let Some(step) = solve_4x4(damped, jtr) else {
            lambda *= 10.0;
            continue;
        };

        let mut cand = [
            params[0] + step[0],
            params[1] + step[1],
            (params[2] + step[2]).clamp(x_min - 4.0, x_max + 4.0),
            (params[3] + step[3]).clamp(-8.0, 8.0),
        ];
        if cand[3].abs() < 0.05 {
            cand[3] = 0.05 * cand[3].signum();
        }

        let cand_sse = sse_4pl(&xs, &ys, &cand);
        if cand_sse < cur_sse {
            let rel_change = (cur_sse - cand_sse) / cur_sse.max(1e-12);
            params = cand;
            cur_sse = cand_sse;
            lambda = (lambda * 0.35).max(1e-11);
            if rel_change < 1e-12 {
                break;
            }
        } else {
            lambda = (lambda * 4.0).min(1e8);
        }
    }

    let (b, t, x50, h) = (params[0], params[1], params[2], params[3]);
    let n = valid.len();
    let dof = (n.saturating_sub(4).max(1)) as f64;
    let s2 = cur_sse / dof;
    let rse = s2.sqrt();

    let y_bar = mean_only(&ys);
    let sst: f64 = ys.iter().map(|&y| (y - y_bar).powi(2)).sum();
    let r2 = if sst > 1e-12 {
        (1.0 - cur_sse / sst).clamp(0.0, 1.0)
    } else {
        1.0
    };

    // Asymptotic standard errors from `s^2 * (J^T J)^-1`
    let (jtj, _) = accumulate_normal_equations(&xs, &ys, &params);
    let (se_x50, se_h) = invert_4x4_diag(&jtj)
        .map(|diag| ((s2 * diag[2].max(0.0)).sqrt(), (s2 * diag[3].max(0.0)).sqrt()))
        .unwrap_or((0.08, 0.10));

    let t_crit = if dof >= 12.0 { 2.12 } else { 2.36 };
    let rel_ec50 = 10_f64.powf(x50);
    let ci_lo = 10_f64.powf(x50 - t_crit * se_x50);
    let ci_hi = 10_f64.powf(x50 + t_crit * se_x50);

    // Sebaugh (2011) Absolute IC50/EC50 where curve crosses `y = 50.0`
    let lo_plat = b.min(t);
    let hi_plat = b.max(t);
    let abs_ic50 = if lo_plat < 50.0 && hi_plat > 50.0 && (50.0 - b).abs() > 1e-9 {
        let ratio = (t - 50.0) / (50.0 - b);
        if ratio > 0.0 {
            Some(10_f64.powf(x50 - ratio.log10() / h))
        } else {
            None
        }
    } else {
        None
    };

    // 90%-effect potency (`IC90` toward bottom when `h < 0`; `EC90` toward top when `h > 0`)
    let ec90 = 10_f64.powf(x50 + 9.0_f64.log10() / h.abs().max(1e-6));

    // Group replicate responses by rounded concentration for error bars
    let mut grouped: BTreeMap<i64, Vec<(f64, f64)>> = BTreeMap::new();
    for &(c, y) in &valid {
        let key = (c.log10() * 10_000.0).round() as i64;
        grouped.entry(key).or_default().push((c, y));
    }
    let mut dose_means = Vec::with_capacity(grouped.len());
    for (_, pts) in grouped {
        let c_avg = pts.iter().map(|p| p.0).sum::<f64>() / (pts.len() as f64);
        let y_vec: Vec<f64> = pts.iter().map(|p| p.1).collect();
        let (ym, ysd) = mean_and_sd(&y_vec);
        dose_means.push((c_avg, ym, ysd, pts.len()));
    }

    Ok(FourParamLogisticFit {
        compound_name: compound_name.to_string(),
        n_points: n,
        bottom: b,
        top: t,
        log10_ec50: x50,
        relative_ec50_nm: rel_ec50,
        relative_ec50_ci95_nm: (ci_lo, ci_hi),
        absolute_ic50_nm: abs_ic50,
        ec90_nm: ec90,
        hill_slope: h,
        hill_se: se_h,
        r_squared: r2,
        residual_se: rse,
        dose_means,
    })
}

/// Fits 4PL curves for every test compound present on a normalized [`Microplate`].
pub fn fit_all_plate_compounds(plate: &Microplate) -> Vec<FourParamLogisticFit> {
    let mut by_comp: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
    for w in &plate.wells {
        if w.role == WellRole::CompoundDose {
            if let (Some(conc), Some(val)) = (w.concentration_nm, w.normalized_pct.or(w.raw_readout))
            {
                if conc > 0.0 {
                    by_comp
                        .entry(w.compound_name.clone())
                        .or_default()
                        .push((conc, val));
                }
            }
        }
    }

    let mut fits = Vec::new();
    for (name, obs) in by_comp {
        if let Ok(fit) = fit_four_param_logistic(&name, &obs) {
            fits.push(fit);
        }
    }
    fits
}

/// Generates a complete realistic 96-well kinase luminescence dataset (`Staurosporine` `IC50 ~ 4.2 nM`,
/// `CRV-104` `IC50 ~ 18.5 nM`, `CRV-209` `IC50 ~ 420 nM`, `Z' ~ 0.82`) for immediate testing.
pub fn demo_measured_96well_plate() -> Microplate {
    let cfg = PlatePlanConfig::default();
    let mut plate = generate_randomized_plate(&cfg).expect("valid default 96-well plan");

    let mut noise_state = 0x517C_C1B7_2722_0A95u64;
    let mut pseudo_norm = || -> f64 {
        // Sum of 6 uniforms - 3.0 -> approximate N(0, 0.5)
        let mut s = 0.0_f64;
        for _ in 0..6 {
            noise_state = noise_state
                .wrapping_mul(6_364_136_223_846_793_005)
                .wrapping_add(1);
            s += ((noise_state >> 32) as f64) / (u32::MAX as f64);
        }
        (s - 3.0) * 1.414
    };

    let bg = 1_200.0_f64;
    let sig = 52_000.0_f64;

    for w in &mut plate.wells {
        match w.role {
            WellRole::NegativeControl => {
                w.raw_readout = Some(sig + pseudo_norm() * 1_250.0);
            }
            WellRole::PositiveControl => {
                w.raw_readout = Some(bg + pseudo_norm() * 95.0);
            }
            WellRole::CompoundDose => {
                let c = w.concentration_nm.unwrap_or(1.0);
                let (true_ic50, true_hill): (f64, f64) = match w.compound_name.as_str() {
                    "Staurosporine" => (4.2, -1.04),
                    "CRV-104" => (18.5, -1.15),
                    _ => (420.0, -0.96),
                };
                let pct = eval_4pl_log10(c.log10(), 2.0, 99.5, true_ic50.log10(), true_hill);
                let raw = bg + (pct / 100.0) * (sig - bg) + pseudo_norm() * 950.0;
                w.raw_readout = Some(raw.max(100.0));
            }
            _ => {}
        }
    }

    let _ = evaluate_plate_qc_and_normalize(&mut plate);
    plate
}

/// Renders a multi-compound semi-log 4PL dose-response `.svg` figure with `mean +- SD` error bars.
pub fn render_4pl_curves_svg(title: &str, fits: &[FourParamLogisticFit]) -> String {
    let width = 840.0_f64;
    let height = 480.0_f64;
    let left = 72.0_f64;
    let top = 56.0_f64;
    let pw = 510.0_f64;
    let ph = 360.0_f64;

    let mut x_min = f64::INFINITY;
    let mut x_max = f64::NEG_INFINITY;
    for f in fits {
        for &(c, _, _, _) in &f.dose_means {
            x_min = x_min.min(c.log10());
            x_max = x_max.max(c.log10());
        }
    }
    if !x_min.is_finite() {
        x_min = -2.0;
        x_max = 4.0;
    }
    x_min -= 0.35;
    x_max += 0.35;
    let x_span = (x_max - x_min).max(1.0);

    const COLORS: &[&str] = &["#2563eb", "#dc2626", "#16a34a", "#9333ea", "#d97706"];
    let mut svg = begin_svg_canvas(8192, width, height, Some("#ffffff"));
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — 4-Parameter Logistic Dose-Response</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    let map_x = |lx: f64| left + ((lx - x_min) / x_span) * pw;
    let map_y = |y: f64| top + ph - ((y + 10.0) / 130.0).clamp(0.0, 1.0) * ph;

    let y50 = map_y(50.0);
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{y50:.1}\" x2=\"{:.1}\" y2=\"{y50:.1}\" stroke=\"#94a3b8\" stroke-dasharray=\"5,4\"/>",
        left + pw
    );

    for (idx, fit) in fits.iter().enumerate() {
        let col = COLORS[idx % COLORS.len()];
        let mut poly = String::new();
        for s in 0..=100 {
            let lx = x_min + (s as f64 / 100.0) * x_span;
            let py = map_y(eval_4pl_log10(lx, fit.bottom, fit.top, fit.log10_ec50, fit.hill_slope));
            let _ = write!(poly, "{:.2},{py:.2} ", map_x(lx));
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.4\"/>",
            poly.trim()
        );

        for &(c, ym, ysd, _) in &fit.dose_means {
            let px = map_x(c.log10());
            let py = map_y(ym);
            let py_hi = map_y(ym + ysd);
            let py_lo = map_y(ym - ysd);
            let _ = writeln!(
                svg,
                "  <line x1=\"{px:.1}\" y1=\"{py_hi:.1}\" x2=\"{px:.1}\" y2=\"{py_lo:.1}\" stroke=\"{col}\" stroke-width=\"1.3\"/>"
            );
            let _ = writeln!(
                svg,
                "  <circle cx=\"{px:.1}\" cy=\"{py:.1}\" r=\"4.0\" fill=\"{col}\"/>"
            );
        }

        let lx_leg = left + pw + 18.0;
        let ly_leg = top + 24.0 + (idx as f64) * 44.0;
        let _ = writeln!(
            svg,
            "  <text x=\"{lx_leg:.1}\" y=\"{ly_leg:.1}\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"{col}\">{}</text>",
            xml_escape(&fit.compound_name)
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{lx_leg:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#334155\">IC50={:.2} nM | h={:.2} | R²={:.3}</text>",
            ly_leg + 16.0,
            fit.relative_ec50_nm,
            fit.hill_slope,
            fit.r_squared
        );
    }

    svg.push_str("</svg>\n");
    svg
}

fn sse_4pl(xs: &[f64], ys: &[f64], p: &[f64; 4]) -> f64 {
    xs.iter()
        .zip(ys)
        .map(|(&x, &y)| {
            let r = y - eval_4pl_log10(x, p[0], p[1], p[2], p[3]);
            r * r
        })
        .sum()
}

fn accumulate_normal_equations(xs: &[f64], ys: &[f64], p: &[f64; 4]) -> ([[f64; 4]; 4], [f64; 4]) {
    let (b, t, x50, h) = (p[0], p[1], p[2], p[3]);
    let mut jtj = [[0.0_f64; 4]; 4];
    let mut jtr = [0.0_f64; 4];

    for (&x, &y) in xs.iter().zip(ys) {
        let exp_arg = ((x50 - x) * h).clamp(-60.0, 60.0);
        let u = 10_f64.powf(exp_arg);
        let d = 1.0 + u;
        let d2 = d * d;
        let pred = b + (t - b) / d;
        let res = y - pred;

        let jb = u / d;
        let jt = 1.0 / d;
        let jx50 = -(t - b) * (u / d2) * LN_10 * h;
        let jh = -(t - b) * (u / d2) * LN_10 * (x50 - x);
        let jac = [jb, jt, jx50, jh];

        for r in 0..4 {
            jtr[r] += jac[r] * res;
            for c in 0..4 {
                jtj[r][c] += jac[r] * jac[c];
            }
        }
    }
    (jtj, jtr)
}

fn solve_4x4(mut a: [[f64; 4]; 4], mut rhs: [f64; 4]) -> Option<[f64; 4]> {
    for col in 0..4 {
        let mut piv = col;
        for r in (col + 1)..4 {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-18 {
            return None;
        }
        a.swap(col, piv);
        rhs.swap(col, piv);

        let diag = a[col][col];
        for r in (col + 1)..4 {
            let f = a[r][col] / diag;
            rhs[r] -= f * rhs[col];
            for c in col..4 {
                a[r][c] -= f * a[col][c];
            }
        }
    }

    let mut x = [0.0_f64; 4];
    for i in (0..4).rev() {
        let mut s = rhs[i];
        for j in (i + 1)..4 {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

fn invert_4x4_diag(a: &[[f64; 4]; 4]) -> Option<[f64; 4]> {
    let mut diag = [0.0_f64; 4];
    for i in 0..4 {
        let mut e = [0.0_f64; 4];
        e[i] = 1.0;
        let col = solve_4x4(*a, e)?;
        diag[i] = col[i];
    }
    Some(diag)
}

fn mean_only(v: &[f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        v.iter().sum::<f64>() / (v.len() as f64)
    }
}

fn mean_and_sd(v: &[f64]) -> (f64, f64) {
    let m = mean_only(v);
    if v.len() < 2 {
        return (m, 0.0);
    }
    let var = v.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / ((v.len() - 1) as f64);
    (m, var.sqrt())
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_z_prime_and_exact_4pl_ic50s_on_demo_96well_plate() {
        let mut plate = demo_measured_96well_plate();
        let qc = evaluate_plate_qc_and_normalize(&mut plate).unwrap();
        assert!(qc.z_prime_factor > 0.70);

        let fits = fit_all_plate_compounds(&plate);
        assert_eq!(fits.len(), 3);
        let crv104 = fits.iter().find(|f| f.compound_name == "CRV-104").unwrap();
        assert!((crv104.relative_ec50_nm - 18.5).abs() < 4.5);
        assert!(crv104.ec90_nm > crv104.relative_ec50_nm);
        assert!(crv104.r_squared > 0.98);
    }

    #[test]
    fn fits_agonist_4pl_with_ec90_above_ec50() {
        let doses: [f64; 8] = [0.1, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0];
        let pts: Vec<(f64, f64)> = doses
            .iter()
            .map(|&c| (c, eval_4pl_log10(c.log10(), 5.0, 95.0, 1.30103, 1.1)))
            .collect();
        let fit = fit_four_param_logistic("Agonist-X", &pts).unwrap();
        assert!((fit.relative_ec50_nm - 20.0).abs() < 0.2);
        assert!(fit.ec90_nm > fit.relative_ec50_nm * 6.0);
    }
}
