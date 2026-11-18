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

//! Tool #15 (`osabt wildlife`) — Closed & open-population capture-recapture engine:
//! - **Chapman (1951) / Seber (1970)** unbiased 2-sample Lincoln–Petersen abundance & log-normal 95% CI
//! - **Otis et al. (1978)** closed-population `M0` (null equal catchability) & Darroch `Mt` (time-varying)
//!   combinatorial profile-likelihood abundance estimators + **Schnabel (1938)** multi-census estimator
//! - **Chao (1987)** *Biometrics* 43(4):783–791 (`https://pubmed.ncbi.nlm.nih.gov/3427163/`) `Mh` unequal-catchability estimator
//! - **Cormack–Jolly–Seber (`CJS`)** — Lebreton et al. (1992) *Ecol. Monogr.* 62(1):67–118
//!   (`https://digitalcommons.usf.edu/kip_articles/7202/`) `m`-array open-population apparent survival (`phi`) & recapture (`p`).

use std::fmt::Write as _;

use crate::svg_util::xml_escape;
use crate::wildlife_occupancy::{
    diag_hessian_se, nelder_mead_minimize, OccupancySuiteReport,
};

/// Unbiased 2-sample Chapman (1951) / Seber (1970) Lincoln–Petersen mark-recapture abundance estimate.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChapmanPetersenEstimate {
    pub marked_sample_1: usize,
    pub captured_sample_2: usize,
    pub recaptured_r: usize,
    pub distinct_observed_mt1: usize,
    pub chapman_n_hat: f64,
    pub seber_se: f64,
    /// Asymmetric log-normal 95% CI bounded strictly above `distinct_observed_mt1`.
    pub ci95: (f64, f64),
}

/// Closed-population abundance estimate for one capture-recapture model (`M0`, `Mt`, `Chao Mh`, or `Schnabel`).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ClosedPopulationEstimate {
    pub model_name: String,
    pub n_hat: f64,
    pub se: f64,
    pub ci95: (f64, f64),
    pub mean_capture_p: f64,
}

/// Complete multi-occasion closed-population capture-recapture summary (`M0`, `Mt`, `Chao Mh`, `Schnabel`).
#[derive(Debug, Clone)]
pub struct ClosedCaptureSuiteReport {
    pub distinct_individuals_mt1: usize,
    pub occasions_k: usize,
    pub total_captures: usize,
    pub captures_per_occasion: Vec<usize>,
    /// Capture-frequency spectrum `f_1..f_K` (`freq_counts[r - 1]` = animals caught `r` times).
    pub frequency_spectrum_f: Vec<usize>,
    pub m0_equal_catchability: ClosedPopulationEstimate,
    pub mt_darroch_time_varying: ClosedPopulationEstimate,
    pub mh_chao1987_heterogeneity: ClosedPopulationEstimate,
    pub schnabel_multicensus: ClosedPopulationEstimate,
}

/// Open-population Cormack–Jolly–Seber (`CJS`, Lebreton et al. 1992) apparent survival & recapture report.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CjsSurvivalReport {
    pub occasions_k: usize,
    pub total_releases: usize,
    pub total_recaptures: usize,
    /// Apparent survival probability `phi` between consecutive occasions (`1 - mortality - permanent emigration`).
    pub phi_hat: f64,
    pub phi_se: f64,
    pub phi_ci95: (f64, f64),
    /// Conditional recapture probability `p` given alive & present in study area.
    pub p_hat: f64,
    pub p_se: f64,
    pub p_ci95: (f64, f64),
    pub log_likelihood: f64,
    pub aicc: f64,
    /// Sufficient-statistic Burnham `m`-array: `(released_R_i, recaptured_m_ij_vec, never_seen_again)`.
    pub m_array: Vec<(usize, Vec<usize>, usize)>,
}

/// Computes Chapman's (1951) unbiased Lincoln–Petersen estimator `N_hat = (M+1)(C+1)/(R+1) - 1`,
/// Seber's (1970) variance, and Chao/Burnham log-normal 95% confidence bounds (`>= M + C - R`).
pub fn calc_chapman_lincoln_petersen(
    marked_m: usize,
    captured_c: usize,
    recaptured_r: usize,
) -> Result<ChapmanPetersenEstimate, String> {
    if marked_m == 0 || captured_c == 0 {
        return Err(String::from(
            "Marked (M) and second-sample captured (C) counts must both be >= 1.",
        ));
    }
    if recaptured_r > marked_m.min(captured_c) {
        return Err(format!(
            "Recaptures R ({recaptured_r}) cannot exceed min(M={marked_m}, C={captured_c})."
        ));
    }

    let m = marked_m as f64;
    let c = captured_c as f64;
    let r = recaptured_r as f64;
    let distinct = marked_m + captured_c - recaptured_r;

    let n_hat = ((m + 1.0) * (c + 1.0)) / (r + 1.0) - 1.0;
    let var = ((m + 1.0) * (c + 1.0) * (m - r) * (c - r))
        / ((r + 1.0).powi(2) * (r + 2.0));
    let se = var.max(0.0).sqrt();
    let ci95 = lognormal_abundance_ci(distinct as f64, n_hat, var);

    Ok(ChapmanPetersenEstimate {
        marked_sample_1: marked_m,
        captured_sample_2: captured_c,
        recaptured_r,
        distinct_observed_mt1: distinct,
        chapman_n_hat: n_hat,
        seber_se: se,
        ci95,
    })
}

/// Parses an individual `0/1` capture-history CSV/TSV (`animal_id,t1,t2,...,tK` or compact string `10101`).
pub fn parse_capture_history_csv(raw_csv: &str) -> Result<(Vec<String>, Vec<Vec<bool>>), String> {
    let mut ids = Vec::new();
    let mut histories: Vec<Vec<bool>> = Vec::new();
    let mut header_checked = false;

    for (line_idx, line) in raw_csv.trim_start_matches('\u{feff}').lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();

        if !header_checked {
            header_checked = true;
            let c0 = cols[0].to_ascii_lowercase();
            if c0.contains("id") || c0.contains("animal") || c0.contains("tag") {
                continue;
            }
        }

        if cols.len() < 2 {
            return Err(format!(
                "Capture-history CSV line {}: expected animal_id plus capture history.",
                line_idx + 1
            ));
        }

        let hist: Vec<bool> = if cols.len() == 2 && cols[1].len() >= 2 {
            // Compact Mark/RMARK capture-history string such as `10110`
            let mut v = Vec::with_capacity(cols[1].len());
            for ch in cols[1].chars() {
                match ch {
                    '1' => v.push(true),
                    '0' => v.push(false),
                    _ => {
                        return Err(format!(
                            "Line {}: invalid character '{}' in capture history string.",
                            line_idx + 1,
                            ch
                        ));
                    }
                }
            }
            v
        } else {
            let mut v = Vec::with_capacity(cols.len() - 1);
            for &tok in &cols[1..] {
                match tok {
                    "1" => v.push(true),
                    "0" => v.push(false),
                    other => {
                        return Err(format!(
                            "Line {}: invalid capture value '{}' (expected 0 or 1).",
                            line_idx + 1,
                            other
                        ));
                    }
                }
            }
            v
        };

        if let Some(first) = histories.first() {
            if hist.len() != first.len() {
                return Err(format!(
                    "Line {}: capture history length ({}) does not match occasion count ({}).",
                    line_idx + 1,
                    hist.len(),
                    first.len()
                ));
            }
        }

        if hist.iter().any(|&b| b) {
            ids.push(cols[0].to_string());
            histories.push(hist);
        }
    }

    if histories.len() < 3 {
        return Err(String::from(
            "Capture-recapture analysis requires at least 3 captured individuals.",
        ));
    }
    if histories[0].len() < 2 {
        return Err(String::from(
            "Capture-recapture analysis requires at least 2 sampling occasions.",
        ));
    }
    Ok((ids, histories))
}

/// Runs closed-population abundance estimation across `M0` (equal catchability), `Mt` (Darroch time-varying),
/// `Chao (1987) Mh` (individual heterogeneity), and `Schnabel (1938)` multi-census models.
pub fn analyze_closed_capture_histories(
    histories: &[Vec<bool>],
) -> Result<ClosedCaptureSuiteReport, String> {
    if histories.len() < 3 || histories[0].len() < 2 {
        return Err(String::from(
            "Closed capture-recapture requires >= 3 individuals across >= 2 occasions.",
        ));
    }

    let mt1 = histories.len();
    let k = histories[0].len();
    let mut n_j = vec![0usize; k];
    let mut freq_f = vec![0usize; k];

    for h in histories {
        let mut c_count = 0usize;
        for (j, &caught) in h.iter().enumerate() {
            if caught {
                n_j[j] += 1;
                c_count += 1;
            }
        }
        if c_count >= 1 && c_count <= k {
            freq_f[c_count - 1] += 1;
        }
    }
    let total_captures: usize = n_j.iter().sum();

    // 1. M0 Equal-Catchability Exact Combinatorial Profile-Likelihood MLE
    let m0_ll = |n_cand: usize| -> f64 {
        let nf = n_cand as f64;
        let kf = k as f64;
        let p = ((total_captures as f64) / (kf * nf)).clamp(1e-12, 1.0 - 1e-12);
        ln_gamma(nf + 1.0) - ln_gamma((n_cand - mt1) as f64 + 1.0)
            + (total_captures as f64) * p.ln()
            + (kf * nf - (total_captures as f64)) * (1.0 - p).ln()
    };
    let (m0_n, m0_se, m0_ci) = profile_discrete_abundance(mt1, &m0_ll);
    let m0_p = (total_captures as f64) / ((k as f64) * m0_n);

    // 2. Mt Darroch Time-Varying Exact Combinatorial Profile-Likelihood MLE
    let min_n_mt = mt1.max(n_j.iter().copied().max().unwrap_or(mt1));
    let mt_ll = |n_cand: usize| -> f64 {
        let nf = n_cand as f64;
        let mut ll = ln_gamma(nf + 1.0) - ln_gamma((n_cand - mt1) as f64 + 1.0);
        for &nj in &n_j {
            let pj = ((nj as f64) / nf).clamp(1e-12, 1.0 - 1e-12);
            ll += (nj as f64) * pj.ln() + (nf - (nj as f64)) * (1.0 - pj).ln();
        }
        ll
    };
    let (mt_n, mt_se, mt_ci) = profile_discrete_abundance(min_n_mt, &mt_ll);
    let mt_p = (total_captures as f64) / ((k as f64) * mt_n);

    // 3. Chao (1987) Biometrics 43:783–791 Finite-Occasion Mh Heterogeneity Lower-Bound Estimator
    let f1 = freq_f[0] as f64;
    let f2 = if k >= 2 { freq_f[1] as f64 } else { 0.0 };
    let k_corr = ((k - 1) as f64) / (k as f64);
    let (mh_extra, mh_var) = if f1 == 0.0 {
        (0.0, 0.0)
    } else if f2 > 0.0 {
        let ex = k_corr * (f1 * f1) / (2.0 * f2);
        let v = f2
            * (0.5 * k_corr * (f1 / f2).powi(2)
                + k_corr.powi(2) * (f1 / f2).powi(3)
                + 0.25 * k_corr.powi(2) * (f1 / f2).powi(4));
        (ex, v.max(0.0))
    } else {
        let ex = k_corr * (f1 * (f1 - 1.0)) / 2.0;
        let v = 0.5 * k_corr * f1 * (f1 - 1.0) + 0.25 * k_corr.powi(2) * f1 * (2.0 * f1 - 1.0).powi(2);
        (ex, v.max(0.0))
    };
    let mh_n = (mt1 as f64) + mh_extra;
    let mh_se = mh_var.sqrt();
    let mh_ci = lognormal_abundance_ci(mt1 as f64, mh_n, mh_var);
    let mh_p = (total_captures as f64) / ((k as f64) * mh_n);

    // 4. Schnabel (1938) Chapman-Modified Multi-Census Estimator
    let mut marked_so_far = vec![false; mt1];
    let mut sum_cm = 0.0_f64;
    let mut sum_recap = 0usize;
    for j in 0..k {
        let m_prior = marked_so_far.iter().filter(|&&b| b).count();
        let mut r_j = 0usize;
        for (i, h) in histories.iter().enumerate() {
            if h[j] {
                if marked_so_far[i] {
                    r_j += 1;
                } else {
                    marked_so_far[i] = true;
                }
            }
        }
        if j >= 1 {
            sum_cm += (n_j[j] * m_prior) as f64;
            sum_recap += r_j;
        }
    }
    let schnabel_n = (sum_cm / ((sum_recap + 1) as f64)).max(mt1 as f64);
    let schnabel_var =
        ((sum_recap as f64) * schnabel_n * schnabel_n) / (((sum_recap + 1) * (sum_recap + 2)) as f64).max(1.0);
    let schnabel_ci = lognormal_abundance_ci(mt1 as f64, schnabel_n, schnabel_var);
    let schnabel_p = (total_captures as f64) / ((k as f64) * schnabel_n);

    Ok(ClosedCaptureSuiteReport {
        distinct_individuals_mt1: mt1,
        occasions_k: k,
        total_captures,
        captures_per_occasion: n_j,
        frequency_spectrum_f: freq_f,
        m0_equal_catchability: ClosedPopulationEstimate {
            model_name: String::from("M0 (Otis et al. 1978 Equal Catchability)"),
            n_hat: m0_n,
            se: m0_se,
            ci95: m0_ci,
            mean_capture_p: m0_p,
        },
        mt_darroch_time_varying: ClosedPopulationEstimate {
            model_name: String::from("Mt (Darroch 1958 Time-Varying Catchability)"),
            n_hat: mt_n,
            se: mt_se,
            ci95: mt_ci,
            mean_capture_p: mt_p,
        },
        mh_chao1987_heterogeneity: ClosedPopulationEstimate {
            model_name: String::from("Mh (Chao 1987 Individual Heterogeneity)"),
            n_hat: mh_n,
            se: mh_se,
            ci95: mh_ci,
            mean_capture_p: mh_p,
        },
        schnabel_multicensus: ClosedPopulationEstimate {
            model_name: String::from("Schnabel (1938 Chapman Multi-Census)"),
            n_hat: schnabel_n,
            se: schnabel_var.sqrt(),
            ci95: schnabel_ci,
            mean_capture_p: schnabel_p,
        },
    })
}

/// Estimates open-population Cormack–Jolly–Seber (`CJS`, Lebreton et al. 1992) apparent survival `phi`
/// and recapture probability `p` via exact Burnham `m`-array multinomial likelihood maximization.
pub fn estimate_cjs_apparent_survival(histories: &[Vec<bool>]) -> Result<CjsSurvivalReport, String> {
    if histories.is_empty() || histories[0].len() < 3 {
        return Err(String::from(
            "Open-population CJS estimation requires at least 3 capture occasions.",
        ));
    }
    let k = histories[0].len();
    let mut r_released = vec![0usize; k - 1];
    let mut m_mat = vec![vec![0usize; k]; k - 1];

    for h in histories {
        let mut prev_occ: Option<usize> = None;
        for (occ, &caught) in h.iter().enumerate() {
            if caught {
                if let Some(pi) = prev_occ {
                    m_mat[pi][occ] += 1;
                }
                if occ < k - 1 {
                    r_released[occ] += 1;
                }
                prev_occ = Some(occ);
            }
        }
    }

    let tot_rel: usize = r_released.iter().sum();
    let tot_recap: usize = m_mat.iter().flat_map(|r| r.iter()).copied().sum();
    if tot_rel == 0 {
        return Err(String::from(
            "No releases prior to the final occasion; CJS survival is unidentifiable.",
        ));
    }

    let mut m_array = Vec::with_capacity(k - 1);
    for i in 0..(k - 1) {
        let rec_sum: usize = m_mat[i][(i + 1)..k].iter().sum();
        let never = r_released[i].saturating_sub(rec_sum);
        m_array.push((r_released[i], m_mat[i][(i + 1)..k].to_vec(), never));
    }

    let neg_ll = |theta: &[f64]| -> f64 {
        if theta[0].abs() > 12.0 || theta[1].abs() > 12.0 {
            return 1e12 + theta[0] * theta[0] + theta[1] * theta[1];
        }
        let phi = 1.0 / (1.0 + (-theta[0]).exp());
        let p = 1.0 / (1.0 + (-theta[1]).exp());

        let mut ll = 0.0_f64;
        for (i, (_, row_rec, never)) in m_array.iter().enumerate() {
            let mut sum_pi = 0.0_f64;
            for (offset, &m_ij) in row_rec.iter().enumerate() {
                let dt = (offset + 1) as i32;
                let pi_ij = phi.powi(dt) * (1.0 - p).powi(dt - 1) * p;
                sum_pi += pi_ij;
                if m_ij > 0 {
                    ll += (m_ij as f64) * pi_ij.max(1e-15).ln();
                }
            }
            let chi_i = (1.0 - sum_pi).max(1e-15);
            if *never > 0 {
                ll += (*never as f64) * chi_i.ln();
            }
            let _ = i;
        }
        -ll
    };

    let (best_th, best_nll) = nelder_mead_minimize(&neg_ll, &[1.1, 0.0], 900, 1e-12);
    let se_logit = diag_hessian_se(&neg_ll, &best_th);

    let logistic = |x: f64| 1.0 / (1.0 + (-x).exp());
    let phi_hat = logistic(best_th[0]);
    let p_hat = logistic(best_th[1]);
    let phi_se = phi_hat * (1.0 - phi_hat) * se_logit[0];
    let p_se = p_hat * (1.0 - p_hat) * se_logit[1];

    let phi_ci95 = (
        logistic(best_th[0] - 1.96 * se_logit[0]),
        logistic(best_th[0] + 1.96 * se_logit[0]),
    );
    let p_ci95 = (
        logistic(best_th[1] - 1.96 * se_logit[1]),
        logistic(best_th[1] + 1.96 * se_logit[1]),
    );

    let ll = -best_nll;
    let eff_n = (tot_rel as f64).max(4.0);
    let aicc = 4.0 - 2.0 * ll + 12.0 / (eff_n - 3.0).max(1.0);

    Ok(CjsSurvivalReport {
        occasions_k: k,
        total_releases: tot_rel,
        total_recaptures: tot_recap,
        phi_hat,
        phi_se,
        phi_ci95,
        p_hat,
        p_se,
        p_ci95,
        log_likelihood: ll,
        aicc,
        m_array,
    })
}

/// Renders a two-panel `.svg` diagnostic figure pairing MacKenzie (2002) cumulative detection
/// `1 - (1 - p)^K` + occupancy `psi` (`left`) with closed-population abundance `N_hat` + 95% CIs (`right`).
pub fn render_wildlife_population_svg(
    occ: &OccupancySuiteReport,
    closed: &ClosedCaptureSuiteReport,
    title: &str,
) -> String {
    let width = 940.0_f64;
    let height = 430.0_f64;
    let p1_x = 68.0_f64;
    let p1_w = 385.0_f64;
    let p2_x = 540.0_f64;
    let p2_w = 350.0_f64;
    let top = 58.0_f64;
    let ph = 305.0_f64;

    let best_occ = &occ.models[occ.best_model_index];
    let mut svg = String::with_capacity(12_288);
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
        "  <text x=\"{p1_x:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — MacKenzie Occupancy &amp; Closed Capture-Recapture Abundance</text>",
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

    // Panel 1: Cumulative detection curve 1 - (1 - p_hat)^k for k = 1..=8 + horizontal dashed 0.95 line
    let y_prob = |pr: f64| top + ph - pr.clamp(0.0, 1.0) * ph;
    let y95 = y_prob(0.95);
    let _ = writeln!(
        svg,
        "  <line x1=\"{p1_x:.1}\" y1=\"{y95:.1}\" x2=\"{:.1}\" y2=\"{y95:.1}\" stroke=\"#dc2626\" stroke-dasharray=\"5,4\" stroke-width=\"1.4\"/>",
        p1_x + p1_w
    );
    let mut poly = String::new();
    for k_vis in 1..=8 {
        let xp = p1_x + ((k_vis - 1) as f64 / 7.0) * p1_w;
        let cp = 1.0 - (1.0 - best_occ.p_hat).powi(k_vis);
        let yp = y_prob(cp);
        let _ = write!(poly, "{xp:.1},{yp:.1} ");
        let _ = writeln!(svg, "  <circle cx=\"{xp:.1}\" cy=\"{yp:.1}\" r=\"4.0\" fill=\"#2563eb\"/>");
    }
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.5\"/>",
        poly.trim()
    );

    // Panel 2: Closed-population N_hat & 95% CI whiskers across M0, Mt, Mh, Schnabel + Mt+1 floor line
    let ests = [
        ("M0", &closed.m0_equal_catchability, "#2563eb"),
        ("Mt", &closed.mt_darroch_time_varying, "#16a34a"),
        ("Chao Mh", &closed.mh_chao1987_heterogeneity, "#9333ea"),
        ("Schnabel", &closed.schnabel_multicensus, "#d97706"),
    ];
    let max_n = ests
        .iter()
        .map(|(_, e, _)| e.ci95.1)
        .fold((closed.distinct_individuals_mt1 as f64) * 1.4, f64::max)
        * 1.08;
    let map_yn = |n_val: f64| top + ph - (n_val / max_n).clamp(0.0, 1.0) * ph;

    let y_mt1 = map_yn(closed.distinct_individuals_mt1 as f64);
    let _ = writeln!(
        svg,
        "  <line x1=\"{p2_x:.1}\" y1=\"{y_mt1:.1}\" x2=\"{:.1}\" y2=\"{y_mt1:.1}\" stroke=\"#64748b\" stroke-dasharray=\"4,4\" stroke-width=\"1.4\"/>",
        p2_x + p2_w
    );

    for (idx, (lbl, est, col)) in ests.iter().enumerate() {
        let xp = p2_x + 48.0 + (idx as f64) * ((p2_w - 96.0) / 3.0);
        let y_lo = map_yn(est.ci95.0);
        let y_hi = map_yn(est.ci95.1);
        let y_pt = map_yn(est.n_hat);
        let _ = writeln!(
            svg,
            "  <line x1=\"{xp:.1}\" y1=\"{y_lo:.1}\" x2=\"{xp:.1}\" y2=\"{y_hi:.1}\" stroke=\"{col}\" stroke-width=\"2.4\"/>"
        );
        let _ = writeln!(svg, "  <circle cx=\"{xp:.1}\" cy=\"{y_pt:.1}\" r=\"5.2\" fill=\"{col}\"/>");
        let _ = writeln!(
            svg,
            "  <text x=\"{xp:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#0f172a\">{lbl}</text>",
            top + ph + 18.0
        );
    }

    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Repeat Visits K — Cumulative Detection 1 - (1 - p)^K (psi={:.2}, naive={:.2})</text>",
        p1_x + 0.5 * p1_w,
        height - 16.0,
        best_occ.psi_hat,
        occ.naive_occupancy
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Closed Abundance Models (95% CI; dashed = M_t+1 = {})</text>",
        p2_x + 0.5 * p2_w,
        height - 16.0,
        closed.distinct_individuals_mt1
    );

    svg.push_str("</svg>\n");
    svg
}

fn profile_discrete_abundance<F: Fn(usize) -> f64>(
    min_n: usize,
    ll_fn: &F,
) -> (f64, f64, (f64, f64)) {
    let max_scan = (min_n * 12).max(min_n + 400);
    let mut best_n = min_n;
    let mut best_ll = f64::NEG_INFINITY;
    let mut vals = Vec::with_capacity(max_scan - min_n + 1);

    for n_cand in min_n..=max_scan {
        let v = ll_fn(n_cand);
        vals.push((n_cand, v));
        if v > best_ll {
            best_ll = v;
            best_n = n_cand;
        }
    }

    // Profile-likelihood 95% CI: `2 * (ll_max - ll(N)) <= 3.84146`
    let cutoff = best_ll - 0.5 * 3.84146;
    let ci_lo = vals
        .iter()
        .find(|(_, v)| *v >= cutoff)
        .map(|(n, _)| *n as f64)
        .unwrap_or(min_n as f64);
    let ci_hi = vals
        .iter()
        .rev()
        .find(|(_, v)| *v >= cutoff)
        .map(|(n, _)| *n as f64)
        .unwrap_or(best_n as f64);
    let approx_se = ((ci_hi - ci_lo) / 3.92).max(0.5);
    (best_n as f64, approx_se, (ci_lo, ci_hi))
}

fn lognormal_abundance_ci(min_observed: f64, n_hat: f64, var: f64) -> (f64, f64) {
    let f0 = (n_hat - min_observed).max(0.0);
    if f0 <= 1e-9 || var <= 1e-12 {
        return (min_observed, n_hat.max(min_observed));
    }
    let c = (1.96 * (1.0 + var / (f0 * f0)).ln().sqrt()).exp();
    (min_observed + f0 / c, min_observed + f0 * c)
}

/// Lanczos approximation to `ln(Gamma(z))` accurate to ~15 digits for `z > 0`.
fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1_259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    let zm1 = z - 1.0;
    let mut x = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        x += c / (zm1 + (i as f64));
    }
    let t = zm1 + 7.5;
    0.5 * (2.0 * std::f64::consts::PI).ln() + (zm1 + 0.5) * t.ln() - t + x.ln()
}


/// Built-in 6-occasion PIT-tag capture-history dataset (`36` Coastal Giant Salamanders
/// *Dicamptodon tenebrosus* in Oak Creek, Benton County, Oregon).
pub fn demo_pnw_capture_history_csv() -> &'static str {
    "\
# Oak Creek (Corvallis, OR) PIT-Tagged Dicamptodon tenebrosus Capture Histories (36 Individuals x 6 Occasions)
animal_id,occ_1,occ_2,occ_3,occ_4,occ_5,occ_6
PIT_001,1,1,0,1,0,1
PIT_002,1,0,1,1,0,0
PIT_003,0,1,1,0,1,0
PIT_004,1,0,0,1,0,1
PIT_005,0,1,0,1,1,0
PIT_006,1,1,0,0,1,0
PIT_007,0,0,1,1,0,1
PIT_008,1,1,1,0,0,0
PIT_009,0,1,1,1,0,0
PIT_010,1,1,0,1,0,0
PIT_011,0,0,1,1,1,0
PIT_012,0,1,1,0,0,1
PIT_013,1,1,0,0,0,0
PIT_014,0,0,0,1,1,0
PIT_015,0,0,1,1,0,1
PIT_016,1,0,0,0,1,0
PIT_017,0,1,1,0,0,0
PIT_018,0,0,0,0,1,1
PIT_019,1,0,0,0,0,0
PIT_020,1,0,0,0,0,0
PIT_021,1,0,0,0,0,0
PIT_022,0,1,0,0,0,0
PIT_023,0,1,0,0,0,0
PIT_024,0,1,0,0,0,0
PIT_025,0,0,1,0,0,0
PIT_026,0,0,1,0,0,0
PIT_027,0,0,1,0,0,0
PIT_028,0,0,0,1,0,0
PIT_029,0,0,0,1,0,0
PIT_030,0,0,0,1,0,0
PIT_031,0,0,0,0,1,0
PIT_032,0,0,0,0,1,0
PIT_033,0,0,0,0,1,0
PIT_034,0,0,0,0,0,1
PIT_035,0,0,0,0,0,1
PIT_036,0,0,0,0,0,1
"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_chapman_closed_m0_mt_mh_and_open_cjs() {
        let cp = calc_chapman_lincoln_petersen(40, 35, 14).unwrap();
        assert!((cp.chapman_n_hat - 97.4).abs() < 0.01);
        assert!(cp.ci95.0 >= (cp.distinct_observed_mt1 as f64));

        let (_, hists) = parse_capture_history_csv(demo_pnw_capture_history_csv()).unwrap();
        let closed = analyze_closed_capture_histories(&hists).unwrap();
        assert_eq!(closed.distinct_individuals_mt1, 36);
        assert!(closed.m0_equal_catchability.n_hat >= 36.0);
        // Individual heterogeneity gives Chao Mh >= M0 when singletons f1 = 18 are abundant
        assert!(closed.mh_chao1987_heterogeneity.n_hat > closed.m0_equal_catchability.n_hat);

        let cjs = estimate_cjs_apparent_survival(&hists).unwrap();
        assert!(cjs.phi_hat > 0.50 && cjs.phi_hat < 0.99);
        assert!(cjs.p_hat > 0.15 && cjs.p_hat < 0.80);
    }
}
