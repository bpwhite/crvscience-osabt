// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #15 (`osabt wildlife`) — Imperfect-detection site-occupancy modeling & survey-design engine:
//! - **Single-Season Occupancy MLE (`Psi(.)p(.)`, `Psi(cov)p(.)`, `Psi(.)p(cov)`)** —
//!   MacKenzie et al. (2002) *Ecology* 83(8):2248–2255 (`https://pubs.usgs.gov/publication/5224176`)
//! - **Empirical Bayes Conditional Site Occupancy & Burnham–Anderson `AICc` Ranking**
//! - **Optimal Survey Effort (`K*`, `S*`) & False-Absence Power Planner** —
//!   MacKenzie & Royle (2005) *J. Appl. Ecol.* 42(6):1105–1114 (`https://pubs.usgs.gov/publication/5224604`).

/// Detection/non-detection survey record for one spatial site across `K` secondary survey occasions.
#[derive(Debug, Clone)]
pub struct SiteDetectionRow {
    pub site_id: String,
    /// `Some(true)` = detected (`1`), `Some(false)` = undetected (`0`), `None` = missing visit (`-` / `NA`).
    pub visits: Vec<Option<bool>>,
    /// Optional site-level habitat covariate (e.g., Old-Growth Canopy Index or Elevation).
    pub covariate: f64,
}

impl SiteDetectionRow {
    pub fn valid_visits(&self) -> usize {
        self.visits.iter().filter(|v| v.is_some()).count()
    }

    pub fn detections(&self) -> usize {
        self.visits
            .iter()
            .filter(|v| matches!(v, Some(true)))
            .count()
    }
}

/// Candidate parameterizations for single-season MacKenzie et al. (2002) occupancy fitting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum OccupancyCovariateModel {
    /// Constant occupancy `psi` and constant per-visit detection `p` (`K = 2` parameters).
    ConstantPsiConstantP,
    /// Logistic habitat covariate on `psi(x)`, constant per-visit detection `p` (`K = 3` parameters).
    CovariatePsiConstantP,
    /// Constant occupancy `psi`, logistic covariate on `p(x)` (`K = 3` parameters).
    ConstantPsiCovariateP,
}

impl OccupancyCovariateModel {
    pub fn label(self) -> &'static str {
        match self {
            Self::ConstantPsiConstantP => "Psi(.) p(.)",
            Self::CovariatePsiConstantP => "Psi(cov) p(.)",
            Self::ConstantPsiCovariateP => "Psi(.) p(cov)",
        }
    }

    pub fn num_params(self) -> usize {
        match self {
            Self::ConstantPsiConstantP => 2,
            Self::CovariatePsiConstantP | Self::ConstantPsiCovariateP => 3,
        }
    }
}

/// Fitted maximum-likelihood occupancy model with asymptotic Hessian SEs and `AICc` diagnostics.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OccupancyModelFit {
    pub model: OccupancyCovariateModel,
    /// Mean/intercept site occupancy probability `psi` at mean covariate (`z_cov = 0`).
    pub psi_hat: f64,
    pub psi_se: f64,
    pub psi_ci95: (f64, f64),
    /// Mean/intercept per-visit detection probability `p` at mean covariate (`z_cov = 0`).
    pub p_hat: f64,
    pub p_se: f64,
    pub p_ci95: (f64, f64),
    /// Standardized logit-slope `beta_1` (on `psi` or `p`, `0.0` for `Psi(.)p(.)`).
    pub covariate_logit_slope: f64,
    pub covariate_slope_se: f64,
    pub log_likelihood: f64,
    pub aicc: f64,
    pub delta_aicc: f64,
    pub akaike_weight: f64,
}

/// Empirical Bayes posterior occupancy status for an individual surveyed site.
#[derive(Debug, Clone)]
pub struct SitePosteriorOccupancy {
    pub site_id: String,
    pub valid_visits_k: usize,
    pub detections_y: usize,
    pub covariate: f64,
    pub fitted_psi_i: f64,
    pub fitted_p_i: f64,
    /// `1.0` if `detections_y > 0`, else `psi_i * (1 - p_i)^k / (psi_i * (1 - p_i)^k + (1 - psi_i))`.
    pub posterior_occupied_prob: f64,
}

/// Complete single-season MacKenzie et al. (2002) occupancy report across candidate models.
#[derive(Debug, Clone)]
pub struct OccupancySuiteReport {
    pub num_sites: usize,
    pub max_occasions: usize,
    pub naive_occupied_sites: usize,
    pub naive_occupancy: f64,
    /// Sum of Empirical Bayes conditional occupancy probabilities across all `S` sites.
    pub finite_sample_occupied_est: f64,
    pub models: Vec<OccupancyModelFit>,
    pub best_model_index: usize,
    pub site_posteriors: Vec<SitePosteriorOccupancy>,
}

/// MacKenzie & Royle (2005) survey-design & false-absence power planner output.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OccupancySurveyPlan {
    pub assumed_psi: f64,
    pub assumed_p: f64,
    pub target_false_absence_alpha: f64,
    /// Minimum visits per site `K_min = ceil(ln(alpha) / ln(1 - p))` so `(1 - p)^K <= alpha`.
    pub min_visits_for_presence_certainty: usize,
    pub cumulative_detection_at_k_min: f64,
    /// Optimal repeat visits `K*` minimizing asymptotic `Var(psi_hat)` per unit survey cost.
    pub optimal_k_for_precision: usize,
    /// Number of sites `S*` required at `K*` to achieve `SE(psi_hat) <= target_psi_se`.
    pub required_sites_for_target_se: usize,
    pub achieved_psi_se: f64,
    pub total_survey_visits: usize,
}

/// Parses a comma- or tab-delimited occupancy detection-history table:
/// `site_id,v1,v2,...,vK[,covariate]` where visits are `1`, `0`, or `-`/`NA` for missing surveys.
pub fn parse_occupancy_csv(raw_csv: &str) -> Result<Vec<SiteDetectionRow>, String> {
    let mut rows = Vec::new();
    let mut header_has_cov = false;
    let mut expected_cols: Option<usize> = None;

    for (line_idx, line) in raw_csv.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let sep = if trimmed.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = trimmed.split(sep).map(|s| s.trim()).collect();

        if expected_cols.is_none() {
            let first_low = cols[0].to_ascii_lowercase();
            if first_low == "site" || first_low == "site_id" {
                if cols.len() < 3 {
                    return Err(String::from(
                        "Occupancy header must have site_id and at least 2 survey columns.",
                    ));
                }
                let last_low = cols.last().copied().unwrap_or("").to_ascii_lowercase();
                header_has_cov = last_low.contains("cov")
                    || last_low.contains("canopy")
                    || last_low.contains("elev")
                    || last_low.contains("habitat");
                expected_cols = Some(cols.len());
                continue;
            }
            expected_cols = Some(cols.len());
        }

        if cols.len() < 3 {
            return Err(format!(
                "Occupancy CSV line {}: expected >= 3 columns, found {}.",
                line_idx + 1,
                cols.len()
            ));
        }

        let visit_end = if header_has_cov {
            cols.len() - 1
        } else {
            // Auto-detect trailing continuous covariate if not `0`, `1`, `-`, `NA`
            let last = cols[cols.len() - 1];
            if matches!(last, "0" | "1" | "-" | "NA" | "na" | ".") {
                cols.len()
            } else {
                cols.len() - 1
            }
        };

        if visit_end <= 1 {
            return Err(format!(
                "Occupancy CSV line {}: no detection-history columns found.",
                line_idx + 1
            ));
        }

        let mut visits = Vec::with_capacity(visit_end - 1);
        for &tok in &cols[1..visit_end] {
            match tok {
                "1" => visits.push(Some(true)),
                "0" => visits.push(Some(false)),
                "-" | "NA" | "na" | "." | "" => visits.push(None),
                other => {
                    return Err(format!(
                        "Occupancy CSV line {}: invalid detection code '{}' (expected 0, 1, or NA).",
                        line_idx + 1,
                        other
                    ));
                }
            }
        }

        let covariate = if visit_end < cols.len() {
            cols[visit_end].parse::<f64>().map_err(|_| {
                format!(
                    "Occupancy CSV line {}: invalid numeric covariate '{}'.",
                    line_idx + 1,
                    cols[visit_end]
                )
            })?
        } else {
            0.0
        };

        rows.push(SiteDetectionRow {
            site_id: cols[0].to_string(),
            visits,
            covariate,
        });
    }

    if rows.len() < 4 {
        return Err(String::from(
            "Occupancy MLE requires at least 4 surveyed sites.",
        ));
    }
    Ok(rows)
}

/// Fits MacKenzie et al. (2002) single-season imperfect-detection occupancy models (`Psi(.)p(.)`,
/// `Psi(cov)p(.)`, `Psi(.)p(cov)`), ranks them via `AICc`, and computes Empirical Bayes site posteriors.
pub fn fit_mackenzie_occupancy_suite(
    sites: &[SiteDetectionRow],
) -> Result<OccupancySuiteReport, String> {
    if sites.len() < 4 {
        return Err(String::from("At least 4 sites are required for occupancy fitting."));
    }

    let s_count = sites.len();
    let max_k = sites.iter().map(|r| r.visits.len()).max().unwrap_or(0);
    let naive_occ_sites = sites.iter().filter(|r| r.detections() > 0).count();
    let naive_occ = (naive_occ_sites as f64) / (s_count as f64);

    // Standardize site covariate to mean 0, SD 1 for numerical stability
    let cov_mean = sites.iter().map(|r| r.covariate).sum::<f64>() / (s_count as f64);
    let cov_var = sites
        .iter()
        .map(|r| (r.covariate - cov_mean).powi(2))
        .sum::<f64>()
        / ((s_count.saturating_sub(1)).max(1) as f64);
    let has_cov_variation = cov_var > 1e-10;
    let cov_sd = cov_var.sqrt().max(1.0);
    let z_covs: Vec<f64> = sites
        .iter()
        .map(|r| (r.covariate - cov_mean) / cov_sd)
        .collect();

    let mut candidates = vec![OccupancyCovariateModel::ConstantPsiConstantP];
    if has_cov_variation {
        candidates.push(OccupancyCovariateModel::CovariatePsiConstantP);
        candidates.push(OccupancyCovariateModel::ConstantPsiCovariateP);
    }

    let mut fits = Vec::with_capacity(candidates.len());
    for cand in candidates {
        fits.push(fit_single_occupancy_model(sites, &z_covs, cand));
    }

    let min_aicc = fits
        .iter()
        .map(|f| f.aicc)
        .fold(f64::INFINITY, f64::min);
    let mut weight_sum = 0.0_f64;
    for fit in &mut fits {
        fit.delta_aicc = (fit.aicc - min_aicc).max(0.0);
        let w = (-0.5 * fit.delta_aicc).exp();
        fit.akaike_weight = w;
        weight_sum += w;
    }
    if weight_sum > 0.0 {
        for fit in &mut fits {
            fit.akaike_weight /= weight_sum;
        }
    }

    let best_idx = fits
        .iter()
        .enumerate()
        .min_by(|(_, a), (_, b)| a.aicc.partial_cmp(&b.aicc).unwrap_or(std::cmp::Ordering::Equal))
        .map(|(i, _)| i)
        .unwrap_or(0);

    let best = &fits[best_idx];
    let b0_psi = logit(best.psi_hat);
    let b0_p = logit(best.p_hat);

    let mut posteriors = Vec::with_capacity(s_count);
    let mut finite_sum = 0.0_f64;

    for (idx, row) in sites.iter().enumerate() {
        let zc = z_covs[idx];
        let psi_i = match best.model {
            OccupancyCovariateModel::CovariatePsiConstantP => {
                logistic(b0_psi + best.covariate_logit_slope * zc)
            }
            _ => best.psi_hat,
        };
        let p_i = match best.model {
            OccupancyCovariateModel::ConstantPsiCovariateP => {
                logistic(b0_p + best.covariate_logit_slope * zc)
            }
            _ => best.p_hat,
        };
        let k_i = row.valid_visits();
        let y_i = row.detections();
        let post = if y_i > 0 {
            1.0
        } else {
            let miss_all = (1.0 - p_i).powi(k_i as i32);
            let num = psi_i * miss_all;
            let den = num + (1.0 - psi_i);
            if den > 0.0 { num / den } else { 0.0 }
        };
        finite_sum += post;
        posteriors.push(SitePosteriorOccupancy {
            site_id: row.site_id.clone(),
            valid_visits_k: k_i,
            detections_y: y_i,
            covariate: row.covariate,
            fitted_psi_i: psi_i,
            fitted_p_i: p_i,
            posterior_occupied_prob: post,
        });
    }

    Ok(OccupancySuiteReport {
        num_sites: s_count,
        max_occasions: max_k,
        naive_occupied_sites: naive_occ_sites,
        naive_occupancy: naive_occ,
        finite_sample_occupied_est: finite_sum,
        models: fits,
        best_model_index: best_idx,
        site_posteriors: posteriors,
    })
}

fn fit_single_occupancy_model(
    sites: &[SiteDetectionRow],
    z_covs: &[f64],
    model: OccupancyCovariateModel,
) -> OccupancyModelFit {
    let n_par = model.num_params();
    let neg_ll = |theta: &[f64]| -> f64 {
        // Penalize extreme logits softly to keep numerical optimization well-conditioned
        if theta.iter().any(|&v| v.abs() > 14.0) {
            return 1e12 + theta.iter().map(|v| v * v).sum::<f64>();
        }
        let mut ll = 0.0_f64;
        for (i, row) in sites.iter().enumerate() {
            let zc = z_covs[i];
            let (psi_i, p_i) = match model {
                OccupancyCovariateModel::ConstantPsiConstantP => {
                    (logistic(theta[0]), logistic(theta[1]))
                }
                OccupancyCovariateModel::CovariatePsiConstantP => {
                    (logistic(theta[0] + theta[2] * zc), logistic(theta[1]))
                }
                OccupancyCovariateModel::ConstantPsiCovariateP => {
                    (logistic(theta[0]), logistic(theta[1] + theta[2] * zc))
                }
            };
            let k_i = row.valid_visits();
            let y_i = row.detections();
            if k_i == 0 {
                continue;
            }
            if y_i > 0 {
                ll += psi_i.max(1e-15).ln()
                    + (y_i as f64) * p_i.max(1e-15).ln()
                    + ((k_i - y_i) as f64) * (1.0 - p_i).max(1e-15).ln();
            } else {
                let p_und = psi_i * (1.0 - p_i).powi(k_i as i32) + (1.0 - psi_i);
                ll += p_und.max(1e-15).ln();
            }
        }
        -ll
    };

    let starts: Vec<Vec<f64>> = if n_par == 2 {
        vec![
            vec![0.2, -0.2],
            vec![0.8, 0.4],
            vec![-0.5, -0.8],
        ]
    } else {
        vec![
            vec![0.2, -0.2, 0.0],
            vec![0.5, 0.0, 0.8],
            vec![0.0, -0.4, -0.8],
        ]
    };

    let mut best_theta = starts[0].clone();
    let mut best_nll = f64::INFINITY;
    for st in &starts {
        let (th, val) = nelder_mead_minimize(&neg_ll, st, 800, 1e-10);
        if val < best_nll {
            best_nll = val;
            best_theta = th;
        }
    }
    // Polish from best solution
    let (best_theta, best_nll) = nelder_mead_minimize(&neg_ll, &best_theta, 900, 1e-12);

    let se_logit = diag_hessian_se(&neg_ll, &best_theta);

    let psi_hat = logistic(best_theta[0]);
    let p_hat = logistic(best_theta[1]);
    // Delta-method SE on probability scale: d(logistic(x))/dx = psi * (1 - psi)
    let psi_se = psi_hat * (1.0 - psi_hat) * se_logit[0];
    let p_se = p_hat * (1.0 - p_hat) * se_logit[1];

    let psi_ci95 = (
        logistic(best_theta[0] - 1.96 * se_logit[0]),
        logistic(best_theta[0] + 1.96 * se_logit[0]),
    );
    let p_ci95 = (
        logistic(best_theta[1] - 1.96 * se_logit[1]),
        logistic(best_theta[1] + 1.96 * se_logit[1]),
    );

    let (slope, slope_se) = if n_par == 3 {
        (best_theta[2], se_logit[2])
    } else {
        (0.0, 0.0)
    };

    let ll = -best_nll;
    let k_f = n_par as f64;
    let n_f = sites.len() as f64;
    let denom = (n_f - k_f - 1.0).max(1.0);
    let aicc = 2.0 * k_f - 2.0 * ll + (2.0 * k_f * (k_f + 1.0)) / denom;

    OccupancyModelFit {
        model,
        psi_hat,
        psi_se,
        psi_ci95,
        p_hat,
        p_se,
        p_ci95,
        covariate_logit_slope: slope,
        covariate_slope_se: slope_se,
        log_likelihood: ll,
        aicc,
        delta_aicc: 0.0,
        akaike_weight: 0.0,
    }
}

/// Computes MacKenzie & Royle (2005) analytical survey-effort recommendations:
/// 1. Minimum repeat visits `K_min` to drive false-absence probability `(1 - p)^K <= target_false_absence_alpha`
/// 2. Optimal repeat visits `K*` minimizing `Var(psi_hat)` per survey-budget dollar (`cost_ratio_site_to_visit`)
/// 3. Required number of sites `S*` at `K*` to attain `SE(psi_hat) <= target_psi_se`.
pub fn plan_occupancy_survey_effort(
    assumed_psi: f64,
    assumed_p: f64,
    target_false_absence_alpha: f64,
    target_psi_se: f64,
    cost_ratio_site_to_visit: f64,
) -> Result<OccupancySurveyPlan, String> {
    if !(0.01..=0.99).contains(&assumed_psi) {
        return Err(String::from("assumed_psi must lie in [0.01, 0.99]."));
    }
    if !(0.02..=0.99).contains(&assumed_p) {
        return Err(String::from("assumed_p must lie in [0.02, 0.99]."));
    }
    let alpha = target_false_absence_alpha.clamp(0.001, 0.50);
    let target_se = target_psi_se.clamp(0.01, 0.25);
    let c_ratio = cost_ratio_site_to_visit.max(0.0);

    let k_min = ((alpha.ln() / (1.0 - assumed_p).ln()).ceil() as usize).max(2);
    let cum_p_kmin = 1.0 - (1.0 - assumed_p).powi(k_min as i32);

    // MacKenzie & Royle (2005) Eq. 2 variance factor per site `V_1(K) = S * Var(psi_hat)`:
    // `V_1(K) = psi * [ (1 - psi) + (1 - p*) / (p* - K * p * (1 - p)^(K - 1)) ]`
    let site_var_factor = |k: usize| -> f64 {
        let p_star = 1.0 - (1.0 - assumed_p).powi(k as i32);
        let denom = (p_star - (k as f64) * assumed_p * (1.0 - assumed_p).powi((k - 1) as i32))
            .max(1e-12);
        assumed_psi * ((1.0 - assumed_psi) + (1.0 - p_star) / denom)
    };

    // Minimize total budget-weighted variance `(c_ratio + K) * V_1(K)` over K in 2..=25
    let mut opt_k = 2usize;
    let mut best_obj = f64::INFINITY;
    for k in 2..=25 {
        let obj = (c_ratio + (k as f64)) * site_var_factor(k);
        if obj < best_obj {
            best_obj = obj;
            opt_k = k;
        }
    }

    let v1_opt = site_var_factor(opt_k);
    let req_s = ((v1_opt / (target_se * target_se)).ceil() as usize).max(10);
    let achieved_se = (v1_opt / (req_s as f64)).sqrt();

    Ok(OccupancySurveyPlan {
        assumed_psi,
        assumed_p,
        target_false_absence_alpha: alpha,
        min_visits_for_presence_certainty: k_min,
        cumulative_detection_at_k_min: cum_p_kmin,
        optimal_k_for_precision: opt_k,
        required_sites_for_target_se: req_s,
        achieved_psi_se: achieved_se,
        total_survey_visits: req_s * opt_k,
    })
}

#[inline]
fn logistic(x: f64) -> f64 {
    if x >= 0.0 {
        1.0 / (1.0 + (-x).exp())
    } else {
        let ex = x.exp();
        ex / (1.0 + ex)
    }
}

#[inline]
fn logit(p: f64) -> f64 {
    let pc = p.clamp(1e-9, 1.0 - 1e-9);
    (pc / (1.0 - pc)).ln()
}

pub(crate) fn nelder_mead_minimize<F: Fn(&[f64]) -> f64>(
    f: &F,
    x0: &[f64],
    max_iter: usize,
    tol: f64,
) -> (Vec<f64>, f64) {
    let n = x0.len();
    let mut simplex: Vec<(Vec<f64>, f64)> = Vec::with_capacity(n + 1);
    simplex.push((x0.to_vec(), f(x0)));
    for i in 0..n {
        let mut pt = x0.to_vec();
        pt[i] += if pt[i].abs() > 0.1 { 0.25 * pt[i].abs() } else { 0.35 };
        let val = f(&pt);
        simplex.push((pt, val));
    }

    for _ in 0..max_iter {
        simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
        if (simplex[n].1 - simplex[0].1).abs() < tol {
            break;
        }
        let mut centroid = vec![0.0; n];
        for item in &simplex[..n] {
            for j in 0..n {
                centroid[j] += item.0[j] / (n as f64);
            }
        }
        let xr: Vec<f64> = (0..n)
            .map(|j| centroid[j] + (centroid[j] - simplex[n].0[j]))
            .collect();
        let fr = f(&xr);
        if fr < simplex[0].1 {
            let xe: Vec<f64> = (0..n)
                .map(|j| centroid[j] + 2.0 * (xr[j] - centroid[j]))
                .collect();
            let fe = f(&xe);
            simplex[n] = if fe < fr { (xe, fe) } else { (xr, fr) };
        } else if fr < simplex[n - 1].1 {
            simplex[n] = (xr, fr);
        } else {
            let xc: Vec<f64> = if fr < simplex[n].1 {
                (0..n).map(|j| centroid[j] + 0.5 * (xr[j] - centroid[j])).collect()
            } else {
                (0..n)
                    .map(|j| centroid[j] + 0.5 * (simplex[n].0[j] - centroid[j]))
                    .collect()
            };
            let fc = f(&xc);
            if fc < simplex[n].1.min(fr) {
                simplex[n] = (xc, fc);
            } else {
                let best = simplex[0].0.clone();
                for item in &mut simplex[1..] {
                    for j in 0..n {
                        item.0[j] = best[j] + 0.5 * (item.0[j] - best[j]);
                    }
                    item.1 = f(&item.0);
                }
            }
        }
    }
    simplex.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal));
    simplex.remove(0)
}

pub(crate) fn diag_hessian_se<F: Fn(&[f64]) -> f64>(neg_ll: &F, theta: &[f64]) -> Vec<f64> {
    let n = theta.len();
    let h = 1e-4_f64;
    let mut hess = vec![vec![0.0_f64; n]; n];
    let f0 = neg_ll(theta);

    for i in 0..n {
        let mut tp = theta.to_vec();
        let mut tm = theta.to_vec();
        tp[i] += h;
        tm[i] -= h;
        hess[i][i] = ((neg_ll(&tp) - 2.0 * f0 + neg_ll(&tm)) / (h * h)).max(1e-6);
        for j in (i + 1)..n {
            let mut tpp = theta.to_vec();
            let mut tpm = theta.to_vec();
            let mut tmp = theta.to_vec();
            let mut tmm = theta.to_vec();
            tpp[i] += h;
            tpp[j] += h;
            tpm[i] += h;
            tpm[j] -= h;
            tmp[i] -= h;
            tmp[j] += h;
            tmm[i] -= h;
            tmm[j] -= h;
            let hij = (neg_ll(&tpp) - neg_ll(&tpm) - neg_ll(&tmp) + neg_ll(&tmm)) / (4.0 * h * h);
            hess[i][j] = hij;
            hess[j][i] = hij;
        }
    }

    invert_small_spd(&hess)
        .map(|inv| inv.iter().enumerate().map(|(i, r)| r[i].max(1e-12).sqrt()).collect())
        .unwrap_or_else(|| (0..n).map(|i| (1.0 / hess[i][i].max(1e-6)).sqrt()).collect())
}

fn invert_small_spd(a: &[Vec<f64>]) -> Option<Vec<Vec<f64>>> {
    let n = a.len();
    let mut aug = vec![vec![0.0_f64; 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n + i] = 1.0;
    }
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if aug[r][col].abs() > aug[piv][col].abs() {
                piv = r;
            }
        }
        if aug[piv][col].abs() < 1e-10 {
            return None;
        }
        aug.swap(col, piv);
        let d = aug[col][col];
        for j in 0..(2 * n) {
            aug[col][j] /= d;
        }
        for r in 0..n {
            if r != col {
                let fac = aug[r][col];
                for j in 0..(2 * n) {
                    aug[r][j] -= fac * aug[col][j];
                }
            }
        }
    }
    let inv: Vec<Vec<f64>> = (0..n).map(|i| aug[i][n..(2 * n)].to_vec()).collect();
    if inv.iter().enumerate().all(|(i, r)| r[i] > 0.0 && r[i].is_finite()) {
        Some(inv)
    } else {
        None
    }
}

/// Built-in 24-site × 5-visit Oregon Coast Range **Red Tree Vole (`Arborimus longicaudus`)**
/// canopy occupancy dataset with old-growth canopy index covariate and realistic missing visits (`NA`).
pub fn demo_pnw_occupancy_csv() -> &'static str {
    "\
# Oregon Coast Range Red Tree Vole (Arborimus longicaudus) Single-Season Canopy Survey (24 Sites x 5 Visits)
site_id,visit_1,visit_2,visit_3,visit_4,visit_5,old_growth_canopy_idx
Siuslaw_01,1,0,1,0,1,0.92
Siuslaw_02,0,1,0,1,0,0.88
Siuslaw_03,1,1,0,0,1,0.95
Siuslaw_04,0,0,1,0,NA,0.84
Siuslaw_05,0,1,1,0,0,0.86
Siuslaw_06,0,0,0,1,0,0.79
MarysPeak_07,1,0,0,1,1,0.91
MarysPeak_08,0,0,0,0,0,0.76
MarysPeak_09,0,1,0,0,1,0.83
MarysPeak_10,1,0,1,1,0,0.94
MarysPeak_11,0,0,1,0,0,0.74
MarysPeak_12,0,0,0,0,NA,0.68
Alsea_13,0,1,0,0,0,0.72
Alsea_14,0,0,0,0,0,0.58
Alsea_15,1,0,0,1,0,0.81
Alsea_16,0,0,0,0,0,0.44
Alsea_17,0,0,0,0,0,0.39
Alsea_18,0,0,1,0,0,0.69
Nestucca_19,0,0,0,0,0,0.35
Nestucca_20,0,0,0,0,0,0.28
Nestucca_21,0,0,0,0,0,0.31
Nestucca_22,0,1,0,0,NA,0.77
Nestucca_23,0,0,0,0,0,0.22
Nestucca_24,0,0,0,0,0,0.19
"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_mackenzie_occupancy_above_naive_and_ranks_covariate_model() {
        let rows = parse_occupancy_csv(demo_pnw_occupancy_csv()).unwrap();
        let rep = fit_mackenzie_occupancy_suite(&rows).unwrap();
        let const_fit = &rep.models[0];

        assert!(const_fit.psi_hat > rep.naive_occupancy);
        assert!(const_fit.p_hat > 0.25 && const_fit.p_hat < 0.75);
        assert!(rep.finite_sample_occupied_est >= (rep.naive_occupied_sites as f64));
        // Positive old-growth canopy association in Psi(cov)p(.)
        assert!(rep.models[1].covariate_logit_slope > 0.5);

        let plan = plan_occupancy_survey_effort(0.60, 0.40, 0.05, 0.06, 2.0).unwrap();
        assert_eq!(plan.min_visits_for_presence_certainty, 6);
        assert!(plan.achieved_psi_se <= 0.06 + 1e-9);
    }
}
