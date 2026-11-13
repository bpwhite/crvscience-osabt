// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #19 (`osabt survival`) — Right-censored time-to-event survival engine:
//! - **Kaplan & Meier (1958)** product-limit survival curve with Greenwood (1926) SE & log-log 95% CIs
//! - **Mantel (1966)** two-group log-rank hypothesis test (`https://pubmed.ncbi.nlm.nih.gov/5910392/`)
//! - **Cox (1972)** proportional-hazards partial-likelihood regression (`https://academic.oup.com/jrsssb/article/34/2/187/7027194`).

use std::collections::BTreeMap;

/// Single right-censored time-to-event observation.
#[derive(Debug, Clone)]
pub struct SurvivalObservation {
    pub subject_id: String,
    pub group: String,
    pub time: f64,
    /// `true` = event occurred, `false` = right-censored.
    pub event: bool,
    pub covariate: f64,
}

/// One step in a Kaplan–Meier product-limit survival curve.
#[derive(Debug, Clone)]
pub struct KaplanMeierStep {
    pub time: f64,
    pub n_at_risk: usize,
    pub n_events: usize,
    pub n_censored: usize,
    pub survival_s: f64,
    pub greenwood_se: f64,
    pub ci95: (f64, f64),
}

/// Complete Kaplan–Meier curve summary for one cohort stratum.
#[derive(Debug, Clone)]
pub struct KaplanMeierCurve {
    pub group: String,
    pub total_n: usize,
    pub total_events: usize,
    pub median_survival_time: Option<f64>,
    pub steps: Vec<KaplanMeierStep>,
}

/// Mantel (1966) two-group log-rank test result.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LogRankTestResult {
    pub group_a: String,
    pub group_b: String,
    pub observed_a: f64,
    pub expected_a: f64,
    pub observed_b: f64,
    pub expected_b: f64,
    pub chi_square: f64,
    pub p_value: f64,
}

/// Cox (1972) single-covariate proportional-hazards partial-likelihood fit (Breslow ties).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CoxPhFit {
    pub covariate_name: String,
    pub n_subjects: usize,
    pub n_events: usize,
    pub beta: f64,
    pub se_beta: f64,
    pub hazard_ratio: f64,
    pub hr_ci95: (f64, f64),
    pub wald_z: f64,
    pub p_value: f64,
    pub partial_log_lik_null: f64,
    pub partial_log_lik_fit: f64,
}

/// Parses a right-censored survival CSV/TSV (`id,group,time,event[,covariate]`).
pub fn parse_survival_csv(raw_csv: &str) -> Result<Vec<SurvivalObservation>, String> {
    let mut out = Vec::new();
    let mut header_seen = false;

    for (line_idx, line) in raw_csv.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();
        if !header_seen {
            header_seen = true;
            let c0 = cols[0].to_ascii_lowercase();
            if c0.contains("id") || c0.contains("subject") || c0.contains("sample") {
                continue;
            }
        }
        if cols.len() < 4 {
            return Err(format!(
                "Survival CSV line {}: expected >= 4 columns (id,group,time,event[,covariate]).",
                line_idx + 1
            ));
        }
        let time = cols[2]
            .parse::<f64>()
            .map_err(|_| format!("Survival line {}: invalid time '{}'.", line_idx + 1, cols[2]))?;
        if time < 0.0 || !time.is_finite() {
            return Err(format!("Survival line {}: time must be finite and >= 0.", line_idx + 1));
        }
        let event = match cols[3] {
            "1" | "true" | "TRUE" | "event" | "dead" => true,
            "0" | "false" | "FALSE" | "censored" | "alive" => false,
            other => {
                return Err(format!(
                    "Survival line {}: invalid event indicator '{}' (expected 0 or 1).",
                    line_idx + 1,
                    other
                ));
            }
        };
        let covariate = if cols.len() >= 5 {
            cols[4].parse::<f64>().unwrap_or(0.0)
        } else {
            0.0
        };
        out.push(SurvivalObservation {
            subject_id: cols[0].to_string(),
            group: cols[1].to_string(),
            time,
            event,
            covariate,
        });
    }

    if out.len() < 4 {
        return Err(String::from("Survival analysis requires at least 4 observations."));
    }
    Ok(out)
}

/// Computes Kaplan–Meier (1958) product-limit survival curves with Greenwood SEs & log-log 95% CIs
/// for every group in `obs`.
pub fn compute_kaplan_meier_curves(obs: &[SurvivalObservation]) -> Vec<KaplanMeierCurve> {
    let mut by_group: BTreeMap<String, Vec<&SurvivalObservation>> = BTreeMap::new();
    for r in obs {
        by_group.entry(r.group.clone()).or_default().push(r);
    }

    let mut curves = Vec::with_capacity(by_group.len());
    for (grp, items) in by_group {
        let total_n = items.len();
        let total_events = items.iter().filter(|r| r.event).count();

        // Group by rounded distinct time points
        let mut by_time: BTreeMap<i64, (f64, usize, usize)> = BTreeMap::new();
        for r in items {
            let key = (r.time * 1_000_000.0).round() as i64;
            let entry = by_time.entry(key).or_insert((r.time, 0usize, 0usize));
            if r.event {
                entry.1 += 1;
            } else {
                entry.2 += 1;
            }
        }

        let mut steps = Vec::with_capacity(by_time.len() + 1);
        steps.push(KaplanMeierStep {
            time: 0.0,
            n_at_risk: total_n,
            n_events: 0,
            n_censored: 0,
            survival_s: 1.0,
            greenwood_se: 0.0,
            ci95: (1.0, 1.0),
        });

        let mut at_risk = total_n;
        let mut s_curr = 1.0_f64;
        let mut gw_sum = 0.0_f64;
        let mut median_t: Option<f64> = None;

        for (_, (t_val, d_i, c_i)) in by_time {
            if at_risk == 0 {
                break;
            }
            if d_i > 0 {
                let surv_step = 1.0 - (d_i as f64) / (at_risk as f64);
                s_curr *= surv_step;
                if at_risk > d_i {
                    gw_sum += (d_i as f64) / ((at_risk * (at_risk - d_i)) as f64);
                }
            }
            let se = s_curr * gw_sum.sqrt();
            let ci95 = if s_curr > 0.0 && s_curr < 1.0 && gw_sum > 0.0 {
                // Kalbfleisch–Prentice log-log CI bounded strictly inside [0, 1]
                let se_loglog = gw_sum.sqrt() / s_curr.ln().abs().max(1e-9);
                let exp_lo = (1.96 * se_loglog).exp();
                let exp_hi = (-1.96 * se_loglog).exp();
                (s_curr.powf(exp_lo).clamp(0.0, 1.0), s_curr.powf(exp_hi).clamp(0.0, 1.0))
            } else {
                (s_curr, s_curr)
            };

            if median_t.is_none() && s_curr <= 0.5 {
                median_t = Some(t_val);
            }

            steps.push(KaplanMeierStep {
                time: t_val,
                n_at_risk: at_risk,
                n_events: d_i,
                n_censored: c_i,
                survival_s: s_curr,
                greenwood_se: se,
                ci95,
            });
            at_risk = at_risk.saturating_sub(d_i + c_i);
        }

        curves.push(KaplanMeierCurve {
            group: grp,
            total_n,
            total_events,
            median_survival_time: median_t,
            steps,
        });
    }
    curves
}

/// Performs Mantel's (1966) two-group log-rank test comparing survival between the first two groups in `obs`.
pub fn compute_log_rank_test(obs: &[SurvivalObservation]) -> Result<LogRankTestResult, String> {
    let mut groups: Vec<String> = Vec::new();
    for r in obs {
        if !groups.contains(&r.group) {
            groups.push(r.group.clone());
        }
    }
    if groups.len() < 2 {
        return Err(String::from("Log-rank test requires at least 2 distinct groups."));
    }
    let ga = &groups[0];
    let gb = &groups[1];

    let mut by_time: BTreeMap<i64, (f64, usize, usize, usize, usize)> = BTreeMap::new();
    let mut n1_rem = 0usize;
    let mut n2_rem = 0usize;

    for r in obs {
        if r.group != *ga && r.group != *gb {
            continue;
        }
        let key = (r.time * 1_000_000.0).round() as i64;
        let e = by_time.entry(key).or_insert((r.time, 0, 0, 0, 0));
        if r.group == *ga {
            n1_rem += 1;
            if r.event {
                e.1 += 1;
            } else {
                e.2 += 1;
            }
        } else {
            n2_rem += 1;
            if r.event {
                e.3 += 1;
            } else {
                e.4 += 1;
            }
        }
    }

    let mut obs_a = 0.0_f64;
    let mut exp_a = 0.0_f64;
    let mut obs_b = 0.0_f64;
    let mut exp_b = 0.0_f64;
    let mut var_sum = 0.0_f64;

    for (_, (_, d1, c1, d2, c2)) in by_time {
        let d_tot = d1 + d2;
        let n_tot = n1_rem + n2_rem;
        if d_tot > 0 && n_tot > 1 {
            let n1f = n1_rem as f64;
            let n2f = n2_rem as f64;
            let nf = n_tot as f64;
            let df = d_tot as f64;

            let e1 = n1f * df / nf;
            let e2 = n2f * df / nf;
            obs_a += d1 as f64;
            obs_b += d2 as f64;
            exp_a += e1;
            exp_b += e2;

            let v = (n1f * n2f * df * (nf - df)) / (nf * nf * (nf - 1.0));
            var_sum += v;
        }
        n1_rem = n1_rem.saturating_sub(d1 + c1);
        n2_rem = n2_rem.saturating_sub(d2 + c2);
    }

    let chi2 = if var_sum > 1e-12 {
        (obs_a - exp_a).powi(2) / var_sum
    } else {
        0.0
    };
    let p_val = chi2_df1_pvalue(chi2);

    Ok(LogRankTestResult {
        group_a: ga.clone(),
        group_b: gb.clone(),
        observed_a: obs_a,
        expected_a: exp_a,
        observed_b: obs_b,
        expected_b: exp_b,
        chi_square: chi2,
        p_value: p_val,
    })
}

/// Fits a Cox (1972) proportional-hazards regression model on `covariate` (or group indicator if `use_group_as_covariate`)
/// via exact Newton–Raphson maximization of the Breslow partial likelihood.
pub fn fit_cox_proportional_hazards(
    obs: &[SurvivalObservation],
    covariate_label: &str,
    use_group_as_covariate: bool,
) -> Result<CoxPhFit, String> {
    if obs.len() < 4 {
        return Err(String::from("Cox PH regression requires >= 4 observations."));
    }
    let first_grp = &obs[0].group;
    let xs: Vec<f64> = obs
        .iter()
        .map(|r| {
            if use_group_as_covariate {
                if r.group == *first_grp { 0.0 } else { 1.0 }
            } else {
                r.covariate
            }
        })
        .collect();

    let n_events = obs.iter().filter(|r| r.event).count();
    if n_events < 2 {
        return Err(String::from("Cox PH regression requires >= 2 uncensored events."));
    }

    let eval_partial = |beta: f64| -> (f64, f64, f64) {
        let mut ll = 0.0_f64;
        let mut grad = 0.0_f64;
        let mut info = 0.0_f64;
        let exps: Vec<f64> = xs.iter().map(|&x| (beta * x).clamp(-30.0, 30.0).exp()).collect();

        for (i, r) in obs.iter().enumerate() {
            if !r.event {
                continue;
            }
            let mut s0 = 0.0_f64;
            let mut s1 = 0.0_f64;
            let mut s2 = 0.0_f64;
            for (j, rj) in obs.iter().enumerate() {
                if rj.time >= r.time {
                    let w = exps[j];
                    s0 += w;
                    s1 += xs[j] * w;
                    s2 += xs[j] * xs[j] * w;
                }
            }
            if s0 > 1e-15 {
                let mean_x = s1 / s0;
                let var_x = (s2 / s0 - mean_x * mean_x).max(0.0);
                ll += beta * xs[i] - s0.ln();
                grad += xs[i] - mean_x;
                info += var_x;
            }
        }
        (ll, grad, info)
    };

    let (ll_null, _, _) = eval_partial(0.0);
    let mut beta = 0.0_f64;
    for _ in 0..30 {
        let (_, grad, info) = eval_partial(beta);
        if info < 1e-12 {
            break;
        }
        let step = (grad / info).clamp(-2.0, 2.0);
        beta += step;
        if step.abs() < 1e-10 {
            break;
        }
    }

    let (ll_fit, _, info_final) = eval_partial(beta);
    let se_beta = if info_final > 1e-12 {
        (1.0 / info_final).sqrt()
    } else {
        0.0
    };
    let hr = beta.exp();
    let hr_ci95 = (
        (beta - 1.96 * se_beta).exp(),
        (beta + 1.96 * se_beta).exp(),
    );
    let wald_z = if se_beta > 1e-12 { beta / se_beta } else { 0.0 };
    let p_val = chi2_df1_pvalue(wald_z * wald_z);

    Ok(CoxPhFit {
        covariate_name: covariate_label.to_string(),
        n_subjects: obs.len(),
        n_events,
        beta,
        se_beta,
        hazard_ratio: hr,
        hr_ci95,
        wald_z,
        p_value: p_val,
        partial_log_lik_null: ll_null,
        partial_log_lik_fit: ll_fit,
    })
}

/// Exact two-sided `P(Chi^2_1 >= x) = erfc(sqrt(x / 2))` via complementary normal tail.
pub(crate) fn chi2_df1_pvalue(chi2: f64) -> f64 {
    if chi2 <= 0.0 || !chi2.is_finite() {
        return 1.0;
    }
    let z = chi2.sqrt();
    let t = 1.0 / (1.0 + 0.231_641_9 * z);
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let pdf = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
    (2.0 * pdf * poly).clamp(0.0, 1.0)
}

/// Built-in 40-fish Pacific Northwest juvenile Chinook salmon (*Oncorhynchus tshawytscha*)
/// 60-day smolt survival cohort (`Wild_Stream` vs `Hatchery_Control` with fork-length covariate).
pub fn demo_survival_cohort_csv() -> &'static str {
    "\
# PNW Juvenile Chinook Salmon 60-Day Smolt Survival Cohort (n=40, Right-Censored)
subject_id,group,time_days,event,fork_len_mm
CHN_01,Wild_Stream,60.0,0,104.2
CHN_02,Wild_Stream,60.0,0,108.5
CHN_03,Wild_Stream,48.0,1,92.1
CHN_04,Wild_Stream,60.0,0,111.0
CHN_05,Wild_Stream,54.0,0,101.4
CHN_06,Wild_Stream,60.0,0,106.8
CHN_07,Wild_Stream,39.0,1,89.4
CHN_08,Wild_Stream,60.0,0,112.5
CHN_09,Wild_Stream,60.0,0,105.0
CHN_10,Wild_Stream,52.0,1,95.2
CHN_11,Wild_Stream,60.0,0,109.1
CHN_12,Wild_Stream,60.0,0,103.6
CHN_13,Wild_Stream,44.0,0,99.8
CHN_14,Wild_Stream,60.0,0,110.3
CHN_15,Wild_Stream,60.0,0,107.4
CHN_16,Wild_Stream,31.0,1,88.0
CHN_17,Wild_Stream,60.0,0,106.2
CHN_18,Wild_Stream,60.0,0,114.0
CHN_19,Wild_Stream,57.0,0,102.9
CHN_20,Wild_Stream,60.0,0,108.0
CHN_21,Hatchery_Control,14.0,1,82.5
CHN_22,Hatchery_Control,18.0,1,84.1
CHN_23,Hatchery_Control,22.0,1,86.0
CHN_24,Hatchery_Control,25.0,1,87.2
CHN_25,Hatchery_Control,29.0,1,89.5
CHN_26,Hatchery_Control,33.0,1,91.0
CHN_27,Hatchery_Control,35.0,0,93.2
CHN_28,Hatchery_Control,38.0,1,90.4
CHN_29,Hatchery_Control,41.0,1,92.8
CHN_30,Hatchery_Control,45.0,1,94.0
CHN_31,Hatchery_Control,47.0,1,93.5
CHN_32,Hatchery_Control,51.0,1,95.1
CHN_33,Hatchery_Control,55.0,0,98.6
CHN_34,Hatchery_Control,58.0,1,96.4
CHN_35,Hatchery_Control,60.0,0,101.0
CHN_36,Hatchery_Control,16.0,1,83.0
CHN_37,Hatchery_Control,27.0,1,88.8
CHN_38,Hatchery_Control,36.0,1,90.9
CHN_39,Hatchery_Control,60.0,0,102.3
CHN_40,Hatchery_Control,49.0,1,94.7
"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_kaplan_meier_log_rank_and_cox_ph() {
        let obs = parse_survival_csv(demo_survival_cohort_csv()).unwrap();
        let curves = compute_kaplan_meier_curves(&obs);
        assert_eq!(curves.len(), 2);

        let lr = compute_log_rank_test(&obs).unwrap();
        assert!(lr.chi_square > 10.0);
        assert!(lr.p_value < 0.002);

        // Larger fork length reduces hazard (negative beta => HR < 1.0)
        let cox = fit_cox_proportional_hazards(&obs, "fork_len_mm", false).unwrap();
        assert!(cox.beta < -0.10);
        assert!(cox.hazard_ratio < 0.90);
        assert!(cox.p_value < 0.001);
    }
}
