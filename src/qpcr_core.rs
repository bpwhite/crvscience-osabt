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

//! Tool #12 (`osabt qpcr`) — Quantitative Real-Time PCR (`RT-qPCR`) analysis engine:
//! - **MIQE Compliance Auditing** — Bustin et al. (2009) *Clin. Chem.* 55(4):611–622
//! - **Efficiency-Corrected Ratio** — Pfaffl (2001) *Nucleic Acids Res.* 29(9):e45
//! - **Comparative `2^-ddCq`** — Livak & Schmittgen (2001) *Methods* 25(4):402–408
//! - **`geNorm` Reference Stability (`M`)** — Vandesompele et al. (2002) *Genome Biol.* 3(7):research0034
//! - **`qBase` Multi-Reference NRQ** — Hellemans et al. (2007) *Genome Biol.* 8(2):R19

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use crate::svg_util::xml_escape;

/// Functional role of a qPCR reaction well on a plate.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum QpcrWellRole {
    Unknown,
    Standard,
    NoTemplateControl,
    NoRtControl,
}

impl QpcrWellRole {
    pub fn from_token(raw: &str) -> Self {
        match raw.trim().to_ascii_lowercase().as_str() {
            "std" | "standard" | "cal" => Self::Standard,
            "ntc" | "no_template" | "water" | "blank" => Self::NoTemplateControl,
            "nort" | "no_rt" | "-rt" | "minus_rt" => Self::NoRtControl,
            _ => Self::Unknown,
        }
    }
}

/// Single parsed qPCR reaction well record.
#[derive(Debug, Clone)]
pub struct QpcrReactionWell {
    pub well_id: String,
    pub sample_id: String,
    pub condition_group: String,
    pub target_gene: String,
    pub role: QpcrWellRole,
    /// Measured quantification cycle (`None` = Undetermined / No amplification within 40 cycles).
    pub cq: Option<f64>,
    /// Starting template copy number or relative dilution quantity (used when `role == Standard`).
    pub template_quantity: Option<f64>,
}

/// Linear regression of `Cq = slope * log10(Quantity) + intercept` for a serial dilution series.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StandardCurveFit {
    pub target_gene: String,
    pub n_points: usize,
    pub slope: f64,
    pub intercept_cq: f64,
    pub r_squared: f64,
    /// Per-cycle amplification factor `E = 10^(-1 / slope)` (`2.0` = 100% doubling).
    pub amplification_factor_e: f64,
    /// Percentage efficiency `(E - 1.0) * 100.0`.
    pub efficiency_pct: f64,
    /// True if `efficiency_pct in [90.0, 110.0]` and `r_squared >= 0.980` per MIQE guidelines.
    pub miqe_pass: bool,
}

/// Vandesompele et al. (2002) `geNorm` pairwise log-ratio reference-gene stability metric `M`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GeNormGeneStability {
    pub gene_name: String,
    /// Mean pairwise standard deviation `M_j` across all candidate internal controls (`M < 0.5` homogeneous).
    pub stability_m: f64,
    pub stable_pass: bool,
}

/// Biological sample-level normalized relative quantity (`NRQ`) for one gene of interest (`GOI`).
#[derive(Debug, Clone)]
pub struct SampleNormalizedExpression {
    pub sample_id: String,
    pub condition_group: String,
    pub target_gene: String,
    pub mean_cq: f64,
    pub tech_sd_cq: f64,
    pub n_tech_reps: usize,
    /// Hellemans et al. (2007) `qBase` multi-reference efficiency-corrected NRQ (calibrator geo-mean = 1.0).
    pub nrq_qbase: f64,
    pub log2_nrq: f64,
    /// Pfaffl (2001) ratio normalized solely against the first reference gene.
    pub pfaffl_single_ref: f64,
    /// Livak & Schmittgen (2001) `2^-ddCq` assuming `E = 2.0` against the geometric mean of reference Cqs.
    pub livak_2ddcq: f64,
}

/// Condition-group differential expression comparison against `calibrator_group`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GroupDifferentialExpression {
    pub target_gene: String,
    pub treatment_group: String,
    pub calibrator_group: String,
    pub n_treatment_bio: usize,
    pub n_calibrator_bio: usize,
    pub fold_change_qbase: f64,
    pub fold_change_ci95: (f64, f64),
    pub log2_fold_change: f64,
    pub log2_fc_se: f64,
    pub pfaffl_ratio: f64,
    pub livak_ratio: f64,
    pub welch_t_stat: f64,
    pub p_value_two_sided: f64,
}

/// Diagnostic status level for Bustin et al. (2009) MIQE checklist auditing.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MiqeFlagSeverity {
    Pass,
    Warning,
    Fail,
}

impl MiqeFlagSeverity {
    pub fn label(self) -> &'static str {
        match self {
            Self::Pass => "PASS",
            Self::Warning => "WARN",
            Self::Fail => "FAIL",
        }
    }
}

/// Single automated MIQE quality-control check item.
#[derive(Debug, Clone)]
pub struct MiqeComplianceCheck {
    pub category: &'static str,
    pub subject: String,
    pub severity: MiqeFlagSeverity,
    pub detail: String,
}

/// Complete RT-qPCR quantification, reference stability, and MIQE compliance report.
#[derive(Debug, Clone)]
pub struct QpcrDatasetAnalysis {
    pub calibrator_group: String,
    pub reference_genes_used: Vec<String>,
    pub standard_curves: Vec<StandardCurveFit>,
    pub genorm_rankings: Vec<GeNormGeneStability>,
    pub sample_nrqs: Vec<SampleNormalizedExpression>,
    pub group_comparisons: Vec<GroupDifferentialExpression>,
    pub miqe_checks: Vec<MiqeComplianceCheck>,
}

/// Parses a comma- or tab-delimited qPCR reaction table with columns:
/// `well,sample,group,target,role,cq,quantity` (`cq` may be `"Undetermined"` / `"NA"` / `""`).
pub fn parse_qpcr_table(raw: &str) -> Result<Vec<QpcrReactionWell>, String> {
    let mut wells = Vec::new();
    let mut header_seen = false;

    for (line_idx, line) in raw.lines().enumerate() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with('#') {
            continue;
        }
        let sep = if trimmed.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = trimmed.split(sep).map(|s| s.trim()).collect();
        if !header_seen {
            header_seen = true;
            if cols[0].eq_ignore_ascii_case("well") || cols[0].eq_ignore_ascii_case("well_id") {
                continue;
            }
        }
        if cols.len() < 6 {
            return Err(format!(
                "Line {}: expected >= 6 columns (well,sample,group,target,role,cq[,quantity]), got {}.",
                line_idx + 1,
                cols.len()
            ));
        }

        let cq_tok = cols[5];
        let cq = if cq_tok.is_empty()
            || cq_tok.eq_ignore_ascii_case("undetermined")
            || cq_tok.eq_ignore_ascii_case("na")
            || cq_tok.eq_ignore_ascii_case("nan")
        {
            None
        } else {
            Some(
                cq_tok
                    .parse::<f64>()
                    .map_err(|_| format!("Line {}: invalid Cq value '{cq_tok}'.", line_idx + 1))?,
            )
        };

        let qty = if cols.len() >= 7 && !cols[6].is_empty() && !cols[6].eq_ignore_ascii_case("na") {
            cols[6].parse::<f64>().ok()
        } else {
            None
        };

        wells.push(QpcrReactionWell {
            well_id: cols[0].to_string(),
            sample_id: cols[1].to_string(),
            condition_group: cols[2].to_string(),
            target_gene: cols[3].to_string(),
            role: QpcrWellRole::from_token(cols[4]),
            cq,
            template_quantity: qty,
        });
    }

    if wells.is_empty() {
        return Err(String::from("No qPCR reaction rows found in input table."));
    }
    Ok(wells)
}

/// Fits log10 standard curves (`Cq = slope * log10(quantity) + intercept`) for every target with `>= 3` `Standard` wells.
pub fn fit_qpcr_standard_curves(wells: &[QpcrReactionWell]) -> Vec<StandardCurveFit> {
    let mut by_gene: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
    for w in wells {
        if w.role == QpcrWellRole::Standard {
            if let (Some(cq), Some(qty)) = (w.cq, w.template_quantity) {
                if qty > 0.0 && cq.is_finite() {
                    by_gene
                        .entry(w.target_gene.clone())
                        .or_default()
                        .push((qty.log10(), cq));
                }
            }
        }
    }

    let mut out = Vec::new();
    for (gene, pts) in by_gene {
        if pts.len() < 3 {
            continue;
        }
        let n = pts.len() as f64;
        let mx = pts.iter().map(|p| p.0).sum::<f64>() / n;
        let my = pts.iter().map(|p| p.1).sum::<f64>() / n;
        let mut sxx = 0.0_f64;
        let mut sxy = 0.0_f64;
        let mut syy = 0.0_f64;
        for &(x, y) in &pts {
            sxx += (x - mx) * (x - mx);
            sxy += (x - mx) * (y - my);
            syy += (y - my) * (y - my);
        }
        if sxx <= 1e-12 || sxy >= -1e-9 {
            continue;
        }
        let slope = sxy / sxx;
        let intercept = my - slope * mx;
        let r2 = if syy > 1e-12 {
            ((sxy * sxy) / (sxx * syy)).clamp(0.0, 1.0)
        } else {
            1.0
        };
        let e_fac = 10_f64.powf(-1.0 / slope);
        let eff_pct = (e_fac - 1.0) * 100.0;
        let pass = (90.0..=110.0).contains(&eff_pct) && r2 >= 0.980;

        out.push(StandardCurveFit {
            target_gene: gene,
            n_points: pts.len(),
            slope,
            intercept_cq: intercept,
            r_squared: r2,
            amplification_factor_e: e_fac,
            efficiency_pct: eff_pct,
            miqe_pass: pass,
        });
    }
    out
}

/// Computes Vandesompele et al. (2002) `geNorm` pairwise expression stability `M` across `candidate_ref_genes`.
pub fn compute_genorm_reference_stability(
    wells: &[QpcrReactionWell],
    candidate_ref_genes: &[String],
    efficiencies: &BTreeMap<String, f64>,
) -> Vec<GeNormGeneStability> {
    if candidate_ref_genes.len() < 2 {
        return candidate_ref_genes
            .iter()
            .map(|g| GeNormGeneStability {
                gene_name: g.clone(),
                stability_m: 0.0,
                stable_pass: true,
            })
            .collect();
    }

    // Aggregate sample mean Cq per (sample_id, gene)
    let mut sample_gene_cqs: BTreeMap<String, BTreeMap<String, Vec<f64>>> = BTreeMap::new();
    for w in wells {
        if w.role == QpcrWellRole::Unknown {
            if let Some(cq) = w.cq {
                if candidate_ref_genes.contains(&w.target_gene) {
                    sample_gene_cqs
                        .entry(w.sample_id.clone())
                        .or_default()
                        .entry(w.target_gene.clone())
                        .or_default()
                        .push(cq);
                }
            }
        }
    }

    // Retain samples that have measurements for all candidate reference genes
    let complete_samples: Vec<BTreeMap<String, f64>> = sample_gene_cqs
        .into_values()
        .filter_map(|gmap| {
            if candidate_ref_genes.iter().all(|g| gmap.contains_key(g)) {
                let mut avg = BTreeMap::new();
                for g in candidate_ref_genes {
                    let v = &gmap[g];
                    avg.insert(g.clone(), v.iter().sum::<f64>() / (v.len() as f64));
                }
                Some(avg)
            } else {
                None
            }
        })
        .collect();

    if complete_samples.len() < 2 {
        return Vec::new();
    }

    // Convert Cq to relative quantities `Q_ij = E_j^(min_Cq_j - Cq_ij)`
    let mut rel_q: BTreeMap<String, Vec<f64>> = BTreeMap::new();
    for g in candidate_ref_genes {
        let e_g = efficiencies.get(g).copied().unwrap_or(2.0);
        let min_cq = complete_samples
            .iter()
            .map(|m| m[g])
            .fold(f64::INFINITY, f64::min);
        let qs: Vec<f64> = complete_samples
            .iter()
            .map(|m| e_g.powf(min_cq - m[g]))
            .collect();
        rel_q.insert(g.clone(), qs);
    }

    let mut active: Vec<String> = candidate_ref_genes.to_vec();
    let mut results = Vec::with_capacity(active.len());

    while active.len() > 2 {
        let m_now: Vec<(usize, String, f64)> = active
            .iter()
            .enumerate()
            .map(|(j, gj)| {
                let qj = &rel_q[gj];
                let mut v_sum = 0.0_f64;
                for (k, gk) in active.iter().enumerate() {
                    if j == k {
                        continue;
                    }
                    let qk = &rel_q[gk];
                    let log_ratios: Vec<f64> = qj
                        .iter()
                        .zip(qk)
                        .map(|(&a, &b)| (a / b.max(1e-30)).log2())
                        .collect();
                    v_sum += sample_sd(&log_ratios);
                }
                (j, gj.clone(), v_sum / ((active.len() - 1) as f64))
            })
            .collect();

        let (worst_idx, worst_name, worst_m) = m_now
            .into_iter()
            .max_by(|a, b| a.2.partial_cmp(&b.2).unwrap())
            .unwrap();
        results.push(GeNormGeneStability {
            gene_name: worst_name,
            stability_m: worst_m,
            stable_pass: worst_m <= 0.50,
        });
        active.remove(worst_idx);
    }

    if active.len() == 2 {
        let q0 = &rel_q[&active[0]];
        let q1 = &rel_q[&active[1]];
        let log_ratios: Vec<f64> = q0
            .iter()
            .zip(q1)
            .map(|(&a, &b)| (a / b.max(1e-30)).log2())
            .collect();
        let best_pair_m = sample_sd(&log_ratios);
        for g in active {
            results.push(GeNormGeneStability {
                gene_name: g,
                stability_m: best_pair_m,
                stable_pass: best_pair_m <= 0.50,
            });
        }
    }

    results.sort_by(|a, b| a.stability_m.partial_cmp(&b.stability_m).unwrap());
    results
}

/// Runs end-to-end RT-qPCR quantification (`qBase` multi-ref, `Pfaffl` 2001, `Livak` `2^-ddCq`),
/// `geNorm` reference stability ranking, and Bustin et al. (2009) MIQE compliance auditing.
pub fn analyze_qpcr_experiment(
    wells: &[QpcrReactionWell],
    reference_genes: &[String],
    calibrator_group: &str,
) -> Result<QpcrDatasetAnalysis, String> {
    if reference_genes.is_empty() {
        return Err(String::from(
            "Specify at least 1 internal reference gene (MIQE recommends >= 2).",
        ));
    }

    let std_curves = fit_qpcr_standard_curves(wells);
    let mut eff_map: BTreeMap<String, f64> = BTreeMap::new();
    for sc in &std_curves {
        eff_map.insert(sc.target_gene.clone(), sc.amplification_factor_e);
    }

    // Also include any gene named with common housekeeping hints in geNorm candidate pool
    let mut all_genes: BTreeSet<String> = BTreeSet::new();
    for w in wells {
        if w.role == QpcrWellRole::Unknown {
            all_genes.insert(w.target_gene.clone());
        }
    }
    let mut genorm_pool: Vec<String> = reference_genes.to_vec();
    for g in &all_genes {
        let up = g.to_ascii_uppercase();
        if !genorm_pool.contains(g)
            && matches!(
                up.as_str(),
                "ACTB" | "GAPDH" | "B2M" | "RPL13A" | "HPRT1" | "TBP" | "18S"
            )
        {
            genorm_pool.push(g.clone());
        }
    }
    let genorm_rankings = compute_genorm_reference_stability(wells, &genorm_pool, &eff_map);

    // Group technical replicates per (sample_id, condition_group, target_gene)
    let mut sample_groups: BTreeMap<String, String> = BTreeMap::new();
    let mut tech_map: BTreeMap<(String, String), Vec<f64>> = BTreeMap::new();
    for w in wells {
        if w.role == QpcrWellRole::Unknown {
            sample_groups.insert(w.sample_id.clone(), w.condition_group.clone());
            if let Some(cq) = w.cq {
                tech_map
                    .entry((w.sample_id.clone(), w.target_gene.clone()))
                    .or_default()
                    .push(cq);
            }
        }
    }

    let cal_samples: Vec<String> = sample_groups
        .iter()
        .filter(|(_, grp)| grp.eq_ignore_ascii_case(calibrator_group))
        .map(|(s, _)| s.clone())
        .collect();
    if cal_samples.is_empty() {
        return Err(format!(
            "No Unknown biological samples matched calibrator group '{calibrator_group}'."
        ));
    }

    // Calibrator group mean Cq per gene
    let mut cal_mean_cq: BTreeMap<String, f64> = BTreeMap::new();
    for g in &all_genes {
        let vals: Vec<f64> = cal_samples
            .iter()
            .filter_map(|s| tech_map.get(&(s.clone(), g.clone())))
            .map(|v| mean_only(v))
            .collect();
        if !vals.is_empty() {
            cal_mean_cq.insert(g.clone(), mean_only(&vals));
        }
    }

    let goi_genes: Vec<String> = all_genes
        .iter()
        .filter(|g| !reference_genes.contains(g) && g.to_ascii_uppercase() != "B2M")
        .cloned()
        .collect();

    let primary_ref = &reference_genes[0];
    let mut sample_nrqs = Vec::new();

    for (sample_id, group_name) in &sample_groups {
        // Compute Hellemans (2007) reference normalization factor NF_s and Livak mean ref Cq
        let mut log_nf_sum = 0.0_f64;
        let mut ref_cq_sum = 0.0_f64;
        let mut cal_ref_cq_sum = 0.0_f64;
        let mut valid_refs = 0usize;

        for rg in reference_genes {
            if let (Some(reps), Some(&c_cal)) =
                (tech_map.get(&(sample_id.clone(), rg.clone())), cal_mean_cq.get(rg))
            {
                let c_s = mean_only(reps);
                let e_r = eff_map.get(rg).copied().unwrap_or(2.0);
                log_nf_sum += (c_cal - c_s) * e_r.ln();
                ref_cq_sum += c_s;
                cal_ref_cq_sum += c_cal;
                valid_refs += 1;
            }
        }
        if valid_refs == 0 {
            continue;
        }
        let nf_s = (log_nf_sum / (valid_refs as f64)).exp();
        let mean_ref_cq_s = ref_cq_sum / (valid_refs as f64);
        let mean_ref_cq_cal = cal_ref_cq_sum / (valid_refs as f64);

        // Primary-ref delta for Pfaffl (2001)
        let pfaffl_denom = if let (Some(r_reps), Some(&r_cal)) = (
            tech_map.get(&(sample_id.clone(), primary_ref.clone())),
            cal_mean_cq.get(primary_ref),
        ) {
            let e_r1 = eff_map.get(primary_ref).copied().unwrap_or(2.0);
            e_r1.powf(r_cal - mean_only(r_reps))
        } else {
            nf_s
        };

        for goi in &goi_genes {
            let Some(t_reps) = tech_map.get(&(sample_id.clone(), goi.clone())) else {
                continue;
            };
            let Some(&t_cal) = cal_mean_cq.get(goi) else {
                continue;
            };
            let c_goi = mean_only(t_reps);
            let sd_goi = sample_sd(t_reps);
            let e_goi = eff_map.get(goi).copied().unwrap_or(2.0);

            let rq_goi = e_goi.powf(t_cal - c_goi);
            let nrq = rq_goi / nf_s.max(1e-30);
            let pfaffl = rq_goi / pfaffl_denom.max(1e-30);
            let dd_cq = (c_goi - mean_ref_cq_s) - (t_cal - mean_ref_cq_cal);
            let livak = 2.0_f64.powf(-dd_cq);

            sample_nrqs.push(SampleNormalizedExpression {
                sample_id: sample_id.clone(),
                condition_group: group_name.clone(),
                target_gene: goi.clone(),
                mean_cq: c_goi,
                tech_sd_cq: sd_goi,
                n_tech_reps: t_reps.len(),
                nrq_qbase: nrq,
                log2_nrq: nrq.log2(),
                pfaffl_single_ref: pfaffl,
                livak_2ddcq: livak,
            });
        }
    }

    // Group-level differential expression (Welch t-test on log2_nrq)
    let distinct_groups: BTreeSet<String> = sample_groups.values().cloned().collect();
    let mut group_comparisons = Vec::new();

    for goi in &goi_genes {
        let cal_pts: Vec<&SampleNormalizedExpression> = sample_nrqs
            .iter()
            .filter(|s| {
                &s.target_gene == goi && s.condition_group.eq_ignore_ascii_case(calibrator_group)
            })
            .collect();
        if cal_pts.is_empty() {
            continue;
        }
        let cal_log2: Vec<f64> = cal_pts.iter().map(|s| s.log2_nrq).collect();
        let cal_pf_l2: Vec<f64> = cal_pts.iter().map(|s| s.pfaffl_single_ref.log2()).collect();
        let cal_lv_l2: Vec<f64> = cal_pts.iter().map(|s| s.livak_2ddcq.log2()).collect();

        for grp in &distinct_groups {
            if grp.eq_ignore_ascii_case(calibrator_group) {
                continue;
            }
            let trt_pts: Vec<&SampleNormalizedExpression> = sample_nrqs
                .iter()
                .filter(|s| &s.target_gene == goi && &s.condition_group == grp)
                .collect();
            if trt_pts.is_empty() {
                continue;
            }
            let trt_log2: Vec<f64> = trt_pts.iter().map(|s| s.log2_nrq).collect();
            let trt_pf_l2: Vec<f64> = trt_pts.iter().map(|s| s.pfaffl_single_ref.log2()).collect();
            let trt_lv_l2: Vec<f64> = trt_pts.iter().map(|s| s.livak_2ddcq.log2()).collect();

            let diff_l2 = mean_only(&trt_log2) - mean_only(&cal_log2);
            let fc_qbase = 2.0_f64.powf(diff_l2);
            let pf_ratio = 2.0_f64.powf(mean_only(&trt_pf_l2) - mean_only(&cal_pf_l2));
            let lv_ratio = 2.0_f64.powf(mean_only(&trt_lv_l2) - mean_only(&cal_lv_l2));

            let (se_diff, df, t_stat, p_val) = welch_t_test(&trt_log2, &cal_log2);
            let t_crit = student_t_crit_95(df);
            let ci_lo = 2.0_f64.powf(diff_l2 - t_crit * se_diff);
            let ci_hi = 2.0_f64.powf(diff_l2 + t_crit * se_diff);

            group_comparisons.push(GroupDifferentialExpression {
                target_gene: goi.clone(),
                treatment_group: grp.clone(),
                calibrator_group: calibrator_group.to_string(),
                n_treatment_bio: trt_log2.len(),
                n_calibrator_bio: cal_log2.len(),
                fold_change_qbase: fc_qbase,
                fold_change_ci95: (ci_lo, ci_hi),
                log2_fold_change: diff_l2,
                log2_fc_se: se_diff,
                pfaffl_ratio: pf_ratio,
                livak_ratio: lv_ratio,
                welch_t_stat: t_stat,
                p_value_two_sided: p_val,
            });
        }
    }

    let miqe_checks = build_miqe_audit(wells, reference_genes, &std_curves, &genorm_rankings, &tech_map);

    Ok(QpcrDatasetAnalysis {
        calibrator_group: calibrator_group.to_string(),
        reference_genes_used: reference_genes.to_vec(),
        standard_curves: std_curves,
        genorm_rankings,
        sample_nrqs,
        group_comparisons,
        miqe_checks,
    })
}

fn build_miqe_audit(
    wells: &[QpcrReactionWell],
    ref_genes: &[String],
    std_curves: &[StandardCurveFit],
    genorm: &[GeNormGeneStability],
    tech_map: &BTreeMap<(String, String), Vec<f64>>,
) -> Vec<MiqeComplianceCheck> {
    let mut checks = Vec::new();

    // 1. Multi-reference normalization check (Bustin et al. 2009 / Vandesompele et al. 2002)
    checks.push(if ref_genes.len() >= 2 {
        MiqeComplianceCheck {
            category: "Normalization",
            subject: ref_genes.join(" + "),
            severity: MiqeFlagSeverity::Pass,
            detail: format!(
                "Normalized against {} validated internal reference genes (MIQE recommends >= 2).",
                ref_genes.len()
            ),
        }
    } else {
        MiqeComplianceCheck {
            category: "Normalization",
            subject: ref_genes.join(""),
            severity: MiqeFlagSeverity::Warning,
            detail: String::from(
                "Single reference gene used; MIQE & Vandesompele (2002) recommend >= 2 validated reference genes.",
            ),
        }
    });

    // 2. Standard curve efficiency & R^2 checks
    for sc in std_curves {
        let sev = if sc.miqe_pass {
            MiqeFlagSeverity::Pass
        } else if (85.0..=115.0).contains(&sc.efficiency_pct) && sc.r_squared >= 0.970 {
            MiqeFlagSeverity::Warning
        } else {
            MiqeFlagSeverity::Fail
        };
        checks.push(MiqeComplianceCheck {
            category: "Assay Efficiency",
            subject: sc.target_gene.clone(),
            severity: sev,
            detail: format!(
                "Slope = {:.3}, Eff = {:.1}% (E = {:.3}), R² = {:.4} [MIQE target: 90–110%, R² >= 0.980].",
                sc.slope, sc.efficiency_pct, sc.amplification_factor_e, sc.r_squared
            ),
        });
    }

    // 3. geNorm reference stability check
    for g in genorm {
        let is_selected = ref_genes.contains(&g.gene_name);
        let sev = if g.stability_m <= 0.50 {
            MiqeFlagSeverity::Pass
        } else if is_selected {
            MiqeFlagSeverity::Fail
        } else {
            MiqeFlagSeverity::Warning
        };
        checks.push(MiqeComplianceCheck {
            category: "geNorm Stability",
            subject: g.gene_name.clone(),
            severity: sev,
            detail: format!(
                "geNorm M = {:.3} ({}) [Target M <= 0.50 for homogeneous samples].",
                g.stability_m,
                if is_selected { "Selected Ref" } else { "Excluded Candidate" }
            ),
        });
    }

    // 4. Technical replicate precision (MIQE: Cq SD <= 0.50, ideally <= 0.20)
    let mut max_tech_sd = 0.0_f64;
    let mut worst_pair = String::from("All wells");
    for ((samp, gene), reps) in tech_map {
        let sd = sample_sd(reps);
        if sd > max_tech_sd {
            max_tech_sd = sd;
            worst_pair = format!("{samp}:{gene}");
        }
    }
    checks.push(MiqeComplianceCheck {
        category: "Replicate Precision",
        subject: worst_pair,
        severity: if max_tech_sd <= 0.30 {
            MiqeFlagSeverity::Pass
        } else if max_tech_sd <= 0.50 {
            MiqeFlagSeverity::Warning
        } else {
            MiqeFlagSeverity::Fail
        },
        detail: format!("Maximum technical replicate Cq SD = {max_tech_sd:.3} cycles [MIQE <= 0.50]."),
    });

    // 5. NTC & No-RT contamination checks (Delta Cq >= 5.0 cycles beyond highest sample Cq)
    for role_check in [QpcrWellRole::NoTemplateControl, QpcrWellRole::NoRtControl] {
        let label = if role_check == QpcrWellRole::NoTemplateControl {
            "NTC Contamination"
        } else {
            "gDNA / No-RT Control"
        };
        for w in wells.iter().filter(|w| w.role == role_check) {
            let max_unk_cq = wells
                .iter()
                .filter(|u| u.role == QpcrWellRole::Unknown && u.target_gene == w.target_gene)
                .filter_map(|u| u.cq)
                .fold(0.0_f64, f64::max);
            match w.cq {
                None => checks.push(MiqeComplianceCheck {
                    category: label,
                    subject: format!("{} ({})", w.target_gene, w.well_id),
                    severity: MiqeFlagSeverity::Pass,
                    detail: String::from("Undetermined (> 40 cycles) — zero background amplification."),
                }),
                Some(ntc_cq) => {
                    let margin = ntc_cq - max_unk_cq;
                    let sev = if margin >= 8.0 || ntc_cq >= 37.0 {
                        MiqeFlagSeverity::Warning
                    } else {
                        MiqeFlagSeverity::Fail
                    };
                    checks.push(MiqeComplianceCheck {
                        category: label,
                        subject: format!("{} ({})", w.target_gene, w.well_id),
                        severity: sev,
                        detail: format!(
                            "Late signal at Cq = {ntc_cq:.2} (ΔCq = {margin:+.1} vs slowest sample {max_unk_cq:.2}); inspect SYBR melt curve for primer-dimer."
                        ),
                    });
                }
            }
        }
    }

    checks
}

/// Renders a side-by-side `.svg` figure comparing `qBase` multi-reference, `Pfaffl (2001)`,
/// and `Livak (2001) 2^-ddCq` log2 fold-changes with 95% confidence intervals.
pub fn render_qpcr_expression_svg(analysis: &QpcrDatasetAnalysis, title: &str) -> String {
    let n_comp = analysis.group_comparisons.len().max(1);
    let width = 840.0_f64;
    let height = 440.0_f64;
    let left = 84.0_f64;
    let top = 62.0_f64;
    let pw = 540.0_f64;
    let ph = 300.0_f64;

    let mut y_max = 2.0_f64;
    let mut y_min = -1.0_f64;
    for c in &analysis.group_comparisons {
        y_max = y_max.max(c.fold_change_ci95.1.log2() + 0.6);
        y_min = y_min.min(c.fold_change_ci95.0.log2() - 0.6);
    }
    let y_span = (y_max - y_min).max(1.0);
    let map_y = |l2: f64| top + ph - ((l2 - y_min) / y_span).clamp(0.0, 1.0) * ph;
    let y0 = map_y(0.0);

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
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Relative Expression (log2 Fold Change vs. {})</text>",
        xml_escape(title),
        xml_escape(&analysis.calibrator_group)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{y0:.1}\" x2=\"{:.1}\" y2=\"{y0:.1}\" stroke=\"#475569\" stroke-width=\"1.4\"/>",
        left + pw
    );

    let group_w = pw / (n_comp as f64);
    let bar_w = (group_w * 0.22).min(38.0);

    for (idx, comp) in analysis.group_comparisons.iter().enumerate() {
        let cx = left + ((idx as f64) + 0.5) * group_w;
        let methods = [
            (comp.log2_fold_change, "#2563eb"),
            (comp.pfaffl_ratio.log2(), "#0d9488"),
            (comp.livak_ratio.log2(), "#9333ea"),
        ];
        for (m_i, &(val_l2, col)) in methods.iter().enumerate() {
            let bx = cx + ((m_i as f64) - 1.5) * (bar_w + 4.0);
            let py = map_y(val_l2);
            let bh = (y0 - py).abs().max(1.5);
            let by = py.min(y0);
            let _ = writeln!(
                svg,
                "  <rect x=\"{bx:.1}\" y=\"{by:.1}\" width=\"{bar_w:.1}\" height=\"{bh:.1}\" fill=\"{col}\" rx=\"2\"/>"
            );
            if m_i == 0 {
                let err_x = bx + 0.5 * bar_w;
                let e_hi = map_y(comp.fold_change_ci95.1.log2());
                let e_lo = map_y(comp.fold_change_ci95.0.log2());
                let _ = writeln!(
                    svg,
                    "  <line x1=\"{err_x:.1}\" y1=\"{e_hi:.1}\" x2=\"{err_x:.1}\" y2=\"{e_lo:.1}\" stroke=\"#0f172a\" stroke-width=\"1.6\"/>"
                );
            }
        }
        let _ = writeln!(
            svg,
            "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11.5\" font-weight=\"700\" fill=\"#0f172a\">{} ({})</text>",
            top + ph + 22.0,
            xml_escape(&comp.target_gene),
            xml_escape(&comp.treatment_group)
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">FC={:.2}x (p={:.2e})</text>",
            top + ph + 38.0,
            comp.fold_change_qbase,
            comp.p_value_two_sided
        );
    }

    let lx = left + pw + 18.0;
    for (li, (lbl, col)) in [
        ("qBase Multi-Ref (2007)", "#2563eb"),
        ("Pfaffl Single-Ref (2001)", "#0d9488"),
        ("Livak 2^-ddCq (2001)", "#9333ea"),
    ]
    .into_iter()
    .enumerate()
    {
        let ly = top + 26.0 + (li as f64) * 28.0;
        let _ = writeln!(
            svg,
            "  <rect x=\"{lx:.1}\" y=\"{:.1}\" width=\"12\" height=\"12\" fill=\"{col}\"/>",
            ly - 10.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{ly:.1}\" font-family=\"sans-serif\" font-size=\"11.5\" fill=\"#1e293b\">{lbl}</text>",
            lx + 18.0
        );
    }

    svg.push_str("</svg>\n");
    svg
}

/// Built-in SYBR Green RT-qPCR dataset (`Control` vs. `HeatShock_42C`, 3 bio reps × 2 tech reps,
/// 5-point 10-fold serial dilution standard curves for `HSP70`, `IL6`, `ACTB`, `GAPDH`, plus unstable `B2M`
/// and `NTC`/`NoRT` controls).
pub fn demo_qpcr_heatshock_plate_csv() -> &'static str {
    "\
well,sample,group,target,role,cq,quantity
A01,Std1,Std,HSP70,Standard,16.12,100000
A02,Std2,Std,HSP70,Standard,19.58,10000
A03,Std3,Std,HSP70,Standard,23.05,1000
A04,Std4,Std,HSP70,Standard,26.51,100
A05,Std5,Std,HSP70,Standard,29.99,10
A06,Std1,Std,IL6,Standard,17.40,100000
A07,Std2,Std,IL6,Standard,20.82,10000
A08,Std3,Std,IL6,Standard,24.25,1000
A09,Std4,Std,IL6,Standard,27.67,100
A10,Std5,Std,IL6,Standard,31.11,10
B01,Std1,Std,ACTB,Standard,14.80,100000
B02,Std2,Std,ACTB,Standard,18.12,10000
B03,Std3,Std,ACTB,Standard,21.45,1000
B04,Std4,Std,ACTB,Standard,24.77,100
B05,Std5,Std,ACTB,Standard,28.10,10
B06,Std1,Std,GAPDH,Standard,15.20,100000
B07,Std2,Std,GAPDH,Standard,18.56,10000
B08,Std3,Std,GAPDH,Standard,21.91,1000
B09,Std4,Std,GAPDH,Standard,25.27,100
B10,Std5,Std,GAPDH,Standard,28.62,10
C01,Ctrl_1,Control,HSP70,Unknown,25.40,NA
C02,Ctrl_1,Control,HSP70,Unknown,25.52,NA
C03,Ctrl_2,Control,HSP70,Unknown,25.61,NA
C04,Ctrl_2,Control,HSP70,Unknown,25.49,NA
C05,Ctrl_3,Control,HSP70,Unknown,25.33,NA
C06,Ctrl_3,Control,HSP70,Unknown,25.45,NA
C07,HS_1,HeatShock_42C,HSP70,Unknown,21.18,NA
C08,HS_1,HeatShock_42C,HSP70,Unknown,21.26,NA
C09,HS_2,HeatShock_42C,HSP70,Unknown,21.05,NA
C10,HS_2,HeatShock_42C,HSP70,Unknown,21.19,NA
C11,HS_3,HeatShock_42C,HSP70,Unknown,21.31,NA
C12,HS_3,HeatShock_42C,HSP70,Unknown,21.21,NA
D01,Ctrl_1,Control,IL6,Unknown,26.80,NA
D02,Ctrl_1,Control,IL6,Unknown,26.92,NA
D03,Ctrl_2,Control,IL6,Unknown,27.01,NA
D04,Ctrl_2,Control,IL6,Unknown,26.89,NA
D05,Ctrl_3,Control,IL6,Unknown,26.74,NA
D06,Ctrl_3,Control,IL6,Unknown,26.84,NA
D07,HS_1,HeatShock_42C,IL6,Unknown,24.62,NA
D08,HS_1,HeatShock_42C,IL6,Unknown,24.74,NA
D09,HS_2,HeatShock_42C,IL6,Unknown,24.55,NA
D10,HS_2,HeatShock_42C,IL6,Unknown,24.67,NA
D11,HS_3,HeatShock_42C,IL6,Unknown,24.79,NA
D12,HS_3,HeatShock_42C,IL6,Unknown,24.71,NA
E01,Ctrl_1,Control,ACTB,Unknown,18.92,NA
E02,Ctrl_1,Control,ACTB,Unknown,19.00,NA
E03,Ctrl_2,Control,ACTB,Unknown,19.08,NA
E04,Ctrl_2,Control,ACTB,Unknown,19.02,NA
E05,Ctrl_3,Control,ACTB,Unknown,18.85,NA
E06,Ctrl_3,Control,ACTB,Unknown,18.93,NA
E07,HS_1,HeatShock_42C,ACTB,Unknown,19.04,NA
E08,HS_1,HeatShock_42C,ACTB,Unknown,19.12,NA
E09,HS_2,HeatShock_42C,ACTB,Unknown,18.96,NA
E10,HS_2,HeatShock_42C,ACTB,Unknown,19.04,NA
E11,HS_3,HeatShock_42C,ACTB,Unknown,19.14,NA
E12,HS_3,HeatShock_42C,ACTB,Unknown,19.06,NA
F01,Ctrl_1,Control,GAPDH,Unknown,19.44,NA
F02,Ctrl_1,Control,GAPDH,Unknown,19.50,NA
F03,Ctrl_2,Control,GAPDH,Unknown,19.59,NA
F04,Ctrl_2,Control,GAPDH,Unknown,19.53,NA
F05,Ctrl_3,Control,GAPDH,Unknown,19.38,NA
F06,Ctrl_3,Control,GAPDH,Unknown,19.42,NA
F07,HS_1,HeatShock_42C,GAPDH,Unknown,19.55,NA
F08,HS_1,HeatShock_42C,GAPDH,Unknown,19.61,NA
F09,HS_2,HeatShock_42C,GAPDH,Unknown,19.47,NA
F10,HS_2,HeatShock_42C,GAPDH,Unknown,19.53,NA
F11,HS_3,HeatShock_42C,GAPDH,Unknown,19.64,NA
F12,HS_3,HeatShock_42C,GAPDH,Unknown,19.58,NA
G01,Ctrl_1,Control,B2M,Unknown,21.10,NA
G03,Ctrl_2,Control,B2M,Unknown,21.25,NA
G05,Ctrl_3,Control,B2M,Unknown,20.95,NA
G07,HS_1,HeatShock_42C,B2M,Unknown,23.15,NA
G09,HS_2,HeatShock_42C,B2M,Unknown,23.40,NA
G11,HS_3,HeatShock_42C,B2M,Unknown,22.95,NA
H01,NTC_HSP70,NTC,HSP70,NTC,Undetermined,NA
H02,NTC_IL6,NTC,IL6,NTC,37.85,NA
H03,NTC_ACTB,NTC,ACTB,NTC,Undetermined,NA
H04,NoRT_Pool,NoRT,HSP70,NoRT,Undetermined,NA
"
}

fn mean_only(v: &[f64]) -> f64 {
    if v.is_empty() {
        0.0
    } else {
        v.iter().sum::<f64>() / (v.len() as f64)
    }
}

fn sample_sd(v: &[f64]) -> f64 {
    if v.len() < 2 {
        return 0.0;
    }
    let m = mean_only(v);
    let var = v.iter().map(|&x| (x - m).powi(2)).sum::<f64>() / ((v.len() - 1) as f64);
    var.sqrt()
}

fn welch_t_test(a: &[f64], b: &[f64]) -> (f64, f64, f64, f64) {
    if a.len() < 2 || b.len() < 2 {
        return (0.0, 2.0, 0.0, 1.0);
    }
    let na = a.len() as f64;
    let nb = b.len() as f64;
    let va = sample_sd(a).powi(2);
    let vb = sample_sd(b).powi(2);
    let se2 = va / na + vb / nb;
    if se2 <= 1e-18 {
        return (0.0, na + nb - 2.0, 0.0, 1.0);
    }
    let se = se2.sqrt();
    let t = (mean_only(a) - mean_only(b)) / se;
    let df_num = se2 * se2;
    let df_den = (va / na).powi(2) / (na - 1.0) + (vb / nb).powi(2) / (nb - 1.0);
    let df = (df_num / df_den.max(1e-18)).max(1.0);
    let p_two = student_t_two_sided_p(t.abs(), df);
    (se, df, t, p_two)
}

fn student_t_crit_95(df: f64) -> f64 {
    if df < 1.5 {
        12.706
    } else if df < 2.5 {
        4.303
    } else if df < 3.5 {
        3.182
    } else if df < 4.5 {
        2.776
    } else if df < 6.5 {
        2.447
    } else if df < 10.5 {
        2.228
    } else {
        1.96 + 2.4 / df
    }
}

/// Accurate Hill/Cornish–Fisher normal-equivalent two-sided Student's `t` p-value (`df >= 1`).
fn student_t_two_sided_p(abs_t: f64, df: f64) -> f64 {
    if abs_t <= 0.0 {
        return 1.0;
    }
    // Peizer–Pratt / Wallace normal transform `z` for Student's t:
    // z = sqrt(df * ln(1 + t^2 / df) * (1 - 1 / (2 * df)))
    let corr = (1.0 - 0.5 / df.max(1.0)).max(0.5);
    let z = (df * (1.0 + (abs_t * abs_t) / df).ln() * corr).sqrt();
    (2.0 * standard_normal_sf(z)).clamp(1e-16, 1.0)
}

fn standard_normal_sf(z: f64) -> f64 {
    // Abramowitz & Stegun 7.1.26 erfc approximation
    let x = z / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.3275911 * x.abs());
    let poly = t
        * (0.254829592
            + t * (-0.284496736 + t * (1.421413741 + t * (-1.453152027 + t * 1.061405429))));
    let erfc_x = poly * (-x * x).exp();
    if z >= 0.0 {
        0.5 * erfc_x
    } else {
        1.0 - 0.5 * erfc_x
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn quantifies_heatshock_plate_genorm_and_miqe_audit() {
        let wells = parse_qpcr_table(demo_qpcr_heatshock_plate_csv()).unwrap();
        let refs = vec![String::from("ACTB"), String::from("GAPDH")];
        let res = analyze_qpcr_experiment(&wells, &refs, "Control").unwrap();

        assert_eq!(res.standard_curves.len(), 4);
        for sc in &res.standard_curves {
            assert!(sc.miqe_pass);
            assert!(sc.r_squared > 0.999);
        }

        // ACTB and GAPDH are highly stable (M < 0.2), whereas heat-regulated B2M is flagged unstable (M > 0.8)
        let b2m = res
            .genorm_rankings
            .iter()
            .find(|g| g.gene_name == "B2M")
            .unwrap();
        assert!(!b2m.stable_pass && b2m.stability_m > 0.80);

        // HSP70 is strongly induced (~15–18x) in HeatShock_42C with p < 0.001
        let hsp = res
            .group_comparisons
            .iter()
            .find(|c| c.target_gene == "HSP70")
            .unwrap();
        assert!(hsp.fold_change_qbase > 13.0 && hsp.fold_change_qbase < 20.0);
        assert!(hsp.p_value_two_sided < 0.001);
    }
}
