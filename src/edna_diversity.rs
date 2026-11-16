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

//! Tool #14 (`osabt edna`) — Biodiversity & community ecology quantification engine:
//! - **Asymptotic Richness (`Chao1`)** — Chao (1984) *Scand. J. Stat.* 11:265–270
//! - **Analytical Hypergeometric Rarefaction (`E[S_m]`) & `PIE`** — Hurlbert (1971) *Ecology* 52:577–586
//! - **Hill Numbers (`^0D, ^1D, ^2D`) & Good–Turing Coverage** — Chao et al. (2014) *Ecol. Monogr.* 84:45–67
//! - **Pairwise Beta Diversity** — Bray–Curtis abundance dissimilarity & Jaccard incidence dissimilarity.

use std::fmt::Write as _;

use crate::edna_fastq::EdnaPipelineReport;
use crate::svg_util::xml_escape;

/// Alpha-diversity, Hill numbers, and sample completeness metrics for one environmental sample.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SampleAlphaDiversity {
    pub sample_name: String,
    pub total_reads_n: usize,
    /// Hill `q = 0`: Observed ASV / species richness `S_obs`.
    pub observed_richness_q0: usize,
    pub singletons_f1: usize,
    pub doubletons_f2: usize,
    /// Chao (1984) bias-corrected lower-bound asymptotic richness `Chao1`.
    pub chao1_asymptotic: f64,
    pub chao1_ci95: (f64, f64),
    /// Hill `q = 1`: Exponential Shannon entropy `exp(H')` (effective number of common species).
    pub hill_shannon_q1: f64,
    /// Hill `q = 2`: Inverse Simpson concentration `1 / sum(p_i^2)` (effective number of dominant species).
    pub hill_simpson_q2: f64,
    /// Pielou's evenness `J' = H' / ln(S_obs)`.
    pub pielou_evenness_j: f64,
    /// Hurlbert (1971) unbiased Probability of Interspecific Encounter `PIE = (N/(N-1)) * (1 - sum p_i^2)`.
    pub hurlbert_pie: f64,
    /// Chao & Jost (2012) / Good–Turing sample completeness `C_hat` in `[0, 1]`.
    pub sample_coverage: f64,
}

/// Analytical Hurlbert (1971) rarefaction curve `E[S_m]` across subsample read depths `m`.
#[derive(Debug, Clone)]
pub struct RarefactionCurve {
    pub sample_name: String,
    pub depths_m: Vec<usize>,
    pub expected_species_sm: Vec<f64>,
}

/// Full multi-site alpha, rarefaction, and beta-diversity summary excluding chimeras & blank contaminants.
#[derive(Debug, Clone)]
pub struct CommunityDiversityReport {
    pub clean_asv_count: usize,
    pub alpha: Vec<SampleAlphaDiversity>,
    pub rarefaction_curves: Vec<RarefactionCurve>,
    pub sample_names: Vec<String>,
    pub bray_curtis_matrix: Vec<Vec<f64>>,
    pub jaccard_matrix: Vec<Vec<f64>>,
}

/// Computes Hurlbert's (1971) exact analytical hypergeometric expected richness `E[S_m]`
/// in a random subsample of `subsample_m` reads drawn without replacement from `counts`:
/// `E[S_m] = sum_{N_i > 0} (1 - C(N - N_i, m) / C(N, m))`.
pub fn hurlbert_expected_richness(counts: &[usize], subsample_m: usize) -> f64 {
    let total_n: usize = counts.iter().sum();
    let s_obs = counts.iter().filter(|&&c| c > 0).count();
    if subsample_m == 0 || total_n == 0 {
        return 0.0;
    }
    if subsample_m >= total_n {
        return s_obs as f64;
    }

    let mut expected = 0.0_f64;
    for &ni in counts {
        if ni == 0 {
            continue;
        }
        if total_n - ni < subsample_m {
            // Guaranteed to draw at least 1 individual of species `i`
            expected += 1.0;
        } else {
            // Exact log-ratio `ln(C(N - Ni, m) / C(N, m)) = sum_{k=0..m-1} ln((N - Ni - k) / (N - k))`
            let mut log_p_miss = 0.0_f64;
            for k in 0..subsample_m {
                let num = (total_n - ni - k) as f64;
                let den = (total_n - k) as f64;
                log_p_miss += (num / den).ln();
            }
            expected += 1.0 - log_p_miss.exp();
        }
    }
    expected
}

/// Computes Chao's (1984 / 1987) bias-corrected `Chao1` richness estimator and asymmetric log-normal 95% CI.
pub fn compute_chao1_with_ci(counts: &[usize]) -> (f64, (f64, f64)) {
    let s_obs = counts.iter().filter(|&&c| c > 0).count() as f64;
    let f1 = counts.iter().filter(|&&c| c == 1).count() as f64;
    let f2 = counts.iter().filter(|&&c| c == 2).count() as f64;

    if s_obs == 0.0 || f1 == 0.0 {
        return (s_obs, (s_obs, s_obs));
    }

    let extra = if f2 > 0.0 {
        (f1 * f1) / (2.0 * f2)
    } else {
        (f1 * (f1 - 1.0)) / 2.0
    };
    let chao1 = s_obs + extra;

    let var = if f2 > 0.0 {
        f2 * (0.5 * (f1 / f2).powi(2) + (f1 / f2).powi(3) + 0.25 * (f1 / f2).powi(4))
    } else {
        0.5 * f1 * (f1 - 1.0) + 0.25 * f1 * (2.0 * f1 - 1.0).powi(2)
            - 0.25 * f1.powi(4) / chao1.max(1.0)
    }
    .max(0.0);

    if extra <= 1e-9 || var <= 1e-12 {
        return (chao1, (s_obs, chao1));
    }

    let c_fac = (1.96 * (1.0 + var / (extra * extra)).ln().sqrt()).exp();
    let ci_lo = s_obs + extra / c_fac;
    let ci_hi = s_obs + extra * c_fac;
    (chao1, (ci_lo, ci_hi))
}

/// Computes alpha diversity (`^0D`, `Chao1`, `^1D`, `^2D`, `PIE`, coverage), Hurlbert rarefaction curves,
/// and pairwise Bray–Curtis & Jaccard beta-diversity across environmental samples in `pipeline`
/// (automatically excluding chimeras, blank-flagged contaminants, and negative control wells).
pub fn analyze_edna_community_diversity(
    pipeline: &EdnaPipelineReport,
    rarefaction_steps: usize,
) -> CommunityDiversityReport {
    let clean_asvs: Vec<_> = pipeline
        .asvs
        .iter()
        .filter(|a| !a.is_chimera && !a.blank_flagged)
        .collect();

    let env_indices: Vec<usize> = (0..pipeline.sample_names.len())
        .filter(|i| !pipeline.blank_sample_indices.contains(i))
        .collect();

    let env_names: Vec<String> = env_indices
        .iter()
        .map(|&i| pipeline.sample_names[i].clone())
        .collect();

    let n_sites = env_indices.len();
    let mut site_vectors: Vec<Vec<usize>> = Vec::with_capacity(n_sites);
    let mut alpha = Vec::with_capacity(n_sites);
    let mut rarefaction_curves = Vec::with_capacity(n_sites);

    let steps = rarefaction_steps.max(8);

    for (pos_s, &s_idx) in env_indices.iter().enumerate() {
        let counts: Vec<usize> = clean_asvs.iter().map(|a| a.counts_by_sample[s_idx]).collect();
        let n_tot: usize = counts.iter().sum();
        let s_obs = counts.iter().filter(|&&c| c > 0).count();
        let f1 = counts.iter().filter(|&&c| c == 1).count();
        let f2 = counts.iter().filter(|&&c| c == 2).count();
        let (chao1, chao1_ci) = compute_chao1_with_ci(&counts);

        let (q1, q2, pielou, pie, cov) = if n_tot > 0 {
            let nf = n_tot as f64;
            let mut h_shannon = 0.0_f64;
            let mut simpson_sum = 0.0_f64;
            for &c in &counts {
                if c > 0 {
                    let p = (c as f64) / nf;
                    h_shannon -= p * p.ln();
                    simpson_sum += p * p;
                }
            }
            let hill1 = h_shannon.exp();
            let hill2 = if simpson_sum > 0.0 {
                1.0 / simpson_sum
            } else {
                0.0
            };
            let j_even = if s_obs > 1 {
                h_shannon / (s_obs as f64).ln()
            } else {
                1.0
            };
            let h_pie = if n_tot > 1 {
                (nf / (nf - 1.0)) * (1.0 - simpson_sum)
            } else {
                0.0
            };
            // Chao & Jost (2012) improved Good–Turing coverage estimator
            let f1_f = f1 as f64;
            let f2_f = f2 as f64;
            let c_hat = if f1 == 0 {
                1.0
            } else if f2 > 0 {
                1.0 - (f1_f / nf) * (((nf - 1.0) * f1_f) / ((nf - 1.0) * f1_f + 2.0 * f2_f))
            } else if n_tot > 1 {
                1.0 - (f1_f / nf) * (((nf - 1.0) * f1_f) / ((nf - 1.0) * f1_f + 2.0))
            } else {
                0.0
            };
            (hill1, hill2, j_even, h_pie, c_hat.clamp(0.0, 1.0))
        } else {
            (0.0, 0.0, 0.0, 0.0, 0.0)
        };

        alpha.push(SampleAlphaDiversity {
            sample_name: env_names[pos_s].clone(),
            total_reads_n: n_tot,
            observed_richness_q0: s_obs,
            singletons_f1: f1,
            doubletons_f2: f2,
            chao1_asymptotic: chao1,
            chao1_ci95: chao1_ci,
            hill_shannon_q1: q1,
            hill_simpson_q2: q2,
            pielou_evenness_j: pielou,
            hurlbert_pie: pie,
            sample_coverage: cov,
        });

        let mut depths_m = Vec::with_capacity(steps + 1);
        let mut expected_sm = Vec::with_capacity(steps + 1);
        depths_m.push(0);
        expected_sm.push(0.0);
        if n_tot > 0 {
            for k in 1..=steps {
                let m = ((k * n_tot) / steps).max(1);
                if depths_m.last().copied() != Some(m) {
                    depths_m.push(m);
                    expected_sm.push(hurlbert_expected_richness(&counts, m));
                }
            }
        }
        rarefaction_curves.push(RarefactionCurve {
            sample_name: env_names[pos_s].clone(),
            depths_m,
            expected_species_sm: expected_sm,
        });

        site_vectors.push(counts);
    }

    let mut bray = vec![vec![0.0_f64; n_sites]; n_sites];
    let mut jacc = vec![vec![0.0_f64; n_sites]; n_sites];

    for i in 0..n_sites {
        for j in (i + 1)..n_sites {
            let va = &site_vectors[i];
            let vb = &site_vectors[j];
            let mut sum_min = 0usize;
            let mut sum_tot = 0usize;
            let mut shared = 0usize;
            let mut union = 0usize;
            for (&ca, &cb) in va.iter().zip(vb) {
                sum_min += ca.min(cb);
                sum_tot += ca + cb;
                if ca > 0 || cb > 0 {
                    union += 1;
                    if ca > 0 && cb > 0 {
                        shared += 1;
                    }
                }
            }
            let bc_val = if sum_tot > 0 {
                1.0 - (2.0 * (sum_min as f64)) / (sum_tot as f64)
            } else {
                0.0
            };
            let jc_val = if union > 0 {
                1.0 - (shared as f64) / (union as f64)
            } else {
                0.0
            };
            bray[i][j] = bc_val;
            bray[j][i] = bc_val;
            jacc[i][j] = jc_val;
            jacc[j][i] = jc_val;
        }
    }

    CommunityDiversityReport {
        clean_asv_count: clean_asvs.len(),
        alpha,
        rarefaction_curves,
        sample_names: env_names,
        bray_curtis_matrix: bray,
        jaccard_matrix: jacc,
    }
}

/// Renders a two-panel `.svg` figure pairing Hurlbert (1971) rarefaction curves (`left`)
/// with Chao et al. (2014) Hill-number diversity profiles (`q = 0, 1, 2`, `right`).
pub fn render_edna_rarefaction_svg(diversity: &CommunityDiversityReport, title: &str) -> String {
    let width = 920.0_f64;
    let height = 440.0_f64;
    let p1_left = 68.0_f64;
    let p1_w = 400.0_f64;
    let p2_left = 545.0_f64;
    let p2_w = 325.0_f64;
    let top = 58.0_f64;
    let ph = 315.0_f64;

    let max_m = diversity
        .rarefaction_curves
        .iter()
        .filter_map(|c| c.depths_m.last().copied())
        .max()
        .unwrap_or(100)
        .max(1) as f64;
    let max_s = diversity
        .alpha
        .iter()
        .map(|a| a.chao1_asymptotic.max(a.observed_richness_q0 as f64))
        .fold(5.0_f64, f64::max)
        * 1.12;

    let map_x1 = |m: usize| p1_left + ((m as f64) / max_m).clamp(0.0, 1.0) * p1_w;
    let map_y = |s: f64| top + ph - (s / max_s).clamp(0.0, 1.0) * ph;

    const COLORS: &[&str] = &["#2563eb", "#16a34a", "#dc2626", "#9333ea", "#d97706"];
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
        "  <text x=\"{p1_left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Hurlbert Rarefaction &amp; Hill Diversity Profiles ({} Clean ASVs)</text>",
        xml_escape(title),
        diversity.clean_asv_count
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p1_left:.1}\" y=\"{top:.1}\" width=\"{p1_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p2_left:.1}\" y=\"{top:.1}\" width=\"{p2_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    for (idx, curve) in diversity.rarefaction_curves.iter().enumerate() {
        let col = COLORS[idx % COLORS.len()];
        let mut poly = String::new();
        for (&m, &sm) in curve.depths_m.iter().zip(&curve.expected_species_sm) {
            let _ = write!(poly, "{:.2},{:.2} ", map_x1(m), map_y(sm));
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.5\"/>",
            poly.trim()
        );

        let a = &diversity.alpha[idx];
        let x_q0 = p2_left + 35.0;
        let x_q1 = p2_left + 0.5 * p2_w;
        let x_q2 = p2_left + p2_w - 35.0;
        let _ = writeln!(
            svg,
            "  <polyline points=\"{x_q0:.1},{:.1} {x_q1:.1},{:.1} {x_q2:.1},{:.1}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.3\"/>",
            map_y(a.observed_richness_q0 as f64),
            map_y(a.hill_shannon_q1),
            map_y(a.hill_simpson_q2)
        );
        for (xp, yp) in [
            (x_q0, map_y(a.observed_richness_q0 as f64)),
            (x_q1, map_y(a.hill_shannon_q1)),
            (x_q2, map_y(a.hill_simpson_q2)),
        ] {
            let _ = writeln!(svg, "  <circle cx=\"{xp:.1}\" cy=\"{yp:.1}\" r=\"4.0\" fill=\"{col}\"/>");
        }
    }

    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Subsampled Sequencing Depth m (Reads)</text>",
        p1_left + 0.5 * p1_w,
        height - 18.0
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Hill Order: q=0 (Richness) → q=1 (exp H') → q=2 (1/D)</text>",
        p2_left + 0.5 * p2_w,
        height - 18.0
    );

    svg.push_str("</svg>\n");
    svg
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::edna_fastq::{
        demo_pnw_12s_edna_dataset, parse_barcode_sheet_csv, parse_fastq_str,
        run_edna_metabarcoding_pipeline, EdnaPipelineConfig,
    };

    #[test]
    fn computes_exact_hurlbert_rarefaction_chao1_and_hill_numbers() {
        let (fq, bc, db) = demo_pnw_12s_edna_dataset();
        let reads = parse_fastq_str(&fq).unwrap();
        let bcs = parse_barcode_sheet_csv(&bc).unwrap();
        let pipe = run_edna_metabarcoding_pipeline(&reads, &bcs, &db, &EdnaPipelineConfig::default());
        let div = analyze_edna_community_diversity(&pipe, 10);

        // Clean ASVs = 6 native PNW fishes/amphibians (excluding 1 chimera and 1 blank-flagged human ASV)
        assert_eq!(div.clean_asv_count, 6);
        assert_eq!(div.alpha.len(), 3);

        for a in &div.alpha {
            assert!(a.chao1_asymptotic >= (a.observed_richness_q0 as f64) - 1e-9);
            assert!(a.observed_richness_q0 as f64 >= a.hill_shannon_q1 - 1e-9);
            assert!(a.hill_shannon_q1 >= a.hill_simpson_q2 - 1e-9);
            assert!(a.sample_coverage > 0.85 && a.sample_coverage <= 1.0);
        }
    }
}
