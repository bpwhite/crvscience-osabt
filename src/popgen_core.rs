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

//! Tool #16 (`osabt popgen`) — Population genetics, neutrality & diploid `VCF` structure engine:
//! - **Nucleotide Diversity (`pi`), Watterson's `theta_W`, & Tajima's `D` Neutrality Test** —
//!   Tajima (1989) *Genetics* 123(3):585–595 (`https://pubmed.ncbi.nlm.nih.gov/2513255/`)
//! - **Pairwise Sequence Divergence (`d_XY`, `d_a`) & Hudson `F_ST`**
//! - **Unbiased Multilocus Diploid `F_IS`, `F_ST`, `F_IT` Variance Components** —
//!   Weir & Cockerham (1984) *Evolution* 38(6):1358–1370 (`https://pubmed.ncbi.nlm.nih.gov/28563791/`)
//! - **Exact Levene–Haldane Hardy–Weinberg Equilibrium (`HWE`) Test** —
//!   Wigginton, Cutler & Abecasis (2005) *Am. J. Hum. Genet.* 76(5):887–893 (`https://pubmed.ncbi.nlm.nih.gov/15789306/`).

use std::collections::BTreeMap;

/// Aligned DNA sequence tagged with individual sample ID and population/watershed label.
#[derive(Debug, Clone)]
pub struct PopSequenceRecord {
    pub sample_id: String,
    pub population: String,
    pub sequence: String,
}

/// Within-population sequence polymorphism and neutrality summary.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PopulationDiversityMetrics {
    pub population: String,
    pub sample_size_n: usize,
    pub sequence_len_bp: usize,
    pub num_haplotypes: usize,
    /// Nei (1987) unbiased haplotype diversity `H_d = (n / (n - 1)) * (1 - sum p_i^2)`.
    pub haplotype_diversity_hd: f64,
    pub segregating_sites_s: usize,
    /// Mean pairwise nucleotide differences `k_hat`.
    pub mean_pairwise_diff_k: f64,
    /// Per-site nucleotide diversity `pi = k_hat / L`.
    pub nucleotide_diversity_pi: f64,
    /// Per-site Watterson's (1975) `theta_W = (S / a_1) / L`.
    pub watterson_theta_w: f64,
    /// Tajima's (1989) `D` statistic (`< 0` excess rare alleles / expansion; `> 0` balancing selection / subdivision).
    pub tajima_d: f64,
    /// Two-sided normal/Beta tail p-value calibrated against Tajima (1989).
    pub tajima_d_pvalue: f64,
}

/// Pairwise sequence differentiation between two populations.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PairwisePopDivergence {
    pub pop_a: String,
    pub pop_b: String,
    /// Nei (1987) mean per-site divergence between sequences from `pop_a` and `pop_b`.
    pub d_xy: f64,
    /// Net nucleotide divergence `d_a = d_XY - 0.5 * (pi_a + pi_b)`.
    pub net_d_a: f64,
    /// Hudson, Slatkin & Maddison (1992) sequence `F_ST = (d_XY - pi_S) / d_XY` clamped to `[0, 1]`.
    pub hudson_fst: f64,
}

/// Full FASTA population-diversity & divergence report.
#[derive(Debug, Clone)]
pub struct PopGenFastaReport {
    pub total_sequences: usize,
    pub alignment_len_bp: usize,
    pub populations: Vec<String>,
    pub global_metrics: PopulationDiversityMetrics,
    pub per_population: Vec<PopulationDiversityMetrics>,
    pub pairwise_divergence: Vec<PairwisePopDivergence>,
}

/// Biallelic diploid `VCF` locus with per-sample genotypes (`Some(0)` = `0/0`, `Some(1)` = `0/1`, `Some(2)` = `1/1`, `None` = `./.`).
#[derive(Debug, Clone)]
pub struct DiploidSnpLocus {
    pub chrom: String,
    pub pos: u64,
    pub snp_id: String,
    pub ref_allele: String,
    pub alt_allele: String,
    pub genotypes: Vec<Option<u8>>,
}

/// Locus-specific Hardy–Weinberg equilibrium (`Wigginton et al. 2005`) & Weir–Cockerham `F_ST` diagnostics.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SnpLocusReport {
    pub snp_id: String,
    pub chrom: String,
    pub pos: u64,
    pub obs_hom_ref: usize,
    pub obs_het: usize,
    pub obs_hom_alt: usize,
    pub minor_allele_freq: f64,
    pub obs_heterozygosity: f64,
    pub exp_heterozygosity: f64,
    pub inbreeding_fis_pooled: f64,
    /// Wigginton, Cutler & Abecasis (2005) exact Levene–Haldane HWE p-value.
    pub hwe_exact_pvalue: f64,
    /// Weir & Cockerham (1984) locus-specific `theta` (`F_ST`).
    pub weir_cockerham_fst: f64,
}

/// Complete diploid `VCF` Weir–Cockerham `F`-statistics and exact `HWE` audit report.
#[derive(Debug, Clone)]
pub struct DiploidVcfReport {
    pub num_samples: usize,
    pub populations: Vec<String>,
    pub num_loci: usize,
    /// Multilocus ratio-of-averages `F_ST = sum(a) / sum(a + b + c)` (Weir & Cockerham 1984 Eq. 10).
    pub multilocus_fst: f64,
    /// Multilocus within-subpopulation inbreeding `F_IS = 1 - sum(c) / sum(b + c)`.
    pub multilocus_fis: f64,
    /// Multilocus total inbreeding `F_IT = 1 - sum(c) / sum(a + b + c)`.
    pub multilocus_fit: f64,
    pub pairwise_wc_fst: Vec<(String, String, f64)>,
    pub loci: Vec<SnpLocusReport>,
}

/// Parses an aligned multi-FASTA file where headers encode `>sample_id|population` (or `pop=Name`).
pub fn parse_popgen_fasta(raw_fasta: &str) -> Result<Vec<PopSequenceRecord>, String> {
    let mut records: Vec<PopSequenceRecord> = Vec::new();
    let mut cur_id = String::new();
    let mut cur_pop = String::new();
    let mut cur_seq = String::new();

    let flush = |id: &mut String,
                 pop: &mut String,
                 seq: &mut String,
                 out: &mut Vec<PopSequenceRecord>| {
        if !id.is_empty() && !seq.is_empty() {
            out.push(PopSequenceRecord {
                sample_id: std::mem::take(id),
                population: if pop.is_empty() {
                    String::from("Panmictic")
                } else {
                    std::mem::take(pop)
                },
                sequence: std::mem::take(seq),
            });
        }
    };

    for line in raw_fasta.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with(';') {
            continue;
        }
        if let Some(hdr) = t.strip_prefix('>') {
            flush(&mut cur_id, &mut cur_pop, &mut cur_seq, &mut records);
            if let Some((left, right)) = hdr.split_once('|') {
                cur_id = left.trim().to_string();
                cur_pop = right
                    .split_whitespace()
                    .next()
                    .unwrap_or("Panmictic")
                    .trim()
                    .to_string();
            } else {
                let mut toks = hdr.split_whitespace();
                cur_id = toks.next().unwrap_or("seq").to_string();
                cur_pop = toks
                    .find_map(|tk| tk.strip_prefix("pop="))
                    .unwrap_or("Panmictic")
                    .to_string();
            }
        } else {
            for ch in t.chars() {
                if !ch.is_ascii_whitespace() {
                    cur_seq.push(ch.to_ascii_uppercase());
                }
            }
        }
    }
    flush(&mut cur_id, &mut cur_pop, &mut cur_seq, &mut records);

    if records.len() < 2 {
        return Err(String::from(
            "Population FASTA requires at least 2 aligned sequences.",
        ));
    }
    let len0 = records[0].sequence.len();
    if len0 == 0 {
        return Err(String::from("Aligned FASTA sequences are empty."));
    }
    for r in &records {
        if r.sequence.len() != len0 {
            return Err(format!(
                "Sequence '{}' length ({}) differs from alignment length ({}).",
                r.sample_id,
                r.sequence.len(),
                len0
            ));
        }
    }
    Ok(records)
}

/// Computes global and per-population nucleotide diversity (`pi`), `theta_W`, `H_d`,
/// Tajima's (1989) `D`, and pairwise `d_XY` / `d_a` / Hudson `F_ST` across populations.
pub fn analyze_popgen_fasta(records: &[PopSequenceRecord]) -> Result<PopGenFastaReport, String> {
    if records.len() < 2 {
        return Err(String::from("Need at least 2 sequences for PopGen analysis."));
    }
    let len_bp = records[0].sequence.len();

    let all_seqs: Vec<&str> = records.iter().map(|r| r.sequence.as_str()).collect();
    let global_metrics = compute_group_diversity("All_Pooled", &all_seqs, len_bp);

    let mut pop_groups: BTreeMap<String, Vec<&str>> = BTreeMap::new();
    for r in records {
        pop_groups
            .entry(r.population.clone())
            .or_default()
            .push(r.sequence.as_str());
    }

    let populations: Vec<String> = pop_groups.keys().cloned().collect();
    let mut per_population = Vec::with_capacity(populations.len());
    for pop in &populations {
        let grp = &pop_groups[pop];
        per_population.push(compute_group_diversity(pop, grp, len_bp));
    }

    let mut pairwise = Vec::new();
    for i in 0..populations.len() {
        for j in (i + 1)..populations.len() {
            let seqs_a = &pop_groups[&populations[i]];
            let seqs_b = &pop_groups[&populations[j]];
            let mut diff_sum = 0usize;
            let mut pair_cnt = 0usize;
            for sa in seqs_a {
                for sb in seqs_b {
                    diff_sum += pairwise_valid_diffs(sa.as_bytes(), sb.as_bytes());
                    pair_cnt += 1;
                }
            }
            let d_xy = if pair_cnt > 0 && len_bp > 0 {
                (diff_sum as f64) / ((pair_cnt * len_bp) as f64)
            } else {
                0.0
            };
            let pi_s = 0.5
                * (per_population[i].nucleotide_diversity_pi
                    + per_population[j].nucleotide_diversity_pi);
            let net_da = d_xy - pi_s;
            let h_fst = if d_xy > 1e-12 {
                ((d_xy - pi_s) / d_xy).clamp(0.0, 1.0)
            } else {
                0.0
            };
            pairwise.push(PairwisePopDivergence {
                pop_a: populations[i].clone(),
                pop_b: populations[j].clone(),
                d_xy,
                net_d_a: net_da,
                hudson_fst: h_fst,
            });
        }
    }

    Ok(PopGenFastaReport {
        total_sequences: records.len(),
        alignment_len_bp: len_bp,
        populations,
        global_metrics,
        per_population,
        pairwise_divergence: pairwise,
    })
}

fn compute_group_diversity(
    label: &str,
    seqs: &[&str],
    len_bp: usize,
) -> PopulationDiversityMetrics {
    let n = seqs.len();
    let mut hap_counts: BTreeMap<&str, usize> = BTreeMap::new();
    for &s in seqs {
        *hap_counts.entry(s).or_insert(0) += 1;
    }
    let num_haps = hap_counts.len();
    let hd = if n > 1 {
        let nf = n as f64;
        let sum_p2: f64 = hap_counts
            .values()
            .map(|&c| {
                let p = (c as f64) / nf;
                p * p
            })
            .sum();
        (nf / (nf - 1.0)) * (1.0 - sum_p2)
    } else {
        0.0
    };

    if n < 2 || len_bp == 0 {
        return PopulationDiversityMetrics {
            population: label.to_string(),
            sample_size_n: n,
            sequence_len_bp: len_bp,
            num_haplotypes: num_haps,
            haplotype_diversity_hd: hd,
            segregating_sites_s: 0,
            mean_pairwise_diff_k: 0.0,
            nucleotide_diversity_pi: 0.0,
            watterson_theta_w: 0.0,
            tajima_d: 0.0,
            tajima_d_pvalue: 1.0,
        };
    }

    let mut seg_s = 0usize;
    for col in 0..len_bp {
        let mut mask = 0u8;
        for s in seqs {
            mask |= match s.as_bytes()[col] {
                b'A' => 1,
                b'C' => 2,
                b'G' => 4,
                b'T' => 8,
                _ => 0,
            };
        }
        if mask.count_ones() >= 2 {
            seg_s += 1;
        }
    }

    let mut total_diffs = 0usize;
    let mut pairs = 0usize;
    for i in 0..n {
        for j in (i + 1)..n {
            total_diffs += pairwise_valid_diffs(seqs[i].as_bytes(), seqs[j].as_bytes());
            pairs += 1;
        }
    }
    let k_hat = (total_diffs as f64) / (pairs as f64);
    let pi = k_hat / (len_bp as f64);

    let mut a1 = 0.0_f64;
    let mut a2 = 0.0_f64;
    for i in 1..n {
        let inv = 1.0 / (i as f64);
        a1 += inv;
        a2 += inv * inv;
    }
    let s_f = seg_s as f64;
    let theta_w_seq = s_f / a1;
    let theta_w_site = theta_w_seq / (len_bp as f64);

    let (tajima_d, pval) = if seg_s > 0 && n >= 4 {
        let nf = n as f64;
        let b1 = (nf + 1.0) / (3.0 * (nf - 1.0));
        let b2 = (2.0 * (nf * nf + nf + 3.0)) / (9.0 * nf * (nf - 1.0));
        let c1 = b1 - 1.0 / a1;
        let c2 = b2 - (nf + 2.0) / (a1 * nf) + a2 / (a1 * a1);
        let e1 = c1 / a1;
        let e2 = c2 / (a1 * a1 + a2);
        let var_d = (e1 * s_f + e2 * s_f * (s_f - 1.0)).max(1e-12);
        let d_val = (k_hat - theta_w_seq) / var_d.sqrt();
        (d_val, approximate_tajima_pvalue(d_val))
    } else {
        (0.0, 1.0)
    };

    PopulationDiversityMetrics {
        population: label.to_string(),
        sample_size_n: n,
        sequence_len_bp: len_bp,
        num_haplotypes: num_haps,
        haplotype_diversity_hd: hd,
        segregating_sites_s: seg_s,
        mean_pairwise_diff_k: k_hat,
        nucleotide_diversity_pi: pi,
        watterson_theta_w: theta_w_site,
        tajima_d,
        tajima_d_pvalue: pval,
    }
}

#[inline]
pub(crate) fn pairwise_valid_diffs(a: &[u8], b: &[u8]) -> usize {
    a.iter()
        .zip(b)
        .filter(|&(&x, &y)| {
            matches!(x, b'A' | b'C' | b'G' | b'T')
                && matches!(y, b'A' | b'C' | b'G' | b'T')
                && x != y
        })
        .count()
}

fn approximate_tajima_pvalue(d: f64) -> f64 {
    let z = d.abs();
    let t = 1.0 / (1.0 + 0.231_641_9 * z);
    let poly = t
        * (0.319_381_530
            + t * (-0.356_563_782
                + t * (1.781_477_937 + t * (-1.821_255_978 + t * 1.330_274_429))));
    let pdf = (-0.5 * z * z).exp() / (2.0 * std::f64::consts::PI).sqrt();
    (2.0 * pdf * poly).clamp(0.0, 1.0)
}

/// Computes Wigginton, Cutler & Abecasis (2005) exact two-sided Levene–Haldane
/// Hardy–Weinberg equilibrium p-value from diploid counts `(obs_hom1, obs_hets, obs_hom2)`.
pub fn calc_exact_hwe_wigginton(obs_hom1: usize, obs_hets: usize, obs_hom2: usize) -> f64 {
    let n = obs_hom1 + obs_hets + obs_hom2;
    if n == 0 {
        return 1.0;
    }
    let a1 = 2 * obs_hom1 + obs_hets;
    let a2 = 2 * obs_hom2 + obs_hets;
    let rare = a1.min(a2);
    let common = 2 * n - rare;
    if rare == 0 {
        return 1.0;
    }
    if obs_hets > rare || (obs_hets % 2) != (rare % 2) {
        return 0.0;
    }

    let mut probs = vec![0.0_f64; rare + 1];
    let mut mid = (rare * common) / (2 * n);
    if (mid % 2) != (rare % 2) {
        mid += 1;
    }
    mid = mid.min(rare);
    probs[mid] = 1.0;
    let mut sum = 1.0_f64;

    // Step downward by 2 from modal heterozygote count `mid`
    let mut curr_hets = mid;
    let mut curr_hom_r = (rare - mid) / 2;
    let mut curr_hom_c = (common - mid) / 2;
    while curr_hets >= 2 {
        let next_p = probs[curr_hets] * ((curr_hets * (curr_hets - 1)) as f64)
            / (4.0 * ((curr_hom_r + 1) * (curr_hom_c + 1)) as f64);
        probs[curr_hets - 2] = next_p;
        sum += next_p;
        curr_hets -= 2;
        curr_hom_r += 1;
        curr_hom_c += 1;
    }

    // Step upward by 2 from modal heterozygote count `mid`
    curr_hets = mid;
    curr_hom_r = (rare - mid) / 2;
    curr_hom_c = (common - mid) / 2;
    while curr_hets + 2 <= rare && curr_hom_r >= 1 && curr_hom_c >= 1 {
        let next_p = probs[curr_hets] * (4.0 * (curr_hom_r * curr_hom_c) as f64)
            / (((curr_hets + 2) * (curr_hets + 1)) as f64);
        probs[curr_hets + 2] = next_p;
        sum += next_p;
        curr_hets += 2;
        curr_hom_r -= 1;
        curr_hom_c -= 1;
    }

    let target = probs[obs_hets];
    let mut p_hwe = 0.0_f64;
    let mut h = rare % 2;
    while h <= rare {
        if probs[h] <= target * (1.0 + 1e-12) {
            p_hwe += probs[h] / sum;
        }
        h += 2;
    }
    p_hwe.clamp(0.0, 1.0)
}

/// Parses a diploid biallelic `VCF` v4.x file where sample columns are named `sample|population`
/// (or `population_sample`). Returns `(sample_populations, snp_loci)`.
pub fn parse_popgen_vcf(raw_vcf: &str) -> Result<(Vec<String>, Vec<DiploidSnpLocus>), String> {
    let mut sample_pops: Vec<String> = Vec::new();
    let mut loci: Vec<DiploidSnpLocus> = Vec::new();

    for (line_idx, line) in raw_vcf.lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with("##") {
            continue;
        }
        if let Some(hdr) = t.strip_prefix("#CHROM") {
            let cols: Vec<&str> = hdr.split_whitespace().collect();
            // Columns after `#CHROM`: POS(0) ID(1) REF(2) ALT(3) QUAL(4) FILTER(5) INFO(6) FORMAT(7) samples(8..)
            if cols.len() < 10 {
                return Err(String::from(
                    "VCF #CHROM header must include FORMAT and at least 2 sample columns.",
                ));
            }
            for &scol in &cols[8..] {
                let pop = if let Some((_, p)) = scol.split_once('|') {
                    p.to_string()
                } else if let Some((p, _)) = scol.rsplit_once('_') {
                    p.to_string()
                } else {
                    String::from("Pop_1")
                };
                sample_pops.push(pop);
            }
            continue;
        }

        let cols: Vec<&str> = t.split_whitespace().collect();
        if cols.len() < 9 + sample_pops.len() || sample_pops.is_empty() {
            return Err(format!(
                "VCF line {}: malformed record or missing #CHROM header.",
                line_idx + 1
            ));
        }
        // Retain biallelic SNP records
        if cols[4].contains(',') {
            continue;
        }
        let pos = cols[1]
            .parse::<u64>()
            .map_err(|_| format!("VCF line {}: invalid POS '{}'.", line_idx + 1, cols[1]))?;

        let mut gts = Vec::with_capacity(sample_pops.len());
        for &gt_field in &cols[9..(9 + sample_pops.len())] {
            let gt_tok = gt_field.split(':').next().unwrap_or("./.");
            let call = match gt_tok {
                "0/0" | "0|0" => Some(0u8),
                "0/1" | "1/0" | "0|1" | "1|0" => Some(1u8),
                "1/1" | "1|1" => Some(2u8),
                _ => None,
            };
            gts.push(call);
        }

        loci.push(DiploidSnpLocus {
            chrom: cols[0].to_string(),
            pos,
            snp_id: cols[2].to_string(),
            ref_allele: cols[3].to_string(),
            alt_allele: cols[4].to_string(),
            genotypes: gts,
        });
    }

    if loci.is_empty() {
        return Err(String::from("VCF contained no biallelic diploid loci."));
    }
    Ok((sample_pops, loci))
}

/// Computes Weir & Cockerham (1984) multilocus & per-locus `F_ST`, `F_IS`, `F_IT` variance components
/// and Wigginton et al. (2005) exact Hardy–Weinberg equilibrium p-values across diploid `VCF` loci.
pub fn analyze_diploid_vcf(
    sample_pops: &[String],
    loci: &[DiploidSnpLocus],
) -> Result<DiploidVcfReport, String> {
    if sample_pops.len() < 2 || loci.is_empty() {
        return Err(String::from("Diploid VCF requires >= 2 samples and >= 1 locus."));
    }

    let mut pop_indices: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for (idx, p) in sample_pops.iter().enumerate() {
        pop_indices.entry(p.clone()).or_default().push(idx);
    }
    let populations: Vec<String> = pop_indices.keys().cloned().collect();
    let idx_groups: Vec<Vec<usize>> = populations
        .iter()
        .map(|p| pop_indices[p].clone())
        .collect();

    let mut sum_a = 0.0_f64;
    let mut sum_b = 0.0_f64;
    let mut sum_c = 0.0_f64;
    let mut locus_reports = Vec::with_capacity(loci.len());

    for loc in loci {
        let mut hom0 = 0usize;
        let mut het = 0usize;
        let mut hom2 = 0usize;
        for &g in &loc.genotypes {
            match g {
                Some(0) => hom0 += 1,
                Some(1) => het += 1,
                Some(2) => hom2 += 1,
                _ => {}
            }
        }
        let n_tot = hom0 + het + hom2;
        let p_alt = if n_tot > 0 {
            ((2 * hom2 + het) as f64) / ((2 * n_tot) as f64)
        } else {
            0.0
        };
        let maf = p_alt.min(1.0 - p_alt);
        let h_obs = if n_tot > 0 {
            (het as f64) / (n_tot as f64)
        } else {
            0.0
        };
        let h_exp = 2.0 * p_alt * (1.0 - p_alt);
        let fis_pool = if h_exp > 1e-12 {
            1.0 - h_obs / h_exp
        } else {
            0.0
        };
        let hwe_p = calc_exact_hwe_wigginton(hom0, het, hom2);

        let (a_l, b_l, c_l) = weir_cockerham_locus_abc(&loc.genotypes, &idx_groups);
        sum_a += a_l;
        sum_b += b_l;
        sum_c += c_l;
        let denom_l = a_l + b_l + c_l;
        let wc_fst = if denom_l > 1e-12 {
            (a_l / denom_l).clamp(-1.0, 1.0)
        } else {
            0.0
        };

        locus_reports.push(SnpLocusReport {
            snp_id: loc.snp_id.clone(),
            chrom: loc.chrom.clone(),
            pos: loc.pos,
            obs_hom_ref: hom0,
            obs_het: het,
            obs_hom_alt: hom2,
            minor_allele_freq: maf,
            obs_heterozygosity: h_obs,
            exp_heterozygosity: h_exp,
            inbreeding_fis_pooled: fis_pool,
            hwe_exact_pvalue: hwe_p,
            weir_cockerham_fst: wc_fst,
        });
    }

    let tot_abc = sum_a + sum_b + sum_c;
    let multilocus_fst = if tot_abc > 1e-12 {
        (sum_a / tot_abc).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let multilocus_fis = if (sum_b + sum_c) > 1e-12 {
        1.0 - sum_c / (sum_b + sum_c)
    } else {
        0.0
    };
    let multilocus_fit = if tot_abc > 1e-12 {
        1.0 - sum_c / tot_abc
    } else {
        0.0
    };

    let mut pairwise_wc = Vec::new();
    for i in 0..populations.len() {
        for j in (i + 1)..populations.len() {
            let sub_groups = vec![idx_groups[i].clone(), idx_groups[j].clone()];
            let mut pa = 0.0_f64;
            let mut pabc = 0.0_f64;
            for loc in loci {
                let (al, bl, cl) = weir_cockerham_locus_abc(&loc.genotypes, &sub_groups);
                pa += al;
                pabc += al + bl + cl;
            }
            let pfst = if pabc > 1e-12 {
                (pa / pabc).clamp(0.0, 1.0)
            } else {
                0.0
            };
            pairwise_wc.push((populations[i].clone(), populations[j].clone(), pfst));
        }
    }

    Ok(DiploidVcfReport {
        num_samples: sample_pops.len(),
        populations,
        num_loci: loci.len(),
        multilocus_fst,
        multilocus_fis,
        multilocus_fit,
        pairwise_wc_fst: pairwise_wc,
        loci: locus_reports,
    })
}

/// Computes Weir & Cockerham (1984) variance components `(a, b, c)` for one biallelic locus.
fn weir_cockerham_locus_abc(genotypes: &[Option<u8>], groups: &[Vec<usize>]) -> (f64, f64, f64) {
    let mut n_vec = Vec::with_capacity(groups.len());
    let mut p_vec = Vec::with_capacity(groups.len());
    let mut h_vec = Vec::with_capacity(groups.len());

    for grp in groups {
        let mut n_i = 0usize;
        let mut alt_alleles = 0usize;
        let mut hets = 0usize;
        for &idx in grp {
            if let Some(gt) = genotypes[idx] {
                n_i += 1;
                alt_alleles += gt as usize;
                if gt == 1 {
                    hets += 1;
                }
            }
        }
        if n_i >= 1 {
            n_vec.push(n_i as f64);
            p_vec.push((alt_alleles as f64) / ((2 * n_i) as f64));
            h_vec.push((hets as f64) / (n_i as f64));
        }
    }

    let r = n_vec.len();
    if r < 2 {
        return (0.0, 0.0, 0.0);
    }
    let rf = r as f64;
    let sum_n: f64 = n_vec.iter().sum();
    let n_bar = sum_n / rf;
    if n_bar <= 1.0 {
        return (0.0, 0.0, 0.0);
    }
    let sum_n2: f64 = n_vec.iter().map(|&ni| ni * ni).sum();
    let n_c = (sum_n - sum_n2 / sum_n) / (rf - 1.0);
    if n_c <= 1e-9 {
        return (0.0, 0.0, 0.0);
    }

    let p_bar: f64 = n_vec
        .iter()
        .zip(&p_vec)
        .map(|(&ni, &pi)| ni * pi)
        .sum::<f64>()
        / sum_n;
    let s2: f64 = n_vec
        .iter()
        .zip(&p_vec)
        .map(|(&ni, &pi)| ni * (pi - p_bar).powi(2))
        .sum::<f64>()
        / ((rf - 1.0) * n_bar);
    let h_bar: f64 = n_vec
        .iter()
        .zip(&h_vec)
        .map(|(&ni, &hi)| ni * hi)
        .sum::<f64>()
        / sum_n;

    let inner = p_bar * (1.0 - p_bar) - ((rf - 1.0) / rf) * s2;
    let a = (n_bar / n_c) * (s2 - (1.0 / (n_bar - 1.0)) * (inner - 0.25 * h_bar));
    let b = (n_bar / (n_bar - 1.0)) * (inner - ((2.0 * n_bar - 1.0) / (4.0 * n_bar)) * h_bar);
    let c = 0.5 * h_bar;
    // Multiply by 2 so `a + b + c` sums over both REF and ALT indicator alleles (cancels in ratios)
    (2.0 * a, 2.0 * b, 2.0 * c)
}

/// Built-in 28-sequence × 60-bp aligned Pacific Northwest `cytb` dataset across 3 watersheds
/// (`Cascade_Headwaters`, `Willamette_Foothills`, `Coast_Range`), exhibiting an ancestral
/// Cascade Median-Joining hub (`MV_1`) and a Coast Range star-expansion (`Tajima's D < -1.5`).
pub fn demo_pnw_popgen_fasta() -> String {
    let anc = b"ATGACCAACATCCGAAAAACCCACCCACTACTAAAAATTATTAACAACTCATTCATTGAC";
    let mutate = |positions: &[(usize, u8)]| -> String {
        let mut s = anc.to_vec();
        for &(p, b) in positions {
            s[p] = b;
        }
        String::from_utf8(s).unwrap()
    };

    // Cascade_Headwaters: 3 deep clades radiating symmetrically from unsampled ancestor `anc` (MV_1)
    let h1 = mutate(&[(3, b'G'), (8, b'T')]);
    let h2 = mutate(&[(18, b'G'), (23, b'T')]);
    let h3 = mutate(&[(33, b'G'), (38, b'C')]);

    // Willamette_Foothills: shares H2 and steps toward Coast Range via H4 & H5
    let h4 = mutate(&[(18, b'G'), (23, b'T'), (44, b'T')]);
    let h5 = mutate(&[(18, b'G'), (23, b'T'), (44, b'T'), (49, b'T')]);

    // Coast_Range: central haplotype H6 (n=8) + 4 singleton 1-bp star-expansion satellites (H7..H10)
    let h6 = mutate(&[(18, b'G'), (23, b'T'), (52, b'C'), (55, b'C')]);
    let h7 = mutate(&[(1, b'C'), (18, b'G'), (23, b'T'), (52, b'C'), (55, b'C')]);
    let h8 = mutate(&[(6, b'G'), (18, b'G'), (23, b'T'), (52, b'C'), (55, b'C')]);
    let h9 = mutate(&[(13, b'A'), (18, b'G'), (23, b'T'), (52, b'C'), (55, b'C')]);
    let h10 = mutate(&[(28, b'C'), (18, b'G'), (23, b'T'), (52, b'C'), (55, b'C')]);

    let mut out = String::with_capacity(4096);
    let mut emit = |prefix: &str, pop: &str, seq: &str, count: usize, start_no: usize| {
        for k in 0..count {
            out.push_str(&format!(">{}_{:02}|{}\n{}\n", prefix, start_no + k, pop, seq));
        }
    };

    emit("CASC", "Cascade_Headwaters", &h1, 3, 1);
    emit("CASC", "Cascade_Headwaters", &h2, 3, 4);
    emit("CASC", "Cascade_Headwaters", &h3, 3, 7);

    emit("WILL", "Willamette_Foothills", &h2, 1, 1);
    emit("WILL", "Willamette_Foothills", &h4, 4, 2);
    emit("WILL", "Willamette_Foothills", &h5, 2, 6);

    emit("COAST", "Coast_Range", &h6, 8, 1);
    emit("COAST", "Coast_Range", &h7, 1, 9);
    emit("COAST", "Coast_Range", &h8, 1, 10);
    emit("COAST", "Coast_Range", &h9, 1, 11);
    emit("COAST", "Coast_Range", &h10, 1, 12);

    out
}

/// Built-in 18-individual × 6-SNP diploid `VCF` v4.2 across `Coast_Range`, `Willamette`, and `Cascade`
/// including neutral HWE-conforming loci, a strongly differentiated adaptive SNP (`F_ST > 0.55`),
/// and a Wahlund/null-allele heterozygote-deficit locus (`HWE p < 0.01`).
pub fn demo_pnw_popgen_vcf() -> &'static str {
    "\
##fileformat=VCFv4.2
##source=OSABT_PopGen_v1.6.0
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS01|Coast\tS02|Coast\tS03|Coast\tS04|Coast\tS05|Coast\tS06|Coast\tS07|Willamette\tS08|Willamette\tS09|Willamette\tS10|Willamette\tS11|Willamette\tS12|Willamette\tS13|Cascade\tS14|Cascade\tS15|Cascade\tS16|Cascade\tS17|Cascade\tS18|Cascade
chr1\t10420\trs_hypoxia_1\tA\tG\t99\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0\t0/1\t0/1\t1/1\t0/1\t0/0\t0/1\t1/1\t1/1\t1/1\t0/1\t1/1\t1/1
chr1\t28910\trs_neutral_2\tC\tT\t99\tPASS\t.\tGT\t0/1\t0/0\t0/1\t1/1\t0/0\t0/1\t0/1\t0/0\t1/1\t0/1\t0/0\t0/1\t0/0\t0/1\t0/1\t1/1\t0/0\t0/1
chr2\t51200\trs_thermal_3\tG\tA\t99\tPASS\t.\tGT\t1/1\t1/1\t0/1\t1/1\t1/1\t1/1\t0/1\t0/1\t0/0\t0/1\t1/1\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t0/0
chr2\t77450\trs_null_hwd_4\tT\tC\t99\tPASS\t.\tGT\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t0/0\t0/0\t0/0\t1/1\t1/1\t1/1\t0/0\t0/0\t0/1\t1/1\t1/1\t1/1
chr3\t14300\trs_neutral_5\tA\tC\t99\tPASS\t.\tGT\t0/0\t0/1\t0/0\t0/1\t0/0\t1/1\t0/0\t0/1\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t0/0\t1/1\t0/0
chr3\t63880\trs_elev_6\tG\tT\t99\tPASS\t.\tGT\t0/0\t0/0\t0/0\t0/0\t0/1\t0/0\t0/1\t0/0\t0/1\t1/1\t0/1\t0/0\t1/1\t0/1\t1/1\t1/1\t1/1\t0/1
"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_tajima_d_weir_cockerham_fst_and_exact_wigginton_hwe() {
        let fasta = demo_pnw_popgen_fasta();
        let recs = parse_popgen_fasta(&fasta).unwrap();
        let rep = analyze_popgen_fasta(&recs).unwrap();

        let coast = rep
            .per_population
            .iter()
            .find(|p| p.population == "Coast_Range")
            .unwrap();
        let casc = rep
            .per_population
            .iter()
            .find(|p| p.population == "Cascade_Headwaters")
            .unwrap();
        // Star-expansion excess of singletons in Coast_Range yields strongly negative Tajima's D
        assert!(coast.tajima_d < -1.4);
        // Deep equal-frequency clades in Cascade_Headwaters yield positive Tajima's D
        assert!(casc.tajima_d > 1.0);

        let (pops, loci) = parse_popgen_vcf(demo_pnw_popgen_vcf()).unwrap();
        let vcf_rep = analyze_diploid_vcf(&pops, &loci).unwrap();
        assert!(vcf_rep.multilocus_fst > 0.20);
        // Locus rs_null_hwd_4 has severe heterozygote deficit -> HWE exact p < 0.005
        let hwd = vcf_rep
            .loci
            .iter()
            .find(|l| l.snp_id == "rs_null_hwd_4")
            .unwrap();
        assert!(hwd.hwe_exact_pvalue < 0.005);
    }
}
