// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #11 (`osabt infoevo`) — Information-Theoretic Model of Evolution (`ITME`),
//! Sequence Logos, Genomic Physical Complexity, Epistasis `MIp`, & Selection Regret Engine:
//! - **Sequence Conservation & Logos (`R_seq`)** — Schneider et al. (1986) *J. Mol. Biol.* 188:415–431;
//!   Schneider & Stephens (1990) *Nucleic Acids Res.* 18(20):6097–6100
//! - **Genomic Physical Complexity (`C`)** — Adami, Ofria & Collier (2000) *PNAS* 97(9):4463–4468
//! - **APC-Corrected Epistasis Mutual Information (`MIp`)** — Dunn, Wahl & Gloor (2008) *Bioinformatics* 24(3):333–340
//! - **Substitutional Load & Information Gain (`I_T <= L_sub`)** — Kimura (1961) *Genet. Res.* 2:127–140

use std::f64::consts::LN_2;

/// Molecular alphabet of a multiple sequence alignment.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BioAlphabet {
    Dna,
    Rna,
    Protein,
}

impl BioAlphabet {
    pub fn symbols(self) -> &'static [u8] {
        match self {
            Self::Dna => b"ACGT",
            Self::Rna => b"ACGU",
            Self::Protein => b"ACDEFGHIKLMNPQRSTVWY",
        }
    }

    pub fn cardinality(self) -> usize {
        self.symbols().len()
    }

    /// Maximum Shannon uncertainty per site `log2(|A|)` (`2.0` bits for DNA/RNA, `4.321928` bits for Protein).
    pub fn max_bits(self) -> f64 {
        (self.cardinality() as f64).log2()
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Dna => "DNA (4-state, 2.00 bits max)",
            Self::Rna => "RNA (4-state, 2.00 bits max)",
            Self::Protein => "Protein (20-state, 4.32 bits max)",
        }
    }

    #[inline]
    pub fn symbol_index(self, b: u8) -> Option<usize> {
        let up = b.to_ascii_uppercase();
        self.symbols().iter().position(|&s| s == up)
    }
}

/// Single aligned sequence in a [`SequenceAlignment`].
#[derive(Debug, Clone)]
pub struct AlignedBioSeq {
    pub id: String,
    pub clade_tag: Option<String>,
    pub residues: Vec<u8>,
}

/// Validated rectangular multiple sequence alignment (`MSA`).
#[derive(Debug, Clone)]
pub struct SequenceAlignment {
    pub title: String,
    pub alphabet: BioAlphabet,
    pub length: usize,
    pub records: Vec<AlignedBioSeq>,
}

/// Per-column Shannon uncertainty and Schneider & Stephens (1990) bit-height stack.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SiteInformationProfile {
    /// 1-based alignment column position.
    pub position_1based: usize,
    pub non_gap_count: usize,
    pub gap_fraction: f64,
    /// Observed Shannon entropy `H(l) = -sum f_b log2(f_b)` in bits.
    pub shannon_bits: f64,
    /// Basharin (1959) / Schneider et al. (1986) small-sample error `e_n = (s - 1) / (2 * N * ln 2)`.
    pub small_sample_error_bits: f64,
    /// Uncorrected Adami per-site information `log2(s) - H(l)` (weighted by `1 - gap_fraction`).
    pub r_sequence_raw_bits: f64,
    /// Small-sample-corrected Schneider `R_seq(l) = max(0, log2(s) - (H(l) + e_n)) * (1 - gap_fraction)`.
    pub r_sequence_bits: f64,
    pub consensus_symbol: char,
    /// Sorted ascending `(symbol, height_in_bits)` for drawing stacked sequence logos (tallest drawn on top).
    pub letter_heights_bits: Vec<(char, f64)> ,
}

/// Locus-wide Schneider `R_sequence` and Adami et al. (2000) Genomic Physical Complexity report.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AlignmentComplexityReport {
    pub title: String,
    pub alphabet: BioAlphabet,
    pub num_sequences: usize,
    pub alignment_length: usize,
    pub max_bits_per_site: f64,
    /// Sum of Schneider `R_seq(l)` across all columns (`bits`).
    pub total_r_sequence_bits: f64,
    /// Adami et al. (2000) physical complexity `C = sum_l (log2(s) - H(l)) * (1 - gap_l)` (`bits`).
    pub adami_physical_complexity_bits: f64,
    /// Fraction of maximum theoretical sequence capacity occupied by environmental information `C / (L * log2(s))`.
    pub complexity_density: f64,
    /// Effective number of strictly invariant monomer equivalents (`C / log2(s)`).
    pub effective_conserved_sites: f64,
    /// Optional Schneider (1986) `R_frequency = log2(genome_size / num_targets)` (`bits`).
    pub r_frequency_target_bits: Option<f64>,
    pub sites: Vec<SiteInformationProfile>,
}

/// Pairwise epistatic / coevolving column pair scored via Dunn, Wahl & Gloor (2008) `MIp`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EpistasisPairScore {
    /// 1-based column index `i`.
    pub pos_i_1based: usize,
    /// 1-based column index `j` (`j > i`).
    pub pos_j_1based: usize,
    /// Raw Shannon mutual information `MI(i, j)` in bits.
    pub mi_bits: f64,
    /// Dunn et al. (2008) Average Product Correction `MI_bar(i) * MI_bar(j) / MI_bar`.
    pub apc_bits: f64,
    /// Phylogeny/entropy-corrected `MIp(i, j) = MI(i, j) - APC(i, j)` in bits.
    pub mip_bits: f64,
    /// Joint-entropy normalized mutual information `MI(i, j) / H(i, j)` in `[0, 1]`.
    pub normalized_mir: f64,
    /// Standardized Z-score of `MIp(i, j)` across all evaluated pairs.
    pub z_score: f64,
    /// Top observed joint state combinations (e.g. `"G-C (25%), A-T (25%)"`).
    pub top_covarying_states: String,
}

/// Full `L x L` mutual information & APC-corrected epistasis analysis.
#[derive(Debug, Clone)]
pub struct EpistasisMatrixReport {
    pub alignment_length: usize,
    pub mean_raw_mi_bits: f64,
    pub mip_std_dev: f64,
    /// Flattened `L x L` `MIp` matrix (`mip_matrix[i * L + j]`).
    pub mip_matrix: Vec<f64>,
    /// Flattened `L x L` raw `MI` matrix (`mi_matrix[i * L + j]`).
    pub mi_matrix: Vec<f64>,
    pub top_pairs: Vec<EpistasisPairScore>,
}

/// Per-site informational divergence between two clades or ecological populations.
#[derive(Debug, Clone)]
pub struct SiteCladeDivergence {
    pub position_1based: usize,
    /// Symmetric Jensen–Shannon divergence `D_JS(P || Q)` in `[0.0, 1.0]` bits.
    pub js_divergence_bits: f64,
    /// Asymmetric Kullback–Leibler divergence `D_KL(P || Q)` (Laplace-smoothed) in bits.
    pub kl_p_from_q_bits: f64,
    pub clade_a_consensus: char,
    pub clade_b_consensus: char,
}

/// Two-population evolutionary divergence summary across an alignment.
#[derive(Debug, Clone)]
pub struct CladeDivergenceReport {
    pub clade_a_name: String,
    pub clade_b_name: String,
    pub clade_a_n: usize,
    pub clade_b_n: usize,
    pub total_js_bits: f64,
    pub mean_js_bits: f64,
    pub total_kl_bits: f64,
    pub sites: Vec<SiteCladeDivergence>,
}

/// Single generation snapshot along a Kimura (1961) replicator selection trajectory.
#[derive(Debug, Clone)]
pub struct SelectionTrajectoryStep {
    pub generation: usize,
    pub fittest_freq: f64,
    pub mean_malthusian_fitness: f64,
    pub population_entropy_bits: f64,
    /// Accumulated Kullback–Leibler genetic information `I_t = D_KL(p_t || p_0)` in bits.
    pub kl_info_gain_bits: f64,
    /// Cumulative substitutional load / fitness regret `L_sub,t` in bits (`I_t <= L_sub,t`).
    pub cumulative_substitutional_load_bits: f64,
    /// Instantaneous conversion efficiency `I_t / L_sub,t` (`-> 1.0` asymptotically as fittest type fixes).
    pub selection_efficiency: f64,
}

/// Summary of Kimura (1961) adaptive genetic information accumulation & substitutional load.
#[derive(Debug, Clone)]
pub struct KimuraSelectionReport {
    pub num_alleles: usize,
    pub generations: usize,
    pub initial_fittest_freq: f64,
    /// Kimura's theoretical limit `-log2(p_0^*)` bits required to isolate the fittest allele.
    pub kimura_limit_bits: f64,
    pub final_info_gain_bits: f64,
    pub final_substitutional_load_bits: f64,
    pub final_efficiency: f64,
    pub steps: Vec<SelectionTrajectoryStep>,
}

/// Parses a multi-FASTA alignment, auto-detects `DNA`, `RNA`, or `Protein`, and extracts optional
/// clade prefixes from headers formatted as `>CladeName|SeqID`.
pub fn parse_fasta_alignment(title: &str, raw_fasta: &str) -> Result<SequenceAlignment, String> {
    let mut records: Vec<(String, Option<String>, Vec<u8>)> = Vec::new();
    let mut cur_header: Option<String> = None;
    let mut cur_seq: Vec<u8> = Vec::new();

    for line in raw_fasta.lines() {
        let trimmed = line.trim();
        if trimmed.is_empty() || trimmed.starts_with(';') {
            continue;
        }
        if let Some(rest) = trimmed.strip_prefix('>') {
            if let Some(hdr) = cur_header.take() {
                let (id, clade) = split_clade_header(&hdr);
                records.push((id, clade, std::mem::take(&mut cur_seq)));
            }
            cur_header = Some(rest.trim().to_string());
        } else {
            for b in trimmed.bytes() {
                if !b.is_ascii_whitespace() {
                    cur_seq.push(b.to_ascii_uppercase());
                }
            }
        }
    }
    if let Some(hdr) = cur_header.take() {
        let (id, clade) = split_clade_header(&hdr);
        records.push((id, clade, cur_seq));
    }

    if records.len() < 2 {
        return Err(String::from(
            "Alignment requires at least 2 FASTA sequences to compute information profiles.",
        ));
    }

    let len0 = records[0].2.len();
    if len0 == 0 {
        return Err(String::from("Aligned FASTA sequences cannot be empty."));
    }
    for (id, _, seq) in &records {
        if seq.len() != len0 {
            return Err(format!(
                "Unaligned FASTA: sequence '{}' has length {}, expected {} (matching '{}').",
                id,
                seq.len(),
                len0,
                records[0].0
            ));
        }
    }

    let alphabet = detect_alphabet(&records);
    let out_records = records
        .into_iter()
        .map(|(id, clade_tag, residues)| AlignedBioSeq {
            id,
            clade_tag,
            residues,
        })
        .collect();

    Ok(SequenceAlignment {
        title: title.to_string(),
        alphabet,
        length: len0,
        records: out_records,
    })
}

fn split_clade_header(hdr: &str) -> (String, Option<String>) {
    let token = hdr.split_whitespace().next().unwrap_or("seq");
    if let Some((clade, id)) = token.split_once('|') {
        if !clade.is_empty() && !id.is_empty() {
            return (token.to_string(), Some(clade.to_string()));
        }
    }
    (token.to_string(), None)
}

fn detect_alphabet(records: &[(String, Option<String>, Vec<u8>)]) -> BioAlphabet {
    let mut total_alpha = 0usize;
    let mut nuc_count = 0usize;
    let mut has_u = false;
    let mut has_protein_only = false;

    for (_, _, seq) in records {
        for &b in seq {
            if b == b'-' || b == b'.' || b == b'?' {
                continue;
            }
            total_alpha += 1;
            match b {
                b'A' | b'C' | b'G' | b'T' => nuc_count += 1,
                b'U' => {
                    nuc_count += 1;
                    has_u = true;
                }
                b'E' | b'F' | b'I' | b'L' | b'P' | b'Q' | b'Z' | b'J' => {
                    has_protein_only = true;
                }
                _ => {}
            }
        }
    }

    if !has_protein_only && total_alpha > 0 && (nuc_count * 100 / total_alpha) >= 90 {
        if has_u {
            BioAlphabet::Rna
        } else {
            BioAlphabet::Dna
        }
    } else {
        BioAlphabet::Protein
    }
}

/// Computes Schneider & Stephens (1990) positional sequence conservation `R_seq(l)`, letter bit stacks,
/// and Adami et al. (2000) genomic physical complexity `C` across `aln`.
pub fn compute_alignment_information(
    aln: &SequenceAlignment,
    apply_small_sample_correction: bool,
    genome_bp_and_site_count: Option<(f64, f64)>,
) -> AlignmentComplexityReport {
    let syms = aln.alphabet.symbols();
    let s_card = syms.len();
    let max_bits = aln.alphabet.max_bits();
    let n_seq = aln.records.len();

    let mut sites = Vec::with_capacity(aln.length);
    let mut total_r_seq = 0.0_f64;
    let mut adami_c = 0.0_f64;

    for col in 0..aln.length {
        let mut counts = vec![0usize; s_card];
        let mut valid_n = 0usize;
        for rec in &aln.records {
            if let Some(idx) = aln.alphabet.symbol_index(rec.residues[col]) {
                counts[idx] += 1;
                valid_n += 1;
            }
        }

        let gap_frac = 1.0 - (valid_n as f64) / (n_seq.max(1) as f64);
        let occupancy = 1.0 - gap_frac;

        let mut h_bits = 0.0_f64;
        let mut best_sym = syms[0] as char;
        let mut best_cnt = 0usize;
        for (i, &cnt) in counts.iter().enumerate() {
            if cnt > best_cnt {
                best_cnt = cnt;
                best_sym = syms[i] as char;
            }
            if cnt > 0 && valid_n > 0 {
                let p = (cnt as f64) / (valid_n as f64);
                h_bits -= p * p.log2();
            }
        }

        let e_n = if valid_n > 0 {
            ((s_card - 1) as f64) / (2.0 * (valid_n as f64) * LN_2)
        } else {
            0.0
        };

        let raw_ic = if valid_n > 0 {
            (max_bits - h_bits).max(0.0) * occupancy
        } else {
            0.0
        };
        let corr_ic = if valid_n > 0 {
            let penalty = if apply_small_sample_correction {
                e_n
            } else {
                0.0
            };
            (max_bits - (h_bits + penalty)).max(0.0) * occupancy
        } else {
            0.0
        };

        total_r_seq += corr_ic;
        adami_c += raw_ic;

        let mut heights = Vec::new();
        if valid_n > 0 && corr_ic > 1e-6 {
            for (i, &cnt) in counts.iter().enumerate() {
                if cnt > 0 {
                    let f_b = (cnt as f64) / (valid_n as f64);
                    heights.push((syms[i] as char, f_b * corr_ic));
                }
            }
            heights.sort_by(|a, b| a.1.partial_cmp(&b.1).unwrap());
        }

        sites.push(SiteInformationProfile {
            position_1based: col + 1,
            non_gap_count: valid_n,
            gap_fraction: gap_frac,
            shannon_bits: h_bits,
            small_sample_error_bits: e_n,
            r_sequence_raw_bits: raw_ic,
            r_sequence_bits: corr_ic,
            consensus_symbol: if valid_n > 0 { best_sym } else { '-' },
            letter_heights_bits: heights,
        });
    }

    let capacity = (aln.length as f64) * max_bits;
    let r_freq = genome_bp_and_site_count.and_then(|(g_bp, gamma)| {
        if g_bp > 0.0 && gamma > 0.0 && g_bp >= gamma {
            Some((g_bp / gamma).log2())
        } else {
            None
        }
    });

    AlignmentComplexityReport {
        title: aln.title.clone(),
        alphabet: aln.alphabet,
        num_sequences: n_seq,
        alignment_length: aln.length,
        max_bits_per_site: max_bits,
        total_r_sequence_bits: total_r_seq,
        adami_physical_complexity_bits: adami_c,
        complexity_density: if capacity > 0.0 {
            adami_c / capacity
        } else {
            0.0
        },
        effective_conserved_sites: if max_bits > 0.0 {
            adami_c / max_bits
        } else {
            0.0
        },
        r_frequency_target_bits: r_freq,
        sites,
    }
}

/// Computes pairwise Shannon Mutual Information `MI(i, j)` and Dunn, Wahl & Gloor (2008)
/// Average Product Corrected `MIp(i, j)` across all column pairs in `aln`.
pub fn compute_epistasis_mip_matrix(
    aln: &SequenceAlignment,
    max_gap_fraction: f64,
    top_k: usize,
) -> EpistasisMatrixReport {
    let l = aln.length;
    let syms = aln.alphabet.symbols();
    let s_card = syms.len();
    let n_seq = aln.records.len().max(1);

    let mut mi_mat = vec![0.0_f64; l * l];
    let mut mir_mat = vec![0.0_f64; l * l];
    let mut joint_summaries: Vec<Vec<String>> = vec![vec![String::new(); l]; l];

    let col_valid: Vec<bool> = (0..l)
        .map(|c| {
            let valid = aln
                .records
                .iter()
                .filter(|r| aln.alphabet.symbol_index(r.residues[c]).is_some())
                .count();
            1.0 - (valid as f64) / (n_seq as f64) <= max_gap_fraction
        })
        .collect();

    for i in 0..l {
        if !col_valid[i] {
            continue;
        }
        for j in (i + 1)..l {
            if !col_valid[j] {
                continue;
            }
            let mut joint = vec![0usize; s_card * s_card];
            let mut pi = vec![0usize; s_card];
            let mut pj = vec![0usize; s_card];
            let mut n_pair = 0usize;

            for rec in &aln.records {
                if let (Some(ai), Some(bj)) = (
                    aln.alphabet.symbol_index(rec.residues[i]),
                    aln.alphabet.symbol_index(rec.residues[j]),
                ) {
                    joint[ai * s_card + bj] += 1;
                    pi[ai] += 1;
                    pj[bj] += 1;
                    n_pair += 1;
                }
            }
            if n_pair < 4 {
                continue;
            }

            let nf = n_pair as f64;
            let mut mi = 0.0_f64;
            let mut h_joint = 0.0_f64;
            let mut state_list: Vec<(usize, char, char)> = Vec::new();

            for a in 0..s_card {
                if pi[a] == 0 {
                    continue;
                }
                let pa = (pi[a] as f64) / nf;
                for b in 0..s_card {
                    let c_ab = joint[a * s_card + b];
                    if c_ab == 0 {
                        continue;
                    }
                    let pb = (pj[b] as f64) / nf;
                    let pab = (c_ab as f64) / nf;
                    mi += pab * (pab / (pa * pb)).log2();
                    h_joint -= pab * pab.log2();
                    state_list.push((c_ab, syms[a] as char, syms[b] as char));
                }
            }
            mi = mi.max(0.0);
            let mir = if h_joint > 1e-12 { mi / h_joint } else { 0.0 };

            state_list.sort_by(|x, y| y.0.cmp(&x.0));
            let summary = state_list
                .iter()
                .take(4)
                .map(|&(cnt, ca, cb)| {
                    format!("{ca}-{cb} ({:.0}%)", 100.0 * (cnt as f64) / nf)
                })
                .collect::<Vec<_>>()
                .join(", ");

            mi_mat[i * l + j] = mi;
            mi_mat[j * l + i] = mi;
            mir_mat[i * l + j] = mir;
            mir_mat[j * l + i] = mir;
            joint_summaries[i][j] = summary;
        }
    }

    let active_cols: Vec<usize> = (0..l).filter(|&c| col_valid[c]).collect();
    let m_active = active_cols.len();
    let mut col_mean_mi = vec![0.0_f64; l];
    let mut pair_sum = 0.0_f64;
    let mut pair_cnt = 0usize;

    if m_active >= 2 {
        for &i in &active_cols {
            let mut s = 0.0_f64;
            for &j in &active_cols {
                if i != j {
                    s += mi_mat[i * l + j];
                }
            }
            col_mean_mi[i] = s / ((m_active - 1) as f64);
        }
        for idx_a in 0..m_active {
            for idx_b in (idx_a + 1)..m_active {
                let i = active_cols[idx_a];
                let j = active_cols[idx_b];
                pair_sum += mi_mat[i * l + j];
                pair_cnt += 1;
            }
        }
    }

    let overall_mean_mi = if pair_cnt > 0 {
        pair_sum / (pair_cnt as f64)
    } else {
        0.0
    };

    let mut mip_mat = vec![0.0_f64; l * l];
    let mut raw_pairs: Vec<(usize, usize, f64, f64, f64, f64)> = Vec::with_capacity(pair_cnt);
    for idx_a in 0..m_active {
        for idx_b in (idx_a + 1)..m_active {
            let i = active_cols[idx_a];
            let j = active_cols[idx_b];
            let mi = mi_mat[i * l + j];
            let apc = if overall_mean_mi > 1e-12 {
                (col_mean_mi[i] * col_mean_mi[j]) / overall_mean_mi
            } else {
                0.0
            };
            let mip = mi - apc;
            mip_mat[i * l + j] = mip;
            mip_mat[j * l + i] = mip;
            raw_pairs.push((i, j, mi, apc, mip, mir_mat[i * l + j]));
        }
    }

    let mip_mean = if !raw_pairs.is_empty() {
        raw_pairs.iter().map(|p| p.4).sum::<f64>() / (raw_pairs.len() as f64)
    } else {
        0.0
    };
    let mip_var = if raw_pairs.len() >= 2 {
        raw_pairs
            .iter()
            .map(|p| (p.4 - mip_mean).powi(2))
            .sum::<f64>()
            / ((raw_pairs.len() - 1) as f64)
    } else {
        0.0
    };
    let mip_sd = mip_var.sqrt().max(1e-12);

    let mut scored: Vec<EpistasisPairScore> = raw_pairs
        .into_iter()
        .map(|(i, j, mi, apc, mip, mir)| EpistasisPairScore {
            pos_i_1based: i + 1,
            pos_j_1based: j + 1,
            mi_bits: mi,
            apc_bits: apc,
            mip_bits: mip,
            normalized_mir: mir,
            z_score: (mip - mip_mean) / mip_sd,
            top_covarying_states: joint_summaries[i][j].clone(),
        })
        .collect();

    scored.sort_by(|a, b| b.mip_bits.partial_cmp(&a.mip_bits).unwrap());
    scored.truncate(top_k.max(1));

    EpistasisMatrixReport {
        alignment_length: l,
        mean_raw_mi_bits: overall_mean_mi,
        mip_std_dev: mip_sd,
        mip_matrix: mip_mat,
        mi_matrix: mi_mat,
        top_pairs: scored,
    }
}

/// Computes per-site symmetric Jensen–Shannon divergence (`D_JS in [0, 1]` bits) and Laplace-smoothed
/// Kullback–Leibler divergence (`D_KL(A || B)`) between two clades/ecotypes in `aln`.
pub fn compute_clade_divergence(
    aln: &SequenceAlignment,
    clade_a_query: &str,
    clade_b_query: &str,
) -> Result<CladeDivergenceReport, String> {
    let (name_a, group_a, name_b, group_b) =
        partition_clades(aln, clade_a_query.trim(), clade_b_query.trim())?;

    let syms = aln.alphabet.symbols();
    let s_card = syms.len();
    let eps = 0.01_f64;

    let mut sites = Vec::with_capacity(aln.length);
    let mut total_js = 0.0_f64;
    let mut total_kl = 0.0_f64;

    for col in 0..aln.length {
        let mut ca = vec![0usize; s_card];
        let mut cb = vec![0usize; s_card];
        let mut na = 0usize;
        let mut nb = 0usize;

        for &idx in &group_a {
            if let Some(si) = aln.alphabet.symbol_index(aln.records[idx].residues[col]) {
                ca[si] += 1;
                na += 1;
            }
        }
        for &idx in &group_b {
            if let Some(si) = aln.alphabet.symbol_index(aln.records[idx].residues[col]) {
                cb[si] += 1;
                nb += 1;
            }
        }

        let cons_a = ca
            .iter()
            .enumerate()
            .max_by_key(|&(_, &c)| c)
            .map(|(i, &c)| if c > 0 { syms[i] as char } else { '-' })
            .unwrap_or('-');
        let cons_b = cb
            .iter()
            .enumerate()
            .max_by_key(|&(_, &c)| c)
            .map(|(i, &c)| if c > 0 { syms[i] as char } else { '-' })
            .unwrap_or('-');

        if na == 0 || nb == 0 {
            sites.push(SiteCladeDivergence {
                position_1based: col + 1,
                js_divergence_bits: 0.0,
                kl_p_from_q_bits: 0.0,
                clade_a_consensus: cons_a,
                clade_b_consensus: cons_b,
            });
            continue;
        }

        let mut js = 0.0_f64;
        let mut kl = 0.0_f64;
        let denom_a_sm = (na as f64) + eps * (s_card as f64);
        let denom_b_sm = (nb as f64) + eps * (s_card as f64);

        for s_i in 0..s_card {
            let pa = (ca[s_i] as f64) / (na as f64);
            let pb = (cb[s_i] as f64) / (nb as f64);
            let m = 0.5 * (pa + pb);
            if pa > 0.0 {
                js += 0.5 * pa * (pa / m).log2();
            }
            if pb > 0.0 {
                js += 0.5 * pb * (pb / m).log2();
            }

            let pa_sm = ((ca[s_i] as f64) + eps) / denom_a_sm;
            let pb_sm = ((cb[s_i] as f64) + eps) / denom_b_sm;
            kl += pa_sm * (pa_sm / pb_sm).log2();
        }

        js = js.clamp(0.0, 1.0);
        kl = kl.max(0.0);
        total_js += js;
        total_kl += kl;

        sites.push(SiteCladeDivergence {
            position_1based: col + 1,
            js_divergence_bits: js,
            kl_p_from_q_bits: kl,
            clade_a_consensus: cons_a,
            clade_b_consensus: cons_b,
        });
    }

    Ok(CladeDivergenceReport {
        clade_a_name: name_a,
        clade_b_name: name_b,
        clade_a_n: group_a.len(),
        clade_b_n: group_b.len(),
        total_js_bits: total_js,
        mean_js_bits: total_js / (aln.length.max(1) as f64),
        total_kl_bits: total_kl,
        sites,
    })
}

fn partition_clades(
    aln: &SequenceAlignment,
    q_a: &str,
    q_b: &str,
) -> Result<(String, Vec<usize>, String, Vec<usize>), String> {
    if !q_a.is_empty() && !q_b.is_empty() {
        let ga: Vec<usize> = aln
            .records
            .iter()
            .enumerate()
            .filter(|(_, r)| {
                r.id.to_ascii_lowercase()
                    .contains(&q_a.to_ascii_lowercase())
            })
            .map(|(i, _)| i)
            .collect();
        let gb: Vec<usize> = aln
            .records
            .iter()
            .enumerate()
            .filter(|(_, r)| {
                r.id.to_ascii_lowercase()
                    .contains(&q_b.to_ascii_lowercase())
            })
            .map(|(i, _)| i)
            .collect();
        if ga.is_empty() || gb.is_empty() {
            return Err(format!(
                "Clade filter matched {} ('{q_a}') and {} ('{q_b}') sequences; both need >= 1.",
                ga.len(),
                gb.len()
            ));
        }
        return Ok((q_a.to_string(), ga, q_b.to_string(), gb));
    }

    let mut tags: Vec<String> = Vec::new();
    for r in &aln.records {
        if let Some(t) = &r.clade_tag {
            if !tags.contains(t) {
                tags.push(t.clone());
            }
        }
    }
    if tags.len() >= 2 {
        let ta = tags[0].clone();
        let tb = tags[1].clone();
        let ga: Vec<usize> = aln
            .records
            .iter()
            .enumerate()
            .filter(|(_, r)| r.clade_tag.as_deref() == Some(ta.as_str()))
            .map(|(i, _)| i)
            .collect();
        let gb: Vec<usize> = aln
            .records
            .iter()
            .enumerate()
            .filter(|(_, r)| r.clade_tag.as_deref() == Some(tb.as_str()))
            .map(|(i, _)| i)
            .collect();
        return Ok((ta, ga, tb, gb));
    }

    let mid = aln.records.len() / 2;
    Ok((
        String::from("Group_1"),
        (0..mid).collect(),
        String::from("Group_2"),
        (mid..aln.records.len()).collect(),
    ))
}

/// Simulates discrete-time asexual/haploid replicator selection (`w_i = exp(dt * m_i)`) and tracks
/// Kimura (1961) genetic information gain `I_t = D_KL(p_t || p_0)` vs. cumulative substitutional load
/// `L_sub,t = sum_{tau < t} log2(W^* / <W>_{p_tau})`.
pub fn simulate_kimura_selection_information(
    initial_freqs: &[f64],
    malthusian_fitnesses: &[f64],
    dt: f64,
    generations: usize,
) -> Result<KimuraSelectionReport, String> {
    if initial_freqs.len() < 2 || initial_freqs.len() != malthusian_fitnesses.len() {
        return Err(String::from(
            "Provide >= 2 matching initial allele frequencies and Malthusian fitnesses.",
        ));
    }
    let f_sum: f64 = initial_freqs.iter().copied().filter(|&x| x > 0.0).sum();
    if f_sum <= 0.0 {
        return Err(String::from("Initial frequencies must contain positive values."));
    }

    let p0: Vec<f64> = initial_freqs.iter().map(|&x| x.max(1e-15) / f_sum).collect();
    let mut p = p0.clone();

    let m_star = malthusian_fitnesses
        .iter()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let p0_star: f64 = p0
        .iter()
        .zip(malthusian_fitnesses)
        .filter(|&(_, &m)| (m - m_star).abs() < 1e-12)
        .map(|(&freq, _)| freq)
        .sum();
    let kimura_limit = -p0_star.log2();

    let w_vec: Vec<f64> = malthusian_fitnesses.iter().map(|&m| (dt * m).exp()).collect();
    let w_star = (dt * m_star).exp();

    let mut cum_load_bits = 0.0_f64;
    let mut steps = Vec::with_capacity(generations + 1);

    for gen in 0..=generations {
        let p_star_now: f64 = p
            .iter()
            .zip(malthusian_fitnesses)
            .filter(|&(_, &m)| (m - m_star).abs() < 1e-12)
            .map(|(&freq, _)| freq)
            .sum();
        let mean_m: f64 = p
            .iter()
            .zip(malthusian_fitnesses)
            .map(|(&pi, &mi)| pi * mi)
            .sum();
        let h_bits: f64 = p
            .iter()
            .filter(|&&pi| pi > 1e-18)
            .map(|&pi| -pi * pi.log2())
            .sum();
        let kl_gain: f64 = p
            .iter()
            .zip(&p0)
            .filter(|&(&pi, _)| pi > 1e-18)
            .map(|(&pi, &p0i)| pi * (pi / p0i).log2())
            .sum::<f64>()
            .max(0.0);

        let eff = if cum_load_bits > 1e-12 {
            (kl_gain / cum_load_bits).clamp(0.0, 1.0)
        } else {
            1.0
        };

        steps.push(SelectionTrajectoryStep {
            generation: gen,
            fittest_freq: p_star_now,
            mean_malthusian_fitness: mean_m,
            population_entropy_bits: h_bits,
            kl_info_gain_bits: kl_gain,
            cumulative_substitutional_load_bits: cum_load_bits,
            selection_efficiency: eff,
        });

        if gen < generations {
            let w_bar: f64 = p.iter().zip(&w_vec).map(|(&pi, &wi)| pi * wi).sum();
            cum_load_bits += (w_star / w_bar.max(1e-30)).log2().max(0.0);
            for (pi, &wi) in p.iter_mut().zip(&w_vec) {
                *pi = (*pi * wi) / w_bar;
            }
        }
    }

    let last = steps.last().unwrap();
    Ok(KimuraSelectionReport {
        num_alleles: p0.len(),
        generations,
        initial_fittest_freq: p0_star,
        kimura_limit_bits: kimura_limit,
        final_info_gain_bits: last.kl_info_gain_bits,
        final_substitutional_load_bits: last.cumulative_substitutional_load_bits,
        final_efficiency: last.selection_efficiency,
        steps,
    })
}

/// Built-in 32-sequence × 24-nt bacterial regulatory + orthogonal `GF(4)` coevolving RNA stem-loop FASTA alignment:
/// - Columns `1..=6`: Invariant Pribnow `-10` promoter box (`TATAAT`, ~11.6 bits `R_seq`)
/// - Columns `7..=8`: Ecotype adaptive discriminator (`Thermal|` = `GC` vs. `Mesic|` = `AT`, exact `1.000` bit `D_JS`)
/// - Columns `9..=10`: Conserved `CA` spacer
/// - Columns `11..=14` paired with `22..=19`: Orthogonal `GF(4)` 4-bp compensatory RNA stem-loop (`11↔22`, `12↔21`, `13↔20`, `14↔19`)
///   where every single stem column has exact `H(l) = 2.000` bits (`0.000` bits on a single-site logo!) yet `2.000` bits of pure
///   Watson–Crick coevolution (`MIp Z > +6.0`)!
/// - Columns `15..=18`: Conserved `GAAA` apical tetraloop; Columns `23..=24`: `TT` terminator.
pub fn demo_ecoli_promoter_and_stemloop_fasta() -> &'static str {
    "\
>Thermal|Thermus_01
TATAATGCCAAAAAGAAATTTTTT
>Thermal|Thermus_02
TATAATGCCAACCGGAAACGGTTT
>Thermal|Thermus_03
TATAATGCCAAGGTGAAAACCTTT
>Thermal|Thermus_04
TATAATGCCAATTCGAAAGAATTT
>Thermal|Thermus_05
TATAATGCCACACCGAAAGGTGTT
>Thermal|Thermus_06
TATAATGCCACCATGAAAATGGTT
>Thermal|Thermus_07
TATAATGCCACGTGGAAACACGTT
>Thermal|Thermus_08
TATAATGCCACTGAGAAATCAGTT
>Thermal|Thermus_09
TATAATGCCAGAGGGAAACCTCTT
>Thermal|Thermus_10
TATAATGCCAGCTAGAAATAGCTT
>Thermal|Thermus_11
TATAATGCCAGGACGAAAGTCCTT
>Thermal|Thermus_12
TATAATGCCAGTCTGAAAAGACTT
>Thermal|Thermus_13
TATAATGCCATATTGAAAAATATT
>Thermal|Thermus_14
TATAATGCCATCGCGAAAGCGATT
>Thermal|Thermus_15
TATAATGCCATGCAGAAATGCATT
>Thermal|Thermus_16
TATAATGCCATTAGGAAACTAATT
>Mesic|Escherichia_01
TATAATATCAAAAAGAAATTTTTT
>Mesic|Escherichia_02
TATAATATCAACCGGAAACGGTTT
>Mesic|Escherichia_03
TATAATATCAAGGTGAAAACCTTT
>Mesic|Escherichia_04
TATAATATCAATTCGAAAGAATTT
>Mesic|Escherichia_05
TATAATATCACACCGAAAGGTGTT
>Mesic|Escherichia_06
TATAATATCACCATGAAAATGGTT
>Mesic|Escherichia_07
TATAATATCACGTGGAAACACGTT
>Mesic|Escherichia_08
TATAATATCACTGAGAAATCAGTT
>Mesic|Escherichia_09
TATAATATCAGAGGGAAACCTCTT
>Mesic|Escherichia_10
TATAATATCAGCTAGAAATAGCTT
>Mesic|Escherichia_11
TATAATATCAGGACGAAAGTCCTT
>Mesic|Escherichia_12
TATAATATCAGTCTGAAAAGACTT
>Mesic|Escherichia_13
TATAATATCATATTGAAAAATATT
>Mesic|Escherichia_14
TATAATATCATCGCGAAAGCGATT
>Mesic|Escherichia_15
TATAATATCATGCAGAAATGCATT
>Mesic|Escherichia_16
TATAATATCATTAGGAAACTAATT
"
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn recovers_pribnow_logo_stemloop_epistasis_clade_jsd_and_kimura_bound() {
        let aln =
            parse_fasta_alignment("Demo_ITME", demo_ecoli_promoter_and_stemloop_fasta()).unwrap();
        assert_eq!(aln.alphabet, BioAlphabet::Dna);
        assert_eq!(aln.records.len(), 32);
        assert_eq!(aln.length, 24);

        let info = compute_alignment_information(&aln, true, Some((4_600_000.0, 4_500.0)));
        assert!(info.sites[0].r_sequence_bits > 1.90 && info.sites[0].r_sequence_bits <= 2.0);
        assert!(info.adami_physical_complexity_bits > 20.0);

        // Epistasis MIp must rank the 4 compensatory RNA stem-loop base-pairs (11-22, 12-21, 13-20, 14-19) #1..#4
        let epi = compute_epistasis_mip_matrix(&aln, 0.2, 8);
        for top in epi.top_pairs.iter().take(4) {
            assert!(
                matches!(
                    (top.pos_i_1based, top.pos_j_1based),
                    (11, 22) | (12, 21) | (13, 20) | (14, 19)
                ),
                "Unexpected top epistatic pair: {}-{}",
                top.pos_i_1based,
                top.pos_j_1based
            );
            assert!(top.z_score > 4.0);
        }

        // Clade JSD between Thermal and Mesic must equal 1.000 bit at adaptive discriminator site 7 ('G' vs 'A')
        let div = compute_clade_divergence(&aln, "", "").unwrap();
        assert!((div.sites[6].js_divergence_bits - 1.0).abs() < 1e-6);

        // Kimura (1961) selection trajectory: p0* = 0.01 -> limit = -log2(0.01) = 6.643856 bits, I_t <= L_sub,t
        let sim = simulate_kimura_selection_information(
            &[0.01, 0.39, 0.60],
            &[0.15, 0.05, 0.0],
            1.0,
            250,
        )
        .unwrap();
        for st in &sim.steps {
            assert!(st.kl_info_gain_bits <= st.cumulative_substitutional_load_bits + 1e-9);
        }
        assert!((sim.final_info_gain_bits - sim.kimura_limit_bits).abs() < 0.01);
        assert!(sim.final_efficiency > 0.99);
    }
}
