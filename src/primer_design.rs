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

//! De novo PCR primer picker, hydrolysis (`TaqMan`) qPCR probe selector,
//! multiplex primer-pool conflict solver, and zero-index FASTA off-target amplicon scanner.

use crate::gibson::{scan_circular_bindings, Strand};
use crate::structure::{evaluate_oligo_qc, score_dimer, OligoQcReport};
use crate::thermo::{recommend_annealing_temp, PolymerasePreset, ReactionConditions};
use crate::{normalize_dna, reverse_complement, OligoError};

/// Configuration constraints for de novo PCR / qPCR primer pair selection.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PrimerDesignConfig {
    /// Primer length search window in nucleotides (default `18..=25`).
    pub primer_len_range: (usize, usize),
    /// Desired PCR product length bounds in base pairs (default `120..=600`).
    pub amplicon_len_range: (usize, usize),
    /// Target primer melting temperature `Tm` (`C`, default `60.0`).
    pub target_tm_c: f64,
    /// Acceptable `Tm` tolerance (`+/- C` around `target_tm_c`, default `3.0`).
    pub tm_tolerance_c: f64,
    /// Maximum allowed `|Tm_fwd - Tm_rev|` (`C`, default `2.0`).
    pub max_delta_tm_c: f64,
    /// Optional 0-based `(start, end)` target region that every amplicon must span.
    pub target_region: Option<(usize, usize)>,
    /// Design an internal hydrolysis (`TaqMan`) probe (`Tm = target_tm + 9 C`, no 5'-G quencher).
    pub include_qpcr_probe: bool,
    /// Polymerase buffer preset governing ion concentrations and `Ta`.
    pub polymerase: PolymerasePreset,
    /// Number of top-ranked candidate pairs to return (default `5`).
    pub num_return: usize,
}

impl Default for PrimerDesignConfig {
    fn default() -> Self {
        Self {
            primer_len_range: (18, 25),
            amplicon_len_range: (120, 600),
            target_tm_c: 60.0,
            tm_tolerance_c: 3.0,
            max_delta_tm_c: 2.5,
            target_region: None,
            include_qpcr_probe: false,
            polymerase: PolymerasePreset::Q5,
            num_return: 5,
        }
    }
}

/// Internal qPCR hydrolysis probe (`5' -> 3'`) avoiding 5'-terminal guanosine quenching.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct QpcrProbe {
    pub sequence: String,
    pub strand: Strand,
    pub start_0: usize,
    pub end_0: usize,
    pub tm_c: f64,
    pub gc_pct: f64,
    pub penalty_score: f64,
}

/// Potential unintended PCR product discovered on a background FASTA sequence set.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OffTargetAmplicon {
    pub record_name: String,
    pub start_0: usize,
    pub end_0: usize,
    pub product_len_bp: usize,
    pub fwd_mismatches: usize,
    pub rev_mismatches: usize,
}

/// Ranked forward + reverse primer pair (with optional qPCR probe and off-target audit).
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PrimerPairCandidate {
    pub rank: usize,
    pub fwd: OligoQcReport,
    pub fwd_start_0: usize,
    pub rev: OligoQcReport,
    /// 0-based inclusive top-strand coordinate of the 3'-most base of the amplicon (`rev` 5' anchor).
    pub rev_end_0: usize,
    pub amplicon_len_bp: usize,
    pub amplicon_gc_pct: f64,
    pub delta_tm_c: f64,
    pub recommended_ta_c: f64,
    pub cross_dimer_3p_dg_kcal: f64,
    pub qpcr_probe: Option<QpcrProbe>,
    pub off_target_count: usize,
    pub composite_penalty: f64,
}

/// Locus target entering a multiplex PCR design run.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiplexTarget {
    pub locus_name: String,
    pub template_sequence: String,
    pub target_region: Option<(usize, usize)>,
}

/// Conflict-checked multiplex primer pool with pairwise cross-dimer matrix.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MultiplexPoolPlan {
    pub selected_pairs: Vec<(String, PrimerPairCandidate)>,
    pub oligo_names: Vec<String>,
    /// Symmetric `2N x 2N` matrix of worst 3'-extensible dimer `Delta G` (`kcal/mol`) at 37 C.
    pub dimer_dg_matrix: Vec<Vec<f64>>,
    pub worst_pool_dimer_dg_kcal: f64,
    pub pooled_ta_c: f64,
}

impl MultiplexPoolPlan {
    /// Exports the entire multiplex primer set as a vendor-ready bulk-order `.tsv`.
    pub fn to_order_tsv(&self) -> String {
        let mut out =
            String::from("Name\tSequence\tLength_nt\tTm_C\tGC_Pct\tExpected_Amplicon_bp\n");
        for (locus, pair) in &self.selected_pairs {
            out.push_str(&format!(
                "{}_F\t{}\t{}\t{:.1}\t{:.1}\t{}\n",
                locus,
                pair.fwd.sequence,
                pair.fwd.length,
                pair.fwd.tm_c,
                pair.fwd.gc_pct,
                pair.amplicon_len_bp
            ));
            out.push_str(&format!(
                "{}_R\t{}\t{}\t{:.1}\t{:.1}\t{}\n",
                locus,
                pair.rev.sequence,
                pair.rev.length,
                pair.rev.tm_c,
                pair.rev.gc_pct,
                pair.amplicon_len_bp
            ));
            if let Some(pr) = &pair.qpcr_probe {
                out.push_str(&format!(
                    "{}_PRB\t{}\t{}\t{:.1}\t{:.1}\t{}\n",
                    locus,
                    pr.sequence,
                    pr.sequence.len(),
                    pr.tm_c,
                    pr.gc_pct,
                    pair.amplicon_len_bp
                ));
            }
        }
        out
    }
}

/// Parses a multi-FASTA text buffer into `(record_id, normalized_sequence)` pairs.
/// Also accepts raw unheaded DNA strings as a single `"Template_1"` record.
pub fn parse_fasta_or_raw(input: &str) -> Result<Vec<(String, String)>, OligoError> {
    let trimmed = input.trim_start_matches('\u{feff}').trim();
    if trimmed.is_empty() {
        return Err(OligoError::EmptySequence);
    }
    if !trimmed.starts_with('>') {
        return Ok(vec![("Template_1".to_string(), normalize_dna(trimmed)?)]);
    }

    let mut records = Vec::new();
    let mut cur_id = String::new();
    let mut cur_seq = String::new();

    for line in trimmed.lines() {
        let l = line.trim();
        if let Some(header) = l.strip_prefix('>') {
            if !cur_id.is_empty() {
                records.push((cur_id.clone(), normalize_dna(&cur_seq)?));
                cur_seq.clear();
            }
            cur_id = header
                .split_whitespace()
                .next()
                .unwrap_or("Seq")
                .to_string();
        } else if !l.is_empty() {
            cur_seq.push_str(l);
        }
    }
    if !cur_id.is_empty() && !cur_seq.is_empty() {
        records.push((cur_id, normalize_dna(&cur_seq)?));
    }
    if records.is_empty() {
        return Err(OligoError::EmptySequence);
    }
    Ok(records)
}

/// Scans a list of background `(record_id, sequence)` templates for unintended PCR amplicons
/// up to `max_product_bp` using up to `max_mismatches` per primer (with exact 3'-OH match).
pub fn scan_offtarget_amplicons(
    fwd_seq: &str,
    rev_seq: &str,
    background_records: &[(String, String)],
    max_product_bp: usize,
    max_mismatches: usize,
    cond: &ReactionConditions,
) -> Result<Vec<OffTargetAmplicon>, OligoError> {
    let oligos = vec![
        ("FWD".to_string(), fwd_seq.to_string()),
        ("REV".to_string(), rev_seq.to_string()),
    ];
    let min_core = fwd_seq.len().min(rev_seq.len()).saturating_sub(2).max(14);
    let mut products = Vec::new();

    for (rec_name, rec_seq) in background_records {
        if rec_seq.len() < min_core {
            continue;
        }
        let hits =
            scan_circular_bindings(rec_seq, false, &oligos, min_core, max_mismatches, cond)?;
        let up_hits: Vec<_> = hits.iter().filter(|h| h.strand == Strand::Forward).collect();
        let dn_hits: Vec<_> = hits.iter().filter(|h| h.strand == Strand::Reverse).collect();

        for u in &up_hits {
            for d in &dn_hits {
                // On a linear top strand, `u.bind_start_0` is 5' start and `d.bind_start_0` is 3' end of product
                if d.bind_start_0 > u.bind_start_0 {
                    let prod_len = d.bind_start_0 - u.bind_start_0 + 1;
                    if prod_len <= max_product_bp {
                        products.push(OffTargetAmplicon {
                            record_name: rec_name.clone(),
                            start_0: u.bind_start_0,
                            end_0: d.bind_start_0,
                            product_len_bp: prod_len,
                            fwd_mismatches: u.mismatches,
                            rev_mismatches: d.mismatches,
                        });
                    }
                }
            }
        }
    }

    Ok(products)
}

struct PrecomputedProbeCand {
    start_0: usize,
    end_excl: usize,
    strand: Strand,
    qc: OligoQcReport,
    base_pen: f64,
}

fn precompute_qpcr_probe_candidates(
    template: &str,
    search_start: usize,
    search_end: usize,
    target_probe_tm_c: f64,
    cond: &ReactionConditions,
) -> Vec<PrecomputedProbeCand> {
    if search_end <= search_start || (search_end - search_start) < 20 {
        return Vec::new();
    }
    let mut cands = Vec::new();
    let max_pos = search_end.saturating_sub(19);
    for pos in search_start..max_pos {
        for plen in 20..=30 {
            if pos + plen > search_end {
                break;
            }
            let slice = &template[pos..(pos + plen)];
            for (strand, cand) in [
                (Strand::Forward, slice.to_string()),
                (Strand::Reverse, reverse_complement(slice)),
            ] {
                // TaqMan rule 1: Never start 5' with G (quenches reporter dye even after cleavage)
                if cand.starts_with('G') {
                    continue;
                }
                // TaqMan rule 2: Avoid >=4 consecutive G bases
                if cand.contains("GGGG") {
                    continue;
                }
                let Ok(qc) = evaluate_oligo_qc(&cand, cond, target_probe_tm_c - 8.0) else {
                    continue;
                };
                let tm_diff = (qc.tm_c - target_probe_tm_c).abs();
                if tm_diff > 4.0 || qc.has_g4_risk {
                    continue;
                }
                let base_pen = tm_diff * 2.5 + qc.penalty_score * 0.5;
                cands.push(PrecomputedProbeCand {
                    start_0: pos,
                    end_excl: pos + plen,
                    strand,
                    qc,
                    base_pen,
                });
            }
        }
    }
    cands.sort_by(|a, b| {
        a.base_pen
            .partial_cmp(&b.base_pen)
            .unwrap_or(core::cmp::Ordering::Equal)
    });
    cands
}

fn pick_internal_qpcr_probe_from_pool(
    pool: &[PrecomputedProbeCand],
    inner_start: usize,
    inner_end: usize,
    fwd_seq: &str,
    rev_seq: &str,
    cond: &ReactionConditions,
) -> Result<Option<QpcrProbe>, OligoError> {
    if inner_end <= inner_start || (inner_end - inner_start) < 20 {
        return Ok(None);
    }
    let mut best: Option<QpcrProbe> = None;
    let mut best_raw_pen = f64::INFINITY;

    for pc in pool {
        if pc.base_pen >= best_raw_pen {
            break;
        }
        if pc.start_0 < inner_start || pc.end_excl > inner_end {
            continue;
        }
        let (_, df) = score_dimer(&pc.qc.sequence, fwd_seq, cond, 37.0)?;
        let (_, dr) = score_dimer(&pc.qc.sequence, rev_seq, cond, 37.0)?;
        let cross_pen = df.map_or(0.0, |d| (-d.delta_g_kcal - 4.0).max(0.0))
            + dr.map_or(0.0, |d| (-d.delta_g_kcal - 4.0).max(0.0));

        let pen = pc.base_pen + cross_pen * 3.0;
        if pen < best_raw_pen {
            best_raw_pen = pen;
            best = Some(QpcrProbe {
                sequence: pc.qc.sequence.clone(),
                strand: pc.strand,
                start_0: pc.start_0,
                end_0: pc.end_excl - 1,
                tm_c: pc.qc.tm_c,
                gc_pct: pc.qc.gc_pct,
                penalty_score: (pen * 100.0).round() / 100.0,
            });
        }
    }

    Ok(best)
}

/// Designs and ranks PCR (or qPCR) primer pairs across `raw_template`, optionally screening
/// against `background_fasta` records to penalize off-target amplicons.
pub fn design_pcr_primers(
    raw_template: &str,
    cfg: &PrimerDesignConfig,
    background_records: Option<&[(String, String)]>,
) -> Result<Vec<PrimerPairCandidate>, OligoError> {
    let template = normalize_dna(raw_template)?;
    let n = template.len();
    let (min_amp, max_amp) = cfg.amplicon_len_range;
    let (min_len, max_len) = cfg.primer_len_range;

    if n < min_amp || min_len < 12 || min_len > max_len {
        return Err(OligoError::SequenceTooShort {
            len: n,
            min: min_amp.max(12),
        });
    }

    let cond = ReactionConditions::for_preset(cfg.polymerase);
    let (roi_start, roi_end) = cfg.target_region.unwrap_or((0, n));
    let fwd_max_start = if cfg.target_region.is_some() {
        roi_start.min(n.saturating_sub(min_amp))
    } else {
        n.saturating_sub(min_amp)
    };
    let fwd_min_start = roi_start.saturating_sub(max_amp);

    // Collect candidate forward oligos `(start_0, OligoQcReport, score)`
    let mut fwd_cands = Vec::new();
    for pos in fwd_min_start..=fwd_max_start {
        for plen in min_len..=max_len {
            if pos + plen > n {
                break;
            }
            let sub = &template[pos..(pos + plen)];
            // Prefer 3' G/C clamp (1-2 G/C in last 2 bases, but not >=4 G/C in last 5)
            let Ok(qc) = evaluate_oligo_qc(sub, &cond, cfg.target_tm_c) else {
                continue;
            };
            let tm_err = (qc.tm_c - cfg.target_tm_c).abs();
            if tm_err > cfg.tm_tolerance_c || qc.has_g4_risk || qc.max_homopolymer >= 5 {
                continue;
            }
            let clamp_bonus = if matches!(sub.as_bytes()[plen - 1], b'G' | b'C') {
                -0.5
            } else {
                0.4
            };
            let score = tm_err * 2.5 + qc.penalty_score + clamp_bonus;
            fwd_cands.push((pos, qc, score));
        }
    }
    fwd_cands.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(core::cmp::Ordering::Equal));
    fwd_cands.truncate(60);

    // Collect candidate reverse oligos `(rev_end_0, OligoQcReport, score)`
    let rev_min_end = if cfg.target_region.is_some() {
        roi_end.max(min_amp - 1).min(n - 1)
    } else {
        min_amp - 1
    };
    let rev_max_end = (roi_end + max_amp).min(n - 1);

    let mut rev_cands = Vec::new();
    for end_0 in rev_min_end..=rev_max_end {
        for plen in min_len..=max_len {
            if end_0 + 1 < plen {
                continue;
            }
            let start_0 = end_0 + 1 - plen;
            let rc = reverse_complement(&template[start_0..=end_0]);
            let Ok(qc) = evaluate_oligo_qc(&rc, &cond, cfg.target_tm_c) else {
                continue;
            };
            let tm_err = (qc.tm_c - cfg.target_tm_c).abs();
            if tm_err > cfg.tm_tolerance_c || qc.has_g4_risk || qc.max_homopolymer >= 5 {
                continue;
            }
            let clamp_bonus = if matches!(rc.as_bytes()[plen - 1], b'G' | b'C') {
                -0.5
            } else {
                0.4
            };
            let score = tm_err * 2.5 + qc.penalty_score + clamp_bonus;
            rev_cands.push((end_0, qc, score));
        }
    }
    rev_cands.sort_by(|a, b| a.2.partial_cmp(&b.2).unwrap_or(core::cmp::Ordering::Equal));
    rev_cands.truncate(60);

    let probe_pool = if cfg.include_qpcr_probe && !fwd_cands.is_empty() && !rev_cands.is_empty() {
        let min_inner_s = fwd_cands
            .iter()
            .map(|(f_start, f_qc, _)| f_start + f_qc.length + 1)
            .min()
            .unwrap_or(0);
        let max_inner_e = rev_cands
            .iter()
            .map(|(r_end, r_qc, _)| (r_end + 1).saturating_sub(r_qc.length))
            .max()
            .unwrap_or(0);
        precompute_qpcr_probe_candidates(
            &template,
            min_inner_s,
            max_inner_e,
            cfg.target_tm_c + 9.0,
            &cond,
        )
    } else {
        Vec::new()
    };

    let mut pairs = Vec::new();
    for (f_start, f_qc, f_score) in &fwd_cands {
        for (r_end, r_qc, r_score) in &rev_cands {
            if *r_end <= f_start + f_qc.length {
                continue;
            }
            let amp_len = r_end - f_start + 1;
            if !(min_amp..=max_amp).contains(&amp_len) {
                continue;
            }
            if let Some((rs, re)) = cfg.target_region {
                if *f_start > rs || *r_end + 1 < re {
                    continue;
                }
            }

            let delta_tm = (f_qc.tm_c - r_qc.tm_c).abs();
            if delta_tm > cfg.max_delta_tm_c {
                continue;
            }

            let (_, cross_3p) = score_dimer(&f_qc.sequence, &r_qc.sequence, &cond, 37.0)?;
            let cross_dg = cross_3p.map(|d| d.delta_g_kcal).unwrap_or(0.0);
            if cross_dg < -6.5 {
                continue;
            }

            let probe = if cfg.include_qpcr_probe {
                let inner_s = f_start + f_qc.length + 1;
                let inner_e = r_end + 1 - r_qc.length;
                let prb = pick_internal_qpcr_probe_from_pool(
                    &probe_pool,
                    inner_s,
                    inner_e,
                    &f_qc.sequence,
                    &r_qc.sequence,
                    &cond,
                )?;
                if prb.is_none() {
                    continue;
                }
                prb
            } else {
                None
            };

            let off_targets = if let Some(bg) = background_records {
                let hits = scan_offtarget_amplicons(
                    &f_qc.sequence,
                    &r_qc.sequence,
                    bg,
                    max_amp * 3,
                    1,
                    &cond,
                )?;
                // Subtract 1 if the primary template itself is inside `bg`
                hits.len().saturating_sub(1)
            } else {
                0
            };

            let ta_rec =
                recommend_annealing_temp(cfg.polymerase, &f_qc.sequence, &r_qc.sequence)?;
            let amp_slice = &template[*f_start..=*r_end];
            let amp_gc = crate::gc_fraction(amp_slice) * 100.0;

            let cross_pen = if cross_dg < -3.5 {
                (-cross_dg - 3.5) * 3.5
            } else {
                0.0
            };
            let probe_pen = probe.as_ref().map_or(0.0, |p| p.penalty_score * 0.4);
            let total_pen = f_score
                + r_score
                + delta_tm * 3.0
                + cross_pen
                + probe_pen
                + (off_targets as f64) * 25.0;

            pairs.push(PrimerPairCandidate {
                rank: 0,
                fwd: f_qc.clone(),
                fwd_start_0: *f_start,
                rev: r_qc.clone(),
                rev_end_0: *r_end,
                amplicon_len_bp: amp_len,
                amplicon_gc_pct: (amp_gc * 10.0).round() / 10.0,
                delta_tm_c: (delta_tm * 100.0).round() / 100.0,
                recommended_ta_c: ta_rec.recommended_ta_c,
                cross_dimer_3p_dg_kcal: cross_dg,
                qpcr_probe: probe,
                off_target_count: off_targets,
                composite_penalty: (total_pen * 100.0).round() / 100.0,
            });
        }
    }

    // Deduplicate overlapping windows (>=80% identical coordinates) so top N pairs are diverse
    pairs.sort_by(|a, b| {
        a.composite_penalty
            .partial_cmp(&b.composite_penalty)
            .unwrap_or(core::cmp::Ordering::Equal)
    });

    let mut diverse: Vec<PrimerPairCandidate> = Vec::new();
    for mut cand in pairs {
        let too_close = diverse.iter().any(|existing| {
            existing.fwd_start_0.abs_diff(cand.fwd_start_0) < 6
                && existing.rev_end_0.abs_diff(cand.rev_end_0) < 6
        });
        if !too_close {
            cand.rank = diverse.len() + 1;
            diverse.push(cand);
            if diverse.len() >= cfg.num_return.max(1) {
                break;
            }
        }
    }

    Ok(diverse)
}

/// Designs a compatible multiplex PCR pool across multiple target loci, enforcing both
/// cross-primer 3' dimer `Delta G >= dg_3p_floor_kcal` (e.g., `-5.5 kcal/mol`) and optional
/// minimum amplicon-length spacing (`min_size_spacing_bp`) for single-lane gel resolution.
pub fn design_multiplex_pool(
    targets: &[MultiplexTarget],
    base_cfg: &PrimerDesignConfig,
    dg_3p_floor_kcal: f64,
    min_size_spacing_bp: usize,
) -> Result<MultiplexPoolPlan, OligoError> {
    if targets.is_empty() {
        return Err(OligoError::EmptySequence);
    }

    let cond = ReactionConditions::for_preset(base_cfg.polymerase);
    let mut per_locus_cands: Vec<(String, Vec<PrimerPairCandidate>)> =
        Vec::with_capacity(targets.len());

    for t in targets {
        let mut local_cfg = base_cfg.clone();
        local_cfg.target_region = t.target_region;
        local_cfg.num_return = 10;
        let cands = design_pcr_primers(&t.template_sequence, &local_cfg, None)?;
        if cands.is_empty() {
            return Err(OligoError::AssemblyDesignFailed {
                junction_index: 0,
                reason: format!("No valid primer candidates found for locus '{}'.", t.locus_name),
            });
        }
        per_locus_cands.push((t.locus_name.clone(), cands));
    }

    // Greedy conflict-minimizing selection across loci
    let mut chosen: Vec<(String, PrimerPairCandidate)> = Vec::with_capacity(targets.len());
    for (locus_name, cands) in per_locus_cands {
        let mut best_choice: Option<(f64, PrimerPairCandidate)> = None;

        for cand in cands {
            let mut pool_pen = cand.composite_penalty;
            for (_, prev) in &chosen {
                if min_size_spacing_bp > 0
                    && cand.amplicon_len_bp.abs_diff(prev.amplicon_len_bp) < min_size_spacing_bp
                {
                    pool_pen += 50.0;
                }
                for s1 in [&cand.fwd.sequence, &cand.rev.sequence] {
                    for s2 in [&prev.fwd.sequence, &prev.rev.sequence] {
                        let (_, d3) = score_dimer(s1, s2, &cond, 37.0)?;
                        if let Some(hit) = d3 {
                            if hit.delta_g_kcal < dg_3p_floor_kcal {
                                pool_pen += (dg_3p_floor_kcal - hit.delta_g_kcal) * 20.0;
                            }
                        }
                    }
                }
            }
            if best_choice
                .as_ref()
                .map_or(true, |(b_pen, _)| pool_pen < *b_pen)
            {
                best_choice = Some((pool_pen, cand));
            }
        }
        let (_, selected_pair) = best_choice.ok_or_else(|| OligoError::AssemblyDesignFailed {
            junction_index: chosen.len(),
            reason: format!("No valid primer pair found for locus '{locus_name}'."),
        })?;
        chosen.push((locus_name, selected_pair));
    }

    // Build full 2N x 2N 3'-extensible dimer matrix for terminal/GUI heatmap inspection
    let mut oligo_names = Vec::with_capacity(chosen.len() * 2);
    let mut oligo_seqs = Vec::with_capacity(chosen.len() * 2);
    let mut min_ta = f64::INFINITY;
    for (locus, pair) in &chosen {
        oligo_names.push(format!("{locus}_F"));
        oligo_seqs.push(pair.fwd.sequence.clone());
        oligo_names.push(format!("{locus}_R"));
        oligo_seqs.push(pair.rev.sequence.clone());
        min_ta = min_ta.min(pair.recommended_ta_c);
    }

    let m = oligo_seqs.len();
    let mut matrix = vec![vec![0.0_f64; m]; m];
    let mut worst_dg = 0.0_f64;
    for i in 0..m {
        for j in i..m {
            let (_, d3) = score_dimer(&oligo_seqs[i], &oligo_seqs[j], &cond, 37.0)?;
            let val = d3.map(|d| d.delta_g_kcal).unwrap_or(0.0);
            matrix[i][j] = val;
            matrix[j][i] = val;
            worst_dg = worst_dg.min(val);
        }
    }

    Ok(MultiplexPoolPlan {
        selected_pairs: chosen,
        oligo_names,
        dimer_dg_matrix: matrix,
        worst_pool_dimer_dg_kcal: worst_dg,
        pooled_ta_c: min_ta,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn designs_ranked_qpcr_primers_and_internal_taqman_probe() {
        let template = "GCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGC\
                        ATGCGATCGCTAGCTAGCTGATCGATCGATCGCTAGCTAGCGATCGATCGTAGCTAGCTGATCGATCGATCGT\
                        CGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCG\
                        TAGCTAGCTGATCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCT";
        let mut cfg = PrimerDesignConfig::default();
        cfg.amplicon_len_range = (100, 240);
        cfg.include_qpcr_probe = true;

        let pairs = design_pcr_primers(template, &cfg, None).unwrap();
        assert!(!pairs.is_empty());
        let top = &pairs[0];
        assert!(top.amplicon_len_bp >= 100 && top.amplicon_len_bp <= 240);
        let probe = top.qpcr_probe.as_ref().expect("Expected internal TaqMan probe");
        assert!(!probe.sequence.starts_with('G'));
        assert!(probe.tm_c > top.fwd.tm_c);

        // UTF-8 BOM-prefixed FASTA and raw DNA parse cleanly
        let bom_fa = parse_fasta_or_raw("\u{feff}>BomSeq\nACGTACGTACGT\n").unwrap();
        assert_eq!(bom_fa, vec![(String::from("BomSeq"), String::from("ACGTACGTACGT"))]);
        let bom_raw = parse_fasta_or_raw("\u{feff}acgtacgt").unwrap();
        assert_eq!(bom_raw, vec![(String::from("Template_1"), String::from("ACGTACGT"))]);
    }
}
