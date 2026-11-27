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

//! Circular plasmid primer-binding discovery and two-stage isothermal Gibson / HiFi
//! multi-fragment assembly design shared by `bench-tui` and `open-plasmid`.

use crate::structure::{evaluate_oligo_qc, fold_hairpin, score_dimer};
use crate::thermo::{
    calc_aligned_duplex, calc_perfect_duplex, recommend_annealing_temp, PolymerasePreset,
    ReactionConditions, SaltCorrectionMethod,
};
use crate::{complement_base, normalize_dna, reverse_complement, OligoError};

/// Template strand orientation (`Forward` = top `5'->3'`, `Reverse` = bottom reverse-complement).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum Strand {
    Forward,
    Reverse,
}

/// Primer binding site discovered on a linear or circular plasmid sequence, supporting
/// origin-wrapping (`N -> 0`) and 5' cloning overhang tails.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CircularBindingSite {
    pub oligo_name: String,
    pub strand: Strand,
    /// 0-based plasmid index of the 5'-most template-bound base of the annealing core.
    pub bind_start_0: usize,
    /// 0-based plasmid index of the 3'-OH terminal base of the primer on the plasmid.
    pub bind_end_0: usize,
    /// Number of template-annealed nucleotides in the 3' core.
    pub bound_len: usize,
    /// Unbound 5' tail sequence (e.g., Gibson homology arm or restriction site pad).
    pub overhang_5p_seq: String,
    /// Template-annealing 3' core sequence (`5' -> 3'`).
    pub core_3p_seq: String,
    /// Number of internal mismatches inside the annealed core.
    pub mismatches: usize,
    /// Salt-corrected `Tm` (`C`) of the 3' template-annealing core (Cycle 1–2 annealing).
    pub core_tm_c: f64,
    /// Salt-corrected `Tm` (`C`) of the full-length primer once incorporated into amplicon.
    pub full_tm_c: f64,
}

/// Named dsDNA fragment (`5' -> 3'` top strand) entering a Gibson / HiFi assembly.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GibsonFragmentInput {
    pub name: String,
    pub sequence: String,
}

/// Configuration rules governing two-stage Gibson / NEBuilder HiFi primer & junction design.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GibsonConfig {
    /// Target `Tm` (`C`) for the 3' template-binding core during PCR amplification.
    pub target_core_tm_c: f64,
    /// Allowed length bounds (`nt`) for the 3' template-binding core.
    pub core_len_range: (usize, usize),
    /// Allowed length bounds (`bp`) for the spliced isothermal dsDNA homology overlap.
    pub overlap_len_range: (usize, usize),
    /// Minimum duplex `Tm` (`C`) required for the overlap at 50 C isothermal incubation.
    pub min_overlap_tm_c: f64,
    /// Maximum (least negative) `Delta G` (`kcal/mol`) at 50 C required for stable overlap annealing.
    pub max_overlap_dg_50c_kcal: f64,
    /// Polymerase preset used to amplify the constituent fragments.
    pub polymerase: PolymerasePreset,
}

impl Default for GibsonConfig {
    fn default() -> Self {
        Self {
            target_core_tm_c: 60.0,
            core_len_range: (18, 32),
            overlap_len_range: (20, 38),
            min_overlap_tm_c: 50.0,
            max_overlap_dg_50c_kcal: -18.0,
            polymerase: PolymerasePreset::Q5,
        }
    }
}

/// Chimeric Gibson PCR primer split into its 5' homology tail and 3' template-binding core.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChimericPrimer {
    pub name: String,
    pub tail_5p_lower: String,
    pub core_3p_upper: String,
    /// Standard cloning convention: lowercase 5' homology tail + UPPERCASE 3' annealing core.
    pub annotated_sequence: String,
    pub full_sequence_upper: String,
    pub length: usize,
    pub core_tm_c: f64,
    pub full_tm_c: f64,
    pub gc_pct: f64,
    pub hairpin_dg_kcal: f64,
    pub self_dimer_3p_dg_kcal: f64,
    pub g4_score: f64,
}

/// Thermodynamic and structural specification for a single fragment-to-fragment junction.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GibsonJunctionPlan {
    pub junction_index: usize,
    pub upstream_fragment: String,
    pub downstream_fragment: String,
    pub overlap_sequence: String,
    pub overlap_len: usize,
    pub overlap_tm_c: f64,
    pub overlap_dg_50c_kcal: f64,
    /// Worst secondary-structure hairpin `Delta G` (`kcal/mol`) at 50 C on T5-exonuclease-exposed 3' ssDNA overhangs.
    pub ssdna_hairpin_dg_50c_kcal: f64,
    pub upstream_rev_primer: ChimericPrimer,
    pub downstream_fwd_primer: ChimericPrimer,
    pub cross_dimer_3p_dg_kcal: f64,
    pub warnings: Vec<String>,
}

/// Per-fragment PCR setup parameters paired with each Gibson amplicon.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FragmentPcrProtocol {
    pub fragment_name: String,
    pub fwd_primer_name: String,
    pub rev_primer_name: String,
    pub amplicon_len_bp: usize,
    pub recommended_ta_c: f64,
    pub extension_sec: u32,
}

/// Complete multi-fragment Gibson / HiFi assembly design output.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GibsonAssemblyPlan {
    pub circular: bool,
    pub assembled_sequence: String,
    pub assembled_length_bp: usize,
    pub junctions: Vec<GibsonJunctionPlan>,
    pub fragment_pcrs: Vec<FragmentPcrProtocol>,
}

impl GibsonAssemblyPlan {
    /// Exports all assembly primers in bulk-order `.tsv` format compatible with IDT/Twist uploads
    /// and `bench-tui` / `open-plasmid` `.crv-oligo.tsv` stock tables.
    pub fn to_order_tsv(&self) -> String {
        let mut out = String::from(
            "Name\tSequence\tAnnotated_5pTail_3pCore\tLength_nt\tCore_Tm_C\tFull_Tm_C\tGC_Pct\n",
        );
        for j in &self.junctions {
            for p in [&j.downstream_fwd_primer, &j.upstream_rev_primer] {
                out.push_str(&format!(
                    "{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\n",
                    p.name,
                    p.full_sequence_upper,
                    p.annotated_sequence,
                    p.length,
                    p.core_tm_c,
                    p.full_tm_c,
                    p.gc_pct
                ));
            }
        }
        out
    }
}

/// Scans a circular (or linear) plasmid sequence for primer binding sites (allowing 5' overhang
/// tails and up to `max_internal_mismatches` inside the annealed 3' core, requiring at least `min_core_len`
/// bound bases and an exact 3'-terminal base match for polymerase extension).
pub fn scan_circular_bindings(
    raw_plasmid_seq: &str,
    circular: bool,
    oligos: &[(String, String)],
    min_core_len: usize,
    max_internal_mismatches: usize,
    cond: &ReactionConditions,
) -> Result<Vec<CircularBindingSite>, OligoError> {
    let plasmid = normalize_dna(raw_plasmid_seq)?;
    let n = plasmid.len();
    let min_core = min_core_len.max(12);
    if n < min_core {
        return Err(OligoError::SequenceTooShort {
            len: n,
            min: min_core,
        });
    }

    let rc_plasmid = reverse_complement(&plasmid);
    let mut hits = Vec::new();

    for (name, raw_oligo) in oligos {
        let oligo = normalize_dna(raw_oligo)?;
        let olen = oligo.len();
        if olen < min_core {
            continue;
        }
        let wrap_pad = if circular { olen.min(n) } else { 0 };
        let fwd_track = format!("{}{}", plasmid, &plasmid[..wrap_pad]);
        let rev_track = format!("{}{}", rc_plasmid, &rc_plasmid[..wrap_pad]);

        for (strand, track) in [
            (Strand::Forward, fwd_track.as_bytes()),
            (Strand::Reverse, rev_track.as_bytes()),
        ] {
            let ob = oligo.as_bytes();
            let max_start = if circular {
                n
            } else {
                n.saturating_sub(min_core) + 1
            };

            for track_3p_pos in (min_core - 1)..(max_start + olen).min(track.len()) {
                let start_limit = track_3p_pos + 1;
                let max_bound = olen.min(start_limit);
                if !circular && track_3p_pos >= n {
                    break;
                }
                if circular && (track_3p_pos - (min_core - 1)) >= n {
                    continue;
                }

                if ob[olen - 1] != track[track_3p_pos] {
                    continue;
                }

                let mut bound = 1usize;
                let mut mm = 0usize;
                while bound < max_bound {
                    let ob_base = ob[olen - 1 - bound];
                    let tr_base = track[track_3p_pos - bound];
                    if ob_base != tr_base {
                        if mm + 1 > max_internal_mismatches {
                            break;
                        }
                        mm += 1;
                    }
                    bound += 1;
                }

                while bound > 1 && ob[olen - bound] != track[track_3p_pos + 1 - bound] {
                    bound -= 1;
                    mm = mm.saturating_sub(1);
                }
                if bound < min_core {
                    continue;
                }
                let anchor_0 = track_3p_pos + 1 - bound;
                if circular && anchor_0 >= n {
                    continue;
                }

                let overhang = oligo[..(olen - bound)].to_string();
                let core = oligo[(olen - bound)..].to_string();
                let template_slice =
                    std::str::from_utf8(&track[anchor_0..=track_3p_pos]).unwrap_or("");
                let target_3to5: String = template_slice
                    .bytes()
                    .map(|b| complement_base(b) as char)
                    .collect();

                let core_tm = calc_aligned_duplex(&core, &target_3to5, cond, 60.0)?
                    .tm_salt_corrected_c;
                let full_tm = calc_perfect_duplex(&oligo, cond, 60.0)?.tm_salt_corrected_c;

                let (bind_start_0, bind_end_0) = match strand {
                    Strand::Forward => (anchor_0 % n, track_3p_pos % n),
                    Strand::Reverse => {
                        let top_5p = (n - 1) - (anchor_0 % n);
                        let top_3p = (n - 1) - (track_3p_pos % n);
                        (top_5p, top_3p)
                    }
                };

                hits.push(CircularBindingSite {
                    oligo_name: name.clone(),
                    strand,
                    bind_start_0,
                    bind_end_0,
                    bound_len: bound,
                    overhang_5p_seq: overhang,
                    core_3p_seq: core,
                    mismatches: mm,
                    core_tm_c: (core_tm * 100.0).round() / 100.0,
                    full_tm_c: (full_tm * 100.0).round() / 100.0,
                });
            }
        }
    }

    Ok(hits)
}

/// Picks the lowest-penalty 3' template-binding core (`5' -> 3'`) starting at the 5' end of `seq`.
fn select_optimal_5p_core(
    seq: &str,
    cfg: &GibsonConfig,
    cond: &ReactionConditions,
) -> Result<(String, f64), OligoError> {
    let (min_l, max_l) = cfg.core_len_range;
    if seq.len() < min_l {
        return Err(OligoError::SequenceTooShort {
            len: seq.len(),
            min: min_l,
        });
    }
    let upper_l = max_l.min(seq.len());
    let mut best: Option<(f64, String, f64)> = None;

    for len in min_l..=upper_l {
        let cand = &seq[..len];
        let qc = evaluate_oligo_qc(cand, cond, cfg.target_core_tm_c)?;
        let tm_err = (qc.tm_c - cfg.target_core_tm_c).abs();
        let gc_clamp_bonus = if matches!(cand.as_bytes()[len - 1], b'G' | b'C') {
            -0.6
        } else {
            0.5
        };
        let score = tm_err * 2.5 + qc.penalty_score * 0.4 + gc_clamp_bonus;
        if best
            .as_ref()
            .map_or(true, |(b_score, _, _)| score < *b_score)
        {
            best = Some((score, cand.to_string(), qc.tm_c));
        }
    }

    let (_, chosen, tm) = best.ok_or(OligoError::AssemblyDesignFailed {
        junction_index: 0,
        reason: String::from("No valid primer candidate within configured length bounds."),
    })?;
    Ok((chosen, tm))
}

/// Designs a complete two-stage Gibson / NEBuilder HiFi assembly across `fragments`
/// (`>= 1` fragment for circular intramolecular closure, `>= 2` fragments for linear assembly).
pub fn design_gibson_assembly(
    fragments: &[GibsonFragmentInput],
    circular: bool,
    cfg: &GibsonConfig,
) -> Result<GibsonAssemblyPlan, OligoError> {
    let min_frags = if circular { 1 } else { 2 };
    if fragments.len() < min_frags {
        return Err(OligoError::AssemblyDesignFailed {
            junction_index: 0,
            reason: format!(
                "At least {min_frags} fragment(s) required when circular = {circular}."
            ),
        });
    }

    let mut clean_frags = Vec::with_capacity(fragments.len());
    for (idx, f) in fragments.iter().enumerate() {
        let norm = normalize_dna(&f.sequence)?;
        if norm.len() < cfg.core_len_range.1 * 2 {
            return Err(OligoError::AssemblyDesignFailed {
                junction_index: idx,
                reason: format!(
                    "Fragment '{}' ({} bp) is shorter than twice the max core primer length ({} bp).",
                    f.name,
                    norm.len(),
                    cfg.core_len_range.1 * 2
                ),
            });
        }
        clean_frags.push((f.name.clone(), norm));
    }

    let pcr_cond = ReactionConditions::for_preset(cfg.polymerase);
    let gibson_50c_cond = ReactionConditions {
        na_mm: 0.0,
        k_mm: 25.0,
        tris_mm: 50.0,
        mg_mm: 10.0,
        dntp_mm: 0.8,
        dmso_pct: 0.0,
        formamide_pct: 0.0,
        strand_a_nm: 50.0,
        strand_b_nm: 50.0,
        salt_method: SaltCorrectionMethod::Owczarzy2008,
    };

    let n_frags = clean_frags.len();
    let n_junctions = if circular { n_frags } else { n_frags - 1 };

    let mut fwd_cores = Vec::with_capacity(n_frags);
    let mut rev_cores = Vec::with_capacity(n_frags);
    for (_, seq) in &clean_frags {
        let (f_core, f_tm) = select_optimal_5p_core(seq, cfg, &pcr_cond)?;
        let rc_seq = reverse_complement(seq);
        let (r_core, r_tm) = select_optimal_5p_core(&rc_seq, cfg, &pcr_cond)?;
        fwd_cores.push((f_core, f_tm));
        rev_cores.push((r_core, r_tm));
    }

    let mut junctions = Vec::with_capacity(n_junctions);
    for j_idx in 0..n_junctions {
        let up_idx = j_idx;
        let dn_idx = (j_idx + 1) % n_frags;
        let (up_name, up_seq) = &clean_frags[up_idx];
        let (dn_name, dn_seq) = &clean_frags[dn_idx];

        let (min_ov, max_ov) = cfg.overlap_len_range;
        let mut chosen_plan: Option<(f64, GibsonJunctionPlan)> = None;

        for ov_len in min_ov..=max_ov {
            let left_half = ov_len / 2;
            let right_half = ov_len - left_half;
            if left_half > up_seq.len() || right_half > dn_seq.len() {
                continue;
            }

            let up_tail_top = &up_seq[(up_seq.len() - left_half)..];
            let dn_head_top = &dn_seq[..right_half];
            let overlap_seq = format!("{up_tail_top}{dn_head_top}");

            let ov_thermo = calc_perfect_duplex(&overlap_seq, &gibson_50c_cond, 50.0)?;
            let ov_rc = reverse_complement(&overlap_seq);
            let hp_top = fold_hairpin(&overlap_seq, &gibson_50c_cond, 50.0)?
                .map(|h| h.delta_g_kcal)
                .unwrap_or(0.0);
            let hp_bot = fold_hairpin(&ov_rc, &gibson_50c_cond, 50.0)?
                .map(|h| h.delta_g_kcal)
                .unwrap_or(0.0);
            let worst_ss_hp = hp_top.min(hp_bot);

            let dn_fwd_tail = up_tail_top.to_ascii_lowercase();
            let (ref dn_core, dn_core_tm) = fwd_cores[dn_idx];
            let dn_full_upper = format!("{up_tail_top}{dn_core}");

            let up_rev_tail_upper = reverse_complement(dn_head_top);
            let up_rev_tail = up_rev_tail_upper.to_ascii_lowercase();
            let (ref up_core, up_core_tm) = rev_cores[up_idx];
            let up_full_upper = format!("{up_rev_tail_upper}{up_core}");

            let dn_qc = evaluate_oligo_qc(&dn_full_upper, &pcr_cond, dn_core_tm)?;
            let up_qc = evaluate_oligo_qc(&up_full_upper, &pcr_cond, up_core_tm)?;
            let (_, cross_3p) = score_dimer(&dn_full_upper, &up_full_upper, &pcr_cond, 37.0)?;
            let cross_dg = cross_3p.map(|d| d.delta_g_kcal).unwrap_or(0.0);

            let mut j_penalty = 0.0_f64;
            let mut warnings = Vec::new();

            if ov_thermo.tm_salt_corrected_c < cfg.min_overlap_tm_c {
                j_penalty += (cfg.min_overlap_tm_c - ov_thermo.tm_salt_corrected_c) * 8.0;
                warnings.push(format!(
                    "Overlap Tm ({:.1} C) below {:.1} C isothermal target.",
                    ov_thermo.tm_salt_corrected_c, cfg.min_overlap_tm_c
                ));
            }
            if ov_thermo.delta_g_at_temp_kcal > cfg.max_overlap_dg_50c_kcal {
                j_penalty += (ov_thermo.delta_g_at_temp_kcal - cfg.max_overlap_dg_50c_kcal) * 5.0;
                warnings.push(format!(
                    "Overlap Delta G_50C ({:.1} kcal/mol) weaker than {:.1} kcal/mol target.",
                    ov_thermo.delta_g_at_temp_kcal, cfg.max_overlap_dg_50c_kcal
                ));
            }
            if worst_ss_hp < -3.0 {
                j_penalty += (-worst_ss_hp) * 4.0;
                warnings.push(format!(
                    "Exonuclease-exposed 3' overhang forms hairpin at 50 C (Delta G = {:.2} kcal/mol).",
                    worst_ss_hp
                ));
            }
            if cross_dg < -5.0 {
                j_penalty += (-cross_dg) * 3.0;
                warnings.push(format!(
                    "Junction primers share 3'-extensible cross-dimer (Delta G = {:.2} kcal/mol).",
                    cross_dg
                ));
            }

            j_penalty +=
                (ov_len as f64) * 0.15 + (dn_qc.penalty_score + up_qc.penalty_score) * 0.25;

            let dn_primer = ChimericPrimer {
                name: format!("J{}_{}_F", j_idx + 1, dn_name),
                annotated_sequence: format!("{dn_fwd_tail}{dn_core}"),
                tail_5p_lower: dn_fwd_tail,
                core_3p_upper: dn_core.clone(),
                full_sequence_upper: dn_full_upper,
                length: dn_qc.length,
                core_tm_c: (dn_core_tm * 100.0).round() / 100.0,
                full_tm_c: dn_qc.tm_c,
                gc_pct: dn_qc.gc_pct,
                hairpin_dg_kcal: dn_qc.hairpin.map(|h| h.delta_g_kcal).unwrap_or(0.0),
                self_dimer_3p_dg_kcal: dn_qc.self_dimer_3p.map(|d| d.delta_g_kcal).unwrap_or(0.0),
                g4_score: dn_qc.g4_hunter_max,
            };

            let up_primer = ChimericPrimer {
                name: format!("J{}_{}_R", j_idx + 1, up_name),
                annotated_sequence: format!("{up_rev_tail}{up_core}"),
                tail_5p_lower: up_rev_tail,
                core_3p_upper: up_core.clone(),
                full_sequence_upper: up_full_upper,
                length: up_qc.length,
                core_tm_c: (up_core_tm * 100.0).round() / 100.0,
                full_tm_c: up_qc.tm_c,
                gc_pct: up_qc.gc_pct,
                hairpin_dg_kcal: up_qc.hairpin.map(|h| h.delta_g_kcal).unwrap_or(0.0),
                self_dimer_3p_dg_kcal: up_qc.self_dimer_3p.map(|d| d.delta_g_kcal).unwrap_or(0.0),
                g4_score: up_qc.g4_hunter_max,
            };

            let plan = GibsonJunctionPlan {
                junction_index: j_idx,
                upstream_fragment: up_name.clone(),
                downstream_fragment: dn_name.clone(),
                overlap_sequence: overlap_seq,
                overlap_len: ov_len,
                overlap_tm_c: (ov_thermo.tm_salt_corrected_c * 100.0).round() / 100.0,
                overlap_dg_50c_kcal: (ov_thermo.delta_g_at_temp_kcal * 100.0).round() / 100.0,
                ssdna_hairpin_dg_50c_kcal: worst_ss_hp,
                upstream_rev_primer: up_primer,
                downstream_fwd_primer: dn_primer,
                cross_dimer_3p_dg_kcal: cross_dg,
                warnings,
            };

            if chosen_plan
                .as_ref()
                .map_or(true, |(best_pen, _)| j_penalty < *best_pen)
            {
                chosen_plan = Some((j_penalty, plan));
            }
        }

        let (_, best_j) = chosen_plan.ok_or_else(|| OligoError::AssemblyDesignFailed {
            junction_index: j_idx,
            reason: format!("No valid overlap window found for junction {up_name} -> {dn_name}."),
        })?;
        junctions.push(best_j);
    }

    let mut fragment_pcrs = Vec::with_capacity(n_frags);
    for f_idx in 0..n_frags {
        let (ref f_name, ref f_seq) = clean_frags[f_idx];
        let (fwd_name, fwd_tail_bp, fwd_core) = if circular {
            let prev_j = &junctions[(f_idx + n_frags - 1) % n_frags];
            (
                prev_j.downstream_fwd_primer.name.clone(),
                prev_j.downstream_fwd_primer.tail_5p_lower.len(),
                prev_j.downstream_fwd_primer.core_3p_upper.clone(),
            )
        } else if f_idx == 0 {
            (
                format!("Linear_{f_name}_5p_F"),
                0usize,
                fwd_cores[0].0.clone(),
            )
        } else {
            let prev_j = &junctions[f_idx - 1];
            (
                prev_j.downstream_fwd_primer.name.clone(),
                prev_j.downstream_fwd_primer.tail_5p_lower.len(),
                prev_j.downstream_fwd_primer.core_3p_upper.clone(),
            )
        };

        let (rev_name, rev_tail_bp, rev_core) = if circular || f_idx + 1 < n_frags {
            let next_j = &junctions[f_idx];
            (
                next_j.upstream_rev_primer.name.clone(),
                next_j.upstream_rev_primer.tail_5p_lower.len(),
                next_j.upstream_rev_primer.core_3p_upper.clone(),
            )
        } else {
            (
                format!("Linear_{f_name}_3p_R"),
                0usize,
                rev_cores[f_idx].0.clone(),
            )
        };

        let ta = recommend_annealing_temp(cfg.polymerase, &fwd_core, &rev_core)?;
        let amp_len = f_seq.len() + fwd_tail_bp + rev_tail_bp;
        let sec_per_kb = match cfg.polymerase {
            PolymerasePreset::Q5 | PolymerasePreset::Phusion => 30.0,
            _ => 60.0,
        };
        let ext_sec = (((amp_len as f64 / 1000.0) * sec_per_kb).ceil() as u32).max(15);

        fragment_pcrs.push(FragmentPcrProtocol {
            fragment_name: f_name.clone(),
            fwd_primer_name: fwd_name,
            rev_primer_name: rev_name,
            amplicon_len_bp: amp_len,
            recommended_ta_c: ta.recommended_ta_c,
            extension_sec: ext_sec,
        });
    }

    let assembled_sequence: String = clean_frags.iter().map(|(_, s)| s.as_str()).collect();
    let assembled_length_bp = assembled_sequence.len();

    Ok(GibsonAssemblyPlan {
        circular,
        assembled_sequence,
        assembled_length_bp,
        junctions,
        fragment_pcrs,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scans_origin_wrapping_circular_primer_with_5p_tail() {
        let plasmid = "GCATGCATCGATCGATCGTAGCTAGCTAGCTGACTGACTGACGTACGTAGCTAGCTAGCGATCG";
        let origin_span = format!("{}{}", &plasmid[54..64], &plasmid[0..10]);
        let tailed_primer = format!("GGATCC{origin_span}");
        let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);

        let hits = scan_circular_bindings(
            plasmid,
            true,
            &[("OriginWrap_F".to_string(), tailed_primer)],
            15,
            1,
            &cond,
        )
        .unwrap();

        let fwd_hit = hits
            .iter()
            .find(|h| h.strand == Strand::Forward && h.bind_start_0 == 54)
            .expect("Expected origin-wrapping forward binding site");
        assert_eq!(fwd_hit.bind_end_0, 9);
        assert_eq!(fwd_hit.bound_len, 20);
        assert_eq!(fwd_hit.overhang_5p_seq, "GGATCC");
        assert_eq!(fwd_hit.mismatches, 0);
    }

    #[test]
    fn designs_seamless_two_part_circular_gibson_plasmid() {
        let backbone = "GCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGC\
                        GATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGTAGCTAGCTAGCGATCGATCGT";
        let gfp_insert = "ATGGTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAAC\
                          GGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATC";

        let frags = vec![
            GibsonFragmentInput {
                name: "pCRV_Backbone".to_string(),
                sequence: backbone.to_string(),
            },
            GibsonFragmentInput {
                name: "eGFP_CDS".to_string(),
                sequence: gfp_insert.to_string(),
            },
        ];

        let plan = design_gibson_assembly(&frags, true, &GibsonConfig::default()).unwrap();
        assert_eq!(plan.junctions.len(), 2);
        assert_eq!(plan.assembled_length_bp, backbone.len() + gfp_insert.len());

        for j in &plan.junctions {
            assert!(j.overlap_len >= 20);
            assert!(j.overlap_tm_c >= 50.0);
            assert!(j.overlap_dg_50c_kcal <= -18.0);
            assert!(plan.to_order_tsv().contains(&j.downstream_fwd_primer.name));
        }
    }
}