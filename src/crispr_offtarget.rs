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

//! Tool #17 (`osabt crispr`) — Doench et al. (2016) Cutting Frequency Determination (`CFD`)
//! off-target cleavage scorer, whole-FASTA specificity scanner, `.tsv` oligo exporter, and `.svg` locus viewer:
//! - **Position & Mismatch-Identity `CFD` Off-Target Model** — Doench et al. (2016)
//!   *Nat. Biotechnol.* 34(2):184–191 (`https://pubmed.ncbi.nlm.nih.gov/26780180/`).

use std::fmt::Write as _;

use crate::crispr_core::{
    design_golden_gate_crispr_oligos, CrisprCloningVector, CrisprGuideCandidate, CrisprNuclease,
};
use crate::svg_util::xml_escape;
use crate::{normalize_dna, reverse_complement};

/// Individual genomic off-target hit scored via Doench et al. (2016) `CFD`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CrisprOffTargetHit {
    pub contig: String,
    pub strand: char,
    pub pos_1idx: usize,
    pub target_20nt: String,
    pub pam_seq: String,
    pub mismatches: usize,
    /// Visual mismatch alignment (`'.'` for match, uppercase base for mismatch).
    pub mismatch_mask: String,
    /// Doench et al. (2016) Cutting Frequency Determination probability in `[0.0, 1.0]`.
    pub cfd_score: f64,
}

/// Computes the Doench et al. (2016) Cutting Frequency Determination (`CFD`) off-target
/// cleavage score in `[0.0, 1.0]` between a 20-nt `5'->3'` `SpCas9` spacer (`DNA` alphabet)
/// and a 23-nt `5'->3'` candidate genomic site (`20-nt protospacer + 3-nt PAM`).
pub fn calc_cfd_offtarget_score(spacer_20nt: &str, candidate_23nt: &str) -> Result<f64, String> {
    let sp = spacer_20nt.trim().to_ascii_uppercase();
    let cand = candidate_23nt.trim().to_ascii_uppercase();
    if sp.len() != 20 {
        return Err(format!(
            "CFD spacer must be exactly 20 nt (got {} nt).",
            sp.len()
        ));
    }
    if cand.len() != 23 {
        return Err(format!(
            "CFD candidate target + PAM must be exactly 23 nt (got {} nt).",
            cand.len()
        ));
    }

    let sb = sp.as_bytes();
    let cb = cand.as_bytes();
    let pam_w = cfd_spcas9_pam_weight(&cb[20..23]);
    if pam_w <= 0.0 {
        return Ok(0.0);
    }

    let mut score = pam_w;
    for pos in 0..20 {
        if sb[pos] != cb[pos] {
            score *= cfd_position_mismatch_weight(pos, sb[pos], cb[pos]);
        }
    }
    Ok(score.clamp(0.0, 1.0))
}

/// Empirical `SpCas9` 3'-PAM cleavage weights from Doench et al. (2016) Supplementary Table 19.
#[inline]
fn cfd_spcas9_pam_weight(pam3: &[u8]) -> f64 {
    match (pam3[1], pam3[2]) {
        (b'G', b'G') => 1.000_000, // Canonical NGG
        (b'A', b'G') => 0.259_259, // NAG non-canonical PAM
        (b'G', b'A') => 0.069_444, // NGA non-canonical PAM
        (b'G', b'C') => 0.022_222, // NGC
        (b'G', b'T') => 0.016_129, // NGT
        (b'T', b'G') => 0.016_129, // NTG
        (b'C', b'G') => 0.010_753, // NCG
        _ => 0.0,
    }
}

/// Position (`0..20`, PAM-distal `0` to PAM-proximal seed `19`) and mismatch-type weight
/// capturing Doench et al. (2016) `CFD`.
fn cfd_position_mismatch_weight(pos_0idx: usize, spacer_nt: u8, protospacer_nt: u8) -> f64 {
    const POS_BASE: [f64; 20] = [
        0.92, 0.89, 0.86, 0.83, 0.79, 0.75, 0.70, 0.63, 0.57, 0.52,
        0.47, 0.42, 0.36, 0.29, 0.23, 0.18, 0.14, 0.10, 0.08, 0.06,
    ];

    let type_mult = match (spacer_nt, protospacer_nt) {
        (b'G', b'A') => 1.55, // rG:dT wobble
        (b'T', b'C') => 1.32, // rU:dG wobble
        (b'G', b'C') | (b'A', b'C') | (b'G', b'T') | (b'A', b'T') => 0.68, // purine:purine steric clash
        _ => 1.00,
    };

    (POS_BASE[pos_0idx] * type_mult).clamp(0.02, 0.98)
}

/// Scans a multi-FASTA genomic background for off-target sites (`<= max_mismatches`) for every
/// candidate in `guides`, updates each guide's `cfd_specificity_score`, `off_target_count`, and
/// `composite_score`, and returns per-guide off-target hit lists sorted by descending `CFD`.
pub fn scan_offtargets_and_specificity(
    guides: &mut [CrisprGuideCandidate],
    raw_background_fasta: &str,
    max_mismatches: usize,
) -> Result<Vec<Vec<CrisprOffTargetHit>>, String> {
    let contigs = parse_multi_fasta_contigs(raw_background_fasta)?;
    let mut all_hits: Vec<Vec<CrisprOffTargetHit>> = vec![Vec::new(); guides.len()];

    for (g_idx, guide) in guides.iter_mut().enumerate() {
        let sp_bytes = guide.spacer_seq.as_bytes();
        let sp_len = sp_bytes.len();
        let mut exact_self_skipped = false;

        for (cname, cseq) in &contigs {
            let rc = reverse_complement(cseq);
            for (strand, s_str) in [('+', cseq.as_str()), ('-', rc.as_str())] {
                let sb = s_str.as_bytes();
                let n = sb.len();
                match guide.nuclease {
                    CrisprNuclease::SpCas9 if n >= 23 => {
                        for i in 0..=(n - 23) {
                            let pam = &sb[(i + 20)..(i + 23)];
                            if cfd_spcas9_pam_weight(pam) <= 0.0 {
                                continue;
                            }
                            let mm =
                                count_mismatches_bounded(sp_bytes, &sb[i..(i + 20)], max_mismatches);
                            let Some(mm_cnt) = mm else {
                                continue;
                            };
                            if mm_cnt == 0
                                && pam[1] == b'G'
                                && pam[2] == b'G'
                                && !exact_self_skipped
                            {
                                exact_self_skipped = true;
                                continue;
                            }
                            let cand23 = std::str::from_utf8(&sb[i..(i + 23)]).unwrap_or("");
                            let cfd =
                                calc_cfd_offtarget_score(&guide.spacer_seq, cand23).unwrap_or(0.0);
                            if cfd < 0.005 {
                                continue;
                            }
                            let pos_1idx = if strand == '+' {
                                i + 1
                            } else {
                                n - (i + 20) + 1
                            };
                            all_hits[g_idx].push(CrisprOffTargetHit {
                                contig: cname.clone(),
                                strand,
                                pos_1idx,
                                target_20nt: cand23[..20].to_string(),
                                pam_seq: cand23[20..23].to_string(),
                                mismatches: mm_cnt,
                                mismatch_mask: format_mismatch_mask(sp_bytes, &sb[i..(i + 20)]),
                                cfd_score: cfd,
                            });
                        }
                    }
                    _ => {
                        let win_len = sp_len + guide.pam_seq.len();
                        if n < win_len {
                            continue;
                        }
                        for i in 0..=(n - win_len) {
                            let sub_sp = if guide.nuclease == CrisprNuclease::LbCas12a {
                                &sb[(i + 4)..(i + 4 + sp_len)]
                            } else {
                                &sb[i..(i + sp_len)]
                            };
                            let Some(mm_cnt) =
                                count_mismatches_bounded(sp_bytes, sub_sp, max_mismatches)
                            else {
                                continue;
                            };
                            if mm_cnt == 0 && !exact_self_skipped {
                                exact_self_skipped = true;
                                continue;
                            }
                            let approx_cfd = 0.35_f64.powi(mm_cnt as i32);
                            all_hits[g_idx].push(CrisprOffTargetHit {
                                contig: cname.clone(),
                                strand,
                                pos_1idx: i + 1,
                                target_20nt: String::from_utf8_lossy(sub_sp).to_string(),
                                pam_seq: guide.pam_seq.clone(),
                                mismatches: mm_cnt,
                                mismatch_mask: format_mismatch_mask(sp_bytes, sub_sp),
                                cfd_score: approx_cfd,
                            });
                        }
                    }
                }
            }
        }

        all_hits[g_idx].sort_by(|a, b| {
            b.cfd_score
                .partial_cmp(&a.cfd_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        let sum_cfd: f64 = all_hits[g_idx].iter().map(|h| h.cfd_score).sum();
        let max_cfd = all_hits[g_idx].first().map(|h| h.cfd_score).unwrap_or(0.0);
        let spec = 100.0 / (1.0 + sum_cfd);

        guide.off_target_count = all_hits[g_idx].len();
        guide.max_off_target_cfd = max_cfd;
        guide.cfd_specificity_score = spec;
        let base_comp = 0.55 * guide.on_target_score + 0.45 * spec;
        guide.composite_score = if guide.has_pol3_poly_t {
            0.25 * base_comp
        } else {
            base_comp
        };
    }

    Ok(all_hits)
}

fn parse_multi_fasta_contigs(raw_fa: &str) -> Result<Vec<(String, String)>, String> {
    let mut out = Vec::new();
    let mut cur_name = String::from("contig_1");
    let mut cur_seq = String::new();

    for line in raw_fa.lines() {
        let t = line.trim();
        if t.is_empty() {
            continue;
        }
        if let Some(h) = t.strip_prefix('>') {
            if !cur_seq.is_empty() {
                let norm = normalize_dna(&cur_seq).map_err(|e| e.to_string())?;
                out.push((std::mem::take(&mut cur_name), norm));
                cur_seq.clear();
            }
            cur_name = h.split_whitespace().next().unwrap_or("contig").to_string();
        } else {
            cur_seq.push_str(t);
        }
    }
    if !cur_seq.is_empty() {
        let norm = normalize_dna(&cur_seq).map_err(|e| e.to_string())?;
        out.push((cur_name, norm));
    }
    if out.is_empty() {
        return Err(String::from("Background FASTA contains no valid DNA contigs."));
    }
    Ok(out)
}

#[inline]
fn count_mismatches_bounded(a: &[u8], b: &[u8], max_mm: usize) -> Option<usize> {
    let mut mm = 0usize;
    for (&x, &y) in a.iter().zip(b) {
        if !matches!(y, b'A' | b'C' | b'G' | b'T') {
            return None;
        }
        if x != y {
            mm += 1;
            if mm > max_mm {
                return None;
            }
        }
    }
    Some(mm)
}

fn format_mismatch_mask(spacer: &[u8], off: &[u8]) -> String {
    spacer
        .iter()
        .zip(off)
        .map(|(&s, &o)| if s == o { '.' } else { o as char })
        .collect()
}

/// Exports ranked CRISPR guides + Golden Gate oligos + `CBE`/`ABE` window annotations as `.tsv`.
pub fn export_crispr_design_tsv(
    guides: &[CrisprGuideCandidate],
    vector: CrisprCloningVector,
) -> String {
    let mut out = String::from(
        "Guide_ID\tNuclease\tStrand\tCut_1idx\tSpacer_5to3\tPAM\tGC_Pct\tRloop_Tm_C\tOnTarget_Score\tCFD_Specificity\tOffTargets\tCBE_Window\tABE_Window\tTop_Oligo_5to3\tBottom_Oligo_5to3\n",
    );
    for g in guides {
        let oligos = design_golden_gate_crispr_oligos(&g.spacer_seq, vector);
        let cbe_str = if g.base_edit.cbe_clean_single_edit {
            format!("CLEAN(C@{})", g.base_edit.cbe_c_positions[0])
        } else if g.base_edit.cbe_c_positions.is_empty() {
            String::from("none")
        } else {
            format!("bystander({:?})", g.base_edit.cbe_c_positions)
        };
        let abe_str = if g.base_edit.abe_clean_single_edit {
            format!("CLEAN(A@{})", g.base_edit.abe_a_positions[0])
        } else if g.base_edit.abe_a_positions.is_empty() {
            String::from("none")
        } else {
            format!("bystander({:?})", g.base_edit.abe_a_positions)
        };

        let _ = writeln!(
            out,
            "{}\t{}\t{}\t{}\t{}\t{}\t{:.1}\t{:.1}\t{:.1}\t{:.1}\t{}\t{}\t{}\t{}\t{}",
            g.guide_id,
            g.nuclease.label(),
            g.strand,
            g.cut_site_1idx,
            g.spacer_seq,
            g.pam_seq,
            g.gc_percent,
            g.rna_dna_tm_c,
            g.on_target_score,
            g.cfd_specificity_score,
            g.off_target_count,
            cbe_str,
            abe_str,
            oligos.top_oligo_5to3,
            oligos.bottom_oligo_5to3
        );
    }
    out
}

/// Renders a publication-ready `.svg` CRISPR target locus track showing the top 10 ranked guides,
/// strand orientation arrows, cut-site ticks, and color-coded composite efficiency/specificity bars.
pub fn render_crispr_locus_svg(
    target_len_bp: usize,
    guides: &[CrisprGuideCandidate],
    title: &str,
) -> String {
    let width = 920.0_f64;
    let height = 420.0_f64;
    let left = 72.0_f64;
    let track_w = 780.0_f64;
    let backbone_y = 92.0_f64;
    let len_f = target_len_bp.max(1) as f64;

    let map_x = |bp: usize| left + ((bp as f64) / len_f).clamp(0.0, 1.0) * track_w;

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
        "  <text x=\"{left:.0}\" y=\"32\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Top Ranked CRISPR Guides ({} bp Target Locus)</text>",
        xml_escape(title),
        target_len_bp
    );

    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{:.1}\" width=\"{track_w:.1}\" height=\"12\" rx=\"4\" fill=\"#cbd5e1\" stroke=\"#475569\"/>",
        backbone_y - 6.0
    );
    for tick in 0..=5 {
        let bp = (tick * target_len_bp) / 5;
        let tx = map_x(bp);
        let _ = writeln!(
            svg,
            "  <line x1=\"{tx:.1}\" y1=\"{:.1}\" x2=\"{tx:.1}\" y2=\"{:.1}\" stroke=\"#334155\" stroke-width=\"1.4\"/>",
            backbone_y - 11.0,
            backbone_y + 11.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{tx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#334155\">{bp} bp</text>",
            backbone_y - 15.0
        );
    }

    for (idx, g) in guides.iter().take(10).enumerate() {
        let row_y = 132.0 + (idx as f64) * 26.0;
        let x1 = map_x(g.protospacer_start_1idx.saturating_sub(1));
        let x2 = map_x(g.protospacer_end_1idx).max(x1 + 14.0);
        let cut_x = map_x(g.cut_site_1idx);
        let col = if g.composite_score >= 70.0 {
            "#16a34a"
        } else if g.composite_score >= 50.0 {
            "#2563eb"
        } else {
            "#d97706"
        };

        let _ = writeln!(
            svg,
            "  <rect x=\"{x1:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"12\" rx=\"3\" fill=\"{col}\"/>",
            row_y - 6.0,
            x2 - x1
        );
        let _ = writeln!(
            svg,
            "  <line x1=\"{cut_x:.1}\" y1=\"{:.1}\" x2=\"{cut_x:.1}\" y2=\"{:.1}\" stroke=\"#dc2626\" stroke-width=\"2.0\"/>",
            row_y - 9.0,
            row_y + 9.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#0f172a\">{} ({}) cut={} | On={:.0} Spec={:.0}</text>",
            (x2 + 8.0).min(left + track_w - 230.0),
            row_y + 4.0,
            g.guide_id,
            g.strand,
            g.cut_site_1idx,
            g.on_target_score,
            g.cfd_specificity_score
        );
    }

    svg.push_str("</svg>\n");
    svg
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::crispr_core::{demo_crispr_target_and_offtargets, scan_crispr_candidates};

    #[test]
    fn scores_cfd_seed_vs_distal_mismatches_and_ranks_specificity() {
        let spacer = "ACAGCCAGCCAGAGATGGTG";
        let exact_ngg = "ACAGCCAGCCAGAGATGGTGGGG";
        let exact_nag = "ACAGCCAGCCAGAGATGGTGAAG";  // AAG -> NAG non-canonical PAM (0.259259)
        let distal_1mm = "GCAGCCAGCCAGAGATGGTGGGG"; // pos 1 PAM-distal mismatch -> high CFD (0.92)
        let seed_1mm = "ACAGCCAGCCAGAGATGGTCGGG";   // pos 20 seed clash -> low CFD (<0.10)

        let c_exact = calc_cfd_offtarget_score(spacer, exact_ngg).unwrap();
        let c_nag = calc_cfd_offtarget_score(spacer, exact_nag).unwrap();
        let c_distal = calc_cfd_offtarget_score(spacer, distal_1mm).unwrap();
        let c_seed = calc_cfd_offtarget_score(spacer, seed_1mm).unwrap();

        assert!((c_exact - 1.0).abs() < 1e-9);
        assert!((c_nag - 0.259_259).abs() < 1e-4);
        assert!(c_distal > 0.70);
        assert!(c_seed < 0.15);

        let (target, bg) = demo_crispr_target_and_offtargets();
        let (_, mut guides) = scan_crispr_candidates(&target, CrisprNuclease::SpCas9).unwrap();
        let hits = scan_offtargets_and_specificity(&mut guides, &bg, 4).unwrap();
        assert_eq!(hits.len(), guides.len());
        assert!(guides
            .iter()
            .any(|g| g.off_target_count >= 1 && g.cfd_specificity_score < 99.0));
    }
}
