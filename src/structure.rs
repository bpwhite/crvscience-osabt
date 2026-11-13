// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)

//! Oligonucleotide secondary-structure and failure-mode screening:
//! - Intramolecular **hairpin** stem-loop thermodynamics (SantaLucia & Hicks, 2004)
//! - Intermolecular **self-dimer** and **cross-dimer** screening with explicit 3'-OH extensibility
//! - **G-quadruplex (`G4Hunter`)** polymerase-stall motif detection (Bedrat et al., 2016)
//! - 3'-end pentamer free-energy (`Delta G^0_37`) GC-clamp auditing

use crate::thermo::{
    calc_perfect_duplex, internal_mismatch_nn, wc_nn_stack, wc_terminal_initiation,
    ReactionConditions, KELVIN_OFFSET,
};
use crate::{complement_base, normalize_dna, OligoError};

/// Most stable intramolecular hairpin stem-loop conformation identified in an oligonucleotide.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct HairpinHit {
    /// Net folding free energy `Delta G` at `eval_temp_c` (`kcal/mol`; negative = spontaneous).
    pub delta_g_kcal: f64,
    /// Unimolecular hairpin melting temperature `Tm` (`C`) where `Delta G = 0`.
    pub tm_c: f64,
    /// Contiguous base-paired stem length (`bp`).
    pub stem_len: usize,
    /// Unpaired apical hairpin loop length (`nt`).
    pub loop_len: usize,
    /// 0-based start coordinate of the 5' stem arm.
    pub stem_5p_start: usize,
    /// Human-readable ASCII stem-loop diagram for terminal/GUI inspection.
    pub ascii_diagram: String,
}

/// Most stable intermolecular antiparallel duplex alignment between two strands.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DimerHit {
    /// Estimated dimerization free energy `Delta G` at `eval_temp_c` (`kcal/mol`).
    pub delta_g_kcal: f64,
    /// Longest contiguous Watson–Crick stem inside the alignment (`bp`).
    pub contiguous_wc_max: usize,
    /// Total paired bases (Watson–Crick + G·T wobble) across the alignment.
    pub paired_bases: usize,
    /// True if at least the last 2 bases of either strand's 3'-OH terminus are paired,
    /// allowing DNA polymerase extension into primer-dimer artifacts.
    pub extensible_3p: bool,
    /// Relative shift of bottom strand (`3' -> 5'`) against top strand (`5' -> 3'`).
    pub offset: isize,
    /// Three-line ASCII alignment (`5'->3'`, `|` WC / `:` G·T bonds, `3'<-5'`).
    pub ascii_alignment: String,
}

/// Comprehensive oligonucleotide quality-control report used by `bench-tui` and `open-plasmid`.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OligoQcReport {
    pub sequence: String,
    pub length: usize,
    pub gc_pct: f64,
    pub tm_c: f64,
    pub three_prime_dg_37_kcal: f64,
    pub hairpin: Option<HairpinHit>,
    pub self_dimer_any: Option<DimerHit>,
    pub self_dimer_3p: Option<DimerHit>,
    pub g4_hunter_max: f64,
    pub has_g4_risk: bool,
    pub max_homopolymer: usize,
    pub penalty_score: f64,
    pub warnings: Vec<String>,
}

/// SantaLucia & Hicks (2004, *Annu. Rev. Biophys.* 33:415, Table 4) DNA hairpin loop
/// initiation free energy `Delta G^0_37` (`kcal/mol`) by unpaired loop length (`nt`).
fn dna_hairpin_loop_dg37(loop_len: usize) -> f64 {
    match loop_len {
        3 => 3.5,
        4 => 3.5,
        5 => 3.3,
        6 => 4.0,
        7 => 4.2,
        8 => 4.3,
        9 => 4.5,
        10 => 4.6,
        11 => 4.8,
        12 => 4.9,
        n if n > 12 => 4.9 + 1.75 * 1.987e-3 * 310.15 * ((n as f64) / 12.0).ln(),
        _ => 6.0,
    }
}

/// Stabilizing DNA triloop and tetraloop sequence bonuses (`kcal/mol`, SantaLucia & Hicks 2004).
fn special_loop_bonus_dg37(closing_5p: u8, loop_seq: &[u8], closing_3p: u8) -> f64 {
    if loop_seq.len() == 4 {
        match (closing_5p, loop_seq, closing_3p) {
            (b'C', b"GAAA", b'G')
            | (b'C', b"GCAA", b'G')
            | (b'C', b"GAGA", b'G')
            | (b'C', b"GTAA", b'G')
            | (b'G', b"GAAA", b'C')
            | (b'G', b"TGCG", b'C') => -1.5,
            (_, b"GAAA", _) | (_, b"GCAA", _) | (_, b"GAGA", _) | (_, b"GTAA", _) => -1.0,
            _ => 0.0,
        }
    } else if loop_seq.len() == 3 {
        match (closing_5p, loop_seq, closing_3p) {
            (b'C', b"GAA", b'G') | (b'G', b"GAA", b'C') => -1.0,
            _ => 0.0,
        }
    } else {
        0.0
    }
}

/// Scans all candidate unimolecular stem-loop conformations (`stem >= 3 bp`, `loop 3..=14 nt`)
/// and returns the lowest-free-energy hairpin at `eval_temp_c` if `Delta G < 0.5 kcal/mol`.
pub fn fold_hairpin(
    raw_seq: &str,
    cond: &ReactionConditions,
    eval_temp_c: f64,
) -> Result<Option<HairpinHit>, OligoError> {
    let seq = normalize_dna(raw_seq)?;
    let b = seq.as_bytes();
    let n = b.len();
    if n < 9 {
        return Ok(None);
    }

    let na_eq = cond.sodium_equivalent_molar().max(1.0e-3);
    let ln_na = na_eq.ln();
    let eval_k = eval_temp_c + KELVIN_OFFSET;

    let mut best: Option<HairpinHit> = None;

    for loop_start in 3..(n - 5) {
        for loop_len in 3..=14 {
            let loop_end = loop_start + loop_len;
            if loop_end + 3 > n {
                break;
            }

            let max_stem = loop_start.min(n - loop_end);
            let mut stem_len = 0usize;
            while stem_len < max_stem {
                let left = b[loop_start - 1 - stem_len];
                let right = b[loop_end + stem_len];
                if complement_base(left) == right {
                    stem_len += 1;
                } else {
                    break;
                }
            }
            if stem_len < 3 {
                continue;
            }

            let s5_start = loop_start - stem_len;
            let mut dh_stem = 0.0;
            let mut ds_stem = 0.0;
            for k in 0..(stem_len - 1) {
                let p1 = b[s5_start + k];
                let p2 = b[s5_start + k + 1];
                if let Some((dh, ds)) = wc_nn_stack(p1, p2) {
                    dh_stem += dh;
                    ds_stem += ds;
                }
            }

            let close_5p = b[loop_start - 1];
            let close_3p = b[loop_end];
            if matches!(close_5p, b'A' | b'T') {
                dh_stem += 0.5;
            }
            if matches!(b[s5_start], b'A' | b'T') {
                dh_stem += 0.5;
            }

            let loop_dg37 = dna_hairpin_loop_dg37(loop_len)
                + special_loop_bonus_dg37(close_5p, &b[loop_start..loop_end], close_3p);
            let ds_loop = -(loop_dg37 * 1000.0) / (37.0 + KELVIN_OFFSET);

            let ds_salt = 0.368 * (stem_len as f64) * ln_na;
            let total_dh = dh_stem;
            let total_ds = ds_stem + ds_loop + ds_salt;

            let dg_eval = total_dh - eval_k * (total_ds * 1.0e-3);
            if dg_eval >= 0.5 || total_ds >= 0.0 {
                continue;
            }

            let tm_c = (total_dh * 1000.0) / total_ds - KELVIN_OFFSET;

            if best.as_ref().map_or(true, |cur| dg_eval < cur.delta_g_kcal) {
                let arm_5p = &seq[s5_start..loop_start];
                let loop_str = &seq[loop_start..loop_end];
                let arm_3p_rev: String =
                    seq[loop_end..(loop_end + stem_len)].chars().rev().collect();
                let bonds = "|".repeat(stem_len);
                let diagram = format!("5'-{arm_5p} ({loop_str})\n   {bonds}\n3'-{arm_3p_rev}");
                best = Some(HairpinHit {
                    delta_g_kcal: (dg_eval * 100.0).round() / 100.0,
                    tm_c: (tm_c * 10.0).round() / 10.0,
                    stem_len,
                    loop_len,
                    stem_5p_start: s5_start,
                    ascii_diagram: diagram,
                });
            }
        }
    }

    Ok(best)
}

#[inline]
fn is_wc_or_wobble(a: u8, b: u8) -> Option<char> {
    if complement_base(a) == b {
        Some('|')
    } else if matches!((a, b), (b'G', b'T') | (b'T', b'G')) {
        Some(':')
    } else {
        None
    }
}

/// Evaluates all antiparallel alignments of `seq_a` (`5' -> 3'`) against `seq_b` (`3' <- 5'`)
/// and returns `(worst_overall_dimer, worst_3p_extensible_dimer)`.
pub fn score_dimer(
    raw_seq_a: &str,
    raw_seq_b: &str,
    cond: &ReactionConditions,
    eval_temp_c: f64,
) -> Result<(Option<DimerHit>, Option<DimerHit>), OligoError> {
    let seq_a = normalize_dna(raw_seq_a)?;
    let seq_b = normalize_dna(raw_seq_b)?;
    let a = seq_a.as_bytes();
    let b_rev: Vec<u8> = seq_b.bytes().rev().collect();

    let na = a.len() as isize;
    let nb = b_rev.len() as isize;
    if na < 4 || nb < 4 {
        return Ok((None, None));
    }

    let na_eq = cond.sodium_equivalent_molar().max(1.0e-3);
    let ln_na = na_eq.ln();
    let eval_k = eval_temp_c + KELVIN_OFFSET;

    let mut worst_any: Option<DimerHit> = None;
    let mut worst_3p: Option<DimerHit> = None;

    for offset in (-(nb - 2))..=(na - 2) {
        let ov_start_a = offset.max(0) as usize;
        let ov_end_a = (offset + nb).min(na) as usize;
        if ov_end_a.saturating_sub(ov_start_a) < 2 {
            continue;
        }

        let mut paired_bases = 0usize;
        let mut cur_wc = 0usize;
        let mut max_wc = 0usize;

        for ia in ov_start_a..ov_end_a {
            let ib = (ia as isize - offset) as usize;
            match is_wc_or_wobble(a[ia], b_rev[ib]) {
                Some('|') => {
                    paired_bases += 1;
                    cur_wc += 1;
                    max_wc = max_wc.max(cur_wc);
                }
                Some(':') => {
                    paired_bases += 1;
                    cur_wc = 0;
                }
                _ => cur_wc = 0,
            }
        }
        if paired_bases < 3 && max_wc < 2 {
            continue;
        }

        let mut best_seg_dg = 0.0_f64;
        let mut seg_i = ov_start_a;
        while seg_i < ov_end_a {
            let ib0 = (seg_i as isize - offset) as usize;
            if is_wc_or_wobble(a[seg_i], b_rev[ib0]).is_none() {
                seg_i += 1;
                continue;
            }
            let mut seg_end = seg_i + 1;
            while seg_end < ov_end_a {
                let ib = (seg_end as isize - offset) as usize;
                if is_wc_or_wobble(a[seg_end], b_rev[ib]).is_some() {
                    seg_end += 1;
                } else if seg_end + 1 < ov_end_a {
                    let ib_next = ((seg_end + 1) as isize - offset) as usize;
                    if complement_base(a[seg_end + 1]) == b_rev[ib_next] {
                        seg_end += 2;
                    } else {
                        break;
                    }
                } else {
                    break;
                }
            }

            let seg_len = seg_end - seg_i;
            if seg_len >= 2 {
                let (dh_5, ds_5) = wc_terminal_initiation(a[seg_i]);
                let (dh_3, ds_3) = wc_terminal_initiation(a[seg_end - 1]);
                let mut dh = dh_5 + dh_3;
                let mut ds = ds_5 + ds_3;

                for k in seg_i..(seg_end - 1) {
                    let kb = (k as isize - offset) as usize;
                    let wc_l = complement_base(a[k]) == b_rev[kb];
                    let wc_r = complement_base(a[k + 1]) == b_rev[kb + 1];
                    if wc_l && wc_r {
                        if let Some((h, s)) = wc_nn_stack(a[k], a[k + 1]) {
                            dh += h;
                            ds += s;
                        }
                    } else if let Some((h, s)) =
                        internal_mismatch_nn(a[k], a[k + 1], b_rev[kb], b_rev[kb + 1])
                    {
                        dh += h;
                        ds += s;
                    } else {
                        dh += 1.8;
                        ds += 2.0;
                    }
                }
                let ds_adj = ds + 0.368 * ((seg_len - 1) as f64) * ln_na;
                let dg = dh - eval_k * (ds_adj * 1.0e-3);
                if dg < best_seg_dg {
                    best_seg_dg = dg;
                }
            }
            seg_i = seg_end;
        }

        if best_seg_dg >= -0.5 {
            continue;
        }

        let a_3p_paired = {
            let ia1 = (na - 1) as usize;
            let ia2 = (na - 2) as usize;
            let ib1 = (na - 1) - offset;
            let ib2 = (na - 2) - offset;
            ib1 >= 0
                && ib1 < nb
                && ib2 >= 0
                && ib2 < nb
                && complement_base(a[ia1]) == b_rev[ib1 as usize]
                && is_wc_or_wobble(a[ia2], b_rev[ib2 as usize]).is_some()
        };
        let b_3p_paired = {
            let ia0 = offset;
            let ia1 = offset + 1;
            ia0 >= 0
                && ia0 < na
                && ia1 >= 0
                && ia1 < na
                && complement_base(a[ia0 as usize]) == b_rev[0]
                && is_wc_or_wobble(a[ia1 as usize], b_rev[1]).is_some()
        };
        let extensible_3p = a_3p_paired || b_3p_paired;

        let ascii = render_dimer_ascii(a, &b_rev, offset);
        let hit = DimerHit {
            delta_g_kcal: (best_seg_dg * 100.0).round() / 100.0,
            contiguous_wc_max: max_wc,
            paired_bases,
            extensible_3p,
            offset,
            ascii_alignment: ascii,
        };

        if worst_any
            .as_ref()
            .map_or(true, |cur| hit.delta_g_kcal < cur.delta_g_kcal)
        {
            worst_any = Some(hit.clone());
        }
        if extensible_3p
            && worst_3p
                .as_ref()
                .map_or(true, |cur| hit.delta_g_kcal < cur.delta_g_kcal)
        {
            worst_3p = Some(hit);
        }
    }

    Ok((worst_any, worst_3p))
}

fn render_dimer_ascii(a: &[u8], b_rev: &[u8], offset: isize) -> String {
    let min_pos = 0isize.min(offset);
    let max_pos = (a.len() as isize).max(offset + b_rev.len() as isize);

    let mut top = String::from("5'-");
    let mut mid = String::from("   ");
    let mut bot = String::from("3'-");

    for pos in min_pos..max_pos {
        let ca = if pos >= 0 && (pos as usize) < a.len() {
            Some(a[pos as usize])
        } else {
            None
        };
        let cb = if pos >= offset && ((pos - offset) as usize) < b_rev.len() {
            Some(b_rev[(pos - offset) as usize])
        } else {
            None
        };

        top.push(ca.map(|x| x as char).unwrap_or(' '));
        bot.push(cb.map(|x| x as char).unwrap_or(' '));
        let bond = match (ca, cb) {
            (Some(x), Some(y)) => is_wc_or_wobble(x, y).unwrap_or(' '),
            _ => ' ',
        };
        mid.push(bond);
    }
    top.push_str("-3'");
    bot.push_str("-5'");
    format!("{top}\n{mid}\n{bot}")
}

/// Computes the maximum absolute sliding-window `G4Hunter` score (Bedrat et al., 2016, *NAR* 44:1746)
/// to detect G-quadruplex (or C-rich i-motif) tracts that stall PCR polymerases (`> 1.5` = high risk).
pub fn g4_hunter_max_window(raw_seq: &str, window: usize) -> f64 {
    let Ok(seq) = normalize_dna(raw_seq) else {
        return 0.0;
    };
    let b = seq.as_bytes();
    if b.is_empty() {
        return 0.0;
    }

    let mut per_base = vec![0.0_f64; b.len()];
    let mut i = 0;
    while i < b.len() {
        if b[i] == b'G' || b[i] == b'C' {
            let target = b[i];
            let mut run = 1usize;
            while i + run < b.len() && b[i + run] == target {
                run += 1;
            }
            let mag = (run.min(4)) as f64;
            let val = if target == b'G' { mag } else { -mag };
            for slot in &mut per_base[i..(i + run)] {
                *slot = val;
            }
            i += run;
        } else {
            i += 1;
        }
    }

    let w = window.min(b.len()).max(1);
    let mut max_abs = 0.0_f64;
    for win in per_base.windows(w) {
        let mean = win.iter().sum::<f64>() / (w as f64);
        max_abs = max_abs.max(mean.abs());
    }
    (max_abs * 100.0).round() / 100.0
}

/// Computes the `Delta G^0_37` (`kcal/mol`) of the 3'-terminal pentamer (`5 bp`) of a primer.
/// Values more negative than `-11.0 kcal/mol` indicate an overly sticky 3' GC clamp prone to mispriming.
pub fn three_prime_pentamer_dg37(raw_seq: &str) -> Result<f64, OligoError> {
    let seq = normalize_dna(raw_seq)?;
    if seq.len() < 5 {
        return Err(OligoError::SequenceTooShort {
            len: seq.len(),
            min: 5,
        });
    }
    let tail = &seq[(seq.len() - 5)..];
    let bytes = tail.as_bytes();
    let (dh_5, ds_5) = wc_terminal_initiation(bytes[0]);
    let (dh_3, ds_3) = wc_terminal_initiation(bytes[4]);
    let mut dh = dh_5 + dh_3;
    let mut ds = ds_5 + ds_3;
    for win in bytes.windows(2) {
        let (h, s) = wc_nn_stack(win[0], win[1]).ok_or(OligoError::InvalidBase {
            byte: win[0] as char,
            pos: 0,
        })?;
        dh += h;
        ds += s;
    }
    Ok(((dh - 310.15 * ds * 1.0e-3) * 100.0).round() / 100.0)
}

fn longest_homopolymer(seq: &str) -> usize {
    let mut best = 0usize;
    let mut cur = 0usize;
    let mut prev = 0u8;
    for &b in seq.as_bytes() {
        if b == prev {
            cur += 1;
        } else {
            prev = b;
            cur = 1;
        }
        best = best.max(cur);
    }
    best
}

/// Runs a complete primer quality-control audit at `anneal_temp_c`.
pub fn evaluate_oligo_qc(
    raw_seq: &str,
    cond: &ReactionConditions,
    anneal_temp_c: f64,
) -> Result<OligoQcReport, OligoError> {
    let seq = normalize_dna(raw_seq)?;
    let seq_len = seq.len();
    if seq_len < 12 {
        return Err(OligoError::SequenceTooShort {
            len: seq_len,
            min: 12,
        });
    }

    let thermo = calc_perfect_duplex(&seq, cond, anneal_temp_c)?;
    let dg_3p = three_prime_pentamer_dg37(&seq)?;
    let hairpin = fold_hairpin(&seq, cond, anneal_temp_c)?;
    let (self_dimer_any, self_dimer_3p) = score_dimer(&seq, &seq, cond, 37.0)?;
    let g4_score = g4_hunter_max_window(&seq, 20);
    let has_g4 = g4_score >= 1.50;
    let homopolymer = longest_homopolymer(&seq);

    let mut penalty = 0.0_f64;
    let mut warnings = Vec::new();

    if !(35.0..=65.0).contains(&thermo.gc_pct) {
        penalty += (thermo.gc_pct - 50.0).abs() * 0.25;
        warnings.push(format!(
            "GC content ({:.1}%) outside 35-65% window.",
            thermo.gc_pct
        ));
    }
    if dg_3p < -11.0 {
        penalty += (-11.0 - dg_3p) * 3.0;
        warnings.push(format!(
            "Sticky 3' pentamer (Delta G_37 = {:.2} kcal/mol < -11.0) risks non-specific priming.",
            dg_3p
        ));
    }
    if let Some(hp) = &hairpin {
        if hp.delta_g_kcal < -3.0 || hp.tm_c > anneal_temp_c - 5.0 {
            penalty += (-hp.delta_g_kcal) * 2.5;
            warnings.push(format!(
                "Stable hairpin (Delta G = {:.2} kcal/mol, Tm = {:.1} C).",
                hp.delta_g_kcal, hp.tm_c
            ));
        }
    }
    if let Some(d3) = &self_dimer_3p {
        if d3.delta_g_kcal < -5.0 {
            penalty += (-d3.delta_g_kcal) * 4.0;
            warnings.push(format!(
                "3'-extensible self-dimer (Delta G_37 = {:.2} kcal/mol).",
                d3.delta_g_kcal
            ));
        }
    }
    if has_g4 {
        penalty += g4_score * 6.0;
        warnings.push(format!(
            "G-quadruplex stall risk detected (G4Hunter = {:.2} >= 1.50).",
            g4_score
        ));
    }
    if homopolymer >= 5 {
        penalty += ((homopolymer - 4) as f64) * 3.0;
        warnings.push(format!("Homopolymer run of {} identical bases.", homopolymer));
    }

    Ok(OligoQcReport {
        sequence: seq,
        length: seq_len,
        gc_pct: (thermo.gc_pct * 10.0).round() / 10.0,
        tm_c: (thermo.tm_salt_corrected_c * 100.0).round() / 100.0,
        three_prime_dg_37_kcal: dg_3p,
        hairpin,
        self_dimer_any,
        self_dimer_3p,
        g4_hunter_max: g4_score,
        has_g4_risk: has_g4,
        max_homopolymer: homopolymer,
        penalty_score: (penalty * 100.0).round() / 100.0,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::thermo::PolymerasePreset;

    #[test]
    fn detects_stable_gnra_tetraloop_hairpin() {
        let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);
        let hp = fold_hairpin("CGCGCGGAAACGCGCG", &cond, 37.0)
            .unwrap()
            .expect("Expected stable hairpin");
        assert!(hp.delta_g_kcal < -5.0);
        assert_eq!(hp.stem_len, 6);
        assert_eq!(hp.loop_len, 4);
        assert!(hp.tm_c > 65.0);
    }

    #[test]
    fn flags_3p_extensible_primer_dimer_and_g4_trap() {
        let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);
        let (_, dim_3p) = score_dimer(
            "ACGTACGTACGTGCATCGAT",
            "ACGTACGTACGTGCATCGAT",
            &cond,
            37.0,
        )
        .unwrap();
        let d3 = dim_3p.expect("Expected 3'-extensible self-dimer");
        assert!(d3.extensible_3p);
        assert!(d3.delta_g_kcal < -7.0);

        let qc = evaluate_oligo_qc("GGGGTAGGGGTAGGGGTAGGGG", &cond, 60.0).unwrap();
        assert!(qc.has_g4_risk);
        assert!(qc.g4_hunter_max > 2.0);
    }
}
