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

//! Tool #17 (`osabt crispr`) — Multi-nuclease CRISPR guide scanner, Rule-Set-2 efficiency model,
//! `CBE`/`ABE` base-editing bystander auditor, and Golden Gate annealed oligo generator:
//! - **On-Target Sequence-Feature Scoring (`SpCas9`, `SaCas9`, `LbCas12a`)** —
//!   Doench et al. (2016) *Nat. Biotechnol.* 34(2):184–191 (`https://pubmed.ncbi.nlm.nih.gov/26780180/`)
//! - **`LbCas12a` (`Cpf1`) `5'-TTTV-3'` Staggered Endonuclease Scanning** —
//!   Zetsche et al. (2015) *Cell* 163(3):759–771 (`https://pubmed.ncbi.nlm.nih.gov/26422227/`)
//! - **Cytosine (`CBE`, `C->T`) & Adenine (`ABE`, `A->G`) Base-Editing Window Auditing** —
//!   Komor et al. (2016) *Nature* 533:420–424 (`https://pubmed.ncbi.nlm.nih.gov/27096365/`) &
//!   Gaudelli et al. (2017) *Nature* 551:464–471 (`https://pmc.ncbi.nlm.nih.gov/articles/PMC5726555/`).

use crate::{gc_fraction, normalize_dna, reverse_complement};

/// Supported RNA-guided CRISPR endonuclease systems.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum CrisprNuclease {
    /// *Streptococcus pyogenes* Cas9: 20-nt spacer + `3'-NGG` PAM, blunt cut `-3` bp upstream of PAM.
    SpCas9,
    /// *Staphylococcus aureus* Cas9: 21-nt spacer + `3'-NNGRRT` PAM (`R = A/G`), compact AAV size.
    SaCas9,
    /// *Lachnospiraceae bacterium* Cas12a (`Cpf1`): `5'-TTTV` PAM (`V = A/C/G`) + 20-nt spacer, staggered cut `+18/+23`.
    LbCas12a,
}

impl CrisprNuclease {
    pub fn label(self) -> &'static str {
        match self {
            Self::SpCas9 => "SpCas9 (20nt + 3'-NGG)",
            Self::SaCas9 => "SaCas9 (21nt + 3'-NNGRRT)",
            Self::LbCas12a => "LbCas12a / Cpf1 (5'-TTTV + 20nt)",
        }
    }

    pub fn spacer_len(self) -> usize {
        match self {
            Self::SpCas9 | Self::LbCas12a => 20,
            Self::SaCas9 => 21,
        }
    }
}

/// Supported Golden Gate sgRNA / crRNA expression vector overhang conventions.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum CrisprCloningVector {
    /// `pSpCas9(BB)-2A-Puro` (`PX459` / `PX330`, `BbsI`): Top `5'-CACCG...`, Bottom `5'-AAAC...C`.
    Px459BbsI,
    /// `lentiCRISPRv2` (`BsmBI` / `Esp3I`): Top `5'-CACCG...`, Bottom `5'-AAAC...C`.
    LentiCrisprV2BsmBI,
    /// `pLbCas12a-crRNA` (`BsaI` direct-repeat overhang): Top `5'-AGAT...`, Bottom `5'-AAAA...`.
    LbCas12aBsaI,
}

impl CrisprCloningVector {
    pub fn label(self) -> &'static str {
        match self {
            Self::Px459BbsI => "PX459 / PX330 (BbsI, U6)",
            Self::LentiCrisprV2BsmBI => "lentiCRISPRv2 (Esp3I/BsmBI, U6)",
            Self::LbCas12aBsaI => "pLbCas12a-crRNA (BsaI)",
        }
    }
}

/// Base-editing window (`positions 4..=8` of `5'->3'` protospacer) target & bystander summary.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BaseEditWindowAudit {
    /// `5`-nt canonical deamination window (`protospacer[3..8]`, 1-indexed `4..=8`).
    pub window_seq_4_to_8: String,
    /// 1-indexed protospacer positions (`4..=8`) containing `C` (`CBE` `C->T` editable).
    pub cbe_c_positions: Vec<usize>,
    /// True if exactly one `C` lies in `4..=8` (zero bystander `C`s) and no `5'-GC-3'` motif inhibits `APOBEC1`.
    pub cbe_clean_single_edit: bool,
    /// 1-indexed protospacer positions (`4..=8`) containing `A` (`ABE` `A->G` editable).
    pub abe_a_positions: Vec<usize>,
    /// True if exactly one `A` lies in `4..=8` (zero bystander `A`s).
    pub abe_clean_single_edit: bool,
}

/// Ready-to-order `5'->3'` annealed oligonucleotide pair for Golden Gate vector cloning.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CrisprGoldenGateOligos {
    pub vector: CrisprCloningVector,
    pub prepended_u6_g: bool,
    pub top_oligo_5to3: String,
    pub bottom_oligo_5to3: String,
}

/// Candidate CRISPR guide identified on a target sequence prior to or paired with off-target scanning.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CrisprGuideCandidate {
    pub guide_id: String,
    pub nuclease: CrisprNuclease,
    /// `'+'` (sense) or `'-'` (antisense).
    pub strand: char,
    /// 1-indexed start coordinate of the protospacer on the `+` target sequence.
    pub protospacer_start_1idx: usize,
    pub protospacer_end_1idx: usize,
    /// 1-indexed double-strand cleavage site on the `+` strand.
    pub cut_site_1idx: usize,
    /// `5'->3'` DNA protospacer sequence (`20` or `21` nt).
    pub spacer_seq: String,
    /// Matched PAM sequence (`3` nt for `SpCas9`, `6` nt for `SaCas9`, `4` nt for `LbCas12a`).
    pub pam_seq: String,
    pub gc_percent: f64,
    /// Sugimoto et al. (1995) RNA:DNA hybrid melting temperature (`°C`) of the spacer R-loop.
    pub rna_dna_tm_c: f64,
    pub has_pol3_poly_t: bool,
    /// Doench et al. (2016) Rule-Set-2-inspired on-target activity score (`0.0..=100.0`).
    pub on_target_score: f64,
    /// Hsu/Doench aggregate `CFD` specificity score (`0.0..=100.0`; `100.0` = zero off-targets).
    pub cfd_specificity_score: f64,
    pub off_target_count: usize,
    pub max_off_target_cfd: f64,
    /// Balanced composite score `0.55 * on_target + 0.45 * cfd_specificity` (penalized if `TTTT`).
    pub composite_score: f64,
    pub base_edit: BaseEditWindowAudit,
}

/// Scans both strands of `raw_target_dna` (`FASTA` or plain sequence) for all valid PAM + protospacer
/// candidates for `nuclease`, computing on-target efficiency, RNA:DNA `Tm`, and `CBE`/`ABE` windows.
pub fn scan_crispr_candidates(
    raw_target_dna: &str,
    nuclease: CrisprNuclease,
) -> Result<(String, Vec<CrisprGuideCandidate>), String> {
    let stripped: String = raw_target_dna
        .lines()
        .filter(|l| !l.trim_start().starts_with('>'))
        .collect::<Vec<_>>()
        .join("");
    let seq = normalize_dna(&stripped).map_err(|e| e.to_string())?;
    if seq.len() < 30 {
        return Err(String::from(
            "Target DNA sequence must be at least 30 bp long for CRISPR scanning.",
        ));
    }

    let fwd = seq.as_bytes();
    let rc_str = reverse_complement(&seq);
    let rev = rc_str.as_bytes();
    let n = fwd.len();

    let mut out = Vec::new();

    for strand in ['+', '-'] {
        let s_bytes = if strand == '+' { fwd } else { rev };
        match nuclease {
            CrisprNuclease::SpCas9 => {
                // 20-nt spacer + 3-nt NGG PAM => window of 23 nt
                for i in 0..=(n.saturating_sub(23)) {
                    let pam = &s_bytes[(i + 20)..(i + 23)];
                    if pam[1] == b'G' && pam[2] == b'G' && is_acgt(&s_bytes[i..(i + 23)]) {
                        let spacer = std::str::from_utf8(&s_bytes[i..(i + 20)]).unwrap();
                        let pam_s = std::str::from_utf8(pam).unwrap();
                        let ctx30 = extract_30mer_context(s_bytes, i);
                        let (p_start, p_end, cut) = if strand == '+' {
                            (i + 1, i + 20, i + 17)
                        } else {
                            // Map RC index `[i..i+20)` back to 1-indexed `+` strand coordinates
                            (n - (i + 20) + 1, n - i, n - (i + 17) + 1)
                        };
                        out.push(build_candidate(
                            out.len() + 1,
                            nuclease,
                            strand,
                            p_start,
                            p_end,
                            cut,
                            spacer,
                            pam_s,
                            ctx30.as_deref(),
                        ));
                    }
                }
            }
            CrisprNuclease::SaCas9 => {
                // 21-nt spacer + 6-nt NNGRRT PAM => window of 27 nt
                for i in 0..=(n.saturating_sub(27)) {
                    let pam = &s_bytes[(i + 21)..(i + 27)];
                    let is_nngrrt = pam[2] == b'G'
                        && matches!(pam[3], b'A' | b'G')
                        && matches!(pam[4], b'A' | b'G')
                        && pam[5] == b'T';
                    if is_nngrrt && is_acgt(&s_bytes[i..(i + 27)]) {
                        let spacer = std::str::from_utf8(&s_bytes[i..(i + 21)]).unwrap();
                        let pam_s = std::str::from_utf8(pam).unwrap();
                        let (p_start, p_end, cut) = if strand == '+' {
                            (i + 1, i + 21, i + 18)
                        } else {
                            (n - (i + 21) + 1, n - i, n - (i + 18) + 1)
                        };
                        out.push(build_candidate(
                            out.len() + 1,
                            nuclease,
                            strand,
                            p_start,
                            p_end,
                            cut,
                            spacer,
                            pam_s,
                            None,
                        ));
                    }
                }
            }
            CrisprNuclease::LbCas12a => {
                // 4-nt 5'-TTTV PAM + 20-nt downstream spacer => window of 24 nt
                for i in 0..=(n.saturating_sub(24)) {
                    let pam = &s_bytes[i..(i + 4)];
                    let is_tttv = pam[0] == b'T'
                        && pam[1] == b'T'
                        && pam[2] == b'T'
                        && matches!(pam[3], b'A' | b'C' | b'G');
                    if is_tttv && is_acgt(&s_bytes[i..(i + 24)]) {
                        let spacer = std::str::from_utf8(&s_bytes[(i + 4)..(i + 24)]).unwrap();
                        let pam_s = std::str::from_utf8(pam).unwrap();
                        let (p_start, p_end, cut) = if strand == '+' {
                            (i + 5, i + 24, i + 22)
                        } else {
                            (n - (i + 24) + 1, n - (i + 4), n - (i + 22) + 1)
                        };
                        out.push(build_candidate(
                            out.len() + 1,
                            nuclease,
                            strand,
                            p_start,
                            p_end,
                            cut,
                            spacer,
                            pam_s,
                            None,
                        ));
                    }
                }
            }
        }
    }

    out.sort_by(|a, b| {
        b.composite_score
            .partial_cmp(&a.composite_score)
            .unwrap_or(std::cmp::Ordering::Equal)
    });
    for (rank, g) in out.iter_mut().enumerate() {
        g.guide_id = format!("gRNA_{:02}", rank + 1);
    }

    Ok((seq, out))
}

#[allow(clippy::too_many_arguments)]
fn build_candidate(
    idx: usize,
    nuclease: CrisprNuclease,
    strand: char,
    p_start: usize,
    p_end: usize,
    cut: usize,
    spacer: &str,
    pam: &str,
    ctx30: Option<&str>,
) -> CrisprGuideCandidate {
    let gc_pct = gc_fraction(spacer) * 100.0;
    let rna_tm = calc_rna_dna_hybrid_tm(spacer);
    let poly_t = spacer.contains("TTTT");
    let on_target = score_sgrna_on_target(spacer, pam, ctx30, nuclease);
    let base_edit = audit_base_editor_window(spacer);
    let composite = if poly_t {
        0.25 * on_target
    } else {
        on_target
    };

    CrisprGuideCandidate {
        guide_id: format!("gRNA_{idx:02}"),
        nuclease,
        strand,
        protospacer_start_1idx: p_start,
        protospacer_end_1idx: p_end,
        cut_site_1idx: cut,
        spacer_seq: spacer.to_string(),
        pam_seq: pam.to_string(),
        gc_percent: gc_pct,
        rna_dna_tm_c: rna_tm,
        has_pol3_poly_t: poly_t,
        on_target_score: on_target,
        cfd_specificity_score: 100.0,
        off_target_count: 0,
        max_off_target_cfd: 0.0,
        composite_score: composite,
        base_edit,
    }
}

/// Computes a sequence-feature on-target activity score (`0.0..=100.0`) grounded in
/// Doench et al. (2016) Rule Set 2 position-specific nucleotide preferences, `GC%` optimality,
/// homopolymer/Pol-III termination filters, and Sugimoto RNA:DNA R-loop `Tm`.
pub fn score_sgrna_on_target(
    spacer: &str,
    pam: &str,
    ctx30: Option<&str>,
    nuclease: CrisprNuclease,
) -> f64 {
    let b = spacer.as_bytes();
    if b.len() < 20 {
        return 0.0;
    }

    let mut linear = 0.25_f64;

    // 1. GC count optimality (Doench 2016: 40%–65% GC peak; steep penalty <30% or >75%)
    let gc = gc_fraction(spacer);
    if (0.40..=0.65).contains(&gc) {
        linear += 0.55;
    } else if (0.35..=0.75).contains(&gc) {
        linear += 0.15;
    } else {
        linear -= 0.75;
    }

    // 2. Position-specific nucleotide weights (Doench 2016 Fig. 2: SpCas9 strongly favors G at pos 20,
    //    C at DSB pos 16–18, A in central positions 9–15, and disfavors T in PAM-proximal 17–20 & C at pos 20)
    match nuclease {
        CrisprNuclease::SpCas9 | CrisprNuclease::SaCas9 => {
            let p20 = b[b.len() - 1];
            linear += match p20 {
                b'G' => 0.68,
                b'A' => 0.12,
                b'C' => -0.32,
                _ => -0.58, // T at pos 20 strongly depleted
            };
            let p17 = b[b.len() - 4]; // Cleavage site (-3/-4 upstream of PAM)
            linear += match p17 {
                b'C' => 0.36,
                b'A' | b'G' => 0.10,
                _ => -0.22,
            };
            // Seed T depletion (positions 17..20)
            let seed_t = b[(b.len() - 4)..].iter().filter(|&&x| x == b'T').count();
            linear -= 0.24 * (seed_t as f64);

            // Mid-spacer purine preference (positions 9..15)
            let mid_ag = b[8..15]
                .iter()
                .filter(|&&x| x == b'A' || x == b'G')
                .count();
            linear += 0.08 * ((mid_ag as f64) - 3.5);

            // Variable PAM nucleotide (`N` of `NGG`): Doench 2016 shows `CGG` / `GGG` vs `TGG`: `TGG` slightly disfavored when preceded by T, `CGG` favored
            if let Some(&pam0) = pam.as_bytes().first() {
                if pam0 == b'C' {
                    linear += 0.18;
                } else if pam0 == b'T' {
                    linear -= 0.14;
                }
            }
            // Downstream +1 base in 30-mer (`pos 27`): G immediately 3' of NGG (`NGGG` extended PAM) inhibits SpCas9!
            if let Some(c30) = ctx30 {
                if c30.len() >= 28 && c30.as_bytes()[27] == b'G' {
                    linear -= 0.38;
                }
            }
        }
        CrisprNuclease::LbCas12a => {
            // Cas12a favors C/G at PAM-proximal seed pos 1..5 immediately downstream of 5'-TTTV
            let seed_cg = b[..5]
                .iter()
                .filter(|&&x| x == b'C' || x == b'G')
                .count();
            linear += 0.16 * ((seed_cg as f64) - 2.0);
            if pam == "TTTC" || pam == "TTTG" {
                linear += 0.18;
            }
        }
    }

    // 3. Pol III U6 terminator (`TTTT`) and 4-mer homopolymer penalties
    if spacer.contains("TTTT") {
        linear -= 1.65;
    }
    if spacer.contains("GGGG") || spacer.contains("CCCC") || spacer.contains("AAAA") {
        linear -= 0.45;
    }

    // 4. R-loop thermodynamic window (Sugimoto RNA:DNA Tm ideally ~48–68 °C)
    let tm = calc_rna_dna_hybrid_tm(spacer);
    if !(42.0..=72.0).contains(&tm) {
        linear -= 0.30;
    }

    let prob = 1.0 / (1.0 + (-linear).exp());
    (100.0 * prob).clamp(1.0, 99.0)
}

/// Computes Sugimoto et al. (1995) *Biochemistry* 34:11211–11216 nearest-neighbor
/// RNA:DNA hybrid melting temperature (`°C`) at `1 M Na+` adjusted to physiological `150 mM` ionic strength (`250 nM` RNP).
pub fn calc_rna_dna_hybrid_tm(dna_protospacer_5to3: &str) -> f64 {
    let b = dna_protospacer_5to3.as_bytes();
    if b.len() < 8 {
        return 0.0;
    }
    let mut dh_kcal = 1.9_f64; // Sugimoto helix initiation enthalpy
    let mut ds_eu = -3.9_f64; // Sugimoto helix initiation entropy (cal / mol / K)

    // Dinucleotide step parameters for 5'-rX rY-3' paired with 3'-dX' dY'-5'
    for win in b.windows(2) {
        let (h, s) = match (win[0], win[1]) {
            (b'A', b'A') => (-7.8, -21.9),
            (b'A', b'C') => (-5.9, -12.3),
            (b'A', b'G') => (-9.1, -23.5),
            (b'A', b'T') => (-8.3, -23.9),
            (b'C', b'A') => (-9.0, -26.1),
            (b'C', b'C') => (-9.3, -23.2),
            (b'C', b'G') => (-16.3, -47.1),
            (b'C', b'T') => (-7.0, -19.7),
            (b'G', b'A') => (-5.5, -13.5),
            (b'G', b'C') => (-8.0, -17.1),
            (b'G', b'G') => (-12.8, -31.9),
            (b'G', b'T') => (-7.8, -21.6),
            (b'T', b'A') => (-7.8, -23.2),
            (b'T', b'C') => (-8.6, -22.9),
            (b'T', b'G') => (-10.4, -28.4),
            (b'T', b'T') => (-11.5, -36.4),
            _ => (-8.0, -22.0),
        };
        dh_kcal += h;
        ds_eu += s;
    }

    let r_gas = 1.9872_f64;
    let c_molar = 2.5e-7_f64;
    let tm_1m = (1000.0 * dh_kcal) / (ds_eu + r_gas * (c_molar / 4.0).ln()) - 273.15;
    // Salt correction from 1.0 M to 0.15 M Na+ (`16.6 * log10(0.15) = -13.68 °C`)
    tm_1m + 16.6 * 0.15_f64.log10()
}

/// Audits the canonical `positions 4..=8` (`1`-indexed from `5'` end of protospacer)
/// deamination window for `CBE` (`C->T`, Komor et al. 2016) and `ABE` (`A->G`, Gaudelli et al. 2017).
pub fn audit_base_editor_window(spacer_5to3: &str) -> BaseEditWindowAudit {
    let b = spacer_5to3.as_bytes();
    if b.len() < 8 {
        return BaseEditWindowAudit {
            window_seq_4_to_8: String::new(),
            cbe_c_positions: Vec::new(),
            cbe_clean_single_edit: false,
            abe_a_positions: Vec::new(),
            abe_clean_single_edit: false,
        };
    }

    let win = String::from_utf8_lossy(&b[3..8]).to_string();
    let mut c_pos = Vec::new();
    let mut a_pos = Vec::new();
    for pos1 in 4..=8 {
        match b[pos1 - 1] {
            b'C' => c_pos.push(pos1),
            b'A' => a_pos.push(pos1),
            _ => {}
        }
    }
    // Rat APOBEC1 (`BE3`/`BE4`) has reduced activity on `5'-GC-3'` motifs; flag clean non-GC single C!
    let cbe_clean = c_pos.len() == 1 && {
        let p = c_pos[0] - 1;
        p == 0 || b[p - 1] != b'G'
    };
    let abe_clean = a_pos.len() == 1;

    BaseEditWindowAudit {
        window_seq_4_to_8: win,
        cbe_c_positions: c_pos,
        cbe_clean_single_edit: cbe_clean,
        abe_a_positions: a_pos,
        abe_clean_single_edit: abe_clean,
    }
}

/// Generates ready-to-anneal `5'->3'` top and bottom Golden Gate oligonucleotides for `vector`.
pub fn design_golden_gate_crispr_oligos(
    spacer_5to3: &str,
    vector: CrisprCloningVector,
) -> CrisprGoldenGateOligos {
    let s = spacer_5to3.to_ascii_uppercase();
    match vector {
        CrisprCloningVector::Px459BbsI | CrisprCloningVector::LentiCrisprV2BsmBI => {
            let needs_g = !s.starts_with('G');
            let core = if needs_g { format!("G{s}") } else { s };
            let rc = reverse_complement(&core);
            CrisprGoldenGateOligos {
                vector,
                prepended_u6_g: needs_g,
                top_oligo_5to3: format!("CACC{core}"),
                bottom_oligo_5to3: format!("AAAC{rc}"),
            }
        }
        CrisprCloningVector::LbCas12aBsaI => {
            let rc = reverse_complement(&s);
            CrisprGoldenGateOligos {
                vector,
                prepended_u6_g: false,
                top_oligo_5to3: format!("AGAT{s}"),
                bottom_oligo_5to3: format!("AAAA{rc}"),
            }
        }
    }
}

fn extract_30mer_context(seq: &[u8], spacer_start: usize) -> Option<String> {
    if spacer_start >= 4 && spacer_start + 26 <= seq.len() {
        let sub = &seq[(spacer_start - 4)..(spacer_start + 26)];
        if is_acgt(sub) {
            return Some(String::from_utf8_lossy(sub).to_string());
        }
    }
    None
}

#[inline]
fn is_acgt(slice: &[u8]) -> bool {
    slice
        .iter()
        .all(|&b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// Built-in 180-bp salmonid *mstnb* (`myostatin-1b`) coding exon target + 360-bp genomic background
/// containing realistic 1-mm PAM-distal, 2-mm seed, and `NAG`-PAM off-target decoys.
pub fn demo_crispr_target_and_offtargets() -> (String, String) {
    let target_fasta = String::from(
        "\
>Oncorhynchus_mstnb_exon1_target
ATGGAGCTGAAACTGATCACCGTGGCCATCCTGTTTGCCTGCAGCCCGGACAGCCAGCCA\
GAGATGGTGGACCAAGCCGACATGGGCCTGGACTGCGACGAGCACTCCACCGAGTCCCGC\
TGCTGCCGCTACCCGCTGACGGTGGACTTCGAGGBFGTTCGGCTGGGACTGGATCATCGC
"
        .replace("BF", "CA"),
    );

    // Background genome includes the exact exon PLUS off-target paralog decoys (`mstna` paralog)
    let background_fasta = format!(
        "{target_fasta}\n\
>Oncorhynchus_mstna_paralog_offtargets\n\
AAACAGCCGGACAGCCAGCCAGAGATGGTGGATCAAGCCGACATGGGCCTGGGCTGCGAC\
GAGCACTCCACCGAGTCCCGCTGCTGCCGCTACCCGCTGACGGTAGACTTCAAGGCATTC\
GGCTGGGACTGGATCATCGCTTTTACCGGACAGCCAGCCAGAGATGGAGGACCAAGCCGAC\n"
    );

    (target_fasta, background_fasta)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn scans_spcas9_cas12a_base_edit_windows_and_golden_gate_oligos() {
        let (target, _) = demo_crispr_target_and_offtargets();
        let (_, sp_guides) = scan_crispr_candidates(&target, CrisprNuclease::SpCas9).unwrap();
        assert!(sp_guides.len() >= 10);
        assert!(sp_guides[0].on_target_score > 55.0);

        let (_, cas12_guides) = scan_crispr_candidates(&target, CrisprNuclease::LbCas12a).unwrap();
        assert!(!cas12_guides.is_empty());

        let oligos = design_golden_gate_crispr_oligos(
            "ACAGCCAGCCAGAGATGGTG",
            CrisprCloningVector::Px459BbsI,
        );
        assert!(oligos.prepended_u6_g);
        assert!(oligos.top_oligo_5to3.starts_with("CACCG"));
        assert!(oligos.bottom_oligo_5to3.starts_with("AAAC"));
        assert!(oligos.bottom_oligo_5to3.ends_with('C'));
    }
}
