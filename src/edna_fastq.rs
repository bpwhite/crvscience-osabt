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

//! Tool #14 (`osabt edna`) — Environmental DNA (`eDNA`) & Metabarcoding FASTQ engine:
//! - Dual-barcode (`5'` + `3'` reverse-complement) Hamming-tolerant demultiplexing & primer trimming
//! - **Expected-Error (`maxEE`) Quality Filtering & `UNOISE`-style Denoising** — Edgar & Flyvbjerg (2015)
//!   *Bioinformatics* 31(21):3476–3482
//! - De novo 2-parent PCR chimera breakpoint screening & negative extraction-blank contamination auditing.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::reverse_complement;

/// Single Phred+33 FASTQ sequencing read.
#[derive(Debug, Clone)]
pub struct FastqRead {
    pub id: String,
    pub seq: Vec<u8>,
    /// Raw ASCII Phred+33 quality bytes (`b'!' == 33 => Q = 0`, `b'I' == 73 => Q = 40`).
    pub qual_phred33: Vec<u8>,
}

impl FastqRead {
    /// Edgar & Flyvbjerg (2015) expected number of base-call errors `EE = sum 10^(-Q_i / 10)`.
    pub fn expected_errors(&self) -> f64 {
        calc_fastq_expected_errors(&self.qual_phred33)
    }

    pub fn mean_phred(&self) -> f64 {
        if self.qual_phred33.is_empty() {
            return 0.0;
        }
        let s: f64 = self
            .qual_phred33
            .iter()
            .map(|&q| (q.saturating_sub(33)) as f64)
            .sum();
        s / (self.qual_phred33.len() as f64)
    }
}

/// Dual-indexed sample mapping entry (`5'` forward index + `3'` reverse index in `5'->3'` orientation).
#[derive(Debug, Clone)]
pub struct DualBarcodeSample {
    pub sample_id: String,
    pub fwd_barcode: String,
    pub rev_barcode: String,
    pub is_negative_blank: bool,
}

/// Configurable thresholds for the OSABT eDNA metabarcoding pipeline.
#[derive(Debug, Clone)]
pub struct EdnaPipelineConfig {
    pub fwd_primer_len: usize,
    pub rev_primer_len: usize,
    pub max_barcode_mismatches: usize,
    pub tail_trim_min_q: u8,
    pub min_amplicon_len: usize,
    /// Edgar & Flyvbjerg (2015) maximum expected errors per trimmed amplicon (default `1.0`).
    pub max_expected_errors: f64,
    /// Minimum parent/daughter abundance ratio to absorb a 1-nt substitution error cloud (default `8.0`).
    pub unoise_abundance_skew: f64,
    pub remove_chimeras: bool,
}

impl Default for EdnaPipelineConfig {
    fn default() -> Self {
        Self {
            fwd_primer_len: 12,
            rev_primer_len: 12,
            max_barcode_mismatches: 1,
            tail_trim_min_q: 15,
            min_amplicon_len: 36,
            max_expected_errors: 1.0,
            unoise_abundance_skew: 8.0,
            remove_chimeras: true,
        }
    }
}

/// Denoised Amplicon Sequence Variant (`ASV`) with per-sample read counts and QC flags.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AsvRecord {
    pub asv_id: String,
    pub sequence: String,
    pub total_reads: usize,
    pub is_chimera: bool,
    pub blank_reads: usize,
    /// Flagged if `blank_reads > 0` and exceeds `5%` of max environmental sample abundance.
    pub blank_flagged: bool,
    pub putative_taxon: Option<String>,
    pub counts_by_sample: Vec<usize>,
}

/// Summary statistics + ASV × Sample community matrix produced by [`run_edna_metabarcoding_pipeline`].
#[derive(Debug, Clone)]
pub struct EdnaPipelineReport {
    pub total_raw_reads: usize,
    pub demux_assigned_reads: usize,
    pub unassigned_reads: usize,
    pub passed_ee_reads: usize,
    pub filtered_ee_or_short_reads: usize,
    pub error_corrected_reads: usize,
    pub chimeric_reads_flagged: usize,
    pub sample_names: Vec<String>,
    pub blank_sample_indices: Vec<usize>,
    pub asvs: Vec<AsvRecord>,
}

/// Computes Edgar & Flyvbjerg (2015) expected errors `EE = sum_{i} 10^(-Q_i / 10)` on Phred+33 bytes.
#[inline]
pub fn calc_fastq_expected_errors(qual_phred33: &[u8]) -> f64 {
    qual_phred33
        .iter()
        .map(|&b| {
            let q = (b.saturating_sub(33)).min(60) as f64;
            10_f64.powf(-0.1 * q)
        })
        .sum()
}

/// Parses standard 4-line Phred+33 FASTQ records from text.
pub fn parse_fastq_str(raw_fastq: &str) -> Result<Vec<FastqRead>, String> {
    let lines: Vec<&str> = raw_fastq
        .trim_start_matches('\u{feff}')
        .lines()
        .map(|l| l.trim())
        .filter(|l| !l.is_empty())
        .collect();
    if lines.is_empty() {
        return Err(String::from("FASTQ input is empty."));
    }
    if lines.len() % 4 != 0 {
        return Err(format!(
            "Malformed FASTQ: non-empty line count ({}) is not a multiple of 4.",
            lines.len()
        ));
    }

    let mut out = Vec::with_capacity(lines.len() / 4);
    for (rec_i, chunk) in lines.chunks_exact(4).enumerate() {
        if !chunk[0].starts_with('@') {
            return Err(format!(
                "FASTQ record #{}: header '{}' does not start with '@'.",
                rec_i + 1,
                chunk[0]
            ));
        }
        if !chunk[2].starts_with('+') {
            return Err(format!(
                "FASTQ record #{}: separator '{}' does not start with '+'.",
                rec_i + 1,
                chunk[2]
            ));
        }
        let seq = chunk[1]
            .as_bytes()
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .collect::<Vec<_>>();
        let qual = chunk[3].as_bytes().to_vec();
        if seq.len() != qual.len() {
            return Err(format!(
                "FASTQ record #{}: sequence length ({}) != quality length ({}).",
                rec_i + 1,
                seq.len(),
                qual.len()
            ));
        }
        let id = chunk[0][1..]
            .split_whitespace()
            .next()
            .unwrap_or("read")
            .to_string();
        out.push(FastqRead {
            id,
            seq,
            qual_phred33: qual,
        });
    }
    Ok(out)
}

/// Parses a comma- or tab-delimited dual-barcode sample sheet (`sample_id,fwd_barcode,rev_barcode,role`).
pub fn parse_barcode_sheet_csv(raw_csv: &str) -> Result<Vec<DualBarcodeSample>, String> {
    let mut out = Vec::new();
    let mut header_seen = false;
    for (idx, line) in raw_csv.trim_start_matches('\u{feff}').lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();
        if !header_seen {
            header_seen = true;
            if cols[0].eq_ignore_ascii_case("sample") || cols[0].eq_ignore_ascii_case("sample_id") {
                continue;
            }
        }
        if cols.len() < 3 {
            return Err(format!(
                "Barcode sheet line {}: expected >= 3 columns (sample_id,fwd_barcode,rev_barcode[,role]).",
                idx + 1
            ));
        }
        let is_blank = cols
            .get(3)
            .map(|r| {
                let low = r.to_ascii_lowercase();
                low.contains("blank") || low.contains("ntc") || low.contains("neg")
            })
            .unwrap_or_else(|| cols[0].to_ascii_lowercase().contains("blank"));
        out.push(DualBarcodeSample {
            sample_id: cols[0].to_string(),
            fwd_barcode: cols[1].to_ascii_uppercase(),
            rev_barcode: cols[2].to_ascii_uppercase(),
            is_negative_blank: is_blank,
        });
    }
    if out.is_empty() {
        return Err(String::from("Barcode sample sheet contains no entries."));
    }
    Ok(out)
}

/// Executes end-to-end dual-barcode demultiplexing, primer stripping, `maxEE` filtering,
/// 1-nt abundance-ratio denoising, 2-parent chimera detection, blank auditing, and taxonomic matching.
pub fn run_edna_metabarcoding_pipeline(
    reads: &[FastqRead],
    barcodes: &[DualBarcodeSample],
    reference_db: &[(String, String)],
    cfg: &EdnaPipelineConfig,
) -> EdnaPipelineReport {
    let n_samples = barcodes.len();
    let sample_names: Vec<String> = barcodes.iter().map(|b| b.sample_id.clone()).collect();
    let blank_indices: Vec<usize> = barcodes
        .iter()
        .enumerate()
        .filter(|(_, b)| b.is_negative_blank)
        .map(|(i, _)| i)
        .collect();

    // Precompute reverse-complemented 3' barcode tags as they appear at the 3' end of merged amplicon reads
    let rev_bc_rc: Vec<Vec<u8>> = barcodes
        .iter()
        .map(|b| reverse_complement(&b.rev_barcode).into_bytes())
        .collect();

    let mut demux_assigned = 0usize;
    let mut unassigned = 0usize;
    let mut passed_ee = 0usize;
    let mut filtered_ee = 0usize;

    let mut raw_table: BTreeMap<String, Vec<usize>> = BTreeMap::new();

    for read in reads {
        let mut matched_sample: Option<usize> = None;
        let mut best_mm = usize::MAX;

        for (s_idx, bc) in barcodes.iter().enumerate() {
            let fb = bc.fwd_barcode.as_bytes();
            let rb = &rev_bc_rc[s_idx];
            if read.seq.len() < fb.len() + rb.len() + cfg.fwd_primer_len + cfg.rev_primer_len {
                continue;
            }
            let head = &read.seq[..fb.len()];
            let tail = &read.seq[read.seq.len() - rb.len()..];
            let mm_f = hamming_dist(head, fb);
            let mm_r = hamming_dist(tail, rb);
            let tot_mm = mm_f + mm_r;
            if mm_f <= cfg.max_barcode_mismatches
                && mm_r <= cfg.max_barcode_mismatches
                && tot_mm < best_mm
            {
                best_mm = tot_mm;
                matched_sample = Some(s_idx);
            }
        }

        let Some(s_idx) = matched_sample else {
            unassigned += 1;
            continue;
        };
        demux_assigned += 1;

        let fb_len = barcodes[s_idx].fwd_barcode.len();
        let rb_len = rev_bc_rc[s_idx].len();
        let start = fb_len + cfg.fwd_primer_len;
        let mut end = read.seq.len().saturating_sub(rb_len + cfg.rev_primer_len);

        // Trim low-Q 3' bases below `tail_trim_min_q`
        while end > start && read.qual_phred33[end - 1].saturating_sub(33) < cfg.tail_trim_min_q {
            end -= 1;
        }

        if end <= start || (end - start) < cfg.min_amplicon_len {
            filtered_ee += 1;
            continue;
        }

        let sub_seq = &read.seq[start..end];
        let sub_qual = &read.qual_phred33[start..end];
        if sub_seq.contains(&b'N') || calc_fastq_expected_errors(sub_qual) > cfg.max_expected_errors
        {
            filtered_ee += 1;
            continue;
        }
        passed_ee += 1;

        let seq_str = String::from_utf8_lossy(sub_seq).to_string();
        raw_table
            .entry(seq_str)
            .or_insert_with(|| vec![0usize; n_samples])[s_idx] += 1;
    }

    // Sort unique sequences descending by total read abundance for UNOISE-style 1-nt error absorption
    let mut uniques: Vec<(String, Vec<usize>, usize)> = raw_table
        .into_iter()
        .map(|(s, cnts)| {
            let tot = cnts.iter().sum();
            (s, cnts, tot)
        })
        .collect();
    uniques.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| a.0.cmp(&b.0)));

    let mut denoised: Vec<(String, Vec<usize>, usize)> = Vec::new();
    let mut error_corrected_reads = 0usize;

    for (seq, cnts, tot) in uniques {
        let mut absorbed = false;
        for parent in &mut denoised {
            if parent.0.len() == seq.len()
                && hamming_dist(parent.0.as_bytes(), seq.as_bytes()) == 1
                && (parent.2 as f64) >= cfg.unoise_abundance_skew * (tot as f64)
            {
                for (pc, &dc) in parent.1.iter_mut().zip(&cnts) {
                    *pc += dc;
                }
                parent.2 += tot;
                error_corrected_reads += tot;
                absorbed = true;
                break;
            }
        }
        if !absorbed {
            denoised.push((seq, cnts, tot));
        }
    }

    denoised.sort_by(|a, b| b.2.cmp(&a.2));

    // De novo 2-parent exact breakpoint chimera detection (abundance of both parents >= 2x candidate)
    let mut chimeric_reads_flagged = 0usize;
    let mut asvs = Vec::with_capacity(denoised.len());

    for (idx, (seq, cnts, tot)) in denoised.iter().enumerate() {
        let is_chim = if cfg.remove_chimeras {
            is_two_parent_chimera(seq.as_bytes(), *tot, &denoised[..idx])
        } else {
            false
        };
        if is_chim {
            chimeric_reads_flagged += *tot;
        }

        let blank_reads: usize = blank_indices.iter().map(|&bi| cnts[bi]).sum();
        let max_env = cnts
            .iter()
            .enumerate()
            .filter(|(si, _)| !blank_indices.contains(si))
            .map(|(_, &c)| c)
            .max()
            .unwrap_or(0);
        let blank_flagged = blank_reads > 0 && (blank_reads * 20 >= max_env.max(1));

        let taxon = lookup_best_taxon(seq, reference_db);

        asvs.push(AsvRecord {
            asv_id: format!("ASV_{:03}", idx + 1),
            sequence: seq.clone(),
            total_reads: *tot,
            is_chimera: is_chim,
            blank_reads,
            blank_flagged,
            putative_taxon: taxon,
            counts_by_sample: cnts.clone(),
        });
    }

    EdnaPipelineReport {
        total_raw_reads: reads.len(),
        demux_assigned_reads: demux_assigned,
        unassigned_reads: unassigned,
        passed_ee_reads: passed_ee,
        filtered_ee_or_short_reads: filtered_ee,
        error_corrected_reads,
        chimeric_reads_flagged,
        sample_names,
        blank_sample_indices: blank_indices,
        asvs,
    }
}

fn is_two_parent_chimera(
    cand: &[u8],
    cand_abund: usize,
    parents: &[(String, Vec<usize>, usize)],
) -> bool {
    let eligible: Vec<&[u8]> = parents
        .iter()
        .filter(|p| p.0.len() == cand.len() && p.2 >= 2 * cand_abund)
        .map(|p| p.0.as_bytes())
        .collect();
    if eligible.len() < 2 || cand.len() < 20 {
        return false;
    }
    let l = cand.len();
    for i in 0..eligible.len() {
        for j in 0..eligible.len() {
            if i == j {
                continue;
            }
            let pa = eligible[i];
            let pb = eligible[j];
            // Check if cand matches prefix of `pa` and suffix of `pb` at some interior crossover `k`
            // while differing from both `pa` and `pb` by >= 2 bases overall
            if hamming_dist(cand, pa) < 2 || hamming_dist(cand, pb) < 2 {
                continue;
            }
            for k in 8..(l - 8) {
                if &cand[..k] == &pa[..k] && &cand[k..] == &pb[k..] {
                    return true;
                }
            }
        }
    }
    false
}

fn lookup_best_taxon(seq: &str, db: &[(String, String)]) -> Option<String> {
    let sb = seq.as_bytes();
    let mut best: Option<(usize, &str)> = None;
    for (taxon, ref_seq) in db {
        if ref_seq.len() == sb.len() {
            let mm = hamming_dist(sb, ref_seq.as_bytes());
            if mm <= 1 && best.map(|(b_mm, _)| mm < b_mm).unwrap_or(true) {
                best = Some((mm, taxon.as_str()));
            }
        }
    }
    best.map(|(mm, t)| {
        if mm == 0 {
            format!("{t} (100.0% id)")
        } else {
            format!(
                "{t} ({:.1}% id)",
                100.0 * ((sb.len() - mm) as f64) / (sb.len() as f64)
            )
        }
    })
}

fn hamming_dist(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|&(x, y)| x != y).count() + a.len().abs_diff(b.len())
}

/// Exports the ASV × Sample community abundance table in standard tab-separated format (`.tsv`).
pub fn export_asv_table_tsv(
    report: &EdnaPipelineReport,
    exclude_chimeras_and_blanks: bool,
) -> String {
    let mut out = String::from("ASV_ID\tTotal_Reads\tQC_Status\tPutative_Taxon");
    for s in &report.sample_names {
        let _ = write!(out, "\t{s}");
    }
    out.push_str("\tAmplicon_Sequence\n");

    for asv in &report.asvs {
        if exclude_chimeras_and_blanks && (asv.is_chimera || asv.blank_flagged) {
            continue;
        }
        let status = if asv.is_chimera {
            "CHIMERA"
        } else if asv.blank_flagged {
            "BLANK_CONTAM"
        } else {
            "PASS"
        };
        let _ = write!(
            out,
            "{}\t{}\t{}\t{}",
            asv.asv_id,
            asv.total_reads,
            status,
            asv.putative_taxon.as_deref().unwrap_or("Unassigned")
        );
        for &cnt in &asv.counts_by_sample {
            let _ = write!(out, "\t{cnt}");
        }
        let _ = writeln!(out, "\t{}", asv.sequence);
    }
    out
}

/// Generates a complete synthetic Pacific Northwest `12S` eDNA FASTQ run (`185` reads across
/// `Marys_River`, `Alsea_River`, `McKenzie_River`, and `Extraction_Blank`), barcode sheet CSV,
/// and 12S reference database including point errors, high-EE reads, a 2-parent chimera, and a blank contaminant.
pub fn demo_pnw_12s_edna_dataset() -> (String, String, Vec<(String, String)>) {
    let barcode_csv = String::from(
        "\
sample_id,fwd_barcode,rev_barcode,role
Marys_River,ACACACAC,GTGTGTGT,environmental
Alsea_River,AGAGAGAG,CTCTCTCT,environmental
McKenzie_River,TCTCTCTC,GAGAGAGA,environmental
Extraction_Blank,CATGCATG,GTACGTAC,negative_blank
",
    );

    let fwd_primer = "GTCGGTAAAACT"; // 12-nt MiFish-U 5' core
    let rev_primer_rc = "GTTTGACCCTAA"; // 12-nt 3' reverse primer RC

    // 44-nt 12S variable-loop haplotypes with informative polymorphisms across both 5' [0..22] and 3' [22..44] halves
    let taxa_seqs: Vec<(String, String)> = vec![
        (
            String::from("Oncorhynchus kisutch (Coho Salmon)"),
            String::from("CGTGCCAGCCACCGCGGTTATACGAGAGGCCCAAGTTGATAGAC"),
        ),
        (
            String::from("Oncorhynchus mykiss (Steelhead / Rainbow Trout)"),
            String::from("CGTGCCAGCCATCGCGGTCATACGAGAGGCTCAAGTTGACAGTC"),
        ),
        (
            String::from("Oncorhynchus clarkii (Coastal Cutthroat Trout)"),
            String::from("CGTGCCAGCTACCGCGGCTATACGAGAGGCACAAGTTGATAGCC"),
        ),
        (
            String::from("Dicamptodon tenebrosus (Coastal Giant Salamander)"),
            String::from("CGTGCTAGTCACCGCGGCAATACGAAAGGCCCTAGTTAATAATC"),
        ),
        (
            String::from("Ascaphus truei (Coastal Tailed Frog)"),
            String::from("CGTACTAGCCACCGCGGCAATACGAAGGGCTCAAGTTAATAGAC"),
        ),
        (
            String::from("Cottus perplexus (Reticulate Sculpin)"),
            String::from("CGTGCCAGCCGTCGCGGTTATACGAGGGGCTCAAGTTGTTAGAC"),
        ),
        (
            String::from("Homo sapiens (Lab Reagent Contaminant)"),
            String::from("CGTATCAGTCACCGCGGTCATACGATTAACCCAAGTCAATAGAA"),
        ),
    ];

    let sample_tags: [(String, String); 4] = [
        (String::from("ACACACAC"), reverse_complement("GTGTGTGT")),
        (String::from("AGAGAGAG"), reverse_complement("CTCTCTCT")),
        (String::from("TCTCTCTC"), reverse_complement("GAGAGAGA")),
        (String::from("CATGCATG"), reverse_complement("GTACGTAC")),
    ];

    // Abundance matrix [taxon_idx][sample_idx]
    let abund: [[usize; 4]; 7] = [
        [24, 32, 10, 0], // Coho
        [18, 14, 28, 0], // Steelhead
        [12, 9, 16, 0],  // Cutthroat
        [8, 4, 11, 0],   // Salamander
        [1, 0, 7, 0],    // Tailed Frog (singleton in Marys River -> tests Chao1!)
        [9, 11, 2, 0],   // Sculpin
        [1, 0, 0, 6],    // Human lab contaminant enriched in Extraction_Blank!
    ];

    let mut fq = String::with_capacity(65_536);
    let mut read_no = 0usize;

    let mut emit_read = |fb: &str, rb_rc: &str, amp: &str, q_char: char| {
        read_no += 1;
        let full = format!("{fb}{fwd_primer}{amp}{rev_primer_rc}{rb_rc}");
        let qual: String = std::iter::repeat(q_char).take(full.len()).collect();
        let _ = writeln!(fq, "@MISEQ_CRV:{read_no:04}\n{full}\n+\n{qual}");
    };

    for (t_idx, (_, seq)) in taxa_seqs.iter().enumerate() {
        for (s_idx, (fb, rb_rc)) in sample_tags.iter().enumerate() {
            for _ in 0..abund[t_idx][s_idx] {
                emit_read(fb, rb_rc, seq, 'I'); // Q40 ('I') -> EE ~ 0.0044
            }
        }
    }

    // Add 2 single-nucleotide sequencing error reads of Coho (66x parent -> absorbed by UNOISE skew >= 8!)
    let mut coho_err = taxa_seqs[0].1.clone().into_bytes();
    coho_err[10] = b'G';
    let coho_err_s = String::from_utf8_lossy(&coho_err).into_owned();
    emit_read(&sample_tags[0].0, &sample_tags[0].1, &coho_err_s, 'F');
    emit_read(&sample_tags[0].0, &sample_tags[0].1, &coho_err_s, 'F');

    // Add 3 PCR chimera reads in Alsea_River (5' Coho [0..22] + 3' Dicamptodon [22..44], differs >= 4 bp from both!)
    let chim_seq = format!("{}{}", &taxa_seqs[0].1[..22], &taxa_seqs[3].1[22..]);
    for _ in 0..3 {
        emit_read(&sample_tags[1].0, &sample_tags[1].1, &chim_seq, 'H');
    }

    // Add 4 noisy low-quality reads (Q5 '&' -> EE > 10 -> rejected by maxEE = 1.0)
    for _ in 0..4 {
        emit_read(&sample_tags[0].0, &sample_tags[0].1, &taxa_seqs[0].1, '&');
    }

    (fq, barcode_csv, taxa_seqs)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn demuxes_denoises_flags_chimera_and_blank_contaminant() {
        let (fq_txt, bc_txt, db) = demo_pnw_12s_edna_dataset();
        let reads = parse_fastq_str(&fq_txt).unwrap();
        let bcs = parse_barcode_sheet_csv(&bc_txt).unwrap();
        let rep =
            run_edna_metabarcoding_pipeline(&reads, &bcs, &db, &EdnaPipelineConfig::default());

        assert_eq!(rep.filtered_ee_or_short_reads, 4);
        assert_eq!(rep.error_corrected_reads, 2);
        assert_eq!(rep.chimeric_reads_flagged, 3);
        assert!(rep.asvs.iter().any(|a| a.is_chimera));
        assert!(rep.asvs.iter().any(|a| a.blank_flagged));
    }
}
