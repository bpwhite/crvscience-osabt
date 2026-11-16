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

//! Zero-lock-in plasmid and sequencing file I/O for `osabt plasmid`:
//! - GenBank flat-file (`.gb` / `.gbk`) reader and writer (preserving circular/linear topology & features)
//! - Clean-room SnapGene binary (`.dna`) packet reader (`0x09` cookie, `0x00` DNA, `0x06` Notes, `0x0A` Features, `0x05` Primers)
//! - Applied Biosystems ABIF (`.ab1`) Sanger chromatogram trace reader (`PBAS` base calls & `PCON` Phred scores)
//! - 6-frame circular Open Reading Frame (`ORF`) scanner, Standard Genetic Code translator, and Sanger trace aligner

use crate::gibson::Strand;
use crate::primer_design::parse_fasta_or_raw;
use crate::{normalize_dna, reverse_complement, OligoError};

/// Annotated sequence feature on a linear or circular plasmid construct.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlasmidFeature {
    pub name: String,
    pub feature_type: String,
    /// 0-based inclusive start coordinate on the top strand.
    /// For circular origin-wrapping features, `start_0 > end_0` (`start_0..len` joined to `0..=end_0`).
    pub start_0: usize,
    /// 0-based inclusive end coordinate on the top strand.
    pub end_0: usize,
    pub strand: Strand,
    pub color_hex: String,
    pub qualifiers: Vec<(String, String)>,
}

impl PlasmidFeature {
    /// Returns the feature's span length in base pairs on a plasmid of `plasmid_len` bp.
    pub fn span_bp(&self, plasmid_len: usize) -> usize {
        if plasmid_len == 0 {
            return 0;
        }
        if self.end_0 >= self.start_0 {
            self.end_0 - self.start_0 + 1
        } else {
            (plasmid_len - self.start_0) + self.end_0 + 1
        }
    }
}

/// Complete annotated vector/construct record shared across `osabt plasmid` CLI and `egui` GUI.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlasmidRecord {
    pub name: String,
    pub definition: String,
    pub accession: String,
    pub circular: bool,
    pub sequence: String,
    pub features: Vec<PlasmidFeature>,
}

impl PlasmidRecord {
    /// Builds a ready-to-explore 840-bp synthetic circular expression demo construct (`pCRV-DemoGFP`)
    /// containing a promoter, RBS, eGFP reporter mini-CDS, poly-His tag, terminator, ori, and unique MCS sites.
    pub fn demo_pcrv_gfp() -> Self {
        // Contains unique EcoRI (GAATTC), BamHI (GGATCC), NcoI (CCATGG), XhoI (CTCGAG), HindIII (AAGCTT), NotI (GCGGCCGC), BsaI (GGTCTC)
        let promoter_rbs =
            "TTGACAGCTAGCTCAGTCCTAGGTATAATGCTAGCGAATTCAAAGGAGGTGATCCATGG";
        let mini_gfp_cds =
            "GTGAGCAAGGGCGAGGAGCTGTTCACCGGGGTGGTGCCCATCCTGGTCGAGCTGGACGGCGACGTAAACGGCCACAAGTTCAGCGTGTCCGGCGAGGGCGAGGGCGATGCCACCTACGGCAAGCTGACCCTGAAGTTCATCTGCACCACCGGCAAGCTGCCCGTGCCCTGGCCCACCCTCGTGACCACCCTGACCTACGGCGTGCAGTGCTTCAGCCGCTACCCCGACCACATGAAGCAGCACGACTTCTTCAAGTCCGCCATGCCCGAAGGCTACGTCCAGGAGCGCACCATCTTCTTCAAGGACGACGGCAACTACAAGACCCGCGCCGAGGTGAAGTTCGAGGGCGACACCCTGGTGAACCGCATCGAGCTGAAGGGCATCGACTTCAAGGAGGACGGCAACATCCTGGGGCACAAGCTGGAGTACAACTACAACAGCCACAACGTCTATATCATGGCCGACAAGCAGAAGAACGGCATCAAGGTGAACTTCAAGATCCGCCACAACATCGAGGACGGCAGCGTGCAGCTCGCCGACCACTACCAGCAGAACACCCCCATCGGCGACGGCCCCGTGCTGCTGCCCGACAACCACTACCTGAGCACCCAGTCCGCCCTGAGCAAAGACCCCAACGAGAAGCGCGATCACATGGTCCTGCTGGAGTTCGTGACCGCCGCCGGGATCACTCTCGGCATGGACGAGCTGTACAAG";
        let his_mcs_term_ori =
            "CTCGAGCACCACCACCACCACCACTGAGGATCCAAGCTTGAGACCGCGGCCGCTACTCAGAGCCAGGCATCAAATAAAACGAAAGGCTCAGTCGAAAGACTGGGCCTTTCGTTTTATCTGTTGTTTGTCGGTGAACGCTCTCTACTAGAGTCACACTGGCTCACCTTCGGGTGGGCCTTTCTGCGTTTATAGGTACCTGACTGCGTTAGCAATTTAACTGTGATAAACTACCGCATTAAAGCTTATCGATGATAAGCTGTCAAACATGAGAATT";

        let seq = format!("{promoter_rbs}{mini_gfp_cds}{his_mcs_term_ori}");
        let p_end = promoter_rbs.len() - 1;
        let cds_start = promoter_rbs.len() - 3; // ATG inside CCATGG (NcoI)
        let cds_end = promoter_rbs.len() + mini_gfp_cds.len() + 24 - 1; // through 6xHis + TGA stop

        Self {
            name: "pCRV_DemoGFP".to_string(),
            definition: "CRV Science Network synthetic GFP cloning & expression demo vector."
                .to_string(),
            accession: "CRV0001".to_string(),
            circular: true,
            sequence: seq.clone(),
            features: vec![
                PlasmidFeature {
                    name: "J23100_Promoter".to_string(),
                    feature_type: "promoter".to_string(),
                    start_0: 0,
                    end_0: 34,
                    strand: Strand::Forward,
                    color_hex: "#f59e0b".to_string(),
                    qualifiers: vec![("label".to_string(), "J23100_Promoter".to_string())],
                },
                PlasmidFeature {
                    name: "RBS_B0034".to_string(),
                    feature_type: "RBS".to_string(),
                    start_0: 41,
                    end_0: p_end,
                    strand: Strand::Forward,
                    color_hex: "#a855f7".to_string(),
                    qualifiers: vec![("label".to_string(), "RBS_B0034".to_string())],
                },
                PlasmidFeature {
                    name: "eGFP_6xHis".to_string(),
                    feature_type: "CDS".to_string(),
                    start_0: cds_start,
                    end_0: cds_end,
                    strand: Strand::Forward,
                    color_hex: "#22c55e".to_string(),
                    qualifiers: vec![
                        ("gene".to_string(), "eGFP".to_string()),
                        ("product".to_string(), "Enhanced GFP + C-term 6xHis".to_string()),
                    ],
                },
                PlasmidFeature {
                    name: "rrnB_T1_Term".to_string(),
                    feature_type: "terminator".to_string(),
                    start_0: cds_end + 25,
                    end_0: cds_end + 150,
                    strand: Strand::Forward,
                    color_hex: "#ef4444".to_string(),
                    qualifiers: vec![("label".to_string(), "rrnB_T1_Term".to_string())],
                },
                PlasmidFeature {
                    name: "Minimal_ori".to_string(),
                    feature_type: "rep_origin".to_string(),
                    start_0: seq.len() - 95,
                    end_0: seq.len() - 8,
                    strand: Strand::Reverse,
                    color_hex: "#3b82f6".to_string(),
                    qualifiers: vec![("label".to_string(), "Minimal_ori".to_string())],
                },
            ],
        }
    }
}

/// Open Reading Frame discovered on a 6-frame scan.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OrfHit {
    pub frame: i8,
    pub strand: Strand,
    pub start_0: usize,
    pub end_0: usize,
    pub length_bp: usize,
    pub aa_sequence: String,
}

/// Extracted Sanger chromatogram base calls and Phred quality values from an ABIF `.ab1` file.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SangerTrace {
    pub sample_name: String,
    pub sequence: String,
    pub qualities: Vec<u8>,
}

/// Base-level discrepancy discovered when aligning a Sanger/FASTA sequencing read to a reference plasmid.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TraceVariant {
    pub plasmid_pos_0: usize,
    pub ref_base: char,
    pub read_base: char,
    pub read_phred_q: u8,
    pub classification: String,
}

/// Alignment summary mapping a verification sequencing trace onto a circular/linear plasmid.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TraceAlignmentResult {
    pub strand: Strand,
    pub plasmid_start_0: usize,
    pub plasmid_end_0: usize,
    pub aligned_len_bp: usize,
    pub identity_pct: f64,
    pub mean_phred_q: f64,
    pub variants: Vec<TraceVariant>,
}

fn default_feature_color(ftype: &str) -> &'static str {
    match ftype {
        "CDS" | "gene" => "#22c55e",
        "promoter" => "#f59e0b",
        "terminator" => "#ef4444",
        "rep_origin" => "#3b82f6",
        "primer_bind" => "#ec4899",
        "RBS" => "#a855f7",
        _ => "#14b8a6",
    }
}

/// Auto-detects SnapGene binary `.dna`, GenBank `.gb`/`.gbk`, or FASTA/raw DNA from raw file bytes.
pub fn load_plasmid_bytes(filename_hint: &str, bytes: &[u8]) -> Result<PlasmidRecord, OligoError> {
    if bytes.is_empty() {
        return Err(OligoError::EmptySequence);
    }
    // Check SnapGene magic packet: byte 0 == 0x09 and bytes 5..13 == b"SnapGene"
    if bytes.len() >= 19 && bytes[0] == 0x09 && &bytes[5..13] == b"SnapGene" {
        return parse_snapgene_dna(bytes);
    }

    let text = std::str::from_utf8(bytes).map_err(|_| OligoError::InvalidBase {
        byte: '?',
        pos: 0,
    })?;
    let trimmed = text.trim_start();
    if trimmed.starts_with("LOCUS") || filename_hint.ends_with(".gb") || filename_hint.ends_with(".gbk") {
        if trimmed.starts_with("LOCUS") {
            return parse_genbank(text);
        }
    }

    let recs = parse_fasta_or_raw(text)?;
    let (name, sequence) = recs.into_iter().next().ok_or(OligoError::EmptySequence)?;
    Ok(PlasmidRecord {
        name: name.clone(),
        definition: format!("Imported sequence {name}"),
        accession: name,
        circular: true,
        sequence,
        features: Vec::new(),
    })
}

/// Parses a standard GenBank (`.gb` / `.gbk`) flat-file into a `PlasmidRecord`.
pub fn parse_genbank(text: &str) -> Result<PlasmidRecord, OligoError> {
    let mut name = String::from("Unnamed_Vector");
    let mut definition = String::new();
    let mut accession = String::from("UNKNOWN");
    let mut circular = true;
    let mut features: Vec<PlasmidFeature> = Vec::new();
    let mut seq_buf = String::new();

    enum Section {
        Header,
        Features,
        Origin,
    }
    let mut sec = Section::Header;

    for raw_line in text.lines() {
        if raw_line.starts_with("//") {
            break;
        }
        if raw_line.starts_with("LOCUS") {
            let tokens: Vec<&str> = raw_line.split_whitespace().collect();
            if tokens.len() >= 2 {
                name = tokens[1].to_string();
                accession = name.clone();
            }
            let lower = raw_line.to_ascii_lowercase();
            if lower.contains("linear") {
                circular = false;
            } else if lower.contains("circular") {
                circular = true;
            }
            continue;
        }
        if let Some(rest) = raw_line.strip_prefix("DEFINITION") {
            definition = rest.trim().to_string();
            continue;
        }
        if let Some(rest) = raw_line.strip_prefix("ACCESSION") {
            if let Some(tok) = rest.split_whitespace().next() {
                accession = tok.to_string();
            }
            continue;
        }
        if raw_line.starts_with("FEATURES") {
            sec = Section::Features;
            continue;
        }
        if raw_line.starts_with("ORIGIN") {
            sec = Section::Origin;
            continue;
        }

        match sec {
            Section::Header => {}
            Section::Features => {
                if raw_line.len() > 5 {
                    let head = raw_line[..21.min(raw_line.len())].trim();
                    let body = if raw_line.len() > 21 {
                        raw_line[21..].trim()
                    } else {
                        ""
                    };
                    if !head.is_empty() && !raw_line.starts_with("      ") {
                        let (strand, start_0, end_0) = parse_gbk_location(body);
                        let color = default_feature_color(head).to_string();
                        features.push(PlasmidFeature {
                            name: head.to_string(),
                            feature_type: head.to_string(),
                            start_0,
                            end_0,
                            strand,
                            color_hex: color,
                            qualifiers: Vec::new(),
                        });
                    } else if let Some(feat) = features.last_mut() {
                        if let Some(q) = body.strip_prefix('/') {
                            if let Some((k, v)) = q.split_once('=') {
                                let clean_v = v.trim_matches('"').to_string();
                                if matches!(k, "label" | "gene" | "product" | "note")
                                    && (feat.name == feat.feature_type || k == "label" || k == "gene")
                                {
                                    feat.name = clean_v.clone();
                                }
                                if (k.eq_ignore_ascii_case("ApEinfo_fwdcolor")
                                    || k.eq_ignore_ascii_case("color"))
                                    && clean_v.starts_with('#')
                                {
                                    feat.color_hex = clean_v.clone();
                                }
                                feat.qualifiers.push((k.to_string(), clean_v));
                            }
                        }
                    }
                }
            }
            Section::Origin => {
                for ch in raw_line.chars() {
                    if ch.is_ascii_alphabetic() {
                        seq_buf.push(ch);
                    }
                }
            }
        }
    }

    let sequence = normalize_dna(&seq_buf)?;
    Ok(PlasmidRecord {
        name,
        definition,
        accession,
        circular,
        sequence,
        features,
    })
}

fn parse_gbk_location(loc: &str) -> (Strand, usize, usize) {
    let mut s = loc.trim();
    let mut strand = Strand::Forward;
    if let Some(inner) = s.strip_prefix("complement(").and_then(|x| x.strip_suffix(')')) {
        strand = Strand::Reverse;
        s = inner.trim();
    }
    if let Some(inner) = s.strip_prefix("join(").and_then(|x| x.strip_suffix(')')) {
        // Origin-spanning join(a..N,1..b) -> start = a-1, end = b-1
        let segs: Vec<&str> = inner.split(',').collect();
        if segs.len() >= 2 {
            let (_, a, _) = parse_gbk_location(segs[0]);
            let (_, _, b) = parse_gbk_location(segs[segs.len() - 1]);
            return (strand, a, b);
        }
    }
    let clean: String = s
        .chars()
        .filter(|c| c.is_ascii_digit() || *c == '.')
        .collect();
    if let Some((a_str, b_str)) = clean.split_once("..") {
        let a = a_str.parse::<usize>().unwrap_or(1).saturating_sub(1);
        let b = b_str.parse::<usize>().unwrap_or(1).saturating_sub(1);
        (strand, a, b)
    } else {
        let pos = clean.parse::<usize>().unwrap_or(1).saturating_sub(1);
        (strand, pos, pos)
    }
}

/// Serializes a `PlasmidRecord` to standard GenBank (`.gbk`) flat-file format.
pub fn write_genbank(rec: &PlasmidRecord) -> String {
    let topo = if rec.circular { "circular" } else { "linear  " };
    let mut out = format!(
        "LOCUS       {:<16} {:>7} bp ds-DNA     {} SYN 11-NOV-2026\n\
         DEFINITION  {}\n\
         ACCESSION   {}\n\
         KEYWORDS    OSABT; CRV Science Network.\n\
         FEATURES             Location/Qualifiers\n",
        rec.name,
        rec.sequence.len(),
        topo,
        if rec.definition.is_empty() {
            &rec.name
        } else {
            &rec.definition
        },
        rec.accession
    );

    let n = rec.sequence.len();
    for f in &rec.features {
        let span = if f.end_0 >= f.start_0 || !rec.circular {
            format!("{}..{}", f.start_0 + 1, f.end_0 + 1)
        } else {
            format!("join({}..{},1..{})", f.start_0 + 1, n, f.end_0 + 1)
        };
        let loc_str = match f.strand {
            Strand::Forward => span,
            Strand::Reverse => format!("complement({span})"),
        };
        out.push_str(&format!("     {:<16}{}\n", f.feature_type, loc_str));
        out.push_str(&format!("                     /label=\"{}\"\n", f.name));
        out.push_str(&format!(
            "                     /ApEinfo_fwdcolor=\"{}\"\n",
            f.color_hex
        ));
        for (k, v) in &f.qualifiers {
            if k != "label" && k != "ApEinfo_fwdcolor" {
                out.push_str(&format!("                     /{k}=\"{v}\"\n"));
            }
        }
    }

    out.push_str("ORIGIN\n");
    let lower_seq = rec.sequence.to_ascii_lowercase();
    let bytes = lower_seq.as_bytes();
    for (row_idx, chunk) in bytes.chunks(60).enumerate() {
        out.push_str(&format!("{:>9}", row_idx * 60 + 1));
        for sub in chunk.chunks(10) {
            out.push(' ');
            out.push_str(std::str::from_utf8(sub).unwrap_or(""));
        }
        out.push('\n');
    }
    out.push_str("//\n");
    out
}

/// Clean-room SnapGene binary (`.dna`) packet reader matching the TLV packet layout:
/// `1-byte type` + `4-byte big-endian length` + `payload`.
pub fn parse_snapgene_dna(bytes: &[u8]) -> Result<PlasmidRecord, OligoError> {
    if bytes.len() < 19 || bytes[0] != 0x09 {
        return Err(OligoError::InvalidConcentration(
            "File does not begin with a valid SnapGene 0x09 cookie packet.",
        ));
    }
    let cookie_len = u32::from_be_bytes([bytes[1], bytes[2], bytes[3], bytes[4]]) as usize;
    if 5 + cookie_len > bytes.len() || &bytes[5..13] != b"SnapGene" {
        return Err(OligoError::InvalidConcentration(
            "Corrupt SnapGene magic cookie header.",
        ));
    }

    let mut cursor = 5 + cookie_len;
    let mut name = String::from("SnapGene_Import");
    let mut definition = String::new();
    let mut accession = String::from("SG_IMPORT");
    let mut circular = true;
    let mut sequence = String::new();
    let mut raw_feature_xmls: Vec<(u8, String)> = Vec::new();

    while cursor + 5 <= bytes.len() {
        let pkt_type = bytes[cursor];
        let pkt_len = u32::from_be_bytes([
            bytes[cursor + 1],
            bytes[cursor + 2],
            bytes[cursor + 3],
            bytes[cursor + 4],
        ]) as usize;
        cursor += 5;
        if cursor + pkt_len > bytes.len() {
            break;
        }
        let payload = &bytes[cursor..(cursor + pkt_len)];
        cursor += pkt_len;

        match pkt_type {
            0x00 if !payload.is_empty() => {
                let flags = payload[0];
                circular = (flags & 0x01) != 0;
                if let Ok(raw_dna) = std::str::from_utf8(&payload[1..]) {
                    sequence = normalize_dna(raw_dna)?;
                }
            }
            0x06 => {
                if let Ok(xml) = std::str::from_utf8(payload) {
                    if let Some(acc) = extract_xml_tag_text(xml, "AccessionNumber") {
                        if !acc.is_empty() {
                            accession = acc.clone();
                            name = acc;
                        }
                    }
                    if let Some(comm) = extract_xml_tag_text(xml, "Comments") {
                        definition = comm;
                    }
                }
            }
            0x05 | 0x0A => {
                if let Ok(xml) = std::str::from_utf8(payload) {
                    raw_feature_xmls.push((pkt_type, xml.to_string()));
                }
            }
            _ => {}
        }
    }

    if sequence.is_empty() {
        return Err(OligoError::EmptySequence);
    }

    let mut features = Vec::new();
    for (pkt_type, xml) in raw_feature_xmls {
        if pkt_type == 0x0A {
            parse_snapgene_features_xml(&xml, sequence.len(), &mut features);
        } else if pkt_type == 0x05 {
            parse_snapgene_primers_xml(&xml, sequence.len(), &mut features);
        }
    }

    Ok(PlasmidRecord {
        name,
        definition,
        accession,
        circular,
        sequence,
        features,
    })
}

fn extract_xml_attr<'a>(tag_header: &'a str, attr: &str) -> Option<&'a str> {
    let pat = format!("{attr}=\"");
    let start = tag_header.find(&pat)? + pat.len();
    let rest = &tag_header[start..];
    let end = rest.find('"')?;
    Some(&rest[..end])
}

fn extract_xml_tag_text(xml: &str, tag: &str) -> Option<String> {
    let open = format!("<{tag}>");
    let close = format!("</{tag}>");
    let s = xml.find(&open)? + open.len();
    let e = xml[s..].find(&close)? + s;
    Some(xml[s..e].trim().to_string())
}

fn parse_snapgene_features_xml(xml: &str, seq_len: usize, out: &mut Vec<PlasmidFeature>) {
    let mut rest = xml;
    while let Some(f_idx) = rest.find("<Feature ") {
        rest = &rest[f_idx..];
        let Some(end_feat) = rest.find("</Feature>") else {
            break;
        };
        let block = &rest[..end_feat];
        let header_end = block.find('>').unwrap_or(block.len());
        let header = &block[..header_end];

        let fname = extract_xml_attr(header, "name")
            .unwrap_or("Feature")
            .to_string();
        let ftype = extract_xml_attr(header, "type")
            .unwrap_or("misc_feature")
            .to_string();
        let dir = extract_xml_attr(header, "directionality").unwrap_or("1");
        let strand = if dir == "2" {
            Strand::Reverse
        } else {
            Strand::Forward
        };

        if let Some(seg_idx) = block.find("<Segment ") {
            let seg = &block[seg_idx..];
            let seg_end = seg.find('>').unwrap_or(seg.len());
            let seg_hdr = &seg[..seg_end];
            let color = extract_xml_attr(seg_hdr, "color")
                .map(|s| s.to_string())
                .unwrap_or_else(|| default_feature_color(&ftype).to_string());
            if let Some(rng) = extract_xml_attr(seg_hdr, "range") {
                if let Some((a_s, b_s)) = rng.split_once('-') {
                    let a = a_s.parse::<usize>().unwrap_or(1).saturating_sub(1) % seq_len.max(1);
                    let b = b_s.parse::<usize>().unwrap_or(1).saturating_sub(1) % seq_len.max(1);
                    out.push(PlasmidFeature {
                        name: fname,
                        feature_type: ftype,
                        start_0: a,
                        end_0: b,
                        strand,
                        color_hex: color,
                        qualifiers: Vec::new(),
                    });
                }
            }
        }
        rest = &rest[(end_feat + "</Feature>".len())..];
    }
}

fn parse_snapgene_primers_xml(xml: &str, seq_len: usize, out: &mut Vec<PlasmidFeature>) {
    let mut rest = xml;
    while let Some(p_idx) = rest.find("<Primer ") {
        rest = &rest[p_idx..];
        let Some(end_p) = rest.find("</Primer>") else {
            break;
        };
        let block = &rest[..end_p];
        let hdr_end = block.find('>').unwrap_or(block.len());
        let pname = extract_xml_attr(&block[..hdr_end], "name")
            .unwrap_or("Primer")
            .to_string();

        if let Some(bs_idx) = block.find("<BindingSite ") {
            let bs = &block[bs_idx..];
            let bs_end = bs.find('>').unwrap_or(bs.len());
            let bs_hdr = &bs[..bs_end];
            let strand = if extract_xml_attr(bs_hdr, "boundStrand") == Some("1") {
                Strand::Reverse
            } else {
                Strand::Forward
            };
            if let Some(loc) = extract_xml_attr(bs_hdr, "location") {
                if let Some((a_s, b_s)) = loc.split_once('-') {
                    let a = a_s.parse::<usize>().unwrap_or(0) % seq_len.max(1);
                    let b = b_s.parse::<usize>().unwrap_or(0) % seq_len.max(1);
                    out.push(PlasmidFeature {
                        name: pname,
                        feature_type: "primer_bind".to_string(),
                        start_0: a,
                        end_0: b,
                        strand,
                        color_hex: "#ec4899".to_string(),
                        qualifiers: Vec::new(),
                    });
                }
            }
        }
        rest = &rest[(end_p + "</Primer>".len())..];
    }
}

/// Parses an Applied Biosystems ABIF (`.ab1`) Sanger chromatogram file, extracting called bases
/// (`PBAS` tag) and Phred quality scores (`PCON` tag).
pub fn parse_ab1_trace(bytes: &[u8]) -> Result<SangerTrace, OligoError> {
    if bytes.len() < 128 || &bytes[0..4] != b"ABIF" {
        return Err(OligoError::InvalidConcentration(
            "Input file is not a valid ABIF (.ab1) chromatogram.",
        ));
    }

    let num_entries =
        u32::from_be_bytes([bytes[18], bytes[19], bytes[20], bytes[21]]) as usize;
    let dir_offset =
        u32::from_be_bytes([bytes[26], bytes[27], bytes[28], bytes[29]]) as usize;

    let mut pbas_1: Option<String> = None;
    let mut pbas_2: Option<String> = None;
    let mut pcon_1: Option<Vec<u8>> = None;
    let mut pcon_2: Option<Vec<u8>> = None;
    let mut sample_name = String::from("AB1_Sample");

    for i in 0..num_entries {
        let off = dir_offset + i * 28;
        if off + 28 > bytes.len() {
            break;
        }
        let tag = &bytes[off..(off + 4)];
        let tag_num = u32::from_be_bytes([
            bytes[off + 4],
            bytes[off + 5],
            bytes[off + 6],
            bytes[off + 7],
        ]);
        let data_size = u32::from_be_bytes([
            bytes[off + 16],
            bytes[off + 17],
            bytes[off + 18],
            bytes[off + 19],
        ]) as usize;

        let slice: &[u8] = if data_size <= 4 {
            &bytes[(off + 20)..(off + 20 + data_size)]
        } else {
            let d_off = u32::from_be_bytes([
                bytes[off + 20],
                bytes[off + 21],
                bytes[off + 22],
                bytes[off + 23],
            ]) as usize;
            if d_off + data_size > bytes.len() {
                continue;
            }
            &bytes[d_off..(d_off + data_size)]
        };

        match tag {
            b"PBAS" => {
                if let Ok(s) = std::str::from_utf8(slice) {
                    let clean: String = s
                        .chars()
                        .filter(|c| c.is_ascii_alphabetic())
                        .collect();
                    if let Ok(norm) = normalize_dna(&clean) {
                        if tag_num == 2 {
                            pbas_2 = Some(norm);
                        } else if tag_num == 1 {
                            pbas_1 = Some(norm);
                        }
                    }
                }
            }
            b"PCON" => {
                if tag_num == 2 {
                    pcon_2 = Some(slice.to_vec());
                } else if tag_num == 1 {
                    pcon_1 = Some(slice.to_vec());
                }
            }
            b"SMPL" if !slice.is_empty() => {
                // Pascal string: byte 0 is length
                let plen = (slice[0] as usize).min(slice.len() - 1);
                if let Ok(name) = std::str::from_utf8(&slice[1..(1 + plen)]) {
                    if !name.trim().is_empty() {
                        sample_name = name.trim().to_string();
                    }
                }
            }
            _ => {}
        }
    }

    let sequence = pbas_2.or(pbas_1).ok_or(OligoError::EmptySequence)?;
    let mut qualities = pcon_2
        .or(pcon_1)
        .unwrap_or_else(|| vec![40u8; sequence.len()]);
    qualities.resize(sequence.len(), 30);

    Ok(SangerTrace {
        sample_name,
        sequence,
        qualities,
    })
}

/// Translates a `5' -> 3'` DNA coding sequence using NCBI Standard Genetic Code (Table 1).
pub fn translate_dna(cds_5to3: &str) -> String {
    let b = cds_5to3.as_bytes();
    let mut aa = String::with_capacity(b.len() / 3);
    for codon in b.chunks_exact(3) {
        let res = match [
            codon[0].to_ascii_uppercase(),
            codon[1].to_ascii_uppercase(),
            codon[2].to_ascii_uppercase(),
        ] {
            [b'T', b'T', b'T'] | [b'T', b'T', b'C'] => 'F',
            [b'T', b'T', b'A'] | [b'T', b'T', b'G'] | [b'C', b'T', _] => 'L',
            [b'A', b'T', b'T'] | [b'A', b'T', b'C'] | [b'A', b'T', b'A'] => 'I',
            [b'A', b'T', b'G'] => 'M',
            [b'G', b'T', _] => 'V',
            [b'T', b'C', _] | [b'A', b'G', b'T'] | [b'A', b'G', b'C'] => 'S',
            [b'C', b'C', _] => 'P',
            [b'A', b'C', _] => 'T',
            [b'G', b'C', _] => 'A',
            [b'T', b'A', b'T'] | [b'T', b'A', b'C'] => 'Y',
            [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A'] => '*',
            [b'C', b'A', b'T'] | [b'C', b'A', b'C'] => 'H',
            [b'C', b'A', b'A'] | [b'C', b'A', b'G'] => 'Q',
            [b'A', b'A', b'T'] | [b'A', b'A', b'C'] => 'N',
            [b'A', b'A', b'A'] | [b'A', b'A', b'G'] => 'K',
            [b'G', b'A', b'T'] | [b'G', b'A', b'C'] => 'D',
            [b'G', b'A', b'A'] | [b'G', b'A', b'G'] => 'E',
            [b'T', b'G', b'T'] | [b'T', b'G', b'C'] => 'C',
            [b'T', b'G', b'G'] => 'W',
            [b'C', b'G', _] | [b'A', b'G', b'A'] | [b'A', b'G', b'G'] => 'R',
            [b'G', b'G', _] => 'G',
            _ => 'X',
        };
        aa.push(res);
    }
    aa
}

/// Scans all 6 reading frames (`+1..+3`, `-1..-3`) across a linear or circular plasmid for
/// `ATG -> Stop` Open Reading Frames of at least `min_aa` amino acids.
pub fn find_orfs(rec: &PlasmidRecord, min_aa: usize) -> Vec<OrfHit> {
    let n = rec.sequence.len();
    let min_bp = min_aa.max(20) * 3;
    if n < min_bp {
        return Vec::new();
    }

    let rc = reverse_complement(&rec.sequence);
    let pad = if rec.circular { n.min(1500) } else { 0 };
    let fwd_track = format!("{}{}", rec.sequence, &rec.sequence[..pad]);
    let rev_track = format!("{}{}", rc, &rc[..pad]);

    let mut orfs = Vec::new();
    for (strand, track) in [
        (Strand::Forward, fwd_track.as_str()),
        (Strand::Reverse, rev_track.as_str()),
    ] {
        let tb = track.as_bytes();
        for frame_off in 0..3usize {
            let mut i = frame_off;
            while i + 3 <= n {
                if &tb[i..(i + 3)] == b"ATG" {
                    let mut j = i + 3;
                    while j + 3 <= tb.len() && (j - i) < n {
                        let c = &tb[j..(j + 3)];
                        if c == b"TAA" || c == b"TAG" || c == b"TGA" {
                            let len_bp = (j + 3) - i;
                            if len_bp >= min_bp {
                                let aa = translate_dna(&track[i..(j + 3)]);
                                let (start_0, end_0) = match strand {
                                    Strand::Forward => (i % n, (j + 2) % n),
                                    Strand::Reverse => {
                                        ((n - 1) - (i % n), (n - 1) - ((j + 2) % n))
                                    }
                                };
                                let signed_frame = match strand {
                                    Strand::Forward => (frame_off as i8) + 1,
                                    Strand::Reverse => -((frame_off as i8) + 1),
                                };
                                orfs.push(OrfHit {
                                    frame: signed_frame,
                                    strand,
                                    start_0,
                                    end_0,
                                    length_bp: len_bp,
                                    aa_sequence: aa,
                                });
                            }
                            break;
                        }
                        j += 3;
                    }
                }
                i += 3;
            }
        }
    }
    orfs.sort_by(|a, b| b.length_bp.cmp(&a.length_bp));
    orfs
}

/// Aligns a high-quality trimmed Sanger `.ab1` or FASTA read onto a circular/linear plasmid,
/// locating the best k-mer seed anchor and classifying high-confidence vs. low-Phred discrepancies.
pub fn align_trace_to_plasmid(
    rec: &PlasmidRecord,
    trace: &SangerTrace,
    min_phred_q: u8,
) -> Option<TraceAlignmentResult> {
    let n = rec.sequence.len();
    if n < 20 || trace.sequence.len() < 20 {
        return None;
    }

    // Quality-trim noisy 5' and 3' ends of Sanger read
    let q = &trace.qualities;
    let trim_s = q
        .iter()
        .position(|&v| v >= min_phred_q)
        .unwrap_or(0)
        .min(trace.sequence.len().saturating_sub(20));
    let trim_e = q
        .iter()
        .rposition(|&v| v >= min_phred_q)
        .map(|idx| idx + 1)
        .unwrap_or(trace.sequence.len())
        .max(trim_s + 20)
        .min(trace.sequence.len());

    let read_f = &trace.sequence[trim_s..trim_e];
    let qual_f = &trace.qualities[trim_s..trim_e];
    let read_r = reverse_complement(read_f);
    let qual_r: Vec<u8> = qual_f.iter().copied().rev().collect();

    let wrap = if rec.circular { read_f.len().min(n) } else { 0 };
    let plasmid_track = format!("{}{}", rec.sequence, &rec.sequence[..wrap]);
    let pb = plasmid_track.as_bytes();

    let mut best: Option<(usize, Strand, usize, &[u8], &[u8])> = None;

    for (strand, r_str, r_qual) in [
        (Strand::Forward, read_f, qual_f),
        (Strand::Reverse, read_r.as_str(), qual_r.as_slice()),
    ] {
        let rb = r_str.as_bytes();
        if rb.len() > pb.len() {
            continue;
        }
        let max_offset = if rec.circular {
            n
        } else {
            n.saturating_sub(rb.len()) + 1
        };

        // Fast 15-mer seed probe at 25%, 50%, 75% of read to find candidate shift
        let seed_len = 15usize.min(rb.len());
        let mut candidate_shifts = Vec::new();
        for frac in [rb.len() / 4, rb.len() / 2, (3 * rb.len()) / 4] {
            let s_pos = frac.min(rb.len() - seed_len);
            let seed = &rb[s_pos..(s_pos + seed_len)];
            for t_pos in s_pos..(max_offset + s_pos).min(pb.len().saturating_sub(seed_len) + 1) {
                if &pb[t_pos..(t_pos + seed_len)] == seed {
                    let shift = t_pos - s_pos;
                    if shift < max_offset && !candidate_shifts.contains(&shift) {
                        candidate_shifts.push(shift);
                    }
                }
            }
        }

        for shift in candidate_shifts {
            if shift + rb.len() > pb.len() {
                continue;
            }
            let matches = rb
                .iter()
                .zip(&pb[shift..(shift + rb.len())])
                .filter(|(a, b)| a == b)
                .count();
            if best.as_ref().map_or(true, |(bm, _, _, _, _)| matches > *bm) {
                best = Some((matches, strand, shift, rb, r_qual));
            }
        }
    }

    let (matches, strand, shift, rb, rq) = best?;
    let aln_len = rb.len();
    let identity_pct = ((matches as f64 / aln_len as f64) * 10_000.0).round() / 100.0;
    let mean_q = rq.iter().map(|&x| x as f64).sum::<f64>() / (aln_len as f64);

    let mut variants = Vec::new();
    for i in 0..aln_len {
        let ref_b = pb[shift + i];
        let rd_b = rb[i];
        if ref_b != rd_b {
            let phred = rq[i];
            let kind = if phred < min_phred_q || rd_b == b'N' {
                "LowQ_Ambiguity"
            } else {
                "HighConf_Substitution"
            };
            variants.push(TraceVariant {
                plasmid_pos_0: (shift + i) % n,
                ref_base: ref_b as char,
                read_base: rd_b as char,
                read_phred_q: phred,
                classification: kind.to_string(),
            });
        }
    }

    Some(TraceAlignmentResult {
        strand,
        plasmid_start_0: shift % n,
        plasmid_end_0: (shift + aln_len - 1) % n,
        aligned_len_bp: aln_len,
        identity_pct,
        mean_phred_q: (mean_q * 10.0).round() / 10.0,
        variants,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_genbank_and_parses_snapgene_packets() {
        let demo = PlasmidRecord::demo_pcrv_gfp();
        let gbk = write_genbank(&demo);
        let parsed = parse_genbank(&gbk).unwrap();
        assert_eq!(parsed.name, "pCRV_DemoGFP");
        assert!(parsed.circular);
        assert_eq!(parsed.sequence, demo.sequence);
        assert_eq!(parsed.features.len(), demo.features.len());

        // Construct synthetic binary SnapGene .dna stream (0x09 Cookie + 0x00 DNA + 0x0A Feature)
        let mut sg = Vec::new();
        sg.push(0x09);
        sg.extend_from_slice(&(14u32.to_be_bytes()));
        sg.extend_from_slice(b"SnapGene\x00\x01\x00\x02\x00\x02");

        let dna = b"\x01ATGAAAGGGCCCTTTAA";
        sg.push(0x00);
        sg.extend_from_slice(&(dna.len() as u32).to_be_bytes());
        sg.extend_from_slice(dna);

        let xml = b"<Features><Feature name=\"MiniORF\" type=\"CDS\" directionality=\"1\"><Segment range=\"1-15\" color=\"#22c55e\"/></Feature></Features>";
        sg.push(0x0A);
        sg.extend_from_slice(&(xml.len() as u32).to_be_bytes());
        sg.extend_from_slice(xml);

        let sg_rec = load_plasmid_bytes("test.dna", &sg).unwrap();
        assert!(sg_rec.circular);
        assert_eq!(sg_rec.sequence, "ATGAAAGGGCCCTTTAA");
        assert_eq!(sg_rec.features[0].name, "MiniORF");
    }
}
