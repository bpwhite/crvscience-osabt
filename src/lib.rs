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

//! # Open Source Ad-hoc Biology Tools (`osabt` `v2.0.0`)
//!
//! Unified GPL-3.0 Rust library and application suite for the **CRV Science Network**
//! (`crvscience.com`), housing 20 zero-paywall wet-lab, cloning, ecological, phylogenetic,
//! epidemiological, etymological, 3D-printing, morphometric, bioassay, bioacoustics,
//! evolutionary information-theory, MIQE RT-qPCR, steady-state enzyme kinetics,
//! environmental DNA (`eDNA`) metabarcoding, wildlife occupancy / capture-recapture,
//! population genetics, CRISPR guide design, flow cytometry, survival / ecotoxicology,
//! and leaf photosynthesis / stream metabolism tools.
//!
//! ## Completed 20-Tool Suite (`v2.0.0` — Tools #1 through #20)
//! - **Tool #1 — `bench` (`osabt bench`):** [`thermo`], [`structure`], [`primer_design`], [`gibson`], [`bench_calc`], [`bench_cli`], [`bench_tui`].
//! - **Tool #2 — `plasmid` (`osabt plasmid`):** [`plasmid_io`], [`plasmid_digest`], [`plasmid_svg`], [`plasmid_gui`], [`plasmid_cli`].
//! - **Tool #3 — `dwc` (`osabt dwc`):** [`dwc_core`], [`dwc_taxonomy`], [`dwc_export`], [`dwc_tui`], [`dwc_cli`].
//! - **Tool #4 — `phylo` (`osabt phylo`):** [`phylo_core`], [`phylo_render`], [`phylo_tui`], [`phylo_cli`].
//! - **Tool #5 — `epi` (`osabt epi`):** [`epi_models`], [`epi_power`], [`epi_tui`], [`epi_cli`].
//! - **Tool #6 — `lexi` (`osabt lexi`):** [`lexi_roots`], [`lexi_engine`], [`lexi_tui`], [`lexi_cli`].
//! - **Tool #7 — `forge3d` (`osabt forge3d`):** [`forge3d_mesh`], [`forge3d_bio`], [`forge3d_gui`], [`forge3d_cli`].
//! - **Tool #8 — `morpho` (`osabt morpho`):** [`morpho_vision`], [`morpho_gpa`], [`morpho_gui`], [`morpho_cli`].
//! - **Tool #9 — `plate` (`osabt plate`):** [`plate_layout`], [`plate_curve`], [`plate_tui`], [`plate_cli`].
//! - **Tool #10 — `ecosound` (`osabt ecosound`):** [`ecosound_dsp`], [`ecosound_slicer`], [`ecosound_gui`], [`ecosound_cli`].
//! - **Tool #11 — `infoevo` (`osabt infoevo`):** [`infoevo_core`], [`infoevo_svg`], [`infoevo_tui`], [`infoevo_cli`].
//! - **Tool #12 — `qpcr` (`osabt qpcr`):** [`qpcr_core`], [`qpcr_melt`], [`qpcr_tui`], [`qpcr_cli`].
//! - **Tool #13 — `kinetics` (`osabt kinetics`):** [`kinetics_core`], [`kinetics_svg`], [`kinetics_gui`], [`kinetics_cli`].
//! - **Tool #14 — `edna` (`osabt edna`):** [`edna_fastq`], [`edna_diversity`], [`edna_tui`], [`edna_cli`].
//! - **Tool #15 — `wildlife` (`osabt wildlife`):** [`wildlife_occupancy`], [`wildlife_markrecap`], [`wildlife_tui`], [`wildlife_cli`].
//! - **Tool #16 — `popgen` (`osabt popgen`):** [`popgen_core`], [`popgen_network`], [`popgen_gui`], [`popgen_cli`].
//! - **Tool #17 — `crispr` (`osabt crispr`):** [`crispr_core`], [`crispr_offtarget`], [`crispr_tui`], [`crispr_cli`].
//! - **Tool #18 — `flow` (`osabt flow`):** [`flow_fcs`], [`flow_gate`], [`flow_gui`], [`flow_cli`].
//! - **Tool #19 — `survival` (`osabt survival`):** [`survival_core`], [`survival_tox`], [`survival_tui`], [`survival_cli`].
//! - **Tool #20 — `flux` (`osabt flux`):**
//!   - [`flux_leaf`]: Farquhar–von Caemmerer–Berry (1980) `FvCB` C3 $A\text{–}C_i$ fitter (`Vcmax`, `Jmax`, `TPU`, `Rd`).
//!   - [`flux_stream`]: Odum (1956) open-channel diel `DO` stream metabolism (`GPP`, `ER`, `NEP`, `K600`) & `.svg` plotter.
//!   - [`flux_gui`] & [`flux_cli`]: Interactive `egui` Flux Studio + `clap` CLI.

pub mod bench_calc;
pub mod crispr_core;
pub mod crispr_offtarget;
pub mod dwc_core;
pub mod dwc_export;
pub mod dwc_taxonomy;
pub mod ecosound_dsp;
pub mod ecosound_slicer;
pub mod edna_diversity;
pub mod edna_fastq;
pub mod epi_models;
pub mod epi_power;
pub mod flow_fcs;
pub mod flow_gate;
pub mod flux_leaf;
pub mod flux_stream;
pub mod forge3d_bio;
pub mod forge3d_mesh;
pub mod gibson;
pub mod infoevo_core;
pub mod infoevo_svg;
pub mod kinetics_core;
pub mod kinetics_svg;
pub mod lexi_engine;
pub mod lexi_roots;
pub mod morpho_gpa;
pub mod morpho_vision;
pub mod phylo_core;
pub mod phylo_render;
pub mod plasmid_digest;
pub mod plasmid_io;
pub mod plasmid_svg;
pub mod plate_curve;
pub mod plate_layout;
pub mod popgen_core;
pub mod popgen_network;
pub mod primer_design;
pub mod qpcr_core;
pub mod qpcr_melt;
pub mod structure;
pub mod survival_core;
pub mod survival_tox;
pub mod svg_util;
pub mod thermo;
pub mod wildlife_markrecap;
pub mod wildlife_occupancy;

#[cfg(feature = "cli")]
pub mod bench_cli;

#[cfg(feature = "cli")]
pub mod crispr_cli;

#[cfg(feature = "cli")]
pub mod dwc_cli;

#[cfg(feature = "cli")]
pub mod ecosound_cli;

#[cfg(feature = "cli")]
pub mod edna_cli;

#[cfg(feature = "cli")]
pub mod epi_cli;

#[cfg(feature = "cli")]
pub mod flow_cli;

#[cfg(feature = "cli")]
pub mod flux_cli;

#[cfg(feature = "cli")]
pub mod forge3d_cli;

#[cfg(feature = "cli")]
pub mod infoevo_cli;

#[cfg(feature = "cli")]
pub mod kinetics_cli;

#[cfg(feature = "cli")]
pub mod lexi_cli;

#[cfg(feature = "cli")]
pub mod morpho_cli;

#[cfg(feature = "cli")]
pub mod phylo_cli;

#[cfg(feature = "cli")]
pub mod plasmid_cli;

#[cfg(feature = "cli")]
pub mod plate_cli;

#[cfg(feature = "cli")]
pub mod popgen_cli;

#[cfg(feature = "cli")]
pub mod qpcr_cli;

#[cfg(feature = "cli")]
pub mod survival_cli;

#[cfg(feature = "cli")]
pub mod wildlife_cli;

#[cfg(feature = "tui")]
pub mod tui_util;

#[cfg(feature = "tui")]
pub mod bench_tui;

#[cfg(feature = "tui")]
pub mod crispr_tui;

#[cfg(feature = "tui")]
pub mod dwc_tui;

#[cfg(feature = "tui")]
pub mod edna_tui;

#[cfg(feature = "tui")]
pub mod epi_tui;

#[cfg(feature = "tui")]
pub mod infoevo_tui;

#[cfg(feature = "tui")]
pub mod lexi_tui;

#[cfg(feature = "tui")]
pub mod phylo_tui;

#[cfg(feature = "tui")]
pub mod plate_tui;

#[cfg(feature = "tui")]
pub mod qpcr_tui;

#[cfg(feature = "tui")]
pub mod survival_tui;

#[cfg(feature = "tui")]
pub mod wildlife_tui;

#[cfg(feature = "gui")]
pub mod ecosound_gui;

#[cfg(feature = "gui")]
pub mod flow_gui;

#[cfg(feature = "gui")]
pub mod flux_gui;

#[cfg(feature = "gui")]
pub mod forge3d_gui;

#[cfg(feature = "gui")]
pub mod kinetics_gui;

#[cfg(feature = "gui")]
pub mod morpho_gui;

#[cfg(feature = "gui")]
pub mod plasmid_gui;

#[cfg(feature = "gui")]
pub mod popgen_gui;

#[cfg(feature = "gui")]
pub mod workbench_gui;

pub use bench_calc::*;
pub use crispr_core::*;
pub use crispr_offtarget::*;
pub use dwc_core::*;
pub use dwc_export::*;
pub use dwc_taxonomy::*;
pub use ecosound_dsp::*;
pub use ecosound_slicer::*;
pub use edna_diversity::*;
pub use edna_fastq::*;
pub use epi_models::*;
pub use epi_power::*;
pub use flow_fcs::*;
pub use flow_gate::*;
pub use flux_leaf::*;
pub use flux_stream::*;
pub use forge3d_bio::*;
pub use forge3d_mesh::*;
pub use gibson::*;
pub use infoevo_core::*;
pub use infoevo_svg::*;
pub use kinetics_core::*;
pub use kinetics_svg::*;
pub use lexi_engine::*;
pub use lexi_roots::*;
pub use morpho_gpa::*;
pub use morpho_vision::*;
pub use phylo_core::*;
pub use phylo_render::*;
pub use plasmid_digest::*;
pub use plasmid_io::*;
pub use plasmid_svg::*;
pub use plate_curve::*;
pub use plate_layout::*;
pub use popgen_core::*;
pub use popgen_network::*;
pub use primer_design::*;
pub use qpcr_core::*;
pub use qpcr_melt::*;
pub use structure::*;
pub use survival_core::*;
pub use survival_tox::*;
pub use thermo::*;
pub use wildlife_markrecap::*;
pub use wildlife_occupancy::*;

use core::fmt;

/// Error type returned across OSABT sequence, thermodynamic, cloning, and wet-lab calculations.
#[derive(Debug, Clone, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum OligoError {
    EmptySequence,
    SequenceTooShort { len: usize, min: usize },
    InvalidBase { byte: char, pos: usize },
    AmbiguousBaseInStrictMode { byte: char, pos: usize },
    InvalidConcentration(&'static str),
    AssemblyDesignFailed { junction_index: usize, reason: String },
}

impl fmt::Display for OligoError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::EmptySequence => write!(f, "DNA sequence is empty after stripping whitespace."),
            Self::SequenceTooShort { len, min } => {
                write!(
                    f,
                    "Sequence length ({len} nt) is shorter than minimum required ({min} nt)."
                )
            }
            Self::InvalidBase { byte, pos } => {
                write!(f, "Invalid IUPAC nucleotide character '{byte}' at position {pos}.")
            }
            Self::AmbiguousBaseInStrictMode { byte, pos } => write!(
                f,
                "Ambiguous IUPAC base '{byte}' at position {pos}; use calc_degenerate_tm_bounds."
            ),
            Self::InvalidConcentration(msg) => write!(f, "{msg}"),
            Self::AssemblyDesignFailed {
                junction_index,
                reason,
            } => write!(
                f,
                "Design failed at locus/junction #{junction_index}: {reason}"
            ),
        }
    }
}

impl std::error::Error for OligoError {}

/// Strips whitespace/digits, converts RNA `U -> T`, uppercases, and validates IUPAC DNA characters.
pub fn normalize_dna(raw: &str) -> Result<String, OligoError> {
    let mut out = String::with_capacity(raw.len());
    for (pos, ch) in raw.chars().enumerate() {
        if ch.is_ascii_whitespace() || ch.is_ascii_digit() || ch == '-' || ch == '_' {
            continue;
        }
        let up = match ch.to_ascii_uppercase() {
            'U' => 'T',
            other => other,
        };
        if matches!(
            up,
            'A' | 'C' | 'G' | 'T' | 'R' | 'Y' | 'S' | 'W' | 'K' | 'M' | 'B' | 'D' | 'H' | 'V' | 'N'
        ) {
            out.push(up);
        } else {
            return Err(OligoError::InvalidBase { byte: ch, pos });
        }
    }
    if out.is_empty() {
        return Err(OligoError::EmptySequence);
    }
    Ok(out)
}

/// Returns true if `seq` contains strictly unambiguous uppercase `A, C, G, T` bases.
#[inline]
pub fn is_unambiguous_dna(seq: &str) -> bool {
    !seq.is_empty()
        && seq
            .bytes()
            .all(|b| matches!(b, b'A' | b'C' | b'G' | b'T'))
}

/// Returns the IUPAC Watson–Crick complement of an uppercase ASCII nucleotide byte.
#[inline]
pub fn complement_base(b: u8) -> u8 {
    match b {
        b'A' => b'T',
        b'T' | b'U' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'S' => b'S',
        b'W' => b'W',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        _ => b'N',
    }
}

/// Computes the `5' -> 3'` reverse complement of a normalized IUPAC DNA string.
pub fn reverse_complement(seq: &str) -> String {
    seq.bytes()
        .rev()
        .map(|b| complement_base(b.to_ascii_uppercase()) as char)
        .collect()
}

/// Computes fractional GC content (`0.0..=1.0`), weighting IUPAC degenerate codes stoichiometrically.
pub fn gc_fraction(seq: &str) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let mut gc_sum = 0.0_f64;
    let mut total = 0.0_f64;
    for b in seq.bytes() {
        total += 1.0;
        gc_sum += match b.to_ascii_uppercase() {
            b'G' | b'C' | b'S' => 1.0,
            b'R' | b'Y' | b'K' | b'M' | b'N' => 0.5,
            b'B' | b'V' => 2.0 / 3.0,
            b'D' | b'H' => 1.0 / 3.0,
            _ => 0.0,
        };
    }
    gc_sum / total
}

/// Returns true if an unambiguous sequence is identical to its own reverse complement.
pub fn is_self_complementary(seq: &str) -> bool {
    let b = seq.as_bytes();
    let n = b.len();
    if n == 0 || n % 2 != 0 {
        return false;
    }
    for i in 0..(n / 2) {
        if complement_base(b[i]) != b[n - 1 - i] {
            return false;
        }
    }
    true
}

/// Expands IUPAC degenerate bases into concrete `A/C/G/T` variants up to `max_variants`.
pub fn expand_iupac(raw_seq: &str, max_variants: usize) -> Result<Vec<String>, OligoError> {
    let seq = normalize_dna(raw_seq)?;
    let cap = max_variants.max(1);
    let mut variants = vec![String::with_capacity(seq.len())];

    for (pos, b) in seq.bytes().enumerate() {
        let choices: &[u8] = match b {
            b'A' => b"A",
            b'C' => b"C",
            b'G' => b"G",
            b'T' => b"T",
            b'R' => b"AG",
            b'Y' => b"CT",
            b'S' => b"GC",
            b'W' => b"AT",
            b'K' => b"GT",
            b'M' => b"AC",
            b'B' => b"CGT",
            b'D' => b"AGT",
            b'H' => b"ACT",
            b'V' => b"ACG",
            b'N' => b"ACGT",
            other => {
                return Err(OligoError::InvalidBase {
                    byte: other as char,
                    pos,
                })
            }
        };

        if variants.len().saturating_mul(choices.len()) > cap {
            return Err(OligoError::InvalidConcentration(
                "IUPAC degeneracy exceeds max_variants limit.",
            ));
        }

        let mut next = Vec::with_capacity(variants.len() * choices.len());
        for prefix in &variants {
            for &ch in choices {
                let mut s = prefix.clone();
                s.push(ch as char);
                next.push(s);
            }
        }
        variants = next;
    }

    Ok(variants)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn integrates_tools_1_through_20() {
        let aci = parse_aci_csv(&demo_fvcb_aci_csv()).unwrap();
        let fit = fit_fvcb_aci_curve(&aci).unwrap();
        assert!(fit.r_squared > 0.999);

        let do_pts = parse_stream_do_csv(&demo_stream_diel_do_csv()).unwrap();
        let stream = fit_odum_stream_metabolism(&do_pts).unwrap();
        assert!(stream.gpp_g_o2_m2_d > 2.0);
    }
}
