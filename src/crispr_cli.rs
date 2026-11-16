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

//! Command-line subcommands (`clap` v4) for Tool #17 (`osabt crispr`):
//! - `osabt crispr tui` — Interactive 4-tab `Ratatui` CRISPR Guide & `CFD` Studio
//! - `osabt crispr design` — Multi-PAM guide scanner, Doench on-target & `CFD` ranker, `.tsv`/`.svg` exporter
//! - `osabt crispr cfd` — Pairwise Doench et al. (2016) `CFD` off-target probability calculator
//! - `osabt crispr base-edit` — Audits `CBE` (`C->T`) and `ABE` (`A->G`) precision zero-bystander windows.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

use crate::crispr_core::{
    demo_crispr_target_and_offtargets, scan_crispr_candidates, CrisprCloningVector, CrisprNuclease,
};
use crate::crispr_offtarget::{
    calc_cfd_offtarget_score, export_crispr_design_tsv, render_crispr_locus_svg,
    scan_offtargets_and_specificity,
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum NucleaseArg {
    Spcas9,
    Sacas9,
    Cas12a,
}

impl From<NucleaseArg> for CrisprNuclease {
    fn from(v: NucleaseArg) -> Self {
        match v {
            NucleaseArg::Spcas9 => Self::SpCas9,
            NucleaseArg::Sacas9 => Self::SaCas9,
            NucleaseArg::Cas12a => Self::LbCas12a,
        }
    }
}

/// Top-level CLI arguments for `osabt crispr`.
#[derive(Debug, Args)]
pub struct CrisprArgs {
    #[command(subcommand)]
    pub command: Option<CrisprSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum CrisprSubcommand {
    /// Launch the interactive 4-tab `Ratatui` CRISPR Guide Designer & `CFD` Scorer.
    Tui {
        /// Optional target FASTA file to scan.
        #[arg(long)]
        target: Option<PathBuf>,

        /// Optional genomic background FASTA for off-target `CFD` specificity scanning.
        #[arg(long)]
        background: Option<PathBuf>,
    },

    /// Scan a target DNA locus for CRISPR guides, score on-target + `CFD` specificity, and export oligos.
    Design {
        #[arg(long)]
        target: Option<PathBuf>,

        #[arg(long)]
        background: Option<PathBuf>,

        /// Endonuclease PAM family (`spcas9`, `sacas9`, or `cas12a`).
        #[arg(long, value_enum, default_value_t = NucleaseArg::Spcas9)]
        nuclease: NucleaseArg,

        /// Optional path to write the Golden Gate oligo + `CFD` table as `.tsv`.
        #[arg(long)]
        tsv_out: Option<PathBuf>,

        /// Optional path to write the ranked CRISPR target locus map as `.svg`.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Compute exact Doench et al. (2016) `CFD` off-target probability for a 20-nt spacer vs. 23-nt target+PAM.
    Cfd {
        /// 20-nt `5'->3'` `SpCas9` guide spacer (`DNA` alphabet).
        #[arg(long)]
        spacer: String,

        /// 23-nt `5'->3'` candidate off-target protospacer + 3-nt `PAM`.
        #[arg(long)]
        offtarget: String,
    },

    /// Screen guides for precision zero-bystander `CBE` (`C->T`, Komor 2016) & `ABE` (`A->G`, Gaudelli 2017) windows.
    BaseEdit {
        #[arg(long)]
        target: Option<PathBuf>,
    },
}

/// Dispatches `osabt crispr` CLI subcommands.
pub fn run_crispr_cli(args: CrisprArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_tui(None, None),
        Some(CrisprSubcommand::Tui { target, background }) => launch_default_tui(target, background),
        Some(CrisprSubcommand::Design {
            target,
            background,
            nuclease,
            tsv_out,
            svg_out,
        }) => {
            let (demo_t, demo_bg) = demo_crispr_target_and_offtargets();
            let t_raw = if let Some(tp) = target {
                fs::read_to_string(&tp)
                    .map_err(|e| format!("Failed reading target '{}': {e}", tp.display()))?
            } else {
                demo_t
            };
            let bg_raw = if let Some(bp) = background {
                fs::read_to_string(&bp)
                    .map_err(|e| format!("Failed reading background '{}': {e}", bp.display()))?
            } else {
                demo_bg
            };

            let nuc = CrisprNuclease::from(nuclease);
            let (seq, mut guides) = scan_crispr_candidates(&t_raw, nuc)?;
            let _ = scan_offtargets_and_specificity(&mut guides, &bg_raw, 4)?;
            guides.sort_by(|a, b| {
                b.composite_score
                    .partial_cmp(&a.composite_score)
                    .unwrap_or(std::cmp::Ordering::Equal)
            });

            println!(
                "=== OSABT CRISPR — {} Guide Ranking ({} Candidates on {} bp Locus) ===",
                nuc.label(),
                guides.len(),
                seq.len()
            );
            for g in guides.iter().take(12) {
                println!(
                    "  {:8} ({}) cut={:3} | {} {} | On={:4.1} Spec={:5.1} (Offs={}) | Comp={:4.1}",
                    g.guide_id,
                    g.strand,
                    g.cut_site_1idx,
                    g.spacer_seq,
                    g.pam_seq,
                    g.on_target_score,
                    g.cfd_specificity_score,
                    g.off_target_count,
                    g.composite_score
                );
            }

            let vec_choice = match nuc {
                CrisprNuclease::LbCas12a => CrisprCloningVector::LbCas12aBsaI,
                _ => CrisprCloningVector::Px459BbsI,
            };
            if let Some(out) = tsv_out {
                let tsv = export_crispr_design_tsv(&guides, vec_choice);
                fs::write(&out, tsv)
                    .map_err(|e| format!("Failed writing TSV '{}': {e}", out.display()))?;
                println!("Wrote CRISPR Oligo & CFD Table -> {}", out.display());
            }
            if let Some(out) = svg_out {
                let svg = render_crispr_locus_svg(seq.len(), &guides, "OSABT CRISPR Guide Track");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", out.display()))?;
                println!("Wrote CRISPR Locus SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(CrisprSubcommand::Cfd { spacer, offtarget }) => {
            let score = calc_cfd_offtarget_score(&spacer, &offtarget)?;
            println!(
                "Doench et al. (2016) CFD Off-Target Score: {:.4} (Spacer={} vs Site+PAM={})",
                score,
                spacer.to_ascii_uppercase(),
                offtarget.to_ascii_uppercase()
            );
            Ok(())
        }
        Some(CrisprSubcommand::BaseEdit { target }) => {
            let (demo_t, _) = demo_crispr_target_and_offtargets();
            let t_raw = if let Some(tp) = target {
                fs::read_to_string(&tp)
                    .map_err(|e| format!("Failed reading target '{}': {e}", tp.display()))?
            } else {
                demo_t
            };
            let (_, guides) = scan_crispr_candidates(&t_raw, CrisprNuclease::SpCas9)?;
            println!("=== OSABT CRISPR — Precision Zero-Bystander Base-Editor Guides (Window 4..8) ===");
            for g in &guides {
                if g.base_edit.cbe_clean_single_edit || g.base_edit.abe_clean_single_edit {
                    println!(
                        "  {:8} ({}) {} {} | Win[4..8]={} | CBE_clean={} ({:?}) | ABE_clean={} ({:?})",
                        g.guide_id,
                        g.strand,
                        g.spacer_seq,
                        g.pam_seq,
                        g.base_edit.window_seq_4_to_8,
                        g.base_edit.cbe_clean_single_edit,
                        g.base_edit.cbe_c_positions,
                        g.base_edit.abe_clean_single_edit,
                        g.base_edit.abe_a_positions
                    );
                }
            }
            Ok(())
        }
    }
}

#[cfg(feature = "tui")]
fn launch_default_tui(
    target: Option<PathBuf>,
    background: Option<PathBuf>,
) -> Result<(), String> {
    crate::crispr_tui::run_crispr_tui(target, background)
}

#[cfg(not(feature = "tui"))]
fn launch_default_tui(
    _target: Option<PathBuf>,
    _background: Option<PathBuf>,
) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `tui` feature; rebuild with `--features tui` or run `osabt crispr design`.",
    ))
}
