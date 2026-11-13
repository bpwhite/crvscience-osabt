// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #14 (`osabt edna`):
//! - `osabt edna tui` — Interactive 4-tab `Ratatui` eDNA Demultiplexer & Rarefaction Studio
//! - `osabt edna demux` — Runs dual-barcode FASTQ demultiplexing, `maxEE` filtering, ASV denoising & `.tsv` export
//! - `osabt edna diversity` — Computes `Chao1` (Chao 1984), Hill numbers (`^0D, ^1D, ^2D`), & Bray–Curtis dissimilarity
//! - `osabt edna rarefy` — Computes Hurlbert (1971) exact hypergeometric rarefaction curves & `.svg` figure

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::edna_diversity::{analyze_edna_community_diversity, render_edna_rarefaction_svg};
use crate::edna_fastq::{
    demo_pnw_12s_edna_dataset, export_asv_table_tsv, parse_barcode_sheet_csv, parse_fastq_str,
    run_edna_metabarcoding_pipeline, EdnaPipelineConfig, EdnaPipelineReport,
};

/// Top-level CLI arguments for `osabt edna`.
#[derive(Debug, Args)]
pub struct EdnaArgs {
    #[command(subcommand)]
    pub command: Option<EdnaSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum EdnaSubcommand {
    /// Launch the interactive 4-tab `Ratatui` eDNA Demultiplexer, ASV & Rarefaction studio.
    Tui {
        /// Optional raw multiplexed Phred+33 `.fastq` file (omit for built-in PNW 12S eDNA dataset).
        #[arg(long)]
        fastq: Option<PathBuf>,

        /// Optional dual-barcode CSV/TSV sample sheet (`sample_id,fwd_barcode,rev_barcode,role`).
        #[arg(long)]
        barcodes: Option<PathBuf>,
    },

    /// Demultiplex dual-indexed FASTQ reads, filter by Edgar & Flyvbjerg (2015) `maxEE`, denoise ASVs, and export `.tsv`.
    Demux {
        #[arg(long)]
        fastq: Option<PathBuf>,

        #[arg(long)]
        barcodes: Option<PathBuf>,

        /// Maximum expected errors (`EE`) allowed per trimmed amplicon read.
        #[arg(long, default_value_t = 1.0)]
        max_ee: f64,

        /// Maximum barcode mismatches allowed per index end.
        #[arg(long, default_value_t = 1)]
        max_mismatches: usize,

        /// Exclude chimeras and blank-flagged contaminants from exported TSV table.
        #[arg(long, default_value_t = false)]
        clean_only: bool,

        /// Optional output path for the ASV × Sample `.tsv` community matrix.
        #[arg(long)]
        tsv_out: Option<PathBuf>,
    },

    /// Compute `Chao1` (Chao 1984), Hill numbers (`^0D, ^1D, ^2D`), Hurlbert `PIE`, coverage, and Bray–Curtis dissimilarity.
    Diversity {
        #[arg(long)]
        fastq: Option<PathBuf>,

        #[arg(long)]
        barcodes: Option<PathBuf>,
    },

    /// Generate exact Hurlbert (1971) hypergeometric rarefaction curves across subsample depths and export `.svg`.
    Rarefy {
        #[arg(long)]
        fastq: Option<PathBuf>,

        #[arg(long)]
        barcodes: Option<PathBuf>,

        /// Number of subsample depth steps per environmental site.
        #[arg(long, default_value_t = 12)]
        steps: usize,

        /// Optional path to export 2-panel Rarefaction + Hill-Profile `.svg`.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

/// Dispatches `osabt edna` CLI subcommands.
pub fn run_edna_cli(args: EdnaArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_tui(None, None),
        Some(EdnaSubcommand::Tui { fastq, barcodes }) => launch_default_tui(fastq, barcodes),
        Some(EdnaSubcommand::Demux {
            fastq,
            barcodes,
            max_ee,
            max_mismatches,
            clean_only,
            tsv_out,
        }) => {
            let cfg = EdnaPipelineConfig {
                max_expected_errors: max_ee,
                max_barcode_mismatches: max_mismatches,
                ..EdnaPipelineConfig::default()
            };
            let rep = run_or_demo_pipeline(fastq.as_ref(), barcodes.as_ref(), &cfg)?;
            println!("=== OSABT eDNA — Metabarcoding Demultiplexing & ASV Summary ===");
            println!(
                "Raw Reads: {} | Demuxed: {} | Passed MaxEE ({:.2}): {} | 1-nt Denoised: {} | Chimeric Reads: {}",
                rep.total_raw_reads,
                rep.demux_assigned_reads,
                max_ee,
                rep.passed_ee_reads,
                rep.error_corrected_reads,
                rep.chimeric_reads_flagged
            );
            let tsv = export_asv_table_tsv(&rep, clean_only);
            if let Some(out) = tsv_out {
                fs::write(&out, &tsv)
                    .map_err(|e| format!("Failed writing ASV TSV '{}': {e}", out.display()))?;
                println!("Wrote ASV Community Table -> {}", out.display());
            } else {
                println!("\n{tsv}");
            }
            Ok(())
        }
        Some(EdnaSubcommand::Diversity { fastq, barcodes }) => {
            let rep = run_or_demo_pipeline(
                fastq.as_ref(),
                barcodes.as_ref(),
                &EdnaPipelineConfig::default(),
            )?;
            let div = analyze_edna_community_diversity(&rep, 10);
            println!(
                "=== OSABT eDNA — Alpha & Hill Biodiversity ({} Clean ASVs) ===",
                div.clean_asv_count
            );
            for a in &div.alpha {
                println!(
                    "  {:15} : N={:3} | ^0D={} | Chao1={:.2} [{:.1}–{:.1}] | ^1D={:.2} | ^2D={:.2} | PIE={:.3} | Cov={:.1}%",
                    a.sample_name,
                    a.total_reads_n,
                    a.observed_richness_q0,
                    a.chao1_asymptotic,
                    a.chao1_ci95.0,
                    a.chao1_ci95.1,
                    a.hill_shannon_q1,
                    a.hill_simpson_q2,
                    a.hurlbert_pie,
                    a.sample_coverage * 100.0
                );
            }
            println!("\n=== Pairwise Bray–Curtis Dissimilarity ===");
            for (i, name) in div.sample_names.iter().enumerate() {
                let row_str = div.bray_curtis_matrix[i]
                    .iter()
                    .map(|v| format!("{v:.3}"))
                    .collect::<Vec<_>>()
                    .join("\t");
                println!("  {name:15}\t{row_str}");
            }
            Ok(())
        }
        Some(EdnaSubcommand::Rarefy {
            fastq,
            barcodes,
            steps,
            svg_out,
        }) => {
            let rep = run_or_demo_pipeline(
                fastq.as_ref(),
                barcodes.as_ref(),
                &EdnaPipelineConfig::default(),
            )?;
            let div = analyze_edna_community_diversity(&rep, steps);
            println!("=== Hurlbert (1971) Analytical Hypergeometric Rarefaction ===");
            for c in &div.rarefaction_curves {
                let pairs = c
                    .depths_m
                    .iter()
                    .zip(&c.expected_species_sm)
                    .map(|(m, s)| format!("m={m}:{s:.2}"))
                    .collect::<Vec<_>>()
                    .join("  ");
                println!("  {:15} : {pairs}", c.sample_name);
            }
            if let Some(out) = svg_out {
                let svg = render_edna_rarefaction_svg(&div, "eDNA Rarefaction & Hill Diversity");
                fs::write(&out, svg).map_err(|e| {
                    format!("Failed writing rarefaction SVG '{}': {e}", out.display())
                })?;
                println!("Wrote Rarefaction & Hill SVG -> {}", out.display());
            }
            Ok(())
        }
    }
}

fn run_or_demo_pipeline(
    fastq_path: Option<&PathBuf>,
    bc_path: Option<&PathBuf>,
    cfg: &EdnaPipelineConfig,
) -> Result<EdnaPipelineReport, String> {
    let (demo_fq, demo_bc, demo_db) = demo_pnw_12s_edna_dataset();
    let reads = if let Some(fp) = fastq_path {
        let raw = fs::read_to_string(fp)
            .map_err(|e| format!("Failed reading FASTQ '{}': {e}", fp.display()))?;
        parse_fastq_str(&raw)?
    } else {
        parse_fastq_str(&demo_fq)?
    };
    let bcs = if let Some(bp) = bc_path {
        let raw = fs::read_to_string(bp)
            .map_err(|e| format!("Failed reading barcode sheet '{}': {e}", bp.display()))?;
        parse_barcode_sheet_csv(&raw)?
    } else {
        parse_barcode_sheet_csv(&demo_bc)?
    };
    Ok(run_edna_metabarcoding_pipeline(&reads, &bcs, &demo_db, cfg))
}

#[cfg(feature = "tui")]
fn launch_default_tui(fastq: Option<PathBuf>, barcodes: Option<PathBuf>) -> Result<(), String> {
    crate::edna_tui::run_edna_tui(fastq, barcodes)
}

#[cfg(not(feature = "tui"))]
fn launch_default_tui(_fastq: Option<PathBuf>, _barcodes: Option<PathBuf>) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `tui` feature; rebuild with `--features tui` or run `osabt edna diversity`.",
    ))
}
