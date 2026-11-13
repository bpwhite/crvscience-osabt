// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #12 (`osabt qpcr`):
//! - `osabt qpcr tui` — Interactive 4-tab `Ratatui` MIQE qPCR & Melt-Curve Studio
//! - `osabt qpcr quantify` — Computes `qBase` (Hellemans 2007), `Pfaffl` (2001), & `Livak` (2001) `2^-ddCq` fold-changes + `.svg`
//! - `osabt qpcr efficiency` — Fits log10 serial-dilution standard curves & Vandesompele (2002) `geNorm` stability `M`
//! - `osabt qpcr melt` — Analyzes SYBR Green `-dF/dT` melt peaks & primer-dimers (Ririe, Rasmussen & Wittwer 1997)
//! - `osabt qpcr miqe` — Runs automated Bustin et al. (2009) MIQE compliance checklist

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::qpcr_core::{
    analyze_qpcr_experiment, demo_qpcr_heatshock_plate_csv, parse_qpcr_table,
    render_qpcr_expression_svg, QpcrReactionWell,
};
use crate::qpcr_melt::{
    analyze_melt_curve, demo_sybr_melt_traces, parse_melt_csv, render_melt_curves_svg, MeltTrace,
};

/// Top-level CLI arguments for `osabt qpcr`.
#[derive(Debug, Args)]
pub struct QpcrArgs {
    #[command(subcommand)]
    pub command: Option<QpcrSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum QpcrSubcommand {
    /// Launch the interactive 4-tab `Ratatui` MIQE qPCR & Melt-Curve terminal studio.
    Tui {
        /// Optional CSV/TSV qPCR reaction table (omit to load built-in heat-shock plate).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Optional CSV/TSV melt-curve trace table.
        #[arg(long)]
        melt: Option<PathBuf>,

        /// Comma-separated reference genes (default: `ACTB,GAPDH`).
        #[arg(long, default_value = "ACTB,GAPDH")]
        refs: String,

        /// Calibrator / control condition group name.
        #[arg(long, default_value = "Control")]
        calibrator: String,
    },

    /// Compute `qBase` multi-reference, `Pfaffl (2001)`, and `Livak (2001) 2^-ddCq` relative expression & Welch `t`-tests.
    Quantify {
        /// Optional CSV/TSV qPCR reaction table (omit to run on built-in heat-shock demo).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Comma-separated internal reference genes.
        #[arg(long, default_value = "ACTB,GAPDH")]
        refs: String,

        /// Calibrator condition group name.
        #[arg(long, default_value = "Control")]
        calibrator: String,

        /// Optional path to export `.svg` fold-change comparison figure.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Fit log10 serial-dilution standard curves (`E = 10^(-1/slope)`) and rank `geNorm` reference stability `M`.
    Efficiency {
        /// Optional CSV/TSV qPCR reaction table.
        #[arg(long)]
        input: Option<PathBuf>,

        /// Comma-separated candidate housekeeping genes for `geNorm` ranking.
        #[arg(long, default_value = "ACTB,GAPDH,B2M")]
        candidates: String,
    },

    /// Call `-dF/dT` SYBR Green melting peaks and audit primer-dimer artifacts (Ririe, Rasmussen & Wittwer 1997).
    Melt {
        /// Optional CSV/TSV melt table (`well,sample,target,temperature_c,fluorescence`).
        #[arg(long)]
        melt: Option<PathBuf>,

        /// Temperature threshold (`°C`) below which peaks are flagged as probable primer-dimers.
        #[arg(long, default_value_t = 77.5)]
        dimer_cutoff_c: f64,

        /// Optional path to write `.svg` `-dF/dT` dissociation curve figure.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Run automated Bustin et al. (2009) MIQE quality-control checklist on a qPCR dataset.
    Miqe {
        /// Optional CSV/TSV qPCR reaction table.
        #[arg(long)]
        input: Option<PathBuf>,

        /// Comma-separated reference genes used for normalization.
        #[arg(long, default_value = "ACTB,GAPDH")]
        refs: String,

        /// Calibrator condition group name.
        #[arg(long, default_value = "Control")]
        calibrator: String,
    },
}

/// Dispatches `osabt qpcr` CLI subcommands.
pub fn run_qpcr_cli(args: QpcrArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_tui(
            None,
            None,
            vec![String::from("ACTB"), String::from("GAPDH")],
            String::from("Control"),
        ),
        Some(QpcrSubcommand::Tui {
            input,
            melt,
            refs,
            calibrator,
        }) => launch_default_tui(input, melt, split_csv_genes(&refs), calibrator),
        Some(QpcrSubcommand::Quantify {
            input,
            refs,
            calibrator,
            svg_out,
        }) => {
            let wells = load_or_demo_wells(input.as_ref())?;
            let ref_list = split_csv_genes(&refs);
            let rep = analyze_qpcr_experiment(&wells, &ref_list, &calibrator)?;

            println!("=== OSABT qPCR — Relative Expression Summary ===");
            println!(
                "Calibrator Group : {} | Reference Panel : {}",
                rep.calibrator_group,
                rep.reference_genes_used.join(" + ")
            );
            for c in &rep.group_comparisons {
                println!(
                    "  {} ({} vs {}): qBase FC = {:.2}x [95% CI {:.2}–{:.2}] | Pfaffl = {:.2}x | 2^-ddCq = {:.2}x | t = {:+.2}, p = {:.3e}",
                    c.target_gene,
                    c.treatment_group,
                    c.calibrator_group,
                    c.fold_change_qbase,
                    c.fold_change_ci95.0,
                    c.fold_change_ci95.1,
                    c.pfaffl_ratio,
                    c.livak_ratio,
                    c.welch_t_stat,
                    c.p_value_two_sided
                );
            }

            if let Some(out) = svg_out {
                let svg = render_qpcr_expression_svg(&rep, "OSABT RT-qPCR Quantification");
                fs::write(&out, svg).map_err(|e| {
                    format!("Failed writing expression SVG '{}': {e}", out.display())
                })?;
                println!("Wrote Expression SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(QpcrSubcommand::Efficiency { input, candidates }) => {
            let wells = load_or_demo_wells(input.as_ref())?;
            let cand_list = split_csv_genes(&candidates);
            let rep = analyze_qpcr_experiment(&wells, &cand_list[..cand_list.len().min(2)], "Control")?;

            println!("=== Serial-Dilution Standard Curves ===");
            for sc in &rep.standard_curves {
                println!(
                    "  {:8} : Slope = {:.3} | Intercept = {:.2} | R² = {:.4} | E = {:.3} ({:.1}%) | MIQE {}",
                    sc.target_gene,
                    sc.slope,
                    sc.intercept_cq,
                    sc.r_squared,
                    sc.amplification_factor_e,
                    sc.efficiency_pct,
                    if sc.miqe_pass { "PASS" } else { "FLAG" }
                );
            }
            println!("\n=== Vandesompele et al. (2002) geNorm Stability Ranking ===");
            for g in &rep.genorm_rankings {
                println!(
                    "  {:8} : geNorm M = {:.4} [{}]",
                    g.gene_name,
                    g.stability_m,
                    if g.stable_pass { "STABLE" } else { "UNSTABLE" }
                );
            }
            Ok(())
        }
        Some(QpcrSubcommand::Melt {
            melt,
            dimer_cutoff_c,
            svg_out,
        }) => {
            let traces = load_or_demo_melt(melt.as_ref())?;
            let results: Vec<_> = traces
                .iter()
                .map(|t| analyze_melt_curve(t, dimer_cutoff_c))
                .collect();
            println!("=== Ririe, Rasmussen & Wittwer (1997) SYBR Melt-Curve Analysis ===");
            for r in &results {
                let tm = r
                    .primary_tm_c
                    .map(|v| format!("{v:.2} °C"))
                    .unwrap_or_else(|| String::from("None"));
                println!(
                    "  Well {} ({} / {}) : Primary Tm = {} | Peaks = {} | {}",
                    r.well_id,
                    r.sample_id,
                    r.target_gene,
                    tm,
                    r.peaks.len(),
                    r.specificity_status
                );
            }
            if let Some(out) = svg_out {
                let svg = render_melt_curves_svg(&results, "SYBR Dissociation Curves");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing melt SVG '{}': {e}", out.display()))?;
                println!("Wrote Melt-Curve SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(QpcrSubcommand::Miqe {
            input,
            refs,
            calibrator,
        }) => {
            let wells = load_or_demo_wells(input.as_ref())?;
            let rep = analyze_qpcr_experiment(&wells, &split_csv_genes(&refs), &calibrator)?;
            println!("=== Bustin et al. (2009) Automated MIQE Compliance Checklist ===");
            for chk in &rep.miqe_checks {
                println!(
                    "  [{:4}] {:20} | {:16} | {}",
                    chk.severity.label(),
                    chk.category,
                    chk.subject,
                    chk.detail
                );
            }
            Ok(())
        }
    }
}

fn split_csv_genes(raw: &str) -> Vec<String> {
    raw.split(',')
        .map(|s| s.trim().to_string())
        .filter(|s| !s.is_empty())
        .collect()
}

fn load_or_demo_wells(path: Option<&PathBuf>) -> Result<Vec<QpcrReactionWell>, String> {
    if let Some(p) = path {
        let raw = fs::read_to_string(p)
            .map_err(|e| format!("Failed reading qPCR table '{}': {e}", p.display()))?;
        parse_qpcr_table(&raw)
    } else {
        parse_qpcr_table(demo_qpcr_heatshock_plate_csv())
    }
}

fn load_or_demo_melt(path: Option<&PathBuf>) -> Result<Vec<MeltTrace>, String> {
    if let Some(p) = path {
        let raw = fs::read_to_string(p)
            .map_err(|e| format!("Failed reading melt CSV '{}': {e}", p.display()))?;
        parse_melt_csv(&raw)
    } else {
        Ok(demo_sybr_melt_traces())
    }
}

#[cfg(feature = "tui")]
fn launch_default_tui(
    input: Option<PathBuf>,
    melt: Option<PathBuf>,
    refs: Vec<String>,
    calibrator: String,
) -> Result<(), String> {
    crate::qpcr_tui::run_qpcr_tui(input, melt, refs, calibrator)
}

#[cfg(not(feature = "tui"))]
fn launch_default_tui(
    _input: Option<PathBuf>,
    _melt: Option<PathBuf>,
    _refs: Vec<String>,
    _calibrator: String,
) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `tui` feature; rebuild with `--features tui` or run `osabt qpcr quantify`.",
    ))
}
