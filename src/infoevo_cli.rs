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

//! Command-line interface (`clap` v4) for Tool #11 (`osabt infoevo`):
//! - `osabt infoevo tui` — Interactive 4-tab `Ratatui` ITME & Sequence Logo Studio
//! - `osabt infoevo analyze` — Computes Schneider `R_seq`, Adami physical complexity `C`, and exports `.svg` Sequence Logo
//! - `osabt infoevo epistasis` — Computes Dunn, Wahl & Gloor (2008) APC-corrected `MIp` coevolution pairs & `.svg` arc map
//! - `osabt infoevo divergence` — Computes two-clade Jensen–Shannon (`D_JS`) & Kullback–Leibler (`D_KL`) divergence per site
//! - `osabt infoevo simulate` — Runs Kimura (1961) adaptive genetic information vs. substitutional load trajectory

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::infoevo_core::{
    compute_alignment_information, compute_clade_divergence, compute_epistasis_mip_matrix,
    demo_ecoli_promoter_and_stemloop_fasta, parse_fasta_alignment,
    simulate_kimura_selection_information, SequenceAlignment,
};
use crate::infoevo_svg::{
    render_epistasis_mip_svg, render_kimura_trajectory_svg, render_sequence_logo_svg,
};

/// Top-level CLI arguments for `osabt infoevo`.
#[derive(Debug, Args)]
pub struct InfoEvoArgs {
    #[command(subcommand)]
    pub command: Option<InfoEvoSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum InfoEvoSubcommand {
    /// Launch the interactive 4-tab `Ratatui` terminal Information-Theoretic Evolution studio.
    Tui {
        /// Optional aligned multi-FASTA file (omit to load built-in Pribnow + RNA stem-loop demo).
        #[arg(long)]
        fasta: Option<PathBuf>,
    },

    /// Compute Schneider & Stephens (1990) `R_seq` logo profile and Adami et al. (2000) physical complexity `C`.
    Analyze {
        /// Aligned multi-FASTA path (omit to analyze built-in demo alignment).
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Disable Basharin/Schneider small-sample bias correction `e_n`.
        #[arg(long, default_value_t = false)]
        no_small_sample_corr: bool,

        /// Optional host genome length (`bp`) to compute Schneider `R_frequency = log2(G / sites)`.
        #[arg(long)]
        genome_bp: Option<f64>,

        /// Optional target binding site count (`gamma`) paired with `--genome-bp`.
        #[arg(long)]
        target_sites: Option<f64>,

        /// Optional path to write publication-grade `.svg` Sequence Logo.
        #[arg(long)]
        logo_svg: Option<PathBuf>,
    },

    /// Rank pairwise coevolving / epistatic sites using Dunn, Wahl & Gloor (2008) APC-corrected `MIp`.
    Epistasis {
        /// Aligned multi-FASTA path (omit to analyze built-in RNA stem-loop demo).
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Maximum allowed gap fraction per column (`0.0..=1.0`).
        #[arg(long, default_value_t = 0.30)]
        max_gap_frac: f64,

        /// Number of top coevolving pairs to report.
        #[arg(long, default_value_t = 15)]
        top: usize,

        /// Optional path to write `.svg` `MIp` coevolution arc & conservation diagram.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Compute per-site Jensen–Shannon (`D_JS`) and Kullback–Leibler (`D_KL`) divergence between two clades.
    Divergence {
        /// Aligned multi-FASTA path (omit to compare `Thermal` vs `Mesic` ecotypes in built-in demo).
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Substring matching Clade A FASTA IDs (leave empty to auto-partition by `Clade|ID` prefix).
        #[arg(long, default_value = "")]
        clade_a: String,

        /// Substring matching Clade B FASTA IDs (leave empty to auto-partition by `Clade|ID` prefix).
        #[arg(long, default_value = "")]
        clade_b: String,
    },

    /// Simulate Kimura (1961) genetic information accumulation `I_t = D_KL(p_t || p_0)` vs. substitutional load `L_sub,t`.
    Simulate {
        /// Initial frequency `p_0*` of the beneficial allele (`0 < p_0* < 1`).
        #[arg(long, default_value_t = 0.01)]
        p0_fittest: f64,

        /// Malthusian selective advantage `s` of the fittest allele over wild-type.
        #[arg(long, default_value_t = 0.12)]
        selection_s: f64,

        /// Number of generations to simulate.
        #[arg(long, default_value_t = 240)]
        generations: usize,

        /// Optional path to write `.svg` Kimura trajectory plot.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

/// Dispatches `osabt infoevo` CLI subcommands.
pub fn run_infoevo_cli(args: InfoEvoArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_tui(None),
        Some(InfoEvoSubcommand::Tui { fasta }) => launch_default_tui(fasta),
        Some(InfoEvoSubcommand::Analyze {
            fasta,
            no_small_sample_corr,
            genome_bp,
            target_sites,
            logo_svg,
        }) => {
            let aln = load_or_demo_fasta(fasta.as_ref())?;
            let r_freq_pair = genome_bp.zip(target_sites);
            let rep = compute_alignment_information(&aln, !no_small_sample_corr, r_freq_pair);

            println!("=== OSABT InfoEvo — Sequence Information & Adami Complexity ===");
            println!(
                "Alignment             : {} ({}, {} seqs × {} cols)",
                rep.title,
                rep.alphabet.label(),
                rep.num_sequences,
                rep.alignment_length
            );
            println!(
                "Schneider R_sequence  : {:.4} bits (small-sample e_n = {})",
                rep.total_r_sequence_bits, !no_small_sample_corr
            );
            println!(
                "Adami Complexity (C)  : {:.4} bits ({:.1}% density, {:.2} effective invariant sites)",
                rep.adami_physical_complexity_bits,
                rep.complexity_density * 100.0,
                rep.effective_conserved_sites
            );
            if let Some(rf) = rep.r_frequency_target_bits {
                println!("Schneider R_frequency : {rf:.4} bits");
            }
            println!("\nPos\tCons\tN\tH(l)\tR_seq(bits)");
            for s in &rep.sites {
                println!(
                    "{}\t{}\t{}\t{:.3}\t{:.3}",
                    s.position_1based,
                    s.consensus_symbol,
                    s.non_gap_count,
                    s.shannon_bits,
                    s.r_sequence_bits
                );
            }

            if let Some(out) = logo_svg {
                let svg = render_sequence_logo_svg(&rep, &aln.title);
                fs::write(&out, svg).map_err(|e| {
                    format!("Failed writing sequence logo SVG '{}': {e}", out.display())
                })?;
                println!("\nWrote Sequence Logo SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(InfoEvoSubcommand::Epistasis {
            fasta,
            max_gap_frac,
            top,
            svg_out,
        }) => {
            let aln = load_or_demo_fasta(fasta.as_ref())?;
            let epi = compute_epistasis_mip_matrix(&aln, max_gap_frac, top);
            println!(
                "=== Dunn, Wahl & Gloor (2008) APC-Corrected Epistasis (`MIp`) — {} ===",
                aln.title
            );
            println!(
                "Mean Raw MI = {:.4} bits | MIp SD = {:.4} bits",
                epi.mean_raw_mi_bits, epi.mip_std_dev
            );
            for (rank, p) in epi.top_pairs.iter().enumerate() {
                println!(
                    "  #{:02}  Pos {:>3} ↔ {:<3} | MI={:.3}b | APC={:.3}b | MIp={:+.3}b | Z={:+.2} | {}",
                    rank + 1,
                    p.pos_i_1based,
                    p.pos_j_1based,
                    p.mi_bits,
                    p.apc_bits,
                    p.mip_bits,
                    p.z_score,
                    p.top_covarying_states
                );
            }

            if let Some(out) = svg_out {
                let info = compute_alignment_information(&aln, true, None);
                let svg = render_epistasis_mip_svg(&epi, &info, &aln.title);
                fs::write(&out, svg).map_err(|e| {
                    format!("Failed writing epistasis SVG '{}': {e}", out.display())
                })?;
                println!("Wrote Epistasis MIp Arc SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(InfoEvoSubcommand::Divergence {
            fasta,
            clade_a,
            clade_b,
        }) => {
            let aln = load_or_demo_fasta(fasta.as_ref())?;
            let div = compute_clade_divergence(&aln, &clade_a, &clade_b)?;
            println!(
                "=== Two-Clade Informational Divergence: {} (n={}) vs {} (n={}) ===",
                div.clade_a_name, div.clade_a_n, div.clade_b_name, div.clade_b_n
            );
            println!(
                "Total D_JS = {:.4} bits (mean {:.4}/site) | Total D_KL(A||B) = {:.4} bits",
                div.total_js_bits, div.mean_js_bits, div.total_kl_bits
            );
            for s in div.sites.iter().filter(|s| s.js_divergence_bits > 0.05) {
                println!(
                    "  Pos {:>3} : {} ('{}') vs {} ('{}') | D_JS = {:.4} bits | D_KL = {:.4} bits",
                    s.position_1based,
                    div.clade_a_name,
                    s.clade_a_consensus,
                    div.clade_b_name,
                    s.clade_b_consensus,
                    s.js_divergence_bits,
                    s.kl_p_from_q_bits
                );
            }
            Ok(())
        }
        Some(InfoEvoSubcommand::Simulate {
            p0_fittest,
            selection_s,
            generations,
            svg_out,
        }) => {
            let p0 = p0_fittest.clamp(1e-6, 0.90);
            let rem = 1.0 - p0;
            let sim = simulate_kimura_selection_information(
                &[p0, 0.4 * rem, 0.6 * rem],
                &[selection_s, 0.4 * selection_s, 0.0],
                1.0,
                generations,
            )?;
            println!("=== Kimura (1961) Genetic Information Accumulation & Regret ===");
            println!(
                "Initial p0* = {:.5} | Kimura Limit -log2(p0*) = {:.4} bits",
                sim.initial_fittest_freq, sim.kimura_limit_bits
            );
            println!(
                "Final KL Info Gain I_T = {:.4} bits | Substitutional Load L_sub,T = {:.4} bits | Efficiency = {:.2}%",
                sim.final_info_gain_bits,
                sim.final_substitutional_load_bits,
                sim.final_efficiency * 100.0
            );
            if let Some(out) = svg_out {
                let svg = render_kimura_trajectory_svg(&sim, "Kimura Selection Trajectory");
                fs::write(&out, svg).map_err(|e| {
                    format!("Failed writing Kimura SVG '{}': {e}", out.display())
                })?;
                println!("Wrote Kimura Trajectory SVG -> {}", out.display());
            }
            Ok(())
        }
    }
}

fn load_or_demo_fasta(path: Option<&PathBuf>) -> Result<SequenceAlignment, String> {
    if let Some(p) = path {
        let raw = fs::read_to_string(p)
            .map_err(|e| format!("Failed reading FASTA '{}': {e}", p.display()))?;
        let name = p
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("Alignment");
        parse_fasta_alignment(name, &raw)
    } else {
        parse_fasta_alignment(
            "Ecoli_Pribnow_RNA_StemLoop_Demo",
            demo_ecoli_promoter_and_stemloop_fasta(),
        )
    }
}

#[cfg(feature = "tui")]
fn launch_default_tui(fasta: Option<PathBuf>) -> Result<(), String> {
    crate::infoevo_tui::run_infoevo_tui(fasta)
}

#[cfg(not(feature = "tui"))]
fn launch_default_tui(_fasta: Option<PathBuf>) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `tui` feature; rebuild with `--features tui` or run `osabt infoevo analyze`.",
    ))
}
