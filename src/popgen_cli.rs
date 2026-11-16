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

//! Command-line subcommands (`clap` v4) for Tool #16 (`osabt popgen`):
//! - `osabt popgen gui` — Launches the interactive `egui` Haplotype Network & `F_ST` Studio
//! - `osabt popgen diversity` — Computes `pi`, `theta_W`, `H_d`, Tajima's `D`, and `d_XY` / Hudson `F_ST`
//! - `osabt popgen vcf` — Computes diploid Weir & Cockerham (1984) `F_IS/F_ST/F_IT` & Wigginton exact `HWE`
//! - `osabt popgen network` — Builds Bandelt et al. (1999) Median-Joining network & exports `.svg`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::popgen_core::{
    analyze_diploid_vcf, analyze_popgen_fasta, demo_pnw_popgen_fasta, demo_pnw_popgen_vcf,
    parse_popgen_fasta, parse_popgen_vcf,
};
use crate::popgen_network::{build_median_joining_network, render_popgen_network_svg};

/// Top-level CLI arguments for `osabt popgen`.
#[derive(Debug, Args)]
pub struct PopGenArgs {
    #[command(subcommand)]
    pub command: Option<PopGenSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum PopGenSubcommand {
    /// Launch the interactive `egui` Median-Joining Haplotype Network & `F_ST` Studio.
    Gui {
        /// Optional aligned population FASTA (`>id|population`).
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Optional diploid biallelic `.vcf` file.
        #[arg(long)]
        vcf: Option<PathBuf>,
    },

    /// Compute nucleotide diversity (`pi`), `theta_W`, `H_d`, Tajima's (1989) `D`, and pairwise `d_XY` / `F_ST`.
    Diversity {
        /// Optional aligned multi-FASTA file (omit for built-in PNW Dicamptodon `cytb` dataset).
        #[arg(long)]
        fasta: Option<PathBuf>,
    },

    /// Audit a diploid `.vcf` for Weir & Cockerham (1984) `F_IS/F_ST/F_IT` and Wigginton (2005) exact `HWE`.
    Vcf {
        /// Optional diploid `.vcf` file (omit for built-in PNW 18-sample adaptive + null-allele VCF).
        #[arg(long)]
        input: Option<PathBuf>,
    },

    /// Construct a Bandelt et al. (1999) Median-Joining haplotype network and export `.svg`.
    Network {
        /// Optional aligned multi-FASTA file.
        #[arg(long)]
        fasta: Option<PathBuf>,

        /// Bandelt homoplasy tolerance parameter `epsilon` (default `0`).
        #[arg(long, default_value_t = 0)]
        epsilon: usize,

        /// Optional output path for the multi-population pie-chart haplotype network `.svg`.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

/// Dispatches `osabt popgen` CLI subcommands.
pub fn run_popgen_cli(args: PopGenArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_gui(None, None),
        Some(PopGenSubcommand::Gui { fasta, vcf }) => launch_default_gui(fasta, vcf),
        Some(PopGenSubcommand::Diversity { fasta }) => {
            let raw = if let Some(fp) = fasta {
                fs::read_to_string(&fp)
                    .map_err(|e| format!("Failed reading FASTA '{}': {e}", fp.display()))?
            } else {
                demo_pnw_popgen_fasta()
            };
            let recs = parse_popgen_fasta(&raw)?;
            let rep = analyze_popgen_fasta(&recs)?;
            println!(
                "=== OSABT PopGen — Sequence Polymorphism & Tajima's D ({} Seqs x {} bp) ===",
                rep.total_sequences, rep.alignment_len_bp
            );
            for m in std::iter::once(&rep.global_metrics).chain(rep.per_population.iter()) {
                println!(
                    "  {:21} | n={:2} | H={:2} (Hd={:.3}) | S={:2} | pi={:.4} | thetaW={:.4} | Tajima D={:+.3} (p={:.4})",
                    m.population,
                    m.sample_size_n,
                    m.num_haplotypes,
                    m.haplotype_diversity_hd,
                    m.segregating_sites_s,
                    m.nucleotide_diversity_pi,
                    m.watterson_theta_w,
                    m.tajima_d,
                    m.tajima_d_pvalue
                );
            }
            println!("\n=== Pairwise Population Divergence (d_XY, net d_a & Hudson F_ST) ===");
            for p in &rep.pairwise_divergence {
                println!(
                    "  {:21} vs {:21} : d_XY={:.4} | d_a={:.4} | F_ST={:.4}",
                    p.pop_a, p.pop_b, p.d_xy, p.net_d_a, p.hudson_fst
                );
            }
            Ok(())
        }
        Some(PopGenSubcommand::Vcf { input }) => {
            let raw = if let Some(vp) = input {
                fs::read_to_string(&vp)
                    .map_err(|e| format!("Failed reading VCF '{}': {e}", vp.display()))?
            } else {
                demo_pnw_popgen_vcf().to_string()
            };
            let (pops, loci) = parse_popgen_vcf(&raw)?;
            let rep = analyze_diploid_vcf(&pops, &loci)?;
            println!(
                "=== OSABT PopGen — Weir & Cockerham (1984) Diploid VCF ({} Samples, {} SNPs) ===",
                rep.num_samples, rep.num_loci
            );
            println!(
                "  Multilocus F_ST = {:.4} | F_IS = {:+.4} | F_IT = {:+.4}",
                rep.multilocus_fst, rep.multilocus_fis, rep.multilocus_fit
            );
            for (pa, pb, pf) in &rep.pairwise_wc_fst {
                println!("  Pairwise WC F_ST ({pa} vs {pb}) = {pf:.4}");
            }
            println!("\n=== Per-Locus Wigginton (2005) Exact HWE & Weir–Cockerham F_ST ===");
            for l in &rep.loci {
                println!(
                    "  {:15} ({}:{}) | MAF={:.3} | H_obs={:.3} H_exp={:.3} | HWE_p={:.4} | F_ST={:.4}",
                    l.snp_id,
                    l.chrom,
                    l.pos,
                    l.minor_allele_freq,
                    l.obs_heterozygosity,
                    l.exp_heterozygosity,
                    l.hwe_exact_pvalue,
                    l.weir_cockerham_fst
                );
            }
            Ok(())
        }
        Some(PopGenSubcommand::Network {
            fasta,
            epsilon,
            svg_out,
        }) => {
            let raw = if let Some(fp) = fasta {
                fs::read_to_string(&fp)
                    .map_err(|e| format!("Failed reading FASTA '{}': {e}", fp.display()))?
            } else {
                demo_pnw_popgen_fasta()
            };
            let recs = parse_popgen_fasta(&raw)?;
            let net = build_median_joining_network(&recs, epsilon)?;
            println!(
                "=== Bandelt et al. (1999) Median-Joining Network (epsilon = {}) ===",
                net.epsilon
            );
            println!(
                "  Sampled Haplotypes: {} | Inferred Median Vectors: {} | Total Steps: {}",
                net.sampled_haplotypes, net.inferred_median_vectors, net.total_network_steps
            );
            for e in &net.edges {
                println!(
                    "  {} <--> {} : {} step(s) [{}]",
                    net.nodes[e.from_idx].node_id,
                    net.nodes[e.to_idx].node_id,
                    e.mutations,
                    e.subst_labels.join(", ")
                );
            }
            if let Some(out) = svg_out {
                let svg = render_popgen_network_svg(&net, "OSABT PopGen Haplotype Network");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing network SVG '{}': {e}", out.display()))?;
                println!("Wrote Median-Joining Network SVG -> {}", out.display());
            }
            Ok(())
        }
    }
}

#[cfg(feature = "gui")]
fn launch_default_gui(
    fasta: Option<PathBuf>,
    vcf: Option<PathBuf>,
) -> Result<(), String> {
    crate::popgen_gui::run_popgen_gui(fasta, vcf)
}

#[cfg(not(feature = "gui"))]
fn launch_default_gui(
    _fasta: Option<PathBuf>,
    _vcf: Option<PathBuf>,
) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `gui` feature; rebuild with `--features gui` or run `osabt popgen diversity`.",
    ))
}