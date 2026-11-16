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

//! Command-line subcommands (`clap` v4) for Tool #4: `osabt phylo`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::phylo_core::{compare_robinson_foulds, parse_phylo_trees, PhyloTree};
use crate::phylo_render::{render_phylo_svg, render_unicode_tree};

#[derive(Debug, Args)]
pub struct PhyloArgs {
    #[command(subcommand)]
    pub command: Option<PhyloSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum PhyloSubcommand {
    /// Launch the interactive SSH/terminal Ratatui tree explorer (`less` for Newick & NEXUS).
    Tui {
        /// Optional `.nwk`, `.tree`, or `#NEXUS` file (loads built-in 14-taxon Vertebrata ML tree if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
    },

    /// Print a Unicode phylogram or cladogram and summary topology statistics to stdout.
    View {
        /// Input `.nwk` or `#NEXUS` file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Render equal-step topological cladogram instead of proportional branch lengths.
        #[arg(long)]
        cladogram: bool,
        /// Collapse internal clades with bootstrap support strictly below this threshold (%).
        #[arg(long)]
        collapse_below: Option<f64>,
        /// Horizontal character width allocated to tree branches.
        #[arg(short, long, default_value_t = 52)]
        width: usize,
    },

    /// Reroot a Newick/NEXUS tree by outgroup taxon or midpoint diameter and export Newick.
    Reroot {
        /// Input `.nwk` or `#NEXUS` file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Outgroup taxon name (mutually compatible or overridden by `--midpoint`).
        #[arg(long)]
        outgroup: Option<String>,
        /// Reroot at the midpoint of the longest leaf-to-leaf patristic path.
        #[arg(long)]
        midpoint: bool,
        /// Ladderize clades after rerooting.
        #[arg(long, default_value_t = true)]
        ladderize: bool,
        /// Output `.nwk` file path (prints to stdout if omitted).
        #[arg(short, long)]
        output: Option<PathBuf>,
    },

    /// Test whether a target list of taxa is monophyletic, paraphyletic, or polyphyletic.
    Monophyly {
        /// Input `.nwk` or `#NEXUS` file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Target taxon names to test (e.g. `Homo_sapiens Pan_troglodytes Mus_musculus`).
        #[arg(required = true)]
        taxa: Vec<String>,
    },

    /// Compute Robinson & Foulds (1981) symmetric split distance and branch-score L1 between two trees.
    RfDiff {
        /// First tree file (`.nwk` or `#NEXUS`).
        #[arg(long)]
        tree_a: PathBuf,
        /// Second tree file (`.nwk` or `#NEXUS`).
        #[arg(long)]
        tree_b: PathBuf,
    },

    /// Export a publication-ready rectangular or circular vector `.svg` phylogram.
    Svg {
        /// Input `.nwk` or `#NEXUS` file (uses built-in demo if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output `.svg` path.
        #[arg(short, long, default_value = "phylo_tree.svg")]
        output: PathBuf,
        /// Render in radial polar / circular layout instead of rectangular phylogram.
        #[arg(long)]
        circular: bool,
        /// Render topological cladogram branches of equal step length.
        #[arg(long)]
        cladogram: bool,
    },
}

pub fn run_phylo_cli(args: PhyloArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "tui")]
            {
                return crate::phylo_tui::run_phylo_tui(None);
            }
            #[cfg(not(feature = "tui"))]
            {
                return Err(String::from(
                    "Build with `--features tui` or run `osabt phylo view`.",
                ));
            }
        }
        Some(PhyloSubcommand::Tui { input }) => {
            #[cfg(feature = "tui")]
            {
                crate::phylo_tui::run_phylo_tui(input)
            }
            #[cfg(not(feature = "tui"))]
            {
                let _ = input;
                Err(String::from("TUI feature is disabled in this build."))
            }
        }
        Some(PhyloSubcommand::View {
            input,
            cladogram,
            collapse_below,
            width,
        }) => {
            let mut tree = load_first_tree_or_demo(input.as_ref())?;
            if let Some(min_sup) = collapse_below {
                tree.collapse_weak_clades(min_sup);
            }
            let st = tree.stats();
            println!(
                "=== {} ({} tips | {} internal | total BL {:.4} | max root-to-tip {:.4}) ===",
                tree.tree_name,
                st.leaf_count,
                st.internal_count,
                st.total_branch_length,
                st.max_root_to_tip
            );
            for line in render_unicode_tree(&tree, !cladogram, width, "") {
                println!("{}", line.text);
            }
            Ok(())
        }
        Some(PhyloSubcommand::Reroot {
            input,
            outgroup,
            midpoint,
            ladderize,
            output,
        }) => {
            let mut tree = load_first_tree_or_demo(input.as_ref())?;
            if midpoint {
                let diam = tree.midpoint_reroot()?;
                eprintln!(
                    "Midpoint rerooted '{}' (max patristic diameter = {:.5}).",
                    tree.tree_name, diam
                );
            } else if let Some(og) = outgroup {
                tree.reroot_by_outgroup(&og)?;
                eprintln!("Rerooted '{}' on outgroup '{og}'.", tree.tree_name);
            } else {
                return Err(String::from(
                    "Specify either `--midpoint` or `--outgroup <TAXON>`.",
                ));
            }
            if ladderize {
                tree.ladderize(true);
            }
            let nwk = tree.to_newick();
            if let Some(out_path) = output {
                fs::write(&out_path, format!("{nwk}\n"))
                    .map_err(|e| format!("Failed writing '{}': {e}", out_path.display()))?;
                println!("Saved rerooted Newick -> {}", out_path.display());
            } else {
                println!("{nwk}");
            }
            Ok(())
        }
        Some(PhyloSubcommand::Monophyly { input, taxa }) => {
            let tree = load_first_tree_or_demo(input.as_ref())?;
            let rep = tree.check_monophyly(&taxa)?;
            println!("Tree           : {}", tree.tree_name);
            println!("Query Taxa     : {}", rep.matched_taxa.join(", "));
            if !rep.missing_query_taxa.is_empty() {
                println!("Missing Taxa   : {}", rep.missing_query_taxa.join(", "));
            }
            println!("Assessment     : {:?}", rep.status);
            println!("MRCA Node ID   : #{}", rep.mrca_node_id);
            if let Some(sup) = rep.mrca_support_pct {
                println!("MRCA Support   : {sup:.1}%");
            }
            if !rep.intruder_taxa.is_empty() {
                println!("Intruder Taxa  : {}", rep.intruder_taxa.join(", "));
            }
            Ok(())
        }
        Some(PhyloSubcommand::RfDiff { tree_a, tree_b }) => {
            let ta = load_first_tree_or_demo(Some(&tree_a))?;
            let tb = load_first_tree_or_demo(Some(&tree_b))?;
            let rf = compare_robinson_foulds(&ta, &tb)?;
            println!("=== Robinson-Foulds Topological Comparison ===");
            println!("Shared Taxa            : {}", rf.shared_taxa_count);
            println!("Non-trivial Splits (A) : {}", rf.splits_tree_a);
            println!("Non-trivial Splits (B) : {}", rf.splits_tree_b);
            println!("Shared Splits          : {}", rf.shared_splits);
            println!(
                "Robinson-Foulds (RF)   : {} / {} (nRF = {:.4})",
                rf.rf_distance, rf.max_rf, rf.normalized_rf
            );
            println!("Branch Score (L1)      : {:.6}", rf.branch_score_l1);
            Ok(())
        }
        Some(PhyloSubcommand::Svg {
            input,
            output,
            circular,
            cladogram,
        }) => {
            let tree = load_first_tree_or_demo(input.as_ref())?;
            let svg = render_phylo_svg(&tree, !cladogram, circular);
            fs::write(&output, svg)
                .map_err(|e| format!("Failed writing SVG '{}': {e}", output.display()))?;
            println!(
                "Exported {} phylogenetic SVG -> {}",
                if circular { "circular" } else { "rectangular" },
                output.display()
            );
            Ok(())
        }
    }
}

fn load_first_tree_or_demo(input: Option<&PathBuf>) -> Result<PhyloTree, String> {
    match input {
        Some(p) => {
            let raw = fs::read_to_string(p)
                .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
            let stem = p
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| String::from("Tree"));
            parse_phylo_trees(&stem, &raw)?
                .into_iter()
                .next()
                .ok_or_else(|| String::from("No trees parsed."))
        }
        None => Ok(PhyloTree::demo_vertebrate_amniotes()),
    }
}
