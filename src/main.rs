// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Open Source Ad-hoc Biology Tools (OSABT v2.0.0) — Unified CLI, TUI & GUI Entrypoint

use std::process;

use clap::{Parser, Subcommand};

use osabt::bench_cli::{run_bench_cli, BenchArgs};
use osabt::crispr_cli::{run_crispr_cli, CrisprArgs};
use osabt::dwc_cli::{run_dwc_cli, DwcArgs};
use osabt::ecosound_cli::{run_ecosound_cli, EcosoundArgs};
use osabt::edna_cli::{run_edna_cli, EdnaArgs};
use osabt::epi_cli::{run_epi_cli, EpiArgs};
use osabt::flow_cli::{run_flow_cli, FlowArgs};
use osabt::flux_cli::{run_flux_cli, FluxArgs};
use osabt::forge3d_cli::{run_forge3d_cli, Forge3dArgs};
use osabt::infoevo_cli::{run_infoevo_cli, InfoEvoArgs};
use osabt::kinetics_cli::{run_kinetics_cli, KineticsArgs};
use osabt::lexi_cli::{run_lexi_cli, LexiArgs};
use osabt::morpho_cli::{run_morpho_cli, MorphoArgs};
use osabt::phylo_cli::{run_phylo_cli, PhyloArgs};
use osabt::plasmid_cli::{run_plasmid_cli, PlasmidArgs};
use osabt::plate_cli::{run_plate_cli, PlateArgs};
use osabt::popgen_cli::{run_popgen_cli, PopGenArgs};
use osabt::qpcr_cli::{run_qpcr_cli, QpcrArgs};
use osabt::survival_cli::{run_survival_cli, SurvivalArgs};
use osabt::wildlife_cli::{run_wildlife_cli, WildlifeArgs};

#[derive(Debug, Parser)]
#[command(
    name = "osabt",
    version,
    author = "Bryan P. White <https://crvscience.com>",
    about = "Open Source Ad-hoc Biology Tools (OSABT v2.0.0) — 20-tool GPL-3.0 suite for the CRV Science Network"
)]
struct OsabtCli {
    #[command(subcommand)]
    tool: Option<OsabtTool>,
}

#[derive(Debug, Subcommand)]
enum OsabtTool {
    /// Tool #1: Modern primer/oligo designer, Gibson/multiplex planner, and wet-lab bench workbench (CLI + Ratatui TUI).
    Bench(BenchArgs),

    /// Tool #2: Vector map, Gibson/Golden Gate cloning, Sanger .ab1 aligner, and virtual gel studio (CLI + egui GUI).
    Plasmid(PlasmidArgs),

    /// Tool #3: Biodiversity & field-survey Darwin Core linter, coordinate scrubber, synonym resolver, and DwC-A packager (CLI + Ratatui TUI).
    Dwc(DwcArgs),

    /// Tool #4: Terminal phylogenetic tree explorer (`less` for Newick & NEXUS) with rerooting, monophyly, RF diff, and SVG export (CLI + Ratatui TUI).
    Phylo(PhyloArgs),

    /// Tool #5: Epidemiology & study-design power calculator (`R_t` renewal, SIR/SEIR, diagnostic Bayes, Cohort & Test-Negative power) (CLI + Ratatui TUI).
    Epi(EpiArgs),

    /// Tool #6: Taxonomic & anatomical Greek/Latin root morpheme decoder, cognate tree viewer, and Anki flashcard generator (CLI + Ratatui TUI).
    Lexi(LexiArgs),

    /// Tool #7: Biomolecule (PDB), Gielis diatom, radiolarian, & CRV Sim voxel-to-manifold `.stl` / `.3mf` 3D print studio (CLI + egui GUI).
    Forge3d(Forge3dArgs),

    /// Tool #8: Automated Petri CFU colony counter, Excess-Green leaf herbivory analyzer, & 2D Generalized Procrustes morphometrics studio (CLI + egui GUI).
    Morpho(MorphoArgs),

    /// Tool #9: 96/384-well edge-buffered randomized assay planner, Zhang `Z'`-factor QC, & 4PL `IC50`/`EC50` curve fitter (CLI + Ratatui TUI).
    Plate(PlateArgs),

    /// Tool #10: Field bioacoustics STFT spectrogram viewer, `ACI`/`NDSI`/`H`/`ADI` index engine, Raven annotator, & `.wav` call slicer (CLI + egui GUI).
    Ecosound(EcosoundArgs),

    /// Tool #11: Information-theoretic evolution (`ITME`), Schneider sequence logos, Adami physical complexity, `MIp` epistasis, & Kimura regret studio (CLI + Ratatui TUI).
    Infoevo(InfoEvoArgs),

    /// Tool #12: MIQE RT-qPCR relative quantification (`qBase` / `Pfaffl` / `2^-ddCq`), `geNorm` reference stability, & SYBR `-dF/dT` melt-curve studio (CLI + Ratatui TUI).
    Qpcr(QpcrArgs),

    /// Tool #13: Steady-state enzyme kinetics (`MM` / `Hill`), global `AICc` reversible `Ki` mode selector, & Morrison (1969) tight-binding studio (CLI + egui GUI).
    Kinetics(KineticsArgs),

    /// Tool #14: eDNA & metabarcoding dual-barcode FASTQ demultiplexer, `maxEE` ASV denoiser, & Chao1 / Hill rarefaction studio (CLI + Ratatui TUI).
    Edna(EdnaArgs),

    /// Tool #15: Imperfect-detection site occupancy (`MacKenzie 2002`), closed mark-recapture (`M0/Mt/Mh`), & `CJS` survival studio (CLI + Ratatui TUI).
    Wildlife(WildlifeArgs),

    /// Tool #16: Population genetics (`Tajima's D`, `d_XY`, diploid `VCF` Weir–Cockerham `F_ST`, exact `HWE`) & Median-Joining haplotype studio (CLI + egui GUI).
    Popgen(PopGenArgs),

    /// Tool #17: CRISPR-Cas9 / Cas12a / Base-Editor guide designer, Doench `CFD` off-target scorer, & Golden Gate oligo studio (CLI + Ratatui TUI).
    Crispr(CrisprArgs),

    /// Tool #18: Flow cytometry `.fcs 3.1` binary parser, `$SPILLOVER` matrix compensation, `Logicle` scaling, & 2D gating studio (CLI + egui GUI).
    Flow(FlowArgs),

    /// Tool #19: Kaplan–Meier survival, Mantel log-rank, Cox proportional hazards (`HR`), & Finney Probit/Logit `LC50` studio (CLI + Ratatui TUI).
    Survival(SurvivalArgs),

    /// Tool #20: Leaf `FvCB` photosynthesis ($A\text{–}C_i$) & Odum (1956) open-channel stream `GPP`/`ER` metabolism studio (CLI + egui GUI).
    Flux(FluxArgs),
}

fn main() {
    let cli = OsabtCli::parse();
    let result = match cli.tool {
        Some(OsabtTool::Bench(args)) => run_bench_cli(args),
        Some(OsabtTool::Plasmid(args)) => run_plasmid_cli(args),
        Some(OsabtTool::Dwc(args)) => run_dwc_cli(args),
        Some(OsabtTool::Phylo(args)) => run_phylo_cli(args),
        Some(OsabtTool::Epi(args)) => run_epi_cli(args),
        Some(OsabtTool::Lexi(args)) => run_lexi_cli(args),
        Some(OsabtTool::Forge3d(args)) => run_forge3d_cli(args),
        Some(OsabtTool::Morpho(args)) => run_morpho_cli(args),
        Some(OsabtTool::Plate(args)) => run_plate_cli(args),
        Some(OsabtTool::Ecosound(args)) => run_ecosound_cli(args),
        Some(OsabtTool::Infoevo(args)) => run_infoevo_cli(args),
        Some(OsabtTool::Qpcr(args)) => run_qpcr_cli(args),
        Some(OsabtTool::Kinetics(args)) => run_kinetics_cli(args),
        Some(OsabtTool::Edna(args)) => run_edna_cli(args),
        Some(OsabtTool::Wildlife(args)) => run_wildlife_cli(args),
        Some(OsabtTool::Popgen(args)) => run_popgen_cli(args),
        Some(OsabtTool::Crispr(args)) => run_crispr_cli(args),
        Some(OsabtTool::Flow(args)) => run_flow_cli(args),
        Some(OsabtTool::Survival(args)) => run_survival_cli(args),
        Some(OsabtTool::Flux(args)) => run_flux_cli(args),
        None => run_bench_cli(BenchArgs { command: None }),
    };

    if let Err(err) = result {
        eprintln!("osabt error: {err}");
        process::exit(1);
    }
}
