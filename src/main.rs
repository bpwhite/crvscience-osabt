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
// Open Source Ad-hoc Biology Tools (OSABT v2.0.0) — Unified CLI, TUI & GUI Entrypoint

#![deny(unsafe_code)]
#![cfg_attr(all(not(debug_assertions), feature = "gui"), windows_subsystem = "windows")]

use std::io::Write as _;
use std::path::PathBuf;
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
    /// Launch the interactive Ratatui terminal workbench (Tool #1) instead of the desktop GUI.
    #[arg(long)]
    tui: bool,

    /// Render an offscreen 1600x1000 GPU screenshot of OSABT Workbench GUI to PNG and exit.
    #[arg(long, visible_alias = "headless-screenshot", value_name = "PNG_PATH")]
    screenshot: Option<PathBuf>,

    /// Select which showcase view to capture (1 = Suite Launcher, 2 = Survival & Ecotox, 3 = Plasmid Studio).
    #[arg(long, value_name = "VIEW", default_value_t = 1)]
    screenshot_view: u8,

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

    /// Tool #20: Leaf `FvCB` photosynthesis (A–Ci) & Odum (1956) open-channel stream `GPP`/`ER` metabolism studio (CLI + egui GUI).
    Flux(FluxArgs),

    /// Launch OSABT Flux Studio GUI (or capture `--screenshot <PNG_PATH>` headlessly).
    Gui {
        /// Optional Leaf A-Ci CSV file path.
        #[arg(long)]
        aci: Option<PathBuf>,
        /// Optional Odum open-channel stream DO CSV file path.
        #[arg(long)]
        stream: Option<PathBuf>,
        /// Render an offscreen 1600x1000 GPU screenshot to PNG and exit without opening a window.
        #[arg(long, value_name = "PNG_PATH")]
        screenshot: Option<PathBuf>,

        /// Select which showcase view to capture (1 = Dual Overview, 2 = Survival & Ecotox, 3 = Kinetics & Plasmid).
        #[arg(long, value_name = "VIEW", default_value_t = 1)]
        screenshot_view: u8,
    },
}

#[cfg(windows)]
fn attach_parent_console_if_cli() {
    if std::env::args_os().len() > 1 {
        extern "system" {
            fn AttachConsole(dwProcessId: u32) -> i32;
        }
        const ATTACH_PARENT_PROCESS: u32 = u32::MAX;
        // SAFETY: Win32 AttachConsole(ATTACH_PARENT_PROCESS) attaches stdout/stderr to the invoking
        // console when CLI arguments are passed to a `windows_subsystem = "windows"` binary.
        #[allow(unsafe_code)]
        unsafe {
            let _ = AttachConsole(ATTACH_PARENT_PROCESS);
        }
    }
}

fn flush_std_streams() {
    let _ = std::io::stdout().flush();
    let _ = std::io::stderr().flush();
}

fn main() {
    #[cfg(windows)]
    attach_parent_console_if_cli();

    let cli = OsabtCli::parse();
    if let Some(out_path) = cli.screenshot {
        #[cfg(feature = "gui")]
        {
            if let Err(err) =
                osabt::workbench_gui::run_workbench_screenshot(&out_path, cli.screenshot_view)
            {
                eprintln!("Error: {err}");
                flush_std_streams();
                process::exit(1);
            }
            flush_std_streams();
            return;
        }
        #[cfg(not(feature = "gui"))]
        {
            let _ = out_path;
            eprintln!("Error: --screenshot requires building with `--features gui`.");
            flush_std_streams();
            process::exit(1);
        }
    }
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
        Some(OsabtTool::Gui {
            aci,
            stream,
            screenshot,
            screenshot_view,
        }) => {
            #[cfg(feature = "gui")]
            {
                if let Some(out_path) = screenshot {
                    osabt::flux_gui::run_flux_screenshot(aci, stream, &out_path, screenshot_view)
                } else {
                    osabt::flux_gui::run_flux_gui(aci, stream)
                }
            }
            #[cfg(not(feature = "gui"))]
            {
                let _ = (aci, stream, screenshot, screenshot_view);
                Err("Built without `gui` feature; rebuild with `--features gui`.".to_string())
            }
        }
        None => {
            if cli.tui {
                run_bench_cli(BenchArgs { command: None })
            } else {
                #[cfg(feature = "gui")]
                {
                    osabt::workbench_gui::run_workbench_gui()
                }
                #[cfg(not(feature = "gui"))]
                {
                    run_bench_cli(BenchArgs { command: None })
                }
            }
        }
    };

    if let Err(err) = result {
        eprintln!("osabt error: {err}");
        flush_std_streams();
        process::exit(1);
    }
    flush_std_streams();
}
