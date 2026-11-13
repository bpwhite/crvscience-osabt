# Open Source Ad-hoc Biology Tools (`osabt` `v2.0.0`)

Unified, zero-paywall **GNU GPL-3.0-or-later** Rust library, CLI, `Ratatui` TUI, and `eframe`/`egui` desktop/WASM GUI suite for the **[CRV Science Network](https://crvscience.com)**.

## Quick Start

```bash
# Build and run all unit & integration tests across the 20-tool suite
cargo test --all-features

# Launch any TUI or GUI tool directly with built-in demo datasets
cargo run -- bench tui
cargo run -- plasmid gui
cargo run -- popgen gui
cargo run -- flux gui
```

### Feature Flags

- `default = ["cli", "tui", "gui"]` — Builds the unified `osabt` binary with all 20 subcommands.
- `cli` — Enables `clap` v4 command-line subcommands (`--no-default-features --features cli` for headless servers).
- `tui` — Enables interactive `Ratatui` + `crossterm` terminal studios.
- `gui` — Enables `eframe` / `egui` native desktop & WASM visual studios.
- `serde` — Derives `Serialize` / `Deserialize` on core data structures and reports.

## Complete 20-Tool Suite (`v2.0.0`)

### Phase 1 (`v1.0.0` — Wet Lab, Cloning, Biodiversity & Bioacoustics)

| # | Subcommand | Interface | Core Capabilities |
| :--- | :--- | :--- | :--- |
| 1 | `osabt bench` | CLI + TUI | SantaLucia nearest-neighbor `Tm`, secondary structure, primer designer, Gibson overlaps, wet-lab calculators |
| 2 | `osabt plasmid` | CLI + GUI | GenBank/FASTA vector I/O, restriction digest, virtual agarose gel, `.svg` plasmid map |
| 3 | `osabt dwc` | CLI + TUI | Darwin Core (`DwC-A`) linter, coordinate scrubber, taxonomic synonym resolver |
| 4 | `osabt phylo` | CLI + TUI | Newick/NEXUS tree explorer, rerooting, monophyly checks, Robinson–Foulds distance, `.svg` cladograms |
| 5 | `osabt epi` | CLI + TUI | Renewal `R_t`, SIR/SEIR ODEs, diagnostic Bayes, cohort & test-negative power |
| 6 | `osabt lexi` | CLI + TUI | Greek/Latin taxonomic & anatomical morpheme decoder, cognate trees, Anki export |
| 7 | `osabt forge3d` | CLI + GUI | PDB biomolecule, Gielis diatom, radiolarian, and voxel-to-manifold `.stl` / `.3mf` 3D print generator |
| 8 | `osabt morpho` | CLI + GUI | Automated Petri CFU counter, Excess-Green leaf herbivory, 2D Generalized Procrustes (`GPA`) morphometrics |
| 9 | `osabt plate` | CLI + TUI | 96/384-well randomized assay planner, Zhang `Z'`-factor QC, 4PL `IC50`/`EC50` fitter |
| 10 | `osabt ecosound` | CLI + GUI | STFT spectrogram viewer, `ACI`/`NDSI`/`H`/`ADI` acoustic indices, Raven annotator, `.wav` call slicer |

### Phase 2 (`v1.5.0` — Molecular Evolution, qPCR, Kinetics, eDNA & Wildlife)

| # | Subcommand | Interface | Core Capabilities |
| :--- | :--- | :--- | :--- |
| 11 | `osabt infoevo` | CLI + TUI | ITME sequence logos (`R_seq`), Adami physical complexity, `MIp` epistasis, Kimura regret |
| 12 | `osabt qpcr` | CLI + TUI | MIQE RT-qPCR relative quantification (`qBase`/`Pfaffl`/`2^-ddCq`), `geNorm` stability, SYBR `-dF/dT` melt curves |
| 13 | `osabt kinetics` | CLI + GUI | Michaelis–Menten/Hill fits, `AICc` reversible `Ki` selector, Morrison (1969) tight-binding solver |
| 14 | `osabt edna` | CLI + TUI | Dual-barcode FASTQ demultiplexer, `maxEE` filter, `UNOISE` 1-nt denoiser, 2-parent chimera audit, `Chao1`/Hill/Hurlbert rarefaction |
| 15 | `osabt wildlife` | CLI + TUI | MacKenzie (2002) occupancy MLE, survey `K*`/`S*` planner, Chapman–Petersen/`M0`/`Mt`/`Chao Mh`/`CJS` mark-recapture |

### Phase 3 (`v2.0.0` — PopGen, CRISPR, Cytometry, Survival & Ecosystem Flux)

| # | Subcommand | Interface | Core Capabilities |
| :--- | :--- | :--- | :--- |
| 16 | `osabt popgen` | CLI + GUI | Tajima's `D`, `pi`/`theta_W`, `d_XY`/Hudson `F_ST`, diploid VCF Weir–Cockerham `F_ST`, Wigginton exact `HWE`, Bandelt Median-Joining networks |
| 17 | `osabt crispr` | CLI + TUI | `SpCas9`/`SaCas9`/`LbCas12a` scanning, Doench Rule Set 2 & `CFD` off-target scoring, `CBE`/`ABE` bystander windows, Golden Gate oligos |
| 18 | `osabt flow` | CLI + GUI | Binary ISAC `FCS 3.1` parser/writer, `$SPILLOVER` inversion, Parks (2006) `Logicle` scaling, polygon/quadrant gating |
| 19 | `osabt survival` | CLI + TUI | Kaplan–Meier product-limit curves, Mantel log-rank test, Cox (1972) PH regression, Finney Probit/Logit acute `LC50` |
| 20 | `osabt flux` | CLI + GUI | Farquhar–von Caemmerer–Berry (`FvCB`) $A\text{–}C_i$ curve fitting (`Vcmax`, `Jmax`, `TPU`, `Rd`) & Odum (1956) diel `DO` stream metabolism |

## License

Copyright (C) 2026 Bryan P. White — CRV Science Network.  
Licensed under the **GNU General Public License v3.0 or later** (`GPL-3.0-or-later`). See [`LICENSE`](LICENSE) and <https://www.gnu.org/licenses/gpl-3.0.txt>.
