# Changelog

All notable releases of **Open Source Ad-hoc Biology Tools (`osabt`)** are documented in this file.

## [2.0.0] — 2026-11-12 (Phase 3 Complete)

### Added
- **Tool #16 (`osabt popgen`):** Tajima's (1989) `D`, nucleotide diversity (`pi`, `theta_W`), `d_XY` / Hudson `F_ST`, diploid VCF Weir–Cockerham (1984) `F_IS`/`F_ST`/`F_IT`, Wigginton et al. (2005) exact `HWE`, and Bandelt et al. (1999) Median-Joining haplotype networks (`egui` GUI + CLI).
- **Tool #17 (`osabt crispr`):** Multi-PAM (`SpCas9`, `SaCas9`, `LbCas12a`) guide scanner, Doench et al. (2016) Rule Set 2 & `CFD` off-target scoring, Sugimoto RNA:DNA R-loop `Tm`, `CBE`/`ABE` bystander window auditing, and Golden Gate oligo synthesis (`Ratatui` TUI + CLI).
- **Tool #18 (`osabt flow`):** Binary ISAC `FCS 3.1` parser & writer (Spidlen et al. 2010), `$SPILLOVER` matrix inversion, Parks et al. (2006) `Logicle` biexponential scaling, and 2D polygon/quadrant gating (`egui` GUI + CLI).
- **Tool #19 (`osabt survival`):** Kaplan–Meier (1958) product-limit curves with Greenwood SE & log-log 95% CIs, Mantel (1966) log-rank test, Cox (1972) proportional-hazards regression, and Finney Probit/Logit acute `LC50` (`Ratatui` TUI + CLI).
- **Tool #20 (`osabt flux`):** Farquhar–von Caemmerer–Berry (1980) `FvCB` $A\text{–}C_i$ photosynthesis fitter (`Vcmax`, `Jmax`, `TPU`, `Rd`) and Odum (1956) open-channel diel `DO` stream metabolism (`GPP`, `ER`, `NEP`, `K600`) (`egui` GUI + CLI).
- Repository packaging: `README.md`, `LICENSE` (`GPL-3.0-or-later`), `.github/workflows/ci.yml`, `CHANGELOG.md`, `CITATION.cff`, and `.gitignore`.

## [1.5.0] — Phase 2 Complete
- **Tool #11 (`osabt infoevo`):** Information-theoretic evolution (`ITME`) sequence logos, Adami physical complexity, `MIp` epistasis, and Kimura regret (`Ratatui` TUI + CLI).
- **Tool #12 (`osabt qpcr`):** MIQE RT-qPCR relative quantification (`qBase`/`Pfaffl`/`2^-ddCq`), `geNorm` stability, and SYBR `-dF/dT` melt curves (`Ratatui` TUI + CLI).
- **Tool #13 (`osabt kinetics`):** Michaelis–Menten/Hill kinetics, `AICc` reversible `Ki` selector, and Morrison (1969) tight-binding solver (`egui` GUI + CLI).
- **Tool #14 (`osabt edna`):** Dual-barcode FASTQ demultiplexer, `maxEE` filter, `UNOISE` 1-nt denoiser, 2-parent chimera & blank auditor, `Chao1`/Hill/Hurlbert rarefaction (`Ratatui` TUI + CLI).
- **Tool #15 (`osabt wildlife`):** MacKenzie (2002) occupancy MLE, survey `K*`/`S*` planner, Chapman–Petersen/`M0`/`Mt`/`Chao Mh`/`CJS` capture-recapture (`Ratatui` TUI + CLI).

## [1.0.0] — Phase 1 Complete
- **Tools #1–#10:** `bench`, `plasmid`, `dwc`, `phylo`, `epi`, `lexi`, `forge3d`, `morpho`, `plate`, and `ecosound`.
