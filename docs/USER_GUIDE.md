# OSABT User Guide (`v2.0.0`)

**Open Source Ad-hoc Biology Tools (`osabt`)** is a pure-Rust, GPL-3.0-or-later workbench combining 20 wet-lab, molecular, evolutionary, kinetic, ecological, and biostatistical tools into a single unified desktop GUI (`eframe`/`egui`), terminal TUI suite (`Ratatui`), command-line interface (`clap`), and library crate.

---

## 1. Installation & Feature Matrix

### Build Requirements
- **Rust toolchain**: `1.75` or newer (`rust-version = "1.75"`).
- **Linux system libraries (for desktop GUI builds only)**:
  ```bash
  sudo apt-get update
  sudo apt-get install -y libxkbcommon-dev libwayland-dev libgtk-3-dev
  ```

### Feature Flags (`Cargo.toml`)
- `default = ["cli", "tui", "gui"]`: Builds the complete 20-tool `osabt` binary with CLI, Ratatui TUIs, and `eframe`/`egui` desktop studios.
- `cli`: Enables `clap` v4 subcommands across all 20 tools (`--no-default-features --features cli` for lightweight headless servers).
- `tui`: Enables `Ratatui` + `crossterm` interactive terminal studios (`--no-default-features --features tui`, automatically enables `cli`).
- `gui`: Enables `eframe`/`egui` desktop studios and offscreen PNG screenshots (`--no-default-features --features gui`).
- `serde`: Derives `serde::Serialize` and `serde::Deserialize` on core report and result structs.

```bash
# Full desktop + TUI + CLI build
cargo build --release

# Headless CLI-only server build
cargo build --release --no-default-features --features cli

# Run the complete unit, headless GUI, integration, and validation test suite
cargo test --all-features --all-targets
```

---

## 2. Launching OSABT

| Mode | Command | Description |
| :--- | :--- | :--- |
| Unified 20-Tool Desktop Workbench | `osabt` | Opens the unified `eframe`/`egui` suite launcher with domain filters, search, navigation rail, 12 integrated panels, and all 8 studios. |
| Terminal Workbench (Tool #1) | `osabt --tui` | Opens the interactive 4-tab Ratatui `bench` TUI directly. |
| Standalone Studio GUI | `osabt <tool> gui` (or `osabt <tool>`) | Launches a dedicated GUI studio (`plasmid`, `forge3d`, `morpho`, `ecosound`, `kinetics`, `popgen`, `flow`, `flux`). |
| Legacy Flux Studio Alias | `osabt gui [--aci <CSV>] [--stream <CSV>]` | Launches `OSABT Flux Studio` directly with optional input CSVs. |
| Terminal TUI Studio | `osabt <tool> tui` (or `osabt <tool>`) | Launches a Ratatui terminal studio (`bench`, `dwc`, `phylo`, `epi`, `lexi`, `plate`, `infoevo`, `qpcr`, `edna`, `wildlife`, `crispr`, `survival`). |
| Offscreen Headless Screenshot | `osabt --screenshot out.png --screenshot-view <1\|2\|3>` | Renders a `1600x1000` PNG of Workbench view `1` (Launcher), `2` (Survival), or `3` (Plasmid) and exits. |

---

## 3. Controls & Keyboard Shortcuts

### Unified Desktop Workbench (`osabt`)
| Shortcut / Control | Scope | Action |
| :--- | :--- | :--- |
| `F1` | Global (when not typing) | Toggle the **About OSABT** modal (version, license, CRV Science links). |
| `Esc` | Global (context-aware) | 1) Surrenders focus if a text box is active; 2) Closes **About OSABT** or paste modal if open; 3) Returns from an active tool to the **Suite Home Launcher**. |
| `Ctrl+F` (or `Cmd+F` on macOS) | Global (when not typing) | Focus the top-bar tool search box (`Search 20 tools...`). |
| `Ctrl+B` (or `Cmd+B` on macOS) | Global (when not typing) | Show or hide the left `#01..#20` tool navigation rail. |
| **Reset Filter** button | Home Launcher empty state | Clears active search query and domain filter when `No matching tools` appears. |

### `osabt forge3d` (3D Print Studio)
| Control / Shortcut | Action |
| :--- | :--- |
| **Left-mouse drag** on 3D viewport | Orbit camera (`yaw` and `pitch`, clamped `[-1.45, +1.45]` rad). |
| **Middle-mouse drag** on 3D viewport | Pan camera horizontally/vertically across the viewport. |
| **Mouse scroll wheel** | Zoom in/out (`0.35x` to `3.50x`). |
| `R` (when not typing) or **Reset Camera** button | Reset orbit angles, pan offset, and zoom to the default print-bed view. |
| **Perspective 3D** checkbox | Toggle perspective projection vs. orthographic CAD projection. |
| **Wireframe** checkbox | Overlay triangle edges over Blinn-Phong + rim + depth-cued shading. |
| Parameter sliders (`Lobes`, `Areolae Pore Pitch`, `Voxel Grid Resolution`, `Max Bed Dimension`) | Heavy voxel remeshing runs once on slider release (`drag_stopped`) or direct numeric entry so dragging stays smooth. |

### `osabt ecosound` (Bioacoustics Studio)
| Control / Shortcut | Action |
| :--- | :--- |
| **Left-mouse drag** on spectrogram | Draw a new manual Raven selection bounding box and measure peak frequency, centroid, and SNR. |
| **Left-mouse click** on box or table row | Select vocalization event (`selected_event_idx`). |
| `Delete` or `Backspace` (when not typing) | Remove the currently selected vocalization box and renumber selections. |
| **Right-mouse click** on box | Delete that bounding box immediately and renumber selections. |

### Ratatui Terminal Studios (`bench`, `dwc`, `phylo`, `epi`, `lexi`, `plate`, `infoevo`, `qpcr`, `edna`, `wildlife`, `crispr`, `survival`)
| Shortcut | Action |
| :--- | :--- |
| `Tab` / `BackTab` or `1`–`4` | Switch between studio tabs. |
| `Up` / `Down` or `j` / `k` | Navigate list/table rows. |
| `Left` / `Right` or `h` / `l` | Adjust active parameter or scroll columns. |
| `F1` or `?` | Display OSABT version & license status line. |
| `q` or `Esc` | Close active inline edit prompt or exit the TUI cleanly (terminal state restored via RAII `TerminalSession`). |

---

## 4. Common CLI Workflows

Every CLI subcommand ships with a built-in realistic biological demo dataset when input flags are omitted, making rapid verification and pipeline scripting straightforward.

### Wet Lab, Cloning & Sequences
```bash
# #1 Primer design + TaqMan probe + vendor TSV
osabt bench design --input template.fasta --amplicon 140..360 --probe --tsv-out primers.tsv
osabt bench calc buffer --system tris --ph 8.0 --molarity 0.05 --liters 1.0 --temp-c 25

# #2 Plasmid map, virtual agarose gel & Golden Gate assembly
osabt plasmid map --input vector.gbk --primers oligos.tsv --out vector_map.svg
osabt plasmid digest --input vector.gbk --enzymes "EcoRI;XhoI;EcoRI,XhoI" --out gel.svg
osabt plasmid goldengate --input parts.fasta --enzyme bsa-i --out-gbk assembled.gbk

# #17 CRISPR multi-PAM design, CFD off-target scoring & base editing
osabt crispr design --nuclease spcas9 --tsv-out guides.tsv --svg-out locus.svg
osabt crispr cfd <20nt-SPACER> <23nt-TARGET+PAM>
```

### Morphology, Bioacoustics & 3D Printing
```bash
# #7 Generate watertight diatom frustule or PDB surface STL/3MF
osabt forge3d diatom -m 6 --pore-pitch 0.28 --resolution 52 --scale-mm 95 -o diatom.3mf
osabt forge3d pdb --input structure.pdb --resolution 48 -o molecule.stl
osabt forge3d inspect-stl diatom.3mf

# #8 CFU counting, leaf herbivory & 2D Procrustes morphometrics
osabt morpho cfu --input plate.bmp --volume-ml 0.1 --dilution 1e-5 --svg-out cfu.svg
osabt morpho leaf --input leaf.bmp --px-per-mm 4.0
osabt morpho gpa --input landmarks.tps --svg-out gpa.svg

# #10 Ecoacoustic indices & Raven selection table + call slicing
osabt ecosound indices --wav dawn_chorus.wav --fft-size 1024
osabt ecosound slice --wav dawn_chorus.wav --out-dir ecosound_exports
```

### Kinetics, qPCR, Plates & Flow Cytometry
```bash
# #9 96/384-well assay layout & 4PL IC50 fitting
osabt plate plan --wells 96 --moat 1 --seed 42 --out-prefix assay_p1
osabt plate fit --input dose_response.csv --svg-out ic50.svg

# #12 MIQE RT-qPCR relative quantification & melt curves
osabt qpcr quantify --input qpcr.csv --refs ACTB,GAPDH --calibrator Control --svg-out fc.svg
osabt qpcr melt --melt melt.csv --dimer-cutoff-c 77.5 --svg-out melt.svg

# #13 Global enzyme inhibition AICc & Morrison tight-binding fit
osabt kinetics inhibit --input matrix.csv --svg-out inhibition.svg
osabt kinetics morrison --input tight_binding.csv --nominal-et-nm 25.0 --substrate-um 70.0 --km-um 35.0

# #18 Flow cytometry FCS 3.1 spillover compensation & Logicle quadrant gates
osabt flow compensate --input sample.fcs --output compensated.fcs
osabt flow gate --input compensated.fcs --x-ch 2 --y-ch 3 --svg-out biaxial.svg
```

### Evolution, Ecology, Epidemiology, Survival & Flux
```bash
# #3 Darwin Core QA & DwC-A packaging
osabt dwc audit --input survey.csv --report-tsv qa.tsv
osabt dwc archive --input survey.csv --out-dir dwca_bundle

# #4 Newick / NEXUS phylogenetics
osabt phylo view --input tree.nwk
osabt phylo svg --input tree.nwk --circular --output cladogram.svg

# #5 Epidemiology Cori R_t & SEIR simulation
osabt epi rt --input daily_cases.csv --si-mean 4.7 --si-sd 2.9 --svg-out rt.svg
osabt epi seir --pop 100000 --r0 2.8 --incubation 3.5 --infectious 5.0 --svg-out seir.svg

# #6 Biological Greek/Latin etymology & Anki deck export
osabt lexi decode "Pseudotsuga menziesii" "Oncorhynchus mykiss"
osabt lexi flashcards --anki-out bio_roots_anki.tsv

# #11 Information-theoretic evolution (R_seq, MIp, Kimura)
osabt infoevo logo --input alignment.fasta --svg-out logo.svg
osabt infoevo epistasis --input alignment.fasta --svg-out mip.svg

# #14 eDNA FASTQ demultiplexing, UNOISE denoising & rarefaction
osabt edna demux --fastq reads.fastq --barcodes sheet.csv --max-ee 1.0 --clean-only --tsv-out asv.tsv
osabt edna rarefy --fastq reads.fastq --barcodes sheet.csv --svg-out rarefaction.svg

# #15 Wildlife occupancy & mark-recapture
osabt wildlife occupancy --input occupancy.csv --svg-out occupancy.svg
osabt wildlife markrecap --input captures.csv

# #16 Population genetics & Median-Joining haplotype networks
osabt popgen diversity --fasta alignment.fasta
osabt popgen vcf --vcf genotypes.vcf
osabt popgen network --fasta alignment.fasta --svg-out haplotype_network.svg

# #19 Kaplan-Meier, Cox PH & acute Probit/Logit LC50
osabt survival km --input cohort.csv
osabt survival lc50 --input bioassay.csv --link probit --svg-out lc50.svg

# #20 Leaf FvCB photosynthesis & Odum stream metabolism
osabt flux aci --input aci_curve.csv --svg-out flux_report.svg
osabt flux stream --input stream_do.csv --svg-out flux_report.svg
```

---

## 5. File Formats & Troubleshooting

- **UTF-8 BOM Support**: All text, CSV/TSV, FASTA, FASTQ, GenBank, PDB, Newick, and `.tps` parsers automatically strip leading UTF-8 BOMs (`\u{feff}`) emitted by Excel and Windows tools.
- **CSV / TSV Auto-Detection**: Tabular parsers detect `,` vs `\t` automatically, ignore blank lines and `#` comments, and handle missing optional columns with documented defaults.
- **Non-Finite Numerical Inputs**: Kinetic curve fitters (`kinetics_core`), leaf `FvCB` evaluation (`flux_leaf`, guarded down to absolute zero `<= -273.15 °C`), and stream DO fitting (`flux_stream`) filter or reject `NaN`/`Inf` values cleanly without panicking.
- **Headless Servers Without Display/GPU**: Build with `cargo build --release --no-default-features --features cli` if X11/Wayland/Vulkan/OpenGL libraries are absent.
- **Windows Console Output**: Because the release GUI binary uses `#![windows_subsystem = "windows"]`, `main.rs` calls `AttachConsole(ATTACH_PARENT_PROCESS)` whenever CLI flags or subcommands are passed so stdout/stderr appear in PowerShell and `cmd.exe`.
