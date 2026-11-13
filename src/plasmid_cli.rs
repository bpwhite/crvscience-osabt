// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #2: `osabt plasmid`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

use crate::gc_fraction;
use crate::gibson::scan_circular_bindings;
use crate::plasmid_digest::{
    design_goldengate_assembly, find_unique_cutters, simulate_restriction_digest, GoldenGateEnzyme,
};
use crate::plasmid_io::{
    align_trace_to_plasmid, find_orfs, load_plasmid_bytes, parse_ab1_trace, write_genbank,
    PlasmidRecord, SangerTrace,
};
use crate::plasmid_svg::{render_plasmid_map_svg, render_virtual_gel_svg, GelLadder, GelLane};
use crate::primer_design::parse_fasta_or_raw;
use crate::thermo::{PolymerasePreset, ReactionConditions};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliGoldenGateEnzyme {
    BsaI,
    BsmBI,
    BbsI,
    SapI,
}

impl From<CliGoldenGateEnzyme> for GoldenGateEnzyme {
    fn from(v: CliGoldenGateEnzyme) -> Self {
        match v {
            CliGoldenGateEnzyme::BsaI => GoldenGateEnzyme::BsaI,
            CliGoldenGateEnzyme::BsmBI => GoldenGateEnzyme::BsmBI,
            CliGoldenGateEnzyme::BbsI => GoldenGateEnzyme::BbsI,
            CliGoldenGateEnzyme::SapI => GoldenGateEnzyme::SapI,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliLadder {
    OneKbPlus,
    HundredBp,
}

impl From<CliLadder> for GelLadder {
    fn from(v: CliLadder) -> Self {
        match v {
            CliLadder::OneKbPlus => GelLadder::OneKbPlus,
            CliLadder::HundredBp => GelLadder::HundredBp,
        }
    }
}

/// Arguments for `osabt plasmid`. Running with no subcommand opens the native `egui` desktop GUI.
#[derive(Debug, Args)]
pub struct PlasmidArgs {
    #[command(subcommand)]
    pub command: Option<PlasmidCommand>,
}

#[derive(Debug, Subcommand)]
pub enum PlasmidCommand {
    /// Launch the interactive desktop vector map, cloning & virtual gel studio (`egui`).
    #[cfg(feature = "gui")]
    Gui {
        /// Optional path to a `.gbk`, `.gb`, `.dna`, or `.fa` plasmid file to open on startup.
        #[arg(short, long)]
        input: Option<PathBuf>,
    },

    /// Inspect a `.gbk`, `.dna`, or `.fa` construct (topology, features, 6-frame ORFs, and unique cutters).
    Inspect {
        /// Path to `.gbk`, `.gb`, `.dna`, or `.fa` file (omit to inspect the built-in `pCRV_DemoGFP` vector).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Optional path to export normalized GenBank (`.gbk`) file (useful for converting `.dna -> .gbk`).
        #[arg(long)]
        export_gbk: Option<PathBuf>,
    },

    /// Render a standalone publication circular plasmid map `.svg` (with optional bound primer stock `.tsv`).
    Map {
        /// Path to `.gbk`, `.dna`, or `.fa` construct (omit to render built-in demo construct).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Optional `.tsv`/`.csv` primer stock (`Name\tSequence`) to auto-bind onto the vector map.
        #[arg(short, long)]
        primers: Option<PathBuf>,
        /// Output `.svg` path.
        #[arg(short, long, default_value = "plasmid_map.svg")]
        out: PathBuf,
        /// Omit CRV Science Network footer attribution on exported SVG.
        #[arg(long)]
        no_attribution: bool,
    },

    /// Simulate single or multi-lane restriction digests and export a virtual agarose gel `.svg`.
    Digest {
        /// Path to `.gbk`, `.dna`, or `.fa` construct (omit to digest built-in demo construct).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Comma-separated enzyme list per lane; separate multiple lanes with `;` (e.g., `"EcoRI;BamHI;EcoRI,XhoI"`).
        #[arg(short, long, default_value = "EcoRI;XhoI;EcoRI,XhoI")]
        enzymes: String,
        /// Agarose percentage (`% w/v`, `0.6..=3.0`).
        #[arg(long, default_value_t = 1.0)]
        agarose: f64,
        /// Reference molecular-weight DNA ladder.
        #[arg(long, value_enum, default_value_t = CliLadder::OneKbPlus)]
        ladder: CliLadder,
        /// Render inverted black-on-white print gel instead of darkroom fluorescence.
        #[arg(long)]
        light_theme: bool,
        /// Optional path to save virtual agarose gel `.svg`.
        #[arg(short, long)]
        out: Option<PathBuf>,
    },

    /// Design a scarless Type IIS Golden Gate assembly with Potapov 4-bp fidelity & silent domestication audit.
    Goldengate {
        /// Multi-FASTA file of ordered dsDNA parts to assemble.
        #[arg(short, long)]
        input: PathBuf,
        #[arg(short, long, value_enum, default_value_t = CliGoldenGateEnzyme::BsaI)]
        enzyme: CliGoldenGateEnzyme,
        /// Assemble linear product rather than circular plasmid.
        #[arg(long)]
        linear: bool,
        /// Optional path to export Golden Gate tailed primers `.tsv`.
        #[arg(long)]
        export_tsv: Option<PathBuf>,
        /// Optional path to export assembled construct `.gbk`.
        #[arg(short, long)]
        out_gbk: Option<PathBuf>,
    },

    /// Align a Sanger `.ab1` chromatogram or FASTA sequencing read against a reference plasmid.
    AlignTrace {
        /// Reference plasmid (`.gbk`, `.dna`, or `.fa`).
        #[arg(short, long)]
        plasmid: PathBuf,
        /// Sanger `.ab1` trace or `.fa` read file.
        #[arg(short, long)]
        trace: PathBuf,
        /// Minimum Phred quality score (`Q`) for high-confidence variant calls.
        #[arg(long, default_value_t = 20)]
        min_q: u8,
    },
}

fn load_or_demo(input: Option<&PathBuf>) -> Result<PlasmidRecord, String> {
    match input {
        Some(p) => {
            let bytes =
                fs::read(p).map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
            load_plasmid_bytes(&p.to_string_lossy(), &bytes).map_err(|e| e.to_string())
        }
        None => Ok(PlasmidRecord::demo_pcrv_gfp()),
    }
}

fn load_oligo_tsv(path: &PathBuf) -> Result<Vec<(String, String)>, String> {
    let txt =
        fs::read_to_string(path).map_err(|e| format!("Failed reading '{}': {e}", path.display()))?;
    let mut out = Vec::new();
    for line in txt.lines() {
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let delim = if l.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = l.split(delim).map(|s| s.trim()).collect();
        if cols.len() >= 2
            && !cols[0].eq_ignore_ascii_case("name")
            && !cols[0].eq_ignore_ascii_case("oligo")
        {
            out.push((cols[0].to_string(), cols[1].to_string()));
        }
    }
    Ok(out)
}

/// Executes `osabt plasmid` CLI commands or launches the interactive `egui` desktop studio.
pub fn run_plasmid_cli(args: PlasmidArgs) -> Result<(), String> {
    let Some(cmd) = args.command else {
        #[cfg(feature = "gui")]
        {
            return crate::plasmid_gui::run_plasmid_gui(None);
        }
        #[cfg(not(feature = "gui"))]
        {
            return Err(
                "GUI feature is disabled in this build; pass a CLI subcommand (e.g., `osabt plasmid --help`)."
                    .to_string(),
            );
        }
    };

    match cmd {
        #[cfg(feature = "gui")]
        PlasmidCommand::Gui { input } => crate::plasmid_gui::run_plasmid_gui(input),

        PlasmidCommand::Inspect { input, export_gbk } => {
            let rec = load_or_demo(input.as_ref())?;
            let uniques = find_unique_cutters(&rec);
            let orfs = find_orfs(&rec, 40);

            println!(
                "OSABT Plasmid Summary — {} ({} bp, {}, {:.1}% GC)\n  Definition: {}\n",
                rec.name,
                rec.sequence.len(),
                if rec.circular { "Circular" } else { "Linear" },
                gc_fraction(&rec.sequence) * 100.0,
                rec.definition
            );
            println!("Annotated Features ({}):", rec.features.len());
            for f in &rec.features {
                println!(
                    "  • {:<18} {:<12} {}..={} ({:?})",
                    f.name,
                    f.feature_type,
                    f.start_0 + 1,
                    f.end_0 + 1,
                    f.strand
                );
            }
            println!("\nUnique Restriction Cutters ({}):", uniques.len());
            let u_strs: Vec<String> = uniques
                .iter()
                .map(|u| format!("{}({})", u.enzyme_name, u.cut_top_0 + 1))
                .collect();
            println!("  {}", u_strs.join(", "));

            println!("\n6-Frame ORFs >= 40 aa ({}):", orfs.len());
            for o in orfs.iter().take(10) {
                println!(
                    "  • Frame {:+2}: {}..={} ({} bp, {} aa)",
                    o.frame,
                    o.start_0 + 1,
                    o.end_0 + 1,
                    o.length_bp,
                    o.aa_sequence.len()
                );
            }

            if let Some(out_p) = export_gbk {
                fs::write(&out_p, write_genbank(&rec))
                    .map_err(|e| format!("Failed writing '{}': {e}", out_p.display()))?;
                println!("\nExported GenBank construct -> {}", out_p.display());
            }
            Ok(())
        }

        PlasmidCommand::Map {
            input,
            primers,
            out,
            no_attribution,
        } => {
            let rec = load_or_demo(input.as_ref())?;
            let bound = if let Some(p_path) = primers {
                let oligos = load_oligo_tsv(&p_path)?;
                let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);
                scan_circular_bindings(&rec.sequence, rec.circular, &oligos, 14, 1, &cond)
                    .map_err(|e| e.to_string())?
            } else {
                Vec::new()
            };

            let svg = render_plasmid_map_svg(&rec, &bound, !no_attribution);
            fs::write(&out, svg)
                .map_err(|e| format!("Failed writing '{}': {e}", out.display()))?;
            println!(
                "Rendered circular plasmid map ({} features, {} bound primers) -> {}",
                rec.features.len(),
                bound.len(),
                out.display()
            );
            Ok(())
        }

        PlasmidCommand::Digest {
            input,
            enzymes,
            agarose,
            ladder,
            light_theme,
            out,
        } => {
            let rec = load_or_demo(input.as_ref())?;
            let mut lanes = Vec::new();

            println!(
                "OSABT Restriction Digest — {} ({} bp, {:.1}% Agarose)\n",
                rec.name,
                rec.sequence.len(),
                agarose
            );
            for lane_spec in enzymes.split(';').map(str::trim).filter(|s| !s.is_empty()) {
                let enz_list: Vec<String> = lane_spec
                    .split(',')
                    .map(|s| s.trim().to_string())
                    .filter(|s| !s.is_empty())
                    .collect();
                let frags = simulate_restriction_digest(&rec, &enz_list)
                    .map_err(|e| e.to_string())?;
                let sizes: Vec<usize> = frags.iter().map(|f| f.length_bp).collect();
                println!("Lane [{}]: {:?} bp", lane_spec, sizes);
                lanes.push(GelLane {
                    label: lane_spec.to_string(),
                    fragment_lengths_bp: sizes,
                });
            }

            if let Some(out_p) = out {
                let svg =
                    render_virtual_gel_svg(ladder.into(), &lanes, agarose, !light_theme, true);
                fs::write(&out_p, svg)
                    .map_err(|e| format!("Failed writing '{}': {e}", out_p.display()))?;
                println!("\nSaved virtual agarose gel SVG -> {}", out_p.display());
            }
            Ok(())
        }

        PlasmidCommand::Goldengate {
            input,
            enzyme,
            linear,
            export_tsv,
            out_gbk,
        } => {
            let txt = fs::read_to_string(&input)
                .map_err(|e| format!("Failed reading '{}': {e}", input.display()))?;
            let parts = parse_fasta_or_raw(&txt).map_err(|e| e.to_string())?;
            let gg = design_goldengate_assembly(&parts, enzyme.into(), !linear)
                .map_err(|e| e.to_string())?;

            println!(
                "OSABT Golden Gate ({:?}, {}) | Predicted 4-bp Overhang Fidelity: {:.1}% | Illegal Internal Sites: {}\n",
                gg.enzyme,
                if gg.circular { "Circular" } else { "Linear" },
                gg.estimated_set_fidelity_pct,
                gg.internal_illegal_sites
            );
            for d in &gg.domestication_edits {
                println!(
                    "  [DOMESTICATE] {} codon at bp {}: {} -> {} (keeps '{}' AA)",
                    d.fragment_name,
                    d.codon_start_0 + 1,
                    d.original_codon,
                    d.silent_codon,
                    d.amino_acid
                );
            }
            for j in &gg.junctions {
                println!(
                    "  J{} ({} -> {}) | Overhang [{}]\n    FWD ({:.1} C): 5'-{}-3'\n    REV ({:.1} C): 5'-{}-3'",
                    j.junction_index + 1,
                    j.upstream_fragment,
                    j.downstream_fragment,
                    j.overhang_seq,
                    j.fwd_core_tm_c,
                    j.downstream_fwd_primer,
                    j.rev_core_tm_c,
                    j.upstream_rev_primer
                );
            }
            if let Some(tsv_p) = export_tsv {
                fs::write(&tsv_p, gg.to_order_tsv())
                    .map_err(|e| format!("Failed writing '{}': {e}", tsv_p.display()))?;
                println!("\nWrote Golden Gate primer order sheet -> {}", tsv_p.display());
            }
            if let Some(gbk_p) = out_gbk {
                let assembled_rec = PlasmidRecord {
                    name: "GoldenGate_Product".to_string(),
                    definition: format!("Assembled by OSABT Golden Gate ({:?})", gg.enzyme),
                    accession: "OSABT_GG".to_string(),
                    circular: gg.circular,
                    sequence: gg.assembled_sequence,
                    features: Vec::new(),
                };
                fs::write(&gbk_p, write_genbank(&assembled_rec))
                    .map_err(|e| format!("Failed writing '{}': {e}", gbk_p.display()))?;
                println!("Wrote assembled GenBank construct -> {}", gbk_p.display());
            }
            Ok(())
        }

        PlasmidCommand::AlignTrace {
            plasmid,
            trace,
            min_q,
        } => {
            let rec = load_or_demo(Some(&plasmid))?;
            let t_bytes = fs::read(&trace)
                .map_err(|e| format!("Failed reading '{}': {e}", trace.display()))?;
            let sanger = if t_bytes.starts_with(b"ABIF") {
                parse_ab1_trace(&t_bytes).map_err(|e| e.to_string())?
            } else {
                let txt = std::str::from_utf8(&t_bytes)
                    .map_err(|_| "Trace is neither ABIF (.ab1) nor UTF-8 FASTA.".to_string())?;
                let (sample_name, sequence) = parse_fasta_or_raw(txt)
                    .map_err(|e| e.to_string())?
                    .into_iter()
                    .next()
                    .ok_or_else(|| "Empty FASTA trace file.".to_string())?;
                let qualities = vec![40u8; sequence.len()];
                SangerTrace {
                    sample_name,
                    sequence,
                    qualities,
                }
            };

            let aln = align_trace_to_plasmid(&rec, &sanger, min_q)
                .ok_or_else(|| "No seed alignment found between trace and plasmid.".to_string())?;
            println!(
                "Aligned '{}' ({:?}) -> {} bp {}..={} ({} bp) | Identity: {:.2}% | Mean Q: {:.1}",
                sanger.sample_name,
                aln.strand,
                rec.name,
                aln.plasmid_start_0 + 1,
                aln.plasmid_end_0 + 1,
                aln.aligned_len_bp,
                aln.identity_pct,
                aln.mean_phred_q
            );
            for v in &aln.variants {
                println!(
                    "  bp {:>5}: {} -> {} (Q{}) [{}]",
                    v.plasmid_pos_0 + 1,
                    v.ref_base,
                    v.read_base,
                    v.read_phred_q,
                    v.classification
                );
            }
            Ok(())
        }
    }
}
