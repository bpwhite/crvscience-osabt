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

//! Command-line subcommands (`clap` v4) for Tool #1: `osabt bench`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

use crate::bench_calc::{
    build_pcr_master_mix, calc_buffer_recipe, calc_c1v1_dilution, calc_molarity_grams,
    calc_oligo_resuspension, calc_serial_dilution, parse_vendor_coa_table, rcf_to_rpm, rpm_to_rcf,
    BufferSystem,
};
use crate::gibson::{design_gibson_assembly, GibsonConfig, GibsonFragmentInput};
use crate::primer_design::{
    design_multiplex_pool, design_pcr_primers, parse_fasta_or_raw, scan_offtarget_amplicons,
    MultiplexTarget, PrimerDesignConfig,
};
use crate::structure::{evaluate_oligo_qc, score_dimer};
use crate::thermo::{
    calc_degenerate_tm_bounds, recommend_annealing_temp, PolymerasePreset, ReactionConditions,
};
use crate::{is_unambiguous_dna, normalize_dna};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliPolymerase {
    Q5,
    Phusion,
    Taq,
    Onetaq,
}

impl From<CliPolymerase> for PolymerasePreset {
    fn from(val: CliPolymerase) -> Self {
        match val {
            CliPolymerase::Q5 => PolymerasePreset::Q5,
            CliPolymerase::Phusion => PolymerasePreset::Phusion,
            CliPolymerase::Taq => PolymerasePreset::Taq,
            CliPolymerase::Onetaq => PolymerasePreset::OneTaq,
        }
    }
}

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliBufferSystem {
    Tris,
    Hepes,
    Mes,
    Phosphate,
}

impl From<CliBufferSystem> for BufferSystem {
    fn from(val: CliBufferSystem) -> Self {
        match val {
            CliBufferSystem::Tris => BufferSystem::Tris,
            CliBufferSystem::Hepes => BufferSystem::Hepes,
            CliBufferSystem::Mes => BufferSystem::Mes,
            CliBufferSystem::Phosphate => BufferSystem::SodiumPhosphate,
        }
    }
}

/// Arguments for `osabt bench`. Running with no subcommand opens the interactive `Ratatui` TUI.
#[derive(Debug, Args)]
pub struct BenchArgs {
    #[command(subcommand)]
    pub command: Option<BenchCommand>,
}

#[derive(Debug, Subcommand)]
pub enum BenchCommand {
    /// Launch the interactive 4-tab terminal workbench (`Ratatui`).
    #[cfg(feature = "tui")]
    Tui,

    /// Design ranked PCR or hydrolysis-probe (`TaqMan`) qPCR primer pairs from a FASTA or raw DNA sequence.
    Design {
        /// Path to template FASTA file, or an inline raw DNA string.
        #[arg(short, long)]
        input: String,
        /// Amplicon length range formatted as `MIN..MAX` (e.g., `120..500`).
        #[arg(long, default_value = "120..500")]
        amplicon: String,
        /// Optional 0-based target region formatted as `START..END` that amplicons must span.
        #[arg(long)]
        region: Option<String>,
        /// Target primer `Tm` in Celsius.
        #[arg(long, default_value_t = 60.0)]
        tm: f64,
        /// Polymerase reaction buffer preset.
        #[arg(short, long, value_enum, default_value_t = CliPolymerase::Q5)]
        polymerase: CliPolymerase,
        /// Also pick an internal hydrolysis (`TaqMan`) qPCR probe (`Tm + 9 C`, no 5'-G).
        #[arg(long)]
        qpcr_probe: bool,
        /// Optional background FASTA file for off-target amplicon / NUMT screening.
        #[arg(long)]
        offtarget: Option<PathBuf>,
        /// Number of top diverse primer pairs to return.
        #[arg(short, long, default_value_t = 5)]
        num: usize,
        /// Optional path to write vendor-ready bulk oligo order `.tsv`.
        #[arg(long)]
        export_tsv: Option<PathBuf>,
    },

    /// Audit an existing primer or forward/reverse pair (thermodynamics, hairpins, 3' dimers, G4Hunter, off-targets).
    Check {
        /// Forward primer sequence (`5' -> 3'`, supports IUPAC degenerate codes).
        #[arg(short, long)]
        fwd: String,
        /// Optional reverse primer sequence (`5' -> 3'`).
        #[arg(short, long)]
        rev: Option<String>,
        /// Polymerase reaction buffer preset.
        #[arg(short, long, value_enum, default_value_t = CliPolymerase::Q5)]
        polymerase: CliPolymerase,
        /// Optional background FASTA file to scan for unintended PCR products.
        #[arg(long)]
        offtarget: Option<PathBuf>,
    },

    /// Design a conflict-free multiplex PCR primer pool across a multi-FASTA file of target loci.
    Multiplex {
        /// Path to multi-FASTA file with one record per target locus.
        #[arg(short, long)]
        input: PathBuf,
        /// Amplicon length range formatted as `MIN..MAX` (e.g., `120..450`).
        #[arg(long, default_value = "120..450")]
        amplicon: String,
        /// Minimum allowed 3'-extensible cross-dimer `Delta G` (`kcal/mol`, e.g., `-5.5`).
        #[arg(long, default_value_t = -5.5)]
        dg_floor: f64,
        /// Minimum bp difference between amplicons for single-lane agarose gel resolution (`0` to disable).
        #[arg(long, default_value_t = 30)]
        size_spacing: usize,
        /// Polymerase preset.
        #[arg(short, long, value_enum, default_value_t = CliPolymerase::Q5)]
        polymerase: CliPolymerase,
        /// Optional path to export vendor order `.tsv`.
        #[arg(long)]
        export_tsv: Option<PathBuf>,
    },

    /// Design two-stage isothermal Gibson / NEBuilder HiFi assembly primers across multi-FASTA fragments.
    Gibson {
        /// Path to multi-FASTA file containing ordered dsDNA fragments (`>= 2` linear, `>= 1` circular).
        #[arg(short, long)]
        input: PathBuf,
        /// Assemble into a closed circular plasmid (`true` by default; pass `--linear` for linear product).
        #[arg(long)]
        linear: bool,
        /// Target 3' core annealing `Tm` (`C`).
        #[arg(long, default_value_t = 60.0)]
        core_tm: f64,
        /// Polymerase preset for fragment amplification.
        #[arg(short, long, value_enum, default_value_t = CliPolymerase::Q5)]
        polymerase: CliPolymerase,
        /// Optional path to export chimeric primer order `.tsv`.
        #[arg(long)]
        export_tsv: Option<PathBuf>,
    },

    /// Calculate lyophilized oligo resuspension volumes (`uL` TE buffer) from inline args or a vendor COA `.csv`/`.tsv`.
    Resuspend {
        /// Optional vendor COA `.csv`/`.tsv` file (`Name,Sequence,nmol`).
        #[arg(short, long)]
        csv: Option<PathBuf>,
        /// Single oligo name (when `--csv` is not provided).
        #[arg(long, default_value = "Oligo_1")]
        name: String,
        /// Single oligo DNA sequence (when `--csv` is not provided).
        #[arg(short, long)]
        seq: Option<String>,
        /// Delivered synthesis yield in nanomoles (`nmol`).
        #[arg(long, default_value_t = 25.0)]
        nmol: f64,
        /// Desired master stock concentration (`uM`, default `100.0`).
        #[arg(long, default_value_t = 100.0)]
        stock_um: f64,
    },

    /// Generate a PCR master-mix pipetting table and thermocycler schedule (prints & optionally writes `.md`).
    Mix {
        #[arg(short, long, value_enum, default_value_t = CliPolymerase::Q5)]
        polymerase: CliPolymerase,
        /// Number of PCR tubes/wells.
        #[arg(short, long, default_value_t = 12)]
        rxns: usize,
        /// Single reaction volume (`uL`).
        #[arg(short, long, default_value_t = 25.0)]
        vol: f64,
        /// Pipetting overage fraction (e.g., `0.10` = +10%).
        #[arg(long, default_value_t = 0.10)]
        overage: f64,
        /// Template DNA volume added per tube (`uL`).
        #[arg(long, default_value_t = 1.0)]
        template_ul: f64,
        /// Annealing temperature `Ta` (`C`).
        #[arg(long, default_value_t = 62.0)]
        ta: f64,
        /// Expected amplicon length (`bp`) for extension time scaling.
        #[arg(long, default_value_t = 500)]
        bp: usize,
        /// Include High-GC Enhancer (Q5) or 3% DMSO.
        #[arg(long)]
        gc_enhancer: bool,
        /// Optional path to write Markdown bench checklist.
        #[arg(short, long)]
        out: Option<PathBuf>,
    },

    /// Wet-lab chemistry & instrument calculators (`buffer`, `molarity`, `dilute`, `serial`, `rcf`).
    Calc {
        #[command(subcommand)]
        calc_cmd: BenchCalcCommand,
    },
}

#[derive(Debug, Subcommand)]
pub enum BenchCalcCommand {
    /// Temperature-corrected Henderson–Hasselbalch dry-salt buffer recipe.
    Buffer {
        #[arg(short, long, value_enum, default_value_t = CliBufferSystem::Tris)]
        system: CliBufferSystem,
        #[arg(long, default_value_t = 8.0)]
        ph: f64,
        /// Total buffer molarity in `M` (e.g., `0.05` for 50 mM).
        #[arg(short, long, default_value_t = 0.05)]
        molarity: f64,
        /// Final buffer volume in liters (`L`).
        #[arg(short, long, default_value_t = 1.0)]
        liters: f64,
        /// Working temperature in Celsius (adjusts effective `pKa`).
        #[arg(short, long, default_value_t = 25.0)]
        temp_c: f64,
    },
    /// Mass (`g`) required to prepare a solution of given molarity (`M`), volume (`L`), and formula weight (`g/mol`).
    Molarity {
        #[arg(short, long)]
        molarity: f64,
        #[arg(short, long)]
        liters: f64,
        #[arg(long)]
        mw: f64,
    },
    /// Single-step `C1 * V1 = C2 * V2` stock dilution solver.
    Dilute {
        #[arg(long)]
        c1: f64,
        #[arg(long)]
        c2: f64,
        #[arg(long)]
        v2: f64,
        #[arg(long, default_value = "uL")]
        unit: String,
    },
    /// Multi-tube serial dilution cascade.
    Serial {
        #[arg(long)]
        start_conc: f64,
        /// Dilution fold per step (e.g., `10.0` for 1:10 serial dilutions).
        #[arg(long, default_value_t = 10.0)]
        fold: f64,
        #[arg(long, default_value_t = 6)]
        steps: usize,
        /// Final usable volume (`uL`) remaining in each tube.
        #[arg(long, default_value_t = 90.0)]
        vol_ul: f64,
    },
    /// Centrifuge rotor `RPM <-> RCF (x g)` converter (pass either `--rpm` or `--rcf`).
    Rcf {
        /// Rotor radius in centimeters (`cm`).
        #[arg(long)]
        radius_cm: f64,
        #[arg(long)]
        rpm: Option<f64>,
        #[arg(long)]
        rcf: Option<f64>,
    },
}

fn parse_span(raw: &str) -> Result<(usize, usize), String> {
    let parts: Vec<&str> = raw.split("..").collect();
    if parts.len() != 2 {
        return Err(format!("Expected range formatted as START..END, got '{raw}'."));
    }
    let a = parts[0]
        .trim()
        .parse::<usize>()
        .map_err(|_| format!("Invalid integer '{}' in range '{raw}'.", parts[0]))?;
    let b = parts[1]
        .trim()
        .trim_start_matches('=')
        .parse::<usize>()
        .map_err(|_| format!("Invalid integer '{}' in range '{raw}'.", parts[1]))?;
    if a > b {
        return Err(format!("Range start ({a}) cannot exceed end ({b})."));
    }
    Ok((a, b))
}

fn read_input_or_inline(input: &str) -> Result<String, String> {
    let p = PathBuf::from(input);
    if p.exists() {
        fs::read_to_string(&p).map_err(|e| format!("Failed reading '{}': {e}", p.display()))
    } else {
        Ok(input.to_string())
    }
}

/// Executes `osabt bench` CLI commands or dispatches to the interactive `Ratatui` TUI.
pub fn run_bench_cli(args: BenchArgs) -> Result<(), String> {
    let Some(cmd) = args.command else {
        #[cfg(feature = "tui")]
        {
            return crate::bench_tui::run_bench_tui().map_err(|e| e.to_string());
        }
        #[cfg(not(feature = "tui"))]
        {
            return Err(
                "TUI feature is disabled in this build; pass a CLI subcommand (e.g., `osabt bench --help`)."
                    .to_string(),
            );
        }
    };

    match cmd {
        #[cfg(feature = "tui")]
        BenchCommand::Tui => crate::bench_tui::run_bench_tui().map_err(|e| e.to_string()),

        BenchCommand::Design {
            input,
            amplicon,
            region,
            tm,
            polymerase,
            qpcr_probe,
            offtarget,
            num,
            export_tsv,
        } => {
            let text = read_input_or_inline(&input)?;
            let records = parse_fasta_or_raw(&text).map_err(|e| e.to_string())?;
            let (rec_name, rec_seq) = &records[0];
            let amp_range = parse_span(&amplicon)?;
            let target_region = match region {
                Some(r) => Some(parse_span(&r)?),
                None => None,
            };
            let bg_records = match offtarget {
                Some(p) => {
                    let bg_txt = fs::read_to_string(&p)
                        .map_err(|e| format!("Failed reading off-target FASTA: {e}"))?;
                    Some(parse_fasta_or_raw(&bg_txt).map_err(|e| e.to_string())?)
                }
                None => None,
            };

            let cfg = PrimerDesignConfig {
                amplicon_len_range: amp_range,
                target_tm_c: tm,
                target_region,
                include_qpcr_probe: qpcr_probe,
                polymerase: polymerase.into(),
                num_return: num,
                ..PrimerDesignConfig::default()
            };

            let cands =
                design_pcr_primers(rec_seq, &cfg, bg_records.as_deref()).map_err(|e| e.to_string())?;
            if cands.is_empty() {
                return Err(
                    "No primer pairs satisfied the requested Tm, amplicon, and secondary-structure gates."
                        .to_string(),
                );
            }

            println!(
                "OSABT Bench — Primer Design Report for '{rec_name}' ({} bp, {:?})\n",
                rec_seq.len(),
                cfg.polymerase
            );
            let mut tsv =
                String::from("Name\tSequence\tLength_nt\tTm_C\tGC_Pct\tAmplicon_bp\tRecommended_Ta_C\n");

            for c in &cands {
                println!(
                    "Rank #{} | Amplicon: {} bp ({}..={}) | Ta: {:.1} C | dTm: {:.2} C | 3'-CrossDimer: {:.2} kcal/mol | Off-Targets: {}",
                    c.rank,
                    c.amplicon_len_bp,
                    c.fwd_start_0 + 1,
                    c.rev_end_0 + 1,
                    c.recommended_ta_c,
                    c.delta_tm_c,
                    c.cross_dimer_3p_dg_kcal,
                    c.off_target_count
                );
                println!(
                    "  FWD ({:>2} nt, Tm {:>5.2} C, GC {:>4.1}%, 3'dG {:>6.2}): 5'-{}-3'",
                    c.fwd.length,
                    c.fwd.tm_c,
                    c.fwd.gc_pct,
                    c.fwd.three_prime_dg_37_kcal,
                    c.fwd.sequence
                );
                println!(
                    "  REV ({:>2} nt, Tm {:>5.2} C, GC {:>4.1}%, 3'dG {:>6.2}): 5'-{}-3'",
                    c.rev.length,
                    c.rev.tm_c,
                    c.rev.gc_pct,
                    c.rev.three_prime_dg_37_kcal,
                    c.rev.sequence
                );
                if let Some(prb) = &c.qpcr_probe {
                    println!(
                        "  PRB ({:>2} nt, Tm {:>5.2} C, GC {:>4.1}%, {:?}): 5'-{}-3'",
                        prb.sequence.len(),
                        prb.tm_c,
                        prb.gc_pct,
                        prb.strand,
                        prb.sequence
                    );
                }
                println!();

                tsv.push_str(&format!(
                    "{rec_name}_P{}_F\t{}\t{}\t{:.2}\t{:.1}\t{}\t{:.1}\n",
                    c.rank,
                    c.fwd.sequence,
                    c.fwd.length,
                    c.fwd.tm_c,
                    c.fwd.gc_pct,
                    c.amplicon_len_bp,
                    c.recommended_ta_c
                ));
                tsv.push_str(&format!(
                    "{rec_name}_P{}_R\t{}\t{}\t{:.2}\t{:.1}\t{}\t{:.1}\n",
                    c.rank,
                    c.rev.sequence,
                    c.rev.length,
                    c.rev.tm_c,
                    c.rev.gc_pct,
                    c.amplicon_len_bp,
                    c.recommended_ta_c
                ));
            }

            if let Some(path) = export_tsv {
                fs::write(&path, tsv)
                    .map_err(|e| format!("Failed writing '{}': {e}", path.display()))?;
                println!("Wrote vendor order sheet -> {}", path.display());
            }
            Ok(())
        }

        BenchCommand::Check {
            fwd,
            rev,
            polymerase,
            offtarget,
        } => {
            let preset: PolymerasePreset = polymerase.into();
            let cond = ReactionConditions::for_preset(preset);
            let norm_f = normalize_dna(&fwd).map_err(|e| e.to_string())?;

            if !is_unambiguous_dna(&norm_f) {
                let (low, high, n_var) =
                    calc_degenerate_tm_bounds(&norm_f, &cond, 256).map_err(|e| e.to_string())?;
                println!(
                    "FWD Degenerate Oligo ({n_var} variants): Tm range {:.2} C .. {:.2} C | GC {:.1}%",
                    low.tm_salt_corrected_c, high.tm_salt_corrected_c, low.gc_pct
                );
                return Ok(());
            }

            let qc_f = evaluate_oligo_qc(&norm_f, &cond, 60.0).map_err(|e| e.to_string())?;
            println!(
                "FWD: 5'-{}-3' ({} nt)\n  Tm: {:.2} C ({:?}, free Mg2+ {:.2} mM) | GC: {:.1}% | 3' Pentamer dG37: {:.2} kcal/mol | G4Hunter: {:.2}",
                qc_f.sequence,
                qc_f.length,
                qc_f.tm_c,
                preset,
                cond.free_mg_molar() * 1000.0,
                qc_f.gc_pct,
                qc_f.three_prime_dg_37_kcal,
                qc_f.g4_hunter_max
            );
            if let Some(hp) = &qc_f.hairpin {
                println!(
                    "  Hairpin: dG = {:.2} kcal/mol, Tm = {:.1} C\n{}",
                    hp.delta_g_kcal, hp.tm_c, hp.ascii_diagram
                );
            }
            if let Some(d3) = &qc_f.self_dimer_3p {
                println!(
                    "  3'-Extensible Self-Dimer: dG37 = {:.2} kcal/mol\n{}",
                    d3.delta_g_kcal, d3.ascii_alignment
                );
            }
            for w in &qc_f.warnings {
                println!("  [WARN] {w}");
            }

            if let Some(raw_r) = rev {
                let norm_r = normalize_dna(&raw_r).map_err(|e| e.to_string())?;
                let qc_r = evaluate_oligo_qc(&norm_r, &cond, 60.0).map_err(|e| e.to_string())?;
                let ta =
                    recommend_annealing_temp(preset, &norm_f, &norm_r).map_err(|e| e.to_string())?;
                println!(
                    "\nREV: 5'-{}-3' ({} nt)\n  Tm: {:.2} C | GC: {:.1}% | 3' Pentamer dG37: {:.2} kcal/mol | G4Hunter: {:.2}",
                    qc_r.sequence,
                    qc_r.length,
                    qc_r.tm_c,
                    qc_r.gc_pct,
                    qc_r.three_prime_dg_37_kcal,
                    qc_r.g4_hunter_max
                );
                println!(
                    "\nPair Summary: dTm = {:.2} C | Recommended {preset:?} Ta = {:.1} C",
                    ta.delta_tm_c, ta.recommended_ta_c
                );
                let (any_cross, d3_cross) =
                    score_dimer(&norm_f, &norm_r, &cond, 37.0).map_err(|e| e.to_string())?;
                if let Some(d) = d3_cross.or(any_cross) {
                    println!(
                        "Worst Cross-Dimer (3'-extensible = {}): dG37 = {:.2} kcal/mol\n{}",
                        d.extensible_3p, d.delta_g_kcal, d.ascii_alignment
                    );
                }
                if let Some(bg_path) = offtarget {
                    let bg_txt = fs::read_to_string(&bg_path)
                        .map_err(|e| format!("Failed reading '{}': {e}", bg_path.display()))?;
                    let recs = parse_fasta_or_raw(&bg_txt).map_err(|e| e.to_string())?;
                    let hits = scan_offtarget_amplicons(&norm_f, &norm_r, &recs, 4000, 1, &cond)
                        .map_err(|e| e.to_string())?;
                    println!("Predicted Amplicons in '{}': {}", bg_path.display(), hits.len());
                    for h in hits {
                        println!(
                            "  - {}: {}..={} ({} bp, mm F={} R={})",
                            h.record_name,
                            h.start_0 + 1,
                            h.end_0 + 1,
                            h.product_len_bp,
                            h.fwd_mismatches,
                            h.rev_mismatches
                        );
                    }
                }
            }
            Ok(())
        }

        BenchCommand::Multiplex {
            input,
            amplicon,
            dg_floor,
            size_spacing,
            polymerase,
            export_tsv,
        } => {
            let txt = fs::read_to_string(&input)
                .map_err(|e| format!("Failed reading '{}': {e}", input.display()))?;
            let recs = parse_fasta_or_raw(&txt).map_err(|e| e.to_string())?;
            let targets: Vec<MultiplexTarget> = recs
                .into_iter()
                .map(|(locus_name, template_sequence)| MultiplexTarget {
                    locus_name,
                    template_sequence,
                    target_region: None,
                })
                .collect();

            let cfg = PrimerDesignConfig {
                amplicon_len_range: parse_span(&amplicon)?,
                polymerase: polymerase.into(),
                ..PrimerDesignConfig::default()
            };
            let pool = design_multiplex_pool(&targets, &cfg, dg_floor, size_spacing)
                .map_err(|e| e.to_string())?;

            println!(
                "OSABT Multiplex Pool ({} loci) | Pooled Ta: {:.1} C | Worst 3'-Dimer: {:.2} kcal/mol\n",
                pool.selected_pairs.len(),
                pool.pooled_ta_c,
                pool.worst_pool_dimer_dg_kcal
            );
            for (locus, pair) in &pool.selected_pairs {
                println!(
                    "  {locus:<18} | {:>4} bp | F: {} ({:.1} C) | R: {} ({:.1} C)",
                    pair.amplicon_len_bp,
                    pair.fwd.sequence,
                    pair.fwd.tm_c,
                    pair.rev.sequence,
                    pair.rev.tm_c
                );
            }
            if let Some(p) = export_tsv {
                fs::write(&p, pool.to_order_tsv())
                    .map_err(|e| format!("Failed writing '{}': {e}", p.display()))?;
                println!("\nWrote multiplex order TSV -> {}", p.display());
            }
            Ok(())
        }

        BenchCommand::Gibson {
            input,
            linear,
            core_tm,
            polymerase,
            export_tsv,
        } => {
            let txt = fs::read_to_string(&input)
                .map_err(|e| format!("Failed reading '{}': {e}", input.display()))?;
            let recs = parse_fasta_or_raw(&txt).map_err(|e| e.to_string())?;
            let frags: Vec<GibsonFragmentInput> = recs
                .into_iter()
                .map(|(name, sequence)| GibsonFragmentInput { name, sequence })
                .collect();
            let cfg = GibsonConfig {
                target_core_tm_c: core_tm,
                polymerase: polymerase.into(),
                ..GibsonConfig::default()
            };
            let plan =
                design_gibson_assembly(&frags, !linear, &cfg).map_err(|e| e.to_string())?;

            println!(
                "OSABT Gibson Assembly ({}, {} bp total, {} junctions)\n",
                if plan.circular { "Circular" } else { "Linear" },
                plan.assembled_length_bp,
                plan.junctions.len()
            );
            for j in &plan.junctions {
                println!(
                    "Junction #{}: {} -> {} | Overlap: {} bp (Tm {:.1} C, dG50 {:.1} kcal/mol)",
                    j.junction_index + 1,
                    j.upstream_fragment,
                    j.downstream_fragment,
                    j.overlap_len,
                    j.overlap_tm_c,
                    j.overlap_dg_50c_kcal
                );
                println!(
                    "  {}: 5'-{}-3' (core Tm {:.1} C, full Tm {:.1} C)",
                    j.downstream_fwd_primer.name,
                    j.downstream_fwd_primer.annotated_sequence,
                    j.downstream_fwd_primer.core_tm_c,
                    j.downstream_fwd_primer.full_tm_c
                );
                println!(
                    "  {}: 5'-{}-3' (core Tm {:.1} C, full Tm {:.1} C)",
                    j.upstream_rev_primer.name,
                    j.upstream_rev_primer.annotated_sequence,
                    j.upstream_rev_primer.core_tm_c,
                    j.upstream_rev_primer.full_tm_c
                );
            }
            println!("\nFragment PCR Setup:");
            for pcr in &plan.fragment_pcrs {
                println!(
                    "  - {}: {} + {} -> {} bp | Ta = {:.1} C | Extend {} s",
                    pcr.fragment_name,
                    pcr.fwd_primer_name,
                    pcr.rev_primer_name,
                    pcr.amplicon_len_bp,
                    pcr.recommended_ta_c,
                    pcr.extension_sec
                );
            }
            if let Some(p) = export_tsv {
                fs::write(&p, plan.to_order_tsv())
                    .map_err(|e| format!("Failed writing '{}': {e}", p.display()))?;
                println!("\nWrote chimeric Gibson oligo sheet -> {}", p.display());
            }
            Ok(())
        }

        BenchCommand::Resuspend {
            csv,
            name,
            seq,
            nmol,
            stock_um,
        } => {
            let plans = if let Some(path) = csv {
                let txt = fs::read_to_string(&path)
                    .map_err(|e| format!("Failed reading '{}': {e}", path.display()))?;
                parse_vendor_coa_table(&txt, stock_um).map_err(|e| e.to_string())?
            } else {
                let s = seq.ok_or_else(|| {
                    "Provide either `--seq <DNA>` or `--csv <vendor_coa.csv>`.".to_string()
                })?;
                vec![calc_oligo_resuspension(&name, &s, nmol, stock_um, 100.0)
                    .map_err(|e| e.to_string())?]
            };

            println!(
                "{:<18} {:>6} {:>9} {:>8} {:>14} {:>22}",
                "Oligo", "Len", "MW(g/mol)", "nmol", "Add TE (uL)", "10 uM Working (100 uL)"
            );
            for p in plans {
                println!(
                    "{:<18} {:>6} {:>9.1} {:>8.1} {:>11.1} uL {:>6.1} uL + {:>5.1} uL H2O",
                    p.name,
                    p.length_nt,
                    p.molecular_weight_g_mol,
                    p.amount_nmol,
                    p.add_buffer_ul,
                    p.working_10um_stock_aliquot_ul,
                    p.working_10um_water_ul
                );
            }
            Ok(())
        }

        BenchCommand::Mix {
            polymerase,
            rxns,
            vol,
            overage,
            template_ul,
            ta,
            bp,
            gc_enhancer,
            out,
        } => {
            let plan = build_pcr_master_mix(
                polymerase.into(),
                rxns,
                vol,
                overage,
                template_ul,
                ta,
                bp,
                gc_enhancer,
            )
            .map_err(|e| e.to_string())?;
            let md = plan.to_markdown_sheet("Standard Assay");
            println!("{md}");
            if let Some(p) = out {
                fs::write(&p, md)
                    .map_err(|e| format!("Failed writing '{}': {e}", p.display()))?;
                println!("Saved Markdown bench checklist -> {}", p.display());
            }
            Ok(())
        }

        BenchCommand::Calc { calc_cmd } => match calc_cmd {
            BenchCalcCommand::Buffer {
                system,
                ph,
                molarity,
                liters,
                temp_c,
            } => {
                let r = calc_buffer_recipe(system.into(), ph, molarity, liters, temp_c)
                    .map_err(|e| e.to_string())?;
                println!(
                    "{:?} Buffer ({:.3} M, {:.2} L) | Target pH {:.2} at {:.1} C (effective pKa = {:.3})\n  - {}: {:.3} g\n  - {}: {:.3} g\n  - Dissolve in ~{:.2} L ddH2O and bring to {:.2} L final volume.",
                    r.system,
                    r.total_molarity_m,
                    r.volume_liters,
                    r.target_ph,
                    r.working_temp_c,
                    r.effective_pka,
                    r.acid_component_name,
                    r.acid_grams,
                    r.base_component_name,
                    r.base_grams,
                    r.volume_liters * 0.8,
                    r.volume_liters
                );
                Ok(())
            }
            BenchCalcCommand::Molarity {
                molarity,
                liters,
                mw,
            } => {
                let g = calc_molarity_grams(molarity, liters, mw).map_err(|e| e.to_string())?;
                println!(
                    "Weigh {:.4} g ({:.2} mg) for {:.4} M in {:.3} L (MW = {:.2} g/mol).",
                    g,
                    g * 1000.0,
                    molarity,
                    liters,
                    mw
                );
                Ok(())
            }
            BenchCalcCommand::Dilute { c1, c2, v2, unit } => {
                let (v1, dil) = calc_c1v1_dilution(c1, c2, v2).map_err(|e| e.to_string())?;
                println!(
                    "Combine {:.3} {unit} stock ({c1}) + {:.3} {unit} diluent -> {:.3} {unit} at {c2}.",
                    v1, dil, v2
                );
                Ok(())
            }
            BenchCalcCommand::Serial {
                start_conc,
                fold,
                steps,
                vol_ul,
            } => {
                let cascade = calc_serial_dilution(start_conc, fold, steps, vol_ul)
                    .map_err(|e| e.to_string())?;
                println!(
                    "Serial 1:{fold:.1} Cascade (Pre-fill each tube with {:.2} uL diluent):",
                    cascade[0].diluent_ul
                );
                for s in cascade {
                    println!(
                        "  Tube #{}: transfer {:.2} uL in -> {:.4e} (final vol {:.1} uL)",
                        s.tube_index, s.transfer_in_ul, s.concentration, s.remaining_final_ul
                    );
                }
                Ok(())
            }
            BenchCalcCommand::Rcf {
                radius_cm,
                rpm,
                rcf,
            } => {
                if let Some(r) = rpm {
                    let g = rpm_to_rcf(r, radius_cm).map_err(|e| e.to_string())?;
                    println!("{r:.0} RPM at r = {radius_cm:.2} cm -> {g:.1} x g");
                } else if let Some(g) = rcf {
                    let r = rcf_to_rpm(g, radius_cm).map_err(|e| e.to_string())?;
                    println!("{g:.1} x g at r = {radius_cm:.2} cm -> {r:.0} RPM");
                } else {
                    return Err("Provide either `--rpm <RPM>` or `--rcf <XG>`.".to_string());
                }
                Ok(())
            }
        },
    }
}
