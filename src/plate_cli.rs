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

//! Command-line subcommands (`clap` v4) for Tool #9: `osabt plate`.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::plate_curve::{
    demo_measured_96well_plate, evaluate_plate_qc_and_normalize, fit_all_plate_compounds,
    fit_four_param_logistic, render_4pl_curves_svg, FourParamLogisticFit,
};
use crate::plate_layout::{
    export_bench_pipetting_markdown, export_plate_pipetting_csv, generate_randomized_plate,
    render_plate_map_svg, CompoundSeriesSpec, PlateFormat, PlatePlanConfig,
};

#[derive(Debug, Args)]
pub struct PlateArgs {
    #[command(subcommand)]
    pub command: Option<PlateSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum PlateSubcommand {
    /// Launch the interactive 4-tab Ratatui 96/384-well assay planner & 4PL curve fitter.
    Tui,

    /// Generate an edge-buffered, Fisher-Yates randomized 96- or 384-well layout, bench protocol, CSV, & SVG map.
    Plan {
        /// Microplate format: `96` or `384`.
        #[arg(long, default_value_t = 96)]
        wells: usize,
        /// Number of outer PBS/media evaporation moat rings (`0..2`).
        #[arg(long, default_value_t = 1)]
        moat: usize,
        /// Deterministic PRNG seed for reproducible inner-well randomization.
        #[arg(long, default_value_t = 2026_1111)]
        seed: u64,
        /// Comma-separated compound names to plate in serial dilution.
        #[arg(long, default_value = "Staurosporine,CRV-104,CRV-209")]
        compounds: String,
        /// Top dose concentration in nM.
        #[arg(long, default_value_t = 2000.0)]
        top_nm: f64,
        /// Serial dilution fold step (`> 1.0`).
        #[arg(long, default_value_t = 4.0)]
        fold: f64,
        /// Number of dilution doses per compound.
        #[arg(long, default_value_t = 8)]
        doses: usize,
        /// Replicate wells per dose.
        #[arg(long, default_value_t = 2)]
        replicates: usize,
        /// Output prefix for `<prefix>_layout.csv`, `<prefix>_protocol.md`, and `<prefix>_map.svg`.
        #[arg(short, long, default_value = "osabt_plate")]
        out_prefix: String,
    },

    /// Fit 4-Parameter Logistic (`4PL`) `IC50`/`EC50` curves from a `compound,conc_nM,response` CSV (or built-in demo).
    Fit {
        /// Optional CSV/TSV file with columns `compound,concentration_nM,response` (uses built-in kinase assay if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Optional `.svg` semi-log 4PL curve figure path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Compute Zhang et al. (1999) `Z'`-factor, `SSMD`, and `%CV` on the built-in or measured assay plate.
    Qc {
        /// Optional `.svg` viability heatmap output path.
        #[arg(long)]
        heatmap_svg: Option<PathBuf>,
    },
}

pub fn run_plate_cli(args: PlateArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "tui")]
            {
                return crate::plate_tui::run_plate_tui();
            }
            #[cfg(not(feature = "tui"))]
            {
                return Err(String::from(
                    "Build with `--features tui` or run `osabt plate plan` / `osabt plate fit`.",
                ));
            }
        }
        Some(PlateSubcommand::Tui) => {
            #[cfg(feature = "tui")]
            {
                crate::plate_tui::run_plate_tui()
            }
            #[cfg(not(feature = "tui"))]
            {
                Err(String::from("TUI feature is disabled in this build."))
            }
        }
        Some(PlateSubcommand::Plan {
            wells,
            moat,
            seed,
            compounds,
            top_nm,
            fold,
            doses,
            replicates,
            out_prefix,
        }) => {
            let format = match wells {
                96 => PlateFormat::Well96,
                384 => PlateFormat::Well384,
                other => return Err(format!("Unsupported plate density '{other}'; choose 96 or 384.")),
            };
            let specs: Vec<CompoundSeriesSpec> = compounds
                .split(',')
                .map(str::trim)
                .filter(|s| !s.is_empty())
                .map(|name| CompoundSeriesSpec {
                    name: name.to_string(),
                    top_conc_nm: top_nm,
                    dilution_fold: fold,
                    num_doses: doses,
                    replicates,
                })
                .collect();

            let cfg = PlatePlanConfig {
                format,
                plate_name: out_prefix.clone(),
                moat_rings: moat,
                negative_control_wells: 6,
                positive_control_wells: 6,
                assay_volume_ul: if wells == 384 { 40.0 } else { 100.0 },
                rng_seed: seed,
                compounds: specs,
            };

            let plate = generate_randomized_plate(&cfg)?;
            let csv_p = format!("{out_prefix}_layout.csv");
            let md_p = format!("{out_prefix}_protocol.md");
            let svg_p = format!("{out_prefix}_map.svg");

            fs::write(&csv_p, export_plate_pipetting_csv(&plate))
                .map_err(|e| format!("Failed writing '{csv_p}': {e}"))?;
            fs::write(&md_p, export_bench_pipetting_markdown(&plate, &cfg))
                .map_err(|e| format!("Failed writing '{md_p}': {e}"))?;
            fs::write(&svg_p, render_plate_map_svg(&plate))
                .map_err(|e| format!("Failed writing '{svg_p}': {e}"))?;

            println!(
                "Generated randomized {} assay plan -> {csv_p}, {md_p}, {svg_p}",
                format.label()
            );
            Ok(())
        }
        Some(PlateSubcommand::Fit { input, svg_out }) => {
            let fits: Vec<FourParamLogisticFit> = match input {
                Some(p) => {
                    let raw = fs::read_to_string(&p)
                        .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                    parse_and_fit_dose_csv(&raw)?
                }
                None => {
                    let plate = demo_measured_96well_plate();
                    fit_all_plate_compounds(&plate)
                }
            };

            println!(
                "compound\tn_pts\trel_IC50_nM\tci95_low_nM\tci95_high_nM\tabs_IC50_nM\tIC90_nM\thill_slope\tbottom\ttop\tr_squared"
            );
            for f in &fits {
                let abs_s = f
                    .absolute_ic50_nm
                    .map(|v| format!("{v:.3}"))
                    .unwrap_or_else(|| String::from("NA"));
                println!(
                    "{}\t{}\t{:.3}\t{:.3}\t{:.3}\t{}\t{:.2}\t{:+.3}\t{:.2}\t{:.2}\t{:.4}",
                    f.compound_name,
                    f.n_points,
                    f.relative_ec50_nm,
                    f.relative_ec50_ci95_nm.0,
                    f.relative_ec50_ci95_nm.1,
                    abs_s,
                    f.ec90_nm,
                    f.hill_slope,
                    f.bottom,
                    f.top,
                    f.r_squared
                );
            }

            if let Some(sp) = svg_out {
                fs::write(&sp, render_4pl_curves_svg("Dose-Response 4PL Analysis", &fits))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                eprintln!("Saved 4PL dose-response SVG -> {}", sp.display());
            }
            Ok(())
        }
        Some(PlateSubcommand::Qc { heatmap_svg }) => {
            let mut plate = demo_measured_96well_plate();
            let qc = evaluate_plate_qc_and_normalize(&mut plate)?;
            println!("=== Zhang et al. (1999) Microplate Assay QC ===");
            println!("Quality Grade          : {}", qc.quality_grade);
            println!("Z'-Factor              : {:.4}", qc.z_prime_factor);
            println!("SSMD                   : {:.3}", qc.ssmd);
            println!("Signal / Background    : {:.2}-fold", qc.signal_to_background);
            println!(
                "Negative Ctrl (n={})    : {:.1} +- {:.1} (CV {:.2}%)",
                qc.neg_ctrl_n, qc.neg_mean, qc.neg_sd, qc.neg_cv_pct
            );
            println!(
                "Positive Ctrl (n={})    : {:.1} +- {:.1} (CV {:.2}%)",
                qc.pos_ctrl_n, qc.pos_mean, qc.pos_sd, qc.pos_cv_pct
            );
            if let Some(sp) = heatmap_svg {
                fs::write(&sp, render_plate_map_svg(&plate))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                println!("Saved heatmap SVG      -> {}", sp.display());
            }
            Ok(())
        }
    }
}

fn parse_and_fit_dose_csv(raw: &str) -> Result<Vec<FourParamLogisticFit>, String> {
    let mut by_comp: BTreeMap<String, Vec<(f64, f64)>> = BTreeMap::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let cols: Vec<&str> = t.split([',', '\t', ';']).map(str::trim).collect();
        if cols.len() < 3 {
            continue;
        }
        let Ok(conc) = cols[1].parse::<f64>() else {
            continue; // skip header row
        };
        let Ok(resp) = cols[2].parse::<f64>() else {
            continue;
        };
        if conc > 0.0 {
            by_comp
                .entry(cols[0].to_string())
                .or_default()
                .push((conc, resp));
        }
    }
    if by_comp.is_empty() {
        return Err(String::from(
            "No valid `compound,conc_nM,response` rows parsed from input file.",
        ));
    }
    let mut out = Vec::new();
    for (name, obs) in by_comp {
        out.push(fit_four_param_logistic(&name, &obs)?);
    }
    Ok(out)
}
