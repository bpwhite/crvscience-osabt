// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #13 (`osabt kinetics`):
//! - `osabt kinetics gui` — Launches the `eframe`/`egui` Enzyme Kinetics & Morrison Studio
//! - `osabt kinetics fit` — Fits Michaelis–Menten, Allosteric Hill, & Haldane Substrate-Inhibition models (`AICc`)
//! - `osabt kinetics inhibit` — Globally fits all 4 reversible inhibition mechanisms (`AICc` + Akaike weights + `.svg`)
//! - `osabt kinetics morrison` — Fits Morrison (1969) tight-binding quadratic `[E]_T` & `Ki` + `.svg`
//! - `osabt kinetics cheng-prusoff` — Converts `IC50` to `Ki` via Cheng & Prusoff (1973)

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

use crate::kinetics_core::{
    convert_ic50_cheng_prusoff, demo_competitive_kinase_matrix_csv, demo_morrison_tight_binding_csv,
    fit_global_reversible_inhibition, fit_morrison_tight_binding, fit_single_substrate_models,
    parse_kinetics_csv, InhibitionMechanism, KineticsPoint,
};
use crate::kinetics_svg::{render_global_inhibition_svg, render_morrison_svg};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliMechanism {
    Competitive,
    Uncompetitive,
    Noncompetitive,
}

impl From<CliMechanism> for InhibitionMechanism {
    fn from(m: CliMechanism) -> Self {
        match m {
            CliMechanism::Competitive => Self::Competitive,
            CliMechanism::Uncompetitive => Self::Uncompetitive,
            CliMechanism::Noncompetitive => Self::PureNoncompetitive,
        }
    }
}

/// Top-level CLI arguments for `osabt kinetics`.
#[derive(Debug, Args)]
pub struct KineticsArgs {
    #[command(subcommand)]
    pub command: Option<KineticsSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum KineticsSubcommand {
    /// Launch the interactive `eframe`/`egui` Enzyme Kinetics & Morrison Tight-Binding desktop studio.
    Gui {
        /// Optional CSV/TSV kinetics file to open on startup.
        #[arg(long)]
        input: Option<PathBuf>,
    },

    /// Fit Michaelis–Menten, Allosteric Hill, and Substrate-Inhibition models ranked by `AICc`.
    Fit {
        /// Optional CSV/TSV `[S], v0` file (omit to fit `[I] = 0` slice of built-in kinase dataset).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Optional active enzyme concentration `[E]_T` (`nM`) to report `kcat` (`s^-1`) and `kcat/Km` (`M^-1 s^-1`).
        #[arg(long)]
        enzyme_nm: Option<f64>,
    },

    /// Globally fit Competitive, Uncompetitive, Pure Noncompetitive, and Mixed-Type inhibition ranked by `AICc`.
    Inhibit {
        /// Optional 3-column CSV/TSV (`substrate_um,inhibitor_nm,velocity`; omit for built-in kinase matrix).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Optional path to write 2-panel Michaelis–Menten + Eadie–Hofstee `.svg` plot.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Fit Morrison's (1969) tight-binding quadratic equation to recover active `[E]_T` and `Ki_app` / intrinsic `Ki`.
    Morrison {
        /// Optional 2-column CSV/TSV (`inhibitor_nm,velocity`; omit for built-in tight-binding titration).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Optional fixed `[E]_T` (`nM`); omit to titrate active `[E]_T` directly from the quadratic bend.
        #[arg(long)]
        fixed_et_nm: Option<f64>,

        /// Optional nominal protein concentration (`nM`) to compute `%` active enzyme purity.
        #[arg(long, default_value_t = 25.0)]
        nominal_et_nm: f64,

        /// Optional competitive substrate concentration `[S]` (`µM`) for intrinsic `Ki` conversion.
        #[arg(long)]
        substrate_um: Option<f64>,

        /// Optional substrate `Km` (`µM`) paired with `--substrate-um`.
        #[arg(long)]
        km_um: Option<f64>,

        /// Optional path to write `.svg` Morrison tight-binding titration plot.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Convert classical `IC50` (`nM`) to `Ki` (`nM`) via Cheng & Prusoff (1973).
    ChengPrusoff {
        #[arg(long, value_enum, default_value_t = CliMechanism::Competitive)]
        mode: CliMechanism,

        #[arg(long)]
        ic50_nm: f64,

        #[arg(long)]
        substrate_um: f64,

        #[arg(long)]
        km_um: f64,
    },
}

/// Dispatches `osabt kinetics` CLI subcommands.
pub fn run_kinetics_cli(args: KineticsArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_gui(None),
        Some(KineticsSubcommand::Gui { input }) => launch_default_gui(input),
        Some(KineticsSubcommand::Fit { input, enzyme_nm }) => {
            let pts = load_or_demo_matrix(input.as_ref())?;
            let fits = fit_single_substrate_models(&pts, enzyme_nm)?;
            println!("=== OSABT Kinetics — Single-Substrate Steady-State Model Selection ===");
            for (idx, f) in fits.iter().enumerate() {
                println!(
                    "  #{} {:28} : Vmax={:.2} | Km/K0.5={:.2} µM | Sec={:?} | kcat={:?} s^-1 | R²={:.5} | AICc={:.2} (w={:.1}%)",
                    idx + 1,
                    f.model.label(),
                    f.vmax,
                    f.km_or_khalf_um,
                    f.secondary_param,
                    f.kcat_per_s,
                    f.r_squared,
                    f.aicc,
                    f.akaike_weight * 100.0
                );
            }
            Ok(())
        }
        Some(KineticsSubcommand::Inhibit { input, svg_out }) => {
            let pts = load_or_demo_matrix(input.as_ref())?;
            let fits = fit_global_reversible_inhibition(&pts)?;
            println!("=== OSABT Kinetics — Global 2D Reversible Inhibition AICc Ranking ===");
            for (idx, f) in fits.iter().enumerate() {
                println!(
                    "  #{} {:32} : Vmax={:.2} | Km={:.2} µM | Ki={:.2} nM | α={:?} | R²={:.5} | ΔAICc={:+.2} | Weight={:.1}%",
                    idx + 1,
                    f.mechanism.label(),
                    f.vmax,
                    f.km_um,
                    f.ki_nm,
                    f.alpha_factor,
                    f.r_squared,
                    f.delta_aicc,
                    f.akaike_weight * 100.0
                );
            }
            if let Some(out) = svg_out {
                let svg = render_global_inhibition_svg(&pts, &fits[0], "Global Reversible Inhibition");
                fs::write(&out, svg).map_err(|e| {
                    format!("Failed writing global inhibition SVG '{}': {e}", out.display())
                })?;
                println!("Wrote Global Inhibition + Eadie–Hofstee SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(KineticsSubcommand::Morrison {
            input,
            fixed_et_nm,
            nominal_et_nm,
            substrate_um,
            km_um,
            svg_out,
        }) => {
            let pts = if let Some(p) = input {
                let raw = fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                parse_kinetics_csv(&raw)?
            } else {
                parse_kinetics_csv(demo_morrison_tight_binding_csv())?
            };
            let obs: Vec<(f64, f64)> = pts.iter().map(|p| (p.substrate_um, p.velocity)).collect();
            let comp_pair = substrate_um.zip(km_um).or(Some((70.0, 35.0)));
            let fit = fit_morrison_tight_binding(&obs, fixed_et_nm, Some(nominal_et_nm), comp_pair)?;

            println!("=== Morrison (1969) Tight-Binding Quadratic Regression ===");
            println!(
                "v0 = {:.2} | Active [E]_T = {:.2} nM ({:.1}% of nominal) | Ki_app = {:.3} nM | Intrinsic Ki = {:?} nM | R² = {:.5}",
                fit.v0_unhibited,
                fit.active_enzyme_et_nm,
                fit.fractional_active_enzyme.unwrap_or(1.0) * 100.0,
                fit.ki_app_nm,
                fit.intrinsic_ki_nm,
                fit.r_squared
            );
            if let Some(out) = svg_out {
                let svg = render_morrison_svg(&obs, &fit, "Morrison Tight-Binding Titration");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing Morrison SVG '{}': {e}", out.display()))?;
                println!("Wrote Morrison SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(KineticsSubcommand::ChengPrusoff {
            mode,
            ic50_nm,
            substrate_um,
            km_um,
        }) => {
            let res = convert_ic50_cheng_prusoff(mode.into(), ic50_nm, substrate_um, km_um);
            println!(
                "Cheng & Prusoff (1973) [{}] : IC50 = {:.2} nM ([S] = {:.1} µM, Km = {:.1} µM) => Ki = {:.3} nM",
                res.mechanism.label(),
                res.ic50_nm,
                res.substrate_um,
                res.km_um,
                res.computed_ki_nm
            );
            Ok(())
        }
    }
}

fn load_or_demo_matrix(path: Option<&PathBuf>) -> Result<Vec<KineticsPoint>, String> {
    if let Some(p) = path {
        let raw = fs::read_to_string(p)
            .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
        parse_kinetics_csv(&raw)
    } else {
        parse_kinetics_csv(demo_competitive_kinase_matrix_csv())
    }
}

#[cfg(feature = "gui")]
fn launch_default_gui(input: Option<PathBuf>) -> Result<(), String> {
    crate::kinetics_gui::run_kinetics_gui(input)
}

#[cfg(not(feature = "gui"))]
fn launch_default_gui(_input: Option<PathBuf>) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `gui` feature; rebuild with `--features gui` or run `osabt kinetics inhibit` / `morrison`.",
    ))
}
