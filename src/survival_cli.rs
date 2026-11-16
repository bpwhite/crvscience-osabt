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

//! Command-line subcommands (`clap` v4) for Tool #19 (`osabt survival`):
//! - `osabt survival tui` — Interactive 4-tab `Ratatui` Kaplan–Meier, Cox PH & `LC50` Studio
//! - `osabt survival km` — Computes Kaplan–Meier curves, median survival, & Mantel (1966) log-rank test
//! - `osabt survival cox` — Fits Cox (1972) proportional-hazards regression (`HR`, 95% CI, Wald `p`)
//! - `osabt survival lc50` — Fits Finney (1971) Probit or Logit acute `LC10/LC50/LC90` & exports `.svg`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

use crate::survival_core::{
    compute_kaplan_meier_curves, compute_log_rank_test, demo_survival_cohort_csv,
    fit_cox_proportional_hazards, parse_survival_csv,
};
use crate::survival_tox::{
    demo_ecotox_lc50_csv, fit_ecotox_lc50, parse_ecotox_csv, render_survival_and_lc50_svg,
    DoseResponseLink,
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum LinkArg {
    Probit,
    Logit,
}

impl From<LinkArg> for DoseResponseLink {
    fn from(v: LinkArg) -> Self {
        match v {
            LinkArg::Probit => Self::Probit,
            LinkArg::Logit => Self::Logit,
        }
    }
}

/// Top-level CLI arguments for `osabt survival`.
#[derive(Debug, Args)]
pub struct SurvivalArgs {
    #[command(subcommand)]
    pub command: Option<SurvivalSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum SurvivalSubcommand {
    /// Launch the interactive 4-tab `Ratatui` Survival & Ecotoxicology `LC50` Studio.
    Tui {
        /// Optional right-censored survival CSV (`id,group,time,event[,covariate]`).
        #[arg(long)]
        survival: Option<PathBuf>,

        /// Optional concentration-mortality bioassay CSV (`concentration,exposed_n,dead_n`).
        #[arg(long)]
        ecotox: Option<PathBuf>,
    },

    /// Compute Kaplan–Meier product-limit survival curves and Mantel (1966) log-rank test.
    Km {
        #[arg(long)]
        input: Option<PathBuf>,
    },

    /// Fit Cox (1972) proportional-hazards partial-likelihood regression (`HR` & 95% CI).
    Cox {
        #[arg(long)]
        input: Option<PathBuf>,

        /// Regress on group indicator instead of numeric covariate column.
        #[arg(long, default_value_t = false)]
        by_group: bool,
    },

    /// Fit Finney Probit or Binomial Logit acute ecotoxicology `LC10/LC50/LC90` and export `.svg`.
    Lc50 {
        #[arg(long)]
        input: Option<PathBuf>,

        #[arg(long, value_enum, default_value_t = LinkArg::Probit)]
        link: LinkArg,

        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

/// Dispatches `osabt survival` CLI subcommands.
pub fn run_survival_cli(args: SurvivalArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_tui(None, None),
        Some(SurvivalSubcommand::Tui { survival, ecotox }) => launch_default_tui(survival, ecotox),
        Some(SurvivalSubcommand::Km { input }) => {
            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading survival CSV '{}': {e}", p.display()))?
            } else {
                demo_survival_cohort_csv().to_string()
            };
            let obs = parse_survival_csv(&raw)?;
            let curves = compute_kaplan_meier_curves(&obs);
            println!("=== OSABT Survival — Kaplan–Meier Product-Limit Curves ===");
            for c in &curves {
                let med = c
                    .median_survival_time
                    .map(|m| format!("{m:.1}"))
                    .unwrap_or_else(|| String::from(">max_t"));
                let fin = c.steps.last().map(|s| s.survival_s).unwrap_or(1.0);
                println!(
                    "  {:18} : n={:2} | events={:2} | median={} | final S(t)={:.3}",
                    c.group, c.total_n, c.total_events, med, fin
                );
            }
            if let Ok(lr) = compute_log_rank_test(&obs) {
                println!(
                    "  Mantel (1966) Log-Rank ({} vs {}): Chi^2 = {:.3}, p = {:.5}",
                    lr.group_a, lr.group_b, lr.chi_square, lr.p_value
                );
            }
            Ok(())
        }
        Some(SurvivalSubcommand::Cox { input, by_group }) => {
            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading survival CSV '{}': {e}", p.display()))?
            } else {
                demo_survival_cohort_csv().to_string()
            };
            let obs = parse_survival_csv(&raw)?;
            let label = if by_group { "group" } else { "covariate" };
            let cox = fit_cox_proportional_hazards(&obs, label, by_group)?;
            println!(
                "=== OSABT Survival — Cox (1972) Proportional Hazards ({}, n={}, events={}) ===",
                cox.covariate_name, cox.n_subjects, cox.n_events
            );
            println!(
                "  beta = {:+.4} ± {:.4} | HR = {:.4} (95% CI [{:.4}, {:.4}]) | Wald z = {:+.3} | p = {:.6}",
                cox.beta,
                cox.se_beta,
                cox.hazard_ratio,
                cox.hr_ci95.0,
                cox.hr_ci95.1,
                cox.wald_z,
                cox.p_value
            );
            Ok(())
        }
        Some(SurvivalSubcommand::Lc50 {
            input,
            link,
            svg_out,
        }) => {
            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading ecotox CSV '{}': {e}", p.display()))?
            } else {
                demo_ecotox_lc50_csv().to_string()
            };
            let pts = parse_ecotox_csv(&raw)?;
            let rep = fit_ecotox_lc50(&pts, DoseResponseLink::from(link))?;
            println!("=== OSABT Survival — {} Acute LC50 ===", rep.link.label());
            println!(
                "  LC10 = {:.3} | LC50 = {:.3} (95% CI [{:.3}, {:.3}]) | LC90 = {:.3} | Pearson Chi^2 = {:.3}",
                rep.lc10, rep.lc50, rep.lc50_ci95.0, rep.lc50_ci95.1, rep.lc90, rep.pearson_chi2
            );
            if let Some(out) = svg_out {
                let obs = parse_survival_csv(demo_survival_cohort_csv())?;
                let curves = compute_kaplan_meier_curves(&obs);
                let lr = compute_log_rank_test(&obs).ok();
                let svg = render_survival_and_lc50_svg(
                    &curves,
                    lr.as_ref(),
                    &pts,
                    &rep,
                    "OSABT Survival & Ecotoxicology",
                );
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", out.display()))?;
                println!("Wrote Survival & LC50 SVG -> {}", out.display());
            }
            Ok(())
        }
    }
}

#[cfg(feature = "tui")]
fn launch_default_tui(
    survival: Option<PathBuf>,
    ecotox: Option<PathBuf>,
) -> Result<(), String> {
    crate::survival_tui::run_survival_tui(survival, ecotox)
}

#[cfg(not(feature = "tui"))]
fn launch_default_tui(
    _survival: Option<PathBuf>,
    _ecotox: Option<PathBuf>,
) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `tui` feature; rebuild with `--features tui` or run `osabt survival km`.",
    ))
}
