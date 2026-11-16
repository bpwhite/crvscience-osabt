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

//! Command-line subcommands (`clap` v4) for Tool #5: `osabt epi`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::epi_models::{
    demo_outbreak_incidence, estimate_cori_rt, evaluate_diagnostic_cascade,
    rogan_gladen_prevalence, simulate_seir_rk4, CoriRtConfig, DiagnosticStepInput, SeirParams,
};
use crate::epi_power::{
    design_case_control_study, design_cohort_study, design_test_negative_study,
    render_power_curve_svg, render_rt_epidemic_svg, render_seir_trajectory_svg,
};

#[derive(Debug, Args)]
pub struct EpiArgs {
    #[command(subcommand)]
    pub command: Option<EpiSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum EpiSubcommand {
    /// Launch the interactive 4-tab Ratatui epidemiology & study-design power studio.
    Tui {
        /// Optional daily incidence CSV/TSV file for the Cori R_t tab.
        #[arg(short, long)]
        input: Option<PathBuf>,
    },

    /// Estimate time-varying reproduction numbers (`R_t`) from daily case counts via Cori et al. (2013).
    Rt {
        /// Optional daily case CSV/TSV (uses built-in 35-day outbreak wave if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Mean serial interval in days.
        #[arg(long, default_value_t = 4.7)]
        si_mean: f64,
        /// Standard deviation of the serial interval in days.
        #[arg(long, default_value_t = 2.9)]
        si_sd: f64,
        /// Sliding estimation window length `tau` in days.
        #[arg(short, long, default_value_t = 7)]
        window: usize,
        /// Optional `.svg` plot output path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Run a 4th-order Runge-Kutta (`RK4`) `SIR` (`--incubation 0`) or `SEIR` compartmental sweep.
    Seir {
        #[arg(long, default_value_t = 100_000.0)]
        pop: f64,
        #[arg(long, default_value_t = 2.8)]
        r0: f64,
        /// Mean latent / incubation duration in days (`0.0` runs classic SIR).
        #[arg(long, default_value_t = 3.5)]
        incubation: f64,
        /// Mean infectious duration in days.
        #[arg(long, default_value_t = 5.0)]
        infectious: f64,
        /// Initial immune/vaccinated proportion (`0.0..1.0`).
        #[arg(long, default_value_t = 0.15)]
        immune: f64,
        /// Simulation duration in days.
        #[arg(long, default_value_t = 160)]
        days: usize,
        /// Optional NPI start day.
        #[arg(long)]
        npi_day: Option<f64>,
        /// Transmission reduction fraction (`0.0..1.0`) starting at `--npi-day`.
        #[arg(long, default_value_t = 0.0)]
        npi_reduction: f64,
        /// Optional `.svg` trajectory plot output path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Sequential Bayesian diagnostic update (`PPV`, `NPV`, `LR+`, `LR-`) and optional Rogan-Gladen prevalence.
    Bayes {
        /// Pre-test probability / baseline prevalence (`0.0..1.0`).
        #[arg(long, default_value_t = 0.08)]
        prior: f64,
        /// Sensitivity of primary assay (`0.0..1.0`).
        #[arg(long, default_value_t = 0.88)]
        se: f64,
        /// Specificity of primary assay (`0.0..1.0`).
        #[arg(long, default_value_t = 0.985)]
        sp: f64,
        /// Whether observed test result is negative (defaults to positive).
        #[arg(long)]
        negative: bool,
        /// Optional field survey positive count for Rogan-Gladen true prevalence estimation.
        #[arg(long)]
        positives: Option<usize>,
        /// Optional field survey total sample size for Rogan-Gladen true prevalence estimation.
        #[arg(long)]
        tested: Option<usize>,
    },

    /// Sample-size & statistical power calculator for Cohort (`RR`), Case-Control (`OR`), and Test-Negative (`VE`) studies.
    Power {
        /// Study design: `tnd` (Test-Negative VE), `cohort` (Relative Risk), or `case-control` (Odds Ratio).
        #[arg(long, default_value = "tnd")]
        design: String,
        /// Control exposure / vaccination coverage (for `tnd`/`case-control`) OR baseline unexposed risk (for `cohort`).
        #[arg(long, default_value_t = 0.50)]
        p0: f64,
        /// Target Vaccine Effectiveness `VE` (`0..1` for `tnd`), Relative Risk `RR` (`cohort`), or Odds Ratio `OR` (`case-control`).
        #[arg(long, default_value_t = 0.60)]
        effect: f64,
        /// Control-to-case or unexposed-to-exposed ratio `k = n0 / n1`.
        #[arg(short = 'k', long, default_value_t = 2.0)]
        ratio: f64,
        /// Two-sided Type I error `alpha`.
        #[arg(long, default_value_t = 0.05)]
        alpha: f64,
        /// Desired statistical power `1 - beta`.
        #[arg(long, default_value_t = 0.80)]
        power: f64,
        /// Optional expected pathogen positivity fraction (`0..1`) among ARI swabs in `tnd` studies.
        #[arg(long)]
        positivity: Option<f64>,
        /// Optional `.svg` power curve output path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

pub fn run_epi_cli(args: EpiArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "tui")]
            {
                return crate::epi_tui::run_epi_tui(None);
            }
            #[cfg(not(feature = "tui"))]
            {
                return Err(String::from(
                    "Build with `--features tui` or run `osabt epi rt` / `osabt epi power`.",
                ));
            }
        }
        Some(EpiSubcommand::Tui { input }) => {
            #[cfg(feature = "tui")]
            {
                crate::epi_tui::run_epi_tui(input)
            }
            #[cfg(not(feature = "tui"))]
            {
                let _ = input;
                Err(String::from("TUI feature is disabled in this build."))
            }
        }
        Some(EpiSubcommand::Rt {
            input,
            si_mean,
            si_sd,
            window,
            svg_out,
        }) => {
            let inc = match input {
                Some(p) => {
                    let raw = fs::read_to_string(&p)
                        .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                    parse_cli_incidence(&raw)?
                }
                None => demo_outbreak_incidence(),
            };
            let cfg = CoriRtConfig {
                serial_interval_mean: si_mean,
                serial_interval_sd: si_sd,
                window_days: window,
                ..CoriRtConfig::default()
            };
            let pts = estimate_cori_rt(&inc, &cfg)?;
            println!(
                "day\tincidence\tlambda_t\trt_mean\trt_median\tci95_low\tci95_high\tprob_rt_gt_1"
            );
            for p in &pts {
                println!(
                    "{}\t{:.1}\t{:.3}\t{}\t{}\t{}\t{}\t{}",
                    p.day_index + 1,
                    p.incidence,
                    p.total_infectiousness_lambda,
                    fmt_opt(p.rt_mean),
                    fmt_opt(p.rt_median),
                    fmt_opt(p.ci95_low),
                    fmt_opt(p.ci95_high),
                    fmt_opt(p.prob_rt_gt_1)
                );
            }
            if let Some(sp) = svg_out {
                fs::write(&sp, render_rt_epidemic_svg("Epidemic Renewal R_t", &pts))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                eprintln!("Saved R_t epidemic SVG -> {}", sp.display());
            }
            Ok(())
        }
        Some(EpiSubcommand::Seir {
            pop,
            r0,
            incubation,
            infectious,
            immune,
            days,
            npi_day,
            npi_reduction,
            svg_out,
        }) => {
            let params = SeirParams {
                population: pop,
                initial_infectious: 10.0,
                initial_exposed: if incubation > 0.0 { 15.0 } else { 0.0 },
                initial_immune_fraction: immune,
                r0,
                incubation_days: incubation,
                infectious_days: infectious,
                total_days: days,
                npi_start_day: npi_day,
                npi_reduction_fraction: npi_reduction,
            };
            let sim = simulate_seir_rk4(&params)?;
            println!("=== RK4 Compartmental Epidemic Summary ===");
            println!(
                "Herd Immunity Threshold : {:.2}%",
                sim.herd_immunity_threshold * 100.0
            );
            println!("Initial Effective Re(0) : {:.3}", sim.initial_effective_re);
            println!(
                "Peak Infectious         : {:.1} on Day {:.0}",
                sim.peak_infectious, sim.peak_day
            );
            println!(
                "Final Attack Fraction   : {:.2}% ({:.0} total infected)",
                sim.final_attack_fraction * 100.0,
                sim.total_infected_persons
            );
            if let Some(sp) = svg_out {
                fs::write(&sp, render_seir_trajectory_svg(&sim))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                println!("Saved SEIR trajectory SVG -> {}", sp.display());
            }
            Ok(())
        }
        Some(EpiSubcommand::Bayes {
            prior,
            se,
            sp,
            negative,
            positives,
            tested,
        }) => {
            let step = DiagnosticStepInput {
                assay_name: String::from("Primary Assay"),
                sensitivity: se,
                specificity: sp,
                observed_positive: !negative,
            };
            let out = &evaluate_diagnostic_cascade(prior, &[step])?[0];
            println!("Pre-Test Probability  : {:.3}%", out.pre_test_prob * 100.0);
            println!(
                "Likelihood Ratios     : LR+ = {:.3} | LR- = {:.4} | Youden J = {:.3}",
                out.lr_plus, out.lr_minus, out.youden_j
            );
            println!(
                "Predictive Values     : PPV = {:.3}% | NPV = {:.3}%",
                out.stage_ppv * 100.0,
                out.stage_npv * 100.0
            );
            println!(
                "Post-Test Probability : {:.3}% ({})",
                out.post_test_prob * 100.0,
                if out.observed_positive {
                    "after T+"
                } else {
                    "after T-"
                }
            );

            if let (Some(pos), Some(tot)) = (positives, tested) {
                let rg = rogan_gladen_prevalence(pos, tot, se, sp)?;
                println!("--- Rogan & Gladen (1978) Prevalence Adjustment ---");
                println!(
                    "Apparent Prevalence   : {:.2}% ({pos}/{tot})",
                    rg.apparent_prevalence * 100.0
                );
                println!(
                    "True Prevalence Est.  : {:.2}% (95% CI {:.2}% .. {:.2}%)",
                    rg.true_prevalence * 100.0,
                    rg.true_ci95_clamped.0 * 100.0,
                    rg.true_ci95_clamped.1 * 100.0
                );
            }
            Ok(())
        }
        Some(EpiSubcommand::Power {
            design,
            p0,
            effect,
            ratio,
            alpha,
            power,
            positivity,
            svg_out,
        }) => {
            let plan = match design.to_ascii_lowercase().as_str() {
                "tnd" | "test-negative" | "ve" => {
                    design_test_negative_study(p0, effect, ratio, alpha, power, positivity)?
                }
                "cohort" | "rct" | "rr" => design_cohort_study(p0, effect, ratio, alpha, power)?,
                "case-control" | "cc" | "or" => {
                    design_case_control_study(p0, effect, ratio, alpha, power)?
                }
                other => {
                    return Err(format!(
                        "Unknown study design '{other}'; use `tnd`, `cohort`, or `case-control`."
                    ))
                }
            };

            println!("=== {} ===", plan.kind.title());
            println!("Target Effect                  : {}", plan.effect_label);
            println!(
                "Proportions                    : p0 = {:.4} | p1 = {:.4} (ratio k = {:.2})",
                plan.p_group_0, plan.p_group_1, plan.allocation_ratio_k
            );
            println!(
                "Uncorrected Sample Size        : n1 = {} | n0 = {} | Total = {}",
                plan.n1_uncorrected,
                plan.n0_uncorrected,
                plan.n1_uncorrected + plan.n0_uncorrected
            );
            println!(
                "Fleiss Continuity-Corrected N  : n1 = {} | n0 = {} | Total = {}",
                plan.n1_fleiss_cc, plan.n0_fleiss_cc, plan.total_fleiss_cc
            );
            if let Some(swabs) = plan.total_ari_swabs_needed {
                println!("Estimated ARI Swabs Needed     : {swabs}");
            }
            if let Some(sp) = svg_out {
                fs::write(&sp, render_power_curve_svg(&plan))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                println!("Saved Power Curve SVG          -> {}", sp.display());
            }
            Ok(())
        }
    }
}

fn fmt_opt(v: Option<f64>) -> String {
    v.map(|x| format!("{x:.3}")).unwrap_or_else(|| String::from("NA"))
}

fn parse_cli_incidence(raw: &str) -> Result<Vec<f64>, String> {
    let mut vals = Vec::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        if let Some(last) = t
            .split([',', '\t', ';', ' '])
            .filter(|s| !s.is_empty())
            .next_back()
        {
            if let Ok(v) = last.parse::<f64>() {
                if v >= 0.0 && v.is_finite() {
                    vals.push(v);
                }
            }
        }
    }
    if vals.len() < 3 {
        return Err(String::from("Expected >= 3 daily incidence numbers."));
    }
    Ok(vals)
}
