// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #15 (`osabt wildlife`):
//! - `osabt wildlife tui` — Interactive 4-tab `Ratatui` Occupancy & Capture-Recapture Studio
//! - `osabt wildlife occupancy` — Fits MacKenzie et al. (2002) `Psi(.)p(.)` & covariate `AICc` models
//! - `osabt wildlife markrecap` — Runs Chapman–Petersen 2-sample or multi-occasion `M0/Mt/Mh/Schnabel` + `CJS` survival
//! - `osabt wildlife plan` — Computes MacKenzie & Royle (2005) optimal survey visits `K*` & required sites `S*`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::wildlife_markrecap::{
    analyze_closed_capture_histories, calc_chapman_lincoln_petersen, demo_pnw_capture_history_csv,
    estimate_cjs_apparent_survival, parse_capture_history_csv, render_wildlife_population_svg,
};
use crate::wildlife_occupancy::{
    demo_pnw_occupancy_csv, fit_mackenzie_occupancy_suite, parse_occupancy_csv,
    plan_occupancy_survey_effort,
};

/// Top-level CLI arguments for `osabt wildlife`.
#[derive(Debug, Args)]
pub struct WildlifeArgs {
    #[command(subcommand)]
    pub command: Option<WildlifeSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum WildlifeSubcommand {
    /// Launch the interactive 4-tab `Ratatui` Wildlife Occupancy & Capture-Recapture Studio.
    Tui {
        /// Optional occupancy detection-history CSV (`site_id,v1,...,vK[,covariate]`).
        #[arg(long)]
        occupancy: Option<PathBuf>,

        /// Optional capture-history CSV (`animal_id,t1,...,tK`).
        #[arg(long)]
        captures: Option<PathBuf>,
    },

    /// Fit MacKenzie et al. (2002) imperfect-detection occupancy models (`Psi(.)p(.)`, `Psi(cov)p(.)`, `Psi(.)p(cov)`).
    Occupancy {
        /// Optional occupancy CSV/TSV file (omit for built-in Oregon Coast Range Red Tree Vole dataset).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Optional path to export 2-panel Occupancy + Capture-Recapture `.svg` figure.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Estimate closed-population abundance (`Chapman–Petersen`, `M0`, `Mt`, `Chao Mh`, `Schnabel`) and open `CJS` survival.
    Markrecap {
        /// Optional multi-occasion capture-history CSV (`animal_id,t1,...,tK`).
        #[arg(long)]
        input: Option<PathBuf>,

        /// Optional 2-sample Lincoln–Petersen marked count `M` (use with `--captured` and `--recaptured`).
        #[arg(long)]
        marked: Option<usize>,

        /// Optional 2-sample Lincoln–Petersen second-sample captured count `C`.
        #[arg(long)]
        captured: Option<usize>,

        /// Optional 2-sample Lincoln–Petersen recaptured count `R`.
        #[arg(long)]
        recaptured: Option<usize>,
    },

    /// Plan repeat visits `K*` and site count `S*` under imperfect detection (MacKenzie & Royle 2005).
    Plan {
        /// Expected site occupancy probability `psi` in `[0.01, 0.99]`.
        #[arg(long, default_value_t = 0.55)]
        psi: f64,

        /// Expected single-visit detection probability `p` in `[0.02, 0.99]`.
        #[arg(long, default_value_t = 0.40)]
        p: f64,

        /// Maximum false-absence probability `(1 - p)^K` at an occupied site.
        #[arg(long, default_value_t = 0.05)]
        alpha: f64,

        /// Target standard error `SE(psi_hat)` for regional occupancy estimation.
        #[arg(long, default_value_t = 0.06)]
        target_se: f64,

        /// Ratio of initial site-setup/travel cost to a single repeat survey cost.
        #[arg(long, default_value_t = 2.0)]
        cost_ratio: f64,
    },
}

/// Dispatches `osabt wildlife` CLI subcommands.
pub fn run_wildlife_cli(args: WildlifeArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_tui(None, None),
        Some(WildlifeSubcommand::Tui {
            occupancy,
            captures,
        }) => launch_default_tui(occupancy, captures),
        Some(WildlifeSubcommand::Occupancy { input, svg_out }) => {
            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading occupancy file '{}': {e}", p.display()))?
            } else {
                demo_pnw_occupancy_csv().to_string()
            };
            let rows = parse_occupancy_csv(&raw)?;
            let rep = fit_mackenzie_occupancy_suite(&rows)?;
            println!(
                "=== OSABT Wildlife — MacKenzie et al. (2002) Site Occupancy ({} Sites, Naive = {:.3}) ===",
                rep.num_sites, rep.naive_occupancy
            );
            println!(
                "Finite-Sample Empirical Bayes Occupied Sites: {:.2} / {}",
                rep.finite_sample_occupied_est, rep.num_sites
            );
            for (i, m) in rep.models.iter().enumerate() {
                let tag = if i == rep.best_model_index { "*" } else { " " };
                println!(
                    " {tag} {:14} | psi={:.3} [{:.2}–{:.2}] | p={:.3} [{:.2}–{:.2}] | slope={:+.2} | AICc={:.2} (w={:.1}%)",
                    m.model.label(),
                    m.psi_hat,
                    m.psi_ci95.0,
                    m.psi_ci95.1,
                    m.p_hat,
                    m.p_ci95.0,
                    m.p_ci95.1,
                    m.covariate_logit_slope,
                    m.aicc,
                    m.akaike_weight * 100.0
                );
            }
            if let Some(out) = svg_out {
                let (_, demo_h) = parse_capture_history_csv(demo_pnw_capture_history_csv())?;
                let closed = analyze_closed_capture_histories(&demo_h)?;
                let svg = render_wildlife_population_svg(&rep, &closed, "OSABT Wildlife Occupancy & Abundance");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", out.display()))?;
                println!("Wrote diagnostic SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(WildlifeSubcommand::Markrecap {
            input,
            marked,
            captured,
            recaptured,
        }) => {
            if let (Some(m), Some(c), Some(r)) = (marked, captured, recaptured) {
                let cp = calc_chapman_lincoln_petersen(m, c, r)?;
                println!("=== Chapman (1951) / Seber (1970) 2-Sample Lincoln–Petersen ===");
                println!(
                    "  M={m}, C={c}, R={r} (Distinct M_t+1={}) => N_hat = {:.2} ± {:.2} (95% CI [{:.1}, {:.1}])",
                    cp.distinct_observed_mt1, cp.chapman_n_hat, cp.seber_se, cp.ci95.0, cp.ci95.1
                );
                if input.is_none() {
                    return Ok(());
                }
            }

            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading capture CSV '{}': {e}", p.display()))?
            } else {
                demo_pnw_capture_history_csv().to_string()
            };
            let (_, hists) = parse_capture_history_csv(&raw)?;
            let closed = analyze_closed_capture_histories(&hists)?;
            println!(
                "=== Closed-Population Capture-Recapture (M_t+1 = {} Individuals, K = {} Occasions) ===",
                closed.distinct_individuals_mt1, closed.occasions_k
            );
            for est in [
                &closed.m0_equal_catchability,
                &closed.mt_darroch_time_varying,
                &closed.mh_chao1987_heterogeneity,
                &closed.schnabel_multicensus,
            ] {
                println!(
                    "  {:42} : N_hat = {:6.2} ± {:5.2}  95% CI [{:5.1}, {:5.1}]  (p={:.3})",
                    est.model_name, est.n_hat, est.se, est.ci95.0, est.ci95.1, est.mean_capture_p
                );
            }

            if closed.occasions_k >= 3 {
                let cjs = estimate_cjs_apparent_survival(&hists)?;
                println!("\n=== Cormack–Jolly–Seber (Lebreton et al. 1992) Open Survival ===");
                println!(
                    "  Apparent Survival phi = {:.3} ± {:.3} [{:.3}–{:.3}] | Recapture p = {:.3} ± {:.3} [{:.3}–{:.3}]",
                    cjs.phi_hat,
                    cjs.phi_se,
                    cjs.phi_ci95.0,
                    cjs.phi_ci95.1,
                    cjs.p_hat,
                    cjs.p_se,
                    cjs.p_ci95.0,
                    cjs.p_ci95.1
                );
            }
            Ok(())
        }
        Some(WildlifeSubcommand::Plan {
            psi,
            p,
            alpha,
            target_se,
            cost_ratio,
        }) => {
            let plan = plan_occupancy_survey_effort(psi, p, alpha, target_se, cost_ratio)?;
            println!("=== MacKenzie & Royle (2005) Occupancy Survey Design ===");
            println!(
                "  Assumed psi={:.2}, p={:.2}, False-Absence alpha={:.3}, Site/Visit Cost={:.1}x",
                plan.assumed_psi,
                plan.assumed_p,
                plan.target_false_absence_alpha,
                cost_ratio
            );
            println!(
                "  Min Visits for Presence Certainty : K_min = {} visits/site (cumulative p* = {:.1}%)",
                plan.min_visits_for_presence_certainty,
                plan.cumulative_detection_at_k_min * 100.0
            );
            println!(
                "  Precision-Optimal Allocation      : S* = {} sites x K* = {} visits ({} total surveys, SE = {:.4})",
                plan.required_sites_for_target_se,
                plan.optimal_k_for_precision,
                plan.total_survey_visits,
                plan.achieved_psi_se
            );
            Ok(())
        }
    }
}

#[cfg(feature = "tui")]
fn launch_default_tui(
    occ_path: Option<PathBuf>,
    cap_path: Option<PathBuf>,
) -> Result<(), String> {
    crate::wildlife_tui::run_wildlife_tui(occ_path, cap_path)
}

#[cfg(not(feature = "tui"))]
fn launch_default_tui(
    _occ_path: Option<PathBuf>,
    _cap_path: Option<PathBuf>,
) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `tui` feature; rebuild with `--features tui` or run `osabt wildlife occupancy`.",
    ))
}
