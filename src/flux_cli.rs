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

//! Command-line subcommands (`clap` v4) for Tool #20 (`osabt flux`):
//! - `osabt flux gui` — Launches the interactive `egui` FvCB Leaf & Stream Metabolism Studio
//! - `osabt flux aci` — Fits Farquhar–von Caemmerer–Berry (1980) `Vcmax`, `Jmax`, `TPU`, `Rd` to $A\text{–}C_i$ CSV
//! - `osabt flux stream` — Fits Odum (1956) open-channel `GPP`, `ER`, `NEP`, `K600` & exports `.svg`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::flux_leaf::{demo_fvcb_aci_csv, fit_fvcb_aci_curve, parse_aci_csv};
use crate::flux_stream::{
    demo_stream_diel_do_csv, fit_odum_stream_metabolism, parse_stream_do_csv, render_flux_svg,
};

/// Top-level CLI arguments for `osabt flux`.
#[derive(Debug, Args)]
pub struct FluxArgs {
    #[command(subcommand)]
    pub command: Option<FluxSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum FluxSubcommand {
    /// Launch the interactive `egui` Leaf `FvCB` Photosynthesis & Stream Metabolism Studio.
    Gui {
        /// Optional $A\text{–}C_i$ leaf gas-exchange CSV (`ci_ubar,a_net[,ppfd,t_leaf_c]`).
        #[arg(long)]
        aci: Option<PathBuf>,

        /// Optional diel dissolved-oxygen stream logger CSV (`hour,do_mg_l,temp_c,par[,depth_m]`).
        #[arg(long)]
        stream: Option<PathBuf>,
    },

    /// Fit Farquhar–von Caemmerer–Berry (1980) `Vcmax`, `Jmax`, `TPU`, and `Rd` to an $A\text{–}C_i$ curve.
    Aci {
        #[arg(long)]
        input: Option<PathBuf>,

        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Fit Odum (1956) open-channel diel dissolved-oxygen `GPP`, `ER`, `NEP`, and `K600`.
    Stream {
        #[arg(long)]
        input: Option<PathBuf>,

        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

/// Dispatches `osabt flux` CLI subcommands.
pub fn run_flux_cli(args: FluxArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_gui(None, None),
        Some(FluxSubcommand::Gui { aci, stream }) => launch_default_gui(aci, stream),
        Some(FluxSubcommand::Aci { input, svg_out }) => {
            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading A-Ci CSV '{}': {e}", p.display()))?
            } else {
                demo_fvcb_aci_csv()
            };
            let pts = parse_aci_csv(&raw)?;
            let fit = fit_fvcb_aci_curve(&pts)?;
            println!(
                "=== OSABT Flux — Farquhar–von Caemmerer–Berry (1980) Leaf A-Ci Fit ({} Points) ===",
                pts.len()
            );
            println!(
                "  Vcmax = {:.2} umol/m2/s | Jmax = {:.2} umol/m2/s | TPU = {:.2} | Rd = {:.2}",
                fit.vcmax, fit.jmax, fit.tpu, fit.rd
            );
            println!(
                "  Jmax/Vcmax = {:.3} | RMSE = {:.3} | R^2 = {:.5}",
                fit.jmax_to_vcmax_ratio, fit.rmse, fit.r_squared
            );
            for ((p, pred), st) in pts.iter().zip(&fit.predicted_a).zip(&fit.limiting_states) {
                println!(
                    "    Ci={:6.1} ubar | A_obs={:6.2} | A_fit={:6.2} | Regime: {}",
                    p.ci_ubar,
                    p.a_net,
                    pred,
                    st.label()
                );
            }
            if let Some(out) = svg_out {
                let do_pts = parse_stream_do_csv(&demo_stream_diel_do_csv())?;
                let stream_rep = fit_odum_stream_metabolism(&do_pts)?;
                let svg = render_flux_svg(&pts, &fit, &do_pts, &stream_rep, "OSABT Flux");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", out.display()))?;
                println!("Wrote Dual-Panel Flux SVG -> {}", out.display());
            }
            Ok(())
        }
        Some(FluxSubcommand::Stream { input, svg_out }) => {
            let raw = if let Some(p) = input {
                fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading stream CSV '{}': {e}", p.display()))?
            } else {
                demo_stream_diel_do_csv()
            };
            let do_pts = parse_stream_do_csv(&raw)?;
            let stream_rep = fit_odum_stream_metabolism(&do_pts)?;
            println!(
                "=== OSABT Flux — Odum (1956) Diel Stream Metabolism ({} Readings) ===",
                do_pts.len()
            );
            println!(
                "  GPP = {:.3} g O2/m2/d | ER = {:.3} g O2/m2/d | NEP = {:+.3} | P/R = {:.3} | K600 = {:.2} d^-1",
                stream_rep.gpp_g_o2_m2_d,
                stream_rep.er_g_o2_m2_d,
                stream_rep.nep_g_o2_m2_d,
                stream_rep.p_to_r_ratio,
                stream_rep.k600_per_day
            );
            println!(
                "  Trophic State: {} | DO RMSE = {:.4} mg/L",
                stream_rep.trophic_state, stream_rep.do_rmse_mg_l
            );
            if let Some(out) = svg_out {
                let aci_pts = parse_aci_csv(&demo_fvcb_aci_csv())?;
                let aci_fit = fit_fvcb_aci_curve(&aci_pts)?;
                let svg = render_flux_svg(&aci_pts, &aci_fit, &do_pts, &stream_rep, "OSABT Flux");
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", out.display()))?;
                println!("Wrote Dual-Panel Flux SVG -> {}", out.display());
            }
            Ok(())
        }
    }
}

#[cfg(feature = "gui")]
fn launch_default_gui(aci: Option<PathBuf>, stream: Option<PathBuf>) -> Result<(), String> {
    crate::flux_gui::run_flux_gui(aci, stream)
}

#[cfg(not(feature = "gui"))]
fn launch_default_gui(_aci: Option<PathBuf>, _stream: Option<PathBuf>) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `gui` feature; rebuild with `--features gui` or run `osabt flux aci`.",
    ))
}
