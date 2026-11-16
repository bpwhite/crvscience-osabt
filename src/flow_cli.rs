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

//! Command-line subcommands (`clap` v4) for Tool #18 (`osabt flow`):
//! - `osabt flow gui` — Launches the interactive `egui` FCS 3.1 Compensation & `Logicle` Gating Studio
//! - `osabt flow inspect` — Prints `FCS 3.0/3.1` channel metadata, `$TOT`, and `$SPILLOVER` matrix
//! - `osabt flow compensate` — Applies `$SPILLOVER` matrix inversion and writes a compensated `.fcs` 3.1 file
//! - `osabt flow gate` — Evaluates `Logicle` quadrant gates (`% Parent`, `MFI`, `rCV%`) and exports `.svg`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::flow_fcs::{
    compensate_fcs_dataset, demo_fcs_31_bytes, parse_fcs_bytes, write_fcs_31_bytes, FcsDataset,
};
use crate::flow_gate::{
    evaluate_quadrant_gates, render_flow_biaxial_svg, LogicleParams, LogicleScale,
};

/// Top-level CLI arguments for `osabt flow`.
#[derive(Debug, Args)]
pub struct FlowArgs {
    #[command(subcommand)]
    pub command: Option<FlowSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum FlowSubcommand {
    /// Launch the interactive `egui` Flow Cytometry `.fcs 3.1` Compensation & `Logicle` Studio.
    Gui {
        /// Optional binary `.fcs` file (omit for built-in 400-event 4-channel `FCS 3.1` demo).
        #[arg(long)]
        input: Option<PathBuf>,
    },

    /// Inspect an `FCS 3.0/3.1` binary file's header, channels (`$PnN/$PnS`), and `$SPILLOVER` matrix.
    Inspect {
        #[arg(long)]
        input: Option<PathBuf>,
    },

    /// Invert `$SPILLOVER` (`S^-1`) and export a compensated binary `FCS 3.1` file.
    Compensate {
        #[arg(long)]
        input: Option<PathBuf>,

        /// Output path for the compensated binary `.fcs` file.
        #[arg(long)]
        fcs_out: Option<PathBuf>,
    },

    /// Evaluate `Logicle` biexponential quadrant gates (`% Parent`, `MFI`, `rCV%`) and export `.svg`.
    Gate {
        #[arg(long)]
        input: Option<PathBuf>,

        /// Zero-based X-axis channel index (default `2`, e.g. `FITC-A`).
        #[arg(long, default_value_t = 2)]
        x_ch: usize,

        /// Zero-based Y-axis channel index (default `3`, e.g. `PE-A`).
        #[arg(long, default_value_t = 3)]
        y_ch: usize,

        /// Logicle display threshold `[0.0..1.0]` for X and Y quadrant dividers.
        #[arg(long, default_value_t = 0.38)]
        threshold: f64,

        /// Optional `.svg` biaxial scatter plot output path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

/// Dispatches `osabt flow` CLI subcommands.
pub fn run_flow_cli(args: FlowArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_gui(None),
        Some(FlowSubcommand::Gui { input }) => launch_default_gui(input),
        Some(FlowSubcommand::Inspect { input }) => {
            let ds = load_or_demo_fcs(input.as_ref())?;
            println!(
                "=== OSABT Flow — {} Header & Channel Summary ({} Events, {} Parameters) ===",
                ds.version,
                ds.num_events,
                ds.channels.len()
            );
            for ch in &ds.channels {
                println!(
                    "  P{}: {:10} | Stain: {:18} | Bits: {} | Range: {:.0}",
                    ch.index_1based, ch.short_name, ch.stain_name, ch.bits, ch.range
                );
            }
            if let Some(sp) = &ds.spillover {
                println!(
                    "  $SPILLOVER ({}x{}): {:?}",
                    sp.channel_names.len(),
                    sp.channel_names.len(),
                    sp.channel_names
                );
            }
            Ok(())
        }
        Some(FlowSubcommand::Compensate { input, fcs_out }) => {
            let ds = load_or_demo_fcs(input.as_ref())?;
            let (comp, inv) = compensate_fcs_dataset(&ds)?;
            println!("=== OSABT Flow — Spectral Spillover Compensation (S^-1) ===");
            for row in &inv {
                let s = row
                    .iter()
                    .map(|v| format!("{v:+.5}"))
                    .collect::<Vec<_>>()
                    .join("  ");
                println!("  [ {s} ]");
            }
            if let Some(out) = fcs_out {
                let bytes = write_fcs_31_bytes(&comp);
                fs::write(&out, bytes)
                    .map_err(|e| format!("Failed writing compensated FCS '{}': {e}", out.display()))?;
                println!("Wrote compensated FCS 3.1 -> {}", out.display());
            }
            Ok(())
        }
        Some(FlowSubcommand::Gate {
            input,
            x_ch,
            y_ch,
            threshold,
            svg_out,
        }) => {
            let raw = load_or_demo_fcs(input.as_ref())?;
            let ds = compensate_fcs_dataset(&raw).map(|(c, _)| c).unwrap_or(raw);
            let scale = LogicleScale::new(LogicleParams::default());
            let quads = evaluate_quadrant_gates(&ds, x_ch, y_ch, threshold, threshold, &scale)?;
            println!(
                "=== OSABT Flow — Logicle Quadrant Gating ({} vs {}, n={}) ===",
                ds.channels.get(x_ch).map(|c| c.short_name.as_str()).unwrap_or("X"),
                ds.channels.get(y_ch).map(|c| c.short_name.as_str()).unwrap_or("Y"),
                ds.num_events
            );
            for q in &quads {
                println!(
                    "  {:20} : {:4} events ({:5.1}%) | X_med={:7.0} | Y_med={:7.0} | rCV_X={:5.1}%",
                    q.gate_name,
                    q.event_count,
                    q.percent_of_parent,
                    q.x_median,
                    q.y_median,
                    q.x_robust_cv_pct
                );
            }
            if let Some(out) = svg_out {
                let svg = render_flow_biaxial_svg(
                    &ds,
                    x_ch,
                    y_ch,
                    threshold,
                    threshold,
                    &scale,
                    "OSABT Flow Biaxial Gate",
                );
                fs::write(&out, svg)
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", out.display()))?;
                println!("Wrote Flow Biaxial SVG -> {}", out.display());
            }
            Ok(())
        }
    }
}

fn load_or_demo_fcs(path: Option<&PathBuf>) -> Result<FcsDataset, String> {
    let bytes = if let Some(p) = path {
        fs::read(p).map_err(|e| format!("Failed reading FCS '{}': {e}", p.display()))?
    } else {
        demo_fcs_31_bytes()
    };
    parse_fcs_bytes(&bytes)
}

#[cfg(feature = "gui")]
fn launch_default_gui(input: Option<PathBuf>) -> Result<(), String> {
    crate::flow_gui::run_flow_gui(input)
}

#[cfg(not(feature = "gui"))]
fn launch_default_gui(_input: Option<PathBuf>) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `gui` feature; rebuild with `--features gui` or run `osabt flow gate`.",
    ))
}
