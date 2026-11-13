// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #8: `osabt morpho`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::morpho_gpa::{
    demo_salmonid_landmarks_tps, parse_tps_landmarks, render_gpa_morphospace_svg,
    run_generalized_procrustes_2d,
};
use crate::morpho_vision::{
    analyze_leaf_herbivory, count_petri_colonies, estimate_cfu_per_ml, render_petri_cfu_svg,
    RgbImage,
};

#[derive(Debug, Args)]
pub struct MorphoArgs {
    #[command(subcommand)]
    pub command: Option<MorphoSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum MorphoSubcommand {
    /// Launch the interactive desktop `egui` Colony Counter, Leaf Herbivory & 2D GPA Studio.
    Gui {
        /// Optional `.bmp`, `.ppm`, `.pgm`, or `.tps` file to open on startup.
        #[arg(short, long)]
        input: Option<PathBuf>,
    },

    /// Automated Petri dish CFU counter with touching-colony distance-peak splitting & `CFU/mL` titer.
    Cfu {
        /// Input `.bmp` or `.ppm`/`.pgm` image (uses built-in 18-CFU plate with doublets if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Detect dark colonies on bright-field agar instead of light colonies on dark background.
        #[arg(long)]
        dark_colonies: bool,
        /// Inner dish circular ROI fraction (`0.50..0.99`) masking Petri rim glare.
        #[arg(long, default_value_t = 0.84)]
        roi: f64,
        /// Plated liquid volume in mL (for `CFU/mL` titer estimation).
        #[arg(long, default_value_t = 0.10)]
        volume_ml: f64,
        /// Dilution factor (e.g. `1e-5` for a `10^-5` dilution plate).
        #[arg(long, default_value_t = 1e-5)]
        dilution: f64,
        /// Optional annotated `.svg` overlay output path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },

    /// Botanical leaf-area & enclosed beetle/caterpillar herbivory quantifier via Excess Green (`ExG`).
    Leaf {
        /// Input `.bmp` or `.ppm`/`.pgm` leaf scan (uses built-in herbivorized leaf if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Spatial scale calibration in pixels per millimeter (`px / mm`).
        #[arg(long, default_value_t = 4.0)]
        px_per_mm: f64,
        /// Optional manual `ExG = 2G - R - B` threshold override.
        #[arg(long)]
        exg_threshold: Option<i32>,
    },

    /// Rohlf & Slice (1990) 2D Generalized Procrustes Analysis (`GPA`) & PCA on `.tps` homologous landmarks.
    Gpa {
        /// Input Rohlf `.tps` landmark file (uses built-in 6-specimen *Oncorhynchus* dataset if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Optional 2-panel Procrustes + `PC1/PC2` morphospace `.svg` output path.
        #[arg(long)]
        svg_out: Option<PathBuf>,
    },
}

pub fn run_morpho_cli(args: MorphoArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "gui")]
            {
                return crate::morpho_gui::run_morpho_gui(None);
            }
            #[cfg(not(feature = "gui"))]
            {
                return Err(String::from(
                    "Build with `--features gui` or run `osabt morpho cfu` / `leaf` / `gpa`.",
                ));
            }
        }
        Some(MorphoSubcommand::Gui { input }) => {
            #[cfg(feature = "gui")]
            {
                crate::morpho_gui::run_morpho_gui(input)
            }
            #[cfg(not(feature = "gui"))]
            {
                let _ = input;
                Err(String::from("GUI feature is disabled in this build."))
            }
        }
        Some(MorphoSubcommand::Cfu {
            input,
            dark_colonies,
            roi,
            volume_ml,
            dilution,
            svg_out,
        }) => {
            let img = match input {
                Some(p) => {
                    let b = fs::read(&p)
                        .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                    RgbImage::from_bytes(&b)?
                }
                None => RgbImage::demo_petri_cfu_plate(),
            };
            let rep = count_petri_colonies(&img, dark_colonies, roi, 8, None);
            let (mle, lo, hi) =
                estimate_cfu_per_ml(&[(rep.total_cfu_count, volume_ml, dilution)])?;

            println!("=== OSABT Petri Dish CFU Report ===");
            println!("Total CFU Count       : {}", rep.total_cfu_count);
            println!("Isolated Colonies     : {}", rep.isolated_colonies);
            println!("Split Touching Peaks  : {}", rep.split_cluster_colonies);
            println!("Otsu Threshold        : {}", rep.threshold_used);
            println!("Mean Colony Diameter  : {:.2} px", rep.mean_colony_diameter_px);
            println!(
                "Estimated Titer       : {:.3e} CFU/mL (95% CI {:.3e} .. {:.3e})",
                mle, lo, hi
            );
            if let Some(sp) = svg_out {
                fs::write(&sp, render_petri_cfu_svg(&img, &rep))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                println!("Saved CFU overlay SVG -> {}", sp.display());
            }
            Ok(())
        }
        Some(MorphoSubcommand::Leaf {
            input,
            px_per_mm,
            exg_threshold,
        }) => {
            let img = match input {
                Some(p) => {
                    let b = fs::read(&p)
                        .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                    RgbImage::from_bytes(&b)?
                }
                None => RgbImage::demo_herbivory_leaf(),
            };
            let rep = analyze_leaf_herbivory(&img, px_per_mm, exg_threshold)?;
            println!("=== OSABT Botanical Leaf Area & Herbivory Report ===");
            println!("Intact Green Lamina   : {:.3} cm² ({} px)", rep.lamina_area_cm2, rep.green_lamina_px);
            println!(
                "Interior Shot-Holes   : {} holes ({:.4} cm², {:.2}% interior herbivory)",
                rep.enclosed_hole_count, rep.enclosed_herbivory_cm2, rep.interior_herbivory_pct
            );
            println!(
                "Convex Hull Envelope  : {:.3} cm² ({:.2}% total marginal + interior defoliation)",
                rep.convex_envelope_cm2, rep.total_defoliation_envelope_pct
            );
            Ok(())
        }
        Some(MorphoSubcommand::Gpa { input, svg_out }) => {
            let tps_txt = match input {
                Some(p) => fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?,
                None => demo_salmonid_landmarks_tps().to_string(),
            };
            let specs = parse_tps_landmarks(&tps_txt)?;
            let gpa = run_generalized_procrustes_2d(&specs)?;

            println!(
                "=== OSABT 2D Generalized Procrustes Analysis ({} specimens, {} landmarks) ===",
                gpa.specimen_ids.len(),
                gpa.num_landmarks
            );
            println!(
                "Tangent PCA Variance  : PC1 = {:.2}% | PC2 = {:.2}%",
                gpa.pc1_variance_pct, gpa.pc2_variance_pct
            );
            println!("specimen_id\tcentroid_size\tprocrustes_dist\tpc1\tpc2");
            for i in 0..gpa.specimen_ids.len() {
                println!(
                    "{}\t{:.4}\t{:.5}\t{:+.5}\t{:+.5}",
                    gpa.specimen_ids[i],
                    gpa.centroid_sizes[i],
                    gpa.procrustes_distances_to_mean[i],
                    gpa.pc12_scores[i][0],
                    gpa.pc12_scores[i][1]
                );
            }
            if let Some(sp) = svg_out {
                fs::write(&sp, render_gpa_morphospace_svg(&gpa))
                    .map_err(|e| format!("Failed writing SVG '{}': {e}", sp.display()))?;
                println!("Saved GPA + PCA SVG   -> {}", sp.display());
            }
            Ok(())
        }
    }
}
