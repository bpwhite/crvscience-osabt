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

//! Command-line subcommands (`clap` v4) for Tool #7: `osabt forge3d`.

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand};

use crate::forge3d_bio::{
    crv_voxel_csv_to_mesh, demo_alpha_helix_pdb, generate_diatom_frustule_mesh,
    generate_radiolarian_mesh, parse_pdb_atoms, pdb_atoms_to_printable_mesh, PdbPrintStyle,
};
use crate::forge3d_mesh::{parse_stl_bytes, TriMesh};

#[derive(Debug, Args)]
pub struct Forge3dArgs {
    #[command(subcommand)]
    pub command: Option<Forge3dSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum Forge3dSubcommand {
    /// Launch the interactive desktop `egui` 3D Biomolecule & Morphology Print Studio.
    Gui {
        /// Optional initial `.pdb` structure file.
        #[arg(short, long)]
        input: Option<PathBuf>,
    },

    /// Convert an RCSB `.pdb` macromolecule (or built-in Alpha-Helix) into a watertight `.stl`, `.3mf`, or `.obj`.
    Pdb {
        /// Optional `.pdb` file (uses built-in 16-residue Alpha-Helix if omitted).
        #[arg(short, long)]
        input: Option<PathBuf>,
        /// Output mesh path (`.stl`, `.3mf`, or `.obj`).
        #[arg(short, long, default_value = "biomolecule.stl")]
        output: PathBuf,
        /// Render C-alpha / Phosphorus backbone worm tube instead of full solvent Gaussian surface.
        #[arg(long)]
        backbone: bool,
        /// Voxel contouring resolution along longest axis (`20..96`).
        #[arg(long, default_value_t = 44)]
        resolution: usize,
        /// Target longest print-bed dimension in millimeters.
        #[arg(long, default_value_t = 90.0)]
        scale_mm: f32,
    },

    /// Generate a Gielis (2003) silica diatom frustule valve with hexagonal areolae pores.
    Diatom {
        /// Output mesh path (`.stl`, `.3mf`, or `.obj`).
        #[arg(short, long, default_value = "diatom_frustule.stl")]
        output: PathBuf,
        /// Rotational symmetry lobes `m` (`4` = quadrate, `5` = pentagonal, `6` = Actinoptychus hex, `8` = octaradiate).
        #[arg(short = 'm', long, default_value_t = 6.0)]
        lobes: f32,
        /// Relative hexagonal pore spacing (`0.18..0.42`).
        #[arg(long, default_value_t = 0.28)]
        pore_pitch: f32,
        /// Voxel contouring resolution (`32..88`).
        #[arg(long, default_value_t = 48)]
        resolution: usize,
        /// Target diameter on print bed in millimeters.
        #[arg(long, default_value_t = 95.0)]
        scale_mm: f32,
    },

    /// Generate a Haeckel-inspired radiolarian perforated cortical test with 12 icosahedral spines.
    Radiolarian {
        /// Output mesh path (`.stl`, `.3mf`, or `.obj`).
        #[arg(short, long, default_value = "radiolarian.stl")]
        output: PathBuf,
        #[arg(long, default_value_t = 46)]
        resolution: usize,
        #[arg(long, default_value_t = 90.0)]
        scale_mm: f32,
    },

    /// Convert a CRV Sim 3D voxel CSV/TSV (`x,y,z,density`) into a smoothed watertight `.stl` / `.3mf`.
    Voxels {
        /// Input voxel CSV/TSV file.
        #[arg(short, long)]
        input: PathBuf,
        /// Output mesh path (`.stl`, `.3mf`, or `.obj`).
        #[arg(short, long, default_value = "crv_voxels.stl")]
        output: PathBuf,
        /// Isosurface density threshold.
        #[arg(long, default_value_t = 0.5)]
        threshold: f32,
        #[arg(long, default_value_t = 100.0)]
        scale_mm: f32,
    },

    /// Audit an existing Binary or ASCII `.stl` file for watertight 2-manifold edges, volume, and PLA/SLA mass.
    InspectStl {
        /// Path to `.stl` file.
        input: PathBuf,
    },
}

pub fn run_forge3d_cli(args: Forge3dArgs) -> Result<(), String> {
    match args.command {
        None => {
            #[cfg(feature = "gui")]
            {
                return crate::forge3d_gui::run_forge3d_gui(None);
            }
            #[cfg(not(feature = "gui"))]
            {
                return Err(String::from(
                    "Build with `--features gui` or run `osabt forge3d pdb` / `osabt forge3d diatom`.",
                ));
            }
        }
        Some(Forge3dSubcommand::Gui { input }) => {
            #[cfg(feature = "gui")]
            {
                crate::forge3d_gui::run_forge3d_gui(input)
            }
            #[cfg(not(feature = "gui"))]
            {
                let _ = input;
                Err(String::from("GUI feature is disabled in this build."))
            }
        }
        Some(Forge3dSubcommand::Pdb {
            input,
            output,
            backbone,
            resolution,
            scale_mm,
        }) => {
            let (name, txt) = match input {
                Some(p) => {
                    let s = fs::read_to_string(&p)
                        .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                    let stem = p
                        .file_stem()
                        .map(|x| x.to_string_lossy().to_string())
                        .unwrap_or_else(|| String::from("PDB_Model"));
                    (stem, s)
                }
                None => (String::from("pCRV_AlphaHelix"), demo_alpha_helix_pdb()),
            };
            let atoms = parse_pdb_atoms(&txt, false)?;
            let style = if backbone {
                PdbPrintStyle::BackboneTube
            } else {
                PdbPrintStyle::GaussianSurface
            };
            let mesh = pdb_atoms_to_printable_mesh(&name, &atoms, style, resolution, scale_mm)?;
            save_and_report_mesh(&mesh, &output)
        }
        Some(Forge3dSubcommand::Diatom {
            output,
            lobes,
            pore_pitch,
            resolution,
            scale_mm,
        }) => {
            let mesh = generate_diatom_frustule_mesh(lobes, pore_pitch, resolution, scale_mm);
            save_and_report_mesh(&mesh, &output)
        }
        Some(Forge3dSubcommand::Radiolarian {
            output,
            resolution,
            scale_mm,
        }) => {
            let mesh = generate_radiolarian_mesh(resolution, scale_mm);
            save_and_report_mesh(&mesh, &output)
        }
        Some(Forge3dSubcommand::Voxels {
            input,
            output,
            threshold,
            scale_mm,
        }) => {
            let csv = fs::read_to_string(&input)
                .map_err(|e| format!("Failed reading '{}': {e}", input.display()))?;
            let stem = input
                .file_stem()
                .map(|x| x.to_string_lossy().to_string())
                .unwrap_or_else(|| String::from("CRV_Voxels"));
            let mesh = crv_voxel_csv_to_mesh(&stem, &csv, threshold, scale_mm)?;
            save_and_report_mesh(&mesh, &output)
        }
        Some(Forge3dSubcommand::InspectStl { input }) => {
            let bytes = fs::read(&input)
                .map_err(|e| format!("Failed reading '{}': {e}", input.display()))?;
            let stem = input
                .file_stem()
                .map(|s| s.to_string_lossy().to_string())
                .unwrap_or_else(|| String::from("Mesh"));
            let mesh = parse_stl_bytes(&stem, &bytes)?;
            print_mesh_audit(&mesh);
            Ok(())
        }
    }
}

fn save_and_report_mesh(mesh: &TriMesh, output: &PathBuf) -> Result<(), String> {
    let ext = output
        .extension()
        .and_then(|s| s.to_str())
        .unwrap_or("stl")
        .to_ascii_lowercase();

    match ext.as_str() {
        "3mf" => fs::write(output, mesh.to_3mf_package())
            .map_err(|e| format!("Failed writing 3MF '{}': {e}", output.display()))?,
        "obj" => fs::write(output, mesh.to_wavefront_obj())
            .map_err(|e| format!("Failed writing OBJ '{}': {e}", output.display()))?,
        "ascii.stl" => fs::write(output, mesh.to_ascii_stl())
            .map_err(|e| format!("Failed writing ASCII STL '{}': {e}", output.display()))?,
        _ => fs::write(output, mesh.to_binary_stl())
            .map_err(|e| format!("Failed writing Binary STL '{}': {e}", output.display()))?,
    }

    print_mesh_audit(mesh);
    println!("Saved 3D-printable mesh -> {}", output.display());
    Ok(())
}

fn print_mesh_audit(mesh: &TriMesh) {
    let r = mesh.inspect_printability();
    println!("=== OSABT Forge3D Slicer Audit: {} ===", mesh.name);
    println!(
        "Vertices / Triangles : {} / {}",
        r.vertex_count, r.triangle_count
    );
    println!(
        "Dimensions (mm)      : {:.2} × {:.2} × {:.2}",
        r.dimensions_mm[0], r.dimensions_mm[1], r.dimensions_mm[2]
    );
    println!(
        "Manifold Status      : {} (bad edges: {}, Euler V-E+F = {})",
        if r.is_watertight_manifold {
            "WATERTIGHT 2-MANIFOLD"
        } else {
            "NON-MANIFOLD / OPEN"
        },
        r.boundary_or_nonmanifold_edges,
        r.euler_characteristic
    );
    println!(
        "Volume & Area        : {:.2} mm³ | {:.2} mm²",
        r.signed_volume_mm3, r.surface_area_mm2
    );
    println!(
        "Material Estimate    : {:.1} g FDM PLA (20% infill) | {:.1} g SLA Resin (solid)",
        r.est_fdm_pla_grams, r.est_sla_resin_grams
    );
}
