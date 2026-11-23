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

//! Pure-Rust cross-platform `eframe` / `egui` desktop 3D Print Studio for Tool #7: `osabt forge3d`.
//! Features an interactive mouse-orbiting/panning Blinn-Phong + rim + depth-cued 3D viewport
//! (perspective & orthographic), build-plate grid & XYZ triad, live manifold & filament mass audit,
//! and 1-click Binary `.stl`, `.3mf`, and `.obj` exports.

use std::fs;
use std::path::PathBuf;

use eframe::egui;
use egui::{Color32, Pos2, Sense, Shape, Stroke, Vec2};

use crate::forge3d_bio::{
    crv_voxel_csv_to_mesh, demo_alpha_helix_pdb, generate_diatom_frustule_mesh,
    generate_radiolarian_mesh, parse_pdb_atoms, pdb_atoms_to_printable_mesh, PdbPrintStyle,
};
use crate::forge3d_mesh::{MeshPrintReport, TriMesh};
use crate::gui_io::{
    draw_status_banner, pick_open_file, pick_save_file, write_bytes_reported, write_text_reported,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum ForgeSourceKind {
    PdbMolecule,
    GielisDiatom,
    Radiolarian,
    CrvVoxelCsv,
}

/// Interactive `egui` application state for `OSABT Forge3D Print Studio`.
pub struct Forge3dStudioApp {
    source_kind: ForgeSourceKind,
    pdb_style: PdbPrintStyle,
    file_path_input: String,
    diatom_lobes: f32,
    diatom_pore_spacing: f32,
    grid_res: usize,
    target_size_mm: f32,
    yaw_rad: f32,
    pitch_rad: f32,
    zoom: f32,
    pan_offset: Vec2,
    use_perspective: bool,
    draw_wireframe: bool,
    mesh_dirty: bool,
    mesh: TriMesh,
    report: MeshPrintReport,
    status_bar: String,
    status_is_error: bool,
    about_open: bool,
}

impl Forge3dStudioApp {
    pub fn new(initial_pdb: Option<PathBuf>) -> Self {
        let mut app = Self {
            source_kind: ForgeSourceKind::GielisDiatom,
            pdb_style: PdbPrintStyle::GaussianSurface,
            file_path_input: initial_pdb
                .as_ref()
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_default(),
            diatom_lobes: 6.0,
            diatom_pore_spacing: 0.28,
            grid_res: 44,
            target_size_mm: 90.0,
            yaw_rad: 0.65,
            pitch_rad: 0.52,
            zoom: 1.0,
            pan_offset: Vec2::ZERO,
            use_perspective: true,
            draw_wireframe: false,
            mesh_dirty: false,
            mesh: TriMesh::new("Init"),
            report: TriMesh::new("Init").inspect_printability(),
            status_bar: String::from(
                "Drag to orbit | Middle-drag to pan | Scroll to zoom | R to reset camera | Export watertight Binary .stl or .3mf.",
            ),
            status_is_error: false,
            about_open: false,
        };
        if initial_pdb.is_some() {
            app.source_kind = ForgeSourceKind::PdbMolecule;
        }
        app.rebuild_mesh();
        app
    }

    /// Resets orbit angles, pan offset, and zoom back to the default isometric print-bed view.
    pub fn reset_camera(&mut self) {
        self.yaw_rad = 0.65;
        self.pitch_rad = 0.52;
        self.zoom = 1.0;
        self.pan_offset = Vec2::ZERO;
    }

    fn rebuild_mesh(&mut self) {
        let mut read_err = false;
        let built = match self.source_kind {
            ForgeSourceKind::PdbMolecule => {
                let pdb_txt = if self.file_path_input.trim().is_empty() {
                    demo_alpha_helix_pdb()
                } else {
                    match fs::read_to_string(&self.file_path_input) {
                        Ok(t) => t,
                        Err(e) => {
                            read_err = true;
                            self.status_is_error = true;
                            self.status_bar = format!("PDB read error: {e} (showing demo helix)");
                            demo_alpha_helix_pdb()
                        }
                    }
                };
                parse_pdb_atoms(&pdb_txt, false).and_then(|atoms| {
                    pdb_atoms_to_printable_mesh(
                        "PDB_Biomolecule",
                        &atoms,
                        self.pdb_style,
                        self.grid_res,
                        self.target_size_mm,
                    )
                })
            }
            ForgeSourceKind::GielisDiatom => Ok(generate_diatom_frustule_mesh(
                self.diatom_lobes,
                self.diatom_pore_spacing,
                self.grid_res,
                self.target_size_mm,
            )),
            ForgeSourceKind::Radiolarian => Ok(generate_radiolarian_mesh(
                self.grid_res,
                self.target_size_mm,
            )),
            ForgeSourceKind::CrvVoxelCsv => match fs::read_to_string(&self.file_path_input) {
                Ok(csv) => crv_voxel_csv_to_mesh("CRV_Sim_Voxels", &csv, 0.5, self.target_size_mm),
                Err(e) => Err(format!("Provide a valid voxel CSV path (`x,y,z,val`): {e}")),
            },
        };

        match built {
            Ok(m) => {
                self.report = m.inspect_printability();
                self.mesh = m;
                if !read_err {
                    self.status_is_error = false;
                }
            }
            Err(err) => {
                self.status_is_error = true;
                self.status_bar = err;
            }
        }
    }
}

/// Launches the native cross-platform `OSABT Forge3D Studio` desktop window.
pub fn run_forge3d_gui(initial_pdb: Option<PathBuf>) -> Result<(), String> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 780.0])
            .with_icon(crate::brand::osabt_icon_data())
            .with_title("OSABT Forge3D — Biomolecule & Morphology-to-STL Studio (GPL-3.0)"),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-forge3d-studio",
        options,
        Box::new(move |_cc| Ok(Box::new(Forge3dStudioApp::new(initial_pdb)))),
    )
    .map_err(|e| format!("Failed starting Forge3D desktop window: {e}"))
}

impl eframe::App for Forge3dStudioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        crate::brand::draw_standalone_footer(ctx, &mut self.about_open);
        self.draw_panels(ctx);
    }
}

impl Forge3dStudioApp {
    fn track_slider_rebuild(resp: &egui::Response, changed: &mut bool, mesh_dirty: &mut bool) {
        if resp.drag_stopped() || (resp.changed() && !resp.dragged()) {
            *changed = true;
            *mesh_dirty = false;
        } else if resp.changed() {
            *mesh_dirty = true;
        }
    }

    /// Renders all Forge3D Studio panels onto `ctx`.
    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        if !ctx.wants_keyboard_input()
            && ctx.input(|i| {
                !i.modifiers.command
                    && !i.modifiers.ctrl
                    && !i.modifiers.alt
                    && i.key_pressed(egui::Key::R)
            })
        {
            self.reset_camera();
        }

        egui::TopBottomPanel::bottom("forge3d_status").show(ctx, |ui| {
            draw_status_banner(ui, &self.status_bar, self.status_is_error);
        });

        egui::SidePanel::left("forge3d_controls")
            .default_width(360.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("Generator & Slicer Parameters");
                ui.separator();

                let mut changed = false;
                ui.label("Biological Morphology Source:");
                changed |= ui
                    .radio_value(
                        &mut self.source_kind,
                        ForgeSourceKind::GielisDiatom,
                        "Gielis (2003) Silica Diatom Frustule",
                    )
                    .changed();
                changed |= ui
                    .radio_value(
                        &mut self.source_kind,
                        ForgeSourceKind::Radiolarian,
                        "Haeckel Radiolarian Endoskeleton",
                    )
                    .changed();
                changed |= ui
                    .radio_value(
                        &mut self.source_kind,
                        ForgeSourceKind::PdbMolecule,
                        "RCSB PDB Macromolecule (.pdb)",
                    )
                    .changed();
                changed |= ui
                    .radio_value(
                        &mut self.source_kind,
                        ForgeSourceKind::CrvVoxelCsv,
                        "CRV Sim 3D Voxel Lattice (.csv)",
                    )
                    .changed();

                ui.separator();
                match self.source_kind {
                    ForgeSourceKind::PdbMolecule => {
                        ui.label("PDB File (blank = built-in Alpha-Helix):");
                        ui.horizontal(|ui| {
                            ui.add(egui::TextEdit::singleline(&mut self.file_path_input).desired_width(210.0));
                            if ui.button("Browse...").clicked() {
                                if let Some(p) = pick_open_file("Open PDB File", "PDB Structure", &["pdb", "ent"]) {
                                    self.file_path_input = p.to_string_lossy().to_string();
                                    changed = true;
                                }
                            }
                        });
                        changed |= ui
                            .radio_value(
                                &mut self.pdb_style,
                                PdbPrintStyle::GaussianSurface,
                                "Solvent Gaussian Surface",
                            )
                            .changed();
                        changed |= ui
                            .radio_value(
                                &mut self.pdb_style,
                                PdbPrintStyle::BackboneTube,
                                "C-Alpha Backbone Worm Tube",
                            )
                            .changed();
                    }
                    ForgeSourceKind::GielisDiatom => {
                        let lobes_resp = ui.add(
                            egui::Slider::new(&mut self.diatom_lobes, 3.0..=12.0)
                                .step_by(1.0)
                                .text("Symmetry Lobes (m)"),
                        );
                        Self::track_slider_rebuild(&lobes_resp, &mut changed, &mut self.mesh_dirty);

                        let pore_resp = ui
                            .add(
                                egui::Slider::new(&mut self.diatom_pore_spacing, 0.18..=0.42)
                                    .text("Areolae Pore Pitch"),
                            )
                            .on_hover_text(
                                "Center-to-center spacing of hexagonal silica areolae pores across the diatom valve (dimensionless units; smaller values yield denser perforation).",
                            );
                        Self::track_slider_rebuild(&pore_resp, &mut changed, &mut self.mesh_dirty);
                    }
                    ForgeSourceKind::Radiolarian => {
                        ui.label("Icosahedral 12-spicule perforated gyroid cortex.");
                    }
                    ForgeSourceKind::CrvVoxelCsv => {
                        ui.label("Voxel CSV Path (`x,y,z,density`):");
                        ui.horizontal(|ui| {
                            ui.add(egui::TextEdit::singleline(&mut self.file_path_input).desired_width(210.0));
                            if ui.button("Browse...").clicked() {
                                if let Some(p) = pick_open_file("Open Voxel CSV", "CSV", &["csv", "tsv"]) {
                                    self.file_path_input = p.to_string_lossy().to_string();
                                    changed = true;
                                }
                            }
                        });
                    }
                }

                let res_resp = ui
                    .add(egui::Slider::new(&mut self.grid_res, 24..=72).text("Voxel Grid Resolution"))
                    .on_hover_text(
                        "3D scalar-field sampling density along the longest axis before Surface Nets extraction (24-72 voxels). Higher values capture finer detail at the cost of slower rebuilds.",
                    );
                Self::track_slider_rebuild(&res_resp, &mut changed, &mut self.mesh_dirty);

                let bed_resp = ui.add(
                    egui::Slider::new(&mut self.target_size_mm, 30.0..=200.0)
                        .text("Max Bed Dimension (mm)"),
                );
                Self::track_slider_rebuild(&bed_resp, &mut changed, &mut self.mesh_dirty);

                if self.mesh_dirty && !ctx.input(|i| i.pointer.any_down()) {
                    changed = true;
                }

                ui.horizontal(|ui| {
                    ui.checkbox(&mut self.draw_wireframe, "Wireframe");
                    ui.checkbox(&mut self.use_perspective, "Perspective 3D");
                });

                ui.horizontal(|ui| {
                    if changed || ui.button("Rebuild Watertight Mesh").clicked() {
                        self.mesh_dirty = false;
                        self.rebuild_mesh();
                    }
                    if ui
                        .button("Reset Camera")
                        .on_hover_text("Reset orbit, pan, and zoom back to default view (R)")
                        .clicked()
                    {
                        self.reset_camera();
                    }
                });

                ui.separator();
                ui.heading("3D Printability & Manifold Audit");
                let r = &self.report;
                let badge = if r.is_watertight_manifold {
                    (Color32::GREEN, "[OK] WATERTIGHT 2-MANIFOLD (Slicer Ready)")
                } else {
                    (Color32::YELLOW, "! Open / Non-manifold Edges Detected")
                };
                ui.colored_label(badge.0, badge.1);
                ui.monospace(format!("Vertices / Faces : {} / {}", r.vertex_count, r.triangle_count));
                ui.monospace(format!("Euler (V-E+F)    : {}", r.euler_characteristic))
                    .on_hover_text(
                        "Topological Euler characteristic (Vertices - Edges + Faces). Closed sphere-like genus-0 shells satisfy V - E + F = 2; genus-g porous meshes satisfy 2 - 2g.",
                    );
                ui.monospace(format!(
                    "Size (X×Y×Z mm)  : {:.1} × {:.1} × {:.1}",
                    r.dimensions_mm[0], r.dimensions_mm[1], r.dimensions_mm[2]
                ));
                ui.monospace(format!("Volume           : {:.1} mm³", r.signed_volume_mm3));
                ui.monospace(format!("Surface Area     : {:.1} mm²", r.surface_area_mm2));
                ui.monospace(format!("Est. FDM PLA     : {:.1} g (20% infill)", r.est_fdm_pla_grams));
                ui.monospace(format!("Est. SLA Resin   : {:.1} g (solid)", r.est_sla_resin_grams));

                ui.separator();
                ui.horizontal(|ui| {
                    if ui.button("Export .stl").clicked() {
                        let def = format!("{}.stl", self.mesh.name);
                        if let Some(path) = pick_save_file("Export Binary STL", &def, "Binary STL", &["stl"]) {
                            write_bytes_reported(
                                &path,
                                &self.mesh.to_binary_stl(),
                                &mut self.status_bar,
                                &mut self.status_is_error,
                                "Saved Binary STL",
                            );
                        }
                    }
                    if ui.button("Export .3mf").clicked() {
                        let def = format!("{}.3mf", self.mesh.name);
                        if let Some(path) = pick_save_file("Export 3MF Package", &def, "3MF Archive", &["3mf"]) {
                            write_bytes_reported(
                                &path,
                                &self.mesh.to_3mf_package(),
                                &mut self.status_bar,
                                &mut self.status_is_error,
                                "Saved 3MF package",
                            );
                        }
                    }
                    if ui.button("Export .obj").clicked() {
                        let def = format!("{}.obj", self.mesh.name);
                        if let Some(path) = pick_save_file("Export Wavefront OBJ", &def, "Wavefront OBJ", &["obj"]) {
                            write_text_reported(
                                &path,
                                &self.mesh.to_wavefront_obj(),
                                &mut self.status_bar,
                                &mut self.status_is_error,
                                "Saved Wavefront OBJ",
                            );
                        }
                    }
                });
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            self.draw_3d_viewport(ui);
        });
    }
}

impl Forge3dStudioApp {
    fn normalize3(v: [f32; 3]) -> [f32; 3] {
        let len = (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]).sqrt().max(1e-9);
        [v[0] / len, v[1] / len, v[2] / len]
    }

    fn draw_3d_viewport(&mut self, ui: &mut egui::Ui) {
        let avail = ui.available_size();
        let (resp, painter) = ui.allocate_painter(avail, Sense::click_and_drag());
        let rect = resp.rect;
        painter.rect_filled(rect, 0.0, Color32::from_rgb(15, 23, 42));
        let clipped = painter.with_clip_rect(rect);

        if resp.dragged_by(egui::PointerButton::Middle) {
            self.pan_offset += resp.drag_delta();
        } else if resp.dragged_by(egui::PointerButton::Primary) {
            let d = resp.drag_delta();
            self.yaw_rad += d.x * 0.012;
            self.pitch_rad = (self.pitch_rad + d.y * 0.012).clamp(-1.45, 1.45);
        }

        let scroll = ui.input(|i| i.smooth_scroll_delta.y);
        if scroll.abs() > 0.1 {
            self.zoom = (self.zoom * (1.0 + scroll * 0.0015)).clamp(0.35, 3.5);
        }

        let (sy, cy) = self.yaw_rad.sin_cos();
        let (sp, cp) = self.pitch_rad.sin_cos();

        let rot_dir = |v: [f32; 3]| -> [f32; 3] {
            let x1 = v[0] * cy - v[1] * sy;
            let y1 = v[0] * sy + v[1] * cy;
            let y2 = y1 * cp - v[2] * sp;
            let z2 = y1 * sp + v[2] * cp;
            [x1, y2, z2]
        };

        if !self.mesh.vertices.is_empty() && !self.mesh.triangles.is_empty() {
            let (bmin, bmax) = self.mesh.bounds();
            let center3 = [
                0.5 * (bmin[0] + bmax[0]),
                0.5 * (bmin[1] + bmax[1]),
                0.5 * (bmin[2] + bmax[2]),
            ];
            let max_dim = self
                .report
                .dimensions_mm
                .iter()
                .copied()
                .fold(1.0_f32, f32::max);
            let px_scale = (rect.width().min(rect.height()) * 0.56 * self.zoom) / max_dim;
            let screen_center = rect.center() + self.pan_offset;
            let cam_dist = (max_dim * 2.6).max(1.0);
            let use_persp = self.use_perspective;

            let rot = |p: [f32; 3]| -> [f32; 3] {
                rot_dir([p[0] - center3[0], p[1] - center3[1], p[2] - center3[2]])
            };
            let proj = |v: [f32; 3]| -> Pos2 {
                let factor = if use_persp {
                    cam_dist / (cam_dist - v[2]).max(max_dim * 0.25)
                } else {
                    1.0
                };
                Pos2::new(
                    screen_center.x + v[0] * px_scale * factor,
                    screen_center.y - v[1] * px_scale * factor,
                )
            };

            let draw_build_plate = |alpha: u8| {
                let half_bed = (max_dim * 0.65).max(self.target_size_mm * 0.52).max(20.0);
                let step = if half_bed > 85.0 { 20.0_f32 } else { 10.0_f32 };
                let n_steps = (half_bed / step).ceil() as i32;
                let extent = (n_steps as f32) * step;
                let z_bed = bmin[2];

                let minor_col = Color32::from_rgba_unmultiplied(40, 54, 78, alpha);
                let major_col = Color32::from_rgba_unmultiplied(88, 108, 138, alpha);
                let x_axis_col = Color32::from_rgba_unmultiplied(239, 68, 68, alpha);
                let y_axis_col = Color32::from_rgba_unmultiplied(34, 197, 94, alpha);

                for i in -n_steps..=n_steps {
                    let coord = (i as f32) * step;
                    let is_border = i.abs() == n_steps;
                    let col_x = if i == 0 {
                        x_axis_col
                    } else if is_border {
                        major_col
                    } else {
                        minor_col
                    };
                    let width_x = if i == 0 || is_border { 1.25_f32 } else { 0.85_f32 };
                    let a_x = proj(rot([-extent, coord, z_bed]));
                    let b_x = proj(rot([extent, coord, z_bed]));
                    clipped.line_segment([a_x, b_x], Stroke::new(width_x, col_x));

                    let col_y = if i == 0 {
                        y_axis_col
                    } else if is_border {
                        major_col
                    } else {
                        minor_col
                    };
                    let width_y = if i == 0 || is_border { 1.25_f32 } else { 0.85_f32 };
                    let a_y = proj(rot([coord, -extent, z_bed]));
                    let b_y = proj(rot([coord, extent, z_bed]));
                    clipped.line_segment([a_y, b_y], Stroke::new(width_y, col_y));
                }
            };

            let bed_above = rot_dir([0.0, 0.0, 1.0])[2] >= 0.0;
            if bed_above {
                draw_build_plate(215);
            }

            let transformed: Vec<[f32; 3]> = self.mesh.vertices.iter().map(|&v| rot(v)).collect();

            // Blinn-Phong specular + rim fill + depth cueing shading
            let key_light = Self::normalize3([0.42, 0.56, 0.71]);
            let rim_light = Self::normalize3([-0.68, 0.24, 0.22]);
            let mut draw_list: Vec<(f32, [Pos2; 3], Color32)> =
                Vec::with_capacity(self.mesh.triangles.len() / 2);

            for &[i0, i1, i2] in &self.mesh.triangles {
                let p0 = transformed[i0 as usize];
                let p1 = transformed[i1 as usize];
                let p2 = transformed[i2 as usize];

                let e1 = [p1[0] - p0[0], p1[1] - p0[1], p1[2] - p0[2]];
                let e2 = [p2[0] - p0[0], p2[1] - p0[1], p2[2] - p0[2]];
                let nx = e1[1] * e2[2] - e1[2] * e2[1];
                let ny = e1[2] * e2[0] - e1[0] * e2[2];
                let nz = e1[0] * e2[1] - e1[1] * e2[0];
                let c = [
                    (p0[0] + p1[0] + p2[0]) * (1.0 / 3.0),
                    (p0[1] + p1[1] + p2[1]) * (1.0 / 3.0),
                    (p0[2] + p1[2] + p2[2]) * (1.0 / 3.0),
                ];

                let view_dir = if use_persp {
                    Self::normalize3([-c[0], -c[1], cam_dist - c[2]])
                } else {
                    [0.0, 0.0, 1.0]
                };
                let face_dot = nx * view_dir[0] + ny * view_dir[1] + nz * view_dir[2];
                if face_dot <= 0.0 {
                    continue; // Back-facing in view space
                }

                let n = Self::normalize3([nx, ny, nz]);
                let ndotl = (n[0] * key_light[0] + n[1] * key_light[1] + n[2] * key_light[2]).max(0.0);
                let half_v = Self::normalize3([
                    key_light[0] + view_dir[0],
                    key_light[1] + view_dir[1],
                    key_light[2] + view_dir[2],
                ]);
                let ndoth = (n[0] * half_v[0] + n[1] * half_v[1] + n[2] * half_v[2]).max(0.0);
                let spec = ndoth.powf(36.0) * 0.45;

                let ndotv = (n[0] * view_dir[0] + n[1] * view_dir[1] + n[2] * view_dir[2]).clamp(0.0, 1.0);
                let fresnel_rim = (1.0 - ndotv).powf(2.4) * 0.26;
                let fill_rim = (n[0] * rim_light[0] + n[1] * rim_light[1] + n[2] * rim_light[2]).max(0.0) * 0.16;

                let depth_t = ((c[2] / (0.65 * max_dim)) * 0.5 + 0.5).clamp(0.0, 1.0);
                let depth_cue = 0.68 + 0.32 * depth_t;

                let diff_term = 0.16 + 0.78 * ndotl;
                let r = ((0.22 * diff_term + 0.45 * (fresnel_rim + fill_rim) + 0.96 * spec)
                    * depth_cue
                    * 255.0)
                    .clamp(0.0, 255.0) as u8;
                let g = ((0.72 * diff_term + 0.78 * (fresnel_rim + fill_rim) + 0.98 * spec)
                    * depth_cue
                    * 255.0)
                    .clamp(0.0, 255.0) as u8;
                let b = ((0.90 * diff_term + 1.00 * (fresnel_rim + fill_rim) + 1.00 * spec)
                    * depth_cue
                    * 255.0)
                    .clamp(0.0, 255.0) as u8;
                let col = Color32::from_rgb(r, g, b);

                let sort_key = if use_persp {
                    let dx = c[0];
                    let dy = c[1];
                    let dz = cam_dist - c[2];
                    -(dx * dx + dy * dy + dz * dz)
                } else {
                    c[2]
                };
                draw_list.push((sort_key, [proj(p0), proj(p1), proj(p2)], col));
            }

            draw_list.sort_by(|a, b| a.0.total_cmp(&b.0));

            let stroke = if self.draw_wireframe {
                Stroke::new(0.6_f32, Color32::from_rgb(15, 23, 42))
            } else {
                Stroke::NONE
            };

            for (_, pts, fill) in draw_list {
                clipped.add(Shape::convex_polygon(pts.to_vec(), fill, stroke));
            }

            if !bed_above {
                draw_build_plate(90);
            }
        }

        // Bottom-left 3D XYZ axis orientation triad HUD
        let triad_origin = Pos2::new(rect.left() + 46.0, rect.bottom() - 46.0);
        let triad_len = 26.0_f32;
        let mut axes = [
            ("X", rot_dir([1.0, 0.0, 0.0]), Color32::from_rgb(248, 113, 113)),
            ("Y", rot_dir([0.0, 1.0, 0.0]), Color32::from_rgb(74, 222, 128)),
            ("Z", rot_dir([0.0, 0.0, 1.0]), Color32::from_rgb(96, 165, 250)),
        ];
        axes.sort_by(|a, b| a.1[2].total_cmp(&b.1[2]));
        clipped.circle_filled(triad_origin, 3.0, Color32::from_rgb(203, 213, 225));
        for (label, dir, col) in axes {
            let tip = Pos2::new(
                triad_origin.x + dir[0] * triad_len,
                triad_origin.y - dir[1] * triad_len,
            );
            clipped.line_segment([triad_origin, tip], Stroke::new(2.0_f32, col));
            clipped.circle_filled(tip, 2.5, col);
            let lbl_pos = Pos2::new(
                triad_origin.x + dir[0] * (triad_len + 10.0),
                triad_origin.y - dir[1] * (triad_len + 10.0),
            );
            clipped.text(
                lbl_pos,
                egui::Align2::CENTER_CENTER,
                label,
                egui::FontId::monospace(11.0),
                col,
            );
        }

        let proj_lbl = if self.use_perspective { "Perspective" } else { "Ortho" };
        clipped.text(
            Pos2::new(rect.left() + 12.0, rect.top() + 10.0),
            egui::Align2::LEFT_TOP,
            format!(
                "{proj_lbl} | Zoom {:.0}% | LMB Orbit | MMB Pan | Scroll Zoom | R Reset",
                self.zoom * 100.0
            ),
            egui::FontId::monospace(11.0),
            Color32::from_rgb(148, 163, 184),
        );
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycles_all_forge_source_kinds_headlessly() {
        let ctx = egui::Context::default();
        let mut app = Forge3dStudioApp::new(None);
        for kind in [
            ForgeSourceKind::GielisDiatom,
            ForgeSourceKind::Radiolarian,
            ForgeSourceKind::PdbMolecule,
            ForgeSourceKind::CrvVoxelCsv,
        ] {
            app.source_kind = kind;
            app.rebuild_mesh();
            for frame in 0..2 {
                let mut raw = egui::RawInput::default();
                raw.screen_rect = Some(egui::Rect::from_min_size(Pos2::ZERO, egui::vec2(1280.0, 800.0)));
                raw.time = Some(frame as f64 / 60.0);
                let out = ctx.run(raw, |c| app.draw_panels(c));
                assert!(!out.shapes.is_empty());
            }
        }
    }

    #[test]
    fn r_hotkey_resets_camera_and_perspective_render_succeeds() {
        let ctx = egui::Context::default();
        let mut app = Forge3dStudioApp::new(None);
        app.yaw_rad = 1.9;
        app.pitch_rad = -0.7;
        app.zoom = 2.4;
        app.pan_offset = Vec2::new(45.0, -30.0);
        app.use_perspective = true;

        let mut raw = egui::RawInput::default();
        raw.screen_rect = Some(egui::Rect::from_min_size(Pos2::ZERO, egui::vec2(1280.0, 800.0)));
        raw.events.push(egui::Event::Key {
            key: egui::Key::R,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        });
        let out = ctx.run(raw, |c| app.draw_panels(c));
        assert!(!out.shapes.is_empty());
        assert!((app.yaw_rad - 0.65).abs() < 1e-6);
        assert!((app.pitch_rad - 0.52).abs() < 1e-6);
        assert!((app.zoom - 1.0).abs() < 1e-6);
        assert_eq!(app.pan_offset, Vec2::ZERO);
    }
}
