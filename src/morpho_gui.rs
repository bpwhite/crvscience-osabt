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

//! Pure-Rust cross-platform `eframe` / `egui` desktop studio for Tool #8: `osabt morpho`.
//! Integrates automated Petri CFU counting (with manual click curation), Excess-Green leaf
//! herbivory thresholding, and interactive `tpsDig`-style 2D homologous landmarking + live GPA.

use std::path::{Path, PathBuf};

use eframe::egui;
use egui::{Color32, ColorImage, Pos2, Sense, Shape, Stroke, TextureHandle, TextureOptions, Vec2};

use crate::gui_io::{
    draw_status_banner, pick_open_file, pick_save_file, read_bytes_reported, read_text_reported,
    write_text_reported,
};
use crate::morpho_gpa::{
    demo_salmonid_landmarks_tps, parse_tps_landmarks, render_gpa_morphospace_svg,
    run_generalized_procrustes_2d, write_tps_landmarks, LandmarkSpecimen2D,
};
use crate::morpho_vision::{
    analyze_leaf_herbivory, count_petri_colonies, estimate_cfu_per_ml, render_petri_cfu_svg,
    DetectedColony, PetriCfuReport, RgbImage,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum MorphoMode {
    PetriCfu,
    LeafHerbivory,
    LandmarksGpa,
}

/// Interactive `egui` application state for `OSABT Morpho Studio`.
pub struct MorphoStudioApp {
    mode: MorphoMode,
    path_input: String,
    petri_img: RgbImage,
    leaf_img: RgbImage,
    texture: Option<TextureHandle>,
    texture_dirty: bool,
    dark_colonies: bool,
    roi_fraction: f64,
    manual_thresh_on: bool,
    manual_thresh: u8,
    plated_vol_ml: f64,
    dilution_exp: i32,
    cfu_report: PetriCfuReport,
    px_per_mm: f64,
    exg_override_on: bool,
    exg_thresh: i32,
    specimens: Vec<LandmarkSpecimen2D>,
    active_specimen_idx: usize,
    status_bar: String,
    status_is_error: bool,
    about_open: bool,
}

impl MorphoStudioApp {
    pub fn new(initial_file: Option<PathBuf>) -> Self {
        let petri_img = RgbImage::demo_petri_cfu_plate();
        let leaf_img = RgbImage::demo_herbivory_leaf();
        let cfu_report = count_petri_colonies(&petri_img, false, 0.84, 8, None);
        let specimens = parse_tps_landmarks(demo_salmonid_landmarks_tps()).unwrap_or_default();

        let mut app = Self {
            mode: MorphoMode::PetriCfu,
            path_input: initial_file
                .as_ref()
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_default(),
            petri_img,
            leaf_img,
            texture: None,
            texture_dirty: true,
            dark_colonies: false,
            roi_fraction: 0.84,
            manual_thresh_on: false,
            manual_thresh: 80,
            plated_vol_ml: 0.10,
            dilution_exp: -5,
            cfu_report,
            px_per_mm: 4.0,
            exg_override_on: false,
            exg_thresh: 45,
            specimens,
            active_specimen_idx: 0,
            status_bar: String::from(
                "Left-click Petri dish to add missed colony, right-click to remove false positive | Switch tabs for Leaf Herbivory & 2D GPA Landmarks.",
            ),
            status_is_error: false,
            about_open: false,
        };

        if let Some(p) = initial_file {
            app.load_from_path(&p);
        }
        app
    }

    fn load_from_path(&mut self, path: &Path) {
        let ext = path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();
        if ext == "tps" {
            if let Some(s) = read_text_reported(path, &mut self.status_bar, &mut self.status_is_error) {
                match parse_tps_landmarks(&s) {
                    Ok(sp) => {
                        self.specimens = sp;
                        self.active_specimen_idx = 0;
                        self.mode = MorphoMode::LandmarksGpa;
                        self.status_is_error = false;
                        self.status_bar = format!("Loaded {} specimens from '{}'.", self.specimens.len(), path.display());
                    }
                    Err(e) => {
                        self.status_is_error = true;
                        self.status_bar = format!("TPS load error: {e}");
                    }
                }
            }
        } else if let Some(bytes) = read_bytes_reported(path, &mut self.status_bar, &mut self.status_is_error) {
            match RgbImage::from_bytes(&bytes) {
                Ok(img) => {
                    if self.mode == MorphoMode::LeafHerbivory {
                        self.leaf_img = img;
                    } else {
                        self.petri_img = img;
                        self.recompute_cfu();
                    }
                    self.texture_dirty = true;
                    self.status_is_error = false;
                    self.status_bar = format!("Loaded image '{}'.", path.display());
                }
                Err(e) => {
                    self.status_is_error = true;
                    self.status_bar = format!("Image parse error: {e}");
                }
            }
        }
    }

    fn recompute_cfu(&mut self) {
        let ov = if self.manual_thresh_on {
            Some(self.manual_thresh)
        } else {
            None
        };
        self.cfu_report = count_petri_colonies(
            &self.petri_img,
            self.dark_colonies,
            self.roi_fraction,
            8,
            ov,
        );
    }

    fn ensure_texture(&mut self, ctx: &egui::Context) {
        if !self.texture_dirty && self.texture.is_some() {
            return;
        }
        let src = match self.mode {
            MorphoMode::LeafHerbivory => &self.leaf_img,
            _ => &self.petri_img,
        };
        let color_img = ColorImage::from_rgb([src.width, src.height], &src.rgb);
        self.texture = Some(ctx.load_texture("morpho_canvas_tex", color_img, TextureOptions::LINEAR));
        self.texture_dirty = false;
    }
}

/// Launches the native cross-platform `OSABT Morpho Studio` desktop GUI.
pub fn run_morpho_gui(initial_file: Option<PathBuf>) -> Result<(), String> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1200.0, 780.0])
            .with_icon(crate::brand::osabt_icon_data())
            .with_title("OSABT Morpho — Colony Counter, Leaf Herbivory & 2D GPA Studio (GPL-3.0)"),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-morpho-studio",
        options,
        Box::new(move |_cc| Ok(Box::new(MorphoStudioApp::new(initial_file)))),
    )
    .map_err(|e| format!("Failed starting Morpho desktop window: {e}"))
}

impl eframe::App for MorphoStudioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        crate::brand::draw_standalone_footer(ctx, &mut self.about_open);
        self.draw_panels(ctx);
    }
}

impl MorphoStudioApp {
    /// Renders all Morpho Studio panels onto `ctx`.
    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("morpho_top").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                if ui
                    .selectable_label(self.mode == MorphoMode::PetriCfu, "Petri CFU Counter")
                    .clicked()
                {
                    self.mode = MorphoMode::PetriCfu;
                    self.texture_dirty = true;
                }
                if ui
                    .selectable_label(
                        self.mode == MorphoMode::LeafHerbivory,
                        "Leaf Area & Herbivory",
                    )
                    .clicked()
                {
                    self.mode = MorphoMode::LeafHerbivory;
                    self.texture_dirty = true;
                }
                if ui
                    .selectable_label(
                        self.mode == MorphoMode::LandmarksGpa,
                        "2D Landmarks & Procrustes GPA",
                    )
                    .clicked()
                {
                    self.mode = MorphoMode::LandmarksGpa;
                }
                ui.separator();
                ui.label("Load (.bmp/.ppm/.tps):");
                ui.add(egui::TextEdit::singleline(&mut self.path_input).desired_width(180.0));
                if ui.button("Browse...").clicked() {
                    if let Some(p) = pick_open_file("Open Morpho Image or TPS", "Image / TPS", &["bmp", "ppm", "tps"]) {
                        self.path_input = p.to_string_lossy().to_string();
                        self.load_from_path(&p);
                    }
                }
                if ui.button("Open").clicked() && !self.path_input.trim().is_empty() {
                    let p = PathBuf::from(self.path_input.trim());
                    self.load_from_path(&p);
                }
            });
        });

        egui::TopBottomPanel::bottom("morpho_bot").show(ctx, |ui| {
            draw_status_banner(ui, &self.status_bar, self.status_is_error);
        });

        egui::SidePanel::left("morpho_side")
            .default_width(370.0)
            .show(ctx, |ui| match self.mode {
                MorphoMode::PetriCfu => self.draw_cfu_sidebar(ui),
                MorphoMode::LeafHerbivory => self.draw_leaf_sidebar(ui),
                MorphoMode::LandmarksGpa => self.draw_gpa_sidebar(ui),
            });

        self.ensure_texture(ctx);
        egui::CentralPanel::default().show(ctx, |ui| match self.mode {
            MorphoMode::PetriCfu => self.draw_cfu_canvas(ui),
            MorphoMode::LeafHerbivory => self.draw_leaf_canvas(ui),
            MorphoMode::LandmarksGpa => self.draw_gpa_canvas(ui),
        });
    }
}

impl MorphoStudioApp {
    fn draw_cfu_sidebar(&mut self, ui: &mut egui::Ui) {
        ui.heading("Automated Petri Dish CFU Counter");
        let mut changed = false;
        changed |= ui
            .checkbox(&mut self.dark_colonies, "Dark Colonies on Bright-Field Agar")
            .changed();
        changed |= ui
            .add(egui::Slider::new(&mut self.roi_fraction, 0.50..=0.98).text("Inner Dish ROI"))
            .changed();
        changed |= ui
            .checkbox(&mut self.manual_thresh_on, "Override Otsu Threshold")
            .changed();
        if self.manual_thresh_on {
            changed |= ui
                .add(egui::Slider::new(&mut self.manual_thresh, 5..=250).text("Threshold"))
                .changed();
        }
        if changed {
            self.recompute_cfu();
        }

        ui.separator();
        ui.add(egui::Slider::new(&mut self.plated_vol_ml, 0.01..=1.0).text("Plated Volume (mL)"));
        ui.add(egui::Slider::new(&mut self.dilution_exp, -9..=0).text("Dilution 10^k"));

        let dil = 10_f64.powi(self.dilution_exp);
        let r = &self.cfu_report;
        ui.separator();
        ui.strong(format!("Detected Total CFUs : {}", r.total_cfu_count));
        ui.monospace(format!("Isolated Singletons : {}", r.isolated_colonies));
        ui.monospace(format!("Split from Clusters : {}", r.split_cluster_colonies));
        ui.monospace(format!("Otsu/Active Thresh  : {}", r.threshold_used));
        ui.monospace(format!("Mean Diameter       : {:.1} px", r.mean_colony_diameter_px));

        if let Ok((mle, lo, hi)) =
            estimate_cfu_per_ml(&[(r.total_cfu_count, self.plated_vol_ml, dil)])
        {
            ui.colored_label(
                Color32::LIGHT_GREEN,
                format!("Est. Titer: {mle:.2e} CFU/mL\n95% CI: [{lo:.2e}, {hi:.2e}]"),
            );
        }

        if ui.button("Export Annotated CFU .svg").clicked() {
            if let Some(path) = pick_save_file("Export Annotated CFU SVG", "petri_cfu_annotated.svg", "SVG Image", &["svg"]) {
                let svg = render_petri_cfu_svg(&self.petri_img, &self.cfu_report);
                write_text_reported(
                    &path,
                    &svg,
                    &mut self.status_bar,
                    &mut self.status_is_error,
                    "Saved Annotated CFU SVG",
                );
            }
        }
    }

    fn draw_cfu_canvas(&mut self, ui: &mut egui::Ui) {
        let Some(tex) = &self.texture else {
            return;
        };
        let scale = 2.2_f32;
        let disp = Vec2::new(
            (self.petri_img.width as f32) * scale,
            (self.petri_img.height as f32) * scale,
        );
        let (resp, painter) = ui.allocate_painter(disp, Sense::click());
        let rect = resp.rect;
        painter.image(
            tex.id(),
            rect,
            egui::Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
            Color32::WHITE,
        );

        // Draw circular ROI mask & detected CFU circles
        let roi_r = (self.petri_img.width.min(self.petri_img.height) as f32)
            * 0.5
            * (self.roi_fraction as f32)
            * scale;
        painter.circle_stroke(
            rect.center(),
            roi_r,
            Stroke::new(1.4_f32, Color32::from_rgb(56, 189, 248)),
        );

        for c in &self.cfu_report.colonies {
            let p = Pos2::new(
                rect.left() + (c.x_px as f32) * scale,
                rect.top() + (c.y_px as f32) * scale,
            );
            let col = if c.split_from_cluster {
                Color32::from_rgb(251, 191, 36)
            } else {
                Color32::from_rgb(34, 197, 94)
            };
            painter.circle_stroke(
                p,
                (c.equivalent_radius_px as f32 * scale).max(4.0),
                Stroke::new(1.8_f32, col),
            );
        }

        // Manual click curation: primary click adds colony, secondary click removes nearest
        if resp.clicked() {
            if let Some(pos) = resp.interact_pointer_pos() {
                let ix = ((pos.x - rect.left()) / scale) as f64;
                let iy = ((pos.y - rect.top()) / scale) as f64;
                self.cfu_report.colonies.push(DetectedColony {
                    x_px: ix,
                    y_px: iy,
                    equivalent_radius_px: 5.0,
                    area_px: 78,
                    split_from_cluster: false,
                });
                self.cfu_report.total_cfu_count = self.cfu_report.colonies.len();
            }
        } else if resp.secondary_clicked() {
            if let Some(pos) = resp.interact_pointer_pos() {
                let ix = ((pos.x - rect.left()) / scale) as f64;
                let iy = ((pos.y - rect.top()) / scale) as f64;
                if let Some((best_i, _)) = self
                    .cfu_report
                    .colonies
                    .iter()
                    .enumerate()
                    .map(|(i, c)| (i, (c.x_px - ix).powi(2) + (c.y_px - iy).powi(2)))
                    .min_by(|a, b| a.1.partial_cmp(&b.1).unwrap_or(std::cmp::Ordering::Equal))
                {
                    self.cfu_report.colonies.remove(best_i);
                    self.cfu_report.total_cfu_count = self.cfu_report.colonies.len();
                }
            }
        }
    }

    fn draw_leaf_sidebar(&mut self, ui: &mut egui::Ui) {
        ui.heading("Excess-Green Leaf Area & Herbivory");
        ui.add(egui::Slider::new(&mut self.px_per_mm, 0.5..=20.0).text("Scale (px / mm)"));
        ui.checkbox(&mut self.exg_override_on, "Manual ExG (2G-R-B) Threshold");
        if self.exg_override_on {
            ui.add(egui::Slider::new(&mut self.exg_thresh, -40..=180).text("ExG Cutoff"));
        }

        let ov = if self.exg_override_on {
            Some(self.exg_thresh)
        } else {
            None
        };
        if let Ok(rep) = analyze_leaf_herbivory(&self.leaf_img, self.px_per_mm, ov) {
            ui.separator();
            ui.monospace(format!("Intact Green Lamina  : {:.2} cm²", rep.lamina_area_cm2));
            ui.monospace(format!(
                "Interior Shot-Holes  : {} holes ({:.3} cm²)",
                rep.enclosed_hole_count, rep.enclosed_herbivory_cm2
            ));
            ui.monospace(format!("Convex Hull Envelope : {:.2} cm²", rep.convex_envelope_cm2));
            ui.colored_label(
                Color32::YELLOW,
                format!("Interior Herbivory   : {:.2}%", rep.interior_herbivory_pct),
            );
            ui.colored_label(
                Color32::LIGHT_RED,
                format!(
                    "Total Envelope Loss  : {:.2}% (incl. marginal bites)",
                    rep.total_defoliation_envelope_pct
                ),
            );
        }
    }

    fn draw_leaf_canvas(&mut self, ui: &mut egui::Ui) {
        let Some(tex) = &self.texture else {
            return;
        };
        let scale = 2.2_f32;
        let disp = Vec2::new(
            (self.leaf_img.width as f32) * scale,
            (self.leaf_img.height as f32) * scale,
        );
        let (resp, painter) = ui.allocate_painter(disp, Sense::hover());
        painter.image(
            tex.id(),
            resp.rect,
            egui::Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
            Color32::WHITE,
        );
    }

    fn draw_gpa_sidebar(&mut self, ui: &mut egui::Ui) {
        ui.heading("2D Homologous Landmarks & GPA");
        ui.label(format!("Loaded Specimens: {}", self.specimens.len()));
        for (idx, sp) in self.specimens.iter().enumerate() {
            ui.selectable_value(
                &mut self.active_specimen_idx,
                idx,
                format!("{} (LM={})", sp.specimen_id, sp.landmarks.len()),
            );
        }

        ui.separator();
        if let Ok(gpa) = run_generalized_procrustes_2d(&self.specimens) {
            ui.monospace(format!("PC1 Variance Explained : {:.1}%", gpa.pc1_variance_pct));
            ui.monospace(format!("PC2 Variance Explained : {:.1}%", gpa.pc2_variance_pct));
            if let Some(&cs) = gpa.centroid_sizes.get(self.active_specimen_idx) {
                ui.monospace(format!("Active Centroid Size   : {cs:.3}"));
            }

            ui.horizontal(|ui| {
                if ui.button("Export .tps").clicked() {
                    if let Some(path) = pick_save_file("Export Landmarks TPS", "osabt_landmarks.tps", "TPS Landmarks", &["tps"]) {
                        write_text_reported(
                            &path,
                            &write_tps_landmarks(&self.specimens),
                            &mut self.status_bar,
                            &mut self.status_is_error,
                            "Saved TPS Landmarks",
                        );
                    }
                }
                if ui.button("Export GPA+PCA .svg").clicked() {
                    if let Some(path) = pick_save_file("Export GPA Morphospace SVG", "osabt_gpa_morphospace.svg", "SVG Image", &["svg"]) {
                        write_text_reported(
                            &path,
                            &render_gpa_morphospace_svg(&gpa),
                            &mut self.status_bar,
                            &mut self.status_is_error,
                            "Saved GPA+PCA SVG",
                        );
                    }
                }
            });
        }
    }

    fn draw_gpa_canvas(&mut self, ui: &mut egui::Ui) {
        let avail = ui.available_size();
        let (resp, painter) = ui.allocate_painter(avail, Sense::click());
        let rect = resp.rect;
        painter.rect_filled(rect, 4.0, Color32::from_rgb(15, 23, 42));

        if let Ok(gpa) = run_generalized_procrustes_2d(&self.specimens) {
            let center = rect.center();
            let radius = rect.width().min(rect.height()) * 0.65;
            let map_pt =
                |p: [f64; 2]| Pos2::new(center.x + (p[0] as f32) * radius, center.y - (p[1] as f32) * radius);

            for (s_idx, shape) in gpa.aligned_shapes.iter().enumerate() {
                let col = if s_idx == self.active_specimen_idx {
                    Color32::from_rgb(250, 204, 21)
                } else {
                    Color32::from_rgba_unmultiplied(125, 211, 252, 160)
                };
                for &pt in shape {
                    painter.circle_filled(map_pt(pt), 4.0, col);
                }
            }

            let mut cons_pts: Vec<Pos2> = gpa.consensus_shape.iter().map(|&p| map_pt(p)).collect();
            if let Some(&first) = cons_pts.first() {
                cons_pts.push(first);
            }
            painter.add(Shape::line(
                cons_pts,
                Stroke::new(2.2_f32, Color32::from_rgb(248, 113, 113)),
            ));
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycles_all_morpho_modes_headlessly() {
        let ctx = egui::Context::default();
        let mut app = MorphoStudioApp::new(None);
        for mode in [
            MorphoMode::PetriCfu,
            MorphoMode::LeafHerbivory,
            MorphoMode::LandmarksGpa,
        ] {
            app.mode = mode;
            app.texture_dirty = true;
            for frame in 0..2 {
                let mut raw = egui::RawInput::default();
                raw.screen_rect = Some(egui::Rect::from_min_size(Pos2::ZERO, Vec2::new(1280.0, 800.0)));
                raw.time = Some(frame as f64 / 60.0);
                let out = ctx.run(raw, |c| app.draw_panels(c));
                assert!(!out.shapes.is_empty());
            }
        }
    }
}
