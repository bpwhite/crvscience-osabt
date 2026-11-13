// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Interactive `eframe`/`egui` Flow Cytometry `.fcs 3.1` Spillover Compensation
//! & `Logicle` Biaxial Gating Studio (`osabt flow`).

use std::fs;
use std::path::PathBuf;

use eframe::egui;

use crate::flow_fcs::{compensate_fcs_dataset, demo_fcs_31_bytes, parse_fcs_bytes, FcsDataset};
use crate::flow_gate::{
    evaluate_polygon_gate, evaluate_quadrant_gates, render_flow_biaxial_svg, GatePopulationStats,
    LogicleParams, LogicleScale, PolygonGate,
};

/// Launches the `eframe`/`egui` native desktop studio for `osabt flow`.
pub fn run_flow_gui(fcs_path: Option<PathBuf>) -> Result<(), String> {
    let app = FlowGuiApp::new(fcs_path)?;
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("OSABT Flow Studio (v1.8.0) — FCS 3.1 Compensation & Logicle Gating")
            .with_inner_size([1140.0, 740.0]),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-flow",
        native_options,
        Box::new(|_cc| Ok(Box::new(app))),
    )
    .map_err(|e| format!("Failed launching egui Flow Studio: {e}"))
}

struct FlowGuiApp {
    raw_ds: FcsDataset,
    comp_ds: Option<FcsDataset>,
    apply_compensation: bool,
    x_ch: usize,
    y_ch: usize,
    x_cut_disp: f64,
    y_cut_disp: f64,
    w_decades: f64,
    scale: LogicleScale,
    quad_stats: Vec<GatePopulationStats>,
    live_gate_stats: Option<GatePopulationStats>,
    status_msg: String,
}

impl FlowGuiApp {
    fn new(fcs_path: Option<PathBuf>) -> Result<Self, String> {
        let bytes = if let Some(fp) = fcs_path {
            fs::read(&fp)
                .map_err(|e| format!("Failed reading FCS file '{}': {e}", fp.display()))?
        } else {
            demo_fcs_31_bytes()
        };
        let raw_ds = parse_fcs_bytes(&bytes)?;
        let comp_ds = compensate_fcs_dataset(&raw_ds).ok().map(|(d, _)| d);
        let apply_comp = comp_ds.is_some();
        let x_ch = if raw_ds.channels.len() >= 4 { 2 } else { 0 };
        let y_ch = if raw_ds.channels.len() >= 4 { 3 } else { 1.min(raw_ds.channels.len().saturating_sub(1)) };

        let params = LogicleParams::default();
        let mut app = Self {
            raw_ds,
            comp_ds,
            apply_compensation: apply_comp,
            x_ch,
            y_ch,
            x_cut_disp: 0.38,
            y_cut_disp: 0.38,
            w_decades: params.w_decades,
            scale: LogicleScale::new(params),
            quad_stats: Vec::new(),
            live_gate_stats: None,
            status_msg: String::new(),
        };
        app.recompute();
        Ok(app)
    }

    fn active_ds(&self) -> &FcsDataset {
        if self.apply_compensation {
            if let Some(c) = &self.comp_ds {
                return c;
            }
        }
        &self.raw_ds
    }

    fn recompute(&mut self) {
        let params = LogicleParams {
            w_decades: self.w_decades,
            ..LogicleParams::default()
        };
        self.scale = LogicleScale::new(params);
        let ds = self.active_ds();
        self.quad_stats = evaluate_quadrant_gates(
            ds,
            self.x_ch,
            self.y_ch,
            self.x_cut_disp,
            self.y_cut_disp,
            &self.scale,
        )
        .unwrap_or_default();

        // Optional Scatter Live-Cell polygon gate when FSC-A (0) & SSC-A (1) exist
        if ds.channels.len() >= 2 {
            let poly = PolygonGate {
                name: String::from("Live_Singlet_Scatter"),
                vertices: vec![(0.45, 0.35), (0.95, 0.35), (0.95, 0.90), (0.45, 0.90)],
            };
            self.live_gate_stats = evaluate_polygon_gate(ds, 0, 1, &poly, &self.scale).ok();
        }
        self.status_msg = format!(
            "Loaded {} ({} events x {} channels) | Compensation: {}",
            ds.version,
            ds.num_events,
            ds.channels.len(),
            if self.apply_compensation { "ON (S^-1 applied)" } else { "RAW" }
        );
    }
}

impl eframe::App for FlowGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("flow_top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("OSABT Flow v1.8.0");
                ui.separator();
                if self.comp_ds.is_some()
                    && ui
                        .checkbox(&mut self.apply_compensation, "Apply $SPILLOVER Compensation")
                        .changed()
                {
                    self.recompute();
                }
                if ui
                    .add(egui::Slider::new(&mut self.w_decades, 0.2..=1.5).text("Logicle W"))
                    .changed()
                {
                    self.recompute();
                }
                if ui.button("Export Biaxial .svg").clicked() {
                    let svg = render_flow_biaxial_svg(
                        self.active_ds(),
                        self.x_ch,
                        self.y_ch,
                        self.x_cut_disp,
                        self.y_cut_disp,
                        &self.scale,
                        "OSABT Flow Cytometry",
                    );
                    if fs::write("flow_biaxial_gate.svg", svg).is_ok() {
                        self.status_msg = String::from("Exported flow_biaxial_gate.svg!");
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("flow_bot").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.columns(2, |cols| {
                let (resp, painter) = cols[0].allocate_painter(
                    egui::vec2(cols[0].available_width().max(480.0), 540.0),
                    egui::Sense::click_and_drag(),
                );
                let rect = resp.rect;
                painter.rect_filled(rect, 6.0, egui::Color32::from_rgb(248, 250, 252));

                if resp.dragged() || resp.clicked() {
                    if let Some(p) = resp.interact_pointer_pos() {
                        self.x_cut_disp =
                            (((p.x - rect.left()) / rect.width()) as f64).clamp(0.05, 0.95);
                        self.y_cut_disp =
                            (((rect.bottom() - p.y) / rect.height()) as f64).clamp(0.05, 0.95);
                        self.recompute();
                    }
                }

                let ds = self.active_ds();
                for ev in &ds.events {
                    let dx = self.scale.transform(ev[self.x_ch]).clamp(0.0, 1.0) as f32;
                    let dy = self.scale.transform(ev[self.y_ch]).clamp(0.0, 1.0) as f32;
                    let px = rect.left() + dx * rect.width();
                    let py = rect.bottom() - dy * rect.height();
                    let col = match (dx as f64 >= self.x_cut_disp, dy as f64 >= self.y_cut_disp) {
                        (false, true) => egui::Color32::from_rgb(147, 51, 234),
                        (true, true) => egui::Color32::from_rgb(220, 38, 38),
                        (false, false) => egui::Color32::from_rgb(37, 99, 235),
                        (true, false) => egui::Color32::from_rgb(22, 163, 74),
                    };
                    painter.circle_filled(egui::pos2(px, py), 2.5, col);
                }

                let qx = rect.left() + (self.x_cut_disp as f32) * rect.width();
                let qy = rect.bottom() - (self.y_cut_disp as f32) * rect.height();
                painter.line_segment(
                    [egui::pos2(qx, rect.top()), egui::pos2(qx, rect.bottom())],
                    egui::Stroke::new(1.5, egui::Color32::from_rgb(15, 23, 42)),
                );
                painter.line_segment(
                    [egui::pos2(rect.left(), qy), egui::pos2(rect.right(), qy)],
                    egui::Stroke::new(1.5, egui::Color32::from_rgb(15, 23, 42)),
                );

                // Right column: Channel selectors + Quadrant statistics
                cols[1].heading("Axes & Quadrant Gate Statistics");
                let ch_names: Vec<String> = self
                    .raw_ds
                    .channels
                    .iter()
                    .map(|c| format!("{} ({})", c.short_name, c.stain_name))
                    .collect();
                let mut ch_changed = false;
                cols[1].horizontal(|ui| {
                    ui.label("X-Axis:");
                    egui::ComboBox::from_id_salt("flow_x")
                        .selected_text(&ch_names[self.x_ch])
                        .show_ui(ui, |ui| {
                            for (i, nm) in ch_names.iter().enumerate() {
                                if ui.selectable_value(&mut self.x_ch, i, nm).changed() {
                                    ch_changed = true;
                                }
                            }
                        });
                });
                cols[1].horizontal(|ui| {
                    ui.label("Y-Axis:");
                    egui::ComboBox::from_id_salt("flow_y")
                        .selected_text(&ch_names[self.y_ch])
                        .show_ui(ui, |ui| {
                            for (i, nm) in ch_names.iter().enumerate() {
                                if ui.selectable_value(&mut self.y_ch, i, nm).changed() {
                                    ch_changed = true;
                                }
                            }
                        });
                });
                if ch_changed {
                    self.recompute();
                }

                cols[1].separator();
                cols[1].label("Click or drag on the scatter plot to move Quadrant crosshairs.");
                egui::Grid::new("flow_quad_grid").striped(true).show(&mut cols[1], |ui| {
                    for h in ["Gate", "Events", "% Parent", "X MFI", "Y MFI", "rCV%"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for q in &self.quad_stats {
                        ui.label(&q.gate_name);
                        ui.label(format!("{}", q.event_count));
                        ui.label(format!("{:.1}%", q.percent_of_parent));
                        ui.label(format!("{:.0}", q.x_median));
                        ui.label(format!("{:.0}", q.y_median));
                        ui.label(format!("{:.1}%", q.x_robust_cv_pct));
                        ui.end_row();
                    }
                });

                if let Some(lg) = &self.live_gate_stats {
                    cols[1].separator();
                    cols[1].strong(format!(
                        "Polygon Gate '{}': {} events ({:.1}% of total)",
                        lg.gate_name, lg.event_count, lg.percent_of_parent
                    ));
                }
            });
        });
    }
}
