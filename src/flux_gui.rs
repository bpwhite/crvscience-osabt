// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Interactive `eframe`/`egui` Leaf `FvCB` Photosynthesis & Stream Open-Channel `GPP`/`ER`
//! Metabolism Studio (`osabt flux`).

use std::fs;
use std::path::PathBuf;

use eframe::egui;

use crate::flux_leaf::{
    demo_fvcb_aci_csv, eval_fvcb_assimilation, fit_fvcb_aci_curve, parse_aci_csv, FvcbAciFit,
    LeafAciPoint,
};
use crate::flux_stream::{
    demo_stream_diel_do_csv, fit_odum_stream_metabolism, parse_stream_do_csv, render_flux_svg,
    StreamDoPoint, StreamMetabolismReport,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FluxView {
    LeafAciFvcb,
    StreamDielDo,
}

/// Launches the `eframe`/`egui` native desktop studio for `osabt flux`.
pub fn run_flux_gui(aci_path: Option<PathBuf>, stream_path: Option<PathBuf>) -> Result<(), String> {
    let app = FluxGuiApp::new(aci_path, stream_path)?;
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("OSABT Flux Studio (v2.0.0) — FvCB Leaf Photosynthesis & Odum Stream Metabolism")
            .with_inner_size([1140.0, 720.0]),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-flux",
        native_options,
        Box::new(|_cc| Ok(Box::new(app))),
    )
    .map_err(|e| format!("Failed launching egui Flux Studio: {e}"))
}

struct FluxGuiApp {
    view: FluxView,
    aci_pts: Vec<LeafAciPoint>,
    aci_fit: FvcbAciFit,
    stream_pts: Vec<StreamDoPoint>,
    stream_rep: StreamMetabolismReport,
    status_msg: String,
}

impl FluxGuiApp {
    fn new(aci_path: Option<PathBuf>, stream_path: Option<PathBuf>) -> Result<Self, String> {
        let aci_raw = if let Some(ap) = aci_path {
            fs::read_to_string(&ap)
                .map_err(|e| format!("Failed reading A-Ci CSV '{}': {e}", ap.display()))?
        } else {
            demo_fvcb_aci_csv()
        };
        let aci_pts = parse_aci_csv(&aci_raw)?;
        let aci_fit = fit_fvcb_aci_curve(&aci_pts)?;

        let stream_raw = if let Some(sp) = stream_path {
            fs::read_to_string(&sp)
                .map_err(|e| format!("Failed reading stream DO CSV '{}': {e}", sp.display()))?
        } else {
            demo_stream_diel_do_csv()
        };
        let stream_pts = parse_stream_do_csv(&stream_raw)?;
        let stream_rep = fit_odum_stream_metabolism(&stream_pts)?;

        let status_msg = format!(
            "FvCB Fit: Vcmax={:.1}, Jmax={:.1}, TPU={:.1} (R^2={:.4}) | Stream: GPP={:.2}, ER={:.2}, NEP={:+.2} g O2/m2/d",
            aci_fit.vcmax,
            aci_fit.jmax,
            aci_fit.tpu,
            aci_fit.r_squared,
            stream_rep.gpp_g_o2_m2_d,
            stream_rep.er_g_o2_m2_d,
            stream_rep.nep_g_o2_m2_d
        );

        Ok(Self {
            view: FluxView::LeafAciFvcb,
            aci_pts,
            aci_fit,
            stream_pts,
            stream_rep,
            status_msg,
        })
    }
}

impl eframe::App for FluxGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("flux_top").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("OSABT Flux v2.0.0");
                ui.separator();
                ui.selectable_value(
                    &mut self.view,
                    FluxView::LeafAciFvcb,
                    "1. FvCB C3 Leaf A-Ci Curve",
                );
                ui.selectable_value(
                    &mut self.view,
                    FluxView::StreamDielDo,
                    "2. Odum Diel DO Stream Metabolism",
                );
                ui.separator();
                if ui.button("Export Dual-Panel .svg").clicked() {
                    let svg = render_flux_svg(
                        &self.aci_pts,
                        &self.aci_fit,
                        &self.stream_pts,
                        &self.stream_rep,
                        "OSABT Ecophysiology & Stream Flux",
                    );
                    if fs::write("osabt_flux_report.svg", svg).is_ok() {
                        self.status_msg = String::from("Exported osabt_flux_report.svg!");
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("flux_bot").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });

        egui::CentralPanel::default().show(ctx, |ui| match self.view {
            FluxView::LeafAciFvcb => self.draw_leaf_panel(ui),
            FluxView::StreamDielDo => self.draw_stream_panel(ui),
        });
    }
}

impl FluxGuiApp {
    fn draw_leaf_panel(&self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            let (resp, painter) = cols[0].allocate_painter(
                egui::vec2(cols[0].available_width().max(480.0), 520.0),
                egui::Sense::hover(),
            );
            let rect = resp.rect;
            painter.rect_filled(rect, 6.0, egui::Color32::from_rgb(248, 250, 252));

            let max_ci = self.aci_pts.last().map(|p| p.ci_ubar).unwrap_or(1600.0) * 1.05;
            let min_a = -2.0_f64;
            let max_a = self
                .aci_pts
                .iter()
                .map(|p| p.a_net)
                .fold(20.0_f64, f64::max)
                * 1.15;

            let map_pt = |ci: f64, a: f64| -> egui::Pos2 {
                let fx = (ci / max_ci).clamp(0.0, 1.0) as f32;
                let fy = ((a - min_a) / (max_a - min_a).max(1.0)).clamp(0.0, 1.0) as f32;
                egui::pos2(rect.left() + fx * rect.width(), rect.bottom() - fy * rect.height())
            };

            let mut curve_pts = Vec::with_capacity(60);
            for k in 1..=60 {
                let ci = 45.0 + ((k as f64) / 60.0) * (max_ci - 45.0);
                let (a_val, _) = eval_fvcb_assimilation(
                    ci,
                    1500.0,
                    25.0,
                    self.aci_fit.vcmax,
                    self.aci_fit.jmax,
                    self.aci_fit.tpu,
                    self.aci_fit.rd,
                );
                curve_pts.push(map_pt(ci, a_val));
            }
            painter.add(egui::Shape::line(
                curve_pts,
                egui::Stroke::new(2.4, egui::Color32::from_rgb(22, 163, 74)),
            ));

            for (p, st) in self.aci_pts.iter().zip(&self.aci_fit.limiting_states) {
                let col = match st {
                    crate::flux_leaf::FvcbLimitationState::RubiscoWc => {
                        egui::Color32::from_rgb(37, 99, 235)
                    }
                    crate::flux_leaf::FvcbLimitationState::RuBpRegenWj => {
                        egui::Color32::from_rgb(217, 119, 6)
                    }
                    crate::flux_leaf::FvcbLimitationState::TriosePhosphateWp => {
                        egui::Color32::from_rgb(220, 38, 38)
                    }
                };
                painter.circle_filled(map_pt(p.ci_ubar, p.a_net), 4.5, col);
            }

            cols[1].heading("Farquhar–von Caemmerer–Berry (1980) Parameters");
            cols[1].separator();
            cols[1].label(format!("Vcmax (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.vcmax));
            cols[1].label(format!("Jmax  (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.jmax));
            cols[1].label(format!("TPU   (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.tpu));
            cols[1].label(format!("Rd    (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.rd));
            cols[1].label(format!("Jmax / Vcmax : {:.3}", self.aci_fit.jmax_to_vcmax_ratio));
            cols[1].label(format!(
                "RMSE = {:.3} | R^2 = {:.5}",
                self.aci_fit.rmse, self.aci_fit.r_squared
            ));
            cols[1].separator();
            cols[1].label("Point Color Key: Blue = Rubisco (Wc), Amber = RuBP Regen (Wj), Red = TPU (Wp)");
        });
    }

    fn draw_stream_panel(&self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            let (resp, painter) = cols[0].allocate_painter(
                egui::vec2(cols[0].available_width().max(480.0), 520.0),
                egui::Sense::hover(),
            );
            let rect = resp.rect;
            painter.rect_filled(rect, 6.0, egui::Color32::from_rgb(248, 250, 252));

            let max_hr = self.stream_pts.last().map(|p| p.hour).unwrap_or(24.0).max(1.0);
            let map_pt = |h: f64, do_v: f64| -> egui::Pos2 {
                let fx = (h / max_hr).clamp(0.0, 1.0) as f32;
                let fy = ((do_v - 6.0) / 7.0).clamp(0.0, 1.0) as f32;
                egui::pos2(rect.left() + fx * rect.width(), rect.bottom() - fy * rect.height())
            };

            let mod_line: Vec<egui::Pos2> = self
                .stream_pts
                .iter()
                .zip(&self.stream_rep.modeled_do_mg_l)
                .map(|(p, &m)| map_pt(p.hour, m))
                .collect();
            let sat_line: Vec<egui::Pos2> = self
                .stream_pts
                .iter()
                .zip(&self.stream_rep.do_sat_mg_l)
                .map(|(p, &s)| map_pt(p.hour, s))
                .collect();
            painter.add(egui::Shape::line(
                sat_line,
                egui::Stroke::new(1.5, egui::Color32::from_rgb(100, 116, 139)),
            ));
            painter.add(egui::Shape::line(
                mod_line,
                egui::Stroke::new(2.4, egui::Color32::from_rgb(37, 99, 235)),
            ));
            for p in &self.stream_pts {
                painter.circle_filled(
                    map_pt(p.hour, p.do_mg_l),
                    3.8,
                    egui::Color32::from_rgb(220, 38, 38),
                );
            }

            cols[1].heading("Odum (1956) Open-Channel Stream Metabolism");
            cols[1].separator();
            cols[1].label(format!(
                "GPP  : {:.3} g O2 m^-2 d^-1",
                self.stream_rep.gpp_g_o2_m2_d
            ));
            cols[1].label(format!(
                "ER   : {:.3} g O2 m^-2 d^-1",
                self.stream_rep.er_g_o2_m2_d
            ));
            cols[1].label(format!(
                "NEP  : {:+.3} g O2 m^-2 d^-1",
                self.stream_rep.nep_g_o2_m2_d
            ));
            cols[1].label(format!("P/R  : {:.3}", self.stream_rep.p_to_r_ratio));
            cols[1].label(format!("K600 : {:.2} d^-1", self.stream_rep.k600_per_day));
            cols[1].label(format!("State: {}", self.stream_rep.trophic_state));
            cols[1].label(format!("RMSE : {:.4} mg/L", self.stream_rep.do_rmse_mg_l));
        });
    }
}
