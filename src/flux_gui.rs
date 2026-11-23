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

//! Interactive `eframe`/`egui` Leaf `FvCB` Photosynthesis & Stream Open-Channel `GPP`/`ER`
//! Metabolism Studio (`osabt flux`).

use std::fs;
use std::io::BufWriter;
use std::path::PathBuf;

use eframe::egui;

use crate::gui_io::{
    draw_status_banner, pick_open_file, pick_save_file, read_text_reported, write_text_reported,
};

use crate::flux_leaf::{
    demo_fvcb_aci_csv, eval_fvcb_assimilation, fit_fvcb_aci_curve, parse_aci_csv, FvcbAciFit,
    LeafAciPoint,
};
use crate::flux_stream::{
    compute_do_y_bounds, demo_stream_diel_do_csv, fit_odum_stream_metabolism, parse_stream_do_csv,
    render_flux_svg, StreamDoPoint, StreamMetabolismReport,
};
use crate::kinetics_core::{
    demo_competitive_kinase_matrix_csv, demo_morrison_tight_binding_csv, eval_inhibition_velocity,
    eval_morrison_velocity, fit_global_reversible_inhibition, fit_morrison_tight_binding,
    parse_kinetics_csv, GlobalInhibitionFit, KineticsPoint, MorrisonTightBindingFit,
};
use crate::plasmid_digest::{
    find_unique_cutters, simulate_restriction_digest, CutSite, DigestFragment,
};
use crate::plasmid_io::PlasmidRecord;
use crate::survival_core::{
    compute_kaplan_meier_curves, compute_log_rank_test, demo_survival_cohort_csv,
    fit_cox_proportional_hazards, parse_survival_csv, CoxPhFit, KaplanMeierCurve,
    LogRankTestResult,
};
use crate::survival_tox::{
    demo_ecotox_lc50_csv, fit_ecotox_lc50, parse_ecotox_csv, DoseResponseLink, EcotoxBioassayPoint,
    EcotoxLc50Report,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum FluxView {
    DualOverview,
    LeafAciFvcb,
    StreamDielDo,
    SurvivalAndEcotox,
    KineticsAndPlasmid,
}

/// Captures a 1600x1000 offscreen GPU screenshot of OSABT Flux Studio without opening a window.
pub fn run_flux_screenshot(
    aci_path: Option<PathBuf>,
    stream_path: Option<PathBuf>,
    out_path: &std::path::Path,
    view_idx: u8,
) -> Result<(), String> {
    let mut app = FluxGuiApp::new(aci_path, stream_path)?;
    app.view = match view_idx {
        1 => FluxView::DualOverview,
        2 => FluxView::SurvivalAndEcotox,
        3 => FluxView::KineticsAndPlasmid,
        other => {
            return Err(format!(
                "Unsupported --screenshot-view {other} (expected 1, 2, or 3)"
            ));
        }
    };
    let ctx = egui::Context::default();
    render_flux_offscreen_png(&ctx, out_path, |c| app.draw_panels(c))
}

/// Launches the `eframe`/`egui` native desktop studio for `osabt flux`.
pub fn run_flux_gui(aci_path: Option<PathBuf>, stream_path: Option<PathBuf>) -> Result<(), String> {
    let app = FluxGuiApp::new(aci_path, stream_path)?;
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_icon(crate::brand::osabt_icon_data())
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

pub struct FluxGuiApp {
    view: FluxView,
    aci_pts: Vec<LeafAciPoint>,
    aci_fit: FvcbAciFit,
    stream_pts: Vec<StreamDoPoint>,
    stream_rep: StreamMetabolismReport,
    km_curves: Vec<KaplanMeierCurve>,
    log_rank: LogRankTestResult,
    cox_fit: CoxPhFit,
    ecotox_pts: Vec<EcotoxBioassayPoint>,
    ecotox_fit: EcotoxLc50Report,
    kin_pts: Vec<KineticsPoint>,
    kin_fits: Vec<GlobalInhibitionFit>,
    morrison_obs: Vec<(f64, f64)>,
    morrison_fit: MorrisonTightBindingFit,
    plasmid: PlasmidRecord,
    unique_cuts: Vec<CutSite>,
    digest_enzymes: Vec<String>,
    digest_frags: Vec<DigestFragment>,
    status_msg: String,
    status_is_error: bool,
    about_open: bool,
}

impl FluxGuiApp {
    pub fn new(aci_path: Option<PathBuf>, stream_path: Option<PathBuf>) -> Result<Self, String> {
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

        let surv_obs = parse_survival_csv(&demo_survival_cohort_csv())?;
        let km_curves = compute_kaplan_meier_curves(&surv_obs);
        let log_rank = compute_log_rank_test(&surv_obs)?;
        let cox_fit = fit_cox_proportional_hazards(&surv_obs, "Group (Hatchery=1)", true)?;

        let ecotox_pts = parse_ecotox_csv(&demo_ecotox_lc50_csv())?;
        let ecotox_fit = fit_ecotox_lc50(&ecotox_pts, DoseResponseLink::Probit)?;

        let kin_pts = parse_kinetics_csv(&demo_competitive_kinase_matrix_csv())?;
        let kin_fits = fit_global_reversible_inhibition(&kin_pts)?;
        let morrison_obs: Vec<(f64, f64)> = demo_morrison_tight_binding_csv()
            .lines()
            .filter_map(|line| {
                let t = line.trim();
                if t.is_empty() || t.starts_with('#') {
                    return None;
                }
                let sep = if t.contains('\t') { '\t' } else { ',' };
                let mut cols = t.split(sep).map(|c| c.trim());
                let x = cols.next()?.parse::<f64>().ok()?;
                let y = cols.next()?.parse::<f64>().ok()?;
                Some((x, y))
            })
            .collect();
        let morrison_fit = fit_morrison_tight_binding(
            &morrison_obs,
            None,
            Some(25.0),
            Some((50.0, 20.0)),
        )?;

        let plasmid = PlasmidRecord::demo_pcrv_gfp();
        let unique_cuts = find_unique_cutters(&plasmid);
        let mut digest_enzymes: Vec<String> = Vec::new();
        for c in &unique_cuts {
            if !digest_enzymes.iter().any(|e| e == &c.enzyme_name) {
                digest_enzymes.push(c.enzyme_name.clone());
            }
            if digest_enzymes.len() >= 3 {
                break;
            }
        }
        let digest_frags =
            simulate_restriction_digest(&plasmid, &digest_enzymes).map_err(|e| format!("{e:?}"))?;

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
            km_curves,
            log_rank,
            cox_fit,
            ecotox_pts,
            ecotox_fit,
            kin_pts,
            kin_fits,
            morrison_obs,
            morrison_fit,
            plasmid,
            unique_cuts,
            digest_enzymes,
            digest_frags,
            status_msg,
            status_is_error: false,
            about_open: false,
        })
    }
}

impl eframe::App for FluxGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        crate::brand::draw_standalone_footer(ctx, &mut self.about_open);
        self.draw_panels(ctx);
    }
}

impl FluxGuiApp {
    /// Selects the dual Leaf FvCB + Odum Stream Metabolism overview tab.
    pub fn set_dual_overview(&mut self) {
        self.view = FluxView::DualOverview;
    }

    /// Renders all Flux Studio panels onto `ctx`.
    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("flux_top").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading("OSABT Flux v2.0.0");
                ui.separator();
                ui.selectable_value(
                    &mut self.view,
                    FluxView::DualOverview,
                    "Overview (FvCB Leaf + Odum Stream)",
                );
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
                ui.selectable_value(
                    &mut self.view,
                    FluxView::SurvivalAndEcotox,
                    "3. Survival & LC50 Ecotox",
                );
                ui.selectable_value(
                    &mut self.view,
                    FluxView::KineticsAndPlasmid,
                    "4. Kinetics & Plasmid Map",
                );
                ui.separator();
                if ui.button("Open A-Ci...").clicked() {
                    if let Some(p) = pick_open_file("Open Leaf A-Ci CSV", "CSV", &["csv", "tsv"]) {
                        if let Some(raw) = read_text_reported(&p, &mut self.status_msg, &mut self.status_is_error) {
                            match parse_aci_csv(&raw).and_then(|pts| fit_fvcb_aci_curve(&pts).map(|f| (pts, f))) {
                                Ok((pts, fit)) => {
                                    self.aci_pts = pts;
                                    self.aci_fit = fit;
                                    self.status_is_error = false;
                                    self.status_msg = format!(
                                        "Loaded A-Ci '{}': Vcmax={:.1}, Jmax={:.1}, TPU={:.1} (R^2={:.4})",
                                        p.display(), self.aci_fit.vcmax, self.aci_fit.jmax, self.aci_fit.tpu, self.aci_fit.r_squared
                                    );
                                }
                                Err(e) => {
                                    self.status_is_error = true;
                                    self.status_msg = format!("A-Ci parse/fit error: {e}");
                                }
                            }
                        }
                    }
                }
                if ui.button("Open Stream DO...").clicked() {
                    if let Some(p) = pick_open_file("Open Stream Diel DO CSV", "CSV", &["csv", "tsv"]) {
                        if let Some(raw) = read_text_reported(&p, &mut self.status_msg, &mut self.status_is_error) {
                            match parse_stream_do_csv(&raw).and_then(|pts| fit_odum_stream_metabolism(&pts).map(|r| (pts, r))) {
                                Ok((pts, rep)) => {
                                    self.stream_pts = pts;
                                    self.stream_rep = rep;
                                    self.status_is_error = false;
                                    self.status_msg = format!(
                                        "Loaded Stream DO '{}': GPP={:.2}, ER={:.2}, NEP={:+.2} g O2/m2/d",
                                        p.display(), self.stream_rep.gpp_g_o2_m2_d, self.stream_rep.er_g_o2_m2_d, self.stream_rep.nep_g_o2_m2_d
                                    );
                                }
                                Err(e) => {
                                    self.status_is_error = true;
                                    self.status_msg = format!("Stream DO parse/fit error: {e}");
                                }
                            }
                        }
                    }
                }
                if ui.button("Export Dual-Panel .svg").clicked() {
                    if let Some(path) = pick_save_file("Export Flux SVG", "osabt_flux_report.svg", "SVG Image", &["svg"]) {
                        let svg = render_flux_svg(
                            &self.aci_pts,
                            &self.aci_fit,
                            &self.stream_pts,
                            &self.stream_rep,
                            "OSABT Ecophysiology & Stream Flux",
                        );
                        write_text_reported(
                            &path,
                            &svg,
                            &mut self.status_msg,
                            &mut self.status_is_error,
                            "Exported Flux SVG",
                        );
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("flux_bot").show(ctx, |ui| {
            if self.status_is_error {
                draw_status_banner(ui, &self.status_msg, true);
            } else {
                match self.view {
                    FluxView::SurvivalAndEcotox => {
                        ui.label(format!(
                            "Kaplan-Meier Log-Rank: chi^2 = {:.3}, p = {:.4} | Cox HR = {:.2} (95% CI {:.2}-{:.2}) | Probit LC50 = {:.2} (95% CI {:.2}-{:.2})",
                            self.log_rank.chi_square,
                            self.log_rank.p_value,
                            self.cox_fit.hazard_ratio,
                            self.cox_fit.hr_ci95.0,
                            self.cox_fit.hr_ci95.1,
                            self.ecotox_fit.lc50,
                            self.ecotox_fit.lc50_ci95.0,
                            self.ecotox_fit.lc50_ci95.1
                        ));
                    }
                    FluxView::KineticsAndPlasmid => {
                        if let Some(top) = self.kin_fits.first() {
                            ui.label(format!(
                                "Top Inhibition Fit: {} (Km={:.1} uM, Ki={:.1} nM, R^2={:.4}) | Morrison [E]_T={:.1} nM, Ki_app={:.1} nM | {} ({} bp, {} digest frags via {})",
                                top.mechanism.label(),
                                top.km_um,
                                top.ki_nm,
                                top.r_squared,
                                self.morrison_fit.active_enzyme_et_nm,
                                self.morrison_fit.ki_app_nm,
                                self.plasmid.name,
                                self.plasmid.sequence.len(),
                                self.digest_frags.len(),
                                self.digest_enzymes.join(" + ")
                            ));
                        } else {
                            draw_status_banner(ui, &self.status_msg, false);
                        }
                    }
                    _ => {
                        draw_status_banner(ui, &self.status_msg, false);
                    }
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| match self.view {
            FluxView::DualOverview => self.draw_dual_overview(ui),
            FluxView::LeafAciFvcb => self.draw_leaf_panel(ui),
            FluxView::StreamDielDo => self.draw_stream_panel(ui),
            FluxView::SurvivalAndEcotox => self.draw_survival_ecotox_panel(ui),
            FluxView::KineticsAndPlasmid => self.draw_kinetics_plasmid_panel(ui),
        });
    }

    fn draw_leaf_plot(&self, ui: &mut egui::Ui, width: f32, height: f32) {
        let (resp, painter) =
            ui.allocate_painter(egui::vec2(width.max(360.0), height), egui::Sense::hover());
        let outer = resp.rect;
        painter.rect_filled(outer, 6.0, egui::Color32::from_rgb(248, 250, 252));
        painter.rect_stroke(
            outer,
            6.0,
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
        );
        let rect = outer.shrink2(egui::vec2(48.0, 32.0));

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

        let grid_stroke = egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(226, 232, 240));
        let axis_color = egui::Color32::from_rgb(71, 85, 105);
        for tick in [0.0, 400.0, 800.0, 1200.0, 1600.0] {
            if tick <= max_ci {
                let p = map_pt(tick, min_a);
                painter.line_segment([egui::pos2(p.x, rect.top()), egui::pos2(p.x, rect.bottom())], grid_stroke);
                painter.text(
                    egui::pos2(p.x, rect.bottom() + 6.0),
                    egui::Align2::CENTER_TOP,
                    format!("{tick:.0}"),
                    egui::FontId::proportional(11.0),
                    axis_color,
                );
            }
        }
        for tick in [0.0, 10.0, 20.0, 30.0, 40.0] {
            if tick <= max_a {
                let p = map_pt(0.0, tick);
                painter.line_segment([egui::pos2(rect.left(), p.y), egui::pos2(rect.right(), p.y)], grid_stroke);
                painter.text(
                    egui::pos2(rect.left() - 6.0, p.y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{tick:.0}"),
                    egui::FontId::proportional(11.0),
                    axis_color,
                );
            }
        }
        painter.text(
            egui::pos2(rect.center().x, outer.top() + 8.0),
            egui::Align2::CENTER_TOP,
            "FvCB C3 Leaf Assimilation vs. Intercellular CO2 (Ci [ubar])",
            egui::FontId::proportional(13.0),
            egui::Color32::from_rgb(15, 23, 42),
        );

        let ppfd_ref = if self.aci_pts.is_empty() {
            1500.0
        } else {
            self.aci_pts.iter().map(|p| p.ppfd).sum::<f64>() / (self.aci_pts.len() as f64)
        };
        let t_leaf_ref = if self.aci_pts.is_empty() {
            25.0
        } else {
            self.aci_pts.iter().map(|p| p.t_leaf_c).sum::<f64>() / (self.aci_pts.len() as f64)
        };
        let mut curve_pts = Vec::with_capacity(60);
        for k in 1..=60 {
            let ci = 45.0 + ((k as f64) / 60.0) * (max_ci - 45.0);
            let (a_val, _) = eval_fvcb_assimilation(
                ci,
                ppfd_ref,
                t_leaf_ref,
                self.aci_fit.vcmax,
                self.aci_fit.jmax,
                self.aci_fit.tpu,
                self.aci_fit.rd,
            );
            curve_pts.push(map_pt(ci, a_val));
        }
        painter.add(egui::Shape::line(
            curve_pts,
            egui::Stroke::new(2.6_f32, egui::Color32::from_rgb(22, 163, 74)),
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
            painter.circle_filled(map_pt(p.ci_ubar, p.a_net), 4.8, col);
        }
    }

    fn draw_leaf_params(&self, ui: &mut egui::Ui) {
        ui.heading("Farquhar-von Caemmerer-Berry (1980) Parameters");
        ui.separator();
        ui.label(format!("Vcmax (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.vcmax));
        ui.label(format!("Jmax  (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.jmax));
        ui.label(format!("TPU   (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.tpu));
        ui.label(format!("Rd    (25 C) : {:.2} umol m^-2 s^-1", self.aci_fit.rd));
        ui.label(format!("Jmax / Vcmax : {:.3}", self.aci_fit.jmax_to_vcmax_ratio));
        ui.label(format!(
            "RMSE = {:.3} | R^2 = {:.5}",
            self.aci_fit.rmse, self.aci_fit.r_squared
        ));
        ui.separator();
        ui.label("Point Color Key: Blue = Rubisco (Wc), Amber = RuBP Regen (Wj), Red = TPU (Wp)");
    }

    fn draw_stream_plot(&self, ui: &mut egui::Ui, width: f32, height: f32) {
        let (resp, painter) =
            ui.allocate_painter(egui::vec2(width.max(360.0), height), egui::Sense::hover());
        let outer = resp.rect;
        painter.rect_filled(outer, 6.0, egui::Color32::from_rgb(248, 250, 252));
        painter.rect_stroke(
            outer,
            6.0,
            egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
        );
        let rect = outer.shrink2(egui::vec2(48.0, 32.0));

        let max_hr = self.stream_pts.last().map(|p| p.hour).unwrap_or(24.0).max(1.0);
        let (min_do, max_do) = compute_do_y_bounds(&self.stream_pts, &self.stream_rep);
        let do_span = (max_do - min_do).max(1e-6);
        let map_pt = |h: f64, do_v: f64| -> egui::Pos2 {
            let fx = (h / max_hr).clamp(0.0, 1.0) as f32;
            let fy = ((do_v - min_do) / do_span).clamp(0.0, 1.0) as f32;
            egui::pos2(rect.left() + fx * rect.width(), rect.bottom() - fy * rect.height())
        };

        let grid_stroke = egui::Stroke::new(0.8_f32, egui::Color32::from_rgb(226, 232, 240));
        let axis_color = egui::Color32::from_rgb(71, 85, 105);
        for hr in [0.0, 6.0, 12.0, 18.0, 24.0] {
            let p = map_pt(hr, min_do);
            painter.line_segment([egui::pos2(p.x, rect.top()), egui::pos2(p.x, rect.bottom())], grid_stroke);
            painter.text(
                egui::pos2(p.x, rect.bottom() + 6.0),
                egui::Align2::CENTER_TOP,
                format!("{hr:.0}h"),
                egui::FontId::proportional(11.0),
                axis_color,
            );
        }
        for k in 0..=4 {
            let do_tick = min_do + (k as f64 / 4.0) * do_span;
            let p = map_pt(0.0, do_tick);
            painter.line_segment([egui::pos2(rect.left(), p.y), egui::pos2(rect.right(), p.y)], grid_stroke);
            painter.text(
                egui::pos2(rect.left() - 6.0, p.y),
                egui::Align2::RIGHT_CENTER,
                format!("{do_tick:.1}"),
                egui::FontId::proportional(11.0),
                axis_color,
            );
        }
        painter.text(
            egui::pos2(rect.center().x, outer.top() + 8.0),
            egui::Align2::CENTER_TOP,
            "Odum Diel Dissolved O2 [mg/L]: Observed (Red), Modeled (Blue), DO Sat (Slate)",
            egui::FontId::proportional(13.0),
            egui::Color32::from_rgb(15, 23, 42),
        );

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
            egui::Stroke::new(1.8_f32, egui::Color32::from_rgb(100, 116, 139)),
        ));
        painter.add(egui::Shape::line(
            mod_line,
            egui::Stroke::new(2.6_f32, egui::Color32::from_rgb(37, 99, 235)),
        ));
        for p in &self.stream_pts {
            painter.circle_filled(
                map_pt(p.hour, p.do_mg_l),
                4.0,
                egui::Color32::from_rgb(220, 38, 38),
            );
        }
    }

    fn draw_stream_params(&self, ui: &mut egui::Ui) {
        ui.heading("Odum (1956) Open-Channel Stream Metabolism");
        ui.separator();
        ui.label(format!("GPP  : {:.3} g O2 m^-2 d^-1", self.stream_rep.gpp_g_o2_m2_d));
        ui.label(format!("ER   : {:.3} g O2 m^-2 d^-1", self.stream_rep.er_g_o2_m2_d));
        ui.label(format!("NEP  : {:+.3} g O2 m^-2 d^-1", self.stream_rep.nep_g_o2_m2_d));
        ui.label(format!("P/R  : {:.3}", self.stream_rep.p_to_r_ratio));
        ui.label(format!("K600 : {:.2} d^-1", self.stream_rep.k600_per_day));
        ui.label(format!("State: {}", self.stream_rep.trophic_state));
        ui.label(format!("RMSE : {:.4} mg/L", self.stream_rep.do_rmse_mg_l));
    }

    fn draw_dual_overview(&self, ui: &mut egui::Ui) {
        let plot_h = (ui.available_height() - 250.0).clamp(360.0, 580.0);
        ui.columns(2, |cols| {
            let w0 = cols[0].available_width();
            self.draw_leaf_plot(&mut cols[0], w0, plot_h);
            cols[0].add_space(10.0);
            self.draw_leaf_params(&mut cols[0]);

            let w1 = cols[1].available_width();
            self.draw_stream_plot(&mut cols[1], w1, plot_h);
            cols[1].add_space(10.0);
            self.draw_stream_params(&mut cols[1]);
        });
    }
}

impl FluxGuiApp {
    fn draw_leaf_panel(&self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            let w = cols[0].available_width().max(480.0);
            self.draw_leaf_plot(&mut cols[0], w, 520.0);
            self.draw_leaf_params(&mut cols[1]);
        });
    }

    fn draw_stream_panel(&self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            let w = cols[0].available_width().max(480.0);
            self.draw_stream_plot(&mut cols[0], w, 520.0);
            self.draw_stream_params(&mut cols[1]);
        });
    }

    fn draw_survival_ecotox_panel(&self, ui: &mut egui::Ui) {
        let avail_h = ui.available_height();
        let plot_h = (avail_h - 210.0).clamp(420.0, 710.0);
        ui.columns(2, |cols| {
            // Left: Kaplan-Meier Survival Curves + 95% CI + Mantel Log-Rank + Cox PH
            let w0 = cols[0].available_width();
            let (resp, painter) =
                cols[0].allocate_painter(egui::vec2(w0.max(420.0), plot_h), egui::Sense::hover());
            let outer = resp.rect;
            painter.rect_filled(outer, 6.0, egui::Color32::from_rgb(248, 250, 252));
            painter.rect_stroke(
                outer,
                6.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
            );
            let p = painter.with_clip_rect(outer);
            let rect = egui::Rect::from_min_max(
                egui::pos2(outer.left() + 56.0, outer.top() + 36.0),
                egui::pos2(outer.right() - 22.0, outer.bottom() - 38.0),
            );
            p.text(
                egui::pos2(rect.center().x, outer.top() + 10.0),
                egui::Align2::CENTER_TOP,
                "Kaplan-Meier Survival Analysis (95% Greenwood CI Bands)",
                egui::FontId::proportional(13.5),
                egui::Color32::from_rgb(15, 23, 42),
            );
            let max_t = self
                .km_curves
                .iter()
                .flat_map(|c| c.steps.iter().map(|s| s.time))
                .fold(10.0_f64, f64::max)
                * 1.05;
            let map_km = |t: f64, s: f64| -> egui::Pos2 {
                let fx = (t / max_t).clamp(0.0, 1.0) as f32;
                let fy = s.clamp(0.0, 1.0) as f32;
                egui::pos2(
                    rect.left() + fx * rect.width(),
                    rect.bottom() - fy * rect.height(),
                )
            };
            let grid_s = egui::Stroke::new(0.9_f32, egui::Color32::from_rgb(226, 232, 240));
            let axis_c = egui::Color32::from_rgb(71, 85, 105);
            for k in 0..=5 {
                let sv = k as f64 * 0.2;
                let pt = map_km(0.0, sv);
                p.line_segment(
                    [egui::pos2(rect.left(), pt.y), egui::pos2(rect.right(), pt.y)],
                    grid_s,
                );
                p.text(
                    egui::pos2(rect.left() - 6.0, pt.y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{sv:.1}"),
                    egui::FontId::proportional(10.5),
                    axis_c,
                );
            }
            for k in 0..=5 {
                let tv = max_t * (k as f64 / 5.0);
                let pt = map_km(tv, 0.0);
                p.line_segment(
                    [egui::pos2(pt.x, rect.top()), egui::pos2(pt.x, rect.bottom())],
                    grid_s,
                );
                p.text(
                    egui::pos2(pt.x, rect.bottom() + 6.0),
                    egui::Align2::CENTER_TOP,
                    format!("{tv:.0}d"),
                    egui::FontId::proportional(10.5),
                    axis_c,
                );
            }
            p.rect_stroke(
                rect,
                0.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
            );
            let y_50 = map_km(0.0, 0.5).y;
            p.line_segment(
                [egui::pos2(rect.left(), y_50), egui::pos2(rect.right(), y_50)],
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
            );
            p.text(
                egui::pos2(rect.left() + 6.0, y_50 - 3.0),
                egui::Align2::LEFT_BOTTOM,
                "Median S(t) = 0.50",
                egui::FontId::proportional(9.5),
                egui::Color32::from_rgb(100, 116, 139),
            );

            let palette = [
                (
                    egui::Color32::from_rgb(37, 99, 235),
                    egui::Color32::from_rgba_premultiplied(37, 99, 235, 28),
                ),
                (
                    egui::Color32::from_rgb(220, 38, 38),
                    egui::Color32::from_rgba_premultiplied(220, 38, 38, 28),
                ),
                (
                    egui::Color32::from_rgb(5, 150, 105),
                    egui::Color32::from_rgba_premultiplied(5, 150, 105, 28),
                ),
            ];
            let leg_box = egui::Rect::from_min_max(
                egui::pos2(rect.left() + 14.0, rect.bottom() - 20.0 - self.km_curves.len() as f32 * 20.0),
                egui::pos2(rect.left() + 300.0, rect.bottom() - 12.0),
            );
            for (ci, curve) in self.km_curves.iter().enumerate() {
                let (line_col, band_col) = palette[ci % palette.len()];
                let mut prev_t = 0.0_f64;
                let mut prev_s = 1.0_f64;
                let mut prev_ci = (1.0_f64, 1.0_f64);
                let mut line_pts = vec![map_km(0.0, 1.0)];
                for step in &curve.steps {
                    let p_top_left = map_km(prev_t, prev_ci.1);
                    let p_bot_right = map_km(step.time, prev_ci.0);
                    let band_rect = egui::Rect::from_min_max(
                        egui::pos2(p_top_left.x.min(p_bot_right.x), p_top_left.y.min(p_bot_right.y)),
                        egui::pos2(p_top_left.x.max(p_bot_right.x), p_top_left.y.max(p_bot_right.y)),
                    );
                    p.rect_filled(band_rect, 0.0, band_col);

                    line_pts.push(map_km(step.time, prev_s));
                    line_pts.push(map_km(step.time, step.survival_s));
                    if step.n_censored > 0 {
                        let cp = map_km(step.time, step.survival_s);
                        p.line_segment(
                            [cp - egui::vec2(3.5, 0.0), cp + egui::vec2(3.5, 0.0)],
                            egui::Stroke::new(1.5_f32, line_col),
                        );
                        p.line_segment(
                            [cp - egui::vec2(0.0, 3.5), cp + egui::vec2(0.0, 3.5)],
                            egui::Stroke::new(1.5_f32, line_col),
                        );
                    }
                    prev_t = step.time;
                    prev_s = step.survival_s;
                    prev_ci = step.ci95;
                }
                p.add(egui::Shape::line(
                    line_pts,
                    egui::Stroke::new(2.5_f32, line_col),
                ));
            }
            p.rect_filled(leg_box, 4.0, egui::Color32::from_rgba_premultiplied(255, 255, 255, 232));
            p.rect_stroke(
                leg_box,
                4.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(203, 213, 225)),
            );
            for (ci, curve) in self.km_curves.iter().enumerate() {
                let (line_col, _) = palette[ci % palette.len()];
                let leg_y = leg_box.top() + 13.0 + ci as f32 * 18.0;
                p.line_segment(
                    [
                        egui::pos2(leg_box.left() + 10.0, leg_y),
                        egui::pos2(leg_box.left() + 30.0, leg_y),
                    ],
                    egui::Stroke::new(2.8_f32, line_col),
                );
                let med_txt = curve
                    .median_survival_time
                    .map(|m| format!("{m:.1}d"))
                    .unwrap_or_else(|| "NR".to_string());
                p.text(
                    egui::pos2(leg_box.left() + 36.0, leg_y),
                    egui::Align2::LEFT_CENTER,
                    format!(
                        "{} (n={}, ev={}, med={})",
                        curve.group, curve.total_n, curve.total_events, med_txt
                    ),
                    egui::FontId::proportional(10.5),
                    egui::Color32::from_rgb(15, 23, 42),
                );
            }

            cols[0].add_space(8.0);
            cols[0].heading("Kaplan-Meier, Mantel Log-Rank & Cox PH");
            cols[0].separator();
            cols[0].label(format!(
                "Mantel Log-Rank ({} vs {}): O_A={:.1} (E={:.2}) | O_B={:.1} (E={:.2})",
                self.log_rank.group_a,
                self.log_rank.group_b,
                self.log_rank.observed_a,
                self.log_rank.expected_a,
                self.log_rank.observed_b,
                self.log_rank.expected_b
            ));
            cols[0].label(format!(
                "Log-Rank Statistic: chi^2 = {:.4}   (p = {:.5})",
                self.log_rank.chi_square, self.log_rank.p_value
            ));
            cols[0].label(format!(
                "Cox PH ({}): beta = {:+.4} +/- {:.4}  |  Wald Z = {:+.3}  (p = {:.5})",
                self.cox_fit.covariate_name,
                self.cox_fit.beta,
                self.cox_fit.se_beta,
                self.cox_fit.wald_z,
                self.cox_fit.p_value
            ));
            cols[0].label(format!(
                "Hazard Ratio (HR) : {:.3}   (95% CI: {:.3} .. {:.3})",
                self.cox_fit.hazard_ratio, self.cox_fit.hr_ci95.0, self.cox_fit.hr_ci95.1
            ));
            cols[0].label(format!(
                "Partial Log-Lik   : null = {:.3}  ->  fit = {:.3}   (events = {} / {})",
                self.cox_fit.partial_log_lik_null,
                self.cox_fit.partial_log_lik_fit,
                self.cox_fit.n_events,
                self.cox_fit.n_subjects
            ));

            // Right: Finney Probit Acute Toxicity LC50 Curve
            let w1 = cols[1].available_width();
            let (resp2, painter2) =
                cols[1].allocate_painter(egui::vec2(w1.max(420.0), plot_h), egui::Sense::hover());
            let outer2 = resp2.rect;
            painter2.rect_filled(outer2, 6.0, egui::Color32::from_rgb(248, 250, 252));
            painter2.rect_stroke(
                outer2,
                6.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
            );
            let p2 = painter2.with_clip_rect(outer2);
            let r2 = egui::Rect::from_min_max(
                egui::pos2(outer2.left() + 56.0, outer2.top() + 36.0),
                egui::pos2(outer2.right() - 24.0, outer2.bottom() - 38.0),
            );
            p2.text(
                egui::pos2(r2.center().x, outer2.top() + 10.0),
                egui::Align2::CENTER_TOP,
                "Ecotoxicology Bioassay: Finney (1971) Probit LC10 / LC50 / LC90 Curve",
                egui::FontId::proportional(13.5),
                egui::Color32::from_rgb(15, 23, 42),
            );
            let min_lx = self
                .ecotox_pts
                .first()
                .map(|x| x.concentration.log10())
                .unwrap_or(-1.0)
                - 0.35;
            let max_lx = self
                .ecotox_pts
                .last()
                .map(|x| x.concentration.log10())
                .unwrap_or(2.0)
                + 0.35;
            let map_tox = |lx: f64, mort: f64| -> egui::Pos2 {
                let fx = ((lx - min_lx) / (max_lx - min_lx).max(0.1)).clamp(0.0, 1.0) as f32;
                let fy = mort.clamp(0.0, 1.0) as f32;
                egui::pos2(r2.left() + fx * r2.width(), r2.bottom() - fy * r2.height())
            };
            for k in 0..=5 {
                let mv = k as f64 * 0.2;
                let pt = map_tox(min_lx, mv);
                p2.line_segment(
                    [egui::pos2(r2.left(), pt.y), egui::pos2(r2.right(), pt.y)],
                    grid_s,
                );
                p2.text(
                    egui::pos2(r2.left() - 6.0, pt.y),
                    egui::Align2::RIGHT_CENTER,
                    format!("{:.0}%", mv * 100.0),
                    egui::FontId::proportional(10.5),
                    axis_c,
                );
            }
            for pt_obs in &self.ecotox_pts {
                let lx = pt_obs.concentration.log10();
                let pt = map_tox(lx, 0.0);
                p2.line_segment(
                    [egui::pos2(pt.x, r2.top()), egui::pos2(pt.x, r2.bottom())],
                    grid_s,
                );
                p2.text(
                    egui::pos2(pt.x, r2.bottom() + 6.0),
                    egui::Align2::CENTER_TOP,
                    format!("{:.1}", pt_obs.concentration),
                    egui::FontId::proportional(10.0),
                    axis_c,
                );
            }
            p2.rect_stroke(
                r2,
                0.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
            );

            // LC50 95% CI vertical band
            let ci_l = map_tox(self.ecotox_fit.lc50_ci95.0.max(1e-9).log10(), 1.0);
            let ci_r = map_tox(self.ecotox_fit.lc50_ci95.1.max(1e-9).log10(), 0.0);
            p2.rect_filled(
                egui::Rect::from_min_max(ci_l, ci_r),
                0.0,
                egui::Color32::from_rgba_premultiplied(5, 150, 105, 24),
            );
            let lc50_pt = map_tox(self.ecotox_fit.lc50.max(1e-9).log10(), 0.5);
            p2.line_segment(
                [egui::pos2(r2.left(), lc50_pt.y), lc50_pt],
                egui::Stroke::new(1.2_f32, egui::Color32::from_rgb(5, 150, 105)),
            );
            p2.line_segment(
                [lc50_pt, egui::pos2(lc50_pt.x, r2.bottom())],
                egui::Stroke::new(1.2_f32, egui::Color32::from_rgb(5, 150, 105)),
            );

            let mut curve_pts = Vec::with_capacity(101);
            for s in 0..=100 {
                let lx = min_lx + (max_lx - min_lx) * (s as f64 / 100.0);
                let z = self.ecotox_fit.intercept_alpha + self.ecotox_fit.slope_beta * lx;
                let prob = probit_norm_cdf(z);
                curve_pts.push(map_tox(lx, prob));
            }
            p2.add(egui::Shape::line(
                curve_pts,
                egui::Stroke::new(2.6_f32, egui::Color32::from_rgb(37, 99, 235)),
            ));

            for pt_obs in &self.ecotox_pts {
                let n = pt_obs.exposed_n.max(1) as f64;
                let prop = pt_obs.dead_n as f64 / n;
                let se = (prop * (1.0 - prop) / n).sqrt();
                let lx = pt_obs.concentration.log10();
                let p_lo = map_tox(lx, (prop - 1.96 * se).max(0.0));
                let p_hi = map_tox(lx, (prop + 1.96 * se).min(1.0));
                p2.line_segment(
                    [p_lo, p_hi],
                    egui::Stroke::new(1.5_f32, egui::Color32::from_rgb(220, 38, 38)),
                );
                p2.circle_filled(
                    map_tox(lx, prop),
                    4.5,
                    egui::Color32::from_rgb(220, 38, 38),
                );
            }
            p2.circle_filled(lc50_pt, 5.5, egui::Color32::from_rgb(5, 150, 105));
            p2.text(
                lc50_pt + egui::vec2(14.0, 12.0),
                egui::Align2::LEFT_TOP,
                format!(
                    "LC50 = {:.2}\n(95% CI: {:.2}-{:.2})",
                    self.ecotox_fit.lc50,
                    self.ecotox_fit.lc50_ci95.0,
                    self.ecotox_fit.lc50_ci95.1
                ),
                egui::FontId::proportional(11.0),
                egui::Color32::from_rgb(6, 95, 70),
            );

            cols[1].add_space(8.0);
            cols[1].heading("Finney Probit LCx Dose-Response Parameters");
            cols[1].separator();
            cols[1].label(format!(
                "Model : {}  (alpha = {:+.3}, beta = {:+.3})",
                self.ecotox_fit.link.label(),
                self.ecotox_fit.intercept_alpha,
                self.ecotox_fit.slope_beta
            ));
            cols[1].label(format!(
                "LC10  : {:.3}   |   LC50 : {:.3} (95% CI: {:.3} .. {:.3})   |   LC90 : {:.3}",
                self.ecotox_fit.lc10,
                self.ecotox_fit.lc50,
                self.ecotox_fit.lc50_ci95.0,
                self.ecotox_fit.lc50_ci95.1,
                self.ecotox_fit.lc90
            ));
            cols[1].label(format!(
                "Log-Likelihood = {:.3}   |   Pearson chi^2 = {:.3}   (n = {} bioassay cohorts)",
                self.ecotox_fit.log_likelihood,
                self.ecotox_fit.pearson_chi2,
                self.ecotox_pts.len()
            ));
            let cohort_str: Vec<String> = self
                .ecotox_pts
                .iter()
                .map(|c| format!("{:.0}:{}/{}", c.concentration, c.dead_n, c.exposed_n))
                .collect();
            cols[1].label(format!("Cohorts (Conc:Dead/N): {}", cohort_str.join(" | ")));
        });
    }

    fn draw_kinetics_plasmid_panel(&self, ui: &mut egui::Ui) {
        let avail_h = ui.available_height();
        let top_h = ((avail_h - 190.0) * 0.55).clamp(240.0, 410.0);
        let bot_h = ((avail_h - 190.0) * 0.45).clamp(200.0, 330.0);
        let grid_s = egui::Stroke::new(0.9_f32, egui::Color32::from_rgb(226, 232, 240));
        let axis_c = egui::Color32::from_rgb(71, 85, 105);

        ui.columns(2, |cols| {
            // Left Column: Top = Global 2D Inhibition V vs [S], Bottom = Morrison Tight-Binding
            let w0 = cols[0].available_width().max(240.0);
            let Some(best_kin) = self.kin_fits.first() else {
                cols[0].label("No global inhibition fit available.");
                return;
            };
            {
                let (resp, painter) =
                    cols[0].allocate_painter(egui::vec2(w0, top_h), egui::Sense::hover());
                let outer = resp.rect;
                painter.rect_filled(outer, 6.0, egui::Color32::from_rgb(248, 250, 252));
                painter.rect_stroke(
                    outer,
                    6.0,
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
                );
                let p = painter.with_clip_rect(outer);
                let rect = egui::Rect::from_min_max(
                    egui::pos2(outer.left() + 52.0, outer.top() + 30.0),
                    egui::pos2(outer.right() - 20.0, outer.bottom() - 32.0),
                );
                p.text(
                    egui::pos2(rect.center().x, outer.top() + 8.0),
                    egui::Align2::CENTER_TOP,
                    format!(
                        "Global 2D Kinase Inhibition (AICc Winner: {})",
                        best_kin.mechanism.label()
                    ),
                    egui::FontId::proportional(12.5),
                    egui::Color32::from_rgb(15, 23, 42),
                );
                let max_s = self
                    .kin_pts
                    .iter()
                    .map(|k| k.substrate_um)
                    .fold(100.0_f64, f64::max)
                    * 1.06;
                let max_v = self
                    .kin_pts
                    .iter()
                    .map(|k| k.velocity)
                    .fold(10.0_f64, f64::max)
                    .max(best_kin.vmax)
                    * 1.1;
                let map_k = |s: f64, v: f64| -> egui::Pos2 {
                    let fx = (s / max_s).clamp(0.0, 1.0) as f32;
                    let fy = (v / max_v).clamp(0.0, 1.0) as f32;
                    egui::pos2(
                        rect.left() + fx * rect.width(),
                        rect.bottom() - fy * rect.height(),
                    )
                };
                for k in 0..=4 {
                    let vv = max_v * (k as f64 / 4.0);
                    let pt = map_k(0.0, vv);
                    p.line_segment(
                        [egui::pos2(rect.left(), pt.y), egui::pos2(rect.right(), pt.y)],
                        grid_s,
                    );
                    p.text(
                        egui::pos2(rect.left() - 6.0, pt.y),
                        egui::Align2::RIGHT_CENTER,
                        format!("{vv:.1}"),
                        egui::FontId::proportional(10.0),
                        axis_c,
                    );
                }
                for k in 0..=5 {
                    let sv = max_s * (k as f64 / 5.0);
                    let pt = map_k(sv, 0.0);
                    p.line_segment(
                        [egui::pos2(pt.x, rect.top()), egui::pos2(pt.x, rect.bottom())],
                        grid_s,
                    );
                    p.text(
                        egui::pos2(pt.x, rect.bottom() + 4.0),
                        egui::Align2::CENTER_TOP,
                        format!("{sv:.0} uM"),
                        egui::FontId::proportional(10.0),
                        axis_c,
                    );
                }
                p.rect_stroke(
                    rect,
                    0.0,
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
                );

                let mut i_levels: Vec<f64> = Vec::new();
                for pt in &self.kin_pts {
                    if !i_levels.iter().any(|&v| (v - pt.inhibitor_nm).abs() < 1e-6) {
                        i_levels.push(pt.inhibitor_nm);
                    }
                }
                i_levels.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
                let cols_i = [
                    egui::Color32::from_rgb(37, 99, 235),
                    egui::Color32::from_rgb(5, 150, 105),
                    egui::Color32::from_rgb(217, 119, 6),
                    egui::Color32::from_rgb(220, 38, 38),
                ];
                let alpha = best_kin.alpha_factor.unwrap_or(1.0);
                for (idx, &inh_nm) in i_levels.iter().enumerate() {
                    let c = cols_i[idx % cols_i.len()];
                    let curve: Vec<egui::Pos2> = (0..=80)
                        .map(|step| {
                            let s = max_s * (step as f64 / 80.0);
                            let v = eval_inhibition_velocity(
                                best_kin.mechanism,
                                s,
                                inh_nm,
                                best_kin.vmax,
                                best_kin.km_um,
                                best_kin.ki_nm,
                                alpha,
                            );
                            map_k(s, v)
                        })
                        .collect();
                    p.add(egui::Shape::line(curve, egui::Stroke::new(2.2_f32, c)));
                    for pt in self
                        .kin_pts
                        .iter()
                        .filter(|k| (k.inhibitor_nm - inh_nm).abs() < 1e-6)
                    {
                        p.circle_filled(map_k(pt.substrate_um, pt.velocity), 3.6, c);
                    }
                    let ly = rect.bottom() - 16.0 - (i_levels.len() - 1 - idx) as f32 * 15.0;
                    p.line_segment(
                        [
                            egui::pos2(rect.right() - 108.0, ly),
                            egui::pos2(rect.right() - 90.0, ly),
                        ],
                        egui::Stroke::new(2.4_f32, c),
                    );
                    p.text(
                        egui::pos2(rect.right() - 84.0, ly),
                        egui::Align2::LEFT_CENTER,
                        format!("[I] = {inh_nm:.0} nM"),
                        egui::FontId::proportional(10.0),
                        egui::Color32::from_rgb(15, 23, 42),
                    );
                }
            }

            cols[0].add_space(6.0);
            {
                let (resp, painter) =
                    cols[0].allocate_painter(egui::vec2(w0, bot_h), egui::Sense::hover());
                let outer = resp.rect;
                painter.rect_filled(outer, 6.0, egui::Color32::from_rgb(248, 250, 252));
                painter.rect_stroke(
                    outer,
                    6.0,
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
                );
                let p = painter.with_clip_rect(outer);
                let rect = egui::Rect::from_min_max(
                    egui::pos2(outer.left() + 52.0, outer.top() + 28.0),
                    egui::pos2(outer.right() - 20.0, outer.bottom() - 30.0),
                );
                p.text(
                    egui::pos2(rect.center().x, outer.top() + 8.0),
                    egui::Align2::CENTER_TOP,
                    "Morrison (1969) Tight-Binding Titration: v_i vs [I] (nM)",
                    egui::FontId::proportional(12.5),
                    egui::Color32::from_rgb(15, 23, 42),
                );
                let max_i = self
                    .morrison_obs
                    .iter()
                    .map(|(x, _)| *x)
                    .fold(50.0_f64, f64::max)
                    * 1.06;
                let max_v0 = self
                    .morrison_obs
                    .iter()
                    .map(|(_, v)| *v)
                    .fold(10.0_f64, f64::max)
                    .max(self.morrison_fit.v0_unhibited)
                    * 1.1;
                let map_m = |inh: f64, v: f64| -> egui::Pos2 {
                    let fx = (inh / max_i).clamp(0.0, 1.0) as f32;
                    let fy = (v / max_v0).clamp(0.0, 1.0) as f32;
                    egui::pos2(
                        rect.left() + fx * rect.width(),
                        rect.bottom() - fy * rect.height(),
                    )
                };
                for k in 0..=4 {
                    let vv = max_v0 * (k as f64 / 4.0);
                    let pt = map_m(0.0, vv);
                    p.line_segment(
                        [egui::pos2(rect.left(), pt.y), egui::pos2(rect.right(), pt.y)],
                        grid_s,
                    );
                    p.text(
                        egui::pos2(rect.left() - 6.0, pt.y),
                        egui::Align2::RIGHT_CENTER,
                        format!("{vv:.1}"),
                        egui::FontId::proportional(9.5),
                        axis_c,
                    );
                }
                for k in 0..=4 {
                    let iv = max_i * (k as f64 / 4.0);
                    let pt = map_m(iv, 0.0);
                    p.line_segment(
                        [egui::pos2(pt.x, rect.top()), egui::pos2(pt.x, rect.bottom())],
                        grid_s,
                    );
                    p.text(
                        egui::pos2(pt.x, rect.bottom() + 4.0),
                        egui::Align2::CENTER_TOP,
                        format!("{iv:.0} nM"),
                        egui::FontId::proportional(9.5),
                        axis_c,
                    );
                }
                p.rect_stroke(
                    rect,
                    0.0,
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
                );
                let curve: Vec<egui::Pos2> = (0..=80)
                    .map(|s| {
                        let inh = max_i * (s as f64 / 80.0);
                        let v = eval_morrison_velocity(
                            inh,
                            self.morrison_fit.v0_unhibited,
                            self.morrison_fit.active_enzyme_et_nm,
                            self.morrison_fit.ki_app_nm,
                        );
                        map_m(inh, v)
                    })
                    .collect();
                p.add(egui::Shape::line(
                    curve,
                    egui::Stroke::new(2.4_f32, egui::Color32::from_rgb(124, 58, 237)),
                ));
                for &(ix, vx) in &self.morrison_obs {
                    p.circle_filled(map_m(ix, vx), 3.8, egui::Color32::from_rgb(220, 38, 38));
                }
                let et_x = map_m(self.morrison_fit.active_enzyme_et_nm, 0.0).x;
                p.line_segment(
                    [egui::pos2(et_x, rect.top()), egui::pos2(et_x, rect.bottom())],
                    egui::Stroke::new(1.2_f32, egui::Color32::from_rgb(5, 150, 105)),
                );
                p.text(
                    egui::pos2(rect.right() - 8.0, rect.top() + 8.0),
                    egui::Align2::RIGHT_TOP,
                    format!(
                        "[E]_T = {:.2} nM | Ki_app = {:.2} nM | R^2 = {:.4}",
                        self.morrison_fit.active_enzyme_et_nm,
                        self.morrison_fit.ki_app_nm,
                        self.morrison_fit.r_squared
                    ),
                    egui::FontId::proportional(10.5),
                    egui::Color32::from_rgb(15, 23, 42),
                );
            }

            cols[0].add_space(6.0);
            cols[0].heading("Global 2D AICc Inhibition Model Selection");
            cols[0].separator();
            for fit in &self.kin_fits {
                cols[0].label(format!(
                    "{:<34} | Vmax={:.1} Km={:.1}uM Ki={:.1}nM | dAICc={:.1} (w={:.3}) R^2={:.4}",
                    fit.mechanism.label(),
                    fit.vmax,
                    fit.km_um,
                    fit.ki_nm,
                    fit.delta_aicc,
                    fit.akaike_weight,
                    fit.r_squared
                ));
            }

            // Right Column: Circular Plasmid Vector Map + Virtual Agarose Gel Digest
            let w1 = cols[1].available_width().max(420.0);
            let viz_h = top_h + bot_h + 6.0;
            let (resp_p, painter_p) =
                cols[1].allocate_painter(egui::vec2(w1, viz_h), egui::Sense::hover());
            let outer_p = resp_p.rect;
            painter_p.rect_filled(outer_p, 6.0, egui::Color32::from_rgb(248, 250, 252));
            painter_p.rect_stroke(
                outer_p,
                6.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
            );
            let pp = painter_p.with_clip_rect(outer_p);
            let map_w = outer_p.width() * 0.59;
            let map_rect = egui::Rect::from_min_max(
                outer_p.left_top(),
                egui::pos2(outer_p.left() + map_w, outer_p.bottom()),
            );
            let gel_rect = egui::Rect::from_min_max(
                egui::pos2(map_rect.right() + 6.0, outer_p.top() + 34.0),
                egui::pos2(outer_p.right() - 12.0, outer_p.bottom() - 16.0),
            );

            // Plasmid Circular Map
            pp.text(
                egui::pos2(map_rect.center().x, outer_p.top() + 10.0),
                egui::Align2::CENTER_TOP,
                format!(
                    "Circular Expression Vector: {} ({} bp)",
                    self.plasmid.name,
                    self.plasmid.sequence.len()
                ),
                egui::FontId::proportional(13.0),
                egui::Color32::from_rgb(15, 23, 42),
            );
            let center = map_rect.center() + egui::vec2(0.0, 24.0);
            let radius = (map_rect.width().min(map_rect.height()) * 0.255).clamp(85.0, 134.0);
            let total_bp = self.plasmid.sequence.len().max(1) as f64;
            pp.circle_stroke(
                center,
                radius,
                egui::Stroke::new(2.0_f32, egui::Color32::from_rgb(100, 116, 139)),
            );
            let bp_to_angle = |bp: f64| -> f64 {
                -std::f64::consts::FRAC_PI_2 + (bp / total_bp) * std::f64::consts::TAU
            };
            for tick_bp in (0..self.plasmid.sequence.len()).step_by(100) {
                let a = bp_to_angle(tick_bp as f64);
                let p0 = center + egui::vec2(a.cos() as f32, a.sin() as f32) * (radius - 6.0);
                let p1 = center + egui::vec2(a.cos() as f32, a.sin() as f32) * (radius + 6.0);
                pp.line_segment(
                    [p0, p1],
                    egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(148, 163, 184)),
                );
                let pt_lbl = center + egui::vec2(a.cos() as f32, a.sin() as f32) * (radius - 20.0);
                pp.text(
                    pt_lbl,
                    egui::Align2::CENTER_CENTER,
                    format!("{tick_bp}"),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(100, 116, 139),
                );
            }
            for (fi, feat) in self.plasmid.features.iter().enumerate() {
                let fcol = parse_hex_color(&feat.color_hex);
                let r_arc = radius + if fi % 2 == 0 { 10.0 } else { 19.0 };
                let end_bp = if feat.end_0 >= feat.start_0 {
                    feat.end_0 as f64
                } else {
                    (self.plasmid.sequence.len() + feat.end_0) as f64
                };
                let steps = 28;
                let arc_pts: Vec<egui::Pos2> = (0..=steps)
                    .map(|k| {
                        let bp = feat.start_0 as f64
                            + (end_bp - feat.start_0 as f64) * (k as f64 / steps as f64);
                        let a = bp_to_angle(bp);
                        center + egui::vec2(a.cos() as f32, a.sin() as f32) * r_arc
                    })
                    .collect();
                pp.add(egui::Shape::line(arc_pts, egui::Stroke::new(6.0_f32, fcol)));
                if end_bp - feat.start_0 as f64 >= 70.0 && feat.name != "Minimal_ori" {
                    let mid_a = bp_to_angle((feat.start_0 as f64 + end_bp) * 0.5);
                    let lbl_pos =
                        center + egui::vec2(mid_a.cos() as f32, mid_a.sin() as f32) * (r_arc + 20.0);
                    pp.text(
                        lbl_pos,
                        egui::Align2::CENTER_CENTER,
                        &feat.name,
                        egui::FontId::proportional(9.8),
                        egui::Color32::from_rgb(15, 23, 42),
                    );
                }
            }
            // Feature Key in Top-Left of map_rect so short elements never collide at 12 o'clock
            for (fi, feat) in self.plasmid.features.iter().enumerate() {
                let fcol = parse_hex_color(&feat.color_hex);
                let fy = outer_p.top() + 38.0 + fi as f32 * 14.0;
                pp.rect_filled(
                    egui::Rect::from_center_size(
                        egui::pos2(outer_p.left() + 18.0, fy),
                        egui::vec2(9.0, 9.0),
                    ),
                    2.0,
                    fcol,
                );
                pp.text(
                    egui::pos2(outer_p.left() + 27.0, fy),
                    egui::Align2::LEFT_CENTER,
                    format!("{} ({} bp)", feat.name, feat.span_bp(self.plasmid.sequence.len())),
                    egui::FontId::proportional(9.5),
                    egui::Color32::from_rgb(51, 65, 85),
                );
            }
            let mut prev_disp_a = -10.0_f64;
            for (ci, cut) in self.unique_cuts.iter().enumerate() {
                let a = bp_to_angle(cut.cut_top_0 as f64);
                let is_used = self.digest_enzymes.iter().any(|e| e == &cut.enzyme_name);
                let ccol = if is_used {
                    egui::Color32::from_rgb(220, 38, 38)
                } else {
                    egui::Color32::from_rgb(71, 85, 105)
                };
                let p0 = center + egui::vec2(a.cos() as f32, a.sin() as f32) * (radius - 2.0);
                let p1 = center + egui::vec2(a.cos() as f32, a.sin() as f32) * (radius + 24.0);
                pp.line_segment(
                    [p0, p1],
                    egui::Stroke::new(if is_used { 2.0_f32 } else { 1.0_f32 }, ccol),
                );
                let disp_a = if ci == 0 {
                    a
                } else {
                    a.max(prev_disp_a + 0.24)
                };
                prev_disp_a = disp_a;
                let moved = (disp_a - a).abs() > 0.02;
                let radial_out = if ci % 2 == 1 && moved {
                    radius + 54.0
                } else {
                    radius + 42.0
                };
                let lp = center + egui::vec2(disp_a.cos() as f32, disp_a.sin() as f32) * radial_out;
                if moved {
                    pp.line_segment([p1, lp], egui::Stroke::new(0.8_f32, ccol));
                }
                let align = if disp_a.cos() > 0.25 {
                    egui::Align2::LEFT_CENTER
                } else if disp_a.cos() < -0.25 {
                    egui::Align2::RIGHT_CENTER
                } else {
                    egui::Align2::CENTER_CENTER
                };
                pp.text(
                    lp,
                    align,
                    format!("{} ({})", cut.enzyme_name, cut.cut_top_0),
                    egui::FontId::monospace(8.5),
                    ccol,
                );
            }
            let gc_count = self
                .plasmid
                .sequence
                .bytes()
                .filter(|b| matches!(b, b'G' | b'C' | b'g' | b'c'))
                .count();
            let gc_pct = 100.0 * gc_count as f64 / total_bp;
            pp.text(
                center + egui::vec2(0.0, -8.0),
                egui::Align2::CENTER_CENTER,
                &self.plasmid.name,
                egui::FontId::proportional(13.0),
                egui::Color32::from_rgb(15, 23, 42),
            );
            pp.text(
                center + egui::vec2(0.0, 10.0),
                egui::Align2::CENTER_CENTER,
                format!("{} bp | GC {:.1}%", self.plasmid.sequence.len(), gc_pct),
                egui::FontId::monospace(10.0),
                egui::Color32::from_rgb(71, 85, 105),
            );

            // Virtual 1.2% Agarose Gel Electrophoresis
            pp.text(
                egui::pos2(gel_rect.center().x, outer_p.top() + 10.0),
                egui::Align2::CENTER_TOP,
                "Virtual 1.2% Agarose Gel Digest",
                egui::FontId::proportional(13.0),
                egui::Color32::from_rgb(15, 23, 42),
            );
            pp.rect_filled(gel_rect, 5.0, egui::Color32::from_rgb(15, 23, 42));
            pp.rect_stroke(
                gel_rect,
                5.0,
                egui::Stroke::new(1.0_f32, egui::Color32::from_rgb(51, 65, 85)),
            );
            let gp = pp.with_clip_rect(gel_rect);
            let lanes_x = [
                gel_rect.left() + gel_rect.width() * 0.25,
                gel_rect.left() + gel_rect.width() * 0.54,
                gel_rect.left() + gel_rect.width() * 0.82,
            ];
            let lane_titles = ["Ladder", "Digest", "Uncut"];
            for (&lx, lt) in lanes_x.iter().zip(lane_titles.iter()) {
                gp.text(
                    egui::pos2(lx, gel_rect.top() + 8.0),
                    egui::Align2::CENTER_TOP,
                    *lt,
                    egui::FontId::monospace(9.5),
                    egui::Color32::from_rgb(148, 163, 184),
                );
                gp.rect_filled(
                    egui::Rect::from_center_size(
                        egui::pos2(lx, gel_rect.top() + 28.0),
                        egui::vec2(34.0, 7.0),
                    ),
                    1.0,
                    egui::Color32::from_rgb(30, 41, 59),
                );
            }
            let bp_to_gel_y = |bp: usize| -> f32 {
                let lmin = (15.0_f64).log10();
                let lmax = (1300.0_f64).log10();
                let lbp = (bp.max(15) as f64).log10().clamp(lmin, lmax);
                let frac = ((lmax - lbp) / (lmax - lmin)) as f32;
                gel_rect.top() + 48.0 + frac * (gel_rect.height() - 68.0)
            };
            for &band_bp in &[1000usize, 800, 600, 500, 400, 300, 200, 100] {
                let by = bp_to_gel_y(band_bp);
                let col = if band_bp == 500 {
                    egui::Color32::from_rgb(56, 189, 248)
                } else {
                    egui::Color32::from_rgb(148, 163, 184)
                };
                gp.rect_filled(
                    egui::Rect::from_center_size(
                        egui::pos2(lanes_x[0], by),
                        egui::vec2(32.0, if band_bp == 500 { 4.2 } else { 2.8 }),
                    ),
                    1.0,
                    col,
                );
                gp.text(
                    egui::pos2(lanes_x[0] - 20.0, by),
                    egui::Align2::RIGHT_CENTER,
                    format!("{band_bp}"),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(148, 163, 184),
                );
            }
            for frag in &self.digest_frags {
                let by = bp_to_gel_y(frag.length_bp);
                gp.rect_filled(
                    egui::Rect::from_center_size(egui::pos2(lanes_x[1], by), egui::vec2(36.0, 4.5)),
                    1.0,
                    egui::Color32::from_rgb(74, 222, 128),
                );
                gp.text(
                    egui::pos2(lanes_x[1] + 21.0, by),
                    egui::Align2::LEFT_CENTER,
                    format!("{}bp", frag.length_bp),
                    egui::FontId::monospace(8.5),
                    egui::Color32::from_rgb(187, 247, 208),
                );
            }
            let uncut_y = bp_to_gel_y(self.plasmid.sequence.len());
            gp.rect_filled(
                egui::Rect::from_center_size(egui::pos2(lanes_x[2], uncut_y), egui::vec2(34.0, 4.8)),
                1.0,
                egui::Color32::from_rgb(250, 204, 21),
            );
            gp.text(
                egui::pos2(lanes_x[2], uncut_y + 10.0),
                egui::Align2::CENTER_TOP,
                format!("{}bp", self.plasmid.sequence.len()),
                egui::FontId::monospace(8.0),
                egui::Color32::from_rgb(253, 224, 71),
            );

            cols[1].add_space(6.0);
            cols[1].heading(format!(
                "Plasmid Features & Restriction Digest ({})",
                self.digest_enzymes.join(" + ")
            ));
            cols[1].separator();
            let feat_line: Vec<String> = self
                .plasmid
                .features
                .iter()
                .map(|f| format!("{} [{}..{}]", f.name, f.start_0, f.end_0))
                .collect();
            cols[1].label(format!("Features : {}", feat_line.join(" | ")));
            let frag_line: Vec<String> = self
                .digest_frags
                .iter()
                .enumerate()
                .map(|(fi, f)| format!("#{} {} bp ({}-{})", fi + 1, f.length_bp, f.left_enzyme, f.right_enzyme))
                .collect();
            cols[1].label(format!("Fragments: {}", frag_line.join(" | ")));
            cols[1].label(format!(
                "Unique Cutters Found: {}",
                self.unique_cuts
                    .iter()
                    .map(|c| format!("{}@{}", c.enzyme_name, c.cut_top_0))
                    .collect::<Vec<_>>()
                    .join(", ")
            ));
        });
    }
}

fn probit_norm_cdf(z: f64) -> f64 {
    let sign = if z < 0.0 { -1.0 } else { 1.0 };
    let x = z.abs() / std::f64::consts::SQRT_2;
    let t = 1.0 / (1.0 + 0.3275911 * x);
    let y = 1.0
        - (((((1.061405429 * t - 1.453152027) * t + 1.421413741) * t - 0.284496736) * t
            + 0.254829592)
            * t)
            * (-x * x).exp();
    0.5 * (1.0 + sign * y)
}

fn parse_hex_color(hex: &str) -> egui::Color32 {
    let s = hex.trim_start_matches('#');
    if s.len() == 6 {
        if let (Ok(r), Ok(g), Ok(b)) = (
            u8::from_str_radix(&s[0..2], 16),
            u8::from_str_radix(&s[2..4], 16),
            u8::from_str_radix(&s[4..6], 16),
        ) {
            return egui::Color32::from_rgb(r, g, b);
        }
    }
    egui::Color32::from_rgb(37, 99, 235)
}

fn block_on_wgpu<F: std::future::Future>(fut: F) -> F::Output {
    struct ThreadWaker(std::thread::Thread);
    impl std::task::Wake for ThreadWaker {
        fn wake(self: std::sync::Arc<Self>) {
            self.0.unpark();
        }
        fn wake_by_ref(self: &std::sync::Arc<Self>) {
            self.0.unpark();
        }
    }
    let waker = std::task::Waker::from(std::sync::Arc::new(ThreadWaker(std::thread::current())));
    let mut fut = std::pin::pin!(fut);
    let mut cx = std::task::Context::from_waker(&waker);
    loop {
        match fut.as_mut().poll(&mut cx) {
            std::task::Poll::Ready(out) => return out,
            std::task::Poll::Pending => std::thread::park(),
        }
    }
}

pub(crate) fn render_flux_offscreen_png(
    ctx: &egui::Context,
    out_path: &std::path::Path,
    mut draw_ui: impl FnMut(&egui::Context),
) -> Result<(), String> {
    use eframe::{egui_wgpu, wgpu};
    const WIDTH: u32 = 1600;
    const HEIGHT: u32 = 1000;
    const BYTES_PER_ROW: u32 = WIDTH * 4;

    let instance = wgpu::Instance::new(wgpu::InstanceDescriptor {
        backends: wgpu::Backends::PRIMARY | wgpu::Backends::GL,
        ..Default::default()
    });
    let adapter = block_on_wgpu(instance.request_adapter(&wgpu::RequestAdapterOptions {
        power_preference: wgpu::PowerPreference::HighPerformance,
        compatible_surface: None,
        force_fallback_adapter: false,
    }))
    .or_else(|| {
        block_on_wgpu(instance.request_adapter(&wgpu::RequestAdapterOptions {
            power_preference: wgpu::PowerPreference::None,
            compatible_surface: None,
            force_fallback_adapter: true,
        }))
    })
    .ok_or_else(|| "No wgpu adapter available for offscreen screenshot".to_string())?;

    let mut limits = wgpu::Limits::default();
    limits.max_texture_dimension_2d = adapter.limits().max_texture_dimension_2d.min(8192);
    let (device, queue) = block_on_wgpu(adapter.request_device(
        &wgpu::DeviceDescriptor {
            label: Some("osabt_offscreen_device"),
            required_features: wgpu::Features::empty(),
            required_limits: limits,
            memory_hints: Default::default(),
        },
        None,
    ))
    .map_err(|e| format!("wgpu request_device failed: {e}"))?;

    let format = wgpu::TextureFormat::Rgba8Unorm;
    let texture = device.create_texture(&wgpu::TextureDescriptor {
        label: Some("osabt_offscreen_target"),
        size: wgpu::Extent3d {
            width: WIDTH,
            height: HEIGHT,
            depth_or_array_layers: 1,
        },
        mip_level_count: 1,
        sample_count: 1,
        dimension: wgpu::TextureDimension::D2,
        format,
        usage: wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::COPY_SRC,
        view_formats: &[],
    });
    let view = texture.create_view(&wgpu::TextureViewDescriptor::default());
    let staging = device.create_buffer(&wgpu::BufferDescriptor {
        label: Some("osabt_offscreen_staging"),
        size: u64::from(BYTES_PER_ROW) * u64::from(HEIGHT),
        usage: wgpu::BufferUsages::COPY_DST | wgpu::BufferUsages::MAP_READ,
        mapped_at_creation: false,
    });

    let mut renderer = egui_wgpu::Renderer::new(&device, format, None, 1, false);
    let screen_desc = egui_wgpu::ScreenDescriptor {
        size_in_pixels: [WIDTH, HEIGHT],
        pixels_per_point: 1.0,
    };

    for frame_idx in 0..4usize {
        let mut raw_input = egui::RawInput::default();
        raw_input.screen_rect = Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(WIDTH as f32, HEIGHT as f32),
        ));
        raw_input
            .viewports
            .entry(egui::ViewportId::ROOT)
            .and_modify(|v| v.native_pixels_per_point = Some(1.0));
        raw_input.time = Some(frame_idx as f64 / 60.0);
        ctx.set_pixels_per_point(1.0);

        let full_output = ctx.run(raw_input, &mut draw_ui);
        let paint_jobs = ctx.tessellate(full_output.shapes, full_output.pixels_per_point);

        for (id, delta) in &full_output.textures_delta.set {
            renderer.update_texture(&device, &queue, *id, delta);
        }

        let mut prep_encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("osabt_offscreen_prep"),
        });
        let user_cmds =
            renderer.update_buffers(&device, &queue, &mut prep_encoder, &paint_jobs, &screen_desc);
        queue.submit(user_cmds.into_iter().chain(std::iter::once(prep_encoder.finish())));

        let mut pass_encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("osabt_offscreen_pass_encoder"),
        });
        {
            let mut rpass = pass_encoder
                .begin_render_pass(&wgpu::RenderPassDescriptor {
                    label: Some("osabt_offscreen_pass"),
                    color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                        view: &view,
                        resolve_target: None,
                        ops: wgpu::Operations {
                            load: wgpu::LoadOp::Clear(wgpu::Color {
                                r: 0.07,
                                g: 0.08,
                                b: 0.11,
                                a: 1.0,
                            }),
                            store: wgpu::StoreOp::Store,
                        },
                    })],
                    depth_stencil_attachment: None,
                    timestamp_writes: None,
                    occlusion_query_set: None,
                })
                .forget_lifetime();
            renderer.render(&mut rpass, &paint_jobs, &screen_desc);
        }

        if frame_idx == 3 {
            pass_encoder.copy_texture_to_buffer(
                wgpu::ImageCopyTexture {
                    texture: &texture,
                    mip_level: 0,
                    origin: wgpu::Origin3d::ZERO,
                    aspect: wgpu::TextureAspect::All,
                },
                wgpu::ImageCopyBuffer {
                    buffer: &staging,
                    layout: wgpu::ImageDataLayout {
                        offset: 0,
                        bytes_per_row: Some(BYTES_PER_ROW),
                        rows_per_image: Some(HEIGHT),
                    },
                },
                wgpu::Extent3d {
                    width: WIDTH,
                    height: HEIGHT,
                    depth_or_array_layers: 1,
                },
            );
        }

        queue.submit(std::iter::once(pass_encoder.finish()));
        for id in &full_output.textures_delta.free {
            renderer.free_texture(id);
        }
        device.poll(wgpu::Maintain::Wait);
    }

    let slice = staging.slice(..);
    let (tx, rx) = std::sync::mpsc::sync_channel(1);
    slice.map_async(wgpu::MapMode::Read, move |res| {
        let _ = tx.send(res);
    });
    device.poll(wgpu::Maintain::Wait);
    rx.recv()
        .map_err(|e| format!("GPU readback channel failed: {e}"))?
        .map_err(|e| format!("GPU buffer map failed: {e}"))?;

    let rgba = slice.get_mapped_range().to_vec();
    staging.unmap();
    write_rgba_png(out_path, WIDTH, HEIGHT, &rgba)
}

fn write_rgba_png(path: &std::path::Path, width: u32, height: u32, rgba: &[u8]) -> Result<(), String> {
    if rgba.len() != (width as usize) * (height as usize) * 4 {
        return Err("RGBA buffer length mismatch".to_string());
    }
    if let Some(parent) = path.parent() {
        if !parent.as_os_str().is_empty() {
            fs::create_dir_all(parent).map_err(|e| format!("mkdir failed: {e}"))?;
        }
    }
    let file = fs::File::create(path)
        .map_err(|e| format!("Failed creating PNG '{}': {e}", path.display()))?;
    let writer = BufWriter::new(file);
    let mut encoder = png::Encoder::new(writer, width, height);
    encoder.set_color(png::ColorType::Rgba);
    encoder.set_depth(png::BitDepth::Eight);
    encoder.set_compression(png::Compression::High);
    let mut png_writer = encoder
        .write_header()
        .map_err(|e| format!("Failed writing PNG header '{}': {e}", path.display()))?;
    png_writer
        .write_image_data(rgba)
        .map_err(|e| format!("Failed writing PNG data '{}': {e}", path.display()))?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cycles_all_flux_views_headlessly() {
        let ctx = egui::Context::default();
        let mut app = FluxGuiApp::new(None, None).unwrap();
        for view in [
            FluxView::DualOverview,
            FluxView::LeafAciFvcb,
            FluxView::StreamDielDo,
            FluxView::SurvivalAndEcotox,
            FluxView::KineticsAndPlasmid,
        ] {
            app.view = view;
            for frame in 0..2 {
                let mut raw = egui::RawInput::default();
                raw.screen_rect = Some(egui::Rect::from_min_size(egui::Pos2::ZERO, egui::vec2(1280.0, 800.0)));
                raw.time = Some(frame as f64 / 60.0);
                let out = ctx.run(raw, |c| app.draw_panels(c));
                assert!(!out.shapes.is_empty());
            }
        }
    }
}
