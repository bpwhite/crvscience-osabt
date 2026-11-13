// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #13 (`osabt kinetics`) — Cross-platform `eframe` / `egui` desktop studio for
//! global 2D reversible inhibition mechanism selection (`AICc` + Eadie–Hofstee diagnostics),
//! Morrison (1969) tight-binding active-site & `Ki` regression, and single-substrate kinetics.

use std::collections::BTreeMap;
use std::fs;
use std::path::PathBuf;

use eframe::egui;
use egui::{Color32, Pos2, Sense, Shape, Stroke, Vec2};

use crate::kinetics_core::{
    convert_ic50_cheng_prusoff, demo_competitive_kinase_matrix_csv, demo_morrison_tight_binding_csv,
    eval_inhibition_velocity, eval_morrison_velocity, eval_single_substrate_velocity,
    fit_global_reversible_inhibition, fit_morrison_tight_binding, fit_single_substrate_models,
    parse_kinetics_csv, GlobalInhibitionFit, InhibitionMechanism, KineticsPoint,
    MorrisonTightBindingFit, SingleSubstrateFit,
};
use crate::kinetics_svg::{render_global_inhibition_svg, render_morrison_svg};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum KineticsMode {
    GlobalInhibition,
    MorrisonTightBinding,
    SingleSubstrateAndChengPrusoff,
}

/// Interactive `egui` application state for `OSABT Kinetics Studio`.
pub struct KineticsStudioApp {
    mode: KineticsMode,
    csv_path_input: String,
    matrix_points: Vec<KineticsPoint>,
    inh_fits: Vec<GlobalInhibitionFit>,
    selected_mech_idx: usize,
    morrison_obs: Vec<(f64, f64)>,
    nominal_et_nm: f64,
    morrison_s_um: f64,
    morrison_km_um: f64,
    fix_morrison_et: bool,
    morrison_fit: MorrisonTightBindingFit,
    enzyme_et_for_kcat_nm: f64,
    single_fits: Vec<SingleSubstrateFit>,
    cp_mech: InhibitionMechanism,
    cp_ic50_nm: f64,
    cp_s_um: f64,
    cp_km_um: f64,
    status_bar: String,
}

impl KineticsStudioApp {
    pub fn new(initial_csv: Option<PathBuf>) -> Self {
        let matrix_points =
            parse_kinetics_csv(demo_competitive_kinase_matrix_csv()).expect("valid demo matrix");
        let inh_fits =
            fit_global_reversible_inhibition(&matrix_points).expect("valid global inhibition fit");
        let single_fits =
            fit_single_substrate_models(&matrix_points, Some(10.0)).expect("valid MM fit");

        let raw_mor =
            parse_kinetics_csv(demo_morrison_tight_binding_csv()).expect("valid Morrison CSV");
        let morrison_obs: Vec<(f64, f64)> = raw_mor
            .iter()
            .map(|p| (p.substrate_um, p.velocity))
            .collect();
        let morrison_fit =
            fit_morrison_tight_binding(&morrison_obs, None, Some(25.0), Some((70.0, 35.0)))
                .expect("valid Morrison fit");

        let mut app = Self {
            mode: KineticsMode::GlobalInhibition,
            csv_path_input: initial_csv
                .as_ref()
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_default(),
            matrix_points,
            inh_fits,
            selected_mech_idx: 0,
            morrison_obs,
            nominal_et_nm: 25.0,
            morrison_s_um: 70.0,
            morrison_km_um: 35.0,
            fix_morrison_et: false,
            morrison_fit,
            enzyme_et_for_kcat_nm: 10.0,
            single_fits,
            cp_mech: InhibitionMechanism::Competitive,
            cp_ic50_nm: 66.0,
            cp_s_um: 70.0,
            cp_km_um: 35.0,
            status_bar: String::from(
                "Global AICc selects Competitive inhibition (Ki = 22 nM) | Switch tabs for Morrison Tight-Binding & Single-Substrate fits.",
            ),
        };

        if let Some(p) = initial_csv {
            app.load_csv_file(&p);
        }
        app
    }

    fn load_csv_file(&mut self, path: &PathBuf) {
        match fs::read_to_string(path).and_then(|s| {
            parse_kinetics_csv(&s)
                .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e))
        }) {
            Ok(pts) => {
                let has_inh = pts.iter().any(|p| p.inhibitor_nm > 0.0);
                if has_inh {
                    if let Ok(fits) = fit_global_reversible_inhibition(&pts) {
                        self.inh_fits = fits;
                        self.selected_mech_idx = 0;
                    }
                    if let Ok(sf) =
                        fit_single_substrate_models(&pts, Some(self.enzyme_et_for_kcat_nm))
                    {
                        self.single_fits = sf;
                    }
                    self.matrix_points = pts;
                    self.mode = KineticsMode::GlobalInhibition;
                    self.status_bar = format!("Loaded 2D kinetics matrix '{}'.", path.display());
                } else {
                    let obs: Vec<(f64, f64)> =
                        pts.iter().map(|p| (p.substrate_um, p.velocity)).collect();
                    if let Ok(mf) = fit_morrison_tight_binding(
                        &obs,
                        if self.fix_morrison_et {
                            Some(self.nominal_et_nm)
                        } else {
                            None
                        },
                        Some(self.nominal_et_nm),
                        Some((self.morrison_s_um, self.morrison_km_um)),
                    ) {
                        self.morrison_obs = obs;
                        self.morrison_fit = mf;
                    }
                    if let Ok(sf) =
                        fit_single_substrate_models(&pts, Some(self.enzyme_et_for_kcat_nm))
                    {
                        self.single_fits = sf;
                    }
                    self.status_bar = format!("Loaded 2-column kinetics curve '{}'.", path.display());
                }
            }
            Err(e) => self.status_bar = format!("CSV load error: {e}"),
        }
    }

    fn refit_morrison(&mut self) {
        let fixed = if self.fix_morrison_et {
            Some(self.nominal_et_nm)
        } else {
            None
        };
        if let Ok(mf) = fit_morrison_tight_binding(
            &self.morrison_obs,
            fixed,
            Some(self.nominal_et_nm),
            Some((self.morrison_s_um, self.morrison_km_um)),
        ) {
            self.morrison_fit = mf;
        }
    }
}

/// Launches the cross-platform `OSABT Kinetics Studio` desktop GUI.
pub fn run_kinetics_gui(initial_csv: Option<PathBuf>) -> Result<(), String> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1240.0, 780.0])
            .with_title(
                "OSABT Kinetics — Global Ki Mode Selection, Eadie-Hofstee & Morrison Tight-Binding Studio (GPL-3.0)",
            ),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-kinetics-studio",
        options,
        Box::new(move |_cc| Ok(Box::new(KineticsStudioApp::new(initial_csv)))),
    )
    .map_err(|e| format!("Failed launching Kinetics GUI: {e}"))
}

impl eframe::App for KineticsStudioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        egui::TopBottomPanel::top("kinetics_top").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading("⚗ OSABT Kinetics");
                ui.separator();
                if ui
                    .selectable_label(
                        self.mode == KineticsMode::GlobalInhibition,
                        "1. Global Reversible Ki (AICc)",
                    )
                    .clicked()
                {
                    self.mode = KineticsMode::GlobalInhibition;
                }
                if ui
                    .selectable_label(
                        self.mode == KineticsMode::MorrisonTightBinding,
                        "2. Morrison Tight-Binding ([E]_T & Ki)",
                    )
                    .clicked()
                {
                    self.mode = KineticsMode::MorrisonTightBinding;
                }
                if ui
                    .selectable_label(
                        self.mode == KineticsMode::SingleSubstrateAndChengPrusoff,
                        "3. MM / Hill / Cheng–Prusoff",
                    )
                    .clicked()
                {
                    self.mode = KineticsMode::SingleSubstrateAndChengPrusoff;
                }
                ui.separator();
                ui.label("CSV:");
                ui.add(egui::TextEdit::singleline(&mut self.csv_path_input).desired_width(185.0));
                if ui.button("Load").clicked() && !self.csv_path_input.trim().is_empty() {
                    let p = PathBuf::from(self.csv_path_input.trim());
                    self.load_csv_file(&p);
                }
                if ui.button("💾 Export SVGs").clicked() {
                    let g_svg = render_global_inhibition_svg(
                        &self.matrix_points,
                        &self.inh_fits[self.selected_mech_idx],
                        "Global Kinase Inhibition",
                    );
                    let m_svg = render_morrison_svg(
                        &self.morrison_obs,
                        &self.morrison_fit,
                        "Morrison Tight-Binding Titration",
                    );
                    if fs::write("kinetics_global_inhibition.svg", g_svg).is_ok()
                        && fs::write("kinetics_morrison_tight_binding.svg", m_svg).is_ok()
                    {
                        self.status_bar = String::from(
                            "Exported kinetics_global_inhibition.svg and kinetics_morrison_tight_binding.svg!",
                        );
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("kinetics_footer").show(ctx, |ui| {
            ui.label(&self.status_bar);
        });

        egui::CentralPanel::default().show(ctx, |ui| match self.mode {
            KineticsMode::GlobalInhibition => self.draw_global_inhibition_view(ui),
            KineticsMode::MorrisonTightBinding => self.draw_morrison_view(ui),
            KineticsMode::SingleSubstrateAndChengPrusoff => self.draw_single_and_cp_view(ui),
        });
    }
}

impl KineticsStudioApp {
    fn draw_global_inhibition_view(&mut self, ui: &mut egui::Ui) {
        let canvas_h = (ui.available_height() * 0.56).max(220.0);
        let (resp, painter) =
            ui.allocate_painter(Vec2::new(ui.available_width(), canvas_h), Sense::hover());
        let rect = resp.rect;
        painter.rect_filled(rect, 4.0, Color32::from_rgb(15, 23, 42));

        let fit = &self.inh_fits[self.selected_mech_idx];
        let pad = 42.0_f32;
        let mid_x = rect.center().x;
        let left_box = egui::Rect::from_min_max(
            Pos2::new(rect.left() + pad, rect.top() + 24.0),
            Pos2::new(mid_x - 24.0, rect.bottom() - 30.0),
        );
        let right_box = egui::Rect::from_min_max(
            Pos2::new(mid_x + 24.0, rect.top() + 24.0),
            Pos2::new(rect.right() - 24.0, rect.bottom() - 30.0),
        );
        painter.rect_stroke(left_box, 2.0, Stroke::new(1.0, Color32::GRAY));
        painter.rect_stroke(right_box, 2.0, Stroke::new(1.0, Color32::GRAY));

        let s_max = self
            .matrix_points
            .iter()
            .map(|p| p.substrate_um)
            .fold(100.0_f64, f64::max)
            * 1.05;
        let v_max = (fit.vmax * 1.08).max(10.0);
        let eh_max = self
            .matrix_points
            .iter()
            .filter(|p| p.substrate_um > 0.0)
            .map(|p| p.velocity / p.substrate_um)
            .fold(0.5_f64, f64::max)
            * 1.12;

        let mut by_i: BTreeMap<i64, Vec<KineticsPoint>> = BTreeMap::new();
        for &pt in &self.matrix_points {
            by_i.entry((pt.inhibitor_nm * 1000.0).round() as i64)
                .or_default()
                .push(pt);
        }

        let palette = [
            Color32::from_rgb(56, 189, 248),
            Color32::from_rgb(74, 222, 128),
            Color32::from_rgb(250, 204, 21),
            Color32::from_rgb(248, 113, 113),
        ];

        for (c_idx, (&ik, grp)) in by_i.iter().enumerate() {
            let inh_nm = (ik as f64) / 1000.0;
            let col = palette[c_idx % palette.len()];
            let mut mm_line = Vec::with_capacity(70);
            let mut eh_line = Vec::with_capacity(70);
            for k in 1..=70 {
                let s = (k as f64 / 70.0) * s_max;
                let vp = eval_inhibition_velocity(
                    fit.mechanism,
                    s,
                    inh_nm,
                    fit.vmax,
                    fit.km_um,
                    fit.ki_nm,
                    fit.alpha_factor.unwrap_or(1.0),
                );
                let px1 = left_box.left() + ((s / s_max) as f32) * left_box.width();
                let py = left_box.bottom() - ((vp / v_max).clamp(0.0, 1.0) as f32) * left_box.height();
                mm_line.push(Pos2::new(px1, py));

                let px2 =
                    right_box.left() + (((vp / s) / eh_max).clamp(0.0, 1.0) as f32) * right_box.width();
                eh_line.push(Pos2::new(px2, py));
            }
            painter.add(Shape::line(mm_line, Stroke::new(2.1, col)));
            painter.add(Shape::line(eh_line, Stroke::new(1.6, col)));

            for pt in grp {
                let py =
                    left_box.bottom() - ((pt.velocity / v_max).clamp(0.0, 1.0) as f32) * left_box.height();
                let px1 = left_box.left() + ((pt.substrate_um / s_max) as f32) * left_box.width();
                painter.circle_filled(Pos2::new(px1, py), 3.8, col);
                if pt.substrate_um > 0.0 {
                    let px2 = right_box.left()
                        + (((pt.velocity / pt.substrate_um) / eh_max).clamp(0.0, 1.0) as f32)
                            * right_box.width();
                    painter.circle_filled(Pos2::new(px2, py), 3.5, col);
                }
            }
        }

        ui.separator();
        ui.strong("Burnham & Anderson (2004) Multimodel AICc Inhibition Ranking (Click row to overlay curve):");
        egui::Grid::new("inh_aicc_grid").striped(true).show(ui, |ui| {
            ui.strong("Rank & Mechanism");
            ui.strong("Vmax");
            ui.strong("Km (µM)");
            ui.strong("Ki (nM)");
            ui.strong("α (Ki'/Ki)");
            ui.strong("R²");
            ui.strong("AICc");
            ui.strong("ΔAICc");
            ui.strong("Akaike Weight");
            ui.end_row();

            for (idx, f) in self.inh_fits.iter().enumerate() {
                if ui
                    .selectable_label(
                        self.selected_mech_idx == idx,
                        format!("#{} {}", idx + 1, f.mechanism.label()),
                    )
                    .clicked()
                {
                    self.selected_mech_idx = idx;
                }
                ui.label(format!("{:.2}", f.vmax));
                ui.label(format!("{:.2}", f.km_um));
                ui.label(format!("{:.2}", f.ki_nm));
                ui.label(
                    f.alpha_factor
                        .map(|a| format!("{a:.2}"))
                        .unwrap_or_else(|| String::from("—")),
                );
                ui.label(format!("{:.5}", f.r_squared));
                ui.label(format!("{:.2}", f.aicc));
                ui.label(format!("{:+.2}", f.delta_aicc));
                ui.strong(format!("{:.1}%", f.akaike_weight * 100.0));
                ui.end_row();
            }
        });
    }

    fn draw_morrison_view(&mut self, ui: &mut egui::Ui) {
        ui.columns(2, |cols| {
            cols[0].heading("Morrison (1969) Tight-Binding Controls");
            let mut changed = false;
            changed |= cols[0]
                .add(
                    egui::Slider::new(&mut self.nominal_et_nm, 1.0..=100.0)
                        .text("Nominal [E] (nM)"),
                )
                .changed();
            changed |= cols[0]
                .add(egui::Slider::new(&mut self.morrison_s_um, 1.0..=500.0).text("Assay [S] (µM)"))
                .changed();
            changed |= cols[0]
                .add(egui::Slider::new(&mut self.morrison_km_um, 1.0..=250.0).text("Enzyme Km (µM)"))
                .changed();
            changed |= cols[0]
                .checkbox(&mut self.fix_morrison_et, "Lock [E]_T to Nominal Concentration")
                .changed();
            if changed {
                self.refit_morrison();
            }

            cols[0].separator();
            let mf = &self.morrison_fit;
            cols[0].label(format!("Uninhibited Velocity v0     : {:.2}", mf.v0_unhibited));
            cols[0].strong(format!(
                "Titrated Active [E]_T       : {:.2} nM ({:.1}% active)",
                mf.active_enzyme_et_nm,
                mf.fractional_active_enzyme.unwrap_or(1.0) * 100.0
            ));
            cols[0].strong(format!("Apparent Ki_app             : {:.3} nM", mf.ki_app_nm));
            if let Some(ki_int) = mf.intrinsic_ki_nm {
                cols[0].strong(format!("Intrinsic Competitive Ki    : {:.3} nM", ki_int));
                let naive_ic50 = 0.5 * mf.active_enzyme_et_nm + mf.ki_app_nm;
                let naive_cp = naive_ic50 / (1.0 + self.morrison_s_um / self.morrison_km_um);
                cols[0].colored_label(
                    Color32::from_rgb(248, 113, 113),
                    format!(
                        "Naive IC50 ({naive_ic50:.2} nM) -> Cheng–Prusoff Ki = {naive_cp:.2} nM ({:.1}x error without Morrison!)",
                        naive_cp / ki_int.max(1e-9)
                    ),
                );
            }
            cols[0].label(format!("Nonlinear R²                : {:.5}", mf.r_squared));

            let avail = cols[1].available_size();
            let (resp, painter) =
                cols[1].allocate_painter(Vec2::new(avail.x, avail.y.max(280.0)), Sense::hover());
            let r = resp.rect.shrink(28.0);
            painter.rect_filled(resp.rect, 4.0, Color32::from_rgb(15, 23, 42));
            painter.rect_stroke(r, 2.0, Stroke::new(1.0, Color32::GRAY));

            let i_max = self
                .morrison_obs
                .iter()
                .map(|p| p.0)
                .fold(60.0_f64, f64::max)
                * 1.05;
            let v_max = (mf.v0_unhibited * 1.08).max(1.0);
            let mut line = Vec::with_capacity(100);
            for k in 0..=100 {
                let it = (k as f64 / 100.0) * i_max;
                let vp = eval_morrison_velocity(
                    it,
                    mf.v0_unhibited,
                    mf.active_enzyme_et_nm,
                    mf.ki_app_nm,
                );
                let px = r.left() + ((it / i_max) as f32) * r.width();
                let py = r.bottom() - ((vp / v_max).clamp(0.0, 1.0) as f32) * r.height();
                line.push(Pos2::new(px, py));
            }
            painter.add(Shape::line(
                line,
                Stroke::new(2.5, Color32::from_rgb(248, 113, 113)),
            ));
            for &(it, vy) in &self.morrison_obs {
                let px = r.left() + ((it / i_max) as f32) * r.width();
                let py = r.bottom() - ((vy / v_max).clamp(0.0, 1.0) as f32) * r.height();
                painter.circle_filled(Pos2::new(px, py), 4.0, Color32::from_rgb(56, 189, 248));
            }
        });
    }

    fn draw_single_and_cp_view(&mut self, ui: &mut egui::Ui) {
        ui.heading("Single-Substrate Rate Laws ([I] = 0 Curve Comparison)");
        if ui
            .add(
                egui::Slider::new(&mut self.enzyme_et_for_kcat_nm, 0.5..=200.0)
                    .text("[E]_T (nM) for kcat"),
            )
            .changed()
        {
            if let Ok(sf) =
                fit_single_substrate_models(&self.matrix_points, Some(self.enzyme_et_for_kcat_nm))
            {
                self.single_fits = sf;
            }
        }

        egui::Grid::new("single_sub_grid").striped(true).show(ui, |ui| {
            ui.strong("Model");
            ui.strong("Vmax");
            ui.strong("Km / K0.5 (µM)");
            ui.strong("n_H / K_si");
            ui.strong("kcat (s⁻¹)");
            ui.strong("kcat/Km (M⁻¹ s⁻¹)");
            ui.strong("AICc Weight");
            ui.end_row();

            for sf in &self.single_fits {
                let _ = eval_single_substrate_velocity(
                    sf.model,
                    50.0,
                    sf.vmax,
                    sf.km_or_khalf_um,
                    sf.secondary_param.unwrap_or(1.0),
                );
                ui.label(sf.model.label());
                ui.label(format!("{:.2}", sf.vmax));
                ui.label(format!("{:.2}", sf.km_or_khalf_um));
                ui.label(
                    sf.secondary_param
                        .map(|v| format!("{v:.2}"))
                        .unwrap_or_else(|| String::from("—")),
                );
                ui.label(
                    sf.kcat_per_s
                        .map(|v| format!("{v:.2}"))
                        .unwrap_or_else(|| String::from("—")),
                );
                ui.label(
                    sf.catalytic_efficiency_m_inv_s_inv
                        .map(|v| format!("{v:.2e}"))
                        .unwrap_or_else(|| String::from("—")),
                );
                ui.strong(format!("{:.1}%", sf.akaike_weight * 100.0));
                ui.end_row();
            }
        });

        ui.separator();
        ui.heading("Cheng & Prusoff (1973) Classical IC50 → Ki Calculator");
        ui.horizontal(|ui| {
            for m in [
                InhibitionMechanism::Competitive,
                InhibitionMechanism::Uncompetitive,
                InhibitionMechanism::PureNoncompetitive,
            ] {
                ui.selectable_value(&mut self.cp_mech, m, m.label());
            }
        });
        ui.add(egui::Slider::new(&mut self.cp_ic50_nm, 0.5..=5_000.0).text("Measured IC50 (nM)"));
        ui.add(egui::Slider::new(&mut self.cp_s_um, 1.0..=1_000.0).text("Assay [S] (µM)"));
        ui.add(egui::Slider::new(&mut self.cp_km_um, 1.0..=500.0).text("Substrate Km (µM)"));

        let conv = convert_ic50_cheng_prusoff(
            self.cp_mech,
            self.cp_ic50_nm,
            self.cp_s_um,
            self.cp_km_um,
        );
        ui.strong(format!(
            "Computed Intrinsic Ki = {:.3} nM (under {})",
            conv.computed_ki_nm,
            conv.mechanism.label()
        ));
    }
}
