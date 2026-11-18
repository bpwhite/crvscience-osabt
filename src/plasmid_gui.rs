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

//! Pure-Rust cross-platform desktop (`eframe` / `egui`) GUI for Tool #2: `osabt plasmid`.
//! Links directly in memory to `thermo`, `structure`, `gibson`, `plasmid_io`, `plasmid_digest`,
//! and `plasmid_svg` with zero IPC serialization overhead.

use std::f32::consts::PI;
use std::path::PathBuf;

use eframe::egui;
use egui::{Color32, Pos2, Rect, Sense, Shape, Stroke, Vec2};

use crate::gibson::{
    design_gibson_assembly, scan_circular_bindings, CircularBindingSite, GibsonConfig,
    GibsonFragmentInput,
};
use crate::gui_io::{
    draw_status_banner, pick_open_file, pick_save_file, read_bytes_reported, write_text_reported,
};
use crate::plasmid_digest::{
    design_goldengate_assembly, find_unique_cutters, simulate_restriction_digest,
    GoldenGateEnzyme, ENZYME_CATALOG,
};
use crate::plasmid_io::{
    align_trace_to_plasmid, find_orfs, load_plasmid_bytes, parse_ab1_trace, translate_dna,
    write_genbank, PlasmidRecord, SangerTrace, TraceAlignmentResult,
};
use crate::plasmid_svg::{
    gel_migration_fraction, layout_circular_callouts, render_plasmid_map_svg,
    render_virtual_gel_svg, CircularCalloutInput, GelLadder, GelLane,
};
use crate::primer_design::parse_fasta_or_raw;
use crate::structure::g4_hunter_max_window;
use crate::thermo::{calc_perfect_duplex, PolymerasePreset, ReactionConditions};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InspectorTab {
    MapAndPrimers,
    DigestAndGel,
    CloningWizard,
    TraceVerify,
}

/// Interactive `egui` state for `OSABT Plasmid Studio`.
pub struct PlasmidStudioApp {
    record: PlasmidRecord,
    load_path_input: String,
    primer_name_input: String,
    primer_seq_input: String,
    stock_oligos: Vec<(String, String)>,
    bound_primers: Vec<CircularBindingSite>,
    selected_span: (usize, usize),
    drag_anchor_bp: Option<usize>,
    selected_feature_idx: Option<usize>,
    active_tab: InspectorTab,
    selected_enzymes: Vec<bool>,
    agarose_pct: f64,
    gg_enzyme: GoldenGateEnzyme,
    insert_name_input: String,
    insert_seq_input: String,
    trace_path_input: String,
    last_trace_result: Option<TraceAlignmentResult>,
    status_bar: String,
    status_is_error: bool,
    about_open: bool,
}

impl PlasmidStudioApp {
    pub fn new(initial_file: Option<PathBuf>) -> Self {
        let mut init_err: Option<String> = None;
        let record = if let Some(p) = &initial_file {
            match std::fs::read(p) {
                Ok(b) => match load_plasmid_bytes(&p.to_string_lossy(), &b) {
                    Ok(rec) => rec,
                    Err(e) => {
                        init_err = Some(format!("Plasmid parse error '{}': {e} (showing demo)", p.display()));
                        PlasmidRecord::demo_pcrv_gfp()
                    }
                },
                Err(e) => {
                    init_err = Some(format!("Read error '{}': {e} (showing demo)", p.display()));
                    PlasmidRecord::demo_pcrv_gfp()
                }
            }
        } else {
            PlasmidRecord::demo_pcrv_gfp()
        };

        let mut selected_enzymes = vec![false; ENZYME_CATALOG.len()];
        for (i, e) in ENZYME_CATALOG.iter().enumerate() {
            if e.name == "EcoRI" || e.name == "XhoI" {
                selected_enzymes[i] = true;
            }
        }

        let default_primer = record.sequence[..22.min(record.sequence.len())].to_string();
        let status_is_error = init_err.is_some();
        let status_bar = init_err.unwrap_or_else(|| {
            String::from(
                "Loaded construct. Drag on the ring to select any bp span, or click a feature arc/label to inspect.",
            )
        });
        let mut app = Self {
            record,
            load_path_input: initial_file
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_default(),
            primer_name_input: String::from("J23100_Fwd"),
            primer_seq_input: default_primer.clone(),
            stock_oligos: vec![(String::from("J23100_Fwd"), default_primer)],
            bound_primers: Vec::new(),
            selected_span: (0, 21),
            drag_anchor_bp: None,
            selected_feature_idx: None,
            active_tab: InspectorTab::MapAndPrimers,
            selected_enzymes,
            agarose_pct: 1.0,
            gg_enzyme: GoldenGateEnzyme::BsaI,
            insert_name_input: String::from("mCherry_Insert"),
            insert_seq_input: String::from(
                "ATGGTGAGCAAGGGCGAGGAGGATAACATGGCCATCATCAAGGAGTTCATGCGCTTCAAGGTGCACATGGAGGGCTCCGTGAACGGCCACGAGTTCGAGATCGAGGGCGAGGGCGAGGGCCGCCCCTACGAGGGCACCCAGACCGCCAAGCTGAAGGTGACCAAGGGTGGCCCCCTGCCCTTCGCCTGGGACATCCTGTCCCCTCAGTTCATGTACGGCTCCAAGGCCTACGTGAAGCACCCCGCCGACATCCCCGACTACTTGAAGCTGTCCTTCCCCGAGGGCTTCAAGTGGGAGCGCGTGATGAACTTCGAGGACGGCGGCGTGGTGACCGTGACCCAGGACTCCTCCCTGCAGGACGGCGAGTTCATCTACAAGGTGAAGCTGCGCGGCACCAACTTCCCCTCCGACGGCCCCGTAATGCAGAAGAAGACCATGGGCTGGGAGGCCTCCTCCGAGCGGATGTACCCCGAGGACGGCGCCCTGAAGGGCGAGATCAAGCAGAGGCTGAAGCTGAAGGACGGCGGCCACTACGACGCTGAGGTCAAGACCACCTACAAGGCCAAGAAGCCCGTGCAGCTGCCCGGCGCCTACAACGTCAACATCAAGTTGGACATCACCTCCCACAACGAGGACTACACCATCGTGGAACAGTACGAACGCGCCGAGGGCCGCCACTCCACCGGCGGCATGGACGAGCTGTACAAGTAA",
            ),
            trace_path_input: String::new(),
            last_trace_result: None,
            status_bar,
            status_is_error,
            about_open: false,
        };
        app.rebind_primers();
        app
    }

    fn load_plasmid_path(&mut self, path: &std::path::Path) {
        if let Some(bytes) = read_bytes_reported(path, &mut self.status_bar, &mut self.status_is_error) {
            match load_plasmid_bytes(&path.to_string_lossy(), &bytes) {
                Ok(rec) => {
                    self.selected_span = (0, 21.min(rec.sequence.len().saturating_sub(1)));
                    self.selected_feature_idx = None;
                    self.drag_anchor_bp = None;
                    self.record = rec;
                    self.rebind_primers();
                    self.status_is_error = false;
                    self.status_bar = format!("Loaded '{}' ({} bp).", self.record.name, self.record.sequence.len());
                }
                Err(e) => {
                    self.status_is_error = true;
                    self.status_bar = format!("Plasmid parse error: {e}");
                }
            }
        }
    }

    fn rebind_primers(&mut self) {
        let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);
        if let Ok(hits) = scan_circular_bindings(
            &self.record.sequence,
            self.record.circular,
            &self.stock_oligos,
            14,
            1,
            &cond,
        ) {
            self.bound_primers = hits;
        }
    }

    fn active_enzyme_names(&self) -> Vec<String> {
        ENZYME_CATALOG
            .iter()
            .zip(&self.selected_enzymes)
            .filter_map(|(e, &on)| if on { Some(e.name.to_string()) } else { None })
            .collect()
    }
}

/// Launches the native cross-platform `OSABT Plasmid Studio` window.
pub fn run_plasmid_gui(initial_file: Option<PathBuf>) -> Result<(), String> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1240.0, 800.0])
            .with_icon(crate::brand::osabt_icon_data())
            .with_title("OSABT Plasmid Studio • CRV Science Network (GPL-3.0)"),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-plasmid-studio",
        options,
        Box::new(move |_cc| Ok(Box::new(PlasmidStudioApp::new(initial_file)))),
    )
    .map_err(|e| format!("Failed to start egui desktop window: {e}"))
}

impl eframe::App for PlasmidStudioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        crate::brand::draw_standalone_footer(ctx, &mut self.about_open);
        self.draw_panels(ctx);
    }
}

impl PlasmidStudioApp {
    /// Renders all Plasmid Studio panels onto `ctx`.
    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("plasmid_top_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading(format!(
                    "{} ({} bp, {})",
                    self.record.name,
                    self.record.sequence.len(),
                    if self.record.circular {
                        "Circular"
                    } else {
                        "Linear"
                    }
                ));
                ui.separator();
                ui.label("File (.gbk / .dna / .fa):");
                ui.add(egui::TextEdit::singleline(&mut self.load_path_input).desired_width(180.0));
                if ui.button("Browse...").clicked() {
                    if let Some(p) = pick_open_file(
                        "Open Plasmid File",
                        "Plasmid / Sequence",
                        &["gbk", "gb", "dna", "fa", "fasta", "fna"],
                    ) {
                        self.load_path_input = p.to_string_lossy().to_string();
                        self.load_plasmid_path(&p);
                    }
                }
                if ui.button("Load").clicked() && !self.load_path_input.trim().is_empty() {
                    let p = PathBuf::from(self.load_path_input.trim());
                    self.load_plasmid_path(&p);
                }
                if ui.button("Demo GFP Vector").clicked() {
                    self.record = PlasmidRecord::demo_pcrv_gfp();
                    self.selected_span = (0, 21);
                    self.selected_feature_idx = None;
                    self.drag_anchor_bp = None;
                    self.rebind_primers();
                    self.status_is_error = false;
                    self.status_bar = String::from("Reset to built-in pCRV_DemoGFP vector.");
                }
                if ui.button("Export .gbk").clicked() {
                    let default_name = format!("{}.gbk", self.record.name);
                    if let Some(path) = pick_save_file("Export GenBank File", &default_name, "GenBank", &["gbk", "gb"]) {
                        write_text_reported(
                            &path,
                            &write_genbank(&self.record),
                            &mut self.status_bar,
                            &mut self.status_is_error,
                            "Exported GenBank",
                        );
                    }
                }
                if ui.button("Export Map .svg").clicked() {
                    let default_name = format!("{}_map.svg", self.record.name);
                    if let Some(path) = pick_save_file("Export Plasmid Map SVG", &default_name, "SVG Image", &["svg"]) {
                        let svg = render_plasmid_map_svg(&self.record, &self.bound_primers, true);
                        write_text_reported(
                            &path,
                            &svg,
                            &mut self.status_bar,
                            &mut self.status_is_error,
                            "Exported Vector Map SVG",
                        );
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("plasmid_bottom_track")
            .min_height(138.0)
            .show(ctx, |ui| {
                self.draw_linear_sequence_hud(ui);
                ui.separator();
                draw_status_banner(ui, &self.status_bar, self.status_is_error);
            });

        egui::SidePanel::right("plasmid_right_inspector")
            .default_width(450.0)
            .resizable(true)
            .show(ctx, |ui| {
                ui.horizontal(|ui| {
                    ui.selectable_value(
                        &mut self.active_tab,
                        InspectorTab::MapAndPrimers,
                        "Map & Primers",
                    );
                    ui.selectable_value(
                        &mut self.active_tab,
                        InspectorTab::DigestAndGel,
                        "Digest & Gel",
                    );
                    ui.selectable_value(
                        &mut self.active_tab,
                        InspectorTab::CloningWizard,
                        "Gibson & Golden Gate",
                    );
                    ui.selectable_value(
                        &mut self.active_tab,
                        InspectorTab::TraceVerify,
                        "Trace (.ab1)",
                    );
                });
                ui.separator();

                egui::ScrollArea::vertical().show(ui, |ui| match self.active_tab {
                    InspectorTab::MapAndPrimers => self.draw_map_and_primers_tab(ui),
                    InspectorTab::DigestAndGel => self.draw_digest_and_gel_tab(ui),
                    InspectorTab::CloningWizard => self.draw_cloning_wizard_tab(ui),
                    InspectorTab::TraceVerify => self.draw_trace_verify_tab(ui),
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            self.draw_circular_plasmid_canvas(ui);
        });
    }
}

impl PlasmidStudioApp {
    fn draw_circular_plasmid_canvas(&mut self, ui: &mut egui::Ui) {
        let avail = ui.available_size();
        let (response, painter) = ui.allocate_painter(avail, Sense::click_and_drag());
        let rect = response.rect;
        let center = rect.center();
        let radius = (rect.width() * 0.21).min(rect.height() * 0.29).max(76.0);
        let n = self.record.sequence.len().max(1);
        let r_feat_band = radius + 16.0;
        let r_feat_outer = radius + 24.0;

        // Precompute unified hemisphere sweep-line callouts for features + unique cutters
        let mut callout_inputs: Vec<CircularCalloutInput> = Vec::new();
        for (idx, feat) in self.record.features.iter().enumerate() {
            let start_a = (feat.start_0 as f64 / n as f64) * 2.0 * std::f64::consts::PI;
            let sweep_a = (feat.span_bp(n) as f64 / n as f64) * 2.0 * std::f64::consts::PI;
            callout_inputs.push(CircularCalloutInput {
                label: feat.name.clone(),
                angle_rad: start_a + sweep_a * 0.5,
                anchor_radius: (r_feat_band + 4.0) as f64,
                is_feature: true,
                color_hex: feat.color_hex.clone(),
                source_index: idx,
            });
        }
        let uniques = find_unique_cutters(&self.record);
        for (u_idx, site) in uniques.iter().take(14).enumerate() {
            let ang = (site.cut_top_0 as f64 / n as f64) * 2.0 * std::f64::consts::PI;
            callout_inputs.push(CircularCalloutInput {
                label: format!("{} ({})", site.enzyme_name, site.cut_top_0 + 1),
                angle_rad: ang,
                anchor_radius: radius as f64,
                is_feature: false,
                color_hex: String::from("#38bdf8"),
                source_index: u_idx,
            });
        }

        let placed_callouts = layout_circular_callouts(
            center.x as f64,
            center.y as f64,
            r_feat_outer as f64,
            (
                (rect.left() + 115.0).min(center.x - 40.0) as f64,
                (rect.top() + 14.0) as f64,
                (rect.right() - 115.0).max(center.x + 40.0) as f64,
                (rect.bottom() - 14.0) as f64,
            ),
            15.5,
            callout_inputs,
        );

        // Pointer drag selection + single-click feature arc/label hit-testing
        let pos_to_bp = |pos: Pos2| -> usize {
            let delta = pos - center;
            let mut angle = delta.y.atan2(delta.x) + PI * 0.5;
            if angle < 0.0 {
                angle += 2.0 * PI;
            }
            (((angle / (2.0 * PI)) * (n as f32)) as usize).min(n - 1)
        };

        if response.drag_started() {
            if let Some(pos) = response.interact_pointer_pos() {
                let d = (pos - center).length();
                if d > radius * 0.45 && d < radius * 1.75 {
                    let bp = pos_to_bp(pos);
                    self.drag_anchor_bp = Some(bp);
                    self.selected_feature_idx = None;
                    self.selected_span = (bp, bp);
                }
            }
        }
        if response.dragged() {
            if let (Some(anchor), Some(pos)) = (self.drag_anchor_bp, response.interact_pointer_pos()) {
                let cur = pos_to_bp(pos);
                let (s0, e0) = if self.record.circular {
                    let fwd_len = if cur >= anchor { cur - anchor } else { n - anchor + cur };
                    let rev_len = n - fwd_len;
                    if fwd_len <= rev_len {
                        (anchor, cur)
                    } else {
                        (cur, anchor)
                    }
                } else {
                    (anchor.min(cur), anchor.max(cur))
                };
                self.selected_span = (s0, e0);
            }
        }
        if response.drag_stopped() {
            if let Some(anchor) = self.drag_anchor_bp.take() {
                if let Some(pos) = response.interact_pointer_pos() {
                    let cur = pos_to_bp(pos);
                    if cur != anchor {
                        let (s0, e0) = self.selected_span;
                        let span_len = if e0 >= s0 { e0 - s0 + 1 } else { (n - s0) + e0 + 1 };
                        self.status_is_error = false;
                        self.status_bar = format!("Selected {} bp ({}..={}).", span_len, s0 + 1, e0 + 1);
                    }
                }
            }
        } else if response.clicked() {
            if let Some(pos) = response.interact_pointer_pos() {
                let mut hit_feat: Option<usize> = None;
                for c in &placed_callouts {
                    if !c.input.is_feature {
                        continue;
                    }
                    let lx = c.p2.0 as f32;
                    let ly = c.p2.1 as f32;
                    let text_size = ui.fonts(|f| {
                        f.layout_no_wrap(
                            c.input.label.clone(),
                            egui::FontId::proportional(11.5),
                            Color32::WHITE,
                        )
                        .size()
                    });
                    let hit_w = text_size.x.max(48.0);
                    let half_h = (text_size.y * 0.5 + 2.0).max(8.5);
                    let lbl_rect = if c.is_right_hemisphere {
                        Rect::from_min_max(
                            Pos2::new(lx + 2.0, ly - half_h),
                            Pos2::new(lx + 6.0 + hit_w, ly + half_h),
                        )
                    } else {
                        Rect::from_min_max(
                            Pos2::new(lx - 6.0 - hit_w, ly - half_h),
                            Pos2::new(lx - 2.0, ly + half_h),
                        )
                    };
                    if lbl_rect.contains(pos) {
                        hit_feat = Some(c.input.source_index);
                        break;
                    }
                }
                if hit_feat.is_none() {
                    let delta = pos - center;
                    let d = delta.length();
                    if (d - r_feat_band).abs() <= 9.0 {
                        let bp = pos_to_bp(pos);
                        for (idx, feat) in self.record.features.iter().enumerate() {
                            let inside = if feat.start_0 <= feat.end_0 {
                                bp >= feat.start_0 && bp <= feat.end_0
                            } else {
                                bp >= feat.start_0 || bp <= feat.end_0
                            };
                            if inside {
                                hit_feat = Some(idx);
                                break;
                            }
                        }
                    }
                }
                if let Some(fidx) = hit_feat {
                    if let Some(feat) = self.record.features.get(fidx) {
                        self.selected_feature_idx = Some(fidx);
                        self.selected_span = (feat.start_0.min(n - 1), feat.end_0.min(n - 1));
                        self.status_is_error = false;
                        self.status_bar = format!(
                            "Feature '{}' [{}] {}..={} ({} bp, {:?})",
                            feat.name,
                            feat.feature_type,
                            feat.start_0 + 1,
                            feat.end_0 + 1,
                            feat.span_bp(n),
                            feat.strand
                        );
                    }
                } else {
                    let delta = pos - center;
                    if delta.length() > radius * 0.45 && delta.length() < radius * 1.45 {
                        let bp = pos_to_bp(pos);
                        self.selected_feature_idx = None;
                        self.selected_span = (bp, bp);
                        self.status_is_error = false;
                        self.status_bar = format!("Cursor at bp {} (drag ring to select span).", bp + 1);
                    }
                }
            }
        }

        // Backbone double-ring
        painter.circle_stroke(center, radius, Stroke::new(3.0_f32, Color32::from_rgb(203, 213, 225)));
        painter.circle_stroke(
            center,
            radius - 5.0,
            Stroke::new(1.2_f32, Color32::from_rgb(100, 116, 139)),
        );

        // Highlight currently selected window arc (supporting origin-wrapping circular spans)
        let (s0, e0) = self.selected_span;
        let sel_bp = if e0 >= s0 { e0 - s0 + 1 } else { (n - s0) + e0 + 1 };
        let pts = arc_polyline(
            center,
            radius - 12.0,
            (s0 as f32 / n as f32) * 2.0 * PI,
            (sel_bp.max(2) as f32 / n as f32) * 2.0 * PI,
            32,
        );
        painter.add(Shape::line(
            pts,
            Stroke::new(6.0_f32, Color32::from_rgb(250, 204, 21)),
        ));

        // Annotated feature arcs
        for (idx, feat) in self.record.features.iter().enumerate() {
            let start_a = (feat.start_0 as f32 / n as f32) * 2.0 * PI;
            let sweep_a = (feat.span_bp(n) as f32 / n as f32) * 2.0 * PI;
            let color = parse_hex_color(&feat.color_hex);
            let pts = arc_polyline(center, r_feat_band, start_a, sweep_a.max(0.05), 36);
            if self.selected_feature_idx == Some(idx) {
                painter.add(Shape::line(pts.clone(), Stroke::new(12.0_f32, Color32::from_rgb(250, 204, 21))));
            }
            painter.add(Shape::line(pts, Stroke::new(8.0_f32, color)));
        }

        // Draw unified relaxed callouts + 3-point elbow leaders
        for c in &placed_callouts {
            let p0 = Pos2::new(c.p0.0 as f32, c.p0.1 as f32);
            let p1 = Pos2::new(c.p1.0 as f32, c.p1.1 as f32);
            let p2 = Pos2::new(c.p2.0 as f32, c.p2.1 as f32);
            let align = if c.is_right_hemisphere {
                egui::Align2::LEFT_CENTER
            } else {
                egui::Align2::RIGHT_CENTER
            };
            let text_dx = if c.is_right_hemisphere { 4.0 } else { -4.0 };
            if c.input.is_feature {
                let col = parse_hex_color(&c.input.color_hex);
                let is_sel = self.selected_feature_idx == Some(c.input.source_index);
                painter.add(Shape::line(
                    vec![p0, p1, p2],
                    Stroke::new(if is_sel { 1.8 } else { 1.1_f32 }, col),
                ));
                painter.text(
                    Pos2::new(p2.x + text_dx, p2.y),
                    align,
                    &c.input.label,
                    egui::FontId::proportional(11.5),
                    if is_sel { Color32::from_rgb(250, 204, 21) } else { col },
                );
            } else {
                let enz_col = Color32::from_rgb(125, 211, 252);
                painter.add(Shape::line(
                    vec![p0, p1, p2],
                    Stroke::new(1.0_f32, Color32::from_rgb(56, 189, 248)),
                ));
                painter.text(
                    Pos2::new(p2.x + text_dx, p2.y),
                    align,
                    &c.input.label,
                    egui::FontId::monospace(10.5),
                    enz_col,
                );
            }
        }

        // Bound primers
        for p in self.bound_primers.iter().take(12) {
            let start_a = (p.bind_start_0 as f32 / n as f32) * 2.0 * PI;
            let sweep_a = ((p.bound_len as f32 / n as f32) * 2.0 * PI).max(0.06);
            let pts = arc_polyline(center, radius - 24.0, start_a, sweep_a, 16);
            painter.add(Shape::line(
                pts,
                Stroke::new(4.0_f32, Color32::from_rgb(244, 114, 182)),
            ));
        }

        // Center nameplate
        painter.text(
            center - Vec2::new(0.0, 10.0),
            egui::Align2::CENTER_CENTER,
            &self.record.name,
            egui::FontId::proportional(16.0),
            Color32::WHITE,
        );
        painter.text(
            center + Vec2::new(0.0, 12.0),
            egui::Align2::CENTER_CENTER,
            format!("{} bp", self.record.sequence.len()),
            egui::FontId::monospace(12.0),
            Color32::LIGHT_GRAY,
        );
    }

    fn draw_linear_sequence_hud(&mut self, ui: &mut egui::Ui) {
        let n = self.record.sequence.len();
        if n == 0 {
            return;
        }
        let s0 = self.selected_span.0.min(n - 1);
        let e0 = self.selected_span.1.min(n - 1);
        let wrapped_sub: String = if e0 >= s0 {
            self.record.sequence[s0..=e0].to_string()
        } else {
            format!("{}{}", &self.record.sequence[s0..], &self.record.sequence[..=e0])
        };
        let span_len = wrapped_sub.len();
        let sub_slice = &wrapped_sub[..wrapped_sub.len().min(500)];
        let aa = translate_dna(sub_slice);
        let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);
        let tm_txt = calc_perfect_duplex(sub_slice, &cond, 60.0)
            .map(|r| format!("Q5 Tm: {:.2} C | GC: {:.1}% | dG37: {:.2} kcal/mol", r.tm_salt_corrected_c, r.gc_pct, r.delta_g_37_kcal))
            .unwrap_or_else(|_| String::from("Select >= 2 bp for Tm"));
        let g4 = g4_hunter_max_window(sub_slice, 20);

        ui.horizontal(|ui| {
            ui.strong(format!(
                "Selection: {}..={} ({} bp)",
                s0 + 1,
                e0 + 1,
                span_len
            ));
            ui.separator();
            ui.label(tm_txt);
            ui.separator();
            ui.label(format!("G4Hunter: {g4:.2}"));
            if ui.button("Use Selection as Primer").clicked() {
                self.primer_seq_input = wrapped_sub.clone();
                self.primer_name_input = format!("Prm_{}_{}", s0 + 1, e0 + 1);
            }
        });

        let preview: String = wrapped_sub.chars().take(78).collect();
        let aa_preview: String = aa.chars().take(60).collect();
        ui.monospace(format!("5' {:>5} {} 3'", s0 + 1, preview));
        ui.monospace(format!("   Frame+1 AA: {aa_preview}"));
    }

    fn draw_map_and_primers_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Freezer Primer Auto-Binder");
        ui.horizontal(|ui| {
            ui.label("Name:");
            ui.add(egui::TextEdit::singleline(&mut self.primer_name_input).desired_width(110.0));
            ui.label("5'->3':");
            ui.add(egui::TextEdit::singleline(&mut self.primer_seq_input).desired_width(180.0));
            if ui.button("Bind").clicked() && !self.primer_seq_input.trim().is_empty() {
                self.stock_oligos.push((
                    self.primer_name_input.trim().to_string(),
                    self.primer_seq_input.trim().to_string(),
                ));
                self.rebind_primers();
            }
        });

        ui.label(format!("Bound Sites on Vector: {}", self.bound_primers.len()));
        for b in &self.bound_primers {
            ui.monospace(format!(
                "• {} ({:?}) {}..={} | core Tm {:.1} C | full Tm {:.1} C | 5'tail '{}'",
                b.oligo_name,
                b.strand,
                b.bind_start_0 + 1,
                b.bind_end_0 + 1,
                b.core_tm_c,
                b.full_tm_c,
                b.overhang_5p_seq
            ));
        }

        ui.separator();
        ui.heading(format!("Annotated Features ({})", self.record.features.len()));
        for f in &self.record.features {
            ui.monospace(format!(
                "• {:<16} [{}] {}..={} ({:?})",
                f.name,
                f.feature_type,
                f.start_0 + 1,
                f.end_0 + 1,
                f.strand
            ));
        }

        ui.separator();
        let orfs = find_orfs(&self.record, 40);
        ui.heading(format!("6-Frame ORFs >= 40 aa ({})", orfs.len()));
        for o in orfs.iter().take(8) {
            ui.monospace(format!(
                "• Frame {:+2} | {}..={} ({} bp, {} aa)",
                o.frame,
                o.start_0 + 1,
                o.end_0 + 1,
                o.length_bp,
                o.aa_sequence.len()
            ));
        }
    }

    fn draw_digest_and_gel_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Restriction Digest & Virtual Agarose Gel");
        ui.add(egui::Slider::new(&mut self.agarose_pct, 0.6..=2.5).text("Agarose % (w/v)"));

        egui::Grid::new("enz_grid").num_columns(4).show(ui, |ui| {
            for (idx, enz) in ENZYME_CATALOG.iter().enumerate() {
                ui.checkbox(&mut self.selected_enzymes[idx], enz.name);
                if idx % 4 == 3 {
                    ui.end_row();
                }
            }
        });

        let active = self.active_enzyme_names();
        let frags = simulate_restriction_digest(&self.record, &active).unwrap_or_default();
        let lengths: Vec<usize> = frags.iter().map(|f| f.length_bp).collect();

        ui.horizontal(|ui| {
            ui.label(format!("Predicted Bands: {}", frags.len()));
            if ui.button("Export Virtual Gel .svg").clicked() {
                let lane_lbl = if active.is_empty() {
                    String::from("Uncut")
                } else {
                    active.join("+")
                };
                let svg = render_virtual_gel_svg(
                    GelLadder::OneKbPlus,
                    &[GelLane {
                        label: lane_lbl,
                        fragment_lengths_bp: lengths.clone(),
                    }],
                    self.agarose_pct,
                    true,
                    true,
                );
                let default_name = format!("{}_gel.svg", self.record.name);
                if let Some(path) = pick_save_file("Export Virtual Gel SVG", &default_name, "SVG Image", &["svg"]) {
                    write_text_reported(
                        &path,
                        &svg,
                        &mut self.status_bar,
                        &mut self.status_is_error,
                        "Saved Virtual Gel SVG",
                    );
                }
            }
        });

        for f in &frags {
            ui.monospace(format!(
                "  Band #{}: {:>6} bp ({} -> {})",
                f.index, f.length_bp, f.left_enzyme, f.right_enzyme
            ));
        }

        // Paint live darkroom SYBR/EtBr gel preview directly inside the egui panel
        let (gel_resp, gp) = ui.allocate_painter(Vec2::new(260.0, 250.0), Sense::hover());
        let grect = gel_resp.rect;
        gp.rect_filled(grect, 6.0, Color32::from_rgb(15, 23, 42));

        // Lane 1: 1 kb Plus Ladder; Lane 2: Digest Sample
        for &(bp, bright) in GelLadder::OneKbPlus.rungs() {
            let y = grect.top()
                + 16.0
                + (gel_migration_fraction(bp, self.agarose_pct) as f32) * (grect.height() - 28.0);
            let col = if bright {
                Color32::from_rgb(253, 224, 71)
            } else {
                Color32::from_rgba_unmultiplied(253, 224, 71, 140)
            };
            gp.line_segment(
                [Pos2::new(grect.left() + 28.0, y), Pos2::new(grect.left() + 98.0, y)],
                Stroke::new(if bright { 3.0 } else { 1.6_f32 }, col),
            );
        }
        for &bp in &lengths {
            let y = grect.top()
                + 16.0
                + (gel_migration_fraction(bp, self.agarose_pct) as f32) * (grect.height() - 28.0);
            gp.line_segment(
                [Pos2::new(grect.left() + 140.0, y), Pos2::new(grect.left() + 220.0, y)],
                Stroke::new(3.2_f32, Color32::from_rgb(134, 239, 172)),
            );
        }
    }

    fn draw_cloning_wizard_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Seamless Insert Cloning (Gibson & Golden Gate)");
        ui.label("Insert Name:");
        ui.text_edit_singleline(&mut self.insert_name_input);
        ui.label("Insert Sequence (5'->3'):");
        ui.add(egui::TextEdit::multiline(&mut self.insert_seq_input).desired_rows(4));

        ui.separator();
        ui.strong("1. Two-Stage Gibson / HiFi Preview");
        let gib_frags = vec![
            GibsonFragmentInput {
                name: self.record.name.clone(),
                sequence: self.record.sequence.clone(),
            },
            GibsonFragmentInput {
                name: self.insert_name_input.clone(),
                sequence: self.insert_seq_input.clone(),
            },
        ];
        if let Ok(gplan) = design_gibson_assembly(&gib_frags, true, &GibsonConfig::default()) {
            ui.label(format!(
                "Assembled Circular Construct: {} bp ({} junctions)",
                gplan.assembled_length_bp,
                gplan.junctions.len()
            ));
            for j in &gplan.junctions {
                ui.monospace(format!(
                    "• J{}: {} bp overlap (Tm {:.1} C, dG50 {:.1} kcal/mol)",
                    j.junction_index + 1,
                    j.overlap_len,
                    j.overlap_tm_c,
                    j.overlap_dg_50c_kcal
                ));
                ui.monospace(format!("  F: {}", j.downstream_fwd_primer.annotated_sequence));
                ui.monospace(format!("  R: {}", j.upstream_rev_primer.annotated_sequence));
            }
            if ui.button("Export Gibson Primers .tsv").clicked() {
                if let Some(path) = pick_save_file("Export Gibson Primers TSV", "osabt_gibson_primers.tsv", "TSV", &["tsv"]) {
                    write_text_reported(
                        &path,
                        &gplan.to_order_tsv(),
                        &mut self.status_bar,
                        &mut self.status_is_error,
                        "Saved Gibson Primers TSV",
                    );
                }
            }
        }

        ui.separator();
        ui.horizontal(|ui| {
            ui.strong("2. Type IIS Golden Gate:");
            ui.selectable_value(&mut self.gg_enzyme, GoldenGateEnzyme::BsaI, "BsaI");
            ui.selectable_value(&mut self.gg_enzyme, GoldenGateEnzyme::BsmBI, "BsmBI");
            ui.selectable_value(&mut self.gg_enzyme, GoldenGateEnzyme::BbsI, "BbsI");
            ui.selectable_value(&mut self.gg_enzyme, GoldenGateEnzyme::SapI, "SapI");
        });
        let gg_parts = vec![
            (self.record.name.clone(), self.record.sequence.clone()),
            (self.insert_name_input.clone(), self.insert_seq_input.clone()),
        ];
        if let Ok(gg) = design_goldengate_assembly(&gg_parts, self.gg_enzyme, true) {
            ui.label(format!(
                "Predicted 4-bp Overhang Set Fidelity: {:.1}% | Internal Illegal Sites: {}",
                gg.estimated_set_fidelity_pct, gg.internal_illegal_sites
            ));
            for d in &gg.domestication_edits {
                ui.colored_label(
                    Color32::YELLOW,
                    format!(
                        "Domesticate {} at bp {}: {} -> {} (preserves '{}' AA)",
                        d.fragment_name,
                        d.codon_start_0 + 1,
                        d.original_codon,
                        d.silent_codon,
                        d.amino_acid
                    ),
                );
            }
            for j in &gg.junctions {
                ui.monospace(format!(
                    "• Overhang [{}] | F ({:.1} C): {}",
                    j.overhang_seq, j.fwd_core_tm_c, j.downstream_fwd_primer
                ));
            }
        }
    }

    fn draw_trace_verify_tab(&mut self, ui: &mut egui::Ui) {
        ui.heading("Sanger (.ab1) & Whole-Plasmid FASTA Aligner");
        ui.horizontal(|ui| {
            ui.label("Trace File:");
            ui.add(egui::TextEdit::singleline(&mut self.trace_path_input).desired_width(180.0));
            if ui.button("Browse...").clicked() {
                if let Some(p) = pick_open_file("Open Sanger Trace or FASTA", "Trace / FASTA", &["ab1", "fa", "fasta", "txt"]) {
                    self.trace_path_input = p.to_string_lossy().to_string();
                }
            }
            if ui.button("Align Read").clicked() && !self.trace_path_input.trim().is_empty() {
                let p = PathBuf::from(self.trace_path_input.trim());
                if let Some(bytes) = read_bytes_reported(&p, &mut self.status_bar, &mut self.status_is_error) {
                    let trace_res: Result<SangerTrace, String> = if bytes.starts_with(b"ABIF") {
                        parse_ab1_trace(&bytes).map_err(|e| format!("AB1 parse error: {e}"))
                    } else {
                        std::str::from_utf8(&bytes)
                            .map_err(|e| format!("Trace file is not valid ABIF or UTF-8 FASTA: {e}"))
                            .and_then(|t| parse_fasta_or_raw(t).map_err(|e| format!("FASTA parse error: {e}")))
                            .and_then(|mut r| {
                                r.pop()
                                    .ok_or_else(|| String::from("No sequence found in trace file."))
                            })
                            .map(|(sample_name, sequence)| {
                                let qualities = vec![40u8; sequence.len()];
                                SangerTrace {
                                    sample_name,
                                    sequence,
                                    qualities,
                                }
                            })
                    };
                    match trace_res {
                        Ok(tr) => {
                            let name = tr.sample_name.clone();
                            self.last_trace_result = align_trace_to_plasmid(&self.record, &tr, 20);
                            if self.last_trace_result.is_some() {
                                self.status_is_error = false;
                                self.status_bar = format!("Aligned trace '{name}'.");
                            } else {
                                self.status_is_error = true;
                                self.status_bar = format!("Trace '{name}' could not align (>= 20 bp anchor) to plasmid.");
                            }
                        }
                        Err(e) => {
                            self.status_is_error = true;
                            self.status_bar = e;
                        }
                    }
                }
            }
        });

        if let Some(aln) = &self.last_trace_result {
            ui.label(format!(
                "Aligned {} bp ({:?}) across {}..={} | Identity: {:.2}% | Mean Phred Q: {:.1}",
                aln.aligned_len_bp,
                aln.strand,
                aln.plasmid_start_0 + 1,
                aln.plasmid_end_0 + 1,
                aln.identity_pct,
                aln.mean_phred_q
            ));
            if aln.variants.is_empty() {
                ui.colored_label(
                    Color32::GREEN,
                    "[OK] 100% sequence concordance across aligned span (zero substitutions).",
                );
            } else {
                for v in &aln.variants {
                    ui.monospace(format!(
                        "  bp {:>5}: Ref {} -> Read {} (Q{}) [{}]",
                        v.plasmid_pos_0 + 1,
                        v.ref_base,
                        v.read_base,
                        v.read_phred_q,
                        v.classification
                    ));
                }
            }
        }
    }
}

fn arc_polyline(
    center: Pos2,
    radius: f32,
    start_angle: f32,
    sweep: f32,
    steps: usize,
) -> Vec<Pos2> {
    let s = steps.max(2);
    (0..=s)
        .map(|i| {
            let t = i as f32 / s as f32;
            let a = start_angle + sweep * t - PI * 0.5;
            center + Vec2::new(a.cos(), a.sin()) * radius
        })
        .collect()
}

fn parse_hex_color(hex: &str) -> Color32 {
    let h = hex.trim_start_matches('#');
    if h.len() == 6 {
        if let (Ok(r), Ok(g), Ok(b)) = (
            u8::from_str_radix(&h[0..2], 16),
            u8::from_str_radix(&h[2..4], 16),
            u8::from_str_radix(&h[4..6], 16),
        ) {
            return Color32::from_rgb(r, g, b);
        }
    }
    Color32::from_rgb(34, 197, 94)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::gibson::Strand;
    use crate::plasmid_io::PlasmidFeature;

    #[test]
    fn cycles_all_inspector_tabs_and_clicks_long_feature_label_beyond_180px() {
        let ctx = egui::Context::default();
        let mut app = PlasmidStudioApp::new(None);
        for tab in [
            InspectorTab::MapAndPrimers,
            InspectorTab::DigestAndGel,
            InspectorTab::CloningWizard,
            InspectorTab::TraceVerify,
        ] {
            app.active_tab = tab;
            for frame in 0..2 {
                let mut raw = egui::RawInput::default();
                raw.screen_rect = Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1280.0, 800.0)));
                raw.time = Some(frame as f64 / 60.0);
                let out = ctx.run(raw, |c| app.draw_panels(c));
                assert!(!out.shapes.is_empty());
            }
        }

        // Synthesize a plasmid with a 37-char feature at 3 o'clock (right hemisphere)
        let seq = "ATGACCATGATTACGCCAAGCTATTTAGGTGACACTATAGAATACTCAAGCTATGCATCCAACGCGTTGGGAGCTCTCCC".repeat(12);
        let n = seq.len();
        let feat = PlasmidFeature {
            name: String::from("CMV_Immediate_Early_Enhancer_Promoter"),
            feature_type: String::from("promoter"),
            start_0: n / 4 - 20,
            end_0: n / 4 + 20,
            strand: Strand::Forward,
            color_hex: String::from("#22c55e"),
            qualifiers: Vec::new(),
        };
        let mut click_app = PlasmidStudioApp::new(None);
        click_app.record = PlasmidRecord {
            name: String::from("LongLabelTest"),
            definition: String::from("Long label click test"),
            accession: String::from("TEST1"),
            circular: true,
            sequence: seq,
            features: vec![feat],
        };
        click_app.selected_feature_idx = None;

        // Frame 0: measure canvas rect and confirm exact galley placement
        let mut canvas_rect = Rect::NOTHING;
        let mut raw0 = egui::RawInput::default();
        raw0.screen_rect = Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1400.0, 860.0)));
        ctx.set_pixels_per_point(1.0);
        let _ = ctx.run(raw0, |c| {
            egui::CentralPanel::default().show(c, |ui| {
                canvas_rect = ui.available_rect_before_wrap();
                click_app.draw_circular_plasmid_canvas(ui);
            });
        });

        let center = canvas_rect.center();
        let radius = (canvas_rect.width() * 0.21)
            .min(canvas_rect.height() * 0.29)
            .max(76.0);
        let r_feat_outer = radius + 24.0;
        let placed = layout_circular_callouts(
            center.x as f64,
            center.y as f64,
            r_feat_outer as f64,
            (
                (canvas_rect.left() + 115.0).min(center.x - 40.0) as f64,
                (canvas_rect.top() + 14.0) as f64,
                (canvas_rect.right() - 115.0).max(center.x + 40.0) as f64,
                (canvas_rect.bottom() - 14.0) as f64,
            ),
            15.5,
            vec![CircularCalloutInput {
                label: String::from("CMV_Immediate_Early_Enhancer_Promoter"),
                angle_rad: std::f64::consts::FRAC_PI_2,
                anchor_radius: (radius + 20.0) as f64,
                is_feature: true,
                color_hex: String::from("#22c55e"),
                source_index: 0,
            }],
        );
        let p2 = Pos2::new(placed[0].p2.0 as f32, placed[0].p2.1 as f32);
        let measured_w = ctx.fonts(|f| {
            f.layout_no_wrap(
                String::from("CMV_Immediate_Early_Enhancer_Promoter"),
                egui::FontId::proportional(11.5),
                Color32::WHITE,
            )
            .size()
            .x
        });
        assert!(measured_w > 200.0, "expected measured_w > 200 px, got {measured_w}");
        let click_pos = Pos2::new(p2.x + 195.0, p2.y);

        for (step, pressed) in [(1usize, true), (2usize, false)] {
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(Rect::from_min_size(Pos2::ZERO, Vec2::new(1400.0, 860.0)));
            raw.time = Some(step as f64 / 60.0);
            raw.events.push(egui::Event::PointerMoved(click_pos));
            raw.events.push(egui::Event::PointerButton {
                pos: click_pos,
                button: egui::PointerButton::Primary,
                pressed,
                modifiers: egui::Modifiers::NONE,
            });
            let _ = ctx.run(raw, |c| {
                egui::CentralPanel::default().show(c, |ui| {
                    click_app.draw_circular_plasmid_canvas(ui);
                });
            });
        }
        assert_eq!(click_app.selected_feature_idx, Some(0));
    }
}
