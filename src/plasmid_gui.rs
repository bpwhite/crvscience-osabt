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
use std::fs;
use std::path::PathBuf;

use eframe::egui;
use egui::{Color32, Pos2, Sense, Shape, Stroke, Vec2};

use crate::gibson::{
    design_gibson_assembly, scan_circular_bindings, CircularBindingSite, GibsonConfig,
    GibsonFragmentInput,
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
    gel_migration_fraction, render_plasmid_map_svg, render_virtual_gel_svg, GelLadder, GelLane,
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
    active_tab: InspectorTab,
    selected_enzymes: Vec<bool>,
    agarose_pct: f64,
    gg_enzyme: GoldenGateEnzyme,
    insert_name_input: String,
    insert_seq_input: String,
    trace_path_input: String,
    last_trace_result: Option<TraceAlignmentResult>,
    status_bar: String,
}

impl PlasmidStudioApp {
    pub fn new(initial_file: Option<PathBuf>) -> Self {
        let record = initial_file
            .as_ref()
            .and_then(|p| fs::read(p).ok().map(|b| (p, b)))
            .and_then(|(p, b)| load_plasmid_bytes(&p.to_string_lossy(), &b).ok())
            .unwrap_or_else(PlasmidRecord::demo_pcrv_gfp);

        let mut selected_enzymes = vec![false; ENZYME_CATALOG.len()];
        for (i, e) in ENZYME_CATALOG.iter().enumerate() {
            if e.name == "EcoRI" || e.name == "XhoI" {
                selected_enzymes[i] = true;
            }
        }

        let default_primer = record.sequence[..22.min(record.sequence.len())].to_string();
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
            status_bar: String::from(
                "Loaded construct. Click anywhere on the circular ring to inspect nucleotides and live SantaLucia/Owczarzy primer Tm.",
            ),
        };
        app.rebind_primers();
        app
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
        egui::TopBottomPanel::top("plasmid_top_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading(format!(
                    "🧬 {} ({} bp, {})",
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
                ui.add(egui::TextEdit::singleline(&mut self.load_path_input).desired_width(220.0));
                if ui.button("Load").clicked() {
                    match fs::read(&self.load_path_input).and_then(|bytes| {
                        load_plasmid_bytes(&self.load_path_input, &bytes)
                            .map_err(|e| std::io::Error::new(std::io::ErrorKind::InvalidData, e.to_string()))
                    }) {
                        Ok(rec) => {
                            self.selected_span = (0, 21.min(rec.sequence.len().saturating_sub(1)));
                            self.record = rec;
                            self.rebind_primers();
                            self.status_bar =
                                format!("Loaded '{}' successfully.", self.load_path_input);
                        }
                        Err(e) => {
                            self.status_bar = format!("Load error: {e}");
                        }
                    }
                }
                if ui.button("Demo GFP Vector").clicked() {
                    self.record = PlasmidRecord::demo_pcrv_gfp();
                    self.rebind_primers();
                    self.status_bar = String::from("Reset to built-in pCRV_DemoGFP vector.");
                }
                if ui.button("Export .gbk").clicked() {
                    let out_name = format!("{}.gbk", self.record.name);
                    match fs::write(&out_name, write_genbank(&self.record)) {
                        Ok(_) => self.status_bar = format!("Exported GenBank -> ./{out_name}"),
                        Err(e) => self.status_bar = format!("Export error: {e}"),
                    }
                }
                if ui.button("Export Map .svg").clicked() {
                    let out_name = format!("{}_map.svg", self.record.name);
                    let svg = render_plasmid_map_svg(&self.record, &self.bound_primers, true);
                    match fs::write(&out_name, svg) {
                        Ok(_) => self.status_bar = format!("Exported Vector Map SVG -> ./{out_name}"),
                        Err(e) => self.status_bar = format!("SVG export error: {e}"),
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("plasmid_bottom_track")
            .min_height(138.0)
            .show(ctx, |ui| {
                self.draw_linear_sequence_hud(ui);
                ui.separator();
                ui.label(&self.status_bar);
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
        let radius = (rect.width().min(rect.height()) * 0.34).max(80.0);
        let n = self.record.sequence.len().max(1);

        // Interactive click/drag selection around circular ring
        if let Some(pos) = response.interact_pointer_pos() {
            let delta = pos - center;
            if delta.length() > radius * 0.45 && delta.length() < radius * 1.45 {
                let mut angle = delta.y.atan2(delta.x) + PI * 0.5;
                if angle < 0.0 {
                    angle += 2.0 * PI;
                }
                let clicked_bp = (((angle / (2.0 * PI)) * (n as f32)) as usize).min(n - 1);
                let win_end = (clicked_bp + 21).min(n - 1);
                self.selected_span = (clicked_bp, win_end);
            }
        }

        // Backbone double-ring
        painter.circle_stroke(center, radius, Stroke::new(3.0_f32, Color32::from_rgb(203, 213, 225)));
        painter.circle_stroke(
            center,
            radius - 5.0,
            Stroke::new(1.2_f32, Color32::from_rgb(100, 116, 139)),
        );

        // Highlight currently selected window arc
        let (s0, e0) = self.selected_span;
        if e0 >= s0 {
            let pts = arc_polyline(
                center,
                radius - 12.0,
                (s0 as f32 / n as f32) * 2.0 * PI,
                ((e0 - s0 + 1).max(2) as f32 / n as f32) * 2.0 * PI,
                24,
            );
            painter.add(Shape::line(
                pts,
                Stroke::new(6.0_f32, Color32::from_rgb(250, 204, 21)),
            ));
        }

        // Annotated features
        for (idx, feat) in self.record.features.iter().enumerate() {
            let start_a = (feat.start_0 as f32 / n as f32) * 2.0 * PI;
            let sweep_a = (feat.span_bp(n) as f32 / n as f32) * 2.0 * PI;
            let color = parse_hex_color(&feat.color_hex);
            let r_feat = radius + 14.0 + ((idx % 2) as f32) * 10.0;
            let pts = arc_polyline(center, r_feat, start_a, sweep_a.max(0.05), 36);
            painter.add(Shape::line(pts, Stroke::new(8.0_f32, color)));

            let mid_a = start_a + sweep_a * 0.5 - PI * 0.5;
            let lbl_pos = center + Vec2::new(mid_a.cos(), mid_a.sin()) * (r_feat + 18.0);
            painter.text(
                lbl_pos,
                egui::Align2::CENTER_CENTER,
                &feat.name,
                egui::FontId::proportional(12.0),
                color,
            );
        }

        // Unique restriction cutters
        for (u_idx, site) in find_unique_cutters(&self.record).iter().take(14).enumerate() {
            let ang = (site.cut_top_0 as f32 / n as f32) * 2.0 * PI - PI * 0.5;
            let dir = Vec2::new(ang.cos(), ang.sin());
            let r_out = radius + 42.0 + ((u_idx % 2) as f32) * 14.0;
            painter.line_segment(
                [center + dir * radius, center + dir * r_out],
                Stroke::new(1.0_f32, Color32::from_rgb(56, 189, 248)),
            );
            painter.text(
                center + dir * (r_out + 10.0),
                egui::Align2::CENTER_CENTER,
                format!("{} ({})", site.enzyme_name, site.cut_top_0 + 1),
                egui::FontId::monospace(10.5),
                Color32::from_rgb(125, 211, 252),
            );
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
            egui::FontId::proportional(17.0),
            Color32::WHITE,
        );
        painter.text(
            center + Vec2::new(0.0, 12.0),
            egui::Align2::CENTER_CENTER,
            format!("{} bp", self.record.sequence.len()),
            egui::FontId::monospace(12.5),
            Color32::LIGHT_GRAY,
        );
    }

    fn draw_linear_sequence_hud(&mut self, ui: &mut egui::Ui) {
        let n = self.record.sequence.len();
        if n == 0 {
            return;
        }
        let s0 = self.selected_span.0.min(n - 1);
        let e0 = self.selected_span.1.min(n - 1).max(s0);
        let sub = &self.record.sequence[s0..=e0];
        let aa = translate_dna(sub);
        let cond = ReactionConditions::for_preset(PolymerasePreset::Q5);
        let tm_txt = calc_perfect_duplex(sub, &cond, 60.0)
            .map(|r| format!("Q5 Tm: {:.2} C | GC: {:.1}% | dG37: {:.2} kcal/mol", r.tm_salt_corrected_c, r.gc_pct, r.delta_g_37_kcal))
            .unwrap_or_else(|_| String::from("Select >= 2 bp for Tm"));
        let g4 = g4_hunter_max_window(sub, 20);

        ui.horizontal(|ui| {
            ui.strong(format!(
                "Selection: {}..={} ({} bp)",
                s0 + 1,
                e0 + 1,
                e0 - s0 + 1
            ));
            ui.separator();
            ui.label(tm_txt);
            ui.separator();
            ui.label(format!("G4Hunter: {g4:.2}"));
            if ui.button("Use Selection as Primer").clicked() {
                self.primer_seq_input = sub.to_string();
                self.primer_name_input = format!("Prm_{}_{}", s0 + 1, e0 + 1);
            }
        });

        let preview_end = (s0 + 78).min(n);
        let top_dna = &self.record.sequence[s0..preview_end];
        ui.monospace(format!("5' {:>5} {} 3'", s0 + 1, top_dna));
        ui.monospace(format!("   Frame+1 AA: {aa}"));
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
                let out_name = format!("{}_gel.svg", self.record.name);
                if fs::write(&out_name, svg).is_ok() {
                    self.status_bar = format!("Saved Virtual Gel SVG -> ./{out_name}");
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
                let _ = fs::write("osabt_gibson_primers.tsv", gplan.to_order_tsv());
                self.status_bar = String::from("Saved ./osabt_gibson_primers.tsv");
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
            ui.add(egui::TextEdit::singleline(&mut self.trace_path_input).desired_width(220.0));
            if ui.button("Align Read").clicked() {
                match fs::read(&self.trace_path_input) {
                    Ok(bytes) => {
                        let trace_opt = if bytes.starts_with(b"ABIF") {
                            parse_ab1_trace(&bytes).ok()
                        } else {
                            std::str::from_utf8(&bytes)
                                .ok()
                                .and_then(|t| parse_fasta_or_raw(t).ok())
                                .and_then(|mut r| r.pop())
                                .map(|(sample_name, sequence)| {
                                    let qualities = vec![40u8; sequence.len()];
                                    SangerTrace {
                                        sample_name,
                                        sequence,
                                        qualities,
                                    }
                                })
                        };
                        if let Some(tr) = trace_opt {
                            self.last_trace_result = align_trace_to_plasmid(&self.record, &tr, 20);
                            self.status_bar = format!("Aligned trace '{}'.", tr.sample_name);
                        }
                    }
                    Err(e) => self.status_bar = format!("Failed reading trace: {e}"),
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
                    "✓ 100% sequence concordance across aligned span (zero substitutions).",
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
