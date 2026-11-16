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

//! Interactive `egui` workbench panels for the 12 core OSABT CLI/TUI tools
//! (`#1 Bench`, `#3 DwC`, `#4 Phylo`, `#5 Epi`, `#6 Lexi`, `#9 Plate`,
//! `#11 InfoEvo`, `#12 qPCR`, `#14 eDNA`, `#15 Wildlife`, `#17 CRISPR`, `#19 Survival`).

use std::collections::HashMap;
use std::fs;

use eframe::egui;
use egui::{Color32, Pos2, Rect, Sense, Shape, Stroke, Vec2};

use crate::bench_calc::{
    build_pcr_master_mix, calc_buffer_recipe, calc_oligo_resuspension, rpm_to_rcf, BufferSystem,
};
use crate::crispr_core::{
    audit_base_editor_window, demo_crispr_target_and_offtargets, design_golden_gate_crispr_oligos,
    scan_crispr_candidates, CrisprCloningVector, CrisprGuideCandidate, CrisprNuclease,
};
use crate::crispr_offtarget::{
    export_crispr_design_tsv, render_crispr_locus_svg, scan_offtargets_and_specificity,
    CrisprOffTargetHit,
};
use crate::dwc_core::{lint_and_clean_dwc, DwcCleanConfig, DwcDataset, LintSeverity};
use crate::dwc_export::{
    export_dwc_geojson, generate_dwca_eml_xml, generate_dwca_meta_xml,
    render_dwc_occurrence_map_svg, write_dwc_delimited,
};
use crate::dwc_taxonomy::{resolve_scientific_name, TaxonBackbone};
use crate::edna_diversity::{
    analyze_edna_community_diversity, render_edna_rarefaction_svg, CommunityDiversityReport,
};
use crate::edna_fastq::{
    demo_pnw_12s_edna_dataset, export_asv_table_tsv, parse_barcode_sheet_csv, parse_fastq_str,
    run_edna_metabarcoding_pipeline, DualBarcodeSample, EdnaPipelineConfig, EdnaPipelineReport,
    FastqRead,
};
use crate::epi_models::{
    demo_outbreak_incidence, estimate_cori_rt, evaluate_diagnostic_cascade,
    rogan_gladen_prevalence, simulate_seir_rk4, CoriRtConfig, DiagnosticStepInput, SeirParams,
};
use crate::epi_power::{
    design_cohort_study, design_test_negative_study, render_power_curve_svg,
    render_rt_epidemic_svg, render_seir_trajectory_svg,
};
use crate::infoevo_core::{
    compute_alignment_information, compute_clade_divergence, compute_epistasis_mip_matrix,
    demo_ecoli_promoter_and_stemloop_fasta, parse_fasta_alignment,
    simulate_kimura_selection_information, AlignmentComplexityReport, CladeDivergenceReport,
    EpistasisMatrixReport, KimuraSelectionReport, SequenceAlignment,
};
use crate::infoevo_svg::{
    render_epistasis_mip_svg, render_kimura_trajectory_svg, render_sequence_logo_svg,
};
use crate::lexi_engine::{
    build_cognate_tree, decode_biological_phrase, demo_etymology_terms, export_anki_flashcards_tsv,
    export_flashcards_markdown, search_morpheme_lexicon, DecodedPhrase,
};
use crate::phylo_core::{parse_phylo_trees, PhyloTree};
use crate::phylo_render::render_phylo_svg;
use crate::plate_curve::{
    demo_measured_96well_plate, evaluate_plate_qc_and_normalize, fit_all_plate_compounds,
    render_4pl_curves_svg,
};
use crate::plate_layout::{
    export_bench_pipetting_markdown, export_plate_pipetting_csv, generate_randomized_plate,
    render_plate_map_svg, Microplate, PlateFormat, PlatePlanConfig, WellRole,
};
use crate::primer_design::{design_pcr_primers, PrimerDesignConfig, PrimerPairCandidate};
use crate::qpcr_core::{
    analyze_qpcr_experiment, demo_qpcr_heatshock_plate_csv, parse_qpcr_table,
    render_qpcr_expression_svg, QpcrDatasetAnalysis, QpcrReactionWell,
};
use crate::qpcr_melt::{
    analyze_melt_curve, demo_sybr_melt_traces, render_melt_curves_svg, MeltCurveResult,
};
use crate::structure::score_dimer;
use crate::survival_core::{
    compute_kaplan_meier_curves, compute_log_rank_test, demo_survival_cohort_csv,
    fit_cox_proportional_hazards, parse_survival_csv, CoxPhFit, KaplanMeierCurve,
    LogRankTestResult,
};
use crate::survival_tox::{
    demo_ecotox_lc50_csv, fit_ecotox_lc50, parse_ecotox_csv, render_survival_and_lc50_svg,
    DoseResponseLink, EcotoxBioassayPoint, EcotoxLc50Report,
};
use crate::thermo::{PolymerasePreset, ReactionConditions};
use crate::wildlife_markrecap::{
    analyze_closed_capture_histories, calc_chapman_lincoln_petersen, demo_pnw_capture_history_csv,
    estimate_cjs_apparent_survival, parse_capture_history_csv, render_wildlife_population_svg,
    ChapmanPetersenEstimate, CjsSurvivalReport, ClosedCaptureSuiteReport,
};
use crate::wildlife_occupancy::{
    demo_pnw_occupancy_csv, fit_mackenzie_occupancy_suite, parse_occupancy_csv,
    plan_occupancy_survey_effort, OccupancySuiteReport, OccupancySurveyPlan,
};

const DEFAULT_COX1_TEMPLATE: &str =
    "GGTCAACAAATCATAAAGATATTGGACTTTATATTTTATTTTTGGAGCTTGAGCTGGAATAGTAGGAACATCTTTAAGTTTATTAATTCGAGCAGAATTAGGAAACCCAGGATCTTTAATTGGAGATGATCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGATTAGTTCCATTAATATTAGGAGCCCCAGATATAGC";

const PALETTE_ACCENTS: [Color32; 5] = [
    Color32::from_rgb(56, 189, 248),
    Color32::from_rgb(52, 211, 153),
    Color32::from_rgb(250, 204, 21),
    Color32::from_rgb(248, 113, 113),
    Color32::from_rgb(192, 132, 252),
];

fn draw_chart_frame(ui: &mut egui::Ui, width: f32, height: f32, title: &str) -> (Rect, egui::Painter) {
    let (resp, painter) = ui.allocate_painter(Vec2::new(width.max(280.0), height.max(180.0)), Sense::hover());
    let outer = resp.rect;
    painter.rect_filled(outer, 6.0, Color32::from_rgb(13, 20, 36));
    painter.rect_stroke(outer, 6.0, Stroke::new(1.0_f32, Color32::from_rgb(51, 65, 85)));
    painter.text(
        Pos2::new(outer.center().x, outer.top() + 8.0),
        egui::Align2::CENTER_TOP,
        title,
        egui::FontId::proportional(13.0),
        Color32::from_rgb(226, 232, 240),
    );
    let inner = Rect::from_min_max(
        Pos2::new(outer.left() + 48.0, outer.top() + 30.0),
        Pos2::new(outer.right() - 18.0, outer.bottom() - 28.0),
    );
    (inner, painter)
}

fn draw_grid_axes(
    painter: &egui::Painter,
    plot: Rect,
    x_ticks: usize,
    y_ticks: usize,
    x_range: (f64, f64),
    y_range: (f64, f64),
    x_fmt: impl Fn(f64) -> String,
    y_fmt: impl Fn(f64) -> String,
) {
    let grid = Stroke::new(0.8_f32, Color32::from_rgb(30, 41, 66));
    let axis = Color32::from_rgb(148, 163, 184);
    for i in 0..=x_ticks {
        let f = i as f32 / x_ticks.max(1) as f32;
        let x = plot.left() + f * plot.width();
        painter.line_segment([Pos2::new(x, plot.top()), Pos2::new(x, plot.bottom())], grid);
        let val = x_range.0 + (x_range.1 - x_range.0) * f64::from(f);
        painter.text(
            Pos2::new(x, plot.bottom() + 4.0),
            egui::Align2::CENTER_TOP,
            x_fmt(val),
            egui::FontId::monospace(10.0),
            axis,
        );
    }
    for i in 0..=y_ticks {
        let f = i as f32 / y_ticks.max(1) as f32;
        let y = plot.bottom() - f * plot.height();
        painter.line_segment([Pos2::new(plot.left(), y), Pos2::new(plot.right(), y)], grid);
        let val = y_range.0 + (y_range.1 - y_range.0) * f64::from(f);
        painter.text(
            Pos2::new(plot.left() - 6.0, y),
            egui::Align2::RIGHT_CENTER,
            y_fmt(val),
            egui::FontId::monospace(10.0),
            axis,
        );
    }
    painter.rect_stroke(plot, 0.0, Stroke::new(1.0_f32, Color32::from_rgb(71, 85, 105)));
}

// ---------------------------------------------------------------------------
// #1 Bench Panel
// ---------------------------------------------------------------------------
pub struct BenchPanelState {
    template_seq: String,
    polymerase: PolymerasePreset,
    target_tm_c: f64,
    include_probe: bool,
    candidates: Vec<PrimerPairCandidate>,
    selected_cand: usize,
    rxn_count: usize,
    buffer_system: BufferSystem,
    buffer_ph: f64,
    buffer_temp_c: f64,
    rotor_rpm: f64,
    oligo_nmol: f64,
    status_msg: String,
}

impl BenchPanelState {
    pub fn new() -> Self {
        let mut s = Self {
            template_seq: DEFAULT_COX1_TEMPLATE.to_string(),
            polymerase: PolymerasePreset::Q5,
            target_tm_c: 60.0,
            include_probe: false,
            candidates: Vec::new(),
            selected_cand: 0,
            rxn_count: 12,
            buffer_system: BufferSystem::Tris,
            buffer_ph: 8.0,
            buffer_temp_c: 25.0,
            rotor_rpm: 12_000.0,
            oligo_nmol: 25.0,
            status_msg: String::from("Ready — adjust template, polymerase, or target Tm to redesign PCR pairs."),
        };
        s.recompute();
        s
    }

    fn recompute(&mut self) {
        let cfg = PrimerDesignConfig {
            amplicon_len_range: (90, 240),
            target_tm_c: self.target_tm_c,
            include_qpcr_probe: self.include_probe,
            polymerase: self.polymerase,
            num_return: 8,
            ..PrimerDesignConfig::default()
        };
        match design_pcr_primers(&self.template_seq, &cfg, None) {
            Ok(list) => {
                self.candidates = list;
                self.selected_cand = 0;
                self.status_msg = format!(
                    "Designed {} primer pairs ({:?}, target Tm {:.1} C).",
                    self.candidates.len(),
                    self.polymerase,
                    self.target_tm_c
                );
            }
            Err(e) => {
                self.candidates.clear();
                self.status_msg = format!("Primer design error: {e}");
            }
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("bench_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("bench_gui_side").default_width(330.0).show(ctx, |ui| {
            egui::ScrollArea::vertical().show(ui, |ui| {
                ui.heading("1. Primer & Master Mix Controls");
                let mut changed = false;
                ui.horizontal_wrapped(|ui| {
                    ui.label("Polymerase:");
                    for p in [
                        PolymerasePreset::Q5,
                        PolymerasePreset::Phusion,
                        PolymerasePreset::Taq,
                        PolymerasePreset::OneTaq,
                    ] {
                        if ui.selectable_value(&mut self.polymerase, p, format!("{p:?}")).changed() {
                            changed = true;
                        }
                    }
                });
                changed |= ui.add(egui::Slider::new(&mut self.target_tm_c, 48.0..=72.0).text("Target Tm (C)")).changed();
                changed |= ui.checkbox(&mut self.include_probe, "Design TaqMan qPCR Probe").changed();
                ui.label("Template DNA (A/C/G/T):");
                if ui.add(egui::TextEdit::multiline(&mut self.template_seq).desired_rows(3).desired_width(f32::INFINITY)).changed() {
                    changed = true;
                }
                if ui.button("Reset Demo COI Template").clicked() {
                    self.template_seq = DEFAULT_COX1_TEMPLATE.to_string();
                    changed = true;
                }
                if changed {
                    self.recompute();
                }

                ui.separator();
                ui.heading("2. Master Mix & Buffers");
                ui.add(egui::Slider::new(&mut self.rxn_count, 1..=96).text("Reactions"));
                ui.horizontal_wrapped(|ui| {
                    ui.label("Buffer:");
                    for (bs, ph, lbl) in [
                        (BufferSystem::Tris, 8.0, "Tris"),
                        (BufferSystem::Hepes, 7.4, "HEPES"),
                        (BufferSystem::Mes, 6.1, "MES"),
                        (BufferSystem::SodiumPhosphate, 7.2, "Na-Phos"),
                    ] {
                        if ui.selectable_label(self.buffer_system == bs, lbl).clicked() {
                            self.buffer_system = bs;
                            self.buffer_ph = ph;
                        }
                    }
                });
                ui.add(egui::Slider::new(&mut self.buffer_ph, 5.5..=9.0).text("Buffer pH"));
                ui.add(egui::Slider::new(&mut self.buffer_temp_c, 4.0..=65.0).text("Temp (C)"));
                ui.add(egui::Slider::new(&mut self.rotor_rpm, 1_000.0..=20_000.0).text("Rotor RPM"));
                ui.add(egui::Slider::new(&mut self.oligo_nmol, 5.0..=100.0).text("Dry Oligo (nmol)"));

                if let Ok(buf) = calc_buffer_recipe(self.buffer_system, self.buffer_ph, 50.0, 1000.0, self.buffer_temp_c) {
                    ui.monospace(format!(
                        "{:?}: {:.2}g acid + {:.2}g base | pKa(T)={:.2}",
                        buf.system,
                        buf.acid_grams,
                        buf.base_grams,
                        buf.effective_pka
                    ));
                }
                if let Ok(rcf) = rpm_to_rcf(self.rotor_rpm, 8.5) {
                    ui.monospace(format!("Rotor r=8.5 cm: {rcf:.0} x g (RCF)"));
                }
                let sample_oligo = self.template_seq.chars().filter(|c| !c.is_whitespace()).take(20).collect::<String>();
                let sample_seq = if sample_oligo.is_empty() { "ATGGCTAGCGATCGATCGAT" } else { sample_oligo.as_str() };
                if let Ok(res) = calc_oligo_resuspension("Oligo", sample_seq, self.oligo_nmol, 100.0, 100.0) {
                    ui.monospace(format!("100 uM stock: add {:.1} uL TE", res.add_buffer_ul));
                }

                ui.separator();
                if ui.button("Export Bench Checklist (.md)").clicked() {
                    let (ta, bp) = self
                        .candidates
                        .get(self.selected_cand)
                        .map(|c| (c.recommended_ta_c, c.amplicon_len_bp))
                        .unwrap_or((62.0, 220));
                    if let Ok(plan) = build_pcr_master_mix(self.polymerase, self.rxn_count, 25.0, 0.10, 1.0, ta, bp, false) {
                        let md = plan.to_markdown_sheet("OSABT Bench Workbench");
                        if fs::write("osabt_bench_sheet.md", md).is_ok() {
                            self.status_msg = String::from("Exported ./osabt_bench_sheet.md");
                        }
                    }
                }
            });
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let w = ui.available_width();
            let (plot, p) = draw_chart_frame(ui, w, 200.0, "Amplicon Binding Map & Primer Tm Profile");
            let seq_len = self.template_seq.len().max(1) as f32;
            p.line_segment(
                [Pos2::new(plot.left(), plot.center().y), Pos2::new(plot.right(), plot.center().y)],
                Stroke::new(3.0_f32, Color32::from_rgb(100, 116, 139)),
            );
            for (i, c) in self.candidates.iter().take(6).enumerate() {
                let y = plot.top() + 20.0 + i as f32 * 22.0;
                let x0 = plot.left() + (c.fwd_start_0 as f32 / seq_len) * plot.width();
                let x1 = plot.left() + (c.rev_end_0 as f32 / seq_len) * plot.width();
                let col = if i == self.selected_cand {
                    Color32::from_rgb(250, 204, 21)
                } else {
                    PALETTE_ACCENTS[i % PALETTE_ACCENTS.len()]
                };
                p.line_segment([Pos2::new(x0, y), Pos2::new(x1, y)], Stroke::new(if i == self.selected_cand { 3.5_f32 } else { 2.0_f32 }, col));
                p.text(
                    Pos2::new(x0, y - 8.0),
                    egui::Align2::LEFT_BOTTOM,
                    format!("#{} {}bp Ta={:.1}C", c.rank, c.amplicon_len_bp, c.recommended_ta_c),
                    egui::FontId::monospace(9.5),
                    col,
                );
            }

            ui.add_space(6.0);
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("bench_primer_grid").striped(true).show(ui, |ui| {
                    for h in ["Rank", "Span", "Amp (bp)", "Fwd Tm", "Rev Tm", "Ta (C)", "3' Dimer dG", "Forward 5'->3'", "Reverse 5'->3'"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for (idx, c) in self.candidates.iter().enumerate() {
                        if ui.selectable_label(self.selected_cand == idx, format!("#{}", c.rank)).clicked() {
                            self.selected_cand = idx;
                        }
                        ui.label(format!("{}..{}", c.fwd_start_0 + 1, c.rev_end_0 + 1));
                        ui.label(format!("{}", c.amplicon_len_bp));
                        ui.label(format!("{:.1}", c.fwd.tm_c));
                        ui.label(format!("{:.1}", c.rev.tm_c));
                        ui.label(format!("{:.1}", c.recommended_ta_c));
                        ui.label(format!("{:.2}", c.cross_dimer_3p_dg_kcal));
                        ui.monospace(&c.fwd.sequence);
                        ui.monospace(&c.rev.sequence);
                        ui.end_row();
                    }
                });
                if let Some(c) = self.candidates.get(self.selected_cand) {
                    ui.separator();
                    let cond = ReactionConditions::for_preset(self.polymerase);
                    if let Ok((any_d, d3)) = score_dimer(&c.fwd.sequence, &c.rev.sequence, &cond, 37.0) {
                        if let Some(d) = d3.or(any_d) {
                            ui.strong("Worst Forward x Reverse Dimer Alignment:");
                            ui.monospace(&d.ascii_alignment);
                        }
                    }
                }
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #3 DwC Panel
// ---------------------------------------------------------------------------
pub struct DwcPanelState {
    dataset: DwcDataset,
    backbone: TaxonBackbone,
    taxon_query: String,
    only_flagged: bool,
    status_msg: String,
}

impl DwcPanelState {
    pub fn new() -> Self {
        let backbone = TaxonBackbone::default_curated();
        let dataset = lint_and_clean_dwc(
            "PNW_Oregon_Demo_Survey",
            DwcDataset::demo_pnw_field_survey_csv(),
            &backbone,
            &DwcCleanConfig::default(),
        )
        .expect("demo DwC parses");
        Self {
            dataset,
            backbone,
            taxon_query: String::from("Pseudotsuga taxifolia (Poir.) Britton"),
            only_flagged: false,
            status_msg: String::from("Darwin Core linter & CoordinateCleaner ready."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("dwc_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        let sum = self.dataset.summary();
        egui::SidePanel::left("dwc_gui_side").default_width(330.0).show(ctx, |ui| {
            ui.heading("Darwin Core Audit Summary");
            ui.monospace(format!("Total Records  : {}", sum.total_records));
            ui.monospace(format!("Clean Records  : {}", sum.clean_records));
            ui.monospace(format!("Georeferenced  : {}", sum.georeferenced_records));
            ui.monospace(format!("Errors / Warns : {} / {}", sum.error_count, sum.warning_count));
            ui.monospace(format!("Unique Taxa    : {}", sum.unique_taxa));
            ui.checkbox(&mut self.only_flagged, "Show only flagged rows");
            ui.separator();
            ui.heading("Taxon Synonym Resolver");
            ui.text_edit_singleline(&mut self.taxon_query);
            let res = resolve_scientific_name(&self.taxon_query, &self.backbone);
            if let Some(m) = &res.matched_entry {
                ui.colored_label(
                    Color32::from_rgb(52, 211, 153),
                    format!("Accepted: {} ({})", m.accepted_name, m.taxonomic_status),
                );
                ui.small(format!("{} > {} > {}", m.kingdom, m.family, m.genus));
            } else {
                ui.colored_label(Color32::from_rgb(250, 204, 21), format!("Match: {:?}", res.match_kind));
            }
            ui.separator();
            if ui.button("Export DwC-A + GeoJSON + SVG").clicked() {
                let out_dir = std::path::PathBuf::from("osabt_dwc_bundle");
                if fs::create_dir_all(&out_dir).is_ok() {
                    let _ = fs::write(out_dir.join("occurrence.csv"), write_dwc_delimited(&self.dataset, ',', true));
                    let _ = fs::write(out_dir.join("meta.xml"), generate_dwca_meta_xml(&self.dataset, "occurrence.csv", Some("eml.xml")));
                    let _ = fs::write(out_dir.join("eml.xml"), generate_dwca_eml_xml(&self.dataset.dataset_name, "CRV Field Team", "CC-BY-4.0"));
                    let _ = fs::write(out_dir.join("occurrences.geojson"), export_dwc_geojson(&self.dataset));
                    let _ = fs::write(out_dir.join("spatial_qa_map.svg"), render_dwc_occurrence_map_svg(&self.dataset));
                    self.status_msg = String::from("Exported DwC-A bundle -> ./osabt_dwc_bundle/");
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let coords: Vec<(f64, f64, bool, &str)> = self
                .dataset
                .records
                .iter()
                .filter_map(|r| {
                    r.valid_lat_lon
                        .map(|(lat, lon)| (lon, lat, r.has_errors || r.has_warnings, r.get("scientificName")))
                })
                .collect();
            let (min_x, max_x) = coords.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p.0), hi.max(p.0)));
            let (min_y, max_y) = coords.iter().fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), p| (lo.min(p.1), hi.max(p.1)));
            let x_rng = if min_x.is_finite() { (min_x - 0.3, max_x + 0.3) } else { (-125.0, -121.0) };
            let y_rng = if min_y.is_finite() { (min_y - 0.3, max_y + 0.3) } else { (42.0, 46.0) };

            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 240.0, "Occurrence Georeference Scatter (Green=Clean, Amber=Flagged)");
            draw_grid_axes(&p, plot, 4, 4, x_rng, y_rng, |v| format!("{v:.2}"), |v| format!("{v:.2}"));
            for (lon, lat, flagged, sp) in coords {
                let fx = ((lon - x_rng.0) / (x_rng.1 - x_rng.0).max(1e-6)).clamp(0.0, 1.0) as f32;
                let fy = ((lat - y_rng.0) / (y_rng.1 - y_rng.0).max(1e-6)).clamp(0.0, 1.0) as f32;
                let pt = Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height());
                let col = if flagged { Color32::from_rgb(250, 204, 21) } else { Color32::from_rgb(52, 211, 153) };
                p.circle_filled(pt, 5.0, col);
                p.text(pt + Vec2::new(6.0, 0.0), egui::Align2::LEFT_CENTER, sp, egui::FontId::monospace(9.5), col);
            }

            ui.add_space(6.0);
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("dwc_records_grid").striped(true).show(ui, |ui| {
                    for h in ["Row", "scientificName", "eventDate", "Lat / Lon", "Status"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for r in &self.dataset.records {
                        if self.only_flagged && !r.has_errors && !r.has_warnings {
                            continue;
                        }
                        ui.label(format!("{}", r.row_number_1));
                        ui.label(r.get("scientificName"));
                        ui.label(r.get("eventDate"));
                        ui.monospace(
                            r.valid_lat_lon
                                .map(|(la, lo)| format!("{la:.4}, {lo:.4}"))
                                .unwrap_or_else(|| String::from("Invalid")),
                        );
                        let (lbl, col) = if r.has_errors {
                            ("ERROR", Color32::from_rgb(248, 113, 113))
                        } else if r.has_warnings {
                            ("WARN", Color32::from_rgb(250, 204, 21))
                        } else {
                            ("CLEAN", Color32::from_rgb(52, 211, 153))
                        };
                        ui.colored_label(col, lbl);
                        ui.end_row();
                    }
                });
                ui.separator();
                for iss in &self.dataset.issues {
                    let col = match iss.severity {
                        LintSeverity::Error => Color32::from_rgb(248, 113, 113),
                        LintSeverity::Warning => Color32::from_rgb(250, 204, 21),
                        LintSeverity::Info => Color32::from_rgb(56, 189, 248),
                    };
                    ui.colored_label(col, format!("Row {} [{}]: {}", iss.row_number_1, iss.column_term, iss.message));
                }
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #4 Phylo Panel
// ---------------------------------------------------------------------------
pub struct PhyloPanelState {
    tree: PhyloTree,
    newick_input: String,
    proportional: bool,
    ladder_asc: bool,
    search: String,
    status_msg: String,
}

impl PhyloPanelState {
    pub fn new() -> Self {
        let tree = PhyloTree::demo_vertebrate_amniotes();
        let newick_input = tree.to_newick();
        Self {
            tree,
            newick_input,
            proportional: true,
            ladder_asc: true,
            search: String::new(),
            status_msg: String::from("Loaded Vertebrate Amniote Phylogeny — click buttons to reroot, ladderize, or parse Newick."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("phylo_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("phylo_gui_side").default_width(310.0).show(ctx, |ui| {
            ui.heading("Phylogeny Controls");
            ui.checkbox(&mut self.proportional, "Proportional Branch Lengths (Phylogram)");
            ui.horizontal(|ui| {
                if ui.button("Midpoint Reroot").clicked() {
                    if let Ok(d) = self.tree.midpoint_reroot() {
                        self.status_msg = format!("Midpoint rerooted (diameter = {d:.4}).");
                    }
                }
                if ui.button("Ladderize").clicked() {
                    self.ladder_asc = !self.ladder_asc;
                    self.tree.ladderize(self.ladder_asc);
                    self.status_msg = format!("Ladderized ({}).", if self.ladder_asc { "asc" } else { "desc" });
                }
            });
            ui.horizontal(|ui| {
                ui.label("Search taxon:");
                ui.text_edit_singleline(&mut self.search);
            });
            ui.separator();
            ui.label("Newick String:");
            ui.add(egui::TextEdit::multiline(&mut self.newick_input).desired_rows(4).desired_width(f32::INFINITY));
            ui.horizontal(|ui| {
                if ui.button("Parse Newick").clicked() {
                    match parse_phylo_trees("CustomTree", &self.newick_input) {
                        Ok(mut list) if !list.is_empty() => {
                            self.tree = list.remove(0);
                            self.status_msg = String::from("Parsed custom Newick tree.");
                        }
                        Ok(_) => self.status_msg = String::from("No trees found."),
                        Err(e) => self.status_msg = format!("Parse error: {e}"),
                    }
                }
                if ui.button("Demo Tree").clicked() {
                    self.tree = PhyloTree::demo_vertebrate_amniotes();
                    self.newick_input = self.tree.to_newick();
                }
            });
            ui.separator();
            if ui.button("Export SVG + Newick").clicked() {
                let s1 = render_phylo_svg(&self.tree, self.proportional, false);
                let s2 = self.tree.to_newick();
                if fs::write("osabt_phylo_tree.svg", s1).is_ok() && fs::write("osabt_phylo_tree.nwk", s2).is_ok() {
                    self.status_msg = String::from("Exported ./osabt_phylo_tree.svg & ./osabt_phylo_tree.nwk");
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), ui.available_height() - 10.0, "Interactive Rectangular Phylogram");
            let mut tips = Vec::new();
            let mut y_map: HashMap<usize, f32> = HashMap::new();
            let mut x_map: HashMap<usize, f32> = HashMap::new();
            self.collect_tips(self.tree.root, &mut tips);
            let n_tips = tips.len().max(1);
            for (i, &nid) in tips.iter().enumerate() {
                let frac = (i as f32 + 0.5) / n_tips as f32;
                y_map.insert(nid, plot.top() + frac * plot.height());
            }
            let max_d = self.compute_depths(self.tree.root, 0.0, &mut x_map).max(1e-6);
            self.assign_internal_y(self.tree.root, &mut y_map);

            let tree_right = plot.right() - 150.0;
            let map_x = |d: f32| plot.left() + (d / max_d) * (tree_right - plot.left()).max(40.0);

            let q = self.search.trim().to_ascii_lowercase();
            for node in &self.tree.nodes {
                let Some(&ny) = y_map.get(&node.id) else { continue };
                let Some(&nx_d) = x_map.get(&node.id) else { continue };
                let nx = map_x(nx_d);
                if let Some(pid) = node.parent {
                    if let (Some(&py), Some(&px_d)) = (y_map.get(&pid), x_map.get(&pid)) {
                        let px = map_x(px_d);
                        p.line_segment([Pos2::new(px, py), Pos2::new(px, ny)], Stroke::new(2.0_f32, Color32::from_rgb(56, 189, 248)));
                        p.line_segment([Pos2::new(px, ny), Pos2::new(nx, ny)], Stroke::new(2.0_f32, Color32::from_rgb(56, 189, 248)));
                    }
                }
                if node.is_leaf() {
                    let hit = !q.is_empty() && node.name.to_ascii_lowercase().contains(&q);
                    let col = if hit { Color32::from_rgb(250, 204, 21) } else { Color32::from_rgb(226, 232, 240) };
                    p.circle_filled(Pos2::new(nx, ny), 3.6, Color32::from_rgb(52, 211, 153));
                    p.text(Pos2::new(nx + 7.0, ny), egui::Align2::LEFT_CENTER, &node.name, egui::FontId::monospace(11.0), col);
                } else if let Some(sup) = node.support_pct {
                    p.text(Pos2::new(nx - 4.0, ny - 6.0), egui::Align2::RIGHT_BOTTOM, format!("{sup:.0}"), egui::FontId::monospace(9.0), Color32::from_rgb(148, 163, 184));
                }
            }
        });
    }

    fn collect_tips(&self, nid: usize, out: &mut Vec<usize>) {
        let node = &self.tree.nodes[nid];
        if node.is_leaf() || node.collapsed {
            out.push(nid);
        } else {
            for &ch in &node.children {
                self.collect_tips(ch, out);
            }
        }
    }

    fn compute_depths(&self, nid: usize, cur: f32, x_map: &mut HashMap<usize, f32>) -> f32 {
        x_map.insert(nid, cur);
        let node = &self.tree.nodes[nid];
        if node.is_leaf() || node.collapsed {
            return cur;
        }
        let mut max_d = cur;
        for &ch in &node.children {
            let step = if self.proportional {
                self.tree.nodes[ch].branch_length.unwrap_or(0.05) as f32
            } else {
                1.0
            };
            max_d = max_d.max(self.compute_depths(ch, cur + step.max(0.002), x_map));
        }
        max_d
    }

    fn assign_internal_y(&self, nid: usize, y_map: &mut HashMap<usize, f32>) -> Option<f32> {
        if let Some(&y) = y_map.get(&nid) {
            return Some(y);
        }
        let node = &self.tree.nodes[nid];
        let ys: Vec<f32> = node.children.iter().filter_map(|&c| self.assign_internal_y(c, y_map)).collect();
        if ys.is_empty() {
            return None;
        }
        let avg = (ys.first().copied().unwrap_or(0.0) + ys.last().copied().unwrap_or(0.0)) * 0.5;
        y_map.insert(nid, avg);
        Some(avg)
    }
}

// ---------------------------------------------------------------------------
// #5 Epi Panel
// ---------------------------------------------------------------------------
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EpiPanelTab {
    RenewalRt,
    SeirSweep,
    DiagnosticBayes,
    StudyPower,
}

pub struct EpiPanelState {
    tab: EpiPanelTab,
    incidence: Vec<f64>,
    rt_cfg: CoriRtConfig,
    seir_params: SeirParams,
    prior_prev: f64,
    pos1: bool,
    pos2: bool,
    tnd_mode: bool,
    p0_control: f64,
    ve_or_rr: f64,
    target_power: f64,
    status_msg: String,
}

impl EpiPanelState {
    pub fn new() -> Self {
        Self {
            tab: EpiPanelTab::RenewalRt,
            incidence: demo_outbreak_incidence(),
            rt_cfg: CoriRtConfig::default(),
            seir_params: SeirParams::default(),
            prior_prev: 0.08,
            pos1: true,
            pos2: true,
            tnd_mode: true,
            p0_control: 0.52,
            ve_or_rr: 0.60,
            target_power: 0.80,
            status_msg: String::from("Epidemiology studio ready — switch tabs for Cori Rt, SEIR, Bayes, and Study Power."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("epi_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("epi_gui_side").default_width(300.0).show(ctx, |ui| {
            ui.heading("Epi Mode & Parameters");
            ui.selectable_value(&mut self.tab, EpiPanelTab::RenewalRt, "1. Cori Renewal Rt");
            ui.selectable_value(&mut self.tab, EpiPanelTab::SeirSweep, "2. RK4 SEIR Simulator");
            ui.selectable_value(&mut self.tab, EpiPanelTab::DiagnosticBayes, "3. Diagnostic Bayes");
            ui.selectable_value(&mut self.tab, EpiPanelTab::StudyPower, "4. Cohort & TND Power");
            ui.separator();
            match self.tab {
                EpiPanelTab::RenewalRt => {
                    ui.add(egui::Slider::new(&mut self.rt_cfg.serial_interval_mean, 1.0..=14.0).text("SI Mean (d)"));
                    ui.add(egui::Slider::new(&mut self.rt_cfg.window_days, 2..=14).text("Window (d)"));
                }
                EpiPanelTab::SeirSweep => {
                    ui.add(egui::Slider::new(&mut self.seir_params.r0, 0.5..=10.0).text("R0"));
                    ui.add(egui::Slider::new(&mut self.seir_params.npi_reduction_fraction, 0.0..=0.90).text("NPI Reduction"));
                }
                EpiPanelTab::DiagnosticBayes => {
                    ui.add(egui::Slider::new(&mut self.prior_prev, 0.005..=0.95).text("Prevalence"));
                    ui.checkbox(&mut self.pos1, "Assay 1 (Rapid Ag) Positive");
                    ui.checkbox(&mut self.pos2, "Assay 2 (RT-qPCR) Positive");
                }
                EpiPanelTab::StudyPower => {
                    ui.checkbox(&mut self.tnd_mode, "Test-Negative Design (vs Cohort)");
                    ui.add(egui::Slider::new(&mut self.p0_control, 0.05..=0.80).text("Control p0"));
                    ui.add(egui::Slider::new(&mut self.ve_or_rr, 0.10..=0.90).text("VE / Effect"));
                    ui.add(egui::Slider::new(&mut self.target_power, 0.60..=0.95).text("Power"));
                }
            }
            ui.separator();
            if ui.button("Export Active .svg").clicked() {
                match self.tab {
                    EpiPanelTab::RenewalRt => {
                        if let Ok(pts) = estimate_cori_rt(&self.incidence, &self.rt_cfg) {
                            let _ = fs::write("epi_rt_curve.svg", render_rt_epidemic_svg("Outbreak Rt", &pts));
                            self.status_msg = String::from("Exported ./epi_rt_curve.svg");
                        }
                    }
                    EpiPanelTab::SeirSweep => {
                        if let Ok(sim) = simulate_seir_rk4(&self.seir_params) {
                            let _ = fs::write("epi_seir_trajectory.svg", render_seir_trajectory_svg(&sim));
                            self.status_msg = String::from("Exported ./epi_seir_trajectory.svg");
                        }
                    }
                    EpiPanelTab::StudyPower => {
                        let plan = if self.tnd_mode {
                            design_test_negative_study(self.p0_control, self.ve_or_rr, 2.0, 0.05, self.target_power, Some(0.25))
                        } else {
                            design_cohort_study(self.p0_control.clamp(0.01, 0.45), (1.0 - self.ve_or_rr).clamp(0.1, 0.95), 2.0, 0.05, self.target_power)
                        };
                        if let Ok(pl) = plan {
                            let _ = fs::write("epi_power_curve.svg", render_power_curve_svg(&pl));
                            self.status_msg = String::from("Exported ./epi_power_curve.svg");
                        }
                    }
                    EpiPanelTab::DiagnosticBayes => {
                        self.status_msg = String::from("Diagnostic Bayes table is live; switch to Rt/SEIR/Power to export SVG.");
                    }
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| match self.tab {
            EpiPanelTab::RenewalRt => {
                let pts = estimate_cori_rt(&self.incidence, &self.rt_cfg).unwrap_or_default();
                let max_d = pts.len().max(1) as f64;
                let max_r = pts.iter().filter_map(|x| x.rt_mean).fold(3.0_f64, f64::max).max(2.0);
                let (plot, p) = draw_chart_frame(ui, ui.available_width(), 320.0, "Cori et al. (2013) Instantaneous Reproduction Number Rt");
                draw_grid_axes(&p, plot, 5, 4, (1.0, max_d), (0.0, max_r), |v| format!("d{v:.0}"), |v| format!("{v:.1}"));
                let y_one = plot.bottom() - ((1.0 / max_r) as f32) * plot.height();
                p.line_segment([Pos2::new(plot.left(), y_one), Pos2::new(plot.right(), y_one)], Stroke::new(1.2_f32, Color32::from_rgb(248, 113, 113)));
                let line: Vec<Pos2> = pts
                    .iter()
                    .filter_map(|pt| {
                        pt.rt_mean.map(|rt| {
                            let fx = (pt.day_index as f32 / (max_d as f32 - 1.0).max(1.0)).clamp(0.0, 1.0);
                            let fy = (rt / max_r).clamp(0.0, 1.0) as f32;
                            Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                        })
                    })
                    .collect();
                p.add(Shape::line(line, Stroke::new(2.4_f32, Color32::from_rgb(56, 189, 248))));
            }
            EpiPanelTab::SeirSweep => {
                if let Ok(sim) = simulate_seir_rk4(&self.seir_params) {
                    let max_d = sim.trajectory.last().map(|t| t.day).unwrap_or(180.0);
                    let pop = sim.params.population.max(1.0);
                    let (plot, p) = draw_chart_frame(ui, ui.available_width(), 320.0, "RK4 SEIR Epidemic Curves (S=Blue, E=Amber, I=Red, R=Green)");
                    draw_grid_axes(&p, plot, 6, 4, (0.0, max_d), (0.0, pop), |v| format!("{v:.0}d"), |v| format!("{:.0}k", v / 1000.0));
                    for (col, ext) in [
                        (Color32::from_rgb(56, 189, 248), (|t: &crate::epi_models::SeirTimePoint| t.susceptible) as fn(&_) -> f64),
                        (Color32::from_rgb(250, 204, 21), |t| t.exposed),
                        (Color32::from_rgb(248, 113, 113), |t| t.infectious),
                        (Color32::from_rgb(52, 211, 153), |t| t.recovered),
                    ] {
                        let pts: Vec<Pos2> = sim.trajectory.iter().map(|tp| {
                            let fx = (tp.day / max_d).clamp(0.0, 1.0) as f32;
                            let fy = (ext(tp) / pop).clamp(0.0, 1.0) as f32;
                            Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                        }).collect();
                        p.add(Shape::line(pts, Stroke::new(2.2_f32, col)));
                    }
                    ui.add_space(6.0);
                    ui.monospace(format!("HIT = {:.1}% | Peak Infectious = {:.0} on Day {:.0} | Final Attack = {:.1}%", sim.herd_immunity_threshold * 100.0, sim.peak_infectious, sim.peak_day, sim.final_attack_fraction * 100.0));
                }
            }
            EpiPanelTab::DiagnosticBayes => {
                let steps = vec![
                    DiagnosticStepInput { assay_name: String::from("Rapid Antigen"), sensitivity: 0.86, specificity: 0.985, observed_positive: self.pos1 },
                    DiagnosticStepInput { assay_name: String::from("Confirmatory qPCR"), sensitivity: 0.975, specificity: 0.995, observed_positive: self.pos2 },
                ];
                ui.heading("Sequential Diagnostic Likelihood Cascade");
                if let Ok(casc) = evaluate_diagnostic_cascade(self.prior_prev, &steps) {
                    egui::Grid::new("epi_bayes_grid").striped(true).show(ui, |ui| {
                        for h in ["Assay", "Result", "LR+", "LR-", "Pre-Test", "Post-Test", "PPV", "NPV"] {
                            ui.strong(h);
                        }
                        ui.end_row();
                        for s in &casc {
                            ui.label(&s.assay_name);
                            ui.label(if s.observed_positive { "POS +" } else { "NEG -" });
                            ui.label(format!("{:.1}", s.lr_plus));
                            ui.label(format!("{:.3}", s.lr_minus));
                            ui.label(format!("{:.2}%", s.pre_test_prob * 100.0));
                            ui.strong(format!("{:.2}%", s.post_test_prob * 100.0));
                            ui.label(format!("{:.1}%", s.stage_ppv * 100.0));
                            ui.label(format!("{:.2}%", s.stage_npv * 100.0));
                            ui.end_row();
                        }
                    });
                }
                if let Ok(rg) = rogan_gladen_prevalence(120, 1000, 0.86, 0.985) {
                    ui.separator();
                    ui.monospace(format!(
                        "Rogan-Gladen (120/1000 positive): Apparent={:.1}% -> True={:.2}% (95% CI [{:.2}%, {:.2}%])",
                        rg.apparent_prevalence * 100.0,
                        rg.true_prevalence * 100.0,
                        rg.true_ci95_clamped.0 * 100.0,
                        rg.true_ci95_clamped.1 * 100.0
                    ));
                }
            }
            EpiPanelTab::StudyPower => {
                let plan = if self.tnd_mode {
                    design_test_negative_study(self.p0_control, self.ve_or_rr, 2.0, 0.05, self.target_power, Some(0.25))
                } else {
                    design_cohort_study(self.p0_control.clamp(0.01, 0.45), (1.0 - self.ve_or_rr).clamp(0.1, 0.95), 2.0, 0.05, self.target_power)
                };
                if let Ok(pl) = plan {
                    let max_n = pl.power_curve.last().map(|x| x.0 as f64).unwrap_or(500.0).max(10.0);
                    let (plot, p) = draw_chart_frame(ui, ui.available_width(), 320.0, pl.kind.title());
                    draw_grid_axes(&p, plot, 5, 5, (0.0, max_n), (0.0, 1.0), |v| format!("{v:.0}"), |v| format!("{:.0}%", v * 100.0));
                    let pts: Vec<Pos2> = pl.power_curve.iter().map(|&(n, pw)| {
                        let fx = (n as f64 / max_n).clamp(0.0, 1.0) as f32;
                        let fy = pw.clamp(0.0, 1.0) as f32;
                        Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                    }).collect();
                    p.add(Shape::line(pts, Stroke::new(2.5_f32, Color32::from_rgb(52, 211, 153))));
                    ui.add_space(6.0);
                    ui.strong(format!("Fleiss CC Sample Size: N1={} + N0={} = {} total", pl.n1_fleiss_cc, pl.n0_fleiss_cc, pl.total_fleiss_cc));
                }
            }
        });
    }
}

// ---------------------------------------------------------------------------
// #6 Lexi Panel
// ---------------------------------------------------------------------------
pub struct LexiPanelState {
    input_phrase: String,
    lexicon_filter: String,
    saved_deck: Vec<DecodedPhrase>,
    quiz_idx: usize,
    quiz_revealed: bool,
    status_msg: String,
}

impl LexiPanelState {
    pub fn new() -> Self {
        let presets = demo_etymology_terms();
        let saved_deck: Vec<DecodedPhrase> = presets.iter().map(|t| decode_biological_phrase(t)).collect();
        Self {
            input_phrase: presets[0].to_string(),
            lexicon_filter: String::new(),
            saved_deck,
            quiz_idx: 0,
            quiz_revealed: false,
            status_msg: String::from("Type any binomial, clade, or anatomical term to segment Greek/Latin morphemes."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("lexi_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("lexi_gui_side").default_width(310.0).show(ctx, |ui| {
            ui.heading("Etymology & Flashcard Deck");
            ui.horizontal_wrapped(|ui| {
                for preset in demo_etymology_terms() {
                    if ui.button(preset).clicked() {
                        self.input_phrase = preset.to_string();
                    }
                }
            });
            ui.separator();
            ui.heading(format!("Study Flashcard Quiz ({} cards)", self.saved_deck.len()));
            if let Some(card) = self.saved_deck.get(self.quiz_idx) {
                ui.strong(&card.raw_input);
                if self.quiz_revealed {
                    ui.colored_label(Color32::from_rgb(52, 211, 153), &card.combined_gloss);
                } else {
                    ui.small("(Click Reveal Meaning)");
                }
                ui.horizontal(|ui| {
                    if ui.button("Reveal/Hide").clicked() {
                        self.quiz_revealed = !self.quiz_revealed;
                    }
                    if ui.button("Next Card").clicked() && !self.saved_deck.is_empty() {
                        self.quiz_idx = (self.quiz_idx + 1) % self.saved_deck.len();
                        self.quiz_revealed = false;
                    }
                });
            }
            ui.separator();
            if ui.button("Add Current Term to Deck").clicked() {
                let d = decode_biological_phrase(&self.input_phrase);
                if !d.raw_input.is_empty() {
                    self.saved_deck.push(d);
                    self.status_msg = format!("Added '{}' (deck = {} cards).", self.input_phrase, self.saved_deck.len());
                }
            }
            if ui.button("Export Anki .tsv + Study .md").clicked() {
                let tsv = export_anki_flashcards_tsv(&self.saved_deck, "OSABT_Etymology");
                let md = export_flashcards_markdown(&self.saved_deck);
                if fs::write("osabt_etymology_anki.tsv", tsv).is_ok() && fs::write("osabt_etymology_study.md", md).is_ok() {
                    self.status_msg = String::from("Exported ./osabt_etymology_anki.tsv & ./osabt_etymology_study.md");
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.strong("Query Term:");
                ui.add(egui::TextEdit::singleline(&mut self.input_phrase).desired_width(360.0));
            });
            let decoded = decode_biological_phrase(&self.input_phrase);
            ui.colored_label(Color32::from_rgb(56, 189, 248), format!("Literal Synthesis: {}", decoded.combined_gloss));

            // Morpheme Painter Diagram
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 120.0, "Morpheme Segmentation Blocks");
            let segs: Vec<_> = decoded.words.iter().flat_map(|w| w.segments.iter()).collect();
            if !segs.is_empty() {
                let bw = (plot.width() / segs.len() as f32).clamp(40.0, 180.0);
                for (i, seg) in segs.iter().enumerate() {
                    let x0 = plot.left() + i as f32 * bw + 4.0;
                    let r = Rect::from_min_max(Pos2::new(x0, plot.top() + 8.0), Pos2::new((x0 + bw - 8.0).min(plot.right()), plot.bottom() - 8.0));
                    let col = PALETTE_ACCENTS[i % PALETTE_ACCENTS.len()];
                    p.rect_stroke(r, 4.0, Stroke::new(1.5_f32, col));
                    p.text(Pos2::new(r.center().x, r.top() + 14.0), egui::Align2::CENTER_CENTER, &seg.surface_slice, egui::FontId::monospace(12.0), col);
                    p.text(Pos2::new(r.center().x, r.bottom() - 14.0), egui::Align2::CENTER_CENTER, &seg.gloss, egui::FontId::proportional(10.0), Color32::from_rgb(203, 213, 225));
                }
            }

            ui.add_space(6.0);
            let cognates = build_cognate_tree(&decoded);
            ui.strong(format!("Cognate Network ({} roots):", cognates.len()));
            for cn in cognates.iter().take(6) {
                ui.label(format!("• {} [{} — {}] '{}': {}", cn.morpheme_id, cn.origin, cn.classical_lemma, cn.gloss, cn.cognates.join(", ")));
            }
            ui.separator();
            ui.horizontal(|ui| {
                ui.label("Search Lexicon:");
                ui.text_edit_singleline(&mut self.lexicon_filter);
            });
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("lexi_root_grid").striped(true).show(ui, |ui| {
                    for h in ["Root", "Origin", "Lemma", "Gloss", "Cognates"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for e in search_morpheme_lexicon(&self.lexicon_filter).into_iter().take(30) {
                        ui.monospace(e.canonical_id);
                        ui.label(e.origin.label());
                        ui.label(e.classical_lemma);
                        ui.label(e.gloss);
                        ui.small(e.cognates.join(", "));
                        ui.end_row();
                    }
                });
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #9 Plate Panel
// ---------------------------------------------------------------------------
pub struct PlatePanelState {
    cfg: PlatePlanConfig,
    plate: Microplate,
    status_msg: String,
}

impl PlatePanelState {
    pub fn new() -> Self {
        let cfg = PlatePlanConfig::default();
        let plate = demo_measured_96well_plate();
        Self {
            cfg,
            plate,
            status_msg: String::from("96/384-Well Edge-Buffered Assay & 4PL IC50/EC50 Studio ready."),
        }
    }

    fn rerandomize(&mut self, seed: u64) {
        self.cfg.rng_seed = seed;
        if let Ok(mut p) = generate_randomized_plate(&self.cfg) {
            let mut st = seed ^ 0xA5A5_5A5A_F0F0_0F0F;
            let mut eps = || -> f64 {
                st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
                (((st >> 32) as f64) / (u32::MAX as f64) - 0.5) * 1400.0
            };
            for w in &mut p.wells {
                match w.role {
                    WellRole::NegativeControl => w.raw_readout = Some(52_000.0 + eps()),
                    WellRole::PositiveControl => w.raw_readout = Some(1_200.0 + eps() * 0.1),
                    WellRole::CompoundDose => {
                        let c = w.concentration_nm.unwrap_or(1.0);
                        let (ic50, hill): (f64, f64) = match w.compound_name.as_str() {
                            "Staurosporine" => (4.2, -1.04),
                            "CRV-104" => (18.5, -1.15),
                            _ => (420.0, -0.96),
                        };
                        let u = 10_f64.powf((ic50.log10() - c.log10()) * hill);
                        let pct = 2.0 + 97.5 / (1.0 + u);
                        w.raw_readout = Some((1_200.0 + (pct / 100.0) * 50_800.0 + eps() * 0.6).max(80.0));
                    }
                    _ => {}
                }
            }
            let _ = evaluate_plate_qc_and_normalize(&mut p);
            self.plate = p;
            self.status_msg = format!("Generated {} randomized plate (seed {}).", self.plate.format.label(), seed);
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("plate_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("plate_gui_side").default_width(300.0).show(ctx, |ui| {
            ui.heading("Plate Format & QC");
            ui.horizontal(|ui| {
                if ui.selectable_label(self.cfg.format == PlateFormat::Well96, "96-Well").clicked() {
                    self.cfg.format = PlateFormat::Well96;
                    self.cfg.assay_volume_ul = 100.0;
                    let s = self.cfg.rng_seed;
                    self.rerandomize(s);
                }
                if ui.selectable_label(self.cfg.format == PlateFormat::Well384, "384-Well").clicked() {
                    self.cfg.format = PlateFormat::Well384;
                    self.cfg.assay_volume_ul = 40.0;
                    let s = self.cfg.rng_seed;
                    self.rerandomize(s);
                }
            });
            if ui.button("Re-Randomize Plate Layout").clicked() {
                let s = self.cfg.rng_seed.wrapping_add(7919);
                self.rerandomize(s);
            }
            ui.separator();
            let mut cp = self.plate.clone();
            if let Ok(qc) = evaluate_plate_qc_and_normalize(&mut cp) {
                ui.strong(qc.quality_grade);
                ui.monospace(format!("Z'-Factor : {:.4}", qc.z_prime_factor));
                ui.monospace(format!("SSMD      : {:.2}", qc.ssmd));
                ui.monospace(format!("S/B Ratio : {:.1}x", qc.signal_to_background));
                ui.monospace(format!("Neg CV    : {:.2}% (n={})", qc.neg_cv_pct, qc.neg_ctrl_n));
                ui.monospace(format!("Pos CV    : {:.2}% (n={})", qc.pos_cv_pct, qc.pos_ctrl_n));
            }
            ui.separator();
            if ui.button("Export CSV + MD + SVGs").clicked() {
                let fits = fit_all_plate_compounds(&self.plate);
                let _ = fs::write("osabt_plate_layout.csv", export_plate_pipetting_csv(&self.plate));
                let _ = fs::write("osabt_bench_pipetting.md", export_bench_pipetting_markdown(&self.plate, &self.cfg));
                let _ = fs::write("osabt_plate_heatmap.svg", render_plate_map_svg(&self.plate));
                let _ = fs::write("osabt_4pl_curves.svg", render_4pl_curves_svg(&self.plate.name, &fits));
                self.status_msg = String::from("Exported plate CSV, MD, heatmap SVG & 4PL curves SVG.");
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let fits = fit_all_plate_compounds(&self.plate);
            ui.columns(2, |cols| {
                let w0 = cols[0].available_width();
                let (plot, p) = draw_chart_frame(&mut cols[0], w0, 280.0, &format!("{} Well Map", self.plate.format.label()));
                let (nr, nc) = self.plate.format.dims();
                let cw = plot.width() / nc as f32;
                let ch = plot.height() / nr as f32;
                let rad = (cw.min(ch) * 0.38).max(2.5);
                for w in &self.plate.wells {
                    let cx = plot.left() + (w.col_idx as f32 + 0.5) * cw;
                    let cy = plot.top() + (w.row_idx as f32 + 0.5) * ch;
                    let col = match w.role {
                        WellRole::BufferMoat => Color32::from_rgb(59, 130, 246),
                        WellRole::NegativeControl => Color32::from_rgb(34, 197, 94),
                        WellRole::PositiveControl => Color32::from_rgb(239, 68, 68),
                        WellRole::CompoundDose => {
                            let pct = (w.normalized_pct.unwrap_or(50.0) / 100.0).clamp(0.0, 1.0) as f32;
                            Color32::from_rgb((250.0 * (1.0 - pct)) as u8 + 40, (200.0 * pct) as u8 + 40, 220)
                        }
                        WellRole::Empty => Color32::from_rgb(51, 65, 85),
                    };
                    p.circle_filled(Pos2::new(cx, cy), rad, col);
                }

                let w1 = cols[1].available_width();
                let (p2, pt2) = draw_chart_frame(&mut cols[1], w1, 280.0, "4PL Dose-Response Curves (% Signal)");
                draw_grid_axes(&pt2, p2, 4, 4, (-1.0, 4.0), (0.0, 110.0), |v| format!("1e{v:.0}"), |v| format!("{v:.0}%"));
                for (fi, fit) in fits.iter().enumerate() {
                    let col = PALETTE_ACCENTS[fi % PALETTE_ACCENTS.len()];
                    let line: Vec<Pos2> = (0..=40).map(|k| {
                        let lx = -1.0 + 5.0 * (k as f64 / 40.0);
                        let resp = fit.bottom + (fit.top - fit.bottom) / (1.0 + 10_f64.powf((fit.log10_ec50 - lx) * fit.hill_slope));
                        let fx = ((lx + 1.0) / 5.0) as f32;
                        let fy = (resp / 110.0).clamp(0.0, 1.0) as f32;
                        Pos2::new(p2.left() + fx * p2.width(), p2.bottom() - fy * p2.height())
                    }).collect();
                    pt2.add(Shape::line(line, Stroke::new(2.2_f32, col)));
                }
            });

            ui.add_space(6.0);
            egui::Grid::new("plate_4pl_grid").striped(true).show(ui, |ui| {
                for h in ["Compound", "Rel IC50 (nM)", "95% CI (nM)", "IC90 (nM)", "Hill Slope", "Bottom/Top", "R^2"] {
                    ui.strong(h);
                }
                ui.end_row();
                for f in &fits {
                    ui.label(&f.compound_name);
                    ui.strong(format!("{:.2}", f.relative_ec50_nm));
                    ui.label(format!("[{:.2}, {:.2}]", f.relative_ec50_ci95_nm.0, f.relative_ec50_ci95_nm.1));
                    ui.label(format!("{:.1}", f.ec90_nm));
                    ui.label(format!("{:+.2}", f.hill_slope));
                    ui.label(format!("{:.1}% / {:.1}%", f.bottom, f.top));
                    ui.label(format!("{:.4}", f.r_squared));
                    ui.end_row();
                }
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #11 InfoEvo Panel
// ---------------------------------------------------------------------------
pub struct InfoEvoPanelState {
    aln: SequenceAlignment,
    small_sample_corr: bool,
    info: AlignmentComplexityReport,
    epistasis: EpistasisMatrixReport,
    divergence: Option<CladeDivergenceReport>,
    p0_fittest: f64,
    sel_coeff_s: f64,
    kimura: KimuraSelectionReport,
    status_msg: String,
}

impl InfoEvoPanelState {
    pub fn new() -> Self {
        let aln = parse_fasta_alignment("Ecoli_Pribnow_RNA_StemLoop_Demo", demo_ecoli_promoter_and_stemloop_fasta()).expect("demo infoevo fasta");
        let info = compute_alignment_information(&aln, true, Some((4_600_000.0, 4_500.0)));
        let epistasis = compute_epistasis_mip_matrix(&aln, 0.30, 20);
        let divergence = compute_clade_divergence(&aln, "", "").ok();
        let p0_fittest = 0.01;
        let sel_coeff_s = 0.12;
        let rem = 1.0 - p0_fittest;
        let kimura = simulate_kimura_selection_information(&[p0_fittest, 0.35 * rem, 0.65 * rem], &[sel_coeff_s, 0.35 * sel_coeff_s, 0.0], 1.0, 240).expect("valid kimura");
        Self {
            aln,
            small_sample_corr: true,
            info,
            epistasis,
            divergence,
            p0_fittest,
            sel_coeff_s,
            kimura,
            status_msg: String::from("Schneider Sequence Logo, Adami Complexity, MIp Epistasis & Kimura Regret ready."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("infoevo_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("infoevo_gui_side").default_width(300.0).show(ctx, |ui| {
            ui.heading("ITME & Logo Controls");
            if ui.checkbox(&mut self.small_sample_corr, "Basharin/Schneider e_n Correction").changed() {
                self.info = compute_alignment_information(&self.aln, self.small_sample_corr, Some((4_600_000.0, 4_500.0)));
            }
            let mut k_chg = false;
            k_chg |= ui.add(egui::Slider::new(&mut self.p0_fittest, 0.001..=0.40).logarithmic(true).text("Kimura p0*")).changed();
            k_chg |= ui.add(egui::Slider::new(&mut self.sel_coeff_s, 0.02..=0.50).text("Selection s")).changed();
            if k_chg {
                let rem = 1.0 - self.p0_fittest;
                if let Ok(k) = simulate_kimura_selection_information(&[self.p0_fittest, 0.35 * rem, 0.65 * rem], &[self.sel_coeff_s, 0.35 * self.sel_coeff_s, 0.0], 1.0, 240) {
                    self.kimura = k;
                }
            }
            ui.separator();
            ui.monospace(format!("R_sequence  : {:.2} bits", self.info.total_r_sequence_bits));
            ui.monospace(format!("Adami C     : {:.2} bits", self.info.adami_physical_complexity_bits));
            ui.monospace(format!("Density     : {:.1}%", self.info.complexity_density * 100.0));
            ui.monospace(format!("Kimura Lim  : {:.2} bits", self.kimura.kimura_limit_bits));
            ui.monospace(format!("Info Gain I : {:.2} bits", self.kimura.final_info_gain_bits));
            ui.monospace(format!("Regret L_sub: {:.2} bits", self.kimura.final_substitutional_load_bits));
            if let Some(div) = &self.divergence {
                ui.monospace(format!("Clade D_JS  : {:.2} bits", div.total_js_bits));
            }
            ui.separator();
            if ui.button("Export 3 .svg Figures").clicked() {
                let _ = fs::write("infoevo_sequence_logo.svg", render_sequence_logo_svg(&self.info, &self.aln.title));
                let _ = fs::write("infoevo_epistasis_mip.svg", render_epistasis_mip_svg(&self.epistasis, &self.info, &self.aln.title));
                let _ = fs::write("infoevo_kimura_selection.svg", render_kimura_trajectory_svg(&self.kimura, "Kimura Adaptive Sweep"));
                self.status_msg = String::from("Exported infoevo_sequence_logo.svg, epistasis_mip.svg, kimura_selection.svg");
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 220.0, "Schneider & Stephens (1990) Sequence Logo (bits)");
            let max_b = self.info.max_bits_per_site.max(2.0) as f32;
            let n = self.info.sites.len().min(48).max(1);
            let col_w = plot.width() / n as f32;
            for (idx, site) in self.info.sites.iter().take(n).enumerate() {
                let x = plot.left() + (idx as f32 + 0.5) * col_w;
                let mut y_bot = plot.bottom();
                for &(ch, hb) in &site.letter_heights_bits {
                    let h = ((hb as f32) / max_b) * plot.height();
                    if h > 2.0 {
                        let col = match ch.to_ascii_uppercase() {
                            'A' => Color32::from_rgb(52, 211, 153),
                            'C' => Color32::from_rgb(56, 189, 248),
                            'G' => Color32::from_rgb(250, 204, 21),
                            'T' | 'U' => Color32::from_rgb(248, 113, 113),
                            _ => Color32::from_rgb(192, 132, 252),
                        };
                        p.text(Pos2::new(x, y_bot - h * 0.5), egui::Align2::CENTER_CENTER, ch.to_string(), egui::FontId::monospace(h.clamp(8.0, 22.0)), col);
                        y_bot -= h;
                    }
                }
            }

            ui.add_space(6.0);
            let max_g = self.kimura.generations.max(1) as f64;
            let max_bits = self.kimura.final_substitutional_load_bits.max(self.kimura.kimura_limit_bits).max(2.0) * 1.1;
            let (p2, pt2) = draw_chart_frame(ui, ui.available_width(), 200.0, "Kimura Selection Information Gain I_t (Green) vs Substitutional Load L_sub (Amber)");
            draw_grid_axes(&pt2, p2, 5, 4, (0.0, max_g), (0.0, max_bits), |v| format!("g{v:.0}"), |v| format!("{v:.1}b"));
            let line_i: Vec<Pos2> = self.kimura.steps.iter().map(|s| Pos2::new(p2.left() + (s.generation as f32 / max_g as f32) * p2.width(), p2.bottom() - (s.kl_info_gain_bits / max_bits).clamp(0.0, 1.0) as f32 * p2.height())).collect();
            let line_l: Vec<Pos2> = self.kimura.steps.iter().map(|s| Pos2::new(p2.left() + (s.generation as f32 / max_g as f32) * p2.width(), p2.bottom() - (s.cumulative_substitutional_load_bits / max_bits).clamp(0.0, 1.0) as f32 * p2.height())).collect();
            pt2.add(Shape::line(line_l, Stroke::new(2.0_f32, Color32::from_rgb(250, 204, 21))));
            pt2.add(Shape::line(line_i, Stroke::new(2.4_f32, Color32::from_rgb(52, 211, 153))));

            ui.add_space(6.0);
            ui.strong("Top APC-Corrected Epistatic Coevolution Pairs (Dunn et al. 2008 MIp):");
            for pair in self.epistasis.top_pairs.iter().take(4) {
                ui.monospace(format!("Sites ({}, {}): MIp={:.3} bits (Z={:+.2}) | {}", pair.pos_i_1based, pair.pos_j_1based, pair.mip_bits, pair.z_score, pair.top_covarying_states));
            }
        });
    }
}

// ---------------------------------------------------------------------------
// #12 qPCR Panel
// ---------------------------------------------------------------------------
pub struct QpcrPanelState {
    wells: Vec<QpcrReactionWell>,
    include_b2m: bool,
    analysis: QpcrDatasetAnalysis,
    melt_results: Vec<MeltCurveResult>,
    status_msg: String,
}

impl QpcrPanelState {
    pub fn new() -> Self {
        let wells = parse_qpcr_table(demo_qpcr_heatshock_plate_csv()).expect("demo qpcr");
        let refs = vec![String::from("ACTB"), String::from("GAPDH")];
        let analysis = analyze_qpcr_experiment(&wells, &refs, "Control_37C").expect("demo analysis");
        let melt_results: Vec<MeltCurveResult> = demo_sybr_melt_traces().iter().map(|t| analyze_melt_curve(t, 77.5)).collect();
        Self {
            wells,
            include_b2m: false,
            analysis,
            melt_results,
            status_msg: String::from("MIQE RT-qPCR & SYBR Melt Studio ready — toggle unstable B2M reference to inspect geNorm M."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("qpcr_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("qpcr_gui_side").default_width(300.0).show(ctx, |ui| {
            ui.heading("Reference Genes & geNorm");
            if ui.checkbox(&mut self.include_b2m, "Include Unstable B2M Reference").changed() {
                let mut refs = vec![String::from("ACTB"), String::from("GAPDH")];
                if self.include_b2m {
                    refs.push(String::from("B2M"));
                }
                if let Ok(rep) = analyze_qpcr_experiment(&self.wells, &refs, "Control_37C") {
                    self.analysis = rep;
                    self.status_msg = format!("Refs updated: {}", refs.join(" + "));
                }
            }
            for g in &self.analysis.genorm_rankings {
                let col = if g.stable_pass { Color32::from_rgb(52, 211, 153) } else { Color32::from_rgb(248, 113, 113) };
                ui.colored_label(col, format!("{} : M = {:.3} ({})", g.gene_name, g.stability_m, if g.stable_pass { "STABLE" } else { "UNSTABLE" }));
            }
            ui.separator();
            ui.heading("Standard Curves");
            for sc in &self.analysis.standard_curves {
                ui.monospace(format!("{}: Eff={:.1}% R2={:.4} [{}]", sc.target_gene, sc.efficiency_pct, sc.r_squared, if sc.miqe_pass { "PASS" } else { "FLAG" }));
            }
            ui.separator();
            if ui.button("Export Expression & Melt .svg").clicked() {
                let _ = fs::write("qpcr_expression.svg", render_qpcr_expression_svg(&self.analysis, "OSABT RT-qPCR Analysis"));
                let _ = fs::write("qpcr_melt_curves.svg", render_melt_curves_svg(&self.melt_results, "OSABT SYBR Dissociation"));
                self.status_msg = String::from("Exported ./qpcr_expression.svg & ./qpcr_melt_curves.svg");
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let max_y = self.melt_results.iter().flat_map(|m| m.neg_df_dt.iter().copied()).fold(1.0_f64, f64::max) * 1.1;
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 240.0, "SYBR Green -dF/dT Dissociation Melt Curves");
            draw_grid_axes(&p, plot, 7, 4, (60.0, 95.0), (0.0, max_y), |v| format!("{v:.0}C"), |v| format!("{v:.0}"));
            for (mi, mr) in self.melt_results.iter().enumerate() {
                let col = PALETTE_ACCENTS[mi % PALETTE_ACCENTS.len()];
                let pts: Vec<Pos2> = mr.temperatures_c.iter().zip(&mr.neg_df_dt).map(|(&t, &df)| {
                    let fx = ((t - 60.0) / 35.0).clamp(0.0, 1.0) as f32;
                    let fy = (df / max_y).clamp(0.0, 1.0) as f32;
                    Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                }).collect();
                p.add(Shape::line(pts, Stroke::new(2.0_f32, col)));
            }

            ui.add_space(6.0);
            egui::Grid::new("qpcr_fc_grid").striped(true).show(ui, |ui| {
                for h in ["Target", "Comparison", "qBase FC (95% CI)", "log2(FC)", "Pfaffl", "Livak", "p-value"] {
                    ui.strong(h);
                }
                ui.end_row();
                for c in &self.analysis.group_comparisons {
                    ui.label(&c.target_gene);
                    ui.label(format!("{} vs {}", c.treatment_group, c.calibrator_group));
                    ui.strong(format!("{:.2}x ({:.2}–{:.2})", c.fold_change_qbase, c.fold_change_ci95.0, c.fold_change_ci95.1));
                    ui.label(format!("{:+.3}", c.log2_fold_change));
                    ui.label(format!("{:.2}x", c.pfaffl_ratio));
                    ui.label(format!("{:.2}x", c.livak_ratio));
                    ui.label(format!("{:.2e}", c.p_value_two_sided));
                    ui.end_row();
                }
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #14 eDNA Panel
// ---------------------------------------------------------------------------
pub struct EdnaPanelState {
    reads: Vec<FastqRead>,
    barcodes: Vec<DualBarcodeSample>,
    ref_db: Vec<(String, String)>,
    cfg: EdnaPipelineConfig,
    pipeline: EdnaPipelineReport,
    diversity: CommunityDiversityReport,
    status_msg: String,
}

impl EdnaPanelState {
    pub fn new() -> Self {
        let (fq, bc, db) = demo_pnw_12s_edna_dataset();
        let reads = parse_fastq_str(&fq).expect("demo fastq");
        let barcodes = parse_barcode_sheet_csv(&bc).expect("demo barcodes");
        let cfg = EdnaPipelineConfig::default();
        let pipeline = run_edna_metabarcoding_pipeline(&reads, &barcodes, &db, &cfg);
        let diversity = analyze_edna_community_diversity(&pipeline, 10);
        Self {
            reads,
            barcodes,
            ref_db: db,
            cfg,
            pipeline,
            diversity,
            status_msg: String::from("eDNA 12S Metabarcoding pipeline ready — adjust maxEE or chimera screening."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("edna_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("edna_gui_side").default_width(300.0).show(ctx, |ui| {
            ui.heading("Pipeline & Denoiser Controls");
            let mut chg = false;
            chg |= ui.add(egui::Slider::new(&mut self.cfg.max_expected_errors, 0.25..=6.0).text("Max Expected Errors")).changed();
            chg |= ui.checkbox(&mut self.cfg.remove_chimeras, "2-Parent Chimera Screen").changed();
            if chg {
                self.pipeline = run_edna_metabarcoding_pipeline(&self.reads, &self.barcodes, &self.ref_db, &self.cfg);
                self.diversity = analyze_edna_community_diversity(&self.pipeline, 10);
                self.status_msg = format!("Re-denoised: {} clean ASVs.", self.diversity.clean_asv_count);
            }
            ui.separator();
            ui.monospace(format!("Raw Reads   : {}", self.pipeline.total_raw_reads));
            ui.monospace(format!("Demux Reads : {}", self.pipeline.demux_assigned_reads));
            ui.monospace(format!("Passed maxEE: {}", self.pipeline.passed_ee_reads));
            ui.monospace(format!("Chimeras    : {}", self.pipeline.chimeric_reads_flagged));
            ui.monospace(format!("Clean ASVs  : {}", self.diversity.clean_asv_count));
            ui.separator();
            if ui.button("Export ASV .tsv + Rarefaction .svg").clicked() {
                let _ = fs::write("edna_asv_table.tsv", export_asv_table_tsv(&self.pipeline, false));
                let _ = fs::write("edna_rarefaction.svg", render_edna_rarefaction_svg(&self.diversity, "PNW 12S eDNA Biodiversity"));
                self.status_msg = String::from("Exported ./edna_asv_table.tsv & ./edna_rarefaction.svg");
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let max_m = self.diversity.rarefaction_curves.iter().flat_map(|c| c.depths_m.iter().copied()).max().unwrap_or(50) as f64;
            let max_s = self.diversity.rarefaction_curves.iter().flat_map(|c| c.expected_species_sm.iter().copied()).fold(5.0_f64, f64::max) * 1.15;
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 240.0, "Hurlbert (1971) Analytical Rarefaction E[S_m]");
            draw_grid_axes(&p, plot, 5, 4, (0.0, max_m), (0.0, max_s), |v| format!("{v:.0}"), |v| format!("{v:.1}"));
            for (ci, curve) in self.diversity.rarefaction_curves.iter().enumerate() {
                let col = PALETTE_ACCENTS[ci % PALETTE_ACCENTS.len()];
                let pts: Vec<Pos2> = curve.depths_m.iter().zip(&curve.expected_species_sm).map(|(&m, &s)| {
                    let fx = (m as f64 / max_m).clamp(0.0, 1.0) as f32;
                    let fy = (s / max_s).clamp(0.0, 1.0) as f32;
                    Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height())
                }).collect();
                p.add(Shape::line(pts, Stroke::new(2.2_f32, col)));
            }

            ui.add_space(6.0);
            egui::Grid::new("edna_alpha_grid").striped(true).show(ui, |ui| {
                for h in ["Site", "Reads", "^0D", "Chao1 (95% CI)", "^1D", "^2D", "Pielou J'", "Coverage"] {
                    ui.strong(h);
                }
                ui.end_row();
                for a in &self.diversity.alpha {
                    ui.label(&a.sample_name);
                    ui.label(format!("{}", a.total_reads_n));
                    ui.label(format!("{}", a.observed_richness_q0));
                    ui.label(format!("{:.1} ({:.1}–{:.1})", a.chao1_asymptotic, a.chao1_ci95.0, a.chao1_ci95.1));
                    ui.label(format!("{:.2}", a.hill_shannon_q1));
                    ui.label(format!("{:.2}", a.hill_simpson_q2));
                    ui.label(format!("{:.3}", a.pielou_evenness_j));
                    ui.label(format!("{:.1}%", a.sample_coverage * 100.0));
                    ui.end_row();
                }
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #15 Wildlife Panel
// ---------------------------------------------------------------------------
pub struct WildlifePanelState {
    occ_report: OccupancySuiteReport,
    closed_report: ClosedCaptureSuiteReport,
    cjs_report: CjsSurvivalReport,
    chapman_demo: ChapmanPetersenEstimate,
    target_se: f64,
    cost_ratio: f64,
    survey_plan: OccupancySurveyPlan,
    status_msg: String,
}

impl WildlifePanelState {
    pub fn new() -> Self {
        let occ_rows = parse_occupancy_csv(demo_pnw_occupancy_csv()).expect("demo occ");
        let occ_report = fit_mackenzie_occupancy_suite(&occ_rows).expect("occ fit");
        let (_, histories) = parse_capture_history_csv(demo_pnw_capture_history_csv()).expect("demo cap");
        let closed_report = analyze_closed_capture_histories(&histories).expect("closed");
        let cjs_report = estimate_cjs_apparent_survival(&histories).expect("cjs");
        let chapman_demo = calc_chapman_lincoln_petersen(40, 35, 14).expect("chapman");
        let best = &occ_report.models[occ_report.best_model_index];
        let target_se = 0.06;
        let cost_ratio = 2.0;
        let survey_plan = plan_occupancy_survey_effort(best.psi_hat.clamp(0.05, 0.95), best.p_hat.clamp(0.05, 0.95), 0.05, target_se, cost_ratio).expect("plan");
        Self {
            occ_report,
            closed_report,
            cjs_report,
            chapman_demo,
            target_se,
            cost_ratio,
            survey_plan,
            status_msg: String::from("MacKenzie Occupancy, Closed M0/Mt/Mh & CJS Survival Studio ready."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("wildlife_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("wildlife_gui_side").default_width(310.0).show(ctx, |ui| {
            ui.heading("Survey Effort Planner");
            let mut chg = false;
            chg |= ui.add(egui::Slider::new(&mut self.target_se, 0.02..=0.15).text("Target SE(psi)")).changed();
            chg |= ui.add(egui::Slider::new(&mut self.cost_ratio, 0.0..=12.0).text("Site/Visit Cost Ratio")).changed();
            if chg {
                let best = &self.occ_report.models[self.occ_report.best_model_index];
                if let Ok(pl) = plan_occupancy_survey_effort(best.psi_hat.clamp(0.05, 0.95), best.p_hat.clamp(0.05, 0.95), 0.05, self.target_se, self.cost_ratio) {
                    self.survey_plan = pl;
                }
            }
            ui.monospace(format!("Optimal Visits K* : {}", self.survey_plan.optimal_k_for_precision));
            ui.monospace(format!("Required Sites S* : {}", self.survey_plan.required_sites_for_target_se));
            ui.monospace(format!("Total Visits      : {}", self.survey_plan.total_survey_visits));
            ui.separator();
            ui.heading("Open CJS Survival");
            ui.monospace(format!("Phi (Survival) : {:.3} ± {:.3}", self.cjs_report.phi_hat, self.cjs_report.phi_se));
            ui.monospace(format!("p (Recapture)  : {:.3} ± {:.3}", self.cjs_report.p_hat, self.cjs_report.p_se));
            ui.monospace(format!("Chapman LP N   : {:.1} ± {:.1}", self.chapman_demo.chapman_n_hat, self.chapman_demo.seber_se));
            ui.separator();
            if ui.button("Export Wildlife .svg").clicked() {
                let svg = render_wildlife_population_svg(&self.occ_report, &self.closed_report, "PNW Occupancy & Mark-Recapture");
                if fs::write("wildlife_population.svg", svg).is_ok() {
                    self.status_msg = String::from("Exported ./wildlife_population.svg");
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let ests = [
                &self.closed_report.m0_equal_catchability,
                &self.closed_report.mt_darroch_time_varying,
                &self.closed_report.mh_chao1987_heterogeneity,
                &self.closed_report.schnabel_multicensus,
            ];
            let max_n = ests.iter().map(|e| e.ci95.1).fold(50.0_f64, f64::max) * 1.1;
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 220.0, "Closed Abundance Estimators N_hat (with 95% CI)");
            draw_grid_axes(&p, plot, 4, 4, (0.0, 4.0), (0.0, max_n), |_| String::new(), |v| format!("{v:.0}"));
            for (i, e) in ests.iter().enumerate() {
                let x = plot.left() + ((i as f32 + 0.5) / 4.0) * plot.width();
                let yn = plot.bottom() - (e.n_hat / max_n).clamp(0.0, 1.0) as f32 * plot.height();
                let ylo = plot.bottom() - (e.ci95.0 / max_n).clamp(0.0, 1.0) as f32 * plot.height();
                let yhi = plot.bottom() - (e.ci95.1 / max_n).clamp(0.0, 1.0) as f32 * plot.height();
                let col = PALETTE_ACCENTS[i % PALETTE_ACCENTS.len()];
                p.line_segment([Pos2::new(x, ylo), Pos2::new(x, yhi)], Stroke::new(2.2_f32, col));
                p.circle_filled(Pos2::new(x, yn), 6.0, col);
                p.text(Pos2::new(x, plot.bottom() + 4.0), egui::Align2::CENTER_TOP, e.model_name.split_whitespace().next().unwrap_or("M"), egui::FontId::monospace(10.0), col);
            }

            ui.add_space(6.0);
            egui::Grid::new("wildlife_occ_grid").striped(true).show(ui, |ui| {
                for h in ["Occupancy Model", "psi_hat (95% CI)", "p_hat (95% CI)", "AICc", "dAICc", "Weight"] {
                    ui.strong(h);
                }
                ui.end_row();
                for m in &self.occ_report.models {
                    ui.label(m.model.label());
                    ui.label(format!("{:.3} ({:.2}–{:.2})", m.psi_hat, m.psi_ci95.0, m.psi_ci95.1));
                    ui.label(format!("{:.3} ({:.2}–{:.2})", m.p_hat, m.p_ci95.0, m.p_ci95.1));
                    ui.label(format!("{:.2}", m.aicc));
                    ui.label(format!("{:.2}", m.delta_aicc));
                    ui.strong(format!("{:.1}%", m.akaike_weight * 100.0));
                    ui.end_row();
                }
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #17 CRISPR Panel
// ---------------------------------------------------------------------------
pub struct CrisprPanelState {
    target_raw: String,
    bg_raw: String,
    target_len_bp: usize,
    nuclease: CrisprNuclease,
    vector: CrisprCloningVector,
    guides: Vec<CrisprGuideCandidate>,
    offtargets: Vec<Vec<CrisprOffTargetHit>>,
    selected_idx: usize,
    status_msg: String,
}

impl CrisprPanelState {
    pub fn new() -> Self {
        let (target_raw, bg_raw) = demo_crispr_target_and_offtargets();
        let mut s = Self {
            target_raw,
            bg_raw,
            target_len_bp: 0,
            nuclease: CrisprNuclease::SpCas9,
            vector: CrisprCloningVector::Px459BbsI,
            guides: Vec::new(),
            offtargets: Vec::new(),
            selected_idx: 0,
            status_msg: String::new(),
        };
        s.rescan();
        s
    }

    fn rescan(&mut self) {
        if let Ok((norm_seq, mut guides)) = scan_crispr_candidates(&self.target_raw, self.nuclease) {
            if let Ok(mut hits) = scan_offtargets_and_specificity(&mut guides, &self.bg_raw, 4) {
                let mut paired: Vec<(CrisprGuideCandidate, Vec<CrisprOffTargetHit>)> = guides.into_iter().zip(hits.drain(..)).collect();
                paired.sort_by(|a, b| b.0.composite_score.partial_cmp(&a.0.composite_score).unwrap_or(std::cmp::Ordering::Equal));
                for (i, (g, _)) in paired.iter_mut().enumerate() {
                    g.guide_id = format!("gRNA_{:02}", i + 1);
                }
                let (gs, hs): (Vec<_>, Vec<_>) = paired.into_iter().unzip();
                self.target_len_bp = norm_seq.len();
                self.guides = gs;
                self.offtargets = hs;
                self.selected_idx = 0;
                self.status_msg = format!("Scored {} {} guides on {} bp target.", self.guides.len(), self.nuclease.label(), self.target_len_bp);
            }
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("crispr_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("crispr_gui_side").default_width(310.0).show(ctx, |ui| {
            ui.heading("CRISPR Nuclease & Vector");
            let mut chg = false;
            for n in [CrisprNuclease::SpCas9, CrisprNuclease::SaCas9, CrisprNuclease::LbCas12a] {
                if ui.selectable_value(&mut self.nuclease, n, n.label()).changed() {
                    self.vector = match n {
                        CrisprNuclease::LbCas12a => CrisprCloningVector::LbCas12aBsaI,
                        _ => CrisprCloningVector::Px459BbsI,
                    };
                    chg = true;
                }
            }
            ui.separator();
            ui.label("Golden Gate Vector:");
            for v in [CrisprCloningVector::Px459BbsI, CrisprCloningVector::LentiCrisprV2BsmBI, CrisprCloningVector::LbCas12aBsaI] {
                ui.selectable_value(&mut self.vector, v, v.label());
            }
            if chg {
                self.rescan();
            }
            if let Some(g) = self.guides.get(self.selected_idx) {
                ui.separator();
                ui.strong(format!("Selected: {} ({} @ {} bp)", g.guide_id, g.strand, g.cut_site_1idx));
                let ol = design_golden_gate_crispr_oligos(&g.spacer_seq, self.vector);
                ui.monospace(format!("Top: {}", ol.top_oligo_5to3));
                ui.monospace(format!("Bot: {}", ol.bottom_oligo_5to3));
                let be = audit_base_editor_window(&g.spacer_seq);
                ui.small(format!("CBE C@4..8: {:?} | ABE A@4..8: {:?}", be.cbe_c_positions, be.abe_a_positions));
            }
            ui.separator();
            if ui.button("Export Guides .tsv + Locus .svg").clicked() {
                let _ = fs::write("crispr_guides.tsv", export_crispr_design_tsv(&self.guides, self.vector));
                let _ = fs::write("crispr_locus.svg", render_crispr_locus_svg(self.target_len_bp, &self.guides, "OSABT CRISPR Locus Guide Map"));
                self.status_msg = String::from("Exported ./crispr_guides.tsv & ./crispr_locus.svg");
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            let len = self.target_len_bp.max(1) as f64;
            let (plot, p) = draw_chart_frame(ui, ui.available_width(), 200.0, "Locus Cut Sites vs Composite Score (0..100)");
            draw_grid_axes(&p, plot, 5, 4, (1.0, len), (0.0, 100.0), |v| format!("{v:.0}bp"), |v| format!("{v:.0}"));
            for (i, g) in self.guides.iter().enumerate() {
                let fx = (g.cut_site_1idx as f64 / len).clamp(0.0, 1.0) as f32;
                let fy = (g.composite_score / 100.0).clamp(0.0, 1.0) as f32;
                let pt = Pos2::new(plot.left() + fx * plot.width(), plot.bottom() - fy * plot.height());
                let col = if i == self.selected_idx { Color32::from_rgb(250, 204, 21) } else { Color32::from_rgb(56, 189, 248) };
                p.line_segment([Pos2::new(pt.x, plot.bottom()), pt], Stroke::new(1.8_f32, col));
                p.circle_filled(pt, if i == self.selected_idx { 6.0 } else { 4.2 }, col);
            }

            ui.add_space(6.0);
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("crispr_guides_grid").striped(true).show(ui, |ui| {
                    for h in ["ID", "Cut", "Str", "Spacer (5'->3')", "PAM", "GC%", "On-Target", "CFD Spec", "Composite"] {
                        ui.strong(h);
                    }
                    ui.end_row();
                    for (i, g) in self.guides.iter().enumerate() {
                        if ui.selectable_label(self.selected_idx == i, &g.guide_id).clicked() {
                            self.selected_idx = i;
                        }
                        ui.label(format!("{}", g.cut_site_1idx));
                        ui.label(format!("{}", g.strand));
                        ui.monospace(&g.spacer_seq);
                        ui.monospace(&g.pam_seq);
                        ui.label(format!("{:.0}%", g.gc_percent));
                        ui.label(format!("{:.1}", g.on_target_score));
                        ui.label(format!("{:.1}", g.cfd_specificity_score));
                        ui.strong(format!("{:.1}", g.composite_score));
                        ui.end_row();
                    }
                });
            });
        });
    }
}

// ---------------------------------------------------------------------------
// #19 Survival Panel
// ---------------------------------------------------------------------------
pub struct SurvivalPanelState {
    curves: Vec<KaplanMeierCurve>,
    log_rank: Option<LogRankTestResult>,
    cox_cov: CoxPhFit,
    cox_grp: CoxPhFit,
    ecotox_pts: Vec<EcotoxBioassayPoint>,
    link: DoseResponseLink,
    lc50_report: EcotoxLc50Report,
    status_msg: String,
}

impl SurvivalPanelState {
    pub fn new() -> Self {
        let obs = parse_survival_csv(demo_survival_cohort_csv()).expect("demo survival");
        let curves = compute_kaplan_meier_curves(&obs);
        let log_rank = compute_log_rank_test(&obs).ok();
        let cox_cov = fit_cox_proportional_hazards(&obs, "Covariate", false).expect("cox cov");
        let cox_grp = fit_cox_proportional_hazards(&obs, "Group", true).expect("cox grp");
        let ecotox_pts = parse_ecotox_csv(demo_ecotox_lc50_csv()).expect("demo tox");
        let link = DoseResponseLink::Probit;
        let lc50_report = fit_ecotox_lc50(&ecotox_pts, link).expect("lc50");
        Self {
            curves,
            log_rank,
            cox_cov,
            cox_grp,
            ecotox_pts,
            link,
            lc50_report,
            status_msg: String::from("Kaplan-Meier, Mantel Log-Rank, Cox PH & Finney Probit/Logit LC50 Studio ready."),
        }
    }

    pub fn draw(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::bottom("survival_gui_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });
        egui::SidePanel::left("survival_gui_side").default_width(305.0).show(ctx, |ui| {
            ui.heading("Survival & Ecotox Models");
            let mut link_changed = false;
            link_changed |= ui.selectable_value(&mut self.link, DoseResponseLink::Probit, "Finney Probit LC50").changed();
            link_changed |= ui.selectable_value(&mut self.link, DoseResponseLink::Logit, "Binomial Logit LC50").changed();
            if link_changed {
                if let Ok(rep) = fit_ecotox_lc50(&self.ecotox_pts, self.link) {
                    self.lc50_report = rep;
                    self.status_msg = format!("{} LC50 = {:.3}", self.link.label(), self.lc50_report.lc50);
                }
            }
            ui.separator();
            ui.heading("Mantel Log-Rank & Cox PH");
            if let Some(lr) = &self.log_rank {
                ui.monospace(format!("Log-Rank chi^2 = {:.3} (p = {:.5})", lr.chi_square, lr.p_value));
            }
            ui.monospace(format!("Group HR = {:.2} (95% CI [{:.2}, {:.2}])", self.cox_grp.hazard_ratio, self.cox_grp.hr_ci95.0, self.cox_grp.hr_ci95.1));
            ui.monospace(format!("Covar HR = {:.2} (p = {:.4})", self.cox_cov.hazard_ratio, self.cox_cov.p_value));
            ui.separator();
            ui.heading("Ecotoxicology Thresholds");
            ui.monospace(format!("LC10 = {:.3}", self.lc50_report.lc10));
            ui.strong(format!("LC50 = {:.3} [{:.2}, {:.2}]", self.lc50_report.lc50, self.lc50_report.lc50_ci95.0, self.lc50_report.lc50_ci95.1));
            ui.monospace(format!("LC90 = {:.3} | chi^2 = {:.2}", self.lc50_report.lc90, self.lc50_report.pearson_chi2));
            ui.separator();
            if ui.button("Export Dual-Panel .svg").clicked() {
                let svg = render_survival_and_lc50_svg(&self.curves, self.log_rank.as_ref(), &self.ecotox_pts, &self.lc50_report, "OSABT Survival & Ecotoxicology");
                if fs::write("survival_lc50.svg", svg).is_ok() {
                    self.status_msg = String::from("Exported ./survival_lc50.svg");
                }
            }
        });

        egui::CentralPanel::default().show(ctx, |ui| {
            ui.columns(2, |cols| {
                let max_t = self.curves.iter().flat_map(|c| c.steps.iter().map(|s| s.time)).fold(10.0_f64, f64::max) * 1.05;
                let w0 = cols[0].available_width();
                let (p_km, km) = draw_chart_frame(&mut cols[0], w0, 320.0, "Kaplan-Meier Survival S(t) (95% Greenwood CI)");
                draw_grid_axes(&km, p_km, 5, 5, (0.0, max_t), (0.0, 1.0), |v| format!("{v:.0}d"), |v| format!("{v:.1}"));
                let map_km = |t: f64, s: f64| Pos2::new(p_km.left() + (t / max_t).clamp(0.0, 1.0) as f32 * p_km.width(), p_km.bottom() - s.clamp(0.0, 1.0) as f32 * p_km.height());
                for (ci, curve) in self.curves.iter().enumerate() {
                    let col = PALETTE_ACCENTS[ci % PALETTE_ACCENTS.len()];
                    let band = Color32::from_rgba_unmultiplied(col.r(), col.g(), col.b(), 30);
                    let mut pts = vec![map_km(0.0, 1.0)];
                    let mut prev_t = 0.0;
                    let mut prev_s = 1.0;
                    let mut prev_ci = (1.0, 1.0);
                    for st in &curve.steps {
                        let tl = map_km(prev_t, prev_ci.1);
                        let br = map_km(st.time, prev_ci.0);
                        km.rect_filled(Rect::from_min_max(Pos2::new(tl.x.min(br.x), tl.y.min(br.y)), Pos2::new(tl.x.max(br.x), tl.y.max(br.y))), 0.0, band);
                        pts.push(map_km(st.time, prev_s));
                        pts.push(map_km(st.time, st.survival_s));
                        prev_t = st.time;
                        prev_s = st.survival_s;
                        prev_ci = st.ci95;
                    }
                    km.add(Shape::line(pts, Stroke::new(2.4_f32, col)));
                }

                let min_lx = self.ecotox_pts.first().map(|x| x.concentration.log10()).unwrap_or(-1.0) - 0.3;
                let max_lx = self.ecotox_pts.last().map(|x| x.concentration.log10()).unwrap_or(2.0) + 0.3;
                let w1 = cols[1].available_width();
                let (p_tx, tx) = draw_chart_frame(&mut cols[1], w1, 320.0, &format!("{} Dose-Mortality Curve", self.link.label()));
                draw_grid_axes(&tx, p_tx, 4, 5, (min_lx, max_lx), (0.0, 1.0), |v| format!("10^{v:.1}"), |v| format!("{:.0}%", v * 100.0));
                let span_lx = (max_lx - min_lx).max(1e-6);
                let curve_pts: Vec<Pos2> = (0..=50).map(|k| {
                    let lx = min_lx + span_lx * (k as f64 / 50.0);
                    let eta = self.lc50_report.intercept_alpha + self.lc50_report.slope_beta * lx;
                    let pr = match self.link {
                        DoseResponseLink::Probit => crate::epi_models::norm_cdf(eta),
                        DoseResponseLink::Logit => 1.0 / (1.0 + (-eta).exp()),
                    };
                    Pos2::new(p_tx.left() + ((lx - min_lx) / span_lx) as f32 * p_tx.width(), p_tx.bottom() - pr.clamp(0.0, 1.0) as f32 * p_tx.height())
                }).collect();
                tx.add(Shape::line(curve_pts, Stroke::new(2.4_f32, Color32::from_rgb(248, 113, 113))));
                for pt in &self.ecotox_pts {
                    let lx = pt.concentration.log10();
                    let mort = pt.dead_n as f64 / pt.exposed_n.max(1) as f64;
                    let pos = Pos2::new(p_tx.left() + ((lx - min_lx) / span_lx) as f32 * p_tx.width(), p_tx.bottom() - mort.clamp(0.0, 1.0) as f32 * p_tx.height());
                    tx.circle_filled(pos, 5.0, Color32::from_rgb(56, 189, 248));
                }
            });

            ui.add_space(6.0);
            egui::Grid::new("survival_ecotox_grid").striped(true).show(ui, |ui| {
                for h in ["Cohort Stratum", "N", "Events", "Median Survival", "Final S(t)"] {
                    ui.strong(h);
                }
                ui.end_row();
                for c in &self.curves {
                    ui.label(&c.group);
                    ui.label(format!("{}", c.total_n));
                    ui.label(format!("{}", c.total_events));
                    ui.label(c.median_survival_time.map(|m| format!("{m:.1} d")).unwrap_or_else(|| String::from(">max")));
                    ui.label(format!("{:.3}", c.steps.last().map(|s| s.survival_s).unwrap_or(1.0)));
                    ui.end_row();
                }
            });
        });
    }
}
