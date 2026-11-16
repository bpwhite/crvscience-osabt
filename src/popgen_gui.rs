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

//! Interactive `eframe`/`egui` desktop/WASM Population Genetics, Weir–Cockerham `F_ST`,
//! & Bandelt Median-Joining Haplotype Network Studio (`osabt popgen`).

use std::fs;
use std::path::PathBuf;

use eframe::egui;

use crate::popgen_core::{
    analyze_diploid_vcf, analyze_popgen_fasta, demo_pnw_popgen_fasta, demo_pnw_popgen_vcf,
    parse_popgen_fasta, parse_popgen_vcf, DiploidVcfReport, PopGenFastaReport,
};
use crate::popgen_network::{
    build_median_joining_network, render_popgen_network_svg, HaplotypeNetwork,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PopGenView {
    HaplotypeNetworkCanvas,
    SequenceDiversityAndFst,
    DiploidVcfStructureAndHwe,
}

/// Launches the `eframe`/`egui` native desktop studio for `osabt popgen`.
pub fn run_popgen_gui(
    fasta_path: Option<PathBuf>,
    vcf_path: Option<PathBuf>,
) -> Result<(), String> {
    let app = PopGenGuiApp::new(fasta_path, vcf_path)?;
    let native_options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_title("OSABT PopGen Studio (v1.6.0) — Median-Joining Networks, Tajima's D & Weir–Cockerham F_ST")
            .with_inner_size([1160.0, 740.0]),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-popgen",
        native_options,
        Box::new(|_cc| Ok(Box::new(app))),
    )
    .map_err(|e| format!("Failed launching egui PopGen Studio: {e}"))
}

pub struct PopGenGuiApp {
    view: PopGenView,
    fasta_text: String,
    vcf_text: String,
    epsilon: usize,
    selected_node: Option<usize>,
    fasta_report: Option<PopGenFastaReport>,
    network: Option<HaplotypeNetwork>,
    vcf_report: Option<DiploidVcfReport>,
    status_msg: String,
}

impl PopGenGuiApp {
    pub fn new(fasta_path: Option<PathBuf>, vcf_path: Option<PathBuf>) -> Result<Self, String> {
        let fasta_text = if let Some(fp) = fasta_path {
            fs::read_to_string(&fp)
                .map_err(|e| format!("Failed reading FASTA '{}': {e}", fp.display()))?
        } else {
            demo_pnw_popgen_fasta()
        };
        let vcf_text = if let Some(vp) = vcf_path {
            fs::read_to_string(&vp)
                .map_err(|e| format!("Failed reading VCF '{}': {e}", vp.display()))?
        } else {
            demo_pnw_popgen_vcf().to_string()
        };

        let mut app = Self {
            view: PopGenView::HaplotypeNetworkCanvas,
            fasta_text,
            vcf_text,
            epsilon: 0,
            selected_node: None,
            fasta_report: None,
            network: None,
            vcf_report: None,
            status_msg: String::new(),
        };
        app.recompute_all();
        Ok(app)
    }

    fn recompute_all(&mut self) {
        match parse_popgen_fasta(&self.fasta_text) {
            Ok(recs) => {
                self.fasta_report = analyze_popgen_fasta(&recs).ok();
                self.network = build_median_joining_network(&recs, self.epsilon).ok();
            }
            Err(e) => {
                self.status_msg = format!("FASTA error: {e}");
                return;
            }
        }
        match parse_popgen_vcf(&self.vcf_text) {
            Ok((pops, loci)) => {
                self.vcf_report = analyze_diploid_vcf(&pops, &loci).ok();
            }
            Err(e) => {
                self.status_msg = format!("VCF error: {e}");
                return;
            }
        }
        if let Some(net) = &self.network {
            self.status_msg = format!(
                "Computed MJN (epsilon={}): {} sampled haplotypes + {} inferred Median Vectors ({} total mutations).",
                self.epsilon,
                net.sampled_haplotypes,
                net.inferred_median_vectors,
                net.total_network_steps
            );
        }
    }
}

impl eframe::App for PopGenGuiApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.draw_panels(ctx);
    }
}

impl PopGenGuiApp {
    /// Renders all PopGen Studio panels onto `ctx`.
    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("popgen_top_bar").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.heading("OSABT PopGen v1.6.0");
                ui.separator();
                ui.selectable_value(
                    &mut self.view,
                    PopGenView::HaplotypeNetworkCanvas,
                    "1. Median-Joining Haplotype Network",
                );
                ui.selectable_value(
                    &mut self.view,
                    PopGenView::SequenceDiversityAndFst,
                    "2. FASTA Diversity (pi, theta_W, Tajima's D & d_XY)",
                );
                ui.selectable_value(
                    &mut self.view,
                    PopGenView::DiploidVcfStructureAndHwe,
                    "3. Diploid VCF (Weir–Cockerham F_ST & Exact HWE)",
                );
                ui.separator();
                if ui
                    .add(egui::Slider::new(&mut self.epsilon, 0..=4).text("MJN epsilon"))
                    .changed()
                {
                    self.recompute_all();
                }
                if ui.button("Export Network .svg").clicked() {
                    if let Some(net) = &self.network {
                        let svg = render_popgen_network_svg(net, "PNW Phylogeography");
                        if fs::write("popgen_haplotype_network.svg", svg).is_ok() {
                            self.status_msg =
                                String::from("Saved popgen_haplotype_network.svg!");
                        }
                    }
                }
            });
        });

        egui::TopBottomPanel::bottom("popgen_status").show(ctx, |ui| {
            ui.label(&self.status_msg);
        });

        egui::CentralPanel::default().show(ctx, |ui| match self.view {
            PopGenView::HaplotypeNetworkCanvas => self.draw_network_canvas(ui),
            PopGenView::SequenceDiversityAndFst => self.draw_diversity_tables(ui),
            PopGenView::DiploidVcfStructureAndHwe => self.draw_vcf_tables(ui),
        });
    }
}

impl PopGenGuiApp {
    fn draw_network_canvas(&mut self, ui: &mut egui::Ui) {
        const COLORS: [egui::Color32; 6] = [
            egui::Color32::from_rgb(37, 99, 235),
            egui::Color32::from_rgb(22, 163, 74),
            egui::Color32::from_rgb(220, 38, 38),
            egui::Color32::from_rgb(147, 51, 234),
            egui::Color32::from_rgb(217, 119, 6),
            egui::Color32::from_rgb(8, 145, 178),
        ];

        let Some(net) = &self.network else {
            ui.label("No haplotype network loaded.");
            return;
        };

        ui.columns(2, |cols| {
            let (resp, painter) = cols[0].allocate_painter(
                egui::vec2(cols[0].available_width().max(520.0), 580.0),
                egui::Sense::click(),
            );
            let rect = resp.rect;
            painter.rect_filled(rect, 6.0, egui::Color32::from_rgb(248, 250, 252));

            let to_screen = |(nx, ny): (f64, f64)| -> egui::Pos2 {
                egui::pos2(
                    rect.left() + (nx as f32) * rect.width(),
                    rect.top() + (ny as f32) * rect.height(),
                )
            };

            for edge in &net.edges {
                let p1 = to_screen(net.nodes[edge.from_idx].pos_xy);
                let p2 = to_screen(net.nodes[edge.to_idx].pos_xy);
                painter.line_segment(
                    [p1, p2],
                    egui::Stroke::new(2.0_f32, egui::Color32::from_rgb(71, 85, 105)),
                );
                let dir = (p2 - p1).normalized();
                let perp = egui::vec2(-dir.y, dir.x);
                let ticks = edge.mutations.min(10);
                for t in 1..=ticks {
                    let frac = 0.35 + 0.30 * ((t as f32) / ((ticks + 1) as f32));
                    let mid = p1 + (p2 - p1) * frac;
                    painter.line_segment(
                        [mid - perp * 5.0, mid + perp * 5.0],
                        egui::Stroke::new(1.8_f32, egui::Color32::from_rgb(220, 38, 38)),
                    );
                }
            }

            let click_pos = if resp.clicked() {
                resp.interact_pointer_pos()
            } else {
                None
            };

            for (idx, node) in net.nodes.iter().enumerate() {
                let center = to_screen(node.pos_xy);
                let rad = if node.is_median_vector {
                    6.5_f32
                } else {
                    9.0 + 4.2 * (node.total_samples.max(1) as f32).sqrt()
                };
                if let Some(cp) = click_pos {
                    if cp.distance(center) <= rad + 4.0 {
                        self.selected_node = Some(idx);
                    }
                }

                if node.is_median_vector {
                    painter.circle_filled(center, rad, egui::Color32::from_rgb(15, 23, 42));
                } else {
                    let dom_pop = node
                        .counts_by_pop
                        .iter()
                        .enumerate()
                        .max_by_key(|&(_, c)| c)
                        .map(|(i, _)| i)
                        .unwrap_or(0);
                    painter.circle_filled(center, rad, COLORS[dom_pop % COLORS.len()]);
                    painter.circle_stroke(
                        center,
                        rad,
                        egui::Stroke::new(1.6_f32, egui::Color32::from_rgb(15, 23, 42)),
                    );
                }
                painter.text(
                    center + egui::vec2(0.0, rad + 10.0),
                    egui::Align2::CENTER_CENTER,
                    &node.node_id,
                    egui::FontId::proportional(12.0),
                    egui::Color32::from_rgb(15, 23, 42),
                );
            }

            cols[1].heading("Network Legend & Node Inspector");
            cols[1].separator();
            for (pi, pop) in net.populations.iter().enumerate() {
                cols[1].colored_label(COLORS[pi % COLORS.len()], format!("■ {pop}"));
            }
            cols[1].label("◆ Black Circle = Bandelt Inferred Median Vector (MV)");
            cols[1].separator();
            if let Some(sel) = self.selected_node.and_then(|i| net.nodes.get(i)) {
                cols[1].strong(format!(
                    "Selected Vertex: {} ({} samples)",
                    sel.node_id, sel.total_samples
                ));
                for (pop, &cnt) in net.populations.iter().zip(&sel.counts_by_pop) {
                    if cnt > 0 {
                        cols[1].label(format!("  • {pop}: {cnt} individuals"));
                    }
                }
                cols[1].monospace(&sel.sequence);
            } else {
                cols[1].label("Click any haplotype circle or Median Vector on the canvas to inspect its watershed composition and sequence.");
            }
        });
    }

    fn draw_diversity_tables(&self, ui: &mut egui::Ui) {
        let Some(rep) = &self.fasta_report else {
            ui.label("No FASTA diversity report available.");
            return;
        };
        ui.heading("Within-Population Polymorphism & Tajima's (1989) Neutrality Test");
        egui::Grid::new("fasta_div_grid").striped(true).show(ui, |ui| {
            for h in [
                "Population",
                "n",
                "Haplotypes",
                "H_d",
                "Seg Sites S",
                "pi (/site)",
                "theta_W (/site)",
                "Tajima's D",
                "p-value",
            ] {
                ui.strong(h);
            }
            ui.end_row();

            for m in std::iter::once(&rep.global_metrics).chain(rep.per_population.iter()) {
                ui.label(&m.population);
                ui.label(format!("{}", m.sample_size_n));
                ui.label(format!("{}", m.num_haplotypes));
                ui.label(format!("{:.3}", m.haplotype_diversity_hd));
                ui.label(format!("{}", m.segregating_sites_s));
                ui.label(format!("{:.4}", m.nucleotide_diversity_pi));
                ui.label(format!("{:.4}", m.watterson_theta_w));
                ui.label(format!("{:+.3}", m.tajima_d));
                ui.label(format!("{:.4}", m.tajima_d_pvalue));
                ui.end_row();
            }
        });

        ui.add_space(16.0);
        ui.heading("Pairwise Population Sequence Divergence (d_XY, net d_a & Hudson F_ST)");
        egui::Grid::new("fasta_fst_grid").striped(true).show(ui, |ui| {
            for h in ["Population A", "Population B", "d_XY", "Net d_a", "Hudson F_ST"] {
                ui.strong(h);
            }
            ui.end_row();
            for p in &rep.pairwise_divergence {
                ui.label(&p.pop_a);
                ui.label(&p.pop_b);
                ui.label(format!("{:.4}", p.d_xy));
                ui.label(format!("{:.4}", p.net_d_a));
                ui.label(format!("{:.4}", p.hudson_fst));
                ui.end_row();
            }
        });
    }

    fn draw_vcf_tables(&self, ui: &mut egui::Ui) {
        let Some(v) = &self.vcf_report else {
            ui.label("No diploid VCF loaded.");
            return;
        };
        ui.heading(format!(
            "Weir & Cockerham (1984) Multilocus Summary ({} Diploid Samples, {} SNPs): F_ST = {:.4} | F_IS = {:+.4} | F_IT = {:+.4}",
            v.num_samples, v.num_loci, v.multilocus_fst, v.multilocus_fis, v.multilocus_fit
        ));
        ui.separator();
        egui::Grid::new("vcf_loci_grid").striped(true).show(ui, |ui| {
            for h in [
                "SNP ID",
                "Chrom:Pos",
                "Genotypes (0/0 : 0/1 : 1/1)",
                "MAF",
                "H_obs",
                "H_exp",
                "Wigginton Exact HWE p",
                "Weir–Cockerham F_ST",
            ] {
                ui.strong(h);
            }
            ui.end_row();
            for l in &v.loci {
                ui.label(&l.snp_id);
                ui.label(format!("{}:{}", l.chrom, l.pos));
                ui.label(format!(
                    "{} : {} : {}",
                    l.obs_hom_ref, l.obs_het, l.obs_hom_alt
                ));
                ui.label(format!("{:.3}", l.minor_allele_freq));
                ui.label(format!("{:.3}", l.obs_heterozygosity));
                ui.label(format!("{:.3}", l.exp_heterozygosity));
                ui.label(format!("{:.4}", l.hwe_exact_pvalue));
                ui.label(format!("{:.4}", l.weir_cockerham_fst));
                ui.end_row();
            }
        });
    }
}
