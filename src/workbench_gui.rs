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

//! Unified 20-Tool Desktop Workbench (`eframe` / `egui`) for Open Source Ad-hoc Biology Tools (`osabt`).

use std::path::Path;

use eframe::egui;
use egui::{Color32, Margin, RichText, Stroke};

use crate::ecosound_gui::EcosoundStudioApp;
use crate::flow_gui::FlowGuiApp;
use crate::flux_gui::{render_flux_offscreen_png, FluxGuiApp};
use crate::forge3d_gui::Forge3dStudioApp;
use crate::kinetics_gui::KineticsStudioApp;
use crate::morpho_gui::MorphoStudioApp;
use crate::plasmid_gui::PlasmidStudioApp;
use crate::popgen_gui::PopGenGuiApp;

#[path = "workbench_panels.rs"]
mod workbench_panels;
use workbench_panels::{
    BenchPanelState, CrisprPanelState, DwcPanelState, EdnaPanelState, EpiPanelState,
    InfoEvoPanelState, LexiPanelState, PhyloPanelState, PlatePanelState, QpcrPanelState,
    SurvivalPanelState, WildlifePanelState,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToolDomain {
    MolecularCloning,
    EvolutionPopGen,
    AssaysKinetics,
    EcologyField,
    EpiBiostats,
}

impl ToolDomain {
    pub const ALL: [Self; 5] = [
        Self::MolecularCloning,
        Self::EvolutionPopGen,
        Self::AssaysKinetics,
        Self::EcologyField,
        Self::EpiBiostats,
    ];

    pub fn label(self) -> &'static str {
        match self {
            Self::MolecularCloning => "Molecular & Cloning",
            Self::EvolutionPopGen => "Evolution & PopGen",
            Self::AssaysKinetics => "Assays & Kinetics",
            Self::EcologyField => "Ecology & Field",
            Self::EpiBiostats => "Epi & Biostatistics",
        }
    }

    pub fn accent(self) -> Color32 {
        match self {
            Self::MolecularCloning => Color32::from_rgb(56, 189, 248),
            Self::EvolutionPopGen => Color32::from_rgb(192, 132, 252),
            Self::AssaysKinetics => Color32::from_rgb(250, 204, 21),
            Self::EcologyField => Color32::from_rgb(52, 211, 153),
            Self::EpiBiostats => Color32::from_rgb(248, 113, 113),
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ToolId {
    Bench,
    Plasmid,
    Dwc,
    Phylo,
    Epi,
    Lexi,
    Forge3d,
    Morpho,
    Plate,
    Ecosound,
    Infoevo,
    Qpcr,
    Kinetics,
    Edna,
    Wildlife,
    Popgen,
    Crispr,
    Flow,
    Survival,
    Flux,
}

struct ToolMeta {
    id: ToolId,
    num: u8,
    short_name: &'static str,
    title: &'static str,
    domain: ToolDomain,
    summary: &'static str,
    exports: &'static str,
    is_dedicated_studio: bool,
}

const TOOL_CATALOG: [ToolMeta; 20] = [
    ToolMeta {
        id: ToolId::Bench,
        num: 1,
        short_name: "Bench",
        title: "Primer, Master Mix & Buffer Bench",
        domain: ToolDomain::MolecularCloning,
        summary: "SantaLucia/Owczarzy primer pairs, 3'-dimer QC, PCR master mix & Henderson-Hasselbalch buffers.",
        exports: ".md checklist",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Plasmid,
        num: 2,
        short_name: "Plasmid",
        title: "Plasmid Map, Digest & Cloning Studio",
        domain: ToolDomain::MolecularCloning,
        summary: "Circular vector canvas, restriction digest, virtual agarose gel, Gibson/Golden Gate & .ab1 traces.",
        exports: ".gbk / .svg",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Dwc,
        num: 3,
        short_name: "DwC",
        title: "Darwin Core Linter & DwC-A Packager",
        domain: ToolDomain::EcologyField,
        summary: "CoordinateCleaner checks, TDWG term verification, taxonomic synonym resolution & occurrence map.",
        exports: "DwC-A / .geojson / .svg",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Phylo,
        num: 4,
        short_name: "Phylo",
        title: "Phylogenetic Tree Explorer",
        domain: ToolDomain::EvolutionPopGen,
        summary: "Newick/NEXUS tree visualization, midpoint & outgroup rerooting, ladderizing & clade search.",
        exports: ".svg / .nwk",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Epi,
        num: 5,
        short_name: "Epi",
        title: "Epidemiology & Study Power",
        domain: ToolDomain::EpiBiostats,
        summary: "Cori renewal R_t, RK4 SIR/SEIR dynamics, Bayesian diagnostic cascades & Fleiss study power.",
        exports: ".svg curves",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Lexi,
        num: 6,
        short_name: "Lexi",
        title: "Taxonomic & Anatomical Etymology",
        domain: ToolDomain::EpiBiostats,
        summary: "Dynamic-programming Greek/Latin morpheme decoder, cognate trees & interactive Anki flashcard quiz.",
        exports: "Anki .tsv / .md",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Forge3d,
        num: 7,
        short_name: "Forge3D",
        title: "Biomolecule & Morphology 3D Print Studio",
        domain: ToolDomain::MolecularCloning,
        summary: "PDB solvent surfaces, Gielis diatoms, radiolarians & voxel lattices with manifold STL/3MF audit.",
        exports: ".stl / .3mf / .obj",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Morpho,
        num: 8,
        short_name: "Morpho",
        title: "CFU Counter, Herbivory & 2D GPA",
        domain: ToolDomain::EcologyField,
        summary: "Automated Petri colony counting, Excess-Green leaf herbivory & Procrustes landmark morphometrics.",
        exports: ".tps / .svg",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Plate,
        num: 9,
        short_name: "Plate",
        title: "96/384-Well Assay & 4PL IC50 Studio",
        domain: ToolDomain::AssaysKinetics,
        summary: "Edge-buffered randomized microplate design, Zhang Z'-factor QC & nonlinear 4PL dose-response.",
        exports: ".csv / .md / .svg",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Ecosound,
        num: 10,
        short_name: "EcoSound",
        title: "Bioacoustics Spectrogram & Call Slicer",
        domain: ToolDomain::EcologyField,
        summary: "Interactive STFT spectrograms, ACI/NDSI/H/ADI biodiversity indices & Raven annotation tables.",
        exports: ".wav clips / Raven .txt / .svg",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Infoevo,
        num: 11,
        short_name: "InfoEvo",
        title: "Sequence Logos, Epistasis & Kimura",
        domain: ToolDomain::EvolutionPopGen,
        summary: "Schneider sequence logos, Adami physical complexity, APC-corrected MIp epistasis & selection regret.",
        exports: "3 x .svg",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Qpcr,
        num: 12,
        short_name: "qPCR",
        title: "MIQE RT-qPCR & Melt-Curve Studio",
        domain: ToolDomain::AssaysKinetics,
        summary: "qBase/Pfaffl/2^-ddCq fold change, geNorm M reference stability, SYBR -dF/dT peaks & MIQE audit.",
        exports: ".svg plots",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Kinetics,
        num: 13,
        short_name: "Kinetics",
        title: "Enzyme Inhibition & Morrison Studio",
        domain: ToolDomain::AssaysKinetics,
        summary: "Global reversible Ki model selection (AICc), Eadie-Hofstee plots & Morrison tight-binding Ki.",
        exports: ".svg plots",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Edna,
        num: 14,
        short_name: "eDNA",
        title: "eDNA Metabarcoding & Hill Diversity",
        domain: ToolDomain::EcologyField,
        summary: "Dual-barcode FASTQ demux, maxEE ASV denoising, chimera screen & Chao1/Hurlbert rarefaction.",
        exports: ".tsv / .svg",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Wildlife,
        num: 15,
        short_name: "Wildlife",
        title: "Occupancy & Mark-Recapture Studio",
        domain: ToolDomain::EcologyField,
        summary: "MacKenzie occupancy AICc, closed M0/Mt/Mh abundance, open CJS survival & effort planner.",
        exports: ".svg plots",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Popgen,
        num: 16,
        short_name: "PopGen",
        title: "Tajima's D, F_ST & Haplotype Networks",
        domain: ToolDomain::EvolutionPopGen,
        summary: "Median-joining haplotype networks, FASTA diversity (pi, theta_W, d_XY) & diploid VCF Weir-Cockerham F_ST.",
        exports: ".svg network",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Crispr,
        num: 17,
        short_name: "CRISPR",
        title: "CRISPR Guide & CFD Off-Target Designer",
        domain: ToolDomain::MolecularCloning,
        summary: "SpCas9/SaCas9/LbCas12a guide ranking, Doench CFD off-target scoring, CBE/ABE windows & Golden Gate oligos.",
        exports: ".tsv / .svg",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Flow,
        num: 18,
        short_name: "Flow",
        title: "Flow Cytometry FCS 3.1 & Logicle Gating",
        domain: ToolDomain::AssaysKinetics,
        summary: "Binary .fcs 3.1 parser, $SPILLOVER matrix compensation, Logicle scaling & interactive quadrant gating.",
        exports: ".svg biaxial",
        is_dedicated_studio: true,
    },
    ToolMeta {
        id: ToolId::Survival,
        num: 19,
        short_name: "Survival",
        title: "Kaplan-Meier, Cox PH & Probit LC50",
        domain: ToolDomain::EpiBiostats,
        summary: "Greenwood 95% CI survival curves, Mantel log-rank test, Cox PH hazard ratios & Finney Probit/Logit LC50.",
        exports: ".svg dual-panel",
        is_dedicated_studio: false,
    },
    ToolMeta {
        id: ToolId::Flux,
        num: 20,
        short_name: "Flux",
        title: "FvCB Photosynthesis & Stream Metabolism",
        domain: ToolDomain::EcologyField,
        summary: "Farquhar-von Caemmerer-Berry C3 A-Ci curve fitting & Odum open-channel stream GPP/ER metabolism.",
        exports: ".svg dual-panel",
        is_dedicated_studio: true,
    },
];

/// Unified desktop application wrapping all 20 OSABT scientific tools.
pub struct OsabtWorkbenchApp {
    active_tool: Option<ToolId>,
    domain_filter: Option<ToolDomain>,
    search_query: String,
    rail_open: bool,
    // 8 dedicated eframe/egui studio apps (lazy-initialized on first launch)
    plasmid: Option<PlasmidStudioApp>,
    forge3d: Option<Forge3dStudioApp>,
    morpho: Option<MorphoStudioApp>,
    ecosound: Option<EcosoundStudioApp>,
    kinetics: Option<KineticsStudioApp>,
    popgen: Option<PopGenGuiApp>,
    flow: Option<FlowGuiApp>,
    flux: Option<FluxGuiApp>,
    // 12 interactive workbench panels (lazy-initialized on first launch)
    bench: Option<BenchPanelState>,
    dwc: Option<DwcPanelState>,
    phylo: Option<PhyloPanelState>,
    epi: Option<EpiPanelState>,
    lexi: Option<LexiPanelState>,
    plate: Option<PlatePanelState>,
    infoevo: Option<InfoEvoPanelState>,
    qpcr: Option<QpcrPanelState>,
    edna: Option<EdnaPanelState>,
    wildlife: Option<WildlifePanelState>,
    crispr: Option<CrisprPanelState>,
    survival: Option<SurvivalPanelState>,
    visuals_applied: bool,
    init_error: Option<String>,
}

impl Default for OsabtWorkbenchApp {
    fn default() -> Self {
        Self::new()
    }
}

impl OsabtWorkbenchApp {
    pub fn new() -> Self {
        Self {
            active_tool: None,
            domain_filter: None,
            search_query: String::new(),
            rail_open: true,
            plasmid: None,
            forge3d: None,
            morpho: None,
            ecosound: None,
            kinetics: None,
            popgen: None,
            flow: None,
            flux: None,
            bench: None,
            dwc: None,
            phylo: None,
            epi: None,
            lexi: None,
            plate: None,
            infoevo: None,
            qpcr: None,
            edna: None,
            wildlife: None,
            crispr: None,
            survival: None,
            visuals_applied: false,
            init_error: None,
        }
    }

    pub fn select_tool(&mut self, tool: ToolId) {
        self.active_tool = Some(tool);
        self.init_error = None;
    }

    fn apply_dark_slate_theme(ctx: &egui::Context) {
        let mut v = egui::Visuals::dark();
        v.panel_fill = Color32::from_rgb(15, 23, 42);
        v.window_fill = Color32::from_rgb(20, 30, 50);
        v.extreme_bg_color = Color32::from_rgb(9, 14, 28);
        v.widgets.noninteractive.bg_fill = Color32::from_rgb(24, 34, 54);
        v.widgets.inactive.bg_fill = Color32::from_rgb(30, 41, 59);
        v.widgets.hovered.bg_fill = Color32::from_rgb(51, 65, 85);
        v.widgets.active.bg_fill = Color32::from_rgb(37, 99, 235);
        v.selection.bg_fill = Color32::from_rgb(30, 64, 175);
        ctx.set_visuals(v);
    }

    fn tool_matches_filter(&self, meta: &ToolMeta) -> bool {
        if let Some(dom) = self.domain_filter {
            if meta.domain != dom {
                return false;
            }
        }
        let q = self.search_query.trim().to_ascii_lowercase();
        if q.is_empty() {
            return true;
        }
        meta.short_name.to_ascii_lowercase().contains(&q)
            || meta.title.to_ascii_lowercase().contains(&q)
            || meta.summary.to_ascii_lowercase().contains(&q)
            || meta.domain.label().to_ascii_lowercase().contains(&q)
            || format!("#{}", meta.num).contains(&q)
    }

    fn draw_suite_header(&mut self, ctx: &egui::Context) {
        egui::TopBottomPanel::top("osabt_suite_header").show(ctx, |ui| {
            ui.add_space(3.0);
            ui.horizontal_wrapped(|ui| {
                let launcher_btn = egui::Button::new(
                    RichText::new("🏠 Suite Launcher")
                        .strong()
                        .color(Color32::WHITE),
                )
                .fill(if self.active_tool.is_none() {
                    Color32::from_rgb(37, 99, 235)
                } else {
                    Color32::from_rgb(30, 41, 59)
                });
                if ui.add(launcher_btn).clicked() {
                    self.active_tool = None;
                }

                let rail_lbl = if self.rail_open { "◧ Hide Rail" } else { "☰ Show Rail" };
                if ui.button(rail_lbl).clicked() {
                    self.rail_open = !self.rail_open;
                }

                ui.separator();
                if ui
                    .selectable_label(self.domain_filter.is_none(), "All (20)")
                    .clicked()
                {
                    self.domain_filter = None;
                }
                for dom in ToolDomain::ALL {
                    let selected = self.domain_filter == Some(dom);
                    let txt = RichText::new(dom.label()).color(if selected {
                        Color32::WHITE
                    } else {
                        dom.accent()
                    });
                    if ui.selectable_label(selected, txt).clicked() {
                        self.domain_filter = if selected { None } else { Some(dom) };
                    }
                }

                ui.separator();
                ui.add(
                    egui::TextEdit::singleline(&mut self.search_query)
                        .hint_text("Search 20 tools...")
                        .desired_width(150.0),
                );

                ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                    let pill_txt = match self.active_tool {
                        Some(id) => {
                            let m = TOOL_CATALOG.iter().find(|t| t.id == id).unwrap();
                            format!("Active: #{} {}", m.num, m.short_name)
                        }
                        None => String::from("OSABT v2.0.0 • 20 Tools"),
                    };
                    egui::Frame::none()
                        .fill(Color32::from_rgb(30, 58, 138))
                        .rounding(10.0)
                        .inner_margin(Margin::symmetric(10.0, 4.0))
                        .show(ui, |ui| {
                            ui.label(RichText::new(pill_txt).small().strong().color(Color32::WHITE));
                        });
                });
            });
            ui.add_space(2.0);
        });
    }

    fn draw_suite_rail(&mut self, ctx: &egui::Context) {
        if !self.rail_open {
            return;
        }
        egui::SidePanel::left("osabt_suite_rail")
            .default_width(205.0)
            .min_width(185.0)
            .max_width(250.0)
            .resizable(true)
            .show(ctx, |ui| {
                ui.add_space(4.0);
                ui.strong(RichText::new("OSABT TOOL SUITE").color(Color32::from_rgb(148, 163, 184)));
                ui.separator();
                egui::ScrollArea::vertical().id_salt("osabt_rail_scroll").show(ui, |ui| {
                    for meta in &TOOL_CATALOG {
                        if !self.tool_matches_filter(meta) {
                            continue;
                        }
                        let is_active = self.active_tool == Some(meta.id);
                        let label = format!("#{:02} {}", meta.num, meta.short_name);
                        let rt = RichText::new(label)
                            .strong()
                            .color(if is_active {
                                Color32::WHITE
                            } else {
                                Color32::from_rgb(203, 213, 225)
                            });
                        if ui
                            .add_sized(
                                [ui.available_width(), 24.0],
                                egui::Button::new(rt).selected(is_active),
                            )
                            .on_hover_text(meta.title)
                            .clicked()
                        {
                            self.select_tool(meta.id);
                        }
                    }
                });
            });
    }

    fn draw_home_launcher(&mut self, ctx: &egui::Context) {
        egui::CentralPanel::default().show(ctx, |ui| {
            egui::ScrollArea::vertical().id_salt("osabt_launcher_scroll").show(ui, |ui| {
                ui.add_space(6.0);
                ui.horizontal(|ui| {
                    ui.heading(
                        RichText::new("Open Source Ad-hoc Biology Tools (OSABT v2.0.0)")
                            .size(22.0)
                            .strong()
                            .color(Color32::from_rgb(241, 245, 249)),
                    );
                    ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                        ui.label(
                            RichText::new("CRV Science Network • GPL-3.0 • Pure-Rust Native Suite")
                                .color(Color32::from_rgb(148, 163, 184)),
                        );
                    });
                });
                ui.label(
                    RichText::new(
                        "Select any domain filter above or click Launch on any of the 20 integrated molecular, evolutionary, kinetic, ecological, and biostatistical workbenches.",
                    )
                    .color(Color32::from_rgb(148, 163, 184)),
                );
                ui.add_space(10.0);

                for dom in ToolDomain::ALL {
                    let matching: Vec<&ToolMeta> = TOOL_CATALOG
                        .iter()
                        .filter(|m| m.domain == dom && self.tool_matches_filter(m))
                        .collect();
                    if matching.is_empty() {
                        continue;
                    }

                    ui.horizontal(|ui| {
                        ui.colored_label(dom.accent(), RichText::new(format!("■ {}", dom.label())).size(15.5).strong());
                        ui.label(RichText::new(format!("({} tools)", matching.len())).small().color(Color32::from_rgb(148, 163, 184)));
                    });
                    ui.add_space(4.0);

                    let avail_w = ui.available_width();
                    let cols = if avail_w >= 1100.0 { 3 } else if avail_w >= 720.0 { 2 } else { 1 };

                    for chunk in matching.chunks(cols) {
                        ui.columns(cols, |col_uis| {
                            for (ci, meta) in chunk.iter().enumerate() {
                                let col_ui = &mut col_uis[ci];
                                let inner_w = (col_ui.available_width() - 28.0).max(180.0);
                                egui::Frame::none()
                                    .fill(Color32::from_rgb(24, 34, 54))
                                    .stroke(Stroke::new(1.0_f32, Color32::from_rgb(51, 65, 85)))
                                    .rounding(8.0)
                                    .inner_margin(Margin::same(12.0))
                                    .show(col_ui, |card| {
                                        card.set_min_width(inner_w);
                                        card.set_max_width(inner_w);
                                        card.set_min_height(114.0);
                                        card.horizontal(|top| {
                                            top.label(
                                                RichText::new(format!("#{} {}", meta.num, meta.short_name))
                                                    .strong()
                                                    .size(14.5)
                                                    .color(dom.accent()),
                                            );
                                            top.with_layout(egui::Layout::right_to_left(egui::Align::Center), |r| {
                                                let tag = if meta.is_dedicated_studio { "Studio" } else { "Workbench" };
                                                r.label(RichText::new(tag).small().color(Color32::from_rgb(148, 163, 184)));
                                            });
                                        });
                                        card.add(
                                            egui::Label::new(
                                                RichText::new(meta.title)
                                                    .strong()
                                                    .color(Color32::from_rgb(241, 245, 249)),
                                            )
                                            .wrap(),
                                        );
                                        card.add_space(2.0);
                                        card.add(
                                            egui::Label::new(
                                                RichText::new(meta.summary)
                                                    .small()
                                                    .color(Color32::from_rgb(148, 163, 184)),
                                            )
                                            .wrap(),
                                        );
                                        card.add_space(6.0);
                                        card.horizontal(|bot| {
                                            bot.label(
                                                RichText::new(format!("Exports: {}", meta.exports))
                                                    .small()
                                                    .color(Color32::from_rgb(100, 116, 139)),
                                            );
                                            bot.with_layout(egui::Layout::right_to_left(egui::Align::Center), |r| {
                                                let btn = egui::Button::new(
                                                    RichText::new("Launch ->").small().strong().color(Color32::WHITE),
                                                )
                                                .fill(Color32::from_rgb(30, 64, 175));
                                                if r.add(btn).clicked() {
                                                    self.select_tool(meta.id);
                                                }
                                            });
                                        });
                                    });
                            }
                        });
                        ui.add_space(8.0);
                    }
                    ui.add_space(6.0);
                }
            });
        });
    }

    fn draw_active_tool(&mut self, ctx: &egui::Context, tool: ToolId) {
        if let Some(err) = &self.init_error {
            egui::CentralPanel::default().show(ctx, |ui| {
                ui.colored_label(Color32::from_rgb(248, 113, 113), err);
            });
            return;
        }
        match tool {
            ToolId::Bench => self.bench.get_or_insert_with(BenchPanelState::new).draw(ctx),
            ToolId::Plasmid => self.plasmid.get_or_insert_with(|| PlasmidStudioApp::new(None)).draw_panels(ctx),
            ToolId::Dwc => self.dwc.get_or_insert_with(DwcPanelState::new).draw(ctx),
            ToolId::Phylo => self.phylo.get_or_insert_with(PhyloPanelState::new).draw(ctx),
            ToolId::Epi => self.epi.get_or_insert_with(EpiPanelState::new).draw(ctx),
            ToolId::Lexi => self.lexi.get_or_insert_with(LexiPanelState::new).draw(ctx),
            ToolId::Forge3d => self.forge3d.get_or_insert_with(|| Forge3dStudioApp::new(None)).draw_panels(ctx),
            ToolId::Morpho => self.morpho.get_or_insert_with(|| MorphoStudioApp::new(None)).draw_panels(ctx),
            ToolId::Plate => self.plate.get_or_insert_with(PlatePanelState::new).draw(ctx),
            ToolId::Ecosound => self.ecosound.get_or_insert_with(|| EcosoundStudioApp::new(None)).draw_panels(ctx),
            ToolId::Infoevo => self.infoevo.get_or_insert_with(InfoEvoPanelState::new).draw(ctx),
            ToolId::Qpcr => self.qpcr.get_or_insert_with(QpcrPanelState::new).draw(ctx),
            ToolId::Kinetics => self.kinetics.get_or_insert_with(|| KineticsStudioApp::new(None)).draw_panels(ctx),
            ToolId::Edna => self.edna.get_or_insert_with(EdnaPanelState::new).draw(ctx),
            ToolId::Wildlife => self.wildlife.get_or_insert_with(WildlifePanelState::new).draw(ctx),
            ToolId::Popgen => {
                if self.popgen.is_none() {
                    match PopGenGuiApp::new(None, None) {
                        Ok(a) => self.popgen = Some(a),
                        Err(e) => { self.init_error = Some(e); return; }
                    }
                }
                if let Some(a) = &mut self.popgen { a.draw_panels(ctx); }
            }
            ToolId::Crispr => self.crispr.get_or_insert_with(CrisprPanelState::new).draw(ctx),
            ToolId::Flow => {
                if self.flow.is_none() {
                    match FlowGuiApp::new(None) {
                        Ok(a) => self.flow = Some(a),
                        Err(e) => { self.init_error = Some(e); return; }
                    }
                }
                if let Some(a) = &mut self.flow { a.draw_panels(ctx); }
            }
            ToolId::Survival => self.survival.get_or_insert_with(SurvivalPanelState::new).draw(ctx),
            ToolId::Flux => {
                if self.flux.is_none() {
                    match FluxGuiApp::new(None, None) {
                        Ok(mut a) => { a.set_dual_overview(); self.flux = Some(a); }
                        Err(e) => { self.init_error = Some(e); return; }
                    }
                }
                if let Some(a) = &mut self.flux { a.draw_panels(ctx); }
            }
        }
    }

    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        if !self.visuals_applied {
            Self::apply_dark_slate_theme(ctx);
            self.visuals_applied = true;
        }
        self.draw_suite_header(ctx);
        self.draw_suite_rail(ctx);
        if let Some(tool) = self.active_tool {
            self.draw_active_tool(ctx, tool);
        } else {
            self.draw_home_launcher(ctx);
        }
    }
}

impl eframe::App for OsabtWorkbenchApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.draw_panels(ctx);
    }
}

/// Launches the unified 20-tool `OSABT Workbench` desktop GUI.
pub fn run_workbench_gui() -> Result<(), String> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1480.0, 920.0])
            .with_min_inner_size([960.0, 640.0])
            .with_title("OSABT v2.0.0 — 20-Tool Biology & Biostatistics Suite (GPL-3.0)"),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-workbench",
        options,
        Box::new(|cc| {
            OsabtWorkbenchApp::apply_dark_slate_theme(&cc.egui_ctx);
            Ok(Box::new(OsabtWorkbenchApp::new()))
        }),
    )
    .map_err(|e| format!("Failed launching OSABT Workbench GUI: {e}"))
}

/// Captures a 1600x1000 offscreen GPU screenshot of OSABT Workbench views 1..=3.
pub fn run_workbench_screenshot(out_path: &Path, view_idx: u8) -> Result<(), String> {
    let mut app = OsabtWorkbenchApp::new();
    match view_idx {
        1 => {
            app.active_tool = None;
            app.rail_open = true;
        }
        2 => {
            app.select_tool(ToolId::Survival);
            app.rail_open = true;
        }
        3 => {
            app.select_tool(ToolId::Plasmid);
            app.rail_open = true;
        }
        other => {
            return Err(format!(
                "Unsupported --screenshot-view {other} (expected 1 = Suite Launcher, 2 = Survival GUI, 3 = Plasmid Studio)"
            ));
        }
    }
    let ctx = egui::Context::default();
    render_flux_offscreen_png(&ctx, out_path, |c| app.draw_panels(c))
}
