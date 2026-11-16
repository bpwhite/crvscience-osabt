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

//! Interactive 4-tab `Ratatui` terminal eDNA & Metabarcoding FASTQ Demultiplexer,
//! ASV Denoiser, & Chao1 / Hill-Number Rarefaction Studio (`osabt edna`).

use std::fs;
use std::io;
use std::path::PathBuf;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::Line;
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::Terminal;

use crate::tui_util::run_tui_session;

use crate::edna_diversity::{
    analyze_edna_community_diversity, render_edna_rarefaction_svg, CommunityDiversityReport,
};
use crate::edna_fastq::{
    demo_pnw_12s_edna_dataset, export_asv_table_tsv, parse_barcode_sheet_csv, parse_fastq_str,
    run_edna_metabarcoding_pipeline, DualBarcodeSample, EdnaPipelineConfig, EdnaPipelineReport,
    FastqRead,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EdnaTab {
    PipelineQc,
    AsvCommunityTable,
    AlphaHillDiversity,
    RarefactionAndBeta,
}

struct EdnaTuiState {
    active_tab: EdnaTab,
    reads: Vec<FastqRead>,
    barcodes: Vec<DualBarcodeSample>,
    ref_db: Vec<(String, String)>,
    cfg: EdnaPipelineConfig,
    pipeline: EdnaPipelineReport,
    diversity: CommunityDiversityReport,
    scroll_offset: usize,
    status_msg: String,
}

impl EdnaTuiState {
    fn new(fastq_path: Option<PathBuf>, barcodes_path: Option<PathBuf>) -> Result<Self, String> {
        let (demo_fq, demo_bc, demo_db) = demo_pnw_12s_edna_dataset();
        let reads = if let Some(fp) = fastq_path {
            let raw = fs::read_to_string(&fp)
                .map_err(|e| format!("Failed reading FASTQ '{}': {e}", fp.display()))?;
            parse_fastq_str(&raw)?
        } else {
            parse_fastq_str(&demo_fq)?
        };

        let barcodes = if let Some(bp) = barcodes_path {
            let raw = fs::read_to_string(&bp)
                .map_err(|e| format!("Failed reading barcode sheet '{}': {e}", bp.display()))?;
            parse_barcode_sheet_csv(&raw)?
        } else {
            parse_barcode_sheet_csv(&demo_bc)?
        };

        let cfg = EdnaPipelineConfig::default();
        let pipeline = run_edna_metabarcoding_pipeline(&reads, &barcodes, &demo_db, &cfg);
        let diversity = analyze_edna_community_diversity(&pipeline, 10);

        Ok(Self {
            active_tab: EdnaTab::PipelineQc,
            reads,
            barcodes,
            ref_db: demo_db,
            cfg,
            pipeline,
            diversity,
            scroll_offset: 0,
            status_msg: String::from(
                "Keys: [1-4/Tab] Views | [Up/Down] Scroll | [+/-] MaxEE Threshold | [c] Toggle Chimera Screen | [e] Export TSV + SVG | [q] Quit",
            ),
        })
    }

    fn rerun(&mut self) {
        self.pipeline =
            run_edna_metabarcoding_pipeline(&self.reads, &self.barcodes, &self.ref_db, &self.cfg);
        self.diversity = analyze_edna_community_diversity(&self.pipeline, 10);
    }

    fn export_files(&mut self) {
        let tsv = export_asv_table_tsv(&self.pipeline, false);
        let svg = render_edna_rarefaction_svg(&self.diversity, "PNW 12S eDNA Biodiversity");
        if fs::write("edna_asv_table.tsv", tsv).is_ok()
            && fs::write("edna_rarefaction.svg", svg).is_ok()
        {
            self.status_msg =
                String::from("Exported edna_asv_table.tsv and edna_rarefaction.svg!");
        } else {
            self.status_msg = String::from("Export failed: check directory write permissions.");
        }
    }
}

/// Launches the interactive 4-tab `Ratatui` terminal studio for `osabt edna`.
pub fn run_edna_tui(
    fastq_path: Option<PathBuf>,
    barcodes_path: Option<PathBuf>,
) -> Result<(), String> {
    let mut state = EdnaTuiState::new(fastq_path, barcodes_path)?;
    run_tui_session(|terminal| run_loop(terminal, &mut state))
}

fn run_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut EdnaTuiState,
) -> Result<(), String> {
    loop {
        terminal
            .draw(|f| draw_ui(f, state))
            .map_err(|e| format!("Terminal draw error: {e}"))?;

        if event::poll(Duration::from_millis(120)).map_err(|e| e.to_string())? {
            if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }
                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                    KeyCode::Char('1') => {
                        state.active_tab = EdnaTab::PipelineQc;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('2') => {
                        state.active_tab = EdnaTab::AsvCommunityTable;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('3') => {
                        state.active_tab = EdnaTab::AlphaHillDiversity;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('4') => {
                        state.active_tab = EdnaTab::RarefactionAndBeta;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            EdnaTab::PipelineQc => EdnaTab::AsvCommunityTable,
                            EdnaTab::AsvCommunityTable => EdnaTab::AlphaHillDiversity,
                            EdnaTab::AlphaHillDiversity => EdnaTab::RarefactionAndBeta,
                            EdnaTab::RarefactionAndBeta => EdnaTab::PipelineQc,
                        };
                        state.scroll_offset = 0;
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        state.scroll_offset = state.scroll_offset.saturating_add(1);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.scroll_offset = state.scroll_offset.saturating_sub(1);
                    }
                    KeyCode::Char('+') | KeyCode::Char('=') => {
                        state.cfg.max_expected_errors =
                            (state.cfg.max_expected_errors + 0.25).min(10.0);
                        state.rerun();
                        state.status_msg = format!(
                            "Max Expected Errors (Edgar & Flyvbjerg 2015) = {:.2}",
                            state.cfg.max_expected_errors
                        );
                    }
                    KeyCode::Char('-') => {
                        state.cfg.max_expected_errors =
                            (state.cfg.max_expected_errors - 0.25).max(0.10);
                        state.rerun();
                        state.status_msg = format!(
                            "Max Expected Errors (Edgar & Flyvbjerg 2015) = {:.2}",
                            state.cfg.max_expected_errors
                        );
                    }
                    KeyCode::Char('c') => {
                        state.cfg.remove_chimeras = !state.cfg.remove_chimeras;
                        state.rerun();
                        state.status_msg = format!(
                            "De novo 2-parent chimera screening: {}",
                            if state.cfg.remove_chimeras { "ON" } else { "OFF" }
                        );
                    }
                    KeyCode::Char('e') => state.export_files(),
                    _ => {}
                }
            }
        }
    }
}

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &EdnaTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Demux & MaxEE QC",
        "[2] ASV Table & Blank Audit",
        "[3] Chao1 & Hill (^0D, ^1D, ^2D)",
        "[4] Hurlbert Rarefaction & Beta",
    ]
    .into_iter()
    .map(Line::from)
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        EdnaTab::PipelineQc => 0,
        EdnaTab::AsvCommunityTable => 1,
        EdnaTab::AlphaHillDiversity => 2,
        EdnaTab::RarefactionAndBeta => 3,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT eDNA Studio (v1.4.0) — {} Raw Reads → {} Clean ASVs across {} Env Sites ",
            state.pipeline.total_raw_reads,
            state.diversity.clean_asv_count,
            state.diversity.sample_names.len()
        )))
        .select(idx)
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        EdnaTab::PipelineQc => draw_qc_tab(f, chunks[1], state),
        EdnaTab::AsvCommunityTable => draw_asv_tab(f, chunks[1], state),
        EdnaTab::AlphaHillDiversity => draw_alpha_tab(f, chunks[1], state),
        EdnaTab::RarefactionAndBeta => draw_rarefaction_beta_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls & Status "))
        .style(Style::default().fg(Color::Gray));
    f.render_widget(footer, chunks[2]);
}

fn draw_qc_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EdnaTuiState) {
    let p = &state.pipeline;
    let lines = vec![
        Line::from("Metabarcoding Read-Retention & Denoising Waterfall:"),
        Line::from(""),
        Line::from(format!("  1. Total Raw FASTQ Reads                 : {}", p.total_raw_reads)),
        Line::from(format!(
            "  2. Dual-Barcode Demux Assigned (mm <= {}) : {} ({} unassigned)",
            state.cfg.max_barcode_mismatches, p.demux_assigned_reads, p.unassigned_reads
        )),
        Line::from(format!(
            "  3. Passed Edgar–Flyvbjerg maxEE <= {:.2}  : {} ({} low-Q / short rejected)",
            state.cfg.max_expected_errors, p.passed_ee_reads, p.filtered_ee_or_short_reads
        )),
        Line::from(format!(
            "  4. 1-nt Error-Cloud Reads Absorbed       : {} (UNOISE skew >= {:.1}x)",
            p.error_corrected_reads, state.cfg.unoise_abundance_skew
        )),
        Line::from(format!(
            "  5. De Novo 2-Parent Chimeric Reads       : {} (screening {})",
            p.chimeric_reads_flagged,
            if state.cfg.remove_chimeras { "ON" } else { "OFF" }
        )),
        Line::from(format!(
            "  6. Total Denoised ASVs (Clean Biological): {} total ({} clean after chimera + blank filter)",
            p.asvs.len(),
            state.diversity.clean_asv_count
        )),
    ];

    let card = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" FASTQ Demultiplexing, Expected-Error Filtering & ASV Denoising Summary "),
        );
    f.render_widget(card, area);
}

fn draw_asv_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EdnaTuiState) {
    let mut cols = vec![
        String::from("ASV"),
        String::from("Status"),
        String::from("Total"),
    ];
    for s in &state.pipeline.sample_names {
        cols.push(s.clone());
    }
    cols.push(String::from("12S Taxonomic Assignment"));

    let hdr = Row::new(cols).style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state
        .pipeline
        .asvs
        .iter()
        .skip(state.scroll_offset)
        .map(|a| {
            let (status, col) = if a.is_chimera {
                ("CHIMERA", Color::LightRed)
            } else if a.blank_flagged {
                ("BLANK_WARN", Color::Yellow)
            } else {
                ("PASS", Color::Green)
            };
            let mut cells = vec![
                Cell::from(a.asv_id.clone()),
                Cell::from(status).style(Style::default().fg(col).add_modifier(Modifier::BOLD)),
                Cell::from(format!("{}", a.total_reads)),
            ];
            for &c in &a.counts_by_sample {
                cells.push(Cell::from(format!("{c}")));
            }
            cells.push(Cell::from(
                a.putative_taxon
                    .clone()
                    .unwrap_or_else(|| String::from("Unassigned")),
            ));
            Row::new(cells)
        });

    let mut widths = vec![
        Constraint::Length(9),
        Constraint::Length(11),
        Constraint::Length(7),
    ];
    for _ in &state.pipeline.sample_names {
        widths.push(Constraint::Length(15));
    }
    widths.push(Constraint::Min(28));

    let tbl = Table::new(rows, widths).header(hdr).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Denoised ASV × Sample Matrix & Extraction-Blank Contamination Audit "),
    );
    f.render_widget(tbl, area);
}

fn draw_alpha_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EdnaTuiState) {
    let hdr = Row::new([
        "Environmental Site",
        "Reads N",
        "^0D (S_obs)",
        "Chao1 (95% CI)",
        "^1D exp(H')",
        "^2D (1/D)",
        "Pielou J'",
        "Hurlbert PIE",
        "Coverage",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state.diversity.alpha.iter().map(|a| {
        Row::new([
            a.sample_name.clone(),
            format!("{}", a.total_reads_n),
            format!("{}", a.observed_richness_q0),
            format!(
                "{:.2} ({:.1}–{:.1})",
                a.chao1_asymptotic, a.chao1_ci95.0, a.chao1_ci95.1
            ),
            format!("{:.2}", a.hill_shannon_q1),
            format!("{:.2}", a.hill_simpson_q2),
            format!("{:.3}", a.pielou_evenness_j),
            format!("{:.3}", a.hurlbert_pie),
            format!("{:.1}%", a.sample_coverage * 100.0),
        ])
    });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(18),
            Constraint::Length(9),
            Constraint::Length(12),
            Constraint::Length(20),
            Constraint::Length(12),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(13),
            Constraint::Min(9),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Chao (1984) Asymptotic Richness & Chao et al. (2014) Hill Numbers "),
    );
    f.render_widget(tbl, area);
}

fn draw_rarefaction_beta_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EdnaTuiState) {
    let split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(54), Constraint::Percentage(46)])
        .split(area);

    let r_hdr = Row::new(["Site", "m=10%", "m=30%", "m=50%", "m=70%", "m=100% (S_obs)"])
        .style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    let r_rows = state.diversity.rarefaction_curves.iter().map(|c| {
        let pick = |frac: f64| -> String {
            let idx = (((c.expected_species_sm.len().saturating_sub(1)) as f64) * frac).round()
                as usize;
            c.expected_species_sm
                .get(idx)
                .map(|v| format!("{v:.2}"))
                .unwrap_or_else(|| String::from("-"))
        };
        Row::new([
            c.sample_name.clone(),
            pick(0.1),
            pick(0.3),
            pick(0.5),
            pick(0.7),
            pick(1.0),
        ])
    });
    let r_tbl = Table::new(
        r_rows,
        [
            Constraint::Length(16),
            Constraint::Length(8),
            Constraint::Length(8),
            Constraint::Length(8),
            Constraint::Length(8),
            Constraint::Min(14),
        ],
    )
    .header(r_hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Hurlbert (1971) Exact Rarefaction E[S_m] "),
    );
    f.render_widget(r_tbl, split[0]);

    let mut b_cols = vec![String::from("Bray–Curtis")];
    for s in &state.diversity.sample_names {
        b_cols.push(s.clone());
    }
    let b_hdr = Row::new(b_cols).style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );
    let b_rows = state
        .diversity
        .sample_names
        .iter()
        .enumerate()
        .map(|(i, name)| {
            let mut cells = vec![name.clone()];
            for &d in &state.diversity.bray_curtis_matrix[i] {
                cells.push(format!("{d:.3}"));
            }
            Row::new(cells)
        });
    let mut b_widths = vec![Constraint::Length(16)];
    for _ in &state.diversity.sample_names {
        b_widths.push(Constraint::Length(15));
    }
    let b_tbl = Table::new(b_rows, b_widths).header(b_hdr).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Pairwise Bray–Curtis Dissimilarity "),
    );
    f.render_widget(b_tbl, split[1]);
}
