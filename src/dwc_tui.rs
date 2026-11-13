// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Interactive `Ratatui` terminal field-survey workbench for Tool #3 (`osabt dwc`).
//! Lets field biologists inspect occurrence rows, step through `CoordinateCleaner` & TDWG
//! diagnostics, test live taxonomic synonym resolution, and export clean DwC-A / SVG bundles.

use std::fs;
use std::io;
use std::path::PathBuf;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use crossterm::execute;
use crossterm::terminal::{
    disable_raw_mode, enable_raw_mode, EnterAlternateScreen, LeaveAlternateScreen,
};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::Terminal;

use crate::dwc_core::{lint_and_clean_dwc, DwcCleanConfig, DwcDataset, LintSeverity};
use crate::dwc_export::{
    export_dwc_geojson, generate_dwca_eml_xml, generate_dwca_meta_xml,
    render_dwc_occurrence_map_svg, write_dwc_delimited,
};
use crate::dwc_taxonomy::{resolve_scientific_name, TaxonBackbone};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum DwcTab {
    Records,
    Diagnostics,
    TaxonLookup,
    SpatialExport,
}

struct DwcTuiState {
    dataset: DwcDataset,
    backbone: TaxonBackbone,
    active_tab: DwcTab,
    selected_row: usize,
    selected_issue: usize,
    taxon_query: String,
    status_msg: String,
}

impl DwcTuiState {
    fn new(initial_csv: Option<PathBuf>, checklist_path: Option<PathBuf>) -> Result<Self, String> {
        let mut backbone = TaxonBackbone::default_curated();
        if let Some(cp) = checklist_path {
            let txt = fs::read_to_string(&cp)
                .map_err(|e| format!("Failed reading checklist '{}': {e}", cp.display()))?;
            backbone.load_checklist_delimited(&txt)?;
        }

        let (name, raw_csv) = match initial_csv {
            Some(p) => {
                let txt = fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading input '{}': {e}", p.display()))?;
                let stem = p
                    .file_stem()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| String::from("Field_Survey"));
                (stem, txt)
            }
            None => (
                String::from("PNW_Oregon_Demo_Survey"),
                DwcDataset::demo_pnw_field_survey_csv().to_string(),
            ),
        };

        let dataset = lint_and_clean_dwc(&name, &raw_csv, &backbone, &DwcCleanConfig::default())?;
        Ok(Self {
            dataset,
            backbone,
            active_tab: DwcTab::Records,
            selected_row: 0,
            selected_issue: 0,
            taxon_query: String::from("Pseudotsuga taxifolia (Poir.) Britton"),
            status_msg: String::from(
                "Tab/1-4: switch views | Up/Down: navigate | [e] Export clean DwC-A + GeoJSON + SVG | [q] Quit",
            ),
        })
    }

    fn export_bundle(&mut self) {
        let out_dir = PathBuf::from(format!("{}_dwca", self.dataset.dataset_name));
        if let Err(e) = fs::create_dir_all(&out_dir) {
            self.status_msg = format!("Export failed creating dir: {e}");
            return;
        }
        let csv_ok = fs::write(
            out_dir.join("occurrence.csv"),
            write_dwc_delimited(&self.dataset, ',', true),
        );
        let meta_ok = fs::write(
            out_dir.join("meta.xml"),
            generate_dwca_meta_xml(&self.dataset, "occurrence.csv", Some("eml.xml")),
        );
        let eml_ok = fs::write(
            out_dir.join("eml.xml"),
            generate_dwca_eml_xml(&self.dataset.dataset_name, "CRV Field Team", "CC-BY-4.0"),
        );
        let geo_ok = fs::write(
            out_dir.join("occurrences.geojson"),
            export_dwc_geojson(&self.dataset),
        );
        let svg_ok = fs::write(
            out_dir.join("spatial_qa_map.svg"),
            render_dwc_occurrence_map_svg(&self.dataset),
        );

        if csv_ok.is_ok() && meta_ok.is_ok() && eml_ok.is_ok() && geo_ok.is_ok() && svg_ok.is_ok() {
            self.status_msg = format!(
                "Exported DwC-A (occurrence.csv, meta.xml, eml.xml, GeoJSON, SVG) -> ./{}/",
                out_dir.display()
            );
        } else {
            self.status_msg = String::from("One or more files failed writing during DwC-A export.");
        }
    }
}

/// Launches the interactive `Ratatui` terminal workbench for `osabt dwc`.
pub fn run_dwc_tui(
    initial_csv: Option<PathBuf>,
    checklist_path: Option<PathBuf>,
) -> Result<(), String> {
    let mut state = DwcTuiState::new(initial_csv, checklist_path)?;

    enable_raw_mode().map_err(|e| format!("Failed enabling terminal raw mode: {e}"))?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)
        .map_err(|e| format!("Failed entering alternate screen: {e}"))?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal =
        Terminal::new(backend).map_err(|e| format!("Failed creating Ratatui terminal: {e}"))?;

    let loop_res = run_event_loop(&mut terminal, &mut state);

    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen);
    let _ = terminal.show_cursor();
    loop_res
}

fn run_event_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut DwcTuiState,
) -> Result<(), String> {
    loop {
        terminal
            .draw(|f| draw_ui(f, state))
            .map_err(|e| format!("Draw error: {e}"))?;

        if event::poll(Duration::from_millis(120)).map_err(|e| e.to_string())? {
            if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }
                if state.active_tab == DwcTab::TaxonLookup {
                    match key.code {
                        KeyCode::Esc => state.active_tab = DwcTab::Records,
                        KeyCode::Tab => state.active_tab = DwcTab::SpatialExport,
                        KeyCode::Backspace => {
                            state.taxon_query.pop();
                        }
                        KeyCode::Char(c) => {
                            state.taxon_query.push(c);
                        }
                        _ => {}
                    }
                    continue;
                }

                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                    KeyCode::Char('1') => state.active_tab = DwcTab::Records,
                    KeyCode::Char('2') => state.active_tab = DwcTab::Diagnostics,
                    KeyCode::Char('3') => state.active_tab = DwcTab::TaxonLookup,
                    KeyCode::Char('4') => state.active_tab = DwcTab::SpatialExport,
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            DwcTab::Records => DwcTab::Diagnostics,
                            DwcTab::Diagnostics => DwcTab::TaxonLookup,
                            DwcTab::TaxonLookup => DwcTab::SpatialExport,
                            DwcTab::SpatialExport => DwcTab::Records,
                        };
                    }
                    KeyCode::Char('e') => state.export_bundle(),
                    KeyCode::Down | KeyCode::Char('j') => match state.active_tab {
                        DwcTab::Records => {
                            if !state.dataset.records.is_empty() {
                                state.selected_row =
                                    (state.selected_row + 1) % state.dataset.records.len();
                            }
                        }
                        DwcTab::Diagnostics => {
                            if !state.dataset.issues.is_empty() {
                                state.selected_issue =
                                    (state.selected_issue + 1) % state.dataset.issues.len();
                            }
                        }
                        _ => {}
                    },
                    KeyCode::Up | KeyCode::Char('k') => match state.active_tab {
                        DwcTab::Records => {
                            if !state.dataset.records.is_empty() {
                                state.selected_row = state
                                    .selected_row
                                    .checked_sub(1)
                                    .unwrap_or(state.dataset.records.len() - 1);
                            }
                        }
                        DwcTab::Diagnostics => {
                            if !state.dataset.issues.is_empty() {
                                state.selected_issue = state
                                    .selected_issue
                                    .checked_sub(1)
                                    .unwrap_or(state.dataset.issues.len() - 1);
                            }
                        }
                        _ => {}
                    },
                    _ => {}
                }
            }
        }
    }
}

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &DwcTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let summary = state.dataset.summary();
    let titles = [
        "[1] Occurrences",
        "[2] QA Diagnostics",
        "[3] Taxon Synonym Lookup",
        "[4] Spatial Grid & DwC-A",
    ]
    .iter()
    .map(|t| Line::from(*t))
    .collect::<Vec<_>>();

    let tab_idx = match state.active_tab {
        DwcTab::Records => 0,
        DwcTab::Diagnostics => 1,
        DwcTab::TaxonLookup => 2,
        DwcTab::SpatialExport => 3,
    };

    let tabs = Tabs::new(titles)
        .select(tab_idx)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT Darwin Core Studio • {} ({} rows | {} auto-fixed | {} errors) ",
            state.dataset.dataset_name,
            summary.total_records,
            summary.auto_fixed_count,
            summary.error_count
        )))
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        DwcTab::Records => draw_records_tab(f, chunks[1], state),
        DwcTab::Diagnostics => draw_diagnostics_tab(f, chunks[1], state),
        DwcTab::TaxonLookup => draw_taxon_tab(f, chunks[1], state),
        DwcTab::SpatialExport => draw_spatial_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls "));
    f.render_widget(footer, chunks[2]);
}

fn draw_records_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &DwcTuiState) {
    let header = Row::new([
        "Row",
        "occurrenceID",
        "scientificName",
        "acceptedNameUsage",
        "eventDate",
        "Lat",
        "Lon",
        "Status",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state.dataset.records.iter().enumerate().map(|(i, r)| {
        let status = if r.has_errors {
            ("ERROR", Color::Red)
        } else if r.has_warnings {
            ("WARN", Color::Yellow)
        } else {
            ("CLEAN", Color::Green)
        };
        let row_style = if i == state.selected_row {
            Style::default()
                .bg(Color::DarkGray)
                .add_modifier(Modifier::BOLD)
        } else {
            Style::default()
        };
        Row::new(vec![
            Cell::from(r.row_number_1.to_string()),
            Cell::from(r.get("occurrenceID").to_string()),
            Cell::from(r.get("scientificName").to_string()),
            Cell::from(r.get("acceptedNameUsage").to_string()),
            Cell::from(r.get("eventDate").to_string()),
            Cell::from(r.get("decimalLatitude").to_string()),
            Cell::from(r.get("decimalLongitude").to_string()),
            Cell::from(Span::styled(status.0, Style::default().fg(status.1))),
        ])
        .style(row_style)
    });

    let widths = [
        Constraint::Length(4),
        Constraint::Length(22),
        Constraint::Length(24),
        Constraint::Length(22),
        Constraint::Length(11),
        Constraint::Length(11),
        Constraint::Length(12),
        Constraint::Length(7),
    ];

    let table = Table::new(rows, widths).header(header).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Scrubbed Simple Darwin Core Records "),
    );
    f.render_widget(table, area);
}

fn draw_diagnostics_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &DwcTuiState) {
    let header = Row::new(["Row", "Sev", "Code", "Term", "Fixed", "Message"]).style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state.dataset.issues.iter().enumerate().map(|(i, iss)| {
        let col = match iss.severity {
            LintSeverity::Error => Color::Red,
            LintSeverity::Warning => Color::Yellow,
            LintSeverity::Info => Color::Cyan,
        };
        let style = if i == state.selected_issue {
            Style::default().bg(Color::DarkGray)
        } else {
            Style::default()
        };
        Row::new(vec![
            Cell::from(iss.row_number_1.to_string()),
            Cell::from(Span::styled(iss.severity.label(), Style::default().fg(col))),
            Cell::from(iss.code.as_str()),
            Cell::from(iss.column_term.clone()),
            Cell::from(if iss.auto_applied { "yes" } else { "no" }),
            Cell::from(iss.message.clone()),
        ])
        .style(style)
    });

    let widths = [
        Constraint::Length(4),
        Constraint::Length(6),
        Constraint::Length(22),
        Constraint::Length(20),
        Constraint::Length(6),
        Constraint::Min(30),
    ];

    let table = Table::new(rows, widths).header(header).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" TDWG & CoordinateCleaner Audit Findings "),
    );
    f.render_widget(table, area);
}

fn draw_taxon_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &DwcTuiState) {
    let res = resolve_scientific_name(&state.taxon_query, &state.backbone);
    let mut lines = vec![
        Line::from(vec![
            Span::styled("Type Scientific Name (Esc to leave input): ", Style::default().fg(Color::Yellow)),
            Span::styled(&state.taxon_query, Style::default().fg(Color::White).add_modifier(Modifier::BOLD)),
        ]),
        Line::from(""),
        Line::from(format!("Canonical Binomial : {}", res.parsed.canonical_name)),
        Line::from(format!("Extracted Author   : {}", if res.parsed.authorship.is_empty() { "(none)" } else { &res.parsed.authorship })),
        Line::from(format!("Inferred Rank      : {}", res.parsed.inferred_rank)),
        Line::from(format!("Match Status       : {:?}", res.match_kind)),
        Line::from(format!("GBIF REST Endpoint : {}", res.gbif_lookup_url)),
    ];

    if let Some(m) = res.matched_entry {
        lines.push(Line::from(""));
        lines.push(Line::from(Span::styled(
            format!(
                "Accepted Taxon: {} {} [GBIF:{}]",
                m.accepted_name, m.authorship, m.gbif_usage_key
            ),
            Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(format!(
            "Hierarchy: {} > {} > {} > {} > {} > {}",
            m.kingdom, m.phylum, m.class_name, m.order, m.family, m.genus
        )));
    }

    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Offline Taxonomic Synonym & Orthographic Resolver "),
        );
    f.render_widget(p, area);
}

fn draw_spatial_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &DwcTuiState) {
    let s = state.dataset.summary();
    let lines = vec![
        Line::from(Span::styled(
            "Darwin Core Archive (DwC-A) & GIS Readiness Summary",
            Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD),
        )),
        Line::from(""),
        Line::from(format!("• Total Occurrence Rows     : {}", s.total_records)),
        Line::from(format!("• Clean Publication Rows    : {}", s.clean_records)),
        Line::from(format!("• Warning Rows              : {}", s.records_with_warnings)),
        Line::from(format!("• Excluded Error Rows       : {}", s.records_with_errors)),
        Line::from(format!("• Valid WGS84 Coordinates   : {}", s.georeferenced_records)),
        Line::from(format!("• Synonyms Auto-Resolved    : {}", s.synonyms_resolved)),
        Line::from(format!("• Distinct Scientific Names : {}", s.unique_taxa)),
        Line::from(""),
        Line::from("Press [e] to write full GBIF IPT-ready bundle:"),
        Line::from("  1. occurrence.csv        (RFC-4180 scrubbed Simple Darwin Core)"),
        Line::from("  2. meta.xml              (TDWG http://rs.tdwg.org/dwc/text/ archive descriptor)"),
        Line::from("  3. eml.xml               (Ecological Metadata Language stub)"),
        Line::from("  4. occurrences.geojson   (RFC-7946 GIS layer for QGIS / Leaflet)"),
        Line::from("  5. spatial_qa_map.svg    (Publication vector distribution & QA map)"),
    ];
    let p = Paragraph::new(lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" DwC-A & Spatial Export "),
    );
    f.render_widget(p, area);
}
