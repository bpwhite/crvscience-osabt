// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Interactive 4-tab `Ratatui` terminal 96/384-Well Assay Planner & 4PL Dose-Response Studio (`osabt plate`).

use std::fs;
use std::io;
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

use crate::plate_curve::{
    demo_measured_96well_plate, evaluate_plate_qc_and_normalize, fit_all_plate_compounds,
    render_4pl_curves_svg,
};
use crate::plate_layout::{
    export_bench_pipetting_markdown, export_plate_pipetting_csv, generate_randomized_plate,
    render_plate_map_svg, Microplate, PlateFormat, PlatePlanConfig, WellRole,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum PlateTab {
    GridMap,
    PipettingSheet,
    AssayQc,
    DoseResponse4Pl,
}

struct PlateTuiState {
    active_tab: PlateTab,
    cfg: PlatePlanConfig,
    plate: Microplate,
    status_msg: String,
}

impl PlateTuiState {
    fn new() -> Self {
        let cfg = PlatePlanConfig::default();
        let plate = demo_measured_96well_plate();
        Self {
            active_tab: PlateTab::GridMap,
            cfg,
            plate,
            status_msg: String::from(
                "1-4/Tab: switch view | [r]: re-seed randomization | [f]: toggle 96/384-well | [e]: export CSV + Markdown + Plate/4PL SVGs | [q]: quit",
            ),
        }
    }

    fn rerandomize(&mut self, new_seed: u64) {
        self.cfg.rng_seed = new_seed;
        if let Ok(mut p) = generate_randomized_plate(&self.cfg) {
            // Populate simulated raw signals so QC & 4PL tabs update live on every re-seed!
            populate_simulated_signals(&mut p, new_seed);
            let _ = evaluate_plate_qc_and_normalize(&mut p);
            self.plate = p;
            self.status_msg = format!(
                "Regenerated {} randomized layout with PRNG seed {}.",
                self.plate.format.label(),
                new_seed
            );
        }
    }

    fn export_all(&mut self) {
        let fits = fit_all_plate_compounds(&self.plate);
        let r1 = fs::write("osabt_plate_layout.csv", export_plate_pipetting_csv(&self.plate));
        let r2 = fs::write(
            "osabt_bench_pipetting.md",
            export_bench_pipetting_markdown(&self.plate, &self.cfg),
        );
        let r3 = fs::write("osabt_plate_heatmap.svg", render_plate_map_svg(&self.plate));
        let r4 = fs::write(
            "osabt_4pl_curves.svg",
            render_4pl_curves_svg(&self.plate.name, &fits),
        );

        if r1.is_ok() && r2.is_ok() && r3.is_ok() && r4.is_ok() {
            self.status_msg = String::from(
                "Exported ./osabt_plate_layout.csv, ./osabt_bench_pipetting.md, ./osabt_plate_heatmap.svg, and ./osabt_4pl_curves.svg",
            );
        } else {
            self.status_msg = String::from("Failed writing one or more plate export files.");
        }
    }
}

fn populate_simulated_signals(plate: &mut Microplate, seed: u64) {
    let mut st = seed ^ 0xA5A5_5A5A_F0F0_0F0F;
    let mut norm_eps = || -> f64 {
        let mut s = 0.0_f64;
        for _ in 0..6 {
            st = st.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
            s += ((st >> 32) as f64) / (u32::MAX as f64);
        }
        (s - 3.0) * 1.414
    };

    let bg = 1_200.0;
    let sig = 52_000.0;
    for w in &mut plate.wells {
        match w.role {
            WellRole::NegativeControl => w.raw_readout = Some(sig + norm_eps() * 1_200.0),
            WellRole::PositiveControl => w.raw_readout = Some(bg + norm_eps() * 90.0),
            WellRole::CompoundDose => {
                let c = w.concentration_nm.unwrap_or(1.0);
                let (ic50, hill) = match w.compound_name.as_str() {
                    "Staurosporine" => (4.2_f64, -1.04_f64),
                    "CRV-104" => (18.5_f64, -1.15_f64),
                    _ => (420.0_f64, -0.96_f64),
                };
                let u = 10_f64.powf((ic50.log10() - c.log10()) * hill);
                let pct = 2.0 + (99.5 - 2.0) / (1.0 + u);
                w.raw_readout = Some((bg + (pct / 100.0) * (sig - bg) + norm_eps() * 900.0).max(80.0));
            }
            _ => {}
        }
    }
}

/// Launches the interactive 4-tab `Ratatui` Microplate & 4PL Dose-Response Studio.
pub fn run_plate_tui() -> Result<(), String> {
    let mut state = PlateTuiState::new();

    enable_raw_mode().map_err(|e| format!("Failed enabling terminal raw mode: {e}"))?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)
        .map_err(|e| format!("Failed entering alternate screen: {e}"))?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal =
        Terminal::new(backend).map_err(|e| format!("Failed creating Ratatui terminal: {e}"))?;

    let res = plate_loop(&mut terminal, &mut state);

    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen);
    let _ = terminal.show_cursor();
    res
}

fn plate_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut PlateTuiState,
) -> Result<(), String> {
    loop {
        terminal
            .draw(|f| draw_plate_ui(f, state))
            .map_err(|e| format!("Draw error: {e}"))?;

        if event::poll(Duration::from_millis(120)).map_err(|e| e.to_string())? {
            if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }
                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                    KeyCode::Char('1') => state.active_tab = PlateTab::GridMap,
                    KeyCode::Char('2') => state.active_tab = PlateTab::PipettingSheet,
                    KeyCode::Char('3') => state.active_tab = PlateTab::AssayQc,
                    KeyCode::Char('4') => state.active_tab = PlateTab::DoseResponse4Pl,
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            PlateTab::GridMap => PlateTab::PipettingSheet,
                            PlateTab::PipettingSheet => PlateTab::AssayQc,
                            PlateTab::AssayQc => PlateTab::DoseResponse4Pl,
                            PlateTab::DoseResponse4Pl => PlateTab::GridMap,
                        };
                    }
                    KeyCode::Char('r') => {
                        let next_seed = state.cfg.rng_seed.wrapping_add(7919);
                        state.rerandomize(next_seed);
                    }
                    KeyCode::Char('f') => {
                        state.cfg.format = match state.cfg.format {
                            PlateFormat::Well96 => {
                                state.cfg.assay_volume_ul = 40.0;
                                PlateFormat::Well384
                            }
                            PlateFormat::Well384 => {
                                state.cfg.assay_volume_ul = 100.0;
                                PlateFormat::Well96
                            }
                        };
                        let s = state.cfg.rng_seed;
                        state.rerandomize(s);
                    }
                    KeyCode::Char('e') => state.export_all(),
                    _ => {}
                }
            }
        }
    }
}

fn draw_plate_ui(f: &mut ratatui::Frame<'_>, state: &PlateTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Randomized Plate Map",
        "[2] Bench Pipetting Sheet",
        "[3] Zhang Z'-Factor QC",
        "[4] 4PL IC50 / EC50 Curves",
    ]
    .iter()
    .map(|t| Line::from(*t))
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        PlateTab::GridMap => 0,
        PlateTab::PipettingSheet => 1,
        PlateTab::AssayQc => 2,
        PlateTab::DoseResponse4Pl => 3,
    };

    let tabs = Tabs::new(titles)
        .select(idx)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT 96/384-Well Assay & 4PL Studio • {} (Seed {}) ",
            state.plate.format.label(),
            state.cfg.rng_seed
        )))
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        PlateTab::GridMap => draw_grid_tab(f, chunks[1], state),
        PlateTab::PipettingSheet => draw_pipette_tab(f, chunks[1], state),
        PlateTab::AssayQc => draw_qc_tab(f, chunks[1], state),
        PlateTab::DoseResponse4Pl => draw_4pl_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls "));
    f.render_widget(footer, chunks[2]);
}

fn draw_grid_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &PlateTuiState) {
    let (nrows, ncols) = state.plate.format.dims();
    let show_cols = ncols.min(12);
    let mut header_cells = vec![Cell::from("Row")];
    for c in 1..=show_cols {
        header_cells.push(Cell::from(format!("{c:02}")));
    }
    let header = Row::new(header_cells).style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let mut rows = Vec::with_capacity(nrows);
    for r in 0..nrows {
        let mut cells = vec![Cell::from(format!(" {} ", (b'A' + r as u8) as char))];
        for c in 0..show_cols {
            if let Some(w) = state.plate.well_at(r, c) {
                let (tag, col) = match w.role {
                    WellRole::BufferMoat => ("MOAT".to_string(), Color::Blue),
                    WellRole::NegativeControl => ("NEG-".to_string(), Color::Green),
                    WellRole::PositiveControl => ("POS+".to_string(), Color::Red),
                    WellRole::CompoundDose => {
                        let pfx: String = w.compound_name.chars().take(4).collect();
                        (pfx, Color::Cyan)
                    }
                    WellRole::Empty => ("----".to_string(), Color::DarkGray),
                };
                cells.push(Cell::from(Span::styled(tag, Style::default().fg(col))));
            }
        }
        rows.push(Row::new(cells));
    }

    let mut widths = vec![Constraint::Length(4)];
    for _ in 0..show_cols {
        widths.push(Constraint::Length(6));
    }

    let tbl = Table::new(rows, widths).header(header).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Edge-Buffered Randomized Well Grid (MOAT = PBS border | NEG- = Vehicle | POS+ = 100% Inhib) "),
    );
    f.render_widget(tbl, area);
}

fn draw_pipette_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &PlateTuiState) {
    let md = export_bench_pipetting_markdown(&state.plate, &state.cfg);
    let lines: Vec<Line<'_>> = md.lines().take(28).map(Line::from).collect();
    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Printable Bench Serial-Dilution & Destination Well Sheet "),
        );
    f.render_widget(p, area);
}

fn draw_qc_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &PlateTuiState) {
    let mut clone_p = state.plate.clone();
    let mut lines = Vec::new();
    if let Ok(qc) = evaluate_plate_qc_and_normalize(&mut clone_p) {
        let grade_col = if qc.z_prime_factor >= 0.5 {
            Color::Green
        } else if qc.z_prime_factor >= 0.0 {
            Color::Yellow
        } else {
            Color::Red
        };
        lines.push(Line::from(Span::styled(
            format!("Assay Quality Assessment: {}", qc.quality_grade),
            Style::default().fg(grade_col).add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(""));
        lines.push(Line::from(format!(
            "• Zhang Z'-Factor                     : {:.4}",
            qc.z_prime_factor
        )));
        lines.push(Line::from(format!(
            "• Strictly Standardized Mean Diff     : {:.2}",
            qc.ssmd
        )));
        lines.push(Line::from(format!(
            "• Signal-to-Background (S/B)          : {:.1}-fold",
            qc.signal_to_background
        )));
        lines.push(Line::from(format!(
            "• Negative Control (n={})            : {:.1} +- {:.1} (CV = {:.2}%)",
            qc.neg_ctrl_n, qc.neg_mean, qc.neg_sd, qc.neg_cv_pct
        )));
        lines.push(Line::from(format!(
            "• Positive Control (n={})            : {:.1} +- {:.1} (CV = {:.2}%)",
            qc.pos_ctrl_n, qc.pos_mean, qc.pos_sd, qc.pos_cv_pct
        )));
    }
    let p = Paragraph::new(lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Zhang, Chung & Oldenburg (1999) High-Throughput Screening QC "),
    );
    f.render_widget(p, area);
}

fn draw_4pl_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &PlateTuiState) {
    let fits = fit_all_plate_compounds(&state.plate);
    let header = Row::new([
        "Compound",
        "Rel IC50 (nM)",
        "95% CI (nM)",
        "Abs IC50 (nM)",
        "IC90 (nM)",
        "Hill Slope",
        "Bottom/Top %",
        "R²",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let rows = fits.iter().map(|fit| {
        let abs_s = fit
            .absolute_ic50_nm
            .map(|v| format!("{v:.2}"))
            .unwrap_or_else(|| String::from("NA"));
        Row::new(vec![
            Cell::from(fit.compound_name.clone()),
            Cell::from(format!("{:.2}", fit.relative_ec50_nm)),
            Cell::from(format!(
                "[{:.2}, {:.2}]",
                fit.relative_ec50_ci95_nm.0, fit.relative_ec50_ci95_nm.1
            )),
            Cell::from(abs_s),
            Cell::from(format!("{:.1}", fit.ec90_nm)),
            Cell::from(format!("{:+.2} +- {:.2}", fit.hill_slope, fit.hill_se)),
            Cell::from(format!("{:.1}% / {:.1}%", fit.bottom, fit.top)),
            Cell::from(format!("{:.4}", fit.r_squared)),
        ])
    });

    let widths = [
        Constraint::Length(15),
        Constraint::Length(14),
        Constraint::Length(20),
        Constraint::Length(14),
        Constraint::Length(12),
        Constraint::Length(15),
        Constraint::Length(15),
        Constraint::Length(8),
    ];

    let tbl = Table::new(rows, widths).header(header).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Levenberg-Marquardt 4-Parameter Logistic (4PL) Nonlinear Regression "),
    );
    f.render_widget(tbl, area);
}
