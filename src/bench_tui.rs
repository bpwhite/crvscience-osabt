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

//! Interactive 4-tab terminal workbench (`Ratatui` + `crossterm`) for `osabt bench`:
//! - **Tab 1 (`[1] Primer Studio`):** Live template editor, polymerase preset switcher, TaqMan probe toggle,
//!   ranked primer table, and ASCII hairpin / 3'-dimer bond inspector.
//! - **Tab 2 (`[2] Multiplex & Stock QC`):** Color-coded `N x N` 3'-extensible dimer `Delta G` matrix.
//! - **Tab 3 (`[3] Master Mix & Cycling`):** Interactive `[x]` bench pipetting checklist & `.md` sheet exporter.
//! - **Tab 4 (`[4] Buffers, Resuspend & RCF`):** Live thermal Henderson–Hasselbalch buffer & rotor calculator.

use std::fs;
use std::io;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::{Frame, Terminal};

use crate::tui_util::{osabt_about_status_line, run_tui_session_io};

use crate::bench_calc::{
    build_pcr_master_mix, calc_buffer_recipe, calc_oligo_resuspension, rpm_to_rcf, BufferSystem,
};
use crate::primer_design::{design_pcr_primers, PrimerDesignConfig, PrimerPairCandidate};
use crate::structure::score_dimer;
use crate::thermo::{PolymerasePreset, ReactionConditions};

const DEFAULT_COX1_TEMPLATE: &str =
    "GGTCAACAAATCATAAAGATATTGGACTTTATATTTTATTTTTGGAGCTTGAGCTGGAATAGTAGGAACATCTTTAAGTTTATTAATTCGAGCAGAATTAGGAAACCCAGGATCTTTAATTGGAGATGATCAAATTTATAATGTAATTGTTACAGCTCATGCTTTTATTATAATTTTTTTTATAGTTATACCTATTATAATTGGAGGATTTGGAAATTGATTAGTTCCATTAATATTAGGAGCCCCAGATATAGC";

struct BenchApp {
    active_tab: usize,
    editing_template: bool,
    template_seq: String,
    polymerase: PolymerasePreset,
    target_tm_c: f64,
    include_probe: bool,
    candidates: Vec<PrimerPairCandidate>,
    selected_cand: usize,
    pool_oligos: Vec<(String, String)>,
    rxn_count: usize,
    checked_reagents: Vec<bool>,
    selected_reagent: usize,
    buffer_system: BufferSystem,
    buffer_ph: f64,
    buffer_temp_c: f64,
    rotor_rpm: f64,
    rotor_radius_cm: f64,
    status_msg: String,
}

impl BenchApp {
    fn new() -> Self {
        let mut app = Self {
            active_tab: 0,
            editing_template: false,
            template_seq: DEFAULT_COX1_TEMPLATE.to_string(),
            polymerase: PolymerasePreset::Q5,
            target_tm_c: 60.0,
            include_probe: false,
            candidates: Vec::new(),
            selected_cand: 0,
            pool_oligos: Vec::new(),
            rxn_count: 12,
            checked_reagents: vec![false; 8],
            selected_reagent: 0,
            buffer_system: BufferSystem::Tris,
            buffer_ph: 8.0,
            buffer_temp_c: 25.0,
            rotor_rpm: 12_000.0,
            rotor_radius_cm: 8.5,
            status_msg: String::from(
                "Keys: [1-4/Tab] Switch Tab | [p] Polymerase | [t] qPCR Probe | [+/-] Target Tm | [a] Add to Pool | [q] Quit",
            ),
        };
        app.recompute_primers();
        if let Some(first) = app.candidates.first() {
            app.pool_oligos.push(("COX1_P1_F".to_string(), first.fwd.sequence.clone()));
            app.pool_oligos.push(("COX1_P1_R".to_string(), first.rev.sequence.clone()));
        }
        app
    }

    fn recompute_primers(&mut self) {
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
                if self.selected_cand >= self.candidates.len() {
                    self.selected_cand = 0;
                }
            }
            Err(e) => {
                self.candidates.clear();
                self.status_msg = format!("Design notice: {e}");
            }
        }
    }

    fn cycle_polymerase(&mut self) {
        self.polymerase = match self.polymerase {
            PolymerasePreset::Q5 => PolymerasePreset::Phusion,
            PolymerasePreset::Phusion => PolymerasePreset::Taq,
            PolymerasePreset::Taq => PolymerasePreset::OneTaq,
            _ => PolymerasePreset::Q5,
        };
        self.recompute_primers();
    }

    fn cycle_buffer(&mut self) {
        self.buffer_system = match self.buffer_system {
            BufferSystem::Tris => {
                self.buffer_ph = 7.4;
                BufferSystem::Hepes
            }
            BufferSystem::Hepes => {
                self.buffer_ph = 6.1;
                BufferSystem::Mes
            }
            BufferSystem::Mes => {
                self.buffer_ph = 7.2;
                BufferSystem::SodiumPhosphate
            }
            BufferSystem::SodiumPhosphate => {
                self.buffer_ph = 8.0;
                BufferSystem::Tris
            }
        };
    }

    fn active_ta_and_bp(&self) -> (f64, usize) {
        self.candidates
            .get(self.selected_cand)
            .map(|c| (c.recommended_ta_c, c.amplicon_len_bp))
            .unwrap_or((62.0, 250))
    }

    fn export_markdown_sheet(&mut self) {
        let (ta, bp) = self.active_ta_and_bp();
        match build_pcr_master_mix(
            self.polymerase,
            self.rxn_count,
            25.0,
            0.10,
            1.0,
            ta,
            bp,
            false,
        ) {
            Ok(plan) => {
                let md = plan.to_markdown_sheet("OSABT Bench Export");
                match fs::write("osabt_bench_sheet.md", md) {
                    Ok(_) => {
                        self.status_msg =
                            String::from("Exported bench checklist -> ./osabt_bench_sheet.md");
                    }
                    Err(e) => {
                        self.status_msg = format!("Export error: {e}");
                    }
                }
            }
            Err(e) => {
                self.status_msg = format!("Export error: {e}");
            }
        }
    }
}

/// Launches the interactive full-screen `Ratatui` application for `osabt bench`.
pub fn run_bench_tui() -> io::Result<()> {
    run_tui_session_io(event_loop)
}

fn event_loop(terminal: &mut Terminal<CrosstermBackend<io::Stdout>>) -> io::Result<()> {
    let mut app = BenchApp::new();

    loop {
        terminal.draw(|f| render_ui(f, &app))?;

        if let Event::Key(key) = event::read()? {
            if key.kind != KeyEventKind::Press {
                continue;
            }

            if app.editing_template {
                match key.code {
                    KeyCode::Esc | KeyCode::Enter => {
                        app.editing_template = false;
                        app.recompute_primers();
                        app.status_msg = String::from("Template updated and primers recalculated.");
                    }
                    KeyCode::Backspace => {
                        app.template_seq.pop();
                    }
                    KeyCode::Char(c) => {
                        let up = c.to_ascii_uppercase();
                        if matches!(up, 'A' | 'C' | 'G' | 'T') {
                            app.template_seq.push(up);
                        }
                    }
                    _ => {}
                }
                continue;
            }

            match key.code {
                KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                KeyCode::F(1) | KeyCode::Char('?') => app.status_msg = osabt_about_status_line(),
                KeyCode::Tab => app.active_tab = (app.active_tab + 1) % 4,
                KeyCode::BackTab => app.active_tab = (app.active_tab + 3) % 4,
                KeyCode::Char('1') => app.active_tab = 0,
                KeyCode::Char('2') => app.active_tab = 1,
                KeyCode::Char('3') => app.active_tab = 2,
                KeyCode::Char('4') => app.active_tab = 3,
                KeyCode::Char('p') => app.cycle_polymerase(),
                KeyCode::Char('e') => app.export_markdown_sheet(),
                _ => match app.active_tab {
                    0 => handle_tab1_keys(&mut app, key.code),
                    1 => handle_tab2_keys(&mut app, key.code),
                    2 => handle_tab3_keys(&mut app, key.code),
                    3 => handle_tab4_keys(&mut app, key.code),
                    _ => {}
                },
            }
        }
    }
}

fn handle_tab1_keys(app: &mut BenchApp, code: KeyCode) {
    match code {
        KeyCode::Up | KeyCode::Char('k') => {
            app.selected_cand = app.selected_cand.saturating_sub(1);
        }
        KeyCode::Down | KeyCode::Char('j') => {
            if !app.candidates.is_empty() {
                app.selected_cand = (app.selected_cand + 1).min(app.candidates.len() - 1);
            }
        }
        KeyCode::Char('t') => {
            app.include_probe = !app.include_probe;
            app.recompute_primers();
        }
        KeyCode::Char('+') | KeyCode::Char('=') => {
            app.target_tm_c = (app.target_tm_c + 0.5).min(72.0);
            app.recompute_primers();
        }
        KeyCode::Char('-') => {
            app.target_tm_c = (app.target_tm_c - 0.5).max(48.0);
            app.recompute_primers();
        }
        KeyCode::Char('i') => {
            app.editing_template = true;
            app.status_msg =
                String::from("EDIT MODE: Type A/C/G/T, Backspace to trim, Enter/Esc to apply.");
        }
        KeyCode::Char('a') => {
            if let Some(c) = app.candidates.get(app.selected_cand) {
                let idx = (app.pool_oligos.len() / 2) + 1;
                app.pool_oligos
                    .push((format!("Set{idx}_F"), c.fwd.sequence.clone()));
                app.pool_oligos
                    .push((format!("Set{idx}_R"), c.rev.sequence.clone()));
                app.status_msg = format!("Added Rank #{} pair to Tab 2 Dimer Matrix.", c.rank);
            }
        }
        _ => {}
    }
}

fn handle_tab2_keys(app: &mut BenchApp, code: KeyCode) {
    if code == KeyCode::Char('c') {
        app.pool_oligos.clear();
        app.status_msg = String::from("Cleared multiplex oligo pool.");
    }
}

fn handle_tab3_keys(app: &mut BenchApp, code: KeyCode) {
    match code {
        KeyCode::Left | KeyCode::Char('h') => {
            app.rxn_count = app.rxn_count.saturating_sub(1).max(1);
        }
        KeyCode::Right | KeyCode::Char('l') => {
            app.rxn_count = (app.rxn_count + 1).min(384);
        }
        KeyCode::Up | KeyCode::Char('k') => {
            app.selected_reagent = app.selected_reagent.saturating_sub(1);
        }
        KeyCode::Down | KeyCode::Char('j') => {
            app.selected_reagent = (app.selected_reagent + 1).min(6);
        }
        KeyCode::Char(' ') => {
            if let Some(slot) = app.checked_reagents.get_mut(app.selected_reagent) {
                *slot = !*slot;
            }
        }
        _ => {}
    }
}

fn handle_tab4_keys(app: &mut BenchApp, code: KeyCode) {
    match code {
        KeyCode::Char('b') => app.cycle_buffer(),
        KeyCode::Left | KeyCode::Char('h') => {
            app.buffer_temp_c = (app.buffer_temp_c - 1.0).max(4.0);
        }
        KeyCode::Right | KeyCode::Char('l') => {
            app.buffer_temp_c = (app.buffer_temp_c + 1.0).min(65.0);
        }
        KeyCode::Up | KeyCode::Char('k') => {
            app.rotor_rpm = (app.rotor_rpm + 500.0).min(30_000.0);
        }
        KeyCode::Down | KeyCode::Char('j') => {
            app.rotor_rpm = (app.rotor_rpm - 500.0).max(500.0);
        }
        _ => {}
    }
}

fn render_ui(f: &mut Frame, app: &BenchApp) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(2),
        ])
        .split(f.area());

    let titles = [
        "[1] Primer Studio",
        "[2] Multiplex & Dimer QC",
        "[3] Master Mix & Cycling",
        "[4] Buffers, Resuspend & RCF",
    ]
    .iter()
    .map(|t| Line::from(Span::styled(*t, Style::default().fg(Color::Cyan))))
    .collect::<Vec<_>>();

    let tabs = Tabs::new(titles)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(format!(
                    " OSABT Bench v{} • CRV Science Network (GPL-3.0) ",
                    env!("CARGO_PKG_VERSION")
                )),
        )
        .select(app.active_tab)
        .highlight_style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match app.active_tab {
        0 => render_primer_studio(f, app, chunks[1]),
        1 => render_dimer_matrix(f, app, chunks[1]),
        2 => render_master_mix(f, app, chunks[1]),
        _ => render_bench_calculators(f, app, chunks[1]),
    }

    let footer = Paragraph::new(Line::from(vec![
        Span::styled(" Status: ", Style::default().fg(Color::Green)),
        Span::raw(&app.status_msg),
    ]));
    f.render_widget(footer, chunks[2]);
}

fn render_primer_studio(f: &mut Frame, app: &BenchApp, area: Rect) {
    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(5), Constraint::Min(8)])
        .split(area);

    let seq_title = format!(
        " Template ({} bp) | Polymerase: {:?} [p] | Target Tm: {:.1} C [+/-] | TaqMan Probe: {} [t] | Edit [i] ",
        app.template_seq.len(),
        app.polymerase,
        app.target_tm_c,
        if app.include_probe { "ON" } else { "OFF" }
    );
    let seq_box = Paragraph::new(app.template_seq.as_str())
        .block(Block::default().borders(Borders::ALL).title(seq_title))
        .wrap(Wrap { trim: true });
    f.render_widget(seq_box, rows[0]);

    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(56), Constraint::Percentage(44)])
        .split(rows[1]);

    let header = Row::new(vec!["#", "Span", "Bp", "Fwd Tm", "Rev Tm", "Ta", "3'Dimer"])
        .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD));

    let table_rows = app.candidates.iter().enumerate().map(|(idx, c)| {
        let style = if idx == app.selected_cand {
            Style::default().bg(Color::DarkGray).fg(Color::White)
        } else {
            Style::default()
        };
        Row::new(vec![
            Cell::from(format!("{}", c.rank)),
            Cell::from(format!("{}..{}", c.fwd_start_0 + 1, c.rev_end_0 + 1)),
            Cell::from(format!("{}", c.amplicon_len_bp)),
            Cell::from(format!("{:.1}", c.fwd.tm_c)),
            Cell::from(format!("{:.1}", c.rev.tm_c)),
            Cell::from(format!("{:.1}", c.recommended_ta_c)),
            Cell::from(format!("{:.1}", c.cross_dimer_3p_dg_kcal)),
        ])
        .style(style)
    });

    let widths = [
        Constraint::Length(3),
        Constraint::Length(11),
        Constraint::Length(5),
        Constraint::Length(7),
        Constraint::Length(7),
        Constraint::Length(6),
        Constraint::Length(8),
    ];
    let table = Table::new(table_rows, widths)
        .header(header)
        .block(Block::default().borders(Borders::ALL).title(" Ranked Primer Pairs ([a] Add to Pool) "));
    f.render_widget(table, cols[0]);

    let detail_text = if let Some(c) = app.candidates.get(app.selected_cand) {
        let cond = ReactionConditions::for_preset(app.polymerase);
        let (any_d, d3) = score_dimer(&c.fwd.sequence, &c.rev.sequence, &cond, 37.0)
            .unwrap_or((None, None));
        let dimer_art = d3
            .or(any_d)
            .map(|d| d.ascii_alignment)
            .unwrap_or_else(|| String::from("No significant cross-dimer (< -0.5 kcal/mol)."));
        let prb_line = c
            .qpcr_probe
            .as_ref()
            .map(|p| format!("PRB: 5'-{}-3' (Tm {:.1} C)\n", p.sequence, p.tm_c))
            .unwrap_or_default();
        format!(
            "FWD: 5'-{}-3'\n  GC {:.1}% | 3'dG {:.2} | G4 {:.2}\n\
             REV: 5'-{}-3'\n  GC {:.1}% | 3'dG {:.2} | G4 {:.2}\n\
             {prb_line}\nWorst Fwd x Rev Alignment:\n{dimer_art}",
            c.fwd.sequence,
            c.fwd.gc_pct,
            c.fwd.three_prime_dg_37_kcal,
            c.fwd.g4_hunter_max,
            c.rev.sequence,
            c.rev.gc_pct,
            c.rev.three_prime_dg_37_kcal,
            c.rev.g4_hunter_max
        )
    } else {
        String::from("No valid primer pairs for current constraints.")
    };

    let detail = Paragraph::new(detail_text)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Thermodynamic & 3'-Dimer Inspector "),
        )
        .wrap(Wrap { trim: false });
    f.render_widget(detail, cols[1]);
}

fn render_dimer_matrix(f: &mut Frame, app: &BenchApp, area: Rect) {
    let cond = ReactionConditions::for_preset(app.polymerase);
    let n = app.pool_oligos.len();
    if n == 0 {
        let empty = Paragraph::new(
            "Pool is empty. Press [1] to select a primer pair and [a] to append it to the multiplex matrix.",
        )
        .block(Block::default().borders(Borders::ALL).title(" Multiplex 3'-Dimer Matrix "));
        f.render_widget(empty, area);
        return;
    }

    let mut hdr_cells = vec![Cell::from("Oligo")];
    for (name, _) in &app.pool_oligos {
        hdr_cells.push(Cell::from(name.as_str()));
    }
    let header = Row::new(hdr_cells)
        .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD));

    let rows = app.pool_oligos.iter().map(|(r_name, r_seq)| {
        let mut cells = vec![Cell::from(r_name.as_str())];
        for (_, c_seq) in &app.pool_oligos {
            let dg = score_dimer(r_seq, c_seq, &cond, 37.0)
                .ok()
                .and_then(|(_, d3)| d3.map(|d| d.delta_g_kcal))
                .unwrap_or(0.0);
            let color = if dg < -5.5 {
                Color::Red
            } else if dg < -3.5 {
                Color::Yellow
            } else {
                Color::Green
            };
            cells.push(Cell::from(format!("{dg:>5.1}")).style(Style::default().fg(color)));
        }
        Row::new(cells)
    });

    let mut widths = vec![Constraint::Length(12)];
    for _ in 0..n {
        widths.push(Constraint::Length(9));
    }

    let table = Table::new(rows, widths).header(header).block(
        Block::default().borders(Borders::ALL).title(
            " 3'-Extensible Dimer Delta G_37 (kcal/mol) — Green > -3.5 | Yellow -3.5..-5.5 | Red < -5.5 ([c] Clear) ",
        ),
    );
    f.render_widget(table, area);
}

fn render_master_mix(f: &mut Frame, app: &BenchApp, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(62), Constraint::Percentage(38)])
        .split(area);

    let (ta, bp) = app.active_ta_and_bp();
    let Ok(plan) = build_pcr_master_mix(
        app.polymerase,
        app.rxn_count,
        25.0,
        0.10,
        1.0,
        ta,
        bp,
        false,
    ) else {
        return;
    };

    let header = Row::new(vec!["Done", "Reagent", "1x (uL)", "MasterMix (uL)"])
        .style(Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD));

    let rows = plan.lines.iter().enumerate().map(|(idx, line)| {
        let checked = app.checked_reagents.get(idx).copied().unwrap_or(false);
        let mark = if checked { "[x]" } else { "[ ]" };
        let style = if idx == app.selected_reagent {
            Style::default().bg(Color::DarkGray).fg(Color::White)
        } else {
            Style::default()
        };
        let mm_txt = if line.add_to_master_mix {
            format!("{:.1}", line.master_mix_ul)
        } else {
            String::from("per tube")
        };
        Row::new(vec![
            Cell::from(mark),
            Cell::from(line.reagent.clone()),
            Cell::from(format!("{:.2}", line.per_rxn_ul)),
            Cell::from(mm_txt),
        ])
        .style(style)
    });

    let widths = [
        Constraint::Length(6),
        Constraint::Min(22),
        Constraint::Length(9),
        Constraint::Length(14),
    ];
    let title = format!(
        " {:?} Master Mix: {} Rxns + 10% ([Left/Right] Rxns, [Space] Check, [e] Export .md) ",
        app.polymerase, app.rxn_count
    );
    let table = Table::new(rows, widths)
        .header(header)
        .block(Block::default().borders(Borders::ALL).title(title));
    f.render_widget(table, cols[0]);

    let mut cyc_txt = format!(
        "Dispense {:.2} uL Master Mix + {:.2} uL DNA per tube.\n\nThermocycler Schedule ({bp} bp):\n",
        plan.dispense_mm_per_tube_ul, plan.template_per_tube_ul
    );
    for (step, temp, dur) in &plan.cycling_steps {
        cyc_txt.push_str(&format!("  • {step}: {temp} ({dur})\n"));
    }
    let cyc_box = Paragraph::new(cyc_txt)
        .block(Block::default().borders(Borders::ALL).title(" Cycling Program "))
        .wrap(Wrap { trim: false });
    f.render_widget(cyc_box, cols[1]);
}

fn render_bench_calculators(f: &mut Frame, app: &BenchApp, area: Rect) {
    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(52), Constraint::Percentage(48)])
        .split(area);

    let buf_txt = match calc_buffer_recipe(
        app.buffer_system,
        app.buffer_ph,
        0.05,
        1.0,
        app.buffer_temp_c,
    ) {
        Ok(r) => format!(
            "System [b]: {:?} (50 mM, 1.00 L)\n\
             Target pH:  {:.2} at {:.1} C ([Left/Right] Temp)\n\
             Effective pKa(T): {:.3}\n\n\
             Dry Salt Masses:\n\
               • {}: {:.3} g\n\
               • {}: {:.3} g\n\n\
             Dissolve in 800 mL ddH2O, verify pH at {:.1} C, bring to 1.00 L.",
            r.system,
            r.target_ph,
            r.working_temp_c,
            r.effective_pka,
            r.acid_component_name,
            r.acid_grams,
            r.base_component_name,
            r.base_grams,
            r.working_temp_c
        ),
        Err(e) => format!("Buffer range notice: {e}"),
    };
    let buf_box = Paragraph::new(buf_txt)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Thermal Henderson-Hasselbalch Buffer Solver "),
        )
        .wrap(Wrap { trim: false });
    f.render_widget(buf_box, cols[0]);

    let g_force = rpm_to_rcf(app.rotor_rpm, app.rotor_radius_cm).unwrap_or(0.0);
    let oligo_demo = app
        .candidates
        .get(app.selected_cand)
        .and_then(|c| calc_oligo_resuspension("Selected_FWD", &c.fwd.sequence, 25.0, 100.0, 100.0).ok());

    let oligo_txt = if let Some(p) = oligo_demo {
        format!(
            "Selected FWD Oligo Hydration (25.0 nmol yield):\n\
               • Sequence: 5'-{}-3' ({} nt)\n\
               • Anhydrous MW: {:.1} g/mol ({:.1} ug)\n\
               • Add {:.1} uL TE -> 100 uM Master Stock\n\
               • Mix {:.1} uL Stock + {:.1} uL H2O -> 100 uL @ 10 uM",
            p.sequence,
            p.length_nt,
            p.molecular_weight_g_mol,
            p.mass_ug,
            p.add_buffer_ul,
            p.working_10um_stock_aliquot_ul,
            p.working_10um_water_ul
        )
    } else {
        String::from("Select a primer in Tab 1 for instant hydration math.")
    };

    let right_txt = format!(
        "Centrifuge Rotor ([Up/Down] RPM):\n\
           • Speed:  {:.0} RPM\n\
           • Radius: {:.1} cm\n\
           • Force:  {:.1} x g\n\n\
         {oligo_txt}",
        app.rotor_rpm, app.rotor_radius_cm, g_force
    );
    let right_box = Paragraph::new(right_txt)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Centrifuge RCF & Oligo Hydration "),
        )
        .wrap(Wrap { trim: false });
    f.render_widget(right_box, cols[1]);
}
