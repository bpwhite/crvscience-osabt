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

//! Interactive 4-tab `Ratatui` terminal Information-Theoretic Evolution (`ITME`),
//! Sequence Logo, Epistasis `MIp`, & Kimura Selection Studio (`osabt infoevo`).

use std::fs;
use std::io;
use std::path::PathBuf;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::Terminal;

use crate::tui_util::run_tui_session;

use crate::infoevo_core::{
    compute_alignment_information, compute_clade_divergence, compute_epistasis_mip_matrix,
    demo_ecoli_promoter_and_stemloop_fasta, parse_fasta_alignment,
    simulate_kimura_selection_information, AlignmentComplexityReport, CladeDivergenceReport,
    EpistasisMatrixReport, KimuraSelectionReport, SequenceAlignment,
};
use crate::infoevo_svg::{
    render_epistasis_mip_svg, render_kimura_trajectory_svg, render_sequence_logo_svg,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum InfoEvoTab {
    SequenceLogo,
    ComplexityAndDivergence,
    EpistasisMip,
    KimuraSelection,
}

struct InfoEvoTuiState {
    active_tab: InfoEvoTab,
    aln: SequenceAlignment,
    small_sample_corr: bool,
    info_report: AlignmentComplexityReport,
    epistasis: EpistasisMatrixReport,
    divergence: Option<CladeDivergenceReport>,
    p0_fittest: f64,
    sel_coeff_s: f64,
    kimura_sim: KimuraSelectionReport,
    scroll_offset: usize,
    status_msg: String,
}

impl InfoEvoTuiState {
    fn new(fasta_path: Option<PathBuf>) -> Result<Self, String> {
        let aln = if let Some(p) = fasta_path {
            let raw = fs::read_to_string(&p)
                .map_err(|e| format!("Failed reading FASTA '{}': {e}", p.display()))?;
            let stem = p
                .file_name()
                .and_then(|s| s.to_str())
                .unwrap_or("Alignment");
            parse_fasta_alignment(stem, &raw)?
        } else {
            parse_fasta_alignment(
                "Ecoli_Pribnow_RNA_StemLoop_Demo",
                demo_ecoli_promoter_and_stemloop_fasta(),
            )?
        };

        let small_sample_corr = true;
        let info_report =
            compute_alignment_information(&aln, small_sample_corr, Some((4_600_000.0, 4_500.0)));
        let epistasis = compute_epistasis_mip_matrix(&aln, 0.30, 20);
        let divergence = compute_clade_divergence(&aln, "", "").ok();

        let p0_fittest = 0.01_f64;
        let sel_coeff_s = 0.12_f64;
        let kimura_sim = Self::build_kimura(p0_fittest, sel_coeff_s);

        Ok(Self {
            active_tab: InfoEvoTab::SequenceLogo,
            aln,
            small_sample_corr,
            info_report,
            epistasis,
            divergence,
            p0_fittest,
            sel_coeff_s,
            kimura_sim,
            scroll_offset: 0,
            status_msg: String::from(
                "Keys: [1-4/Tab] Views | [Up/Down] Scroll | [c] Small-Sample e_n | [+/-] Kimura p0* | [[/]] Selection s | [e] Export SVGs | [q] Quit",
            ),
        })
    }

    fn build_kimura(p0_star: f64, s: f64) -> KimuraSelectionReport {
        let p0 = p0_star.clamp(1e-4, 0.50);
        let rem = 1.0 - p0;
        simulate_kimura_selection_information(
            &[p0, 0.35 * rem, 0.65 * rem],
            &[s, 0.35 * s, 0.0],
            1.0,
            240,
        )
        .expect("valid Kimura parameters")
    }

    fn recompute_info(&mut self) {
        self.info_report = compute_alignment_information(
            &self.aln,
            self.small_sample_corr,
            Some((4_600_000.0, 4_500.0)),
        );
    }

    fn export_all_svgs(&mut self) {
        let logo = render_sequence_logo_svg(&self.info_report, &self.aln.title);
        let arcs = render_epistasis_mip_svg(&self.epistasis, &self.info_report, &self.aln.title);
        let kim = render_kimura_trajectory_svg(&self.kimura_sim, "Kimura Adaptive Sweep");

        let r1 = fs::write("infoevo_sequence_logo.svg", logo);
        let r2 = fs::write("infoevo_epistasis_mip.svg", arcs);
        let r3 = fs::write("infoevo_kimura_selection.svg", kim);
        if r1.is_ok() && r2.is_ok() && r3.is_ok() {
            self.status_msg = String::from(
                "Exported infoevo_sequence_logo.svg, infoevo_epistasis_mip.svg, and infoevo_kimura_selection.svg!",
            );
        } else {
            self.status_msg = String::from("SVG export failed: check write permissions.");
        }
    }
}

/// Launches the interactive 4-tab `Ratatui` terminal studio for `osabt infoevo`.
pub fn run_infoevo_tui(fasta_path: Option<PathBuf>) -> Result<(), String> {
    let mut state = InfoEvoTuiState::new(fasta_path)?;
    run_tui_session(|terminal| run_event_loop(terminal, &mut state))
}

fn run_event_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut InfoEvoTuiState,
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
                        state.active_tab = InfoEvoTab::SequenceLogo;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('2') => {
                        state.active_tab = InfoEvoTab::ComplexityAndDivergence;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('3') => {
                        state.active_tab = InfoEvoTab::EpistasisMip;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('4') => {
                        state.active_tab = InfoEvoTab::KimuraSelection;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            InfoEvoTab::SequenceLogo => InfoEvoTab::ComplexityAndDivergence,
                            InfoEvoTab::ComplexityAndDivergence => InfoEvoTab::EpistasisMip,
                            InfoEvoTab::EpistasisMip => InfoEvoTab::KimuraSelection,
                            InfoEvoTab::KimuraSelection => InfoEvoTab::SequenceLogo,
                        };
                        state.scroll_offset = 0;
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        state.scroll_offset = state.scroll_offset.saturating_add(1);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.scroll_offset = state.scroll_offset.saturating_sub(1);
                    }
                    KeyCode::Char('c') => {
                        state.small_sample_corr = !state.small_sample_corr;
                        state.recompute_info();
                        state.status_msg = format!(
                            "Small-sample Basharin/Schneider correction e_n: {}",
                            if state.small_sample_corr { "ON" } else { "OFF" }
                        );
                    }
                    KeyCode::Char('+') | KeyCode::Char('=') => {
                        state.p0_fittest = (state.p0_fittest * 1.5).min(0.45);
                        state.kimura_sim =
                            InfoEvoTuiState::build_kimura(state.p0_fittest, state.sel_coeff_s);
                    }
                    KeyCode::Char('-') => {
                        state.p0_fittest = (state.p0_fittest / 1.5).max(0.0005);
                        state.kimura_sim =
                            InfoEvoTuiState::build_kimura(state.p0_fittest, state.sel_coeff_s);
                    }
                    KeyCode::Char(']') => {
                        state.sel_coeff_s = (state.sel_coeff_s + 0.02).min(0.60);
                        state.kimura_sim =
                            InfoEvoTuiState::build_kimura(state.p0_fittest, state.sel_coeff_s);
                    }
                    KeyCode::Char('[') => {
                        state.sel_coeff_s = (state.sel_coeff_s - 0.02).max(0.02);
                        state.kimura_sim =
                            InfoEvoTuiState::build_kimura(state.p0_fittest, state.sel_coeff_s);
                    }
                    KeyCode::Char('e') => state.export_all_svgs(),
                    _ => {}
                }
            }
        }
    }
}

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &InfoEvoTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Sequence Logo & R_seq",
        "[2] Adami Complexity & Clade JSD",
        "[3] Epistasis MIp Coevolution",
        "[4] Kimura Selection Regret (ITME)",
    ]
    .into_iter()
    .map(Line::from)
    .collect::<Vec<_>>();

    let tab_idx = match state.active_tab {
        InfoEvoTab::SequenceLogo => 0,
        InfoEvoTab::ComplexityAndDivergence => 1,
        InfoEvoTab::EpistasisMip => 2,
        InfoEvoTab::KimuraSelection => 3,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT InfoEvo (v1.1.0) — {} ({} seqs × {} cols) ",
            state.aln.title, state.info_report.num_sequences, state.info_report.alignment_length
        )))
        .select(tab_idx)
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        InfoEvoTab::SequenceLogo => draw_logo_tab(f, chunks[1], state),
        InfoEvoTab::ComplexityAndDivergence => draw_complexity_tab(f, chunks[1], state),
        InfoEvoTab::EpistasisMip => draw_epistasis_tab(f, chunks[1], state),
        InfoEvoTab::KimuraSelection => draw_kimura_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Status & Controls "))
        .style(Style::default().fg(Color::Gray));
    f.render_widget(footer, chunks[2]);
}

fn sym_ratatui_color(ch: char) -> Color {
    match ch.to_ascii_uppercase() {
        'A' => Color::Green,
        'C' => Color::Blue,
        'G' => Color::Yellow,
        'T' | 'U' => Color::Red,
        _ => Color::Cyan,
    }
}

fn draw_logo_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &InfoEvoTuiState) {
    let split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(9), Constraint::Min(6)])
        .split(area);

    // Unicode terminal sequence logo across up to 64 columns
    let max_cols = state.info_report.sites.len().min(60);
    let max_b = state.info_report.max_bits_per_site.max(2.0);
    let mut logo_lines: Vec<Line<'static>> = Vec::with_capacity(7);

    for row_level in (0..5).rev() {
        let thresh = ((row_level as f64) + 0.35) / 5.0 * max_b;
        let mut spans = vec![Span::styled(
            format!("{:>4.1}b │ ", ((row_level + 1) as f64 / 5.0) * max_b),
            Style::default().fg(Color::DarkGray),
        )];
        for s in state.info_report.sites.iter().take(max_cols) {
            if s.r_sequence_bits >= thresh {
                spans.push(Span::styled(
                    format!("{}", s.consensus_symbol),
                    Style::default()
                        .fg(sym_ratatui_color(s.consensus_symbol))
                        .add_modifier(Modifier::BOLD),
                ));
            } else {
                spans.push(Span::raw(" "));
            }
        }
        logo_lines.push(Line::from(spans));
    }
    let mut ruler = String::from("  pos │ ");
    for i in 1..=max_cols {
        ruler.push(char::from_digit((i % 10) as u32, 10).unwrap_or(' '));
    }
    logo_lines.push(Line::from(Span::styled(
        ruler,
        Style::default().fg(Color::Gray),
    )));

    let logo_box = Paragraph::new(logo_lines).block(Block::default().borders(Borders::ALL).title(
        format!(
            " Schneider & Stephens (1990) Terminal Logo — Total R_seq = {:.2} bits (e_n {}) ",
            state.info_report.total_r_sequence_bits,
            if state.small_sample_corr { "ON" } else { "OFF" }
        ),
    ));
    f.render_widget(logo_box, split[0]);

    let header = Row::new([
        "Pos",
        "Consensus",
        "Non-Gap N",
        "H(l) (bits)",
        "e_n (bits)",
        "R_seq (bits)",
        "Stack Breakdown (bits)",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state
        .info_report
        .sites
        .iter()
        .skip(state.scroll_offset)
        .map(|s| {
            let stack_str = s
                .letter_heights_bits
                .iter()
                .rev()
                .map(|&(ch, hb)| format!("{ch}:{hb:.2}"))
                .collect::<Vec<_>>()
                .join(" ");
            Row::new([
                Cell::from(format!("{}", s.position_1based)),
                Cell::from(format!("{}", s.consensus_symbol))
                    .style(Style::default().fg(sym_ratatui_color(s.consensus_symbol))),
                Cell::from(format!("{}", s.non_gap_count)),
                Cell::from(format!("{:.3}", s.shannon_bits)),
                Cell::from(format!("{:.3}", s.small_sample_error_bits)),
                Cell::from(format!("{:.3}", s.r_sequence_bits)),
                Cell::from(stack_str),
            ])
        });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(6),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(12),
            Constraint::Length(11),
            Constraint::Length(13),
            Constraint::Min(24),
        ],
    )
    .header(header)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Per-Site Positional Entropy & Information Profile "),
    );
    f.render_widget(tbl, split[1]);
}

fn draw_complexity_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &InfoEvoTuiState) {
    let split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(44), Constraint::Percentage(56)])
        .split(area);

    let r_freq_txt = state
        .info_report
        .r_frequency_target_bits
        .map(|v| format!("{v:.2} bits (E. coli ~4.6 Mb / 4,500 sites)"))
        .unwrap_or_else(|| String::from("N/A"));

    let mut summary_lines = vec![
        Line::from(Span::styled(
            "Adami et al. (2000) Genomic Physical Complexity",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        )),
        Line::from(""),
        Line::from(format!(
            "Alphabet                 : {}",
            state.info_report.alphabet.label()
        )),
        Line::from(format!(
            "Max Locus Capacity       : {:.2} bits",
            (state.info_report.alignment_length as f64) * state.info_report.max_bits_per_site
        )),
        Line::from(format!(
            "Adami Complexity (C)     : {:.3} bits",
            state.info_report.adami_physical_complexity_bits
        )),
        Line::from(format!(
            "Schneider R_sequence     : {:.3} bits",
            state.info_report.total_r_sequence_bits
        )),
        Line::from(format!(
            "R_frequency Benchmark    : {}",
            r_freq_txt
        )),
        Line::from(format!(
            "Complexity Density       : {:.1}%",
            state.info_report.complexity_density * 100.0
        )),
        Line::from(format!(
            "Effective Invariant Pos  : {:.2} / {} sites",
            state.info_report.effective_conserved_sites, state.info_report.alignment_length
        )),
    ];

    if let Some(div) = &state.divergence {
        summary_lines.push(Line::from(""));
        summary_lines.push(Line::from(Span::styled(
            "Two-Clade Informational Divergence",
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        )));
        summary_lines.push(Line::from(format!(
            "Clades Compared          : {} (n={}) vs {} (n={})",
            div.clade_a_name, div.clade_a_n, div.clade_b_name, div.clade_b_n
        )));
        summary_lines.push(Line::from(format!(
            "Total Jensen–Shannon D_JS: {:.3} bits (mean {:.3}/site)",
            div.total_js_bits, div.mean_js_bits
        )));
        summary_lines.push(Line::from(format!(
            "Total KL D_KL(A || B)    : {:.3} bits",
            div.total_kl_bits
        )));
    }

    let left_card = Paragraph::new(summary_lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Physical Complexity & Maxwell's Demon Summary "),
        );
    f.render_widget(left_card, split[0]);

    if let Some(div) = &state.divergence {
        let hdr = Row::new(["Pos", &div.clade_a_name, &div.clade_b_name, "D_JS (bits)", "D_KL (bits)"])
            .style(
                Style::default()
                    .fg(Color::Yellow)
                    .add_modifier(Modifier::BOLD),
            );
        let rows = div.sites.iter().skip(state.scroll_offset).map(|s| {
            let st = if s.js_divergence_bits >= 0.50 {
                Style::default().fg(Color::LightRed).add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };
            Row::new([
                format!("{}", s.position_1based),
                format!("{}", s.clade_a_consensus),
                format!("{}", s.clade_b_consensus),
                format!("{:.4}", s.js_divergence_bits),
                format!("{:.4}", s.kl_p_from_q_bits),
            ])
            .style(st)
        });
        let tbl = Table::new(
            rows,
            [
                Constraint::Length(6),
                Constraint::Length(12),
                Constraint::Length(12),
                Constraint::Length(13),
                Constraint::Length(13),
            ],
        )
        .header(hdr)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Per-Site Clade Jensen–Shannon & KL Divergence "),
        );
        f.render_widget(tbl, split[1]);
    }
}

fn draw_epistasis_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &InfoEvoTuiState) {
    let hdr = Row::new([
        "Rank",
        "Pair (i–j)",
        "Raw MI (b)",
        "APC (b)",
        "MIp (bits)",
        "MIr",
        "Z-Score",
        "Top Covarying States",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state
        .epistasis
        .top_pairs
        .iter()
        .enumerate()
        .skip(state.scroll_offset)
        .map(|(idx, p)| {
            let st = if p.z_score >= 3.0 {
                Style::default()
                    .fg(Color::LightGreen)
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };
            Row::new([
                format!("#{}", idx + 1),
                format!("{} ↔ {}", p.pos_i_1based, p.pos_j_1based),
                format!("{:.3}", p.mi_bits),
                format!("{:.3}", p.apc_bits),
                format!("{:.3}", p.mip_bits),
                format!("{:.3}", p.normalized_mir),
                format!("{:+.2}", p.z_score),
                p.top_covarying_states.clone(),
            ])
            .style(st)
        });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(6),
            Constraint::Length(11),
            Constraint::Length(11),
            Constraint::Length(9),
            Constraint::Length(11),
            Constraint::Length(7),
            Constraint::Length(9),
            Constraint::Min(28),
        ],
    )
    .header(hdr)
    .block(Block::default().borders(Borders::ALL).title(format!(
        " Dunn, Wahl & Gloor (2008) APC-Corrected Epistasis Ranking (Mean MI = {:.3} bits) ",
        state.epistasis.mean_raw_mi_bits
    )));
    f.render_widget(tbl, area);
}

fn draw_kimura_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &InfoEvoTuiState) {
    let split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(42), Constraint::Percentage(58)])
        .split(area);

    let sim = &state.kimura_sim;
    let info_txt = vec![
        Line::from(Span::styled(
            "Kimura (1961) Genetic Information & Regret",
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        )),
        Line::from(""),
        Line::from(format!("Initial Fittest Freq p0* : {:.5}", sim.initial_fittest_freq)),
        Line::from(format!("Selection Advantage s    : {:.3}", state.sel_coeff_s)),
        Line::from(format!(
            "Kimura Limit -log2(p0*)  : {:.4} bits",
            sim.kimura_limit_bits
        )),
        Line::from(format!(
            "Final KL Gain I_T        : {:.4} bits",
            sim.final_info_gain_bits
        )),
        Line::from(format!(
            "Cumulative Load L_sub,T  : {:.4} bits",
            sim.final_substitutional_load_bits
        )),
        Line::from(format!(
            "Conversion Efficiency    : {:.2}%",
            sim.final_efficiency * 100.0
        )),
        Line::from(""),
        Line::from(Span::styled(
            "Theorem Check: I_t <= L_sub,t holds at every generation, and lim(I_T) = lim(L_sub,T) = -log2(p0*).",
            Style::default().fg(Color::Green),
        )),
    ];
    let left_p = Paragraph::new(info_txt)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Parameters (+/- p0*, [/] s) "),
        );
    f.render_widget(left_p, split[0]);

    let hdr = Row::new([
        "Gen t",
        "p_t*",
        "H(p_t) (b)",
        "I_t Gain (b)",
        "L_sub,t (b)",
        "Eff (%)",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );
    let rows = sim
        .steps
        .iter()
        .step_by(10)
        .map(|st| {
            Row::new([
                format!("{}", st.generation),
                format!("{:.4}", st.fittest_freq),
                format!("{:.3}", st.population_entropy_bits),
                format!("{:.3}", st.kl_info_gain_bits),
                format!("{:.3}", st.cumulative_substitutional_load_bits),
                format!("{:.1}%", st.selection_efficiency * 100.0),
            ])
        });
    let tbl = Table::new(
        rows,
        [
            Constraint::Length(7),
            Constraint::Length(8),
            Constraint::Length(11),
            Constraint::Length(13),
            Constraint::Length(12),
            Constraint::Length(9),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Selection Trajectory (Every 10 Generations) "),
    );
    f.render_widget(tbl, split[1]);
}
