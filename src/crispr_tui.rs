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

//! Interactive 4-tab `Ratatui` terminal CRISPR-Cas9 / Cas12a / Base-Editor Guide Designer
//! & Doench `CFD` Off-Target Scorer (`osabt crispr`).

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

use crate::crispr_core::{
    demo_crispr_target_and_offtargets, design_golden_gate_crispr_oligos, scan_crispr_candidates,
    CrisprCloningVector, CrisprGuideCandidate, CrisprNuclease,
};
use crate::crispr_offtarget::{
    export_crispr_design_tsv, render_crispr_locus_svg, scan_offtargets_and_specificity,
    CrisprOffTargetHit,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum CrisprTab {
    RankedGuides,
    OffTargetCfdInspector,
    BaseEditorWindows,
    GoldenGateOligos,
}

struct CrisprTuiState {
    active_tab: CrisprTab,
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

impl CrisprTuiState {
    fn new(target_path: Option<PathBuf>, bg_path: Option<PathBuf>) -> Result<Self, String> {
        let (demo_t, demo_bg) = demo_crispr_target_and_offtargets();
        let target_raw = if let Some(tp) = target_path {
            fs::read_to_string(&tp)
                .map_err(|e| format!("Failed reading target FASTA '{}': {e}", tp.display()))?
        } else {
            demo_t
        };
        let bg_raw = if let Some(bp) = bg_path {
            fs::read_to_string(&bp)
                .map_err(|e| format!("Failed reading background FASTA '{}': {e}", bp.display()))?
        } else {
            demo_bg
        };

        let mut state = Self {
            active_tab: CrisprTab::RankedGuides,
            target_raw,
            bg_raw,
            target_len_bp: 0,
            nuclease: CrisprNuclease::SpCas9,
            vector: CrisprCloningVector::Px459BbsI,
            guides: Vec::new(),
            offtargets: Vec::new(),
            selected_idx: 0,
            status_msg: String::from(
                "Keys: [1-4/Tab] Views | [Up/Down] Select Guide | [n] Cycle Nuclease (SpCas9/SaCas9/LbCas12a) | [v] Vector | [e] Export TSV+SVG | [q] Quit",
            ),
        };
        state.rescan()?;
        Ok(state)
    }

    fn rescan(&mut self) -> Result<(), String> {
        let (norm_seq, mut guides) = scan_crispr_candidates(&self.target_raw, self.nuclease)?;
        let mut hits = scan_offtargets_and_specificity(&mut guides, &self.bg_raw, 4)?;

        // Pair-sort guides + their off-target hit vectors by descending composite_score
        let mut paired: Vec<(CrisprGuideCandidate, Vec<CrisprOffTargetHit>)> =
            guides.into_iter().zip(hits.drain(..)).collect();
        paired.sort_by(|a, b| {
            b.0.composite_score
                .partial_cmp(&a.0.composite_score)
                .unwrap_or(std::cmp::Ordering::Equal)
        });
        for (i, (g, _)) in paired.iter_mut().enumerate() {
            g.guide_id = format!("gRNA_{:02}", i + 1);
        }

        let (g_sorted, h_sorted): (Vec<_>, Vec<_>) = paired.into_iter().unzip();
        self.target_len_bp = norm_seq.len();
        self.guides = g_sorted;
        self.offtargets = h_sorted;
        self.selected_idx = 0;
        Ok(())
    }
}

/// Launches the interactive 4-tab `Ratatui` terminal studio for `osabt crispr`.
pub fn run_crispr_tui(
    target_path: Option<PathBuf>,
    bg_path: Option<PathBuf>,
) -> Result<(), String> {
    let mut state = CrisprTuiState::new(target_path, bg_path)?;
    run_tui_session(|terminal| run_loop(terminal, &mut state))
}

fn run_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut CrisprTuiState,
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
                    KeyCode::Char('1') => state.active_tab = CrisprTab::RankedGuides,
                    KeyCode::Char('2') => state.active_tab = CrisprTab::OffTargetCfdInspector,
                    KeyCode::Char('3') => state.active_tab = CrisprTab::BaseEditorWindows,
                    KeyCode::Char('4') => state.active_tab = CrisprTab::GoldenGateOligos,
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            CrisprTab::RankedGuides => CrisprTab::OffTargetCfdInspector,
                            CrisprTab::OffTargetCfdInspector => CrisprTab::BaseEditorWindows,
                            CrisprTab::BaseEditorWindows => CrisprTab::GoldenGateOligos,
                            CrisprTab::GoldenGateOligos => CrisprTab::RankedGuides,
                        };
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        if !state.guides.is_empty() {
                            state.selected_idx =
                                (state.selected_idx + 1).min(state.guides.len() - 1);
                        }
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.selected_idx = state.selected_idx.saturating_sub(1);
                    }
                    KeyCode::Char('n') => {
                        state.nuclease = match state.nuclease {
                            CrisprNuclease::SpCas9 => CrisprNuclease::SaCas9,
                            CrisprNuclease::SaCas9 => CrisprNuclease::LbCas12a,
                            CrisprNuclease::LbCas12a => CrisprNuclease::SpCas9,
                        };
                        state.vector = match state.nuclease {
                            CrisprNuclease::LbCas12a => CrisprCloningVector::LbCas12aBsaI,
                            _ => CrisprCloningVector::Px459BbsI,
                        };
                        let _ = state.rescan();
                        state.status_msg = format!(
                            "Switched nuclease to {} -> {} candidate guides.",
                            state.nuclease.label(),
                            state.guides.len()
                        );
                    }
                    KeyCode::Char('v') => {
                        state.vector = match state.vector {
                            CrisprCloningVector::Px459BbsI => {
                                CrisprCloningVector::LentiCrisprV2BsmBI
                            }
                            CrisprCloningVector::LentiCrisprV2BsmBI => {
                                CrisprCloningVector::LbCas12aBsaI
                            }
                            CrisprCloningVector::LbCas12aBsaI => CrisprCloningVector::Px459BbsI,
                        };
                        state.status_msg =
                            format!("Cloning overhang vector: {}", state.vector.label());
                    }
                    KeyCode::Char('e') => {
                        let tsv = export_crispr_design_tsv(&state.guides, state.vector);
                        let svg = render_crispr_locus_svg(
                            state.target_len_bp,
                            &state.guides,
                            "OSABT CRISPR Locus Guide Map",
                        );
                        if fs::write("crispr_guides.tsv", tsv).is_ok()
                            && fs::write("crispr_locus.svg", svg).is_ok()
                        {
                            state.status_msg =
                                String::from("Exported crispr_guides.tsv and crispr_locus.svg!");
                        } else {
                            state.status_msg =
                                String::from("Export failed: check write permissions.");
                        }
                    }
                    _ => {}
                }
            }
        }
    }
}

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &CrisprTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Ranked Guides (Rule Set 2 + CFD)",
        "[2] Off-Target CFD Mismatch Inspector",
        "[3] Base-Editor (CBE/ABE) Bystanders",
        "[4] Golden Gate Oligo Order Sheet",
    ]
    .into_iter()
    .map(Line::from)
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        CrisprTab::RankedGuides => 0,
        CrisprTab::OffTargetCfdInspector => 1,
        CrisprTab::BaseEditorWindows => 2,
        CrisprTab::GoldenGateOligos => 3,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT CRISPR Studio (v1.7.0) — {} | {} Guides on {} bp Target ",
            state.nuclease.label(),
            state.guides.len(),
            state.target_len_bp
        )))
        .select(idx)
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        CrisprTab::RankedGuides => draw_guides_tab(f, chunks[1], state),
        CrisprTab::OffTargetCfdInspector => draw_offtarget_tab(f, chunks[1], state),
        CrisprTab::BaseEditorWindows => draw_base_edit_tab(f, chunks[1], state),
        CrisprTab::GoldenGateOligos => draw_oligos_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls & Status "))
        .style(Style::default().fg(Color::Gray));
    f.render_widget(footer, chunks[2]);
}

fn draw_guides_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &CrisprTuiState) {
    let hdr = Row::new([
        "Sel",
        "ID",
        "Str",
        "Cut",
        "Spacer (5'->3')",
        "PAM",
        "GC%",
        "R-Loop Tm",
        "On-Target",
        "CFD Spec",
        "Offs",
        "Composite",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );

    let start = state.selected_idx.saturating_sub(8);
    let rows = state
        .guides
        .iter()
        .enumerate()
        .skip(start)
        .map(|(i, g)| {
            let is_sel = i == state.selected_idx;
            let style = if is_sel {
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD)
            } else {
                Style::default()
            };
            Row::new([
                Cell::from(if is_sel { ">" } else { " " }),
                Cell::from(g.guide_id.clone()),
                Cell::from(format!("{}", g.strand)),
                Cell::from(format!("{}", g.cut_site_1idx)),
                Cell::from(g.spacer_seq.clone()),
                Cell::from(g.pam_seq.clone()),
                Cell::from(format!("{:.0}%", g.gc_percent)),
                Cell::from(format!("{:.1}C", g.rna_dna_tm_c)),
                Cell::from(format!("{:.1}", g.on_target_score)),
                Cell::from(format!("{:.1}", g.cfd_specificity_score)),
                Cell::from(format!("{}", g.off_target_count)),
                Cell::from(format!("{:.1}", g.composite_score)),
            ])
            .style(style)
        });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(4),
            Constraint::Length(9),
            Constraint::Length(4),
            Constraint::Length(6),
            Constraint::Length(23),
            Constraint::Length(8),
            Constraint::Length(6),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Length(6),
            Constraint::Min(10),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Doench et al. (2016) On-Target Efficiency & CFD Specificity Ranking "),
    );
    f.render_widget(tbl, area);
}

fn draw_offtarget_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &CrisprTuiState) {
    let Some(g) = state.guides.get(state.selected_idx) else {
        return;
    };
    let hits = &state.offtargets[state.selected_idx];

    let mut lines = vec![
        Line::from(format!(
            "Selected Guide: {} ({} strand, cut @ {} bp) | Spacer = 5'-{}-3' + PAM {} | Aggregate Specificity = {:.1}/100",
            g.guide_id,
            g.strand,
            g.cut_site_1idx,
            g.spacer_seq,
            g.pam_seq,
            g.cfd_specificity_score
        )),
        Line::from(""),
    ];

    if hits.is_empty() {
        lines.push(Line::from(
            "  Zero off-target hits (<= 4 mismatches) detected across the background FASTA!",
        ));
    } else {
        lines.push(Line::from(
            "  Contig                         Pos(Str)   MM   Alignment vs Spacer    PAM    Doench CFD",
        ));
        for h in hits.iter().take(14) {
            lines.push(Line::from(format!(
                "  {:28} {:6}({})   {:2}   {:20}   {:4}   {:.4}",
                h.contig, h.pos_1idx, h.strand, h.mismatches, h.mismatch_mask, h.pam_seq, h.cfd_score
            )));
        }
    }

    let card = Paragraph::new(lines).wrap(Wrap { trim: false }).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Doench et al. (2016) Cutting Frequency Determination (CFD) Off-Target Alignments "),
    );
    f.render_widget(card, area);
}

fn draw_base_edit_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &CrisprTuiState) {
    let hdr = Row::new([
        "Guide",
        "Spacer (5'->3')",
        "Window [4..8]",
        "CBE (C->T, Komor 2016)",
        "ABE (A->G, Gaudelli 2017)",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state.guides.iter().skip(state.selected_idx.saturating_sub(6)).map(|g| {
        let cbe = if g.base_edit.cbe_clean_single_edit {
            format!("PRECISION 1x C @ pos {}", g.base_edit.cbe_c_positions[0])
        } else if g.base_edit.cbe_c_positions.is_empty() {
            String::from("No C in window")
        } else {
            format!("Bystander C @ {:?}", g.base_edit.cbe_c_positions)
        };
        let abe = if g.base_edit.abe_clean_single_edit {
            format!("PRECISION 1x A @ pos {}", g.base_edit.abe_a_positions[0])
        } else if g.base_edit.abe_a_positions.is_empty() {
            String::from("No A in window")
        } else {
            format!("Bystander A @ {:?}", g.base_edit.abe_a_positions)
        };
        Row::new([
            g.guide_id.clone(),
            g.spacer_seq.clone(),
            g.base_edit.window_seq_4_to_8.clone(),
            cbe,
            abe,
        ])
    });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(9),
            Constraint::Length(23),
            Constraint::Length(15),
            Constraint::Length(28),
            Constraint::Min(28),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Cytosine (CBE) & Adenine (ABE) Deamination Window & Bystander Audit "),
    );
    f.render_widget(tbl, area);
}

fn draw_oligos_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &CrisprTuiState) {
    let hdr = Row::new([
        "Guide",
        "+G?",
        "Top Annealed Oligo (5'->3')",
        "Bottom Annealed Oligo (5'->3')",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );
    let rows = state.guides.iter().take(14).map(|g| {
        let ol = design_golden_gate_crispr_oligos(&g.spacer_seq, state.vector);
        Row::new([
            g.guide_id.clone(),
            String::from(if ol.prepended_u6_g { "Yes" } else { "No" }),
            ol.top_oligo_5to3,
            ol.bottom_oligo_5to3,
        ])
    });
    let tbl = Table::new(
        rows,
        [
            Constraint::Length(9),
            Constraint::Length(6),
            Constraint::Length(32),
            Constraint::Min(32),
        ],
    )
    .header(hdr)
    .block(Block::default().borders(Borders::ALL).title(format!(
        " Golden Gate Cloning Oligo Order Sheet — Vector: {} ([v] to switch) ",
        state.vector.label()
    )));
    f.render_widget(tbl, area);
}
