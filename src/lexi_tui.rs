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

//! Interactive 4-tab `Ratatui` terminal biological & anatomical etymology studio for Tool #6: `osabt lexi`.

use std::fs;
use std::io;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::Terminal;

use crate::tui_util::run_tui_session;

use crate::lexi_engine::{
    build_cognate_tree, decode_biological_phrase, demo_etymology_terms, export_anki_flashcards_tsv,
    export_flashcards_markdown, search_morpheme_lexicon, DecodedPhrase,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum LexiTab {
    Decoder,
    CognateTree,
    RootBrowser,
    FlashcardQuiz,
}

struct LexiTuiState {
    active_tab: LexiTab,
    input_phrase: String,
    editing_input: bool,
    lexicon_filter: String,
    saved_deck: Vec<DecodedPhrase>,
    quiz_idx: usize,
    quiz_revealed: bool,
    preset_idx: usize,
    status_msg: String,
}

impl LexiTuiState {
    fn new(initial_phrase: Option<String>) -> Self {
        let presets = demo_etymology_terms();
        let saved_deck: Vec<DecodedPhrase> = presets
            .iter()
            .map(|t| decode_biological_phrase(t))
            .collect();
        let start_phrase = initial_phrase.unwrap_or_else(|| String::from(presets[0]));

        Self {
            active_tab: LexiTab::Decoder,
            input_phrase: start_phrase,
            editing_input: false,
            lexicon_filter: String::new(),
            saved_deck,
            quiz_idx: 0,
            quiz_revealed: false,
            preset_idx: 0,
            status_msg: String::from(
                "Tab/1-4: switch view | [i]: type term | [n]: next demo taxon | [a]: add to deck | [e]: export Anki TSV + Markdown | [q]: quit",
            ),
        }
    }

    fn export_deck(&mut self) {
        let anki = export_anki_flashcards_tsv(&self.saved_deck, "OSABT_Biology_Etymology");
        let md = export_flashcards_markdown(&self.saved_deck);
        let r1 = fs::write("osabt_etymology_anki.tsv", anki);
        let r2 = fs::write("osabt_etymology_study.md", md);
        if r1.is_ok() && r2.is_ok() {
            self.status_msg = format!(
                "Exported {} cards -> ./osabt_etymology_anki.tsv & ./osabt_etymology_study.md",
                self.saved_deck.len()
            );
        } else {
            self.status_msg = String::from("Failed writing flashcard output files.");
        }
    }
}

/// Launches the interactive 4-tab `Ratatui` biological etymology & study flashcard studio.
pub fn run_lexi_tui(initial_phrase: Option<String>) -> Result<(), String> {
    let mut state = LexiTuiState::new(initial_phrase);
    run_tui_session(|terminal| lexi_loop(terminal, &mut state))
}

fn lexi_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut LexiTuiState,
) -> Result<(), String> {
    loop {
        terminal
            .draw(|f| draw_lexi_ui(f, state))
            .map_err(|e| format!("Draw error: {e}"))?;

        if event::poll(Duration::from_millis(120)).map_err(|e| e.to_string())? {
            if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }

                if state.editing_input {
                    match key.code {
                        KeyCode::Enter | KeyCode::Esc => {
                            state.editing_input = false;
                            state.status_msg = String::from(
                                "Decoded input. Press [a] to add to study deck or [e] to export Anki TSV.",
                            );
                        }
                        KeyCode::Backspace => {
                            if state.active_tab == LexiTab::RootBrowser {
                                state.lexicon_filter.pop();
                            } else {
                                state.input_phrase.pop();
                            }
                        }
                        KeyCode::Char(c) => {
                            if state.active_tab == LexiTab::RootBrowser {
                                state.lexicon_filter.push(c);
                            } else {
                                state.input_phrase.push(c);
                            }
                        }
                        _ => {}
                    }
                    continue;
                }

                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                    KeyCode::Char('1') => state.active_tab = LexiTab::Decoder,
                    KeyCode::Char('2') => state.active_tab = LexiTab::CognateTree,
                    KeyCode::Char('3') => state.active_tab = LexiTab::RootBrowser,
                    KeyCode::Char('4') => state.active_tab = LexiTab::FlashcardQuiz,
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            LexiTab::Decoder => LexiTab::CognateTree,
                            LexiTab::CognateTree => LexiTab::RootBrowser,
                            LexiTab::RootBrowser => LexiTab::FlashcardQuiz,
                            LexiTab::FlashcardQuiz => LexiTab::Decoder,
                        };
                    }
                    KeyCode::Char('i') | KeyCode::Char('/') => {
                        state.editing_input = true;
                        state.status_msg =
                            String::from("Type scientific/anatomical term (Enter/Esc when done):");
                    }
                    KeyCode::Char('n') => {
                        let presets = demo_etymology_terms();
                        state.preset_idx = (state.preset_idx + 1) % presets.len();
                        state.input_phrase = presets[state.preset_idx].to_string();
                        if state.active_tab == LexiTab::FlashcardQuiz && !state.saved_deck.is_empty()
                        {
                            state.quiz_idx = (state.quiz_idx + 1) % state.saved_deck.len();
                            state.quiz_revealed = false;
                        }
                    }
                    KeyCode::Char(' ') => {
                        if state.active_tab == LexiTab::FlashcardQuiz {
                            state.quiz_revealed = !state.quiz_revealed;
                        }
                    }
                    KeyCode::Char('a') => {
                        let dec = decode_biological_phrase(&state.input_phrase);
                        if !dec.raw_input.is_empty() {
                            state.saved_deck.push(dec);
                            state.status_msg = format!(
                                "Added '{}' to flashcard study deck ({} total cards).",
                                state.input_phrase,
                                state.saved_deck.len()
                            );
                        }
                    }
                    KeyCode::Char('e') => state.export_deck(),
                    _ => {}
                }
            }
        }
    }
}

fn draw_lexi_ui(f: &mut ratatui::Frame<'_>, state: &LexiTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Morpheme Decoder",
        "[2] Cognate Tree",
        "[3] Greek/Latin Lexicon",
        "[4] Study Quiz & Anki Deck",
    ]
    .iter()
    .map(|t| Line::from(*t))
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        LexiTab::Decoder => 0,
        LexiTab::CognateTree => 1,
        LexiTab::RootBrowser => 2,
        LexiTab::FlashcardQuiz => 3,
    };

    let tabs = Tabs::new(titles)
        .select(idx)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT Taxonomic & Anatomical Etymology Studio • Study Deck: {} cards ",
            state.saved_deck.len()
        )))
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        LexiTab::Decoder => draw_decoder_tab(f, chunks[1], state),
        LexiTab::CognateTree => draw_cognate_tab(f, chunks[1], state),
        LexiTab::RootBrowser => draw_lexicon_tab(f, chunks[1], state),
        LexiTab::FlashcardQuiz => draw_quiz_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls "));
    f.render_widget(footer, chunks[2]);
}

fn draw_decoder_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &LexiTuiState) {
    let decoded = decode_biological_phrase(&state.input_phrase);
    let split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(5), Constraint::Min(8)])
        .split(area);

    let prompt_style = if state.editing_input {
        Style::default()
            .fg(Color::Green)
            .add_modifier(Modifier::BOLD)
    } else {
        Style::default()
            .fg(Color::White)
            .add_modifier(Modifier::BOLD)
    };
    let top_lines = vec![
        Line::from(vec![
            Span::styled("Query Term ([i] edit, [n] next preset): ", Style::default().fg(Color::Yellow)),
            Span::styled(&state.input_phrase, prompt_style),
        ]),
        Line::from(vec![
            Span::styled("Literal Synthesis: ", Style::default().fg(Color::Cyan)),
            Span::raw(&decoded.combined_gloss),
        ]),
    ];
    let top_box = Paragraph::new(top_lines)
        .wrap(Wrap { trim: false })
        .block(Block::default().borders(Borders::ALL).title(" Input & Literal Meaning "));
    f.render_widget(top_box, split[0]);

    let header = Row::new(["Word", "Slice", "Role", "Origin", "Classical Lemma", "Meaning"])
        .style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );

    let mut rows = Vec::new();
    for w in &decoded.words {
        for seg in &w.segments {
            rows.push(Row::new(vec![
                Cell::from(w.original_word.clone()),
                Cell::from(seg.surface_slice.clone()),
                Cell::from(seg.role_label.clone()),
                Cell::from(seg.origin_label.clone()),
                Cell::from(seg.classical_lemma.clone()),
                Cell::from(seg.gloss.clone()),
            ]));
        }
    }

    let widths = [
        Constraint::Length(18),
        Constraint::Length(12),
        Constraint::Length(10),
        Constraint::Length(18),
        Constraint::Length(26),
        Constraint::Min(24),
    ];

    let tbl = Table::new(rows, widths).header(header).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Dynamic-Programming Morpheme Segmentation "),
    );
    f.render_widget(tbl, split[1]);
}

fn draw_cognate_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &LexiTuiState) {
    let decoded = decode_biological_phrase(&state.input_phrase);
    let tree = build_cognate_tree(&decoded);
    let mut lines = vec![
        Line::from(Span::styled(
            format!("Cognate Network for: {}", decoded.raw_input),
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        )),
        Line::from(""),
    ];

    for node in &tree {
        lines.push(Line::from(Span::styled(
            format!(
                "┌─ Root `{}` [{} — {}] = \"{}\"",
                node.morpheme_id, node.origin, node.classical_lemma, node.gloss
            ),
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        )));
        for (idx, cog) in node.cognates.iter().enumerate() {
            let branch = if idx + 1 == node.cognates.len() {
                "└──"
            } else {
                "├──"
            };
            lines.push(Line::from(format!("   {branch} {cog}")));
        }
        lines.push(Line::from(""));
    }

    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Shared Morpheme Cognate Tree "),
        );
    f.render_widget(p, area);
}

fn draw_lexicon_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &LexiTuiState) {
    let hits = search_morpheme_lexicon(&state.lexicon_filter);
    let header = Row::new(["Morpheme", "Role", "Origin", "Classical Lemma", "Gloss", "Cognates"])
        .style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );

    let rows = hits.iter().take(22).map(|e| {
        Row::new(vec![
            Cell::from(e.canonical_id),
            Cell::from(e.role.label()),
            Cell::from(e.origin.label()),
            Cell::from(e.classical_lemma),
            Cell::from(e.gloss),
            Cell::from(e.cognates.join(", ")),
        ])
    });

    let widths = [
        Constraint::Length(11),
        Constraint::Length(10),
        Constraint::Length(17),
        Constraint::Length(24),
        Constraint::Length(28),
        Constraint::Min(24),
    ];

    let title = format!(
        " Offline Greek/Latin Lexicon ({} matches | press [/] to filter: '{}') ",
        hits.len(),
        state.lexicon_filter
    );
    let tbl = Table::new(rows, widths)
        .header(header)
        .block(Block::default().borders(Borders::ALL).title(title));
    f.render_widget(tbl, area);
}

fn draw_quiz_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &LexiTuiState) {
    let mut lines = Vec::new();
    if let Some(card) = state.saved_deck.get(state.quiz_idx) {
        lines.push(Line::from(Span::styled(
            format!(
                "Flashcard [{}/{}] — Press [Space] to flip, [n] for next card, [e] to export Anki TSV",
                state.quiz_idx + 1,
                state.saved_deck.len()
            ),
            Style::default().fg(Color::Yellow),
        )));
        lines.push(Line::from(""));
        lines.push(Line::from(Span::styled(
            format!("FRONT:  {}", card.raw_input),
            Style::default()
                .fg(Color::White)
                .add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(""));
        if state.quiz_revealed {
            lines.push(Line::from(Span::styled(
                format!("BACK :  {}", card.combined_gloss),
                Style::default()
                    .fg(Color::Green)
                    .add_modifier(Modifier::BOLD),
            )));
            for w in &card.words {
                lines.push(Line::from(format!("  • {} ({})", w.original_word, w.hyphenated)));
                for seg in &w.segments {
                    if !seg.is_connecting_vowel {
                        lines.push(Line::from(format!(
                            "      - `{}` [{}]: {}",
                            seg.surface_slice, seg.classical_lemma, seg.gloss
                        )));
                    }
                }
            }
        } else {
            lines.push(Line::from(Span::styled(
                "[Hidden — press Spacebar to reveal morpheme breakdown]",
                Style::default().fg(Color::DarkGray),
            )));
        }
    }

    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Active-Recall Etymology Flashcard Deck "),
        );
    f.render_widget(p, area);
}
