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

//! Interactive `less`-style SSH/terminal phylogenetic explorer (`Ratatui`) for Tool #4: `osabt phylo`.
//! Supports scrolling large Newick/NEXUS trees, `/`-searching taxa, collapsing/expanding clades,
//! outgroup & midpoint rerooting, ladderizing, cladogram/phylogram toggling, and `.svg`/`.nwk` export.

use std::fs;
use std::io;
use std::path::PathBuf;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::{Line, Span};
use ratatui::widgets::{Block, Borders, Paragraph, Wrap};
use ratatui::Terminal;

use crate::tui_util::run_tui_session;

use crate::phylo_core::{parse_phylo_trees, PhyloTree};
use crate::phylo_render::{render_phylo_svg, render_unicode_tree, RenderedTreeLine};

struct PhyloTuiState {
    trees: Vec<PhyloTree>,
    tree_idx: usize,
    cursor_line: usize,
    scroll_offset: usize,
    proportional_lengths: bool,
    ladder_asc: bool,
    search_mode: bool,
    search_query: String,
    rendered: Vec<RenderedTreeLine>,
    status_msg: String,
}

impl PhyloTuiState {
    fn new(input: Option<PathBuf>) -> Result<Self, String> {
        let trees = match input {
            Some(p) => {
                let raw = fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                let stem = p
                    .file_stem()
                    .map(|s| s.to_string_lossy().to_string())
                    .unwrap_or_else(|| String::from("Tree"));
                parse_phylo_trees(&stem, &raw)?
            }
            None => vec![PhyloTree::demo_vertebrate_amniotes()],
        };

        let mut st = Self {
            trees,
            tree_idx: 0,
            cursor_line: 0,
            scroll_offset: 0,
            proportional_lengths: true,
            ladder_asc: true,
            search_mode: false,
            search_query: String::new(),
            rendered: Vec::new(),
            status_msg: String::from(
                "j/k: move | Space/c: collapse clade | r: reroot | m: midpoint | l: ladderize | p: phylogram/cladogram | /: search | e: export SVG+NWK | q: quit",
            ),
        };
        st.rebuild_render(50);
        Ok(st)
    }

    fn current_tree(&self) -> &PhyloTree {
        &self.trees[self.tree_idx]
    }

    fn current_tree_mut(&mut self) -> &mut PhyloTree {
        &mut self.trees[self.tree_idx]
    }

    fn rebuild_render(&mut self, branch_cols: usize) {
        let q = self.search_query.clone();
        let prop = self.proportional_lengths;
        self.rendered = render_unicode_tree(self.current_tree(), prop, branch_cols, &q);
        if !self.rendered.is_empty() {
            self.cursor_line = self.cursor_line.min(self.rendered.len() - 1);
        } else {
            self.cursor_line = 0;
        }
    }

    fn selected_node_id(&self) -> Option<usize> {
        self.rendered.get(self.cursor_line).and_then(|r| r.node_id)
    }

    fn export_current(&mut self) {
        let name = self.current_tree().tree_name.clone();
        let rect_svg = render_phylo_svg(self.current_tree(), self.proportional_lengths, false);
        let circ_svg = render_phylo_svg(self.current_tree(), self.proportional_lengths, true);
        let nwk = self.current_tree().to_newick();

        let r1 = fs::write(format!("{name}_rect.svg"), rect_svg);
        let r2 = fs::write(format!("{name}_circular.svg"), circ_svg);
        let r3 = fs::write(format!("{name}_edited.nwk"), nwk);

        if r1.is_ok() && r2.is_ok() && r3.is_ok() {
            self.status_msg = format!(
                "Exported ./{name}_rect.svg, ./{name}_circular.svg, and ./{name}_edited.nwk"
            );
        } else {
            self.status_msg = String::from("Export error writing .svg/.nwk files.");
        }
    }
}

/// Launches the interactive SSH/terminal `Ratatui` tree explorer.
pub fn run_phylo_tui(input: Option<PathBuf>) -> Result<(), String> {
    let mut state = PhyloTuiState::new(input)?;
    run_tui_session(|terminal| event_loop(terminal, &mut state))
}

fn event_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut PhyloTuiState,
) -> Result<(), String> {
    loop {
        let term_w = terminal.size().map(|r| r.width as usize).unwrap_or(110);
        let branch_cols = term_w.saturating_sub(68).clamp(24, 84);
        state.rebuild_render(branch_cols);

        terminal
            .draw(|f| draw_phylo_screen(f, state))
            .map_err(|e| format!("Draw error: {e}"))?;

        if event::poll(Duration::from_millis(120)).map_err(|e| e.to_string())? {
            if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }

                if state.search_mode {
                    match key.code {
                        KeyCode::Enter | KeyCode::Esc => {
                            state.search_mode = false;
                            // Jump to first search hit if available
                            if let Some(hit_idx) =
                                state.rendered.iter().position(|l| l.is_search_hit)
                            {
                                state.cursor_line = hit_idx;
                            }
                        }
                        KeyCode::Backspace => {
                            state.search_query.pop();
                        }
                        KeyCode::Char(c) => {
                            state.search_query.push(c);
                        }
                        _ => {}
                    }
                    continue;
                }

                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                    KeyCode::Char('/') => {
                        state.search_mode = true;
                        state.status_msg =
                            String::from("Search taxon substring (Enter/Esc to finish):");
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        if state.cursor_line + 1 < state.rendered.len() {
                            state.cursor_line += 1;
                        }
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.cursor_line = state.cursor_line.saturating_sub(1);
                    }
                    KeyCode::Char(' ') | KeyCode::Char('c') => {
                        if let Some(nid) = state.selected_node_id() {
                            let target = if state.current_tree().nodes[nid].is_leaf() {
                                state.current_tree().nodes[nid].parent.unwrap_or(nid)
                            } else {
                                nid
                            };
                            let collapsed = state.current_tree_mut().toggle_collapse(target);
                            state.status_msg = if collapsed {
                                format!("Collapsed clade at node #{target}.")
                            } else {
                                format!("Expanded clade at node #{target}.")
                            };
                        }
                    }
                    KeyCode::Char('r') => {
                        if let Some(nid) = state.selected_node_id() {
                            match state.current_tree_mut().reroot_at_node(nid, 0.5) {
                                Ok(()) => {
                                    state.status_msg = format!("Rerooted tree on branch #{nid}.")
                                }
                                Err(e) => state.status_msg = format!("Reroot notice: {e}"),
                            }
                        }
                    }
                    KeyCode::Char('m') => match state.current_tree_mut().midpoint_reroot() {
                        Ok(d) => {
                            state.status_msg = format!(
                                "Midpoint rerooted at half max patristic diameter ({:.4} subst/site).",
                                d * 0.5
                            );
                        }
                        Err(e) => state.status_msg = format!("Midpoint error: {e}"),
                    },
                    KeyCode::Char('l') => {
                        state.ladder_asc = !state.ladder_asc;
                        let asc = state.ladder_asc;
                        state.current_tree_mut().ladderize(asc);
                        state.status_msg = format!(
                            "Ladderized tree ({}).",
                            if asc { "ascending" } else { "descending" }
                        );
                    }
                    KeyCode::Char('p') => {
                        state.proportional_lengths = !state.proportional_lengths;
                    }
                    KeyCode::Char('[') => {
                        if state.trees.len() > 1 {
                            state.tree_idx = state
                                .tree_idx
                                .checked_sub(1)
                                .unwrap_or(state.trees.len() - 1);
                            state.cursor_line = 0;
                        }
                    }
                    KeyCode::Char(']') => {
                        if state.trees.len() > 1 {
                            state.tree_idx = (state.tree_idx + 1) % state.trees.len();
                            state.cursor_line = 0;
                        }
                    }
                    KeyCode::Char('e') => state.export_current(),
                    _ => {}
                }
            }
        }
    }
}

fn draw_phylo_screen(f: &mut ratatui::Frame<'_>, state: &mut PhyloTuiState) {
    let rows = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(10),
            Constraint::Length(3),
        ])
        .split(f.area());

    let st = state.current_tree().stats();
    let header_txt = format!(
        " Tree [{}/{}]: {} | Tips: {} | Internal: {} | Total BL: {:.4} | Max Height: {:.4} | Mode: {} ",
        state.tree_idx + 1,
        state.trees.len(),
        state.current_tree().tree_name,
        st.leaf_count,
        st.internal_count,
        st.total_branch_length,
        st.max_root_to_tip,
        if state.proportional_lengths {
            "Phylogram"
        } else {
            "Cladogram"
        }
    );
    let header = Paragraph::new(header_txt)
        .style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        )
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" OSABT Terminal Phylogenetic Explorer "),
        );
    f.render_widget(header, rows[0]);

    let cols = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(68), Constraint::Percentage(32)])
        .split(rows[1]);

    draw_tree_pane(f, cols[0], state);
    draw_inspector_pane(f, cols[1], state);

    let footer_text = if state.search_mode {
        format!("/{}", state.search_query)
    } else {
        state.status_msg.clone()
    };
    let footer = Paragraph::new(footer_text)
        .block(Block::default().borders(Borders::ALL).title(" Controls "));
    f.render_widget(footer, rows[2]);
}

fn draw_tree_pane(f: &mut ratatui::Frame<'_>, area: Rect, state: &mut PhyloTuiState) {
    let inner_h = area.height.saturating_sub(2) as usize;
    if inner_h > 0 {
        if state.cursor_line < state.scroll_offset {
            state.scroll_offset = state.cursor_line;
        } else if state.cursor_line >= state.scroll_offset + inner_h {
            state.scroll_offset = state.cursor_line + 1 - inner_h;
        }
    }

    let visible_lines: Vec<Line<'_>> = state
        .rendered
        .iter()
        .enumerate()
        .skip(state.scroll_offset)
        .take(inner_h.max(1))
        .map(|(idx, rline)| {
            let is_cur = idx == state.cursor_line;
            let mut style = Style::default().fg(Color::White);
            if rline.is_search_hit {
                style = style.fg(Color::Yellow).add_modifier(Modifier::BOLD);
            }
            if is_cur {
                style = style.bg(Color::DarkGray).add_modifier(Modifier::BOLD);
            }
            let prefix = if is_cur { "▶ " } else { "  " };
            Line::from(Span::styled(format!("{prefix}{}", rline.text), style))
        })
        .collect();

    let p = Paragraph::new(visible_lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Tree Canvas "),
    );
    f.render_widget(p, area);
}

fn draw_inspector_pane(f: &mut ratatui::Frame<'_>, area: Rect, state: &PhyloTuiState) {
    let mut lines = Vec::new();
    if let Some(nid) = state.selected_node_id() {
        let tree = state.current_tree();
        let node = &tree.nodes[nid];
        let tips = tree.leaf_names_under(nid);

        lines.push(Line::from(Span::styled(
            format!("Selected Node #{nid}"),
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(format!(
            "Type       : {}",
            if node.is_leaf() {
                "Terminal Taxon (Tip)"
            } else if nid == tree.root {
                "Root Clade"
            } else {
                "Internal Clade"
            }
        )));
        if !node.name.is_empty() {
            lines.push(Line::from(format!("Label      : {}", node.name)));
        }
        lines.push(Line::from(format!(
            "Branch Len : {}",
            node.branch_length
                .map(|b| format!("{b:.6}"))
                .unwrap_or_else(|| String::from("(none)"))
        )));
        lines.push(Line::from(format!(
            "Support    : {}",
            node.support_pct
                .map(|s| format!("{s:.1}%"))
                .unwrap_or_else(|| String::from("(none)"))
        )));
        lines.push(Line::from(format!("Descendant Tips ({}):", tips.len())));
        for t in tips.iter().take(12) {
            lines.push(Line::from(format!("  • {}", t.replace('_', " "))));
        }
        if tips.len() > 12 {
            lines.push(Line::from(format!("  ... +{} more", tips.len() - 12)));
        }
    } else {
        lines.push(Line::from("Move j/k onto a branch or tip row to inspect clade details."));
    }

    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Clade Inspector "),
        );
    f.render_widget(p, area);
}
