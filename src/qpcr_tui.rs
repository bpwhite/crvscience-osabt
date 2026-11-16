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

//! Interactive 4-tab `Ratatui` terminal MIQE RT-qPCR Quantification, `geNorm`/`qBase`,
//! & SYBR `-dF/dT` Melt-Curve Studio (`osabt qpcr`).

use std::fs;
use std::io;
use std::path::PathBuf;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::Line;
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs};
use ratatui::Terminal;

use crate::tui_util::run_tui_session;

use crate::qpcr_core::{
    analyze_qpcr_experiment, demo_qpcr_heatshock_plate_csv, parse_qpcr_table,
    render_qpcr_expression_svg, MiqeFlagSeverity, QpcrDatasetAnalysis, QpcrReactionWell,
};
use crate::qpcr_melt::{
    analyze_melt_curve, demo_sybr_melt_traces, parse_melt_csv, render_melt_curves_svg,
    MeltCurveResult,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum QpcrTab {
    ExpressionFoldChange,
    EfficiencyAndGeNorm,
    MeltCurves,
    MiqeAudit,
}

struct QpcrTuiState {
    active_tab: QpcrTab,
    wells: Vec<QpcrReactionWell>,
    calibrator_group: String,
    include_unstable_b2m: bool,
    ref_genes: Vec<String>,
    analysis: QpcrDatasetAnalysis,
    melt_results: Vec<MeltCurveResult>,
    scroll_offset: usize,
    status_msg: String,
}

impl QpcrTuiState {
    fn new(
        cq_path: Option<PathBuf>,
        melt_path: Option<PathBuf>,
        refs: Vec<String>,
        calibrator: String,
    ) -> Result<Self, String> {
        let wells = if let Some(p) = cq_path {
            let raw = fs::read_to_string(&p)
                .map_err(|e| format!("Failed reading qPCR table '{}': {e}", p.display()))?;
            parse_qpcr_table(&raw)?
        } else {
            parse_qpcr_table(demo_qpcr_heatshock_plate_csv())?
        };

        let traces = if let Some(mp) = melt_path {
            let raw = fs::read_to_string(&mp)
                .map_err(|e| format!("Failed reading melt CSV '{}': {e}", mp.display()))?;
            parse_melt_csv(&raw)?
        } else {
            demo_sybr_melt_traces()
        };
        let melt_results: Vec<MeltCurveResult> = traces
            .iter()
            .map(|t| analyze_melt_curve(t, 77.5))
            .collect();

        let ref_genes = if refs.is_empty() {
            vec![String::from("ACTB"), String::from("GAPDH")]
        } else {
            refs
        };
        let analysis = analyze_qpcr_experiment(&wells, &ref_genes, &calibrator)?;

        Ok(Self {
            active_tab: QpcrTab::ExpressionFoldChange,
            wells,
            calibrator_group: calibrator,
            include_unstable_b2m: false,
            ref_genes,
            analysis,
            melt_results,
            scroll_offset: 0,
            status_msg: String::from(
                "Keys: [1-4/Tab] Views | [Up/Down] Scroll | [b] Toggle Unstable B2M Ref | [e] Export Expression & Melt SVGs | [q] Quit",
            ),
        })
    }

    fn toggle_b2m_ref(&mut self) {
        self.include_unstable_b2m = !self.include_unstable_b2m;
        let mut next_refs = vec![String::from("ACTB"), String::from("GAPDH")];
        if self.include_unstable_b2m {
            next_refs.push(String::from("B2M"));
        }
        if let Ok(rep) = analyze_qpcr_experiment(&self.wells, &next_refs, &self.calibrator_group) {
            self.ref_genes = next_refs;
            self.analysis = rep;
            self.status_msg = format!(
                "Active Reference Panel: {} (watch geNorm & MIQE tab!)",
                self.ref_genes.join(" + ")
            );
        }
    }

    fn export_svgs(&mut self) {
        let expr_svg = render_qpcr_expression_svg(&self.analysis, "OSABT RT-qPCR Analysis");
        let melt_svg = render_melt_curves_svg(&self.melt_results, "OSABT SYBR Dissociation");
        if fs::write("qpcr_expression.svg", expr_svg).is_ok()
            && fs::write("qpcr_melt_curves.svg", melt_svg).is_ok()
        {
            self.status_msg =
                String::from("Exported qpcr_expression.svg and qpcr_melt_curves.svg!");
        } else {
            self.status_msg = String::from("SVG export failed: check directory permissions.");
        }
    }
}

/// Launches the interactive 4-tab `Ratatui` terminal studio for `osabt qpcr`.
pub fn run_qpcr_tui(
    cq_path: Option<PathBuf>,
    melt_path: Option<PathBuf>,
    refs: Vec<String>,
    calibrator: String,
) -> Result<(), String> {
    let mut state = QpcrTuiState::new(cq_path, melt_path, refs, calibrator)?;
    run_tui_session(|terminal| run_loop(terminal, &mut state))
}

fn run_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut QpcrTuiState,
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
                        state.active_tab = QpcrTab::ExpressionFoldChange;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('2') => {
                        state.active_tab = QpcrTab::EfficiencyAndGeNorm;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('3') => {
                        state.active_tab = QpcrTab::MeltCurves;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('4') => {
                        state.active_tab = QpcrTab::MiqeAudit;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            QpcrTab::ExpressionFoldChange => QpcrTab::EfficiencyAndGeNorm,
                            QpcrTab::EfficiencyAndGeNorm => QpcrTab::MeltCurves,
                            QpcrTab::MeltCurves => QpcrTab::MiqeAudit,
                            QpcrTab::MiqeAudit => QpcrTab::ExpressionFoldChange,
                        };
                        state.scroll_offset = 0;
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        state.scroll_offset = state.scroll_offset.saturating_add(1);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.scroll_offset = state.scroll_offset.saturating_sub(1);
                    }
                    KeyCode::Char('b') => state.toggle_b2m_ref(),
                    KeyCode::Char('e') => state.export_svgs(),
                    _ => {}
                }
            }
        }
    }
}

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &QpcrTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] qBase / Pfaffl / 2^-ddCq",
        "[2] Std Curves & geNorm M",
        "[3] SYBR -dF/dT Melt Peaks",
        "[4] Bustin (2009) MIQE Audit",
    ]
    .into_iter()
    .map(Line::from)
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        QpcrTab::ExpressionFoldChange => 0,
        QpcrTab::EfficiencyAndGeNorm => 1,
        QpcrTab::MeltCurves => 2,
        QpcrTab::MiqeAudit => 3,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT qPCR Studio (v1.2.0) — Refs: [{}] | Calibrator: {} ",
            state.ref_genes.join(", "),
            state.calibrator_group
        )))
        .select(idx)
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        QpcrTab::ExpressionFoldChange => draw_expression_tab(f, chunks[1], state),
        QpcrTab::EfficiencyAndGeNorm => draw_efficiency_genorm_tab(f, chunks[1], state),
        QpcrTab::MeltCurves => draw_melt_tab(f, chunks[1], state),
        QpcrTab::MiqeAudit => draw_miqe_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls & Status "))
        .style(Style::default().fg(Color::Gray));
    f.render_widget(footer, chunks[2]);
}

fn draw_expression_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &QpcrTuiState) {
    let split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(8), Constraint::Min(6)])
        .split(area);

    let grp_hdr = Row::new([
        "Target",
        "Comparison",
        "qBase FC (95% CI)",
        "log2(FC)",
        "Pfaffl FC",
        "Livak 2^-ddCq",
        "Welch t",
        "p-value",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );

    let grp_rows = state.analysis.group_comparisons.iter().map(|c| {
        Row::new([
            c.target_gene.clone(),
            format!("{} vs {}", c.treatment_group, c.calibrator_group),
            format!(
                "{:.2}x ({:.2}–{:.2})",
                c.fold_change_qbase, c.fold_change_ci95.0, c.fold_change_ci95.1
            ),
            format!("{:+.3}", c.log2_fold_change),
            format!("{:.2}x", c.pfaffl_ratio),
            format!("{:.2}x", c.livak_ratio),
            format!("{:+.2}", c.welch_t_stat),
            format!("{:.2e}", c.p_value_two_sided),
        ])
    });

    let grp_tbl = Table::new(
        grp_rows,
        [
            Constraint::Length(9),
            Constraint::Length(22),
            Constraint::Length(21),
            Constraint::Length(10),
            Constraint::Length(11),
            Constraint::Length(14),
            Constraint::Length(9),
            Constraint::Min(10),
        ],
    )
    .header(grp_hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Group Differential Expression (Hellemans 2007 qBase vs. Pfaffl 2001 vs. Livak 2001) "),
    );
    f.render_widget(grp_tbl, split[0]);

    let samp_hdr = Row::new([
        "Sample",
        "Group",
        "Target",
        "Mean Cq ± SD (n)",
        "qBase NRQ",
        "log2(NRQ)",
        "Pfaffl Ratio",
        "2^-ddCq",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );

    let samp_rows = state
        .analysis
        .sample_nrqs
        .iter()
        .skip(state.scroll_offset)
        .map(|s| {
            Row::new([
                s.sample_id.clone(),
                s.condition_group.clone(),
                s.target_gene.clone(),
                format!("{:.2} ± {:.3} ({})", s.mean_cq, s.tech_sd_cq, s.n_tech_reps),
                format!("{:.3}", s.nrq_qbase),
                format!("{:+.3}", s.log2_nrq),
                format!("{:.3}", s.pfaffl_single_ref),
                format!("{:.3}", s.livak_2ddcq),
            ])
        });

    let samp_tbl = Table::new(
        samp_rows,
        [
            Constraint::Length(10),
            Constraint::Length(15),
            Constraint::Length(9),
            Constraint::Length(19),
            Constraint::Length(11),
            Constraint::Length(11),
            Constraint::Length(13),
            Constraint::Min(9),
        ],
    )
    .header(samp_hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Biological Replicate Normalized Relative Quantities (NRQ) "),
    );
    f.render_widget(samp_tbl, split[1]);
}

fn draw_efficiency_genorm_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &QpcrTuiState) {
    let split = Layout::default()
        .direction(Direction::Horizontal)
        .constraints([Constraint::Percentage(58), Constraint::Percentage(42)])
        .split(area);

    let sc_hdr = Row::new(["Target", "n", "Slope", "Intercept", "R²", "Factor E", "% Eff", "MIQE"])
        .style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    let sc_rows = state.analysis.standard_curves.iter().map(|sc| {
        let st = if sc.miqe_pass {
            Style::default().fg(Color::Green)
        } else {
            Style::default().fg(Color::LightRed)
        };
        Row::new([
            sc.target_gene.clone(),
            format!("{}", sc.n_points),
            format!("{:.3}", sc.slope),
            format!("{:.2}", sc.intercept_cq),
            format!("{:.4}", sc.r_squared),
            format!("{:.3}", sc.amplification_factor_e),
            format!("{:.1}%", sc.efficiency_pct),
            String::from(if sc.miqe_pass { "PASS" } else { "FLAG" }),
        ])
        .style(st)
    });
    let sc_tbl = Table::new(
        sc_rows,
        [
            Constraint::Length(8),
            Constraint::Length(4),
            Constraint::Length(8),
            Constraint::Length(10),
            Constraint::Length(8),
            Constraint::Length(9),
            Constraint::Length(8),
            Constraint::Min(6),
        ],
    )
    .header(sc_hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Serial-Dilution Standard Curves (E = 10^(-1/slope)) "),
    );
    f.render_widget(sc_tbl, split[0]);

    let gn_hdr = Row::new(["Candidate Gene", "geNorm M", "Status (M <= 0.50)"])
        .style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );
    let gn_rows = state.analysis.genorm_rankings.iter().map(|g| {
        let st = if g.stable_pass {
            Style::default().fg(Color::Green)
        } else {
            Style::default()
                .fg(Color::LightRed)
                .add_modifier(Modifier::BOLD)
        };
        Row::new([
            g.gene_name.clone(),
            format!("{:.4}", g.stability_m),
            String::from(if g.stable_pass {
                "STABLE (PASS)"
            } else {
                "UNSTABLE (REJECT)"
            }),
        ])
        .style(st)
    });
    let gn_tbl = Table::new(
        gn_rows,
        [
            Constraint::Length(15),
            Constraint::Length(11),
            Constraint::Min(18),
        ],
    )
    .header(gn_hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Vandesompele et al. (2002) geNorm Stability "),
    );
    f.render_widget(gn_tbl, split[1]);
}

fn draw_melt_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &QpcrTuiState) {
    let hdr = Row::new([
        "Well",
        "Sample",
        "Target",
        "Primary Tm (°C)",
        "FWHM (°C)",
        "Peak -dF/dT",
        "# Peaks",
        "Wittwer (1997) Specificity Diagnosis",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );

    let rows = state.melt_results.iter().map(|m| {
        let has_dimer = m.peaks.iter().any(|p| p.is_primer_dimer_suspect);
        let st = if !has_dimer && m.peaks.len() == 1 {
            Style::default().fg(Color::Green)
        } else {
            Style::default()
                .fg(Color::LightRed)
                .add_modifier(Modifier::BOLD)
        };
        let (tm_s, fwhm_s, pk_s) = if let Some(p0) = m.peaks.first() {
            (
                format!("{:.2}", p0.tm_celsius),
                format!("{:.2}", p0.fwhm_celsius),
                format!("{:.1}", p0.peak_neg_df_dt),
            )
        } else {
            (String::from("-"), String::from("-"), String::from("-"))
        };
        Row::new([
            m.well_id.clone(),
            m.sample_id.clone(),
            m.target_gene.clone(),
            tm_s,
            fwhm_s,
            pk_s,
            format!("{}", m.peaks.len()),
            String::from(m.specificity_status),
        ])
        .style(st)
    });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(6),
            Constraint::Length(11),
            Constraint::Length(8),
            Constraint::Length(16),
            Constraint::Length(10),
            Constraint::Length(12),
            Constraint::Length(8),
            Constraint::Min(32),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" SYBR Green I Dissociation (-dF/dT) Amplicon & Primer-Dimer Audit "),
    );
    f.render_widget(tbl, area);
}

fn draw_miqe_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &QpcrTuiState) {
    let hdr = Row::new(["Status", "MIQE Category", "Subject", "Diagnostic Finding"])
        .style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );

    let rows = state
        .analysis
        .miqe_checks
        .iter()
        .skip(state.scroll_offset)
        .map(|chk| {
            let col = match chk.severity {
                MiqeFlagSeverity::Pass => Color::Green,
                MiqeFlagSeverity::Warning => Color::Yellow,
                MiqeFlagSeverity::Fail => Color::LightRed,
            };
            Row::new([
                Cell::from(chk.severity.label())
                    .style(Style::default().fg(col).add_modifier(Modifier::BOLD)),
                Cell::from(chk.category),
                Cell::from(chk.subject.clone()),
                Cell::from(chk.detail.clone()),
            ])
        });

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(7),
            Constraint::Length(20),
            Constraint::Length(18),
            Constraint::Min(40),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Bustin et al. (2009) Automated MIQE Compliance Checklist "),
    );
    f.render_widget(tbl, area);
}
