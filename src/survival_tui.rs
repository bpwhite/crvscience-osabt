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

//! Interactive 4-tab `Ratatui` terminal Kaplan–Meier, Log-Rank, Cox Proportional Hazards,
//! & Acute Ecotoxicology `LC50` Studio (`osabt survival`).

use std::fs;
use std::io;
use std::path::PathBuf;
use std::time::Duration;

use crossterm::event::{self, Event, KeyCode, KeyEventKind};
use ratatui::backend::CrosstermBackend;
use ratatui::layout::{Constraint, Direction, Layout, Rect};
use ratatui::style::{Color, Modifier, Style};
use ratatui::text::Line;
use ratatui::widgets::{Block, Borders, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::Terminal;

use crate::tui_util::run_tui_session;

use crate::survival_core::{
    compute_kaplan_meier_curves, compute_log_rank_test, demo_survival_cohort_csv,
    fit_cox_proportional_hazards, parse_survival_csv, CoxPhFit, KaplanMeierCurve,
    LogRankTestResult,
};
use crate::survival_tox::{
    demo_ecotox_lc50_csv, fit_ecotox_lc50, parse_ecotox_csv, render_survival_and_lc50_svg,
    DoseResponseLink, EcotoxBioassayPoint, EcotoxLc50Report,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum SurvivalTab {
    KaplanMeier,
    CoxPhRegression,
    EcotoxLc50,
    DoseResponseTable,
}

struct SurvivalTuiState {
    active_tab: SurvivalTab,
    curves: Vec<KaplanMeierCurve>,
    log_rank: Option<LogRankTestResult>,
    cox_cov: CoxPhFit,
    cox_grp: CoxPhFit,
    ecotox_pts: Vec<EcotoxBioassayPoint>,
    link: DoseResponseLink,
    lc50_report: EcotoxLc50Report,
    scroll_offset: usize,
    status_msg: String,
}

impl SurvivalTuiState {
    fn new(surv_path: Option<PathBuf>, tox_path: Option<PathBuf>) -> Result<Self, String> {
        let surv_raw = if let Some(sp) = surv_path {
            fs::read_to_string(&sp)
                .map_err(|e| format!("Failed reading survival CSV '{}': {e}", sp.display()))?
        } else {
            demo_survival_cohort_csv().to_string()
        };
        let obs = parse_survival_csv(&surv_raw)?;
        let curves = compute_kaplan_meier_curves(&obs);
        let log_rank = compute_log_rank_test(&obs).ok();
        let cox_cov = fit_cox_proportional_hazards(&obs, "Covariate", false)?;
        let cox_grp = fit_cox_proportional_hazards(&obs, "Group_Indicator", true)?;

        let tox_raw = if let Some(tp) = tox_path {
            fs::read_to_string(&tp)
                .map_err(|e| format!("Failed reading ecotox CSV '{}': {e}", tp.display()))?
        } else {
            demo_ecotox_lc50_csv().to_string()
        };
        let ecotox_pts = parse_ecotox_csv(&tox_raw)?;
        let link = DoseResponseLink::Probit;
        let lc50_report = fit_ecotox_lc50(&ecotox_pts, link)?;

        Ok(Self {
            active_tab: SurvivalTab::KaplanMeier,
            curves,
            log_rank,
            cox_cov,
            cox_grp,
            ecotox_pts,
            link,
            lc50_report,
            scroll_offset: 0,
            status_msg: String::from(
                "Keys: [1-4/Tab] Views | [Up/Down] Scroll | [l] Toggle Probit/Logit LC50 | [e] Export .svg | [q] Quit",
            ),
        })
    }
}

/// Launches the interactive 4-tab `Ratatui` terminal studio for `osabt survival`.
pub fn run_survival_tui(
    surv_path: Option<PathBuf>,
    tox_path: Option<PathBuf>,
) -> Result<(), String> {
    let mut state = SurvivalTuiState::new(surv_path, tox_path)?;
    run_tui_session(|terminal| run_loop(terminal, &mut state))
}

fn run_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut SurvivalTuiState,
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
                        state.active_tab = SurvivalTab::KaplanMeier;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('2') => {
                        state.active_tab = SurvivalTab::CoxPhRegression;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('3') => {
                        state.active_tab = SurvivalTab::EcotoxLc50;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Char('4') => {
                        state.active_tab = SurvivalTab::DoseResponseTable;
                        state.scroll_offset = 0;
                    }
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            SurvivalTab::KaplanMeier => SurvivalTab::CoxPhRegression,
                            SurvivalTab::CoxPhRegression => SurvivalTab::EcotoxLc50,
                            SurvivalTab::EcotoxLc50 => SurvivalTab::DoseResponseTable,
                            SurvivalTab::DoseResponseTable => SurvivalTab::KaplanMeier,
                        };
                        state.scroll_offset = 0;
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        state.scroll_offset = state.scroll_offset.saturating_add(1);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.scroll_offset = state.scroll_offset.saturating_sub(1);
                    }
                    KeyCode::Char('l') => {
                        state.link = match state.link {
                            DoseResponseLink::Probit => DoseResponseLink::Logit,
                            DoseResponseLink::Logit => DoseResponseLink::Probit,
                        };
                        if let Ok(rep) = fit_ecotox_lc50(&state.ecotox_pts, state.link) {
                            state.lc50_report = rep;
                            state.status_msg = format!(
                                "Switched to {} -> LC50 = {:.2} (95% CI [{:.2}, {:.2}])",
                                state.link.label(),
                                state.lc50_report.lc50,
                                state.lc50_report.lc50_ci95.0,
                                state.lc50_report.lc50_ci95.1
                            );
                        }
                    }
                    KeyCode::Char('e') => {
                        let svg = render_survival_and_lc50_svg(
                            &state.curves,
                            state.log_rank.as_ref(),
                            &state.ecotox_pts,
                            &state.lc50_report,
                            "OSABT Survival & Ecotoxicology",
                        );
                        if fs::write("survival_lc50.svg", svg).is_ok() {
                            state.status_msg = String::from("Exported survival_lc50.svg!");
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

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &SurvivalTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Kaplan-Meier & Log-Rank",
        "[2] Cox PH Regression (HR)",
        "[3] Ecotox LC50 (Probit/Logit)",
        "[4] Bioassay Dose Table",
    ]
    .into_iter()
    .map(Line::from)
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        SurvivalTab::KaplanMeier => 0,
        SurvivalTab::CoxPhRegression => 1,
        SurvivalTab::EcotoxLc50 => 2,
        SurvivalTab::DoseResponseTable => 3,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT Survival Studio (v1.9.0) — {} Cohort Strata | {} LC50 = {:.2} ",
            state.curves.len(),
            state.link.label(),
            state.lc50_report.lc50
        )))
        .select(idx)
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        SurvivalTab::KaplanMeier => draw_km_tab(f, chunks[1], state),
        SurvivalTab::CoxPhRegression => draw_cox_tab(f, chunks[1], state),
        SurvivalTab::EcotoxLc50 => draw_lc50_tab(f, chunks[1], state),
        SurvivalTab::DoseResponseTable => draw_dose_table_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls & Status "))
        .style(Style::default().fg(Color::Gray));
    f.render_widget(footer, chunks[2]);
}

fn draw_km_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &SurvivalTuiState) {
    let split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(6), Constraint::Min(8)])
        .split(area);

    let mut summary_lines = Vec::new();
    for c in &state.curves {
        let med_s = c
            .median_survival_time
            .map(|m| format!("{m:.1}"))
            .unwrap_or_else(|| String::from(">max_t"));
        let final_s = c.steps.last().map(|s| s.survival_s).unwrap_or(1.0);
        summary_lines.push(Line::from(format!(
            "  Group {:18} : n={:2} | Events={:2} | Median Survival={} | Final S(t)={:.3}",
            c.group, c.total_n, c.total_events, med_s, final_s
        )));
    }
    if let Some(lr) = &state.log_rank {
        summary_lines.push(Line::from(format!(
            "  Mantel (1966) Log-Rank ({} vs {}): Chi^2 = {:.3}, p = {:.5}",
            lr.group_a, lr.group_b, lr.chi_square, lr.p_value
        )));
    }
    let top_card = Paragraph::new(summary_lines).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Kaplan–Meier Stratum Summary & Mantel Log-Rank Test "),
    );
    f.render_widget(top_card, split[0]);

    let hdr = Row::new([
        "Group",
        "Time t",
        "At Risk n_i",
        "Events d_i",
        "Censored",
        "S(t)",
        "Greenwood SE",
        "95% Log-Log CI",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );

    let all_rows: Vec<Row<'_>> = state
        .curves
        .iter()
        .flat_map(|c| {
            c.steps.iter().map(move |st| {
                Row::new([
                    c.group.clone(),
                    format!("{:.1}", st.time),
                    format!("{}", st.n_at_risk),
                    format!("{}", st.n_events),
                    format!("{}", st.n_censored),
                    format!("{:.3}", st.survival_s),
                    format!("{:.3}", st.greenwood_se),
                    format!("[{:.3}, {:.3}]", st.ci95.0, st.ci95.1),
                ])
            })
        })
        .skip(state.scroll_offset)
        .collect();

    let tbl = Table::new(
        all_rows,
        [
            Constraint::Length(18),
            Constraint::Length(9),
            Constraint::Length(12),
            Constraint::Length(11),
            Constraint::Length(10),
            Constraint::Length(8),
            Constraint::Length(13),
            Constraint::Min(18),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Product-Limit Survival Schedule (Scroll with Up/Down) "),
    );
    f.render_widget(tbl, split[1]);
}

fn draw_cox_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &SurvivalTuiState) {
    let lines = vec![
        Line::from("Cox (1972) Proportional-Hazards Partial-Likelihood Regression (Breslow Ties):"),
        Line::from(""),
        Line::from(format!(
            "  1. Continuous Covariate Model ({}, n={}, events={}):",
            state.cox_cov.covariate_name, state.cox_cov.n_subjects, state.cox_cov.n_events
        )),
        Line::from(format!(
            "     beta = {:+.4} ± {:.4}  |  Hazard Ratio HR = {:.4} (95% CI [{:.4}, {:.4}])",
            state.cox_cov.beta,
            state.cox_cov.se_beta,
            state.cox_cov.hazard_ratio,
            state.cox_cov.hr_ci95.0,
            state.cox_cov.hr_ci95.1
        )),
        Line::from(format!(
            "     Wald z = {:+.3}  |  p = {:.6}  |  Partial LogLik: {:.2} (Null {:.2})",
            state.cox_cov.wald_z,
            state.cox_cov.p_value,
            state.cox_cov.partial_log_lik_fit,
            state.cox_cov.partial_log_lik_null
        )),
        Line::from(""),
        Line::from("  2. Binary Group Indicator Model:"),
        Line::from(format!(
            "     beta = {:+.4} ± {:.4}  |  Hazard Ratio HR = {:.3} (95% CI [{:.3}, {:.3}])",
            state.cox_grp.beta,
            state.cox_grp.se_beta,
            state.cox_grp.hazard_ratio,
            state.cox_grp.hr_ci95.0,
            state.cox_grp.hr_ci95.1
        )),
        Line::from(format!(
            "     Wald z = {:+.3}  |  p = {:.6}",
            state.cox_grp.wald_z, state.cox_grp.p_value
        )),
    ];
    let card = Paragraph::new(lines).wrap(Wrap { trim: false }).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Cox (1972) Proportional Hazards & Hazard Ratio Inference "),
    );
    f.render_widget(card, area);
}

fn draw_lc50_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &SurvivalTuiState) {
    let r = &state.lc50_report;
    let lines = vec![
        Line::from(format!(
            "Acute Ecotoxicology Dose-Response Model: {} (press [l] to toggle Probit/Logit)",
            r.link.label()
        )),
        Line::from(""),
        Line::from(format!("  Intercept alpha      : {:+.4}", r.intercept_alpha)),
        Line::from(format!("  Log10(C) Slope beta  : {:+.4}", r.slope_beta)),
        Line::from(format!("  LC10 (10% Mortality) : {:.3}", r.lc10)),
        Line::from(format!(
            "  LC50 (50% Mortality) : {:.3}  (95% CI [{:.3}, {:.3}])",
            r.lc50, r.lc50_ci95.0, r.lc50_ci95.1
        )),
        Line::from(format!("  LC90 (90% Mortality) : {:.3}", r.lc90)),
        Line::from(format!(
            "  Binomial LogLik = {:.2}  |  Pearson Chi^2 = {:.3}",
            r.log_likelihood, r.pearson_chi2
        )),
    ];
    let card = Paragraph::new(lines).wrap(Wrap { trim: false }).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Finney (1971) Probit & Logit Acute LC50 Estimator "),
    );
    f.render_widget(card, area);
}

fn draw_dose_table_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &SurvivalTuiState) {
    let hdr = Row::new(["Concentration", "Exposed n", "Dead k", "Observed Mortality %"])
        .style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );
    let rows = state.ecotox_pts.iter().map(|pt| {
        Row::new([
            format!("{:.2}", pt.concentration),
            format!("{}", pt.exposed_n),
            format!("{}", pt.dead_n),
            format!("{:.1}%", 100.0 * (pt.dead_n as f64) / (pt.exposed_n as f64)),
        ])
    });
    let tbl = Table::new(
        rows,
        [
            Constraint::Length(16),
            Constraint::Length(12),
            Constraint::Length(10),
            Constraint::Min(22),
        ],
    )
    .header(hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Acute Bioassay Concentration-Mortality Schedule "),
    );
    f.render_widget(tbl, area);
}
