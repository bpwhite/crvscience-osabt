// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Interactive 4-tab `Ratatui` terminal Wildlife Capture-Recapture & Imperfect-Detection
//! Occupancy Studio (`osabt wildlife`).

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
use ratatui::text::Line;
use ratatui::widgets::{Block, Borders, Cell, Paragraph, Row, Table, Tabs, Wrap};
use ratatui::Terminal;

use crate::wildlife_markrecap::{
    analyze_closed_capture_histories, calc_chapman_lincoln_petersen, demo_pnw_capture_history_csv,
    estimate_cjs_apparent_survival, parse_capture_history_csv, render_wildlife_population_svg,
    ChapmanPetersenEstimate, CjsSurvivalReport, ClosedCaptureSuiteReport,
};
use crate::wildlife_occupancy::{
    demo_pnw_occupancy_csv, fit_mackenzie_occupancy_suite, parse_occupancy_csv,
    plan_occupancy_survey_effort, OccupancySuiteReport, OccupancySurveyPlan,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum WildlifeTab {
    OccupancyMle,
    ClosedMarkRecap,
    OpenCjsSurvival,
    SurveyEffortPlanner,
}

struct WildlifeTuiState {
    active_tab: WildlifeTab,
    occ_report: OccupancySuiteReport,
    closed_report: ClosedCaptureSuiteReport,
    cjs_report: CjsSurvivalReport,
    chapman_demo: ChapmanPetersenEstimate,
    target_se: f64,
    cost_ratio: f64,
    survey_plan: OccupancySurveyPlan,
    scroll_offset: usize,
    status_msg: String,
}

impl WildlifeTuiState {
    fn new(occ_path: Option<PathBuf>, cap_path: Option<PathBuf>) -> Result<Self, String> {
        let occ_raw = if let Some(op) = occ_path {
            fs::read_to_string(&op)
                .map_err(|e| format!("Failed reading occupancy CSV '{}': {e}", op.display()))?
        } else {
            demo_pnw_occupancy_csv().to_string()
        };
        let occ_rows = parse_occupancy_csv(&occ_raw)?;
        let occ_report = fit_mackenzie_occupancy_suite(&occ_rows)?;

        let cap_raw = if let Some(cp) = cap_path {
            fs::read_to_string(&cp)
                .map_err(|e| format!("Failed reading capture-history CSV '{}': {e}", cp.display()))?
        } else {
            demo_pnw_capture_history_csv().to_string()
        };
        let (_, histories) = parse_capture_history_csv(&cap_raw)?;
        let closed_report = analyze_closed_capture_histories(&histories)?;
        let cjs_report = estimate_cjs_apparent_survival(&histories)?;
        let chapman_demo = calc_chapman_lincoln_petersen(40, 35, 14)?;

        let best = &occ_report.models[occ_report.best_model_index];
        let target_se = 0.06_f64;
        let cost_ratio = 2.0_f64;
        let survey_plan = plan_occupancy_survey_effort(
            best.psi_hat.clamp(0.05, 0.95),
            best.p_hat.clamp(0.05, 0.95),
            0.05,
            target_se,
            cost_ratio,
        )?;

        Ok(Self {
            active_tab: WildlifeTab::OccupancyMle,
            occ_report,
            closed_report,
            cjs_report,
            chapman_demo,
            target_se,
            cost_ratio,
            survey_plan,
            scroll_offset: 0,
            status_msg: String::from(
                "Keys: [1-4/Tab] Views | [Up/Down] Scroll Sites | [+/-] Target SE(psi) | [[/]] Site/Visit Cost | [e] Export SVG | [q] Quit",
            ),
        })
    }

    fn recompute_plan(&mut self) {
        let best = &self.occ_report.models[self.occ_report.best_model_index];
        if let Ok(p) = plan_occupancy_survey_effort(
            best.psi_hat.clamp(0.05, 0.95),
            best.p_hat.clamp(0.05, 0.95),
            0.05,
            self.target_se,
            self.cost_ratio,
        ) {
            self.survey_plan = p;
        }
    }
}

/// Launches the interactive 4-tab `Ratatui` terminal studio for `osabt wildlife`.
pub fn run_wildlife_tui(
    occ_path: Option<PathBuf>,
    cap_path: Option<PathBuf>,
) -> Result<(), String> {
    let mut state = WildlifeTuiState::new(occ_path, cap_path)?;

    enable_raw_mode().map_err(|e| format!("Failed enabling terminal raw mode: {e}"))?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)
        .map_err(|e| format!("Failed entering alternate screen: {e}"))?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal =
        Terminal::new(backend).map_err(|e| format!("Failed creating Ratatui terminal: {e}"))?;

    let res = run_loop(&mut terminal, &mut state);

    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen);
    let _ = terminal.show_cursor();

    res
}

fn run_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut WildlifeTuiState,
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
                    KeyCode::Char('1') => state.active_tab = WildlifeTab::OccupancyMle,
                    KeyCode::Char('2') => state.active_tab = WildlifeTab::ClosedMarkRecap,
                    KeyCode::Char('3') => state.active_tab = WildlifeTab::OpenCjsSurvival,
                    KeyCode::Char('4') => state.active_tab = WildlifeTab::SurveyEffortPlanner,
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            WildlifeTab::OccupancyMle => WildlifeTab::ClosedMarkRecap,
                            WildlifeTab::ClosedMarkRecap => WildlifeTab::OpenCjsSurvival,
                            WildlifeTab::OpenCjsSurvival => WildlifeTab::SurveyEffortPlanner,
                            WildlifeTab::SurveyEffortPlanner => WildlifeTab::OccupancyMle,
                        };
                    }
                    KeyCode::Down | KeyCode::Char('j') => {
                        state.scroll_offset = state.scroll_offset.saturating_add(1);
                    }
                    KeyCode::Up | KeyCode::Char('k') => {
                        state.scroll_offset = state.scroll_offset.saturating_sub(1);
                    }
                    KeyCode::Char('+') | KeyCode::Char('=') => {
                        state.target_se = (state.target_se + 0.01).min(0.15);
                        state.recompute_plan();
                        state.status_msg = format!(
                            "Updated target SE(psi) = {:.2} -> S* = {} sites x K* = {} visits",
                            state.target_se,
                            state.survey_plan.required_sites_for_target_se,
                            state.survey_plan.optimal_k_for_precision
                        );
                    }
                    KeyCode::Char('-') => {
                        state.target_se = (state.target_se - 0.01).max(0.02);
                        state.recompute_plan();
                        state.status_msg = format!(
                            "Updated target SE(psi) = {:.2} -> S* = {} sites x K* = {} visits",
                            state.target_se,
                            state.survey_plan.required_sites_for_target_se,
                            state.survey_plan.optimal_k_for_precision
                        );
                    }
                    KeyCode::Char(']') => {
                        state.cost_ratio = (state.cost_ratio + 0.5).min(20.0);
                        state.recompute_plan();
                        state.status_msg = format!(
                            "Updated Site/Visit Cost Ratio = {:.1} -> Optimal K* = {}",
                            state.cost_ratio, state.survey_plan.optimal_k_for_precision
                        );
                    }
                    KeyCode::Char('[') => {
                        state.cost_ratio = (state.cost_ratio - 0.5).max(0.0);
                        state.recompute_plan();
                        state.status_msg = format!(
                            "Updated Site/Visit Cost Ratio = {:.1} -> Optimal K* = {}",
                            state.cost_ratio, state.survey_plan.optimal_k_for_precision
                        );
                    }
                    KeyCode::Char('e') => {
                        let svg = render_wildlife_population_svg(
                            &state.occ_report,
                            &state.closed_report,
                            "PNW Imperfect-Detection Occupancy & Capture-Recapture",
                        );
                        if fs::write("wildlife_population.svg", svg).is_ok() {
                            state.status_msg = String::from("Exported wildlife_population.svg!");
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

fn draw_ui(f: &mut ratatui::Frame<'_>, state: &WildlifeTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] MacKenzie Occupancy & AICc",
        "[2] Closed Mark-Recapture (M0/Mt/Mh)",
        "[3] Open CJS Survival (Phi & p)",
        "[4] Survey Effort Planner (K* & S*)",
    ]
    .into_iter()
    .map(Line::from)
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        WildlifeTab::OccupancyMle => 0,
        WildlifeTab::ClosedMarkRecap => 1,
        WildlifeTab::OpenCjsSurvival => 2,
        WildlifeTab::SurveyEffortPlanner => 3,
    };

    let tabs = Tabs::new(titles)
        .block(Block::default().borders(Borders::ALL).title(format!(
            " OSABT Wildlife Studio (v1.5.0) — {} Sites (Naive={:.2}) | M_t+1={} Tagged Individuals ",
            state.occ_report.num_sites,
            state.occ_report.naive_occupancy,
            state.closed_report.distinct_individuals_mt1
        )))
        .select(idx)
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        WildlifeTab::OccupancyMle => draw_occupancy_tab(f, chunks[1], state),
        WildlifeTab::ClosedMarkRecap => draw_closed_tab(f, chunks[1], state),
        WildlifeTab::OpenCjsSurvival => draw_cjs_tab(f, chunks[1], state),
        WildlifeTab::SurveyEffortPlanner => draw_planner_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls & Status "))
        .style(Style::default().fg(Color::Gray));
    f.render_widget(footer, chunks[2]);
}

fn draw_occupancy_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &WildlifeTuiState) {
    let split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(8), Constraint::Min(8)])
        .split(area);

    let m_hdr = Row::new([
        "Model",
        "psi_hat (95% CI)",
        "p_hat (95% CI)",
        "Cov Slope (SE)",
        "LogLik",
        "AICc",
        "dAICc",
        "Weight",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );
    let m_rows = state.occ_report.models.iter().enumerate().map(|(i, m)| {
        let style = if i == state.occ_report.best_model_index {
            Style::default().fg(Color::Green).add_modifier(Modifier::BOLD)
        } else {
            Style::default()
        };
        Row::new([
            Cell::from(m.model.label()),
            Cell::from(format!(
                "{:.3} ({:.2}–{:.2})",
                m.psi_hat, m.psi_ci95.0, m.psi_ci95.1
            )),
            Cell::from(format!(
                "{:.3} ({:.2}–{:.2})",
                m.p_hat, m.p_ci95.0, m.p_ci95.1
            )),
            Cell::from(format!(
                "{:+.2} ({:.2})",
                m.covariate_logit_slope, m.covariate_slope_se
            )),
            Cell::from(format!("{:.2}", m.log_likelihood)),
            Cell::from(format!("{:.2}", m.aicc)),
            Cell::from(format!("{:.2}", m.delta_aicc)),
            Cell::from(format!("{:.1}%", m.akaike_weight * 100.0)),
        ])
        .style(style)
    });

    let m_tbl = Table::new(
        m_rows,
        [
            Constraint::Length(16),
            Constraint::Length(19),
            Constraint::Length(19),
            Constraint::Length(16),
            Constraint::Length(10),
            Constraint::Length(9),
            Constraint::Length(8),
            Constraint::Min(9),
        ],
    )
    .header(m_hdr)
    .block(Block::default().borders(Borders::ALL).title(format!(
        " MacKenzie et al. (2002) Occupancy AICc Table — Naive={}/{} ({:.1}%), Finite-Sample Est={:.2} Sites ",
        state.occ_report.naive_occupied_sites,
        state.occ_report.num_sites,
        state.occ_report.naive_occupancy * 100.0,
        state.occ_report.finite_sample_occupied_est
    )));
    f.render_widget(m_tbl, split[0]);

    let s_hdr = Row::new([
        "Site ID",
        "Visits K",
        "Detections y",
        "Covariate",
        "psi_i",
        "p_i",
        "Pr(z_i = 1 | history)",
    ])
    .style(
        Style::default()
            .fg(Color::Yellow)
            .add_modifier(Modifier::BOLD),
    );
    let s_rows = state
        .occ_report
        .site_posteriors
        .iter()
        .skip(state.scroll_offset)
        .map(|s| {
            Row::new([
                s.site_id.clone(),
                format!("{}", s.valid_visits_k),
                format!("{}", s.detections_y),
                format!("{:.2}", s.covariate),
                format!("{:.3}", s.fitted_psi_i),
                format!("{:.3}", s.fitted_p_i),
                format!("{:.3}", s.posterior_occupied_prob),
            ])
        });
    let s_tbl = Table::new(
        s_rows,
        [
            Constraint::Length(16),
            Constraint::Length(10),
            Constraint::Length(14),
            Constraint::Length(12),
            Constraint::Length(10),
            Constraint::Length(10),
            Constraint::Min(22),
        ],
    )
    .header(s_hdr)
    .block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Empirical Bayes Conditional Site Occupancy Posteriors "),
    );
    f.render_widget(s_tbl, split[1]);
}

fn draw_closed_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &WildlifeTuiState) {
    let c = &state.closed_report;
    let cp = &state.chapman_demo;
    let ests = [
        &c.m0_equal_catchability,
        &c.mt_darroch_time_varying,
        &c.mh_chao1987_heterogeneity,
        &c.schnabel_multicensus,
    ];
    let hdr = Row::new([
        "Closed-Population Estimator",
        "N_hat",
        "SE(N_hat)",
        "95% Confidence Interval",
        "Mean Capture p",
    ])
    .style(
        Style::default()
            .fg(Color::Cyan)
            .add_modifier(Modifier::BOLD),
    );
    let rows = ests.iter().map(|e| {
        Row::new([
            e.model_name.clone(),
            format!("{:.2}", e.n_hat),
            format!("{:.2}", e.se),
            format!("[{:.1}, {:.1}]", e.ci95.0, e.ci95.1),
            format!("{:.3}", e.mean_capture_p),
        ])
    });

    let split = Layout::default()
        .direction(Direction::Vertical)
        .constraints([Constraint::Length(9), Constraint::Min(6)])
        .split(area);

    let tbl = Table::new(
        rows,
        [
            Constraint::Length(42),
            Constraint::Length(10),
            Constraint::Length(11),
            Constraint::Length(24),
            Constraint::Min(15),
        ],
    )
    .header(hdr)
    .block(Block::default().borders(Borders::ALL).title(format!(
        " Otis et al. (1978) & Chao (1987) Closed Abundance — M_t+1 = {} Individuals across K = {} Occasions ",
        c.distinct_individuals_mt1, c.occasions_k
    )));
    f.render_widget(tbl, split[0]);

    let freqs = c
        .frequency_spectrum_f
        .iter()
        .enumerate()
        .map(|(i, &cnt)| format!("f_{}={}", i + 1, cnt))
        .collect::<Vec<_>>()
        .join("  ");
    let info = vec![
        Line::from(format!("  Capture-Frequency Spectrum : {freqs}")),
        Line::from(format!(
            "  Captures per Occasion n_j  : {:?}",
            c.captures_per_occasion
        )),
        Line::from(""),
        Line::from(format!(
            "  Reference 2-Sample Chapman (1951) Lincoln–Petersen (M={}, C={}, R={}): N_hat = {:.2} ± {:.2} (95% CI [{:.1}, {:.1}])",
            cp.marked_sample_1,
            cp.captured_sample_2,
            cp.recaptured_r,
            cp.chapman_n_hat,
            cp.seber_se,
            cp.ci95.0,
            cp.ci95.1
        )),
    ];
    let p_card = Paragraph::new(info).wrap(Wrap { trim: false }).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Capture-Frequency Diagnostics & Chapman–Petersen Reference "),
    );
    f.render_widget(p_card, split[1]);
}

fn draw_cjs_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &WildlifeTuiState) {
    let cjs = &state.cjs_report;
    let mut lines = vec![
        Line::from("Cormack–Jolly–Seber (Lebreton et al. 1992) Open-Population Survival & Recapture:"),
        Line::from(""),
        Line::from(format!(
            "  Apparent Survival   phi_hat = {:.3} ± {:.3}   (95% CI [{:.3}, {:.3}])",
            cjs.phi_hat, cjs.phi_se, cjs.phi_ci95.0, cjs.phi_ci95.1
        )),
        Line::from(format!(
            "  Recapture Prob      p_hat   = {:.3} ± {:.3}   (95% CI [{:.3}, {:.3}])",
            cjs.p_hat, cjs.p_se, cjs.p_ci95.0, cjs.p_ci95.1
        )),
        Line::from(format!(
            "  Releases R_tot = {} | First Recaptures = {} | LogLik = {:.2} | AICc = {:.2}",
            cjs.total_releases, cjs.total_recaptures, cjs.log_likelihood, cjs.aicc
        )),
        Line::from(""),
        Line::from("  Burnham Sufficient-Statistic m-Array (Release Occasion i -> First Recapture j):"),
    ];
    for (i, (rel, rec_vec, never)) in cjs.m_array.iter().enumerate() {
        lines.push(Line::from(format!(
            "    Occasion {}: Released R_{} = {:2} | Recaptured m_ij = {:?} | Never Seen Again = {}",
            i + 1,
            i + 1,
            rel,
            rec_vec,
            never
        )));
    }
    let card = Paragraph::new(lines).wrap(Wrap { trim: false }).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" Open-Population CJS Apparent Survival & Burnham m-Array "),
    );
    f.render_widget(card, area);
}

fn draw_planner_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &WildlifeTuiState) {
    let pl = &state.survey_plan;
    let lines = vec![
        Line::from("MacKenzie & Royle (2005) Imperfect-Detection Survey Effort & Power Allocation:"),
        Line::from(""),
        Line::from(format!(
            "  Assumed True Occupancy psi           : {:.3}",
            pl.assumed_psi
        )),
        Line::from(format!(
            "  Assumed Per-Visit Detection p        : {:.3}",
            pl.assumed_p
        )),
        Line::from(format!(
            "  Site-Setup / Single-Visit Cost Ratio : {:.1}x  (adjust with [ / ])",
            state.cost_ratio
        )),
        Line::from(""),
        Line::from(format!(
            "  1. Presence-Absence Certainty (False-Absence alpha <= {:.1}%): K_min = {} visits/site (p* = {:.1}%)",
            pl.target_false_absence_alpha * 100.0,
            pl.min_visits_for_presence_certainty,
            pl.cumulative_detection_at_k_min * 100.0
        )),
        Line::from(format!(
            "  2. Variance-Optimal Repeat Visits per Site                 : K* = {} visits/site",
            pl.optimal_k_for_precision
        )),
        Line::from(format!(
            "  3. Required Sites S* for SE(psi_hat) <= {:.3} (+/- keys)    : S* = {} sites ({} total visits, achieved SE = {:.4})",
            state.target_se,
            pl.required_sites_for_target_se,
            pl.total_survey_visits,
            pl.achieved_psi_se
        )),
    ];
    let card = Paragraph::new(lines).wrap(Wrap { trim: false }).block(
        Block::default()
            .borders(Borders::ALL)
            .title(" MacKenzie & Royle (2005) Survey Design & Effort Planner "),
    );
    f.render_widget(card, area);
}
