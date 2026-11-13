// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Interactive 4-tab `Ratatui` terminal epidemiology & biostatistics studio for Tool #5: `osabt epi`.

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

use crate::epi_models::{
    demo_outbreak_incidence, estimate_cori_rt, evaluate_diagnostic_cascade,
    rogan_gladen_prevalence, simulate_seir_rk4, CoriRtConfig, DiagnosticStepInput, SeirParams,
};
use crate::epi_power::{
    design_cohort_study, design_test_negative_study, render_power_curve_svg,
    render_rt_epidemic_svg, render_seir_trajectory_svg,
};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum EpiTab {
    RenewalRt,
    SeirSweep,
    DiagnosticBayes,
    StudyPower,
}

struct EpiTuiState {
    active_tab: EpiTab,
    incidence: Vec<f64>,
    rt_cfg: CoriRtConfig,
    seir_params: SeirParams,
    prior_prev: f64,
    se1: f64,
    sp1: f64,
    pos1: bool,
    se2: f64,
    sp2: f64,
    pos2: bool,
    tnd_mode: bool,
    p0_control: f64,
    ve_or_rr: f64,
    ratio_k: f64,
    target_power: f64,
    status_msg: String,
}

impl EpiTuiState {
    fn new(incidence_csv: Option<PathBuf>) -> Result<Self, String> {
        let incidence = match incidence_csv {
            Some(p) => {
                let raw = fs::read_to_string(&p)
                    .map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
                parse_incidence_series(&raw)?
            }
            None => demo_outbreak_incidence(),
        };

        Ok(Self {
            active_tab: EpiTab::RenewalRt,
            incidence,
            rt_cfg: CoriRtConfig::default(),
            seir_params: SeirParams::default(),
            prior_prev: 0.08,
            se1: 0.86,
            sp1: 0.985,
            pos1: true,
            se2: 0.975,
            sp2: 0.995,
            pos2: true,
            tnd_mode: true,
            p0_control: 0.52,
            ve_or_rr: 0.60,
            ratio_k: 2.0,
            target_power: 0.80,
            status_msg: String::from(
                "1-4/Tab: switch tool | Left/Right: adjust primary parameter | Up/Down: secondary parameter | m: toggle TND/Cohort | e: export SVG | q: quit",
            ),
        })
    }

    fn export_active_svg(&mut self) {
        match self.active_tab {
            EpiTab::RenewalRt => {
                if let Ok(pts) = estimate_cori_rt(&self.incidence, &self.rt_cfg) {
                    let svg = render_rt_epidemic_svg("Outbreak Surveillance", &pts);
                    if fs::write("epi_rt_curve.svg", svg).is_ok() {
                        self.status_msg = String::from("Exported ./epi_rt_curve.svg");
                    }
                }
            }
            EpiTab::SeirSweep => {
                if let Ok(sim) = simulate_seir_rk4(&self.seir_params) {
                    let svg = render_seir_trajectory_svg(&sim);
                    if fs::write("epi_seir_trajectory.svg", svg).is_ok() {
                        self.status_msg = String::from("Exported ./epi_seir_trajectory.svg");
                    }
                }
            }
            EpiTab::DiagnosticBayes => {
                self.status_msg = String::from(
                    "Diagnostic Bayes summary displayed live; switch to Rt, SEIR, or Power tab to export .svg curves.",
                );
            }
            EpiTab::StudyPower => {
                let plan_res = if self.tnd_mode {
                    design_test_negative_study(
                        self.p0_control,
                        self.ve_or_rr.clamp(0.05, 0.95),
                        self.ratio_k,
                        0.05,
                        self.target_power,
                        Some(0.25),
                    )
                } else {
                    design_cohort_study(
                        self.p0_control.clamp(0.01, 0.45),
                        (1.0 - self.ve_or_rr).clamp(0.10, 0.95),
                        self.ratio_k,
                        0.05,
                        self.target_power,
                    )
                };
                if let Ok(plan) = plan_res {
                    let svg = render_power_curve_svg(&plan);
                    if fs::write("epi_power_curve.svg", svg).is_ok() {
                        self.status_msg = String::from("Exported ./epi_power_curve.svg");
                    }
                }
            }
        }
    }
}

/// Parses daily case counts from whitespace, comma, or multi-column CSV text (uses last numeric column per row).
pub fn parse_incidence_series(raw: &str) -> Result<Vec<f64>, String> {
    let mut vals = Vec::new();
    for line in raw.lines() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let tokens: Vec<&str> = t.split([',', '\t', ';', ' ']).filter(|s| !s.is_empty()).collect();
        if let Some(last) = tokens.last() {
            if let Ok(v) = last.parse::<f64>() {
                if v >= 0.0 && v.is_finite() {
                    vals.push(v);
                }
            }
        }
    }
    if vals.len() < 3 {
        return Err(String::from(
            "Could not parse >= 3 non-negative daily incidence counts from input file.",
        ));
    }
    Ok(vals)
}

/// Launches the interactive 4-tab `Ratatui` terminal workbench for `osabt epi`.
pub fn run_epi_tui(incidence_csv: Option<PathBuf>) -> Result<(), String> {
    let mut state = EpiTuiState::new(incidence_csv)?;

    enable_raw_mode().map_err(|e| format!("Failed enabling terminal raw mode: {e}"))?;
    let mut stdout = io::stdout();
    execute!(stdout, EnterAlternateScreen)
        .map_err(|e| format!("Failed entering alternate screen: {e}"))?;
    let backend = CrosstermBackend::new(stdout);
    let mut terminal =
        Terminal::new(backend).map_err(|e| format!("Failed creating Ratatui terminal: {e}"))?;

    let res = epi_loop(&mut terminal, &mut state);

    let _ = disable_raw_mode();
    let _ = execute!(terminal.backend_mut(), LeaveAlternateScreen);
    let _ = terminal.show_cursor();
    res
}

fn epi_loop(
    terminal: &mut Terminal<CrosstermBackend<io::Stdout>>,
    state: &mut EpiTuiState,
) -> Result<(), String> {
    loop {
        terminal
            .draw(|f| draw_epi_ui(f, state))
            .map_err(|e| format!("Draw error: {e}"))?;

        if event::poll(Duration::from_millis(120)).map_err(|e| e.to_string())? {
            if let Event::Key(key) = event::read().map_err(|e| e.to_string())? {
                if key.kind != KeyEventKind::Press {
                    continue;
                }
                match key.code {
                    KeyCode::Char('q') | KeyCode::Esc => return Ok(()),
                    KeyCode::Char('1') => state.active_tab = EpiTab::RenewalRt,
                    KeyCode::Char('2') => state.active_tab = EpiTab::SeirSweep,
                    KeyCode::Char('3') => state.active_tab = EpiTab::DiagnosticBayes,
                    KeyCode::Char('4') => state.active_tab = EpiTab::StudyPower,
                    KeyCode::Tab => {
                        state.active_tab = match state.active_tab {
                            EpiTab::RenewalRt => EpiTab::SeirSweep,
                            EpiTab::SeirSweep => EpiTab::DiagnosticBayes,
                            EpiTab::DiagnosticBayes => EpiTab::StudyPower,
                            EpiTab::StudyPower => EpiTab::RenewalRt,
                        };
                    }
                    KeyCode::Char('e') => state.export_active_svg(),
                    KeyCode::Char('m') => {
                        if state.active_tab == EpiTab::DiagnosticBayes {
                            state.pos2 = !state.pos2;
                        } else if state.active_tab == EpiTab::StudyPower {
                            state.tnd_mode = !state.tnd_mode;
                        }
                    }
                    KeyCode::Right => match state.active_tab {
                        EpiTab::RenewalRt => {
                            state.rt_cfg.serial_interval_mean =
                                (state.rt_cfg.serial_interval_mean + 0.3).min(15.0);
                        }
                        EpiTab::SeirSweep => {
                            state.seir_params.r0 = (state.seir_params.r0 + 0.15).min(12.0);
                        }
                        EpiTab::DiagnosticBayes => {
                            state.prior_prev = (state.prior_prev + 0.02).min(0.95);
                        }
                        EpiTab::StudyPower => {
                            state.ve_or_rr = (state.ve_or_rr + 0.05).min(0.90);
                        }
                    },
                    KeyCode::Left => match state.active_tab {
                        EpiTab::RenewalRt => {
                            state.rt_cfg.serial_interval_mean =
                                (state.rt_cfg.serial_interval_mean - 0.3).max(1.0);
                        }
                        EpiTab::SeirSweep => {
                            state.seir_params.r0 = (state.seir_params.r0 - 0.15).max(0.4);
                        }
                        EpiTab::DiagnosticBayes => {
                            state.prior_prev = (state.prior_prev - 0.02).max(0.005);
                        }
                        EpiTab::StudyPower => {
                            state.ve_or_rr = (state.ve_or_rr - 0.05).max(0.10);
                        }
                    },
                    KeyCode::Up => match state.active_tab {
                        EpiTab::RenewalRt => {
                            state.rt_cfg.window_days = (state.rt_cfg.window_days + 1).min(14);
                        }
                        EpiTab::SeirSweep => {
                            state.seir_params.npi_reduction_fraction =
                                (state.seir_params.npi_reduction_fraction + 0.05).min(0.90);
                        }
                        EpiTab::DiagnosticBayes => {
                            state.pos1 = !state.pos1;
                        }
                        EpiTab::StudyPower => {
                            state.ratio_k = (state.ratio_k + 0.5).min(5.0);
                        }
                    },
                    KeyCode::Down => match state.active_tab {
                        EpiTab::RenewalRt => {
                            state.rt_cfg.window_days = state.rt_cfg.window_days.saturating_sub(1).max(2);
                        }
                        EpiTab::SeirSweep => {
                            state.seir_params.npi_reduction_fraction =
                                (state.seir_params.npi_reduction_fraction - 0.05).max(0.0);
                        }
                        EpiTab::DiagnosticBayes => {
                            state.pos1 = !state.pos1;
                        }
                        EpiTab::StudyPower => {
                            state.ratio_k = (state.ratio_k - 0.5).max(0.5);
                        }
                    },
                    _ => {}
                }
            }
        }
    }
}

fn draw_epi_ui(f: &mut ratatui::Frame<'_>, state: &EpiTuiState) {
    let chunks = Layout::default()
        .direction(Direction::Vertical)
        .constraints([
            Constraint::Length(3),
            Constraint::Min(12),
            Constraint::Length(3),
        ])
        .split(f.area());

    let titles = [
        "[1] Cori Renewal Rt",
        "[2] RK4 SIR / SEIR",
        "[3] Diagnostic Bayes & Rogan-Gladen",
        "[4] Cohort & Test-Negative Power",
    ]
    .iter()
    .map(|t| Line::from(*t))
    .collect::<Vec<_>>();

    let idx = match state.active_tab {
        EpiTab::RenewalRt => 0,
        EpiTab::SeirSweep => 1,
        EpiTab::DiagnosticBayes => 2,
        EpiTab::StudyPower => 3,
    };

    let tabs = Tabs::new(titles)
        .select(idx)
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" OSABT Epidemiology & Study-Design Power Studio "),
        )
        .highlight_style(
            Style::default()
                .fg(Color::Cyan)
                .add_modifier(Modifier::BOLD),
        );
    f.render_widget(tabs, chunks[0]);

    match state.active_tab {
        EpiTab::RenewalRt => draw_rt_tab(f, chunks[1], state),
        EpiTab::SeirSweep => draw_seir_tab(f, chunks[1], state),
        EpiTab::DiagnosticBayes => draw_bayes_tab(f, chunks[1], state),
        EpiTab::StudyPower => draw_power_tab(f, chunks[1], state),
    }

    let footer = Paragraph::new(state.status_msg.as_str())
        .block(Block::default().borders(Borders::ALL).title(" Controls "));
    f.render_widget(footer, chunks[2]);
}

fn draw_rt_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EpiTuiState) {
    let pts = estimate_cori_rt(&state.incidence, &state.rt_cfg).unwrap_or_default();
    let header = Row::new(["Day", "Cases I_t", "Lambda_t", "Mean R_t", "95% CrI", "Pr(Rt>1)"])
        .style(
            Style::default()
                .fg(Color::Yellow)
                .add_modifier(Modifier::BOLD),
        );

    let rows = pts.iter().skip(state.rt_cfg.window_days.saturating_sub(1)).take(18).map(|p| {
        let rt_str = p.rt_mean.map(|v| format!("{v:.2}")).unwrap_or_else(|| "-".into());
        let cri_str = match (p.ci95_low, p.ci95_high) {
            (Some(l), Some(h)) => format!("[{l:.2}, {h:.2}]"),
            _ => "-".into(),
        };
        let pr_str = p
            .prob_rt_gt_1
            .map(|pr| format!("{:.1}%", pr * 100.0))
            .unwrap_or_else(|| "-".into());
        let col = match p.rt_mean {
            Some(v) if v >= 1.0 => Color::Red,
            Some(_) => Color::Green,
            None => Color::Gray,
        };
        Row::new(vec![
            Cell::from(format!("d{:02}", p.day_index + 1)),
            Cell::from(format!("{:.0}", p.incidence)),
            Cell::from(format!("{:.2}", p.total_infectiousness_lambda)),
            Cell::from(Span::styled(rt_str, Style::default().fg(col).add_modifier(Modifier::BOLD))),
            Cell::from(cri_str),
            Cell::from(pr_str),
        ])
    });

    let widths = [
        Constraint::Length(6),
        Constraint::Length(10),
        Constraint::Length(11),
        Constraint::Length(10),
        Constraint::Length(16),
        Constraint::Length(10),
    ];

    let title = format!(
        " Cori et al. (2013) Renewal R_t • SI Mean={:.1}d (Left/Right) | Window tau={}d (Up/Down) ",
        state.rt_cfg.serial_interval_mean, state.rt_cfg.window_days
    );
    let table = Table::new(rows, widths)
        .header(header)
        .block(Block::default().borders(Borders::ALL).title(title));
    f.render_widget(table, area);
}

fn draw_seir_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EpiTuiState) {
    let mut lines = Vec::new();
    if let Ok(sim) = simulate_seir_rk4(&state.seir_params) {
        lines.push(Line::from(Span::styled(
            format!(
                "Parameters: R0 = {:.2} (Left/Right) | NPI Reduction = {:.0}% at day 35 (Up/Down) | Baseline Immune = {:.0}%",
                sim.params.r0,
                sim.params.npi_reduction_fraction * 100.0,
                sim.params.initial_immune_fraction * 100.0
            ),
            Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(""));
        lines.push(Line::from(format!(
            "• Herd Immunity Threshold (1 - 1/R0) : {:.1}%",
            sim.herd_immunity_threshold * 100.0
        )));
        lines.push(Line::from(format!(
            "• Initial Effective Re(0)            : {:.2}",
            sim.initial_effective_re
        )));
        lines.push(Line::from(format!(
            "• Peak Infectious Prevalence         : {:.0} persons on Day {:.0}",
            sim.peak_infectious, sim.peak_day
        )));
        lines.push(Line::from(format!(
            "• Cumulative Epidemic Attack Rate    : {:.1}% ({:.0} / {:.0})",
            sim.final_attack_fraction * 100.0,
            sim.total_infected_persons,
            sim.params.population
        )));
        lines.push(Line::from(""));
        lines.push(Line::from("Sampled Trajectory (every 20 days):"));
        for tp in sim.trajectory.iter().step_by(20) {
            lines.push(Line::from(format!(
                "  Day {:>3.0} | S={:>7.0}  E={:>6.0}  I={:>6.0}  R={:>7.0} | Re={:.2}",
                tp.day, tp.susceptible, tp.exposed, tp.infectious, tp.recovered, tp.effective_re
            )));
        }
    }
    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" RK4 Deterministic SIR / SEIR Simulator "),
        );
    f.render_widget(p, area);
}

fn draw_bayes_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EpiTuiState) {
    let steps = vec![
        DiagnosticStepInput {
            assay_name: String::from("Assay 1 (Rapid Antigen)"),
            sensitivity: state.se1,
            specificity: state.sp1,
            observed_positive: state.pos1,
        },
        DiagnosticStepInput {
            assay_name: String::from("Assay 2 (Confirmatory RT-qPCR)"),
            sensitivity: state.se2,
            specificity: state.sp2,
            observed_positive: state.pos2,
        },
    ];

    let mut lines = vec![
        Line::from(Span::styled(
            format!(
                "Pre-Test Prevalence = {:.1}% (Left/Right) | Up/Down toggles Assay 1 | [m] toggles Assay 2",
                state.prior_prev * 100.0
            ),
            Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
        )),
        Line::from(""),
    ];

    if let Ok(cascade) = evaluate_diagnostic_cascade(state.prior_prev, &steps) {
        for s in &cascade {
            lines.push(Line::from(format!(
                "• {} [{}] | LR+={:.1}, LR-={:.3} | Pre {:.2}% -> Post {:.2}% (Single-stage PPV {:.1}%, NPV {:.2}%)",
                s.assay_name,
                if s.observed_positive { "POS +" } else { "NEG -" },
                s.lr_plus,
                s.lr_minus,
                s.pre_test_prob * 100.0,
                s.post_test_prob * 100.0,
                s.stage_ppv * 100.0,
                s.stage_npv * 100.0
            )));
        }
    }

    if let Ok(rg) = rogan_gladen_prevalence(120, 1000, state.se1, state.sp1) {
        lines.push(Line::from(""));
        lines.push(Line::from(Span::styled(
            "Rogan & Gladen (1978) Field Serosurvey Example (120 / 1000 positive on Assay 1):",
            Style::default().fg(Color::Cyan).add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(format!(
            "  Apparent Prevalence : {:.2}% (95% Wilson CI {:.2}% .. {:.2}%)",
            rg.apparent_prevalence * 100.0,
            rg.apparent_ci95.0 * 100.0,
            rg.apparent_ci95.1 * 100.0
        )));
        lines.push(Line::from(format!(
            "  True Prevalence     : {:.2}% (95% RG CI {:.2}% .. {:.2}% | Youden J = {:.3})",
            rg.true_prevalence * 100.0,
            rg.true_ci95_clamped.0 * 100.0,
            rg.true_ci95_clamped.1 * 100.0,
            rg.youden_j
        )));
    }

    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Sequential Diagnostic Bayes & Rogan-Gladen Prevalence "),
        );
    f.render_widget(p, area);
}

fn draw_power_tab(f: &mut ratatui::Frame<'_>, area: Rect, state: &EpiTuiState) {
    let plan_res = if state.tnd_mode {
        design_test_negative_study(
            state.p0_control,
            state.ve_or_rr.clamp(0.05, 0.95),
            state.ratio_k,
            0.05,
            state.target_power,
            Some(0.25),
        )
    } else {
        design_cohort_study(
            state.p0_control.clamp(0.01, 0.45),
            (1.0 - state.ve_or_rr).clamp(0.10, 0.95),
            state.ratio_k,
            0.05,
            state.target_power,
        )
    };

    let mut lines = Vec::new();
    if let Ok(plan) = plan_res {
        lines.push(Line::from(Span::styled(
            format!(
                "{} • {} | Allocation Ratio k={:.1} (Up/Down) | [m] Toggle TND/Cohort",
                plan.kind.title(),
                plan.effect_label,
                plan.allocation_ratio_k
            ),
            Style::default().fg(Color::Yellow).add_modifier(Modifier::BOLD),
        )));
        lines.push(Line::from(""));
        lines.push(Line::from(format!(
            "• Group 0 Proportion (Control / Unexposed) : {:.2}%",
            plan.p_group_0 * 100.0
        )));
        lines.push(Line::from(format!(
            "• Group 1 Proportion (Case / Exposed)      : {:.2}%",
            plan.p_group_1 * 100.0
        )));
        lines.push(Line::from(format!(
            "• Uncorrected Sample Size                  : n1 = {} | n0 = {} (Total = {})",
            plan.n1_uncorrected,
            plan.n0_uncorrected,
            plan.n1_uncorrected + plan.n0_uncorrected
        )));
        lines.push(Line::from(Span::styled(
            format!(
                "• Fleiss (1980) Continuity-Corrected N     : n1 = {} | n0 = {} (Total = {})",
                plan.n1_fleiss_cc, plan.n0_fleiss_cc, plan.total_fleiss_cc
            ),
            Style::default().fg(Color::Green).add_modifier(Modifier::BOLD),
        )));
        if let Some(swabs) = plan.total_ari_swabs_needed {
            lines.push(Line::from(format!(
                "• Total Symptomatic ARI Swabs (@25% pos)   : {swabs} clinic encounters"
            )));
        }
        lines.push(Line::from(""));
        lines.push(Line::from("Power Curve Checkpoints (press [e] to export .svg):"));
        for &(tot_n, pwr) in plan.power_curve.iter().step_by(5) {
            lines.push(Line::from(format!(
                "  Total N = {:>5} -> Power = {:>5.1}%",
                tot_n,
                pwr * 100.0
            )));
        }
    }

    let p = Paragraph::new(lines)
        .wrap(Wrap { trim: false })
        .block(
            Block::default()
                .borders(Borders::ALL)
                .title(" Fleiss Cohort & Jackson-Nelson Test-Negative Design Power "),
        );
    f.render_widget(p, area);
}
