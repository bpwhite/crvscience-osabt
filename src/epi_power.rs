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

//! Tool #5 (`osabt epi`) — Epidemiological study-design power & sample-size calculators
//! (Fleiss, Tytun & Ury 1980 continuity-corrected Cohort / RCT two-proportion design;
//! Unmatched Case–Control odds-ratio design; and Jackson & Nelson 2013 Test-Negative Design
//! `VE = 1 - OR` vaccine effectiveness calculator) plus standalone vector `.svg` chart exporters.

use std::fmt::Write as _;

use crate::epi_models::{norm_cdf, norm_ppf, RtEstimatePoint, SeirSimulationResult};
use crate::svg_util::{begin_svg_canvas, xml_escape};

/// Supported epidemiological observational and interventional study designs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum EpiStudyKind {
    CohortRelativeRisk,
    CaseControlOddsRatio,
    TestNegativeVaccineEffectiveness,
}

impl EpiStudyKind {
    pub fn title(self) -> &'static str {
        match self {
            Self::CohortRelativeRisk => "Two-Arm Cohort / RCT (Relative Risk)",
            Self::CaseControlOddsRatio => "Unmatched Case-Control (Odds Ratio)",
            Self::TestNegativeVaccineEffectiveness => {
                "Test-Negative Design (Vaccine Effectiveness VE = 1 - OR)"
            }
        }
    }
}

/// Complete sample-size & power-curve solution for a two-group binary-outcome study.
#[derive(Debug, Clone)]
pub struct StudyPowerPlan {
    pub kind: EpiStudyKind,
    /// Control / unexposed proportion `p_0` (risk in unexposed for Cohort; exposure/vaccine coverage in controls for Case-Control/TND).
    pub p_group_0: f64,
    /// Exposed / case proportion `p_1` (risk in exposed for Cohort; exposure/vaccine coverage in cases for Case-Control/TND).
    pub p_group_1: f64,
    /// Ratio `k = n_0 / n_1` (controls per case, or unexposed per exposed).
    pub allocation_ratio_k: f64,
    pub alpha_two_sided: f64,
    pub target_power: f64,
    /// Effect metric (`RR`, `OR`, or `VE = 1 - OR`).
    pub effect_label: String,
    /// Uncorrected asymptotic normal sample size in Group 1 (`n_exposed` or `n_cases`).
    pub n1_uncorrected: usize,
    /// Uncorrected asymptotic normal sample size in Group 0 (`n_unexposed` or `n_controls`).
    pub n0_uncorrected: usize,
    /// Fleiss, Tytun & Ury (1980) continuity-corrected sample size in Group 1.
    pub n1_fleiss_cc: usize,
    /// Fleiss, Tytun & Ury (1980) continuity-corrected sample size in Group 0.
    pub n0_fleiss_cc: usize,
    pub total_fleiss_cc: usize,
    /// For Test-Negative Designs: estimated total symptomatic patients to swab given expected pathogen positivity.
    pub total_ari_swabs_needed: Option<usize>,
    /// `(total_sample_size, statistical_power_0_to_1)` curve across `0.2x .. 2.0x` required N.
    pub power_curve: Vec<(usize, f64)>,
}

/// Core Fleiss, Tytun & Ury (1980) sample-size solver comparing independent proportions `p_0` and `p_1`
/// with group size ratio `k = n_0 / n_1`, two-sided Type I error `alpha`, and desired `power`.
pub fn solve_two_proportion_sample_size(
    kind: EpiStudyKind,
    p_0: f64,
    p_1: f64,
    ratio_k: f64,
    alpha_two_sided: f64,
    target_power: f64,
    effect_label: String,
    ari_positivity: Option<f64>,
) -> Result<StudyPowerPlan, String> {
    if !(0.0001..=0.9999).contains(&p_0) || !(0.0001..=0.9999).contains(&p_1) {
        return Err(String::from("Group proportions p0 and p1 must lie strictly within (0, 1)."));
    }
    if (p_0 - p_1).abs() < 1e-6 {
        return Err(String::from(
            "Group proportions p0 and p1 are identical (null effect); sample size is infinite.",
        ));
    }
    if ratio_k <= 0.0 || !(1e-6..0.5).contains(&alpha_two_sided) || !(0.05..0.9999).contains(&target_power) {
        return Err(String::from(
            "Require ratio k > 0, alpha in (0, 0.5), and target_power in (0.05, 0.9999).",
        ));
    }

    let z_alpha = norm_ppf(1.0 - alpha_two_sided * 0.5);
    let z_beta = norm_ppf(target_power);
    let diff = (p_1 - p_0).abs();

    let p_bar = (p_1 + ratio_k * p_0) / (1.0 + ratio_k);
    let q_bar = 1.0 - p_bar;
    let term_null = z_alpha * ((1.0 + 1.0 / ratio_k) * p_bar * q_bar).sqrt();
    let term_alt = z_beta * (p_1 * (1.0 - p_1) + (p_0 * (1.0 - p_0)) / ratio_k).sqrt();

    let n1_prime = ((term_null + term_alt) / diff).powi(2);
    // Fleiss, Tytun & Ury (1980) Biometrics 36:343-346 Yates continuity correction:
    let cc_factor = 1.0 + (1.0 + 2.0 * (ratio_k + 1.0) / (n1_prime * ratio_k * diff)).sqrt();
    let n1_cc_f = (n1_prime / 4.0) * cc_factor.powi(2);

    let n1_unc = n1_prime.ceil() as usize;
    let n0_unc = (n1_prime * ratio_k).ceil() as usize;
    let n1_cc = n1_cc_f.ceil() as usize;
    let n0_cc = ((n1_cc as f64) * ratio_k).ceil() as usize;
    let total_cc = n1_cc + n0_cc;

    let swabs = ari_positivity.and_then(|pos| {
        if pos > 0.0 && pos < 1.0 {
            let need_for_cases = ((n1_cc as f64) / pos).ceil() as usize;
            let need_for_ctrls = ((n0_cc as f64) / (1.0 - pos)).ceil() as usize;
            Some(need_for_cases.max(need_for_ctrls))
        } else {
            None
        }
    });

    let mut curve = Vec::with_capacity(25);
    for step in 1..=25 {
        let scale = 0.16 * (step as f64) / 2.0; // 0.08x .. 2.0x of n1_cc
        let cand_n1 = ((n1_cc as f64) * scale).round().max(4.0) as usize;
        let cand_n0 = ((cand_n1 as f64) * ratio_k).round().max(4.0) as usize;
        let pwr = power_at_sample_size(p_0, p_1, cand_n0, cand_n1, alpha_two_sided, true);
        curve.push((cand_n0 + cand_n1, pwr));
    }

    Ok(StudyPowerPlan {
        kind,
        p_group_0: p_0,
        p_group_1: p_1,
        allocation_ratio_k: ratio_k,
        alpha_two_sided,
        target_power,
        effect_label,
        n1_uncorrected: n1_unc,
        n0_uncorrected: n0_unc,
        n1_fleiss_cc: n1_cc,
        n0_fleiss_cc: n0_cc,
        total_fleiss_cc: total_cc,
        total_ari_swabs_needed: swabs,
        power_curve: curve,
    })
}

/// Computes two-sided statistical power (`0.0..=1.0`) for independent group sizes `n_0` and `n_1`.
pub fn power_at_sample_size(
    p_0: f64,
    p_1: f64,
    n_0: usize,
    n_1: usize,
    alpha_two_sided: f64,
    continuity_corrected: bool,
) -> f64 {
    let f0 = n_0.max(1) as f64;
    let f1 = n_1.max(1) as f64;
    let diff = (p_1 - p_0).abs();
    let cc = if continuity_corrected {
        0.5 * (1.0 / f0 + 1.0 / f1)
    } else {
        0.0
    };
    let eff_diff = (diff - cc).max(0.0);
    if eff_diff <= 0.0 {
        return alpha_two_sided;
    }

    let z_alpha = norm_ppf(1.0 - alpha_two_sided * 0.5);
    let p_bar = (f0 * p_0 + f1 * p_1) / (f0 + f1);
    let se_null = (p_bar * (1.0 - p_bar) * (1.0 / f0 + 1.0 / f1)).sqrt();
    let se_alt = ((p_0 * (1.0 - p_0)) / f0 + (p_1 * (1.0 - p_1)) / f1).sqrt();
    if se_alt <= 1e-15 {
        return 1.0;
    }

    let z_power = (eff_diff - z_alpha * se_null) / se_alt;
    norm_cdf(z_power).clamp(0.0, 1.0)
}

/// Designs a Two-Arm Cohort / RCT study parameterized by baseline unexposed risk `p_unexposed`
/// and Relative Risk `relative_risk` (`p_exposed = p_unexposed * relative_risk`).
pub fn design_cohort_study(
    p_unexposed: f64,
    relative_risk: f64,
    unexposed_per_exposed_k: f64,
    alpha_two_sided: f64,
    target_power: f64,
) -> Result<StudyPowerPlan, String> {
    let p_1 = p_unexposed * relative_risk;
    if !(0.0001..=0.9999).contains(&p_1) {
        return Err(format!(
            "Exposed risk p1 = p0 * RR = {p_1:.4} must lie in (0, 1)."
        ));
    }
    solve_two_proportion_sample_size(
        EpiStudyKind::CohortRelativeRisk,
        p_unexposed,
        p_1,
        unexposed_per_exposed_k,
        alpha_two_sided,
        target_power,
        format!("RR = {relative_risk:.3}"),
        None,
    )
}

/// Designs an unmatched Case–Control study parameterized by control exposure prevalence `p_control_exp`
/// and target `odds_ratio`.
pub fn design_case_control_study(
    p_control_exposed: f64,
    odds_ratio: f64,
    controls_per_case_k: f64,
    alpha_two_sided: f64,
    target_power: f64,
) -> Result<StudyPowerPlan, String> {
    if odds_ratio <= 0.0 || !(0.0001..=0.9999).contains(&p_control_exposed) {
        return Err(String::from(
            "Require control exposure in (0, 1) and odds_ratio > 0.",
        ));
    }
    let odds0 = p_control_exposed / (1.0 - p_control_exposed);
    let odds1 = odds0 * odds_ratio;
    let p_case_exposed = odds1 / (1.0 + odds1);
    solve_two_proportion_sample_size(
        EpiStudyKind::CaseControlOddsRatio,
        p_control_exposed,
        p_case_exposed,
        controls_per_case_k,
        alpha_two_sided,
        target_power,
        format!("OR = {odds_ratio:.3}"),
        None,
    )
}

/// Designs a Jackson & Nelson (2013) Test-Negative Design (`TND`) vaccine effectiveness study
/// given population/test-negative vaccination coverage `vax_coverage_controls`, target
/// `vaccine_effectiveness` (`0..1`, where `OR = 1 - VE`), control-to-case ratio `k`, and optional
/// expected laboratory positivity among acute illness swabs (`expected_positivity`).
pub fn design_test_negative_study(
    vax_coverage_controls: f64,
    vaccine_effectiveness: f64,
    controls_per_case_k: f64,
    alpha_two_sided: f64,
    target_power: f64,
    expected_positivity: Option<f64>,
) -> Result<StudyPowerPlan, String> {
    if !(0.01..=0.99).contains(&vaccine_effectiveness) {
        return Err(String::from(
            "Target vaccine effectiveness (VE) must lie in 0.01..=0.99.",
        ));
    }
    let or = 1.0 - vaccine_effectiveness;
    let odds0 = vax_coverage_controls / (1.0 - vax_coverage_controls);
    let p_vax_cases = (odds0 * or) / (1.0 + odds0 * or);

    solve_two_proportion_sample_size(
        EpiStudyKind::TestNegativeVaccineEffectiveness,
        vax_coverage_controls,
        p_vax_cases,
        controls_per_case_k,
        alpha_two_sided,
        target_power,
        format!("VE = {:.1}% (OR = {or:.3})", vaccine_effectiveness * 100.0),
        expected_positivity,
    )
}

// ============================================================================
// 5. Publication Vector `.svg` Exporters for Rt, SEIR, and Power Curves
// ============================================================================

/// Renders a two-panel Cori et al. (2013) epidemic incidence + `R_t` (with 95% CrI ribbon and `R_t = 1` threshold) `.svg`.
pub fn render_rt_epidemic_svg(title: &str, pts: &[RtEstimatePoint]) -> String {
    let width = 820.0_f64;
    let height = 540.0_f64;
    let left = 68.0_f64;
    let pw = 716.0_f64;
    let n = pts.len().max(2);

    let max_inc = pts
        .iter()
        .map(|p| p.incidence)
        .fold(1.0_f64, f64::max)
        * 1.12;
    let max_rt = pts
        .iter()
        .filter_map(|p| p.ci95_high)
        .fold(2.5_f64, f64::max)
        .min(8.0)
        * 1.1;

    let mut svg = String::with_capacity(8192);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#ffffff\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Incidence &amp; Cori Renewal R_t</text>",
        xml_escape(title)
    );

    // Top panel: Daily Incidence bars (`y = 48..230`)
    let top_y = 48.0_f64;
    let top_h = 182.0_f64;
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top_y:.1}\" width=\"{pw:.1}\" height=\"{top_h:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let bar_w = (pw / (n as f64) * 0.72).max(2.0);
    for (i, p) in pts.iter().enumerate() {
        let x = left + ((i as f64 + 0.14) / (n as f64)) * pw;
        let bh = (p.incidence / max_inc) * top_h;
        let y = top_y + top_h - bh;
        let _ = writeln!(
            svg,
            "  <rect x=\"{x:.2}\" y=\"{y:.2}\" width=\"{bar_w:.2}\" height=\"{bh:.2}\" fill=\"#0284c7\"/>"
        );
    }

    // Bottom panel: R_t posterior mean + 95% CrI polygon + R_t = 1 dashed line (`y = 278..490`)
    let bot_y = 278.0_f64;
    let bot_h = 212.0_f64;
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{bot_y:.1}\" width=\"{pw:.1}\" height=\"{bot_h:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let y_rt1 = bot_y + bot_h - (1.0 / max_rt).clamp(0.0, 1.0) * bot_h;
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{y_rt1:.2}\" x2=\"{:.1}\" y2=\"{y_rt1:.2}\" stroke=\"#dc2626\" stroke-width=\"1.5\" stroke-dasharray=\"5,4\"/>",
        left + pw
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#dc2626\">R=1.0</text>",
        left - 6.0,
        y_rt1 + 4.0
    );

    let valid: Vec<(f64, f64, f64, f64)> = pts
        .iter()
        .filter_map(|p| {
            let m = p.rt_mean?;
            let lo = p.ci95_low?;
            let hi = p.ci95_high?.min(max_rt);
            let x = left + ((p.day_index as f64 + 0.5) / (n as f64)) * pw;
            Some((
                x,
                bot_y + bot_h - (m.min(max_rt) / max_rt) * bot_h,
                bot_y + bot_h - (lo.min(max_rt) / max_rt) * bot_h,
                bot_y + bot_h - (hi / max_rt) * bot_h,
            ))
        })
        .collect();

    if valid.len() >= 2 {
        let mut poly = String::new();
        for &(x, _, _, y_hi) in &valid {
            let _ = write!(poly, "{x:.2},{y_hi:.2} ");
        }
        for &(x, _, y_lo, _) in valid.iter().rev() {
            let _ = write!(poly, "{x:.2},{y_lo:.2} ");
        }
        let _ = writeln!(
            svg,
            "  <polygon points=\"{}\" fill=\"#bae6fd\" fill-opacity=\"0.65\"/>",
            poly.trim()
        );

        let mut line_pts = String::new();
        for &(x, ym, _, _) in &valid {
            let _ = write!(line_pts, "{x:.2},{ym:.2} ");
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#0369a1\" stroke-width=\"2.4\"/>",
            line_pts.trim()
        );
    }

    svg.push_str("</svg>\n");
    svg
}

/// Renders an `SIR`/`SEIR` compartmental trajectory `.svg` showing `S(t)`, `E(t)`, `I(t)`, and `R(t)`.
pub fn render_seir_trajectory_svg(sim: &SeirSimulationResult) -> String {
    let width = 820.0_f64;
    let height = 460.0_f64;
    let left = 72.0_f64;
    let top = 54.0_f64;
    let pw = 600.0_f64;
    let ph = 350.0_f64;
    let pop = sim.params.population.max(1.0);
    let max_d = (sim.params.total_days.max(1)) as f64;

    let mut svg = String::with_capacity(8192);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#ffffff\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">SEIR Compartmental Sweep (R0 = {:.2} | Peak Day {:.0} | Final Attack {:.1}%)</text>",
        sim.params.r0,
        sim.peak_day,
        sim.final_attack_fraction * 100.0
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    let series: [(&str, &str, fn(&crate::epi_models::SeirTimePoint) -> f64); 4] = [
        ("Susceptible (S)", "#2563eb", |p| p.susceptible),
        ("Exposed (E)", "#d97706", |p| p.exposed),
        ("Infectious (I)", "#dc2626", |p| p.infectious),
        ("Recovered/Immune (R)", "#16a34a", |p| p.recovered),
    ];

    for (idx, (lbl, col, getter)) in series.iter().enumerate() {
        let mut pts = String::new();
        for p in &sim.trajectory {
            let x = left + (p.day / max_d) * pw;
            let y = top + ph - (getter(p) / pop).clamp(0.0, 1.0) * ph;
            let _ = write!(pts, "{x:.2},{y:.2} ");
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.2\"/>",
            pts.trim()
        );
        let ly = top + 24.0 + (idx as f64) * 24.0;
        let lx = left + pw + 16.0;
        let _ = writeln!(
            svg,
            "  <line x1=\"{lx:.1}\" y1=\"{ly:.1}\" x2=\"{:.1}\" y2=\"{ly:.1}\" stroke=\"{col}\" stroke-width=\"3\"/>",
            lx + 18.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"11.5\" fill=\"#1e293b\">{lbl}</text>",
            lx + 24.0,
            ly + 4.0
        );
    }

    svg.push_str("</svg>\n");
    svg
}

/// Renders a statistical power-by-sample-size curve `.svg` highlighting the Fleiss target `N` and `1 - beta` threshold.
pub fn render_power_curve_svg(plan: &StudyPowerPlan) -> String {
    let width = 780.0_f64;
    let height = 440.0_f64;
    let left = 68.0_f64;
    let top = 56.0_f64;
    let pw = 668.0_f64;
    let ph = 330.0_f64;

    let max_n = plan
        .power_curve
        .last()
        .map(|p| p.0 as f64)
        .unwrap_or(plan.total_fleiss_cc as f64)
        .max(10.0);

    let mut svg = begin_svg_canvas(4096, width, height, Some("#ffffff"));
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"14.5\" font-weight=\"700\" fill=\"#0f172a\">{} — {} (Required N = {})</text>",
        xml_escape(plan.kind.title()),
        xml_escape(&plan.effect_label),
        plan.total_fleiss_cc
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    let y_target = top + ph - plan.target_power * ph;
    let x_target = left + ((plan.total_fleiss_cc as f64) / max_n).clamp(0.0, 1.0) * pw;
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{y_target:.1}\" x2=\"{:.1}\" y2=\"{y_target:.1}\" stroke=\"#16a34a\" stroke-dasharray=\"4,4\"/>",
        left + pw
    );
    let _ = writeln!(
        svg,
        "  <line x1=\"{x_target:.1}\" y1=\"{top:.1}\" x2=\"{x_target:.1}\" y2=\"{:.1}\" stroke=\"#dc2626\" stroke-dasharray=\"4,4\"/>",
        top + ph
    );

    let mut poly = String::new();
    for &(n, pwr) in &plan.power_curve {
        let x = left + ((n as f64) / max_n) * pw;
        let y = top + ph - pwr.clamp(0.0, 1.0) * ph;
        let _ = write!(poly, "{x:.2},{y:.2} ");
    }
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#4f46e5\" stroke-width=\"2.6\"/>",
        poly.trim()
    );
    let _ = writeln!(
        svg,
        "  <circle cx=\"{x_target:.2}\" cy=\"{y_target:.2}\" r=\"5\" fill=\"#dc2626\"/>"
    );

    svg.push_str("</svg>\n");
    svg
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn computes_fleiss_cohort_and_jackson_nelson_tnd_sample_sizes() {
        // Classic Fleiss benchmark: p0 = 0.20, RR = 0.50 (p1 = 0.10), k = 1, alpha = 0.05, power = 0.80
        // Uncorrected n per arm ~ 199, Fleiss CC ~ 219 per arm (438 total)
        let coh = design_cohort_study(0.20, 0.50, 1.0, 0.05, 0.80).unwrap();
        assert!((coh.n1_uncorrected as isize - 199).abs() <= 2);
        assert!((coh.n1_fleiss_cc as isize - 219).abs() <= 2);

        // Test-Negative Design: 50% vax coverage in controls, VE = 60% (OR = 0.40 -> p_case = 0.2857), k = 2 controls/case
        let tnd = design_test_negative_study(0.50, 0.60, 2.0, 0.05, 0.80, Some(0.25)).unwrap();
        assert!(tnd.n1_fleiss_cc > 50 && tnd.n1_fleiss_cc < 130);
        assert!(tnd.total_ari_swabs_needed.unwrap() > tnd.total_fleiss_cc);
    }
}
