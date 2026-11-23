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

//! Tool #20 (`osabt flux`) — Odum (1956) Diel Dissolved-Oxygen Stream Metabolism Engine
//! & Two-Panel `.svg` Ecophysiology / Watershed Plotter:
//! - **Primary Production in Flowing Waters** — Odum (1956) *Limnol. Oceanogr.* 1(2):102–117
//!   (`https://www.emergysociety.com/wp-content/uploads/OdumHT.1956.Primary-production-in-flowing-waters.-Limnol-and-Oceanogr.V.1.pdf`)
//! - Estimates daily `GPP`, `ER`, `NEP = GPP - ER` (`g O2 m^-2 d^-1`), `P/R` ratio, and `K600` (`d^-1`).

use std::fmt::Write as _;

use crate::flux_leaf::{eval_fvcb_regime_rates, FvcbAciFit, FvcbLimitationState, LeafAciPoint};
use crate::svg_util::{begin_svg_canvas, format_polyline_points, xml_escape};
use crate::wildlife_occupancy::nelder_mead_minimize;

/// Single hourly/sub-daily stream dissolved-oxygen (`DO`) and irradiance observation.
#[derive(Debug, Clone)]
pub struct StreamDoPoint {
    pub hour: f64,
    pub do_mg_l: f64,
    pub temp_c: f64,
    pub par_umol_m2_s: f64,
    pub depth_m: f64,
}

/// Fitted Odum (1956) open-channel daily stream metabolism summary.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct StreamMetabolismReport {
    /// Gross Primary Production (`g O2 m^-2 d^-1`, `>= 0`).
    pub gpp_g_o2_m2_d: f64,
    /// Ecosystem Respiration (`g O2 m^-2 d^-1`, positive magnitude of oxygen consumption).
    pub er_g_o2_m2_d: f64,
    /// Net Ecosystem Production `NEP = GPP - ER` (`> 0` autotrophic, `< 0` heterotrophic).
    pub nep_g_o2_m2_d: f64,
    pub p_to_r_ratio: f64,
    /// Standardized gas-exchange reaeration coefficient `K600` (`d^-1`).
    pub k600_per_day: f64,
    pub trophic_state: String,
    pub do_rmse_mg_l: f64,
    pub modeled_do_mg_l: Vec<f64>,
    pub do_sat_mg_l: Vec<f64>,
}

/// Temperature dependence model for daily Ecosystem Respiration (`ER`) in Odum stream metabolism.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum ErTempScaling {
    /// Constant temperature-independent `ER` (`Q10 = 1.0`, classic Odum 1956).
    Constant,
    /// Exponential `Q10` scaling relative to `ref_temp_c`: `ER(T) = ER_ref * Q10^((T - T_ref) / 10)`.
    Q10 { q10: f64, ref_temp_c: f64 },
    /// Arrhenius activation-energy scaling relative to `ref_temp_c` (`ea_kj_mol` in `kJ mol^-1`).
    Arrhenius { ea_kj_mol: f64, ref_temp_c: f64 },
}

/// Optional configuration for [`fit_odum_stream_metabolism_with_config`] supporting
/// barometric/elevation `O2` saturation adjustment and `ER` temperature scaling.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct StreamMetabolismConfig {
    /// Ambient barometric pressure in atmospheres (`1.0 atm = 101.325 kPa` at sea level).
    pub pressure_atm: f64,
    /// Ecosystem respiration temperature scaling (`Constant`, `Q10`, or `Arrhenius`).
    pub er_temp_scaling: ErTempScaling,
}

impl Default for StreamMetabolismConfig {
    fn default() -> Self {
        Self {
            pressure_atm: 1.0,
            er_temp_scaling: ErTempScaling::Constant,
        }
    }
}

impl StreamMetabolismConfig {
    /// Constructs a configuration using standard U.S. Atmosphere barometric pressure at `elevation_m`.
    pub fn at_elevation_m(elevation_m: f64, er_temp_scaling: ErTempScaling) -> Self {
        Self {
            pressure_atm: barometric_pressure_atm_at_elevation(elevation_m),
            er_temp_scaling,
        }
    }
}

/// Computes standard atmospheric barometric pressure (`atm`, relative to sea level `P0 = 1.0 atm`)
/// as a function of stream elevation (`m` above sea level).
#[inline]
pub fn barometric_pressure_atm_at_elevation(elevation_m: f64) -> f64 {
    let h = elevation_m.clamp(-500.0, 11_000.0);
    (1.0 - 2.255_77e-5 * h).powf(5.255_88)
}

/// Computes Benson & Krause (1984) freshwater dissolved-oxygen saturation (`mg L^-1` at `1 atm`).
#[inline]
pub fn calc_do_saturation_mg_l(temp_c: f64) -> f64 {
    let tk = temp_c.clamp(0.0, 40.0) + 273.15;
    let ln_c = -139.344_11 + 1.575_701e5 / tk - 6.642_308e7 / (tk * tk)
        + 1.243_800e10 / (tk * tk * tk)
        - 8.621_949e11 / (tk * tk * tk * tk);
    ln_c.exp()
}

/// Computes Benson & Krause (1984) freshwater dissolved-oxygen saturation (`mg L^-1`) corrected
/// for barometric pressure (`P / P0`, where `P0 = 1.0 atm`).
#[inline]
pub fn calc_do_saturation_mg_l_at_pressure(temp_c: f64, pressure_atm: f64) -> f64 {
    calc_do_saturation_mg_l(temp_c) * pressure_atm.max(0.0)
}

/// Computes Benson & Krause (1984) freshwater dissolved-oxygen saturation (`mg L^-1`) corrected
/// for elevation (`m` above sea level) via [`barometric_pressure_atm_at_elevation`].
#[inline]
pub fn calc_do_saturation_mg_l_at_elevation(temp_c: f64, elevation_m: f64) -> f64 {
    calc_do_saturation_mg_l_at_pressure(temp_c, barometric_pressure_atm_at_elevation(elevation_m))
}

/// Scales reference Ecosystem Respiration `er_ref` to `temp_c` using a `Q10` temperature coefficient.
#[inline]
pub fn scale_er_q10(er_ref: f64, temp_c: f64, ref_temp_c: f64, q10: f64) -> f64 {
    if !er_ref.is_finite() || !temp_c.is_finite() || !ref_temp_c.is_finite() || !q10.is_finite() {
        return er_ref;
    }
    er_ref * q10.max(0.1).powf((temp_c - ref_temp_c) / 10.0)
}

/// Scales reference Ecosystem Respiration `er_ref` to `temp_c` using an Arrhenius activation energy `ea_kj_mol`.
#[inline]
pub fn scale_er_arrhenius(er_ref: f64, temp_c: f64, ref_temp_c: f64, ea_kj_mol: f64) -> f64 {
    if !er_ref.is_finite() || !temp_c.is_finite() || !ref_temp_c.is_finite() || !ea_kj_mol.is_finite() {
        return er_ref;
    }
    let r_gas = 0.008_314_5_f64;
    let tk = (temp_c + 273.15).max(1.0);
    let t_ref_k = (ref_temp_c + 273.15).max(1.0);
    er_ref * ((ea_kj_mol / r_gas) * (1.0 / t_ref_k - 1.0 / tk)).exp()
}

/// Evaluates instantaneous Ecosystem Respiration `ER(temp_c)` under the specified [`ErTempScaling`] law.
#[inline]
pub fn scale_er_temp(scaling: ErTempScaling, er_ref: f64, temp_c: f64) -> f64 {
    match scaling {
        ErTempScaling::Constant => er_ref,
        ErTempScaling::Q10 { q10, ref_temp_c } => scale_er_q10(er_ref, temp_c, ref_temp_c, q10),
        ErTempScaling::Arrhenius {
            ea_kj_mol,
            ref_temp_c,
        } => scale_er_arrhenius(er_ref, temp_c, ref_temp_c, ea_kj_mol),
    }
}

/// Wanninkhof (1992) freshwater `O2` Schmidt number `Sc(T)` for scaling `K600` -> `K(T)`.
#[inline]
fn schmidt_o2(temp_c: f64) -> f64 {
    let t = temp_c.clamp(0.0, 40.0);
    1800.6 - 120.10 * t + 3.7818 * t * t - 0.047_608 * t * t * t
}

/// Parses a diel stream dissolved-oxygen CSV/TSV (`hour,do_mg_l,temp_c,par[,depth_m]`).
pub fn parse_stream_do_csv(raw_csv: &str) -> Result<Vec<StreamDoPoint>, String> {
    let mut pts = Vec::new();
    let mut hdr = false;
    for (idx, line) in raw_csv.trim_start_matches('\u{feff}').lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();
        if !hdr {
            hdr = true;
            let c0 = cols[0].to_ascii_lowercase();
            if c0.contains("hour") || c0.contains("time") {
                continue;
            }
        }
        if cols.len() < 4 {
            return Err(format!(
                "Stream DO CSV line {}: expected >= 4 columns (hour,do_mg_l,temp_c,par[,depth_m]).",
                idx + 1
            ));
        }
        let hour = cols[0]
            .parse::<f64>()
            .map_err(|_| format!("Stream line {}: invalid hour '{}'.", idx + 1, cols[0]))?;
        let do_val = cols[1]
            .parse::<f64>()
            .map_err(|_| format!("Stream line {}: invalid DO '{}'.", idx + 1, cols[1]))?;
        let temp = cols[2]
            .parse::<f64>()
            .map_err(|_| format!("Stream line {}: invalid temp '{}'.", idx + 1, cols[2]))?;
        let par = cols[3]
            .parse::<f64>()
            .map_err(|_| format!("Stream line {}: invalid PAR '{}'.", idx + 1, cols[3]))?;
        let depth = cols.get(4).and_then(|s| s.parse::<f64>().ok()).unwrap_or(0.35);

        if !hour.is_finite()
            || !do_val.is_finite()
            || !temp.is_finite()
            || !par.is_finite()
            || !depth.is_finite()
        {
            return Err(format!(
                "Stream line {}: non-finite numeric value encountered.",
                idx + 1
            ));
        }

        pts.push(StreamDoPoint {
            hour,
            do_mg_l: do_val,
            temp_c: temp,
            par_umol_m2_s: par.max(0.0),
            depth_m: depth.max(0.05),
        });
    }
    if pts.len() < 8 {
        return Err(String::from(
            "Stream metabolism fitting requires at least 8 diel dissolved-oxygen measurements.",
        ));
    }
    pts.sort_by(|a, b| a.hour.total_cmp(&b.hour));
    Ok(pts)
}

/// Fits Odum's (1956) open-channel oxygen balance (`GPP`, `ER`, `K600`) to diel `DO` observations
/// at sea-level pressure (`1 atm`) and constant `ER` via Crank–Nicolson/midpoint ODE integration
/// and Nelder–Mead SSE minimization.
pub fn fit_odum_stream_metabolism(points: &[StreamDoPoint]) -> Result<StreamMetabolismReport, String> {
    fit_odum_stream_metabolism_with_config(points, &StreamMetabolismConfig::default())
}

/// Fits Odum's (1956) open-channel oxygen balance (`GPP`, `ER`, `K600`) with configurable
/// barometric pressure (`pressure_atm`) and Ecosystem Respiration temperature scaling (`Q10` / `Arrhenius`).
pub fn fit_odum_stream_metabolism_with_config(
    points: &[StreamDoPoint],
    cfg: &StreamMetabolismConfig,
) -> Result<StreamMetabolismReport, String> {
    if points.len() < 8 {
        return Err(String::from("At least 8 diel points are required."));
    }
    if points.iter().any(|p| {
        !p.hour.is_finite()
            || !p.do_mg_l.is_finite()
            || !p.temp_c.is_finite()
            || !p.par_umol_m2_s.is_finite()
            || !p.depth_m.is_finite()
    }) {
        return Err(String::from(
            "Stream metabolism points contain non-finite (NaN/Inf) values.",
        ));
    }
    if !cfg.pressure_atm.is_finite() || cfg.pressure_atm <= 0.0 {
        return Err(String::from("Stream pressure_atm must be finite and > 0."));
    }

    let mut sum_par_dt = 0.0_f64;
    let mut total_days = 0.0_f64;
    for i in 1..points.len() {
        let dt_days = ((points[i].hour - points[i - 1].hour) / 24.0).max(1e-4);
        let par_mid = 0.5 * (points[i].par_umol_m2_s + points[i - 1].par_umol_m2_s);
        sum_par_dt += par_mid * dt_days;
        total_days += dt_days;
    }
    let par_norm = sum_par_dt.max(1e-6);
    let do_sat: Vec<f64> = points
        .iter()
        .map(|p| calc_do_saturation_mg_l_at_pressure(p.temp_c, cfg.pressure_atm))
        .collect();

    let simulate_do = |gpp: f64, er_ref: f64, k600: f64| -> (Vec<f64>, f64) {
        let mut pred = Vec::with_capacity(points.len());
        pred.push(points[0].do_mg_l);
        let mut er_integral = 0.0_f64;
        for i in 1..points.len() {
            let dt_days = ((points[i].hour - points[i - 1].hour) / 24.0).max(1e-4);
            let t_mid = 0.5 * (points[i].temp_c + points[i - 1].temp_c);
            let z_mid = 0.5 * (points[i].depth_m + points[i - 1].depth_m);
            let par_mid = 0.5 * (points[i].par_umol_m2_s + points[i - 1].par_umol_m2_s);
            let sat_mid = 0.5 * (do_sat[i] + do_sat[i - 1]);

            let k_t = k600 * (schmidt_o2(t_mid) / 600.0).powf(-0.5);
            let er_t = scale_er_temp(cfg.er_temp_scaling, er_ref, t_mid).max(0.0);
            er_integral += er_t * dt_days;
            let bio_flux = (gpp * (par_mid / par_norm) - er_t) / z_mid;
            // Semi-implicit step for unconditional stability across large K600
            let prev = pred[i - 1];
            let next = (prev + dt_days * (bio_flux + k_t * sat_mid)) / (1.0 + dt_days * k_t);
            pred.push(next);
        }
        let er_daily = er_integral / total_days.max(1e-6);
        (pred, er_daily)
    };

    let sse = |th: &[f64]| -> f64 {
        let (gpp, er_ref, k600) = (th[0], th[1], th[2]);
        if gpp < 0.0 || gpp > 40.0 || er_ref < 0.05 || er_ref > 40.0 || k600 < 0.2 || k600 > 60.0 {
            return 1e12 + gpp * gpp + er_ref * er_ref + k600 * k600;
        }
        let (pred, _) = simulate_do(gpp, er_ref, k600);
        pred.iter()
            .zip(points)
            .map(|(&pr, ob)| (pr - ob.do_mg_l).powi(2))
            .sum()
    };

    let starts = [[3.0, 4.5, 6.0], [1.5, 2.5, 10.0], [6.0, 5.0, 4.0]];
    let mut best_th = starts[0].to_vec();
    let mut best_sse = f64::INFINITY;
    for st in &starts {
        let (th, val) = nelder_mead_minimize(&sse, st, 900, 1e-9);
        if val < best_sse {
            best_sse = val;
            best_th = th;
        }
    }
    let (best_th, best_sse) = nelder_mead_minimize(&sse, &best_th, 1100, 1e-11);
    let (gpp, er_ref, k600) = (best_th[0], best_th[1], best_th[2]);
    let (modeled, er) = simulate_do(gpp, er_ref, k600);
    let nep = gpp - er;
    let pr = gpp / er.max(1e-6);
    let trophic = if nep > 0.15 {
        String::from("Net Autotrophic (P/R > 1)")
    } else if nep < -0.15 {
        String::from("Net Heterotrophic (P/R < 1)")
    } else {
        String::from("Metabolically Balanced (P/R ~ 1)")
    };

    Ok(StreamMetabolismReport {
        gpp_g_o2_m2_d: gpp,
        er_g_o2_m2_d: er,
        nep_g_o2_m2_d: nep,
        p_to_r_ratio: pr,
        k600_per_day: k600,
        trophic_state: trophic,
        do_rmse_mg_l: (best_sse / (points.len() as f64)).sqrt(),
        modeled_do_mg_l: modeled,
        do_sat_mg_l: do_sat,
    })
}

/// Computes auto-scaled `(min_do, max_do)` Y-axis bounds (`mg L^-1`) across observed, modeled,
/// and saturation dissolved-oxygen series with padding and a degenerate constant-data guard.
pub fn compute_do_y_bounds(
    do_pts: &[StreamDoPoint],
    stream_rep: &StreamMetabolismReport,
) -> (f64, f64) {
    let mut raw_min = f64::INFINITY;
    let mut raw_max = f64::NEG_INFINITY;
    for v in do_pts
        .iter()
        .map(|p| p.do_mg_l)
        .chain(stream_rep.modeled_do_mg_l.iter().copied())
        .chain(stream_rep.do_sat_mg_l.iter().copied())
    {
        if v.is_finite() {
            raw_min = raw_min.min(v);
            raw_max = raw_max.max(v);
        }
    }
    if !raw_min.is_finite() || !raw_max.is_finite() {
        return (0.0, 14.0);
    }
    let span = raw_max - raw_min;
    if span < 1e-6 {
        let pad = (raw_min.abs() * 0.15).max(0.5);
        let lo = if raw_min >= 0.0 {
            (raw_min - pad).max(0.0)
        } else {
            raw_min - pad
        };
        let hi = raw_max + pad;
        return (lo, if hi > lo + 1e-3 { hi } else { lo + 1.0 });
    }
    let pad = (span * 0.10).max(0.25);
    let lo = if raw_min >= 0.0 {
        (raw_min - pad).max(0.0)
    } else {
        raw_min - pad
    };
    let hi = raw_max + pad;
    (lo, if hi > lo + 1e-3 { hi } else { lo + 1.0 })
}

/// Renders a two-panel `.svg` figure pairing the leaf `FvCB` $A\text{–}C_i$ curve (`left`,
/// including `Wc`, `Wj`, and `Wp` limitation curves) with the diel stream dissolved-oxygen
/// metabolism curve (`right`, with auto-scaled DO bounds).
pub fn render_flux_svg(
    aci_pts: &[LeafAciPoint],
    aci_fit: &FvcbAciFit,
    do_pts: &[StreamDoPoint],
    stream_rep: &StreamMetabolismReport,
    title: &str,
) -> String {
    let width = 960.0_f64;
    let height = 472.0_f64;
    let p1_x = 72.0_f64;
    let p1_w = 390.0_f64;
    let p2_x = 548.0_f64;
    let p2_w = 368.0_f64;
    let top = 74.0_f64;
    let ph = 296.0_f64;

    let mut svg = begin_svg_canvas(24_576, width, height, Some("#ffffff"));
    let _ = writeln!(
        svg,
        "  <text x=\"{p1_x:.0}\" y=\"26\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — FvCB Leaf A-Ci &amp; Odum Stream DO Metabolism</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p1_x:.1}\" y=\"{top:.1}\" width=\"{p1_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#94a3b8\"/>"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p2_x:.1}\" y=\"{top:.1}\" width=\"{p2_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#94a3b8\"/>"
    );

    // Panel 1: FvCB A-Ci curve + Wc / Wj / Wp limitation curves
    let raw_max_ci = aci_pts
        .iter()
        .map(|p| p.ci_ubar)
        .filter(|v| v.is_finite())
        .fold(0.0_f64, f64::max);
    let max_ci = (if raw_max_ci > 0.0 { raw_max_ci } else { 1600.0 }) * 1.05;
    let mut raw_min_a = f64::INFINITY;
    let mut raw_max_a = f64::NEG_INFINITY;
    for v in aci_pts
        .iter()
        .map(|p| p.a_net)
        .chain(aci_fit.predicted_a.iter().copied())
    {
        if v.is_finite() {
            raw_min_a = raw_min_a.min(v);
            raw_max_a = raw_max_a.max(v);
        }
    }
    let min_a = if raw_min_a.is_finite() {
        raw_min_a.min(-2.0) - 0.5
    } else {
        -2.5
    };
    let max_a = if raw_max_a.is_finite() {
        (raw_max_a * 1.15).max(min_a + 5.0)
    } else {
        25.0
    };
    let map_x1 = |ci: f64| p1_x + (ci / max_ci).clamp(0.0, 1.0) * p1_w;
    let map_y1 = |a: f64| top + ph - ((a - min_a) / (max_a - min_a).max(1.0)).clamp(0.0, 1.0) * ph;

    // Panel 1 gridlines, ticks, and numeric axis labels
    for k in 0..=4 {
        let frac = (k as f64) / 4.0;
        let ci_tick = frac * max_ci;
        let tx = map_x1(ci_tick);
        if k > 0 && k < 4 {
            let _ = writeln!(
                svg,
                "  <line x1=\"{tx:.1}\" y1=\"{top:.1}\" x2=\"{tx:.1}\" y2=\"{:.1}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\"/>",
                top + ph
            );
        }
        let _ = writeln!(
            svg,
            "  <line x1=\"{tx:.1}\" y1=\"{:.1}\" x2=\"{tx:.1}\" y2=\"{:.1}\" stroke=\"#64748b\" stroke-width=\"1.0\"/>",
            top + ph,
            top + ph + 5.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{tx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#475569\">{ci_tick:.0}</text>",
            top + ph + 18.0
        );

        let a_tick = min_a + frac * (max_a - min_a);
        let ty = map_y1(a_tick);
        if k > 0 && k < 4 {
            let _ = writeln!(
                svg,
                "  <line x1=\"{p1_x:.1}\" y1=\"{ty:.1}\" x2=\"{:.1}\" y2=\"{ty:.1}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\"/>",
                p1_x + p1_w
            );
        }
        let _ = writeln!(
            svg,
            "  <line x1=\"{:.1}\" y1=\"{ty:.1}\" x2=\"{p1_x:.1}\" y2=\"{ty:.1}\" stroke=\"#64748b\" stroke-width=\"1.0\"/>",
            p1_x - 5.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#475569\">{a_tick:.1}</text>",
            p1_x - 8.0,
            ty + 3.5
        );
    }

    let ppfd_ref = if aci_pts.is_empty() {
        1500.0
    } else {
        aci_pts.iter().map(|p| p.ppfd).sum::<f64>() / (aci_pts.len() as f64)
    };
    let t_leaf_ref = if aci_pts.is_empty() {
        25.0
    } else {
        aci_pts.iter().map(|p| p.t_leaf_c).sum::<f64>() / (aci_pts.len() as f64)
    };

    let mut wc_pts: Vec<(f64, f64)> = Vec::with_capacity(81);
    let mut wj_pts: Vec<(f64, f64)> = Vec::with_capacity(81);
    let mut wp_pts: Vec<(f64, f64)> = Vec::with_capacity(81);
    let mut net_pts: Vec<(f64, f64)> = Vec::with_capacity(81);

    let steps = 80usize;
    let mut prev_wc: Option<(f64, f64)> = None;
    for k in 0..=steps {
        let ci = 35.0 + ((k as f64) / (steps as f64)) * (max_ci * 0.97 - 35.0).max(50.0);
        let (a_wc, a_wj, a_wp, a_net, _) = eval_fvcb_regime_rates(
            ci,
            ppfd_ref,
            t_leaf_ref,
            aci_fit.vcmax,
            aci_fit.jmax,
            aci_fit.tpu,
            aci_fit.rd,
        );

        if a_wc.is_finite() {
            if a_wc >= min_a && a_wc <= max_a {
                wc_pts.push((map_x1(ci), map_y1(a_wc)));
            } else if a_wc > max_a {
                if let Some((p_ci, p_a)) = prev_wc {
                    if p_a <= max_a && (a_wc - p_a).abs() > 1e-9 {
                        let frac = ((max_a - p_a) / (a_wc - p_a)).clamp(0.0, 1.0);
                        let cross_ci = p_ci + frac * (ci - p_ci);
                        wc_pts.push((map_x1(cross_ci), map_y1(max_a)));
                    }
                }
            }
            prev_wc = Some((ci, a_wc));
        }

        if a_wj.is_finite() && a_wj >= min_a && a_wj <= max_a {
            wj_pts.push((map_x1(ci), map_y1(a_wj)));
        }
        if a_wp.is_finite() && a_wp >= min_a && a_wp <= max_a {
            wp_pts.push((map_x1(ci), map_y1(a_wp)));
        }
        if a_net.is_finite() {
            net_pts.push((map_x1(ci), map_y1(a_net)));
        }
    }

    if wc_pts.len() >= 2 {
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#2563eb\" stroke-dasharray=\"5,3\" stroke-width=\"1.5\"/>",
            format_polyline_points(&wc_pts)
        );
    }
    if wj_pts.len() >= 2 {
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#d97706\" stroke-dasharray=\"4,3\" stroke-width=\"1.5\"/>",
            format_polyline_points(&wj_pts)
        );
    }
    if wp_pts.len() >= 2 {
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#dc2626\" stroke-dasharray=\"2,3\" stroke-width=\"1.5\"/>",
            format_polyline_points(&wp_pts)
        );
    }
    if net_pts.len() >= 2 {
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#16a34a\" stroke-width=\"2.4\"/>",
            format_polyline_points(&net_pts)
        );
    }

    for (idx, p) in aci_pts.iter().enumerate() {
        let pt_fill = match aci_fit.limiting_states.get(idx) {
            Some(FvcbLimitationState::RubiscoWc) => "#2563eb",
            Some(FvcbLimitationState::RuBpRegenWj) => "#d97706",
            Some(FvcbLimitationState::TriosePhosphateWp) => "#dc2626",
            None => "#0f172a",
        };
        let _ = writeln!(
            svg,
            "  <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"4.0\" fill=\"{pt_fill}\" stroke=\"#ffffff\" stroke-width=\"1.0\"/>",
            map_x1(p.ci_ubar),
            map_y1(p.a_net)
        );
    }

    // Panel 1 legend & Y-axis label
    let leg1_y = 58.0_f64;
    let _ = writeln!(
        svg,
        "  <text x=\"{p1_x:.1}\" y=\"44\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"#0f172a\">Leaf FvCB A-Ci Curve</text>"
    );
    for (dx, col, dash, lbl) in [
        (0.0, "#16a34a", "", "A_net"),
        (70.0, "#2563eb", " stroke-dasharray=\"5,3\"", "Wc"),
        (124.0, "#d97706", " stroke-dasharray=\"4,3\"", "Wj"),
        (178.0, "#dc2626", " stroke-dasharray=\"2,3\"", "Wp"),
    ] {
        let lx = p1_x + dx;
        let _ = writeln!(
            svg,
            "  <line x1=\"{lx:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{col}\" stroke-width=\"2.0\"{dash}/>",
            leg1_y - 3.5,
            lx + 18.0,
            leg1_y - 3.5
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{leg1_y:.1}\" font-family=\"sans-serif\" font-size=\"10.5\" fill=\"#334155\">{lbl}</text>",
            lx + 22.0
        );
    }
    let _ = writeln!(
        svg,
        "  <text x=\"18\" y=\"{:.1}\" text-anchor=\"middle\" transform=\"rotate(-90 18 {:.1})\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#334155\">A_n (umol m^-2 s^-1)</text>",
        top + 0.5 * ph,
        top + 0.5 * ph
    );

    // Panel 2: Diel DO observed vs modeled + DO saturation (auto-scaled)
    let min_hr = do_pts
        .iter()
        .map(|p| p.hour)
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let max_hr = do_pts
        .iter()
        .map(|p| p.hour)
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    let (min_hr, max_hr) = if min_hr.is_finite() && max_hr.is_finite() {
        if (max_hr - min_hr) < 1e-6 {
            (min_hr.min(0.0), (max_hr + 1.0).max(min_hr + 1.0))
        } else {
            (min_hr.min(0.0), max_hr)
        }
    } else {
        (0.0, 24.0)
    };
    let (min_do, max_do) = compute_do_y_bounds(do_pts, stream_rep);
    let map_x2 = |h: f64| {
        p2_x + ((h - min_hr) / (max_hr - min_hr).max(1.0)).clamp(0.0, 1.0) * p2_w
    };
    let map_y2 = |d: f64| {
        top + ph - ((d - min_do) / (max_do - min_do).max(1e-6)).clamp(0.0, 1.0) * ph
    };

    // Panel 2 gridlines, ticks, and numeric axis labels
    for k in 0..=4 {
        let frac = (k as f64) / 4.0;
        let hr_tick = min_hr + frac * (max_hr - min_hr);
        let tx = map_x2(hr_tick);
        if k > 0 && k < 4 {
            let _ = writeln!(
                svg,
                "  <line x1=\"{tx:.1}\" y1=\"{top:.1}\" x2=\"{tx:.1}\" y2=\"{:.1}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\"/>",
                top + ph
            );
        }
        let _ = writeln!(
            svg,
            "  <line x1=\"{tx:.1}\" y1=\"{:.1}\" x2=\"{tx:.1}\" y2=\"{:.1}\" stroke=\"#64748b\" stroke-width=\"1.0\"/>",
            top + ph,
            top + ph + 5.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{tx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#475569\">{hr_tick:.0}h</text>",
            top + ph + 18.0
        );

        let do_tick = min_do + frac * (max_do - min_do);
        let ty = map_y2(do_tick);
        if k > 0 && k < 4 {
            let _ = writeln!(
                svg,
                "  <line x1=\"{p2_x:.1}\" y1=\"{ty:.1}\" x2=\"{:.1}\" y2=\"{ty:.1}\" stroke=\"#e2e8f0\" stroke-width=\"0.8\"/>",
                p2_x + p2_w
            );
        }
        let _ = writeln!(
            svg,
            "  <line x1=\"{:.1}\" y1=\"{ty:.1}\" x2=\"{p2_x:.1}\" y2=\"{ty:.1}\" stroke=\"#64748b\" stroke-width=\"1.0\"/>",
            p2_x - 5.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"10.5\" fill=\"#475569\">{do_tick:.1}</text>",
            p2_x - 8.0,
            ty + 3.5
        );
    }

    let mut mod_pts: Vec<(f64, f64)> = Vec::with_capacity(do_pts.len());
    let mut sat_pts: Vec<(f64, f64)> = Vec::with_capacity(do_pts.len());
    for (i, p) in do_pts.iter().enumerate() {
        let xp = map_x2(p.hour);
        if let Some(&m_do) = stream_rep.modeled_do_mg_l.get(i) {
            if m_do.is_finite() {
                mod_pts.push((xp, map_y2(m_do)));
            }
        }
        if let Some(&s_do) = stream_rep.do_sat_mg_l.get(i) {
            if s_do.is_finite() {
                sat_pts.push((xp, map_y2(s_do)));
            }
        }
        if p.do_mg_l.is_finite() {
            let _ = writeln!(
                svg,
                "  <circle cx=\"{xp:.2}\" cy=\"{:.2}\" r=\"3.5\" fill=\"#dc2626\"/>",
                map_y2(p.do_mg_l)
            );
        }
    }
    if sat_pts.len() >= 2 {
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#64748b\" stroke-dasharray=\"4,4\" stroke-width=\"1.5\"/>",
            format_polyline_points(&sat_pts)
        );
    }
    if mod_pts.len() >= 2 {
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.4\"/>",
            format_polyline_points(&mod_pts)
        );
    }

    // Panel 2 legend & Y-axis label
    let _ = writeln!(
        svg,
        "  <text x=\"{p2_x:.1}\" y=\"44\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"#0f172a\">Odum Diel Stream DO Metabolism</text>"
    );
    let _ = writeln!(
        svg,
        "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"3.5\" fill=\"#dc2626\"/>",
        p2_x + 5.0,
        leg1_y - 3.5
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{leg1_y:.1}\" font-family=\"sans-serif\" font-size=\"10.5\" fill=\"#334155\">Obs DO</text>",
        p2_x + 12.0
    );
    for (dx, col, dash, lbl) in [
        (72.0, "#2563eb", "", "Modeled"),
        (152.0, "#64748b", " stroke-dasharray=\"4,4\"", "DO Sat"),
    ] {
        let lx = p2_x + dx;
        let _ = writeln!(
            svg,
            "  <line x1=\"{lx:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"{col}\" stroke-width=\"2.0\"{dash}/>",
            leg1_y - 3.5,
            lx + 18.0,
            leg1_y - 3.5
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{leg1_y:.1}\" font-family=\"sans-serif\" font-size=\"10.5\" fill=\"#334155\">{lbl}</text>",
            lx + 22.0
        );
    }
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" transform=\"rotate(-90 {:.1} {:.1})\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#334155\">DO (mg/L)</text>",
        p2_x - 44.0,
        top + 0.5 * ph,
        p2_x - 44.0,
        top + 0.5 * ph
    );

    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Ci (ubar) — Vcmax={:.1}, Jmax={:.1}, TPU={:.1}</text>",
        p1_x + 0.5 * p1_w,
        height - 16.0,
        aci_fit.vcmax,
        aci_fit.jmax,
        aci_fit.tpu
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">Hour — GPP={:.2}, ER={:.2}, NEP={:+.2} g O2/m2/d</text>",
        p2_x + 0.5 * p2_w,
        height - 16.0,
        stream_rep.gpp_g_o2_m2_d,
        stream_rep.er_g_o2_m2_d,
        stream_rep.nep_g_o2_m2_d
    );

    svg.push_str("</svg>\n");
    svg
}


/// Built-in 25-hour diel dissolved-oxygen record from a forested Pacific Northwest stream
/// (`Oak Creek, Oregon`, `depth = 0.35 m`, true `GPP = 2.40`, `ER = 4.10`, `K600 = 7.5 d^-1`).
pub fn demo_stream_diel_do_csv() -> String {
    let mut out = String::from(
        "# Oak Creek (Corvallis, OR) 24-h Diel Dissolved Oxygen Logger (depth=0.35 m)\n\
         hour,do_mg_l,temp_c,par_umol_m2_s,depth_m\n",
    );
    let true_gpp = 2.40_f64;
    let true_er = 4.10_f64;
    let true_k600 = 7.50_f64;
    let depth = 0.35_f64;

    let mut hours = Vec::with_capacity(25);
    let mut temps = Vec::with_capacity(25);
    let mut pars = Vec::with_capacity(25);
    for h in 0..=24 {
        let hf = h as f64;
        let temp = 12.2 + 2.3 * ((std::f64::consts::PI * (hf - 9.0)) / 12.0).sin();
        let par = if (6..=18).contains(&h) {
            850.0 * ((std::f64::consts::PI * (hf - 6.0)) / 12.0).sin()
        } else {
            0.0
        };
        hours.push(hf);
        temps.push(temp);
        pars.push(par);
    }

    let mut sum_par_dt = 0.0_f64;
    for i in 1..=24 {
        sum_par_dt += 0.5 * (pars[i] + pars[i - 1]) * (1.0 / 24.0);
    }

    let mut do_curr = calc_do_saturation_mg_l(temps[0]) - 0.55;
    for i in 0..=24 {
        if i > 0 {
            let dt = 1.0 / 24.0;
            let t_mid = 0.5 * (temps[i] + temps[i - 1]);
            let par_mid = 0.5 * (pars[i] + pars[i - 1]);
            let sat_mid = 0.5 * (calc_do_saturation_mg_l(temps[i]) + calc_do_saturation_mg_l(temps[i - 1]));
            let k_t = true_k600 * (schmidt_o2(t_mid) / 600.0).powf(-0.5);
            let bio = (true_gpp * (par_mid / sum_par_dt) - true_er) / depth;
            do_curr = (do_curr + dt * (bio + k_t * sat_mid)) / (1.0 + dt * k_t);
        }
        out.push_str(&format!(
            "{:.1},{:.3},{:.2},{:.1},{:.2}\n",
            hours[i], do_curr, temps[i], pars[i], depth
        ));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_odum_stream_gpp_er_and_k600_accurately() {
        let pts = parse_stream_do_csv(&demo_stream_diel_do_csv()).unwrap();
        let rep = fit_odum_stream_metabolism(&pts).unwrap();
        assert!((rep.gpp_g_o2_m2_d - 2.40).abs() < 0.25);
        assert!((rep.er_g_o2_m2_d - 4.10).abs() < 0.25);
        assert!((rep.k600_per_day - 7.50).abs() < 0.80);
        assert!(rep.nep_g_o2_m2_d < -1.0);
    }

    #[test]
    fn pressure_elevation_q10_and_nan_safety_checks() {
        assert!((calc_do_saturation_mg_l_at_pressure(15.0, 1.0) - calc_do_saturation_mg_l(15.0)).abs() < 1e-12);
        let p_2000 = barometric_pressure_atm_at_elevation(2000.0);
        assert!(p_2000 > 0.75 && p_2000 < 0.82);
        let sat_sea = calc_do_saturation_mg_l(15.0);
        let sat_2000 = calc_do_saturation_mg_l_at_elevation(15.0, 2000.0);
        assert!(sat_2000 < sat_sea * 0.85);

        assert!((scale_er_q10(4.0, 25.0, 15.0, 2.0) - 8.0).abs() < 1e-12);
        let er_cold = scale_er_arrhenius(4.0, 10.0, 20.0, 45.0);
        let er_warm = scale_er_arrhenius(4.0, 25.0, 20.0, 45.0);
        assert!(er_cold < 4.0 && er_warm > 4.0);

        let pts = parse_stream_do_csv(&demo_stream_diel_do_csv()).unwrap();
        let cfg = StreamMetabolismConfig::at_elevation_m(
            1200.0,
            ErTempScaling::Q10 {
                q10: 2.0,
                ref_temp_c: 12.2,
            },
        );
        let rep_cfg = fit_odum_stream_metabolism_with_config(&pts, &cfg).unwrap();
        assert!(rep_cfg.gpp_g_o2_m2_d.is_finite() && rep_cfg.er_g_o2_m2_d.is_finite());

        // NaN safety in CSV parser and metabolism solver
        let bad_csv = "hour,do_mg_l,temp_c,par,depth_m\nNaN,9.0,12.0,0.0,0.35\n1,9.0,12.0,0.0,0.35\n";
        assert!(parse_stream_do_csv(bad_csv).is_err());

        let mut nan_pts = pts.clone();
        nan_pts[3].hour = f64::NAN;
        assert!(fit_odum_stream_metabolism(&nan_pts).is_err());
    }

    #[test]
    fn render_flux_svg_auto_scales_low_do_and_plots_wc_wj_wp_axes_and_legend() {
        let aci_pts = crate::flux_leaf::parse_aci_csv(&crate::flux_leaf::demo_fvcb_aci_csv()).unwrap();
        let aci_fit = crate::flux_leaf::fit_fvcb_aci_curve(&aci_pts).unwrap();

        let low_do_pts: Vec<StreamDoPoint> = (0..12)
            .map(|h| StreamDoPoint {
                hour: (h * 2) as f64,
                do_mg_l: 3.0,
                temp_c: 22.0,
                par_umol_m2_s: if (3..=9).contains(&h) { 500.0 } else { 0.0 },
                depth_m: 0.4,
            })
            .collect();
        let low_rep = StreamMetabolismReport {
            gpp_g_o2_m2_d: 1.1,
            er_g_o2_m2_d: 3.4,
            nep_g_o2_m2_d: -2.3,
            p_to_r_ratio: 1.1 / 3.4,
            k600_per_day: 4.5,
            trophic_state: String::from("Net Heterotrophic (P/R < 1)"),
            do_rmse_mg_l: 0.05,
            modeled_do_mg_l: vec![3.0; low_do_pts.len()],
            do_sat_mg_l: vec![3.0; low_do_pts.len()],
        };

        let svg = render_flux_svg(&aci_pts, &aci_fit, &low_do_pts, &low_rep, "Hypoxic Stream");
        assert!(svg.starts_with("<svg ") && svg.trim_end().ends_with("</svg>"));
        assert!(svg.contains("<line "), "Expected axis ticks and gridlines");
        assert!(svg.contains(">A_net<"));
        assert!(svg.contains(">Wc<"));
        assert!(svg.contains(">Wj<"));
        assert!(svg.contains(">Wp<"));
        assert!(svg.contains(">Obs DO<"));
        assert!(svg.contains(">Modeled<"));
        assert!(svg.contains(">DO Sat<"));
        assert!(svg.matches("<polyline ").count() >= 5);

        // Verify DO = 3.0 mg/L circles sit centered inside Panel 2 (top=74, ph=296 => mid y = 222.0, NOT bottom 370.0)
        assert!(
            svg.contains("cy=\"222.00\" r=\"3.5\" fill=\"#dc2626\""),
            "Expected auto-scaled DO=3.0 circles near panel center (y=222.00), svg:\n{svg}"
        );
    }
}
