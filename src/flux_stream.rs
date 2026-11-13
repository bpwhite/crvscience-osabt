// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #20 (`osabt flux`) — Odum (1956) Diel Dissolved-Oxygen Stream Metabolism Engine
//! & Two-Panel `.svg` Ecophysiology / Watershed Plotter:
//! - **Primary Production in Flowing Waters** — Odum (1956) *Limnol. Oceanogr.* 1(2):102–117
//!   (`https://www.emergysociety.com/wp-content/uploads/OdumHT.1956.Primary-production-in-flowing-waters.-Limnol-and-Oceanogr.V.1.pdf`)
//! - Estimates daily `GPP`, `ER`, `NEP = GPP - ER` (`g O2 m^-2 d^-1`), `P/R` ratio, and `K600` (`d^-1`).

use std::fmt::Write as _;

use crate::flux_leaf::{eval_fvcb_assimilation, FvcbAciFit, LeafAciPoint};
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

/// Computes Benson & Krause (1984) freshwater dissolved-oxygen saturation (`mg L^-1` at `1 atm`).
#[inline]
pub fn calc_do_saturation_mg_l(temp_c: f64) -> f64 {
    let tk = temp_c.clamp(0.0, 40.0) + 273.15;
    let ln_c = -139.344_11 + 1.575_701e5 / tk - 6.642_308e7 / (tk * tk)
        + 1.243_800e10 / (tk * tk * tk)
        - 8.621_949e11 / (tk * tk * tk * tk);
    ln_c.exp()
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
    for (idx, line) in raw_csv.lines().enumerate() {
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
    pts.sort_by(|a, b| a.hour.partial_cmp(&b.hour).unwrap());
    Ok(pts)
}

/// Fits Odum's (1956) open-channel oxygen balance (`GPP`, `ER`, `K600`) to diel `DO` observations
/// via Crank–Nicolson/midpoint ODE integration and Nelder–Mead SSE minimization.
pub fn fit_odum_stream_metabolism(points: &[StreamDoPoint]) -> Result<StreamMetabolismReport, String> {
    if points.len() < 8 {
        return Err(String::from("At least 8 diel points are required."));
    }
    let mut sum_par_dt = 0.0_f64;
    for i in 1..points.len() {
        let dt_days = ((points[i].hour - points[i - 1].hour) / 24.0).max(1e-4);
        let par_mid = 0.5 * (points[i].par_umol_m2_s + points[i - 1].par_umol_m2_s);
        sum_par_dt += par_mid * dt_days;
    }
    let par_norm = sum_par_dt.max(1e-6);
    let do_sat: Vec<f64> = points.iter().map(|p| calc_do_saturation_mg_l(p.temp_c)).collect();

    let simulate_do = |gpp: f64, er: f64, k600: f64| -> Vec<f64> {
        let mut pred = Vec::with_capacity(points.len());
        pred.push(points[0].do_mg_l);
        for i in 1..points.len() {
            let dt_days = ((points[i].hour - points[i - 1].hour) / 24.0).max(1e-4);
            let t_mid = 0.5 * (points[i].temp_c + points[i - 1].temp_c);
            let z_mid = 0.5 * (points[i].depth_m + points[i - 1].depth_m);
            let par_mid = 0.5 * (points[i].par_umol_m2_s + points[i - 1].par_umol_m2_s);
            let sat_mid = 0.5 * (do_sat[i] + do_sat[i - 1]);

            let k_t = k600 * (schmidt_o2(t_mid) / 600.0).powf(-0.5);
            let bio_flux = (gpp * (par_mid / par_norm) - er) / z_mid;
            // Semi-implicit step for unconditional stability across large K600
            let prev = pred[i - 1];
            let next = (prev + dt_days * (bio_flux + k_t * sat_mid)) / (1.0 + dt_days * k_t);
            pred.push(next);
        }
        pred
    };

    let sse = |th: &[f64]| -> f64 {
        let (gpp, er, k600) = (th[0], th[1], th[2]);
        if gpp < 0.0 || gpp > 40.0 || er < 0.05 || er > 40.0 || k600 < 0.2 || k600 > 60.0 {
            return 1e12 + gpp * gpp + er * er + k600 * k600;
        }
        let pred = simulate_do(gpp, er, k600);
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
    let (gpp, er, k600) = (best_th[0], best_th[1], best_th[2]);
    let modeled = simulate_do(gpp, er, k600);
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

/// Renders a two-panel `.svg` figure pairing the leaf `FvCB` $A\text{–}C_i$ curve (`left`)
/// with the diel stream dissolved-oxygen metabolism curve (`right`).
pub fn render_flux_svg(
    aci_pts: &[LeafAciPoint],
    aci_fit: &FvcbAciFit,
    do_pts: &[StreamDoPoint],
    stream_rep: &StreamMetabolismReport,
    title: &str,
) -> String {
    let width = 940.0_f64;
    let height = 430.0_f64;
    let p1_x = 68.0_f64;
    let p1_w = 390.0_f64;
    let p2_x = 538.0_f64;
    let p2_w = 352.0_f64;
    let top = 58.0_f64;
    let ph = 305.0_f64;

    let mut svg = String::with_capacity(14_336);
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
        "  <text x=\"{p1_x:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — FvCB Leaf A-Ci &amp; Odum Stream DO Metabolism</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p1_x:.1}\" y=\"{top:.1}\" width=\"{p1_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{p2_x:.1}\" y=\"{top:.1}\" width=\"{p2_w:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    // Panel 1: FvCB A-Ci curve
    let max_ci = aci_pts.last().map(|p| p.ci_ubar).unwrap_or(1600.0) * 1.05;
    let min_a = aci_pts.iter().map(|p| p.a_net).fold(0.0_f64, f64::min).min(-2.0);
    let max_a = aci_pts.iter().map(|p| p.a_net).fold(20.0_f64, f64::max) * 1.12;
    let map_x1 = |ci: f64| p1_x + (ci / max_ci).clamp(0.0, 1.0) * p1_w;
    let map_y1 = |a: f64| top + ph - ((a - min_a) / (max_a - min_a).max(1.0)).clamp(0.0, 1.0) * ph;

    let mut aci_poly = String::new();
    for k in 1..=60 {
        let ci = 45.0 + ((k as f64) / 60.0) * (max_ci - 45.0);
        let (a_val, _) = eval_fvcb_assimilation(
            ci,
            1500.0,
            25.0,
            aci_fit.vcmax,
            aci_fit.jmax,
            aci_fit.tpu,
            aci_fit.rd,
        );
        let _ = write!(aci_poly, "{:.1},{:.1} ", map_x1(ci), map_y1(a_val));
    }
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#16a34a\" stroke-width=\"2.4\"/>",
        aci_poly.trim()
    );
    for p in aci_pts {
        let _ = writeln!(
            svg,
            "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"4.0\" fill=\"#0f172a\"/>",
            map_x1(p.ci_ubar),
            map_y1(p.a_net)
        );
    }

    // Panel 2: Diel DO observed vs modeled + DO saturation
    let max_hr = do_pts.last().map(|p| p.hour).unwrap_or(24.0).max(1.0);
    let min_do = 6.0_f64;
    let max_do = 13.0_f64;
    let map_x2 = |h: f64| p2_x + (h / max_hr).clamp(0.0, 1.0) * p2_w;
    let map_y2 = |d: f64| top + ph - ((d - min_do) / (max_do - min_do)).clamp(0.0, 1.0) * ph;

    let mut mod_poly = String::new();
    let mut sat_poly = String::new();
    for (i, p) in do_pts.iter().enumerate() {
        let xp = map_x2(p.hour);
        let _ = write!(mod_poly, "{xp:.1},{:.1} ", map_y2(stream_rep.modeled_do_mg_l[i]));
        let _ = write!(sat_poly, "{xp:.1},{:.1} ", map_y2(stream_rep.do_sat_mg_l[i]));
        let _ = writeln!(
            svg,
            "  <circle cx=\"{xp:.1}\" cy=\"{:.1}\" r=\"3.5\" fill=\"#dc2626\"/>",
            map_y2(p.do_mg_l)
        );
    }
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#64748b\" stroke-dasharray=\"4,4\" stroke-width=\"1.5\"/>",
        sat_poly.trim()
    );
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.4\"/>",
        mod_poly.trim()
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

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
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
}
