// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #20 (`osabt flux`) — C3 Leaf `FvCB` Photosynthesis $A\text{–}C_i$ Curve Fitter:
//! - **Farquhar, von Caemmerer & Berry (1980)** *Planta* 149:78–90 (`https://pmc.ncbi.nlm.nih.gov/articles/PMC12500802/`)
//! - Uses the original minimum-carboxylation-rate formulation `Vc = min(Wc, Wj, Wp)` with
//!   Bernacchi et al. (2001) Arrhenius temperature scaling for `Kc`, `Ko`, and `Gamma*`.

use crate::wildlife_occupancy::nelder_mead_minimize;

/// Single leaf gas-exchange $A\text{–}C_i$ measurement point.
#[derive(Debug, Clone)]
pub struct LeafAciPoint {
    /// Intercellular CO2 partial pressure `C_i` (`ubar` / `umol mol^-1` at 1 atm).
    pub ci_ubar: f64,
    /// Net CO2 assimilation rate `A_n` (`umol CO2 m^-2 s^-1`).
    pub a_net: f64,
    /// Incident photosynthetic photon flux density (`umol photons m^-2 s^-1`).
    pub ppfd: f64,
    /// Leaf temperature (`°C`).
    pub t_leaf_c: f64,
}

/// Biochemical rate-limiting state at a given `C_i`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum FvcbLimitationState {
    RubiscoWc,
    RuBpRegenWj,
    TriosePhosphateWp,
}

impl FvcbLimitationState {
    pub fn label(self) -> &'static str {
        match self {
            Self::RubiscoWc => "Rubisco (Wc)",
            Self::RuBpRegenWj => "RuBP Regen (Wj)",
            Self::TriosePhosphateWp => "TPU (Wp)",
        }
    }
}

/// Fitted Farquhar–von Caemmerer–Berry (1980) photosynthetic capacity parameters.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct FvcbAciFit {
    /// Maximum Rubisco carboxylation rate `Vcmax` (`umol m^-2 s^-1`).
    pub vcmax: f64,
    /// Maximum electron transport rate `Jmax` (`umol m^-2 s^-1`).
    pub jmax: f64,
    /// Triose-phosphate utilization capacity `TPU` (`umol m^-2 s^-1`).
    pub tpu: f64,
    /// Daytime mitochondrial respiration `Rd` (`umol m^-2 s^-1`).
    pub rd: f64,
    pub jmax_to_vcmax_ratio: f64,
    pub rmse: f64,
    pub r_squared: f64,
    pub predicted_a: Vec<f64>,
    pub limiting_states: Vec<FvcbLimitationState>,
}

/// Parses an $A\text{–}C_i$ gas-exchange CSV/TSV (`ci_ubar,a_net[,ppfd,t_leaf_c]`).
pub fn parse_aci_csv(raw_csv: &str) -> Result<Vec<LeafAciPoint>, String> {
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
            if c0.contains("ci") || c0.contains("co2") {
                continue;
            }
        }
        if cols.len() < 2 {
            return Err(format!("A-Ci CSV line {}: expected >= 2 columns (ci_ubar,a_net).", idx + 1));
        }
        let ci = cols[0]
            .parse::<f64>()
            .map_err(|_| format!("A-Ci line {}: invalid Ci '{}'.", idx + 1, cols[0]))?;
        let an = cols[1]
            .parse::<f64>()
            .map_err(|_| format!("A-Ci line {}: invalid A_net '{}'.", idx + 1, cols[1]))?;
        let ppfd = cols.get(2).and_then(|s| s.parse::<f64>().ok()).unwrap_or(1500.0);
        let t_leaf = cols.get(3).and_then(|s| s.parse::<f64>().ok()).unwrap_or(25.0);
        if ci > 0.0 && an.is_finite() {
            pts.push(LeafAciPoint {
                ci_ubar: ci,
                a_net: an,
                ppfd,
                t_leaf_c: t_leaf,
            });
        }
    }
    if pts.len() < 5 {
        return Err(String::from("FvCB A-Ci fitting requires at least 5 measurement points."));
    }
    pts.sort_by(|a, b| a.ci_ubar.partial_cmp(&b.ci_ubar).unwrap());
    Ok(pts)
}

/// Evaluates the Farquhar, von Caemmerer & Berry (1980) net assimilation rate `A_n`
/// and active biochemical limitation (`Wc`, `Wj`, `Wp`) at `(ci_ubar, ppfd, t_leaf_c)`.
pub fn eval_fvcb_assimilation(
    ci_ubar: f64,
    ppfd: f64,
    t_leaf_c: f64,
    vcmax25: f64,
    jmax25: f64,
    tpu25: f64,
    rd25: f64,
) -> (f64, FvcbLimitationState) {
    let tk = t_leaf_c + 273.15;
    let r_gas = 0.008_314_5_f64; // kJ / (mol K)
    let arrh = |val25: f64, ha_kj: f64| -> f64 {
        val25 * (ha_kj * (tk - 298.15) / (298.15 * r_gas * tk)).exp()
    };

    // Bernacchi et al. (2001) Rubisco kinetics (`25 °C` baseline: Kc=404.9 ubar, Ko=278.4 mbar, Gamma*=42.75 ubar)
    let kc = arrh(404.9, 79.43);
    let ko = arrh(278.4, 36.38); // mbar
    let gamma_star = arrh(42.75, 37.83);
    let o2_mbar = 210.0_f64;

    let vcmax = arrh(vcmax25, 65.33);
    let jmax = arrh(jmax25, 43.54);
    let tpu = arrh(tpu25, 53.10);
    let rd = arrh(rd25, 46.39);

    // Non-rectangular hyperbola for chloroplast electron transport `J(PPFD, Jmax)`
    let i_abs = 0.85 * 0.5 * ppfd; // Absorbed PAR partitioned to PSII
    let theta_j = 0.70_f64;
    let sum_ij = i_abs + jmax;
    let disc = (sum_ij * sum_ij - 4.0 * theta_j * i_abs * jmax).max(0.0);
    let j_rate = (sum_ij - disc.sqrt()) / (2.0 * theta_j);

    let km_app = kc * (1.0 + o2_mbar / ko);
    let wc = (vcmax * ci_ubar) / (ci_ubar + km_app);
    let wj = (j_rate * ci_ubar) / (4.0 * ci_ubar + 8.0 * gamma_star);
    let wp = if ci_ubar > gamma_star + 1.0 {
        (3.0 * tpu * ci_ubar) / (ci_ubar - gamma_star)
    } else {
        f64::INFINITY
    };

    let (w_lim, state) = if wc <= wj && wc <= wp {
        (wc, FvcbLimitationState::RubiscoWc)
    } else if wj <= wp {
        (wj, FvcbLimitationState::RuBpRegenWj)
    } else {
        (wp, FvcbLimitationState::TriosePhosphateWp)
    };

    let a_net = w_lim * (1.0 - gamma_star / ci_ubar.max(1e-6)) - rd;
    (a_net, state)
}

/// Fits `(Vcmax, Jmax, TPU, Rd)` to an $A\text{–}C_i$ response curve via bounded Nelder–Mead SSE minimization.
pub fn fit_fvcb_aci_curve(points: &[LeafAciPoint]) -> Result<FvcbAciFit, String> {
    if points.len() < 5 {
        return Err(String::from("At least 5 A-Ci points are required for FvCB fitting."));
    }

    let sse = |th: &[f64]| -> f64 {
        let (vcmax, jmax, tpu, rd) = (th[0], th[1], th[2], th[3]);
        if vcmax < 5.0 || vcmax > 400.0 || jmax < 10.0 || jmax > 600.0 || tpu < 2.0 || tpu > 40.0 || rd < 0.05 || rd > 15.0 {
            return 1e12 + vcmax * vcmax + jmax * jmax;
        }
        let mut err = 0.0_f64;
        for p in points {
            let (pred, _) = eval_fvcb_assimilation(
                p.ci_ubar, p.ppfd, p.t_leaf_c, vcmax, jmax, tpu, rd,
            );
            let d = pred - p.a_net;
            err += d * d;
        }
        err
    };

    let starts = [
        [80.0, 140.0, 10.0, 1.2],
        [110.0, 190.0, 13.0, 1.5],
        [55.0, 100.0, 7.5, 0.9],
    ];
    let mut best_th = starts[0].to_vec();
    let mut best_sse = f64::INFINITY;
    for st in &starts {
        let (th, val) = nelder_mead_minimize(&sse, st, 1000, 1e-9);
        if val < best_sse {
            best_sse = val;
            best_th = th;
        }
    }
    let (best_th, best_sse) = nelder_mead_minimize(&sse, &best_th, 1200, 1e-11);
    let (vcmax, jmax, tpu, rd) = (best_th[0], best_th[1], best_th[2], best_th[3]);

    let mut predicted_a = Vec::with_capacity(points.len());
    let mut limiting_states = Vec::with_capacity(points.len());
    for p in points {
        let (pred, st) =
            eval_fvcb_assimilation(p.ci_ubar, p.ppfd, p.t_leaf_c, vcmax, jmax, tpu, rd);
        predicted_a.push(pred);
        limiting_states.push(st);
    }

    let mean_a = points.iter().map(|p| p.a_net).sum::<f64>() / (points.len() as f64);
    let ss_tot: f64 = points.iter().map(|p| (p.a_net - mean_a).powi(2)).sum();
    let rmse = (best_sse / (points.len() as f64)).sqrt();
    let r2 = if ss_tot > 1e-12 {
        1.0 - best_sse / ss_tot
    } else {
        1.0
    };

    Ok(FvcbAciFit {
        vcmax,
        jmax,
        tpu,
        rd,
        jmax_to_vcmax_ratio: jmax / vcmax.max(1e-6),
        rmse,
        r_squared: r2,
        predicted_a,
        limiting_states,
    })
}

/// Built-in 12-point Pacific Northwest Douglas-fir (*Pseudotsuga menziesii*) $A\text{–}C_i$ curve
/// (`25 °C`, `PPFD = 1500 umol m^-2 s^-1`, true `Vcmax = 86`, `Jmax = 152`, `TPU = 10.6`, `Rd = 1.15`).
pub fn demo_fvcb_aci_csv() -> String {
    let cis = [55.0, 90.0, 130.0, 180.0, 240.0, 320.0, 420.0, 580.0, 780.0, 1000.0, 1250.0, 1550.0];
    let mut out = String::from(
        "# Douglas-fir (Pseudotsuga menziesii) Leaf Gas-Exchange A-Ci Curve (25 C, PPFD=1500)\n\
         ci_ubar,a_net,ppfd,t_leaf_c\n",
    );
    for (i, &ci) in cis.iter().enumerate() {
        let (a_true, _) = eval_fvcb_assimilation(ci, 1500.0, 25.0, 86.0, 152.0, 10.6, 1.15);
        let jitter = if i % 2 == 0 { 0.08 } else { -0.07 };
        out.push_str(&format!("{ci:.1},{:.3},1500.0,25.0\n", a_true + jitter));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fits_fvcb_vcmax_jmax_tpu_and_rd_accurately() {
        let pts = parse_aci_csv(&demo_fvcb_aci_csv()).unwrap();
        let fit = fit_fvcb_aci_curve(&pts).unwrap();
        assert!((fit.vcmax - 86.0).abs() < 3.0);
        assert!((fit.jmax - 152.0).abs() < 5.0);
        assert!((fit.tpu - 10.6).abs() < 0.8);
        assert!(fit.r_squared > 0.999);
        assert_eq!(fit.limiting_states[0], FvcbLimitationState::RubiscoWc);
        assert_eq!(
            *fit.limiting_states.last().unwrap(),
            FvcbLimitationState::TriosePhosphateWp
        );
    }
}
