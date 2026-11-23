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

//! Tool #13 (`osabt kinetics`) — Steady-state enzyme kinetics, global reversible inhibition
//! mechanism discrimination (`AICc` + Akaike weights), Cheng–Prusoff conversion, and
//! Morrison tight-binding quadratic solver:
//! - **Michaelis–Menten Kinetics** — Michaelis & Menten (1913) / Johnson & Goody (2011) *Biochemistry* 50:8264–8269
//! - **Tight-Binding Quadratic Kinetics** — Morrison (1969) *Biochim. Biophys. Acta* 185:269–286;
//!   Williams & Morrison (1979) *Methods Enzymol.* 63:437–467
//! - **`Ki`–`IC50` Relationships** — Cheng & Prusoff (1973) *Biochem. Pharmacol.* 22:3099–3108
//! - **Multimodel `AICc` Selection** — Burnham & Anderson (2004) *Sociol. Methods Res.* 33:261–304

/// Single steady-state initial velocity observation `([S] µM, [I] nM, v0)`.
#[derive(Debug, Clone, Copy)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct KineticsPoint {
    pub substrate_um: f64,
    pub inhibitor_nm: f64,
    pub velocity: f64,
}

/// Candidate single-substrate steady-state rate models (`[I] == 0`).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SingleSubstrateModel {
    MichaelisMenten,
    AllostericHill,
    SubstrateInhibition,
}

impl SingleSubstrateModel {
    pub fn label(self) -> &'static str {
        match self {
            Self::MichaelisMenten => "Michaelis–Menten",
            Self::AllostericHill => "Allosteric Hill (Cooperative)",
            Self::SubstrateInhibition => "Haldane Substrate Inhibition",
        }
    }
}

/// Nonlinear least-squares fit + `AICc` score for a single-substrate `[S] -> v0` curve.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SingleSubstrateFit {
    pub model: SingleSubstrateModel,
    pub n_points: usize,
    pub vmax: f64,
    pub km_or_khalf_um: f64,
    /// Hill coefficient `n_H` (for `AllostericHill`) or substrate-inhibition constant `K_si` (`µM`).
    pub secondary_param: Option<f64>,
    /// Turnover number `kcat` (`s^-1`) when total active enzyme concentration `[E]_T` (`nM`) is supplied
    /// and velocities are in `µM/min` (`kcat = (Vmax_uM_min / 60) / (Et_nM * 1e-3)`).
    pub kcat_per_s: Option<f64>,
    /// Specificity constant `kcat / Km` (`M^-1 s^-1`).
    pub catalytic_efficiency_m_inv_s_inv: Option<f64>,
    pub sse: f64,
    pub r_squared: f64,
    pub aicc: f64,
    pub akaike_weight: f64,
}

/// Classical reversible enzyme inhibition modalities discriminated via global 2D fitting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum InhibitionMechanism {
    Competitive,
    Uncompetitive,
    PureNoncompetitive,
    MixedType,
}

impl InhibitionMechanism {
    pub fn label(self) -> &'static str {
        match self {
            Self::Competitive => "Competitive (binds free E)",
            Self::Uncompetitive => "Uncompetitive (binds ES complex)",
            Self::PureNoncompetitive => "Pure Noncompetitive (Ki = Ki')",
            Self::MixedType => "Mixed-Type (Ki != α·Ki)",
        }
    }
}

/// Global 2D `([S], [I]) -> v0` nonlinear fit ranked by Burnham & Anderson (2004) `AICc`.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GlobalInhibitionFit {
    pub mechanism: InhibitionMechanism,
    pub n_points: usize,
    pub num_params_k: usize,
    pub vmax: f64,
    pub km_um: f64,
    /// Dissociation constant `Ki` (`nM`) for free enzyme `E` (or `Ki'` for `Uncompetitive`).
    pub ki_nm: f64,
    /// Cooperativity factor `alpha` (`Ki' = alpha * Ki`) for `MixedType`.
    pub alpha_factor: Option<f64>,
    pub ki_prime_nm: Option<f64>,
    pub sse: f64,
    pub r_squared: f64,
    pub aicc: f64,
    pub delta_aicc: f64,
    pub akaike_weight: f64,
}

/// Morrison (1969) / Williams & Morrison (1979) tight-binding quadratic fit (`[E]_T ~ Ki_app`).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MorrisonTightBindingFit {
    pub n_points: usize,
    pub v0_unhibited: f64,
    pub active_enzyme_et_nm: f64,
    pub ki_app_nm: f64,
    /// Intrinsic `Ki` (`nM`) when competitive `[S]` and `Km` (`µM`) are provided (`Ki = Ki_app / (1 + [S]/Km)`).
    pub intrinsic_ki_nm: Option<f64>,
    /// Fractional active enzyme purity `[E]_T / [E]_nominal` when nominal protein concentration is known.
    pub fractional_active_enzyme: Option<f64>,
    pub sse: f64,
    pub r_squared: f64,
}

/// Cheng & Prusoff (1973) analytical conversion from `IC50` (`nM`) to `Ki` (`nM`).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ChengPrusoffConversion {
    pub mechanism: InhibitionMechanism,
    pub ic50_nm: f64,
    pub substrate_um: f64,
    pub km_um: f64,
    pub computed_ki_nm: f64,
}

/// Parses a CSV or TSV kinetics table (`substrate_um,inhibitor_nm,velocity` or 2-column `concentration,velocity`).
pub fn parse_kinetics_csv(raw: &str) -> Result<Vec<KineticsPoint>, String> {
    let mut pts = Vec::new();
    let mut header_checked = false;

    for (idx, line) in raw.trim_start_matches('\u{feff}').lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();
        if !header_checked {
            header_checked = true;
            if cols[0]
                .chars()
                .next()
                .map(|c| c.is_ascii_alphabetic() || c == '[')
                .unwrap_or(false)
            {
                continue;
            }
        }
        if cols.len() == 2 {
            let s: f64 = cols[0]
                .parse()
                .map_err(|_| format!("Invalid concentration on line {}.", idx + 1))?;
            let v: f64 = cols[1]
                .parse()
                .map_err(|_| format!("Invalid velocity on line {}.", idx + 1))?;
            pts.push(KineticsPoint {
                substrate_um: s,
                inhibitor_nm: 0.0,
                velocity: v,
            });
        } else if cols.len() >= 3 {
            let s: f64 = cols[0]
                .parse()
                .map_err(|_| format!("Invalid substrate_um on line {}.", idx + 1))?;
            let i: f64 = cols[1]
                .parse()
                .map_err(|_| format!("Invalid inhibitor_nm on line {}.", idx + 1))?;
            let v: f64 = cols[2]
                .parse()
                .map_err(|_| format!("Invalid velocity on line {}.", idx + 1))?;
            pts.push(KineticsPoint {
                substrate_um: s,
                inhibitor_nm: i,
                velocity: v,
            });
        } else {
            return Err(format!(
                "Line {}: expected 2 or 3 numeric columns.",
                idx + 1
            ));
        }
    }

    if pts.len() < 4 {
        return Err(String::from("Kinetics fitting requires at least 4 data points."));
    }
    Ok(pts)
}

/// Evaluates initial velocity `v0` under the specified reversible inhibition mechanism.
#[inline]
pub fn eval_inhibition_velocity(
    mech: InhibitionMechanism,
    s_um: f64,
    i_nm: f64,
    vmax: f64,
    km_um: f64,
    ki_nm: f64,
    alpha: f64,
) -> f64 {
    let s = s_um.max(0.0);
    let i = i_nm.max(0.0);
    let ki = ki_nm.max(1e-12);
    match mech {
        InhibitionMechanism::Competitive => (vmax * s) / (km_um * (1.0 + i / ki) + s).max(1e-18),
        InhibitionMechanism::Uncompetitive => (vmax * s) / (km_um + s * (1.0 + i / ki)).max(1e-18),
        InhibitionMechanism::PureNoncompetitive => {
            (vmax * s) / ((km_um + s) * (1.0 + i / ki)).max(1e-18)
        }
        InhibitionMechanism::MixedType => {
            let ki_p = (alpha.max(1e-6) * ki).max(1e-12);
            (vmax * s) / (km_um * (1.0 + i / ki) + s * (1.0 + i / ki_p)).max(1e-18)
        }
    }
}

/// Evaluates single-substrate steady-state models (`MichaelisMenten`, `AllostericHill`, `SubstrateInhibition`).
#[inline]
pub fn eval_single_substrate_velocity(
    model: SingleSubstrateModel,
    s_um: f64,
    vmax: f64,
    km_um: f64,
    sec: f64,
) -> f64 {
    let s = s_um.max(0.0);
    if s == 0.0 {
        return 0.0;
    }
    match model {
        SingleSubstrateModel::MichaelisMenten => (vmax * s) / (km_um + s).max(1e-18),
        SingleSubstrateModel::AllostericHill => {
            let nh = sec.clamp(0.15, 8.0);
            let sn = s.powf(nh);
            let kn = km_um.max(1e-12).powf(nh);
            (vmax * sn) / (kn + sn).max(1e-18)
        }
        SingleSubstrateModel::SubstrateInhibition => {
            let ksi = sec.max(1e-6);
            (vmax * s) / (km_um + s * (1.0 + s / ksi)).max(1e-18)
        }
    }
}

/// Fits Michaelis–Menten, Allosteric Hill, and Haldane Substrate-Inhibition models to zero-inhibitor
/// (`[I] == 0`) points and ranks them via `AICc` and Akaike weights.
pub fn fit_single_substrate_models(
    points: &[KineticsPoint],
    enzyme_conc_nm: Option<f64>,
) -> Result<Vec<SingleSubstrateFit>, String> {
    let uninh: Vec<KineticsPoint> = points
        .iter()
        .copied()
        .filter(|p| {
            p.inhibitor_nm.is_finite()
                && p.inhibitor_nm.abs() <= 1e-9
                && p.substrate_um.is_finite()
                && p.substrate_um >= 0.0
                && p.velocity.is_finite()
        })
        .collect();
    let data = if uninh.len() >= 5 {
        uninh
    } else {
        points
            .iter()
            .copied()
            .filter(|p| p.substrate_um.is_finite() && p.substrate_um >= 0.0 && p.velocity.is_finite())
            .collect()
    };
    if data.len() < 5 {
        return Err(String::from(
            "Single-substrate model comparison requires >= 5 finite velocity points.",
        ));
    }

    let vmax_init = data
        .iter()
        .map(|p| p.velocity)
        .fold(1.0_f64, f64::max)
        * 1.15;
    let half_v = 0.5 * vmax_init;
    let km_init = data
        .iter()
        .min_by(|a, b| (a.velocity - half_v).abs().total_cmp(&(b.velocity - half_v).abs()))
        .map(|p| p.substrate_um.max(1.0))
        .unwrap_or(25.0);
    let max_s = data
        .iter()
        .map(|p| p.substrate_um)
        .fold(100.0_f64, f64::max);

    let ys: Vec<f64> = data.iter().map(|p| p.velocity).collect();
    let sst = total_sum_of_squares(&ys);
    let n = data.len();

    // 1. Michaelis–Menten (p = 2, K = 3)
    let mm_theta = lm_optimize_log_params(
        &[vmax_init.ln(), km_init.ln()],
        &ys,
        |th, idx| {
            eval_single_substrate_velocity(
                SingleSubstrateModel::MichaelisMenten,
                data[idx].substrate_um,
                th[0].exp(),
                th[1].exp(),
                1.0,
            )
        },
    );
    let (mm_vmax, mm_km) = (mm_theta[0].exp(), mm_theta[1].exp());
    let mm_sse: f64 = data
        .iter()
        .map(|p| {
            let r = p.velocity
                - eval_single_substrate_velocity(
                    SingleSubstrateModel::MichaelisMenten,
                    p.substrate_um,
                    mm_vmax,
                    mm_km,
                    1.0,
                );
            r * r
        })
        .sum();

    // 2. Allosteric Hill (p = 3, K = 4)
    let hill_theta = lm_optimize_log_params(
        &[mm_vmax.ln(), mm_km.ln(), 0.0],
        &ys,
        |th, idx| {
            eval_single_substrate_velocity(
                SingleSubstrateModel::AllostericHill,
                data[idx].substrate_um,
                th[0].exp(),
                th[1].exp(),
                th[2].exp().clamp(0.15, 8.0),
            )
        },
    );
    let (h_vmax, h_km, h_n) = (
        hill_theta[0].exp(),
        hill_theta[1].exp(),
        hill_theta[2].exp().clamp(0.15, 8.0),
    );
    let h_sse: f64 = data
        .iter()
        .map(|p| {
            let r = p.velocity
                - eval_single_substrate_velocity(
                    SingleSubstrateModel::AllostericHill,
                    p.substrate_um,
                    h_vmax,
                    h_km,
                    h_n,
                );
            r * r
        })
        .sum();

    // 3. Substrate Inhibition (p = 3, K = 4)
    let si_theta = lm_optimize_log_params(
        &[mm_vmax.ln(), mm_km.ln(), (max_s * 3.0).ln()],
        &ys,
        |th, idx| {
            eval_single_substrate_velocity(
                SingleSubstrateModel::SubstrateInhibition,
                data[idx].substrate_um,
                th[0].exp(),
                th[1].exp(),
                th[2].exp(),
            )
        },
    );
    let (si_vmax, si_km, si_ksi) = (si_theta[0].exp(), si_theta[1].exp(), si_theta[2].exp());
    let si_sse: f64 = data
        .iter()
        .map(|p| {
            let r = p.velocity
                - eval_single_substrate_velocity(
                    SingleSubstrateModel::SubstrateInhibition,
                    p.substrate_um,
                    si_vmax,
                    si_km,
                    si_ksi,
                );
            r * r
        })
        .sum();

    let build_fit =
        |model: SingleSubstrateModel, vmax: f64, km: f64, sec: Option<f64>, sse: f64, k_par: usize| {
            let (kcat, cat_eff) = if let Some(et_nm) = enzyme_conc_nm {
                if et_nm > 0.0 {
                    // Vmax (µM/min) -> µM/s = Vmax/60; [E]_T (µM) = Et_nM * 1e-3
                    let kc = (vmax / 60.0) / (et_nm * 1e-3);
                    // kcat / Km (M^-1 s^-1) where Km is in µM = Km * 1e-6 M
                    let ce = kc / (km * 1e-6).max(1e-18);
                    (Some(kc), Some(ce))
                } else {
                    (None, None)
                }
            } else {
                (None, None)
            };
            SingleSubstrateFit {
                model,
                n_points: n,
                vmax,
                km_or_khalf_um: km,
                secondary_param: sec,
                kcat_per_s: kcat,
                catalytic_efficiency_m_inv_s_inv: cat_eff,
                sse,
                r_squared: if sst > 1e-12 {
                    (1.0 - sse / sst).clamp(0.0, 1.0)
                } else {
                    1.0
                },
                aicc: compute_aicc(n, sse, k_par),
                akaike_weight: 0.0,
            }
        };

    let mut fits = vec![
        build_fit(SingleSubstrateModel::MichaelisMenten, mm_vmax, mm_km, None, mm_sse, 3),
        build_fit(
            SingleSubstrateModel::AllostericHill,
            h_vmax,
            h_km,
            Some(h_n),
            h_sse,
            4,
        ),
        build_fit(
            SingleSubstrateModel::SubstrateInhibition,
            si_vmax,
            si_km,
            Some(si_ksi),
            si_sse,
            4,
        ),
    ];

    fits.sort_by(|a, b| a.aicc.total_cmp(&b.aicc));
    let min_aicc = fits[0].aicc;
    let w_sum: f64 = fits
        .iter()
        .map(|f| (-0.5 * (f.aicc - min_aicc)).exp())
        .sum::<f64>()
        .max(1e-30);
    for f in &mut fits {
        f.akaike_weight = (-0.5 * (f.aicc - min_aicc)).exp() / w_sum;
    }

    Ok(fits)
}

/// Globally fits all 4 classical reversible inhibition mechanisms (`Competitive`, `Uncompetitive`,
/// `PureNoncompetitive`, `MixedType`) across a 2D `([S], [I]) -> v0` matrix and ranks them via `AICc`.
pub fn fit_global_reversible_inhibition(
    points: &[KineticsPoint],
) -> Result<Vec<GlobalInhibitionFit>, String> {
    if points.len() < 8 {
        return Err(String::from(
            "Global 2D inhibition mechanism selection requires >= 8 ([S], [I], v0) observations.",
        ));
    }
    let has_inh = points.iter().any(|p| p.inhibitor_nm > 0.0);
    if !has_inh {
        return Err(String::from(
            "Dataset contains no positive inhibitor ([I] > 0) points for Ki estimation.",
        ));
    }

    let ys: Vec<f64> = points.iter().map(|p| p.velocity).collect();
    let sst = total_sum_of_squares(&ys);
    let n = points.len();

    let vmax0 = points
        .iter()
        .map(|p| p.velocity)
        .fold(1.0_f64, f64::max)
        * 1.15;
    let km0 = median_positive(points.iter().map(|p| p.substrate_um).collect()).max(1.0);
    let ki0 = median_positive(
        points
            .iter()
            .filter(|p| p.inhibitor_nm > 0.0)
            .map(|p| p.inhibitor_nm)
            .collect(),
    )
    .max(1.0);

    let mut fits = Vec::with_capacity(4);
    for mech in [
        InhibitionMechanism::Competitive,
        InhibitionMechanism::Uncompetitive,
        InhibitionMechanism::PureNoncompetitive,
    ] {
        let th = lm_optimize_log_params(&[vmax0.ln(), km0.ln(), ki0.ln()], &ys, |p, idx| {
            eval_inhibition_velocity(
                mech,
                points[idx].substrate_um,
                points[idx].inhibitor_nm,
                p[0].exp(),
                p[1].exp(),
                p[2].exp(),
                1.0,
            )
        });
        let (vmax, km, ki) = (th[0].exp(), th[1].exp(), th[2].exp());
        let sse: f64 = points
            .iter()
            .map(|pt| {
                let r = pt.velocity
                    - eval_inhibition_velocity(
                        mech,
                        pt.substrate_um,
                        pt.inhibitor_nm,
                        vmax,
                        km,
                        ki,
                        1.0,
                    );
                r * r
            })
            .sum();
        fits.push(GlobalInhibitionFit {
            mechanism: mech,
            n_points: n,
            num_params_k: 4,
            vmax,
            km_um: km,
            ki_nm: ki,
            alpha_factor: None,
            ki_prime_nm: None,
            sse,
            r_squared: (1.0 - sse / sst.max(1e-12)).clamp(0.0, 1.0),
            aicc: compute_aicc(n, sse, 4),
            delta_aicc: 0.0,
            akaike_weight: 0.0,
        });
    }

    // Seed Mixed-Type from best 3-parameter mechanism so global minimum is guaranteed
    let best_3p = fits
        .iter()
        .min_by(|a, b| a.sse.total_cmp(&b.sse))
        .unwrap();
    let init_alpha = match best_3p.mechanism {
        InhibitionMechanism::Competitive => 20.0_f64,
        InhibitionMechanism::Uncompetitive => 0.05_f64,
        _ => 1.0_f64,
    };
    let init_ki_mixed = if best_3p.mechanism == InhibitionMechanism::Uncompetitive {
        best_3p.ki_nm / init_alpha
    } else {
        best_3p.ki_nm
    };
    let m_th = lm_optimize_log_params(
        &[
            best_3p.vmax.ln(),
            best_3p.km_um.ln(),
            init_ki_mixed.ln(),
            init_alpha.ln(),
        ],
        &ys,
        |p, idx| {
            eval_inhibition_velocity(
                InhibitionMechanism::MixedType,
                points[idx].substrate_um,
                points[idx].inhibitor_nm,
                p[0].exp(),
                p[1].exp(),
                p[2].exp(),
                p[3].exp().clamp(1e-4, 1e5),
            )
        },
    );
    let (mv, mkm, mki, malpha) = (
        m_th[0].exp(),
        m_th[1].exp(),
        m_th[2].exp(),
        m_th[3].exp().clamp(1e-4, 1e5),
    );
    let msse: f64 = points
        .iter()
        .map(|pt| {
            let r = pt.velocity
                - eval_inhibition_velocity(
                    InhibitionMechanism::MixedType,
                    pt.substrate_um,
                    pt.inhibitor_nm,
                    mv,
                    mkm,
                    mki,
                    malpha,
                );
            r * r
        })
        .sum();
    fits.push(GlobalInhibitionFit {
        mechanism: InhibitionMechanism::MixedType,
        n_points: n,
        num_params_k: 5,
        vmax: mv,
        km_um: mkm,
        ki_nm: mki,
        alpha_factor: Some(malpha),
        ki_prime_nm: Some(malpha * mki),
        sse: msse,
        r_squared: (1.0 - msse / sst.max(1e-12)).clamp(0.0, 1.0),
        aicc: compute_aicc(n, msse, 5),
        delta_aicc: 0.0,
        akaike_weight: 0.0,
    });

    fits.sort_by(|a, b| a.aicc.total_cmp(&b.aicc));
    let best_aicc = fits[0].aicc;
    let w_sum: f64 = fits
        .iter()
        .map(|f| (-0.5 * (f.aicc - best_aicc)).exp())
        .sum::<f64>()
        .max(1e-30);
    for f in &mut fits {
        f.delta_aicc = f.aicc - best_aicc;
        f.akaike_weight = (-0.5 * f.delta_aicc).exp() / w_sum;
    }

    Ok(fits)
}

/// Evaluates Morrison's (1969) / Williams & Morrison's (1979) quadratic rate law for tight-binding inhibition:
/// `v_i = v_0 * (1 - ((Et + It + Ki_app) - sqrt((Et + It + Ki_app)^2 - 4 * Et * It)) / (2 * Et))`.
#[inline]
pub fn eval_morrison_velocity(i_nm: f64, v0: f64, et_nm: f64, ki_app_nm: f64) -> f64 {
    let it = i_nm.max(0.0);
    let et = et_nm.max(1e-9);
    let kia = ki_app_nm.max(1e-9);
    let b = et + it + kia;
    let disc = (b * b - 4.0 * et * it).max(0.0).sqrt();
    let frac_bound = ((b - disc) / (2.0 * et)).clamp(0.0, 1.0);
    v0 * (1.0 - frac_bound)
}

/// Fits Morrison's (1969) tight-binding quadratic model to `(inhibitor_nm, velocity)` observations,
/// either jointly estimating active enzyme concentration `[E]_T` (active-site titration) or holding `[E]_T` fixed.
pub fn fit_morrison_tight_binding(
    observations: &[(f64, f64)],
    fixed_et_nm: Option<f64>,
    nominal_et_nm: Option<f64>,
    competitive_s_and_km_um: Option<(f64, f64)>,
) -> Result<MorrisonTightBindingFit, String> {
    if observations.len() < 5 {
        return Err(String::from(
            "Morrison tight-binding regression requires >= 5 ([I], v_i) points.",
        ));
    }

    let ys: Vec<f64> = observations.iter().map(|p| p.1).collect();
    let sst = total_sum_of_squares(&ys);
    let v0_init = observations
        .iter()
        .filter(|p| p.0 <= 1e-9)
        .map(|p| p.1)
        .fold(0.0_f64, f64::max)
        .max(ys.iter().copied().fold(1.0_f64, f64::max));

    let half_v = 0.5 * v0_init;
    let ic50_approx = observations
        .iter()
        .min_by(|a, b| (a.1 - half_v).abs().total_cmp(&(b.1 - half_v).abs()))
        .map(|p| p.0.max(0.5))
        .unwrap_or(10.0);

    let (v0_fit, et_fit, kia_fit) = if let Some(et_fixed) = fixed_et_nm {
        let kia0 = (ic50_approx - 0.5 * et_fixed).max(0.05);
        let th = lm_optimize_log_params(&[v0_init.ln(), kia0.ln()], &ys, |p, idx| {
            eval_morrison_velocity(observations[idx].0, p[0].exp(), et_fixed, p[1].exp())
        });
        (th[0].exp(), et_fixed, th[1].exp())
    } else {
        let et0 = nominal_et_nm.unwrap_or(1.4 * ic50_approx).max(0.2);
        let kia0 = (0.15 * ic50_approx).max(0.05);
        let th = lm_optimize_log_params(
            &[v0_init.ln(), et0.ln(), kia0.ln()],
            &ys,
            |p, idx| {
                eval_morrison_velocity(
                    observations[idx].0,
                    p[0].exp(),
                    p[1].exp(),
                    p[2].exp(),
                )
            },
        );
        (th[0].exp(), th[1].exp(), th[2].exp())
    };

    let sse: f64 = observations
        .iter()
        .map(|&(it, vy)| {
            let r = vy - eval_morrison_velocity(it, v0_fit, et_fit, kia_fit);
            r * r
        })
        .sum();

    let intrinsic_ki = competitive_s_and_km_um.and_then(|(s_um, km_um)| {
        if km_um > 0.0 && s_um >= 0.0 {
            Some(kia_fit / (1.0 + s_um / km_um))
        } else {
            None
        }
    });
    let frac_active = nominal_et_nm.and_then(|nom| if nom > 0.0 { Some(et_fit / nom) } else { None });

    Ok(MorrisonTightBindingFit {
        n_points: observations.len(),
        v0_unhibited: v0_fit,
        active_enzyme_et_nm: et_fit,
        ki_app_nm: kia_fit,
        intrinsic_ki_nm: intrinsic_ki,
        fractional_active_enzyme: frac_active,
        sse,
        r_squared: (1.0 - sse / sst.max(1e-12)).clamp(0.0, 1.0),
    })
}

/// Converts a measured `IC50` (`nM`) into an inhibition constant `Ki` (`nM`) using
/// Cheng & Prusoff (1973) equations for Classical (`[E]_T << IC50`) inhibitors.
pub fn convert_ic50_cheng_prusoff(
    mechanism: InhibitionMechanism,
    ic50_nm: f64,
    substrate_um: f64,
    km_um: f64,
) -> ChengPrusoffConversion {
    let s = substrate_um.max(0.0);
    let km = km_um.max(1e-9);
    let ki = match mechanism {
        InhibitionMechanism::Competitive => ic50_nm / (1.0 + s / km),
        InhibitionMechanism::Uncompetitive => ic50_nm / (1.0 + km / s.max(1e-9)),
        InhibitionMechanism::PureNoncompetitive | InhibitionMechanism::MixedType => ic50_nm,
    };
    ChengPrusoffConversion {
        mechanism,
        ic50_nm,
        substrate_um: s,
        km_um: km,
        computed_ki_nm: ki,
    }
}

/// Built-in 28-point global 2D kinase kinetics matrix (`[ATP] = 5..500 µM`, `[I] = 0, 15, 45, 135 nM`,
/// true `Vmax = 120.0 µM/min`, `Km = 35.0 µM`, Competitive `Ki = 22.0 nM`).
pub fn demo_competitive_kinase_matrix_csv() -> &'static str {
    "\
substrate_um,inhibitor_nm,velocity
5.0,0.0,15.018
10.0,0.0,26.619
25.0,0.0,50.071
50.0,0.0,70.566
100.0,0.0,88.803
250.0,0.0,105.446
500.0,0.0,111.970
5.0,15.0,9.401
10.0,15.0,17.438
25.0,15.0,35.714
50.0,15.0,55.209
100.0,15.0,75.468
250.0,15.0,97.093
500.0,15.0,107.519
5.0,45.0,5.367
10.0,45.0,10.304
25.0,45.0,22.800
50.0,45.0,38.268
100.0,45.0,58.190
250.0,45.0,84.014
500.0,45.0,98.937
5.0,135.0,2.357
10.0,135.0,4.611
25.0,135.0,10.935
50.0,135.0,20.004
100.0,135.0,34.282
250.0,135.0,60.127
500.0,135.0,79.889
"
}

/// Built-in 12-point Morrison (1969) tight-binding titration dataset (`nominal [E] = 25 nM`,
/// true active `[E]_T = 20.0 nM` (`80%` active), `v0 = 95.0 RFU/s`, `[S] = 70 µM`, `Km = 35 µM`,
/// true `Ki = 0.45 nM` -> `Ki_app = 1.35 nM`).
pub fn demo_morrison_tight_binding_csv() -> &'static str {
    "\
inhibitor_nm,velocity
0.0,94.952
2.5,83.870
5.0,73.258
7.5,62.700
10.0,52.578
12.5,43.316
15.0,34.864
17.5,27.526
20.0,21.681
25.0,13.863
35.0,7.174
60.0,3.051
"
}

fn compute_aicc(n: usize, sse: f64, k_params: usize) -> f64 {
    let nf = n as f64;
    let kf = k_params as f64;
    let mse = (sse / nf.max(1.0)).max(1e-30);
    let aic = nf * mse.ln() + 2.0 * kf;
    if n > k_params + 1 {
        aic + (2.0 * kf * (kf + 1.0)) / ((n - k_params - 1) as f64)
    } else {
        aic + 100.0
    }
}

fn total_sum_of_squares(ys: &[f64]) -> f64 {
    if ys.is_empty() {
        return 0.0;
    }
    let m = ys.iter().sum::<f64>() / (ys.len() as f64);
    ys.iter().map(|&y| (y - m).powi(2)).sum()
}

fn median_positive(mut vals: Vec<f64>) -> f64 {
    vals.retain(|&x| x > 0.0 && x.is_finite());
    if vals.is_empty() {
        return 10.0;
    }
    vals.sort_by(|a, b| a.total_cmp(b));
    vals[vals.len() / 2]
}

/// Damped Levenberg–Marquardt optimizer in unconstrained log-parameter space (`param_k = exp(theta_k)`).
fn lm_optimize_log_params<F>(init_theta: &[f64], ys: &[f64], model_fn: F) -> Vec<f64>
where
    F: Fn(&[f64], usize) -> f64,
{
    let p_dim = init_theta.len();
    let n = ys.len();
    let mut theta = init_theta.to_vec();
    let mut lambda = 1e-3_f64;

    let calc_sse = |th: &[f64]| -> f64 {
        (0..n)
            .map(|i| {
                let r = ys[i] - model_fn(th, i);
                r * r
            })
            .sum()
    };

    let mut cur_sse = calc_sse(&theta);
    let eps = 1e-6_f64;

    for _iter in 0..140 {
        let mut jtj = vec![vec![0.0_f64; p_dim]; p_dim];
        let mut jtr = vec![0.0_f64; p_dim];

        for i in 0..n {
            let f0 = model_fn(&theta, i);
            let res = ys[i] - f0;
            let mut jac = vec![0.0_f64; p_dim];
            for k in 0..p_dim {
                let mut th_p = theta.clone();
                let mut th_m = theta.clone();
                th_p[k] += eps;
                th_m[k] -= eps;
                jac[k] = (model_fn(&th_p, i) - model_fn(&th_m, i)) / (2.0 * eps);
            }
            for r in 0..p_dim {
                jtr[r] += jac[r] * res;
                for c in 0..p_dim {
                    jtj[r][c] += jac[r] * jac[c];
                }
            }
        }

        let mut damped = jtj;
        for d in 0..p_dim {
            damped[d][d] += lambda * damped[d][d].max(1e-6);
        }

        let Some(step) = solve_dense_linear(damped, jtr) else {
            lambda = (lambda * 10.0).min(1e9);
            continue;
        };

        let cand: Vec<f64> = theta
            .iter()
            .zip(&step)
            .map(|(&t, &s)| (t + s.clamp(-2.5, 2.5)).clamp(-18.0, 18.0))
            .collect();
        let cand_sse = calc_sse(&cand);

        if cand_sse < cur_sse {
            let rel = (cur_sse - cand_sse) / cur_sse.max(1e-15);
            theta = cand;
            cur_sse = cand_sse;
            lambda = (lambda * 0.35).max(1e-11);
            if rel < 1e-12 {
                break;
            }
        } else {
            lambda = (lambda * 4.0).min(1e9);
        }
    }

    theta
}

fn solve_dense_linear(mut a: Vec<Vec<f64>>, mut b: Vec<f64>) -> Option<Vec<f64>> {
    let m = b.len();
    for col in 0..m {
        let mut piv = col;
        for r in (col + 1)..m {
            if a[r][col].abs() > a[piv][col].abs() {
                piv = r;
            }
        }
        if a[piv][col].abs() < 1e-18 {
            return None;
        }
        a.swap(col, piv);
        b.swap(col, piv);

        let diag = a[col][col];
        for r in (col + 1)..m {
            let f = a[r][col] / diag;
            b[r] -= f * b[col];
            for c in col..m {
                a[r][c] -= f * a[col][c];
            }
        }
    }

    let mut x = vec![0.0_f64; m];
    for i in (0..m).rev() {
        let mut s = b[i];
        for j in (i + 1)..m {
            s -= a[i][j] * x[j];
        }
        x[i] = s / a[i][i];
    }
    Some(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn global_aicc_identifies_competitive_inhibition_and_recovers_ki() {
        let pts = parse_kinetics_csv(demo_competitive_kinase_matrix_csv()).unwrap();
        let fits = fit_global_reversible_inhibition(&pts).unwrap();
        assert_eq!(fits[0].mechanism, InhibitionMechanism::Competitive);
        assert!(fits[0].akaike_weight > 0.65);
        assert!((fits[0].vmax - 120.0).abs() < 2.5);
        assert!((fits[0].km_um - 35.0).abs() < 2.5);
        assert!((fits[0].ki_nm - 22.0).abs() < 1.5);
    }

    #[test]
    fn morrison_quadratic_recovers_active_enzyme_concentration_and_intrinsic_ki() {
        let raw_pts = parse_kinetics_csv(demo_morrison_tight_binding_csv()).unwrap();
        let obs: Vec<(f64, f64)> = raw_pts.iter().map(|p| (p.substrate_um, p.velocity)).collect();
        let fit =
            fit_morrison_tight_binding(&obs, None, Some(25.0), Some((70.0, 35.0))).unwrap();
        assert!((fit.active_enzyme_et_nm - 20.0).abs() < 1.2);
        assert!((fit.ki_app_nm - 1.35).abs() < 0.25);
        assert!((fit.intrinsic_ki_nm.unwrap() - 0.45).abs() < 0.10);
        assert!(fit.r_squared > 0.999);
    }

    #[test]
    fn single_substrate_fit_filters_nan_and_inf_in_fallback_branch() {
        // Fewer than 5 uninhibited points forces fallback to `points` slice,
        // where NaN and Inf velocities must be filtered rather than panicking in sorting.
        let pts = vec![
            KineticsPoint { substrate_um: 5.0, inhibitor_nm: 15.0, velocity: 9.4 },
            KineticsPoint { substrate_um: 10.0, inhibitor_nm: 15.0, velocity: 17.4 },
            KineticsPoint { substrate_um: 25.0, inhibitor_nm: 15.0, velocity: 35.7 },
            KineticsPoint { substrate_um: 50.0, inhibitor_nm: 15.0, velocity: 55.2 },
            KineticsPoint { substrate_um: 100.0, inhibitor_nm: 15.0, velocity: 75.5 },
            KineticsPoint { substrate_um: 250.0, inhibitor_nm: 15.0, velocity: 97.1 },
            KineticsPoint { substrate_um: 75.0, inhibitor_nm: 15.0, velocity: f64::NAN },
            KineticsPoint { substrate_um: 150.0, inhibitor_nm: 15.0, velocity: f64::INFINITY },
        ];
        let fits = fit_single_substrate_models(&pts, None).expect("NaN/Inf filtered safely");
        assert_eq!(fits.len(), 3);
        assert!(fits.iter().all(|f| f.vmax.is_finite() && f.aicc.is_finite()));

        let too_few = vec![
            KineticsPoint { substrate_um: 5.0, inhibitor_nm: 15.0, velocity: 9.4 },
            KineticsPoint { substrate_um: 10.0, inhibitor_nm: 15.0, velocity: f64::NAN },
            KineticsPoint { substrate_um: 25.0, inhibitor_nm: 15.0, velocity: f64::NAN },
            KineticsPoint { substrate_um: 50.0, inhibitor_nm: 15.0, velocity: 55.2 },
            KineticsPoint { substrate_um: 100.0, inhibitor_nm: 15.0, velocity: f64::INFINITY },
        ];
        assert!(fit_single_substrate_models(&too_few, None).is_err());
    }
}
