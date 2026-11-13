// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #5 (`osabt epi`) — Epidemiological mathematical & statistical models:
//! 1. High-precision pure-Rust special functions (`ln_gamma`, regularized incomplete Gamma
//!    `gamma_p`/`gamma_q`, Newton-refined Gamma quantiles, normal CDF `norm_cdf`, and Acklam's
//!    inverse normal quantile `norm_ppf`).
//! 2. Cori et al. (2013) `EpiEstim` Bayesian Gamma–Poisson renewal estimation of time-varying
//!    reproduction numbers (`R_t`) with discretized serial interval distributions.
//! 3. 4th-order Runge–Kutta (`RK4`) deterministic `SIR` and `SEIR` compartmental transmission
//!    simulators with vaccination/pre-immunity and NPI transmission-reduction interventions.
//! 4. Sequential multi-assay Bayesian diagnostic updater (`LR+`, `LR-`, `PPV`, `NPV`) and
//!    Rogan & Gladen (1978) misclassification-adjusted true prevalence estimator.

use std::f64::consts::{PI, SQRT_2};

// ============================================================================
// 1. Pure-Rust Statistical Special Functions
// ============================================================================

/// Lanczos `g = 7` natural log-Gamma function `ln(Gamma(z))` for `z > 0` (~15 digits accuracy).
pub fn ln_gamma(z: f64) -> f64 {
    const COEFFS: [f64; 9] = [
        0.999_999_999_999_809_9,
        676.520_368_121_885_1,
        -1259.139_216_722_402_8,
        771.323_428_777_653_1,
        -176.615_029_162_140_6,
        12.507_343_278_686_905,
        -0.138_571_095_265_720_12,
        9.984_369_578_019_572e-6,
        1.505_632_735_149_311_6e-7,
    ];
    if z < 0.5 {
        return PI.ln() - (PI * z).sin().ln() - ln_gamma(1.0 - z);
    }
    let z_m1 = z - 1.0;
    let mut x = COEFFS[0];
    for (i, &c) in COEFFS.iter().enumerate().skip(1) {
        x += c / (z_m1 + (i as f64));
    }
    let t = z_m1 + 7.5;
    0.5 * (2.0 * PI).ln() + (z_m1 + 0.5) * t.ln() - t + x.ln()
}

/// Lower regularized incomplete Gamma function `P(a, x) = gamma(a, x) / Gamma(a)` for `a > 0, x >= 0`.
pub fn gamma_p(a: f64, x: f64) -> f64 {
    if x <= 0.0 || a <= 0.0 {
        return 0.0;
    }
    if x < a + 1.0 {
        gamma_series(a, x)
    } else {
        1.0 - gamma_continued_fraction(a, x)
    }
}

/// Upper regularized incomplete Gamma function `Q(a, x) = 1 - P(a, x)`.
pub fn gamma_q(a: f64, x: f64) -> f64 {
    if x <= 0.0 {
        return 1.0;
    }
    if x < a + 1.0 {
        1.0 - gamma_series(a, x)
    } else {
        gamma_continued_fraction(a, x)
    }
}

fn gamma_series(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let mut ap = a;
    let mut sum = 1.0 / a;
    let mut del = sum;
    for _ in 0..200 {
        ap += 1.0;
        del *= x / ap;
        sum += del;
        if del.abs() < sum.abs() * 1e-14 {
            break;
        }
    }
    (sum * (-x + a * x.ln() - gln).exp()).clamp(0.0, 1.0)
}

fn gamma_continued_fraction(a: f64, x: f64) -> f64 {
    let gln = ln_gamma(a);
    let fpmin = 1e-300_f64;
    let mut b = x + 1.0 - a;
    let mut c = 1.0 / fpmin;
    let mut d = 1.0 / b.max(fpmin);
    let mut h = d;
    for i in 1..=200 {
        let fi = i as f64;
        let an = -fi * (fi - a);
        b += 2.0;
        d = an * d + b;
        if d.abs() < fpmin {
            d = fpmin;
        }
        c = b + an / c;
        if c.abs() < fpmin {
            c = fpmin;
        }
        d = 1.0 / d;
        let del = d * c;
        h *= del;
        if (del - 1.0).abs() < 1e-14 {
            break;
        }
    }
    ((-x + a * x.ln() - gln).exp() * h).clamp(0.0, 1.0)
}

/// Inverse CDF (quantile) of `Gamma(shape = a, rate = b)` at probability `p in (0, 1)`.
/// Uses Wilson–Hilferty cube-root initialization followed by safeguarded Newton–Raphson steps.
pub fn gamma_quantile_shape_rate(shape_a: f64, rate_b: f64, p: f64) -> f64 {
    if p <= 0.0 || shape_a <= 0.0 || rate_b <= 0.0 {
        return 0.0;
    }
    if p >= 1.0 {
        return f64::INFINITY;
    }

    // Solve for unit-rate `x` such that `gamma_p(shape_a, x) == p`, then return `x / rate_b`
    let mut x = if shape_a >= 1.0 {
        let z = norm_ppf(p);
        let wh = 1.0 - 1.0 / (9.0 * shape_a) + z / (3.0 * shape_a.sqrt());
        (shape_a * wh.max(0.05).powi(3)).max(1e-9)
    } else {
        (p.powf(1.0 / shape_a) * (ln_gamma(shape_a + 1.0)).exp()).max(1e-12)
    };

    let gln = ln_gamma(shape_a);
    let mut lo = 0.0_f64;
    let mut hi = (shape_a + 40.0 * shape_a.sqrt().max(1.0)).max(x * 4.0);
    while gamma_p(shape_a, hi) < p {
        hi *= 2.0;
    }

    for _ in 0..24 {
        let cdf = gamma_p(shape_a, x);
        let err = cdf - p;
        if err.abs() < 1e-12 {
            break;
        }
        if err < 0.0 {
            lo = lo.max(x);
        } else {
            hi = hi.min(x);
        }
        let log_pdf = (shape_a - 1.0) * x.ln() - x - gln;
        let pdf = log_pdf.exp();
        let next = if pdf > 1e-30 { x - err / pdf } else { 0.5 * (lo + hi) };
        x = if next > lo && next < hi {
            next
        } else {
            0.5 * (lo + hi)
        };
    }

    x / rate_b
}

/// Standard normal cumulative distribution function `Phi(z)`.
pub fn norm_cdf(z: f64) -> f64 {
    0.5 * (1.0 + erf_approx(z / SQRT_2))
}

fn erf_approx(x: f64) -> f64 {
    // Abramowitz & Stegun 7.1.26 / rational approximation (max abs error <= 1.5e-7; enhanced with 9-term Chebyshev erfc)
    let z = x.abs();
    let t = 1.0 / (1.0 + 0.5 * z);
    let ans = t
        * (-z * z - 1.265_512_23
            + t * (1.000_023_68
                + t * (0.374_091_96
                    + t * (0.096_784_18
                        + t * (-0.186_288_06
                            + t * (0.278_868_07
                                + t * (-1.135_203_98
                                    + t * (1.488_515_87
                                        + t * (-0.822_152_23 + t * 0.170_872_77)))))))))
            .exp();
    if x >= 0.0 {
        1.0 - ans
    } else {
        ans - 1.0
    }
}

/// Peter J. Acklam's rational approximation to the inverse standard normal CDF `Phi^-1(p)`
/// (maximum relative error `< 1.15e-9` across `p in (0, 1)`).
pub fn norm_ppf(p: f64) -> f64 {
    const A: [f64; 6] = [
        -3.969_683_028_665_376e1,
        2.209_460_984_245_205e2,
        -2.759_285_104_469_687e2,
        1.383_577_518_672_69e2,
        -3.066_479_806_614_716e1,
        2.506_628_277_459_239,
    ];
    const B: [f64; 5] = [
        -5.447_609_879_822_406e1,
        1.615_858_368_580_409e2,
        -1.556_989_798_598_866e2,
        6.680_131_188_771_972e1,
        -1.328_068_155_288_572e1,
    ];
    const C: [f64; 6] = [
        -7.784_894_002_430_293e-3,
        -3.223_964_580_411_365e-1,
        -2.400_758_277_161_838,
        -2.549_732_539_343_734,
        4.374_664_141_464_968,
        2.938_163_982_698_783,
    ];
    const D: [f64; 4] = [
        7.784_695_709_041_462e-3,
        3.224_671_290_700_398e-1,
        2.445_134_137_142_996,
        3.754_408_661_907_416,
    ];

    let p_clamped = p.clamp(1e-15, 1.0 - 1e-15);
    let p_low = 0.02425;
    let p_high = 1.0 - p_low;

    if p_clamped < p_low {
        let q = (-2.0 * p_clamped.ln()).sqrt();
        (((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    } else if p_clamped <= p_high {
        let q = p_clamped - 0.5;
        let r = q * q;
        (((((A[0] * r + A[1]) * r + A[2]) * r + A[3]) * r + A[4]) * r + A[5]) * q
            / (((((B[0] * r + B[1]) * r + B[2]) * r + B[3]) * r + B[4]) * r + 1.0)
    } else {
        let q = (-2.0 * (1.0 - p_clamped).ln()).sqrt();
        -(((((C[0] * q + C[1]) * q + C[2]) * q + C[3]) * q + C[4]) * q + C[5])
            / ((((D[0] * q + D[1]) * q + D[2]) * q + D[3]) * q + 1.0)
    }
}

// ============================================================================
// 2. Cori et al. (2013) EpiEstim Renewal R_t Estimation
// ============================================================================

/// Parameters for Cori et al. (2013) instantaneous reproduction number (`R_t`) estimation.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct CoriRtConfig {
    /// Mean of the positive serial interval distribution (days, must be `> 1.0`).
    pub serial_interval_mean: f64,
    /// Standard deviation of the serial interval distribution (days, `> 0.0`).
    pub serial_interval_sd: f64,
    /// Length `tau` of the trailing estimation window in days (standard `EpiEstim` default = `7`).
    pub window_days: usize,
    /// Gamma prior shape `a_0` (default `1.0` for mean `5.0`, SD `5.0`).
    pub prior_shape: f64,
    /// Gamma prior rate `b_0` (default `0.2` for mean `5.0`, SD `5.0`).
    pub prior_rate: f64,
}

impl Default for CoriRtConfig {
    fn default() -> Self {
        Self {
            serial_interval_mean: 4.7,
            serial_interval_sd: 2.9,
            window_days: 7,
            prior_shape: 1.0,
            prior_rate: 0.2,
        }
    }
}

/// Time-step estimate of `R_t` produced by [`estimate_cori_rt`].
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct RtEstimatePoint {
    pub day_index: usize,
    pub incidence: f64,
    pub total_infectiousness_lambda: f64,
    pub rt_mean: Option<f64>,
    pub rt_median: Option<f64>,
    pub ci95_low: Option<f64>,
    pub ci95_high: Option<f64>,
    /// Posterior probability `Pr(R_t > 1.0 | I_0..t)` indicating active epidemic growth.
    pub prob_rt_gt_1: Option<f64>,
}

/// Discretizes a shifted Gamma serial interval (`mean > 0, sd > 0`) over integer lags `u = 1..=max_lag`
/// so `w_0 = 0` and `sum_{u=1..=max_lag} w_u = 1.0`.
pub fn discretize_serial_interval(mean: f64, sd: f64, max_lag: usize) -> Vec<f64> {
    let mu = mean.max(0.5);
    let sigma = sd.max(0.1);
    let shape = (mu * mu) / (sigma * sigma);
    let scale = (sigma * sigma) / mu;

    let n = max_lag.max(1);
    let mut w = vec![0.0_f64; n + 1];
    let mut sum = 0.0_f64;
    for (u, slot) in w.iter_mut().enumerate().skip(1) {
        let hi = ((u as f64) + 0.5) / scale;
        let lo = (((u as f64) - 0.5).max(0.0)) / scale;
        let mass = (gamma_p(shape, hi) - gamma_p(shape, lo)).max(0.0);
        *slot = mass;
        sum += mass;
    }
    if sum > 0.0 {
        for slot in w.iter_mut().skip(1) {
            *slot /= sum;
        }
    }
    w
}

/// Estimates time-varying `R_t` from daily incidence counts using Cori et al. (2013)
/// analytical Gamma-Poisson conjugate renewal posteriors.
pub fn estimate_cori_rt(
    daily_incidence: &[f64],
    cfg: &CoriRtConfig,
) -> Result<Vec<RtEstimatePoint>, String> {
    if daily_incidence.len() < 3 {
        return Err(String::from(
            "Incidence time series must contain at least 3 daily observations.",
        ));
    }
    if cfg.serial_interval_mean <= 0.5 || cfg.serial_interval_sd <= 0.0 {
        return Err(String::from(
            "Serial interval mean (> 0.5) and SD (> 0.0) must be positive.",
        ));
    }

    let n = daily_incidence.len();
    let tau = cfg.window_days.max(1).min(n);
    let w = discretize_serial_interval(cfg.serial_interval_mean, cfg.serial_interval_sd, n);

    let mut lambdas = vec![0.0_f64; n];
    for t in 1..n {
        let mut lam = 0.0_f64;
        for u in 1..=t {
            lam += daily_incidence[t - u].max(0.0) * w[u];
        }
        lambdas[t] = lam;
    }

    let mut out = Vec::with_capacity(n);
    for t in 0..n {
        let inc = daily_incidence[t].max(0.0);
        let lam = lambdas[t];
        if t + 1 < tau || t == 0 {
            out.push(RtEstimatePoint {
                day_index: t,
                incidence: inc,
                total_infectiousness_lambda: lam,
                rt_mean: None,
                rt_median: None,
                ci95_low: None,
                ci95_high: None,
                prob_rt_gt_1: None,
            });
            continue;
        }

        let start = t + 1 - tau;
        let sum_i: f64 = daily_incidence[start..=t].iter().map(|v| v.max(0.0)).sum();
        let sum_lam: f64 = lambdas[start..=t].iter().sum();

        if sum_lam <= 1e-9 {
            out.push(RtEstimatePoint {
                day_index: t,
                incidence: inc,
                total_infectiousness_lambda: lam,
                rt_mean: None,
                rt_median: None,
                ci95_low: None,
                ci95_high: None,
                prob_rt_gt_1: None,
            });
            continue;
        }

        let a_post = cfg.prior_shape + sum_i;
        let b_post = cfg.prior_rate + sum_lam;
        let mean = a_post / b_post;
        let median = gamma_quantile_shape_rate(a_post, b_post, 0.50);
        let lo = gamma_quantile_shape_rate(a_post, b_post, 0.025);
        let hi = gamma_quantile_shape_rate(a_post, b_post, 0.975);
        // Pr(R_t > 1) = Q(a_post, b_post * 1.0)
        let prob_gt_1 = gamma_q(a_post, b_post);

        out.push(RtEstimatePoint {
            day_index: t,
            incidence: inc,
            total_infectiousness_lambda: lam,
            rt_mean: Some(mean),
            rt_median: Some(median),
            ci95_low: Some(lo),
            ci95_high: Some(hi),
            prob_rt_gt_1: Some(prob_gt_1),
        });
    }

    Ok(out)
}

/// Built-in 35-day respiratory outbreak daily case series illustrating growth (`R_t ~ 2.2`),
/// intervention inflection around day 15, and post-control decline (`R_t ~ 0.65`).
pub fn demo_outbreak_incidence() -> Vec<f64> {
    vec![
        2.0, 3.0, 4.0, 6.0, 9.0, 13.0, 18.0, 26.0, 35.0, 48.0, 64.0, 82.0, 105.0, 128.0, 146.0,
        152.0, 148.0, 136.0, 121.0, 106.0, 91.0, 78.0, 66.0, 55.0, 46.0, 38.0, 31.0, 25.0, 20.0,
        16.0, 13.0, 10.0, 8.0, 6.0, 5.0,
    ]
}

// ============================================================================
// 3. Deterministic RK4 SIR / SEIR Compartmental Transmission Simulator
// ============================================================================

/// Configuration for an RK4 `SIR` (`incubation_days <= 0`) or `SEIR` (`incubation_days > 0`) run.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SeirParams {
    pub population: f64,
    pub initial_infectious: f64,
    pub initial_exposed: f64,
    /// Baseline immune / vaccinated proportion in `0.0..1.0`.
    pub initial_immune_fraction: f64,
    /// Basic reproduction number `R_0` in a wholly susceptible population.
    pub r0: f64,
    /// Mean latent / incubation period `1 / sigma` in days (`0.0` reduces cleanly to classic `SIR`).
    pub incubation_days: f64,
    /// Mean infectious duration `1 / gamma` in days.
    pub infectious_days: f64,
    /// Total simulation horizon in days.
    pub total_days: usize,
    /// Optional non-pharmaceutical intervention (NPI) start day (`None` = unmitigated).
    pub npi_start_day: Option<f64>,
    /// Fractional reduction in transmission `beta` (`0.0..=1.0`) once `t >= npi_start_day`.
    pub npi_reduction_fraction: f64,
}

impl Default for SeirParams {
    fn default() -> Self {
        Self {
            population: 100_000.0,
            initial_infectious: 10.0,
            initial_exposed: 15.0,
            initial_immune_fraction: 0.15,
            r0: 2.8,
            incubation_days: 3.5,
            infectious_days: 5.0,
            total_days: 160,
            npi_start_day: Some(35.0),
            npi_reduction_fraction: 0.40,
        }
    }
}

/// Daily state snapshot of the compartmental trajectory.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SeirTimePoint {
    pub day: f64,
    pub susceptible: f64,
    pub exposed: f64,
    pub infectious: f64,
    pub recovered: f64,
    pub daily_new_infections: f64,
    pub effective_re: f64,
}

/// Complete trajectory + epidemic summary metrics returned by [`simulate_seir_rk4`].
#[derive(Debug, Clone)]
pub struct SeirSimulationResult {
    pub params: SeirParams,
    pub herd_immunity_threshold: f64,
    pub initial_effective_re: f64,
    pub peak_infectious: f64,
    pub peak_day: f64,
    pub final_attack_fraction: f64,
    pub total_infected_persons: f64,
    pub trajectory: Vec<SeirTimePoint>,
}

/// Integrates the deterministic `SIR` or `SEIR` ODE system via 4th-order Runge–Kutta (`dt = 0.1` day).
pub fn simulate_seir_rk4(params: &SeirParams) -> Result<SeirSimulationResult, String> {
    if params.population <= 0.0 || params.infectious_days <= 0.0 || params.r0 < 0.0 {
        return Err(String::from(
            "Population (>0), infectious_days (>0), and R0 (>=0) must be valid.",
        ));
    }

    let n_pop = params.population;
    let gamma = 1.0 / params.infectious_days;
    let beta0 = params.r0 * gamma;
    let is_seir = params.incubation_days > 1e-6;
    let sigma = if is_seir {
        1.0 / params.incubation_days
    } else {
        0.0
    };

    let init_r = (n_pop * params.initial_immune_fraction.clamp(0.0, 0.999)).min(n_pop);
    let init_i = params.initial_infectious.max(0.0).min(n_pop - init_r);
    let init_e = if is_seir {
        params
            .initial_exposed
            .max(0.0)
            .min(n_pop - init_r - init_i)
    } else {
        0.0
    };
    let init_s = (n_pop - init_r - init_i - init_e).max(0.0);

    let substeps = 10usize;
    let dt = 1.0 / (substeps as f64);
    let mut state = [init_s, init_e, init_i, init_r];
    let mut trajectory = Vec::with_capacity(params.total_days + 1);

    let beta_at = |t: f64| -> f64 {
        match params.npi_start_day {
            Some(t_npi) if t >= t_npi => {
                beta0 * (1.0 - params.npi_reduction_fraction.clamp(0.0, 1.0))
            }
            _ => beta0,
        }
    };

    let deriv = |t: f64, y: &[f64; 4]| -> [f64; 4] {
        let (s, e, i, _r) = (y[0], y[1], y[2], y[3]);
        let b = beta_at(t);
        let foi = b * s * i / n_pop;
        if is_seir {
            [-foi, foi - sigma * e, sigma * e - gamma * i, gamma * i]
        } else {
            [-foi, 0.0, foi - gamma * i, gamma * i]
        }
    };

    let mut peak_i = init_i;
    let mut peak_day = 0.0_f64;

    for day in 0..=params.total_days {
        let t_now = day as f64;
        let re_now = (beta_at(t_now) / gamma) * (state[0] / n_pop);
        if state[2] > peak_i {
            peak_i = state[2];
            peak_day = t_now;
        }

        let s_before = state[0];
        if day < params.total_days {
            for step in 0..substeps {
                let t = t_now + (step as f64) * dt;
                let k1 = deriv(t, &state);
                let y2 = [
                    state[0] + 0.5 * dt * k1[0],
                    state[1] + 0.5 * dt * k1[1],
                    state[2] + 0.5 * dt * k1[2],
                    state[3] + 0.5 * dt * k1[3],
                ];
                let k2 = deriv(t + 0.5 * dt, &y2);
                let y3 = [
                    state[0] + 0.5 * dt * k2[0],
                    state[1] + 0.5 * dt * k2[1],
                    state[2] + 0.5 * dt * k2[2],
                    state[3] + 0.5 * dt * k2[3],
                ];
                let k3 = deriv(t + 0.5 * dt, &y3);
                let y4 = [
                    state[0] + dt * k3[0],
                    state[1] + dt * k3[1],
                    state[2] + dt * k3[2],
                    state[3] + dt * k3[3],
                ];
                let k4 = deriv(t + dt, &y4);
                for c in 0..4 {
                    state[c] = (state[c]
                        + (dt / 6.0) * (k1[c] + 2.0 * k2[c] + 2.0 * k3[c] + k4[c]))
                        .max(0.0);
                }
            }
        }
        let new_inf = (s_before - state[0]).max(0.0);

        trajectory.push(SeirTimePoint {
            day: t_now,
            susceptible: if day == params.total_days {
                state[0]
            } else {
                s_before
            },
            exposed: state[1],
            infectious: state[2],
            recovered: state[3],
            daily_new_infections: new_inf,
            effective_re: re_now,
        });
    }

    let final_s = state[0];
    let total_newly_infected = (init_s - final_s).max(0.0) + init_e + init_i;
    let hit = if params.r0 > 1.0 {
        1.0 - 1.0 / params.r0
    } else {
        0.0
    };

    Ok(SeirSimulationResult {
        params: params.clone(),
        herd_immunity_threshold: hit,
        initial_effective_re: params.r0 * (init_s / n_pop),
        peak_infectious: peak_i,
        peak_day,
        final_attack_fraction: total_newly_infected / n_pop,
        total_infected_persons: total_newly_infected,
        trajectory,
    })
}

// ============================================================================
// 4. Diagnostic Bayesian Cascade & Rogan-Gladen (1978) Prevalence
// ============================================================================

/// Specification of a single diagnostic test step in a sequential testing cascade.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct DiagnosticStepInput {
    pub assay_name: String,
    pub sensitivity: f64,
    pub specificity: f64,
    pub observed_positive: bool,
}

/// Output of a single sequential Bayesian update step.
#[derive(Debug, Clone)]
pub struct DiagnosticStepResult {
    pub assay_name: String,
    pub observed_positive: bool,
    pub lr_plus: f64,
    pub lr_minus: f64,
    pub applied_lr: f64,
    pub youden_j: f64,
    pub pre_test_prob: f64,
    pub post_test_prob: f64,
    /// Single-step PPV (`Pr(D+ | T+)`) at this stage's pre-test probability.
    pub stage_ppv: f64,
    /// Single-step NPV (`Pr(D- | T-)`) at this stage's pre-test probability.
    pub stage_npv: f64,
}

/// Evaluates a 1-or-more step sequential Bayesian diagnostic cascade starting from `prior_prevalence`.
pub fn evaluate_diagnostic_cascade(
    prior_prevalence: f64,
    steps: &[DiagnosticStepInput],
) -> Result<Vec<DiagnosticStepResult>, String> {
    if !(0.0..=1.0).contains(&prior_prevalence) {
        return Err(String::from("Prior prevalence must lie in 0.0..=1.0."));
    }
    let mut cur_p = prior_prevalence.clamp(1e-12, 1.0 - 1e-12);
    let mut results = Vec::with_capacity(steps.len());

    for st in steps {
        if !(0.0..=1.0).contains(&st.sensitivity) || !(0.0..=1.0).contains(&st.specificity) {
            return Err(format!(
                "Assay '{}' sensitivity and specificity must lie in 0.0..=1.0.",
                st.assay_name
            ));
        }
        let se = st.sensitivity.clamp(1e-9, 1.0 - 1e-9);
        let sp = st.specificity.clamp(1e-9, 1.0 - 1e-9);
        let lr_plus = se / (1.0 - sp);
        let lr_minus = (1.0 - se) / sp;
        let applied_lr = if st.observed_positive {
            lr_plus
        } else {
            lr_minus
        };

        let pre_odds = cur_p / (1.0 - cur_p);
        let ppv = (se * cur_p) / (se * cur_p + (1.0 - sp) * (1.0 - cur_p));
        let npv = (sp * (1.0 - cur_p)) / ((1.0 - se) * cur_p + sp * (1.0 - cur_p));

        let post_odds = pre_odds * applied_lr;
        let post_p = post_odds / (1.0 + post_odds);

        results.push(DiagnosticStepResult {
            assay_name: st.assay_name.clone(),
            observed_positive: st.observed_positive,
            lr_plus,
            lr_minus,
            applied_lr,
            youden_j: st.sensitivity + st.specificity - 1.0,
            pre_test_prob: cur_p,
            post_test_prob: post_p,
            stage_ppv: ppv,
            stage_npv: npv,
        });
        cur_p = post_p;
    }

    Ok(results)
}

/// Result of Rogan & Gladen (1978) imperfect-assay true prevalence estimation.
#[derive(Debug, Clone)]
pub struct RoganGladenResult {
    pub positive_count: usize,
    pub sample_size: usize,
    pub apparent_prevalence: f64,
    pub apparent_ci95: (f64, f64),
    pub true_prevalence: f64,
    pub true_ci95_clamped: (f64, f64),
    pub youden_j: f64,
}

/// Estimates true disease prevalence from `positives` out of `total_tested` screened subjects
/// using Rogan & Gladen (1978) with Wilson score confidence limits mapped through `(p + Sp - 1) / J`.
pub fn rogan_gladen_prevalence(
    positives: usize,
    total_tested: usize,
    sensitivity: f64,
    specificity: f64,
) -> Result<RoganGladenResult, String> {
    if total_tested == 0 || positives > total_tested {
        return Err(String::from("Require total_tested > 0 and positives <= total_tested."));
    }
    let j = sensitivity + specificity - 1.0;
    if j <= 1e-6 {
        return Err(String::from(
            "Rogan-Gladen estimator requires informative assay with Se + Sp > 1.0.",
        ));
    }

    let n = total_tested as f64;
    let p_app = (positives as f64) / n;
    let z = 1.959_963_984_540_054; // Phi^-1(0.975)
    let z2 = z * z;
    let denom = 1.0 + z2 / n;
    let center = (p_app + z2 / (2.0 * n)) / denom;
    let half = (z / denom) * ((p_app * (1.0 - p_app) / n) + (z2 / (4.0 * n * n))).sqrt();
    let app_lo = (center - half).clamp(0.0, 1.0);
    let app_hi = (center + half).clamp(0.0, 1.0);

    let map_rg = |pa: f64| ((pa + specificity - 1.0) / j).clamp(0.0, 1.0);
    Ok(RoganGladenResult {
        positive_count: positives,
        sample_size: total_tested,
        apparent_prevalence: p_app,
        apparent_ci95: (app_lo, app_hi),
        true_prevalence: map_rg(p_app),
        true_ci95_clamped: (map_rg(app_lo), map_rg(app_hi)),
        youden_j: j,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn verifies_gamma_quantiles_cori_rt_seir_and_rogan_gladen() {
        // Gamma(shape=10, rate=2) has mean 5.0; median ~ 4.834357
        let q50 = gamma_quantile_shape_rate(10.0, 2.0, 0.5);
        assert!((q50 - 4.834_357).abs() < 1e-4);

        let rt = estimate_cori_rt(&demo_outbreak_incidence(), &CoriRtConfig::default()).unwrap();
        // Early day 10 should have R_t > 1.5; late day 30 should have R_t < 1.0
        assert!(rt[10].rt_mean.unwrap() > 1.5);
        assert!(rt[30].rt_mean.unwrap() < 1.0);

        let seir = simulate_seir_rk4(&SeirParams::default()).unwrap();
        assert!(seir.peak_infectious > 500.0);
        assert!((seir.herd_immunity_threshold - (1.0 - 1.0 / 2.8)).abs() < 1e-9);

        // Rogan-Gladen: 140/1000 apparent (14%), Se=0.90, Sp=0.95 -> (0.14 + 0.95 - 1)/0.85 = 0.09/0.85 = 0.105882
        let rg = rogan_gladen_prevalence(140, 1000, 0.90, 0.95).unwrap();
        assert!((rg.true_prevalence - (0.09 / 0.85)).abs() < 1e-6);
    }
}
