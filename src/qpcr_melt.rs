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

//! Tool #12 (`osabt qpcr`) — SYBR Green I / EvaGreen DNA Melting-Curve (`-dF/dT`)
//! Savitzky–Golay derivative peak caller & primer-dimer specificity auditor, grounded in
//! Ririe, Rasmussen & Wittwer (1997) *Anal. Biochem.* 245(2):154–160.

use std::collections::BTreeMap;
use std::fmt::Write as _;

use crate::svg_util::xml_escape;

/// Raw thermal dissociation fluorescence trace (`60 °C .. 95 °C`) for one qPCR well.
#[derive(Debug, Clone)]
pub struct MeltTrace {
    pub well_id: String,
    pub sample_id: String,
    pub target_gene: String,
    pub temperatures_c: Vec<f64>,
    pub fluorescence_rfu: Vec<f64>,
}

/// Detected local `-dF/dT` melting peak with parabolic sub-step `Tm` and FWHM width.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MeltPeakCall {
    pub tm_celsius: f64,
    pub peak_neg_df_dt: f64,
    pub fwhm_celsius: f64,
    pub relative_height_pct: f64,
    pub is_primer_dimer_suspect: bool,
}

/// Processed derivative melt curve + amplicon specificity classification for one well.
#[derive(Debug, Clone)]
pub struct MeltCurveResult {
    pub well_id: String,
    pub sample_id: String,
    pub target_gene: String,
    pub temperatures_c: Vec<f64>,
    pub neg_df_dt: Vec<f64>,
    pub primary_tm_c: Option<f64>,
    pub peaks: Vec<MeltPeakCall>,
    pub specificity_status: &'static str,
}

/// Parses a long-format melt CSV/TSV (`well,sample,target,temperature_c,fluorescence`).
pub fn parse_melt_csv(raw: &str) -> Result<Vec<MeltTrace>, String> {
    let mut grouped: BTreeMap<(String, String, String), Vec<(f64, f64)>> = BTreeMap::new();
    let mut header_seen = false;

    for (idx, line) in raw.trim_start_matches('\u{feff}').lines().enumerate() {
        let t = line.trim();
        if t.is_empty() || t.starts_with('#') {
            continue;
        }
        let sep = if t.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = t.split(sep).map(|s| s.trim()).collect();
        if !header_seen {
            header_seen = true;
            if cols[0].eq_ignore_ascii_case("well") || cols[0].eq_ignore_ascii_case("well_id") {
                continue;
            }
        }
        if cols.len() < 5 {
            return Err(format!(
                "Melt CSV line {}: expected 5 columns (well,sample,target,temperature_c,fluorescence).",
                idx + 1
            ));
        }
        let temp: f64 = cols[3]
            .parse()
            .map_err(|_| format!("Invalid temperature '{}' on line {}.", cols[3], idx + 1))?;
        let fl: f64 = cols[4]
            .parse()
            .map_err(|_| format!("Invalid fluorescence '{}' on line {}.", cols[4], idx + 1))?;
        grouped
            .entry((
                cols[0].to_string(),
                cols[1].to_string(),
                cols[2].to_string(),
            ))
            .or_default()
            .push((temp, fl));
    }

    let mut traces = Vec::with_capacity(grouped.len());
    for ((well_id, sample_id, target_gene), mut pts) in grouped {
        pts.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
        traces.push(MeltTrace {
            well_id,
            sample_id,
            target_gene,
            temperatures_c: pts.iter().map(|p| p.0).collect(),
            fluorescence_rfu: pts.iter().map(|p| p.1).collect(),
        });
    }
    if traces.is_empty() {
        return Err(String::from("No melt-curve traces found."));
    }
    Ok(traces)
}

/// Computes Savitzky–Golay smoothed `-dF/dT` melting derivative, calls primary & secondary peaks,
/// and flags primer-dimer (`Tm < dimer_cutoff_c`, default `77.5 °C`) or multimodal artifacts.
pub fn analyze_melt_curve(trace: &MeltTrace, dimer_cutoff_c: f64) -> MeltCurveResult {
    let n = trace.temperatures_c.len().min(trace.fluorescence_rfu.len());
    if n < 7 {
        return MeltCurveResult {
            well_id: trace.well_id.clone(),
            sample_id: trace.sample_id.clone(),
            target_gene: trace.target_gene.clone(),
            temperatures_c: trace.temperatures_c.clone(),
            neg_df_dt: vec![0.0; n],
            primary_tm_c: None,
            peaks: Vec::new(),
            specificity_status: "INSUFFICIENT POINTS",
        };
    }

    // 5-point Savitzky–Golay quadratic smoothing (`[-3, 12, 17, 12, -3] / 35`) on raw fluorescence
    let f_raw = &trace.fluorescence_rfu[..n];
    let mut f_smooth = f_raw.to_vec();
    for i in 2..(n - 2) {
        f_smooth[i] = (-3.0 * f_raw[i - 2]
            + 12.0 * f_raw[i - 1]
            + 17.0 * f_raw[i]
            + 12.0 * f_raw[i + 1]
            - 3.0 * f_raw[i + 2])
            / 35.0;
    }

    // Central-difference negative derivative `-dF/dT`
    let temps = &trace.temperatures_c[..n];
    let mut deriv = vec![0.0_f64; n];
    for i in 1..(n - 1) {
        let dt = (temps[i + 1] - temps[i - 1]).max(1e-6);
        deriv[i] = (-(f_smooth[i + 1] - f_smooth[i - 1]) / dt).max(0.0);
    }
    deriv[0] = deriv[1];
    deriv[n - 1] = deriv[n - 2];

    let max_d = deriv.iter().copied().fold(0.0_f64, f64::max);
    if max_d < 5.0 {
        return MeltCurveResult {
            well_id: trace.well_id.clone(),
            sample_id: trace.sample_id.clone(),
            target_gene: trace.target_gene.clone(),
            temperatures_c: temps.to_vec(),
            neg_df_dt: deriv,
            primary_tm_c: None,
            peaks: Vec::new(),
            specificity_status: "FLAT / NO MELT PEAK",
        };
    }

    let mut raw_peaks: Vec<(usize, f64, f64)> = Vec::new();
    for i in 2..(n - 2) {
        let y0 = deriv[i];
        if y0 >= 0.12 * max_d
            && y0 >= deriv[i - 1]
            && y0 > deriv[i + 1]
            && y0 >= deriv[i - 2]
            && y0 >= deriv[i + 2]
        {
            // Parabolic interpolation around `(i - 1, i, i + 1)`
            let ym = deriv[i - 1];
            let yp = deriv[i + 1];
            let denom = 2.0 * (2.0 * y0 - ym - yp);
            let step = 0.5 * (temps[i + 1] - temps[i - 1]);
            let offset = if denom.abs() > 1e-12 {
                ((ym - yp) / denom).clamp(-0.5, 0.5)
            } else {
                0.0
            };
            let tm_interp = temps[i] + offset * step;
            raw_peaks.push((i, tm_interp, y0));
        }
    }

    raw_peaks.sort_by(|a, b| b.2.partial_cmp(&a.2).unwrap());
    let mut peaks = Vec::with_capacity(raw_peaks.len());
    for &(idx, tm_c, pk_h) in &raw_peaks {
        let half = 0.5 * pk_h;
        let mut l = idx;
        while l > 0 && deriv[l] > half {
            l -= 1;
        }
        let mut r = idx;
        while r + 1 < n && deriv[r] > half {
            r += 1;
        }
        let fwhm = (temps[r] - temps[l]).max(0.2);
        peaks.push(MeltPeakCall {
            tm_celsius: tm_c,
            peak_neg_df_dt: pk_h,
            fwhm_celsius: fwhm,
            relative_height_pct: (pk_h / max_d) * 100.0,
            is_primer_dimer_suspect: tm_c < dimer_cutoff_c,
        });
    }

    let primary_tm = peaks.first().map(|p| p.tm_celsius);
    let has_dimer = peaks.iter().any(|p| p.is_primer_dimer_suspect);
    let status = if peaks.len() == 1 && !has_dimer {
        "CLEAN SINGLE AMPLICON"
    } else if has_dimer && peaks.len() == 1 {
        "PRIMER-DIMER ARTIFACT (Low Tm)"
    } else if has_dimer {
        "AMPLICON + PRIMER-DIMER SHOULDER"
    } else {
        "MULTIMODAL / NON-SPECIFIC PRODUCT"
    };

    MeltCurveResult {
        well_id: trace.well_id.clone(),
        sample_id: trace.sample_id.clone(),
        target_gene: trace.target_gene.clone(),
        temperatures_c: temps.to_vec(),
        neg_df_dt: deriv,
        primary_tm_c: primary_tm,
        peaks,
        specificity_status: status,
    }
}

/// Generates 4 realistic SYBR Green I dissociation curves (`65.0 .. 95.0 °C` in `0.2 °C` steps):
/// - `HSP70 (HS_1)` — Single sharp amplicon peak at `85.4 °C`
/// - `ACTB (Ctrl_1)` — Single sharp reference peak at `87.2 °C`
/// - `IL6 (HS_1)` — Clean target amplicon peak at `83.6 °C`
/// - `IL6 (NTC_H02)` — Classic `73.8 °C` primer-dimer melt peak explaining the late `Cq = 37.85` NTC signal!
pub fn demo_sybr_melt_traces() -> Vec<MeltTrace> {
    let temps: Vec<f64> = (0..=150).map(|k| 65.0 + 0.2 * (k as f64)).collect();

    let synth = |well: &str, samp: &str, gene: &str, transitions: &[(f64, f64, f64)]| -> MeltTrace {
        let mut fl = Vec::with_capacity(temps.len());
        for &t in &temps {
            // Linear thermal quenching envelope + Boltzmann sigmoid dsDNA melting drops
            let quench = 1.0 - 0.0085 * (t - 65.0);
            let mut sig = 45.0;
            for &(amp, tm, width) in transitions {
                let ds_frac = 1.0 / (1.0 + ((t - tm) / width).exp());
                sig += amp * quench.max(0.25) * ds_frac;
            }
            fl.push(sig);
        }
        MeltTrace {
            well_id: well.to_string(),
            sample_id: samp.to_string(),
            target_gene: gene.to_string(),
            temperatures_c: temps.clone(),
            fluorescence_rfu: fl,
        }
    };

    vec![
        synth("C07", "HS_1", "HSP70", &[(1_850.0, 85.4, 0.52)]),
        synth("E01", "Ctrl_1", "ACTB", &[(1_920.0, 87.2, 0.50)]),
        synth("D07", "HS_1", "IL6", &[(1_640.0, 83.6, 0.55)]),
        synth("H02", "NTC_IL6", "IL6", &[(620.0, 73.8, 0.68)]),
    ]
}

/// Renders a multi-trace `-dF/dT` SYBR Green dissociation plot (`.svg`) with vertical primer-dimer
/// threshold shading (`< 77.5 °C`) and annotated `Tm` callouts.
pub fn render_melt_curves_svg(results: &[MeltCurveResult], title: &str) -> String {
    let width = 840.0_f64;
    let height = 440.0_f64;
    let left = 74.0_f64;
    let top = 56.0_f64;
    let pw = 530.0_f64;
    let ph = 320.0_f64;

    let t_min = 65.0_f64;
    let t_max = 95.0_f64;
    let mut y_max = 50.0_f64;
    for r in results {
        for &d in &r.neg_df_dt {
            y_max = y_max.max(d * 1.12);
        }
    }

    let map_x = |tc: f64| left + ((tc - t_min) / (t_max - t_min)).clamp(0.0, 1.0) * pw;
    let map_y = |v: f64| top + ph - (v / y_max).clamp(0.0, 1.0) * ph;

    const COLORS: &[&str] = &["#2563eb", "#16a34a", "#9333ea", "#dc2626", "#d97706"];
    let mut svg = String::with_capacity(12_288);
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
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — SYBR Green Derivative Melt Curves (-dF/dT)</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    // Primer-dimer caution zone (< 77.5 °C)
    let x_dimer = map_x(77.5);
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{:.1}\" height=\"{ph:.1}\" fill=\"#fee2e2\" fill-opacity=\"0.45\"/>",
        x_dimer - left
    );
    let _ = writeln!(
        svg,
        "  <line x1=\"{x_dimer:.1}\" y1=\"{top:.1}\" x2=\"{x_dimer:.1}\" y2=\"{:.1}\" stroke=\"#ef4444\" stroke-dasharray=\"5,4\"/>",
        top + ph
    );

    for (idx, res) in results.iter().enumerate() {
        let col = COLORS[idx % COLORS.len()];
        let mut poly = String::new();
        for (&tc, &dv) in res.temperatures_c.iter().zip(&res.neg_df_dt) {
            let _ = write!(poly, "{:.2},{:.2} ", map_x(tc), map_y(dv));
        }
        let _ = writeln!(
            svg,
            "  <polyline points=\"{}\" fill=\"none\" stroke=\"{col}\" stroke-width=\"2.3\"/>",
            poly.trim()
        );

        let lx = left + pw + 16.0;
        let ly = top + 24.0 + (idx as f64) * 48.0;
        let tm_txt = res
            .primary_tm_c
            .map(|t| format!("Tm = {t:.2} °C"))
            .unwrap_or_else(|| String::from("No Peak"));
        let _ = writeln!(
            svg,
            "  <text x=\"{lx:.1}\" y=\"{ly:.1}\" font-family=\"sans-serif\" font-size=\"11.5\" font-weight=\"700\" fill=\"{col}\">{} ({})</text>",
            xml_escape(&res.target_gene),
            xml_escape(&res.sample_id)
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#334155\">{tm_txt} | {}</text>",
            ly + 15.0,
            xml_escape(res.specificity_status)
        );
    }

    svg.push_str("</svg>\n");
    svg
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn calls_exact_amplicon_and_primer_dimer_melt_peaks() {
        let traces = demo_sybr_melt_traces();
        let hsp = analyze_melt_curve(&traces[0], 77.5);
        assert_eq!(hsp.specificity_status, "CLEAN SINGLE AMPLICON");
        assert!((hsp.primary_tm_c.unwrap() - 85.4).abs() < 0.35);

        let ntc_dimer = analyze_melt_curve(&traces[3], 77.5);
        assert_eq!(ntc_dimer.specificity_status, "PRIMER-DIMER ARTIFACT (Low Tm)");
        assert!((ntc_dimer.primary_tm_c.unwrap() - 73.8).abs() < 0.35);
    }
}
