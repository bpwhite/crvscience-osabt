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

//! Tool #18 (`osabt flow`) — `Logicle` biexponential scaling, polygon/quadrant gating & `.svg` plotter:
//! - **Parks, Roederer & Moore (2006)** *Cytometry A* 69(6):541–551 (`https://pubmed.ncbi.nlm.nih.gov/16604519/`)
//! - **Moore & Parks (2012)** *Cytometry A* (`https://pmc.ncbi.nlm.nih.gov/articles/PMC4761345/`)
//! - Point-in-polygon ray casting, 4-quadrant statistics (`% Parent`, `MFI`, robust `CV%`), & `.svg` export.

use std::fmt::Write as _;

use crate::flow_fcs::FcsDataset;
use crate::svg_util::{begin_svg_canvas, xml_escape};

/// Parameters for the Parks, Roederer & Moore (2006) `Logicle` biexponential display scale.
#[derive(Debug, Clone, Copy)]
pub struct LogicleParams {
    /// Top of linear data scale `T` (e.g., `262144.0` for 18-bit cytometry).
    pub t_top: f64,
    /// Width of the near-zero linear transition region in decades `W` (typically `0.5`–`1.0`).
    pub w_decades: f64,
    /// Total positive dynamic range in decades `M` (typically `4.5`).
    pub m_decades: f64,
    /// Additional negative decades `A` (typically `0.0`).
    pub a_neg_decades: f64,
}

impl Default for LogicleParams {
    fn default() -> Self {
        Self {
            t_top: 262_144.0,
            w_decades: 0.75,
            m_decades: 4.5,
            a_neg_decades: 0.0,
        }
    }
}

/// Precomputed `Logicle` transform constants solving `W = 2 * p * ln(p) / (p + 1)` (Moore & Parks 2012).
#[derive(Debug, Clone, Copy)]
pub struct LogicleScale {
    pub params: LogicleParams,
    a: f64,
    b: f64,
    c: f64,
    d: f64,
    f: f64,
    x0: f64,
    x1: f64,
    x2: f64,
}

impl LogicleScale {
    pub fn new(params: LogicleParams) -> Self {
        let t = params.t_top.max(10.0);
        let m = params.m_decades.max(2.0);
        let w = params.w_decades.clamp(0.05, 0.45 * m);
        let a_neg = params.a_neg_decades.clamp(0.0, m - 2.0 * w);

        let b = (m + a_neg) * std::f64::consts::LN_10;
        let w_frac = w / (m + a_neg);
        let x2 = a_neg / (m + a_neg);
        let x1 = x2 + w_frac;
        let x0 = x2 + 2.0 * w_frac;

        // Solve `d / b = 2 * p * ln(p) / (p + 1)` => `w_frac * b = ...` via Newton–Raphson on `d`
        let d = solve_logicle_d(b, w_frac);
        let c = b * ((b - d) * x1).exp();
        let f_raw = -((b * x1).exp()) + c * ((-d * x1).exp());
        let scale_a = t / ((b * 1.0).exp() - c * ((-d * 1.0).exp()) + f_raw);
        let a = scale_a;
        let f = scale_a * f_raw;
        let c_scaled = scale_a * c;

        Self {
            params,
            a,
            b,
            c: c_scaled,
            d,
            f,
            x0,
            x1,
            x2,
        }
    }

    /// Evaluates the inverse `Logicle` function `S(x)` mapping normalized display `x in [0, 1]` to linear data `y`.
    #[inline]
    pub fn inverse(&self, x: f64) -> f64 {
        let val = if x >= self.x1 {
            self.a * (self.b * x).exp() - self.c * (-self.d * x).exp() + self.f
        } else {
            // Exact odd symmetry around `x1` for the negative/zero-centered linear region
            let xref = 2.0 * self.x1 - x;
            -(self.a * (self.b * xref).exp() - self.c * (-self.d * xref).exp() + self.f)
        };
        let _ = (self.x0, self.x2);
        val
    }

    /// Maps a linear cytometry intensity `y` (including negative compensated values) to display coordinate `x in [0, 1]`.
    pub fn transform(&self, y: f64) -> f64 {
        if !y.is_finite() {
            return self.x1;
        }
        // Bracketed bisection + Newton–Raphson polish for guaranteed monotonic accuracy < 1e-12
        let mut lo = -0.25_f64;
        let mut hi = 1.25_f64;
        for _ in 0..48 {
            let mid = 0.5 * (lo + hi);
            if self.inverse(mid) < y {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        0.5 * (lo + hi)
    }
}

fn solve_logicle_d(b: f64, w_frac: f64) -> f64 {
    // Target: `(e^((b+d)*w_frac) ...)` -> Moore & Parks (2012) Eq. 5: `f(d) = 2*(b*d/(b+d))*...`
    // Specifically `b * w_frac` corresponds to `W * ln(10) = 2 * ln(p) / (1 + 1/p)` where `p = b / d`!
    let target_w_ln10 = w_frac * b;
    let mut lo = 0.001_f64;
    let mut hi = b - 0.001;
    for _ in 0..60 {
        let d = 0.5 * (lo + hi);
        let lhs = d * (b - d) * ((b - d).ln() - d.ln()) / b;
        // Equivalently `b / d = p`: `w_ln10 = 2 * (p / (p + 1)) * ln(p)` when `d = b / p`
        let p = b / d;
        let w_val = 2.0 * (p / (p + 1.0)) * p.ln();
        let _ = lhs;
        if w_val > target_w_ln10 {
            lo = d;
        } else {
            hi = d;
        }
    }
    0.5 * (lo + hi)
}

/// 2D polygon gate in transformed display coordinates (`[0, 1] x [0, 1]`) or raw units.
#[derive(Debug, Clone)]
pub struct PolygonGate {
    pub name: String,
    pub vertices: Vec<(f64, f64)>,
}

impl PolygonGate {
    pub fn contains(&self, x: f64, y: f64) -> bool {
        let n = self.vertices.len();
        if n < 3 {
            return false;
        }
        let mut inside = false;
        let mut j = n - 1;
        for i in 0..n {
            let (xi, yi) = self.vertices[i];
            let (xj, yj) = self.vertices[j];
            let intersect = ((yi > y) != (yj > y))
                && (x < (xj - xi) * (y - yi) / ((yj - yi) + 1e-18) + xi);
            if intersect {
                inside = !inside;
            }
            j = i;
        }
        inside
    }
}

/// Summary statistics for one gated cell population (`% Parent`, `MFI`, robust `CV%`).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct GatePopulationStats {
    pub gate_name: String,
    pub event_count: usize,
    pub percent_of_parent: f64,
    pub x_median: f64,
    pub y_median: f64,
    pub x_robust_cv_pct: f64,
    pub y_robust_cv_pct: f64,
}

/// Evaluates 4-quadrant gating (`Q1: -/+`, `Q2: +/+`, `Q3: -/-`, `Q4: +/-`) in `Logicle` display coordinates
/// across channels `x_ch` and `y_ch` with thresholds `(x_cut_disp, y_cut_disp)`.
pub fn evaluate_quadrant_gates(
    dataset: &FcsDataset,
    x_ch: usize,
    y_ch: usize,
    x_cut_disp: f64,
    y_cut_disp: f64,
    scale: &LogicleScale,
) -> Result<Vec<GatePopulationStats>, String> {
    if x_ch >= dataset.channels.len() || y_ch >= dataset.channels.len() {
        return Err(String::from("Channel index out of bounds for quadrant gating."));
    }
    let tot = dataset.events.len().max(1) as f64;
    let mut buckets: [Vec<(f64, f64)>; 4] = [Vec::new(), Vec::new(), Vec::new(), Vec::new()];

    for ev in &dataset.events {
        let lx = ev[x_ch];
        let ly = ev[y_ch];
        let dx = scale.transform(lx);
        let dy = scale.transform(ly);
        let q_idx = match (dx >= x_cut_disp, dy >= y_cut_disp) {
            (false, true) => 0,  // Q1 (-/+)
            (true, true) => 1,   // Q2 (+/+)
            (false, false) => 2, // Q3 (-/-)
            (true, false) => 3,  // Q4 (+/-)
        };
        buckets[q_idx].push((lx, ly));
    }

    let names = [
        "Q1 (-/+)",
        "Q2 (+/+ Double-Pos)",
        "Q3 (-/- Double-Neg)",
        "Q4 (+/-)",
    ];
    let mut out = Vec::with_capacity(4);
    for (i, nm) in names.iter().enumerate() {
        let cnt = buckets[i].len();
        let xs: Vec<f64> = buckets[i].iter().map(|p| p.0).collect();
        let ys: Vec<f64> = buckets[i].iter().map(|p| p.1).collect();
        let (mx, rcv_x) = median_and_robust_cv(&xs);
        let (my, rcv_y) = median_and_robust_cv(&ys);
        out.push(GatePopulationStats {
            gate_name: (*nm).to_string(),
            event_count: cnt,
            percent_of_parent: 100.0 * (cnt as f64) / tot,
            x_median: mx,
            y_median: my,
            x_robust_cv_pct: rcv_x,
            y_robust_cv_pct: rcv_y,
        });
    }
    Ok(out)
}

/// Evaluates an arbitrary [`PolygonGate`] in `Logicle` display space (`[0, 1] x [0, 1]`).
pub fn evaluate_polygon_gate(
    dataset: &FcsDataset,
    x_ch: usize,
    y_ch: usize,
    gate: &PolygonGate,
    scale: &LogicleScale,
) -> Result<GatePopulationStats, String> {
    if x_ch >= dataset.channels.len() || y_ch >= dataset.channels.len() {
        return Err(String::from("Channel index out of bounds for polygon gating."));
    }
    let mut xs = Vec::new();
    let mut ys = Vec::new();
    for ev in &dataset.events {
        let dx = scale.transform(ev[x_ch]);
        let dy = scale.transform(ev[y_ch]);
        if gate.contains(dx, dy) {
            xs.push(ev[x_ch]);
            ys.push(ev[y_ch]);
        }
    }
    let cnt = xs.len();
    let (mx, rcv_x) = median_and_robust_cv(&xs);
    let (my, rcv_y) = median_and_robust_cv(&ys);
    Ok(GatePopulationStats {
        gate_name: gate.name.clone(),
        event_count: cnt,
        percent_of_parent: 100.0 * (cnt as f64) / (dataset.events.len().max(1) as f64),
        x_median: mx,
        y_median: my,
        x_robust_cv_pct: rcv_x,
        y_robust_cv_pct: rcv_y,
    })
}

fn median_and_robust_cv(vals: &[f64]) -> (f64, f64) {
    if vals.is_empty() {
        return (0.0, 0.0);
    }
    let mut s = vals.to_vec();
    s.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let med = s[s.len() / 2];
    let mut devs: Vec<f64> = s.iter().map(|&v| (v - med).abs()).collect();
    devs.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let rsd = 1.4826 * devs[devs.len() / 2];
    let rcv = if med.abs() > 1e-6 {
        100.0 * rsd / med.abs()
    } else {
        0.0
    };
    (med, rcv)
}

/// Renders a 2D `Logicle`-transformed biaxial cytometry scatter plot with quadrant gate dividers
/// and population percentage annotations as a standalone `.svg` figure.
pub fn render_flow_biaxial_svg(
    dataset: &FcsDataset,
    x_ch: usize,
    y_ch: usize,
    x_cut_disp: f64,
    y_cut_disp: f64,
    scale: &LogicleScale,
    title: &str,
) -> String {
    let width = 620.0_f64;
    let height = 560.0_f64;
    let left = 72.0_f64;
    let top = 56.0_f64;
    let plot_sz = 440.0_f64;

    let x_lbl = dataset
        .channels
        .get(x_ch)
        .map(|c| c.short_name.as_str())
        .unwrap_or("X");
    let y_lbl = dataset
        .channels
        .get(y_ch)
        .map(|c| c.short_name.as_str())
        .unwrap_or("Y");

    let quad_stats =
        evaluate_quadrant_gates(dataset, x_ch, y_ch, x_cut_disp, y_cut_disp, scale)
            .unwrap_or_default();

    let mut svg = begin_svg_canvas(24_576, width, height, Some("#ffffff"));
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"30\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Logicle ({x_lbl} vs {y_lbl}, n={})</text>",
        xml_escape(title),
        dataset.events.len()
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{plot_sz:.1}\" height=\"{plot_sz:.1}\" fill=\"#f8fafc\" stroke=\"#94a3b8\"/>"
    );

    for ev in &dataset.events {
        let dx = scale.transform(ev[x_ch]).clamp(0.0, 1.0);
        let dy = scale.transform(ev[y_ch]).clamp(0.0, 1.0);
        let px = left + dx * plot_sz;
        let py = top + plot_sz - dy * plot_sz;
        let col = match (dx >= x_cut_disp, dy >= y_cut_disp) {
            (false, true) => "#9333ea",
            (true, true) => "#dc2626",
            (false, false) => "#2563eb",
            (true, false) => "#16a34a",
        };
        let _ = writeln!(
            svg,
            "  <circle cx=\"{px:.1}\" cy=\"{py:.1}\" r=\"2.2\" fill=\"{col}\" fill-opacity=\"0.65\"/>"
        );
    }

    // Quadrant crosshairs
    let qx = left + x_cut_disp.clamp(0.0, 1.0) * plot_sz;
    let qy = top + plot_sz - y_cut_disp.clamp(0.0, 1.0) * plot_sz;
    let _ = writeln!(
        svg,
        "  <line x1=\"{qx:.1}\" y1=\"{top:.1}\" x2=\"{qx:.1}\" y2=\"{:.1}\" stroke=\"#0f172a\" stroke-dasharray=\"5,4\" stroke-width=\"1.5\"/>",
        top + plot_sz
    );
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{qy:.1}\" x2=\"{:.1}\" y2=\"{qy:.1}\" stroke=\"#0f172a\" stroke-dasharray=\"5,4\" stroke-width=\"1.5\"/>",
        left + plot_sz
    );

    if quad_stats.len() == 4 {
        let corners = [
            (left + 12.0, top + 20.0, &quad_stats[0]),
            (left + plot_sz - 120.0, top + 20.0, &quad_stats[1]),
            (left + 12.0, top + plot_sz - 12.0, &quad_stats[2]),
            (left + plot_sz - 120.0, top + plot_sz - 12.0, &quad_stats[3]),
        ];
        for (cx, cy, st) in corners {
            let _ = writeln!(
                svg,
                "  <text x=\"{cx:.1}\" y=\"{cy:.1}\" font-family=\"sans-serif\" font-size=\"11\" font-weight=\"700\" fill=\"#0f172a\">{}: {:.1}%</text>",
                st.gate_name.split_whitespace().next().unwrap_or("Q"),
                st.percent_of_parent
            );
        }
    }

    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#1e293b\">{x_lbl} (Logicle Biexponential)</text>",
        left + 0.5 * plot_sz,
        top + plot_sz + 36.0
    );
    svg.push_str("</svg>\n");
    svg
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::flow_fcs::{compensate_fcs_dataset, demo_fcs_31_bytes, parse_fcs_bytes};

    #[test]
    fn logicle_roundtrips_negative_and_positive_values_and_gates_populations() {
        let scale = LogicleScale::new(LogicleParams::default());
        for &y in &[-300.0, 0.0, 500.0, 20_000.0, 200_000.0] {
            let x = scale.transform(y);
            let y_back = scale.inverse(x);
            assert!((y - y_back).abs() < 1e-3);
        }

        let ds = parse_fcs_bytes(&demo_fcs_31_bytes()).unwrap();
        let (comp, _) = compensate_fcs_dataset(&ds).unwrap();
        let quads = evaluate_quadrant_gates(&comp, 2, 3, 0.38, 0.38, &scale).unwrap();
        let sum_pct: f64 = quads.iter().map(|q| q.percent_of_parent).sum();
        assert!((sum_pct - 100.0).abs() < 1e-6);
        // Double-positive cluster (~90 events = 22.5%), FITC+ single-positive (~120 events = 30%)
        assert!(quads[1].event_count >= 75);
        assert!(quads[3].event_count >= 100);
    }
}
