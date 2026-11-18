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

//! Publication-quality vector `.svg` generators for `osabt plasmid`:
//! - Circular plasmid vector maps (`render_plasmid_map_svg`) with directional feature arcs,
//!   GC-content track, unique restriction site callouts, and bound primer arrows
//! - Multi-lane virtual agarose gel electrophoresis figures (`render_virtual_gel_svg`)

use std::f64::consts::PI;

use crate::gibson::{CircularBindingSite, Strand};
use crate::plasmid_digest::find_unique_cutters;
use crate::plasmid_io::PlasmidRecord;
use crate::gc_fraction;
use crate::svg_util::{format_polyline_points, xml_escape};

/// Input callout descriptor for circular plasmid hemisphere sweep-line relaxation.
#[derive(Debug, Clone)]
pub struct CircularCalloutInput {
    pub label: String,
    /// Clockwise angle in `[0, 2*PI)` from 12 o'clock.
    pub angle_rad: f64,
    /// Radius where the leader line anchors (backbone or feature outer edge).
    pub anchor_radius: f64,
    pub is_feature: bool,
    pub color_hex: String,
    pub source_index: usize,
}

/// Relaxed circular callout output with 3-point elbow leader coordinates.
#[derive(Debug, Clone)]
pub struct CircularCalloutPlaced {
    pub input: CircularCalloutInput,
    /// Anchor point on backbone or feature arc.
    pub p0: (f64, f64),
    /// Radial elbow point outside feature tracks.
    pub p1: (f64, f64),
    /// Relaxed horizontal label anchor point.
    pub p2: (f64, f64),
    pub is_right_hemisphere: bool,
}

/// Solves collision-free callout positions across Left and Right hemispheres using
/// sweep-line vertical relaxation, ring-exclusion horizontal push-out, and 3-point elbow leaders.
pub fn layout_circular_callouts(
    cx: f64,
    cy: f64,
    r_feat_outer: f64,
    canvas_bounds: (f64, f64, f64, f64),
    min_gap_y: f64,
    inputs: Vec<CircularCalloutInput>,
) -> Vec<CircularCalloutPlaced> {
    let (min_x, min_y, max_x, max_y) = canvas_bounds;
    let y_lo = (min_y + 12.0).min(cy);
    let y_hi = (max_y - 12.0).max(y_lo + 10.0);
    let r_elbow = r_feat_outer + 16.0;
    let r_ring_excl = r_feat_outer + 16.0;

    let mut right_items: Vec<(f64, CircularCalloutInput)> = Vec::new();
    let mut left_items: Vec<(f64, CircularCalloutInput)> = Vec::new();

    for inp in inputs {
        let ang = inp.angle_rad.rem_euclid(2.0 * PI);
        let theta = ang - PI * 0.5;
        let y_nat = cy + r_elbow * theta.sin();
        if theta.cos() >= -1e-9 {
            right_items.push((y_nat, inp));
        } else {
            left_items.push((y_nat, inp));
        }
    }

    right_items.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.1.angle_rad.partial_cmp(&b.1.angle_rad).unwrap_or(std::cmp::Ordering::Equal))
    });
    left_items.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| b.1.angle_rad.partial_cmp(&a.1.angle_rad).unwrap_or(std::cmp::Ordering::Equal))
    });

    let mut placed = Vec::with_capacity(right_items.len() + left_items.len());
    for (is_right, items) in [(true, right_items), (false, left_items)] {
        let count = items.len();
        if count == 0 {
            continue;
        }
        let avail_h = (y_hi - y_lo).max(1.0);
        let gap = if count > 1 {
            min_gap_y.min(avail_h / ((count - 1) as f64))
        } else {
            min_gap_y
        };
        let mut ys: Vec<f64> = items.iter().map(|(y, _)| y.clamp(y_lo, y_hi)).collect();
        if count > 1 {
            for i in 1..count {
                ys[i] = ys[i].max(ys[i - 1] + gap);
            }
            if ys[count - 1] > y_hi {
                ys[count - 1] = y_hi;
            }
            for i in (0..count - 1).rev() {
                ys[i] = ys[i].min(ys[i + 1] - gap);
            }
            if ys[0] < y_lo {
                ys[0] = y_lo;
                for i in 1..count {
                    ys[i] = ys[i].max(ys[i - 1] + gap).min(y_hi);
                }
            }
        }

        for ((_, inp), y_rel) in items.into_iter().zip(ys) {
            let ang = inp.angle_rad.rem_euclid(2.0 * PI);
            let p0 = polar_xy(cx, cy, inp.anchor_radius, ang);
            let p1 = polar_xy(cx, cy, r_elbow, ang);
            let dy = (y_rel - cy).abs();
            let ring_clear_dx = if dy < r_ring_excl {
                (r_ring_excl * r_ring_excl - dy * dy).sqrt() + 14.0
            } else {
                18.0
            };
            let elbow_dx = (p1.0 - cx).abs() + 10.0;
            let dx = ring_clear_dx.max(elbow_dx);
            let label_x = if is_right {
                let min_right_x = cx + ring_clear_dx;
                (cx + dx).min((max_x - 8.0).max(min_right_x))
            } else {
                let max_left_x = cx - ring_clear_dx;
                (cx - dx).max((min_x + 8.0).min(max_left_x))
            };
            placed.push(CircularCalloutPlaced {
                input: inp,
                p0,
                p1,
                p2: (label_x, y_rel),
                is_right_hemisphere: is_right,
            });
        }
    }

    placed
}

/// DNA molecular-weight ladder standard for virtual agarose gel electrophoresis.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum GelLadder {
    /// Standard 1 kb Plus ladder (`100 bp .. 10,000 bp`; bright reference bands at `500, 1000, 3000 bp`).
    OneKbPlus,
    /// Low-molecular-weight 100 bp ladder (`100 bp .. 1,500 bp`; bright bands at `500, 1000 bp`).
    HundredBp,
}

impl GelLadder {
    /// Returns `(fragment_bp, is_bright_reference_band)` rungs.
    pub fn rungs(self) -> &'static [(usize, bool)] {
        match self {
            Self::OneKbPlus => &[
                (10_000, false),
                (8_000, false),
                (6_000, false),
                (5_000, false),
                (4_000, false),
                (3_000, true),
                (2_000, false),
                (1_500, false),
                (1_000, true),
                (850, false),
                (650, false),
                (500, true),
                (400, false),
                (300, false),
                (200, false),
                (100, false),
            ],
            Self::HundredBp => &[
                (1_500, false),
                (1_200, false),
                (1_000, true),
                (900, false),
                (800, false),
                (700, false),
                (600, false),
                (500, true),
                (400, false),
                (300, false),
                (200, false),
                (100, false),
            ],
        }
    }
}

/// Sample lane rendered on a virtual agarose gel.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct GelLane {
    pub label: String,
    pub fragment_lengths_bp: Vec<usize>,
}

/// Computes normalized migration distance (`0.0` at loading well -> `1.0` at dye front)
/// for a linear dsDNA fragment of `bp` base pairs in an agarose gel of `agarose_pct` (`0.6..=3.0%`).
pub fn gel_migration_fraction(bp: usize, agarose_pct: f64) -> f64 {
    let ag = agarose_pct.clamp(0.5, 3.0);
    let log_bp = (bp.max(40) as f64).log10();
    // Higher agarose % compresses high-MW bands near the well and resolves 100-800 bp cleanly
    let log_max = (4.35 - 0.28 * (ag - 1.0)).max(3.4);
    let log_min = (1.75 - 0.15 * (ag - 1.0)).max(1.3);
    let raw = (log_max - log_bp) / (log_max - log_min);
    // Slight Ogston/reptation S-curve compression near high-MW limit
    raw.clamp(0.03, 0.96)
}


fn polar_xy(cx: f64, cy: f64, radius: f64, angle_rad: f64) -> (f64, f64) {
    // 0 bp sits at 12 o'clock (-PI/2), advancing clockwise
    let theta = angle_rad - PI * 0.5;
    (cx + radius * theta.cos(), cy + radius * theta.sin())
}

fn svg_arc_band_path(
    cx: f64,
    cy: f64,
    r_inner: f64,
    r_outer: f64,
    start_angle: f64,
    sweep_angle: f64,
    strand: Strand,
) -> String {
    let clamped_sweep = sweep_angle.clamp(0.03, 2.0 * PI - 0.01);
    let r_mid = 0.5 * (r_inner + r_outer);
    let head_rad = (10.0 / r_mid).min(clamped_sweep * 0.45);

    let (body_s, body_e, tip_angle) = match strand {
        Strand::Forward => (
            start_angle,
            start_angle + clamped_sweep - head_rad,
            start_angle + clamped_sweep,
        ),
        Strand::Reverse => (
            start_angle + head_rad,
            start_angle + clamped_sweep,
            start_angle,
        ),
    };

    let body_span = (body_e - body_s).max(0.005);
    let large_arc = if body_span > PI { 1 } else { 0 };

    let (x_out_s, y_out_s) = polar_xy(cx, cy, r_outer, body_s);
    let (x_out_e, y_out_e) = polar_xy(cx, cy, r_outer, body_e);
    let (x_in_e, y_in_e) = polar_xy(cx, cy, r_inner, body_e);
    let (x_in_s, y_in_s) = polar_xy(cx, cy, r_inner, body_s);
    let (x_tip, y_tip) = polar_xy(cx, cy, r_mid, tip_angle);

    match strand {
        Strand::Forward => {
            let (xf1, yf1) = polar_xy(cx, cy, r_outer + 3.0, body_e);
            let (xf2, yf2) = polar_xy(cx, cy, r_inner - 3.0, body_e);
            format!(
                "M {x_out_s:.2} {y_out_s:.2} \
                 A {r_outer:.2} {r_outer:.2} 0 {large_arc} 1 {x_out_e:.2} {y_out_e:.2} \
                 L {xf1:.2} {yf1:.2} L {x_tip:.2} {y_tip:.2} L {xf2:.2} {yf2:.2} \
                 L {x_in_e:.2} {y_in_e:.2} \
                 A {r_inner:.2} {r_inner:.2} 0 {large_arc} 0 {x_in_s:.2} {y_in_s:.2} Z"
            )
        }
        Strand::Reverse => {
            let (xr1, yr1) = polar_xy(cx, cy, r_outer + 3.0, body_s);
            let (xr2, yr2) = polar_xy(cx, cy, r_inner - 3.0, body_s);
            format!(
                "M {x_tip:.2} {y_tip:.2} L {xr1:.2} {yr1:.2} L {x_out_s:.2} {y_out_s:.2} \
                 A {r_outer:.2} {r_outer:.2} 0 {large_arc} 1 {x_out_e:.2} {y_out_e:.2} \
                 L {x_in_e:.2} {y_in_e:.2} \
                 A {r_inner:.2} {r_inner:.2} 0 {large_arc} 0 {x_in_s:.2} {y_in_s:.2} \
                 L {xr2:.2} {yr2:.2} Z"
            )
        }
    }
}

/// Renders a standalone, publication-ready circular plasmid vector map `.svg` showing annotated
/// features, sliding-window GC% inner ring, unique restriction cutters, and bound primers.
pub fn render_plasmid_map_svg(
    rec: &PlasmidRecord,
    primers: &[CircularBindingSite],
    include_attribution: bool,
) -> String {
    let width = 820.0_f64;
    let height = 660.0_f64;
    let cx = 410.0_f64;
    let cy = 325.0_f64;
    let r_backbone = 162.0_f64;
    let r_gc_base = 116.0_f64;
    let n = rec.sequence.len().max(1);

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"#ffffff\"/>\n\
         <style>\n\
           .title {{ font-family: system-ui, -apple-system, sans-serif; font-size: 18px; font-weight: 700; fill: #0f172a; }}\n\
           .sub   {{ font-family: ui-monospace, monospace; font-size: 12px; fill: #475569; }}\n\
           .feat  {{ font-family: system-ui, sans-serif; font-size: 11px; font-weight: 600; fill: #1e293b; }}\n\
           .enz   {{ font-family: ui-monospace, monospace; font-size: 11px; fill: #0369a1; }}\n\
           .prm   {{ font-family: ui-monospace, monospace; font-size: 10px; fill: #be185d; }}\n\
           .foot  {{ font-family: system-ui, sans-serif; font-size: 10px; fill: #64748b; }}\n\
         </style>\n"
    );

    // 1. Inner GC% deviation polyline ring
    let bins = 120usize.min(n);
    let step = (n / bins).max(1);
    let mean_gc = gc_fraction(&rec.sequence);
    let mut gc_pts = String::new();
    for b_idx in 0..=bins {
        let pos = (b_idx * step) % n;
        let win_end = (pos + step).min(n);
        let local_gc = gc_fraction(&rec.sequence[pos..win_end.max(pos + 1)]);
        let r = r_gc_base + (local_gc - mean_gc) * 38.0;
        let ang = (pos as f64 / n as f64) * 2.0 * PI;
        let (gx, gy) = polar_xy(cx, cy, r, ang);
        if b_idx == 0 {
            gc_pts.push_str(&format!("M {gx:.2} {gy:.2} "));
        } else {
            gc_pts.push_str(&format!("L {gx:.2} {gy:.2} "));
        }
    }
    svg.push_str(&format!(
        "<circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{r_gc_base:.1}\" fill=\"none\" stroke=\"#e2e8f0\" stroke-dasharray=\"3 3\"/>\n\
         <path d=\"{gc_pts}Z\" fill=\"none\" stroke=\"#94a3b8\" stroke-width=\"1.2\"/>\n"
    ));

    // 2. Main dsDNA circular backbone
    svg.push_str(&format!(
        "<circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#1e293b\" stroke-width=\"3.0\"/>\n\
         <circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{:.1}\" fill=\"none\" stroke=\"#64748b\" stroke-width=\"1.2\"/>\n",
        r_backbone,
        r_backbone - 4.0
    ));

    // 3. Coordinate ruler ticks (4 quadrants)
    for q in 0..4 {
        let bp_pos = (q * n) / 4;
        let ang = (bp_pos as f64 / n as f64) * 2.0 * PI;
        let (x1, y1) = polar_xy(cx, cy, r_backbone - 10.0, ang);
        let (x2, y2) = polar_xy(cx, cy, r_backbone - 4.0, ang);
        let (xt, yt) = polar_xy(cx, cy, r_backbone - 22.0, ang);
        svg.push_str(&format!(
            "<line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#475569\" stroke-width=\"1.2\"/>\n\
             <text class=\"sub\" x=\"{xt:.1}\" y=\"{yt:.1}\" text-anchor=\"middle\" dominant-baseline=\"central\">{}</text>\n",
            bp_pos + 1
        ));
    }

    // 4. Annotated directional feature arcs + unified callout pool
    let mut callout_inputs: Vec<CircularCalloutInput> = Vec::new();
    for (idx, feat) in rec.features.iter().enumerate() {
        let start_ang = (feat.start_0 as f64 / n as f64) * 2.0 * PI;
        let span_bp = feat.span_bp(n);
        let sweep = (span_bp as f64 / n as f64) * 2.0 * PI;
        let d = svg_arc_band_path(
            cx,
            cy,
            r_backbone + 6.0,
            r_backbone + 22.0,
            start_ang,
            sweep,
            feat.strand,
        );
        svg.push_str(&format!(
            "<path d=\"{d}\" fill=\"{}\" fill-opacity=\"0.88\" stroke=\"#0f172a\" stroke-width=\"0.8\"/>\n",
            xml_escape(&feat.color_hex),
        ));
        callout_inputs.push(CircularCalloutInput {
            label: feat.name.clone(),
            angle_rad: start_ang + sweep * 0.5,
            anchor_radius: r_backbone + 22.0,
            is_feature: true,
            color_hex: feat.color_hex.clone(),
            source_index: idx,
        });
    }

    // 5. Unique restriction cutters added to same unified callout pool
    let uniques = find_unique_cutters(rec);
    for (u_idx, site) in uniques.iter().take(16).enumerate() {
        let ang = (site.cut_top_0 as f64 / n as f64) * 2.0 * PI;
        callout_inputs.push(CircularCalloutInput {
            label: format!("{} ({})", site.enzyme_name, site.cut_top_0 + 1),
            angle_rad: ang,
            anchor_radius: r_backbone,
            is_feature: false,
            color_hex: String::from("#0284c7"),
            source_index: u_idx,
        });
    }

    let placed_callouts = layout_circular_callouts(
        cx,
        cy,
        r_backbone + 24.0,
        (120.0, 22.0, width - 120.0, height - 30.0),
        16.0,
        callout_inputs,
    );
    for c in &placed_callouts {
        let pts = format_polyline_points(&[c.p0, c.p1, c.p2]);
        let anchor = if c.is_right_hemisphere { "start" } else { "end" };
        let dx = if c.is_right_hemisphere { 4.0 } else { -4.0 };
        if c.input.is_feature {
            svg.push_str(&format!(
                "<polyline fill=\"none\" stroke=\"#64748b\" stroke-width=\"1.0\" points=\"{pts}\"/>\n\
                 <text class=\"feat\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"{anchor}\" dominant-baseline=\"central\">{}</text>\n",
                c.p2.0 + dx,
                c.p2.1,
                xml_escape(&c.input.label),
            ));
        } else {
            svg.push_str(&format!(
                "<polyline fill=\"none\" stroke=\"#0284c7\" stroke-width=\"1.0\" stroke-dasharray=\"2 2\" points=\"{pts}\"/>\n\
                 <text class=\"enz\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"{anchor}\" dominant-baseline=\"central\">{}</text>\n",
                c.p2.0 + dx,
                c.p2.1,
                xml_escape(&c.input.label),
            ));
        }
    }

    // 6. Bound primer arrows (with 5' tail indicator)
    for (p_idx, p) in primers.iter().take(12).enumerate() {
        let start_ang = (p.bind_start_0 as f64 / n as f64) * 2.0 * PI;
        let sweep = ((p.bound_len.max(12) as f64 / n as f64) * 2.0 * PI).max(0.06);
        let r_p = r_backbone - 34.0 - ((p_idx % 2) as f64) * 10.0;
        let d = svg_arc_band_path(cx, cy, r_p - 4.0, r_p + 4.0, start_ang, sweep, p.strand);
        let (tx, ty) = polar_xy(cx, cy, r_p - 12.0, start_ang + sweep * 0.5);
        svg.push_str(&format!(
            "<path d=\"{d}\" fill=\"#ec4899\"/>\n\
             <text class=\"prm\" x=\"{tx:.1}\" y=\"{ty:.1}\" text-anchor=\"middle\" dominant-baseline=\"central\">{}</text>\n",
            xml_escape(&p.oligo_name)
        ));
    }

    // 7. Central plasmid nameplate
    let topo = if rec.circular { "Circular" } else { "Linear" };
    svg.push_str(&format!(
        "<text class=\"title\" x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>\n\
         <text class=\"sub\" x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{} bp • {} • {:.1}% GC</text>\n",
        cy - 8.0,
        xml_escape(&rec.name),
        cy + 14.0,
        rec.sequence.len(),
        topo,
        mean_gc * 100.0
    ));

    if include_attribution {
        svg.push_str(&format!(
            "<text class=\"foot\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\">Generated by OSABT (osabt plasmid) • crvscience.com • GPL-3.0</text>\n",
            width - 18.0,
            height - 14.0
        ));
    }
    svg.push_str("</svg>\n");
    svg
}

/// Renders a multi-lane virtual agarose gel `.svg` next to a reference molecular-weight ladder.
pub fn render_virtual_gel_svg(
    ladder: GelLadder,
    lanes: &[GelLane],
    agarose_pct: f64,
    darkroom_theme: bool,
    include_attribution: bool,
) -> String {
    let total_lanes = 1 + lanes.len();
    let lane_w = 92.0_f64;
    let left_margin = 86.0_f64;
    let right_margin = 36.0_f64;
    let top_margin = 64.0_f64;
    let gel_h = 420.0_f64;
    let width = (left_margin + (total_lanes as f64) * lane_w + right_margin).max(360.0);
    let height = top_margin + gel_h + 48.0;

    let (bg, gel_bg, band_color, txt_color, well_stroke) = if darkroom_theme {
        ("#090d16", "#111827", "#fde047", "#e2e8f0", "#475569")
    } else {
        ("#ffffff", "#f1f5f9", "#0f172a", "#0f172a", "#94a3b8")
    };

    let mut svg = format!(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">\n\
         <rect width=\"100%\" height=\"100%\" fill=\"{bg}\"/>\n\
         <rect x=\"{:.1}\" y=\"{top_margin:.1}\" width=\"{:.1}\" height=\"{gel_h:.1}\" rx=\"8\" fill=\"{gel_bg}\" stroke=\"{well_stroke}\" stroke-width=\"1.2\"/>\n\
         <style>\n\
           .hdr  {{ font-family: system-ui, sans-serif; font-size: 12px; font-weight: 600; fill: {txt_color}; }}\n\
           .bp   {{ font-family: ui-monospace, monospace; font-size: 10px; fill: {txt_color}; }}\n\
           .foot {{ font-family: system-ui, sans-serif; font-size: 10px; fill: #94a3b8; }}\n\
         </style>\n",
        left_margin - 12.0,
        (total_lanes as f64) * lane_w + 24.0
    );

    svg.push_str(&format!(
        "<text class=\"hdr\" x=\"{left_margin:.1}\" y=\"26\">Virtual {:.1}% Agarose Gel Electrophoresis</text>\n",
        agarose_pct
    ));

    // Lane 0: Reference Ladder
    let lad_x = left_margin + 12.0;
    let band_w = lane_w - 24.0;
    let ladder_title = match ladder {
        GelLadder::OneKbPlus => "1 kb Plus",
        GelLadder::HundredBp => "100 bp",
    };
    svg.push_str(&format!(
        "<text class=\"hdr\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{ladder_title}</text>\n\
         <rect x=\"{lad_x:.1}\" y=\"{:.1}\" width=\"{band_w:.1}\" height=\"6\" fill=\"none\" stroke=\"{well_stroke}\"/>\n",
        lad_x + band_w * 0.5,
        top_margin - 10.0,
        top_margin + 8.0
    ));

    for &(bp, bright) in ladder.rungs() {
        let frac = gel_migration_fraction(bp, agarose_pct);
        let y = top_margin + 18.0 + frac * (gel_h - 32.0);
        let h = if bright { 3.6 } else { 2.0 };
        let op = if bright { 0.96 } else { 0.62 };
        svg.push_str(&format!(
            "<rect x=\"{lad_x:.1}\" y=\"{:.2}\" width=\"{band_w:.1}\" height=\"{h:.1}\" rx=\"1\" fill=\"{band_color}\" fill-opacity=\"{op:.2}\"/>\n",
            y - h * 0.5
        ));
        if bright || bp == 10_000 || bp == 2_000 || bp == 200 || bp == 100 {
            svg.push_str(&format!(
                "<text class=\"bp\" x=\"{:.1}\" y=\"{y:.2}\" text-anchor=\"end\" dominant-baseline=\"central\">{bp} bp</text>\n",
                left_margin - 18.0
            ));
        }
    }

    // Sample Digest / PCR Lanes
    for (idx, lane) in lanes.iter().enumerate() {
        let lx = left_margin + ((idx + 1) as f64) * lane_w + 12.0;
        let max_bp = lane.fragment_lengths_bp.iter().copied().max().unwrap_or(1).max(1);

        svg.push_str(&format!(
            "<text class=\"hdr\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\">{}</text>\n\
             <rect x=\"{lx:.1}\" y=\"{:.1}\" width=\"{band_w:.1}\" height=\"6\" fill=\"none\" stroke=\"{well_stroke}\"/>\n",
            lx + band_w * 0.5,
            top_margin - 10.0,
            xml_escape(&lane.label),
            top_margin + 8.0
        ));

        for &bp in &lane.fragment_lengths_bp {
            let frac = gel_migration_fraction(bp, agarose_pct);
            let y = top_margin + 18.0 + frac * (gel_h - 32.0);
            // Stoichiometric equimolar dsDNA brightness scales with fragment mass (`bp / max_bp`)
            let rel_mass = ((bp as f64) / (max_bp as f64)).sqrt().clamp(0.35, 1.0);
            let h = 2.0 + 2.4 * rel_mass;
            let op = 0.45 + 0.50 * rel_mass;
            svg.push_str(&format!(
                "<rect x=\"{lx:.1}\" y=\"{:.2}\" width=\"{band_w:.1}\" height=\"{h:.2}\" rx=\"1.2\" fill=\"{band_color}\" fill-opacity=\"{op:.2}\"/>\n\
                 <text class=\"bp\" x=\"{:.1}\" y=\"{:.2}\" text-anchor=\"middle\">{bp}</text>\n",
                y - h * 0.5,
                lx + band_w * 0.5,
                y - h * 0.5 - 4.0
            ));
        }
    }

    if include_attribution {
        svg.push_str(&format!(
            "<text class=\"foot\" x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\">Generated by OSABT (osabt plasmid) • crvscience.com • GPL-3.0</text>\n",
            width - 14.0,
            height - 12.0
        ));
    }
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_valid_vector_map_and_agarose_gel_svgs() {
        let rec = PlasmidRecord::demo_pcrv_gfp();
        let map_svg = render_plasmid_map_svg(&rec, &[], true);
        assert!(map_svg.starts_with("<svg"));
        assert!(map_svg.contains("pCRV_DemoGFP"));
        assert!(map_svg.contains("EcoRI"));

        let gel_svg = render_virtual_gel_svg(
            GelLadder::OneKbPlus,
            &[GelLane {
                label: "EcoRI+XhoI".to_string(),
                fragment_lengths_bp: vec![742, 165],
            }],
            1.0,
            true,
            true,
        );
        assert!(gel_svg.starts_with("<svg"));
        assert!(gel_svg.contains("1 kb Plus"));
        assert!(gel_svg.contains("EcoRI+XhoI"));
    }

    #[test]
    fn circular_callout_layout_separates_hemisphere_labels_and_clears_ring() {
        let rec = PlasmidRecord::demo_pcrv_gfp();
        let n = rec.sequence.len().max(1);
        let cx = 410.0;
        let cy = 325.0;
        let r_backbone = 162.0;
        let r_feat_outer = r_backbone + 24.0;
        let mut inputs = Vec::new();
        for (i, feat) in rec.features.iter().enumerate() {
            let start_ang = (feat.start_0 as f64 / n as f64) * 2.0 * PI;
            let sweep = (feat.span_bp(n) as f64 / n as f64) * 2.0 * PI;
            inputs.push(CircularCalloutInput {
                label: feat.name.clone(),
                angle_rad: start_ang + sweep * 0.5,
                anchor_radius: r_backbone + 22.0,
                is_feature: true,
                color_hex: feat.color_hex.clone(),
                source_index: i,
            });
        }
        for (i, site) in find_unique_cutters(&rec).iter().take(16).enumerate() {
            inputs.push(CircularCalloutInput {
                label: format!("{} ({})", site.enzyme_name, site.cut_top_0 + 1),
                angle_rad: (site.cut_top_0 as f64 / n as f64) * 2.0 * PI,
                anchor_radius: r_backbone,
                is_feature: false,
                color_hex: String::from("#0284c7"),
                source_index: i,
            });
        }

        let placed = layout_circular_callouts(
            cx,
            cy,
            r_feat_outer,
            (120.0, 22.0, 700.0, 630.0),
            16.0,
            inputs,
        );
        assert!(!placed.is_empty());

        for hemisphere in [true, false] {
            let mut ys: Vec<f64> = placed
                .iter()
                .filter(|c| c.is_right_hemisphere == hemisphere)
                .map(|c| c.p2.1)
                .collect();
            ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for win in ys.windows(2) {
                assert!(
                    win[1] - win[0] >= 15.9,
                    "hemisphere={hemisphere} labels collided: {} vs {}",
                    win[0],
                    win[1]
                );
            }
        }
        for c in &placed {
            let dist = ((c.p2.0 - cx).powi(2) + (c.p2.1 - cy).powi(2)).sqrt();
            assert!(
                dist >= r_feat_outer + 15.0,
                "label '{}' intersects ring: dist={dist:.1} < {:.1}",
                c.input.label,
                r_feat_outer + 15.0
            );
        }

        // Dense synthetic clustered input near 12 o'clock and 9 o'clock
        let dense: Vec<CircularCalloutInput> = (0..22)
            .map(|i| CircularCalloutInput {
                label: format!("Site_{i}"),
                angle_rad: if i < 11 { 0.03 * (i as f64) } else { 4.5 + 0.03 * ((i - 11) as f64) },
                anchor_radius: r_backbone,
                is_feature: i % 3 == 0,
                color_hex: String::from("#0284c7"),
                source_index: i,
            })
            .collect();
        let dense_placed = layout_circular_callouts(
            cx,
            cy,
            r_feat_outer,
            (120.0, 22.0, 700.0, 630.0),
            16.0,
            dense,
        );
        for hemisphere in [true, false] {
            let mut ys: Vec<f64> = dense_placed
                .iter()
                .filter(|c| c.is_right_hemisphere == hemisphere)
                .map(|c| c.p2.1)
                .collect();
            ys.sort_by(|a, b| a.partial_cmp(b).unwrap());
            for win in ys.windows(2) {
                assert!(win[1] - win[0] >= 15.9);
            }
        }

        // Narrow-bounds ring clearance regression check
        let narrow_placed = layout_circular_callouts(
            cx,
            cy,
            r_feat_outer,
            (250.0, 22.0, 570.0, 630.0),
            16.0,
            vec![
                CircularCalloutInput {
                    label: String::from("Right3"),
                    angle_rad: PI * 0.5,
                    anchor_radius: r_backbone,
                    is_feature: true,
                    color_hex: String::from("#22c55e"),
                    source_index: 0,
                },
                CircularCalloutInput {
                    label: String::from("Left9"),
                    angle_rad: PI * 1.5,
                    anchor_radius: r_backbone,
                    is_feature: true,
                    color_hex: String::from("#38bdf8"),
                    source_index: 1,
                },
            ],
        );
        for c in &narrow_placed {
            let dist = ((c.p2.0 - cx).powi(2) + (c.p2.1 - cy).powi(2)).sqrt();
            assert!(
                dist >= r_feat_outer + 16.0,
                "narrow bounds pulled label '{}' inside ring: dist={dist:.2}",
                c.input.label
            );
        }
    }
}
