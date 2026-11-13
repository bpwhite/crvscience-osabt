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

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
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
    let width = 760.0_f64;
    let height = 640.0_f64;
    let cx = 380.0_f64;
    let cy = 310.0_f64;
    let r_backbone = 175.0_f64;
    let r_gc_base = 128.0_f64;
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
        let r = r_gc_base + (local_gc - mean_gc) * 42.0;
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

    // 4. Annotated directional feature arcs + radial labels
    for feat in &rec.features {
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
        let mid_ang = start_ang + sweep * 0.5;
        let (lx, ly) = polar_xy(cx, cy, r_backbone + 38.0, mid_ang);
        let anchor = if lx >= cx { "start" } else { "end" };
        svg.push_str(&format!(
            "<path d=\"{d}\" fill=\"{}\" fill-opacity=\"0.88\" stroke=\"#0f172a\" stroke-width=\"0.8\"/>\n\
             <text class=\"feat\" x=\"{lx:.1}\" y=\"{ly:.1}\" text-anchor=\"{anchor}\" dominant-baseline=\"central\">{}</text>\n",
            xml_escape(&feat.color_hex),
            xml_escape(&feat.name)
        ));
    }

    // 5. Unique restriction cutter callouts
    let uniques = find_unique_cutters(rec);
    for (u_idx, site) in uniques.iter().take(16).enumerate() {
        let ang = (site.cut_top_0 as f64 / n as f64) * 2.0 * PI;
        let r_tick_out = r_backbone + 56.0 + ((u_idx % 2) as f64) * 14.0;
        let (x1, y1) = polar_xy(cx, cy, r_backbone, ang);
        let (x2, y2) = polar_xy(cx, cy, r_tick_out, ang);
        let anchor = if x2 >= cx { "start" } else { "end" };
        let dx = if x2 >= cx { 4.0 } else { -4.0 };
        svg.push_str(&format!(
            "<line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#0284c7\" stroke-width=\"1.0\" stroke-dasharray=\"2 2\"/>\n\
             <text class=\"enz\" x=\"{:.1}\" y=\"{y2:.1}\" text-anchor=\"{anchor}\" dominant-baseline=\"central\">{} ({})</text>\n",
            x2 + dx,
            xml_escape(&site.enzyme_name),
            site.cut_top_0 + 1
        ));
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
            top_margin + 8.0,
            xml_escape(&lane.label)
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
}
