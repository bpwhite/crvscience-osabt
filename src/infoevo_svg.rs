// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #11 (`osabt infoevo`) — Publication-grade vector `.svg` visualizers for:
//! 1. **Schneider & Stephens (1990) Sequence Logos** (`0..2` bits DNA/RNA or `0..4.32` bits Protein
//!    with vertically scaled glyph stacks proportional to `f_b(l) * R_seq(l)`).
//! 2. **Dunn, Wahl & Gloor (2008) `MIp` Epistasis Arc & Contact Map** overlaying pairwise
//!    coevolving interactions onto positional conservation `R_seq(l)`.
//! 3. **Kimura (1961) Adaptive Information Gain vs. Substitutional Load** (`I_t` vs. `L_sub,t`).

use std::fmt::Write as _;

use crate::infoevo_core::{
    AlignmentComplexityReport, EpistasisMatrixReport, KimuraSelectionReport,
};

/// Standard nucleotide & amino-acid biochemistry color palette for sequence logos.
pub fn logo_symbol_color(ch: char) -> &'static str {
    match ch.to_ascii_uppercase() {
        'A' => "#16a34a", // Adenine / Alanine — green
        'C' => "#2563eb", // Cytosine / Cysteine — blue
        'G' => "#d97706", // Guanine / Glycine — amber/orange
        'T' | 'U' => "#dc2626", // Thymine / Uracil — red
        'D' | 'E' => "#e11d48", // Acidic — crimson
        'K' | 'R' | 'H' => "#1d4ed8", // Basic — royal blue
        'F' | 'W' | 'Y' => "#7c3aed", // Aromatic — violet
        'N' | 'Q' | 'S' => "#0d9488", // Polar — teal
        _ => "#475569",   // Hydrophobic / other — slate
    }
}

/// Renders a vector `.svg` **Schneider & Stephens (1990) Sequence Logo** using vertical `<g transform>`
/// scaling so every residue glyph exactly spans its information height `h_b(l) = f_b(l) * R_seq(l)` in bits.
pub fn render_sequence_logo_svg(report: &AlignmentComplexityReport, title: &str) -> String {
    let n_pos = report.alignment_length.max(1);
    let col_w = (720.0 / (n_pos as f64)).clamp(16.0, 38.0);
    let left = 72.0_f64;
    let top = 58.0_f64;
    let plot_w = col_w * (n_pos as f64);
    let plot_h = 250.0_f64;
    let width = left + plot_w + 48.0;
    let height = top + plot_h + 76.0;
    let max_bits = report.max_bits_per_site.max(2.0);

    let mut svg = String::with_capacity(16_384);
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
        "  <text x=\"{left:.0}\" y=\"26\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Sequence Logo (R_seq = {:.2} bits | Adami C = {:.2} bits)</text>",
        xml_escape(title),
        report.total_r_sequence_bits,
        report.adami_physical_complexity_bits
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"44\" font-family=\"sans-serif\" font-size=\"11.5\" fill=\"#475569\">{} | n = {} sequences | L = {} columns | Complexity density = {:.1}%</text>",
        report.alphabet.label(),
        report.num_sequences,
        report.alignment_length,
        report.complexity_density * 100.0
    );

    // Y-axis gridlines in bits
    let y_base = top + plot_h;
    let max_int_tick = max_bits.floor() as usize;
    for b_tick in 0..=max_int_tick {
        let b_val = b_tick as f64;
        let py = y_base - (b_val / max_bits) * plot_h;
        let _ = writeln!(
            svg,
            "  <line x1=\"{left:.1}\" y1=\"{py:.1}\" x2=\"{:.1}\" y2=\"{py:.1}\" stroke=\"#e2e8f0\" stroke-width=\"1\"/>",
            left + plot_w
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"11\" fill=\"#334155\">{b_tick}.0</text>",
            left - 8.0,
            py + 4.0
        );
    }
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{top:.1}\" x2=\"{left:.1}\" y2=\"{y_base:.1}\" stroke=\"#334155\" stroke-width=\"1.4\"/>"
    );
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{y_base:.1}\" x2=\"{:.1}\" y2=\"{y_base:.1}\" stroke=\"#334155\" stroke-width=\"1.4\"/>",
        left + plot_w
    );
    let _ = writeln!(
        svg,
        "  <text x=\"20\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"600\" fill=\"#1e293b\" transform=\"rotate(-90 20 {:.1})\">Information (bits)</text>",
        top + 0.5 * plot_h,
        top + 0.5 * plot_h
    );

    // Stacked glyphs per column
    for (idx, site) in report.sites.iter().enumerate() {
        let cx = left + ((idx as f64) + 0.5) * col_w;
        let mut cur_y_bottom = y_base;

        for &(sym, h_bits) in &site.letter_heights_bits {
            let box_h = (h_bits / max_bits) * plot_h;
            if box_h < 0.8 {
                cur_y_bottom -= box_h;
                continue;
            }
            // Reference capital glyph cap-height in 100-unit font is ~72 units
            let scale_y = box_h / 72.0;
            let scale_x = (col_w * 0.88) / 62.0;
            let col = logo_symbol_color(sym);
            let _ = writeln!(
                svg,
                "  <g transform=\"translate({cx:.2},{cur_y_bottom:.2}) scale({scale_x:.3},{scale_y:.3})\"><text x=\"0\" y=\"0\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"100\" font-weight=\"800\" fill=\"{col}\">{sym}</text></g>"
            );
            cur_y_bottom -= box_h;
        }

        if n_pos <= 40 || (site.position_1based % 5 == 0) || site.position_1based == 1 {
            let _ = writeln!(
                svg,
                "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10\" fill=\"#475569\">{}</text>",
                y_base + 16.0,
                site.position_1based
            );
        }
    }

    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#334155\">Alignment Position (5' → 3' / N → C)</text>",
        left + 0.5 * plot_w,
        height - 18.0
    );

    svg.push_str("</svg>\n");
    svg
}

/// Renders an `.svg` **Epistasis / Coevolution Arc & Conservation Diagram** showing top Dunn et al. (2008)
/// `MIp` covarying pairs as weighted Bézier arcs arching over the positional `R_seq(l)` conservation bar track.
pub fn render_epistasis_mip_svg(
    epistasis: &EpistasisMatrixReport,
    info: &AlignmentComplexityReport,
    title: &str,
) -> String {
    let n_pos = info.alignment_length.max(1);
    let col_w = (720.0 / (n_pos as f64)).clamp(16.0, 34.0);
    let left = 68.0_f64;
    let plot_w = col_w * (n_pos as f64);
    let arc_base_y = 240.0_f64;
    let bar_h_max = 95.0_f64;
    let width = left + plot_w + 48.0;
    let height = 400.0_f64;

    let max_mip = epistasis
        .top_pairs
        .first()
        .map(|p| p.mip_bits.max(0.05))
        .unwrap_or(1.0);

    let mut svg = String::with_capacity(12_288);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#0f172a\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#f8fafc\">{} — APC-Corrected Epistasis Arc Map (Dunn et al. 2008 MIp)</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"46\" font-family=\"sans-serif\" font-size=\"11.5\" fill=\"#94a3b8\">Top arcs = compensatory coevolution (MIp Z &gt; 2.0)  |  Lower bars = single-site conservation R_seq (bits)</text>"
    );

    // Draw coevolving arcs above `arc_base_y`
    for pair in epistasis.top_pairs.iter().filter(|p| p.mip_bits > 0.02) {
        let x1 = left + ((pair.pos_i_1based - 1) as f64 + 0.5) * col_w;
        let x2 = left + ((pair.pos_j_1based - 1) as f64 + 0.5) * col_w;
        let mid_x = 0.5 * (x1 + x2);
        let span_ratio = ((pair.pos_j_1based - pair.pos_i_1based) as f64) / (n_pos as f64);
        let ctrl_y = (arc_base_y - (45.0 + 215.0 * span_ratio)).max(52.0);
        let rel = (pair.mip_bits / max_mip).clamp(0.15, 1.0);
        let stroke_w = 1.2 + 3.4 * rel;
        let col = if pair.z_score >= 3.0 {
            "#38bdf8"
        } else {
            "#a78bfa"
        };
        let _ = writeln!(
            svg,
            "  <path d=\"M {x1:.1} {arc_base_y:.1} Q {mid_x:.1} {ctrl_y:.1} {x2:.1} {arc_base_y:.1}\" fill=\"none\" stroke=\"{col}\" stroke-opacity=\"{:.2}\" stroke-width=\"{stroke_w:.2}\"/>",
            0.35 + 0.60 * rel
        );
    }

    // Draw R_seq bars downward from `arc_base_y + 16`
    let bar_top = arc_base_y + 18.0;
    let max_b = info.max_bits_per_site.max(2.0);
    for (idx, s) in info.sites.iter().enumerate() {
        let x = left + (idx as f64 + 0.12) * col_w;
        let bw = col_w * 0.76;
        let bh = ((s.r_sequence_bits / max_b).clamp(0.0, 1.0)) * bar_h_max;
        let _ = writeln!(
            svg,
            "  <rect x=\"{x:.1}\" y=\"{bar_top:.1}\" width=\"{bw:.1}\" height=\"{bh:.1}\" fill=\"#22c55e\" rx=\"2\"/>"
        );
        let cx = left + (idx as f64 + 0.5) * col_w;
        let _ = writeln!(
            svg,
            "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"10\" fill=\"#e2e8f0\">{}</text>",
            arc_base_y + 12.0,
            s.consensus_symbol
        );
        if n_pos <= 36 || s.position_1based % 5 == 0 || s.position_1based == 1 {
            let _ = writeln!(
                svg,
                "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"9.5\" fill=\"#94a3b8\">{}</text>",
                bar_top + bar_h_max + 15.0,
                s.position_1based
            );
        }
    }

    svg.push_str("</svg>\n");
    svg
}

/// Renders an `.svg` plot of Kimura (1961) accumulated genetic information `I_t = D_KL(p_t || p_0)`
/// alongside cumulative substitutional load `L_sub,t` and the asymptotic `-log2(p_0^*)` bound.
pub fn render_kimura_trajectory_svg(sim: &KimuraSelectionReport, title: &str) -> String {
    let width = 820.0_f64;
    let height = 430.0_f64;
    let left = 72.0_f64;
    let top = 56.0_f64;
    let pw = 540.0_f64;
    let ph = 315.0_f64;

    let y_max = (sim
        .final_substitutional_load_bits
        .max(sim.kimura_limit_bits)
        * 1.12)
        .max(1.0);
    let g_max = (sim.generations.max(1)) as f64;

    let map_x = |g: usize| left + ((g as f64) / g_max) * pw;
    let map_y = |b: f64| top + ph - (b / y_max).clamp(0.0, 1.0) * ph;

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
        "  <text x=\"{left:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Kimura (1961) Information Gain vs. Substitutional Load</text>",
        xml_escape(title)
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );

    let y_lim = map_y(sim.kimura_limit_bits);
    let _ = writeln!(
        svg,
        "  <line x1=\"{left:.1}\" y1=\"{y_lim:.1}\" x2=\"{:.1}\" y2=\"{y_lim:.1}\" stroke=\"#64748b\" stroke-width=\"1.4\" stroke-dasharray=\"6,4\"/>",
        left + pw
    );

    let mut poly_info = String::new();
    let mut poly_load = String::new();
    let stride = (sim.steps.len() / 160).max(1);
    for (idx, st) in sim.steps.iter().enumerate() {
        if idx % stride == 0 || idx + 1 == sim.steps.len() {
            let x = map_x(st.generation);
            let _ = write!(poly_info, "{x:.2},{:.2} ", map_y(st.kl_info_gain_bits));
            let _ = write!(
                poly_load,
                "{x:.2},{:.2} ",
                map_y(st.cumulative_substitutional_load_bits)
            );
        }
    }

    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#dc2626\" stroke-width=\"2.4\"/>",
        poly_load.trim()
    );
    let _ = writeln!(
        svg,
        "  <polyline points=\"{}\" fill=\"none\" stroke=\"#2563eb\" stroke-width=\"2.6\"/>",
        poly_info.trim()
    );

    let lx = left + pw + 18.0;
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"#2563eb\">KL Info Gain I_T: {:.3} bits</text>",
        top + 28.0,
        sim.final_info_gain_bits
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"#dc2626\">Subst. Load L_sub: {:.3} bits</text>",
        top + 52.0,
        sim.final_substitutional_load_bits
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#475569\">Kimura Limit: {:.3} bits</text>",
        top + 76.0,
        sim.kimura_limit_bits
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"11.5\" fill=\"#0f172a\">Efficiency: {:.1}%</text>",
        top + 100.0,
        sim.final_efficiency * 100.0
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

#[cfg(test)]
mod tests {
    use super::*;
    use crate::infoevo_core::{
        compute_alignment_information, compute_epistasis_mip_matrix,
        demo_ecoli_promoter_and_stemloop_fasta, parse_fasta_alignment,
        simulate_kimura_selection_information,
    };

    #[test]
    fn renders_valid_infoevo_logo_epistasis_and_kimura_svgs_with_escaping() {
        let aln = parse_fasta_alignment("demo", demo_ecoli_promoter_and_stemloop_fasta())
            .expect("demo FASTA parses");
        let info = compute_alignment_information(&aln, true, None);
        let mip = compute_epistasis_mip_matrix(&aln, 0.25, 8);
        let kimura = simulate_kimura_selection_information(&[0.5, 0.3, 0.2], &[0.0, 0.4, 0.8], 1.0, 20)
            .expect("Kimura trajectory succeeds");

        let logo_svg = render_sequence_logo_svg(&info, "Logo <Promoter>");
        let mip_svg = render_epistasis_mip_svg(&mip, &info, "MIp <Matrix>");
        let kim_svg = render_kimura_trajectory_svg(&kimura, "Kimura <Gain & Load>");

        for svg in [&logo_svg, &mip_svg, &kim_svg] {
            assert!(svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>"));
        }
        assert!(logo_svg.contains("Logo &lt;Promoter&gt;"));
        assert!(mip_svg.contains("MIp &lt;Matrix&gt;"));
        assert!(kim_svg.contains("Kimura &lt;Gain &amp; Load&gt;"));
    }
}
