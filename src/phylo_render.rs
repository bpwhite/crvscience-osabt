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

//! Tool #4 (`osabt phylo`) — Unicode box-drawing terminal renderer (for SSH/CLI and `Ratatui`
//! interactive navigation) and publication-grade rectangular & circular vector `.svg` phylogram
//! generator.

use std::collections::HashMap;
use std::f64::consts::PI;
use std::fmt::Write as _;

use crate::phylo_core::PhyloTree;
use crate::svg_util::xml_escape;

/// A single rendered terminal line paired with the arena node ID it represents for TUI cursor tracking.
#[derive(Debug, Clone)]
pub struct RenderedTreeLine {
    pub node_id: Option<usize>,
    pub text: String,
    pub is_search_hit: bool,
}

/// Renders `tree` using crisp Unicode box-drawing characters (`┌`, `├`, `└`, `─`, `│`).
/// When `proportional_lengths` is true, horizontal stem segments scale proportionally to branch length.
pub fn render_unicode_tree(
    tree: &PhyloTree,
    proportional_lengths: bool,
    target_branch_cols: usize,
    search_query: &str,
) -> Vec<RenderedTreeLine> {
    let mut visible_tips: Vec<usize> = Vec::new();
    collect_visible_tips(tree, tree.root, &mut visible_tips);
    if visible_tips.is_empty() {
        return Vec::new();
    }

    // Assign integer vertical row coordinates: visible tips on even rows `0, 2, 4, ...`
    // and internal nodes at the midpoint `(first_child_y + last_child_y) / 2`
    let mut y_coord: HashMap<usize, usize> = HashMap::new();
    for (idx, &tip_id) in visible_tips.iter().enumerate() {
        y_coord.insert(tip_id, idx * 2);
    }
    assign_internal_y(tree, tree.root, &mut y_coord);

    // Compute horizontal column coordinates
    let mut x_coord: HashMap<usize, usize> = HashMap::new();
    let max_depth_bl = max_visible_depth(tree, tree.root, proportional_lengths).max(1e-9);
    let max_cols = target_branch_cols.clamp(18, 120);
    assign_x_coords(
        tree,
        tree.root,
        0.0,
        proportional_lengths,
        max_depth_bl,
        max_cols,
        &mut x_coord,
    );

    let grid_rows = visible_tips.len() * 2 - 1;
    let grid_cols = max_cols + 4;
    let mut grid = vec![vec![' '; grid_cols]; grid_rows];
    let mut row_node: Vec<Option<usize>> = vec![None; grid_rows];

    draw_node_to_grid(tree, tree.root, &x_coord, &y_coord, &mut grid, &mut row_node);

    let q_low = search_query.trim().to_ascii_lowercase();
    let mut lines = Vec::with_capacity(grid_rows + 2);

    for r in 0..grid_rows {
        let mut s: String = grid[r].iter().collect();
        while s.ends_with(' ') {
            s.pop();
        }
        let mut is_hit = false;

        if let Some(nid) = row_node[r] {
            let node = &tree.nodes[nid];
            if node.is_leaf() {
                let display_name = node.name.replace('_', " ");
                if !q_low.is_empty() && display_name.to_ascii_lowercase().contains(&q_low) {
                    is_hit = true;
                }
                let _ = write!(s, " {}", node.name);
                if let Some(bl) = node.branch_length {
                    let _ = write!(s, " ({bl:.4})");
                }
            } else if node.collapsed {
                let tips = tree.leaf_names_under(nid);
                let first = tips.first().map(String::as_str).unwrap_or("?");
                let last = tips.last().map(String::as_str).unwrap_or("?");
                if !q_low.is_empty()
                    && tips.iter().any(|t| t.to_ascii_lowercase().contains(&q_low))
                {
                    is_hit = true;
                }
                let _ = write!(s, "▸ [Collapsed: {} taxa | {first} .. {last}]", tips.len());
            } else if let Some(sup) = node.support_pct {
                let _ = write!(s, "[{sup:.0}]");
            }
        }

        lines.push(RenderedTreeLine {
            node_id: row_node[r],
            text: s,
            is_search_hit: is_hit,
        });
    }

    if proportional_lengths {
        let scale_unit = nice_scale_unit(max_depth_bl * 0.25);
        let bar_cols = (((scale_unit / max_depth_bl) * (max_cols as f64)).round() as usize).max(4);
        let bar = format!(
            "├{}┤ {:.3} substitutions/site",
            "─".repeat(bar_cols.saturating_sub(2)),
            scale_unit
        );
        lines.push(RenderedTreeLine {
            node_id: None,
            text: String::new(),
            is_search_hit: false,
        });
        lines.push(RenderedTreeLine {
            node_id: None,
            text: bar,
            is_search_hit: false,
        });
    }

    lines
}

fn collect_visible_tips(tree: &PhyloTree, u: usize, out: &mut Vec<usize>) {
    let n = &tree.nodes[u];
    if n.is_leaf() || n.collapsed {
        out.push(u);
    } else {
        for &ch in &n.children {
            collect_visible_tips(tree, ch, out);
        }
    }
}

fn assign_internal_y(tree: &PhyloTree, u: usize, y_map: &mut HashMap<usize, usize>) -> usize {
    if let Some(&y) = y_map.get(&u) {
        return y;
    }
    let n = &tree.nodes[u];
    let mut first_y = 0usize;
    let mut last_y = 0usize;
    for (idx, &ch) in n.children.iter().enumerate() {
        let cy = assign_internal_y(tree, ch, y_map);
        if idx == 0 {
            first_y = cy;
        }
        last_y = cy;
    }
    let mid = (first_y + last_y) / 2;
    y_map.insert(u, mid);
    mid
}

fn max_visible_depth(tree: &PhyloTree, u: usize, proportional: bool) -> f64 {
    let n = &tree.nodes[u];
    let step = if u == tree.root {
        0.0
    } else if proportional {
        n.branch_length.unwrap_or(0.05).max(0.005)
    } else {
        1.0
    };
    if n.is_leaf() || n.collapsed {
        step
    } else {
        let max_child = n
            .children
            .iter()
            .map(|&c| max_visible_depth(tree, c, proportional))
            .fold(0.0_f64, f64::max);
        step + max_child
    }
}

fn assign_x_coords(
    tree: &PhyloTree,
    u: usize,
    accum: f64,
    proportional: bool,
    max_depth: f64,
    max_cols: usize,
    x_map: &mut HashMap<usize, usize>,
) {
    let n = &tree.nodes[u];
    let next_accum = if u == tree.root {
        0.0
    } else if proportional {
        accum + n.branch_length.unwrap_or(0.05).max(0.005)
    } else {
        accum + 1.0
    };

    let mut col = ((next_accum / max_depth) * (max_cols as f64)).round() as usize;
    if let Some(p) = n.parent {
        let px = x_map.get(&p).copied().unwrap_or(0);
        col = col.max(px + 3).min(max_cols + 2);
    }
    x_map.insert(u, col);

    if !n.is_leaf() && !n.collapsed {
        for &ch in &n.children {
            assign_x_coords(tree, ch, next_accum, proportional, max_depth, max_cols, x_map);
        }
    }
}

fn draw_node_to_grid(
    tree: &PhyloTree,
    u: usize,
    x_map: &HashMap<usize, usize>,
    y_map: &HashMap<usize, usize>,
    grid: &mut [Vec<char>],
    row_node: &mut [Option<usize>],
) {
    let ux = x_map[&u];
    let uy = y_map[&u];
    if row_node[uy].is_none() || tree.nodes[u].is_leaf() || tree.nodes[u].collapsed {
        row_node[uy] = Some(u);
    }

    let node = &tree.nodes[u];
    if node.is_leaf() || node.collapsed {
        return;
    }

    let child_ys: Vec<usize> = node.children.iter().map(|c| y_map[c]).collect();
    if let (Some(&top_y), Some(&bot_y)) = (child_ys.first(), child_ys.last()) {
        for r in top_y..=bot_y {
            if ux < grid[r].len() {
                grid[r][ux] = '│';
            }
        }
        for (idx, &ch) in node.children.iter().enumerate() {
            let cy = y_map[&ch];
            let cx = x_map[&ch];
            let corner = if top_y == bot_y {
                '─'
            } else if idx == 0 {
                '┌'
            } else if idx + 1 == node.children.len() {
                '└'
            } else {
                '├'
            };
            if ux < grid[cy].len() {
                grid[cy][ux] = corner;
            }
            for c_col in (ux + 1)..cx.min(grid[cy].len()) {
                grid[cy][c_col] = '─';
            }
            draw_node_to_grid(tree, ch, x_map, y_map, grid, row_node);
        }
    }
}

/// Renders a standalone vector `.svg` phylogenetic figure in either rectangular (`circular = false`)
/// or radial polar (`circular = true`) layout.
pub fn render_phylo_svg(
    tree: &PhyloTree,
    proportional_lengths: bool,
    circular_layout: bool,
) -> String {
    if circular_layout {
        render_circular_phylo_svg(tree, proportional_lengths)
    } else {
        render_rectangular_phylo_svg(tree, proportional_lengths)
    }
}

fn render_rectangular_phylo_svg(tree: &PhyloTree, proportional: bool) -> String {
    let mut tips = Vec::new();
    collect_visible_tips(tree, tree.root, &mut tips);
    let n_tips = tips.len().max(1);

    let row_h = 28.0_f64;
    let width = 880.0_f64;
    let height = (n_tips as f64 * row_h + 104.0).max(260.0);
    let left_pad = 44.0_f64;
    let right_pad = 290.0_f64;
    let top_pad = 54.0_f64;
    let plot_w = width - left_pad - right_pad;

    let mut y_px: HashMap<usize, f64> = HashMap::new();
    for (idx, &tid) in tips.iter().enumerate() {
        y_px.insert(tid, top_pad + (idx as f64) * row_h);
    }
    assign_svg_y(tree, tree.root, &mut y_px);

    let max_d = max_visible_depth(tree, tree.root, proportional).max(1e-9);
    let mut x_px: HashMap<usize, f64> = HashMap::new();
    assign_svg_x(
        tree,
        tree.root,
        0.0,
        proportional,
        max_d,
        left_pad,
        plot_w,
        &mut x_px,
    );

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
        "  <text x=\"{left_pad:.0}\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} ({} tips)</text>",
        xml_escape(&tree.tree_name),
        tree.stats().leaf_count
    );

    emit_rectangular_edges_svg(tree, tree.root, &x_px, &y_px, &mut svg);

    if proportional {
        let unit = nice_scale_unit(max_d * 0.2);
        let bar_w = ((unit / max_d) * plot_w).clamp(24.0, plot_w * 0.5);
        let by = height - 24.0;
        let _ = writeln!(
            svg,
            "  <line x1=\"{left_pad:.1}\" y1=\"{by:.1}\" x2=\"{:.1}\" y2=\"{by:.1}\" stroke=\"#0f172a\" stroke-width=\"2\"/>",
            left_pad + bar_w
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#334155\">{unit:.3} substitutions/site</text>",
            left_pad + bar_w + 8.0,
            by + 4.0
        );
    }

    svg.push_str("</svg>\n");
    svg
}

fn assign_svg_y(tree: &PhyloTree, u: usize, y_px: &mut HashMap<usize, f64>) -> f64 {
    if let Some(&y) = y_px.get(&u) {
        return y;
    }
    let n = &tree.nodes[u];
    let first = assign_svg_y(tree, n.children[0], y_px);
    let mut last = first;
    for &ch in n.children.iter().skip(1) {
        last = assign_svg_y(tree, ch, y_px);
    }
    let mid = 0.5 * (first + last);
    y_px.insert(u, mid);
    mid
}

#[allow(clippy::too_many_arguments)]
fn assign_svg_x(
    tree: &PhyloTree,
    u: usize,
    accum: f64,
    proportional: bool,
    max_d: f64,
    left_pad: f64,
    plot_w: f64,
    x_px: &mut HashMap<usize, f64>,
) {
    let n = &tree.nodes[u];
    let next = if u == tree.root {
        0.0
    } else if proportional {
        accum + n.branch_length.unwrap_or(0.05).max(0.005)
    } else {
        accum + 1.0
    };
    x_px.insert(u, left_pad + (next / max_d) * plot_w);
    if !n.is_leaf() && !n.collapsed {
        for &ch in &n.children {
            assign_svg_x(tree, ch, next, proportional, max_d, left_pad, plot_w, x_px);
        }
    }
}

fn emit_rectangular_edges_svg(
    tree: &PhyloTree,
    u: usize,
    x_px: &HashMap<usize, f64>,
    y_px: &HashMap<usize, f64>,
    svg: &mut String,
) {
    let node = &tree.nodes[u];
    let ux = x_px[&u];
    let uy = y_px[&u];

    if node.is_leaf() {
        let clean_lbl = xml_escape(&node.name.replace('_', " "));
        let _ = writeln!(
            svg,
            "  <circle cx=\"{ux:.1}\" cy=\"{uy:.1}\" r=\"2.6\" fill=\"#0284c7\"/>"
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-style=\"italic\" font-size=\"12\" fill=\"#0f172a\">{clean_lbl}</text>",
            ux + 8.0,
            uy + 4.0
        );
        return;
    }

    if node.collapsed {
        let n_desc = tree.leaf_ids_under(u).len();
        let _ = writeln!(
            svg,
            "  <polygon points=\"{ux:.1},{uy:.1} {:.1},{:.1} {:.1},{:.1}\" fill=\"#cbd5e1\" stroke=\"#475569\" stroke-width=\"1.2\"/>",
            ux + 36.0,
            uy - 9.0,
            ux + 36.0,
            uy + 9.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"11.5\" fill=\"#334155\">[Collapsed clade: {n_desc} taxa]</text>",
            ux + 44.0,
            uy + 4.0
        );
        return;
    }

    if let (Some(first_ch), Some(last_ch)) = (node.children.first(), node.children.last()) {
        let first_cy = y_px[first_ch];
        let last_cy = y_px[last_ch];
        let _ = writeln!(
            svg,
            "  <line x1=\"{ux:.1}\" y1=\"{first_cy:.1}\" x2=\"{ux:.1}\" y2=\"{last_cy:.1}\" stroke=\"#334155\" stroke-width=\"1.6\"/>"
        );
    }

    if let Some(sup) = node.support_pct {
        let col = if sup >= 95.0 {
            "#16a34a"
        } else if sup >= 70.0 {
            "#d97706"
        } else {
            "#64748b"
        };
        let _ = writeln!(
            svg,
            "  <circle cx=\"{ux:.1}\" cy=\"{uy:.1}\" r=\"3.2\" fill=\"{col}\"/>"
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"9.5\" fill=\"{col}\">{sup:.0}</text>",
            ux - 4.0,
            uy - 4.0
        );
    }

    for &ch in &node.children {
        let cx = x_px[&ch];
        let cy = y_px[&ch];
        let _ = writeln!(
            svg,
            "  <line x1=\"{ux:.1}\" y1=\"{cy:.1}\" x2=\"{cx:.1}\" y2=\"{cy:.1}\" stroke=\"#334155\" stroke-width=\"1.6\"/>"
        );
        emit_rectangular_edges_svg(tree, ch, x_px, y_px, svg);
    }
}

fn render_circular_phylo_svg(tree: &PhyloTree, proportional: bool) -> String {
    let mut tips = Vec::new();
    collect_visible_tips(tree, tree.root, &mut tips);
    let n_tips = tips.len().max(1);

    let size = 780.0_f64;
    let cx = size * 0.5;
    let cy = size * 0.52;
    let max_r = 225.0_f64;

    let mut theta: HashMap<usize, f64> = HashMap::new();
    let span = 2.0 * PI * 0.92;
    for (idx, &tid) in tips.iter().enumerate() {
        let ang = -PI * 0.5 + (idx as f64 / (n_tips.max(2) - 1) as f64) * span;
        theta.insert(tid, ang);
    }
    assign_svg_y(tree, tree.root, &mut theta);

    let max_d = max_visible_depth(tree, tree.root, proportional).max(1e-9);
    let mut rad: HashMap<usize, f64> = HashMap::new();
    assign_svg_x(
        tree,
        tree.root,
        0.0,
        proportional,
        max_d,
        18.0,
        max_r,
        &mut rad,
    );

    let mut svg = String::with_capacity(8192);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {size:.0} {size:.0}\" width=\"{size:.0}\" height=\"{size:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{size:.0}\" height=\"{size:.0}\" fill=\"#ffffff\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"24\" y=\"30\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Circular Phylogram ({} tips)</text>",
        xml_escape(&tree.tree_name),
        tree.stats().leaf_count
    );

    emit_circular_edges_svg(tree, tree.root, cx, cy, &rad, &theta, &mut svg);
    svg.push_str("</svg>\n");
    svg
}

fn emit_circular_edges_svg(
    tree: &PhyloTree,
    u: usize,
    cx: f64,
    cy: f64,
    rad: &HashMap<usize, f64>,
    theta: &HashMap<usize, f64>,
    svg: &mut String,
) {
    let node = &tree.nodes[u];
    let ur = rad[&u];
    let ua = theta[&u];
    let ux = cx + ur * ua.cos();
    let uy = cy + ur * ua.sin();

    if node.is_leaf() || node.collapsed {
        let lbl = if node.collapsed {
            format!("[{} taxa]", tree.leaf_ids_under(u).len())
        } else {
            node.name.replace('_', " ")
        };
        let lx = cx + (ur + 8.0) * ua.cos();
        let ly = cy + (ur + 8.0) * ua.sin();
        let deg = ua * (180.0 / PI);
        let _ = writeln!(
            svg,
            "  <circle cx=\"{ux:.1}\" cy=\"{uy:.1}\" r=\"2.4\" fill=\"#0284c7\"/>"
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{lx:.1}\" y=\"{ly:.1}\" transform=\"rotate({deg:.1},{lx:.1},{ly:.1})\" font-family=\"sans-serif\" font-style=\"italic\" font-size=\"10.5\" fill=\"#0f172a\">{}</text>",
            xml_escape(&lbl)
        );
        return;
    }

    if let (Some(first_ch), Some(last_ch)) = (node.children.first(), node.children.last()) {
        let a0 = theta[first_ch];
        let a1 = theta[last_ch];
        let x0 = cx + ur * a0.cos();
        let y0 = cy + ur * a0.sin();
        let x1 = cx + ur * a1.cos();
        let y1 = cy + ur * a1.sin();
        let large_arc = usize::from((a1 - a0).abs() > PI);
        let _ = writeln!(
            svg,
            "  <path d=\"M {x0:.2} {y0:.2} A {ur:.2} {ur:.2} 0 {large_arc} 1 {x1:.2} {y1:.2}\" fill=\"none\" stroke=\"#334155\" stroke-width=\"1.4\"/>"
        );
    }

    for &ch in &node.children {
        let cr = rad[&ch];
        let ca = theta[&ch];
        let sx = cx + ur * ca.cos();
        let sy = cy + ur * ca.sin();
        let ex = cx + cr * ca.cos();
        let ey = cy + cr * ca.sin();
        let _ = writeln!(
            svg,
            "  <line x1=\"{sx:.2}\" y1=\"{sy:.2}\" x2=\"{ex:.2}\" y2=\"{ey:.2}\" stroke=\"#334155\" stroke-width=\"1.4\"/>"
        );
        emit_circular_edges_svg(tree, ch, cx, cy, rad, theta, svg);
    }
}

fn nice_scale_unit(raw: f64) -> f64 {
    if raw <= 0.0 {
        return 0.01;
    }
    let exp = raw.log10().floor();
    let base = 10_f64.powf(exp);
    let norm = raw / base;
    let step = if norm < 1.5 {
        1.0
    } else if norm < 3.5 {
        2.0
    } else if norm < 7.5 {
        5.0
    } else {
        10.0
    };
    step * base
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn renders_unicode_phylogram_and_both_svg_layouts() {
        let tree = PhyloTree::demo_vertebrate_amniotes();
        let lines = render_unicode_tree(&tree, true, 48, "Homo");
        assert!(lines.iter().any(|l| l.is_search_hit));

        let rect_svg = render_phylo_svg(&tree, true, false);
        assert!(rect_svg.starts_with("<svg"));
        assert!(rect_svg.contains("Homo sapiens"));

        let circ_svg = render_phylo_svg(&tree, true, true);
        assert!(circ_svg.starts_with("<svg"));
        assert!(circ_svg.trim_end().ends_with("</svg>"));
    }
}
