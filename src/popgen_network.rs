// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #16 (`osabt popgen`) — Median-Joining (`MJN`) haplotype network builder & `.svg` studio:
//! - **Median-Joining Network (`epsilon >= 0`)** — Bandelt, Forster & Röhl (1999)
//!   *Mol. Biol. Evol.* 16(1):37–48 (`https://pubmed.ncbi.nlm.nih.gov/10331250/`)
//! - **Fruchterman–Reingold 2D Planar Force Layout** + multi-population pie-chart `.svg` export
//!   with perpendicular branch mutation hatch marks.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::Write as _;

use crate::popgen_core::{pairwise_valid_diffs, PopSequenceRecord};

/// Haplotype or inferred Bandelt median-vector vertex in a Median-Joining network.
#[derive(Debug, Clone)]
pub struct HaplotypeNode {
    pub node_id: String,
    pub sequence: String,
    pub is_median_vector: bool,
    pub total_samples: usize,
    /// Sample counts stratified across `network.populations`.
    pub counts_by_pop: Vec<usize>,
    /// Normalized 2D layout coordinates `(x, y)` in `[0.0, 1.0]`.
    pub pos_xy: (f64, f64),
}

/// Undirected parsimony link between two haplotypes / median vectors with mutation count.
#[derive(Debug, Clone)]
pub struct HaplotypeEdge {
    pub from_idx: usize,
    pub to_idx: usize,
    pub mutations: usize,
    /// Human-readable 1-indexed nucleotide substitutions (`e.g. "12:C>T"`).
    pub subst_labels: Vec<String>,
}

/// Complete Bandelt et al. (1999) Median-Joining haplotype network with 2D layout.
#[derive(Debug, Clone)]
pub struct HaplotypeNetwork {
    pub populations: Vec<String>,
    pub epsilon: usize,
    pub sampled_haplotypes: usize,
    pub inferred_median_vectors: usize,
    pub total_network_steps: usize,
    pub nodes: Vec<HaplotypeNode>,
    pub edges: Vec<HaplotypeEdge>,
}

/// Constructs a Bandelt, Forster & Röhl (1999) Median-Joining (`MJN`) haplotype network
/// from aligned `PopSequenceRecord`s with homoplasy parameter `epsilon` (`default = 0`).
pub fn build_median_joining_network(
    records: &[PopSequenceRecord],
    epsilon: usize,
) -> Result<HaplotypeNetwork, String> {
    if records.len() < 2 {
        return Err(String::from(
            "Median-Joining network requires at least 2 aligned sequences.",
        ));
    }

    let mut pop_set: BTreeSet<String> = BTreeSet::new();
    for r in records {
        pop_set.insert(r.population.clone());
    }
    let populations: Vec<String> = pop_set.into_iter().collect();
    let pop_index: BTreeMap<&str, usize> = populations
        .iter()
        .enumerate()
        .map(|(i, p)| (p.as_str(), i))
        .collect();

    // Collapse identical sequences into sampled haplotypes H_1..H_k ordered by abundance
    let mut hap_map: BTreeMap<String, Vec<usize>> = BTreeMap::new();
    for r in records {
        let pi = pop_index[r.population.as_str()];
        hap_map
            .entry(r.sequence.clone())
            .or_insert_with(|| vec![0usize; populations.len()])[pi] += 1;
    }

    let mut raw_haps: Vec<(String, Vec<usize>, usize)> = hap_map
        .into_iter()
        .map(|(seq, cnts)| {
            let tot = cnts.iter().sum();
            (seq, cnts, tot)
        })
        .collect();
    raw_haps.sort_by(|a, b| b.2.cmp(&a.2).then_with(|| a.0.cmp(&b.0)));

    let sampled_count = raw_haps.len();
    let mut nodes: Vec<HaplotypeNode> = raw_haps
        .into_iter()
        .enumerate()
        .map(|(idx, (seq, cnts, tot))| HaplotypeNode {
            node_id: format!("H_{}", idx + 1),
            sequence: seq,
            is_median_vector: false,
            total_samples: tot,
            counts_by_pop: cnts,
            pos_xy: (0.5, 0.5),
        })
        .collect();

    // Iterative Bandelt et al. (1999) Median-Joining loop (bounded to <= 12 median vectors for parsimony)
    let mut known_seqs: BTreeSet<String> = nodes.iter().map(|n| n.sequence.clone()).collect();
    let mut mv_counter = 0usize;

    for _iter in 0..8 {
        if nodes.len() < 3 || mv_counter >= 12 {
            break;
        }
        let msn_edges = compute_epsilon_msn(&nodes, epsilon);
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); nodes.len()];
        for &(u, v, _) in &msn_edges {
            adj[u].push(v);
            adj[v].push(u);
        }

        let mut added_any = false;
        let cur_len = nodes.len();
        'triplet: for center in 0..cur_len {
            if adj[center].len() < 2 {
                continue;
            }
            for ai in 0..adj[center].len() {
                for bi in (ai + 1)..adj[center].len() {
                    let u = center;
                    let v = adj[center][ai];
                    let w = adj[center][bi];
                    let med = majority_median_sequence(
                        nodes[u].sequence.as_bytes(),
                        nodes[v].sequence.as_bytes(),
                        nodes[w].sequence.as_bytes(),
                    );
                    if known_seqs.contains(&med) {
                        continue;
                    }
                    let duv = pairwise_valid_diffs(
                        nodes[u].sequence.as_bytes(),
                        nodes[v].sequence.as_bytes(),
                    );
                    let duw = pairwise_valid_diffs(
                        nodes[u].sequence.as_bytes(),
                        nodes[w].sequence.as_bytes(),
                    );
                    let dvw = pairwise_valid_diffs(
                        nodes[v].sequence.as_bytes(),
                        nodes[w].sequence.as_bytes(),
                    );
                    // Parsimony cost of spanning {u, v, w} without median = sum of 2 shortest pairwise links
                    let mut pair_d = [duv, duw, dvw];
                    pair_d.sort_unstable();
                    let two_link_cost = pair_d[0] + pair_d[1];

                    let m_cost = pairwise_valid_diffs(med.as_bytes(), nodes[u].sequence.as_bytes())
                        + pairwise_valid_diffs(med.as_bytes(), nodes[v].sequence.as_bytes())
                        + pairwise_valid_diffs(med.as_bytes(), nodes[w].sequence.as_bytes());

                    if m_cost + epsilon < two_link_cost + epsilon && m_cost < two_link_cost {
                        mv_counter += 1;
                        known_seqs.insert(med.clone());
                        nodes.push(HaplotypeNode {
                            node_id: format!("MV_{mv_counter}"),
                            sequence: med,
                            is_median_vector: true,
                            total_samples: 0,
                            counts_by_pop: vec![0usize; populations.len()],
                            pos_xy: (0.5, 0.5),
                        });
                        added_any = true;
                        if mv_counter >= 12 {
                            break 'triplet;
                        }
                    }
                }
            }
        }
        if !added_any {
            break;
        }
    }

    // Prune obsolete median vectors with degree <= 2 in the final epsilon-MSN
    loop {
        let msn = compute_epsilon_msn(&nodes, epsilon);
        let mut deg = vec![0usize; nodes.len()];
        for &(u, v, _) in &msn {
            deg[u] += 1;
            deg[v] += 1;
        }
        if let Some(rm_idx) = (0..nodes.len())
            .rev()
            .find(|&i| nodes[i].is_median_vector && deg[i] <= 2)
        {
            nodes.remove(rm_idx);
        } else {
            break;
        }
    }

    let final_msn = compute_epsilon_msn(&nodes, epsilon);
    let mut edges = Vec::with_capacity(final_msn.len());
    let mut total_steps = 0usize;
    for (u, v, d) in final_msn {
        total_steps += d;
        let subst_labels = list_substitutions(
            nodes[u].sequence.as_bytes(),
            nodes[v].sequence.as_bytes(),
        );
        edges.push(HaplotypeEdge {
            from_idx: u,
            to_idx: v,
            mutations: d,
            subst_labels,
        });
    }

    layout_haplotype_network_2d(&mut nodes, &edges);

    let mv_count = nodes.iter().filter(|n| n.is_median_vector).count();
    Ok(HaplotypeNetwork {
        populations,
        epsilon,
        sampled_haplotypes: sampled_count,
        inferred_median_vectors: mv_count,
        total_network_steps: total_steps,
        nodes,
        edges,
    })
}

/// Computes the Bandelt/Excoffier `epsilon`-Minimum Spanning Network (`epsilon`-MSN) via Kruskal distance thresholds.
fn compute_epsilon_msn(nodes: &[HaplotypeNode], epsilon: usize) -> Vec<(usize, usize, usize)> {
    let n = nodes.len();
    if n <= 1 {
        return Vec::new();
    }
    let mut all_pairs: Vec<(usize, usize, usize)> = Vec::with_capacity(n * (n - 1) / 2);
    for i in 0..n {
        for j in (i + 1)..n {
            let d = pairwise_valid_diffs(
                nodes[i].sequence.as_bytes(),
                nodes[j].sequence.as_bytes(),
            );
            all_pairs.push((d, i, j));
        }
    }
    all_pairs.sort_unstable_by_key(|&(d, i, j)| (d, i, j));

    // Find the bottleneck MST max-edge distance first
    let mut parent: Vec<usize> = (0..n).collect();
    fn find(p: &mut [usize], mut x: usize) -> usize {
        while p[x] != x {
            p[x] = p[p[x]];
            x = p[x];
        }
        x
    }

    let mut mst_edges = Vec::new();
    let mut max_mst_d = 0usize;
    let mut comps = n;
    for &(d, u, v) in &all_pairs {
        let ru = find(&mut parent, u);
        let rv = find(&mut parent, v);
        if ru != rv {
            parent[ru] = rv;
            mst_edges.push((u, v, d));
            max_mst_d = max_mst_d.max(d);
            comps -= 1;
            if comps == 1 {
                break;
            }
        }
    }

    // Retain all MST edges PLUS any alternative parsimonious ties within `d <= min_connecting_step + epsilon`
    let mut parent2: Vec<usize> = (0..n).collect();
    let mut selected: BTreeSet<(usize, usize, usize)> = BTreeSet::new();
    let mut idx = 0usize;
    let mut comps2 = n;
    let mut threshold_after_conn: Option<usize> = None;

    while idx < all_pairs.len() {
        let cur_d = all_pairs[idx].0;
        if let Some(lim) = threshold_after_conn {
            if cur_d > lim {
                break;
            }
        }
        let start = idx;
        while idx < all_pairs.len() && all_pairs[idx].0 == cur_d {
            idx += 1;
        }
        // Check which edges at distance `cur_d` connect distinct components (or tie within epsilon)
        let mut merged_at_d = false;
        for &(_, u, v) in &all_pairs[start..idx] {
            if find(&mut parent2, u) != find(&mut parent2, v) || threshold_after_conn.is_some() {
                selected.insert((u, v, cur_d));
                merged_at_d = true;
            }
        }
        if merged_at_d {
            for &(_, u, v) in &all_pairs[start..idx] {
                let ru = find(&mut parent2, u);
                let rv = find(&mut parent2, v);
                if ru != rv {
                    parent2[ru] = rv;
                    comps2 -= 1;
                }
            }
            if comps2 == 1 && threshold_after_conn.is_none() {
                threshold_after_conn = Some(cur_d + epsilon);
            }
        }
    }

    for e in mst_edges {
        selected.insert(e);
    }
    selected.into_iter().collect()
}

fn majority_median_sequence(a: &[u8], b: &[u8], c: &[u8]) -> String {
    let mut out = String::with_capacity(a.len());
    for i in 0..a.len() {
        let ch = if a[i] == b[i] || a[i] == c[i] {
            a[i]
        } else if b[i] == c[i] {
            b[i]
        } else {
            a[i]
        };
        out.push(ch as char);
    }
    out
}

fn list_substitutions(a: &[u8], b: &[u8]) -> Vec<String> {
    let mut out = Vec::new();
    for (idx, (&x, &y)) in a.iter().zip(b).enumerate() {
        if x != y {
            out.push(format!("{}:{}>{}", idx + 1, x as char, y as char));
        }
    }
    out
}

fn layout_haplotype_network_2d(nodes: &mut [HaplotypeNode], edges: &[HaplotypeEdge]) {
    let n = nodes.len();
    if n == 0 {
        return;
    }
    if n == 1 {
        nodes[0].pos_xy = (0.5, 0.5);
        return;
    }

    // Deterministic circular initialization followed by Fruchterman–Reingold spring-electrical relaxation
    let mut pos: Vec<(f64, f64)> = (0..n)
        .map(|i| {
            let ang = 2.0 * std::f64::consts::PI * (i as f64) / (n as f64);
            (0.5 + 0.35 * ang.cos(), 0.5 + 0.35 * ang.sin())
        })
        .collect();

    let k_ideal = (0.65 / (n as f64).sqrt()).max(0.10);
    for step in 0..220 {
        let temp = 0.12 * (1.0 - (step as f64) / 220.0).max(0.05);
        let mut disp = vec![(0.0_f64, 0.0_f64); n];

        for i in 0..n {
            for j in (i + 1)..n {
                let dx = pos[i].0 - pos[j].0;
                let dy = pos[i].1 - pos[j].1;
                let dist = (dx * dx + dy * dy).sqrt().max(0.01);
                let rep = (k_ideal * k_ideal) / dist;
                let rx = (dx / dist) * rep;
                let ry = (dy / dist) * rep;
                disp[i].0 += rx;
                disp[i].1 += ry;
                disp[j].0 -= rx;
                disp[j].1 -= ry;
            }
        }

        for e in edges {
            let u = e.from_idx;
            let v = e.to_idx;
            let dx = pos[u].0 - pos[v].0;
            let dy = pos[u].1 - pos[v].1;
            let dist = (dx * dx + dy * dy).sqrt().max(0.005);
            let target_len = k_ideal * (0.75 + 0.25 * (e.mutations.min(6) as f64));
            let attr = (dist - target_len) * 0.45;
            let ax = (dx / dist) * attr;
            let ay = (dy / dist) * attr;
            disp[u].0 -= ax;
            disp[u].1 -= ay;
            disp[v].0 += ax;
            disp[v].1 += ay;
        }

        for i in 0..n {
            let mag = (disp[i].0 * disp[i].0 + disp[i].1 * disp[i].1).sqrt().max(1e-9);
            let step_len = mag.min(temp);
            pos[i].0 += (disp[i].0 / mag) * step_len;
            pos[i].1 += (disp[i].1 / mag) * step_len;
        }
    }

    let min_x = pos.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
    let max_x = pos.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
    let min_y = pos.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
    let max_y = pos.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
    let span_x = (max_x - min_x).max(1e-6);
    let span_y = (max_y - min_y).max(1e-6);

    for (node, &(px, py)) in nodes.iter_mut().zip(&pos) {
        node.pos_xy = (
            0.12 + 0.76 * ((px - min_x) / span_x),
            0.14 + 0.72 * ((py - min_y) / span_y),
        );
    }
}

/// Renders a publication-ready `.svg` Median-Joining haplotype network with multi-population
/// pie-chart nodes, black diamonds for inferred median vectors, and mutation hatch ticks.
pub fn render_popgen_network_svg(network: &HaplotypeNetwork, title: &str) -> String {
    let width = 900.0_f64;
    let height = 560.0_f64;
    const PALETTE: &[&str] = &[
        "#2563eb", "#16a34a", "#dc2626", "#9333ea", "#d97706", "#0891b2",
    ];

    let map_xy = |(nx, ny): (f64, f64)| (60.0 + nx * 640.0, 65.0 + ny * 440.0);

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
        "  <text x=\"42\" y=\"32\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — Bandelt et al. (1999) Median-Joining Network ({} Haplotypes, {} MVs, {} Steps)</text>",
        xml_escape(title),
        network.sampled_haplotypes,
        network.inferred_median_vectors,
        network.total_network_steps
    );

    // Draw edges + perpendicular mutation hatch marks
    for e in &network.edges {
        let (x1, y1) = map_xy(network.nodes[e.from_idx].pos_xy);
        let (x2, y2) = map_xy(network.nodes[e.to_idx].pos_xy);
        let _ = writeln!(
            svg,
            "  <line x1=\"{x1:.1}\" y1=\"{y1:.1}\" x2=\"{x2:.1}\" y2=\"{y2:.1}\" stroke=\"#475569\" stroke-width=\"2.0\"/>"
        );
        let dx = x2 - x1;
        let dy = y2 - y1;
        let len = (dx * dx + dy * dy).sqrt().max(1.0);
        let nx = -dy / len;
        let ny = dx / len;
        let ticks = e.mutations.min(12);
        for t in 1..=ticks {
            let frac = 0.35 + 0.30 * ((t as f64) / ((ticks + 1) as f64));
            let tx = x1 + frac * dx;
            let ty = y1 + frac * dy;
            let _ = writeln!(
                svg,
                "  <line x1=\"{:.1}\" y1=\"{:.1}\" x2=\"{:.1}\" y2=\"{:.1}\" stroke=\"#dc2626\" stroke-width=\"1.8\"/>",
                tx - 5.0 * nx,
                ty - 5.0 * ny,
                tx + 5.0 * nx,
                ty + 5.0 * ny
            );
        }
    }

    // Draw haplotype pie-chart nodes & median-vector diamonds
    for node in &network.nodes {
        let (cx, cy) = map_xy(node.pos_xy);
        if node.is_median_vector {
            let _ = writeln!(
                svg,
                "  <polygon points=\"{cx:.1},{:.1} {:.1},{cy:.1} {cx:.1},{:.1} {:.1},{cy:.1}\" fill=\"#0f172a\" stroke=\"#ffffff\" stroke-width=\"1.5\"/>",
                cy - 7.0,
                cx + 7.0,
                cy + 7.0,
                cx - 7.0
            );
            let _ = writeln!(
                svg,
                "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"10\" fill=\"#334155\">{}</text>",
                cy - 11.0,
                node.node_id
            );
        } else {
            let rad = 9.0 + 4.5 * (node.total_samples.max(1) as f64).sqrt();
            let active_pops: Vec<(usize, usize)> = node
                .counts_by_pop
                .iter()
                .copied()
                .enumerate()
                .filter(|&(_, c)| c > 0)
                .collect();

            if active_pops.len() <= 1 {
                let p_idx = active_pops.first().map(|&(i, _)| i).unwrap_or(0);
                let col = PALETTE[p_idx % PALETTE.len()];
                let _ = writeln!(
                    svg,
                    "  <circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{rad:.1}\" fill=\"{col}\" stroke=\"#0f172a\" stroke-width=\"1.6\"/>"
                );
            } else {
                let tot = node.total_samples as f64;
                let mut cur_ang = -0.5 * std::f64::consts::PI;
                for (p_idx, cnt) in active_pops {
                    let sweep = 2.0 * std::f64::consts::PI * ((cnt as f64) / tot);
                    let end_ang = cur_ang + sweep;
                    let sx = cx + rad * cur_ang.cos();
                    let sy = cy + rad * cur_ang.sin();
                    let ex = cx + rad * end_ang.cos();
                    let ey = cy + rad * end_ang.sin();
                    let large_arc = if sweep > std::f64::consts::PI { 1 } else { 0 };
                    let col = PALETTE[p_idx % PALETTE.len()];
                    let _ = writeln!(
                        svg,
                        "  <path d=\"M {cx:.2} {cy:.2} L {sx:.2} {sy:.2} A {rad:.2} {rad:.2} 0 {large_arc} 1 {ex:.2} {ey:.2} Z\" fill=\"{col}\" stroke=\"#0f172a\" stroke-width=\"1.2\"/>"
                    );
                    cur_ang = end_ang;
                }
            }
            let _ = writeln!(
                svg,
                "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11\" font-weight=\"700\" fill=\"#0f172a\">{} (n={})</text>",
                cy + rad + 14.0,
                node.node_id,
                node.total_samples
            );
        }
    }

    // Legend box on right side
    let lx = 730.0_f64;
    let mut ly = 85.0_f64;
    let _ = writeln!(
        svg,
        "  <text x=\"{lx:.0}\" y=\"{ly:.0}\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"#0f172a\">Populations</text>"
    );
    for (idx, pop) in network.populations.iter().enumerate() {
        ly += 24.0;
        let col = PALETTE[idx % PALETTE.len()];
        let _ = writeln!(
            svg,
            "  <rect x=\"{lx:.0}\" y=\"{:.0}\" width=\"13\" height=\"13\" rx=\"2\" fill=\"{col}\"/>",
            ly - 11.0
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.0}\" y=\"{ly:.0}\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#1e293b\">{}</text>",
            lx + 20.0,
            xml_escape(pop)
        );
    }
    ly += 28.0;
    let _ = writeln!(
        svg,
        "  <polygon points=\"{:.1},{:.1} {:.1},{:.1} {:.1},{:.1} {:.1},{:.1}\" fill=\"#0f172a\"/>",
        lx + 6.5,
        ly - 11.0,
        lx + 13.0,
        ly - 4.5,
        lx + 6.5,
        ly + 2.0,
        lx,
        ly - 4.5
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{:.0}\" y=\"{ly:.0}\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#1e293b\">Median Vector (MV)</text>",
        lx + 20.0
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

    #[test]
    fn median_joining_infers_missing_steiner_node_and_renders_svg() {
        let records = vec![
            PopSequenceRecord {
                sample_id: "s1".into(),
                population: "North".into(),
                sequence: "GAAAAA".into(),
            },
            PopSequenceRecord {
                sample_id: "s2".into(),
                population: "North".into(),
                sequence: "GAAAAA".into(),
            },
            PopSequenceRecord {
                sample_id: "s3".into(),
                population: "South".into(),
                sequence: "ACAAAA".into(),
            },
            PopSequenceRecord {
                sample_id: "s4".into(),
                population: "West".into(),
                sequence: "AAATAA".into(),
            },
        ];

        let net = build_median_joining_network(&records, 0).expect("MJN succeeds");
        assert_eq!(net.sampled_haplotypes, 3);
        assert!(
            net.inferred_median_vectors >= 1 && net.nodes.iter().any(|n| n.is_median_vector),
            "expected inferred median vector for 3-way star split"
        );
        let observed_total: usize = net
            .nodes
            .iter()
            .filter(|n| !n.is_median_vector)
            .map(|n| n.total_samples)
            .sum();
        assert_eq!(observed_total, records.len());
        assert!(net
            .edges
            .iter()
            .all(|e| e.from_idx < net.nodes.len() && e.to_idx < net.nodes.len() && e.mutations > 0));

        let svg = render_popgen_network_svg(&net, "Steiner <Median> Network");
        assert!(svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>"));
        assert!(svg.contains("Steiner &lt;Median&gt; Network"));
    }
}
