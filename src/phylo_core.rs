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

//! Tool #4 (`osabt phylo`) — Arena-allocated phylogenetic tree engine supporting
//! Newick (`.nwk`), NHX (`[&&NHX:...]`), and Maddison et al. (1997) NEXUS (`BEGIN TREES;` +
//! `TRANSLATE`) I/O, outgroup & midpoint rerooting, ladderizing, clade collapsing, MRCA
//! monophyly testing, and Robinson & Foulds (1981) symmetric bipartition distances.

use std::collections::{BTreeMap, BTreeSet, HashMap};
use std::fmt::Write as _;

/// Single node inside an arena-allocated [`PhyloTree`].
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PhyloNode {
    pub id: usize,
    pub parent: Option<usize>,
    pub children: Vec<usize>,
    /// Taxon name (for leaves) or clade/internal label.
    pub name: String,
    /// Branch length leading from `parent` to this node (substitutions/site or time).
    pub branch_length: Option<f64>,
    /// Bootstrap, UFBoot, or posterior support (`0.0..=100.0` or `0.0..=1.0` normalized to `%`).
    pub support_pct: Option<f64>,
    /// Optional NHX / BEAST / MrBayes bracketed key-value annotations.
    pub annotations: BTreeMap<String, String>,
    /// Interactive viewer state: when true on an internal node, descendants are folded into a wedge.
    pub collapsed: bool,
}

impl PhyloNode {
    pub fn is_leaf(&self) -> bool {
        self.children.is_empty()
    }
}

/// Arena-based rooted or unrooted phylogenetic tree.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PhyloTree {
    pub tree_name: String,
    pub nodes: Vec<PhyloNode>,
    pub root: usize,
}

/// Summary tree metrics for CLI and TUI headers.
#[derive(Debug, Clone, Default)]
pub struct PhyloTreeStats {
    pub leaf_count: usize,
    pub internal_count: usize,
    pub polytomy_count: usize,
    pub collapsed_clades: usize,
    pub total_branch_length: f64,
    pub max_root_to_tip: f64,
    pub min_root_to_tip: f64,
    pub mean_support_pct: Option<f64>,
}

/// Classification of a target taxon set under MRCA monophyly testing.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum MonophylyStatus {
    Monophyletic,
    Paraphyletic,
    Polyphyletic,
}

/// Detailed result of testing whether a set of query taxa forms a monophyletic clade.
#[derive(Debug, Clone)]
pub struct MonophylyReport {
    pub status: MonophylyStatus,
    pub mrca_node_id: usize,
    pub matched_taxa: Vec<String>,
    pub missing_query_taxa: Vec<String>,
    /// Non-query taxa nested inside the MRCA clade (empty iff `Monophyletic`).
    pub intruder_taxa: Vec<String>,
    pub mrca_support_pct: Option<f64>,
    pub mrca_stem_length: Option<f64>,
}

/// Result of comparing two phylogenetic trees via Robinson & Foulds (1981) bipartition splits.
#[derive(Debug, Clone)]
pub struct RobinsonFouldsReport {
    pub shared_taxa_count: usize,
    pub splits_tree_a: usize,
    pub splits_tree_b: usize,
    pub shared_splits: usize,
    /// Unweighted symmetric difference `|A \ B| + |B \ A|`.
    pub rf_distance: usize,
    /// Maximum possible RF distance `splits_tree_a + splits_tree_b`.
    pub max_rf: usize,
    /// Normalized RF distance in `0.0..=1.0`.
    pub normalized_rf: f64,
    /// Kuhner–Felsenstein branch-length L1 difference across non-trivial bipartitions.
    pub branch_score_l1: f64,
    pub unique_to_a: Vec<Vec<String>>,
    pub unique_to_b: Vec<Vec<String>>,
}

impl PhyloTree {
    /// Creates an empty single-root tree.
    pub fn new(tree_name: impl Into<String>) -> Self {
        let root_node = PhyloNode {
            id: 0,
            parent: None,
            children: Vec::new(),
            name: String::new(),
            branch_length: None,
            support_pct: None,
            annotations: BTreeMap::new(),
            collapsed: false,
        };
        Self {
            tree_name: tree_name.into(),
            nodes: vec![root_node],
            root: 0,
        }
    }

    /// Built-in 14-taxon Vertebrate + Tunicate outgroup maximum-likelihood phylogram
    /// with realistic substitutions/site branch lengths and UFBoot support percentages.
    pub fn demo_vertebrate_amniotes() -> Self {
        let nwk = "(\
            Ciona_intestinalis:0.380,\
            (\
                (\
                    Danio_rerio:0.185,\
                    (Oncorhynchus_mykiss:0.092,Oncorhynchus_tshawytscha:0.088)99:0.074\
                )96:0.115,\
                (\
                    (Xenopus_tropicalis:0.148,Rana_cascadae:0.132)95:0.089,\
                    (\
                        (\
                            (Gallus_gallus:0.078,Taeniopygia_guttata:0.084)100:0.062,\
                            Alligator_mississippiensis:0.091\
                        )98:0.054,\
                        (\
                            Ornithorhynchus_anatinus:0.126,\
                            (\
                                (Mus_musculus:0.089,Rattus_norvegicus:0.082)100:0.071,\
                                (Homo_sapiens:0.031,Pan_troglodytes:0.034)100:0.064\
                            )97:0.058\
                        )94:0.049\
                    )100:0.083\
                )93:0.076\
            )100:0.160\
        )Vertebrata_Root;";
        let mut t = parse_newick("Vertebrata_14Taxa_ML", nwk).expect("valid built-in Newick");
        t.ladderize(true);
        t
    }

    /// Returns all reachable node IDs in preorder (`node` before `children`).
    pub fn preorder_ids(&self, start: usize) -> Vec<usize> {
        let mut out = Vec::with_capacity(self.nodes.len());
        let mut stack = vec![start];
        while let Some(u) = stack.pop() {
            out.push(u);
            for &ch in self.nodes[u].children.iter().rev() {
                stack.push(ch);
            }
        }
        out
    }

    /// Returns all reachable node IDs in postorder (`children` before `node`).
    pub fn postorder_ids(&self, start: usize) -> Vec<usize> {
        let mut pre = self.preorder_ids(start);
        pre.reverse();
        pre
    }

    /// Collects all descendant leaf IDs beneath `node_id` (ignoring `collapsed` flag).
    pub fn leaf_ids_under(&self, node_id: usize) -> Vec<usize> {
        let mut leaves = Vec::new();
        for u in self.preorder_ids(node_id) {
            if self.nodes[u].is_leaf() {
                leaves.push(u);
            }
        }
        leaves
    }

    /// Collects all descendant leaf names beneath `node_id`.
    pub fn leaf_names_under(&self, node_id: usize) -> Vec<String> {
        self.leaf_ids_under(node_id)
            .into_iter()
            .map(|i| self.nodes[i].name.clone())
            .collect()
    }

    /// Finds the first node (preferring leaves) matching `query` case-insensitively (exact or `_`/` ` normalized).
    pub fn find_taxon_node(&self, query: &str) -> Option<usize> {
        let norm = |s: &str| s.trim().replace('_', " ").to_ascii_lowercase();
        let q = norm(query);
        for u in self.preorder_ids(self.root) {
            if self.nodes[u].is_leaf() && norm(&self.nodes[u].name) == q {
                return Some(u);
            }
        }
        for u in self.preorder_ids(self.root) {
            if norm(&self.nodes[u].name).contains(&q) {
                return Some(u);
            }
        }
        None
    }

    /// Computes summary topology & branch-length statistics across all reachable nodes.
    pub fn stats(&self) -> PhyloTreeStats {
        let mut st = PhyloTreeStats {
            min_root_to_tip: f64::INFINITY,
            ..Default::default()
        };
        let mut depth_map = vec![0.0_f64; self.nodes.len()];
        let mut sup_sum = 0.0_f64;
        let mut sup_cnt = 0usize;

        for u in self.preorder_ids(self.root) {
            let node = &self.nodes[u];
            if u != self.root {
                let bl = node.branch_length.unwrap_or(0.0).max(0.0);
                st.total_branch_length += bl;
                let p_depth = node.parent.map(|p| depth_map[p]).unwrap_or(0.0);
                depth_map[u] = p_depth + bl;
            }
            if node.is_leaf() {
                st.leaf_count += 1;
                st.max_root_to_tip = st.max_root_to_tip.max(depth_map[u]);
                st.min_root_to_tip = st.min_root_to_tip.min(depth_map[u]);
            } else {
                st.internal_count += 1;
                if node.children.len() > 2 {
                    st.polytomy_count += 1;
                }
                if node.collapsed {
                    st.collapsed_clades += 1;
                }
                if let Some(s) = node.support_pct {
                    sup_sum += s;
                    sup_cnt += 1;
                }
            }
        }

        if st.leaf_count == 0 {
            st.min_root_to_tip = 0.0;
        }
        if sup_cnt > 0 {
            st.mean_support_pct = Some(sup_sum / (sup_cnt as f64));
        }
        st
    }

    /// Ladderizes every internal node by sorting children according to clade leaf count
    /// (`ascending = true` places smaller clades at the top, larger clades at the bottom).
    pub fn ladderize(&mut self, ascending: bool) {
        let mut tip_counts = vec![0usize; self.nodes.len()];
        let mut min_label = vec![String::new(); self.nodes.len()];

        for u in self.postorder_ids(self.root) {
            if self.nodes[u].is_leaf() {
                tip_counts[u] = 1;
                min_label[u] = self.nodes[u].name.clone();
            } else {
                let mut sum = 0usize;
                let mut smallest = String::from("~");
                for &ch in &self.nodes[u].children {
                    sum += tip_counts[ch];
                    if min_label[ch] < smallest {
                        smallest = min_label[ch].clone();
                    }
                }
                tip_counts[u] = sum;
                min_label[u] = smallest;

                self.nodes[u].children.sort_by(|&a, &b| {
                    let cmp = tip_counts[a]
                        .cmp(&tip_counts[b])
                        .then_with(|| min_label[a].cmp(&min_label[b]));
                    if ascending {
                        cmp
                    } else {
                        cmp.reverse()
                    }
                });
            }
        }
    }

    /// Toggles collapsed state on an internal node. Returns the new `collapsed` boolean.
    pub fn toggle_collapse(&mut self, node_id: usize) -> bool {
        if let Some(n) = self.nodes.get_mut(node_id) {
            if !n.is_leaf() && node_id != self.root {
                n.collapsed = !n.collapsed;
                return n.collapsed;
            }
        }
        false
    }

    /// Collapses all internal clades whose bootstrap support falls strictly below `min_support_pct`.
    pub fn collapse_weak_clades(&mut self, min_support_pct: f64) -> usize {
        let mut count = 0usize;
        for u in self.preorder_ids(self.root) {
            if u != self.root && !self.nodes[u].is_leaf() {
                if let Some(sup) = self.nodes[u].support_pct {
                    if sup < min_support_pct {
                        self.nodes[u].collapsed = true;
                        count += 1;
                    }
                }
            }
        }
        count
    }

    /// Computes the Most Recent Common Ancestor (MRCA) node ID of a slice of node IDs.
    pub fn mrca_of_nodes(&self, targets: &[usize]) -> Option<usize> {
        let first = *targets.first()?;
        let mut common_ancestors = self.ancestors_to_root(first);
        for &t in targets.iter().skip(1) {
            let anc_t: BTreeSet<usize> = self.ancestors_to_root(t).into_iter().collect();
            common_ancestors.retain(|a| anc_t.contains(a));
        }
        common_ancestors.first().copied()
    }

    fn ancestors_to_root(&self, mut u: usize) -> Vec<usize> {
        let mut chain = vec![u];
        while let Some(p) = self.nodes[u].parent {
            chain.push(p);
            u = p;
        }
        chain
    }

    /// Evaluates whether `query_taxa` form a monophyletic, paraphyletic, or polyphyletic group.
    pub fn check_monophyly(&self, query_taxa: &[String]) -> Result<MonophylyReport, String> {
        let mut matched_ids = Vec::new();
        let mut matched_names = Vec::new();
        let mut missing = Vec::new();

        for q in query_taxa {
            if let Some(id) = self.find_taxon_node(q) {
                if !matched_ids.contains(&id) {
                    matched_ids.push(id);
                    matched_names.push(self.nodes[id].name.clone());
                }
            } else {
                missing.push(q.clone());
            }
        }

        if matched_ids.is_empty() {
            return Err(String::from("None of the requested taxa were found in the tree."));
        }

        let mrca = self
            .mrca_of_nodes(&matched_ids)
            .ok_or_else(|| String::from("Failed computing MRCA across disconnected tree."))?;
        let clade_leaves = self.leaf_ids_under(mrca);
        let matched_set: BTreeSet<usize> = matched_ids.iter().copied().collect();

        let intruders: Vec<String> = clade_leaves
            .iter()
            .filter(|&&lid| !matched_set.contains(&lid))
            .map(|&lid| self.nodes[lid].name.clone())
            .collect();

        let status = if intruders.is_empty() {
            MonophylyStatus::Monophyletic
        } else {
            // Distinguish paraphyly (single contiguous grade excluding 1..2 subclades) vs polyphyly:
            // Count minimal monophyletic target subclades required to cover `matched_set`.
            let subclades = self.count_minimal_monophyletic_clusters(mrca, &matched_set);
            if subclades <= 2 && intruders.len() <= matched_set.len() {
                MonophylyStatus::Paraphyletic
            } else {
                MonophylyStatus::Polyphyletic
            }
        };

        Ok(MonophylyReport {
            status,
            mrca_node_id: mrca,
            matched_taxa: matched_names,
            missing_query_taxa: missing,
            intruder_taxa: intruders,
            mrca_support_pct: self.nodes[mrca].support_pct,
            mrca_stem_length: self.nodes[mrca].branch_length,
        })
    }

    fn count_minimal_monophyletic_clusters(
        &self,
        node: usize,
        targets: &BTreeSet<usize>,
    ) -> usize {
        let leaves = self.leaf_ids_under(node);
        let target_hits = leaves.iter().filter(|l| targets.contains(l)).count();
        if target_hits == 0 {
            return 0;
        }
        if target_hits == leaves.len() {
            return 1;
        }
        self.nodes[node]
            .children
            .iter()
            .map(|&ch| self.count_minimal_monophyletic_clusters(ch, targets))
            .sum()
    }

    /// Reroots the tree along the edge leading into `target_node`, placing the new root at
    /// `fraction_from_parent` (`0.5` splits the target branch in half).
    pub fn reroot_at_node(&mut self, target_node: usize, fraction_from_parent: f64) -> Result<(), String> {
        if target_node == self.root {
            return Ok(());
        }
        let Some(orig_parent) = self.nodes[target_node].parent else {
            return Err(String::from("Target node has no parent branch to split."));
        };

        let frac = fraction_from_parent.clamp(0.0, 1.0);
        let orig_bl = self.nodes[target_node].branch_length.unwrap_or(0.0);
        let orig_sup = self.nodes[target_node].support_pct;

        // Collect path from `orig_parent` up to `self.root`: [orig_parent, ..., old_root]
        let mut path = Vec::new();
        let mut cur = orig_parent;
        loop {
            path.push(cur);
            if cur == self.root {
                break;
            }
            cur = self.nodes[cur]
                .parent
                .ok_or_else(|| String::from("Broken parent chain in tree."))?;
        }

        // Detach `target_node` from `orig_parent`
        self.nodes[orig_parent].children.retain(|&c| c != target_node);

        // Reverse directed edges along `path` from `orig_parent` up to `old_root`
        for win in path.windows(2).rev() {
            let child = win[0];
            let par = win[1];
            let edge_bl = self.nodes[child].branch_length;
            let edge_sup = self.nodes[child].support_pct;

            self.nodes[par].children.retain(|&c| c != child);
            self.nodes[child].children.push(par);
            self.nodes[par].parent = Some(child);
            self.nodes[par].branch_length = edge_bl;
            self.nodes[par].support_pct = edge_sup;
        }

        // Suppress degree-2 old root if it now has only 1 remaining child (`len() == 1` after reversal)
        let old_root = self.root;
        if self.nodes[old_root].children.len() == 1 && self.nodes[old_root].parent.is_some() {
            let sole_child = self.nodes[old_root].children[0];
            let new_par_of_old_root = self.nodes[old_root].parent.unwrap();
            let merged_bl = self.nodes[old_root].branch_length.unwrap_or(0.0)
                + self.nodes[sole_child].branch_length.unwrap_or(0.0);
            let has_bl = self.nodes[old_root].branch_length.is_some()
                || self.nodes[sole_child].branch_length.is_some();

            for slot in &mut self.nodes[new_par_of_old_root].children {
                if *slot == old_root {
                    *slot = sole_child;
                }
            }
            self.nodes[sole_child].parent = Some(new_par_of_old_root);
            if has_bl {
                self.nodes[sole_child].branch_length = Some(merged_bl);
            }
            self.nodes[old_root].children.clear();
            self.nodes[old_root].parent = None;
        }

        // Create new bifurcating root node between `target_node` and `orig_parent`
        let new_root_id = self.nodes.len();
        self.nodes.push(PhyloNode {
            id: new_root_id,
            parent: None,
            children: vec![target_node, orig_parent],
            name: String::from("Root"),
            branch_length: None,
            support_pct: None,
            annotations: BTreeMap::new(),
            collapsed: false,
        });

        self.nodes[target_node].parent = Some(new_root_id);
        self.nodes[orig_parent].parent = Some(new_root_id);

        if self.nodes[target_node].branch_length.is_some() {
            self.nodes[target_node].branch_length = Some(orig_bl * (1.0 - frac));
            self.nodes[orig_parent].branch_length = Some(orig_bl * frac);
        } else {
            self.nodes[orig_parent].branch_length = None;
        }
        self.nodes[orig_parent].support_pct = orig_sup;

        self.root = new_root_id;
        Ok(())
    }

    /// Reroots the tree using `outgroup_taxon` (matched case-insensitively).
    pub fn reroot_by_outgroup(&mut self, outgroup_taxon: &str) -> Result<(), String> {
        let tid = self
            .find_taxon_node(outgroup_taxon)
            .ok_or_else(|| format!("Outgroup taxon '{outgroup_taxon}' not found in tree."))?;
        self.reroot_at_node(tid, 0.5)
    }

    /// Midpoint-reroots the tree at the exact half-distance of the maximum leaf-to-leaf patristic path.
    pub fn midpoint_reroot(&mut self) -> Result<f64, String> {
        let leaves = self.leaf_ids_under(self.root);
        if leaves.len() < 3 {
            return Err(String::from("Midpoint rooting requires >= 3 tips."));
        }

        let mut best_pair = (leaves[0], leaves[1], 0.0_f64);
        for (i, &u) in leaves.iter().enumerate() {
            for &v in leaves.iter().skip(i + 1) {
                let d = self.patristic_distance(u, v);
                if d > best_pair.2 {
                    best_pair = (u, v, d);
                }
            }
        }

        let (u, v, max_dist) = best_pair;
        if max_dist <= 0.0 {
            return Err(String::from(
                "Tree lacks positive branch lengths required for midpoint rooting.",
            ));
        }

        let half = max_dist * 0.5;
        let path_edges = self.directed_path_between(u, v)?;
        let mut accum = 0.0_f64;

        for (from_node, _to_node, child_with_bl, bl) in path_edges {
            if accum + bl >= half && bl > 0.0 {
                let offset_from_from = half - accum;
                // `child_with_bl` is whichever endpoint (`from_node` or `to_node`) is currently the child of the other
                let frac_from_parent = if child_with_bl == from_node {
                    // Moving from child (`from_node`) up toward parent (`to_node`): distance from parent is `bl - offset`
                    ((bl - offset_from_from) / bl).clamp(0.0, 1.0)
                } else {
                    // Moving from parent (`from_node`) down toward child (`to_node`): distance from parent is `offset`
                    (offset_from_from / bl).clamp(0.0, 1.0)
                };
                self.reroot_at_node(child_with_bl, frac_from_parent)?;
                return Ok(max_dist);
            }
            accum += bl;
        }

        Ok(max_dist)
    }

    /// Computes patristic branch-length distance between any two nodes `u` and `v`.
    pub fn patristic_distance(&self, u: usize, v: usize) -> f64 {
        let Some(mrca) = self.mrca_of_nodes(&[u, v]) else {
            return 0.0;
        };
        self.dist_to_ancestor(u, mrca) + self.dist_to_ancestor(v, mrca)
    }

    fn dist_to_ancestor(&self, mut cur: usize, anc: usize) -> f64 {
        let mut sum = 0.0_f64;
        while cur != anc {
            sum += self.nodes[cur].branch_length.unwrap_or(0.0).max(0.0);
            match self.nodes[cur].parent {
                Some(p) => cur = p,
                None => break,
            }
        }
        sum
    }

    fn directed_path_between(
        &self,
        u: usize,
        v: usize,
    ) -> Result<Vec<(usize, usize, usize, f64)>, String> {
        let mrca = self
            .mrca_of_nodes(&[u, v])
            .ok_or_else(|| String::from("Nodes do not share an ancestor."))?;
        let mut edges = Vec::new();

        let mut cur = u;
        while cur != mrca {
            let p = self.nodes[cur].parent.unwrap();
            let bl = self.nodes[cur].branch_length.unwrap_or(0.0).max(0.0);
            edges.push((cur, p, cur, bl));
            cur = p;
        }

        let mut down = Vec::new();
        let mut cur_v = v;
        while cur_v != mrca {
            let p = self.nodes[cur_v].parent.unwrap();
            let bl = self.nodes[cur_v].branch_length.unwrap_or(0.0).max(0.0);
            down.push((p, cur_v, cur_v, bl));
            cur_v = p;
        }
        down.reverse();
        edges.extend(down);
        Ok(edges)
    }

    /// Serializes the reachable tree from `self.root` into standard Newick format.
    pub fn to_newick(&self) -> String {
        let mut buf = String::with_capacity(self.nodes.len() * 24);
        self.write_newick_rec(self.root, &mut buf);
        buf.push(';');
        buf
    }

    fn write_newick_rec(&self, u: usize, buf: &mut String) {
        let node = &self.nodes[u];
        if !node.children.is_empty() {
            buf.push('(');
            for (idx, &ch) in node.children.iter().enumerate() {
                if idx > 0 {
                    buf.push(',');
                }
                self.write_newick_rec(ch, buf);
            }
            buf.push(')');
        }

        if node.is_leaf() {
            push_newick_label(buf, &node.name);
        } else if let Some(sup) = node.support_pct {
            if (sup.fract()).abs() < 1e-4 {
                let _ = write!(buf, "{sup:.0}");
            } else {
                let _ = write!(buf, "{sup:.1}");
            }
        } else if !node.name.is_empty() && u != self.root {
            push_newick_label(buf, &node.name);
        }

        if let Some(bl) = node.branch_length {
            let _ = write!(buf, ":{bl:.6}");
        }
    }
}

fn push_newick_label(buf: &mut String, label: &str) {
    if label.is_empty() {
        return;
    }
    let needs_quotes = label
        .chars()
        .any(|c| matches!(c, '(' | ')' | ',' | ':' | ';' | '[' | ']' | ' '));
    if needs_quotes {
        buf.push('\'');
        for ch in label.chars() {
            if ch == '\'' {
                buf.push_str("''");
            } else {
                buf.push(ch);
            }
        }
        buf.push('\'');
    } else {
        buf.push_str(label);
    }
}

/// Auto-detects Maddison et al. (1997) NEXUS (`#NEXUS`) vs multi-line Newick (`.nwk`) files
/// and parses all trees contained in `text`.
pub fn parse_phylo_trees(default_prefix: &str, text: &str) -> Result<Vec<PhyloTree>, String> {
    let trimmed = text.trim_start_matches('\u{feff}').trim();
    if trimmed.is_empty() {
        return Err(String::from("Tree file is empty."));
    }
    if trimmed
        .get(0..6)
        .map(|s| s.eq_ignore_ascii_case("#NEXUS"))
        .unwrap_or(false)
    {
        return parse_nexus_trees(trimmed);
    }

    // Parse 1 or more semicolon-terminated Newick statements
    let mut trees = Vec::new();
    for (idx, chunk) in split_newick_statements(trimmed).into_iter().enumerate() {
        let name = if idx == 0 {
            default_prefix.to_string()
        } else {
            format!("{default_prefix}_{}", idx + 1)
        };
        trees.push(parse_newick(&name, &format!("{chunk};"))?);
    }
    if trees.is_empty() {
        return Err(String::from("No valid Newick tree statements found."));
    }
    Ok(trees)
}

fn split_newick_statements(s: &str) -> Vec<String> {
    let mut out = Vec::new();
    let mut cur = String::new();
    let mut in_quote = false;
    let mut bracket_depth = 0usize;

    for ch in s.chars() {
        if in_quote {
            if ch == '\'' {
                in_quote = false;
            }
            cur.push(ch);
        } else if bracket_depth > 0 {
            if ch == '[' {
                bracket_depth += 1;
            } else if ch == ']' {
                bracket_depth = bracket_depth.saturating_sub(1);
            }
            cur.push(ch);
        } else {
            match ch {
                '\'' => {
                    in_quote = true;
                    cur.push(ch);
                }
                '[' => {
                    bracket_depth = 1;
                    cur.push(ch);
                }
                ';' => {
                    let t = cur.trim();
                    if !t.is_empty() {
                        out.push(t.to_string());
                    }
                    cur.clear();
                }
                _ => cur.push(ch),
            }
        }
    }
    let tail = cur.trim();
    if !tail.is_empty() && tail.starts_with('(') {
        out.push(tail.to_string());
    }
    out
}

/// Parses a Maddison, Swofford & Maddison (1997) `#NEXUS` file (`BEGIN TREES;`, optional `TRANSLATE` map).
pub fn parse_nexus_trees(text: &str) -> Result<Vec<PhyloTree>, String> {
    let text = text.trim_start_matches('\u{feff}');
    let upper = text.to_ascii_uppercase();
    let begin_idx = upper
        .find("BEGIN TREES;")
        .or_else(|| upper.find("BEGIN TREES"))
        .ok_or_else(|| String::from("NEXUS file does not contain a `BEGIN TREES;` block."))?;
    let block_rest = &text[begin_idx..];
    let block_upper = block_rest.to_ascii_uppercase();
    let end_idx = block_upper
        .find("END;")
        .or_else(|| block_upper.find("ENDBLOCK;"))
        .unwrap_or(block_rest.len());
    let trees_block = &block_rest[..end_idx];

    let mut translate_map: HashMap<String, String> = HashMap::new();
    let mut trees = Vec::new();

    for stmt in split_newick_statements(trees_block) {
        let s_trim = stmt.trim();
        let s_up = s_trim.to_ascii_uppercase();
        if s_up.starts_with("BEGIN ") {
            continue;
        }
        if s_up.starts_with("TRANSLATE") {
            let body = s_trim["TRANSLATE".len()..].trim();
            for pair in body.split(',') {
                let p = pair.trim();
                if p.is_empty() {
                    continue;
                }
                let mut toks = p.splitn(2, char::is_whitespace);
                if let (Some(k), Some(v)) = (toks.next(), toks.next()) {
                    let clean_v = v.trim().trim_matches('\'').trim_matches('"').to_string();
                    translate_map.insert(k.trim().to_string(), clean_v);
                }
            }
            continue;
        }
        if s_up.starts_with("TREE ") || s_up.starts_with("UTREE ") {
            if let Some((lhs, rhs)) = s_trim.split_once('=') {
                let tname = lhs
                    .split_whitespace()
                    .nth(1)
                    .unwrap_or("NexusTree")
                    .trim_start_matches('*')
                    .trim_matches('\'')
                    .to_string();
                let mut nwk_rhs = rhs.trim();
                // Strip leading NEXUS rooting comment `[&R]` or `[&U]` before `(`
                while nwk_rhs.starts_with('[') {
                    if let Some(close) = nwk_rhs.find(']') {
                        nwk_rhs = nwk_rhs[(close + 1)..].trim();
                    } else {
                        break;
                    }
                }
                let mut tree = parse_newick(&tname, &format!("{nwk_rhs};"))?;
                if !translate_map.is_empty() {
                    for node in &mut tree.nodes {
                        if node.is_leaf() {
                            if let Some(mapped) = translate_map.get(&node.name) {
                                node.name = mapped.clone();
                            }
                        }
                    }
                }
                trees.push(tree);
            }
        }
    }

    if trees.is_empty() {
        return Err(String::from("No `TREE name = (...)` definitions found in NEXUS TREES block."));
    }
    Ok(trees)
}

/// Recursive-descent parser for a single Newick / NHX tree string.
pub fn parse_newick(tree_name: &str, newick: &str) -> Result<PhyloTree, String> {
    let bytes: Vec<char> = newick.trim_start_matches('\u{feff}').trim().chars().collect();
    let mut pos = 0usize;
    let mut nodes: Vec<PhyloNode> = Vec::new();

    let root = parse_subtree(&bytes, &mut pos, &mut nodes, None)?;
    skip_ws(&bytes, &mut pos);
    if pos < bytes.len() && bytes[pos] == ';' {
        let _ = pos;
    }

    Ok(PhyloTree {
        tree_name: tree_name.to_string(),
        nodes,
        root,
    })
}

fn parse_subtree(
    chars: &[char],
    pos: &mut usize,
    nodes: &mut Vec<PhyloNode>,
    parent: Option<usize>,
) -> Result<usize, String> {
    skip_ws(chars, pos);
    let my_id = nodes.len();
    nodes.push(PhyloNode {
        id: my_id,
        parent,
        children: Vec::new(),
        name: String::new(),
        branch_length: None,
        support_pct: None,
        annotations: BTreeMap::new(),
        collapsed: false,
    });

    if *pos < chars.len() && chars[*pos] == '(' {
        *pos += 1;
        loop {
            let child_id = parse_subtree(chars, pos, nodes, Some(my_id))?;
            nodes[my_id].children.push(child_id);
            skip_ws(chars, pos);
            if *pos < chars.len() && chars[*pos] == ',' {
                *pos += 1;
                continue;
            } else if *pos < chars.len() && chars[*pos] == ')' {
                *pos += 1;
                break;
            } else {
                return Err(format!(
                    "Malformed Newick at char {}: expected ',' or ')'.",
                    *pos
                ));
            }
        }
    }

    skip_ws(chars, pos);
    let label = parse_label(chars, pos);
    let mut annotations = BTreeMap::new();
    parse_optional_comments(chars, pos, &mut annotations);

    let mut branch_length: Option<f64> = None;
    skip_ws(chars, pos);
    if *pos < chars.len() && chars[*pos] == ':' {
        *pos += 1;
        skip_ws(chars, pos);
        parse_optional_comments(chars, pos, &mut annotations);
        let num_str = parse_unquoted_token(chars, pos);
        if !num_str.is_empty() {
            let val: f64 = num_str
                .parse()
                .map_err(|_| format!("Invalid branch length '{num_str}' in Newick."))?;
            branch_length = Some(val);
        }
        parse_optional_comments(chars, pos, &mut annotations);
    }

    let is_leaf = nodes[my_id].children.is_empty();
    let mut support_pct: Option<f64> = None;
    let mut final_name = label;

    if !is_leaf && !final_name.is_empty() {
        // Handle IQ-TREE `SH-aLRT/UFBoot` like `98.4/100` or pure numeric bootstrap `99` or posterior `0.98`
        let numeric_candidate = final_name
            .rsplit('/')
            .next()
            .unwrap_or(&final_name)
            .trim();
        if let Ok(v) = numeric_candidate.parse::<f64>() {
            support_pct = Some(if v <= 1.0 && numeric_candidate.contains('.') {
                v * 100.0
            } else {
                v
            });
            if final_name.parse::<f64>().is_ok() || final_name.contains('/') {
                final_name.clear();
            }
        }
    }

    if support_pct.is_none() {
        if let Some(b_str) = annotations
            .get("B")
            .or_else(|| annotations.get("posterior"))
            .or_else(|| annotations.get("bootstrap"))
        {
            if let Ok(v) = b_str.parse::<f64>() {
                support_pct = Some(if v <= 1.0 && b_str.contains('.') {
                    v * 100.0
                } else {
                    v
                });
            }
        }
    }

    nodes[my_id].name = final_name;
    nodes[my_id].branch_length = branch_length;
    nodes[my_id].support_pct = support_pct;
    nodes[my_id].annotations = annotations;
    Ok(my_id)
}

fn skip_ws(chars: &[char], pos: &mut usize) {
    while *pos < chars.len() && chars[*pos].is_whitespace() {
        *pos += 1;
    }
}

fn parse_label(chars: &[char], pos: &mut usize) -> String {
    skip_ws(chars, pos);
    if *pos >= chars.len() {
        return String::new();
    }
    if chars[*pos] == '\'' || chars[*pos] == '"' {
        let quote = chars[*pos];
        *pos += 1;
        let mut out = String::new();
        while *pos < chars.len() {
            let c = chars[*pos];
            *pos += 1;
            if c == quote {
                if *pos < chars.len() && chars[*pos] == quote {
                    out.push(quote);
                    *pos += 1;
                } else {
                    break;
                }
            } else {
                out.push(c);
            }
        }
        return out;
    }
    parse_unquoted_token(chars, pos)
}

fn parse_unquoted_token(chars: &[char], pos: &mut usize) -> String {
    let mut out = String::new();
    while *pos < chars.len() {
        let c = chars[*pos];
        if c.is_whitespace() || matches!(c, '(' | ')' | ',' | ':' | ';' | '[' | ']') {
            break;
        }
        out.push(c);
        *pos += 1;
    }
    out
}

fn parse_optional_comments(
    chars: &[char],
    pos: &mut usize,
    annotations: &mut BTreeMap<String, String>,
) {
    loop {
        skip_ws(chars, pos);
        if *pos >= chars.len() || chars[*pos] != '[' {
            break;
        }
        *pos += 1;
        let mut raw = String::new();
        let mut depth = 1usize;
        while *pos < chars.len() && depth > 0 {
            let c = chars[*pos];
            *pos += 1;
            if c == '[' {
                depth += 1;
            } else if c == ']' {
                depth -= 1;
                if depth == 0 {
                    break;
                }
            }
            raw.push(c);
        }

        // Parse NHX `[&&NHX:S=Human:B=100]` or Beast `[&posterior=0.99,rate=0.01]`
        let trimmed = raw
            .trim()
            .trim_start_matches("&&NHX:")
            .trim_start_matches('&');
        let sep = if raw.contains("&&NHX:") { ':' } else { ',' };
        for kv in trimmed.split(sep) {
            if let Some((k, v)) = kv.split_once('=') {
                annotations.insert(k.trim().to_string(), v.trim().to_string());
            }
        }
    }
}

/// Computes the Robinson & Foulds (1981) unrooted symmetric bipartition distance and
/// Kuhner–Felsenstein branch-length L1 score across the shared leaf set of `tree_a` and `tree_b`.
pub fn compare_robinson_foulds(
    tree_a: &PhyloTree,
    tree_b: &PhyloTree,
) -> Result<RobinsonFouldsReport, String> {
    let set_a: BTreeSet<String> = tree_a.leaf_names_under(tree_a.root).into_iter().collect();
    let set_b: BTreeSet<String> = tree_b.leaf_names_under(tree_b.root).into_iter().collect();
    let shared: BTreeSet<String> = set_a.intersection(&set_b).cloned().collect();

    if shared.len() < 4 {
        return Err(format!(
            "Robinson-Foulds comparison requires >= 4 shared taxa (found {}).",
            shared.len()
        ));
    }

    let anchor = shared.iter().next().unwrap().clone();
    let splits_a = extract_canonical_splits(tree_a, &shared, &anchor);
    let splits_b = extract_canonical_splits(tree_b, &shared, &anchor);

    let mut shared_count = 0usize;
    let mut unique_a = Vec::new();
    let mut unique_b = Vec::new();
    let mut l1_sum = 0.0_f64;

    for (sp, &bl_a) in &splits_a {
        if let Some(&bl_b) = splits_b.get(sp) {
            shared_count += 1;
            l1_sum += (bl_a - bl_b).abs();
        } else {
            unique_a.push(sp.iter().cloned().collect());
            l1_sum += bl_a;
        }
    }
    for (sp, &bl_b) in &splits_b {
        if !splits_a.contains_key(sp) {
            unique_b.push(sp.iter().cloned().collect());
            l1_sum += bl_b;
        }
    }

    let rf = unique_a.len() + unique_b.len();
    let max_rf = splits_a.len() + splits_b.len();
    let norm_rf = if max_rf > 0 {
        (rf as f64) / (max_rf as f64)
    } else {
        0.0
    };

    Ok(RobinsonFouldsReport {
        shared_taxa_count: shared.len(),
        splits_tree_a: splits_a.len(),
        splits_tree_b: splits_b.len(),
        shared_splits: shared_count,
        rf_distance: rf,
        max_rf,
        normalized_rf: norm_rf,
        branch_score_l1: l1_sum,
        unique_to_a: unique_a,
        unique_to_b: unique_b,
    })
}

fn extract_canonical_splits(
    tree: &PhyloTree,
    shared: &BTreeSet<String>,
    anchor_taxon: &str,
) -> BTreeMap<BTreeSet<String>, f64> {
    let mut splits: BTreeMap<BTreeSet<String>, f64> = BTreeMap::new();
    let n_shared = shared.len();

    for u in tree.postorder_ids(tree.root) {
        if u == tree.root {
            continue;
        }
        let mut side: BTreeSet<String> = tree
            .leaf_names_under(u)
            .into_iter()
            .filter(|t| shared.contains(t))
            .collect();

        if side.len() <= 1 || side.len() >= n_shared - 1 {
            continue; // Trivial terminal split
        }

        // Canonicalize unrooted split: always choose the partition half NOT containing `anchor_taxon`
        if side.contains(anchor_taxon) {
            side = shared.difference(&side).cloned().collect();
        }
        let bl = tree.nodes[u].branch_length.unwrap_or(0.0).max(0.0);
        // If root is bifurcating, its two child edges represent the SAME unrooted split — sum their lengths!
        *splits.entry(side).or_insert(0.0) += bl;
    }

    splits
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parses_nexus_reroot_invariance_and_monophyly() {
        let mut demo = PhyloTree::demo_vertebrate_amniotes();
        let st = demo.stats();
        assert_eq!(st.leaf_count, 14);

        // Amniotes should be strictly monophyletic
        let amniotes = vec![
            "Homo_sapiens".to_string(),
            "Pan_troglodytes".to_string(),
            "Mus_musculus".to_string(),
            "Rattus_norvegicus".to_string(),
            "Ornithorhynchus_anatinus".to_string(),
            "Gallus_gallus".to_string(),
            "Taeniopygia_guttata".to_string(),
            "Alligator_mississippiensis".to_string(),
        ];
        let mono = demo.check_monophyly(&amniotes).unwrap();
        assert_eq!(mono.status, MonophylyStatus::Monophyletic);

        // Reptiles without birds (`Alligator` + `Xenopus`?) or `Gallus` excluded -> Paraphyletic/Polyphyletic
        let mut copy = demo.clone();
        copy.midpoint_reroot().unwrap();

        // Unrooted Robinson-Foulds distance must be 0 after midpoint rerooting!
        let rf = compare_robinson_foulds(&demo, &copy).unwrap();
        assert_eq!(rf.rf_distance, 0);

        // NEXUS with TRANSLATE block
        let nex = "#NEXUS\nBEGIN TREES;\n  TRANSLATE\n    1 Homo_sapiens,\n    2 Pan_troglodytes,\n    3 Mus_musculus,\n    4 Gallus_gallus;\n  TREE t1 = [&R] ((1:0.03,2:0.03)100:0.07,(3:0.1,4:0.2):0.05);\nEND;\n";
        let parsed = parse_phylo_trees("nex", nex).unwrap();
        assert_eq!(parsed.len(), 1);
        assert_eq!(parsed[0].stats().leaf_count, 4);
        assert!(parsed[0].find_taxon_node("Homo_sapiens").is_some());
        let _ = &mut demo;
    }
}
