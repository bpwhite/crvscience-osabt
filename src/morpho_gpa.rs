// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #8 (`osabt morpho`) — 2D Geometric Morphometrics engine:
//! 1. Rohlf `tpsDig` `.tps` landmark file reader & writer (`LM=`, `IMAGE=`, `ID=`, `SCALE=`).
//! 2. Rohlf & Slice (1990) Generalized Procrustes Analysis (`GPA`) removing non-shape translation,
//!    isometric Centroid Size (`CS`), and rigid planar rotation without reflection.
//! 3. Cyclic Jacobi symmetric eigendecomposition Principal Component Analysis (`PCA`) on tangent
//!    shape coordinates plus publication `.svg` Procrustes superimposition & PCA morphospace plots.

use std::fmt::Write as _;

/// Single specimen's 2D homologous landmark configuration.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct LandmarkSpecimen2D {
    pub specimen_id: String,
    pub image_path: String,
    pub scale_factor: Option<f64>,
    pub landmarks: Vec<[f64; 2]>,
}

/// Complete output of Rohlf & Slice (1990) 2D Generalized Procrustes Analysis + tangent-space PCA.
#[derive(Debug, Clone)]
pub struct GpaAnalysis2D {
    pub specimen_ids: Vec<String>,
    pub num_landmarks: usize,
    /// Calibrated Centroid Size (`CS`) for each specimen.
    pub centroid_sizes: Vec<f64>,
    /// Consensus (mean) Procrustes shape (`num_landmarks` `[x, y]` pairs).
    pub consensus_shape: Vec<[f64; 2]>,
    /// Superimposed Procrustes coordinates for every specimen (`[specimen][landmark]`).
    pub aligned_shapes: Vec<Vec<[f64; 2]>>,
    /// Procrustes root-mean-square distance of each specimen from `consensus_shape`.
    pub procrustes_distances_to_mean: Vec<f64>,
    /// `(PC1_score, PC2_score)` for each specimen in Procrustes tangent morphospace.
    pub pc12_scores: Vec<[f64; 2]>,
    pub pc1_variance_pct: f64,
    pub pc2_variance_pct: f64,
}

/// Parses a Rohlf `tpsDig` `.tps` file containing one or more `LM=K` landmark blocks.
pub fn parse_tps_landmarks(tps_text: &str) -> Result<Vec<LandmarkSpecimen2D>, String> {
    let mut specimens = Vec::new();
    let lines: Vec<&str> = tps_text.lines().collect();
    let mut i = 0usize;

    while i < lines.len() {
        let line = lines[i].trim();
        if line.is_empty() || line.starts_with('#') {
            i += 1;
            continue;
        }
        let up = line.to_ascii_uppercase();
        if let Some(rest) = up.strip_prefix("LM=").or_else(|| up.strip_prefix("LM2=")) {
            let k: usize = rest
                .trim()
                .parse()
                .map_err(|_| format!("Invalid landmark count in '{line}'"))?;
            i += 1;
            let mut pts = Vec::with_capacity(k);
            while pts.len() < k && i < lines.len() {
                let pline = lines[i].trim();
                i += 1;
                if pline.is_empty() {
                    continue;
                }
                let mut toks = pline.split_whitespace();
                let x: f64 = toks
                    .next()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| format!("Malformed X landmark coordinate: '{pline}'"))?;
                let y: f64 = toks
                    .next()
                    .and_then(|s| s.parse().ok())
                    .ok_or_else(|| format!("Malformed Y landmark coordinate: '{pline}'"))?;
                pts.push([x, y]);
            }
            if pts.len() != k {
                return Err(format!("Expected {k} landmarks, found {}.", pts.len()));
            }

            let mut image_path = String::new();
            let mut specimen_id = format!("Specimen_{}", specimens.len() + 1);
            let mut scale_factor = None;

            while i < lines.len() {
                let mline = lines[i].trim();
                let mup = mline.to_ascii_uppercase();
                if mup.starts_with("LM=") || mup.starts_with("LM2=") {
                    break;
                }
                i += 1;
                if let Some( val) = mline.split_once('=').and_then(|(k_str, v)| {
                    if k_str.trim().eq_ignore_ascii_case("IMAGE") {
                        Some(v.trim())
                    } else {
                        None
                    }
                }) {
                    image_path = val.to_string();
                } else if let Some(val) = mline.split_once('=').and_then(|(k_str, v)| {
                    if k_str.trim().eq_ignore_ascii_case("ID") {
                        Some(v.trim())
                    } else {
                        None
                    }
                }) {
                    specimen_id = val.to_string();
                } else if let Some(val) = mline.split_once('=').and_then(|(k_str, v)| {
                    if k_str.trim().eq_ignore_ascii_case("SCALE") {
                        v.trim().parse::<f64>().ok()
                    } else {
                        None
                    }
                }) {
                    scale_factor = Some(val);
                }
            }

            specimens.push(LandmarkSpecimen2D {
                specimen_id,
                image_path,
                scale_factor,
                landmarks: pts,
            });
        } else {
            i += 1;
        }
    }

    if specimens.is_empty() {
        return Err(String::from("No `LM=` blocks found in .tps input."));
    }
    Ok(specimens)
}

/// Serializes homologous 2D landmark specimens into standard Rohlf `tpsDig` `.tps` format.
pub fn write_tps_landmarks(specimens: &[LandmarkSpecimen2D]) -> String {
    let mut out = String::new();
    for sp in specimens {
        let _ = writeln!(out, "LM={}", sp.landmarks.len());
        for &[x, y] in &sp.landmarks {
            let _ = writeln!(out, "{x:.5} {y:.5}");
        }
        if !sp.image_path.is_empty() {
            let _ = writeln!(out, "IMAGE={}", sp.image_path);
        }
        let _ = writeln!(out, "ID={}", sp.specimen_id);
        if let Some(sc) = sp.scale_factor {
            let _ = writeln!(out, "SCALE={sc:.6}");
        }
        out.push('\n');
    }
    out
}

/// Built-in 6-specimen, 8-landmark salmonid (*Oncorhynchus*) morphometrics `.tps` dataset
/// spanning 3 stream-resident slender morphs and 3 deep-bodied lacustrine morphs across
/// arbitrary camera positions, scales, and orientations.
pub fn demo_salmonid_landmarks_tps() -> &'static str {
    "\
LM=8
10.0 50.0
28.0 61.0
62.0 64.0
98.0 56.0
118.0 50.0
98.0 44.0
62.0 36.0
28.0 39.0
IMAGE=stream_01.bmp
ID=Stream_Morph_01
SCALE=0.250000

LM=8
21.0 80.5
40.2 92.0
76.5 95.1
114.0 86.2
135.2 80.1
113.8 73.9
76.1 65.2
39.8 68.8
IMAGE=stream_02.bmp
ID=Stream_Morph_02
SCALE=0.250000

LM=8
50.0 10.2
39.1 28.1
36.2 62.1
44.1 98.2
50.1 118.1
56.0 97.9
63.9 61.9
61.1 28.0
IMAGE=stream_03_rotated.bmp
ID=Stream_Morph_03
SCALE=0.250000

LM=8
10.0 50.0
27.0 67.5
60.0 74.0
96.0 60.0
115.0 50.0
96.0 40.0
60.0 26.0
27.0 32.5
IMAGE=lake_01.bmp
ID=Lake_Morph_01
SCALE=0.250000

LM=8
15.0 45.0
33.0 63.2
67.0 70.1
103.5 55.4
123.0 45.0
103.2 34.8
66.8 20.2
32.8 27.1
IMAGE=lake_02.bmp
ID=Lake_Morph_02
SCALE=0.250000

LM=8
20.0 60.0
55.0 96.0
122.0 109.0
193.0 81.0
231.0 60.0
193.0 39.0
122.0 11.0
55.0 24.0
IMAGE=lake_03_2x_zoom.bmp
ID=Lake_Morph_03
SCALE=0.125000
"
}

/// Runs Rohlf & Slice (1990) 2D Generalized Procrustes Analysis (`GPA`) + Jacobi PCA across `specimens`.
pub fn run_generalized_procrustes_2d(
    specimens: &[LandmarkSpecimen2D],
) -> Result<GpaAnalysis2D, String> {
    if specimens.len() < 2 {
        return Err(String::from("GPA requires at least 2 landmarked specimens."));
    }
    let k = specimens[0].landmarks.len();
    if k < 3 {
        return Err(String::from("GPA requires at least 3 homologous 2D landmarks per specimen."));
    }
    for sp in specimens {
        if sp.landmarks.len() != k {
            return Err(format!(
                "Specimen '{}' has {} landmarks; expected {k}.",
                sp.specimen_id,
                sp.landmarks.len()
            ));
        }
    }

    let m = specimens.len();
    let mut centroid_sizes = Vec::with_capacity(m);
    let mut shapes: Vec<Vec<[f64; 2]>> = Vec::with_capacity(m);

    // Step 1: Translate each configuration to origin and normalize to unit Centroid Size (`CS = 1`)
    for sp in specimens {
        let scale = sp.scale_factor.unwrap_or(1.0);
        let mx = sp.landmarks.iter().map(|p| p[0]).sum::<f64>() / (k as f64);
        let my = sp.landmarks.iter().map(|p| p[1]).sum::<f64>() / (k as f64);

        let mut centered = Vec::with_capacity(k);
        let mut ss = 0.0_f64;
        for &[x, y] in &sp.landmarks {
            let cx = (x - mx) * scale;
            let cy = (y - my) * scale;
            ss += cx * cx + cy * cy;
            centered.push([cx, cy]);
        }
        let cs = ss.sqrt().max(1e-12);
        for p in &mut centered {
            p[0] /= cs;
            p[1] /= cs;
        }
        centroid_sizes.push(cs);
        shapes.push(centered);
    }

    // Step 2: Iterative orthogonal planar Procrustes alignment to the evolving consensus shape
    let mut consensus = shapes[0].clone();
    for _iter in 0..50 {
        for s in &mut shapes {
            *s = rotate_to_reference_2d(s, &consensus);
        }
        let mut next_mean = vec![[0.0_f64; 2]; k];
        for s in &shapes {
            for (j, p) in s.iter().enumerate() {
                next_mean[j][0] += p[0] / (m as f64);
                next_mean[j][1] += p[1] / (m as f64);
            }
        }
        // Renormalize consensus to unit Centroid Size
        let norm = next_mean
            .iter()
            .map(|p| p[0] * p[0] + p[1] * p[1])
            .sum::<f64>()
            .sqrt()
            .max(1e-12);
        let mut delta = 0.0_f64;
        for j in 0..k {
            next_mean[j][0] /= norm;
            next_mean[j][1] /= norm;
            delta += (next_mean[j][0] - consensus[j][0]).powi(2)
                + (next_mean[j][1] - consensus[j][1]).powi(2);
        }
        consensus = next_mean;
        if delta < 1e-20 {
            break;
        }
    }

    let dists: Vec<f64> = shapes
        .iter()
        .map(|s| procrustes_distance_2d(s, &consensus))
        .collect();

    // Step 3: Tangent-space covariance matrix (`dim = 2 * k`) + cyclic Jacobi PCA
    let dim = 2 * k;
    let mut resid = vec![vec![0.0_f64; dim]; m];
    for (i, s) in shapes.iter().enumerate() {
        for j in 0..k {
            resid[i][2 * j] = s[j][0] - consensus[j][0];
            resid[i][2 * j + 1] = s[j][1] - consensus[j][1];
        }
    }

    let denom = (m.saturating_sub(1).max(1)) as f64;
    let mut cov = vec![vec![0.0_f64; dim]; dim];
    for a in 0..dim {
        for b in a..dim {
            let sum: f64 = (0..m).map(|i| resid[i][a] * resid[i][b]).sum();
            let v = sum / denom;
            cov[a][b] = v;
            cov[b][a] = v;
        }
    }

    let (eigvals, eigvecs) = jacobi_eigen_symmetric(&cov);
    let tot_var: f64 = eigvals.iter().map(|v| v.max(0.0)).sum::<f64>().max(1e-18);
    let pc1_pct = (eigvals[0].max(0.0) / tot_var) * 100.0;
    let pc2_pct = (eigvals.get(1).copied().unwrap_or(0.0).max(0.0) / tot_var) * 100.0;

    let mut pc12 = Vec::with_capacity(m);
    for row in &resid {
        let s1: f64 = (0..dim).map(|d| row[d] * eigvecs[d][0]).sum();
        let s2: f64 = if dim > 1 {
            (0..dim).map(|d| row[d] * eigvecs[d][1]).sum()
        } else {
            0.0
        };
        pc12.push([s1, s2]);
    }

    Ok(GpaAnalysis2D {
        specimen_ids: specimens.iter().map(|s| s.specimen_id.clone()).collect(),
        num_landmarks: k,
        centroid_sizes,
        consensus_shape: consensus,
        aligned_shapes: shapes,
        procrustes_distances_to_mean: dists,
        pc12_scores: pc12,
        pc1_variance_pct: pc1_pct,
        pc2_variance_pct: pc2_pct,
    })
}

/// Optimal proper 2D orthogonal rotation (`SO(2)`, no reflection) aligning `target` onto `reference`.
pub fn rotate_to_reference_2d(target: &[[f64; 2]], reference: &[[f64; 2]]) -> Vec<[f64; 2]> {
    let mut dot_sum = 0.0_f64;
    let mut cross_sum = 0.0_f64;
    for (t, r) in target.iter().zip(reference.iter()) {
        dot_sum += t[0] * r[0] + t[1] * r[1];
        cross_sum += t[0] * r[1] - t[1] * r[0];
    }
    let norm = (dot_sum * dot_sum + cross_sum * cross_sum).sqrt().max(1e-15);
    let cos_t = dot_sum / norm;
    let sin_t = cross_sum / norm;

    target
        .iter()
        .map(|p| [p[0] * cos_t - p[1] * sin_t, p[0] * sin_t + p[1] * cos_t])
        .collect()
}

pub fn procrustes_distance_2d(a: &[[f64; 2]], b: &[[f64; 2]]) -> f64 {
    a.iter()
        .zip(b.iter())
        .map(|(p, q)| (p[0] - q[0]).powi(2) + (p[1] - q[1]).powi(2))
        .sum::<f64>()
        .sqrt()
}

fn jacobi_eigen_symmetric(mat: &[Vec<f64>]) -> (Vec<f64>, Vec<Vec<f64>>) {
    let n = mat.len();
    let mut a = mat.to_vec();
    let mut v = vec![vec![0.0_f64; n]; n];
    for (i, row) in v.iter_mut().enumerate() {
        row[i] = 1.0;
    }

    for _sweep in 0..80 {
        let mut max_off = 0.0_f64;
        let (mut p, mut q) = (0usize, 1usize.min(n.saturating_sub(1)));
        for (i, row) in a.iter().enumerate() {
            for (j, &val) in row.iter().enumerate().skip(i + 1) {
                if val.abs() > max_off {
                    max_off = val.abs();
                    p = i;
                    q = j;
                }
            }
        }
        if max_off < 1e-14 || p == q {
            break;
        }

        let app = a[p][p];
        let aqq = a[q][q];
        let apq = a[p][q];
        let tau = (aqq - app) / (2.0 * apq);
        let t = if tau >= 0.0 {
            1.0 / (tau + (1.0 + tau * tau).sqrt())
        } else {
            -1.0 / (-tau + (1.0 + tau * tau).sqrt())
        };
        let c = 1.0 / (1.0 + t * t).sqrt();
        let s = t * c;

        a[p][p] = app - t * apq;
        a[q][q] = aqq + t * apq;
        a[p][q] = 0.0;
        a[q][p] = 0.0;

        for r in 0..n {
            if r != p && r != q {
                let arp = a[r][p];
                let arq = a[r][q];
                a[r][p] = c * arp - s * arq;
                a[p][r] = a[r][p];
                a[r][q] = s * arp + c * arq;
                a[q][r] = a[r][q];
            }
            let vrp = v[r][p];
            let vrq = v[r][q];
            v[r][p] = c * vrp - s * vrq;
            v[r][q] = s * vrp + c * vrq;
        }
    }

    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| a[j][j].partial_cmp(&a[i][i]).unwrap_or(std::cmp::Ordering::Equal));

    let vals: Vec<f64> = order.iter().map(|&idx| a[idx][idx]).collect();
    let mut vecs = vec![vec![0.0_f64; n]; n];
    for r in 0..n {
        for (col_out, &col_src) in order.iter().enumerate() {
            vecs[r][col_out] = v[r][col_src];
        }
    }
    (vals, vecs)
}

/// Renders a two-panel publication `.svg` displaying (Left) Superimposed Procrustes Landmarks &
/// Consensus Wireframe and (Right) `PC1` vs `PC2` Tangent Morphospace Scatter Plot.
pub fn render_gpa_morphospace_svg(gpa: &GpaAnalysis2D) -> String {
    let width = 880.0_f64;
    let height = 440.0_f64;
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
        "  <text x=\"28\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">Generalized Procrustes Superimposition &amp; Morphospace PCA ({} specimens, {} landmarks)</text>",
        gpa.specimen_ids.len(),
        gpa.num_landmarks
    );

    // Panel A: Procrustes superimposition (`x = 30..410, y = 50..400`)
    let _ = writeln!(
        svg,
        "  <rect x=\"30\" y=\"50\" width=\"380\" height=\"350\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let map_proc = |p: [f64; 2]| -> (f64, f64) { (220.0 + p[0] * 260.0, 225.0 - p[1] * 260.0) };

    for shape in &gpa.aligned_shapes {
        for &pt in shape {
            let (px, py) = map_proc(pt);
            let _ = writeln!(
                svg,
                "  <circle cx=\"{px:.2}\" cy=\"{py:.2}\" r=\"3.2\" fill=\"#93c5fd\" fill-opacity=\"0.85\"/>"
            );
        }
    }

    let mut poly = String::new();
    for &pt in &gpa.consensus_shape {
        let (px, py) = map_proc(pt);
        let _ = write!(poly, "{px:.2},{py:.2} ");
    }
    let _ = writeln!(
        svg,
        "  <polygon points=\"{}\" fill=\"none\" stroke=\"#0f172a\" stroke-width=\"2.0\"/>",
        poly.trim()
    );
    for (idx, &pt) in gpa.consensus_shape.iter().enumerate() {
        let (px, py) = map_proc(pt);
        let _ = writeln!(
            svg,
            "  <circle cx=\"{px:.2}\" cy=\"{py:.2}\" r=\"4.5\" fill=\"#dc2626\"/>"
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10\" fill=\"#0f172a\">LM{}</text>",
            px + 6.0,
            py - 4.0,
            idx + 1
        );
    }

    // Panel B: PC1 vs PC2 (`x = 460..850, y = 50..400`)
    let _ = writeln!(
        svg,
        "  <rect x=\"460\" y=\"50\" width=\"390\" height=\"350\" fill=\"#f8fafc\" stroke=\"#cbd5e1\"/>"
    );
    let max_pc = gpa
        .pc12_scores
        .iter()
        .flat_map(|s| [s[0].abs(), s[1].abs()])
        .fold(0.05_f64, f64::max)
        * 1.25;

    for (sp_id, &[pc1, pc2]) in gpa.specimen_ids.iter().zip(&gpa.pc12_scores) {
        let sx = 655.0 + (pc1 / max_pc) * 165.0;
        let sy = 225.0 - (pc2 / max_pc) * 145.0;
        let _ = writeln!(
            svg,
            "  <circle cx=\"{sx:.2}\" cy=\"{sy:.2}\" r=\"5.2\" fill=\"#2563eb\"/>"
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"sans-serif\" font-size=\"10.5\" fill=\"#1e293b\">{}</text>",
            sx + 7.0,
            sy + 3.5,
            xml_escape(sp_id)
        );
    }

    let _ = writeln!(
        svg,
        "  <text x=\"655\" y=\"422\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"11.5\" fill=\"#334155\">PC1 ({:.1}% var) vs PC2 ({:.1}% var)</text>",
        gpa.pc1_variance_pct, gpa.pc2_variance_pct
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
    fn verifies_procrustes_rotation_scale_invariance_and_pca_separation() {
        let specs = parse_tps_landmarks(demo_salmonid_landmarks_tps()).unwrap();
        assert_eq!(specs.len(), 6);

        let gpa = run_generalized_procrustes_2d(&specs).unwrap();
        // Specimen 3 (`Stream_Morph_03`) was rotated 90 deg from `Stream_Morph_01`; Procrustes distance must be near 0!
        let d_stream_rot = procrustes_distance_2d(&gpa.aligned_shapes[0], &gpa.aligned_shapes[2]);
        assert!(d_stream_rot < 0.02);

        // Specimen 6 (`Lake_Morph_03`) was 2x zoomed with `SCALE=0.125` (half of 0.25); Centroid Size must match Lake_01!
        assert!((gpa.centroid_sizes[3] - gpa.centroid_sizes[5]).abs() < 1.0);

        // PC1 should capture the major Stream vs Lake body-depth axis (> 85% variance)
        assert!(gpa.pc1_variance_pct > 85.0);
    }
}
