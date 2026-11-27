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

//! Tool #7 (`osabt forge3d`) — Watertight 3D triangle-mesh engine (`TriMesh`), Surface Nets
//! scalar-field isosurface mesher, Taubin `lambda/mu` volume-preserving smoother, manifold
//! half-edge verifier, filament mass calculator, and Binary `.stl`, ASCII `.stl`, Wavefront `.obj`,
//! and native ZIP-packaged 3D Manufacturing Format (`.3mf`) exporters.

use std::collections::HashMap;
use std::fmt::Write as _;

use crate::svg_util::xml_escape;

/// Indexed 3D triangle mesh in millimeter print coordinates.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct TriMesh {
    pub name: String,
    pub vertices: Vec<[f32; 3]>,
    pub triangles: Vec<[u32; 3]>,
}

/// Summary geometry, manifold health, and 3D-printer material metrics for a [`TriMesh`].
#[derive(Debug, Clone)]
pub struct MeshPrintReport {
    pub vertex_count: usize,
    pub triangle_count: usize,
    pub bounding_min_mm: [f32; 3],
    pub bounding_max_mm: [f32; 3],
    pub dimensions_mm: [f32; 3],
    pub signed_volume_mm3: f64,
    pub surface_area_mm2: f64,
    /// True iff every directed edge `(u, v)` appears once and has an exact opposite `(v, u)`.
    pub is_watertight_manifold: bool,
    pub boundary_or_nonmanifold_edges: usize,
    pub euler_characteristic: i64,
    /// Estimated filament/resin mass assuming 20% FDM infill + 2 perimeters (`PLA ~ 1.24 g/cm^3`) vs solid SLA resin (`1.12 g/cm^3`).
    pub est_fdm_pla_grams: f64,
    pub est_sla_resin_grams: f64,
}

/// Dense 3D scalar field `values[z * ny * nx + y * nx + x]` sampled on a regular grid.
#[derive(Debug, Clone)]
pub struct ScalarGrid3D {
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
    pub origin: [f32; 3],
    pub voxel_size: f32,
    pub values: Vec<f32>,
}

impl ScalarGrid3D {
    pub fn new(nx: usize, ny: usize, nz: usize, origin: [f32; 3], voxel_size: f32, fill: f32) -> Self {
        Self {
            nx,
            ny,
            nz,
            origin,
            voxel_size,
            values: vec![fill; nx * ny * nz],
        }
    }

    #[inline]
    pub fn idx(&self, x: usize, y: usize, z: usize) -> usize {
        (z * self.ny + y) * self.nx + x
    }

    #[inline]
    pub fn get(&self, x: usize, y: usize, z: usize) -> f32 {
        self.values[self.idx(x, y, z)]
    }

    #[inline]
    pub fn set(&mut self, x: usize, y: usize, z: usize, v: f32) {
        let i = self.idx(x, y, z);
        self.values[i] = v;
    }

    #[inline]
    pub fn world_pos(&self, x: f32, y: f32, z: f32) -> [f32; 3] {
        [
            self.origin[0] + x * self.voxel_size,
            self.origin[1] + y * self.voxel_size,
            self.origin[2] + z * self.voxel_size,
        ]
    }
}

impl TriMesh {
    pub fn new(name: impl Into<String>) -> Self {
        Self {
            name: name.into(),
            vertices: Vec::new(),
            triangles: Vec::new(),
        }
    }

    /// Appends all vertices and faces of `other` into `self`.
    pub fn append_mesh(&mut self, other: &TriMesh) {
        let offset = self.vertices.len() as u32;
        self.vertices.extend_from_slice(&other.vertices);
        for &[a, b, c] in &other.triangles {
            self.triangles.push([a + offset, b + offset, c + offset]);
        }
    }

    /// Centers XY at `(0, 0)`, seats `Z_min = 0.0` flat on the print bed, and scales the longest
    /// axis to `target_max_extent_mm`.
    pub fn normalize_to_print_bed(&mut self, target_max_extent_mm: f32) {
        if self.vertices.is_empty() || target_max_extent_mm <= 0.0 {
            return;
        }
        let (bmin, bmax) = self.bounds();
        let dx = (bmax[0] - bmin[0]).max(1e-6);
        let dy = (bmax[1] - bmin[1]).max(1e-6);
        let dz = (bmax[2] - bmin[2]).max(1e-6);
        let max_span = dx.max(dy).max(dz);
        let scale = target_max_extent_mm / max_span;

        let cx = 0.5 * (bmin[0] + bmax[0]);
        let cy = 0.5 * (bmin[1] + bmax[1]);
        let z0 = bmin[2];

        for v in &mut self.vertices {
            v[0] = (v[0] - cx) * scale;
            v[1] = (v[1] - cy) * scale;
            v[2] = (v[2] - z0) * scale;
        }
    }

    pub fn bounds(&self) -> ([f32; 3], [f32; 3]) {
        if self.vertices.is_empty() {
            return ([0.0; 3], [0.0; 3]);
        }
        let mut mn = [f32::INFINITY; 3];
        let mut mx = [f32::NEG_INFINITY; 3];
        for v in &self.vertices {
            for k in 0..3 {
                mn[k] = mn[k].min(v[k]);
                mx[k] = mx[k].max(v[k]);
            }
        }
        (mn, mx)
    }

    /// Applies Taubin (`lambda = +0.50, mu = -0.53`) two-step Laplacian smoothing, removing
    /// voxel staircase faceting without shrinking biological structures.
    pub fn taubin_smooth(&mut self, iterations: usize) {
        if self.vertices.is_empty() || self.triangles.is_empty() {
            return;
        }
        let n = self.vertices.len();
        let mut adj: Vec<Vec<usize>> = vec![Vec::new(); n];
        for &[a, b, c] in &self.triangles {
            let (ia, ib, ic) = (a as usize, b as usize, c as usize);
            for (u, v) in [(ia, ib), (ib, ic), (ic, ia)] {
                if !adj[u].contains(&v) {
                    adj[u].push(v);
                }
                if !adj[v].contains(&u) {
                    adj[v].push(u);
                }
            }
        }

        for _ in 0..iterations {
            self.laplacian_step(&adj, 0.50);
            self.laplacian_step(&adj, -0.53);
        }
    }

    fn laplacian_step(&mut self, adj: &[Vec<usize>], factor: f32) {
        let prev = self.vertices.clone();
        for (i, nbrs) in adj.iter().enumerate() {
            if nbrs.is_empty() {
                continue;
            }
            let mut sum = [0.0_f32; 3];
            for &nb in nbrs {
                sum[0] += prev[nb][0];
                sum[1] += prev[nb][1];
                sum[2] += prev[nb][2];
            }
            let inv = 1.0 / (nbrs.len() as f32);
            for k in 0..3 {
                let avg = sum[k] * inv;
                self.vertices[i][k] = prev[i][k] + factor * (avg - prev[i][k]);
            }
        }
    }

    /// Computes exact signed divergence volume, surface area, half-edge watertightness, Euler
    /// characteristic `V - E + F`, and PLA/SLA print mass estimates.
    pub fn inspect_printability(&self) -> MeshPrintReport {
        let (bmin, bmax) = self.bounds();
        let dims = [bmax[0] - bmin[0], bmax[1] - bmin[1], bmax[2] - bmin[2]];

        let mut vol6 = 0.0_f64;
        let mut area = 0.0_f64;
        let mut half_edges: HashMap<(u32, u32), usize> =
            HashMap::with_capacity(self.triangles.len() * 3);

        for &[i0, i1, i2] in &self.triangles {
            let v0 = to_f64(self.vertices[i0 as usize]);
            let v1 = to_f64(self.vertices[i1 as usize]);
            let v2 = to_f64(self.vertices[i2 as usize]);

            let c = cross(v1, v2);
            vol6 += dot(v0, c);

            let e1 = sub(v1, v0);
            let e2 = sub(v2, v0);
            let n = cross(e1, e2);
            area += 0.5 * dot(n, n).sqrt();

            for (u, v) in [(i0, i1), (i1, i2), (i2, i0)] {
                *half_edges.entry((u, v)).or_insert(0) += 1;
            }
        }

        let mut bad_edges = 0usize;
        let mut undirected_edges = 0usize;
        for (&(u, v), &cnt) in &half_edges {
            let rev_cnt = half_edges.get(&(v, u)).copied().unwrap_or(0);
            if cnt != 1 || rev_cnt != 1 {
                bad_edges += 1;
            }
            if u < v {
                undirected_edges += 1;
            } else if rev_cnt == 0 {
                undirected_edges += 1;
            }
        }

        let signed_vol_mm3 = vol6 / 6.0;
        let abs_cm3 = signed_vol_mm3.abs() / 1000.0;
        let shell_cm3 = (area * 0.8) / 1000.0;
        let fdm_eff_cm3 =
            (shell_cm3 + (abs_cm3 - shell_cm3).max(0.0) * 0.20).min(abs_cm3.max(shell_cm3));

        let euler = (self.vertices.len() as i64) - (undirected_edges as i64)
            + (self.triangles.len() as i64);

        MeshPrintReport {
            vertex_count: self.vertices.len(),
            triangle_count: self.triangles.len(),
            bounding_min_mm: bmin,
            bounding_max_mm: bmax,
            dimensions_mm: dims,
            signed_volume_mm3: signed_vol_mm3,
            surface_area_mm2: area,
            is_watertight_manifold: !self.triangles.is_empty() && bad_edges == 0,
            boundary_or_nonmanifold_edges: bad_edges,
            euler_characteristic: euler,
            est_fdm_pla_grams: fdm_eff_cm3 * 1.24,
            est_sla_resin_grams: abs_cm3 * 1.12,
        }
    }

    /// Serializes the mesh to standard IEEE-754 little-endian Binary `.stl` (`84 + 50 * N` bytes).
    pub fn to_binary_stl(&self) -> Vec<u8> {
        let mut out = Vec::with_capacity(84 + self.triangles.len() * 50);
        let mut header = [b' '; 80];
        let title = format!("OSABT forge3d (GPL-3.0) - {} (mm)", self.name);
        let hb = title.as_bytes();
        let copy_len = hb.len().min(80);
        header[..copy_len].copy_from_slice(&hb[..copy_len]);
        out.extend_from_slice(&header);
        out.extend_from_slice(&(self.triangles.len() as u32).to_le_bytes());

        for &[i0, i1, i2] in &self.triangles {
            let v0 = self.vertices[i0 as usize];
            let v1 = self.vertices[i1 as usize];
            let v2 = self.vertices[i2 as usize];
            let n = unit_normal_f32(v0, v1, v2);
            for comp in n {
                out.extend_from_slice(&comp.to_le_bytes());
            }
            for v in [v0, v1, v2] {
                for comp in v {
                    out.extend_from_slice(&comp.to_le_bytes());
                }
            }
            out.extend_from_slice(&0u16.to_le_bytes());
        }
        out
    }

    /// Serializes the mesh to human-readable ASCII `.stl`.
    pub fn to_ascii_stl(&self) -> String {
        let mut s = String::with_capacity(self.triangles.len() * 180 + 64);
        let safe_name: String = self
            .name
            .chars()
            .map(|c| if c.is_ascii_alphanumeric() { c } else { '_' })
            .collect();
        let _ = writeln!(s, "solid {safe_name}");
        for &[i0, i1, i2] in &self.triangles {
            let v0 = self.vertices[i0 as usize];
            let v1 = self.vertices[i1 as usize];
            let v2 = self.vertices[i2 as usize];
            let n = unit_normal_f32(v0, v1, v2);
            let _ = writeln!(s, "  facet normal {:.6e} {:.6e} {:.6e}", n[0], n[1], n[2]);
            s.push_str("    outer loop\n");
            for v in [v0, v1, v2] {
                let _ = writeln!(s, "      vertex {:.6e} {:.6e} {:.6e}", v[0], v[1], v[2]);
            }
            s.push_str("    endloop\n  endfacet\n");
        }
        let _ = writeln!(s, "endsolid {safe_name}");
        s
    }

    /// Serializes the mesh to Wavefront `.obj` (1-indexed vertices).
    pub fn to_wavefront_obj(&self) -> String {
        let mut s = String::with_capacity(self.vertices.len() * 36 + self.triangles.len() * 24);
        let _ = writeln!(s, "# OSABT forge3d (GPL-3.0) — {}", self.name);
        for v in &self.vertices {
            let _ = writeln!(s, "v {:.5} {:.5} {:.5}", v[0], v[1], v[2]);
        }
        for &[a, b, c] in &self.triangles {
            let _ = writeln!(s, "f {} {} {}", a + 1, b + 1, c + 1);
        }
        s
    }

    /// Packages a genuine Open Packaging Conventions `.3mf` archive (ZIP `STORE` container with
    /// IEEE 802.3 CRC-32) containing `[Content_Types].xml`, `_rels/.rels`, and `3D/3dmodel.model`
    /// conforming to the 3MF Core Specification (`http://schemas.microsoft.com/3dmanufacturing/core/2015/02`).
    pub fn to_3mf_package(&self) -> Vec<u8> {
        let content_types = b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<Types xmlns=\"http://schemas.openxmlformats.org/package/2006/content-types\">\n\
  <Default Extension=\"rels\" ContentType=\"application/vnd.openxmlformats-package.relationships+xml\"/>\n\
  <Default Extension=\"model\" ContentType=\"application/vnd.ms-package.3dmanufacturing-3dmodel+xml\"/>\n\
</Types>\n";

        let rels = b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n\
<Relationships xmlns=\"http://schemas.openxmlformats.org/package/2006/relationships\">\n\
  <Relationship Target=\"/3D/3dmodel.model\" Id=\"rel0\" Type=\"http://schemas.microsoft.com/3dmanufacturing/2013/01/3dmodel\"/>\n\
</Relationships>\n";

        let mut model_xml =
            String::with_capacity(self.vertices.len() * 48 + self.triangles.len() * 48 + 512);
        model_xml.push_str("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
        model_xml.push_str("<model unit=\"millimeter\" xml:lang=\"en-US\" xmlns=\"http://schemas.microsoft.com/3dmanufacturing/core/2015/02\">\n");
        let _ = writeln!(
            model_xml,
            "  <metadata name=\"Title\">{}</metadata>",
            xml_escape(&self.name)
        );
        model_xml.push_str(
            "  <resources>\n    <object id=\"1\" type=\"model\">\n      <mesh>\n        <vertices>\n",
        );
        for v in &self.vertices {
            let _ = writeln!(
                model_xml,
                "          <vertex x=\"{:.5}\" y=\"{:.5}\" z=\"{:.5}\"/>",
                v[0], v[1], v[2]
            );
        }
        model_xml.push_str("        </vertices>\n        <triangles>\n");
        for &[v1, v2, v3] in &self.triangles {
            let _ = writeln!(
                model_xml,
                "          <triangle v1=\"{v1}\" v2=\"{v2}\" v3=\"{v3}\"/>"
            );
        }
        model_xml.push_str("        </triangles>\n      </mesh>\n    </object>\n  </resources>\n");
        model_xml.push_str("  <build>\n    <item objectid=\"1\"/>\n  </build>\n</model>\n");

        build_uncompressed_zip(&[
            ("[Content_Types].xml", content_types),
            ("_rels/.rels", rels),
            ("3D/3dmodel.model", model_xml.as_bytes()),
        ])
    }
}

/// Parses a Binary or ASCII `.stl` byte buffer into a deduplicated [`TriMesh`] for inspection.
pub fn parse_stl_bytes(name: &str, bytes: &[u8]) -> Result<TriMesh, String> {
    let bytes_nobom = bytes.strip_prefix(b"\xEF\xBB\xBF").unwrap_or(bytes);
    if bytes.len() >= 84 && !bytes_nobom.starts_with(b"solid ") {
        let n_tri = u32::from_le_bytes([bytes[80], bytes[81], bytes[82], bytes[83]]) as usize;
        if 84 + n_tri * 50 <= bytes.len() {
            let mut mesh = TriMesh::new(name);
            let mut vmap: HashMap<[u32; 3], u32> = HashMap::new();
            for i in 0..n_tri {
                let off = 84 + i * 50 + 12;
                let mut tri = [0u32; 3];
                for (vi, slot) in tri.iter_mut().enumerate() {
                    let voff = off + vi * 12;
                    let x = f32::from_le_bytes(
                        bytes
                            .get(voff..voff + 4)
                            .and_then(|s| s.try_into().ok())
                            .ok_or_else(|| String::from("Truncated binary STL vertex X."))?,
                    );
                    let y = f32::from_le_bytes(
                        bytes
                            .get(voff + 4..voff + 8)
                            .and_then(|s| s.try_into().ok())
                            .ok_or_else(|| String::from("Truncated binary STL vertex Y."))?,
                    );
                    let z = f32::from_le_bytes(
                        bytes
                            .get(voff + 8..voff + 12)
                            .and_then(|s| s.try_into().ok())
                            .ok_or_else(|| String::from("Truncated binary STL vertex Z."))?,
                    );
                    let key = [x.to_bits(), y.to_bits(), z.to_bits()];
                    let idx = *vmap.entry(key).or_insert_with(|| {
                        let id = mesh.vertices.len() as u32;
                        mesh.vertices.push([x, y, z]);
                        id
                    });
                    *slot = idx;
                }
                mesh.triangles.push(tri);
            }
            return Ok(mesh);
        }
    }

    let text = std::str::from_utf8(bytes).map_err(|_| String::from("Invalid STL file."))?;
    let mut mesh = TriMesh::new(name);
    let mut vmap: HashMap<[u32; 3], u32> = HashMap::new();
    let mut cur_vs: Vec<u32> = Vec::with_capacity(3);
    for line in text.trim_start_matches('\u{feff}').lines() {
        let t = line.trim();
        if let Some(rest) = t.strip_prefix("vertex") {
            let nums: Vec<f32> = rest
                .split_whitespace()
                .filter_map(|s| s.parse::<f32>().ok())
                .collect();
            if nums.len() == 3 {
                let key = [nums[0].to_bits(), nums[1].to_bits(), nums[2].to_bits()];
                let idx = *vmap.entry(key).or_insert_with(|| {
                    let id = mesh.vertices.len() as u32;
                    mesh.vertices.push([nums[0], nums[1], nums[2]]);
                    id
                });
                cur_vs.push(idx);
                if cur_vs.len() == 3 {
                    mesh.triangles.push([cur_vs[0], cur_vs[1], cur_vs[2]]);
                    cur_vs.clear();
                }
            }
        }
    }
    if mesh.triangles.is_empty() {
        return Err(String::from("No triangles found in STL buffer."));
    }
    Ok(mesh)
}

/// Extracts a guaranteed watertight, outward-oriented 2-manifold isosurface at `iso_level`
/// from `grid` (`inside` defined as `value >= iso_level`) via Naive Surface Nets (Gibson 1998).
pub fn extract_isosurface_surface_nets(
    name: &str,
    grid: &ScalarGrid3D,
    iso_level: f32,
) -> TriMesh {
    let mut mesh = TriMesh::new(name);
    if grid.nx < 3 || grid.ny < 3 || grid.nz < 3 {
        return mesh;
    }

    const CUBE_EDGES: [(usize, usize); 12] = [
        (0, 1),
        (2, 3),
        (4, 5),
        (6, 7),
        (0, 2),
        (1, 3),
        (4, 6),
        (5, 7),
        (0, 4),
        (1, 5),
        (2, 6),
        (3, 7),
    ];

    let cx_max = grid.nx - 1;
    let cy_max = grid.ny - 1;
    let cz_max = grid.nz - 1;
    let cell_idx = |cx: usize, cy: usize, cz: usize| -> usize { (cz * cy_max + cy) * cx_max + cx };

    // Keep outer grid boundary vertices strictly outside (`false`) so every crossing surface
    // closes cleanly inside the domain, and resolve any 2x2 face-checkerboard sign configurations
    // so Surface Nets emits a strictly watertight 2-manifold (`bad_edges == 0`).
    let mut ins: Vec<bool> = vec![false; grid.nx * grid.ny * grid.nz];
    for z in 1..(grid.nz - 1) {
        for y in 1..(grid.ny - 1) {
            for x in 1..(grid.nx - 1) {
                ins[grid.idx(x, y, z)] = grid.get(x, y, z) >= iso_level;
            }
        }
    }

    loop {
        let mut changed = false;
        // XY faces
        for z in 1..(grid.nz - 1) {
            for y in 1..(grid.ny - 2) {
                for x in 1..(grid.nx - 2) {
                    let i00 = ins[grid.idx(x, y, z)];
                    let i10 = ins[grid.idx(x + 1, y, z)];
                    let i01 = ins[grid.idx(x, y + 1, z)];
                    let i11 = ins[grid.idx(x + 1, y + 1, z)];
                    if i00 == i11 && i10 == i01 && i00 != i10 {
                        let (ax, ay) = if !i00 {
                            if grid.get(x, y, z) >= grid.get(x + 1, y + 1, z) {
                                (x, y)
                            } else {
                                (x + 1, y + 1)
                            }
                        } else if grid.get(x + 1, y, z) >= grid.get(x, y + 1, z) {
                            (x + 1, y)
                        } else {
                            (x, y + 1)
                        };
                        let idx = grid.idx(ax, ay, z);
                        ins[idx] = true;
                        changed = true;
                    }
                }
            }
        }
        // XZ faces
        for y in 1..(grid.ny - 1) {
            for z in 1..(grid.nz - 2) {
                for x in 1..(grid.nx - 2) {
                    let i00 = ins[grid.idx(x, y, z)];
                    let i10 = ins[grid.idx(x + 1, y, z)];
                    let i01 = ins[grid.idx(x, y, z + 1)];
                    let i11 = ins[grid.idx(x + 1, y, z + 1)];
                    if i00 == i11 && i10 == i01 && i00 != i10 {
                        let (ax, az) = if !i00 {
                            if grid.get(x, y, z) >= grid.get(x + 1, y, z + 1) {
                                (x, z)
                            } else {
                                (x + 1, z + 1)
                            }
                        } else if grid.get(x + 1, y, z) >= grid.get(x, y, z + 1) {
                            (x + 1, z)
                        } else {
                            (x, z + 1)
                        };
                        let idx = grid.idx(ax, y, az);
                        ins[idx] = true;
                        changed = true;
                    }
                }
            }
        }
        // YZ faces
        for x in 1..(grid.nx - 1) {
            for z in 1..(grid.nz - 2) {
                for y in 1..(grid.ny - 2) {
                    let i00 = ins[grid.idx(x, y, z)];
                    let i10 = ins[grid.idx(x, y + 1, z)];
                    let i01 = ins[grid.idx(x, y, z + 1)];
                    let i11 = ins[grid.idx(x, y + 1, z + 1)];
                    if i00 == i11 && i10 == i01 && i00 != i10 {
                        let (ay, az) = if !i00 {
                            if grid.get(x, y, z) >= grid.get(x, y + 1, z + 1) {
                                (y, z)
                            } else {
                                (y + 1, z + 1)
                            }
                        } else if grid.get(x, y + 1, z) >= grid.get(x, y, z + 1) {
                            (y + 1, z)
                        } else {
                            (y, z + 1)
                        };
                        let idx = grid.idx(x, ay, az);
                        ins[idx] = true;
                        changed = true;
                    }
                }
            }
        }
        if !changed {
            break;
        }
    }

    let mut cell_vert: Vec<u32> = vec![u32::MAX; cx_max * cy_max * cz_max];

    for cz in 0..cz_max {
        for cy in 0..cy_max {
            for cx in 0..cx_max {
                let mut corner_vals = [0.0_f32; 8];
                let mut corner_pos = [[0.0_f32; 3]; 8];
                let mut mask = 0u8;
                for c in 0..8 {
                    let dx = c & 1;
                    let dy = (c >> 1) & 1;
                    let dz = (c >> 2) & 1;
                    let val = grid.get(cx + dx, cy + dy, cz + dz);
                    corner_vals[c] = val;
                    corner_pos[c] = [(cx + dx) as f32, (cy + dy) as f32, (cz + dz) as f32];
                    if ins[grid.idx(cx + dx, cy + dy, cz + dz)] {
                        mask |= 1 << c;
                    }
                }
                if mask == 0 || mask == 0xFF {
                    continue;
                }

                let mut sum = [0.0_f32; 3];
                let mut crossings = 0usize;
                for &(a, b) in &CUBE_EDGES {
                    let in_a = (mask & (1 << a)) != 0;
                    let in_b = (mask & (1 << b)) != 0;
                    if in_a != in_b {
                        let va = corner_vals[a];
                        let vb = corner_vals[b];
                        let denom = vb - va;
                        let t = if denom.abs() > 1e-9 {
                            ((iso_level - va) / denom).clamp(0.05, 0.95)
                        } else {
                            0.5
                        };
                        for k in 0..3 {
                            sum[k] += corner_pos[a][k] + t * (corner_pos[b][k] - corner_pos[a][k]);
                        }
                        crossings += 1;
                    }
                }

                let inv = 1.0 / (crossings.max(1) as f32);
                let wpos = grid.world_pos(sum[0] * inv, sum[1] * inv, sum[2] * inv);
                let vid = mesh.vertices.len() as u32;
                mesh.vertices.push(wpos);
                cell_vert[cell_idx(cx, cy, cz)] = vid;
            }
        }
    }

    // X-directed edges
    for z in 1..cz_max {
        for y in 1..cy_max {
            for x in 0..cx_max {
                let in0 = ins[grid.idx(x, y, z)];
                let in1 = ins[grid.idx(x + 1, y, z)];
                if in0 != in1 {
                    let c0 = cell_vert[cell_idx(x, y - 1, z - 1)];
                    let c1 = cell_vert[cell_idx(x, y, z - 1)];
                    let c2 = cell_vert[cell_idx(x, y, z)];
                    let c3 = cell_vert[cell_idx(x, y - 1, z)];
                    emit_quad(&mut mesh.triangles, c0, c1, c2, c3, in0);
                }
            }
        }
    }

    // Y-directed edges
    for z in 1..cz_max {
        for y in 0..cy_max {
            for x in 1..cx_max {
                let in0 = ins[grid.idx(x, y, z)];
                let in1 = ins[grid.idx(x, y + 1, z)];
                if in0 != in1 {
                    let c0 = cell_vert[cell_idx(x - 1, y, z - 1)];
                    let c1 = cell_vert[cell_idx(x - 1, y, z)];
                    let c2 = cell_vert[cell_idx(x, y, z)];
                    let c3 = cell_vert[cell_idx(x, y, z - 1)];
                    emit_quad(&mut mesh.triangles, c0, c1, c2, c3, in0);
                }
            }
        }
    }

    // Z-directed edges
    for z in 0..cz_max {
        for y in 1..cy_max {
            for x in 1..cx_max {
                let in0 = ins[grid.idx(x, y, z)];
                let in1 = ins[grid.idx(x, y, z + 1)];
                if in0 != in1 {
                    let c0 = cell_vert[cell_idx(x - 1, y - 1, z)];
                    let c1 = cell_vert[cell_idx(x, y - 1, z)];
                    let c2 = cell_vert[cell_idx(x, y, z)];
                    let c3 = cell_vert[cell_idx(x - 1, y, z)];
                    emit_quad(&mut mesh.triangles, c0, c1, c2, c3, in0);
                }
            }
        }
    }

    mesh
}

#[inline]
fn emit_quad(tris: &mut Vec<[u32; 3]>, c0: u32, c1: u32, c2: u32, c3: u32, positive_normal: bool) {
    if c0 == u32::MAX || c1 == u32::MAX || c2 == u32::MAX || c3 == u32::MAX {
        return;
    }
    if positive_normal {
        tris.push([c0, c1, c2]);
        tris.push([c0, c2, c3]);
    } else {
        tris.push([c0, c2, c1]);
        tris.push([c0, c3, c2]);
    }
}

fn unit_normal_f32(v0: [f32; 3], v1: [f32; 3], v2: [f32; 3]) -> [f32; 3] {
    let e1 = [v1[0] - v0[0], v1[1] - v0[1], v1[2] - v0[2]];
    let e2 = [v2[0] - v0[0], v2[1] - v0[1], v2[2] - v0[2]];
    let nx = e1[1] * e2[2] - e1[2] * e2[1];
    let ny = e1[2] * e2[0] - e1[0] * e2[2];
    let nz = e1[0] * e2[1] - e1[1] * e2[0];
    let len = (nx * nx + ny * ny + nz * nz).sqrt();
    if len > 1e-12 {
        [nx / len, ny / len, nz / len]
    } else {
        [0.0, 0.0, 1.0]
    }
}

fn to_f64(v: [f32; 3]) -> [f64; 3] {
    [v[0] as f64, v[1] as f64, v[2] as f64]
}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn crc32_ieee(bytes: &[u8]) -> u32 {
    let mut crc = 0xFFFF_FFFFu32;
    for &b in bytes {
        crc ^= b as u32;
        for _ in 0..8 {
            if (crc & 1) != 0 {
                crc = (crc >> 1) ^ 0xEDB8_8320;
            } else {
                crc >>= 1;
            }
        }
    }
    !crc
}

fn build_uncompressed_zip(files: &[(&str, &[u8])]) -> Vec<u8> {
    let mut archive = Vec::new();
    let mut central_dir = Vec::new();

    for &(name, data) in files {
        let local_offset = archive.len() as u32;
        let crc = crc32_ieee(data);
        let size = data.len() as u32;
        let name_bytes = name.as_bytes();
        let name_len = name_bytes.len() as u16;

        archive.extend_from_slice(&0x0403_4b50u32.to_le_bytes());
        archive.extend_from_slice(&20u16.to_le_bytes());
        archive.extend_from_slice(&0u16.to_le_bytes());
        archive.extend_from_slice(&0u16.to_le_bytes());
        archive.extend_from_slice(&0u16.to_le_bytes());
        archive.extend_from_slice(&0u16.to_le_bytes());
        archive.extend_from_slice(&crc.to_le_bytes());
        archive.extend_from_slice(&size.to_le_bytes());
        archive.extend_from_slice(&size.to_le_bytes());
        archive.extend_from_slice(&name_len.to_le_bytes());
        archive.extend_from_slice(&0u16.to_le_bytes());
        archive.extend_from_slice(name_bytes);
        archive.extend_from_slice(data);

        central_dir.extend_from_slice(&0x0201_4b50u32.to_le_bytes());
        central_dir.extend_from_slice(&20u16.to_le_bytes());
        central_dir.extend_from_slice(&20u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&crc.to_le_bytes());
        central_dir.extend_from_slice(&size.to_le_bytes());
        central_dir.extend_from_slice(&size.to_le_bytes());
        central_dir.extend_from_slice(&name_len.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u16.to_le_bytes());
        central_dir.extend_from_slice(&0u32.to_le_bytes());
        central_dir.extend_from_slice(&local_offset.to_le_bytes());
        central_dir.extend_from_slice(name_bytes);
    }

    let cd_offset = archive.len() as u32;
    let cd_size = central_dir.len() as u32;
    let count = files.len() as u16;
    archive.extend_from_slice(&central_dir);

    archive.extend_from_slice(&0x0605_4b50u32.to_le_bytes());
    archive.extend_from_slice(&0u16.to_le_bytes());
    archive.extend_from_slice(&0u16.to_le_bytes());
    archive.extend_from_slice(&count.to_le_bytes());
    archive.extend_from_slice(&count.to_le_bytes());
    archive.extend_from_slice(&cd_size.to_le_bytes());
    archive.extend_from_slice(&cd_offset.to_le_bytes());
    archive.extend_from_slice(&0u16.to_le_bytes());
    archive
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn extracts_watertight_positive_volume_sphere_and_exports_stl_3mf() {
        let mut grid = ScalarGrid3D::new(22, 22, 22, [-11.0, -11.0, -11.0], 1.0, 0.0);
        for z in 1..21 {
            for y in 1..21 {
                for x in 1..21 {
                    let p = grid.world_pos(x as f32, y as f32, z as f32);
                    let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                    grid.set(x, y, z, 8.0 - r);
                }
            }
        }
        let mut mesh = extract_isosurface_surface_nets("TestSphere", &grid, 0.0);
        mesh.taubin_smooth(2);
        let rep = mesh.inspect_printability();
        assert!(rep.is_watertight_manifold);
        assert_eq!(rep.boundary_or_nonmanifold_edges, 0);
        assert_eq!(rep.euler_characteristic, 2);
        assert!(rep.signed_volume_mm3 > 1500.0);

        let bin = mesh.to_binary_stl();
        assert_eq!(bin.len(), 84 + mesh.triangles.len() * 50);
        let parsed = parse_stl_bytes("Reloaded", &bin).unwrap();
        assert_eq!(parsed.triangles.len(), mesh.triangles.len());

        let pkg_3mf = mesh.to_3mf_package();
        assert_eq!(&pkg_3mf[0..4], b"PK\x03\x04");
    }
}
