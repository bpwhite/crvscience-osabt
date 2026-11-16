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

//! Tool #7 (`osabt forge3d`) — Biological 3D structure generators:
//! 1. wwPDB v3.3 (`ATOM`/`HETATM`) molecular parser with Bondi (1964) van der Waals atomic radii
//!    and watertight 3D-printable `BackboneTube` & `GaussianSurface` meshers.
//! 2. Gielis (2003) superformula silica diatom frustule generator (with marginal rim & areolae pores).
//! 3. Haeckel-inspired radiolarian cortical shell & icosahedral radial spicule generator.
//! 4. CRV Sim 3D voxel scalar/lattice importer (`x,y,z,density` CSV/TSV) to watertight manifold mesh.

use std::fmt::Write as _;

use crate::forge3d_mesh::{extract_isosurface_surface_nets, ScalarGrid3D, TriMesh};

/// Single parsed atom from an RCSB / wwPDB `ATOM` or `HETATM` record.
#[derive(Debug, Clone)]
pub struct PdbAtom {
    pub serial: u32,
    pub atom_name: String,
    pub res_name: String,
    pub chain_id: char,
    pub res_seq: i32,
    pub pos_angstrom: [f32; 3],
    pub element: String,
    pub vdw_radius_angstrom: f32,
}

/// 3D print representation styles for PDB macromolecular structures.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PdbPrintStyle {
    /// Smooth 3D-printable solvent-excluded Gaussian molecular surface across heavy atoms.
    GaussianSurface,
    /// Sturdy 3D-printable C-alpha / Phosphorus backbone worm tube.
    BackboneTube,
}

/// Returns the Bondi (1964) *J. Phys. Chem.* 68:441–451 van der Waals radius (in Å) for an element.
pub fn bondi_vdw_radius(element: &str) -> f32 {
    match element.trim().to_ascii_uppercase().as_str() {
        "H" => 1.20,
        "C" => 1.70,
        "N" => 1.55,
        "O" => 1.52,
        "F" => 1.47,
        "P" => 1.80,
        "S" => 1.80,
        "CL" => 1.75,
        "MG" => 1.73,
        "CA" => 2.31,
        "ZN" => 1.39,
        "FE" => 2.00,
        _ => 1.70,
    }
}

/// Parses wwPDB v3.3 `ATOM` and `HETATM` coordinate lines (excluding waters `HOH` unless requested).
/// Reads only the first NMR `MODEL` if multiple models are present.
pub fn parse_pdb_atoms(pdb_text: &str, include_hetatm: bool) -> Result<Vec<PdbAtom>, String> {
    let mut atoms = Vec::new();
    for line in pdb_text.lines() {
        if line.starts_with("ENDMDL") && !atoms.is_empty() {
            break;
        }
        let is_atom = line.starts_with("ATOM  ");
        let is_het = line.starts_with("HETATM");
        if !is_atom && !(include_hetatm && is_het) {
            continue;
        }
        if line.len() < 54 {
            continue;
        }
        let res_name = line.get(17..20).unwrap_or("UNK").trim().to_string();
        if res_name == "HOH" || res_name == "WAT" {
            continue;
        }
        let alt_loc = line.as_bytes().get(16).copied().unwrap_or(b' ') as char;
        if alt_loc != ' ' && alt_loc != 'A' {
            continue;
        }

        let serial = line
            .get(6..11)
            .and_then(|s| s.trim().parse::<u32>().ok())
            .unwrap_or((atoms.len() + 1) as u32);
        let atom_name = line.get(12..16).unwrap_or("CA").trim().to_string();
        let chain_id = line.as_bytes().get(21).copied().unwrap_or(b'A') as char;
        let res_seq = line
            .get(22..26)
            .and_then(|s| s.trim().parse::<i32>().ok())
            .unwrap_or(1);

        let x: f32 = line[30..38]
            .trim()
            .parse()
            .map_err(|_| format!("Invalid PDB X coordinate on line: {line}"))?;
        let y: f32 = line[38..46]
            .trim()
            .parse()
            .map_err(|_| format!("Invalid PDB Y coordinate on line: {line}"))?;
        let z: f32 = line[46..54]
            .trim()
            .parse()
            .map_err(|_| format!("Invalid PDB Z coordinate on line: {line}"))?;

        let elem = line
            .get(76..78)
            .map(str::trim)
            .filter(|s| !s.is_empty())
            .map(String::from)
            .unwrap_or_else(|| {
                atom_name
                    .chars()
                    .find(|c| c.is_ascii_alphabetic())
                    .map(|c| c.to_string())
                    .unwrap_or_else(|| String::from("C"))
            });

        atoms.push(PdbAtom {
            serial,
            atom_name,
            res_name,
            chain_id,
            res_seq,
            pos_angstrom: [x, y, z],
            vdw_radius_angstrom: bondi_vdw_radius(&elem),
            element: elem,
        });
    }

    if atoms.is_empty() {
        return Err(String::from("No valid ATOM records found in PDB input."));
    }
    Ok(atoms)
}

/// Built-in 16-residue canonical right-handed alpha-helix (`3.6` residues/turn, `1.5 Å` rise/res)
/// PDB structure (`pCRV_AlphaHelix`) for instant 3D molecular mesh generation.
pub fn demo_alpha_helix_pdb() -> String {
    let mut pdb = String::from("HEADER    DE NOVO PROTEIN                         11-NOV-26   CRV1\n");
    pdb.push_str("TITLE     OSABT BUILT-IN AMPHIPATHIC ALPHA HELIX DEMO\n");
    let residues = [
        "ALA", "LEU", "GLU", "LYS", "ALA", "ILE", "ALA", "GLN",
        "LEU", "GLU", "LYS", "ALA", "VAL", "ALA", "LYS", "GLY",
    ];
    let mut serial = 1u32;
    for (i, res) in residues.iter().enumerate() {
        let idx = i as f32;
        let theta = idx * (100.0_f32.to_radians());
        let z_ca = idx * 1.50;
        let ca = [2.30 * theta.cos(), 2.30 * theta.sin(), z_ca];
        let n = [1.55 * (theta - 0.45).cos(), 1.55 * (theta - 0.45).sin(), z_ca - 0.75];
        let c = [1.75 * (theta + 0.48).cos(), 1.75 * (theta + 0.48).sin(), z_ca + 0.72];
        let o = [2.05 * (theta + 0.52).cos(), 2.05 * (theta + 0.52).sin(), z_ca + 1.92];
        let cb = [3.65 * (theta + 0.12).cos(), 3.65 * (theta + 0.12).sin(), z_ca + 0.25];

        for (aname, el, p) in [
            ("N", "N", n),
            ("CA", "C", ca),
            ("C", "C", c),
            ("O", "O", o),
            ("CB", "C", cb),
        ] {
            if *res == "GLY" && aname == "CB" {
                continue;
            }
            let _ = writeln!(
                pdb,
                "ATOM  {:>5}  {:<3} {:>3} A{:>4}    {:>8.3}{:>8.3}{:>8.3}  1.00 20.00          {:>2}",
                serial,
                aname,
                res,
                i + 1,
                p[0],
                p[1],
                p[2],
                el
            );
            serial += 1;
        }
    }
    pdb.push_str("END\n");
    pdb
}

/// Converts parsed PDB atoms into a watertight, Taubin-smoothed [`TriMesh`] scaled to `target_mm`.
pub fn pdb_atoms_to_printable_mesh(
    name: &str,
    atoms: &[PdbAtom],
    style: PdbPrintStyle,
    grid_resolution: usize,
    target_mm: f32,
) -> Result<TriMesh, String> {
    if atoms.is_empty() {
        return Err(String::from("Cannot mesh empty atom list."));
    }

    let res = grid_resolution.clamp(20, 96);
    let mut mn = [f32::INFINITY; 3];
    let mut mx = [f32::NEG_INFINITY; 3];
    for a in atoms {
        for k in 0..3 {
            mn[k] = mn[k].min(a.pos_angstrom[k]);
            mx[k] = mx[k].max(a.pos_angstrom[k]);
        }
    }

    let pad = 4.5_f32;
    for k in 0..3 {
        mn[k] -= pad;
        mx[k] += pad;
    }
    let max_span = (mx[0] - mn[0]).max(mx[1] - mn[1]).max(mx[2] - mn[2]).max(4.0);
    let vox = max_span / ((res - 3) as f32);
    let nx = (((mx[0] - mn[0]) / vox).ceil() as usize + 3).max(8);
    let ny = (((mx[1] - mn[1]) / vox).ceil() as usize + 3).max(8);
    let nz = (((mx[2] - mn[2]) / vox).ceil() as usize + 3).max(8);

    let mut grid = ScalarGrid3D::new(nx, ny, nz, mn, vox, 0.0);

    match style {
        PdbPrintStyle::GaussianSurface => {
            for a in atoms {
                if a.element == "H" {
                    continue;
                }
                splat_sphere_gaussian(&mut grid, a.pos_angstrom, a.vdw_radius_angstrom * 1.25);
            }
        }
        PdbPrintStyle::BackboneTube => {
            let trace: Vec<[f32; 3]> = atoms
                .iter()
                .filter(|a| a.atom_name == "CA" || a.atom_name == "P")
                .map(|a| a.pos_angstrom)
                .collect();
            let pts = if trace.len() >= 2 {
                trace
            } else {
                atoms.iter().map(|a| a.pos_angstrom).collect()
            };
            for win in pts.windows(2) {
                let d = dist3(win[0], win[1]);
                if d <= 5.5 {
                    let steps = ((d / (vox * 0.5)).ceil() as usize).max(2);
                    for s in 0..=steps {
                        let t = (s as f32) / (steps as f32);
                        let p = [
                            win[0][0] + t * (win[1][0] - win[0][0]),
                            win[0][1] + t * (win[1][1] - win[0][1]),
                            win[0][2] + t * (win[1][2] - win[0][2]),
                        ];
                        splat_sphere_gaussian(&mut grid, p, 1.85);
                    }
                }
            }
        }
    }

    zero_grid_boundary(&mut grid);
    let mut mesh = extract_isosurface_surface_nets(name, &grid, 0.45);
    mesh.taubin_smooth(3);
    mesh.normalize_to_print_bed(target_mm);
    Ok(mesh)
}

fn splat_sphere_gaussian(grid: &mut ScalarGrid3D, center: [f32; 3], radius: f32) {
    let r_cut = radius * 2.1;
    let ix0 =
        (((center[0] - r_cut - grid.origin[0]) / grid.voxel_size).floor() as isize).max(1) as usize;
    let iy0 =
        (((center[1] - r_cut - grid.origin[1]) / grid.voxel_size).floor() as isize).max(1) as usize;
    let iz0 =
        (((center[2] - r_cut - grid.origin[2]) / grid.voxel_size).floor() as isize).max(1) as usize;

    let ix1 = (((center[0] + r_cut - grid.origin[0]) / grid.voxel_size).ceil() as usize)
        .min(grid.nx.saturating_sub(2));
    let iy1 = (((center[1] + r_cut - grid.origin[1]) / grid.voxel_size).ceil() as usize)
        .min(grid.ny.saturating_sub(2));
    let iz1 = (((center[2] + r_cut - grid.origin[2]) / grid.voxel_size).ceil() as usize)
        .min(grid.nz.saturating_sub(2));

    let inv_r2 = 1.0 / (radius * radius).max(1e-4);
    for z in iz0..=iz1 {
        for y in iy0..=iy1 {
            for x in ix0..=ix1 {
                let p = grid.world_pos(x as f32, y as f32, z as f32);
                let dx = p[0] - center[0];
                let dy = p[1] - center[1];
                let dz = p[2] - center[2];
                let d2 = dx * dx + dy * dy + dz * dz;
                let contrib = (-0.80 * d2 * inv_r2).exp();
                if contrib > 0.01 {
                    let cur = grid.get(x, y, z);
                    grid.set(x, y, z, cur + contrib);
                }
            }
        }
    }
}

/// Evaluates Johan Gielis's (2003) *Am. J. Bot.* 90:333–338 biological superformula radius `r(phi)`.
#[inline]
pub fn gielis_superformula(phi: f32, m: f32, a: f32, b: f32, n1: f32, n2: f32, n3: f32) -> f32 {
    let t1 = ((m * phi * 0.25).cos() / a).abs().powf(n2);
    let t2 = ((m * phi * 0.25).sin() / b).abs().powf(n3);
    let sum = (t1 + t2).max(1e-9);
    sum.powf(-1.0 / n1.max(0.1))
}

/// Generates a watertight 3D-printable silica **Diatom Frustule** valve (centric or multipolar
/// Gielis outline with raised marginal rim and hexagonal areolae pores).
pub fn generate_diatom_frustule_mesh(
    symmetry_lobes_m: f32,
    pore_spacing: f32,
    grid_res: usize,
    target_mm: f32,
) -> TriMesh {
    let n = grid_res.clamp(32, 88);
    let nz = (n / 2).max(22);
    let vox = 2.4 / (n as f32);
    let mut grid = ScalarGrid3D::new(n, n, nz, [-1.2, -1.2, -0.50], vox, -1.0);

    let sqrt3 = 3.0_f32.sqrt();
    let sp = pore_spacing.clamp(0.14, 0.45);
    let min_half = vox * 1.45;

    for z in 1..(nz - 1) {
        for y in 1..(n - 1) {
            for x in 1..(n - 1) {
                let p = grid.world_pos(x as f32, y as f32, z as f32);
                let rho = (p[0] * p[0] + p[1] * p[1]).sqrt();
                let phi = p[1].atan2(p[0]);
                let r_edge =
                    gielis_superformula(phi, symmetry_lobes_m, 1.0, 1.0, 2.5, 2.2, 2.2) * 0.92;

                if rho > r_edge + 0.08 {
                    continue;
                }

                let dome_z = 0.14 * (1.0 - (rho / r_edge.max(0.1)).powi(2)).max(0.0);
                let is_rim = rho > r_edge * 0.84 && rho <= r_edge;
                let half_thick = if is_rim {
                    0.12_f32.max(min_half * 1.3)
                } else {
                    0.065_f32.max(min_half)
                };
                let dz = (p[2] - dome_z).abs();

                let mut sdf = (r_edge - rho).min(half_thick - dz);

                if rho > 0.18 && rho < r_edge * 0.80 {
                    let row = (p[1] / (sp * 0.5 * sqrt3)).round();
                    let x_shift = if (row as i32).rem_euclid(2) != 0 {
                        0.5 * sp
                    } else {
                        0.0
                    };
                    let cx = ((p[0] - x_shift) / sp).round() * sp + x_shift;
                    let cy = row * sp * 0.5 * sqrt3;
                    let d_pore = ((p[0] - cx).powi(2) + (p[1] - cy).powi(2)).sqrt();
                    let pore_rad = sp * 0.30;
                    sdf = sdf.min(d_pore - pore_rad);
                }

                grid.set(x, y, z, sdf);
            }
        }
    }

    zero_negative_boundary(&mut grid);
    let mut mesh = extract_isosurface_surface_nets("Diatom_Frustule", &grid, 0.0);
    mesh.taubin_smooth(2);
    mesh.normalize_to_print_bed(target_mm);
    mesh
}

/// Generates a watertight Haeckel-inspired **Radiolarian** perforated cortical test with 12
/// icosahedral radial spines suitable for SLA or FDM printing.
pub fn generate_radiolarian_mesh(grid_res: usize, target_mm: f32) -> TriMesh {
    let n = grid_res.clamp(34, 84);
    let vox = 2.8 / (n as f32);
    let mut grid = ScalarGrid3D::new(n, n, n, [-1.4, -1.4, -1.4], vox, -1.0);

    let gr = 0.5 * (1.0 + 5.0_f32.sqrt());
    let raw_dirs: [[f32; 3]; 6] = [
        [0.0, 1.0, gr],
        [0.0, -1.0, gr],
        [1.0, gr, 0.0],
        [-1.0, gr, 0.0],
        [gr, 0.0, 1.0],
        [-gr, 0.0, 1.0],
    ];
    let norm_dirs: Vec<[f32; 3]> = raw_dirs
        .iter()
        .map(|d| {
            let l = (d[0] * d[0] + d[1] * d[1] + d[2] * d[2]).sqrt();
            [d[0] / l, d[1] / l, d[2] / l]
        })
        .collect();

    for z in 1..(n - 1) {
        for y in 1..(n - 1) {
            for x in 1..(n - 1) {
                let p = grid.world_pos(x as f32, y as f32, z as f32);
                let r = (p[0] * p[0] + p[1] * p[1] + p[2] * p[2]).sqrt();
                if r > 1.25 {
                    continue;
                }

                let shell_sdf = 0.085 - (r - 0.70).abs();
                let gyroid = (p[0] * 9.5).sin() * (p[1] * 9.5).cos()
                    + (p[1] * 9.5).sin() * (p[2] * 9.5).cos()
                    + (p[2] * 9.5).sin() * (p[0] * 9.5).cos();
                let perforated_shell = shell_sdf.min(gyroid.abs() - 0.34);

                let mut spine_sdf = -1.0_f32;
                if (0.55..=1.18).contains(&r) {
                    for d in &norm_dirs {
                        let proj = (p[0] * d[0] + p[1] * d[1] + p[2] * d[2]).abs();
                        let perp2 = (r * r - proj * proj).max(0.0);
                        let taper = 0.095 * (1.28 - 0.65 * r).max(0.35);
                        spine_sdf = spine_sdf.max(taper - perp2.sqrt());
                    }
                }

                grid.set(x, y, z, perforated_shell.max(spine_sdf));
            }
        }
    }

    zero_negative_boundary(&mut grid);
    let mut mesh = extract_isosurface_surface_nets("Radiolarian_Endoskeleton", &grid, 0.0);
    mesh.taubin_smooth(2);
    mesh.normalize_to_print_bed(target_mm);
    mesh
}

/// Converts a CRV Sim 3D voxel scalar CSV/TSV (`x,y,z,value` integer lattice rows) into a
/// watertight, Taubin-smoothed manifold [`TriMesh`].
pub fn crv_voxel_csv_to_mesh(
    name: &str,
    voxel_csv: &str,
    threshold: f32,
    target_mm: f32,
) -> Result<TriMesh, String> {
    let mut samples: Vec<(usize, usize, usize, f32)> = Vec::new();
    let (mut mx, mut my, mut mz) = (0usize, 0usize, 0usize);

    for line in voxel_csv.lines() {
        let t = line.trim();
        if t.is_empty()
            || t.starts_with('#')
            || t.chars().next().map(|c| c.is_ascii_alphabetic()).unwrap_or(false)
        {
            continue;
        }
        let cols: Vec<&str> = t.split([',', '\t', ' ']).filter(|s| !s.is_empty()).collect();
        if cols.len() < 3 {
            continue;
        }
        let x: usize = cols[0].parse().map_err(|_| format!("Bad voxel X in '{t}'"))?;
        let y: usize = cols[1].parse().map_err(|_| format!("Bad voxel Y in '{t}'"))?;
        let z: usize = cols[2].parse().map_err(|_| format!("Bad voxel Z in '{t}'"))?;
        let val: f32 = cols
            .get(3)
            .and_then(|s| s.parse::<f32>().ok())
            .unwrap_or(1.0);
        mx = mx.max(x);
        my = my.max(y);
        mz = mz.max(z);
        samples.push((x, y, z, val));
    }

    if samples.is_empty() {
        return Err(String::from("No voxel coordinates parsed from CSV/TSV input."));
    }

    let mut grid = ScalarGrid3D::new(
        mx + 5,
        my + 5,
        mz + 5,
        [0.0, 0.0, 0.0],
        1.0,
        threshold - 1.0,
    );
    for (x, y, z, v) in samples {
        grid.set(x + 2, y + 2, z + 2, v);
    }

    let mut mesh = extract_isosurface_surface_nets(name, &grid, threshold);
    mesh.taubin_smooth(2);
    mesh.normalize_to_print_bed(target_mm);
    Ok(mesh)
}

fn zero_grid_boundary(grid: &mut ScalarGrid3D) {
    set_boundary(grid, 0.0);
}

fn zero_negative_boundary(grid: &mut ScalarGrid3D) {
    set_boundary(grid, -1.0);
}

fn set_boundary(grid: &mut ScalarGrid3D, outside_val: f32) {
    let (nx, ny, nz) = (grid.nx, grid.ny, grid.nz);
    for z in 0..nz {
        for y in 0..ny {
            for x in 0..nx {
                if x == 0 || y == 0 || z == 0 || x + 1 == nx || y + 1 == ny || z + 1 == nz {
                    grid.set(x, y, z, outside_val);
                }
            }
        }
    }
}

fn dist3(a: [f32; 3], b: [f32; 3]) -> f32 {
    let dx = a[0] - b[0];
    let dy = a[1] - b[1];
    let dz = a[2] - b[2];
    (dx * dx + dy * dy + dz * dz).sqrt()
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn generates_watertight_pdb_diatom_radiolarian_and_voxel_meshes() {
        let atoms = parse_pdb_atoms(&demo_alpha_helix_pdb(), false).unwrap();
        assert!(atoms.len() >= 70);

        let pdb_mesh =
            pdb_atoms_to_printable_mesh("Helix", &atoms, PdbPrintStyle::GaussianSurface, 28, 80.0)
                .unwrap();
        let rep = pdb_mesh.inspect_printability();
        assert!(rep.is_watertight_manifold);
        assert!(rep.signed_volume_mm3 > 100.0);

        let diatom = generate_diatom_frustule_mesh(6.0, 0.28, 36, 90.0);
        assert!(diatom.inspect_printability().is_watertight_manifold);

        let rad = generate_radiolarian_mesh(34, 85.0);
        assert!(rad.inspect_printability().is_watertight_manifold);
    }
}
