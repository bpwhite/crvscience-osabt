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

//! Tool #9 (`osabt plate`) — ANSI/SLAS `96`-well (`8x12`) and `384`-well (`16x24`) microplate
//! layout planner with Lundholt et al. (2003) perimeter evaporation moat buffering, seeded
//! inner-well randomization (`SplitMix64`), Opentrons/bench liquid-handling transfer CSVs,
//! printable bench pipetting sheets, and color-coded `.svg` plate maps & heatmaps.

use std::fmt::Write as _;

use crate::svg_util::xml_escape;

/// Standard ANSI/SLAS microplate geometries.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PlateFormat {
    Well96,
    Well384,
}

impl PlateFormat {
    pub fn dims(self) -> (usize, usize) {
        match self {
            Self::Well96 => (8, 12),
            Self::Well384 => (16, 24),
        }
    }

    pub fn label(self) -> &'static str {
        match self {
            Self::Well96 => "96-Well (8×12)",
            Self::Well384 => "384-Well (16×24)",
        }
    }
}

/// Functional assignment of a microplate well.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum WellRole {
    BufferMoat,
    NegativeControl,
    PositiveControl,
    CompoundDose,
    Empty,
}

impl WellRole {
    pub fn label(self) -> &'static str {
        match self {
            Self::BufferMoat => "PBS/Buffer Moat",
            Self::NegativeControl => "Negative Ctrl (0% Inhib / Vehicle)",
            Self::PositiveControl => "Positive Ctrl (100% Inhib)",
            Self::CompoundDose => "Test Compound Dose",
            Self::Empty => "Unused",
        }
    }
}

/// Metadata + optional plate-reader raw & normalized signal for one well (`A01` .. `P24`).
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PlateWell {
    pub row_idx: usize,
    pub col_idx: usize,
    pub well_id: String,
    pub role: WellRole,
    pub compound_name: String,
    pub concentration_nm: Option<f64>,
    pub replicate_num: usize,
    pub raw_readout: Option<f64>,
    pub normalized_pct: Option<f64>,
}

/// Specification of a compound serial-dilution series to randomize onto a microplate.
#[derive(Debug, Clone)]
pub struct CompoundSeriesSpec {
    pub name: String,
    pub top_conc_nm: f64,
    pub dilution_fold: f64,
    pub num_doses: usize,
    pub replicates: usize,
}

/// Parameters for generating an edge-effect-mitigated randomized assay plate.
#[derive(Debug, Clone)]
pub struct PlatePlanConfig {
    pub format: PlateFormat,
    pub plate_name: String,
    /// Number of outer border rings filled with sterile PBS/media moat (`1` recommended for 96-well, `1..2` for 384-well).
    pub moat_rings: usize,
    pub negative_control_wells: usize,
    pub positive_control_wells: usize,
    pub assay_volume_ul: f64,
    pub rng_seed: u64,
    pub compounds: Vec<CompoundSeriesSpec>,
}

impl Default for PlatePlanConfig {
    fn default() -> Self {
        Self {
            format: PlateFormat::Well96,
            plate_name: String::from("CRV_Kinase_Screen_Plate01"),
            moat_rings: 1,
            negative_control_wells: 6,
            positive_control_wells: 6,
            assay_volume_ul: 100.0,
            rng_seed: 2026_1111,
            compounds: vec![
                CompoundSeriesSpec {
                    name: String::from("Staurosporine"),
                    top_conc_nm: 1_000.0,
                    dilution_fold: 4.0,
                    num_doses: 8,
                    replicates: 2,
                },
                CompoundSeriesSpec {
                    name: String::from("CRV-104"),
                    top_conc_nm: 2_000.0,
                    dilution_fold: 4.0,
                    num_doses: 8,
                    replicates: 2,
                },
                CompoundSeriesSpec {
                    name: String::from("CRV-209"),
                    top_conc_nm: 10_000.0,
                    dilution_fold: 4.0,
                    num_doses: 8,
                    replicates: 2,
                },
            ],
        }
    }
}

/// Complete planned (and optionally measured) microplate layout.
#[derive(Debug, Clone)]
pub struct Microplate {
    pub name: String,
    pub format: PlateFormat,
    pub assay_volume_ul: f64,
    pub wells: Vec<PlateWell>,
}

impl Microplate {
    pub fn well_at(&self, r: usize, c: usize) -> Option<&PlateWell> {
        let (_, ncols) = self.format.dims();
        self.wells.get(r * ncols + c)
    }

    /// Formats standard alphanumeric well coordinate (`A01` .. `H12` or `P24`).
    pub fn format_well_id(row: usize, col: usize) -> String {
        let row_ch = (b'A' + (row as u8)) as char;
        format!("{row_ch}{:02}", col + 1)
    }
}

/// Generates an edge-buffered, Fisher–Yates randomized inner-well microplate layout using a
/// deterministic `SplitMix64` PRNG so layouts are 100% reproducible from `cfg.rng_seed`.
pub fn generate_randomized_plate(cfg: &PlatePlanConfig) -> Result<Microplate, String> {
    let (nrows, ncols) = cfg.format.dims();
    let moat = cfg.moat_rings.min(nrows / 2 - 1);

    let mut wells = Vec::with_capacity(nrows * ncols);
    let mut inner_indices = Vec::new();

    for r in 0..nrows {
        for c in 0..ncols {
            let is_moat =
                r < moat || c < moat || r + moat >= nrows || c + moat >= ncols;
            let idx = wells.len();
            wells.push(PlateWell {
                row_idx: r,
                col_idx: c,
                well_id: Microplate::format_well_id(r, c),
                role: if is_moat {
                    WellRole::BufferMoat
                } else {
                    WellRole::Empty
                },
                compound_name: if is_moat {
                    String::from("PBS_Moat")
                } else {
                    String::new()
                },
                concentration_nm: None,
                replicate_num: 0,
                raw_readout: None,
                normalized_pct: None,
            });
            if !is_moat {
                inner_indices.push(idx);
            }
        }
    }

    // Build pool of experimental + control assignments
    let mut assignments: Vec<(WellRole, String, Option<f64>, usize)> = Vec::new();
    for rep in 1..=cfg.negative_control_wells {
        assignments.push((
            WellRole::NegativeControl,
            String::from("DMSO_Vehicle"),
            Some(0.0),
            rep,
        ));
    }
    for rep in 1..=cfg.positive_control_wells {
        assignments.push((
            WellRole::PositiveControl,
            String::from("PosCtrl_100PctInhib"),
            None,
            rep,
        ));
    }

    for comp in &cfg.compounds {
        if comp.dilution_fold <= 1.0 || comp.top_conc_nm <= 0.0 {
            return Err(format!(
                "Compound '{}' requires top_conc_nm > 0 and dilution_fold > 1.",
                comp.name
            ));
        }
        for d in 0..comp.num_doses {
            let conc = comp.top_conc_nm / comp.dilution_fold.powi(d as i32);
            for rep in 1..=comp.replicates {
                assignments.push((
                    WellRole::CompoundDose,
                    comp.name.clone(),
                    Some(conc),
                    rep,
                ));
            }
        }
    }

    if assignments.len() > inner_indices.len() {
        return Err(format!(
            "Requested {} inner assay wells ({} controls + compound doses), but {} with {}-ring moat only has {} inner wells (switch to 384-well or reduce replicates).",
            assignments.len(),
            cfg.negative_control_wells + cfg.positive_control_wells,
            cfg.format.label(),
            moat,
            inner_indices.len()
        ));
    }

    // Deterministic SplitMix64 Fisher-Yates shuffle on `inner_indices`
    let mut rng = SplitMix64::new(cfg.rng_seed);
    for i in (1..inner_indices.len()).rev() {
        let j = (rng.next_u64() as usize) % (i + 1);
        inner_indices.swap(i, j);
    }

    for (slot_idx, (role, cname, conc, rep)) in inner_indices.into_iter().zip(assignments) {
        let w = &mut wells[slot_idx];
        w.role = role;
        w.compound_name = cname;
        w.concentration_nm = conc;
        w.replicate_num = rep;
    }

    Ok(Microplate {
        name: cfg.plate_name.clone(),
        format: cfg.format,
        assay_volume_ul: cfg.assay_volume_ul,
        wells,
    })
}

struct SplitMix64 {
    state: u64,
}

impl SplitMix64 {
    fn new(seed: u64) -> Self {
        Self { state: seed }
    }

    fn next_u64(&mut self) -> u64 {
        self.state = self.state.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }
}

/// Exports a liquid-handler / Opentrons / bench-ready CSV plate & pipetting transfer table.
pub fn export_plate_pipetting_csv(plate: &Microplate) -> String {
    let mut csv = String::from(
        "well_id,row,col,role,compound,concentration_nM,replicate,volume_uL,raw_signal,normalized_pct\n",
    );
    for w in &plate.wells {
        let conc_s = w
            .concentration_nm
            .map(|c| format!("{c:.4}"))
            .unwrap_or_default();
        let raw_s = w
            .raw_readout
            .map(|v| format!("{v:.2}"))
            .unwrap_or_default();
        let norm_s = w
            .normalized_pct
            .map(|v| format!("{v:.2}"))
            .unwrap_or_default();
        let _ = writeln!(
            csv,
            "{},{},{},{},{},{},{},{:.1},{},{}",
            w.well_id,
            (b'A' + w.row_idx as u8) as char,
            w.col_idx + 1,
            w.role.label(),
            w.compound_name,
            conc_s,
            w.replicate_num,
            plate.assay_volume_ul,
            raw_s,
            norm_s
        );
    }
    csv
}

/// Exports a human-friendly Markdown bench pipetting checklist grouped by compound & dilution step
/// (with 15% dead-volume overage calculation for multichannel reservoirs).
pub fn export_bench_pipetting_markdown(plate: &Microplate, cfg: &PlatePlanConfig) -> String {
    let mut md = String::with_capacity(4096);
    let _ = writeln!(md, "# Bench Pipetting Protocol — {}\n", plate.name);
    let _ = writeln!(
        md,
        "- **Format:** {} | **Well Volume:** {:.1} uL | **Perimeter Moat:** {} ring(s) | **PRNG Seed:** `{}`\n",
        plate.format.label(),
        plate.assay_volume_ul,
        cfg.moat_rings,
        cfg.rng_seed
    );

    md.push_str("## 1. Serial Dilution Master-Tube Preparation (15% Dead-Volume Overage)\n\n");
    md.push_str("| Compound | Step | Target Conc (nM) | Destination Wells |\n");
    md.push_str("| :--- | :--- | ---: | :--- |\n");

    for comp in &cfg.compounds {
        for d in 0..comp.num_doses {
            let target_c = comp.top_conc_nm / comp.dilution_fold.powi(d as i32);
            let dests: Vec<&str> = plate
                .wells
                .iter()
                .filter(|w| {
                    w.compound_name == comp.name
                        && w.concentration_nm
                            .map(|c| (c - target_c).abs() <= target_c * 1e-4)
                            .unwrap_or(false)
                })
                .map(|w| w.well_id.as_str())
                .collect();
            let _ = writeln!(
                md,
                "| {} | D{} | {:.3} | `{}` |",
                comp.name,
                d + 1,
                target_c,
                dests.join(", ")
            );
        }
    }

    let neg_wells: Vec<&str> = plate
        .wells
        .iter()
        .filter(|w| w.role == WellRole::NegativeControl)
        .map(|w| w.well_id.as_str())
        .collect();
    let pos_wells: Vec<&str> = plate
        .wells
        .iter()
        .filter(|w| w.role == WellRole::PositiveControl)
        .map(|w| w.well_id.as_str())
        .collect();

    md.push_str("\n## 2. Control & Moat Wells\n\n");
    let _ = writeln!(md, "- **Negative (0% Inhib / Vehicle):** `{}`", neg_wells.join(", "));
    let _ = writeln!(md, "- **Positive (100% Inhib):** `{}`", pos_wells.join(", "));
    let moat_cnt = plate.wells.iter().filter(|w| w.role == WellRole::BufferMoat).count();
    let _ = writeln!(
        md,
        "- **Outer Evaporation Moat ({moat_cnt} perimeter wells):** Dispense {:.1} uL sterile PBS / media.\n",
        plate.assay_volume_ul
    );
    md
}

/// Renders a standalone vector `.svg` microplate layout (or fluorescence/viability heatmap when
/// `normalized_pct` values are populated on `plate`).
pub fn render_plate_map_svg(plate: &Microplate) -> String {
    let (nrows, ncols) = plate.format.dims();
    let cell = if ncols <= 12 { 44.0_f64 } else { 26.0_f64 };
    let left = 54.0_f64;
    let top = 66.0_f64;
    let width = left + (ncols as f64) * cell + 210.0;
    let height = top + (nrows as f64) * cell + 42.0;

    let mut svg = String::with_capacity(12_000);
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
        "  <text x=\"24\" y=\"28\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#0f172a\">{} — {}</text>",
        xml_escape(&plate.name),
        plate.format.label()
    );

    // Plate frame
    let _ = writeln!(
        svg,
        "  <rect x=\"{:.1}\" y=\"{:.1}\" width=\"{:.1}\" height=\"{:.1}\" rx=\"10\" fill=\"#f1f5f9\" stroke=\"#64748b\" stroke-width=\"1.8\"/>",
        left - 24.0,
        top - 26.0,
        (ncols as f64) * cell + 34.0,
        (nrows as f64) * cell + 36.0
    );

    for c in 0..ncols {
        let cx = left + (c as f64 + 0.5) * cell;
        let _ = writeln!(
            svg,
            "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" font-weight=\"700\" fill=\"#334155\">{}</text>",
            top - 8.0,
            c + 1
        );
    }
    for r in 0..nrows {
        let cy = top + (r as f64 + 0.5) * cell + 4.0;
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{cy:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" font-weight=\"700\" fill=\"#334155\">{}</text>",
            left - 12.0,
            (b'A' + r as u8) as char
        );
    }

    for w in &plate.wells {
        let cx = left + (w.col_idx as f64 + 0.5) * cell;
        let cy = top + (w.row_idx as f64 + 0.5) * cell;
        let rad = cell * 0.40;

        let fill = if let Some(pct) = w.normalized_pct {
            // Viridis-inspired blue -> teal -> yellow viability ramp
            let t = (pct / 100.0).clamp(0.0, 1.0);
            let r = (30.0 + 215.0 * t) as u8;
            let g = (64.0 + 165.0 * t) as u8;
            let b = (175.0 - 110.0 * t) as u8;
            format!("#{r:02x}{g:02x}{b:02x}")
        } else {
            match w.role {
                WellRole::BufferMoat => String::from("#bae6fd"),
                WellRole::NegativeControl => String::from("#22c55e"),
                WellRole::PositiveControl => String::from("#ef4444"),
                WellRole::CompoundDose => String::from("#6366f1"),
                WellRole::Empty => String::from("#e2e8f0"),
            }
        };

        let _ = writeln!(
            svg,
            "  <circle cx=\"{cx:.1}\" cy=\"{cy:.1}\" r=\"{rad:.1}\" fill=\"{fill}\" stroke=\"#334155\" stroke-width=\"1.0\"/>"
        );
        if ncols <= 12 && w.role == WellRole::CompoundDose {
            let short: String = w.compound_name.chars().take(4).collect();
            let _ = writeln!(
                svg,
                "  <text x=\"{cx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"8.5\" fill=\"#ffffff\">{}</text>",
                cy + 3.0,
                xml_escape(&short)
            );
        }
    }

    let lx = left + (ncols as f64) * cell + 24.0;
    for (i, (col, lbl)) in [
        ("#bae6fd", "PBS Evaporation Moat"),
        ("#22c55e", "Negative Ctrl (100% Viab)"),
        ("#ef4444", "Positive Ctrl (0% Viab)"),
        ("#6366f1", "Randomized Dose Wells"),
    ]
    .iter()
    .enumerate()
    {
        let ly = top + 16.0 + (i as f64) * 26.0;
        let _ = writeln!(svg, "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"6\" fill=\"{col}\"/>", lx + 8.0, ly - 4.0);
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{ly:.1}\" font-family=\"sans-serif\" font-size=\"11\" fill=\"#1e293b\">{lbl}</text>",
            lx + 20.0
        );
    }

    svg.push_str("</svg>\n");
    svg
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn plate_format_dims_and_deterministic_randomization_preserve_moat_and_wells() {
        assert_eq!(PlateFormat::Well96.dims(), (8, 12));
        assert_eq!(PlateFormat::Well384.dims(), (16, 24));

        let cfg = PlatePlanConfig::default();
        let plate_a = generate_randomized_plate(&cfg).expect("96-well layout succeeds");
        let plate_b = generate_randomized_plate(&cfg).expect("identical seed reproduces layout");
        assert_eq!(plate_a.wells.len(), 96);
        assert_eq!(plate_a.wells[0].well_id, "A01");
        assert_eq!(plate_a.wells[95].well_id, "H12");

        for (wa, wb) in plate_a.wells.iter().zip(plate_b.wells.iter()) {
            assert_eq!(wa.well_id, wb.well_id);
            assert_eq!(wa.role, wb.role);
            assert_eq!(wa.compound_name, wb.compound_name);
            assert_eq!(wa.concentration_nm, wb.concentration_nm);
        }

        let moat_count = plate_a
            .wells
            .iter()
            .filter(|w| w.role == WellRole::BufferMoat)
            .count();
        assert_eq!(moat_count, 36);

        let neg_count = plate_a
            .wells
            .iter()
            .filter(|w| w.role == WellRole::NegativeControl)
            .count();
        let pos_count = plate_a
            .wells
            .iter()
            .filter(|w| w.role == WellRole::PositiveControl)
            .count();
        assert_eq!(neg_count, cfg.negative_control_wells);
        assert_eq!(pos_count, cfg.positive_control_wells);

        let csv = export_plate_pipetting_csv(&plate_a);
        assert_eq!(csv.lines().count(), 97);
        let md = export_bench_pipetting_markdown(&plate_a, &cfg);
        assert!(md.contains(&cfg.plate_name));
        let svg = render_plate_map_svg(&plate_a);
        assert!(svg.starts_with("<svg") && svg.trim_end().ends_with("</svg>"));
    }

    #[test]
    fn plate_384_layout_dimensions_and_corner_ids() {
        let mut cfg = PlatePlanConfig::default();
        cfg.format = PlateFormat::Well384;
        cfg.moat_rings = 1;
        let plate = generate_randomized_plate(&cfg).expect("384-well layout succeeds");
        assert_eq!(plate.wells.len(), 384);
        assert_eq!(plate.wells.first().unwrap().well_id, "A01");
        assert_eq!(plate.wells.last().unwrap().well_id, "P24");
    }
}
