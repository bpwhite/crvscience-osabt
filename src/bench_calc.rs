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

//! Wet-lab bench calculators powering `osabt bench`:
//! - Lyophilized oligonucleotide resuspension (`nmol` / `OD260` / anhydrous `MW` -> `uL` TE buffer)
//! - PCR master-mix reagent tables (with pipetting overage) & 3-step/touchdown thermocycler blocks
//! - Molarity (`g = M * V * MW`), `C1 * V1 = C2 * V2` dilutions, and fold serial dilution series
//! - Temperature-corrected Henderson–Hasselbalch biological buffers (`Tris`, `HEPES`, `MES`, `Phosphate`)
//! - Centrifuge rotor `RCF (x g) <-> RPM` conversions

use crate::thermo::PolymerasePreset;
use crate::{normalize_dna, OligoError};

/// Physical properties and stock hydration recipe for a synthesized DNA oligonucleotide.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct OligoResuspensionPlan {
    pub name: String,
    pub sequence: String,
    pub length_nt: usize,
    /// Anhydrous unphosphorylated single-stranded molecular weight (`g/mol`).
    pub molecular_weight_g_mol: f64,
    /// Estimated molar extinction coefficient at 260 nm (`L / (mol * cm)`).
    pub extinction_coeff_260: f64,
    /// Delivered synthesis yield in nanomoles (`nmol`).
    pub amount_nmol: f64,
    /// Delivered mass in micrograms (`ug`).
    pub mass_ug: f64,
    /// Target master stock concentration in micromolar (`uM`, typically `100.0`).
    pub target_stock_um: f64,
    /// Volume of TE / nuclease-free water (`uL`) to add to lyophilized tube.
    pub add_buffer_ul: f64,
    /// Recipe to prepare `working_vol_ul` (e.g., `100 uL`) of `10 uM` bench working stock.
    pub working_10um_stock_aliquot_ul: f64,
    pub working_10um_water_ul: f64,
}

/// Computes anhydrous ssDNA molecular weight (`g/mol`), `epsilon_260`, and hydration volume (`uL`).
pub fn calc_oligo_resuspension(
    name: &str,
    raw_seq: &str,
    amount_nmol: f64,
    target_stock_um: f64,
    working_final_vol_ul: f64,
) -> Result<OligoResuspensionPlan, OligoError> {
    if amount_nmol <= 0.0 || target_stock_um <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Oligo nmol yield and target stock uM must both be positive.",
        ));
    }
    let seq = normalize_dna(raw_seq)?;
    let mut mw = -61.96_f64; // Remove terminal phosphate HPO3 + H
    let mut eps = 0.0_f64;

    for b in seq.bytes() {
        let (m, e) = match b {
            b'A' => (313.21, 15_400.0),
            b'C' => (289.18, 7_400.0),
            b'G' => (329.21, 11_500.0),
            b'T' => (304.20, 8_700.0),
            _ => (308.95, 10_750.0), // IUPAC degenerate average
        };
        mw += m;
        eps += e * 0.90; // Hypochromicity-adjusted ssDNA extinction
    }
    let mw = mw.max(200.0);
    let mass_ug = amount_nmol * mw * 1.0e-3;
    // C (uM) = nmol / V (mL) = (nmol * 1000) / V (uL) => V (uL) = (nmol * 1000) / target_stock_um
    let add_buffer_ul = (amount_nmol * 1000.0) / target_stock_um;

    let work_v = working_final_vol_ul.max(10.0);
    let aliquot_ul = if target_stock_um >= 10.0 {
        (10.0 * work_v) / target_stock_um
    } else {
        work_v
    };

    Ok(OligoResuspensionPlan {
        name: name.to_string(),
        sequence: seq.clone(),
        length_nt: seq.len(),
        molecular_weight_g_mol: (mw * 100.0).round() / 100.0,
        extinction_coeff_260: eps.round(),
        amount_nmol: (amount_nmol * 100.0).round() / 100.0,
        mass_ug: (mass_ug * 100.0).round() / 100.0,
        target_stock_um,
        add_buffer_ul: (add_buffer_ul * 10.0).round() / 10.0,
        working_10um_stock_aliquot_ul: (aliquot_ul * 10.0).round() / 10.0,
        working_10um_water_ul: ((work_v - aliquot_ul).max(0.0) * 10.0).round() / 10.0,
    })
}

/// Parses a vendor COA `.csv` or `.tsv` (`Name,Sequence,nmol` columns) into hydration plans.
pub fn parse_vendor_coa_table(
    table_text: &str,
    target_stock_um: f64,
) -> Result<Vec<OligoResuspensionPlan>, OligoError> {
    let mut plans = Vec::new();
    for line in table_text.lines() {
        let l = line.trim();
        if l.is_empty() || l.starts_with('#') {
            continue;
        }
        let delim = if l.contains('\t') { '\t' } else { ',' };
        let cols: Vec<&str> = l.split(delim).map(|s| s.trim().trim_matches('"')).collect();
        if cols.len() < 3 {
            continue;
        }
        if cols[0].eq_ignore_ascii_case("name") || cols[0].eq_ignore_ascii_case("oligo") {
            continue;
        }
        if let Ok(nmol) = cols[2].parse::<f64>() {
            plans.push(calc_oligo_resuspension(
                cols[0],
                cols[1],
                nmol,
                target_stock_um,
                100.0,
            )?);
        }
    }
    Ok(plans)
}

/// Single reagent line inside a PCR master-mix pipetting table.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct MasterMixReagentLine {
    pub reagent: String,
    pub stock_label: String,
    pub final_label: String,
    pub per_rxn_ul: f64,
    pub master_mix_ul: f64,
    pub add_to_master_mix: bool,
}

/// Complete PCR master-mix recipe and thermocycler program ready for bench export.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct PcrMasterMixPlan {
    pub polymerase: PolymerasePreset,
    pub num_reactions: usize,
    pub overage_fraction: f64,
    pub rxn_volume_ul: f64,
    pub dispense_mm_per_tube_ul: f64,
    pub template_per_tube_ul: f64,
    pub lines: Vec<MasterMixReagentLine>,
    pub cycling_steps: Vec<(String, String, String)>,
}

impl PcrMasterMixPlan {
    /// Renders a printable Markdown wet-lab bench sheet with interactive checkboxes.
    pub fn to_markdown_sheet(&self, assay_title: &str) -> String {
        let mut md = format!(
            "# PCR Bench Sheet — {assay_title}\n\n\
             - **Polymerase Preset:** `{:?}` | **Reactions:** `{}` (+{:.0}% overage) | **Rxn Volume:** `{:.1} uL`\n\n\
             ## 1. Master Mix Pipetting Checklist\n\n\
             | Done | Reagent | Stock | Final | 1x Rxn (uL) | Master Mix (uL) |\n\
             | :--- | :--- | :--- | :--- | ---: | ---: |\n",
            self.polymerase,
            self.num_reactions,
            self.overage_fraction * 100.0,
            self.rxn_volume_ul
        );
        for line in &self.lines {
            let mm_col = if line.add_to_master_mix {
                format!("{:.1}", line.master_mix_ul)
            } else {
                String::from("*(add per tube)*")
            };
            md.push_str(&format!(
                "| [ ] | {} | {} | {} | {:.2} | {} |\n",
                line.reagent, line.stock_label, line.final_label, line.per_rxn_ul, mm_col
            ));
        }
        md.push_str(&format!(
            "\n- **Dispense:** `{:.2} uL` Master Mix per tube + `{:.2} uL` DNA template.\n\n\
             ## 2. Thermocycler Program\n\n\
             | Step | Temperature | Duration |\n\
             | :--- | :--- | :--- |\n",
            self.dispense_mm_per_tube_ul, self.template_per_tube_ul
        ));
        for (step, temp, dur) in &self.cycling_steps {
            md.push_str(&format!("| {step} | {temp} | {dur} |\n"));
        }
        md.push_str("\n*Generated by OSABT (`osabt bench`) • CRV Science Network (GPL-3.0)*\n");
        md
    }
}

/// Builds a complete PCR master mix and cycling program for `num_reactions` tubes.
pub fn build_pcr_master_mix(
    polymerase: PolymerasePreset,
    num_reactions: usize,
    rxn_volume_ul: f64,
    overage_fraction: f64,
    template_ul_per_rxn: f64,
    anneal_temp_c: f64,
    amplicon_bp: usize,
    with_gc_enhancer: bool,
) -> Result<PcrMasterMixPlan, OligoError> {
    if num_reactions == 0 || rxn_volume_ul < 5.0 || template_ul_per_rxn >= rxn_volume_ul * 0.5 {
        return Err(OligoError::InvalidConcentration(
            "Invalid reaction count, volume (< 5 uL), or excessive template volume.",
        ));
    }

    let mult = (num_reactions as f64) * (1.0 + overage_fraction.max(0.0));
    let (buf_x, buf_name, primer_final_um, denat_c, init_sec, step_denat_sec, sec_per_kb, ext_c) =
        match polymerase {
            PolymerasePreset::Q5 => (5.0, "5X Q5 Reaction Buffer", 0.5, 98, 30, 10, 30.0, 72),
            PolymerasePreset::Phusion => (5.0, "5X Phusion HF Buffer", 0.5, 98, 30, 10, 30.0, 72),
            PolymerasePreset::OneTaq => {
                (5.0, "5X OneTaq Standard Buffer", 0.2, 94, 30, 20, 60.0, 68)
            }
            PolymerasePreset::Taq | PolymerasePreset::UvMeltBenchmark => {
                (10.0, "10X Standard Taq Buffer", 0.2, 95, 30, 20, 60.0, 68)
            }
        };

    let buf_ul = rxn_volume_ul / buf_x;
    let dntp_ul = rxn_volume_ul * (0.2 / 10.0); // 10 mM each dNTP stock -> 200 uM each final
    let primer_ul = rxn_volume_ul * (primer_final_um / 10.0); // 10 uM working stock
    let enhancer_ul = if with_gc_enhancer && matches!(polymerase, PolymerasePreset::Q5) {
        rxn_volume_ul / 5.0
    } else if with_gc_enhancer {
        rxn_volume_ul * 0.03 // 3% DMSO
    } else {
        0.0
    };
    let pol_ul = match polymerase {
        PolymerasePreset::Q5 | PolymerasePreset::Phusion => rxn_volume_ul * (0.25 / 25.0),
        _ => rxn_volume_ul * (0.125 / 25.0),
    };

    let fixed_ul =
        buf_ul + dntp_ul + 2.0 * primer_ul + enhancer_ul + pol_ul + template_ul_per_rxn;
    if fixed_ul > rxn_volume_ul {
        return Err(OligoError::InvalidConcentration(
            "Component volumes exceed total reaction volume.",
        ));
    }
    let water_ul = rxn_volume_ul - fixed_ul;

    let mut lines = vec![
        MasterMixReagentLine {
            reagent: "Nuclease-Free Water".to_string(),
            stock_label: "—".to_string(),
            final_label: "to volume".to_string(),
            per_rxn_ul: water_ul,
            master_mix_ul: water_ul * mult,
            add_to_master_mix: true,
        },
        MasterMixReagentLine {
            reagent: buf_name.to_string(),
            stock_label: format!("{buf_x:.0}X"),
            final_label: "1X".to_string(),
            per_rxn_ul: buf_ul,
            master_mix_ul: buf_ul * mult,
            add_to_master_mix: true,
        },
        MasterMixReagentLine {
            reagent: "dNTP Mix".to_string(),
            stock_label: "10 mM each".to_string(),
            final_label: "200 uM each".to_string(),
            per_rxn_ul: dntp_ul,
            master_mix_ul: dntp_ul * mult,
            add_to_master_mix: true,
        },
        MasterMixReagentLine {
            reagent: "Forward Primer".to_string(),
            stock_label: "10 uM".to_string(),
            final_label: format!("{primer_final_um:.2} uM"),
            per_rxn_ul: primer_ul,
            master_mix_ul: primer_ul * mult,
            add_to_master_mix: true,
        },
        MasterMixReagentLine {
            reagent: "Reverse Primer".to_string(),
            stock_label: "10 uM".to_string(),
            final_label: format!("{primer_final_um:.2} uM"),
            per_rxn_ul: primer_ul,
            master_mix_ul: primer_ul * mult,
            add_to_master_mix: true,
        },
    ];

    if enhancer_ul > 0.0 {
        let label = if matches!(polymerase, PolymerasePreset::Q5) {
            ("5X Q5 High GC Enhancer", "5X", "1X")
        } else {
            ("DMSO", "100%", "3% v/v")
        };
        lines.push(MasterMixReagentLine {
            reagent: label.0.to_string(),
            stock_label: label.1.to_string(),
            final_label: label.2.to_string(),
            per_rxn_ul: enhancer_ul,
            master_mix_ul: enhancer_ul * mult,
            add_to_master_mix: true,
        });
    }

    lines.push(MasterMixReagentLine {
        reagent: format!("{polymerase:?} DNA Polymerase"),
        stock_label: "2 U/uL".to_string(),
        final_label: "0.02 U/uL".to_string(),
        per_rxn_ul: pol_ul,
        master_mix_ul: pol_ul * mult,
        add_to_master_mix: true,
    });
    lines.push(MasterMixReagentLine {
        reagent: "Template DNA".to_string(),
        stock_label: "1-10 ng/uL".to_string(),
        final_label: "variable".to_string(),
        per_rxn_ul: template_ul_per_rxn,
        master_mix_ul: template_ul_per_rxn,
        add_to_master_mix: false,
    });

    let ext_sec = (((amplicon_bp.max(100) as f64 / 1000.0) * sec_per_kb).ceil() as u32).max(15);
    let cycling_steps = vec![
        (
            "1. Initial Denaturation".to_string(),
            format!("{denat_c} C"),
            format!("{init_sec} s"),
        ),
        (
            "2. 30-35x Denature".to_string(),
            format!("{denat_c} C"),
            format!("{step_denat_sec} s"),
        ),
        (
            "3. 30-35x Anneal".to_string(),
            format!("{anneal_temp_c:.1} C"),
            "20 s".to_string(),
        ),
        (
            "4. 30-35x Extend".to_string(),
            format!("{ext_c} C"),
            format!("{ext_sec} s"),
        ),
        (
            "5. Final Extension".to_string(),
            format!("{ext_c} C"),
            "120 s".to_string(),
        ),
        (
            "6. Hold".to_string(),
            "4-10 C".to_string(),
            "infinite".to_string(),
        ),
    ];

    Ok(PcrMasterMixPlan {
        polymerase,
        num_reactions,
        overage_fraction,
        rxn_volume_ul,
        dispense_mm_per_tube_ul: rxn_volume_ul - template_ul_per_rxn,
        template_per_tube_ul: template_ul_per_rxn,
        lines,
        cycling_steps,
    })
}

/// Supported biological buffer conjugate acid/base systems with `dpKa/dT` thermal coefficients.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum BufferSystem {
    /// Tris-HCl + Tris Base (`pKa_25 = 8.06`, `dpKa/dT = -0.028 / C`).
    Tris,
    /// HEPES Free Acid + HEPES Na-Salt (`pKa_25 = 7.48`, `dpKa/dT = -0.014 / C`).
    Hepes,
    /// MES Free Acid + MES Na-Salt (`pKa_25 = 6.15`, `dpKa/dT = -0.011 / C`).
    Mes,
    /// Sodium Phosphate (`NaH2PO4*H2O` + `Na2HPO4`, `pKa2_25 = 6.86` at 0.1 M ionic strength, `dpKa/dT = -0.0028 / C`).
    SodiumPhosphate,
}

/// Two-component dry-salt Henderson–Hasselbalch buffer recipe at a specified working temperature.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BufferRecipe {
    pub system: BufferSystem,
    pub target_ph: f64,
    pub working_temp_c: f64,
    pub effective_pka: f64,
    pub total_molarity_m: f64,
    pub volume_liters: f64,
    pub acid_component_name: String,
    pub acid_grams: f64,
    pub base_component_name: String,
    pub base_grams: f64,
}

/// Computes exact dry-chemical masses (`g`) of conjugate acid and conjugate base required
/// to hit `target_ph` at `working_temp_c` (accounting for thermal `pKa` drift).
pub fn calc_buffer_recipe(
    system: BufferSystem,
    target_ph: f64,
    total_molarity_m: f64,
    volume_liters: f64,
    working_temp_c: f64,
) -> Result<BufferRecipe, OligoError> {
    if total_molarity_m <= 0.0 || volume_liters <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Buffer molarity and volume must be greater than zero.",
        ));
    }

    let (pka_25, dpka_dt, acid_name, mw_acid, base_name, mw_base) = match system {
        BufferSystem::Tris => (
            8.06,
            -0.028,
            "Tris-HCl",
            157.60,
            "Tris Base",
            121.14,
        ),
        BufferSystem::Hepes => (
            7.48,
            -0.014,
            "HEPES Free Acid",
            238.30,
            "HEPES Sodium Salt",
            260.29,
        ),
        BufferSystem::Mes => (
            6.15,
            -0.011,
            "MES Free Acid",
            195.24,
            "MES Sodium Salt",
            217.22,
        ),
        BufferSystem::SodiumPhosphate => (
            6.86,
            -0.0028,
            "NaH2PO4 * H2O (Monobasic)",
            137.99,
            "Na2HPO4 Anhydrous (Dibasic)",
            141.96,
        ),
    };

    let eff_pka = pka_25 + dpka_dt * (working_temp_c - 25.0);
    if (target_ph - eff_pka).abs() > 1.5 {
        return Err(OligoError::InvalidConcentration(
            "Target pH is outside the +/- 1.5 pH unit buffering capacity of this salt pair at the target temperature.",
        ));
    }

    // Henderson-Hasselbalch: pH = pKa + log10([Base]/[Acid]) => ratio = 10^(pH - pKa)
    let ratio = 10.0_f64.powf(target_ph - eff_pka);
    let total_moles = total_molarity_m * volume_liters;
    let moles_acid = total_moles / (1.0 + ratio);
    let moles_base = total_moles - moles_acid;

    Ok(BufferRecipe {
        system,
        target_ph,
        working_temp_c,
        effective_pka: (eff_pka * 1000.0).round() / 1000.0,
        total_molarity_m,
        volume_liters,
        acid_component_name: acid_name.to_string(),
        acid_grams: (moles_acid * mw_acid * 1000.0).round() / 1000.0,
        base_component_name: base_name.to_string(),
        base_grams: (moles_base * mw_base * 1000.0).round() / 1000.0,
    })
}

/// Computes dry solute mass in grams (`g = M * V_L * MW_g_mol`).
pub fn calc_molarity_grams(
    molarity_m: f64,
    volume_liters: f64,
    formula_weight_g_mol: f64,
) -> Result<f64, OligoError> {
    if molarity_m <= 0.0 || volume_liters <= 0.0 || formula_weight_g_mol <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Molarity, volume, and formula weight must all be positive.",
        ));
    }
    Ok(molarity_m * volume_liters * formula_weight_g_mol)
}

/// Solves `C1 * V1 = C2 * V2` returning `(stock_volume_v1, diluent_volume)` in identical volume units.
pub fn calc_c1v1_dilution(
    stock_c1: f64,
    target_c2: f64,
    final_v2: f64,
) -> Result<(f64, f64), OligoError> {
    if stock_c1 <= 0.0 || target_c2 <= 0.0 || final_v2 <= 0.0 || target_c2 > stock_c1 {
        return Err(OligoError::InvalidConcentration(
            "Stock C1 must be >= target C2, and all concentrations/volumes must be positive.",
        ));
    }
    let v1 = (target_c2 * final_v2) / stock_c1;
    Ok((v1, final_v2 - v1))
}

/// Step entry in a serial dilution cascade.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct SerialDilutionStep {
    pub tube_index: usize,
    pub concentration: f64,
    pub transfer_in_ul: f64,
    pub diluent_ul: f64,
    pub remaining_final_ul: f64,
}

/// Generates an `n_steps` serial dilution cascade with constant `dilution_factor` (e.g., `10.0` for 1:10)
/// leaving `final_usable_vol_ul` in every tube.
pub fn calc_serial_dilution(
    start_conc: f64,
    dilution_factor: f64,
    n_steps: usize,
    final_usable_vol_ul: f64,
) -> Result<Vec<SerialDilutionStep>, OligoError> {
    if start_conc <= 0.0 || dilution_factor <= 1.0 || n_steps == 0 || final_usable_vol_ul <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Dilution factor must be > 1.0, and concentration, steps, and usable volume must be > 0.",
        ));
    }
    // Each tube initially receives `diluent_ul + transfer_ul`, then passes `transfer_ul` onward (except last tube):
    // V_total = final_usable_vol_ul * (dilution_factor / (dilution_factor - 1.0))
    let transfer_ul = final_usable_vol_ul / (dilution_factor - 1.0);
    let diluent_ul = final_usable_vol_ul;
    let mut steps = Vec::with_capacity(n_steps);
    let mut conc = start_conc;

    for idx in 1..=n_steps {
        conc /= dilution_factor;
        let rem = if idx == n_steps {
            diluent_ul + transfer_ul
        } else {
            final_usable_vol_ul
        };
        steps.push(SerialDilutionStep {
            tube_index: idx,
            concentration: conc,
            transfer_in_ul: (transfer_ul * 100.0).round() / 100.0,
            diluent_ul: (diluent_ul * 100.0).round() / 100.0,
            remaining_final_ul: (rem * 100.0).round() / 100.0,
        });
    }
    Ok(steps)
}

/// Converts rotor speed `RPM` and radius `r_cm` into relative centrifugal force (`x g`).
#[inline]
pub fn rpm_to_rcf(rpm: f64, radius_cm: f64) -> Result<f64, OligoError> {
    if rpm < 0.0 || radius_cm <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Rotor radius must be > 0 cm and RPM >= 0.",
        ));
    }
    Ok((1.118e-5 * radius_cm * rpm * rpm * 10.0).round() / 10.0)
}

/// Converts desired relative centrifugal force (`x g`) and radius `r_cm` into rotor `RPM`.
#[inline]
pub fn rcf_to_rpm(rcf_g: f64, radius_cm: f64) -> Result<f64, OligoError> {
    if rcf_g < 0.0 || radius_cm <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Rotor radius must be > 0 cm and RCF >= 0.",
        ));
    }
    Ok(((rcf_g / (1.118e-5 * radius_cm)).sqrt()).round())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn resuspends_oligo_and_shifts_tris_pka_with_temperature() {
        // 28.4 nmol oligo -> 284.0 uL for 100 uM stock
        let plan =
            calc_oligo_resuspension("COX1_F", "GGTCAACAAATCATAAAGATATTGG", 28.4, 100.0, 100.0)
                .unwrap();
        assert!((plan.add_buffer_ul - 284.0).abs() < 1e-6);
        assert!((plan.working_10um_stock_aliquot_ul - 10.0).abs() < 1e-6);

        // Tris buffer at 37 C has lower effective pKa (~7.724) than at 25 C (8.06), requiring more Tris Base
        let r_25 = calc_buffer_recipe(BufferSystem::Tris, 8.0, 0.05, 1.0, 25.0).unwrap();
        let r_37 = calc_buffer_recipe(BufferSystem::Tris, 8.0, 0.05, 1.0, 37.0).unwrap();
        assert!(r_37.base_grams > r_25.base_grams);
        assert!(r_37.acid_grams < r_25.acid_grams);

        // RCF <-> RPM round-trip at r = 8.5 cm, 12,000 RPM -> ~13,684 x g
        let g = rpm_to_rcf(12_000.0, 8.5).unwrap();
        assert!((g - 13_684.3).abs() < 1.0);
        let rpm = rcf_to_rpm(g, 8.5).unwrap();
        assert!((rpm - 12_000.0).abs() <= 1.0);
    }
}
