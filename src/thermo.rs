// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)

//! Nearest-neighbor DNA duplex thermodynamics, mismatch penalties, and empirical
//! monovalent/divalent cation (`Mg2+` + `dNTP`) melting temperature corrections.
//!
//! ### Primary Literature References
//! - **Unified Watson–Crick NN stacks:** SantaLucia, J. Jr. (1998). *PNAS* 95(4), 1460–1465.
//! - **Internal mismatch thermodynamics:** Allawi, H. T. & SantaLucia, J. Jr. (1997–1998).
//!   *Biochemistry* 36, 10581–10594; *Biochemistry* 37, 2170–2179 & 9435–9444; *NAR* 26, 2694–2701.
//! - **Monovalent salt correction:** Owczarzy, R. et al. (2004). *Biochemistry* 43(12), 3537–3554.
//! - **Divalent `Mg2+` + `dNTP` chelation & ion competition:** Owczarzy, R. et al. (2008).
//!   *Biochemistry* 47(19), 5336–5353.
//! - **High-salt polymerase adjustment (Phusion):** Schildkraut, C. & Lifson, S. (1965).
//!   *Biopolymers* 3(2), 195–208 (as parameterized in the NEB Tm Calculator).

use crate::{
    complement_base, expand_iupac, gc_fraction, is_self_complementary, is_unambiguous_dna,
    normalize_dna, OligoError,
};

/// Universal gas constant in `cal / (mol * K)` matching SantaLucia (1998).
pub const GAS_CONST_CAL_MOL_K: f64 = 1.987;

/// Absolute zero offset in Celsius.
pub const KELVIN_OFFSET: f64 = 273.15;

/// Association constant `Ka` (`M^-1`) for `Mg2+:dNTP` 1:1 chelation in PCR buffers
/// (Owczarzy et al., 2008).
pub const MG_DNTP_KA_PER_M: f64 = 3.0e4;

/// Salt correction formulation applied to the 1 M `Na+` nearest-neighbor baseline `Tm`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SaltCorrectionMethod {
    /// Owczarzy et al. (2008) three-regime `Mg2+` vs. monovalent competition model
    /// with quadratic `dNTP` chelation (`R < 0.22`, `0.22 <= R < 6.0`, `R >= 6.0`).
    Owczarzy2008,
    /// Owczarzy et al. (2004) GC-dependent `1/Tm` monovalent cation model.
    Owczarzy2004,
    /// SantaLucia (1998) entropy adjustment `0.368 * (N - 1) * ln([Na+_eq])`.
    SantaLucia1998,
    /// Schildkraut & Lifson (1965) linear `16.6 * log10([Na+_eq])` used for Phusion buffers.
    Schildkraut1965,
}

/// Built-in PCR polymerase reaction presets matching standard vendor buffer formulations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum PolymerasePreset {
    /// Q5 High-Fidelity DNA Polymerase (500 nM primer, 2.0 mM `Mg2+`, 0.8 mM total dNTPs).
    Q5,
    /// Phusion High-Fidelity DNA Polymerase (500 nM primer, Schildkraut salt model, `Ta = Tm_min + 3 C`).
    Phusion,
    /// Standard Taq DNA Polymerase (200 nM primer, 1.5 mM `Mg2+`, 0.8 mM total dNTPs, `Ta = Tm_min - 5 C`).
    Taq,
    /// OneTaq DNA Polymerase (200 nM primer, 1.8 mM `Mg2+`, 0.8 mM total dNTPs, `Ta = Tm_min - 3 C`).
    OneTaq,
    /// Equimolar UV optical melting benchmark (25 nM + 25 nM strands, 50 mM `Na+`, 10 mM Tris, 1.5 mM `Mg2+`).
    UvMeltBenchmark,
}

/// Complete ionic, co-solvent, and strand-concentration specification for a hybridization reaction.
#[derive(Debug, Clone, Copy, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ReactionConditions {
    /// Monovalent `Na+` concentration in millimolar (`mM`).
    pub na_mm: f64,
    /// Monovalent `K+` concentration in millimolar (`mM`).
    pub k_mm: f64,
    /// `Tris+` buffer concentration in millimolar (`mM`); contributes `0.5 * [Tris]` to `[Mon+]`.
    pub tris_mm: f64,
    /// Total divalent `Mg2+` concentration in millimolar (`mM`) prior to `dNTP` chelation.
    pub mg_mm: f64,
    /// Total `dNTP` concentration in millimolar (`mM`, e.g., `0.8` for `0.2 mM` each of dATP/dCTP/dGTP/dTTP).
    pub dntp_mm: f64,
    /// Co-solvent DMSO percentage (`% v/v`, lowers `Tm` by ~0.60 C per 1% DMSO).
    pub dmso_pct: f64,
    /// Co-solvent formamide percentage (`% v/v`, lowers `Tm` by ~0.65 C per 1% formamide).
    pub formamide_pct: f64,
    /// Concentration of the excess strand (e.g., PCR primer) in nanomolar (`nM`).
    pub strand_a_nm: f64,
    /// Concentration of the limiting strand (e.g., `0.0` for picomolar PCR template, or equal to `strand_a_nm` for UV melts) in `nM`.
    pub strand_b_nm: f64,
    /// Salt correction formulation.
    pub salt_method: SaltCorrectionMethod,
}

impl Default for ReactionConditions {
    fn default() -> Self {
        Self::for_preset(PolymerasePreset::Q5)
    }
}

impl ReactionConditions {
    /// Returns calibrated reaction conditions for a standard polymerase or UV melt preset.
    pub fn for_preset(preset: PolymerasePreset) -> Self {
        match preset {
            PolymerasePreset::Q5 => Self {
                na_mm: 0.0,
                k_mm: 50.0,
                tris_mm: 25.0,
                mg_mm: 2.0,
                dntp_mm: 0.8,
                dmso_pct: 0.0,
                formamide_pct: 0.0,
                strand_a_nm: 500.0,
                strand_b_nm: 0.0,
                salt_method: SaltCorrectionMethod::Owczarzy2008,
            },
            PolymerasePreset::Phusion => Self {
                na_mm: 0.0,
                k_mm: 50.0,
                tris_mm: 20.0,
                mg_mm: 1.5,
                dntp_mm: 0.8,
                dmso_pct: 0.0,
                formamide_pct: 0.0,
                strand_a_nm: 500.0,
                strand_b_nm: 0.0,
                salt_method: SaltCorrectionMethod::Schildkraut1965,
            },
            PolymerasePreset::Taq => Self {
                na_mm: 0.0,
                k_mm: 50.0,
                tris_mm: 10.0,
                mg_mm: 1.5,
                dntp_mm: 0.8,
                dmso_pct: 0.0,
                formamide_pct: 0.0,
                strand_a_nm: 200.0,
                strand_b_nm: 0.0,
                salt_method: SaltCorrectionMethod::Owczarzy2008,
            },
            PolymerasePreset::OneTaq => Self {
                na_mm: 0.0,
                k_mm: 60.0,
                tris_mm: 20.0,
                mg_mm: 1.8,
                dntp_mm: 0.8,
                dmso_pct: 0.0,
                formamide_pct: 0.0,
                strand_a_nm: 200.0,
                strand_b_nm: 0.0,
                salt_method: SaltCorrectionMethod::Owczarzy2008,
            },
            PolymerasePreset::UvMeltBenchmark => Self {
                na_mm: 50.0,
                k_mm: 0.0,
                tris_mm: 10.0,
                mg_mm: 1.5,
                dntp_mm: 0.0,
                dmso_pct: 0.0,
                formamide_pct: 0.0,
                strand_a_nm: 25.0,
                strand_b_nm: 25.0,
                salt_method: SaltCorrectionMethod::Owczarzy2008,
            },
        }
    }

    /// Total monovalent cation concentration `[Mon+] = [Na+] + [K+] + 0.5 * [Tris+]` in molar (`M`).
    pub fn monovalent_molar(&self) -> f64 {
        (self.na_mm + self.k_mm + 0.5 * self.tris_mm).max(0.0) * 1.0e-3
    }

    /// Computes free unchelated `[Mg2+]` in molar (`M`) after 1:1 `Mg:dNTP` binding (`Ka = 3.0e4 M^-1`,
    /// Owczarzy et al., 2008).
    pub fn free_mg_molar(&self) -> f64 {
        let mg = self.mg_mm.max(0.0) * 1.0e-3;
        if mg <= 0.0 {
            return 0.0;
        }
        let dntp = self.dntp_mm.max(0.0) * 1.0e-3;
        if dntp <= 0.0 {
            return mg;
        }
        let term = MG_DNTP_KA_PER_M * dntp - MG_DNTP_KA_PER_M * mg + 1.0;
        let disc = (term * term + 4.0 * MG_DNTP_KA_PER_M * mg).sqrt();
        ((-term + disc) / (2.0 * MG_DNTP_KA_PER_M)).max(0.0)
    }

    /// Sodium-equivalent monovalent concentration in molar (`M`) according to von Ahsen et al. (2001),
    /// used when `SantaLucia1998` or `Schildkraut1965` is selected in the presence of `Mg2+`.
    pub fn sodium_equivalent_molar(&self) -> f64 {
        let mut mon_mm = (self.na_mm + self.k_mm + 0.5 * self.tris_mm).max(0.0);
        if self.mg_mm > self.dntp_mm && self.mg_mm > 0.0 {
            mon_mm += 120.0 * (self.mg_mm - self.dntp_mm.max(0.0)).sqrt();
        }
        mon_mm * 1.0e-3
    }

    /// Validates concentrations and returns the effective mass-action strand factor `k` (`M`)
    /// per SantaLucia (1998) Eq. 3: `C_A - C_B / 2` for non-self-complementary duplexes,
    /// or `C_A` for self-complementary palindromes.
    pub fn effective_strand_molar(&self, self_comp: bool) -> Result<f64, OligoError> {
        if self.strand_a_nm <= 0.0 && self.strand_b_nm <= 0.0 {
            return Err(OligoError::InvalidConcentration(
                "At least one strand concentration must be greater than 0 nM.",
            ));
        }
        if self.strand_a_nm < 0.0 || self.strand_b_nm < 0.0 {
            return Err(OligoError::InvalidConcentration(
                "Strand concentrations cannot be negative.",
            ));
        }
        let (ca, cb) = if self.strand_a_nm >= self.strand_b_nm {
            (self.strand_a_nm, self.strand_b_nm)
        } else {
            (self.strand_b_nm, self.strand_a_nm)
        };
        let eff_nm = if self_comp { ca } else { ca - 0.5 * cb };
        if eff_nm <= 0.0 {
            return Err(OligoError::InvalidConcentration(
                "Effective strand concentration must be positive.",
            ));
        }
        Ok(eff_nm * 1.0e-9)
    }
}

/// Thermodynamic output for a DNA duplex or primer–template hybrid.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct ThermodynamicResult {
    /// Standard enthalpy change `Delta H^0` in `kcal/mol` (1 M `Na+` reference state).
    pub delta_h_kcal: f64,
    /// Standard entropy change `Delta S^0` in `cal/(mol * K)` (1 M `Na+` reference state).
    pub delta_s_cal_k: f64,
    /// Gibbs free energy `Delta G^0_37` at 37.0 C (1 M `Na+`) in `kcal/mol`.
    pub delta_g_37_kcal: f64,
    /// Salt-adjusted Gibbs free energy `Delta G_T` at `eval_temp_c` in `kcal/mol`.
    pub delta_g_at_temp_kcal: f64,
    /// Uncorrected melting temperature at 1 M `Na+` in Celsius.
    pub tm_1m_na_c: f64,
    /// Final salt- and co-solvent-corrected melting temperature in Celsius.
    pub tm_salt_corrected_c: f64,
    /// Free unchelated `[Mg2+]` in `mM`.
    pub free_mg_mm: f64,
    /// Fractional GC percentage (`0.0..=100.0`).
    pub gc_pct: f64,
    /// True if the oligonucleotide is an exact Watson–Crick palindrome.
    pub is_self_complementary: bool,
}

/// Paired primer `Tm` and enzyme-specific PCR annealing temperature (`Ta`) recommendation.
#[derive(Debug, Clone, PartialEq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct AnnealingRecommendation {
    pub preset: PolymerasePreset,
    pub fwd_tm_c: f64,
    pub rev_tm_c: f64,
    pub delta_tm_c: f64,
    pub recommended_ta_c: f64,
}

/// SantaLucia (1998) unified Watson–Crick dinucleotide step `(dH_kcal_mol, dS_cal_mol_k)`.
/// `b1, b2` are consecutive `5' -> 3'` bases on the top strand paired with exact complements.
#[inline]
pub fn wc_nn_stack(b1: u8, b2: u8) -> Option<(f64, f64)> {
    match (b1, b2) {
        (b'A', b'A') | (b'T', b'T') => Some((-7.9, -22.2)),
        (b'A', b'T') => Some((-7.2, -20.4)),
        (b'T', b'A') => Some((-7.2, -21.3)),
        (b'C', b'A') | (b'T', b'G') => Some((-8.5, -22.7)),
        (b'G', b'T') | (b'A', b'C') => Some((-8.4, -22.4)),
        (b'C', b'T') | (b'A', b'G') => Some((-7.8, -21.0)),
        (b'G', b'A') | (b'T', b'C') => Some((-8.2, -22.2)),
        (b'C', b'G') => Some((-10.6, -27.2)),
        (b'G', b'C') => Some((-9.8, -24.4)),
        (b'G', b'G') | (b'C', b'C') => Some((-8.0, -19.9)),
        _ => None,
    }
}

/// SantaLucia (1998) terminal base-pair initiation `(dH_kcal_mol, dS_cal_mol_k)`.
#[inline]
pub fn wc_terminal_initiation(base: u8) -> (f64, f64) {
    match base {
        b'G' | b'C' => (0.1, -2.8),
        _ => (2.3, 4.1),
    }
}

/// Allawi & SantaLucia (1997–1998) & Peyret et al. (1999) internal single-mismatch
/// nearest-neighbor step `(dH_kcal_mol, dS_cal_mol_k)` for `5'-x1 x2-3' / 3'-y1 y2-5'`.
pub fn internal_mismatch_nn(x1: u8, x2: u8, y1: u8, y2: u8) -> Option<(f64, f64)> {
    if let Some((dh, ds)) = lookup_imm_directed(x1, x2, y1, y2) {
        return Some((dh, ds));
    }
    // Dinucleotide symmetry: 5'-x1 x2-3' / 3'-y1 y2-5' == 5'-y2 y1-3' / 3'-x2 x1-5'
    lookup_imm_directed(y2, y1, x2, x1)
}

fn lookup_imm_directed(x1: u8, x2: u8, y1: u8, y2: u8) -> Option<(f64, f64)> {
    match ([x1, x2], [y1, y2]) {
        // G-T wobble family (Allawi & SantaLucia 1997, Biochemistry 36:10581)
        ([b'A', b'G'], [b'T', b'T']) => Some((1.0, 0.9)),
        ([b'A', b'T'], [b'T', b'G']) => Some((-2.5, -8.3)),
        ([b'C', b'G'], [b'G', b'T']) => Some((-4.1, -11.7)),
        ([b'C', b'T'], [b'G', b'G']) => Some((-2.8, -8.0)),
        ([b'G', b'G'], [b'C', b'T']) => Some((3.3, 10.4)),
        ([b'G', b'G'], [b'T', b'T']) => Some((5.8, 16.3)),
        ([b'G', b'T'], [b'C', b'G']) => Some((-4.4, -12.3)),
        ([b'G', b'T'], [b'T', b'G']) => Some((4.1, 9.5)),
        ([b'T', b'G'], [b'A', b'T']) => Some((-0.1, -1.7)),
        ([b'T', b'G'], [b'G', b'T']) => Some((-1.4, -6.2)),
        ([b'T', b'T'], [b'A', b'G']) => Some((-1.3, -5.3)),
        // G-A mismatch family (Allawi & SantaLucia 1998, Biochemistry 37:2170)
        ([b'A', b'A'], [b'T', b'G']) => Some((-0.6, -2.3)),
        ([b'A', b'G'], [b'T', b'A']) => Some((-0.7, -2.3)),
        ([b'C', b'A'], [b'G', b'G']) => Some((-0.7, -2.3)),
        ([b'C', b'G'], [b'G', b'A']) => Some((-4.0, -13.2)),
        ([b'G', b'A'], [b'C', b'G']) => Some((-0.6, -1.0)),
        ([b'G', b'G'], [b'C', b'A']) => Some((0.5, 3.2)),
        ([b'T', b'A'], [b'A', b'G']) => Some((0.7, 0.7)),
        ([b'T', b'G'], [b'A', b'A']) => Some((3.0, 7.4)),
        // C-T mismatch family (Allawi & SantaLucia 1998, NAR 26:2694)
        ([b'A', b'C'], [b'T', b'T']) => Some((0.7, 0.2)),
        ([b'A', b'T'], [b'T', b'C']) => Some((-1.2, -6.2)),
        ([b'C', b'C'], [b'G', b'T']) => Some((-0.8, -4.5)),
        ([b'C', b'T'], [b'G', b'C']) => Some((-1.5, -6.1)),
        ([b'G', b'C'], [b'C', b'T']) => Some((2.3, 5.4)),
        ([b'G', b'T'], [b'C', b'C']) => Some((5.2, 13.5)),
        ([b'T', b'C'], [b'A', b'T']) => Some((1.2, 0.7)),
        ([b'T', b'T'], [b'A', b'C']) => Some((1.0, 0.7)),
        // A-C mismatch family (Allawi & SantaLucia 1998, Biochemistry 37:9435)
        ([b'A', b'A'], [b'T', b'C']) => Some((2.3, 4.6)),
        ([b'A', b'C'], [b'T', b'A']) => Some((5.3, 14.6)),
        ([b'C', b'A'], [b'G', b'C']) => Some((1.9, 3.7)),
        ([b'C', b'C'], [b'G', b'A']) => Some((0.6, -0.6)),
        ([b'G', b'A'], [b'C', b'C']) => Some((5.2, 14.2)),
        ([b'G', b'C'], [b'C', b'A']) => Some((-0.7, -3.8)),
        ([b'T', b'A'], [b'A', b'C']) => Some((3.4, 8.0)),
        ([b'T', b'C'], [b'A', b'A']) => Some((7.6, 20.2)),
        // Homopurine / homopyrimidine AA, CC, GG, TT mismatches (Peyret et al. 1999)
        ([b'A', b'A'], [b'T', b'A']) => Some((1.2, 1.7)),
        ([b'C', b'A'], [b'G', b'A']) => Some((-0.9, -4.2)),
        ([b'G', b'A'], [b'C', b'A']) => Some((-2.9, -9.8)),
        ([b'T', b'A'], [b'A', b'A']) => Some((4.7, 12.9)),
        ([b'A', b'C'], [b'T', b'C']) => Some((0.0, -4.4)),
        ([b'C', b'C'], [b'G', b'C']) => Some((-1.5, -7.2)),
        ([b'G', b'C'], [b'C', b'C']) => Some((3.6, 8.9)),
        ([b'T', b'C'], [b'A', b'C']) => Some((6.1, 16.4)),
        ([b'A', b'G'], [b'T', b'G']) => Some((-3.1, -9.5)),
        ([b'C', b'G'], [b'G', b'G']) => Some((-4.9, -15.3)),
        ([b'G', b'G'], [b'C', b'G']) => Some((-6.0, -15.8)),
        ([b'T', b'G'], [b'A', b'G']) => Some((1.6, 3.6)),
        ([b'A', b'T'], [b'T', b'T']) => Some((-2.7, -10.8)),
        ([b'C', b'T'], [b'G', b'T']) => Some((-5.0, -15.8)),
        ([b'G', b'T'], [b'C', b'T']) => Some((-2.2, -8.4)),
        ([b'T', b'T'], [b'A', b'T']) => Some((0.2, -1.5)),
        _ => None,
    }
}

/// Computes `1/Tm` correction (`K^-1`) under Owczarzy et al. (2004) monovalent equation.
fn owczarzy_2004_inv_tm_corr(f_gc: f64, mon_molar: f64) -> Result<f64, OligoError> {
    if mon_molar <= 0.0 {
        return Err(OligoError::InvalidConcentration(
            "Monovalent cation concentration must be > 0 M when Mg2+ is absent.",
        ));
    }
    let ln_mon = mon_molar.ln();
    Ok((4.29 * f_gc - 3.95) * 1.0e-5 * ln_mon + 9.40e-6 * ln_mon * ln_mon)
}

/// Computes `1/Tm` correction (`K^-1`) under Owczarzy et al. (2008) divalent + monovalent competition.
fn owczarzy_2008_inv_tm_corr(
    f_gc: f64,
    n_bp: usize,
    mon_molar: f64,
    free_mg_molar: f64,
) -> Result<f64, OligoError> {
    if free_mg_molar <= 0.0 {
        return owczarzy_2004_inv_tm_corr(f_gc, mon_molar);
    }

    let mut a = 3.92;
    let b = -0.911;
    let c = 6.26;
    let mut d = 1.42;
    let e = -48.2;
    let f = 52.5;
    let mut g = 8.31;

    if mon_molar > 0.0 {
        let ratio = free_mg_molar.sqrt() / mon_molar;
        if ratio < 0.22 {
            // Monovalent-dominant regime (Owczarzy et al. 2008, Eq. 4)
            return owczarzy_2004_inv_tm_corr(f_gc, mon_molar);
        } else if ratio < 6.0 {
            // Competitive Mg2+ / Mon+ regime (Owczarzy et al. 2008, Eqs. 18–20)
            let ln_mon = mon_molar.ln();
            a = 3.92 * (0.843 - 0.352 * mon_molar.sqrt() * ln_mon);
            d = 1.42 * (1.279 - 4.03e-3 * ln_mon - 8.03e-3 * ln_mon * ln_mon);
            g = 8.31 * (0.486 - 0.258 * ln_mon + 5.25e-3 * ln_mon * ln_mon * ln_mon);
        }
    }

    let ln_mg = free_mg_molar.ln();
    let inv_phosphates = 1.0 / (2.0 * ((n_bp.saturating_sub(1)).max(1) as f64));
    let corr = (a
        + b * ln_mg
        + f_gc * (c + d * ln_mg)
        + inv_phosphates * (e + f * ln_mg + g * ln_mg * ln_mg))
        * 1.0e-5;
    Ok(corr)
}

/// Evaluates `Tm` and Gibbs free energy for a perfectly paired DNA oligonucleotide (`5' -> 3'`).
pub fn calc_perfect_duplex(
    raw_seq: &str,
    cond: &ReactionConditions,
    eval_temp_c: f64,
) -> Result<ThermodynamicResult, OligoError> {
    let seq = normalize_dna(raw_seq)?;
    if seq.len() < 2 {
        return Err(OligoError::SequenceTooShort {
            len: seq.len(),
            min: 2,
        });
    }
    if !is_unambiguous_dna(&seq) {
        for (idx, ch) in seq.chars().enumerate() {
            if !matches!(ch, 'A' | 'C' | 'G' | 'T') {
                return Err(OligoError::AmbiguousBaseInStrictMode {
                    byte: ch,
                    pos: idx,
                });
            }
        }
    }

    let bytes = seq.as_bytes();
    let mut dh_kcal = 0.0;
    let mut ds_cal_k = 0.0;

    // Terminal initiation terms (5' and 3' ends)
    let (dh_5, ds_5) = wc_terminal_initiation(bytes[0]);
    let (dh_3, ds_3) = wc_terminal_initiation(bytes[bytes.len() - 1]);
    dh_kcal += dh_5 + dh_3;
    ds_cal_k += ds_5 + ds_3;

    // Nearest-neighbor propagation steps
    for window in bytes.windows(2) {
        let (dh, ds) = wc_nn_stack(window[0], window[1]).ok_or(OligoError::InvalidBase {
            byte: window[1] as char,
            pos: 0,
        })?;
        dh_kcal += dh;
        ds_cal_k += ds;
    }

    let self_comp = is_self_complementary(&seq);
    if self_comp {
        ds_cal_k -= 1.4;
    }

    finalize_thermodynamics(&seq, dh_kcal, ds_cal_k, self_comp, cond, eval_temp_c)
}

/// Evaluates `Tm` and `Delta G` for an aligned duplex (`seq_5to3` paired antiparallel against
/// `target_3to5` of identical length) containing Watson–Crick pairs and isolated internal mismatches.
pub fn calc_aligned_duplex(
    raw_seq_5to3: &str,
    raw_target_3to5: &str,
    cond: &ReactionConditions,
    eval_temp_c: f64,
) -> Result<ThermodynamicResult, OligoError> {
    let seq = normalize_dna(raw_seq_5to3)?;
    let target = normalize_dna(raw_target_3to5)?;
    if seq.len() < 2 || seq.len() != target.len() {
        return Err(OligoError::SequenceTooShort {
            len: seq.len().min(target.len()),
            min: 2,
        });
    }

    let s = seq.as_bytes();
    let t = target.as_bytes();
    let n = s.len();

    let mut dh_kcal = 0.0;
    let mut ds_cal_k = 0.0;

    // Terminal initiation on ends
    let (dh_5, ds_5) = wc_terminal_initiation(s[0]);
    let (dh_3, ds_3) = wc_terminal_initiation(s[n - 1]);
    dh_kcal += dh_5 + dh_3;
    ds_cal_k += ds_5 + ds_3;

    for i in 0..(n - 1) {
        let wc_left = complement_base(s[i]) == t[i];
        let wc_right = complement_base(s[i + 1]) == t[i + 1];
        if wc_left && wc_right {
            let (dh, ds) = wc_nn_stack(s[i], s[i + 1]).ok_or(OligoError::InvalidBase {
                byte: s[i] as char,
                pos: i,
            })?;
            dh_kcal += dh;
            ds_cal_k += ds;
        } else if let Some((dh, ds)) = internal_mismatch_nn(s[i], s[i + 1], t[i], t[i + 1]) {
            dh_kcal += dh;
            ds_cal_k += ds;
        } else {
            // Tandem adjacent mismatches lack closed-form trimer SVD parameters; apply empirical
            // destabilizing loop penalty (+1.8 kcal/mol, +2.0 cal/mol*K per tandem mismatch step).
            dh_kcal += 1.8;
            ds_cal_k += 2.0;
        }
    }

    finalize_thermodynamics(&seq, dh_kcal, ds_cal_k, false, cond, eval_temp_c)
}

fn finalize_thermodynamics(
    seq: &str,
    dh_kcal: f64,
    ds_cal_k: f64,
    self_comp: bool,
    cond: &ReactionConditions,
    eval_temp_c: f64,
) -> Result<ThermodynamicResult, OligoError> {
    let k_molar = cond.effective_strand_molar(self_comp)?;
    let f_gc = gc_fraction(seq);
    let free_mg_m = cond.free_mg_molar();

    let dg_37_kcal = dh_kcal - (37.0 + KELVIN_OFFSET) * (ds_cal_k * 1.0e-3);
    let denom_1m = ds_cal_k + GAS_CONST_CAL_MOL_K * k_molar.ln();
    let tm_1m_k = (dh_kcal * 1000.0) / denom_1m;
    let tm_1m_c = tm_1m_k - KELVIN_OFFSET;

    let (tm_salt_c, effective_ds_cal_k) = match cond.salt_method {
        SaltCorrectionMethod::SantaLucia1998 => {
            let na_eq = cond.sodium_equivalent_molar();
            if na_eq <= 0.0 {
                return Err(OligoError::InvalidConcentration(
                    "Sodium-equivalent concentration must be > 0 M.",
                ));
            }
            let ds_corr = 0.368 * ((seq.len() - 1) as f64) * na_eq.ln();
            let ds_adj = ds_cal_k + ds_corr;
            let tm_c = (dh_kcal * 1000.0) / (ds_adj + GAS_CONST_CAL_MOL_K * k_molar.ln())
                - KELVIN_OFFSET;
            (tm_c, ds_adj)
        }
        SaltCorrectionMethod::Schildkraut1965 => {
            let na_eq = cond.sodium_equivalent_molar();
            if na_eq <= 0.0 {
                return Err(OligoError::InvalidConcentration(
                    "Sodium-equivalent concentration must be > 0 M.",
                ));
            }
            let tm_corr = 16.6 * na_eq.log10();
            let tm_c = tm_1m_c + tm_corr;
            let implied_ds =
                (dh_kcal * 1000.0) / (tm_c + KELVIN_OFFSET) - GAS_CONST_CAL_MOL_K * k_molar.ln();
            (tm_c, implied_ds)
        }
        SaltCorrectionMethod::Owczarzy2004 => {
            let inv_corr = owczarzy_2004_inv_tm_corr(f_gc, cond.monovalent_molar())?;
            let tm_k = 1.0 / ((1.0 / tm_1m_k) + inv_corr);
            let implied_ds = (dh_kcal * 1000.0) / tm_k - GAS_CONST_CAL_MOL_K * k_molar.ln();
            (tm_k - KELVIN_OFFSET, implied_ds)
        }
        SaltCorrectionMethod::Owczarzy2008 => {
            let inv_corr = owczarzy_2008_inv_tm_corr(
                f_gc,
                seq.len(),
                cond.monovalent_molar(),
                free_mg_m,
            )?;
            let tm_k = 1.0 / ((1.0 / tm_1m_k) + inv_corr);
            let implied_ds = (dh_kcal * 1000.0) / tm_k - GAS_CONST_CAL_MOL_K * k_molar.ln();
            (tm_k - KELVIN_OFFSET, implied_ds)
        }
    };

    // Apply DMSO (von Ahsen 2001: -0.60 C per % DMSO) and formamide (-0.65 C per % formamide)
    let tm_final_c = tm_salt_c - 0.60 * cond.dmso_pct.max(0.0) - 0.65 * cond.formamide_pct.max(0.0);
    let eval_k = eval_temp_c + KELVIN_OFFSET;
    let dg_eval_kcal = dh_kcal - eval_k * (effective_ds_cal_k * 1.0e-3);

    Ok(ThermodynamicResult {
        delta_h_kcal: dh_kcal,
        delta_s_cal_k: ds_cal_k,
        delta_g_37_kcal: dg_37_kcal,
        delta_g_at_temp_kcal: dg_eval_kcal,
        tm_1m_na_c: tm_1m_c,
        tm_salt_corrected_c: tm_final_c,
        free_mg_mm: free_mg_m * 1.0e3,
        gc_pct: f_gc * 100.0,
        is_self_complementary: self_comp,
    })
}

/// Expands IUPAC degenerate bases (up to `max_variants`, e.g., `64`) and returns
/// `(lowest_tm_variant, highest_tm_variant, expanded_count)`.
pub fn calc_degenerate_tm_bounds(
    raw_seq: &str,
    cond: &ReactionConditions,
    max_variants: usize,
) -> Result<(ThermodynamicResult, ThermodynamicResult, usize), OligoError> {
    let variants = expand_iupac(raw_seq, max_variants)?;
    let count = variants.len();
    let mut min_res: Option<ThermodynamicResult> = None;
    let mut max_res: Option<ThermodynamicResult> = None;

    for v in &variants {
        let res = calc_perfect_duplex(v, cond, 37.0)?;
        if min_res
            .as_ref()
            .map_or(true, |m| res.tm_salt_corrected_c < m.tm_salt_corrected_c)
        {
            min_res = Some(res.clone());
        }
        if max_res
            .as_ref()
            .map_or(true, |m| res.tm_salt_corrected_c > m.tm_salt_corrected_c)
        {
            max_res = Some(res);
        }
    }

    Ok((min_res.unwrap(), max_res.unwrap(), count))
}

/// Computes forward/reverse primer `Tm` values and recommends the optimal PCR annealing
/// temperature `Ta` under the chosen polymerase buffer preset.
pub fn recommend_annealing_temp(
    preset: PolymerasePreset,
    fwd_seq: &str,
    rev_seq: &str,
) -> Result<AnnealingRecommendation, OligoError> {
    let cond = ReactionConditions::for_preset(preset);
    let fwd = calc_perfect_duplex(fwd_seq, &cond, 60.0)?;
    let rev = calc_perfect_duplex(rev_seq, &cond, 60.0)?;
    let min_tm = fwd.tm_salt_corrected_c.min(rev.tm_salt_corrected_c);
    let delta_tm = (fwd.tm_salt_corrected_c - rev.tm_salt_corrected_c).abs();

    let min_len = normalize_dna(fwd_seq)?.len().min(normalize_dna(rev_seq)?.len());
    let ta_raw = match preset {
        PolymerasePreset::Q5 => {
            if min_len >= 20 {
                (min_tm + 1.0).min(72.0)
            } else {
                min_tm.min(72.0)
            }
        }
        PolymerasePreset::Phusion => {
            if min_len > 20 {
                (min_tm + 3.0).min(72.0)
            } else {
                min_tm.min(72.0)
            }
        }
        PolymerasePreset::Taq | PolymerasePreset::UvMeltBenchmark => min_tm - 5.0,
        PolymerasePreset::OneTaq => min_tm - 3.0,
    };

    Ok(AnnealingRecommendation {
        preset,
        fwd_tm_c: fwd.tm_salt_corrected_c,
        rev_tm_c: rev.tm_salt_corrected_c,
        delta_tm_c: delta_tm,
        recommended_ta_c: (ta_raw * 10.0).round() / 10.0,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn matches_biopython_owczarzy_2008_doctests() {
        // Biopython Bio.SeqUtils.MeltingTemp doctest benchmark:
        // seq = 'CGTTCCAAAGATGTGGGCATGAGCTTAC', dnac1=25, dnac2=25
        // Na=50, Tris=10, Mg=1.5, dNTPs=0.0, saltcorr=7 -> 66.81 C
        // Na=50, Tris=10, Mg=1.5, dNTPs=0.6, saltcorr=7 -> 66.04 C
        let seq = "CGTTCCAAAGATGTGGGCATGAGCTTAC";
        let mut cond = ReactionConditions::for_preset(PolymerasePreset::UvMeltBenchmark);
        let res_no_dntp = calc_perfect_duplex(seq, &cond, 37.0).unwrap();
        assert!(
            (res_no_dntp.tm_salt_corrected_c - 66.81).abs() < 0.05,
            "Expected ~66.81 C, got {:.4} C",
            res_no_dntp.tm_salt_corrected_c
        );

        cond.dntp_mm = 0.6;
        let res_with_dntp = calc_perfect_duplex(seq, &cond, 37.0).unwrap();
        assert!(
            (res_with_dntp.tm_salt_corrected_c - 66.04).abs() < 0.05,
            "Expected ~66.04 C, got {:.4} C",
            res_with_dntp.tm_salt_corrected_c
        );
    }

    #[test]
    fn matches_biopython_monovalent_and_mismatch_doctests() {
        // Biopython doctest: Na=50, Mg=0, saltcorr=6/7 -> 59.78 C; saltcorr=5 -> 60.32 C
        let seq = "CGTTCCAAAGATGTGGGCATGAGCTTAC";
        let mut cond = ReactionConditions {
            na_mm: 50.0,
            k_mm: 0.0,
            tris_mm: 0.0,
            mg_mm: 0.0,
            dntp_mm: 0.0,
            dmso_pct: 0.0,
            formamide_pct: 0.0,
            strand_a_nm: 25.0,
            strand_b_nm: 25.0,
            salt_method: SaltCorrectionMethod::Owczarzy2004,
        };
        let ow04 = calc_perfect_duplex(seq, &cond, 37.0).unwrap();
        assert!((ow04.tm_salt_corrected_c - 59.78).abs() < 0.05);

        cond.salt_method = SaltCorrectionMethod::SantaLucia1998;
        let sl98 = calc_perfect_duplex(seq, &cond, 37.0).unwrap();
        assert!((sl98.tm_salt_corrected_c - 60.32).abs() < 0.05);

        // Single internal A-C mismatch at index 6: c_seq = 'GCAAGGcTTCTACACCCGTACTCGAATG' -> 55.39 C
        let mm = calc_aligned_duplex(seq, "GCAAGGCTTCTACACCCGTACTCGAATG", &cond, 37.0).unwrap();
        assert!(
            (mm.tm_salt_corrected_c - 55.39).abs() < 0.05,
            "Expected ~55.39 C, got {:.4} C",
            mm.tm_salt_corrected_c
        );
    }

    #[test]
    fn expands_degenerate_coi_primer_bounds() {
        let cond = ReactionConditions::for_preset(PolymerasePreset::Taq);
        let (low, high, count) =
            calc_degenerate_tm_bounds("GGTCAACAAATCATAAAGATATTGGYAC", &cond, 64).unwrap();
        assert_eq!(count, 2);
        assert!(high.tm_salt_corrected_c > low.tm_salt_corrected_c);
    }
}
