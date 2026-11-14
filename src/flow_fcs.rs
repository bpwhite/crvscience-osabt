// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #18 (`osabt flow`) — Binary `ISAC FCS 3.1` list-mode parser, writer, and spectral
//! spillover compensation engine:
//! - **Data File Standard for Flow Cytometry, version `FCS 3.1`** — Spidlen et al. (2010)
//!   *Cytometry A* 77(1):97–100 (`https://pubmed.ncbi.nlm.nih.gov/19937951/`)
//! - Supports `f32` (`F`), `f64` (`D`), and unsigned integer (`I`) list-mode binary payloads in
//!   little-endian (`1,2,3,4`) or big-endian (`4,3,2,1`) byte order, plus `$SPILLOVER` matrix inversion.

use std::collections::BTreeMap;
use std::fmt::Write as _;

/// Channel metadata parsed from `$PnN`, `$PnS`, `$PnB`, and `$PnR` keywords.
#[derive(Debug, Clone)]
pub struct FcsChannelInfo {
    pub index_1based: usize,
    pub short_name: String,
    pub stain_name: String,
    pub bits: usize,
    pub range: f64,
}

/// Square spectral spillover matrix `$SPILLOVER` (`S[emitter_row][detector_col]`).
#[derive(Debug, Clone)]
pub struct SpilloverMatrix {
    pub channel_names: Vec<String>,
    /// Row-major `n x n` spillover matrix `S` (`1.0` on diagonal).
    pub matrix: Vec<Vec<f64>>,
}

/// Complete parsed `FCS 3.0 / 3.1` list-mode flow cytometry dataset.
#[derive(Debug, Clone)]
pub struct FcsDataset {
    pub version: String,
    pub num_events: usize,
    pub channels: Vec<FcsChannelInfo>,
    /// Row-major event matrix: `events[event_idx][channel_idx]`.
    pub events: Vec<Vec<f64>>,
    pub spillover: Option<SpilloverMatrix>,
    pub keywords: BTreeMap<String, String>,
}

/// Parses a raw binary `ISAC FCS 3.0 / 3.1` byte slice into [`FcsDataset`].
pub fn parse_fcs_bytes(bytes: &[u8]) -> Result<FcsDataset, String> {
    if bytes.len() < 58 {
        return Err(String::from("File too small to be a valid FCS 3.0/3.1 file (< 58 bytes)."));
    }
    let version = std::str::from_utf8(&bytes[0..6])
        .map_err(|_| String::from("Invalid ASCII header in FCS file."))?
        .trim()
        .to_string();
    if !version.starts_with("FCS") {
        return Err(format!("Unsupported file signature '{version}' (expected FCS3.0 or FCS3.1)."));
    }

    let parse_offset = |slice: &[u8], label: &str| -> Result<usize, String> {
        let s = std::str::from_utf8(slice)
            .map_err(|_| format!("Invalid ASCII in FCS header offset {label}."))?
            .trim();
        if s.is_empty() {
            return Ok(0);
        }
        s.parse::<usize>()
            .map_err(|_| format!("Failed parsing FCS header offset {label} ('{s}')."))
    };

    let text_start = parse_offset(&bytes[10..18], "TEXT_START")?;
    let text_end = parse_offset(&bytes[18..26], "TEXT_END")?;
    let mut data_start = parse_offset(&bytes[26..34], "DATA_START")?;
    let mut data_end = parse_offset(&bytes[34..42], "DATA_END")?;

    if text_start < 58 || text_end < text_start || text_end >= bytes.len() {
        return Err(format!(
            "Invalid FCS TEXT segment bounds ({text_start}..={text_end}) for {} bytes.",
            bytes.len()
        ));
    }

    let text_bytes = &bytes[text_start..=text_end];
    let keywords = parse_fcs_text_segment(text_bytes)?;

    // FCS 3.1 allows large >99MB data offsets in $BEGINDATA / $ENDDATA keywords when header fields are 0
    if data_start == 0 {
        if let Some(bd) = keywords.get("$BEGINDATA") {
            data_start = bd.trim().parse::<usize>().unwrap_or(0);
        }
    }
    if data_end == 0 {
        if let Some(ed) = keywords.get("$ENDDATA") {
            data_end = ed.trim().parse::<usize>().unwrap_or(0);
        }
    }

    let par: usize = keywords
        .get("$PAR")
        .ok_or_else(|| String::from("FCS TEXT segment missing required $PAR keyword."))?
        .trim()
        .parse()
        .map_err(|_| String::from("Invalid $PAR integer in FCS TEXT segment."))?;
    let tot: usize = keywords
        .get("$TOT")
        .ok_or_else(|| String::from("FCS TEXT segment missing required $TOT keyword."))?
        .trim()
        .parse()
        .map_err(|_| String::from("Invalid $TOT integer in FCS TEXT segment."))?;

    let datatype = keywords
        .get("$DATATYPE")
        .map(|s| s.trim().to_ascii_uppercase())
        .unwrap_or_else(|| String::from("F"));
    let little_endian = keywords
        .get("$BYTEORD")
        .map(|s| s.trim().starts_with("1,2"))
        .unwrap_or(true);

    let mut channels = Vec::with_capacity(par);
    for p in 1..=par {
        let short_name = keywords
            .get(&format!("$P{p}N"))
            .cloned()
            .unwrap_or_else(|| format!("P{p}"));
        let stain_name = keywords
            .get(&format!("$P{p}S"))
            .cloned()
            .unwrap_or_else(|| short_name.clone());
        let bits = keywords
            .get(&format!("$P{p}B"))
            .and_then(|s| s.trim().parse::<usize>().ok())
            .unwrap_or(32);
        let range = keywords
            .get(&format!("$P{p}R"))
            .and_then(|s| s.trim().parse::<f64>().ok())
            .unwrap_or(262_144.0);
        channels.push(FcsChannelInfo {
            index_1based: p,
            short_name,
            stain_name,
            bits,
            range,
        });
    }

    if data_start >= bytes.len() || data_end < data_start || data_end >= bytes.len() {
        return Err(format!(
            "Invalid FCS DATA segment bounds ({data_start}..={data_end})."
        ));
    }
    let payload = &bytes[data_start..=data_end];
    let mut cursor = 0usize;
    let mut events = Vec::with_capacity(tot);

    for ev in 0..tot {
        let mut row = Vec::with_capacity(par);
        for ch in &channels {
            let val = match datatype.as_str() {
                "F" => {
                    if cursor + 4 > payload.len() {
                        return Err(format!("Unexpected EOF reading f32 at event {ev}."));
                    }
                    let arr: [u8; 4] = payload[cursor..cursor + 4].try_into().unwrap();
                    cursor += 4;
                    let v = if little_endian {
                        f32::from_le_bytes(arr)
                    } else {
                        f32::from_be_bytes(arr)
                    };
                    v as f64
                }
                "D" => {
                    if cursor + 8 > payload.len() {
                        return Err(format!("Unexpected EOF reading f64 at event {ev}."));
                    }
                    let arr: [u8; 8] = payload[cursor..cursor + 8].try_into().unwrap();
                    cursor += 8;
                    if little_endian {
                        f64::from_le_bytes(arr)
                    } else {
                        f64::from_be_bytes(arr)
                    }
                }
                "I" => {
                    let n_bytes = (ch.bits / 8).clamp(1, 8);
                    if cursor + n_bytes > payload.len() {
                        return Err(format!("Unexpected EOF reading int at event {ev}."));
                    }
                    let mut acc = 0u64;
                    if little_endian {
                        for b_i in (0..n_bytes).rev() {
                            acc = (acc << 8) | (payload[cursor + b_i] as u64);
                        }
                    } else {
                        for b_i in 0..n_bytes {
                            acc = (acc << 8) | (payload[cursor + b_i] as u64);
                        }
                    }
                    cursor += n_bytes;
                    acc as f64
                }
                other => {
                    return Err(format!("Unsupported FCS $DATATYPE '{other}'."));
                }
            };
            row.push(val);
        }
        events.push(row);
    }

    let spillover = keywords
        .get("$SPILLOVER")
        .or_else(|| keywords.get("$SPILL"))
        .and_then(|s| parse_spillover_keyword(s).ok());

    Ok(FcsDataset {
        version,
        num_events: tot,
        channels,
        events,
        spillover,
        keywords,
    })
}

fn parse_fcs_text_segment(text_bytes: &[u8]) -> Result<BTreeMap<String, String>, String> {
    if text_bytes.len() < 3 {
        return Err(String::from("FCS TEXT segment is empty."));
    }
    let delim = text_bytes[0] as char;
    let raw = String::from_utf8_lossy(&text_bytes[1..]);
    let parts: Vec<&str> = raw.split(delim).collect();
    let mut map = BTreeMap::new();
    let mut i = 0usize;
    while i + 1 < parts.len() {
        let key = parts[i].trim().to_ascii_uppercase();
        let val = parts[i + 1].trim().to_string();
        if !key.is_empty() {
            map.insert(key, val);
        }
        i += 2;
    }
    Ok(map)
}

fn parse_spillover_keyword(raw: &str) -> Result<SpilloverMatrix, String> {
    let toks: Vec<&str> = raw.split(',').map(|s| s.trim()).collect();
    if toks.is_empty() {
        return Err(String::from("Empty $SPILLOVER string."));
    }
    let n: usize = toks[0]
        .parse()
        .map_err(|_| String::from("Invalid $SPILLOVER matrix dimension."))?;
    if toks.len() < 1 + n + n * n {
        return Err(format!(
            "$SPILLOVER expected {} tokens for {n}x{n} matrix, found {}.",
            1 + n + n * n,
            toks.len()
        ));
    }
    let channel_names: Vec<String> = toks[1..(1 + n)].iter().map(|s| s.to_string()).collect();
    let mut matrix = vec![vec![0.0_f64; n]; n];
    let mut k = 1 + n;
    for r in 0..n {
        for c in 0..n {
            matrix[r][c] = toks[k]
                .parse::<f64>()
                .map_err(|_| format!("Invalid $SPILLOVER entry '{}'.", toks[k]))?;
            k += 1;
        }
    }
    Ok(SpilloverMatrix {
        channel_names,
        matrix,
    })
}

/// Applies spectral spillover compensation (`x_comp = S^-1 * x_raw`) to matching channels
/// in `dataset` and returns the compensated `FcsDataset` plus the computed inverse matrix.
pub fn compensate_fcs_dataset(dataset: &FcsDataset) -> Result<(FcsDataset, Vec<Vec<f64>>), String> {
    let Some(spill) = &dataset.spillover else {
        return Err(String::from(
            "FCS dataset does not contain a $SPILLOVER or $SPILL compensation matrix.",
        ));
    };
    let n = spill.channel_names.len();
    let mut ch_indices = Vec::with_capacity(n);
    for cname in &spill.channel_names {
        let idx = dataset
            .channels
            .iter()
            .position(|c| c.short_name == *cname)
            .ok_or_else(|| format!("Spillover channel '{cname}' not found in FCS $PnN list."))?;
        ch_indices.push(idx);
    }

    // In FCS 3.1, spillover matrix S is stored such that observed = S^T * true (row = emitter, col = detector).
    // Transpose S so `A = S^T` satisfies `observed = A * true`, then invert `A` to get `true = A^-1 * observed`.
    let mut a_mat = vec![vec![0.0_f64; n]; n];
    for r in 0..n {
        for c in 0..n {
            a_mat[c][r] = spill.matrix[r][c];
        }
    }
    let comp_inv = invert_square_matrix(&a_mat)?;

    let mut out = dataset.clone();
    for ev in &mut out.events {
        let raw_sub: Vec<f64> = ch_indices.iter().map(|&ci| ev[ci]).collect();
        for (r, &dest_idx) in ch_indices.iter().enumerate() {
            let mut sum = 0.0_f64;
            for c in 0..n {
                sum += comp_inv[r][c] * raw_sub[c];
            }
            ev[dest_idx] = sum;
        }
    }
    Ok((out, comp_inv))
}

fn invert_square_matrix(a: &[Vec<f64>]) -> Result<Vec<Vec<f64>>, String> {
    let n = a.len();
    let mut aug = vec![vec![0.0_f64; 2 * n]; n];
    for i in 0..n {
        for j in 0..n {
            aug[i][j] = a[i][j];
        }
        aug[i][n + i] = 1.0;
    }
    for col in 0..n {
        let mut piv = col;
        for r in (col + 1)..n {
            if aug[r][col].abs() > aug[piv][col].abs() {
                piv = r;
            }
        }
        if aug[piv][col].abs() < 1e-12 {
            return Err(String::from("Spillover matrix is singular and cannot be inverted."));
        }
        aug.swap(col, piv);
        let d = aug[col][col];
        for j in 0..(2 * n) {
            aug[col][j] /= d;
        }
        for r in 0..n {
            if r != col {
                let fac = aug[r][col];
                for j in 0..(2 * n) {
                    aug[r][j] -= fac * aug[col][j];
                }
            }
        }
    }
    Ok((0..n).map(|i| aug[i][n..(2 * n)].to_vec()).collect())
}

/// Serializes an [`FcsDataset`] into a standards-compliant binary `FCS3.1` (`f32` little-endian) buffer.
pub fn write_fcs_31_bytes(dataset: &FcsDataset) -> Vec<u8> {
    let mut text = String::from("/");
    let _ = write!(
        text,
        "$BEGINANALYSIS/0/$ENDANALYSIS/0/$BEGINSTEXT/0/$ENDSTEXT/0/$BYTEORD/1,2,3,4/$DATATYPE/F/$MODE/L/$NEXTDATA/0/$PAR/{}/$TOT/{}/",
        dataset.channels.len(),
        dataset.events.len()
    );
    for ch in &dataset.channels {
        let p = ch.index_1based;
        let _ = write!(
            text,
            "$P{p}N/{}/$P{p}S/{}/$P{p}B/32/$P{p}R/{:.0}/$P{p}E/0,0/",
            ch.short_name, ch.stain_name, ch.range
        );
    }
    if let Some(sp) = &dataset.spillover {
        let mut sp_str = format!("{}", sp.channel_names.len());
        for nm in &sp.channel_names {
            let _ = write!(sp_str, ",{nm}");
        }
        for row in &sp.matrix {
            for &v in row {
                let _ = write!(sp_str, ",{v:.6}");
            }
        }
        let _ = write!(text, "$SPILLOVER/{sp_str}/");
    }

    let text_start = 256usize;
    let text_bytes = text.as_bytes();
    let text_end = text_start + text_bytes.len() - 1;
    let data_start = text_end + 1;
    let data_len = dataset.events.len() * dataset.channels.len() * 4;
    let data_end = data_start + data_len.saturating_sub(1);

    let mut out = vec![b' '; data_start + data_len];
    let hdr = format!(
        "FCS3.1    {:>8}{:>8}{:>8}{:>8}{:>8}{:>8}",
        text_start, text_end, data_start, data_end, 0, 0
    );
    out[..58].copy_from_slice(&hdr.as_bytes()[..58]);
    out[text_start..=text_end].copy_from_slice(text_bytes);

    let mut cursor = data_start;
    for ev in &dataset.events {
        for &val in ev {
            let b = (val as f32).to_le_bytes();
            out[cursor..cursor + 4].copy_from_slice(&b);
            cursor += 4;
        }
    }
    out
}

/// Generates a synthetic 4-channel (`FSC-A`, `SSC-A`, `FITC-A`, `PE-A`), 400-event binary `FCS 3.1`
/// dataset with 14% `FITC->PE` and 4% `PE->FITC` spectral spillover and 4 distinct cell/debris clusters.
pub fn demo_fcs_31_bytes() -> Vec<u8> {
    let channels = vec![
        FcsChannelInfo {
            index_1based: 1,
            short_name: String::from("FSC-A"),
            stain_name: String::from("Forward Scatter"),
            bits: 32,
            range: 262_144.0,
        },
        FcsChannelInfo {
            index_1based: 2,
            short_name: String::from("SSC-A"),
            stain_name: String::from("Side Scatter"),
            bits: 32,
            range: 262_144.0,
        },
        FcsChannelInfo {
            index_1based: 3,
            short_name: String::from("FITC-A"),
            stain_name: String::from("GFP-Reporter"),
            bits: 32,
            range: 262_144.0,
        },
        FcsChannelInfo {
            index_1based: 4,
            short_name: String::from("PE-A"),
            stain_name: String::from("mCherry-Reporter"),
            bits: 32,
            range: 262_144.0,
        },
    ];

    let spillover = SpilloverMatrix {
        channel_names: vec![String::from("FITC-A"), String::from("PE-A")],
        matrix: vec![vec![1.0, 0.14], vec![0.04, 1.0]],
    };

    // Deterministic LCG PRNG for reproducible cytometry clusters
    let mut seed = 0x5EED_F10C_u64.wrapping_add(141_592);
    let mut next_unit = || -> f64 {
        seed = seed.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        ((seed >> 33) as f64) / ((1u64 << 31) as f64)
    };
    let mut next_norm = || -> f64 {
        let u1 = next_unit().max(1e-9);
        let u2 = next_unit();
        (-2.0 * u1.ln()).sqrt() * (2.0 * std::f64::consts::PI * u2).cos()
    };

    let mut events = Vec::with_capacity(400);
    for i in 0..400 {
        let (fsc, ssc, true_fitc, true_pe) = if i < 60 {
            // Debris / dead-cell low scatter
            (
                22_000.0 + 4_500.0 * next_norm(),
                18_000.0 + 3_500.0 * next_norm(),
                40.0 + 90.0 * next_norm(),
                30.0 + 85.0 * next_norm(),
            )
        } else if i < 190 {
            // Double-negative live cells (includes negative compensated noise tail!)
            (
                115_000.0 + 14_000.0 * next_norm(),
                48_000.0 + 7_500.0 * next_norm(),
                25.0 + 140.0 * next_norm(),
                20.0 + 135.0 * next_norm(),
            )
        } else if i < 310 {
            // FITC+ / mCherry- single-positive population
            (
                122_000.0 + 15_000.0 * next_norm(),
                54_000.0 + 8_000.0 * next_norm(),
                18_500.0 + 3_400.0 * next_norm(),
                35.0 + 145.0 * next_norm(),
            )
        } else {
            // Double-positive FITC+ / PE+ population
            (
                132_000.0 + 15_000.0 * next_norm(),
                66_000.0 + 9_000.0 * next_norm(),
                21_000.0 + 3_800.0 * next_norm(),
                14_500.0 + 2_900.0 * next_norm(),
            )
        };

        // Apply forward spectral spillover (`obs = S^T * true`)
        let obs_fitc = true_fitc + 0.04 * true_pe;
        let obs_pe = 0.14 * true_fitc + true_pe;
        events.push(vec![fsc.max(100.0), ssc.max(100.0), obs_fitc, obs_pe]);
    }

    let ds = FcsDataset {
        version: String::from("FCS3.1"),
        num_events: events.len(),
        channels,
        events,
        spillover: Some(spillover),
        keywords: BTreeMap::new(),
    };
    write_fcs_31_bytes(&ds)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn roundtrips_fcs31_binary_and_inverts_spillover() {
        let raw = demo_fcs_31_bytes();
        let ds = parse_fcs_bytes(&raw).unwrap();
        assert_eq!(ds.version, "FCS3.1");
        assert_eq!(ds.num_events, 400);
        assert_eq!(ds.channels.len(), 4);

        let (comp, _) = compensate_fcs_dataset(&ds).unwrap();
        // FITC+ single-positive events (190..310) had +14% spillover into PE-A; after compensation PE-A median ~ 35
        let mut pe_sp: Vec<f64> = comp.events[190..310].iter().map(|r| r[3]).collect();
        pe_sp.sort_by(|a, b| a.partial_cmp(b).unwrap());
        let med_pe = pe_sp[pe_sp.len() / 2];
        assert!(med_pe.abs() < 250.0);
    }
}
