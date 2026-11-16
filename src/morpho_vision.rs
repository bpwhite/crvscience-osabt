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

//! Tool #8 (`osabt morpho`) — Pure-Rust biological computer-vision engine:
//! 1. Zero-dependency `RgbImage` buffer with uncompressed 24/32-bit `.bmp` and Netpbm `.ppm`/`.pgm`
//!    (`P2`/`P3`/`P5`/`P6`) readers and writers.
//! 2. Otsu (1979) maximum between-class variance histogram thresholding.
//! 3. Automated Petri dish Colony-Forming Unit (`CFU`) detector with circular inner-dish ROI mask,
//!    8-connected component labeling, Chamfer `(3, 4)` Euclidean distance-transform local-peak
//!    splitting for touching doublet/triplet colonies, and Poisson MLE serial-dilution `CFU/mL`.
//! 4. Woebbecke et al. (1995) Excess Green (`ExG = 2G - R - B`) leaf-lamina segmenter, border
//!    flood-fill enclosed herbivory shot-hole detector, and Andrew's monotone-chain convex hull
//!    envelope estimator.

use std::collections::VecDeque;
use std::f64::consts::PI;
use std::fmt::Write as _;

/// Packed 8-bit RGB raster (`3 * width * height` bytes in row-major top-to-bottom order).
#[derive(Debug, Clone)]
pub struct RgbImage {
    pub width: usize,
    pub height: usize,
    pub rgb: Vec<u8>,
}

impl RgbImage {
    pub fn new(width: usize, height: usize, fill: [u8; 3]) -> Self {
        let mut rgb = Vec::with_capacity(width * height * 3);
        for _ in 0..(width * height) {
            rgb.extend_from_slice(&fill);
        }
        Self { width, height, rgb }
    }

    #[inline]
    pub fn pixel(&self, x: usize, y: usize) -> [u8; 3] {
        let i = (y * self.width + x) * 3;
        [self.rgb[i], self.rgb[i + 1], self.rgb[i + 2]]
    }

    #[inline]
    pub fn set_pixel(&mut self, x: usize, y: usize, color: [u8; 3]) {
        if x < self.width && y < self.height {
            let i = (y * self.width + x) * 3;
            self.rgb[i] = color[0];
            self.rgb[i + 1] = color[1];
            self.rgb[i + 2] = color[2];
        }
    }

    /// Loads either a Windows 24/32-bit uncompressed `.bmp` or Netpbm `P2`/`P3`/`P5`/`P6` `.ppm`/`.pgm`.
    pub fn from_bytes(bytes: &[u8]) -> Result<Self, String> {
        if bytes.len() >= 54 && &bytes[0..2] == b"BM" {
            return parse_bmp_rgb(bytes);
        }
        if bytes.len() >= 2 && bytes[0] == b'P' && matches!(bytes[1], b'2' | b'3' | b'5' | b'6') {
            return parse_pnm_rgb(bytes);
        }
        Err(String::from(
            "Unsupported raster header; provide uncompressed 24/32-bit .bmp or Netpbm .ppm/.pgm.",
        ))
    }

    /// Encodes the image as an uncompressed 24-bit Windows `.bmp` (bottom-up padded scanlines).
    pub fn to_bmp_bytes(&self) -> Vec<u8> {
        let row_raw = self.width * 3;
        let pad = (4 - (row_raw % 4)) % 4;
        let stride = row_raw + pad;
        let pixel_bytes = stride * self.height;
        let file_size = 54 + pixel_bytes;

        let mut out = Vec::with_capacity(file_size);
        out.extend_from_slice(b"BM");
        out.extend_from_slice(&(file_size as u32).to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out.extend_from_slice(&54u32.to_le_bytes());

        // BITMAPINFOHEADER (40 bytes)
        out.extend_from_slice(&40u32.to_le_bytes());
        out.extend_from_slice(&(self.width as i32).to_le_bytes());
        out.extend_from_slice(&(self.height as i32).to_le_bytes());
        out.extend_from_slice(&1u16.to_le_bytes());
        out.extend_from_slice(&24u16.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out.extend_from_slice(&(pixel_bytes as u32).to_le_bytes());
        out.extend_from_slice(&2835i32.to_le_bytes());
        out.extend_from_slice(&2835i32.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());
        out.extend_from_slice(&0u32.to_le_bytes());

        for y in (0..self.height).rev() {
            for x in 0..self.width {
                let [r, g, b] = self.pixel(x, y);
                out.extend_from_slice(&[b, g, r]);
            }
            for _ in 0..pad {
                out.push(0);
            }
        }
        out
    }

    /// Synthesizes a realistic `240 x 240` dark-field LB agar Petri dish containing 14 isolated
    /// colonies and 2 touching doublet pairs (`18` total true CFUs) plus outer plastic dish rim glare.
    pub fn demo_petri_cfu_plate() -> Self {
        let (w, h) = (240usize, 240usize);
        let mut img = Self::new(w, h, [18, 22, 28]);
        let (cx, cy) = (120.0_f64, 120.0_f64);

        for y in 0..h {
            for x in 0..w {
                let d = (((x as f64) - cx).powi(2) + ((y as f64) - cy).powi(2)).sqrt();
                if d <= 104.0 {
                    img.set_pixel(x, y, [32, 42, 46]); // Agar interior
                } else if (104.0..=110.0).contains(&d) {
                    img.set_pixel(x, y, [195, 202, 210]); // Bright plastic rim glare (must be masked!)
                }
            }
        }

        // 14 isolated single colonies + 2 touching doublets = 18 total CFUs
        let centers: &[(f64, f64, f64)] = &[
            (80.0, 74.0, 5.2),
            (112.0, 62.0, 4.8),
            (150.0, 76.0, 5.4),
            (68.0, 108.0, 4.9),
            (98.0, 102.0, 5.1),
            (134.0, 98.0, 4.6),
            (168.0, 112.0, 5.0),
            (74.0, 144.0, 5.3),
            (108.0, 140.0, 4.7),
            (144.0, 148.0, 5.2),
            (172.0, 142.0, 4.5),
            (92.0, 174.0, 5.0),
            (126.0, 178.0, 4.9),
            (154.0, 172.0, 4.8),
            // Touching doublet #1
            (116.0, 118.0, 5.6),
            (124.5, 120.0, 5.5),
            // Touching doublet #2
            (86.0, 124.0, 5.4),
            (93.8, 127.0, 5.3),
        ];

        for &(px, py, rad) in centers {
            let x0 = (px - rad - 2.0).max(0.0) as usize;
            let x1 = ((px + rad + 2.0) as usize).min(w - 1);
            let y0 = (py - rad - 2.0).max(0.0) as usize;
            let y1 = ((py + rad + 2.0) as usize).min(h - 1);
            for y in y0..=y1 {
                for x in x0..=x1 {
                    let dist = (((x as f64) - px).powi(2) + ((y as f64) - py).powi(2)).sqrt();
                    if dist <= rad {
                        let profile = (1.0 - (dist / (rad + 0.8)).powi(2)).clamp(0.2, 1.0);
                        let cur = img.pixel(x, y);
                        let r = (45.0 + 200.0 * profile) as u8;
                        let g = (48.0 + 192.0 * profile) as u8;
                        let b = (40.0 + 160.0 * profile) as u8;
                        img.set_pixel(x, y, [cur[0].max(r), cur[1].max(g), cur[2].max(b)]);
                    }
                }
            }
        }
        img
    }

    /// Synthesizes a `260 x 200` backlit botanical leaf image (`px_per_mm = 4.0`) containing
    /// 3 enclosed interior beetle shot-holes and 1 marginal caterpillar notch.
    pub fn demo_herbivory_leaf() -> Self {
        let (w, h) = (260usize, 200usize);
        let mut img = Self::new(w, h, [242, 244, 246]);

        let enclosed_holes: &[(f64, f64, f64)] = &[
            (112.0, 86.0, 8.5),
            (148.0, 108.0, 10.0),
            (126.0, 122.0, 6.5),
        ];
        let marginal_bite = (165.0_f64, 58.0_f64, 15.0_f64);

        for y in 0..h {
            for x in 0..w {
                // Ovate-acuminate lamina centered around (130, 100), half-length 86 px, half-width 52 px
                let u = ((x as f64) - 130.0) / 86.0;
                let v = ((y as f64) - 100.0) / (52.0 * (1.0 - 0.28 * u));
                if u.abs() <= 1.0 && (u * u + v * v) <= 1.0 {
                    // Check marginal caterpillar bite
                    let db = (((x as f64) - marginal_bite.0).powi(2)
                        + ((y as f64) - marginal_bite.1).powi(2))
                    .sqrt();
                    if db <= marginal_bite.2 {
                        continue;
                    }
                    // Check enclosed beetle shot-holes
                    let in_hole = enclosed_holes.iter().any(|&(hx, hy, hr)| {
                        (((x as f64) - hx).powi(2) + ((y as f64) - hy).powi(2)).sqrt() <= hr
                    });
                    if !in_hole {
                        img.set_pixel(x, y, [46, 154, 62]);
                    }
                }
            }
        }
        img
    }
}

fn parse_bmp_rgb(bytes: &[u8]) -> Result<RgbImage, String> {
    let data_offset = u32::from_le_bytes(bytes[10..14].try_into().unwrap()) as usize;
    let width_i = i32::from_le_bytes(bytes[18..22].try_into().unwrap());
    let height_i = i32::from_le_bytes(bytes[22..26].try_into().unwrap());
    let bpp = u16::from_le_bytes(bytes[28..30].try_into().unwrap());
    let compression = u32::from_le_bytes(bytes[30..34].try_into().unwrap());

    if width_i <= 0 || height_i == 0 || (bpp != 24 && bpp != 32) || compression != 0 {
        return Err(String::from(
            "Only uncompressed (BI_RGB) 24-bit or 32-bit BMP images are supported.",
        ));
    }

    let width = width_i as usize;
    let top_down = height_i < 0;
    let height = height_i.unsigned_abs() as usize;
    let bytes_per_px = (bpp / 8) as usize;
    let stride = (width * bytes_per_px + 3) & !3;

    if data_offset + stride * height > bytes.len() {
        return Err(String::from("Truncated BMP pixel payload."));
    }

    let mut img = RgbImage::new(width, height, [0, 0, 0]);
    for row in 0..height {
        let dst_y = if top_down { row } else { height - 1 - row };
        let src_off = data_offset + row * stride;
        for x in 0..width {
            let p = src_off + x * bytes_per_px;
            let (b, g, r) = (bytes[p], bytes[p + 1], bytes[p + 2]);
            img.set_pixel(x, dst_y, [r, g, b]);
        }
    }
    Ok(img)
}

fn parse_pnm_rgb(bytes: &[u8]) -> Result<RgbImage, String> {
    let magic = bytes[1];
    let mut idx = 2usize;

    fn next_token(buf: &[u8], idx: &mut usize) -> Option<String> {
        while *idx < buf.len() {
            if buf[*idx] == b'#' {
                while *idx < buf.len() && buf[*idx] != b'\n' {
                    *idx += 1;
                }
            } else if buf[*idx].is_ascii_whitespace() {
                *idx += 1;
            } else {
                break;
            }
        }
        if *idx >= buf.len() {
            return None;
        }
        let start = *idx;
        while *idx < buf.len() && !buf[*idx].is_ascii_whitespace() && buf[*idx] != b'#' {
            *idx += 1;
        }
        std::str::from_utf8(&buf[start..*idx]).ok().map(String::from)
    }

    let width: usize = next_token(bytes, &mut idx)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| String::from("Invalid PNM width."))?;
    let height: usize = next_token(bytes, &mut idx)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| String::from("Invalid PNM height."))?;
    let maxval: u32 = next_token(bytes, &mut idx)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| String::from("Invalid PNM maxval."))?;
    if maxval == 0 {
        return Err(String::from("PNM maxval must be > 0."));
    }
    if idx < bytes.len() && bytes[idx].is_ascii_whitespace() {
        idx += 1; // Consume single whitespace byte before binary payload
    }

    let mut img = RgbImage::new(width, height, [0, 0, 0]);
    match magic {
        b'6' => {
            let need = width * height * 3;
            if idx + need > bytes.len() {
                return Err(String::from("Truncated binary PPM (P6) data."));
            }
            img.rgb.copy_from_slice(&bytes[idx..(idx + need)]);
        }
        b'5' => {
            let need = width * height;
            if idx + need > bytes.len() {
                return Err(String::from("Truncated binary PGM (P5) data."));
            }
            for i in 0..need {
                let g = bytes[idx + i];
                img.rgb[i * 3] = g;
                img.rgb[i * 3 + 1] = g;
                img.rgb[i * 3 + 2] = g;
            }
        }
        b'3' => {
            for y in 0..height {
                for x in 0..width {
                    let r = next_token(bytes, &mut idx).and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                    let g = next_token(bytes, &mut idx).and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                    let b = next_token(bytes, &mut idx).and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                    img.set_pixel(
                        x,
                        y,
                        [
                            ((r * 255) / maxval) as u8,
                            ((g * 255) / maxval) as u8,
                            ((b * 255) / maxval) as u8,
                        ],
                    );
                }
            }
        }
        b'2' => {
            for y in 0..height {
                for x in 0..width {
                    let g = next_token(bytes, &mut idx).and_then(|s| s.parse::<u32>().ok()).unwrap_or(0);
                    let v = ((g * 255) / maxval) as u8;
                    img.set_pixel(x, y, [v, v, v]);
                }
            }
        }
        _ => unreachable!(),
    }
    Ok(img)
}

/// Computes Otsu's (1979) exact optimal gray-level threshold (`0..=255`) maximizing between-class variance.
pub fn otsu_threshold_u8(samples: &[u8]) -> u8 {
    if samples.is_empty() {
        return 128;
    }
    let mut hist = [0usize; 256];
    for &v in samples {
        hist[v as usize] += 1;
    }
    let total = samples.len() as f64;
    let sum_all: f64 = hist
        .iter()
        .enumerate()
        .map(|(i, &c)| (i as f64) * (c as f64))
        .sum();

    let mut w_bg = 0.0_f64;
    let mut sum_bg = 0.0_f64;
    let mut best_var = -1.0_f64;
    let mut best_t = 128u8;

    for (t, &cnt) in hist.iter().enumerate() {
        w_bg += cnt as f64;
        if w_bg == 0.0 {
            continue;
        }
        let w_fg = total - w_bg;
        if w_fg == 0.0 {
            break;
        }
        sum_bg += (t as f64) * (cnt as f64);
        let m_bg = sum_bg / w_bg;
        let m_fg = (sum_all - sum_bg) / w_fg;
        let between = w_bg * w_fg * (m_bg - m_fg).powi(2);
        if between > best_var {
            best_var = between;
            best_t = t as u8;
        }
    }
    best_t
}

/// Detected Colony-Forming Unit (`CFU`) centroid and equivalent radius.
#[derive(Debug, Clone)]
pub struct DetectedColony {
    pub x_px: f64,
    pub y_px: f64,
    pub equivalent_radius_px: f64,
    pub area_px: usize,
    pub split_from_cluster: bool,
}

/// Complete Petri dish CFU counting output.
#[derive(Debug, Clone)]
pub struct PetriCfuReport {
    pub total_cfu_count: usize,
    pub isolated_colonies: usize,
    pub split_cluster_colonies: usize,
    pub threshold_used: u8,
    pub mean_colony_diameter_px: f64,
    pub colonies: Vec<DetectedColony>,
}

/// Counts bacterial/yeast colonies inside a circular Petri dish ROI (`roi_radius_fraction` of half-min-dim,
/// e.g. `0.84` to mask rim glare), splitting touching colonies via Chamfer `(3, 4)` distance peaks.
pub fn count_petri_colonies(
    img: &RgbImage,
    dark_colonies_on_light_agar: bool,
    roi_radius_fraction: f64,
    min_colony_area_px: usize,
    threshold_override: Option<u8>,
) -> PetriCfuReport {
    let (w, h) = (img.width, img.height);
    let (cx, cy) = ((w as f64) * 0.5, (h as f64) * 0.5);
    let max_r = (w.min(h) as f64) * 0.5 * roi_radius_fraction.clamp(0.3, 0.99);
    let max_r2 = max_r * max_r;

    let mut gray = vec![0u8; w * h];
    let mut in_roi = vec![false; w * h];
    let mut roi_vals = Vec::with_capacity(w * h);

    for y in 0..h {
        for x in 0..w {
            let idx = y * w + x;
            let d2 = ((x as f64) - cx).powi(2) + ((y as f64) - cy).powi(2);
            if d2 <= max_r2 {
                in_roi[idx] = true;
                let [r, g, b] = img.pixel(x, y);
                let luma = ((77u32 * r as u32 + 150 * g as u32 + 29 * b as u32) >> 8) as u8;
                let sig = if dark_colonies_on_light_agar {
                    255 - luma
                } else {
                    luma
                };
                gray[idx] = sig;
                roi_vals.push(sig);
            }
        }
    }

    let thresh = threshold_override.unwrap_or_else(|| otsu_threshold_u8(&roi_vals).max(25));
    let mut fg = vec![false; w * h];
    for i in 0..(w * h) {
        if in_roi[i] && gray[i] > thresh {
            fg[i] = true;
        }
    }

    // Chamfer (3, 4) integer distance transform (scaled by 3 so 1 orthogonal pixel = 3 units)
    let mut dist = vec![0i32; w * h];
    let inf = 1_000_000i32;
    for i in 0..(w * h) {
        dist[i] = if fg[i] { inf } else { 0 };
    }
    for y in 1..(h.saturating_sub(1)) {
        for x in 1..(w.saturating_sub(1)) {
            let i = y * w + x;
            if dist[i] > 0 {
                let m = (dist[(y - 1) * w + x] + 3)
                    .min(dist[y * w + (x - 1)] + 3)
                    .min(dist[(y - 1) * w + (x - 1)] + 4)
                    .min(dist[(y - 1) * w + (x + 1)] + 4);
                dist[i] = dist[i].min(m);
            }
        }
    }
    for y in (1..(h.saturating_sub(1))).rev() {
        for x in (1..(w.saturating_sub(1))).rev() {
            let i = y * w + x;
            if dist[i] > 0 {
                let m = (dist[(y + 1) * w + x] + 3)
                    .min(dist[y * w + (x + 1)] + 3)
                    .min(dist[(y + 1) * w + (x + 1)] + 4)
                    .min(dist[(y + 1) * w + (x - 1)] + 4);
                dist[i] = dist[i].min(m);
            }
        }
    }

    // Connected components + local distance maxima peak splitting
    let mut visited = vec![false; w * h];
    let mut colonies = Vec::new();
    let mut isolated_cnt = 0usize;
    let mut split_cnt = 0usize;

    for y in 1..(h.saturating_sub(1)) {
        for x in 1..(w.saturating_sub(1)) {
            let start = y * w + x;
            if !fg[start] || visited[start] {
                continue;
            }

            let mut comp = Vec::new();
            let mut q = VecDeque::new();
            visited[start] = true;
            q.push_back((x, y));

            while let Some((cx_i, cy_i)) = q.pop_front() {
                comp.push((cx_i, cy_i));
                for dy in -1isize..=1 {
                    for dx in -1isize..=1 {
                        if dx == 0 && dy == 0 {
                            continue;
                        }
                        let nx = (cx_i as isize + dx) as usize;
                        let ny = (cy_i as isize + dy) as usize;
                        if nx < w && ny < h {
                            let ni = ny * w + nx;
                            if fg[ni] && !visited[ni] {
                                visited[ni] = true;
                                q.push_back((nx, ny));
                            }
                        }
                    }
                }
            }

            if comp.len() < min_colony_area_px.max(4) {
                continue;
            }

            // Find local maxima in `dist` separated by >= 5 pixels (`dist >= 9`, i.e., >= 3 px radius)
            let mut candidates: Vec<(usize, usize, i32)> = comp
                .iter()
                .filter_map(|&(px, py)| {
                    let d0 = dist[py * w + px];
                    if d0 < 9 {
                        return None;
                    }
                    for oy in -2isize..=2 {
                        for ox in -2isize..=2 {
                            let tx = (px as isize + ox).clamp(0, (w - 1) as isize) as usize;
                            let ty = (py as isize + oy).clamp(0, (h - 1) as isize) as usize;
                            if dist[ty * w + tx] > d0 {
                                return None;
                            }
                        }
                    }
                    Some((px, py, d0))
                })
                .collect();

            candidates.sort_by(|a, b| b.2.cmp(&a.2));
            let mut peaks: Vec<(usize, usize, i32)> = Vec::new();
            for (px, py, d0) in candidates {
                let min_sep2 = 30.0_f64; // ~5.5 px minimum center-to-center distance
                let too_close = peaks.iter().any(|&(kx, ky, _)| {
                    ((px as f64) - (kx as f64)).powi(2) + ((py as f64) - (ky as f64)).powi(2)
                        < min_sep2
                });
                if !too_close {
                    peaks.push((px, py, d0));
                }
            }

            if peaks.len() <= 1 {
                let sum_x: f64 = comp.iter().map(|p| p.0 as f64).sum();
                let sum_y: f64 = comp.iter().map(|p| p.1 as f64).sum();
                let area = comp.len();
                isolated_cnt += 1;
                colonies.push(DetectedColony {
                    x_px: sum_x / (area as f64),
                    y_px: sum_y / (area as f64),
                    equivalent_radius_px: ((area as f64) / PI).sqrt(),
                    area_px: area,
                    split_from_cluster: false,
                });
            } else {
                let sub_area = (comp.len() / peaks.len()).max(1);
                split_cnt += peaks.len();
                for (px, py, d0) in peaks {
                    colonies.push(DetectedColony {
                        x_px: px as f64,
                        y_px: py as f64,
                        equivalent_radius_px: (d0 as f64) / 3.0,
                        area_px: sub_area,
                        split_from_cluster: true,
                    });
                }
            }
        }
    }

    let mean_diam = if colonies.is_empty() {
        0.0
    } else {
        colonies
            .iter()
            .map(|c| 2.0 * c.equivalent_radius_px)
            .sum::<f64>()
            / (colonies.len() as f64)
    };

    PetriCfuReport {
        total_cfu_count: colonies.len(),
        isolated_colonies: isolated_cnt,
        split_cluster_colonies: split_cnt,
        threshold_used: thresh,
        mean_colony_diameter_px: mean_diam,
        colonies,
    }
}

/// Computes Poisson Maximum-Likelihood `CFU/mL` and 95% confidence bounds across one or more
/// `(colony_count, plated_volume_ml, dilution_factor)` observations (e.g. `dilution_factor = 1e-5`).
pub fn estimate_cfu_per_ml(observations: &[(usize, f64, f64)]) -> Result<(f64, f64, f64), String> {
    let mut sum_c = 0usize;
    let mut sum_eff_vol = 0.0_f64;
    for &(cnt, vol_ml, dil) in observations {
        if vol_ml <= 0.0 || dil <= 0.0 {
            return Err(String::from(
                "Plated volume (mL) and dilution factor (e.g. 1e-5) must be > 0.",
            ));
        }
        sum_c += cnt;
        sum_eff_vol += vol_ml * dil;
    }
    if sum_eff_vol <= 0.0 {
        return Err(String::from("At least one valid plate observation is required."));
    }

    let c = sum_c as f64;
    let mle = c / sum_eff_vol;
    // Byers / score-interval Poisson 95% bounds on total count `c`
    let z = 1.959_963_984_540_054;
    let lo_c = if sum_c == 0 {
        0.0
    } else {
        c * (1.0 - 1.0 / (9.0 * c) - z / (3.0 * c.sqrt())).powi(3)
    };
    let hi_c = (c + 1.0)
        * (1.0 - 1.0 / (9.0 * (c + 1.0)) + z / (3.0 * (c + 1.0).sqrt())).powi(3);

    Ok((mle, lo_c / sum_eff_vol, hi_c / sum_eff_vol))
}

/// Metrics produced by [`analyze_leaf_herbivory`].
#[derive(Debug, Clone)]
pub struct LeafHerbivoryReport {
    pub green_lamina_px: usize,
    pub enclosed_hole_px: usize,
    pub enclosed_hole_count: usize,
    pub convex_hull_px: f64,
    pub lamina_area_cm2: f64,
    pub enclosed_herbivory_cm2: f64,
    /// Estimated total intact leaf area (`lamina + enclosed_holes`, or convex hull envelope).
    pub convex_envelope_cm2: f64,
    /// `enclosed_hole_px / (green_lamina_px + enclosed_hole_px) * 100`.
    pub interior_herbivory_pct: f64,
    /// `(convex_hull_px - green_lamina_px) / convex_hull_px * 100` capturing both marginal notches + interior holes.
    pub total_defoliation_envelope_pct: f64,
}

/// Segments green leaf tissue via Woebbecke et al. (1995) Excess Green (`ExG = 2G - R - B`),
/// isolates interior non-border herbivory holes via 4-connected exterior flood fill, and computes
/// convex-hull marginal + interior defoliation calibrated by `px_per_mm`.
pub fn analyze_leaf_herbivory(
    img: &RgbImage,
    px_per_mm: f64,
    exg_threshold_override: Option<i32>,
) -> Result<LeafHerbivoryReport, String> {
    if px_per_mm <= 0.0 {
        return Err(String::from("Calibration px_per_mm must be > 0."));
    }
    let (w, h) = (img.width, img.height);
    if w < 4 || h < 4 {
        return Err(String::from("Leaf image must be at least 4x4 pixels."));
    }

    let mut exg_shifted = Vec::with_capacity(w * h);
    let mut raw_exg = Vec::with_capacity(w * h);
    for y in 0..h {
        for x in 0..w {
            let [r, g, b] = img.pixel(x, y);
            let exg = 2 * (g as i32) - (r as i32) - (b as i32);
            raw_exg.push(exg);
            // Map `-128..=382` into `0..=255` for automatic Otsu threshold selection
            exg_shifted.push(((exg + 128) / 2).clamp(0, 255) as u8);
        }
    }

    let exg_thresh = exg_threshold_override.unwrap_or_else(|| {
        let t_u8 = otsu_threshold_u8(&exg_shifted) as i32;
        (t_u8 * 2 - 128).max(20)
    });

    let mut is_leaf = vec![false; w * h];
    let mut leaf_pts: Vec<(f64, f64)> = Vec::new();
    for y in 0..h {
        for x in 0..w {
            let i = y * w + x;
            if raw_exg[i] >= exg_thresh {
                is_leaf[i] = true;
                leaf_pts.push((x as f64, y as f64));
            }
        }
    }

    // Flood-fill exterior non-leaf background starting from all 4 image borders
    let mut exterior_bg = vec![false; w * h];
    let mut q = VecDeque::new();
    for x in 0..w {
        for y in [0, h - 1] {
            let i = y * w + x;
            if !is_leaf[i] && !exterior_bg[i] {
                exterior_bg[i] = true;
                q.push_back((x, y));
            }
        }
    }
    for y in 0..h {
        for x in [0, w - 1] {
            let i = y * w + x;
            if !is_leaf[i] && !exterior_bg[i] {
                exterior_bg[i] = true;
                q.push_back((x, y));
            }
        }
    }

    while let Some((cx, cy)) = q.pop_front() {
        for (dx, dy) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
            let nx = cx as isize + dx;
            let ny = cy as isize + dy;
            if nx >= 0 && ny >= 0 && (nx as usize) < w && (ny as usize) < h {
                let (ux, uy) = (nx as usize, ny as usize);
                let ni = uy * w + ux;
                if !is_leaf[ni] && !exterior_bg[ni] {
                    exterior_bg[ni] = true;
                    q.push_back((ux, uy));
                }
            }
        }
    }

    // Any pixel with `!is_leaf[i] && !exterior_bg[i]` is an enclosed interior herbivory hole!
    let mut hole_visited = vec![false; w * h];
    let mut hole_px = 0usize;
    let mut hole_count = 0usize;

    for y in 1..(h - 1) {
        for x in 1..(w - 1) {
            let i = y * w + x;
            if !is_leaf[i] && !exterior_bg[i] && !hole_visited[i] {
                let mut sz = 0usize;
                let mut hq = VecDeque::new();
                hole_visited[i] = true;
                hq.push_back((x, y));
                while let Some((hx, hy)) = hq.pop_front() {
                    sz += 1;
                    for (dx, dy) in [(-1isize, 0isize), (1, 0), (0, -1), (0, 1)] {
                        let ux = (hx as isize + dx) as usize;
                        let uy = (hy as isize + dy) as usize;
                        let ni = uy * w + ux;
                        if !is_leaf[ni] && !exterior_bg[ni] && !hole_visited[ni] {
                            hole_visited[ni] = true;
                            hq.push_back((ux, uy));
                        }
                    }
                }
                if sz >= 4 {
                    hole_px += sz;
                    hole_count += 1;
                }
            }
        }
    }

    let lamina_px = leaf_pts.len();
    let hull_px = convex_hull_area(&mut leaf_pts).max((lamina_px + hole_px) as f64);
    let px_per_cm2 = (px_per_mm * 10.0).powi(2);

    let int_denom = (lamina_px + hole_px).max(1) as f64;
    let int_pct = (hole_px as f64 / int_denom) * 100.0;
    let env_pct = if hull_px > 0.0 {
        ((hull_px - (lamina_px as f64)).max(0.0) / hull_px) * 100.0
    } else {
        0.0
    };

    Ok(LeafHerbivoryReport {
        green_lamina_px: lamina_px,
        enclosed_hole_px: hole_px,
        enclosed_hole_count: hole_count,
        convex_hull_px: hull_px,
        lamina_area_cm2: (lamina_px as f64) / px_per_cm2,
        enclosed_herbivory_cm2: (hole_px as f64) / px_per_cm2,
        convex_envelope_cm2: hull_px / px_per_cm2,
        interior_herbivory_pct: int_pct,
        total_defoliation_envelope_pct: env_pct,
    })
}

fn convex_hull_area(pts: &mut [(f64, f64)]) -> f64 {
    if pts.len() < 3 {
        return 0.0;
    }
    pts.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap()
            .then_with(|| a.1.partial_cmp(&b.1).unwrap())
    });

    let cross2 = |o: (f64, f64), a: (f64, f64), b: (f64, f64)| -> f64 {
        (a.0 - o.0) * (b.1 - o.1) - (a.1 - o.1) * (b.0 - o.0)
    };

    let mut lower: Vec<(f64, f64)> = Vec::new();
    for &p in pts.iter() {
        while lower.len() >= 2 && cross2(lower[lower.len() - 2], lower[lower.len() - 1], p) <= 0.0 {
            lower.pop();
        }
        lower.push(p);
    }

    let mut upper: Vec<(f64, f64)> = Vec::new();
    for &p in pts.iter().rev() {
        while upper.len() >= 2 && cross2(upper[upper.len() - 2], upper[upper.len() - 1], p) <= 0.0 {
            upper.pop();
        }
        upper.push(p);
    }
    lower.pop();
    upper.pop();
    lower.extend(upper);

    let mut area2 = 0.0_f64;
    for i in 0..lower.len() {
        let j = (i + 1) % lower.len();
        area2 += lower[i].0 * lower[j].1 - lower[j].0 * lower[i].1;
    }
    0.5 * area2.abs()
}

/// Renders an annotated vector `.svg` overlay of detected Petri dish CFUs.
pub fn render_petri_cfu_svg(img: &RgbImage, rep: &PetriCfuReport) -> String {
    let mut svg = String::with_capacity(4096);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {} {}\" width=\"{}\" height=\"{}\">",
        img.width, img.height, img.width, img.height
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{}\" height=\"{}\" fill=\"#0f172a\"/>",
        img.width, img.height
    );
    let _ = writeln!(
        svg,
        "  <circle cx=\"{:.1}\" cy=\"{:.1}\" r=\"{:.1}\" fill=\"#1e293b\" stroke=\"#64748b\" stroke-width=\"2\"/>",
        img.width as f64 * 0.5,
        img.height as f64 * 0.5,
        (img.width.min(img.height) as f64) * 0.42
    );
    for c in &rep.colonies {
        let stroke = if c.split_from_cluster {
            "#f59e0b"
        } else {
            "#22c55e"
        };
        let _ = writeln!(
            svg,
            "  <circle cx=\"{:.2}\" cy=\"{:.2}\" r=\"{:.2}\" fill=\"none\" stroke=\"{stroke}\" stroke-width=\"1.6\"/>",
            c.x_px,
            c.y_px,
            c.equivalent_radius_px.max(2.5)
        );
    }
    let _ = writeln!(
        svg,
        "  <text x=\"10\" y=\"18\" font-family=\"sans-serif\" font-size=\"12\" font-weight=\"700\" fill=\"#f8fafc\">Total CFUs: {} ({} split from touching clusters)</text>",
        rep.total_cfu_count, rep.split_cluster_colonies
    );
    svg.push_str("</svg>\n");
    svg
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_bmp_splits_touching_cfus_and_measures_leaf_holes() {
        let plate = RgbImage::demo_petri_cfu_plate();
        let bmp = plate.to_bmp_bytes();
        let reloaded = RgbImage::from_bytes(&bmp).unwrap();
        assert_eq!(reloaded.width, 240);

        let rep = count_petri_colonies(&reloaded, false, 0.84, 8, None);
        // 14 isolated + 2 doublets = 18 true colonies
        assert_eq!(rep.total_cfu_count, 18);
        assert!(rep.split_cluster_colonies >= 4);

        let leaf = RgbImage::demo_herbivory_leaf();
        let lrep = analyze_leaf_herbivory(&leaf, 4.0, None).unwrap();
        assert_eq!(lrep.enclosed_hole_count, 3);
        assert!(lrep.interior_herbivory_pct > 2.0 && lrep.interior_herbivory_pct < 15.0);
        assert!(lrep.total_defoliation_envelope_pct > lrep.interior_herbivory_pct);
    }
}
