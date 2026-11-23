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

//! Shared application branding, embedded window icon decoding, and About / Copyright UI
//! for OSABT (`CRV Science Network`).

use std::sync::Arc;

use eframe::egui;
use egui::{Color32, RichText};

const OSABT_ICO_BYTES: &[u8] = include_bytes!("../assets/osabt.ico");

/// Decodes the largest 32-bit BGRA BMP-DIB entry from `assets/osabt.ico` into top-down RGBA `IconData`.
pub fn osabt_icon_data() -> Arc<egui::IconData> {
    if let Some(decoded) = decode_largest_ico_rgba(OSABT_ICO_BYTES) {
        return Arc::new(decoded);
    }
    // Minimal deterministic fallback (never panics)
    let w = 16u32;
    let h = 16u32;
    let mut rgba = vec![0u8; (w * h * 4) as usize];
    for px in rgba.chunks_exact_mut(4) {
        px[0] = 15;
        px[1] = 23;
        px[2] = 42;
        px[3] = 255;
    }
    Arc::new(egui::IconData {
        rgba,
        width: w,
        height: h,
    })
}

fn read_u16_le(bytes: &[u8], off: usize) -> Option<u16> {
    let s = bytes.get(off..off + 2)?;
    Some(u16::from_le_bytes([s[0], s[1]]))
}

fn read_u32_le(bytes: &[u8], off: usize) -> Option<u32> {
    let s = bytes.get(off..off + 4)?;
    Some(u32::from_le_bytes([s[0], s[1], s[2], s[3]]))
}

fn read_i32_le(bytes: &[u8], off: usize) -> Option<i32> {
    let s = bytes.get(off..off + 4)?;
    Some(i32::from_le_bytes([s[0], s[1], s[2], s[3]]))
}

fn decode_largest_ico_rgba(ico: &[u8]) -> Option<egui::IconData> {
    if ico.len() < 6 || read_u16_le(ico, 2)? != 1 {
        return None;
    }
    let count = read_u16_le(ico, 4)? as usize;
    let mut best: Option<(u32, u32, usize, usize)> = None;
    for i in 0..count {
        let dir_off = 6 + i * 16;
        if dir_off + 16 > ico.len() {
            break;
        }
        let w = match ico[dir_off] {
            0 => 256u32,
            v => u32::from(v),
        };
        let h = match ico[dir_off + 1] {
            0 => 256u32,
            v => u32::from(v),
        };
        let img_size = read_u32_le(ico, dir_off + 8)? as usize;
        let img_off = read_u32_le(ico, dir_off + 12)? as usize;
        if img_off.checked_add(img_size)? > ico.len() {
            continue;
        }
        if best.map(|(bw, bh, _, _)| w * h > bw * bh).unwrap_or(true) {
            best = Some((w, h, img_off, img_size));
        }
    }

    let (w, h, img_off, img_size) = best?;
    let dib = ico.get(img_off..img_off + img_size)?;
    if dib.len() < 40 {
        return None;
    }
    let header_sz = read_u32_le(dib, 0)? as usize;
    let dib_w = read_i32_le(dib, 4)?;
    let dib_h = read_i32_le(dib, 8)?;
    let bpp = read_u16_le(dib, 14)?;
    if dib_w <= 0 || dib_h <= 0 || bpp != 32 {
        return None;
    }
    let width = u32::try_from(dib_w).ok().unwrap_or(w);
    let height = u32::try_from(dib_h / 2).ok().filter(|&v| v > 0).unwrap_or(h);
    let pixel_bytes = (width as usize) * (height as usize) * 4;
    let xor = dib.get(header_sz..header_sz + pixel_bytes)?;

    let mut rgba = vec![0u8; pixel_bytes];
    let row_bytes = (width as usize) * 4;
    for y in 0..(height as usize) {
        let src_row = (height as usize - 1) - y;
        let src = &xor[src_row * row_bytes..(src_row + 1) * row_bytes];
        let dst = &mut rgba[y * row_bytes..(y + 1) * row_bytes];
        for x in 0..(width as usize) {
            let s = x * 4;
            dst[s] = src[s + 2];     // R
            dst[s + 1] = src[s + 1]; // G
            dst[s + 2] = src[s];     // B
            dst[s + 3] = src[s + 3]; // A
        }
    }

    Some(egui::IconData {
        rgba,
        width,
        height,
    })
}

/// Renders the interactive `About OSABT` modal window when `*open` is true.
pub fn draw_about_window(ctx: &egui::Context, open: &mut bool) {
    if !ctx.wants_keyboard_input()
        && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::F1))
    {
        *open = !*open;
    }
    if *open && ctx.input_mut(|i| i.consume_key(egui::Modifiers::NONE, egui::Key::Escape)) {
        *open = false;
    }
    if !*open {
        return;
    }
    egui::Window::new("About OSABT")
        .open(open)
        .collapsible(false)
        .resizable(false)
        .anchor(egui::Align2::CENTER_CENTER, egui::Vec2::ZERO)
        .show(ctx, |ui| {
            ui.set_min_width(380.0);
            ui.heading(
                RichText::new("OSABT — Open Source Ad-hoc Biology Tools")
                    .strong()
                    .color(Color32::from_rgb(56, 189, 248)),
            );
            ui.label(format!("Version {}", env!("CARGO_PKG_VERSION")));
            ui.add_space(6.0);
            ui.label("Unified 20-tool pure-Rust wet-lab, cloning, ecological, evolutionary, and biostatistical suite.");
            ui.add_space(6.0);
            ui.strong("Copyright CRV Science Network");
            ui.label("License: GNU General Public License v3.0 or later (GPL-3.0-or-later)");
            ui.add_space(8.0);
            ui.horizontal(|ui| {
                ui.hyperlink_to("CRV Science Network", "https://crvscience.com");
                ui.label("|");
                ui.hyperlink_to("CRV Science Apps", "https://crvscience.com/apps");
            });
        });
}

/// Renders a compact bottom footer with `"Copyright CRV Science Network"` and `"About"`
/// for standalone studio GUIs (`osabt plasmid`, `osabt flux`, etc.).
pub fn draw_standalone_footer(ctx: &egui::Context, about_open: &mut bool) {
    egui::TopBottomPanel::bottom("osabt_standalone_copyright_footer").show(ctx, |ui| {
        ui.horizontal(|ui| {
            ui.label(
                RichText::new("Copyright CRV Science Network • GPL-3.0-or-later")
                    .small()
                    .color(Color32::from_rgb(148, 163, 184)),
            );
            ui.with_layout(egui::Layout::right_to_left(egui::Align::Center), |ui| {
                if ui.small_button("About").clicked() {
                    *about_open = !*about_open;
                }
            });
        });
    });
    draw_about_window(ctx, about_open);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn icon_decodes_valid_rgba() {
        let icon = osabt_icon_data();
        assert_eq!(icon.width, 64);
        assert_eq!(icon.height, 64);
        assert_eq!(icon.rgba.len(), 64 * 64 * 4);
        let has_opaque = icon.rgba.chunks_exact(4).any(|p| p[3] == 255);
        let has_transparent = icon.rgba.chunks_exact(4).any(|p| p[3] == 0);
        assert!(has_opaque && has_transparent);
    }

    #[test]
    fn about_window_toggles_with_f1_and_closes_with_escape() {
        let ctx = egui::Context::default();
        let mut open = false;
        let run_key = |ctx: &egui::Context, open: &mut bool, key: egui::Key| {
            let mut raw = egui::RawInput::default();
            raw.screen_rect = Some(egui::Rect::from_min_size(
                egui::Pos2::ZERO,
                egui::vec2(1024.0, 768.0),
            ));
            raw.events.push(egui::Event::Key {
                key,
                physical_key: None,
                pressed: true,
                repeat: false,
                modifiers: egui::Modifiers::NONE,
            });
            let _ = ctx.run(raw, |c| draw_about_window(c, open));
        };

        run_key(&ctx, &mut open, egui::Key::F1);
        assert!(open, "F1 should open About OSABT modal");

        run_key(&ctx, &mut open, egui::Key::Escape);
        assert!(!open, "Escape should close About OSABT modal");

        run_key(&ctx, &mut open, egui::Key::F1);
        assert!(open, "F1 should reopen About OSABT modal");

        run_key(&ctx, &mut open, egui::Key::F1);
        assert!(!open, "Second F1 should toggle About OSABT modal closed");
    }
}
