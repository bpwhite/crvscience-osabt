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

//! Shared native OS file dialogs (`rfd`), surfaced file I/O error reporting, and reusable
//! clipboard/paste import modal helpers for all OSABT GUI panels and studio windows.

use std::fs;
use std::path::{Path, PathBuf};

use eframe::egui;
use egui::{Color32, RichText};

/// Opens a native OS file picker for loading scientific datasets.
pub fn pick_open_file(title: &str, filter_name: &str, exts: &[&str]) -> Option<PathBuf> {
    let mut dlg = rfd::FileDialog::new().set_title(title);
    if !exts.is_empty() {
        dlg = dlg.add_filter(filter_name, exts);
    }
    dlg.pick_file()
}

/// Opens a native OS file save dialog with a suggested filename.
pub fn pick_save_file(
    title: &str,
    default_name: &str,
    filter_name: &str,
    exts: &[&str],
) -> Option<PathBuf> {
    let mut dlg = rfd::FileDialog::new()
        .set_title(title)
        .set_file_name(default_name);
    if !exts.is_empty() {
        dlg = dlg.add_filter(filter_name, exts);
    }
    dlg.save_file()
}

/// Opens a native OS folder picker for multi-file bundle exports.
pub fn pick_folder(title: &str) -> Option<PathBuf> {
    rfd::FileDialog::new().set_title(title).pick_folder()
}

/// Reads UTF-8 text from `path` and surfaces any I/O failure into `status` + `is_error`.
pub fn read_text_reported(path: &Path, status: &mut String, is_error: &mut bool) -> Option<String> {
    match fs::read_to_string(path) {
        Ok(text) => {
            *is_error = false;
            Some(text.trim_start_matches('\u{feff}').to_string())
        }
        Err(err) => {
            *is_error = true;
            *status = format!("Read error '{}': {err}", path.display());
            None
        }
    }
}

/// Reads raw bytes from `path` and surfaces any I/O failure into `status` + `is_error`.
pub fn read_bytes_reported(path: &Path, status: &mut String, is_error: &mut bool) -> Option<Vec<u8>> {
    match fs::read(path) {
        Ok(bytes) => {
            *is_error = false;
            Some(bytes)
        }
        Err(err) => {
            *is_error = true;
            *status = format!("Read error '{}': {err}", path.display());
            None
        }
    }
}

/// Writes UTF-8 text to `path`, updating `status` and `is_error` on success or failure.
pub fn write_text_reported(
    path: &Path,
    content: &str,
    status: &mut String,
    is_error: &mut bool,
    ok_label: &str,
) -> bool {
    write_bytes_reported(path, content.as_bytes(), status, is_error, ok_label)
}

/// Writes raw bytes to `path`, updating `status` and `is_error` on success or failure.
pub fn write_bytes_reported(
    path: &Path,
    bytes: &[u8],
    status: &mut String,
    is_error: &mut bool,
    ok_label: &str,
) -> bool {
    match fs::write(path, bytes) {
        Ok(()) => {
            *is_error = false;
            *status = format!("{ok_label} -> {}", path.display());
            true
        }
        Err(err) => {
            *is_error = true;
            *status = format!("Write error '{}': {err}", path.display());
            false
        }
    }
}

/// Renders a status message colored red on errors and light slate on normal status updates.
pub fn draw_status_banner(ui: &mut egui::Ui, status_msg: &str, is_error: bool) {
    let col = if is_error {
        Color32::from_rgb(248, 113, 113)
    } else {
        Color32::from_rgb(203, 213, 225)
    };
    ui.label(RichText::new(status_msg).color(col));
}

/// Reusable modal state for pasting CSV/TSV/FASTA/FASTQ/Newick user data into any panel.
#[derive(Debug, Clone, Default)]
pub struct PasteModalState {
    pub open: bool,
    pub title: String,
    pub hint: String,
    pub buffer: String,
    pub target_slot: usize,
}

impl PasteModalState {
    pub fn open_for(&mut self, slot: usize, title: &str, hint: &str, prefill: &str) {
        self.open = true;
        self.target_slot = slot;
        self.title = title.to_string();
        self.hint = hint.to_string();
        self.buffer = prefill.to_string();
    }
}

/// Draws `modal` if `modal.open` is true and returns `Some((target_slot, pasted_text))` when applied.
pub fn show_paste_modal(
    ctx: &egui::Context,
    modal: &mut PasteModalState,
) -> Option<(usize, String)> {
    if !modal.open {
        return None;
    }
    if ctx.input(|i| i.key_pressed(egui::Key::Escape)) {
        modal.open = false;
        return None;
    }
    let mut applied: Option<(usize, String)> = None;
    let mut open = modal.open;
    let title = if modal.title.is_empty() {
        "Paste Input Data".to_string()
    } else {
        modal.title.clone()
    };

    egui::Window::new(title)
        .open(&mut open)
        .collapsible(false)
        .resizable(true)
        .default_width(540.0)
        .anchor(egui::Align2::CENTER_CENTER, egui::Vec2::ZERO)
        .show(ctx, |ui| {
            if !modal.hint.is_empty() {
                ui.label(RichText::new(&modal.hint).small().color(Color32::from_rgb(148, 163, 184)));
                ui.add_space(4.0);
            }
            egui::ScrollArea::vertical()
                .max_height(260.0)
                .show(ui, |ui| {
                    ui.add(
                        egui::TextEdit::multiline(&mut modal.buffer)
                            .code_editor()
                            .desired_rows(12)
                            .desired_width(f32::INFINITY),
                    );
                });
            ui.add_space(6.0);
            ui.horizontal(|ui| {
                if ui.button("Apply Data").clicked() {
                    applied = Some((modal.target_slot, modal.buffer.clone()));
                }
                if ui.button("Cancel").clicked() {
                    modal.open = false;
                }
            });
        });

    if applied.is_some() {
        modal.open = false;
    } else {
        modal.open = open;
    }
    applied
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn write_text_reported_handles_success_and_error() {
        let tmp = std::env::temp_dir().join(format!("osabt_gui_io_test_{}", std::process::id()));
        let _ = fs::create_dir_all(&tmp);
        let ok_path = tmp.join("sample.txt");

        let mut status = String::new();
        let mut is_err = true;
        let ok = write_text_reported(&ok_path, "hello osabt", &mut status, &mut is_err, "Saved sample");
        assert!(ok);
        assert!(!is_err);
        assert!(status.contains("Saved sample"));
        assert_eq!(fs::read_to_string(&ok_path).unwrap(), "hello osabt");

        let bad_path = tmp.join("nonexistent_subdir_12345").join("nested.txt");
        let ok2 = write_text_reported(&bad_path, "fail", &mut status, &mut is_err, "Saved");
        assert!(!ok2);
        assert!(is_err);
        assert!(status.starts_with("Write error"));

        let _ = fs::remove_dir_all(&tmp);
    }

    #[test]
    fn paste_modal_closes_on_escape() {
        let ctx = egui::Context::default();
        let mut modal = PasteModalState::default();
        modal.open_for(0, "Paste Test", "Hint", "abc");
        assert!(modal.open);

        let mut raw = egui::RawInput::default();
        raw.screen_rect = Some(egui::Rect::from_min_size(
            egui::Pos2::ZERO,
            egui::vec2(1024.0, 768.0),
        ));
        raw.events.push(egui::Event::Key {
            key: egui::Key::Escape,
            physical_key: None,
            pressed: true,
            repeat: false,
            modifiers: egui::Modifiers::NONE,
        });
        let mut applied = None;
        let _ = ctx.run(raw, |c| {
            applied = show_paste_modal(c, &mut modal);
        });
        assert!(applied.is_none());
        assert!(!modal.open, "Escape should close paste modal");
    }
}
