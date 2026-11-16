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

//! Tool #10 (`osabt ecosound`) — GPU-accelerated cross-platform `eframe` / `egui` desktop studio
//! for STFT spectrogram visualization, interactive drag-to-box vocalization annotation,
//! frequency-band SNR call slicing, Raven Pro Selection Table export, and live `ACI`/`NDSI`/`H`/`ADI`
//! ecoacoustic biodiversity monitoring.

use std::fs;
use std::path::PathBuf;

use eframe::egui;
use egui::{
    Color32, ColorImage, Pos2, Rect, Sense, Stroke, TextureHandle, TextureOptions, Vec2,
};

use crate::ecosound_dsp::{
    compute_ecoacoustic_indices, compute_stft_spectrogram, demo_dawn_chorus_recording,
    demo_ultrasonic_bat_recording, parse_wav_bytes, AcousticSpectrogram, EcoacousticIndices,
    SpectrogramConfig, WavRecording, WindowFunction,
};
use crate::ecosound_slicer::{
    detect_calls_in_band, export_raven_selection_tsv, render_spectrogram_svg, slice_call_wav_clips,
    BioacousticCallEvent, CallSlicerConfig, SpectrogramPalette,
};

/// Interactive `egui` application state for `OSABT EcoSound Studio`.
pub struct EcosoundStudioApp {
    wav_path_input: String,
    export_dir_input: String,
    recording: WavRecording,
    spec_cfg: SpectrogramConfig,
    spectrogram: AcousticSpectrogram,
    indices: EcoacousticIndices,
    slicer_cfg: CallSlicerConfig,
    events: Vec<BioacousticCallEvent>,
    palette: SpectrogramPalette,
    texture: Option<TextureHandle>,
    texture_dirty: bool,
    drag_start_pos: Option<Pos2>,
    selected_event_idx: Option<usize>,
    status_bar: String,
}

impl EcosoundStudioApp {
    pub fn new(initial_wav: Option<PathBuf>) -> Self {
        let recording = demo_dawn_chorus_recording();
        let spec_cfg = SpectrogramConfig {
            fft_size: 1024,
            hop_size: 256,
            window: WindowFunction::Hann,
            min_freq_hz: 0.0,
            max_freq_hz: 10_000.0,
            floor_dbfs: -85.0,
        };
        let spectrogram =
            compute_stft_spectrogram(&recording, &spec_cfg).expect("valid demo STFT");
        let indices = compute_ecoacoustic_indices(
            &recording,
            &spectrogram,
            (1_000.0, 2_000.0),
            (2_000.0, 11_000.0),
        );
        let slicer_cfg = CallSlicerConfig::default();
        let events = detect_calls_in_band(&spectrogram, &slicer_cfg);

        let mut app = Self {
            wav_path_input: initial_wav
                .as_ref()
                .map(|p| p.to_string_lossy().to_string())
                .unwrap_or_default(),
            export_dir_input: String::from("ecosound_exports"),
            recording,
            spec_cfg,
            spectrogram,
            indices,
            slicer_cfg,
            events,
            palette: SpectrogramPalette::Magma,
            texture: None,
            texture_dirty: true,
            drag_start_pos: None,
            selected_event_idx: Some(0),
            status_bar: String::from(
                "Drag on spectrogram to draw custom Raven selection box | Right-click box to delete | Adjust band & click Auto-Detect.",
            ),
        };

        if let Some(path) = initial_wav {
            app.load_wav_file(&path);
        }
        app
    }

    fn load_wav_file(&mut self, path: &PathBuf) {
        match fs::read(path) {
            Ok(bytes) => {
                let title = path
                    .file_name()
                    .and_then(|n| n.to_str())
                    .unwrap_or("recording.wav");
                match parse_wav_bytes(title, &bytes) {
                    Ok(rec) => {
                        let nyq = rec.nyquist_hz();
                        self.spec_cfg.max_freq_hz = nyq;
                        self.slicer_cfg.low_freq_hz = (0.08 * nyq).max(200.0);
                        self.slicer_cfg.high_freq_hz = (0.85 * nyq).max(1_000.0);
                        self.recording = rec;
                        self.recompute_spectrogram_and_indices(true);
                        self.status_bar = format!(
                            "Loaded '{}' ({:.2} s @ {} Hz, {} selections).",
                            self.recording.title,
                            self.recording.duration_seconds(),
                            self.recording.sample_rate_hz,
                            self.events.len()
                        );
                    }
                    Err(e) => self.status_bar = format!("WAV parse error: {e}"),
                }
            }
            Err(e) => self.status_bar = format!("Failed reading '{}': {e}", path.display()),
        }
    }

    fn recompute_spectrogram_and_indices(&mut self, rerun_detector: bool) {
        if let Ok(spec) = compute_stft_spectrogram(&self.recording, &self.spec_cfg) {
            let nyq = self.recording.nyquist_hz();
            let bio_band = if nyq > 30_000.0 {
                (15_000.0, nyq)
            } else {
                (2_000.0, nyq.min(11_000.0))
            };
            let anthro_band = if nyq > 30_000.0 {
                (1_000.0, 15_000.0)
            } else {
                (1_000.0, 2_000.0)
            };
            self.indices =
                compute_ecoacoustic_indices(&self.recording, &spec, anthro_band, bio_band);
            if rerun_detector {
                self.events = detect_calls_in_band(&spec, &self.slicer_cfg);
                self.selected_event_idx = if self.events.is_empty() {
                    None
                } else {
                    Some(0)
                };
            }
            self.spectrogram = spec;
            self.texture_dirty = true;
        }
    }

    fn ensure_spectrogram_texture(&mut self, ctx: &egui::Context) {
        if !self.texture_dirty && self.texture.is_some() {
            return;
        }
        let w = self.spectrogram.num_frames.max(1);
        let h = self.spectrogram.num_bins.max(1);
        let dyn_db = (-self.spectrogram.floor_dbfs).max(20.0);
        let mut rgba = vec![0u8; w * h * 4];

        for fx in 0..w {
            for by in 0..h {
                let img_row = h - 1 - by;
                let db = self.spectrogram.dbfs_at(fx, by);
                let u = ((db - self.spectrogram.floor_dbfs) / dyn_db).clamp(0.0, 1.0);
                let [r, g, b] = self.palette.map_unit_rgb(u);
                let off = (img_row * w + fx) * 4;
                rgba[off] = r;
                rgba[off + 1] = g;
                rgba[off + 2] = b;
                rgba[off + 3] = 255;
            }
        }

        let cimg = ColorImage::from_rgba_unmultiplied([w, h], &rgba);
        self.texture = Some(ctx.load_texture(
            "osabt_ecosound_spectrogram",
            cimg,
            TextureOptions::LINEAR,
        ));
        self.texture_dirty = false;
    }

    fn measure_manual_box(&self, t0: f64, t1: f64, f0: f64, f1: f64) -> BioacousticCallEvent {
        let begin_t = t0.min(t1).max(0.0);
        let end_t = t0.max(t1).max(begin_t + 0.01);
        let lo_f = f0.min(f1).max(0.0);
        let hi_f = f0.max(f1).max(lo_f + 25.0);

        let mut peak_mag = 0.0_f64;
        let mut peak_f = 0.5 * (lo_f + hi_f);
        let mut peak_db = self.spectrogram.floor_dbfs;
        let mut w_f = 0.0_f64;
        let mut w_s = 0.0_f64;

        for (fr, &ft) in self.spectrogram.frame_times_s.iter().enumerate() {
            if ft < begin_t || ft > end_t {
                continue;
            }
            for (b, &bf) in self.spectrogram.bin_freqs_hz.iter().enumerate() {
                if bf < lo_f || bf > hi_f {
                    continue;
                }
                let m = self.spectrogram.mag_at(fr, b);
                let db = self.spectrogram.dbfs_at(fr, b);
                w_f += m * m * bf;
                w_s += m * m;
                if m > peak_mag {
                    peak_mag = m;
                    peak_f = bf;
                    peak_db = db;
                }
            }
        }

        BioacousticCallEvent {
            selection_id: self.events.len() + 1,
            begin_time_s: begin_t,
            end_time_s: end_t,
            duration_ms: (end_t - begin_t) * 1000.0,
            low_freq_hz: lo_f,
            high_freq_hz: hi_f,
            peak_freq_hz: peak_f,
            spectral_centroid_hz: if w_s > 1e-18 { w_f / w_s } else { peak_f },
            peak_dbfs: peak_db,
            snr_db: (peak_db - (self.spectrogram.floor_dbfs + 25.0)).max(0.0),
            fm_slope_hz_per_s: 0.0,
            annotation: self.slicer_cfg.species_label.clone(),
        }
    }

    fn renumber_selections(&mut self) {
        for (i, ev) in self.events.iter_mut().enumerate() {
            ev.selection_id = i + 1;
        }
    }
}

/// Launches the native cross-platform `OSABT EcoSound Studio` desktop GUI.
pub fn run_ecosound_gui(initial_wav: Option<PathBuf>) -> Result<(), String> {
    let options = eframe::NativeOptions {
        viewport: egui::ViewportBuilder::default()
            .with_inner_size([1280.0, 800.0])
            .with_title(
                "OSABT EcoSound — Field Bioacoustics Spectrogram, Ecoacoustic Indices & Call Slicer (GPL-3.0)",
            ),
        ..Default::default()
    };
    eframe::run_native(
        "osabt-ecosound-studio",
        options,
        Box::new(move |_cc| Ok(Box::new(EcosoundStudioApp::new(initial_wav)))),
    )
    .map_err(|e| format!("Failed starting EcoSound desktop window: {e}"))
}

impl eframe::App for EcosoundStudioApp {
    fn update(&mut self, ctx: &egui::Context, _frame: &mut eframe::Frame) {
        self.draw_panels(ctx);
    }
}

impl EcosoundStudioApp {
    /// Renders all EcoSound Studio panels onto `ctx`.
    pub fn draw_panels(&mut self, ctx: &egui::Context) {
        self.ensure_spectrogram_texture(ctx);

        egui::TopBottomPanel::top("ecosound_top_bar").show(ctx, |ui| {
            ui.horizontal_wrapped(|ui| {
                ui.heading("🦉 OSABT EcoSound");
                ui.separator();
                if ui.button("🌲 Demo: PNW Dawn Chorus (44.1 kHz)").clicked() {
                    self.recording = demo_dawn_chorus_recording();
                    self.spec_cfg.max_freq_hz = 10_000.0;
                    self.slicer_cfg = CallSlicerConfig {
                        species_label: String::from("PNW Riparian Avian/Anuran"),
                        ..CallSlicerConfig::default()
                    };
                    self.recompute_spectrogram_and_indices(true);
                }
                if ui.button("🦇 Demo: Myotis Ultrasonic (192 kHz)").clicked() {
                    self.recording = demo_ultrasonic_bat_recording();
                    self.spec_cfg.max_freq_hz = 96_000.0;
                    self.slicer_cfg.low_freq_hz = 32_000.0;
                    self.slicer_cfg.high_freq_hz = 90_000.0;
                    self.slicer_cfg.snr_threshold_db = 12.0;
                    self.slicer_cfg.min_duration_ms = 18.0;
                    self.slicer_cfg.merge_gap_ms = 12.0;
                    self.slicer_cfg.species_label = String::from("Myotis lucifugus FM pulse");
                    self.recompute_spectrogram_and_indices(true);
                }
                ui.separator();
                ui.label("WAV File:");
                ui.add(egui::TextEdit::singleline(&mut self.wav_path_input).desired_width(210.0));
                if ui.button("Open .wav").clicked() && !self.wav_path_input.trim().is_empty() {
                    let p = PathBuf::from(self.wav_path_input.trim());
                    self.load_wav_file(&p);
                }
            });
        });

        egui::TopBottomPanel::bottom("ecosound_status_footer").show(ctx, |ui| {
            ui.horizontal(|ui| {
                ui.label(&self.status_bar);
            });
        });

        egui::SidePanel::left("ecosound_controls")
            .default_width(315.0)
            .show(ctx, |ui| {
                egui::ScrollArea::vertical().show(ui, |ui| {
                    ui.heading("1. STFT Spectrogram");
                    ui.label(format!(
                        "{} ({:.2} s @ {} Hz)",
                        self.recording.title,
                        self.recording.duration_seconds(),
                        self.recording.sample_rate_hz
                    ));
                    if let Some(cmt) = &self.recording.comment_metadata {
                        ui.small(cmt);
                    }

                    let mut stft_changed = false;
                    ui.horizontal(|ui| {
                        ui.label("FFT Size:");
                        for sz in [512usize, 1024, 2048] {
                            if ui
                                .selectable_label(self.spec_cfg.fft_size == sz, format!("{sz}"))
                                .clicked()
                            {
                                self.spec_cfg.fft_size = sz;
                                self.spec_cfg.hop_size = sz / 4;
                                stft_changed = true;
                            }
                        }
                    });

                    egui::ComboBox::from_label("Window")
                        .selected_text(self.spec_cfg.window.label())
                        .show_ui(ui, |ui| {
                            for w in [
                                WindowFunction::Hann,
                                WindowFunction::Hamming,
                                WindowFunction::BlackmanHarris,
                            ] {
                                if ui
                                    .selectable_value(&mut self.spec_cfg.window, w, w.label())
                                    .changed()
                                {
                                    stft_changed = true;
                                }
                            }
                        });

                    egui::ComboBox::from_label("Palette")
                        .selected_text(self.palette.label())
                        .show_ui(ui, |ui| {
                            for p in [
                                SpectrogramPalette::Magma,
                                SpectrogramPalette::Viridis,
                                SpectrogramPalette::Inferno,
                                SpectrogramPalette::GrayscaleInv,
                            ] {
                                if ui
                                    .selectable_value(&mut self.palette, p, p.label())
                                    .changed()
                                {
                                    self.texture_dirty = true;
                                }
                            }
                        });

                    let nyq = self.recording.nyquist_hz();
                    if ui
                        .add(
                            egui::Slider::new(&mut self.spec_cfg.max_freq_hz, 1_000.0..=nyq)
                                .text("Max View (Hz)"),
                        )
                        .changed()
                    {
                        stft_changed = true;
                    }
                    if ui
                        .add(
                            egui::Slider::new(&mut self.spec_cfg.floor_dbfs, -110.0..=-45.0)
                                .text("Floor (dBFS)"),
                        )
                        .changed()
                    {
                        stft_changed = true;
                    }
                    if stft_changed {
                        self.recompute_spectrogram_and_indices(false);
                    }

                    ui.separator();
                    ui.heading("2. Frequency-Band Call Slicer");
                    ui.add(
                        egui::Slider::new(&mut self.slicer_cfg.low_freq_hz, 50.0..=nyq)
                            .text("Band Low (Hz)"),
                    );
                    ui.add(
                        egui::Slider::new(&mut self.slicer_cfg.high_freq_hz, 100.0..=nyq)
                            .text("Band High (Hz)"),
                    );
                    ui.add(
                        egui::Slider::new(&mut self.slicer_cfg.snr_threshold_db, 4.0..=30.0)
                            .text("SNR Trigger (dB)"),
                    );
                    ui.add(
                        egui::Slider::new(&mut self.slicer_cfg.min_duration_ms, 10.0..=600.0)
                            .text("Min Call (ms)"),
                    );
                    ui.checkbox(
                        &mut self.slicer_cfg.bandpass_clips,
                        "Butterworth Bandpass Sliced .wav Clips",
                    );
                    ui.horizontal(|ui| {
                        ui.label("Label:");
                        ui.text_edit_singleline(&mut self.slicer_cfg.species_label);
                    });

                    ui.horizontal(|ui| {
                        if ui.button("🔍 Auto-Detect Band Calls").clicked() {
                            self.events =
                                detect_calls_in_band(&self.spectrogram, &self.slicer_cfg);
                            self.selected_event_idx =
                                if self.events.is_empty() { None } else { Some(0) };
                            self.status_bar = format!(
                                "Detected {} vocalizations in [{:.0}, {:.0}] Hz.",
                                self.events.len(),
                                self.slicer_cfg.low_freq_hz,
                                self.slicer_cfg.high_freq_hz
                            );
                        }
                        if ui.button("Clear Boxes").clicked() {
                            self.events.clear();
                            self.selected_event_idx = None;
                        }
                    });

                    ui.separator();
                    ui.heading("3. Ecoacoustic Biodiversity Indices");
                    egui::Grid::new("eco_indices_grid")
                        .striped(true)
                        .show(ui, |ui| {
                            ui.label("ACI (Pieretti 2011):");
                            ui.strong(format!(
                                "{:.1} ({:.1}/min)",
                                self.indices.aci_total, self.indices.aci_per_minute
                            ));
                            ui.end_row();

                            ui.label("NDSI (Kasten 2012):");
                            ui.strong(format!("{:+.3}", self.indices.ndsi));
                            ui.end_row();

                            ui.label("Entropy H (Sueur 2008):");
                            ui.strong(format!(
                                "{:.3} (Ht={:.2}, Hf={:.2})",
                                self.indices.acoustic_entropy_h,
                                self.indices.temporal_entropy_ht,
                                self.indices.spectral_entropy_hf
                            ));
                            ui.end_row();

                            ui.label("ADI / AEI (Villanueva 2011):");
                            ui.strong(format!(
                                "{:.3} / {:.3}",
                                self.indices.adi_shannon, self.indices.aei_gini
                            ));
                            ui.end_row();
                        });

                    ui.separator();
                    ui.heading("4. Export Clips, Raven TSV & SVG");
                    ui.horizontal(|ui| {
                        ui.label("Out Dir:");
                        ui.text_edit_singleline(&mut self.export_dir_input);
                    });
                    if ui
                        .button("💾 Export Sliced .wav Clips + Raven TSV + SVG")
                        .clicked()
                    {
                        let out_dir = PathBuf::from(self.export_dir_input.trim());
                        if let Err(e) = fs::create_dir_all(&out_dir) {
                            self.status_bar = format!("Failed creating output dir: {e}");
                        } else {
                            let raven =
                                export_raven_selection_tsv(&self.events, &self.recording.title);
                            let raven_path = out_dir.join("selections.Raven.1.selections.txt");
                            let _ = fs::write(&raven_path, raven);

                            let svg = render_spectrogram_svg(
                                &self.spectrogram,
                                &self.events,
                                Some((
                                    self.slicer_cfg.low_freq_hz,
                                    self.slicer_cfg.high_freq_hz,
                                )),
                                self.palette,
                                &self.recording.title,
                            );
                            let _ = fs::write(out_dir.join("spectrogram_annotated.svg"), svg);

                            let clips = slice_call_wav_clips(
                                &self.recording,
                                &self.events,
                                &self.slicer_cfg,
                                "call",
                            );
                            for (fname, bytes) in &clips {
                                let _ = fs::write(out_dir.join(fname), bytes);
                            }
                            self.status_bar = format!(
                                "Exported {} .wav call clips, Raven selection table, and annotated .svg to '{}'.",
                                clips.len(),
                                out_dir.display()
                            );
                        }
                    }
                });
            });

        egui::CentralPanel::default().show(ctx, |ui| {
            let total_h = ui.available_height();
            let canvas_h = (total_h * 0.66).max(240.0);

            let (resp, painter) = ui.allocate_painter(
                Vec2::new(ui.available_width(), canvas_h),
                Sense::click_and_drag(),
            );
            let rect = resp.rect;
            painter.rect_filled(rect, 4.0, Color32::from_rgb(12, 18, 30));

            if let Some(tex) = &self.texture {
                painter.image(
                    tex.id(),
                    rect,
                    Rect::from_min_max(Pos2::ZERO, Pos2::new(1.0, 1.0)),
                    Color32::WHITE,
                );
            }

            let dur = self.recording.duration_seconds().max(1e-6);
            let f_min = *self.spectrogram.bin_freqs_hz.first().unwrap_or(&0.0);
            let f_max = *self
                .spectrogram
                .bin_freqs_hz
                .last()
                .unwrap_or(&self.recording.nyquist_hz());
            let f_span = (f_max - f_min).max(1.0);

            let map_x = |t: f64| rect.left() + ((t / dur).clamp(0.0, 1.0) as f32) * rect.width();
            let map_y = |f: f64| {
                rect.bottom() - (((f - f_min) / f_span).clamp(0.0, 1.0) as f32) * rect.height()
            };
            let unmap_t = |px: f32| (((px - rect.left()) / rect.width()).clamp(0.0, 1.0) as f64) * dur;
            let unmap_f = |py: f32| {
                f_min + (((rect.bottom() - py) / rect.height()).clamp(0.0, 1.0) as f64) * f_span
            };

            for f_band in [self.slicer_cfg.low_freq_hz, self.slicer_cfg.high_freq_hz] {
                if f_band >= f_min && f_band <= f_max {
                    let y = map_y(f_band);
                    painter.line_segment(
                        [Pos2::new(rect.left(), y), Pos2::new(rect.right(), y)],
                        Stroke::new(1.3_f32, Color32::from_rgba_unmultiplied(56, 189, 248, 200)),
                    );
                }
            }

            let mut right_clicked_delete: Option<usize> = None;
            for (idx, ev) in self.events.iter().enumerate() {
                let b_rect = Rect::from_min_max(
                    Pos2::new(map_x(ev.begin_time_s), map_y(ev.high_freq_hz)),
                    Pos2::new(map_x(ev.end_time_s), map_y(ev.low_freq_hz)),
                );
                let is_sel = self.selected_event_idx == Some(idx);
                let stroke_col = if is_sel {
                    Color32::from_rgb(250, 204, 21)
                } else {
                    Color32::from_rgb(34, 211, 238)
                };
                painter.rect_stroke(b_rect, 2.0, Stroke::new(if is_sel { 2.4 } else { 1.6_f32 }, stroke_col));
                painter.text(
                    Pos2::new(b_rect.left() + 3.0, (b_rect.top() - 13.0).max(rect.top() + 2.0)),
                    egui::Align2::LEFT_TOP,
                    format!("#{} {:.1}kHz", ev.selection_id, ev.peak_freq_hz / 1000.0),
                    egui::FontId::monospace(11.0),
                    stroke_col,
                );

                if resp.clicked() {
                    if let Some(pos) = resp.interact_pointer_pos() {
                        if b_rect.contains(pos) {
                            self.selected_event_idx = Some(idx);
                        }
                    }
                }
                if resp.secondary_clicked() {
                    if let Some(pos) = resp.interact_pointer_pos() {
                        if b_rect.contains(pos) {
                            right_clicked_delete = Some(idx);
                        }
                    }
                }
            }

            if let Some(del_i) = right_clicked_delete {
                self.events.remove(del_i);
                self.renumber_selections();
                self.selected_event_idx = None;
            }

            if resp.drag_started() {
                self.drag_start_pos = resp.interact_pointer_pos();
            }
            if let (Some(p0), Some(p1)) = (self.drag_start_pos, resp.interact_pointer_pos()) {
                if resp.dragged() {
                    let preview = Rect::from_two_pos(p0, p1);
                    painter.rect_stroke(
                        preview,
                        2.0,
                        Stroke::new(2.0_f32, Color32::from_rgb(74, 222, 128)),
                    );
                }
            }
            if resp.drag_stopped() {
                if let (Some(p0), Some(p1)) = (self.drag_start_pos.take(), resp.interact_pointer_pos())
                {
                    if (p1.x - p0.x).abs() > 6.0 && (p1.y - p0.y).abs() > 6.0 {
                        let new_ev =
                            self.measure_manual_box(unmap_t(p0.x), unmap_t(p1.x), unmap_f(p0.y), unmap_f(p1.y));
                        self.events.push(new_ev);
                        self.selected_event_idx = Some(self.events.len() - 1);
                    }
                }
            }

            if let Some(hover) = resp.hover_pos() {
                let ht = unmap_t(hover.x);
                let hf = unmap_f(hover.y);
                self.status_bar = format!(
                    "Cursor: t = {ht:.3} s | f = {hf:.0} Hz | Selections: {} (Drag to add box, right-click box to remove)",
                    self.events.len()
                );
            }

            ui.separator();
            ui.strong(format!(
                "Raven Pro Selection Table Preview ({} vocalizations)",
                self.events.len()
            ));
            egui::ScrollArea::vertical().show(ui, |ui| {
                egui::Grid::new("raven_selections_table")
                    .striped(true)
                    .min_col_width(68.0)
                    .show(ui, |ui| {
                        ui.strong("Sel #");
                        ui.strong("Begin (s)");
                        ui.strong("End (s)");
                        ui.strong("Dur (ms)");
                        ui.strong("Low–High (Hz)");
                        ui.strong("Peak (Hz)");
                        ui.strong("SNR (dB)");
                        ui.strong("FM Slope");
                        ui.strong("Annotation");
                        ui.end_row();

                        for (idx, ev) in self.events.iter_mut().enumerate() {
                            if ui
                                .selectable_label(
                                    self.selected_event_idx == Some(idx),
                                    format!("#{}", ev.selection_id),
                                )
                                .clicked()
                            {
                                self.selected_event_idx = Some(idx);
                            }
                            ui.label(format!("{:.3}", ev.begin_time_s));
                            ui.label(format!("{:.3}", ev.end_time_s));
                            ui.label(format!("{:.0}", ev.duration_ms));
                            ui.label(format!("{:.0}–{:.0}", ev.low_freq_hz, ev.high_freq_hz));
                            ui.label(format!("{:.0}", ev.peak_freq_hz));
                            ui.label(format!("{:.1}", ev.snr_db));
                            ui.label(format!("{:+.0} Hz/s", ev.fm_slope_hz_per_s));
                            ui.text_edit_singleline(&mut ev.annotation);
                            ui.end_row();
                        }
                    });
            });
        });
    }
}
