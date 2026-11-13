// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #10 (`osabt ecosound`) — Adaptive SNR frequency-band bioacoustic call detector,
//! individual `.wav` call-clip slicer (with optional Butterworth bandpass & cosine micro-fades),
//! Cornell Lab of Ornithology **Raven Pro 1.6** Selection Table `.tsv` exporter, and
//! publication-grade `.svg` spectrogram + bounding-box renderer.

use std::f64::consts::PI;
use std::fmt::Write as _;

use crate::ecosound_dsp::{
    apply_bandpass_biquad, encode_wav_pcm16, AcousticSpectrogram, WavRecording,
};

/// Color palettes for spectrogram rendering in `.svg` and `egui` viewports.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum SpectrogramPalette {
    Viridis,
    Magma,
    Inferno,
    GrayscaleInv,
}

impl SpectrogramPalette {
    pub fn label(self) -> &'static str {
        match self {
            Self::Viridis => "Viridis",
            Self::Magma => "Magma",
            Self::Inferno => "Inferno",
            Self::GrayscaleInv => "Inverted Grayscale (Print)",
        }
    }

    /// Maps normalized intensity `t in [0.0, 1.0]` to sRGB `[R, G, B]`.
    pub fn map_unit_rgb(self, t: f64) -> [u8; 3] {
        let u = t.clamp(0.0, 1.0);
        let stops: &[[f64; 3]] = match self {
            Self::Viridis => &[
                [68.0, 1.0, 84.0],
                [59.0, 82.0, 139.0],
                [33.0, 145.0, 140.0],
                [94.0, 201.0, 98.0],
                [253.0, 231.0, 37.0],
            ],
            Self::Magma => &[
                [0.0, 0.0, 4.0],
                [81.0, 18.0, 124.0],
                [183.0, 55.0, 121.0],
                [252.0, 137.0, 97.0],
                [252.0, 253.0, 191.0],
            ],
            Self::Inferno => &[
                [0.0, 0.0, 4.0],
                [87.0, 16.0, 110.0],
                [188.0, 55.0, 84.0],
                [249.0, 142.0, 9.0],
                [252.0, 255.0, 164.0],
            ],
            Self::GrayscaleInv => &[
                [250.0, 250.0, 252.0],
                [195.0, 200.0, 208.0],
                [130.0, 138.0, 150.0],
                [65.0, 72.0, 84.0],
                [12.0, 16.0, 24.0],
            ],
        };
        let scaled = u * ((stops.len() - 1) as f64);
        let idx = (scaled.floor() as usize).min(stops.len() - 2);
        let frac = scaled - (idx as f64);
        let a = stops[idx];
        let b = stops[idx + 1];
        [
            (a[0] + frac * (b[0] - a[0])).round() as u8,
            (a[1] + frac * (b[1] - a[1])).round() as u8,
            (a[2] + frac * (b[2] - a[2])).round() as u8,
        ]
    }
}

/// Configuration for frequency-band bioacoustic call segmentation & `.wav` slicing.
#[derive(Debug, Clone)]
pub struct CallSlicerConfig {
    /// Lower edge of target frequency band (`Hz`).
    pub low_freq_hz: f64,
    /// Upper edge of target frequency band (`Hz`).
    pub high_freq_hz: f64,
    /// Detection threshold (`dB`) above the target band's median background noise floor.
    pub snr_threshold_db: f64,
    /// Minimum contiguous call duration (`ms`) required to retain an event.
    pub min_duration_ms: f64,
    /// Maximum inter-syllable gap (`ms`) bridged into a single phrase/selection.
    pub merge_gap_ms: f64,
    /// Context padding (`ms`) prepended before each sliced `.wav` clip.
    pub pad_before_ms: f64,
    /// Context padding (`ms`) appended after each sliced `.wav` clip.
    pub pad_after_ms: f64,
    /// Apply Butterworth bandpass `[low_freq_hz, high_freq_hz]` when exporting sliced `.wav` clips.
    pub bandpass_clips: bool,
    /// Default taxon or call-type label populated into Raven `Annotation` column.
    pub species_label: String,
}

impl Default for CallSlicerConfig {
    fn default() -> Self {
        Self {
            low_freq_hz: 1_800.0,
            high_freq_hz: 8_000.0,
            snr_threshold_db: 10.0,
            min_duration_ms: 120.0,
            merge_gap_ms: 145.0,
            pad_before_ms: 40.0,
            pad_after_ms: 60.0,
            bandpass_clips: false,
            species_label: String::from("Avian/Anuran Vocalization"),
        }
    }
}

/// Detected or manually curated bioacoustic vocalization bounding box + spectral contour summary.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct BioacousticCallEvent {
    pub selection_id: usize,
    pub begin_time_s: f64,
    pub end_time_s: f64,
    pub duration_ms: f64,
    pub low_freq_hz: f64,
    pub high_freq_hz: f64,
    pub peak_freq_hz: f64,
    pub spectral_centroid_hz: f64,
    pub peak_dbfs: f64,
    pub snr_db: f64,
    /// Linear regression slope (`Hz/s`) of framewise peak frequency across the call (`>0` upsweep, `<0` downsweep).
    pub fm_slope_hz_per_s: f64,
    pub annotation: String,
}

/// Detects vocalizations inside `[cfg.low_freq_hz, cfg.high_freq_hz]` using an adaptive
/// median-floor SNR detector, gap merging, `-14 dB` spectral bounding-box tightening, and FM slope estimation.
pub fn detect_calls_in_band(
    spec: &AcousticSpectrogram,
    cfg: &CallSlicerConfig,
) -> Vec<BioacousticCallEvent> {
    if spec.num_frames == 0 || spec.num_bins == 0 {
        return Vec::new();
    }
    let f_lo = cfg.low_freq_hz.min(cfg.high_freq_hz);
    let f_hi = cfg.low_freq_hz.max(cfg.high_freq_hz);

    let band_bins: Vec<usize> = spec
        .bin_freqs_hz
        .iter()
        .enumerate()
        .filter(|&(_, &f)| f >= f_lo && f <= f_hi)
        .map(|(i, _)| i)
        .collect();
    if band_bins.is_empty() {
        return Vec::new();
    }

    let mut frame_band_db = Vec::with_capacity(spec.num_frames);
    for t in 0..spec.num_frames {
        let mut energy = 0.0_f64;
        for &b in &band_bins {
            let m = spec.mag_at(t, b);
            energy += m * m;
        }
        let rms_mag = (energy / (band_bins.len() as f64)).sqrt();
        let db = if rms_mag > 1e-12 {
            20.0 * rms_mag.log10()
        } else {
            spec.floor_dbfs
        };
        frame_band_db.push(db);
    }

    let mut sorted_db = frame_band_db.clone();
    sorted_db.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let noise_floor_db = sorted_db[sorted_db.len() / 2];
    let trigger_db = noise_floor_db + cfg.snr_threshold_db;

    let dt_s = (spec.hop_size as f64) / (spec.sample_rate_hz as f64);
    let half_win_s = 0.5 * (spec.fft_size as f64) / (spec.sample_rate_hz as f64);
    let max_gap_frames = ((cfg.merge_gap_ms / 1000.0) / dt_s).round() as usize;

    let mut spans: Vec<(usize, usize)> = Vec::new();
    let mut t = 0usize;
    while t < spec.num_frames {
        if frame_band_db[t] >= trigger_db {
            let start = t;
            let mut end = t;
            while end + 1 < spec.num_frames && frame_band_db[end + 1] >= trigger_db - 2.5 {
                end += 1;
            }
            if let Some(prev) = spans.last_mut() {
                if start <= prev.1 + max_gap_frames + 1 {
                    prev.1 = end;
                } else {
                    spans.push((start, end));
                }
            } else {
                spans.push((start, end));
            }
            t = end + 1;
        } else {
            t += 1;
        }
    }

    let mut events = Vec::new();
    for (s_frame, e_frame) in spans {
        let t_begin = (spec.frame_times_s[s_frame] - half_win_s).max(0.0);
        let t_end = spec.frame_times_s[e_frame] + half_win_s;
        let dur_ms = (t_end - t_begin) * 1000.0;
        if dur_ms < cfg.min_duration_ms {
            continue;
        }

        let mut peak_mag = 0.0_f64;
        let mut peak_bin = band_bins[0];
        let mut peak_dbfs = spec.floor_dbfs;
        let mut weighted_freq = 0.0_f64;
        let mut weight_sum = 0.0_f64;
        let mut bin_max_db = vec![spec.floor_dbfs; band_bins.len()];

        let mut ridge_t = Vec::with_capacity(e_frame - s_frame + 1);
        let mut ridge_f = Vec::with_capacity(e_frame - s_frame + 1);

        for fr in s_frame..=e_frame {
            let mut fr_best_mag = 0.0_f64;
            let mut fr_best_freq = spec.bin_freqs_hz[band_bins[0]];
            for (k_idx, &b) in band_bins.iter().enumerate() {
                let m = spec.mag_at(fr, b);
                let db = spec.dbfs_at(fr, b);
                let p = m * m;
                weighted_freq += p * spec.bin_freqs_hz[b];
                weight_sum += p;
                if db > bin_max_db[k_idx] {
                    bin_max_db[k_idx] = db;
                }
                if m > peak_mag {
                    peak_mag = m;
                    peak_bin = b;
                    peak_dbfs = db;
                }
                if m > fr_best_mag {
                    fr_best_mag = m;
                    fr_best_freq = spec.bin_freqs_hz[b];
                }
            }
            if frame_band_db[fr] >= noise_floor_db + 0.5 * cfg.snr_threshold_db {
                ridge_t.push(spec.frame_times_s[fr]);
                ridge_f.push(fr_best_freq);
            }
        }

        let cutoff_db = peak_dbfs - 14.0;
        let mut tight_lo = spec.bin_freqs_hz[peak_bin];
        let mut tight_hi = spec.bin_freqs_hz[peak_bin];
        for (k_idx, &b) in band_bins.iter().enumerate() {
            if bin_max_db[k_idx] >= cutoff_db {
                let f = spec.bin_freqs_hz[b];
                tight_lo = tight_lo.min(f);
                tight_hi = tight_hi.max(f);
            }
        }
        let df = if spec.bin_freqs_hz.len() >= 2 {
            spec.bin_freqs_hz[1] - spec.bin_freqs_hz[0]
        } else {
            50.0
        };
        let box_lo = (tight_lo - df).max(f_lo);
        let box_hi = (tight_hi + df).min(f_hi).max(box_lo + df);

        let centroid_hz = if weight_sum > 1e-18 {
            weighted_freq / weight_sum
        } else {
            spec.bin_freqs_hz[peak_bin]
        };

        let fm_slope = linear_slope(&ridge_t, &ridge_f);
        let max_band_db = (s_frame..=e_frame)
            .map(|fr| frame_band_db[fr])
            .fold(spec.floor_dbfs, f64::max);
        let snr_db = (max_band_db - noise_floor_db).max(0.0);

        events.push(BioacousticCallEvent {
            selection_id: events.len() + 1,
            begin_time_s: t_begin,
            end_time_s: t_end,
            duration_ms: dur_ms,
            low_freq_hz: box_lo,
            high_freq_hz: box_hi,
            peak_freq_hz: spec.bin_freqs_hz[peak_bin],
            spectral_centroid_hz: centroid_hz,
            peak_dbfs,
            snr_db,
            fm_slope_hz_per_s: fm_slope,
            annotation: cfg.species_label.clone(),
        });
    }

    events
}

fn linear_slope(xs: &[f64], ys: &[f64]) -> f64 {
    if xs.len() < 2 || xs.len() != ys.len() {
        return 0.0;
    }
    let n = xs.len() as f64;
    let mx = xs.iter().sum::<f64>() / n;
    let my = ys.iter().sum::<f64>() / n;
    let mut num = 0.0_f64;
    let mut den = 0.0_f64;
    for (&x, &y) in xs.iter().zip(ys) {
        num += (x - mx) * (y - my);
        den += (x - mx) * (x - mx);
    }
    if den > 1e-12 {
        num / den
    } else {
        0.0
    }
}

/// Slices individual `.wav` audio clips for each [`BioacousticCallEvent`] with configurable padding,
/// optional Butterworth frequency-band isolation, and 3 ms cosine micro-fades.
/// Returns `(clip_filename, wav_bytes)` pairs.
pub fn slice_call_wav_clips(
    rec: &WavRecording,
    events: &[BioacousticCallEvent],
    cfg: &CallSlicerConfig,
    file_prefix: &str,
) -> Vec<(String, Vec<u8>)> {
    let sr = rec.sample_rate_hz as f64;
    let total = rec.samples_mono.len();
    let fade_samples = ((0.003 * sr).round() as usize).max(1);
    let mut out = Vec::with_capacity(events.len());

    let source_audio: Vec<f32> = if cfg.bandpass_clips {
        apply_bandpass_biquad(
            &rec.samples_mono,
            rec.sample_rate_hz,
            cfg.low_freq_hz.min(cfg.high_freq_hz),
            cfg.low_freq_hz.max(cfg.high_freq_hz),
        )
    } else {
        rec.samples_mono.clone()
    };

    for ev in events {
        let t0 = (ev.begin_time_s - cfg.pad_before_ms / 1000.0).max(0.0);
        let t1 = ev.end_time_s + cfg.pad_after_ms / 1000.0;
        let i0 = ((t0 * sr).floor() as usize).min(total);
        let i1 = ((t1 * sr).ceil() as usize).min(total);
        if i1 <= i0 + 8 {
            continue;
        }

        let mut clip = source_audio[i0..i1].to_vec();
        let f_len = fade_samples.min(clip.len() / 2);
        for k in 0..f_len {
            let gain = 0.5 * (1.0 - (PI * (k as f64) / (f_len as f64)).cos()) as f32;
            clip[k] *= gain;
            let tail = clip.len() - 1 - k;
            clip[tail] *= gain;
        }

        let fname = format!(
            "{}_sel{:03}_{:.0}ms_{:.0}Hz.wav",
            sanitize_stem(file_prefix),
            ev.selection_id,
            ev.begin_time_s * 1000.0,
            ev.peak_freq_hz
        );
        let provenance = format!(
            "OSABT ecosound slice #{:03} from {} [{:.3}-{:.3}s, {:.0}-{:.0} Hz, Peak={:.0} Hz, SNR={:.1} dB] {}",
            ev.selection_id,
            rec.title,
            ev.begin_time_s,
            ev.end_time_s,
            ev.low_freq_hz,
            ev.high_freq_hz,
            ev.peak_freq_hz,
            ev.snr_db,
            ev.annotation
        );
        let wav_bytes = encode_wav_pcm16(rec.sample_rate_hz, &clip, Some(&provenance));
        out.push((fname, wav_bytes));
    }

    out
}

/// Exports detected or curated [`BioacousticCallEvent`] boxes as a **Cornell Lab of Ornithology
/// Raven Pro 1.6** Selection Table (`.Selection.1.selections.txt` / `.tsv`), preserving Raven's
/// mandatory 7 leftmost columns (`Selection`, `View`, `Channel`, `Begin Time (s)`, `End Time (s)`,
/// `Low Freq (Hz)`, `High Freq (Hz)`).
pub fn export_raven_selection_tsv(
    events: &[BioacousticCallEvent],
    begin_file_name: &str,
) -> String {
    let mut tsv = String::with_capacity(256 + events.len() * 140);
    tsv.push_str("Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)\tLow Freq (Hz)\tHigh Freq (Hz)\tBegin File\tPeak Freq (Hz)\tDelta Time (s)\tSNR (dB)\tFM Slope (Hz/s)\tAnnotation\n");
    for ev in events {
        let clean_note = ev.annotation.replace(['\t', '\r', '\n'], " ");
        let clean_file = begin_file_name.replace(['\t', '\r', '\n'], "_");
        let _ = writeln!(
            tsv,
            "{}\tSpectrogram 1\t1\t{:.4}\t{:.4}\t{:.1}\t{:.1}\t{}\t{:.1}\t{:.4}\t{:.2}\t{:.1}\t{}",
            ev.selection_id,
            ev.begin_time_s,
            ev.end_time_s,
            ev.low_freq_hz,
            ev.high_freq_hz,
            clean_file,
            ev.peak_freq_hz,
            ev.end_time_s - ev.begin_time_s,
            ev.snr_db,
            ev.fm_slope_hz_per_s,
            clean_note
        );
    }
    tsv
}

/// Renders a publication-ready `.svg` spectrogram with optional target frequency-band cursors
/// and annotated call-selection bounding boxes.
pub fn render_spectrogram_svg(
    spec: &AcousticSpectrogram,
    events: &[BioacousticCallEvent],
    active_band_hz: Option<(f64, f64)>,
    palette: SpectrogramPalette,
    title: &str,
) -> String {
    let width = 920.0_f64;
    let height = 460.0_f64;
    let left = 74.0_f64;
    let top = 52.0_f64;
    let pw = 760.0_f64;
    let ph = 340.0_f64;

    let nx = spec.num_frames.clamp(1, 220);
    let ny = spec.num_bins.clamp(1, 110);
    let cw = pw / (nx as f64);
    let ch = ph / (ny as f64);

    let f_min = *spec.bin_freqs_hz.first().unwrap_or(&0.0);
    let f_max = *spec.bin_freqs_hz.last().unwrap_or(&22_050.0);
    let f_span = (f_max - f_min).max(1.0);

    let half_w = 0.5 * (spec.fft_size as f64) / (spec.sample_rate_hz as f64);
    let t_min = 0.0_f64;
    let t_max = spec
        .frame_times_s
        .last()
        .map(|&t| t + half_w)
        .unwrap_or(1.0)
        .max(1e-6);

    let dyn_db = (-spec.floor_dbfs).max(20.0);

    let mut svg = String::with_capacity(65_536);
    let _ = writeln!(
        svg,
        "<svg xmlns=\"http://www.w3.org/2000/svg\" viewBox=\"0 0 {width:.0} {height:.0}\" width=\"{width:.0}\" height=\"{height:.0}\">"
    );
    let _ = writeln!(
        svg,
        "  <rect x=\"0\" y=\"0\" width=\"{width:.0}\" height=\"{height:.0}\" fill=\"#0f172a\"/>"
    );
    let _ = writeln!(
        svg,
        "  <text x=\"{left:.0}\" y=\"30\" font-family=\"sans-serif\" font-size=\"15\" font-weight=\"700\" fill=\"#f8fafc\">{} — STFT Bioacoustic Spectrogram ({} selections)</text>",
        xml_escape(title),
        events.len()
    );

    for gx in 0..nx {
        let f0 = (gx * spec.num_frames) / nx;
        let f1 = (((gx + 1) * spec.num_frames) / nx).max(f0 + 1).min(spec.num_frames);
        let x = left + (gx as f64) * cw;
        for gy in 0..ny {
            let b0 = (gy * spec.num_bins) / ny;
            let b1 = (((gy + 1) * spec.num_bins) / ny).max(b0 + 1).min(spec.num_bins);
            let mut best_db = spec.floor_dbfs;
            for fr in f0..f1 {
                for b in b0..b1 {
                    let db = spec.dbfs_at(fr, b);
                    if db > best_db {
                        best_db = db;
                    }
                }
            }
            let unit = ((best_db - spec.floor_dbfs) / dyn_db).clamp(0.0, 1.0);
            let [r, g, b_col] = palette.map_unit_rgb(unit);
            let y = top + ph - ((gy + 1) as f64) * ch;
            let _ = writeln!(
                svg,
                "  <rect x=\"{x:.2}\" y=\"{y:.2}\" width=\"{:.2}\" height=\"{:.2}\" fill=\"#{r:02x}{g:02x}{b_col:02x}\"/>",
                cw + 0.35,
                ch + 0.35
            );
        }
    }

    let map_t = |t: f64| left + ((t - t_min) / (t_max - t_min)).clamp(0.0, 1.0) * pw;
    let map_f = |f: f64| top + ph - ((f - f_min) / f_span).clamp(0.0, 1.0) * ph;

    if let Some((b_lo, b_hi)) = active_band_hz {
        for f_line in [b_lo, b_hi] {
            if f_line >= f_min && f_line <= f_max {
                let py = map_f(f_line);
                let _ = writeln!(
                    svg,
                    "  <line x1=\"{left:.1}\" y1=\"{py:.1}\" x2=\"{:.1}\" y2=\"{py:.1}\" stroke=\"#38bdf8\" stroke-width=\"1.4\" stroke-dasharray=\"6,4\"/>",
                    left + pw
                );
            }
        }
    }

    for ev in events {
        let x0 = map_t(ev.begin_time_s);
        let x1 = map_t(ev.end_time_s).max(x0 + 2.0);
        let y_hi = map_f(ev.high_freq_hz);
        let y_lo = map_f(ev.low_freq_hz).max(y_hi + 3.0);
        let _ = writeln!(
            svg,
            "  <rect x=\"{x0:.1}\" y=\"{y_hi:.1}\" width=\"{:.1}\" height=\"{:.1}\" fill=\"#22d3ee\" fill-opacity=\"0.12\" stroke=\"#22d3ee\" stroke-width=\"1.8\"/>",
            x1 - x0,
            y_lo - y_hi
        );
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" font-family=\"monospace\" font-size=\"10.5\" font-weight=\"700\" fill=\"#67e8f9\">#{} {:.1}kHz</text>",
            x0 + 2.0,
            (y_hi - 4.0).max(top + 11.0),
            ev.selection_id,
            ev.peak_freq_hz / 1000.0
        );
    }

    let _ = writeln!(
        svg,
        "  <rect x=\"{left:.1}\" y=\"{top:.1}\" width=\"{pw:.1}\" height=\"{ph:.1}\" fill=\"none\" stroke=\"#94a3b8\" stroke-width=\"1.2\"/>"
    );
    for tick in 0..=5 {
        let u = (tick as f64) / 5.0;
        let tx = left + u * pw;
        let t_val = t_min + u * (t_max - t_min);
        let _ = writeln!(
            svg,
            "  <text x=\"{tx:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"monospace\" font-size=\"11\" fill=\"#cbd5e1\">{t_val:.2}s</text>",
            top + ph + 20.0
        );

        let fy = top + ph - u * ph;
        let f_val_khz = (f_min + u * f_span) / 1000.0;
        let _ = writeln!(
            svg,
            "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"end\" font-family=\"monospace\" font-size=\"11\" fill=\"#cbd5e1\">{f_val_khz:.1}k</text>",
            left - 8.0,
            fy + 4.0
        );
    }
    let _ = writeln!(
        svg,
        "  <text x=\"{:.1}\" y=\"{:.1}\" text-anchor=\"middle\" font-family=\"sans-serif\" font-size=\"12\" fill=\"#e2e8f0\">Time (s)  |  Palette: {}  |  Floor: {:.0} dBFS</text>",
        left + 0.5 * pw,
        height - 16.0,
        palette.label(),
        spec.floor_dbfs
    );

    svg.push_str("</svg>\n");
    svg
}

fn sanitize_stem(raw: &str) -> String {
    let s: String = raw
        .chars()
        .map(|c| {
            if c.is_ascii_alphanumeric() || c == '-' || c == '_' {
                c
            } else {
                '_'
            }
        })
        .collect();
    if s.is_empty() {
        String::from("clip")
    } else {
        s
    }
}

fn xml_escape(s: &str) -> String {
    s.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
        .replace('\'', "&apos;")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::ecosound_dsp::{
        compute_stft_spectrogram, demo_dawn_chorus_recording, parse_wav_bytes, SpectrogramConfig,
    };

    #[test]
    fn detects_all_five_phrases_in_demo_dawn_chorus_and_exports_raven_and_wav_clips() {
        let rec = demo_dawn_chorus_recording();
        let spec = compute_stft_spectrogram(&rec, &SpectrogramConfig::default()).unwrap();
        let cfg = CallSlicerConfig::default();
        let events = detect_calls_in_band(&spec, &cfg);
        assert_eq!(events.len(), 5);

        let raven = export_raven_selection_tsv(&events, &rec.title);
        assert!(raven.starts_with(
            "Selection\tView\tChannel\tBegin Time (s)\tEnd Time (s)\tLow Freq (Hz)\tHigh Freq (Hz)"
        ));

        let clips = slice_call_wav_clips(&rec, &events, &cfg, "dawn");
        assert_eq!(clips.len(), 5);
        let first_clip = parse_wav_bytes(&clips[0].0, &clips[0].1).unwrap();
        assert!(first_clip.duration_seconds() > 0.45);
    }
}
