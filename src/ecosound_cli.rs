// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Command-line subcommands (`clap` v4) for Tool #10 (`osabt ecosound`):
//! - `osabt ecosound gui` — Launches the GPU-accelerated `eframe`/`egui` spectrogram & call-slicer studio
//! - `osabt ecosound indices` — Computes `ACI` (Pieretti 2011), `NDSI` (Kasten 2012), `H` (Sueur 2008), & `ADI`/`AEI` (Villanueva-Rivera 2011)
//! - `osabt ecosound slice` — Segments vocalizations in `[low_hz, high_hz]`, exports sliced `.wav` clips, Raven Pro `.tsv`, and `.svg` spectrogram
//! - `osabt ecosound demo` — Writes synthetic 44.1 kHz PNW dawn-chorus or 192 kHz ultrasonic bat `.wav` recordings to disk

use std::fs;
use std::path::PathBuf;

use clap::{Args, Subcommand, ValueEnum};

use crate::ecosound_dsp::{
    compute_ecoacoustic_indices, compute_stft_spectrogram, demo_dawn_chorus_recording,
    demo_ultrasonic_bat_recording, encode_wav_pcm16, parse_wav_bytes, SpectrogramConfig,
    WavRecording, WindowFunction,
};
use crate::ecosound_slicer::{
    detect_calls_in_band, export_raven_selection_tsv, render_spectrogram_svg, slice_call_wav_clips,
    CallSlicerConfig, SpectrogramPalette,
};

#[derive(Debug, Clone, Copy, ValueEnum)]
pub enum CliPalette {
    Magma,
    Viridis,
    Inferno,
    Gray,
}

impl From<CliPalette> for SpectrogramPalette {
    fn from(v: CliPalette) -> Self {
        match v {
            CliPalette::Magma => Self::Magma,
            CliPalette::Viridis => Self::Viridis,
            CliPalette::Inferno => Self::Inferno,
            CliPalette::Gray => Self::GrayscaleInv,
        }
    }
}

/// Top-level CLI arguments for `osabt ecosound`.
#[derive(Debug, Args)]
pub struct EcosoundArgs {
    #[command(subcommand)]
    pub command: Option<EcosoundSubcommand>,
}

#[derive(Debug, Subcommand)]
pub enum EcosoundSubcommand {
    /// Launch the interactive `eframe`/`egui` bioacoustic spectrogram & call-slicer desktop GUI.
    Gui {
        /// Optional input `.wav` file to open on startup (defaults to PNW dawn-chorus demo).
        #[arg(long)]
        wav: Option<PathBuf>,
    },

    /// Compute ecoacoustic biodiversity indices (`ACI`, `NDSI`, `H = Ht * Hf`, `ADI`, `AEI`) on a `.wav` file or built-in soundscape.
    Indices {
        /// Input `.wav` file (omit to evaluate built-in 44.1 kHz PNW dawn chorus or `--ultrasonic`).
        #[arg(long)]
        wav: Option<PathBuf>,

        /// Use the built-in 192 kHz *Myotis* ultrasonic bat echolocation demo when `--wav` is omitted.
        #[arg(long, default_value_t = false)]
        ultrasonic: bool,

        /// STFT FFT window size (power of two, e.g. `512`, `1024`, `2048`).
        #[arg(long, default_value_t = 1024)]
        fft_size: usize,

        /// Anthrophony lower frequency bound (`Hz`) for Kasten et al. (2012) `NDSI`.
        #[arg(long, default_value_t = 1_000.0)]
        anthro_min_hz: f64,

        /// Anthrophony upper frequency bound (`Hz`) for Kasten et al. (2012) `NDSI`.
        #[arg(long, default_value_t = 2_000.0)]
        anthro_max_hz: f64,

        /// Biophony lower frequency bound (`Hz`) for `NDSI` & `ACI`.
        #[arg(long, default_value_t = 2_000.0)]
        bio_min_hz: f64,

        /// Biophony upper frequency bound (`Hz`) for `NDSI` & `ACI`.
        #[arg(long, default_value_t = 11_000.0)]
        bio_max_hz: f64,
    },

    /// Detect vocalizations in a target frequency band and export Raven Pro `.tsv` selections, `.svg` spectrogram, and sliced `.wav` clips.
    Slice {
        /// Input `.wav` recording (omit to run on built-in PNW dawn-chorus or `--ultrasonic` demo).
        #[arg(long)]
        wav: Option<PathBuf>,

        /// Use built-in 192 kHz ultrasonic *Myotis* demo when `--wav` is omitted.
        #[arg(long, default_value_t = false)]
        ultrasonic: bool,

        /// Target band lower frequency (`Hz`).
        #[arg(long, default_value_t = 1_800.0)]
        low_hz: f64,

        /// Target band upper frequency (`Hz`).
        #[arg(long, default_value_t = 8_000.0)]
        high_hz: f64,

        /// Adaptive SNR trigger (`dB`) above band median noise floor.
        #[arg(long, default_value_t = 12.0)]
        snr_db: f64,

        /// Minimum call duration (`ms`).
        #[arg(long, default_value_t = 100.0)]
        min_ms: f64,

        /// Maximum inter-syllable gap (`ms`) merged into one phrase.
        #[arg(long, default_value_t = 95.0)]
        merge_gap_ms: f64,

        /// Apply zero-phase Butterworth bandpass filter to exported `.wav` clips.
        #[arg(long, default_value_t = false)]
        bandpass: bool,

        /// Taxon or vocalization label written to Raven `Annotation` column.
        #[arg(long, default_value = "Target Vocalization")]
        species: String,

        /// Optional directory to write individual sliced `.wav` call clips.
        #[arg(long)]
        clips_dir: Option<PathBuf>,

        /// Optional path to write Cornell Raven Pro 1.6 Selection Table (`.txt`/`.tsv`).
        #[arg(long)]
        raven_out: Option<PathBuf>,

        /// Optional path to write annotated `.svg` spectrogram.
        #[arg(long)]
        svg_out: Option<PathBuf>,

        /// Color palette for `.svg` spectrogram.
        #[arg(long, value_enum, default_value_t = CliPalette::Magma)]
        palette: CliPalette,
    },

    /// Write a built-in synthetic `.wav` field recording (44.1 kHz dawn chorus or 192 kHz ultrasonic bat pass) to disk.
    Demo {
        /// Destination `.wav` file path.
        #[arg(long)]
        out: PathBuf,

        /// Generate the 192 kHz ultrasonic *Myotis* echolocation demo instead of the 44.1 kHz dawn chorus.
        #[arg(long, default_value_t = false)]
        ultrasonic: bool,
    },
}

/// Executes the `osabt ecosound` CLI workflow.
pub fn run_ecosound_cli(args: EcosoundArgs) -> Result<(), String> {
    match args.command {
        None => launch_default_gui(None),
        Some(EcosoundSubcommand::Gui { wav }) => launch_default_gui(wav),
        Some(EcosoundSubcommand::Indices {
            wav,
            ultrasonic,
            fft_size,
            anthro_min_hz,
            anthro_max_hz,
            bio_min_hz,
            bio_max_hz,
        }) => {
            let rec = load_or_demo_wav(wav.as_ref(), ultrasonic)?;
            let spec = compute_stft_spectrogram(
                &rec,
                &SpectrogramConfig {
                    fft_size,
                    hop_size: (fft_size / 4).max(1),
                    window: WindowFunction::Hann,
                    ..SpectrogramConfig::default()
                },
            )?;
            let idx = compute_ecoacoustic_indices(
                &rec,
                &spec,
                (anthro_min_hz, anthro_max_hz),
                (bio_min_hz, bio_max_hz),
            );

            println!("=== OSABT EcoSound — Soundscape Biodiversity Indices ===");
            println!(
                "Recording              : {} ({:.2} s @ {} Hz, {} ch)",
                rec.title, idx.duration_seconds, idx.sample_rate_hz, rec.channels
            );
            if let Some(cmt) = &rec.comment_metadata {
                println!("Metadata (ICMT)        : {cmt}");
            }
            println!(
                "ACI (Pieretti 2011)    : {:.2} total  |  {:.2} / min",
                idx.aci_total, idx.aci_per_minute
            );
            println!(
                "NDSI (Kasten 2012)     : {:+.4}  (Biophony PSD={:.3e}, Anthrophony PSD={:.3e})",
                idx.ndsi, idx.biophony_psd, idx.anthrophony_psd
            );
            println!(
                "Entropy H (Sueur 2008) : {:.4}  (Ht={:.4} × Hf={:.4})",
                idx.acoustic_entropy_h, idx.temporal_entropy_ht, idx.spectral_entropy_hf
            );
            println!(
                "ADI / AEI (Villanueva) : ADI={:.4} (Shannon)  |  AEI={:.4} (Gini)",
                idx.adi_shannon, idx.aei_gini
            );
            Ok(())
        }
        Some(EcosoundSubcommand::Slice {
            wav,
            ultrasonic,
            low_hz,
            high_hz,
            snr_db,
            min_ms,
            merge_gap_ms,
            bandpass,
            species,
            clips_dir,
            raven_out,
            svg_out,
            palette,
        }) => {
            let rec = load_or_demo_wav(wav.as_ref(), ultrasonic)?;
            let spec = compute_stft_spectrogram(&rec, &SpectrogramConfig::default())?;
            let cfg = CallSlicerConfig {
                low_freq_hz: low_hz,
                high_freq_hz: high_hz,
                snr_threshold_db: snr_db,
                min_duration_ms: min_ms,
                merge_gap_ms,
                bandpass_clips: bandpass,
                species_label: species,
                ..CallSlicerConfig::default()
            };

            let events = detect_calls_in_band(&spec, &cfg);
            println!(
                "Detected {} vocalizations in [{:.0}, {:.0}] Hz on '{}':",
                events.len(),
                low_hz,
                high_hz,
                rec.title
            );
            for ev in &events {
                println!(
                    "  #{:03}  {:.3}–{:.3} s ({:.0} ms) | {:.0}–{:.0} Hz | Peak {:.0} Hz | SNR {:.1} dB | Slope {:+.0} Hz/s",
                    ev.selection_id,
                    ev.begin_time_s,
                    ev.end_time_s,
                    ev.duration_ms,
                    ev.low_freq_hz,
                    ev.high_freq_hz,
                    ev.peak_freq_hz,
                    ev.snr_db,
                    ev.fm_slope_hz_per_s
                );
            }

            if let Some(rp) = raven_out {
                let tsv = export_raven_selection_tsv(&events, &rec.title);
                fs::write(&rp, tsv)
                    .map_err(|e| format!("Failed writing Raven table '{}': {e}", rp.display()))?;
                println!("Wrote Raven Pro Selection Table -> {}", rp.display());
            }

            if let Some(sp) = svg_out {
                let svg = render_spectrogram_svg(
                    &spec,
                    &events,
                    Some((low_hz, high_hz)),
                    palette.into(),
                    &rec.title,
                );
                fs::write(&sp, svg)
                    .map_err(|e| format!("Failed writing spectrogram SVG '{}': {e}", sp.display()))?;
                println!("Wrote annotated spectrogram SVG -> {}", sp.display());
            }

            if let Some(dir) = clips_dir {
                fs::create_dir_all(&dir)
                    .map_err(|e| format!("Failed creating '{}': {e}", dir.display()))?;
                let clips = slice_call_wav_clips(&rec, &events, &cfg, "ecosound");
                for (fname, bytes) in &clips {
                    let p = dir.join(fname);
                    fs::write(&p, bytes)
                        .map_err(|e| format!("Failed writing clip '{}': {e}", p.display()))?;
                }
                println!("Exported {} sliced .wav call clips -> {}", clips.len(), dir.display());
            }

            Ok(())
        }
        Some(EcosoundSubcommand::Demo { out, ultrasonic }) => {
            let rec = if ultrasonic {
                demo_ultrasonic_bat_recording()
            } else {
                demo_dawn_chorus_recording()
            };
            let bytes = encode_wav_pcm16(
                rec.sample_rate_hz,
                &rec.samples_mono,
                rec.comment_metadata.as_deref(),
            );
            fs::write(&out, bytes)
                .map_err(|e| format!("Failed writing demo WAV '{}': {e}", out.display()))?;
            println!(
                "Wrote synthetic bioacoustic recording '{}' ({:.2} s @ {} Hz) -> {}",
                rec.title,
                rec.duration_seconds(),
                rec.sample_rate_hz,
                out.display()
            );
            Ok(())
        }
    }
}

fn load_or_demo_wav(path: Option<&PathBuf>, ultrasonic: bool) -> Result<WavRecording, String> {
    if let Some(p) = path {
        let bytes = fs::read(p).map_err(|e| format!("Failed reading '{}': {e}", p.display()))?;
        let name = p
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or("recording.wav");
        parse_wav_bytes(name, &bytes)
    } else if ultrasonic {
        Ok(demo_ultrasonic_bat_recording())
    } else {
        Ok(demo_dawn_chorus_recording())
    }
}

#[cfg(feature = "gui")]
fn launch_default_gui(wav: Option<PathBuf>) -> Result<(), String> {
    crate::ecosound_gui::run_ecosound_gui(wav)
}

#[cfg(not(feature = "gui"))]
fn launch_default_gui(_wav: Option<PathBuf>) -> Result<(), String> {
    Err(String::from(
        "OSABT was compiled without the `gui` feature; rebuild with `--features gui` or use `osabt ecosound indices` / `slice`.",
    ))
}
