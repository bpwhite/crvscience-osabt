// SPDX-License-Identifier: GPL-3.0-or-later
// Copyright (C) 2026 Bryan P. White — CRV Science Network (https://crvscience.com)
// Part of Open Source Ad-hoc Biology Tools (OSABT)

//! Tool #10 (`osabt ecosound`) — Zero-dependency RIFF/WAVE (`.wav` + AudioMoth `LIST`/`INFO`)
//! parser & encoder, Cooley–Tukey radix-2 FFT, Short-Time Fourier Transform (`STFT`) spectrogram
//! engine, Butterworth biquad bandpass filter, and peer-reviewed soundscape ecology indices:
//! - **ACI** — Pieretti, Farina & Morri (2011) *Ecological Indicators* 11(3):868–873
//! - **NDSI** — Kasten, Gage, Fox & Joo (2012) *Ecological Informatics* 7(1):50–67
//! - **Acoustic Entropy (`H = Ht * Hf`)** — Sueur et al. (2008) *PLOS ONE* 3(12):e4065
//! - **ADI & AEI** — Villanueva-Rivera et al. (2011) *Landscape Ecology* 26(9):1233–1246

use std::f64::consts::PI;

/// Decoded mono-normalized waveform + header metadata from a RIFF/WAVE recording.
#[derive(Debug, Clone)]
pub struct WavRecording {
    pub title: String,
    pub sample_rate_hz: u32,
    pub channels: u16,
    pub bits_per_sample: u16,
    /// Mono-summed waveform normalized to `[-1.0, 1.0]`.
    pub samples_mono: Vec<f32>,
    /// Optional `IART` (e.g. AudioMoth device ID) from `LIST`/`INFO` chunk.
    pub artist_metadata: Option<String>,
    /// Optional `ICMT` (e.g. AudioMoth deployment timestamp/gain/battery) from `LIST`/`INFO` chunk.
    pub comment_metadata: Option<String>,
}

impl WavRecording {
    pub fn duration_seconds(&self) -> f64 {
        if self.sample_rate_hz == 0 {
            0.0
        } else {
            (self.samples_mono.len() as f64) / (self.sample_rate_hz as f64)
        }
    }

    pub fn nyquist_hz(&self) -> f64 {
        0.5 * (self.sample_rate_hz as f64)
    }
}

/// Spectral windowing functions for STFT frame tapering.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub enum WindowFunction {
    Hann,
    Hamming,
    BlackmanHarris,
}

impl WindowFunction {
    pub fn label(self) -> &'static str {
        match self {
            Self::Hann => "Hann",
            Self::Hamming => "Hamming",
            Self::BlackmanHarris => "Blackman-Harris (4-term)",
        }
    }

    pub fn weights(self, n: usize) -> Vec<f64> {
        if n <= 1 {
            return vec![1.0; n];
        }
        let denom = (n - 1) as f64;
        (0..n)
            .map(|i| {
                let phase = 2.0 * PI * (i as f64) / denom;
                match self {
                    Self::Hann => 0.5 * (1.0 - phase.cos()),
                    Self::Hamming => 0.54 - 0.46 * phase.cos(),
                    Self::BlackmanHarris => {
                        0.35875 - 0.48829 * phase.cos() + 0.14128 * (2.0 * phase).cos()
                            - 0.01168 * (3.0 * phase).cos()
                    }
                }
            })
            .collect()
    }
}

/// Parameters controlling STFT spectrogram computation.
#[derive(Debug, Clone)]
pub struct SpectrogramConfig {
    /// Power-of-two FFT frame size (e.g. `512`, `1024`, `2048`).
    pub fft_size: usize,
    /// Hop size in samples between consecutive frames (e.g. `256` or `512`).
    pub hop_size: usize,
    pub window: WindowFunction,
    pub min_freq_hz: f64,
    /// Upper frequency display bound (`0.0` = auto-set to Nyquist).
    pub max_freq_hz: f64,
    /// Dynamic range floor in `dBFS` (default `-90.0`).
    pub floor_dbfs: f64,
}

impl Default for SpectrogramConfig {
    fn default() -> Self {
        Self {
            fft_size: 1024,
            hop_size: 256,
            window: WindowFunction::Hann,
            min_freq_hz: 0.0,
            max_freq_hz: 0.0,
            floor_dbfs: -90.0,
        }
    }
}

/// Time-frequency STFT spectrogram matrix (`row-major: frame_idx * num_bins + bin_idx`).
#[derive(Debug, Clone)]
pub struct AcousticSpectrogram {
    pub sample_rate_hz: u32,
    pub fft_size: usize,
    pub hop_size: usize,
    pub num_frames: usize,
    pub num_bins: usize,
    pub frame_times_s: Vec<f64>,
    pub bin_freqs_hz: Vec<f64>,
    /// One-sided linear amplitude spectrum per `[frame, bin]`.
    pub mag_linear: Vec<f64>,
    /// One-sided power level in `dBFS` (`20 * log10(mag)`, clamped to `floor_dbfs..=0.0`).
    pub power_dbfs: Vec<f64>,
    pub floor_dbfs: f64,
}

impl AcousticSpectrogram {
    #[inline]
    pub fn dbfs_at(&self, frame: usize, bin: usize) -> f64 {
        self.power_dbfs[frame * self.num_bins + bin]
    }

    #[inline]
    pub fn mag_at(&self, frame: usize, bin: usize) -> f64 {
        self.mag_linear[frame * self.num_bins + bin]
    }
}

/// Soundscape ecology & rapid biodiversity survey indices computed on a recording.
#[derive(Debug, Clone)]
#[cfg_attr(feature = "serde", derive(serde::Serialize, serde::Deserialize))]
pub struct EcoacousticIndices {
    pub duration_seconds: f64,
    pub sample_rate_hz: u32,
    /// Pieretti, Farina & Morri (2011) total Acoustic Complexity Index across biophony bins.
    pub aci_total: f64,
    /// Duration-standardized ACI per minute (`aci_total * 60 / duration_seconds`).
    pub aci_per_minute: f64,
    /// Kasten et al. (2012) Normalized Difference Soundscape Index `(biophony - anthrophony) / (biophony + anthrophony)` in `[-1, +1]`.
    pub ndsi: f64,
    pub anthrophony_psd: f64,
    pub biophony_psd: f64,
    /// Sueur et al. (2008) temporal amplitude-envelope Shannon entropy `Ht` in `[0, 1]`.
    pub temporal_entropy_ht: f64,
    /// Sueur et al. (2008) spectral mean-energy Shannon entropy `Hf` in `[0, 1]`.
    pub spectral_entropy_hf: f64,
    /// Sueur et al. (2008) composite Acoustic Entropy Index `H = Ht * Hf` in `[0, 1]`.
    pub acoustic_entropy_h: f64,
    /// Villanueva-Rivera et al. (2011) Acoustic Diversity Index (Shannon diversity over 1 kHz bins above `-50 dBFS`).
    pub adi_shannon: f64,
    /// Villanueva-Rivera et al. (2011) Acoustic Evenness Index (Gini coefficient over 1 kHz bins above `-50 dBFS`).
    pub aei_gini: f64,
}

/// Parses a RIFF/WAVE byte buffer (`8`/`16`/`24`/`32`-bit PCM or `32`-bit IEEE float, mono or multi-channel)
/// including optional AudioMoth / Wildlife Acoustics `LIST`/`INFO` metadata tags.
pub fn parse_wav_bytes(title: &str, bytes: &[u8]) -> Result<WavRecording, String> {
    if bytes.len() < 44 || &bytes[0..4] != b"RIFF" || &bytes[8..12] != b"WAVE" {
        return Err(String::from("Input is not a valid RIFF/WAVE (.wav) audio file."));
    }

    let mut pos = 12usize;
    let mut format_tag = 0u16;
    let mut channels = 0u16;
    let mut sample_rate = 0u32;
    let mut bits_per_sample = 0u16;
    let mut data_slice: Option<&[u8]> = None;
    let mut artist_meta: Option<String> = None;
    let mut comment_meta: Option<String> = None;

    while pos + 8 <= bytes.len() {
        let chunk_id = &bytes[pos..pos + 4];
        let chunk_len = u32::from_le_bytes([
            bytes[pos + 4],
            bytes[pos + 5],
            bytes[pos + 6],
            bytes[pos + 7],
        ]) as usize;
        let payload_start = pos + 8;
        let payload_end = payload_start.saturating_add(chunk_len).min(bytes.len());
        let payload = &bytes[payload_start..payload_end];

        match chunk_id {
            b"fmt " => {
                if payload.len() < 16 {
                    return Err(String::from("WAV fmt chunk is truncated (< 16 bytes)."));
                }
                let raw_tag = u16::from_le_bytes([payload[0], payload[1]]);
                channels = u16::from_le_bytes([payload[2], payload[3]]);
                sample_rate =
                    u32::from_le_bytes([payload[4], payload[5], payload[6], payload[7]]);
                bits_per_sample = u16::from_le_bytes([payload[14], payload[15]]);
                if raw_tag == 0xFFFE && payload.len() >= 26 {
                    format_tag = u16::from_le_bytes([payload[24], payload[25]]);
                } else {
                    format_tag = raw_tag;
                }
            }
            b"data" => {
                data_slice = Some(payload);
            }
            b"LIST" => {
                if payload.len() >= 4 && &payload[0..4] == b"INFO" {
                    let mut ip = 4usize;
                    while ip + 8 <= payload.len() {
                        let sub_id = &payload[ip..ip + 4];
                        let sub_len = u32::from_le_bytes([
                            payload[ip + 4],
                            payload[ip + 5],
                            payload[ip + 6],
                            payload[ip + 7],
                        ]) as usize;
                        let s0 = ip + 8;
                        let s1 = s0.saturating_add(sub_len).min(payload.len());
                        let txt = String::from_utf8_lossy(&payload[s0..s1])
                            .trim_matches(char::from(0))
                            .trim()
                            .to_string();
                        if !txt.is_empty() {
                            if sub_id == b"IART" {
                                artist_meta = Some(txt);
                            } else if sub_id == b"ICMT" {
                                comment_meta = Some(txt);
                            }
                        }
                        ip = s1 + (sub_len % 2);
                    }
                }
            }
            _ => {}
        }

        pos = payload_end + (chunk_len % 2);
    }

    if channels == 0 || sample_rate == 0 {
        return Err(String::from("WAV file missing valid fmt chunk."));
    }
    let raw_data = data_slice.ok_or_else(|| String::from("WAV file missing data chunk."))?;

    let bytes_per_sample = (bits_per_sample as usize) / 8;
    if bytes_per_sample == 0 {
        return Err(format!("Unsupported bits_per_sample = {bits_per_sample}."));
    }
    let frame_stride = bytes_per_sample * (channels as usize);
    let total_frames = raw_data.len() / frame_stride;
    let mut samples_mono = Vec::with_capacity(total_frames);
    let inv_ch = 1.0_f32 / (channels as f32);

    for f in 0..total_frames {
        let base = f * frame_stride;
        let mut acc = 0.0_f32;
        for ch in 0..(channels as usize) {
            let off = base + ch * bytes_per_sample;
            let s: f32 = match (format_tag, bits_per_sample) {
                (1, 8) => ((raw_data[off] as f32) - 128.0) / 128.0,
                (1, 16) => {
                    let v = i16::from_le_bytes([raw_data[off], raw_data[off + 1]]);
                    (v as f32) / 32768.0
                }
                (1, 24) => {
                    let b0 = raw_data[off] as i32;
                    let b1 = (raw_data[off + 1] as i32) << 8;
                    let b2 = ((raw_data[off + 2] as i8) as i32) << 16;
                    ((b0 | b1 | b2) as f32) / 8_388_608.0
                }
                (1, 32) => {
                    let v = i32::from_le_bytes([
                        raw_data[off],
                        raw_data[off + 1],
                        raw_data[off + 2],
                        raw_data[off + 3],
                    ]);
                    (v as f32) / 2_147_483_648.0
                }
                (3, 32) => f32::from_le_bytes([
                    raw_data[off],
                    raw_data[off + 1],
                    raw_data[off + 2],
                    raw_data[off + 3],
                ])
                .clamp(-1.0, 1.0),
                _ => {
                    return Err(format!(
                        "Unsupported WAV encoding (format_tag={format_tag}, bits={bits_per_sample})."
                    ));
                }
            };
            acc += s;
        }
        samples_mono.push(acc * inv_ch);
    }

    Ok(WavRecording {
        title: title.to_string(),
        sample_rate_hz: sample_rate,
        channels,
        bits_per_sample,
        samples_mono,
        artist_metadata: artist_meta,
        comment_metadata: comment_meta,
    })
}

/// Encodes mono normalized `[-1.0, 1.0]` samples into a canonical 16-bit PCM RIFF/WAVE `.wav` buffer,
/// optionally embedding an `ICMT` provenance comment in a `LIST`/`INFO` metadata chunk.
pub fn encode_wav_pcm16(sample_rate_hz: u32, samples: &[f32], comment: Option<&str>) -> Vec<u8> {
    let data_bytes_len = samples.len() * 2;
    let mut info_chunk = Vec::new();
    if let Some(cmt) = comment {
        if !cmt.is_empty() {
            let mut c_bytes = cmt.as_bytes().to_vec();
            c_bytes.push(0);
            let c_len = c_bytes.len() as u32;
            info_chunk.extend_from_slice(b"INFO");
            info_chunk.extend_from_slice(b"ICMT");
            info_chunk.extend_from_slice(&c_len.to_le_bytes());
            info_chunk.extend_from_slice(&c_bytes);
            if c_bytes.len() % 2 == 1 {
                info_chunk.push(0);
            }
        }
    }

    let list_total = if info_chunk.is_empty() {
        0usize
    } else {
        8 + info_chunk.len()
    };
    let riff_size = (36 + data_bytes_len + list_total) as u32;
    let mut out = Vec::with_capacity(8 + riff_size as usize);

    out.extend_from_slice(b"RIFF");
    out.extend_from_slice(&riff_size.to_le_bytes());
    out.extend_from_slice(b"WAVE");

    out.extend_from_slice(b"fmt ");
    out.extend_from_slice(&16u32.to_le_bytes());
    out.extend_from_slice(&1u16.to_le_bytes());
    out.extend_from_slice(&1u16.to_le_bytes());
    out.extend_from_slice(&sample_rate_hz.to_le_bytes());
    let byte_rate = sample_rate_hz * 2;
    out.extend_from_slice(&byte_rate.to_le_bytes());
    out.extend_from_slice(&2u16.to_le_bytes());
    out.extend_from_slice(&16u16.to_le_bytes());

    if !info_chunk.is_empty() {
        out.extend_from_slice(b"LIST");
        out.extend_from_slice(&(info_chunk.len() as u32).to_le_bytes());
        out.extend_from_slice(&info_chunk);
    }

    out.extend_from_slice(b"data");
    out.extend_from_slice(&(data_bytes_len as u32).to_le_bytes());
    for &s in samples {
        let q = (s.clamp(-1.0, 1.0) * 32767.0).round() as i16;
        out.extend_from_slice(&q.to_le_bytes());
    }
    out
}

/// In-place Cooley–Tukey radix-2 Fast Fourier Transform (`re.len()` must be a power of two).
pub fn fft_radix2_inplace(re: &mut [f64], im: &mut [f64], inverse: bool) {
    let n = re.len();
    assert_eq!(n, im.len(), "re and im slices must have equal length");
    assert!(n.is_power_of_two(), "FFT size must be a power of two");
    if n <= 1 {
        return;
    }

    let mut j = 0usize;
    for i in 1..n {
        let mut bit = n >> 1;
        while j & bit != 0 {
            j ^= bit;
            bit >>= 1;
        }
        j ^= bit;
        if i < j {
            re.swap(i, j);
            im.swap(i, j);
        }
    }

    let sign = if inverse { 1.0 } else { -1.0 };
    let mut len = 2usize;
    while len <= n {
        let half = len / 2;
        let ang = sign * 2.0 * PI / (len as f64);
        let wlen_re = ang.cos();
        let wlen_im = ang.sin();

        let mut start = 0usize;
        while start < n {
            let mut w_re = 1.0_f64;
            let mut w_im = 0.0_f64;
            for k in 0..half {
                let u_re = re[start + k];
                let u_im = im[start + k];
                let v_re = re[start + k + half] * w_re - im[start + k + half] * w_im;
                let v_im = re[start + k + half] * w_im + im[start + k + half] * w_re;

                re[start + k] = u_re + v_re;
                im[start + k] = u_im + v_im;
                re[start + k + half] = u_re - v_re;
                im[start + k + half] = u_im - v_im;

                let next_w_re = w_re * wlen_re - w_im * wlen_im;
                w_im = w_re * wlen_im + w_im * wlen_re;
                w_re = next_w_re;
            }
            start += len;
        }
        len <<= 1;
    }

    if inverse {
        let inv_n = 1.0 / (n as f64);
        for k in 0..n {
            re[k] *= inv_n;
            im[k] *= inv_n;
        }
    }
}

/// Computes a one-sided Short-Time Fourier Transform (`STFT`) spectrogram in linear magnitude and `dBFS`.
pub fn compute_stft_spectrogram(
    rec: &WavRecording,
    cfg: &SpectrogramConfig,
) -> Result<AcousticSpectrogram, String> {
    let n = cfg.fft_size;
    if !n.is_power_of_two() || n < 64 {
        return Err(format!(
            "fft_size ({n}) must be a power of two >= 64 (e.g. 512, 1024, 2048)."
        ));
    }
    let hop = cfg.hop_size.max(1);
    if rec.samples_mono.len() < n {
        return Err(format!(
            "Recording has {} samples, shorter than fft_size ({n}).",
            rec.samples_mono.len()
        ));
    }

    let sr = rec.sample_rate_hz as f64;
    let nyq = 0.5 * sr;
    let max_f = if cfg.max_freq_hz <= 0.0 || cfg.max_freq_hz > nyq {
        nyq
    } else {
        cfg.max_freq_hz
    };
    let min_f = cfg.min_freq_hz.clamp(0.0, max_f);

    let full_bins = n / 2 + 1;
    let df = sr / (n as f64);
    let b_start = ((min_f / df).floor() as usize).min(full_bins - 1);
    let b_end = ((max_f / df).ceil() as usize).clamp(b_start + 1, full_bins - 1);
    let num_bins = b_end - b_start + 1;
    let bin_freqs_hz: Vec<f64> = (b_start..=b_end).map(|k| (k as f64) * df).collect();

    let win = cfg.window.weights(n);
    let win_sum: f64 = win.iter().sum::<f64>().max(1e-9);
    let norm = 2.0 / win_sum;

    let num_frames = (rec.samples_mono.len() - n) / hop + 1;
    let mut frame_times_s = Vec::with_capacity(num_frames);
    let mut mag_linear = Vec::with_capacity(num_frames * num_bins);
    let mut power_dbfs = Vec::with_capacity(num_frames * num_bins);

    let mut re = vec![0.0_f64; n];
    let mut im = vec![0.0_f64; n];

    for f_idx in 0..num_frames {
        let s0 = f_idx * hop;
        let center_s = ((s0 as f64) + 0.5 * (n as f64)) / sr;
        frame_times_s.push(center_s);

        for i in 0..n {
            re[i] = (rec.samples_mono[s0 + i] as f64) * win[i];
            im[i] = 0.0;
        }
        fft_radix2_inplace(&mut re, &mut im, false);

        for k in b_start..=b_end {
            let mag = (re[k] * re[k] + im[k] * im[k]).sqrt() * norm;
            let db = if mag > 1e-12 {
                (20.0 * mag.log10()).clamp(cfg.floor_dbfs, 0.0)
            } else {
                cfg.floor_dbfs
            };
            mag_linear.push(mag);
            power_dbfs.push(db);
        }
    }

    Ok(AcousticSpectrogram {
        sample_rate_hz: rec.sample_rate_hz,
        fft_size: n,
        hop_size: hop,
        num_frames,
        num_bins,
        frame_times_s,
        bin_freqs_hz,
        mag_linear,
        power_dbfs,
        floor_dbfs: cfg.floor_dbfs,
    })
}

/// Applies a cascaded 2nd-order Butterworth high-pass + low-pass zero-phase (`filtfilt`-style forward-backward)
/// biquad bandpass filter to isolate vocalizations in `[low_hz, high_hz]`.
pub fn apply_bandpass_biquad(
    samples: &[f32],
    sample_rate_hz: u32,
    low_hz: f64,
    high_hz: f64,
) -> Vec<f32> {
    if samples.is_empty() || sample_rate_hz == 0 {
        return Vec::new();
    }
    let sr = sample_rate_hz as f64;
    let nyq = 0.495 * sr;
    let mut buf: Vec<f64> = samples.iter().map(|&x| x as f64).collect();

    if low_hz > 5.0 && low_hz < nyq {
        let coeffs = butterworth_biquad_coeffs(low_hz / sr, true);
        run_biquad_in_place(&mut buf, coeffs);
        buf.reverse();
        run_biquad_in_place(&mut buf, coeffs);
        buf.reverse();
    }

    if high_hz > low_hz && high_hz < nyq {
        let coeffs = butterworth_biquad_coeffs(high_hz / sr, false);
        run_biquad_in_place(&mut buf, coeffs);
        buf.reverse();
        run_biquad_in_place(&mut buf, coeffs);
        buf.reverse();
    }

    buf.into_iter()
        .map(|v| (v.clamp(-1.0, 1.0)) as f32)
        .collect()
}

fn butterworth_biquad_coeffs(norm_fc: f64, highpass: bool) -> [f64; 5] {
    let w0 = 2.0 * PI * norm_fc.clamp(1e-4, 0.495);
    let cos_w0 = w0.cos();
    let sin_w0 = w0.sin();
    let alpha = sin_w0 / std::f64::consts::SQRT_2;
    let a0 = 1.0 + alpha;
    let (b0, b1, b2) = if highpass {
        let c = 0.5 * (1.0 + cos_w0);
        (c / a0, -(1.0 + cos_w0) / a0, c / a0)
    } else {
        let c = 0.5 * (1.0 - cos_w0);
        (c / a0, (1.0 - cos_w0) / a0, c / a0)
    };
    let a1 = (-2.0 * cos_w0) / a0;
    let a2 = (1.0 - alpha) / a0;
    [b0, b1, b2, a1, a2]
}

fn run_biquad_in_place(buf: &mut [f64], c: [f64; 5]) {
    let [b0, b1, b2, a1, a2] = c;
    let (mut x1, mut x2, mut y1, mut y2) = (0.0, 0.0, 0.0, 0.0);
    for sample in buf.iter_mut() {
        let x0 = *sample;
        let y0 = b0 * x0 + b1 * x1 + b2 * x2 - a1 * y1 - a2 * y2;
        x2 = x1;
        x1 = x0;
        y2 = y1;
        y1 = y0;
        *sample = y0;
    }
}

/// Computes `ACI` (Pieretti et al. 2011), `NDSI` (Kasten et al. 2012), `H = Ht * Hf` (Sueur et al. 2008),
/// and `ADI`/`AEI` (Villanueva-Rivera et al. 2011) on `rec` and its `STFT` `spec`.
pub fn compute_ecoacoustic_indices(
    rec: &WavRecording,
    spec: &AcousticSpectrogram,
    anthro_band_hz: (f64, f64),
    bio_band_hz: (f64, f64),
) -> EcoacousticIndices {
    let dur = rec.duration_seconds().max(1e-6);
    let nyq = rec.nyquist_hz();
    let bio_lo = bio_band_hz.0.max(0.0);
    let bio_hi = if bio_band_hz.1 <= bio_lo || bio_band_hz.1 > nyq {
        nyq
    } else {
        bio_band_hz.1
    };

    let dt = (spec.hop_size as f64) / (spec.sample_rate_hz as f64);
    let clump_frames = ((5.0 / dt).round() as usize).max(2).min(spec.num_frames.max(2));
    let mut aci_total = 0.0_f64;
    if spec.num_frames >= 2 {
        for (b, &f_hz) in spec.bin_freqs_hz.iter().enumerate() {
            if f_hz < bio_lo || f_hz > bio_hi {
                continue;
            }
            let mut f0 = 0usize;
            while f0 + 1 < spec.num_frames {
                let f1 = (f0 + clump_frames).min(spec.num_frames);
                let mut diff_sum = 0.0_f64;
                let mut amp_sum = 0.0_f64;
                for t in f0..f1 {
                    let cur = spec.mag_at(t, b);
                    amp_sum += cur;
                    if t + 1 < f1 {
                        diff_sum += (cur - spec.mag_at(t + 1, b)).abs();
                    }
                }
                if amp_sum > 1e-12 {
                    aci_total += diff_sum / amp_sum;
                }
                f0 = f1;
            }
        }
    }
    let aci_per_minute = aci_total * (60.0 / dur);

    let mut anthro_sum = 0.0_f64;
    let mut anthro_bins = 0usize;
    let mut bio_sum = 0.0_f64;
    let mut bio_bins = 0usize;

    for (b, &f_hz) in spec.bin_freqs_hz.iter().enumerate() {
        let mut mean_psd = 0.0_f64;
        for t in 0..spec.num_frames {
            let m = spec.mag_at(t, b);
            mean_psd += m * m;
        }
        mean_psd /= spec.num_frames.max(1) as f64;

        if f_hz >= anthro_band_hz.0 && f_hz < anthro_band_hz.1 {
            anthro_sum += mean_psd;
            anthro_bins += 1;
        } else if f_hz >= bio_lo && f_hz <= bio_hi {
            bio_sum += mean_psd;
            bio_bins += 1;
        }
    }

    let anthro_psd = if anthro_bins > 0 {
        anthro_sum / (anthro_bins as f64)
    } else {
        0.0
    };
    let bio_psd = if bio_bins > 0 {
        bio_sum / (bio_bins as f64)
    } else {
        0.0
    };
    let ndsi = if anthro_psd + bio_psd > 1e-18 {
        ((bio_psd - anthro_psd) / (bio_psd + anthro_psd)).clamp(-1.0, 1.0)
    } else {
        0.0
    };

    let env_win = ((0.005 * (rec.sample_rate_hz as f64)).round() as usize).max(1);
    let n_env = rec.samples_mono.len() / env_win;
    let temporal_entropy_ht = if n_env >= 2 {
        let mut env = Vec::with_capacity(n_env);
        let mut env_sum = 0.0_f64;
        for c in 0..n_env {
            let slice = &rec.samples_mono[c * env_win..(c + 1) * env_win];
            let rms = (slice.iter().map(|&s| (s as f64) * (s as f64)).sum::<f64>()
                / (env_win as f64))
                .sqrt();
            env.push(rms);
            env_sum += rms;
        }
        if env_sum > 1e-12 {
            let h_num: f64 = env
                .iter()
                .filter(|&&e| e > 1e-15)
                .map(|&e| {
                    let p = e / env_sum;
                    -p * p.ln()
                })
                .sum();
            (h_num / (n_env as f64).ln()).clamp(0.0, 1.0)
        } else {
            0.0
        }
    } else {
        0.0
    };

    let spectral_entropy_hf = if spec.num_bins > 2 {
        let mut mean_spec = Vec::with_capacity(spec.num_bins - 1);
        let mut spec_sum = 0.0_f64;
        for b in 1..spec.num_bins {
            let mut m_acc = 0.0_f64;
            for t in 0..spec.num_frames {
                m_acc += spec.mag_at(t, b);
            }
            let m_avg = m_acc / (spec.num_frames.max(1) as f64);
            mean_spec.push(m_avg);
            spec_sum += m_avg;
        }
        if spec_sum > 1e-12 {
            let h_num: f64 = mean_spec
                .iter()
                .filter(|&&v| v > 1e-15)
                .map(|&v| {
                    let p = v / spec_sum;
                    -p * p.ln()
                })
                .sum();
            (h_num / (mean_spec.len() as f64).ln()).clamp(0.0, 1.0)
        } else {
            0.0
        }
    } else {
        0.0
    };
    let acoustic_entropy_h = temporal_entropy_ht * spectral_entropy_hf;

    let top_adi_hz = nyq.min(10_000.0).max(1_000.0);
    let band_w = top_adi_hz / 10.0;
    let mut props = [0.0_f64; 10];
    for band_idx in 0..10 {
        let lo = (band_idx as f64) * band_w;
        let hi = lo + band_w;
        let mut total_cells = 0usize;
        let mut active_cells = 0usize;
        for (b, &f_hz) in spec.bin_freqs_hz.iter().enumerate() {
            if f_hz >= lo && f_hz < hi {
                for t in 0..spec.num_frames {
                    total_cells += 1;
                    if spec.dbfs_at(t, b) > -50.0 {
                        active_cells += 1;
                    }
                }
            }
        }
        props[band_idx] = if total_cells > 0 {
            (active_cells as f64) / (total_cells as f64)
        } else {
            0.0
        };
    }

    let prop_sum: f64 = props.iter().sum();
    let adi_shannon = if prop_sum > 1e-12 {
        props
            .iter()
            .filter(|&&v| v > 1e-12)
            .map(|&v| {
                let p = v / prop_sum;
                -p * p.ln()
            })
            .sum()
    } else {
        0.0
    };

    let mut sorted_p = props;
    sorted_p.sort_by(|a, b| a.partial_cmp(b).unwrap());
    let aei_gini = if prop_sum > 1e-12 {
        let n_b = sorted_p.len() as f64;
        let weighted: f64 = sorted_p
            .iter()
            .enumerate()
            .map(|(i, &v)| ((i + 1) as f64) * v)
            .sum();
        ((2.0 * weighted) / (n_b * prop_sum) - (n_b + 1.0) / n_b).clamp(0.0, 1.0)
    } else {
        1.0
    };

    EcoacousticIndices {
        duration_seconds: dur,
        sample_rate_hz: rec.sample_rate_hz,
        aci_total,
        aci_per_minute,
        ndsi,
        anthrophony_psd: anthro_psd,
        biophony_psd: bio_psd,
        temporal_entropy_ht,
        spectral_entropy_hf,
        acoustic_entropy_h,
        adi_shannon,
        aei_gini,
    }
}

/// Synthesizes a realistic 5.0-second `44,100 Hz` Pacific Northwest riparian dawn-chorus recording
/// featuring distinct avian & anuran vocalizations for immediate inspection and call slicing:
/// - `0.35–0.90 s`: *Poecile atricapillus* (Black-capped Chickadee) *fee-bee* whistle (`3.95 -> 3.45 kHz`)
/// - `1.30–2.05 s`: *Catharus ustulatus* (Swainson's Thrush) upward spiraling flute phrase (`2.5 -> 5.7 kHz`)
/// - `2.45–3.05 s`: *Pseudacris regilla* (Pacific Chorus Frog) pulsed harmonic call (`2.15 kHz` + `4.30 kHz`)
/// - `3.45–4.05 s`: *Poecile atricapillus* second *fee-bee* phrase (`3.92 -> 3.42 kHz`)
/// - `4.30–4.78 s`: *Setophaga townsendi* (Townsend's Warbler) high-frequency FM trill (`5.9 -> 7.2 kHz`)
pub fn demo_dawn_chorus_recording() -> WavRecording {
    let sr = 44_100u32;
    let dur_s = 5.0_f64;
    let n = (dur_s * (sr as f64)) as usize;
    let mut samples = vec![0.0_f32; n];

    let mut rng = 0xC0DE_CAFE_2026_1010u64;
    let mut next_noise = || -> f64 {
        rng = rng.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = rng;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^= z >> 31;
        ((z >> 11) as f64) / ((1u64 << 53) as f64) * 2.0 - 1.0
    };

    for (i, s) in samples.iter_mut().enumerate() {
        let t = (i as f64) / (sr as f64);
        let bg = 0.0022 * next_noise() + 0.0015 * (2.0 * PI * 1_420.0 * t).sin();
        *s = bg as f32;
    }

    let mut add_chirp = |t0: f64, t1: f64, f0: f64, f1: f64, amp: f64, am_hz: f64, h2_rel: f64| {
        let i0 = ((t0 * sr as f64).round() as usize).min(n);
        let i1 = ((t1 * sr as f64).round() as usize).min(n);
        if i1 <= i0 + 4 {
            return;
        }
        let span = (i1 - i0) as f64;
        let mut phase = 0.0_f64;
        for idx in i0..i1 {
            let u = ((idx - i0) as f64) / span;
            let env = (PI * u).sin().powf(0.65);
            let am = if am_hz > 0.0 {
                0.55 + 0.45 * (2.0 * PI * am_hz * ((idx - i0) as f64 / sr as f64)).sin()
            } else {
                1.0
            };
            let f_inst = f0 + (f1 - f0) * u;
            phase += 2.0 * PI * f_inst / (sr as f64);
            let sig = amp * env * am * (phase.sin() + h2_rel * (2.0 * phase).sin());
            samples[idx] = (samples[idx] + sig as f32).clamp(-0.98, 0.98);
        }
    };

    add_chirp(0.35, 0.60, 3_960.0, 3_910.0, 0.42, 0.0, 0.08);
    add_chirp(0.64, 0.90, 3_460.0, 3_430.0, 0.40, 0.0, 0.08);

    add_chirp(1.30, 1.46, 2_550.0, 3_150.0, 0.36, 0.0, 0.12);
    add_chirp(1.49, 1.66, 3_250.0, 4_100.0, 0.40, 18.0, 0.10);
    add_chirp(1.69, 1.86, 4_200.0, 5_050.0, 0.38, 22.0, 0.08);
    add_chirp(1.88, 2.05, 5_050.0, 5_750.0, 0.32, 26.0, 0.05);

    add_chirp(2.45, 3.05, 2_120.0, 2_480.0, 0.44, 34.0, 0.45);

    add_chirp(3.45, 3.70, 3_920.0, 3_880.0, 0.39, 0.0, 0.08);
    add_chirp(3.75, 4.02, 3_420.0, 3_390.0, 0.37, 0.0, 0.08);

    add_chirp(4.30, 4.78, 5_900.0, 7_250.0, 0.35, 42.0, 0.04);

    WavRecording {
        title: String::from("CRV_MarysRiver_DawnChorus_44kHz.wav"),
        sample_rate_hz: sr,
        channels: 1,
        bits_per_sample: 16,
        samples_mono: samples,
        artist_metadata: Some(String::from("AudioMoth 1.2.0 — Corvallis OR")),
        comment_metadata: Some(String::from(
            "Recorded at 2026-05-18 05:42:00 PDT by AudioMoth (Gain=Medium, 44.1 kHz)",
        )),
    }
}

/// Synthesizes a 1.5-second `192,000 Hz` full-spectrum ultrasonic chiropteran recording
/// containing 4 *Myotis lucifugus* downward FM echolocation pulses (`82 kHz -> 39 kHz`).
pub fn demo_ultrasonic_bat_recording() -> WavRecording {
    let sr = 192_000u32;
    let dur_s = 1.5_f64;
    let n = (dur_s * (sr as f64)) as usize;
    let mut samples = vec![0.0_f32; n];

    let mut rng = 0xBA75_EC40_2026_1920u64;
    for s in &mut samples {
        rng = rng.wrapping_mul(6_364_136_223_846_793_005).wrapping_add(1);
        let u = ((rng >> 32) as f64) / (u32::MAX as f64) * 2.0 - 1.0;
        *s = (0.002 * u) as f32;
    }

    for &t_start in &[0.18_f64, 0.49, 0.83, 1.16] {
        let pulse_dur = 0.035;
        let i0 = (t_start * sr as f64) as usize;
        let i1 = ((t_start + pulse_dur) * sr as f64).min(n as f64) as usize;
        let mut phase = 0.0_f64;
        for idx in i0..i1 {
            let u = ((idx - i0) as f64) / ((i1 - i0) as f64);
            let f_hz = 82_000.0 / (1.0 + 1.10 * u);
            let env = (PI * u).sin();
            phase += 2.0 * PI * f_hz / (sr as f64);
            samples[idx] = (samples[idx] + (0.52 * env * phase.sin()) as f32).clamp(-0.98, 0.98);
        }
    }

    WavRecording {
        title: String::from("CRV_Myotis_Ultrasonic_192kHz.wav"),
        sample_rate_hz: sr,
        channels: 1,
        bits_per_sample: 16,
        samples_mono: samples,
        artist_metadata: Some(String::from("AudioMoth 1.2.0 Ultrasonic")),
        comment_metadata: Some(String::from("192 kHz full-spectrum chiropteran pass")),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn wav_pcm16_roundtrip_and_fft_parseval() {
        let demo = demo_dawn_chorus_recording();
        let bytes = encode_wav_pcm16(
            demo.sample_rate_hz,
            &demo.samples_mono[..4096],
            demo.comment_metadata.as_deref(),
        );
        let parsed = parse_wav_bytes("roundtrip.wav", &bytes).unwrap();
        assert_eq!(parsed.sample_rate_hz, 44_100);
        assert_eq!(parsed.samples_mono.len(), 4096);
        assert!(parsed.comment_metadata.unwrap().contains("AudioMoth"));

        let mut re: Vec<f64> = (0..256).map(|i| (0.2 * i as f64).sin()).collect();
        let orig = re.clone();
        let mut im = vec![0.0; 256];
        fft_radix2_inplace(&mut re, &mut im, false);
        fft_radix2_inplace(&mut re, &mut im, true);
        for (a, b) in orig.iter().zip(&re) {
            assert!((a - b).abs() < 1e-11);
        }
    }

    #[test]
    fn ecoacoustic_indices_within_theoretical_bounds() {
        let demo = demo_dawn_chorus_recording();
        let spec = compute_stft_spectrogram(&demo, &SpectrogramConfig::default()).unwrap();
        let idx = compute_ecoacoustic_indices(&demo, &spec, (1_000.0, 2_000.0), (2_000.0, 11_000.0));
        assert!(idx.aci_total > 10.0);
        assert!(idx.ndsi > 0.50 && idx.ndsi <= 1.0);
        assert!(idx.acoustic_entropy_h > 0.05 && idx.acoustic_entropy_h <= 1.0);
        assert!(idx.adi_shannon > 0.5);
    }
}
