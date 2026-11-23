# OSABT Scientific Methods & Algorithms (`v2.0.0`)

Copyright (C) 2026 Bryan P. White — CRV Science Network (`https://crvscience.com`).  
Licensed under the **GNU General Public License v3.0 or later** (`GPL-3.0-or-later`).

This reference documents the mathematical formulations, physical units, boundary guards, and peer-reviewed citations implemented in the `osabt` (`v2.0.0`) pure-Rust library crate.

---

## Tool #1 — `bench` (Oligonucleotide Thermodynamics & Wet-Lab Chemistry)

- **Modules**: `thermo`, `structure`, `primer_design`, `gibson`, `bench_calc`
- **Unified Nearest-Neighbor Duplex Free Energy** (*SantaLucia 1998*):
  ```text
  ΔH° = ΔH°_init + sum_i ΔH°(step_i) + ΔH°_sym
  ΔS° = ΔS°_init + sum_i ΔS°(step_i) + ΔS°_sym
  ΔG°(37 °C) = ΔH° - 310.15 * ΔS°
  Tm(1 M Na+) = ΔH° / (ΔS° + R * ln(C_T / x)) - 273.15   (x = 4 nonself, x = 1 self-comp)
  ```
- **Monovalent & Divalent Salt Correction** (*Owczarzy et al. 2004, 2008*):
  Effective free `Mg2+` accounts for `[dNTP]` chelation (`Mg_free` via quadratic stoichiometric binding). Depending on the ratio `R = sqrt(Mg_free) / Monovalent`, Owczarzy et al. (2008) smoothly switches among monovalent-dominated, intermediate, and divalent-dominated `1/Tm` correction regimes. Phusion high-salt polymerase presets apply Schildkraut & Lifson (1965) scaling.
- **Secondary Structure & G-Quadruplex Screening** (*SantaLucia & Hicks 2004; Bedrat et al. 2016*):
  Intramolecular hairpin stems, 3'-OH extensible self/cross-dimers (`ΔG°_37`), and sliding-window `G4Hunter` polymerase-stall propensity scores.
- **Henderson-Hasselbalch Thermal Buffer Scaling** (`bench_calc`):
  `pKa(T) = pKa(25 °C) + (dpKa/dT) * (T - 25)` for Tris, HEPES, MES, and Sodium Phosphate; conjugate acid/base mass balance via `[Base]/[Acid] = 10^(pH - pKa(T))`.

---

## Tool #2 — `plasmid` (Vector I/O, Digests, Golden Gate & Trace Alignment)

- **Modules**: `plasmid_io`, `plasmid_digest`, `plasmid_svg`
- **Circular Wrap-Around Restriction Mapping & Electrophoresis**:
  Recognizes palindromic and Type IIS sites across circular origin boundaries (`0..L`), simulates agarose gel migration (`log10(bp)` piecewise-linear mobility against 1 kb Plus / 100 bp ladders).
- **Golden Gate 4-bp Overhang Fidelity** (*Potapov et al. 2018*):
  Scores multi-part Type IIS (`BsaI`, `BsmBI`, `BbsI`, `SapI`) junction sets against empirical T4 DNA ligase mismatch ligation frequencies and identifies synonymous silent codon substitutions to domesticate internal cut sites.

---

## Tool #3 — `dwc` (Biodiversity & Darwin Core QA)

- **Modules**: `dwc_core`, `dwc_taxonomy`, `dwc_export`
- **Spatial & Taxonomic Scrubbing** (*Zizka et al. 2019* `CoordinateCleaner` rules):
  Normalizes DMS/DDM/DD coordinates, flags `(0, 0)` Null Island, swapped latitude/longitude, out-of-range coordinates (`|lat| > 90`, `|lon| > 180`), and future/pre-1500 `eventDate` ISO-8601 entries; resolves binomial authorship and synonyms into GBIF IPT-ready Darwin Core Archives (`occurrence.csv`, `meta.xml`, `eml.xml`, GeoJSON, SVG).

---

## Tool #4 — `phylo` (Phylogenetics & Tree Comparison)

- **Modules**: `phylo_core`, `phylo_render`
- **Arena Tree Operations**:
  Parses Newick, NHX, and NEXUS (`TRANSLATE` + `BEGIN TREES`) trees; performs outgroup and longest-patristic-diameter midpoint rerooting; classifies target taxon sets as monophyletic, paraphyletic, or polyphyletic via MRCA clade coverage.
- **Robinson-Foulds Distance** (*Robinson & Foulds 1981*):
  Unrooted bipartition split sets `S(A)` and `S(B)` yield symmetric topological distance `RF = |S(A) Δ S(B)|` plus normalized `RF` and branch-length `L1` score.

---

## Tool #5 — `epi` (Epidemiology & Study-Design Power)

- **Modules**: `epi_models`, `epi_power`
- **Cori Time-Varying Reproduction Number `R_t`** (*Cori et al. 2013*):
  Discretizes a shifted Gamma serial interval `w_s` (`mean`, `sd`) and computes exact conjugate Gamma posterior `(a_t, b_t)` over sliding window `tau` (`a_post = a_prior + sum I_u`, `b_post = b_prior + sum Lambda_u`, `Lambda_u = sum_{s=1}^u I_{u-s} w_s`) with 95% credible intervals via Newton-refined regularized incomplete Gamma quantiles.
- **Compartmental `SIR` / `SEIR` ODEs**:
  4th-order Runge-Kutta (`RK4`) integration supporting initial immunity and step NPI transmission reduction.
- **Diagnostic Bayes & Study Power**:
  `LR+ = Se / (1 - Sp)`, `LR- = (1 - Se) / Sp`, Rogan-Gladen (1978) true prevalence `(p_obs + Sp - 1) / (Se + Sp - 1)`, Fleiss et al. (1980) continuity-corrected cohort/RCT power, and Jackson & Nelson (2013) Test-Negative Design (`VE = 1 - OR`) sample size.

---

## Tool #6 — `lexi` (Taxonomic & Anatomical Etymology)

- **Modules**: `lexi_roots`, `lexi_engine`
- **Viterbi-Style Morpheme Lattice**:
  Dynamic programming over Classical Greek, Latin, and ICZN/ICNafp taxonomic suffixes minimizes segmentation penalty across stems, connecting vowels, and affixes to reconstruct cognate trees and Anki flashcard decks.

---

## Tool #7 — `forge3d` (Biomolecule & Morphology 3D Mesh Engine)

- **Modules**: `forge3d_mesh`, `forge3d_bio`, `forge3d_gui`
- **PDB Solvent Gaussian Surface & Backbone Tube**:
  Assigns Bondi (1964) van der Waals radii (`r_vdW`, Å) by element and splats Gaussian density `sum exp(-0.80 * d^2 / (1.25 * r_vdW)^2)` on a regular 3D scalar grid.
- **Gielis Biological Superformula** (*Gielis 2003*):
  ```text
  r(φ) = ( |cos(m * φ / 4) / a|^n2 + |sin(m * φ / 4) / b|^n3 )^(-1 / n1)
  ```
  Modulates silica diatom frustule symmetry lobes `m` with a domed valve profile and hexagonal areolae perforations.
- **Surface Nets Extraction, Taubin Smoothing & Printability Audit**:
  Extracts closed isosurfaces via Surface Nets, applies shrinkage-free Taubin (1995) two-pass Laplacian smoothing (`λ = +0.50`, `μ = -0.53`), normalizes to `Z_min = 0` on the print bed, and computes:
  - Exact divergence signed volume: `V = (1 / 6) * sum(v0 · (v1 × v2))` (`mm³`)
  - Half-edge 2-manifold watertightness (`count(u -> v) == 1` and `count(v -> u) == 1` for every directed edge)
  - Topological Euler characteristic: `χ = V - E + F` (`= 2` for a single genus-0 closed shell; `2 - 2g` for genus `g` porous shells)
  - Estimated FDM PLA (`1.24 g/cm³`, 20% infill + 0.8 mm shell) and SLA resin (`1.12 g/cm³`, solid) mass
- **3D Viewport Rendering**:
  Supports perspective projection (`cam_dist / (cam_dist - z)`) and orthographic CAD mode, Blinn-Phong specular (`(N · H)^36`) + Fresnel rim + secondary fill + depth cueing, a world-space `Z_min` build-plate grid, and an orientation-locked XYZ axis triad.

---

## Tool #8 — `morpho` (Colony Counting, Leaf Herbivory & 2D GPA)

- **Modules**: `morpho_vision`, `morpho_gpa`
- **Automated CFU & Leaf Herbivory**:
  Otsu (1979) variance-maximizing thresholding, `(3, 4)` Chamfer Euclidean distance transform local-peak watershed splitting for touching colonies, Poisson MLE `CFU/mL`, and Woebbecke et al. (1995) Excess Green (`ExG = 2G - R - B`) + border flood-fill enclosed hole area estimation.
- **2D Generalized Procrustes Analysis (`GPA`)** (*Rohlf & Slice 1990*):
  Centers each specimen's landmark configuration at the origin, normalizes by isometric Centroid Size `CS = sqrt(sum ||x_i - x_bar||^2)`, aligns via Procrustes optimal planar rotation `atan2(sum(x' y - y' x), sum(x' x + y' y))` without reflection, and runs cyclic Jacobi PCA on tangent shape residuals.

---

## Tool #9 — `plate` (96/384-Well Assay Design & 4PL Curves)

- **Modules**: `plate_layout`, `plate_curve`
- **Microplate QC & 4PL Regression**:
  Generates edge-moat buffered (`Lundholt et al. 2003`) `SplitMix64`-randomized layouts; computes Zhang et al. (1999) `Z' = 1 - 3*(σ_pos + σ_neg) / |μ_pos - μ_neg|` and `SSMD`; fits 4-Parameter Logistic (`4PL`) dose-response models (`Sebaugh 2011`) via damped Levenberg-Marquardt with analytical Jacobians:
  ```text
  y(x) = Bottom + (Top - Bottom) / (1 + (IC50 / x)^Hill)
  ```

---

## Tool #10 — `ecosound` (Bioacoustics STFT & Ecoacoustic Indices)

- **Modules**: `ecosound_dsp`, `ecosound_slicer`, `ecosound_gui`
- **STFT & Zero-Phase Filtering**:
  Cooley-Tukey radix-2 FFT with Hann, Hamming, or 4-term Blackman-Harris windows (`dBFS = 20 * log10(mag)`); zero-phase forward-backward 2nd-order Butterworth bandpass filtering (`filtfilt` style).
- **Soundscape Ecology Indices**:
  - **ACI** (*Pieretti, Farina & Morri 2011*): Sum across biophony bins and 5 s frame clumps of `sum(|I_k - I_{k+1}|) / sum(I_k)`, plus duration-normalized `ACI/min`. Range `>= 0` (unitless).
  - **NDSI** (*Kasten et al. 2012*): `(PSD_bio - PSD_anthro) / (PSD_bio + PSD_anthro)` in `[-1.0, +1.0]` (unitless).
  - **Acoustic Entropy `H`** (*Sueur et al. 2008*): `H = H_t * H_f` in `[0.0, 1.0]` (unitless), where `H_t` is normalized Shannon entropy of the 5 ms RMS temporal envelope and `H_f` is normalized Shannon entropy of the mean amplitude spectrum.
  - **ADI & AEI** (*Villanueva-Rivera et al. 2011*): Activity fraction `p_b` of STFT cells `> -50 dBFS` across ten 1 kHz sub-bands up to `min(Nyquist, 10 kHz)`: `ADI = -sum(p_b' * ln(p_b'))` (`>= 0` nats) and `AEI` = Gini coefficient across `p_b` in `[0.0, 1.0]`.
- **Call Slicer & Raven Pro Export**:
  Adaptive median-floor SNR detector with gap merging, `-14 dB` spectral bounding-box tightening, FM slope linear regression (`Hz/s`), 3 ms cosine-faded `.wav` clip slicing, and 13-column Raven Pro 1.6 `.tsv` selection tables.

---

## Tool #11 — `infoevo` (Information-Theoretic Evolution)

- **Modules**: `infoevo_core`, `infoevo_svg`
- **Sequence Conservation & Complexity**:
  Schneider & Stephens (1990) `R_seq(l) = log2(A) - H_obs(l) - e_n` (bits) with Basharin small-sample correction; Adami, Ofria & Collier (2000) physical complexity `C = sum R_seq(l)`.
- **Epistasis & Kimura Selection Bound**:
  Dunn, Wahl & Gloor (2008) Average Product Corrected Mutual Information `MIp(i, j) = MI(i, j) - (MI_bar(i) * MI_bar(j)) / MI_bar_global`; clade divergence `D_JS` / `D_KL`; Kimura (1961) genetic information gain `I_t = D_KL(p_t || p_0) <= L_sub,t`.

---

## Tool #12 — `qpcr` (MIQE RT-qPCR & Melt Curves)

- **Modules**: `qpcr_core`, `qpcr_melt`
- **Relative Quantification**:
  Efficiency-calibrated `E = 10^(-1/slope)`, Livak & Schmittgen (2001) `2^-ΔΔCq`, Pfaffl (2001) single-reference ratio, Vandesompele et al. (2002) `geNorm` pairwise log-ratio SD stability `M`, and Hellemans et al. (2007) `qBase` geometric-mean multi-reference Normalized Relative Quantities (`NRQ`) with Welch `t`-tests and Bustin et al. (2009) MIQE checks.
- **SYBR Melt Curves** (*Ririe, Rasmussen & Wittwer 1997*):
  Savitzky-Golay smoothed `-dF/dT` peak calling and low-Tm primer-dimer flagging.

---

## Tool #13 — `kinetics` (Enzyme Inhibition & Morrison Tight Binding)

- **Modules**: `kinetics_core`, `kinetics_svg`
- **Single-Substrate & Global Reversible Inhibition Models**:
  Filters non-finite (`NaN`/`Inf`) observations prior to parameter initialization and log-space damped Levenberg-Marquardt nonlinear least squares. Compares Michaelis-Menten, Allosteric Hill, and Haldane Substrate Inhibition (`[I] = 0`) and global 2D `(Competitive, Uncompetitive, PureNoncompetitive, MixedType)` inhibition via Burnham & Anderson (2004) small-sample `AICc` and Akaike weights `w_i`.
- **Morrison Tight-Binding Quadratic** (*Morrison 1969; Williams & Morrison 1979*):
  ```text
  v_i = v_0 * (1 - ((E_T + I_T + Ki_app) - sqrt((E_T + I_T + Ki_app)^2 - 4 * E_T * I_T)) / (2 * E_T))
  ```
  Jointly estimates active-site enzyme concentration `E_T` (`nM`) and apparent `Ki_app` (`nM`), plus competitive intrinsic `Ki = Ki_app / (1 + [S] / Km)` and Cheng & Prusoff (1973) `IC50 -> Ki` conversions.

---

## Tool #14 — `edna` (Metabarcoding FASTQ & Biodiversity)

- **Modules**: `edna_fastq`, `edna_diversity`
- **Amplicon Pipeline**:
  Hamming-tolerant dual-barcode demultiplexing, Edgar & Flyvbjerg (2015) expected-error (`EE = sum 10^(-Q/10)`) filtering and `UNOISE` 1-nt abundance-ratio denoising (`alpha = 2.0`), de novo 2-parent chimera check, and negative extraction-blank contaminant auditing.
- **Diversity Estimators**:
  Chao (1984) `Chao1`, Hurlbert (1971) exact hypergeometric rarefaction `E[S_m] = sum_i [1 - C(N - N_i, m) / C(N, m)]` and `PIE`, Chao et al. (2014) Hill numbers (`^0D`, `^1D = exp(H')`, `^2D = 1 / sum p_i^2`) and Good-Turing coverage, plus Bray-Curtis and Jaccard beta dissimilarity.

---

## Tool #15 — `wildlife` (Site Occupancy & Capture-Recapture)

- **Modules**: `wildlife_occupancy`, `wildlife_markrecap`
- **Imperfect-Detection Occupancy** (*MacKenzie et al. 2002, 2005*):
  Nelder-Mead binomial/detection-history MLE for `Psi(.)p(.)`, `Psi(cov)p(.)`, and `Psi(.)p(cov)` ranked by `AICc`, empirical Bayes conditional occupancy, and survey effort planning (`K*`, `S*`).
- **Closed & Open Capture-Recapture**:
  Chapman (1951)/Seber (1970) 2-sample estimator, Otis et al. (1978) `M0` / Darroch `Mt` profile likelihoods, Schnabel (1938), Chao (1987) `Mh` heterogeneity estimator, and Lebreton et al. (1992) Cormack-Jolly-Seber (`CJS`) `m`-array open survival (`phi`) and recapture (`p`).

---

## Tool #16 — `popgen` (Population Genetics & Haplotype Networks)

- **Modules**: `popgen_core`, `popgen_network`
- **FASTA Diversity & VCF Variance Components**:
  Tajima (1989) `pi`, Watterson's `theta_W`, and Tajima's `D`; pairwise `d_XY`, `d_a`, and Hudson `F_ST`; diploid VCF Weir & Cockerham (1984) `F_IS`, `F_ST`, `F_IT`; Wigginton et al. (2005) exact Levene-Haldane `HWE` p-values; Bandelt et al. (1999) Median-Joining haplotype network construction with Fruchterman-Reingold planar layout.

---

## Tool #17 — `crispr` (Guide Design & CFD Off-Target Scoring)

- **Modules**: `crispr_core`, `crispr_offtarget`
- **Multi-Nuclease & Base-Editing Audit**:
  PAM scanning for `SpCas9` (`3'-NGG`), `SaCas9` (`3'-NNGRRT`), and `LbCas12a` (`5'-TTTV`), Doench et al. (2016) Rule Set 2 on-target scoring and position/mismatch `CFD` off-target cleavage matrix, Sugimoto RNA:DNA R-loop `Tm`, and Komor et al. (2016) `CBE` (`C -> T`) / Gaudelli et al. (2017) `ABE` (`A -> G`) bystander windows.

---

## Tool #18 — `flow` (ISAC FCS 3.1 Cytometry & Logicle Gating)

- **Modules**: `flow_fcs`, `flow_gate`
- **Binary I/O, Spillover & Biexponential Transform**:
  ISAC `FCS 3.0/3.1` binary parser/writer (`F`/`D`/`I` data types, big/little-endian) with CRC/offset validation (*Spidlen et al. 2010*); exact `$SPILLOVER` matrix inversion (`S^-1`) via Gaussian elimination; Parks, Roederer & Moore (2006) / Moore & Parks (2012) `Logicle` biexponential inverse (`S(y; T, W, M, A)` solved via monotonic Newton-Raphson); ray-casting 2D polygon and 4-quadrant gating (`% Parent`, `MFI`, robust `CV%`).

---

## Tool #19 — `survival` (Kaplan-Meier, Cox PH & Acute Ecotoxicology)

- **Modules**: `survival_core`, `survival_tox`
- **Time-to-Event & Dose-Mortality Models**:
  Kaplan & Meier (1958) product-limit survival `S(t)` with Greenwood (1926) SE and log-log 95% CIs; Mantel (1966) log-rank `χ²` statistic; Cox (1972) Breslow partial-likelihood Newton-Raphson regression (`HR` and Wald 95% CI); Finney (1971) Probit and Binomial Logit `LC10` / `LC50` / `LC90` MLE with Fieller 95% fiducial limits.

---

## Tool #20 — `flux` (Leaf FvCB Photosynthesis & Odum Stream Metabolism)

- **Modules**: `flux_leaf`, `flux_stream`
- **C3 Leaf FvCB Photosynthesis** (*Farquhar, von Caemmerer & Berry 1980; Bernacchi et al. 2001, 2003; Medlyn et al. 2002; Harley et al. 1992; von Caemmerer 2000*):
  Net CO2 assimilation `A_n` (`umol m^-2 s^-1`) is evaluated across Rubisco-limited (`Wc`), RuBP-regeneration-limited (`Wj`), and triose-phosphate-utilization-limited (`Wp`) rates:
  ```text
  T_K = T_leaf_C + 273.15   (guarded cleanly for T_leaf_C <= -273.15 °C and non-finite inputs)
  Arrhenius(k25, Ha) = k25 * exp(Ha * (T_K - 298.15) / (298.15 * R * T_K)),   R = 0.0083145 kJ/(mol*K)
  Peaked(k25, Ha, Hd, dS) = Arrhenius(k25, Ha) * (1 + exp((298.15*dS - Hd)/(298.15*R)))
                                               / (1 + exp((T_K*dS - Hd)/(T_K*R)))
  ```
  - Rubisco kinetics (`Bernacchi et al. 2001`): `Kc(25) = 404.9 ubar (Ha = 79.43)`, `Ko(25) = 278.4 mbar (Ha = 36.38)`, `Γ*(25) = 42.75 ubar (Ha = 37.83)`, `[O2] = 210 mbar`, `Rd Ha = 46.39`.
  - Peaked high-temperature deactivation (`Medlyn 2002 / Bernacchi 2003 / Harley 1992`): `Vcmax (Ha = 65.33, Hd = 200.0, dS = 0.650)`, `Jmax (Ha = 43.54, Hd = 200.0, dS = 0.650)`, `TPU (Ha = 53.10, Hd = 201.8, dS = 0.650)` — exact at `25 °C` (`deactivation factor = 1.0`), physiologically declining above `35–45 °C`.
  - PSII electron transport with `(1 - f)` (`f = 0.15`) spectral-quality factor:
    ```text
    I_abs = 0.85 * 0.5 * (1 - 0.15) * PPFD
    J = ((I_abs + Jmax) - sqrt((I_abs + Jmax)^2 - 4 * 0.70 * I_abs * Jmax)) / (2 * 0.70)
    Wc = Vcmax * Ci / (Ci + Kc * (1 + O2 / Ko))
    Wj = J * Ci / (4 * Ci + 8 * Γ*)
    Wp = 3 * TPU * Ci / (Ci - Γ*)      (for Ci > Γ* + 1.0)
    A_n = min(Wc, Wj, Wp) * (1 - Γ* / Ci) - Rd
    ```
  Parameters `(Vcmax25, Jmax25, TPU25, Rd25)` are fitted via multi-start bounded Nelder-Mead SSE minimization.
- **Odum Open-Channel Diel Stream Metabolism** (*Odum 1956; Benson & Krause 1984; Wanninkhof 1992*):
  - **Dissolved-Oxygen Saturation** (`mg L^-1`, Benson & Krause 1984 at `1 atm`) with optional barometric/elevation `P / P0` pressure correction (`P(h) = (1 - 2.25577e-5 * h)^5.25588` atm).
  - **Reaeration & Respiration Scaling**:
    `K(T) = K600 * (Sc_O2(T) / 600)^(-0.5)` via Wanninkhof (1992) freshwater `O2` Schmidt polynomial `Sc_O2(T)`. Optional Ecosystem Respiration temperature scaling (`ErTempScaling::Constant`, `Q10`: `ER_ref * Q10^((T - T_ref)/10)`, or `Arrhenius`: `ER_ref * exp((Ea / R) * (1/T_ref_K - 1/T_K))`).
  - **Semi-Implicit ODE Step & Reporting**:
    ```text
    C_{i} = (C_{i-1} + dt * ((GPP * (PAR_mid / sum(PAR*dt)) - ER(T_mid)) / z_mid + K_T * C_sat,mid)) / (1 + dt * K_T)
    NEP = GPP - ER_daily   (g O2 m^-2 d^-1)
    ```
    Dual-panel SVG rendering (`render_flux_svg` using shared `svg_util` primitives) automatically scales Y-axis bounds to observed/modeled/saturation `DO` ranges and overlays `Wc`, `Wj`, `Wp`, and `A_net` curves with ticks, gridlines, and legends.
