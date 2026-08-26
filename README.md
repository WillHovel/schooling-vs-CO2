# Kinemetrix — Fish Swimming Kinematics Toolkit

A MATLAB toolkit for analysing fish swimming kinematics from digitised landmark data (DeepLabCut or similar). Given a CSV of tracked body points, the pipeline aligns each midline to the swimming axis, then extracts beat frequency, amplitude envelopes, propulsive wavelength, curvature, speed, stride length, Strouhal number, fin/girdle kinematics, stance–swing duty factor, pectoral phase, and school-level metrics (polarization, angle-to-flow, inter-individual distance).

Built on the midline-alignment method of **Castro-Santos & Goerig (2017)**, with an auto-orientation extension and numerous accuracy fixes (see the `CHANGE NOTE` comments throughout the source). Every metric below is documented with its formula and the exact code block that computes it.

> **Metric reference at a glance** — the sections below catalogue every number the program outputs:
> [1. Midline kinematics](#1-midline-kinematics) · [2. Body & gait metrics](#2-body--gait-metrics) · [3. Fin kinematics](#3-fin-kinematics) · [4. Girdle kinematics](#4-girdle-kinematics) · [5. Stance–swing (duty factor)](#5-stanceswing-duty-factor) · [6. Pectoral phase](#6-pectoral-phase) · [7. School metrics](#7-school-metrics) · [8. Per-cycle frequencies](#8-per-cycle-frequencies) · [Preprocessing](#preprocessing) · [CSV export](#csv-export)

---

## Pipeline at a glance

```
CSV file
   │
   ├─ detectFormat()                FishKinematicsApp.m   — classify A/B/C/D/E or DLC multi-animal
   │
   ├─ load_*()                      — points into struct array  fp  (RAW camera frame)
   │     ├─ filter_dlc_jumps()      (optional) NaN-out DLC tracking jumps
   │     ├─ remap_axes()            (optional, 3D)  swap/negate X/Y/Z columns
   │     └─ [school metrics stash]  fp saved RAW before any transform
   │
   ├─ transform_fish()              — per-frame rotate + translate + BL-normalize
   │
   ├─ compute_kinematics()          — FFT analysis  →  kine
   ├─ compute_body_extended()       — speed / stride / Strouhal / angles  →  ext
   │
   ├─ (Fin tab)    compute_fin_kinematics()  → compute_girdle_kinematics()  → compute_stance_swing()
   ├─ (School tab) compute_polarization() / compute_angle_to_flow() / compute_distance_between_individuals()
   │
   └─ collect_*_row()               — stage one flat row per metric set for CSV export
```

### The GUI tabs

| Tab | Purpose | Metrics computed by |
|---|---|---|
| **Kinematics** (main) | Load a trial, run midline analysis, browse per-fish results | `transform_fish`, `compute_kinematics`, `compute_body_extended` |
| **Fin Analysis (3D)** | Fin root→tip Euler angles, fin beat frequency, girdle protraction, stance/swing | `compute_fin_kinematics`, `compute_girdle_kinematics`, `compute_stance_swing` |
| **Pectoral Phase (2D)** | Dual-camera right-vs-left pectoral phase classification | `compute_pect_phase` (in `load_fish_points_named.m`) |
| **Batch Processing** | Run the midline pipeline over a folder of CSVs, per-file flow speed from filename tokens | same as main tab, headless |
| **School Metrics** | Multi-fish polarization, angle-to-flow, nearest-neighbour distance | the three `compute_*` group functions |

### Data formats (auto-detected by `detectFormat()`, FishKinematicsApp.m ≈ line 650)

| Format | Detection | Files |
|---|---|---|
| **A — DLC multi-animal** | native DLC export, row 2 = `individuals` (probed with `readcell` **before** `readtable`, which would coerce the text row to NaN) | 4-row header (`scorer`/`individuals`/`bodyparts`/`coords`), N animals × M bodyparts |
| **A — DLC single** | `Fish1_P1_x` style indexed columns | one animal per file, multi-animal columns also loadable |
| **B — Named** | `snout_X`, `peduncle_Y`, … columns | landmark-named CSVs |
| **C — Numbered** | `pt1_X … pt12_X` (12 anatomical points) | numbered-point CSVs |
| **D — Dual-camera** | `pt1_cam1_X …` (cam1 = lateral, cam2 = ventral) | pectoral-phase trials |
| **E — CURVES** | pre-transformed `X/Y` in BL | files already in body frame (`pre_transformed = true` skips `transform_fish`) |

### Units & conventions

- **BL** = body length. By default the per-frame distance between the two endpoint landmarks (in raw units); override with a known body length via `bl_override` / the GUI's body-length field (recommended for Format C, whose endpoints are fin bases, not snout/tail).
- **Speeds** in BL/s. Angles in degrees. Frequencies in Hz.
- **Flow sign convention** (`compute_body_extended`): *positive* flow = swimming **against/upstream** (through-water speed = ground speed + flow); *negative* = with the current.
- **Headings** follow `atan2d` (counter-clockwise positive); pectoral "phase shift" is wrapped `mod(…, 360)`.
- NaN everywhere means *not computable* — the code deliberately returns NaN rather than a fabricated value for degenerate input (all-NaN, zero variance, too few points). See the `CHANGE NOTE` comments for the specific bugs this guards against.

---

# Metric reference

Every metric the program outputs, grouped by the compute function that produces it. Line numbers are approximate (as of 2026-08).

## 1. Midline kinematics

**File:** `compute_kinematics.m` — **GUI:** Kinematics tab.
Runs on `transform_fish` output. Interpolates each frame's n-point midline to `N_OUT = 200` stations `s_norm = 0…1` via zero-padded FFT interpolation (§1, lines 43–107), after restoring the world-frame lateral displacement `Yc_raw` (adds back the per-frame fitted midline's rotation so the rigid-body recoil is later removable by spatial-mean subtraction).

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `head_TBF` | Hz | Dominant temporal frequency of the cleaned world-frame lateral displacement at the head station. `Yc_raw(:,1)` is first **linearly detrended** to remove the rigid-body translation ramp of the swimming path (a mean-removal alone leaves a sawtooth whose 1/f spectrum buries the beat), then `dominant_freq` takes the peak of `power = (2/N)·|FFT|²` above `min_freq`. | §4 lines 127–185; helper `dominant_freq` 425–468 |
| `tail_TBF` | Hz | Same FFT, but on the **camera-space tail Y relative to the per-frame centroid of all tracked points** — removes common-mode sway without inheriting the per-frame line fit's wobble (raw camera Y is sway-dominated; cleaned Y can be line-fit-artifact-dominated; centroid-relative measured the true beat in validated swim *and* walk trials). Falls back to `Yc_raw(:,end)` for pre-transformed data. | §4 lines 155–185 |
| `spline_freq_Hz` | Hz | Comparison value: dominant temporal frequency at **each** of the 200 interpolated stations, summarized as the **median over stations with amplitude ≥ 10 % of the max station amplitude** (quiet near-head stations don't skew it). | §4b lines 198–221 |
| `headZ_TBF` / `tailZ_TBF` | Hz | Dorso-ventral beat (3-D only). Head measured at station 2 — station 1 is identically zero after BL normalization. | §4 lines 187–196 |
| `amp_mean` / `amp_std` | BL | Per-station amplitude envelope (1×200): temporal mean and SD of the **absolute** lateral displacement `|Y_interp|` at each station — i.e. the mean half-amplitude of the traveling wave along the body. | §2 lines 109–116; `amplitude_stats` 355–375 |
| `headAmp` / `tailAmp` | BL | Envelope values averaged over the head window (s ≤ 0.05) and tail window (s ≥ 0.95). | `amplitude_stats` 366–367 |
| `headTailAmpRatio` | – | `headAmp / tailAmp`. | line 115 |
| `minAmp`/`maxAmp` + `minAmpLoc`/`maxAmpLoc` | BL / s | Extrema of the amplitude envelope and their body positions. | `amplitude_stats` 369–374 |
| `ampZ_mean` / `ampZ_std`, `headAmpZ` … | BL | Dorso-ventral equivalents of the whole amplitude set (3-D only). | §3 lines 117–125 |
| `wavelength` | BL | Propulsive wavelength via a **two-stage complex-amplitude traveling-wave fit**: the per-station complex amplitudes at the beat frequency are corrected for the known "s-linear contamination" (the transform's line-fit wobble, §`slinear_contamination` 736–759), then a complex exponential model is fit to the offset-corrected amplitudes and the wavelength read from the phase gradient. This replaced an envelope-FFT method that always returned ≈1.005 BL regardless of input. | §5 lines 223–238; `spatial_wavelength` 469–718 |
| `wave_speed_BL_s` | BL/s | `wavelength × tail_TBF` (TBF-based by design; head/tail/spline frequencies are exported alongside for comparison). | line 238 |
| `curv_mean` / `curv_std` | 1/BL | Per-station mean/SD curvature profile (1×200): three-point geometric curvature via the circumcircle radius, `κ = 1/R`, `R = ABC/(4·area)` of the triangle formed by stations k−lag, k, k+lag with `lag = N_OUT/40` (≈2.5 % BL each side), averaged across frames. | §6 lines 240–244; `curvature_stats` 378–422 |
| `maxCurv` / `maxCurvLoc` | 1/BL, s | Peak of the mean curvature profile and its body position. | `curvature_stats` 416–421 |
| `curv3d_mean` … `maxCurv3DLoc` | 1/BL, s | 3-D curvature equivalents (same circumcircle method with Z included; 3-D only). | §7 lines 246–253 |

**Parameters:** `fps` (frame rate), `min_freq` (Hz floor for every FFT peak — set it from a visual beat count for the trial; too high a floor makes the FFT lock onto a harmonic).

## 2. Body & gait metrics

**File:** `compute_body_extended.m` — **GUI:** Kinematics tab (runs after `compute_kinematics`; needs `kine` for TBF). Also takes an optional **flow speed** (BL/s) and, for 3-D data, an optional left/right landmark pair for roll.

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `body_angle_deg` + mean/SD/range | deg | Per-frame heading: slope of the least-squares line through the **middle points** (same fit `transform_fish` uses), `atan2d`-style via `atand(slope)`, then **unwrapped mod 180°** so each step is forced into (−90°, 90°] — a fish turning past ±90° no longer produces a fake 180° jump. Mean/SD/range are relative to the first valid frame. Left NaN for pre-transformed (CURVES) data, where heading no longer exists. | §1 lines 107–171 |
| `angular_velocity_deg_s` | deg/s | Frame-to-frame `diff(body_angle_unwrapped) × fps`. | line 162 |
| `speed_BL_s` + mean/SD/peak | BL/s | Forward speed of the body **centroid**: frame-to-frame centroid displacement (raw units) ÷ this frame's body length × fps. Falls back to the previous frame's BL when the current one is missing; pre-transformed data uses BL = 1. | §2 lines 173–217 |
| `speed_through_water_BL_s` + mean/SD/peak | BL/s | Ground speed corrected for flume flow: `+flow` when swimming against, `abs(ground + flow)` when with the current. Assumes the path is aligned with the flow axis (standard station-holding assumption). | §2b lines 219–246 |
| `stride_length_BL` | BL | Mean distance travelled per tail-beat: `mean_speed / tail_TBF`. | §3 lines 248–256 |
| `tail_amp_pp_BL` | BL | Peak-to-peak tail lateral excursion. Measured as the tail's **perpendicular distance from a smoothed body axis** (line orientation + mid-point centroid averaged over ≈1.5 beat periods, with the orientation smoothed as `exp(i·2α)` to survive ±180° ambiguity), then linearly detrended and `A_pp = 2·√2·std(yt)`. The smoothing matters: the per-frame fitted axis itself oscillates with the beat and inflates the apparent excursion ≈50 %. | §3b lines 266–326 |
| `strouhal` | – | `St = tail_TBF × tail_amp_pp_BL / U`, with U = mean **through-water** speed when flow is set (physically correct reference in a flume), else mean ground speed. | §3b lines 328–333 |
| `head_pitch_deg` + mean/SD/range | deg | Elevation of the head-to-next-point vector vs horizontal, `atand(Δz/horiz_dist)`, + = head raised (3-D only). | §4 lines 338–363 |
| `head_Z_raw` | BL | Head Z relative to its first valid frame's Z (3-D only). | §4 lines 354–355 |
| `roll_deg` + mean/SD/range | deg | Body roll from a left/right landmark pair: `atan2d(zR−zL, yR−yL)` (3-D + supplied pair only; `roll_available` flags whether a pair was found). | §5 lines 365–395 |

## 3. Fin kinematics

**File:** `compute_fin_kinematics.m` — **GUI:** Fin Analysis (3D) tab.
Loads two landmarks (root, tip — each may be a compound `'Rpectbase+Lpectbase'` average) and computes the fin vector tip−root.

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `fin_length`, `mean_length`, `std_length` | raw | `‖tip − root‖` per frame, mean/SD over valid frames. | lines 100–101, 169–171 |
| `yaw` / `pitch` / `roll` (+ mean/SD/range) | deg | Fin-vector Euler angles: `yaw = atan2d(vy,vx)`, `pitch = atan2d(−vz,√(vx²+vy²))`, `roll = atan2d(vz,vy)`; means/SDs/ranges over valid frames. | lines 103–106, 192–197 |
| `step_dist` / `cum_dist` / `total_dist` | raw | Tip displacement per frame / cumulative / total tip travel over the trial. | lines 108–119 |
| `tip_speed` + mean/SD/peak | raw/s | `step_dist × fps`. (Frame 1 excluded from stats — it has no previous frame, so its speed is 0 by construction.) | lines 119, 150–155 |
| `d_yaw` / `d_pitch` / `d_roll` | deg/s | Angular rates: `diff(angle) × fps`. | lines 121–124 |
| `ang_vel` + mean/SD/peak | deg/s | Total angular speed `√(d_yaw² + d_pitch² + d_roll²)`. | line 125, 202–204 |
| `fin_freq_Hz` | Hz | **Fin-beat frequency** from the temporal FFT of the yaw signal (dominant peak above `min_freq`; NaN-safe: degenerate input → NaN). | lines 127–134 |
| `fin_freq_pitch_Hz` | Hz | Same from the pitch signal — use whichever matches your visual beat count (flap- vs sweep-dominated fins differ). | lines 132–134 |
| `stride_duration_s` | s | `1 / fin_freq_Hz`. | lines 137–144 |
| `stride_length_BL` | BL | `body_speed_BL_s / fin_freq_Hz` — forward distance per fin-beat (NaN unless you pass `ext.mean_speed_BL_s` in). | lines 141–143 |
| `n_valid`, `pct_valid` | –, % | Valid frames / % of the *tracked* window (padding rows excluded). | lines 208–213 |

## 4. Girdle kinematics

**File:** `compute_girdle_kinematics.m` — **GUI:** Fin Analysis (3D) tab.
Expresses a girdle-base landmark in the **same body-relative frame as the midline** (via `apply_body_transform` + the midline's `transform_params`), capturing girdle rotation about the body — a signal genuinely separate from fin yaw/pitch.

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `X` / `Y` / `Z` | BL | Girdle point in body-relative coordinates (X fore-aft: protraction → smaller X; Y lateral splay). | projection lines 114–115 |
| `protraction_range_BL` | BL | Fore-aft excursion range of X: `range(X)`. | lines 129–130 |
| `lateral_range_BL` | BL | Side-to-side excursion range of Y. | lines 129–130 |
| `girdle_freq_Hz` | Hz | Dominant FFT frequency of the fore-aft (X) signal — the girdle's own protraction/retraction beat. | lines 132–136 |
| `n_valid` / `pct_valid` | –, % | Valid frames / % of tracked window. | lines 117–121 |

## 5. Stance–swing (duty factor)

**File:** `compute_stance_swing.m` — **GUI:** Fin Analysis (3D) tab.
> **Important:** "contact" here is a **kinematic proxy**, not a force measurement — a frame is STANCE when the fin tip's speed < 25 % (default) of the trial's own 95th-percentile tip speed, SWING otherwise. If you have real contact-time data, treat that as ground truth and this as a secondary check.

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `threshold_used` | raw/s | `0.25 × p95(tip_speed)` — the actual speed threshold. | lines 92–96 |
| `is_stance` | – | Per-frame logical classification (zero-by-construction speeds at tracking gaps count as *unknown*, never stance). | lines 95–96 |
| `mean_contact_time_s` / `mean_swing_time_s` | s | Mean bout duration per phase, after bridging data gaps ≤ 2 frames (default) and discarding bouts < 3 frames (noise filter). | lines 98–134 |
| `duty_factor` | – | `mean_contact / (mean_contact + mean_swing)`. | line 133 |
| `n_cycles` | – | Complete stance→swing→stance cycles found. | line 136 |
| `pct_valid_consecutive` | % | Frames with usable frame-to-frame tip speed. **Below ≈70 % this estimate is likely unreliable** (tracking gaps fragment real bouts) — check this before trusting `duty_factor`. | lines 76–84 |

## 6. Pectoral phase

**File:** `compute_pect_phase` (local function in `load_fish_points_named.m`, lines 227–320) — **GUI:** Pectoral Phase (2D) tab (display: `displayPectPhase`, FishKinematicsApp.m line 2139).
For Format D dual-camera files: picks the best-tracked camera per fin, then compares the phase of the right (pt2) vs left (pt12) pectoral tip signals.

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `phase_shift_deg` | deg | Phase difference between the fins at the dominant FFT bin of the right-fin signal: `mod(phase12 − phase2, 360)`, where each phase is `angle(FFT(signal))` at that bin. | lines 279–292 |
| `phase_right_deg` / `phase_left_deg` | deg | Individual fin phases at the dominant bin. | lines 289–290 |
| `classification` | – | `In-phase` (\|shift\| ≤ 45° or ≥ 315°), `Antiphase` (135–225°), `Intermediate` (otherwise). | lines 294–300 |
| `peak_lag_frames` | frames | Lag (in frames) of the peak of the normalized cross-correlation between the two fin signals (positive = left lags right). | lines 302–304 |
| `dom_freq_norm` | cycles/frame | Dominant (beat) frequency of the right-fin signal in normalized units. | line 310 |
| `n_valid` | – | Frames where both fins were tracked (needs ≥ 20, else `N/A`). | lines 266–274 |

## 7. School metrics

**School Metrics tab** (GUI) — computed from the **RAW, untransformed** struct array (stashed as `app.school_fp` before `transform_fish`), because each fish's private post-transform frame would destroy the shared spatial relationships these metrics need. Loads native DLC multi-animal exports directly (Format A), or any multi-fish file via the tab's own Browse.

### 7a. Polarization — `compute_polarization.m`

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `polarization` | 0–1 | Per-frame school alignment: `P(f) = |mean_k exp(i·heading_k(f))|` — the length of the average unit heading vector (headings = peduncle→snout vectors, `atan2d`). 0 = headings random/opposed, 1 = all facing the same way; NaN when < 2 fish tracked that frame. | lines 99–107 |
| `mean_polarization` / `std_polarization` | 0–1 | Mean/SD over valid frames. | lines 109–115 |
| `heading_deg` / `heading_to_flow_deg` | deg | Per-fish per-frame headings, raw and re-centered on the flow axis. | lines 94–96, 117 |

### 7b. Angle to flow — `compute_angle_to_flow.m`

Per-fish (not group) heading relative to a flow direction (`flowAxisDeg`, default 0 = left-to-right).

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `angle_to_flow_deg` | deg | Heading re-centered on the flow axis, wrapped to (−180°, 180°]. | line 70 |
| `mean_angle_to_flow_deg` | deg | **Circular** mean: direction of the average unit heading vector, `atan2d(Σsin, Σcos)` — an arithmetic mean of wrapped angles is garbage when headings straddle ±180°. | lines 75–76 |
| `std_angle_to_flow_deg` | deg | Circular SD: `√(−2·ln R)·180/π`, where R = mean resultant length. | lines 77–81 |
| `range_angle_to_flow_deg` | deg | Wrap-safe angular span measured around the circular mean (179/−179 → ≈2°, not 358°). | lines 83–86 |
| `n_valid_frames` | – | Frames with both landmarks tracked. | line 97 |

### 7c. Distance between individuals — `compute_distance_between_individuals.m`

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `pairwise_dist` | raw × cm_per_unit | Snout-to-snout Euclidean distance matrix per frame (symmetric, NaN where untracked). | lines 71–82 |
| `mean_nn_dist` | raw × cm_per_unit | Per-frame **mean nearest-neighbour** distance: each fish's distance to its single closest neighbour, averaged over present fish. | lines 84–86 |
| `mean_pairwise_dist` | raw × cm_per_unit | Per-frame mean of **all** pairwise distances. | lines 88–90 |
| `overall_mean_nn_dist` / `overall_mean_pairwise_dist` | raw × cm_per_unit | Trial summaries over valid frames. | lines 97–98 |
| `cm_per_unit` | – | Units conversion factor (default 1.0 — distances are in *raw coordinate units* unless you supply your pixel-to-cm calibration). | line 99 |

## 8. Per-cycle frequencies

**File:** `compute_cycle_frequencies.m` — cycle-by-cycle beat frequency instead of one trial-averaged FFT number.

| Metric | Units | Meaning & how it is computed | Code |
|---|---|---|---|
| `freqs_Hz` / `periods_s` | Hz / s | One frequency **per detected cycle**: positive-going zero-crossings of the detrended signal, sub-frame interpolated, `f = 1/period`. | lines 84–111 |
| `cycle_start_frame` / `cycle_end_frame` | frames | Fractional (interpolated) frame index of each cycle boundary. | lines 107–113 |
| `n_cycles` | – | Number of cycles surviving the `min_freq`/`max_freq` gates. | line 114 |
| `mean_freq_Hz` / `std_freq_Hz` | Hz | Mean/SD of the per-cycle frequencies. | lines 117–123 |
| `detrend_window_frames` | frames | Moving-average baseline window actually used (n/5 heuristic) — the first thing to adjust if cycle counts look wrong. | lines 70–82 |

> **Accuracy caveat (from the function's own docstring):** zero-crossing counting is sensitive to slow drift in body-relative Y. On a validated trial it returned 0.21 Hz (1 cycle) where the FFT-based value, a full power-spectrum check, and a visual beat count all agreed on 0.30 Hz. Always sanity-check `n_cycles`/`mean_freq_Hz` against a visual count or `tail_TBF` — or pass `ref_freq_Hz` and the function warns automatically when the two disagree by > 25 %.

---

# Preprocessing

Steps applied between loading and the compute functions; they shape every metric above.

| Step | File / function | What it does | Code |
|---|---|---|---|
| **Midline alignment** | `transform_fish.m` | Per frame: LS line through the **middle points** (all except endpoints); rotate about the y-intercept by `θ = 2π − atan(slope)`; **auto-orient** — whichever *endpoint* has the smaller rotated X is the head (fixes zero-body-length failures when point order is reversed); normalize `X = (x' − head)/bl` (0=head→1=tail), `Y = (y' − a)/bl`, `Z = (z − z_head)/bl`. Exports `transform_params` (a, θ, bl, x1, sign_flip) used downstream. | docstring 26–85; main loop 87 onward (auto-orientation ≈ 182) |
| **Axis remapping** | `remap_axes()` in FishKinematicsApp.m | Permute/negate the 3 dims of `.points` per the GUI's top-X/top-Y/vertical dropdowns (identity when unmapped; skipped for CURVES and 2-D). | lines 2227–2264 |
| **DLC jump filter** | `filter_dlc_jumps.m` | NaN-out tracking jumps: per point, deviation from the local median of neighbours (13-frame window) exceeding `threshold_frac × body_length` (GUI default 0.5) is replaced with NaN — preventing DLC identity swaps from fabricating giant speeds/amplitudes. | whole file |
| **Body-length override** | `transform_fish(…, bl_override)` / GUI body-length field | Known body length replaces the per-frame endpoint distance for BL normalization (needed for Format C, whose endpoints are fin bases). | transform_fish 8–16 |
| **Batch filename parsing** | `parse_filename_tokens()` in FishKinematicsApp.m | Extracts `{speed}`, `{bodylength}`, `{id}` tokens from each batch filename (e.g. `WSBS_*_{speed}BL*_{id}xyzpts`); parsed speed feeds the flow correction (converted to BL/s via `{bodylength}` or the panel body length). | lines 2831–2860 |

---

# CSV export

Runs stage one flat row per metric set (`app.csv_rows`); the export writes them all to one CSV (columns unioned across row types, `row_type` tags each row). The collectors live in FishKinematicsApp.m:

| Row type | Staged by | Notable fields |
|---|---|---|
| `kinematics` | `collect_kine_rows()` (line 2270) | provenance + point identity; `head_TBF_Hz`, `tail_TBF_Hz`, `spline_freq_Hz`, `head/tail_TBF_Z_Hz`; `head/tail_amp_Y_BL`, `head_tail_amp_ratio`, `min/max_amp_Y_BL(+loc)`; Z-amp equivalents; `wavelength_BL`, `wave_speed_BL_s`; `max_curv_XY(+loc)`, `max_curv_3D(+loc)`; extended body: `body_angle_mean/std/range_deg`, `speed_mean/std/peak_BL_s`, `stride_length_BL`, `head_pitch_*`, `roll_*`, `flow_BL_s`, `flow_orientation`, `speed_through_water_*`, `tail_amp_pp_BL`, `strouhal` |
| `fin_kinematics` | `collect_fin_row()` (line 2418) | `fin_mean/std_length`, `fin_total_tip_dist`, `fin_mean/std/peak_speed`, `fin_mean/std/range_{yaw,pitch,roll}_deg`, `fin_mean/std/peak_ang_vel`, `fin_freq_Hz`, `fin_freq_pitch_Hz`, `fin_stride_length_BL`, `fin_stride_duration_s`, `fin_n_valid_frames`, `fin_pct_valid` |
| `girdle_kinematics` | `collect_girdle_row()` (line 2490) | `girdle_point`, `girdle_protraction_range_BL`, `girdle_lateral_range_BL`, `girdle_freq_Hz`, `girdle_n_valid_frames`, `girdle_pct_valid` |
| `stance_swing_estimate` | `collect_stance_row()` (line 2520) | `stance_threshold_used`, `stance_mean_contact_time_s`, `stance_mean_swing_time_s`, `stance_duty_factor`, `stance_n_cycles`, `stance_pct_valid_consecutive`, `stance_is_estimate_not_measured = TRUE` |
| `school_multi_fish` | `collect_polarization_row()` / `collect_angle_to_flow_row()` / `collect_distance_row()` | mean/std polarization; per-fish circular angle-to-flow stats; `mean_nearest_neighbor_dist`, `mean_pairwise_dist`, `cm_per_unit_is_default` |

---

# Repository structure

```
Kinemetrix/
├── FishKinematicsApp.m                 % GUI — all tabs, format detection, CSV collectors
├── load_fish_points.m                  % indexed CSVs (Fish1_P1_x …)
├── load_fish_points_named.m            % named CSVs + dual-camera + compute_pect_phase
├── load_fish_points_dlc_multianimal.m  % native DLC multi-animal exports (struct array)
├── load_fish_curves.m                  % CURVES (Format E) pre-transformed data
├── transform_fish.m                    % midline alignment + BL normalization
├── apply_body_transform.m              % shared per-frame transform (used by girdle)
├── compute_kinematics.m                % FFT analysis: TBF, amplitudes, wavelength, curvature
├── compute_body_extended.m             % speed, stride, Strouhal, body/head/roll angles, flow
├── compute_fin_kinematics.m            % fin Euler angles, fin beat, fin stride
├── compute_girdle_kinematics.m         % girdle protraction/retraction
├── compute_stance_swing.m              % duty-factor estimate from tip speed
├── compute_cycle_frequencies.m         % per-cycle beat frequencies
├── compute_polarization.m              % school alignment metric
├── compute_angle_to_flow.m             % per-fish circular heading vs flow
├── compute_distance_between_individuals.m  % NN + pairwise distances
├── filter_dlc_jumps.m                  % DLC identity-swap jump filter
├── tests/                              % headless ground-truth test suite (12/12 passing)
│   ├── run_all_tests.m                 %   run: matlab -batch "cd('…/tests'); run_all_tests"
│   ├── synth_fish.m                    %   synthetic fish with exact camera geometry
│   └── test_*.m                        %   one file per module
├── Analysis/                           % analysis datasets (e.g. GlassCatfish DLC multi-animal)
└── Demos/                              % demo scripts
```

---

# Test suite

`tests/` holds a synthetic ground-truth suite: `synth_fish.m` generates fish with exactly known camera geometry, kinematics, and body lengths, and each `test_*.m` verifies the corresponding module reproduces them. Run headlessly:

```matlab
matlab -batch "cd('C:/Users/willh/Desktop/fish_analysis_v2/tests'); run_all_tests"
```

Currently **12/12 passing**. The suite has caught real production bugs (head_TBF translation-ramp detrend, dominant-frequency sidelobe leakage, per-cycle detrending) — any change to a metric should keep it green.

---

# Requirements

- MATLAB (recent; uses `readtable` with `VariableNamingRule`, `readcell`, `movmean`, `prctile`). No toolboxes required.

---

# Citation

Midline-alignment method:

> Castro-Santos, T. & Goerig, E. (2017). *Transformer.m* — MATLAB function for aligning fish midlines to the swimming axis.

---

# License

MIT License. See `LICENSE` for details.
