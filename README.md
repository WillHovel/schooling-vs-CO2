# Kinemetrix: Fish Swimming Kinematics Toolkit

A MATLAB toolkit for analysing fish swimming kinematics from digitised landmark data (DeepLabCut or similar). Given a CSV of tracked body points, the pipeline aligns each midline to the swimming axis, then extracts beat frequency, amplitude envelopes, propulsive wavelength, curvature, speed, stride length, Strouhal number, fin/girdle kinematics, stance–swing duty factor, pectoral phase, and school-level metrics (polarization, angle-to-flow, inter-individual distance).

Built on the midline-alignment method of **Castro-Santos & Goerig (2017)**, with an auto-orientation extension and numerous accuracy fixes (see the `CHANGE NOTE` comments throughout the source). Every metric below is documented in plain language, with its formula and the exact code block that computes it.

> **Metric reference at a glance**: the sections below catalogue every number the program outputs:
> [1. Midline kinematics](#1-midline-kinematics) · [2. Body & gait metrics](#2-body--gait-metrics) · [3. Fin kinematics](#3-fin-kinematics) · [4. Girdle kinematics](#4-girdle-kinematics) · [5. Stance–swing (duty factor)](#5-stanceswing-duty-factor) · [6. Pectoral phase](#6-pectoral-phase) · [7. School metrics](#7-school-metrics) · [8. Per-cycle frequencies](#8-per-cycle-frequencies) · [Preprocessing](#preprocessing) · [CSV export](#csv-export)

---

## Pipeline at a glance

```
CSV file
   │
   ├─ detectFormat()                FishKinematicsApp.m  : classify A/B/C/D/E or DLC multi-animal
   │
   ├─ load_*()                     : points into struct array  fp  (RAW camera frame)
   │     ├─ filter_dlc_jumps()      (optional) NaN-out DLC tracking jumps
   │     ├─ remap_axes()            (optional, 3D)  swap/negate X/Y/Z columns
   │     └─ [school metrics stash]  fp saved RAW before any transform
   │
   ├─ transform_fish()             : per-frame rotate + translate + BL-normalize
   │
   ├─ compute_kinematics()         : FFT analysis  →  kine
   ├─ compute_body_extended()      : speed / stride / Strouhal / angles  →  ext
   │
   ├─ (Fin tab)    compute_fin_kinematics()  → compute_girdle_kinematics()  → compute_stance_swing()
   ├─ (School tab) compute_polarization() / compute_angle_to_flow() / compute_distance_between_individuals()
   │
   └─ collect_*_row()              : stage one flat row per metric set for CSV export
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
| **A: DLC multi-animal** | native DLC export, row 2 = `individuals` (probed with `readcell` **before** `readtable`, which would coerce the text row to NaN) | 4-row header (`scorer`/`individuals`/`bodyparts`/`coords`), N animals × M bodyparts |
| **A: DLC single** | `Fish1_P1_x` style indexed columns | one animal per file, multi-animal columns also loadable |
| **B: Named** | `snout_X`, `peduncle_Y`, … columns | landmark-named CSVs |
| **C: Numbered** | `pt1_X … pt12_X` (12 anatomical points) | numbered-point CSVs |
| **D: Dual-camera** | `pt1_cam1_X …` (cam1 = lateral, cam2 = ventral) | pectoral-phase trials |
| **E: CURVES** | pre-transformed `X/Y` in BL | files already in body frame (`pre_transformed = true` skips `transform_fish`) |

### Units & conventions

- **BL** = body length. By default the per-frame distance between the two endpoint landmarks (in raw units); override with a known body length via `bl_override` / the GUI's body-length field (recommended for Format C, whose endpoints are fin bases, not snout/tail).
- **Speeds** in BL/s. Angles in degrees. Frequencies in Hz.
- **Flow sign convention** (`compute_body_extended.m`): *positive* flow = swimming **against/upstream** (through-water speed = ground speed + flow); *negative* = with the current.
- **Headings** follow `atan2d` (counter-clockwise positive); pectoral "phase shift" is wrapped `mod(…, 360)`.
- **NaN** everywhere means *not computable*: the code deliberately returns NaN rather than a fabricated value for degenerate input (all-NaN, zero variance, too few points). See the `CHANGE NOTE` comments for the specific bugs this guards against.

---

# Metric reference

Every metric the program outputs, grouped by the compute function that produces it. The **Code** column gives the exact file and line range (approximate, as of 2026-08). Throughout, the body is sampled at **200 evenly spaced "stations"** running from the head (station 0) to the tail (station 1), so "at station 0.9" means "90 % of the way from snout to tail".

## 1. Midline kinematics

**File:** `compute_kinematics.m`, **GUI:** Kinematics tab.
Runs after `transform_fish` aligns the fish. For every frame it resamples the tracked midline onto the 200 body stations, cleans the side-to-side motion (removing whole-body sway and the line-fit's own wobble), then extracts beat frequencies, amplitudes, wavelength, and curvature.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `head_TBF` | Hz | **Tail-beat frequency measured at the head**: how many side-to-side beat cycles per second the front of the body completes. This is the same number you'd get by counting beats on video and dividing by trial length, but computed automatically: the program finds the dominant oscillation rate of the head's side-to-side motion, after removing the fish's steady forward drift (which would otherwise masquerade as a slow beat). | `compute_kinematics.m` §4, lines 127–185 (helper `dominant_freq`, lines 425–468) |
| `tail_TBF` | Hz | **Tail-beat frequency at the tail, the classic beat rate** used for stride length and Strouhal number. Computed from the tail's side-to-side motion **relative to the body's midline** (subtracting the whole-body sway that dominates the raw camera signal, while avoiding the line-fit's own oscillation). Measured correctly in validated swimming *and* walking trials; prefer this over `head_TBF`. | `compute_kinematics.m` §4, lines 155–185 |
| `spline_freq_Hz` | Hz | **A second, whole-midline estimate of beat frequency**, the median beat rate across all 200 body stations that show clear oscillation. Use it as a consistency check against `head_TBF`/`tail_TBF`. | `compute_kinematics.m` §4b, lines 198–221 |
| `headZ_TBF` / `tailZ_TBF` | Hz | **Vertical (up-down) beat frequency** at head and tail, the rate of dorso-ventral oscillation in 3-D recordings (e.g. the body's vertical wiggle during swimming). | `compute_kinematics.m` §4, lines 187–196 |
| `amp_mean` / `amp_std` | BL | **The amplitude envelope, how far each body station swings side to side, averaged over the trial** (in body lengths). `amp_mean` is the mean half-amplitude at each of the 200 stations (typically small at the head, growing toward the tail); `amp_std` is how much that amplitude varies between beats. | `compute_kinematics.m` §2, lines 109–116 (helper `amplitude_stats`, lines 355–375) |
| `headAmp` / `tailAmp` | BL | **Average side-to-side excursion of the head region and the tail region** (stations within 5 % of each end). Tail amplitude is the classic "tail-beat amplitude" used in swimming studies. | `compute_kinematics.m` (helper `amplitude_stats`, lines 366–367) |
| `headTailAmpRatio` | – | **Head amplitude ÷ tail amplitude.** Near 0 = most motion in the tail (typical swimming). Higher values mean the head swings too, a stiff-bodied gait or a measurement problem worth checking. | `compute_kinematics.m` line 115 |
| `minAmp`/`maxAmp` + `minAmpLoc`/`maxAmpLoc` | BL / s | **The smallest and largest amplitudes along the body, and where they occur** (location given as body position 0 = head to 1 = tail). The minimum is usually the wave's "node"; the maximum usually sits near the tail. | `compute_kinematics.m` (helper `amplitude_stats`, lines 369–374) |
| `ampZ_mean` / `ampZ_std`, `headAmpZ` … | BL | **The same amplitude set, but for vertical (up-down) motion** instead of side-to-side (3-D recordings only). | `compute_kinematics.m` §3, lines 117–125 |
| `wavelength` | BL | **Propulsive wavelength, the distance between successive bends of the traveling wave along the body**, in body lengths. One full wave of bending per body length ≈ a wavelength of 1. Longer wavelengths mean stiffer, more undulatory swimming; shorter wavelengths mean more flexible, eel-like motion. Computed by fitting a traveling-wave model to the body's motion at the beat frequency (replacing an older method that always returned ≈1 body length regardless of the animal). | `compute_kinematics.m` §5, lines 223–238 (helper `spatial_wavelength`, lines 469–718) |
| `wave_speed_BL_s` | BL/s | **How fast the bend travels down the body**: wavelength × tail-beat frequency. | `compute_kinematics.m` line 238 |
| `curv_mean` / `curv_std` | 1/BL | **The curvature profile, how sharply each body station bends, averaged over the trial** (bend per body length). High values = tightly curved (e.g. near the tail in sharp turns); low values = straight. | `compute_kinematics.m` §6, lines 240–244 (helper `curvature_stats`, lines 378–422) |
| `maxCurv` / `maxCurvLoc` | 1/BL, s | **The sharpest bend anywhere along the body and where it occurs**, useful for identifying where the body bends most during a gait. | `compute_kinematics.m` (helper `curvature_stats`, lines 416–421) |
| `curv3d_mean` … `maxCurv3DLoc` | 1/BL, s | **The same curvature set including vertical (up-down) bending**, true 3-D curvature (3-D recordings only). | `compute_kinematics.m` §7, lines 246–253 |

**Parameters:** `fps` (frame rate of the video), `min_freq` (the lowest beat frequency the program will consider, set it from a visual beat count for your trial; too high a floor makes the frequency finder lock onto a harmonic instead of the true beat).

## 2. Body & gait metrics

**File:** `compute_body_extended.m`, **GUI:** Kinematics tab (runs after `compute_kinematics`; needs its beat frequencies). Also takes an optional **flow speed** (BL/s, e.g. flume current) and, for 3-D data, an optional left/right landmark pair for roll.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `body_angle_deg` + mean/SD/range | deg | **The fish's heading, the direction its body axis points in the video frame.** The trace is "unwrapped" so a fish that slowly turns past ±90° shows a continuous heading instead of a fake 180° jump. Mean/SD/range summarize the average heading and how much it wanders over the trial. | `compute_body_extended.m` §1, lines 107–171 |
| `angular_velocity_deg_s` | deg/s | **Turn rate, how fast the heading changes** from frame to frame. | `compute_body_extended.m` line 162 |
| `speed_BL_s` + mean/SD/peak | BL/s | **Forward swimming speed**, how many body lengths the fish's centre travels per second (mean, variability, and top speed). Computed from the frame-to-frame movement of the body centroid, divided by the measured body length so it is comparable across fish sizes. | `compute_body_extended.m` §2, lines 173–217 |
| `speed_through_water_BL_s` + mean/SD/peak | BL/s | **Speed relative to the water, not the tank**, ground speed corrected for flume current (flow added when swimming upstream, subtracted when swept downstream). This is the physically meaningful speed for flume studies. | `compute_body_extended.m` §2b, lines 219–246 |
| `stride_length_BL` | BL | **How far the fish travels during one complete tail beat**: forward speed ÷ tail-beat frequency. A key measure of gait: each beat should push the fish forward roughly a consistent distance. | `compute_body_extended.m` §3, lines 248–256 |
| `tail_amp_pp_BL` | BL | **Peak-to-peak tail amplitude, total side-to-side excursion of the tail tip across a full beat** (from maximum left to maximum right, ≈ twice the half-amplitude you'd measure by eye). Measured as the tail's distance from a smoothed version of the body axis so the fitted axis's own wobble doesn't inflate it. | `compute_body_extended.m` §3b, lines 266–326 |
| `strouhal` | – | **Strouhal number, a dimensionless measure of swimming efficiency**: (tail-beat frequency × tail amplitude) ÷ swimming speed. Efficient cruising swimmers typically fall in the 0.2–0.4 range; much higher/lower suggests a different gait or an unusual speed. | `compute_body_extended.m` §3b, lines 328–333 |
| `head_pitch_deg` + mean/SD/range | deg | **Head pitch, how much the head tilts up or down relative to horizontal** (+ = head raised). Only for 3-D recordings. | `compute_body_extended.m` §4, lines 338–363 |
| `head_Z_raw` | BL | **Vertical position of the head relative to where it started**, e.g. diving or rising within the tank (3-D only). | `compute_body_extended.m` §4, lines 354–355 |
| `roll_deg` + mean/SD/range | deg | **Body roll, rotation about the long axis** (e.g. banking into turns), measured from a left/right landmark pair (3-D only; `roll_available` flags whether the pair was found). | `compute_body_extended.m` §5, lines 365–395 |

## 3. Fin kinematics

**File:** `compute_fin_kinematics.m`, **GUI:** Fin Analysis (3D) tab.
Loads two landmarks: the fin **root** (base, e.g. `Rpectbase`) and the fin **tip** (e.g. `Rpecttip`), either of which may be a left+right average (`'Rpectbase+Lpectbase'`), and tracks the fin as a vector from root to tip.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `fin_length`, `mean_length`, `std_length` | raw | **Fin size, the root-to-tip distance each frame**, with its mean and variability across the trial. | `compute_fin_kinematics.m` lines 100–101, 169–171 |
| `yaw` / `pitch` / `roll` (+ mean/SD/range) | deg | **Fin orientation relative to the body**: yaw = fore-aft sweep angle, pitch = up-down tilt, roll = rotation about the fin's own long axis. Mean/SD/range describe the fin's resting orientation and how much it moves. | `compute_fin_kinematics.m` lines 103–106, 192–197 |
| `step_dist` / `cum_dist` / `total_dist` | raw | **How far the fin tip travels**, per frame, cumulatively, and in total over the trial. | `compute_fin_kinematics.m` lines 108–119 |
| `tip_speed` + mean/SD/peak | raw/s | **Speed of the fin tip** (per-frame travel × frame rate). | `compute_fin_kinematics.m` lines 119, 150–155 |
| `d_yaw` / `d_pitch` / `d_roll` | deg/s | **How fast the fin is sweeping** (angular velocity about each of the three axes). | `compute_fin_kinematics.m` lines 121–124 |
| `ang_vel` + mean/SD/peak | deg/s | **Total angular speed of the fin**, all three axes combined, the overall "how fast is it flapping". | `compute_fin_kinematics.m` line 125, 202–204 |
| `fin_freq_Hz` | Hz | **Fin-beat frequency, how many full fin beats per second**, from the dominant oscillation of the sweep (yaw) angle. The fin's analogue of tail-beat frequency. | `compute_fin_kinematics.m` lines 127–134 |
| `fin_freq_pitch_Hz` | Hz | **Fin-beat frequency from the flap (pitch) signal instead**, use whichever matches your visual beat count (fins that sweep vs fins that flap differ). | `compute_fin_kinematics.m` lines 132–134 |
| `stride_duration_s` | s | **Length of one fin-beat cycle**: 1 ÷ fin-beat frequency. | `compute_fin_kinematics.m` lines 137–144 |
| `stride_length_BL` | BL | **Body distance travelled per fin beat**: body speed ÷ fin-beat frequency (requires the body speed to be supplied, e.g. `mean_speed_BL_s` from `compute_body_extended.m`). | `compute_fin_kinematics.m` lines 141–143 |
| `n_valid`, `pct_valid` | –, % | **How much of the trial the fin was actually tracked**, valid frames and % of the tracked window (padding rows excluded). Low % = gaps in tracking. | `compute_fin_kinematics.m` lines 208–213 |

## 4. Girdle kinematics

**File:** `compute_girdle_kinematics.m`, **GUI:** Fin Analysis (3D) tab.
Expresses a girdle-base landmark (e.g. `RPectBase`) in the **same body-based coordinate frame as the midline**, so it measures how the girdle itself moves relative to the body, a signal genuinely separate from the fin's own sweep.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `X` / `Y` / `Z` | BL | **Girdle position in body-based coordinates** (X = fore-aft along the body axis, Y = side-to-side). Protraction (girdle swinging forward) moves X toward 0; retraction moves it toward the tail. | `compute_girdle_kinematics.m` lines 114–115 |
| `protraction_range_BL` | BL | **How far the girdle travels fore-aft** (protraction–retraction excursion) in body lengths. | `compute_girdle_kinematics.m` lines 129–130 |
| `lateral_range_BL` | BL | **Side-to-side excursion of the girdle.** | `compute_girdle_kinematics.m` lines 129–130 |
| `girdle_freq_Hz` | Hz | **The girdle's own protraction/retraction beat rate**, the dominant oscillation frequency of its fore-aft motion. | `compute_girdle_kinematics.m` lines 132–136 |
| `n_valid` / `pct_valid` | –, % | **Tracking coverage** for the girdle point (% of the tracked window). | `compute_girdle_kinematics.m` lines 117–121 |

## 5. Stance–swing (duty factor)

**File:** `compute_stance_swing.m`, **GUI:** Fin Analysis (3D) tab.
> **Important:** "contact" here is a **kinematic proxy**, not a force measurement: a frame is called STANCE (fin planted) when the fin tip's speed drops below 25 % (default) of the trial's own top-end (95th-percentile) tip speed, and SWING (repositioning) otherwise. If you have real contact-time data (e.g. from a pressure mat), treat that as ground truth and this as a secondary check.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `threshold_used` | raw/s | **The actual speed cutoff used** to decide "planted" vs "moving" (25 % of the trial's 95th-percentile tip speed), reported so you know exactly where the line was drawn. | `compute_stance_swing.m` lines 92–96 |
| `is_stance` | – | **Per-frame flag: is the fin planted (slow) or swinging (moving)?** Tracking gaps count as *unknown*, never as planted. | `compute_stance_swing.m` lines 95–96 |
| `mean_contact_time_s` / `mean_swing_time_s` | s | **Average time the fin stays planted per step, and average time it spends moving** to the next position (after bridging short tracking gaps and discarding very short blips). | `compute_stance_swing.m` lines 98–134 |
| `duty_factor` | – | **The fraction of each step cycle the fin spends planted**: 0.5 = equal plant and move time, higher = more time planted (a "heavier" gait). The classic duty factor. | `compute_stance_swing.m` line 133 |
| `n_cycles` | – | **Number of complete stance→swing→stance cycles detected** in the trial. | `compute_stance_swing.m` line 136 |
| `pct_valid_consecutive` | % | **Coverage check, % of frames with usable frame-to-frame speed.** Below ≈70 % this estimate is likely unreliable (tracking gaps fragment the real bouts); check this before trusting `duty_factor`. | `compute_stance_swing.m` lines 76–84 |

## 6. Pectoral phase

**File:** `compute_pect_phase` (a local function inside `load_fish_points_named.m`, lines 227–320), **GUI:** Pectoral Phase (2D) tab (results shown by `displayPectPhase`, FishKinematicsApp.m line 2139).
For Format D dual-camera files: picks the best-tracked camera for each fin, then compares the beat timing of the right pectoral (pt2) vs left pectoral (pt12) tip signals.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `phase_shift_deg` | deg | **The timing difference between the right and left pectoral beats**, in degrees of the beat cycle: 0° = both fins move together (in-phase), 180° = exactly opposite (antiphase). Computed from the phase of each fin's signal at its dominant beat frequency. | `load_fish_points_named.m` lines 279–292 |
| `phase_right_deg` / `phase_left_deg` | deg | **Each fin's individual phase** at the beat frequency, the raw numbers behind the shift. | `load_fish_points_named.m` lines 289–290 |
| `classification` | – | **Plain-language verdict**: `In-phase` (shift ≤ 45° or ≥ 315°), `Antiphase` (135–225°), or `Intermediate` (anything else). | `load_fish_points_named.m` lines 294–300 |
| `peak_lag_frames` | frames | **How many frames one fin trails the other**, from the peak of the cross-correlation between the two signals (positive = left lags right). | `load_fish_points_named.m` lines 302–304 |
| `dom_freq_norm` | cycles/frame | **The pectoral beat rate** (the dominant frequency, in cycles per frame). | `load_fish_points_named.m` line 310 |
| `n_valid` | – | **Frames where both fins were tracked** (needs ≥ 20, otherwise reported as `N/A`). | `load_fish_points_named.m` lines 266–274 |

## 7. School metrics

**School Metrics tab** (GUI), computed from the **RAW, untransformed** tracking data (stashed before `transform_fish`), because each fish's individual body-normalized frame would destroy the shared spatial relationships these metrics need. Loads native DLC multi-animal exports directly (Format A), or any multi-fish file via the tab's own Browse.

### 7a. Polarization: `compute_polarization.m`

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `polarization` | 0–1 | **How aligned the fish are with each other**, the core schooling metric: 1 = every fish facing exactly the same way, 0 = directions random or opposed. Computed per frame as the length of the average unit heading vector (each fish's heading = its peduncle→snout direction), a standard circular-statistics measure that ignores *which* way they face and only measures *how agreed* they are. NaN when fewer than 2 fish are tracked that frame. | `compute_polarization.m` lines 99–107 |
| `mean_polarization` / `std_polarization` | 0–1 | **Average alignment across the trial**, and how much it fluctuates. | `compute_polarization.m` lines 109–115 |
| `heading_deg` / `heading_to_flow_deg` | deg | **Each fish's heading**, raw and re-centered so 0 = aligned with the flow direction you entered. | `compute_polarization.m` lines 94–96, 117 |

### 7b. Angle to flow: `compute_angle_to_flow.m`

Per-fish (not group) heading relative to a flow direction (`flowAxisDeg`, default 0 = left-to-right).

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `angle_to_flow_deg` | deg | **Each fish's heading relative to the flow**, 0 = aligned with the flow axis, ±180 = pointing the opposite way. | `compute_angle_to_flow.m` line 70 |
| `mean_angle_to_flow_deg` | deg | **Average heading relative to flow, computed "circularly"**, the direction of the average unit heading vector. This matters: a plain arithmetic average of angles straddling ±180° (e.g. fish at +179° and −179°, both heading into the flow) would wrongly average to ~0° ("pointing every which way"). | `compute_angle_to_flow.m` lines 75–76 |
| `std_angle_to_flow_deg` | deg | **Circular standard deviation**, how much the heading varies around its average (also wrap-safe). | `compute_angle_to_flow.m` lines 77–81 |
| `range_angle_to_flow_deg` | deg | **The total angular spread of the headings, measured around the circular mean**, e.g. +179/−179 gives ≈2°, not 358°. | `compute_angle_to_flow.m` lines 83–86 |
| `n_valid_frames` | – | **Frames where both landmarks were tracked** for this fish. | `compute_angle_to_flow.m` line 97 |

### 7c. Distance between individuals: `compute_distance_between_individuals.m`

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `pairwise_dist` | raw × cm_per_unit | **Snout-to-snout distance between every pair of fish, per frame** (symmetric; NaN where a fish's snout isn't tracked). | `compute_distance_between_individuals.m` lines 71–82 |
| `mean_nn_dist` | raw × cm_per_unit | **Mean nearest-neighbour distance, how close each fish is to its single closest neighbour**, averaged over the fish present that frame. The classic "how tightly are they schooling" number. | `compute_distance_between_individuals.m` lines 84–86 |
| `mean_pairwise_dist` | raw × cm_per_unit | **Mean of ALL pairwise distances** that frame (not just nearest neighbours), a broader measure of group spread. | `compute_distance_between_individuals.m` lines 88–90 |
| `overall_mean_nn_dist` / `overall_mean_pairwise_dist` | raw × cm_per_unit | **Trial-wide averages** of the two measures above. | `compute_distance_between_individuals.m` lines 97–98 |
| `cm_per_unit` | – | **Unit conversion factor** (default 1.0, distances are in your raw coordinate units unless you supply a pixel-to-cm calibration). | `compute_distance_between_individuals.m` line 99 |

## 8. Per-cycle frequencies

**File:** `compute_cycle_frequencies.m`, cycle-by-cycle beat frequency instead of one trial-averaged number.

| Metric | Units | What it means & how it is computed | Code |
|---|---|---|---|
| `freqs_Hz` / `periods_s` | Hz / s | **One frequency per individual beat cycle** (rather than one average for the whole trial), lets you see whether the beat rate is steady or varies cycle to cycle. Each cycle is detected as a rising crossing of the signal's zero level, and its period = time between crossings. | `compute_cycle_frequencies.m` lines 84–111 |
| `cycle_start_frame` / `cycle_end_frame` | frames | **The frame at which each cycle begins and ends** (interpolated to sub-frame precision). | `compute_cycle_frequencies.m` lines 107–113 |
| `n_cycles` | – | **Number of complete cycles detected** (after applying any min/max frequency gates you set). | `compute_cycle_frequencies.m` line 114 |
| `mean_freq_Hz` / `std_freq_Hz` | Hz | **Average cycle frequency and its variability.** | `compute_cycle_frequencies.m` lines 117–123 |
| `detrend_window_frames` | frames | **An internal smoothing window** used to remove slow drift before counting cycles, reported because it's the first thing to adjust if your cycle counts look wrong. | `compute_cycle_frequencies.m` lines 70–82 |

> **Accuracy caveat (from the function's own documentation):** counting cycles by zero-crossings is sensitive to slow drift in the body-relative signal. On one validated trial it returned 0.21 Hz (1 cycle) where the FFT-based value, a full power-spectrum check, and a visual beat count all agreed on 0.30 Hz. Always sanity-check `n_cycles`/`mean_freq_Hz` against a visual beat count or `tail_TBF`, or pass `ref_freq_Hz` and the function warns automatically when the two disagree by more than 25 %.

---

# Preprocessing

Steps applied between loading and the compute functions; they shape every metric above.

| Step | File / function | What it does | Code |
|---|---|---|---|
| **Midline alignment** | `transform_fish.m` | Straightens each frame so the body axis lies along X with the head at 0, in body-length units, the coordinate frame every midline metric lives in. Per frame: fits a line through the middle points, rotates the fish to it, auto-detects which end is the head, and normalizes by body length. | `transform_fish.m` docstring lines 26–85; main loop line 87 onward (auto-orientation ≈ line 182) |
| **Axis remapping** | `remap_axes()` in FishKinematicsApp.m | Swaps/negates the X/Y/Z columns per the GUI's axis dropdowns (e.g. if the camera's "up" is the CSV's Z, tell the app and it fixes the frame before analysis). | `FishKinematicsApp.m` lines 2227–2264 |
| **DLC jump filter** | `filter_dlc_jumps.m` | Removes tracking jumps: if a point jumps more than a threshold × body length away from its neighbours' local median (13-frame window), that frame is marked untracked, preventing DLC identity swaps from producing giant fake speeds/amplitudes. | `filter_dlc_jumps.m` (whole file) |
| **Body-length override** | `transform_fish(…, bl_override)` / GUI body-length field | Uses a known body length instead of the measured endpoint distance for all BL-normalized metrics, recommended for Format C, whose endpoints are fin bases, not the true snout/tail. | `transform_fish.m` lines 8–16 |
| **Batch filename parsing** | `parse_filename_tokens()` in FishKinematicsApp.m | Reads `{speed}`, `{bodylength}`, `{id}` values out of each batch filename (e.g. `WSBS_*_{speed}BL*_{id}xyzpts`) and uses the parsed speed as that file's flume flow for the through-water speed correction. | `FishKinematicsApp.m` lines 2831–2860 |

---

# CSV export

Each analysis run stages one flat row per metric set (`app.csv_rows`); the export writes them all to one CSV (columns unioned across row types, with `row_type` tagging each row). The collectors all live in FishKinematicsApp.m:

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
├── FishKinematicsApp.m                 % GUI, all tabs, format detection, CSV collectors
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

Currently **12/12 passing**. The suite has caught real production bugs (head_TBF translation-ramp detrend, dominant-frequency sidelobe leakage, per-cycle detrending); any change to a metric should keep it green.

---

# Requirements

- MATLAB (recent; uses `readtable` with `VariableNamingRule`, `readcell`, `movmean`, `prctile`). No toolboxes required.

---

# Citation

Midline-alignment method:

> Castro-Santos, T. & Goerig, E. (2017). *Transformer.m*, MATLAB function for aligning fish midlines to the swimming axis.

---

# License

MIT License. See `LICENSE` for details.
