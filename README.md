# Fish Swimming Kinematics Toolkit (Kinemetrix)

A MATLAB toolkit for analysing fish swimming and walking kinematics from 2-D/3-D body-tracking data (e.g. DeepLabCut or DLTdv exports). Given a CSV of digitised landmark coordinates, the pipeline aligns each midline to the swimming axis and extracts tail-beat frequency, lateral amplitude, propulsive wavelength, body curvature, body angle, speed, stride length, head elevation, fin kinematics, girdle protraction/retraction, and an experimental duty-factor estimate.

Developed originally for *Polypterus* locomotion research, extended for frogfish and walking-shark (punting) locomotion, and designed to work with any fish species and any number of tracked body points.

---

## Features

- **Flexible CSV loading** — supports generic indexed columns (`Fish1_P1_x/y/z`), named landmark columns (`snout_X`, `peduncle_Y`, …), dual-camera 2-D columns (`pt1_cam1_X` …), and pre-transformed CURVES/DLTdv exports (`s60mm`, `s70mm`, … station format). Pair-averaging (e.g. left + right pectoral base) is built in.
- **Midline alignment** — rotates each frame so the mean swimming axis is parallel to X. **Head/tail orientation is auto-detected per frame** (via mirroring) rather than assumed from point order, so a dataset whose digitized point order happens to run tail-to-head in the raw coordinate frame is handled automatically — no need to manually reverse column order in your source data.
- **2-D and 3-D support** — all functions auto-detect whether Z data are present and compute dorso-ventral kinematics alongside lateral ones.
- **NaN-safe throughout** — missing/degenerate data (occluded landmarks, untracked frames, zero valid frames) propagate as `NaN` with a clear warning, rather than silently producing a plausible-looking but fabricated number. If a metric is `NaN`, the underlying data genuinely didn't support computing it.
- **FFT interpolation** — upsamples sparse midlines (typically 5–21 points) to 200 evenly-spaced positions along the body for smooth envelope and curvature profiles.
- **Body kinematic outputs** — tail-beat frequency (head and tail separately), lateral and dorso-ventral amplitude envelopes, head/tail amplitude ratio, propulsive wavelength, curvature profile, body angle (heading), forward speed, stride length, head elevation (pitch), and roll (when a paired left/right landmark is available).
- **Fin kinematic outputs** — 3-D fin yaw/pitch/roll angles, fin length, tip speed/distance traveled, angular velocity, fin-beat frequency, and fin stride length/duration.
- **Girdle kinematics** — protraction/retraction range and oscillation frequency of a pectoral or pelvic girdle base point, expressed in the same body-relative reference frame as the midline.
- **Duty factor / contact time (experimental)** — a kinematic *estimate* of stance/swing phase and duty factor from fin-tip velocity thresholding, for walking/punting gaits. This is a proxy, not a substitute for measured (e.g. force-plate) contact-time data — treat it as exploratory.
- **Batch processing (GUI)** — process many CSVs of the same format in one pass via the app's Batch Processing tab: computes and stages results for every file without rendering plots per file, then exports one combined CSV.
- **Ready-made figures** — publication-style plots (midline overlays, amplitude envelope, curvature profile, FFT power spectra, fin trajectories, pectoral phase) generated automatically by the demo scripts and the GUI.
- **Interactive GUI** (`FishKinematicsApp.m`) — load, analyse, visualize, and export kinematics, fin analysis, girdle analysis, and duty-factor estimates without writing any code, including batch mode for large datasets.

---

## Repository structure

```
├── load_fish_points.m           % Load indexed CSV  (Fish1_P1_x …, or pt1_X … numbered format)
├── load_fish_points_named.m     % Load named CSV (snout_X …) or dual-camera (pt1_cam1_X …), with point selection
├── load_fish_curves.m           % Load pre-transformed CURVES/DLTdv station-format CSV or XLS
├── transform_fish.m             % Rotate, translate, and auto-orient midlines to the swimming axis
├── apply_body_transform.m       % Project any OTHER tracked point (fin root, girdle marker) into the
│                                 %   same body-relative frame transform_fish computed for the midline
├── compute_kinematics.m         % FFT analysis → body kinematic struct (TBF, amplitude, wavelength, curvature)
├── compute_body_extended.m      % Body angle, speed, stride length, head elevation, roll
├── compute_fin_kinematics.m     % 3-D fin angle/speed/frequency/stride kinematics from root+tip landmarks
├── compute_girdle_kinematics.m  % Girdle protraction/retraction range and frequency
├── compute_stance_swing.m       % Experimental duty-factor / contact-time estimate from fin-tip velocity
├── FishKinematicsApp.m          % Interactive GUI wrapper, incl. batch processing
├── demo_kinematics.m            % End-to-end demo with 4 figures
├── demo_transform.m             % Before/after midline alignment demo
└── data/
    ├── polypterus_poly1_data.csv
    ├── polypterus_poly2_data.csv
    └── polypterus_poly4_data.csv
```

---

## Quick start

### Option A — named landmark CSV (e.g. the included *Polypterus* data)

```matlab
% Step 1: discover available landmarks
load_fish_points_named('polypterus_poly1_data.csv');

% Step 2: load the points you want, in head-to-tail order
fp = load_fish_points_named('polypterus_poly1_data.csv', ...
       {'snout', {'Rpectbase','Lpectbase'}, 'peduncle', 'caudaltip'}, ...
       [1 2 3 4]);

% Step 3: align midlines to the swimming axis (auto-orients head/tail — see Methods)
fp = transform_fish(fp);

% Step 4: compute body kinematics (100 fps, ignore frequencies below 0.5 Hz)
kine = compute_kinematics(fp, 100, 0.5);

% Step 5: compute extended body kinematics (body angle, speed, stride length, head elevation)
% Pass a {left,right} point-name pair as the 4th argument if you have one, for roll.
ext = compute_body_extended(fp, 100, kine, {});

% Step 6: inspect results
fprintf('Tail-beat frequency: %.2f Hz\n', kine.tail_TBF);
fprintf('Tail amplitude:      %.4f BL\n', kine.tailAmp);
fprintf('Wavelength:          %.4f BL\n', kine.wavelength);
fprintf('Speed:                %.4f BL/s\n', ext.mean_speed_BL_s);
fprintf('Stride length:        %.4f BL/cycle\n', ext.stride_length_BL);
```

### Option B — indexed CSV (`Fish1_P1_x/y/z` or `pt1_X/Y/Z` columns)

```matlab
fp   = load_fish_points('data.csv');
fp   = transform_fish(fp);
kine = compute_kinematics(fp, 100, 0.5);
```

### Option C — pre-transformed CURVES/DLTdv export (`.xls`/.csv station format)

```matlab
% CURVES data is already in body-relative coordinates — skip transform_fish.
fp   = load_fish_curves('trial_CURVES.xls');
kine = compute_kinematics(fp, 30, 0.5);   % use the ACTUAL digitizing fps, not a guess
```

> **Setting `min_freq` correctly matters.** If it's set too high relative to the true beat frequency, the FFT can lock onto a harmonic (e.g. 2x the real frequency) instead of the fundamental. Do a quick visual beat count on a few seconds of video and set `min_freq` comfortably below that, not just at a generic default.

### Option D — fin kinematics from root/tip landmarks

```matlab
fin = compute_fin_kinematics('trial.csv', 'Rpectbase', 'Rpecttip', 100, 0.5, ext.mean_speed_BL_s);
fprintf('Fin-beat frequency: %.2f Hz\n', fin.fin_freq_Hz);
fprintf('Fin stride length:  %.4f BL/cycle\n', fin.stride_length_BL);
```

### Option E — girdle protraction/retraction

```matlab
% Requires fp from transform_fish() (needs its .transform_params field).
girdle = compute_girdle_kinematics('trial.csv', 'Rpectbase', fp, 100, 0.5);
fprintf('Girdle protraction range: %.4f BL\n', girdle.protraction_range_BL);
```

### Option F — duty factor / contact time (experimental — walking/punting gaits only)

```matlab
fin    = compute_fin_kinematics('trial.csv', 'Rpectbase', 'Rpecttip', 100, 0.5);
stance = compute_stance_swing(fin);
fprintf('Duty factor (estimate): %.3f  (coverage: %.1f%%)\n', ...
        stance.duty_factor, stance.pct_valid_consecutive);
% Check stance.pct_valid_consecutive -- below ~70%% means tracking gaps are
% likely fragmenting real stance/swing bouts. Treat low-coverage estimates
% as rough, and prefer measured (e.g. force-plate) contact-time data if
% you have it.
```

### Option G — full demo with figures

```matlab
% Edit CSV_FILE and FPS at the top of the script, then run:
demo_kinematics
```

### Option H — interactive GUI (single file or batch)

```matlab
FishKinematicsApp
```

Load a file, set FPS/min frequency, select points, and hit **Load & Analyse** for body kinematics (with extended metrics shown automatically). Use the **Fin Analysis** tab for fin/girdle/duty-factor analysis on the same file. For many files of the same format, use the **Batch Processing** tab: select all files at once, confirm midline point names, optionally enable fin/girdle/duty-factor analysis (using whichever Root/Tip points are set on the Fin Analysis tab), and run — results stage into the same CSV export as single-file mode, so one export covers the whole batch.

---

## Input CSV formats

### Named landmarks (Format B — recommended for new datasets)

One row per frame, columns named `<landmark>_X`, `<landmark>_Y`, `<landmark>_Z`. Z is optional. Missing values should be `NaN` or a recognized missing-value token (e.g. `NA`) — fully-occluded landmarks (100% missing across the whole file) are detected and handled cleanly rather than crashing.

```
snout_X, snout_Y, snout_Z, peduncle_X, peduncle_Y, peduncle_Z, caudaltip_X, …
-2.09,    4.97,   -7.97,    5.46,       3.74,       -3.71,      7.37, …
```

### Indexed landmarks (Format A/C — DeepLabCut-style multi-animal, or numbered `pt1_X…`)

Columns: `frame`, `FishN_Pk_x`, `FishN_Pk_y` (and optionally `FishN_Pk_z`); or `pt1_X`, `pt1_Y`, `pt1_Z`, `pt2_X`, … with predefined anatomical labels for `pt1`-`pt12`. Multiple fish in one file are loaded as a struct array.

```
frame, Fish1_P1_x, Fish1_P1_y, Fish1_P2_x, Fish1_P2_y, …
0,     312.4,      187.2,      345.1,       190.8, …
```

### Dual-camera 2-D (Format D)

Columns: `pt1_cam1_X`, `pt1_cam1_Y`, `pt1_cam2_X`, `pt1_cam2_Y`, … (`cam1` = lateral view, `cam2` = ventral view). Includes automatic pectoral fin phase (in-phase / antiphase / intermediate) analysis between `pt2` and `pt12`.

### CURVES / DLTdv station export (Format E)

Row 1: body station positions (e.g. every 10mm from snout). Subsequent rows: per-frame `[X, Y]` pairs per station, already in body-relative coordinates. Loaded via `load_fish_curves.m` — **skip `transform_fish`** for this format, it's already aligned.

> **A note on this format's X column:** in some CURVES exports, the X value has been observed to correlate strongly with elapsed time/frame progression rather than being a clean per-frame 0->1 head-to-tail position (i.e. it behaves more like cumulative swim distance than instantaneous body position). This does not affect `head_TBF`/`tail_TBF` (which only use the Y column), but can affect `wavelength`/`maxCurv` if present. If those values look implausible for a CURVES-loaded trial, verify what your specific export's X column actually encodes before trusting spatial metrics.

---

## Kinematic outputs

### `compute_kinematics` (body kinematics)

| Field | Description |
|---|---|
| `head_TBF` / `tail_TBF` | Tail-beat frequency at head / tail (Hz) -- `NaN` if undetectable, never a fabricated value |
| `headAmp` / `tailAmp` | Mean lateral half-amplitude at head / tail (body lengths, BL) |
| `headTailAmpRatio` | Ratio of head to tail amplitude |
| `minAmp` / `minAmpLoc` | Minimum amplitude and its body position (0-1) |
| `maxAmp` / `maxAmpLoc` | Maximum amplitude and its body position (0-1) |
| `wavelength` | Dominant propulsive wavelength (BL) -- resolution is limited to roughly `1/(N_OUT x ds)` per FFT bin; treat as coarse (e.g. "~1 BL" vs "~0.5 BL"), not a fine continuous measurement, unless upsampled |
| `maxCurv` / `maxCurvLoc` | Peak mean curvature and its body position |
| `amp_mean` / `amp_std` | Full amplitude envelope (1 x 200) |
| `curv_mean` / `curv_std` | Full curvature profile (1 x 200) |
| `X_interp` / `Y_interp` | Interpolated midlines (nFrames x 200) |
| `headZ_TBF` / `tailZ_TBF` | Dorso-ventral beat frequency (3-D only) |
| `ampZ_mean` / `ampZ_std` | Dorso-ventral amplitude envelope (3-D only) |
| `n_frames_with_data` | Number of frames with complete X/Y data -- check this if other fields are unexpectedly `NaN` |

### `compute_body_extended` (body angle, speed, stride, head elevation, roll)

| Field | Description |
|---|---|
| `body_angle_deg` / `mean_body_angle_deg` etc. | Per-frame and summary heading angle of the fitted body axis, in the raw/world coordinate frame (turning behavior, not lateral undulation) |
| `angular_velocity_deg_s` | Turning rate (deg/s) |
| `speed_BL_s` / `mean_speed_BL_s` etc. | Forward speed of the body centroid (BL/s) |
| `stride_length_BL` | Mean forward distance traveled per tail-beat cycle (needs `kine.tail_TBF`) |
| `head_pitch_deg` / `mean_head_pitch_deg` etc. | Head elevation angle (3-D only) |
| `head_Z_raw` | Raw vertical head position relative to its own first valid frame (3-D only) |
| `roll_available` / `roll_deg` / `mean_roll_deg` etc. | Roll angle -- **only computed if a `roll_pair` of two left/right landmark names is supplied**; a single midline cannot determine roll on its own. Returns `roll_available = false` and `NaN` otherwise. |

### `compute_fin_kinematics` (3-D fin kinematics)

| Field | Description |
|---|---|
| `yaw` / `pitch` / `roll` | Fin vector Euler angles (deg) over time -- see in-file docstring for exact definitions |
| `fin_length` / `mean_length` / `std_length` | Root-to-tip distance |
| `tip_speed` / `total_dist` / `cum_dist` | Fin tip travel |
| `ang_vel` / `mean_ang_vel` / `peak_ang_vel` | Angular velocity of the fin vector |
| `fin_freq_Hz` | Dominant fin-beat frequency, from the yaw signal |
| `fin_freq_pitch_Hz` | Same, from the pitch signal -- compare the two against a visual beat count if your fin's dominant motion is a flap rather than a sweep |
| `stride_length_BL` | Forward body distance traveled per fin-beat cycle (needs a body speed argument) |
| `stride_duration_s` | Period of one fin-beat cycle |

### `compute_girdle_kinematics` (girdle protraction/retraction)

| Field | Description |
|---|---|
| `X` / `Y` / `Z` | Girdle position in the SAME body-relative frame as the midline |
| `protraction_range_BL` | Fore-aft excursion range of the girdle base itself |
| `lateral_range_BL` | Side-to-side excursion range |
| `girdle_freq_Hz` | Dominant fore-aft oscillation frequency |
| `pct_valid` | % of frames with usable data -- treat low-coverage results as rough estimates |

### `compute_stance_swing` (experimental duty factor)

| Field | Description |
|---|---|
| `is_stance` | Per-frame logical stance/swing classification (velocity-threshold proxy) |
| `mean_contact_time_s` / `mean_swing_time_s` | Mean bout durations |
| `duty_factor` | `contact / (contact + swing)` |
| `n_cycles` | Number of complete stance->swing->stance cycles detected |
| `pct_valid_consecutive` | % of frames with usable frame-to-frame tip speed -- **below ~70% means bouts are likely fragmented by tracking gaps, not real transitions**; treat the estimate as rough and cross-check against any measured contact-time data you have |

---

## Methods

### Midline alignment (`transform_fish`)

For each frame, a least-squares line is fit through the **middle points** (excluding head and tail) in the XY plane. The rotation angle `theta = 2*pi - atan(slope)` is applied about the y-intercept.

**Head/tail orientation is then auto-detected per frame**, rather than assumed from point order: body length is computed as `|x_r(end) - x_r(1)|` (the distance between the two rotated endpoints), which is always non-negative. If point-end happens to rotate to a smaller X than point 1 (which can happen depending on a dataset's raw coordinate handedness -- not necessarily an error in your data), the frame is mirrored about point 1 so that point 1 always maps to `X = 0` and point-end always maps to `X = 1`. A pure mirror preserves all pairwise distances, so amplitude/curvature/wavelength are unaffected -- only which end is labeled "head" changes. This means point order in your source CSV never needs to be manually reversed to match a particular coordinate convention.

Frames where either endpoint is missing, or where a frame's body length collapses to zero/degenerate, are left as `NaN` rather than producing a fabricated position. `transform_fish` reports `pct_frames_valid` and `n_frames_reversed` per animal, and warns if an animal ends up with 0% or under 20% valid frames.

This follows the general approach described in Castro-Santos & Goerig (2017), extended with the auto-orientation step above.

### FFT interpolation (`compute_kinematics`)

Sparse midlines (typically 5-21 landmarks) are upsampled to 200 body positions via zero-padded FFT interpolation. Where the number of input points exceeds 200, cubic spline interpolation is used instead.

### Beat frequency

The dominant tail-beat (or fin-beat, or girdle-oscillation) frequency is extracted from the temporal FFT of the relevant position time series. Only frequencies above `min_freq` are considered, to avoid DC or very-low-frequency drift being selected -- but setting `min_freq` too high relative to the true beat frequency risks the FFT locking onto a harmonic instead of the fundamental. There is no universally "correct" default; set it based on a visual beat count for your specific trial/species.

Degenerate input (all-NaN, or zero-variance after gap-filling) returns `NaN`, not a fabricated frequency -- earlier versions of this toolkit could silently substitute an all-zero signal for missing data and report the lowest allowed FFT bin as if it were real, which looked plausible but wasn't. This has been fixed throughout `compute_kinematics`, `compute_fin_kinematics`, and `compute_girdle_kinematics`.

### Curvature

Three-point geometric (Menger) curvature is computed at each of the 200 interpolated body positions using a configurable lag (default: 5% of body length), then averaged across frames.

### Girdle kinematics (`apply_body_transform` + `compute_girdle_kinematics`)

A fin root/girdle point that wasn't part of the original midline point set can be projected into the exact same body-relative reference frame the midline uses, via `apply_body_transform`, which replays each frame's saved rotation/translation/mirroring parameters (`transform_fish`'s `transform_params` output) against the new point's raw coordinates. This lets girdle protraction/retraction be measured in body-relative units directly comparable to the midline's own kinematics.

### Duty factor estimate (`compute_stance_swing`)

There is no force/contact sensor involved -- stance (fin planted) vs. swing (repositioning) is estimated purely from fin-tip speed, classifying a frame as "stance" when tip speed drops below a threshold fraction (default 25%) of that trial's own 95th-percentile tip speed. This is a standard proxy in locomotion kinematics but is **not the same measurement** as true contact time from a force plate or frame-by-frame visual scoring -- treat measured data as ground truth where you have it, and this as a secondary/exploratory cross-check. The function reports its own reliability (`pct_valid_consecutive`) and warns when tracking gaps make bout detection unreliable.

---

## Known limitations / things to verify on new data

- **`wavelength` has coarse resolution.** With `N_OUT = 200` spanning exactly 1 body length, the FFT's frequency resolution is `~fps_spatial/N` -- practically, this metric distinguishes "about 1 BL" from "about 0.5 BL," not fine differences in between. Spatial upsampling would improve this if needed.
- **`min_freq` is trial-specific, not a safe universal default.** Always sanity-check the returned TBF/fin-freq against a quick visual beat count before trusting it, especially for slow swimmers.
- **Roll requires paired left/right landmarks.** A single midline or a single off-axis marker cannot determine roll on its own; `compute_body_extended`'s `roll_pair` argument needs two named points. A single asymmetric marker (e.g. one dorsal point) can still be analysed via `compute_fin_kinematics`'s `.roll` field, treating a body-centerline reference as the "root" and the marker as the "tip" -- this gives a pseudo-roll proxy, not a fully general roll measurement.
- **Duty factor is a kinematic estimate**, not a measured contact time -- see above.
- **100%-occluded landmarks** (e.g. a fin never visible from the filming angle) are detected and produce clean `NaN` output with a warning, rather than a crash -- but obviously can't be analysed.

---

## Requirements

- MATLAB R2019b or later (uses `readtable`/`readmatrix` with `VariableNamingRule`, `isnan`, `fft`, `polyfit`, `uicheckbox`/`uilistbox` for the GUI)
- No additional toolboxes required

---

## Citation

If you use this toolkit, please cite the midline-alignment method:

> Castro-Santos, T. & Goerig, E. (2017). *Transformer.m* -- MATLAB function for aligning fish midlines to the swimming axis.

---

## License

MIT License. See `LICENSE` for details.
