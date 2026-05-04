# Fish Swimming Kinematics Toolkit

A MATLAB toolkit for analysing fish swimming kinematics from 3-D body-tracking data. Given a CSV of digitised landmark coordinates (e.g. from DeepLabCut), the pipeline aligns each midline to the swimming axis and extracts tail-beat frequency, lateral amplitude, propulsive wavelength, and body curvature via FFT interpolation.

Developed for *Polypterus* locomotion research, but designed to work with any fish species and any number of tracked body points.

---

## Features

- **Flexible CSV loading** — supports both generic indexed columns (`Fish1_P1_x/y/z`) and named landmark columns (`snout_X`, `peduncle_Y`, …). Pair-averaging (e.g. left + right pectoral base) is built in.
- **Midline alignment** — rotates each frame so the mean swimming axis is parallel to X, with the snout at the origin (follows Castro-Santos & Goerig 2017).
- **2-D and 3-D support** — all functions auto-detect whether Z data are present and compute dorso-ventral kinematics alongside lateral ones.
- **FFT interpolation** — upsamples sparse midlines (typically 5–10 points) to 200 evenly-spaced positions along the body for smooth envelope and curvature profiles.
- **Kinematic outputs** — tail-beat frequency (head and tail separately), lateral and dorso-ventral amplitude envelopes, head/tail amplitude ratio, propulsive wavelength, and curvature profile.
- **Ready-made figures** — four publication-style plots (midline overlays, amplitude envelope, curvature profile, FFT power spectra) generated automatically.

---

## Repository structure

```
├── load_fish_points.m        % Load indexed CSV  (Fish1_P1_x …)
├── load_fish_points_named.m  % Load named CSV    (snout_X …) with point selection
├── transform_fish.m          % Rotate & translate midlines
├── compute_kinematics.m      % FFT analysis → kinematic struct
├── FishKinematicsApp.m       % Interactive GUI wrapper
├── demo_kinematics.m         % End-to-end demo with 4 figures
├── demo_transform.m          % Before/after midline alignment demo
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

% Step 3: align midlines to the swimming axis
fp = transform_fish(fp);

% Step 4: compute kinematics (100 fps, ignore frequencies below 0.5 Hz)
kine = compute_kinematics(fp, 100, 0.5);

% Step 5: inspect results
fprintf('Tail-beat frequency: %.2f Hz\n', kine.tail_TBF);
fprintf('Tail amplitude:      %.4f BL\n', kine.tailAmp);
fprintf('Wavelength:          %.4f BL\n', kine.wavelength);
```

### Option B — indexed CSV (`Fish1_P1_x/y/z` columns)

```matlab
fp   = load_fish_points('data.csv');
fp   = transform_fish(fp);
kine = compute_kinematics(fp, 100, 0.5);
```

### Option C — full demo with figures

```matlab
% Edit CSV_FILE and FPS at the top of the script, then run:
demo_kinematics
```

---

## Input CSV format

### Named landmarks (recommended for new datasets)

One row per frame, columns named `<landmark>_X`, `<landmark>_Y`, `<landmark>_Z`. Z is optional. Missing values should be `NaN`.

```
snout_X, snout_Y, snout_Z, peduncle_X, peduncle_Y, peduncle_Z, caudaltip_X, …
-2.09,    4.97,   -7.97,    5.46,       3.74,       -3.71,      7.37, …
```

### Indexed landmarks (DeepLabCut-style multi-animal)

Columns: `frame`, `FishN_Pk_x`, `FishN_Pk_y` (and optionally `FishN_Pk_z`). Multiple fish in one file are loaded as a struct array.

```
frame, Fish1_P1_x, Fish1_P1_y, Fish1_P2_x, Fish1_P2_y, …
0,     312.4,      187.2,      345.1,       190.8, …
```

---

## Kinematic outputs

`compute_kinematics` returns a struct (one element per animal) with the following fields.

| Field | Description |
|---|---|
| `head_TBF` / `tail_TBF` | Tail-beat frequency at head / tail (Hz) |
| `headAmp` / `tailAmp` | Mean lateral half-amplitude at head / tail (body lengths, BL) |
| `headTailAmpRatio` | Ratio of head to tail amplitude |
| `minAmp` / `minAmpLoc` | Minimum amplitude and its body position (0–1) |
| `maxAmp` / `maxAmpLoc` | Maximum amplitude and its body position (0–1) |
| `wavelength` | Dominant propulsive wavelength (BL) |
| `maxCurv` / `maxCurvLoc` | Peak mean curvature and its body position |
| `amp_mean` / `amp_std` | Full amplitude envelope (1 × 200) |
| `curv_mean` / `curv_std` | Full curvature profile (1 × 200) |
| `X_interp` / `Y_interp` | Interpolated midlines (nFrames × 200) |
| `headZ_TBF` / `tailZ_TBF` | Dorso-ventral beat frequency (3-D only) |
| `ampZ_mean` / `ampZ_std` | Dorso-ventral amplitude envelope (3-D only) |

---

## Methods

### Midline alignment (`transform_fish`)

For each frame, a least-squares line is fit through the **middle points** (excluding head and tail) in the XY plane. The rotation angle `θ = 2π − atan(slope)` is applied about the y-intercept, then the midline is translated so the head sits at X = 0. This follows the approach described in Castro-Santos & Goerig (2017).

### FFT interpolation (`compute_kinematics`)

Sparse midlines (typically 5–10 landmarks) are upsampled to 200 body positions via zero-padded FFT interpolation. Where the number of input points exceeds 200, cubic spline interpolation is used instead.

### Beat frequency

The dominant tail-beat frequency is extracted from the temporal FFT of the head (or tail) Y-position time series. Only frequencies above `min_freq` are considered, preventing DC or very-low-frequency noise from being selected.

### Curvature

Three-point geometric curvature is computed at each of the 200 interpolated body positions using a configurable lag (default: 5% of body length), then averaged across frames.

---

## Requirements

- MATLAB R2019b or later (uses `readtable` with `VariableNamingRule`, `isnan`, `fft`, `polyfit`)
- No additional toolboxes required

---

## Citation

If you use this toolkit, please cite the midline-alignment method:

> Castro-Santos, T. & Goerig, E. (2017). *Transformer.m* — MATLAB function for aligning fish midlines to the swimming axis.

---

## License

MIT License. See `LICENSE` for details.
