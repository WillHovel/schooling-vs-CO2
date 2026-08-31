#!/usr/bin/env python3
"""
Generate realistic, fabricated fish tracking CSVs for the Kinemetrix toolkit.

One CSV (or set of CSVs) per GUI tab, all built from the same travelling-wave
fish model used by the synthetic test suite (tests/synth_fish.m), but dressed
up with realistic landmark anatomy, fin beating, schooling, and dual-camera
views so they read naturally in a demonstration / tutorial video.

Outputs (written into this script's directory):

  demo_kinematics.csv                 Kinematics tab     (Format A: DLC single, 3D)
  demo_fin.csv                        Fin Analysis tab   (Format B: named, 3D)
  demo_pectoral_inphase.csv           Pectoral Phase tab (Format D: dual-camera)
  demo_pectoral_antiphase.csv         Pectoral Phase tab (Format D: dual-camera)
  demo_school_3fish.csv               School Metrics tab (Format A: DLC multi-animal)
  batch_trials/WSBS_Swim_*_Fish*.csv  Batch Processing tab (Format B: named, 3D)

Run:  python generate_demo_data.py
"""

import csv
import math
import os
import random

OUT = os.path.dirname(os.path.abspath(__file__))
BATCH = os.path.join(OUT, "batch_trials")
os.makedirs(BATCH, exist_ok=True)

random.seed(42)


# ---------------------------------------------------------------------------
# Coordinate helpers (mirror the transform in tests/synth_fish.m)
# ---------------------------------------------------------------------------
def rot_xy(xb, yb, heading):
    x = math.cos(heading) * xb - math.sin(heading) * yb
    y = math.sin(heading) * xb + math.cos(heading) * yb
    return x, y


def swim_center(t, U, heading, c0):
    cx = c0[0] + U * math.cos(heading) * t
    cy = c0[1] + U * math.sin(heading) * t
    return cx, cy


def fmt(v, nd=4):
    return f"{v:.{nd}f}"


def jitter():
    return random.gauss(0.0, 0.15)


def write_csv(path, header, rows):
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(header)
        w.writerows(rows)
    print(f"  wrote {path}  ({len(rows)} rows x {len(header)} cols)")


# ---------------------------------------------------------------------------
# Shared fish parameters
# ---------------------------------------------------------------------------
FPS = 100
L = 120.0          # body length in camera units (px)
TAIL_FREQ = 5.0    # tail-beat frequency, Hz
LAMBDA = 1.0       # propulsive wavelength, BL
HEADING = math.radians(-12.0)   # swim up-and-to-the-right
C0 = (220.0, 480.0)


def amp_fn(s):
    # amplitude envelope in camera units: grows head -> tail (0.015 -> 0.045 BL)
    # small enough that Strouhal ~ 0.3 (the efficient-cruising sweet spot)
    return L * (0.015 + 0.03 * s)


def body_xy(t, s):
    xb = -s * L
    yb = amp_fn(s) * math.sin(2 * math.pi * (TAIL_FREQ * t - s / LAMBDA))
    return xb, yb


def body_z(t, s):
    # gentle dorso-ventral wave so the data is genuinely 3D but well behaved
    return L * 0.03 * math.sin(2 * math.pi * (TAIL_FREQ * t - s / LAMBDA) + 0.7)


def add_jitter(x, y, z=None):
    x += jitter()
    y += jitter()
    if z is not None:
        z += jitter()
        return x, y, z
    return x, y


# ---------------------------------------------------------------------------
# 1. Kinematics tab  --  Format A (DLC single), 3D, Fish1_P1_x/y/z ...
# ---------------------------------------------------------------------------
def gen_kinematics():
    print("Kinematics tab (Format A DLC single, 3D):")
    n_frames = 600
    n_pts = 13
    s = [i / (n_pts - 1) for i in range(n_pts)]
    U = 1.5 * L  # 1.5 BL/s ground speed

    cols = ["frame"]
    for p in range(1, n_pts + 1):
        cols += [f"Fish1_P{p}_x", f"Fish1_P{p}_y", f"Fish1_P{p}_z"]

    rows = []
    for f in range(n_frames):
        t = f / FPS
        cx, cy = swim_center(t, U, HEADING, C0)
        row = [f]
        for si in s:
            xb, yb = body_xy(t, si)
            zb = body_z(t, si)
            x, y = rot_xy(xb, yb, HEADING)
            x, y, z = add_jitter(cx + x, cy + y, zb)
            row += [fmt(x), fmt(y), fmt(z)]
        rows.append(row)
    write_csv(os.path.join(OUT, "demo_kinematics.csv"), cols, rows)


# ---------------------------------------------------------------------------
# 2. Fin Analysis tab  --  Format B (named), 3D
#    Columns ordered so the fin root/tip dropdowns default to Rpectbase /
#    Rpecttip. Midline landmarks (snout..caudaltip) selected in the main tab.
# ---------------------------------------------------------------------------
def gen_fin():
    print("Fin Analysis tab (Format B named, 3D):")
    n_frames = 600
    U = 1.5 * L
    FF = 3.0  # pectoral fin-beat frequency, Hz

    # station positions along the body for each landmark
    s_root = 0.18
    s_mid = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]  # snout..caudaltip

    names = [
        "Rpectbase", "Rpecttip", "Lpectbase", "Lpecttip",
        "Rpelvicbase", "Rpelvictip", "Lpelvicbase", "Lpelvictip",
        "snout", "mid1", "mid2", "mid3", "peduncle", "caudaltip",
    ]
    cols = []
    for n in names:
        cols += [f"{n}_X", f"{n}_Y", f"{n}_Z"]

    def fin_root_tip(sr, side, phase):
        """side = +1 (right) or -1 (left); phase shifts the beat."""
        rows = []
        for f in range(n_frames):
            t = f / FPS
            cx, cy = swim_center(t, U, HEADING, C0)
            # root: attached to body at station sr, offset laterally
            xb_r, yb_r = body_xy(t, sr)
            xb_r += 0.0
            yb_r += side * 0.06 * L
            zb_r = body_z(t, sr) + 0.02 * L
            # tip: sweeps fore-aft (x) and flaps up-down (z) around the root.
            # A tanh-sharpened sine dwells near the extremes (slow -> stance)
            # and transitions quickly (fast -> swing), so the stance/swing
            # duty factor is meaningful rather than NaN.
            sharp = math.tanh(2.5 * math.sin(2 * math.pi * FF * t + phase))
            xb_t = xb_r - 0.08 * L * sharp
            yb_t = yb_r + side * (0.20 * L + 0.02 * L * sharp)
            zb_t = zb_r - 0.12 * L * sharp
            xr, yr = rot_xy(xb_r, yb_r, HEADING)
            xt, yt = rot_xy(xb_t, yb_t, HEADING)
            xr, yr, zr = add_jitter(cx + xr, cy + yr, zb_r)
            xt, yt, zt = add_jitter(cx + xt, cy + yt, zb_t)
            rows.append((xr, yr, zr, xt, yt, zt))
        return rows

    R_root = fin_root_tip(s_root, +1, 0.0)
    L_root = fin_root_tip(s_root, -1, math.pi)

    data = {n: [] for n in names}
    for f in range(n_frames):
        t = f / FPS
        cx, cy = swim_center(t, U, HEADING, C0)
        # pectoral fins
        r = R_root[f]
        l = L_root[f]
        data["Rpectbase"].append([fmt(r[0]), fmt(r[1]), fmt(r[2])])
        data["Rpecttip"].append([fmt(r[3]), fmt(r[4]), fmt(r[5])])
        data["Lpectbase"].append([fmt(l[0]), fmt(l[1]), fmt(l[2])])
        data["Lpecttip"].append([fmt(l[3]), fmt(l[4]), fmt(l[5])])
        # pelvic fins: small, slow, near station 0.55
        for side, base, tip in ((+1, "Rpelvicbase", "Rpelvictip"), (-1, "Lpelvicbase", "Lpelvictip")):
            xb_b, yb_b = body_xy(t, 0.55)
            yb_b += side * 0.05 * L
            zb_b = body_z(t, 0.55)
            xb_t = xb_b + side * 0.06 * L * math.sin(2 * math.pi * 2.0 * t)
            yb_t = yb_b + side * 0.10 * L
            zb_t = zb_b - 0.06 * L * math.sin(2 * math.pi * 2.0 * t + 0.5)
            xb, yb = rot_xy(xb_b, yb_b, HEADING)
            xt_, yt_ = rot_xy(xb_t, yb_t, HEADING)
            xb, yb, zb = add_jitter(cx + xb, cy + yb, zb_b)
            xt_, yt_, zt_ = add_jitter(cx + xt_, cy + yt_, zb_t)
            data[base].append([fmt(xb), fmt(yb), fmt(zb)])
            data[tip].append([fmt(xt_), fmt(yt_), fmt(zt_)])
        # midline
        for nm, si in zip(["snout", "mid1", "mid2", "mid3", "peduncle", "caudaltip"], s_mid):
            xb, yb = body_xy(t, si)
            zb = body_z(t, si)
            xb, yb = rot_xy(xb, yb, HEADING)
            xb, yb, zb = add_jitter(cx + xb, cy + yb, zb)
            data[nm].append([fmt(xb), fmt(yb), fmt(zb)])

    rows = []
    for f in range(n_frames):
        row = []
        for n in names:
            row += data[n][f]
        rows.append(row)
    write_csv(os.path.join(OUT, "demo_fin.csv"), cols, rows)


# ---------------------------------------------------------------------------
# 3. Pectoral Phase tab  --  Format D (dual-camera), 2D
#    pt2 (right pectoral tip) vs pt12 (left pectoral) Y signals beat at 3 Hz.
#    phase=0 -> In-phase, phase=pi -> Antiphase.
# ---------------------------------------------------------------------------
def gen_pectoral(phase_rad, fname):
    print(f"Pectoral Phase tab (Format D dual-camera, 2D): {os.path.basename(fname)}")
    n_frames = 400
    FP = 3.0
    A = 20.0

    pts = list(range(1, 13))
    cams = [1, 2]
    cols = []
    for p in pts:
        for c in cams:
            cols += [f"pt{p}_cam{c}_X", f"pt{p}_cam{c}_Y"]

    # nominal resting (x, y) for each point, body centred ~ (400, 400)
    resting = {
        1: (350, 430), 2: (360, 455), 3: (180, 400), 4: (80, 400),
        5: (395, 390), 6: (250, 415), 7: (245, 440), 8: (210, 408),
        9: (205, 435), 10: (230, 360), 11: (395, 410), 12: (350, 370),
    }

    rows = []
    for f in range(n_frames):
        t = f / FPS
        row = []
        for p in pts:
            for c in cams:
                rx, ry = resting[p]
                if p == 2:   # right pectoral tip
                    if c == 2:  # ventral view: lateral beat is the Y axis
                        x = rx + 8 * math.sin(2 * math.pi * FP * t)
                        y = ry + A * math.sin(2 * math.pi * FP * t)
                    else:       # lateral view: smaller vertical beat
                        x = rx + 8 * math.sin(2 * math.pi * FP * t)
                        y = ry + 6 * math.sin(2 * math.pi * FP * t)
                elif p == 12:  # left pectoral
                    if c == 2:
                        x = rx + 8 * math.sin(2 * math.pi * FP * t + phase_rad)
                        y = ry + A * math.sin(2 * math.pi * FP * t + phase_rad)
                    else:
                        x = rx + 8 * math.sin(2 * math.pi * FP * t + phase_rad)
                        y = ry + 6 * math.sin(2 * math.pi * FP * t + phase_rad)
                else:           # body points: gentle sway
                    x = rx + 1.5 * math.sin(2 * math.pi * FP * t)
                    y = ry + 1.0 * math.cos(2 * math.pi * FP * t)
                x += jitter()
                y += jitter()
                row += [fmt(x), fmt(y)]
        rows.append(row)
    write_csv(fname, cols, rows)


# ---------------------------------------------------------------------------
# 4. School Metrics tab  --  Format A (DLC multi-animal), 2D
# ---------------------------------------------------------------------------
def gen_school():
    print("School Metrics tab (Format A DLC multi-animal):")
    n_frames = 400
    U = 1.0 * L
    bodyparts = ["snout", "mid1", "mid2", "mid3", "peduncle"]
    n_indiv = 3
    # per-fish (x0, y0, heading-jitter phase) -> all swim roughly +x
    fish_cfg = [
        (100.0, 200.0, 0.0),
        (60.0, 260.0, 1.3),
        (40.0, 150.0, 2.1),
    ]
    # body-frame offsets along heading axis (positive = ahead)
    offsets = {"snout": 0.05 * L, "mid1": -0.15 * L, "mid2": -0.35 * L,
               "mid3": -0.55 * L, "peduncle": -0.85 * L}
    scorer = "DLC_resnet50_demoSchoolJul1shuffle1_100000"

    # build 4-row header
    header = []
    header.append(["scorer"] + [scorer] * (len(bodyparts) * 3 * n_indiv))
    header.append(["individuals"] + [f"individual{i}" for i in range(1, n_indiv + 1)
                                     for _ in range(len(bodyparts) * 3)])
    header.append(["bodyparts"] + [bp for _ in range(n_indiv) for bp in bodyparts
                                   for _ in range(3)])
    header.append(["coords"] + ["x", "y", "likelihood"] * (len(bodyparts) * n_indiv))

    rows = []
    for f in range(n_frames):
        t = f / FPS
        row = [f]
        for (x0, y0, ph) in fish_cfg:
            h = math.radians(3.0 * math.sin(2 * math.pi * 0.15 * t + ph))  # +/-3 deg wander
            cx = x0 + U * t
            cy = y0 + 1.5 * math.sin(2 * math.pi * 0.1 * t + ph)
            for bp in bodyparts:
                xb = offsets[bp]
                yb = amp_fn(abs(offsets[bp]) / L) * 0.3 * math.sin(2 * math.pi * (TAIL_FREQ * t - abs(offsets[bp]) / L))
                x, y = rot_xy(xb, yb, h)
                x, y = add_jitter(cx + x, cy + y)
                lik = fmt(random.uniform(0.95, 1.0))
                row += [fmt(x), fmt(y), lik]
        rows.append(row)

    with open(os.path.join(OUT, "demo_school_3fish.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        for hrow in header:
            w.writerow(hrow)
        w.writerows(rows)
    print(f"  wrote {os.path.join(OUT, 'demo_school_3fish.csv')}  ({len(rows)} rows x {len(header[3])} cols)")


# ---------------------------------------------------------------------------
# 5. Batch Processing tab  --  Format B (named), 3D, tokenised filenames
# ---------------------------------------------------------------------------
def gen_batch():
    print("Batch Processing tab (Format B named, 3D):")
    midline = ["snout", "mid1", "mid2", "mid3", "mid4", "peduncle", "caudaltip"]
    s_vals = [0.0, 0.16, 0.33, 0.5, 0.66, 0.83, 1.0]
    cols = []
    for n in midline:
        cols += [f"{n}_X", f"{n}_Y", f"{n}_Z"]

    # (flow speed token, tail freq, ground speed) -> fish swims against flow
    trials = [
        (0.5, 4.5, 1.0),
        (1.0, 5.0, 1.2),
        (1.5, 5.5, 1.5),
        (2.0, 6.0, 1.8),
        (2.5, 6.5, 2.0),
    ]
    n_frames = 400

    for i, (flow, freq, Ubl) in enumerate(trials, start=1):
        fname = os.path.join(BATCH, f"WSBS_Swim_{flow}BL_Fish{i:02d}xyzpts.csv")
        U = Ubl * L
        rows = []
        for f in range(n_frames):
            t = f / FPS
            cx, cy = swim_center(t, U, HEADING, C0)
            row = []
            for nm, si in zip(midline, s_vals):
                xb = -si * L
                yb = amp_fn(si) * math.sin(2 * math.pi * (freq * t - si / LAMBDA))
                zb = L * 0.03 * math.sin(2 * math.pi * (freq * t - si / LAMBDA) + 0.7)
                x, y = rot_xy(xb, yb, HEADING)
                x, y, z = add_jitter(cx + x, cy + y, zb)
                row += [fmt(x), fmt(y), fmt(z)]
            rows.append(row)
        write_csv(fname, cols, rows)


if __name__ == "__main__":
    print("Generating Kinemetrix demo data...\n")
    gen_kinematics()
    gen_fin()
    gen_pectoral(0.0, os.path.join(OUT, "demo_pectoral_inphase.csv"))
    gen_pectoral(math.pi, os.path.join(OUT, "demo_pectoral_antiphase.csv"))
    gen_school()
    gen_batch()
    print("\nDone.")
