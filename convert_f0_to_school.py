#!/usr/bin/env python3
"""
Convert a wide reshaped DeepLabCut export into the FishN_Pk_<coord> convention
that Kinemetrix's School Metrics tab accepts via load_fish_points.

The coworker's F0 file uses columns like:  index, ind1_snout_x, ind1_snout_y,
ind1_snout_z, ind1_midline1_x, ...  (10 individuals x 10 bodyparts x 3 coords).
load_fish_points only recognises  frame + FishN_Pk_[xyz]  (or ptN_[XY]),
so this script renames:
    index            -> frame
    ind<N>_<bp>_<c>  -> Fish<N>_P<k>_<c>

Bodypart -> point-number mapping (order in the source header):
    snout           -> P1     caudal_fin_tip  -> P6
    midline1        -> P2     pect_fin_tip_R  -> P7
    midline2        -> P3     pect_fin_tip_L  -> P8
    midline3        -> P4     pect_fin_base_R -> P9
    peduncle        -> P5     pect_fin_base_L -> P10

NA values are kept as-is (readtable treats them as NaN).
After converting, open the School Metrics tab and pick snout=P1, peduncle=P5.

Run:  python convert_f0_to_school.py
"""

import csv
import os
import re

HERE = os.path.dirname(os.path.abspath(__file__))
SRC = os.path.join(HERE, "Data", "F0_7.4_3_0.5BLS_15_AUG_2026_000000_DLC_3D.csv")
DST = os.path.join(HERE, "Data", "F0_7.4_3_0.5BLS_15_AUG_2026_000000_DLC_3D_school.csv")

# bodyparts in the exact order they appear per individual in the source header
BODYPARTS = [
    "snout",
    "midline1",
    "midline2",
    "midline3",
    "peduncle",
    "caudal_fin_tip",
    "pect_fin_tip_R",
    "pect_fin_tip_L",
    "pect_fin_base_R",
    "pect_fin_base_L",
]

COORDS = ["x", "y", "z"]
FRAME_SRC = "index"

IND_RE = re.compile(r"^ind(\d+)_")


def main():
    with open(SRC, newline="") as fh:
        reader = csv.reader(fh)
        header = next(reader)
        rows = list(reader)

    col_index = {name: i for i, name in enumerate(header)}

    # discover individuals present in the header (sorted numerically)
    ind_numbers = sorted(
        {int(m.group(1)) for h in header for m in [IND_RE.match(h)] if m}
    )

    out_header = ["frame"]
    for n in ind_numbers:
        for k, _bp in enumerate(BODYPARTS, start=1):
            for c in COORDS:
                out_header.append(f"Fish{n}_P{k}_{c}")

    missing_src = []
    for n in ind_numbers:
        for k, bp in enumerate(BODYPARTS, start=1):
            for c in COORDS:
                src = f"ind{n}_{bp}_{c}"
                if src not in col_index:
                    missing_src.append(src)
    if missing_src:
        raise SystemExit(f"Missing source columns: {missing_src[:5]} ...")

    def num(v):
        # readtable treats 'NA' as text (-> cell column); 'NaN' is a numeric
        # token, so writing it keeps every coordinate column a double.
        return "NaN" if v.strip().upper() in ("NA", "N/A") else v

    out_rows = []
    for row in rows:
        out = [row[col_index[FRAME_SRC]]]
        for n in ind_numbers:
            for k, bp in enumerate(BODYPARTS, start=1):
                for c in COORDS:
                    out.append(num(row[col_index[f"ind{n}_{bp}_{c}"]]))
        out_rows.append(out)

    with open(DST, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(out_header)
        w.writerows(out_rows)

    print(f"Converted {len(rows)} rows -> {len(out_header)} columns")
    print(f"  individuals: ind{ind_numbers[0]} .. ind{ind_numbers[-1]}")
    print(f"  output: {DST}")
    print("  School Metrics tab: snout=P1, peduncle=P5")


if __name__ == "__main__":
    main()
