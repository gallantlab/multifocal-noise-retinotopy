#!/usr/bin/env python3
"""Render the multifocal m-sequence design matrix to design_matrix.png.

Reads msequence_design.json (written by make_msequence_movie.py) and draws a
data-only heatmap: one row per m-sequence state (top = first state -> bottom =
last), one column per wedge regressor (8 columns). A cell is "on" (red) when
that wedge is red during that state. No axes/margins, so the image's vertical
extent maps linearly to time, letting the web page overlay an aligned cursor.
"""

import os
import json
import numpy as np
from PIL import Image

HERE = os.path.dirname(os.path.abspath(__file__))

COL_W = 44          # px per regressor column
DIVIDER = 1         # px between columns
ROW_H = 8           # px per state row (height = L * ROW_H)
ON = np.array([200, 40, 40], dtype=np.uint8)        # wedge on
OFF = np.array([245, 245, 245], dtype=np.uint8)     # wedge off
LINE = np.array([170, 170, 170], dtype=np.uint8)    # dividers


def main():
    with open(os.path.join(HERE, "msequence_design.json")) as f:
        d = json.load(f)
    design = np.array(d["design"], dtype=bool)   # (L, n_wedges)
    L, n_wedges = design.shape

    width = n_wedges * COL_W + (n_wedges + 1) * DIVIDER
    height = L * ROW_H
    img = np.empty((height, width, 3), dtype=np.uint8)
    img[:] = LINE

    for k in range(n_wedges):
        x0 = DIVIDER + k * (COL_W + DIVIDER)
        col = design[:, k]                                   # (L,)
        cells = np.where(col[:, None], ON, OFF)              # (L, 3)
        col_px = np.repeat(cells, ROW_H, axis=0)             # (height, 3)
        img[:, x0:x0 + COL_W, :] = col_px[:, None, :]

    out = os.path.join(HERE, "design_matrix.png")
    Image.fromarray(img, "RGB").save(out)
    print(f"Wrote {out}  ({width}x{height}, {L} states x {n_wedges} regressors)")


if __name__ == "__main__":
    main()
