#!/usr/bin/env python3
"""Generate a retinotopic polar-angle wedge movie as a PNG sequence.

512x512, 30 Hz. A 45-degree red pie wedge (from a diameter-512 pie split into
8 wedges) is shown for 4.1 s at each of 8 angular positions, sweeping CCW with
wedge 1 spanning 0-45 degrees (0 deg = 3 o'clock). Frames are written to
frames/frame_NNNNN.png on a [128,128,128] gray background.
"""

import os
import numpy as np
from PIL import Image

# --- Parameters ---
SIZE = 512                      # frame width/height in pixels
RADIUS = SIZE / 2              # pie radius (diameter 512)
N_WEDGES = 8                   # number of pie wedges
WEDGE_DEG = 360.0 / N_WEDGES   # 45 degrees per wedge
FPS = 30
SECONDS_PER_WEDGE = 4.1
FRAMES_PER_WEDGE = round(FPS * SECONDS_PER_WEDGE)  # 123
BG = np.array([128, 128, 128], dtype=np.uint8)
RED = np.array([255, 0, 0], dtype=np.uint8)
OUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "frames")


def main():
    os.makedirs(OUT_DIR, exist_ok=True)

    # Pixel-center coordinates so the disc is symmetric about the image center.
    cx = cy = SIZE / 2.0
    xs = np.arange(SIZE) + 0.5
    ys = np.arange(SIZE) + 0.5
    gx, gy = np.meshgrid(xs - cx, ys - cy)   # gx, gy shape (SIZE, SIZE)

    radius = np.hypot(gx, gy)
    disc = radius <= RADIUS

    # Angle 0 deg at 3 o'clock, increasing CCW (upward on screen). Image rows
    # increase downward, so negate gy. Map to [0, 360).
    angle = np.degrees(np.arctan2(-gy, gx)) % 360.0

    frame_index = 0
    for w in range(N_WEDGES):
        lo = w * WEDGE_DEG
        hi = (w + 1) * WEDGE_DEG
        wedge_mask = disc & (angle >= lo) & (angle < hi)

        frame = np.empty((SIZE, SIZE, 3), dtype=np.uint8)
        frame[:] = BG
        frame[wedge_mask] = RED

        img = Image.fromarray(frame, mode="RGB")
        for _ in range(FRAMES_PER_WEDGE):
            img.save(os.path.join(OUT_DIR, f"frame_{frame_index:05d}.png"))
            frame_index += 1

    print(f"Wrote {frame_index} frames to {OUT_DIR}")
    print(f"  {N_WEDGES} wedges x {FRAMES_PER_WEDGE} frames/wedge")


if __name__ == "__main__":
    main()
