#!/usr/bin/env python3
"""Multifocal m-sequence retinotopy stimulus as a PNG sequence.

A single binary maximum-length sequence (m-sequence, length L = 2^n - 1) drives
all 8 pie wedges: wedge k uses a distinct circular shift of the sequence. In
each state a wedge is red when its shifted m-sequence value is 1, gray when 0;
several wedges can be on at once, and all-zero states render as a blank gray
frame. Each state is held for 4.1 s (123 frames @ 30 Hz).

Outputs:
  frames/frame_NNNNN.png   - the movie (L * 123 frames), 512x512, RGB
  wedge_1.png .. wedge_8.png - single-wedge reference images (for the web montage)
  msequence_design.json    - design array + metadata (for the matrix renderer)
  msequence_design.js      - same payload as window.MSEQ (for the web page, file:// safe)
"""

import os
import json
import numpy as np
from PIL import Image

# --- Stimulus geometry (consistent with make_movie.py) ---
SIZE = 512
RADIUS = SIZE / 2.0
N_WEDGES = 8
WEDGE_DEG = 360.0 / N_WEDGES
BG = np.array([128, 128, 128], dtype=np.uint8)
RED = np.array([255, 0, 0], dtype=np.uint8)

# --- m-sequence / timing ---
N = 6                                   # LFSR order -> L = 2^N - 1 = 63 states
TAPS = [6, 1]                           # primitive poly x^6 + x + 1
FPS = 30
SECONDS_PER_STATE = 4.1
FRAMES_PER_STATE = round(FPS * SECONDS_PER_STATE)   # 123

HERE = os.path.dirname(os.path.abspath(__file__))
FRAME_DIR = os.path.join(HERE, "frames")


def m_sequence(n, taps, seed=1):
    """Fibonacci LFSR -> binary m-sequence of length 2^n - 1 (verified maximal)."""
    period = (1 << n) - 1
    state = seed
    seq = []
    for _ in range(period):
        seq.append(state & 1)
        fb = 0
        for t in taps:
            fb ^= (state >> (t - 1)) & 1
        state = (state >> 1) | (fb << (n - 1))
    # Verify it is a true maximum-length sequence.
    assert sum(seq) == (1 << (n - 1)), f"expected {1 << (n-1)} ones, got {sum(seq)}"
    assert len(seq) == period
    return seq


def wedge_masks():
    """Boolean pixel mask for each of the 8 wedges (disc & angular sector)."""
    cx = cy = SIZE / 2.0
    coords = np.arange(SIZE) + 0.5
    gx, gy = np.meshgrid(coords - cx, coords - cy)
    disc = np.hypot(gx, gy) <= RADIUS
    angle = np.degrees(np.arctan2(-gy, gx)) % 360.0   # 0 deg at 3 o'clock, CCW
    masks = []
    for k in range(N_WEDGES):
        lo, hi = k * WEDGE_DEG, (k + 1) * WEDGE_DEG
        masks.append(disc & (angle >= lo) & (angle < hi))
    return masks


def main():
    os.makedirs(FRAME_DIR, exist_ok=True)
    seq = m_sequence(N, TAPS)
    L = len(seq)

    # Distinct, ~evenly spaced circular shifts -> near-orthogonal regressors.
    shifts = [round(k * L / N_WEDGES) for k in range(N_WEDGES)]
    assert len(set(shifts)) == N_WEDGES, f"shifts not distinct: {shifts}"

    # design[s, k] = 1 when wedge k is on during state s
    design = np.array([[seq[(s + shifts[k]) % L] for k in range(N_WEDGES)]
                       for s in range(L)], dtype=int)

    masks = wedge_masks()

    # Single-wedge reference images for the montage.
    for k, m in enumerate(masks):
        ref = np.empty((SIZE, SIZE, 3), dtype=np.uint8)
        ref[:] = BG
        ref[m] = RED
        Image.fromarray(ref, "RGB").save(os.path.join(HERE, f"wedge_{k+1}.png"))

    # Movie frames: union of all wedges that are on in each state.
    frame_idx = 0
    for s in range(L):
        frame = np.empty((SIZE, SIZE, 3), dtype=np.uint8)
        frame[:] = BG
        for k in range(N_WEDGES):
            if design[s, k]:
                frame[masks[k]] = RED
        img = Image.fromarray(frame, "RGB")
        for _ in range(FRAMES_PER_STATE):
            img.save(os.path.join(FRAME_DIR, f"frame_{frame_idx:05d}.png"))
            frame_idx += 1

    payload = {
        "n": N, "taps": TAPS, "L": L, "n_wedges": N_WEDGES,
        "fps": FPS, "seconds_per_state": SECONDS_PER_STATE,
        "frames_per_state": FRAMES_PER_STATE, "total_frames": frame_idx,
        "shifts": shifts,
        "design": design.tolist(),
    }
    with open(os.path.join(HERE, "msequence_design.json"), "w") as f:
        json.dump(payload, f)
    with open(os.path.join(HERE, "msequence_design.js"), "w") as f:
        f.write("window.MSEQ = " + json.dumps(payload) + ";\n")

    print(f"m-sequence n={N} taps={TAPS} -> L={L} states, shifts={shifts}")
    print(f"Wrote {frame_idx} frames ({L} states x {FRAMES_PER_STATE}) to {FRAME_DIR}")
    print(f"Wrote 8 wedge refs, msequence_design.json, msequence_design.js")


if __name__ == "__main__":
    main()
