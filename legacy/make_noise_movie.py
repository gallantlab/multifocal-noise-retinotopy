#!/usr/bin/env python3
"""Fill the multifocal m-sequence wedges with bandpass 1/f spatiotemporal noise.

Each m-sequence state gets an independent 512x512x123 block of noise whose
amplitude spectrum is 1/f in BOTH space and time, bandpassed to:
  - spatial:  4 .. 128 cyc/image (isotropic, radial)
  - temporal: ~0.25 .. 15 cyc/sec (15 = Nyquist at 30 Hz; 0.244 = 1 cyc / 4.1 s window)
The noise (grayscale, zero-mean about 128) is shown only inside the wedges that
the m-sequence turns on for that state; everything else stays gray (128).

Set N_STATES_TO_GENERATE to control how many states to render (short test runs
vs the full L=63). Also writes temporal/spatial spectrum plots and noise_meta.js
for the web viewer. Reuses msequence_design.json (the m-seq design) and the
existing solid wedge_*.png references unchanged.
"""

import os
import json
import numpy as np
from PIL import Image
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FRAME_DIR = os.path.join(HERE, "frames")

# ---- knobs ----
N_STATES_TO_GENERATE = 5          # bump to 63 for the full movie
SIZE = 512
FPS = 30
FRAMES_PER_STATE = 123
SF_LO, SF_HI = 4.0, 128.0         # spatial band, cyc/image
ORIENT_TAPS = [6, 4, 3, 1]        # 2nd m-seq (x^6+x^4+x^3+x+1): per-wedge orientation
TF_HI = 15.0                      # temporal high (Nyquist at 30 Hz)
TF_LO = FPS / FRAMES_PER_STATE    # ~0.2439 cyc/sec = 1 cycle per state window
GAIN = 38.0                       # gray-levels per std (controls contrast/clipping)
FADE_FRAMES = 10                  # contrast fade-in / fade-out at each state boundary
BASE_SEED = 1234

BG = 128
RADIUS = SIZE / 2.0
N_WEDGES = 8
WEDGE_DEG = 360.0 / N_WEDGES


def m_sequence(n, taps, seed=1):
    """Fibonacci LFSR -> binary maximum-length sequence of length 2^n - 1."""
    period = (1 << n) - 1
    state, seq = seed, []
    for _ in range(period):
        seq.append(state & 1)
        fb = 0
        for t in taps:
            fb ^= (state >> (t - 1)) & 1
        state = (state >> 1) | (fb << (n - 1))
    assert sum(seq) == (1 << (n - 1)), "not a maximum-length sequence"
    return seq


def orientation_mask(FX, FY, orientation):
    """Infinitely narrow orientation band: keep a single Fourier line.

    Vertical image orientation (vertical bars vary horizontally) -> energy on the
    fy=0 line. Horizontal -> fx=0 line. Zero everywhere else => the carrier
    contains exactly one orientation and nothing else.
    """
    if orientation == "vertical":
        return (FY == 0)
    if orientation == "horizontal":
        return (FX == 0)
    raise ValueError(f"unknown orientation {orientation!r}")


def build_filter(orientation):
    """Separable 1/f bandpass amplitude filter matching rfftn(axes=(0,1,2))."""
    fx = np.fft.fftfreq(SIZE, d=1.0) * SIZE          # cyc/image, axis 1
    fy = np.fft.fftfreq(SIZE, d=1.0) * SIZE          # cyc/image, axis 0
    FX, FY = np.meshgrid(fx, fy)
    fr = np.hypot(FX, FY)
    with np.errstate(divide="ignore"):
        a_sp = np.where((fr >= SF_LO) & (fr <= SF_HI), 1.0 / fr, 0.0)  # (512,512)
    a_sp *= orientation_mask(FX, FY, orientation)

    ft = np.fft.rfftfreq(FRAMES_PER_STATE, d=1.0 / FPS)               # cyc/sec, axis 2
    with np.errstate(divide="ignore"):
        a_t = np.where((ft >= TF_LO * 0.999) & (ft <= TF_HI), 1.0 / ft, 0.0)
    a_t[0] = 0.0                                                       # kill DC

    return a_sp[:, :, None] * a_t[None, None, :]                       # (512,512,62)


def wedge_masks():
    cx = cy = SIZE / 2.0
    c = np.arange(SIZE) + 0.5
    gx, gy = np.meshgrid(c - cx, c - cy)
    disc = np.hypot(gx, gy) <= RADIUS
    angle = np.degrees(np.arctan2(-gy, gx)) % 360.0
    return [disc & (angle >= k * WEDGE_DEG) & (angle < (k + 1) * WEDGE_DEG)
            for k in range(N_WEDGES)]


def make_noise(seed, H, want_lum=False):
    """One state's colored bandpass 1/f noise.

    Each RGB channel is an INDEPENDENT noise field pushed through the given
    1/f orientation filter, so every Fourier component carries an independent random
    [r,g,b] -> randomly multicolored, while each channel keeps the exact spectrum.
    Returns (rgb_uint8 (512,512,123,3), luminance_float or None).
    """
    rgb = np.empty((SIZE, SIZE, FRAMES_PER_STATE, 3), dtype=np.uint8)
    lum = np.zeros((SIZE, SIZE, FRAMES_PER_STATE)) if want_lum else None
    for c in range(3):
        rng = np.random.default_rng(seed * 3 + c)
        w = rng.standard_normal((SIZE, SIZE, FRAMES_PER_STATE))
        spec = np.fft.rfftn(w, axes=(0, 1, 2)) * H
        noise = np.fft.irfftn(spec, s=(SIZE, SIZE, FRAMES_PER_STATE), axes=(0, 1, 2))
        noise /= noise.std()
        rgb[..., c] = np.clip(BG + GAIN * noise, 0, 255).astype(np.uint8)
        if want_lum:
            lum += noise / 3.0
    return rgb, lum


STYLE = dict(fg="#e8eaed", muted="#9aa0a8", panel="#24272d",
             axbg="#1a1c20", data="#2bd1ff", ideal="#e02828", band="#e0282822")


def _dark(ax):
    ax.set_facecolor(STYLE["axbg"])
    for s in ax.spines.values():
        s.set_color(STYLE["muted"])
    ax.tick_params(colors=STYLE["muted"])
    ax.xaxis.label.set_color(STYLE["fg"]); ax.yaxis.label.set_color(STYLE["fg"])
    ax.title.set_color(STYLE["fg"])
    ax.grid(True, which="both", color="#ffffff12", lw=0.6)


def plot_temporal(noise):
    ft = np.fft.rfftfreq(FRAMES_PER_STATE, d=1.0 / FPS)
    power = (np.abs(np.fft.rfft(noise, axis=2)) ** 2).mean(axis=(0, 1))
    band = (ft >= TF_LO * 0.999) & (ft <= TF_HI) & (ft > 0)
    K = np.median(power[band] * ft[band] ** 2)                  # ideal: power = K/f^2

    fig, ax = plt.subplots(figsize=(4.2, 3.3), dpi=130)
    fig.patch.set_facecolor(STYLE["panel"])
    ax.axvspan(TF_LO, TF_HI, color=STYLE["ideal"], alpha=0.10)
    ax.loglog(ft[1:], power[1:], color=STYLE["data"], lw=1.6, label="measured")
    ax.loglog(ft[band], K / ft[band] ** 2, "--", color=STYLE["ideal"], lw=1.4,
              label="ideal 1/f²")
    for f in (TF_LO, TF_HI):
        ax.axvline(f, color=STYLE["muted"], lw=0.8, ls=":")
    ax.set_xlabel("temporal frequency (cyc/s)"); ax.set_ylabel("power")
    ax.set_title("Temporal spectrum"); _dark(ax)
    ax.legend(facecolor=STYLE["panel"], edgecolor=STYLE["muted"],
              labelcolor=STYLE["fg"], fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "temporal_spectrum.png"))
    plt.close(fig)


def time_avg_power2d(noise):
    p2d = np.zeros((SIZE, SIZE))
    for t in range(FRAMES_PER_STATE):
        p2d += np.abs(np.fft.fft2(noise[:, :, t])) ** 2
    return p2d / FRAMES_PER_STATE


def plot_spatial(noise):
    p2d = time_avg_power2d(noise)
    f = np.fft.fftfreq(SIZE, d=1.0) * SIZE                  # cyc/image
    # Profile along fy=0 (the vertical-orientation axis); the horizontal carrier
    # has an identical radial 1/f profile on the fx=0 axis.
    line = p2d[0, :]
    pos = f > 0
    fpos, ppos = f[pos], line[pos]
    band = (fpos >= SF_LO) & (fpos <= SF_HI)
    K = np.median(ppos[band] * fpos[band] ** 2)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.3), dpi=130)
    fig.patch.set_facecolor(STYLE["panel"])
    ax1.axvspan(SF_LO, SF_HI, color=STYLE["ideal"], alpha=0.10)
    ax1.loglog(fpos, np.maximum(ppos, 1e-30), color=STYLE["data"], lw=1.6, label="measured")
    ax1.loglog(fpos[band], K / fpos[band] ** 2, "--", color=STYLE["ideal"], lw=1.4,
               label="ideal 1/f²")
    for v in (SF_LO, SF_HI):
        ax1.axvline(v, color=STYLE["muted"], lw=0.8, ls=":")
    ax1.set_ylim(K / SF_HI ** 2 / 1e3, K / SF_LO ** 2 * 1e2)
    ax1.set_xlabel("spatial frequency (cyc/image)"); ax1.set_ylabel("power")
    ax1.set_title("Spatial spectrum (along vertical axis)"); _dark(ax1)
    ax1.legend(facecolor=STYLE["panel"], edgecolor=STYLE["muted"],
               labelcolor=STYLE["fg"], fontsize=8)

    shown = np.log10(np.fft.fftshift(p2d) + 1e-9)
    ext = [-SIZE // 2, SIZE // 2, -SIZE // 2, SIZE // 2]
    ax2.imshow(shown, extent=ext, cmap="magma", origin="lower")
    ax2.set_xlim(-160, 160); ax2.set_ylim(-160, 160)
    ax2.set_xlabel("fx (cyc/image)"); ax2.set_ylabel("fy (cyc/image)")
    ax2.set_title("2D power (log) — cross = vertical + horizontal"); _dark(ax2)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "spatial_spectrum.png"))
    plt.close(fig)


def plot_orientation(noise):
    """Angular power distribution -> verify all energy is at one orientation."""
    p2d = time_avg_power2d(noise)
    f = np.fft.fftfreq(SIZE, d=1.0) * SIZE
    FX, FY = np.meshgrid(f, f)
    fr = np.hypot(FX, FY)
    band = (fr >= SF_LO) & (fr <= SF_HI)
    # image orientation = Fourier wavevector angle + 90 deg, folded to [0,180)
    orient = (np.degrees(np.arctan2(FY, FX)) + 90.0) % 180.0
    obins = np.arange(0, 181, 1.0)
    idx = np.clip(np.digitize(orient[band], obins) - 1, 0, len(obins) - 2)
    power = np.zeros(len(obins) - 1)
    np.add.at(power, idx, p2d[band])
    centers = obins[:-1] + 0.5
    power = power / power.max()

    fig, ax = plt.subplots(figsize=(4.2, 3.3), dpi=130)
    fig.patch.set_facecolor(STYLE["panel"])
    ax.fill_between(centers, 0, power, color=STYLE["data"], alpha=0.5)
    ax.plot(centers, power, color=STYLE["data"], lw=1.6)
    ax.axvline(0, color="#ffb028", lw=1.0, ls="--", label="horizontal (0°)")
    ax.axvline(90, color=STYLE["ideal"], lw=1.0, ls="--", label="vertical (90°)")
    ax.axvline(180, color="#ffb028", lw=1.0, ls="--")
    ax.set_xlim(0, 180); ax.set_xticks([0, 45, 90, 135, 180])
    ax.set_xlabel("image orientation (deg)  0=horiz, 90=vert")
    ax.set_ylabel("normalized power")
    ax.set_title("Orientation spectrum"); _dark(ax)
    ax.legend(facecolor=STYLE["panel"], edgecolor=STYLE["muted"],
              labelcolor=STYLE["fg"], fontsize=8)
    fig.tight_layout(); fig.savefig(os.path.join(HERE, "orientation_spectrum.png"))
    plt.close(fig)


def render_orientation_matrix(orient, design, out):
    """L x 8 image: vertical=cyan, horizontal=orange; off-wedges dimmed."""
    COL_W, DIV, ROW_H = 44, 1, 8
    VERT = np.array([43, 150, 255]); HORIZ = np.array([255, 150, 40])
    LINE = np.array([170, 170, 170]); DIM = 0.32
    L = orient.shape[0]
    width = N_WEDGES * COL_W + (N_WEDGES + 1) * DIV
    img = np.empty((L * ROW_H, width, 3), dtype=np.uint8)
    img[:] = LINE
    for k in range(N_WEDGES):
        x0 = DIV + k * (COL_W + DIV)
        cells = np.where(orient[:, k:k+1] == 0, VERT, HORIZ).astype(float)   # (L,3)
        off = design[:, k] == 0
        cells[off] = BG + (cells[off] - BG) * DIM                            # fade off-wedges
        col = np.repeat(cells.astype(np.uint8), ROW_H, axis=0)               # (L*ROW_H,3)
        img[:, x0:x0 + COL_W, :] = col[:, None, :]
    Image.fromarray(img, "RGB").save(out)


def main():
    with open(os.path.join(HERE, "msequence_design.json")) as f:
        design = np.array(json.load(f)["design"], dtype=int)   # (L, 8) on/off
    L = design.shape[0]
    n_gen = min(N_STATES_TO_GENERATE, L)

    # Second, independent m-sequence -> per-wedge orientation (0=vertical, 1=horizontal).
    oseq = m_sequence(6, ORIENT_TAPS)
    assert len(oseq) == L
    oshifts = [round(k * L / N_WEDGES) for k in range(N_WEDGES)]
    orient = np.array([[oseq[(s + oshifts[k]) % L] for k in range(N_WEDGES)]
                       for s in range(L)], dtype=int)           # (L, 8)

    H_vert = build_filter("vertical")
    H_horiz = build_filter("horizontal")
    masks = wedge_masks()

    # Per-state contrast envelope: fade in from gray, hold, fade out to gray.
    env = np.ones(FRAMES_PER_STATE)
    F = min(FADE_FRAMES, FRAMES_PER_STATE // 2)
    ramp = np.linspace(0.0, 1.0, F, endpoint=False)   # 0 .. (F-1)/F
    env[:F] = ramp
    env[FRAMES_PER_STATE - F:] = ramp[::-1]

    os.makedirs(FRAME_DIR, exist_ok=True)
    for f in os.listdir(FRAME_DIR):
        if f.endswith(".png"):
            os.remove(os.path.join(FRAME_DIR, f))

    spec_v = spec_h = None      # one V and one H carrier from state 0, for the spectra
    frame_idx = 0
    for s in range(n_gen):
        on = [k for k in range(N_WEDGES) if design[s, k]]
        buf = np.full((SIZE, SIZE, FRAMES_PER_STATE, 3), BG, dtype=np.uint8)
        for k in on:
            is_v = orient[s, k] == 0
            Hk = H_vert if is_v else H_horiz
            need = (s == 0) and ((is_v and spec_v is None) or (not is_v and spec_h is None))
            # Unique seed per (state, wedge) => every wedge is INDEPENDENT noise,
            # so same-orientation neighbors show a discontinuity at the boundary.
            rgb, lum = make_noise(BASE_SEED + s * N_WEDGES + k, Hk, want_lum=need)
            if need:
                if is_v: spec_v = lum
                else: spec_h = lum
            buf[masks[k]] = rgb[masks[k]]
        for t in range(FRAMES_PER_STATE):
            fr = buf[:, :, t, :].astype(np.int16)
            out = (BG + env[t] * (fr - BG)).clip(0, 255).astype(np.uint8)
            Image.fromarray(out, "RGB").save(
                os.path.join(FRAME_DIR, f"frame_{frame_idx:05d}.png"))
            frame_idx += 1
        print(f"  state {s}: " + (", ".join(
            f"{k+1}{'V' if orient[s,k]==0 else 'H'}" for k in on) or "blank"))

    # combine a V and an H carrier so the orientation plot shows both spikes
    spectra_noise = (spec_v if spec_v is not None else 0) + \
                    (spec_h if spec_h is not None else 0)

    plot_temporal(spectra_noise)
    plot_spatial(spectra_noise)
    plot_orientation(spectra_noise)
    render_orientation_matrix(orient, design, os.path.join(HERE, "orientation_matrix.png"))

    meta = {
        "generated_states": n_gen, "total_states": L,
        "frames_per_state": FRAMES_PER_STATE, "fps": FPS,
        "sf_band": [SF_LO, SF_HI], "tf_band": [round(TF_LO, 3), TF_HI],
        "color": "random rgb per Fourier component", "gain": GAIN,
        "fade_frames": F,
        "orient_taps": ORIENT_TAPS, "orient_shifts": oshifts,
        "orient_map": {"0": "vertical", "1": "horizontal"},
        "orient_design": orient.tolist(),
    }
    with open(os.path.join(HERE, "noise_meta.js"), "w") as f:
        f.write("window.NOISE = " + json.dumps(meta) + ";\n")

    print(f"Generated {frame_idx} noise frames ({n_gen}/{L} states); per-wedge V/H m-seq taps={ORIENT_TAPS}.")
    print("Wrote spectra, orientation_matrix.png, noise_meta.js")


if __name__ == "__main__":
    main()
