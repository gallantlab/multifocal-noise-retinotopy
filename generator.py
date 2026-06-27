#!/usr/bin/env python3
"""Generate multifocal m-sequence retinotopy noise movies.

The stimulus is a disc split into ``N`` regions. In **wedges** geometry the disc
is split into equal angular sectors (polar-angle mapping); in **rings** geometry
it is split into ``N`` concentric annuli (eccentricity mapping). A binary
maximum-length sequence (m-sequence) turns each region on/off across 63 *states*;
a second, independent m-sequence sets each on-region's noise *orientation*. Each
on-region is filled with bandpass 1/f spatiotemporal noise (multicolored or
grayscale), confined to a single orientation, and independent from every other
region. States fade in/out to the background, and the movie can be padded with
full-screen isotropic noise.

The single public entry point is :func:`generate_movie`, which writes a PNG
sequence to ``frames/`` plus spectra plots, design/orientation matrices, and
``movie_meta.json`` for the web viewer. See ``README.md`` for the experimental
caveats (temporal Nyquist, off-axis orientation bandwidth, K-ary balance, etc.).

Conventions
-----------
- Image arrays are ``(row=y, col=x, ...)`` ``uint8`` RGB unless noted.
- Spatial frequency is in **cyc/width** (= integer FFT index for a width-px field).
- Temporal frequency is in **Hz**; image orientation in **degrees** (0=horizontal
  bars, 90=vertical bars).
- "1/f" means amplitude ∝ 1/f (power ∝ 1/f²).
"""

from __future__ import annotations

import json
import math
import os
from collections.abc import Callable

import matplotlib
matplotlib.use("Agg")
import matplotlib.cm as cm
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

# --- paths ---
HERE = os.path.dirname(os.path.abspath(__file__))
FRAME_DIR = os.path.join(HERE, "frames")

# --- fixed design constants ---
BG = 128                       # mid-gray background / zero-contrast level
GAIN = 38.0                    # gray-levels per noise std (sets contrast; ~0.05% clip at 1/f)
BASE_SEED = 1234               # base RNG seed (fixed -> reproducible movies)
LFSR_N = 6                     # m-sequence order -> L = 2**6 - 1 = 63 states
ONOFF_TAPS = [6, 1]            # primitive poly x^6+x+1     : on/off m-sequence
ORIENT_TAPS = [6, 4, 3, 1]     # primitive poly x^6+x^4+x^3+x+1: orientation m-sequence
DEMO_STATES = 5                # states rendered in "demo" mode
WEDGE_HALFWIDTH = 4.0          # angular half-width (deg) of off-axis orientation wedges
RING_FOVEA_FRAC = 0.05         # foveal radius offset (frac of max radius) for log ring spacing
RING_SPACINGS = ("log", "equal_width", "equal_area")

#: Default parameters (mirror INTERFACE.txt). ``generate_movie`` fills missing keys.
#: The ``n_wedges`` / ``wedge_*`` keys are generic *region* parameters that apply to
#: both geometries (a "region" is an angular wedge or a concentric ring); the names
#: are kept for backward compatibility with previously saved movies and the UI.
DEFAULTS = {
    "geometry": "wedges",                                # wedges | rings
    "width": 512, "n_wedges": 8, "wedge_rotation": 22.5, "wedge_sec": 4.1, "fps": 30,
    "ring_spacing": "log",                               # log | equal_width | equal_area (rings only)
    "color": "color",                                    # color | bw
    "tf_shape": "1/f", "tf_lo": 0.5, "tf_hi": 15.0,      # 1/f | flat
    "sf_shape": "1/f", "sf_lo": 2.0, "sf_hi": 128.0,     # cyc/width (cycles across the movie width)
    "n_orientations": 2,
    "background": "gray",                                # gray | random | oriented
    "fade_frames": 5, "pad_sec": 2.0,
    "fixation": "on",                                    # off | on (central 2x2 spot)
    "mode": "demo",                                      # demo | full
}

FIX_BLOCK_SEC = 0.5            # fixation spot re-colors every this many seconds
FIX_SIZE = 2                   # fixation spot side length (pixels)


class Cancelled(Exception):
    """Raised inside :func:`generate_movie` when the UI requests cancellation."""


# ============================================================ m-sequence design
def m_sequence(order: int, taps: list[int], seed: int = 1) -> list[int]:
    """Binary maximum-length sequence via a Fibonacci LFSR.

    Parameters
    ----------
    order : int
        Register length ``n``; the sequence has period ``2**n - 1``.
    taps : list[int]
        1-indexed feedback tap positions of a primitive polynomial.
    seed : int
        Non-zero initial register state.

    Returns
    -------
    list[int]
        The 0/1 sequence of length ``2**order - 1``. Asserts maximal length
        (exactly ``2**(order-1)`` ones).
    """
    period = (1 << order) - 1
    state, seq = seed, []
    for _ in range(period):
        seq.append(state & 1)
        feedback = 0
        for t in taps:
            feedback ^= (state >> (t - 1)) & 1
        state = (state >> 1) | (feedback << (order - 1))
    assert sum(seq) == (1 << (order - 1)), "taps are not a primitive (maximal) polynomial"
    return seq


def build_design(n_regions: int, n_orientations: int
                 ) -> tuple[np.ndarray, np.ndarray, list[float], int]:
    """Build the on/off and orientation designs from two m-sequences.

    Both designs use a distinct circular shift of an m-sequence per region so the
    region regressors are near-orthogonal. The orientation index is decoded from
    ``ceil(log2 K)`` consecutive bits of the orientation sequence (``mod K``) --
    a pragmatic K-ary decoder, exact only for ``K == 2`` (see README).

    Parameters
    ----------
    n_regions, n_orientations : int

    Returns
    -------
    design : np.ndarray, shape (L, n_regions), int
        1 where a region is on in a given state.
    orient_index : np.ndarray, shape (L, n_regions), int in [0, K)
        Orientation index per region per state.
    angles : list[float]
        The ``K`` equally spaced image orientations in degrees (``i*180/K``).
    L : int
        Number of states (``2**LFSR_N - 1`` = 63).
    """
    onoff = m_sequence(LFSR_N, ONOFF_TAPS)
    oseq = m_sequence(LFSR_N, ORIENT_TAPS)
    L = len(onoff)
    shifts = [round(k * L / n_regions) for k in range(n_regions)]
    n_bits = max(1, math.ceil(math.log2(n_orientations))) if n_orientations > 1 else 1

    def orient_at(state: int, region: int) -> int:
        bits = 0
        for j in range(n_bits):
            bits |= oseq[(state + shifts[region] + j) % L] << j
        return bits % n_orientations

    design = np.array([[onoff[(s + shifts[k]) % L] for k in range(n_regions)]
                       for s in range(L)], dtype=int)
    orient_index = np.array([[orient_at(s, k) for k in range(n_regions)]
                             for s in range(L)], dtype=int)
    angles = [i * 180.0 / n_orientations for i in range(n_orientations)]
    return design, orient_index, angles, L


def orthogonal_orientation(present_deg: list[float]) -> float:
    """Orientation (deg, in [0, 180)) maximally orthogonal to all given orientations.

    Orientation is periodic mod 180, so the maximally-orthogonal angle is the
    midpoint of the largest empty arc among ``present_deg`` on that circle -- the
    angle whose *minimum* angular distance to every present orientation is as
    large as possible. With one orientation this is exactly perpendicular
    (``theta + 90``); with several it is the best compromise (e.g. 45 deg when
    both 0 and 90 are present). Empty input returns 0.
    """
    pts = sorted({a % 180.0 for a in present_deg})
    if not pts:
        return 0.0
    if len(pts) == 1:
        return (pts[0] + 90.0) % 180.0
    best_mid, best_gap = 0.0, -1.0
    for i in range(len(pts)):
        lo = pts[i]
        hi = pts[i + 1] if i + 1 < len(pts) else pts[0] + 180.0   # wrap mod 180
        gap = hi - lo
        if gap > best_gap:
            best_gap, best_mid = gap, (lo + gap / 2.0) % 180.0
    return best_mid


# ============================================================ Fourier filters
def orientation_mask(FX: np.ndarray, FY: np.ndarray, theta_deg: float,
                     halfwidth: float) -> np.ndarray:
    """Boolean mask selecting one image orientation in the 2-D Fourier plane.

    Axis-aligned orientations (0 / 90 deg) are kept as the *exact* infinitely
    narrow Fourier line; an off-axis orientation cannot fall on the integer FFT
    grid, so it is approximated by a thin angular wedge of the given half-width.

    Parameters
    ----------
    FX, FY : np.ndarray
        Spatial-frequency grids (cyc/width), from ``np.meshgrid``.
    theta_deg : float
        Target image orientation (0=horizontal bars, 90=vertical bars).
    halfwidth : float
        Angular half-width (deg) used for off-axis orientations.
    """
    t = theta_deg % 180.0
    if abs(t) < 1e-6 or abs(t - 180.0) < 1e-6:
        return FX == 0                       # horizontal bars -> energy on fy axis
    if abs(t - 90.0) < 1e-6:
        return FY == 0                       # vertical bars   -> energy on fx axis
    pixel_orient = (np.degrees(np.arctan2(FY, FX)) + 90.0) % 180.0
    dist = np.abs(((pixel_orient - t + 90.0) % 180.0) - 90.0)
    return dist <= halfwidth


def spatial_filter(width: int, sf_lo: float, sf_hi: float, shape: str,
                   theta: float | None = None,
                   halfwidth: float = WEDGE_HALFWIDTH) -> np.ndarray:
    """2-D spatial amplitude filter (``(width, width)``).

    A radial bandpass in ``[sf_lo, sf_hi]`` cyc/width, with ``1/f`` or ``flat``
    amplitude, optionally restricted to a single orientation ``theta`` (deg).
    ``theta=None`` gives an isotropic filter (all orientations).
    """
    f = np.fft.fftfreq(width, d=1.0) * width          # cyc/width
    FX, FY = np.meshgrid(f, f)
    fr = np.hypot(FX, FY)
    band = (fr >= sf_lo) & (fr <= sf_hi)
    with np.errstate(divide="ignore"):
        amp = (1.0 / fr) if shape == "1/f" else np.ones_like(fr)
    a_sp = np.where(band, amp, 0.0)
    a_sp[~np.isfinite(a_sp)] = 0.0                    # guard fr == 0
    if theta is not None:
        a_sp = a_sp * orientation_mask(FX, FY, theta, halfwidth)
    return a_sp


def temporal_filter(n_frames: int, fps: int, tf_lo: float, tf_hi: float,
                    shape: str) -> np.ndarray:
    """1-D temporal amplitude filter for ``np.fft.rfft`` (length ``n_frames//2+1``).

    Bandpass in ``[tf_lo, tf_hi]`` Hz with ``1/f`` or ``flat`` amplitude; DC is
    zeroed. Note the realized band is limited by the frame rate (Nyquist =
    ``fps/2``) and window (lowest bin = ``fps/n_frames``).
    """
    ft = np.fft.rfftfreq(n_frames, d=1.0 / fps)
    band = (ft >= tf_lo * 0.999) & (ft <= tf_hi)
    with np.errstate(divide="ignore"):
        amp = (1.0 / np.where(ft == 0, 1, ft)) if shape == "1/f" else np.ones_like(ft)
    a_t = np.where(band, amp, 0.0)
    a_t[0] = 0.0
    return a_t


# ============================================================ noise synthesis
def make_noise(seed: int, H: np.ndarray, width: int, n_frames: int, color: str,
               want_lum: bool = False) -> tuple[np.ndarray, np.ndarray | None]:
    """Render one colored/grayscale bandpass noise block.

    White noise is shaped by the separable filter ``H`` in the 3-D Fourier domain
    (``rfftn`` over x, y, t). For ``color`` each RGB channel is an *independent*
    field through the same ``H`` (so each Fourier component carries a random RGB);
    for ``bw`` a single field is replicated to all channels.

    Parameters
    ----------
    seed : int
        Base RNG seed; channel ``c`` uses ``seed*3 + c`` (disjoint streams).
    H : np.ndarray, shape (width, width, n_frames//2+1)
        Separable spatial x temporal amplitude filter.
    color : str
        ``"color"`` (3 independent channels) or ``"bw"`` (1 channel).
    want_lum : bool
        If True, also return the mean-of-channels luminance (for spectra plots).

    Returns
    -------
    rgb : np.ndarray, shape (width, width, n_frames, 3), uint8
    lum : np.ndarray | None, shape (width, width, n_frames), float
    """
    rgb = np.empty((width, width, n_frames, 3), dtype=np.uint8)
    n_channels = 3 if color == "color" else 1
    fields = []
    for c in range(n_channels):
        rng = np.random.default_rng(seed * 3 + c)
        white = rng.standard_normal((width, width, n_frames))
        noise = np.fft.irfftn(np.fft.rfftn(white, axes=(0, 1, 2)) * H,
                              s=(width, width, n_frames), axes=(0, 1, 2))
        std = noise.std()
        noise /= std if std > 0 else 1.0              # guard an all-zero (empty-band) filter
        fields.append(noise)
    for c in range(3):                                # bw replicates the single field
        rgb[..., c] = np.clip(BG + GAIN * fields[c % n_channels], 0, 255).astype(np.uint8)
    lum = np.mean(fields, axis=0) if want_lum else None
    return rgb, lum


# ============================================================ geometry
def _disc_grid(width: int) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Centered coordinate grids and the unit-disc mask for a ``width × width`` field.

    Returns ``(gx, gy, r, disc)``: pixel offsets from the center (``gx`` rightward,
    ``gy`` downward), radius ``r = hypot(gx, gy)``, and the boolean disc mask
    ``r <= width/2``. Pixel centers are sampled at ``i + 0.5``.
    """
    center = width / 2.0
    coord = np.arange(width) + 0.5
    gx, gy = np.meshgrid(coord - center, coord - center)
    r = np.hypot(gx, gy)
    return gx, gy, r, r <= width / 2.0


def wedge_masks(width: int, n_wedges: int, rotation: float = 0.0) -> list[np.ndarray]:
    """Boolean pixel mask per wedge: inside the disc and in the wedge's sector.

    Angles are measured from 3 o'clock, counter-clockwise; wedge ``k`` spans
    ``[k, k+1) * 360/n_wedges`` degrees, with all boundaries offset
    counter-clockwise by ``rotation`` degrees (e.g. ``360/(2·n_wedges)`` puts the
    divisions halfway between the default ones).
    """
    gx, gy, _, disc = _disc_grid(width)
    angle = (np.degrees(np.arctan2(-gy, gx)) - rotation) % 360.0
    step = 360.0 / n_wedges
    return [disc & (angle >= k * step) & (angle < (k + 1) * step) for k in range(n_wedges)]


def ring_boundaries(n_rings: int, radius: float, spacing: str) -> np.ndarray:
    """Radii (length ``n_rings+1``) delimiting ``n_rings`` concentric annuli of the disc.

    Boundaries run from ``0`` (fovea) to ``radius`` (disc edge). The spacing sets
    how eccentricity is sampled:

    - ``"equal_width"``: equal radial thickness, ``r_k = radius * k/n``.
    - ``"equal_area"``: equal screen area per ring, ``r_k = radius * sqrt(k/n)``
      (rings thin outward).
    - ``"log"``: equal spacing in ``log(r + r0)`` with a small foveal offset
      ``r0 = RING_FOVEA_FRAC * radius``, so rings are thin near the fovea and thick
      in the periphery (~equal cortical area; the retinotopy standard). Still
      spans ``[0, radius]`` exactly because of the offset.
    """
    k = np.arange(n_rings + 1)
    if spacing == "equal_width":
        return radius * k / n_rings
    if spacing == "equal_area":
        return radius * np.sqrt(k / n_rings)
    r0 = radius * RING_FOVEA_FRAC                          # log (default)
    return r0 * ((radius + r0) / r0) ** (k / n_rings) - r0


def ring_masks(width: int, n_rings: int, spacing: str = "log") -> list[np.ndarray]:
    """Boolean pixel mask per ring: inside the disc and within the ring's annulus.

    Ring ``0`` is the innermost (foveal) disc; ring ``n_rings-1`` reaches the disc
    edge. ``spacing`` is one of :data:`RING_SPACINGS` (see :func:`ring_boundaries`).
    """
    _, _, r, disc = _disc_grid(width)
    radius = width / 2.0
    if spacing not in RING_SPACINGS:
        spacing = "log"
    b = ring_boundaries(n_rings, radius, spacing)
    masks = []
    for k in range(n_rings):
        if k == n_rings - 1:                              # outer ring: disc handles the edge
            masks.append(disc & (r >= b[k]))
        else:
            masks.append(disc & (r >= b[k]) & (r < b[k + 1]))
    return masks


def fade_envelope(n_frames: int, fade_frames: int) -> np.ndarray:
    """Per-state contrast envelope: ramp 0->1, hold at 1, ramp 1->0 (length ``n_frames``)."""
    env = np.ones(n_frames)
    F = min(fade_frames, n_frames // 2)
    if F > 0:
        ramp = np.linspace(0.0, 1.0, F, endpoint=False)
        env[:F] = ramp
        env[n_frames - F:] = ramp[::-1]
    return env


# ============================================================ figures
def render_design_matrix(design: np.ndarray, out: str, row_h: int = 8, col_w: int = 44) -> None:
    """Render the on/off design matrix (states x regions; red=on, white=off)."""
    L, N = design.shape
    on, off, line = np.array([200, 40, 40]), np.array([245, 245, 245]), np.array([170, 170, 170])
    img = np.empty((L * row_h, N * col_w + (N + 1), 3), np.uint8)
    img[:] = line
    for k in range(N):
        x0 = 1 + k * (col_w + 1)
        cells = np.where(design[:, k:k + 1] == 1, on, off).astype(np.uint8)
        img[:, x0:x0 + col_w, :] = np.repeat(cells, row_h, axis=0)[:, None, :]
    Image.fromarray(img, "RGB").save(out)


def render_orientation_matrix(orient_index: np.ndarray, design: np.ndarray, n_orient: int,
                              out: str, row_h: int = 8, col_w: int = 44) -> None:
    """Render the orientation matrix (states x regions; one hue per orientation, off-regions dimmed)."""
    L, N = orient_index.shape
    colors = np.array([cm.hsv(i / max(n_orient, 1))[:3] for i in range(n_orient)]) * 255
    line, dim = np.array([170, 170, 170]), 0.32
    img = np.empty((L * row_h, N * col_w + (N + 1), 3), np.uint8)
    img[:] = line
    for k in range(N):
        x0 = 1 + k * (col_w + 1)
        cells = colors[orient_index[:, k]].astype(float)
        cells[design[:, k] == 0] = BG + (cells[design[:, k] == 0] - BG) * dim
        img[:, x0:x0 + col_w, :] = np.repeat(cells.astype(np.uint8), row_h, axis=0)[:, None, :]
    Image.fromarray(img, "RGB").save(out)


# matplotlib dark theme to match the viewer
_STYLE = dict(fg="#e8eaed", muted="#9aa0a8", panel="#24272d", axbg="#1a1c20",
              data="#2bd1ff", ideal="#e02828")


def _dark(ax) -> None:
    """Apply the dark viewer theme to an Axes."""
    ax.set_facecolor(_STYLE["axbg"])
    for spine in ax.spines.values():
        spine.set_color(_STYLE["muted"])
    ax.tick_params(colors=_STYLE["muted"])
    ax.xaxis.label.set_color(_STYLE["fg"])
    ax.yaxis.label.set_color(_STYLE["fg"])
    ax.title.set_color(_STYLE["fg"])
    ax.grid(True, which="both", color="#ffffff12", lw=0.6)


def _power2d(noise: np.ndarray) -> np.ndarray:
    """Time-averaged 2-D power spectrum of a ``(width, width, n_frames)`` field."""
    width, _, n_frames = noise.shape
    power = np.zeros((width, width))
    for t in range(n_frames):
        power += np.abs(np.fft.fft2(noise[:, :, t])) ** 2
    return power / n_frames


def plot_temporal(noise: np.ndarray, fps: int, tf_lo: float, tf_hi: float,
                  shape: str, out: str) -> None:
    """Plot the measured temporal power spectrum vs. the ideal (1/f² or flat)."""
    n_frames = noise.shape[2]
    ft = np.fft.rfftfreq(n_frames, d=1.0 / fps)
    power = (np.abs(np.fft.rfft(noise, axis=2)) ** 2).mean(axis=(0, 1))
    band = (ft >= tf_lo * 0.999) & (ft <= tf_hi) & (ft > 0)
    fig, ax = plt.subplots(figsize=(4.2, 3.3), dpi=130)
    fig.patch.set_facecolor(_STYLE["panel"])
    ax.axvspan(tf_lo, tf_hi, color=_STYLE["ideal"], alpha=0.10)
    ax.loglog(ft[1:], np.maximum(power[1:], 1e-30), color=_STYLE["data"], lw=1.6, label="measured")
    if band.any():                                    # skip the ideal fit if the band is empty
        if shape == "1/f":
            K = np.median(power[band] * ft[band] ** 2)
            ax.loglog(ft[band], K / ft[band] ** 2, "--", color=_STYLE["ideal"], lw=1.4, label="ideal 1/f²")
        else:
            K = np.median(power[band])
            ax.loglog(ft[band], np.full(band.sum(), K), "--", color=_STYLE["ideal"], lw=1.4, label="flat")
    ax.set_xlabel("temporal frequency (cyc/s)")
    ax.set_ylabel("power")
    ax.set_title("Temporal spectrum")
    _dark(ax)
    ax.legend(facecolor=_STYLE["panel"], edgecolor=_STYLE["muted"], labelcolor=_STYLE["fg"], fontsize=8)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_spatial(noise: np.ndarray, width: int, sf_lo: float, sf_hi: float,
                 shape: str, out: str) -> None:
    """Plot the radially-averaged spatial power spectrum and the 2-D power map."""
    power = _power2d(noise)
    f = np.fft.fftfreq(width, d=1.0) * width
    FX, FY = np.meshgrid(f, f)
    fr = np.round(np.hypot(FX, FY)).astype(int).ravel()
    radial = np.bincount(fr, weights=power.ravel()) / np.maximum(np.bincount(fr), 1)
    frs = np.arange(len(radial))
    band = (frs >= sf_lo) & (frs <= sf_hi)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(7.0, 3.3), dpi=130)
    fig.patch.set_facecolor(_STYLE["panel"])
    ax1.axvspan(sf_lo, sf_hi, color=_STYLE["ideal"], alpha=0.10)
    ax1.loglog(frs[1:width // 2], np.maximum(radial[1:width // 2], 1e-30),
               color=_STYLE["data"], lw=1.6, label="measured")
    if band.any():                                    # skip the ideal fit if the band is empty
        if shape == "1/f":
            K = np.median(radial[band] * frs[band] ** 2)
            ax1.loglog(frs[band], K / frs[band] ** 2, "--", color=_STYLE["ideal"], lw=1.4, label="ideal 1/f²")
        else:
            K = np.median(radial[band])
            ax1.loglog(frs[band], np.full(band.sum(), K), "--", color=_STYLE["ideal"], lw=1.4, label="flat")
    ax1.set_xlabel("spatial frequency (cyc/width)")
    ax1.set_ylabel("power")
    ax1.set_title("Spatial spectrum (radial)")
    _dark(ax1)
    ax1.legend(facecolor=_STYLE["panel"], edgecolor=_STYLE["muted"], labelcolor=_STYLE["fg"], fontsize=8)

    lim = min(width // 2, int(sf_hi * 1.2) + 8)
    ax2.imshow(np.log10(np.fft.fftshift(power) + 1e-9),
               extent=[-width // 2, width // 2, -width // 2, width // 2], cmap="magma", origin="lower")
    ax2.set_xlim(-lim, lim)
    ax2.set_ylim(-lim, lim)
    ax2.set_xlabel("fx (cyc/width)")
    ax2.set_ylabel("fy (cyc/width)")
    ax2.set_title("2D power (log)")
    _dark(ax2)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_orientation(noise: np.ndarray, width: int, sf_lo: float, sf_hi: float,
                     angles: list[float], out: str) -> None:
    """Plot in-band power vs. image orientation, marking the target angles."""
    power = _power2d(noise)
    f = np.fft.fftfreq(width, d=1.0) * width
    FX, FY = np.meshgrid(f, f)
    fr = np.hypot(FX, FY)
    band = (fr >= sf_lo) & (fr <= sf_hi)
    orient = (np.degrees(np.arctan2(FY, FX)) + 90.0) % 180.0
    bins = np.arange(0, 181, 1.0)
    idx = np.clip(np.digitize(orient[band], bins) - 1, 0, len(bins) - 2)
    binned = np.zeros(len(bins) - 1)
    np.add.at(binned, idx, power[band])
    binned = binned / binned.max() if binned.max() > 0 else binned

    fig, ax = plt.subplots(figsize=(4.2, 3.3), dpi=130)
    fig.patch.set_facecolor(_STYLE["panel"])
    centers = bins[:-1] + 0.5
    ax.fill_between(centers, 0, binned, color=_STYLE["data"], alpha=0.5)
    ax.plot(centers, binned, color=_STYLE["data"], lw=1.6)
    for a in angles:
        ax.axvline(a % 180, color=_STYLE["ideal"], lw=1.0, ls="--")
    ax.set_xlim(0, 180)
    ax.set_xticks([0, 45, 90, 135, 180])
    ax.set_xlabel("image orientation (deg)  0=horiz, 90=vert")
    ax.set_ylabel("normalized power")
    ax.set_title(f"Orientation spectrum ({len(angles)})")
    _dark(ax)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


# ============================================================ top-level
def generate_movie(params: dict,
                   progress_cb: Callable[[int, int], None] | None = None,
                   cancel_cb: Callable[[], bool] | None = None) -> dict:
    """Render a full movie and its companion artifacts.

    Parameters
    ----------
    params : dict
        Overrides for :data:`DEFAULTS` (see INTERFACE.txt / README).
    progress_cb : callable(done:int, total:int) | None
        Called periodically with frame progress.
    cancel_cb : callable() -> bool | None
        Polled at each frame and before each heavy FFT; if it returns True,
        generation aborts with :class:`Cancelled`.

    Returns
    -------
    dict
        The metadata also written to ``movie_meta.json``.

    Side effects
    ------------
    Clears and repopulates ``frames/``; writes the spectra/matrix PNGs and the
    metadata files into the project directory.
    """
    p = {**DEFAULTS, **(params or {})}
    p["geometry"] = "rings" if p.get("geometry") == "rings" else "wedges"
    p["ring_spacing"] = str(p.get("ring_spacing", "log")).replace(" ", "_").replace("-", "_")
    if p["ring_spacing"] not in RING_SPACINGS:
        p["ring_spacing"] = "log"
    geometry = p["geometry"]
    W, fps = int(p["width"]), int(p["fps"])
    # Regions are capped to [1, L]: the design has L = 2**LFSR_N - 1 = 63 states, so more
    # than L regions cannot get distinct m-sequence shifts (they would collide into
    # duplicate, collinear regressors); fewer than 1 is degenerate. Mirrors the UI.
    N = max(1, min(int(p["n_wedges"]), (1 << LFSR_N) - 1))
    p["n_wedges"] = N
    K = max(1, int(p["n_orientations"]))
    frames_per_state = max(1, round(fps * float(p["wedge_sec"])))
    fade, pad = int(p["fade_frames"]), round(float(p["pad_sec"]) * fps)
    color, bg_mode = p["color"], p["background"]
    sf_lo, sf_hi, sf_shape = float(p["sf_lo"]), float(p["sf_hi"]), p["sf_shape"]
    tf_lo, tf_hi, tf_shape = float(p["tf_lo"]), float(p["tf_hi"]), p["tf_shape"]

    def check_cancel():
        if cancel_cb and cancel_cb():
            raise Cancelled()

    # --- design, geometry, filters ---
    design, orient_index, angles, L = build_design(N, K)
    # optional region subset: excluded regions are forced always-off (background).
    # The disc geometry and m-sequence are unchanged; only which regions are ever
    # shown is masked. wedge_mask is a length-N list of booleans (default all on).
    wedge_mask = p.get("wedge_mask")
    if not isinstance(wedge_mask, list) or len(wedge_mask) != N:
        wedge_mask = [True] * N
    wedge_mask = [bool(x) for x in wedge_mask]
    p["wedge_mask"] = wedge_mask
    for k in range(N):
        if not wedge_mask[k]:
            design[:, k] = 0
    # Background orientation for EVERY state (not just rendered ones), so the full
    # design is recorded. Deterministic from the design; only used when
    # background == "oriented" (the orthogonal-to-wedges angle, per state).
    bg_orient = [orthogonal_orientation([angles[orient_index[s, k]]
                                         for k in range(N) if design[s, k]])
                 for s in range(L)] if bg_mode == "oriented" else []
    n_states = min(DEMO_STATES if p["mode"] == "demo" else L, L)
    masks = (ring_masks(W, N, p["ring_spacing"]) if geometry == "rings"
             else wedge_masks(W, N, float(p["wedge_rotation"])))
    env = fade_envelope(frames_per_state, fade)

    # off-axis orientation wedges narrow as K grows, to avoid overlap
    halfwidth = min(WEDGE_HALFWIDTH, 90.0 / K / 2.0) if K > 2 else WEDGE_HALFWIDTH
    a_t_state = temporal_filter(frames_per_state, fps, tf_lo, tf_hi, tf_shape)
    H_orient = [spatial_filter(W, sf_lo, sf_hi, sf_shape, angles[i], halfwidth)[:, :, None] * a_t_state
                for i in range(K)]
    H_iso_state = spatial_filter(W, sf_lo, sf_hi, sf_shape)[:, :, None] * a_t_state
    H_iso_pad = None
    if pad > 0:
        H_iso_pad = spatial_filter(W, sf_lo, sf_hi, sf_shape)[:, :, None] * \
            temporal_filter(pad, fps, tf_lo, tf_hi, tf_shape)

    # --- output frame writer ---
    os.makedirs(FRAME_DIR, exist_ok=True)
    for fn in os.listdir(FRAME_DIR):
        if fn.endswith(".png"):
            os.remove(os.path.join(FRAME_DIR, fn))
    total = 2 * pad + n_states * frames_per_state
    saved = [0]

    # central fixation spot: a FIX_SIZE x FIX_SIZE square recolored at random every
    # FIX_BLOCK_SEC; the colour sequence is seeded (reproducible) and recorded.
    fixation = p["fixation"] == "on"
    fix_block = max(1, round(FIX_BLOCK_SEC * fps))
    fix_colors = None
    if fixation:
        n_blocks = (total + fix_block - 1) // fix_block
        fix_colors = np.random.default_rng(BASE_SEED + 4242).integers(
            0, 256, (n_blocks, 3), dtype=np.uint8)

    def emit(frame: np.ndarray) -> None:
        check_cancel()
        u8 = frame if frame.dtype == np.uint8 else np.clip(frame, 0, 255).astype(np.uint8)
        u8 = np.ascontiguousarray(u8)                      # own copy before overlay
        if fixation:
            c = W // 2
            u8[c - FIX_SIZE // 2:c + (FIX_SIZE + 1) // 2,
               c - FIX_SIZE // 2:c + (FIX_SIZE + 1) // 2, :] = fix_colors[saved[0] // fix_block]
        Image.fromarray(u8, "RGB").save(os.path.join(FRAME_DIR, f"frame_{saved[0]:05d}.png"))
        saved[0] += 1
        if progress_cb and saved[0] % 5 == 0:
            progress_cb(saved[0], total)

    def emit_padding(seed: int) -> None:
        if pad <= 0:
            return
        rgb, _ = make_noise(seed, H_iso_pad, W, pad, color)   # full-screen isotropic noise
        for t in range(pad):
            emit(rgb[:, :, t, :])

    # --- render: padding, states, padding ---
    emit_padding(BASE_SEED + 7000)
    movie_start = saved[0]

    # one luminance carrier per orientation, sampled from the first rendered state that
    # has any on-regions (a subset mask can leave the early states blank).
    spectra_lum: dict[int, np.ndarray] = {}
    capture_state = None
    for s in range(n_states):
        check_cancel()
        on = [k for k in range(N) if design[s, k]]
        if capture_state is None and on:
            capture_state = s
        background = None
        if bg_mode == "random":
            background, _ = make_noise(BASE_SEED + 5000 + s, H_iso_state, W, frames_per_state, color)
        elif bg_mode == "oriented":
            # single-orientation background, orthogonal to this state's wedges
            H_bg = spatial_filter(W, sf_lo, sf_hi, sf_shape, bg_orient[s], halfwidth)[:, :, None] * a_t_state
            background, _ = make_noise(BASE_SEED + 5000 + s, H_bg, W, frames_per_state, color)
        carriers = {}
        for k in on:
            check_cancel()
            orient = orient_index[s, k]
            capture = s == capture_state and orient not in spectra_lum
            rgb, lum = make_noise(BASE_SEED + s * N + k, H_orient[orient], W, frames_per_state,
                                  color, want_lum=capture)
            carriers[k] = rgb
            if capture:
                spectra_lum[orient] = lum
        for t in range(frames_per_state):
            base = (background[:, :, t, :].astype(float) if background is not None
                    else np.full((W, W, 3), BG, float))
            frame = base.copy()
            for k in on:                          # blend each on-region toward the background by env
                carrier = carriers[k][:, :, t, :].astype(float)
                frame[masks[k]] = (base + env[t] * (carrier - base))[masks[k]]
            emit(frame)
    movie_end = saved[0]
    emit_padding(BASE_SEED + 8000)

    # --- spectra (sampled from the first non-empty state's carriers) ---
    spectra_files = ["temporal_spectrum.png", "spatial_spectrum.png", "orientation_spectrum.png"]
    if spectra_lum:
        spectra = sum(spectra_lum.values())       # combine orientations -> multiple spikes
        plot_temporal(spectra, fps, tf_lo, tf_hi, tf_shape, os.path.join(HERE, spectra_files[0]))
        plot_spatial(spectra, W, sf_lo, sf_hi, sf_shape, os.path.join(HERE, spectra_files[1]))
        plot_orientation(spectra, W, sf_lo, sf_hi, [angles[i] for i in spectra_lum],
                         os.path.join(HERE, spectra_files[2]))
    else:
        for fn in spectra_files:                  # no on-regions anywhere: drop stale spectra
            stale = os.path.join(HERE, fn)
            if os.path.exists(stale):
                os.remove(stale)
    render_design_matrix(design, os.path.join(HERE, "design_matrix.png"))
    render_orientation_matrix(orient_index, design, K, os.path.join(HERE, "orientation_matrix.png"))

    # --- metadata for the viewer ---
    meta = {
        "geometry": geometry, "ring_spacing": p["ring_spacing"],
        "width": W, "n_wedges": N, "wedge_rotation": float(p["wedge_rotation"]),
        "fps": fps, "wedge_sec": float(p["wedge_sec"]),
        "frames_per_state": frames_per_state, "generated_states": n_states, "total_states": L,
        "pad_frames": pad, "movie_start": movie_start, "movie_end": movie_end,
        "total_frames": saved[0], "fade_frames": min(fade, frames_per_state // 2),
        "color": color, "background": bg_mode, "bg_orient": bg_orient,
        "fixation": p["fixation"],
        "fixation_block_frames": fix_block if fixation else 0,
        "fixation_colors": fix_colors.tolist() if fixation else [],
        "sf_band": [sf_lo, sf_hi], "sf_shape": sf_shape,
        "tf_band": [tf_lo, tf_hi], "tf_shape": tf_shape,
        "n_orientations": K, "orient_angles": angles,
        "wedge_mask": wedge_mask,
        "design": design.tolist(), "orient_design": orient_index.tolist(),
        "params": p,
    }
    with open(os.path.join(HERE, "movie_meta.json"), "w") as f:
        json.dump(meta, f)
    if progress_cb:
        progress_cb(total, total)
    return meta


if __name__ == "__main__":
    import time
    t0 = time.time()
    meta = generate_movie({"mode": "demo"},
                          progress_cb=lambda d, t: print(f"\r  {d}/{t} frames", end="", flush=True))
    print(f"\nDone in {time.time() - t0:.0f}s: {meta['total_frames']} frames "
          f"({meta['generated_states']} states, {meta['n_wedges']} wedges, "
          f"{meta['n_orientations']} orientations).")
