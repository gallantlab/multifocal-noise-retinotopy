# Multifocal m-sequence retinotopy noise stimulus

A generator and web UI for building **spatiotemporal noise movies** for
retinotopic mapping and voxelwise encoding-model experiments. The stimulus is a
pie of wedges; each wedge is independently turned on/off by an m-sequence and
filled with bandpass **1/f noise** of a specific **orientation** (set by a second
m-sequence). Movies are written as a sequence of PNG frames.

![viewer](docs/screenshot.png)

---

## Quick start

Requirements: Python 3 with `numpy`, `Pillow`, `matplotlib` (no other deps).

```bash
python3 server.py          # starts a local server (stdlib only)
# open http://localhost:8000  in a browser
```

Set parameters in the left panel, click **Generate movie**, and the viewer shows
the result (movie, design matrices, spectra). **Cancel** stops a long run.

> The page must be opened through the server (`http://localhost:8000`), **not**
> as a `file://` — generation needs the backend. Opening the file directly shows
> a message telling you this.

The movie is written to `frames/frame_00000.png …`. Metadata and plots are
written alongside (see [Outputs](#outputs)).

---

## What the stimulus is

- The display is a disc (diameter = movie width) split into **N equal wedges**.
- An **on/off m-sequence** turns each wedge on or off in each *state* (multifocal
  design): each wedge follows a distinct circular shift of one binary
  maximum-length sequence, so the wedge regressors are near-orthogonal.
- Each *on* wedge is filled with **multicolored bandpass 1/f noise**, confined to
  a single **orientation** chosen by a **second, independent m-sequence**.
- The noise is **independent per wedge** (its own random seed), so neighboring
  wedges show a discontinuity at their shared boundary.
- Each state lasts `wedge_sec` seconds and **fades** in/out to the background.
- The movie is optionally **padded** before and after with full-screen isotropic
  noise.

A *state* = one configuration of (which wedges are on) × (each on-wedge's
orientation). There are **63 states** (a length-2⁶−1 m-sequence). "Demo" renders
the first 5 states; "full" renders all 63.

---

## Files

| File | Role |
|------|------|
| `generator.py` | The parameterized generator: `generate_movie(params, progress_cb, cancel_cb)`. Builds both m-sequences, renders frames, spectra, matrices, and `movie_meta.json`. |
| `server.py` | Local HTTP server (stdlib). Serves the UI and assets; `POST /generate`, `GET /status`, `POST /cancel`. |
| `index.html` | Web UI (parameter form + Generate/Cancel) and the viewer (movie player, on/off & orientation matrices, spectra). Reads `movie_meta.json`. |
| `INTERFACE.txt` | The parameter spec (name, units, default) the UI is built from. |
| `make_movie.py`, `make_msequence_movie.py`, `make_noise_movie.py`, `generate_design_matrix.py` | **Legacy / incremental** development scripts (rotating wedge → m-sequence → noise → orientation → …). Superseded by `generator.py`; kept for reference. |

---

## Parameters

From `INTERFACE.txt` (defaults shown). All are exposed in the UI.

| Parameter | Units | Default | Notes |
|-----------|-------|---------|-------|
| movie width | pixels | 512 | Frame is `width × width`. |
| # wedges | – | 8 | Equal angular wedges of the disc. |
| wedges shown | – | all on | Per-wedge include/exclude toggles (with all/none/alternate presets). Excluded wedges are always off; the geometry and m-sequence are unchanged. |
| wedge rotation | degrees | 22.5 | Offsets all wedge boundaries counter-clockwise (`360/(2·N)` = 22.5° for 8 wedges puts the divisions halfway between the cardinal axes). Geometry only. |
| wedge duration | sec | 4.1 | Length of each state → `frames_per_state = round(fps·dur)`. |
| frame rate | Hz | 30 | Render/display rate. |
| color / BW | – | color | `color` = 3 independent RGB channels; `bw` = grayscale. |
| TF shape | 1/f \| flat | 1/f | Temporal amplitude spectrum. |
| lowest / highest TF | Hz | 0.5 / 15 | Temporal passband. |
| SF shape | 1/f \| flat | 1/f | Spatial amplitude spectrum. |
| lowest / highest SF | cyc/width | 2 / 128 | Spatial passband (cycles per movie width). |
| # orientations | – | 2 | Equally spaced image orientations `i·180/K`. |
| background | gray \| random \| oriented | gray | `random` = isotropic noise behind the wedges; `oriented` = single-orientation noise whose angle is, each state, as orthogonal as possible to that state's wedge orientations. |
| fade frames | frames | 5 | Per-state fade in/out. |
| padding | sec | 2 | Full-screen isotropic noise before & after the movie. |
| fixation spot | off \| on | off | A 2×2-px square at the display center that switches to a random color every 0.5 s (baked into every frame; color sequence recorded in metadata). |
| demo / full | – | demo | `demo` = first 5 states, `full` = all 63. |

"1/f" means **amplitude ∝ 1/f** (power ∝ 1/f²).

---

## UI reference (every control)

### Parameters panel

- **movie width (pix)** — side length of the square frame (`width × width`). Sets
  the spatial resolution and the spatial Nyquist (`width/2` cyc/width). Larger
  widths mean sharper noise but slower generation and larger files. Default 512.

- **\# wedges** — number of equal angular sectors the disc is divided into. Each
  wedge becomes one on/off regressor (with m-sequence shift `round(k·63/N)`).
  Changing this rebuilds the "wedges shown" toggles (all on). Default 8.

- **wedges shown** — one toggle chip per wedge (red = shown, dark = hidden), with
  **all / none / alt** preset links (`alt` = every other wedge). Hidden wedges are
  forced always-off and never display noise; the disc geometry and m-sequence are
  unchanged. Use this to present, say, only every other wedge.

- **wedge rotation (deg)** — rotates all wedge boundaries counter-clockwise by this
  many degrees. Purely geometric — it moves where the divisions fall without
  touching the m-sequence, noise, or orientations. e.g. `22.5` (= 360/(2·8)) puts
  the divisions halfway between the cardinal axes. Default 22.5.

- **wedge duration (sec)** — how long each of the 63 states is shown. Sets
  `frames_per_state = round(fps · duration)` and the **lowest temporal frequency**
  the noise can contain within a state (`1/duration`, ≈ 0.24 Hz at 4.1 s). Default 4.1.

- **frame rate (hz)** — render/display rate. Caps the temporal band at the Nyquist
  `fps/2` (15 Hz at 30 Hz). Default 30.

- **color / BW** — `color` draws three independent noise fields (one per RGB
  channel) through the same filter, so each Fourier component gets a random colour;
  `bw` uses one grayscale field. Default color.

- **temporal frequency → shape** — `1/f` (amplitude ∝ 1/f) or `flat` (equal
  amplitude in band) temporal spectrum. Default 1/f.
- **lowest / highest TF (hz)** — temporal passband edges. Keep `highest ≤ fps/2`
  (no aliasing) and `lowest ≥ 1/wedge_sec` (representable in a state). Defaults 0.5 / 15.

- **spatial frequency → shape** — `1/f` or `flat` spatial (radial) spectrum. Default 1/f.
- **lowest / highest SF (cyc/width)** — spatial passband in cycles per movie width.
  Keep within `[~1, width/2]`. Defaults 2 / 128.

- **\# orientations** — number `K` of equally spaced image orientations (`i·180/K`)
  the wedge/background textures can take. `K=2` → {0°, 90°}; `K=4` → add 45°/135°.
  Each on-wedge's orientation is chosen by the orientation m-sequence. **Note:** 0°
  and 90° are exact infinitely-narrow orientation bands; other angles use a thin
  Fourier wedge (see caveats). Default 2.

- **background** — what fills the disc behind/around the wedges:
  - `gray` — flat mid-gray.
  - `random` — full-frame **isotropic** 1/f noise (all orientations).
  - `oriented` — full-frame **single-orientation** 1/f noise whose angle is chosen,
    *each state*, to be as orthogonal as possible to that state's on-wedge
    orientations (the midpoint of the largest empty arc on the orientation circle;
    exactly perpendicular when the wedges share one orientation, 45° when both 0°
    and 90° are present). The angle changes from state to state and is recorded for
    all 63 states in `bg_orient` (see Outputs). Default gray.

- **fade between states (frames)** — length of the per-state contrast ramp. Each
  state fades its wedges in from / out to the background over this many frames (so a
  gray background returns to gray between states; an oriented background stays).
  Default 5.

- **padding before/after (sec)** — seconds of full-screen **isotropic** 1/f noise
  prepended and appended to the movie (e.g. for fMRI run lead-in/out). Default 2.

- **fixation spot** — `on` overlays a **2×2-pixel square at the exact center** of
  every frame (padding included) that switches to a new **random color every
  0.5 s**. Drawn on top at full contrast as a fixation target; the seeded
  per-block color sequence is saved to `fixation_colors` in the metadata. Default off.

- **demo / full** — `demo` renders only the first 5 states (quick iteration);
  `full` renders all 63. Default demo.

### Buttons

- **Generate movie** — POSTs the current parameters to the server, which renders in
  the background; the progress bar tracks frames and the viewer reloads when done.
- **Cancel** — appears while generating; aborts the run within ~one wedge's compute.
  A cancelled run leaves `frames/` partial (just regenerate).

### Viewer

- **Movie** — canvas player (Pause/Play), looping the rendered frames at the true
  frame rate. The pill shows `states / preview / pad`; the readout shows the current
  frame, time, state (or "padding"), and the on-wedges with their orientations.
- **On/off m-sequence** — the design matrix (states top→bottom × wedges W1…WN; red =
  on). A cyan cursor and a time axis track playback; active wedge columns highlight.
- **Orientation m-sequence** — orientation index per state×wedge, one hue per
  orientation (legend lists the angles), off-wedges dimmed.
- **Noise spectra** — measured orientation, spatial (radial + 2-D), and temporal
  spectra of the generated noise vs. the ideal, for verification. Sampled from
  state 0, so they show the orientations present in that state.

---

## Important stimulus-generation issues

These are real constraints of the method — read before designing an experiment.

### 1. Temporal frequency is limited by the frame rate (Nyquist)
The highest representable temporal frequency is **`fps/2`**. At 30 Hz that is
**15 Hz**. Requesting a higher `highest TF` will alias and the realized spectrum
will not match the request. The UI does **not** clamp this — keep
`highest TF ≤ fps/2`.

### 2. The lowest temporal frequency is set by the state duration
Each state is an **independent** noise block of length `frames_per_state`. The
lowest temporal frequency it can contain is **`1 / wedge_sec`** (= `fps /
frames_per_state`). For a 4.1 s state that is ≈ **0.244 Hz**. A `lowest TF`
below this can't be represented within a state (the default 0.5 Hz is safe).
There is no temporal continuity *across* states.

### 3. Orientation filtering and the discrete FFT grid
Orientation is imposed by keeping only the Fourier components along one line
through the origin (an infinitely narrow orientation band):

- **0° and 90°** are axis-aligned (`fx=0` / `fy=0`), so they are kept **exactly** —
  a true infinitely-narrow line.
- **Off-axis orientations** (e.g. 45° and 135° when `# orientations = 4`) cannot
  be an exact line, because an arbitrary-angle line misses the integer FFT grid.
  They instead use a **thin angular wedge** (a few degrees of orientation
  bandwidth).

**Consequence:** with `# orientations = 4`, the 0°/90° orientations are spectrally
purer than the 45°/135° ones. With `# orientations = 2` (default), both are exact.

### 4. "# orientations" does not change the m-sequence — it changes bit decoding
There is **one** binary orientation m-sequence (length 63). For `K` orientations,
each wedge reads a window of **`ceil(log2 K)` bits** of that sequence (at its own
shift) and takes the value `mod K`:

- `K = 2` → read 1 bit → {0°, 90°}.
- `K = 4` → read 2 bits → {0°, 45°, 90°, 135°}.

This is a pragmatic decoder, **not** a true K-ary (GF(K)) maximum-length
sequence, so the ideal balance/decorrelation guarantees do **not** hold exactly
for `K > 2`. Measured balance over all 63×8 assignments:

| K | counts per orientation | ideal even |
|---|------------------------|-----------|
| 2 | 248 / 256              | 252 |
| 4 | 120 / 128 / 128 / 128  | 126 |

Powers of two (2, 4, 8) stay close to balanced; non-powers (e.g. 3) are more
biased because `mod K` folds the bit values unevenly. A proper K-ary m-sequence
and grid-aligned off-axis orientations are possible but not yet implemented.

### 5. Two decorrelated m-sequences (location vs. orientation)
On/off uses primitive polynomial `x⁶+x+1` (taps `[6,1]`); orientation uses
`x⁶+x⁴+x³+x+1` (taps `[6,4,3,1]`). Both have length 63. They are essentially
uncorrelated (per-wedge correlation ≈ −0.02), so **wedge location and orientation
are separately estimable** in an encoding model.

### 6. Spatial frequency units and limits
`cyc/width` = cycles per movie width, which equals the integer FFT index when the
field is `width` pixels. The spatial Nyquist is **`width/2`** cyc/width (256 for a
512-px movie). Keep the SF passband within `[~1, width/2]`.

### 7. Independent per-wedge noise → boundary discontinuities
Each wedge gets its own seed, so adjacent wedges (even of the same orientation)
are independent and show a visible seam at the boundary. This is intentional.

### 7a. Wedge subsets
"Wedges shown" masks the design: excluded wedges are forced always-off
(`design[:, k] = 0`) and never display noise, while the disc geometry and the
m-sequence shifts of the remaining wedges are unchanged. So "every other wedge"
keeps the 8-wedge layout but only fills wedges 1, 3, 5, 7. The excluded columns
appear empty in the on/off matrix and dimmed in the orientation matrix.

"Wedge rotation" is purely geometric: it offsets where the wedge boundaries fall
(`angle - rotation`) without touching the m-sequence, noise, or orientations.

### 8. Color vs. BW
`color` draws **three independent noise fields** through the same filter, so every
Fourier component carries an independent random RGB (randomly multicolored) while
each channel keeps the exact target spectrum. `bw` uses a single grayscale field.

### 9. Background and padding noise
`random` background and the start/end padding are **isotropic** 1/f noise (all
orientations present), i.e. the same parameters but **without** orientation
filtering. The `oriented` background is single-orientation noise whose angle is
chosen *per state* to be maximally orthogonal to that state's wedge orientations;
when a state mixes orientations the background takes the best-compromise angle
(e.g. 45° for {0°, 90°}), which may itself be off-axis (a thin Fourier wedge). The
per-state background orientations are saved in `bg_orient` for the full design.

### 10. Fades are to background, not crossfades
Each state fades its contrast in from / out to the background over `fade frames`.
States do **not** overlap (no crossfade), so between states the screen returns to
the background (gray, or the background noise).

### 11. Contrast clipping
Noise is mapped `128 + 38·z` and clipped to `[0,255]` (≈0.05% of pixels clip for
the default 1/f settings). Flat spectra or very wide bands concentrate more
energy and can clip more, slightly distorting the realized spectrum. Lower the
gain in `generator.py` (`GAIN`) if clipping matters for your stimulus.

### 12. Frame layout, scale, and cancellation
- Layout: `[pad] [state 0] … [state n] [pad]`. `movie_meta.json` records
  `pad_frames`, `frames_per_state`, `generated_states`, etc. for the frame→state
  mapping.
- A **full 512-px movie is large** (thousands of frames, hundreds of MB) and can
  take many minutes. Use **demo** while iterating.
- **Cancel** is cooperative: generation checks a flag at each frame and before
  each per-wedge FFT, so it stops within roughly one wedge's compute. A cancelled
  run leaves `frames/` partially written — just regenerate.
- The viewer preloads at most ~900 frames; longer movies are truncated **in the
  preview only** (all frames are still written to disk).

---

## Outputs

| File | Contents |
|------|----------|
| `frames/frame_NNNNN.png` | The movie, one PNG per frame (`width × width`, RGB). |
| `movie_meta.json` / `movie_meta.js` | The **complete design + parameters**, enough to reconstruct the stimulus: `design` (states × wedges on/off), `orient_design` (orientation index per state × wedge), `orient_angles`, `wedge_mask`, `wedge_rotation`, `bg_orient` (background orientation for **every** state in oriented mode), `fixation_colors` (the per-0.5 s spot colors) + `fixation_block_frames`, timing, bands, and the echoed `params`. All design arrays cover all 63 states even in demo mode. |
| `design_matrix.png` | On/off m-sequence (states × wedges). |
| `orientation_matrix.png` | Orientation m-sequence (states × wedges; one color per orientation, faded where off). |
| `temporal_spectrum.png`, `spatial_spectrum.png`, `orientation_spectrum.png` | Measured spectra of the generated noise vs. the ideal, for verification. The spectra are sampled from state 0, so they show the orientations present in that state. |

---

## Notes

- The number of states is fixed at 63 (a 6-bit m-sequence). Demo = first 5.
- Random seeds are fixed (`BASE_SEED` in `generator.py`), so a given parameter set
  reproduces the same movie.
- Earlier development happened in stages (rotating wedge → m-sequence → 1/f noise
  → vertical-only orientation → multicolor → per-wedge V/H → fades → UI); the
  `make_*.py` scripts capture those steps and are superseded by `generator.py`.
