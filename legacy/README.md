# Legacy / development scripts

These are the **incremental development stages** of the stimulus, kept for
reference. They are **superseded by `../generator.py`**, which is parameterized
and used by the web UI. They are not needed to run the tool and are not
maintained.

Development order:

1. `make_movie.py` — rotating solid-red pie wedge (the original sanity check).
2. `make_msequence_movie.py` — multifocal m-sequence design (wedges on/off).
3. `make_noise_movie.py` — wedges filled with bandpass 1/f spatiotemporal noise
   (then vertical-only orientation, then multicolor, then per-wedge V/H, then fades).
4. `generate_design_matrix.py` — renders the on/off design matrix PNG.

Each writes PNG frames into `../frames/` and assorted PNG/JSON artifacts into the
project root. See `../README.md` for the current pipeline.
