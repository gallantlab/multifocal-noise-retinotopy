#!/usr/bin/env python3
"""Local web server for the retinotopy movie-maker UI.

Serves ``index.html`` and the generated assets, and drives generation in a
background thread so the browser stays responsive:

==================  ======================================================
Route               Behaviour
==================  ======================================================
``GET  /``          the viewer (``index.html``)
``GET  /status``    current job ``{state, done, total, version, error}``
``POST /generate``  start a job with the posted JSON parameters
``POST /cancel``    request cancellation of the running job
==================  ======================================================

``state`` is one of ``idle | running | done | cancelled | error``. The client
polls ``/status`` and reloads ``movie_meta.json`` once ``state == "done"``
(``version`` increments per completed job, for cache-busting).

Run::

    python3 server.py        # then open http://localhost:8000
"""

import http.server
import json
import os
import shutil
import socketserver
import subprocess
import threading
import traceback
from functools import partial

import generator

PORT = 8000

# --- output saving / optional mp4 encoding ---------------------------------
# The renderer writes PNG frames to ``frames/`` (scratch for the viewer). After a
# job finishes, save_outputs() copies the metadata into a user-named output folder
# and -- if the ``encode_mp4`` toggle is on -- encodes those frames to an mp4 with
# ffmpeg. The output folder is git-ignored (see .gitignore).
DEFAULT_OUTPUT_DIR = "output"   # folder (under the app dir) for saved movies + meta
# mp4 quality presets -> (libx264 CRF, pixel format). "lossless" uses yuv444p so it
# is mathematically lossless for color too; the lossy presets use yuv420p (max
# compatibility; chroma subsampling is irrelevant for bw stimuli). Noise compresses
# poorly, so lossless is large (~tens of GB) but faithful; compressed visibly alters
# the noise texture.
MP4_QUALITY = {
    "lossless":      {"crf": 0,  "pix_fmt": "yuv444p"},
    "near-lossless": {"crf": 15, "pix_fmt": "yuv420p"},
    "compressed":    {"crf": 23, "pix_fmt": "yuv420p"},
}
DEFAULT_MP4_QUALITY = "lossless"

# Shared job state, guarded by LOCK. CANCEL signals the worker to stop.
JOB = {"state": "idle", "done": 0, "total": 0, "version": 0, "error": ""}
LOCK = threading.Lock()
CANCEL = threading.Event()


def _safe_dirname(name: str) -> str:
    """One safe folder name under the app dir (no path traversal / absolute paths)."""
    return os.path.basename(str(name or "").strip()) or DEFAULT_OUTPUT_DIR


def save_outputs(params: dict, meta: dict) -> None:
    """Save the movie's metadata into the chosen output folder and, when the
    ``encode_mp4`` toggle is on, encode the rendered frames to ``<output_dir>/movie.mp4``.

    Frames stay in ``frames/`` (the viewer reads them); the named output folder is
    the durable deliverable. mp4 encoding needs ffmpeg on PATH; if it's missing or
    fails, the job still succeeds (the frames remain) and a note is logged.
    """
    out_dir = os.path.join(generator.HERE, _safe_dirname(params.get("output_dir", DEFAULT_OUTPUT_DIR)))
    os.makedirs(out_dir, exist_ok=True)

    meta_src = os.path.join(generator.HERE, "movie_meta.json")   # save the design alongside the movie
    if os.path.exists(meta_src):
        shutil.copy2(meta_src, os.path.join(out_dir, "movie_meta.json"))

    if str(params.get("encode_mp4", "off")).lower() not in ("on", "true", "1", "yes"):
        return
    if shutil.which("ffmpeg") is None:
        print("[encode] ffmpeg not found on PATH -- skipping mp4 (frames are in frames/).")
        return

    fps = int(meta.get("fps", 30))
    quality = str(params.get("mp4_quality", DEFAULT_MP4_QUALITY)).lower()
    preset = MP4_QUALITY.get(quality, MP4_QUALITY[DEFAULT_MP4_QUALITY])
    crf = int(params.get("mp4_crf", preset["crf"]))   # explicit mp4_crf still overrides
    pix_fmt = preset["pix_fmt"]
    out_mp4 = os.path.join(out_dir, "movie.mp4")
    cmd = ["ffmpeg", "-y", "-framerate", str(fps),
           "-i", os.path.join(generator.FRAME_DIR, "frame_%05d.png"),
           "-c:v", "libx264", "-preset", "medium", "-crf", str(crf),
           "-pix_fmt", pix_fmt, out_mp4]
    print(f"[encode] {out_mp4}  (fps={fps}, quality={quality}, crf={crf}, {pix_fmt}) ...")
    try:
        subprocess.run(cmd, check=True, capture_output=True)
        print(f"[encode] wrote {out_mp4}")
    except subprocess.CalledProcessError as exc:                 # non-fatal: frames still exist
        print("[encode] ffmpeg failed:\n", exc.stderr.decode(errors="replace")[-2000:])


def run_job(params: dict) -> None:
    """Worker thread: run generation, recording progress and final state in JOB."""
    def on_progress(done: int, total: int) -> None:
        with LOCK:
            JOB["done"], JOB["total"] = done, total

    try:
        meta = generator.generate_movie(params, progress_cb=on_progress, cancel_cb=CANCEL.is_set)
        save_outputs(params, meta)                    # named output folder + optional mp4 encode
        with LOCK:
            JOB["state"] = "done"
            JOB["version"] += 1
    except generator.Cancelled:
        with LOCK:
            JOB["state"] = "cancelled"
    except Exception as exc:                          # surface failures to the UI
        traceback.print_exc()
        with LOCK:
            JOB["state"] = "error"
            JOB["error"] = f"{type(exc).__name__}: {exc}"


class Handler(http.server.SimpleHTTPRequestHandler):
    """Static file serving (no-cache) plus the /status, /generate, /cancel routes."""

    def _json(self, obj: dict, code: int = 200) -> None:
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()                            # adds the no-cache header (see below)
        self.wfile.write(body)

    def end_headers(self) -> None:
        # never cache assets, so regenerated frames/metadata are always fresh
        self.send_header("Cache-Control", "no-store, max-age=0")
        super().end_headers()

    def do_GET(self) -> None:
        if self.path.split("?")[0] == "/status":
            with LOCK:
                self._json(dict(JOB))
            return
        if self.path in ("/", ""):
            self.path = "/index.html"
        super().do_GET()

    def do_POST(self) -> None:
        route = self.path.split("?")[0]
        if route == "/cancel":
            with LOCK:
                running = JOB["state"] == "running"
            if running:
                CANCEL.set()
            self._json({"ok": running})
            return
        if route != "/generate":
            self._json({"error": "unknown endpoint"}, 404)
            return
        length = int(self.headers.get("Content-Length", 0))
        try:
            params = json.loads(self.rfile.read(length) or b"{}")
        except json.JSONDecodeError as exc:
            self._json({"error": f"bad json: {exc}"}, 400)
            return
        with LOCK:
            if JOB["state"] == "running":
                self._json({"error": "already running"}, 409)
                return
            JOB.update(state="running", done=0, total=0, error="")
            CANCEL.clear()                            # clear under the lock so a /cancel can't be lost
        threading.Thread(target=run_job, args=(params,), daemon=True).start()
        self._json({"ok": True})

    def log_message(self, fmt, *args):                # quieter logs (drop /status polls)
        if "/status" not in (args[0] if args else ""):
            super().log_message(fmt, *args)


def main() -> None:
    socketserver.TCPServer.allow_reuse_address = True
    handler = partial(Handler, directory=generator.HERE)
    with socketserver.ThreadingTCPServer(("", PORT), handler) as httpd:
        print(f"Retinotopy movie-maker running at http://localhost:{PORT}/")
        print("Open that URL, set parameters, and click Generate.  Ctrl-C to stop.")
        try:
            httpd.serve_forever()
        except KeyboardInterrupt:
            print("\nstopped.")


if __name__ == "__main__":
    main()
