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
import socketserver
import threading
from functools import partial

import generator

PORT = 8000

# Shared job state, guarded by LOCK. CANCEL signals the worker to stop.
JOB = {"state": "idle", "done": 0, "total": 0, "version": 0, "error": ""}
LOCK = threading.Lock()
CANCEL = threading.Event()


def run_job(params: dict) -> None:
    """Worker thread: run generation, recording progress and final state in JOB."""
    def on_progress(done: int, total: int) -> None:
        with LOCK:
            JOB["done"], JOB["total"] = done, total

    try:
        generator.generate_movie(params, progress_cb=on_progress, cancel_cb=CANCEL.is_set)
        with LOCK:
            JOB["state"] = "done"
            JOB["version"] += 1
    except generator.Cancelled:
        with LOCK:
            JOB["state"] = "cancelled"
    except Exception as exc:                          # surface failures to the UI
        import traceback
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
        self.send_header("Cache-Control", "no-store")
        self.end_headers()
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
        CANCEL.clear()
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
