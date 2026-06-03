#!/usr/bin/env python3
"""Local HTTP bridge: browser → RANS fast pipeline → forces + flow field.

The browser can't run the GPU solver, so launch this ONCE in a normal (non-sandboxed)
interactive shell — it runs the pipeline in that shell context, where the solver works:

    <compute venv>/bin/python /home/qiqi/flexcompute/flexfoil/rans/rans_server.py

Vite proxies /api/rans → http://localhost:8077 (see flexfoil-ui/vite.config.ts).

POST /api/rans/sweep   body: { elements:[{name,contour:[[x,y]…]}…], alphas:[…], max_pseudo? }
                       resp: { jobId, n }   (starts an unsteady-as-steady sweep job)
GET  /api/rans/sweep/status?job=ID
                       resp: { points:[{alpha,CL,CD,L_over_D,flowField}], done, error, n }
                       Points stream in as the single unsteady run marches through each α,
                       so the UI updates live; if a late α stalls/diverges, the converged
                       points are still returned.
GET  /api/rans/health  → {ok:true}
GET  /api/rans/last-flowfield → most recent solve's flow field (LIC debugging)
"""
import csv
import json
import shutil
import sys
import threading
import time
import traceback
from http.server import BaseHTTPRequestHandler, ThreadingHTTPServer
from pathlib import Path
from urllib.parse import urlparse, parse_qs

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rans.pipeline import run               # noqa: E402
from rans import solve as _solve            # noqa: E402
from rans import flowfield as _ff           # noqa: E402
from rans.env import make_env               # noqa: E402

PORT = 8077
RUN_DIR = Path("/tmp/rans_server")
RUN_DIR.mkdir(parents=True, exist_ok=True)
_lock = threading.Lock()       # serialize GPU work (one solver at a time)
_counter = [0]
_env, _find = make_env(None)   # locate the compute install once

# Accuracy-oriented defaults (overridable per request). vs the old fast-preview
# (growth 1.4 / hwall 0.02 / farfield r22 n120 ≈ 7.6k cells) this is growth 1.2 /
# hwall 0.006 / farfield r50 n240 ≈ 48k cells — finer surface + BL gradation + a
# farther boundary, the level needed for trustworthy CL/CD (~6× slower solve).
DEFAULT_FARFIELD = {"type": "circle", "center": [0.5, -0.1], "radius": 50, "n": 240}
DEFAULT_MESH = {"span": 0.1, "nspan": 1, "yplus": 1, "growth": 1.2, "hwall": 0.006, "hmax": 1.5}

# In-flight / finished sweep jobs: id -> {points, done, error, n}.
_jobs: dict[str, dict] = {}
_jobs_lock = threading.Lock()


def build_case(payload: dict) -> dict:
    mesh = {**DEFAULT_MESH, **(payload.get("mesh") or {})}
    return {
        "elements": payload["elements"],
        "farfield": {**DEFAULT_FARFIELD, **(payload.get("farfield") or {})},
        "flow": {"reynolds": payload.get("reynolds", 1.0e7), "mach": payload.get("mach", 0.2),
                 "alpha_deg": payload.get("alpha", 0.0), "temperature": 288.15},
        "mesh": mesh,
        "solver": {"max_steps": payload.get("steps", 1000)},
        "fast": True,
    }


# --- progressive α-sweep (background job) ---

def _cl_cd(out: Path, k: int) -> tuple[float, float]:
    """Converged (CL, CD) at the end of physical step k, from total_forces_v2.csv. Only
    called once step k's slice exists, so its rows are present (last match = converged)."""
    rows = list(csv.reader(open(out / "total_forces_v2.csv")))
    hdr = [h.strip() for h in rows[0]]
    iCL, iCD = hdr.index("CL"), hdr.index("CD")
    for r in reversed(rows[1:]):
        c = [x.strip() for x in r if x.strip()]
        if len(c) > iCD and int(float(c[0])) == k:
            return float(c[iCL]), float(c[iCD])
    raise ValueError(f"no converged row for physical step {k}")


def _extract_point(out: Path, alphas: list[float], k: int) -> dict:
    cl, cd = _cl_cd(out, k)
    return {"alpha": alphas[k], "CL": cl, "CD": cd, "L_over_D": cl / cd if cd else None,
            "flowField": _ff.extract_flow_mesh(out, step=k + 1)}


def _set_job(job_id: str, **kw) -> None:
    with _jobs_lock:
        _jobs[job_id].update(kw)


def _scan(out: Path, alphas: list[float], done: set[int], job_id: str) -> None:
    """Append each α whose per-step slice has been written (= that step converged)."""
    for k in range(len(alphas)):
        if k not in done and (out / f"slice_centerSpan_time_{k + 1}_proc0.vtu").exists():
            done.add(k)
            with _jobs_lock:
                _jobs[job_id]["points"].append(_extract_point(out, alphas, k))
            print(f"[rans-server]   α={alphas[k]}° ({len(done)}/{len(alphas)})", flush=True)


def _sweep_worker(job_id: str, payload: dict, job_dir: Path) -> None:
    alphas = [float(a) for a in payload["alphas"]]
    n = len(alphas)
    out = job_dir / "out"
    t0 = time.time()
    try:
        with _lock:                              # serialize the GPU portion
            # mesh + preprocess only → the unsteady+UDD Flow360.json (no solve yet)
            run(job_dir / "case.json", out, solve=False, fast=True, alpha_sweep=alphas,
                sweep_max_pseudo=int(payload.get("max_pseudo", 3000)))
            _set_job(job_id, stage="solving")

            # solve in a thread; append each α as its per-step slice lands
            holder: dict = {}

            def _run_solver():
                try:
                    _solve.run_solver(out, _find, _env, gpu=int(payload.get("gpu", 0)))
                except Exception as e:           # noqa: BLE001 — divergence aborts the solver
                    holder["error"] = str(e)

            th = threading.Thread(target=_run_solver, daemon=True)
            th.start()
            done: set[int] = set()
            while th.is_alive():
                _scan(out, alphas, done, job_id)
                time.sleep(1.0)
            th.join()
            _scan(out, alphas, done, job_id)     # catch the final step

        wall = round(time.time() - t0, 1)
        err = (f"solver stopped after {len(done)}/{n} α (α≥{alphas[len(done)]}° stalled/diverged)"
               if "error" in holder and len(done) < n else None)
        _set_job(job_id, done=True, wall_s=wall, error=err)
        print(f"[rans-server]   sweep {len(done)}/{n} ({wall}s)", flush=True)
    except Exception as e:                       # noqa: BLE001
        traceback.print_exc()
        _set_job(job_id, done=True, error=str(e))


def start_sweep(payload: dict) -> dict:
    _counter[0] += 1
    job_id = f"sweep{_counter[0]}"
    job_dir = RUN_DIR / job_id
    # The counter resets on restart, so a reused dir can hold a prior run's per-step
    # slices — _scan would report them as "done" instantly. Start clean.
    shutil.rmtree(job_dir, ignore_errors=True)
    job_dir.mkdir(parents=True, exist_ok=True)
    (job_dir / "case.json").write_text(json.dumps(build_case(payload)))
    n = len(payload["alphas"])
    with _jobs_lock:
        _jobs[job_id] = {"points": [], "done": False, "error": None, "n": n, "stage": "meshing"}
    threading.Thread(target=_sweep_worker, args=(job_id, payload, job_dir), daemon=True).start()
    return {"jobId": job_id, "n": n}


class Handler(BaseHTTPRequestHandler):
    def _send(self, code: int, obj: dict) -> None:
        body = json.dumps(obj).encode()
        self.send_response(code)
        self.send_header("Content-Type", "application/json")
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Content-Length", str(len(body)))
        self.end_headers()
        self.wfile.write(body)

    def do_OPTIONS(self) -> None:
        self.send_response(204)
        self.send_header("Access-Control-Allow-Origin", "*")
        self.send_header("Access-Control-Allow-Headers", "Content-Type")
        self.send_header("Access-Control-Allow-Methods", "POST, GET, OPTIONS")
        self.end_headers()

    def do_GET(self) -> None:
        u = urlparse(self.path)
        if u.path == "/api/rans/health":
            self._send(200, {"ok": True})
        elif u.path == "/api/rans/sweep/status":
            job_id = (parse_qs(u.query).get("job") or [""])[0]
            with _jobs_lock:
                job = _jobs.get(job_id)
                snapshot = dict(job) if job else None
            self._send(200, snapshot) if snapshot else self._send(404, {"error": "unknown job"})
        elif u.path == "/api/rans/last-flowfield":
            ffs = sorted(RUN_DIR.glob("*/out/flow_field*.json"), key=lambda p: p.stat().st_mtime)
            if not ffs:
                self._send(404, {"error": "no saved flow field yet — run a case first"})
                return
            self._send(200, {"flowField": json.loads(ffs[-1].read_text()), "source": str(ffs[-1])})
        else:
            self._send(404, {"error": "not found"})

    def do_POST(self) -> None:
        if urlparse(self.path).path != "/api/rans/sweep":
            self._send(404, {"error": "not found"})
            return
        try:
            payload = json.loads(self.rfile.read(int(self.headers.get("Content-Length", 0))))
            print(f"[rans-server] sweep alphas={payload.get('alphas')} "
                  f"elements={len(payload.get('elements', []))}", flush=True)
            self._send(200, start_sweep(payload))
        except Exception as e:  # noqa: BLE001
            traceback.print_exc()
            self._send(500, {"error": str(e)})

    def log_message(self, *args) -> None:  # quiet default access log
        pass


if __name__ == "__main__":
    print(f"[rans-server] listening on http://localhost:{PORT}  "
          f"(POST /api/rans/sweep)", flush=True)
    ThreadingHTTPServer(("127.0.0.1", PORT), Handler).serve_forever()
