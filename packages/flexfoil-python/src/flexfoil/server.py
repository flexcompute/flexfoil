"""Local web server — serves the flexfoil frontend and a REST API backed by SQLite.

The frontend connects to this server's API when running in local mode,
instead of using its built-in sql.js + IndexedDB storage.
"""

from __future__ import annotations

import asyncio
import json
import logging
import os
import re
import time
import webbrowser
from functools import lru_cache
from pathlib import Path
from typing import AsyncGenerator
from urllib.request import urlopen, Request as UrlRequest
from urllib.error import URLError

from starlette.applications import Starlette
from starlette.exceptions import HTTPException
from starlette.middleware import Middleware
from starlette.middleware.cors import CORSMiddleware
from starlette.requests import Request
from starlette.responses import JSONResponse, Response, StreamingResponse
from starlette.routing import Mount, Route
from starlette.staticfiles import StaticFiles

from flexfoil.database import RunDatabase

logger = logging.getLogger("flexfoil.server")

_db: RunDatabase | None = None
_sse_subscribers: list[asyncio.Queue] = []

STATIC_DIR = Path(__file__).parent / "_static"


def _get_db() -> RunDatabase:
    assert _db is not None, "Database not initialized"
    return _db


# ---------------------------------------------------------------------------
# SSE — Server-Sent Events for live updates
# ---------------------------------------------------------------------------

async def _broadcast(event: str, data: dict) -> None:
    payload = f"event: {event}\ndata: {json.dumps(data)}\n\n"
    dead: list[asyncio.Queue] = []
    for q in _sse_subscribers:
        try:
            q.put_nowait(payload)
        except asyncio.QueueFull:
            dead.append(q)
    for q in dead:
        _sse_subscribers.remove(q)


async def sse_endpoint(request: Request) -> StreamingResponse:
    queue: asyncio.Queue[str] = asyncio.Queue(maxsize=64)
    _sse_subscribers.append(queue)

    async def event_stream() -> AsyncGenerator[str, None]:
        try:
            yield "event: connected\ndata: {}\n\n"
            while True:
                msg = await queue.get()
                yield msg
        except asyncio.CancelledError:
            pass
        finally:
            if queue in _sse_subscribers:
                _sse_subscribers.remove(queue)

    return StreamingResponse(
        event_stream(),
        media_type="text/event-stream",
        headers={"Cache-Control": "no-cache", "X-Accel-Buffering": "no"},
    )


# ---------------------------------------------------------------------------
# API routes
# ---------------------------------------------------------------------------

async def health(request: Request) -> JSONResponse:
    db = _get_db()
    return JSONResponse({
        "status": "ok",
        "db_path": str(db.path),
        "run_count": db.row_count(),
    })


async def list_runs(request: Request) -> JSONResponse:
    db = _get_db()
    name = request.query_params.get("airfoil_name")
    limit = int(request.query_params.get("limit", "1000"))
    offset = int(request.query_params.get("offset", "0"))
    rows = db.query_runs(airfoil_name=name, limit=limit, offset=offset)
    return JSONResponse(rows)


async def get_run(request: Request) -> JSONResponse:
    db = _get_db()
    run_id = int(request.path_params["id"])
    row = db.get_run(run_id)
    if row is None:
        return JSONResponse({"error": "Not found"}, status_code=404)
    return JSONResponse(row)


async def create_run(request: Request) -> JSONResponse:
    db = _get_db()
    body = await request.json()
    run_id = db.insert_run(**body)
    await _broadcast("run_added", {"id": run_id})
    return JSONResponse({"id": run_id}, status_code=201)


async def delete_runs(request: Request) -> JSONResponse:
    db = _get_db()
    count = db.delete_all_runs()
    await _broadcast("runs_cleared", {})
    return JSONResponse({"deleted": count})


async def list_airfoils(request: Request) -> JSONResponse:
    db = _get_db()
    return JSONResponse(db.list_airfoils())


async def save_airfoil(request: Request) -> JSONResponse:
    db = _get_db()
    body = await request.json()
    aid = db.save_airfoil(
        name=body["name"],
        coordinates_json=body["coordinates_json"],
        n_panels=body.get("n_panels"),
    )
    return JSONResponse({"id": aid}, status_code=201)


async def export_db(request: Request) -> Response:
    db = _get_db()
    data = db.export_bytes()
    return Response(
        content=data,
        media_type="application/octet-stream",
        headers={"Content-Disposition": "attachment; filename=flexfoil-runs.sqlite"},
    )


async def import_db(request: Request) -> JSONResponse:
    db = _get_db()
    body = await request.body()
    db.import_bytes(body)
    await _broadcast("db_imported", {})
    return JSONResponse({"status": "ok"})


async def run_count(request: Request) -> JSONResponse:
    db = _get_db()
    return JSONResponse({"count": db.row_count()})


async def lookup_cache(request: Request) -> JSONResponse:
    db = _get_db()
    body = await request.json()
    row = db.lookup_cache(
        airfoil_hash=body["airfoil_hash"],
        alpha=body["alpha"],
        reynolds=body["reynolds"],
        mach=body["mach"],
        ncrit=body["ncrit"],
        n_panels=body["n_panels"],
        max_iter=body["max_iter"],
    )
    if row is None:
        return JSONResponse(None, status_code=200)
    return JSONResponse(row)


# ---------------------------------------------------------------------------
# UIUC Selig database proxy — avoids CORS issues for the frontend
# ---------------------------------------------------------------------------

_UIUC_BASE = "https://m-selig.ae.illinois.edu/ads/coord/"
_SAFE_FILENAME = re.compile(r"^[\w\-]+\.dat$")


@lru_cache(maxsize=256)
def _fetch_uiuc_dat(filename: str) -> str:
    req = UrlRequest(
        _UIUC_BASE + filename,
        headers={"User-Agent": "flexfoil/1.0"},
    )
    with urlopen(req, timeout=15) as resp:
        return resp.read().decode("utf-8", errors="replace")


async def uiuc_proxy(request: Request) -> Response:
    filename = request.path_params["filename"]
    if not _SAFE_FILENAME.match(filename):
        return JSONResponse({"error": "Invalid filename"}, status_code=400)
    try:
        text = await asyncio.to_thread(_fetch_uiuc_dat, filename)
    except URLError as exc:
        return JSONResponse({"error": str(exc)}, status_code=502)
    except Exception as exc:
        return JSONResponse({"error": str(exc)}, status_code=500)
    return Response(text, media_type="text/plain")


# ---------------------------------------------------------------------------
# App factory
# ---------------------------------------------------------------------------

_LOCAL_META_TAG = '<meta name="flexfoil-local" content="true">'
_index_html_cache: str | None = None


def _make_index_response() -> Response:
    """Return index.html with the local-mode meta tag injected."""
    global _index_html_cache
    if _index_html_cache is None:
        index_path = STATIC_DIR / "index.html"
        raw = index_path.read_text()
        _index_html_cache = raw.replace("<head>", f"<head>\n    {_LOCAL_META_TAG}", 1)
    return Response(_index_html_cache, media_type="text/html")


class _SPAStaticFiles(StaticFiles):
    """StaticFiles with SPA fallback — unknown paths serve index.html."""

    async def get_response(self, path: str, scope):
        try:
            return await super().get_response(path, scope)
        except HTTPException as ex:
            if ex.status_code == 404:
                return _make_index_response()
            raise


# ─── CFD endpoints ──────────────────────────────────────────────────────

async def cfd_mesh(request: Request) -> JSONResponse:
    """Generate a structured O-mesh around an airfoil."""
    from flexfoil._rustfoil import cfd_generate_mesh as _cfd_mesh

    body = await request.json()
    coords = body.get("coordinates", [])
    ni = int(body.get("ni", 128))
    nj = int(body.get("nj", 64))
    far_field = float(body.get("far_field", 15.0))
    ds0 = float(body.get("ds0", 0.001))

    # Flatten [{x,y},...] to [x0,y0,x1,y1,...] if needed
    if coords and isinstance(coords[0], dict):
        flat = []
        for p in coords:
            flat.append(p["x"])
            flat.append(p["y"])
        coords = flat

    result = _cfd_mesh(coords, ni, nj, far_field, ds0)
    return JSONResponse({
        "x": result["x"],
        "y": result["y"],
        "ni": result["ni"],
        "nj": result["nj"],
    })


async def cfd_naca_mesh(request: Request) -> JSONResponse:
    """Generate a NACA 4-digit airfoil and mesh it in one call."""
    from flexfoil._rustfoil import generate_naca4, cfd_generate_mesh as _cfd_mesh

    body = await request.json()
    code = int(body.get("naca", 12))  # NACA designation as integer
    n_pts = int(body.get("n_points", 128))
    ni = int(body.get("ni", 128))
    nj = int(body.get("nj", 64))
    far_field = float(body.get("far_field", 15.0))
    ds0 = float(body.get("ds0", 0.001))

    # Generate airfoil
    foil = generate_naca4(code, n_pts)
    coords_flat = []
    for x, y in foil:
        coords_flat.append(x)
        coords_flat.append(y)

    result = _cfd_mesh(coords_flat, ni, nj, far_field, ds0)
    return JSONResponse({
        "x": result["x"],
        "y": result["y"],
        "ni": result["ni"],
        "nj": result["nj"],
        "airfoil": foil,
    })


async def cfd_page(request: Request) -> Response:
    """Serve the CFD mesh visualization page."""
    html = _CFD_HTML
    return Response(content=html, media_type="text/html")


def _build_routes() -> list:
    routes = [
        Route("/api/health", health),
        Route("/api/cfd/mesh", cfd_mesh, methods=["POST"]),
        Route("/api/cfd/naca-mesh", cfd_naca_mesh, methods=["POST"]),
        Route("/cfd", cfd_page),
        Route("/api/runs", list_runs, methods=["GET"]),
        Route("/api/runs", create_run, methods=["POST"]),
        Route("/api/runs", delete_runs, methods=["DELETE"]),
        Route("/api/runs/count", run_count),
        Route("/api/runs/lookup", lookup_cache, methods=["POST"]),
        Route("/api/runs/{id:int}", get_run),
        Route("/api/airfoils", list_airfoils, methods=["GET"]),
        Route("/api/airfoils", save_airfoil, methods=["POST"]),
        Route("/api/db/export", export_db),
        Route("/api/db/import", import_db, methods=["POST"]),
        Route("/api/uiuc-proxy/{filename:path}", uiuc_proxy),
        Route("/api/events", sse_endpoint),
    ]
    if STATIC_DIR.is_dir():
        routes.append(Mount("/", app=_SPAStaticFiles(directory=str(STATIC_DIR), html=False)))
    return routes


def create_app(db_path: str | None = None) -> Starlette:
    global _db
    _db = RunDatabase(db_path)

    return Starlette(
        routes=_build_routes(),
        middleware=[
            Middleware(
                CORSMiddleware,
                allow_origins=["*"],
                allow_methods=["*"],
                allow_headers=["*"],
            ),
        ],
    )


def _find_free_port(host: str, preferred: int, max_attempts: int = 20) -> int:
    """Return *preferred* if available, otherwise probe upward until a free port is found."""
    import socket

    for offset in range(max_attempts):
        candidate = preferred + offset
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            try:
                s.bind((host, candidate))
                return candidate
            except OSError:
                continue
    # Last resort: let the OS pick
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind((host, 0))
        return s.getsockname()[1]


def run_server(
    *,
    host: str = "127.0.0.1",
    port: int = 8420,
    open_browser: bool = True,
    db_path: str | None = None,
) -> None:
    """Start the local flexfoil server (blocking)."""
    import uvicorn

    port = _find_free_port(host, port)
    app = create_app(db_path=db_path)
    url = f"http://{host}:{port}"

    has_static = STATIC_DIR.is_dir()
    if has_static:
        logger.info("Serving frontend from %s", STATIC_DIR)
    else:
        logger.warning(
            "Frontend assets not found at %s — API-only mode. "
            "Run `npm run build` in flexfoil-ui/ and copy dist/ to %s",
            STATIC_DIR,
            STATIC_DIR,
        )

    if open_browser and has_static:
        import threading
        def _open():
            time.sleep(1.0)
            webbrowser.open(url)
        threading.Thread(target=_open, daemon=True).start()

    print(f"flexfoil server running at {url}")
    print(f"  Database: {_db.path if _db else 'N/A'}")
    if has_static:
        print(f"  Web UI:   {url}")
    print(f"  API:      {url}/api/health")
    print()

    print(f"  CFD Mesh: {url}/cfd")
    print()

    uvicorn.run(app, host=host, port=port, log_level="info")


# ─── Inline CFD mesh visualization frontend ─────────────────────────────

_CFD_HTML = r"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>FlexFoil CFD Mesh</title>
<style>
  * { margin: 0; padding: 0; box-sizing: border-box; }
  body { font-family: system-ui, -apple-system, sans-serif; background: #1a1a2e; color: #e0e0e0; }
  .header { background: #16213e; padding: 12px 24px; display: flex; align-items: center; gap: 16px; border-bottom: 1px solid #0f3460; }
  .header h1 { font-size: 18px; font-weight: 600; color: #00d4aa; }
  .header span { font-size: 13px; color: #8888aa; }
  .container { display: flex; height: calc(100vh - 48px); }
  .sidebar { width: 280px; padding: 16px; background: #16213e; border-right: 1px solid #0f3460; overflow-y: auto; flex-shrink: 0; }
  .canvas-wrap { flex: 1; position: relative; }
  canvas { width: 100%; height: 100%; display: block; }
  fieldset { border: 1px solid #0f3460; border-radius: 6px; padding: 10px; margin-bottom: 12px; }
  legend { font-size: 11px; color: #00d4aa; text-transform: uppercase; letter-spacing: 0.5px; padding: 0 6px; }
  label { display: flex; justify-content: space-between; align-items: center; font-size: 13px; margin: 4px 0; }
  input, select { background: #0a0a1a; border: 1px solid #333; color: #e0e0e0; padding: 4px 8px; border-radius: 4px; width: 90px; font-size: 13px; }
  select { width: 100%; }
  .btn { display: block; width: 100%; padding: 10px; border: none; border-radius: 6px; font-size: 14px; font-weight: 600; cursor: pointer; margin-top: 8px; }
  .btn-primary { background: #00d4aa; color: #1a1a2e; }
  .btn-primary:hover { background: #00e6bb; }
  .btn-primary:disabled { background: #335; color: #666; cursor: not-allowed; }
  .status { font-size: 12px; color: #8888aa; padding: 8px 0; text-align: center; }
  .status.ok { color: #00d4aa; }
  .status.err { color: #ff4466; }
  .info { font-size: 12px; color: #8888aa; margin-top: 8px; font-family: monospace; line-height: 1.6; }
  .info span { color: #00d4aa; }
  .zoom-hint { position: absolute; bottom: 8px; right: 12px; font-size: 11px; color: #555; }
</style>
</head>
<body>
<div class="header">
  <h1>FlexFoil CFD Mesh</h1>
  <span>Structured O-grid mesh generation</span>
</div>
<div class="container">
  <div class="sidebar">
    <fieldset>
      <legend>Airfoil</legend>
      <label>NACA <input id="naca" type="text" value="0012" style="width:70px"></label>
      <label>Points <input id="npts" type="number" value="128" min="32" max="512" step="16" style="width:70px"></label>
    </fieldset>
    <fieldset>
      <legend>Grid</legend>
      <label>ni (circum.) <input id="ni" type="number" value="128" min="32" max="512" step="16"></label>
      <label>nj (radial) <input id="nj" type="number" value="64" min="16" max="256" step="8"></label>
      <label>Far-field <input id="ff" type="number" value="15" min="5" max="50" step="1"></label>
      <label>ds0 <input id="ds0" type="number" value="0.001" min="0.0001" max="0.1" step="0.0001"></label>
    </fieldset>
    <fieldset>
      <legend>Display</legend>
      <label>Show grid lines <input id="showGrid" type="checkbox" checked style="width:auto"></label>
      <label>Show airfoil <input id="showAirfoil" type="checkbox" checked style="width:auto"></label>
      <label>Zoom level
        <select id="zoomLevel">
          <option value="full">Full domain</option>
          <option value="near" selected>Near field</option>
          <option value="wall">Wall region</option>
        </select>
      </label>
    </fieldset>
    <button class="btn btn-primary" id="genBtn" onclick="generate()">Generate Mesh</button>
    <div id="status" class="status">Ready</div>
    <div id="info" class="info"></div>
  </div>
  <div class="canvas-wrap">
    <canvas id="canvas"></canvas>
    <div class="zoom-hint">Scroll to zoom, drag to pan</div>
  </div>
</div>
<script>
const $ = id => document.getElementById(id);
let meshData = null;
let viewX = 0.5, viewY = 0, viewScale = 1;
let dragging = false, dragX = 0, dragY = 0;

async function generate() {
  const btn = $('genBtn');
  const status = $('status');
  btn.disabled = true;
  status.className = 'status';
  status.textContent = 'Generating mesh...';

  try {
    const t0 = performance.now();
    const resp = await fetch('/api/cfd/naca-mesh', {
      method: 'POST',
      headers: {'Content-Type': 'application/json'},
      body: JSON.stringify({
        naca: parseInt($('naca').value),
        n_points: +$('npts').value,
        ni: +$('ni').value,
        nj: +$('nj').value,
        far_field: +$('ff').value,
        ds0: +$('ds0').value,
      }),
    });
    const data = await resp.json();
    const dt = ((performance.now() - t0) / 1000).toFixed(2);

    meshData = data;

    // Auto-set view
    const zoom = $('zoomLevel').value;
    autoView(zoom);

    status.className = 'status ok';
    status.textContent = `Mesh generated in ${dt}s`;

    const ni = data.ni, nj = data.nj;
    $('info').innerHTML =
      `Grid: <span>${ni} x ${nj}</span> = <span>${(ni*nj).toLocaleString()}</span> cells<br>` +
      `Points: <span>${data.x.length.toLocaleString()}</span><br>` +
      `Memory: <span>${((data.x.length * 8) / 1024).toFixed(0)} KB</span>`;

    draw();
  } catch (e) {
    status.className = 'status err';
    status.textContent = 'Error: ' + e.message;
  }
  btn.disabled = false;
}

function autoView(zoom) {
  if (!meshData) return;
  const {x, y, ni, nj} = meshData;
  if (zoom === 'wall') {
    viewX = 0.5; viewY = 0; viewScale = 8;
  } else if (zoom === 'near') {
    viewX = 0.5; viewY = 0; viewScale = 2.5;
  } else {
    // Fit full domain
    let xmin=Infinity, xmax=-Infinity, ymin=Infinity, ymax=-Infinity;
    for (let k = 0; k < x.length; k++) {
      if (x[k]<xmin) xmin=x[k]; if (x[k]>xmax) xmax=x[k];
      if (y[k]<ymin) ymin=y[k]; if (y[k]>ymax) ymax=y[k];
    }
    viewX = (xmin+xmax)/2; viewY = (ymin+ymax)/2;
    const canvas = $('canvas');
    const range = Math.max(xmax-xmin, (ymax-ymin)*canvas.width/canvas.height);
    viewScale = canvas.width / range * 0.9;
  }
}

function draw() {
  const canvas = $('canvas');
  const ctx = canvas.getContext('2d');
  const dpr = window.devicePixelRatio || 1;
  canvas.width = canvas.clientWidth * dpr;
  canvas.height = canvas.clientHeight * dpr;
  ctx.setTransform(dpr, 0, 0, dpr, 0, 0);

  const W = canvas.clientWidth, H = canvas.clientHeight;
  ctx.fillStyle = '#111118';
  ctx.fillRect(0, 0, W, H);

  if (!meshData) return;

  const {x, y, ni, nj, airfoil} = meshData;
  const showGrid = $('showGrid').checked;
  const showAirfoil = $('showAirfoil').checked;

  // Transform: world → screen
  const cx = W / 2, cy = H / 2;
  const sx = (wx) => cx + (wx - viewX) * viewScale;
  const sy = (wy) => cy - (wy - viewY) * viewScale;

  if (showGrid) {
    // j-lines (circumferential) — teal
    ctx.strokeStyle = 'rgba(0, 210, 170, 1.0)';
    ctx.lineWidth = 0.7;
    for (let j = 0; j < nj; j++) {
      ctx.beginPath();
      for (let i = 0; i <= ni; i++) {
        const ii = i % ni;
        const idx = j * ni + ii;
        const px = sx(x[idx]), py = sy(y[idx]);
        if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }

    // i-lines (radial) — blue
    ctx.strokeStyle = 'rgba(60, 160, 220, 1.0)';
    ctx.lineWidth = 0.7;
    for (let i = 0; i < ni; i++) {
      ctx.beginPath();
      for (let j = 0; j < nj; j++) {
        const idx = j * ni + i;
        const px = sx(x[idx]), py = sy(y[idx]);
        if (j === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.stroke();
    }
  }

  // Airfoil surface (j=0 wall boundary)
  ctx.strokeStyle = '#ffffff';
  ctx.lineWidth = 2.5;
  ctx.beginPath();
  for (let i = 0; i <= ni; i++) {
    const ii = i % ni;
    const px = sx(x[ii]), py = sy(y[ii]);
    if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
  }
  ctx.stroke();

  // Original airfoil points
  if (showAirfoil && airfoil) {
    ctx.fillStyle = '#ff8855';
    for (const [ax, ay] of airfoil) {
      const px = sx(ax), py = sy(ay);
      ctx.fillRect(px - 2, py - 2, 4, 4);
    }
  }
}

// Zoom / Pan
const canvas = $('canvas');
canvas.addEventListener('wheel', e => {
  e.preventDefault();
  const factor = e.deltaY > 0 ? 0.9 : 1.1;
  viewScale *= factor;
  draw();
});
canvas.addEventListener('mousedown', e => { dragging = true; dragX = e.clientX; dragY = e.clientY; });
canvas.addEventListener('mousemove', e => {
  if (!dragging) return;
  viewX -= (e.clientX - dragX) / viewScale;
  viewY += (e.clientY - dragY) / viewScale;
  dragX = e.clientX; dragY = e.clientY;
  draw();
});
canvas.addEventListener('mouseup', () => dragging = false);
canvas.addEventListener('mouseleave', () => dragging = false);

$('zoomLevel').addEventListener('change', () => { autoView($('zoomLevel').value); draw(); });
$('showGrid').addEventListener('change', draw);
$('showAirfoil').addEventListener('change', draw);

window.addEventListener('resize', draw);

// Auto-generate on load
generate();
</script>
</body>
</html>
"""
