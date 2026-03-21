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
# RANS endpoints — Flow360 cloud CFD
# ---------------------------------------------------------------------------

_rans_jobs: dict[str, dict] = {}  # case_id → {status, result, ...}


async def rans_submit(request: Request) -> JSONResponse:
    """Submit a RANS case to Flow360 (runs in background)."""
    try:
        from flexfoil.rans.flow360 import check_auth, run_rans
    except ImportError:
        return JSONResponse(
            {"error": "flow360client not installed. Run: pip install flexfoil[rans]"},
            status_code=501,
        )

    body = await request.json()
    coords_json = body.get("coordinates_json", "[]")
    coords = [(p["x"], p["y"]) for p in json.loads(coords_json)]
    alpha = float(body.get("alpha", 0.0))
    Re = float(body.get("reynolds", 1e6))
    mach = float(body.get("mach", 0.2))
    airfoil_name = body.get("airfoil_name", "airfoil")

    if not check_auth():
        return JSONResponse(
            {"error": "Flow360 credentials not configured"},
            status_code=401,
        )

    # Generate a tracking ID
    import uuid
    job_id = str(uuid.uuid4())[:8]
    _rans_jobs[job_id] = {"status": "submitted", "result": None}

    async def _background():
        try:
            async def progress(status, frac):
                _rans_jobs[job_id]["status"] = status
                await _broadcast("rans_status", {
                    "job_id": job_id,
                    "status": status,
                    "progress": frac,
                })

            # Run in thread pool (blocking I/O)
            def _sync_progress(status, frac):
                _rans_jobs[job_id]["status"] = status

            result = await asyncio.to_thread(
                run_rans,
                coords,
                alpha=alpha,
                Re=Re,
                mach=mach,
                airfoil_name=airfoil_name,
                on_progress=_sync_progress,
            )

            _rans_jobs[job_id]["status"] = "complete" if result.success else "failed"
            _rans_jobs[job_id]["result"] = {
                "cl": result.cl,
                "cd": result.cd,
                "cm": result.cm,
                "alpha": result.alpha,
                "reynolds": result.reynolds,
                "mach": result.mach,
                "converged": result.converged,
                "success": result.success,
                "error": result.error,
                "case_id": result.case_id,
                "cd_pressure": result.cd_pressure,
                "cd_friction": result.cd_friction,
            }

            await _broadcast("rans_status", {
                "job_id": job_id,
                "status": _rans_jobs[job_id]["status"],
                "progress": 1.0,
                "result": _rans_jobs[job_id]["result"],
            })

        except Exception as e:
            _rans_jobs[job_id]["status"] = "failed"
            _rans_jobs[job_id]["result"] = {"error": str(e), "success": False}
            await _broadcast("rans_status", {
                "job_id": job_id,
                "status": "failed",
                "error": str(e),
            })

    asyncio.create_task(_background())
    return JSONResponse({"job_id": job_id, "status": "submitted"}, status_code=202)


async def rans_status(request: Request) -> JSONResponse:
    """Get status of a RANS job."""
    job_id = request.path_params["job_id"]
    job = _rans_jobs.get(job_id)
    if job is None:
        return JSONResponse({"error": "Job not found"}, status_code=404)
    return JSONResponse({"job_id": job_id, **job})


async def rans_result(request: Request) -> JSONResponse:
    """Get result of a completed RANS job."""
    job_id = request.path_params["job_id"]
    job = _rans_jobs.get(job_id)
    if job is None:
        return JSONResponse({"error": "Job not found"}, status_code=404)
    if job["result"] is None:
        return JSONResponse({"error": "Not complete yet", "status": job["status"]}, status_code=202)
    return JSONResponse(job["result"])


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


def _build_routes() -> list:
    routes = [
        Route("/api/health", health),
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
        Route("/api/rans/submit", rans_submit, methods=["POST"]),
        Route("/api/rans/status/{job_id:str}", rans_status),
        Route("/api/rans/result/{job_id:str}", rans_result),
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

    uvicorn.run(app, host=host, port=port, log_level="info")
