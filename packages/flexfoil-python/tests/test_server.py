"""Tests for the local API server."""

import pytest


@pytest.fixture(autouse=True)
def _use_tmp_db(tmp_path, monkeypatch):
    monkeypatch.setenv("FLEXFOIL_DATA_DIR", str(tmp_path))
    import flexfoil.database
    flexfoil.database._default_db = None


def _make_client(tmp_path):
    from starlette.testclient import TestClient
    from flexfoil.server import create_app
    app = create_app(db_path=str(tmp_path / "test.db"))
    return TestClient(app)


_SAMPLE_RUN = {
    "airfoil_name": "NACA 2412",
    "airfoil_hash": "testhash",
    "alpha": 5.0,
    "reynolds": 1e6,
    "mach": 0.0,
    "ncrit": 9.0,
    "n_panels": 160,
    "max_iter": 100,
            "cache_key": "current-engine-settings",
    "cl": 0.92,
    "cd": 0.008,
    "cm": -0.05,
    "converged": True,
    "iterations": 15,
    "residual": 1e-6,
    "x_tr_upper": 0.12,
    "x_tr_lower": 0.45,
    "solver_mode": "viscous",
    "success": True,
    "error": None,
    "coordinates_json": None,
    "panels_json": None,
}


# ---------------------------------------------------------------------------
# Health
# ---------------------------------------------------------------------------

class TestHealth:
    def test_health(self, tmp_path):
        client = _make_client(tmp_path)
        resp = client.get("/api/health")
        assert resp.status_code == 200
        data = resp.json()
        assert data["status"] == "ok"
        assert data["run_count"] == 0


# ---------------------------------------------------------------------------
# Runs CRUD
# ---------------------------------------------------------------------------

class TestRunsAPI:
    def test_create_run(self, tmp_path):
        client = _make_client(tmp_path)
        resp = client.post("/api/runs", json=_SAMPLE_RUN)
        assert resp.status_code == 201
        assert resp.json()["id"] > 0

    def test_list_runs(self, tmp_path):
        client = _make_client(tmp_path)
        client.post("/api/runs", json=_SAMPLE_RUN)
        resp = client.get("/api/runs")
        assert resp.status_code == 200
        rows = resp.json()
        assert len(rows) == 1
        assert rows[0]["airfoil_name"] == "NACA 2412"

    def test_get_run_by_id(self, tmp_path):
        client = _make_client(tmp_path)
        rid = client.post("/api/runs", json=_SAMPLE_RUN).json()["id"]
        resp = client.get(f"/api/runs/{rid}")
        assert resp.status_code == 200
        assert resp.json()["cl"] == 0.92

    def test_get_run_not_found(self, tmp_path):
        client = _make_client(tmp_path)
        resp = client.get("/api/runs/99999")
        assert resp.status_code == 404

    def test_delete_all_runs(self, tmp_path):
        client = _make_client(tmp_path)
        client.post("/api/runs", json=_SAMPLE_RUN)
        resp = client.delete("/api/runs")
        assert resp.status_code == 200
        assert client.get("/api/runs/count").json()["count"] == 0

    def test_run_count(self, tmp_path):
        client = _make_client(tmp_path)
        assert client.get("/api/runs/count").json()["count"] == 0
        client.post("/api/runs", json=_SAMPLE_RUN)
        assert client.get("/api/runs/count").json()["count"] == 1


# ---------------------------------------------------------------------------
# Cache lookup
# ---------------------------------------------------------------------------

class TestCacheLookup:
    def test_cache_hit(self, tmp_path):
        client = _make_client(tmp_path)
        client.post("/api/runs", json=_SAMPLE_RUN)
        resp = client.post("/api/runs/lookup", json={
            "airfoil_hash": "testhash",
            "alpha": 5.0,
            "reynolds": 1e6,
            "mach": 0.0,
            "ncrit": 9.0,
            "n_panels": 160,
            "max_iter": 100,
            "cache_key": "current-engine-settings",
        })
        assert resp.status_code == 200
        assert resp.json() is not None
        assert resp.json()["cl"] == 0.92

    def test_cache_miss(self, tmp_path):
        client = _make_client(tmp_path)
        resp = client.post("/api/runs/lookup", json={
            "airfoil_hash": "missing",
            "alpha": 0.0,
            "reynolds": 1e6,
            "mach": 0.0,
            "ncrit": 9.0,
            "n_panels": 160,
            "max_iter": 100,
            "cache_key": "current-engine-settings",
        })
        assert resp.status_code == 200
        assert resp.json() is None


# ---------------------------------------------------------------------------
# Airfoils API
# ---------------------------------------------------------------------------

class TestAirfoilsAPI:
    def test_save_and_list(self, tmp_path):
        client = _make_client(tmp_path)
        resp = client.post("/api/airfoils", json={
            "name": "NACA 2412",
            "coordinates_json": '[{"x":1,"y":0}]',
            "n_panels": 160,
        })
        assert resp.status_code == 201
        resp = client.get("/api/airfoils")
        assert len(resp.json()) == 1


# ---------------------------------------------------------------------------
# DB export/import
# ---------------------------------------------------------------------------

class TestDbExportImport:
    def test_export_returns_bytes(self, tmp_path):
        client = _make_client(tmp_path)
        client.post("/api/runs", json=_SAMPLE_RUN)
        resp = client.get("/api/db/export")
        assert resp.status_code == 200
        assert len(resp.content) > 0

    def test_import_replaces_data(self, tmp_path):
        client = _make_client(tmp_path)
        client.post("/api/runs", json=_SAMPLE_RUN)
        export_data = client.get("/api/db/export").content

        client.delete("/api/runs")
        assert client.get("/api/runs/count").json()["count"] == 0

        resp = client.post("/api/db/import", content=export_data)
        assert resp.status_code == 200
        assert client.get("/api/runs/count").json()["count"] == 1
