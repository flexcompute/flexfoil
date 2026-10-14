"""Tests for the local SQLite database."""

import sqlite3

import pytest
from flexfoil.database import RunDatabase


@pytest.fixture
def db(tmp_path):
    return RunDatabase(tmp_path / "test.db")


def _insert_sample(db, *, airfoil_name="NACA 2412", airfoil_hash="abc123", alpha=5.0, **kw):
    defaults = dict(
        reynolds=1e6, mach=0.0, ncrit=9.0, n_panels=160, max_iter=100,
        cl=0.92, cd=0.008, cm=-0.05, converged=True, iterations=15,
        residual=1e-6, x_tr_upper=0.12, x_tr_lower=0.45,
        solver_mode="viscous", success=True, cache_key="current-engine-settings",
    )
    defaults.update(kw)
    return db.insert_run(airfoil_name=airfoil_name, airfoil_hash=airfoil_hash, alpha=alpha, **defaults)


# ---------------------------------------------------------------------------
# Basic CRUD
# ---------------------------------------------------------------------------

class TestInsertAndQuery:
    def test_insert_returns_positive_id(self, db):
        rid = _insert_sample(db)
        assert rid > 0

    def test_query_returns_inserted_row(self, db):
        _insert_sample(db)
        rows = db.query_all_runs()
        assert len(rows) == 1
        assert rows[0]["airfoil_name"] == "NACA 2412"
        assert rows[0]["cl"] == 0.92

    def test_get_run_by_id(self, db):
        rid = _insert_sample(db)
        row = db.get_run(rid)
        assert row is not None
        assert row["id"] == rid

    def test_get_run_missing_id(self, db):
        assert db.get_run(9999) is None

    def test_multiple_inserts(self, db):
        _insert_sample(db, alpha=0.0, airfoil_hash="h0")
        _insert_sample(db, alpha=5.0, airfoil_hash="h5")
        assert db.row_count() == 2

    def test_preserves_repeat_evaluations(self, db):
        _insert_sample(db)
        _insert_sample(db)
        assert db.row_count() == 2

    def test_insert_different_alpha_not_duplicate(self, db):
        _insert_sample(db, alpha=0.0)
        _insert_sample(db, alpha=5.0)
        assert db.row_count() == 2


# ---------------------------------------------------------------------------
# Cache lookup
# ---------------------------------------------------------------------------

class TestCacheLookup:
    def test_cache_hit(self, db):
        _insert_sample(db, airfoil_hash="hash0012", alpha=0.0)
        hit = db.lookup_cache("hash0012", 0.0, 1e6, 0.0, 9.0, 160, 100, "current-engine-settings")
        assert hit is not None
        assert hit["cl"] == 0.92

    def test_cache_miss_different_alpha(self, db):
        _insert_sample(db, airfoil_hash="hash0012", alpha=0.0)
        miss = db.lookup_cache("hash0012", 5.0, 1e6, 0.0, 9.0, 160, 100, "current-engine-settings")
        assert miss is None

    def test_cache_miss_different_re(self, db):
        _insert_sample(db, airfoil_hash="hash0012", alpha=0.0)
        miss = db.lookup_cache("hash0012", 0.0, 2e6, 0.0, 9.0, 160, 100, "current-engine-settings")
        assert miss is None


# ---------------------------------------------------------------------------
# Filtered query
# ---------------------------------------------------------------------------

class TestQueryRuns:
    def test_filter_by_airfoil_name(self, db):
        _insert_sample(db, airfoil_name="A", airfoil_hash="hA")
        _insert_sample(db, airfoil_name="B", airfoil_hash="hB")
        rows = db.query_runs(airfoil_name="A")
        assert len(rows) == 1
        assert rows[0]["airfoil_name"] == "A"

    def test_limit_and_offset(self, db):
        for i in range(5):
            _insert_sample(db, alpha=float(i), airfoil_hash=f"h{i}")
        rows = db.query_runs(limit=2, offset=1)
        assert len(rows) == 2


# ---------------------------------------------------------------------------
# Delete
# ---------------------------------------------------------------------------

class TestDelete:
    def test_delete_all(self, db):
        _insert_sample(db)
        count = db.delete_all_runs()
        assert count == 1
        assert db.row_count() == 0

    def test_delete_empty_db(self, db):
        count = db.delete_all_runs()
        assert count == 0


# ---------------------------------------------------------------------------
# Export / import
# ---------------------------------------------------------------------------

class TestExportImport:
    def test_round_trip(self, db, tmp_path):
        _insert_sample(db)
        data = db.export_bytes()
        assert len(data) > 0

        db2 = RunDatabase(tmp_path / "imported.db")
        db2.import_bytes(data)
        rows = db2.query_all_runs()
        assert len(rows) == 1
        assert rows[0]["airfoil_name"] == "NACA 2412"

    def test_export_empty_db(self, db):
        data = db.export_bytes()
        assert len(data) > 0


# ---------------------------------------------------------------------------
# Airfoils table
# ---------------------------------------------------------------------------

class TestAirfoils:
    def test_save_and_list(self, db):
        db.save_airfoil("NACA 2412", '[{"x":1,"y":0},{"x":0,"y":0}]', 160)
        foils = db.list_airfoils()
        assert len(foils) == 1
        assert foils[0]["name"] == "NACA 2412"

    def test_get_airfoil(self, db):
        db.save_airfoil("NACA 2412", '[{"x":1,"y":0}]', 160)
        got = db.get_airfoil("NACA 2412")
        assert got is not None
        assert got["n_panels"] == 160

    def test_get_airfoil_missing(self, db):
        assert db.get_airfoil("nonexistent") is None

    def test_save_replaces_existing(self, db):
        db.save_airfoil("test", '[]', 80)
        db.save_airfoil("test", '[]', 160)
        foils = db.list_airfoils()
        assert len(foils) == 1
        assert foils[0]["n_panels"] == 160


# ---------------------------------------------------------------------------
# Schema migration
# ---------------------------------------------------------------------------

class TestSchemaMigration:
    def test_opens_old_db_without_flaps_json(self, tmp_path):
        """Simulate an older DB that was created without flaps_json."""
        old_path = tmp_path / "old.db"
        conn = sqlite3.connect(str(old_path))
        conn.executescript("""
            CREATE TABLE runs (
              id INTEGER PRIMARY KEY AUTOINCREMENT,
              airfoil_name TEXT NOT NULL,
              airfoil_hash TEXT NOT NULL,
              alpha REAL NOT NULL,
              reynolds REAL NOT NULL,
              mach REAL NOT NULL,
              ncrit REAL NOT NULL,
              n_panels INTEGER NOT NULL,
              max_iter INTEGER NOT NULL,
              cl REAL, cd REAL, cm REAL,
              converged INTEGER NOT NULL DEFAULT 0,
              iterations INTEGER, residual REAL,
              x_tr_upper REAL, x_tr_lower REAL,
              success INTEGER NOT NULL DEFAULT 0,
              error TEXT,
              coordinates_json TEXT,
              panels_json TEXT,
              created_at TEXT DEFAULT (datetime('now')),
              session_id TEXT
            );
            CREATE UNIQUE INDEX idx_cache_key
              ON runs(airfoil_hash, alpha, reynolds, mach, ncrit, n_panels, max_iter);
        """)
        conn.execute(
            "INSERT INTO runs (airfoil_name, airfoil_hash, alpha, reynolds, mach, ncrit, n_panels, max_iter, cl, success) "
            "VALUES ('test', 'h1', 0.0, 1e6, 0.0, 9.0, 160, 100, 0.5, 1)"
        )
        conn.commit()
        conn.close()

        db = RunDatabase(old_path)
        rows = db.query_all_runs()
        assert len(rows) == 1
        assert rows[0]["cl"] == 0.5

        _insert_sample(db, airfoil_hash="new")
        assert db.row_count() == 2

    def test_opens_old_db_without_solver_mode(self, tmp_path):
        """Simulate an even older DB missing solver_mode and flaps_json."""
        old_path = tmp_path / "ancient.db"
        conn = sqlite3.connect(str(old_path))
        conn.executescript("""
            CREATE TABLE runs (
              id INTEGER PRIMARY KEY AUTOINCREMENT,
              airfoil_name TEXT NOT NULL,
              airfoil_hash TEXT NOT NULL,
              alpha REAL NOT NULL,
              reynolds REAL NOT NULL,
              mach REAL NOT NULL,
              ncrit REAL NOT NULL,
              n_panels INTEGER NOT NULL,
              max_iter INTEGER NOT NULL,
              cl REAL, cd REAL, cm REAL,
              converged INTEGER NOT NULL DEFAULT 0,
              iterations INTEGER, residual REAL,
              x_tr_upper REAL, x_tr_lower REAL,
              success INTEGER NOT NULL DEFAULT 0,
              error TEXT,
              created_at TEXT DEFAULT (datetime('now')),
              session_id TEXT
            );
            CREATE UNIQUE INDEX idx_cache_key
              ON runs(airfoil_hash, alpha, reynolds, mach, ncrit, n_panels, max_iter);
        """)
        conn.commit()
        conn.close()

        db = RunDatabase(old_path)
        _insert_sample(db, airfoil_hash="migrated")
        rows = db.query_all_runs()
        assert len(rows) == 1
        assert rows[0]["solver_mode"] == "viscous"


# ---------------------------------------------------------------------------
# Pruning
# ---------------------------------------------------------------------------

class TestPruning:
    def test_prune_does_not_delete_under_limit(self, db):
        for i in range(5):
            _insert_sample(db, alpha=float(i), airfoil_hash=f"h{i}")
        assert db.row_count() == 5

class TestCacheHistory:
    def test_failed_attempt_and_legacy_run_cannot_hide_corrected_result(self, db):
        legacy = _insert_sample(db, cache_key=None, converged=False, cl=0.15)
        failed = _insert_sample(db, converged=False, cl=0.2)
        assert db.lookup_cache('abc123', 5, 1e6, 0, 9, 160, 100, 'current-engine-settings') is None
        fresh = _insert_sample(db, cl=0.4432)
        assert db.lookup_cache('abc123', 5, 1e6, 0, 9, 160, 100, 'current-engine-settings')['id'] == fresh
        assert db.get_run(legacy)['cl'] == 0.15
        assert db.get_run(failed)['cl'] == 0.2
        assert db.row_count() == 3
        assert db.lookup_cache('abc123', 5, 1e6, 0, 9, 160, 100) is None
        assert db.lookup_cache('abc123', 5, 1e6, 0, 9, 160, 100, 'other-engine-or-settings') is None

    def test_import_migrates_legacy_unique_index_without_losing_history(self, db, tmp_path):
        old_path = tmp_path / 'legacy.db'
        old = RunDatabase(old_path)
        legacy = _insert_sample(old, cache_key=None)
        old._conn.execute('CREATE UNIQUE INDEX idx_cache_key ON runs(airfoil_hash,alpha,reynolds,mach,ncrit,n_panels,max_iter)')
        old._conn.execute('DROP INDEX idx_cache_lookup')
        old._conn.execute('ALTER TABLE runs DROP COLUMN cache_key')
        data = old.export_bytes()
        old.close()
        db.import_bytes(data)
        assert db.get_run(legacy)['cl'] == 0.92
        assert db.get_run(legacy)['cache_key'] is None
        fresh = _insert_sample(db, cl=0.4432)
        assert db.row_count() == 2
        assert db.lookup_cache('abc123', 5, 1e6, 0, 9, 160, 100, 'current-engine-settings')['id'] == fresh
        restored = RunDatabase(tmp_path / 'roundtrip.db')
        restored.import_bytes(db.export_bytes())
        assert restored.row_count() == 2
        assert restored.lookup_cache('abc123', 5, 1e6, 0, 9, 160, 100, 'current-engine-settings')['cl'] == 0.4432
        restored.close()


def test_invalid_import_keeps_existing_history(db):
    saved = _insert_sample(db)
    with pytest.raises(sqlite3.DatabaseError):
        db.import_bytes(b"not a SQLite database")
    assert db.get_run(saved)["cl"] == 0.92
    assert db.row_count() == 1
