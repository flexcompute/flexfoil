"""Local SQLite database for solver run history and caching.

The schema is intentionally compatible with the browser's sql.js database
so that .sqlite export/import works seamlessly between Python and the web UI.
"""

from __future__ import annotations

import os
import sqlite3
from pathlib import Path
from tempfile import TemporaryDirectory

_SCHEMA = """\
CREATE TABLE IF NOT EXISTS runs (
  id            INTEGER PRIMARY KEY AUTOINCREMENT,
  airfoil_name  TEXT NOT NULL,
  airfoil_hash  TEXT NOT NULL,
  alpha         REAL NOT NULL,
  reynolds      REAL NOT NULL,
  mach          REAL NOT NULL,
  ncrit         REAL NOT NULL,
  n_panels      INTEGER NOT NULL,
  max_iter      INTEGER NOT NULL,
  cl            REAL,
  cd            REAL,
  cm            REAL,
  converged     INTEGER NOT NULL DEFAULT 0,
  iterations    INTEGER,
  residual      REAL,
  x_tr_upper    REAL,
  x_tr_lower    REAL,
  solver_mode   TEXT NOT NULL DEFAULT 'viscous',
  success       INTEGER NOT NULL DEFAULT 0,
  error         TEXT,
  coordinates_json TEXT,
  panels_json   TEXT,
  flaps_json    TEXT,
  created_at    TEXT DEFAULT (datetime('now')),
  session_id    TEXT
);

DROP INDEX IF EXISTS idx_cache_key;

CREATE TABLE IF NOT EXISTS airfoils (
  id              INTEGER PRIMARY KEY AUTOINCREMENT,
  name            TEXT NOT NULL UNIQUE,
  coordinates_json TEXT NOT NULL,
  n_panels        INTEGER,
  created_at      TEXT DEFAULT (datetime('now'))
);
"""

MAX_ROWS = 50_000

_REQUIRED_COLUMNS: list[tuple[str, str]] = [
    ("solver_mode", "TEXT NOT NULL DEFAULT 'viscous'"),
    ("coordinates_json", "TEXT"),
    ("panels_json", "TEXT"),
    ("flaps_json", "TEXT"),
    ("cache_key", "TEXT"),
]


def default_db_path() -> Path:
    base = Path(os.environ.get("FLEXFOIL_DATA_DIR", "~/.flexfoil")).expanduser()
    base.mkdir(parents=True, exist_ok=True)
    return base / "runs.db"


class RunDatabase:
    """Thin wrapper around a local SQLite database."""

    def __init__(self, path: str | Path | None = None):
        self.path = Path(path) if path else default_db_path()
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._conn = sqlite3.connect(str(self.path), check_same_thread=False)
        self._conn.execute("PRAGMA journal_mode=WAL")
        self._conn.execute("PRAGMA foreign_keys=ON")
        self._conn.executescript(_SCHEMA)
        self._ensure_columns()
        self._conn.commit()

    def _ensure_columns(self) -> None:
        """Add any columns that are missing from an older schema."""
        cur = self._conn.execute("PRAGMA table_info('runs')")
        existing = {row[1] for row in cur.fetchall()}
        for col_name, col_def in _REQUIRED_COLUMNS:
            if col_name not in existing:
                self._conn.execute(
                    f"ALTER TABLE runs ADD COLUMN {col_name} {col_def}"
                )

        self._conn.execute("CREATE INDEX IF NOT EXISTS idx_cache_lookup ON runs(cache_key, airfoil_hash, alpha)")

    def insert_run(
        self,
        *,
        airfoil_name: str,
        airfoil_hash: str,
        alpha: float,
        reynolds: float,
        mach: float,
        ncrit: float,
        n_panels: int,
        max_iter: int,
        cl: float | None = None,
        cd: float | None = None,
        cm: float | None = None,
        converged: bool = False,
        iterations: int | None = None,
        residual: float | None = None,
        x_tr_upper: float | None = None,
        x_tr_lower: float | None = None,
        solver_mode: str = "viscous",
        success: bool = False,
        error: str | None = None,
        coordinates_json: str | None = None,
        panels_json: str | None = None,
        flaps_json: str | None = None,
        session_id: str | None = None,
        cache_key: str | None = None,
    ) -> int:
        cur = self._conn.execute(
            """INSERT INTO runs
               (airfoil_name, airfoil_hash, alpha, reynolds, mach, ncrit,
                n_panels, max_iter, cl, cd, cm, converged, iterations, residual,
                x_tr_upper, x_tr_lower, solver_mode, success, error,
                coordinates_json, panels_json, flaps_json, session_id, cache_key)
             VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)""",
            (
                airfoil_name,
                airfoil_hash,
                alpha,
                reynolds,
                mach,
                ncrit,
                n_panels,
                max_iter,
                cl,
                cd,
                cm,
                1 if converged else 0,
                iterations,
                residual,
                x_tr_upper,
                x_tr_lower,
                solver_mode,
                1 if success else 0,
                error,
                coordinates_json,
                panels_json,
                flaps_json,
                session_id,
                cache_key,
            ),
        )
        self._conn.commit()
        self._prune()
        return cur.lastrowid or -1

    def lookup_cache(
        self,
        airfoil_hash: str,
        alpha: float,
        reynolds: float,
        mach: float,
        ncrit: float,
        n_panels: int,
        max_iter: int,
        cache_key: str | None = None,
    ) -> dict | None:
        """Reuse only successful, converged runs with known matching solver/settings provenance."""
        if not cache_key:
            return None
        cur = self._conn.execute(
            """SELECT * FROM runs
               WHERE airfoil_hash=? AND alpha=? AND reynolds=? AND mach=?
                 AND ncrit=? AND n_panels=? AND max_iter=?
                 AND cache_key=? AND success=1 AND converged=1
               ORDER BY id DESC LIMIT 1""",
            (airfoil_hash, alpha, reynolds, mach, ncrit, n_panels, max_iter, cache_key),
        )
        row = cur.fetchone()
        if row is None:
            return None
        return self._row_to_dict(cur.description, row)

    def query_all_runs(self) -> list[dict]:
        cur = self._conn.execute("SELECT * FROM runs ORDER BY id DESC")
        return [self._row_to_dict(cur.description, row) for row in cur.fetchall()]

    def query_runs(
        self,
        *,
        airfoil_name: str | None = None,
        limit: int = 1000,
        offset: int = 0,
    ) -> list[dict]:
        query = "SELECT * FROM runs"
        params: list = []
        if airfoil_name:
            query += " WHERE airfoil_name = ?"
            params.append(airfoil_name)
        query += " ORDER BY id DESC LIMIT ? OFFSET ?"
        params.extend([limit, offset])
        cur = self._conn.execute(query, params)
        return [self._row_to_dict(cur.description, row) for row in cur.fetchall()]

    def get_run(self, run_id: int) -> dict | None:
        cur = self._conn.execute("SELECT * FROM runs WHERE id = ?", (run_id,))
        row = cur.fetchone()
        if row is None:
            return None
        return self._row_to_dict(cur.description, row)

    def delete_all_runs(self) -> int:
        cur = self._conn.execute("DELETE FROM runs")
        self._conn.commit()
        return cur.rowcount

    def row_count(self) -> int:
        cur = self._conn.execute("SELECT COUNT(*) FROM runs")
        return cur.fetchone()[0]

    # -- Airfoils table --

    def save_airfoil(self, name: str, coordinates_json: str, n_panels: int | None = None) -> int:
        cur = self._conn.execute(
            """INSERT OR REPLACE INTO airfoils (name, coordinates_json, n_panels)
               VALUES (?, ?, ?)""",
            (name, coordinates_json, n_panels),
        )
        self._conn.commit()
        return cur.lastrowid or -1

    def list_airfoils(self) -> list[dict]:
        cur = self._conn.execute("SELECT * FROM airfoils ORDER BY name")
        return [self._row_to_dict(cur.description, row) for row in cur.fetchall()]

    def get_airfoil(self, name: str) -> dict | None:
        cur = self._conn.execute("SELECT * FROM airfoils WHERE name = ?", (name,))
        row = cur.fetchone()
        if row is None:
            return None
        return self._row_to_dict(cur.description, row)

    # -- Export / import --

    def export_bytes(self) -> bytes:
        """Export the entire database as raw bytes (for .sqlite download)."""
        self._conn.commit()
        self._conn.execute("PRAGMA wal_checkpoint(FULL)")
        with open(self.path, "rb") as f:
            return f.read()

    def import_bytes(self, data: bytes) -> None:
        """Replace the database with imported bytes."""
        # Validate and migrate a separate copy before replacing the live history.
        with TemporaryDirectory(dir=self.path.parent) as directory:
            candidate_path = Path(directory) / "import.db"
            candidate_path.write_bytes(data)
            candidate = RunDatabase(candidate_path)
            try:
                candidate._conn.backup(self._conn)
                self._conn.commit()
            finally:
                candidate.close()

    def close(self) -> None:
        self._conn.close()

    # -- Internal --

    def _prune(self) -> None:
        count = self.row_count()
        if count <= MAX_ROWS:
            return
        excess = count - MAX_ROWS
        self._conn.execute(
            "DELETE FROM runs WHERE id IN (SELECT id FROM runs ORDER BY id ASC LIMIT ?)",
            (excess,),
        )
        self._conn.commit()

    @staticmethod
    def _row_to_dict(description, row) -> dict:
        return {desc[0]: val for desc, val in zip(description, row)}


_default_db: RunDatabase | None = None


def get_database(path: str | None = None) -> RunDatabase:
    """Get (or create) the module-level default database instance."""
    global _default_db
    if _default_db is None:
        _default_db = RunDatabase(path)
    return _default_db
