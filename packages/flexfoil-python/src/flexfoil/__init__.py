"""
flexfoil — Airfoil analysis in Python with a local web UI.

    import flexfoil

    foil = flexfoil.naca("2412")
    result = foil.solve(alpha=5.0, Re=1e6)
    print(result.cl, result.cd)

    foil.polar(alpha=(-5, 15, 0.5), Re=1e6).plot()

    flexfoil.serve()   # opens the web UI at localhost:8420
"""

from flexfoil.airfoil import Airfoil, BLResult, SolveResult
from flexfoil.database import RunDatabase, get_database
from flexfoil.polar import PolarResult

try:
    from flexfoil.rans import RANSPolarResult, RANSResult
except ImportError:
    pass  # flow360client not installed

__all__ = [
    "Airfoil",
    "BLResult",
    "SolveResult",
    "RunDatabase",
    "PolarResult",
    "naca",
    "load",
    "from_coordinates",
    "runs",
    "serve",
    "get_database",
]


def naca(designation: str, *, n_panels: int = 160) -> Airfoil:
    """Create an airfoil from a NACA 4-digit designation (e.g. '2412')."""
    return Airfoil.from_naca(designation, n_panels=n_panels)


def load(path: str, *, n_panels: int = 160) -> Airfoil:
    """Load an airfoil from a Selig/Lednicer .dat file."""
    return Airfoil.from_dat(path, n_panels=n_panels)


def from_coordinates(
    x: "list[float]", y: "list[float]", *, name: str = "custom", n_panels: int = 160
) -> Airfoil:
    """Create an airfoil from raw x, y coordinate arrays."""
    return Airfoil.from_xy(x, y, name=name, n_panels=n_panels)


def runs():
    """Return all cached solver runs as a list of dicts (or DataFrame if pandas available)."""
    db = get_database()
    rows = db.query_all_runs()
    try:
        import pandas as pd
        return pd.DataFrame(rows)
    except ImportError:
        return rows


def serve(
    *,
    port: int = 8420,
    host: str = "127.0.0.1",
    open_browser: bool = True,
    db_path: str | None = None,
):
    """Launch the local flexfoil web UI and API server."""
    from flexfoil.server import run_server
    run_server(host=host, port=port, open_browser=open_browser, db_path=db_path)
