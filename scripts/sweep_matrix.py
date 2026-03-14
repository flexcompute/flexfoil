from __future__ import annotations

from pathlib import Path


def get_sweep_cases(workspace: Path) -> list[dict]:
    """
    Standard large comparison matrix for XFOIL vs faithful Rust reproduction.

    `no_repanel=True` is only used for the existing pre-paneled NACA references
    where we already have trusted coordinate parity. All other foils are loaded
    from raw `.dat` files and repaneled on the Rust side to XFOIL-style spacing.
    """
    return [
        {
            "foil": "naca0012",
            "file": workspace / "naca0012_xfoil_paneled.dat",
            "reynolds": [5e5, 1e6, 3e6, 6e6],
            "no_repanel": True,
        },
        {
            "foil": "naca2412",
            "file": workspace / "naca2412_xfoil_paneled.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": True,
        },
        {
            "foil": "naca4412",
            "file": workspace / "naca4412_xfoil_paneled.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": True,
        },
        {
            "foil": "dae11",
            "file": workspace / "Xfoil" / "runs" / "dae11.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "dae21",
            "file": workspace / "Xfoil" / "runs" / "dae21.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "dae31",
            "file": workspace / "Xfoil" / "runs" / "dae31.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "dae51",
            "file": workspace / "Xfoil" / "runs" / "dae51.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "e387",
            "file": workspace / "Xfoil" / "runs" / "e387.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "la203",
            "file": workspace / "Xfoil" / "runs" / "la203.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "la203t",
            "file": workspace / "Xfoil" / "runs" / "la203t.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
        {
            "foil": "lnv109a",
            "file": workspace / "Xfoil" / "runs" / "lnv109a.dat",
            "reynolds": [5e5, 1e6, 3e6],
            "no_repanel": False,
        },
    ]


def case_tuples(workspace: Path) -> list[tuple[str, str]]:
    tuples: list[tuple[str, str]] = []
    for case in get_sweep_cases(workspace):
        for re in case["reynolds"]:
            re_str = f"re{re:.0e}".replace("+", "").replace(".", "")
            tuples.append((case["foil"], re_str))
    return tuples
