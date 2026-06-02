"""Turn a CaseConfig into the mesher's ``contours.txt`` input.

Assembles a global point list and closed curves: the farfield outer loop (CCW)
plus each element wall loop (wound CW, opposite the farfield, so the constrained
Delaunay carve treats them as holes). Emits the Spalding metric parameters the
C++ mesher uses to build the a-priori anisotropic boundary-layer metric.
"""
from __future__ import annotations

import math
from pathlib import Path

from .config import CaseConfig


def _signed_area(pts: list[list[float]]) -> float:
    a = 0.0
    n = len(pts)
    for i in range(n):
        x0, y0 = pts[i]
        x1, y1 = pts[(i + 1) % n]
        a += x0 * y1 - x1 * y0
    return 0.5 * a


def _farfield_loop(cfg: CaseConfig) -> list[list[float]]:
    ff = cfg.farfield
    if ff.type != "circle":
        raise ValueError(f"unsupported farfield type: {ff.type!r}")
    cx, cy = ff.center
    return [[cx + ff.radius * math.cos(2 * math.pi * i / ff.n),
             cy + ff.radius * math.sin(2 * math.pi * i / ff.n)] for i in range(ff.n)]


def write_contours(cfg: CaseConfig, path: str | Path) -> dict:
    """Write contours.txt. Returns a small summary (element names, counts)."""
    points: list[list[float]] = []
    curves: list[tuple[int, list[int]]] = []   # (is_wall, node indices, closed)

    def add_loop(pts: list[list[float]], is_wall: bool, ccw_wanted: bool):
        # normalize winding: farfield CCW, walls CW
        if (_signed_area(pts) > 0) != ccw_wanted:
            pts = pts[::-1]
        base = len(points)
        points.extend(pts)
        idx = list(range(base, base + len(pts))) + [base]   # close the loop
        curves.append((1 if is_wall else 0, idx))

    add_loop(_farfield_loop(cfg), is_wall=False, ccw_wanted=True)
    for el in cfg.elements:
        add_loop([list(p) for p in el.contour], is_wall=el.is_wall, ccw_wanted=False)

    m = cfg.mesh
    h0 = cfg.wall_h0()
    with open(path, "w") as f:
        f.write(f"H0 {h0:.8e}\nGROWTH {m.growth}\nHWALL {m.hwall}\nHMAX {m.hmax}\n")
        f.write(f"NPOINTS {len(points)}\n")
        for x, y in points:
            f.write(f"{x:.10f} {y:.10f}\n")
        f.write(f"NCURVES {len(curves)}\n")
        for is_wall, idx in curves:
            f.write(f"{is_wall} {len(idx)} " + " ".join(map(str, idx)) + "\n")

    return {
        "h0": h0,
        "n_points": len(points),
        "n_curves": len(curves),
        "elements": [e.name for e in cfg.elements],
        "wall_names": [e.name for e in cfg.elements if e.is_wall],
    }
