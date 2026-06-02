"""Input schema for a quasi-2D multi-element RANS case.

A case is fully described by a contours JSON: the deployed element polylines (the
global-frame output of the webUI's ``buildConfiguration()`` / ``toPoints()``),
the farfield, the flow conditions, and the meshing knobs. Everything downstream
(mesh, Flow360 case, solve) is derived from this.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path


@dataclass
class Element:
    """One closed element contour (global frame). ``contour`` is an ordered list
    of [x, y] points; the loop is closed automatically (no need to repeat the
    first point). ``is_wall`` marks it as a viscous wall (vs. the farfield)."""
    name: str
    contour: list[list[float]]
    is_wall: bool = True


@dataclass
class Farfield:
    """Outer boundary. Currently a circle (good enough for the PoC); the radius
    is in chords."""
    type: str = "circle"
    center: list[float] = field(default_factory=lambda: [0.5, -0.1])
    radius: float = 50.0
    n: int = 240


@dataclass
class Flow:
    reynolds: float = 1.0e7          # based on chord = 1
    mach: float = 0.2
    alpha_deg: float = 0.0
    temperature: float = 288.15      # K


@dataclass
class Mesh:
    """Meshing knobs. ``span``/``nspan`` set the quasi-2D extrusion (span is
    independent of the solution; 1 cell is adequate for Flow360 — see the span
    study). The Spalding knobs drive the a-priori anisotropic boundary-layer
    metric; ``h0`` is derived from the Reynolds number when left at 0."""
    span: float = 0.1
    nspan: int = 1
    yplus: float = 1.0
    growth: float = 1.2
    hwall: float = 0.004
    hmax: float = 3.0
    h0: float = 0.0                  # 0 -> derive from Re (flat-plate y+ estimate)


@dataclass
class Solver:
    max_steps: int = 5000


@dataclass
class CaseConfig:
    elements: list[Element]
    farfield: Farfield = field(default_factory=Farfield)
    flow: Flow = field(default_factory=Flow)
    mesh: Mesh = field(default_factory=Mesh)
    solver: Solver = field(default_factory=Solver)

    # ---- derived ----
    def wall_h0(self) -> float:
        """Wall-normal spacing for y+ = ``mesh.yplus`` from a flat-plate estimate
        (Cf = 0.026/Re^(1/7)), in chord units. Honors an explicit mesh.h0."""
        if self.mesh.h0 > 0:
            return self.mesh.h0
        Re = self.flow.reynolds
        Cf = 0.026 / Re ** (1.0 / 7.0)
        utau_over_U = math.sqrt(Cf / 2.0)
        return self.mesh.yplus / (Re * utau_over_U)

    def mu_ref(self) -> float:
        """Flow360 muRef = Mach / Re_per_chord (chord is the mesh unit)."""
        return self.flow.mach / self.flow.reynolds

    # ---- (de)serialization ----
    @staticmethod
    def load(path: str | Path) -> "CaseConfig":
        d = json.loads(Path(path).read_text())
        return CaseConfig(
            elements=[Element(**e) for e in d["elements"]],
            farfield=Farfield(**d.get("farfield", {})),
            flow=Flow(**d.get("flow", {})),
            mesh=Mesh(**d.get("mesh", {})),
            solver=Solver(**d.get("solver", {})),
        )

    def dump(self, path: str | Path) -> None:
        Path(path).write_text(json.dumps(asdict(self), indent=2))
