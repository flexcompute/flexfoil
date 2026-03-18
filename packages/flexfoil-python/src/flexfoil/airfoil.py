"""Airfoil class — the core object users interact with."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

from flexfoil._rustfoil import (
    analyze_faithful,
    analyze_inviscid,
    deflect_flap,
    generate_naca4,
    parse_dat_file,
    repanel_xfoil,
)

if TYPE_CHECKING:
    from flexfoil.polar import PolarResult


@dataclass
class SolveResult:
    """Result of a single-point aerodynamic analysis."""

    cl: float
    cd: float
    cm: float
    converged: bool
    iterations: int
    residual: float
    x_tr_upper: float
    x_tr_lower: float
    alpha: float
    reynolds: float
    mach: float
    ncrit: float
    success: bool
    error: str | None = None

    @property
    def ld(self) -> float | None:
        """Lift-to-drag ratio."""
        if self.cd and abs(self.cd) > 1e-10:
            return self.cl / self.cd
        return None

    def __repr__(self) -> str:
        if not self.success:
            return f"SolveResult(success=False, error={self.error!r})"
        conv = "converged" if self.converged else "NOT converged"
        return (
            f"SolveResult(α={self.alpha:.2f}°, Re={self.reynolds:.0e}, "
            f"CL={self.cl:.4f}, CD={self.cd:.5f}, CM={self.cm:.4f}, {conv})"
        )


class Airfoil:
    """An airfoil with coordinates, paneling, and solver methods."""

    def __init__(
        self,
        name: str,
        raw_coords: list[tuple[float, float]],
        panel_coords: list[tuple[float, float]],
    ):
        self.name = name
        self.raw_coords = raw_coords
        self.panel_coords = panel_coords
        self._hash: str | None = None

    @classmethod
    def from_naca(cls, designation: str, *, n_panels: int = 160) -> Airfoil:
        """Create from a NACA 4-digit string (e.g. '2412', '0012')."""
        if len(designation) < 4:
            designation = designation.zfill(4)
        naca_int = int(designation)
        raw = generate_naca4(naca_int)
        flat = []
        for x, y in raw:
            flat.extend([x, y])
        paneled = repanel_xfoil(flat, n_panels)
        return cls(f"NACA {designation}", raw, paneled)

    @classmethod
    def from_dat(cls, path: str, *, n_panels: int = 160) -> Airfoil:
        """Load from a Selig/Lednicer .dat file."""
        import os
        raw = parse_dat_file(path)
        if not raw:
            raise ValueError(f"No valid coordinates found in {path}")
        flat = []
        for x, y in raw:
            flat.extend([x, y])
        paneled = repanel_xfoil(flat, n_panels)
        name = os.path.splitext(os.path.basename(path))[0]
        return cls(name, raw, paneled)

    @classmethod
    def from_xy(
        cls,
        x: list[float],
        y: list[float],
        *,
        name: str = "custom",
        n_panels: int = 160,
    ) -> Airfoil:
        """Create from raw x, y arrays."""
        if len(x) != len(y):
            raise ValueError("x and y must have the same length")
        raw = list(zip(x, y))
        flat = []
        for xi, yi in raw:
            flat.extend([xi, yi])
        paneled = repanel_xfoil(flat, n_panels)
        return cls(name, raw, paneled)

    def with_flap(
        self,
        hinge_x: float = 0.75,
        deflection: float = 10.0,
        *,
        hinge_y_frac: float = 0.5,
        n_panels: int | None = None,
    ) -> Airfoil:
        """Return a new Airfoil with a plain flap deflected.

        Parameters
        ----------
        hinge_x : hinge location as fraction of chord (0-1)
        deflection : deflection angle in degrees (positive = down)
        hinge_y_frac : vertical hinge position (0 = lower surface, 1 = upper)
        n_panels : repanel count (default: same as current)
        """
        flat = self._flat_panels()
        flapped = deflect_flap(flat, hinge_x, deflection, hinge_y_frac)
        if not flapped:
            raise ValueError("Flap deflection failed (bad geometry?)")
        n = n_panels or self.n_panels
        flat_flapped = [v for x, y in flapped for v in (x, y)]
        paneled = repanel_xfoil(flat_flapped, n)
        suffix = f"+flap({hinge_x:.0%}, {deflection:+.1f}°)"
        return Airfoil(f"{self.name} {suffix}", flapped, paneled)

    @property
    def n_panels(self) -> int:
        """Number of panels. For open-TE bodies, n_panels == n_nodes."""
        return len(self.panel_coords)

    @property
    def hash(self) -> str:
        """SHA-256 of canonical panel coordinates (for cache keys)."""
        if self._hash is None:
            canonical = "|".join(
                f"{x:.8f},{y:.8f}" for x, y in self.panel_coords
            )
            self._hash = hashlib.sha256(canonical.encode()).hexdigest()[:16]
        return self._hash

    def _flat_panels(self) -> list[float]:
        flat = []
        for x, y in self.panel_coords:
            flat.extend([x, y])
        return flat

    def solve(
        self,
        alpha: float = 0.0,
        *,
        Re: float = 1e6,
        mach: float = 0.0,
        ncrit: float = 9.0,
        max_iter: int = 100,
        viscous: bool = True,
        store: bool = True,
    ) -> SolveResult:
        """Run aerodynamic analysis at a single operating point.

        Parameters
        ----------
        alpha : angle of attack in degrees
        Re : Reynolds number (ignored for inviscid)
        mach : Mach number
        ncrit : e^N transition criterion
        max_iter : maximum Newton iterations
        viscous : if False, run inviscid panel method only
        store : if True, cache result in the local database
        """
        coords = self._flat_panels()

        if viscous:
            raw = analyze_faithful(coords, alpha, Re, mach, ncrit, max_iter)
            result = SolveResult(
                cl=raw.get("cl", 0.0),
                cd=raw.get("cd", 0.0),
                cm=raw.get("cm", 0.0),
                converged=raw.get("converged", False),
                iterations=raw.get("iterations", 0),
                residual=raw.get("residual", 0.0),
                x_tr_upper=raw.get("x_tr_upper", 1.0),
                x_tr_lower=raw.get("x_tr_lower", 1.0),
                alpha=alpha,
                reynolds=Re,
                mach=mach,
                ncrit=ncrit,
                success=raw.get("success", False),
                error=raw.get("error"),
            )
        else:
            raw = analyze_inviscid(coords, alpha)
            result = SolveResult(
                cl=raw.get("cl", 0.0),
                cd=0.0,
                cm=raw.get("cm", 0.0),
                converged=True,
                iterations=0,
                residual=0.0,
                x_tr_upper=1.0,
                x_tr_lower=1.0,
                alpha=alpha,
                reynolds=0.0,
                mach=0.0,
                ncrit=0.0,
                success=raw.get("success", False),
                error=raw.get("error"),
            )

        if store and result.success:
            self._store_run(result, viscous=viscous, max_iter=max_iter)

        return result

    def polar(
        self,
        alpha: tuple[float, float, float] = (-5, 15, 0.5),
        *,
        Re: float = 1e6,
        mach: float = 0.0,
        ncrit: float = 9.0,
        max_iter: int = 100,
        viscous: bool = True,
        store: bool = True,
        parallel: bool = True,
    ) -> PolarResult:
        """Run a polar sweep over a range of angles of attack.

        Parameters
        ----------
        alpha : (start, end, step) in degrees
        Re, mach, ncrit, max_iter : solver parameters
        viscous : if False, inviscid only
        store : if True, cache each point in the local database
        parallel : if True (default), solve all alphas in parallel via rayon
        """
        from flexfoil.polar import PolarResult

        import numpy as np

        start, end, step = alpha
        alphas = np.arange(start, end + step * 0.5, step)

        if parallel:
            results = self._polar_batch(
                [float(a) for a in alphas],
                Re=Re, mach=mach, ncrit=ncrit, max_iter=max_iter,
                viscous=viscous, store=store,
            )
        else:
            results = [
                self.solve(
                    float(a), Re=Re, mach=mach, ncrit=ncrit,
                    max_iter=max_iter, viscous=viscous, store=store,
                )
                for a in alphas
            ]

        return PolarResult(
            airfoil_name=self.name,
            reynolds=Re,
            mach=mach,
            ncrit=ncrit,
            results=results,
        )

    def _polar_batch(
        self,
        alphas: list[float],
        *,
        Re: float,
        mach: float,
        ncrit: float,
        max_iter: int,
        viscous: bool,
        store: bool,
    ) -> list[SolveResult]:
        """Solve all alphas in a single parallel Rust call."""
        coords = self._flat_panels()

        if viscous:
            from flexfoil._rustfoil import analyze_faithful_batch

            raw_list = analyze_faithful_batch(coords, alphas, Re, mach, ncrit, max_iter)
            results = [
                SolveResult(
                    cl=raw.get("cl", 0.0),
                    cd=raw.get("cd", 0.0),
                    cm=raw.get("cm", 0.0),
                    converged=raw.get("converged", False),
                    iterations=raw.get("iterations", 0),
                    residual=raw.get("residual", 0.0),
                    x_tr_upper=raw.get("x_tr_upper", 1.0),
                    x_tr_lower=raw.get("x_tr_lower", 1.0),
                    alpha=a,
                    reynolds=Re,
                    mach=mach,
                    ncrit=ncrit,
                    success=raw.get("success", False),
                    error=raw.get("error"),
                )
                for a, raw in zip(alphas, raw_list)
            ]
        else:
            from flexfoil._rustfoil import analyze_inviscid_batch

            raw_list = analyze_inviscid_batch(coords, alphas)
            results = [
                SolveResult(
                    cl=raw.get("cl", 0.0),
                    cd=0.0,
                    cm=raw.get("cm", 0.0),
                    converged=True,
                    iterations=0,
                    residual=0.0,
                    x_tr_upper=1.0,
                    x_tr_lower=1.0,
                    alpha=a,
                    reynolds=0.0,
                    mach=0.0,
                    ncrit=0.0,
                    success=raw.get("success", False),
                    error=raw.get("error"),
                )
                for a, raw in zip(alphas, raw_list)
            ]

        if store:
            for r in results:
                if r.success:
                    self._store_run(r, viscous=viscous, max_iter=max_iter)

        return results

    def _store_run(self, result: SolveResult, *, viscous: bool, max_iter: int) -> None:
        """Insert a run into the local database."""
        import json
        from flexfoil.database import get_database

        coords_json = json.dumps([{"x": x, "y": y} for x, y in self.raw_coords])
        panels_json = json.dumps([{"x": x, "y": y} for x, y in self.panel_coords])

        db = get_database()
        db.insert_run(
            airfoil_name=self.name,
            airfoil_hash=self.hash,
            alpha=result.alpha,
            reynolds=result.reynolds,
            mach=result.mach,
            ncrit=result.ncrit,
            n_panels=self.n_panels,
            max_iter=max_iter,
            cl=result.cl,
            cd=result.cd,
            cm=result.cm,
            converged=result.converged,
            iterations=result.iterations,
            residual=result.residual,
            x_tr_upper=result.x_tr_upper,
            x_tr_lower=result.x_tr_lower,
            solver_mode="viscous" if viscous else "inviscid",
            success=result.success,
            error=result.error,
            coordinates_json=coords_json,
            panels_json=panels_json,
        )

    def __repr__(self) -> str:
        return f"Airfoil({self.name!r}, n_panels={self.n_panels})"
