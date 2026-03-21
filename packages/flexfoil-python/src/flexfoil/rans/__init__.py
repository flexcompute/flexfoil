"""RANS CFD analysis via Flow360 — pseudo-2D airfoil simulations."""

from __future__ import annotations

from dataclasses import dataclass, field


@dataclass
class RANSResult:
    """Result of a RANS CFD analysis via Flow360."""

    cl: float
    cd: float
    cm: float
    alpha: float
    reynolds: float
    mach: float
    converged: bool
    success: bool
    error: str | None = None

    # Force breakdown
    cd_pressure: float | None = None
    cd_friction: float | None = None

    # Flow360 identifiers
    case_id: str | None = None
    mesh_id: str | None = None
    wall_time_s: float | None = None

    # Surface distributions
    cp: list[float] | None = None
    cp_x: list[float] | None = None
    cf: list[float] | None = None
    cf_x: list[float] | None = None

    @property
    def ld(self) -> float | None:
        """Lift-to-drag ratio."""
        if self.cd and abs(self.cd) > 1e-10:
            return self.cl / self.cd
        return None

    def __repr__(self) -> str:
        if not self.success:
            return f"RANSResult(success=False, error={self.error!r})"
        conv = "converged" if self.converged else "NOT converged"
        return (
            f"RANSResult(α={self.alpha:.2f}°, Re={self.reynolds:.0e}, M={self.mach:.2f}, "
            f"CL={self.cl:.4f}, CD={self.cd:.5f}, CM={self.cm:.4f}, {conv})"
        )


@dataclass
class RANSPolarResult:
    """Result of a RANS polar sweep via Flow360."""

    airfoil_name: str
    reynolds: float
    mach: float
    results: list[RANSResult] = field(default_factory=list)

    @property
    def converged(self) -> list[RANSResult]:
        return [r for r in self.results if r.converged]

    @property
    def alpha(self) -> list[float]:
        return [r.alpha for r in self.converged]

    @property
    def cl(self) -> list[float]:
        return [r.cl for r in self.converged]

    @property
    def cd(self) -> list[float]:
        return [r.cd for r in self.converged]

    @property
    def cm(self) -> list[float]:
        return [r.cm for r in self.converged]

    @property
    def cl_max(self) -> float | None:
        vals = self.cl
        return max(vals) if vals else None

    @property
    def cd_min(self) -> float | None:
        vals = self.cd
        return min(vals) if vals else None

    @property
    def ld_max(self) -> float | None:
        lds = [r.ld for r in self.converged if r.ld is not None]
        return max(lds) if lds else None

    def to_dict(self) -> dict:
        return {
            "airfoil_name": self.airfoil_name,
            "reynolds": self.reynolds,
            "mach": self.mach,
            "alpha": self.alpha,
            "cl": self.cl,
            "cd": self.cd,
            "cm": self.cm,
        }

    def __repr__(self) -> str:
        n = len(self.converged)
        total = len(self.results)
        return (
            f"RANSPolarResult({self.airfoil_name}, Re={self.reynolds:.0e}, M={self.mach:.2f}, "
            f"{n}/{total} converged)"
        )
