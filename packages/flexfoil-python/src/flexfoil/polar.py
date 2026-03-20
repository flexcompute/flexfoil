"""Polar sweep result container with plotting, export, and aggregate statistics."""

from __future__ import annotations

import statistics
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from flexfoil.airfoil import SolveResult


def _iqr_filter(values: list[float], multiplier: float = 1.5) -> list[int]:
    """Return indices of values that pass the 1.5×IQR fence."""
    if len(values) < 4:
        return list(range(len(values)))
    s = sorted(values)
    n = len(s)
    q1 = s[int(n * 0.25)]
    q3 = s[int(n * 0.75)]
    iqr = q3 - q1
    if iqr == 0:
        return list(range(len(values)))
    lo, hi = q1 - multiplier * iqr, q3 + multiplier * iqr
    return [i for i, v in enumerate(values) if lo <= v <= hi]


@dataclass
class PolarResult:
    """Container for a polar sweep (CL, CD, CM vs alpha)."""

    airfoil_name: str
    reynolds: float
    mach: float
    ncrit: float
    results: list[SolveResult] = field(default_factory=list)

    # ── Column accessors (converged points only) ──────────────

    @property
    def converged(self) -> list[SolveResult]:
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
    def ld(self) -> list[float | None]:
        return [r.ld for r in self.converged]

    # ── Aggregate properties ──────────────────────────────────

    @property
    def cl_max(self) -> float | None:
        """Maximum lift coefficient from converged points."""
        vals = self.cl
        return max(vals) if vals else None

    @property
    def cl_min(self) -> float | None:
        """Minimum lift coefficient from converged points."""
        vals = self.cl
        return min(vals) if vals else None

    @property
    def cd_min(self) -> float | None:
        """Minimum drag coefficient from converged points."""
        vals = self.cd
        return min(vals) if vals else None

    @property
    def ld_max(self) -> float | None:
        """Maximum lift-to-drag ratio from converged points."""
        finite = [v for v in self.ld if v is not None]
        return max(finite) if finite else None

    @property
    def alpha_stall(self) -> float | None:
        """Angle of attack at CL_max."""
        return self.argmax("cl", "alpha")

    @property
    def alpha_at_ld_max(self) -> float | None:
        """Angle of attack at L/D_max."""
        return self.argmax("ld", "alpha")

    @property
    def alpha_at_cd_min(self) -> float | None:
        """Angle of attack at CD_min."""
        return self.argmin("cd", "alpha")

    # ── Generic column aggregation ────────────────────────────

    def _get_column(
        self, col: str, *, filter_outliers: bool = False,
    ) -> tuple[list[float], list[SolveResult]]:
        """Retrieve a numeric column from converged results.

        Returns (values, corresponding_results) with None entries removed.
        When *filter_outliers* is True, applies IQR filtering.
        """
        results = self.converged
        raw = [getattr(r, col) for r in results]
        pairs = [(v, r) for v, r in zip(raw, results) if v is not None]
        if not pairs:
            return [], []

        vals, res = zip(*pairs)
        vals_list = list(vals)
        res_list = list(res)

        if filter_outliers and len(vals_list) >= 4:
            keep = _iqr_filter(vals_list)
            vals_list = [vals_list[i] for i in keep]
            res_list = [res_list[i] for i in keep]

        return vals_list, res_list

    def column_max(
        self, col: str, *, filter_outliers: bool = False,
    ) -> float | None:
        """Maximum value of any numeric column (cl, cd, cm, ld, alpha)."""
        vals, _ = self._get_column(col, filter_outliers=filter_outliers)
        return max(vals) if vals else None

    def column_min(
        self, col: str, *, filter_outliers: bool = False,
    ) -> float | None:
        """Minimum value of any numeric column."""
        vals, _ = self._get_column(col, filter_outliers=filter_outliers)
        return min(vals) if vals else None

    def column_mean(
        self, col: str, *, filter_outliers: bool = False,
    ) -> float | None:
        """Arithmetic mean of any numeric column."""
        vals, _ = self._get_column(col, filter_outliers=filter_outliers)
        return statistics.mean(vals) if vals else None

    def column_median(
        self, col: str, *, filter_outliers: bool = False,
    ) -> float | None:
        """Median of any numeric column."""
        vals, _ = self._get_column(col, filter_outliers=filter_outliers)
        return statistics.median(vals) if vals else None

    def column_stdev(
        self, col: str, *, filter_outliers: bool = False,
    ) -> float | None:
        """Sample standard deviation of any numeric column."""
        vals, _ = self._get_column(col, filter_outliers=filter_outliers)
        return statistics.stdev(vals) if len(vals) >= 2 else None

    def argmax(
        self,
        target: str,
        return_col: str,
        *,
        filter_outliers: bool = False,
    ) -> float | None:
        """Value of *return_col* at the row where *target* is maximized.

        Example: ``polar.argmax('cl', 'alpha')`` returns alpha_stall.
        """
        vals, res = self._get_column(target, filter_outliers=filter_outliers)
        if not vals:
            return None
        idx = max(range(len(vals)), key=lambda i: vals[i])
        return getattr(res[idx], return_col, None)

    def argmin(
        self,
        target: str,
        return_col: str,
        *,
        filter_outliers: bool = False,
    ) -> float | None:
        """Value of *return_col* at the row where *target* is minimized.

        Example: ``polar.argmin('cd', 'alpha')`` returns alpha at CD_min.
        """
        vals, res = self._get_column(target, filter_outliers=filter_outliers)
        if not vals:
            return None
        idx = min(range(len(vals)), key=lambda i: vals[i])
        return getattr(res[idx], return_col, None)

    # ── Export ─────────────────────────────────────────────────

    def to_dict(self, *, summary: bool = False) -> dict:
        d = {
            "alpha": self.alpha,
            "cl": self.cl,
            "cd": self.cd,
            "cm": self.cm,
            "ld": self.ld,
        }
        if summary:
            d["_summary"] = {
                "cl_max": self.cl_max,
                "cl_min": self.cl_min,
                "cd_min": self.cd_min,
                "ld_max": self.ld_max,
                "alpha_stall": self.alpha_stall,
                "alpha_at_ld_max": self.alpha_at_ld_max,
                "alpha_at_cd_min": self.alpha_at_cd_min,
            }
        return d

    def to_dataframe(self, *, summary: bool = False):
        """Return a pandas DataFrame (requires pandas)."""
        import pandas as pd

        df = pd.DataFrame(self.to_dict())
        if summary:
            df.attrs["cl_max"] = self.cl_max
            df.attrs["cl_min"] = self.cl_min
            df.attrs["cd_min"] = self.cd_min
            df.attrs["ld_max"] = self.ld_max
            df.attrs["alpha_stall"] = self.alpha_stall
            df.attrs["alpha_at_ld_max"] = self.alpha_at_ld_max
            df.attrs["alpha_at_cd_min"] = self.alpha_at_cd_min
        return df

    def summary(self) -> dict[str, float | None]:
        """Return a dict of all aggregate statistics."""
        return {
            "cl_max": self.cl_max,
            "cl_min": self.cl_min,
            "cd_min": self.cd_min,
            "ld_max": self.ld_max,
            "alpha_stall": self.alpha_stall,
            "alpha_at_ld_max": self.alpha_at_ld_max,
            "alpha_at_cd_min": self.alpha_at_cd_min,
        }

    # ── Plotting ──────────────────────────────────────────────

    def plot(self, *, show: bool = True, backend: str = "plotly"):
        """Plot CL vs alpha, CD vs alpha, CL vs CD, and CM vs alpha.

        Parameters
        ----------
        show : if True, display the figure immediately
        backend : "plotly" (default) or "matplotlib"
        """
        if backend == "matplotlib":
            return self._plot_matplotlib(show=show)
        return self._plot_plotly(show=show)

    def _plot_plotly(self, *, show: bool = True):
        from plotly.subplots import make_subplots
        import plotly.graph_objects as go

        a, cl, cd, cm = self.alpha, self.cl, self.cd, self.cm
        title = f"{self.airfoil_name}  Re={self.reynolds:.2e}  M={self.mach}"

        fig = make_subplots(
            rows=2, cols=2,
            subplot_titles=("CL vs α", "CD vs α", "CL vs CD", "CM vs α"),
            horizontal_spacing=0.10,
            vertical_spacing=0.12,
        )

        line_style = dict(color="#636EFA")
        marker_style = dict(size=5)

        fig.add_trace(go.Scatter(
            x=a, y=cl, mode="lines+markers",
            line=line_style, marker=marker_style, showlegend=False,
        ), row=1, col=1)

        fig.add_trace(go.Scatter(
            x=a, y=cd, mode="lines+markers",
            line=line_style, marker=marker_style, showlegend=False,
        ), row=1, col=2)

        fig.add_trace(go.Scatter(
            x=cd, y=cl, mode="lines+markers",
            line=line_style, marker=marker_style, showlegend=False,
        ), row=2, col=1)

        fig.add_trace(go.Scatter(
            x=a, y=cm, mode="lines+markers",
            line=line_style, marker=marker_style, showlegend=False,
        ), row=2, col=2)

        fig.update_xaxes(title_text="α (°)", row=1, col=1)
        fig.update_yaxes(title_text="CL", row=1, col=1)
        fig.update_xaxes(title_text="α (°)", row=1, col=2)
        fig.update_yaxes(title_text="CD", row=1, col=2)
        fig.update_xaxes(title_text="CD", row=2, col=1)
        fig.update_yaxes(title_text="CL", row=2, col=1)
        fig.update_xaxes(title_text="α (°)", row=2, col=2)
        fig.update_yaxes(title_text="CM", row=2, col=2)

        fig.update_layout(
            title_text=title,
            height=700,
            width=900,
            template="plotly_white",
        )

        if show:
            fig.show()
        return fig

    def _plot_matplotlib(self, *, show: bool = True):
        import matplotlib.pyplot as plt

        fig, axes = plt.subplots(2, 2, figsize=(10, 8))
        fig.suptitle(
            f"{self.airfoil_name}  Re={self.reynolds:.2e}  M={self.mach}",
            fontsize=13,
        )

        a, cl, cd, cm = self.alpha, self.cl, self.cd, self.cm

        axes[0, 0].plot(a, cl, ".-")
        axes[0, 0].set_xlabel("α (°)")
        axes[0, 0].set_ylabel("CL")
        axes[0, 0].grid(True, alpha=0.3)

        axes[0, 1].plot(a, cd, ".-")
        axes[0, 1].set_xlabel("α (°)")
        axes[0, 1].set_ylabel("CD")
        axes[0, 1].grid(True, alpha=0.3)

        axes[1, 0].plot(cd, cl, ".-")
        axes[1, 0].set_xlabel("CD")
        axes[1, 0].set_ylabel("CL")
        axes[1, 0].grid(True, alpha=0.3)

        axes[1, 1].plot(a, cm, ".-")
        axes[1, 1].set_xlabel("α (°)")
        axes[1, 1].set_ylabel("CM")
        axes[1, 1].grid(True, alpha=0.3)

        fig.tight_layout()
        if show:
            plt.show()
        return fig

    def __repr__(self) -> str:
        n_conv = len(self.converged)
        n_total = len(self.results)
        parts = [
            f"PolarResult({self.airfoil_name!r}, Re={self.reynolds:.0e}, "
            f"{n_conv}/{n_total} converged",
        ]
        if self.cl_max is not None:
            parts.append(f", CLmax={self.cl_max:.4f}")
        if self.alpha_stall is not None:
            parts.append(f", α_stall={self.alpha_stall:.1f}°")
        if self.ld_max is not None:
            parts.append(f", L/D_max={self.ld_max:.1f}")
        parts.append(")")
        return "".join(parts)
