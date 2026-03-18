"""Polar sweep result container with plotting and export."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from flexfoil.airfoil import SolveResult


@dataclass
class PolarResult:
    """Container for a polar sweep (CL, CD, CM vs alpha)."""

    airfoil_name: str
    reynolds: float
    mach: float
    ncrit: float
    results: list[SolveResult] = field(default_factory=list)

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

    def to_dict(self) -> dict:
        return {
            "alpha": self.alpha,
            "cl": self.cl,
            "cd": self.cd,
            "cm": self.cm,
            "ld": self.ld,
        }

    def to_dataframe(self):
        """Return a pandas DataFrame (requires pandas)."""
        import pandas as pd
        return pd.DataFrame(self.to_dict())

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
        return (
            f"PolarResult({self.airfoil_name!r}, Re={self.reynolds:.0e}, "
            f"{n_conv}/{n_total} converged)"
        )
