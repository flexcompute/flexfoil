"""
Matrix sweep — sweep alpha x Reynolds number for a flapped airfoil.

Demonstrates the full parameter matrix capability: every combination
of alpha and Re is solved, with results stored in the local database.

    python examples/12_matrix_sweep.py
"""

import time
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import flexfoil

foil = flexfoil.naca("2412")
flapped = foil.with_flap(hinge_x=0.75, deflection=10)

alphas = (-2, 12, 2.0)
reynolds = [2e5, 5e5, 1e6, 3e6]

fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=("CL vs α", "CL vs CD"),
)

t0 = time.time()
for Re in reynolds:
    polar = flapped.polar(alpha=alphas, Re=Re)
    label = f"Re = {Re:.0e}"

    fig.add_trace(go.Scatter(
        x=polar.alpha, y=polar.cl,
        mode="lines+markers", name=label, marker=dict(size=5),
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=polar.cd, y=polar.cl,
        mode="lines+markers", name=label, showlegend=False, marker=dict(size=5),
    ), row=1, col=2)

    print(f"  {label}: {len(polar.converged)}/{len(polar.results)} converged")

elapsed = time.time() - t0
total = len(reynolds) * len(polar.results)
print(f"\n{total} solves in {elapsed:.1f}s")

fig.update_xaxes(title_text="α (°)", row=1, col=1)
fig.update_yaxes(title_text="CL", row=1, col=1)
fig.update_xaxes(title_text="CD", row=1, col=2)
fig.update_yaxes(title_text="CL", row=1, col=2)
fig.update_layout(
    title=f"{flapped.name} — Matrix sweep",
    height=500, width=1000, template="plotly_white",
)
fig.show()

print(f"\nAll runs cached — open the web UI to explore:")
print(f"  flexfoil serve")
