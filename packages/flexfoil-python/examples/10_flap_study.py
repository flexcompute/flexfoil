"""
Flap study — how does flap deflection affect CL and CD?

    python examples/10_flap_study.py
"""

import plotly.graph_objects as go
from plotly.subplots import make_subplots
import flexfoil

foil = flexfoil.naca("2412")

deflections = [0, 5, 10, 15, 20]
alpha_range = (-4, 14, 1.0)

fig = make_subplots(
    rows=1, cols=2,
    subplot_titles=("CL vs α", "CL vs CD (drag polar)"),
)

for defl in deflections:
    if defl == 0:
        f = foil
    else:
        f = foil.with_flap(hinge_x=0.75, deflection=defl)

    polar = f.polar(alpha=alpha_range, Re=1e6)
    label = f"δ = {defl}°"

    fig.add_trace(go.Scatter(
        x=polar.alpha, y=polar.cl,
        mode="lines+markers", name=label, marker=dict(size=4),
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=polar.cd, y=polar.cl,
        mode="lines+markers", name=label, showlegend=False, marker=dict(size=4),
    ), row=1, col=2)

    cl_max = max(polar.cl) if polar.cl else 0
    print(f"  {label:>8}: CL_max = {cl_max:.3f}")

fig.update_xaxes(title_text="α (°)", row=1, col=1)
fig.update_yaxes(title_text="CL", row=1, col=1)
fig.update_xaxes(title_text="CD", row=1, col=2)
fig.update_yaxes(title_text="CL", row=1, col=2)
fig.update_layout(
    title=f"{foil.name} — Flap at 75% chord, Re=1e6",
    height=500, width=1000, template="plotly_white",
)
fig.show()
