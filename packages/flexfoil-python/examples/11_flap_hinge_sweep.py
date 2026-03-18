"""
Flap hinge position sweep — find the best hinge location.

Sweeps hinge x/c from 60% to 85% at a fixed 10° flap deflection.

    python examples/11_flap_hinge_sweep.py
"""

import plotly.graph_objects as go
import flexfoil

foil = flexfoil.naca("2412")
alpha = 5.0
Re = 1e6

hinge_positions = [0.60, 0.65, 0.70, 0.75, 0.80, 0.85]
cls, cds, lds = [], [], []

print(f"{'hinge x/c':>10}  {'CL':>8}  {'CD':>8}  {'L/D':>8}")
print("-" * 42)

for hx in hinge_positions:
    flapped = foil.with_flap(hinge_x=hx, deflection=10.0)
    r = flapped.solve(alpha=alpha, Re=Re, store=False)
    cls.append(r.cl)
    cds.append(r.cd)
    lds.append(r.ld or 0)
    print(f"{hx:10.2f}  {r.cl:8.4f}  {r.cd:8.5f}  {r.ld or 0:8.1f}")

fig = go.Figure()
fig.add_trace(go.Scatter(x=hinge_positions, y=lds, mode="lines+markers", name="L/D"))
fig.update_layout(
    title=f"{foil.name} — 10° flap, α={alpha}°, Re={Re:.0e}",
    xaxis_title="Hinge x/c",
    yaxis_title="L/D",
    template="plotly_white",
)
fig.show()
