#!/usr/bin/env python3
"""
Interactive comparison dashboard: RustFoil vs XFOIL
Requires: pip install plotly dash pandas
"""

import json
from pathlib import Path
import numpy as np

# Required packages: pip install plotly dash
import plotly.graph_objects as go
from plotly.subplots import make_subplots
import dash
from dash import dcc, html, Input, Output

# Load data
DATA_DIR = Path(__file__).parent / "comparison_data"

def load_data():
    with open(DATA_DIR / "xfoil_reference.json") as f:
        xfoil_data = json.load(f)
    with open(DATA_DIR / "rustfoil_data.json") as f:
        rustfoil_data = json.load(f)
    return xfoil_data, rustfoil_data

xfoil_data, rustfoil_data = load_data()

# Extract available foils and alphas
FOILS = list(rustfoil_data.keys())
ALPHAS = sorted([int(a) for a in list(rustfoil_data[FOILS[0]]["alphas"].keys())])

# Calculate lift coefficients from XFOIL (integrate Cp)
def calc_cl_from_cp(x, y, cp):
    """
    Calculate CL by integrating Cp around airfoil.
    Data is ordered: upper TE → LE → lower TE (open contour)
    CL = ∫ Cp dx
    """
    n = len(x)
    cl = 0.0
    for i in range(n - 1):
        dx = x[i+1] - x[i]
        cp_avg = 0.5 * (cp[i] + cp[i+1])
        cl += cp_avg * dx
    return cl

# Build polar data
def build_polar_data():
    """Build CL vs alpha data for all foils."""
    polars = {}
    for foil in FOILS:
        polars[foil] = {
            "xfoil": {"alpha": [], "cl": []},
            "rustfoil": {"alpha": [], "cl": []}
        }
        
        for alpha in ALPHAS:
            alpha_str = str(alpha)
            
            # XFOIL
            if foil in xfoil_data and alpha_str in xfoil_data[foil]["alphas"]:
                xf = xfoil_data[foil]["alphas"][alpha_str]
                cl = calc_cl_from_cp(xf["x"], xf["y"], xf["cp"])
                polars[foil]["xfoil"]["alpha"].append(alpha)
                polars[foil]["xfoil"]["cl"].append(cl)
            
            # RustFoil
            if foil in rustfoil_data and alpha_str in rustfoil_data[foil]["alphas"]:
                rf = rustfoil_data[foil]["alphas"][alpha_str]
                polars[foil]["rustfoil"]["alpha"].append(alpha)
                polars[foil]["rustfoil"]["cl"].append(rf["cl"])
    
    return polars

polars = build_polar_data()

# Create Dash app
app = dash.Dash(__name__)

app.layout = html.Div([
    html.H1("RustFoil vs XFOIL Comparison", style={"textAlign": "center"}),
    
    html.Div([
        html.Div([
            html.Label("Airfoil:"),
            dcc.Dropdown(
                id="foil-dropdown",
                options=[{"label": f.upper(), "value": f} for f in FOILS],
                value=FOILS[0],
                style={"width": "200px"}
            ),
        ], style={"display": "inline-block", "marginRight": "20px"}),
        
        html.Div([
            html.Label("Angle of Attack:"),
            dcc.Dropdown(
                id="alpha-dropdown",
                options=[{"label": f"{a}°", "value": a} for a in ALPHAS],
                value=0,
                style={"width": "100px"}
            ),
        ], style={"display": "inline-block", "marginRight": "20px"}),
        
        html.Div([
            html.Label("View:"),
            dcc.RadioItems(
                id="view-radio",
                options=[
                    {"label": "Distributions", "value": "dist"},
                    {"label": "Polar (CL vs α)", "value": "polar"},
                    {"label": "All Foils Polar", "value": "all_polar"},
                ],
                value="dist",
                inline=True
            ),
        ], style={"display": "inline-block"}),
    ], style={"textAlign": "center", "marginBottom": "20px"}),
    
    dcc.Graph(id="main-graph", style={"height": "80vh"}),
    
    html.Div(id="stats-output", style={"textAlign": "center", "marginTop": "10px", "fontFamily": "monospace"})
])

@app.callback(
    [Output("main-graph", "figure"), Output("stats-output", "children")],
    [Input("foil-dropdown", "value"), Input("alpha-dropdown", "value"), Input("view-radio", "value")]
)
def update_graph(foil, alpha, view):
    if view == "polar":
        return create_polar_plot(foil), ""
    elif view == "all_polar":
        return create_all_polar_plot(), ""
    else:
        return create_distribution_plot(foil, alpha)

def create_distribution_plot(foil, alpha):
    """Create Cp and Gamma distribution plots."""
    alpha_str = str(alpha)
    
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            f"Pressure Coefficient (Cp) - {foil.upper()} α={alpha}°",
            f"Vortex Strength (γ) vs x - {foil.upper()} α={alpha}°",
            f"Airfoil Geometry - {foil.upper()}",
            f"γ vs Arc Length (s) - {foil.upper()} α={alpha}°"
        ),
        vertical_spacing=0.12,
        horizontal_spacing=0.08
    )
    
    stats = []
    
    # Get data
    xf = xfoil_data.get(foil, {}).get("alphas", {}).get(alpha_str, {})
    rf = rustfoil_data.get(foil, {}).get("alphas", {}).get(alpha_str, {})
    
    if not xf or not rf:
        return fig, "No data available"
    
    # Cp vs x (row 1, col 1) - inverted y-axis for Cp
    fig.add_trace(go.Scatter(
        x=xf["x"], y=xf["cp"], mode="lines", name="XFOIL Cp",
        line=dict(color="blue", width=2)
    ), row=1, col=1)
    fig.add_trace(go.Scatter(
        x=rf["x"], y=rf["cp"], mode="lines", name="RustFoil Cp",
        line=dict(color="red", width=2, dash="dash")
    ), row=1, col=1)
    fig.update_yaxes(autorange="reversed", title_text="Cp", row=1, col=1)
    fig.update_xaxes(title_text="x/c", row=1, col=1)
    
    # Gamma vs x (row 1, col 2)
    fig.add_trace(go.Scatter(
        x=xf["x"], y=xf["gamma"], mode="lines", name="XFOIL γ",
        line=dict(color="blue", width=2), showlegend=False
    ), row=1, col=2)
    fig.add_trace(go.Scatter(
        x=rf["x"], y=rf["gamma"], mode="lines", name="RustFoil γ",
        line=dict(color="red", width=2, dash="dash"), showlegend=False
    ), row=1, col=2)
    fig.update_yaxes(title_text="γ (Ue/V∞)", row=1, col=2)
    fig.update_xaxes(title_text="x/c", row=1, col=2)
    
    # Airfoil geometry (row 2, col 1)
    fig.add_trace(go.Scatter(
        x=xf["x"], y=xf["y"], mode="lines", name="Geometry",
        line=dict(color="black", width=2), showlegend=False, fill="toself", fillcolor="rgba(200,200,200,0.3)"
    ), row=2, col=1)
    fig.update_yaxes(scaleanchor="x3", scaleratio=1, title_text="y/c", row=2, col=1)
    fig.update_xaxes(title_text="x/c", row=2, col=1)
    
    # Gamma vs s (row 2, col 2)
    fig.add_trace(go.Scatter(
        x=xf["s"], y=xf["gamma"], mode="lines", name="XFOIL γ(s)",
        line=dict(color="blue", width=2), showlegend=False
    ), row=2, col=2)
    fig.add_trace(go.Scatter(
        x=rf["s"], y=rf["gamma"], mode="lines", name="RustFoil γ(s)",
        line=dict(color="red", width=2, dash="dash"), showlegend=False
    ), row=2, col=2)
    fig.update_yaxes(title_text="γ (Ue/V∞)", row=2, col=2)
    fig.update_xaxes(title_text="Arc length s", row=2, col=2)
    
    fig.update_layout(
        legend=dict(x=0.5, y=1.02, xanchor="center", orientation="h"),
        margin=dict(t=80, b=40, l=60, r=40)
    )
    
    # Calculate statistics
    gamma_err = np.array(rf["gamma"]) - np.array(xf["gamma"])
    cp_err = np.array(rf["cp"]) - np.array(xf["cp"])
    
    cl_xfoil = calc_cl_from_cp(xf["x"], xf["y"], xf["cp"])
    cl_rust = rf["cl"]
    
    stats_text = (
        f"γ max error: {np.max(np.abs(gamma_err)):.2e} | "
        f"γ RMS: {np.sqrt(np.mean(gamma_err**2)):.2e} | "
        f"Cp max error: {np.max(np.abs(cp_err)):.2e} | "
        f"CL XFOIL: {cl_xfoil:.4f} | CL RustFoil: {cl_rust:.4f} | "
        f"CL error: {abs(cl_rust - cl_xfoil):.2e}"
    )
    
    return fig, stats_text

def create_polar_plot(foil):
    """Create CL vs alpha polar plot for single foil."""
    fig = go.Figure()
    
    polar = polars[foil]
    
    fig.add_trace(go.Scatter(
        x=polar["xfoil"]["alpha"], y=polar["xfoil"]["cl"],
        mode="lines+markers", name="XFOIL",
        line=dict(color="blue", width=2), marker=dict(size=8)
    ))
    fig.add_trace(go.Scatter(
        x=polar["rustfoil"]["alpha"], y=polar["rustfoil"]["cl"],
        mode="lines+markers", name="RustFoil",
        line=dict(color="red", width=2, dash="dash"), marker=dict(size=8, symbol="x")
    ))
    
    # Add thin airfoil theory line
    alpha_theory = np.linspace(-6, 12, 50)
    cl_theory = 2 * np.pi * np.deg2rad(alpha_theory)
    fig.add_trace(go.Scatter(
        x=alpha_theory, y=cl_theory,
        mode="lines", name="2πα (thin airfoil)",
        line=dict(color="gray", width=1, dash="dot")
    ))
    
    # Calculate lift curve slope
    xf_alpha = np.array(polar["xfoil"]["alpha"])
    xf_cl = np.array(polar["xfoil"]["cl"])
    rf_alpha = np.array(polar["rustfoil"]["alpha"])
    rf_cl = np.array(polar["rustfoil"]["cl"])
    
    if len(xf_alpha) > 1:
        xf_slope = np.polyfit(np.deg2rad(xf_alpha), xf_cl, 1)[0]
        rf_slope = np.polyfit(np.deg2rad(rf_alpha), rf_cl, 1)[0]
        
        fig.add_annotation(
            x=0.02, y=0.98, xref="paper", yref="paper",
            text=f"Lift curve slope:<br>XFOIL: {xf_slope:.4f}/rad<br>RustFoil: {rf_slope:.4f}/rad<br>Theory (2π): {2*np.pi:.4f}/rad",
            showarrow=False, align="left",
            bgcolor="rgba(255,255,255,0.8)", bordercolor="black", borderwidth=1
        )
    
    fig.update_layout(
        title=f"Lift Polar - {foil.upper()}",
        xaxis_title="Angle of Attack (°)",
        yaxis_title="Lift Coefficient (CL)",
        legend=dict(x=0.7, y=0.1),
        hovermode="x unified"
    )
    
    return fig

def create_all_polar_plot():
    """Create CL vs alpha for all foils."""
    fig = go.Figure()
    
    colors = {"naca0012": "blue", "naca2412": "green", "naca4412": "red"}
    
    for foil in FOILS:
        polar = polars[foil]
        color = colors.get(foil, "gray")
        
        fig.add_trace(go.Scatter(
            x=polar["xfoil"]["alpha"], y=polar["xfoil"]["cl"],
            mode="lines+markers", name=f"{foil.upper()} XFOIL",
            line=dict(color=color, width=2), marker=dict(size=6)
        ))
        fig.add_trace(go.Scatter(
            x=polar["rustfoil"]["alpha"], y=polar["rustfoil"]["cl"],
            mode="lines+markers", name=f"{foil.upper()} RustFoil",
            line=dict(color=color, width=2, dash="dash"), marker=dict(size=6, symbol="x")
        ))
    
    fig.update_layout(
        title="Lift Polars - All Foils Comparison",
        xaxis_title="Angle of Attack (°)",
        yaxis_title="Lift Coefficient (CL)",
        legend=dict(x=0.01, y=0.99),
        hovermode="x unified"
    )
    
    return fig

if __name__ == "__main__":
    print("\n" + "="*60)
    print("RustFoil vs XFOIL Comparison Dashboard")
    print("="*60)
    print("\nOpen http://127.0.0.1:8050 in your browser")
    print("Press Ctrl+C to stop\n")
    
    app.run(debug=True, port=8050)
