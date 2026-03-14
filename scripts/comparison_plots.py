#!/usr/bin/env python3
"""
Generate static comparison plots: RustFoil vs XFOIL
Creates an HTML file with interactive Plotly plots (no server required)
"""

import json
from pathlib import Path
import numpy as np

try:
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots
except ImportError:
    print("Installing plotly...")
    import subprocess
    subprocess.run(["pip", "install", "plotly"])
    import plotly.graph_objects as go
    from plotly.subplots import make_subplots

# Load data
DATA_DIR = Path(__file__).parent / "comparison_data"

with open(DATA_DIR / "xfoil_reference.json") as f:
    xfoil_data = json.load(f)
with open(DATA_DIR / "rustfoil_data.json") as f:
    rustfoil_data = json.load(f)

FOILS = list(rustfoil_data.keys())
ALPHAS = sorted([int(a) for a in list(rustfoil_data[FOILS[0]]["alphas"].keys())])

def calc_cl_from_cp(x, y, cp, alpha_deg=0.0):
    """
    Calculate CL by integrating Cp around airfoil.
    
    Uses XFOIL's method: closed contour + wind-axis rotation
    CL = ∮ Cp * dx_wind, where dx_wind = dx*cos(α) + dy*sin(α)
    """
    import math
    n = len(x)
    alpha_rad = math.radians(alpha_deg)
    cosa, sina = math.cos(alpha_rad), math.sin(alpha_rad)
    
    cl = 0.0
    for i in range(n):
        ip = (i + 1) % n  # Closed contour
        dx = x[ip] - x[i]
        dy = y[ip] - y[i]
        dx_wind = dx * cosa + dy * sina  # Wind-axis rotation
        cp_avg = 0.5 * (cp[i] + cp[ip])
        cl += cp_avg * dx_wind
    return cl

# Create comprehensive comparison figure
def create_full_comparison():
    """Create a comprehensive multi-tab comparison."""
    
    # Build summary statistics
    print("\n" + "="*80)
    print("RUSTFOIL vs XFOIL COMPARISON SUMMARY")
    print("="*80)
    
    all_stats = []
    
    for foil in FOILS:
        print(f"\n{foil.upper()}")
        print("-" * 40)
        print(f"{'α':>6} {'CL_xfoil':>10} {'CL_rust':>10} {'CL_err':>10} {'γ_max_err':>12} {'Cp_max_err':>12}")
        
        for alpha in ALPHAS:
            alpha_str = str(alpha)
            xf = xfoil_data[foil]["alphas"].get(alpha_str, {})
            rf = rustfoil_data[foil]["alphas"].get(alpha_str, {})
            
            if not xf or not rf:
                continue
            
            cl_xfoil = calc_cl_from_cp(xf["x"], xf["y"], xf["cp"], alpha)
            cl_rust = rf["cl"]
            cl_err = abs(cl_rust - cl_xfoil)
            
            gamma_err = np.array(rf["gamma"]) - np.array(xf["gamma"])
            gamma_max = np.max(np.abs(gamma_err))
            
            cp_err = np.array(rf["cp"]) - np.array(xf["cp"])
            cp_max = np.max(np.abs(cp_err))
            
            print(f"{alpha:>6}° {cl_xfoil:>10.4f} {cl_rust:>10.4f} {cl_err:>10.2e} {gamma_max:>12.2e} {cp_max:>12.2e}")
            
            all_stats.append({
                "foil": foil, "alpha": alpha,
                "cl_xfoil": cl_xfoil, "cl_rust": cl_rust,
                "gamma_max_err": gamma_max, "cp_max_err": cp_max
            })
    
    # Create HTML with multiple figures
    figures = []
    
    # 1. All polars comparison
    fig_polar = go.Figure()
    colors = {"naca0012": "blue", "naca2412": "green", "naca4412": "red"}
    
    for foil in FOILS:
        color = colors.get(foil, "gray")
        
        xf_alpha, xf_cl = [], []
        rf_alpha, rf_cl = [], []
        
        for alpha in ALPHAS:
            alpha_str = str(alpha)
            xf = xfoil_data[foil]["alphas"].get(alpha_str, {})
            rf = rustfoil_data[foil]["alphas"].get(alpha_str, {})
            
            if xf:
                xf_alpha.append(alpha)
                xf_cl.append(calc_cl_from_cp(xf["x"], xf["y"], xf["cp"], alpha))
            if rf:
                rf_alpha.append(alpha)
                rf_cl.append(rf["cl"])
        
        fig_polar.add_trace(go.Scatter(
            x=xf_alpha, y=xf_cl, mode="lines+markers",
            name=f"{foil.upper()} XFOIL", line=dict(color=color, width=2),
            marker=dict(size=8)
        ))
        fig_polar.add_trace(go.Scatter(
            x=rf_alpha, y=rf_cl, mode="lines+markers",
            name=f"{foil.upper()} RustFoil", line=dict(color=color, width=2, dash="dash"),
            marker=dict(size=8, symbol="x")
        ))
    
    fig_polar.update_layout(
        title="Lift Polars: RustFoil vs XFOIL (All Foils)",
        xaxis_title="Angle of Attack (°)",
        yaxis_title="Lift Coefficient (CL)",
        height=600
    )
    figures.append(("Lift Polars", fig_polar))
    
    # 2. Cp comparison at different alphas for each foil
    for foil in FOILS:
        fig_cp = make_subplots(
            rows=2, cols=4,
            subplot_titles=[f"α={a}°" for a in ALPHAS],
            vertical_spacing=0.12
        )
        
        for idx, alpha in enumerate(ALPHAS):
            row = idx // 4 + 1
            col = idx % 4 + 1
            alpha_str = str(alpha)
            
            xf = xfoil_data[foil]["alphas"].get(alpha_str, {})
            rf = rustfoil_data[foil]["alphas"].get(alpha_str, {})
            
            if xf and rf:
                fig_cp.add_trace(go.Scatter(
                    x=xf["x"], y=xf["cp"], mode="lines", name="XFOIL",
                    line=dict(color="blue", width=1.5),
                    showlegend=(idx == 0)
                ), row=row, col=col)
                fig_cp.add_trace(go.Scatter(
                    x=rf["x"], y=rf["cp"], mode="lines", name="RustFoil",
                    line=dict(color="red", width=1.5, dash="dash"),
                    showlegend=(idx == 0)
                ), row=row, col=col)
                
                fig_cp.update_yaxes(autorange="reversed", row=row, col=col)
        
        fig_cp.update_layout(
            title=f"Pressure Distributions: {foil.upper()}",
            height=600
        )
        figures.append((f"Cp - {foil.upper()}", fig_cp))
    
    # 3. Gamma comparison
    for foil in FOILS:
        fig_gamma = make_subplots(
            rows=2, cols=4,
            subplot_titles=[f"α={a}°" for a in ALPHAS],
            vertical_spacing=0.12
        )
        
        for idx, alpha in enumerate(ALPHAS):
            row = idx // 4 + 1
            col = idx % 4 + 1
            alpha_str = str(alpha)
            
            xf = xfoil_data[foil]["alphas"].get(alpha_str, {})
            rf = rustfoil_data[foil]["alphas"].get(alpha_str, {})
            
            if xf and rf:
                fig_gamma.add_trace(go.Scatter(
                    x=xf["s"], y=xf["gamma"], mode="lines", name="XFOIL",
                    line=dict(color="blue", width=1.5),
                    showlegend=(idx == 0)
                ), row=row, col=col)
                fig_gamma.add_trace(go.Scatter(
                    x=rf["s"], y=rf["gamma"], mode="lines", name="RustFoil",
                    line=dict(color="red", width=1.5, dash="dash"),
                    showlegend=(idx == 0)
                ), row=row, col=col)
        
        fig_gamma.update_layout(
            title=f"Vortex Strength (γ) vs Arc Length: {foil.upper()}",
            height=600
        )
        figures.append((f"Gamma - {foil.upper()}", fig_gamma))
    
    # 4. Error analysis
    fig_error = make_subplots(
        rows=1, cols=2,
        subplot_titles=["Max Gamma Error", "Max Cp Error"]
    )
    
    for foil in FOILS:
        foil_stats = [s for s in all_stats if s["foil"] == foil]
        alphas = [s["alpha"] for s in foil_stats]
        gamma_errs = [s["gamma_max_err"] for s in foil_stats]
        cp_errs = [s["cp_max_err"] for s in foil_stats]
        
        color = colors.get(foil, "gray")
        
        fig_error.add_trace(go.Scatter(
            x=alphas, y=gamma_errs, mode="lines+markers",
            name=foil.upper(), line=dict(color=color),
            showlegend=True
        ), row=1, col=1)
        fig_error.add_trace(go.Scatter(
            x=alphas, y=cp_errs, mode="lines+markers",
            name=foil.upper(), line=dict(color=color),
            showlegend=False
        ), row=1, col=2)
    
    fig_error.update_layout(
        title="Error Analysis: RustFoil vs XFOIL",
        height=400
    )
    fig_error.update_yaxes(type="log", row=1, col=1)
    fig_error.update_yaxes(type="log", row=1, col=2)
    fig_error.update_xaxes(title_text="Angle of Attack (°)", row=1, col=1)
    fig_error.update_xaxes(title_text="Angle of Attack (°)", row=1, col=2)
    figures.append(("Error Analysis", fig_error))
    
    return figures

def save_html(figures, output_path):
    """Save figures to a single HTML file with tabs."""
    
    html_content = """
<!DOCTYPE html>
<html>
<head>
    <title>RustFoil vs XFOIL Comparison</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        body { font-family: Arial, sans-serif; margin: 20px; }
        h1 { text-align: center; }
        .tab-container { margin: 20px 0; }
        .tab-button { 
            padding: 10px 20px; margin: 2px; cursor: pointer;
            background: #f0f0f0; border: 1px solid #ccc; border-radius: 4px;
        }
        .tab-button:hover { background: #e0e0e0; }
        .tab-button.active { background: #007bff; color: white; }
        .tab-content { display: none; }
        .tab-content.active { display: block; }
        .stats { 
            font-family: monospace; background: #f5f5f5; 
            padding: 10px; margin: 10px 0; border-radius: 4px;
        }
    </style>
</head>
<body>
    <h1>RustFoil vs XFOIL Inviscid Solver Comparison</h1>
    <p style="text-align:center">3 NACA airfoils × 8 angles of attack (-4° to +10°)</p>
    
    <div class="tab-container">
"""
    
    for idx, (name, _) in enumerate(figures):
        active = "active" if idx == 0 else ""
        html_content += f'        <button class="tab-button {active}" onclick="showTab({idx})">{name}</button>\n'
    
    html_content += "    </div>\n\n"
    
    for idx, (name, fig) in enumerate(figures):
        active = "active" if idx == 0 else ""
        plot_html = fig.to_html(full_html=False, include_plotlyjs=False)
        html_content += f'    <div id="tab{idx}" class="tab-content {active}">\n'
        html_content += f'        {plot_html}\n'
        html_content += f'    </div>\n\n'
    
    html_content += """
    <script>
        function showTab(idx) {
            // Hide all tabs
            document.querySelectorAll('.tab-content').forEach(t => t.classList.remove('active'));
            document.querySelectorAll('.tab-button').forEach(b => b.classList.remove('active'));
            // Show selected tab
            document.getElementById('tab' + idx).classList.add('active');
            document.querySelectorAll('.tab-button')[idx].classList.add('active');
            // Resize plotly charts
            window.dispatchEvent(new Event('resize'));
        }
    </script>
</body>
</html>
"""
    
    with open(output_path, "w") as f:
        f.write(html_content)
    
    print(f"\nSaved to {output_path}")

if __name__ == "__main__":
    figures = create_full_comparison()
    output_path = DATA_DIR / "comparison_report.html"
    save_html(figures, output_path)
    print(f"\nOpen {output_path} in your browser to view interactive plots")
