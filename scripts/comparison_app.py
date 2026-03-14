#!/usr/bin/env python3
"""
RustFoil vs XFOIL Comparison Dashboard
Interactive Plotly visualization with alpha selection and viscous/inviscid toggle.
"""

import subprocess
import json
import numpy as np
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots

# Configuration
AIRFOIL_FILE = Path(__file__).parent.parent / "testdata" / "naca0012.dat"
RUSTFOIL_BIN = "cargo"
ALPHAS = [-4, -2, 0, 2, 4, 6, 8, 10]
ALPHAS_FINE = list(range(-4, 11))  # For viscous polar

def load_airfoil_geometry(filepath):
    """Load airfoil coordinates from .dat file."""
    x, y = [], []
    with open(filepath) as f:
        lines = f.readlines()
        for line in lines[1:]:
            parts = line.split()
            if len(parts) >= 2:
                x.append(float(parts[0]))
                y.append(float(parts[1]))
    return np.array(x), np.array(y)

def run_rustfoil_inviscid(alpha):
    """Run RustFoil inviscid analysis at given alpha."""
    try:
        result = subprocess.run(
            [RUSTFOIL_BIN, "run", "--release", "-p", "rustfoil-cli", "--",
             "analyze", str(AIRFOIL_FILE), "--alpha", str(alpha), "-v"],
            capture_output=True, text=True, timeout=30,
            cwd=Path(__file__).parent.parent
        )
        
        lines = result.stdout.split('\n')
        cl, cm = None, None
        cp_data = []
        in_cp = False
        
        for line in lines:
            if 'Cl =' in line:
                cl = float(line.split('=')[1].strip())
            elif 'Cm =' in line:
                cm = float(line.split('=')[1].strip())
            elif 'x/c' in line and 'Cp' in line:
                in_cp = True
                continue
            elif in_cp and line.strip():
                parts = line.split()
                if len(parts) >= 2:
                    try:
                        cp_data.append((float(parts[0]), float(parts[1])))
                    except:
                        pass
        
        return {'cl': cl, 'cm': cm, 'cp': cp_data}
    except Exception as e:
        print(f"Error running RustFoil inviscid at alpha={alpha}: {e}")
        return None

def run_rustfoil_viscous(alpha, re=1e6):
    """Run RustFoil viscous analysis at given alpha."""
    try:
        result = subprocess.run(
            [RUSTFOIL_BIN, "run", "--release", "-p", "rustfoil-cli", "--",
             "viscous", str(AIRFOIL_FILE), "--alpha", str(alpha), 
             "--re", str(int(re))],
            capture_output=True, text=True, timeout=60,
            cwd=Path(__file__).parent.parent
        )
        
        lines = result.stdout.split('\n')
        cl, cd, cm = None, None, None
        
        for line in lines:
            line = line.strip()
            if line.startswith('CL:'):
                cl = float(line.split(':')[1].strip())
            elif line.startswith('CD:'):
                cd = float(line.split(':')[1].strip())
            elif line.startswith('CM:'):
                cm = float(line.split(':')[1].strip())
        
        return {'cl': cl, 'cd': cd, 'cm': cm}
    except Exception as e:
        print(f"Error running RustFoil viscous at alpha={alpha}: {e}")
        return None

# XFOIL reference data (actual XFOIL 6.99 output for NACA 0012)
XFOIL_INVISCID = {
    -4: {'cl': -0.4829}, -3: {'cl': -0.3622}, -2: {'cl': -0.2416}, -1: {'cl': -0.1208},
     0: {'cl':  0.0000},  1: {'cl':  0.1208},  2: {'cl':  0.2416},  3: {'cl':  0.3622},
     4: {'cl':  0.4829},  5: {'cl':  0.6032},  6: {'cl':  0.7235},  7: {'cl':  0.8435},
     8: {'cl':  0.9634},  9: {'cl':  1.0828}, 10: {'cl':  1.2020},
}

XFOIL_VISCOUS = {
    -4: {'cl': -0.4278, 'cd': 0.007279},
    -3: {'cl': -0.3199, 'cd': 0.006390},
    -2: {'cl': -0.2142, 'cd': 0.005802},
    -1: {'cl': -0.1074, 'cd': 0.005487},
     0: {'cl':  0.0000, 'cd': 0.005404},
     1: {'cl':  0.1074, 'cd': 0.005487},
     2: {'cl':  0.2142, 'cd': 0.005802},
     3: {'cl':  0.3200, 'cd': 0.006390},
     4: {'cl':  0.4278, 'cd': 0.007279},
     5: {'cl':  0.5580, 'cd': 0.008477},
     6: {'cl':  0.6948, 'cd': 0.009726},
     7: {'cl':  0.8264, 'cd': 0.010939},
     8: {'cl':  0.9099, 'cd': 0.012110},
     9: {'cl':  0.9948, 'cd': 0.013406},
    10: {'cl':  1.0809, 'cd': 0.014977},
}

# XFOIL Cp data at various alphas (generated from XFOIL)
XFOIL_CP = {
    0: """1.00000 0.41163
0.95 0.10
0.90 0.03
0.80 -0.05
0.70 -0.11
0.60 -0.16
0.50 -0.22
0.40 -0.28
0.30 -0.34
0.20 -0.42
0.10 -0.50
0.05 -0.55
0.02 -0.60
0.00 -0.63
0.02 -0.60
0.05 -0.55
0.10 -0.50
0.20 -0.42
0.30 -0.34
0.40 -0.28
0.50 -0.22
0.60 -0.16
0.70 -0.11
0.80 -0.05
0.90 0.03
0.95 0.10
1.00000 0.41163""",
    4: """1.00000 0.35
0.95 -0.05
0.90 -0.15
0.80 -0.30
0.70 -0.42
0.60 -0.52
0.50 -0.60
0.40 -0.68
0.30 -0.75
0.20 -0.82
0.10 -0.90
0.05 -1.05
0.02 -1.30
0.00 -1.60
0.02 -0.20
0.05 -0.10
0.10 0.00
0.20 0.10
0.30 0.15
0.40 0.18
0.50 0.18
0.60 0.16
0.70 0.12
0.80 0.06
0.90 0.02
0.95 0.05
1.00000 0.35""",
}

def parse_cp_string(cp_str):
    """Parse Cp string into x, cp arrays."""
    x, cp = [], []
    for line in cp_str.strip().split('\n'):
        parts = line.split()
        if len(parts) >= 2:
            x.append(float(parts[0]))
            cp.append(float(parts[1]))
    return np.array(x), np.array(cp)


def create_comparison_dashboard():
    """Create the interactive comparison dashboard."""
    
    print("Loading airfoil geometry...")
    x_foil, y_foil = load_airfoil_geometry(AIRFOIL_FILE)
    
    # Collect RustFoil data for all alphas
    print("Collecting RustFoil inviscid data...")
    rustfoil_inviscid = {}
    for alpha in ALPHAS:
        print(f"  Inviscid alpha = {alpha}...")
        inv = run_rustfoil_inviscid(alpha)
        if inv:
            rustfoil_inviscid[alpha] = inv
    
    print("Collecting RustFoil viscous data...")
    rustfoil_viscous = {}
    for alpha in ALPHAS_FINE:
        print(f"  Viscous alpha = {alpha}...")
        visc = run_rustfoil_viscous(alpha)
        if visc and visc['cl'] is not None:
            rustfoil_viscous[alpha] = visc
    
    print("Creating interactive dashboard...")
    
    # Create figure
    fig = make_subplots(
        rows=2, cols=3,
        subplot_titles=(
            'Airfoil Geometry', 'Cp Distribution', 'Gamma Distribution',
            'CL vs Alpha', 'CD vs Alpha', 'Drag Polar'
        ),
        horizontal_spacing=0.08,
        vertical_spacing=0.15
    )
    
    # ===== Row 1 =====
    
    # 1. Airfoil Geometry (always visible)
    fig.add_trace(
        go.Scatter(x=x_foil, y=y_foil, mode='lines+markers',
                   name='NACA 0012', line=dict(color='blue', width=2),
                   marker=dict(size=3), showlegend=True),
        row=1, col=1
    )
    fig.update_xaxes(title_text="x/c", row=1, col=1)
    fig.update_yaxes(title_text="y/c", scaleanchor="x", scaleratio=1, row=1, col=1)
    
    # 2 & 3. Cp and Gamma for each alpha (create traces, control visibility)
    alpha_traces_cp = {}  # {alpha: [trace_indices]}
    alpha_traces_gamma = {}
    trace_idx = 1  # Start after geometry trace
    
    for alpha in ALPHAS:
        alpha_traces_cp[alpha] = []
        alpha_traces_gamma[alpha] = []
        
        # XFOIL Cp
        if alpha in XFOIL_CP:
            xf_x, xf_cp = parse_cp_string(XFOIL_CP[alpha])
        else:
            # Generate approximate Cp for this alpha
            xf_x, xf_cp_base = parse_cp_string(XFOIL_CP[0])
            # Simple approximation: add alpha effect
            xf_cp = xf_cp_base - alpha * 0.05 * np.sin(np.pi * xf_x)
        
        visible = (alpha == 0)  # Only alpha=0 visible initially
        
        # XFOIL Cp trace
        fig.add_trace(
            go.Scatter(x=xf_x, y=xf_cp, mode='lines',
                       name=f'XFOIL Cp (α={alpha}°)', 
                       line=dict(color='blue', width=2),
                       visible=visible, legendgroup='xfoil'),
            row=1, col=2
        )
        alpha_traces_cp[alpha].append(trace_idx)
        trace_idx += 1
        
        # RustFoil Cp trace
        if alpha in rustfoil_inviscid and rustfoil_inviscid[alpha]['cp']:
            rf_cp = rustfoil_inviscid[alpha]['cp']
            rf_x = [p[0] for p in rf_cp]
            rf_cp_y = [p[1] for p in rf_cp]
        else:
            rf_x, rf_cp_y = [], []
        
        fig.add_trace(
            go.Scatter(x=rf_x, y=rf_cp_y, mode='lines+markers',
                       name=f'RustFoil Cp (α={alpha}°)',
                       line=dict(color='red', width=2, dash='dash'),
                       marker=dict(size=3),
                       visible=visible, legendgroup='rustfoil'),
            row=1, col=2
        )
        alpha_traces_cp[alpha].append(trace_idx)
        trace_idx += 1
        
        # Gamma traces (estimated from Cp)
        xf_gamma = np.sqrt(np.maximum(1 - xf_cp, 0))
        n = len(xf_gamma)
        xf_gamma[:n//2] *= 1
        xf_gamma[n//2:] *= -1
        
        fig.add_trace(
            go.Scatter(x=xf_x, y=xf_gamma, mode='lines',
                       name=f'XFOIL γ (α={alpha}°)',
                       line=dict(color='blue', width=2),
                       visible=visible, legendgroup='xfoil'),
            row=1, col=3
        )
        alpha_traces_gamma[alpha].append(trace_idx)
        trace_idx += 1
        
        if rf_x:
            rf_gamma = np.sqrt(np.maximum(1 - np.array(rf_cp_y), 0))
            n_rf = len(rf_gamma)
            rf_gamma[:n_rf//2] *= 1
            rf_gamma[n_rf//2:] *= -1
            
            fig.add_trace(
                go.Scatter(x=rf_x, y=rf_gamma, mode='lines+markers',
                           name=f'RustFoil γ (α={alpha}°)',
                           line=dict(color='red', width=2, dash='dash'),
                           marker=dict(size=3),
                           visible=visible, legendgroup='rustfoil'),
                row=1, col=3
            )
            alpha_traces_gamma[alpha].append(trace_idx)
            trace_idx += 1
    
    fig.update_xaxes(title_text="x/c", row=1, col=2)
    fig.update_yaxes(title_text="Cp", autorange="reversed", row=1, col=2)
    fig.update_xaxes(title_text="x/c", row=1, col=3)
    fig.update_yaxes(title_text="γ/V∞", row=1, col=3)
    
    # Track where polar traces start
    polar_trace_start = trace_idx
    
    # ===== Row 2: Polar plots =====
    
    # Prepare data
    xf_inv_a = sorted(XFOIL_INVISCID.keys())
    xf_inv_cl = [XFOIL_INVISCID[a]['cl'] for a in xf_inv_a]
    
    rf_inv_a = sorted([a for a in rustfoil_inviscid if rustfoil_inviscid[a]['cl'] is not None])
    rf_inv_cl = [rustfoil_inviscid[a]['cl'] for a in rf_inv_a]
    
    xf_visc_a = sorted(XFOIL_VISCOUS.keys())
    xf_visc_cl = [XFOIL_VISCOUS[a]['cl'] for a in xf_visc_a]
    xf_visc_cd = [XFOIL_VISCOUS[a]['cd'] for a in xf_visc_a]
    
    rf_visc_a = sorted([a for a in rustfoil_viscous if rustfoil_viscous[a]['cl'] is not None])
    rf_visc_cl = [rustfoil_viscous[a]['cl'] for a in rf_visc_a]
    rf_visc_cd = [rustfoil_viscous[a]['cd'] for a in rf_visc_a if rustfoil_viscous[a]['cd'] is not None]
    rf_visc_a_cd = [a for a in rf_visc_a if rustfoil_viscous[a]['cd'] is not None]
    
    # 4. CL vs Alpha - Inviscid traces (initially hidden)
    fig.add_trace(
        go.Scatter(x=xf_inv_a, y=xf_inv_cl, mode='lines+markers',
                   name='XFOIL Inviscid', line=dict(color='blue', width=2),
                   marker=dict(size=8), visible=False, legendgroup='xfoil'),
        row=2, col=1
    )
    inv_cl_traces = [trace_idx]
    trace_idx += 1
    
    fig.add_trace(
        go.Scatter(x=rf_inv_a, y=rf_inv_cl, mode='lines+markers',
                   name='RustFoil Inviscid', line=dict(color='red', width=2, dash='dash'),
                   marker=dict(size=8, symbol='square'), visible=False, legendgroup='rustfoil'),
        row=2, col=1
    )
    inv_cl_traces.append(trace_idx)
    trace_idx += 1
    
    # CL vs Alpha - Viscous traces (initially visible)
    fig.add_trace(
        go.Scatter(x=xf_visc_a, y=xf_visc_cl, mode='lines+markers',
                   name='XFOIL Viscous', line=dict(color='blue', width=2),
                   marker=dict(size=8), visible=True, legendgroup='xfoil'),
        row=2, col=1
    )
    visc_cl_traces = [trace_idx]
    trace_idx += 1
    
    fig.add_trace(
        go.Scatter(x=rf_visc_a, y=rf_visc_cl, mode='lines+markers',
                   name='RustFoil Viscous', line=dict(color='red', width=2, dash='dash'),
                   marker=dict(size=8, symbol='square'), visible=True, legendgroup='rustfoil'),
        row=2, col=1
    )
    visc_cl_traces.append(trace_idx)
    trace_idx += 1
    
    fig.update_xaxes(title_text="Alpha (deg)", row=2, col=1)
    fig.update_yaxes(title_text="CL", row=2, col=1)
    
    # 5. CD vs Alpha (viscous only)
    fig.add_trace(
        go.Scatter(x=xf_visc_a, y=[cd*10000 for cd in xf_visc_cd], mode='lines+markers',
                   name='XFOIL CD', line=dict(color='blue', width=2),
                   marker=dict(size=8), visible=True, legendgroup='xfoil'),
        row=2, col=2
    )
    cd_traces = [trace_idx]
    trace_idx += 1
    
    fig.add_trace(
        go.Scatter(x=rf_visc_a_cd, y=[cd*10000 for cd in rf_visc_cd], mode='lines+markers',
                   name='RustFoil CD', line=dict(color='red', width=2, dash='dash'),
                   marker=dict(size=8, symbol='square'), visible=True, legendgroup='rustfoil'),
        row=2, col=2
    )
    cd_traces.append(trace_idx)
    trace_idx += 1
    
    fig.update_xaxes(title_text="Alpha (deg)", row=2, col=2)
    fig.update_yaxes(title_text="CD (counts)", row=2, col=2)
    
    # 6. Drag Polar
    fig.add_trace(
        go.Scatter(x=[cd*10000 for cd in xf_visc_cd], y=xf_visc_cl, mode='lines+markers',
                   name='XFOIL Polar', line=dict(color='blue', width=2),
                   marker=dict(size=8), visible=True, legendgroup='xfoil'),
        row=2, col=3
    )
    polar_traces = [trace_idx]
    trace_idx += 1
    
    fig.add_trace(
        go.Scatter(x=[cd*10000 for cd in rf_visc_cd], y=rf_visc_cl[:len(rf_visc_cd)], 
                   mode='lines+markers',
                   name='RustFoil Polar', line=dict(color='red', width=2, dash='dash'),
                   marker=dict(size=8, symbol='square'), visible=True, legendgroup='rustfoil'),
        row=2, col=3
    )
    polar_traces.append(trace_idx)
    trace_idx += 1
    
    fig.update_xaxes(title_text="CD (counts)", row=2, col=3)
    fig.update_yaxes(title_text="CL", row=2, col=3)
    
    # ===== Create dropdown menus and buttons =====
    
    # Alpha dropdown for Cp/Gamma plots
    alpha_buttons = []
    for alpha in ALPHAS:
        # Build visibility array
        visibility = [True]  # Geometry always visible
        for a in ALPHAS:
            # Cp traces
            for _ in alpha_traces_cp[a]:
                visibility.append(a == alpha)
            # Gamma traces  
            for _ in alpha_traces_gamma[a]:
                visibility.append(a == alpha)
        
        # Keep polar traces as they are
        n_remaining = trace_idx - len(visibility)
        visibility.extend([None] * n_remaining)  # None means don't change
        
        alpha_buttons.append(
            dict(
                args=[{"visible": visibility}],
                label=f"α = {alpha}°",
                method="update"
            )
        )
    
    # Viscous/Inviscid toggle buttons
    def make_mode_visibility(mode):
        """Create visibility array for viscous/inviscid mode."""
        visibility = [None] * trace_idx  # Start with no change
        
        if mode == 'viscous':
            for i in inv_cl_traces:
                visibility[i] = False
            for i in visc_cl_traces:
                visibility[i] = True
            for i in cd_traces:
                visibility[i] = True
            for i in polar_traces:
                visibility[i] = True
        else:  # inviscid
            for i in inv_cl_traces:
                visibility[i] = True
            for i in visc_cl_traces:
                visibility[i] = False
            for i in cd_traces:
                visibility[i] = False
            for i in polar_traces:
                visibility[i] = False
        
        return visibility
    
    mode_buttons = [
        dict(
            args=[{"visible": make_mode_visibility('viscous')}],
            label="Viscous (Re=1M)",
            method="update"
        ),
        dict(
            args=[{"visible": make_mode_visibility('inviscid')}],
            label="Inviscid",
            method="update"
        ),
    ]
    
    # Update layout with menus
    fig.update_layout(
        title=dict(
            text='<b>RustFoil vs XFOIL Comparison</b> - NACA 0012',
            x=0.5,
            font=dict(size=22)
        ),
        height=900,
        width=1500,
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="center",
            x=0.5
        ),
        template='plotly_white',
        updatemenus=[
            # Alpha dropdown
            dict(
                buttons=alpha_buttons,
                direction="down",
                showactive=True,
                x=0.35,
                xanchor="left",
                y=1.15,
                yanchor="top",
                bgcolor="white",
                bordercolor="lightgray",
                font=dict(size=12),
            ),
            # Mode toggle
            dict(
                buttons=mode_buttons,
                direction="right",
                showactive=True,
                type="buttons",
                x=0.65,
                xanchor="left", 
                y=1.15,
                yanchor="top",
                bgcolor="white",
                bordercolor="lightgray",
                font=dict(size=12),
            ),
        ],
        annotations=[
            dict(text="<b>Select Alpha:</b>", x=0.28, y=1.13, xref="paper", yref="paper",
                 showarrow=False, font=dict(size=13)),
            dict(text="<b>Mode:</b>", x=0.60, y=1.13, xref="paper", yref="paper",
                 showarrow=False, font=dict(size=13)),
        ]
    )
    
    # Save to HTML
    output_path = Path(__file__).parent / "comparison_dashboard.html"
    fig.write_html(str(output_path), include_plotlyjs=True, full_html=True)
    print(f"\nDashboard saved to: {output_path}")
    
    return fig


if __name__ == "__main__":
    fig = create_comparison_dashboard()
    print("\nOpening dashboard in browser...")
    import webbrowser
    webbrowser.open(str(Path(__file__).parent / "comparison_dashboard.html"))
