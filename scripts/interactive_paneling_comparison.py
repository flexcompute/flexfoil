#!/usr/bin/env python3
"""
Interactive comparison of XFOIL and RustFoil paneling using ACTUAL dumped coordinates.

This script loads real paneled coordinates from both XFOIL and RustFoil
to provide a true comparison.
"""

import numpy as np
import plotly.graph_objects as go
from plotly.subplots import make_subplots
from pathlib import Path
import webbrowser


def load_dat(filepath):
    """Load airfoil coordinates from .dat file."""
    coords = []
    with open(filepath) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    coords.append((float(parts[0]), float(parts[1])))
                except:
                    pass
    return np.array(coords)


def panel_lengths(coords):
    """Compute panel lengths."""
    diff = np.diff(coords, axis=0)
    return np.sqrt(diff[:, 0]**2 + diff[:, 1]**2)


def find_le_idx(coords):
    """Find leading edge index (minimum x)."""
    return np.argmin(coords[:, 0])


def create_comparison():
    workspace = Path(__file__).parent.parent
    
    foils = [
        ('NACA 0012', 'naca0012'),
        ('NACA 2412', 'naca2412'),
        ('NACA 4412', 'naca4412'),
    ]
    
    n_foils = len(foils)
    
    # Create figure
    fig = make_subplots(
        rows=4, cols=n_foils,
        subplot_titles=[f[0] for f in foils] + 
                       ['Panel Lengths' for _ in foils] +
                       ['Position Error' for _ in foils] +
                       ['LE Region Detail' for _ in foils],
        vertical_spacing=0.07,
        horizontal_spacing=0.06,
        row_heights=[0.35, 0.2, 0.2, 0.25],
        specs=[[{"type": "scatter"}]*n_foils,
               [{"type": "scatter"}]*n_foils,
               [{"type": "scatter"}]*n_foils,
               [{"type": "scatter"}]*n_foils],
    )
    
    colors = {
        'XFOIL': '#1f77b4',
        'RustFoil': '#ff7f0e',
        'error': '#d62728',
    }
    
    summary = []
    
    for col, (label, foil_id) in enumerate(foils, 1):
        # Load data
        xfoil_path = workspace / f'{foil_id}_paneled_xfoil.dat'
        rustfoil_path = workspace / f'{foil_id}_rustfoil_paneled.dat'
        
        if not xfoil_path.exists() or not rustfoil_path.exists():
            print(f"Skipping {foil_id}: files not found")
            continue
        
        xfoil = load_dat(xfoil_path)
        rustfoil = load_dat(rustfoil_path)
        
        xf_le = find_le_idx(xfoil)
        rf_le = find_le_idx(rustfoil)
        
        # Calculate errors
        pos_error = np.sqrt((xfoil[:, 0] - rustfoil[:, 0])**2 + 
                            (xfoil[:, 1] - rustfoil[:, 1])**2)
        rms_error = np.sqrt(np.mean(pos_error**2))
        max_error = pos_error.max()
        
        summary.append({
            'foil': label,
            'xfoil_le': xf_le,
            'rustfoil_le': rf_le,
            'le_diff': rf_le - xf_le,
            'rms': rms_error,
            'max': max_error,
        })
        
        # Row 1: Full geometry
        fig.add_trace(
            go.Scatter(
                x=xfoil[:, 0], y=xfoil[:, 1],
                mode='lines+markers',
                name='XFOIL',
                marker=dict(size=4, color=colors['XFOIL']),
                line=dict(width=1.5, color=colors['XFOIL']),
                legendgroup='xfoil',
                showlegend=(col == 1),
                hovertemplate='idx=%{customdata}<br>x=%{x:.6f}<br>y=%{y:.6f}<extra>XFOIL</extra>',
                customdata=list(range(len(xfoil))),
            ),
            row=1, col=col
        )
        
        fig.add_trace(
            go.Scatter(
                x=rustfoil[:, 0], y=rustfoil[:, 1],
                mode='lines+markers',
                name='RustFoil',
                marker=dict(size=4, color=colors['RustFoil'], symbol='x'),
                line=dict(width=1.5, dash='dash', color=colors['RustFoil']),
                legendgroup='rustfoil',
                showlegend=(col == 1),
                hovertemplate='idx=%{customdata}<br>x=%{x:.6f}<br>y=%{y:.6f}<extra>RustFoil</extra>',
                customdata=list(range(len(rustfoil))),
            ),
            row=1, col=col
        )
        
        # Mark LE points
        fig.add_trace(
            go.Scatter(
                x=[xfoil[xf_le, 0]], y=[xfoil[xf_le, 1]],
                mode='markers',
                name=f'XFOIL LE (idx={xf_le})',
                marker=dict(size=12, color=colors['XFOIL'], symbol='star'),
                showlegend=False,
                hovertemplate=f'XFOIL LE<br>idx={xf_le}<br>x=%{{x:.6f}}<br>y=%{{y:.6f}}<extra></extra>',
            ),
            row=1, col=col
        )
        
        fig.add_trace(
            go.Scatter(
                x=[rustfoil[rf_le, 0]], y=[rustfoil[rf_le, 1]],
                mode='markers',
                name=f'RustFoil LE (idx={rf_le})',
                marker=dict(size=12, color=colors['RustFoil'], symbol='star'),
                showlegend=False,
                hovertemplate=f'RustFoil LE<br>idx={rf_le}<br>x=%{{x:.6f}}<br>y=%{{y:.6f}}<extra></extra>',
            ),
            row=1, col=col
        )
        
        # Row 2: Panel lengths
        xf_len = panel_lengths(xfoil) * 1000  # millicord
        rf_len = panel_lengths(rustfoil) * 1000
        
        fig.add_trace(
            go.Scatter(
                x=list(range(len(xf_len))), y=xf_len,
                mode='lines',
                name='XFOIL',
                line=dict(width=2, color=colors['XFOIL']),
                legendgroup='xfoil',
                showlegend=False,
                hovertemplate='Panel %{x}<br>Length: %{y:.3f} mc<extra>XFOIL</extra>',
            ),
            row=2, col=col
        )
        
        fig.add_trace(
            go.Scatter(
                x=list(range(len(rf_len))), y=rf_len,
                mode='lines',
                name='RustFoil',
                line=dict(width=2, dash='dash', color=colors['RustFoil']),
                legendgroup='rustfoil',
                showlegend=False,
                hovertemplate='Panel %{x}<br>Length: %{y:.3f} mc<extra>RustFoil</extra>',
            ),
            row=2, col=col
        )
        
        # Row 3: Position error
        fig.add_trace(
            go.Scatter(
                x=list(range(len(pos_error))), y=pos_error,
                mode='lines',
                name='Error',
                line=dict(width=2, color=colors['error']),
                fill='tozeroy',
                fillcolor='rgba(214, 39, 40, 0.3)',
                legendgroup='error',
                showlegend=(col == 1),
                hovertemplate='Node %{x}<br>Error: %{y:.2e}<extra></extra>',
            ),
            row=3, col=col
        )
        
        # Add RMS and LE diff annotation
        fig.add_annotation(
            x=len(pos_error)//2, 
            y=max_error * 0.8,
            text=f'RMS: {rms_error:.2e}<br>LE diff: {rf_le - xf_le:+d}',
            showarrow=False,
            font=dict(size=10),
            bgcolor='rgba(255,255,255,0.8)',
            row=3, col=col
        )
        
        # Row 4: LE region zoom
        le_range = 10
        start = max(0, min(xf_le, rf_le) - le_range)
        end = min(len(xfoil), max(xf_le, rf_le) + le_range + 1)
        
        fig.add_trace(
            go.Scatter(
                x=xfoil[start:end, 0], y=xfoil[start:end, 1],
                mode='lines+markers',
                name='XFOIL',
                marker=dict(size=6, color=colors['XFOIL']),
                line=dict(width=2, color=colors['XFOIL']),
                legendgroup='xfoil',
                showlegend=False,
                hovertemplate='idx=%{customdata}<br>x=%{x:.8f}<br>y=%{y:.8f}<extra>XFOIL</extra>',
                customdata=list(range(start, end)),
            ),
            row=4, col=col
        )
        
        fig.add_trace(
            go.Scatter(
                x=rustfoil[start:end, 0], y=rustfoil[start:end, 1],
                mode='lines+markers',
                name='RustFoil',
                marker=dict(size=6, color=colors['RustFoil'], symbol='x'),
                line=dict(width=2, dash='dash', color=colors['RustFoil']),
                legendgroup='rustfoil',
                showlegend=False,
                hovertemplate='idx=%{customdata}<br>x=%{x:.8f}<br>y=%{y:.8f}<extra>RustFoil</extra>',
                customdata=list(range(start, end)),
            ),
            row=4, col=col
        )
        
        # Mark LE in zoom
        fig.add_trace(
            go.Scatter(
                x=[xfoil[xf_le, 0]], y=[xfoil[xf_le, 1]],
                mode='markers+text',
                text=[f'XF LE={xf_le}'],
                textposition='top center',
                marker=dict(size=10, color=colors['XFOIL'], symbol='star'),
                showlegend=False,
            ),
            row=4, col=col
        )
        
        fig.add_trace(
            go.Scatter(
                x=[rustfoil[rf_le, 0]], y=[rustfoil[rf_le, 1]],
                mode='markers+text',
                text=[f'RF LE={rf_le}'],
                textposition='bottom center',
                marker=dict(size=10, color=colors['RustFoil'], symbol='star'),
                showlegend=False,
            ),
            row=4, col=col
        )
    
    # Update layout
    fig.update_layout(
        title=dict(
            text='<b>XFOIL vs RustFoil Paneling: Actual Coordinates Comparison</b><br>' +
                 '<sup>Key finding: LE index differs for cambered airfoils (NACA 2412: +2, NACA 4412: +3)</sup>',
            x=0.5,
            font=dict(size=18)
        ),
        height=1100,
        width=450 * n_foils,
        showlegend=True,
        legend=dict(
            orientation='h',
            yanchor='bottom',
            y=1.01,
            xanchor='center',
            x=0.5
        ),
        hovermode='closest',
    )
    
    # Update axes
    for col in range(1, n_foils + 1):
        fig.update_xaxes(title_text='x/c', row=1, col=col)
        fig.update_yaxes(title_text='y/c', scaleanchor=f'x{col if col > 1 else ""}', row=1, col=col)
        
        fig.update_xaxes(title_text='Panel Index', row=2, col=col)
        fig.update_yaxes(title_text='Length (mc)', row=2, col=col)
        
        fig.update_xaxes(title_text='Node Index', row=3, col=col)
        fig.update_yaxes(title_text='Position Error', type='log', row=3, col=col)
        
        fig.update_xaxes(title_text='x/c (LE zoom)', row=4, col=col)
        fig.update_yaxes(title_text='y/c', scaleanchor=f'x{col + n_foils*3 if col > 1 else n_foils*3 + 1}', row=4, col=col)
    
    # Save
    output_path = workspace / 'scripts' / 'paneling_comparison_actual.html'
    fig.write_html(str(output_path), include_plotlyjs='cdn')
    
    # Print summary
    print("\n" + "="*70)
    print("PANELING COMPARISON SUMMARY (Actual XFOIL vs RustFoil coordinates)")
    print("="*70)
    print(f"{'Foil':<12} {'XFOIL LE':<10} {'RustFoil LE':<12} {'LE Diff':<10} {'RMS Error':<12} {'Max Error'}")
    print("-"*70)
    for s in summary:
        print(f"{s['foil']:<12} {s['xfoil_le']:<10} {s['rustfoil_le']:<12} {s['le_diff']:+d}{'':8} {s['rms']:<12.2e} {s['max']:.2e}")
    
    print("="*70)
    print(f"\nSaved interactive plot to: {output_path}")
    
    # Open in browser
    webbrowser.open(f'file://{output_path}')
    
    return summary


if __name__ == '__main__':
    create_comparison()
