#!/usr/bin/env python3
"""
Visualization Script for Newton Iteration Comparison

Creates plots showing:
- Residual evolution (XFOIL vs RustFoil)
- BL variable profiles at divergence iteration
- Matrix element differences heatmap
- Station-by-station error progression

Usage:
    python visualize_divergence.py --xfoil-trace xfoil_alpha_2.json \
                                   --rustfoil-trace rustfoil_alpha_2.json \
                                   --output divergence_plots/
"""

import argparse
import sys
from pathlib import Path
from typing import Dict, List, Optional, Any
import json

# Try to import matplotlib
try:
    import matplotlib.pyplot as plt
    import matplotlib.colors as mcolors
    import numpy as np
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False
    print("Warning: matplotlib not available. Install with: pip install matplotlib")

# Add parent directory to path for imports
sys.path.insert(0, str(Path(__file__).parent))

from utils import (
    load_trace, extract_events, calculate_rel_error,
    get_iteration_events, summarize_trace
)


def plot_residual_evolution(
    xfoil_data: Dict, rustfoil_data: Dict, 
    output_dir: Path, alpha: float
):
    """Plot residual (RMSBL) evolution across iterations."""
    if not HAS_MATPLOTLIB:
        return
    
    # Extract residuals from NEWTON_ITER events
    xf_residuals = []
    rf_residuals = []
    
    xf_events = xfoil_data.get('events', [])
    rf_events = rustfoil_data.get('events', [])
    
    # XFOIL: look for VISCAL_RESULT or NEWTON_CONVERGE
    for e in xf_events:
        if e.get('subroutine') == 'VISCAL_RESULT':
            xf_residuals.append((e.get('iteration', 0), e.get('rms_residual', 0)))
        elif e.get('subroutine') == 'NEWTON_CONVERGE':
            xf_residuals.append((e.get('iteration', 0), e.get('rmsbl', 0)))
    
    # RustFoil: look for NEWTON_ITER
    for e in rf_events:
        if e.get('subroutine') == 'NEWTON_ITER':
            rf_residuals.append((e.get('iteration', 0), e.get('residual', 0)))
    
    if not xf_residuals and not rf_residuals:
        print("  No residual data found for plotting")
        return
    
    fig, ax = plt.subplots(figsize=(10, 6))
    
    if xf_residuals:
        iters, resids = zip(*sorted(xf_residuals))
        ax.semilogy(iters, resids, 'b-o', label='XFOIL', markersize=4)
    
    if rf_residuals:
        iters, resids = zip(*sorted(rf_residuals))
        ax.semilogy(iters, resids, 'r-s', label='RustFoil', markersize=4)
    
    ax.set_xlabel('Newton Iteration')
    ax.set_ylabel('RMSBL (log scale)')
    ax.set_title(f'Residual Evolution - Alpha = {alpha}°')
    ax.legend()
    ax.grid(True, alpha=0.3)
    
    output_path = output_dir / f'residual_evolution_a{int(alpha)}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_bl_profiles(
    xfoil_data: Dict, rustfoil_data: Dict,
    output_dir: Path, alpha: float, iteration: int = 0
):
    """Plot BL variable profiles comparing XFOIL and RustFoil."""
    if not HAS_MATPLOTLIB:
        return
    
    # Get FULL_ITER events for the specified iteration
    xf_full = extract_events(xfoil_data, 'FULL_ITER', iteration)
    rf_full = extract_events(rustfoil_data, 'FULL_ITER', iteration)
    
    if not xf_full or not rf_full:
        # Try BL_STATE_SUMMARY
        xf_full = extract_events(xfoil_data, 'BL_STATE', iteration)
        rf_full = extract_events(rustfoil_data, 'BL_STATE_SUMMARY', iteration)
        
        if not xf_full or not rf_full:
            print(f"  No BL data found for iteration {iteration}")
            return
    
    xf = xf_full[0] if isinstance(xf_full, list) else xf_full
    rf = rf_full[0] if isinstance(rf_full, list) else rf_full
    
    # Extract upper surface data
    xf_upper = xf.get('bl_upper', [])
    rf_upper = rf.get('bl_upper', rf.get('upper_stations', []))
    
    if not xf_upper or not rf_upper:
        print(f"  No BL upper surface data for iteration {iteration}")
        return
    
    # Create figure with subplots for different variables
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))
    
    variables = [
        ('theta', 'Momentum Thickness θ'),
        ('dstar', 'Displacement Thickness δ*'),
        ('ue', 'Edge Velocity Ue'),
        ('hk', 'Shape Factor Hk'),
    ]
    
    for ax, (var, title) in zip(axes.flat, variables):
        xf_x = [s.get('x', 0) for s in xf_upper]
        xf_y = [s.get(var, 0) for s in xf_upper]
        rf_x = [s.get('x', 0) for s in rf_upper]
        rf_y = [s.get(var, 0) for s in rf_upper]
        
        ax.plot(xf_x, xf_y, 'b-o', label='XFOIL', markersize=3)
        ax.plot(rf_x, rf_y, 'r-s', label='RustFoil', markersize=3)
        
        ax.set_xlabel('x/c')
        ax.set_ylabel(var)
        ax.set_title(title)
        ax.legend(fontsize=8)
        ax.grid(True, alpha=0.3)
    
    fig.suptitle(f'BL Profiles - Alpha = {alpha}°, Iteration {iteration}')
    plt.tight_layout()
    
    output_path = output_dir / f'bl_profiles_a{int(alpha)}_iter{iteration}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_error_heatmap(
    xfoil_data: Dict, rustfoil_data: Dict,
    output_dir: Path, alpha: float
):
    """Plot heatmap of relative errors across iterations and stations."""
    if not HAS_MATPLOTLIB:
        return
    
    # Collect error data: rows = iterations, cols = stations
    max_iter = 10
    max_station = 50
    
    errors = np.zeros((max_iter, max_station))
    
    for iteration in range(max_iter):
        xf_full = extract_events(xfoil_data, 'FULL_ITER', iteration)
        rf_full = extract_events(rustfoil_data, 'FULL_ITER', iteration)
        
        if not xf_full or not rf_full:
            continue
        
        xf = xf_full[0]
        rf = rf_full[0]
        
        xf_upper = xf.get('bl_upper', [])
        rf_upper = rf.get('bl_upper', [])
        
        for xf_st, rf_st in zip(xf_upper, rf_upper):
            ibl = xf_st.get('ibl', 0)
            if ibl >= max_station:
                continue
            
            # Use theta as the comparison variable
            xf_theta = xf_st.get('theta', 0)
            rf_theta = rf_st.get('theta', 0)
            
            if xf_theta != 0:
                rel_err = abs(rf_theta - xf_theta) / abs(xf_theta)
                errors[iteration, ibl] = rel_err
    
    if errors.sum() == 0:
        print("  No error data collected for heatmap")
        return
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    # Log scale for errors, with clipping
    errors_log = np.log10(errors + 1e-10)
    errors_log = np.clip(errors_log, -6, 0)  # Clip to 1e-6 to 1
    
    im = ax.imshow(errors_log, aspect='auto', cmap='RdYlGn_r',
                   vmin=-4, vmax=0)
    
    ax.set_xlabel('Station Index')
    ax.set_ylabel('Newton Iteration')
    ax.set_title(f'Relative Error in θ (log10) - Alpha = {alpha}°')
    
    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label('log10(relative error)')
    
    output_path = output_dir / f'error_heatmap_a{int(alpha)}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def plot_station_errors(
    xfoil_data: Dict, rustfoil_data: Dict,
    output_dir: Path, alpha: float, iteration: int = 1
):
    """Plot station-by-station relative errors for multiple variables."""
    if not HAS_MATPLOTLIB:
        return
    
    xf_full = extract_events(xfoil_data, 'FULL_ITER', iteration)
    rf_full = extract_events(rustfoil_data, 'FULL_ITER', iteration)
    
    if not xf_full or not rf_full:
        print(f"  No data for iteration {iteration}")
        return
    
    xf = xf_full[0]
    rf = rf_full[0]
    
    xf_upper = xf.get('bl_upper', [])
    rf_upper = rf.get('bl_upper', [])
    
    if not xf_upper or not rf_upper:
        return
    
    # Collect errors
    variables = ['theta', 'dstar', 'ue', 'mass']
    errors = {var: [] for var in variables}
    stations = []
    
    for xf_st, rf_st in zip(xf_upper, rf_upper):
        ibl = xf_st.get('ibl', 0)
        stations.append(ibl)
        
        for var in variables:
            xf_val = xf_st.get(var, 0)
            rf_val = rf_st.get(var, 0)
            
            if abs(xf_val) > 1e-10:
                rel_err = (rf_val - xf_val) / xf_val * 100
            else:
                rel_err = 0
            
            errors[var].append(rel_err)
    
    fig, ax = plt.subplots(figsize=(12, 6))
    
    colors = ['blue', 'red', 'green', 'orange']
    for var, color in zip(variables, colors):
        ax.plot(stations, errors[var], '-o', color=color, label=var, markersize=3)
    
    ax.axhline(y=0, color='k', linestyle='--', alpha=0.3)
    ax.axhline(y=1, color='r', linestyle=':', alpha=0.3, label='1% threshold')
    ax.axhline(y=-1, color='r', linestyle=':', alpha=0.3)
    
    ax.set_xlabel('Station Index')
    ax.set_ylabel('Relative Error (%)')
    ax.set_title(f'Station-wise Errors - Alpha = {alpha}°, Iteration {iteration}')
    ax.legend(fontsize=8)
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-10, 10)
    
    output_path = output_dir / f'station_errors_a{int(alpha)}_iter{iteration}.png'
    plt.savefig(output_path, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Saved: {output_path}")


def main():
    parser = argparse.ArgumentParser(
        description="Visualize Newton iteration divergence"
    )
    parser.add_argument(
        "--xfoil-trace", type=Path, required=True,
        help="Path to XFOIL trace JSON"
    )
    parser.add_argument(
        "--rustfoil-trace", type=Path, required=True,
        help="Path to RustFoil trace JSON"
    )
    parser.add_argument(
        "--output", type=Path, default=Path("divergence_plots"),
        help="Output directory for plots"
    )
    parser.add_argument(
        "--alpha", type=float, default=2.0,
        help="Alpha value for labeling"
    )
    parser.add_argument(
        "--iteration", type=int, default=1,
        help="Iteration to analyze for detailed plots"
    )
    
    args = parser.parse_args()
    
    if not HAS_MATPLOTLIB:
        print("ERROR: matplotlib is required for visualization")
        print("Install with: pip install matplotlib numpy")
        return 1
    
    print(f"Loading traces...")
    
    try:
        xfoil_data = load_trace(args.xfoil_trace)
        rustfoil_data = load_trace(args.rustfoil_trace)
    except Exception as e:
        print(f"Error loading traces: {e}")
        return 1
    
    # Create output directory
    args.output.mkdir(parents=True, exist_ok=True)
    
    print(f"\nGenerating plots for alpha = {args.alpha}°...")
    
    # Generate all plots
    plot_residual_evolution(xfoil_data, rustfoil_data, args.output, args.alpha)
    plot_bl_profiles(xfoil_data, rustfoil_data, args.output, args.alpha, args.iteration)
    plot_error_heatmap(xfoil_data, rustfoil_data, args.output, args.alpha)
    plot_station_errors(xfoil_data, rustfoil_data, args.output, args.alpha, args.iteration)
    
    print(f"\nPlots saved to: {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
