#!/usr/bin/env python3
"""
Compare DIJ influence matrices between XFOIL and RustFoil.

The DIJ matrix represents dUe_i / d(delta_star * Ue)_j - the influence of
mass defect changes at station j on edge velocity at station i.

This script:
1. Loads XFOIL and RustFoil debug JSON files
2. Extracts FULLDIJ events containing the DIJ matrix
3. Compares element-by-element:
   - Computes max/mean absolute difference
   - Computes relative differences
   - Finds first element (i,j) where difference exceeds threshold
4. Generates heatmap visualization (optional with --plot)

Usage:
    python scripts/compare_dij.py xfoil_debug.json rustfoil_debug.json
    python scripts/compare_dij.py xfoil_debug.json rustfoil_debug.json --tolerance 0.01 --plot
"""

import json
import argparse
import numpy as np
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Any


def load_debug_events(json_path: Path) -> List[Dict[str, Any]]:
    """Load events from debug JSON file."""
    with open(json_path) as f:
        data = json.load(f)
    
    if isinstance(data, dict) and 'events' in data:
        return data['events']
    elif isinstance(data, list):
        return data
    else:
        raise ValueError(f"Unexpected JSON structure in {json_path}")


def extract_dij_matrix(events: List[Dict], source_name: str = "") -> Optional[Tuple[int, np.ndarray]]:
    """
    Extract DIJ matrix from FULLDIJ event.
    
    Returns:
        Tuple of (nsys, dij_matrix) or None if not found
    """
    # Look for FULLDIJ event (RustFoil format)
    for event in events:
        if event.get('subroutine') == 'FULLDIJ':
            nsys = event.get('nsys', 0)
            dij_flat = event.get('dij', [])
            if nsys > 0 and len(dij_flat) == nsys * nsys:
                dij = np.array(dij_flat).reshape((nsys, nsys))
                return nsys, dij
    
    # XFOIL doesn't log full DIJ in the same format, but we can 
    # reconstruct from QDCALC samples or the txt dump
    for event in events:
        if event.get('subroutine') == 'QDCALC':
            # QDCALC contains samples but not full matrix
            print(f"  [{source_name}] Found QDCALC event with samples, but not full DIJ")
            print(f"    n_total: {event.get('n_total', 'N/A')}")
            print(f"    diagonal_sample: {event.get('DIJ_diagonal_sample', [])[:5]}...")
            return None
    
    return None


def load_xfoil_dij_txt(txt_path: Path) -> Optional[Tuple[int, np.ndarray]]:
    """
    Load DIJ matrix from XFOIL text dump file.
    
    Format:
        N_PANELS= 160
        N_WAKE=  23
        # DIJ diagonal (first 20)
        ...
        # Full DIJ matrix (optional)
    """
    if not txt_path.exists():
        return None
    
    with open(txt_path) as f:
        lines = f.readlines()
    
    n_panels = 0
    n_wake = 0
    
    # Parse header
    for line in lines:
        if line.startswith('N_PANELS='):
            n_panels = int(line.split('=')[1].strip())
        elif line.startswith('N_WAKE='):
            n_wake = int(line.split('=')[1].strip())
    
    if n_panels == 0:
        return None
    
    n_total = n_panels + n_wake
    
    # Try to find and parse full matrix data
    # For now, we extract what diagonal and sample data we have
    dij = np.zeros((n_total, n_total))
    
    mode = None
    for line in lines:
        line = line.strip()
        if line.startswith('#'):
            if 'diagonal' in line.lower():
                mode = 'diagonal'
            elif 'row 1' in line.lower():
                mode = 'row1'
            elif 'row 80' in line.lower():
                mode = 'row80'
            elif 'stag' in line.lower():
                mode = 'stag'
            elif 'full' in line.lower():
                mode = 'full'
            continue
        
        parts = line.split()
        if not parts:
            continue
        
        try:
            if mode == 'diagonal' and len(parts) == 2:
                i = int(parts[0]) - 1  # 1-indexed to 0-indexed
                val = float(parts[1])
                if 0 <= i < n_total:
                    dij[i, i] = val
            elif mode == 'row1' and len(parts) == 2:
                j = int(parts[0]) - 1
                val = float(parts[1])
                if 0 <= j < n_total:
                    dij[0, j] = val
            elif mode == 'row80' and len(parts) == 2:
                j = int(parts[0]) - 1
                val = float(parts[1])
                if 0 <= j < n_total and 79 < n_total:
                    dij[79, j] = val
            elif mode == 'stag' and len(parts) == 3:
                i = int(parts[0]) - 1
                j = int(parts[1]) - 1
                val = float(parts[2])
                if 0 <= i < n_total and 0 <= j < n_total:
                    dij[i, j] = val
            elif mode == 'full' and len(parts) == 3:
                i = int(parts[0]) - 1
                j = int(parts[1]) - 1
                val = float(parts[2])
                if 0 <= i < n_total and 0 <= j < n_total:
                    dij[i, j] = val
        except (ValueError, IndexError):
            continue
    
    return n_total, dij


def compare_matrices(
    mat1: np.ndarray,
    mat2: np.ndarray,
    name1: str,
    name2: str,
    tolerance: float = 0.01
) -> Dict[str, Any]:
    """
    Compare two matrices element-by-element.
    
    Returns dict with comparison statistics.
    """
    # Ensure same shape (pad smaller matrix with zeros if needed)
    n1, m1 = mat1.shape
    n2, m2 = mat2.shape
    n = max(n1, n2)
    m = max(m1, m2)
    
    if n1 != n2 or m1 != m2:
        print(f"  Warning: Matrix size mismatch - {name1}: {mat1.shape}, {name2}: {mat2.shape}")
        print(f"  Comparing common region: ({min(n1, n2)}, {min(m1, m2)})")
        n = min(n1, n2)
        m = min(m1, m2)
        mat1 = mat1[:n, :m]
        mat2 = mat2[:n, :m]
    
    diff = mat1 - mat2
    abs_diff = np.abs(diff)
    
    # Compute statistics
    max_abs_diff = np.max(abs_diff)
    mean_abs_diff = np.mean(abs_diff)
    
    # Find first element exceeding tolerance
    first_exceed = None
    for i in range(n):
        for j in range(m):
            ref = max(abs(mat1[i, j]), abs(mat2[i, j]), 1e-10)
            rel_diff = abs_diff[i, j] / ref
            if rel_diff > tolerance:
                first_exceed = (i, j, mat1[i, j], mat2[i, j], abs_diff[i, j], rel_diff)
                break
        if first_exceed:
            break
    
    # Relative differences (avoid division by zero)
    with np.errstate(divide='ignore', invalid='ignore'):
        ref_max = np.maximum(np.abs(mat1), np.abs(mat2))
        ref_max = np.where(ref_max < 1e-20, 1.0, ref_max)
        rel_diff = abs_diff / ref_max
    
    max_rel_diff = np.max(rel_diff)
    
    # Diagonal analysis
    diag1 = np.diag(mat1)
    diag2 = np.diag(mat2)
    diag_diff = np.abs(diag1 - diag2)
    
    # Symmetry analysis
    sym_err1 = np.max(np.abs(mat1 - mat1.T))
    sym_err2 = np.max(np.abs(mat2 - mat2.T))
    
    return {
        'shape': (n, m),
        'max_abs_diff': max_abs_diff,
        'mean_abs_diff': mean_abs_diff,
        'max_rel_diff': max_rel_diff,
        'first_exceed': first_exceed,
        'diag_max_diff': np.max(diag_diff),
        'diag_mean_diff': np.mean(diag_diff),
        'symmetry_err_1': sym_err1,
        'symmetry_err_2': sym_err2,
        'diff_matrix': diff,
        'mat1': mat1,
        'mat2': mat2,
    }


def print_comparison_report(result: Dict[str, Any], name1: str, name2: str, tolerance: float):
    """Print detailed comparison report."""
    print("\n" + "=" * 70)
    print("DIJ MATRIX COMPARISON")
    print("=" * 70)
    
    print(f"\nMatrix size: {result['shape'][0]} x {result['shape'][1]}")
    print(f"Tolerance: {tolerance}")
    
    print("\n--- Global Statistics ---")
    print(f"  Max absolute difference:  {result['max_abs_diff']:.6e}")
    print(f"  Mean absolute difference: {result['mean_abs_diff']:.6e}")
    print(f"  Max relative difference:  {result['max_rel_diff']:.2%}")
    
    print("\n--- Diagonal Analysis ---")
    print(f"  Max diagonal difference:  {result['diag_max_diff']:.6e}")
    print(f"  Mean diagonal difference: {result['diag_mean_diff']:.6e}")
    
    print("\n--- Symmetry Check ---")
    print(f"  {name1} symmetry error: {result['symmetry_err_1']:.6e}")
    print(f"  {name2} symmetry error: {result['symmetry_err_2']:.6e}")
    
    if result['first_exceed']:
        i, j, v1, v2, abs_d, rel_d = result['first_exceed']
        print(f"\n--- First Exceedance (tol={tolerance}) ---")
        print(f"  Location: ({i}, {j})")
        print(f"  {name1}[{i},{j}] = {v1:.6e}")
        print(f"  {name2}[{i},{j}] = {v2:.6e}")
        print(f"  Absolute diff: {abs_d:.6e}")
        print(f"  Relative diff: {rel_d:.2%}")
    else:
        print(f"\n  All elements within tolerance ({tolerance})")
    
    # Sample diagonal values
    mat1 = result['mat1']
    mat2 = result['mat2']
    n = min(10, result['shape'][0])
    
    print(f"\n--- Diagonal Sample (first {n}) ---")
    print(f"{'i':>4} | {name1:>14} | {name2:>14} | {'Diff':>12} | {'Rel%':>8}")
    print("-" * 60)
    for i in range(n):
        v1 = mat1[i, i]
        v2 = mat2[i, i]
        diff = abs(v1 - v2)
        ref = max(abs(v1), abs(v2), 1e-10)
        rel = diff / ref * 100
        print(f"{i:4d} | {v1:14.6e} | {v2:14.6e} | {diff:12.4e} | {rel:7.2f}%")
    
    # Sample off-diagonal (row 1)
    print(f"\n--- Row 0 Sample (first {n} cols) ---")
    print(f"{'j':>4} | {name1:>14} | {name2:>14} | {'Diff':>12} | {'Rel%':>8}")
    print("-" * 60)
    for j in range(n):
        v1 = mat1[0, j]
        v2 = mat2[0, j]
        diff = abs(v1 - v2)
        ref = max(abs(v1), abs(v2), 1e-10)
        rel = diff / ref * 100
        print(f"{j:4d} | {v1:14.6e} | {v2:14.6e} | {diff:12.4e} | {rel:7.2f}%")
    
    # Stagnation region (around n/2)
    mid = result['shape'][0] // 2
    print(f"\n--- Stagnation Region (rows {mid-2}:{mid+2}, cols {mid-2}:{mid+2}) ---")
    for i in range(mid - 2, mid + 3):
        if 0 <= i < result['shape'][0]:
            row_str = f"Row {i:3d}: "
            for j in range(mid - 2, mid + 3):
                if 0 <= j < result['shape'][1]:
                    diff = abs(mat1[i, j] - mat2[i, j])
                    row_str += f"{diff:10.2e} "
            print(row_str)


def plot_comparison(result: Dict[str, Any], name1: str, name2: str, output_path: Optional[Path] = None):
    """Generate heatmap visualization of DIJ comparison."""
    try:
        import matplotlib.pyplot as plt
        from matplotlib.colors import LogNorm, SymLogNorm
    except ImportError:
        print("Warning: matplotlib not available, skipping plot")
        return
    
    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    
    mat1 = result['mat1']
    mat2 = result['mat2']
    diff = result['diff_matrix']
    
    # Determine color scale
    vmax1 = max(np.abs(mat1).max(), np.abs(mat2).max())
    vmin1 = -vmax1
    
    # Plot mat1
    ax = axes[0, 0]
    im = ax.imshow(mat1, aspect='auto', cmap='RdBu_r', vmin=vmin1, vmax=vmax1)
    ax.set_title(f'{name1} DIJ Matrix')
    ax.set_xlabel('Column (j)')
    ax.set_ylabel('Row (i)')
    plt.colorbar(im, ax=ax, label='DIJ value')
    
    # Plot mat2
    ax = axes[0, 1]
    im = ax.imshow(mat2, aspect='auto', cmap='RdBu_r', vmin=vmin1, vmax=vmax1)
    ax.set_title(f'{name2} DIJ Matrix')
    ax.set_xlabel('Column (j)')
    ax.set_ylabel('Row (i)')
    plt.colorbar(im, ax=ax, label='DIJ value')
    
    # Plot absolute difference (log scale)
    ax = axes[1, 0]
    abs_diff = np.abs(diff)
    abs_diff = np.where(abs_diff < 1e-20, 1e-20, abs_diff)  # Floor for log scale
    im = ax.imshow(abs_diff, aspect='auto', cmap='hot', norm=LogNorm(vmin=1e-10, vmax=abs_diff.max()))
    ax.set_title('Absolute Difference (log scale)')
    ax.set_xlabel('Column (j)')
    ax.set_ylabel('Row (i)')
    plt.colorbar(im, ax=ax, label='|DIJ1 - DIJ2|')
    
    # Plot relative difference
    ax = axes[1, 1]
    with np.errstate(divide='ignore', invalid='ignore'):
        ref = np.maximum(np.abs(mat1), np.abs(mat2))
        ref = np.where(ref < 1e-20, 1.0, ref)
        rel_diff = abs_diff / ref
    rel_diff = np.clip(rel_diff, 1e-10, 10)  # Clip for visualization
    im = ax.imshow(rel_diff, aspect='auto', cmap='viridis', norm=LogNorm(vmin=1e-6, vmax=1))
    ax.set_title('Relative Difference')
    ax.set_xlabel('Column (j)')
    ax.set_ylabel('Row (i)')
    plt.colorbar(im, ax=ax, label='Relative diff')
    
    plt.tight_layout()
    
    if output_path:
        plt.savefig(output_path, dpi=150, bbox_inches='tight')
        print(f"\nPlot saved to: {output_path}")
    else:
        plt.show()


def main():
    parser = argparse.ArgumentParser(
        description='Compare DIJ influence matrices between XFOIL and RustFoil',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/compare_dij.py xfoil_debug.json rustfoil_debug.json
    python scripts/compare_dij.py xfoil_debug.json rustfoil_debug.json --tolerance 0.01 --plot
    python scripts/compare_dij.py --xfoil-txt Xfoil-instrumented/bin/xfoil_dij_dump.txt rustfoil_debug.json
        """
    )
    parser.add_argument('xfoil_json', nargs='?', help='Path to XFOIL debug JSON file')
    parser.add_argument('rustfoil_json', help='Path to RustFoil debug JSON file')
    parser.add_argument('--xfoil-txt', help='Path to XFOIL DIJ text dump (alternative to JSON)')
    parser.add_argument('--tolerance', '-t', type=float, default=0.01,
                        help='Relative tolerance for element comparison (default: 0.01)')
    parser.add_argument('--plot', '-p', action='store_true',
                        help='Generate heatmap visualization')
    parser.add_argument('--output', '-o', help='Output path for plot (default: show interactively)')
    
    args = parser.parse_args()
    
    # Load RustFoil DIJ
    print(f"Loading RustFoil debug: {args.rustfoil_json}")
    rf_events = load_debug_events(Path(args.rustfoil_json))
    rf_dij = extract_dij_matrix(rf_events, "RustFoil")
    
    if rf_dij is None:
        print("ERROR: No FULLDIJ event found in RustFoil debug output")
        print("Ensure RustFoil is run with RUSTFOIL_DEBUG=<path> to capture DIJ matrix")
        return 1
    
    rf_nsys, rf_mat = rf_dij
    print(f"  RustFoil DIJ: {rf_nsys} x {rf_nsys}")
    
    # Load XFOIL DIJ
    xf_mat = None
    xf_nsys = 0
    
    if args.xfoil_txt:
        print(f"Loading XFOIL DIJ from text: {args.xfoil_txt}")
        xf_dij = load_xfoil_dij_txt(Path(args.xfoil_txt))
        if xf_dij:
            xf_nsys, xf_mat = xf_dij
            print(f"  XFOIL DIJ: {xf_nsys} x {xf_nsys} (partial from txt)")
    elif args.xfoil_json:
        print(f"Loading XFOIL debug: {args.xfoil_json}")
        xf_events = load_debug_events(Path(args.xfoil_json))
        xf_dij = extract_dij_matrix(xf_events, "XFOIL")
        if xf_dij:
            xf_nsys, xf_mat = xf_dij
            print(f"  XFOIL DIJ: {xf_nsys} x {xf_nsys}")
    
    if xf_mat is None:
        print("\nWarning: No XFOIL DIJ matrix available")
        print("Options:")
        print("  1. Use --xfoil-txt to load from xfoil_dij_dump.txt")
        print("  2. Instrument XFOIL to dump full DIJ matrix")
        print("\nShowing RustFoil DIJ statistics only:")
        print(f"  Shape: {rf_mat.shape}")
        print(f"  Diagonal range: [{np.diag(rf_mat).min():.4e}, {np.diag(rf_mat).max():.4e}]")
        print(f"  Symmetry error: {np.max(np.abs(rf_mat - rf_mat.T)):.4e}")
        return 0
    
    # Compare matrices
    result = compare_matrices(xf_mat, rf_mat, "XFOIL", "RustFoil", args.tolerance)
    print_comparison_report(result, "XFOIL", "RustFoil", args.tolerance)
    
    # Plot if requested
    if args.plot:
        output_path = Path(args.output) if args.output else None
        plot_comparison(result, "XFOIL", "RustFoil", output_path)
    
    # Return success if within tolerance
    if result['first_exceed'] is None:
        print("\n✓ DIJ matrices match within tolerance")
        return 0
    else:
        print("\n✗ DIJ matrices differ beyond tolerance")
        return 1


if __name__ == '__main__':
    exit(main())
