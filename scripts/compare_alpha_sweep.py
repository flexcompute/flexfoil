#!/usr/bin/env python3
"""
Compare RustFoil and XFOIL alpha sweeps and identify discrepancies.
"""
import json
from pathlib import Path
import numpy as np

def load_data():
    """Load both RustFoil and XFOIL data."""
    rf_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/rustfoil_alpha_sweep.json")
    xf_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/xfoil_reference.json")
    
    with open(rf_file) as f:
        rustfoil_data = json.load(f)
    
    with open(xf_file) as f:
        xfoil_data = json.load(f)
    
    return rustfoil_data, xfoil_data

def match_by_alpha(rf_results, xf_results):
    """Match RustFoil and XFOIL results by alpha."""
    rf_dict = {r['alpha']: r for r in rf_results if r['success']}
    xf_dict = {r['alpha']: r for r in xf_results}
    
    common_alphas = sorted(set(rf_dict.keys()) & set(xf_dict.keys()))
    
    comparisons = []
    for alpha in common_alphas:
        rf = rf_dict[alpha]
        xf = xf_dict[alpha]
        
        cl_error = None
        cd_error = None
        cl_pct_error = None
        cd_pct_error = None
        
        if rf['cl'] is not None and xf['cl'] is not None and xf['cl'] != 0:
            cl_error = rf['cl'] - xf['cl']
            cl_pct_error = 100 * cl_error / xf['cl']
        
        if rf['cd'] is not None and xf['cd'] is not None and xf['cd'] != 0:
            cd_error = rf['cd'] - xf['cd']
            cd_pct_error = 100 * cd_error / xf['cd']
        
        comparisons.append({
            'alpha': alpha,
            'rustfoil_cl': rf['cl'],
            'xfoil_cl': xf['cl'],
            'cl_error': cl_error,
            'cl_pct_error': cl_pct_error,
            'rustfoil_cd': rf['cd'],
            'xfoil_cd': xf['cd'],
            'cd_error': cd_error,
            'cd_pct_error': cd_pct_error,
            'converged': rf['converged']
        })
    
    return comparisons

def analyze_discrepancies(comparisons):
    """Identify patterns in discrepancies."""
    print("\n" + "="*80)
    print("ALPHA SWEEP COMPARISON: RustFoil vs XFOIL")
    print("="*80)
    print(f"{'Alpha':>6} {'RF CL':>10} {'XF CL':>10} {'ΔCL':>10} {'%Err':>8} "
          f"{'RF CD':>10} {'XF CD':>10} {'ΔCD':>10} {'%Err':>8}")
    print("-"*100)
    
    for c in comparisons:
        cl_err_str = f"{c['cl_pct_error']:>7.1f}%" if c['cl_pct_error'] is not None else "    N/A"
        cd_err_str = f"{c['cd_pct_error']:>7.1f}%" if c['cd_pct_error'] is not None else "    N/A"
        
        print(f"{c['alpha']:>6.1f} {c['rustfoil_cl']:>10.4f} {c['xfoil_cl']:>10.4f} "
              f"{c['cl_error']:>10.4f} {cl_err_str:>8} "
              f"{c['rustfoil_cd']:>10.6f} {c['xfoil_cd']:>10.6f} "
              f"{c['cd_error']:>10.6f} {cd_err_str:>8}")
    
    # Statistics
    print("\n" + "="*80)
    print("DISCREPANCY ANALYSIS")
    print("="*80)
    
    cl_errors = [c['cl_pct_error'] for c in comparisons if c['cl_pct_error'] is not None]
    cd_errors = [c['cd_pct_error'] for c in comparisons if c['cd_pct_error'] is not None]
    
    if cl_errors:
        print(f"\nCL Error Statistics:")
        print(f"  Mean error: {np.mean(cl_errors):>7.2f}%")
        print(f"  Std dev:    {np.std(cl_errors):>7.2f}%")
        print(f"  RMS error:  {np.sqrt(np.mean(np.array(cl_errors)**2)):>7.2f}%")
        print(f"  Max error:  {np.max(np.abs(cl_errors)):>7.2f}% at α={comparisons[np.argmax(np.abs(cl_errors))]['alpha']:.1f}°")
    
    if cd_errors:
        print(f"\nCD Error Statistics:")
        print(f"  Mean error: {np.mean(cd_errors):>7.2f}%")
        print(f"  Std dev:    {np.std(cd_errors):>7.2f}%")
        print(f"  RMS error:  {np.sqrt(np.mean(np.array(cd_errors)**2)):>7.2f}%")
        print(f"  Max error:  {np.max(np.abs(cd_errors)):>7.2f}% at α={comparisons[np.argmax(np.abs(cd_errors))]['alpha']:.1f}°")
    
    # Identify regions of concern
    print(f"\nREGIONS OF CONCERN:")
    print("-"*80)
    
    # Low alpha (near zero)
    low_alpha = [c for c in comparisons if abs(c['alpha']) <= 2]
    if low_alpha:
        cl_errs = [abs(c['cl_pct_error']) for c in low_alpha if c['cl_pct_error'] is not None]
        if cl_errs:
            print(f"\nLow α region (|α| ≤ 2°):")
            print(f"  Mean |CL error|: {np.mean(cl_errs):.2f}%")
            for c in low_alpha:
                if c['cl_pct_error'] and abs(c['cl_pct_error']) > 10:
                    print(f"    α={c['alpha']:>5.1f}°: CL error = {c['cl_pct_error']:>6.1f}%")
    
    # Mid alpha (3-8 degrees)
    mid_alpha = [c for c in comparisons if 3 <= c['alpha'] <= 8]
    if mid_alpha:
        cl_errs = [abs(c['cl_pct_error']) for c in mid_alpha if c['cl_pct_error'] is not None]
        if cl_errs:
            print(f"\nMid α region (3° ≤ α ≤ 8°):")
            print(f"  Mean |CL error|: {np.mean(cl_errs):.2f}%")
            for c in mid_alpha:
                if c['cl_pct_error'] and abs(c['cl_pct_error']) > 15:
                    print(f"    α={c['alpha']:>5.1f}°: CL error = {c['cl_pct_error']:>6.1f}%")
    
    # High alpha (>10 degrees)
    high_alpha = [c for c in comparisons if c['alpha'] > 10]
    if high_alpha:
        cl_errs = [abs(c['cl_pct_error']) for c in high_alpha if c['cl_pct_error'] is not None]
        if cl_errs:
            print(f"\nHigh α region (α > 10°):")
            print(f"  Mean |CL error|: {np.mean(cl_errs):.2f}%")
            for c in high_alpha:
                if c['cl_pct_error'] and abs(c['cl_pct_error']) > 20:
                    print(f"    α={c['alpha']:>5.1f}°: CL error = {c['cl_pct_error']:>6.1f}%")
    
    # Drag anomalies
    print(f"\nDrag Anomalies (CD error > 100%):")
    for c in comparisons:
        if c['cd_pct_error'] and abs(c['cd_pct_error']) > 100:
            print(f"  α={c['alpha']:>5.1f}°: CD={c['rustfoil_cd']:.6f} (XF: {c['xfoil_cd']:.6f}), error={c['cd_pct_error']:>7.1f}%")
    
    # Systematic bias
    print(f"\nSystematic Bias:")
    if cl_errors:
        positive_errors = len([e for e in cl_errors if e > 0])
        print(f"  CL: {positive_errors}/{len(cl_errors)} points are high ({100*positive_errors/len(cl_errors):.1f}%)")
        if positive_errors > 0.7 * len(cl_errors):
            print(f"  → RustFoil systematically OVER-predicts lift")
        elif positive_errors < 0.3 * len(cl_errors):
            print(f"  → RustFoil systematically UNDER-predicts lift")
    
    # Save comparison
    output_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/alpha_sweep_comparison.json")
    with open(output_file, 'w') as f:
        json.dump({
            'airfoil': 'naca0012',
            'reynolds': 3e6,
            'comparisons': comparisons,
            'statistics': {
                'cl_mean_error': float(np.mean(cl_errors)) if cl_errors else None,
                'cl_rms_error': float(np.sqrt(np.mean(np.array(cl_errors)**2))) if cl_errors else None,
                'cd_mean_error': float(np.mean(cd_errors)) if cd_errors else None,
                'cd_rms_error': float(np.sqrt(np.mean(np.array(cd_errors)**2))) if cd_errors else None,
            }
        }, f, indent=2)
    
    print(f"\nDetailed comparison saved to {output_file}")

def main():
    rustfoil_data, xfoil_data = load_data()
    comparisons = match_by_alpha(rustfoil_data['results'], xfoil_data['results'])
    analyze_discrepancies(comparisons)

if __name__ == '__main__':
    main()
