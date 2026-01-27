#!/usr/bin/env python3
"""
Final comprehensive comparison of RustFoil vs actual XFOIL data.
"""
import json
from pathlib import Path
import numpy as np

def main():
    # Load data
    with open("/Users/harry/flexfoil-boundary-layer/comparison_results/rustfoil_alpha_sweep.json") as f:
        rustfoil_data = json.load(f)
    
    with open("/Users/harry/flexfoil-boundary-layer/comparison_results/xfoil_actual_data.json") as f:
        xfoil_data = json.load(f)
    
    # Match by alpha
    rf_dict = {r['alpha']: r for r in rustfoil_data['results'] if r['success']}
    xf_dict = {r['alpha']: r for r in xfoil_data['results']}
    
    comparisons = []
    for alpha in sorted(set(rf_dict.keys()) & set(xf_dict.keys())):
        rf = rf_dict[alpha]
        xf = xf_dict[alpha]
        
        cl_error = rf['cl'] - xf['cl']
        cd_error = rf['cd'] - xf['cd']
        cl_pct = 100 * cl_error / abs(xf['cl']) if xf['cl'] != 0 else None
        cd_pct = 100 * cd_error / xf['cd'] if xf['cd'] != 0 else None
        
        comparisons.append({
            'alpha': alpha,
            'rf_cl': rf['cl'],
            'xf_cl': xf['cl'],
            'cl_error': cl_error,
            'cl_pct': cl_pct,
            'rf_cd': rf['cd'],
            'xf_cd': xf['cd'],
            'cd_error': cd_error,
            'cd_pct': cd_pct
        })
    
    # Print comparison table
    print("\n" + "="*100)
    print(f"{'COMPREHENSIVE COMPARISON: RustFoil vs XFOIL':^100}")
    print("="*100)
    print(f"{'Alpha':>6} │ {'RustFoil':^22} │ {'XFOIL':^22} │ {'Error':^22} │ {'% Error':^18}")
    print(f"{'(deg)':>6} │ {'CL':>10} {'CD':>10} │ {'CL':>10} {'CD':>10} │ {'ΔCL':>10} {'ΔCD':>10} │ {'CL%':>8} {'CD%':>8}")
    print("-"*100)
    
    for c in comparisons:
        cl_pct_str = f"{c['cl_pct']:>7.1f}%" if c['cl_pct'] is not None else "    N/A"
        cd_pct_str = f"{c['cd_pct']:>7.1f}%" if c['cd_pct'] is not None else "    N/A"
        
        # Color coding for errors
        cl_flag = "⚠" if c['cl_pct'] and abs(c['cl_pct']) > 20 else " "
        cd_flag = "⚠" if c['cd_pct'] and abs(c['cd_pct']) > 200 else " "
        
        print(f"{c['alpha']:>6.1f} │ {c['rf_cl']:>10.4f} {c['rf_cd']:>10.6f} │ "
              f"{c['xf_cl']:>10.4f} {c['xf_cd']:>10.6f} │ "
              f"{c['cl_error']:>10.4f} {c['cd_error']:>10.6f} │ "
              f"{cl_pct_str:>8} {cd_pct_str:>8} {cl_flag}{cd_flag}")
    
    # Statistics
    print("\n" + "="*100)
    print("STATISTICAL ANALYSIS")
    print("="*100)
    
    cl_errors = [c['cl_pct'] for c in comparisons if c['cl_pct'] is not None and abs(c['alpha']) >= 1]
    cd_errors = [c['cd_pct'] for c in comparisons if c['cd_pct'] is not None]
    
    print(f"\n{'Metric':<20} {'Mean':>12} {'Std Dev':>12} {'RMS':>12} {'Max |Error|':>15}")
    print("-"*75)
    if cl_errors:
        print(f"{'CL Error (%)':<20} {np.mean(cl_errors):>12.2f} {np.std(cl_errors):>12.2f} "
              f"{np.sqrt(np.mean(np.array(cl_errors)**2)):>12.2f} "
              f"{np.max(np.abs(cl_errors)):>15.2f}")
    if cd_errors:
        print(f"{'CD Error (%)':<20} {np.mean(cd_errors):>12.2f} {np.std(cd_errors):>12.2f} "
              f"{np.sqrt(np.mean(np.array(cd_errors)**2)):>12.2f} "
              f"{np.max(np.abs(cd_errors)):>15.2f}")
    
    # Regional analysis
    print(f"\n{'REGIONAL ANALYSIS':-^75}")
    
    regions = [
        ("Low α (1° ≤ |α| ≤ 3°)", lambda a: 1 <= abs(a) <= 3),
        ("Mid α (4° ≤ α ≤ 8°)", lambda a: 4 <= a <= 8),
        ("High α (9° ≤ α ≤ 12°)", lambda a: 9 <= a <= 12),
        ("Very High α (|α| > 12°)", lambda a: abs(a) > 12),
    ]
    
    for name, condition in regions:
        region_comps = [c for c in comparisons if condition(c['alpha'])]
        if region_comps:
            cl_errs = [c['cl_pct'] for c in region_comps if c['cl_pct'] is not None]
            cd_errs = [c['cd_pct'] for c in region_comps if c['cd_pct'] is not None]
            
            print(f"\n{name}:")
            if cl_errs:
                print(f"  CL: {len(cl_errs)} points, mean error = {np.mean(cl_errs):>6.1f}%, "
                      f"RMS = {np.sqrt(np.mean(np.array(cl_errs)**2)):>6.1f}%")
            if cd_errs:
                print(f"  CD: {len(cd_errs)} points, mean error = {np.mean(cd_errs):>6.1f}%, "
                      f"RMS = {np.sqrt(np.mean(np.array(cd_errs)**2)):>6.1f}%")
    
    # Key findings
    print(f"\n{'KEY FINDINGS':-^75}")
    
    # Systematic bias
    cl_positive = len([c for c in comparisons if c['cl_pct'] and c['cl_pct'] > 0 and abs(c['alpha']) >= 1])
    cl_total = len([c for c in comparisons if c['cl_pct'] is not None and abs(c['alpha']) >= 1])
    
    print(f"\n1. Systematic Bias:")
    print(f"   • {cl_positive}/{cl_total} points have CL > XFOIL ({100*cl_positive/cl_total:.0f}%)")
    if cl_positive > 0.7 * cl_total:
        print(f"   → RustFoil OVER-predicts lift across most of the range")
    
    # Worst cases
    worst_cl = max(comparisons, key=lambda c: abs(c['cl_pct']) if c['cl_pct'] else 0)
    worst_cd = max(comparisons, key=lambda c: abs(c['cd_pct']) if c['cd_pct'] else 0)
    
    print(f"\n2. Worst Discrepancies:")
    print(f"   • CL: α={worst_cl['alpha']:.1f}°, error={worst_cl['cl_pct']:.1f}%")
    print(f"   • CD: α={worst_cd['alpha']:.1f}°, error={worst_cd['cd_pct']:.1f}%")
    
    # Best region
    best_region = min(regions, key=lambda r: np.mean([abs(c['cl_pct']) for c in comparisons 
                                                       if r[1](c['alpha']) and c['cl_pct']]) 
                      if [c for c in comparisons if r[1](c['alpha']) and c['cl_pct']] else float('inf'))
    print(f"\n3. Best Agreement:")
    print(f"   • {best_region[0]}")
    
    # Anomalies
    print(f"\n4. Anomalies (CD error > 500%):")
    anomalies = [c for c in comparisons if c['cd_pct'] and abs(c['cd_pct']) > 500]
    for c in anomalies:
        print(f"   • α={c['alpha']:>5.1f}°: CD={c['rf_cd']:.6f} (XFOIL: {c['xf_cd']:.6f})")
    
    # Save
    output_file = Path("/Users/harry/flexfoil-boundary-layer/comparison_results/final_comparison.json")
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
    
    print(f"\n{'='*100}")
    print(f"✓ Full comparison saved to {output_file}")
    print("="*100)

if __name__ == '__main__':
    main()
