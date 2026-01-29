#!/usr/bin/env python3
"""
Systematic comparison of XFOIL vs RustFoil traces to find first point of divergence.
Checks every intermediate calculation stage in order.
"""

import json
import sys
import numpy as np
from pathlib import Path

def load_trace(path):
    """Load trace file and index events by subroutine."""
    with open(path) as f:
        data = json.load(f)
    return data.get('events', [])

def get_events_by_type(events, subroutine):
    """Get all events of a specific type."""
    return [e for e in events if e.get('subroutine') == subroutine]

def compare_arrays(name, xf_arr, rf_arr, tolerance=0.01, max_show=5):
    """Compare two arrays and return divergence info."""
    if len(xf_arr) != len(rf_arr):
        return {
            'match': False,
            'reason': f'Length mismatch: XFOIL={len(xf_arr)}, RustFoil={len(rf_arr)}',
            'first_diff_idx': 0
        }
    
    xf = np.array(xf_arr)
    rf = np.array(rf_arr)
    
    # Compute relative error where possible
    with np.errstate(divide='ignore', invalid='ignore'):
        rel_err = np.abs(xf - rf) / np.maximum(np.abs(xf), 1e-10)
    
    # Find first significant divergence
    mask = rel_err > tolerance
    if not mask.any():
        return {
            'match': True,
            'max_rel_err': float(np.nanmax(rel_err)),
            'mean_rel_err': float(np.nanmean(rel_err))
        }
    
    first_idx = int(np.argmax(mask))
    
    # Get worst offenders
    worst_indices = np.argsort(rel_err)[-max_show:][::-1]
    
    return {
        'match': False,
        'first_diff_idx': first_idx,
        'first_diff_val': {
            'xfoil': float(xf[first_idx]),
            'rustfoil': float(rf[first_idx]),
            'rel_err': float(rel_err[first_idx])
        },
        'max_rel_err': float(np.nanmax(rel_err)),
        'mean_rel_err': float(np.nanmean(rel_err)),
        'worst': [{
            'idx': int(i),
            'xfoil': float(xf[i]),
            'rustfoil': float(rf[i]),
            'rel_err': float(rel_err[i])
        } for i in worst_indices if not np.isnan(rel_err[i])]
    }

def compare_inviscid(xf_events, rf_events):
    """Compare inviscid solution (gamma, Cp)."""
    print("\n" + "="*80)
    print("STAGE 1: INVISCID SOLUTION")
    print("="*80)
    
    # Get full gamma from both
    xf_gamma_events = get_events_by_type(xf_events, 'FULLGAMMA')
    rf_gamma_events = get_events_by_type(rf_events, 'FULL_INVISCID')
    
    if not xf_gamma_events:
        print("WARNING: No FULLGAMMA event in XFOIL trace")
        return None
    if not rf_gamma_events:
        print("WARNING: No FULL_INVISCID event in RustFoil trace")
        return None
    
    xf_gamma = xf_gamma_events[0].get('gamma', [])
    rf_gamma = rf_gamma_events[0].get('gamma', [])
    
    print(f"\nGamma comparison (XFOIL n={len(xf_gamma)}, RustFoil n={len(rf_gamma)}):")
    result = compare_arrays('gamma', xf_gamma, rf_gamma, tolerance=0.01)
    
    if result['match']:
        print(f"  ✓ MATCH - max rel err: {result['max_rel_err']:.4%}")
    else:
        print(f"  ✗ DIVERGENCE at index {result['first_diff_idx']}")
        print(f"    XFOIL:    {result['first_diff_val']['xfoil']:.6f}")
        print(f"    RustFoil: {result['first_diff_val']['rustfoil']:.6f}")
        print(f"    Rel err:  {result['first_diff_val']['rel_err']:.4%}")
        return {'stage': 'INVISCID_GAMMA', 'details': result}
    
    # Compare Cp
    xf_cp = xf_gamma_events[0].get('cp', [])
    rf_cp = rf_gamma_events[0].get('cp', [])
    
    print(f"\nCp comparison:")
    cp_result = compare_arrays('cp', xf_cp, rf_cp, tolerance=0.01)
    
    if cp_result['match']:
        print(f"  ✓ MATCH - max rel err: {cp_result['max_rel_err']:.4%}")
    else:
        print(f"  ✗ DIVERGENCE at index {cp_result['first_diff_idx']}")
        return {'stage': 'INVISCID_CP', 'details': cp_result}
    
    return None  # No divergence found

def compare_initial_bl(xf_events, rf_events):
    """Compare initial BL state."""
    print("\n" + "="*80)
    print("STAGE 2: INITIAL BL STATE")
    print("="*80)
    
    # Get BL_INIT events
    xf_init = get_events_by_type(xf_events, 'BL_INIT')
    rf_init = get_events_by_type(rf_events, 'BL_INIT')
    
    if not xf_init:
        print("WARNING: No BL_INIT in XFOIL - checking BLVAR iter=0")
        xf_init = [e for e in xf_events if e.get('subroutine') == 'BLVAR' and e.get('iteration', 1) == 0]
    
    if not rf_init:
        print("WARNING: No BL_INIT in RustFoil")
        return None
    
    print(f"\nFound {len(xf_init)} XFOIL init events, {len(rf_init)} RustFoil init events")
    
    # Compare first station state
    # This is complex due to different event structures - skip for now
    return None

def compare_bl_iterations(xf_events, rf_events):
    """Compare BL state across iterations."""
    print("\n" + "="*80)
    print("STAGE 3: BL ITERATION COMPARISON")
    print("="*80)
    
    # Get full BL state events
    xf_bl = get_events_by_type(xf_events, 'FULL_BL_STATE')
    rf_bl = get_events_by_type(rf_events, 'FULL_BL_STATE')
    
    if not xf_bl:
        print("WARNING: No FULL_BL_STATE in XFOIL trace")
        return None
    if not rf_bl:
        print("WARNING: No FULL_BL_STATE in RustFoil trace")
        return None
    
    print(f"\nFound {len(xf_bl)} XFOIL iterations, {len(rf_bl)} RustFoil iterations")
    
    # Compare iteration 1 (first real iteration)
    for iter_num in range(min(len(xf_bl), len(rf_bl))):
        xf_state = xf_bl[iter_num]
        rf_state = rf_bl[iter_num]
        
        print(f"\n--- Iteration {iter_num + 1} ---")
        
        # Compare upper surface theta
        xf_theta_upper = xf_state.get('upper_theta', [])
        rf_upper = rf_state.get('upper', {})
        rf_theta_upper = rf_upper.get('theta', [])
        
        if xf_theta_upper and rf_theta_upper:
            result = compare_arrays('upper_theta', xf_theta_upper, rf_theta_upper, tolerance=0.05)
            if result['match']:
                print(f"  Upper θ: ✓ max err {result['max_rel_err']:.2%}")
            else:
                print(f"  Upper θ: ✗ DIVERGENCE at station {result['first_diff_idx']}")
                print(f"    XFOIL:    {result['first_diff_val']['xfoil']:.6e}")
                print(f"    RustFoil: {result['first_diff_val']['rustfoil']:.6e}")
                print(f"    Rel err:  {result['first_diff_val']['rel_err']:.2%}")
                return {'stage': f'BL_ITER{iter_num+1}_UPPER_THETA', 'iteration': iter_num+1, 'details': result}
        
        # Compare upper surface Ue
        xf_ue_upper = xf_state.get('upper_Ue', [])
        rf_ue_upper = rf_upper.get('ue', [])
        
        if xf_ue_upper and rf_ue_upper:
            result = compare_arrays('upper_Ue', xf_ue_upper, rf_ue_upper, tolerance=0.02)
            if result['match']:
                print(f"  Upper Ue: ✓ max err {result['max_rel_err']:.2%}")
            else:
                print(f"  Upper Ue: ✗ DIVERGENCE at station {result['first_diff_idx']}")
                return {'stage': f'BL_ITER{iter_num+1}_UPPER_UE', 'iteration': iter_num+1, 'details': result}
        
        # Compare upper surface Hk
        xf_hk_upper = xf_state.get('upper_Hk', [])
        rf_hk_upper = rf_upper.get('hk', [])
        
        if xf_hk_upper and rf_hk_upper:
            result = compare_arrays('upper_Hk', xf_hk_upper, rf_hk_upper, tolerance=0.05)
            if result['match']:
                print(f"  Upper Hk: ✓ max err {result['max_rel_err']:.2%}")
            else:
                print(f"  Upper Hk: ✗ DIVERGENCE at station {result['first_diff_idx']}")
                return {'stage': f'BL_ITER{iter_num+1}_UPPER_HK', 'iteration': iter_num+1, 'details': result}
        
        # Lower surface
        xf_theta_lower = xf_state.get('lower_theta', [])
        rf_lower = rf_state.get('lower', {})
        rf_theta_lower = rf_lower.get('theta', [])
        
        if xf_theta_lower and rf_theta_lower:
            result = compare_arrays('lower_theta', xf_theta_lower, rf_theta_lower, tolerance=0.05)
            if result['match']:
                print(f"  Lower θ: ✓ max err {result['max_rel_err']:.2%}")
            else:
                print(f"  Lower θ: ✗ DIVERGENCE at station {result['first_diff_idx']}")
                return {'stage': f'BL_ITER{iter_num+1}_LOWER_THETA', 'iteration': iter_num+1, 'details': result}
    
    return None

def compare_transition(xf_events, rf_events):
    """Compare transition prediction."""
    print("\n" + "="*80)
    print("STAGE 4: TRANSITION PREDICTION")
    print("="*80)
    
    # Get transition events
    xf_trans = get_events_by_type(xf_events, 'TRANSITION')
    rf_trans = get_events_by_type(rf_events, 'TRCHEK2_FINAL')
    
    print(f"\nFound {len(xf_trans)} XFOIL transition events, {len(rf_trans)} RustFoil transition events")
    
    # Compare final transition locations
    if xf_trans and rf_trans:
        # Get last transition for each surface
        xf_upper = [e for e in xf_trans if e.get('side', e.get('ISIDE')) == 1]
        xf_lower = [e for e in xf_trans if e.get('side', e.get('ISIDE')) == 2]
        rf_upper = [e for e in rf_trans if e.get('side') == 1]
        rf_lower = [e for e in rf_trans if e.get('side') == 2]
        
        if xf_upper and rf_upper:
            xf_xtr = xf_upper[-1].get('x_tr', xf_upper[-1].get('XTRAN'))
            rf_xtr = rf_upper[-1].get('xt_final')
            if xf_xtr and rf_xtr:
                rel_err = abs(xf_xtr - rf_xtr) / max(abs(xf_xtr), 0.01)
                print(f"\n  Upper x_tr: XFOIL={xf_xtr:.4f}, RustFoil={rf_xtr:.4f}, err={rel_err:.1%}")
                if rel_err > 0.1:
                    return {'stage': 'TRANSITION_UPPER', 'xfoil': xf_xtr, 'rustfoil': rf_xtr, 'rel_err': rel_err}
    
    return None

def compare_forces(xf_events, rf_events):
    """Compare final forces."""
    print("\n" + "="*80)
    print("STAGE 5: FORCES (CL, CD)")
    print("="*80)
    
    xf_cl = get_events_by_type(xf_events, 'CL_DETAIL')
    rf_cl = get_events_by_type(rf_events, 'CL_DETAIL')
    
    if xf_cl and rf_cl:
        xf_val = xf_cl[-1].get('cl')
        rf_val = rf_cl[-1].get('cl')
        if xf_val is not None and rf_val is not None:
            rel_err = abs(xf_val - rf_val) / max(abs(xf_val), 0.01)
            print(f"\n  CL: XFOIL={xf_val:.4f}, RustFoil={rf_val:.4f}, err={rel_err:.1%}")

def main():
    if len(sys.argv) < 3:
        print("Usage: python find_divergence.py <xfoil_trace.json> <rustfoil_trace.json>")
        sys.exit(1)
    
    xf_path = sys.argv[1]
    rf_path = sys.argv[2]
    
    print(f"Loading XFOIL trace: {xf_path}")
    xf_events = load_trace(xf_path)
    print(f"  {len(xf_events)} events")
    
    print(f"Loading RustFoil trace: {rf_path}")
    rf_events = load_trace(rf_path)
    print(f"  {len(rf_events)} events")
    
    # Run comparison stages in order
    divergence = None
    
    divergence = divergence or compare_inviscid(xf_events, rf_events)
    divergence = divergence or compare_initial_bl(xf_events, rf_events)
    divergence = divergence or compare_bl_iterations(xf_events, rf_events)
    divergence = divergence or compare_transition(xf_events, rf_events)
    compare_forces(xf_events, rf_events)
    
    print("\n" + "="*80)
    if divergence:
        print(f"FIRST DIVERGENCE FOUND: {divergence['stage']}")
        print("="*80)
        print(json.dumps(divergence, indent=2))
    else:
        print("NO SIGNIFICANT DIVERGENCE FOUND")
        print("="*80)

if __name__ == '__main__':
    main()
