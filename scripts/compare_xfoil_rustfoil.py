#!/usr/bin/env python3
"""
Compare XFOIL instrumented output with RustFoil.

Usage:
    python compare_xfoil_rustfoil.py xfoil_debug.json
    python compare_xfoil_rustfoil.py xfoil_debug.json --phase 1
    python compare_xfoil_rustfoil.py xfoil_debug.json --phase 3 --station 10
"""

import json
import argparse
import numpy as np
from pathlib import Path


def load_xfoil_events(json_path):
    """Load events from xfoil_debug.json."""
    with open(json_path) as f:
        return json.load(f)['events']


def filter_events(events, subroutine, iteration=None, side=None, ibl=None):
    """Filter events by criteria."""
    result = [e for e in events if e['subroutine'] == subroutine]
    if iteration is not None:
        result = [e for e in result if e.get('iteration') == iteration]
    if side is not None:
        result = [e for e in result if e.get('side') == side]
    if ibl is not None:
        result = [e for e in result if e.get('ibl') == ibl]
    return result


def summarize_events(events):
    """Print summary of event types."""
    counts = {}
    for e in events:
        sub = e['subroutine']
        counts[sub] = counts.get(sub, 0) + 1
    
    print("=" * 60)
    print("Event Summary")
    print("=" * 60)
    for k, v in sorted(counts.items(), key=lambda x: -x[1]):
        print(f"  {k:20s}: {v:6d}")
    print(f"  {'TOTAL':20s}: {len(events):6d}")
    print()


def phase1_closures(events):
    """Phase 1: Compare closure functions."""
    print("\n" + "=" * 60)
    print("Phase 1: Closure Function Values")
    print("=" * 60)
    
    # HKIN examples
    print("\n--- HKIN (H -> Hk) ---")
    hkin_events = filter_events(events, 'HKIN')[:10]
    print(f"{'H':>10s} {'Msq':>10s} {'Hk':>10s} {'Hk_H':>10s} {'Hk_Msq':>10s}")
    for ev in hkin_events:
        print(f"{ev['H']:10.4f} {ev['Msq']:10.4f} {ev['Hk']:10.4f} {ev['Hk_H']:10.4f} {ev['Hk_Msq']:10.4f}")
    
    # HSL examples
    print("\n--- HSL (Laminar Hs) ---")
    hsl_events = filter_events(events, 'HSL')[:10]
    print(f"{'Hk':>10s} {'Rtheta':>12s} {'Hs':>10s} {'Hs_Hk':>10s}")
    for ev in hsl_events:
        print(f"{ev['Hk']:10.4f} {ev['Rtheta']:12.2f} {ev['Hs']:10.4f} {ev['Hs_Hk']:10.4f}")
    
    # CFL examples
    print("\n--- CFL (Laminar Cf) ---")
    cfl_events = filter_events(events, 'CFL')[:10]
    print(f"{'Hk':>10s} {'Rtheta':>12s} {'Cf':>12s} {'Cf_Hk':>12s}")
    for ev in cfl_events:
        print(f"{ev['Hk']:10.4f} {ev['Rtheta']:12.2f} {ev['Cf']:12.6f} {ev['Cf_Hk']:12.6f}")
    
    # DAMPL examples
    print("\n--- DAMPL (Amplification Rate) ---")
    dampl_events = filter_events(events, 'DAMPL')[:10]
    print(f"{'Hk':>10s} {'theta':>12s} {'Rtheta':>12s} {'Ax':>12s}")
    for ev in dampl_events:
        print(f"{ev['Hk']:10.4f} {ev['theta']:12.6f} {ev['Rtheta']:12.2f} {ev['Ax']:12.4f}")


def phase2_blvar(events, station=None):
    """Phase 2: BLVAR secondary variables."""
    print("\n" + "=" * 60)
    print("Phase 2: BLVAR Secondary Variables")
    print("=" * 60)
    
    blvar_events = filter_events(events, 'BLVAR', iteration=1)
    if station is not None:
        blvar_events = [e for e in blvar_events if e.get('ibl') == station]
    
    for ev in blvar_events[:20]:
        print(f"\n--- Side {ev.get('side', '?')}, Station {ev.get('ibl', '?')}, Type {ev.get('flow_type', '?')} ---")
        inp = ev.get('input', {})
        out = ev.get('output', {})
        
        print("  Inputs:")
        print(f"    x        = {inp.get('x', 0):12.6e}")
        print(f"    u (Ue)   = {inp.get('u', 0):12.6e}")
        print(f"    theta    = {inp.get('theta', 0):12.6e}")
        print(f"    delta*   = {inp.get('delta_star', 0):12.6e}")
        print(f"    ctau     = {inp.get('ctau', 0):12.6e}")
        
        print("  Outputs:")
        print(f"    H        = {out.get('H', 0):12.6f}")
        print(f"    Hk       = {out.get('Hk', 0):12.6f}")
        print(f"    Hs       = {out.get('Hs', 0):12.6f}")
        print(f"    Cf       = {out.get('Cf', 0):12.6e}")
        print(f"    Cd       = {out.get('Cd', 0):12.6e}")
        print(f"    Rtheta   = {out.get('Rtheta', 0):12.2f}")


def phase3_jacobians(
    events,
    station=None,
    rustfoil_debug_path=None,
    xfoil_side=None,
    rust_side=None,
    rust_ibl=None,
):
    """Phase 3: BLDIF Jacobian analysis."""
    print("\n" + "=" * 60)
    print("Phase 3: BLDIF Jacobian Analysis")
    print("=" * 60)
    
    bldif_events = filter_events(events, 'BLDIF', iteration=1)
    if xfoil_side is not None:
        bldif_events = [e for e in bldif_events if e.get('side') == xfoil_side]
    if station is not None:
        bldif_events = [e for e in bldif_events if e.get('ibl') == station]

    # Keep only the final Newton iteration per station when available
    if bldif_events:
        grouped = {}
        for ev in bldif_events:
            key = (ev.get('side'), ev.get('ibl'), ev.get('flow_type'))
            iter_val = ev.get('newton_iter', ev.get('iteration', 0))
            if key not in grouped or iter_val >= grouped[key].get('newton_iter', grouped[key].get('iteration', 0)):
                grouped[key] = ev
        bldif_events = list(grouped.values())

    rust_events = []
    if rustfoil_debug_path and Path(rustfoil_debug_path).exists():
        with open(rustfoil_debug_path) as f:
            rust_events = json.load(f).get('events', [])
    rust_bldif = [e for e in rust_events if e.get('subroutine') == 'BLDIF']
    rust_map = {}
    for ev in rust_bldif:
        key = (ev.get('side'), ev.get('ibl'), ev.get('flow_type'))
        rust_map[key] = ev
    
    if rust_bldif:
        print(f"{'Side':>4} {'IBL':>4} {'Type':>4} {'Iter':>4} {'Max|dVS1|':>12} {'Max|dVS2|':>12} {'Max|dVSREZ|':>12}")

    for ev in bldif_events[:10]:
        print(f"\n--- Side {ev.get('side', '?')}, Station {ev.get('ibl', '?')}, Type {ev.get('flow_type', '?')} ---")
        
        vs1 = np.array(ev['VS1'])
        vs2 = np.array(ev['VS2'])
        vsrez = np.array(ev['VSREZ'])

        if rust_bldif:
            x_side = ev.get('side')
            x_ibl = ev.get('ibl')
            x_type = ev.get('flow_type')
            r_side = rust_side if rust_side is not None else x_side
            r_ibl = rust_ibl if rust_ibl is not None else x_ibl
            rust = rust_map.get((r_side, r_ibl, x_type))
            if rust:
                r_vs1 = np.array(rust['VS1'])
                r_vs2 = np.array(rust['VS2'])
                r_vsrez = np.array(rust['VSREZ'])
                max_vs1 = np.max(np.abs(vs1 - r_vs1))
                max_vs2 = np.max(np.abs(vs2 - r_vs2))
                max_vsr = np.max(np.abs(vsrez - r_vsrez))
                iter_val = ev.get('newton_iter', ev.get('iteration', 0))
                print(f"{x_side:>4} {x_ibl:>4} {x_type:>4} {iter_val:>4} {max_vs1:12.4e} {max_vs2:12.4e} {max_vsr:12.4e}")
            else:
                print("  RustFoil BLDIF not found for matching side/ibl/flow_type")
        
        print(f"  VS1 (upstream station derivatives):")
        print(f"    Shape: {vs1.shape}")
        for i in range(min(3, vs1.shape[0])):
            print(f"    Row {i}: [{', '.join(f'{v:12.4e}' for v in vs1[i])}]")
        
        print(f"  VS2 (current station derivatives):")
        for i in range(min(3, vs2.shape[0])):
            print(f"    Row {i}: [{', '.join(f'{v:12.4e}' for v in vs2[i])}]")
        
        print(f"  VSREZ (residuals): [{', '.join(f'{v:12.4e}' for v in vsrez)}]")
        
        # Check for problematic patterns
        warnings = []
        
        # Check column 2 (delta* derivatives) - was previously zero in RustFoil
        if vs1.shape[1] > 2:
            col2_vs1 = vs1[:3, 2] if vs1.shape[0] >= 3 else vs1[:, 2]
            col2_vs2 = vs2[:3, 2] if vs2.shape[0] >= 3 else vs2[:, 2]
            if np.allclose(col2_vs1, 0, atol=1e-20):
                warnings.append("VS1 column 2 (δ* derivs) is all zeros!")
            if np.allclose(col2_vs2, 0, atol=1e-20):
                warnings.append("VS2 column 2 (δ* derivs) is all zeros!")
        
        # Check 3x3 VA block determinant
        if vs1.shape[0] >= 3 and vs1.shape[1] >= 3:
            va = vs1[:3, :3]
            try:
                det = np.linalg.det(va)
                if abs(det) < 1e-10:
                    warnings.append(f"Near-singular VA block! det={det:.2e}")
                else:
                    print(f"  VA block det: {det:.4e} (OK)")
            except:
                warnings.append("Could not compute VA determinant")
        
        for w in warnings:
            print(f"  ⚠️  WARNING: {w}")


def phase4_march(events):
    """Phase 4: Initial march state."""
    print("\n" + "=" * 60)
    print("Phase 4: MRCHUE Initial March")
    print("=" * 60)
    
    mrchue_events = filter_events(events, 'MRCHUE')
    
    # Group by side
    side1 = [e for e in mrchue_events if e.get('side') == 1]
    side2 = [e for e in mrchue_events if e.get('side') == 2]
    
    for side, side_events in [(1, side1), (2, side2)]:
        print(f"\n--- Side {side} ({len(side_events)} stations) ---")
        print(f"{'IBL':>4s} {'x':>12s} {'Ue':>12s} {'theta':>12s} {'delta*':>12s} {'Hk':>8s} {'Cf':>12s}")
        for ev in side_events[:15]:
            print(f"{ev.get('ibl', 0):4d} {ev['x']:12.6f} {ev['Ue']:12.6f} "
                  f"{ev['theta']:12.6e} {ev['delta_star']:12.6e} "
                  f"{ev['Hk']:8.4f} {ev['Cf']:12.6e}")


def phase5_iteration(events):
    """Phase 5: Newton iteration convergence."""
    print("\n" + "=" * 60)
    print("Phase 5: Iteration Convergence")
    print("=" * 60)
    
    viscal_results = filter_events(events, 'VISCAL_RESULT')
    
    print(f"{'Iter':>4s} {'RMS':>12s} {'Max':>12s} {'CL':>10s} {'CD':>12s} {'CM':>10s}")
    for ev in viscal_results:
        print(f"{ev['iteration']:4d} {ev['rms_residual']:12.4e} {ev['max_residual']:12.4e} "
              f"{ev['CL']:10.6f} {ev['CD']:12.6e} {ev['CM']:10.6f}")
    
    # UPDATE statistics per iteration
    print("\n--- UPDATE deltas (first iteration) ---")
    update_events = filter_events(events, 'UPDATE', iteration=1)[:10]
    print(f"{'Side':>4s} {'IBL':>4s} {'Δctau':>12s} {'Δtheta':>12s} {'Δmass':>12s} {'ΔUe':>12s} {'RLX':>8s}")
    for ev in update_events:
        print(f"{ev['side']:4d} {ev['ibl']:4d} {ev['delta_ctau']:12.4e} {ev['delta_theta']:12.4e} "
              f"{ev['delta_mass']:12.4e} {ev['delta_Ue']:12.4e} {ev['relaxation']:8.4f}")


def phase6_dij(events):
    """Phase 6: DIJ influence matrix."""
    print("\n" + "=" * 60)
    print("Phase 6: QDCALC DIJ Matrix")
    print("=" * 60)
    
    qdcalc_events = filter_events(events, 'QDCALC')
    
    for ev in qdcalc_events:
        print(f"  n_airfoil: {ev['n_airfoil']}")
        print(f"  n_wake:    {ev['n_wake']}")
        print(f"  n_total:   {ev['n_total']}")
        
        diag = ev.get('DIJ_diagonal_sample', [])
        row1 = ev.get('DIJ_row1_sample', [])
        
        print(f"  Diagonal (first 10): {[f'{v:.4e}' for v in diag]}")
        print(f"  Row 1 (first 10):    {[f'{v:.4e}' for v in row1]}")


def phase8_trchek2(
    events,
    rustfoil_debug_path=None,
    xfoil_side=None,
    rust_side=None,
    ibl=None,
    rust_ibl=None,
):
    """Phase 8: TRCHEK2 transition-point comparison."""
    print("\n" + "=" * 60)
    print("Phase 8: TRCHEK2 Transition Comparison")
    print("=" * 60)

    tr_final = filter_events(events, 'TRCHEK2_FINAL')
    tr_iter = filter_events(events, 'TRCHEK2_ITER')

    if not tr_final:
        print("No TRCHEK2_FINAL events found in XFOIL data")
        return

    # Index the latest TRCHEK2_ITER by (side, ibl) and by call_id
    tr_iter_map = {}
    tr_iter_by_call = {}
    for ev in tr_iter:
        side = ev.get('side')
        if xfoil_side is not None and side != xfoil_side:
            continue
        if ibl is not None and ev.get('ibl') != ibl:
            continue
        key = (side, ev.get('ibl'))
        if key not in tr_iter_map or ev.get('trchek_iter', 0) >= tr_iter_map[key].get('trchek_iter', 0):
            tr_iter_map[key] = ev
        call_key = (side, ev.get('ibl'), ev.get('call_id'))
        if call_key not in tr_iter_by_call or ev.get('trchek_iter', 0) >= tr_iter_by_call[call_key].get('trchek_iter', 0):
            tr_iter_by_call[call_key] = ev

    # Load RustFoil debug if available
    rust_events = []
    if rustfoil_debug_path and Path(rustfoil_debug_path).exists():
        with open(rustfoil_debug_path) as f:
            rust_events = json.load(f).get('events', [])

    rust_tr_final = [e for e in rust_events if e.get('subroutine') == 'TRCHEK2_FINAL']
    rust_tr_iter = [e for e in rust_events if e.get('subroutine') == 'TRCHEK2_ITER']
    rust_map = {(e.get('side'), e.get('ibl')): e for e in rust_tr_final}
    rust_iter_map = {}
    rust_iter_by_call = {}
    for ev in rust_tr_iter:
        side = ev.get('side')
        if rust_side is not None and side != rust_side:
            continue
        if rust_ibl is not None and ev.get('ibl') != rust_ibl:
            continue
        key = (side, ev.get('ibl'))
        if key not in rust_iter_map or ev.get('trchek_iter', 0) >= rust_iter_map[key].get('trchek_iter', 0):
            rust_iter_map[key] = ev
        call_key = (side, ev.get('ibl'), ev.get('call_id'))
        if call_key not in rust_iter_by_call or ev.get('trchek_iter', 0) >= rust_iter_by_call[call_key].get('trchek_iter', 0):
            rust_iter_by_call[call_key] = ev

    xfoil_mrchue = [e for e in events if e.get('subroutine') == 'MRCHUE_ITER']
    rust_mrchue = [e for e in rust_events if e.get('subroutine') == 'MRCHUE_ITER']
    xfoil_mrchue_map = {}
    for ev in xfoil_mrchue:
        side = ev.get('side')
        if xfoil_side is not None and side != xfoil_side:
            continue
        if ibl is not None and ev.get('ibl') != ibl:
            continue
        key = (side, ev.get('ibl'))
        if key not in xfoil_mrchue_map or ev.get('newton_iter', 0) >= xfoil_mrchue_map[key].get('newton_iter', 0):
            xfoil_mrchue_map[key] = ev

    rust_mrchue_map = {}
    for ev in rust_mrchue:
        side = ev.get('side')
        if rust_side is not None and side != rust_side:
            continue
        if rust_ibl is not None and ev.get('ibl') != rust_ibl:
            continue
        key = (side, ev.get('ibl'))
        iter_val = ev.get('newton_iter', ev.get('iter', ev.get('iteration', 0)))
        if key not in rust_mrchue_map or iter_val >= rust_mrchue_map[key].get('newton_iter', rust_mrchue_map[key].get('iter', 0)):
            rust_mrchue_map[key] = ev

    header = [
        "Side", "IBL", "XT_XFOIL", "XT_RUST", "TT_XFOIL", "TT_RUST", "Err%_TT"
    ]
    print(f"{header[0]:>4} {header[1]:>4} {header[2]:>12} {header[3]:>12} {header[4]:>12} {header[5]:>12} {header[6]:>8}")

    if ibl is not None:
        compare_events = [
            e for e in tr_final
            if (xfoil_side is None or e.get('side') == xfoil_side) and e.get('ibl') == ibl
        ]
        if compare_events:
            compare_events = [
                max(compare_events, key=lambda ev: ev.get('trchek_iter', ev.get('iteration', 0)))
            ]
    elif rust_tr_final:
        compare_keys = sorted(rust_map.keys())
        tr_final_map = {}
        tr_final_by_global = {}
        for ev in tr_final:
            key = (ev.get('side'), ev.get('ibl'))
            tr_final_map[key] = ev
            gkey = (ev.get('side'), ev.get('ibl'), ev.get('global_iter'))
            tr_final_by_global[gkey] = ev
        compare_events = []
        for key in compare_keys:
            rust_ev = rust_map.get(key, {})
            g_iter = rust_ev.get('global_iter')
            if g_iter is not None:
                pick = tr_final_by_global.get((key[0], key[1], g_iter))
                if pick is not None:
                    compare_events.append(pick)
                    continue
            if tr_final_map.get(key):
                compare_events.append(tr_final_map.get(key))
    else:
        compare_events = tr_final[:30]

    def pick_iter_event(iter_events, xt_final):
        if not iter_events:
            return {}
        # Prefer events closest to final XT, then highest trchek_iter
        best = None
        best_key = None
        for it in iter_events:
            xt_val = it.get('xt')
            if xt_val is None or xt_final is None:
                continue
            key = (abs(xt_val - xt_final), -it.get('trchek_iter', 0))
            if best_key is None or key < best_key:
                best_key = key
                best = it
        return best or max(iter_events, key=lambda it: it.get('trchek_iter', 0))

    for ev in compare_events:
        side = ev.get('side')
        ibl_val = ev.get('ibl')
        if xfoil_side is not None and side != xfoil_side:
            continue
        if ibl is not None and ibl_val != ibl:
            continue
        key = (side, ibl_val)

        xt_xfoil = ev.get('xt_final', ev.get('xt'))
        iter_candidates = [e for e in tr_iter if e.get('side') == side and e.get('ibl') == ibl_val]
        iter_ev = tr_iter_by_call.get((side, ibl_val, ev.get('call_id')))
        if not iter_ev:
            iter_ev = pick_iter_event(iter_candidates, xt_xfoil) or tr_iter_map.get(key, {})
        wf1 = iter_ev.get('wf1')
        wf2 = iter_ev.get('wf2')
        t1 = iter_ev.get('T1')
        t2 = iter_ev.get('T2')
        tt_xfoil = None
        if wf1 is not None and wf2 is not None and t1 is not None and t2 is not None:
            tt_xfoil = t1 * wf1 + t2 * wf2

        rust_ibl_val = rust_ibl if rust_ibl is not None else ibl_val
        rust_key = (rust_side if rust_side is not None else side, rust_ibl_val)
        rust = rust_map.get(rust_key, {})
        xt_rust = rust.get('xt_final')
        rust_iter_candidates = [e for e in rust_tr_iter if e.get('side') == rust_key[0] and e.get('ibl') == rust_key[1]]
        rust_iter_ev = rust_iter_by_call.get((rust_key[0], rust_key[1], rust.get('call_id')))
        if not rust_iter_ev:
            rust_iter_ev = pick_iter_event(rust_iter_candidates, xt_rust) or rust_iter_map.get(rust_key, {})
        t1_r = rust_iter_ev.get('T1')
        t2_r = rust_iter_ev.get('T2')
        wf1_r_iter = rust_iter_ev.get('wf1')
        wf2_r_iter = rust_iter_ev.get('wf2')
        wf1_r_final = rust.get('wf1')
        wf2_r_final = rust.get('wf2')
        tt_rust = None
        if wf1_r_final is not None and wf2_r_final is not None and t1_r is not None and t2_r is not None:
            tt_rust = t1_r * wf1_r_final + t2_r * wf2_r_final
        elif wf1_r_iter is not None and wf2_r_iter is not None and t1_r is not None and t2_r is not None:
            tt_rust = t1_r * wf1_r_iter + t2_r * wf2_r_iter
        elif rust.get('tt') is not None:
            tt_rust = rust.get('tt')

        err_tt = 0.0
        if tt_xfoil and tt_rust and abs(tt_xfoil) > 1e-15:
            err_tt = abs(tt_xfoil - tt_rust) / abs(tt_xfoil) * 100.0

        print(f"{side:>4} {ibl_val:>4} {xt_xfoil:12.6e} {xt_rust if xt_rust is not None else 0.0:12.6e} "
              f"{tt_xfoil if tt_xfoil is not None else 0.0:12.6e} {tt_rust if tt_rust is not None else 0.0:12.6e} "
              f"{err_tt:8.2f}")
        if wf1 is not None or wf1_r_final is not None or wf1_r_iter is not None:
            print(
                f"  WF XFOIL: {wf1} / {wf2} | "
                f"WF Rust(final): {wf1_r_final} / {wf2_r_final} | "
                f"WF Rust(iter): {wf1_r_iter} / {wf2_r_iter}"
            )

        if ibl is not None or rust_ibl is not None:
            x_mrchue = xfoil_mrchue_map.get(key, {})
            r_mrchue = rust_mrchue_map.get(rust_key, {})
            x_theta2 = None
            r_theta2 = None
            if isinstance(x_mrchue.get('output'), dict):
                x_theta2 = x_mrchue['output'].get('theta')
            if isinstance(r_mrchue.get('output'), dict):
                r_theta2 = r_mrchue['output'].get('theta')
            if r_theta2 is None:
                r_theta2 = r_mrchue.get('theta')
            err_theta2 = 0.0
            if x_theta2 and r_theta2 and abs(x_theta2) > 1e-15:
                err_theta2 = abs(x_theta2 - r_theta2) / abs(x_theta2) * 100.0

            x_ue = x_mrchue.get('Ue')
            r_ue = r_mrchue.get('Ue')
            err_ue = 0.0
            if x_ue and r_ue and abs(x_ue) > 1e-15:
                err_ue = abs(x_ue - r_ue) / abs(x_ue) * 100.0

            print("\nTransition station detail:")
            print(f"  XFOIL theta2: {x_theta2} | Rust theta2: {r_theta2} | Err%: {err_theta2:.2f}")
            print(f"  XFOIL Ue: {x_ue} | Rust Ue: {r_ue} | Err%: {err_ue:.2f}")

    if rust_tr_final:
        print("\nRustFoil TRCHEK2_FINAL extra fields (last 5):")
        for ev in rust_tr_final[-5:]:
            print(f"  side={ev.get('side')} ibl={ev.get('ibl')} "
                  f"wf1={ev.get('wf1')} wf2={ev.get('wf2')} tt={ev.get('tt')} "
                  f"dt={ev.get('dt')} ut={ev.get('ut')} Hk_t={ev.get('Hk_t')} "
                  f"Rt_t={ev.get('Rt_t')} St={ev.get('St')} Cq_t={ev.get('Cq_t')}")


def phase9_trdif(
    events,
    rustfoil_debug_path=None,
    xfoil_side=None,
    rust_side=None,
    ibl=None,
    rust_ibl=None,
):
    """Phase 9: TRDIF matrix comparison."""
    print("\n" + "=" * 60)
    print("Phase 9: TRDIF Matrix Comparison")
    print("=" * 60)

    trdif_events = [e for e in events if e.get('subroutine') == 'TRDIF']
    if xfoil_side is not None:
        trdif_events = [e for e in trdif_events if e.get('side') == xfoil_side]
    if ibl is not None:
        trdif_events = [e for e in trdif_events if e.get('ibl') == ibl]

    if not trdif_events:
        print("No TRDIF events found in XFOIL data")
        return

    # Keep only final Newton iteration per station
    grouped = {}
    for ev in trdif_events:
        key = (ev.get('side'), ev.get('ibl'))
        iter_val = ev.get('newton_iter', ev.get('iteration', 0))
        if key not in grouped or iter_val >= grouped[key].get('newton_iter', grouped[key].get('iteration', 0)):
            grouped[key] = ev
    trdif_events = list(grouped.values())

    rust_events = []
    if rustfoil_debug_path and Path(rustfoil_debug_path).exists():
        with open(rustfoil_debug_path) as f:
            rust_events = json.load(f).get('events', [])
    rust_trdif = [e for e in rust_events if e.get('subroutine') == 'TRDIF']
    rust_map = {}
    for ev in rust_trdif:
        key = (ev.get('side'), ev.get('ibl'))
        iter_val = ev.get('iter', ev.get('iteration', 0))
        if key not in rust_map or iter_val >= rust_map[key].get('iter', rust_map[key].get('iteration', 0)):
            rust_map[key] = ev

    print(f"{'Side':>4} {'IBL':>4} {'Iter':>4} {'Max|dVS1|':>12} {'Max|dVS2|':>12} {'Max|dVSREZ|':>12}")

    for ev in trdif_events:
        x_side = ev.get('side')
        x_ibl = ev.get('ibl')
        r_side = rust_side if rust_side is not None else x_side
        r_ibl = rust_ibl if rust_ibl is not None else x_ibl
        rust = rust_map.get((r_side, r_ibl))
        if not rust:
            print(f"{x_side:>4} {x_ibl:>4} {'?':>4} {'-':>12} {'-':>12} {'-':>12}")
            continue
        vs1 = np.array(ev['VS1'])
        vs2 = np.array(ev['VS2'])
        vsrez = np.array(ev['VSREZ'])
        r_vs1 = np.array(rust['VS1'])
        r_vs2 = np.array(rust['VS2'])
        r_vsrez = np.array(rust['VSREZ'])
        max_vs1 = np.max(np.abs(vs1 - r_vs1))
        max_vs2 = np.max(np.abs(vs2 - r_vs2))
        max_vsr = np.max(np.abs(vsrez - r_vsrez))
        iter_val = ev.get('newton_iter', ev.get('iteration', 0))
        print(f"{x_side:>4} {x_ibl:>4} {iter_val:>4} {max_vs1:12.4e} {max_vs2:12.4e} {max_vsr:12.4e}")

        if ibl is not None or rust_ibl is not None:
            def top_diffs(label, a, b, n=5):
                diffs = np.abs(a - b)
                flat_idx = np.argsort(diffs, axis=None)[-n:][::-1]
                print(f"  Top {n} {label} diffs:")
                for idx in flat_idx:
                    r, c = np.unravel_index(idx, diffs.shape)
                    print(f"    [{r},{c}] xfoil={a[r,c]:.6e} rust={b[r,c]:.6e} diff={diffs[r,c]:.6e}")

            top_diffs("VS1", vs1, r_vs1)
            top_diffs("VS2", vs2, r_vs2)
            diffs = np.abs(vsrez - r_vsrez)
            order = np.argsort(diffs)[::-1][:4]
            print("  VSREZ diffs:")
            for i in order:
                print(f"    [{i}] xfoil={vsrez[i]:.6e} rust={r_vsrez[i]:.6e} diff={diffs[i]:.6e}")


def phase10_trdif_derivs(
    events,
    rustfoil_debug_path=None,
    xfoil_side=None,
    rust_side=None,
    ibl=None,
    rust_ibl=None,
):
    """Phase 10: TRDIF derivative comparison."""
    print("\n" + "=" * 60)
    print("Phase 10: TRDIF_DERIVS Comparison")
    print("=" * 60)

    deriv_events = [e for e in events if e.get('subroutine') == 'TRDIF_DERIVS']
    if xfoil_side is not None:
        deriv_events = [e for e in deriv_events if e.get('side') == xfoil_side]
    if ibl is not None:
        deriv_events = [e for e in deriv_events if e.get('ibl') == ibl]

    if not deriv_events:
        print("No TRDIF_DERIVS events found in XFOIL data")
        return

    # Final Newton iteration per station
    grouped = {}
    for ev in deriv_events:
        key = (ev.get('side'), ev.get('ibl'))
        iter_val = ev.get('newton_iter', ev.get('iteration', 0))
        if key not in grouped or iter_val >= grouped[key].get('newton_iter', grouped[key].get('iteration', 0)):
            grouped[key] = ev
    deriv_events = list(grouped.values())

    rust_events = []
    if rustfoil_debug_path and Path(rustfoil_debug_path).exists():
        with open(rustfoil_debug_path) as f:
            rust_events = json.load(f).get('events', [])
    rust_derivs = [e for e in rust_events if e.get('subroutine') == 'TRDIF_DERIVS']
    rust_map = {}
    for ev in rust_derivs:
        key = (ev.get('side'), ev.get('ibl'))
        iter_val = ev.get('iter', ev.get('iteration', 0))
        if key not in rust_map or iter_val >= rust_map[key].get('iter', rust_map[key].get('iteration', 0)):
            rust_map[key] = ev

    ignore_keys = {
        'subroutine', 'call_id', 'side', 'ibl', 'iteration', 'newton_iter', 'iter'
    }

    for ev in deriv_events:
        x_side = ev.get('side')
        x_ibl = ev.get('ibl')
        r_side = rust_side if rust_side is not None else x_side
        r_ibl = rust_ibl if rust_ibl is not None else x_ibl
        rust = rust_map.get((r_side, r_ibl))
        if not rust:
            print(f"No Rust TRDIF_DERIVS for side={r_side} ibl={r_ibl}")
            continue

        diffs = []
        for k, v in ev.items():
            if k in ignore_keys:
                continue
            if k not in rust:
                continue
            rv = rust[k]
            if not isinstance(v, (int, float)) or not isinstance(rv, (int, float)):
                continue
            diff = rv - v
            rel = abs(diff) / abs(v) if abs(v) > 1e-12 else abs(diff)
            diffs.append((k, v, rv, diff, rel))

        diffs.sort(key=lambda x: abs(x[3]), reverse=True)
        print(f"\n--- Side {x_side}, IBL {x_ibl} (Rust {r_side}/{r_ibl}) ---")
        print(f"{'Key':>8} {'XFOIL':>14} {'Rust':>14} {'Diff':>14} {'Rel':>12}")
        for k, v, rv, diff, rel in diffs[:10]:
            print(f"{k:>8} {v:14.6e} {rv:14.6e} {diff:14.6e} {rel:12.4e}")


def export_test_vectors(events, output_dir):
    """Export test vectors for RustFoil unit tests."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)
    
    # Export closure function test cases
    hkin_events = filter_events(events, 'HKIN')[:50]
    with open(output_dir / 'hkin_test_vectors.json', 'w') as f:
        json.dump(hkin_events, f, indent=2)
    
    hsl_events = filter_events(events, 'HSL')[:50]
    with open(output_dir / 'hsl_test_vectors.json', 'w') as f:
        json.dump(hsl_events, f, indent=2)
    
    cfl_events = filter_events(events, 'CFL')[:50]
    with open(output_dir / 'cfl_test_vectors.json', 'w') as f:
        json.dump(cfl_events, f, indent=2)
    
    blvar_events = filter_events(events, 'BLVAR', iteration=1)[:50]
    with open(output_dir / 'blvar_test_vectors.json', 'w') as f:
        json.dump(blvar_events, f, indent=2)
    
    # Export DAMPL - mix of zero and non-zero Ax values
    all_dampl = filter_events(events, 'DAMPL')
    dampl_zero = [e for e in all_dampl if e.get('Ax', 0) == 0][:25]
    dampl_nonzero = [e for e in all_dampl if e.get('Ax', 0) != 0][:25]
    dampl_events = dampl_zero + dampl_nonzero
    with open(output_dir / 'dampl_test_vectors.json', 'w') as f:
        json.dump(dampl_events, f, indent=2)
    
    print(f"Test vectors exported to {output_dir}")


def phase_mrchue_stations(xfoil_events, rustfoil_debug_path=None):
    """
    Compare theta/dstar evolution station by station.
    
    This compares XFOIL's MRCHUE_ITER final values with RustFoil's march results.
    """
    print("=" * 70)
    print("MRCHUE Station-by-Station Comparison")
    print("=" * 70)
    
    # Group MRCHUE_ITER events by side and ibl
    mrchue_events = [e for e in xfoil_events if e.get('subroutine') == 'MRCHUE_ITER']
    
    if not mrchue_events:
        print("No MRCHUE_ITER events found in XFOIL data")
        return
    
    # Organize by side -> ibl -> final iteration
    stations = {}
    for e in mrchue_events:
        side = e['side']
        ibl = e['ibl']
        key = (side, ibl)
        
        # Keep only the final (converged) iteration for each station
        if key not in stations or e['newton_iter'] > stations[key]['newton_iter']:
            stations[key] = e
    
    # Load RustFoil debug if available
    rust_stations = {}
    rust_data = None
    if rustfoil_debug_path and Path(rustfoil_debug_path).exists():
        with open(rustfoil_debug_path) as f:
            rust_data = json.load(f)
            # Extract best (lowest dmax) MRCHUE_ITER per station
            for e in rust_data.get('events', []):
                if e.get('subroutine') == 'MRCHUE_ITER':
                    side = e['side']
                    ibl = e['ibl']
                    key = (side, ibl)
                    dmax = e.get('dmax')
                    if key not in rust_stations:
                        rust_stations[key] = e
                    else:
                        prev_dmax = rust_stations[key].get('dmax')
                        if dmax is not None and (prev_dmax is None or dmax < prev_dmax):
                            rust_stations[key] = e

            # Optional arc-length sanity check to catch mismatched geometry/paneling.
            xfoil_mrchue = {(e['side'], e['ibl']): e for e in xfoil_events if e.get('subroutine') == 'MRCHUE'}
            rust_mrchue = {
                (e['side'], e['ibl']): e
                for e in rust_data.get('events', [])
                if e.get('subroutine') == 'MRCHUE'
            }
            dx_samples = []
            for key, x_ev in xfoil_mrchue.items():
                r_ev = rust_mrchue.get(key)
                if not r_ev:
                    continue
                x_xfoil = x_ev.get('x')
                x_rust = r_ev.get('x')
                if x_xfoil is None or x_rust is None:
                    continue
                dx_samples.append(x_rust - x_xfoil)
            if dx_samples:
                max_dx = max(dx_samples, key=lambda v: abs(v))
                if abs(max_dx) > 1e-3:
                    print(
                        f\"\\n[warn] MRCHUE arc-length mismatch detected (max |dx|={abs(max_dx):.4e}). \"
                        \"XFOIL and RustFoil may be using different geometry/paneling.\"\
                    )
    
    for side in [1, 2]:
        side_name = "Upper" if side == 1 else "Lower"
        side_stations = {k: v for k, v in stations.items() if k[0] == side}
        
        if not side_stations:
            continue
        
        print(f"\n=== {side_name} Surface (Side {side}) ===")
        
        # Header depends on whether we have RustFoil data
        if rust_stations:
            print(f"{'IBL':>4} {'X':>10} {'XFOIL_θ':>12} {'RUST_θ':>12} {'Err%':>8} {'XFOIL_Ue':>12} {'RUST_Ue':>12} {'Err%':>8}")
        else:
            print(f"{'IBL':>4} {'X':>10} {'θ_init':>12} {'θ_final':>12} {'N_iter':>6}")
        
        for (s, ibl), e in sorted(side_stations.items()):
            x = e['x']
            theta_in = e['input']['theta']
            theta_out = e['output']['theta']
            n_iter = e['newton_iter']
            
            if rust_stations and (side, ibl) in rust_stations:
                rust = rust_stations[(side, ibl)]
                rust_theta = None
                if isinstance(rust.get('output'), dict):
                    rust_theta = rust['output'].get('theta')
                if rust_theta is None:
                    rust_theta = rust.get('theta', 0.0)
                err = abs(theta_out - rust_theta) / theta_out * 100 if theta_out != 0 else 0
                xfoil_ue = e.get('Ue', 0.0)
                rust_ue = rust.get('Ue', 0.0)
                ue_err = abs(xfoil_ue - rust_ue) / xfoil_ue * 100 if xfoil_ue != 0 else 0
                print(f"{ibl:4d} {x:10.4f} {theta_out:12.4e} {rust_theta:12.4e} {err:8.2f} {xfoil_ue:12.6f} {rust_ue:12.6f} {ue_err:8.2f}")
            else:
                print(f"{ibl:4d} {x:10.4f} {theta_in:12.4e} {theta_out:12.4e} {n_iter:6d}")
        
        # Summary stats
        thetas = [e['output']['theta'] for (s, ibl), e in side_stations.items()]
        print(f"\n  Stations: {len(thetas)}")
        print(f"  θ range: [{min(thetas):.4e}, {max(thetas):.4e}]")
        
        # Check for transition (where flow type changes)
        for (s, ibl), e in sorted(side_stations.items()):
            if e.get('converged') and e.get('dmax', 1.0) < 1e-4:
                # Found converged station
                pass


def main():
    parser = argparse.ArgumentParser(description='Compare XFOIL instrumented output with RustFoil')
    parser.add_argument('xfoil_json', help='Path to xfoil_debug.json')
    parser.add_argument('--phase', type=int, default=0,
                        help='Comparison phase (0=summary, 1=closures, 2=blvar, 3=jacobians, 4=march, 5=iteration, 6=dij, 7=mrchue, 8=trchek2, 9=trdif, 10=trdif_derivs)')
    parser.add_argument('--station', type=int, default=None,
                        help='Filter to specific station number')
    parser.add_argument('--export', type=str, default=None,
                        help='Export test vectors to directory')
    parser.add_argument('--rustfoil', type=str, default=None,
                        help='Path to RustFoil debug JSON for comparison')
    parser.add_argument('--xfoil-side', type=int, default=None,
                        help='Override XFOIL side for phase 8 (e.g., 1=upper, 2=lower)')
    parser.add_argument('--rust-side', type=int, default=None,
                        help='Override Rust side for phase 8 (e.g., 0=single-surface debug)')
    parser.add_argument('--rust-ibl', type=int, default=None,
                        help='Override Rust IBL for phase 8 when numbering differs')
    args = parser.parse_args()
    
    print(f"Loading {args.xfoil_json}...")
    events = load_xfoil_events(args.xfoil_json)
    
    summarize_events(events)
    
    if args.export:
        export_test_vectors(events, args.export)
        return
    
    if args.phase == 0:
        print("Use --phase 1-8 for detailed comparison")
        print("  1: Closure functions (HKIN, HSL, CFL, DAMPL)")
        print("  2: BLVAR secondary variables")
        print("  3: BLDIF Jacobian matrices")
        print("  4: MRCHUE initial march")
        print("  5: Newton iteration convergence")
        print("  6: QDCALC DIJ matrix")
        print("  7: MRCHUE station-by-station theta/dstar")
        print("  8: TRCHEK2 transition comparison")
        print("  9: TRDIF matrix comparison")
        print(" 10: TRDIF_DERIVS comparison")
    elif args.phase == 1:
        phase1_closures(events)
    elif args.phase == 2:
        phase2_blvar(events, args.station)
    elif args.phase == 3:
        phase3_jacobians(
            events,
            station=args.station,
            rustfoil_debug_path=args.rustfoil,
            xfoil_side=args.xfoil_side,
            rust_side=args.rust_side,
            rust_ibl=args.rust_ibl,
        )
    elif args.phase == 4:
        phase4_march(events)
    elif args.phase == 5:
        phase5_iteration(events)
    elif args.phase == 6:
        phase6_dij(events)
    elif args.phase == 7:
        phase_mrchue_stations(events, args.rustfoil)
    elif args.phase == 8:
        phase8_trchek2(
            events,
            rustfoil_debug_path=args.rustfoil,
            xfoil_side=args.xfoil_side,
            rust_side=args.rust_side,
            ibl=args.station,
            rust_ibl=args.rust_ibl,
        )
    elif args.phase == 9:
        phase9_trdif(
            events,
            rustfoil_debug_path=args.rustfoil,
            xfoil_side=args.xfoil_side,
            rust_side=args.rust_side,
            ibl=args.station,
            rust_ibl=args.rust_ibl,
        )
    elif args.phase == 10:
        phase10_trdif_derivs(
            events,
            rustfoil_debug_path=args.rustfoil,
            xfoil_side=args.xfoil_side,
            rust_side=args.rust_side,
            ibl=args.station,
            rust_ibl=args.rust_ibl,
        )


if __name__ == '__main__':
    main()
