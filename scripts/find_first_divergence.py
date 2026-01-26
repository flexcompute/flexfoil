#!/usr/bin/env python3
"""
Master Divergence Finder: Identifies the first point where RustFoil diverges from XFOIL.

This script runs comparisons in order of the computation pipeline and reports
the first stage/station where divergence exceeds the specified tolerance.

Pipeline stages:
  1. Inviscid: gamma, Cp, CL_inv (FULLGAMMA events)
  2. Setup: stagnation point, surface extraction, station counts (IBLPAN, STFIND)
  3. DIJ Matrix: influence coefficient matrix (FULLDIJ events)
  4. Initial March: station-by-station MRCHUE results
  5. Newton Iterations: SETBL, BLDIF, BLSOLV, UPDATE comparisons
  6. Forces: CD, CL, CM final results (VISCAL_RESULT)

Usage:
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json --tolerance 0.01
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json --stop-on-first --verbose
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json --stage 4
"""

import argparse
import json
import sys
from dataclasses import dataclass, field
from enum import Enum
from pathlib import Path
from typing import Any, Optional

import numpy as np


# =============================================================================
# Data Structures
# =============================================================================

class StageStatus(Enum):
    PASS = "PASS"
    FAIL = "FAIL"
    SKIP = "SKIP"  # No data available for comparison


@dataclass
class DivergenceInfo:
    """Information about a single divergence point."""
    quantity: str
    index: Optional[int] = None
    location: Optional[str] = None  # e.g., "x=0.24" or "station 38"
    xfoil_value: Any = None
    rustfoil_value: Any = None
    abs_diff: float = 0.0
    rel_diff: float = 0.0
    extra_info: dict = field(default_factory=dict)


@dataclass
class StageResult:
    """Result of comparing a single pipeline stage."""
    stage_num: int
    stage_name: str
    status: StageStatus
    summary: str = ""
    max_diff_pct: float = 0.0
    max_diff_location: str = ""
    divergences: list = field(default_factory=list)
    details: dict = field(default_factory=dict)
    recommendation: str = ""


# =============================================================================
# Utility Functions
# =============================================================================

def load_json(path: Path) -> dict:
    """Load and parse a debug JSON file."""
    with open(path) as f:
        return json.load(f)


def get_events(data: dict) -> list:
    """Extract events list from debug JSON."""
    if isinstance(data, dict) and 'events' in data:
        return data['events']
    elif isinstance(data, list):
        return data
    return []


def find_events(events: list, subroutine: str, **filters) -> list:
    """Find all events matching subroutine name and optional filters."""
    results = []
    for e in events:
        if e.get('subroutine') != subroutine:
            continue
        match = True
        for key, val in filters.items():
            if e.get(key) != val:
                match = False
                break
        if match:
            results.append(e)
    return results


def find_event(events: list, subroutine: str, **filters) -> Optional[dict]:
    """Find first event matching subroutine name and optional filters."""
    matches = find_events(events, subroutine, **filters)
    return matches[0] if matches else None


def find_last_event(events: list, subroutine: str, **filters) -> Optional[dict]:
    """Find last event matching subroutine name and optional filters."""
    matches = find_events(events, subroutine, **filters)
    return matches[-1] if matches else None


def compute_rel_diff(val1: float, val2: float, eps: float = 1e-12) -> float:
    """Compute relative difference between two values."""
    abs_diff = abs(val1 - val2)
    ref = max(abs(val1), abs(val2), eps)
    return abs_diff / ref


def compute_array_diff(arr1: np.ndarray, arr2: np.ndarray) -> dict:
    """Compute statistics for difference between two arrays."""
    if len(arr1) == 0 or len(arr2) == 0:
        return {'error': 'Empty array'}
    
    n = min(len(arr1), len(arr2))
    arr1, arr2 = arr1[:n], arr2[:n]
    
    abs_diff = np.abs(arr1 - arr2)
    
    with np.errstate(divide='ignore', invalid='ignore'):
        ref = np.maximum(np.abs(arr1), 1e-12)
        rel_diff = abs_diff / ref
        rel_diff = np.where(np.isinf(rel_diff), np.nan, rel_diff)
    
    max_idx = int(np.argmax(abs_diff))
    max_rel_idx = int(np.nanargmax(rel_diff)) if not np.all(np.isnan(rel_diff)) else max_idx
    
    return {
        'n_points': n,
        'max_abs': float(np.max(abs_diff)),
        'max_abs_idx': max_idx,
        'mean_abs': float(np.mean(abs_diff)),
        'max_rel': float(np.nanmax(rel_diff)) if not np.all(np.isnan(rel_diff)) else float('inf'),
        'max_rel_idx': max_rel_idx,
        'mean_rel': float(np.nanmean(rel_diff)) if not np.all(np.isnan(rel_diff)) else float('inf'),
    }


def find_first_exceed(arr1: np.ndarray, arr2: np.ndarray, tol: float) -> Optional[int]:
    """Find first index where relative difference exceeds tolerance."""
    if len(arr1) != len(arr2):
        return None
    
    for i in range(len(arr1)):
        rel = compute_rel_diff(arr1[i], arr2[i])
        if rel > tol:
            return i
    return None


# =============================================================================
# Stage 1: Inviscid Comparison
# =============================================================================

def compare_inviscid(xf_events: list, rf_events: list, tol: float, verbose: bool) -> StageResult:
    """Compare inviscid solver results: gamma, Cp, CL_inv."""
    result = StageResult(1, "Inviscid", StageStatus.SKIP)
    
    # Look for FULLGAMMA events
    xf_gamma = find_last_event(xf_events, 'FULLGAMMA')
    rf_gamma = find_last_event(rf_events, 'FULLGAMMA')
    
    if not xf_gamma and not rf_gamma:
        result.summary = "No FULLGAMMA events found in either file"
        return result
    
    if not xf_gamma or not rf_gamma:
        result.status = StageStatus.SKIP
        result.summary = f"FULLGAMMA events: XFOIL={'Yes' if xf_gamma else 'No'}, RustFoil={'Yes' if rf_gamma else 'No'}"
        return result
    
    result.status = StageStatus.PASS
    result.details = {'compared': []}
    
    # Compare gamma distribution
    if 'gamma' in xf_gamma and 'gamma' in rf_gamma:
        xf_g = np.array(xf_gamma['gamma'])
        rf_g = np.array(rf_gamma['gamma'])
        stats = compute_array_diff(xf_g, rf_g)
        result.details['gamma'] = stats
        result.details['compared'].append('gamma')
        
        if stats['max_rel'] > tol:
            result.status = StageStatus.FAIL
            result.divergences.append(DivergenceInfo(
                quantity='gamma',
                index=stats['max_rel_idx'],
                xfoil_value=float(xf_g[stats['max_rel_idx']]),
                rustfoil_value=float(rf_g[stats['max_rel_idx']]),
                rel_diff=stats['max_rel'],
            ))
        
        if stats['max_rel'] * 100 > result.max_diff_pct:
            result.max_diff_pct = stats['max_rel'] * 100
            result.max_diff_location = f"station {stats['max_rel_idx']}"
    
    # Compare Cp distribution
    if 'cp' in xf_gamma and 'cp' in rf_gamma:
        xf_cp = np.array(xf_gamma['cp'])
        rf_cp = np.array(rf_gamma['cp'])
        stats = compute_array_diff(xf_cp, rf_cp)
        result.details['cp'] = stats
        result.details['compared'].append('Cp')
        
        if stats['max_rel'] > tol:
            result.status = StageStatus.FAIL
            result.divergences.append(DivergenceInfo(
                quantity='Cp',
                index=stats['max_rel_idx'],
                xfoil_value=float(xf_cp[stats['max_rel_idx']]),
                rustfoil_value=float(rf_cp[stats['max_rel_idx']]),
                rel_diff=stats['max_rel'],
            ))
    
    # Compare CL_inv
    xf_cl = xf_gamma.get('cl_inv')
    rf_cl = rf_gamma.get('cl_inv')
    if xf_cl is not None and rf_cl is not None:
        rel = compute_rel_diff(xf_cl, rf_cl)
        result.details['cl_inv'] = {'xfoil': xf_cl, 'rustfoil': rf_cl, 'rel_diff': rel}
        result.details['compared'].append('CL_inv')
        
        if rel > tol:
            result.status = StageStatus.FAIL
            result.divergences.append(DivergenceInfo(
                quantity='CL_inv',
                xfoil_value=xf_cl,
                rustfoil_value=rf_cl,
                rel_diff=rel,
            ))
    
    # Build summary
    parts = []
    for name in result.details.get('compared', []):
        if name in result.details and isinstance(result.details[name], dict):
            if 'max_rel' in result.details[name]:
                pct = result.details[name]['max_rel'] * 100
                idx = result.details[name].get('max_rel_idx', '?')
                parts.append(f"{name}: {pct:.2f}% at station {idx}")
            elif 'rel_diff' in result.details[name]:
                pct = result.details[name]['rel_diff'] * 100
                parts.append(f"{name}: {pct:.2f}%")
    
    result.summary = "; ".join(parts) if parts else "Comparison complete"
    
    if result.status == StageStatus.FAIL:
        result.recommendation = "Check inviscid solver: panel influence matrix, gamma solution"
    
    return result


# =============================================================================
# Stage 2: Setup Comparison
# =============================================================================

def compare_setup(xf_events: list, rf_events: list, tol: float, verbose: bool) -> StageResult:
    """Compare setup: stagnation point, surface extraction, station counts."""
    result = StageResult(2, "Setup", StageStatus.SKIP)
    result.details = {}
    
    # Look for IBLPAN/STFIND events (stagnation point)
    xf_stag = find_event(xf_events, 'STFIND') or find_event(xf_events, 'IBLPAN')
    rf_stag = find_event(rf_events, 'STFIND') or find_event(rf_events, 'IBLPAN')
    
    # Also check VISCAL for station counts
    xf_viscal = find_event(xf_events, 'VISCAL')
    rf_viscal = find_event(rf_events, 'VISCAL')
    
    if xf_stag and rf_stag:
        result.status = StageStatus.PASS
        
        # Compare stagnation point
        xf_ist = xf_stag.get('IST') or xf_stag.get('ist')
        rf_ist = rf_stag.get('IST') or rf_stag.get('ist')
        
        if xf_ist is not None and rf_ist is not None:
            result.details['stagnation'] = {'xfoil': xf_ist, 'rustfoil': rf_ist}
            if xf_ist != rf_ist:
                result.status = StageStatus.FAIL
                result.divergences.append(DivergenceInfo(
                    quantity='stagnation_index',
                    xfoil_value=xf_ist,
                    rustfoil_value=rf_ist,
                    abs_diff=abs(xf_ist - rf_ist),
                ))
                result.recommendation = "Check stagnation point finding in setup.rs"
    
    # Count stations from BLVAR events (look at all iterations, not just iteration=1)
    # XFOIL uses side=1 for upper, side=2 for lower
    # Count unique ibl values per side across all events
    xf_blvar_events = find_events(xf_events, 'BLVAR')
    rf_blvar_events = find_events(rf_events, 'BLVAR')
    
    # Get unique (side, ibl) pairs
    xf_side1_stations = set(e.get('ibl') for e in xf_blvar_events if e.get('side') == 1 and e.get('ibl') is not None)
    xf_side2_stations = set(e.get('ibl') for e in xf_blvar_events if e.get('side') == 2 and e.get('ibl') is not None)
    rf_side1_stations = set(e.get('ibl') for e in rf_blvar_events if e.get('side') == 1 and e.get('ibl') is not None)
    rf_side2_stations = set(e.get('ibl') for e in rf_blvar_events if e.get('side') == 2 and e.get('ibl') is not None)
    
    xf_upper_stations = len(xf_side1_stations)
    xf_lower_stations = len(xf_side2_stations)
    rf_upper_stations = len(rf_side1_stations)
    rf_lower_stations = len(rf_side2_stations)
    
    # Also count total BLVAR events for comparison
    result.details['blvar_event_counts'] = {
        'xfoil': len(xf_blvar_events),
        'rustfoil': len(rf_blvar_events),
    }
    
    if xf_upper_stations > 0 or rf_upper_stations > 0 or xf_lower_stations > 0 or rf_lower_stations > 0:
        result.status = StageStatus.PASS if result.status != StageStatus.FAIL else result.status
        result.details['station_counts'] = {
            'xfoil_upper': xf_upper_stations,
            'xfoil_lower': xf_lower_stations,
            'rustfoil_upper': rf_upper_stations,
            'rustfoil_lower': rf_lower_stations,
        }
        
        # Only flag as mismatch if both have data and they differ significantly
        # Allow for some flexibility since logging may differ
        if (xf_upper_stations > 0 and rf_upper_stations > 0 and 
            abs(xf_upper_stations - rf_upper_stations) > 5):
            if result.status != StageStatus.FAIL:
                result.status = StageStatus.FAIL
            result.divergences.append(DivergenceInfo(
                quantity='station_counts',
                extra_info=result.details['station_counts'],
            ))
            result.recommendation = "Check surface extraction and paneling in setup.rs"
    
    # Build summary
    parts = []
    if 'stagnation' in result.details:
        xf_s = result.details['stagnation']['xfoil']
        rf_s = result.details['stagnation']['rustfoil']
        if xf_s == rf_s:
            parts.append(f"Stagnation: IST={xf_s} (both)")
        else:
            parts.append(f"Stagnation: XFOIL={xf_s}, RustFoil={rf_s} MISMATCH")
    
    if 'station_counts' in result.details:
        sc = result.details['station_counts']
        parts.append(f"Stations: XF upper={sc['xfoil_upper']}/lower={sc['xfoil_lower']}, "
                    f"RF upper={sc['rustfoil_upper']}/lower={sc['rustfoil_lower']}")
    
    if 'blvar_event_counts' in result.details:
        bc = result.details['blvar_event_counts']
        parts.append(f"BLVAR events: XF={bc['xfoil']}, RF={bc['rustfoil']}")
    
    result.summary = "; ".join(parts) if parts else "Limited setup data available"
    
    return result


# =============================================================================
# Stage 3: DIJ Matrix Comparison
# =============================================================================

def compare_dij(xf_events: list, rf_events: list, tol: float, verbose: bool) -> StageResult:
    """Compare DIJ influence coefficient matrix."""
    result = StageResult(3, "DIJ Matrix", StageStatus.SKIP)
    
    # Look for FULLDIJ events
    xf_dij = find_last_event(xf_events, 'FULLDIJ')
    rf_dij = find_last_event(rf_events, 'FULLDIJ')
    
    if not xf_dij and not rf_dij:
        result.summary = "No FULLDIJ events found"
        return result
    
    if not xf_dij or not rf_dij:
        result.status = StageStatus.SKIP
        result.summary = f"FULLDIJ: XFOIL={'Yes' if xf_dij else 'No'}, RustFoil={'Yes' if rf_dij else 'No'}"
        return result
    
    result.status = StageStatus.PASS
    
    # Extract DIJ matrices
    xf_nsys = xf_dij.get('nsys', 0)
    rf_nsys = rf_dij.get('nsys', 0)
    xf_mat = np.array(xf_dij.get('dij', []))
    rf_mat = np.array(rf_dij.get('dij', []))
    
    if xf_nsys > 0 and xf_mat.size == xf_nsys * xf_nsys:
        xf_mat = xf_mat.reshape((xf_nsys, xf_nsys))
    if rf_nsys > 0 and rf_mat.size == rf_nsys * rf_nsys:
        rf_mat = rf_mat.reshape((rf_nsys, rf_nsys))
    
    result.details['nsys'] = {'xfoil': xf_nsys, 'rustfoil': rf_nsys}
    
    if xf_nsys != rf_nsys:
        result.status = StageStatus.FAIL
        result.divergences.append(DivergenceInfo(
            quantity='nsys',
            xfoil_value=xf_nsys,
            rustfoil_value=rf_nsys,
        ))
        result.summary = f"DIJ size mismatch: XFOIL={xf_nsys}x{xf_nsys}, RustFoil={rf_nsys}x{rf_nsys}"
        result.recommendation = "Check station count / panel setup"
        return result
    
    if xf_mat.shape != rf_mat.shape:
        result.summary = "DIJ matrix shape mismatch"
        return result
    
    # Compare element-by-element
    diff = np.abs(xf_mat - rf_mat)
    with np.errstate(divide='ignore', invalid='ignore'):
        ref = np.maximum(np.abs(xf_mat), np.abs(rf_mat))
        ref = np.where(ref < 1e-20, 1.0, ref)
        rel_diff = diff / ref
    
    max_abs = float(np.max(diff))
    max_rel = float(np.max(rel_diff))
    max_idx = np.unravel_index(np.argmax(rel_diff), rel_diff.shape)
    
    result.details['max_abs_diff'] = max_abs
    result.details['max_rel_diff'] = max_rel
    result.details['max_diff_location'] = (int(max_idx[0]), int(max_idx[1]))
    result.max_diff_pct = max_rel * 100
    result.max_diff_location = f"({max_idx[0]}, {max_idx[1]})"
    
    # Find first element exceeding tolerance
    first_exceed = None
    for i in range(xf_mat.shape[0]):
        for j in range(xf_mat.shape[1]):
            if rel_diff[i, j] > tol:
                first_exceed = (i, j, xf_mat[i, j], rf_mat[i, j], rel_diff[i, j])
                break
        if first_exceed:
            break
    
    if first_exceed:
        result.status = StageStatus.FAIL
        i, j, xf_v, rf_v, rd = first_exceed
        result.divergences.append(DivergenceInfo(
            quantity=f'DIJ[{i},{j}]',
            index=i * xf_nsys + j,
            location=f"({i}, {j})",
            xfoil_value=xf_v,
            rustfoil_value=rf_v,
            rel_diff=rd,
        ))
        result.summary = f"First DIJ divergence at ({i}, {j}): {rd*100:.2f}%"
        result.recommendation = "Check dij.rs: influence coefficient calculation"
    else:
        result.summary = f"Max DIJ diff: {max_rel*100:.2f}% at ({max_idx[0]}, {max_idx[1]})"
    
    return result


# =============================================================================
# Stage 4: Initial March (MRCHUE) Comparison
# =============================================================================

def compare_march(xf_events: list, rf_events: list, tol: float, verbose: bool) -> StageResult:
    """Compare station-by-station march results (MRCHUE)."""
    result = StageResult(4, "Initial March (MRCHUE)", StageStatus.SKIP)
    
    # Get MRCHUE events (RustFoil) and BLVAR events (both)
    rf_mrchue = find_events(rf_events, 'MRCHUE')
    rf_blvar = find_events(rf_events, 'BLVAR')
    xf_blvar = find_events(xf_events, 'BLVAR')
    
    if not rf_mrchue and not rf_blvar:
        result.summary = "No MRCHUE or BLVAR events in RustFoil"
        return result
    
    if not xf_blvar:
        result.summary = "No BLVAR events in XFOIL"
        return result
    
    def get_x(e: dict) -> Optional[float]:
        """Extract x coordinate from event."""
        x = e.get('x')
        if x is None:
            inp = e.get('input', {})
            x = inp.get('x')
        return x
    
    # Index XFOIL BLVAR by (side, ibl) - use the first occurrence
    xf_by_key = {}
    xf_by_x_list = []  # List of (x, side, event) for proximity matching
    for e in xf_blvar:
        key = (e.get('side'), e.get('ibl'))
        if key not in xf_by_key and key[0] is not None and key[1] is not None:
            xf_by_key[key] = e
        
        x = get_x(e)
        side = e.get('side')
        if x is not None and side is not None:
            xf_by_x_list.append((x, side, e))
    
    # Index RustFoil: prefer MRCHUE for march comparison, fall back to BLVAR
    rf_by_key = {}
    rf_by_x_list = []
    rf_source = rf_mrchue if rf_mrchue else rf_blvar
    for e in rf_source:
        key = (e.get('side'), e.get('ibl'))
        if key not in rf_by_key and key[0] is not None and key[1] is not None:
            rf_by_key[key] = e
        
        x = get_x(e)
        side = e.get('side')
        if x is not None and side is not None:
            rf_by_x_list.append((x, side, e))
    
    if not xf_by_key or not rf_by_key:
        result.summary = f"Insufficient march data: XF keys={len(xf_by_key)}, RF keys={len(rf_by_key)}"
        return result
    
    result.status = StageStatus.PASS
    result.details = {'stations_compared': 0, 'stations_passed': 0, 'quantities': {}}
    result.details['xfoil_stations'] = len(xf_by_key)
    result.details['rustfoil_stations'] = len(rf_by_key)
    
    # Try matching by (side, ibl) first
    common_keys = sorted(set(xf_by_key.keys()) & set(rf_by_key.keys()))
    
    # If no common keys, try matching by x-coordinate with proximity
    matched_pairs = []
    if not common_keys and xf_by_x_list and rf_by_x_list:
        # Match by closest x-coordinate (ignoring side since conventions differ)
        # Use a tolerance of 0.001 (0.1% chord) for matching
        x_tol = 0.001
        
        # Sort both by x
        xf_sorted = sorted(xf_by_x_list, key=lambda t: t[0])
        rf_sorted = sorted(rf_by_x_list, key=lambda t: t[0])
        
        # For each RF station, find closest XF station
        used_xf_idx = set()
        for rf_x, rf_side, rf_e in rf_sorted:
            best_idx = None
            best_diff = float('inf')
            for i, (xf_x, xf_side, xf_e) in enumerate(xf_sorted):
                if i in used_xf_idx:
                    continue
                diff = abs(rf_x - xf_x)
                if diff < best_diff and diff < x_tol:
                    best_diff = diff
                    best_idx = i
            
            if best_idx is not None:
                used_xf_idx.add(best_idx)
                xf_x, xf_side, xf_e = xf_sorted[best_idx]
                matched_pairs.append((xf_e, rf_e, f"x≈{rf_x:.4f}"))
        
        result.details['matched_by'] = 'x-coordinate (proximity)'
        result.details['x_matched_count'] = len(matched_pairs)
    elif common_keys:
        for key in common_keys:
            matched_pairs.append((xf_by_key[key], rf_by_key[key], f"side={key[0]}, ibl={key[1]}"))
        result.details['matched_by'] = 'side/ibl'
    
    result.details['stations_compared'] = len(matched_pairs)
    
    if not matched_pairs:
        # Show sample x ranges to help diagnose
        xf_x_sample = sorted(set(round(t[0], 4) for t in xf_by_x_list[:20]))[:5] if xf_by_x_list else []
        rf_x_sample = sorted(set(round(t[0], 4) for t in rf_by_x_list[:20]))[:5] if rf_by_x_list else []
        result.summary = (f"No common stations found "
                         f"(XF x={xf_x_sample}, RF x={rf_x_sample})")
        return result
    
    first_diverge = None
    max_diff_station = None
    max_diff_value = 0.0
    
    # Define paths for extracting values from both XFOIL BLVAR and RustFoil MRCHUE/BLVAR
    quantities_to_compare = [
        ('theta', ['input.theta', 'theta']),
        ('Hk', ['output.Hk', 'Hk']),
        ('Cf', ['output.Cf', 'Cf']),
        ('delta_star', ['input.delta_star', 'delta_star']),
        ('ctau', ['input.ctau', 'ctau']),
        ('Ue', ['input.u', 'Ue', 'U']),
    ]
    
    def get_value(d: dict, paths: list):
        """Get value trying multiple paths."""
        for path in paths:
            parts = path.split('.')
            val = d
            for p in parts:
                if isinstance(val, dict):
                    val = val.get(p)
                else:
                    val = None
                    break
            if val is not None:
                return val
        return None
    
    for xf_e, rf_e, location_str in matched_pairs:
        station_diverged = False
        
        for name, paths in quantities_to_compare:
            xf_val = get_value(xf_e, paths)
            rf_val = get_value(rf_e, paths)
            
            if xf_val is None or rf_val is None:
                continue
            
            rel = compute_rel_diff(xf_val, rf_val)
            
            if rel > max_diff_value:
                max_diff_value = rel
                max_diff_station = (location_str, name, xf_val, rf_val, rel)
            
            if rel > tol:
                station_diverged = True
                if first_diverge is None:
                    x_val = get_value(xf_e, ['input.x', 'x'])
                    first_diverge = {
                        'location': location_str,
                        'quantity': name,
                        'xfoil': xf_val,
                        'rustfoil': rf_val,
                        'rel_diff': rel,
                        'x': x_val,
                    }
        
        if not station_diverged:
            result.details['stations_passed'] += 1
    
    if first_diverge:
        result.status = StageStatus.FAIL
        result.divergences.append(DivergenceInfo(
            quantity=first_diverge['quantity'],
            location=first_diverge['location'],
            xfoil_value=first_diverge['xfoil'],
            rustfoil_value=first_diverge['rustfoil'],
            rel_diff=first_diverge['rel_diff'],
        ))
        result.summary = (f"Station {first_diverge['location']}: "
                         f"{first_diverge['quantity']} diff {first_diverge['rel_diff']*100:.2f}%")
        result.max_diff_pct = first_diverge['rel_diff'] * 100
        result.max_diff_location = first_diverge['location']
        
        # Provide specific recommendation
        q = first_diverge['quantity']
        if q in ('theta', 'delta_star'):
            result.recommendation = "Check momentum/displacement thickness equations in march.rs"
        elif q == 'Hk':
            result.recommendation = "Check shape parameter calculation in equations.rs (HKIN)"
        elif q == 'Cf':
            result.recommendation = "Check skin friction correlation in equations.rs (CFL)"
        elif q == 'ctau':
            result.recommendation = "Check turbulent shear stress in march.rs, transition handling"
        else:
            result.recommendation = f"Check {q} calculation in march.rs"
    else:
        match_type = result.details.get('matched_by', 'key')
        result.summary = (f"Compared {result.details['stations_compared']} stations (by {match_type}); "
                         f"{result.details['stations_passed']} passed; "
                         f"max diff: {max_diff_value*100:.4f}%")
        if max_diff_station:
            loc, name, _, _, rel = max_diff_station
            result.max_diff_pct = rel * 100
            result.max_diff_location = f"{loc}, {name}"
    
    return result


# =============================================================================
# Stage 5: Newton Iterations Comparison
# =============================================================================

def compare_newton(xf_events: list, rf_events: list, tol: float, verbose: bool) -> StageResult:
    """Compare Newton iteration: SETBL, BLDIF, BLSOLV, UPDATE."""
    result = StageResult(5, "Newton Iterations", StageStatus.SKIP)
    
    # Check for SETBL events (XFOIL) and BLDIF events (both)
    xf_setbl = find_events(xf_events, 'SETBL')
    xf_bldif = find_events(xf_events, 'BLDIF')
    rf_bldif = find_events(rf_events, 'BLDIF')
    
    # Also check BLSOLV and UPDATE
    xf_blsolv = find_events(xf_events, 'BLSOLV')
    xf_update = find_events(xf_events, 'UPDATE')
    
    result.details = {
        'xfoil_bldif_count': len(xf_bldif),
        'rustfoil_bldif_count': len(rf_bldif),
        'xfoil_setbl_count': len(xf_setbl),
        'xfoil_blsolv_count': len(xf_blsolv),
        'xfoil_update_count': len(xf_update),
    }
    
    if not xf_bldif or not rf_bldif:
        result.summary = f"Insufficient BLDIF events: XF={len(xf_bldif)}, RF={len(rf_bldif)}"
        return result
    
    result.status = StageStatus.PASS
    
    # Compare BLDIF events (Jacobians and residuals)
    # Index by (side, ibl) only - iteration may not match between implementations
    xf_by_key = {}
    for e in xf_bldif:
        key = (e.get('side'), e.get('ibl'))
        if key not in xf_by_key and key[0] is not None and key[1] is not None:
            xf_by_key[key] = e
    
    rf_by_key = {}
    for e in rf_bldif:
        key = (e.get('side'), e.get('ibl'))
        if key not in rf_by_key and key[0] is not None and key[1] is not None:
            rf_by_key[key] = e
    
    common_keys = sorted(set(xf_by_key.keys()) & set(rf_by_key.keys()))
    result.details['bldif_compared'] = len(common_keys)
    result.details['xfoil_unique_stations'] = len(xf_by_key)
    result.details['rustfoil_unique_stations'] = len(rf_by_key)
    
    if not common_keys:
        result.summary = (f"No common BLDIF stations found "
                         f"(XF: {list(xf_by_key.keys())[:3]}..., RF: {list(rf_by_key.keys())[:3]}...)")
        return result
    
    first_diverge = None
    max_vsrez_diff = 0.0
    max_vsrez_location = None
    compared_count = 0
    
    for key in common_keys:
        xf_e = xf_by_key[key]
        rf_e = rf_by_key[key]
        
        side, ibl = key
        
        # Compare VSREZ (residuals)
        xf_vsrez = np.array(xf_e.get('VSREZ', []))
        rf_vsrez = np.array(rf_e.get('VSREZ', []))
        
        if len(xf_vsrez) > 0 and len(rf_vsrez) > 0:
            compared_count += 1
            n = min(len(xf_vsrez), len(rf_vsrez))
            diff = np.abs(xf_vsrez[:n] - rf_vsrez[:n])
            max_diff = float(np.max(diff))
            
            # Use absolute diff for residuals (they're already small)
            if max_diff > max_vsrez_diff:
                max_vsrez_diff = max_diff
                max_vsrez_location = key
            
            # Check relative for non-tiny values
            for i in range(n):
                if abs(xf_vsrez[i]) > 1e-10:
                    rel = compute_rel_diff(xf_vsrez[i], rf_vsrez[i])
                    if rel > tol and first_diverge is None:
                        first_diverge = {
                            'side': side,
                            'ibl': ibl,
                            'component': f'VSREZ[{i}]',
                            'xfoil': float(xf_vsrez[i]),
                            'rustfoil': float(rf_vsrez[i]),
                            'rel_diff': rel,
                        }
        
        # Compare VS1, VS2 Jacobians
        for vs_name in ['VS1', 'VS2']:
            xf_vs = np.array(xf_e.get(vs_name, []))
            rf_vs = np.array(rf_e.get(vs_name, []))
            
            if xf_vs.size > 0 and rf_vs.size > 0:
                xf_vs = xf_vs.flatten()
                rf_vs = rf_vs.flatten()
                n = min(len(xf_vs), len(rf_vs))
                
                for i in range(n):
                    if abs(xf_vs[i]) > 1e-10:
                        rel = compute_rel_diff(xf_vs[i], rf_vs[i])
                        if rel > tol and first_diverge is None:
                            first_diverge = {
                                'side': side,
                                'ibl': ibl,
                                'component': f'{vs_name}[{i}]',
                                'xfoil': float(xf_vs[i]),
                                'rustfoil': float(rf_vs[i]),
                                'rel_diff': rel,
                            }
    
    if first_diverge:
        result.status = StageStatus.FAIL
        result.divergences.append(DivergenceInfo(
            quantity=first_diverge['component'],
            index=first_diverge['ibl'],
            location=f"side={first_diverge['side']}, ibl={first_diverge['ibl']}",
            xfoil_value=first_diverge['xfoil'],
            rustfoil_value=first_diverge['rustfoil'],
            rel_diff=first_diverge['rel_diff'],
        ))
        result.summary = (f"Station side={first_diverge['side']}, ibl={first_diverge['ibl']}: "
                         f"{first_diverge['component']} diff {first_diverge['rel_diff']*100:.2f}%")
        result.max_diff_pct = first_diverge['rel_diff'] * 100
        result.max_diff_location = f"side={first_diverge['side']}, ibl={first_diverge['ibl']}"
        result.recommendation = "Check BLDIF equation setup in newton.rs"
    else:
        result.summary = f"Compared {compared_count} BLDIF events at {len(common_keys)} stations; max VSREZ diff: {max_vsrez_diff:.2e}"
        if max_vsrez_location:
            result.max_diff_location = f"side={max_vsrez_location[0]}, ibl={max_vsrez_location[1]}"
    
    return result


# =============================================================================
# Stage 6: Forces Comparison
# =============================================================================

def compare_forces(xf_events: list, rf_events: list, tol: float, verbose: bool) -> StageResult:
    """Compare final force results: CD, CL, CM."""
    result = StageResult(6, "Forces", StageStatus.SKIP)
    
    # Look for VISCAL_RESULT or VISCOUS_FINAL
    xf_result = find_last_event(xf_events, 'VISCAL_RESULT') or find_last_event(xf_events, 'VISCOUS_FINAL')
    rf_result = find_last_event(rf_events, 'VISCAL_RESULT') or find_last_event(rf_events, 'VISCOUS_FINAL')
    
    if not xf_result and not rf_result:
        result.summary = "No VISCAL_RESULT events found"
        return result
    
    if not xf_result or not rf_result:
        result.status = StageStatus.SKIP
        result.summary = f"VISCAL_RESULT: XFOIL={'Yes' if xf_result else 'No'}, RustFoil={'Yes' if rf_result else 'No'}"
        return result
    
    result.status = StageStatus.PASS
    result.details = {}
    
    force_quantities = [
        ('CL', 'CL', 'Lift coefficient'),
        ('CD', 'CD', 'Drag coefficient'),
        ('CDf', 'CDf', 'Friction drag'),
        ('CDp', 'CDp', 'Pressure drag'),
        ('CM', 'CM', 'Pitching moment'),
    ]
    
    first_diverge = None
    
    for key, alt_key, desc in force_quantities:
        xf_val = xf_result.get(key) or xf_result.get(alt_key)
        rf_val = rf_result.get(key) or rf_result.get(alt_key)
        
        if xf_val is not None and rf_val is not None:
            rel = compute_rel_diff(xf_val, rf_val)
            result.details[key] = {
                'xfoil': xf_val,
                'rustfoil': rf_val,
                'rel_diff': rel,
            }
            
            if rel > tol and first_diverge is None:
                first_diverge = {
                    'quantity': key,
                    'description': desc,
                    'xfoil': xf_val,
                    'rustfoil': rf_val,
                    'rel_diff': rel,
                }
            
            if rel * 100 > result.max_diff_pct:
                result.max_diff_pct = rel * 100
                result.max_diff_location = key
    
    if first_diverge:
        result.status = StageStatus.FAIL
        result.divergences.append(DivergenceInfo(
            quantity=first_diverge['quantity'],
            xfoil_value=first_diverge['xfoil'],
            rustfoil_value=first_diverge['rustfoil'],
            rel_diff=first_diverge['rel_diff'],
        ))
        result.summary = (f"{first_diverge['quantity']} ({first_diverge['description']}): "
                         f"XFOIL={first_diverge['xfoil']:.6f}, "
                         f"RustFoil={first_diverge['rustfoil']:.6f}, "
                         f"diff={first_diverge['rel_diff']*100:.2f}%")
        result.recommendation = "Check force integration in forces.rs"
    else:
        parts = []
        for key in ['CL', 'CD']:
            if key in result.details:
                d = result.details[key]
                parts.append(f"{key}: {d['rel_diff']*100:.2f}%")
        result.summary = "; ".join(parts) if parts else "Force comparison complete"
    
    return result


# =============================================================================
# Main Comparison Runner
# =============================================================================

def run_all_stages(
    xf_events: list,
    rf_events: list,
    tolerance: float,
    stop_on_first: bool,
    verbose: bool,
    stage_filter: Optional[int] = None,
) -> list[StageResult]:
    """Run all comparison stages and return results."""
    
    stages = [
        (1, compare_inviscid),
        (2, compare_setup),
        (3, compare_dij),
        (4, compare_march),
        (5, compare_newton),
        (6, compare_forces),
    ]
    
    results = []
    first_fail = None
    
    for stage_num, compare_fn in stages:
        if stage_filter is not None and stage_num != stage_filter:
            continue
        
        result = compare_fn(xf_events, rf_events, tolerance, verbose)
        results.append(result)
        
        if result.status == StageStatus.FAIL and first_fail is None:
            first_fail = result
            if stop_on_first:
                break
    
    return results


def print_report(results: list[StageResult], tolerance: float, verbose: bool):
    """Print the divergence analysis report."""
    print("\n" + "=" * 70)
    print(" DIVERGENCE ANALYSIS REPORT")
    print("=" * 70)
    print(f" Tolerance: {tolerance*100:.2f}%")
    print()
    
    first_fail = None
    
    for r in results:
        status_str = {
            StageStatus.PASS: "PASS",
            StageStatus.FAIL: "FAIL",
            StageStatus.SKIP: "SKIP",
        }[r.status]
        
        marker = ""
        if r.status == StageStatus.FAIL and first_fail is None:
            first_fail = r
            marker = " ← FIRST DIVERGENCE"
        
        print(f"Stage {r.stage_num} ({r.stage_name}): {status_str}{marker}")
        
        if r.summary:
            # Indent summary lines
            for line in r.summary.split('\n'):
                print(f"  {line}")
        
        if verbose and r.details:
            for key, val in r.details.items():
                if isinstance(val, dict) and 'max_rel' in val:
                    print(f"    {key}: max_rel={val['max_rel']*100:.4f}% at idx {val.get('max_rel_idx', '?')}")
        
        if r.status == StageStatus.FAIL and r.divergences:
            print()
            for div in r.divergences:
                loc_str = f" ({div.location})" if div.location else ""
                print(f"    {div.quantity}{loc_str}:")
                print(f"      XFOIL:    {div.xfoil_value}")
                print(f"      RustFoil: {div.rustfoil_value}")
                if div.rel_diff > 0:
                    print(f"      Diff:     {div.rel_diff*100:.4f}%")
        
        print()
    
    print("-" * 70)
    
    if first_fail:
        print(f"RESULT: DIVERGENCE FOUND at Stage {first_fail.stage_num} ({first_fail.stage_name})")
        if first_fail.recommendation:
            print(f"\nRecommendation: {first_fail.recommendation}")
        return 1
    else:
        all_pass = all(r.status in (StageStatus.PASS, StageStatus.SKIP) for r in results)
        if all_pass:
            print("RESULT: ALL STAGES PASSED")
            return 0
        else:
            print("RESULT: Some stages skipped (missing data)")
            return 0


def main():
    parser = argparse.ArgumentParser(
        description="Master divergence finder: identifies first point where RustFoil diverges from XFOIL",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json --tolerance 0.01
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json --stop-on-first --verbose
    python scripts/find_first_divergence.py xfoil_debug.json rustfoil_debug.json --stage 4

Stages:
    1: Inviscid (gamma, Cp, CL_inv)
    2: Setup (stagnation, station counts)
    3: DIJ Matrix (influence coefficients)
    4: Initial March (MRCHUE station-by-station)
    5: Newton Iterations (SETBL, BLDIF, BLSOLV)
    6: Forces (CD, CL, CM)
        """
    )
    parser.add_argument("xfoil_json", type=Path, help="Path to XFOIL debug JSON file")
    parser.add_argument("rustfoil_json", type=Path, help="Path to RustFoil debug JSON file")
    parser.add_argument(
        "--tolerance", "-t",
        type=float,
        default=0.01,
        help="Relative difference threshold for divergence (default: 0.01 = 1%%)"
    )
    parser.add_argument(
        "--stop-on-first", "-s",
        action="store_true",
        help="Stop at first divergence found"
    )
    parser.add_argument(
        "--verbose", "-v",
        action="store_true",
        help="Show detailed output per stage"
    )
    parser.add_argument(
        "--stage",
        type=int,
        choices=[1, 2, 3, 4, 5, 6],
        help="Run only specific stage (1-6)"
    )
    
    args = parser.parse_args()
    
    # Validate input files
    if not args.xfoil_json.exists():
        print(f"Error: XFOIL JSON file not found: {args.xfoil_json}")
        return 1
    if not args.rustfoil_json.exists():
        print(f"Error: RustFoil JSON file not found: {args.rustfoil_json}")
        return 1
    
    # Load data
    print(f"Loading XFOIL debug: {args.xfoil_json}")
    xf_data = load_json(args.xfoil_json)
    xf_events = get_events(xf_data)
    print(f"  {len(xf_events)} events")
    
    print(f"Loading RustFoil debug: {args.rustfoil_json}")
    rf_data = load_json(args.rustfoil_json)
    rf_events = get_events(rf_data)
    print(f"  {len(rf_events)} events")
    
    # Run comparisons
    results = run_all_stages(
        xf_events,
        rf_events,
        args.tolerance,
        args.stop_on_first,
        args.verbose,
        args.stage,
    )
    
    # Print report and return exit code
    return print_report(results, args.tolerance, args.verbose)


if __name__ == "__main__":
    sys.exit(main())
