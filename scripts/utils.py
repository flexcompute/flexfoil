"""
Shared utility functions for XFOIL-RustFoil comparison scripts.

This module provides common functionality for:
- Loading JSON debug traces from both XFOIL and RustFoil
- Extracting events by subroutine and iteration
- Calculating relative errors
- Matching stations between solvers
"""

import json
from pathlib import Path
from typing import Any, Dict, List, Optional, Tuple, Union
import re


def load_trace(filepath: Union[str, Path]) -> Dict[str, Any]:
    """
    Load a debug JSON trace file with error recovery for malformed JSON.
    
    Handles common issues with XFOIL Fortran-generated JSON:
    - Trailing commas
    - Incomplete objects
    - Non-standard number formats
    - Debug print statements mixed in with JSON
    
    Args:
        filepath: Path to JSON debug trace file
        
    Returns:
        Parsed JSON data as a dictionary
        
    Raises:
        FileNotFoundError: If file doesn't exist
        json.JSONDecodeError: If JSON is unrecoverable
    """
    filepath = Path(filepath)
    
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Try standard parsing first
    try:
        return json.loads(content)
    except json.JSONDecodeError:
        pass
    
    # Clean up XFOIL debug output issues
    
    # Remove non-JSON debug print lines (like "SETBL ENTRY", "DUE: is=...")
    lines = content.split('\n')
    cleaned_lines = []
    for line in lines:
        stripped = line.strip()
        # Skip empty lines in certain contexts
        if not stripped:
            cleaned_lines.append(line)
            continue
        # Skip lines that look like debug prints
        if stripped.startswith('DUE:'):
            continue
        if 'ENTRY' in stripped and '{' not in stripped:
            continue
        # Skip other non-JSON lines
        if stripped and not stripped.startswith(('{', '}', '[', ']', '"', ',')):
            if '=' in stripped and ':' not in stripped:
                # Likely a debug print like "var= value"
                continue
        cleaned_lines.append(line)
    content = '\n'.join(cleaned_lines)
    
    # Remove trailing commas before ] or }
    content = re.sub(r',\s*([}\]])', r'\1', content)
    
    # Remove double commas
    content = re.sub(r',\s*,', ',', content)
    
    # Fix incomplete last event (truncated file)
    if not content.rstrip().endswith('}'):
        # Find last complete event
        last_brace = content.rfind('}')
        if last_brace > 0:
            content = content[:last_brace+1]
            # Ensure we close the events array
            if '"events"' in content and not content.rstrip().endswith(']}'):
                content = content + ']}'
    
    try:
        return json.loads(content)
    except json.JSONDecodeError as e:
        # Last resort: try to extract just the events array
        match = re.search(r'"events"\s*:\s*\[', content)
        if match:
            # Find matching close bracket
            start = match.end() - 1
            depth = 0
            end = start
            for i, c in enumerate(content[start:]):
                if c == '[':
                    depth += 1
                elif c == ']':
                    depth -= 1
                    if depth == 0:
                        end = start + i + 1
                        break
            
            events_str = content[start:end]
            try:
                events = json.loads(events_str)
                return {'events': events}
            except:
                pass
        
        raise e


def extract_events(
    data: Dict[str, Any],
    subroutine: str,
    iteration: Optional[int] = None,
    side: Optional[int] = None,
) -> List[Dict[str, Any]]:
    """
    Extract events by subroutine name and optional filters.
    
    Args:
        data: Parsed JSON trace data with 'events' key
        subroutine: Name of subroutine to filter (e.g., 'BLVAR', 'SETBL')
        iteration: Optional iteration number to filter
        side: Optional surface side to filter (0=upper, 1=lower or 1=upper, 2=lower)
        
    Returns:
        List of matching events
    """
    events = data.get('events', [])
    if isinstance(data, list):
        events = data
    
    result = []
    for event in events:
        if event.get('subroutine') != subroutine:
            continue
        if iteration is not None:
            event_iter = event.get('iteration')
            if event_iter is None or event_iter != iteration:
                continue
        if side is not None:
            event_side = event.get('side') or event.get('is')
            if event_side is not None and event_side != side:
                continue
        result.append(event)
    
    return result


def calculate_rel_error(val1: float, val2: float, eps: float = 1e-10) -> float:
    """
    Calculate relative error between two values.
    
    For near-zero values, returns absolute difference instead.
    
    Args:
        val1: First value (typically XFOIL reference)
        val2: Second value (typically RustFoil)
        eps: Threshold below which to use absolute error
        
    Returns:
        Relative error as a fraction (not percentage)
    """
    if abs(val1) > eps:
        return abs(val2 - val1) / abs(val1)
    else:
        return abs(val2 - val1)


def format_rel_error(val1: float, val2: float, eps: float = 1e-10) -> str:
    """
    Format relative error as a percentage string.
    
    Args:
        val1: First value (typically XFOIL reference)
        val2: Second value (typically RustFoil)
        eps: Threshold for near-zero handling
        
    Returns:
        Formatted string like "+5.2%" or "-0.1%"
    """
    if abs(val1) > eps:
        err = (val2 - val1) / abs(val1) * 100
        return f"{err:+.1f}%"
    else:
        diff = val2 - val1
        return f"Δ={diff:.2e}"


def match_stations(
    xfoil_events: List[Dict],
    rustfoil_events: List[Dict],
    key: str = 'ibl',
) -> List[Tuple[Dict, Dict]]:
    """
    Match events from XFOIL and RustFoil by station index.
    
    Args:
        xfoil_events: List of XFOIL events
        rustfoil_events: List of RustFoil events  
        key: Key to match on (default 'ibl')
        
    Returns:
        List of (xfoil_event, rustfoil_event) tuples for matched stations
    """
    # Build lookup for RustFoil events
    rf_lookup = {}
    for event in rustfoil_events:
        idx = event.get(key)
        if idx is not None:
            rf_lookup[idx] = event
    
    # Match XFOIL events
    matched = []
    for xf_event in xfoil_events:
        idx = xf_event.get(key)
        if idx is not None and idx in rf_lookup:
            matched.append((xf_event, rf_lookup[idx]))
    
    return matched


def extract_value(event: Dict, field: str) -> Optional[float]:
    """
    Extract a numeric value from an event, handling different JSON structures.
    
    RustFoil uses nested input/output dicts, XFOIL uses flat structure.
    
    Args:
        event: Event dictionary
        field: Field name to extract
        
    Returns:
        Extracted float value or None if not found
    """
    # Direct field access
    if field in event:
        return event[field]
    
    # RustFoil nested structure
    if 'input' in event and field in event['input']:
        return event['input'][field]
    if 'output' in event and field in event['output']:
        return event['output'][field]
    
    # Try data subfield
    if 'data' in event:
        data = event['data']
        if field in data:
            return data[field]
        if 'input' in data and field in data['input']:
            return data['input'][field]
        if 'output' in data and field in data['output']:
            return data['output'][field]
    
    return None


def find_first_divergence(
    comparisons: List[Dict],
    variable: str,
    threshold: float = 0.01,
) -> Optional[Dict]:
    """
    Find the first comparison where a variable exceeds the threshold.
    
    Args:
        comparisons: List of comparison dictionaries with station data
        variable: Variable name to check (e.g., 'theta', 'dstar')
        threshold: Relative error threshold (default 1%)
        
    Returns:
        First comparison dict exceeding threshold, or None
    """
    for comp in comparisons:
        rel_key = f'{variable}_rel'
        if rel_key in comp:
            if comp[rel_key] is not None and comp[rel_key] > threshold:
                return comp
    return None


def compare_matrices(
    xf_matrix: List[List[float]],
    rf_matrix: List[List[float]],
    threshold: float = 0.05,
) -> Dict[str, Any]:
    """
    Compare two matrices element-by-element.
    
    Args:
        xf_matrix: XFOIL matrix (list of rows)
        rf_matrix: RustFoil matrix (list of rows)
        threshold: Relative error threshold
        
    Returns:
        Dictionary with comparison results:
        - max_diff: Maximum absolute difference
        - max_rel: Maximum relative error
        - first_exceed: (row, col) of first threshold exceedance
        - n_exceed: Count of elements exceeding threshold
    """
    max_diff = 0.0
    max_rel = 0.0
    first_exceed = None
    n_exceed = 0
    
    for i, (xf_row, rf_row) in enumerate(zip(xf_matrix, rf_matrix)):
        for j, (xf_val, rf_val) in enumerate(zip(xf_row, rf_row)):
            diff = abs(xf_val - rf_val)
            rel = calculate_rel_error(xf_val, rf_val)
            
            if diff > max_diff:
                max_diff = diff
            if rel > max_rel:
                max_rel = rel
            
            if rel > threshold:
                n_exceed += 1
                if first_exceed is None:
                    first_exceed = (i, j)
    
    return {
        'max_diff': max_diff,
        'max_rel': max_rel,
        'first_exceed': first_exceed,
        'n_exceed': n_exceed,
    }


def get_iteration_events(data: Dict, iteration: int) -> Dict[str, List[Dict]]:
    """
    Get all events for a specific iteration, grouped by subroutine.
    
    Args:
        data: Parsed JSON trace data
        iteration: Iteration number to extract
        
    Returns:
        Dictionary mapping subroutine names to lists of events
    """
    events = data.get('events', [])
    by_subroutine = {}
    
    for event in events:
        if event.get('iteration') == iteration:
            sub = event.get('subroutine', 'UNKNOWN')
            if sub not in by_subroutine:
                by_subroutine[sub] = []
            by_subroutine[sub].append(event)
    
    return by_subroutine


def summarize_trace(data: Dict) -> Dict[str, Any]:
    """
    Generate a summary of a debug trace file.
    
    Args:
        data: Parsed JSON trace data
        
    Returns:
        Summary dictionary with event counts and key metrics
    """
    events = data.get('events', [])
    
    # Count events by subroutine
    by_subroutine = {}
    iterations = set()
    
    for event in events:
        sub = event.get('subroutine', 'UNKNOWN')
        by_subroutine[sub] = by_subroutine.get(sub, 0) + 1
        
        it = event.get('iteration')
        if it is not None:
            iterations.add(it)
    
    # Find convergence events
    converged = False
    final_rmsbl = None
    for event in events:
        if event.get('subroutine') == 'NEWTON_CONVERGE':
            if event.get('converged'):
                converged = True
            final_rmsbl = event.get('rmsbl')
    
    return {
        'n_events': len(events),
        'n_iterations': len(iterations),
        'max_iteration': max(iterations) if iterations else 0,
        'by_subroutine': by_subroutine,
        'converged': converged,
        'final_rmsbl': final_rmsbl,
    }


# Convenience exports
__all__ = [
    'load_trace',
    'extract_events',
    'calculate_rel_error',
    'format_rel_error',
    'match_stations',
    'extract_value',
    'find_first_divergence',
    'compare_matrices',
    'get_iteration_events',
    'summarize_trace',
]
