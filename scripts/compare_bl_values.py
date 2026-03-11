#!/usr/bin/env python3
"""Compare BL intermediate values between XFOIL and RustFoil."""

import json
import re
import sys

def parse_xfoil_debug(filepath):
    """Parse XFOIL debug JSON, handling potential format issues."""
    with open(filepath, 'r') as f:
        content = f.read()
    
    # Fix common JSON issues
    # Remove trailing commas before }
    content = re.sub(r',\s*}', '}', content)
    content = re.sub(r',\s*]', ']', content)
    
    # Remove non-JSON lines like "SETBL ENTRY"
    lines = content.split('\n')
    cleaned_lines = []
    for line in lines:
        if 'ENTRY' not in line and 'DUE:' not in line:
            cleaned_lines.append(line)
    content = '\n'.join(cleaned_lines)
    
    try:
        data = json.loads(content)
        return data.get('events', data) if isinstance(data, dict) else data
    except json.JSONDecodeError as e:
        print(f"JSON parse error at position {e.pos}: {e.msg}")
        # Try to find the error location
        start = max(0, e.pos - 100)
        end = min(len(content), e.pos + 100)
        print(f"Context: ...{content[start:end]}...")
        return []

def main():
    xfoil_path = "/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/bin/xfoil_debug.json"
    
    events = parse_xfoil_debug(xfoil_path)
    if not events:
        print("Failed to parse XFOIL debug file")
        return
    
    # Find converged iteration
    iterations = set()
    for e in events:
        if isinstance(e, dict) and 'iteration' in e:
            iterations.add(e['iteration'])
    
    max_iter = max(iterations) if iterations else 0
    print(f"XFOIL converged at iteration {max_iter}")
    
    # Get BL_STATE at final iteration
    bl_states = [e for e in events 
                 if isinstance(e, dict) 
                 and e.get('subroutine') == 'BL_STATE' 
                 and e.get('iteration') == max_iter]
    
    print(f"\n=== XFOIL BL State at Iteration {max_iter} ===")
    print("IBL | Side |    x/c     |    theta     |    dstar     |     Ue      |  mass")
    print("-" * 85)
    
    for s in sorted(bl_states, key=lambda x: (x.get('is',0), x.get('ibl',0))):
        ibl = s.get('ibl', 0)
        side = s.get('is', 0)
        xssi = s.get('xssi', 0)
        thet = s.get('thet', 0)
        dstr = s.get('dstr', 0)
        uedg = s.get('uedg', 0)
        mass = s.get('mass', 0)
        print(f" {ibl:3d} |  {side}   | {xssi:10.6f} | {thet:12.6e} | {dstr:12.6e} | {uedg:10.6f} | {mass:12.6e}")
    
    # Get UESET data
    ueset = [e for e in events 
             if isinstance(e, dict) 
             and e.get('subroutine') == 'UESET' 
             and e.get('iteration') == max_iter]
    
    if ueset:
        u = ueset[-1]
        upper = u.get('upper_surface', {})
        lower = u.get('lower_surface', {})
        
        print("\n=== XFOIL Upper Surface Ue (first 10 stations) ===")
        ue_upper = upper.get('ue_after', [])[:10]
        for i, ue in enumerate(ue_upper):
            print(f"  Station {i+2}: Ue = {ue:.6f}")
        
        print("\n=== XFOIL Upper Surface theta (first 10 stations) ===")
        theta_upper = upper.get('theta', [])[:10] if 'theta' in upper else []
        for i, th in enumerate(theta_upper):
            print(f"  Station {i+2}: theta = {th:.6e}")

if __name__ == '__main__':
    main()
