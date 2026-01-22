# RustFoil Solver Debugging Plan

## Objective

Step through the viscous solver comparing XFOIL vs RustFoil outputs at each stage to identify where values diverge.

## Prerequisites

We have:
- XFOIL instrumented output: `Xfoil-instrumented/bin/xfoil_debug.json` (~25,000 events)
- RustFoil debug module: `rustfoil_bl::debug` (same event format)
- Events are ordered by `call_id` (monotonically increasing)

## Phase 1: Generate RustFoil Debug Output

### Step 1.1: Add Debug Initialization to CLI

Modify `rustfoil-cli` to enable debug output:

```rust
// In run_viscous() add:
rustfoil_bl::init_debug("rustfoil_debug.json");

// ... run solver ...

rustfoil_bl::finalize_debug();
```

### Step 1.2: Instrument Missing Events

Currently instrumented:
- [x] VISCAL start
- [x] VISCAL_RESULT per iteration
- [x] BLVAR (via blvar_debug)
- [x] BLDIF (via bldif_debug)
- [x] MRCHUE (via march_fixed_ue_debug)
- [x] QDCALC (via build_dij_matrix_debug)

**Need to add:**
- [ ] UPDATE - Newton deltas applied at each station
- [ ] BLSOLV - Newton system solve

### Step 1.3: Run RustFoil with Debug

```bash
cargo run --release -p rustfoil-cli -- viscous testdata/naca0012_repaneled.dat \
    --alpha 4.0 --re 3e6 --debug rustfoil_debug.json
```

---

## Phase 2: Side-by-Side Comparison Script

Create `scripts/debug_compare.py`:

```python
#!/usr/bin/env python3
"""
Side-by-side comparison of XFOIL and RustFoil debug events.

Aligns events by type and station, highlights divergences.
"""

import json
import argparse
from collections import defaultdict
from pathlib import Path

def load_events(path):
    with open(path) as f:
        data = json.load(f)
    return data.get('events', data)

def group_by_type(events):
    """Group events by subroutine name."""
    groups = defaultdict(list)
    for ev in events:
        groups[ev['subroutine']].append(ev)
    return groups

def compare_mrchue(xfoil_events, rust_events, tolerance=0.01):
    """Compare initial march profiles."""
    print("\n" + "="*70)
    print("MRCHUE (Initial March) Comparison")
    print("="*70)
    
    # Group by station
    xf_by_station = {(e['side'], e['ibl']): e for e in xfoil_events}
    rf_by_station = {(e.get('side', 1), e['ibl']): e for e in rust_events}
    
    all_stations = sorted(set(xf_by_station.keys()) | set(rf_by_station.keys()))
    
    print(f"{'Station':<10} {'Field':<10} {'XFOIL':>14} {'RustFoil':>14} {'Diff%':>10} {'Status':<6}")
    print("-" * 70)
    
    errors = []
    for station in all_stations[:20]:  # First 20 stations
        xf = xf_by_station.get(station, {})
        rf = rf_by_station.get(station, {})
        
        if not xf or not rf:
            print(f"{station!s:<10} {'MISSING':<10} {'---':>14} {'---':>14}")
            continue
            
        for field in ['theta', 'delta_star', 'Hk', 'Cf', 'Ue']:
            xf_val = xf.get(field, 0)
            rf_val = rf.get(field, 0)
            
            if abs(xf_val) > 1e-15:
                diff_pct = abs(xf_val - rf_val) / abs(xf_val) * 100
            else:
                diff_pct = 0 if abs(rf_val) < 1e-15 else 100
                
            status = "✓" if diff_pct < tolerance * 100 else "✗"
            if diff_pct >= tolerance * 100:
                errors.append((station, field, xf_val, rf_val, diff_pct))
            
            print(f"{station!s:<10} {field:<10} {xf_val:>14.6e} {rf_val:>14.6e} {diff_pct:>9.2f}% {status:<6}")
    
    return errors

def compare_blvar(xfoil_events, rust_events, iteration=1, tolerance=0.01):
    """Compare BLVAR secondary variables for a specific iteration."""
    print(f"\n" + "="*70)
    print(f"BLVAR (Secondary Variables) - Iteration {iteration}")
    print("="*70)
    
    # Filter by iteration
    xf_iter = [e for e in xfoil_events if e.get('iteration') == iteration]
    rf_iter = [e for e in rust_events if e.get('iteration') == iteration]
    
    # Group by station
    xf_by_station = {(e['side'], e['ibl']): e for e in xf_iter}
    rf_by_station = {(e.get('side', 1), e['ibl']): e for e in rf_iter}
    
    errors = []
    fields = ['H', 'Hk', 'Hs', 'Cf', 'Cd', 'Rtheta']
    
    for station in sorted(xf_by_station.keys())[:10]:
        xf = xf_by_station.get(station, {})
        rf = rf_by_station.get(station, {})
        
        if not rf:
            print(f"Station {station}: MISSING in RustFoil")
            continue
            
        xf_out = xf.get('output', {})
        rf_out = rf.get('output', {})
        
        print(f"\nStation {station}:")
        for field in fields:
            xf_val = xf_out.get(field, 0)
            rf_val = rf_out.get(field, 0)
            
            if abs(xf_val) > 1e-15:
                diff_pct = abs(xf_val - rf_val) / abs(xf_val) * 100
            else:
                diff_pct = 0
                
            status = "✓" if diff_pct < tolerance * 100 else "✗"
            print(f"  {field:<8} XFOIL={xf_val:>12.6e}  Rust={rf_val:>12.6e}  diff={diff_pct:>6.2f}% {status}")
            
            if diff_pct >= tolerance * 100:
                errors.append((station, field, xf_val, rf_val, diff_pct))
    
    return errors

def compare_update(xfoil_events, rust_events, iteration=1, tolerance=0.10):
    """Compare UPDATE deltas - this is where divergence likely starts."""
    print(f"\n" + "="*70)
    print(f"UPDATE (Newton Deltas) - Iteration {iteration}")
    print("="*70)
    
    xf_iter = [e for e in xfoil_events if e.get('iteration') == iteration]
    rf_iter = [e for e in rust_events if e.get('iteration') == iteration]
    
    print(f"XFOIL has {len(xf_iter)} UPDATE events, RustFoil has {len(rf_iter)}")
    
    if not xf_iter or not rf_iter:
        print("Cannot compare - missing events")
        return []
    
    # Compare first few stations
    for i, (xf, rf) in enumerate(zip(xf_iter[:10], rf_iter[:10])):
        station = (xf.get('side', 1), xf.get('ibl', i))
        print(f"\nStation {station}:")
        
        for field in ['delta_ctau', 'delta_theta', 'delta_mass', 'delta_Ue', 'relaxation']:
            xf_val = xf.get(field, 0)
            rf_val = rf.get(field, 0)
            
            print(f"  {field:<12} XFOIL={xf_val:>12.6e}  Rust={rf_val:>12.6e}")
    
    return []

def compare_convergence(xfoil_events, rust_events):
    """Compare iteration convergence history."""
    print(f"\n" + "="*70)
    print("VISCAL_RESULT (Convergence History)")
    print("="*70)
    
    print(f"{'Iter':<6} {'XFOIL RMS':>12} {'Rust RMS':>12} {'XFOIL CL':>10} {'Rust CL':>10} {'XFOIL CD':>12} {'Rust CD':>12}")
    print("-" * 80)
    
    for xf in xfoil_events:
        it = xf.get('iteration', 0)
        rf = next((e for e in rust_events if e.get('iteration') == it), None)
        
        xf_rms = xf.get('rms_residual', 0)
        xf_cl = xf.get('CL', 0)
        xf_cd = xf.get('CD', 0)
        
        rf_rms = rf.get('rms_residual', 0) if rf else 0
        rf_cl = rf.get('CL', 0) if rf else 0
        rf_cd = rf.get('CD', 0) if rf else 0
        
        print(f"{it:<6} {xf_rms:>12.4e} {rf_rms:>12.4e} {xf_cl:>10.4f} {rf_cl:>10.4f} {xf_cd:>12.6e} {rf_cd:>12.6e}")

def main():
    parser = argparse.ArgumentParser(description="Compare XFOIL and RustFoil debug output")
    parser.add_argument('xfoil_json', help='Path to XFOIL debug JSON')
    parser.add_argument('rust_json', help='Path to RustFoil debug JSON')
    parser.add_argument('--stage', choices=['mrchue', 'blvar', 'update', 'convergence', 'all'],
                       default='all', help='Which stage to compare')
    parser.add_argument('--iteration', type=int, default=1, help='Iteration to compare')
    args = parser.parse_args()
    
    print(f"Loading XFOIL: {args.xfoil_json}")
    xf_events = load_events(args.xfoil_json)
    print(f"  {len(xf_events)} events")
    
    print(f"Loading RustFoil: {args.rust_json}")
    rf_events = load_events(args.rust_json)
    print(f"  {len(rf_events)} events")
    
    xf_groups = group_by_type(xf_events)
    rf_groups = group_by_type(rf_events)
    
    print("\nEvent counts:")
    all_types = set(xf_groups.keys()) | set(rf_groups.keys())
    for t in sorted(all_types):
        print(f"  {t:<15} XFOIL={len(xf_groups.get(t, [])): >5}  Rust={len(rf_groups.get(t, [])): >5}")
    
    if args.stage in ['mrchue', 'all']:
        compare_mrchue(xf_groups.get('MRCHUE', []), rf_groups.get('MRCHUE', []))
    
    if args.stage in ['blvar', 'all']:
        compare_blvar(xf_groups.get('BLVAR', []), rf_groups.get('BLVAR', []), args.iteration)
    
    if args.stage in ['update', 'all']:
        compare_update(xf_groups.get('UPDATE', []), rf_groups.get('UPDATE', []), args.iteration)
    
    if args.stage in ['convergence', 'all']:
        compare_convergence(xf_groups.get('VISCAL_RESULT', []), rf_groups.get('VISCAL_RESULT', []))

if __name__ == '__main__':
    main()
```

---

## Phase 3: Systematic Debugging Sequence

### Stage 1: Initial March (MRCHUE)

**First divergence point to check.** If march outputs differ, all subsequent calculations will be wrong.

```bash
python scripts/debug_compare.py \
    Xfoil-instrumented/bin/xfoil_debug.json \
    rustfoil_debug.json \
    --stage mrchue
```

**Expected:** θ, δ*, Hk should match within 1%.

**If they don't match:**
- Check `march_fixed_ue()` implementation
- Verify stagnation point initialization
- Compare arc length computation

### Stage 2: Secondary Variables (BLVAR) - Iteration 1

**Check BLVAR at first Newton iteration.**

```bash
python scripts/debug_compare.py ... --stage blvar --iteration 1
```

**Expected:** H, Hk, Hs, Cf, Cd should match within 1%.

**If they don't match:**
- The closure functions are verified correct (unit tests pass)
- Issue is likely in how `blvar()` is called or inputs differ
- Check: Are Ue, θ, δ* inputs identical?

### Stage 3: Newton Deltas (UPDATE) - Iteration 1

**Critical check.** This is where XFOIL's UPDATE applies Newton corrections.

```bash
python scripts/debug_compare.py ... --stage update --iteration 1
```

**Fields to check:**
- `delta_ctau` - Change in shear stress coefficient
- `delta_theta` - Change in momentum thickness
- `delta_mass` - Change in mass defect (δ* × Ue)
- `delta_Ue` - Change in edge velocity
- `relaxation` - Under-relaxation factor

**If deltas differ:**
- Issue is in Newton system (build or solve)
- Check `BlNewtonSystem::build()`
- Check `solve_bl_system()` output

### Stage 4: Edge Velocity Update (UESET)

**Check if Ue is updated correctly from mass defect changes.**

After UPDATE, XFOIL calls UESET which computes:
```
Ue_new = Ue_inviscid + DIJ × Δ(mass_defect)
```

**If Ue doesn't update correctly:**
- Check `set_edge_velocities()` implementation
- Verify DIJ matrix diagonal/off-diagonal structure

### Stage 5: Convergence History

```bash
python scripts/debug_compare.py ... --stage convergence
```

**Expected XFOIL pattern:**
```
Iter    RMS         CL        CD
1    9.5e-02    0.465    0.00603
2    4.2e-02    0.450    0.00623
3    1.2e-02    0.443    0.00618  ← Converged
```

---

## Phase 4: Deep Dive Debugging

Once the divergence point is identified, add more detailed logging.

### 4.1: Newton System Matrices

Add debug output for VA, VB blocks:

```rust
// In newton.rs build()
if rustfoil_bl::is_debug_active() {
    rustfoil_bl::add_event(DebugEvent::newton_block(
        iteration, i, 
        system.va[i], 
        system.vb[i],
        system.rhs[i]
    ));
}
```

### 4.2: Matrix Solution Check

After `solve_bl_system()`, verify residual:
```rust
// Compute A*x - b residual
for i in 1..n {
    let va_x = multiply_3x3_vec(&system.va[i], &solution[i]);
    let vb_xprev = if i > 1 { multiply_3x3_vec(&system.vb[i], &solution[i-1]) } else { [0.0; 3] };
    let residual = [
        va_x[0] + vb_xprev[0] - system.rhs[i][0],
        va_x[1] + vb_xprev[1] - system.rhs[i][1],
        va_x[2] + vb_xprev[2] - system.rhs[i][2],
    ];
    // Log if |residual| > 1e-10
}
```

### 4.3: DIJ Matrix Verification

Compare DIJ diagonal terms with XFOIL's QDCALC output:

```python
# In compare script
def compare_dij(xfoil_dij, rust_dij):
    print("DIJ Diagonal comparison:")
    for i, (xf, rf) in enumerate(zip(xfoil_dij, rust_dij)):
        diff = abs(xf - rf) / abs(xf) if xf != 0 else 0
        print(f"  [{i}] XFOIL={xf:.6e}  Rust={rf:.6e}  diff={diff*100:.2f}%")
```

---

## Phase 5: Likely Bug Locations

Based on failure pattern (CD=321, no convergence), probable issues:

### 5.1: DIJ Matrix (`dij.rs`)

XFOIL's QDCALC builds influence matrix from panel geometry.
Current RustFoil implementation may be:
- Missing curvature terms
- Wrong sign convention
- Incorrect self-influence (diagonal)

### 5.2: Edge Velocity Update (`update.rs`)

`set_edge_velocities()` may:
- Apply wrong delta to Ue
- Not account for mass defect definition (δ* vs δ*×Ue)
- Miss inviscid baseline restoration

### 5.3: Newton System Build (`newton.rs`)

`BlNewtonSystem::build()` may:
- Not include Ue coupling column (column 4 of VS1/VS2)
- Wrong residual scaling
- Missing wake coupling terms

---

## Execution Checklist

1. [x] Add `--debug` flag to CLI viscous command
2. [x] Instrument UPDATE events in `update_stations()`
3. [x] Run RustFoil with debug output
4. [x] Run `debug_compare.py --stage mrchue`
5. [x] Identify first divergence point → **FOUND: MRCHUE setup**

---

## DIAGNOSIS RESULTS (2026-01-21)

### Root Cause Identified: Setup Before March

The **first divergence point is in the initial setup**, before `march_fixed_ue` even runs:

| Field | XFOIL | RustFoil | Issue |
|-------|-------|----------|-------|
| x (station 2) | 0.0011 | 0.0079 | Wrong arc length reference |
| Ue (station 2) | 0.061 | 0.885 | Not starting from stagnation |
| θ (station 2) | 2.3e-05 | 7.9e-02 | ~3,500x wrong |

### What XFOIL Does

1. **Finds stagnation point** on airfoil
2. **Splits into upper/lower surfaces** from stagnation
3. **Computes arc length from stagnation** (x=0 at stagnation)
4. **Edge velocity starts at 0** and increases (stagnation acceleration)
5. **Marches each surface** separately

### What RustFoil Does Wrong

1. Uses **total arc length** from first panel (not from stagnation)
2. Uses **raw gamma values** as Ue (not stagnation-referenced)
3. Doesn't split surfaces properly
4. θ initialization is wrong (starts ~0.1 instead of ~2e-5)

### Fix Required

The fix is in `run_viscous_analysis()` in `rustfoil-cli/src/main.rs`:

```rust
// Current (WRONG):
let arc_lengths = compute_arc_lengths(&node_x, &node_y);  // From LE
let ue_inviscid = inv_solution.gamma.clone();  // Raw gamma

// Should be:
// 1. Find stagnation index
let stag_idx = find_stagnation(&ue_inviscid);

// 2. Split into upper/lower from stagnation
let (upper_x, upper_ue) = extract_surface(stag_idx, n, &node_x, &node_y, &ue_inviscid, true);
let (lower_x, lower_ue) = extract_surface(stag_idx, n, &node_x, &node_y, &ue_inviscid, false);

// 3. Arc length from stagnation (x=0 at stagnation)
let upper_arc = compute_arc_from_stagnation(&upper_x, &upper_y);
let lower_arc = compute_arc_from_stagnation(&lower_x, &lower_y);

// 4. March each surface separately
let upper_result = march_fixed_ue(&upper_arc, &upper_ue, re, msq, &config);
let lower_result = march_fixed_ue(&lower_arc, &lower_ue, re, msq, &config);
```

### Next Steps

1. [ ] Implement `extract_surface()` to split airfoil at stagnation
2. [ ] Compute arc lengths from stagnation for each surface
3. [ ] Ensure Ue is properly scaled (edge velocity, not circulation)
4. [ ] Re-run comparison to verify MRCHUE matches

---

## Success Criteria

Debugging is complete when:
- [ ] MRCHUE outputs match within 1%
- [ ] BLVAR outputs match within 1% per iteration
- [ ] UPDATE deltas match within 10% (sign and magnitude)
- [ ] Convergence history follows same pattern (RMS decreasing)
- [ ] Final CL, CD within 1% of XFOIL
