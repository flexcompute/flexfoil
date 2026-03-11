---
created: 2026-01-28
source: flexfoil-boundary-layer
tags: [type/plan, status/active, topic/debugging, topic/xfoil-matching]
---

# Newton Iteration Brute Force Comparison Plan

## Context
After extensive debugging of the [[rustfoil-bl]] viscous solver, we've narrowed down the issue to the Newton iteration not converging to the same values as XFOIL. The DUI (mass defect velocity correction) is 32% smaller than XFOIL's, but the root cause remains unclear despite:
- Same formula
- Same VTI signs
- Similar DIJ values (within 4%)
- Same mass values

This plan outlines a systematic brute-force comparison approach.

## Objective
Find the exact point of divergence between XFOIL and RustFoil by comparing every intermediate value at every Newton iteration across multiple angles of attack.

## Test Cases

| Alpha | Expected Behavior |
|-------|-------------------|
| 0° | Symmetric, simplest case |
| 2° | Current test case, small asymmetry |
| 4° | Moderate lift |
| 8° | Higher loading, approaching stall |

## Data to Dump at Each Newton Iteration

### 1. Global State (per iteration)

```
RMSBL / rms_change - convergence metric
ITRAN(1), ITRAN(2) - transition station indices
x_tr_upper, x_tr_lower - transition locations
```

### 2. Per-Station BL Variables (every 10th station)

```
IBL, IS, x, theta, dstar, Ue, mass, Hk, Cf, Ctau, Ampl
```

### 3. Newton System Components (every 10th station)

```
IBL, IS:
  - VS1[3x5] (upstream Jacobian)
  - VS2[3x5] (downstream Jacobian)  
  - VSREZ[3] (raw residuals before forced changes)
  - VDEL[3] (residuals after forced changes)
  - DUE (Ue mismatch)
```

### 4. DIJ-Related Values (once per iteration)

```
For stations 1-5 on each surface:
  - DUI_upper, DUI_lower, DUI_total
  - ue_inviscid, ue_from_mass, ue_current
```

### 5. Solution Deltas (per iteration, every 10th station)

```
IBL, IS:
  - delta_ctau, delta_theta, delta_mass
  - RLX (relaxation factor applied)
```

## Implementation Plan

### Phase 1: Instrumentation

1. **XFOIL**: Create `DBGFULLITER` subroutine that dumps all above data to JSON
2. **RustFoil**: Add matching debug output with identical format
3. **Output format**: One JSON file per (alpha, solver) combination

### Phase 2: Comparison Script

```python
# compare_newton_iterations.py

def compare_iterations(xfoil_json, rustfoil_json):
    """
    Load XFOIL and RustFoil JSON for same alpha
    For each iteration:
      - Compare global state
      - Find first station where values diverge > 1%
      - Report: "Iteration X, Station Y, Variable Z diverged: XFOIL=A, RF=B"
    Generate summary: "First divergence at alpha=X, iter=Y, station=Z"
    """
    pass
```

### Phase 3: Execution

```bash
# Run XFOIL at each alpha, dump traces
for alpha in 0 2 4 8; do
  ./xfoil_instrumented < run_alpha_$alpha.inp
done

# Run RustFoil at each alpha, dump traces  
for alpha in 0 2 4 8; do
  RUSTFOIL_FULL_TRACE=1 cargo test --release ... -- alpha_$alpha
done

# Compare
python3 compare_newton_iterations.py
```

## Expected Output Format

```
=== Alpha = 0° ===
Iteration 1: Match through station 80
Iteration 2: DIVERGENCE at station 45 (lower)
  theta: XFOIL=3.21e-5, RF=3.45e-5 (+7.5%)
  First cause: DUE at station 2 differs by 12%

=== Alpha = 2° ===  
Iteration 1: DIVERGENCE at station 3
  VDEL[1]: XFOIL=-12.4, RF=-8.9 (-28%)
  Root cause: DUI sum differs
```

## Deliverables

1. `Xfoil-instrumented/src/xfoil_debug.f` - new DBGFULLITER subroutine
2. `crates/rustfoil-coupling/src/debug_trace.rs` - matching RustFoil output
3. `scripts/compare_newton_iterations.py` - comparison tool
4. `docs/debugging/newton_comparison_results.md` - findings

## Time Estimate

| Phase | Estimate |
|-------|----------|
| Phase 1 (Instrumentation) | 2-3 hours |
| Phase 2 (Comparison script) | 1 hour |
| Phase 3 (Run & analyze) | 1-2 hours |
| **Total** | **4-6 hours** |

## Current Status

**2026-01-28**: Plan created. Ready to begin Phase 1.

## Related

- [[rustfoil-xfoil-methodology-comparison]]
- [[XFOIL_VISCOUS_ARCHITECTURE]]
- [[newton-vi-coupling-plan]]
