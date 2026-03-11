# Laminar Separation Bubble Investigation Summary

## Problem Statement
At α=10° on NACA 0012:
- **XFOIL**: Max H = 7.73 at x=0.0617 (Laminar Separation Bubble)
- **RustFoil**: Max H = 3.35 at x=0.0617 (No LSB captured)

This difference explains why RustFoil doesn't predict stall - the separation physics isn't being captured.

## Investigation Findings

### 1. htarg Upper Bound Fix (Implemented)
- Removed artificial upper bound on `htarg` that was capping at 4.8
- Now matches XFOIL's `HTARG = MAX(HTARG, HMAX)` (lower bound only)
- **Result**: Minor improvement (H: 3.08 → 3.35)

### 2. Closure Relations (Verified)
- Hs, Cf, DI formulas match XFOIL exactly
- No formula differences found at high H values

### 3. Edge Velocity Comparison
- Ue values **identical** up to x=0.0588
- Diverge at x=0.0617: XFOIL Ue=2.538, RustFoil Ue=2.487
- XFOIL's inverse mode adjusts Ue, RustFoil stays in direct mode

### 4. HLMAX Threshold Experiment
**Key finding**: Lowering HLMAX from 3.8 to 3.2:
- **Initial march** achieves H=5.52 at x=0.065 (close to XFOIL's 7.73)
- **Global Newton** converges back to H~3.3

This reveals the real issue: **VI coupling during global iteration suppresses the LSB**

## Root Cause Analysis

### The Two-Phase Problem

**Phase 1: Initial March**
- With lower HLMAX, inverse mode activates earlier
- htarg formula produces H=5.5 (good LSB physics)

**Phase 2: Global Newton Iteration**
- VI coupling adjusts Ue based on δ* changes
- Modified Ue has weaker adverse pressure gradient
- Re-march produces lower H (~3.3)
- LSB "disappears" during convergence

### Why XFOIL Preserves the LSB
XFOIL's global Newton includes:
1. "Forced changes" (DDS, DUE terms) that maintain BL-inviscid consistency
2. More sophisticated Ue update that preserves separation characteristics
3. Different relaxation strategy during global iteration

## Technical Details

### BLDIF Residual Comparison at x=0.0617
| Parameter | XFOIL | RustFoil |
|-----------|-------|----------|
| VSREZ[1] (θ correction) | +0.20 | +0.11 |
| VSREZ[2] (δ* correction) | -0.063 | -0.003 |

XFOIL produces 10× larger δ* corrections, driving H growth.

### H Values at Critical Station (x=0.065)
| Stage | XFOIL | RustFoil (hlmax=3.8) | RustFoil (hlmax=3.2) |
|-------|-------|---------------------|---------------------|
| Initial march | 7.73 | 3.35 | 5.52 |
| After global Newton | 7.73 | 3.35 | 3.35 |

## Conclusions

1. **The marching equations can produce high H** when inverse mode activates early (proven by hlmax=3.2 experiment)

2. **The global VI coupling suppresses the LSB** by:
   - Adjusting Ue to reduce adverse pressure gradient
   - Not preserving the separation-induced changes

3. **CL error source**: Without the LSB:
   - Flow appears attached at high α
   - No CL reduction from separation
   - RustFoil overpredicts CL at high angles

## Recommendations

### Immediate Fix (Try First)
Lower HLMAX to 3.2-3.5 AND modify global Newton to preserve separation:
```rust
// In viscal.rs global Newton loop:
// Don't fully relax Ue changes that would suppress separation
// Or: skip re-march when LSB is detected
```

### Medium-Term
1. Implement XFOIL's "forced change" mechanism in global Newton
2. Add LSB detection with special handling during VI coupling
3. Use different relaxation for stations in inverse mode

### Long-Term
1. Full port of XFOIL's SETBL/BLSOLVE global Newton system
2. Proper handling of transition in separated shear layer
3. Turbulent reattachment modeling

## Files Modified
- `crates/rustfoil-coupling/src/march.rs`: Removed htarg upper bound

## Test Commands
```bash
# Run comparison at α=10°
python3 scripts/compare_stall_physics.py 10

# Run with debug trace
RUSTFOIL_TEST_ALPHA=10 RUSTFOIL_FULL_TRACE=1 cargo test --release -p rustfoil-solver \
  --test xfoil_viscous_comparison -- test_newton_iteration_trace --nocapture
```
