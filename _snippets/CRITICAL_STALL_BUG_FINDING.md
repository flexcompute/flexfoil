---
created: 2026-01-28
type: critical-bug-report
tags: [bug/critical, boundary-layer, stall-prediction]
priority: HIGHEST
---

# CRITICAL: RustFoil Does Not Re-March BL During Global Iterations

## Summary

RustFoil fails to predict stall because it marches the boundary layer ONCE, then uses Jacobian-based Newton updates for subsequent global iterations. XFOIL re-marches the BL at EVERY global iteration, allowing nonlinear effects like H > htmax at the trailing edge.

## Evidence

### XFOIL Behavior at α=15° (NACA 0012)
```
MRCHUE: Inverse mode at  93     Hk =   2.500  (TE during marching)
MRCHUE: Inverse mode at  94     Hk =   2.500
...
MRCHUE: Inverse mode at  99     Hk =   2.500

Final converged state:
  TE (x=1.0): H = 3.033, Cf = 0.000034 (separating!)
```

XFOIL's BL march uses inverse mode with htarg=2.5 at the TE. But after multiple global viscous-inviscid iterations, the final H reaches 3.033 (exceeds htmax!).

### RustFoil Behavior at α=15°
```
ibl=93 (x=0.9986, TE): Appears in MRCHUE events ONCE
  Final state: Hk = 2.378, Cf = 0.000294 (attached)
  
Global Newton iterations: 20
  But MRCHUE events show each station only ONCE!
```

RustFoil marches once, then does 20 global Newton iterations using only Jacobian updates. H stays at 2.378.

### Code Confirmation

**RustFoil `solve_viscous_two_surfaces` (viscal.rs:584):**

```rust
// Line 627: March ONCE before Newton loop
let upper_result = march_surface(&upper_arc, upper_ue, re, msq, &march_config, 1);
let lower_result = march_surface(&lower_arc, lower_ue, re, msq, &march_config, 2);

// Lines 809-1000: Newton iteration loop
for iter in 0..max_newton_iter {
    // Update secondary variables
    blvar(station, flow_type, msq, re);  // Line 865
    
    // Build Jacobian system
    global_system.build_global_system(...);  // Line 939
    
    // Solve and update
    // ... NO march_surface call!
}
```

## Why This Matters

1. **Nonlinear Effects:** The relationship between Ue and H is highly nonlinear, especially near separation. Small changes in Ue can cause large H changes.

2. **Inverse Mode Limiting:** During marching, inverse mode limits Hk to htmax. But as Ue evolves through viscous-inviscid coupling, the CONVERGED H can exceed htmax.

3. **TE Separation:** At high alpha, TE separation occurs when H ~ 3.0 for turbulent flow. If we can't exceed htmax=2.5, we can't capture this!

4. **Stall Prediction:** Stall occurs when TE separation causes a large drop in CL. Without proper H growth, no stall is predicted.

## The Fix

**Option 1: Re-march at each iteration (XFOIL-like)**
```rust
for iter in 0..max_newton_iter {
    // Re-march both surfaces with updated Ue
    let upper_result = march_surface(&upper_arc, &updated_ue_upper, re, msq, &march_config, 1);
    let lower_result = march_surface(&lower_arc, &updated_ue_lower, re, msq, &march_config, 2);
    
    // Copy results back to stations
    // ...
    
    // Then build Jacobian for next correction
    global_system.build_global_system(...);
}
```

**Pros:**
- Matches XFOIL exactly
- Captures all nonlinear effects
- Guaranteed correct

**Cons:**
- Slower (march is expensive)
- May need relaxation tuning

**Option 2: Hybrid approach**
- March every N iterations (e.g., N=3)
- Use Jacobian updates in between
- Balance accuracy vs performance

## Impact

**Without fix:**
- ✅ Pre-stall CL within 5-10% of XFOIL
- ❌ No stall prediction at any alpha
- ❌ CD overpredicted by 30-50% at high alpha
- ❌ H capped at ~2.4, never reaches separation threshold

**With fix:**
- ✅ Should match XFOIL CL polar including stall
- ✅ Proper TE separation detection
- ✅ Correct CD at high alpha

## Recommended Action

1. Implement Option 1 (full re-march) first
2. Validate against XFOIL at α=0°, 10°, 15°, 20°
3. If performance is acceptable, ship it
4. If too slow, try Option 2 (march every 3 iterations)

## Related Files

- `crates/rustfoil-solver/src/viscous/viscal.rs` - Main viscal loop
- `crates/rustfoil-coupling/src/march.rs` - BL marching logic
- `_snippets/2026-01-28.md` - Daily investigation log
- `_snippets/BL_Investigation_Plan.md` - Investigation plan

---

**Status:** Root cause identified, fix ready to implement
**Assigned:** AI Assistant
**Priority:** CRITICAL - blocks stall prediction entirely
