# Re-March Fix Test Report: α=15°

## Executive Summary

**Status: ⚠️  PARTIAL SUCCESS**

The re-marching mechanism **IS WORKING** (stations appear multiple times per iteration), but the solution is **NOT CONVERGING PROPERLY**:
- ✅ Re-march is happening (ibl=93 appears 9 times)
- ❌ H is decreasing instead of growing (2.386 → 2.202)
- ❌ CL is 16% too high (1.797 vs XFOIL's 1.553)
- ❌ Newton solver shows instability in late iterations

---

## 1. Re-March Verification

### Station Coverage (Upper Surface)
- **Total unique stations**: 99
- **Trailing edge station**: ibl=93 (x≈0.9986)
- **TE appears**: 9 times

**✅ PASS**: The re-march is working! Each station appears once per global Newton iteration.

### Sample Station Counts
```
ibl=93  (x≈0.9986): 9 occurrences  [TE]
ibl=92  (x≈0.9820): 9 occurrences
ibl=91  (x≈0.9652): 9 occurrences
ibl=88  (x≈0.9149): 9 occurrences
ibl=83  (x≈0.8308): 9 occurrences
```

All stations show 9 occurrences, confirming the boundary layer is being re-marched at each Newton iteration.

---

## 2. H Evolution at Trailing Edge

### Table: H at ibl=93 (x≈0.9986) Across Iterations

| Occurrence | H (Hk) | Cf       | Ue     | theta    | dstar    |
|------------|--------|----------|--------|----------|----------|
| 1          | 2.386  | 0.000289 | 1.0207 | 0.009836 | 0.023473 |
| 2          | 2.340  | 0.000318 | 1.0207 | 0.009684 | 0.022664 |
| 3          | 2.351  | 0.000312 | 1.0174 | 0.009688 | 0.022771 |
| 4          | 2.370  | 0.000300 | 1.0138 | 0.009708 | 0.023007 |
| 5          | 2.476  | 0.000240 | 1.0087 | 0.009923 | 0.024570 |
| 6          | 2.432  | 0.000264 | 1.0090 | 0.009849 | 0.023950 |
| 7          | 2.133  | 0.000477 | 1.0398 | 0.009170 | 0.019559 |
| 8          | 2.166  | 0.000448 | 1.0335 | 0.009254 | 0.020048 |
| 9          | 2.202  | 0.000419 | 1.0274 | 0.009334 | 0.020549 |

### H Statistics
- **First H**: 2.386
- **Last H**: 2.202
- **Min H**: 2.133
- **Max H**: 2.476
- **Range**: 0.343
- **Growth**: -0.185 (-7.8%)

**❌ FAIL**: H is **decreasing** instead of growing! Expected H to grow from ~2.4 → ~3.0, but it decreased by 7.8%.

### Observations
- Iterations 1-6: H oscillates around 2.3-2.5
- Iteration 7-9: H drops sharply to ~2.1-2.2
- The sharp drop at iteration 7 coincides with large Newton residuals in the test output

---

## 3. Final Forces Comparison

### RustFoil vs XFOIL at α=15°

| Metric | RustFoil | XFOIL | Error |
|--------|----------|-------|-------|
| **CL** | 1.7969   | 1.5528 | +15.7% |
| **CD** | 0.02614  | 0.01926 | +35.7% |
| **H (upper TE)** | 2.202 | 2.294 | -4.0% |
| **H (lower TE)** | 1.352 | 3.033 | -55.4% |

**❌ FAIL**: CL is 16% too high, suggesting the solution has not converged to the correct stall state.

---

## 4. Convergence Analysis

From test output line 994:
```
[DEBUG Newton] iter 7 relaxation used: 1.000000e-2 (requested: 1.000000e0), rms_change: 1.198909e2
```

### Issues Observed:
1. **Iteration 7 instability**: Extremely large RMS change (119.9), relaxation reduced to 0.01
2. **Ue changes are erratic**: Upper ibl=20 shows due=-3.541 (line 984)
3. **Delta-star warnings**: Multiple "Reduced rlx" messages due to large dn3 changes
4. **Early termination**: Only 8 Newton iterations before stopping

### Likely Root Causes:
- **Inviscid Ue updates are too aggressive**: The coupling between viscous and inviscid is unstable
- **Missing damping**: No damping on Ue updates when H is near stall
- **Incorrect stagnation handling**: The test shows "x_tr upper = 0.0088" which seems wrong

---

## 5. Comparison to Expected Behavior

### Expected (from XFOIL)
At α=15° near stall:
- CL ≈ 1.55
- Upper surface TE: H ≈ 2.3 (attached flow)
- Lower surface TE: H ≈ 3.0 (approaching separation)
- Smooth convergence in ~9 iterations

### Actual (RustFoil)
- CL = 1.80 (too high, over-predicting lift)
- Upper surface TE: H oscillating 2.1-2.5
- Lower surface TE: H = 1.35 (way too low, suggests error)
- Newton solver is unstable, large relaxation reductions

---

## 6. Success Criteria Assessment

| Criterion | Status | Details |
|-----------|--------|---------|
| ibl=93 appears 8-20 times | ✅ PASS | 9 occurrences |
| H grows from ~2.4 → ~3.0 | ❌ FAIL | H decreased -7.8% |
| Max H ≥ 2.9 | ❌ FAIL | Max H = 2.476 |
| CL within 10% of XFOIL | ❌ FAIL | +15.7% error |

**Overall: ❌ FAIL** - While re-marching works, the solution is not converging correctly.

---

## 7. Next Steps (Recommended)

### Immediate Fixes Needed:
1. **Investigate Ue update instability** (line 984-993 in test output)
   - The due values are huge: -3.54 for upper, +3.29 for lower
   - This suggests the viscous-inviscid coupling is broken

2. **Check stagnation point handling**
   - Test output shows "x_tr upper = 0.0088" which is suspiciously low
   - May be confusing transition with stagnation

3. **Add damping for high-α cases**
   - Implement adaptive relaxation based on H values
   - Reduce Ue updates when H > 2.5

4. **Verify lower surface calculation**
   - Lower surface H=1.35 at TE is wrong (should be ~3.0)
   - This is a major error independent of the re-march fix

### Debugging Strategy:
1. Test at lower α (8°-10°) where flow is attached
2. Compare station-by-station H values to XFOIL
3. Add more debug output for Ue updates
4. Check if the issue is in `compute_new_ue` or `apply_global_updates`

---

## 8. Conclusion

The **critical stall bug fix is partially working**:
- ✅ Re-marching mechanism is functional
- ❌ Solution convergence is broken
- ❌ CL prediction is 16% too high
- ❌ H evolution is incorrect (decreasing instead of growing)

**Root cause appears to be in the viscous-inviscid coupling**, not the re-march mechanism itself. The large Ue updates (line 984-993) and early relaxation reductions suggest the Newton solver is fighting against itself.

**Recommendation**: Fix the Ue update stability issues before retesting the stall prediction capability.
