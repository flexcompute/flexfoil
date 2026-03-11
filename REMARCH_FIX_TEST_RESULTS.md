# Re-March Fix Test Results

## Executive Summary

**Status: ⚠️  PARTIAL SUCCESS - Critical Issues Found**

### What Works ✅
- **Re-marching mechanism is functional**: Stations appear multiple times per iteration (8-21 times)
- **March surface is being called**: The boundary layer is re-computed each Newton iteration

### What's Broken ❌
- **H values at TE are completely wrong**: Off by 100-150% compared to XFOIL
- **CL is over-predicted by 16% at α=15°**: 1.797 vs 1.553
- **Newton solver shows instability**: Large Ue updates, excessive relaxation reductions
- **Lower surface H is wrong**: 1.35 vs XFOIL's 3.03 at α=15°

---

## Test Results: α=15° (Near Stall)

### Re-March Verification
- **Trailing edge station**: ibl=93 (x≈0.9986)
- **Occurrences**: 9 times ✅
- **Conclusion**: Re-marching IS working

### H Evolution at TE (ibl=93)

| Iteration | H (Hk) | Cf       | Ue     | Status |
|-----------|--------|----------|--------|--------|
| 1         | 2.386  | 0.000289 | 1.0207 | Start  |
| 2         | 2.340  | 0.000318 | 1.0207 |        |
| 3         | 2.351  | 0.000312 | 1.0174 |        |
| 4         | 2.370  | 0.000300 | 1.0138 |        |
| 5         | 2.476  | 0.000240 | 1.0087 | Peak   |
| 6         | 2.432  | 0.000264 | 1.0090 |        |
| 7         | 2.133  | 0.000477 | 1.0398 | Drop!  |
| 8         | 2.166  | 0.000448 | 1.0335 |        |
| 9         | 2.202  | 0.000419 | 1.0274 | End    |

**Problem**: H decreased by 7.8% instead of growing! Expected H to grow to ~3.0.

### Forces Comparison

| Metric            | RustFoil | XFOIL  | Error    |
|-------------------|----------|--------|----------|
| CL                | 1.7969   | 1.5528 | **+15.7%** ❌ |
| CD                | 0.02614  | 0.01926 | +35.7%   |
| H (upper TE)      | 2.202    | 2.294  | -4.0%    |
| H (lower TE)      | 1.352    | 3.033  | **-55.4%** ❌ |
| x_tr (upper)      | 0.0088   | 0.0104 | -15.4%   |

**Critical Issues**:
1. CL is 16% too high - solution not converged
2. Lower surface H is completely wrong (1.35 vs 3.03)
3. Upper transition location seems wrong (0.0088 is suspiciously low)

### Newton Convergence Issues

From test output (iteration 7):
```
[DEBUG Newton] iter 7 relaxation used: 1.000000e-2 (requested: 1.000000e0), rms_change: 1.198909e2
[DEBUG Ue UPDATE] upper ibl=20: new_u=0.010000, current_u=3.551501, due=-3.541501e0
```

**Problems**:
- RMS change exploded to 119.9 at iteration 7
- Ue updates are huge: due=-3.54 (change of 350%!)
- Relaxation reduced to 0.01 (99% damping)
- Multiple "Reduced rlx" warnings due to delta_star changes

---

## Test Results: α=10° (Attached Flow)

### Re-March Verification
- **Trailing edge station**: ibl=90 (x≈0.9962)
- **Occurrences**: 21 times ✅
- **Conclusion**: Re-marching works better here (more iterations)

### H Evolution at TE (ibl=90)

| Iteration | H (Hk) | Trend |
|-----------|--------|-------|
| 1         | 1.736  | Start |
| 5         | 1.772  | ↑     |
| 10        | 1.802  | ↑     |
| 15        | 1.816  | ↑     |
| 21        | 1.827  | End   |

**Good**: H is growing monotonically (+5.2%)

### Forces Comparison

| Metric            | RustFoil | XFOIL  | Error      |
|-------------------|----------|--------|------------|
| CL                | 1.0716   | 1.1176 | **-4.1%** ⚠️ |
| CD                | 0.01342  | 0.01125 | +19.3%     |
| H (upper TE)      | 1.827    | 4.912  | **-62.8%** ❌ |
| H (lower TE)      | 1.168    | 1.971  | **-40.8%** ❌ |
| x_tr (upper)      | 0.0167   | 0.0186 | -10.2%     |

**Critical Finding**: RustFoil's H at TE is **completely wrong**!
- Upper TE: 1.827 vs XFOIL's 4.912 (off by 169%!)
- This is NOT a convergence issue - the solution converged well at α=10°
- This suggests a fundamental error in the H calculation or boundary layer integration

---

## Analysis: What's Going Wrong?

### 1. H at TE is Systematically Too Low

**α=10°**: RustFoil H=1.83, XFOIL H=4.91 (63% error)
**α=15°**: RustFoil H=2.20, XFOIL H=2.29 (4% error, but only after wrong convergence)

**Possible Causes**:
- Boundary layer integration errors near the trailing edge
- Incorrect application of the Kutta condition
- Displacement thickness calculation error
- Shape parameter closure equations error

### 2. Viscous-Inviscid Coupling is Unstable

At α=15°, iteration 7 shows:
```
due=-3.541 at ibl=20 (99% change!)
rlx reduced to 0.01 due to delta_star changes
```

**Possible Causes**:
- `compute_new_ue` is too aggressive with updates
- No damping for high-α or near-separation cases
- Jacobian is inaccurate when H is large

### 3. Lower Surface H is Always Wrong

**α=10°**: RustFoil H=1.17, XFOIL H=1.97 (41% error)
**α=15°**: RustFoil H=1.35, XFOIL H=3.03 (55% error)

**This is independent of convergence issues!**

### 4. Re-March is Working But Exposes Other Bugs

The re-march mechanism itself is functional:
- Stations appear multiple times ✅
- `march_surface` is being called ✅
- Results are copied back to stations ✅

BUT the underlying BL calculations are wrong, so re-marching just iterates towards the wrong answer!

---

## Root Cause Hypothesis

Based on the evidence, the most likely root cause is **incorrect boundary layer integration near the trailing edge**:

1. **H grows correctly in the middle of the airfoil** (α=10°: H grows from 1.7 → 1.8)
2. **But fails near the TE** (should be 4.9, got 1.8)
3. **Lower surface is consistently wrong** (always 30-55% too low)
4. **CL errors correlate with H errors** (α=10° has -4% CL error, α=15° has +16% CL error)

**Likely Issues**:
- Trailing edge panel handling in `march_surface`
- Wake initialization using wrong H values
- Displacement thickness incorrectly calculated at TE
- Kutta condition application error

---

## Success Criteria Assessment

| Criterion                           | α=10° | α=15° | Overall |
|-------------------------------------|-------|-------|---------|
| ibl=TE appears 8-20 times           | ✅ 21 | ✅ 9  | ✅ PASS |
| H grows (not decreases)             | ✅ +5.2% | ❌ -7.8% | ⚠️  MIXED |
| H reaches expected value            | ❌ 1.83 vs 4.91 | ❌ 2.20 vs 2.29 | ❌ FAIL |
| CL within 10% of XFOIL              | ✅ -4.1% | ❌ +15.7% | ⚠️  MIXED |
| Final convergence                   | ✅ Good | ❌ Unstable | ⚠️  MIXED |

**Overall Verdict**: ❌ **FAIL** - Re-march works but exposes fundamental BL errors

---

## Recommended Next Steps

### Immediate Priority (Blocking):
1. **Investigate TE boundary layer calculation**
   - Compare RustFoil vs XFOIL station-by-station from x=0.9 to 1.0
   - Check displacement thickness, momentum thickness, and H at each station
   - Verify wake initialization

2. **Debug lower surface entirely**
   - Lower surface H is wrong at ALL angles tested
   - Check if it's a sign error, indexing error, or closure equation error

3. **Test at lower angles first** (α=0°, 4°)
   - Eliminate separation/transition complications
   - Focus on attached flow accuracy

### Medium Priority:
4. **Fix Ue update stability**
   - Add adaptive damping based on H values
   - Limit maximum Ue change per iteration
   - Better relaxation strategy for near-stall cases

5. **Add station-by-station comparison test**
   - Create test that compares RustFoil vs XFOIL at each station
   - Identify exactly where H calculation goes wrong

### Lower Priority:
6. **Re-test stall prediction** (after fixes above)
   - Only test stall after TE boundary layer is correct
   - Current results are meaningless due to underlying errors

---

## Conclusion

The **critical stall bug fix (re-marching at each Newton iteration) is technically working**, but it has exposed much more serious bugs in the boundary layer calculations:

1. ✅ Re-marching mechanism works correctly
2. ❌ H at trailing edge is systematically wrong (30-170% error)
3. ❌ Lower surface calculations are completely broken
4. ❌ Viscous-inviscid coupling is unstable at high α
5. ❌ Cannot test stall prediction until basic BL accuracy is fixed

**The test reveals that the original stall bug was masking these fundamental errors.** Now that we're re-marching, the errors compound and prevent convergence.

**Action Required**: Fix TE boundary layer calculations and lower surface errors before proceeding with stall testing.

---

## Appendix: Debug Commands Used

```bash
# Run RustFoil at α=15°
RUSTFOIL_TEST_ALPHA=15 RUSTFOIL_FULL_TRACE=1 cargo test --package rustfoil-solver \
  --test xfoil_viscous_comparison test_newton_iteration_trace -- --nocapture

# Run RustFoil at α=10°
RUSTFOIL_TEST_ALPHA=10 RUSTFOIL_FULL_TRACE=1 cargo test --package rustfoil-solver \
  --test xfoil_viscous_comparison test_newton_iteration_trace -- --nocapture

# Run XFOIL for comparison
cd Xfoil-instrumented/bin
./xfoil_instrumented << EOF
NACA 0012
OPER
VISC 3e6
ALFA 15
DUMP /tmp/xfoil_a15_verify.txt

QUIT
EOF

# Analyze traces
python3 scripts/analyze_remarch_fix.py
```

## Appendix: Key Files

- Test output: `test_remarch_a15.txt` (in workspace root)
- RustFoil trace: `crates/rustfoil-solver/traces/rustfoil_new/rustfoil_alpha_{10,15}.json`
- XFOIL dump: `/tmp/xfoil_a{10,15}_verify.txt`
- Analysis script: `scripts/analyze_remarch_fix.py`
- This report: `REMARCH_FIX_TEST_RESULTS.md`
