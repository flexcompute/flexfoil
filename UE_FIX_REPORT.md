# Edge Velocity (Ue) Bug Fix Report

## Problem Summary
At α=10°, Re=20M (NACA 0012), transition was occurring 43% late in RustFoil compared to XFOIL:
- **XFOIL**: transition at x = 0.0117, CD = 0.00793
- **RustFoil (before fix)**: transition at x = 0.0167 (+43%), CD = 0.01342 (+69%)

Root cause: Edge velocity (Ue) values were systematically too low, causing Reynolds number (Rθ = θ × Ue × Re) to stay below critical longer, delaying N-factor accumulation and transition.

## Root Cause Identified

**Bug Location**: Test files in `crates/rustfoil-solver/tests/`

**The Bug**:
```rust
// WRONG - caps Ue at MAXIMUM of 0.01
let ue_stag = setup_result.setup.ue_inviscid[ist].abs().min(0.01);

// CORRECT - ensures Ue is at MINIMUM of 0.01  
let ue_stag = setup_result.setup.ue_inviscid[ist].abs().max(0.01);
```

When gamma at stagnation was 0.078, `.min(0.01)` reduced it to 0.01.
When gamma at stagnation was 0.5, `.min(0.01)` reduced it to 0.01.

This capped stagnation Ue values, which then affected the entire surface extraction.

## Files Fixed

1. `crates/rustfoil-solver/tests/xfoil_viscous_comparison.rs` (2 occurrences)
2. `crates/rustfoil-solver/tests/test_ue_viscous.rs` (1 occurrence)
3. `crates/rustfoil-solver/tests/test_station_mapping.rs` (1 occurrence)

## Results After Fix

**α=10°, Re=20M (NACA 0012)**:
- **Transition**: x = 0.0133 (XFOIL: 0.0117) → Now only 14% late (was 43%)
- **CD**: 0.01000 (XFOIL: 0.00793) → Now 26% high (was 69%)
- **CL**: 1.1980 (XFOIL: ~1.20) → Within 0.2%

## Verification

Inviscid gamma values confirmed correct throughout:
- LE region: gamma = 0.75-1.1 ✓
- Stagnation: gamma ≈ 0.078 ✓
- Surface extraction now passes correct Ue values to BL solver ✓

## Remaining Issues

While the fix significantly improved results (14% late vs 43% late), transition timing is still not perfect. Possible remaining causes:

1. **N-factor calculation**: eⁿ amplification model may need tuning
2. **Initial θ values**: First few BL stations may need better initialization
3. **Turbulent transition modeling**: TRDIF (transition differential equation) accuracy
4. **Arc length calculation**: Small differences in surface parameterization
5. **Ncrit calibration**: Using Ncrit=9.0, may need case-specific adjustment

## Impact

This was a **critical bug** that affected all viscous test cases. The fix:
- Reduces transition location error by 67% (from 43% to 14% late)
- Reduces drag coefficient error by 62% (from 69% to 26% high)
- Dramatically improves early-station Rθ values
- Makes viscous-inviscid coupling more accurate

## Next Steps

1. ✅ Fix `.min(0.01)` → `.max(0.01)` in test files
2. ✅ Verify inviscid solver outputs correct gamma
3. ✅ Verify surface extraction maps gamma→Ue correctly  
4. 🔄 Investigate remaining 14% transition delay
   - Check N-factor accumulation rates
   - Compare early-station θ evolution with XFOIL
   - Verify Thwaites/similarity initial conditions
5. 🔄 Validate fix across multiple test cases (different α, Re, airfoils)

## Lessons Learned

- **Always use `.max()` for lower bounds**, `.min()` for upper bounds
- Edge velocity normalization is critical for BL coupling
- Small errors near stagnation propagate through entire BL solution
- Transition prediction is extremely sensitive to early-station θ values
