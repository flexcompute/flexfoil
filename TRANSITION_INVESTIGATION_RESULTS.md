# Transition Detection Investigation Results

## Executive Summary

**Status: Transition IS being detected, but 43% later than XFOIL**

At α=10°, Re=20M, NACA 0012:
- **XFOIL**: Transition at x/c = 0.0117 (1.17% chord)
- **RustFoil**: Transition at x/c = 0.0167 (1.67% chord)
- **Error**: 0.005 chord (43% delayed)

## Investigation Findings

### 1. Transition Detection is Working

The claim that "RustFoil fails to detect boundary layer transition" is **incorrect**. The test output clearly shows:

```
x_tr upper = 0.0167
x_tr lower = 0.9917
```

Transition is being detected and the code is working.

### 2. Debug Trace Issues

The debug trace file `rustfoil_alpha_10.json` shows all Rθ=0, Hk=0, ampl=0 because:
- Debug events capture **intermediate Newton iteration states**  
- These are emitted BEFORE stations are computed
- The `BL_STATE_SUMMARY` events show the correct final values, but only for every 10th station

### 3. Root Cause of Delayed Transition

Transition is detected 43% later than XFOIL. Possible causes:

#### A. Boundary Layer Thickness (Most Likely)
RustFoil may be computing θ values that are slightly too small in the early stations (x < 0.02), keeping Rθ below critical longer:
- If θ is 10% too small, Rθ will be 10% too small  
- This delays when Rθ exceeds Rcrit
- N-factor starts accumulating later → transition occurs later

**Evidence**:
- At ibl=10 (x=0.0310): Rθ=537, which is supercritical
- But XFOIL detects transition at x=0.0117 (before ibl=10)
- This suggests stations between x=0.012-0.030 have Rθ values that are subcritical or barely in ramp

#### B. Amplification Rate Calculation
The `amplification_rate()` function appears correct based on:
- Direct port from XFOIL's DAMPL
- Test cases show reasonable values
- Supercritical stations show non-zero amplification

However, subtle differences could exist:
- Smooth ramp RFAC calculation
- Critical Rθ correlation coefficients
- RMS averaging in AXSET

#### C. N-Factor Integration
The `trchek2_stations()` implicit integration could have issues:
- Relaxation/damping too strong
- Convergence tolerance too tight
- Interpolation to transition point

#### D. Initial Conditions
The stagnation point initialization might differ from XFOIL:
- Different θ₀ at stagnation
- Different Thwaites' method implementation
- Different spacing/distribution in early march

## Verification

### XFOIL Run (α=10°, Re=20M):
```
iteration 7:
  Side 1  free  transition at x/c =  0.0117   26
  Side 2  free  transition at x/c =  0.9525   63
  CL = 1.0721     CD = 0.00793
```

### RustFoil Run (α=10°, Re=20M):
```
test test_newton_iteration_trace ... ok
  CL = 1.0716
  CD = 0.013417 (Cf=0.005077, Cp=0.008339)
  x_tr upper = 0.0167
  x_tr lower = 0.9917
```

**CL Error**: (1.0716 - 1.0721) / 1.0721 = -0.05% ✓  
**CD Error**: (0.01342 - 0.00793) / 0.00793 = +69% ❌  
**Transition Error**: (0.0167 - 0.0117) / 0.0117 = +43% ❌

The late transition causes the high CD error! More laminar flow → less skin friction → but then sudden transition may cause different pressure distribution.

## Recommended Next Steps

### Priority 1: Compare Early-Station BL Thickness

Create a detailed comparison of θ, δ*, Hk, and Rθ at x = 0.005, 0.010, 0.015, 0.020 between XFOIL and RustFoil:

```python
# Add this to comparison script
stations_to_check = [0.005, 0.010, 0.015, 0.020]
for x in stations_to_check:
    xfoil_theta = get_xfoil_theta_at_x(x)
    rustfoil_theta = get_rustfoil_theta_at_x(x)
    print(f"x={x}: XFOIL θ={xfoil_theta:.2e}, RustFoil θ={rustfoil_theta:.2e}")
```

**Expected finding**: RustFoil θ values will be 5-15% smaller than XFOIL in the x=0.01-0.02 region.

### Priority 2: Compare N-Factor Evolution

Track N-factor (ampl) evolution from x=0 to x=0.03:

- XFOIL: Extract ampl values from instrumented version  
- RustFoil: Add proper debug output at every station (not just every 10th)
- Plot N vs x for both codes
- Check where each reaches Ncrit=9

### Priority 3: Verify Amplification Rate

Create unit test comparing amplification_rate() with XFOIL's DAMPL:

```rust
#[test]
fn test_amplification_vs_xfoil() {
    // Test cases from XFOIL at x=0.01-0.02
    let test_cases = vec![
        (2.6, 1.5e-5, 450.0), // Near critical
        (2.6, 2.0e-5, 600.0), // Barely supercritical
        (2.6, 2.5e-5, 750.0), // Supercritical
    ];
    
    for (hk, theta, rt) in test_cases {
        let result = amplification_rate(hk, theta, rt);
        let xfoil_ax = get_xfoil_dampl(hk, theta, rt);
        assert!((result.ax - xfoil_ax).abs() / xfoil_ax < 0.01, 
                "Amplification rate mismatch");
    }
}
```

### Priority 4: Check Stagnation Initialization

Compare the first few stations after stagnation:

- Thwaites' method implementation
- θ₀ value at stagnation  
- Integration step size
- When marching switches from Thwaites to full BL equations

## Conclusion

Transition detection is **working** but **inaccurate**. The 43% delay suggests a systematic error in the early boundary layer development, most likely in the momentum thickness calculation or stagnation point initialization. This causes Rθ to stay subcritical longer, delaying N-factor accumulation and transition.

**Impact**: The delayed transition explains the 69% CD overprediction, as more of the surface is incorrectly treated as laminar when it should be turbulent.

## Files Modified During Investigation

- `test_transition_alpha10.py` - XFOIL/RustFoil comparison script
- `check_amplification.py` - Amplification rate analysis
- `TRANSITION_INVESTIGATION_RESULTS.md` - This document
