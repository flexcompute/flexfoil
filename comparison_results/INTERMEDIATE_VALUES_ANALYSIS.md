# Intermediate Values Discrepancy Analysis

**Purpose**: Identify which intermediate computational values differ between RustFoil and XFOIL, causing the remaining CL/CD discrepancies.

## Final Results Summary (After All Fixes)

### CL Accuracy
- **Mean error**: 4.0% (was -24% before fixes)
- **RMS error**: 12.7%
- **Best region**: 3°-8° with 9-10% error
- **Worst case**: α=-7° with 29% error

### CD Accuracy  
- **Mean error**: 124% (was 3000% before fixes)
- **RMS error**: 419%
- **Issues**: Many angles report CD=0 (likely CDP calculation issue)
- **Anomaly remaining**: α=-7° with CD=0.0446 (should be 0.0021) → 2000% error

##Key Remaining Discrepancies

### 1. **CL Systematically High** (3°-15°)

All angles from 3° to 15° show CL 9-16% too high. This suggests a systematic issue in one of:

**A. Boundary Layer Thickness**
- **Hypothesis**: BL is still too thin (low theta → low mass defect → high lift)
- **Test**: Compare theta at x=0.5, 0.9, 1.0 for α=4°
- **Files**: `crates/rustfoil-bl/src/equations.rs` (BLDIF), closures

**B. Transition Location**
- **Known**: RustFoil transitions ~5% earlier than XFOIL
- **Impact**: Earlier transition → thicker turbulent BL → more drag, less lift
- **But**: This should REDUCE lift, not increase it
- **Conflict**: Suggests other factor is dominating

**C. Edge Velocity (UESET)**
- **Hypothesis**: DIJ coupling gives slightly higher Ue than XFOIL
- **Test**: Compare Ue at station 1, 10, 50 for α=4°
- **Files**: `crates/rustfoil-coupling/src/global_newton.rs` (compute_ue_from_mass_both)

**D. Force Calculation (CL Integration)**
- **Current**: Uses circulation-based CL = ∫(γ_upper - γ_lower) ds
- **XFOIL**: Uses multiple methods, cross-validates
- **Test**: Compare Cp distribution and verify integration

### 2. **CD = 0 at Multiple Angles**

RustFoil reports CD = 0 at 6 angles: [-3°, 0°, 3°, 11°, 12°, 13°]

**Hypothesis**: Pressure drag (CDP) calculation returns zero in certain conditions

**Possible Causes**:
- Integration bounds wrong (missing portions of airfoil)
- Cp = 1 - Ue² produces zero when integrated (cancellation error)
- Conditional logic skips CDP calculation
- Wake-based CDP calculation fails

**Investigation**:
File: `crates/rustfoil-solver/src/viscous/forces.rs` or `viscal.rs`

Check:
- How is `cd_pressure` computed?
- When can it be exactly zero?
- Is there a fallback to wake-based CD?

### 3. **CD Anomaly at α=-7°**

CD = 0.0446 (should be 0.0021) → 2066% error

**Unique to this angle**, suggesting:
- Separation onset
- Wake march instability
- Transition prediction issue
- Newton convergence failure

**Investigation**: Run with debug and check:
```bash
cargo run --release -- viscous naca0012_xfoil_paneled.dat --alpha=-7 --re=3000000 --no-repanel --debug=debug_a-7.json
```

Check for:
- Separation events
- Large theta values at TE
- Wake station anomalies
- Convergence issues

### 4. **Low Alpha Scatter** (|α| ≤ 2°)

Errors range from -21% to +20% with large variance.

**Hypothesis**: Small numerical errors dominate at low α where forces are small

**Contributing factors**:
- Stagnation point positioning sensitivity
- Small Ue values amplify relative errors
- CL/CD approaching zero makes % errors large

## Investigation Priority

### Priority 1: CL Systematic High Bias (3°-15°)
**Impact**: Affects 80% of useful operating range
**Magnitude**: 9-16% error
**Approach**: Compare theta, Ue, and force integration at α=4°, 8°, 12°

### Priority 2: CD = 0 Cases  
**Impact**: 6 angles report zero pressure drag
**Magnitude**: CD should be O(0.001) but is 0
**Approach**: Debug CDP calculation logic, check integration bounds

### Priority 3: α=-7° CD Anomaly
**Impact**: Single angle outlier
**Magnitude**: 2000% error
**Approach**: Check for separation, wake instability, or convergence failure

### Priority 4: Low α Scatter
**Impact**: 5 angles with |α| ≤ 2°
**Magnitude**: High variance but small absolute errors
**Approach**: Improve numerical precision at small values

## Recommended Comparison Points

To identify root causes, compare at α=4° (good CL match, 9% error):

| Intermediate Value | Location | File | Method |
|-------------------|----------|------|--------|
| **Inviscid Ue** | All panels | debug trace | FULL_INVISCID event |
| **BL theta** | Stations 1, 10, 50, TE | debug trace | MRCHUE events |
| **BL delta_star** | Same stations | debug trace | MRCHUE events |
| **Mass defect** | Same stations | Computed | theta * delta_star * Ue |
| **DUI** | Stations 1, 10, 50 | debug trace | UESET event (if available) |
| **Ue after VI** | Same stations | debug trace | Final MRCHUE iteration |
| **Cp distribution** | All panels | Computed | 1 - Ue² |
| **Transition x/c** | Upper & lower | CLI output | x_tr values |

## Next Steps

1. **Generate detailed comparison at α=4°**:
   - Extract MRCHUE events from both RustFoil and XFOIL
   - Match stations by x/c
   - Compare theta, H, Ue, Cf at each station
   - Identify where divergence begins

2. **Debug CD=0 cases**:
   - Add debug output to CDP calculation
   - Check what triggers zero CDP
   - Verify integration limits and Cp values

3. **Fix identified issues**:
   - If theta is consistently low, adjust laminar/turbulent closures
   - If Ue is wrong, check DIJ or UESET calculation
   - If CDP calculation has bugs, fix the integration

4. **Rerun sweep** and verify improvements

## Tools Available

- `scripts/compare_bl_stations.py` - Station-by-station BL comparison
- `scripts/compare_traces.py` - Event-by-event trace comparison
- Debug traces in `/Users/harry/flexfoil-boundary-layer/debug_*.json`
- XFOIL traces in `traces/xfoil/naca0012/re3e06/alpha_*.json`
