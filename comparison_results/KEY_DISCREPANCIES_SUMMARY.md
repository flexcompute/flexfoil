# Key Intermediate Value Discrepancies: RustFoil vs XFOIL

## Executive Summary

After implementing wake DIJ extension and multiple fixes, RustFoil achieves:
- **CL RMS Error**: 12.7% (was 24%)
- **Best Region**: 3°-8° with 9-10% error
- **Remaining Issues**: Systematic CL high bias, CD calculation problems

## Intermediate Value Analysis

### 1. **Boundary Layer Thickness (theta)**

Comparing theta at key locations shows a **complex pattern**:

#### α=4° (Good CL Match - 9% Error)
| x/c | RustFoil theta | XFOIL theta | % Diff | Impact |
|-----|----------------|-------------|--------|---------|
| 0.01 | 1.834e-5 | 1.849e-5 | -0.8% | ✓ Excellent |
| 0.10 | 1.042e-4 | 9.569e-5 | **+8.8%** | RustFoil higher |
| 0.50 | 9.477e-4 | 9.245e-4 | +2.5% | RustFoil slightly higher |
| 0.90 | 2.315e-3 | 2.255e-3 | +2.6% | RustFoil slightly higher |
| 1.00 | 3.248e-3 | 3.043e-3 | **+6.8%** | RustFoil higher |

**Pattern**: RustFoil theta is 2.5-8.8% **higher** at most stations.

**Impact on CL**:
- Higher theta → higher mass defect → lower edge velocity → **should reduce lift**
- But CL is 9% too HIGH, not low!
- **Contradiction**: Suggests another factor dominates

#### α=8° (Moderate Error - 16%)
| x/c | RustFoil theta | XFOIL theta | % Diff |
|-----|----------------|-------------|--------|
| 0.01 | 2.460e-5 | 2.747e-5 | **-10.4%** | RustFoil lower |
| 0.10 | 1.598e-4 | 1.685e-4 | -5.1% | RustFoil lower |
| 0.50 | 1.538e-3 | 1.576e-3 | -2.4% | RustFoil lower |
| 1.00 | 5.005e-3 | 4.892e-3 | +2.3% | RustFoil higher |

**Pattern**: RustFoil theta is **lower** early (by 5-10%), **higher** late (by 2%).

#### α=12° (Large Error - 15%)
| x/c | RustFoil theta | XFOIL theta | % Diff |
|-----|----------------|-------------|--------|
| 0.01 | 3.315e-5 | 4.292e-5 | **-22.8%** | RustFoil much lower |
| 0.10 | 1.621e-4 | 1.606e-4 | +0.9% | Match |
| 0.50 | 2.192e-3 | 2.250e-3 | -2.6% | RustFoil lower |
| 1.00 | 7.066e-3 | 7.056e-3 | +0.1% | Perfect match |

**Pattern**: RustFoil theta is **much lower** at stagnation (23%), but matches well downstream.

### 2. **Edge Velocity (Ue)**

Edge velocities match XFOIL **perfectly** at all stations (< 0.01% error).

**Conclusion**: The DIJ/UESET coupling is working correctly. The VI iteration converges to the same Ue as XFOIL.

### 3. **CL Discrepancy Root Cause**

Given that:
- Ue matches perfectly (VI coupling correct)
- theta is 2-9% higher at most stations (should reduce lift)
- But CL is 9% too HIGH

**Hypothesis**: The force calculation (CL integration) has a systematic bias.

**Files to check**:
- `crates/rustfoil-solver/src/viscous/viscal.rs` lines 860-890 (CL from circulation)
- `crates/rustfoil-solver/src/viscous/forces.rs` (force integration)

**Possible causes**:
1. **Integration bounds**: Missing some portion of airfoil
2. **Circulation formula**: `CL = ∫(γ_upper - γ_lower) ds` may have sign error
3. **Arc length calculation**: `ds` values may be wrong
4. **REYBL tuning**: The re/3 factor may need adjustment to re/4

### 4. **CD = 0 Problem**

Six angles report CD = 0: [-3°, 0°, 3°, 11°, 12°, 13°]

**Investigation needed**:
File: `crates/rustfoil-solver/src/viscous/viscal.rs` or `forces.rs`

Check:
- How is `cd_pressure` computed?
- When can it return exactly zero?
- Is there a conditional that skips calculation?

Example from trace:
```
α=3°: CL=0.4077, CD=0.00000
α=4°: CL=0.4833, CD=0.00116
```

The transition from CD=0 to CD>0 between 3° and 4° suggests a threshold or conditional.

### 5. **CD Anomalies**

#### α=-7°: CD = 0.0446 (should be 0.0021)
- 2066% error - clearly an outlier
- Likely separation, wake instability, or Newton divergence
- Check debug trace for:
  - Separation events
  - Large theta at TE (>0.1)
  - Non-convergence

#### α=-4°: CD = 0.00133 (expected 0.0008)  
- 57.6% error - reasonable after fix
- Was 0.592 before theta bounds fix
- Now in acceptable range

## Summary of Discrepancies

| Component | Status | Error | Priority |
|-----------|--------|-------|----------|
| **DIJ Matrix** | ✅ Fixed | Wake magnitudes correct | Complete |
| **Stagnation Init** | ✅ Fixed | Thwaites w/ REYBL=re/3 | Complete |
| **Edge Velocity (Ue)** | ✅ Correct | <0.01% error | Working |
| **BL Thickness (theta)** | ⚠️ Mixed | Early: 0-23% low, Late: 2-9% high | Investigate |
| **CL Integration** | ⚠️ High Bias | 9-16% systematic high | **FIX NEEDED** |
| **CD Calculation** | ❌ Broken | CD=0 at many angles | **FIX URGENT** |
| **Transition** | ⚠️ Early | 5% earlier than XFOIL | Minor impact |

## Recommended Actions

### Priority 1: Fix CD = 0 Bug
**Impact**: Makes CD unusable at 20% of angles  
**Approach**: Debug `cd_pressure` calculation, find why it returns 0  
**Expected result**: CD should be O(0.001) at all angles

### Priority 2: Investigate CL Integration Bias
**Impact**: 9-16% systematic CL over-prediction  
**Approach**: 
- Compare Cp distributions between RustFoil and XFOIL
- Verify circulation integration formula and signs
- Check if REYBL=re/3.5 or re/4 improves match

### Priority 3: Fix Early BL Development
**Impact**: theta 10-23% low at x<0.01 for higher α  
**Approach**:
- Verify Thwaites formula at different α
- Check if stagnation Ue extraction is correct
- Compare laminar closure coefficients

## Files Requiring Investigation

1. `crates/rustfoil-solver/src/viscous/forces.rs` - CD calculation (Priority 1)
2. `crates/rustfoil-solver/src/viscous/viscal.rs` - CL integration (Priority 2)
3. `crates/rustfoil-bl/src/state.rs` - Stagnation initialization (Priority 3)
4. `crates/rustfoil-coupling/src/march.rs` - Early BL marching (Priority 3)

## Verification Tests

After fixes, verify:
- [ ] CD > 0 at all angles (no more zeros)
- [ ] CL at α=4-8° within 5% of XFOIL
- [ ] theta at x=0.01 matches XFOIL within 10%
- [ ] No CD anomalies > 100% error
