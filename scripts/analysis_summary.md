# RustFoil vs XFOIL Systematic Comparison Analysis

## Key Findings

### 1. Edge Velocity (Ue) - Generally Good
- **Match within 1-3%** at most stations
- **Exception at TE**: RustFoil Ue ~5% higher (0.80 vs 0.77)
- This doesn't explain the large CD errors

### 2. Momentum Thickness (θ) - Critical Issue!
| Alpha | Location | RF θ | XF θ | Error | Impact |
|-------|----------|------|------|-------|--------|
| 4° | Upper TE | 3.99e-3 | 4.76e-3 | **-16%** | Lowers CD |
| 4° | Lower TE | 1.43e-3 | 8.02e-3 | **-82%** | Lowers CD |
| 6° | Upper TE | 2.76e-3 | 7.69e-3 | **-64%** | Lowers CD |
| 6° | Lower TE | 1.30e-3 | 9.01e-3 | **-86%** | Lowers CD |

**θ at TE is 15-86% LOW in RustFoil!**
Since CD ∝ 2θ_TE (Squire-Young), this directly explains the CD being 25-55% low.

### 3. Shape Factor (H) - Severe Issues
At certain stations, H spikes to unrealistic values:
- Alpha=0°, x=0.7: RF H=3.2 vs XF H=1.7 (**+90%**)
- Alpha=6°, x=0.9: RF H=7.8 vs XF H=1.5 (**+400%!**)

H > 4 indicates laminar separation, but XFOIL shows turbulent values (H ≈ 1.5).
This suggests **transition is happening too late or not at all** in RustFoil.

### 4. Skin Friction (Cf) - Mixed
- Near LE: Generally good match
- Mid-chord: 10-30% low (consistent with H issues)
- Near TE/separation: **Negative Cf** in RustFoil where XFOIL has positive

Negative Cf indicates RustFoil thinks flow is separated when it isn't.

## Error Trend by Alpha

| Alpha | CL Error | CD Error | Root Cause |
|-------|----------|----------|------------|
| 0° | N/A (asymmetry) | -25% | Stagnation discretization |
| 2° | -22% | -30% | Newton not converging properly |
| 4° | **-5%** | -36% | θ_TE too low |
| 6° | **+0.2%** | -55% | θ_TE much too low |
| 8° | +10% | +9% | Different convergence state |

## Primary Issues to Fix

### Issue 1: θ at TE is too low
**Why**: The march isn't accumulating enough momentum thickness by the trailing edge.
**Check**: Compare θ at intermediate stations - are we losing θ somewhere?

At x=0.1: θ matches within 10% (good)
At x=0.5: θ matches within 10% (good)  
At TE: θ is 15-86% LOW (bad!)

The θ is being lost somewhere between x=0.5 and TE.

### Issue 2: H spikes indicate wrong flow regime
**Why**: Flow appears to stay laminar (H > 3) where it should be turbulent (H ≈ 1.5)
**Check**: Transition locations

From the output at α=4°:
- Upper: x_tr = 0.21 (reasonable for α=4°)
- But H at x=0.3 is 3.1 (should be ~1.7 if turbulent)

The transition happens at x=0.21, but H doesn't drop - something wrong with turbulent closure.

### Issue 3: Lower surface TE θ is catastrophically low
At α=4°:
- RF θ_lower_TE = 1.43e-3
- XF θ_lower_TE = 8.02e-3
- Error: **-82%!**

The lower surface isn't growing θ properly near the TE.

## Recommendations

1. **Debug θ growth**: Add logging to see θ at every station, compare with XFOIL
2. **Check turbulent closures**: After transition, H should drop to ~1.5, but we see H > 3
3. **Check Cf computation**: Negative Cf shouldn't occur in attached flow
4. **Investigate wake**: XFOIL's TE station (side 2 in wake) has x > 1, may need to adjust θ_TE extraction
