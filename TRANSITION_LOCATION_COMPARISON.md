# Transition Location Comparison: RustFoil vs XFOIL

## Summary

Comparison of laminar-turbulent transition locations between RustFoil and XFOIL for NACA 0012 at Re = 3×10⁶.

## Comparison Table

| Alpha | RF x_tr_upper | XF x_tr_upper | Err% | RF x_tr_lower | XF x_tr_lower | Err% |
|-------|---------------|---------------|------|---------------|---------------|------|
| 0°    | 0.6648        | 0.5133        | 29.51| 0.6648        | 0.5133        | 29.51|
| 2°    | 0.3061        | 0.3212        | -4.73| 0.8186        | 0.7024        | 16.54|
| 4°    | 0.0065        | 0.1475        | -95.57| 0.8868        | 0.8704        | 1.88 |
| 6°    | 0.0039        | 0.0580        | -93.32| 0.9549        | 0.9684        | -1.40|
| 8°    | 0.0028        | 0.0281        | -89.99| 0.9719        | 0.9953        | -2.35|

## Key Findings

### Critical Issue: Premature Transition on Upper Surface

At higher angles of attack (α ≥ 4°), RustFoil transitions **much earlier** on the upper surface than XFOIL:

- **α = 4°**: RustFoil transitions at x = 0.0065 vs XFOIL at x = 0.1475 (**-95.57% error**)
- **α = 6°**: RustFoil transitions at x = 0.0039 vs XFOIL at x = 0.0580 (**-93.32% error**)
- **α = 8°**: RustFoil transitions at x = 0.0028 vs XFOIL at x = 0.0281 (**-89.99% error**)

RustFoil transitions almost immediately (x_tr ≈ 0.003-0.007) while XFOIL transitions much later (x_tr ≈ 0.03-0.15).

### Impact on Boundary Layer and Lift

**Earlier transition → longer turbulent region → thicker boundary layer → reduced lift**

1. **Thicker Boundary Layer**: 
   - RustFoil has a much longer turbulent region on the upper surface
   - Turbulent boundary layers are thicker than laminar ones
   - This explains the BL thickness over-prediction

2. **Reduced Suction (Lower Pressure)**:
   - Thicker BL reduces the effective camber/curvature
   - Lower pressure recovery on upper surface
   - Reduced pressure difference between upper and lower surfaces

3. **Lower Lift Coefficient**:
   - Reduced suction on upper surface directly reduces lift
   - This explains the CL under-prediction at higher angles

### Lower Surface Behavior

The lower surface shows much better agreement:
- Differences are typically < 2% for α ≥ 4°
- At α = 2°, RustFoil transitions slightly later (16.54% error)
- This suggests the transition model works better for favorable pressure gradients

### Zero Angle of Attack Anomaly

At α = 0°, RustFoil transitions **later** than XFOIL (29.51% error):
- RF: x_tr = 0.6648 (both surfaces)
- XF: x_tr = 0.5133 (both surfaces)
- This is opposite to the high-angle behavior
- May indicate different transition criteria or N-factor calculation

## Detailed Analysis by Angle

### α = 0°
- **CL = 0.0219**
- RustFoil transitions later on both surfaces
- More laminar flow → thinner BL → increased suction
- However, at zero angle, lift is minimal, so impact is small

### α = 2°
- **CL = 0.2334**
- Upper: RF transitions slightly earlier (0.3061 vs 0.3212, -4.73%)
- Lower: RF transitions later (0.8186 vs 0.7024, +16.54%)
- Small differences, minimal impact

### α = 4°
- **CL = 0.4655**
- Upper: RF transitions **much earlier** (0.0065 vs 0.1475, -95.57%)
  - RF has ~23× longer turbulent region
  - Significant impact on BL thickness and lift
- Lower: Good agreement (0.8868 vs 0.8704, +1.88%)

### α = 6°
- **CL = 0.7814**
- Upper: RF transitions **much earlier** (0.0039 vs 0.0580, -93.32%)
  - RF has ~15× longer turbulent region
- Lower: Good agreement (0.9549 vs 0.9684, -1.40%)

### α = 8°
- **CL = 1.1639**
- Upper: RF transitions **much earlier** (0.0028 vs 0.0281, -89.99%)
  - RF has ~10× longer turbulent region
- Lower: Good agreement (0.9719 vs 0.9953, -2.35%)

## Root Cause Hypothesis

The premature transition on the upper surface at higher angles suggests:

1. **Transition Model Issue**: 
   - N-factor amplification calculation may be incorrect
   - Transition criteria (Ncrit = 9.0) may be applied incorrectly
   - Amplification rate integration may have errors

2. **Pressure Gradient Sensitivity**:
   - Transition model may be too sensitive to adverse pressure gradients
   - At higher angles, upper surface has stronger adverse gradient
   - May trigger transition too early

3. **Edge Velocity Calculation**:
   - If edge velocities are incorrect, amplification rates will be wrong
   - This could cause premature transition detection

## Recommendations

1. **Investigate Transition Model**:
   - Compare N-factor evolution between RustFoil and XFOIL
   - Check amplification rate calculations
   - Verify transition criteria application

2. **Check Edge Velocities**:
   - Compare edge velocity distributions at transition locations
   - Verify inviscid solution accuracy

3. **Examine Transition Detection Logic**:
   - Review how transition is detected in RustFoil
   - Compare with XFOIL's TRANSITION subroutine behavior

4. **Fix Transition Model**:
   - Once root cause is identified, fix the premature transition
   - This should improve both BL thickness and CL predictions

## Conclusion

The premature transition on the upper surface at higher angles of attack is a **critical issue** that explains:
- ✅ BL thickness over-prediction
- ✅ CL under-prediction
- ✅ Reduced pressure recovery on upper surface

Fixing the transition model should significantly improve agreement with XFOIL.
