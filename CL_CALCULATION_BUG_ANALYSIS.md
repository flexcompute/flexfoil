# CL Calculation Bug Analysis: RustFoil vs XFOIL at α=4°

## Problem Summary
- **RustFoil CL**: 0.4655
- **XFOIL CL**: 0.4424
- **Error**: +5.2% (RustFoil overpredicts)
- **Ue match**: Perfect (0.000003% error)

Since Ue matches perfectly, the CL error must come from the CL calculation chain.

## XFOIL's CLCALC Formula (xfoil.f line ~1088-1168)

```fortran
DO I=1, N
  IP = I+1
  IF(I.EQ.N) IP = 1
  
  CGINC = 1.0 - (GAM(IP)/QINF)**2
  CPG2 = CGINC/(BETA + BFAC*CGINC)  ! Compressible Cp
  
  DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA  ! Wind-axis projection
  
  AG = 0.5*(CPG2 + CPG1)  ! Average Cp
  
  CL = CL + DX * AG  ! Accumulate lift
ENDDO
```

For incompressible flow (MINF=0, QINF=1):
- `CGINC = 1.0 - GAM²`
- `CPG = CGINC` (since BETA=1, BFAC=0)
- So: **Cp = 1 - gamma²**

## RustFoil's clcalc Implementation (viscal.rs line 60-99)

```rust
pub fn clcalc(panel_x: &[f64], panel_y: &[f64], gamma: &[f64], alpha: f64) -> f64 {
    let ca = alpha.cos();
    let sa = alpha.sin();
    
    let mut cpg1 = 1.0 - gamma[0] * gamma[0];  // ✓ Correct Cp formula
    let mut cl = 0.0;
    
    for i in 0..n {
        let ip = (i + 1) % n;
        let cpg2 = 1.0 - gamma[ip] * gamma[ip];  // ✓ Correct Cp formula
        
        let dx_wind = (panel_x[ip] - panel_x[i]) * ca + (panel_y[ip] - panel_y[i]) * sa;  // ✓ Correct wind-axis projection
        
        let cp_avg = 0.5 * (cpg1 + cpg2);  // ✓ Correct averaging
        
        cl += cp_avg * dx_wind;  // ✓ Correct integration
        
        cpg1 = cpg2;
    }
    
    cl
}
```

**The clcalc function itself is CORRECT** - it matches XFOIL exactly.

## The Bug: construct_gamma_from_stations (viscal.rs line 113-171)

The bug is in how gamma is constructed from BL station edge velocities:

### Current (WRONG) Implementation:
```rust
// Upper surface: VTI = +1
gamma[idx] = station.u.abs();  // ❌ BUG: Using .abs()

// Lower surface: VTI = -1  
gamma[idx] = -station.u.abs();  // ❌ BUG: Using .abs()
```

### XFOIL's Correct Formula:
According to XFOIL's QVFUE (xpanel.f:1597) and GAMQV (xpanel.f:1633):
```fortran
QVIS(I) = VTI(IBL,IS) * UEDG(IBL,IS)
GAM(I) = QVIS(I)
```

So:
- **Upper surface (VTI=+1)**: `gamma = +1 * Ue = +Ue`
- **Lower surface (VTI=-1)**: `gamma = -1 * Ue = -Ue`

### The Issue:
If `station.u` is always positive (magnitude), then:
- Upper: `gamma = +|Ue|` ✓ (correct if Ue > 0)
- Lower: `gamma = -|Ue|` ✓ (correct if Ue > 0)

But if `station.u` can be negative OR if there's a sign convention issue, then using `.abs()` loses information.

### Correct Fix:
```rust
// Upper surface: VTI = +1
gamma[idx] = station.u;  // Use signed value directly

// Lower surface: VTI = -1
gamma[idx] = -station.u;  // Flip sign for lower surface
```

**However**, if `station.u` is always stored as positive magnitude, then the current code should work. The bug might be elsewhere.

## Alternative Bug: Missing Panels or Integration Weights

If gamma construction is correct, the bug could be:

1. **Missing panels**: Some panels not included in the sum
2. **Wrong integration weights**: dx_wind computed incorrectly
3. **Sign error**: Wrong sign in accumulation

## Investigation Steps

1. **Check gamma values**: Compare RustFoil gamma array with XFOIL GAM array
2. **Check Cp values**: Compare RustFoil Cp = 1 - gamma² with XFOIL CPG
3. **Check dx_wind**: Compare wind-axis projections panel-by-panel
4. **Check integration**: Verify all panels are included in the sum

## Expected Fix

The most likely bug is in `construct_gamma_from_stations`:
- Remove `.abs()` calls
- Use signed `station.u` directly with VTI sign convention
- Upper: `gamma = +station.u`
- Lower: `gamma = -station.u`

This ensures gamma has the correct sign for Cp = 1 - gamma² calculation.

**Note**: If `station.u` is always stored as positive magnitude, then removing `.abs()` won't change the result. However, it's still the correct implementation according to XFOIL's formula: `GAM = VTI * UEDG`.

## Additional Checks

1. **Unfilled panels**: Check if all panels have corresponding BL stations. Unfilled panels will have gamma=0, Cp=1, which may contribute incorrectly to CL.

2. **Panel ordering**: Verify that panel indices match XFOIL's ordering (TE → upper → LE → lower → TE).

3. **Sign consistency**: Ensure that the sign convention for upper/lower surfaces matches XFOIL's VTI convention throughout the codebase.

## Line-by-Line Comparison

### XFOIL CLCALC (xfoil.f:1088-1168)
```fortran
DO I=1, N
  IP = I+1
  IF(I.EQ.N) IP = 1
  
  CGINC = 1.0 - (GAM(IP)/QINF)**2      ! Cp = 1 - gamma²
  CPG2 = CGINC/(BETA + BFAC*CGINC)     ! Compressible correction
  
  DX = (X(IP) - X(I))*CA + (Y(IP) - Y(I))*SA  ! Wind-axis projection
  
  AG = 0.5*(CPG2 + CPG1)                ! Average Cp
  
  CL = CL + DX * AG                     ! Accumulate
ENDDO
```

### RustFoil clcalc (viscal.rs:60-99)
```rust
for i in 0..n {
    let ip = (i + 1) % n;  // ✓ Wrap-around correct
    
    let cpg2 = 1.0 - gamma[ip] * gamma[ip];  // ✓ Cp formula correct
    
    let dx_wind = (panel_x[ip] - panel_x[i]) * ca + (panel_y[ip] - panel_y[i]) * sa;  // ✓ Wind-axis correct
    
    let cp_avg = 0.5 * (cpg1 + cpg2);  // ✓ Averaging correct
    
    cl += cp_avg * dx_wind;  // ✓ Integration correct
}
```

**The clcalc function matches XFOIL exactly.** The bug must be in `construct_gamma_from_stations`.
