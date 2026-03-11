# Newton System Matrix Comparison: XFOIL vs RustFoil

## Test Case
- Airfoil: NACA 0012
- α = 4°
- Re = 3×10⁶
- Mach = 0.0
- Ncrit = 9.0

## Aerodynamic Results Comparison

| Metric | XFOIL | RustFoil | Difference |
|--------|-------|----------|------------|
| CL | 0.4424 | 0.4161 | -5.9% |
| CD | 0.00618 | 0.00763 | +23.5% |
| CM | 0.0014 | 0.0000 | - |
| x_tr (Upper) | 0.1475 | 0.1376 | -6.7% |
| x_tr (Lower) | 0.8704 | 0.8384 | -3.7% |
| Converged | Yes (5 iter) | Yes (30 iter) | - |

## Matrix Comparison (Iteration 1)

### VA Matrix (Downstream Jacobian) - EXCELLENT MATCH

| Station | Element | XFOIL | RustFoil | Diff % |
|---------|---------|-------|----------|--------|
| SIMI (X:2/R:1) | [1][1] | -2.852e5 | -2.854e5 | 0.1% |
| SIMI (X:2/R:1) | [2][1] | 3.788e5 | 3.791e5 | 0.1% |
| X:3 / R:2 | [1][1] | -1.170e5 | -1.160e5 | 0.9% |
| X:3 / R:2 | [2][1] | 2.750e5 | 2.730e5 | 0.7% |
| X:4 / R:3 | [1][1] | -3.520e4 | -3.490e4 | 0.9% |
| X:5 / R:4 | [1][1] | -6.330e3 | -6.330e3 | 0.0% |

### VB Matrix (Upstream Jacobian) - EXCELLENT MATCH

| Station | Element | XFOIL | RustFoil | Diff % |
|---------|---------|-------|----------|--------|
| X:3 / R:2 | [0][0] | -1.0 | -1.0 | 0.0% |
| X:3 / R:2 | [1][1] | -2.370e5 | -2.350e5 | 0.8% |
| X:3 / R:2 | [2][1] | 2.070e5 | 2.050e5 | 1.0% |
| X:4 / R:3 | [1][1] | -1.460e5 | -1.460e5 | 0.0% |
| X:5 / R:4 | [1][1] | -1.160e5 | -1.160e5 | 0.0% |

### VDEL Residuals (After Forced Changes) - SIMILAR PATTERNS

| Station | Element | XFOIL | RustFoil | Comment |
|---------|---------|-------|----------|---------|
| SIMI | [1] res_mom | -42.3 | -33.1 | 22% smaller |
| SIMI | [2] res_shape | 31.5 | 27.7 | 12% smaller |
| X:3/R:2 | [1] res_mom | -32.7 | -34.3 | 5% larger |
| X:3/R:2 | [2] res_shape | 23.5 | 22.9 | 2% smaller |

## Key Finding: SIMI Combining Bug Fixed

### The Bug
XFOIL combines the Jacobian matrices at SIMI (similarity) stations:
```fortran
VS2(k,j) = VS1(k,j) + VS2(k,j)
VS1(k,j) = 0
```

This **doubles** the Jacobian because at SIMI, both upstream (VS1) and downstream (VS2) 
refer to the same station, so their derivatives are identical.

### Before Fix
- RustFoil VA[1][1] at SIMI = -1.427e5 (half of XFOIL)
- Newton iteration: **DIVERGENT**

### After Fix
- RustFoil VA[1][1] at SIMI = -2.854e5 (matches XFOIL)
- Newton iteration: **CONVERGENT**

## Implementation Change

In `global_newton.rs`, line ~647:

```rust
// Before (WRONG):
self.va[iv][eq][1] = jacobian.vs2[eq][1];

// After (CORRECT):
self.va[iv][eq][1] = jacobian.vs1[eq][1] + jacobian.vs2[eq][1];
```

## Remaining Differences

The 23% CD overprediction is likely due to:

1. **Forced Change (DUE) Computation**: The VDEL residuals differ by ~20% at SIMI, 
   suggesting differences in how edge velocity mismatches are computed

2. **March Initialization**: RustFoil may use slightly different initial BL states 
   from the march phase

3. **SIMI Residuals**: We commented out zeroing SIMI residuals, but XFOIL may handle
   them differently

4. **Convergence Tolerance**: RustFoil takes 30 iterations vs XFOIL's 5, suggesting
   relaxation or tolerance differences

## Next Steps

1. Compare DUE (edge velocity mismatch) computation between XFOIL and RustFoil
2. Verify march initialization produces same theta/delta_star as XFOIL
3. Check if SIMI residual handling affects convergence rate
4. Compare relaxation factor computation (RMSBL formula)
