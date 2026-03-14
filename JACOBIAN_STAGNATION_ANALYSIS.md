# Jacobian Matrix Stagnation Analysis

## Problem
Residual stagnates at ~0.64 in RustFoil's Newton iteration, indicating inaccurate Jacobian computation.

## Root Cause Identified

### Bug in `build_vm_global()` - Line 626

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:626`

**Current (WRONG) code**:
```rust
let u1_m_j = if ibl_i > 1 && panel_im1 < dij.nrows() && panel_j < dij.ncols() {
    -self.vti[iv - 1].max(self.vti[iv]) * self.vti[jv] * dij[(panel_im1, panel_j)]
```

**Problem**: The `.max(self.vti[iv])` is incorrect. It should use the VTI of the upstream station directly.

**Correct code**:
```rust
let u1_m_j = if ibl_i > 1 && panel_im1 < dij.nrows() && panel_j < dij.ncols() {
    -self.vti[iv - 1] * self.vti[jv] * dij[(panel_im1, panel_j)]
```

**Why this matters**:
- XFOIL's formula: `U1_M(JV) = -VTI(IBL-1,IS)*VTI(JBL,JS)*DIJ(I-1,J)`
- The upstream station's VTI must be used, not a max of both stations
- This bug causes incorrect VM matrix entries, leading to inaccurate Jacobian
- An inaccurate Jacobian causes Newton iteration to stagnate

## Verification

Compare with:
1. **XFOIL reference** (`Xfoil-instrumented/src/xbl.f:206`):
   ```fortran
   U2_M(JV) = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
   ```
   Then at end of loop (line 481): `U1_M(JV) = U2_M(JV)` for next iteration

2. **Single-surface version** (`crates/rustfoil-coupling/src/newton.rs:691`):
   ```rust
   -self.vti[i - 1] * self.vti[j] * dij[(panel_im1, panel_j)]
   ```
   This correctly uses `vti[i - 1]` directly.

## Impact

- **VM matrix entries are wrong** for all stations where `ibl_i > 1`
- **Jacobian accuracy degraded**, especially for cross-surface coupling
- **Newton convergence stalls** because updates are computed from wrong derivatives
- **Residual stagnates** at ~0.64 instead of converging to zero

## Other Checks Performed

1. ✅ **VA/VB assembly**: Correctly extracts columns 0,1 from VS1/VS2
2. ✅ **VS1/VS2 storage**: Correctly stores delta_star and Ue derivatives
3. ✅ **VM loop structure**: Correctly loops over both surfaces (jv from 1 to nsys)
4. ✅ **Forced changes**: Correctly adds DUE and DDS to residuals
5. ✅ **Cross-surface coupling**: VM matrix includes both surfaces correctly
6. ❌ **U1_M computation**: BUG - uses `.max()` instead of direct VTI value

## Fix

Replace line 626 in `global_newton.rs`:
```rust
// BEFORE (WRONG):
-self.vti[iv - 1].max(self.vti[iv]) * self.vti[jv] * dij[(panel_im1, panel_j)]

// AFTER (CORRECT):
-self.vti[iv - 1] * self.vti[jv] * dij[(panel_im1, panel_j)]
```

This matches XFOIL's behavior and the single-surface implementation.
