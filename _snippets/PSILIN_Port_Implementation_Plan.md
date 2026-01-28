---
created: 2026-01-28
project_source: "flexfoil-boundary-layer"
tags: [status/active, type/implementation-plan, topic/aerodynamics]
---

# PSILIN Port Implementation Plan

## Context

The RustFoil viscous solver has a 32% DUI (mass defect velocity correction) discrepancy compared to XFOIL. Root cause analysis traced this to the **DQDM computation** in the wake-airfoil block of the DIJ matrix.

### The Problem

RustFoil's `compute_source_tangent_influence()` uses a **simple point-source approximation**:

```rust
// Current (WRONG) implementation
let factor = ds / (2.0 * PI * r2);
dqdm[j] = factor * (dx * tx + dy * ty);
```

**Result:** DQDM values are ~100x smaller with wrong signs compared to XFOIL.

| Panel | XFOIL DQDM | RustFoil DQDM | Ratio |
|-------|------------|---------------|-------|
| 0 | 32.3 | 0.08 | 0.003x |
| 1 | -20.5 | 0.20 | -0.01x |
| 2 | -5.71 | 0.10 | -0.02x |

### Inputs Verified Correct

- Wake point coordinates: ✅ Match
- Wake point normals: ✅ Match
- Panel count (n=160): ✅ Match
- Panel coordinates: ✅ Match to 8 decimal places

**Conclusion:** The algorithm is wrong, not the inputs.

---

## XFOIL's PSILIN Algorithm

XFOIL computes DQDM using full **panel source integrals** with:

1. **Half-panel decomposition** (1-0 and 0-2 segments)
2. **Quadratic source distribution** (SSUM, SDIF terms)
3. **Proper geometric derivatives** (PSNI, PDNI)
4. **Three-point stencil** (JM, JO, JP for 1-0; JO, JP, JQ for 0-2)

### Key XFOIL Variables

| Variable | Description |
|----------|-------------|
| `X1, X2, YY` | Local panel coordinates (field point relative to panel) |
| `G1, G2, G0` | log(r²) at panel endpoints and midpoint |
| `T1, T2, T0` | atan2 angles at panel endpoints and midpoint |
| `DSIO, DSIM, DSIP` | Inverse panel lengths |
| `PSUM, PDIF` | Symmetric/antisymmetric streamfunction contributions |
| `PSNI, PDNI` | Normal derivatives of PSUM, PDIF |
| `DQDM[j]` | Output: dQtan/dSigma for panel j |

### Algorithm Flow (SIGLIN=.TRUE.)

```
For each panel JO = 1 to N-1:
    JP = JO + 1
    JM = JO - 1  (or JO if JO=1)
    JQ = JP + 1  (or JP if JO=N-1)
    
    1. Compute local coordinates (X1, X2, YY)
    2. Compute log/atan terms (G1, G2, T1, T2)
    3. Set up midpoint (X0, G0, T0)
    
    // Half-panel 1-0:
    4. PSUM = X0*(T0-APAN) - X1*(T1-APAN) + 0.5*YY*(G1-G0)
    5. PDIF = quadratic correction
    6. PSNI, PDNI = normal derivatives
    7. DQDM[JM] += QOPI*(-PSNI*DSIM + PDNI*DSIM)
       DQDM[JO] += QOPI*(-PSNI*DSIO - PDNI*DSIO)
       DQDM[JP] += QOPI*( PSNI*(DSIO+DSIM) + PDNI*(DSIO-DSIM))
    
    // Half-panel 0-2:
    8. Similar computation with X0, X2, G0, G2, T0, T2
    9. DQDM[JO] += ...
       DQDM[JP] += ...
       DQDM[JQ] += ...
```

---

## Implementation Status

| Phase | Description | Status |
|-------|-------------|--------|
| 1 | Create `psilin_with_dqdm()` function | ✅ COMPLETE |
| 2 | Port core panel integral | ✅ COMPLETE |
| 3 | Replace in `build_dij_with_wake()` | ✅ COMPLETE |
| 4 | Add indirect effects (CIJ coupling) | ⏳ PENDING |

### Verification Results (2026-01-28)

DQDM values now match XFOIL:
- Wake point 161: Max rel error **0.73%**, all 20 values within 5%
- Wake point 162: Max rel error **0.19%**, all 20 values within 5%

**Before fix:** Values were ~100x smaller with wrong signs (ratios 0.003 to -0.13)
**After fix:** Ratios ~1.0000, avg error 0.04%

**Remaining issue:** CD still 23% overpredicted. May be due to missing indirect effects or other issues.

---

## Implementation Plan

### Phase 1: Create `psilin_with_siglin()` Function [COMPLETE]

**File:** `crates/rustfoil-inviscid/src/influence.rs`

**Signature:**
```rust
pub struct PsilinResult {
    pub psi: f64,           // Streamfunction
    pub psi_ni: f64,        // dPsi/dni
    pub dzdg: Vec<f64>,     // dPsi/dGamma
    pub dzdm: Vec<f64>,     // dPsi/dSigma (mass source)
    pub dqdg: Vec<f64>,     // dQtan/dGamma
    pub dqdm: Vec<f64>,     // dQtan/dSigma (THE KEY OUTPUT)
}

pub fn psilin_with_siglin(
    geom: &AirfoilGeometry,
    field_idx: usize,       // I in XFOIL (can be > N for wake)
    xi: f64, yi: f64,       // Field point
    nxi: f64, nyi: f64,     // Field point normal
    gam: &[f64],            // Vorticity distribution (optional)
    sig: &[f64],            // Source distribution (optional)
) -> PsilinResult
```

### Phase 2: Port Core Panel Integral

**Steps:**

1. **Local coordinate transformation:**
   ```rust
   let sx = (x[jp] - x[jo]) * dsio;
   let sy = (y[jp] - y[jo]) * dsio;
   let x1 = sx * rx1 + sy * ry1;  // tangent projection
   let x2 = sx * rx2 + sy * ry2;
   let yy = sx * ry1 - sy * rx1;  // normal projection
   ```

2. **Log/atan computation with reflection handling:**
   ```rust
   let sgn = if io >= 1 && io <= n { 1.0 } else { yy.signum() };
   let g1 = if io != jo && rs1 > 0.0 { rs1.ln() } else { 0.0 };
   let t1 = if io != jo && rs1 > 0.0 {
       (sgn * x1).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI
   } else { 0.0 };
   ```

3. **Half-panel 1-0 contribution:**
   ```rust
   let x0 = 0.5 * (x1 + x2);
   let rs0 = x0 * x0 + yy * yy;
   let g0 = rs0.ln();
   let t0 = (sgn * x0).atan2(sgn * yy) + (0.5 - 0.5 * sgn) * PI;
   
   let dxinv = 1.0 / (x1 - x0);
   let psum = x0 * (t0 - apan) - x1 * (t1 - apan) + 0.5 * yy * (g1 - g0);
   let pdif = ((x1 + x0) * psum + rs1 * (t1 - apan) - rs0 * (t0 - apan)
              + (x0 - x1) * yy) * dxinv;
   
   // Normal derivatives
   let psx1 = -(t1 - apan);
   let psx0 = t0 - apan;
   let psyy = 0.5 * (g1 - g0);
   
   let pdx1 = ((x1 + x0) * psx1 + psum + 2.0 * x1 * (t1 - apan) - pdif) * dxinv;
   let pdx0 = ((x1 + x0) * psx0 + psum - 2.0 * x0 * (t0 - apan) + pdif) * dxinv;
   let pdyy = ((x1 + x0) * psyy + 2.0 * (x0 - x1 + yy * (t1 - t0))) * dxinv;
   
   let psni = psx1 * x1i + psx0 * (x1i + x2i) * 0.5 + psyy * yyi;
   let pdni = pdx1 * x1i + pdx0 * (x1i + x2i) * 0.5 + pdyy * yyi;
   
   // Accumulate DQDM
   dqdm[jm] += QOPI * (-psni * dsim + pdni * dsim);
   dqdm[jo] += QOPI * (-psni * dsio - pdni * dsio);
   dqdm[jp] += QOPI * (psni * (dsio + dsim) + pdni * (dsio - dsim));
   ```

4. **Half-panel 0-2 contribution:** (similar structure)

### Phase 3: Replace in `build_dij_with_wake()`

**File:** `crates/rustfoil-inviscid/src/system.rs`

```rust
// BEFORE (wrong):
let dqdm_airfoil = compute_source_tangent_influence(&self.geom, xi, yi, tx, ty);

// AFTER (correct):
let result = psilin_with_siglin(
    &self.geom,
    i_global,  // wake point index
    xi, yi,
    nxi, nyi,
    &[], &[],  // no gam/sig needed, just computing influence coefficients
);
let dqdm_airfoil = result.dqdm;
```

### Phase 4: Add Indirect Effects

XFOIL adds indirect gamma effects to wake-airfoil DIJ (lines 1254-1276):

```fortran
C---- add on effect of all sources on airfoil vorticity which effects wake Qtan
DO 80 I=N+1, N+NW
  IW = I-N
  DO 810 J=1, N
    SUM = 0.
    DO 8100 K=1, N
      SUM = SUM + CIJ(IW,K)*DIJ(K,J)
    CONTINUE
    DIJ(I,J) = DIJ(I,J) + SUM
  CONTINUE
CONTINUE
```

This requires also computing DQDG (dQtan/dGamma) in psilin, then:

```rust
// CIJ[iw][k] = DQDG[k] for wake point iw
// DIJ[i][j] += sum_k(CIJ[iw][k] * DIJ[k][j])
for j in 0..n {
    let mut sum = 0.0;
    for k in 0..n {
        sum += cij[iw][k] * dij[(k, j)];
    }
    dij[(i_global, j)] += sum;
}
```

---

## Testing Strategy

### 1. Unit Test: DQDM Values

```rust
#[test]
fn test_dqdm_matches_xfoil() {
    // Load XFOIL reference from xfoil_dqdm.json
    // Compare first 20 DQDM values for wake point 1
    // Tolerance: 1e-6 relative error
}
```

### 2. Integration Test: DIJ Matrix

```rust
#[test]
fn test_dij_wake_airfoil_block() {
    // Compare DIJ[N+1, 1..N] between XFOIL and RustFoil
    // This is the wake-airfoil block
}
```

### 3. End-to-End: VDEL Residuals

```rust
#[test]
fn test_vdel_matches_xfoil() {
    // Run full viscous solve
    // Compare VDEL values at iteration 1
    // Should match within 5%
}
```

---

## Files to Modify

| File | Changes |
|------|---------|
| `crates/rustfoil-inviscid/src/influence.rs` | Add `psilin_with_siglin()` (~200 lines) |
| `crates/rustfoil-inviscid/src/system.rs` | Replace `compute_source_tangent_influence()`, add indirect effects |
| `crates/rustfoil-inviscid/src/lib.rs` | Export new function |

## Estimated Complexity

- **psilin_with_siglin():** ~200 lines (core panel integral math)
- **Indirect effects:** ~30 lines
- **Tests:** ~100 lines

**Total:** ~330 lines of new code

---

## Reference Files

- **XFOIL PSILIN:** `Xfoil-instrumented/src/xpanel.f` lines 99-420
- **XFOIL QDCALC:** `Xfoil-instrumented/src/xpanel.f` lines 1164-1290
- **Current RustFoil:** `crates/rustfoil-inviscid/src/system.rs` lines 550-584
- **Debug comparison:** `docs/debugging/newton_comparison_results.md`

---

## Related Notes

- [[Newton_BruteForce_Comparison_Plan]] - Original comparison methodology
- [[rustfoil-xfoil-methodology-comparison]] - High-level methodology differences

