# UESET Comparison: XFOIL vs RustFoil

## Overview

This document compares XFOIL's `UESET` subroutine with RustFoil's `compute_ue_from_mass_both` function. Both implement the viscous-inviscid coupling that computes edge velocities from mass defect using the DIJ influence matrix.

## XFOIL Implementation (xpanel.f)

### Source Code

```fortran
SUBROUTINE UESET
C---------------------------------------------------------
C     Sets Ue from inviscid Ue plus all source influence
C---------------------------------------------------------
      INCLUDE 'XFOIL.INC'

      DO 1 IS=1, 2
        DO 10 IBL=2, NBL(IS)
          I = IPAN(IBL,IS)

          DUI = 0.
          DO 100 JS=1, 2
            DO 1000 JBL=2, NBL(JS)
              J  = IPAN(JBL,JS)
              UE_M = -VTI(IBL,IS)*VTI(JBL,JS)*DIJ(I,J)
              DUI = DUI + UE_M*MASS(JBL,JS)
 1000       CONTINUE
  100     CONTINUE

          UEDG(IBL,IS) = UINV(IBL,IS) + DUI

   10   CONTINUE
    1 CONTINUE

      RETURN
      END
```

### Formula

```
UEDG(IBL,IS) = UINV(IBL,IS) + Σ_JS Σ_JBL [ -VTI(IBL,IS) * VTI(JBL,JS) * DIJ(I,J) * MASS(JBL,JS) ]
```

Where:
- `I = IPAN(IBL,IS)` - panel index for current BL station
- `J = IPAN(JBL,JS)` - panel index for influence station
- `UINV` - inviscid edge velocity
- `MASS` - mass defect (δ* × Ue)
- `DIJ` - source influence matrix

### VTI Sign Convention

Set in `IBLPAN` subroutine (xpanel.f):

| Surface | IS | VTI Value |
|---------|-----|-----------|
| Upper (top) | 1 | +1.0 |
| Lower (bottom) | 2 | -1.0 |
| Wake (both) | 1,2 | -1.0 (lower), +1.0 (upper plotting) |

### Station Ranges

- **Outer loop:** `IBL = 2 to NBL(IS)` — **skips stagnation (IBL=1)**
- **Inner loop:** `JBL = 2 to NBL(JS)` — **skips stagnation (JBL=1)**
- **Surfaces:** Both IS=1 (upper) and IS=2 (lower) are summed

### Panel Index Mapping (IPAN)

From `IBLPAN`:
- **Upper surface (IS=1):** Stations go from stagnation (IST) toward panel 1, so `IPAN(IBL,1) = IST, IST-1, IST-2, ...`
- **Lower surface (IS=2):** Stations go from stagnation (IST+1) toward panel N, so `IPAN(IBL,2) = IST+1, IST+2, ...`

---

## RustFoil Implementation (global_newton.rs)

### Source Code

```rust
fn compute_ue_from_mass_both(
    &self,
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    ue_inviscid: &[f64],
    dij: &DMatrix<f64>,
    surface: usize,
) -> Vec<f64> {
    let stations = if surface == 0 { upper_stations } else { lower_stations };
    let n = stations.len();
    let mut ue = vec![0.0; n];

    for i in 0..n {
        let panel_i = stations[i].panel_idx;
        let vti_i = if surface == 0 { 1.0 } else { -1.0 };

        let mut dui = 0.0;

        // Sum over upper surface (VTI = +1)
        // XFOIL starts from JBL=2 (first BL station, skipping stagnation)
        for j in 1..upper_stations.len() {
            let panel_j = upper_stations[j].panel_idx;
            let vti_j = 1.0;  // Upper surface

            if panel_i < dij.nrows() && panel_j < dij.ncols() {
                let ue_m = -vti_i * vti_j * dij[(panel_i, panel_j)];
                dui += ue_m * upper_stations[j].mass_defect;
            }
        }

        // Sum over lower surface (VTI = -1)
        // XFOIL starts from JBL=2 (first BL station, skipping stagnation)
        for j in 1..lower_stations.len() {
            let panel_j = lower_stations[j].panel_idx;
            let vti_j = -1.0;  // Lower surface

            if panel_i < dij.nrows() && panel_j < dij.ncols() {
                let ue_m = -vti_i * vti_j * dij[(panel_i, panel_j)];
                dui += ue_m * lower_stations[j].mass_defect;
            }
        }

        ue[i] = if i < ue_inviscid.len() {
            ue_inviscid[i] + dui
        } else {
            dui
        };
    }

    ue
}
```

### Formula

```
ue[i] = ue_inviscid[i] + Σ_upper(j=1..n) [ -vti_i * vti_j * dij[(panel_i, panel_j)] * mass[j] ]
                       + Σ_lower(j=1..n) [ -vti_i * vti_j * dij[(panel_i, panel_j)] * mass[j] ]
```

### VTI Sign Convention

| Surface | Index | vti_i Value | vti_j Value |
|---------|-------|-------------|-------------|
| Upper | surface=0 | +1.0 | +1.0 |
| Lower | surface=1 | -1.0 | -1.0 |

### Station Ranges

- **Outer loop:** `i = 0..n` — **INCLUDES stagnation (i=0)**
- **Inner loop:** `j = 1..n` — **skips stagnation (j=0)**
- **Surfaces:** Both upper and lower are summed (cross-surface coupling)

### Panel Index Mapping

From `initialize_surface_stations_with_panel_idx`:
- **Upper surface:** `panel_idx = stagnation_idx - (i - 1)` for stations after stagnation
- **Lower surface:** `panel_idx = stagnation_idx + i` for stations after stagnation

---

## Key Differences

### 1. Outer Loop Station Range ⚠️

| Aspect | XFOIL | RustFoil |
|--------|-------|----------|
| Outer loop start | IBL=2 (skips stagnation) | i=0 (includes stagnation) |
| Effect | Does NOT compute Ue for stagnation station | Computes Ue for ALL stations including stagnation |

**Impact:** RustFoil computes Ue for the stagnation station (i=0), while XFOIL leaves `UEDG(1,IS)` unchanged. This is a **potential issue** if:
- The stagnation Ue should remain at its initialized value
- Subsequent code assumes stagnation Ue is special

### 2. Inner Loop Station Range ✓

| Aspect | XFOIL | RustFoil |
|--------|-------|----------|
| Inner loop start | JBL=2 | j=1 |
| Effect | Skips stagnation in sum | Skips stagnation in sum |

**These are equivalent** — both skip the stagnation station (JBL=1 or j=0) when computing the mass defect contribution.

### 3. VTI Sign Convention ✓

| Surface | XFOIL | RustFoil |
|---------|-------|----------|
| Upper | VTI = +1 | vti = +1 |
| Lower | VTI = -1 | vti = -1 |

**Match!** Both implementations use the same sign convention.

### 4. DIJ Indexing ✓

| Aspect | XFOIL | RustFoil |
|--------|-------|----------|
| Row index | I = IPAN(IBL,IS) | panel_i = stations[i].panel_idx |
| Col index | J = IPAN(JBL,JS) | panel_j = stations[j].panel_idx |

**Match!** Both use panel indices from the station-to-panel mapping.

### 5. Cross-Surface Coupling ✓

| Aspect | XFOIL | RustFoil |
|--------|-------|----------|
| Surfaces summed | Both (JS=1,2) | Both (upper + lower) |

**Match!** Both implementations sum contributions from both surfaces.

---

## Recommended Fixes

### Fix 1: Skip Stagnation Station in Outer Loop

The outer loop should skip stagnation to match XFOIL:

**Current (incorrect):**
```rust
for i in 0..n {
    // ...
}
```

**Recommended (match XFOIL):**
```rust
for i in 1..n {  // Skip i=0 (stagnation)
    // ...
}
// Keep ue[0] = ue_inviscid[0] (stagnation unchanged)
```

However, this may require adjusting how the result array is indexed. A safer fix:

```rust
fn compute_ue_from_mass_both(...) -> Vec<f64> {
    let stations = if surface == 0 { upper_stations } else { lower_stations };
    let n = stations.len();
    let mut ue = vec![0.0; n];

    // Station 0 (stagnation): keep inviscid value unchanged
    if n > 0 && !ue_inviscid.is_empty() {
        ue[0] = ue_inviscid[0];
    }

    // Stations 1..n: compute from mass defect (XFOIL IBL=2 to NBL)
    for i in 1..n {
        let panel_i = stations[i].panel_idx;
        let vti_i = if surface == 0 { 1.0 } else { -1.0 };

        let mut dui = 0.0;

        // Sum over upper surface (j=1..n, skipping stagnation)
        for j in 1..upper_stations.len() {
            // ... (unchanged)
        }

        // Sum over lower surface (j=1..n, skipping stagnation)
        for j in 1..lower_stations.len() {
            // ... (unchanged)
        }

        ue[i] = if i < ue_inviscid.len() {
            ue_inviscid[i] + dui
        } else {
            dui
        };
    }

    ue
}
```

---

## Summary

| Aspect | Status | Notes |
|--------|--------|-------|
| VTI signs | ✓ Match | +1 upper, -1 lower |
| DIJ indexing | ✓ Match | Via panel_idx/IPAN |
| Inner loop (sum) | ✓ Match | Both skip stagnation |
| Cross-surface | ✓ Match | Both sum both surfaces |
| **Outer loop** | ⚠️ **Differ** | RustFoil includes stagnation, XFOIL skips |

The main difference is that RustFoil computes Ue for the stagnation station (i=0) while XFOIL does not (starts at IBL=2). This should be fixed to match XFOIL's behavior exactly.
