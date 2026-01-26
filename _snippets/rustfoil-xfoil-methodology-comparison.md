# RustFoil vs XFOIL: Methodology Comparison

## Executive Summary

RustFoil implements XFOIL's viscous-inviscid coupling algorithm with high fidelity. The core methodology matches XFOIL's approach: linear vorticity panel method, DIJ mass defect coupling, integral boundary layer equations, and global Newton iteration. Current accuracy: **CL error ~7.7%**, **CD within 1.45x XFOIL**.

---

## 1. Inviscid Panel Method

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Method** | Linear vorticity stream function | Linear vorticity stream function | ✅ |
| **Influence** | PSILIN (xpanel.f:99) | `psilin()` (influence.rs) | ✅ |
| **System** | (N+1)×(N+1) with Kutta row | (N+1)×(N+1) with Kutta row | ✅ |
| **Factorization** | LUDCMP/BAKSUB | nalgebra LU | ✅ |
| **Alpha superposition** | γ = cos(α)·γ₀ + sin(α)·γ₉₀ | Same | ✅ |

**Notes:**
- Both use Drela's linear vorticity formulation where γ represents surface velocity
- Panel midpoint evaluation for influence coefficients
- Kutta condition: γ_TE_upper + γ_TE_lower = 0

---

## 2. DIJ Matrix (Mass Defect Coupling)

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Location** | QDCALC (xpanel.f:1149) | `build_dij_matrix()` (dij.rs) | ✅ |
| **Formula** | DIJ = AIJ⁻¹ × BIJ | Same via back_substitute | ✅ |
| **BIJ source** | PSILIN with source panels | `psilin_with_sources()` | ✅ |
| **Application** | ΔUe = Σ DIJ × Δmass | Same in `compute_new_ue_via_dij()` | ✅ |

**Physical meaning:**
```
DIJ[i,j] = ∂Ue_i / ∂mass_j
mass = Ue × δ* (mass defect)
```

---

## 3. Boundary Layer Equations

### 3.1 BLVAR (Variable Computation)

| Variable | XFOIL (xblsys.f:784) | RustFoil (equations.rs:248) | Match |
|----------|---------------------|---------------------------|-------|
| H = δ*/θ | ✅ | ✅ | ✅ |
| Hk (kinematic) | `hkin()` | `hkin()` | ✅ |
| Hc (density) | `hct()` | `density_shape()` | ✅ |
| Hs (energy) | `hs_laminar/turbulent` | Same | ✅ |
| Us (slip) | 0.5·Hs·(1-(Hk-1)/(GBCON·H)) | Same | ✅ |
| Cf | `cfl()`, `cft()` | `cf_laminar()`, `cf_turbulent()` | ✅ |
| Dissipation | DI = f(Cf, Us, Hs) | Same formulas | ✅ |
| DE (BL thickness) | Green's correlation | Same | ✅ |

### 3.2 BLDIF (Residual Equations)

| Equation | XFOIL (xblsys.f:1552) | RustFoil (equations.rs:834) | Match |
|----------|----------------------|---------------------------|-------|
| **Amplification** | dn/dx = AX (Orr-Sommerfeld) | Same via `amplification_rate()` | ✅ |
| **Shear-lag** | d√Cτ/dx = (CQ - √Cτ)/SCC | Same | ✅ |
| **Momentum** | dθ/dx = Cf/2 - (H+2-M²)θ·d(lnUe)/dx | Same | ✅ |
| **Shape** | d(lnHs)/dx = (Cf/2-DI)/θ - (2Hc/Hs+1-H)·d(lnUe)/dx | Same | ✅ |

**Jacobian structure:** Both produce 3×5 blocks (VS1, VS2) for upstream/downstream stations.

---

## 4. Transition Prediction

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Method** | eⁿ envelope | Same | ✅ |
| **Ncrit default** | 9.0 | 9.0 | ✅ |
| **AX correlation** | `AXSET` | `amplification_rate()` | ✅ |
| **TRCHEK** | Implicit integration | `trchek2_stations()` | ✅ |
| **Ctau init** | 0.03 at transition | Same | ✅ |

---

## 5. Newton Iteration

### 5.1 System Construction

| Aspect | XFOIL (SETBL) | RustFoil | Match |
|--------|---------------|----------|-------|
| **Block structure** | 3×3 tridiagonal (VA, VB) | Same | ✅ |
| **VM matrix** | Mass influence on BL eqns | `build_with_vm_full()` | ✅ |
| **Forced changes** | DUE, DDS for nonlinearity | Same concept | ✅ |
| **Direct march first** | MRCHDU before SETBL | `march_surface()` first | ✅ |

**VM formula (both):**
```
VM[k,i,j] = VS2[k,3]·D2_M[j] + VS2[k,4]·U2_M[j] 
          + VS1[k,3]·D1_M[j] + VS1[k,4]·U1_M[j]
```

### 5.2 Solver

| Aspect | XFOIL (BLSOLV) | RustFoil | Match |
|--------|----------------|----------|-------|
| **Algorithm** | Forward elim + back sub | Same in `solve_blsolv_xfoil()` | ✅ |
| **VACCEL sparsity** | Yes (1e-3 threshold) | Yes | ✅ |
| **Row 3 special** | Mass equation pivot | Same | ✅ |

### 5.3 Update (Relaxation)

| Aspect | XFOIL (UPDATE) | RustFoil | Match |
|--------|----------------|----------|-------|
| **New Ue first** | UNEW = UINV + DIJ×(mass+Δmass) | `compute_new_ue_via_dij()` | ✅ |
| **Relaxation limits** | DHI=1.5, DLO=-0.5 | Same | ✅ |
| **Global RLX** | Single factor all vars | Same | ✅ |
| **Mass update** | mass = Ue × δ* (nonlinear) | Same | ✅ |
| **Ctau clamp** | [1e-7, 0.25] | [1e-7, 0.25] | ✅ |
| **Ue clamp** | None explicit | **[0.01, 5.0]** | ⚠️ Added |

---

## 6. Direct/Inverse Mode Switching (MRCHUE)

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Direct mode** | Fixed Ue, solve θ, δ* | Same | ✅ |
| **Inverse trigger** | Hk ≥ HLMAX (3.8 lam, 2.5 turb) | Same | ✅ |
| **Inverse target** | Htarg = f(dx/θ) | Same | ✅ |
| **Mode in iteration** | 25 Newton iters | Same | ✅ |

---

## 7. Stagnation Point Handling

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Initial location** | STFIND (γ sign change) | `find_stagnation()` (min |Ue|) | ⚠️ Different |
| **During iteration** | STMOVE relocates | **Not implemented** | ❌ Missing |
| **BL array shift** | Shifts upstream/downstream | N/A | ❌ Missing |

**Impact:** CL errors at low α (especially α=1°: 31.6% error) due to fixed stagnation point.

---

## 8. Force Computation

### 8.1 Lift (CL)

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Method** | Pressure integration ∮Cp·dx | **Circulation** 2∮Ue·ds | ⚠️ Different |
| **Compressibility** | Karman-Tsien | Not applied | ⚠️ |

### 8.2 Drag (CD)

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Friction** | ∫τ·dx where τ=Cf·Ue²/2 | ∫Cf·dx | ⚠️ Missing Ue² |
| **Pressure** | Squire-Young at wake end | Squire-Young at TE | ⚠️ No wake march |
| **Formula** | 2θ·Ue^((5+H)/2) | Same | ✅ |

---

## 9. Wake Modeling

| Aspect | XFOIL | RustFoil | Match |
|--------|-------|----------|-------|
| **Wake march** | Full BL march in wake | Simplified extrapolation | ❌ |
| **Wake stations** | From potential flow | Analytical growth model | ⚠️ |
| **TE combination** | θ_wake = θ_upper + θ_lower | Same | ✅ |
| **Wake H decay** | Marched to H→1 | Exponential model | ⚠️ |

---

## 10. Key Differences Summary

### Implemented Correctly ✅
1. Panel method (linear vorticity)
2. DIJ matrix computation
3. Boundary layer equations (BLVAR, BLDIF)
4. Transition (eⁿ method)
5. Newton system (VM matrix, BLSOLV)
6. Direct/inverse switching
7. Basic relaxation strategy

### Partial/Different ⚠️
1. **CL computation**: Circulation vs pressure integration
2. **Friction drag**: Missing Ue² weighting
3. **Wake model**: Simplified vs full march
4. **Stagnation initial**: Min |Ue| vs γ sign change
5. **Ue safeguards**: Added clamps [0.01, 5.0] not in XFOIL

### Missing ❌
1. **STMOVE**: Stagnation point relocation during iteration
2. **Wake march**: Full BL solution in wake
3. **Compressibility**: Karman-Tsien corrections

---

## 11. Current Accuracy

| Metric | Current | Target |
|--------|---------|--------|
| CL error (avg) | 7.7% | < 5% |
| CL error (max) | 31.6% (α=1°) | < 15% |
| CD ratio (avg) | 1.45x | < 1.2x |
| CD ratio (max) | 1.80x | < 1.5x |

---

## 12. Recommended Improvements (Priority Order)

### High Priority
1. **Implement STMOVE** - Will fix low-α CL errors (31.6% at α=1°)
2. **Add wake march** - Will improve CD accuracy
3. **Friction drag with Ue²** - Currently omitting velocity weighting

### Medium Priority
4. **Pressure-integration CL** - More accurate than circulation
5. **Karman-Tsien compressibility** - For M > 0.3

### Low Priority
6. **Stagnation by γ sign** - Minor accuracy improvement
7. **Newton forced changes refinement** - DUE/DDS terms

---

## Appendix: XFOIL Subroutine Mapping

| RustFoil | XFOIL | File | Status |
|----------|-------|------|--------|
| `psilin()` | PSILIN | xpanel.f:99 | ✅ |
| `build_dij_matrix()` | QDCALC | xpanel.f:1149 | ✅ |
| `blvar()` | BLVAR | xblsys.f:784 | ✅ |
| `bldif()` | BLDIF | xblsys.f:1552 | ✅ |
| `trchek2_stations()` | TRCHEK | xbl.f:1050 | ✅ |
| `march_surface()` | MRCHUE | xbl.f:549 | ✅ |
| `build_with_vm_full()` | SETBL | xbl.f:21 | ✅ |
| `solve_blsolv_xfoil()` | BLSOLV | xsolve.f:283 | ✅ |
| `update_xfoil_style()` | UPDATE | xbl.f:1316 | ✅ |
| N/A | STMOVE | xpanel.f:1628 | ❌ |
| N/A | Wake march | xbl.f | ❌ |
