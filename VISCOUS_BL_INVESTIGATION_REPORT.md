# Viscous Boundary Layer Investigation Report

**Date**: 2026-03-08  
**Branch**: `feature/boundary-layer`  
**Scope**: Full comparison of Rustfoil BL vs XFOIL reference implementation

---

## Executive Summary

The Rustfoil inviscid solver appears correct. The BL closure relations (Hkin, Hs, Cf, CD, Hc, DAMPL) are **formula-identical** to XFOIL. The problems are concentrated in three areas:

1. **The wake is solved with simplified empirical physics instead of the full BL equations** (Critical)
2. **The TE junction coupling (VZ block, D1_M redefinition) is fundamentally wrong** (Critical)
3. **The global Newton solver is missing XFOIL's second RHS column for α/CL sensitivity** (Critical)
4. Several medium-severity Jacobian completeness issues and update relaxation differences

These issues explain the observed symptoms: CL 5-16% wrong, CD pressure 200%+ wrong, Newton instability at high α, and incorrect H at the trailing edge.

---

## Part 1: What's Working

### 1.1 Inviscid Solver (VERIFIED CORRECT)

The `rustfoil-inviscid` crate implements the XFOIL panel method faithfully:

- Linear vorticity panels with stream-function formulation
- Kutta condition: γ₀ + γₙ₋₁ = 0
- PSILIN influence coefficients with correct SGN branch-cut handling
- TE panel treatment with PSIG/PGAM and HOPI
- CLCALC wind-axis force integration
- DIJ matrix construction (airfoil-airfoil, airfoil-wake, wake-wake)

**Note**: There is a separate inviscid implementation in `rustfoil-solver/src/inviscid/mod.rs` that may have TE sign differences. The production path uses `rustfoil-inviscid` via `setup_from_body()`, which is correct.

### 1.2 BL Closure Relations (VERIFIED IDENTICAL)

All 9 closure functions match XFOIL exactly:

| Closure | Rust File | Verdict |
|---------|-----------|---------|
| Hkin (Whitfield) | `closures/hkin.rs` | Identical |
| Hs laminar | `closures/hs.rs:40-81` | Identical |
| Hs turbulent | `closures/hs.rs:107-198` | Identical |
| Cf laminar (Falkner-Skan) | `closures/cf.rs:41-65` | Identical |
| Cf turbulent (Coles) | `closures/cf.rs:85-135` | Identical |
| DI laminar | `closures/dissipation.rs:108-127` | Identical |
| DI turbulent | `closures/dissipation.rs:159-176` | Identical |
| DI wake (DILW) | `closures/dissipation.rs:208-230` | Identical |
| Hc (density shape) | `closures/density.rs:30-36` | Identical |
| DAMPL (amplification) | `closures/transition.rs:77-165` | Identical |
| AXSET (averaged amplification) | `closures/transition.rs:261-414` | Identical |

### 1.3 BLPAR Constants (VERIFIED IDENTICAL)

All constants in `constants.rs` match `BLPAR.INC` exactly (SCCON=5.6, GACON=6.70, GBCON=0.75, GCCON=18.0, DLCON=0.9, CTRCON=1.8, CTRCEX=3.3, DUXCON=1.0, CTCON=0.5/(GACON²·GBCON)).

### 1.4 BLVAR Secondary Variables (MOSTLY CORRECT)

The `blvar()` function in `equations.rs:248-806` correctly computes: Hk clamping, Us, CQ (equilibrium shear coefficient), Cf selection (max of laminar/turbulent), DI with DFAC correction, DE boundary layer thickness. First BL station matches XFOIL to machine precision in trace comparisons.

### 1.5 BLDIF Residuals (CORRECT)

The three residual equations match XFOIL:
- Row 1 laminar: `AMPL2 - AMPL1 - AX·(X2-X1) = 0`
- Row 1 turbulent: shear-lag `SCC·(CQA-SA·ALD)·DXI - DEA·2·SLOG + ...`
- Row 2 momentum: `TLOG + BTMP·ULOG - XLOG·0.5·CFX = 0`
- Row 3 shape: `HLOG + BTMP·ULOG + XLOG·(0.5·CFX - DIX) = 0`

### 1.6 BL Marching (CORRECT FOR SINGLE SURFACE)

The `march_surface()` function correctly marches the BL from stagnation, handling transition via TRCHEK2. Transition locations are within ~6% of XFOIL for the cases checked.

---

## Part 2: Critical Bugs

### BUG 1 (CRITICAL): Wake Uses Simplified Empirical Physics

**Location**: `crates/rustfoil-coupling/src/wake.rs:273-349`

**Problem**: The wake is solved with a simplified momentum conservation + empirical dissipation model:

```rust
let theta_conserved = prev.theta * ue_ratio * ue_ratio;
let dissipation_rate = 0.15;
let dissipation_factor = (-dissipation_rate * dx).exp();
station.theta = theta_conserved * dissipation_factor;
```

**XFOIL** marches the wake using the **full BL equations** (BLSYS + MRCHUE with `WAKE=.TRUE.`), solving the momentum integral, shape parameter, and shear-lag equations with proper Newton iteration at each station. The wake edge velocity comes from the inviscid panel solution (not an empirical recovery formula).

**Impact**: 
- Wake θ, δ*, and H evolution are all wrong
- CD computed via Squire-Young formula is directly affected (observed: pressure CD 200%+ too high)
- The empirical `dissipation_rate = 0.15` and `wake_edge_velocity()` function have no theoretical basis

**Fix**: Replace `solve_wake_station()` with actual BLDIF(FlowType::Wake) Newton iteration at each wake station, using the inviscid wake panel velocities for Ue.

---

### BUG 2 (CRITICAL): VZ Block at TE Uses Wrong Derivatives

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:1007-1010`

**Problem**: 

```rust
for k in 0..3 {
    self.vz[k][0] = self.vs1_delta[jv_lower_wake][k];
    self.vz[k][1] = self.vs1_ue[jv_lower_wake][k];
}
```

This stores δ* and Ue derivatives from the lower wake station.

**XFOIL** (`xbl.f:455-468`) computes VZ using **chain-rule derivatives** through the TE combination:

```fortran
VZ(1,1) = VS1(1,1)*CTE_CTE1
VZ(1,2) = VS1(1,1)*CTE_TTE1 + VS1(1,2)*TTE_TTE1
VB(1,1,IV) = VS1(1,1)*CTE_CTE2
VB(1,2,IV) = VS1(1,1)*CTE_TTE2 + VS1(1,2)*TTE_TTE2
```

Where:
- `CTE_CTE1 = θ_upper / TTE` (∂ctau_wake/∂ctau_upper)
- `CTE_TTE1 = (ctau_upper - CTE)/TTE` (∂ctau_wake/∂theta_upper)
- `TTE_TTE1 = 1.0` (∂theta_wake/∂theta_upper)

The VB block is also **redefined** at the first wake station to use derivatives w.r.t. the lower TE station.

**Impact**: The entire TE junction coupling is broken. The Newton system cannot correctly propagate changes from the upper surface through the TE combination into the wake.

**Fix**: Implement the full TESYS chain-rule transformation for VZ and redefine VB at the first wake station.

---

### BUG 3 (CRITICAL): D1_M Not Redefined at First Wake Station

**Location**: `crates/rustfoil-coupling/src/global_newton.rs` (build_vm_global)

**Problem**: In XFOIL (`xbl.f:266-274`), the `D1_M(JV)` vector is redefined at the first wake station because the "upstream" δ* is the combined TE δ* which depends on **both** surfaces' mass defects:

```fortran
DO 35 JS=1, 2
  DO 350 JBL=2, NBL(JS)
    JV = ISYS(JBL,JS)
    D1_M(JV) = DTE_UTE1*UTE1_M(JV) + DTE_UTE2*UTE2_M(JV)
  350 CONTINUE
35 CONTINUE
D1_M(JVTE1) = D1_M(JVTE1) + DTE_MTE1
D1_M(JVTE2) = D1_M(JVTE2) + DTE_MTE2
```

Rustfoil does not perform this redefinition. The VM matrix at the first wake station uses standard single-surface D1_M instead of the combined-TE D1_M.

**Impact**: VM matrix entries for the first wake station are wrong, missing cross-surface coupling through the TE.

**Fix**: Add the D1_M redefinition at the first wake station in `build_vm_global()`.

---

### BUG 4 (CRITICAL): Missing Second RHS Column in BLSOLV

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:1280` (solve_global_system)

**Problem**: XFOIL's BLSOLV (`xsolve.f:332-493`) processes **two RHS columns** simultaneously: `VDEL(k,1,IV)` (residuals) and `VDEL(k,2,IV)` (α/CL sensitivity). Both columns go through the entire forward sweep and back-substitution.

Rustfoil only has a single RHS column `vdel[iv][k]`. The second column that carries the α-sensitivity is completely absent.

In XFOIL, the UPDATE routine uses `VDEL(3,2,JV)` to compute `DUI_AC` and `U_AC`, which are then used to compute the `DAC` update for the global α/CL coupling variable.

**Impact**: The solver cannot compute the α/CL update that XFOIL uses for outer viscous-inviscid coupling convergence. This may cause the iteration to oscillate or fail to converge.

**Fix**: Add the second RHS column to the BLSOLV implementation. Populate it with the α-sensitivity (dFk/dα) terms as XFOIL does.

---

## Part 3: High-Severity Bugs

### BUG 5 (HIGH): Forced Changes at First Wake Station Wrong

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:1088-1167`

**Problem**: XFOIL (`xbl.f:277-279`) computes forced changes at the first wake station specially:

```fortran
DUE1 = 0.
DDS1 = DTE_UTE1*(UEDG(IBLTE(1),1) - USAV(IBLTE(1),1))
     + DTE_UTE2*(UEDG(IBLTE(2),2) - USAV(IBLTE(2),2))
```

The forced DUE1 is **zero** and DDS1 is a combination of **both** surfaces' Ue mismatches. Rustfoil uses the standard single-surface formula for all stations including the wake.

**Impact**: Forced changes at TE wake are wrong, affecting convergence.

---

### BUG 6 (HIGH): VTI Index Cross-Surface at Lower Surface Start

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:883`

**Problem**: When computing `U1_M` for the first lower surface station (`iv = n_upper`), the code uses `self.vti[iv - 1]` which indexes the last upper surface station's VTI instead of the lower surface's stagnation point VTI.

XFOIL uses `VTI(IBL-1, IS)` which stays on the same surface.

**Impact**: VM matrix sign error near the lower surface start, corrupting the Newton system near stagnation.

---

### BUG 7 (HIGH): Wake Edge Velocity From Empirical Formula

**Location**: `crates/rustfoil-coupling/src/wake.rs:138-172`

**Problem**: Wake Ue is computed from an empirical recovery formula instead of from the inviscid panel solution. XFOIL uses the actual panel-method-computed velocities at wake panels.

```rust
// Empirical formula - NOT from panel method
let recovery_factor = 1.0 - (-recovery_rate * x_far).exp();
ue_te + (1.0 - ue_te) * recovery_factor
```

**Impact**: Wrong wake velocity profile directly affects wake BL evolution and thus CD computation.

**Fix**: Use the inviscid panel solution's wake velocities (available from the DIJ setup which already computes wake influence coefficients).

---

## Part 4: Medium-Severity Issues

### BUG 8 (MEDIUM): HWA Wake Displacement Term Missing from BTMP

**Location**: `crates/rustfoil-bl/src/equations.rs:1433,1585`

**Problem**: In XFOIL's BLDIF, the momentum and shape equation BTMP coefficients include `HWA = 0.5*(DW1/T1 + DW2/T2)`:

```fortran
BTMP_momentum = HA + 2.0 - MA + HWA
BTMP_shape = 2.0*HCA/HSA + 1.0 - HA - HWA
```

Rustfoil omits HWA:

```rust
let btmp = ha + 2.0 - ma_avg;          // Missing + HWA
let btmp_shape = 2.0 * hca / hsa + 1.0 - ha;  // Missing - HWA
```

**Impact**: Affects wake and TE region accuracy. For thin airfoils at moderate angles, effect is small but non-negligible.

---

### BUG 9 (MEDIUM): Shear-Lag Jacobian Missing UPW_U and DE_U Terms

**Location**: `crates/rustfoil-bl/src/equations.rs:1332,1392`

**Problem**: In XFOIL's BLDIF, the turbulent shear-lag Jacobian row (row 1) includes:

```fortran
VS1(1,4) = Z_U1 + Z_UPW*UPW_U1 + Z_DE1*DE1_U1 + Z_US1*US1_U1
```

And the UQ derivative includes:

```fortran
UQ_U1 = (1-UPW)*(UQ_CFA*CF1_U1 + UQ_HKA*HK1_U1) + UQ_UPW*UPW_U1
```

Rustfoil is missing the `Z_UPW*UPW_U1`, `Z_DE1*DE1_U1`, and `UQ_HKA*HK1_U1` terms. While `UPW_U` and `DE_U` are small in most conditions, they contribute to Newton convergence robustness.

**Impact**: Slower convergence, particularly near transition where UPW changes rapidly.

---

### BUG 10 (MEDIUM): VACCEL Uses Single Fixed Threshold

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:1288`

**Problem**: XFOIL scales VACCEL for rows 2 and 3 by `2.0/(S(N)-S(1))`:

```fortran
VACC1 = VACCEL
VACC2 = VACCEL * 2.0 / (S(N) - S(1))
VACC3 = VACCEL * 2.0 / (S(N) - S(1))
```

Rustfoil uses a single `vaccel = 0.005` for all rows.

**Impact**: Rows 2 and 3 may drop too many VM entries, reducing accuracy of the mass-influence coupling.

---

### BUG 11 (MEDIUM): Non-Equivalent Update Relaxation Logic

**Location**: `crates/rustfoil-coupling/src/global_newton.rs` (apply_global_updates)

**Problem**: XFOIL uses signed limit checks:

```fortran
IF(RDN1 .GT. DHI) RLX = DHI/DN1
IF(RDN1 .LT. DLO) RLX = DLO/DN1
```

Rustfoil uses `abs()` checks that don't preserve the signed limit behavior. Additionally, extra non-XFOIL safeguards (theta clamp 0.5×-2×, ctau min 1e-6, Ue min 0.01, H clamp 1.0-20.0) interfere with Newton convergence.

**Impact**: Relaxation can be too aggressive or too conservative, affecting convergence speed and stability.

---

### BUG 12 (MEDIUM): SIMI Station Residuals Zeroed

**Location**: `crates/rustfoil-coupling/src/global_newton.rs:636-639`

**Problem**: Residuals for the first station (SIMI) are explicitly zeroed:

```rust
if i == 1 {
    residuals.res_mom = 0.0;
    residuals.res_shape = 0.0;
}
```

XFOIL does NOT zero the SIMI residuals. This removes the forced-change driving term.

**Impact**: Stagnation region may not converge properly.

---

### BUG 13 (MEDIUM): Hardcoded Ncrit in AXSET

**Location**: `crates/rustfoil-bl/src/equations.rs:1098`

**Problem**: `axset_full` is called with `ncrit = 9.0` hardcoded instead of using the user-supplied value.

**Impact**: Transition prediction always uses Ncrit=9 regardless of configuration.

---

## Part 5: Compressibility Gaps (LOW priority for incompressible)

These are correctly noted as simplified in the code and only matter for M > 0.3:

- **M_U = 0** hardcoded (`equations.rs:293`): velocity-Mach coupling derivative
- **Z_MA·M_U Jacobian terms** missing in momentum equation (`equations.rs:1547`)
- **BLPRV compressible velocity correction** not implemented (U = UEI × compressibility factor)

---

## Part 6: Trace Data Comparison (NACA 0012, Re=3e6, α=4°)

### Final Results

| Quantity | XFOIL | Rustfoil | Difference |
|----------|-------|----------|------------|
| CL | 0.4424 | 0.3737 | -15.5% |
| CD total | 0.006183 | 0.007162 | +15.8% |
| CD friction | 0.004579 | 0.004595 | +0.3% |
| CD pressure | 0.000845 | 0.002567 | **+204%** |
| Transition upper (x/c) | 0.1475 | 0.1567 | +6.2% |
| Transition lower (x/c) | 0.8704 | 0.8383 | -3.7% |

### Key Observations

1. **Friction drag matches within 0.3%** -- confirms closure relations are correct
2. **Pressure drag is 3× too high** -- confirms wake/TE coupling is broken
3. **CL 15.5% too low** -- Newton iteration not fully converging the VI coupling
4. **Transition locations within 4-6%** -- BL march and e^N model are working reasonably

### Station-Level Match

| Station | Quantity | XFOIL | Rustfoil | Match? |
|---------|----------|-------|----------|--------|
| ibl=2 (x~0.001) | theta | 2.1332e-5 | 2.1332e-5 | Exact |
| ibl=2 (x~0.001) | Hk | 2.200 | 2.200 | Exact |
| ibl=2 (x~0.001) | Cf | 0.19224 | 0.19224 | Exact |
| Upper TE | theta | 0.003790 | 0.003838 | +1.3% |
| Upper TE | H | 1.675 | 1.616 | -3.5% |
| Lower TE | theta | 0.000785 | 0.000922 | **+17.5%** |

### NACA 0012 Symmetry Check

At α=0° for a symmetric airfoil, Rustfoil reports CL = -0.047 instead of 0.000. This indicates a numerical asymmetry in the VI coupling, likely from paneling or stagnation point handling.

---

## Part 7: Previously Identified and Fixed Issues

These have been addressed and should not be re-investigated:

| Issue | Status |
|-------|--------|
| CD reported as total instead of pressure | Fixed (forces.rs) |
| Ue min used `.min(0.01)` instead of `.max(0.01)` | Fixed |
| U1_M used `.max(self.vti[iv])` instead of `self.vti[iv-1]` | Fixed |
| Newton relaxation missing Ue and δ* checks | Fixed |
| Re-march at each iteration | Implemented |

---

## Part 8: Recommended Fix Priority

### Phase 1: Critical (will fix convergence and CD)

1. **Replace simplified wake with full BLDIF wake marching** (Bug 1)
   - Use inviscid wake panel velocities for Ue (Bug 7)
   - Implement proper Newton solve at each wake station
   
2. **Fix VZ block with TESYS chain-rule derivatives** (Bug 2)
   - Implement CTE_CTE1, CTE_TTE1, TTE_TTE1 etc.
   - Redefine VB at first wake station

3. **Add D1_M redefinition at first wake station** (Bug 3)
   - Compute DTE_UTE1, DTE_UTE2 derivatives
   - Loop over both surfaces to build combined D1_M

4. **Add second RHS column to BLSOLV** (Bug 4)
   - Carry α-sensitivity through forward sweep and back-substitution
   - Implement DAC computation in UPDATE

### Phase 2: High (will improve accuracy and robustness)

5. Fix forced changes at first wake station (Bug 5)
6. Fix VTI cross-surface index (Bug 6)
7. Add HWA term to BTMP (Bug 8)

### Phase 3: Medium (polish)

8. Complete shear-lag Jacobian (Bug 9)
9. Scale VACCEL per row (Bug 10)
10. Fix update relaxation to match XFOIL signed logic (Bug 11)
11. Remove SIMI residual zeroing (Bug 12)
12. Use config Ncrit in AXSET (Bug 13)

---

## Appendix: Architecture Reference

### Crate Structure

```
rustfoil-bl/          BL equations, closures, state
rustfoil-inviscid/    Panel method (XFOIL-faithful)
rustfoil-core/        Geometry, splines, NACA generation
rustfoil-coupling/    VI coupling, Newton, marching, wake
rustfoil-solver/      Top-level solver, viscal, forces
rustfoil-cli/         Command-line interface
```

### XFOIL ↔ Rust File Mapping

| XFOIL | Rust |
|-------|------|
| xblsys.f BLVAR | rustfoil-bl/equations.rs `blvar()` |
| xblsys.f BLDIF | rustfoil-bl/equations.rs `bldif()` |
| xblsys.f TRDIF | rustfoil-bl/equations.rs `trdif()` |
| xblsys.f TRCHEK2 | rustfoil-bl/closures/transition.rs `trchek2()` |
| xbl.f MRCHUE | rustfoil-coupling/march.rs `march_surface()` |
| xbl.f SETBL | rustfoil-coupling/global_newton.rs `build_system()` |
| xbl.f TESYS | rustfoil-coupling/wake.rs `combine_te_for_wake()` |
| xbl.f UPDATE | rustfoil-coupling/global_newton.rs `apply_global_updates()` |
| xsolve.f BLSOLV | rustfoil-coupling/solve.rs + global_newton.rs |
| xoper.f VISCAL | rustfoil-solver/viscal.rs `solve_viscous_two_surfaces()` |
| BLPAR.INC | rustfoil-bl/constants.rs |
