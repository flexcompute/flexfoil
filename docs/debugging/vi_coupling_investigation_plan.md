# VI Coupling Investigation Plan

## Objective
Identify and fix differences between RustFoil and XFOIL's viscous-inviscid (VI) coupling mechanism that cause the observed CL errors:
- Low α (0-4°): RustFoil CL ~10% higher than XFOIL
- High α (6-8°): RustFoil CL ~5% lower than XFOIL
- α=0: RustFoil gives CL=-0.047 (should be 0 for symmetric airfoil)

## Background

### VI Coupling Flow
```
┌─────────────────────────────────────────────────────────────────┐
│                    XFOIL VI Coupling Flow                        │
├─────────────────────────────────────────────────────────────────┤
│ 1. MRCHDU: March BL with current Ue, compute θ, δ*, mass        │
│ 2. UESET:  Ue_new = Ue_inv + DIJ * mass_defect                  │
│ 3. SETBL:  Build global Newton system with Ue-mass coupling     │
│ 4. BLSOLV: Solve block-tridiagonal system for deltas            │
│ 5. UPDATE: Apply relaxed updates to BL vars and Ue              │
│ 6. Repeat until converged                                        │
└─────────────────────────────────────────────────────────────────┘
```

### Key Equations
- Mass defect: `mass = Ue * δ*`
- Ue update: `Ue[i] = Ue_inv[i] + Σⱼ(-VTI[i]*VTI[j]*DIJ[i,j] * mass[j])`
- Sensitivity: `dUe/d(mass) = -VTI[i]*VTI[j]*DIJ[i,j]`

## Test Matrix

| Foil | Reynolds | Alpha Range | Purpose |
|------|----------|-------------|---------|
| NACA 0012 | 3e6 | -4, 0, 4, 8 | Symmetric baseline |
| NACA 0012 | 1e6, 6e6 | 0, 4, 8 | Re sensitivity |
| NACA 2412 | 3e6 | -4, 0, 4, 8 | Cambered validation |

## Investigation Steps

### Phase 1: DIJ Matrix Comparison
Compare the DIJ influence coefficient matrix between XFOIL and RustFoil.

**Key questions:**
- Are DIJ values identical for the same paneling?
- How is DIJ computed in each solver?

**Traces needed:**
- XFOIL: Enable `DBGFULLDIJ` in QDCALC
- RustFoil: Add DIJ dump to inviscid solver

### Phase 2: UESET Comparison (Ue from mass defect)
Compare how each solver computes Ue from mass defect.

**Key variables to compare:**
- `UINV[i]` - inviscid edge velocity
- `MASS[i]` - mass defect (Ue * δ*)
- `UEDG[i]` - updated edge velocity after VI coupling

**Traces needed:**
- XFOIL: `DBGUESET` dumps these
- RustFoil: Add `compute_ue_from_mass_both()` debug output

### Phase 3: Newton System Comparison
Compare the global Newton system matrices at iteration 1.

**Key matrices:**
- `VA[4,4]` - diagonal block (current station)
- `VB[4,4]` - off-diagonal block (upstream station)
- `VM[4,nsys]` - VI coupling (mass defect influence)
- `VDEL[4,2]` - solution deltas

**Traces needed:**
- XFOIL: `DBGSETBLSYSTEM` dumps full system
- RustFoil: `emit_setbl_debug()` already dumps samples

### Phase 4: Update Comparison
Compare how Newton deltas are applied.

**Key variables:**
- `DUI` - Ue change from mass update
- Relaxation factor
- Final `UEDG` after update

**Traces needed:**
- XFOIL: `DBGUPDATEDETAIL` dumps before/after
- RustFoil: Existing debug prints in `apply_global_updates()`

### Phase 5: Cross-Surface Coupling
Check if upper/lower surface coupling differs.

**Key questions:**
- Does RustFoil's separate upper/lower handling match XFOIL's unified arrays?
- Are VTI signs consistent?
- Is trailing edge coupling correct?

## Comparison Scripts

### 1. `compare_dij.py`
Compare DIJ matrices from both solvers.

### 2. `compare_ueset.py`
Compare UINV, MASS, UEDG at each VI iteration.

### 3. `compare_newton_system.py`
Compare VA, VB, VM, VDEL at iteration 1.

### 4. `compare_ue_updates.py`
Compare Ue updates at each global Newton iteration.

## Findings (2026-01-26)

### Phase 1: DIJ Matrix Comparison - COMPLETED
**Status**: DIJ matrices match within 2-5%, which is acceptable.
- RustFoil DIJ[0,0] = -49.51 vs XFOIL = -51.41 (+3.7%)
- RustFoil DIJ[0,1] = +39.02 vs XFOIL = +40.56 (-3.8%)

### Root Cause Identified: Wake Panel Coverage in DIJ

**CRITICAL FINDING**: XFOIL and RustFoil have different DIJ matrix sizes:
- **XFOIL DIJ**: 183×183 (160 airfoil panels + 23 wake panels)
- **RustFoil DIJ**: 160×160 (airfoil panels only)

This causes two problems:

1. **Wake stations cannot contribute to UESET in RustFoil**: In XFOIL, wake BL stations 
   (x > chord length) have panel indices N+1 to N+NW and contribute to the mass defect 
   coupling via DIJ. RustFoil's DIJ doesn't cover wake panels, so these contributions 
   are missing.

2. **Panel index clamping for upper surface**: RustFoil assigns panel_idx using 
   `saturating_sub(stagnation_idx, station_idx)`, which clamps negative results to 0. 
   For stations beyond the panel count, this creates duplicate panel indices at 0, 
   causing incorrect DIJ contributions.

### DUI Analysis at Station 1 (First Iteration)

| Solver | Upper Contrib | Lower Contrib | Total DUI |
|--------|--------------|---------------|-----------|
| XFOIL | ~-0.17 | ~+0.04 | ~-0.17 |
| RustFoil (original) | -0.49 | +0.045 | -0.45 |
| RustFoil (with wake skip) | +0.05 | -0.03 | +0.02 |

The original RustFoil DUI was 2.6× too large (negative) due to clamped panel indices.
Skipping wake stations overcorrects and makes DUI small positive.

### Required Fix (Architectural)

To properly match XFOIL's VI coupling, RustFoil needs:

1. **Extend DIJ matrix to include wake panels**: The inviscid solver should compute 
   DIJ for both airfoil and wake panels (size = N + N_wake).

2. **Set up wake panel indices correctly**: Wake BL stations should have panel indices
   from N+1 to N+N_wake, not clamped to airfoil panel range.

3. **Include wake stations in UESET sum**: Once DIJ covers wake panels, the UESET
   loop should include all BL stations (not just airfoil stations).

### Workaround (Current State)

Until the architectural fix is implemented, the code retains the original behavior
which has excessive DUI magnitude. This causes:
- CL ~10% higher than XFOIL at low alpha
- CL ~5% lower than XFOIL at high alpha
- Non-zero CL at alpha=0 for symmetric airfoils

## Success Criteria

1. DIJ matrices match to <1e-6
2. UINV matches (should already - from inviscid match)
3. MASS defect matches at each station
4. Ue updates follow same pattern
5. CL at α=0 is exactly 0 for symmetric airfoil
6. CL matches XFOIL within 1% across test matrix

## Files to Modify

### XFOIL Instrumentation
- `Xfoil-instrumented/src/xpanel.f` - DIJ dump in QDCALC
- `Xfoil-instrumented/src/xbl.f` - Enhanced UESET/SETBL dumps

### RustFoil Debug Output
- `crates/rustfoil-coupling/src/global_newton.rs` - Add VI coupling traces
- `crates/rustfoil-inviscid/src/solver.rs` - Add DIJ dump

### Comparison Scripts
- `scripts/compare_dij.py`
- `scripts/compare_ueset.py`
- `scripts/compare_newton_system.py`
- `scripts/compare_ue_updates.py`
- `scripts/generate_vi_traces.py` - Generate all traces for test matrix
