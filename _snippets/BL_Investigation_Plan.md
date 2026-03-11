---
created: 2026-01-28
source: flexfoil-boundary-layer
tags: [type/investigation, status/active, topic/boundary-layer]
---

# Boundary Layer Investigation Plan

## Current Status

After the PSILIN/DQDM fix:
- ✅ CL at α=0° is correct (both symmetric and cambered airfoils)
- ✅ CL slope in linear region is reasonable (±5-8%)
- ❌ **No stall prediction** - CL continues to increase linearly beyond physical limits
- ❌ **CD overprediction** - 30-50% higher than XFOIL across the polar
- ❌ **Drag polar shape wrong** - Missing the characteristic stall "hook"

## Root Cause Hypothesis

The missing stall indicates RustFoil's BL solver is **not detecting or handling flow separation**. In XFOIL, stall occurs when:

1. Shape factor H increases (adverse pressure gradient)
2. When H exceeds threshold (~2.5-4.0), flow separates
3. XFOIL switches to **inverse mode** (prescribe H, solve for Ue)
4. Separation limits lift and increases pressure drag

---

## Phase 1 Results (Completed)

### Bug 1: Delayed Transition (CRITICAL)
**Symptom**: XFOIL transitions at x/c=0.0186, RustFoil at x/c=0.065
**Cause**: At the leading edge, Rθ is very low (~40-100), well below the critical value (~200-300) needed for amplification

**Evidence**:
- ibl=6 (x=0.018): Rθ = 42, no amplification
- ibl=24 (x=0.059): Rθ = 296, Hk = 5.43 (already separated!)
- Transition occurs at ibl=26 (x=0.065), too late

**Impact**: By the time transition occurs, the laminar BL has already separated, producing non-physical Hk values (7-8)

**Possible fixes**:
1. Check if XFOIL uses a different critical Rθ correlation
2. Investigate if there's a separation-induced transition mechanism
3. Lower Ncrit for high-alpha cases

### Bug 2: Inverse Mode htarg Formula (MAJOR)
**Symptom**: In inverse mode, htarg = 5.6 instead of expected 3.8
**Cause**: Formula `htarg = prev.hk + 0.03 * dx / theta` produces huge values when θ is tiny

**Evidence**:
- At ibl=24: prev.hk=3.05, dx=0.003, theta=3.5e-5
- htarg = 3.05 + 0.03 * 0.003 / 3.5e-5 = 3.05 + 2.57 = **5.62**
- Inverse mode converges to htarg=5.4 (wrong!)

**Impact**: Laminar separation produces Hk=7-8 instead of being limited to ~3.8

---

## Investigation Phases

### Phase 1: Instrument Shape Factor (H) Evolution
**Goal**: Determine if H is increasing correctly as α increases

**Tasks**:
1. Add debug output for H at each BL station during march
2. Compare H distribution at α=10° (pre-stall) vs α=15° (near stall)
3. Plot H vs x/c for both surfaces at various α
4. Compare against XFOIL H values if available

**Key Questions**:
- Is H increasing near the trailing edge at high α?
- Does H ever exceed 2.5 (turbulent separation threshold)?
- Is H similar to XFOIL in the attached flow region?

**Files to modify**:
- `crates/rustfoil-bl/src/debug.rs` - Add H tracking event
- `crates/rustfoil-coupling/src/march.rs` - Emit H values

### Phase 2: Check Separation Detection
**Goal**: Verify if/how RustFoil detects separation

**Tasks**:
1. Search for separation criteria in XFOIL (`SEPCHK`, `MRCHDU`)
2. Find RustFoil equivalent (if exists)
3. Compare threshold values (H_sep, Cf criteria)
4. Check if `is_separated` flag is ever set in stations

**XFOIL separation indicators**:
- `Cf < 0` (negative skin friction → separation)
- `H > H_SEP` (shape factor threshold)
- `TFORCE=.TRUE.` triggers transition/separation handling

**Files to examine**:
- `crates/rustfoil-bl/src/equations.rs` - BL equations
- `crates/rustfoil-coupling/src/march.rs` - Marching logic
- `Xfoil-instrumented/src/xbl.f` - XFOIL BL routines

### Phase 3: Investigate Inverse Mode
**Goal**: Understand XFOIL's inverse BL mode and check if RustFoil has it

**Background**:
In XFOIL, when separation is detected:
- Direct mode: prescribe Ue → solve for θ, δ*, H
- **Inverse mode**: prescribe H → solve for Ue

Inverse mode prevents the non-physical solution of continuing to accelerate the BL.

**Tasks**:
1. Identify XFOIL inverse mode trigger in `MRCHUE`/`MRCHDU`
2. Search RustFoil for inverse mode implementation
3. If missing, plan implementation

**Key XFOIL subroutines**:
- `MRCHUE` - Direct marching (prescribed Ue)
- `MRCHDU` - Inverse marching (prescribed H or δ*)
- `UICALC` - Computes inverse Ue

### Phase 4: Transition Model Verification
**Goal**: Ensure eN transition model matches XFOIL

**Tasks**:
1. Compare amplification factor (N) growth rate
2. Verify critical N value (typically N_crit = 9)
3. Check transition location (x_tr) against XFOIL
4. Verify transition triggers turbulent closures correctly

**Files**:
- `crates/rustfoil-bl/src/closures/` - Closure relations
- `crates/rustfoil-bl/src/equations.rs` - Amplification equation

### Phase 5: CD Component Analysis
**Goal**: Identify which drag component is overpredicted

**Tasks**:
1. Break down CD into components:
   - Cf (friction drag) - integrate τ_wall
   - Cp (pressure drag) - from momentum deficit
   - Squire-Young wake contribution
2. Compare each component to XFOIL
3. Identify which is overpredicted

**Current observation**: CD_pressure seems too high → likely separation-related

---

## Diagnostic Scripts to Create

### 1. `compare_shape_factor.py`
```python
# Compare H distribution along x/c
# Input: RustFoil JSON trace, XFOIL dump
# Output: H vs x/c plot for upper/lower surfaces
```

### 2. `compare_transition.py`
```python
# Compare transition location
# Check N-factor growth
# Verify x_tr matches
```

### 3. `separation_analysis.py`
```python
# Track where H > 2.5
# Track where Cf approaches 0
# Identify separation onset
```

---

## Priority Order

1. ~~**Phase 1** (Instrument H)~~ - ✅ COMPLETED
2. ~~**Phase 2** (Separation detection)~~ - ✅ COMPLETED (inverse mode exists but wasn't limiting correctly)
3. **Phase 3** (Inverse mode) - ⚠️ Fixed htarg computation, but TE separation still not occurring
4. **Phase 4** (Transition) - ⚠️ Added separation-induced transition, but was too aggressive
5. **Phase 5** (CD components) - Quantify remaining errors

### Current Status (2026-01-28)

**Fixes implemented:**
- htarg gradient clamping (prevents runaway to Hk=7+)
- Separation-induced transition (threshold 4.3)

**Remaining issue:**
The turbulent BL at high alpha doesn't show the H and Cf evolution that leads to TE separation:
- XFOIL: H=3.03, Cf=0.000034 at TE (stall)
- RustFoil: H=2.50, Cf=0.000227 at TE (attached)

**Next investigation:** Why doesn't the turbulent BL develop high H near the TE?

---

## Phase 6: Re-March Fix Implementation (2026-01-28)

### Status: IMPLEMENTED BUT EXPOSED DEEPER BUGS ⚠️

**What was fixed:**
- Modified `viscal.rs` to re-march BL at each Newton iteration (lines 827-910)
- Stations now appear multiple times in debug trace (re-march confirmed working)
- Code compiles and runs successfully

**What was discovered:**
The re-march fix revealed that the single-march approach was **converging to a wrong answer**:

| Issue | Evidence |
|-------|----------|
| **TE H calculation** | Off by 170% at α=10° (RustFoil: 1.83, XFOIL: 4.91) |
| **Lower surface broken** | H wrong by 40-55% at all angles tested |
| **VI coupling unstable** | RMS explodes to 120 at α=15° iteration 7 |
| **Compounding errors** | Re-marching amplifies underlying bugs, causing divergence |

**Conclusion:**
✅ Root cause diagnosis was CORRECT - must re-march for stall prediction  
❌ Implementation exposed that basic BL accuracy is broken  
⚠️ Cannot proceed to stall testing until fundamental accuracy is fixed

**Required before continuing:**
1. Debug trailing edge BL integration (why H off by 170%)
2. Fix lower surface calculations
3. Test at simple cases (α=0°, 4°) to isolate bugs
4. Add detailed station-by-station XFOIL comparison
5. Fix Newton solver stability at high alpha

See `REMARCH_FIX_TEST_RESULTS.md` for complete test data.

---

## Expected Outcomes

After this investigation, we should understand:
1. Whether H is evolving correctly
2. If/where separation is detected
3. Whether inverse mode is implemented
4. What specifically causes the CD overprediction

The most likely fix is **implementing inverse BL mode**, which allows the solver to handle separated flow correctly.

---

## References

- [[PSILIN_Port_Implementation_Plan]] - Previous DIJ fix
- XFOIL Theory Doc: Section 5 (Viscous Formulation)
- Drela, M. "XFOIL: An Analysis and Design System for Low Reynolds Number Airfoils"
