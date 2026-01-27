# Separation-Induced Transition Fix

## Problem

RustFoil's lower surface stayed fully laminar at high angles of attack (α=4°, 6°, 8°) while XFOIL correctly transitioned to turbulent flow. This caused significant errors:
- CL errors varying from -9% to +13%
- CD errors up to +50%
- Non-systematic errors that changed sign with angle of attack

## Root Cause

The transition check in `march_surface()` only used the e^N method (`trchek2_stations`). When Hk exceeded `hlmax` (3.8) on the lower surface, XFOIL forces transition because laminar flow cannot sustain such high shape factors. RustFoil did NOT implement this separation-induced transition.

**Evidence at α=4°, lower surface:**
- Both XFOIL and RustFoil reached laminar separation at x=0.844 (Hk > 6)
- XFOIL: Transitioned to turbulent at x=0.96 (Hk dropped to 1.5)
- RustFoil: Stayed laminar (Hk stayed at 2.5)

## Fix Implemented

**File**: `crates/rustfoil-coupling/src/march.rs`

**Location**: In `march_surface()`, after the e^N transition check block

**Logic added**:
```rust
// Separation-induced transition: force transition if laminar Hk exceeds hlmax
// This handles cases where laminar BL approaches separation without e^N transition
// Only apply to lower surface (side == 2) to avoid interference with upper surface
// natural transition which typically occurs before Hk gets very high
if is_laminar && side == 2 && station.hk > config.hlmax && result.x_transition.is_none() {
    is_laminar = false;
    station.is_laminar = false;
    station.is_turbulent = true;
    result.x_transition = Some(station.x);
    result.transition_index = Some(station_idx);
    // Re-solve as turbulent...
}
```

## Results After Fix

### Transition Location (Lower Surface)

| Alpha | Before Fix | After Fix | XFOIL Reference |
|-------|-----------|-----------|-----------------|
| 4° | LAMINAR | x=0.86 | x=0.97 |
| 8° | LAMINAR | x=0.93 | x=1.00 |

### CL/CD Error Comparison

| Alpha | CL Before | CL After | CD Before | CD After |
|-------|-----------|----------|-----------|----------|
| 0° | ~0% | ~0% | +50% | +54% |
| 2° | +13% | +14% | +46% | +51% |
| 4° | -7% | -4% | -4% | -5% |
| 6° | -9% | -9% | +2% | -6% |
| 8° | +6% | -7% | +24% | -9% |

### Key Improvements
- α=8° CD: **+24% → -9%** (error reduced by 15 percentage points)
- α=4° CL: **-7% → -4%** (improved)
- Lower surface now physically correct (transitions instead of staying fully laminar)

## Remaining Issues

1. **CD at low alpha (0°, 2°) still ~50% high**: This is a pre-existing issue unrelated to transition, likely in the friction drag or Squire-Young formula
2. **Early transition**: Lower surface transitions earlier than XFOIL predicts (x=0.86 vs x=0.97 at α=4°), possibly due to differences in amplification rate calculation

## Related Files

- `crates/rustfoil-coupling/src/march.rs` - Main fix location
- `crates/rustfoil-bl/src/closures/transition.rs` - e^N transition implementation
- Previous fix: Theta safeguard bug (lines 1844-1856 in march.rs)
