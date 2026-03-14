# Turbulent BL Solver Diagnostic Plan

## Problem Statement

After transition, RustFoil's turbulent boundary layer solution diverges from XFOIL:
- Hk collapses to ~1.05 instead of settling at ~1.5
- theta errors grow large post-transition
- The solution is non-physical

## Root Cause Analysis

### How XFOIL handles the laminar→turbulent transition

At station 32 (first turbulent station, x=0.1408):

1. **Initial state**: Inherits laminar solution (δ*=1.196e-3, Hk=5.4)
2. **Direct mode check**: Computes `HKTEST = (δ* + rlx*vsrez[3])/(θ + rlx*vsrez[2])`
3. **Inverse mode trigger**: Since `HKTEST > htmax (2.5)`, switches to inverse mode
4. **Inverse constraint**: Sets `HTARG = HK1 - 0.15*(X2-X1)/T1` ≈ 2.5 (clamped to htmax)
5. **Newton solve with constraint**: VS2 row 3 = `[0, HK2_T2, HK2_D2, HK2_U2]`, residual = `HTARG - HK2`
6. **Result**: Over 6 iterations, δ* drops from 1.2e-3 → 6.1e-4, Hk settles at 2.5

### RustFoil's current broken behavior

1. **Artificial reset**: On transition detection, sets `delta_star = htmax * theta` BEFORE solving
2. **Skip inverse mode**: This hack bypasses the natural inverse mode mechanism
3. **Wrong initial state**: The turbulent solver starts from an artificial state, not the physical laminar values

## Diagnostic Steps

### Step 1: Trace first turbulent station (station 32)
- Add detailed iteration logging to `newton_solve_station` for station 32
- Log: initial (θ, δ*, ctau), each iteration (vsrez, relaxation, mode), final state
- Compare against XFOIL's iteration trace

### Step 2: Verify inverse mode is being triggered correctly
- Check that `HKTEST > htmax` triggers inverse mode
- Verify the inverse mode constraint equation matches XFOIL:
  - `VS2[3] = [0, HK2_T2, HK2_D2, HK2_U2]`
  - `VSREZ[3] = HTARG - HK2`

### Step 3: Verify HTARG calculation for turbulent flow
XFOIL's formula (from MRCHUE lines 725-730):
```fortran
C----------- turbulent case: relatively fast decrease in Hk downstream
             HTARG = HK1 - 0.15*(X2-X1)/T1
C---------- limit specified Hk to something reasonable
             HTARG = MAX( HTARG , HMAX )
```

RustFoil must use the same formula when in inverse mode for turbulent stations.

### Step 4: Remove the artificial delta* reset
The hack at lines 624-628 in march.rs should be removed:
```rust
// REMOVE THIS:
turb_prev.delta_star = config.htmax * turb_prev.theta;
turb_prev.hk = config.htmax;
```

Instead, pass the actual laminar solution to the turbulent Newton solve and let inverse mode do its job.

### Step 5: Verify ctau evolution via shear-lag equation
- Check that `vsrez[0]` (ctau update) is being computed correctly for turbulent flow
- Verify ctau clamps: `0.0000001 ≤ ctau ≤ 0.30`
- Trace ctau from initial 0.03 → final values at each turbulent station

### Step 6: Verify turbulent closures (Cf, Hs, Cd)
- Compare `cf_turbulent`, `hs_turbulent`, `dissipation_turbulent` outputs against XFOIL
- These should already be correct (tested for laminar), but verify for turbulent Hk range

## Implementation Order

1. **Remove artificial delta* reset** → let inverse mode handle transition
2. **Fix HTARG calculation** → use XFOIL's turbulent formula `HK1 - 0.15*(X2-X1)/T1`
3. **Add diagnostic tracing** → verify iteration-by-iteration match
4. **Verify shear-lag equation** → ensure ctau evolves correctly
5. **Test and validate** → compare against XFOIL turbulent stations 32-45

## Success Criteria

- Station 32: θ within 5% of XFOIL, Hk=2.5 (matching htmax clamp)
- Station 33-35: Hk decreasing toward ~1.5
- Station 36+: Hk stable at ~1.47-1.48
- ctau evolving from 0.03 → ~0.05 → ~0.04 (matching XFOIL pattern)

## XFOIL Reference Data (Station 32)

| Iteration | θ | δ* | dmax |
|-----------|---|----|----|
| 2 | 2.36e-4 | 1.07e-3 | 0.94 |
| 3 | 2.54e-4 | 8.97e-4 | 0.57 |
| 4 | 2.54e-4 | 6.36e-4 | 0.29 |
| 5 | 2.44e-4 | 6.09e-4 | 0.07 |
| 6 | 2.44e-4 | 6.10e-4 | 0.001 |
| 7 | 2.44e-4 | 6.10e-4 | 3e-6 (converged) |

Final: θ=2.44e-4, δ*=6.10e-4, Hk=2.50, ctau=0.0546
