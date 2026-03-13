# XFOIL vs RustFoil Remaining Architecture Differences

## Executive Summary

RustFoil is no longer dominated by obviously wrong local BL closures. The remaining gaps are now mostly structural:

- The wake march is much closer to XFOIL than before, because the dedicated reduced wake Newton path has been replaced by the shared station solver.
- Canonical row/state ownership is materially better than before: wake geometry, `uinv`, `uinv_a`, operating sensitivity, and stagnation rollback are now more canonical-state-driven.
- The largest remaining architectural difference is the adapter boundary between the canonical XFOIL-like state and the global Newton/update layer, which still operates on split `Vec<BlStation>` views plus side arrays instead of operating directly on canonical rows/state.

Current broad status from the latest full sweep:

- `cargo test -p rustfoil-solver test_full_polar_comparison -- --nocapture`
- `CD ratio: min=0.11x, max=1.11x, avg=0.57x`
- `CL error: min=0.4%, max=22.1%, avg=7.4%`

Current targeted status for the sharp-TE debug case:

- `NACA0012 @ Re=1e6, alpha=2`
- `CL = 0.2225`
- `CD = 0.002663`

## What Has Already Been Architecturally Fixed

- Wake seeding is no longer driven by a heuristic recovery law. The pure solver now uses XFOIL-like wake geometry and inviscid wake-velocity ownership in `crates/rustfoil-solver/src/viscous/setup.rs`.
- Wake marching no longer uses a separate reduced Newton owner. `crates/rustfoil-coupling/src/wake.rs` now delegates wake marching into the shared station Newton machinery in `crates/rustfoil-coupling/src/march.rs`.
- Canonical rows now own more of the state that XFOIL rows own:
  - wake `y` geometry
  - `uinv`
  - `uinv_a`
  - row-derived wake arrays
- Newton rollback is now canonical-state-driven rather than juggling split `upper_ue_inv`, `lower_ue_inv`, and a shadow stagnation owner in parallel.

These fixes removed several real structural mismatches and measurably improved parity, especially in the wake/drag path.

## Remaining Architecture Differences

### 1. Global Newton/update still runs on split station views, not canonical rows

Files:

- `crates/rustfoil-solver/src/viscous/viscal.rs`
- `crates/rustfoil-coupling/src/global_newton.rs`

Current RustFoil architecture:

- `viscal` repeatedly projects canonical state out to `upper_stations` / `lower_stations`
- `GlobalNewtonSystem::build_global_system()` consumes `&[BlStation]`
- `apply_global_updates()` mutates `&mut [BlStation]`
- the canonical state is then refreshed from those split station vectors again

Why this still differs from XFOIL:

- In XFOIL, the canonical owner is the panel/BL array state itself. `SETBL`, `BLSOLV`, `UPDATE`, `QVFUE`, `GAMQV`, and `STMOVE` all operate on that one owner.
- In RustFoil, canonical ownership has improved, but the core Newton/update path still treats split station vectors as the primary mutable representation.

Why it matters:

- This is now the biggest remaining adapter boundary in the codebase.
- Any mismatch in station-view projection or writeback can still alter branch selection, especially around wake rows, stagnation, and transition ownership.

Leverage: `high`

### 2. Global Newton still depends on external side arrays for `Ue` bookkeeping

Files:

- `crates/rustfoil-coupling/src/global_newton.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`

Current RustFoil architecture:

- `build_global_system()` takes `ue_current_upper/lower`, `ue_inviscid_upper/lower`, `ue_from_mass_upper/lower`, and `ue_operating_upper/lower` as external slices
- `apply_global_updates()` and `preview_global_update_ue()` likewise depend on external `Ue` side arrays

Why this still differs from XFOIL:

- XFOIL’s update path uses its canonical array state directly for `USAV`, `UINV`, `UEDG`, `QVIS`, `GAM`, and operating sensitivity
- RustFoil still passes a collection of side arrays through the Newton/update interface rather than letting the canonical state own the full `Ue` story

Why it matters:

- It keeps the coupling layer more like a functional transform over slices than a direct mutation of canonical XFOIL-like state
- This is the main reason the state owner is still “transitional” rather than authoritative

Leverage: `high`

### 3. Stagnation-derivative setup still uses a local proxy instead of canonical panel arrays

Files:

- `crates/rustfoil-solver/src/viscous/viscal.rs`

Current RustFoil architecture:

- `viscal` computes `sst_go` / `sst_gp` from a temporary two-point proxy built from first upper/lower post-stagnation stations

Why this still differs from XFOIL:

- XFOIL computes stagnation derivatives from its canonical panel-space circulation/velocity arrays during `STFIND`
- RustFoil now stores canonical stagnation metadata, but this derivative setup is still reconstructed locally in `viscal`

Why it matters:

- It is still a non-canonical preconditioner for the global Newton solve
- Small stagnation derivative mismatches can alter leading-edge coupling behavior and relaxation behavior in difficult cases

Leverage: `medium`

### 4. Canonical state still does not own the full wake-panel inviscid basis

Files:

- `crates/rustfoil-solver/src/viscous/state.rs`
- `crates/rustfoil-solver/src/viscous/setup.rs`
- `crates/rustfoil-xfoil/src/state.rs`
- `crates/rustfoil-xfoil/src/wake_panel.rs`

Current RustFoil architecture:

- canonical state owns `wake_x`, `wake_y`, `wake_s`, `wake_uedg`, and `wake_mass`
- it does not own wake normals, wake panel angles, wake `qinvu_0`, wake `qinvu_90`, wake `qinv`, or wake `qinv_a`

Why this still differs from XFOIL:

- `XfoilState` owns the full wake-panel basis and derived wake inviscid arrays
- the pure solver’s canonical owner is still a reduced wake owner compared with XFOIL’s canonical wake state

Why it matters:

- repeated wake rebuilds are still less panel-native than XFOIL
- this limits how faithfully RustFoil can mirror `XYWAKE`, `QWCALC`, `QISET`, and wake-array reuse across the full operating-point loop

Leverage: `medium`

### 5. Friction drag still includes non-XFOIL filtering/clamping heuristics

Files:

- `crates/rustfoil-solver/src/viscous/forces.rs`

Current RustFoil architecture:

- `compute_friction_drag()` and `compute_row_friction_drag()` still:
  - skip low-`R_theta` rows
  - clamp `Cf`
  - skip numerically small/invalid rows

Why this still differs from XFOIL:

- XFOIL’s `CDCALC` integrates `TAU` directly from `IBL=3..IBLTE`
- it does not apply these Rust-side safety filters as part of the nominal architecture

Why it matters:

- This is less likely than the Newton/update boundary to be the dominant current error source
- but it is still a structural difference in the final force path

Leverage: `medium`

### 6. Legacy single-surface force helpers remain non-faithful utilities

Files:

- `crates/rustfoil-solver/src/viscous/forces.rs`

Current RustFoil architecture:

- `compute_forces()` and some older helper paths still contain simplified or fallback logic that is not XFOIL-faithful

Why this still differs from XFOIL:

- XFOIL’s canonical operating-point path is two-surface plus converged wake plus `CLCALC/CDCALC`
- these helper paths still exist as alternate utility logic

Why it matters:

- They are not the main current parity path
- but they still leave architectural ambiguity in the codebase

Leverage: `low` for current parity, `medium` for codebase clarity

### 7. Canonical state is still explicitly transitional rather than the sole owner

Files:

- `crates/rustfoil-solver/src/viscous/state.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`

Current RustFoil architecture:

- `XfoilLikeViscousState` is the main owner for more data than before
- but its own comments and call patterns still reflect an additive/transitional design
- `viscal` still orchestrates repeated projection to/from station views

Why this still differs from XFOIL:

- `XfoilState` in the faithful path is not transitional; it is the canonical mutable owner used by the full loop

Why it matters:

- This is the umbrella structural gap behind several smaller adapter mismatches
- the remaining parity work becomes easier once the canonical state is truly authoritative end-to-end

Leverage: `high`, but broad in scope

## Highest-Leverage Next Steps

If the goal is to continue matching architecture rather than tuning outputs, the next work should be ordered like this:

1. Collapse the Newton/update adapter boundary so `build_global_system()` and `apply_global_updates()` operate on canonical rows/state more directly.
2. Move `Ue` side arrays (`ue_current`, `ue_inviscid`, `ue_from_mass`, `ue_operating`) behind canonical-state-owned views instead of passing them around as free slices.
3. Make stagnation derivative setup canonical/panel-native rather than locally reconstructed in `viscal`.
4. Expand canonical wake ownership to include the wake-panel inviscid basis (`wake_qinvu_0`, `wake_qinvu_90`, `wake_qinv`, `wake_qinv_a`, normals/angles).
5. Remove the remaining force-path safety heuristics only after the canonical Newton/update owner is more direct and stable.

## Bottom Line

RustFoil is now much closer to XFOIL in the wake march and row/state ownership than it was before. The remaining differences are no longer “mystery physics” problems; they are mostly concentrated in one architectural seam:

- the global Newton/update process still treats split `BlStation` vectors and external `Ue` slices as the working representation,
- while XFOIL works directly on one canonical owner state throughout the operating-point loop.

That is the main remaining architectural frontier.
