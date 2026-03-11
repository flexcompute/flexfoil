# XFOIL vs RustFoil Architecture Audit

## Executive Summary

RustFoil is no longer failing because the basic local boundary-layer physics are wildly different from XFOIL. The strongest evidence from the current repository points to a narrower conclusion:

- Early-station surface marching is already close to XFOIL.
- Transition detection in the stored trace gates is already close to XFOIL at the tested attached cases.
- The largest remaining parity gaps are architectural: outer coupling, TE/wake handling, panel-space gamma bookkeeping, and operating-variable sensitivity.

The current live solver state is materially better than earlier in the project:

- Full polar check: `CL avg error = 7.2%`, `CL max error = 41.8%`
- Full polar check: `CD ratio min/max/avg = 0.62x / 1.19x / 1.00x`
- `cargo test -p rustfoil-solver test_full_polar_comparison -- --nocapture` passes

The key conclusion is that the next parity work should focus less on isolated closure tweaks and more on matching XFOIL's process flow in four places:

1. outer Newton coupling / operating-variable sensitivity
2. full TE-to-wake process matching
3. panel-array-native stagnation relocation path
4. final drag from converged wake state rather than synthetic wake reconstruction

## Top-Level Pipeline Comparison

```mermaid
flowchart TD
    subgraph xfoilFlow [XFOILFlow]
        xViscal[VISCAL]
        xWake[XYWAKE_QWCALC_QISET]
        xStag[STFIND_IBLPAN_XICALC_IBLSYS_UICALC]
        xSetbl[SETBL]
        xSolve[BLSOLV_UPDATE]
        xPanel[QVFUE_GAMQV_STMOVE]
        xForces[CLCALC_CDCALC]
        xViscal --> xWake --> xStag --> xSetbl --> xSolve --> xPanel --> xForces
    end
    subgraph rustFlow [RustFoilFlow]
        rSetup[setup_from_body]
        rInit[initialize_surface_and_wake_stations]
        rMarch[march_surface]
        rNewton[march_mixed_du_GlobalNewtonSystem]
        rStmove[apply_stmove_like_xfoil]
        rForces[clcalc_compute_forces_syntheticWakeCd]
        rSetup --> rInit --> rMarch --> rNewton --> rStmove --> rForces
    end
```

## Side-by-Side Architecture

| Stage | XFOIL | RustFoil | Notes |
| --- | --- | --- | --- |
| Geometry / panel setup | `xfoil.f`, `xpanel.f` | `rustfoil-core`, `rustfoil-inviscid/src/geometry.rs` | Broadly comparable responsibility split |
| Inviscid base solve | `GGCALC`, `QISET`, `QDCALC` | `FactorizedSystem::solve_alpha`, `build_dij_with_wake` | RustFoil already has wake columns in `DIJ` |
| Stagnation / BL indexing | `STFIND`, `IBLPAN`, `XICALC`, `IBLSYS`, `UICALC` | `setup.rs`, `extract_surface_xfoil`, `find_stagnation_with_derivs` | RustFoil still reconstructs some panel-space data from split stations |
| Initial BL march | `MRCHUE` | `march_surface` | Close in intent; trace gates suggest early attached-region parity is good |
| Mixed remarch | `MRCHDU` | `march_mixed_du` | RustFoil currently skips wake stations here |
| Transition interval | `TRCHEK2`, `TRDIF` | `trchek2_full`, `trdif_full` | Much closer than before; not the dominant remaining system-level gap |
| TE join / first wake | `TESYS` + `VZ` chain-rule coupling | `apply_te_wake_interface`, `build_vz_block` | RustFoil has a reduced algebraic version, not the full XFOIL transformation |
| Global Newton solve | `SETBL` + `BLSOLV` | `build_global_system` + `solve_global_system` + `apply_global_updates` | Main structural mismatch remains the missing operating-variable sensitivity column |
| Stagnation relocation | `QVFUE -> GAMQV -> STMOVE` | `construct_gamma_from_stations*` -> `apply_stmove_like_xfoil` | RustFoil is sequence-similar but not panel-array-native |
| Final forces | `CLCALC`, `CDCALC` | `clcalc`, `compute_forces_two_surfaces`, synthetic wake CD | Main drag mismatch source is still architectural |

## Source-Level Call Order

### XFOIL

`VISCAL` in [`Xfoil-instrumented/src/xoper.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xoper.f) performs the canonical operating-point loop:

```2893:3039:/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xoper.f
SUBROUTINE VISCAL(NITER1)
...
CALL XYWAKE
CALL QWCALC
CALL QISET
...
CALL STFIND
CALL IBLPAN
CALL XICALC
CALL IBLSYS
CALL UICALC
...
CALL SETBL
CALL BLSOLV
CALL UPDATE
...
CALL QVFUE
CALL GAMQV
CALL STMOVE
CALL CLCALC(...)
CALL CDCALC
```

### RustFoil

`solve_viscous_two_surfaces()` in [`crates/rustfoil-solver/src/viscous/viscal.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs) is the closest Rust analogue:

```910:1186:/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs
pub fn solve_viscous_two_surfaces(
...
let upper_result = march_surface(&upper_arc, upper_ue, re, msq, &march_config, 1);
let lower_result = march_surface(&lower_arc, lower_ue, re, msq, &march_config, 2);
...
for iter in 0..max_newton_iter {
    ...
    let mut global_system =
        GlobalNewtonSystem::new(n_upper, n_lower, iblte_upper, iblte_lower);
    ...
    let current_gamma =
        construct_gamma_from_stations_qvfue(upper_stations, lower_stations, panel_x.len());
    ...
```

and later in the same function:

```1221:1586:/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs
march_mixed_du(upper_stations, re, msq, &march_config, 1);
march_mixed_du(lower_stations, re, msq, &march_config, 2);
...
global_system.build_global_system(...)
let deltas = solve_global_system(&mut global_system);
let update_result = apply_global_updates(...)
...
if let Some(new_ist) = apply_stmove_like_xfoil(...)
```

This is broadly the right shape, but several important XFOIL steps are either reduced or reconstructed differently.

## State Ownership and Dataflow Differences

### 1. Operating-variable sensitivity: `VDEL` mismatch

This is the biggest architectural gap.

RustFoil stores a single residual/solution vector per station:

```97:99:/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/global_newton.rs
/// RHS residual / solution vector (3 per station, 2 for RHS/solution)
/// VDEL[IV] = [res_third, res_mom, res_shape]
pub vdel: Vec<[f64; 3]>,
```

XFOIL fills two `VDEL` columns at assembly time: one for the residual and one for the operating-variable sensitivity path:

```368:378:/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xbl.f
IF(LALFA) THEN
 VDEL(1,2,IV) = VSR(1)*RE_CLMR + VSM(1)*MSQ_CLMR
ELSE
 VDEL(1,2,IV) = ...
ENDIF
VDEL(1,1,IV) = VSREZ(1) + ...
```

Why it matters:

- XFOIL is not just solving for local BL corrections. It is carrying the operating-variable sensitivity through the same assembled system.
- RustFoil can currently converge station states, but it cannot reproduce the same outer alpha/CL coupling architecture.
- This is a strong candidate for remaining asymmetry and robustness differences near small alpha and at difficult operating points.

### 2. Stagnation relocation still depends on reconstructed panel gamma

RustFoil reconstructs panel `gamma` from split BL stations:

```130:168:/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs
fn construct_gamma_from_stations_inner(
    upper_stations: &[BlStation],
    lower_stations: &[BlStation],
    n_panels: usize,
    lower_overwrites: bool,
) -> Vec<f64> {
    let mut gamma = vec![0.0; n_panels];
    ...
    gamma[idx] = station.u;
    ...
    gamma[idx] = -station.u;
```

XFOIL instead updates panel-space arrays directly through `QVFUE` and `GAMQV`, then calls `STMOVE` on those canonical arrays.

Why it matters:

- RustFoil’s stagnation relocation is process-similar but not dataflow-identical.
- Overwrite order and missing/unfilled panel behavior matter around symmetry and the stagnation neighborhood.
- This directly lines up with the remaining `-1 deg` lift asymmetry.

### 3. Wake ownership is split across two incompatible paths

RustFoil carries wake stations inside the lower-side coupled solve, but it does not trust them for final drag:

```1817:1868:/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs
// The lower surface stations already include wake stations ...
// Use those directly for Squire-Young CD ...
...
// Always use synthetic wake march from the combined TE state for CD computation.
...
let wake_stations_synth = march_wake(&wake_initial, &wake_x, &wake_ue, re, msq);
let cd_wake = compute_cd_from_wake(&wake_stations_synth, forces.cd_friction);
```

Why it matters:

- RustFoil currently has two wake concepts: the Newton-coupled wake and the synthetic drag wake.
- XFOIL’s `CDCALC` uses the actual converged wake state.
- As long as drag comes from a synthetic wake, RustFoil cannot be process-identical even if the surface BL is accurate.

## Validated Runtime Evidence

## What Already Matches Well

Two trace-based parity gates pass today:

- `cargo test -p rustfoil-solver test_alpha4_upper_station_windows_and_transition_gate -- --nocapture`
- `cargo test -p rustfoil-solver test_alpha10_alpha12_upper_station_windows_keep_ue_locked -- --nocapture`

These tests verify that:

- early upper-surface `Ue` remains locked to the stored XFOIL traces
- the `alpha=4` upper transition interval still lands on the same station/window
- attached-region upper-surface `Ue` remains close through `alpha=10` and `alpha=12`

Interpretation:

- RustFoil’s early marching and attached-region `Ue` propagation are not the main remaining problem.
- The biggest remaining differences are therefore downstream of those early gates: outer coupling, TE/wake handling, or panel-space bookkeeping.

## Live Full Polar

Current live full-polar check:

- `cargo test -p rustfoil-solver test_full_polar_comparison -- --nocapture`
- result: pass
- summary: `CD ratio min/max/avg = 0.62x / 1.19x / 1.00x`
- summary: `CL error min/max/avg = 0.0% / 41.8% / 7.2%`

Interpretation:

- Drag is now numerically decent overall, but still not process-identical because the final drag path is synthetic.
- Lift still has one dominant outlier regime rather than broad failure.

## First-Divergence Sweep

The structural sweep in `test_naca2412_first_divergence_sweep` is especially useful because it labels the earliest stage where the two solvers diverge for each alpha:

- Negative high alpha: often `Wake`
- Mid alpha attached / mildly loaded cases: often `Transition` or `Wake`
- Highest positive case in the sweep: `Relocation`

The emitted summary from the current run was:

```text
alpha    stage
-15      Wake
-12      Wake
-10      Wake
 -8      Transition
 -6      Transition
 -4      Transition
 -2      Transition
  0      Wake
  2      Wake
  4      Wake
  6      Transition
  8      Transition
 10      Wake
 12      Wake
 15      Relocation
```

Interpretation:

- This lines up with the source-level conclusions.
- The solver is not first diverging at geometry or very early `Ue` setup.
- The dominant divergence points are exactly the subsystems that remain architecturally non-identical: transition-process handling, wake evolution, and relocation/panel gamma bookkeeping.

## Severity-Ranked Mismatches

### Critical

#### 1. Missing outer operating-variable column in the global viscous solve

- XFOIL reference: [`Xfoil-instrumented/src/xbl.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xbl.f), [`Xfoil-instrumented/src/xsolve.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xsolve.f)
- RustFoil reference: [`crates/rustfoil-coupling/src/global_newton.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/global_newton.rs)
- Difference: XFOIL assembles both residual and operating-variable sensitivity into `VDEL`; RustFoil only solves the residual side.
- Why it matters: this changes the architecture of the outer viscous-inviscid correction loop.
- Fix order: before deeper tuning of small-alpha asymmetry.

#### 2. Wake is not fully marched through `MRCHDU`

- XFOIL reference: [`Xfoil-instrumented/src/xbl.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xbl.f)
- RustFoil reference: [`crates/rustfoil-coupling/src/march.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/march.rs)
- Difference: RustFoil explicitly skips wake stations in `march_mixed_du()`.
- Why it matters: the wake history is not evolving through the same remarch process as XFOIL.
- Fix order: alongside TE/wake parity work.

```2518:2528:/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/march.rs
for i in (start_idx + 1)..n {
    // Skip wake stations — they need TESYS (cross-surface join) which
    // requires both upper and lower TE states.
    if stations[i].is_wake {
        result.stations.push(stations[i].clone());
        continue;
    }
```

#### 3. Final drag still comes from a synthetic wake

- XFOIL reference: `CDCALC` in [`Xfoil-instrumented/src/xoper.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xoper.f) / [`Xfoil-instrumented/src/xfoil.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xfoil.f)
- RustFoil reference: [`crates/rustfoil-solver/src/viscous/viscal.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs), [`crates/rustfoil-coupling/src/wake.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/wake.rs)
- Difference: RustFoil discards the converged wake for drag and synthesizes a new one.
- Why it matters: this is a process mismatch by definition, even when numerical drag looks acceptable.
- Fix order: before declaring wake parity done.

### High

#### 4. Reduced `TESYS` / first-wake coupling

- XFOIL reference: [`Xfoil-instrumented/src/xbl.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xbl.f)
- RustFoil reference: [`crates/rustfoil-coupling/src/global_newton.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/global_newton.rs)
- Difference: RustFoil uses hand-built algebraic TE/wake coupling rather than the full XFOIL chain-rule structure.
- Why it matters: upper/lower TE sensitivity into the first wake station is still structurally weaker than in XFOIL.
- Fix order: pair with wake remarch work.

#### 5. Stagnation relocation is sequence-similar but not array-native

- XFOIL reference: `QVFUE -> GAMQV -> STMOVE` in [`Xfoil-instrumented/src/xoper.f`](/Users/harry/flexfoil-boundary-layer/Xfoil-instrumented/src/xoper.f)
- RustFoil reference: [`crates/rustfoil-solver/src/viscous/viscal.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-solver/src/viscous/viscal.rs)
- Difference: RustFoil rebuilds `gamma` from split BL stations and then estimates stagnation derivatives from that reconstruction.
- Why it matters: this is the most plausible remaining source of the persistent `-1 deg` lift asymmetry.
- Fix order: immediately after or in parallel with operating-variable sensitivity work.

#### 6. Wake geometry and recovery law remain synthetic

- XFOIL reference: wake panels, `WGAP`, and wake velocities from `XYWAKE`, `QWCALC`, `UICALC`
- RustFoil reference: [`crates/rustfoil-coupling/src/wake.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/wake.rs)
- Difference: RustFoil uses geometric wake spacing and empirical `wake_edge_velocity()`.
- Why it matters: this affects both wake growth and drag extrapolation.
- Fix order: after converged-wake drag ownership is corrected.

### Medium

#### 7. Transition kernel is close, but transition ownership is still split

- XFOIL reference: `TRCHEK2`, `TRDIF`, `SETBL`
- RustFoil reference: [`crates/rustfoil-coupling/src/march.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/march.rs), [`crates/rustfoil-coupling/src/global_newton.rs`](/Users/harry/flexfoil-boundary-layer/crates/rustfoil-coupling/src/global_newton.rs)
- Difference: local transition handling is now much closer, but the transition interval is still being fed by a solver architecture that differs upstream and downstream.
- Why it matters: this likely explains why trace gates pass in attached windows while some global small-alpha asymmetry remains.
- Fix order: after the larger architecture mismatches above.

## Recommended Fix Order

1. Extend the global coupled solve to carry the operating-variable sensitivity column like XFOIL `VDEL(*,2,IV)`.
2. Unify TE/wake handling so `MRCHDU`, `TESYS`, and final drag all use the same converged wake state.
3. Replace reconstructed stagnation-driving `gamma` with a panel-array-native `QVFUE/GAM` flow wherever possible.
4. Re-run the small-alpha symmetric cases and only then revisit residual transition asymmetry.
5. After those changes, tighten remaining wake geometry details such as `ANTE/WGAP` ownership and wake-spacing/velocity extraction.

## Bottom Line

RustFoil is no longer "mystically different" from XFOIL everywhere. The architecture audit points to a smaller set of high-leverage mismatches:

- missing outer coupling dimension
- incomplete TE/wake process matching
- panel-space stagnation bookkeeping differences
- synthetic final wake for drag

That is good news, because it means the parity problem is becoming narrower and more mechanical. The next wins are likely to come from matching XFOIL's dataflow, not from inventing new local closure tweaks.
