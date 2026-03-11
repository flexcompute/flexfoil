# XFOIL State Parity Implementation Plan

## Goal

Make RustFoil's viscous/inviscid coupling architecture match XFOIL's state model as closely as possible, so parity work is no longer fighting architectural translation gaps.

The target end state is:

- one canonical viscous state object that owns the XFOIL-equivalent arrays
- split upper/lower station views derived from that canonical state, not co-owning it
- one panel-space circulation path for `QVIS/GAM/GAM_A`
- one wake ownership path
- one `UPDATE` path with live operating-variable correction
- one force-evaluation path

This file is intended as an execution backlog for another agent.

## Current Architectural Gaps

### 1. State ownership is still station-first

RustFoil still treats split `upper_stations` / `lower_stations` as the primary iterate in:

- `crates/rustfoil-solver/src/viscous/viscal.rs`
- `crates/rustfoil-solver/src/viscous/setup.rs`
- `crates/rustfoil-coupling/src/march.rs`

Panel arrays such as `QVIS/GAM` are still reconstructed views in:

- `crates/rustfoil-solver/src/viscous/circulation.rs`

XFOIL is the opposite: panel and BL arrays coexist as one canonical shared state.

### 2. Operating-variable plumbing exists, but is not yet the real XFOIL mode system

RustFoil now carries a second RHS through:

- `crates/rustfoil-coupling/src/global_newton.rs`

But the top-level solver still runs fixed-alpha with:

- `crates/rustfoil-solver/src/viscous/viscal.rs`

and `dac = 0.0`.

XFOIL's `UPDATE` uses the second RHS as a live control variable path.

### 3. Wake ownership is still split

Wake geometry, wake-aware DIJ, BL wake initialization, wake remarch, and drag ownership are still spread across:

- `crates/rustfoil-inviscid/src/system.rs`
- `crates/rustfoil-solver/src/viscous/setup.rs`
- `crates/rustfoil-coupling/src/wake.rs`
- `crates/rustfoil-coupling/src/march.rs`
- `crates/rustfoil-solver/src/viscous/forces.rs`

XFOIL centralizes wake behavior much more tightly through the panel layer and shared BL state.

### 4. Marching and stagnation relocation still reconstruct state

RustFoil's `STMOVE`-like flow still rebuilds surface templates and copies state around in:

- `crates/rustfoil-solver/src/viscous/viscal.rs`

XFOIL instead shifts the canonical BL arrays in place.

## Reference Files

### XFOIL references

- `Xfoil-instrumented/src/xpanel.f`
- `Xfoil-instrumented/src/xbl.f`
- `Xfoil-instrumented/src/xfoil.f`
- `Xfoil-instrumented/src/xblsys.f`

### RustFoil primary files

- `crates/rustfoil-solver/src/viscous/setup.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`
- `crates/rustfoil-solver/src/viscous/circulation.rs`
- `crates/rustfoil-solver/src/viscous/forces.rs`
- `crates/rustfoil-coupling/src/global_newton.rs`
- `crates/rustfoil-coupling/src/march.rs`
- `crates/rustfoil-coupling/src/wake.rs`
- `crates/rustfoil-inviscid/src/system.rs`
- `crates/rustfoil-bl/src/equations.rs`
- `crates/rustfoil-bl/src/closures/transition.rs`

## Execution Rules

- Do not try to complete this in one giant refactor.
- After each phase, re-run the focused parity tests before moving on.
- Prefer structure-first parity before output-first tuning.
- Do not preserve redundant fallback paths unless they are temporarily needed to keep the code compiling during a phase.
- If a temporary adapter is introduced, mark it clearly and delete it in a later phase.

## Phase 0: Lock the Test Harness

### Objective

Make sure there are regression gates for the architectural transition.

### Tasks

1. Keep and use these existing tests as required gates:
   - `cargo test -p rustfoil-solver test_full_polar_comparison --test xfoil_viscous_comparison -- --nocapture`
   - `cargo test -p rustfoil-solver test_full_sweep_cd_cdp_cm_gate --test xfoil_viscous_comparison -- --nocapture`
   - `cargo test -p rustfoil-solver test_naca2412_first_divergence_sweep --test structural_viscous_parity -- --nocapture`
   - `cargo test -p rustfoil-coupling test_apply_global_updates_uses_operating_correction -- --nocapture`

2. Add structural snapshot tests for:
   - panel `QVIS/GAM` ownership near stagnation
   - `IST/SST` and first BL rows after setup
   - first wake station row and final wake row
   - `DAC`, `UNEW`, `U_AC`, and `RLX` once the mode system is active

### Files

- `crates/rustfoil-solver/tests/xfoil_viscous_comparison.rs`
- `crates/rustfoil-solver/tests/structural_viscous_parity.rs`
- `crates/rustfoil-coupling/src/global_newton.rs`

### Exit Criteria

- All baseline parity tests are green before starting the deeper refactor.

## Phase 1: Introduce a Canonical XFOIL-Style State Object

### Objective

Create one authoritative state owner for the viscous loop.

### Design Target

Add a new module and struct, e.g.:

- `crates/rustfoil-solver/src/viscous/state.rs`
- `XfoilLikeViscousState`

This state should own XFOIL-equivalent arrays and metadata:

- stagnation metadata: `ist`, `sst`, `sst_go`, `sst_gp`
- BL topology: `nbl_upper`, `nbl_lower`, `iblte_upper`, `iblte_lower`
- surface maps: `ipan`, `vti`, `isys`
- inviscid panel arrays: `qinv`, `qinv_a`, maybe `uinv`
- viscous panel arrays: `qvis`, `gam`, `gam_a`
- BL arrays by `(ibl, is)`: `theta`, `dstr`, `ctau_or_ampl`, `uedg`, `mass`
- wake arrays and wake geometry
- transition metadata
- operating-mode/update scratch: `u_new`, `u_ac`, `dac`, `rlx`

### Tasks

1. Create the canonical state type.
2. Add indexing helpers so code can work in XFOIL-like `(ibl, is)` coordinates.
3. Add adapter methods to derive temporary `BlStation` views from the canonical state.
4. Add adapter methods to write `BlStation` edits back into the canonical state.

### Files

- New: `crates/rustfoil-solver/src/viscous/state.rs`
- Update exports from the viscous module root
- Temporary adapters in `crates/rustfoil-solver/src/viscous/viscal.rs`

### Exit Criteria

- The solver can instantiate and carry a canonical state object without changing physics yet.

## Phase 2: Move Setup to Populate Canonical State Directly

### Objective

Stop setup from producing independently-owned split station vectors as the primary representation.

### Tasks

1. Refactor setup flow to populate canonical arrays directly:
   - `STFIND` equivalent
   - `IPAN/VTI` equivalent
   - `XICALC/UICALC` equivalent seed data
   - wake append on lower side only

2. Convert these functions from "owners" to "initializers/helpers":
   - `setup_from_body()`
   - `setup_from_coords()`
   - `extract_surface_xfoil()`
   - `initialize_surface_stations_with_panel_idx()`
   - `initialize_wake_bl_stations()`

3. Keep split station construction only as an adapter for subsystems not yet ported.

### Files

- `crates/rustfoil-solver/src/viscous/setup.rs`
- `crates/rustfoil-inviscid/src/system.rs`

### Exit Criteria

- Setup creates canonical XFOIL-like state first.
- Split station vectors are derived, not authored, by setup.

## Phase 3: Make Panel Arrays Primary Throughout the Loop

### Objective

Make `QVIS/GAM/GAM_A` continuously owned, not reconstructed ad hoc.

### Tasks

1. Move panel array ownership into canonical state.
2. Replace on-demand projection as primary logic in:
   - `project_panel_qvis_from_two_surfaces()`
   - `project_panel_gamma_from_two_surfaces()`

3. Keep projections only as temporary adapters if needed.
4. Make `STMOVE` and `CLCALC` consume canonical `gam`.
5. Remove or demote alternate CL/circulation fallback paths.

### Files

- `crates/rustfoil-solver/src/viscous/circulation.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`
- `crates/rustfoil-solver/src/viscous/forces.rs`

### Exit Criteria

- `QVIS/GAM` exist continuously in state during the Newton loop.
- `CLCALC` reads canonical panel `gam`.

## Phase 4: Port `STMOVE` as a True State Mutation

### Objective

Stop rebuilding surface templates and copying station state around.

### Tasks

1. Refactor `apply_stmove_like_xfoil()` to mutate canonical state in place.
2. Match XFOIL flow:
   - update `IST/SST`
   - rebuild `IPAN/VTI`
   - refresh inviscid split values from panel arrays
   - shift BL rows in place
   - preserve wake and transition rows consistently

3. Remove temporary station-template rebuild logic after parity is proven.

### Files

- `crates/rustfoil-solver/src/viscous/viscal.rs`
- new or existing helpers in `crates/rustfoil-solver/src/viscous/state.rs`

### Exit Criteria

- Stagnation relocation no longer reconstructs the iterate from station vectors.

## Phase 5: Make `UPDATE` a Literal XFOIL-Style Mode-Dependent Path

### Objective

Turn the second RHS plumbing into the actual XFOIL operating-variable architecture.

### Tasks

1. Add an operating mode abstraction to the viscous config:
   - prescribed alpha
   - prescribed CL

2. Teach the canonical state/update path to own:
   - `u_new`
   - `u_ac`
   - `dac`
   - mode-dependent operating column meaning

3. Refactor `build_operating_rhs()` so the second RHS matches XFOIL mode semantics, not always alpha sensitivity.
4. Activate live `dac` in the top-level viscous loop.
5. Add `DALMAX/DCLMAX` style limits before the main `DN1..DN4` relaxation logic.

### Files

- `crates/rustfoil-coupling/src/global_newton.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`
- `crates/rustfoil-solver/src/viscous/config.rs`

### Exit Criteria

- The fixed-CL branch exists and uses real `DAC`.
- The fixed-alpha branch matches XFOIL's operating-variable logic shape.

## Phase 6: Unify Wake Ownership

### Objective

Ensure there is only one wake geometry/state source.

### Tasks

1. Make wake geometry, wake inviscid velocity, and wake-aware `DIJ` come from one canonical wake representation.
2. Match XFOIL's first wake node semantics:
   - first wake node constraints
   - first wake DIJ row behavior
   - first wake `Ue` ownership

3. Refactor or remove duplicated synthetic wake setup paths in:
   - `generate_wake_positions()`
   - `wake_edge_velocity()`
   - `march_wake()` fallback use

4. Ensure wake BL state is authored from canonical wake data, not re-created from separate heuristics.

### Files

- `crates/rustfoil-inviscid/src/system.rs`
- `crates/rustfoil-coupling/src/wake.rs`
- `crates/rustfoil-solver/src/viscous/setup.rs`
- `crates/rustfoil-coupling/src/march.rs`
- `crates/rustfoil-solver/src/viscous/forces.rs`

### Exit Criteria

- DIJ wake, marched wake, and drag wake all refer to the same wake representation.

## Phase 7: Port Marching to Canonical `(ibl, is)` State

### Objective

Move marching logic off independent `Vec<BlStation>` ownership.

### Tasks

1. Refactor `march_surface()` and `march_mixed_du()` to operate on canonical BL rows.
2. Use XFOIL-like indexing and first-wake semantics directly.
3. Make transition logic read/write canonical row state, not detached station structs.
4. Keep temporary `BlStation` wrappers only where required by lower-level equation APIs.

### Files

- `crates/rustfoil-coupling/src/march.rs`
- `crates/rustfoil-bl/src/equations.rs`
- `crates/rustfoil-bl/src/closures/transition.rs`

### Exit Criteria

- Marching is state-centric, not vector-centric.

## Phase 8: Collapse Force Ownership to One Path

### Objective

Make force computation follow XFOIL's singular architecture.

### Tasks

1. Make `CL` and `CM` come only from canonical `CLCALC`-style integration on panel `gam`.
2. Make `CD` come only from the converged wake/BL end state.
3. Remove synthetic wake drag fallback paths once parity is stable.
4. Ensure `compute_forces_two_surfaces()` is either simplified heavily or replaced by a canonical-state force module.

### Files

- `crates/rustfoil-solver/src/viscous/forces.rs`
- `crates/rustfoil-solver/src/viscous/viscal.rs`

### Exit Criteria

- One authoritative force path remains.

## Phase 9: Delete Transitional Adapters

### Objective

Remove temporary scaffolding introduced to bridge old and new architectures.

### Candidates for deletion or demotion

- station-first projection helpers once no longer needed as primary logic
- duplicate CL fallback logic
- synthetic wake fallback paths
- any rebuild-from-template stagnation move path
- any duplicated lower-surface sign conversion paths

### Exit Criteria

- The architecture reads like a direct structured Rust port of XFOIL's state flow.

## Recommended Task Breakdown For Another Agent

### Task A: canonical state skeleton

Create `XfoilLikeViscousState`, indexing helpers, and adapter views.

### Task B: setup migration

Refactor setup to populate canonical state directly, preserving current tests.

### Task C: panel ownership migration

Move `QVIS/GAM` and `CLCALC` to canonical ownership.

### Task D: `STMOVE` in-place mutation

Replace template rebuild/copy logic with in-place state mutation.

### Task E: live mode-dependent `UPDATE`

Implement true `DAC` behavior and operating-mode switching.

### Task F: wake ownership unification

Make wake geometry, wake BL state, and wake drag share one state source.

### Task G: marching migration

Move remarch/direct march to canonical `(ibl, is)` state.

### Task H: cleanup and adapter removal

Delete transitional scaffolding after parity tests are stable.

## Validation Gates After Each Major Phase

Run:

```sh
cargo test -p rustfoil-coupling test_apply_global_updates_uses_operating_correction -- --nocapture
cargo test -p rustfoil-solver test_full_polar_comparison --test xfoil_viscous_comparison -- --nocapture
cargo test -p rustfoil-solver test_full_sweep_cd_cdp_cm_gate --test xfoil_viscous_comparison -- --nocapture
cargo test -p rustfoil-solver test_naca2412_first_divergence_sweep --test structural_viscous_parity -- --nocapture
```

And after substantive edits:

```sh
cargo test -p rustfoil-solver test_initialize_surface_stations_with_panel_idx_virtual_anchor_zero_velocity -- --nocapture
cargo test -p rustfoil-solver test_build_panel_gamma_respects_write_order -- --nocapture
```

## Success Criteria

This project is done when:

- canonical state ownership mirrors XFOIL's architecture
- `STMOVE`, `QVFUE`, `GAMQV`, `UPDATE`, and wake handling operate on that state directly
- no synthetic/fallback architecture remains in the main solve path
- `CL` and `CD` parity improvements come from the same state architecture as XFOIL, not from patching translated views

