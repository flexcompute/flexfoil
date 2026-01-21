# RustFoil Viscous BL Implementation Tasks

## Overview

This folder contains 18 self-contained task files for implementing XFOIL's viscous-inviscid coupling in Rust. Each task is designed to be completed in a single Cursor chat session.

## Quick Start: Copy-Paste Prompts

### Wave 1 (1 chat)
```
@docs/tasks/TASK_01_TESTKIT.md
Implement this task completely. Create rustfoil-testkit with FORTRAN compilation utilities.
Verify with: cd crates/rustfoil-testkit/fortran && make
```

### Wave 2 (1 chat)
```
@docs/tasks/TASK_02_BL_CONSTANTS.md
Implement this task completely. Create rustfoil-bl with BLPAR constants.
Verify with: cargo build -p rustfoil-bl
```

### Wave 3 (7 PARALLEL chats)
Run these in separate Cursor windows simultaneously:

| Chat | Prompt |
|------|--------|
| 3A | `@docs/tasks/TASK_03_HKIN.md` - Implement HKIN closure |
| 3B | `@docs/tasks/TASK_04_HS.md` - Implement HSL/HST closures |
| 3C | `@docs/tasks/TASK_05_CF.md` - Implement CFL/CFT closures |
| 3D | `@docs/tasks/TASK_06_DISSIPATION.md` - Implement DIL/DIT/DILW |
| 3E | `@docs/tasks/TASK_07_DENSITY.md` - Implement HCT closure |
| 3F | `@docs/tasks/TASK_08_TRANSITION.md` - Implement DAMPL |
| 3G | `@docs/tasks/TASK_09_STATE.md` - Implement BlStation |

### Wave 4+ (sequential, 1 chat each)
```
@docs/tasks/TASK_10_EQUATIONS.md - Implement BLVAR/BLDIF
@docs/tasks/TASK_11_COUPLING_DIJ.md - Create rustfoil-coupling + QDCALC
@docs/tasks/TASK_12_COUPLING_NEWTON.md - Implement BLSYS
@docs/tasks/TASK_13_COUPLING_SOLVE.md - Implement BLSOLV
@docs/tasks/TASK_14_COUPLING_MARCH.md - Implement MRCHUE/MRCHDU
@docs/tasks/TASK_15_COUPLING_UPDATE.md - Implement UPDATE/UESET
@docs/tasks/TASK_16_VISCAL.md - Implement main solver
@docs/tasks/TASK_17_CLI.md - Add CLI commands
@docs/tasks/TASK_18_INTEGRATION_TESTS.md - Integration tests
```

For each, add: "Implement this task completely. Verify with the command in the task file."

---

## Full Prompts

See the Docusaurus page for complete copy-paste prompts:
`/docs-site/docs/rustfoil-implementation-plan.mdx`

Or view at: http://localhost:3000/docs/rustfoil-implementation-plan

---

## How to Use (Detailed)

1. **Start a new chat** for each task
2. **Reference the task file**: `@docs/tasks/TASK_XX_NAME.md`
3. **Tell Claude**: "Implement this task completely"
4. **Verify** the deliverables before moving to next task

## Task List

| # | File | Description | Dependencies |
|---|------|-------------|--------------|
| 01 | TASK_01_TESTKIT.md | Create rustfoil-testkit with FORTRAN harness | None |
| 02 | TASK_02_BL_CONSTANTS.md | Create rustfoil-bl with BLPAR constants | 01 |
| 03 | TASK_03_HKIN.md | Implement HKIN closure | 01, 02 |
| 04 | TASK_04_HS.md | Implement HSL/HST closures | 01-03 |
| 05 | TASK_05_CF.md | Implement CFL/CFT closures | 01-04 |
| 06 | TASK_06_DISSIPATION.md | Implement DIL/DIT/DILW closures | 01-05 |
| 07 | TASK_07_DENSITY.md | Implement HCT closure | 01-06 |
| 08 | TASK_08_TRANSITION.md | Implement DAMPL/TRCHEK2 | 01-07 |
| 09 | TASK_09_STATE.md | Implement BlStation struct | 01-08 |
| 10 | TASK_10_EQUATIONS.md | Implement BLVAR/BLDIF | 01-09 |
| 11 | TASK_11_COUPLING_DIJ.md | Create rustfoil-coupling with QDCALC | 01-10 |
| 12 | TASK_12_COUPLING_NEWTON.md | Implement BLSYS Newton builder | 11 |
| 13 | TASK_13_COUPLING_SOLVE.md | Implement BLSOLV solver | 12 |
| 14 | TASK_14_COUPLING_MARCH.md | Implement MRCHUE/MRCHDU | 13 |
| 15 | TASK_15_COUPLING_UPDATE.md | Implement UPDATE/UESET | 14 |
| 16 | TASK_16_VISCAL.md | Implement VISCAL main loop | 15 |
| 17 | TASK_17_CLI.md | Add viscous CLI commands | 16 |
| 18 | TASK_18_INTEGRATION_TESTS.md | Integration tests vs XFOIL | 17 |

## Critical Context for All Tasks

### File Locations
- **Rust workspace**: `/Users/harry/flexfoil-boundary-layer/crates/`
- **FORTRAN source**: `/Users/harry/flexfoil-boundary-layer/Xfoil/src/`
- **gfortran**: `/opt/homebrew/bin/gfortran`

### Key FORTRAN Files
- `xblsys.f` - Closure functions (HKIN, HSL, HST, CFL, CFT, etc.)
- `xbl.f` - BL marching (SETBL, MRCHUE, MRCHDU, UPDATE)
- `xsolve.f` - Block solver (BLSOLV)
- `xoper.f` - Main loop (VISCAL)
- `xpanel.f` - Panel influence (QDCALC, UESET)
- `BLPAR.INC` - Constants
- `XBL.INC` - State variables

### Testing Strategy
Every Rust function must be unit tested against FORTRAN:
1. Compile FORTRAN test harness with gfortran
2. Generate JSON reference data
3. Rust tests load JSON and compare with tolerance 1e-10

### Test Documentation Requirements (CRITICAL)

**Every completed task MUST document its tests** in `docs-site/docs/bl-implementation-progress.mdx`:

1. **Test command**: The exact `cargo test -p {crate}` command
2. **Result summary**: "X passed, Y failed"
3. **Test documentation** grouped by module/purpose. For EACH test explain:
   - **Why it's needed**: What problem does this test prevent? What failure mode does it catch?
   - **What it proves**: The specific assertion and what passing demonstrates

4. **FORTRAN reference data** (if applicable):
   - **Why we need it**: Validation strategy (ground-truth oracle)
   - **Parameter coverage rationale**: Physical meaning of test ranges
   - **Outputs captured**: What values/derivatives and why they matter
   - **Validation**: How data was verified before depending on it

Example entry:
```markdown
#### Test Results

**Command:** `cargo test -p rustfoil-bl`  
**Result:** 8 passed, 0 failed

##### Closure Tests (`closures::hkin` module)

HKIN converts H to kinematic Hk. Errors propagate to all downstream closures.

| Test | Why It's Needed | What It Proves |
|------|-----------------|----------------|
| `test_hkin_incompressible` | Sanity check: at M=0, correction vanishes | `hkin(H, 0).hk == H` within `1e-10` |
| `test_hkin_matches_fortran` | Ensures exact XFOIL compatibility | All 400 cases match within `1e-10` |

#### FORTRAN Reference Data

**Why we need this:** Ground-truth oracle ensures bit-exact compatibility with XFOIL.
```

## Quick Start

```bash
# Task 01: Create test infrastructure
# In Cursor: @docs/tasks/TASK_01_TESTKIT.md
# "Implement this task"

# Verify
cd /Users/harry/flexfoil-boundary-layer/crates/rustfoil-testkit/fortran
make run
```

## Estimated Effort

- Tasks 01-10 (closures + equations): ~2000 lines Rust, ~300 lines FORTRAN
- Tasks 11-15 (coupling): ~1200 lines Rust
- Tasks 16-18 (solver + CLI + tests): ~800 lines Rust

Total: ~4000 lines of Rust

## Troubleshooting

If a chat gets stuck:
1. Save progress (commit any completed code)
2. Start fresh chat
3. Reference the same task file
4. Continue from where you left off

If tests fail:
1. Check FORTRAN reference matches XFOIL behavior
2. Verify constants match BLPAR.INC
3. Check derivative calculations numerically

---

## Wave 3 Summary: Closure Functions (Tasks 03-09)

Tasks 03-09 form the **closure layer** of the boundary layer solver. These implement the empirical correlations that close the integral boundary layer equations. All seven tasks can be run in **parallel** since they have no interdependencies.

### Task 03: HKIN - Kinematic Shape Factor Transformation

**Purpose:** Converts shape factor H (δ*/θ) to kinematic shape factor Hk with compressibility correction.

**FORTRAN Source:** `xblsys.f` line 2276

**Key Formula:**
- `Hk = (H - 0.29·M²) / (1 + 0.113·M²)`
- At M=0: `Hk = H` (incompressible limit)

**Outputs:** `hk`, `∂Hk/∂H`, `∂Hk/∂M²`

**Why it matters:** Hk is the input to all other closures. Errors here propagate everywhere.

---

### Task 04: HS - Energy Shape Factor Closures

**Purpose:** Compute energy shape factor H* (or Hs) for laminar and turbulent flow.

**FORTRAN Source:** `xblsys.f` HSL (line 2327), HST (line 2388)

**Key Formulas:**
- **Laminar (HSL):** Piecewise correlation depending on whether Hk < 4.35 or Hk ≥ 4.35
- **Turbulent (HST):** Complex correlation with Rθ dependence (~90 lines)

**Outputs:** `hs`, `∂Hs/∂Hk`, `∂Hs/∂Rθ`, `∂Hs/∂M²`

**Why it matters:** Hs appears in the energy integral equation and dissipation calculations.

---

### Task 05: CF - Skin Friction Closures

**Purpose:** Compute skin friction coefficient Cf for laminar and turbulent boundary layers.

**FORTRAN Source:** `xblsys.f` CFL (line 2354), CFT (line 2483)

**Key Formulas:**
- **Laminar (CFL):** El Hady correlation - piecewise depending on Hk vs 5.5
- **Turbulent (CFT):** Power law with compressibility correction

**Outputs:** `cf`, `∂Cf/∂Hk`, `∂Cf/∂Rθ`, `∂Cf/∂M²`

**Why it matters:** Cf directly gives wall shear stress (drag). Also enters the momentum integral equation.

---

### Task 06: Dissipation Closures ✅ COMPLETE

**Purpose:** Compute dissipation coefficient 2CD/H* for laminar, turbulent, and wake regions.

**FORTRAN Source:** `xblsys.f` DIL (line 2290), DIT (line 2375), DILW (line 2308)

**Key Formulas:**
- **Laminar (DIL):** Uses `HDCON = 5·0.0003/(1 + 0.1·DLCON)` with Hk³ dependence
- **Turbulent (DIT):** Takes (Hs, Us, Cf, St) as inputs
- **Wake (DILW):** Includes internal HSL helper for laminar H* calculation

**Outputs:** `di`, `∂Di/∂Hk`, `∂Di/∂Rθ` (plus additional derivatives for turbulent)

**Why it matters:** Dissipation appears in the energy integral equation (kinetic energy loss).

---

### Task 07: HCT - Density Shape Factor

**Purpose:** Compute density thickness shape factor Hc for compressible flows.

**FORTRAN Source:** `xblsys.f` HCT (line 2514)

**Key Formula (Whitfield correlation):**
- `Hc = M² · (0.064/(Hk - 0.8) + 0.251)`
- At M=0: `Hc = 0` (incompressible limit)

**Outputs:** `hc`, `∂Hc/∂Hk`, `∂Hc/∂M²`

**Why it matters:** Hc accounts for density variations in compressible boundary layers.

---

### Task 08: Transition Prediction

**Purpose:** Implement e^n method for laminar-turbulent transition prediction.

**FORTRAN Source:** `xblsys.f` DAMPL (line 1981, ~100 lines)

**Key Components:**
- **DAMPL:** Computes amplification rate dN/dRθ using Drela-Giles correlation for Tollmien-Schlichting waves
- **TRCHEK2:** Checks if N ≥ Ncrit (typically Ncrit = 9 for free flight)

**Outputs:** `ax` (dN/dRθ), `∂ax/∂Hk`, `∂ax/∂θ`, `∂ax/∂Rθ`

**Why it matters:** Determines where the boundary layer transitions from laminar to turbulent. Critical for drag prediction accuracy.

---

### Task 09: BlStation State Structure

**Purpose:** Define the Rust struct that holds all boundary layer state at one station, matching XFOIL's `XBL.INC`.

**FORTRAN Source:** `XBL.INC` (73+ state variables)

**Key Components:**
- **Primary variables** (Newton unknowns): x, Ue, θ, δ*, √Cτ, N
- **Secondary variables** (computed by closures): H, Hk, Hs, Hc, Rθ, Cf, CD
- **Mode flags:** is_laminar, is_wake, is_turbulent
- **BlDerivatives struct:** All partial derivatives for Jacobian construction

**Special initializers:**
- `BlStation::new()` - Default values
- `BlStation::stagnation(ue, re)` - Hiemenz stagnation point solution

**Why it matters:** This struct is passed throughout the solver. Clean design here simplifies all downstream tasks.

---

### Wave 3 Summary Table

| Task | Module | Functions | Complexity | Status |
|------|--------|-----------|------------|--------|
| 03 | `closures::hkin` | `hkin` | Simple (5 lines) | Pending |
| 04 | `closures::hs` | `hs_laminar`, `hs_turbulent` | Medium/Complex | Pending |
| 05 | `closures::cf` | `cf_laminar`, `cf_turbulent` | Medium | Pending |
| 06 | `closures::dissipation` | `dissipation_laminar`, `dissipation_turbulent`, `dissipation_wake` | Medium | ✅ Complete |
| 07 | `closures::density` | `density_shape` | Simple (3 lines) | Pending |
| 08 | `transition` | `amplification_rate`, `check_transition` | Complex (~100 lines) | Pending |
| 09 | `state` | `BlStation`, `BlDerivatives` | Medium (struct design) | Pending |

### Dependencies After Wave 3

Once Wave 3 is complete, Task 10 (`BLVAR`/`BLDIF`) can proceed. It uses all the closures to compute secondary variables and integral equation residuals.
