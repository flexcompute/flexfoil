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
- **Rust workspace**: `/Users/harry/flexfoil/crates/`
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

## Quick Start

```bash
# Task 01: Create test infrastructure
# In Cursor: @docs/tasks/TASK_01_TESTKIT.md
# "Implement this task"

# Verify
cd /Users/harry/flexfoil/crates/rustfoil-testkit/fortran
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
