---
title: Multi-Element Architecture & Decisions
sidebar_label: Multi-Element (ADR)
sidebar_position: 20
description: Architectural decisions for extending FlexFoil from single-airfoil to coupled multi-element high-lift configurations
keywords: [multi-element, high-lift, slat, flap, vane, coupled, panel method]
---

# Multi-Element Architecture & Decisions

Living record of the decisions behind extending FlexFoil from single-airfoil to
**coupled multi-element** (slat / main / vane / flap) aerodynamics. Future tasks
should read this before touching the multi-element path.

**Status:** in progress. Phase 1 (coupled inviscid core) implemented & validated.

## Guiding principles

1. **Do not make changes irrelevant to the multi-element goal.** No refactor for
   its own sake.
2. **Leave every touched area better organized than we found it** (boy-scout
   rule), but only where our own changes already land.

## Goal

A web-based design tool for a multi-element high-lift airfoil (reference: the
3-element STOL config — coved LS(1)-0417 main + NACA vane + NACA aft-flap on a
flap track). Two result panes: a drag polar and a Cp distribution. The geometry
engine and flap-track kinematics already exist in `flexfoil-ui/src/highlift/`.

## Decisions

### D1 — Build on the "faithful" solver stack, not the legacy one
The workspace has **two parallel solver lineages**:
- **Faithful (production):** `rustfoil-inviscid` (XFOIL-exact PSILIN panel
  method) → `rustfoil-bl` (closures) → `rustfoil-coupling` (Newton VII) →
  `rustfoil-xfoil` (orchestration). This is what the UI's `analyze_airfoil_faithful`
  uses.
- **Legacy:** `rustfoil-solver` carries its **own** `inviscid/influence.rs` and
  `viscous/viscal.rs`, plus the flow-visualization code (`smoke.rs`,
  `velocity.rs`). It also defines duplicate `FlowConditions` / `InviscidSolution`.

All multi-element work goes in the **faithful** stack. We deliberately **do not**
touch or extend `rustfoil-solver`'s duplicate inviscid/types — sidestepping that
debt rather than expanding it. (Retiring the legacy stack and unifying the
duplicated types are real but **off-goal** refactors; deferred.)

### D2 — Fidelity roadmap
- **(a) Coupled inviscid** — combined influence matrix + one Kutta condition per
  body. Correct interacting Cp/Cl. No drag. *(Phase 1 — done.)*
- **(b) Hybrid BL** — march each element's boundary layer on its *coupled* edge
  velocity (no wake/element confluence). Approximate drag → restores the drag
  polar. *(Phase 3.)*
- **(c) Full viscous multi-element** — confluent BL + wake-on-element coupling +
  global Newton (MSES-class). Research-grade; only if (b) proves insufficient.

A drag polar needs viscosity, so (b) is required for the headline pane.

### D3 — First cut is bodies-only (no wake panels)
Phase 1 omits wake panels. Consequence: trailing-downwash effects between
elements are not fully modeled inviscidly; near-field potential coupling is.

### D4 — Coupled inviscid system structure
Implemented in `crates/rustfoil-inviscid/src/multibody.rs`:
- Unknowns: `N_total` vortex strengths γ (all bodies) + `K` stream-function
  constants ψ₀ (one per body). System size `N_total + K`.
- One **Kutta condition per body**: γ(upper TE) + γ(lower TE) = 0.
- One **ψ₀ per body** (each body is its own streamline).
- Cross-body influence: evaluate body B's panels at a control point on body A via
  `psilin` with an **out-of-range node index** (the external/wake convention
  already used for wake points in the DIJ code). Self-blocks keep the
  node-index singularity handling.
- Per-body sharp-TE bisector rows preserved as in single-body GGCALC.
- **Reduces exactly to single-body when K=1** (guaranteed and tested).

### D5 — Validation strategy (no multi-element reference exists)
Neither XFOIL nor the bundled `mfoil` does multi-element (both single-element;
`mfoil`'s "gap" code is TE-gap/wake-gap, not slot gaps). So the repo's
"match-XFOIL-Fortran" parity culture does not apply here. Anchors used
(`crates/rustfoil-inviscid/tests/multibody.rs`):
1. **Single-body equivalence** — K=1 multi matches single-body to 1e-9.
2. **Far-apart recovery** — bodies 50c apart recover isolated Cl (<1e-3).
3. **Monotone decay** — interaction strength falls monotonically with separation.

External checks: an analytic biplane interference factor (optional, inviscid);
**30P-30N** experimental/published data (the real check, but *viscous* — belongs
after Phase 3).

### D6 — Overlap / geometry validity: rely on the human, no geometric check
The tool is human-in-the-loop and the UI **draws all elements**, so a contour
intersection (e.g. low-deploy nesting inside the cove) is **visually obvious** to
the user. A geometric pre-check would be defensive complexity guarding a state
the user can already see and avoid — so we **omit it** (YAGNI / wu-wei).

Rationale that this is safe without a gate:
- True intersection leaves the AIC non-singular and returns finite-but-wrong
  numbers, which the user correlates with the visibly-bad geometry.
- The rare exact-degenerate case trips the existing `SingularMatrix` error.
- The Web Worker isolates the solve, so even NaN/garbage merely fails to plot —
  it never crashes the UI.

No intersection test, no minimum-gap threshold. **Revisit only if** a headless /
batch multi-element API (no canvas, no human) becomes a goal — there, overlap is
not visually obvious and a documented precondition or check would be warranted.

### D7 — WASM / UI boundary (Phase 2)
- Add **one** `analyze_multi_element` WASM entry with a coherent
  geometry-in (per-element contours) / results-out (per-element + total Cp &
  forces) DTO + small marshaling helper — instead of cloning the per-function
  flat-array pattern. This is the only on-path organizational improvement we
  actively pursue (boy-scout, #3 from the repo review).
- Swap `flexfoil-ui/src/highlift/solve.ts` from N independent solves to one
  coupled solve; this removes the interim **effective-α** hack (free cleanup of
  code we're already touching).
- Until Phase 3, the drag-polar pane shows Cl-vs-α (inviscid has no drag).

## Code pointers
- Coupled inviscid solver: `crates/rustfoil-inviscid/src/multibody.rs`
  (`build_and_factorize_multi`, `FactorizedMultiSystem::solve_alpha`,
  `MultiInviscidSolution`).
- Validation: `crates/rustfoil-inviscid/tests/multibody.rs`.
- Geometry engine + UI: `flexfoil-ui/src/highlift/` (`geometry.ts`,
  `solve.ts` seam, `solveWorker.ts`, `preview.ts`).

## Explicitly deferred (off-goal)
Retiring the legacy `rustfoil-solver` stack; unifying duplicate
`FlowConditions`/`InviscidSolution`; splitting inverse design out of
`rustfoil-xfoil`; relocating `flap.rs` / `mfoil/`. Revisit only if a future
multi-element task genuinely lands in those areas.
