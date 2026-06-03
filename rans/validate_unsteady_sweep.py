#!/usr/bin/env python3
"""Validate the "unsteady-as-steady" α-sweep against independent steady runs.

Builds the mesh ONCE from a case.json, then on that same mesh runs:
  • one UNSTEADY case that marches through ``--alphas`` (one physical step per α,
    a huge time step, ``--max-pseudo`` pseudo-iters/step, a UDD driving alphaAngle), and
  • one STEADY case per α (the trusted reference).
Then it compares the per-physical-step CL/CD of the unsteady run to the steady runs.

This answers the two open questions before we wire the sweep into the server/UI:
  (1) does the solver honor the UDD ``alphaAngle`` control?  (2) is huge-Δt unsteady
  numerically equivalent to steady?

Build + preprocess run anywhere (self-test). The SOLVE needs a non-sandboxed shell:

    ~/flexcompute/compute/.venv/bin/python rans/validate_unsteady_sweep.py \
        --case /tmp/rans_server/run1/case.json --alphas 0,4 --solve \
        2>&1 | tee /tmp/unsteady_val/run.log

Without --solve it only builds the cases (validates the Flow360 JSON is accepted).
"""
from __future__ import annotations

import argparse
import json
import shutil
import sys
import time
from functools import partial
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rans import case as _case          # noqa: E402
from rans import contours as _contours  # noqa: E402
from rans import mesh as _mesh          # noqa: E402
from rans.config import CaseConfig      # noqa: E402
from rans.env import make_env           # noqa: E402

MESH = "mesh.cgns"


def boundary_lists(cfg: CaseConfig):
    elem = [e.name for e in cfg.elements]
    walls = [f"fluid/{e.name}" for e in cfg.elements if e.is_wall]
    bnd = (["fluid/farfield"] + [f"fluid/{n}" for n in elem]
           + ["fluid/symmetry1", "fluid/symmetry2"])
    return walls, bnd


def build_mesh(cfg: CaseConfig, out: Path, mesher_bin: Path, env, find) -> None:
    out.mkdir(parents=True, exist_ok=True)
    _contours.write_contours(cfg, out / "contours.txt")
    _mesh.run_mesher(out / "contours.txt", out / "mesh2d.vtk", mesher_bin, env, cwd=out)
    _mesh.write_volume_msh(out / "mesh2d.vtk", cfg, out / "mesh.msh")
    _mesh.gmsh_to_cgns(out / "mesh.msh", out / MESH, find("flow360gmshtocgns"), env)


def prep_and_solve(workdir: Path, mesh_src: Path, cfg, walls, bnd, find, env,
                   *, sim_builder, solve: bool, gpu: int) -> dict:
    workdir.mkdir(parents=True, exist_ok=True)
    shutil.copy(mesh_src, workdir / MESH)
    t = {}
    s = time.perf_counter()
    _case.preprocess(workdir, MESH, find, env, cfg=cfg, wall_names=walls,
                     boundary_names=bnd, timings=t, sim_builder=sim_builder)
    t["preprocess_total"] = round(time.perf_counter() - s, 2)
    if not solve:
        return {"timing": t}
    from rans import solve as _solve
    s = time.perf_counter()
    _solve.run_solver(workdir, find, env, gpu=gpu)
    t["solve"] = round(time.perf_counter() - s, 2)
    return {"timing": t, "per_step": _solve.extract_forces_per_physical_step(workdir)}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", required=True, help="case.json (geometry + flow conditions)")
    ap.add_argument("--out", default="/tmp/unsteady_val")
    ap.add_argument("--alphas", default="0,4", help="comma list, evenly spaced")
    ap.add_argument("--dt", type=float, default=1.0e6, help="unsteady physical step size [s]")
    ap.add_argument("--max-pseudo", type=int, default=2000)
    ap.add_argument("--steady-steps", type=int, default=2000, help="steady max_steps (reference)")
    ap.add_argument("--solve", action="store_true")
    ap.add_argument("--gpu", type=int, default=0)
    ap.add_argument("--compute-root", default=None)
    args = ap.parse_args()

    alphas = [float(a) for a in args.alphas.split(",")]
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    env, find = make_env(args.compute_root)
    mesher_bin = Path(__file__).resolve().parent / "bin" / "mesh2d"

    cfg = CaseConfig.load(args.case)
    walls, bnd = boundary_lists(cfg)
    print(f"[val] alphas={alphas} dt={args.dt:g}s max_pseudo={args.max_pseudo} "
          f"steady_steps={args.steady_steps} solve={args.solve}", flush=True)

    print("[val] building mesh once …", flush=True)
    build_mesh(cfg, out / "mesh", mesher_bin, env, find)
    mesh_src = out / "mesh" / MESH

    # --- unsteady α sweep ---
    print("[val] unsteady sweep: preprocess" + ("+solve" if args.solve else ""), flush=True)
    sweep_builder = partial(_case.build_alpha_sweep_simulation_json,
                            alphas_deg=alphas, step_size_s=args.dt,
                            max_pseudo_steps=args.max_pseudo)
    cfg.flow.alpha_deg = alphas[0]
    uns = prep_and_solve(out / "unsteady", mesh_src, cfg, walls, bnd, find, env,
                         sim_builder=sweep_builder, solve=args.solve, gpu=args.gpu)

    # --- steady reference per α ---
    steady = {}
    for a in alphas:
        print(f"[val] steady α={a}: preprocess" + ("+solve" if args.solve else ""), flush=True)
        c = CaseConfig.load(args.case)
        c.flow.alpha_deg = a
        c.solver.max_steps = args.steady_steps
        r = prep_and_solve(out / f"steady_a{a:g}", mesh_src, c, *boundary_lists(c), find, env,
                           sim_builder=_case.build_simulation_json, solve=args.solve, gpu=args.gpu)
        steady[a] = r

    result = {"alphas": alphas, "dt": args.dt, "max_pseudo": args.max_pseudo,
              "unsteady": uns, "steady": {str(a): steady[a] for a in alphas}}
    (out / "validation.json").write_text(json.dumps(result, indent=2, default=str))

    if args.solve:
        print("\n=== unsteady-vs-steady (same mesh) ===", flush=True)
        print(f"{'α':>6} {'CL_uns':>10} {'CL_steady':>11} {'ΔCL':>9} "
              f"{'CD_uns':>10} {'CD_steady':>11} {'ΔCD':>9}", flush=True)
        per = {p["physical_step"]: p for p in uns.get("per_step", [])}
        for i, a in enumerate(alphas):
            pu = per.get(i, {})
            ps = steady[a].get("per_step", [{}])[-1]
            clu, cls = pu.get("CL"), ps.get("CL")
            cdu, cds = pu.get("CD"), ps.get("CD")
            dcl = (clu - cls) if (clu is not None and cls is not None) else float("nan")
            dcd = (cdu - cds) if (cdu is not None and cds is not None) else float("nan")
            print(f"{a:6.1f} {clu!s:>10} {cls!s:>11} {dcl:9.4f} "
                  f"{cdu!s:>10} {cds!s:>11} {dcd:9.4f}", flush=True)
        steady_solve_s = {a: r["timing"].get("solve") for a, r in steady.items()}
        print(f"\n[val] timings: unsteady {uns['timing']} | steady_solve_s {steady_solve_s}", flush=True)
    else:
        print("[val] build OK (no solve). Flow360 JSON accepted for all cases.", flush=True)
    print(f"[val] wrote {out/'validation.json'}", flush=True)


if __name__ == "__main__":
    main()
