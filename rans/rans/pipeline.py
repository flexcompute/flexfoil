"""End-to-end orchestrator: contours JSON -> mesh -> Flow360 case -> (solve).

    from rans.pipeline import run
    run("case.json", "out/")            # build mesh + case (in-session OK)
    run("case.json", "out/", solve=True)  # also solve (needs a non-sandboxed shell)

Everything up to and including ``preprocess`` runs anywhere; ``solve`` needs a
normal interactive shell (see rans.solve).
"""
from __future__ import annotations

import json
import time
from pathlib import Path

from . import case as _case
from . import contours as _contours
from . import mesh as _mesh
from .config import CaseConfig
from .env import make_env

MESH_NAME = "mesh.cgns"
# Bump when case.py's simulation-JSON builders change in a way the cache key wouldn't
# otherwise capture (e.g. CFL controller, solver tolerances), to invalidate stale caches.
_CASE_BUILD_VERSION = 2


def run(config_path: str | Path, outdir: str | Path, *, solve: bool = False,
        compute_root: str | Path | None = None, gpu: int = 0, fast: bool = False,
        cache_root: str | Path = "/tmp/rans_cache",
        flow_field: bool = True,
        alpha_sweep: list[float] | None = None,
        sweep_step_size_s: float = 1.0e6, sweep_max_pseudo: int = 3000) -> dict:
    """``fast`` enables the fast-iteration path: cache the mesh-independent SDK case
    JSONs (skip the ~10 s flow360 imports on a cache hit, keyed by flow conditions +
    boundary structure + steps) and skip auto-visualization. Pair it with a coarse
    mesh + fewer steps in the config for a ~15 s turnaround.

    ``alpha_sweep`` (a list of α in degrees, evenly spaced) builds the unsteady-as-steady
    sweep case instead of a steady one: one physical step per α, huge Δt ⇒ steady per
    step, a UDD driving the freestream angle (warm-started from the previous α). The
    per-α solve + extraction is progressive, driven by the caller (rans_server)."""
    cfg = CaseConfig.load(config_path)
    if alpha_sweep:
        cfg.flow.alpha_deg = alpha_sweep[0]
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    env, find = make_env(compute_root)
    pkg_root = Path(__file__).resolve().parent.parent
    mesher_bin = pkg_root / "bin" / "mesh2d"

    summary: dict = {"outdir": str(out)}
    timing: dict = {}
    clk = time.perf_counter
    t_all = clk()

    # 1. contours -> mesher input
    s = clk()
    summary["contours"] = _contours.write_contours(cfg, out / "contours.txt")
    timing["contours_txt"] = round(clk() - s, 3)

    # 2. 2D anisotropic mesh (cwd=out so mesher side artifacts land there)
    s = clk()
    _mesh.run_mesher(out / "contours.txt", out / "mesh2d.vtk", mesher_bin, env, cwd=out)
    timing["mesh2d"] = round(clk() - s, 3)

    # 3. extrude to quasi-2D volume + name patches
    s = clk()
    summary["mesh"] = _mesh.write_volume_msh(out / "mesh2d.vtk", cfg, out / "mesh.msh")
    timing["extrude"] = round(clk() - s, 3)

    # 4. gmsh -> cgns
    s = clk()
    _mesh.gmsh_to_cgns(out / "mesh.msh", out / MESH_NAME, find("flow360gmshtocgns"), env)
    timing["gmsh_to_cgns"] = round(clk() - s, 3)

    # 5. boundary names follow the CGNS 'fluid/<patch>' convention
    elem_names = [e.name for e in cfg.elements]
    wall_names = [f"fluid/{e.name}" for e in cfg.elements if e.is_wall]
    boundary_names = (["fluid/farfield"] + [f"fluid/{n}" for n in elem_names]
                      + ["fluid/symmetry1", "fluid/symmetry2"])

    # SDK case-JSON cache key: everything that determines the (mesh-independent)
    # simulation.json / Flow360_processed.json — flow conditions, span, steps, boundaries.
    sdk_cache_dir = None
    if fast:
        import hashlib
        key_src = json.dumps({
            "mach": cfg.flow.mach, "re": cfg.flow.reynolds, "alpha": cfg.flow.alpha_deg,
            "T": cfg.flow.temperature, "span": cfg.mesh.span, "steps": cfg.solver.max_steps,
            "bnd": boundary_names,
            # the sweep produces a different (unsteady + UDD) case JSON, so key on it too
            "sweep": [alpha_sweep, sweep_step_size_s, sweep_max_pseudo] if alpha_sweep else None,
            # bump when the case-JSON builders change (e.g. CFL controller) so stale
            # cached SDK JSONs aren't reused across code changes
            "build_v": _CASE_BUILD_VERSION,
        }, sort_keys=True)
        sdk_cache_dir = Path(cache_root) / hashlib.md5(key_src.encode()).hexdigest()[:16]

    # 6. preprocessing chain -> Flow360.json (builds the SDK case JSONs; sub-timings recorded)
    sim_builder = _case.build_simulation_json
    if alpha_sweep:
        from functools import partial
        sim_builder = partial(_case.build_alpha_sweep_simulation_json, alphas_deg=alpha_sweep,
                              step_size_s=sweep_step_size_s, max_pseudo_steps=sweep_max_pseudo)
    pre: dict = {}
    s = clk()
    summary["flow360_json"] = _case.preprocess(
        out, MESH_NAME, find, env, cfg=cfg, wall_names=wall_names, boundary_names=boundary_names,
        timings=pre, sdk_cache_dir=sdk_cache_dir, sim_builder=sim_builder)
    timing["preprocess"] = {"total": round(clk() - s, 3), **pre}

    # 7. solve (optional; non-sandboxed shell only). A sweep's per-α solve + extraction
    #    is driven by the caller (the progressive server path); here ``solve`` is the
    #    single steady solve that also writes the LIC flow field.
    if solve:
        from . import solve as _solve
        from . import flowfield as _ff
        s = clk()
        _solve.run_solver(out, find, env, gpu=gpu)
        timing["solve"] = round(clk() - s, 3)
        s = clk()
        summary["forces"] = _solve.extract_forces(out)
        timing["extract"] = round(clk() - s, 3)
        if flow_field:
            s = clk()
            (out / "flow_field.json").write_text(json.dumps(_ff.extract_flow_mesh(out)))
            summary["flow_field"] = str(out / "flow_field.json")
            timing["flow_field"] = round(clk() - s, 3)

    timing["total"] = round(clk() - t_all, 3)
    summary["timing_s"] = timing
    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary
