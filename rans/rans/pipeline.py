"""End-to-end orchestrator: contours JSON -> mesh -> Flow360 case -> (solve).

    from rans.pipeline import run
    run("case.json", "out/")            # build mesh + case (in-session OK)
    run("case.json", "out/", solve=True)  # also solve (needs a non-sandboxed shell)

Everything up to and including ``preprocess`` runs anywhere; ``solve`` needs a
normal interactive shell (see rans.solve).
"""
from __future__ import annotations

import json
from pathlib import Path

from . import case as _case
from . import contours as _contours
from . import mesh as _mesh
from .config import CaseConfig
from .env import make_env

MESH_NAME = "mesh.cgns"


def run(config_path: str | Path, outdir: str | Path, *, solve: bool = False,
        compute_root: str | Path | None = None, gpu: int = 0) -> dict:
    cfg = CaseConfig.load(config_path)
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    env, find = make_env(compute_root)
    pkg_root = Path(__file__).resolve().parent.parent
    mesher_bin = pkg_root / "bin" / "mesh2d"

    summary: dict = {"outdir": str(out)}

    # 1. contours -> mesher input
    summary["contours"] = _contours.write_contours(cfg, out / "contours.txt")

    # 2. 2D anisotropic mesh (cwd=out so mesher side artifacts land there)
    _mesh.run_mesher(out / "contours.txt", out / "mesh2d.vtk", mesher_bin, env, cwd=out)

    # 3. extrude to quasi-2D volume + name patches
    summary["mesh"] = _mesh.write_volume_msh(out / "mesh2d.vtk", cfg, out / "mesh.msh")

    # 4. gmsh -> cgns
    _mesh.gmsh_to_cgns(out / "mesh.msh", out / MESH_NAME, find("flow360gmshtocgns"), env)

    # 5. simulation.json (boundary names follow the CGNS 'fluid/<patch>' convention)
    elem_names = [e.name for e in cfg.elements]
    wall_names = [f"fluid/{e.name}" for e in cfg.elements if e.is_wall]
    boundary_names = (["fluid/farfield"] + [f"fluid/{n}" for n in elem_names]
                      + ["fluid/symmetry1", "fluid/symmetry2"])
    _case.build_simulation_json(cfg, wall_names, boundary_names, out / "simulation.json")

    # 6. preprocessing chain -> Flow360.json
    summary["flow360_json"] = _case.preprocess(out, MESH_NAME, find, env)

    # 7. solve (optional; non-sandboxed shell only)
    if solve:
        from . import solve as _solve
        _solve.run_solver(out, find, env, gpu=gpu)
        summary["forces"] = _solve.extract_forces(out)

    (out / "summary.json").write_text(json.dumps(summary, indent=2))
    return summary
