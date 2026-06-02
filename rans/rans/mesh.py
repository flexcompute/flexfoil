"""Meshing: contours.txt -> 2D anisotropic mesh (C++ mesher) -> quasi-2D volume.

The 2D in-plane mesh comes from the in-house CavityBasedMesher driven by the
a-priori Spalding metric (the ``mesh2d`` binary). It is then extruded one or more
cells in the span and oriented so streamwise=x, lift=z, span=y (a rigid +90 deg
rotation about x — preserves prism orientation). Boundary edges are classified
against the original element/farfield curves to name the wall patches.
"""
from __future__ import annotations

import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np

from .config import CaseConfig
from .contours import _farfield_loop


def run_mesher(contours_txt: str | Path, vtk_out: str | Path,
               mesher_bin: str | Path, env: dict, cwd: str | Path | None = None) -> None:
    """Run the C++ mesher: contours.txt -> 2D triangle mesh (legacy VTK).
    ``cwd`` keeps the mesher's side artifacts (e.g. *.meshb) out of the package dir."""
    subprocess.run([str(mesher_bin), str(Path(contours_txt).resolve()),
                    str(Path(vtk_out).resolve())],
                   check=True, env=env, capture_output=True, text=True,
                   cwd=str(cwd) if cwd else None)


def _read_vtk(path: str | Path):
    toks = Path(path).read_text().split("\n")
    P, T, i = [], [], 0
    while i < len(toks):
        w = toks[i].split()
        if w and w[0] == "POINTS":
            n = int(w[1]); vals = []; i += 1
            while len(vals) < 3 * n:
                vals += toks[i].split(); i += 1
            P = np.array(vals, float).reshape(-1, 3)[:, :2]; continue
        if w and w[0] == "CELLS":
            m = int(w[1]); i += 1
            for _ in range(m):
                c = toks[i].split(); T.append([int(c[1]), int(c[2]), int(c[3])]); i += 1
            continue
        i += 1
    return np.array(P), np.array(T)


def _named_curves(cfg: CaseConfig):
    """[(name, Nx2 points)] for the farfield + each element, for edge classification."""
    curves = [("farfield", np.array(_farfield_loop(cfg)))]
    for el in cfg.elements:
        curves.append((el.name, np.array(el.contour)))
    return curves


def write_volume_msh(vtk_2d: str | Path, cfg: CaseConfig, msh_out: str | Path) -> dict:
    """Extrude the 2D mesh to a quasi-2D Gmsh v2.2 volume mesh with named patches."""
    P, T = _read_vtk(vtk_2d)
    N = len(P)
    span, nspan = cfg.mesh.span, cfg.mesh.nspan

    curves = _named_curves(cfg)
    cnames = [c[0] for c in curves]
    ckd = [c[1] for c in curves]

    def nearest_curve(mx, my):
        best, bd = -1, 1e30
        for ci, pts in enumerate(ckd):
            d = np.min((pts[:, 0] - mx) ** 2 + (pts[:, 1] - my) ** 2)
            if d < bd:
                bd, best = d, ci
        return best

    ec = defaultdict(int)
    for a, b, c in T:
        for e in [(a, b), (b, c), (c, a)]:
            ec[tuple(sorted(e))] += 1
    bedges = [e for e, n in ec.items() if n == 1]
    edge_tag = {}
    for e in bedges:
        mx, my = 0.5 * (P[e[0]] + P[e[1]])
        edge_tag[e] = nearest_curve(mx, my)

    # patch tags: fluid (3D) = 1; 2D patches = walls/farfield + the two span faces
    wall_names = [c for c in cnames]                      # farfield + elements
    patch_names = wall_names + ["symmetry1", "symmetry2"]
    surf_tag = {n: i + 2 for i, n in enumerate(patch_names)}
    phys = [(2, surf_tag[n], n) for n in patch_names] + [(3, 1, "fluid")]

    NL = nspan + 1

    def nid(L, k):
        return L * N + k + 1

    elems, eid = [], 1

    def emit(s):
        nonlocal eid
        elems.append(f"{eid} {s}"); eid += 1

    with open(msh_out, "w") as f:
        f.write("$MeshFormat\n2.2 0 8\n$EndMeshFormat\n")
        f.write(f"$PhysicalNames\n{len(phys)}\n")
        for dim, tag, name in phys:
            f.write(f'{dim} {tag} "{name}"\n')
        f.write("$EndPhysicalNames\n")
        f.write(f"$Nodes\n{NL * N}\n")
        for L in range(NL):
            y = -span * L / nspan
            for k in range(N):
                f.write(f"{nid(L, k)} {P[k,0]:.16g} {y:.16g} {P[k,1]:.16g}\n")
        f.write("$EndNodes\n")
        # 2D first: span end planes (tris) then side walls (quads)
        for a, b, c in T:
            emit(f"2 2 {surf_tag['symmetry1']} {surf_tag['symmetry1']} {nid(0,a)} {nid(0,b)} {nid(0,c)}")
        for a, b, c in T:
            emit(f"2 2 {surf_tag['symmetry2']} {surf_tag['symmetry2']} {nid(nspan,a)} {nid(nspan,b)} {nid(nspan,c)}")
        for e in bedges:
            a, b = e; tg = surf_tag[cnames[edge_tag[e]]]
            for L in range(nspan):
                emit(f"3 2 {tg} {tg} {nid(L,a)} {nid(L,b)} {nid(L+1,b)} {nid(L+1,a)}")
        # 3D prisms
        for L in range(nspan):
            for a, b, c in T:
                emit(f"6 2 1 1 {nid(L,a)} {nid(L,b)} {nid(L,c)} {nid(L+1,a)} {nid(L+1,b)} {nid(L+1,c)}")
        f.write(f"$Elements\n{len(elems)}\n")
        f.write("\n".join(elems) + "\n$EndElements\n")

    return {
        "n_points_2d": int(N),
        "n_tris_2d": int(len(T)),
        "n_prisms": int(nspan * len(T)),
        "boundary_edges": {cnames[k]: v for k, v in Counter(edge_tag.values()).items()},
        "wall_patches": [e.name for e in cfg.elements if e.is_wall],
    }


def gmsh_to_cgns(msh: str | Path, cgns: str | Path, gmshtocgns_bin: str | Path,
                 env: dict) -> None:
    subprocess.run([str(gmshtocgns_bin), str(msh), "-o", str(cgns)],
                   check=True, env=env, capture_output=True, text=True)
