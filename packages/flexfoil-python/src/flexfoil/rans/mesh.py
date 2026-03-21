"""Mesh generation for pseudo-2D airfoil RANS via Flow360.

Two approaches:
1. **CSM geometry** (preferred): Generate a CSM file from airfoil coordinates,
   upload to Flow360, and let its automated mesher handle surface + volume meshing.
   This produces high-quality BL meshes with proper wake refinement.

2. **Direct UGRID** (fallback): Generate a structured C-grid mesh in pure numpy
   and write it as UGRID binary. No external dependencies but mesh quality is limited.

No external meshing dependencies — pure numpy + string generation.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def estimate_first_cell_height(
    Re: float, chord: float = 1.0, y_plus: float = 1.0
) -> float:
    """Estimate first cell height for a target y+ using flat-plate correlation."""
    cf = 0.058 * Re ** (-0.2)
    u_tau_norm = np.sqrt(cf / 2.0)
    y1 = y_plus * chord / (Re * u_tau_norm)
    return float(y1)


# ---------------------------------------------------------------------------
# CSM geometry generation (for Flow360 automated meshing)
# ---------------------------------------------------------------------------

def generate_csm(
    coords: list[tuple[float, float]],
    *,
    span: float = 0.01,
    chord: float = 1.0,
) -> str:
    """Generate an OpenCSM (.csm) file for a quasi-2D airfoil extrusion.

    The CSM file defines the airfoil as a spline sketch in the x-z plane,
    then extrudes it along the y-axis by `span`. This produces a thin slab
    geometry suitable for pseudo-2D RANS with symmetry BCs.

    Parameters
    ----------
    coords : list of (x, y) tuples
        Airfoil coordinates in Selig ordering (upper TE → LE → lower TE).
        The y-coordinate becomes z in the CSM file (x-z plane).
    span : float
        Extrusion depth along y-axis.
    chord : float
        Reference chord for normalization.

    Returns
    -------
    str
        CSM file contents.
    """
    pts = np.array(coords, dtype=np.float64)
    n = len(pts)

    if n < 20:
        raise ValueError(f"Need at least 20 points, got {n}")

    # Find LE (leftmost point) to split upper/lower
    le_idx = np.argmin(pts[:, 0])

    # Upper surface: TE → LE (indices 0..le_idx)
    upper = pts[:le_idx + 1]
    # Lower surface: LE → TE (indices le_idx..end)
    lower = pts[le_idx:]

    # TE points
    te_upper = upper[0]
    te_lower = lower[-1]

    lines = []
    lines.append(f"# FlexFoil quasi-2D airfoil geometry")
    lines.append(f"# {n} surface points, span={span}")
    lines.append("")

    # Upper surface sketch: TE → LE (x decreasing)
    # CSM coordinates: x=airfoil_x, y=0, z=airfoil_y
    lines.append(f"skbeg {te_upper[0]:.10f} 0 {te_upper[1]:.10f} 0")
    for i in range(1, len(upper)):
        x, z = upper[i]
        lines.append(f"spline {x:.10f} 0 {z:.10f}")
    lines.append("skend 0")

    # Lower surface sketch: LE → TE (x increasing)
    lines.append(f"skbeg {lower[0][0]:.10f} 0 {lower[0][1]:.10f} 0")
    for i in range(1, len(lower)):
        x, z = lower[i]
        lines.append(f"spline {x:.10f} 0 {z:.10f}")
    lines.append("skend 0")

    # TE closure: connect lower TE to upper TE
    lines.append(f"skbeg {te_lower[0]:.10f} 0 {te_lower[1]:.10f} 0")
    lines.append(f"linseg {te_upper[0]:.10f} 0 {te_upper[1]:.10f}")
    lines.append("skend 0")

    # Combine the three sketches into a single closed cross-section
    lines.append("combine")

    # Extrude along y-axis
    lines.append(f"extrude 0 {span:.10f} 0")

    # Tag boundaries
    lines.append('select face')
    lines.append('attribute groupName $airfoil')

    lines.append("")
    lines.append("end")

    return "\n".join(lines) + "\n"


def write_csm(
    path: str | Path,
    coords: list[tuple[float, float]],
    *,
    span: float = 0.01,
) -> Path:
    """Write a CSM geometry file for an airfoil."""
    path = Path(path)
    csm_content = generate_csm(coords, span=span)
    path.write_text(csm_content)
    return path


# ---------------------------------------------------------------------------
# Surface mesh configuration (for Flow360 automated meshing)
# ---------------------------------------------------------------------------

def build_surface_mesh_config(
    *,
    max_edge_length: float = 0.05,
    curvature_resolution: float = 15.0,
    growth_rate: float = 1.2,
) -> dict:
    """Build surface mesh configuration for Flow360's automated mesher."""
    return {
        "maxEdgeLength": max_edge_length,
        "curvatureResolutionAngle": curvature_resolution,
        "growthRate": growth_rate,
    }


def build_volume_mesh_config(
    *,
    first_layer_thickness: float | None = None,
    Re: float = 1e6,
    growth_rate: float = 1.15,
    n_layers: int = 40,
    farfield_type: str = "quasi-3d",
) -> dict:
    """Build volume mesh configuration for Flow360's automated mesher."""
    if first_layer_thickness is None:
        first_layer_thickness = estimate_first_cell_height(Re)

    return {
        "firstLayerThickness": first_layer_thickness,
        "growthRate": growth_rate,
        "volume": {
            "firstLayerThickness": first_layer_thickness,
            "growthRate": growth_rate,
        },
    }


# ---------------------------------------------------------------------------
# Direct C-grid mesh generation (fallback when automated meshing unavailable)
# ---------------------------------------------------------------------------

def _compute_normals(curve: np.ndarray) -> np.ndarray:
    """Compute unit outward normals along an open curve."""
    tangents = np.zeros_like(curve)
    tangents[1:-1] = curve[2:] - curve[:-2]
    tangents[0] = curve[1] - curve[0]
    tangents[-1] = curve[-1] - curve[-2]

    lengths = np.linalg.norm(tangents, axis=1, keepdims=True)
    lengths = np.maximum(lengths, 1e-14)
    tangents = tangents / lengths

    normals = np.column_stack([tangents[:, 1], -tangents[:, 0]])

    centroid = curve.mean(axis=0)
    outward = curve - centroid
    dots = np.sum(normals * outward, axis=1)
    if np.sum(dots < 0) > np.sum(dots > 0):
        normals = -normals

    return normals


def generate_airfoil_mesh(
    coords: list[tuple[float, float]],
    *,
    n_normal: int = 80,
    n_wake: int = 40,
    first_cell_height: float | None = None,
    Re: float = 1e6,
    growth_rate: float = 1.1,
    farfield_radius: float = 50.0,
    wake_length: float = 5.0,
    chord: float = 1.0,
    y_plus: float = 1.0,
) -> dict:
    """Generate a structured C-grid mesh around a 2D airfoil.

    The C-grid splits the airfoil at the TE and extends a wake cut downstream.
    """
    surface = np.array(coords, dtype=np.float64)
    n_pts = len(surface)
    if n_pts < 20:
        raise ValueError(f"Need at least 20 surface points, got {n_pts}")

    le_idx = np.argmin(surface[:, 0])
    upper = surface[:le_idx + 1]
    lower = surface[le_idx:]

    te_upper = upper[0].copy()
    te_lower = lower[-1].copy()

    if first_cell_height is None:
        first_cell_height = estimate_first_cell_height(Re, chord, y_plus)

    # Wake points
    wake_dx = np.zeros(n_wake)
    wake_first = chord * 0.01
    for i in range(n_wake):
        wake_dx[i] = wake_first * (1.2 ** i)
    wake_x_offsets = np.cumsum(wake_dx)
    if wake_x_offsets[-1] > 0:
        wake_x_offsets *= (wake_length * chord) / wake_x_offsets[-1]

    wake_upper_pts = np.column_stack([te_upper[0] + wake_x_offsets, np.full(n_wake, te_upper[1])])
    wake_lower_pts = np.column_stack([te_lower[0] + wake_x_offsets, np.full(n_wake, te_lower[1])])

    lower_reversed = lower[::-1]
    wake_lower_reversed = wake_lower_pts[::-1]

    c_boundary = np.vstack([
        wake_lower_reversed,
        lower_reversed[1:],
        upper[1:],
        wake_upper_pts,
    ])
    n_c = len(c_boundary)

    normals = _compute_normals(c_boundary)
    for i in range(n_wake):
        normals[i] = [0.0, -1.0]
    for i in range(n_c - n_wake, n_c):
        normals[i] = [0.0, 1.0]

    n_layers = n_normal + 1
    heights = np.zeros(n_layers)
    for i in range(1, n_layers):
        heights[i] = heights[i - 1] + first_cell_height * growth_rate ** (i - 1)

    max_h = heights[-1]
    target = farfield_radius * chord
    if max_h < target:
        scale = target / max_h
        t = np.linspace(0, 1, n_layers)
        blend = t ** 1.5
        heights = heights * (1.0 - blend) + heights * scale * blend
    elif max_h > target:
        heights *= target / max_h

    nodes_2d = np.zeros((n_c * n_layers, 2))
    nodes_2d[:n_c] = c_boundary
    for j in range(1, n_layers):
        nodes_2d[j * n_c: (j + 1) * n_c] = c_boundary + heights[j] * normals

    return {
        "nodes_2d": nodes_2d,
        "n_c": n_c,
        "n_layers": n_layers,
        "n_wake": n_wake,
        "c_boundary": c_boundary,
        "first_cell_height": first_cell_height,
    }


def extrude_to_3d(mesh_2d: dict, *, span: float = 0.01) -> dict:
    """Extrude a 2D C-grid one cell deep in y. Airfoil in x-z plane."""
    nodes_2d = mesh_2d["nodes_2d"]
    n_c = mesh_2d["n_c"]
    n_layers = mesh_2d["n_layers"]
    n_nodes_2d = len(nodes_2d)
    n_wake = mesh_2d["n_wake"]

    nodes_y0 = np.column_stack([nodes_2d[:, 0], np.zeros(n_nodes_2d), nodes_2d[:, 1]])
    nodes_y1 = np.column_stack([nodes_2d[:, 0], np.full(n_nodes_2d, span), nodes_2d[:, 1]])
    nodes = np.vstack([nodes_y0, nodes_y1])

    n_cells_c = n_c - 1
    n_cells_normal = n_layers - 1

    hexes = np.zeros((n_cells_c * n_cells_normal, 8), dtype=np.int32)
    idx = 0
    for j in range(n_cells_normal):
        for i in range(n_cells_c):
            n0 = j * n_c + i
            n1 = j * n_c + (i + 1)
            n2 = (j + 1) * n_c + (i + 1)
            n3 = (j + 1) * n_c + i
            n4, n5, n6, n7 = n0 + n_nodes_2d, n1 + n_nodes_2d, n2 + n_nodes_2d, n3 + n_nodes_2d

            p0, p1, p3, p4 = nodes_y0[n0], nodes_y0[n1], nodes_y0[n3], nodes_y1[n0]
            vol = np.dot(p1 - p0, np.cross(p3 - p0, p4 - p0))
            if vol > 0:
                hexes[idx] = [n0, n1, n2, n3, n4, n5, n6, n7]
            else:
                hexes[idx] = [n0, n3, n2, n1, n4, n7, n6, n5]
            idx += 1

    hexes = hexes[:idx] + 1

    # Boundary faces
    airfoil_start = n_wake
    airfoil_end = n_c - n_wake
    wall_quads = []
    for i in range(airfoil_start, airfoil_end - 1):
        wall_quads.append([i, i + n_nodes_2d, (i + 1) + n_nodes_2d, i + 1])

    farfield_quads = []
    j = n_cells_normal
    for i in range(n_cells_c):
        n0, n1 = j * n_c + i, j * n_c + (i + 1)
        farfield_quads.append([n0, n1, n1 + n_nodes_2d, n0 + n_nodes_2d])

    wake_end1 = []
    wake_end2 = []
    for j in range(n_cells_normal):
        n0, n3 = j * n_c, (j + 1) * n_c
        wake_end1.append([n0, n0 + n_nodes_2d, n3 + n_nodes_2d, n3])
        n0r, n3r = j * n_c + (n_c - 1), (j + 1) * n_c + (n_c - 1)
        wake_end2.append([n0r, n3r, n3r + n_nodes_2d, n0r + n_nodes_2d])

    sym_y0, sym_y1 = [], []
    for j in range(n_cells_normal):
        for i in range(n_cells_c):
            n0 = j * n_c + i
            n1, n2, n3 = n0 + 1, (j + 1) * n_c + (i + 1), (j + 1) * n_c + i
            sym_y0.append([n0, n3, n2, n1])
            sym_y1.append([n0 + n_nodes_2d, n1 + n_nodes_2d, n2 + n_nodes_2d, n3 + n_nodes_2d])

    boundary_quads = {
        "wall": np.array(wall_quads, dtype=np.int32) + 1,
        "farfield": np.array(farfield_quads, dtype=np.int32) + 1,
        "symmetry_y0": np.array(sym_y0, dtype=np.int32) + 1,
        "symmetry_y1": np.array(sym_y1, dtype=np.int32) + 1,
        "wake": np.vstack([np.array(wake_end1, dtype=np.int32), np.array(wake_end2, dtype=np.int32)]) + 1,
    }

    boundary_ids = {"wall": 1, "farfield": 2, "symmetry_y0": 3, "symmetry_y1": 4, "wake": 5}

    return {
        "nodes": nodes, "hexes": hexes,
        "boundary_quads": boundary_quads, "boundary_ids": boundary_ids,
        "n_nodes": len(nodes), "n_hexes": len(hexes),
    }


# ---------------------------------------------------------------------------
# UGRID writer
# ---------------------------------------------------------------------------

def write_ugrid(path: str | Path, mesh_3d: dict) -> None:
    """Write a 3D hex mesh in AFLR3/UGRID big-endian binary format (.b8.ugrid)."""
    path = Path(path)
    nodes, hexes = mesh_3d["nodes"], mesh_3d["hexes"]
    boundary_quads, boundary_ids = mesh_3d["boundary_quads"], mesh_3d["boundary_ids"]

    all_bquads, all_tags = [], []
    for name, quads in boundary_quads.items():
        all_bquads.append(quads)
        all_tags.extend([boundary_ids[name]] * len(quads))
    all_bquads = np.vstack(all_bquads) if all_bquads else np.zeros((0, 4), dtype=np.int32)

    with open(path, "wb") as f:
        f.write(struct.pack(">7i", len(nodes), 0, len(all_bquads), 0, 0, 0, len(hexes)))
        for node in nodes:
            f.write(struct.pack(">3d", *node))
        for quad in all_bquads:
            f.write(struct.pack(">4i", *quad))
        for tag in all_tags:
            f.write(struct.pack(">i", tag))
        for hex_elem in hexes:
            f.write(struct.pack(">8i", *hex_elem))


def write_mapbc(path: str | Path, mesh_3d: dict) -> None:
    """Write .mapbc boundary condition mapping file."""
    path = Path(path)
    ids = mesh_3d["boundary_ids"]
    bc_map = {"wall": 4000, "farfield": 3000, "symmetry_y0": 5000, "symmetry_y1": 5000, "wake": 3000}
    lines = [str(len(ids))]
    for name, tag in sorted(ids.items(), key=lambda x: x[1]):
        lines.append(f"{tag} {bc_map.get(name, 0)} {name}")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Top-level pipelines
# ---------------------------------------------------------------------------

def generate_and_write_csm(
    coords: list[tuple[float, float]],
    output_dir: str | Path,
    *,
    span: float = 0.01,
    mesh_name: str = "airfoil",
) -> Path:
    """Generate a CSM geometry file from airfoil coords.

    Returns path to the CSM file.
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    csm_path = output_dir / f"{mesh_name}.csm"
    write_csm(csm_path, coords, span=span)
    return csm_path


def generate_and_write_mesh(
    coords: list[tuple[float, float]],
    output_dir: str | Path,
    *,
    Re: float = 1e6,
    n_normal: int = 80,
    n_wake: int = 40,
    growth_rate: float = 1.1,
    farfield_radius: float = 50.0,
    wake_length: float = 5.0,
    span: float = 0.01,
    mesh_name: str = "airfoil",
) -> tuple[Path, Path]:
    """Full pipeline: airfoil coords → UGRID mesh files (C-grid fallback)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh_2d = generate_airfoil_mesh(
        coords, n_normal=n_normal, n_wake=n_wake, Re=Re,
        growth_rate=growth_rate, farfield_radius=farfield_radius,
        wake_length=wake_length,
    )
    mesh_3d = extrude_to_3d(mesh_2d, span=span)

    ugrid_path = output_dir / f"{mesh_name}.b8.ugrid"
    mapbc_path = output_dir / f"{mesh_name}.mapbc"
    write_ugrid(ugrid_path, mesh_3d)
    write_mapbc(mapbc_path, mesh_3d)

    return ugrid_path, mapbc_path
