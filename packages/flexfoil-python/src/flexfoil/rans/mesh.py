"""Structured C-grid mesh generation for pseudo-2D airfoil RANS.

Generates a single-cell-deep hex mesh from 2D airfoil coordinates using a
C-grid topology. The airfoil contour is split at the trailing edge and a
wake cut extends downstream, giving the mesh lines a clean exit path.

This avoids the degenerate trailing-edge cells and poor outer-field coverage
that plague O-grid approaches for airfoils.

No external meshing dependencies — pure numpy.
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


def _compute_normals(curve: np.ndarray) -> np.ndarray:
    """Compute unit outward normals along an open curve.

    Uses central differences for interior points, one-sided at endpoints.
    Outward direction is determined by the curve's winding sense.
    """
    n = len(curve)
    tangents = np.zeros_like(curve)
    tangents[1:-1] = curve[2:] - curve[:-2]
    tangents[0] = curve[1] - curve[0]
    tangents[-1] = curve[-1] - curve[-2]

    lengths = np.linalg.norm(tangents, axis=1, keepdims=True)
    lengths = np.maximum(lengths, 1e-14)
    tangents = tangents / lengths

    # Rotate 90° to get normal: (tx, ty) → (ty, -tx)
    normals = np.column_stack([tangents[:, 1], -tangents[:, 0]])

    # Ensure normals point outward (away from contour centroid)
    centroid = curve.mean(axis=0)
    outward = curve - centroid
    dots = np.sum(normals * outward, axis=1)
    if np.sum(dots < 0) > np.sum(dots > 0):
        normals = -normals

    return normals


# ---------------------------------------------------------------------------
# C-grid generation
# ---------------------------------------------------------------------------

def generate_airfoil_mesh(
    coords: list[tuple[float, float]] | list[list[tuple[float, float]]],
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

    The C-grid topology:
    - The inner boundary is: wake_lower → TE_lower → LE → TE_upper → wake_upper
    - Mesh lines march outward from this C-shaped boundary
    - The wake cut extends downstream from the TE
    - No degenerate cells at the trailing edge

    Parameters
    ----------
    coords : list of (x, y) or list of lists for multi-body
        Airfoil coordinates in Selig ordering (upper TE → LE → lower TE).
    n_normal : int
        Number of cells in the wall-normal direction.
    n_wake : int
        Number of streamwise cells in the wake region downstream of TE.
    first_cell_height : float or None
        First cell height (auto-estimated from Re if None).
    Re : float
        Reynolds number.
    growth_rate : float
        Geometric growth rate for boundary-layer cells.
    farfield_radius : float
        Farfield distance in chord lengths.
    wake_length : float
        Wake extent downstream of TE in chord lengths.
    chord : float
        Reference chord length.
    y_plus : float
        Target y+ for first cell.

    Returns
    -------
    dict with mesh data including C-grid topology info.
    """
    # Normalize input
    if (len(coords) > 0
            and isinstance(coords[0], (tuple, list))
            and len(coords[0]) == 2
            and isinstance(coords[0][0], (int, float))):
        surface = np.array(coords, dtype=np.float64)
    else:
        # Multi-body: concatenate
        surface = np.vstack([np.array(b, dtype=np.float64) for b in coords])

    n_pts = len(surface)
    if n_pts < 20:
        raise ValueError(f"Need at least 20 surface points, got {n_pts}")

    # Selig ordering: upper TE → LE → lower TE
    # Find LE (leftmost point)
    le_idx = np.argmin(surface[:, 0])

    # Split into upper and lower surfaces
    # Upper: indices 0..le_idx (TE → LE, x decreasing)
    # Lower: indices le_idx..end (LE → TE, x increasing)
    upper = surface[:le_idx + 1]  # TE_upper → LE
    lower = surface[le_idx:]      # LE → TE_lower

    # TE points
    te_upper = upper[0].copy()
    te_lower = lower[-1].copy()
    te_mid = 0.5 * (te_upper + te_lower)

    if first_cell_height is None:
        first_cell_height = estimate_first_cell_height(Re, chord, y_plus)

    # --- Build the C-shaped inner boundary ---
    # Order: wake_lower (far→TE) + lower_reversed (TE→LE) + upper (LE→TE) + wake_upper (TE→far)

    # Wake points: extend downstream from TE
    # Geometric spacing in the wake, coarser than BL
    wake_dx = np.zeros(n_wake)
    wake_first = chord * 0.01  # first wake cell ~1% chord
    for i in range(n_wake):
        wake_dx[i] = wake_first * (1.2 ** i)
    wake_x_offsets = np.cumsum(wake_dx)
    # Scale to fit wake_length
    if wake_x_offsets[-1] > 0:
        wake_x_offsets *= (wake_length * chord) / wake_x_offsets[-1]

    wake_upper_pts = np.column_stack([
        te_upper[0] + wake_x_offsets,
        np.full(n_wake, te_upper[1]),
    ])
    wake_lower_pts = np.column_stack([
        te_lower[0] + wake_x_offsets,
        np.full(n_wake, te_lower[1]),
    ])

    # Reverse lower surface so it goes TE → LE
    lower_reversed = lower[::-1]  # now TE_lower → LE

    # C-boundary: wake_lower_reversed → lower(TE→LE) → upper(LE→TE) → wake_upper
    wake_lower_reversed = wake_lower_pts[::-1]  # far → TE
    c_boundary = np.vstack([
        wake_lower_reversed,   # far downstream → TE_lower
        lower_reversed[1:],    # TE_lower → LE (skip duplicate TE point)
        upper[1:],             # LE → TE_upper (skip duplicate LE point)
        wake_upper_pts,        # TE_upper → far downstream
    ])

    n_c = len(c_boundary)

    # --- Compute normals along the C-boundary ---
    normals = _compute_normals(c_boundary)

    # For wake points, force normals to be purely vertical (±y)
    # Wake lower normals should point downward (-y), wake upper should point upward (+y)
    n_wake_lower = n_wake  # first n_wake points are wake_lower_reversed
    n_wake_upper = n_wake  # last n_wake points are wake_upper

    # Wake lower: normal should point in -y direction
    for i in range(n_wake_lower):
        normals[i] = [0.0, -1.0]

    # Wake upper: normal should point in +y direction
    for i in range(n_c - n_wake_upper, n_c):
        normals[i] = [0.0, 1.0]

    # --- Generate radial layers ---
    n_layers = n_normal + 1
    heights = np.zeros(n_layers)
    for i in range(1, n_layers):
        heights[i] = heights[i - 1] + first_cell_height * growth_rate ** (i - 1)

    # Scale heights so the outermost layer reaches farfield_radius
    max_height = heights[-1]
    target = farfield_radius * chord
    if max_height < target:
        # Use a blending approach: geometric near wall, stretched outer
        scale = target / max_height
        # Quadratic blend: inner layers keep spacing, outer layers stretch
        t = np.linspace(0, 1, n_layers)
        blend = t ** 1.5  # gentle blend
        heights = heights * (1.0 - blend) + heights * scale * blend
    elif max_height > target:
        heights *= target / max_height

    # --- Build 2D mesh nodes ---
    nodes_2d = np.zeros((n_c * n_layers, 2))
    nodes_2d[:n_c] = c_boundary  # layer 0 = C-boundary

    for j in range(1, n_layers):
        h = heights[j]
        layer = c_boundary + h * normals
        nodes_2d[j * n_c: (j + 1) * n_c] = layer

    return {
        "nodes_2d": nodes_2d,
        "n_c": n_c,  # points along the C-boundary
        "n_layers": n_layers,
        "n_wake": n_wake,
        "c_boundary": c_boundary,
        "first_cell_height": first_cell_height,
        "le_idx_in_c": n_wake + len(lower_reversed) - 1,  # LE position in C-boundary
    }


# ---------------------------------------------------------------------------
# 3D extrusion
# ---------------------------------------------------------------------------

def extrude_to_3d(mesh_2d: dict, *, span: float = 0.01) -> dict:
    """Extrude a 2D C-grid mesh one cell deep in the y-direction.

    The airfoil lies in the x-z plane (chord along x, thickness along z).
    Span (extrusion) is along y. This matches Flow360's alphaAngle convention.

    The C-grid has an open topology: the two ends of the C (wake_lower_far and
    wake_upper_far) are NOT connected. This means the mesh has n_c-1 cells
    along the C direction per layer.
    """
    nodes_2d = mesh_2d["nodes_2d"]
    n_c = mesh_2d["n_c"]
    n_layers = mesh_2d["n_layers"]
    n_nodes_2d = len(nodes_2d)

    # 3D nodes: (x_2d, y_span, z_2d_y)
    nodes_y0 = np.column_stack([nodes_2d[:, 0], np.zeros(n_nodes_2d), nodes_2d[:, 1]])
    nodes_y1 = np.column_stack([nodes_2d[:, 0], np.full(n_nodes_2d, span), nodes_2d[:, 1]])
    nodes = np.vstack([nodes_y0, nodes_y1])

    # Hex elements: (n_c - 1) cells around C × (n_layers - 1) cells radially
    n_cells_c = n_c - 1  # open C: no wraparound
    n_cells_normal = n_layers - 1
    n_hexes = n_cells_c * n_cells_normal

    hexes = np.zeros((n_hexes, 8), dtype=np.int32)
    idx = 0

    for j in range(n_cells_normal):
        for i in range(n_cells_c):
            n0 = j * n_c + i
            n1 = j * n_c + (i + 1)
            n2 = (j + 1) * n_c + (i + 1)
            n3 = (j + 1) * n_c + i
            n4 = n0 + n_nodes_2d
            n5 = n1 + n_nodes_2d
            n6 = n2 + n_nodes_2d
            n7 = n3 + n_nodes_2d

            # Check winding via scalar triple product
            p0, p1, p3, p4 = nodes_y0[n0], nodes_y0[n1], nodes_y0[n3], nodes_y1[n0]
            vol = np.dot(p1 - p0, np.cross(p3 - p0, p4 - p0))

            if vol > 0:
                hexes[idx] = [n0, n1, n2, n3, n4, n5, n6, n7]
            else:
                hexes[idx] = [n0, n3, n2, n1, n4, n7, n6, n5]
            idx += 1

    hexes = hexes[:idx] + 1  # 1-based for UGRID

    # --- Boundary faces ---
    n_wake = mesh_2d["n_wake"]

    # 1. Wall: inner ring (j=0), but ONLY the airfoil portion, not the wake
    # C-boundary indices: [0..n_wake-1] = wake_lower, [n_wake..n_c-n_wake-1] = airfoil, [n_c-n_wake..n_c-1] = wake_upper
    airfoil_start = n_wake
    airfoil_end = n_c - n_wake
    wall_quads = []
    for i in range(airfoil_start, airfoil_end - 1):
        n0, n1 = i, i + 1
        wall_quads.append([n0, n0 + n_nodes_2d, n1 + n_nodes_2d, n1])

    # 2. Farfield: outer ring (j = n_layers-1), ALL cells
    farfield_quads = []
    j = n_cells_normal
    for i in range(n_cells_c):
        n0 = j * n_c + i
        n1 = j * n_c + (i + 1)
        farfield_quads.append([n0, n1, n1 + n_nodes_2d, n0 + n_nodes_2d])

    # 3. Wake cut: the two open ends of the C-grid, stacked radially
    # Wake end 1 (i=0 face): all layers, first point of each layer
    wake_end1_quads = []
    for j in range(n_cells_normal):
        n0 = j * n_c + 0  # first point
        n3 = (j + 1) * n_c + 0
        wake_end1_quads.append([n0, n0 + n_nodes_2d, n3 + n_nodes_2d, n3])

    # Wake end 2 (i=n_c-1 face): all layers, last point of each layer
    wake_end2_quads = []
    for j in range(n_cells_normal):
        n0 = j * n_c + (n_c - 1)  # last point
        n3 = (j + 1) * n_c + (n_c - 1)
        wake_end2_quads.append([n0, n3, n3 + n_nodes_2d, n0 + n_nodes_2d])

    # 4. Symmetry faces (y=0 and y=span): ALL 2D quad cells
    sym_y0_quads = []
    sym_y1_quads = []
    for j in range(n_cells_normal):
        for i in range(n_cells_c):
            n0 = j * n_c + i
            n1 = j * n_c + (i + 1)
            n2 = (j + 1) * n_c + (i + 1)
            n3 = (j + 1) * n_c + i
            sym_y0_quads.append([n0, n3, n2, n1])
            sym_y1_quads.append([n0 + n_nodes_2d, n1 + n_nodes_2d,
                                 n2 + n_nodes_2d, n3 + n_nodes_2d])

    boundary_quads = {
        "wall": np.array(wall_quads, dtype=np.int32) + 1,
        "farfield": np.array(farfield_quads, dtype=np.int32) + 1,
        "symmetry_y0": np.array(sym_y0_quads, dtype=np.int32) + 1,
        "symmetry_y1": np.array(sym_y1_quads, dtype=np.int32) + 1,
    }

    # Wake cut faces get Freestream BC (outflow)
    wake_faces = np.vstack([
        np.array(wake_end1_quads, dtype=np.int32),
        np.array(wake_end2_quads, dtype=np.int32),
    ]) + 1
    boundary_quads["wake"] = wake_faces

    boundary_ids = {
        "wall": 1,
        "farfield": 2,
        "symmetry_y0": 3,
        "symmetry_y1": 4,
        "wake": 5,
    }

    return {
        "nodes": nodes,
        "hexes": hexes,
        "boundary_quads": boundary_quads,
        "boundary_ids": boundary_ids,
        "n_nodes": len(nodes),
        "n_hexes": len(hexes),
    }


# ---------------------------------------------------------------------------
# UGRID writer
# ---------------------------------------------------------------------------

def write_ugrid(path: str | Path, mesh_3d: dict) -> None:
    """Write a 3D hex mesh in AFLR3/UGRID big-endian binary format (.b8.ugrid)."""
    path = Path(path)
    nodes = mesh_3d["nodes"]
    hexes = mesh_3d["hexes"]
    boundary_quads = mesh_3d["boundary_quads"]
    boundary_ids = mesh_3d["boundary_ids"]

    all_bquads = []
    all_tags = []
    for name, quads in boundary_quads.items():
        tag = boundary_ids[name]
        all_bquads.append(quads)
        all_tags.extend([tag] * len(quads))
    all_bquads = np.vstack(all_bquads) if all_bquads else np.zeros((0, 4), dtype=np.int32)
    all_tags = np.array(all_tags, dtype=np.int32)

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
    boundary_ids = mesh_3d["boundary_ids"]
    bc_type_map = {
        "wall": 4000,
        "farfield": 3000,
        "symmetry_y0": 5000,
        "symmetry_y1": 5000,
        "wake": 3000,  # Freestream on wake cut (outflow)
    }
    lines = [str(len(boundary_ids))]
    for name, tag in sorted(boundary_ids.items(), key=lambda x: x[1]):
        lines.append(f"{tag} {bc_type_map.get(name, 0)} {name}")
    path.write_text("\n".join(lines) + "\n")


# ---------------------------------------------------------------------------
# Top-level pipeline
# ---------------------------------------------------------------------------

def generate_and_write_mesh(
    coords: list[tuple[float, float]] | list[list[tuple[float, float]]],
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
    """Full pipeline: airfoil coords → UGRID mesh files."""
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
