"""Structured O-grid mesh generation for pseudo-2D airfoil RANS.

Generates a single-cell-deep hex mesh from 2D airfoil coordinates using a
hybrid normal-offset / TFI O-grid. The inner BL layers use normal-offset for
proper wall resolution; the outer layers blend to transfinite interpolation
(surface → farfield circle) for guaranteed cell validity.

Supports single-body and multi-body (slat + main + flap) configurations.

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
    """Compute unit outward normals along a closed curve."""
    n = len(curve)
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


def _make_farfield_ring(surface: np.ndarray, center: np.ndarray, radius: float) -> np.ndarray:
    """Create a farfield circle with points matched to surface point directions."""
    directions = surface - center
    angles = np.arctan2(directions[:, 1], directions[:, 0])
    ring = np.column_stack([
        center[0] + radius * np.cos(angles),
        center[1] + radius * np.sin(angles),
    ])
    return ring


# ---------------------------------------------------------------------------
# O-grid generation
# ---------------------------------------------------------------------------

def generate_airfoil_mesh(
    coords: list[tuple[float, float]] | list[list[tuple[float, float]]],
    *,
    n_normal: int = 80,
    first_cell_height: float | None = None,
    Re: float = 1e6,
    growth_rate: float = 1.1,
    farfield_radius: float = 50.0,
    chord: float = 1.0,
    y_plus: float = 1.0,
) -> dict:
    """Generate a structured O-grid mesh around one or more 2D airfoil bodies.

    Uses a hybrid approach: normal-offset for the inner boundary-layer region,
    smoothly blending to transfinite interpolation (TFI) towards the farfield
    circle. This gives proper y+ wall resolution while guaranteeing positive
    cell volumes everywhere.

    Parameters
    ----------
    coords : list of (x, y), or list of lists for multi-body
        Single body: list of (x, y) in Selig ordering.
        Multi-body: list of bodies, each a list of (x, y) in Selig ordering,
                    ordered upstream to downstream (e.g. [slat, main, flap]).
    n_normal : int
        Number of cells in the wall-normal direction.
    first_cell_height : float or None
        First cell height. Auto-estimated from Re if None.
    Re : float
        Reynolds number.
    growth_rate : float
        Geometric growth rate for BL layers.
    farfield_radius : float
        Farfield distance in chord lengths.
    chord : float
        Reference chord length.
    y_plus : float
        Target y+.

    Returns
    -------
    dict with mesh data.
    """
    # Normalize input: single body → list of one body
    if (len(coords) > 0
            and isinstance(coords[0], (tuple, list))
            and len(coords[0]) == 2
            and isinstance(coords[0][0], (int, float))):
        bodies = [np.array(coords, dtype=np.float64)]
    else:
        bodies = [np.array(b, dtype=np.float64) for b in coords]

    # Concatenate all bodies into one contour
    # For multi-body, the contour goes: body1 upper→LE→lower, body2 upper→LE→lower, ...
    surface = np.vstack(bodies)
    n_surface = len(surface)

    if n_surface < 20:
        raise ValueError(f"Need at least 20 surface points, got {n_surface}")

    # Ensure closed contour
    if np.linalg.norm(surface[0] - surface[-1]) > 1e-10:
        surface = np.vstack([surface, surface[0]])
        n_surface = len(surface)

    # Track body ranges
    body_ranges = []
    offset = 0
    for body in bodies:
        body_ranges.append((offset, offset + len(body)))
        offset += len(body)

    if first_cell_height is None:
        first_cell_height = estimate_first_cell_height(Re, chord, y_plus)

    center = np.array([0.5 * chord, 0.0])
    normals = _compute_normals(surface)
    farfield = _make_farfield_ring(surface, center, farfield_radius * chord)

    # Generate radial layer heights with geometric growth
    n_layers = n_normal + 1
    heights = np.zeros(n_layers)
    for i in range(1, n_layers):
        heights[i] = heights[i - 1] + first_cell_height * growth_rate ** (i - 1)

    total_height = heights[-1]
    s = heights / total_height if total_height > 0 else np.linspace(0, 1, n_layers)

    # Hybrid normal-offset / TFI
    # Inner layers use normal-offset for BL resolution; outer layers blend to TFI.
    # Transition at ~10% chord to keep more of the near-field properly resolved.
    bl_transition = 0.10 * chord

    nodes_2d = np.zeros((n_surface * n_layers, 2))
    nodes_2d[:n_surface] = surface

    for j in range(1, n_layers):
        h = heights[j]

        normal_layer = surface + h * normals
        tfi_layer = (1.0 - s[j]) * surface + s[j] * farfield

        if h <= bl_transition:
            blend = 0.0
        elif h >= bl_transition * 4.0:
            blend = 1.0
        else:
            t = (h - bl_transition) / (bl_transition * 3.0)
            blend = 3 * t * t - 2 * t * t * t  # smoothstep

        layer = (1.0 - blend) * normal_layer + blend * tfi_layer
        nodes_2d[j * n_surface: (j + 1) * n_surface] = layer

    return {
        "nodes_2d": nodes_2d,
        "n_surface": n_surface,
        "n_layers": n_layers,
        "surface_coords": surface,
        "first_cell_height": first_cell_height,
        "body_ranges": body_ranges,
    }


# ---------------------------------------------------------------------------
# 3D extrusion
# ---------------------------------------------------------------------------

def extrude_to_3d(mesh_2d: dict, *, span: float = 0.01) -> dict:
    """Extrude a 2D O-grid mesh one cell deep in the y-direction.

    The airfoil lies in the x-z plane (chord along x, thickness along z).
    Span (extrusion) is along y. This matches Flow360's alphaAngle convention.
    """
    nodes_2d = mesh_2d["nodes_2d"]
    n_surface = mesh_2d["n_surface"]
    n_layers = mesh_2d["n_layers"]
    n_nodes_2d = len(nodes_2d)

    # 3D nodes: (x, y_span, z=2D_y)
    nodes_y0 = np.column_stack([nodes_2d[:, 0], np.zeros(n_nodes_2d), nodes_2d[:, 1]])
    nodes_y1 = np.column_stack([nodes_2d[:, 0], np.full(n_nodes_2d, span), nodes_2d[:, 1]])
    nodes = np.vstack([nodes_y0, nodes_y1])

    # Hex elements
    n_cells_circ = n_surface - 1  # closed contour
    n_cells_normal = n_layers - 1
    n_hexes = n_cells_circ * n_cells_normal

    hexes = np.zeros((n_hexes, 8), dtype=np.int32)
    idx = 0

    for j in range(n_cells_normal):
        for i in range(n_cells_circ):
            i_next = (i + 1) % (n_surface - 1)  # wrap for closed contour

            n0 = j * n_surface + i
            n1 = j * n_surface + i_next
            n2 = (j + 1) * n_surface + i_next
            n3 = (j + 1) * n_surface + i
            n4 = n0 + n_nodes_2d
            n5 = n1 + n_nodes_2d
            n6 = n2 + n_nodes_2d
            n7 = n3 + n_nodes_2d

            # Check winding
            p0, p1, p3, p4 = nodes_y0[n0], nodes_y0[n1], nodes_y0[n3], nodes_y1[n0]
            vol = np.dot(p1 - p0, np.cross(p3 - p0, p4 - p0))

            if vol > 0:
                hexes[idx] = [n0, n1, n2, n3, n4, n5, n6, n7]
            else:
                hexes[idx] = [n0, n3, n2, n1, n4, n7, n6, n5]
            idx += 1

    hexes = hexes[:idx] + 1  # 1-based

    # Boundary faces
    # 1. Wall: inner ring (j=0)
    wall_quads = []
    for i in range(n_cells_circ):
        i_next = (i + 1) % (n_surface - 1)
        n0, n1 = i, i_next
        wall_quads.append([n0, n0 + n_nodes_2d, n1 + n_nodes_2d, n1])

    # 2. Farfield: outer ring (j = n_layers-1)
    farfield_quads = []
    j = n_cells_normal
    for i in range(n_cells_circ):
        i_next = (i + 1) % (n_surface - 1)
        n0 = j * n_surface + i
        n1 = j * n_surface + i_next
        farfield_quads.append([n0, n1, n1 + n_nodes_2d, n0 + n_nodes_2d])

    # 3/4. Symmetry faces (y=0 and y=span)
    sym_y0_quads = []
    sym_y1_quads = []
    for j in range(n_cells_normal):
        for i in range(n_cells_circ):
            i_next = (i + 1) % (n_surface - 1)
            n0 = j * n_surface + i
            n1 = j * n_surface + i_next
            n2 = (j + 1) * n_surface + i_next
            n3 = (j + 1) * n_surface + i
            sym_y0_quads.append([n0, n3, n2, n1])
            sym_y1_quads.append([n0 + n_nodes_2d, n1 + n_nodes_2d,
                                 n2 + n_nodes_2d, n3 + n_nodes_2d])

    boundary_quads = {
        "wall": np.array(wall_quads, dtype=np.int32) + 1,
        "farfield": np.array(farfield_quads, dtype=np.int32) + 1,
        "symmetry_y0": np.array(sym_y0_quads, dtype=np.int32) + 1,
        "symmetry_y1": np.array(sym_y1_quads, dtype=np.int32) + 1,
    }

    boundary_ids = {
        "wall": 1,
        "farfield": 2,
        "symmetry_y0": 3,
        "symmetry_y1": 4,
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
        "wall": 4000, "farfield": 3000,
        "symmetry_y0": 5000, "symmetry_y1": 5000,
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
    growth_rate: float = 1.1,
    farfield_radius: float = 50.0,
    span: float = 0.01,
    mesh_name: str = "airfoil",
) -> tuple[Path, Path]:
    """Full pipeline: airfoil coords → UGRID mesh files."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh_2d = generate_airfoil_mesh(
        coords, n_normal=n_normal, Re=Re,
        growth_rate=growth_rate, farfield_radius=farfield_radius,
    )
    mesh_3d = extrude_to_3d(mesh_2d, span=span)

    ugrid_path = output_dir / f"{mesh_name}.b8.ugrid"
    mapbc_path = output_dir / f"{mesh_name}.mapbc"

    write_ugrid(ugrid_path, mesh_3d)
    write_mapbc(mapbc_path, mesh_3d)

    return ugrid_path, mapbc_path
