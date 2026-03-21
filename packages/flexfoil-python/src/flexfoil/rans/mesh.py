"""Structured mesh generation for pseudo-2D airfoil RANS.

Generates a single-cell-deep hex mesh from 2D airfoil coordinates:
  airfoil coords → structured O-grid (TFI) → extrude to 3D hexes → write UGRID.

Uses Transfinite Interpolation between the airfoil surface and a farfield circle
to guarantee positive cell volumes.

No external meshing dependencies — pure numpy.
"""

from __future__ import annotations

import struct
from pathlib import Path

import numpy as np


def estimate_first_cell_height(
    Re: float, chord: float = 1.0, y_plus: float = 1.0
) -> float:
    """Estimate first cell height for a target y+ using flat-plate correlation.

    Uses the Schlichting skin-friction formula:
        Cf = 0.058 * Re^(-0.2)
        u_tau = sqrt(Cf / 2) * U_inf
        y = y+ * nu / u_tau

    In nondimensional form (chord = 1):
        y1 = y+ * chord / (Re * sqrt(Cf / 2))
    """
    cf = 0.058 * Re ** (-0.2)
    u_tau_norm = np.sqrt(cf / 2.0)
    y1 = y_plus * chord / (Re * u_tau_norm)
    return float(y1)


def _compute_normals(surface: np.ndarray) -> np.ndarray:
    """Compute unit outward normals at each point of a closed airfoil contour."""
    n = len(surface)
    tangents = np.zeros_like(surface)
    tangents[1:-1] = surface[2:] - surface[:-2]
    tangents[0] = surface[1] - surface[0]
    tangents[-1] = surface[-1] - surface[-2]

    lengths = np.linalg.norm(tangents, axis=1, keepdims=True)
    lengths = np.maximum(lengths, 1e-14)
    tangents = tangents / lengths

    # Rotate 90° to get normal: (tx, ty) → (ty, -tx)
    normals = np.column_stack([tangents[:, 1], -tangents[:, 0]])

    # Ensure outward: normals should point away from centroid
    centroid = surface.mean(axis=0)
    outward_check = surface - centroid
    dots = np.sum(normals * outward_check, axis=1)
    if np.sum(dots < 0) > np.sum(dots > 0):
        normals = -normals

    return normals


def _make_farfield_ring(surface: np.ndarray, center: np.ndarray, radius: float) -> np.ndarray:
    """Create a farfield circle matching the angular distribution of surface points.

    Each farfield point is placed on a circle in the direction from center to
    the corresponding surface point, ensuring smooth mesh lines.
    """
    n_pts = len(surface)
    directions = surface - center
    angles = np.arctan2(directions[:, 1], directions[:, 0])
    ring = np.column_stack([
        center[0] + radius * np.cos(angles),
        center[1] + radius * np.sin(angles),
    ])
    return ring


def generate_airfoil_mesh(
    coords: list[tuple[float, float]],
    *,
    n_normal: int = 64,
    first_cell_height: float | None = None,
    Re: float = 1e6,
    growth_rate: float = 1.15,
    farfield_radius: float = 100.0,
    chord: float = 1.0,
    y_plus: float = 1.0,
) -> dict:
    """Generate a structured O-grid mesh around a 2D airfoil using TFI.

    Uses transfinite interpolation between the airfoil surface and a farfield
    circle, with geometric stretching in the wall-normal direction for
    boundary-layer resolution.

    Parameters
    ----------
    coords : list of (x, y) tuples
        Airfoil coordinates, Selig ordering (upper TE → LE → lower TE).
    n_normal : int
        Number of cells in the wall-normal direction.
    first_cell_height : float or None
        Height of the first cell off the wall. If None, estimated from Re and y_plus.
    Re : float
        Reynolds number (used to estimate first_cell_height if not given).
    growth_rate : float
        Geometric growth rate for cell layers.
    farfield_radius : float
        Farfield distance in chord lengths.
    chord : float
        Chord length.
    y_plus : float
        Target y+ for the first cell.

    Returns
    -------
    dict with keys:
        'nodes_2d': (N, 2) array of all 2D mesh nodes
        'n_surface': number of surface points
        'n_layers': number of layers (n_normal + 1)
        'surface_coords': (N_surface, 2)
        'first_cell_height': float
    """
    surface = np.array(coords, dtype=np.float64)
    n_surface = len(surface)

    if n_surface < 20:
        raise ValueError(f"Need at least 20 surface points, got {n_surface}")

    # Ensure the contour is closed (last point == first point)
    if np.linalg.norm(surface[0] - surface[-1]) > 1e-10:
        surface = np.vstack([surface, surface[0]])
        n_surface = len(surface)

    # Estimate first cell height
    if first_cell_height is None:
        first_cell_height = estimate_first_cell_height(Re, chord, y_plus)

    center = np.array([0.5 * chord, 0.0])
    normals = _compute_normals(surface)

    # Create farfield ring matched to surface point angles
    farfield = _make_farfield_ring(surface, center, farfield_radius * chord)

    # Build layers using a hybrid approach:
    # - Inner layers (BL region): normal-offset from surface with geometric growth
    # - Outer layers: blend to TFI (surface→farfield) for guaranteed validity
    #
    # The transition height is set to ~5% chord — within this zone, normal-offset
    # is safe because the BL thickness is small relative to curvature. Beyond it,
    # we blend to TFI to avoid TE crossing issues.
    n_layers = n_normal + 1
    bl_transition = 0.05 * chord  # offset distance where we start blending to TFI

    # Compute cumulative wall-normal distances with geometric growth
    heights = np.zeros(n_layers)
    for i in range(1, n_layers):
        heights[i] = heights[i - 1] + first_cell_height * growth_rate ** (i - 1)

    # Normalize heights to [0, 1] for TFI parameter
    total_height = heights[-1]
    s = heights / total_height if total_height > 0 else np.linspace(0, 1, n_layers)

    nodes_2d = np.zeros((n_surface * n_layers, 2))
    nodes_2d[:n_surface] = surface  # layer 0 = surface

    for j in range(1, n_layers):
        h = heights[j]

        # Normal-offset layer (ideal for BL resolution)
        normal_layer = surface + h * normals

        # TFI layer (guaranteed valid, but poor BL resolution)
        tfi_layer = (1.0 - s[j]) * surface + s[j] * farfield

        # Blending weight: 0 = pure normal-offset, 1 = pure TFI
        # Smooth transition between bl_transition and 3*bl_transition
        if h <= bl_transition:
            blend = 0.0
        elif h >= bl_transition * 3.0:
            blend = 1.0
        else:
            t = (h - bl_transition) / (bl_transition * 2.0)
            blend = 3 * t * t - 2 * t * t * t  # smoothstep

        layer = (1.0 - blend) * normal_layer + blend * tfi_layer
        nodes_2d[j * n_surface: (j + 1) * n_surface] = layer

    return {
        "nodes_2d": nodes_2d,
        "n_surface": n_surface,
        "n_layers": n_layers,
        "surface_coords": surface,
        "first_cell_height": first_cell_height,
    }


def extrude_to_3d(
    mesh_2d: dict,
    *,
    span: float = 0.01,
) -> dict:
    """Extrude a 2D O-grid mesh one cell deep in the z-direction.

    Returns a dict with:
        'nodes': (N, 3) float64 — all 3D node coordinates
        'hexes': (M, 8) int32 — hex element connectivity (1-based for UGRID)
        'boundary_quads': dict mapping boundary name → (K, 4) int32 connectivity
        'boundary_ids': dict mapping boundary name → integer BC tag
    """
    nodes_2d = mesh_2d["nodes_2d"]
    n_surface = mesh_2d["n_surface"]
    n_layers = mesh_2d["n_layers"]
    n_nodes_2d = len(nodes_2d)

    # 3D nodes: airfoil in x-z plane, span in y-direction
    # Flow360 alpha rotates freestream in x-z plane, so the airfoil
    # chord must be along x, thickness along z, span along y.
    # 2D coords are (x, y_2d) → 3D coords are (x, y_span, z=y_2d)
    nodes_y0 = np.column_stack([nodes_2d[:, 0], np.zeros(n_nodes_2d), nodes_2d[:, 1]])
    nodes_y1 = np.column_stack([nodes_2d[:, 0], np.full(n_nodes_2d, span), nodes_2d[:, 1]])
    nodes = np.vstack([nodes_y0, nodes_y1])

    # Build hex elements
    # Closed contour: n_cells_circ = n_surface - 1 (last point == first)
    n_cells_circ = n_surface - 1
    n_cells_normal = n_layers - 1
    n_hexes = n_cells_circ * n_cells_normal

    hexes = np.zeros((n_hexes, 8), dtype=np.int32)
    idx = 0

    for j in range(n_cells_normal):
        for i in range(n_cells_circ):
            i_next = i + 1
            # Wrap around for closed contour
            if i_next >= n_surface - 1:
                i_next = 0

            # Bottom face (z=0): quad nodes
            n0 = j * n_surface + i
            n1 = j * n_surface + i_next
            n2 = (j + 1) * n_surface + i_next
            n3 = (j + 1) * n_surface + i

            # Top face (z=span)
            n4 = n0 + n_nodes_2d
            n5 = n1 + n_nodes_2d
            n6 = n2 + n_nodes_2d
            n7 = n3 + n_nodes_2d

            # Determine correct winding for positive volume
            # Extrusion is in +y direction; check using scalar triple product
            p0 = nodes_y0[n0]
            p1 = nodes_y0[n1]
            p3 = nodes_y0[n3]
            p4 = nodes_y1[n0]
            d1 = p1 - p0
            d2 = p3 - p0
            d3 = p4 - p0
            cross_z = np.dot(d1, np.cross(d2, d3))

            if cross_z > 0:
                # Bottom face is CCW (correct for +z extrusion)
                hexes[idx] = [n0, n1, n2, n3, n4, n5, n6, n7]
            else:
                # Reverse winding
                hexes[idx] = [n0, n3, n2, n1, n4, n7, n6, n5]
            idx += 1

    hexes = hexes[:idx]

    # Convert to 1-based indexing for UGRID
    hexes_1based = hexes + 1

    # Build boundary face quads

    # 1. Airfoil wall: inner ring (j=0)
    wall_quads = []
    for i in range(n_cells_circ):
        i_next = i + 1
        if i_next >= n_surface - 1:
            i_next = 0
        n0 = i
        n1 = i_next
        n4 = n0 + n_nodes_2d
        n5 = n1 + n_nodes_2d
        wall_quads.append([n0, n4, n5, n1])

    # 2. Farfield: outer ring (j = n_layers-1)
    farfield_quads = []
    j = n_cells_normal
    for i in range(n_cells_circ):
        i_next = i + 1
        if i_next >= n_surface - 1:
            i_next = 0
        n0 = j * n_surface + i
        n1 = j * n_surface + i_next
        n4 = n0 + n_nodes_2d
        n5 = n1 + n_nodes_2d
        farfield_quads.append([n0, n1, n5, n4])

    # 3. Symmetry z=0: all quads at z=0 plane
    sym_z0_quads = []
    for j in range(n_cells_normal):
        for i in range(n_cells_circ):
            i_next = i + 1
            if i_next >= n_surface - 1:
                i_next = 0
            n0 = j * n_surface + i
            n1 = j * n_surface + i_next
            n2 = (j + 1) * n_surface + i_next
            n3 = (j + 1) * n_surface + i
            sym_z0_quads.append([n0, n3, n2, n1])

    # 4. Symmetry z=span: all quads at z=span plane
    sym_z1_quads = []
    for j in range(n_cells_normal):
        for i in range(n_cells_circ):
            i_next = i + 1
            if i_next >= n_surface - 1:
                i_next = 0
            n0 = j * n_surface + i + n_nodes_2d
            n1 = j * n_surface + i_next + n_nodes_2d
            n2 = (j + 1) * n_surface + i_next + n_nodes_2d
            n3 = (j + 1) * n_surface + i + n_nodes_2d
            sym_z1_quads.append([n0, n1, n2, n3])

    boundary_quads = {
        "wall": np.array(wall_quads, dtype=np.int32) + 1,
        "farfield": np.array(farfield_quads, dtype=np.int32) + 1,
        "symmetry_z0": np.array(sym_z0_quads, dtype=np.int32) + 1,
        "symmetry_z1": np.array(sym_z1_quads, dtype=np.int32) + 1,
    }

    boundary_ids = {
        "wall": 1,
        "farfield": 2,
        "symmetry_z0": 3,
        "symmetry_z1": 4,
    }

    return {
        "nodes": nodes,
        "hexes": hexes_1based,
        "boundary_quads": boundary_quads,
        "boundary_ids": boundary_ids,
        "n_nodes": len(nodes),
        "n_hexes": len(hexes_1based),
    }


def write_ugrid(path: str | Path, mesh_3d: dict) -> None:
    """Write a 3D hex mesh in AFLR3/UGRID binary format (.b8.ugrid).

    Format: big-endian 64-bit floats, 32-bit ints.
    """
    path = Path(path)
    nodes = mesh_3d["nodes"]
    hexes = mesh_3d["hexes"]
    boundary_quads = mesh_3d["boundary_quads"]
    boundary_ids = mesh_3d["boundary_ids"]

    n_nodes = len(nodes)
    n_tris = 0
    n_tets = 0
    n_pyramids = 0
    n_prisms = 0
    n_hexes = len(hexes)

    # Collect all boundary quads with tags
    all_bquads = []
    all_bquad_tags = []
    for name, quads in boundary_quads.items():
        tag = boundary_ids[name]
        all_bquads.append(quads)
        all_bquad_tags.extend([tag] * len(quads))
    all_bquads = np.vstack(all_bquads) if all_bquads else np.zeros((0, 4), dtype=np.int32)
    all_bquad_tags = np.array(all_bquad_tags, dtype=np.int32)
    n_surf_quads = len(all_bquads)

    with open(path, "wb") as f:
        # Header
        f.write(struct.pack(">7i", n_nodes, n_tris, n_surf_quads, n_tets, n_pyramids, n_prisms, n_hexes))

        # Node coordinates
        for node in nodes:
            f.write(struct.pack(">3d", *node))

        # Surface quads (4 node indices each, 1-based)
        for quad in all_bquads:
            f.write(struct.pack(">4i", *quad))

        # Surface quad tags
        for tag in all_bquad_tags:
            f.write(struct.pack(">i", tag))

        # Volume hexes (8 node indices each, 1-based)
        for hex_elem in hexes:
            f.write(struct.pack(">8i", *hex_elem))


def write_mapbc(path: str | Path, mesh_3d: dict) -> None:
    """Write .mapbc boundary condition mapping file."""
    path = Path(path)
    boundary_ids = mesh_3d["boundary_ids"]

    bc_type_map = {
        "wall": 4000,
        "farfield": 3000,
        "symmetry_z0": 5000,
        "symmetry_z1": 5000,
    }

    n_boundaries = len(boundary_ids)
    lines = [str(n_boundaries)]
    for name, tag in sorted(boundary_ids.items(), key=lambda x: x[1]):
        bc_type = bc_type_map.get(name, 0)
        lines.append(f"{tag} {bc_type} {name}")

    path.write_text("\n".join(lines) + "\n")


def generate_and_write_mesh(
    coords: list[tuple[float, float]],
    output_dir: str | Path,
    *,
    Re: float = 1e6,
    n_normal: int = 64,
    growth_rate: float = 1.15,
    farfield_radius: float = 100.0,
    span: float = 0.01,
    mesh_name: str = "airfoil",
) -> tuple[Path, Path]:
    """Full pipeline: airfoil coords → UGRID mesh files.

    Returns (ugrid_path, mapbc_path).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    mesh_2d = generate_airfoil_mesh(
        coords,
        n_normal=n_normal,
        Re=Re,
        growth_rate=growth_rate,
        farfield_radius=farfield_radius,
    )

    mesh_3d = extrude_to_3d(mesh_2d, span=span)

    ugrid_path = output_dir / f"{mesh_name}.b8.ugrid"
    mapbc_path = output_dir / f"{mesh_name}.mapbc"

    write_ugrid(ugrid_path, mesh_3d)
    write_mapbc(mapbc_path, mesh_3d)

    return ugrid_path, mapbc_path
