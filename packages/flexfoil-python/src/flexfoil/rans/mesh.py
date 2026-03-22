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

    # Tag airfoil surface faces (lateral faces of the extrusion)
    # After extrude, the body has 5 faces:
    #   - 3 lateral faces (upper surface, lower surface, TE) → airfoil wall
    #   - 2 end caps (y=0 and y=span) → these get deleted by AutomatedFarfield
    # We tag only the lateral faces as "airfoil"
    lines.append("select face")
    lines.append("attribute groupName $airfoil")
    lines.append("attribute faceName $airfoil")

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
    """Full pipeline: airfoil coords → UGRID mesh files.

    Uses gmsh if available (high-quality BL mesh), falls back to the
    algebraic C-grid generator otherwise.
    """
    try:
        return generate_and_write_mesh_gmsh(
            coords, output_dir, Re=Re, growth_rate=growth_rate,
            farfield_radius=farfield_radius, span=span, mesh_name=mesh_name,
        )
    except ImportError:
        pass

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


# ---------------------------------------------------------------------------
# gmsh-based mesh generation (preferred)
# ---------------------------------------------------------------------------

def generate_and_write_mesh_gmsh(
    coords: list[tuple[float, float]],
    output_dir: str | Path,
    *,
    Re: float = 1e6,
    growth_rate: float = 1.15,
    farfield_radius: float = 50.0,
    span: float = 0.01,
    mesh_name: str = "airfoil",
    y_plus: float = 1.0,
    bl_thickness: float = 0.2,
    surface_size_min: float = 0.01,
    surface_size_max: float = 5.0,
) -> tuple[Path, Path]:
    """Generate a pseudo-2D airfoil mesh using gmsh.

    Uses gmsh's BoundaryLayer field for proper BL meshing, then manually
    extrudes one cell deep in the spanwise direction. Produces hex + prism
    elements with no degenerate cells.

    Requires ``pip install gmsh``.

    Returns (ugrid_path, mapbc_path).
    """
    import gmsh

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    first_cell = estimate_first_cell_height(Re, y_plus=y_plus)

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("airfoil")

    # Airfoil spline
    pts = [gmsh.model.geo.addPoint(x, y, 0, surface_size_min) for x, y in coords]
    spline = gmsh.model.geo.addBSpline(pts + [pts[0]])
    airfoil_loop = gmsh.model.geo.addCurveLoop([spline])

    # Farfield circle
    R = farfield_radius
    cx, cy = 0.5, 0.0
    center = gmsh.model.geo.addPoint(cx, cy, 0)
    p1 = gmsh.model.geo.addPoint(cx + R, cy, 0, surface_size_max)
    p2 = gmsh.model.geo.addPoint(cx, cy + R, 0, surface_size_max)
    p3 = gmsh.model.geo.addPoint(cx - R, cy, 0, surface_size_max)
    p4 = gmsh.model.geo.addPoint(cx, cy - R, 0, surface_size_max)
    a1 = gmsh.model.geo.addCircleArc(p1, center, p2)
    a2 = gmsh.model.geo.addCircleArc(p2, center, p3)
    a3 = gmsh.model.geo.addCircleArc(p3, center, p4)
    a4 = gmsh.model.geo.addCircleArc(p4, center, p1)
    farfield_loop = gmsh.model.geo.addCurveLoop([a1, a2, a3, a4])

    surf = gmsh.model.geo.addPlaneSurface([farfield_loop, airfoil_loop])
    gmsh.model.geo.synchronize()

    # BL field
    bl = gmsh.model.mesh.field.add("BoundaryLayer")
    gmsh.model.mesh.field.setNumbers(bl, "CurvesList", [spline])
    gmsh.model.mesh.field.setNumber(bl, "Size", first_cell)
    gmsh.model.mesh.field.setNumber(bl, "Ratio", growth_rate)
    gmsh.model.mesh.field.setNumber(bl, "Thickness", bl_thickness)
    gmsh.model.mesh.field.setNumber(bl, "Quads", 1)
    gmsh.model.mesh.field.setAsBoundaryLayer(bl)

    # Size control
    dist = gmsh.model.mesh.field.add("Distance")
    gmsh.model.mesh.field.setNumbers(dist, "CurvesList", [spline])
    thresh = gmsh.model.mesh.field.add("Threshold")
    gmsh.model.mesh.field.setNumber(thresh, "InField", dist)
    gmsh.model.mesh.field.setNumber(thresh, "SizeMin", surface_size_min)
    gmsh.model.mesh.field.setNumber(thresh, "SizeMax", surface_size_max)
    gmsh.model.mesh.field.setNumber(thresh, "DistMin", bl_thickness)
    gmsh.model.mesh.field.setNumber(thresh, "DistMax", farfield_radius * 0.6)
    gmsh.model.mesh.field.setAsBackgroundMesh(thresh)

    gmsh.option.setNumber("Mesh.RecombineAll", 1)
    gmsh.option.setNumber("Mesh.Algorithm", 8)

    gmsh.model.mesh.generate(2)

    # Extract 2D mesh
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    n_nodes_2d = len(node_tags)
    coords_2d = node_coords.reshape(-1, 3)[:, :2]

    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}

    # Get elements
    elem_types, elem_tags, elem_nodes = gmsh.model.mesh.getElements(dim=2)
    quad_conn, tri_conn = None, None
    for et, _, enodes in zip(elem_types, elem_tags, elem_nodes):
        props = gmsh.model.mesh.getElementProperties(et)
        if "Quad" in props[0]:
            quad_conn = enodes.reshape(-1, 4)
        elif "Tri" in props[0]:
            tri_conn = enodes.reshape(-1, 3)

    # Boundary edges
    wall_edges = []
    et, _, enodes = gmsh.model.mesh.getElements(dim=1, tag=spline)
    for en in enodes:
        wall_edges.append(en.reshape(-1, 2))
    wall_edges = np.vstack(wall_edges) if wall_edges else np.zeros((0, 2), dtype=int)

    ff_edges = []
    for arc in [a1, a2, a3, a4]:
        et, _, enodes = gmsh.model.mesh.getElements(dim=1, tag=arc)
        for en in enodes:
            ff_edges.append(en.reshape(-1, 2))
    ff_edges = np.vstack(ff_edges) if ff_edges else np.zeros((0, 2), dtype=int)

    gmsh.finalize()

    # === EXTRUDE TO 3D ===
    # Nodes: airfoil in x-z plane, span in y
    nodes_y0 = np.column_stack([coords_2d[:, 0], np.zeros(n_nodes_2d), coords_2d[:, 1]])
    nodes_y1 = np.column_stack([coords_2d[:, 0], np.full(n_nodes_2d, span), coords_2d[:, 1]])
    all_nodes = np.vstack([nodes_y0, nodes_y1])

    # Hex elements from quads — check orientation so volume is positive
    hexes = []
    if quad_conn is not None:
        for q in quad_conn:
            idx = [tag_to_idx[int(t)] for t in q]
            n0, n1, n2, n3 = idx
            # Check hex volume via cross product
            p0, p1, p3 = nodes_y0[n0], nodes_y0[n1], nodes_y0[n3]
            p4 = nodes_y1[n0]
            vol = np.dot(p1 - p0, np.cross(p3 - p0, p4 - p0))
            if vol > 0:
                hexes.append([n0+1, n1+1, n2+1, n3+1,
                             n0+n_nodes_2d+1, n1+n_nodes_2d+1, n2+n_nodes_2d+1, n3+n_nodes_2d+1])
            else:
                # Flip winding to get positive volume
                hexes.append([n0+1, n3+1, n2+1, n1+1,
                             n0+n_nodes_2d+1, n3+n_nodes_2d+1, n2+n_nodes_2d+1, n1+n_nodes_2d+1])

    # Prism elements from triangles — check orientation
    prisms = []
    if tri_conn is not None:
        for t in tri_conn:
            idx = [tag_to_idx[int(tag)] for tag in t]
            n0, n1, n2 = idx
            p0, p1, p2 = nodes_y0[n0], nodes_y0[n1], nodes_y0[n2]
            p3 = nodes_y1[n0]
            vol = np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0))
            if vol > 0:
                prisms.append([n0+1, n1+1, n2+1, n0+n_nodes_2d+1, n1+n_nodes_2d+1, n2+n_nodes_2d+1])
            else:
                prisms.append([n0+1, n2+1, n1+1, n0+n_nodes_2d+1, n2+n_nodes_2d+1, n1+n_nodes_2d+1])

    hexes = np.array(hexes, dtype=np.int32) if hexes else np.zeros((0, 8), dtype=np.int32)
    prisms = np.array(prisms, dtype=np.int32) if prisms else np.zeros((0, 6), dtype=np.int32)

    # Boundary faces
    wall_quads = np.array([
        [tag_to_idx[int(e[0])]+1, tag_to_idx[int(e[1])]+1,
         tag_to_idx[int(e[1])]+n_nodes_2d+1, tag_to_idx[int(e[0])]+n_nodes_2d+1]
        for e in wall_edges
    ], dtype=np.int32)

    ff_quads = np.array([
        [tag_to_idx[int(e[0])]+1, tag_to_idx[int(e[1])]+1,
         tag_to_idx[int(e[1])]+n_nodes_2d+1, tag_to_idx[int(e[0])]+n_nodes_2d+1]
        for e in ff_edges
    ], dtype=np.int32)

    # Symmetry planes: all 2D quads at y=0 and y=span
    sym_y0_quads, sym_y0_tris = [], []
    sym_y1_quads, sym_y1_tris = [], []
    if quad_conn is not None:
        for q in quad_conn:
            idx = [tag_to_idx[int(t)]+1 for t in q]
            sym_y0_quads.append(idx)
            sym_y1_quads.append([i + n_nodes_2d for i in idx])
    if tri_conn is not None:
        for t in tri_conn:
            idx = [tag_to_idx[int(tag)]+1 for tag in t]
            sym_y0_tris.append(idx)
            sym_y1_tris.append([i + n_nodes_2d for i in idx])

    # === WRITE UGRID ===
    # Collect all boundary tris and quads with tags
    all_bnd_tris, all_tri_tags = [], []
    for tris, tag in [(sym_y0_tris, 3), (sym_y1_tris, 4)]:
        for t in tris:
            all_bnd_tris.append(t)
            all_tri_tags.append(tag)

    all_bnd_quads, all_quad_tags = [], []
    for quads, tag in [
        (wall_quads.tolist(), 1),
        (ff_quads.tolist(), 2),
        (sym_y0_quads, 3),
        (sym_y1_quads, 4),
    ]:
        for q in quads:
            all_bnd_quads.append(q)
            all_quad_tags.append(tag)

    ugrid_path = output_dir / f"{mesh_name}.b8.ugrid"
    mapbc_path = output_dir / f"{mesh_name}.mapbc"

    with open(ugrid_path, "wb") as f:
        f.write(struct.pack(">7i",
            len(all_nodes), len(all_bnd_tris), len(all_bnd_quads),
            0, 0, len(prisms), len(hexes)))
        for node in all_nodes:
            f.write(struct.pack(">3d", *node))
        for tri in all_bnd_tris:
            f.write(struct.pack(">3i", *tri))
        for quad in all_bnd_quads:
            f.write(struct.pack(">4i", *quad))
        for tag in all_tri_tags:
            f.write(struct.pack(">i", tag))
        for tag in all_quad_tags:
            f.write(struct.pack(">i", tag))
        for prism in prisms:
            f.write(struct.pack(">6i", *prism))
        for hex_elem in hexes:
            f.write(struct.pack(">8i", *hex_elem))

    mapbc_path.write_text("4\n1 4000 wall\n2 3000 farfield\n3 5000 symmetry_y0\n4 5000 symmetry_y1\n")

    return ugrid_path, mapbc_path
