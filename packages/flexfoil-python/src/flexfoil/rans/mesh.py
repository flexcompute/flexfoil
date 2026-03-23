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

    # Build layer spacing: geometric near wall, stretched to reach farfield
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

    # Normalize heights to [0, 1] for TFI parameter
    s = heights / heights[-1]  # s[0]=0 (surface), s[-1]=1 (farfield)

    # Build C-shaped outer boundary that preserves the wake opening.
    # The outer boundary is a semicircle from the lower wake tip, around
    # the front, to the upper wake tip — with straight segments along the
    # wake exit on both sides. This keeps the wake cut open at the farfield.
    R = farfield_radius * chord
    cx = 0.5 * chord  # circle center at mid-chord

    outer = np.zeros_like(c_boundary)

    # Wake tip positions at farfield distance
    wake_tip_x = c_boundary[0, 0]  # x of outermost wake point (lower)
    wake_lower_y = -R  # lower wake at farfield radius below
    wake_upper_y = R   # upper wake at farfield radius above

    for i in range(n_c):
        pt = c_boundary[i]

        if i < n_wake:
            # Lower wake: go straight down from wake point to farfield
            outer[i] = [pt[0], wake_lower_y]
        elif i >= n_c - n_wake:
            # Upper wake: go straight up from wake point to farfield
            outer[i] = [pt[0], wake_upper_y]
        else:
            # Airfoil portion: map to semicircle in front of wake
            # Parameterize along the airfoil contour
            i_airfoil = i - n_wake  # 0 = lower TE, n_airfoil-1 = upper TE
            n_airfoil = n_c - 2 * n_wake
            t_airfoil = i_airfoil / max(n_airfoil - 1, 1)  # 0 to 1

            # Map to angle: -pi/2 (lower TE) → -pi (LE) → -3pi/2 (upper TE)
            # i.e. the semicircle wrapping around the front
            angle = -np.pi / 2 - t_airfoil * np.pi

            outer[i] = [cx + R * np.cos(angle), R * np.sin(angle)]

    # TFI: blend inner (s=0) and outer (s=1) boundaries.
    # Convex combination guarantees no cell crossing.
    nodes_2d = np.zeros((n_c * n_layers, 2))
    for j in range(n_layers):
        nodes_2d[j * n_c: (j + 1) * n_c] = (1.0 - s[j]) * c_boundary + s[j] * outer

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
    span: float = 0.01,
    mesh_name: str = "airfoil",
    y_plus: float = 1.0,
    growth_rate: float = 1.12,
    n_bl_layers: int = 60,
    farfield_height: float = 20.0,
    wake_length: float = 10.0,
    airfoil_mesh_size: float = 0.005,
    farfield_mesh_size: float = 0.5,
    **kwargs,
) -> tuple[Path, Path]:
    """Generate a structured C-grid mesh using our forked gmshairfoil2d.

    Uses the gmshairfoil2d library (forked to preserve open trailing edges)
    to generate a 2D structured C-grid, then extrudes one cell deep in z
    and writes as UGRID binary.

    The airfoil lives in the x-y plane (gmsh convention). Flow360 uses
    betaAngle (not alphaAngle) for angle of attack with this orientation.

    Requires ``pip install gmsh``.

    Returns (ugrid_path, mapbc_path).
    """
    import gmsh
    from collections import Counter

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    first_cell = estimate_first_cell_height(Re, y_plus=y_plus)

    # Write airfoil coordinates to .dat file
    dat_path = output_dir / f"{mesh_name}.dat"
    with open(dat_path, "w") as f:
        f.write(f"{mesh_name}\n")
        for x, y in coords:
            f.write(f"  {x:.8f}  {y:.8f}\n")

    # Generate 2D structured C-grid using the forked gmshairfoil2d classes directly.
    # We bypass gmsh_main() and read_airfoil_from_file() to avoid their
    # point reordering issues — the coords are already in correct Selig order.
    from flexfoil.rans.gmshairfoil2d.geometry_def import AirfoilSpline, CType

    gmsh.initialize()
    gmsh.model.add("flexfoil_airfoil")
    gmsh.option.setNumber("General.Terminal", 0)

    # Create point cloud in (x, y, 0) format — preserving Selig ordering
    point_cloud = [(x, y, 0) for x, y in coords]

    airfoil = AirfoilSpline(point_cloud, airfoil_mesh_size, name="airfoil")

    dx_wake = 10.0  # wake extension length
    dy = farfield_height * 2  # total domain height

    mesh_obj = CType(
        airfoil, dx_wake, dy,
        airfoil_mesh_size, first_cell, growth_rate, aoa=0,
    )
    mesh_obj.define_bc()

    # Add physical surface for the fluid domain
    if mesh_obj.surfaces:
        ps = gmsh.model.addPhysicalGroup(2, mesh_obj.surfaces)
        gmsh.model.setPhysicalName(2, ps, "fluid")

    gmsh.model.geo.synchronize()

    # Write the 2D geo as .geo_unrolled, then re-open with extrusion appended.
    # This ensures gmsh handles shared block-interface nodes correctly
    # (direct extrude API creates duplicate nodes at block boundaries).
    geo_2d_path = output_dir / "airfoil_2d.geo_unrolled"
    gmsh.write(str(geo_2d_path))
    gmsh.finalize()

    # Append extrusion and re-generate as 3D
    geo_3d_path = output_dir / "airfoil_3d.geo"
    geo_text = geo_2d_path.read_text()
    geo_3d_path.write_text(
        geo_text + f"\nExtrude {{0, 0, {span}}} {{ Surface{{:}}; Layers{{1}}; Recombine; }}\n"
    )

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.open(str(geo_3d_path))
    gmsh.model.mesh.generate(3)

    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    all_coords = node_coords.reshape(-1, 3)
    n_nodes = len(node_tags)
    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}
    new_idx = {int(t): i + 1 for i, t in enumerate(node_tags)}

    hex_tags, hex_nodes = gmsh.model.mesh.getElementsByType(5)
    hex_conn = hex_nodes.reshape(-1, 8)
    n_hexes = len(hex_tags)

    gmsh.finalize()

    # Extract boundary faces from hex connectivity (guaranteed watertight)
    hex_face_defs = [
        (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6),
        (3, 0, 4, 7), (0, 3, 2, 1), (4, 5, 6, 7),
    ]
    face_count = Counter()
    face_to_nodes = {}
    for h in hex_conn:
        for fd in hex_face_defs:
            fk = tuple(sorted([int(h[i]) for i in fd]))
            face_count[fk] += 1
            if fk not in face_to_nodes:
                face_to_nodes[fk] = [int(h[i]) for i in fd]

    boundary_faces = [face_to_nodes[k] for k, c in face_count.items() if c == 1]

    # Classify boundary faces by position
    wall_q, ff_q, sym0_q, sym1_q = [], [], [], []
    for face in boundary_faces:
        indices = [tag_to_idx[n] for n in face]
        pts = all_coords[indices]
        z_range = pts[:, 2].max() - pts[:, 2].min()
        z_mean = pts[:, 2].mean()

        if z_range < 1e-8:
            # Flat in z → symmetry plane
            if abs(z_mean) < 1e-8:
                sym0_q.append(face)
            elif abs(z_mean - span) < 1e-8:
                sym1_q.append(face)
        else:
            # Lateral face → wall or farfield
            x_vals, y_vals = pts[:, 0], pts[:, 1]
            if x_vals.min() >= -0.01 and x_vals.max() <= 1.06 and abs(y_vals).max() < 0.1:
                wall_q.append(face)
            else:
                ff_q.append(face)

    total = len(wall_q) + len(ff_q) + len(sym0_q) + len(sym1_q)
    if total != len(boundary_faces):
        raise RuntimeError(
            f"Boundary classification mismatch: {total} classified vs "
            f"{len(boundary_faces)} extracted"
        )

    # Write UGRID
    ugrid_path = output_dir / f"{mesh_name}.b8.ugrid"
    mapbc_path = output_dir / f"{mesh_name}.mapbc"

    all_bnd_quads = wall_q + ff_q + sym0_q + sym1_q
    all_quad_tags = ([1] * len(wall_q) + [2] * len(ff_q) +
                     [3] * len(sym0_q) + [4] * len(sym1_q))

    with open(ugrid_path, "wb") as f:
        f.write(struct.pack(">7i", n_nodes, 0, len(all_bnd_quads), 0, 0, 0, n_hexes))
        for c in all_coords:
            f.write(struct.pack(">3d", *c))
        for q in all_bnd_quads:
            f.write(struct.pack(">4i", *[new_idx[n] for n in q]))
        for t in all_quad_tags:
            f.write(struct.pack(">i", t))
        for h in hex_conn:
            f.write(struct.pack(">8i", *[new_idx[int(n)] for n in h]))

    mapbc_path.write_text(
        "4\n1 4000 wall\n2 3000 farfield\n3 5000 symmetry_y0\n4 5000 symmetry_y1\n"
    )

    return ugrid_path, mapbc_path
