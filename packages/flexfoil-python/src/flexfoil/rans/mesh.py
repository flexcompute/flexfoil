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
    farfield_radius: float = 15.0,
    span: float = 0.01,
    mesh_name: str = "airfoil",
    y_plus: float = 1.0,
    n_airfoil: int = 200,
    n_normal: int = 125,
    n_wake: int = 150,
    n_le: int = 100,
    wake_length: float = 25.0,
    te_cell_size: float = 0.01,
    le_length: float = 0.05,
) -> tuple[Path, Path]:
    """Generate a C-block structured mesh using gmsh transfinite meshing.

    Uses a 5-block C-grid topology following the standard approach
    (wuFoil / NASA CFL3D validation grids):

      - Inlet block: semicircle around the leading edge region
      - Upper block: upper surface (LE split → TE) to farfield
      - Lower block: lower surface (LE split → TE) to farfield
      - Upper wake: TE to downstream exit (upper half)
      - Lower wake: TE to downstream exit (lower half)

    The airfoil is split into 3 segments: upper aft, leading edge, lower aft.
    The LE gets a separate block with "Bump" distribution for good LE resolution.
    The semicircle is centered at the LE split point.

    Requires ``pip install gmsh``.

    Returns (ugrid_path, mapbc_path).
    """
    import gmsh

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    first_cell = estimate_first_cell_height(Re, y_plus=y_plus)
    pts = np.array(coords, dtype=np.float64)

    # Split airfoil into 3 segments: upper aft, LE, lower aft
    # LE region = points with x < le_length
    le_idx = int(np.argmin(pts[:, 0]))

    # Find the split points where x crosses le_length
    # Upper: going from TE (x=1) toward LE (x=0), find where x < le_length
    upper_all = pts[:le_idx + 1]  # TE → LE
    lower_all = pts[le_idx:]      # LE → TE

    # Split upper at le_length
    upper_aft = []
    upper_le = []
    for i, (x, y) in enumerate(upper_all):
        if x > le_length:
            upper_aft.append((x, y))
        else:
            upper_le.append((x, y))
    # Include the split point in both segments
    if upper_aft and upper_le:
        upper_aft.append(upper_le[0])

    # Split lower at le_length
    lower_le = []
    lower_aft = []
    for i, (x, y) in enumerate(lower_all):
        if x <= le_length:
            lower_le.append((x, y))
        else:
            lower_aft.append((x, y))
    if lower_le and lower_aft:
        lower_aft.insert(0, lower_le[-1])

    R = farfield_radius
    center = (le_length, 0.0)

    # Growth factors (matching wuFoil's formulas)
    bl_growth = 1.0 / np.exp(np.log(first_cell) / (n_normal - 1))
    te_growth = (te_cell_size / 0.1) ** (1.0 / max(n_airfoil - 1, 1))
    wake_growth = 1.0 / np.exp(np.log(te_cell_size) / (n_wake - 1))

    # Detect open vs closed TE
    te_upper = upper_aft[0]
    te_lower = lower_aft[-1]
    te_gap = np.sqrt((te_upper[0] - te_lower[0])**2 + (te_upper[1] - te_lower[1])**2)
    open_te = te_gap > 1e-8

    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    gmsh.model.add("airfoil_cblock")
    geo = gmsh.model.geo

    # === AIRFOIL POINTS & SPLINES ===
    # For open-TE airfoils: both splines pass through their true TE coordinates
    # then converge to a shared midpoint. This preserves the airfoil shape
    # while giving a clean single-point C-grid topology.

    if open_te:
        te_mid = ((te_upper[0] + te_lower[0]) / 2.0,
                  (te_upper[1] + te_lower[1]) / 2.0)
        p_te_mid = geo.addPoint(te_mid[0], te_mid[1], 0)

        # Upper: TE_mid → TE_upper → ... → LE_split
        af_upper_pts = [p_te_mid]
        af_upper_pts.append(geo.addPoint(te_upper[0], te_upper[1], 0))
        for x, y in upper_aft[1:]:  # skip first (TE) since we added it explicitly
            af_upper_pts.append(geo.addPoint(x, y, 0))

        # LE region
        af_le_pts = [geo.addPoint(x, y, 0) for x, y in upper_le]
        for x, y in lower_le[1:]:
            af_le_pts.append(geo.addPoint(x, y, 0))

        # Lower: LE_split → ... → TE_lower → TE_mid
        af_lower_pts = []
        for x, y in lower_aft[:-1]:
            af_lower_pts.append(geo.addPoint(x, y, 0))
        af_lower_pts.append(geo.addPoint(te_lower[0], te_lower[1], 0))
        af_lower_pts.append(p_te_mid)  # shared TE midpoint
    else:
        # Closed TE: standard wuFoil approach
        af_upper_pts = [geo.addPoint(x, y, 0) for x, y in upper_aft]
        af_le_pts = [geo.addPoint(x, y, 0) for x, y in upper_le]
        for x, y in lower_le[1:]:
            af_le_pts.append(geo.addPoint(x, y, 0))
        af_lower_pts = [geo.addPoint(x, y, 0) for x, y in lower_aft[:-1]]
        af_lower_pts.append(af_upper_pts[0])  # share TE point
        p_te_mid = af_upper_pts[0]

    # Key shared points
    p_te = p_te_mid if open_te else af_upper_pts[0]
    p_top = af_upper_pts[-1]         # upper LE split
    p_bottom = af_lower_pts[0]       # lower LE split
    af_le_pts[0] = p_top             # share the LE split point
    af_le_pts[-1] = p_bottom         # share the LE split point

    c_upper = geo.addBSpline(af_upper_pts)
    c_le = geo.addBSpline(af_le_pts)
    c_lower = geo.addBSpline(af_lower_pts)

    # === FARFIELD & WAKE (unified 5-block topology) ===
    # Both open and closed TE now share a single p_te point.
    te_x = upper_aft[0][0]
    p_center = geo.addPoint(center[0], center[1], 0)
    p_inlet_top = geo.addPoint(center[0], R, 0)
    p_inlet_bottom = geo.addPoint(center[0], -R, 0)
    p_top_te = geo.addPoint(te_x, R, 0)
    p_bottom_te = geo.addPoint(te_x, -R, 0)

    wake_x = te_x + wake_length
    p_wake_top = geo.addPoint(wake_x, R, 0)
    p_wake_bottom = geo.addPoint(wake_x, -R, 0)
    p_wake_center = geo.addPoint(wake_x, 0.0, 0)

    c_inlet = geo.addCircleArc(p_inlet_top, p_center, p_inlet_bottom)
    c_top_to_inlet = geo.addLine(p_top, p_inlet_top)
    c_bottom_to_inlet = geo.addLine(p_inlet_bottom, p_bottom)
    c_te_to_top = geo.addLine(p_top_te, p_te)
    c_te_to_bottom = geo.addLine(p_te, p_bottom_te)
    c_top_line = geo.addLine(p_top_te, p_inlet_top)
    c_bottom_line = geo.addLine(p_inlet_bottom, p_bottom_te)
    c_wake_center = geo.addLine(p_te, p_wake_center)
    c_outlet_top = geo.addLine(p_wake_center, p_wake_top)
    c_outlet_bottom = geo.addLine(p_wake_bottom, p_wake_center)
    c_wake_top_line = geo.addLine(p_top_te, p_wake_top)
    c_wake_bottom_line = geo.addLine(p_bottom_te, p_wake_bottom)

    # 5 blocks
    inlet_loop = geo.addCurveLoop([c_inlet, c_bottom_to_inlet, -c_le, c_top_to_inlet])
    s_inlet = geo.addPlaneSurface([inlet_loop])
    top_loop = geo.addCurveLoop([-c_upper, -c_top_to_inlet, c_top_line, -c_te_to_top])
    s_top = geo.addPlaneSurface([top_loop])
    bottom_loop = geo.addCurveLoop([-c_bottom_to_inlet, -c_lower, -c_te_to_bottom, c_bottom_line])
    s_bottom = geo.addPlaneSurface([bottom_loop])
    wake_top_loop = geo.addCurveLoop([c_te_to_top, c_wake_center, c_outlet_top, -c_wake_top_line])
    s_wake_top = geo.addPlaneSurface([wake_top_loop])
    wake_bottom_loop = geo.addCurveLoop([-c_wake_center, c_outlet_bottom, c_wake_bottom_line, c_te_to_bottom])
    s_wake_bottom = geo.addPlaneSurface([wake_bottom_loop])

    all_surfaces = [s_inlet, s_top, s_bottom, s_wake_top, s_wake_bottom]
    geo.synchronize()

    mesh = gmsh.model.mesh
    mesh.setTransfiniteCurve(c_inlet, n_le, "Bump", -0.1)
    mesh.setTransfiniteCurve(c_le, n_le)
    mesh.setTransfiniteCurve(c_top_to_inlet, n_normal, "Progression", bl_growth)
    mesh.setTransfiniteCurve(c_bottom_to_inlet, n_normal, "Progression", -bl_growth)
    mesh.setTransfiniteCurve(c_te_to_top, n_normal, "Progression", -bl_growth)
    mesh.setTransfiniteCurve(c_upper, n_airfoil, "Progression", -te_growth)
    mesh.setTransfiniteCurve(c_top_line, n_airfoil, "Progression", -te_growth)
    mesh.setTransfiniteCurve(c_te_to_bottom, n_normal, "Progression", bl_growth)
    mesh.setTransfiniteCurve(c_lower, n_airfoil, "Progression", te_growth)
    mesh.setTransfiniteCurve(c_bottom_line, n_airfoil, "Progression", -te_growth)
    mesh.setTransfiniteCurve(c_wake_top_line, n_wake, "Progression", -wake_growth)
    mesh.setTransfiniteCurve(c_wake_center, n_wake, "Progression", wake_growth)
    mesh.setTransfiniteCurve(c_outlet_top, n_normal, "Progression", bl_growth)
    mesh.setTransfiniteCurve(c_wake_bottom_line, n_wake, "Progression", wake_growth)
    mesh.setTransfiniteCurve(c_outlet_bottom, n_normal, "Progression", -bl_growth)

    for s in all_surfaces:
        mesh.setTransfiniteSurface(s)
        mesh.setRecombine(2, s)

    gmsh.model.mesh.generate(2)

    wall_curves = [c_upper, c_lower, c_le]
    ff_curves = [c_inlet, c_top_line, c_bottom_line,
                 c_wake_top_line, c_wake_bottom_line,
                 c_outlet_top, c_outlet_bottom]

    gmsh.model.addPhysicalGroup(1, wall_curves, tag=1, name="wall")
    gmsh.model.addPhysicalGroup(1, ff_curves, tag=2, name="farfield")
    gmsh.model.addPhysicalGroup(2, all_surfaces, tag=10, name="fluid")

    # === EXTRACT MESH ===
    node_tags, node_coords, _ = gmsh.model.mesh.getNodes()
    n_nodes_2d = len(node_tags)
    coords_2d = node_coords.reshape(-1, 3)[:, :2]
    tag_to_idx = {int(t): i for i, t in enumerate(node_tags)}

    # Quads only (all recombined)
    elem_types, elem_tags_list, elem_nodes_list = gmsh.model.mesh.getElements(dim=2)
    quad_conn = None
    for et, _, enodes in zip(elem_types, elem_tags_list, elem_nodes_list):
        props = gmsh.model.mesh.getElementProperties(et)
        if "Quadrangle" in props[0] or "Quad" in props[0]:
            quad_conn = enodes.reshape(-1, 4)

    if quad_conn is None:
        gmsh.finalize()
        raise RuntimeError("gmsh did not produce quad elements — transfinite meshing failed")

    # Wall edges — use the wall_curves/ff_curves set by topology branch above
    wall_edges = []
    for curve in wall_curves:
        _, _, enodes = gmsh.model.mesh.getElements(dim=1, tag=curve)
        for en in enodes:
            wall_edges.append(en.reshape(-1, 2))
    wall_edges = np.vstack(wall_edges)

    # Farfield + outlet edges
    ff_edges = []
    for curve in ff_curves:
        _, _, enodes = gmsh.model.mesh.getElements(dim=1, tag=curve)
        for en in enodes:
            ff_edges.append(en.reshape(-1, 2))
    ff_edges = np.vstack(ff_edges)

    gmsh.finalize()

    # === EXTRUDE TO 3D (airfoil in x-z plane, span in y) ===
    nodes_y0 = np.column_stack([coords_2d[:, 0], np.zeros(n_nodes_2d), coords_2d[:, 1]])
    nodes_y1 = np.column_stack([coords_2d[:, 0], np.full(n_nodes_2d, span), coords_2d[:, 1]])
    all_nodes = np.vstack([nodes_y0, nodes_y1])

    # Hex elements from quads
    hexes = []
    for q in quad_conn:
        idx = [tag_to_idx[int(t)] for t in q]
        n0, n1, n2, n3 = idx
        p0, p1, p3, p4 = nodes_y0[n0], nodes_y0[n1], nodes_y0[n3], nodes_y1[n0]
        vol = np.dot(p1 - p0, np.cross(p3 - p0, p4 - p0))
        if vol > 0:
            hexes.append([n0+1, n1+1, n2+1, n3+1,
                         n0+n_nodes_2d+1, n1+n_nodes_2d+1, n2+n_nodes_2d+1, n3+n_nodes_2d+1])
        else:
            hexes.append([n0+1, n3+1, n2+1, n1+1,
                         n0+n_nodes_2d+1, n3+n_nodes_2d+1, n2+n_nodes_2d+1, n1+n_nodes_2d+1])
    hexes = np.array(hexes, dtype=np.int32)

    # Boundary quads from edges
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
    sym_y0_quads = [[tag_to_idx[int(t)]+1 for t in q] for q in quad_conn]
    sym_y1_quads = [[tag_to_idx[int(t)]+n_nodes_2d+1 for t in q] for q in quad_conn]

    # === WRITE UGRID ===
    all_bnd_quads, all_quad_tags = [], []
    for quads, tag in [(wall_quads.tolist(), 1), (ff_quads.tolist(), 2),
                       (sym_y0_quads, 3), (sym_y1_quads, 4)]:
        for q in quads:
            all_bnd_quads.append(q)
            all_quad_tags.append(tag)

    # === WATERTIGHT CHECK ===
    # Every hex face that's not a boundary quad must be shared by exactly 2 hexes.
    # Boundary faces must be shared by exactly 1 hex.
    from collections import Counter
    def _face_key(nodes):
        return tuple(sorted(nodes))

    all_hex_faces = Counter()
    for h in hexes:
        # 6 faces of a hex (each is 4 nodes)
        faces = [
            (h[0], h[1], h[2], h[3]),  # bottom
            (h[4], h[5], h[6], h[7]),  # top
            (h[0], h[1], h[5], h[4]),  # front
            (h[1], h[2], h[6], h[5]),  # right
            (h[2], h[3], h[7], h[6]),  # back
            (h[3], h[0], h[4], h[7]),  # left
        ]
        for f in faces:
            all_hex_faces[_face_key(f)] += 1

    bnd_face_keys = set()
    for q in all_bnd_quads:
        bnd_face_keys.add(_face_key(q))

    # Interior faces should appear exactly 2 times
    # Boundary faces should appear exactly 1 time
    n_interior_bad = 0
    n_boundary_bad = 0
    n_unmatched_bnd = 0
    for fk, count in all_hex_faces.items():
        if fk in bnd_face_keys:
            if count != 1:
                n_boundary_bad += 1
        else:
            if count != 2:
                n_interior_bad += 1

    for fk in bnd_face_keys:
        if fk not in all_hex_faces:
            n_unmatched_bnd += 1

    if n_interior_bad > 0 or n_boundary_bad > 0 or n_unmatched_bnd > 0:
        raise RuntimeError(
            f"Mesh is NOT watertight: {n_interior_bad} bad interior faces, "
            f"{n_boundary_bad} bad boundary faces, {n_unmatched_bnd} unmatched boundary faces"
        )

    ugrid_path = output_dir / f"{mesh_name}.b8.ugrid"
    mapbc_path = output_dir / f"{mesh_name}.mapbc"

    with open(ugrid_path, "wb") as f:
        f.write(struct.pack(">7i", len(all_nodes), 0, len(all_bnd_quads), 0, 0, 0, len(hexes)))
        for node in all_nodes:
            f.write(struct.pack(">3d", *node))
        for quad in all_bnd_quads:
            f.write(struct.pack(">4i", *quad))
        for tag in all_quad_tags:
            f.write(struct.pack(">i", tag))
        for hex_elem in hexes:
            f.write(struct.pack(">8i", *hex_elem))

    mapbc_path.write_text("4\n1 4000 wall\n2 3000 farfield\n3 5000 symmetry_y0\n4 5000 symmetry_y1\n")

    return ugrid_path, mapbc_path
