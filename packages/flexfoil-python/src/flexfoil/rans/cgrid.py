"""Structured C-grid mesh generator for 2D airfoils.

Generates a single-block structured quad mesh using:
1. Normal extrusion from the airfoil surface with geometric BL stretching
2. Transfinite interpolation (TFI) to blend near-wall grid with farfield
3. Optional Laplace smoothing for improved orthogonality

Handles both open (blunt) and closed (sharp) trailing edges.
No external meshing dependencies — pure numpy.

References:
- Thompson, Thames, Mastin (1974), J. Comp. Physics 15, 299-319
- Gordon & Hall (1973), Transfinite Interpolation (TFI)
- Construct2D (Fortran, GPL3) for algorithm reference

Grid topology (C-grid):
    Inner boundary (j=0): wake_lower → lower_surface → upper_surface → wake_upper
    Outer boundary (j=N): straight_lower → semicircle → straight_upper
    i runs along the C-shape, j runs wall-normal (surface to farfield).
"""

from __future__ import annotations

import numpy as np
from scipy.interpolate import interp1d


def generate_cgrid(
    coords: list[tuple[float, float]],
    *,
    n_normal: int = 100,
    n_wake: int = 60,
    first_cell_height: float | None = None,
    Re: float = 1e6,
    growth_rate: float = 1.08,
    farfield_radius: float = 50.0,
    wake_length: float = 10.0,
    y_plus: float = 1.0,
    smooth_iterations: int = 0,
) -> dict:
    """Generate a structured C-grid around a 2D airfoil.

    Parameters
    ----------
    coords : list of (x, y)
        Airfoil coordinates in Selig ordering (upper TE → LE → lower TE).
    n_normal : int
        Number of cells in the wall-normal direction.
    n_wake : int
        Number of cells along each wake segment.
    first_cell_height : float or None
        First cell height off the wall. Auto-estimated from Re if None.
    Re : float
        Reynolds number (for auto first_cell_height).
    growth_rate : float
        Geometric growth rate for BL layers.
    farfield_radius : float
        Farfield distance in chord lengths.
    wake_length : float
        Wake extension downstream of TE in chord lengths.
    y_plus : float
        Target y+ for first cell estimation.
    smooth_iterations : int
        Number of Laplace smoothing iterations (0 = skip).

    Returns
    -------
    dict with:
        'X': (n_i, n_j) array of x-coordinates
        'Y': (n_i, n_j) array of y-coordinates
        'n_i': number of points along the C-curve
        'n_j': number of points in the wall-normal direction
        'n_wake': number of wake points on each side
        'n_surface': number of airfoil surface points
    """
    surface = np.array(coords, dtype=np.float64)
    n_pts = len(surface)

    if n_pts < 20:
        raise ValueError(f"Need at least 20 surface points, got {n_pts}")

    # Find LE and TE
    le_idx = int(np.argmin(surface[:, 0]))
    upper = surface[:le_idx + 1]  # TE_upper → LE
    lower = surface[le_idx:]      # LE → TE_lower

    te_upper = upper[0].copy()
    te_lower = lower[-1].copy()

    # Estimate first cell height if needed
    if first_cell_height is None:
        cf = 0.058 * Re ** (-0.2)
        u_tau_norm = np.sqrt(cf / 2.0)
        first_cell_height = y_plus / (Re * u_tau_norm)

    # -------------------------------------------------------------------------
    # Step 1: Build the inner boundary (C-curve)
    # Order: wake_lower(reversed) → lower_surface(reversed) → upper_surface → wake_upper
    # This traces the C from bottom-right, around the LE, to top-right
    # -------------------------------------------------------------------------

    # Wake lines extend downstream from TE
    wake_spacing = np.zeros(n_wake)
    wake_first = 0.01  # first wake cell matches TE spacing
    for i in range(n_wake):
        wake_spacing[i] = wake_first * (1.2 ** i)
    wake_x_offsets = np.cumsum(wake_spacing)
    if wake_x_offsets[-1] > 0:
        wake_x_offsets *= wake_length / wake_x_offsets[-1]

    # Upper wake: extends from te_upper
    wake_upper = np.column_stack([
        te_upper[0] + wake_x_offsets,
        np.full(n_wake, te_upper[1]),
    ])

    # Lower wake: extends from te_lower
    wake_lower = np.column_stack([
        te_lower[0] + wake_x_offsets,
        np.full(n_wake, te_lower[1]),
    ])

    # Build the C-curve: wake_lower(reversed) → lower(reversed) → upper → wake_upper
    lower_reversed = lower[::-1]  # TE_lower → LE
    wake_lower_reversed = wake_lower[::-1]  # far wake → TE_lower

    inner = np.vstack([
        wake_lower_reversed,       # far wake bottom → TE_lower
        lower_reversed[1:],        # TE_lower → LE (skip duplicate at TE_lower)
        upper[1:],                 # LE → TE_upper (skip duplicate at LE)
        wake_upper,                # TE_upper → far wake top
    ])
    n_i = len(inner)

    # -------------------------------------------------------------------------
    # Step 2: Compute outward normals along the inner boundary
    # -------------------------------------------------------------------------
    normals = _compute_normals(inner)

    # Force wake normals to point straight up/down (no lateral component)
    for i in range(n_wake):
        normals[i] = [0.0, -1.0]        # lower wake: point down
    for i in range(n_i - n_wake, n_i):
        normals[i] = [0.0, 1.0]         # upper wake: point up

    # -------------------------------------------------------------------------
    # Step 3: Build the outer boundary (C-shape at farfield)
    # -------------------------------------------------------------------------
    R = farfield_radius
    outer = np.zeros_like(inner)

    # The outer boundary mirrors the inner boundary's topology:
    # - Wake sections: straight lines at ±R from the wake
    # - Airfoil section: semicircle at radius R centered at (0.5, 0)
    cx = 0.5  # circle center x (mid-chord)

    # The outer boundary must follow the SAME direction as the inner boundary.
    # Inner goes: lower wake far → TE_lower → LE → TE_upper → upper wake far
    # So y goes: te_lower.y → negative(lower surface) → LE → positive(upper surface) → te_upper.y
    # The outer boundary should trace the same path at radius R:
    # y = -R (below TE) → semicircle around front → y = +R (above TE)

    for i in range(n_i):
        if i < n_wake:
            # Lower wake: straight line at y = te_lower.y, extending to x = farfield
            # But we want the outermost point to be at the farfield radius below
            # Project the inner wake point outward
            outer[i] = [inner[i, 0], -R]
        elif i >= n_i - n_wake:
            # Upper wake: straight line at y = te_upper.y, extending to farfield
            outer[i] = [inner[i, 0], R]
        else:
            # Airfoil portion: map to semicircle
            # The inner boundary goes from lower TE (y<0) through LE (x_min)
            # to upper TE (y>0). Map this to a semicircle going the same way:
            # bottom (y=-R) → left (x=-R) → top (y=+R)
            i_airfoil = i - n_wake
            n_airfoil = n_i - 2 * n_wake

            # Use the inner boundary's angular position relative to the center
            # to determine the outer boundary position on the semicircle.
            # This ensures the mapping is monotonic and follows the same direction.
            pt = inner[i]
            angle_inner = np.arctan2(pt[1], pt[0] - cx)

            # Scale to the semicircle: maintain the same angle, just at radius R
            outer[i] = [cx + R * np.cos(angle_inner), R * np.sin(angle_inner)]

    # -------------------------------------------------------------------------
    # Step 4: Build radial stretching (geometric near wall, blend to farfield)
    # -------------------------------------------------------------------------
    n_j = n_normal + 1

    # Geometric heights
    heights = np.zeros(n_j)
    for j in range(1, n_j):
        heights[j] = heights[j - 1] + first_cell_height * growth_rate ** (j - 1)

    # Scale to reach farfield
    max_h = heights[-1]
    target = R
    if max_h < target:
        # Blend: keep geometric near wall, stretch outer layers
        t = np.linspace(0, 1, n_j)
        blend = t ** 1.5  # quadratic blend favoring near-wall
        heights = heights * (1 - blend) + heights * (target / max_h) * blend
    elif max_h > target:
        heights *= target / max_h

    # Normalize to [0, 1]
    s = heights / heights[-1]

    # -------------------------------------------------------------------------
    # Step 5: Fill the grid using normal extrusion + TFI blending
    # -------------------------------------------------------------------------
    # Phase 1 (j < j_blend): Pure normal extrusion from the surface.
    #   Guarantees orthogonality at the wall and avoids cell crossing.
    # Phase 2 (j >= j_blend): Blend the extruded grid toward the outer boundary
    #   using a smooth transition. The blend factor goes from 0 (pure extrusion)
    #   to 1 (pure outer boundary) over the remaining layers.

    X = np.zeros((n_i, n_j))
    Y = np.zeros((n_i, n_j))

    # Determine blend start: where extrusion height reaches ~50% of farfield
    blend_fraction = 0.5
    j_blend = 0
    for j in range(n_j):
        if heights[j] > blend_fraction * R:
            j_blend = j
            break
    if j_blend == 0:
        j_blend = int(0.7 * n_j)  # fallback

    # Phase 1: Normal extrusion
    for j in range(j_blend + 1):
        X[:, j] = inner[:, 0] + heights[j] * normals[:, 0]
        Y[:, j] = inner[:, 1] + heights[j] * normals[:, 1]

    # Phase 2: Blend from extruded position to outer boundary
    extruded_at_blend = np.column_stack([X[:, j_blend], Y[:, j_blend]])
    for j in range(j_blend + 1, n_j):
        # Blend factor: 0 at j_blend, 1 at j=n_j-1
        t_blend = (j - j_blend) / (n_j - 1 - j_blend)
        # Smooth blend (cubic ease)
        t_smooth = t_blend * t_blend * (3 - 2 * t_blend)
        X[:, j] = (1 - t_smooth) * extruded_at_blend[:, 0] + t_smooth * outer[:, 0]
        Y[:, j] = (1 - t_smooth) * extruded_at_blend[:, 1] + t_smooth * outer[:, 1]

    # -------------------------------------------------------------------------
    # Step 6: Optional Laplace smoothing
    # -------------------------------------------------------------------------
    if smooth_iterations > 0:
        X, Y = _laplace_smooth(X, Y, n_wake, smooth_iterations)

    return {
        'X': X,
        'Y': Y,
        'n_i': n_i,
        'n_j': n_j,
        'n_wake': n_wake,
        'n_surface': n_pts,
        'first_cell_height': first_cell_height,
    }


def _compute_normals(curve: np.ndarray) -> np.ndarray:
    """Compute unit outward normals along an open curve."""
    n = len(curve)
    tangents = np.zeros_like(curve)

    # Central differences for interior, one-sided at endpoints
    tangents[1:-1] = curve[2:] - curve[:-2]
    tangents[0] = curve[1] - curve[0]
    tangents[-1] = curve[-1] - curve[-2]

    lengths = np.linalg.norm(tangents, axis=1, keepdims=True)
    lengths = np.maximum(lengths, 1e-14)
    tangents = tangents / lengths

    # Rotate 90° to get normal
    normals = np.column_stack([tangents[:, 1], -tangents[:, 0]])

    # Check orientation: normals should point away from centroid
    centroid = curve.mean(axis=0)
    outward = curve - centroid
    dots = np.sum(normals * outward, axis=1)
    if np.sum(dots < 0) > np.sum(dots > 0):
        normals = -normals

    return normals


def _laplace_smooth(X, Y, n_wake, iterations, omega=1.0):
    """Laplace smoothing of interior grid points.

    Fixes all boundary points (j=0, j=N-1, i=0, i=M-1) and smooths interior.
    """
    n_i, n_j = X.shape

    for _it in range(iterations):
        for j in range(1, n_j - 1):
            for i in range(1, n_i - 1):
                X[i, j] = (1 - omega) * X[i, j] + omega * 0.25 * (
                    X[i+1, j] + X[i-1, j] + X[i, j+1] + X[i, j-1]
                )
                Y[i, j] = (1 - omega) * Y[i, j] + omega * 0.25 * (
                    Y[i+1, j] + Y[i-1, j] + Y[i, j+1] + Y[i, j-1]
                )

    return X, Y


def cgrid_to_3d(
    grid: dict,
    *,
    span: float = 0.01,
) -> dict:
    """Extrude a 2D C-grid one cell deep in the z-direction.

    Returns a dict with nodes, hex connectivity, and boundary faces,
    ready for write_ugrid().
    """
    X, Y = grid['X'], grid['Y']
    n_i, n_j = grid['n_i'], grid['n_j']
    n_wake = grid['n_wake']

    # Create 3D nodes: z=0 and z=span
    n_2d = n_i * n_j
    nodes_z0 = np.column_stack([X.ravel(), Y.ravel(), np.zeros(n_2d)])
    nodes_z1 = np.column_stack([X.ravel(), Y.ravel(), np.full(n_2d, span)])
    nodes = np.vstack([nodes_z0, nodes_z1])

    def idx(i, j, layer=0):
        """Node index in the flat array (1-based for UGRID)."""
        return layer * n_2d + j * n_i + i + 1

    # Hex elements: one per quad cell
    n_cells_i = n_i - 1
    n_cells_j = n_j - 1
    hexes = []
    for j in range(n_cells_j):
        for i in range(n_cells_i):
            # Bottom face (z=0): CCW when viewed from below
            n0 = idx(i, j, 0)
            n1 = idx(i+1, j, 0)
            n2 = idx(i+1, j+1, 0)
            n3 = idx(i, j+1, 0)
            # Top face (z=span)
            n4 = idx(i, j, 1)
            n5 = idx(i+1, j, 1)
            n6 = idx(i+1, j+1, 1)
            n7 = idx(i, j+1, 1)
            # Check volume and fix orientation if needed
            p0, p1, p3, p4 = nodes[n0-1], nodes[n1-1], nodes[n3-1], nodes[n4-1]
            vol = np.dot(p1 - p0, np.cross(p3 - p0, p4 - p0))
            if vol >= 0:
                hexes.append([n0, n1, n2, n3, n4, n5, n6, n7])
            else:
                hexes.append([n0, n3, n2, n1, n4, n7, n6, n5])

    hexes = np.array(hexes, dtype=np.int32)

    # Boundary faces
    # Wall: j=0, excluding wake segments
    wall_start = n_wake  # first airfoil cell index
    wall_end = n_i - n_wake - 1  # last airfoil cell index
    wall_quads = []
    for i in range(wall_start, wall_end):
        wall_quads.append([idx(i, 0, 0), idx(i, 0, 1), idx(i+1, 0, 1), idx(i+1, 0, 0)])

    # Farfield: j=n_j-1 (outermost ring)
    farfield_quads = []
    for i in range(n_cells_i):
        j = n_cells_j
        farfield_quads.append([idx(i, j, 0), idx(i+1, j, 0), idx(i+1, j, 1), idx(i, j, 1)])

    # Wake: i=0 and i=n_i-1 boundaries (where the C-grid opens)
    wake_quads = []
    for j in range(n_cells_j):
        # i=0 (lower wake exit)
        wake_quads.append([idx(0, j, 0), idx(0, j, 1), idx(0, j+1, 1), idx(0, j+1, 0)])
        # i=n_i-1 (upper wake exit)
        wake_quads.append([idx(n_i-1, j, 0), idx(n_i-1, j+1, 0), idx(n_i-1, j+1, 1), idx(n_i-1, j, 1)])

    # Symmetry z=0 and z=span
    sym_z0_quads = []
    sym_z1_quads = []
    for j in range(n_cells_j):
        for i in range(n_cells_i):
            sym_z0_quads.append([idx(i, j, 0), idx(i, j+1, 0), idx(i+1, j+1, 0), idx(i+1, j, 0)])
            sym_z1_quads.append([idx(i, j, 1), idx(i+1, j, 1), idx(i+1, j+1, 1), idx(i, j+1, 1)])

    boundary_quads = {
        'wall': np.array(wall_quads, dtype=np.int32),
        'farfield': np.array(farfield_quads + wake_quads, dtype=np.int32),
        'symmetry_y0': np.array(sym_z0_quads, dtype=np.int32),
        'symmetry_y1': np.array(sym_z1_quads, dtype=np.int32),
    }

    boundary_ids = {'wall': 1, 'farfield': 2, 'symmetry_y0': 3, 'symmetry_y1': 4}

    return {
        'nodes': nodes,
        'hexes': hexes,
        'boundary_quads': boundary_quads,
        'boundary_ids': boundary_ids,
        'n_nodes': len(nodes),
        'n_hexes': len(hexes),
    }
