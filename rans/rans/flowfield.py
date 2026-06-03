"""Export the RANS center-span slice as a flow mesh for the web UI's LIC renderer.

The solver writes ``slice_centerSpan_proc0.vtu`` — the 2D field on the mid-span plane
(mesh frame: x=streamwise, y=span, z=lift). We export it as a triangle mesh with
per-vertex velocity + a Mach scalar, projected to UI coordinates (x=mesh-x, y=mesh-z;
u=vx, v=vz). This is exactly the input the ported flow360 two-pass LIC shader needs
(pass 1 renders the mesh with the per-vertex flow vector; see flexfoil-ui LIC layer).
"""
from __future__ import annotations

from pathlib import Path

import vtk
from vtk.util.numpy_support import vtk_to_numpy
import numpy as np


# Clip the exported slice to the airfoil view (+ margin) so we don't ship the whole
# far-field disk. Slightly larger than the UI view so LIC streaks have context.
CLIP_BOUNDS = (-0.4, 1.6, -0.7, 0.5)


def extract_flow_mesh(workdir: str | Path, clip: tuple[float, float, float, float] = CLIP_BOUNDS,
                      step: int | None = None) -> dict:
    """Return the center-span slice as ``{points, tris, u, v, mach, machRange, bounds}``
    with flat arrays (points/vel interleaved x,y; tris flat i,j,k), projected to UI
    coords and clipped to ``clip``. Triangle mesh, per-vertex velocity — ready for the
    LIC pass-1 render.

    ``step`` selects an unsteady physical step's slice (1-based, e.g. an α-sweep point ⇒
    ``slice_centerSpan_time_{step}_proc0.vtu``); ``None`` reads the final/steady slice."""
    workdir = Path(workdir)
    name = "slice_centerSpan_proc0.vtu" if step is None else f"slice_centerSpan_time_{step}_proc0.vtu"
    r = vtk.vtkXMLUnstructuredGridReader()
    r.SetFileName(str(workdir / name))
    r.Update()
    # surface + triangulate (robust if any non-triangle cells appear)
    geo = vtk.vtkGeometryFilter()
    geo.SetInputConnection(r.GetOutputPort())
    tf = vtk.vtkTriangleFilter()
    tf.SetInputConnection(geo.GetOutputPort())
    tf.Update()
    g = tf.GetOutput()

    pts = vtk_to_numpy(g.GetPoints().GetData())               # (N,3): x, y_span, z
    vel = vtk_to_numpy(g.GetPointData().GetArray("velocity")) # (N,3): vx, vy, vz
    mach = vtk_to_numpy(g.GetPointData().GetArray("Mach"))     # (N,)
    tris = vtk_to_numpy(g.GetPolys().GetConnectivityArray()).reshape(-1, 3)

    # UI projection: (x, z) for position, (vx, vz) for velocity
    px, py = pts[:, 0], pts[:, 2]
    u, v = vel[:, 0], vel[:, 2]

    # Clip to the view region: keep triangles whose centroid is inside `clip`,
    # then drop unused points and reindex.
    xmin, xmax, ymin, ymax = clip
    cx = px[tris].mean(axis=1)
    cy = py[tris].mean(axis=1)
    keep = (cx >= xmin) & (cx <= xmax) & (cy >= ymin) & (cy <= ymax)
    tris = tris[keep]
    used = np.unique(tris)
    remap = np.full(len(px), -1, dtype=np.int64)
    remap[used] = np.arange(len(used))
    tris = remap[tris]
    px, py, u, v, mach = px[used], py[used], u[used], v[used], mach[used]

    points = np.column_stack([px, py]).round(5).ravel().tolist()
    uv = np.column_stack([u, v]).round(5).ravel().tolist()
    return {
        "points": points,                                     # flat [x0,y0,x1,y1,...]
        "tris": tris.astype(np.int32).ravel().tolist(),       # flat [i0,j0,k0,...]
        "vel": uv,                                            # flat [u0,v0,u1,v1,...]
        "mach": np.round(mach, 5).tolist(),
        "machRange": [float(mach.min()), float(mach.max())],
        "bounds": [float(px.min()), float(px.max()), float(py.min()), float(py.max())],
        "nPoints": int(len(px)), "nTris": int(len(tris)),
    }
