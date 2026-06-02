# flexfoil/rans — quasi-2D multi-element RANS pipeline

Turns a **contours JSON** (deployed element polylines + flow conditions) into a
Flow360 RANS solution: CL/CD plus ParaView fields. This is the viscous engine
that slots in alongside the inviscid coupled solver behind the high-lift design
tool.

```
contours JSON ──▶ Spalding metric ──▶ in-house 2D mesh ──▶ extrude (quasi-2D)
            ──▶ CGNS ──▶ Flow360 case ──▶ GPU solve ──▶ CL/CD + Cp/volume.pvtu
```

## Layout

```
rans/
├── rans/                 python package
│   ├── config.py         CaseConfig schema (the contours JSON) + Spalding h0
│   ├── contours.py       config → contours.txt (mesher input)
│   ├── mesh.py           mesh2d → extrude (x=streamwise, z=lift, y=span) → CGNS
│   ├── case.py           simulation.json (flow360 API) + preprocessing → Flow360.json
│   ├── solve.py          run solver (+postprocessor IPC) + extract forces
│   ├── env.py            locate the compute install; resolve tools
│   └── pipeline.py       run(config, outdir, solve=…)  ← entry point
├── mesher/               mesh2d.cpp (in-house CavityBasedMesher + Spalding) + build.sh
├── bin/mesh2d            prebuilt mesher binary
├── examples/
│   └── highlift_deploy012.json   the 3-element high-lift config (deploy 0.12)
└── run_case.py           CLI
```

## Usage

Run with the **compute venv's python** (it has the version-matched `flow360` SDK):

```bash
PY=/home/qiqi/flexcompute/compute/.venv/bin/python

# Build mesh + Flow360 case (runs anywhere, incl. a sandbox):
$PY run_case.py examples/highlift_deploy012.json -o out/

# Build AND solve — run from a normal interactive shell (see caveat):
$PY run_case.py examples/highlift_deploy012.json -o out/ --solve
```

Or programmatically:

```python
from rans import run
summary = run("examples/highlift_deploy012.json", "out/", solve=True)
print(summary["forces"])   # {'CL':…, 'CD':…, 'L_over_D':…, 'paraview':{…}}
```

## Contours JSON

```jsonc
{
  "elements": [                       // deployed, global frame (buildConfiguration output)
    {"name": "main", "contour": [[x,y], …], "is_wall": true},
    {"name": "vane", "contour": [[x,y], …]},
    {"name": "flap", "contour": [[x,y], …]}
  ],
  "farfield": {"type": "circle", "center": [0.5,-0.1], "radius": 50, "n": 240},
  "flow":   {"reynolds": 1e7, "mach": 0.2, "alpha_deg": 0.0, "temperature": 288.15},
  "mesh":   {"span": 0.1, "nspan": 1, "yplus": 1.0,
             "growth": 1.2, "hwall": 0.004, "hmax": 3.0, "h0": 0.0},  // h0=0 → from Re
  "solver": {"max_steps": 5000}
}
```

Element names become the wall patch names (`fluid/<name>`). The reference area is
set to `chord·span` so CL/CD are span-normalized and span-independent (1 cell is
adequate — see the span study). Lift is on **z**, drag on **x**, span on **y**.

## Web-tool integration

The clean seam mirrors `flexfoil-ui/src/highlift/solve.ts`: the UI's
`buildConfiguration()` / `toPoints()` already produce the deployed contours that
go straight into `elements[].contour`. A backend serializes that geometry + the
chosen flow conditions to a contours JSON, calls `run(..., solve=True)`, and reads
`summary["forces"]` (CL/CD) and `summary["forces"]["paraview"]` (field files) —
the RANS analogue of `analyzeConfiguration()`.

## ⚠️ Solving needs a non-sandboxed shell

Meshing and case construction run anywhere. **`Flow360Solver` (GPU/MPI) must run
in a normal interactive shell** — inside a restricted agent sandbox it is killed
before it starts. Build the case in-session if you like, then run the solve (or the
whole thing with `--solve`) from a real shell.

## Configuration

- `FLOW360_COMPUTE_ROOT` (default `/home/qiqi/flexcompute/compute`) — the compute
  install (`install/release/bin` binaries, `.venv` Python tools).
- Rebuild the mesher only if the compute libraries change: `mesher/build.sh`.
