# flexfoil

Airfoil analysis in Python -- XFOIL-faithful viscous solver with a local web UI.

```python
import flexfoil

foil = flexfoil.naca("2412")
result = foil.solve(alpha=5.0, Re=1e6)
print(result)
# SolveResult(α=5.00°, Re=1e+06, CL=0.8094, CD=0.00775, CM=-0.0540, converged)
```

## Install

```bash
pip install flexfoil
```

Pre-built wheels are available for:
- macOS (Apple Silicon and Intel)
- Linux x86_64
- Windows x86_64

Plotly is included by default for interactive plots. Optional extras:

```bash
pip install "flexfoil[all]"        # server + matplotlib + pandas
pip install "flexfoil[matplotlib]" # matplotlib for polar.plot(backend="matplotlib")
pip install "flexfoil[dataframe]"  # pandas for polar.to_dataframe()
pip install "flexfoil[server]"     # starlette + uvicorn for flexfoil.serve()
```

## What is this?

flexfoil wraps the [RustFoil](https://foil.flexcompute.com/flexfoil/) solver
(an XFOIL-faithful Rust reimplementation) in native Python bindings via
[PyO3](https://pyo3.rs). Every solve is cached in a local SQLite database
(`~/.flexfoil/runs.db`). The same web UI from
[foil.flexcompute.com](https://foil.flexcompute.com/flexfoil/) can be launched
locally with `flexfoil.serve()`, reading from the same database. **No data
leaves your machine** unless you explicitly export it.

## Quick start

### Single-point solve

```python
import flexfoil

foil = flexfoil.naca("2412")
result = foil.solve(alpha=5.0, Re=1e6)

print(result.cl)          # 0.8094
print(result.cd)          # 0.00775
print(result.cm)          # -0.0540
print(result.converged)   # True
print(result.ld)          # 104.5
print(result.x_tr_upper)  # 0.315  (transition location, upper surface)
```

### Polar sweep

Polar sweeps are **parallelized by default** using all available CPU cores
(via Rust's rayon thread pool). A 41-point polar runs ~3x faster than sequential.

```python
polar = foil.polar(alpha=(-5, 15, 0.5), Re=1e6)
print(polar)
# PolarResult('NACA 2412', Re=1e+06, 40/41 converged)

# Interactive plotly figure (default): CL-α, CD-α, CL-CD, CM-α
polar.plot()

# Or use matplotlib
polar.plot(backend="matplotlib")

# Sequential mode (for debugging or progress output)
polar = foil.polar(alpha=(-5, 15, 0.5), Re=1e6, parallel=False)

# Export to pandas
df = polar.to_dataframe()
df.to_csv("polar.csv", index=False)
```

### Compare multiple airfoils

```python
import matplotlib.pyplot as plt
import flexfoil

for naca in ["0012", "2412", "4412"]:
    foil = flexfoil.naca(naca)
    polar = foil.polar(alpha=(-4, 14, 1.0), Re=1e6)
    plt.plot(polar.alpha, polar.cl, ".-", label=foil.name)

plt.xlabel("α (°)")
plt.ylabel("CL")
plt.legend()
plt.grid(True, alpha=0.3)
plt.show()
```

### Load a .dat file

```python
foil = flexfoil.load("e387.dat")   # Selig or Lednicer format
result = foil.solve(alpha=4.0, Re=2e5)
```

### Custom coordinates

```python
foil = flexfoil.from_coordinates(x_list, y_list, name="my shape")
result = foil.solve(alpha=3.0, Re=1e6)
```

### Inviscid analysis

```python
result = foil.solve(alpha=5.0, viscous=False)
# result.cd == 0.0 (no drag in potential flow)
```

## Launch the web UI

```python
flexfoil.serve()
# Opens the full flexfoil web app in your browser.
# Defaults to port 8420; automatically picks the next free port if busy.
```

Or from the command line:

```bash
flexfoil serve
flexfoil serve --port 9000 --no-browser
```

The web UI reads from the same local SQLite database as the Python API.
Runs you solve from the browser (via the built-in WASM solver) are also
written to the shared database, so `flexfoil.runs()` in Python sees them
and vice versa.

## Query the run database

Every solve (from Python or the web UI) is cached in `~/.flexfoil/runs.db`.

```python
# All runs as a pandas DataFrame (or list[dict] if pandas not installed)
df = flexfoil.runs()
print(f"{len(df)} runs cached")

# Filter by airfoil
naca_runs = df[df.airfoil_name == "NACA 2412"]
```

## CLI

```
flexfoil solve 2412 -a 5 -r 1e6    # quick solve from the terminal
flexfoil serve                       # launch the web UI
flexfoil info                        # show config and DB location
```

## API reference

### Top-level functions

| Function | Description |
| --- | --- |
| `flexfoil.naca(designation, n_panels=160)` | Create airfoil from NACA 4-digit string |
| `flexfoil.load(path, n_panels=160)` | Load from `.dat` file |
| `flexfoil.from_coordinates(x, y, name, n_panels)` | Create from raw x/y arrays |
| `flexfoil.runs()` | All cached runs (DataFrame or list) |
| `flexfoil.serve(port, host, open_browser)` | Launch local web UI + API server |

### `Airfoil`

| Property / Method | Description |
| --- | --- |
| `.name` | Airfoil name |
| `.n_panels` | Number of panel nodes |
| `.raw_coords` | Original coordinates `list[(x, y)]` |
| `.panel_coords` | Repaneled coordinates `list[(x, y)]` |
| `.hash` | SHA-256 hash of panel coords (cache key) |
| `.solve(alpha, Re, mach, ncrit, max_iter, viscous, store)` | Single-point analysis |
| `.polar(alpha, Re, mach, ncrit, max_iter, viscous, store, parallel)` | Sweep over alpha range (parallel by default) |

### `SolveResult`

| Field | Type | Description |
| --- | --- | --- |
| `.cl` | `float` | Lift coefficient |
| `.cd` | `float` | Drag coefficient |
| `.cm` | `float` | Moment coefficient (quarter-chord) |
| `.converged` | `bool` | Whether the Newton solve converged |
| `.iterations` | `int` | Newton iterations used |
| `.residual` | `float` | Final Newton residual |
| `.x_tr_upper` | `float` | Transition x/c, upper surface |
| `.x_tr_lower` | `float` | Transition x/c, lower surface |
| `.alpha` | `float` | Angle of attack (degrees) |
| `.reynolds` | `float` | Reynolds number |
| `.mach` | `float` | Mach number |
| `.ncrit` | `float` | e^N transition criterion |
| `.ld` | `float \| None` | Lift-to-drag ratio |
| `.success` | `bool` | Overall success flag |
| `.error` | `str \| None` | Error message if failed |

### `PolarResult`

| Property / Method | Description |
| --- | --- |
| `.alpha` | `list[float]` — angles (converged only) |
| `.cl` | `list[float]` — lift coefficients |
| `.cd` | `list[float]` — drag coefficients |
| `.cm` | `list[float]` — moment coefficients |
| `.ld` | `list[float]` — lift-to-drag ratios |
| `.converged` | `list[SolveResult]` — converged results only |
| `.results` | `list[SolveResult]` — all results |
| `.to_dict()` | Export as `dict` |
| `.to_dataframe()` | Export as `pandas.DataFrame` |
| `.plot(show=True, backend="plotly")` | 4-panel figure (plotly default, or `"matplotlib"`) |

### `RunDatabase`

| Method | Description |
| --- | --- |
| `.insert_run(...)` | Insert a solver run |
| `.lookup_cache(...)` | Cache lookup by (hash, alpha, Re, ...) |
| `.query_all_runs()` | All runs as `list[dict]` |
| `.query_runs(airfoil_name, limit, offset)` | Filtered query |
| `.row_count()` | Number of cached runs |
| `.delete_all_runs()` | Clear all runs |
| `.save_airfoil(name, coords_json)` | Save a named airfoil |
| `.list_airfoils()` | List saved airfoils |
| `.export_bytes()` | Export SQLite as bytes |
| `.import_bytes(data)` | Import SQLite from bytes |

## Configuration

| Environment variable | Default | Description |
| --- | --- | --- |
| `FLEXFOIL_DATA_DIR` | `~/.flexfoil` | Directory for `runs.db` |

## How it works

```
┌──────────────────────────────────┐     ┌──────────────────────────┐
│  Python                          │     │  Browser                 │
│                                  │     │                          │
│  import flexfoil                 │     │  localhost:8420           │
│  foil.solve() ──► Rust solver    │     │  React SPA               │
│        │         (PyO3)          │     │  WASM solver ──┐         │
│        ▼                         │     │                │         │
│  ~/.flexfoil/runs.db ◄──────────────── REST API ◄──────┘         │
│  (SQLite, WAL mode)              │     │                          │
│                                  │     │  SSE live updates ◄──┐  │
│  flexfoil.serve() ───► Starlette ──────►                     │  │
│                        uvicorn   │     │  (runs appear as     │  │
│                                  │     │   Python solves them) │  │
│  foil.polar() ──► insert rows ───────────────────────────────┘  │
└──────────────────────────────────┘     └──────────────────────────┘
```

The Rust solver called from Python is the exact same code as the WASM solver
in the browser — both are compiled from the `rustfoil-xfoil` crate. Results
are byte-identical.

## Examples

See the [`examples/`](examples/) directory:

| Script | What it does |
| --- | --- |
| `01_quickstart.py` | Single-point solve |
| `02_polar_sweep.py` | Full polar with table + plot |
| `03_compare_airfoils.py` | Overlay multiple NACA foils |
| `04_reynolds_sweep.py` | Re effect on drag polar |
| `05_dat_file.py` | Load from `.dat` file |
| `06_pandas_export.py` | Export to CSV via pandas |
| `07_inviscid_vs_viscous.py` | Inviscid vs viscous comparison |
| `08_batch_matrix.py` | Batch sweep: airfoils x Re x alpha |
| `09_custom_coordinates.py` | Build airfoil from x, y arrays |

## Supported platforms

| Platform | Status |
| --- | --- |
| macOS Apple Silicon (arm64) | Pre-built wheel |
| macOS Intel (x86_64) | Pre-built wheel |
| Linux x86_64 | Pre-built wheel |
| Windows x86_64 | Pre-built wheel |
| Linux aarch64 (ARM) | Build from source (requires Rust toolchain) |

## Development

```bash
cd packages/flexfoil-python

# Create a venv
python -m venv .venv
source .venv/bin/activate  # or .venv\Scripts\activate on Windows

# Install Rust toolchain and maturin
curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
pip install maturin

# Install dependencies
pip install starlette 'uvicorn[standard]' matplotlib pandas pytest

# Build and install in development mode
maturin develop

# Run tests (28 tests)
pytest tests/ -v

# Bundle the web UI (optional, for `flexfoil serve`)
./build_frontend.sh
```
