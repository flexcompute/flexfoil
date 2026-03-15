#!/usr/bin/env python3
"""Interactive boundary-layer explorer: hover on the airfoil to see Cp, Cf, δ*, θ, H.

Generates a self-contained Plotly HTML that renders:
  1. Airfoil shape with δ* displacement-thickness envelope (RustFoil + XFOIL)
  2. Cp distribution (inverted y-axis)
  3. Cf distribution (skin friction)
  4. H  distribution (shape factor)

Every surface point carries a rich hover tooltip with all BL quantities.

Usage:
    python scripts/visualize_bl_interactive.py                          # NACA 0012 at α=5° Re=3M
    python scripts/visualize_bl_interactive.py --naca 2412 --alpha 8    # NACA 2412 at α=8°
    python scripts/visualize_bl_interactive.py --dat my_foil.dat -a 3   # custom airfoil
    python scripts/visualize_bl_interactive.py --alphas -2,0,2,4,6,8    # multi-alpha sweep
"""

from __future__ import annotations

import argparse
import json
import math
import os
import re as regex
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

SCRIPT_DIR = Path(__file__).resolve().parent
PROJECT_DIR = SCRIPT_DIR.parent
XFOIL_PATH = PROJECT_DIR / "Xfoil" / "bin" / "xfoil"
OUTPUT_DIR = PROJECT_DIR / "comparison_results" / "bl_explorer"


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SurfaceBL:
    """Boundary-layer quantities along one surface (upper or lower)."""
    x: list[float] = field(default_factory=list)     # x/c coordinate
    y: list[float] = field(default_factory=list)     # y/c coordinate
    s: list[float] = field(default_factory=list)     # arc length
    ue: list[float] = field(default_factory=list)    # Ue/Vinf
    dstar: list[float] = field(default_factory=list) # displacement thickness
    theta: list[float] = field(default_factory=list) # momentum thickness
    cf: list[float] = field(default_factory=list)    # skin friction
    h: list[float] = field(default_factory=list)     # shape factor H = δ*/θ

    @property
    def cp(self) -> list[float]:
        """Pressure coefficient from edge velocity (Bernoulli, incompressible)."""
        return [1.0 - ue_val ** 2 for ue_val in self.ue]


@dataclass
class CaseData:
    """All data for one operating-point comparison."""
    foil_name: str
    alpha: float
    re: float
    # XFOIL results
    xfoil_cl: Optional[float] = None
    xfoil_cd: Optional[float] = None
    xfoil_upper: Optional[SurfaceBL] = None
    xfoil_lower: Optional[SurfaceBL] = None
    # RustFoil results
    rustfoil_cl: Optional[float] = None
    rustfoil_cd: Optional[float] = None
    rustfoil_upper: Optional[SurfaceBL] = None
    rustfoil_lower: Optional[SurfaceBL] = None
    # Airfoil geometry (full panel loop from .dat file)
    geom_x: Optional[list[float]] = None
    geom_y: Optional[list[float]] = None


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------

def load_airfoil_dat(path: Path) -> tuple[list[float], list[float]]:
    """Load Selig-format airfoil coordinates, return (x_list, y_list)."""
    xs, ys = [], []
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) >= 2:
                try:
                    xs.append(float(parts[0]))
                    ys.append(float(parts[1]))
                except ValueError:
                    continue
    return xs, ys


def split_upper_lower(xs: list[float], ys: list[float]):
    """Split full-loop airfoil into upper and lower surfaces.

    Assumes Selig ordering: TE -> upper LE -> lower TE.
    The leading-edge point is the one with minimum x.
    Returns (upper_x, upper_y, lower_x, lower_y) both ordered LE -> TE.
    """
    le_idx = int(np.argmin(xs))
    # Upper: indices 0..le_idx  (TE -> LE), reverse to get LE -> TE
    ux = list(reversed(xs[: le_idx + 1]))
    uy = list(reversed(ys[: le_idx + 1]))
    # Lower: indices le_idx..end (LE -> TE)
    lx = list(xs[le_idx:])
    ly = list(ys[le_idx:])
    return ux, uy, lx, ly


def cumulative_arc_length(xs: list[float], ys: list[float]) -> list[float]:
    s = [0.0]
    for i in range(1, len(xs)):
        ds = math.hypot(xs[i] - xs[i - 1], ys[i] - ys[i - 1])
        s.append(s[-1] + ds)
    return s


# ---------------------------------------------------------------------------
# XFOIL runner
# ---------------------------------------------------------------------------

def generate_naca_dat(naca: str) -> Path:
    """Use XFOIL to write a .dat file for a NACA airfoil, return the path."""
    dat_path = Path(tempfile.gettempdir()) / f"naca{naca}.dat"
    script = f"NACA {naca}\nPSAV {dat_path}\n\nQUIT\n"
    subprocess.run(
        [str(XFOIL_PATH)], input=script,
        capture_output=True, text=True, timeout=10,
    )
    if not dat_path.exists():
        raise RuntimeError(f"XFOIL failed to generate {dat_path}")
    return dat_path


def run_xfoil_bl(dat_path: Path, alpha: float, re: float, naca: Optional[str] = None) -> dict:
    """Run XFOIL viscous + DUMP to get BL data. Returns dict with cl, cd, bl."""
    dump_file = f"/tmp/_xfoil_dump_{os.getpid()}.dat"
    if os.path.exists(dump_file):
        os.remove(dump_file)

    load_cmd = f"NACA {naca}" if naca else f"LOAD {dat_path}"
    script = (
        f"{load_cmd}\nPANE\nOPER\n"
        f"VISC {re:.0f}\nITER 200\nALFA {alpha}\n"
        f"DUMP {dump_file}\n\nQUIT\n"
    )
    result = subprocess.run(
        [str(XFOIL_PATH)], input=script,
        capture_output=True, text=True, timeout=30,
    )
    cl, cd = None, None
    for line in result.stdout.split("\n"):
        m = regex.search(r"CL\s*=\s*([-\d.]+)", line)
        if m:
            cl = float(m.group(1))
        m = regex.search(r"CD\s*=\s*([-\d.]+)", line)
        if m:
            cd = float(m.group(1))

    bl = _parse_xfoil_dump(dump_file) if os.path.exists(dump_file) else None
    if os.path.exists(dump_file):
        os.remove(dump_file)
    return {"cl": cl, "cd": cd, "bl": bl}


def _parse_xfoil_dump(filename: str) -> dict[str, SurfaceBL]:
    """Parse XFOIL DUMP file into upper/lower SurfaceBL objects.

    The dump is a continuous loop (TE -> upper -> LE -> lower -> TE) with
    columns: s, x, y, Ue/Vinf, Dstar, Theta, Cf, H [, ...extra].
    We split upper/lower at the minimum-x point (leading edge).
    """
    all_rows: list[tuple[float, ...]] = []
    with open(filename) as f:
        for line in f:
            stripped = line.strip()
            if not stripped or stripped.startswith("#"):
                continue
            # Skip any text-header lines
            if "Ue" in stripped and "Dstar" in stripped:
                continue
            if "Side" in stripped or "UPPER" in stripped.upper() or "LOWER" in stripped.upper():
                continue
            parts = stripped.split()
            if len(parts) >= 8:
                try:
                    vals = tuple(float(p) for p in parts[:8])
                    all_rows.append(vals)
                except ValueError:
                    continue

    if not all_rows:
        return {"upper": SurfaceBL(), "lower": SurfaceBL()}

    # Find LE (minimum x) to split upper / lower
    xs = [r[1] for r in all_rows]
    le_idx = int(np.argmin(xs))

    def rows_to_surface(rows: list[tuple[float, ...]]) -> SurfaceBL:
        bl = SurfaceBL()
        for r in rows:
            bl.s.append(r[0])
            bl.x.append(r[1])
            bl.y.append(r[2])
            bl.ue.append(r[3])
            bl.dstar.append(r[4])
            bl.theta.append(r[5])
            bl.cf.append(r[6])
            bl.h.append(r[7])
        return bl

    # Upper surface: rows 0..le_idx (TE -> LE), reverse to get LE -> TE
    upper = rows_to_surface(list(reversed(all_rows[: le_idx + 1])))
    # Lower surface: rows le_idx..end (LE -> TE), strip wake (x > 1.0005)
    lower_rows = [r for r in all_rows[le_idx:] if r[1] <= 1.0005]
    lower = rows_to_surface(lower_rows)

    return {"upper": upper, "lower": lower}


# ---------------------------------------------------------------------------
# RustFoil runner
# ---------------------------------------------------------------------------

def run_rustfoil_bl(dat_path: Path, alpha: float, re: float) -> Optional[dict]:
    """Run RustFoil viscous with --dump-bl. Returns parsed JSON or None."""
    result = subprocess.run(
        [
            "cargo", "run", "--release", "-q", "-p", "rustfoil-cli", "--",
            "viscous", str(dat_path),
            f"--alpha={alpha}", f"--re={re:.0f}",
            "--format", "json", "--dump-bl",
        ],
        capture_output=True, text=True, timeout=120,
        cwd=str(PROJECT_DIR),
    )
    if result.returncode != 0:
        print(f"  RustFoil failed (exit {result.returncode})", file=sys.stderr)
        if result.stderr:
            for line in result.stderr.strip().split("\n")[-3:]:
                print(f"    {line}", file=sys.stderr)
        return None
    try:
        return json.loads(result.stdout)
    except json.JSONDecodeError:
        print("  RustFoil returned non-JSON output", file=sys.stderr)
        return None


def _rustfoil_bl_to_surface(
    bl_side: dict,
    geom_x: list[float],
    geom_y: list[float],
) -> SurfaceBL:
    """Convert RustFoil bl_state side dict to SurfaceBL with geometry coords.

    RustFoil bl_state gives arc-length x; we map to geometry (x, y) by
    computing arc lengths along the geometry surface and interpolating.
    """
    rf_arc = np.array(bl_side.get("x", []))
    rf_theta = np.array(bl_side.get("theta", []))
    rf_dstar = np.array(bl_side.get("delta_star", []))
    rf_ue = np.array(bl_side.get("ue", []))
    rf_hk = np.array(bl_side.get("hk", []))
    rf_cf = np.array(bl_side.get("cf", []))

    if len(rf_arc) == 0 or len(geom_x) < 2:
        return SurfaceBL()

    geom_arc = np.array(cumulative_arc_length(geom_x, geom_y))
    gx = np.array(geom_x)
    gy = np.array(geom_y)

    # Interpolate geometry coords at RustFoil arc-length positions
    interp_x = np.interp(rf_arc, geom_arc, gx)
    interp_y = np.interp(rf_arc, geom_arc, gy)

    # Compute H from δ*/θ (avoiding division by zero)
    h_vals = np.where(rf_theta > 1e-15, rf_dstar / rf_theta, 0.0)

    return SurfaceBL(
        x=interp_x.tolist(),
        y=interp_y.tolist(),
        s=rf_arc.tolist(),
        ue=rf_ue.tolist(),
        dstar=rf_dstar.tolist(),
        theta=rf_theta.tolist(),
        cf=rf_cf.tolist(),
        h=h_vals.tolist(),
    )


# ---------------------------------------------------------------------------
# Collect all data for one case
# ---------------------------------------------------------------------------

def collect_case(
    foil_name: str,
    dat_path: Path,
    alpha: float,
    re: float,
    naca: Optional[str] = None,
) -> CaseData:
    """Run both solvers and assemble a CaseData."""
    case = CaseData(foil_name=foil_name, alpha=alpha, re=re)

    # Geometry
    gx, gy = load_airfoil_dat(dat_path)
    case.geom_x, case.geom_y = gx, gy
    ux, uy, lx, ly = split_upper_lower(gx, gy)

    # --- XFOIL ---
    print(f"  XFOIL  α={alpha:+.1f}° ... ", end="", flush=True)
    xf = run_xfoil_bl(dat_path, alpha, re, naca=naca)
    case.xfoil_cl, case.xfoil_cd = xf["cl"], xf["cd"]
    if xf["bl"]:
        case.xfoil_upper = xf["bl"]["upper"]
        case.xfoil_lower = xf["bl"]["lower"]
    print(f"CL={case.xfoil_cl}, CD={case.xfoil_cd}")

    # --- RustFoil ---
    print(f"  Rust   α={alpha:+.1f}° ... ", end="", flush=True)
    rf = run_rustfoil_bl(dat_path, alpha, re)
    if rf:
        case.rustfoil_cl = rf.get("cl")
        case.rustfoil_cd = rf.get("cd")
        bl = rf.get("bl_state")
        if bl:
            case.rustfoil_upper = _rustfoil_bl_to_surface(bl.get("upper", {}), ux, uy)
            case.rustfoil_lower = _rustfoil_bl_to_surface(bl.get("lower", {}), lx, ly)
        print(f"CL={case.rustfoil_cl}, CD={case.rustfoil_cd}")
    else:
        print("FAILED")

    return case


# ---------------------------------------------------------------------------
# HTML generation
# ---------------------------------------------------------------------------

def _surface_to_json(surf: Optional[SurfaceBL]) -> Optional[dict]:
    if surf is None or not surf.x:
        return None
    return {
        "x": surf.x, "y": surf.y,
        "cp": surf.cp,
        "cf": surf.cf,
        "dstar": surf.dstar,
        "theta": surf.theta,
        "h": surf.h,
        "ue": surf.ue,
    }


def _case_to_json(case: CaseData) -> dict:
    return {
        "foil": case.foil_name,
        "alpha": case.alpha,
        "re": case.re,
        "xfoil_cl": case.xfoil_cl,
        "xfoil_cd": case.xfoil_cd,
        "rustfoil_cl": case.rustfoil_cl,
        "rustfoil_cd": case.rustfoil_cd,
        "geom_x": case.geom_x,
        "geom_y": case.geom_y,
        "xfoil_upper": _surface_to_json(case.xfoil_upper),
        "xfoil_lower": _surface_to_json(case.xfoil_lower),
        "rustfoil_upper": _surface_to_json(case.rustfoil_upper),
        "rustfoil_lower": _surface_to_json(case.rustfoil_lower),
    }


def write_html(cases: list[CaseData], output_path: Path) -> None:
    """Write self-contained interactive HTML with Plotly."""
    data_json = json.dumps(
        [_case_to_json(c) for c in cases],
        separators=(",", ":"),
    )

    html = _HTML_TEMPLATE.replace("__DATA_JSON__", data_json)
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        f.write(html)
    print(f"\nWrote {output_path}")


# ---------------------------------------------------------------------------
# HTML template (Plotly CDN, dark theme matching full_sweep_plot.html)
# ---------------------------------------------------------------------------

_HTML_TEMPLATE = r"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Boundary-Layer Explorer</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root {
      --bg:#08101f; --panel:#101a36; --panel2:#18254a;
      --text:#edf2ff; --muted:#9bacdf;
      --grid:rgba(155,172,223,.16);
      --xfoil:#7dd3fc; --xfoil2:#38bdf8;
      --rust:#f472b6; --rust2:#ec4899;
      --upper-opacity:1; --lower-opacity:.55;
    }
    *{box-sizing:border-box}
    body{
      margin:0;
      font-family:ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif;
      color:var(--text);
      background:radial-gradient(circle at top left,rgba(125,211,252,.12),transparent 28%),
                 radial-gradient(circle at top right,rgba(244,114,182,.08),transparent 24%),
                 linear-gradient(180deg,#09101e,#08101f);
    }
    .wrap{max-width:1500px;margin:0 auto;padding:28px 24px 48px}
    h1{margin:0 0 4px;font-size:28px;font-weight:760;letter-spacing:-.03em}
    .sub{margin:0 0 20px;color:var(--muted);font-size:14px;line-height:1.5;max-width:1000px}
    /* Controls row */
    .controls{display:flex;gap:14px;flex-wrap:wrap;align-items:end;margin-bottom:18px}
    .field{background:linear-gradient(180deg,rgba(24,37,74,.92),rgba(16,26,54,.94));
           border:1px solid rgba(155,172,223,.14);border-radius:14px;padding:12px 14px;
           box-shadow:0 8px 28px rgba(0,0,0,.18)}
    label{display:block;font-size:11px;font-weight:700;letter-spacing:.08em;
          text-transform:uppercase;color:var(--muted);margin-bottom:8px}
    select{width:100%;min-width:260px;padding:10px 12px;border-radius:10px;
           border:1px solid rgba(155,172,223,.18);background:rgba(10,15,30,.8);
           color:var(--text);font-size:14px}
    /* Summary cards */
    .summary{display:grid;grid-template-columns:repeat(auto-fit,minmax(170px,1fr));gap:14px;margin-bottom:18px}
    .card{background:linear-gradient(180deg,rgba(24,37,74,.92),rgba(16,26,54,.94));
          border:1px solid rgba(155,172,223,.14);border-radius:14px;padding:14px 16px;
          box-shadow:0 8px 28px rgba(0,0,0,.18)}
    .eyebrow{font-size:11px;letter-spacing:.08em;text-transform:uppercase;color:var(--muted);margin-bottom:8px}
    .val{font-size:22px;font-weight:760;letter-spacing:-.03em}
    .detail{color:var(--muted);font-size:12px;margin-top:4px}
    /* Plots */
    .plot-row{display:grid;gap:14px;margin-bottom:14px}
    .plot-row.top{grid-template-columns:1fr}
    .plot-row.bottom{grid-template-columns:1fr 1fr 1fr}
    .plot-card{background:linear-gradient(180deg,rgba(24,37,74,.92),rgba(16,26,54,.94));
               border:1px solid rgba(155,172,223,.14);border-radius:14px;padding:10px 10px 4px;
               box-shadow:0 8px 28px rgba(0,0,0,.18)}
    .plot-main{height:420px}
    .plot-sm{height:300px}
    /* Legend strip */
    .legend-strip{display:flex;gap:24px;justify-content:center;margin:-4px 0 10px;font-size:13px;color:var(--muted)}
    .legend-strip span{display:flex;align-items:center;gap:6px}
    .legend-strip .swatch{width:18px;height:4px;border-radius:2px}
    @media(max-width:1100px){.plot-row.bottom{grid-template-columns:1fr}}
  </style>
</head>
<body>
<div class="wrap">
  <h1>Boundary-Layer Explorer</h1>
  <p class="sub">
    Hover on the airfoil surface to inspect local Cp, Cf, &delta;*, &theta;, and H.
    The shaded envelope shows displacement thickness (&delta;*) offset normal to the surface.
  </p>

  <div class="controls">
    <div class="field">
      <label for="caseSelect">Operating Point</label>
      <select id="caseSelect"></select>
    </div>
  </div>

  <div class="summary" id="cards">
    <div class="card"><div class="eyebrow">XFOIL C<sub>L</sub></div><div class="val" id="xCl">-</div></div>
    <div class="card"><div class="eyebrow">RustFoil C<sub>L</sub></div><div class="val" id="rCl">-</div></div>
    <div class="card"><div class="eyebrow">XFOIL C<sub>D</sub></div><div class="val" id="xCd">-</div></div>
    <div class="card"><div class="eyebrow">RustFoil C<sub>D</sub></div><div class="val" id="rCd">-</div></div>
    <div class="card"><div class="eyebrow">CL Error</div><div class="val" id="clErr">-</div><div class="detail" id="clErrDetail"></div></div>
  </div>

  <div class="legend-strip">
    <span><span class="swatch" style="background:var(--xfoil)"></span> XFOIL Upper</span>
    <span><span class="swatch" style="background:var(--xfoil);opacity:.55"></span> XFOIL Lower</span>
    <span><span class="swatch" style="background:var(--rust)"></span> RustFoil Upper</span>
    <span><span class="swatch" style="background:var(--rust);opacity:.55"></span> RustFoil Lower</span>
  </div>

  <div class="plot-row top">
    <div class="plot-card"><div id="airfoilPlot" class="plot-main"></div></div>
  </div>
  <div class="plot-row bottom">
    <div class="plot-card"><div id="cpPlot" class="plot-sm"></div></div>
    <div class="plot-card"><div id="cfPlot" class="plot-sm"></div></div>
    <div class="plot-card"><div id="hPlot"  class="plot-sm"></div></div>
  </div>
</div>

<script>
const DATA = __DATA_JSON__;

const plotCfg = {responsive:true, displaylogo:false,
  modeBarButtonsToRemove:["lasso2d","select2d","autoScale2d"]};

const baseLayout = {
  paper_bgcolor:"rgba(0,0,0,0)", plot_bgcolor:"rgba(0,0,0,0)",
  font:{color:"#edf2ff",size:12},
  margin:{l:52,r:16,t:36,b:42},
  xaxis:{gridcolor:"rgba(155,172,223,.16)",zerolinecolor:"rgba(155,172,223,.2)"},
  yaxis:{gridcolor:"rgba(155,172,223,.16)",zerolinecolor:"rgba(155,172,223,.2)"},
  legend:{orientation:"h",yanchor:"bottom",y:1.02,xanchor:"right",x:1,font:{size:11}},
  hovermode:"closest"
};

function fmtRe(re){ return "Re=" + (re>=1e6 ? (re/1e6).toFixed(0)+"M" : (re/1e3).toFixed(0)+"k"); }
function fmtCase(c){ return c.foil.toUpperCase() + "  α=" + c.alpha.toFixed(1) + "°  " + fmtRe(c.re); }

// Hover template shared by all BL surface traces
const hoverBL =
  "<b>%{meta}</b><br>" +
  "x/c = %{x:.4f}<br>" +
  "Cp = %{customdata[0]:.4f}<br>" +
  "Cf×10³ = %{customdata[1]:.3f}<br>" +
  "δ*×10³ = %{customdata[2]:.3f}<br>" +
  "θ×10³ = %{customdata[3]:.3f}<br>" +
  "H = %{customdata[4]:.2f}<br>" +
  "Ue/V∞ = %{customdata[5]:.4f}" +
  "<extra></extra>";

function customData(s){
  // pack [cp, cf*1e3, dstar*1e3, theta*1e3, H, ue]
  return s.x.map((_,i) => [
    s.cp[i], s.cf[i]*1e3, s.dstar[i]*1e3, s.theta[i]*1e3, s.h[i], s.ue[i]
  ]);
}

function envelopeTrace(surf, scale, color, name){
  // Compute δ* envelope offset normal to the surface.
  // Approximate local normal from finite differences.
  const n = surf.x.length;
  if(n < 3) return null;
  const ex = [], ey = [];
  for(let i = 0; i < n; i++){
    const ip = Math.min(i+1, n-1), im = Math.max(i-1, 0);
    const dx = surf.x[ip] - surf.x[im];
    const dy = surf.y[ip] - surf.y[im];
    const len = Math.sqrt(dx*dx + dy*dy) || 1e-12;
    // outward normal: rotate tangent 90° (upper: +, lower: −)
    const nx = -dy/len * scale;
    const ny =  dx/len * scale;
    ex.push(surf.x[i] + nx * surf.dstar[i]);
    ey.push(surf.y[i] + ny * surf.dstar[i]);
  }
  // closed fill: surface → envelope reversed → back to start
  const fx = surf.x.concat(ex.slice().reverse(), [surf.x[0]]);
  const fy = surf.y.concat(ey.slice().reverse(), [surf.y[0]]);
  return {
    x:fx, y:fy, mode:"lines", fill:"toself",
    fillcolor:color.replace(")",",0.12)").replace("rgb","rgba"),
    line:{color:color.replace(")",",0.25)").replace("rgb","rgba"), width:0},
    name:name, showlegend:false, hoverinfo:"skip"
  };
}

function renderCase(idx){
  const c = DATA[idx];

  // Cards
  document.getElementById("xCl").textContent = c.xfoil_cl != null ? c.xfoil_cl.toFixed(4) : "—";
  document.getElementById("rCl").textContent = c.rustfoil_cl != null ? c.rustfoil_cl.toFixed(4) : "—";
  document.getElementById("xCd").textContent = c.xfoil_cd != null ? c.xfoil_cd.toFixed(5) : "—";
  document.getElementById("rCd").textContent = c.rustfoil_cd != null ? c.rustfoil_cd.toFixed(5) : "—";
  if(c.xfoil_cl != null && c.rustfoil_cl != null && Math.abs(c.xfoil_cl) > 1e-9){
    const err = Math.abs((c.rustfoil_cl - c.xfoil_cl) / c.xfoil_cl * 100);
    document.getElementById("clErr").textContent = err.toFixed(2) + "%";
    document.getElementById("clErrDetail").textContent =
      "ΔCL = " + (c.rustfoil_cl - c.xfoil_cl).toFixed(5);
  } else {
    document.getElementById("clErr").textContent = "—";
    document.getElementById("clErrDetail").textContent = "";
  }

  // --- Airfoil plot ---
  const airTraces = [];

  // Geometry outline (faint)
  if(c.geom_x){
    airTraces.push({
      x:c.geom_x, y:c.geom_y, mode:"lines",
      line:{color:"rgba(155,172,223,.25)",width:1},
      showlegend:false, hoverinfo:"skip"
    });
  }

  // δ* envelopes
  const envColors = {
    xfoil_upper:"rgb(125,211,252)", xfoil_lower:"rgb(125,211,252)",
    rustfoil_upper:"rgb(244,114,182)", rustfoil_lower:"rgb(244,114,182)"
  };
  for(const key of ["xfoil_upper","xfoil_lower","rustfoil_upper","rustfoil_lower"]){
    const surf = c[key];
    if(!surf) continue;
    const isUpper = key.endsWith("upper");
    const env = envelopeTrace(surf, isUpper ? 1 : -1, envColors[key],
      key.replace("_"," ").replace("xfoil","XFOIL").replace("rustfoil","RustFoil") + " δ*");
    if(env) airTraces.push(env);
  }

  // Interactive surface lines (carry hover)
  function addSurfTrace(surf, color, opacity, name){
    if(!surf || !surf.x.length) return;
    airTraces.push({
      x:surf.x, y:surf.y, mode:"lines+markers",
      marker:{size:4, color:color, opacity:opacity},
      line:{color:color, width:2.5, opacity:opacity},
      customdata:customData(surf), meta:name,
      hovertemplate:hoverBL,
      name:name, showlegend:false
    });
  }
  addSurfTrace(c.xfoil_upper,   "#7dd3fc", 1,    "XFOIL Upper");
  addSurfTrace(c.xfoil_lower,   "#7dd3fc", 0.55, "XFOIL Lower");
  addSurfTrace(c.rustfoil_upper, "#f472b6", 1,    "RustFoil Upper");
  addSurfTrace(c.rustfoil_lower, "#f472b6", 0.55, "RustFoil Lower");

  Plotly.newPlot("airfoilPlot", airTraces, {
    ...baseLayout,
    title:{text:"Airfoil + δ* Envelope", x:0.02, xanchor:"left"},
    xaxis:{...baseLayout.xaxis, title:"x/c", range:[-0.02, 1.05]},
    yaxis:{...baseLayout.yaxis, title:"y/c", scaleanchor:"x", scaleratio:1},
    margin:{l:52,r:16,t:36,b:42}
  }, plotCfg);

  // --- Cp plot ---
  const cpTraces = [];
  function addCpTrace(surf, color, opacity, name){
    if(!surf || !surf.x.length) return;
    cpTraces.push({
      x:surf.x, y:surf.cp, mode:"lines+markers",
      marker:{size:3, color:color, opacity:opacity},
      line:{color:color, width:2, opacity:opacity},
      customdata:customData(surf), meta:name,
      hovertemplate:hoverBL, name:name, showlegend:false
    });
  }
  addCpTrace(c.xfoil_upper,   "#7dd3fc", 1,    "XFOIL Upper");
  addCpTrace(c.xfoil_lower,   "#7dd3fc", 0.55, "XFOIL Lower");
  addCpTrace(c.rustfoil_upper, "#f472b6", 1,    "RustFoil Upper");
  addCpTrace(c.rustfoil_lower, "#f472b6", 0.55, "RustFoil Lower");

  Plotly.newPlot("cpPlot", cpTraces, {
    ...baseLayout,
    title:{text:"Cp Distribution", x:0.02, xanchor:"left"},
    xaxis:{...baseLayout.xaxis, title:"x/c", range:[-0.02,1.05]},
    yaxis:{...baseLayout.yaxis, title:"Cp", autorange:"reversed"},
    margin:{l:52,r:12,t:36,b:42}
  }, plotCfg);

  // --- Cf plot ---
  const cfTraces = [];
  function addCfTrace(surf, color, opacity, name){
    if(!surf || !surf.x.length) return;
    cfTraces.push({
      x:surf.x, y:surf.cf.map(v=>v*1e3), mode:"lines+markers",
      marker:{size:3, color:color, opacity:opacity},
      line:{color:color, width:2, opacity:opacity},
      customdata:customData(surf), meta:name,
      hovertemplate:hoverBL, name:name, showlegend:false
    });
  }
  addCfTrace(c.xfoil_upper,   "#7dd3fc", 1,    "XFOIL Upper");
  addCfTrace(c.xfoil_lower,   "#7dd3fc", 0.55, "XFOIL Lower");
  addCfTrace(c.rustfoil_upper, "#f472b6", 1,    "RustFoil Upper");
  addCfTrace(c.rustfoil_lower, "#f472b6", 0.55, "RustFoil Lower");

  // Cf = 0 separation line
  cfTraces.push({
    x:[0,1], y:[0,0], mode:"lines",
    line:{color:"rgba(251,113,133,.5)", width:1.5, dash:"dot"},
    showlegend:false, hoverinfo:"skip"
  });

  Plotly.newPlot("cfPlot", cfTraces, {
    ...baseLayout,
    title:{text:"Cf × 10³ (Cf < 0 → separated)", x:0.02, xanchor:"left"},
    xaxis:{...baseLayout.xaxis, title:"x/c", range:[-0.02,1.05]},
    yaxis:{...baseLayout.yaxis, title:"Cf × 10³", range:[-5,25]},
    margin:{l:52,r:12,t:36,b:42}
  }, plotCfg);

  // --- H plot ---
  const hTraces = [];
  function addHTrace(surf, color, opacity, name){
    if(!surf || !surf.x.length) return;
    hTraces.push({
      x:surf.x, y:surf.h, mode:"lines+markers",
      marker:{size:3, color:color, opacity:opacity},
      line:{color:color, width:2, opacity:opacity},
      customdata:customData(surf), meta:name,
      hovertemplate:hoverBL, name:name, showlegend:false
    });
  }
  addHTrace(c.xfoil_upper,   "#7dd3fc", 1,    "XFOIL Upper");
  addHTrace(c.xfoil_lower,   "#7dd3fc", 0.55, "XFOIL Lower");
  addHTrace(c.rustfoil_upper, "#f472b6", 1,    "RustFoil Upper");
  addHTrace(c.rustfoil_lower, "#f472b6", 0.55, "RustFoil Lower");

  // H = 4 separation threshold
  hTraces.push({
    x:[0,1], y:[4,4], mode:"lines",
    line:{color:"rgba(251,113,133,.5)", width:1.5, dash:"dot"},
    showlegend:false, hoverinfo:"skip"
  });

  Plotly.newPlot("hPlot", hTraces, {
    ...baseLayout,
    title:{text:"Shape Factor H (H > 4 → separation)", x:0.02, xanchor:"left"},
    xaxis:{...baseLayout.xaxis, title:"x/c", range:[-0.02,1.05]},
    yaxis:{...baseLayout.yaxis, title:"H", range:[0, 10]},
    margin:{l:52,r:12,t:36,b:42}
  }, plotCfg);
}

// Populate selector
const sel = document.getElementById("caseSelect");
DATA.forEach((c, i) => {
  const opt = document.createElement("option");
  opt.value = i; opt.textContent = fmtCase(c);
  sel.appendChild(opt);
});
sel.addEventListener("change", e => renderCase(+e.target.value));
if(DATA.length) renderCase(0);
</script>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Interactive BL explorer: hover on the airfoil to see Cp, Cf, δ*, θ, H",
    )
    grp = parser.add_mutually_exclusive_group()
    grp.add_argument("--naca", default="0012", help="NACA 4-digit code (default: 0012)")
    grp.add_argument("--dat", type=Path, help="Path to an airfoil .dat file")
    parser.add_argument("-a", "--alpha", type=float, default=None,
                        help="Single angle of attack (default: 5)")
    parser.add_argument("--alphas", type=str, default=None,
                        help="Comma-separated list of alphas (e.g. -2,0,2,4,6,8)")
    parser.add_argument("--re", type=float, default=3e6, help="Reynolds number (default: 3e6)")
    parser.add_argument("-o", "--output", type=Path, default=None,
                        help="Output HTML path (default: comparison_results/bl_explorer/bl_explorer.html)")
    args = parser.parse_args()

    # Resolve alpha list
    if args.alphas:
        alphas = [float(a.strip()) for a in args.alphas.split(",")]
    elif args.alpha is not None:
        alphas = [args.alpha]
    else:
        alphas = [5.0]

    # Resolve airfoil
    if args.dat:
        dat_path = args.dat.resolve()
        foil_name = dat_path.stem
    else:
        foil_name = f"NACA {args.naca}"
        print(f"Generating {foil_name} geometry via XFOIL ...")
        dat_path = generate_naca_dat(args.naca)

    output_path = args.output or (OUTPUT_DIR / "bl_explorer.html")

    print(f"\n{'='*60}")
    print(f"  BL Explorer: {foil_name}")
    print(f"  Re = {args.re:.2e}    α = {', '.join(f'{a:+.1f}' for a in alphas)}°")
    print(f"{'='*60}\n")

    naca_code = args.naca if not args.dat else None
    cases: list[CaseData] = []
    for alpha in alphas:
        case = collect_case(foil_name, dat_path, alpha, args.re, naca=naca_code)
        cases.append(case)

    write_html(cases, output_path)
    print(f"\nOpen in your browser:\n  file://{output_path.resolve()}")


if __name__ == "__main__":
    main()
