#!/usr/bin/env python3
"""
Rebuild and compare the full RustFoil vs XFOIL viscous sweep matrix.

This is the standard end-to-end script intended to be run after a fix:

    python scripts/run_full_sweep_comparison.py

It will:
1. Regenerate XFOIL traces for the standard 5-case matrix
2. Regenerate RustFoil traces for the same matrix
3. Build a machine-readable comparison summary
4. Build a self-contained Plotly HTML dashboard

By default, outputs go to timestamped directories so each run is preserved.
Use --tag latest (or another stable tag) if you want predictable paths.

IMPORTANT:
This script compares against the faithful rustfoil-xfoil side path only.
It is intentionally not able to invoke the legacy rustfoil viscous solver.
For large matrices, compact summary-only traces are used by default.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from statistics import mean
from typing import Any

from sweep_matrix import case_tuples


@dataclass
class SweepPaths:
    xfoil_dir: Path
    rustfoil_dir: Path
    comparison_dir: Path
    summary_json: Path
    plot_html: Path
    markdown_report: Path


def load_json(path: Path) -> dict[str, Any]:
    with open(path) as f:
        return json.load(f)


def run_command(command: list[str], cwd: Path) -> None:
    print()
    print("$", " ".join(command))
    result = subprocess.run(command, cwd=cwd)
    if result.returncode != 0:
        raise SystemExit(result.returncode)


def run_commands_parallel(commands: list[tuple[str, list[str]]], cwd: Path) -> None:
    procs: list[tuple[str, subprocess.Popen[str]]] = []
    try:
        for label, command in commands:
            print()
            print(f"$ ({label})", " ".join(command))
            proc = subprocess.Popen(command, cwd=cwd, text=True)
            procs.append((label, proc))

        while procs:
            still_running: list[tuple[str, subprocess.Popen[str]]] = []
            for label, proc in procs:
                code = proc.poll()
                if code is None:
                    still_running.append((label, proc))
                    continue
                if code != 0:
                    for other_label, other_proc in still_running:
                        if other_proc.poll() is None:
                            other_proc.terminate()
                    raise SystemExit(f"{label} failed with exit code {code}")
            procs = still_running
            if procs:
                time.sleep(0.5)
    finally:
        for _, proc in procs:
            if proc.poll() is None:
                proc.terminate()


def extract_trace_metadata(path: Path) -> dict[str, Any]:
    data = load_json(path)
    meta = data.get("metadata", {})
    cl = meta.get("cl")
    cd = meta.get("cd")
    converged = meta.get("converged")
    solver_path = meta.get("solver_path")
    if cl is not None and cd is not None:
        return {
            "cl": cl,
            "cd": cd,
            "converged": True if converged is None else converged,
            "solver_path": solver_path,
        }

    for event in reversed(data.get("events", [])):
        subroutine = event.get("subroutine")
        if subroutine in {"VISCOUS_FINAL", "VISCAL_RESULT"}:
            return {
                "cl": event.get("CL", event.get("cl")),
                "cd": event.get("CD", event.get("cd")),
                "converged": event.get("converged", True),
                "solver_path": solver_path,
            }

    return {"cl": None, "cd": None, "converged": False, "solver_path": solver_path}


def build_case_points(
    xfoil_dir: Path,
    rustfoil_dir: Path,
    foil: str,
    re_dir: str,
    expected_rust_solver_path: str,
) -> list[dict[str, Any]]:
    xfoil_summary = load_json(xfoil_dir / foil / re_dir / "summary.json")
    rustfoil_summary = load_json(rustfoil_dir / foil / re_dir / "summary.json")

    xfoil_by_alpha = {round(row["alpha"], 6): row for row in xfoil_summary["results"]}
    points: list[dict[str, Any]] = []

    for rust_row in rustfoil_summary["results"]:
        alpha = round(rust_row["alpha"], 6)
        xfoil_row = xfoil_by_alpha.get(alpha)
        if xfoil_row is None:
            continue

        x_meta = extract_trace_metadata(Path(xfoil_row["output_file"]))
        r_meta = extract_trace_metadata(Path(rust_row["output_file"]))
        if None in (x_meta["cl"], x_meta["cd"], r_meta["cl"], r_meta["cd"]):
            continue
        if r_meta.get("solver_path") != expected_rust_solver_path:
            raise RuntimeError(
                f"Unexpected Rust trace solver path for {foil} {re_dir} alpha={alpha:+.1f}: "
                f"{r_meta.get('solver_path')!r} != {expected_rust_solver_path!r}"
            )

        xfoil_cl = x_meta["cl"]
        rustfoil_cl = r_meta["cl"]
        xfoil_cd = x_meta["cd"]
        rustfoil_cd = r_meta["cd"]
        if not all(
            math.isfinite(value)
            for value in (xfoil_cl, rustfoil_cl, xfoil_cd, rustfoil_cd)
        ):
            continue
        cl_error_pct = (
            abs((rustfoil_cl - xfoil_cl) / xfoil_cl) * 100.0
            if abs(xfoil_cl) > 1.0e-9
            else abs(rustfoil_cl - xfoil_cl) * 100.0
        )
        cd_ratio = rustfoil_cd / xfoil_cd if xfoil_cd else None

        points.append(
            {
                "alpha": alpha,
                "xfoil_cl": xfoil_cl,
                "rustfoil_cl": rustfoil_cl,
                "xfoil_cd": xfoil_cd,
                "rustfoil_cd": rustfoil_cd,
                "cl_error_pct": cl_error_pct,
                "cd_ratio": cd_ratio,
                "xfoil_converged": bool(x_meta["converged"]),
                "rustfoil_converged": bool(r_meta["converged"]),
            }
        )

    return points


def summarize_case(case: dict[str, Any]) -> dict[str, Any]:
    if case.get("error"):
        return {
            "foil": case["foil"],
            "re": case["re"],
            "error": case["error"],
            "point_count": 0,
            "rustfoil_converged_points": 0,
        }

    points = case["points"]
    if not points:
        return {
            "foil": case["foil"],
            "re": case["re"],
            "error": "No comparable points",
            "point_count": 0,
            "rustfoil_converged_points": 0,
        }
    cl_errors = [point["cl_error_pct"] for point in points]
    cd_ratios = [point["cd_ratio"] for point in points if point["cd_ratio"] is not None]
    worst_cl = max(points, key=lambda point: point["cl_error_pct"])
    worst_cd = max(points, key=lambda point: point["cd_ratio"])
    best_cd = min(points, key=lambda point: point["cd_ratio"])

    return {
        "foil": case["foil"],
        "re": case["re"],
        "point_count": len(points),
        "rustfoil_converged_points": sum(1 for point in points if point["rustfoil_converged"]),
        "avg_cl_error_pct": mean(cl_errors),
        "max_cl_error_pct": max(cl_errors),
        "avg_cd_ratio": mean(cd_ratios),
        "min_cd_ratio": min(cd_ratios),
        "max_cd_ratio": max(cd_ratios),
        "worst_cl_point": worst_cl,
        "worst_cd_point": worst_cd,
        "best_cd_point": best_cd,
    }


def build_summary(cases: list[dict[str, Any]]) -> dict[str, Any]:
    case_summaries = [summarize_case(case) for case in cases]
    all_points = [point for case in cases for point in case.get("points", [])]
    valid_cd = [point["cd_ratio"] for point in all_points if point["cd_ratio"] is not None]
    return {
        "generated_at": datetime.now().isoformat(),
        "cases": cases,
        "case_summaries": case_summaries,
        "overall": {
            "point_count": len(all_points),
            "rustfoil_converged_points": sum(1 for point in all_points if point["rustfoil_converged"]),
            "avg_cl_error_pct": mean(point["cl_error_pct"] for point in all_points) if all_points else None,
            "max_cl_error_pct": max(point["cl_error_pct"] for point in all_points) if all_points else None,
            "avg_cd_ratio": mean(valid_cd) if valid_cd else None,
            "min_cd_ratio": min(valid_cd) if valid_cd else None,
            "max_cd_ratio": max(valid_cd) if valid_cd else None,
        },
    }


def write_markdown_report(summary: dict[str, Any], output_path: Path) -> None:
    overall = summary["overall"]
    lines = [
        "# Full Sweep Comparison",
        "",
        f"- Generated: `{summary['generated_at']}`",
        f"- Rust solver path: `{summary['rust_solver_path']}`",
        f"- Total points: `{overall['point_count']}`",
        f"- RustFoil converged points: `{overall['rustfoil_converged_points']}`",
        f"- Avg CL error: `{overall['avg_cl_error_pct']:.2f}%`",
        f"- Max CL error: `{overall['max_cl_error_pct']:.2f}%`",
        f"- Avg CD ratio: `{overall['avg_cd_ratio']:.2f}x`",
        f"- CD ratio range: `{overall['min_cd_ratio']:.2f}x .. {overall['max_cd_ratio']:.2f}x`",
        "",
        "## Cases",
        "",
        "| Foil | Re | Points | RF Converged | Avg CL Err | Max CL Err | Avg CD Ratio | CD Range |",
        "|---|---:|---:|---:|---:|---:|---:|---|",
    ]

    for case in summary["case_summaries"]:
        if case.get("error"):
            lines.append(
                f"| `{case['foil']}` | `{case['re']}` | 0 | 0 | - | - | - | `{case['error']}` |"
            )
        else:
            lines.append(
                f"| `{case['foil']}` | `{case['re']}` | {case['point_count']} | "
                f"{case['rustfoil_converged_points']} | {case['avg_cl_error_pct']:.2f}% | "
                f"{case['max_cl_error_pct']:.2f}% | {case['avg_cd_ratio']:.2f}x | "
                f"{case['min_cd_ratio']:.2f}x .. {case['max_cd_ratio']:.2f}x |"
            )

    output_path.write_text("\n".join(lines) + "\n")


def write_plot_html(cases: list[dict[str, Any]], output_path: Path, rust_solver_path: str) -> None:
    plot_cases = [case for case in cases if case.get("points")]
    html = f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>RustFoil vs XFOIL Full Sweep Comparison</title>
  <script src="https://cdn.plot.ly/plotly-2.35.2.min.js"></script>
  <style>
    :root {{
      --bg:#08101f; --panel:#101a36; --panel2:#18254a; --text:#edf2ff; --muted:#9bacdf; --grid:rgba(155,172,223,.16);
    }}
    * {{ box-sizing:border-box; }}
    body {{
      margin:0; font-family:ui-sans-serif,system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif; color:var(--text);
      background:radial-gradient(circle at top left, rgba(125,211,252,.12), transparent 28%),
                 radial-gradient(circle at top right, rgba(244,114,182,.08), transparent 24%),
                 linear-gradient(180deg,#09101e 0%,#08101f 100%);
    }}
    .wrap {{ max-width:1500px; margin:0 auto; padding:28px 24px 48px; }}
    h1 {{ margin:0 0 8px; font-size:30px; font-weight:760; letter-spacing:-.03em; }}
    .sub {{ margin:0 0 24px; color:var(--muted); font-size:15px; line-height:1.5; max-width:1000px; }}
    .controls,.summary,.plot-grid {{ display:grid; gap:16px; }}
    .controls {{ grid-template-columns:minmax(280px,420px) 1fr; align-items:end; margin-bottom:18px; }}
    .field,.card,.plot-card {{ background:linear-gradient(180deg, rgba(24,37,74,.92), rgba(16,26,54,.94)); border:1px solid rgba(155,172,223,.14); border-radius:16px; box-shadow:0 12px 40px rgba(0,0,0,.22); }}
    .field {{ padding:14px 16px; }}
    .note {{ padding:14px 16px; color:var(--muted); font-size:14px; line-height:1.5; }}
    label {{ display:block; font-size:12px; font-weight:700; letter-spacing:.08em; text-transform:uppercase; color:var(--muted); margin-bottom:10px; }}
    select {{ width:100%; padding:12px 14px; border-radius:12px; border:1px solid rgba(155,172,223,.18); background:rgba(10,15,30,.8); color:var(--text); font-size:15px; }}
    .summary {{ grid-template-columns:repeat(5,minmax(0,1fr)); margin-bottom:18px; }}
    .card {{ padding:16px 18px; min-height:108px; }}
    .eyebrow {{ font-size:12px; letter-spacing:.08em; text-transform:uppercase; color:var(--muted); margin-bottom:10px; }}
    .value {{ font-size:26px; font-weight:760; letter-spacing:-.03em; margin-bottom:6px; }}
    .detail {{ color:var(--muted); font-size:13px; line-height:1.45; }}
    .plot-grid {{ grid-template-columns:1fr 1fr; }}
    .plot-card {{ padding:12px 12px 6px; }}
    .plot {{ height:360px; }}
    .dot-row {{ display:none; justify-content:center; gap:8px; margin:-8px 0 16px; }}
    .dot {{ width:8px; height:8px; border-radius:50%; background:rgba(155,172,223,.25); transition:background .2s, transform .2s; }}
    .dot.active {{ background:#7dd3fc; transform:scale(1.3); }}
    .plot-tabs {{ display:none; gap:8px; margin-bottom:12px; overflow-x:auto; -webkit-overflow-scrolling:touch; }}
    .plot-tab {{ padding:10px 18px; border-radius:10px; font-size:13px; font-weight:600; letter-spacing:.02em; background:rgba(24,37,74,.6); border:1px solid rgba(155,172,223,.14); color:var(--muted); cursor:pointer; white-space:nowrap; transition:all .2s; }}
    .plot-tab.active {{ background:rgba(125,211,252,.15); border-color:rgba(125,211,252,.4); color:var(--text); }}
    @media (max-width:1100px) {{ .summary,.plot-grid,.controls {{ grid-template-columns:1fr; }} }}
    @media (max-width:600px) {{
      .wrap {{ padding:16px 12px 32px; }}
      h1 {{ font-size:22px; }}
      .sub {{ font-size:13px; margin-bottom:16px; }}
      select {{ min-height:44px; font-size:16px; }}
      .controls {{ grid-template-columns:1fr; gap:10px; margin-bottom:12px; }}
      .summary {{ display:flex; overflow-x:auto; scroll-snap-type:x mandatory; -webkit-overflow-scrolling:touch; gap:12px; padding-bottom:4px; margin-bottom:8px; scrollbar-width:none; }}
      .summary::-webkit-scrollbar {{ display:none; }}
      .card {{ min-width:78vw; flex-shrink:0; scroll-snap-align:start; min-height:auto; padding:14px 16px; }}
      .value {{ font-size:22px; }}
      .dot-row {{ display:flex; }}
      .plot-tabs {{ display:flex; }}
      .plot-grid {{ display:block; }}
      .plot-card {{ display:none; }}
      .plot-card.active {{ display:block; }}
      .plot {{ height:300px; }}
    }}
  </style>
</head>
<body>
  <div class="wrap">
    <h1>RustFoil vs XFOIL Full Sweep Comparison</h1>
    <p class="sub">
      Fresh rebuild from raw trace directories. Rust side path: <strong>{rust_solver_path}</strong>.
      Use this script output as the standard post-fix regression view.
      The lower panels make CL and drag regressions obvious even when raw CL/CD overlays look superficially close.
    </p>
    <div class="controls">
      <div class="field">
        <label for="caseSelect">Sweep Case</label>
        <select id="caseSelect"></select>
      </div>
      <div class="field note" id="caseNote"></div>
    </div>
    <div class="summary">
      <div class="card"><div class="eyebrow">Avg CL Error</div><div class="value" id="avgClErr">-</div><div class="detail" id="avgClErrDetail"></div></div>
      <div class="card"><div class="eyebrow">Worst CL Point</div><div class="value" id="worstCl">-</div><div class="detail" id="worstClDetail"></div></div>
      <div class="card"><div class="eyebrow">Avg CD Ratio</div><div class="value" id="avgCdRatio">-</div><div class="detail" id="avgCdRatioDetail"></div></div>
      <div class="card"><div class="eyebrow">CD Ratio Range</div><div class="value" id="cdRange">-</div><div class="detail" id="cdRangeDetail"></div></div>
      <div class="card"><div class="eyebrow">Convergence</div><div class="value" id="pointCount">-</div><div class="detail" id="pointCountDetail"></div></div>
    </div>
    <div class="dot-row">
      <div class="dot active"></div><div class="dot"></div><div class="dot"></div><div class="dot"></div><div class="dot"></div>
    </div>
    <nav class="plot-tabs">
      <div class="plot-tab active" data-target="clPlot">CL vs Alpha</div>
      <div class="plot-tab" data-target="cdPlot">CD vs Alpha</div>
      <div class="plot-tab" data-target="clErrPlot">CL Error %</div>
      <div class="plot-tab" data-target="cdRatioPlot">CD Ratio</div>
    </nav>
    <div class="plot-grid">
      <div class="plot-card active" data-plot="clPlot"><div id="clPlot" class="plot"></div></div>
      <div class="plot-card" data-plot="cdPlot"><div id="cdPlot" class="plot"></div></div>
      <div class="plot-card" data-plot="clErrPlot"><div id="clErrPlot" class="plot"></div></div>
      <div class="plot-card" data-plot="cdRatioPlot"><div id="cdRatioPlot" class="plot"></div></div>
    </div>
  </div>
  <script>
    const DATA = {json.dumps(plot_cases, separators=(",", ":"))};
    const caseSelect = document.getElementById("caseSelect");
    const caseNote = document.getElementById("caseNote");
    const isMobile = () => window.matchMedia("(max-width:600px)").matches;
    const plotLayoutBase = {{
      paper_bgcolor: "rgba(0,0,0,0)",
      plot_bgcolor: "rgba(0,0,0,0)",
      margin: isMobile() ? {{ l: 44, r: 12, t: 36, b: 40 }} : {{ l: 58, r: 20, t: 48, b: 48 }},
      font: {{ color: "#edf2ff", size: 13 }},
      xaxis: {{ gridcolor: "rgba(155,172,223,.16)", zerolinecolor: "rgba(155,172,223,.2)", title: "Alpha (deg)" }},
      yaxis: {{ gridcolor: "rgba(155,172,223,.16)", zerolinecolor: "rgba(155,172,223,.2)" }},
      legend: {{ orientation: "h", yanchor: "bottom", y: 1.02, xanchor: "right", x: 1 }}
    }};
    const plotConfig = {{ responsive: true, displaylogo: false, modeBarButtonsToRemove: ["lasso2d", "select2d", "autoScale2d"] }};
    function fmtRe(re) {{ return re.replace("re", "Re=").replace("e0", "e"); }}
    function fmtCaseLabel(item) {{ return `${{item.foil.toUpperCase()}} • ${{fmtRe(item.re)}}`; }}
    function avg(arr) {{ return arr.reduce((a, b) => a + b, 0) / arr.length; }}
    DATA.forEach((item, index) => {{
      const option = document.createElement("option");
      option.value = String(index);
      option.textContent = fmtCaseLabel(item);
      caseSelect.appendChild(option);
    }});
    function updateCards(item) {{
      const points = item.points;
      const clErrors = points.map(p => p.cl_error_pct);
      const cdRatios = points.map(p => p.cd_ratio);
      const worstCl = points.reduce((best, p) => p.cl_error_pct > best.cl_error_pct ? p : best, points[0]);
      const minCd = points.reduce((best, p) => p.cd_ratio < best.cd_ratio ? p : best, points[0]);
      const maxCd = points.reduce((best, p) => p.cd_ratio > best.cd_ratio ? p : best, points[0]);
      const rustConverged = points.filter(p => p.rustfoil_converged).length;
      document.getElementById("avgClErr").textContent = `${{avg(clErrors).toFixed(2)}}%`;
      document.getElementById("avgClErrDetail").textContent = `Min ${{Math.min(...clErrors).toFixed(2)}}%, max ${{Math.max(...clErrors).toFixed(2)}}%`;
      document.getElementById("worstCl").textContent = `${{worstCl.alpha.toFixed(0)}} deg`;
      document.getElementById("worstClDetail").textContent = `XFOIL CL ${{worstCl.xfoil_cl.toFixed(4)}} vs RustFoil ${{worstCl.rustfoil_cl.toFixed(4)}} (${{worstCl.cl_error_pct.toFixed(2)}}%)`;
      document.getElementById("avgCdRatio").textContent = `${{avg(cdRatios).toFixed(2)}}x`;
      document.getElementById("avgCdRatioDetail").textContent = "Large values here usually indicate a real drag blow-up or non-converged state.";
      document.getElementById("cdRange").textContent = `${{minCd.cd_ratio.toFixed(2)}}x .. ${{maxCd.cd_ratio.toFixed(2)}}x`;
      document.getElementById("cdRangeDetail").textContent = `Min at ${{minCd.alpha.toFixed(0)}} deg, max at ${{maxCd.alpha.toFixed(0)}} deg`;
      document.getElementById("pointCount").textContent = `${{rustConverged}} / ${{points.length}}`;
      document.getElementById("pointCountDetail").textContent = "RustFoil converged points in this run";
      caseNote.innerHTML = `<strong>${{fmtCaseLabel(item)}}</strong><br>Standard post-fix full-sweep comparison. Dashed parity line at 1.0 in the CD-ratio panel.`;
    }}
    function renderPlots(item) {{
      const alpha = item.points.map(p => p.alpha);
      const xfoilCl = item.points.map(p => p.xfoil_cl);
      const rustfoilCl = item.points.map(p => p.rustfoil_cl);
      const xfoilCd = item.points.map(p => p.xfoil_cd);
      const rustfoilCd = item.points.map(p => p.rustfoil_cd);
      const clError = item.points.map(p => p.cl_error_pct);
      const cdRatio = item.points.map(p => p.cd_ratio);
      Plotly.newPlot("clPlot", [
        {{ x: alpha, y: xfoilCl, mode: "lines+markers", name: "XFOIL CL", line: {{ color: "#7dd3fc", width: 3 }}, marker: {{ size: 6 }} }},
        {{ x: alpha, y: rustfoilCl, mode: "lines+markers", name: "RustFoil CL", line: {{ color: "#f472b6", width: 3 }}, marker: {{ size: 6 }} }}
      ], {{ ...plotLayoutBase, title: {{ text: "CL vs Alpha", x: 0.02, xanchor: "left" }}, yaxis: {{ ...plotLayoutBase.yaxis, title: "CL" }} }}, plotConfig);
      Plotly.newPlot("cdPlot", [
        {{ x: alpha, y: xfoilCd, mode: "lines+markers", name: "XFOIL CD", line: {{ color: "#7dd3fc", width: 3 }}, marker: {{ size: 6 }} }},
        {{ x: alpha, y: rustfoilCd, mode: "lines+markers", name: "RustFoil CD", line: {{ color: "#f472b6", width: 3 }}, marker: {{ size: 6 }} }}
      ], {{ ...plotLayoutBase, title: {{ text: "CD vs Alpha", x: 0.02, xanchor: "left" }}, yaxis: {{ ...plotLayoutBase.yaxis, title: "CD", type: "log" }} }}, plotConfig);
      Plotly.newPlot("clErrPlot", [
        {{ x: alpha, y: clError, mode: "lines+markers", name: "CL Error %", line: {{ color: "#fbbf24", width: 3 }}, marker: {{ size: 6 }} }}
      ], {{ ...plotLayoutBase, title: {{ text: "CL Error %", x: 0.02, xanchor: "left" }}, yaxis: {{ ...plotLayoutBase.yaxis, title: "Absolute % error" }}, showlegend: false }}, plotConfig);
      Plotly.newPlot("cdRatioPlot", [
        {{ x: alpha, y: cdRatio, mode: "lines+markers", name: "CD Ratio", line: {{ color: "#34d399", width: 3 }}, marker: {{ size: 6 }} }},
        {{ x: [alpha[0], alpha[alpha.length - 1]], y: [1, 1], mode: "lines", name: "Parity", line: {{ color: "#fb7185", width: 2, dash: "dash" }} }}
      ], {{ ...plotLayoutBase, title: {{ text: "RustFoil / XFOIL CD Ratio", x: 0.02, xanchor: "left" }}, yaxis: {{ ...plotLayoutBase.yaxis, title: "CD ratio", type: "log" }} }}, plotConfig);
    }}
    function renderCase(index) {{
      const item = DATA[index];
      updateCards(item);
      renderPlots(item);
    }}
    caseSelect.addEventListener("change", (event) => renderCase(Number(event.target.value)));
    if (DATA.length > 0) {{
      renderCase(0);
    }}
    if (isMobile()) {{
      const dots = document.querySelectorAll(".dot-row .dot");
      const cards = document.querySelectorAll(".summary .card");
      const observer = new IntersectionObserver((entries) => {{
        entries.forEach(entry => {{
          if (entry.isIntersecting) {{
            const idx = [...cards].indexOf(entry.target);
            dots.forEach((d, i) => d.classList.toggle("active", i === idx));
          }}
        }});
      }}, {{ root: document.querySelector(".summary"), threshold: 0.6 }});
      cards.forEach(card => observer.observe(card));
    }}
    document.querySelectorAll(".plot-tab").forEach(tab => {{
      tab.addEventListener("click", () => {{
        document.querySelectorAll(".plot-tab").forEach(t => t.classList.remove("active"));
        document.querySelectorAll(".plot-card").forEach(c => c.classList.remove("active"));
        tab.classList.add("active");
        const target = tab.dataset.target;
        document.querySelector(`.plot-card[data-plot="${{target}}"]`).classList.add("active");
        Plotly.Plots.resize(document.getElementById(target));
      }});
    }});
  </script>
</body>
</html>
"""
    output_path.write_text(html)


def resolve_paths(workspace: Path, tag: str) -> SweepPaths:
    xfoil_dir = workspace / "traces" / f"xfoil_{tag}"
    rustfoil_dir = workspace / "traces" / f"rustfoil_faithful_{tag}"
    comparison_dir = workspace / "comparison_results" / tag
    return SweepPaths(
        xfoil_dir=xfoil_dir,
        rustfoil_dir=rustfoil_dir,
        comparison_dir=comparison_dir,
        summary_json=comparison_dir / "comparison_summary.json",
        plot_html=comparison_dir / "full_sweep_plot.html",
        markdown_report=comparison_dir / "README.md",
    )


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run the standard full RustFoil/XFOIL comparison sweep")
    parser.add_argument(
        "--tag",
        default=datetime.now().strftime("rebuilt_%Y%m%d_%H%M%S"),
        help="Tag used for output directories (default: timestamped rebuilt_*)",
    )
    parser.add_argument(
        "--skip-xfoil",
        action="store_true",
        help="Reuse an existing XFOIL output directory for this tag",
    )
    parser.add_argument(
        "--skip-rustfoil",
        action="store_true",
        help="Reuse an existing RustFoil output directory for this tag",
    )
    parser.add_argument(
        "--alpha-start",
        type=float,
        default=-15.0,
        help="Starting alpha for both generators",
    )
    parser.add_argument(
        "--alpha-end",
        type=float,
        default=25.0,
        help="Ending alpha for both generators",
    )
    parser.add_argument(
        "--alpha-step",
        type=float,
        default=1.0,
        help="Alpha step for both generators",
    )
    parser.add_argument(
        "--jobs",
        type=int,
        default=4,
        help="Parallel jobs for each generator (default: 4)",
    )
    parser.add_argument(
        "--full-debug-traces",
        action="store_true",
        help="Persist full per-alpha debug payloads instead of compact summary-only traces",
    )
    return parser.parse_args()


def main() -> int:
    args = parse_args()
    workspace = Path(__file__).parent.parent
    paths = resolve_paths(workspace, args.tag)
    paths.comparison_dir.mkdir(parents=True, exist_ok=True)
    rust_solver_path = "faithful-xfoil-side"
    cases = case_tuples(workspace)

    print("=" * 72)
    print(" FULL SWEEP COMPARISON")
    print("=" * 72)
    print(f"  Tag:            {args.tag}")
    print(f"  XFOIL traces:   {paths.xfoil_dir}")
    print(f"  RustFoil traces:{paths.rustfoil_dir}")
    print(f"  Comparison dir: {paths.comparison_dir}")
    print(f"  Rust path:      {rust_solver_path}")
    print(f"  Trace mode:     {'full-debug' if args.full_debug_traces else 'summary-only'}")

    xfoil_command = [
        sys.executable,
        "scripts/generate_xfoil_traces.py",
        "--output-dir",
        str(paths.xfoil_dir.relative_to(workspace)),
        "--alpha-start",
        str(args.alpha_start),
        "--alpha-end",
        str(args.alpha_end),
        "--alpha-step",
        str(args.alpha_step),
        "--jobs",
        str(args.jobs),
        *(["--summary-only"] if not args.full_debug_traces else []),
    ]
    rustfoil_command = [
        sys.executable,
        "scripts/generate_rustfoil_traces.py",
        "--output-dir",
        str(paths.rustfoil_dir.relative_to(workspace)),
        "--alpha-start",
        str(args.alpha_start),
        "--alpha-end",
        str(args.alpha_end),
        "--alpha-step",
        str(args.alpha_step),
        "--jobs",
        str(args.jobs),
        *(["--summary-only"] if not args.full_debug_traces else []),
    ]

    if args.skip_xfoil and not paths.xfoil_dir.exists():
        raise SystemExit(f"Missing reused XFOIL directory: {paths.xfoil_dir}")
    if args.skip_rustfoil and not paths.rustfoil_dir.exists():
        raise SystemExit(f"Missing reused RustFoil directory: {paths.rustfoil_dir}")

    commands_to_run: list[tuple[str, list[str]]] = []
    if not args.skip_xfoil:
        commands_to_run.append(("xfoil", xfoil_command))
    if not args.skip_rustfoil:
        commands_to_run.append(("rustfoil-faithful", rustfoil_command))

    if len(commands_to_run) == 1:
        _, command = commands_to_run[0]
        run_command(command, workspace)
    elif len(commands_to_run) == 2:
        run_commands_parallel(commands_to_run, workspace)

    case_entries = []
    for foil, re_dir in cases:
        x_summary = paths.xfoil_dir / foil / re_dir / "summary.json"
        r_summary = paths.rustfoil_dir / foil / re_dir / "summary.json"
        if not x_summary.exists():
            case_entries.append({"foil": foil, "re": re_dir, "points": [], "error": "Missing XFOIL summary"})
            continue
        if not r_summary.exists():
            case_entries.append({"foil": foil, "re": re_dir, "points": [], "error": "Missing RustFoil summary"})
            continue
        try:
            points = build_case_points(
                paths.xfoil_dir,
                paths.rustfoil_dir,
                foil,
                re_dir,
                rust_solver_path,
            )
            case_entries.append({"foil": foil, "re": re_dir, "points": points})
        except Exception as exc:
            case_entries.append({"foil": foil, "re": re_dir, "points": [], "error": str(exc)})
    summary = build_summary(case_entries)
    summary["rust_solver_path"] = rust_solver_path

    with open(paths.summary_json, "w") as f:
        json.dump(summary, f, indent=2)
    write_markdown_report(summary, paths.markdown_report)
    write_plot_html(case_entries, paths.plot_html, rust_solver_path)

    overall = summary["overall"]
    print()
    print("=" * 72)
    print(" RESULTS")
    print("=" * 72)
    print(f"  Total points:              {overall['point_count']}")
    print(f"  RustFoil converged points: {overall['rustfoil_converged_points']}")
    print(f"  Avg CL error:              {overall['avg_cl_error_pct']:.2f}%")
    print(f"  Max CL error:              {overall['max_cl_error_pct']:.2f}%")
    print(f"  Avg CD ratio:              {overall['avg_cd_ratio']:.2f}x")
    print(f"  CD ratio range:            {overall['min_cd_ratio']:.2f}x .. {overall['max_cd_ratio']:.2f}x")
    print()
    print(f"  Summary JSON: {paths.summary_json}")
    print(f"  Plot HTML:    {paths.plot_html}")
    print(f"  Report:       {paths.markdown_report}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
