"""Run the Flow360 solve and extract results.

NOTE: Flow360Solver is a GPU/MPI binary. It must run in a normal (non-sandboxed)
shell — launch the pipeline with ``--solve`` from an interactive session, not from
inside a restricted agent sandbox (the solver gets killed there). Meshing and case
construction have no such restriction.

The solver requires a columnar-data postprocessor co-process that owns the IPC
socket (``ipc_control.sock``); this module launches it, waits for the socket, runs
the solver under ``mpirun -np 1``, then shuts the postprocessor down.
"""
from __future__ import annotations

import csv
import subprocess
import time
from pathlib import Path


def run_solver(workdir: str | Path, find, env: dict, *,
               gpu: int = 0, mpirun: str = "/usr/bin/mpirun", timeout: int = 7200) -> None:
    """``find`` resolves a tool name to its absolute path (see rans.env.make_env)."""
    workdir = Path(workdir)
    senv = dict(env)
    senv["CUDA_VISIBLE_DEVICES"] = str(gpu)
    senv["OMP_NUM_THREADS"] = "1"

    sock = workdir / "ipc_control.sock"
    if sock.exists():
        sock.unlink()

    pp_log = open(workdir / "postprocessor.log", "w")
    pp = subprocess.Popen(
        [find("columnarDataProcessor.py"), "--asyncMode",
         "--inputSimulationJson", "simulation.json",
         "--columnarDataProcessorJson", "columnar.json"],
        cwd=str(workdir), env=senv, stdout=pp_log, stderr=subprocess.STDOUT)
    try:
        for _ in range(40):
            if sock.exists():
                break
            time.sleep(0.5)
        with open(workdir / "solver.log", "w") as slog:
            subprocess.run([mpirun, "-np", "1", "--bind-to", "none",
                            find("Flow360Solver")],
                           cwd=str(workdir), env=senv, stdout=slog,
                           stderr=subprocess.STDOUT, check=True, timeout=timeout)
    finally:
        pp.terminate()
        try:
            pp.wait(timeout=10)
        except subprocess.TimeoutExpired:
            pp.kill()
        pp_log.close()


def extract_forces(workdir: str | Path) -> dict:
    """Final CL/CD/(L/D) from total_forces_v2.csv + a listing of result files."""
    workdir = Path(workdir)
    rows = list(csv.reader(open(workdir / "total_forces_v2.csv")))
    hdr = [h.strip() for h in rows[0]]
    iCL, iCD = hdr.index("CL"), hdr.index("CD")
    data = [[c.strip() for c in r if c.strip() != ""] for r in rows[1:] if len(r) > 6]
    last = data[-1]
    CL, CD = float(last[iCL]), float(last[iCD])
    return {
        "step": int(last[1]),
        "CL": CL,
        "CD": CD,
        "L_over_D": CL / CD if CD else float("nan"),
        "paraview": {
            "volume": str(workdir / "volume.pvtu"),
            "surfaces": sorted(str(p) for p in workdir.glob("surface_fluid_*.pvtu")),
            "center_slice": str(workdir / "slice_centerSpan.pvtu"),
        },
    }
