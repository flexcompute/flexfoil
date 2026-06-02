"""Locate the Flow360 compute install and build the subprocess environment.

Override the install location with FLOW360_COMPUTE_ROOT (default matches this
machine). The compute venv ships the Python tools (convert/generate/preprocess,
the flow360 SDK) in its bin; the binaries (solver, mesh tools) live in
install/release/bin.
"""
from __future__ import annotations

import os
from pathlib import Path

DEFAULT_COMPUTE_ROOT = Path(
    os.environ.get("FLOW360_COMPUTE_ROOT", "/home/qiqi/flexcompute/compute"))


def resolve(compute_root: str | Path | None = None):
    root = Path(compute_root or DEFAULT_COMPUTE_ROOT)
    release = root / "install" / "release"
    venv = root / ".venv"
    return root, release, venv, release / "bin"


def _site_packages(release: Path) -> str:
    hits = sorted(release.glob("lib/python3.*/site-packages"))
    return str(hits[0]) if hits else str(release / "lib")


def make_env(compute_root: str | Path | None = None):
    """Return (env dict, find) for running the compute tools as subprocesses.

    The C++ binaries live in install/release/bin; the Python tools (convert/
    generate/preprocess, columnarDataProcessor) live in the venv's bin. ``find``
    resolves a tool name to its absolute path across both locations.
    """
    root, release, venv, bindir = resolve(compute_root)
    venv_bin = venv / "bin"
    env = dict(os.environ)
    env["VIRTUAL_ENV"] = str(venv)
    env["PATH"] = f"{bindir}:{venv_bin}:" + env.get("PATH", "")
    env["LD_LIBRARY_PATH"] = f"{release}/lib:" + env.get("LD_LIBRARY_PATH", "")
    env["PYTHONPATH"] = _site_packages(release) + ":" + env.get("PYTHONPATH", "")
    env["OMP_NUM_THREADS"] = "1"
    env["FLOW360_SUPPRESS_BETA_WARNING"] = "1"   # quiet the SDK beta banner in subprocesses

    def find(name: str) -> str:
        for d in (bindir, venv_bin):
            p = d / name
            if p.exists():
                return str(p)
        raise FileNotFoundError(f"tool {name!r} not found in {bindir} or {venv_bin}")

    return env, find
