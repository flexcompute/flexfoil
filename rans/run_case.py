#!/usr/bin/env python3
"""CLI for the quasi-2D multi-element RANS pipeline.

    # build mesh + Flow360 case (runs anywhere):
    python run_case.py examples/highlift_deploy012.json -o out/

    # also solve (run from a normal, non-sandboxed interactive shell):
    python run_case.py examples/highlift_deploy012.json -o out/ --solve

Run with the compute venv's python so the flow360 SDK is importable, e.g.
    /home/qiqi/flexcompute/compute/.venv/bin/python run_case.py ...
"""
import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from rans.pipeline import run  # noqa: E402


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("config", help="contours JSON describing the case")
    ap.add_argument("-o", "--out", required=True, help="output directory")
    ap.add_argument("--solve", action="store_true",
                    help="run the GPU solver (non-sandboxed shell only)")
    ap.add_argument("--compute-root", default=None,
                    help="Flow360 compute install root (else $FLOW360_COMPUTE_ROOT)")
    ap.add_argument("--gpu", type=int, default=0, help="CUDA device for the solve")
    ap.add_argument("--fast", action="store_true",
                    help="fast-iteration mode: cache the SDK case JSONs + skip auto-vis")
    args = ap.parse_args()

    # `fast` may also be set by a top-level "fast": true in the case JSON, so the
    # daemon (which calls run_case.py unchanged) can request it per-job.
    fast = args.fast or bool(json.loads(Path(args.config).read_text()).get("fast", False))
    summary = run(args.config, args.out, solve=args.solve,
                  compute_root=args.compute_root, gpu=args.gpu, fast=fast)
    print(json.dumps(summary, indent=2))
    if "forces" in summary:
        f = summary["forces"]
        print(f"\nCL={f['CL']:.4f}  CD={f['CD']:.5f}  L/D={f['L_over_D']:.2f}  (step {f['step']})")
    else:
        print(f"\nCase built in {args.out}. Solve from an interactive shell with --solve.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
