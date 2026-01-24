#!/usr/bin/env python3
"""
Compare TRCHEK2 intermediate and final states between XFOIL and RustFoil.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path


def load_events(path: Path):
    with open(path) as f:
        data = json.load(f)
    return data.get("events", data)


def index_last(events, key_fields, iter_field):
    indexed = {}
    for ev in events:
        key = tuple(ev.get(k) for k in key_fields)
        if key not in indexed or ev.get(iter_field, 0) >= indexed[key].get(iter_field, 0):
            indexed[key] = ev
    return indexed


def pct_diff(a, b):
    if a is None or b is None:
        return None
    if abs(a) < 1e-15:
        return 0.0 if abs(b) < 1e-15 else 100.0
    return abs(a - b) / abs(a) * 100.0


def main():
    parser = argparse.ArgumentParser(description="Compare TRCHEK2 intermediate/final dumps")
    parser.add_argument("xfoil_json", help="Path to XFOIL debug JSON")
    parser.add_argument("rust_json", help="Path to RustFoil debug JSON")
    parser.add_argument("--side", type=int, default=None, help="Filter by side (1/2)")
    parser.add_argument("--limit", type=int, default=50, help="Max stations to print")
    parser.add_argument("--window", type=int, default=2, help="MRCHUE window radius")
    args = parser.parse_args()

    xf_events = load_events(Path(args.xfoil_json))
    rf_events = load_events(Path(args.rust_json))

    xf_iter = [e for e in xf_events if e.get("subroutine") == "TRCHEK2_ITER"]
    xf_final = [e for e in xf_events if e.get("subroutine") == "TRCHEK2_FINAL"]
    xf_mrchue = [e for e in xf_events if e.get("subroutine") == "MRCHUE"]
    rf_iter = [e for e in rf_events if e.get("subroutine") == "TRCHEK2_ITER"]
    rf_final = [e for e in rf_events if e.get("subroutine") == "TRCHEK2_FINAL"]
    rf_mrchue = [e for e in rf_events if e.get("subroutine") == "MRCHUE"]

    xf_iter_map = index_last(xf_iter, ["side", "ibl"], "trchek_iter")
    rf_iter_map = index_last(rf_iter, ["side", "ibl"], "trchek_iter")
    xf_final_map = {(e.get("side"), e.get("ibl")): e for e in xf_final}
    rf_final_map = {(e.get("side"), e.get("ibl")): e for e in rf_final}

    keys = sorted(set(xf_iter_map.keys()) | set(rf_iter_map.keys()))
    if args.side is not None:
        keys = [k for k in keys if k[0] == args.side]

    print("=" * 90)
    print("TRCHEK2 INTERMEDIATE (last iter per station)")
    print("=" * 90)
    header = "Side IBL  XT_xf      XT_rf      TT_xf      TT_rf    ErrTT%  WF2_xf  WF2_rf"
    print(header)
    count = 0
    for key in keys:
        if count >= args.limit:
            break
        xf = xf_iter_map.get(key)
        rf = rf_iter_map.get(key)
        if not xf or not rf:
            continue

        tt_xf = None
        if xf.get("T1") is not None and xf.get("T2") is not None:
            tt_xf = xf.get("T1") * xf.get("wf1") + xf.get("T2") * xf.get("wf2")
        tt_rf = None
        if rf.get("T1") is not None and rf.get("T2") is not None:
            tt_rf = rf.get("T1") * rf.get("wf1") + rf.get("T2") * rf.get("wf2")

        err_tt = pct_diff(tt_xf, tt_rf)
        print(
            f"{key[0]:4d} {key[1]:3d} "
            f"{xf.get('xt', 0.0):10.6f} {rf.get('xt', 0.0):10.6f} "
            f"{(tt_xf or 0.0):10.6e} {(tt_rf or 0.0):10.6e} "
            f"{(err_tt or 0.0):7.2f} "
            f"{xf.get('wf2', 0.0):7.3f} {rf.get('wf2', 0.0):7.3f}"
        )
        count += 1

    # Detailed comparison for matching stations
    print("\nDetailed TRCHEK2_ITER inputs (matching stations):")
    detail_keys = [k for k in keys if k in xf_iter_map and k in rf_iter_map]
    for key in detail_keys[: min(len(detail_keys), args.limit)]:
        xf = xf_iter_map[key]
        rf = rf_iter_map[key]
        print(f"\n-- Side {key[0]} IBL {key[1]} --")
        for field in [
            "x1", "x2", "ampl1", "ampl2", "ax", "residual",
            "wf1", "wf2", "xt", "Hk1", "Hk2", "Rt1", "Rt2", "T1", "T2",
        ]:
            xf_val = xf.get(field)
            rf_val = rf.get(field)
            diff = pct_diff(xf_val, rf_val)
            if diff is None:
                diff_str = "n/a"
            else:
                diff_str = f"{diff:6.2f}%"
            print(f"  {field:8s} xf={xf_val!s:>14} rf={rf_val!s:>14}  diff={diff_str}")

    print("\n" + "=" * 90)
    print("MRCHUE station window comparison (output state)")
    print("=" * 90)
    xf_mrchue_map = {(e.get("side"), e.get("ibl")): e for e in xf_mrchue}
    rf_mrchue_map = {(e.get("side"), e.get("ibl")): e for e in rf_mrchue}

    for key in detail_keys[: min(len(detail_keys), args.limit)]:
        side, ibl = key
        print(f"\n-- Side {side} IBL {ibl} (window ±{args.window}) --")
        print("  IBL    x_xf       x_rf      theta_xf   theta_rf   err%   Hk_xf   Hk_rf   err%")
        for j in range(ibl - args.window, ibl + args.window + 1):
            xf = xf_mrchue_map.get((side, j), {})
            rf = rf_mrchue_map.get((side, j), {})
            x_xf = xf.get("x")
            x_rf = rf.get("x")
            th_xf = xf.get("theta")
            th_rf = rf.get("theta")
            hk_xf = xf.get("Hk")
            hk_rf = rf.get("Hk")
            err_th = pct_diff(th_xf, th_rf)
            err_hk = pct_diff(hk_xf, hk_rf)
            print(
                f"  {j:3d} "
                f"{(x_xf or 0.0):10.6f} {(x_rf or 0.0):10.6f} "
                f"{(th_xf or 0.0):10.6e} {(th_rf or 0.0):10.6e} "
                f"{(err_th or 0.0):6.2f} "
                f"{(hk_xf or 0.0):8.4f} {(hk_rf or 0.0):8.4f} "
                f"{(err_hk or 0.0):6.2f}"
            )

    print("\n" + "=" * 90)
    print("TRCHEK2 FINAL (transition stations only)")
    print("=" * 90)
    final_keys = sorted(set(xf_final_map.keys()) | set(rf_final_map.keys()))
    if args.side is not None:
        final_keys = [k for k in final_keys if k[0] == args.side]

    print("Side IBL  XT_xf      XT_rf      AX_xf      AX_rf    ErrAX%  Niter_xf Niter_rf")
    for key in final_keys:
        xf = xf_final_map.get(key)
        rf = rf_final_map.get(key)
        if not xf or not rf:
            continue
        err_ax = pct_diff(xf.get("ax_final"), rf.get("ax_final"))
        print(
            f"{key[0]:4d} {key[1]:3d} "
            f"{xf.get('xt_final', 0.0):10.6f} {rf.get('xt_final', 0.0):10.6f} "
            f"{xf.get('ax_final', 0.0):10.6e} {rf.get('ax_final', 0.0):10.6e} "
            f"{(err_ax or 0.0):7.2f} "
            f"{xf.get('n_iterations', 0):9d} {rf.get('n_iterations', 0):9d}"
        )

    print("\nRust-only fields at transition (first 5):")
    rust_only = [e for e in rf_final if e.get("transition")]
    for ev in rust_only[:5]:
        print(
            f"  side={ev.get('side')} ibl={ev.get('ibl')} "
            f"tt={ev.get('tt')} dt={ev.get('dt')} ut={ev.get('ut')} "
            f"Hk_t={ev.get('Hk_t')} Rt_t={ev.get('Rt_t')} St={ev.get('St')} Cq_t={ev.get('Cq_t')}"
        )


if __name__ == "__main__":
    main()
