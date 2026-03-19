"""
Compare XFOIL 6.99 (Fortran) vs Rustfoil for NACA 0012 with flap deflections.

Runs alpha sweeps at a range of control-surface deflections (TE-up and TE-down)
and prints a side-by-side convergence/accuracy comparison.
"""

import subprocess
import os
import sys
import re
import tempfile

# ---------------------------------------------------------------------------
# Config
# ---------------------------------------------------------------------------
XFOIL_BIN = os.path.join(os.path.dirname(__file__), "..", "Xfoil", "bin", "xfoil")
RE = 1e6
NCRIT = 9.0
MACH = 0.0
HINGE_X = 0.75
ALPHAS = (-4, 12, 1.0)
DEFLECTIONS = [-10, -5, 0, 5, 10, 15, 20]

# ---------------------------------------------------------------------------
# XFOIL driver (stdin piping)
# ---------------------------------------------------------------------------

def run_xfoil(deflection_deg: float, alpha_start: float, alpha_end: float,
              alpha_step: float) -> list[dict]:
    """Drive XFOIL via subprocess, return list of {alpha, cl, cd, cm}."""
    with tempfile.TemporaryDirectory(dir="/tmp") as tmpdir:
        polar_file = os.path.join(tmpdir, "p.dat")

        cmds = []
        cmds.append("PLOP")
        cmds.append("G")                   # toggle graphics off
        cmds.append("")                    # exit PLOP
        cmds.append("NACA 0012")
        if abs(deflection_deg) > 0.01:
            cmds.append("GDES")
            cmds.append("FLAP")
            cmds.append(f"{HINGE_X}")       # x/c hinge
            cmds.append("0")               # y/c hinge (midline)
            cmds.append(f"{deflection_deg}")
            cmds.append("X")               # execute
            cmds.append("")                 # back to main
            cmds.append("PANE")            # repanel
        cmds.append("OPER")
        cmds.append(f"VISC {RE:.0f}")
        cmds.append(f"MACH {MACH}")
        cmds.append(f"VPAR")
        cmds.append(f"N {NCRIT}")
        cmds.append("")
        cmds.append(f"PACC")
        cmds.append(polar_file)
        cmds.append("")                    # no dump file
        cmds.append(f"ASEQ {alpha_start} {alpha_end} {alpha_step}")
        cmds.append("")
        cmds.append("QUIT")

        stdin_text = "\n".join(cmds) + "\n"

        env = {k: v for k, v in os.environ.items()}
        env["DISPLAY"] = ""

        try:
            proc = subprocess.run(
                [XFOIL_BIN],
                input=stdin_text,
                capture_output=True,
                text=True,
                timeout=120,
                env=env,
            )
        except subprocess.TimeoutExpired:
            return []

        results = []
        if os.path.exists(polar_file):
            with open(polar_file) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("-") or line.startswith("alpha"):
                        continue
                    parts = line.split()
                    if len(parts) >= 4:
                        try:
                            alpha = float(parts[0])
                            cl = float(parts[1])
                            cd = float(parts[2])
                            cm = float(parts[4]) if len(parts) > 4 else 0.0
                            results.append(dict(alpha=alpha, cl=cl, cd=cd, cm=cm))
                        except (ValueError, IndexError):
                            continue
        return results


# ---------------------------------------------------------------------------
# Rustfoil driver
# ---------------------------------------------------------------------------

def run_rustfoil(deflection_deg: float, alpha_start: float, alpha_end: float,
                 alpha_step: float) -> list[dict]:
    """Run flexfoil (Rustfoil) and return list of {alpha, cl, cd, cm}."""
    import flexfoil

    foil = flexfoil.naca("0012")
    if abs(deflection_deg) > 0.01:
        foil = foil.with_flap(hinge_x=HINGE_X, deflection=deflection_deg)

    polar = foil.polar(
        alpha=(alpha_start, alpha_end, alpha_step),
        Re=RE,
        mach=MACH,
        ncrit=NCRIT,
    )

    results = []
    for r in polar.results:
        entry = dict(alpha=r.alpha, cl=r.cl, cd=r.cd, cm=r.cm,
                     converged=r.converged)
        results.append(entry)
    return results


# ---------------------------------------------------------------------------
# Comparison
# ---------------------------------------------------------------------------

def compare(deflection: float) -> tuple[int, int, int]:
    a_start, a_end, a_step = ALPHAS

    xfoil_pts = run_xfoil(deflection, a_start, a_end, a_step)
    rust_pts = run_rustfoil(deflection, a_start, a_end, a_step)

    xfoil_map = {round(p["alpha"], 1): p for p in xfoil_pts}
    rust_map = {round(p["alpha"], 1): p for p in rust_pts}

    all_alphas = sorted(set(xfoil_map.keys()) | set(rust_map.keys()))

    xfoil_conv = len(xfoil_pts)
    rust_conv = sum(1 for p in rust_pts if p.get("converged", True))
    total = len(all_alphas)

    print(f"\n{'='*78}")
    print(f"  NACA 0012  |  Flap δ = {deflection:+.0f}°  |  Re = {RE:.0e}  |  hinge x/c = {HINGE_X}")
    print(f"  XFOIL converged: {xfoil_conv}/{total}   Rustfoil converged: {rust_conv}/{total}")
    print(f"{'='*78}")
    print(f"  {'α':>5}  │ {'CL_xf':>8} {'CL_rf':>8} {'ΔCL':>7} │ "
          f"{'CD_xf':>8} {'CD_rf':>8} {'ΔCD':>7} │ "
          f"{'CM_xf':>8} {'CM_rf':>8} │ {'note':>6}")
    print(f"  {'─'*5}──┼─{'─'*8}─{'─'*8}─{'─'*7}─┼─"
          f"{'─'*8}─{'─'*8}─{'─'*7}─┼─"
          f"{'─'*8}─{'─'*8}─┼─{'─'*6}")

    cl_diffs = []
    cd_diffs = []

    for a in all_alphas:
        xp = xfoil_map.get(a)
        rp = rust_map.get(a)
        note = ""

        if xp and rp and rp.get("converged", True):
            dcl = rp["cl"] - xp["cl"]
            dcd = rp["cd"] - xp["cd"]
            cl_diffs.append(abs(dcl))
            cd_diffs.append(abs(dcd))
            print(f"  {a:5.1f}  │ {xp['cl']:8.4f} {rp['cl']:8.4f} {dcl:+7.4f} │ "
                  f"{xp['cd']:8.5f} {rp['cd']:8.5f} {dcd:+7.5f} │ "
                  f"{xp['cm']:8.4f} {rp['cm']:8.4f} │ {'':>6}")
        elif xp and (not rp or not rp.get("converged", True)):
            note = "RF ✗"
            cl_rf = f"{rp['cl']:8.4f}" if rp else "     N/A"
            cd_rf = f"{rp['cd']:8.5f}" if rp else "      N/A"
            cm_rf = f"{rp['cm']:8.4f}" if rp else "     N/A"
            print(f"  {a:5.1f}  │ {xp['cl']:8.4f} {cl_rf} {'':>7} │ "
                  f"{xp['cd']:8.5f} {cd_rf} {'':>7} │ "
                  f"{xp['cm']:8.4f} {cm_rf} │ {note:>6}")
        elif rp and rp.get("converged", True) and not xp:
            note = "XF ✗"
            print(f"  {a:5.1f}  │ {'     N/A':>8} {rp['cl']:8.4f} {'':>7} │ "
                  f"{'      N/A':>8} {rp['cd']:8.5f} {'':>7} │ "
                  f"{'     N/A':>8} {rp['cm']:8.4f} │ {note:>6}")
        else:
            note = "both✗"
            print(f"  {a:5.1f}  │ {'     N/A':>8} {'     N/A':>8} {'':>7} │ "
                  f"{'      N/A':>8} {'      N/A':>8} {'':>7} │ "
                  f"{'     N/A':>8} {'     N/A':>8} │ {note:>6}")

    if cl_diffs:
        avg_dcl = sum(cl_diffs) / len(cl_diffs)
        max_dcl = max(cl_diffs)
        avg_dcd = sum(cd_diffs) / len(cd_diffs)
        max_dcd = max(cd_diffs)
        print(f"  {'─'*74}")
        print(f"  ΔCL  avg={avg_dcl:.4f}  max={max_dcl:.4f}   │   "
              f"ΔCD  avg={avg_dcd:.5f}  max={max_dcd:.5f}")

    return xfoil_conv, rust_conv, total


def main():
    print("╔══════════════════════════════════════════════════════════════════╗")
    print("║  XFOIL 6.99 vs Rustfoil — NACA 0012 Flap Deflection Comparison ║")
    print("╚══════════════════════════════════════════════════════════════════╝")

    summary = []
    for defl in DEFLECTIONS:
        xc, rc, total = compare(defl)
        summary.append((defl, xc, rc, total))

    print("\n\n" + "="*78)
    print("  CONVERGENCE SUMMARY")
    print("="*78)
    print(f"  {'δ (°)':>6} │ {'XFOIL conv':>12} │ {'Rustfoil conv':>14} │ {'match':>6}")
    print(f"  {'─'*6}──┼─{'─'*12}──┼─{'─'*14}──┼─{'─'*6}")

    for defl, xc, rc, total in summary:
        print(f"  {defl:+6.0f} │ {xc:>5}/{total:<5} │ {rc:>7}/{total:<5}  │ "
              f"{'OK' if xc == rc else f'Δ={rc-xc:+d}'}")


if __name__ == "__main__":
    main()
