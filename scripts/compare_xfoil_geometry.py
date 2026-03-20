"""
Diagnostic: Use XFOIL's own FLAP+PANE geometry, feed it to Rustfoil's solver.

This isolates whether convergence failures are due to:
  (a) Our flap geometry pipeline (deflect_flap + repanel_xfoil), or
  (b) Our BL solver itself

Procedure for each deflection:
  1. XFOIL: NACA 0012 → GDES FLAP → PANE → save paneled coords → run polar
  2. Rustfoil with OUR geometry: naca("0012").with_flap(...)  → run polar
  3. Rustfoil with XFOIL geometry: load XFOIL's saved coords → run polar
"""

import subprocess
import os
import sys
import tempfile

XFOIL_BIN = os.path.join(os.path.dirname(__file__), "..", "Xfoil", "bin", "xfoil")
RE = 1e6
NCRIT = 9.0
HINGE_X = 0.75
ALPHAS = (-4, 12, 1.0)
DEFLECTIONS = [-10, -5, 5, 10, 15, 20]


def run_xfoil_save_geometry_and_polar(deflection_deg, alpha_start, alpha_end, alpha_step):
    """Run XFOIL: deflect flap, save paneled coords, run polar. Returns (coords, polar_pts)."""
    polar_file = f"/tmp/xf_polar_{deflection_deg:+.0f}.dat"
    coord_file = f"/tmp/xf_coords_{deflection_deg:+.0f}.dat"

    for f in [polar_file, coord_file]:
        if os.path.exists(f):
            os.remove(f)

    cmds = ["PLOP", "G", "", "NACA 0012"]
    if abs(deflection_deg) > 0.01:
        cmds += ["GDES", "FLAP", f"{HINGE_X}", "0", f"{deflection_deg}", "X", "", "PANE"]
    cmds += [f"PSAV {coord_file}",
             "OPER", f"VISC {RE:.0f}", "VPAR", f"N {NCRIT}", "",
             "PACC", polar_file, "",
             f"ASEQ {alpha_start} {alpha_end} {alpha_step}", "", "QUIT"]

    env = dict(os.environ)
    env["DISPLAY"] = ""
    subprocess.run([XFOIL_BIN], input="\n".join(cmds) + "\n",
                   capture_output=True, text=True, timeout=120, env=env)

    coords = []
    if os.path.exists(coord_file):
        with open(coord_file) as f:
            for line in f:
                parts = line.split()
                if len(parts) == 2:
                    try:
                        coords.append((float(parts[0]), float(parts[1])))
                    except ValueError:
                        pass

    polar_pts = []
    if os.path.exists(polar_file):
        with open(polar_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("-") or line.startswith("a"):
                    continue
                parts = line.split()
                if len(parts) >= 5:
                    try:
                        polar_pts.append(dict(
                            alpha=float(parts[0]), cl=float(parts[1]),
                            cd=float(parts[2]), cm=float(parts[4])
                        ))
                    except (ValueError, IndexError):
                        pass

    return coords, polar_pts


def run_rustfoil_our_geometry(deflection_deg, alpha_start, alpha_end, alpha_step):
    import flexfoil
    foil = flexfoil.naca("0012")
    if abs(deflection_deg) > 0.01:
        foil = foil.with_flap(hinge_x=HINGE_X, deflection=deflection_deg)
    polar = foil.polar(alpha=(alpha_start, alpha_end, alpha_step), Re=RE, ncrit=NCRIT)
    return [(r.alpha, r.cl, r.cd, r.cm, r.converged) for r in polar.results]


def run_rustfoil_xfoil_geometry(coords, alpha_start, alpha_end, alpha_step):
    import flexfoil
    if not coords:
        return []
    x = [c[0] for c in coords]
    y = [c[1] for c in coords]
    foil = flexfoil.Airfoil.from_xy(x, y, name="XFOIL-paneled")
    polar = foil.polar(alpha=(alpha_start, alpha_end, alpha_step), Re=RE, ncrit=NCRIT)
    return [(r.alpha, r.cl, r.cd, r.cm, r.converged) for r in polar.results]


def main():
    a_start, a_end, a_step = ALPHAS
    total = int((a_end - a_start) / a_step) + 1

    print("╔════════════════════════════════════════════════════════════════════════╗")
    print("║  Geometry Pipeline Diagnostic: XFOIL geom vs Rustfoil geom vs Solver ║")
    print("╚════════════════════════════════════════════════════════════════════════╝")
    print()

    summary = []

    for defl in DEFLECTIONS:
        print(f"  Running δ = {defl:+d}° ...", end=" ", flush=True)

        xf_coords, xf_polar = run_xfoil_save_geometry_and_polar(defl, a_start, a_end, a_step)
        rf_ours = run_rustfoil_our_geometry(defl, a_start, a_end, a_step)
        rf_xfgeom = run_rustfoil_xfoil_geometry(xf_coords, a_start, a_end, a_step)

        xf_conv = len(xf_polar)
        rf_ours_conv = sum(1 for r in rf_ours if r[4])
        rf_xfgeom_conv = sum(1 for r in rf_xfgeom if r[4])

        summary.append((defl, xf_conv, rf_ours_conv, rf_xfgeom_conv, total,
                         len(xf_coords)))
        print(f"XFOIL={xf_conv}/{total}  RF(ours)={rf_ours_conv}/{total}  "
              f"RF(XF geom)={rf_xfgeom_conv}/{total}  [{len(xf_coords)} pts]")

    print()
    print("=" * 78)
    print("  CONVERGENCE SUMMARY")
    print("=" * 78)
    print(f"  {'δ':>4} │ {'XFOIL':>7} │ {'RF ours':>9} │ {'RF+XF geom':>11} │ {'diagnosis':>20}")
    print(f"  {'─'*4}──┼─{'─'*7}──┼─{'─'*9}──┼─{'─'*11}──┼─{'─'*20}")

    for defl, xf, rf_o, rf_xf, tot, npts in summary:
        if rf_xf >= xf - 1 and rf_o >= xf - 1:
            diag = "OK"
        elif rf_xf >= xf - 1 and rf_o < xf - 2:
            diag = "GEOMETRY PIPELINE"
        elif rf_xf < xf - 2:
            diag = "SOLVER"
        else:
            diag = "mixed"

        print(f"  {defl:+4d} │ {xf:>3}/{tot:<3} │ {rf_o:>5}/{tot:<3} │ {rf_xf:>7}/{tot:<3} │ {diag:>20}")


if __name__ == "__main__":
    main()
