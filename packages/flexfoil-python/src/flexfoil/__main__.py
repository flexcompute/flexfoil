"""CLI entry point: ``python -m flexfoil`` or ``flexfoil`` (after pip install)."""

from __future__ import annotations

import argparse
import sys


def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        prog="flexfoil",
        description="Airfoil analysis — XFOIL-faithful solver with a local web UI",
    )
    sub = parser.add_subparsers(dest="command")

    # -- serve --
    serve_p = sub.add_parser("serve", help="Launch the local web UI and API server")
    serve_p.add_argument("--port", type=int, default=8420)
    serve_p.add_argument("--host", default="127.0.0.1")
    serve_p.add_argument("--no-browser", action="store_true")
    serve_p.add_argument("--db", default=None, help="Path to SQLite database")

    # -- info --
    sub.add_parser("info", help="Show flexfoil configuration info")

    # -- solve --
    solve_p = sub.add_parser("solve", help="Quick single-point solve from the CLI")
    solve_p.add_argument("airfoil", help="NACA designation or .dat file path")
    solve_p.add_argument("-a", "--alpha", type=float, default=0.0)
    solve_p.add_argument("-r", "--Re", type=float, default=1e6)
    solve_p.add_argument("-m", "--mach", type=float, default=0.0)

    args = parser.parse_args(argv)

    if args.command == "serve":
        from flexfoil import serve
        serve(
            port=args.port,
            host=args.host,
            open_browser=not args.no_browser,
            db_path=args.db,
        )

    elif args.command == "info":
        _print_info()

    elif args.command == "solve":
        _run_solve(args)

    else:
        parser.print_help()


def _print_info() -> None:
    from flexfoil.database import default_db_path

    print("flexfoil v1.1.3")
    print(f"  Database: {default_db_path()}")
    try:
        from flexfoil._rustfoil import analyze_faithful
        print("  Solver:   rustfoil (native Rust via PyO3)")
    except ImportError:
        print("  Solver:   NOT AVAILABLE (Rust extension not built)")


def _run_solve(args) -> None:
    import flexfoil

    name: str = args.airfoil
    if name.endswith(".dat"):
        foil = flexfoil.load(name)
    else:
        foil = flexfoil.naca(name)

    result = foil.solve(alpha=args.alpha, Re=args.Re, mach=args.mach)
    print(result)

    if result.success:
        print(f"  L/D = {result.ld:.2f}" if result.ld else "  L/D = N/A")


if __name__ == "__main__":
    main()
