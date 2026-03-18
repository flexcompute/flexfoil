"""Tests for the high-level Python API (Airfoil, solve, polar)."""

import pytest


@pytest.fixture(autouse=True)
def _use_tmp_db(tmp_path, monkeypatch):
    """Direct the database to a temporary directory for test isolation."""
    monkeypatch.setenv("FLEXFOIL_DATA_DIR", str(tmp_path))
    import flexfoil.database
    flexfoil.database._default_db = None


# ---------------------------------------------------------------------------
# Airfoil creation
# ---------------------------------------------------------------------------

class TestAirfoilCreation:
    def test_from_naca(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        assert foil.name == "NACA 2412"
        assert foil.n_panels == 160
        assert len(foil.panel_coords) == 160

    def test_from_naca_short(self):
        import flexfoil
        foil = flexfoil.naca("12")
        assert foil.name == "NACA 0012"
        assert foil.n_panels == 160

    def test_from_naca_custom_panels(self):
        import flexfoil
        foil = flexfoil.naca("2412", n_panels=80)
        assert foil.n_panels == 80

    def test_hash_deterministic(self):
        import flexfoil
        f1 = flexfoil.naca("2412")
        f2 = flexfoil.naca("2412")
        assert f1.hash == f2.hash

    def test_hash_differs_across_foils(self):
        import flexfoil
        f1 = flexfoil.naca("2412")
        f2 = flexfoil.naca("0012")
        assert f1.hash != f2.hash

    def test_from_coordinates(self):
        import flexfoil
        base = flexfoil.naca("0012")
        x = [p[0] for p in base.raw_coords]
        y = [p[1] for p in base.raw_coords]
        foil = flexfoil.from_coordinates(x, y, name="custom")
        assert foil.name == "custom"
        assert foil.n_panels > 0

    def test_from_coordinates_mismatched_lengths(self):
        import flexfoil
        with pytest.raises(ValueError, match="same length"):
            flexfoil.from_coordinates([0, 1], [0], name="bad")

    def test_repr(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        assert "NACA 2412" in repr(foil)
        assert "n_panels=160" in repr(foil)


# ---------------------------------------------------------------------------
# Single-point solve
# ---------------------------------------------------------------------------

class TestSolve:
    def test_viscous_solve(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        result = foil.solve(alpha=5.0, Re=1e6)
        assert result.success
        assert result.converged
        assert result.cl > 0.5
        assert 0 < result.cd < 0.1

    def test_inviscid_solve(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        result = foil.solve(alpha=5.0, viscous=False)
        assert result.success
        assert result.cl > 0.5
        assert result.cd == 0.0

    def test_solve_result_ld(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        result = foil.solve(alpha=5.0, Re=1e6)
        assert result.ld is not None
        assert result.ld > 10

    def test_inviscid_ld_is_none(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        result = foil.solve(alpha=5.0, viscous=False)
        assert result.ld is None

    def test_solve_result_repr(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        result = foil.solve(alpha=5.0, Re=1e6)
        r = repr(result)
        assert "CL=" in r
        assert "CD=" in r
        assert "converged" in r

    def test_solve_at_zero_alpha_symmetric(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        result = foil.solve(alpha=0.0, Re=1e6)
        assert result.success
        assert abs(result.cl) < 0.01

    def test_solve_stores_in_db(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        foil.solve(alpha=0.0, Re=1e6)
        runs = flexfoil.runs()
        assert len(runs) >= 1

    def test_solve_store_false(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        foil.solve(alpha=0.0, Re=1e6, store=False)
        runs = flexfoil.runs()
        assert len(runs) == 0

    def test_solve_different_re(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        r1 = foil.solve(alpha=5.0, Re=1e5)
        r2 = foil.solve(alpha=5.0, Re=3e6)
        assert r1.cd > r2.cd  # lower Re => higher drag


# ---------------------------------------------------------------------------
# Polar sweep
# ---------------------------------------------------------------------------

class TestPolar:
    def test_polar_sweep(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 5, 2.5), Re=1e6)
        assert len(polar.results) == 3
        assert len(polar.converged) >= 2

    def test_polar_properties(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        assert len(polar.alpha) == len(polar.cl) == len(polar.cd) == len(polar.cm)

    def test_polar_to_dict(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        d = polar.to_dict()
        assert set(d.keys()) == {"alpha", "cl", "cd", "cm", "ld"}

    def test_polar_repr(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        r = repr(polar)
        assert "NACA 2412" in r
        assert "converged" in r

    def test_polar_inviscid(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6, viscous=False)
        assert all(cd == 0.0 for cd in polar.cd)

    def test_polar_plot_plotly(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        fig = polar.plot(show=False)
        assert fig is not None

    def test_polar_plot_matplotlib(self):
        import flexfoil
        import matplotlib
        matplotlib.use("Agg")
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        fig = polar.plot(show=False, backend="matplotlib")
        assert fig is not None


# ---------------------------------------------------------------------------
# Pandas integration
# ---------------------------------------------------------------------------

class TestPandas:
    def test_polar_to_dataframe(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        df = polar.to_dataframe()
        assert "cl" in df.columns
        assert "cd" in df.columns
        assert len(df) >= 2

    def test_runs_returns_dataframe(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        foil.solve(alpha=5.0, Re=1e6)
        df = flexfoil.runs()
        import pandas as pd
        assert isinstance(df, pd.DataFrame)
        assert len(df) >= 1


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

class TestCli:
    def test_info(self, capsys):
        from flexfoil.__main__ import main
        main(["info"])
        out = capsys.readouterr().out
        assert "flexfoil" in out
        assert "Database" in out

    def test_solve_naca(self, capsys):
        from flexfoil.__main__ import main
        main(["solve", "2412", "-a", "5", "-r", "1e6"])
        out = capsys.readouterr().out
        assert "CL=" in out

    def test_no_args_prints_help(self, capsys):
        from flexfoil.__main__ import main
        main([])
        out = capsys.readouterr().out
        assert "usage" in out.lower() or "flexfoil" in out.lower()


# ---------------------------------------------------------------------------
# Parallel polar
# ---------------------------------------------------------------------------

class TestParallelPolar:
    def test_parallel_default(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        polar = foil.polar(alpha=(0, 4, 2.0), Re=1e6)
        assert len(polar.converged) >= 2

    def test_parallel_matches_sequential(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        par = foil.polar(alpha=(0, 4, 2.0), Re=1e6, parallel=True, store=False)
        seq = foil.polar(alpha=(0, 4, 2.0), Re=1e6, parallel=False, store=False)
        assert len(par.converged) == len(seq.converged)
        for p, s in zip(par.converged, seq.converged):
            assert abs(p.cl - s.cl) < 1e-10
            assert abs(p.cd - s.cd) < 1e-10

    def test_parallel_inviscid(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        polar = foil.polar(alpha=(0, 6, 2.0), Re=1e6, viscous=False, parallel=True)
        assert all(cd == 0.0 for cd in polar.cd)

    def test_parallel_stores_in_db(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        foil.polar(alpha=(0, 2, 1.0), Re=1e6, parallel=True)
        runs = flexfoil.runs()
        assert len(runs) >= 3


# ---------------------------------------------------------------------------
# Flap deflection
# ---------------------------------------------------------------------------

class TestFlap:
    def test_with_flap_returns_new_airfoil(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        flapped = foil.with_flap(hinge_x=0.75, deflection=10)
        assert "flap" in flapped.name
        assert flapped.n_panels > 0
        assert flapped is not foil

    def test_flap_increases_cl(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        flapped = foil.with_flap(hinge_x=0.75, deflection=10)
        r_clean = foil.solve(alpha=5, Re=1e6, store=False)
        r_flap = flapped.solve(alpha=5, Re=1e6, store=False)
        assert r_flap.cl > r_clean.cl

    def test_zero_deflection_similar_to_clean(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        flapped = foil.with_flap(hinge_x=0.75, deflection=0.0)
        r_clean = foil.solve(alpha=3, Re=1e6, store=False)
        r_flap = flapped.solve(alpha=3, Re=1e6, store=False)
        assert abs(r_flap.cl - r_clean.cl) < 0.05

    def test_negative_deflection(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        flapped = foil.with_flap(hinge_x=0.75, deflection=-5)
        r_clean = foil.solve(alpha=5, Re=1e6, store=False)
        r_flap = flapped.solve(alpha=5, Re=1e6, store=False)
        assert r_flap.cl < r_clean.cl

    def test_flap_different_hinge(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        f60 = foil.with_flap(hinge_x=0.60, deflection=10)
        f80 = foil.with_flap(hinge_x=0.80, deflection=10)
        assert f60.hash != f80.hash

    def test_flap_polar(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        flapped = foil.with_flap(hinge_x=0.75, deflection=10)
        polar = flapped.polar(alpha=(0, 4, 2.0), Re=1e6)
        assert len(polar.converged) >= 2
        assert all(cl > 0 for cl in polar.cl)

    def test_flap_repr(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        flapped = foil.with_flap(hinge_x=0.75, deflection=10)
        r = repr(flapped)
        assert "flap" in r
        assert "NACA 2412" in r


# ---------------------------------------------------------------------------
# Matrix sweep (multi-Re)
# ---------------------------------------------------------------------------

class TestMatrixSweep:
    def test_multi_re_sweep(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        results = {}
        for Re in [2e5, 1e6]:
            polar = foil.polar(alpha=(0, 4, 2.0), Re=Re, store=False)
            results[Re] = polar
        assert results[2e5].cd[0] > results[1e6].cd[0]

    def test_flap_x_re_matrix(self):
        import flexfoil
        foil = flexfoil.naca("0012")
        flapped = foil.with_flap(hinge_x=0.75, deflection=5)
        for Re in [5e5, 1e6]:
            polar = flapped.polar(alpha=(0, 4, 2.0), Re=Re, store=False)
            assert len(polar.converged) >= 2

    def test_matrix_results_stored(self):
        import flexfoil
        foil = flexfoil.naca("2412")
        for Re in [5e5, 1e6]:
            foil.polar(alpha=(0, 2, 1.0), Re=Re)
        runs = flexfoil.runs()
        assert len(runs) >= 6
