"""Tests for the Rust solver bindings."""

import pytest
from flexfoil._rustfoil import (
    analyze_faithful,
    analyze_inviscid,
    generate_naca4,
    repanel_xfoil,
    parse_dat_file,
)


# ---------------------------------------------------------------------------
# NACA generation
# ---------------------------------------------------------------------------

class TestNacaGeneration:
    def test_generates_points(self):
        pts = generate_naca4(2412)
        assert len(pts) > 100
        assert all(isinstance(p, tuple) and len(p) == 2 for p in pts)

    def test_symmetric_naca0012(self):
        pts = generate_naca4(12)
        xs, ys = zip(*pts)
        max_y = max(abs(y) for y in ys)
        assert 0.04 < max_y < 0.08

    def test_naca0012_symmetric_cl(self):
        """Symmetric foil at alpha=0 should have CL ~ 0."""
        pts = generate_naca4(12)
        flat = [v for x, y in pts for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_inviscid(coords, 0.0)
        assert result["success"]
        assert abs(result["cl"]) < 0.01

    def test_different_designations(self):
        for des in [12, 2412, 4412, 23012]:
            pts = generate_naca4(des)
            assert len(pts) > 50, f"NACA {des} should generate points"

    def test_custom_nside(self):
        pts_default = generate_naca4(2412)
        pts_small = generate_naca4(2412, 40)
        assert len(pts_small) < len(pts_default)

    def test_leading_edge_near_origin(self):
        pts = generate_naca4(12)
        xs = [p[0] for p in pts]
        assert min(xs) < 0.01


# ---------------------------------------------------------------------------
# Repaneling
# ---------------------------------------------------------------------------

class TestRepanel:
    def test_repanel_160(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        assert len(paneled) == 160

    def test_repanel_80(self):
        raw = generate_naca4(12)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 80)
        assert len(paneled) == 80

    def test_repanel_custom_params(self):
        raw = generate_naca4(12)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 80, 1.0, 0.15, 0.667)
        assert len(paneled) == 80

    def test_repanel_preserves_chord(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        xs = [p[0] for p in paneled]
        assert max(xs) > 0.99
        assert min(xs) < 0.01

    def test_empty_input(self):
        assert repanel_xfoil([], 160) == []

    def test_too_few_points(self):
        assert repanel_xfoil([1.0, 0.0], 160) == []


# ---------------------------------------------------------------------------
# Inviscid solver
# ---------------------------------------------------------------------------

class TestInviscidSolver:
    def test_naca0012_zero_alpha(self):
        raw = generate_naca4(12)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_inviscid(coords, 0.0)
        assert result["success"]
        assert abs(result["cl"]) < 0.01

    def test_naca2412_positive_lift(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_inviscid(coords, 5.0)
        assert result["success"]
        assert result["cl"] > 0.5

    def test_returns_cp(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_inviscid(coords, 3.0)
        assert result["success"]
        assert len(result["cp"]) > 0
        assert len(result["cp_x"]) > 0

    def test_invalid_coords(self):
        result = analyze_inviscid([1.0, 0.0], 0.0)
        assert not result["success"]


# ---------------------------------------------------------------------------
# Faithful (viscous) solver
# ---------------------------------------------------------------------------

class TestFaithfulSolver:
    def test_naca2412_viscous(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100)
        assert result["success"]
        assert result["converged"]
        assert result["cl"] > 0.5
        assert 0.0 < result["cd"] < 0.1
        assert result["iterations"] > 0

    def test_viscous_returns_transition(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100)
        assert result["success"]
        assert 0.0 < result["x_tr_upper"] < 1.0
        assert 0.0 < result["x_tr_lower"] <= 1.0

    def test_viscous_cd_greater_than_zero(self):
        raw = generate_naca4(12)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_faithful(coords, 3.0, 1e6, 0.0, 9.0, 100)
        assert result["success"]
        assert result["cd"] > 0.0

    def test_different_reynolds(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        r_low = analyze_faithful(coords, 5.0, 1e5, 0.0, 9.0, 100)
        r_high = analyze_faithful(coords, 5.0, 3e6, 0.0, 9.0, 100)
        assert r_low["success"] and r_high["success"]
        assert r_low["cd"] > r_high["cd"]

    def test_invalid_coords(self):
        result = analyze_faithful([1.0, 0.0], 0.0, 1e6, 0.0, 9.0, 100)
        assert not result["success"]
        assert result["error"] is not None

    def test_default_parameters(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        coords = [v for x, y in paneled for v in (x, y)]
        result = analyze_faithful(coords, 0.0)
        assert result["success"]


# ---------------------------------------------------------------------------
# .dat file parsing
# ---------------------------------------------------------------------------

class TestParseDat:
    def test_parse_nonexistent_file(self):
        with pytest.raises(Exception):
            parse_dat_file("/nonexistent/path.dat")


# ---------------------------------------------------------------------------
# Batch (parallel) solvers
# ---------------------------------------------------------------------------

class TestFaithfulBatch:
    def _make_coords(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        return [v for x, y in paneled for v in (x, y)]

    def test_batch_returns_correct_count(self):
        from flexfoil._rustfoil import analyze_faithful_batch
        coords = self._make_coords()
        results = analyze_faithful_batch(coords, [0.0, 3.0, 6.0])
        assert len(results) == 3

    def test_batch_results_match_single(self):
        from flexfoil._rustfoil import analyze_faithful_batch
        coords = self._make_coords()
        alphas = [0.0, 5.0, 10.0]

        batch = analyze_faithful_batch(coords, alphas, 1e6, 0.0, 9.0, 100)
        singles = [analyze_faithful(coords, a, 1e6, 0.0, 9.0, 100) for a in alphas]

        for b, s in zip(batch, singles):
            assert b["success"] == s["success"]
            if b["success"]:
                assert abs(b["cl"] - s["cl"]) < 1e-10
                assert abs(b["cd"] - s["cd"]) < 1e-10

    def test_batch_empty_alphas(self):
        from flexfoil._rustfoil import analyze_faithful_batch
        coords = self._make_coords()
        results = analyze_faithful_batch(coords, [])
        assert results == []

    def test_batch_invalid_coords(self):
        from flexfoil._rustfoil import analyze_faithful_batch
        results = analyze_faithful_batch([1.0, 0.0], [0.0, 5.0])
        assert len(results) == 2
        assert not results[0]["success"]
        assert not results[1]["success"]


class TestDeflectFlap:
    def _make_flat(self):
        raw = generate_naca4(2412)
        return [v for x, y in raw for v in (x, y)]

    def test_returns_points(self):
        from flexfoil._rustfoil import deflect_flap
        flat = self._make_flat()
        result = deflect_flap(flat, 0.75, 10.0)
        assert len(result) > 50

    def test_zero_deflection_preserves_shape(self):
        from flexfoil._rustfoil import deflect_flap
        flat = self._make_flat()
        result = deflect_flap(flat, 0.75, 0.0)
        original = list(zip(flat[::2], flat[1::2]))
        assert len(result) == len(original)
        for (rx, ry), (ox, oy) in zip(result, original):
            assert abs(rx - ox) < 1e-8
            assert abs(ry - oy) < 1e-8

    def test_positive_deflection_moves_te_down(self):
        from flexfoil._rustfoil import deflect_flap
        flat = self._make_flat()
        original_te_y = flat[1]
        result = deflect_flap(flat, 0.75, 15.0)
        flapped_te_y = result[0][1]
        assert flapped_te_y < original_te_y

    def test_invalid_coords(self):
        from flexfoil._rustfoil import deflect_flap
        assert deflect_flap([1.0, 0.0], 0.75, 10.0) == []


class TestBLDistribution:
    def _make_coords(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        return [v for x, y in paneled for v in (x, y)]

    def test_returns_bl_data(self):
        from flexfoil._rustfoil import get_bl_distribution
        coords = self._make_coords()
        result = get_bl_distribution(coords, 5.0, 1e6, 0.0, 9.0, 100)
        assert result["success"]
        assert result["converged"]
        assert len(result["x_upper"]) > 0
        assert len(result["x_lower"]) > 0
        assert len(result["cf_upper"]) == len(result["x_upper"])
        assert len(result["delta_star_upper"]) == len(result["x_upper"])
        assert len(result["theta_upper"]) == len(result["x_upper"])
        assert len(result["h_upper"]) == len(result["x_upper"])
        assert len(result["ue_upper"]) == len(result["x_upper"])

    def test_bl_transition_locations(self):
        from flexfoil._rustfoil import get_bl_distribution
        coords = self._make_coords()
        result = get_bl_distribution(coords, 5.0, 1e6, 0.0, 9.0, 100)
        assert result["success"]
        assert 0.0 < result["x_tr_upper"] < 1.0

    def test_bl_invalid_coords(self):
        from flexfoil._rustfoil import get_bl_distribution
        result = get_bl_distribution([1.0, 0.0], 0.0, 1e6, 0.0, 9.0, 100)
        assert not result["success"]


class TestPythonAPI:
    """Tests for the Python-level Airfoil API enhancements."""

    def test_solve_returns_cp(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("2412")
        r = foil.solve(5.0, viscous=False)
        assert r.success
        assert r.cp is not None
        assert r.cp_x is not None
        assert len(r.cp) > 0
        assert len(r.cp_x) > 0

    def test_solve_viscous_returns_cp(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("2412")
        r = foil.solve(5.0, Re=1e6)
        assert r.success
        assert r.cp is not None
        assert len(r.cp) > 0

    def test_polar_explicit_alphas(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("0012")
        polar = foil.polar(alpha=[0.0, 5.0, 10.0])
        assert len(polar.results) == 3

    def test_polar_matrix_sweep(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("0012")
        polars = foil.polar(alpha=[0.0, 5.0], Re=[1e5, 1e6])
        assert isinstance(polars, list)
        assert len(polars) == 2
        assert polars[0].reynolds == 1e5
        assert polars[1].reynolds == 1e6

    def test_polar_matrix_re_mach(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("0012")
        polars = foil.polar(alpha=[0.0, 5.0], Re=[1e5, 1e6], mach=[0.0, 0.3])
        assert len(polars) == 4

    def test_polar_scalar_returns_single(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("0012")
        polar = foil.polar(alpha=[0.0, 5.0], Re=1e6)
        from flexfoil.polar import PolarResult
        assert isinstance(polar, PolarResult)

    def test_bl_distribution(self):
        from flexfoil import Airfoil
        foil = Airfoil.from_naca("2412")
        bl = foil.bl_distribution(5.0, Re=1e6)
        assert bl.success
        assert bl.converged
        assert len(bl.x_upper) > 0
        assert len(bl.cf_upper) == len(bl.x_upper)
        assert len(bl.delta_star_upper) == len(bl.x_upper)
        assert len(bl.theta_upper) == len(bl.x_upper)
        assert len(bl.h_upper) == len(bl.x_upper)
        assert len(bl.ue_upper) == len(bl.x_upper)
        assert 0.0 < bl.x_tr_upper < 1.0


class TestReType:
    """Tests for polar Type 2/3 (variable Re) at the Rust binding level."""

    def _make_coords(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        return [v for x, y in paneled for v in (x, y)]

    def test_type1_reynolds_eff_equals_nominal(self):
        coords = self._make_coords()
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=1)
        assert result["success"]
        assert result["reynolds_eff"] == 1e6

    def test_type2_reynolds_eff_differs(self):
        coords = self._make_coords()
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=2)
        assert result["success"]
        assert result["reynolds_eff"] != 1e6

    def test_type3_reynolds_eff_differs(self):
        coords = self._make_coords()
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=3)
        assert result["success"]
        assert result["reynolds_eff"] != 1e6

    def test_type2_formula(self):
        """Re_eff should equal Re / sqrt(|CL|) for Type 2."""
        import math
        coords = self._make_coords()
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=2)
        assert result["success"]
        expected = 1e6 / math.sqrt(abs(result["cl"]))
        assert abs(result["reynolds_eff"] - expected) < 1.0

    def test_type3_formula(self):
        """Re_eff should equal Re / |CL| for Type 3."""
        coords = self._make_coords()
        result = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=3)
        assert result["success"]
        expected = 1e6 / abs(result["cl"])
        assert abs(result["reynolds_eff"] - expected) < 1.0

    def test_type2_changes_cd(self):
        """Type 2 should produce different CD than Type 1 at the same nominal Re."""
        coords = self._make_coords()
        r1 = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=1)
        r2 = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=2)
        assert r1["success"] and r2["success"]
        assert r1["cd"] != r2["cd"]

    def test_default_re_type_is_type1(self):
        """Omitting re_type should behave as Type 1."""
        coords = self._make_coords()
        r_default = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100)
        r_type1 = analyze_faithful(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=1)
        assert r_default["reynolds_eff"] == r_type1["reynolds_eff"]
        assert abs(r_default["cl"] - r_type1["cl"]) < 1e-10

    def test_batch_type2_reynolds_eff(self):
        from flexfoil._rustfoil import analyze_faithful_batch
        coords = self._make_coords()
        results = analyze_faithful_batch(coords, [3.0, 5.0, 7.0], 1e6, 0.0, 9.0, 100, re_type=2)
        for r in results:
            assert r["success"]
            assert r["reynolds_eff"] != 1e6

    def test_batch_type2_matches_single(self):
        from flexfoil._rustfoil import analyze_faithful_batch
        coords = self._make_coords()
        alphas = [3.0, 5.0]
        batch = analyze_faithful_batch(coords, alphas, 1e6, 0.0, 9.0, 100, re_type=2)
        singles = [analyze_faithful(coords, a, 1e6, 0.0, 9.0, 100, re_type=2) for a in alphas]
        for b, s in zip(batch, singles):
            assert abs(b["cl"] - s["cl"]) < 1e-10
            assert abs(b["reynolds_eff"] - s["reynolds_eff"]) < 1e-10

    def test_bl_distribution_type2(self):
        from flexfoil._rustfoil import get_bl_distribution
        coords = self._make_coords()
        result = get_bl_distribution(coords, 5.0, 1e6, 0.0, 9.0, 100, re_type=2)
        assert result["success"]
        assert result["converged"]


class TestInviscidBatch:
    def _make_coords(self):
        raw = generate_naca4(2412)
        flat = [v for x, y in raw for v in (x, y)]
        paneled = repanel_xfoil(flat, 160)
        return [v for x, y in paneled for v in (x, y)]

    def test_batch_returns_correct_count(self):
        from flexfoil._rustfoil import analyze_inviscid_batch
        coords = self._make_coords()
        results = analyze_inviscid_batch(coords, [0.0, 3.0, 6.0])
        assert len(results) == 3

    def test_batch_results_match_single(self):
        from flexfoil._rustfoil import analyze_inviscid_batch
        coords = self._make_coords()
        alphas = [0.0, 5.0, 10.0]

        batch = analyze_inviscid_batch(coords, alphas)
        singles = [analyze_inviscid(coords, a) for a in alphas]

        for b, s in zip(batch, singles):
            assert b["success"] == s["success"]
            if b["success"]:
                assert abs(b["cl"] - s["cl"]) < 1e-10

    def test_batch_invalid_coords(self):
        from flexfoil._rustfoil import analyze_inviscid_batch
        results = analyze_inviscid_batch([1.0, 0.0], [0.0, 5.0])
        assert len(results) == 2
        assert not results[0]["success"]
