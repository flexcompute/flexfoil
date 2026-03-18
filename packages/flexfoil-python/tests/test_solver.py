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
