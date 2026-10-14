"""Complete polar exports preserve coverage and recorded assumptions."""
from flexfoil.airfoil import Airfoil, SolveResult
from flexfoil.polar import PolarResult


def point(alpha, success=True, converged=True):
    return SolveResult(1, .01, -.1, converged, 3, .001, .4, .7, alpha,
                       1e6, 0, 9, success, None if success else "singular system")


def test_complete_export_keeps_failed_and_unconverged_points():
    polar = PolarResult("test", 1e6, 0, 9,
                        [point(3), point(1, False, False), point(2, True, False)],
                        geometry_hash="abc", solver_version="test", max_iter=10,
                        viscous=True, re_type=2, xstrip_upper=.2, xstrip_lower=.6)
    assert polar.to_dict()["alpha"] == [3]
    data = polar.to_dict(include_failed=True, summary=True)
    assert data["alpha"] == [3, 1, 2]
    assert data["cl"] == [1, None, 1]
    assert data["success"] == [True, False, True]
    assert data["converged"] == [True, False, False]
    assert data["error"][1] == "singular system"
    assert data["xstrip_upper"] == [.2] * 3
    assert data["geometry_hash"] == ["abc"] * 3
    assert data["_summary"]["cl_max"] == 1
    frame = polar.to_dataframe(include_failed=True, summary=True)
    assert frame["alpha"].tolist() == [3, 1, 2]
    assert frame["cl"].isna().tolist() == [False, True, False]
    assert frame.attrs["cl_max"] == 1


def test_inviscid_complete_export_records_requested_conditions(monkeypatch):
    def raw_for(alpha):
        if alpha == 4:
            return {"success": False, "error": "singular system"}
        return {"success": True, "cl": 0.5, "cm": -0.05}

    monkeypatch.setattr("flexfoil.airfoil.analyze_inviscid", lambda _coords, alpha: raw_for(alpha))
    monkeypatch.setattr(
        "flexfoil._rustfoil.analyze_inviscid_batch",
        lambda _coords, alphas: [raw_for(alpha) for alpha in alphas],
    )
    foil = Airfoil("test", [], [])
    for parallel in (False, True):
        polar = foil.polar(
            alpha=[2, 4], Re=2e6, mach=0.2, ncrit=7,
            viscous=False, store=False, parallel=parallel,
        )
        assert [(row.reynolds, row.mach, row.ncrit) for row in polar.results] == [(2e6, 0.2, 7)] * 2
        data = polar.to_dict(include_failed=True)
        assert data["alpha"] == [2, 4]
        assert data["cl"] == [0.5, None]
        assert data["reynolds"] == [2e6, 2e6]
        assert data["mach"] == [0.2, 0.2]
        assert data["ncrit"] == [7, 7]


def test_empty_complete_export():
    polar = PolarResult("test", 1e6, 0, 9)
    assert polar.to_dict(include_failed=True)["alpha"] == []
    assert polar.to_dataframe(include_failed=True).empty


def test_generated_polar_records_conditions(monkeypatch):
    from types import SimpleNamespace
    import flexfoil.airfoil
    monkeypatch.setattr('importlib.metadata.distribution', lambda _: SimpleNamespace(version='1.4.2', locate_file=lambda _: flexfoil.airfoil.__file__))
    foil = Airfoil("test", [(1., 0.), (0., .1), (1., -.1)], [])
    monkeypatch.setattr(foil, "_polar_batch", lambda *args, **kwargs: [point(0)])
    polar = foil.polar(alpha=[0], Re=1e6, re_type=2, max_iter=42,
                       xstrip_upper=.2, xstrip_lower=.6, store=False)
    data = polar.to_dict(include_failed=True)
    assert data["max_iter"] == [42]
    assert data["re_type"] == [2]
    assert data["xstrip_lower"] == [.6]
    assert data["geometry_hash"] == [foil.hash]
    assert data["solver_version"][0] == "1.4.2"
