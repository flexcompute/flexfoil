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


def test_empty_complete_export():
    polar = PolarResult("test", 1e6, 0, 9)
    assert polar.to_dict(include_failed=True)["alpha"] == []
    assert polar.to_dataframe(include_failed=True).empty


def test_generated_polar_records_conditions(monkeypatch):
    foil = Airfoil("test", [(1., 0.), (0., .1), (1., -.1)], [])
    monkeypatch.setattr(foil, "_polar_batch", lambda *args, **kwargs: [point(0)])
    polar = foil.polar(alpha=[0], Re=1e6, re_type=2, max_iter=42,
                       xstrip_upper=.2, xstrip_lower=.6, store=False)
    data = polar.to_dict(include_failed=True)
    assert data["max_iter"] == [42]
    assert data["re_type"] == [2]
    assert data["xstrip_lower"] == [.6]
    assert data["geometry_hash"] == [foil.hash]
    assert data["solver_version"][0]
