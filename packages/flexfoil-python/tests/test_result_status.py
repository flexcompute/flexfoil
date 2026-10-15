"""Failed solves must never count as usable polar points."""

import pytest

from flexfoil.airfoil import Airfoil, SolveResult
from flexfoil.polar import PolarResult


@pytest.mark.parametrize("parallel", [False, True])
def test_failed_inviscid_sweep_retains_failure(monkeypatch, parallel):
    failure = {"success": False, "error": "singular system"}
    monkeypatch.setattr("flexfoil.airfoil.analyze_inviscid", lambda *args: failure)
    monkeypatch.setattr("flexfoil._rustfoil.analyze_inviscid_batch", lambda *args: [failure])
    foil = Airfoil("test", [], [])
    polar = foil.polar(alpha=[0], viscous=False, store=False, parallel=parallel)
    assert len(polar.results) == 1
    assert not polar.results[0].success
    assert not polar.results[0].converged
    assert polar.results[0].error == "singular system"
    assert polar.converged == []
    assert polar.cl_max is None
    assert polar.to_dict()["alpha"] == []


def test_polar_rejects_inconsistent_success_flag():
    failed = SolveResult(99, 0, 0, True, 0, 0, 1, 1, 0, 0, 0, 0, False)
    polar = PolarResult("test", 0, 0, 0, [failed])
    assert polar.converged == []
    assert polar.cl_max is None
