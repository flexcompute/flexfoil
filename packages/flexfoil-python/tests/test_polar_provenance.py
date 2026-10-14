"""Source/vendored polar calls retain unknown provenance without changing solver results."""
from importlib.metadata import PackageNotFoundError
from types import SimpleNamespace
import pytest
import flexfoil

@pytest.mark.parametrize('parallel', [False, True])
def test_polar_without_distribution_metadata(monkeypatch, parallel):
    def missing(_name):
        raise PackageNotFoundError('flexfoil')
    monkeypatch.setattr('importlib.metadata.distribution', missing)
    result = flexfoil.naca('0012').polar(alpha=[0, 4], viscous=False, store=False, parallel=parallel)
    assert result.solver_version is None
    assert len(result.results) == 2
    assert all(point.success for point in result.results)
    assert result.results[1].cl > result.results[0].cl


def test_shadowed_wheel_metadata_is_not_claimed(monkeypatch, tmp_path):
    monkeypatch.setattr('importlib.metadata.distribution', lambda _: SimpleNamespace(version='0.0.1-stale', locate_file=lambda _: tmp_path/'other-package/airfoil.py'))
    result = flexfoil.naca('0012').polar(alpha=[0], viscous=False, store=False)
    assert result.solver_version is None
    assert result.results[0].success


def test_matching_installed_distribution_is_retained(monkeypatch):
    import flexfoil.airfoil
    monkeypatch.setattr('importlib.metadata.distribution', lambda _: SimpleNamespace(version='1.4.2', locate_file=lambda _: flexfoil.airfoil.__file__))
    result = flexfoil.naca('0012').polar(alpha=[0], viscous=False, store=False)
    assert result.solver_version == '1.4.2'
