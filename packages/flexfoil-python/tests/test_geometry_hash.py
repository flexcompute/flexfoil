"""Tests for placement-aware geometry hashing (the run-cache key).

The hash keys ``runs.idx_cache_key``, so a collision between two different
geometries is a silently wrong answer. These tests pin (a) backward
compatibility with the coordinate-only digest that existing local databases
are keyed on, and (b) the exact element/placement encoding, which must stay
byte-for-byte identical to the browser implementation in
``flexfoil-ui/src/lib/airfoilHash.ts``.
"""

from __future__ import annotations

import importlib.util
import sys
import types
from pathlib import Path

_AIRFOIL_PY = Path(__file__).resolve().parents[1] / "src" / "flexfoil" / "airfoil.py"

_RUSTFOIL_NAMES = (
    "analyze_faithful",
    "analyze_inviscid",
    "deflect_flap",
    "generate_naca4",
    "get_bl_distribution",
    "parse_dat_file",
    "repanel_xfoil",
)


def _load_airfoil_module() -> types.ModuleType:
    """Load ``src/flexfoil/airfoil.py`` without the compiled extension.

    The hashing code is pure Python, but ``airfoil.py`` imports the
    maturin-built ``flexfoil._rustfoil`` at module scope, which is not present
    in a fresh checkout. Stub that import for the duration of the load and then
    restore ``sys.modules`` so that no other test in the session sees the stub.
    """
    sentinel = object()
    previous = sys.modules.get("flexfoil._rustfoil", sentinel)
    if previous is sentinel:
        stub = types.ModuleType("flexfoil._rustfoil")
        for name in _RUSTFOIL_NAMES:
            setattr(stub, name, _unavailable)
        sys.modules["flexfoil._rustfoil"] = stub

    spec = importlib.util.spec_from_file_location(
        "flexfoil_airfoil_source", _AIRFOIL_PY
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    try:
        spec.loader.exec_module(module)
    finally:
        if previous is sentinel:
            sys.modules.pop("flexfoil._rustfoil", None)
        else:
            sys.modules["flexfoil._rustfoil"] = previous
    return module


def _unavailable(*_args, **_kwargs):  # pragma: no cover - never called here
    raise RuntimeError("flexfoil._rustfoil is not built in this checkout")


_mod = _load_airfoil_module()

Airfoil = _mod.Airfoil
GeometryElement = _mod.GeometryElement
Placement = _mod.Placement
canonical_geometry = _mod.canonical_geometry
geometry_hash = _mod.geometry_hash


# ---------------------------------------------------------------------------
# Shared cross-language fixture.
#
# The identical fixture and the identical expected canonical strings appear in
# flexfoil-ui/src/lib/airfoilHash.test.ts. If either side's canonicalisation
# drifts, one of the two pinned strings stops matching.
# ---------------------------------------------------------------------------

MAIN = [(0.0, 0.0), (1.0, 0.0), (0.5, 0.06), (0.25, -0.03)]
FLAP = [(0.8, 0.0), (1.2, 0.0), (1.0, 0.02)]

MAIN_COORDS = (
    "0.00000000,0.00000000|1.00000000,0.00000000|"
    "0.50000000,0.06000000|0.25000000,-0.03000000"
)
FLAP_COORDS = "0.80000000,0.00000000|1.20000000,0.00000000|1.00000000,0.02000000"

# Digest of MAIN alone, as produced by the coordinate-only implementation that
# shipped in flexfoil 1.1.6. Every run already cached in a user's
# ~/.flexfoil/runs.db is keyed on this; it must never change.
MAIN_LEGACY_DIGEST = "cadb528511d42bd7"

PLACEMENT = Placement(
    rotation=25.0,
    pivot=(0.8, 0.01),
    translation=(0.02, -0.05),
    scale=1.0,
)

# The placement/element encoding is the part that MUST agree byte-for-byte with
# the TypeScript implementation. This exact literal is asserted on both sides.
PLACEMENT_SUFFIX = "@25.00000000,0.80000000,0.01000000,0.02000000,-0.05000000,1.00000000"
ELEMENT_SEP = "#"


def _hash(coords, placement=None) -> str:
    return geometry_hash([GeometryElement(coords, placement)])


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    def test_legacy_single_element_digest_is_pinned(self):
        assert _hash(MAIN) == MAIN_LEGACY_DIGEST

    def test_single_unplaced_element_canonicalises_to_coordinates_alone(self):
        assert canonical_geometry([GeometryElement(MAIN)]) == MAIN_COORDS

    def test_airfoil_hash_matches_legacy_digest(self):
        assert Airfoil("t", [], MAIN).hash == MAIN_LEGACY_DIGEST

    def test_identity_placement_hashes_as_absent(self):
        assert _hash(MAIN, Placement()) == MAIN_LEGACY_DIGEST
        assert Placement().is_identity

    def test_pivot_only_placement_is_identity(self):
        # A pivot with no rotation/translation/scale moves nothing.
        pivot_only = Placement(pivot=(0.75, 0.02))
        assert pivot_only.is_identity
        assert _hash(MAIN, pivot_only) == MAIN_LEGACY_DIGEST

    def test_airfoil_placement_defaults_to_none(self):
        assert Airfoil("t", [], MAIN).placement is None

    def test_negative_zero_collapses(self):
        # Python's %.8f would emit "-0.00000000" where JS toFixed emits
        # "0.00000000"; the formatter collapses -0.0 so the two agree.
        assert _hash(MAIN, Placement(translation=(-0.0, 0.1))) == _hash(
            MAIN, Placement(translation=(0.0, 0.1))
        )

    def test_digest_length_unchanged(self):
        assert len(_hash(MAIN)) == 16
        assert len(_hash(MAIN, PLACEMENT)) == 16


# ---------------------------------------------------------------------------
# Placement sensitivity — the bug this hash exists to prevent
# ---------------------------------------------------------------------------

class TestPlacementSensitivity:
    def test_translation_changes_the_hash(self):
        a = _hash(MAIN, Placement(translation=(0.01, 0.0)))
        b = _hash(MAIN, Placement(translation=(0.02, 0.0)))
        assert a != b
        assert a != MAIN_LEGACY_DIGEST

    def test_rotation_changes_the_hash(self):
        a = _hash(MAIN, Placement(rotation=10.0, pivot=(0.75, 0.0)))
        b = _hash(MAIN, Placement(rotation=20.0, pivot=(0.75, 0.0)))
        assert a != b

    def test_pivot_changes_the_hash_when_rotated(self):
        a = _hash(MAIN, Placement(rotation=10.0, pivot=(0.70, 0.0)))
        b = _hash(MAIN, Placement(rotation=10.0, pivot=(0.75, 0.0)))
        assert a != b

    def test_scale_changes_the_hash(self):
        assert _hash(MAIN, Placement(scale=1.0000001)) != _hash(
            MAIN, Placement(scale=0.3)
        )

    def test_sub_precision_differences_share_a_key(self):
        # Documented precision limit: finer than 1e-8 is the same cache key.
        assert _hash(MAIN, Placement(translation=(0.1, 0.0))) == _hash(
            MAIN, Placement(translation=(0.1 + 1e-12, 0.0))
        )

    def test_placed_airfoil_differs_from_unplaced(self):
        placed = Airfoil("t", [], MAIN, placement=Placement(translation=(0.0, 0.05)))
        assert placed.hash != MAIN_LEGACY_DIGEST


# ---------------------------------------------------------------------------
# Element count / order, and cross-language encoding
# ---------------------------------------------------------------------------

class TestAssembly:
    def test_placement_suffix_matches_the_typescript_encoding(self):
        canonical = canonical_geometry(
            [GeometryElement(MAIN), GeometryElement(FLAP, PLACEMENT)]
        )
        assert canonical == MAIN_COORDS + ELEMENT_SEP + FLAP_COORDS + PLACEMENT_SUFFIX

    def test_element_count_matters(self):
        one = geometry_hash([GeometryElement(MAIN)])
        two = geometry_hash([GeometryElement(MAIN), GeometryElement(MAIN)])
        assert one != two

    def test_element_order_matters(self):
        ab = geometry_hash([GeometryElement(MAIN), GeometryElement(FLAP)])
        ba = geometry_hash([GeometryElement(FLAP), GeometryElement(MAIN)])
        assert ab != ba

    def test_same_shapes_different_gap_do_not_collide(self):
        closed = geometry_hash([
            GeometryElement(MAIN),
            GeometryElement(FLAP, Placement(translation=(0.0, -0.01))),
        ])
        opened = geometry_hash([
            GeometryElement(MAIN),
            GeometryElement(FLAP, Placement(translation=(0.0, -0.03))),
        ])
        assert closed != opened

    def test_repeated_calls_are_stable(self):
        elements = [GeometryElement(MAIN), GeometryElement(FLAP, PLACEMENT)]
        assert geometry_hash(elements) == geometry_hash(elements)

    def test_empty_assembly_does_not_raise(self):
        assert len(geometry_hash([])) == 16
