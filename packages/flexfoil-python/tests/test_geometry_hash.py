"""Tests for placement-aware geometry hashing (the run-cache key).

The hash keys ``runs.idx_cache_key``, so a collision between two different
geometries returns the wrong answer for one of them. These tests pin

* backward compatibility with the coordinate-only digest that existing local
  databases are keyed on (the frozen legacy digest, which is deliberately not
  comparable with the browser's legacy digest — see the module notes in
  ``src/flexfoil/airfoil.py``), and
* the shared assembly digest, whose exact hex literals also appear in
  ``flexfoil-ui/src/lib/airfoilHash.test.ts``. Those literals are the
  cross-language contract: if either implementation's canonicalisation drifts,
  one of the two suites fails.
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
is_legacy_single_element = _mod.is_legacy_single_element
shared_canonical_geometry = _mod.shared_canonical_geometry


# ---------------------------------------------------------------------------
# Shared cross-language fixtures.
#
# The identical fixtures, canonical strings and digest literals appear in
# flexfoil-ui/src/lib/airfoilHash.test.ts.
# ---------------------------------------------------------------------------

MAIN = [(0.0, 0.0), (1.0, 0.0), (0.5, 0.06), (0.25, -0.03)]
FLAP = [(0.8, 0.0), (1.2, 0.0), (1.0, 0.02)]

# Legacy encoding joins nodes with "|"; the shared encoding joins them with ";"
# (which happens to be what the browser's frozen legacy encoding also uses).
MAIN_LEGACY_COORDS = (
    "0.00000000,0.00000000|1.00000000,0.00000000|"
    "0.50000000,0.06000000|0.25000000,-0.03000000"
)
MAIN_SHARED_COORDS = (
    "0.00000000,0.00000000;1.00000000,0.00000000;"
    "0.50000000,0.06000000;0.25000000,-0.03000000"
)
FLAP_SHARED_COORDS = "0.80000000,0.00000000;1.20000000,0.00000000;1.00000000,0.02000000"

# Digest of MAIN alone, as produced by the coordinate-only implementation that
# shipped in flexfoil 1.1.6. Every run already cached in a user's
# ~/.flexfoil/runs.db is keyed on this; it must never change. It is not
# comparable with the browser's legacy digest for the same coordinates.
MAIN_LEGACY_DIGEST = "cadb528511d42bd7"

PLACEMENT = Placement(
    rotation=25.0,
    pivot=(0.8, 0.01),
    translation=(0.02, -0.05),
    scale=1.0,
)

SHARED_SCHEMA = "ffgeom1:"
PLACEMENT_SUFFIX = "@25.00000000,0.80000000,0.01000000,0.02000000,-0.05000000,1.00000000"
ELEMENT_SEP = "#"

# Vector 1 — a two-element assembly with a non-identity placement on the second
# element. Pinned byte-for-byte against TypeScript.
SHARED_ASSEMBLY_CANONICAL = (
    SHARED_SCHEMA + MAIN_SHARED_COORDS + ELEMENT_SEP + FLAP_SHARED_COORDS + PLACEMENT_SUFFIX
)
SHARED_ASSEMBLY_DIGEST = (
    "cbec85f2a62ddd0613c392aec6b3ecc751665d6096a2206e67bed2e2b1c5e53f"
)

# Vector 2 — float-formatting edge cases, which is where the two languages
# diverge if the formatter is not shared: negative zero in a coordinate and in a
# placement field, and 0.001953125 (an odd multiple of 1/512, the smallest
# family of doubles that lands exactly on a half at the eighth decimal, where
# Python's round-half-to-even and JavaScript's round-half-away-from-zero
# disagree). Both must render 0.00195313.
EDGE = [(-0.0, 0.0), (0.001953125, -0.001953125), (0.5, -0.0)]
EDGE_PLACEMENT = Placement(
    rotation=-0.0,
    pivot=(0.0, 0.0),
    translation=(0.001953125, -0.0),
    scale=1.001953125,
)
SHARED_EDGE_CANONICAL = (
    SHARED_SCHEMA
    + "0.00000000,0.00000000;0.00195313,-0.00195313;0.50000000,0.00000000"
    + "@0.00000000,0.00000000,0.00000000,0.00195313,0.00000000,1.00195313"
)
SHARED_EDGE_DIGEST = (
    "3b3c8bf81354ec30ee2fe3f151771c7ed2791eff9a70d76cb2080b5807ed9c98"
)


def _hash(coords, placement=None) -> str:
    return geometry_hash([GeometryElement(coords, placement)])


# ---------------------------------------------------------------------------
# Backward compatibility
# ---------------------------------------------------------------------------

class TestBackwardCompatibility:
    def test_legacy_single_element_digest_is_pinned(self):
        assert _hash(MAIN) == MAIN_LEGACY_DIGEST

    def test_single_unplaced_element_canonicalises_to_coordinates_alone(self):
        assert canonical_geometry([GeometryElement(MAIN, None)]) == MAIN_LEGACY_COORDS

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

    def test_legacy_digest_length_unchanged(self):
        assert len(_hash(MAIN)) == 16

    def test_only_the_lone_unplaced_element_routes_to_legacy(self):
        assert is_legacy_single_element([GeometryElement(MAIN, None)])
        assert is_legacy_single_element([GeometryElement(MAIN, Placement())])
        assert not is_legacy_single_element([GeometryElement(MAIN, PLACEMENT)])
        assert not is_legacy_single_element(
            [GeometryElement(MAIN, None), GeometryElement(FLAP, None)]
        )

    def test_shared_form_cannot_collide_with_the_legacy_form(self):
        # The schema tag keeps the two canonical spaces disjoint, so one
        # geometry always maps to exactly one digest.
        assert shared_canonical_geometry([GeometryElement(MAIN, None)]).startswith(
            SHARED_SCHEMA
        )
        assert not canonical_geometry([GeometryElement(MAIN, None)]).startswith(
            SHARED_SCHEMA
        )


# ---------------------------------------------------------------------------
# Placement sensitivity — the collision this hash exists to prevent
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
# Element count and order
# ---------------------------------------------------------------------------

class TestAssembly:
    def test_element_count_matters(self):
        one = geometry_hash([GeometryElement(MAIN, None)])
        two = geometry_hash([GeometryElement(MAIN, None), GeometryElement(MAIN, None)])
        assert one != two

    def test_element_order_matters(self):
        ab = geometry_hash([GeometryElement(MAIN, None), GeometryElement(FLAP, None)])
        ba = geometry_hash([GeometryElement(FLAP, None), GeometryElement(MAIN, None)])
        assert ab != ba

    def test_same_shapes_different_gap_do_not_collide(self):
        closed = geometry_hash([
            GeometryElement(MAIN, None),
            GeometryElement(FLAP, Placement(translation=(0.0, -0.01))),
        ])
        opened = geometry_hash([
            GeometryElement(MAIN, None),
            GeometryElement(FLAP, Placement(translation=(0.0, -0.03))),
        ])
        assert closed != opened

    def test_repeated_calls_are_stable(self):
        elements = [GeometryElement(MAIN, None), GeometryElement(FLAP, PLACEMENT)]
        assert geometry_hash(elements) == geometry_hash(elements)

    def test_empty_assembly_does_not_raise(self):
        assert len(geometry_hash([])) == 64


# ---------------------------------------------------------------------------
# Shared cross-language vectors. These digests are pinned as the same literal
# hex strings in flexfoil-ui/src/lib/airfoilHash.test.ts.
# ---------------------------------------------------------------------------

class TestSharedVectors:
    ASSEMBLY = [GeometryElement(MAIN, None), GeometryElement(FLAP, PLACEMENT)]
    EDGE_ELEMENTS = [GeometryElement(EDGE, EDGE_PLACEMENT)]

    def test_assembly_vector_canonicalises_as_typescript_does(self):
        assert shared_canonical_geometry(self.ASSEMBLY) == SHARED_ASSEMBLY_CANONICAL
        assert canonical_geometry(self.ASSEMBLY) == SHARED_ASSEMBLY_CANONICAL

    def test_assembly_vector_digest_is_shared_with_typescript(self):
        assert geometry_hash(self.ASSEMBLY) == SHARED_ASSEMBLY_DIGEST

    def test_float_formatting_vector_canonicalises_as_typescript_does(self):
        assert canonical_geometry(self.EDGE_ELEMENTS) == SHARED_EDGE_CANONICAL

    def test_float_formatting_vector_digest_is_shared_with_typescript(self):
        assert geometry_hash(self.EDGE_ELEMENTS) == SHARED_EDGE_DIGEST

    def test_signed_zero_normalises_in_coordinates_and_placement_fields(self):
        neg_coord = geometry_hash([GeometryElement([(-0.0, -0.0)], PLACEMENT)])
        pos_coord = geometry_hash([GeometryElement([(0.0, 0.0)], PLACEMENT)])
        assert neg_coord == pos_coord

        neg_field = _hash(MAIN, Placement(translation=(-0.0, 0.1)))
        pos_field = _hash(MAIN, Placement(translation=(0.0, 0.1)))
        assert neg_field == pos_field

    def test_non_zero_negative_rounding_to_zero_keeps_its_sign(self):
        # Both languages emit "-0.00000000" here; only exact -0.0 is normalised.
        canonical = shared_canonical_geometry(
            [GeometryElement([(-1e-12, 0.0)], PLACEMENT)]
        )
        assert "-0.00000000,0.00000000" in canonical

    def test_shared_digest_is_the_full_sha256(self):
        assert len(geometry_hash(self.ASSEMBLY)) == 64
