import math

import pytest

from flexfoil.airfoil import BLResult


def test_bl_summary_locates_negative_shear_without_inventing_bubble_bounds():
    bl = BLResult(True, True, x_tr_upper=.3, x_upper=[.1, .4, .8],
                  h_upper=[2.6, 4., 2.], cf_upper=[.001, -.002, 0.])
    data = bl.summary()
    assert data["max_h_upper"] == 4.
    assert data["min_cf_upper"] == -.002
    assert data["negative_cf_x_upper"] == [.4]
    assert data["x_tr_upper"] == .3
    assert data["max_h_lower"] is None
    assert data["negative_cf_x_lower"] is None


@pytest.mark.parametrize("success,converged", [(False, False), (True, False), (False, True)])
def test_invalid_status_does_not_produce_summary_numbers(success, converged):
    bl = BLResult(success, converged, x_upper=[.1], h_upper=[2.6], cf_upper=[.001])
    assert bl.summary()["max_h_upper"] is None
    assert bl.summary()["negative_cf_x_upper"] is None


@pytest.mark.parametrize("h,cf", [([math.nan], [.001]), ([2.6], [math.inf]), ([2.6, 3.], [.001])])
def test_invalid_or_mismatched_station_arrays_remain_unknown(h, cf):
    data = BLResult(True, True, x_upper=[.1], h_upper=h, cf_upper=cf).summary()
    assert data["max_h_upper"] is None
    assert data["min_cf_upper"] is None
