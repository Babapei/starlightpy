import numpy as np

from starlightpy.extinction import get_extinction_curve


def test_ccm_is_unity_near_v_band():
    q = get_extinction_curve([5500.0], law="CCM", r_v=3.1)
    assert abs(q[0] - 1.0) < 0.02


def test_more_extinction_in_the_blue():
    q = get_extinction_curve([4000.0, 5500.0, 8000.0], law="CCM")
    assert q[0] > q[1] > q[2]


def test_calzetti_and_gordon_are_finite():
    wave = np.linspace(3000, 8000, 50)
    assert np.all(np.isfinite(get_extinction_curve(wave, law="CAL")))
    assert np.all(np.isfinite(get_extinction_curve(wave, law="GD1")))
