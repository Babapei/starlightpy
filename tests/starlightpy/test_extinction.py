from starlightpy.extinction import get_extinction_curve


def test_ccm_near_v_band():
    q = get_extinction_curve([5500.0], law="CCM", r_v=3.1)
    assert abs(q[0] - 1.0) < 0.02


def test_bluer_is_more_extincted():
    q = get_extinction_curve([4000.0, 5500.0, 8000.0], law="CCM")
    assert q[0] > q[1] > q[2]
