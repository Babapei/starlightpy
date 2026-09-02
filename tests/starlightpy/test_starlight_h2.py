import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _h2_noiseless():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    config = FitConfig(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
    )
    flux, err = mock_observation(
        wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None
    )
    return wave, bases, true_x, true_av, config, flux, err


def test_h2_chi2_slice_covers_true_av():
    wave, bases, true_x, true_av, config, flux, err = _h2_noiseless()
    config.error_method = "chi2_slice"
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.errors is not None
    assert "a_v" in result.errors
    hw = float(result.errors["a_v"])
    assert np.isfinite(hw) and hw >= 0.0
    assert abs(result.a_v - true_av) <= hw + 1e-9 or abs(result.a_v - true_av) <= 0.10
    assert true_av <= result.a_v + hw + 1e-9
    assert true_av >= result.a_v - hw - 1e-9


def test_h2_default_errors_are_none():
    wave, bases, _true_x, _true_av, config, flux, err = _h2_noiseless()
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.errors is None


def test_h2_repeat_x_scatter():
    wave, bases, true_x, true_av, config, flux, err = _h2_noiseless()
    flux, err = mock_observation(
        wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=30.0, rng=np.random.default_rng(1)
    )
    config.error_method = "repeat"
    config.n_repeat = 5
    config.repeat_seed = 2
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.errors is not None
    xs = np.asarray(result.errors["x"])
    assert xs.shape == (3,)
    assert np.all(xs >= 0)
    assert np.any(xs > 0)
