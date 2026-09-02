import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _noiseless_mock():
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


def test_g4_estimate_error_recovers():
    wave, bases, true_x, true_av, config, flux, _err = _noiseless_mock()
    config.estimate_error = True
    with pytest.warns(UserWarning, match="not a true"):
        result = fit_spectrum(wave, flux, None, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_g4_missing_error_without_estimate_raises():
    wave, bases, _true_x, _true_av, config, flux, _err = _noiseless_mock()
    config.estimate_error = False
    with pytest.raises(ValueError, match="error is required"):
        fit_spectrum(wave, flux, None, bases, config=config)


def test_g4_real_error_not_overwritten():
    wave, bases, true_x, true_av, config, flux, err = _noiseless_mock()
    with_flag = FitConfig(
        a_v_bounds=config.a_v_bounds,
        a_v_step=config.a_v_step,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
        estimate_error=True,
    )
    without_flag = FitConfig(
        a_v_bounds=config.a_v_bounds,
        a_v_step=config.a_v_step,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
        estimate_error=False,
    )
    with pytest.warns(UserWarning):
        r_flag = fit_spectrum(wave, flux, err, bases, config=with_flag)
        r_plain = fit_spectrum(wave, flux, err, bases, config=without_flag)
    assert r_flag.chi2 == pytest.approx(r_plain.chi2)
    np.testing.assert_allclose(r_flag.x_fraction, r_plain.x_fraction)
    np.testing.assert_allclose(r_flag.x_fraction, true_x, atol=0.12)
    assert abs(r_flag.a_v - true_av) <= 0.10
