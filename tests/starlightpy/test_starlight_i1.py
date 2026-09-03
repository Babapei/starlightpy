import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _i1_setup():
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
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None)
    return wave, bases, true_x, true_av, config, flux, err


def test_i1_wavelength_and_error_match_fit_grid():
    wave, bases, true_x, true_av, config, flux, err = _i1_setup()
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.wavelength is not None
    assert result.error is not None
    np.testing.assert_allclose(result.wavelength, wave)
    assert result.error.shape == wave.shape
    assert result.model.shape == wave.shape
    assert result.good is not None and result.good.shape == wave.shape
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_i1_error_is_observed_flux_units():
    wave, bases, _x, _av, config, flux, err = _i1_setup()
    flux = flux * 10.0
    err = err * 10.0
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.error is not None
    np.testing.assert_allclose(result.error, err)
    assert result.obs_scale > 2.0
    assert not np.allclose(result.error, err / result.obs_scale)


def test_i1_estimated_error_is_positive_and_same_length():
    wave, bases, true_x, true_av, config, flux, _err = _i1_setup()
    config.estimate_error = True
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, None, bases, config=config)
    assert result.error is not None
    assert result.error.shape == wave.shape
    assert np.all(np.isfinite(result.error))
    assert np.all(result.error > 0)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_i1_g1_fields_still_present():
    wave, bases, _x, _av, config, flux, err = _i1_setup()
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.good is not None and result.good.shape == wave.shape
    assert result.obs_scale > 0
    assert result.config is not None
    assert result.dropped is not None and result.dropped.shape == (bases.shape[1],)
    assert not np.any(result.dropped)
