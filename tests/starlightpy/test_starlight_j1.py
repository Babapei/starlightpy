import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _j1_setup():
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


def test_j1_flux_matches_model_grid():
    wave, bases, true_x, true_av, config, flux, err = _j1_setup()
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.flux is not None
    assert result.model is not None
    assert result.wavelength is not None
    assert result.flux.shape == result.model.shape
    assert result.flux.shape == result.wavelength.shape
    np.testing.assert_allclose(result.flux, flux)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_j1_flux_is_observed_units_not_normalized():
    wave, bases, _x, _av, config, flux, err = _j1_setup()
    flux = flux * 10.0
    err = err * 10.0
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.flux is not None
    np.testing.assert_allclose(result.flux, flux)
    assert result.obs_scale > 2.0
    assert not np.allclose(result.flux, flux / result.obs_scale)


def test_j1_wavelength_and_error_still_present():
    wave, bases, _x, _av, config, flux, err = _j1_setup()
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.wavelength is not None
    assert result.error is not None
    np.testing.assert_allclose(result.wavelength, wave)
    np.testing.assert_allclose(result.error, err)
