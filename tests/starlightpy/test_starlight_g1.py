import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _g1_setup():
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


def test_g1_result_fields_on_noiseless_recovery():
    wave, bases, true_x, true_av, config, flux, err = _g1_setup()
    with pytest.warns(UserWarning, match="velocity dispersion is not fitted"):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.good is not None and result.good.shape == wave.shape
    assert bool(np.all(result.good))
    assert result.obs_scale > 0
    assert result.config is not None
    assert result.config.norm_window == config.norm_window
    assert result.config is not config
    assert result.dropped is not None and result.dropped.shape == (bases.shape[1],)
    assert not np.any(result.dropped)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_g1_nan_flux_raises():
    wave, bases, _x, _av, config, flux, err = _g1_setup()
    flux = flux.copy()
    flux[10] = np.nan
    with pytest.raises(ValueError, match="non-finite"):
        fit_spectrum(wave, flux, err, bases, config=config)


def test_g1_norm_window_outside_raises():
    wave, bases, _x, _av, config, flux, err = _g1_setup()
    config.norm_window = (100.0, 200.0)
    with pytest.raises(ValueError, match="Normalization window"):
        fit_spectrum(wave, flux, err, bases, config=config)


def test_g1_ex0_sets_dropped():
    wave, bases, _x, _av, config, flux, err = _g1_setup()
    true_x = np.array([0.55, 0.40, 0.05])
    config.x_min_keep = 0.10
    flux, err = mock_observation(wave, bases, true_x, 0.30, 0.0, 0.0, config=config, snr=None)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.dropped is not None
    assert result.dropped[2]
    assert not result.dropped[0]
    assert result.x[2] == 0.0
