import numpy as np
import pytest

from starlightpy import FitConfig, air_to_vacuum, fit_spectrum, resample_to, to_rest_frame, vacuum_to_air
from starlightpy.simulate import default_absorption_bases, mock_observation


def _rest_mock():
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
        redshift=0.0,
    )
    flux, err = mock_observation(
        wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None
    )
    return wave, bases, true_x, true_av, config, flux, err


def test_air_vacuum_roundtrip():
    wave = np.linspace(3700.0, 7000.0, 200)
    air = vacuum_to_air(wave)
    back = air_to_vacuum(air)
    np.testing.assert_allclose(back, wave, rtol=1e-6, atol=0.0)


def test_g2_to_rest_frame_then_fit_recovers():
    z = 0.02
    wave, bases, true_x, true_av, config, flux_rest, err_rest = _rest_mock()
    wave_obs = wave * (1.0 + z)
    flux_obs = flux_rest / (1.0 + z)
    err_obs = err_rest / (1.0 + z)
    wr, fr, er = to_rest_frame(wave_obs, flux_obs, z, err_obs)
    flux_fit = resample_to(wr, fr, wave)
    err_fit = resample_to(wr, er, wave)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux_fit, err_fit, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_g2_skipping_rest_frame_is_worse():
    z = 0.02
    wave, bases, true_x, true_av, config, flux_rest, err_rest = _rest_mock()
    wave_obs = wave * (1.0 + z)
    flux_obs = flux_rest / (1.0 + z)
    err_obs = err_rest / (1.0 + z)
    wr, fr, er = to_rest_frame(wave_obs, flux_obs, z, err_obs)
    flux_ok = resample_to(wr, fr, wave)
    err_ok = resample_to(wr, er, wave)
    flux_wrong = resample_to(wave_obs, flux_obs, wave)
    err_wrong = resample_to(wave_obs, err_obs, wave)
    with pytest.warns(UserWarning):
        good = fit_spectrum(wave, flux_ok, err_ok, bases, config=config)
        bad = fit_spectrum(wave, flux_wrong, err_wrong, bases, config=config)
    assert bad.chi2 > good.chi2 * 5.0
    assert abs(bad.a_v - true_av) > 0.10 or np.max(np.abs(bad.x_fraction - true_x)) > 0.12


def test_g2_fit_spectrum_redshift_config():
    z = 0.02
    wave, bases, true_x, true_av, config, flux_rest, err_rest = _rest_mock()
    flux_obs = resample_to(wave, flux_rest, wave / (1.0 + z)) / (1.0 + z)
    err_obs = resample_to(wave, err_rest, wave / (1.0 + z)) / (1.0 + z)
    config.redshift = z
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux_obs, err_obs, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_g2_wave_frame_air_recovers():
    wave, bases, true_x, true_av, config, flux_rest, err_rest = _rest_mock()
    flux_air = resample_to(wave, flux_rest, vacuum_to_air(wave))
    err_air = resample_to(wave, err_rest, vacuum_to_air(wave))
    config.wave_frame = "air"
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux_air, err_air, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)
