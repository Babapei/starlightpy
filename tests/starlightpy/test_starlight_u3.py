import numpy as np

from starlightpy import (
    FitConfig,
    default_absorption_bases,
    fit_spectrum,
    mock_observation,
    resample_to,
    to_rest_frame,
)


def _truth():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    true_v = 80.0
    true_s = 130.0
    inner = FitConfig(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        search_kinematics=False,
        v0_kms=0.0,
        sigma_kms=0.0,
    )
    flux, err = mock_observation(
        wave, bases, true_x, true_av, true_v, true_s, config=inner, snr=None
    )
    return wave, bases, true_x, true_av, true_v, true_s, flux, err


def _search_config(**kwargs) -> FitConfig:
    return FitConfig(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        search_kinematics=True,
        refine_kinematics=True,
        v_bounds=(50.0, 125.0),
        v_step=25.0,
        sigma_bounds=(100.0, 175.0),
        sigma_step=25.0,
        **kwargs,
    )


def _assert_recovered(result, true_x, true_av, true_v, true_s):
    assert abs(result.a_v - true_av) <= 0.10
    assert abs(result.v0_kms - true_v) <= 50.0
    assert abs(result.sigma_kms - true_s) <= 50.0
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_u3_redshift_config_on_rest_grid_recovers():
    wave, bases, true_x, true_av, true_v, true_s, flux_rest, err_rest = _truth()
    z = 0.02
    flux_obs = resample_to(wave, flux_rest, wave / (1.0 + z)) / (1.0 + z)
    err_obs = resample_to(wave, err_rest, wave / (1.0 + z)) / (1.0 + z)
    config = _search_config(redshift=z, wave_frame="vacuum")
    result = fit_spectrum(wave, flux_obs, err_obs, bases, config=config)
    _assert_recovered(result, true_x, true_av, true_v, true_s)


def test_u3_to_rest_frame_then_redshift_zero_recovers():
    wave, bases, true_x, true_av, true_v, true_s, flux_rest, err_rest = _truth()
    z = 0.02
    wave_obs = wave * (1.0 + z)
    flux_obs = flux_rest / (1.0 + z)
    err_obs = err_rest / (1.0 + z)
    wave_r, flux_r, err_r = to_rest_frame(wave_obs, flux_obs, z, err_obs)
    flux_fit = resample_to(wave_r, flux_r, wave)
    err_fit = resample_to(wave_r, err_r, wave)
    config = _search_config(redshift=0.0)
    result = fit_spectrum(wave, flux_fit, err_fit, bases, config=config)
    _assert_recovered(result, true_x, true_av, true_v, true_s)
