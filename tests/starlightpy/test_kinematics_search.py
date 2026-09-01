import numpy as np

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _b2_config() -> FitConfig:
    return FitConfig(
        search_kinematics=True,
        v_bounds=(0.0, 200.0),
        v_step=50.0,
        sigma_bounds=(50.0, 200.0),
        sigma_step=50.0,
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
    )


def _b2_truth():
    wave = np.arange(3800.0, 5601.0, 1.0)
    bases = default_absorption_bases(wave)
    x = np.array([0.50, 0.35, 0.15])
    return wave, bases, x, 0.30, 100.0, 150.0


def test_search_kinematics_tiny_grid_does_not_raise():
    wave = np.arange(4010.0, 4065.0)
    flux = np.ones(wave.size)
    err = np.ones(wave.size)
    bases = np.column_stack([flux, flux * 0.8])
    config = FitConfig(
        search_kinematics=True,
        v_bounds=(0.0, 0.0),
        v_step=50.0,
        sigma_bounds=(50.0, 50.0),
        sigma_step=50.0,
        a_v_bounds=(0.0, 0.1),
        a_v_step=0.1,
    )
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.v0_kms == 0.0
    assert result.sigma_kms == 50.0


def test_search_kinematics_noiseless_recovers_on_grid():
    wave, bases, true_x, true_av, true_v, true_sig = _b2_truth()
    config = _b2_config()
    flux, err = mock_observation(
        wave, bases, true_x, true_av, true_v, true_sig, config=config, snr=None
    )
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.v0_kms - true_v) <= 50.0
    assert abs(result.sigma_kms - true_sig) <= 50.0
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_search_kinematics_snr30_recovers_loosely():
    wave, bases, true_x, true_av, true_v, true_sig = _b2_truth()
    config = _b2_config()
    flux, err = mock_observation(
        wave,
        bases,
        true_x,
        true_av,
        true_v,
        true_sig,
        config=config,
        snr=30.0,
        rng=np.random.default_rng(1),
    )
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.v0_kms - true_v) <= 100.0
    assert abs(result.sigma_kms - true_sig) <= 100.0
    assert abs(result.a_v - true_av) <= 0.20
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.20)
