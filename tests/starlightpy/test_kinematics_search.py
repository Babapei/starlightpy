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


def _assert_off_grid(value, low, high, step):
    n = round((value - low) / step)
    node = low + n * step
    assert abs(value - node) > 1e-6, f"{value} landed on a grid node"


def test_b3_off_grid_recovers_v80_sigma130():
    """PLAN B3: truth not on the kinematic grid; noiseless, tolerances no looser than B2."""
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av, true_v, true_sig = 0.30, 80.0, 130.0
    config = FitConfig(
        search_kinematics=True,
        v_bounds=(50.0, 125.0),
        v_step=25.0,
        sigma_bounds=(100.0, 175.0),
        sigma_step=25.0,
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
    )
    _assert_off_grid(true_v, *config.v_bounds, config.v_step)
    _assert_off_grid(true_sig, *config.sigma_bounds, config.sigma_step)
    flux, err = mock_observation(
        wave, bases, true_x, true_av, true_v, true_sig, config=config, snr=None
    )
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.v0_kms - true_v) <= 25.0
    assert abs(result.sigma_kms - true_sig) <= 25.0
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_b3_off_grid_second_mix():
    """PLAN B3: a second (x, A_V, v, σ) set, still off-grid and noiseless."""
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.20, 0.55, 0.25])
    true_av, true_v, true_sig = 0.35, 120.0, 90.0
    config = FitConfig(
        search_kinematics=True,
        v_bounds=(50.0, 175.0),
        v_step=25.0,
        sigma_bounds=(50.0, 150.0),
        sigma_step=25.0,
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
    )
    _assert_off_grid(true_v, *config.v_bounds, config.v_step)
    _assert_off_grid(true_sig, *config.sigma_bounds, config.sigma_step)
    flux, err = mock_observation(
        wave, bases, true_x, true_av, true_v, true_sig, config=config, snr=None
    )
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.v0_kms - true_v) <= 25.0
    assert abs(result.sigma_kms - true_sig) <= 25.0
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)
