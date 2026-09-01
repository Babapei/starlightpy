import numpy as np

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _d_config(**kwargs) -> FitConfig:
    cfg = dict(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
    )
    cfg.update(kwargs)
    return FitConfig(**cfg)


def test_nsigma_clip_ignores_spike_and_recovers():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    config = _d_config(clip_nsigma=3.0)
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None)
    spike = int(np.argmin(np.abs(wave - 4500.0)))
    assert not (4010.0 <= wave[spike] <= 4060.0)
    flux = flux.copy()
    flux[spike] = float(np.median(flux)) * 40.0
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.n_clipped >= 1
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_ex0_drops_tiny_component():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.55, 0.40, 0.05])
    true_av = 0.30
    config = _d_config(x_min_keep=0.10)
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None)
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert result.x[2] == 0.0
    expected = np.array([0.55, 0.40, 0.0]) / 0.95
    np.testing.assert_allclose(result.x_fraction, expected, atol=0.12)
    assert abs(result.a_v - true_av) <= 0.10
    assert result.n_clipped == 0
