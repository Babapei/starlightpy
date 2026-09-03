import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import default_absorption_bases, mock_observation


def _h4_setup(true_ayv: float):
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.20
    young = np.array([1.0, 0.0, 0.0])
    config = FitConfig(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
        fit_ayv=True,
        a_yv_bounds=(0.0, 0.8),
        a_yv_step=0.2,
    )
    flux, err = mock_observation(
        wave,
        bases,
        true_x,
        true_av,
        0.0,
        0.0,
        config=config,
        snr=None,
        a_yv=true_ayv,
        young_flags=young,
    )
    return wave, bases, true_x, true_av, young, config, flux, err


def test_h4_recovers_ayv_when_enabled():
    wave, bases, true_x, true_av, young, config, flux, err = _h4_setup(0.40)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config, young_flags=young)
    assert abs(result.a_v - true_av) <= 0.10
    assert abs(result.a_yv - 0.40) <= 0.15
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)
    assert np.all(result.x >= 0)


def test_h4_shared_av_is_worse_on_ayv_truth():
    wave, bases, true_x, true_av, young, config, flux, err = _h4_setup(0.40)
    with pytest.warns(UserWarning):
        with_ayv = fit_spectrum(wave, flux, err, bases, config=config, young_flags=young)
    off = FitConfig(
        a_v_bounds=config.a_v_bounds,
        a_v_step=config.a_v_step,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
        fit_ayv=False,
    )
    with pytest.warns(UserWarning):
        shared = fit_spectrum(wave, flux, err, bases, config=off)
    assert shared.a_yv == 0.0
    assert shared.chi2 > with_ayv.chi2 + 1e-8 or abs(shared.a_v - true_av) > 0.10


def test_h4_zero_ayv_stays_near_zero():
    wave, bases, true_x, true_av, young, config, flux, err = _h4_setup(0.0)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config, young_flags=young)
    assert abs(result.a_yv) <= config.a_yv_step + 1e-12
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_h4_fit_ayv_requires_young_flags():
    wave, bases, _x, _av, _young, config, flux, err = _h4_setup(0.40)
    with pytest.raises(ValueError, match="young_flags"):
        fit_spectrum(wave, flux, err, bases, config=config)
