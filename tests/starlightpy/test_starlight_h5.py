import sys

import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum, get_extinction_curve
from starlightpy.simulate import default_absorption_bases, mock_observation


def _h5_config(**kwargs) -> FitConfig:
    cfg = dict(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
    )
    cfg.update(kwargs)
    return FitConfig(**cfg)


def test_h5_dust_f99_recovers():
    pytest.importorskip("dust_extinction")
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    config = _h5_config(law="dust:F99")
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_h5_builtin_ccm_is_worse_on_f99_truth():
    pytest.importorskip("dust_extinction")
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    f99 = _h5_config(law="dust:F99")
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=f99, snr=None)
    with pytest.warns(UserWarning):
        matched = fit_spectrum(wave, flux, err, bases, config=f99)
    ccm = _h5_config(law="CCM")
    with pytest.warns(UserWarning):
        mismatched = fit_spectrum(wave, flux, err, bases, config=ccm)
    assert mismatched.chi2 > matched.chi2 + 1e-8 or abs(mismatched.a_v - true_av) > 0.05


def test_h5_default_ccm_still_recovers():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    config = _h5_config()
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    assert config.law == "CCM"
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_h5_missing_package_is_explicit(monkeypatch):
    for key in list(sys.modules):
        if key == "dust_extinction" or key.startswith("dust_extinction."):
            monkeypatch.setitem(sys.modules, key, None)
    with pytest.raises((ImportError, ModuleNotFoundError, ValueError), match="dust_extinction"):
        get_extinction_curve([5000.0], law="dust:F99")


def test_h5_unknown_dust_model_raises():
    pytest.importorskip("dust_extinction")
    with pytest.raises(ValueError, match="Unknown dust_extinction model"):
        get_extinction_curve([5000.0], law="dust:NOTALAW")
