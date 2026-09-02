import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.optimize import metropolis_anneal
from starlightpy.simulate import default_absorption_bases, mock_observation


def _off_grid_setup():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av, true_v, true_sig = 0.30, 80.0, 130.0
    config = FitConfig(
        search_kinematics=True,
        refine_kinematics=False,
        v_bounds=(50.0, 125.0),
        v_step=25.0,
        sigma_bounds=(100.0, 175.0),
        sigma_step=25.0,
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
    )
    flux, err = mock_observation(
        wave, bases, true_x, true_av, true_v, true_sig, config=config, snr=None
    )
    return wave, bases, true_x, true_av, true_v, true_sig, config, flux, err


def test_g5_coarse_stays_on_nearby_node():
    wave, bases, true_x, true_av, true_v, true_sig, config, flux, err = _off_grid_setup()
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.v0_kms - true_v) <= 25.0
    assert abs(result.sigma_kms - true_sig) <= 25.0
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_g5_refine_closer_or_chi2_not_worse():
    wave, bases, true_x, true_av, true_v, true_sig, config, flux, err = _off_grid_setup()
    coarse = fit_spectrum(wave, flux, err, bases, config=config)
    config.refine_kinematics = True
    refined = fit_spectrum(wave, flux, err, bases, config=config)
    dv_c = abs(coarse.v0_kms - true_v)
    ds_c = abs(coarse.sigma_kms - true_sig)
    dv_r = abs(refined.v0_kms - true_v)
    ds_r = abs(refined.sigma_kms - true_sig)
    closer = (dv_r < dv_c) and (ds_r < ds_c)
    chi2_ok = refined.chi2 <= coarse.chi2 + 1e-8
    assert closer or chi2_ok
    assert abs(refined.a_v - true_av) <= 0.10
    np.testing.assert_allclose(refined.x_fraction, true_x, atol=0.12)


def test_g5_metropolis_still_closed():
    with pytest.raises(NotImplementedError):
        metropolis_anneal()
