import numpy as np
import pytest

from starlightpy import (
    FitConfig,
    apply_mask,
    fit_spectrum,
    match_instrumental_fwhm,
    optical_emission_mask_regions,
)
from starlightpy.simulate import default_absorption_bases, mock_observation


def _noiseless_mock(dw=1.0):
    wave = np.arange(3800.0, 5601.0, dw)
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
    flux, err = mock_observation(
        wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None
    )
    return wave, bases, true_x, true_av, config, flux, err


def test_g3_lsf_match_recovers():
    wave, bases, true_x, true_av, config, flux, err = _noiseless_mock()
    flux4 = match_instrumental_fwhm(wave, flux, fwhm_data=4.0, fwhm_template=2.0)
    config.fwhm_data = 4.0
    config.fwhm_template = 2.0
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux4, err, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_g3_skipping_lsf_is_worse():
    wave, bases, true_x, true_av, config, flux, err = _noiseless_mock()
    flux4 = match_instrumental_fwhm(wave, flux, fwhm_data=4.0, fwhm_template=2.0)
    matched = FitConfig(
        a_v_bounds=config.a_v_bounds,
        a_v_step=config.a_v_step,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
        fwhm_data=4.0,
        fwhm_template=2.0,
    )
    unmatched = FitConfig(
        a_v_bounds=config.a_v_bounds,
        a_v_step=config.a_v_step,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
    )
    with pytest.warns(UserWarning):
        good = fit_spectrum(wave, flux4, err, bases, config=matched)
        bad = fit_spectrum(wave, flux4, err, bases, config=unmatched)
    assert bad.chi2 > good.chi2
    assert abs(good.a_v - true_av) <= 0.10
    np.testing.assert_allclose(good.x_fraction, true_x, atol=0.12)


def test_g3_cannot_deconvolve():
    wave = np.arange(3800.0, 4001.0, 1.0)
    y = np.exp(-0.5 * ((wave - 3900.0) / 4.0) ** 2)
    out = match_instrumental_fwhm(wave, y, fwhm_data=2.0, fwhm_template=4.0)
    np.testing.assert_array_equal(out, y)


def test_g3_halpha_masked():
    wave = np.arange(6500.0, 6621.0, 1.0)
    good = apply_mask(wave, optical_emission_mask_regions())
    i_ha = int(np.argmin(np.abs(wave - 6563.0)))
    assert not bool(good[i_ha])
