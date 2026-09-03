import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum
from starlightpy.simulate import absorption_template, mock_observation


def _degenerate_six(wave):
    cols = [
        absorption_template(wave, [4861.0, 5172.0], line_depth=0.45),
        absorption_template(wave, [4861.0, 5172.0], line_depth=0.42),
        absorption_template(wave, [4861.0, 5172.0], line_depth=0.48),
        absorption_template(wave, [4920.0, 5270.0], line_depth=0.35),
        absorption_template(wave, [4920.0, 5270.0], line_depth=0.32),
        absorption_template(wave, [4920.0, 5270.0], line_depth=0.38),
    ]
    return np.column_stack(cols)


def _h1_setup():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = _degenerate_six(wave)
    true_x = np.array([0.40, 0.0, 0.0, 0.60, 0.0, 0.0])
    true_av = 0.30
    ages = np.array([1.0e7, 1.0e7, 1.0e7, 1.0e10, 1.0e10, 1.0e10])
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
    return wave, bases, true_x, true_av, ages, config, flux, err


def test_h1_unregularized_is_sparse_inside_bins():
    wave, bases, true_x, true_av, ages, config, flux, err = _h1_setup()
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config, template_ages=ages)
    young = result.x_fraction[:3]
    old = result.x_fraction[3:]
    assert abs(young.sum() - 0.40) <= 0.12
    assert abs(old.sum() - 0.60) <= 0.12
    assert young.max() - young.min() > 0.08
    assert np.all(result.x >= 0)
    assert abs(result.a_v - true_av) <= 0.10


def test_h1_age_bins_share_light():
    wave, bases, true_x, true_av, ages, config, flux, err = _h1_setup()
    config.regularize_x = "age_bins"
    config.age_bin_edges = (0.0, 1.0e8, 1.0e11)
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config, template_ages=ages)
    young = result.x_fraction[:3]
    old = result.x_fraction[3:]
    assert abs(young.sum() - 0.40) <= 0.12
    assert abs(old.sum() - 0.60) <= 0.12
    assert young.max() - young.min() <= 0.05
    assert old.max() - old.min() <= 0.05
    assert np.all(result.x >= 0)
    assert abs(result.a_v - true_av) <= 0.10


def test_h1_smooth_age_reduces_spread():
    wave, bases, true_x, true_av, ages, config, flux, err = _h1_setup()
    with pytest.warns(UserWarning):
        plain = fit_spectrum(wave, flux, err, bases, config=config, template_ages=ages)
    config.regularize_x = "smooth_age"
    config.regularize_strength = 1.0
    with pytest.warns(UserWarning):
        smooth = fit_spectrum(wave, flux, err, bases, config=config, template_ages=ages)
    spread_plain = plain.x_fraction[:3].max() - plain.x_fraction[:3].min()
    spread_smooth = smooth.x_fraction[:3].max() - smooth.x_fraction[:3].min()
    assert spread_smooth < spread_plain
    assert np.all(smooth.x >= 0)
    assert abs(smooth.x_fraction[:3].sum() - 0.40) <= 0.12
    assert abs(smooth.a_v - true_av) <= 0.10


def test_h1_missing_ages_raises():
    wave, bases, _true_x, _true_av, _ages, config, flux, err = _h1_setup()
    config.regularize_x = "age_bins"
    with pytest.raises(ValueError, match="template_ages"):
        fit_spectrum(wave, flux, err, bases, config=config, template_ages=None)
