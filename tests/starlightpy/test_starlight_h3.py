import numpy as np
import pytest

from starlightpy import FitConfig, apply_losvd, fit_spectrum, losvd_design
from starlightpy.extinction import get_extinction_curve
from starlightpy.model import build_model, normalize_at, normalize_bases
from starlightpy.simulate import absorption_template


def _edge_bases(wave):
    return np.column_stack(
        [
            absorption_template(wave, [3810.0, 3850.0]),
            absorption_template(wave, [3825.0, 3880.0], line_depth=0.35),
            absorption_template(wave, [3840.0, 3910.0], line_width=8.0, line_depth=0.40),
        ]
    )


def _padded_mock():
    wave = np.arange(3800.0, 5201.0, 1.0)
    bases = _edge_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    true_av = 0.30
    v0, sig = 100.0, 150.0
    config = FitConfig(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        v0_kms=v0,
        sigma_kms=sig,
        search_kinematics=False,
        pad_losvd=True,
        losvd_oversample=2,
    )
    bases_n = normalize_bases(wave, bases, config.norm_window)
    q = get_extinction_curve(wave, law=config.law, r_v=config.r_v)
    win = (wave >= config.norm_window[0]) & (wave <= config.norm_window[1])
    q0 = float(np.median(q[win]))
    mixed = build_model(true_x, bases_n, q, true_av, q_lambda0=q0)
    model = apply_losvd(wave, mixed, v0, sig, pad=True, oversample=2)
    _, scale = normalize_at(wave, model, config.norm_window)
    flux = model / scale
    err = np.full_like(flux, max(0.01 * float(np.median(np.abs(flux[win]))), 1e-8))
    return wave, bases, true_x, true_av, config, flux, err


def test_h3_padded_losvd_recovers():
    wave, bases, true_x, true_av, config, flux, err = _padded_mock()
    result = fit_spectrum(wave, flux, err, bases, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)


def test_h3_padding_chi2_not_worse_than_default():
    wave, bases, true_x, true_av, config, flux, err = _padded_mock()
    padded = fit_spectrum(wave, flux, err, bases, config=config)
    plain = FitConfig(
        a_v_bounds=config.a_v_bounds,
        a_v_step=config.a_v_step,
        v0_kms=config.v0_kms,
        sigma_kms=config.sigma_kms,
        search_kinematics=False,
    )
    unpadded = fit_spectrum(wave, flux, err, bases, config=plain)
    assert padded.chi2 <= unpadded.chi2 + 1e-8


def test_h3_losvd_design_precompute_matches():
    wave = np.arange(4000.0, 4101.0, 1.0)
    bases = np.column_stack(
        [
            absorption_template(wave, [4015.0, 4055.0]),
            absorption_template(wave, [4020.0, 4065.0], line_depth=0.35),
        ]
    )
    a = losvd_design(wave, bases, 80.0, 120.0, pad=True, oversample=2)
    b = losvd_design(wave, bases, 80.0, 120.0, pad=True, oversample=2)
    np.testing.assert_allclose(a, b)
    from starlightpy.kinematics import LosvdDesignCache

    cache = LosvdDesignCache()
    c1 = cache.design(wave, bases, 80.0, 120.0, pad=True, oversample=2)
    c2 = cache.design(wave, bases, 80.0, 120.0, pad=True, oversample=2)
    np.testing.assert_allclose(c1, c2)
    np.testing.assert_allclose(c1, a)
    assert len(cache._store) == 1
