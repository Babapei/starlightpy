from pathlib import Path

import numpy as np
import pytest

from starlightpy import (
    FitConfig,
    combine_good,
    fit_spectrum,
    light_to_mass,
    light_weighted_age,
    load_fit_result,
    load_spectrum,
    optical_emission_mask_regions,
    save_fit_result,
)
from starlightpy.simulate import default_absorption_bases, mock_observation


def test_i3_cxt_mask_fit_products_save_roundtrip(tmp_path: Path):
    wave = np.arange(3800.0, 5601.0, 2.0)
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
    flux, err = mock_observation(wave, bases, true_x, true_av, 0.0, 0.0, config=config, snr=None)

    cxt = tmp_path / "galaxy.cxt"
    lines = [f"{wl:.4f} {fl:.8f} {e:.8f} 0\n" for wl, fl, e in zip(wave, flux, err)]
    cxt.write_text("".join(lines), encoding="utf-8")

    wave_in, flux_in, err_in, flags = load_spectrum(cxt)
    good = combine_good(wave_in, flags=flags, mask_regions=optical_emission_mask_regions())
    assert good.shape == wave_in.shape
    assert not bool(np.all(good))

    with pytest.warns(UserWarning):
        result = fit_spectrum(wave_in, flux_in, err_in, bases, mask=good, config=config)
    assert abs(result.a_v - true_av) <= 0.10
    np.testing.assert_allclose(result.x_fraction, true_x, atol=0.12)

    mass, mu = light_to_mass(result.x_fraction, np.array([1.0, 2.0, 3.0]))
    age = light_weighted_age(result.x_fraction, np.array([1.0e7, 1.0e8, 1.0e10]))
    assert mass.shape == true_x.shape
    assert np.isclose(mu.sum(), 1.0)
    assert age > 0

    out_path = tmp_path / "fit.npz"
    save_fit_result(out_path, result)
    loaded = load_fit_result(out_path)
    assert loaded.wavelength is not None
    assert loaded.wavelength.shape == loaded.model.shape
    np.testing.assert_allclose(loaded.wavelength, result.wavelength)
    assert list(tmp_path.glob("*.out")) == []


def test_i3_light_to_mass_requires_ml():
    with pytest.raises(ValueError, match="mass_to_light"):
        light_to_mass(np.array([0.50, 0.35, 0.15]), mass_to_light=None)
