import numpy as np
import pytest

ppxf = pytest.importorskip("ppxf")

from easyppxf import fit_spectrum


def _absorption_template(wave: np.ndarray, line_centers, width: float = 8.0) -> np.ndarray:
    spec = np.ones_like(wave)
    for center in line_centers:
        spec -= 0.4 * np.exp(-0.5 * ((wave - center) / width) ** 2)
    return np.clip(spec, 0.2, None)


def test_u2_same_wavelength_does_not_hit_ppxf_assert():
    wave = np.arange(4800.0, 5400.0, 1.0)
    templates = np.column_stack(
        [
            _absorption_template(wave, [4861.0, 5175.0]),
            _absorption_template(wave, [5000.0, 5270.0]),
        ]
    )
    rng = np.random.default_rng(1)
    flux = templates[:, 0] + rng.normal(0.0, 0.005, size=wave.size)
    err = np.full_like(flux, 0.005)
    result = fit_spectrum(
        wave,
        flux,
        templates,
        wave,
        error=err,
        mask_emission=False,
        start=(0.0, 80.0),
        degree=0,
        moments=2,
    )
    assert np.isfinite(result.velocity)
    assert np.isfinite(result.sigma)
    assert result.bestfit.shape == result.galaxy.shape


def test_u2_short_templates_raise_value_error():
    wave_gal = np.arange(4800.0, 5400.0, 1.0)
    wave_temp = np.arange(5100.0, 5111.0, 1.0)
    templates = _absorption_template(wave_temp, [5105.0])
    flux = np.ones_like(wave_gal)
    err = np.full_like(flux, 0.01)
    with pytest.raises(ValueError, match=r"wavelength|2900|km/s"):
        fit_spectrum(
            wave_gal,
            flux,
            templates,
            wave_temp,
            error=err,
            mask_emission=False,
            degree=0,
        )
