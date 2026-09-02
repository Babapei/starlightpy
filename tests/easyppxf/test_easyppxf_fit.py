import numpy as np
import pytest

ppxf = pytest.importorskip("ppxf")

from easyppxf import fit_spectrum


def _absorption_template(wave: np.ndarray, line_centers, width: float = 8.0) -> np.ndarray:
    spec = np.ones_like(wave)
    for center in line_centers:
        spec -= 0.4 * np.exp(-0.5 * ((wave - center) / width) ** 2)
    return np.clip(spec, 0.2, None)


def test_ppxf_recovers_near_zero_velocity_on_matching_template():
    wave_gal = np.arange(4800.0, 5400.0, 1.0)
    wave_temp = np.arange(4300.0, 5900.0, 1.0)
    templates = np.column_stack(
        [
            _absorption_template(wave_temp, [4861.0, 5175.0]),
            _absorption_template(wave_temp, [5000.0, 5270.0]),
            _absorption_template(wave_temp, [4920.0, 5320.0]),
        ]
    )
    true = np.interp(wave_gal, wave_temp, templates[:, 0])
    rng = np.random.default_rng(0)
    flux = true + rng.normal(0.0, 0.005, size=true.size)
    err = np.full_like(flux, 0.005)

    result = fit_spectrum(
        wave_gal,
        flux,
        templates,
        wave_temp,
        error=err,
        mask_emission=False,
        start=(0.0, 80.0),
        degree=0,
        moments=2,
    )
    assert abs(result.velocity) < 80.0
    assert 10.0 < result.sigma < 250.0
    assert result.bestfit.shape == result.galaxy.shape
    assert "Cappellari" in result.cite
