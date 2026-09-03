import json
from pathlib import Path

import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum, load_fit_result, save_fit_result
from starlightpy.simulate import default_absorption_bases, mock_observation


def _j2_fit():
    wave = np.arange(3800.0, 5601.0, 2.0)
    bases = default_absorption_bases(wave)
    true_x = np.array([0.50, 0.35, 0.15])
    config = FitConfig(
        a_v_bounds=(0.0, 0.6),
        a_v_step=0.1,
        v0_kms=0.0,
        sigma_kms=0.0,
        search_kinematics=False,
    )
    flux, err = mock_observation(wave, bases, true_x, 0.30, 0.0, 0.0, config=config, snr=None)
    with pytest.warns(UserWarning):
        return fit_spectrum(wave, flux, err, bases, config=config)


def test_j2_npz_roundtrip_keeps_flux(tmp_path: Path):
    result = _j2_fit()
    path = tmp_path / "result.npz"
    save_fit_result(path, result)
    loaded = load_fit_result(path)
    assert loaded.flux is not None and result.flux is not None
    np.testing.assert_allclose(loaded.flux, result.flux)
    np.testing.assert_allclose(loaded.model, result.model)
    np.testing.assert_allclose(loaded.wavelength, result.wavelength)
    assert loaded.flux.shape == loaded.model.shape


def test_j2_json_roundtrip_keeps_flux(tmp_path: Path):
    result = _j2_fit()
    path = tmp_path / "result.json"
    save_fit_result(path, result)
    loaded = load_fit_result(path)
    assert loaded.flux is not None
    np.testing.assert_allclose(loaded.flux, result.flux)
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["version"] == 3


def test_j2_v2_payload_without_flux_loads(tmp_path: Path):
    result = _j2_fit()
    path = tmp_path / "old.json"
    save_fit_result(path, result)
    data = json.loads(path.read_text(encoding="utf-8"))
    data["version"] = 2
    data.pop("flux", None)
    path.write_text(json.dumps(data), encoding="utf-8")
    loaded = load_fit_result(path)
    assert loaded.flux is None
    np.testing.assert_allclose(loaded.model, result.model)
    np.testing.assert_allclose(loaded.wavelength, result.wavelength)


def test_j2_rejects_fortran_out(tmp_path: Path):
    result = _j2_fit()
    with pytest.raises(ValueError, match=r"\.out"):
        save_fit_result(tmp_path / "galaxy.out", result)
