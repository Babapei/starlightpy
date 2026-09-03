import json
from pathlib import Path

import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum, load_fit_result, save_fit_result
from starlightpy.simulate import default_absorption_bases, mock_observation


def _i2_fit():
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
    with pytest.warns(UserWarning):
        result = fit_spectrum(wave, flux, err, bases, config=config)
    return result


def _assert_wave_error(original, loaded):
    assert loaded.wavelength is not None and original.wavelength is not None
    assert loaded.error is not None and original.error is not None
    np.testing.assert_allclose(loaded.wavelength, original.wavelength)
    np.testing.assert_allclose(loaded.error, original.error)
    np.testing.assert_allclose(loaded.model, original.model)
    assert loaded.wavelength.shape == loaded.model.shape
    assert loaded.error.shape == loaded.model.shape


def test_i2_npz_roundtrip_keeps_wavelength_and_error(tmp_path: Path):
    result = _i2_fit()
    path = tmp_path / "result.npz"
    save_fit_result(path, result)
    loaded = load_fit_result(path)
    _assert_wave_error(result, loaded)


def test_i2_json_roundtrip_keeps_wavelength_and_error(tmp_path: Path):
    result = _i2_fit()
    path = tmp_path / "result.json"
    save_fit_result(path, result)
    loaded = load_fit_result(path)
    _assert_wave_error(result, loaded)
    data = json.loads(path.read_text(encoding="utf-8"))
    assert data["version"] == 2


def test_i2_v1_payload_without_wavelength_loads(tmp_path: Path):
    result = _i2_fit()
    path = tmp_path / "old.json"
    save_fit_result(path, result)
    data = json.loads(path.read_text(encoding="utf-8"))
    data["version"] = 1
    data.pop("wavelength", None)
    data.pop("error", None)
    path.write_text(json.dumps(data), encoding="utf-8")
    loaded = load_fit_result(path)
    assert loaded.wavelength is None
    assert loaded.error is None
    np.testing.assert_allclose(loaded.x, result.x)
    np.testing.assert_allclose(loaded.model, result.model)


def test_i2_rejects_fortran_out(tmp_path: Path):
    result = _i2_fit()
    with pytest.raises(ValueError, match=r"\.out"):
        save_fit_result(tmp_path / "galaxy.out", result)
    with pytest.raises(ValueError, match=r"\.out"):
        load_fit_result(tmp_path / "galaxy.out")
