from pathlib import Path

import numpy as np
import pytest

from starlightpy import FitConfig, fit_spectrum, load_fit_result, save_fit_result
from starlightpy.simulate import default_absorption_bases, mock_observation


def _h6_fit():
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


def _assert_roundtrip(original, loaded):
    np.testing.assert_allclose(loaded.x, original.x)
    np.testing.assert_allclose(loaded.x_fraction, original.x_fraction)
    assert loaded.a_v == original.a_v
    assert loaded.a_yv == original.a_yv
    assert loaded.v0_kms == original.v0_kms
    assert loaded.sigma_kms == original.sigma_kms
    assert loaded.chi2 == original.chi2
    np.testing.assert_allclose(loaded.model, original.model)
    assert loaded.model.shape == original.model.shape
    assert loaded.good is not None and original.good is not None
    assert loaded.good.shape == original.good.shape
    np.testing.assert_array_equal(loaded.good, original.good)
    assert loaded.config is not None and original.config is not None
    assert loaded.config.law == original.config.law
    assert loaded.config.a_v_bounds == original.config.a_v_bounds


def test_h6_npz_roundtrip(tmp_path: Path):
    result = _h6_fit()
    path = tmp_path / "result.npz"
    save_fit_result(path, result)
    loaded = load_fit_result(path)
    _assert_roundtrip(result, loaded)


def test_h6_json_roundtrip(tmp_path: Path):
    result = _h6_fit()
    path = tmp_path / "result.json"
    save_fit_result(path, result)
    loaded = load_fit_result(path)
    _assert_roundtrip(result, loaded)


def test_h6_fit_does_not_write_out(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.chdir(tmp_path)
    _h6_fit()
    assert list(tmp_path.glob("*.out")) == []


def test_h6_rejects_fortran_out(tmp_path: Path):
    result = _h6_fit()
    with pytest.raises(ValueError, match=r"\.out"):
        save_fit_result(tmp_path / "galaxy.out", result)
    with pytest.raises(ValueError, match=r"\.out"):
        load_fit_result(tmp_path / "galaxy.out")
