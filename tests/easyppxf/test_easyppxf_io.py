from pathlib import Path

import numpy as np
import pytest

from easyppxf.io import load_sdss_fits, load_spectrum


def test_load_two_and_three_columns(tmp_path: Path):
    path = tmp_path / "s.txt"
    path.write_text("# wl flux err\n4000 1.0 0.1\n4001 1.1\n")
    wave, flux, err = load_spectrum(path)
    assert list(wave) == [4000.0, 4001.0]
    assert list(flux) == [1.0, 1.1]
    assert err[0] == 0.1
    assert err[1] == 1.0


def test_easyppxf_load_sdss_fits_is_independent(tmp_path: Path):
    pytest.importorskip("astropy")
    import easyppxf.io as ppxf_io
    import starlightpy.io as star_io

    assert ppxf_io.load_sdss_fits is not star_io.load_sdss_fits

    from astropy.table import Table

    wave = np.arange(5000.0, 5010.0)
    flux = np.ones(wave.size)
    table = Table([flux, np.log10(wave)], names=("flux", "loglam"))
    path = tmp_path / "spec.fits"
    table.write(path, format="fits")
    out_wave, out_flux, out_err, good = load_sdss_fits(path)
    np.testing.assert_allclose(out_wave, wave)
    np.testing.assert_allclose(out_flux, flux)
    assert np.all(good)
    assert out_err.shape == flux.shape
