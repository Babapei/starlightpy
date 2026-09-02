import gzip
from pathlib import Path

import numpy as np
import pytest

from starlightpy.io import (
    combine_good,
    good_from_flags,
    iter_ascii_spectra,
    load_base_master,
    load_base_spectra,
    load_cxt,
    load_mask,
    load_sdss_fits,
    load_spectrum,
)


def test_load_spectrum_and_mask(tmp_path: Path):
    spec = tmp_path / "g.cxt"
    spec.write_text("3500 1.0 0.1 0\n3501 1.1 0.1 0\n3502 1.2 0.1 2\n")
    wave, flux, err, flags = load_spectrum(spec)
    assert wave.tolist() == [3500.0, 3501.0, 3502.0]
    assert flags[-1] == 2
    assert good_from_flags(flags).tolist() == [True, True, False]

    mask = tmp_path / "m.mask"
    mask.write_text("1\n3500 3500.5 0.0\n")
    good = combine_good(wave, flags=flags, mask_regions=load_mask(mask))
    assert good.tolist() == [False, True, False]


def test_cxt_skips_npix_line(tmp_path: Path):
    spec = tmp_path / "g.cxt"
    spec.write_text("3\n3500 1.0 0.1 0\n3501 1.1 0.1 0\n3502 1.2 0.1 2\n")
    wave, flux, err, flags = load_cxt(spec)
    assert wave.tolist() == [3500.0, 3501.0, 3502.0]
    assert flags[-1] == 2


def test_load_gzipped_cxt_and_bases(tmp_path: Path):
    cxt = tmp_path / "g.cxt.gz"
    with gzip.open(cxt, "wt", encoding="utf-8") as handle:
        handle.write("2\n4000 1.0 0.1 0\n4001 1.1 0.1 0\n")
    wave, flux, err, flags = load_spectrum(cxt)
    assert wave.tolist() == [4000.0, 4001.0]

    with gzip.open(tmp_path / "a.spec.gz", "wt", encoding="utf-8") as handle:
        handle.write("4000 1.0\n4001 2.0\n")
    with gzip.open(tmp_path / "b.spec.gz", "wt", encoding="utf-8") as handle:
        handle.write("4000 3.0\n4001 4.0\n")
    master = tmp_path / "Base.N"
    master.write_text(
        "2\n"
        "a.spec 1e9 0.02 young 1.0 0 0.0\n"
        "b.spec 1e10 0.02 old 1.0 0 0.0\n"
    )
    bases = load_base_master(master)
    wave_b, matrix = load_base_spectra(tmp_path, bases)
    assert wave_b.tolist() == [4000.0, 4001.0]
    assert matrix[1, 1] == 4.0


def test_iter_ascii_spectra(tmp_path: Path):
    (tmp_path / "a.cxt").write_text("4000 1.0 0.1 0\n")
    (tmp_path / "b.cxt").write_text("4001 2.0 0.1 0\n")
    (tmp_path / "skip.txt").write_text("1 1 1 0\n")
    items = list(iter_ascii_spectra(tmp_path, pattern="*.cxt"))
    assert len(items) == 2
    names = [path.name for path, _spec in items]
    assert names == ["a.cxt", "b.cxt"]
    assert items[1][1][1][0] == 2.0


def test_load_base(tmp_path: Path):
    (tmp_path / "a.spec").write_text("4000 1.0\n4001 2.0\n")
    (tmp_path / "b.spec").write_text("4000 3.0\n4001 4.0\n")
    master = tmp_path / "Base.N"
    master.write_text(
        "2\n"
        "a.spec 1e9 0.02 young 1.0 0 0.0\n"
        "b.spec 1e10 0.02 old 1.0 0 0.0\n"
    )
    bases = load_base_master(master)
    wave, matrix = load_base_spectra(tmp_path, bases)
    assert wave.tolist() == [4000.0, 4001.0]
    assert matrix[0, 0] == 1.0
    assert matrix[1, 1] == 4.0


def test_starlightpy_sources_do_not_import_easyppxf():
    root = Path(__file__).resolve().parents[2] / "src" / "starlightpy"
    for path in root.glob("*.py"):
        text = path.read_text(encoding="utf-8")
        assert "import easyppxf" not in text
        assert "from easyppxf" not in text


def test_load_sdss_fits_table(tmp_path: Path):
    pytest.importorskip("astropy")
    from astropy.table import Table

    wave = np.arange(4000.0, 4010.0)
    flux = np.linspace(1.0, 2.0, wave.size)
    ivar = np.full(wave.size, 100.0)
    ivar[1] = 0.0
    and_mask = np.zeros(wave.size, dtype=np.int32)
    and_mask[2] = 1
    table = Table(
        [flux, np.log10(wave), ivar, and_mask],
        names=("flux", "loglam", "ivar", "and_mask"),
    )
    path = tmp_path / "spec.fits"
    table.write(path, format="fits")
    out_wave, out_flux, out_err, good = load_sdss_fits(path)
    np.testing.assert_allclose(out_wave, wave)
    np.testing.assert_allclose(out_flux, flux)
    assert not good[1]
    assert not good[2]
    assert good[0]
    np.testing.assert_allclose(out_err[0], 0.1)
