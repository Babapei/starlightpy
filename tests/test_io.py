from pathlib import Path

import numpy as np

from starlightpy.io import apply_mask, load_base_master, load_base_spectra, load_mask, load_spectrum


def test_load_spectrum_and_mask(tmp_path: Path):
    spec = tmp_path / "g.cxt"
    spec.write_text("3500 1.0 0.1 0\n3501 1.1 0.1 0\n3502 1.2 0.1 2\n")
    wave, flux, err, flags = load_spectrum(spec)
    assert wave.tolist() == [3500.0, 3501.0, 3502.0]
    assert flags[-1] == 2

    mask = tmp_path / "m.mask"
    mask.write_text("1\n3500 3500.5 0.0\n")
    good = apply_mask(wave, load_mask(mask))
    assert good.tolist() == [False, True, True]


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
    np.testing.assert_allclose(matrix[:, 0], [1.0, 2.0])
    np.testing.assert_allclose(matrix[:, 1], [3.0, 4.0])
