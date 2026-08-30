from pathlib import Path

from easyppxf.io import load_spectrum


def test_load_two_and_three_columns(tmp_path: Path):
    path = tmp_path / "s.txt"
    path.write_text("# wl flux err\n4000 1.0 0.1\n4001 1.1\n")
    wave, flux, err = load_spectrum(path)
    assert list(wave) == [4000.0, 4001.0]
    assert list(flux) == [1.0, 1.1]
    assert err[0] == 0.1
    assert err[1] == 1.0
