from pathlib import Path


def test_u1_readme_file_path_and_redshift_contract():
    text = (Path(__file__).resolve().parents[2] / "README.md").read_text(encoding="utf-8")
    assert "to_rest_frame" in text
    assert "load_sdss_fits" in text or "load_spectrum" in text
    assert "redshift=0" in text or "redshift=0.0" in text
    assert "静止系" in text
