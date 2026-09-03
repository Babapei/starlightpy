from pathlib import Path


def test_j3_readme_says_loaded_flux_and_model_can_plot():
    text = (Path(__file__).resolve().parents[2] / "README.md").read_text(encoding="utf-8")
    assert "loaded.flux" in text
    assert "loaded.model" in text
    assert "loaded.wavelength" in text
