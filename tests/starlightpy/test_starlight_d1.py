from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parents[2]


def test_d1_develop_page_keeps_working_readme():
    text = (ROOT / "docs" / "DEVELOP.md").read_text(encoding="utf-8")
    assert "PLAN.md" in text
    assert "阶段 I" in text
    assert "按" in text and "开发" in text


def test_d1_product_readme_install_and_example_figure():
    text = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "git+https://github.com/Babapei/starlight.git" in text
    assert "default_absorption_bases" in text
    assert "to_rest_frame" in text
    assert "loaded.flux" in text
    assert "docs/figures/ngc3522_emiles.png" in text
    assert "不是" in text and ("标准" in text or "真理" in text)
    assert (ROOT / "docs" / "figures" / "ngc3522_emiles.png").is_file()
    assert (ROOT / "examples" / "fit_ngc3522_emiles.py").is_file()


def test_d1_product_readme_synthetic_snippet_runs():
    """The five-minute example in README must execute, not just look copy-pasteable."""
    code = r"""
import numpy as np
from starlightpy import (
    FitConfig,
    apply_mask,
    default_absorption_bases,
    fit_spectrum,
    load_fit_result,
    mock_observation,
    optical_emission_mask_regions,
    save_fit_result,
)
wave = np.arange(3800.0, 5601.0, 2.0)
bases = default_absorption_bases(wave)
flux, error = mock_observation(
    wave, bases, [0.50, 0.35, 0.15], 0.30, 80.0, 130.0,
    config=FitConfig(search_kinematics=False),
    snr=40.0,
)
config = FitConfig(
    search_kinematics=True,
    refine_kinematics=True,
    a_v_bounds=(0.0, 0.6),
    a_v_step=0.1,
    v_bounds=(50.0, 125.0),
    v_step=25.0,
    sigma_bounds=(100.0, 175.0),
    sigma_step=25.0,
)
good = apply_mask(wave, optical_emission_mask_regions())
result = fit_spectrum(wave, flux, error, bases, mask=good, config=config)
assert result.flux is not None and result.model is not None
"""
    subprocess.check_call([sys.executable, "-c", code], cwd=str(ROOT))
