from pathlib import Path


def test_g7_readme_usage_and_checklist():
    readme = Path(__file__).resolve().parents[2] / "README.md"
    text = readme.read_text(encoding="utf-8")
    assert "from starlightpy import" in text
    assert "fit_spectrum" in text
    for needle in ("静止系", "波长", "分辨率", "发射线", "简并"):
        assert needle in text, f"README missing checklist item: {needle}"
