from pathlib import Path


def test_i5_readme_lists_optional_flags_and_save():
    text = (Path(__file__).resolve().parents[2] / "README.md").read_text(encoding="utf-8")
    required = (
        "regularize_x",
        "error_method",
        "pad_losvd",
        "fit_ayv",
        "dust:F99",
        "save_fit_result",
        ".[dust]",
        ".[fits]",
    )
    missing = [token for token in required if token not in text]
    assert missing == []
    assert "1.0.0" in text
    develop = (Path(__file__).resolve().parents[2] / "docs" / "DEVELOP.md").read_text(
        encoding="utf-8"
    )
    assert "阶段 I" in develop
