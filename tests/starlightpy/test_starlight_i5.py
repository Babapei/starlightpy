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
    assert "阶段 I" in text
    assert "1.0.0" in text
