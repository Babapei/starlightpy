from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]


def test_r1_plan_reviewed_never_and_later_items():
    text = (ROOT / "docs" / "PLAN.md").read_text(encoding="utf-8")
    assert "审查记录" in text
    assert "改口径" in text
    assert "以后可做" in text
    assert "不是自动开工" in text
    for needle in ("O1", "O2", "O7", "O8", "O10", "波长相关 LSF"):
        assert needle in text, f"PLAN missing later item: {needle}"
    assert "退火" in text
    assert "backend=" in text
    assert "pad_losvd" in text
    assert "同一 `fit()` 里拟合发射线" not in text or "不再当禁令原文" in text


def test_r1_does_not_open_forbidden_config_names():
    config = (ROOT / "src" / "starlightpy" / "config.py").read_text(encoding="utf-8")
    assert "search_redshift" not in config
    assert "anneal_x" not in config
    assert "fit_emission" not in config
    assert "pad_losvd: bool = False" in config
