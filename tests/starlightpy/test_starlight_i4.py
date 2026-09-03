from pathlib import Path

import starlightpy


def test_i4_package_version_is_1_0_0():
    assert starlightpy.__version__ == "1.0.0"


def test_i4_ci_workflow_installs_dev_and_runs_pytest():
    path = Path(__file__).resolve().parents[2] / ".github" / "workflows" / "tests.yml"
    text = path.read_text(encoding="utf-8")
    assert 'pip install -e ".[dev]"' in text
    assert "pytest" in text
