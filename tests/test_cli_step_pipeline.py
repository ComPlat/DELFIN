"""Tests for delfin-step and delfin-pipeline CLI entry points."""

from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any, Optional

import pytest

from delfin.tools._base import StepAdapter
from delfin.tools._types import StepResult, StepStatus


# ── Fake adapter (same as test_tools_layer.py) ────────────────────────

class FakeAdapter(StepAdapter):
    name = "fake_step"
    description = "Fake adapter for CLI testing"
    produces_geometry = True

    def validate_params(self, **kwargs: Any) -> None:
        if kwargs.get("fail_validation"):
            raise ValueError("validation failed on purpose")

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        if kwargs.get("fail"):
            return self._make_result(self.name, StepStatus.FAILED, work_dir, start, error="forced failure")
        out_xyz = work_dir / "output.xyz"
        out_xyz.write_text("1\nfake\nH  0.0  0.0  0.0\n")
        data = {"cores": cores, **{k: v for k, v in kwargs.items() if not k.startswith("_")}}
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start,
                                 geometry=out_xyz, data=data)


class FakeAnalysis(StepAdapter):
    name = "fake_analysis"
    description = "Fake analysis"
    produces_geometry = False

    def execute(self, work_dir: Path, *, geometry: Optional[Path] = None,
                cores: int = 1, **kwargs: Any) -> StepResult:
        start = time.monotonic()
        return self._make_result(self.name, StepStatus.SUCCESS, work_dir, start,
                                 data={"result": 42})


@pytest.fixture(autouse=True)
def _register_fakes(monkeypatch):
    from delfin.tools import _registry
    old_reg = dict(_registry._REGISTRY)
    old_disc = _registry._DISCOVERED
    _registry._DISCOVERED = True
    _registry._REGISTRY.clear()
    _registry.register(FakeAdapter())
    _registry.register(FakeAnalysis())
    yield
    _registry._REGISTRY.clear()
    _registry._REGISTRY.update(old_reg)
    _registry._DISCOVERED = old_disc


# ======================================================================
#  delfin-step CLI
# ======================================================================

class TestCliStep:
    def test_list(self, capsys):
        from delfin.cli_step import main
        rc = main(["--list"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "fake_step" in out
        assert "fake_analysis" in out

    def test_no_args_shows_help(self, capsys):
        from delfin.cli_step import main
        rc = main([])
        assert rc == 1

    def test_run_step(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        from delfin.cli_step import main
        rc = main(["fake_step", "--geometry", str(xyz), "--work-dir", str(tmp_path / "work")])
        assert rc == 0
        out = capsys.readouterr().out
        assert "[OK]" in out

    def test_run_step_json(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        from delfin.cli_step import main
        rc = main(["fake_step", "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--json"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["ok"] is True
        assert data["step_name"] == "fake_step"

    def test_extra_kwargs(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        from delfin.cli_step import main
        rc = main(["fake_step", "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"),
                    "--json", "--charge=2", "--method", "B3LYP"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["data"]["charge"] == 2
        assert data["data"]["method"] == "B3LYP"

    def test_unknown_step(self, capsys):
        from delfin.cli_step import main
        rc = main(["nonexistent"])
        assert rc == 1

    def test_failed_step(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        from delfin.cli_step import main
        rc = main(["fake_step", "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--fail"])
        assert rc == 1
        out = capsys.readouterr().out
        assert "[FAILED]" in out


# ======================================================================
#  delfin-pipeline CLI
# ======================================================================

class TestCliPipeline:
    def _write_yaml(self, tmp_path, content):
        p = tmp_path / "pipeline.yaml"
        p.write_text(content)
        return str(p)

    def test_dry_run(self, tmp_path, capsys):
        yaml_path = self._write_yaml(tmp_path, """
name: test
steps:
  - step: fake_step
  - step: fake_analysis
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--dry-run"])
        assert rc == 0
        out = capsys.readouterr().out
        assert "Dry run" in out
        assert "2" in out  # 2 steps

    def test_run_pipeline(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: test
steps:
  - step: fake_step
  - step: fake_analysis
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work")])
        assert rc == 0
        out = capsys.readouterr().out
        assert "OK" in out

    def test_run_with_defaults(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: test
defaults:
  charge: 2
steps:
  - step: fake_step
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--json"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["ok"]
        assert data["steps"][0]["data"]["charge"] == 2

    def test_run_with_branches(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: test
steps:
  - step: fake_step
branches:
  alpha:
    - step: fake_step
    - step: fake_analysis
  beta:
    - step: fake_analysis
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--json"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["ok"]
        assert "alpha" in data["branches"]
        assert "beta" in data["branches"]

    def test_template_with_params(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: tpl
template: true
steps:
  - step: fake_step
    charge: "{charge}"
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"),
                    "--json", "--param", "charge=2"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["steps"][0]["data"]["charge"] == 2

    def test_missing_yaml(self, capsys):
        from delfin.cli_pipeline import main
        rc = main(["/no/such/file.yaml"])
        assert rc == 1
        err = capsys.readouterr().err
        assert "not found" in err

    def test_invalid_yaml(self, tmp_path, capsys):
        yaml_path = self._write_yaml(tmp_path, "just a string")
        from delfin.cli_pipeline import main
        rc = main([yaml_path])
        assert rc == 1

    def test_yaml_missing_name(self, tmp_path, capsys):
        yaml_path = self._write_yaml(tmp_path, """
steps:
  - step: fake_step
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path])
        assert rc == 1
        err = capsys.readouterr().err
        assert "name" in err

    # -- YAML flow control tests --

    def test_yaml_retry_step(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: retry_test
steps:
  - step: fake_step
    type: retry
    max_attempts: 2
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--json"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["ok"]

    def test_yaml_checkpoint_step(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: cp_test
steps:
  - step: fake_step
  - type: checkpoint
  - step: fake_analysis
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--json"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["ok"]
        assert len(data["steps"]) == 3

    def test_yaml_if_step(self, tmp_path, capsys):
        xyz = tmp_path / "mol.xyz"
        xyz.write_text("1\ntest\nH  0.0  0.0  0.0\n")
        yaml_path = self._write_yaml(tmp_path, """
name: if_test
steps:
  - step: fake_step
  - step: fake_analysis
    type: if
    condition: "last.ok"
""")
        from delfin.cli_pipeline import main
        rc = main([yaml_path, "--geometry", str(xyz),
                    "--work-dir", str(tmp_path / "work"), "--json"])
        assert rc == 0
        data = json.loads(capsys.readouterr().out)
        assert data["ok"]
        assert len(data["steps"]) == 2
