"""SW2a phase 2: swarm.yaml save/load + Pre-PR display (red control).

The tab must offer save/load of the spec through SW1's own functions
(write_example_spec + load_spec) and show the Pre-PR check as the
COMMAND for the user (merge_plan) -- never executed by the tab.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.swarm import merge as _merge
from delfin.agent.swarm.spec import SpecError, load_spec
from delfin.dashboard import tab_swarm

VALID_SPEC = """\
goal: test goal
packages:
  - name: pkg1
    owned: ["delfin/dashboard/tab_swarm.py"]
    task: build
    roles:
      builder: true
"""


def test_save_spec_uses_sw1_writer(tmp_path: Path) -> None:
    target = tmp_path / "swarm.yaml"
    returned = tab_swarm.save_spec(VALID_SPEC, target)
    assert returned == target
    assert target.exists()
    # The file round-trips through SW1's own loader unchanged.
    spec = load_spec(target)
    assert spec["goal"] == "test goal"
    assert spec["packages"][0]["name"] == "pkg1"


def test_save_spec_refuses_overwrite_without_force(tmp_path: Path) -> None:
    target = tmp_path / "swarm.yaml"
    tab_swarm.save_spec(VALID_SPEC, target)
    with pytest.raises(SpecError):
        tab_swarm.save_spec(VALID_SPEC, target)


def test_load_spec_text_roundtrip(tmp_path: Path) -> None:
    target = tmp_path / "swarm.yaml"
    target.write_text(VALID_SPEC, encoding="utf-8")
    text = tab_swarm.load_spec_text(target)
    assert text == VALID_SPEC


def test_load_spec_text_reports_missing_file(tmp_path: Path) -> None:
    with pytest.raises(FileNotFoundError):
        tab_swarm.load_spec_text(tmp_path / "absent.yaml")


def test_pre_pr_display_is_the_command_not_an_execution(
        tmp_path: Path) -> None:
    plan = tab_swarm.pre_pr_display("pkg1", base="main")
    # Comes from SW1's merge_plan: commands as TEXT, executed always False.
    assert plan["executed"] is False
    assert plan["merge_command"] == "git merge swarm/pkg1-builder"
    assert plan["pre_pr_check"].startswith("delfin-agent prcheck")
    assert plan["pre_pr_check_status"] == "not_run"


def test_pre_pr_display_rejects_unsafe_package() -> None:
    with pytest.raises(_merge.TriageError):
        tab_swarm.pre_pr_display("bad name; rm -rf")
