"""Package SW1 phase 2b: the `delfin-agent swarm` CLI surface + status.

Control: tests/test_sw1_cli.py, 5 cases, all 5 red on the previous
commit (ModuleNotFoundError: delfin.agent.cli_swarm).
"""

from __future__ import annotations

import json
import textwrap

import pytest

from delfin.agent import cli_swarm
from delfin.agent.swarm import status as swarm_status
from delfin.agent.swarm.spec import SpecError


@pytest.fixture()
def swarm_repo(tmp_path, monkeypatch):
    (tmp_path / "delfin").mkdir()
    (tmp_path / "swarm.yaml").write_text(textwrap.dedent("""
        goal: Ship it
        packages:
          - name: core
            owned: ["delfin/**"]
            task: Build core
    """), encoding="utf-8")
    monkeypatch.chdir(tmp_path)
    return tmp_path


def test_init_writes_example_and_refuses_overwrite(swarm_repo, capsys):
    rc = cli_swarm.main(["init"])
    assert rc == 0
    assert (swarm_repo / "swarm.yaml").exists()
    # The shipped example itself must validate.
    from delfin.agent.swarm.launch import EXAMPLE_SPEC
    from delfin.agent.swarm.spec import load_spec, validate_spec
    import yaml
    spec = yaml.safe_load(EXAMPLE_SPEC)
    validate_spec(spec)
    rc = cli_swarm.main(["init"])
    assert rc == 1  # refuses to overwrite


def test_plan_validates_and_lists_members(swarm_repo, capsys):
    rc = cli_swarm.main(["plan", "Ship it"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "core-builder" in out
    assert "never bypass" in out


def test_plan_rejects_bad_spec_with_prose_error(swarm_repo, capsys):
    (swarm_repo / "swarm.yaml").write_text("goal: x\npackages: []\n",
                                           encoding="utf-8")
    rc = cli_swarm.main(["plan", "Ship it"])
    assert rc == 1
    err = capsys.readouterr().err
    assert "ERROR" in err and "packages" in err


def test_start_dry_run_starts_nothing(swarm_repo, capsys):
    rc = cli_swarm.main(["start", "--dry-run"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "Nothing started" in out


def test_status_json_matches_schema(swarm_repo, capsys, monkeypatch):
    monkeypatch.setattr(swarm_status, "_event_subscribers", [])
    rc = cli_swarm.main(["status", "--json"])
    assert rc == 0
    report = json.loads(capsys.readouterr().out)
    assert swarm_status.validate_report(report) == []
    assert report["packages"][0]["name"] == "core"
    assert report["packages"][0]["members"][0]["name"] == "core-builder"
