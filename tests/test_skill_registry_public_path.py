"""Integration control: the evidence gate on the PUBLIC call path.

Calls publish_skill/pull_skills exactly the way delfin's dashboard
(tab_agent /skills push|pull) calls them: keyword connection args from
the transfer settings, workspace from the repo dir, the transport
replaced by the existing _default_run stand-in (no real transfer).
Red on the previous commit (publish without evidence succeeded there).
"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from delfin.agent import skill_registry as sr  # noqa: E402


_REPO = Path(__file__).resolve().parents[1]
_CONN = {"host": "login", "user": "grp", "remote_path": "/archive"}


def _stamped_green_run() -> list[dict]:
    from delfin.agent import evidence_freshness as ef
    run = {"command": "tests/test_skill_registry.py", "exit_code": 0,
           "status": "ok", "passed": 1, "failed": 0}
    return [ef.stamp(run, _REPO)]


def _ok_run(calls):
    def run(cmd):
        calls.append(cmd)
        return 0, "", ""
    return run


def test_public_path_publish_without_evidence_refused():
    """The dashboard's /skills push call shape (no evidence argument):
    must now refuse -- nothing leaves the machine."""
    calls: list = []
    ok, msg = sr.publish_skill("casscf-setup", workspace=_REPO,
                               run_fn=_ok_run(calls), **_CONN)
    assert ok is False and "evidence" in msg.lower()
    assert calls == []


def test_public_path_publish_with_evidence_builds_ssh_commands():
    """Same call shape plus a valid evidence record: the mkdir + rsync
    commands are built through the proven ssh_transfer_jobs builders
    (stand-in transport, no real transfer)."""
    calls: list = []
    ok, msg = sr.publish_skill(
        "casscf-setup", workspace=_REPO,
        evidence=[{"kind": "test",
                   "ref": "tests/test_skill_registry.py::test_public_path_publish_with_evidence_builds_ssh_commands"}],
        runs=_stamped_green_run(),
        run_fn=_ok_run(calls), **_CONN)
    assert ok is True, msg
    assert msg.endswith("AGENT_SKILLS/casscf-setup.md")
    assert len(calls) == 2
    joined = " ".join(" ".join(map(str, c)) for c in calls)
    assert "ssh" in joined or "rsync" in joined
    assert "casscf-setup.md" in joined


def test_public_path_publish_chemistry_needs_calc_or_test():
    """Domain 'chemistry' (the dashboard's chemistry skills): a VALID
    recipe evidence alone still cannot attest it."""
    calls: list = []
    ok, msg = sr.publish_skill("casscf-setup", workspace=_REPO,
                               domain="chemistry",
                               evidence=[{"kind": "recipe",
                                          "ref": "workspace"}],
                               run_fn=_ok_run(calls), **_CONN)
    assert ok is False and "chemistry" in msg.lower()
    assert calls == []


def test_public_path_pull_yields_proposals_not_installs(tmp_path,
                                                        monkeypatch):
    """The dashboard's /skills pull call shape: downloaded skills must
    arrive as proposals, never in the active skill dir."""
    def run(cmd):
        target = Path(cmd[-1])
        (target / "delta.md").write_text("# Delta")
        return 0, "", ""

    proposed: list[str] = []

    def fake_propose(name, text):
        proposed.append(name)
        return "proposed", name

    monkeypatch.setattr(sr, "_propose_pulled_skill", fake_propose,
                        raising=False)
    ok, results = sr.pull_skills("", dest_dir=tmp_path, run_fn=run,
                                 **_CONN)
    assert ok is True
    assert proposed == ["delta"]
    assert all(status == "proposed" for _, status in results)
    # nothing installed into any active dir
    assert not (tmp_path / "delta.md").is_file()
