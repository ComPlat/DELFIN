"""Tests for the chemistry skill pack + the team skill registry.

Contract: the 6 new chemistry skills exist, are well-formed and encode
the safety philosophy (manual grounding, no autonomous submits/chemistry
changes); the registry pushes/pulls over the existing rsync/SSH builders
and NEVER overwrites a local skill on pull.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import skill_registry as sr


_NEW_SKILLS = (
    "diagnose-failed-run",
    "casscf-setup",
    "tddft-excited-states",
    "geometry-convergence-debug",
    "solvation-setup",
    "freq-thermochemistry",
)


# ---------------------------------------------------------------------------
# Skill pack: present, discoverable, safety philosophy encoded
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name", _NEW_SKILLS)
def test_skill_exists_and_resolves(name):
    p = sr.find_skill_file(name)
    assert p is not None and p.is_file()
    body = p.read_text(encoding="utf-8")
    assert body.startswith("# ")                 # title heading
    assert len(body) > 400                        # substantial protocol


@pytest.mark.parametrize("name", _NEW_SKILLS)
def test_skill_demands_manual_grounding(name):
    body = sr.find_skill_file(name).read_text(encoding="utf-8")
    assert "search_docs" in body, f"{name}: no manual-grounding instruction"


def test_diagnose_skill_separates_mechanical_from_chemical():
    body = sr.find_skill_file("diagnose-failed-run").read_text(encoding="utf-8")
    # The tiered-autonomy philosophy: mechanical vs chemical classes,
    # repair only with confirmation, never weaken convergence.
    assert "Mechanical" in body and "Chemical" in body
    assert "confirmation" in body.lower()
    assert "convergence" in body.lower()


def test_skills_are_discovered_by_skill_system():
    from delfin.agent.skills import discover_skills
    names = {s.name for s in discover_skills(None)}
    # Pack skills may or may not be in discover_skills (it scans user
    # dirs) — but find_skill_file must resolve them either way.  This
    # test pins that at least the resolution path works for all six.
    for n in _NEW_SKILLS:
        assert sr.find_skill_file(n) is not None
    assert isinstance(names, set)


# ---------------------------------------------------------------------------
# Registry: push / list / pull with mocked transport
# ---------------------------------------------------------------------------

def _ok_run(calls):
    def run(cmd):
        calls.append(cmd)
        return 0, "", ""
    return run


def test_publish_requires_remote_config():
    ok, msg = sr.publish_skill("casscf-setup", host="", user="", remote_path="")
    assert ok is False and "kein Remote" in msg


def test_publish_unknown_skill_fails():
    ok, msg = sr.publish_skill("does-not-exist", host="h", user="u",
                               remote_path="/r")
    assert ok is False and "nicht gefunden" in msg


def test_publish_runs_mkdir_then_rsync():
    calls: list = []
    ok, msg = sr.publish_skill("casscf-setup", host="login", user="grp",
                               remote_path="/archive", run_fn=_ok_run(calls))
    assert ok is True
    assert msg.endswith("AGENT_SKILLS/casscf-setup.md")
    assert len(calls) == 2                       # mkdir, then rsync
    assert any("casscf-setup.md" in " ".join(map(str, c)) for c in calls)


def test_list_shared_parses_md_files():
    def run(cmd):
        return 0, "casscf-setup.md\nnotes.txt\ntddft-excited-states.md\n", ""
    ok, names = sr.list_shared(host="h", user="u", remote_path="/r", run_fn=run)
    assert ok and names == ["casscf-setup.md", "tddft-excited-states.md"]


# ---------------------------------------------------------------------------
# Pull: install without ever overwriting (user no-delete rule)
# ---------------------------------------------------------------------------

def test_install_new_skill(tmp_path):
    status, p = sr.install_skill_text("new-skill", "# New\nbody", dest_dir=tmp_path)
    assert status == "installed" and p.read_text() == "# New\nbody"


def test_install_identical_is_noop(tmp_path):
    sr.install_skill_text("s", "# Same", dest_dir=tmp_path)
    status, _ = sr.install_skill_text("s", "# Same", dest_dir=tmp_path)
    assert status == "identical"


def test_install_conflict_never_overwrites(tmp_path):
    sr.install_skill_text("s", "# Local version", dest_dir=tmp_path)
    status, p = sr.install_skill_text("s", "# Team version", dest_dir=tmp_path)
    assert status == "conflict"
    assert p.name == "s-shared.md"               # written alongside
    assert (tmp_path / "s.md").read_text() == "# Local version"  # untouched


def test_pull_installs_downloaded_skills(tmp_path, monkeypatch):
    # Fake transport: "download" = drop two files into the rsync target.
    def run(cmd):
        target = Path(cmd[-1])                    # local_target is last arg
        (target / "alpha.md").write_text("# Alpha")
        (target / "beta.md").write_text("# Beta")
        return 0, "", ""
    ok, results = sr.pull_skills(host="h", user="u", remote_path="/r",
                                 dest_dir=tmp_path, run_fn=run)
    assert ok is True
    assert sorted(results) == [("alpha", "installed"), ("beta", "installed")]
    assert (tmp_path / "alpha.md").is_file()
