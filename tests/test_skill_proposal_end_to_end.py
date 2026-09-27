"""End-to-end: propose -> show -> accept -> discover -> skill tool.

The public call paths, not just the store: a proposal is created the way
skill_learning will create it, previewed the way the CLI previews it,
accepted the way a person accepts it, and then loaded through the very
executor the `skill` tool runs -- plus the REPL listing. A skill that
only works through the store but not through the caller was a real wave-4
failure; this is the gate against repeating it.
"""
import sys
import types

import pytest

from delfin.agent import skill_proposals as sp

BODY = ("# Draft CSV tables\n\n"
        "Use plot_energy_distribution with plot_type=bar_by_method\n"
        "when the folders span more than one functional.\n")


@pytest.fixture
def home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    mod = types.ModuleType("delfin.agent.skill_safety")
    mod.check = lambda text: []
    monkeypatch.setitem(sys.modules, "delfin.agent.skill_safety", mod)
    # ``from . import skill_safety`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "skill_safety",
                        mod, raising=False)
    ev = types.ModuleType("delfin.agent.evidence")
    ev.verify_evidence = lambda e, **kw: (True, "verified")
    monkeypatch.setitem(sys.modules, "delfin.agent.evidence", ev)
    # ``from . import evidence`` reads the package attribute
    # first; once the real module was imported by another test, a
    # sys.modules entry alone would be bypassed.
    monkeypatch.setattr(__import__("delfin.agent").agent, "evidence",
                        ev, raising=False)
    return tmp_path


def _run_cli(argv, monkeypatch):
    from delfin.agent import cli
    ns = cli.build_parser().parse_args(argv)

    class FakeOut:
        parts = []
        def write(self, s):
            self.parts.append(s)

    out = FakeOut()
    err = FakeOut()
    monkeypatch.setattr(sys, "stdout", out)
    monkeypatch.setattr(sys, "stderr", err)
    code = ns.func(ns)
    return code, "".join(out.parts), "".join(err.parts)


def test_full_path_propose_to_skill_tool(home, monkeypatch):
    # 1. propose, the way skill_learning (package 2) will call it
    p = sp.propose("draft-csv-tables", BODY,
                   evidence=[sp.Evidence(
                       kind="test", ref="tests/test_skill_proposals.py",
                       detail="store green on 6eb3d906")],
                   source="session end")
    assert p.status == "pending"

    # 2. invisible before the decision
    from delfin.agent.skills import discover_skills
    assert "draft-csv-tables" not in [s.name for s in discover_skills()]

    # 3. the human previews through the CLI and accepts
    code, shown, _ = _run_cli(["skills", "show", p.name], monkeypatch)
    assert code == 0 and "bar_by_method" in shown
    code, _, _ = _run_cli(["skills", "accept", p.name], monkeypatch)
    assert code == 0

    # 4. now it is a live skill
    found = {s.name: s for s in discover_skills()}
    assert "draft-csv-tables" in found
    assert found["draft-csv-tables"].body == BODY.strip()

    # 5. the `skill` tool executor finds and renders it
    from delfin.agent import api_client
    # _execute_skill lives on _DocToolExecutor -- the class behind the
    # `skill` tool dispatch at api_client.py:11701.
    client = api_client._DocToolExecutor.__new__(
        api_client._DocToolExecutor)
    perms = types.SimpleNamespace(workspace=home,
                                  skip_skill_discovery=False)
    result = client._execute_skill({"name": "draft-csv-tables"}, perms)
    import json
    payload = json.loads(result)
    assert payload.get("status") == "ok"
    assert payload["skill"] == "draft-csv-tables"
    assert "bar_by_method" in payload["content"]

    # 6. and a proposal that was never accepted still cannot be executed
    ghost = sp.propose("ghost-skill", "# Ghost\nnever accepted",
                       evidence=[sp.Evidence(kind="test", ref="t")],
                       source="test")
    result = client._execute_skill({"name": ghost.name}, perms)
    payload = json.loads(result)
    assert "not found" in payload.get("error", "")
