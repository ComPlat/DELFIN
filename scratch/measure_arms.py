"""Measure both arms against one set of cases, in each arm's own tree.

    PYTHONPATH=. python measure_arms.py <tree>

The module's real path is printed first, so the answer cannot come from an
installed copy of DELFIN instead of the tree under measurement.

Two test-file shapes matter and are measured separately, because an arm that
only handles one of them blocks or waves through every member:

  * EXISTING -- the test file is on the base and red there. This is the shape
    both arms' own fixtures use.
  * NEW -- the test file does not exist on the base. This is every wave
    member's normal case (DELIVERY_RULES 4: "a NEW test must be red on the
    baseline", which for a new file means it is not there yet).

Cases, and the verdict the package requires:

  A existing test, stage row + a probe log that touched the change  ACCEPT
  B existing test, stage row + a probe log that touched nothing     refuse
  C existing test, NO stage row, only a member-written probe log    refuse
  D existing test, a model name in the PR body the squash publishes refuse
  E NEW test file, stage row + a probe log that touched the change  ACCEPT
  F existing test, no decision table at all                         refuse
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
import tempfile
from pathlib import Path

PY = sys.executable
AGENT = "424242+delfin-agent@users.noreply.github.com"
PERSON = "person@example.org"

TABLE = ("| item | package | decision | test | real run | reason |\n"
         "|---|---|---|---|---|---|\n"
         "| x | g1 | landed | tests/test_new.py | probe1 | fine |\n")

TOUCHED = json.dumps({"type": "tool_use", "name": "edit_file",
                      "input": {"path": "tests/test_new.py"}})
CHATTED = json.dumps({"type": "text", "text": "I read the code and answered."})


def git(*a, cwd, env=None):
    r = subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True,
                       env={**os.environ, **(env or {})})
    assert r.returncode == 0, (a, r.stderr)
    return r.stdout


def make_repo(tmp: Path, existing_test: bool) -> Path:
    w = tmp / "wt"
    (w / "scratch").mkdir(parents=True)
    (w / "tests").mkdir()
    env = {"GIT_AUTHOR_NAME": "p", "GIT_AUTHOR_EMAIL": PERSON,
           "GIT_COMMITTER_NAME": "p", "GIT_COMMITTER_EMAIL": PERSON,
           "GIT_CONFIG_GLOBAL": str(tmp / "gitconfig"), "HOME": str(tmp)}
    git("init", "-q", "-b", "main", cwd=w, env=env)
    (w / "scratch" / ".gitignore").write_text("*\n")
    if existing_test:
        (w / "tests" / "test_new.py").write_text("def test_new():\n    assert False\n")
    (w / "base.txt").write_text("base\n")
    git("add", "-A", cwd=w, env=env)
    git("commit", "-qm", "base", cwd=w, env=env)
    git("update-ref", "refs/remotes/origin/main", "HEAD", cwd=w, env=env)
    (w / "tests" / "test_new.py").write_text("def test_new():\n    assert True\n")
    (w / "code.py").write_text("VALUE = 1\n")
    git("add", "-A", cwd=w, env=env)
    git("commit", "-qm", "the change and its test", cwd=w, env=env)
    return w


def build(L, wave, case, tmp):
    """The repository one case presents to the gate."""
    w = make_repo(tmp, existing_test=case["existing"])
    body = "# a title\n\n"
    if case["table"]:
        body += TABLE + "\n"
    body += case["extra"]
    if case.get("trailer"):
        body += ("\nCo-Authored-By: DELFIN-Agent "
                 f"<{AGENT}>\n")
    (w / "scratch" / "PR.md").write_text(body)
    if case["log"] is not None:
        d = w / "scratch" / "realrun"
        d.mkdir()
        (d / "probe1-branch.jsonl").write_text(case["log"] + "\n")
    if case["stage"]:
        (tmp / "w54t").mkdir(exist_ok=True)
        wave.root = lambda: tmp
        wave.record("w54t", "realrun", "g1", rows=[
            {"probe": "probe1", "side": "base", "stop": "end_turn", "tools": 2, "error": ""},
            {"probe": "probe1", "side": "branch", "stop": "end_turn", "tools": 4,
             "error": ""}])
    return w


CASES = [
    dict(name="A existing,row+touched", existing=True, table=True, stage=True,
         log=TOUCHED, extra="", want="ACCEPT", trailer=True),
    dict(name="B existing,row+chatted", existing=True, table=True, stage=True,
         log=CHATTED, extra="", want="refuse", trailer=True),
    dict(name="C existing,member log only", existing=True, table=True, stage=False,
         log=TOUCHED, extra="", want="refuse", trailer=True),
    dict(name="D existing,model in body", existing=True, table=True, stage=True,
         log=TOUCHED, extra="built with kit.glm-5.3-flash", want="refuse", trailer=True),
    dict(name="E NEW test,row+touched", existing=False, table=True, stage=True,
         log=TOUCHED, extra="", want="ACCEPT", trailer=True),
    dict(name="F existing,no table", existing=True, table=False, stage=True,
         log=TOUCHED, extra="just prose, no table", want="refuse", trailer=True),
]


def main(tree: str) -> int:
    sys.path.insert(0, os.path.abspath(tree))
    import delfin
    print(f"# module: {os.path.realpath(delfin.__file__)}")
    from delfin.agent.swarm import land as L
    from delfin.agent.swarm import wave
    L.identity_email = lambda: PERSON
    for n in ("coauthor_email", "_coauthor_email"):
        if hasattr(L, n):
            setattr(L, n, lambda: AGENT)
    L.PYCI = PY
    L.review_gate = lambda wt, run, tid, agent: (True, "stub", {})

    def gate(w):
        """Drive the PRODUCT entry point: land(), which delfin-agent swarm wave
        land calls (cli_swarm). check() alone never composes the squash message,
        so a probe that only calls check() cannot see that rule at all."""
        root = w.parent
        (root / "w54t").mkdir(exist_ok=True)
        (root / "w54t" / "state.json").write_text(json.dumps({
            "repo": str(w), "land": "pr",
            "agents": {"g1": {"worktree": str(w), "branch": "wave/w54t-g1",
                              "task": {"id": "g1", "base": "origin/main"}}}}))
        L.ROOT = root
        out = []
        recs = L.land(["w54t"], dry_run=True, log=lambda m: out.append(str(m)))
        if not recs:
            return False, "land() considered no member: " + " | ".join(out[-3:])
        r = recs[0]
        return r["outcome"] in ("ready", "opened"), f"{r['outcome']}: {r['why']}"

    wrong = 0
    for c in CASES:
        try:
            with tempfile.TemporaryDirectory() as t:
                tmp = Path(t)
                w = build(L, wave, c, tmp)
                ok, why = gate(w)
            got = "ACCEPT" if ok else "refuse"
            flag = "" if got == c["want"] else f"  <-- WRONG, wanted {c['want']}"
            wrong += bool(flag)
            print(f"{c['name']}\t{got}\t{str(why)[:130]}{flag}")
        except Exception as exc:
            print(f"{c['name']}\tERROR\t{type(exc).__name__}: {exc}")
            wrong += 1
    print(f"# wrong verdicts: {wrong}/{len(CASES)}")
    return wrong


if __name__ == "__main__":
    raise SystemExit(1 if main(sys.argv[1]) else 0)
