"""Does the gate accept scratch/PR.md in the format the wave itself mandates?

DELIVERY_RULES 17 says: "First line: the title. Then one table row per item
(item, decision, commit, reason in one line)". The gate's proof check reads a
decision-table column called "package". If it accepts only that word, then
every branch written to the wave's own documented rule is refused for having
no proof table -- a false refusal on the house format.

    PYTHONPATH=. python scratch/probe_house_format.py <tree>
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
TOUCHED = json.dumps({"type": "tool_use", "name": "edit_file",
                      "input": {"path": "tests/test_new.py"}})

#: Rule 17's own four columns, plus the two the package asks a package row to
#: carry. The header word is "item", which is the wave's word for it.
TABLE_ITEM = ("| item | decision | test | real run | commit | reason |\n"
              "|---|---|---|---|---|---|\n"
              "| g1 | landed | tests/test_new.py | probe1 | abc1234 | measured |\n")
TABLE_PACKAGE = ("| package | test | real run | decision |\n"
                 "|---|---|---|---|\n"
                 "| g1 | `tests/test_new.py` | probe1 | landed |\n")


def git(*a, cwd, env=None):
    r = subprocess.run(["git", *a], cwd=cwd, capture_output=True, text=True,
                       env={**os.environ, **(env or {})})
    assert r.returncode == 0, (a, r.stderr)
    return r.stdout


def repo(tmp: Path) -> Path:
    w = tmp / "wt"
    (w / "scratch").mkdir(parents=True)
    (w / "tests").mkdir()
    env = {"GIT_AUTHOR_NAME": "p", "GIT_AUTHOR_EMAIL": PERSON,
           "GIT_COMMITTER_NAME": "p", "GIT_COMMITTER_EMAIL": PERSON,
           "GIT_CONFIG_GLOBAL": str(tmp / "gitconfig"), "HOME": str(tmp)}
    git("init", "-q", "-b", "main", cwd=w, env=env)
    (w / "scratch" / ".gitignore").write_text("*\n")
    (w / "tests" / "test_new.py").write_text("def test_new():\n    assert False\n")
    (w / "base.txt").write_text("x\n")
    git("add", "-A", cwd=w, env=env)
    git("commit", "-qm", "base", cwd=w, env=env)
    git("update-ref", "refs/remotes/origin/main", "HEAD", cwd=w, env=env)
    (w / "tests" / "test_new.py").write_text("def test_new():\n    assert True\n")
    git("add", "-A", cwd=w, env=env)
    git("commit", "-qm", "change", cwd=w, env=env)
    return w


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
    bad = 0
    for label, table in (("item (rule 17's word)", TABLE_ITEM),
                         ("package", TABLE_PACKAGE)):
        with tempfile.TemporaryDirectory() as t:
            tmp = Path(t)
            w = repo(tmp)
            (w / "scratch" / "PR.md").write_text("# a title\n\n" + table)
            (tmp / "w54t").mkdir()
            wave.root = lambda: tmp
            wave.record("w54t", "realrun", "g1", rows=[
                {"probe": "probe1", "side": "base", "stop": "end_turn", "tools": 2, "error": ""},
                {"probe": "probe1", "side": "branch", "stop": "end_turn", "tools": 4,
                 "error": ""}])
            d = w / "scratch" / "realrun"
            d.mkdir()
            (d / "probe1-branch.jsonl").write_text(TOUCHED + "\n")
            try:
                ok, why = L.check(w, "origin/main", run="w54t", tid="g1")
            except TypeError:
                ok, why = L.check(w, "origin/main")
            print(f"{label}: {'ACCEPTED' if ok else 'refused'} -- {str(why)[:150]}")
            if not ok:
                bad += 1
    return bad


if __name__ == "__main__":
    raise SystemExit(1 if main(sys.argv[1]) else 0)
