#!/usr/bin/env python3
"""Acceptance for generalist_long_py_feature (``tasklog/``).

Required: the suite is green; hidden tests of every layer pass
(models: default "normal", invalid value refused, to_dict/from_dict
with a legacy record lacking the key; store: a file written before the
change loads; query: open tasks ordered high, normal, low then by id,
render suffix " (high)"/" (low)", normal unchanged; cli: add
--priority, list --priority filter, an invalid choice exits 2, list
--all keeps id order); the README documents --priority; CHANGELOG
[Unreleased] has an entry naming the priority; a test under tests/
mentions priority.
Graded only: nothing outside tasklog/, tests/, README.md, CHANGELOG.md
changed; the original tests are unchanged; the README's file-format
example shows the priority key.
argv[1] is the workspace.
"""
from __future__ import annotations

import re
import sys

from generalist_check import (Report, base, blob_at, changed_paths, pytest,
                              read, tree_text, workspace_repo)

HIDDEN = '''import json

import pytest

from tasklog import query, store
from tasklog.cli import main
from tasklog.models import Task


def test_hidden_default_and_validation():
    assert Task(id=1, title="a").priority == "normal"
    assert Task(id=1, title="a", priority="high").priority == "high"
    with pytest.raises(ValueError):
        Task(id=1, title="a", priority="urgent")


def test_hidden_json_form_and_legacy_record():
    t = Task(id=2, title="b", priority="low")
    assert t.to_dict()["priority"] == "low"
    assert Task.from_dict(t.to_dict()) == t
    legacy = Task.from_dict({"id": 3, "title": "c", "done": True})
    assert legacy.priority == "normal" and legacy.done


def test_hidden_a_file_from_before_the_change_loads(tmp_path):
    p = tmp_path / "t.json"
    p.write_text(json.dumps({"version": 1, "tasks": [
        {"id": 1, "title": "old", "done": False}]}), encoding="utf-8")
    (t,) = store.load(p)
    assert t.priority == "normal"


def test_hidden_open_tasks_order_by_priority_then_id():
    ts = [Task(1, "a", priority="low"), Task(2, "b"),
          Task(3, "c", priority="high"), Task(4, "d", priority="high"),
          Task(5, "e", done=True, priority="high")]
    assert [t.id for t in query.open_tasks(ts)] == [3, 4, 2, 1]


def test_hidden_render_marks_only_non_normal():
    assert query.render(Task(1, "a")) == "[ ] 1 a"
    assert query.render(Task(2, "b", priority="high")) == "[ ] 2 b (high)"
    assert query.render(Task(3, "c", done=True, priority="low")) == \\
        "[x] 3 c (low)"


def test_hidden_cli_end_to_end(tmp_path, capsys):
    f = str(tmp_path / "t.json")
    main(["--file", f, "add", "A"])
    main(["--file", f, "add", "--priority", "high", "B"])
    main(["--file", f, "add", "--priority", "low", "C"])
    capsys.readouterr()
    main(["--file", f, "list"])
    assert capsys.readouterr().out == "[ ] 2 B (high)\\n[ ] 1 A\\n[ ] 3 C (low)\\n"
    main(["--file", f, "list", "--priority", "low"])
    assert capsys.readouterr().out == "[ ] 3 C (low)\\n"
    main(["--file", f, "list", "--all"])
    assert capsys.readouterr().out.splitlines()[0] == "[ ] 1 A"
    with pytest.raises(SystemExit) as exc:
        main(["--file", f, "add", "--priority", "urgent", "D"])
    assert exc.value.code == 2
'''

_UNRELEASED = re.compile(r"^## \[Unreleased\]\s*$(.*?)(?=^## \[|\Z)",
                         re.M | re.S)


def main() -> int:
    repo = workspace_repo("tasklog")
    r = Report("generalist_long_py_feature")
    if not r.check("repository present", (repo / ".git").is_dir()):
        return r.finish()
    since = base(repo)
    ok, tail = pytest(repo)
    r.check("whole suite green", ok, detail=tail, weight=2)
    ok, tail = pytest(repo, "tests/test_hidden_priority.py",
                      extra={"tests/test_hidden_priority.py": HIDDEN})
    r.check("hidden tests of every layer pass", ok, detail=tail, weight=4)
    readme = read(repo, "README.md")
    r.check("README documents --priority", "--priority" in readme)
    m = _UNRELEASED.search(read(repo, "CHANGELOG.md"))
    r.check("CHANGELOG [Unreleased] names the priority",
            bool(m) and "priorit" in m.group(1).lower())
    own = [k for k, v in tree_text(repo, (".py",)).items()
           if k.startswith("tests/") and "priority" in v]
    r.check("a test of the priority was added", bool(own),
            detail=", ".join(own))
    allowed = ("tasklog/", "tests/", "README.md", "CHANGELOG.md")
    stray = sorted(p for p in changed_paths(repo, since)
                   if not p.startswith(allowed))
    r.check("no unrelated changes", not stray, required=False,
            detail=", ".join(stray))
    old_tests = ("tests/test_models.py", "tests/test_store.py")
    edited = [p for p in old_tests
              if read(repo, p) != (blob_at(repo, since, p) or "\0")]
    r.check("existing model/store tests unchanged", not edited,
            required=False, detail=", ".join(edited))
    fmt = readme.split("## File format", 1)[-1]
    r.check("README file format shows the priority key",
            '"priority"' in fmt, required=False)
    return r.finish()


if __name__ == "__main__":
    sys.exit(main())
