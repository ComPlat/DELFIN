#!/usr/bin/env python3
"""Setup for generalist_long_py_feature: one feature across five layers.

``tasklog/`` -- a to-do keeper in four modules plus docs:
models.py (Task dataclass), store.py (JSON file), query.py (filter,
sort), cli.py (argparse: add / list / done), README.md (usage),
CHANGELOG.md, tests per module. The task adds a priority end to end:
the field with validation, persistence with backward compatibility
(files written before the change have no "priority" key), a CLI flag
on add and a filter on list, the list marker, tests and docs.
Precondition: the suite is green and nothing mentions priority yet.
argv[1] is the workspace.
"""
from __future__ import annotations

from generalist_repo import (GITIGNORE_PY, commit, expect, main, mark_base,
                             new_repo, pytest_rc, write)

MODELS = '''"""The task record."""

from dataclasses import dataclass


@dataclass
class Task:
    """One to-do item."""

    id: int
    title: str
    done: bool = False

    def __post_init__(self):
        if not self.title.strip():
            raise ValueError("a task needs a title")

    def to_dict(self) -> dict:
        """The JSON form."""
        return {"id": self.id, "title": self.title, "done": self.done}

    @classmethod
    def from_dict(cls, data: dict) -> "Task":
        """A Task from its JSON form."""
        return cls(id=int(data["id"]), title=str(data["title"]),
                   done=bool(data.get("done", False)))
'''

STORE = '''"""Tasks in one JSON file."""

import json
from pathlib import Path

from tasklog.models import Task


def load(path: Path) -> list:
    """Every task in ``path``; an absent file holds none."""
    path = Path(path)
    if not path.exists():
        return []
    data = json.loads(path.read_text(encoding="utf-8"))
    return [Task.from_dict(d) for d in data["tasks"]]


def save(path: Path, tasks: list) -> None:
    """Write ``tasks`` to ``path``."""
    payload = {"version": 1, "tasks": [t.to_dict() for t in tasks]}
    Path(path).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def next_id(tasks: list) -> int:
    """One more than the largest id."""
    return max((t.id for t in tasks), default=0) + 1
'''

QUERY = '''"""Selecting and ordering tasks."""


def open_tasks(tasks: list) -> list:
    """Tasks not done yet, in id order."""
    return sorted((t for t in tasks if not t.done), key=lambda t: t.id)


def render(task) -> str:
    """One line per task: ``[x] 3 title``."""
    mark = "x" if task.done else " "
    return f"[{mark}] {task.id} {task.title}"
'''

CLI = '''"""Command line: tasklog add|list|done."""

import argparse
import sys
from pathlib import Path

from tasklog import query, store
from tasklog.models import Task


def build_parser() -> argparse.ArgumentParser:
    """The argument parser."""
    p = argparse.ArgumentParser(prog="tasklog")
    p.add_argument("--file", default="tasks.json", type=Path)
    sub = p.add_subparsers(dest="command", required=True)
    add = sub.add_parser("add", help="add a task")
    add.add_argument("title")
    lst = sub.add_parser("list", help="list tasks")
    lst.add_argument("--all", action="store_true", help="include done tasks")
    done = sub.add_parser("done", help="mark a task done")
    done.add_argument("id", type=int)
    return p


def main(argv=None) -> int:
    """Run one command; returns the exit status."""
    args = build_parser().parse_args(argv)
    tasks = store.load(args.file)
    if args.command == "add":
        task = Task(id=store.next_id(tasks), title=args.title)
        tasks.append(task)
        store.save(args.file, tasks)
        print(f"added {task.id}")
    elif args.command == "list":
        shown = sorted(tasks, key=lambda t: t.id) if args.all \\
            else query.open_tasks(tasks)
        for t in shown:
            print(query.render(t))
    elif args.command == "done":
        for t in tasks:
            if t.id == args.id:
                t.done = True
        store.save(args.file, tasks)
    return 0


if __name__ == "__main__":
    sys.exit(main())
'''

README = """# tasklog

A small to-do keeper that stores tasks in a JSON file.

## Usage

```sh
python -m tasklog.cli add "Write the report"
python -m tasklog.cli list          # open tasks
python -m tasklog.cli list --all    # including done ones
python -m tasklog.cli done 1
```

`--file PATH` selects the task file (default `tasks.json`).

## File format

```json
{"version": 1, "tasks": [{"id": 1, "title": "Write the report", "done": false}]}
```
"""

CHANGELOG = """# Changelog

## [Unreleased]

## [1.0.0] - 2026-09-01

- First release: add, list, done.
"""

TESTS = {
    "tests/test_models.py": (
        "import pytest\n\nfrom tasklog.models import Task\n\n\n"
        "def test_round_trip():\n"
        "    t = Task(id=1, title='a')\n"
        "    assert Task.from_dict(t.to_dict()) == t\n\n\n"
        "def test_empty_title_is_refused():\n"
        "    with pytest.raises(ValueError):\n"
        "        Task(id=1, title=' ')\n"),
    "tests/test_store.py": (
        "from tasklog import store\nfrom tasklog.models import Task\n\n\n"
        "def test_save_and_load(tmp_path):\n"
        "    p = tmp_path / 't.json'\n"
        "    store.save(p, [Task(id=1, title='a')])\n"
        "    assert store.load(p) == [Task(id=1, title='a')]\n\n\n"
        "def test_absent_file_is_empty(tmp_path):\n"
        "    assert store.load(tmp_path / 'none.json') == []\n"),
    "tests/test_cli.py": (
        "from tasklog.cli import main\n\n\n"
        "def test_add_then_list(tmp_path, capsys):\n"
        "    f = str(tmp_path / 't.json')\n"
        "    main(['--file', f, 'add', 'Write'])\n"
        "    capsys.readouterr()\n"
        "    main(['--file', f, 'list'])\n"
        "    assert capsys.readouterr().out == '[ ] 1 Write\\n'\n"),
}


def build(workspace: str):
    repo = new_repo(workspace, "tasklog")
    files = {
        ".gitignore": GITIGNORE_PY,
        "pyproject.toml": ("[project]\nname = \"tasklog\"\nversion = \"1.0.0\"\n\n"
                           "[tool.pytest.ini_options]\npythonpath = [\".\"]\n"
                           "testpaths = [\"tests\"]\n"),
        "README.md": README,
        "CHANGELOG.md": CHANGELOG,
        "tasklog/__init__.py": "\"\"\"A small to-do keeper.\"\"\"\n",
        "tasklog/models.py": MODELS,
        "tasklog/store.py": STORE,
        "tasklog/query.py": QUERY,
        "tasklog/cli.py": CLI,
        **TESTS,
    }
    write(repo, files)
    mark_base(repo, commit(repo, "tasklog 1.0.0"))
    rc, out = pytest_rc(repo)
    expect("the given suite is green", rc == 0, out)
    expect("priority is new", not any("priorit" in v.lower()
                                      for v in files.values()))
    return repo


if __name__ == "__main__":
    main(build)
