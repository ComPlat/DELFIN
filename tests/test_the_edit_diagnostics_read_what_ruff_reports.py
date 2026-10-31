"""The diagnostics after an edit report what the pinned ruff reports.

Input: delfin/agent/diagnostics.py, driven with a stand-in ``ruff`` on PATH
that answers the way ruff 0.16.3 (the CI pin, .github/ci/constraints-ci.txt)
does:
- an unknown or removed rule in ``--select`` (E999 was removed) is exit 2,
  "ruff failed" on stderr, nothing on stdout;
- a syntax error is a finding with code ``invalid-syntax`` and severity
  ``error``;
- ``--stdin-filename P -`` checks the text on stdin as if it were P.

Measured before this change with the real ruff 0.16.3: a file with an
undefined name and an unused import gave no finding and no note, because
the selected E999 made every run exit 2 with an empty stdout, which was
read as an empty list.

The pre-edit text was checked from a copy written next to the file
(``.<name>.lsp1-old.py``); the harness wrote it with Path.write_text, which
follows a symlink, so a link planted under that name in a repository made
the harness overwrite its target -- outside the workspace, outside the
cage -- with the file's previous text, and then removed the link. The
pre-image now goes to ruff on stdin, and nothing is written.

pyright is no longer run: the harness started it on the workspace without
the cage, its findings of the old text were never filtered, and its start
alone costs seconds per edit.
"""
from __future__ import annotations

import os
import stat

import pytest

from delfin.agent import diagnostics as diag

FAKE_RUFF = r"""#!/usr/bin/env python3
import json, re, sys

args = sys.argv[1:]
select = args[args.index("--select") + 1].split(",") if "--select" in args else []
if "E999" in select:
    sys.stderr.write("ruff failed\n  Cause: Rule `E999` was removed and "
                     "cannot be selected.\n")
    sys.exit(2)
if "--stdin-filename" in args:
    content = sys.stdin.read()
else:
    with open(args[-1], encoding="utf-8") as fh:
        content = fh.read()
out = []
for i, line in enumerate(content.splitlines(), start=1):
    if line.rstrip().endswith("(:"):
        out.append({"code": "invalid-syntax", "severity": "error",
                    "message": "Expected a parameter",
                    "location": {"row": i, "column": 7}})
    for name in ("alpha", "beta"):
        if re.search(rf"\b{name}\(\)", line):
            out.append({"code": "F821", "message": f"Undefined name `{name}`",
                        "location": {"row": i, "column": line.index(name) + 1}})
print(json.dumps(out))
sys.exit(1 if out else 0)
"""


@pytest.fixture()
def ruff_016(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    exe = bin_dir / "ruff"
    exe.write_text(FAKE_RUFF, encoding="utf-8")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", str(bin_dir) + os.pathsep
                       + os.environ.get("PATH", ""))
    diag.reset_unavailability_note()
    ws = tmp_path / "ws"
    ws.mkdir()
    return ws


def test_a_new_undefined_name_is_reported(ruff_016):
    p = ruff_016 / "m.py"
    p.write_text("def go():\n    alpha()\n    beta()\n")
    found = diag.check(p, ruff_016, old_text="def go():\n    alpha()\n")
    assert [(d.code, d.line) for d in found] == [("F821", 3)]


def test_a_syntax_error_is_an_error(ruff_016):
    p = ruff_016 / "m.py"
    p.write_text("def f(:\n    pass\n")
    found = diag.check(p, ruff_016)
    assert [(d.code, d.severity) for d in found] == [("invalid-syntax",
                                                     "error")]


def test_a_failed_ruff_run_is_a_note_not_a_clean_file(ruff_016, monkeypatch):
    monkeypatch.setattr(diag, "_SELECT", "F,E999")
    p = ruff_016 / "m.py"
    p.write_text("def go():\n    beta()\n")
    result = diag.check_detailed(p, ruff_016)
    assert result.diagnostics == []
    assert result.note.startswith("diagnostics unavailable: ruff failed")


def test_the_pre_edit_text_is_never_written_beside_the_file(ruff_016,
                                                            tmp_path):
    victim = tmp_path / "outside.txt"
    victim.write_text("USER DATA\n")
    p = ruff_016 / "x.py"
    p.write_text("def go():\n    beta()\n")
    os.symlink(victim, ruff_016 / ".x.py.lsp1-old.py")
    found = diag.check(p, ruff_016, old_text="def go():\n    alpha()\n")
    assert victim.read_text() == "USER DATA\n"
    assert os.path.islink(ruff_016 / ".x.py.lsp1-old.py")
    assert [d.code for d in found] == ["F821"]
    assert sorted(f.name for f in ruff_016.iterdir()) == [
        ".x.py.lsp1-old.py", "x.py"]


def test_pyright_is_not_started(ruff_016, tmp_path, monkeypatch):
    marker = tmp_path / "pyright-ran"
    exe = tmp_path / "bin" / "pyright"
    exe.write_text(f"#!/bin/sh\ntouch '{marker}'\necho '{{}}'\n")
    exe.chmod(exe.stat().st_mode | stat.S_IXUSR)
    p = ruff_016 / "m.py"
    p.write_text("x = 1\n")
    diag.check(p, ruff_016)
    assert not marker.exists()
