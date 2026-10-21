"""A file git reads as configuration is asked about before it is written.

git runs programs named by files in the work tree: a ``filter=`` or
``diff=`` attribute in ``.gitattributes`` selects a driver, a ``.git``
gitfile points at a directory whose ``config`` names filters, hooks and
fsmonitor, and a ``HEAD`` next to ``objects/`` and ``refs/`` makes a
directory a bare repository with its own ``config``. ``git status`` and
``git add`` are auto-allowed, so a program named that way ran without a
dialog once its file existed. ``.git/**`` is already refused; these are
the files that reach the same configuration from outside it.

The tests assert on the gate's decision only; nothing here runs git.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


class _Asker:
    def __init__(self, answer=False):
        self.answer = answer
        self.asked = []
        self.last_refusal_reason = ""
        self.last_timed_out = False

    def callback(self, name, args, preview):
        self.asked.append(preview)
        return self.answer


def _perms(tmp_path, mode="acceptEdits", asker=None):
    return A.KitToolPermissions(
        workspace=str(tmp_path), mode=mode,
        confirm_callback=asker.callback if asker else None)


def _write(perms, path):
    return A._doc_executor._run_permission_gate(
        "write_file", {"path": path, "content": "x\n"}, perms)


GIT_CONFIG_FILES = [
    ".gitattributes",
    "sub/.gitattributes",
    ".gitmodules",
    ".git",
    "sub/.git",
    "HEAD",
    "sub/HEAD",
    ".GITATTRIBUTES",
]


@pytest.mark.parametrize("path", GIT_CONFIG_FILES)
def test_acceptedits_asks_before_writing_it(tmp_path, path):
    asker = _Asker(answer=True)
    assert _write(_perms(tmp_path, asker=asker), path) is None
    assert len(asker.asked) == 1, path
    assert "git" in asker.asked[0].lower()


@pytest.mark.parametrize("path", GIT_CONFIG_FILES)
def test_with_nobody_to_ask_it_is_refused(tmp_path, path):
    err = _write(_perms(tmp_path, mode="bypassPermissions"), path)
    assert err and "refused" in err, path


@pytest.mark.parametrize("cmd", [
    "printf '* filter=x\\n' > .gitattributes",
    "echo 'gitdir: g' > .git",
    "cp attrs sub/.gitattributes",
])
def test_the_shell_route_is_asked_too(tmp_path, cmd):
    (tmp_path / "sub").mkdir()
    (tmp_path / "attrs").write_text("x\n")
    asker = _Asker(answer=False)
    perms = _perms(tmp_path, asker=asker)
    out = A._doc_executor._gate_bash_write_targets(cmd, {}, perms)
    assert len(asker.asked) == 1, cmd
    assert out and "error" in json.loads(out)


@pytest.mark.parametrize("path", [
    ".gitignore", "sub/.gitignore", "README.md", "config.yaml",
    "refs/paper.bib", "objects/mol.xyz", "HEADER.txt",
])
def test_ordinary_files_are_not_asked_about(tmp_path, path):
    asker = _Asker(answer=True)
    assert _write(_perms(tmp_path, asker=asker), path) is None
    assert asker.asked == [], path
