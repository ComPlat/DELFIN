"""A field the model sent must reach the tool, or the call is refused.

The failure this guards is not an error: a tool that DROPS a field it did
not recognise runs on its defaults and answers with real data about a
different subject. ``api_client._missing_required_argument``'s own
docstring calls that the worst thing a scientific tool can hand back
("not a refusal, not an error: plausible data about the wrong subject").

Every test here drives the real entry point -- ``_doc_executor.execute``
on a real workspace with a real permissions object -- not a helper.

Measured on the base commit ``cdadd9efb`` before this file existed:

* ``grep_file(pattern="def", file_path="calc")`` answered with files from
  OUTSIDE ``calc`` (it searched the whole workspace): the alias
  ``file_path`` was translated for the schema check but never written back
  under ``path``, and ``grep_file`` reads ``arguments["path"]``.
* ``list_files(glob="*.py")`` listed every file in the workspace:
  ``glob`` is not a schema field, was silently passed over, and the call
  ran with ``pattern`` absent.
* ``grep_file(pattern="def", pth="calc")`` searched the whole workspace:
  the misspelt field was passed over, exactly as ``glob`` was.

The first two are fixed by ``patches/api_client.alias-writeback.patch``,
which this branch SHIPS UNAPPLIED: ``delfin/agent/api_client.py`` is part
of the agent's own safety layer and the Self-Modification Guard refuses an
in-session edit of it. Their tests below are therefore xfail (strict) with
that reason, and they turn to pass the moment the patch is applied --
``tests/test_the_alias_patch_is_the_documented_fix`` pins the pair.

The third is fixed here and now, in ``delfin/agent/tool_args.py``.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import api_client as A


@pytest.fixture
def ws(tmp_path):
    """A workspace where the SUBJECT of a search gives a different answer.

    ``calc/orca.out`` holds the pattern and lives in a subdirectory, so a
    search that reaches ``calc`` and one that misses it cannot be confused;
    ``notes.md`` is a file with no match at all.
    """
    (tmp_path / "calc").mkdir()
    (tmp_path / "sub").mkdir()
    (tmp_path / "hello.py").write_text("def main():\n    return 1\n",
                                       encoding="utf-8")
    (tmp_path / "notes.md").write_text("# Notes\n\ntext\n", encoding="utf-8")
    (tmp_path / "sub" / "deep.py").write_text("def other():\n    return 2\n",
                                              encoding="utf-8")
    (tmp_path / "calc" / "orca.out").write_text("def stub():\n",
                                                encoding="utf-8")
    perms = A.KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    return tmp_path, perms


def _call(ws_fixture, tool: str, arguments: dict) -> str:
    _root, perms = ws_fixture
    return str(A._doc_executor.execute(tool, dict(arguments), perms))


def _error(out: str) -> str:
    try:
        parsed = json.loads(out)
    except (ValueError, TypeError):
        return ""
    return str(parsed.get("error") or "") if isinstance(parsed, dict) else ""


# ---------------------------------------------------------------------------
# The alias must reach the field the schema names
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("alias", ["file_path", "filename"])
@pytest.mark.xfail(
    strict=True,
    reason="needs patches/api_client.alias-writeback.patch, which ships "
           "unapplied: delfin/agent/api_client.py is the agent's own "
           "safety layer and the Self-Modification Guard refuses the edit "
           "in-session. RED on base for this reason, and green the moment "
           "the patch is applied.")
def test_a_path_alias_reaches_the_tool_it_names(ws, alias):
    """Red on base: the whole workspace was searched, not ``calc``."""
    out = _call(ws, "grep_file", {"pattern": "def", alias: "calc"})
    assert not _error(out), _error(out)
    assert "calc/orca.out" in out, out
    assert "hello.py" not in out, (
        f"grep_file({alias}='calc') answered with a file outside calc: "
        f"{out}")


def test_the_documented_field_still_answers_the_same_way(ws):
    """Control: `path` was already correct, and must stay correct."""
    out = _call(ws, "grep_file", {"pattern": "def", "path": "calc"})
    assert not _error(out), _error(out)
    assert "calc/orca.out" in out, out
    assert "hello.py" not in out, out


# ---------------------------------------------------------------------------
# A synonym the table names must be accepted as that field
# ---------------------------------------------------------------------------

@pytest.mark.xfail(
    strict=True,
    reason="needs patches/api_client.alias-writeback.patch (see the "
           "module docstring): the alias table lives in "
           "delfin/agent/api_client.py, which the Self-Modification Guard "
           "protects. RED on base for this reason.")
def test_glob_is_accepted_as_the_pattern_it_names(ws):
    """Red on base: ``glob`` was passed over and ALL files were listed."""
    out = _call(ws, "list_files", {"glob": "*.py"})
    assert not _error(out), _error(out)
    assert "hello.py" in out, out
    assert "notes.md" not in out, (
        f"list_files(glob='*.py') listed a file that does not match: {out}")
    assert "orca.out" not in out, out


def test_the_documented_pattern_still_answers_the_same_way(ws):
    """Control: `pattern` was already correct."""
    out = _call(ws, "list_files", {"pattern": "*.py"})
    assert "hello.py" in out and "notes.md" not in out, out


# ---------------------------------------------------------------------------
# A stray field that stands in for an absent property is REFUSED
# ---------------------------------------------------------------------------

def test_a_typoed_field_is_refused_and_named(ws):
    """Red on base: `pth` was passed over, so the whole workspace was
    searched and real content came back about the wrong subject."""
    out = _call(ws, "grep_file", {"pattern": "def", "pth": "calc"})
    err = _error(out)
    assert err, f"a stray field ran anyway and answered: {out}"
    assert "pth" in err and "path" in err, err
    # The refusal must say what to do next, not only what was wrong.
    assert "Re-send" in err or "re-send" in err, err


def test_the_refused_call_never_reaches_the_tool(ws):
    """The refusal is the point: no answer about the wrong subject."""
    out = _call(ws, "grep_file", {"pattern": "def", "pth": "calc"})
    assert "orca.out" not in out and "hello.py" not in out, (
        f"a refused call still returned matches: {out}")


def test_a_stray_field_that_names_nothing_does_not_refuse(ws):
    """The narrowness has to be pinned, or the rule drifts into "refuse
    every extra field" the first time someone widens it.

    `zzzqqq` resembles no schema field, so there is no name to offer and
    nothing to confuse: the call is complete without it.
    """
    out = _call(ws, "read_file", {"path": "hello.py", "zzzqqq": 1})
    assert not _error(out), _error(out)
    assert "def main" in out, out


# ---------------------------------------------------------------------------
# The patch this branch could not apply is a tracked fact, not a hole
# ---------------------------------------------------------------------------

def test_the_alias_patch_is_shipped_and_still_applies():
    """delfin/agent/api_client.py is the agent's own safety layer, so the
    alias-writeback fix could not be applied in-session. It ships as a
    unified diff against the base; this pins that it is present, that its
    hunks are the two the docstring names, and that it still applies --
    a patch that no longer applies is a fix that quietly stopped existing.
    """
    import subprocess
    from pathlib import Path

    patch = (Path(__file__).resolve().parents[1] / "patches"
             / "api_client.alias-writeback.patch")
    assert patch.is_file(), f"the promised patch is missing: {patch}"
    text = patch.read_text(encoding="utf-8")
    assert "+++ b/delfin/agent/api_client.py" in text, (
        "the patch must target the real module, not the scratch copy")
    assert "/home/" not in text and "/tmp/" not in text, (
        "a patch carrying an absolute home or temp path cannot apply "
        "anywhere else")
    # Both halves of the fix must be in it: the synonym entries and the
    # write-back into the caller's dict.
    assert '"pattern": ("glob", "regex")' in text, text
    assert "arguments[canonical] = arguments[alias]" in text, text
    repo = Path(__file__).resolve().parents[1]
    got = subprocess.run(
        ["git", "apply", "--check", str(patch)],
        cwd=repo, capture_output=True, text=True)
    # A tree that has ALREADY had the patch applied reports "does not
    # apply"; that is success from this test's point of view, so it is
    # told apart from a patch that is genuinely broken.
    if got.returncode != 0:
        assert "already exists" in got.stderr or "patch does not apply" in \
            got.stderr, got.stderr
        assert '"pattern": ("glob", "regex")' in (
            repo / "delfin" / "agent" / "api_client.py").read_text(
                encoding="utf-8"), (
            "the patch does not apply AND the module does not carry the "
            f"fix -- the fix is gone. git said: {got.stderr.strip()}")


# ---------------------------------------------------------------------------
# The tolerance the codebase decided on must survive
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("extra", [{"chunks": 4}, {"note": "why"},
                                   {"cache": True}])
def test_a_harmless_extra_field_still_passes(ws, extra):
    """Control, green on base and after: an extra field that names nothing
    in the schema is noise beside a complete call, and refusing it would
    be a regression (api_client.py:2433 says so in as many words)."""
    args = {"path": "hello.py", **extra}
    out = _call(ws, "read_file", args)
    assert not _error(out), _error(out)
    assert "def main" in out, out


def test_an_extra_field_beside_the_field_it_resembles_still_passes(ws):
    """Control: the stray is redundant -- `path` is already there.

    Refusing this would be the defect in reverse: the call says everything
    it needs to, and the model only repeated itself under another name.
    """
    out = _call(ws, "read_file", {"path": "hello.py", "filepath": "hello.py"})
    assert not _error(out), _error(out)
    assert "def main" in out, out
