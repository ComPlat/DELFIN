"""How a tool call with bad arguments ends.

Written (runde3-s6) as the "before" of wiring ``tool_args.check`` into
the executor, and kept as its contract: the tolerances that stayed
(the `file_path` alias, extra fields, the quoted integer) still pass;
only the shape check changed (see the multi_edit case at the end).

A real permissions object on a tmp workspace is used so the file tools
actually run and their argument handling is what answers, not the
sandbox refusal. The gated tools (bash) still go through their gate
with these arguments, which is itself part of the "before".
"""

import json
from pathlib import Path

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions


@pytest.fixture
def ws(tmp_path: Path):
    (tmp_path / "a.txt").write_text("one\ntwo\nthree\n", encoding="utf-8")
    perms = KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "tool-args-characterization"
    return tmp_path, perms


def _call(ws_fixture, name, args):
    tmp_path, perms = ws_fixture
    # execute() returns a plain string -- file content, prose, or JSON
    # with an "error" key; never a promise about which.
    return A._DocToolExecutor().execute(name, args, perms)


def _is_error(out: str) -> bool:
    try:
        return "error" in json.loads(out)
    except (json.JSONDecodeError, TypeError):
        return False


def _error_text(out: str) -> str:
    return json.loads(out).get("error", "")


def test_wrong_field_name_is_silently_accepted_today(ws):
    # `file_path` is not in read_file's schema. Today the alias table
    # resolves it and the read succeeds: no schema check exists.
    out = _call(ws, "read_file", {"file_path": "a.txt"})
    assert not _is_error(out), out
    assert "one" in out


def test_missing_required_field_gives_plain_guard_message(ws):
    out = _call(ws, "read_file", {})
    assert _is_error(out)
    assert "path" in _error_text(out)
    # The current guard message explains the call, not the schema:
    assert "required" in _error_text(out)


def test_string_where_schema_says_integer_is_tolerated(ws):
    # `limit` is typed integer in the schema; a string is coerced
    # silently by _as_int.
    out = _call(ws, "read_file", {"path": "a.txt", "limit": "2"})
    assert not _is_error(out), out


def test_unknown_extra_field_is_silently_ignored(ws):
    out = _call(ws, "read_file", {"path": "a.txt", "chunks": 4})
    assert not _is_error(out), out


def test_bad_regex_pattern_gives_engine_not_schema_error(ws):
    # A schema-valid but broken pattern fails inside the regex engine,
    # not with an argument-schema message. The error does not even
    # name the offending field.
    out = _call(ws, "grep_file", {"pattern": "["})
    assert _is_error(out)
    assert "Invalid regex" in _error_text(out)
    assert "schema" not in _error_text(out).lower()


def test_multi_edit_without_edits_key_is_caught_by_the_guard(ws):
    # `edits` is required; the existing missing-argument guard names it
    # (message is about the call, phrased as prose, not a schema report
    # -- no types, no field list). But the guard only checks required
    # KEYS: a broken edit shape inside `edits` is not its business.
    out = _call(ws, "multi_edit", {"path": "a.txt"})
    assert _is_error(out)
    assert "edits" in _error_text(out)


def test_multi_edit_with_wrongly_typed_edits_names_the_shape(ws):
    # Was: the failure happened later, inside the executor, in terms that
    # did not say the expected shape. Since tool_args is wired in, the
    # schema answers first.
    out = _call(ws, "multi_edit", {"path": "a.txt", "edits": "not a list"})
    assert _is_error(out), out
    assert "array" in _error_text(out).lower()
