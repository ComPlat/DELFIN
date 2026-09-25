"""tool_args on the execute path (api_client._malformed_arguments).

Wired LAST in the refusal chain, as _missing_required_argument is: a
role or mode refusal outranks "your call was malformed" (the reason is
written at the chain). Before any dialog, so a broken call never costs
the user a question. Aliases (`file_path` for `path`) keep working and
extra fields pass; an unknown field is named when it explains a
missing one.
"""

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions

@pytest.fixture
def ws(tmp_path):
    (tmp_path / "a.txt").write_text("one\ntwo\n", encoding="utf-8")
    perms = KitToolPermissions(workspace=str(tmp_path))
    perms.mode = "acceptEdits"
    perms.task_session_id = "tool-args-integration"
    return tmp_path, perms


def _call(ws_fixture, name, args):
    _, perms = ws_fixture
    return A._DocToolExecutor().execute(name, args, perms)


def _error_text(out: str) -> str:
    try:
        return json.loads(out).get("error", "")
    except (json.JSONDecodeError, TypeError):
        return ""


def test_a_misspelt_field_is_refused_with_a_hint(ws):
    out = _call(ws, "read_file", {"pth": "a.txt"})
    err = _error_text(out)
    assert "did you mean 'path'?" in err
    assert "two" not in out           # nothing was read (file: one/two)


def test_an_alias_still_works(ws):
    # _ARG_ALIASES is documented tolerance for weak models; kept.
    out = _call(ws, "read_file", {"file_path": "a.txt"})
    assert "one" in out


def test_wrongly_typed_edits_is_refused_with_shape(ws):
    out = _call(ws, "multi_edit", {"path": "a.txt", "edits": "not a list"})
    assert "array" in _error_text(out)


# NOTE: a quoted integer is deliberately NOT tested here. The
# executor's _as_int coerces "2" today, so "read two lines" would
# pass before AND after the wiring -- it measures nothing. What the
# wiring adds for repairs is the Result.repairs record, which has no
# channel to the caller until api_client builds one.


def test_a_broken_call_asks_nothing(ws):
    asked = []
    perms = ws[1]
    perms.mode = "default"            # a write would ask here
    perms.confirm_callback = lambda *a: (asked.append(a), True)[1]
    out = _call(ws, "write_file", {"path": "b.txt", "content": ["x"]})
    assert "string" in _error_text(out)
    assert not asked
    assert not (ws[0] / "b.txt").exists()
