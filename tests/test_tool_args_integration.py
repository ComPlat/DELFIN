"""The integration the operator will wire: bad arguments are refused
by tool_args BEFORE the permission gate answers.

Strict xfail on purpose: today these fail because nothing calls
tool_args.check on the execute path (see the characterization in
test_tool_args_current_behavior.py). When the check is wired in, the
tests start passing -- strict turns that into an XPASS failure, which
is the built-in reminder to remove the marker in the same commit that
does the wiring.
"""

import json

import pytest

from delfin.agent import api_client as A
from delfin.agent.api_client import KitToolPermissions

pytestmark = pytest.mark.xfail(
    strict=True,
    reason="tool_args.check is not yet called on the execute path; "
    "remove this marker in the commit that wires it in",
)


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


def test_wrong_field_name_is_refused_with_hint(ws):
    out = _call(ws, "read_file", {"file_path": "a.txt"})
    err = _error_text(out)
    assert "file_path" in err
    assert "did you mean 'path'?" in err
    # and the file was NOT read: the alias table did not run
    assert "one" not in out


def test_wrongly_typed_edits_is_refused_with_shape(ws):
    out = _call(ws, "multi_edit", {"path": "a.txt", "edits": "not a list"})
    assert "array" in _error_text(out)


# NOTE: a quoted integer is deliberately NOT tested here. The
# executor's _as_int coerces "2" today, so "read two lines" would
# pass before AND after the wiring -- it measures nothing. What the
# wiring adds for repairs is the Result.repairs record, which has no
# channel to the caller until api_client builds one.


def test_broken_call_never_reaches_the_permission_gate(ws):
    # A broken call must be refused by the argument check BEFORE any
    # permission logic could ask a question -- the pre_tool_hook is the
    # earliest gate signal, so it must not see the broken call.
    seen = []

    def hook(name, args):
        seen.append(name)

    perms = ws[1]
    perms.pre_tool_hook = hook
    _call(ws, "read_file", {"file_path": "a.txt"})
    assert "read_file" not in seen
