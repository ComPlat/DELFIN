"""The dashboard's approval row can answer a shell command with a change.

Input: an open shell approval in the dashboard (``KitConfirmBroker``), the
row's "Change…" button, the edited command and "Run changed command".
Output and semantics, the same as in the terminal and the window
(``terminal_confirm``, one grant store for both brokers):

- the request is refused as asked; the gate tells the model the command
  was changed and gives it word for word;
- the changed command runs once, exactly as written, for the same tool:
  the model's next call with it is not asked; a second one, a different
  one, or one behind a guard or widening question is asked as usual;
- an empty, unchanged or control-character change is refused, the row
  says why, and the request stays open;
- a command holding a secret is offered from an empty field.
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import terminal_confirm as tc

widgets = pytest.importorskip("ipywidgets")


@pytest.fixture
def broker(tmp_path, monkeypatch):
    from delfin.agent import attention, audit_log
    from delfin.agent.kit_confirm import KitConfirmBroker
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    monkeypatch.setattr(attention, "emit_attention", lambda *a, **kw: "")
    records = []
    monkeypatch.setattr(audit_log, "append", lambda rec, *a, **kw: records.append(rec))
    b = KitConfirmBroker(default_timeout_s=10.0)
    b.records = records
    return b


def _ask(broker, command, tool="bash", preview=None):
    out = {}

    def run():
        out["ok"] = broker.callback(tool, {"command": command},
                                    f"$ {command}" if preview is None else preview)
        out["changed_to"] = broker.last_changed_to

    t = threading.Thread(target=run, daemon=True)
    t.start()
    for _ in range(300):
        if broker._pending or "ok" in out:
            break
        time.sleep(0.01)
    return t, out


def _buttons(box):
    found = []
    for child in getattr(box, "children", ()):
        if isinstance(child, widgets.Button):
            found.append(child)
        found.extend(_buttons(child))
    return found


def _find(box, kind):
    for child in getattr(box, "children", ()):
        if isinstance(child, kind):
            return child
        inner = _find(child, kind)
        if inner is not None:
            return inner
    return None


def _change_in_row(broker, text):
    req = list(broker._pending)[0]
    row = broker._build_request_row(req, widgets)
    by_name = {b.description: b for b in _buttons(row)}
    by_name["Change…"].click()
    field = _find(row, widgets.Textarea)
    prefilled = field.value
    field.value = text
    by_name["Run changed command"].click()
    return req, row, prefilled


def test_the_row_changes_the_command_and_it_runs_once(broker):
    t, out = _ask(broker, "rm -r results")
    req, _row, prefilled = _change_in_row(broker, "rm -r results/old")
    t.join(5)
    assert prefilled == "rm -r results"
    assert out["ok"] is False and out["changed_to"] == "rm -r results/old"
    # The changed command: not asked, and recorded.
    assert broker.callback("bash", {"command": "rm -r results/old"},
                           "$ rm -r results/old") is True
    assert not broker._pending
    assert any(r.get("by") == "changed" for r in broker.records)
    # Once.
    t, out = _ask(broker, "rm -r results/old")
    assert broker._pending, "the same command a second time is asked"
    broker.resolve(list(broker._pending)[0], False)
    t.join(5)


@pytest.mark.parametrize("command, preview", [
    ("rm -r results/old -f", None),
    ("rm -r results/old", f"{A._WIDENING_BANNER}\n$ rm -r results/old"),
    ("rm -r results/old", f"{tc._GUARD_MARKER}\n$ rm -r results/old"),
])
def test_anything_but_the_plain_same_call_is_asked(broker, command, preview):
    t, _ = _ask(broker, "rm -r results")
    _change_in_row(broker, "rm -r results/old")
    t.join(5)
    t, out = _ask(broker, command, preview=preview)
    assert broker._pending, "asked as usual"
    broker.resolve(list(broker._pending)[0], False)
    t.join(5)


@pytest.mark.parametrize("text", ["", "rm -r results", "ls\x1b[2J"])
def test_a_change_that_is_no_command_leaves_the_request_open(broker, text):
    t, out = _ask(broker, "rm -r results")
    req, row, _ = _change_in_row(broker, text)
    assert req in broker._pending and req.decision is None
    status = [w.value for w in row.children if isinstance(w, widgets.HTML)]
    assert any("not run" in v for v in status), status
    broker.resolve(req, False)
    t.join(5)
    assert out["changed_to"] == ""


def test_a_command_with_a_secret_starts_from_an_empty_field(broker):
    key = "sk-" + "a1B2c3D4e5F6g7H8i9J0k1L2m3N4o5P6"
    t, _ = _ask(broker, f"curl -H 'Authorization: Bearer {key}' x")
    req = list(broker._pending)[0]
    row = broker._build_request_row(req, widgets)
    {b.description: b for b in _buttons(row)}["Change…"].click()
    assert _find(row, widgets.Textarea).value == ""
    broker.resolve(req, False)
    t.join(5)


def test_a_write_has_no_change_button(broker):
    from delfin.agent.kit_confirm import _ConfirmRequest
    req = _ConfirmRequest(seq=1, tool_name="write_file", args={"path": "a"}, preview="+x")
    row = broker._build_request_row(req, widgets)
    assert "Change…" not in {b.description for b in _buttons(row)}


def test_the_model_reads_the_change(broker, tmp_path):
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=broker.callback)

    def change_when_asked():
        for _ in range(300):
            if broker._pending:
                break
            time.sleep(0.01)
        _change_in_row(broker, "printf changed > out.txt")

    t = threading.Thread(target=change_when_asked, daemon=True)
    t.start()
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": "some-unlisted-command -x"}, perms)
    t.join(5)
    assert err and "changed the command" in err and "printf changed > out.txt" in err
    assert A._doc_executor._run_permission_gate(
        "bash", {"command": "printf changed > out.txt"}, perms) is None
