"""An approval of a shell command can be answered with a changed command.

Input: an open shell approval and the user's own version of the command
(key ``c``; ``permission.reply`` with ``choice: "c"`` and ``command``).
Output and semantics:

- the command as asked is refused, and the model reads that the user
  changed it, and the changed command, word for word;
- the changed command is approved once, exactly as written, for the same
  tool: the model's next call with it is not asked again. A second call,
  a call that differs by one character, another tool, or a request that
  carries a guard, calculation, outside-read or widening marker is asked
  as usual;
- the change is a pre-answer to the question only. Every check that runs
  before a question (deny rules, sandbox, isolation) still runs on the
  changed command;
- an empty, unchanged or control-character command, a command with any
  other key, and typed ``text`` are refused and change nothing.
"""

from __future__ import annotations

import threading
import time

import pytest

from delfin.agent import api_client as A
from delfin.agent import terminal_confirm as tc
from delfin.agent.service.approvals import InvalidAnswer, ProtocolBroker


@pytest.fixture
def broker(tmp_path, monkeypatch):
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    from delfin.agent import audit_log
    records = []
    monkeypatch.setattr(audit_log, "append", lambda rec, *a, **kw: records.append(rec))
    events = []
    b = ProtocolBroker(lambda name, **f: events.append((name, f)))
    b.events = events
    b.records = records
    return b


def _asked(broker, after=0):
    for _ in range(200):
        asked = [f for n, f in broker.events[after:] if n == "permission.asked"]
        if asked:
            return asked[0]
        time.sleep(0.01)
    raise AssertionError("nothing was asked")


def _ask_and_change(broker, command, changed, tool="bash"):
    """Ask for `command` on a worker thread; answer it with `changed`."""
    out = {}

    def ask():
        out["ok"] = broker.callback(tool, {"command": command}, command)
        out["changed_to"] = broker.last_changed_to

    t = threading.Thread(target=ask)
    t.start()
    rid = _asked(broker)["requestId"]
    out["reply"] = broker.answer(rid, choice="c", command=changed)
    t.join(5)
    return out


def _asks_again(broker, tool, command, preview=None):
    """True when the call opens a new question (it is then refused)."""
    before = len(broker.events)
    out = {}
    t = threading.Thread(target=lambda: out.setdefault(
        "ok", broker.callback(tool, {"command": command},
                              command if preview is None else preview)))
    t.start()
    try:
        rid = _asked(broker, before)["requestId"]
    except AssertionError:
        t.join(5)
        return False
    broker.answer(rid, choice="n")
    t.join(5)
    return True


def test_a_shell_approval_offers_to_change_the_command():
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                            args={"command": "rm -r results"}, preview="rm -r results")
    assert "c" in {o.key for o in tc.options_for(req)}
    write = tc.ConfirmRequest(kind=tc.CONFIRM, tool="write_file",
                              args={"path": "a"}, preview="+x")
    assert "c" not in {o.key for o in tc.options_for(write)}


def test_the_asked_command_is_refused_and_the_change_named(broker):
    out = _ask_and_change(broker, "rm -r results", "rm -r results/old")
    assert out["reply"][0] is True
    assert out["ok"] is False
    assert out["changed_to"] == "rm -r results/old"


def test_the_changed_command_runs_once_without_a_new_question(broker):
    _ask_and_change(broker, "rm -r results", "rm -r results/old")
    before = len(broker.events)
    assert broker.callback("bash", {"command": "rm -r results/old"},
                           "rm -r results/old") is True
    assert not [n for n, _ in broker.events[before:] if n == "permission.asked"]
    assert any(r.get("by") == "changed" and r.get("answer") == "approved"
               for r in broker.records)
    # Once: the same command again is asked.
    assert _asks_again(broker, "bash", "rm -r results/old")


@pytest.mark.parametrize("tool, command, preview", [
    ("bash", "rm -r results/old ", None),                   # one character more
    ("bash_background", "rm -r results/old", None),         # another tool
    ("bash", "rm -r results/old", f"{tc._GUARD_MARKER}\nrm -r results/old"),
    ("bash", "rm -r results/old", f"{tc._CALC_MARKER}\nrm -r results/old"),
    ("bash", "rm -r results/old", f"{tc._OUTSIDE_MARKER}\nrm -r results/old"),
    ("bash", "rm -r results/old", f"{A._WIDENING_BANNER}\nrm -r results/old"),
])
def test_anything_but_the_plain_same_call_is_asked(broker, tool, command, preview):
    _ask_and_change(broker, "rm -r results", "rm -r results/old")
    assert _asks_again(broker, tool, command, preview)


def test_the_grant_expires(broker, monkeypatch):
    _ask_and_change(broker, "rm -r results", "rm -r results/old")
    later = time.monotonic() + tc.CHANGE_GRANT_S + 1
    monkeypatch.setattr(tc.time, "monotonic", lambda: later)
    assert _asks_again(broker, "bash", "rm -r results/old")


@pytest.mark.parametrize("changed", ["", "   ", "rm -r results",
                                     "ls\x1b[2J", "ls\x00", "x" * 4001])
def test_a_change_that_is_no_command_is_refused(broker, changed):
    out = {}
    t = threading.Thread(target=lambda: out.setdefault(
        "ok", broker.callback("bash", {"command": "rm -r results"}, "rm -r results")))
    t.start()
    rid = _asked(broker)["requestId"]
    with pytest.raises(InvalidAnswer):
        broker.answer(rid, choice="c", command=changed)
    assert broker.pending(), "a refused change must leave the question open"
    broker.answer(rid, choice="n")
    t.join(5)
    assert out["ok"] is False and broker.last_changed_to == ""


def test_a_command_with_another_key_or_as_text_is_refused(broker):
    out = {}
    t = threading.Thread(target=lambda: out.setdefault(
        "ok", broker.callback("bash", {"command": "rm -r results"}, "rm -r results")))
    t.start()
    rid = _asked(broker)["requestId"]
    with pytest.raises(InvalidAnswer):
        broker.answer(rid, choice="y", command="ls")
    with pytest.raises(InvalidAnswer):
        broker.answer(rid, choice="c", text="ls")
    broker.answer(rid, choice="n")
    t.join(5)
    assert broker.callback.__self__ is broker
    assert not broker._changed


def test_the_model_reads_the_change_and_must_run_it_as_written(broker, tmp_path, monkeypatch):
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=broker.callback)

    def answer_when_asked():
        rid = _asked(broker)["requestId"]
        broker.answer(rid, choice="c", command="printf changed > out.txt")

    t = threading.Thread(target=answer_when_asked)
    t.start()
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": "some-unlisted-command -x"}, perms)
    t.join(5)
    assert err and "changed the command" in err
    assert "printf changed > out.txt" in err
    # The changed command passes the gate without a new question.
    before = len(broker.events)
    assert A._doc_executor._run_permission_gate(
        "bash", {"command": "printf changed > out.txt"}, perms) is None
    assert not [n for n, _ in broker.events[before:] if n == "permission.asked"]


def test_a_change_never_passes_a_check_that_runs_before_the_question(broker, tmp_path):
    perms = A.KitToolPermissions(workspace=str(tmp_path), mode="default",
                                 confirm_callback=broker.callback)
    _ask_and_change(broker, "some-unlisted-command -x", "rm -rf /")
    err = A._doc_executor._run_permission_gate("bash", {"command": "rm -rf /"}, perms)
    assert err, "a hard-refused command must stay refused after a change"


class _Keys:
    """The terminal's reader: hands out the chunks given, then silence."""

    active = True

    def __init__(self, chunks):
        self._chunks = list(chunks)

    def read_ready(self, timeout):
        return self._chunks.pop(0) if self._chunks else ""


class _Engine:
    token_usage = {"input": 0, "output": 0}

    def get_status(self):
        return {}

    def request_stop(self):
        pass


def test_the_terminal_dialog_edits_the_command_on_one_line(tmp_path, monkeypatch):
    import io

    from delfin.agent import repl
    from delfin.agent import repl_render as rr
    monkeypatch.setattr(tc, "_PENDING_DIR", tmp_path / "pending")
    broker = tc.TerminalConfirmBroker(timeout_s=5.0)
    agent = repl.TerminalAgent(_Engine(), out=io.StringIO(), err=io.StringIO(),
                               broker=broker)
    agent.transcript.theme = rr.Theme(enabled=False)
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                            args={"command": "rm -r results"},
                            preview="$ rm -r results")
    # c, the prefilled line edited (one character off, "/old" on), enter.
    # The dialog reads on until it is answered: run it beside the test,
    # and refuse the request if the keys never answered it.
    t = threading.Thread(target=agent._answer_request,
                         args=(req, _Keys(["c", "\x7f", "s/old", "\r"])),
                         daemon=True)
    t.start()
    t.join(5)
    if t.is_alive():
        broker.resolve(req, False)
        t.join(5)
    assert getattr(req, "changed_to", "") == "rm -r results/old"
    assert req.resolved and req.decision is False
    assert broker.callback("bash", {"command": "rm -r results/old"},
                           "$ rm -r results/old") is True


def test_a_window_gets_the_command_to_change_unless_it_holds_a_secret(broker):
    plain = tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                              args={"command": "rm -r results"}, preview="$ rm -r results")
    assert broker.describe(plain)["command"] == "rm -r results"
    key = "sk-" + "a1B2c3D4e5F6g7H8i9J0k1L2m3N4o5P6"
    req = tc.ConfirmRequest(kind=tc.CONFIRM, tool="bash",
                            args={"command": f"curl -H 'Authorization: Bearer {key}' x"},
                            preview="$ curl ...")
    fields = broker.describe(req)
    # Neither the secret nor its mask is offered to edit (the mask would run).
    assert "command" not in fields
    assert key not in str(fields)
    assert "c" in {o["key"] for o in fields["options"]}
    write = tc.ConfirmRequest(kind=tc.CONFIRM, tool="write_file", args={"path": "a"}, preview="+x")
    assert "command" not in broker.describe(write)
