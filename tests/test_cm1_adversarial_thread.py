"""CM1 reviewer tests for operator hook-build item 2: the RPC serve loop
runs in a daemon THREAD while the agent thread is blocked in the bash
call. Whatever execute() relies on must survive being called from that
thread: the same permissions object, no cross-talk with the blocked main
thread, and (for the hook) a copied context if contextvars are involved.

This file pins the code_mode side (server + Connection carry the exact
permissions object into execute, from a non-main thread). The real
KitToolExecutor + contextvars wiring is the operator's build; the px
comment fix (_dispatch_under_rules -> execute) is nacht-s9's.
"""
import contextvars
import json
import threading

from delfin.agent import code_mode


class _ThreadRecordingExecutor:
    """Records the calling thread and the permissions identity it saw."""

    def __init__(self, marker):
        self.marker = marker
        self.saw_thread = None
        self.saw_permissions = None

    def execute(self, name, arguments, permissions):
        self.saw_thread = threading.current_thread()
        self.saw_permissions = permissions
        return json.dumps({"result": "ok"})


def _serve_one_in_thread(server, payload):
    """handle_payload driven exactly like the hook drives it: from a
    plain daemon thread, while the 'agent thread' (main here) waits."""
    box = {}

    def _run():
        box["reply"] = server.handle_payload(payload)

    thread = threading.Thread(target=_run, daemon=True)
    thread.start()
    thread.join(timeout=10.0)
    return box["reply"]


def test_execute_from_serve_thread_sees_the_same_permissions_object():
    """The permissions object must survive into execute() when the call
    arrives from the serve thread — identity, not a copy: a copy would
    let a permission flip land between the two views."""
    marker = object()
    ex = _ThreadRecordingExecutor(marker)
    server = code_mode.RpcServer(ex, permissions=marker)
    reply = _serve_one_in_thread(
        server, json.dumps({"tool": "get_calc_info", "arguments": {}}))
    assert reply == json.dumps({"result": {"result": "ok"}}), reply
    assert ex.saw_thread is not threading.main_thread(), (
        "call did not come from a serve thread")
    assert ex.saw_permissions is marker, (
        "execute() saw a different permissions object than the server "
        "was constructed with")


def test_contextvars_do_not_leak_or_vanish_across_the_serve_thread():
    """Documents the hazard the operator named for the hook: a contextvar
    set in the agent thread is NOT visible inside a plain daemon thread
    (threading.Thread does not inherit context). The hook must copy the
    context (contextvars.copy_context().run) if anything downstream reads
    one — pinned here so the contract is written down next to the code
    that must respect it."""
    var = contextvars.ContextVar("cm1_review_probe", default=None)

    def _read():
        _read.value_in_thread = var.get()

    _read.value_in_thread = "unset"
    var.set("main-thread-value")
    thread = threading.Thread(target=_read, daemon=True)
    thread.start()
    thread.join(timeout=10.0)
    assert _read.value_in_thread is None, (
        "a plain thread unexpectedly inherited the caller's context "
        "value — the hazard this documents does not reproduce here")
