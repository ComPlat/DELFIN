"""Bug 065503: a plain "Erlauben (1×)" on an outside-workspace access grants
the directory for THIS session (live perms, non-persisted) so the agent stops
re-prompting per file. 'Dauerhaft' still persists; bash 1× grants nothing.
"""

from __future__ import annotations

import threading
import time

import pytest

widgets = pytest.importorskip("ipywidgets")

from delfin.agent.kit_confirm import KitConfirmBroker


def _broker_with_recorder():
    broker = KitConfirmBroker(default_timeout_s=5.0)
    calls: list[tuple] = []

    def _persist(kind, value):
        calls.append((kind, value))
        return True, "ok"

    broker.set_persist_callback(_persist)
    return broker, calls


def _run_request(broker, tool, args):
    out = {}

    def _r():
        out["ok"] = broker.callback(tool, args, preview="")

    t = threading.Thread(target=_r, daemon=True)
    t.start()
    for _ in range(100):                       # wait until it's pending
        if broker._pending:
            break
        time.sleep(0.01)
    req = list(broker._pending)[0]
    return t, req, out


def _buttons(broker, req):
    row = broker._build_request_row(req, widgets)
    hbox = row.children[-1]                     # VBox -> [..., HBox(buttons)]
    approve, approve_persist, deny = hbox.children
    return approve, approve_persist, deny


def test_plain_allow_grants_dir_for_session(tmp_path):
    broker, calls = _broker_with_recorder()
    f = tmp_path / "data" / "mod.py"
    f.parent.mkdir(parents=True)
    f.write_text("x")
    t, req, out = _run_request(broker, "read_file", {"path": str(f)})
    approve, _persist_btn, _deny = _buttons(broker, req)
    approve.click()                            # "Erlauben (1×)"
    t.join(timeout=5)

    parent = str(f.parent.resolve())
    assert out["ok"] is True
    # session grant fired with the parent dir, NOT a persisted extra_dir
    assert ("extra_dir_session", parent) in calls
    assert all(k != "extra_dir" for k, _ in calls)


def test_dauerhaft_persists_and_does_not_double_grant(tmp_path):
    broker, calls = _broker_with_recorder()
    f = tmp_path / "docs" / "readme.md"
    f.parent.mkdir(parents=True)
    f.write_text("x")
    t, req, out = _run_request(broker, "read_file", {"path": str(f)})
    _approve, approve_persist, _deny = _buttons(broker, req)
    approve_persist.click()                    # "Erlauben + Dauerhaft"
    t.join(timeout=5)

    parent = str(f.parent.resolve())
    assert out["ok"] is True
    assert ("extra_dir", parent) in calls          # persisted
    assert all(k != "extra_dir_session" for k, _ in calls)   # no double grant


def test_bash_allow_once_grants_no_directory(tmp_path):
    broker, calls = _broker_with_recorder()
    t, req, out = _run_request(broker, "bash", {"command": "ls -la /etc"})
    approve, _persist_btn, _deny = _buttons(broker, req)
    approve.click()                            # 1× on a bash command
    t.join(timeout=5)

    assert out["ok"] is True
    # A bash command grants no directory — only file-access requests do.
    assert all(k not in ("extra_dir_session", "extra_dir") for k, _ in calls)


def test_deny_grants_nothing(tmp_path):
    broker, calls = _broker_with_recorder()
    f = tmp_path / "secret" / "k.py"
    f.parent.mkdir(parents=True)
    f.write_text("x")
    t, req, out = _run_request(broker, "read_file", {"path": str(f)})
    _approve, _persist_btn, deny = _buttons(broker, req)
    deny.click()
    t.join(timeout=5)

    assert out["ok"] is False
    assert calls == []


def test_allow_permanent_hidden_when_no_persist_hook():
    """The 'Allow + Permanent' button is only shown when something is actually
    persistable AND a persist hook is wired. Otherwise it used to render
    DISABLED — a dead button that did nothing on click and left the dialog up
    (user bug 2026-06-25: 'Dauerhaft geht garnicht, Fenster verschwindet
    nicht'). Now only acting buttons are shown, and any click dismisses."""
    broker = KitConfirmBroker(default_timeout_s=5.0)   # NO persist callback
    t, req, out = _run_request(broker, "bash", {"command": "ls -la /etc"})
    row = broker._build_request_row(req, widgets)
    labels = [b.description for b in row.children[-1].children]
    assert "Allow + Permanent" not in labels      # dead button is gone
    assert any("Allow" in lbl for lbl in labels)  # plain allow stays
    assert any("Deny" in lbl for lbl in labels)
    # The plain Allow still resolves AND removes the row from the panel.
    row.children[-1].children[0].click()
    t.join(timeout=5)
    assert out["ok"] is True
    assert req not in list(broker._pending)        # dialog dismissed


# --- Adversarial-review fix: don't auto-grant sensitive system/secret dirs ---

@pytest.mark.parametrize("path,grantable", [
    ("/tmp/project", True),
    ("/home/u/ka_xn0397/Porpoise", True),
    ("/etc", False), ("/etc/hostname", False),   # system root + child
    ("/root", False), ("/", False), ("/usr/lib", False), ("/var/log", False),
    ("/home/u/.ssh", False),                     # secret dir
])
def test_is_grantable_session_dir(path, grantable):
    from delfin.agent.kit_confirm import _is_grantable_session_dir
    assert _is_grantable_session_dir(path) is grantable


def test_approving_a_system_file_does_not_grant_its_dir(tmp_path):
    # Approving a read of /etc/hostname must NOT grant all of /etc for the
    # session (the single read still proceeds; we just don't broaden).
    broker, calls = _broker_with_recorder()
    t, req, out = _run_request(broker, "read_file", {"path": "/etc/hostname"})
    approve, _persist_btn, _deny = _buttons(broker, req)
    approve.click()
    t.join(timeout=5)
    assert out["ok"] is True                         # the read is allowed
    assert all(k != "extra_dir_session" for k, _ in calls)   # but no broad grant
