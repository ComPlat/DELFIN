"""Permission-UX hardening: approval timeout is not a denial, question
timeouts are visible, and allowlist generalisation can't turn one click
into a blanket grant for destructive commands.
"""

from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import kit_settings as ks


# ---------------------------------------------------------------------------
# Timeout vs. denial in the confirm broker
# ---------------------------------------------------------------------------

def test_broker_timeout_sets_flag_and_denies():
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=0.05)
    assert broker.callback("bash", {"command": "ls"}, "$ ls") is False
    assert broker.last_timed_out is True


def test_broker_explicit_deny_clears_flag():
    import threading
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=5.0)

    def _deny_soon():
        import time
        for _ in range(100):
            with broker._lock:
                if broker._pending:
                    broker._pending[0].decision = False
                    broker._pending[0].event.set()
                    return
            time.sleep(0.01)

    t = threading.Thread(target=_deny_soon)
    t.start()
    assert broker.callback("bash", {"command": "ls"}, "$ ls") is False
    t.join()
    assert broker.last_timed_out is False


def test_bash_gate_reports_timeout_distinctly(tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent.api_client import KitToolPermissions
    from delfin.agent.kit_confirm import KitConfirmBroker
    broker = KitConfirmBroker(default_timeout_s=0.05)
    perms = KitToolPermissions(workspace=tmp_path, mode="default")
    perms.confirm_callback = broker.callback
    err = A._doc_executor._run_permission_gate(
        "bash", {"command": "sleep 1"}, perms)
    assert err is not None
    assert "TIMED OUT" in err
    assert "NOT a denial" in err


# ---------------------------------------------------------------------------
# Allowlist generalisation guard
# ---------------------------------------------------------------------------

def test_safe_heads_still_generalize():
    assert ks.suggest_pattern_for_command("git status") == r"^\s*git\s+status\b"
    assert ks.suggest_pattern_for_command("pytest -x tests/") == r"^\s*pytest\b"


def test_destructive_heads_never_generalize_bare():
    pat = ks.suggest_pattern_for_command("rm build/tmp.txt")
    import re
    assert re.match(pat, "rm build/tmp.txt")
    assert not re.match(pat, "rm -rf /home")
    for cmd in ("curl -O http://x/y", "ssh host", "dd if=/dev/zero of=x",
                "bash script.sh"):
        p = ks.suggest_pattern_for_command(cmd)
        head = cmd.split()[0]
        assert p != rf"^\s*{head}\b", cmd


def test_path_invoked_binary_never_generalizes_bare():
    pat = ks.suggest_pattern_for_command("/usr/bin/foo --danger")
    import re
    assert re.match(pat, "/usr/bin/foo --danger")
    assert not re.match(pat, "/usr/bin/foo --other")


def test_persist_pattern_rejects_invalid_regex(tmp_path):
    with pytest.raises(ValueError, match="invalid regex"):
        ks.persist_pattern(r"^\s*(unclosed", kind="allow",
                           user_path=tmp_path / "settings.json")


# ---------------------------------------------------------------------------
# ask_user_question surfaces timed_out
# ---------------------------------------------------------------------------

def test_ask_user_question_surfaces_timeout(tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    perms.ask_user_callback = lambda payload: {
        "answers": [], "timed_out": True}
    out = json.loads(A._doc_executor._execute_ask_user_question(
        {"question": "Proceed?", "options": [
            {"label": "yes"}, {"label": "no"}]}, perms))
    assert out["answers"] == []
    assert out["timed_out"] is True
    assert "Do not re-ask" in out["note"]


def test_ask_user_question_real_answer_has_no_timeout_note(tmp_path):
    from delfin.agent import api_client as A
    from delfin.agent.api_client import KitToolPermissions
    perms = KitToolPermissions(workspace=tmp_path, mode="bypassPermissions")
    perms.ask_user_callback = lambda payload: {"answers": ["yes"]}
    out = json.loads(A._doc_executor._execute_ask_user_question(
        {"question": "Proceed?", "options": [
            {"label": "yes"}, {"label": "no"}]}, perms))
    assert out["answers"] == ["yes"]
    assert "timed_out" not in out
