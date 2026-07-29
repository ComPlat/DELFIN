"""Self-preservation: process-kill commands must not hit the host stack.

Field case 20260729-112242: the agent ran ``pkill -f "voila.*8866"`` to
stop "old Voila instances" and killed the user's own DELFIN dashboard
(hosted under Voila on the same port), losing the whole unpersisted turn.
"""

from __future__ import annotations

import json

from delfin.agent.api_client import (
    _DocToolExecutor,
    KitToolPermissions,
    _kill_targets_host_process,
)


_ANC = [
    (4242, "python3", "/usr/bin/python3 -m ipykernel_launcher -f kernel.json"),
    (4100, "voila", "/pfs/home/user/.venv/bin/voila dashboard.ipynb "
                    "--port=8866 --no-browser"),
    (3999, "bash", "-bash"),
]


def test_field_case_pkill_pattern_is_blocked():
    reason = _kill_targets_host_process(
        'pkill -f "voila.*8866" || true', ancestry=_ANC)
    assert "voila" in reason and "host process" in reason


def test_killall_by_ancestor_name_is_blocked():
    assert _kill_targets_host_process("killall voila", ancestry=_ANC)
    assert _kill_targets_host_process("killall -9 voila", ancestry=_ANC)


def test_kill_by_ancestor_pid_is_blocked():
    assert _kill_targets_host_process("kill -9 4100", ancestry=_ANC)


def test_unrelated_kills_are_allowed():
    assert _kill_targets_host_process(
        'pkill -f "my_test_server"', ancestry=_ANC) == ""
    assert _kill_targets_host_process("killall orca", ancestry=_ANC) == ""
    assert _kill_targets_host_process("kill 99999", ancestry=_ANC) == ""


def test_pkill_without_dash_f_matches_process_name_not_cmdline():
    # Without -f, pkill matches the NAME; 'voila' hits the ancestor name.
    assert _kill_targets_host_process("pkill voila", ancestry=_ANC)
    # '8866' appears only in the cmdline, not the name -> no name match.
    assert _kill_targets_host_process("pkill 8866", ancestry=_ANC) == ""


def test_no_ancestry_falls_back_to_host_stack_names():
    assert _kill_targets_host_process('pkill -f "voila"', ancestry=[])
    assert _kill_targets_host_process("killall jupyter-lab", ancestry=[])
    assert _kill_targets_host_process(
        'pkill -f "my_test_server"', ancestry=[]) == ""


def test_non_kill_commands_pass_through():
    assert _kill_targets_host_process("ls -la && pytest -q", ancestry=_ANC) == ""
    assert _kill_targets_host_process(
        "echo 'do not pkill anything'", ancestry=_ANC) == ""


def test_chained_kill_segment_is_still_caught():
    cmd = "cd game && curl localhost:8866; pkill -f 'voila.*8866'"
    assert _kill_targets_host_process(cmd, ancestry=_ANC)


def test_bad_regex_pattern_never_raises():
    assert _kill_targets_host_process(
        "pkill -f '[unclosed'", ancestry=_ANC) == ""


def test_executor_blocks_with_guidance(tmp_path, monkeypatch):
    import delfin.agent.api_client as ac
    monkeypatch.setattr(ac, "_host_process_ancestry", lambda: _ANC)
    ws = tmp_path / "ws"
    ws.mkdir()
    ex = _DocToolExecutor()
    perms = KitToolPermissions(workspace=str(ws))
    perms.mode = "bypassPermissions"      # protection holds in every mode
    out = json.loads(ex.execute(
        "bash", {"command": 'pkill -f "voila.*8866"'}, perms))
    assert "blocked" in out.get("error", "")
    assert "bash_kill" in out["error"]
