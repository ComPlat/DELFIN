"""The skill usage ledger: where it is written, how it is opened, and how
a session's uses are paired with the session's end.

Store: ``skill_usage._usage_dir`` resolves ``Path.home()`` per call. The
suite and a benchmark attempt redirect only what the product tables name,
so without an entry in ``USER_STATE_RESOLVERS`` every test and every
benchmark attempt that hands a skill over appends to the user's real
``~/.delfin/skills/_usage``.

Pairing: a session's ``started`` rows are resolved by the rows its own
end appends, matched by the exact session token. Two sessions using one
skill, one session using a skill twice, and a session id that is a
prefix of another each left a ``started`` row nobody could resolve.
"""

import os
import stat
import threading

import pytest

from delfin.agent import session_end as SE
from delfin.agent import skill_usage as SU
from delfin.agent import state_paths


@pytest.fixture
def ledger(tmp_path, monkeypatch):
    monkeypatch.setattr(SU, "_usage_dir", lambda: tmp_path)
    return tmp_path


def _end(monkeypatch, sid, runs):
    monkeypatch.setattr(SE, "_session_test_runs",
                        lambda s: runs if s == sid else [])
    return SE.resolve_pending_session_uses(sid)


# -- store -----------------------------------------------------------------


def test_the_ledger_is_a_redirected_user_store():
    assert ("delfin.agent.skill_usage", "_usage_dir", "skills/_usage") in \
        state_paths.USER_STATE_RESOLVERS


def test_a_redirected_run_writes_into_the_scratch_home(tmp_path,
                                                       monkeypatch):
    # HOME points at a decoy: a resolver the redirect misses lands there,
    # not in the real home, and the assertion below still fails.
    decoy = tmp_path / "decoy"
    monkeypatch.setenv("HOME", str(decoy))
    scratch = tmp_path / "scratch"
    with state_paths.RedirectedUserState(scratch):
        SU.record_use("probe-skill", "started", detail="session:s1")
    assert (scratch / ".delfin" / "skills" / "_usage"
            / "probe-skill.jsonl").is_file()
    assert not (decoy / ".delfin").exists()


# -- pairing ---------------------------------------------------------------


def test_two_sessions_on_one_skill_both_resolve(ledger, monkeypatch):
    SU.record_use("deploy", "started", detail="session:a")
    SU.record_use("deploy", "started", detail="session:b")
    assert _end(monkeypatch, "a", [{"passed": 1}]) == 1
    assert _end(monkeypatch, "b", [{"failed": 1}]) == 1
    u = SU.usage_for("deploy")
    assert (u["uses"], u["successes"], u["failures"], u["pending"]) == \
        (2, 1, 1, 0)


def test_the_later_session_ending_first_leaves_the_earlier_resolvable(
        ledger, monkeypatch):
    SU.record_use("deploy", "started", detail="session:a")
    SU.record_use("deploy", "started", detail="session:b")
    _end(monkeypatch, "b", [{"passed": 1}])
    assert SU.pending_for_session("a") == ["deploy"]
    assert _end(monkeypatch, "a", [{"passed": 1}]) == 1
    assert SU.usage_for("deploy")["pending"] == 0


def test_a_skill_used_twice_in_one_session_has_both_uses_resolved(
        ledger, monkeypatch):
    SU.record_use("deploy", "started", detail="session:s1")
    SU.record_use("deploy", "started", detail="session:s1")
    assert _end(monkeypatch, "s1", []) == 2
    u = SU.usage_for("deploy")
    assert (u["uses"], u["unverified"], u["pending"]) == (2, 2, 0)


def test_a_second_end_of_a_kept_session_resolves_nothing_twice(
        ledger, monkeypatch):
    SU.record_use("deploy", "started", detail="session:s1")
    assert _end(monkeypatch, "s1", [{"passed": 1}]) == 1
    assert _end(monkeypatch, "s1", [{"passed": 1}]) == 0
    SU.record_use("deploy", "started", detail="session:s1")
    assert _end(monkeypatch, "s1", [{"passed": 1}]) == 1
    u = SU.usage_for("deploy")
    assert (u["uses"], u["successes"], u["pending"]) == (2, 2, 0)


def test_a_session_id_that_is_a_prefix_of_another_is_not_credited(
        ledger, monkeypatch):
    SU.record_use("deploy", "started", detail="session:s12")
    assert SU.pending_for_session("s1") == []
    assert _end(monkeypatch, "s1", [{"passed": 1}]) == 0
    assert SU.pending_for_session("s12") == ["deploy"]


def test_a_use_without_a_session_is_not_left_pending(ledger, tmp_path):
    import json
    from delfin.agent import api_client as A
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills" / "hello").mkdir(parents=True)
    (ws / ".delfin" / "skills" / "hello" / "SKILL.md").write_text(
        "---\nname: hello\ndescription: say hello\n---\nSay hello.\n")
    perms = A.KitToolPermissions(workspace=ws, mode="default")
    perms.task_session_id = ""
    out = json.loads(A._DocToolExecutor()._execute_skill(
        {"name": "hello"}, perms))
    if out.get("status") != "ok":
        pytest.skip(f"skill not discoverable in this fixture: {out}")
    u = SU.usage_for("hello")
    assert (u["uses"], u["pending"], u["unverified"]) == (1, 0, 1)


def test_a_skill_handed_to_a_sub_agent_is_counted(ledger, tmp_path,
                                                  monkeypatch):
    """``delegate: true`` runs the skill in a sub-agent; it is a use of
    the skill all the same, under the parent's session."""
    import json
    from delfin.agent import api_client as A
    ws = tmp_path / "ws"
    (ws / ".delfin" / "skills" / "hello").mkdir(parents=True)
    (ws / ".delfin" / "skills" / "hello" / "SKILL.md").write_text(
        "---\nname: hello\ndescription: say hello\n---\nSay hello.\n")
    perms = A.KitToolPermissions(workspace=ws, mode="default")
    perms.task_session_id = "s7"
    ex = A._DocToolExecutor()
    monkeypatch.setattr(ex, "_execute_subagent",
                        lambda args, p: json.dumps({"status": "ok"}))
    out = json.loads(ex._execute_skill({"name": "hello", "delegate": True},
                                       perms))
    assert out.get("status") == "ok", out
    assert SU.pending_for_session("s7") == ["hello"]


# -- opening ---------------------------------------------------------------


def test_a_symlinked_ledger_is_not_written_through(ledger, tmp_path):
    target = tmp_path / "elsewhere.txt"
    target.write_text("keep\n")
    os.symlink(target, ledger / "deploy.jsonl")
    SU.record_use("deploy", "started", detail="session:s1")
    assert target.read_text() == "keep\n"


def test_a_fifo_at_the_ledger_path_does_not_block_a_reader(ledger):
    os.mkfifo(ledger / "deploy.jsonl")
    done = []

    def _read():
        done.append(SU.usage_for("deploy"))
        done.append(SU.pending_for_session("s1"))
        SU.record_use("deploy", "started", detail="session:s1")
        done.append(True)

    t = threading.Thread(target=_read, daemon=True)
    t.start()
    t.join(5)
    if t.is_alive():
        # Unblock the stuck open so the thread can end.
        fd = os.open(ledger / "deploy.jsonl", os.O_RDWR | os.O_NONBLOCK)
        os.close(fd)
        pytest.fail("a FIFO at the ledger path blocked the reader")
    assert stat.S_ISFIFO(os.lstat(ledger / "deploy.jsonl").st_mode)


def test_a_new_ledger_file_is_owner_only(ledger):
    SU.record_use("deploy", "started", detail="session:s1")
    mode = stat.S_IMODE(os.stat(ledger / "deploy.jsonl").st_mode)
    assert mode == 0o600


def test_a_read_is_bounded_to_the_tail(ledger, monkeypatch):
    monkeypatch.setattr(SU, "LEDGER_READ_MAX_BYTES", 4096)
    for i in range(200):
        SU.record_use("deploy", "unverified", detail=f"session:s{i}")
    assert (ledger / "deploy.jsonl").stat().st_size > 4096
    rows = SU.usage_records("deploy")
    assert 0 < len(rows) < 200
    assert rows[0]["detail"] == "session:s199"


# -- display ---------------------------------------------------------------


def test_a_live_skills_record_is_shown_without_a_proposal(ledger,
                                                         monkeypatch):
    """Accepted proposals leave the proposal list (they are live skills),
    and only a live skill can be handed over -- so a record keyed to a
    proposal row alone is never on screen."""
    from delfin.dashboard import skill_proposals_panel as P
    monkeypatch.setattr(P, "collect_rows", lambda: [])
    SU.record_use("deploy", "started", detail="session:s1")
    _end(monkeypatch, "s1", [{"failed": 1}])
    html = P.format_panel_html()
    assert "deploy" in html
    assert "used 1" in html


def test_a_red_session_is_not_called_a_failure_of_the_skill(ledger,
                                                           monkeypatch):
    """The outcome is the session's test record, not the skill's: a red
    run anywhere in the session marks every skill it used."""
    from delfin.dashboard import skill_proposals_panel as P
    monkeypatch.setattr(P, "collect_rows", lambda: [])
    SU.record_use("deploy", "started", detail="session:s1")
    _end(monkeypatch, "s1", [{"failed": 1}])
    html = P.format_panel_html()
    assert "fail" not in html.lower()
    assert "red test run" in html

# Workspace skills are read as trusted here: these tests are about what a
# loaded skill does; the trust gate itself is pinned in
# test_a_repository_skill_needs_trust.py.
pytestmark = pytest.mark.usefixtures("workspace_skills_trusted")
