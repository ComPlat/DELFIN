"""The TUI opens no session by mistake (Claude Code's conventions):

* ``delfin-agent tui --resume <id>`` with an id nobody knows is an error
  that names the close ones -- never a new session under that name
  (``claude --resume`` says "No conversation found"); a unique prefix of a
  known id stands for it, as in ``service end``;
* a worker asked to resume an unknown id refuses (chat would say "not
  found, starting fresh" and start one under that id);
* a session that was never asked anything is not saved when it ends, so
  every window opened and closed does not leave an "(untitled)" behind.
"""

from __future__ import annotations

import argparse
import json

import pytest

from delfin.agent import session_store
from delfin.agent.service import cli_service, worker

from tui_worker_client import WorkerClient


@pytest.fixture
def store(tmp_path, monkeypatch):
    """A session store of its own, with two saved sessions."""
    d = tmp_path / "agent_sessions"
    d.mkdir()
    monkeypatch.setattr(session_store, "_SESSIONS_DIR", d)
    for sid, title in (("0f3a9c21aa77", "relax Fe complex"), ("0f3b1111bb88", "benzene scan")):
        (d / f"{sid}.json").write_text(json.dumps({"session_id": sid, "title": title,
                                                   "chat_messages": [{"role": "user", "content": title}]}))
    return d


HELLO = {"pid": 1, "sessions": [{"sessionId": "live5555cc99", "title": "running one", "state": "idle"}]}


# -- resolve_resume ---------------------------------------------------------------------

def test_a_known_id_is_opened_as_given(store):
    assert cli_service.resolve_resume("0f3a9c21aa77", HELLO) == ("0f3a9c21aa77", "")
    assert cli_service.resolve_resume("live5555cc99", HELLO) == ("live5555cc99", "")


def test_a_unique_prefix_stands_for_its_id(store):
    assert cli_service.resolve_resume("0f3a9", HELLO) == ("0f3a9c21aa77", "")
    assert cli_service.resolve_resume("live", HELLO) == ("live5555cc99", "")


def test_a_prefix_too_short_stands_for_no_id(store):
    # Shorter than MIN_PREFIX: an unknown id, even where only one id starts so.
    sid, why = cli_service.resolve_resume("liv", HELLO)
    assert sid == "" and "No session found with id liv" in why


def test_two_sessions_under_one_prefix_are_both_named(store):
    (store / "0f3a9c21ffff.json").write_text(json.dumps({"session_id": "0f3a9c21ffff", "title": "x"}))
    sid, why = cli_service.resolve_resume("0f3a9c21", HELLO)
    assert sid == ""
    assert "0f3a9c21aa77" in why and "0f3a9c21ffff" in why and "give more of the id" in why


def test_an_unknown_id_is_an_error_naming_the_close_ones(store):
    sid, why = cli_service.resolve_resume("0f3a9c21aa78", HELLO)
    assert sid == ""
    assert "No session found with id 0f3a9c21aa78" in why
    assert "Did you mean" in why and "0f3a9c21aa77" in why and "relax Fe complex" in why
    assert "--resume` (no id) opens the picker" in why


def test_an_unknown_id_with_nothing_close_says_where_to_look(store):
    sid, why = cli_service.resolve_resume("doesnotexist", HELLO)
    assert sid == ""
    assert "No session found with id doesnotexist." in why and "Did you mean" not in why
    assert "delfin-agent sessions" in why


def test_a_title_typed_for_an_id_names_its_session(store):
    sid, why = cli_service.resolve_resume("benzene", HELLO)
    assert sid == "" and "0f3b1111bb88" in why


def test_a_path_is_never_looked_up(store, tmp_path):
    (tmp_path / "outside.json").write_text(json.dumps({"session_id": "outside"}))
    sid, why = cli_service.resolve_resume("../outside", HELLO)
    assert sid == "" and "No session found" in why


# -- the launcher -----------------------------------------------------------------------

def _launch(monkeypatch, argv):
    from delfin.agent import cli
    from delfin.agent.service import launcher

    handed: list[dict] = []

    class Exec(Exception):
        pass

    def execvpe(path, args, env):
        handed.append(dict(env))
        raise Exec()

    monkeypatch.setattr(cli_service, "_security_preamble", lambda: None)
    monkeypatch.setattr(cli_service, "_has_terminal", lambda: True)
    monkeypatch.setattr(launcher, "ensure_service", lambda: HELLO)
    monkeypatch.setattr(launcher, "tui_command", lambda: ["/bin/true"])
    monkeypatch.setattr(launcher, "service_command", lambda: ["svc"])
    monkeypatch.setattr(launcher, "other_host_notes", lambda: [])
    monkeypatch.setattr(cli_service.os, "execvpe", execvpe)
    try:
        code = cli_service.cmd_tui(cli.build_parser().parse_args(argv))
    except Exec:
        code = None
    return code, handed


def test_tui_resume_with_an_unknown_id_exits_2_and_opens_no_window(store, monkeypatch, tmp_path, capsys):
    code, handed = _launch(monkeypatch, ["tui", "--cwd", str(tmp_path), "--resume", "doesnotexist"])
    assert code == 2 and handed == []
    err = capsys.readouterr().err
    assert "No session found with id doesnotexist" in err


def test_tui_resume_with_a_prefix_hands_the_whole_id_to_the_window(store, monkeypatch, tmp_path):
    code, handed = _launch(monkeypatch, ["tui", "--cwd", str(tmp_path), "--resume", "0f3a9"])
    assert code is None
    assert json.loads(handed[-1]["DELFIN_TUI_START"])["resume"] == "0f3a9c21aa77"


def test_tui_without_resume_is_not_checked(store, monkeypatch, tmp_path):
    code, handed = _launch(monkeypatch, ["tui", "--cwd", str(tmp_path)])
    assert code is None and "resume" not in json.loads(handed[-1]["DELFIN_TUI_START"])


# -- the worker -------------------------------------------------------------------------

def test_the_worker_resumes_only_a_saved_session(store):
    assert worker.resumable("0f3a9c21aa77")
    assert worker.resumable("latest")                    # chat's own word
    assert not worker.resumable("doesnotexist")
    assert not worker.resumable("../agent_sessions/0f3a9c21aa77")
    assert not worker.resumable("")


def test_a_session_the_service_ran_is_resumable_without_a_save(store, tmp_path, monkeypatch):
    # Put to sleep before anything was asked of it: not saved, but its
    # event log is there, and it wakes under its id.
    monkeypatch.setattr(worker, "log_path", lambda sid: tmp_path / "sessions" / f"{sid}.jsonl")
    (tmp_path / "sessions").mkdir()
    (tmp_path / "sessions" / "asleep000001.jsonl").write_text("{}\n")
    assert worker.resumable("asleep000001")
    assert not worker.resumable("asleep000002")


def test_the_worker_refuses_to_start_an_unknown_id_as_a_new_session(store, tmp_path):
    with pytest.raises(worker.SessionRefused, match="No saved session doesnotexist"):
        worker.open_real_session(argparse.Namespace(resume="doesnotexist"), tmp_path)


def _saved_files(w: WorkerClient, sid: str):
    return [p for root in (w.home, w.scratch) for p in root.rglob(f"{sid}.json")]


def test_a_session_never_asked_anything_is_not_saved(tmp_path):
    with WorkerClient(tmp_path, "--fake-save") as w:
        sid = w.init_result["session"]["sessionId"]
    assert _saved_files(w, sid) == []


def test_a_session_that_was_asked_is_saved(tmp_path):
    with WorkerClient(tmp_path, "--fake-save") as w:
        w.turn("hello")
        sid = w.init_result["session"]["sessionId"]
    assert _saved_files(w, sid), "the session that was asked something was not saved"


# -- no terminal ------------------------------------------------------------------------

def test_tui_without_a_terminal_exits_2_and_names_print_mode(store, monkeypatch, tmp_path, capsys):
    """As Claude Code: no interactive window in a pipe; -p is the way there."""
    from delfin.agent import cli
    from delfin.agent.service import launcher

    monkeypatch.setattr(cli_service, "_has_terminal", lambda: False)
    monkeypatch.setattr(cli_service, "_security_preamble", lambda: None)
    monkeypatch.setattr(launcher, "ensure_service", lambda: HELLO)
    monkeypatch.setattr(launcher, "tui_command", lambda: ["/bin/true"])
    monkeypatch.setattr(cli_service.os, "execvpe",
                        lambda *a: pytest.fail("the window was started without a terminal"))
    code = cli_service.cmd_tui(cli.build_parser().parse_args(["tui", "--cwd", str(tmp_path)]))
    assert code == 2
    assert "needs a terminal" in capsys.readouterr().err


def test_a_closed_stdin_counts_as_no_terminal(monkeypatch):
    class Closed:
        def isatty(self):
            raise ValueError("I/O operation on closed file")

    monkeypatch.setattr(cli_service.sys, "stdin", Closed())
    assert cli_service._has_terminal() is False
