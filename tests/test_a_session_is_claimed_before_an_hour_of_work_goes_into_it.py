"""Sessions: claimed up front, scoped to this directory, and findable.

THE LOCK. `save_session` refuses to overwrite a session another live
process is writing, and the headless saver swallows that refusal into one
WARN line on stderr. That is right for a scheduled run and wrong for a
person: they work for an hour and discover at the end that nothing was
written. Claiming at second zero costs nothing and can say what to do
instead.

THE SCOPE. `run --session latest` resolved through `resume_latest()`,
which asks the whole machine rather than this directory — so continuing
project A could pick up project B's conversation and answer out of it.
`latest_session(workspace=…)` exists for exactly this and its own
docstring says so.

THE TITLE. Every session the headless path ever wrote passed
`chat_messages=[]`, and the store derives the title from the first user
message — so they are all "Untitled" and none of them can be found by
`session search`.
"""

from __future__ import annotations

import argparse

from delfin.agent import cli as agent_cli
from delfin.agent import launch_guard as lg


# ---------------------------------------------------------------------------
# The lock
# ---------------------------------------------------------------------------

def test_a_session_another_process_is_writing_is_refused_at_the_start(
        monkeypatch, capsys):
    from delfin.agent import session_store as ss

    def _locked(sid):
        raise ss.SessionLockedError(sid, 4242)

    monkeypatch.setattr(ss, "acquire_session_lock", _locked)
    assert agent_cli._claim_session("sid-x") is False

    err = capsys.readouterr().err
    assert "4242" in err
    assert "--fork-session" in err, "a refusal has to say what to do instead"
    assert "--new" in err


def test_a_free_session_is_claimed_and_released(monkeypatch):
    from delfin.agent import session_store as ss
    taken: list[str] = []
    monkeypatch.setattr(ss, "acquire_session_lock", lambda sid: taken.append(sid))
    registered: list = []
    monkeypatch.setattr("atexit.register",
                        lambda fn, *a: registered.append((fn, a)))

    assert agent_cli._claim_session("sid-y") is True
    assert taken == ["sid-y"]
    assert registered and registered[0][1] == ("sid-y",), (
        "the lock has to be given back when the process ends")


def test_a_lock_that_cannot_be_taken_at_all_does_not_block_the_start(
        monkeypatch):
    """A missing state directory is not another writer."""
    from delfin.agent import session_store as ss

    def _boom(sid):
        raise OSError("no such directory")

    monkeypatch.setattr(ss, "acquire_session_lock", _boom)
    assert agent_cli._claim_session("sid-z") is True


def test_a_fresh_session_needs_no_claim():
    assert agent_cli._claim_session("") is True


# ---------------------------------------------------------------------------
# The scope
# ---------------------------------------------------------------------------

def test_latest_means_this_directorys_last_conversation(monkeypatch, tmp_path):
    from delfin.agent import session_store as ss
    asked: dict = {}

    def _latest(workspace=None):
        asked["workspace"] = workspace
        return {"session_id": "here-1"}

    monkeypatch.setattr(ss, "latest_session", _latest)
    monkeypatch.setattr(ss, "resume_latest",
                        lambda **kw: {"session_id": "somewhere-else"})

    restored: dict = {}

    class _Engine:
        session_id = ""

        def restore_state(self, data):
            restored.update(data)
            return None

    args = argparse.Namespace(session="latest", mode="solo",
                              cwd=str(tmp_path), verbose=False)
    agent_cli._resume_or_create(_Engine(), args)

    assert asked["workspace"] == str(tmp_path)
    assert restored.get("session_id") == "here-1", (
        "resume_latest asks the whole machine; this must ask the project")


def test_the_machine_wide_fallback_still_exists(monkeypatch, tmp_path):
    """A session saved before workspaces were recorded must stay reachable."""
    from delfin.agent import session_store as ss
    monkeypatch.setattr(ss, "latest_session", lambda workspace=None: None)
    monkeypatch.setattr(ss, "resume_latest",
                        lambda **kw: {"session_id": "old-one"})

    restored: dict = {}

    class _Engine:
        session_id = ""

        def restore_state(self, data):
            restored.update(data)
            return None

    agent_cli._resume_or_create(
        _Engine(), argparse.Namespace(session="latest", mode="solo",
                                      cwd=str(tmp_path), verbose=False))
    assert restored.get("session_id") == "old-one"


# ---------------------------------------------------------------------------
# The title
# ---------------------------------------------------------------------------

def test_a_saved_session_carries_what_was_asked():
    class _Engine:
        messages = [
            {"role": "system", "content": "you are an agent"},
            {"role": "user", "content": "fix the failing test in calc.py"},
            {"role": "assistant", "content": "I will look."},
        ]

    rows = agent_cli._display_messages(_Engine())
    assert rows[0] == {"role": "user",
                       "content": "fix the failing test in calc.py"}
    assert all(r["role"] in ("user", "assistant") for r in rows)


def test_a_content_block_list_becomes_text_not_a_crash():
    """The store slices the first user message to build a title.

    A content list would make that slice a LIST slice, and .replace on it
    takes the whole save down.
    """
    class _Engine:
        messages = [{"role": "user", "content": [
            {"type": "text", "text": "check the sheet"},
            {"type": "tool_result", "content": "..."},
        ]}]

    rows = agent_cli._display_messages(_Engine())
    assert rows == [{"role": "user", "content": "check the sheet"}]
    assert isinstance(rows[0]["content"], str)


def test_an_engine_with_no_history_saves_nothing_odd():
    class _Engine:
        messages = None

    assert agent_cli._display_messages(_Engine()) == []


# ---------------------------------------------------------------------------
# The flags
# ---------------------------------------------------------------------------

def test_the_session_flags_parse_the_way_they_read():
    p = agent_cli.build_parser()

    args = p.parse_args(["chat", "-c"])
    assert args.continue_session is True
    assert args.resume is None

    args = p.parse_args(["chat", "-r"])
    assert args.resume == "", "a bare -r means: show me the list"

    args = p.parse_args(["chat", "-r", "abc123"])
    assert args.resume == "abc123"

    args = p.parse_args(["chat", "--new"])
    assert args.new_session is True

    args = p.parse_args(["chat", "-r", "abc", "--fork-session"])
    assert args.fork_session is True

    args = p.parse_args(["chat", "-n", "the migration"])
    assert args.session_name == "the migration"


def test_a_fresh_session_is_the_default():
    args = agent_cli.build_parser().parse_args(["chat"])
    assert args.continue_session is False
    assert args.resume is None, (
        "a terminal opened in a directory must not silently inherit a "
        "conversation")


def test_forking_leaves_the_original_alone(monkeypatch, tmp_path, capsys):
    from delfin.agent import session_store as ss
    monkeypatch.setattr(ss, "fork_session", lambda sid: {"session_id": "copy-1"})
    monkeypatch.setattr(ss, "acquire_session_lock", lambda sid: None)
    monkeypatch.setattr(ss, "consume_crash_recovery_note", lambda sid: None)
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda e, a: "")

    class _Engine:
        session_id = "copy-1"

    args = argparse.Namespace(resume="orig-1", continue_session=False,
                              new_session=False, fork_session=True,
                              session="", mode="solo", cwd=str(tmp_path),
                              verbose=False)
    assert agent_cli._open_session(_Engine(), args, tmp_path) is True
    assert args.session == "copy-1"
    assert "untouched" in capsys.readouterr().err


def test_a_crash_note_from_last_time_is_handed_over(monkeypatch, tmp_path,
                                                    capsys):
    from delfin.agent import session_store as ss
    monkeypatch.setattr(ss, "acquire_session_lock", lambda sid: None)
    monkeypatch.setattr(ss, "consume_crash_recovery_note",
                        lambda sid: "the last turn was cut short mid-tool")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda e, a: "")

    class _Engine:
        session_id = "sid-1"

    args = argparse.Namespace(resume="sid-1", continue_session=False,
                              new_session=False, fork_session=False,
                              session="", mode="solo", cwd=str(tmp_path),
                              verbose=False)
    agent_cli._open_session(_Engine(), args, tmp_path)
    assert "cut short" in capsys.readouterr().err


# ---------------------------------------------------------------------------
# The banner's own facts
# ---------------------------------------------------------------------------

def test_the_first_changed_file_keeps_its_first_letter(monkeypatch, tmp_path):
    """`git status --porcelain` puts the status in columns 1-2.

    A convenience .strip() on the whole output ate the leading space of the
    FIRST line and shifted its path by one character, so a banner
    truthfully counting two modified files named one of them
    `elfin/agent/cli.py`. Only the first line was ever wrong, which is why
    it reads as a typo rather than a parse.
    """
    calls: list[tuple] = []

    class _Proc:
        returncode = 0

        def __init__(self, out):
            self.stdout = out

    def _run(cmd, **kw):
        calls.append(tuple(cmd))
        if cmd[1:] == ["rev-parse", "--show-toplevel"]:
            return _Proc(str(tmp_path) + "\n")
        if cmd[1:] == ["rev-parse", "--abbrev-ref", "HEAD"]:
            return _Proc("main\n")
        if cmd[1:] == ["status", "--porcelain"]:
            return _Proc(" M delfin/agent/cli.py\n M delfin/agent/repl.py\n"
                         "?? notes.txt\n")
        return _Proc("")

    monkeypatch.setattr("subprocess.run", _run)
    info = lg._git_info(tmp_path)
    assert info.dirty == ("delfin/agent/cli.py", "delfin/agent/repl.py",
                          "notes.txt")
    assert info.branch == "main"


def test_a_renamed_file_is_named_by_where_it_ended_up(monkeypatch, tmp_path):
    class _Proc:
        returncode = 0

        def __init__(self, out):
            self.stdout = out

    def _run(cmd, **kw):
        if cmd[1:] == ["rev-parse", "--show-toplevel"]:
            return _Proc(str(tmp_path) + "\n")
        if cmd[1:] == ["status", "--porcelain"]:
            return _Proc("R  old_name.py -> new_name.py\n")
        return _Proc("")

    monkeypatch.setattr("subprocess.run", _run)
    assert lg._git_info(tmp_path).dirty == ("new_name.py",)
