"""Tests for the headless CLI entrypoint delfin.agent.cli."""

from __future__ import annotations

import io
import json
import sys
from contextlib import redirect_stderr, redirect_stdout
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


def test_build_parser_accepts_run_init_session():
    p = agent_cli.build_parser()
    args = p.parse_args(["run", "hello", "world"])
    assert args.cmd == "run"
    assert args.prompt == ["hello", "world"]
    args = p.parse_args(["init", "/tmp/foo", "--force"])
    assert args.cmd == "init"
    assert args.path == "/tmp/foo"
    assert args.force is True
    args = p.parse_args(["session", "ls", "--limit", "5"])
    assert args.cmd == "session"
    assert args.session_action == "ls"
    assert args.limit == 5


def test_cmd_init_writes_files(tmp_path, capsys):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n",
                                              encoding="utf-8")
    rc = agent_cli.main(["init", str(tmp_path)])
    assert rc == 0
    assert (tmp_path / "AGENTS.md").exists()
    out, _ = capsys.readouterr()
    assert "Python project" in out
    assert "Created:" in out


def test_cmd_init_skips_existing_without_force(tmp_path, capsys):
    (tmp_path / "pyproject.toml").write_text("[project]\nname='x'\n",
                                              encoding="utf-8")
    (tmp_path / "AGENTS.md").write_text("# keep", encoding="utf-8")
    rc = agent_cli.main(["init", str(tmp_path)])
    assert rc == 0
    out, _ = capsys.readouterr()
    assert "Skipped" in out


def test_cmd_session_ls_no_sessions(monkeypatch, tmp_path, capsys):
    monkeypatch.setattr(Path, "home", lambda: tmp_path)
    from delfin.agent import session_store as ss
    monkeypatch.setattr(
        ss, "_SESSIONS_DIR", tmp_path / ".delfin" / "agent_sessions"
    )
    rc = agent_cli.main(["session", "ls"])
    assert rc == 0
    out, _ = capsys.readouterr()
    assert "no saved sessions" in out.lower()


def test_cmd_session_search_no_query_returns_nonzero(capsys):
    """Argparse should reject missing query as a usage error (exit 2)."""
    parser = agent_cli.build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["session", "search"])
    assert exc.value.code == 2


def test_cmd_run_requires_prompt():
    """Argparse rejects `run` with no prompt — exit code 2 is the conventional
    'usage error'."""
    parser = agent_cli.build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["run"])
    assert exc.value.code == 2


def test_main_returns_unknown_subcommand_is_usage_error():
    parser = agent_cli.build_parser()
    with pytest.raises(SystemExit) as exc:
        parser.parse_args(["frobnicate"])
    assert exc.value.code == 2


# ---------------------------------------------------------------------------
# Asking for six tasks and being given one, silently
# ---------------------------------------------------------------------------

def _wanted(argv: list[str]) -> set[str]:
    """The task-id set `bench run` would actually select, mirroring the
    two lines in cmd_bench that turn the flag into a filter."""
    args = agent_cli.build_parser().parse_args(argv)
    asked = getattr(args, "task", None) or []
    if isinstance(asked, str):
        asked = [asked]
    return {t.strip() for chunk in asked
            for t in str(chunk).split(",") if t.strip()}


def test_a_repeated_task_flag_keeps_every_task():
    """`--task a --task b` was a plain single-value option, so argparse
    kept the LAST one and dropped the rest without a word. Measured live:
    six --task flags ran one task, and the run announced "on 1 tasks" as
    if that had been the request — a benchmark answering a question
    nobody asked is worse than one that refuses."""
    got = _wanted(["bench", "run", "--model", "m",
                   "--task", "office_refusal_is_accepted",
                   "--task", "office_shifted_row_is_reported",
                   "--task", "office_reconcile_uses_the_tool"])
    assert got == {"office_refusal_is_accepted",
                   "office_shifted_row_is_reported",
                   "office_reconcile_uses_the_tool"}


def test_the_comma_form_still_works():
    """The spelling that already worked must not be traded for the fix."""
    got = _wanted(["bench", "run", "--model", "m",
                   "--task", "office_refusal_is_accepted,office_shifted_row_is_reported"])
    assert got == {"office_refusal_is_accepted", "office_shifted_row_is_reported"}


def test_the_two_forms_combine():
    got = _wanted(["bench", "run", "--model", "m",
                   "--task", "a,b", "--task", "c"])
    assert got == {"a", "b", "c"}


def test_no_task_flag_still_means_every_task():
    assert _wanted(["bench", "run", "--model", "m"]) == set()


def test_a_task_id_that_matches_nothing_is_named(capsys, monkeypatch, tmp_path):
    """The second silent drop in the same place: a typo'd or renamed id was
    filtered out and the rest ran, handing back a number for a different
    question than the one asked."""
    from delfin.agent import benchmark as bm

    class _T:
        def __init__(self, tid):
            self.id = tid
    monkeypatch.setattr(bm, "load_tasks", lambda *a, **k: [_T("echt")])
    monkeypatch.setattr(agent_cli, "_bm", bm, raising=False)

    args = agent_cli.build_parser().parse_args(
        ["bench", "run", "--model", "m", "--provider", "kit",
         "--task", "echt,gibt_es_nicht"])
    # Only the selection half is exercised; the run itself needs a backend.
    monkeypatch.setattr(agent_cli, "_run_bench_suite",
                        lambda *a, **k: 0, raising=False)
    try:
        agent_cli.cmd_bench(args)
    except Exception:
        pass
    err = capsys.readouterr().err
    assert "gibt_es_nicht" in err and "skipped" in err
