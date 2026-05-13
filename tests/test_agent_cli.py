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
