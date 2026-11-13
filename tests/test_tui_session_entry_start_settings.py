"""Entering a session from the terminal as ``claude`` does:

    delfin-agent tui                a new session (unchanged)
    delfin-agent tui --continue     the folder's latest session
    delfin-agent tui --resume <id>  that session
    delfin-agent tui --resume       the resume picker first

The launcher hands the flags to the window in DELFIN_TUI_START
(``cli_service._start_settings``); a bare ``--resume`` arrives as
``RESUME_PICK`` and becomes ``pick`` -- never a session id the service
would try to resume.
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pytest

from delfin.agent.service import cli_service

REPO = Path(__file__).resolve().parents[1]


def _args(tmp_path: Path, **given) -> argparse.Namespace:
    return argparse.Namespace(cwd=str(tmp_path), **given)


def test_no_flag_starts_a_new_session_as_before(tmp_path):
    settings = cli_service._start_settings(_args(tmp_path))
    assert settings == {"cwd": str(tmp_path)}


def test_keep_starts_the_windows_sessions_kept(tmp_path):
    """`delfin-agent tui --keep`: the sessions this window starts keep
    running after it closes (Keep session on). Without it they end with
    the window (ctrl+c twice); /keep on keeps one."""
    assert cli_service._start_settings(_args(tmp_path, keep=True)) == {
        "cwd": str(tmp_path), "keep": True}
    from delfin.agent import cli
    assert cli.build_parser().parse_args(["tui", "--keep"]).keep is True
    assert cli.build_parser().parse_args(["tui"]).keep is False


def test_continue_and_resume_with_an_id_are_handed_over(tmp_path):
    assert cli_service._start_settings(_args(tmp_path, continue_session=True)) == {
        "cwd": str(tmp_path), "continue": True}
    assert cli_service._start_settings(_args(tmp_path, session="0f3a9c21")) == {
        "cwd": str(tmp_path), "resume": "0f3a9c21"}


def test_a_bare_resume_opens_the_picker_and_resumes_nothing(tmp_path):
    settings = cli_service._start_settings(_args(tmp_path, session=cli_service.RESUME_PICK))
    assert settings == {"cwd": str(tmp_path), "pick": True}


def test_the_marker_is_the_window_s_own():
    entry_ts = REPO / "tui" / "src" / "entry.ts"
    if not (REPO / "tui").is_dir():
        # A copy of the Python side only (test-on-slurm): CI has the window.
        pytest.skip("no tui/ in this copy of the repository")
    source = entry_ts.read_text(encoding="utf-8")
    m = re.search(r'export const RESUME_PICK = "([^"]+)"', source)
    assert m, "entry.ts names RESUME_PICK"
    assert m.group(1) == cli_service.RESUME_PICK


def test_the_parser_takes_continue_and_resume_with_an_id():
    from delfin.agent.cli import build_parser

    args = build_parser().parse_args(["tui", "--continue"])
    assert args.continue_session is True and not args.session
    args = build_parser().parse_args(["tui", "--resume", "0f3a9c21"])
    assert args.session == "0f3a9c21"
    args = build_parser().parse_args(["tui"])
    assert not args.session and not args.continue_session


@pytest.mark.xfail(reason="needs the protected cli.py change: --resume nargs='?' "
                          "const=cli_service.RESUME_PICK (see the OPERATOR_A note)",
                   strict=False)
def test_the_parser_takes_a_bare_resume():
    from delfin.agent.cli import build_parser

    try:
        args = build_parser().parse_args(["tui", "--resume"])
    except SystemExit:
        pytest.fail("a bare --resume is refused by the parser")
    assert args.session == cli_service.RESUME_PICK
