"""`delfin-agent audit last` shows the newest records through the one reader.

Before: ``audit_log.read_last_n`` existed, was tested directly, and
nothing in the product called it — two consolidations (this file and
#849) rewired it through ``_scan`` as "the one reader's public face",
but the operator who wanted to SEE the newest audit records had no
command. Mounted as ``delfin-agent audit last [N]``: the records print
as one-liners, and the integrity verdict on the SAME bytes is printed
FIRST, so a tampered log is named before its records are shown.
"""
from __future__ import annotations

import pytest

from delfin.agent import audit_log, cli


@pytest.fixture
def log_home(tmp_path, monkeypatch):
    monkeypatch.setenv("HOME", str(tmp_path))
    return tmp_path


def _write_log(path, n=3, sid="s-audit"):
    for i in range(n):
        audit_log.append(
            audit_log.make_record(tool="write_file", decision="ok",
                                  path=f"file{i}.py", session_id=sid),
            log_path=path)
    return path


def test_audit_last_prints_the_newest_records(log_home, capsys):
    path = log_home / "audit.log"
    _write_log(path, 3)
    code = cli.main(["audit", "last", "2", "--log", str(path)])
    assert code == 0
    out = capsys.readouterr().out
    assert "file2.py" in out          # newest is included
    assert "file0.py" not in out      # older than the count is not
    assert "NOT intact" not in out    # an intact log warns about nothing


def test_audit_last_names_a_tampered_log_before_its_records(
        log_home, capsys):
    """An edited line means the verdict comes first."""
    path = log_home / "audit.log"
    _write_log(path, 3)
    lines = path.read_text().splitlines()
    lines[1] = lines[1].replace("file1.py", "fileX.py")
    path.write_text("\n".join(lines) + "\n")
    code = cli.main(["audit", "last", "10", "--log", str(path)])
    assert code == 0
    out = capsys.readouterr().out
    warn_at = out.find("NOT intact")
    assert warn_at != -1, out
    assert warn_at < out.find("file"), "the warning must come first"


def test_audit_last_on_an_empty_log_is_a_clear_line(log_home, capsys):
    code = cli.main(["audit", "last", "5", "--log",
                     str(log_home / "audit.log")])
    assert code == 0
    assert "no audit records found" in capsys.readouterr().out


def test_audit_verify_still_works(log_home, capsys):
    """The existing subcommand is untouched by the new one."""
    path = log_home / "audit.log"
    _write_log(path, 2)
    code = cli.main(["audit", "verify", "--log", str(path)])
    assert code == 0
    assert "intact" in capsys.readouterr().out
