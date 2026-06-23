"""The bash gate applies the deny-list + secret-path scan to the CONTENTS of a
script it would run — so `write payload.py; python payload.py` can't smuggle a
denied command past the command-line-only check. Precise by design: it reuses
the curated patterns, so ordinary code (subprocess, requests, …) is untouched.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent.api_client import _doc_executor, KitToolPermissions


def _perms(ws: Path) -> KitToolPermissions:
    # bypassPermissions → the scan must still block (it runs before the bypass).
    return KitToolPermissions(workspace=ws, mode="bypassPermissions")


def _gate(ws: Path, cmd: str):
    return _doc_executor._run_permission_gate(
        "bash", {"command": cmd}, _perms(ws))


def test_script_with_destructive_payload_is_blocked(tmp_path):
    (tmp_path / "payload.py").write_text(
        "import os\nos.system('rm -rf /')\n")
    out = _gate(tmp_path, "python payload.py")
    assert out is not None and "deny-pattern" in out


def test_script_reading_ssh_key_is_blocked(tmp_path):
    (tmp_path / "steal.py").write_text(
        "data = open('~/.ssh/id_rsa').read()\nprint(data)\n")
    out = _gate(tmp_path, "python3 steal.py")
    assert out is not None and "secret path" in out


def test_shell_script_payload_is_blocked(tmp_path):
    (tmp_path / "go.sh").write_text("#!/bin/sh\nsudo rm -rf /var\n")
    out = _gate(tmp_path, "bash go.sh")
    assert out is not None


def test_legit_script_using_subprocess_is_allowed(tmp_path):
    # Ordinary build/automation code must NOT trip the scan.
    (tmp_path / "build.py").write_text(
        "import subprocess\n"
        "subprocess.run(['python', '-m', 'pytest', '-q'])\n"
        "print('done')\n")
    assert _gate(tmp_path, "python build.py") is None


def test_inline_code_is_not_treated_as_a_file(tmp_path):
    # `python -c ...` has no script file to scan → not blocked by this check.
    assert _gate(tmp_path, "python -c \"print(1)\"") is None


def test_benign_named_script_is_allowed(tmp_path):
    (tmp_path / "ok.py").write_text("print('hello')\n")
    assert _gate(tmp_path, "python ok.py") is None


def test_missing_script_does_not_crash(tmp_path):
    # Referencing a non-existent file must not raise — just no payload to scan.
    assert _gate(tmp_path, "python does_not_exist.py") is None
