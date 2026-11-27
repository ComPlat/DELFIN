"""A Pi that already runs the one-account layout moves to the per-person variant.

``pi_bootstrap.sh`` with ``DELFIN_PI_ACCOUNTS=1`` sets the shared root up on
a fresh box only: where ``~delfin/delfin`` is a directory it stops. With
``DELFIN_PI_MIGRATE=1`` it first runs ``scripts/pi/pi_accounts_migrate.sh``,
which creates the shared root closed (0700) before anything moves, moves the
old installation aside (kept whole), keeps the agent's workspace in its own
home (``~/delfin-work``), makes ``agent.env`` root's copy in the shared root,
switches the launcher's people mode in the agent account off and warns,
loudly, about every key and every person's home still in that account.
``DELFIN_PI_MIGRATE=back`` reverses each step.

The migration script needs no root (renames inside the agent's home and the
shared root), so it runs for real here on a tree laid out like the lab Pi's.
A full run as root -- bootstrap, deploy, services, a person's login -- on a
Debian system that ran the one-account layout is recorded in the pull request.
"""

from __future__ import annotations

import os
import pathlib
import subprocess

from conftest import child_env

ROOT = pathlib.Path(__file__).resolve().parents[1]
PI = ROOT / "scripts" / "pi"
MIGRATE = PI / "pi_accounts_migrate.sh"


def _box(tmp_path: pathlib.Path):
    """An agent home as on the lab Pi (2026-11-26): installation, workspace,
    people mode with two persons, a launcher key and an open invitation."""
    home = tmp_path / "home" / "delfin"
    inst = home / "delfin"
    for rel, text in {
        "src/delfin/__init__.py": "__version__ = '1.4.0'\n",
        ".venv/bin/python": "#!/bin/sh\n",
        "src.prev/README": "previous deploy\n",
        "work/session-1/notes.md": "agent work\n",
        "agent.env": "export DELFIN_AGENT_SANDBOX=0\n",
        "DEPLOYED_VERSION": "a54b758c6\n",
    }.items():
        path = inst / rel
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(text)
    (home / ".config" / "delfin").mkdir(parents=True)
    (home / ".config" / "delfin" / "launcher-people.conf").write_text("# people mode\n")
    for person in ("hmaxi", "sx4535"):
        (home / "people" / person / ".delfin").mkdir(parents=True)
    (home / ".ssh").mkdir()
    (home / ".ssh" / "authorized_keys").write_text(
        "ssh-ed25519 AAAAdeploy delfin-pi-deploy\n"
        "restrict,port-forwarding,pty ssh-ed25519 AAAAsilvia delfin-launcher-silvia\n"
        'restrict,command="/x --enroll nicole",expiry-time="202611280900Z" ssh-ed25519 AAAAinv delfin-invite-nicole\n'
    )
    shared = tmp_path / "opt" / "delfin"
    conf = tmp_path / "etc" / "accounts.conf"
    conf.parent.mkdir()
    return home, shared, conf


def _run(direction, home, shared, conf):
    tmp = shared.parent.parent
    env = child_env(tmp)
    env["DELFIN_ACCOUNTS_CONF"] = str(conf)
    return subprocess.run(
        ["bash", str(MIGRATE), direction, str(home), str(shared)],
        capture_output=True, text=True, env=env, timeout=60,
    )


def _tree(base: pathlib.Path):
    out = {}
    for path in sorted(base.rglob("*")):
        if path.is_file() and not path.is_symlink():
            out[str(path.relative_to(base))] = path.read_text()
    return out


def test_the_old_installation_is_moved_aside_whole_and_the_workspace_stays_the_agents(tmp_path):
    home, shared, conf = _box(tmp_path)
    before = _tree(home / "delfin")
    result = _run("to", home, shared, conf)
    assert result.returncode == 0, result.stderr
    kept = home / "delfin.one-account"
    assert not (home / "delfin").exists()
    # Everything but the workspace stays where it was, unchanged.
    after = _tree(kept)
    assert {k: v for k, v in before.items() if not k.startswith("work/")} == {
        k: v for k, v in after.items() if k != ".moved-to-shared-root"
    }
    # The workspace is the agent's data: in its home, not in the shared root.
    assert (home / "delfin-work" / "session-1" / "notes.md").read_text() == "agent work\n"
    assert not (shared / "work").exists()
    # agent.env: a copy of the shared root's own, not a link and not shared mode.
    assert (shared / "agent.env").read_text() == "export DELFIN_AGENT_SANDBOX=0\n"
    assert not (shared / "agent.env").is_symlink()
    assert oct((shared / "agent.env").stat().st_mode & 0o777) == "0o644"
    assert (kept / "agent.env").exists()
    # The shared root was created closed; the bootstrap opens it afterwards.
    assert oct(shared.stat().st_mode & 0o777) == "0o700"


def test_people_mode_is_switched_off_and_every_old_way_in_is_named_not_removed(tmp_path):
    home, shared, conf = _box(tmp_path)
    keys_before = (home / ".ssh" / "authorized_keys").read_text()
    result = _run("to", home, shared, conf)
    assert result.returncode == 0, result.stderr
    cfg = home / ".config" / "delfin"
    assert not (cfg / "launcher-people.conf").exists()
    assert (cfg / "launcher-people.conf.off-accounts").read_text() == "# people mode\n"
    out = result.stdout
    assert "WARNING: the agent account still holds the persons of the old layout." in out
    assert "delfin-launcher-silvia" in out and "delfin-invite-nicole" in out
    assert "delfin-pi-deploy" not in out
    assert "sudo delfin-accounts add NAME" in out
    assert str(home / "people" / "hmaxi") in out and str(home / "people" / "sx4535") in out
    # Access is the admin's decision: the file is exactly as it was.
    assert (home / ".ssh" / "authorized_keys").read_text() == keys_before
    assert (home / "people" / "sx4535" / ".delfin").is_dir()


def test_moving_again_after_the_link_exists_changes_nothing(tmp_path):
    home, shared, conf = _box(tmp_path)
    assert _run("to", home, shared, conf).returncode == 0
    (home / "delfin").symlink_to(shared)          # what the bootstrap does next
    snapshot = _tree(tmp_path)
    again = _run("to", home, shared, conf)
    assert again.returncode == 0
    assert "already a link" in again.stdout
    assert _tree(tmp_path) == snapshot


def test_a_move_that_would_overwrite_something_is_refused_before_anything(tmp_path):
    home, shared, conf = _box(tmp_path)
    (home / "delfin.one-account").mkdir()
    snapshot = _tree(tmp_path)
    result = _run("to", home, shared, conf)
    assert result.returncode == 1 and "Nothing was changed" in result.stderr
    assert _tree(tmp_path) == snapshot

    home2, shared2, conf2 = _box(tmp_path / "b")
    shared2.mkdir(parents=True)
    (shared2 / "something").write_text("x")
    snapshot = _tree(tmp_path / "b")
    result = _run("to", home2, shared2, conf2)
    assert result.returncode == 1 and "not an empty directory" in result.stderr
    assert _tree(tmp_path / "b") == snapshot

    # A link planted where the shared root goes is refused too.
    home3, shared3, conf3 = _box(tmp_path / "c")
    shared3.parent.mkdir(parents=True)
    elsewhere = tmp_path / "c" / "elsewhere"
    elsewhere.mkdir()
    shared3.symlink_to(elsewhere)
    result = _run("to", home3, shared3, conf3)
    assert result.returncode == 1 and (home3 / "delfin").is_dir()


def test_an_agent_env_link_does_not_carry_another_file_into_the_shared_root(tmp_path):
    home, shared, conf = _box(tmp_path)
    secret = tmp_path / "secret"
    secret.write_text("not for the persons\n")
    env_file = home / "delfin" / "agent.env"
    env_file.unlink()
    env_file.symlink_to(secret)
    # As root the copy runs with the agent account's rights; here, as the test
    # user, the rule that matters is that the shared copy is a plain file the
    # script wrote, never the link itself.
    result = _run("to", home, shared, conf)
    assert result.returncode == 0, result.stderr
    assert not (shared / "agent.env").is_symlink()
    script = MIGRATE.read_text(encoding="utf-8")
    assert 'as_owner cat -- "$kept/agent.env" > "$shared/agent.env.tmp"' in script
    assert 'install -d -m 700 "$shared"' in script and "mkdir -p \"$shared\"" not in script
    assert script.index('install -d -m 700 "$shared"') < script.index('as_owner mv "$current" "$kept"')


def test_back_returns_the_box_exactly_as_it_was(tmp_path):
    home, shared, conf = _box(tmp_path)
    before = _tree(home)
    assert _run("to", home, shared, conf).returncode == 0
    (home / "delfin").symlink_to(shared)
    (shared / "work").symlink_to(home / "delfin-work")   # the bootstrap's link
    conf.write_text("DELFIN_SHARED_ROOT=/opt/delfin\n")
    result = _run("back", home, shared, conf)
    assert result.returncode == 0, result.stderr
    assert (home / "delfin").is_dir() and not (home / "delfin").is_symlink()
    assert not (home / "delfin.one-account").exists()
    assert not (home / "delfin-work").exists()
    assert _tree(home) == before
    assert not conf.exists()
    assert "deploy_pi.sh" in result.stdout


def test_back_refuses_a_link_to_somewhere_else_or_a_real_directory(tmp_path):
    home, shared, conf = _box(tmp_path)
    assert _run("to", home, shared, conf).returncode == 0
    elsewhere = tmp_path / "elsewhere"
    elsewhere.mkdir()
    (home / "delfin").symlink_to(elsewhere)
    snapshot = _tree(tmp_path)
    result = _run("back", home, shared, conf)
    assert result.returncode == 1 and "Nothing was changed" in result.stderr
    assert _tree(tmp_path) == snapshot

    (home / "delfin").unlink()
    (home / "delfin").mkdir()
    result = _run("back", home, shared, conf)
    assert result.returncode == 1 and "is a directory" in result.stderr


def test_the_bootstrap_moves_only_when_asked_and_goes_back_before_installing_anything():
    boot = (PI / "pi_bootstrap.sh").read_text(encoding="utf-8")
    deploy = (PI / "deploy_pi.sh").read_text(encoding="utf-8")
    # The upload carries the script, and the bootstrap insists on it.
    assert deploy.count("pi_accounts_migrate.sh") >= 2
    assert "delfin_accounts.sh delfin_people.sh delfin_deploy.sh pi_accounts_migrate.sh" in boot
    # Without DELFIN_PI_MIGRATE=1 an existing directory still stops the run.
    assert 'if [[ "$MIGRATE" != 1 ]]; then' in boot
    assert 'pi_accounts_migrate.sh" to "$AGENT_HOME" "$SHARED_ROOT"' in boot
    # The way back runs before apt, the user and the venv are touched.
    assert boot.index('pi_accounts_migrate.sh" back') < boot.index("apt-get update")
    # A move to the variant without the variant is refused.
    assert "set DELFIN_PI_ACCOUNTS=1 too" in boot


def test_the_migration_script_is_in_the_bootstrap_folder_and_executable_text():
    assert MIGRATE.is_file()
    assert MIGRATE.read_bytes().startswith(b"#!/usr/bin/env bash\n")
    assert b"\r\n" not in MIGRATE.read_bytes()
    assert os.access(MIGRATE, os.R_OK)
