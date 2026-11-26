"""One account per person: the agent account cannot change what a person runs.

The agent account (``delfin``) runs the agent. With one Linux account per
person (``scripts/pi``, ``DELFIN_PI_ACCOUNTS=1``) every person's login shell
loads the shared installation: ``agent.env``, the venv's ``activate`` and
``settings-defaults.json`` (``delfin_accounts.sh shell-block``), and every
person's DELFIN runs the code in it. Anything there the agent account could
write would be a way from the agent into every person's account.

So: the shared root and the two setting files are root's, the build is the
install account's (no login, not the agent), the deploy is root's narrow step
``delfin-deploy`` driven by the admin, the agent account's part of a deploy
reads only, every person's defaults come from root's ``/etc/delfin`` and never
from the agent account's settings, and root never reads a file a person
controls with root's rights. A real run as root on a Debian system is recorded
in the pull request; here the scripts' rules and the parts that run without
root are checked.
"""

from __future__ import annotations

import json
import os
import pathlib
import re
import subprocess

from conftest import child_env

ROOT = pathlib.Path(__file__).resolve().parents[1]
PI = ROOT / "scripts" / "pi"
BOOT = (PI / "pi_bootstrap.sh").read_text(encoding="utf-8")
DEPLOY = (PI / "delfin_deploy.sh").read_text(encoding="utf-8")
INSTALL = (PI / "pi_install.sh").read_text(encoding="utf-8")
ACCOUNTS = (PI / "delfin_accounts.sh").read_text(encoding="utf-8")
DEPLOY_PI = (PI / "deploy_pi.sh").read_text(encoding="utf-8")


def _accounts_block() -> str:
    start = BOOT.index('if [[ "$ACCOUNTS" == "1" ]]; then\n    echo "==> one account per person')
    return BOOT[start:BOOT.index("\nfi\n", start)]


def test_the_shared_root_and_what_a_login_shell_loads_are_not_the_agent_accounts():
    block = _accounts_block()
    assert 'install -d -m 750 -o root -g delfin-shared "$SHARED_ROOT"' in block
    # Nothing under the shared root is created for the agent account.
    for line in block.splitlines():
        if "$SHARED_ROOT" in line and ("install -d" in line or "chown" in line):
            assert '"$AGENT_USER"' not in line and "$AGENT_USER:" not in line, line
    assert 'chown root:root "$SHARED_ROOT/agent.env"' in block
    assert 'chmod 644 "$SHARED_ROOT/agent.env"' in block
    # The build's directories are the install account's, a system account
    # without a login that is not the agent.
    assert 'install -d -m 755 -o "$INSTALL_USER" -g "$INSTALL_USER" "$SHARED_ROOT/$d"' in block
    assert "--shell /usr/sbin/nologin" in block
    assert 'INSTALL_USER="${DELFIN_PI_INSTALL_USER:-delfin-install}"' in BOOT
    # The venv and its Python: built as the install account with root's uv.
    assert 'runuser -u "$INSTALL_USER"' in BOOT and '"$UV_SHARED"' in BOOT


def test_a_persons_login_shell_loads_only_files_of_root_and_the_install_account():
    start = ACCOUNTS.index("shell_block() {")
    func = ACCOUNTS[start:ACCOUNTS.index("\n}\n", start)]
    loaded = set(re.findall(r"'\$root/([A-Za-z0-9_./-]+)'", func))
    assert loaded == {"agent.env", ".venv/bin/activate", "settings-defaults.json"}
    # agent.env and settings-defaults.json: root writes them (deploy step),
    # the venv is the install account's (bootstrap, shared-build).
    assert 'chown root:root "$f"' in DEPLOY and 'chmod 644 "$f"' in DEPLOY
    assert 'chown root:root "$ROOT/settings-defaults.json"' in DEPLOY


def test_the_deploy_builds_as_the_install_account_and_never_runs_the_tree_as_root():
    # Every run of the deployed installer goes through one of the two
    # unprivileged accounts; root runs nothing from the tree.
    for line in DEPLOY.splitlines():
        if "pi_install.sh" in line and not line.lstrip().startswith("#"):
            assert line.lstrip().startswith(("as_installer ", "as_agent ")), line
    assert "as_installer tar -xzf - -C \"$ROOT/src.new\" < \"$tarball\"" in DEPLOY
    assert 'pi_install.sh" shared-build' in DEPLOY and 'pi_install.sh" account' in DEPLOY
    assert '[[ "$INSTALLER" != "$AGENT" ]]' in DEPLOY
    # The settings defaults are computed with the system python in isolated
    # mode, from root's files: no venv, no import of DELFIN as root.
    assert "python3 -I - /etc/delfin/pi-agent.conf /etc/delfin/pi-tabs.conf" in DEPLOY
    assert "load_settings" not in DEPLOY and "user_settings" not in DEPLOY
    # A deploy ends with the measured boundary before the agent's step.
    assert DEPLOY.index("    check_boundary\n") < DEPLOY.index('as_agent bash "$ROOT/src/scripts/pi/pi_install.sh" account')
    assert '-user "$AGENT" -print -quit' in DEPLOY


def test_the_agent_account_cannot_deploy_the_shared_installation(tmp_path):
    start = INSTALL.index('MODE="${1:-}"')
    head = INSTALL[start:INSTALL.index('VENV_PY="$ROOT/.venv/bin/python"', start)]
    home = tmp_path / "home"
    home.mkdir()
    (tmp_path / "opt").mkdir()
    (home / "delfin").symlink_to(tmp_path / "opt")
    for mode in ("install", "rollback"):
        r = subprocess.run(["bash", "-c", head + "echo REACHED", "_", mode],
                           env=dict(child_env(tmp_path), HOME=str(home)),
                           capture_output=True, text=True, timeout=30)
        assert r.returncode == 1 and "REACHED" not in r.stdout, (mode, r.stdout)
        assert "deploy_pi.sh --shared" in r.stderr
    # pi_install.sh no longer writes the persons' defaults at all.
    assert "settings-defaults.json" not in INSTALL.replace(
        "delfin-deploy writes $ROOT/settings-defaults.json", "")
    assert "shared-build)" in INSTALL and "is the install account's" in INSTALL


def test_every_persons_defaults_come_from_roots_files_not_the_agents_settings(tmp_path):
    start = DEPLOY.index("root_settings_defaults() {")
    py = DEPLOY[DEPLOY.index("<<'EOF'\n", start) + len("<<'EOF'\n"):DEPLOY.index("\nEOF\n", start)]
    agent_conf = tmp_path / "pi-agent.conf"
    agent_conf.write_text("# box defaults\nmode=solo\nprovider=kit\npermission_profile=ask_all\neffort=\n")
    tabs_conf = tmp_path / "pi-tabs.conf"
    tabs_conf.write_text("# shown\nagent home\n")
    dashboard = tmp_path / "dashboard.py"
    dashboard.write_text("TABS = [{'id': 'home'}, {'id': 'agent'}, {'id': 'archive'}, {'id': 'settings'}]\n")
    out = tmp_path / "settings-defaults.json"
    # The agent account's own settings say Bypass: they must not reach a person.
    agent_home = tmp_path / "agent-home"
    agent_home.mkdir()
    (agent_home / ".delfin_settings.json").write_text(json.dumps({"agent": {"permission_profile": "all_free"}}))
    r = subprocess.run(["python3", "-I", "-", str(agent_conf), str(tabs_conf), str(dashboard), str(out)],
                       input=py, capture_output=True, text=True, timeout=30,
                       env=dict(child_env(tmp_path), HOME=str(agent_home)))
    assert r.returncode == 0, r.stderr
    data = json.loads(out.read_text())
    assert data["agent"] == {"mode": "solo", "provider": "kit", "permission_profile": "ask_all", "effort": None}
    assert data["ui"]["tabs"] == {"order": ["agent", "home", "archive"], "hidden": ["archive"]}
    assert oct(out.stat().st_mode & 0o777) == "0o644"
    # A link planted at the target is replaced, never written through.
    victim = tmp_path / "victim"
    victim.write_text("unchanged")
    out.unlink()
    (tmp_path / "settings-defaults.json.tmp").symlink_to(victim)
    r = subprocess.run(["python3", "-I", "-", str(agent_conf), str(tabs_conf), str(dashboard), str(out)],
                       input=py, capture_output=True, text=True, timeout=30, env=child_env(tmp_path))
    assert victim.read_text() == "unchanged"


def test_root_never_reads_a_persons_file_with_roots_rights():
    start = ACCOUNTS.index("invite() {")
    func = ACCOUNTS[start:ACCOUNTS.index("\n}\n", start)]
    # The package is root's copy, not the shared tree's.
    assert 'package="$LAUNCHER_ZIP"' in func and "$root/src" not in func
    assert 'LAUNCHER_ZIP="${DELFIN_ACCOUNTS_LAUNCHER_ZIP:-/usr/local/lib/delfin/delfin-windows-launcher.zip}"' in ACCOUNTS
    # The invitation in the person's home is read and removed as the person.
    assert 'run_as "$name" cat -- "$home/invites/DELFIN-$name.zip" > "$tmp"' in func
    assert 'run_as "$name" rm -f -- "$home/invites/DELFIN-$name.zip"' in func
    assert 'install -m 600 -o "$admin"' not in func
    # list reads each person's authorized_keys as that person.
    start = ACCOUNTS.index("list() {")
    func = ACCOUNTS[start:ACCOUNTS.index("\n}\n", start)]
    for line in func.splitlines():
        if "authorized_keys" in line:
            assert 'run_as "$u"' in line, line
    # add: the agent account's keep default is tested as that account.
    assert 'run_as "$deployer" test -f "$dhome/.config/delfin/launcher-keep.conf"' in ACCOUNTS


def test_the_admin_deploys_the_shared_installation_through_root():
    assert "--shared)          action=shared" in DEPLOY_PI
    assert "/usr/local/sbin/delfin-deploy install '$version' '$remote/src.tgz'" in DEPLOY_PI
    # The upload lands in a mktemp directory of the admin's, not a fixed name.
    assert "d=\"$(mktemp -d)\" && cat > \"$d/src.tgz\"" in DEPLOY_PI
    # The installer staged for one-account deploys no longer lives in
    # ~/delfin (root's shared root in the variant).
    assert "~/delfin/pi_install.sh" not in DEPLOY_PI
    assert "install -D -m 755 -o root -g root \"$HERE/delfin_deploy.sh\" /usr/local/sbin/delfin-deploy" in BOOT
