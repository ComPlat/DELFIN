"""One Linux account per person on the Pi, as an opt-in beside people mode.

People mode (the default) puts everyone into the agent user's account under
``~/people/<name>``: one Unix user, an overview and no boundary. The accounts
variant (``pi_bootstrap.sh`` with ``DELFIN_PI_ACCOUNTS=1``) shares the
installation read only from ``/opt/delfin`` and gives each person an account
of their own through ``delfin-accounts`` (``scripts/pi/delfin_accounts.sh``):
created without a password, let through sshd's allow-list, reached by an
invitation that leads into that account, locked or removed on its own.

The account operations need root (useradd, sshd); here the real script runs
with those few system commands replaced by stubs that keep a passwd, group
and shadow file of their own, so its own logic -- the allow-list it writes,
what lock and unlock change, which accounts it refuses -- runs for real. The
launcher half runs for real too: the bootstrap the Windows app starts finds
the shared installation from a person's login shell. A full run as root on a
Debian system (bootstrap, deploy, add, invitation, login, lock, remove) is
recorded in the pull request.
"""

from __future__ import annotations

import json
import os
import pathlib
import pwd
import shlex
import subprocess
import sys
import zipfile

from conftest import child_env

ROOT = pathlib.Path(__file__).resolve().parents[1]
PI = ROOT / "scripts" / "pi"
ACCOUNTS = PI / "delfin_accounts.sh"
PEOPLE = PI / "delfin_people.sh"
LAUNCHER = ROOT / "tools" / "windows-launcher"

# --------------------------------------------------------------------------
# A small user database behind stub commands. Formats are the system's:
# passwd name:x:uid:gid:gecos:home:shell, group name:x:gid:members,
# shadow name:pw:...:expire(8th field):.
STUBS = {
    "getent": r'''#!/usr/bin/env bash
db="$1"; key="${2:-}"
[[ -n "$key" ]] || { cat "$FAKE/$db"; exit 0; }
grep -m1 "^$key:" "$FAKE/$db" || exit 2
''',
    "id": r'''#!/usr/bin/env bash
case "$1" in
  -un) echo "${USER:-root}" ;;
  -gn) /usr/bin/id -gn "${2:-}" ;;
  -*) /usr/bin/id "$@" ;;
  *) grep -q "^$1:" "$FAKE/passwd" ;;
esac
''',
    "groupadd": r'''#!/usr/bin/env bash
echo "${@: -1}:x:900:" >> "$FAKE/group"
''',
    "useradd": r'''#!/usr/bin/env bash
name="${@: -1}"; groups=""; echo "useradd $*" >> "$FAKE/log"
while [[ $# -gt 1 ]]; do [[ "$1" == --groups ]] && groups="$2"; shift; done
home="$FAKE/home/$name"; mkdir -p "$home"
printf '# skel\n' > "$home/.bashrc"
echo "$name:x:1001:1001:DELFIN person:$home:/bin/bash" >> "$FAKE/passwd"
echo "$name:x:1001:" >> "$FAKE/group"
echo "$name:!:20000:0:99999:7:::" >> "$FAKE/shadow"
for g in ${groups//,/ }; do
  sed -i -E "s/^($g:x:[0-9]+:)(.*)$/\1\2,$name/; s/^($g:x:[0-9]+:),/\1/" "$FAKE/group"
done
''',
    "usermod": r'''#!/usr/bin/env bash
if [[ "$1" == -aG ]]; then
  grep -q "^$2:.*[:,]$3\(,\|$\)" "$FAKE/group" && exit 0
  sed -i -E "s/^($2:x:[0-9]+:)(.*)$/\1\2,$3/; s/^($2:x:[0-9]+:),/\1/" "$FAKE/group"
  echo "usermod $*" >> "$FAKE/log"; exit 0
fi
[[ "$1" == --expiredate ]] || exit 1
awk -F: -v OFS=: -v u="$3" -v e="$2" '$1==u{$8=e}1' "$FAKE/shadow" > "$FAKE/shadow.new"
mv "$FAKE/shadow.new" "$FAKE/shadow"; echo "usermod $*" >> "$FAKE/log"
''',
    "chage": r'''#!/usr/bin/env bash
awk -F: -v OFS=: -v u="$3" '$1==u{$8=""}1' "$FAKE/shadow" > "$FAKE/shadow.new"
mv "$FAKE/shadow.new" "$FAKE/shadow"; echo "chage $*" >> "$FAKE/log"
''',
    # As the person: the clean environment run_as builds, with the stubs
    # still on PATH (id -un answers USER) and the fake host keys.
    "runuser": r'''#!/usr/bin/env bash
shift 3  # -u NAME --
[[ "$1" == env && "$2" == -i ]] && shift 2
args=()
for a in "$@"; do [[ "$a" == PATH=* ]] || args+=("$a"); done
exec env -i "PATH=$STUB_BIN:/usr/bin:/bin" "DELFIN_HOST_KEYS=$DELFIN_HOST_KEYS" "${args[@]}"
''',
    "sshd": r'''#!/usr/bin/env bash
echo "sshd $*" >> "$FAKE/log"
[[ ! -f "$FAKE/sshd-refuses" ]]
''',
    "systemctl": '#!/usr/bin/env bash\necho "systemctl $*" >> "$FAKE/log"\n',
    "loginctl": '#!/usr/bin/env bash\necho "loginctl $*" >> "$FAKE/log"\n',
    "pkill": '#!/usr/bin/env bash\necho "pkill $*" >> "$FAKE/log"\n',
    "pgrep": '#!/usr/bin/env bash\necho 0\n',
}


def _box(tmp_path: pathlib.Path):
    """A fake box: user database, hardened sshd file, accounts.conf, the
    root-owned delfin-people, and delfin_accounts.sh pointed at them."""
    fake = tmp_path / "fake"
    (fake / "home").mkdir(parents=True)
    me = pwd.getpwuid(os.getuid()).pw_name
    admin_home = fake / "home" / "pi"
    admin_home.mkdir()
    (fake / "passwd").write_text(
        f"pi:x:1000:1000::{admin_home}:/bin/bash\n"
        f"delfin:x:1002:1002::{fake / 'home' / 'delfin'}:/bin/bash\n"
        f"{me}:x:1003:1003::{admin_home}:/bin/bash\n")
    (fake / "group").write_text("pi:x:1000:\ndelfin:x:1002:\ndelfin-shared:x:901:delfin\n")
    (fake / "shadow").write_text("pi:*:20000:0:99999:7:::\ndelfin:*:20000:0:99999:7:::\n")
    (fake / "log").write_text("")
    sshd_conf = fake / "40-delfin-hardening.conf"
    sshd_conf.write_text("# DELFIN agent box\nMatch all\nPasswordAuthentication no\n"
                         "MaxAuthTries 3\nAllowUsers pi delfin\n")
    shared = tmp_path / "opt" / "delfin"
    (shared / "src" / "tools").mkdir(parents=True)
    # Invitations are built from root's copy (delfin-deploy), never the tree.
    lib = fake / "usr-local-lib-delfin"
    lib.mkdir()
    with zipfile.ZipFile(lib / "delfin-windows-launcher.zip", "w") as z:
        z.writestr("windows-launcher/Install.cmd", "@echo off\r\n")
    conf = fake / "accounts.conf"
    conf.write_text(f"DELFIN_SHARED_ROOT={shared}\nDELFIN_DEPLOY_USER=delfin\n")
    host_keys = fake / "etc-ssh"
    host_keys.mkdir()
    subprocess.run(["ssh-keygen", "-q", "-t", "ed25519", "-N", "", "-f",
                    str(host_keys / "ssh_host_ed25519_key")], check=True)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in STUBS.items():
        (bin_dir / name).write_text(body)
        (bin_dir / name).chmod(0o755)
    # The script, pointed at the fake box. Only the constants and the root
    # check change; everything else runs as written.
    text = ACCOUNTS.read_text(encoding="utf-8")
    for old, new in (
        ('SSHD_CONF="/etc/ssh/sshd_config.d/40-delfin-hardening.conf"', f'SSHD_CONF="{sshd_conf}"'),
        ('PEOPLE_TOOL="/usr/local/lib/delfin/delfin-people"', f'PEOPLE_TOOL="{bin_dir / "delfin-people"}"'),
        ('ARCHIVE_DIR="/var/backups/delfin-people"', f'ARCHIVE_DIR="{fake / "archive"}"'),
        ('[[ $EUID -eq 0 ]] || die "run with sudo"', ': stubbed-root-check'),
        # Files here belong to the test user, not root.
        ('[[ "$(stat -c %U "$package")" == root ]]', '[[ -n stubbed-owner-check ]]'),
    ):
        assert text.count(old) == 1, old
        text = text.replace(old, new)
    script = bin_dir / "delfin-accounts"
    script.write_text(text)
    script.chmod(0o755)
    (bin_dir / "delfin-people").write_text(PEOPLE.read_text(encoding="utf-8"))
    (bin_dir / "delfin-people").chmod(0o755)
    env = dict(os.environ, FAKE=str(fake), PATH=f"{bin_dir}:/usr/bin:/bin", STUB_BIN=str(bin_dir),
               DELFIN_ACCOUNTS_CONF=str(conf), SUDO_USER=me,
               DELFIN_ACCOUNTS_LAUNCHER_ZIP=str(lib / "delfin-windows-launcher.zip"),
               DELFIN_HOST_KEYS=str(host_keys / "ssh_host_*_key.pub"))
    return fake, script, env, sshd_conf, shared


def _run(script, env, *args):
    return subprocess.run(["bash", str(script), *args], env=env, capture_output=True,
                          text=True, timeout=60)


def _allow(sshd_conf):
    return [line for line in sshd_conf.read_text().splitlines() if line.startswith("AllowUsers")]


# --------------------------------------------------------------------------
def test_adding_a_person_makes_an_account_without_password_and_lets_it_through_sshd(tmp_path):
    fake, script, env, sshd_conf, shared = _box(tmp_path)
    r = _run(script, env, "add", "anna", "--no-invite")
    assert r.returncode == 0, r.stderr + r.stdout
    log = (fake / "log").read_text()
    assert "useradd --create-home --shell /bin/bash --user-group --groups delfin-people" in log
    assert "delfin-people:x:900:anna" in (fake / "group").read_text()
    # No password ever set: the shadow entry stays locked ('!'), keys only.
    assert "anna:!:" in (fake / "shadow").read_text()
    assert _allow(sshd_conf) == ["AllowUsers pi delfin anna"]
    assert "sshd -t" in log and "systemctl reload ssh" in log
    home = fake / "home" / "anna"
    assert oct(home.stat().st_mode & 0o777) == "0o700"
    bashrc = (home / ".bashrc").read_text()
    assert bashrc.startswith("# skel\n"), "the account's own .bashrc is kept, the block appended"
    assert f". '{shared}/.venv/bin/activate'" in bashrc
    assert oct((home / ".ssh" / "authorized_keys").stat().st_mode & 0o777) == "0o600"


def test_adding_again_or_an_account_that_is_not_a_person_changes_nothing(tmp_path):
    fake, script, env, sshd_conf, _ = _box(tmp_path)
    assert _run(script, env, "add", "anna", "--no-invite").returncode == 0
    before = (fake / "passwd").read_text()
    # Run again: the setup is checked, nothing is replaced, no new account.
    again = _run(script, env, "add", "anna")
    assert again.returncode == 0 and "already a person" in again.stdout, again.stderr
    assert not (fake / "home" / "pi" / "invites").exists(), "no new invitation on a re-run"
    for other in ("pi", "delfin"):
        r = _run(script, env, "add", other, "--no-invite")
        assert r.returncode != 0 and "not a DELFIN person" in r.stderr
        for verb in ("lock", "remove", "invite"):
            r = _run(script, env, verb, other)
            assert r.returncode != 0 and "not a DELFIN person" in r.stderr, (verb, other)
    assert (fake / "passwd").read_text() == before
    assert _allow(sshd_conf) == ["AllowUsers pi delfin anna"]


def test_names_that_are_not_plain_user_names_are_refused_before_anything(tmp_path):
    fake, script, env, _, _ = _box(tmp_path)
    for bad in ("-rf", "Anna", "1anna", "an na", "anna/..", "a" * 33, ""):
        r = _run(script, env, "add", bad)
        assert r.returncode != 0, bad
    assert "useradd" not in (fake / "log").read_text()


def test_locking_one_person_takes_only_that_one_out_and_unlocking_brings_it_back(tmp_path):
    fake, script, env, sshd_conf, _ = _box(tmp_path)
    for name in ("anna", "ben"):
        assert _run(script, env, "add", name, "--no-invite").returncode == 0
    assert _allow(sshd_conf) == ["AllowUsers pi delfin anna ben"]
    r = _run(script, env, "lock", "anna")
    assert r.returncode == 0, r.stderr
    log = (fake / "log").read_text()
    # Expired, not just a locked password: sshd's PAM account check refuses
    # an expired account even with a valid key.
    assert "usermod --expiredate 1 anna" in log
    assert "loginctl terminate-user anna" in log and "pkill -KILL -u anna" in log
    assert _allow(sshd_conf) == ["AllowUsers pi delfin ben"]
    listing = _run(script, env, "list").stdout
    assert any(line.split()[:2] == ["anna", "locked"] for line in listing.splitlines())
    assert any(line.split()[:2] == ["ben", "active"] for line in listing.splitlines())
    assert _run(script, env, "invite", "anna").returncode != 0
    assert _run(script, env, "unlock", "anna").returncode == 0
    assert _allow(sshd_conf) == ["AllowUsers pi delfin anna ben"]


def test_a_removed_person_leaves_the_allow_list(tmp_path):
    # Found in the real run: after userdel the name was no member any more
    # and stayed in AllowUsers as if it were an account like pi.
    fake, script, env, sshd_conf, _ = _box(tmp_path)
    for name in ("anna", "ben"):
        assert _run(script, env, "add", name, "--no-invite").returncode == 0
    # userdel as the system does it: the account and its membership go.
    (fake / "passwd").write_text("".join(l for l in (fake / "passwd").read_text().splitlines(True)
                                         if not l.startswith("anna:")))
    (fake / "group").write_text((fake / "group").read_text().replace(":anna,ben", ":ben"))
    r = _run(script, env, "unlock", "ben")  # any verb that re-syncs the list
    assert r.returncode == 0, r.stderr
    assert _allow(sshd_conf) == ["AllowUsers pi delfin ben"]


def test_an_allow_list_sshd_refuses_is_put_back(tmp_path):
    fake, script, env, sshd_conf, _ = _box(tmp_path)
    before = sshd_conf.read_text()
    (fake / "sshd-refuses").write_text("")
    r = _run(script, env, "add", "anna", "--no-invite")
    assert r.returncode != 0 and "previous one is back" in r.stderr
    assert sshd_conf.read_text() == before
    assert not list(sshd_conf.parent.glob("40-delfin-hardening.conf.*"))


def test_the_invitation_leads_into_the_persons_own_account(tmp_path):
    fake, script, env, _, _ = _box(tmp_path)
    r = _run(script, env, "add", "anna", "--host", "box.example", "--hours", "6")
    assert r.returncode == 0, r.stderr + r.stdout
    admin_zip = fake / "home" / "pi" / "invites" / "DELFIN-anna.zip"
    assert admin_zip.is_file(), r.stdout
    assert oct(admin_zip.stat().st_mode & 0o777) == "0o600"
    with zipfile.ZipFile(admin_zip) as z:
        invite = json.loads(z.read("windows-launcher/invite.json"))
        assert "windows-launcher/invite_key" in z.namelist()
    assert invite["target"] == "anna@box.example"
    assert invite["host_keys"] and invite["host_keys"][0].startswith("ssh-ed25519 ")
    # The key's line is in ANNA's authorized_keys, forced to enroll, and the
    # forced command is the root-owned delfin-people, not a copy she owns.
    home = fake / "home" / "anna"
    lines = [line for line in (home / ".ssh" / "authorized_keys").read_text().splitlines()
             if line.endswith(" delfin-invite-anna")]
    assert len(lines) == 1
    assert lines[0].startswith(f'restrict,command="{script.parent / "delfin-people"} --enroll anna",expiry-time="')
    assert not (home / "invites" / "DELFIN-anna.zip").exists(), "no copy is left in the person's home"


def test_delfin_people_forwards_the_person_verbs_and_says_when_the_variant_is_missing(tmp_path):
    text = PEOPLE.read_text(encoding="utf-8")
    for verb in ("--add-person", "--invite-person", "--lock-person", "--unlock-person",
                 "--remove-person", "--persons"):
        assert verb in text
    if not pathlib.Path("/usr/local/sbin/delfin-accounts").exists():
        r = subprocess.run(["bash", str(PEOPLE), "--add-person", "anna"], capture_output=True,
                           text=True, timeout=30, env=dict(os.environ, HOME=str(tmp_path)))
        assert r.returncode == 1 and "DELFIN_PI_ACCOUNTS=1" in r.stderr


# --------------------------------------------------------------------------
# The Windows launcher, from a person's login shell.
def _shared_install(shared: pathlib.Path, tmp_path: pathlib.Path) -> pathlib.Path:
    """/opt/delfin as pi_install.sh leaves it: src/ (the package), .venv/ a
    real venv whose site-packages points at src (the compat editable)."""
    src = shared / "src"
    for package in ("delfin", "delfin/agent", "delfin/dashboard", "voila"):
        (src / package).mkdir(parents=True, exist_ok=True)
        (src / package / "__init__.py").write_text("")
    (src / "delfin/agent/where.py").write_text("")
    (src / "delfin/dashboard/session.py").write_text("")
    (src / "delfin/cli_voila.py").write_text("def main(argv):\n    print('--strict-port')\n")
    venv = shared / ".venv"
    subprocess.run([sys.executable, "-m", "venv", "--without-pip", str(venv)], check=True,
                   env=child_env(tmp_path))
    site = subprocess.run([str(venv / "bin" / "python"), "-c",
                           "import sysconfig; print(sysconfig.get_paths()['purelib'])"],
                          check=True, capture_output=True, text=True, env=child_env(tmp_path)).stdout.strip()
    pathlib.Path(site, "__editable__.delfin.pth").write_text(f"{src}\n")
    (shared / "agent.env").write_text("export DELFIN_AGENT_SANDBOX=auto\n")
    (shared / "settings-defaults.json").write_text('{"agent": {"mode": "solo"}}\n')
    return venv / "bin" / "python"


def test_the_launcher_finds_the_shared_install_from_a_persons_login_shell(tmp_path):
    shared = tmp_path / "opt" / "delfin"
    python = _shared_install(shared, tmp_path)
    home = tmp_path / "home" / "anna"
    home.mkdir(parents=True)
    # Debian's skeleton: .profile reads .bashrc, which stops early unless
    # the shell is interactive. The launcher's login shell is both.
    (home / ".profile").write_text('if [ -n "$BASH_VERSION" ] && [ -f "$HOME/.bashrc" ]; then . "$HOME/.bashrc"; fi\n')
    block = subprocess.run(["bash", str(ACCOUNTS), "shell-block", str(shared)],
                           check=True, capture_output=True, text=True).stdout
    (home / ".bashrc").write_text("case $- in *i*) ;; *) return;; esac\n" + block)
    payload = tmp_path / "payload.py"
    payload.write_text('import os, sys\nprint("CHOSEN=" + sys.executable)\n'
                       'print("DEFAULTS=" + os.environ.get("DELFIN_SETTINGS_DEFAULTS", ""))\n'
                       'print("SANDBOX=" + os.environ.get("DELFIN_AGENT_SANDBOX", ""))\n')
    env = {k: v for k, v in child_env(tmp_path).items()
           if k not in ("VIRTUAL_ENV", "CONDA_PREFIX", "DELFIN_REPO", "PYTHONPATH", "DELFIN_SETTINGS_DEFAULTS",
                        "DELFIN_AGENT_SANDBOX", "BASH_ENV", "ENV")}
    env.update(PATH="/usr/bin:/bin", HOME=str(home))
    # As the app runs it (DELFIN.ps1, Get-RemoteCommand): an interactive
    # login shell that sources the bootstrap with the arguments as data.
    command = f". {shlex.quote(str(LAUNCHER / 'remote_bootstrap.sh'))}"
    r = subprocess.run(["bash", "-ilc", command, "delfin", "0", "", "", str(payload), "dashboard"],
                       env=env, capture_output=True, text=True, timeout=60, cwd=home)
    assert r.returncode == 0, r.stderr + r.stdout
    assert f"CHOSEN={python}" in r.stdout
    assert f"DEFAULTS={shared}/settings-defaults.json" in r.stdout
    assert "SANDBOX=auto" in r.stdout


# --------------------------------------------------------------------------
# The host scripts.
def test_the_hardening_keeps_every_unlocked_person_in_the_allow_list(tmp_path):
    text = (PI / "pi_harden_host.sh").read_text(encoding="utf-8")
    assert "AllowUsers $LOGIN_USER $AGENT_USER$PEOPLE_ALLOW\n" in text
    start = text.index("people_allowed() {")
    func = text[start:text.index("\n}\n", start) + 3]
    fake = tmp_path / "fake"
    fake.mkdir()
    (fake / "group").write_text("delfin-people:x:900:anna,ben,cleo\n")
    (fake / "shadow").write_text("anna:!:1:0:99999:7:::\nben:!:1:0:99999:7::1:\ncleo:!:1:0:99999:7::99999:\n")
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "getent").write_text(STUBS["getent"])
    (bin_dir / "getent").chmod(0o755)
    r = subprocess.run(["bash", "-c", func + 'printf "[%s]" "$(people_allowed)"'],
                       env=dict(os.environ, FAKE=str(fake), PATH=f"{bin_dir}:/usr/bin:/bin"),
                       capture_output=True, text=True, timeout=30)
    assert r.stdout == "[ anna cleo]", r.stderr
    # Without the group: the list stays exactly "pi delfin".
    (fake / "group").write_text("")
    r = subprocess.run(["bash", "-c", func + 'printf "[%s]" "$(people_allowed)"'],
                       env=dict(os.environ, FAKE=str(fake), PATH=f"{bin_dir}:/usr/bin:/bin"),
                       capture_output=True, text=True, timeout=30)
    assert r.stdout == "[]"


def test_the_shared_install_is_opt_in_and_resolves_to_the_shared_path():
    boot = (PI / "pi_bootstrap.sh").read_text(encoding="utf-8")
    assert 'ACCOUNTS="${DELFIN_PI_ACCOUNTS:-0}"' in boot
    # Root-owned copies: never run as root from the agent user's tree.
    assert "install -D -m 755 -o root -g root \"$HERE/delfin_accounts.sh\" /usr/local/sbin/delfin-accounts" in boot
    assert 'UV_PYTHON_INSTALL_DIR="$SHARED_ROOT/python"' in boot
    # Root's, open to the readers' group, not to every account on the box.
    assert 'install -d -m 750 -o root -g delfin-shared "$SHARED_ROOT"' in boot
    install = (PI / "pi_install.sh").read_text(encoding="utf-8")
    assert 'ROOT="$(readlink -f "$HOME/delfin"' in install
    deploy = (PI / "deploy_pi.sh").read_text(encoding="utf-8")
    assert "delfin_accounts.sh delfin_people.sh delfin_deploy.sh" in deploy


def test_pi_install_takes_a_linked_root_as_shared_and_a_plain_one_as_before(tmp_path):
    install = (PI / "pi_install.sh").read_text(encoding="utf-8")
    start = install.index('MODE="${1:-}"')
    head = install[start:install.index('VENV_PY="$ROOT/.venv/bin/python"', start)]
    home = tmp_path / "home"
    home.mkdir()
    (tmp_path / "opt").mkdir()
    probe = head + 'printf "%s %s %s" "$ROOT" "$SHARED" "$(umask)"'
    (home / "delfin").mkdir()
    r = subprocess.run(["bash", "-c", "umask 077; " + probe], env=dict(os.environ, HOME=str(home)),
                       capture_output=True, text=True, timeout=30)
    assert r.stdout == f"{os.path.realpath(home)}/delfin 0 0077"
    (home / "delfin").rmdir()
    (home / "delfin").symlink_to(tmp_path / "opt")
    # The agent account's part of a shared deploy (pi_install.sh account).
    r = subprocess.run(["bash", "-c", "umask 077; " + probe, "_", "account"],
                       env=dict(os.environ, HOME=str(home)), capture_output=True, text=True, timeout=30)
    assert r.stdout == f"{os.path.realpath(tmp_path / 'opt')} 1 0022"
