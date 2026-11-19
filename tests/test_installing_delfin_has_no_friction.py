"""A first install on a fresh account: what is missing is said first, a failed
step says how to go on, and the end says how to start DELFIN.

The installer is run here as a new user runs it -- an empty HOME, a PATH
holding only the host commands each test names -- with a stand-in for the
Python it builds the venv from: ``-m venv`` and ``-m pip`` are answered
without the network, everything else (reading the tool list, writing the
settings, the summary) runs the real interpreter on this checkout. curl is a
stand-in too where a test needs one, so nothing here reaches the network.
"""

import json
import os
import pathlib
import shutil
import subprocess
import sys

import pytest

REPO = pathlib.Path(__file__).resolve().parents[1]
INSTALL_SH = REPO / "install.sh"
INSTALLER = REPO / "delfin" / "installers" / "install_delfin.sh"
ENTRY_SCRIPTS = (
    INSTALL_SH,
    INSTALLER,
    REPO / "delfin" / "installers" / "install_delfin_bwu.sh",
    REPO / "scripts" / "install_delfin_bwu.sh",
    REPO / "scripts" / "verify_delfin_bwu.sh",
)

# What a minimal Linux account has, without curl, wget, a compiler, make,
# an MPI or a Python of its own.
HOST_COMMANDS = (
    "bash", "sh", "env", "uname", "dirname", "basename", "awk", "sed", "grep",
    "cut", "sort", "ls", "mkdir", "date", "readlink", "nproc", "tr", "head",
    "tail", "df", "mv", "rm", "cp", "cat", "tar", "chmod", "tee", "id",
    "hostname", "mktemp", "touch", "wc", "sleep", "find", "xargs", "stat",
)


def _host_bin(tmp_path, *, without=(), stubs=None):
    """A directory with the host commands (links) and the stubs given."""
    bindir = tmp_path / "hostbin"
    bindir.mkdir(exist_ok=True)
    for name in HOST_COMMANDS:
        if name in without:
            continue
        found = shutil.which(name)
        if found and not (bindir / name).exists():
            (bindir / name).symlink_to(found)
    for name, text in (stubs or {}).items():
        path = bindir / name
        if path.exists() or path.is_symlink():
            path.unlink()
        path.write_text(text)
        path.chmod(0o755)
    return bindir


def _stand_in_python(tmp_path):
    """A Python 3.11 that makes venvs and pip-installs without the network."""
    bash = shutil.which("bash")
    path = tmp_path / "standin" / "python3.11"
    path.parent.mkdir(exist_ok=True)
    path.write_text(f"""#!{bash}
real={sys.executable!s}
if [ "${{1:-}}" = "-c" ]; then
  case "${{2:-}}" in *version_info*exit*|*exit*version_info*) exit "${{STANDIN_TOO_NEW:-0}}" ;; esac
fi
case "${{1:-}} ${{2:-}}" in
  "-m venv") mkdir -p "$3/bin" && cp "$0" "$3/bin/python" && exit 0 ;;
  "-m pip") printf 'pip %s\\n' "$*" >> "$STANDIN_LOG"; exit "${{STANDIN_PIP_EXIT:-0}}" ;;
esac
exec "$real" "$@"
""")
    path.chmod(0o755)
    if not (path.parent / "python3").is_symlink():
        (path.parent / "python3").symlink_to(path)
    return path


def _recording(name):
    return f'#!/bin/sh\nprintf "{name} %s\\n" "$*" >> "$STANDIN_LOG"\nexit 0\n'


def _checkout(tmp_path):
    """A directory the installer takes for a DELFIN checkout."""
    root = tmp_path / "checkout"
    root.mkdir(exist_ok=True)
    (root / "pyproject.toml").write_text('[project]\nname = "delfin-complat"\n')
    return root


def _run(tmp_path, *args, script=INSTALLER, bindir=None, env=None, home=None):
    home = home or tmp_path / "home"
    home.mkdir(exist_ok=True)
    bindir = bindir or _host_bin(tmp_path, stubs={"curl": _recording("curl")})
    standin = _stand_in_python(tmp_path)
    full_env = {
        "HOME": str(home),
        "PATH": os.pathsep.join([str(standin.parent), str(bindir)]),
        "PYTHONPATH": str(REPO),
        "STANDIN_LOG": str(tmp_path / "standin.log"),
        "LANG": "C.UTF-8",
    }
    full_env.update(env or {})
    return subprocess.run(["bash", str(script), *args], capture_output=True, text=True,
                          timeout=300, env=full_env, stdin=subprocess.DEVNULL), home


@pytest.mark.parametrize("script", ENTRY_SCRIPTS, ids=lambda p: p.name)
def test_every_entry_script_runs_when_started_with_sh(script):
    """`sh install.sh` is dash on Debian and Ubuntu; set -o pipefail stopped it
    at its first line with "Illegal option -o pipefail"."""
    sh = shutil.which("sh")
    done = subprocess.run([sh, str(script), "--help"], capture_output=True, text=True,
                          timeout=60, stdin=subprocess.DEVNULL)
    assert done.returncode == 0, done.stderr
    assert "Illegal option" not in done.stderr
    assert "Usage" in done.stdout


def test_no_curl_and_no_wget_is_said_before_anything_changes(tmp_path):
    bindir = _host_bin(tmp_path, without=("curl", "wget"))
    prefix = tmp_path / "software"
    done, home = _run(tmp_path, "--profile", "core", "--no-openmpi", "--no-orca",
                      "--prefix", str(prefix), "--repo", str(_checkout(tmp_path)),
                      bindir=bindir)

    assert done.returncode != 0
    assert "curl" in done.stderr and "wget" in done.stderr
    assert not prefix.exists(), "the prefix was created before the check refused"
    assert list(home.iterdir()) == [], "HOME was written before the check refused"
    assert not (tmp_path / "checkout" / ".venv").exists()


def test_a_prefix_that_cannot_be_created_is_refused_up_front(tmp_path):
    blocker = tmp_path / "a-file"
    blocker.write_text("")
    done, home = _run(tmp_path, "--profile", "core", "--no-openmpi", "--no-orca",
                      "--prefix", str(blocker / "software"),
                      "--repo", str(_checkout(tmp_path)))

    assert done.returncode != 0
    assert "--prefix" in done.stderr
    assert str(blocker) in done.stderr
    assert list(home.iterdir()) == []


def test_a_missing_compiler_is_named_before_the_long_steps(tmp_path):
    """Building OpenMPI needs a C compiler and make; finding that out after the
    Python install is twenty minutes too late."""
    bindir = _host_bin(tmp_path, without=("gcc", "cc", "make"),
                       stubs={"curl": _recording("curl")})
    done, _ = _run(tmp_path, "--dry-run", "--profile", "core", "--no-orca",
                   "--repo", str(_checkout(tmp_path)), bindir=bindir)

    assert done.returncode == 0, done.stderr
    assert "C compiler" in done.stderr
    assert "--no-openmpi" in done.stderr


def test_little_free_disk_space_is_said_up_front(tmp_path):
    tiny_df = ('#!/bin/sh\n'
               'echo "Filesystem 1024-blocks Used Available Capacity Mounted on"\n'
               'echo "/dev/x 1000000 900000 100000 90% /"\n')
    bindir = _host_bin(tmp_path, stubs={"curl": _recording("curl"), "df": tiny_df})
    done, _ = _run(tmp_path, "--dry-run", "--profile", "core", "--no-openmpi", "--no-orca",
                   "--repo", str(_checkout(tmp_path)), bindir=bindir)

    assert done.returncode == 0, done.stderr
    assert "free" in done.stderr and "GB" in done.stderr


def _settings(home):
    return json.loads((home / ".delfin_settings.json").read_text())


def test_a_fresh_install_ends_with_how_to_start_and_runs_again_cleanly(tmp_path):
    from delfin.user_settings import DEFAULT_SETTINGS

    checkout = _checkout(tmp_path)
    prefix = tmp_path / "software"
    home = tmp_path / "home"
    home.mkdir()
    old_activation = "source /elsewhere/DELFIN/.venv/bin/activate"
    (home / ".bashrc").write_text(f"{old_activation}\n")
    args = ("--profile", "core", "--no-openmpi", "--no-orca",
            "--prefix", str(prefix), "--repo", str(checkout))

    first, _ = _run(tmp_path, *args, home=home)

    assert first.returncode == 0, first.stdout + first.stderr
    out = first.stdout
    assert "delfin-voila" in out
    assert "delfin-agent" in out
    assert f"bash {checkout}/install.sh --update" in out
    assert "ssh -L" in out, "the summary does not say how to reach a dashboard on a login node"
    assert old_activation in out + first.stderr, \
        "an older DELFIN activated from ~/.bashrc is not named"
    logs = sorted((prefix / "logs").glob("delfin-install-*.log"))
    assert len(logs) == 1, "no log of the install was kept"
    assert "delfin-voila" in logs[0].read_text()
    assert str(logs[0]) in out
    # The dashboard's own tab defaults stay in charge: the installer's settings
    # write keeps ui.tabs as shipped.
    saved_tabs = (_settings(home).get("ui") or {}).get("tabs", DEFAULT_SETTINGS["ui"]["tabs"])
    assert saved_tabs == DEFAULT_SETTINGS["ui"]["tabs"]

    second, _ = _run(tmp_path, *args, home=home)

    assert second.returncode == 0, second.stdout + second.stderr
    rc = (home / ".bashrc").read_text()
    assert rc.count(".delfin_env.sh") == 1, rc
    assert not list(checkout.glob(".venv.before-*")), "the second run moved the venv aside"
    assert len(list((prefix / "logs").glob("delfin-install-*.log"))) >= 1


def test_a_failed_step_names_itself_and_how_to_go_on(tmp_path):
    checkout = _checkout(tmp_path)
    prefix = tmp_path / "software"
    done, _ = _run(tmp_path, "--profile", "core", "--no-openmpi", "--no-orca",
                   "--prefix", str(prefix), "--repo", str(checkout),
                   env={"STANDIN_PIP_EXIT": "1"})

    assert done.returncode != 0
    err = done.stderr
    assert "installing DELFIN" in err, err
    assert "run the same command again" in err, err
    logs = sorted((prefix / "logs").glob("delfin-install-*.log"))
    assert logs and str(logs[0]) in err
    assert "delfin-voila" not in done.stdout, "a failed install printed the start instructions"


def test_a_failed_micromamba_download_says_what_to_do(tmp_path):
    """No Python 3.10/3.11 here and the download of micromamba fails (no
    connection, a proxy): the message used to be only "pass --python PATH"."""
    failing_curl = ('#!/bin/sh\necho "curl: (7) Failed to connect to github.com port 443" >&2\n'
                    'exit 7\n')
    bindir = _host_bin(tmp_path, stubs={"curl": failing_curl})
    prefix = tmp_path / "software"
    done, _ = _run(tmp_path, "--profile", "core", "--no-openmpi", "--no-orca",
                   "--prefix", str(prefix), "--repo", str(_checkout(tmp_path)),
                   bindir=bindir, env={"STANDIN_TOO_NEW": "1"})

    assert done.returncode != 0
    err = done.stderr
    assert "micromamba could not be downloaded" in err, err
    assert "https_proxy" in err
    assert "--python PATH" in err
    assert "run the same command again" in err
    bindir_mm = prefix / "micromamba" / "bin"
    assert not (bindir_mm / "micromamba").exists()
    assert not (bindir_mm.exists() and list(bindir_mm.glob(".micromamba.partial.*")))
    assert not (tmp_path / "checkout" / ".venv").exists()
