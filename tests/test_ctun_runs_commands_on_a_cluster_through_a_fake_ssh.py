"""The parts of ctun's Python API that DELFIN's remote layer is built on.

Input: ``pyproject.toml``; ``cluster_tunnel.api`` with a ctun config in a
temporary directory and a fake ``ssh`` first on PATH. The fake answers
``ssh -O check`` from a marker file that stands for a logged-in master,
records every argv it is given, and runs the remote command line with a local
``sh`` in place of the cluster's login shell. Output: assertions.

cluster-tunnel (ctun) is a core dependency: it keeps the logged-in ssh
connection to a cluster that DELFIN sends calculations through. These tests
pin what the remote layer relies on, so that a ctun release that changes it
fails here and not on a cluster: no master means LoginRequired (ctun's exit
code 10), a command travels as an argv list that the remote shell cannot
re-read, results come back captured with their exit code, and a cluster whose
host key is pinned runs every ssh against that key alone. ctun's state
(config, sockets, sessions, command log) lives in the temporary directory,
never in the user's own ~/.config or ~/.cache.
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import pytest
import tomllib

from cluster_tunnel import api

REPO = Path(__file__).resolve().parents[1]

_FAKE_SSH = """\
#!{python}
import json, os, subprocess, sys

state = os.environ["FAKE_SSH_STATE"]
args = sys.argv[1:]
with open(os.path.join(state, "argv.jsonl"), "a", encoding="utf-8") as log:
    log.write(json.dumps(args) + "\\n")
if "-O" in args:
    if args[args.index("-O") + 1] == "check" and os.path.exists(os.path.join(state, "master")):
        print("Master running (pid=4242)", file=sys.stderr)
        sys.exit(0)
    print("Control socket connect(sock): No such file or directory", file=sys.stderr)
    sys.exit(255)
# ssh hands the words after the destination to the login shell on the
# cluster as one command line; a local sh stands in for that shell.
sys.exit(subprocess.run(["sh", "-c", args[-1]], cwd=os.path.join(state, "home")).returncode)
"""


@pytest.fixture
def cluster(tmp_path, monkeypatch):
    """A ctun config with one cluster, ``fake``, reached through the fake ssh."""
    state = tmp_path / "ssh"
    (state / "home").mkdir(parents=True)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    fake = bin_dir / "ssh"
    fake.write_text(_FAKE_SSH.format(python=sys.executable), encoding="utf-8")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:/usr/bin:/bin")
    monkeypatch.setenv("FAKE_SSH_STATE", str(state))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    config = tmp_path / "ctun.yaml"
    config.write_text("clusters:\n  fake:\n    host: fake.example\n    user: me\n",
                      encoding="utf-8")

    class Cluster:
        path = config
        home = state / "home"

        @staticmethod
        def log_in():
            (state / "master").touch()

        @staticmethod
        def argvs() -> list[list[str]]:
            log = state / "argv.jsonl"
            if not log.exists():
                return []
            return [json.loads(line) for line in log.read_text(encoding="utf-8").splitlines()]

    return Cluster


def test_cluster_tunnel_is_a_core_dependency_from_the_release_with_the_api():
    """0.10.0: ``cluster_tunnel.api`` with host-key pinning (0.7.0), a
    login with given credentials (the window's masked field, 0.8.0), and
    the presets and add/remove calls of the window's add form; the ceiling
    is the next minor, since a 0.x minor release may change the API."""
    project = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))["project"]
    ctun = [r for r in project["dependencies"]
            if re.split(r"[<>=!~\[; ]", r, maxsplit=1)[0].lower() == "cluster-tunnel"]
    assert ctun == ["cluster-tunnel>=0.10,<0.11"]


def test_without_a_master_a_command_needs_a_login_first(cluster):
    assert api.is_live("fake", config_path=cluster.path) is False
    with pytest.raises(api.LoginRequired) as caught:
        api.run("fake", ["true"], capture=True, config_path=cluster.path)
    assert caught.value.exit_code == 10
    # Refused before anything but the liveness check reached ssh.
    assert cluster.argvs() and all("-O" in argv for argv in cluster.argvs())


def test_a_command_comes_back_with_its_output_and_exit_code(cluster):
    cluster.log_in()
    assert api.is_live("fake", config_path=cluster.path) is True
    result = api.run("fake", ["sh", "-c", "echo out; echo err >&2; exit 3"],
                     capture=True, config_path=cluster.path)
    assert (result.exit_code, result.stdout, result.stderr) == (3, "out\n", "err\n")
    assert not result.ok


def test_the_remote_shell_sees_each_word_as_one_argument(cluster):
    """The remote layer will pass file and job names through; none of them
    may be read as shell syntax on the login node."""
    cluster.log_in()
    words = ["a b", "$HOME", "x; touch pwned", "`touch pwned2`", "'q'"]
    result = api.run("fake", ["printf", "%s\\n", *words], capture=True,
                     config_path=cluster.path)
    assert result.ok, result.stderr
    assert result.stdout.splitlines() == words
    assert not list(cluster.home.glob("pwned*"))


def test_an_empty_command_is_refused_instead_of_opening_a_login_shell(cluster):
    cluster.log_in()
    with pytest.raises(ValueError):
        api.run("fake", [], capture=True, config_path=cluster.path)


def test_a_pinned_cluster_runs_every_ssh_against_its_own_key_file(cluster, tmp_path):
    """DELFIN will use a cluster only once its host key is pinned; this is
    what the pin has to mean for every ssh that ctun starts."""
    known_hosts = tmp_path / "known_hosts"
    cluster.path.write_text(
        "clusters:\n  fake:\n    host: fake.example\n    user: me\n"
        "    strict_host_key_checking: 'yes'\n"
        f"    known_hosts_file: {known_hosts}\n",
        encoding="utf-8")
    cluster.log_in()
    assert api.run("fake", ["true"], capture=True, config_path=cluster.path).ok
    argvs = cluster.argvs()
    assert len(argvs) >= 2, argvs  # the liveness check and the command
    for argv in argvs:
        assert "StrictHostKeyChecking=yes" in argv, argv
        assert f'UserKnownHostsFile="{known_hosts}"' in argv, argv
        assert "GlobalKnownHostsFile=/dev/null" in argv, argv


def test_credentials_for_a_login_never_print():
    creds = api.Credentials("s3cret-pw", "424242")
    assert "s3cret-pw" not in repr(creds) and "424242" not in repr(creds)


def test_the_errors_the_remote_layer_maps_keep_their_exit_codes():
    assert (api.LoginRequired.exit_code, api.BudgetExhausted.exit_code,
            api.BudgetGuardError.exit_code) == (10, 11, 12)
    for error in (api.LoginRequired, api.BudgetExhausted, api.BudgetGuardError,
                  api.ConfigError, api.UnknownClusterError, api.HostKeyError):
        assert issubclass(error, api.ClusterTunnelError), error


def test_ctun_keeps_its_state_in_the_temporary_directories(cluster, tmp_path):
    cluster.log_in()
    api.run("fake", ["true"], capture=True, config_path=cluster.path)
    assert any((tmp_path / "cache").rglob("*")), "ctun wrote nothing under XDG_CACHE_HOME"
