"""No agent's cage reaches a pairing request, a door or activation records of anyone.

Review finding: a caged agent could take another person's invitation link.
(a) The cages masked ``.delfin/messenger_pairing`` only under the agent's
own HOME, not under ``~/people/*`` nor the account's home inside a
person's cage, and the home secret masks came BEFORE the workspace bind,
so a workspace at $HOME re-exposed every one of them. (b) The bridge
handed the link to any door socket of this uid, also one an agent process
listened on.

Output, each with its control:

* both bubblewrap builders mask, after every bind, the home secrets and,
  in people mode, every ``~/people/*/.delfin/messenger_pairing``,
  ``.../messenger_sessions`` and ``~/people/*/.run``, plus the account's
  own pairing and activation directories inside a person's cage; Landlock
  hides the same paths;
* a request whose door is held by an agent process -- or by any process
  other than the one the request names (pid and start time) -- gets
  nothing delivered; the invitation is withdrawn and audited. The asking
  process's own door gets the link (control).
"""
from __future__ import annotations

import os
import subprocess
import sys
import time
from pathlib import Path

import pytest

from _messenger_people_rig import Account
from conftest import requires_a_working_sandbox
from delfin.agent import api_client as A
from delfin.agent import sandbox, security_events
from delfin.messenger import pairing, store


@pytest.fixture
def account_home(tmp_path, monkeypatch):
    home = tmp_path / "acct"
    (home / ".config" / "delfin").mkdir(parents=True)
    (home / ".config" / "delfin" / "launcher-people.conf").write_text("on\n")
    (home / ".ssh").mkdir()
    for d in (home, home / "people" / "anna", home / "people" / "ben"):
        (d / ".delfin" / "messenger_pairing").mkdir(parents=True)
        (d / ".delfin" / "messenger_sessions").mkdir(parents=True)
        (d / ".run").mkdir(parents=True, exist_ok=True)
    return home


def _secrets(home, own=None):
    """The paths every cage must mask: the cage owner's own home secrets
    (``.ssh`` as the sample) and everyone's pairing and activation
    directories and runtime directories."""
    own = own or home
    out = [own / ".ssh", home / ".delfin" / "messenger_pairing", home / ".delfin" / "messenger_sessions"]
    for name in ("anna", "ben"):
        p = home / "people" / name
        out += [p / ".delfin" / "messenger_pairing", p / ".delfin" / "messenger_sessions", p / ".run"]
    return [str(x.resolve()) for x in out]


def _masked_after_binds(argv, path, bound):
    pairs = list(zip(argv, argv[1:]))
    masks = [i for i, (a, b) in enumerate(pairs) if a == "--tmpfs" and b == path]
    binds = [i for i, (a, b) in enumerate(pairs) if a in ("--bind", "--bind-try") and b == bound]
    return bool(masks) and (not binds or max(masks) > max(binds))


@pytest.mark.parametrize("who", ["owner", "anna"])
def test_both_bwrap_builders_mask_every_persons_doors_after_the_workspace_bind(
        account_home, monkeypatch, who):
    home = account_home if who == "owner" else account_home / "people" / "anna"
    (home / ".ssh").mkdir(exist_ok=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    monkeypatch.setattr(A.shutil, "which", lambda name, *a, **k: f"/usr/bin/{name}")
    # run_agent_command's cage with the repository at $HOME
    bound = str(home.resolve())
    argv = sandbox._bwrap_argv("true", home, allow_network=False)
    for path in _secrets(account_home, home):
        assert _masked_after_binds(argv, path, bound), (path, argv)
    # the shell tool's cage with the account's home granted as a write root
    ws = home / "ws"
    ws.mkdir()
    perms = A.KitToolPermissions(workspace=ws.resolve(), mode="bypassPermissions",
                                 extra_workspace_dirs=(account_home / "people", home / ".delfin"))
    argv2 = A._bash_isolation_argv("true", ws, perms, mode="bwrap")
    for path in _secrets(account_home, home):
        for root in (account_home / "people", home / ".delfin"):
            assert _masked_after_binds(argv2, path, str(root.resolve())), (path, argv2)


def test_landlock_hides_every_persons_doors(account_home, monkeypatch):
    anna = account_home / "people" / "anna"
    (anna / ".ssh").mkdir(exist_ok=True)
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: anna))
    hides = A._home_secret_paths()
    for path in _secrets(account_home, anna):
        assert path in hides, path


@requires_a_working_sandbox
def test_a_caged_command_sees_no_pairing_request_of_another_person(account_home, monkeypatch, tmp_path):
    (account_home / "people" / "ben" / ".delfin" / "messenger_pairing" / "req-x.json").write_text("{}")
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: account_home / "people" / "anna"))
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setenv("DELFIN_AGENT_SANDBOX", "bwrap")
    res = sandbox.run_agent_command(
        f"ls -A {account_home}/people/ben/.delfin/messenger_pairing | wc -l", account_home)
    assert res.stdout.strip().splitlines()[-1] == "0"


_LISTENER = r"""
import json, os, socket, sys
path = sys.argv[1]
s = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
s.bind(path); os.chmod(path, 0o600); s.listen(1); s.settimeout(8)
print("ready", flush=True)
try:
    c, _ = s.accept()
    c.settimeout(2)
    data = c.recv(65536)
    c.sendall(b'{"ok": true}\n')
    print("GOT", data.decode(errors="replace"), flush=True)
except OSError:
    print("NOTHING", flush=True)
"""


@pytest.fixture
def acct(tmp_path, monkeypatch):
    security_events.clear()
    yield Account(tmp_path, monkeypatch)
    security_events.clear()


def _forged_door(acct, person="anna", marker=True):
    from delfin.agent.approval_guard import ENV_MARKER
    path = os.path.join(str(acct.runtime[person]), "pair-evil.sock")
    env = dict(os.environ)
    if marker:
        env[ENV_MARKER] = "1"
    proc = subprocess.Popen([sys.executable, "-c", _LISTENER, path], env=env,
                            stdout=subprocess.PIPE, text=True)
    assert proc.stdout.readline().strip() == "ready"
    return proc, path


def test_an_agent_held_door_gets_nothing_and_the_invitation_is_withdrawn(acct):
    proc, path = _forged_door(acct)
    rid = pairing.submit("pair", {"requester": "x", "nonce_sha": "a" * 64, "door": path,
                                  "pid": proc.pid, "start": pairing._start_time(proc.pid)},
                         directory=acct.pairing_dir("anna"))
    acct.round()
    out, _ = proc.communicate(timeout=15)
    assert "GOT" not in out
    res = pairing.result(rid, directory=acct.pairing_dir("anna"))
    assert res["state"] == "refused" and acct.daemon.dropped and not store.invitations()
    assert "messenger_invitation_refused" in [e.kind for e in security_events.recent()]


def test_a_door_of_another_process_than_the_one_named_gets_nothing(acct):
    proc, path = _forged_door(acct, marker=False)
    rid = pairing.submit("pair", {"requester": "x", "nonce_sha": "a" * 64, "door": path,
                                  "pid": os.getpid(), "start": pairing._start_time(os.getpid())},
                         directory=acct.pairing_dir("anna"))
    acct.round()
    out, _ = proc.communicate(timeout=15)
    assert "GOT" not in out
    assert pairing.result(rid, directory=acct.pairing_dir("anna"))["state"] == "refused"


def test_the_asking_process_gets_its_link_control(acct):
    rid = pairing.begin("session:p1", directory=acct.pairing_dir("anna"))
    acct.round()
    assert pairing.read(rid, "session:p1", directory=acct.pairing_dir("anna"))["link"] == acct.daemon.links[0]
    time.sleep(0)
