"""The service lists DELFIN's clusters and hands a typed login to ctun, and only to ctun.

Input: a private agent service (tests/tui_service.py) whose HOME holds
ctun's config (one pinned cluster, ``kcist``, and an unpinned one) and
DELFIN's ``remote`` settings; a fake ``ssh`` first on its PATH that acts as
the login node: asked to open a master it prompts for the password on the
pty ctun gives it, and answers ``ssh -O check`` from a marker the right
password creates. Output: ``clusters.list`` and ``cluster.login`` results,
every event the window received, and every file under the test's
directories.

Semantics:

* each cluster of remote.clusters is listed with its state: refused (with
  its problems), logged_out, logged_in; and with what a login asks for;
* a login with the right password opens the tunnel through ctun's pty login,
  and the row then says logged_in; a wrong one is answered with
  login_failed and a message that does not hold it;
* the typed password is in no event any window received, and in no file the
  service, its workers or ctun wrote -- logs, event logs, sessions included;
* a connection from an agent's own process may list the clusters but not log
  in (or out), and its attempt starts no ssh.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent.service import protocol
from delfin.agent.service.client import ServiceError
from tui_service import REPO, ServiceProcess

PASSWORD = "pw-3f9c-not-anywhere-else"
WRONG = "wrong-7d21-not-anywhere-else"

_FAKE_SSH = """\
#!{python}
# The login node, as ctun's ssh sees it. -M: prompt for the password on the
# pty and open the "master" (a marker) only for the right one. -O check:
# the marker. -O exit: remove it. Anything else: run it with sh in $HOME.
import os, subprocess, sys
state = os.environ["FAKE_STATE"]
master = os.path.join(state, "master")
args = sys.argv[1:]
with open(os.path.join(state, "calls"), "a", encoding="utf-8") as log:
    log.write(("-M " if "-M" in args else "") + ("-O " + args[args.index("-O") + 1]
              if "-O" in args else "") + "\\n")
if "-O" in args:
    op = args[args.index("-O") + 1]
    if op == "exit":
        if os.path.exists(master):
            os.remove(master)
        sys.exit(0)
    sys.exit(0 if os.path.exists(master) else 255)
if "-M" in args:
    sys.stdout.write("me@login.example's password: ")
    sys.stdout.flush()
    got = sys.stdin.readline().rstrip("\\r\\n")
    if got == os.environ["FAKE_PASSWORD"]:
        open(master, "w").close()
        sys.exit(0)
    sys.stdout.write("Permission denied, please try again.\\n")
    sys.exit(255)
WITH_VALUE = set("bcDEeFIiJLlmOoPpQRSWw")
i = 0
while i < len(args) and args[i].startswith("-"):
    i += 2 if len(args[i]) == 2 and args[i][1] in WITH_VALUE else 1
sys.exit(subprocess.run(["sh", "-c", " ".join(args[i + 1:])], cwd=os.environ["HOME"]).returncode)
"""


@pytest.fixture
def service(tmp_path):
    svc = ServiceProcess(tmp_path, start=False)
    home = svc.home
    ctun_dir = home / ".config" / "cluster-tunnel"
    ctun_dir.mkdir(parents=True)
    (ctun_dir / "config.yaml").write_text(
        "clusters:\n"
        "  kcist:\n    host: login.example\n    user: me\n"
        "    strict_host_key_checking: 'yes'\n"
        f"    known_hosts_file: {ctun_dir / 'known_hosts'}\n"
        "  open:\n    host: open.example\n    user: me\n", encoding="utf-8")
    (home / ".delfin_settings.json").write_text(json.dumps(
        {"remote": {"clusters": {"kcist": {"default_partition": "small"}, "open": {}}}}),
        encoding="utf-8")
    state = tmp_path / "ssh-state"
    state.mkdir()
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "ssh").write_text(_FAKE_SSH.format(python=sys.executable), encoding="utf-8")
    (bin_dir / "ssh").chmod(0o755)
    svc.env.update(PATH=f"{bin_dir}:{svc.env.get('PATH', '/usr/bin:/bin')}",
                   XDG_CONFIG_HOME=str(home / ".config"), XDG_CACHE_HOME=str(home / ".cache"),
                   FAKE_STATE=str(state), FAKE_PASSWORD=PASSWORD)
    svc.env.pop("CTUN_CONFIG", None)
    svc.state = state
    svc.start()
    try:
        yield svc
    finally:
        svc.stop()


def _rows(c) -> dict:
    return {r["name"]: r for r in c.request("clusters.list", timeout=60)["clusters"]}


def _files_holding(root: Path, secret: str) -> list[str]:
    needle = secret.encode()
    hits = []
    for path in root.rglob("*"):
        try:
            if path.is_file() and not path.is_symlink() and needle in path.read_bytes():
                hits.append(str(path))
        except OSError:
            continue
    return hits


def test_each_cluster_is_listed_with_its_state(service):
    c = service.client()
    rows = _rows(c)
    assert rows["kcist"]["state"] == "logged_out"
    assert (rows["kcist"]["requiresPassword"], rows["kcist"]["requiresOtp"]) == (True, False)
    assert rows["kcist"]["defaultPartition"] == "small"
    assert rows["kcist"]["target"] == "me@login.example"
    assert rows["open"]["state"] == "refused"
    assert any("hostkey --pin" in p for p in rows["open"]["problems"])


def test_a_typed_password_opens_the_tunnel_and_goes_nowhere_else(service):
    c = service.client()
    sid = service.session(c)          # a running session: its worker must see nothing
    wrong = c.request("cluster.login", {"name": "kcist", "password": WRONG}, timeout=120)
    assert (wrong["ok"], wrong["code"]) == (False, "login_failed")
    assert WRONG not in json.dumps(wrong)
    assert _rows(c)["kcist"]["state"] == "logged_out"

    done = c.request("cluster.login", {"name": "kcist", "password": PASSWORD}, timeout=120)
    assert done == {"ok": True, "opened": True, "code": None, "error": None}
    assert _rows(c)["kcist"]["state"] == "logged_in"
    again = c.request("cluster.login", {"name": "kcist", "password": PASSWORD}, timeout=120)
    assert (again["ok"], again["opened"]) == (True, False)

    c.request("session.status", {"sessionId": sid}, timeout=60)
    seen = json.dumps(c.events)
    assert PASSWORD not in seen and WRONG not in seen
    for secret in (PASSWORD, WRONG):
        assert _files_holding(service.tmp, secret) == [], secret
        assert _files_holding(Path(service.runtime), secret) == [], secret

    out = c.request("cluster.logout", {"name": "kcist"}, timeout=60)
    assert (out["ok"], out["closed"]) == (True, True)
    assert _rows(c)["kcist"]["state"] == "logged_out"


def test_an_unusable_cluster_is_not_logged_in(service):
    c = service.client()
    got = c.request("cluster.login", {"name": "open", "password": PASSWORD}, timeout=60)
    assert (got["ok"], got["code"]) == (False, "cluster_not_usable")
    assert not (service.state / "calls").exists() or "-M" not in (service.state / "calls").read_text()


def test_an_agents_own_process_may_look_but_not_log_in_or_read_details(service):
    probe = r"""
import json, sys
sys.path.insert(0, sys.argv[2])
from delfin.agent.service.client import ServiceClient, ServiceError
c = ServiceClient(sys.argv[1])
out = {"listed": len(c.request("clusters.list", timeout=60)["clusters"])}
for method, params in (("cluster.login", {"name": "kcist", "password": sys.argv[3]}),
                       ("cluster.logout", {"name": "kcist"}),
                       ("cluster.detail", {"name": "kcist"}),
                       ("job.detail", {"source": "cluster", "cluster": "kcist", "id": "1"})):
    try:
        c.request(method, params, timeout=60)
        out[method] = "done"
    except ServiceError as exc:
        out[method] = exc.code
print(json.dumps(out))
"""
    from conftest import child_env
    env = child_env(service.tmp)
    env.update(DELFIN_RUNTIME_DIR=service.env["DELFIN_RUNTIME_DIR"],
               PYTHONPATH=service.env["PYTHONPATH"], DELFIN_AGENT_SESSION="12345")
    run = subprocess.run([sys.executable, "-c", probe, str(service.socket), str(REPO), PASSWORD],
                         env=env, capture_output=True, text=True, timeout=120)
    got = json.loads(run.stdout.strip().splitlines()[-1])
    assert got == {"listed": 2, "cluster.login": protocol.NOT_ALLOWED,
                   "cluster.logout": protocol.NOT_ALLOWED, "cluster.detail": protocol.NOT_ALLOWED,
                   "job.detail": protocol.NOT_ALLOWED}
    calls = (service.state / "calls").read_text() if (service.state / "calls").exists() else ""
    assert "-M" not in calls, "the refused login started ssh"


def test_the_login_request_is_the_services_own_and_never_a_workers():
    from delfin.agent.service import supervisor
    assert "cluster.login" in supervisor.ACTING and "cluster.logout" in supervisor.ACTING
    assert "cluster.login" not in supervisor.FORWARDED
    assert "cluster.login" in protocol.SERVICE_REQUESTS
    assert supervisor.SERVICE_VERSION >= 4
