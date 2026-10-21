"""A refused dashboard login is recorded as a labelled security event.

Input: a wrong token -- on the login form, in the URL, in the
``Authorization`` header -- sent to a real ``delfin.cli_voila`` server
(``--no-browser``, HOME in tmp_path, no desktop session so nothing is
shown on a screen). Output: a WARNING line in the server log and a durable
``security_alert`` in the attention inbox under the server's HOME, at most
one record per minute however many refusals arrive, naming the client and
the path and never the token that was presented.

Unit tests drive the rate limit with a fake clock.
"""

from __future__ import annotations

import http.client
import json
import socket
import subprocess
import sys
import time
import urllib.parse
from pathlib import Path

import pytest
from conftest import child_env

from delfin.agent import security_events
from delfin.dashboard import login_session

TOKEN = "refused-login-test-token-0123456789abcdef0123"
WRONG = "not-the-token-but-also-long-enough-0123456789"
REPO = Path(__file__).resolve().parents[1]


class Clock:
    def __init__(self):
        self.now = 5_000.0

    def __call__(self):
        return self.now


def _limiter(clock):
    seen = []
    limiter = login_session.RefusedLogins(
        window=60.0, clock=clock,
        record=lambda kind, tool, detail, blocked: seen.append(
            (kind, tool, detail, blocked)))
    return limiter, seen


def test_the_first_refusal_is_recorded_as_a_blocked_security_event():
    limiter, seen = _limiter(Clock())
    assert limiter.refused("wrong token in the request", "127.0.0.1",
                           "/api/status?token=" + WRONG)
    [(kind, tool, detail, blocked)] = seen
    assert kind == login_session.REFUSED_KIND and tool == "dashboard"
    assert blocked is True
    assert "127.0.0.1" in detail and "/api/status" in detail
    assert WRONG not in detail and "token=" not in detail


def test_a_flood_of_refusals_is_one_record_per_window_with_the_count():
    clock = Clock()
    limiter, seen = _limiter(clock)
    for _ in range(500):
        limiter.refused("wrong token on the login form", "10.0.0.7", "/login")
        clock.now += 0.1                     # 500 in 50 s
    assert len(seen) == 1
    clock.now += 60
    limiter.refused("wrong token on the login form", "10.0.0.7", "/login")
    assert len(seen) == 2
    assert "499 more" in seen[1][2]


def test_the_kind_has_a_label_in_the_containment_panel():
    assert login_session.REFUSED_KIND in security_events.known_kinds()


# -- a real server ---------------------------------------------------------------

def _free_port() -> int:
    with socket.socket() as candidate:
        candidate.bind(("127.0.0.1", 0))
        return int(candidate.getsockname()[1])


def _request(port, method, path, headers=None, body=None):
    conn = http.client.HTTPConnection("127.0.0.1", port, timeout=15)
    try:
        conn.request(method, path, body=body, headers=headers or {})
        response = conn.getresponse()
        data = response.read()
        return response.status, response.getheaders(), data
    finally:
        conn.close()


@pytest.fixture(scope="module")
def server(tmp_path_factory):
    base = tmp_path_factory.mktemp("refused-login")
    root = base / "root"
    root.mkdir()
    env = child_env(base)
    for name in list(env):
        if name.startswith("VSCODE_") or name in (
                "BROWSER", "TERM_PROGRAM", "DISPLAY", "WAYLAND_DISPLAY",
                "DBUS_SESSION_BUS_ADDRESS"):
            env.pop(name)
    env.update({
        "PYTHONPATH": str(REPO),
        "DELFIN_VOILA_ROOT_DIR": str(root),
        "DELFIN_DASHBOARD_STAY_UP": "1",
        "DELFIN_NO_BROWSER": "1",
    })
    port = _free_port()
    log = base / "server.log"
    with log.open("wb") as output:
        proc = subprocess.Popen(
            [sys.executable, "-m", "delfin.cli_voila", "--port", str(port),
             "--ip", "127.0.0.1", "--token", TOKEN, "--no-browser"],
            cwd=str(REPO), env=env, stdout=output, stderr=subprocess.STDOUT,
            start_new_session=True)
    try:
        deadline = time.time() + 90
        while time.time() < deadline:
            if proc.poll() is not None:
                pytest.fail("server exited:\n" + log.read_text(errors="replace")[-2000:])
            try:
                status, _, _ = _request(port, "GET", "/api/status",
                                        {"Authorization": f"token {TOKEN}"})
                if status == 200:
                    break
            except OSError:
                pass
            time.sleep(0.3)
        else:
            pytest.fail("server did not answer")
        yield {"port": port, "log": log, "home": Path(env["HOME"])}
    finally:
        if proc.poll() is None:
            proc.terminate()
            try:
                proc.wait(timeout=20)
            except subprocess.TimeoutExpired:
                proc.kill()
                proc.wait()


def _wait_for(predicate, seconds=10.0):
    deadline = time.time() + seconds
    while time.time() < deadline:
        if predicate():
            return True
        time.sleep(0.2)
    return predicate()


def test_wrong_tokens_reach_the_log_and_the_inbox_once(server):
    port = server["port"]
    for _ in range(20):
        status, _, _ = _request(port, "GET", f"/api/status?token={WRONG}")
        assert status == 403
    _request(port, "GET", "/api/status",
             {"Authorization": f"token {WRONG}"})
    status, headers, _ = _request(port, "GET", "/login")
    xsrf = ""
    for name, value in headers:
        if name.lower() == "set-cookie" and value.startswith("_xsrf="):
            xsrf = value.split(";", 1)[0].split("=", 1)[1]
    body = urllib.parse.urlencode({"password": WRONG, "_xsrf": xsrf})
    status, _, _ = _request(port, "POST", "/login", {
        "Content-Type": "application/x-www-form-urlencoded",
        "Cookie": f"_xsrf={xsrf}"}, body)
    assert status == 401

    log = server["log"]
    assert _wait_for(lambda: "dashboard login refused" in log.read_text(
        errors="replace")), log.read_text(errors="replace")[-2000:]
    text = log.read_text(errors="replace")
    assert text.count("dashboard login refused") == 1, text[-3000:]
    assert WRONG not in text

    inbox = server["home"] / ".delfin" / "attention_inbox.jsonl"
    assert _wait_for(inbox.exists), sorted(
        str(p) for p in server["home"].rglob("*"))[:40]
    events = [json.loads(line) for line in inbox.read_text().splitlines()
              if line.strip()]
    alerts = [e for e in events if e.get("kind") == "security_alert"
              and login_session.REFUSED_KIND in e.get("detail", "")]
    assert len(alerts) == 1, events
    assert WRONG not in json.dumps(events)


def test_the_right_token_is_not_a_refusal(server):
    before = server["log"].read_text(errors="replace").count("login refused")
    status, _, _ = _request(server["port"], "GET", f"/api/status?token={TOKEN}")
    assert status == 200
    after = server["log"].read_text(errors="replace").count("login refused")
    assert after == before
