"""'merk dir das' through the REAL front door: the first remember call wins.

The engine layer has its own test
(test_a_first_remember_succeeds_outside_a_swarm_run.py). This one drives
delfin-agent chat on a pty against the scripted endpoint in
agent_terminal_probes/fake_model, and the scripted model makes the call
the live failure showed (marker do-remember: remember, type=finding,
outside a swarm run). What is asserted:

* the turn reaches idle;
* ONE model request carried the typed phrase -- a second remember call
  to fix the type would be a second request;
* exactly one project memory file exists under the scratch state, with
  the text and the type project in it.
"""

from __future__ import annotations

import json

import pytest

pytest.importorskip("pexpect")
pytest.importorskip("termios")

from agent_terminal_probes import fake_model  # noqa: E402
from agent_terminal_probes.pty_chat import Chat, close_all  # noqa: E402

PHRASE = "do-remember Kaffee ist wichtig, merk dir das."
MEMTEXT = "Der Anwender ist ein schneller Kaffee-Trinker."


def _loopback_reachable(port: int) -> bool:
    """Whether THIS process can connect to 127.0.0.1:<port> at all.

    The chat child inherits the network of the process that spawns it. On
    a host whose sandbox already holds every loopback connect (e.g. an
    egress-caged shell), the child can never reach the endpoint, and the
    failure is environmental, not DELFIN's: every pty test fails the same
    way there (measured 2026-11-13: direct urllib connect -> Permission
    denied; the existing pty suites fail identically on unmodified code).
    """
    import socket
    try:
        with socket.create_connection(("127.0.0.1", port), timeout=3):
            return True
    except OSError:
        return False


def test_the_first_remember_call_ends_the_turn_and_stores_the_memory(
        tmp_path):
    opened: list[Chat] = []
    servers: list = []
    try:
        # The endpoint must exist BEFORE the settings name its port: the
        # egress proxy re-reads agent.sandbox_network.allowed_ports per
        # request, but the port is only knowable once the endpoint is
        # bound. The child's HOME is <root>/home, and load_settings reads
        # <HOME>/.delfin_settings.json (Path.home() is not redirected by
        # DELFIN_SCRATCH_STATE; the memory sinks are, which is where the
        # store is asserted below).
        model_log = tmp_path / "s" / "model.log"
        server, port = fake_model.serve(log_path=str(model_log))
        servers.append(server)
        if not _loopback_reachable(port):
            pytest.skip("this host's sandbox holds loopback connects, so "
                        "the pty child cannot reach a local endpoint")
        home = tmp_path / "s" / "home"
        home.mkdir(parents=True, exist_ok=True)
        (home / ".delfin_settings.json").write_text(json.dumps({
            "agent": {"sandbox_network": {
                "mode": "proxy", "allowed_domains": [],
                "allowed_ports": [port], "upload_domains": []}}}) + "\n")
        # The sandbox shell exports http(s)_proxy (the harness egress
        # proxy). httpx would send the loopback endpoint call through it,
        # which refuses. NO_PROXY makes httpx connect directly; empty
        # proxy vars would NOT help (httpx treats them as a broken proxy,
        # measured: APIConnectionError).
        session = Chat(tmp_path / "s", cols=100, rows=30,
                       args=("--permission-mode", "acceptEdits",
                             "--no-isolate"),
                       server_port=port,
                       env={"FAKE_MODEL_LOG": str(model_log),
                            "no_proxy": "127.0.0.1,localhost",
                            "NO_PROXY": "127.0.0.1,localhost"})
        opened.append(session)
        assert session.wait_for("idle", 60) == "idle", (
            "\n".join(session.visible()))
        assert session.turn(PHRASE) == "idle", (
            "\n".join(session.visible()))
        assert session.model_log.exists(), (
            "the turn reached no endpoint:\n" + "\n".join(session.visible()))
        asked = [json.loads(line) for line
                 in session.model_log.read_text().splitlines()]
        rounds = [r for r in asked if PHRASE in (r.get("last_user") or "")]
        calls = [r for r in rounds if r.get("plan", {}).get("tool") == "remember"]
        # ONE remember call: the first one succeeds, so the tool result does
        # not send the model back for a second try (the request after it is
        # the model's closing text).
        assert len(calls) == 1, (
            "expected one remember call, got "
            + json.dumps([r.get("plan") for r in rounds]))
        assert calls[0]["plan"] == {"tool": "remember",
                                    "args": {"text": MEMTEXT,
                                             "type": "finding"}}, (
            json.dumps(calls[0]["plan"]))
        assert all("tool" not in r.get("plan", {}) for r in rounds[1:]), (
            json.dumps([r.get("plan") for r in rounds]))
        # One project memory under the redirected scratch state, whatever
        # slug the store derives from the scratch workspace.
        stores = [p for p in (tmp_path / "s" / "scratch" / ".delfin")
                  .glob("projects/*/memory/*.md") if p.name != "MEMORY.md"]
        assert len(stores) == 1, (
            "expected exactly one memory file, got: "
            + ", ".join(str(p) for p in stores))
        body = stores[0].read_text(encoding="utf-8")
        assert MEMTEXT in body, body[:400]
        assert "type: project" in body, body[:400]
    finally:
        close_all(opened)
        for server in servers:
            try:
                server.shutdown()
            except Exception:
                pass
