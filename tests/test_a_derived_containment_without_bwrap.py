"""What happens to a DERIVED containment when the host has no bwrap.

Measured first (see the first case): with the derived isolation switched
on and bubblewrap unusable, today's start() refuses the server outright
-- the same refusal a user-written containment gets. On a host without
bwrap or without user namespaces (macOS, some clusters) that takes out
``delfin-tools`` and ``delfin-ops`` entirely, and delfin-ops carries the
prescribed first step of every ORCA question.

The decision these cases pin (task 9; the DEFAULT itself is not flipped
here -- ``agent.mcp_isolation`` stays opt-in):

* a DECLARED containment still refuses. The user wrote the containment
  down; running without it is not a lesser version of what they asked
  for. Unchanged, and still covered by
  test_a_server_is_contained_at_its_launch.py.
* a DERIVED containment falls back to the uncontained start with a
  visible note. The roots were inferred, not asked for, and the module
  already settled this class of question for absent directories: "an
  inference that stops the server would be the derivation making
  policy" (delfin_roots). On a bwrap-less host the uncontained start is
  also not a new hole -- it is exactly what runs there today with the
  setting off; the fallback only keeps parity with the status quo while
  the note says so out loud. The listing follows: a server that will
  start uncontained is shown with no containment, so /mcp and the
  banner keep telling the truth.
"""
import pytest

from delfin.agent import mcp_client, mcp_isolation


def _derived_iso(tmp_path):
    home = tmp_path / "home"
    home.mkdir()
    return mcp_isolation.delfin_roots(settings={}, home=home,
                                      workspace=tmp_path)


def _builtin_server(iso):
    name, cfg = sorted(mcp_client._BUILTIN_SERVERS.items())[0]
    return mcp_client.MCPServer(
        name=name, command=cfg["command"], args=list(cfg["args"]),
        isolation=iso)


def test_the_derived_roots_are_marked_as_derived(tmp_path):
    """The distinction this decision turns on: roots DELFIN inferred vs
    roots the user wrote down. ``parse_isolation`` speaks for the user,
    ``delfin_roots`` for the machine."""
    iso = _derived_iso(tmp_path)
    assert iso is not None and iso.derived is True

    root = tmp_path / "project"
    root.mkdir()
    declared = mcp_isolation.parse_isolation({"roots": [str(root)]})
    assert declared is not None and declared.derived is False


def test_a_derived_server_falls_back_to_the_uncontained_start(
        monkeypatch, tmp_path):
    """The chosen behaviour, and the measurement it replaces. Measured
    on the code BEFORE it (see the commit message): with bwrap unusable
    the derived isolation refused the start exactly like a declared
    one -- on a host without bwrap that would take out both built-in
    servers, and delfin-ops carries the prescribed first step of every
    ORCA question. The old decision was the module's own logic applied
    one step too far: delfin_roots already drops (not refuses) a
    directory the user does not keep, because an inference must not
    make policy. The same reasoning carries the no-bwrap host: the
    derived roots were inferred, not asked for, and on a host that
    cannot honour them the server starts as it always ran there --
    which is not a new hole, it is exactly what runs on such a host
    today with the setting off. The note keeps the fallback from being
    a silent one."""
    monkeypatch.setattr(mcp_isolation, "bwrap_functional", lambda: False)
    server = _builtin_server(_derived_iso(tmp_path))
    server.start()

    assert server.proc is not None, server.last_error
    assert server.containment_note, "the fallback was silent"
    assert "bubblewrap" in server.containment_note
    assert "uncontained" in server.containment_note
    server.stop()


def test_a_declared_server_still_refuses_without_bwrap(
        monkeypatch, tmp_path):
    """No regression on the user's own words: a declared containment
    refuses, as it always has. The fallback is for what DELFIN inferred,
    not for what someone wrote down."""
    monkeypatch.setattr(mcp_isolation, "bwrap_functional", lambda: False)
    root = tmp_path / "project"
    root.mkdir()
    server = mcp_client.MCPServer(
        name="docs", command="/bin/echo", args=["hi"],
        isolation=mcp_isolation.parse_isolation({"roots": [str(root)]}))
    server.start()

    assert server.proc is None
    assert server.containment_note == ""
    assert "bubblewrap" in server.last_error


def test_the_listing_admits_the_fallback(tmp_path, monkeypatch):
    """A server that will start uncontained must not be listed as
    contained -- /mcp and the banner read this row, and a row that says
    "(rw)" about a process with no namespace is the config-file shape of
    the defect this module exists to close."""
    monkeypatch.setattr(mcp_isolation, "bwrap_functional", lambda: False)
    monkeypatch.setattr(mcp_isolation, "builtin_isolation_enabled",
                        lambda *a, **k: True)
    monkeypatch.setattr(mcp_isolation, "delfin_roots",
                        lambda **k: mcp_isolation.Isolation(
                            ("/somewhere", ()), (), derived=True))
    rows = {r["name"]: r for r in mcp_client.effective_servers(None)}
    for name in mcp_client._BUILTIN_SERVERS:
        assert rows[name]["isolation"] == "", rows[name]
