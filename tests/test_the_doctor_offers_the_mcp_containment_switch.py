"""The MCP row names the switch that contains DELFIN's own servers.

`agent.mcp_isolation = "builtin"` exists and is opt-in by its author's
decision (the roots are derived; it graduates to a default after a real
session). The doctor said "2 without declared roots (outside the shell's
isolation)" and left the user to find the switch. Now the row carries it
as a proposal `/fix` can apply with approval -- only where the loose
servers are the built-ins the switch covers, and only while it is off.
"""

from __future__ import annotations

import pytest

from delfin.agent import doctor as D


@pytest.fixture
def mcp(monkeypatch):
    def _make(configs, *, enabled=False):
        from delfin.agent import mcp_client as MC
        monkeypatch.setattr(MC, "_load_configs", lambda ws: configs)
        from delfin.agent import mcp_isolation as MI
        monkeypatch.setattr(MI, "builtin_isolation_enabled",
                            lambda settings=None: enabled)
        rows = D._check_mcp({"workspace": "", "fast": True})
        assert len(rows) == 1
        return rows[0]
    return _make


def test_loose_builtins_with_the_switch_off_get_the_proposal(mcp):
    row = mcp({"delfin-tools": {"command": "python", "args": ["-m", "x"]},
               "delfin-ops": {"command": "python", "args": ["-m", "y"]}})
    assert row["status"] == "PASS", row
    assert "without declared roots" in row["detail"]
    assert tuple(row.get("setting") or ()) == ("agent.mcp_isolation", "builtin"), row
    assert "mcp_isolation" in row.get("fix", "")


def test_with_the_switch_on_nothing_is_proposed(mcp):
    row = mcp({"delfin-tools": {"command": "python", "args": ["-m", "x"]}},
              enabled=True)
    assert not row.get("setting"), row


def test_a_third_party_server_is_not_offered_the_builtin_switch(mcp):
    """The switch covers DELFIN's own servers; proposing it for somebody
    else's would promise containment the switch does not give."""
    row = mcp({"other": {"command": "node", "args": ["srv.js"]}})
    assert "without declared roots" in row["detail"]
    assert not row.get("setting"), row


def test_a_declared_server_is_not_loose(mcp, tmp_path):
    """A declared root has to exist: a root that is not there contains
    nothing, and the registry rightly reads such an entry as undeclared."""
    # Declared roots sit at the top of the entry, as parse_isolation reads
    # them -- not under an "isolation" key, which is the "off" switch.
    row = mcp({"delfin-tools": {"command": "python", "args": ["-m", "x"],
                                "roots": [str(tmp_path)]}})
    assert "without declared roots" not in row["detail"]
    assert not row.get("setting")


def test_the_row_survives_normalisation():
    """`_normalise` is the trust boundary for every row; a setting it
    rejects would silently become a row with no proposal."""
    row = D._normalise({"check": "mcp servers", "status": "PASS",
                        "detail": "x", "fix": "set it",
                        "setting": ("agent.mcp_isolation", "builtin")}, "mcp")
    assert tuple(row.get("setting") or ()) == ("agent.mcp_isolation", "builtin"), row
