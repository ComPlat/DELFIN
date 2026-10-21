"""The agent service behind ``delfin-agent tui``.

The terminal UI is a TypeScript program; the agent stays here, in Python.
This package is everything on the Python side of that boundary:

* ``protocol``  -- the wire format both sides speak (JSON-RPC 2.0, one
  JSON object per line) and the single definition of every message;
* ``worker``    -- one agent session driven over that protocol;
* ``approvals`` -- the session's approval questions, answered remotely;
* ``transcript`` -- what a turn shows, as structured events;
* ``eventlog``  -- the durable record a window replays when it attaches;
* ``fake_engine`` -- scripted turns for tests and UI work, no model.

docs/TUI_REWRITE.md has the design and docs/TUI_PROTOCOL.md the messages.
"""
