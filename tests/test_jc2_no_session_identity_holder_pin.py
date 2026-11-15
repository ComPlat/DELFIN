"""Pin: nothing in delfin.tools.* consults a module-global for session
identity (deny-pin, both directions).

The session id reaches a chain via the per-request `_meta` of the
tools/call request (JC2 target 2, operator route).  An earlier design
used a module-global holder on session_store; it was stripped because a
dormant identity holder is exactly the dead code a later change
mistakes for a working channel.  This pin fires if a future commit
reintroduces that pattern into the tools/agent modules this package
owns.

Control: committed WITH the strip — it is green after the strip and
would have been red against ca8b2f5ff (the holder existed then).
"""

from __future__ import annotations

import ast

_FORBIDDEN_HOLDER_NAMES = {
    "current_agent_session_id",
    "note_agent_session_id",
}

_MODULES = [
    "delfin.tools.chain_state",
    "delfin.tools.chain_index",
    "delfin.tools.pipeline",
    "delfin.cli_chain",
    "delfin.tools.platform",
    "delfin.tools.mcp_server",
]


def _source(module_name: str) -> str:
    import delfin
    from pathlib import Path

    root = Path(delfin.__path__[0])
    rel = module_name.replace(".", "/") + ".py"
    return (root.parent / rel).read_text(encoding="utf-8")


def test_no_tools_module_defines_a_session_identity_holder():
    for mod_name in _MODULES:
        tree = ast.parse(_source(mod_name))
        for node in ast.walk(tree):
            names = []
            if isinstance(node, ast.Name):
                names = [node.id]
            elif isinstance(node, ast.Attribute):
                names = [node.attr]
            assert not (_FORBIDDEN_HOLDER_NAMES & set(names)), (
                f"{mod_name} consults or defines a session-identity "
                f"module-global ({sorted(set(names) & _FORBIDDEN_HOLDER_NAMES)})"
                " — the session id travels per request in _meta "
                "(delfin/agent_session_id); a global holder is the dead "
                "channel pattern this package stripped")


def test_chain_index_records_only_explicit_session_ids():
    """The one legitimate session-id source in the tools tree is the
    explicit kwarg/`_meta` parameter: the index records what the caller
    passed and nothing else — no environment discovery."""
    from delfin.tools import chain_index

    import inspect

    src = inspect.getsource(chain_index.record)
    assert "os.environ" not in src and "environ" not in src, (
        "chain_index.record must never discover the session id from "
        "the environment")
