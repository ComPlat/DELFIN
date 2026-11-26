"""The metrics row measures what the product writes, and refuses to guess.

Every check here exists because the first version of this feature passed its
own tests while being wrong in a way a reader would never have seen. The
defects were measured on 2026-11-25, each through a real entry point, and
each is named in the test that now holds it:

1-2. (The quality half now reads ``bench_history``; its checks live in
   ``tests/test_a_commit_row_measures_every_agent.py``.)
3. **A root the scan could not read reported a clean bill of health.**
   ``unmounted: 0, within_ceiling: true`` for a directory that does not
   exist, which is a row about nothing that reads as a row about everything.
4. **One tree was measured with two instruments.** The ops-server half was
   imported from the installed package instead of read from the tree being
   measured, so a pull request's own ops-tool change cancelled out between
   the two sides of the comparison.
5. **Ten agents were missing from every row.** ``role_names`` listed the
   delegation presets only, so ``dashboard_agent``, ``office_agent`` and the
   engine's other session roles -- the most restricted roles in the product
   -- had no surface measured at all.
6. **``delfin-agent metrics`` crashed.** A bare call raised
   ``AttributeError: 'Namespace' object has no attribute 'func'``. On a
   machine with no model configured it dies earlier, in ``cli.py``, with
   ``SystemExit`` -- so the red proof of ``test_bare_metrics_prints_its_actions``
   on the base names ``SystemExit``, not the ``AttributeError``. Both are the
   same defect: a bare command that is not a chat prompt. The property this
   suite pins is the machine-independent one -- ``metrics`` reaches its own
   subcommand instead of being routed into a chat turn.
7. **An unreadable source was counted as zero.** The ops-server import and
   all three role sources were wrapped in a bare ``except Exception`` that
   returned an empty set, so a failed import published "0 tools" and a
   shrunken agent count with nothing in the row saying it had not looked.
"""

from __future__ import annotations

import pytest

import builtins


from delfin.agent import cli


def _main(argv):
    return cli.main(argv)


# ---------------------------------------------------------------------------
# (6) the command itself answers when asked with no action
# ---------------------------------------------------------------------------

def test_bare_metrics_prints_its_actions(capsys):
    """The failing input: a bare call, which raised AttributeError.

    A person who types the command to see what it does must get the list,
    not a traceback, and not a 60-second measurement they did not ask for.
    """
    rc = _main(["metrics"])
    out = capsys.readouterr().out
    assert rc == 0, rc
    for action in ("record", "show", "compare"):
        assert action in out, out
    assert "within noise" in out, out


# ---------------------------------------------------------------------------
# (3) an unreadable root is not a count of zero
# ---------------------------------------------------------------------------

def test_an_unreadable_root_is_not_a_count_of_zero():
    """The failing input: a path that does not exist.

    ``unmounted: 0`` with ``within_ceiling: true`` is the most dangerous row
    this command can print, because it answers "nothing is unreachable" about
    a tree it never opened.
    """
    from delfin.agent.commit_metrics import reachability

    row = reachability(__import__("pathlib").Path("/nonexistent/metrics-tree"))
    assert row["unmounted"] is None, row
    assert row["within_ceiling"] is None, row
    assert row["callables_defined"] == 0, row
    assert "not a count of zero" in row["error"], row


# ---------------------------------------------------------------------------
# (4) both halves of the MCP surface come from the tree being measured
# ---------------------------------------------------------------------------

def test_the_ops_half_is_read_from_the_measured_tree(tmp_path):
    """The failing input: a tree whose ops server registers two tools.

    Read from the tree, the answer is those two names. Imported from the
    installed package, it is the 102 this checkout ships -- and then a pull
    request that adds one is counted on both sides and cancels out.
    """
    from delfin.agent.commit_metrics import mcp_tool_names

    server = tmp_path / "delfin" / "ops_server" / "server.py"
    server.parent.mkdir(parents=True)
    server.write_text(
        "def _register(mcp):\n"
        "    mcp.tool(name='only_one')(lambda: 1)\n"
        "    mcp.tool(name=\"only_two\")(lambda: 2)\n"
        "    return mcp\n", encoding="utf-8")
    names = mcp_tool_names(tmp_path)["delfin-ops"]
    assert names == frozenset({"only_one", "only_two"}), names


def test_the_shipped_ops_surface_is_still_read_in_full():
    """The control case: the real tree still reports its whole surface.

    Without this, the test above would pass on a parser that returned the
    two synthetic names because it had stopped finding anything real.
    """
    from pathlib import Path

    from delfin.agent.commit_metrics import mcp_tool_names

    names = mcp_tool_names(Path(__file__).resolve().parents[1])["delfin-ops"]
    assert len(names) > 50, len(names)


def test_the_operating_knowledge_layer_is_counted():
    """The failing input: a glob of ``knowledge/*.md`` only.

    Ten operating notes live one level down and reach an agent's turn, so the
    old count reported that layer as zero -- the same "absent reads as zero"
    defect as the unreadable root, in a smaller place.
    """
    from pathlib import Path

    from delfin.agent.commit_metrics import _packs_row

    row = _packs_row(Path(__file__).resolve().parents[1])
    assert row["knowledge_operating"] > 0, row
    assert row["skills"] > 0 and row["agents"] > 0, row


# ---------------------------------------------------------------------------
# (5) every role that can run a session is in the row
# ---------------------------------------------------------------------------

def test_every_role_the_product_restricts_is_measured():
    """The failing input: the preset registry alone.

    A role with its own allow-list runs a session and is offered a surface
    that differs from the generalist's. Measured 2026-11-25: the presets
    alone left ten of them out, including both of the tightest.
    """
    from delfin.agent.api_client import _ROLE_EXEC_ALLOWLIST
    from delfin.agent.commit_metrics import role_names
    from delfin.agent.engine import _ROLE_TOOL_WHITELIST

    measured = set(role_names())
    missing = (set(_ROLE_EXEC_ALLOWLIST) | set(_ROLE_TOOL_WHITELIST)) - measured
    assert not missing, sorted(missing)


# ---------------------------------------------------------------------------
# (7) an unreadable source is named, never counted as zero
# ---------------------------------------------------------------------------

def _block_ops_import(monkeypatch):
    """Make ``delfin.ops_server.server`` fail to import, for one test.

    The failing input is a broken import, which happens for real when a
    dependency of that server (the MCP SDK) is absent.
    """
    real_import = builtins.__import__

    def _import(name, *a, **kw):
        if name == "delfin.ops_server.server" or (
                name == "delfin.ops_server" and "server" in (a[2] or ())):
            raise ImportError("blocked for this test")
        return real_import(name, *a, **kw)

    monkeypatch.setattr(builtins, "__import__", _import)


# Scans the whole repository (30-100 s): in the slow suite, so a fast
# CI shard stays within its 30-minute limit (2026-11-26).
@pytest.mark.slow
def test_an_unreadable_ops_surface_is_named_not_counted_as_zero(tmp_path):
    """The failing input: a tree with no ops server at all.

    The reachable form of the defect: ``build_row`` always passes a tree, so
    the AST path is the one the product takes. On the base a tree without
    ``delfin/ops_server/server.py`` yielded ``delfin-ops: 0`` -- a number a
    reader compares between commits, published as if it had been measured,
    and indistinguishable from a server that registers no tools.
    """
    from delfin.agent import commit_metrics

    empty_tree = tmp_path / "tree"
    (empty_tree / "delfin" / "agent").mkdir(parents=True)
    row = commit_metrics.build_row(tree_root=empty_tree)
    surface = row.get("mcp_surface")
    if surface is None:
        raise AssertionError(
            "the row cannot report that the ops server was unreadable: a "
            "missing server is published as a server exposing no tools")
    assert surface["complete"] is False, surface
    assert surface["errors"].get("delfin-ops"), surface
    assert "ops server" in surface["errors"]["delfin-ops"], surface


def test_a_readable_ops_surface_reports_no_error():
    """The control case: the shipped tree, whose ops server imports.

    It passes on the base and on the fix, so the test above measures the
    failure path rather than the ordinary one.
    """
    from pathlib import Path

    from delfin.agent import commit_metrics

    names = commit_metrics._ops_server_names(Path(__file__).resolve().parents[1])
    assert names, "the shipped ops server exposes no tools"


# Scans the whole repository (30-100 s): in the slow suite, so a fast
# CI shard stays within its 30-minute limit (2026-11-26).
@pytest.mark.slow
def test_the_agent_count_says_when_a_role_source_was_unreadable(monkeypatch):
    """The failing input: a role source whose import fails.

    The row's agent count is a claim about *every* agent. A source that
    cannot be imported used to shrink that count with nothing in the row
    saying so, so "30 agents" could mean "30 of 41, and eleven had no
    import". Asserted on the ROW -- the object a reader compares between
    commits -- without touching a helper's name on the base.
    """
    from pathlib import Path

    from delfin.agent import commit_metrics

    def _boom():
        raise ImportError("engine is not importable here")

    sources_fn = getattr(commit_metrics, "_role_sources", None)
    if sources_fn is None:
        raise AssertionError(
            "no way to ask which role sources this row can read: the row "
            "cannot say that a source was unreadable")
    monkeypatch.setattr(commit_metrics, "_role_sources", lambda: [
        {"source": "a source that works", "names": lambda: ["alpha"]},
        {"source": "a broken source", "names": _boom},
    ])
    row = commit_metrics.build_row(tree_root=Path(__file__).resolve().parents[1])
    sources = row.get("agent_sources")
    assert sources is not None, "the row does not carry agent_sources"
    assert sources["complete"] is False, sources
    assert sources["read"] == {"a source that works": 1}, sources
    assert "ImportError" in sources["unreadable"]["a broken source"], sources


def test_every_role_source_of_the_shipped_tree_is_readable():
    """The control case for the field above, on the real tree.

    Passes on both versions; it fails if a rename makes a source silently
    unimportable, which is the defect the field exists to surface.
    """
    from delfin.agent import commit_metrics

    row = commit_metrics.role_sources_row()
    assert row["complete"] is True, row["unreadable"]
    assert row["read"], row
