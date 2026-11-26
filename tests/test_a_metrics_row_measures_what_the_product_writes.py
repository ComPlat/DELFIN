"""The metrics row measures what the product writes, and refuses to guess.

Every check here exists because the first version of this feature passed its
own tests while being wrong in a way a reader would never have seen. The
defects were measured on 2026-11-25, each through a real entry point, and
each is named in the test that now holds it:

1. **The quality half had never read a real run.** ``quality_row`` understood
   only an invented ``{"agent": {...}}`` shape that nothing in the product
   writes, so a file from ``delfin-agent bench run`` recorded ``agents: {}``
   -- and ``agents: {}`` reads like "quality was fine". The failing input is
   exactly what ``benchmark.write_run`` produces.
2. **A zero spread was read as certainty.** Three identical repeats are a
   coincidence, not a measurement of no noise; ``0.93/0.93/0.93`` against
   ``0.92/0.92/0.92`` was called ``worse``.
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

import builtins
import json

import pytest

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
    assert "not a result" in out, out


# ---------------------------------------------------------------------------
# (1) the quality half reads a real benchmark run file
# ---------------------------------------------------------------------------

def _write_bench_file(path, qualities, profile="kit.model-a", unmeasured=False):
    """One run file, in the exact shape ``benchmark.write_run`` writes.

    The first line is the ``run_header`` -- provenance, not a result -- and
    every further line is one task's result.
    """
    lines = [json.dumps({"run_header": {"model": profile,
                                        "settings": {"repeats": 1}}})]
    for index, quality in enumerate(qualities):
        lines.append(json.dumps({
            "task_id": f"t{index}",
            "task_class": "chem",
            "model": profile,
            "profile_name": profile,
            "mode": "solo",
            "success": quality >= 50,
            "quality_0_100": quality,
            "unmeasured": unmeasured,
        }))
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def test_one_real_run_file_is_read_at_all(tmp_path):
    """The failing input: ONE file in the shape ``bench run`` writes.

    Sharper than the three-file test above: this one has no comma-path or
    repeat counting in the way, so it fails on the defect itself -- the row
    recorded ``agents: {}`` for a real run file, which reads as "no agent had
    a quality" and as "quality was not bad" at the same time.
    """
    bench = _write_bench_file(tmp_path / "run.jsonl", [90, 80], profile="kit.m")
    history = tmp_path / "h"
    rc = _main(["metrics", "record", "--tree", ".", "--bench-runs", str(bench),
                "--history", str(history), "--json"])
    assert rc == 0, rc
    quality = json.loads(
        (history / "rows.jsonl").read_text().splitlines()[-1])["quality"]
    assert "kit.m" in quality["agents"], quality
    stats = quality["agents"]["kit.m"]
    assert stats["mean"] == pytest.approx(85.0), stats
    # One file is one repeat, and one repeat is not a result.
    assert stats["is_result"] is False, stats


def test_three_real_run_files_are_three_repeats(tmp_path):
    """The failing input: a real run file, which gave ``agents: {}``.

    Three files, because one ``bench run`` invocation is one repeat and the
    noise rule needs three. The check is on the ROW the command records, so
    the count and the file-reading are measured together.
    """
    files = []
    for index, qualities in enumerate(([90, 80], [60, 70], [100, 75])):
        files.append(_write_bench_file(tmp_path / f"run{index}.jsonl",
                                       qualities))
    history = tmp_path / "h"
    rc = _main(["metrics", "record", "--tree", ".", "--bench-runs",
                ",".join(str(f) for f in files), "--history", str(history),
                "--json"])
    assert rc == 0, rc
    row = json.loads((history / "rows.jsonl").read_text().splitlines()[-1])
    quality = row["quality"]
    assert quality["grouped_by"] == "profile_name", quality
    stats = quality["agents"]["kit.model-a"]
    # 85, 65, 87.5 -- one sample per FILE, not one per task.
    assert stats["n"] == 3, stats
    assert stats["is_result"] is True, stats
    assert stats["min"] == pytest.approx(65.0), stats
    assert stats["max"] == pytest.approx(87.5), stats


def test_a_run_file_of_unmeasured_rows_is_not_a_score_of_zero(tmp_path):
    """The failing input: booking a run that never reached the model as 0.

    ``aggregate_replicates`` excludes unmeasured replicates for exactly this
    reason: a zero lowers a baseline permanently and silently.
    """
    files = [_write_bench_file(tmp_path / f"u{i}.jsonl", [90, 80],
                               profile="kit.model-b", unmeasured=True)
             for i in range(3)]
    history = tmp_path / "h"
    rc = _main(["metrics", "record", "--tree", ".", "--bench-runs",
                ",".join(str(f) for f in files), "--history", str(history),
                "--json"])
    assert rc == 0, rc
    quality = json.loads(
        (history / "rows.jsonl").read_text().splitlines()[-1])["quality"]
    assert "kit.model-b" not in quality["agents"], quality
    assert "not a quality of zero" in quality["note"], quality


# ---------------------------------------------------------------------------
# (3) a clean bill of health needs something to have been read
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
# (2) three identical repeats are not certainty
# ---------------------------------------------------------------------------

def test_coincidentally_identical_repeats_do_not_grant_a_direction():
    """The failing input: 0.93 vs 0.92 over three identical repeats each.

    Both sides have ``stderr == 0``. The old rule read that as "any
    difference is outside the spread" and said ``worse`` about a 0.01 move
    that three coincidences cannot distinguish from noise.
    """
    from delfin.agent.commit_metrics import compare_number, spread

    old = spread([0.93, 0.93, 0.93])
    new = spread([0.92, 0.92, 0.92])
    result = compare_number("quality.generalist", old, new)
    assert result["verdict"] == "within noise", result
    assert result["threshold"] == pytest.approx(1.0 / 3.0), result


def test_a_large_move_on_identical_repeats_is_still_called():
    """The mirror: the floor must not swallow a real move.

    Three repeats each, 0.9 down to 0.1. If the floor above were implemented
    as "zero spread means within noise", this would fail and the whole rule
    would be unable to report a regression.
    """
    from delfin.agent.commit_metrics import compare_number, spread

    result = compare_number("quality.generalist",
                            spread([0.9, 0.9, 0.9]), spread([0.1, 0.1, 0.1]))
    assert result["verdict"] == "worse", result


def test_a_row_with_no_spread_is_not_compared():
    """The failing input: a mean and an n, with ``stderr`` absent.

    An absent spread used to read as 0.0 and take the same path as agreement
    between repeats. It is not a measurement of zero noise; it is the absence
    of the measurement.
    """
    from delfin.agent.commit_metrics import compare_number

    no_spread = {"n": 3, "mean": 0.9, "min": 0.9, "max": 0.9, "stderr": None}
    other = {"n": 3, "mean": 0.5, "min": 0.5, "max": 0.5, "stderr": 0.0}
    result = compare_number("quality.generalist", no_spread, other)
    assert result["verdict"] == "not a result", result


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
