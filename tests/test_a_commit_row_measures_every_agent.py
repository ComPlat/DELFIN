"""A commit row measures every agent, and calls a difference noise only when it is.

Three claims are pinned here, each with the input that would make it fail
named in the test that makes it:

1. **The row names every agent separately.** One count for "the agent"
   answers a question about no agent that runs: the generalist is offered 97
   tools and ``workflow-operator`` 16, and they are not interchangeable. The
   failing input is a row that reports only the generalist, or that reports a
   specialist as having the generalist's surface.

2. **The offered count includes the MCP servers.** Measured on 41e15fe9f,
   ``api_client.role_tool_surface_report`` reports ``xyz-editor`` as 5 tools
   because it reads the built-in catalogue only; the session offers it 23,
   the other 18 arriving from ``delfin-ops``. The failing input is a count
   taken from the catalogue alone -- which is exactly what the older function
   returns, so this asserts the two numbers differ rather than merely exist.

3. **Quality is read from the kept benchmark rounds, per capability tag.**
   The rounds are written through ``bench_history``'s own writers and read
   back by ``metrics show``/``compare``; a difference is ``within noise``
   unless ``bench_ab.newcombe_interval`` excludes zero, and rounds of
   different models are not compared.

Every check drives ``delfin.agent.cli.main(argv)`` -- parser, dispatch, exit
code, stdout.
"""

from __future__ import annotations

import json

import pytest

from delfin.agent import cli


def _main(argv):
    return cli.main(argv)


@pytest.fixture(scope="module")
def recorded(tmp_path_factory):
    """One real row, measured once through the CLI, for the whole module.

    The scan walks the whole tree twice (product and suite), so measuring it
    per test would pay the same 30 s several times for the same answer.
    """
    history = tmp_path_factory.mktemp("metrics_history")
    rc = _main(["metrics", "record", "--tree", ".",
                "--history", str(history), "--json"])
    assert rc == 0, "metrics record failed"
    rows = (history / "rows.jsonl").read_text(encoding="utf-8").splitlines()
    assert len(rows) == 1, rows
    return json.loads(rows[0]), history


# --------------------------------------------------------------------------
# (1) every agent, separately
# --------------------------------------------------------------------------

# Scans the whole repository (30-100 s): in the slow suite, so a fast
# CI shard stays within its 30-minute limit (2026-11-26).
@pytest.mark.slow
def test_the_row_carries_one_surface_per_agent(recorded):
    """The failing input: a row whose `agents` holds only the generalist."""
    row, _ = recorded
    agents = row["agents"]
    assert "(generalist)" in agents, sorted(agents)
    specialists = [a for a in agents if a != "(generalist)"]
    # The preset registry ships this many; a row that measured only the
    # generalist would report zero here and fail.
    assert len(specialists) >= 10, sorted(agents)


def test_each_agent_carries_its_own_numbers(recorded):
    """The failing input: a row that copies the generalist's surface onto
    every specialist, which is what a single shared count does."""
    row, _ = recorded
    agents = row["agents"]
    generalist = agents["(generalist)"]
    narrow = agents["workflow-operator"]
    assert narrow["count"] != generalist["count"], (
        "workflow-operator is reported with the generalist's surface, so the "
        "row did not measure it separately")
    # Both halves are always reported, so a reader can see which moved.
    for entry in agents.values():
        assert "catalogue_count" in entry and "mcp_count" in entry
        assert entry["count"] == entry["catalogue_count"] + entry["mcp_count"]


def test_a_specialist_gets_its_prompt_size_and_it_differs(recorded):
    """The failing input: a row that reports the generalist's prompt size for
    every agent -- the specialist prompts are a fifth of its length."""
    row, _ = recorded
    agents = row["agents"]
    generalist = agents["(generalist)"]["prompt_chars"]
    reviewer = agents["chemistry-reviewer"]["prompt_chars"]
    assert reviewer < generalist, (reviewer, generalist)


# --------------------------------------------------------------------------
# (2) the offered count includes what the MCP servers expose
# --------------------------------------------------------------------------

def test_the_offered_count_reaches_past_the_builtin_catalogue(recorded):
    """The failing input: a count taken from the catalogue alone.

    The comparison is against the older function's own answer, measured on
    this tree, so the test fails if either number drifts into the other --
    and it fails if the MCP half is dropped, which is the defect it exists
    for.
    """
    from delfin.agent.api_client import role_tool_surface_report

    row, _ = recorded
    editor = row["agents"]["xyz-editor"]
    catalogue_only = role_tool_surface_report(["xyz-editor"])["xyz-editor"]
    assert editor["catalogue_count"] == catalogue_only["count"], (
        "the catalogue half no longer agrees with the product's own report")
    assert editor["mcp_count"] > 0, (
        "xyz-editor declares 18 geometry tools and all of them live on the "
        "delfin-ops server; an mcp_count of 0 means the MCP half was lost")
    assert editor["count"] > catalogue_only["count"], (
        "the offered count equals the catalogue-only count, so the row would "
        "report the structure specialists as offered almost nothing")


def test_a_declared_tool_is_reachable_somewhere(recorded):
    """The failing input: a row whose MCP half lists names no server exposes.

    Pinned because a count is not a capability: every MCP name in the row
    must be one the server module actually registers.
    """
    from delfin.agent.commit_metrics import mcp_tool_names
    from pathlib import Path

    servers = mcp_tool_names(Path("."))
    assert servers["delfin-tools"], "the tools server registered nothing"
    assert servers["delfin-ops"], "the ops server exposed nothing"
    row, _ = recorded
    known = {f"mcp__{s}__{n}" for s, names in servers.items() for n in names}
    for agent, entry in row["agents"].items():
        unknown = [n for n in entry["mcp_names"] if n not in known]
        assert not unknown, (agent, unknown)


# --------------------------------------------------------------------------
# (3) quality comes from the kept benchmark rounds, per capability tag
# --------------------------------------------------------------------------

def _keep_round(root, round_id, commit, outcomes, *, model="m", started=1.0):
    """One round through ``bench_history``'s own writers.

    *outcomes* is ``[(capability, success), ...]``; the capability goes into
    ``task_class``, the field the benchmark writer stamps.
    """
    from delfin.agent import bench_history

    bench_history.write_round(
        {"round_id": round_id, "commit": commit, "model": model,
         "seed_mode": "practice", "task_set_version": "v1",
         "started_at": started},
        history_root=root)
    bench_history.write_instances(
        round_id,
        [{"task_id": f"t{i}", "task_class": cap, "ran": True,
          "unmeasured": False, "success": ok}
         for i, (cap, ok) in enumerate(outcomes)],
        history_root=root, started_at=started)


def _two_rows(history, old="aaaaaaa", new="bbbbbbb"):
    history.mkdir(parents=True, exist_ok=True)
    with (history / "rows.jsonl").open("w", encoding="utf-8") as fh:
        for commit, unmounted in ((old, 9), (new, 8)):
            fh.write(json.dumps({"commit": commit, "reachability":
                                 {"unmounted": unmounted}}) + "\n")


def test_compare_reads_kept_rounds_per_capability(tmp_path, capsys):
    """The failing input: a clear drop in one capability, flat in another.

    Three rounds per commit; ``chem`` goes from 10/10 to 0/10 per round and
    ``code`` stays 5/10. The first must be ``worse``, the second
    ``within noise``, and the reachability count 9 -> 8 ``better``.
    """
    bench = tmp_path / "bench"
    for k in range(3):
        _keep_round(bench, f"o{k}", "aaaaaaa",
                    [("chem", True)] * 10 + [("code", i < 5) for i in range(10)],
                    started=float(k + 1))
        _keep_round(bench, f"n{k}", "bbbbbbb",
                    [("chem", False)] * 10 + [("code", i < 5) for i in range(10)],
                    started=float(k + 10))
    history = tmp_path / "rows"
    _two_rows(history)
    rc = _main(["metrics", "compare", "--old", "aaaaaaa", "--new", "bbbbbbb",
                "--history", str(history), "--bench-history", str(bench),
                "--json"])
    assert rc == 0
    result = json.loads(capsys.readouterr().out)
    by = {c["number"]: c for c in result["comparisons"]}
    assert by["reachability.unmounted"]["verdict"] == "better"
    assert by["quality.capability:chem"]["verdict"] == "worse", by
    assert by["quality.capability:code"]["verdict"] == "within noise", by
    assert by["quality.all"]["old_rounds"] == 3


def test_one_task_per_side_claims_nothing(tmp_path, capsys):
    """The failing input: 1/1 against 0/1 -- the means moved by 1.0."""
    bench = tmp_path / "bench"
    _keep_round(bench, "o", "aaaaaaa", [("chem", True)])
    _keep_round(bench, "n", "bbbbbbb", [("chem", False)], started=2.0)
    history = tmp_path / "rows"
    _two_rows(history)
    rc = _main(["metrics", "compare", "--old", "aaaaaaa", "--new", "bbbbbbb",
                "--history", str(history), "--bench-history", str(bench),
                "--json"])
    assert rc == 0
    by = {c["number"]: c
          for c in json.loads(capsys.readouterr().out)["comparisons"]}
    assert by["quality.all"]["verdict"] == "within noise", by


def test_rounds_of_different_models_are_not_compared(tmp_path, capsys):
    bench = tmp_path / "bench"
    _keep_round(bench, "o", "aaaaaaa", [("chem", True)] * 10, model="m1")
    _keep_round(bench, "n", "bbbbbbb", [("chem", False)] * 10, model="m2",
                started=2.0)
    history = tmp_path / "rows"
    _two_rows(history)
    rc = _main(["metrics", "compare", "--old", "aaaaaaa", "--new", "bbbbbbb",
                "--history", str(history), "--bench-history", str(bench),
                "--json"])
    assert rc == 0
    by = {c["number"]: c
          for c in json.loads(capsys.readouterr().out)["comparisons"]}
    assert by["quality.all"]["verdict"] == "not comparable", by


def test_show_carries_per_capability_quality_with_spread(tmp_path, capsys):
    """``show`` reads the same rounds, and reports the per-round spread."""
    bench = tmp_path / "bench"
    _keep_round(bench, "a", "aaaaaaa", [("chem", True), ("chem", True)])
    _keep_round(bench, "b", "aaaaaaa", [("chem", True), ("chem", False)],
                started=2.0)
    history = tmp_path / "rows"
    _two_rows(history)
    rc = _main(["metrics", "show", "--commit", "aaaaaaa", "--history",
                str(history), "--bench-history", str(bench), "--json"])
    assert rc == 0
    quality = json.loads(capsys.readouterr().out)["quality"]
    chem = quality["numbers"]["capability:chem"]
    assert quality["rounds"] == 2
    assert (chem["passed"], chem["n"], chem["min"], chem["max"]) == (3, 4, 0.5, 1.0)
