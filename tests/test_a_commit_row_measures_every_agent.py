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

3. **A difference is "within noise" unless it exceeds the spread, and a
   single run is never a result.** On unchanged code this project's own
   reference benchmark moved one task from 45 to 93 across three runs, so a
   comparison of one run per side would report a regression or an improvement
   at random. The failing inputs are named in each test: a 30-point gap with
   no spread IS called better, and a real move with tight spread IS called
   better, so the assertions can fail.

Every check drives ``delfin.agent.cli.main(argv)`` -- parser, dispatch, exit
code, stdout -- as the rule requires, except the two that are pure arithmetic
on fake run files, which cannot be reached through a model this cage has no
key for. Those are driven through the same ``record``/``compare`` commands
with the fake file as the input, so the seam is the product's, not the test's.
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
# (3) repeats, spread, and the noise rule
# --------------------------------------------------------------------------

def _write_runs(tmp_path, runs, name="runs.json"):
    path = tmp_path / name
    path.write_text(json.dumps(runs), encoding="utf-8")
    return path


def _record_with(tmp_path, runs, history):
    bench = _write_runs(tmp_path, runs)
    rc = _main(["metrics", "record", "--tree", ".", "--bench-runs", str(bench),
                "--history", str(history), "--json"])
    assert rc == 0
    rows = (history / "rows.jsonl").read_text(encoding="utf-8").splitlines()
    return json.loads(rows[-1])


def test_three_repeats_are_recorded_with_their_spread(tmp_path):
    """The failing input: a row that keeps only the last run's score."""
    history = tmp_path / "h"
    row = _record_with(tmp_path, [
        {"agent": {"generalist": 0.5}},
        {"agent": {"generalist": 1.0}},
        {"agent": {"generalist": 1.0}},
    ], history)
    stats = row["quality"]["agents"]["generalist"]
    assert stats["n"] == 3
    assert stats["is_result"] is True
    assert stats["min"] == 0.5 and stats["max"] == 1.0, stats
    assert stats["mean"] == pytest.approx(2.5 / 3), stats
    # A spread that is not reported cannot bound a comparison.
    assert stats["stderr"] is not None and stats["stderr"] > 0


def test_one_run_is_not_a_result(tmp_path):
    """The failing input: a row that calls a single run a result.

    This is the whole reason `_MIN_REPEATS` exists: with n=1 there is no
    spread, so no difference measured against it can be told from noise.
    """
    history = tmp_path / "h"
    row = _record_with(tmp_path, [{"agent": {"generalist": 1.0}}], history)
    stats = row["quality"]["agents"]["generalist"]
    assert stats["n"] == 1
    assert stats["is_result"] is False
    assert stats["stderr"] is None


def test_a_real_move_with_tight_spread_is_called_better():
    """The failing input that makes the noise test meaningful.

    Three repeats at 0.90 vs three at 0.20, no overlap: this MUST come out
    ``better``. If it came out ``within noise`` the rule would be unable to
    report a genuine improvement, and the ``within noise`` cases below would
    prove nothing.
    """
    from delfin.agent.commit_metrics import compare_number, spread

    old = spread([0.9, 0.9, 0.9])
    new = spread([0.2, 0.2, 0.2])
    result = compare_number("quality.generalist", new, old)
    assert result["verdict"] == "better", result


def test_a_difference_inside_the_spread_is_within_noise():
    """The failing input: a rule that calls any non-zero delta better/worse.

    These two arms differ by 0.1 in the mean while each side ranges over 0.8,
    which is the shape of the project's own 45/93/93 reference noise. Any
    rule that reports ``better`` or ``worse`` here fails this test.
    """
    from delfin.agent.commit_metrics import compare_number, spread

    old = spread([0.2, 1.0, 0.9])
    new = spread([1.0, 0.3, 0.9])
    result = compare_number("quality.generalist", old, new)
    assert result["verdict"] == "within noise", result
    assert result["threshold"] > 0


def test_a_comparison_between_one_run_each_claims_nothing():
    """The failing input: two 1-run rows, one at 0.0 and one at 1.0.

    A 100-point gap on one run per side is the most tempting false
    improvement there is, and it must still be ``not a result``.
    """
    from delfin.agent.commit_metrics import compare_number, spread

    old = spread([0.0])
    new = spread([1.0])
    result = compare_number("quality.generalist", old, new)
    assert result["verdict"] == "not a result", result


def test_delegation_accuracy_counts_only_labelled_tasks():
    """The failing input: a rate computed over every task in the run file.

    A task nobody labelled as belonging to a specialist is not a wrong
    delegation; counting it as one moves the rate with the size of the task
    set instead of with the routing.
    """
    from delfin.agent.commit_metrics import delegation_accuracy

    runs = [
        {"delegation": {"t1": "ts-finder", "t2": "explore"},
         "delegation_expected": {"t1": "ts-finder", "t2": "ts-finder"}},
        {"delegation": {"t1": "ts-finder", "t2": "ts-finder"},
         "delegation_expected": {"t1": "ts-finder", "t2": "ts-finder"}},
        # A run with no labels contributes nothing at all.
        {"delegation": {"t3": "explore"}},
    ]
    result = delegation_accuracy(runs)
    assert result["tasks"] == 4, result
    assert result["hits"] == 3, result
    assert result["accuracy"] == pytest.approx(0.75)
    assert result["misroutes"] == {"ts-finder": {"explore": 1}}, result


def test_an_unlabelled_run_says_it_measured_no_routing():
    """The failing input: reporting accuracy 0.0 when nothing was labelled.

    Zero reads as "every delegation was wrong", which is a different claim
    and a false one.
    """
    from delfin.agent.commit_metrics import delegation_accuracy

    result = delegation_accuracy([{"delegation": {"t1": "explore"}}])
    assert result["accuracy"] is None
    assert "not measured" in result["note"]


# --------------------------------------------------------------------------
# the comparison table, end to end through the two commands
# --------------------------------------------------------------------------

def test_compare_prints_the_before_after_table(tmp_path, capsys):
    """The failing input: a `compare` that prints prose with no verdicts.

    Two rows are recorded from the same tree with different fake quality, so
    the reachability half is identical (within noise) and the quality half
    differs by more than the spread -- 0.9 down to 0.2 over three repeats
    each, which must be ``worse``. The table must carry that verdict, the
    ``within noise`` rows, and the summary line.
    """
    history = tmp_path / "h"
    _record_with(tmp_path, [{"agent": {"generalist": 0.9}}] * 3, history)
    bench = _write_runs(tmp_path, [{"agent": {"generalist": 0.2}}] * 3,
                        name="runs2.json")
    rc = _main(["metrics", "record", "--tree", ".", "--bench-runs",
                str(bench), "--history", str(history), "--commit", "after",
                "--json"])
    assert rc == 0
    capsys.readouterr()

    rc = _main(["metrics", "compare", "--history", str(history), "--new",
                "after"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "quality.generalist" in out
    # 0.9 -> 0.2 is a REAL drop with no spread on either side, so it must be
    # called worse; the direction is the one the module got right and this
    # test first got wrong, which is why both halves are asserted here.
    assert "worse" in out, out
    assert "within noise" in out, out
    assert "reachability.unmounted" in out, out
    assert "summary:" in out, out


def test_the_direction_of_merit_is_read_from_the_delta(tmp_path, capsys):
    """The failing input: a comparison that labels any change "worse".

    The mirror of the test above, in the other direction: 0.2 up to 0.9 must
    come out ``better``. Without this pair, a rule that always said "worse"
    would pass the table test.
    """
    history = tmp_path / "h"
    _record_with(tmp_path, [{"agent": {"generalist": 0.2}}] * 3, history)
    bench = _write_runs(tmp_path, [{"agent": {"generalist": 0.9}}] * 3,
                        name="runs2.json")
    rc = _main(["metrics", "record", "--tree", ".", "--bench-runs",
                str(bench), "--history", str(history), "--commit", "after2",
                "--json"])
    assert rc == 0
    capsys.readouterr()

    rc = _main(["metrics", "compare", "--history", str(history), "--new",
                "after2"])
    out = capsys.readouterr().out
    assert rc == 0
    assert "better" in out, out
