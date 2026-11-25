"""A real run says whether it touched the change it was meant to measure.

Input: ``realrun.run(run, tid)`` on a toy repository whose ``delfin/agent/cli``
is a stand-in for the agent's one-shot entry point (the same module and the
same command line the product uses, so the probe child is started exactly the
way a member's real run is started) and whose ``delfin/agent/inner`` holds one
function the branch edits. Two probes: one whose prompt makes the child call
the edited function, one that calls a different one. No model, no network, no
provider key: the child answers for itself.

Output pinned, per probe and side: the tools the model was OFFERED in its first
turn, the tools it CALLED, and ``reached the change`` -- yes / no / unknown --
written into ``scratch/realrun/summary.md``, plus the stage's resume note that
says plainly when no probe reached the change.

Why it exists: on 2026-11-25 three probes of one package never executed the
changed code, and the report found that out only by reading transcripts. A run
that cannot see the change measured nothing, and nothing in the summary said
so. The same defect could not be caught by a test that checks the column exists
alone, so each of the three states the column can hold is produced here by a
run that really ends in it -- an ``unknown`` by a child that dies before its
observer can report, a ``no`` by a probe that answers without entering the
edited function, a ``yes`` by one that enters it.

The observer mechanism is measured, not assumed: the offered tool list is read
from the request frame, and DELFIN issues that request from a thread it starts
for the purpose, so a profile installed with ``sys.setprofile`` alone would see
nothing and every probe would answer "no tools were offered". The child below
reproduces that shape (``_spawn_create``) for the same reason.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

from delfin.agent.swarm import realrun, stages, wave

# The toy repository: a base commit, and a branch that edits ONE line inside
# ``inner.touched``. ``inner.untouched`` is a decoy of the same shape, so a
# probe that calls it enters the package and still misses the change -- which
# is the case that used to look like evidence.
_INNER_BASE = '''"""The package a probe either reaches or does not."""


def touched():
    return 1


def untouched():
    return 2


class Holder:
    def nested(self):
        def helper():
            return 3
        return helper()
'''

_INNER_BRANCH = '''"""The package a probe either reaches or does not."""


def touched():
    return 99   # the one change on the branch


def untouched():
    return 2


class Holder:
    def nested(self):
        def helper():
            return 3
        return helper()
'''

#: The child DELFIN's real-run code starts: ``python -P -m delfin.agent.cli
#: chat -p ...``. It answers for itself, and it is written to fail loudly if the
#: observer changes what the run does -- the observer must never be able to
#: steer a probe.
_CHILD = r'''
import json, os, sys, threading

TOOLS = [{"function": {"name": "bash"}}, {"function": {"name": "read_file"}},
         {"function": {"name": "write_file"}}]


class _Completions:
    def create(self, **kwargs):
        # The provider client's request: this frame is where the whole offered
        # tool list exists, exactly as in the real OpenAIClient.
        return "resp"


class _Chat:
    completions = _Completions()


class Client:
    chat = _Chat()


def emit(obj):
    print(json.dumps(obj), flush=True)


def main(prompt):
    if prompt.startswith("die"):
        os._exit(3)                      # before the observer can report
    from delfin.agent import inner
    called = []
    if "touched" in prompt and not prompt.startswith("do-untouched"):
        inner.touched(); called.append("bash")
    else:
        inner.untouched(); called.append("read_file")
    client = Client()
    kwargs = {"model": "toy-model", "messages": [{"role": "user", "content": prompt}],
              "tools": TOOLS}
    box = []
    # DELFIN issues its request from a thread started for that purpose; the
    # profile has to reach it, so this shape is not decoration.
    t = threading.Thread(target=lambda: box.append(client.chat.completions.create(**kwargs)))
    t.start(); t.join()
    for name in called:
        emit({"type": "tool_use", "name": name, "input": {}})
    emit({"type": "result", "stop_reason": "done", "input_tokens": 11, "output_tokens": 5,
          "tool_calls": [{"name": n} for n in called], "text": "toy answer: " + prompt})


if __name__ == "__main__":
    a = sys.argv[1:]
    main(a[a.index("-p") + 1] if "-p" in a else "")
'''


def _dead_pid() -> int:
    """A pid no process owns: a member whose session has ended.

    ``advance`` refuses to touch a member whose pid is alive, so a fixture that
    means "this session is over" cannot borrow the test's own pid -- it would
    make the stage do nothing, and a test that asserts on the note it never
    received would pass by saying nothing. A number no pid can reach is used
    rather than a reaped child's pid, which the kernel may hand to some later
    process and turn the same fixture back into a live member.
    """
    return 2 ** 30


def _git(repo: Path, *args) -> None:
    subprocess.run(["git", *args], cwd=str(repo), check=True,
                   env=dict(os.environ, GIT_AUTHOR_NAME="t", GIT_AUTHOR_EMAIL="t@example.org",
                            GIT_COMMITTER_NAME="t", GIT_COMMITTER_EMAIL="t@example.org"))


def _toy_run(tmp_path, monkeypatch, probes, change_body=_INNER_BRANCH):
    """A wave whose member's branch edits ``inner.touched``, and its real runs.

    Returns ``(result, summary_text)``. The member is never started: its
    worktree and its ``scratch/realrun.json`` are written directly, which is
    all ``realrun.run`` reads. ``change_body=None`` changes only a text file,
    for the case where the branch holds no Python change at all -- a Python
    comment is NOT that case, because editing module level changes the code the
    import runs.
    """
    repo = tmp_path / "repo"
    (repo / "delfin" / "agent").mkdir(parents=True)
    (repo / "delfin" / "__init__.py").write_text("")
    (repo / "delfin" / "agent" / "__init__.py").write_text("")
    (repo / "delfin" / "agent" / "cli.py").write_text(_CHILD)
    (repo / "delfin" / "agent" / "inner.py").write_text(_INNER_BASE)
    (repo / "README").write_text("toy\n")
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-qm", "base")

    wt = tmp_path / "wt"
    _git(repo.parent / "repo", "clone", "-q", "--branch", "main", str(repo), str(wt))
    _git(wt, "checkout", "-q", "-b", "change")
    if change_body is None:
        (wt / "README").write_text("toy\n\none line the branch added\n")
    else:
        (wt / "delfin" / "agent" / "inner.py").write_text(change_body)
    _git(wt, "add", "-A")
    _git(wt, "commit", "-qm", "edit inner.touched")
    (wt / "scratch").mkdir(exist_ok=True)
    (wt / "scratch" / "realrun.json").write_text(json.dumps(
        {"probes": probes, "provider": "ollama", "model": "toy-model"}))

    root = tmp_path / "waves"
    run = "toy"
    d = root / run
    d.mkdir(parents=True)
    cfg = dict(wave._DEFAULTS, waves_root=str(root), min_available_gb=0, max_load_per_core=1e9)
    monkeypatch.setattr(wave, "settings", lambda: dict(cfg))
    monkeypatch.setattr(wave, "_scoped", lambda cmd, s: cmd)
    (d / "state.json").write_text(json.dumps({"repo": str(repo), "agents": {
        "a1": {"worktree": str(wt), "log": str(d / "a1.log"), "pid": os.getpid(),
               "task": {"id": "a1", "brief": "b", "base": "main", "provider": "ollama",
                        "model": "toy-model"}}}}))
    res = realrun.run(run, "a1")
    return res, (wt / "scratch" / "realrun" / "summary.md").read_text()


_PROBES = [{"name": "onpath", "prompt": "do-touched", "turns": 3},
           {"name": "offpath", "prompt": "do-untouched", "turns": 3}]


def _row(res, probe, side):
    for r in res["rows"]:
        if r["probe"] == probe and r["side"] == side:
            return r
    raise AssertionError(f"no row for {probe}/{side}: "
                         + str([(x["probe"], x["side"]) for x in res["rows"]]))


# -- the real run, end to end ------------------------------------------------------

def test_the_summary_says_which_probe_reached_the_change(tmp_path, monkeypatch):
    """The column, and the two opposite answers in it.

    Red on the base for the reason it names: there is no ``reached the change``
    column, and no row says anything about the changed function.
    """
    res, text = _toy_run(tmp_path, monkeypatch, _PROBES)
    assert "reached the change" in text, text[:1500]
    assert _row(res, "onpath", "branch")["reached"] == "yes"
    assert _row(res, "offpath", "branch")["reached"] == "no"
    # The word in the table is the word in the row, per probe: one aggregate
    # for the member would hide that half the probes measured nothing.
    for probe, word in (("onpath", "yes"), ("offpath", "no")):
        line = [l for l in text.splitlines() if l.startswith(f"| {probe} | branch |")]
        assert len(line) == 1, line
        assert f"| {word} |" in line[0], line[0]
    assert res["reached"] == "yes"


def test_a_probe_that_reached_it_names_the_changed_function(tmp_path, monkeypatch):
    """``yes`` says WHICH function it ran -- a bare yes proves nothing."""
    res, text = _toy_run(tmp_path, monkeypatch, _PROBES)
    hit = _row(res, "onpath", "branch")
    assert "delfin.agent.inner.touched" in hit["reached_reason"], hit["reached_reason"]
    miss = _row(res, "offpath", "branch")
    assert "delfin.agent.inner.touched" not in miss["reached_reason"]
    assert "delfin.agent.inner.touched" in text


def test_the_offered_tools_are_the_ones_the_request_carried(tmp_path, monkeypatch):
    """What the model was GIVEN in its first turn.

    The request goes out from a thread; without ``threading.setprofile`` this
    reads as an empty list, so the assertion below is the measurement that the
    observer reaches that thread. The child offers three tools and calls one.
    """
    res, text = _toy_run(tmp_path, monkeypatch, _PROBES)
    row = _row(res, "onpath", "branch")
    assert row["offered"] == ["bash", "read_file", "write_file"], row["offered"]
    assert row["called_tools"] == ["bash"], row["called_tools"]
    assert _row(res, "offpath", "branch")["called_tools"] == ["read_file"]


def test_a_run_that_died_early_is_unknown_and_not_no(tmp_path, monkeypatch):
    """The third state, reached by a child that dies before it can report.

    ``no`` claims the change was not reached, which a run that never got that
    far cannot support.
    """
    res, text = _toy_run(tmp_path, monkeypatch, [{"name": "dead", "prompt": "die now",
                                                  "turns": 2}])
    row = _row(res, "dead", "branch")
    assert row["reached"] == "unknown", row
    assert "observer wrote nothing" in row["reached_reason"]
    assert res["reached"] == "unknown"
    assert "| unknown |" in text


def test_a_base_with_no_python_change_is_unknown(tmp_path, monkeypatch):
    """A diff that changes no Python file cannot be reached or missed.

    Control-shaped: it names the reason instead of answering ``no``, because
    ``no`` would tell the member its probe missed something that does not
    exist. A Python comment is not this case -- editing module level changes the
    code the import runs, which is what the last test in this file pins.
    """
    res, text = _toy_run(tmp_path, monkeypatch, _PROBES, change_body=None)
    assert res["changed_functions"] == 0, res["changed_functions"]
    assert _row(res, "onpath", "branch")["reached"] == "unknown"
    assert "no Python" in _row(res, "onpath", "branch")["reached_reason"]


def test_the_observer_does_not_change_the_run_it_watches(tmp_path, monkeypatch):
    """Control on both sides of the change: watching is not steering.

    The stop reason, the seconds, and the final answer are what they were
    before the column existed -- the pre-existing columns still carry the run.
    """
    res, text = _toy_run(tmp_path, monkeypatch, _PROBES)
    row = _row(res, "onpath", "branch")
    assert row["stop"] == "done" and row["tools"] == 1
    assert row["in"] == 11 and row["out"] == 5
    assert "toy answer: do-touched" in text
    assert row["error"] in ("", None)


# -- what the diff says changed ----------------------------------------------------

def _toy_repo(tmp_path, base, head):
    """A two-commit repo whose second commit turns ``base`` into ``head``.

    Same content on both sides is allowed and leaves the branch with no Python
    change to find: git refuses an empty commit, so the commit is skipped
    rather than faked with an unrelated file.
    """
    repo = tmp_path / "diffrepo"
    (repo / "pkg").mkdir(parents=True)
    (repo / "pkg" / "__init__.py").write_text("")
    (repo / "pkg" / "mod.py").write_text(base)
    (repo / "notes.txt").write_text("x = 1\n")
    _git(repo, "init", "-q", "-b", "main")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-qm", "base")
    _git(repo, "checkout", "-q", "-b", "change")
    (repo / "pkg" / "mod.py").write_text(head)
    _git(repo, "add", "-A")
    if subprocess.run(["git", "diff", "--cached", "--quiet"], cwd=str(repo)).returncode != 0:
        _git(repo, "commit", "-qm", "change")
    return repo


_BASE = ('VALUE = 1\n\n\ndef helper():\n    return 1\n\n\nclass C:\n    def m(self):\n'
         '        def deep():\n            return 2\n        return deep()\n')
_HEAD_M = _BASE.replace("    return 1\n", "    return 2\n")
_HEAD_NESTED = _BASE.replace("            return 2\n", "            return 3\n")
_HEAD_MODULE = _BASE.replace("VALUE = 1\n", "VALUE = 2\n")


def test_the_diff_names_the_changed_function_by_qualified_name(tmp_path):
    repo = _toy_repo(tmp_path, _BASE, _HEAD_M)
    names = set(realrun.changed_functions(repo, "main", "HEAD")["names"])
    assert names == {"pkg.mod.helper"}, names


def test_a_change_inside_a_nested_function_names_that_function(tmp_path):
    """The qualified name must be the one the trace reports, or nothing matches.

    ``co_qualname`` joins a function nested in a function with ``.<locals>.``,
    which is measured here rather than assumed: a name built with a plain dot
    never matches a frame and every probe would answer ``no``.
    """
    repo = _toy_repo(tmp_path, _BASE, _HEAD_NESTED)
    names = set(realrun.changed_functions(repo, "main", "HEAD")["names"])
    assert "pkg.mod.C.m.<locals>.deep" in names, names


def test_a_module_level_change_is_a_module_level_name(tmp_path):
    repo = _toy_repo(tmp_path, _BASE, _HEAD_MODULE)
    names = set(realrun.changed_functions(repo, "main", "HEAD")["names"])
    assert names == {"pkg.mod.<module>"}, names


def test_a_non_python_change_names_nothing_and_says_so(tmp_path):
    repo = _toy_repo(tmp_path, _BASE, _BASE)
    _git(repo, "checkout", "-q", "change")
    (repo / "notes.txt").write_text("x = 2\n")
    _git(repo, "add", "-A")
    _git(repo, "commit", "-qm", "only text")
    got = realrun.changed_functions(repo, "main", "HEAD")
    assert got["names"] == {} and "no .py file" in got["note"], got


def test_the_changed_names_reach_the_child_before_it_runs(tmp_path, monkeypatch):
    """The names are handed to the observer, not only used after the run.

    Without this the child could not report a hit at all: the parent's set of
    names and the child's copy must be the same set.
    """
    d, env = realrun._instrument(tmp_path, "branch", {"pkg.mod.helper": "pkg/mod.py"})
    assert json.loads((d / "names.json").read_text()) == ["pkg.mod.helper"]
    assert env["DELFIN_PROBE_NAMES"].endswith("names.json")
    assert Path(env["DELFIN_PROBE_OUT"]).name == "trace.json"
    assert "sys.setprofile" in (d / "sitecustomize.py").read_text()
    assert "threading.setprofile" in (d / "sitecustomize.py").read_text()


def test_a_stale_trace_is_not_read_as_a_reached_change(tmp_path):
    """An old trace file must not survive into the next run's answer."""
    d, env = realrun._instrument(tmp_path, "branch", {"pkg.mod.helper": "pkg/mod.py"})
    Path(env["DELFIN_PROBE_OUT"]).write_text(json.dumps(
        {"offered": ["bash"], "funcs": ["pkg.mod.helper"]}))
    realrun._instrument(tmp_path, "branch", {"pkg.mod.helper": "pkg/mod.py"})
    assert not Path(env["DELFIN_PROBE_OUT"]).exists()


def test_the_summary_reads_the_tools_a_run_called(tmp_path):
    """``called`` from one log, both sources: emitted events and the result line.

    Control: the counts and the stop reason are what they were before this
    change, so the new column did not displace the old accounting.
    """
    log = tmp_path / "one.jsonl"
    log.write_text("\n".join([
        json.dumps({"type": "tool_use", "name": "bash", "input": {}}),
        json.dumps({"type": "notice", "text": "not a call"}),
        json.dumps({"type": "result", "stop_reason": "done", "input_tokens": 1,
                    "output_tokens": 2, "tool_calls": [{"name": "read_file"},
                                                       {"name": "bash"}],
                    "text": "answer"}),
    ]) + "\n")
    got = realrun._summary(log)
    assert got["tools"] == 1 and got["stop"] == "done"
    assert got["called"] == ["bash", "read_file"], got["called"]


def test_the_column_aggregate_counts_only_the_branch_side():
    """The base has no change to reach, so its row must not decide the word."""
    assert realrun._reached_overall(
        [{"side": "base", "reached": "yes"}, {"side": "branch", "reached": "no"}]) == "no"
    assert realrun._reached_overall(
        [{"side": "branch", "reached": "no"}, {"side": "branch", "reached": "yes"}]) == "yes"
    # A probe that cannot tell outranks a clean miss: "we do not know" is not
    # the same claim as "it did not run".
    assert realrun._reached_overall(
        [{"side": "branch", "reached": "no"},
         {"side": "branch", "reached": "unknown"}]) == "unknown"


# -- the stage's note --------------------------------------------------------------

@pytest.fixture
def _waves(tmp_path, monkeypatch):
    """A run dir with one ended member whose real runs have been read."""
    d = tmp_path / "waves" / "note"
    d.mkdir(parents=True)
    cfg = dict(wave._DEFAULTS, waves_root=str(tmp_path / "waves"), min_available_gb=0,
               max_load_per_core=1e9)
    monkeypatch.setattr(wave, "settings", lambda: dict(cfg))
    (d / "state.json").write_text(json.dumps({"repo": str(tmp_path), "agents": {
        "a1": {"worktree": str(tmp_path / "wt"), "log": str(d / "a1.log"),
               "pid": _dead_pid(), "task": {"id": "a1", "brief": "b", "base": "main"}}}}))
    (d / "tasks.json").write_text(json.dumps([{"id": "a1", "brief": "b", "base": "main"}]))
    (d / "a1.log").write_text(json.dumps({"type": "result", "stop_reason": "done"}) + "\n")
    # What routes the member to the real runs at all: without it the stage asks
    # for a request instead of making one, and the note under test is never
    # the note the member receives.
    (tmp_path / "wt" / "scratch").mkdir(parents=True)
    (tmp_path / "wt" / "scratch" / "realrun.json").write_text(json.dumps(
        {"probes": [{"name": "p", "prompt": "show it"}]}))
    return d


def _note_for(_waves, monkeypatch, reached):
    """The note the stage hands the member, driven through ``advance``."""
    seen = []
    monkeypatch.setattr(wave, "resume", lambda run, tid, note="", **k: seen.append(note))
    monkeypatch.setattr(realrun, "run",
                        lambda run, tid: {"rows": [], "reached": reached, "summary": "x"})
    did = stages.advance("note", start_real_run=lambda run, tid, busy:
                         stages._real_run_then_resume(run, tid, busy))
    assert [x["did"] for x in did] == ["real run"], did
    assert seen, "the member was never resumed with a note"
    return seen[0]


def test_the_note_says_plainly_when_no_probe_reached_the_change(_waves, monkeypatch):
    note = _note_for(_waves, monkeypatch, "no")
    assert "NONE OF THEM REACHED YOUR CHANGE" in note
    assert "reached the change" in note
    # It must not invite a second request the stage will never run: the
    # ``.realrun`` marker makes the member final at its next ending.
    assert "Do not write scratch/realrun.json again" in note


def test_the_note_asks_for_the_column_to_be_read_when_a_probe_reached_it(_waves, monkeypatch):
    note = _note_for(_waves, monkeypatch, "yes")
    assert "NONE OF THEM REACHED YOUR CHANGE" not in note
    assert "reached the change" in note
    assert "Do not write scratch/realrun.json again" in note


def test_the_note_does_not_accuse_when_the_run_cannot_tell(_waves, monkeypatch):
    """``unknown`` gets the plain note: the change was not shown to be missed."""
    note = _note_for(_waves, monkeypatch, "unknown")
    assert "NONE OF THEM REACHED YOUR CHANGE" not in note


def _note_for_a_writer_that_records_no_column(_waves, monkeypatch):
    """What the stage says when the summary it was handed has no column.

    The writer that ran need not be the branch's: the orchestrator imports the
    installed DELFIN, so a member's real runs can be written by a copy that
    predates the column. Its result then carries no ``reached`` key, and the
    plain note would tell the member to read a column the table does not have.
    """
    seen = []
    monkeypatch.setattr(wave, "resume", lambda run, tid, note="", **k: seen.append(note))
    monkeypatch.setattr(realrun, "run", lambda run, tid: {"rows": [], "summary": "x"})
    did = stages.advance("note", start_real_run=lambda run, tid, busy:
                         stages._real_run_then_resume(run, tid, busy))
    assert [x["did"] for x in did] == ["real run"], did
    assert seen, "the member was never resumed with a note"
    return seen[0]


def test_the_note_says_the_summary_cannot_answer_when_it_has_no_column(_waves, monkeypatch):
    """No ``reached`` key: do not point at a column that is not in the table.

    Red on the base for the reason it names: there ``realrun.run`` returns no
    ``reached`` key either, so the stage hands the plain note, which says "Each
    row says whether that run reached the change: read that column" -- a column
    the delivered summary may not carry. The note is the only instruction the
    member gets, so a member following it looks for evidence that was never
    written and can read the base/branch difference as proof of its change.
    """
    note = _note_for_a_writer_that_records_no_column(_waves, monkeypatch)
    assert "THE SUMMARY DOES NOT SAY WHETHER ANY RUN REACHED YOUR CHANGE" in note
    # The false instruction is what must be gone, word for word.
    assert "read that column" not in note
    assert "Each row says whether that run reached the change" not in note
    # Not an accusation either: this is "cannot tell", not "none reached it".
    assert "NONE OF THEM REACHED YOUR CHANGE" not in note
    # Still no second request the stage would never run.
    assert "Do not write scratch/realrun.json again" in note
