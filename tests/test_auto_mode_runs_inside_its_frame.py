"""Auto mode: the harness starts the next iteration, inside a frame it enforces.

Input: an AutoController for one session (workspace, session id, an
injectable clock, check runner and notifier), a stand-in engine (its
permission object, token usage, messages, a stop recorder), and the turns
it is told about. Output: the iteration prompts it hands the surface, the
verdicts its veto gives, the commits, reverts and log lines it makes, and
the stop report.

No model runs here: the "turn" is whatever the test writes into the
workspace between ``next_prompt`` and ``finish_iteration``.
"""

from __future__ import annotations

import json
import subprocess
from pathlib import Path
from types import SimpleNamespace

import pytest

from delfin.agent import auto_mode as am
from delfin.agent import loops
from delfin.agent.harness_rows import HarnessText


# ---------------------------------------------------------------------------
# fixtures
# ---------------------------------------------------------------------------

class Clock:
    def __init__(self):
        self.t = 1_000_000.0

    def __call__(self):
        return self.t

    def advance(self, s):
        self.t += s


class Engine:
    def __init__(self, mode="acceptEdits"):
        self.kit_permissions = SimpleNamespace(mode=mode, auto_veto=None)
        self.token_usage = {"input": 0, "output": 0}
        self.messages = []
        self.session_id = "s1"
        self.stops = 0
        self.pending_turn_continuation = ""

    def request_stop(self):
        self.stops += 1

    def spend(self, n):
        self.token_usage["input"] += n


def _git(ws, *args):
    return subprocess.run(["git", "-C", str(ws), *args], capture_output=True,
                          text=True, check=True).stdout.strip()


@pytest.fixture()
def home(tmp_path, monkeypatch):
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: h))
    return h


@pytest.fixture()
def repo(tmp_path, home):
    ws = tmp_path / "ws"
    ws.mkdir()
    _git(ws, "init", "-q", "-b", "main")
    _git(ws, "config", "user.name", "t")
    _git(ws, "config", "user.email", "t@example.org")
    _git(ws, "config", "commit.gpgsign", "false")
    (ws / "results").mkdir()
    (ws / "results" / "table.csv").write_text("job,energy\n", encoding="utf-8")
    (ws / "scripts").mkdir()
    (ws / "scripts" / "check_table.py").write_text("print('missing=3')\n",
                                                   encoding="utf-8")
    (ws / "tests").mkdir()
    (ws / "tests" / "test_table.py").write_text(
        "def test_rows():\n    assert rows() == 3\n", encoding="utf-8")
    _git(ws, "add", "-A")
    _git(ws, "commit", "-q", "-m", "init")
    return ws


class Checks:
    """The check command's results, one per call, last one repeated."""

    def __init__(self, *results):
        self.results = list(results) or [(0, "")]
        self.calls = []

    def __call__(self, engine, command):
        self.calls.append(command)
        r = self.results.pop(0) if len(self.results) > 1 else self.results[0]
        return r


def _ctl(ws, clock=None, checks=None, **settings):
    notes = []
    ctl = am.AutoController(
        ws, session_id="s1", clock=clock or Clock(),
        run_check=checks or Checks(), notify=notes.append,
        settings={"idle_seconds": 60, "loop_gap_seconds": 5, **settings})
    ctl._notes = notes
    return ctl


LOOP = """---
goal: Every finished job has its energy in results/table.csv
trigger: manual
check: python scripts/check_table.py
metric: 'missing=(\\d+)'
better: lower
budget: {{iterations: {its}, wall: 2h, tokens: 50000}}
stop: {{plateau: {plateau}, repeat_failures: 2}}
touch: [results/]
---
1. Add the next missing row to results/table.csv.
"""


def _loop(ws, its=5, plateau=3, name="table"):
    d = ws / ".delfin" / "loops"
    d.mkdir(parents=True, exist_ok=True)
    (d / f"{name}.md").write_text(LOOP.format(its=its, plateau=plateau),
                                  encoding="utf-8")
    spec = loops.load(name, ws)
    assert spec.errors == []
    loops.confirm(spec)
    return loops.load(name, ws)


def _iterate(ctl, eng, clock, write=None, answer="LOOP: progress row added"):
    clock.advance(10)
    prompt = ctl.next_prompt(eng)
    assert isinstance(prompt, HarnessText), ctl.status_line()
    if write:
        for rel, text in write.items():
            p = ctl.workspace / rel
            p.parent.mkdir(parents=True, exist_ok=True)
            p.write_text(text, encoding="utf-8")
    return ctl.finish_iteration(eng, answer)


# ---------------------------------------------------------------------------
# enabling: the posture is taken as it is, never widened
# ---------------------------------------------------------------------------

def test_auto_does_not_start_under_bypass(repo):
    ctl = _ctl(repo)
    msg = ctl.enable(Engine(mode="bypassPermissions"))
    assert not ctl.enabled and "bypass" in msg.lower()


def test_auto_needs_the_native_gate(repo):
    ctl = _ctl(repo)
    eng = Engine()
    eng.kit_permissions = None
    assert "gate" in ctl.enable(eng).lower() and not ctl.enabled


def test_auto_leaves_the_default_branch_for_its_own(repo):
    ctl = _ctl(repo)
    ctl.enable(Engine())
    assert ctl.enabled
    branch = _git(repo, "branch", "--show-current")
    assert branch.startswith("auto/") and branch != "main"


def test_a_dirty_default_branch_refuses(repo):
    (repo / "results" / "table.csv").write_text("x\n", encoding="utf-8")
    ctl = _ctl(repo)
    msg = ctl.enable(Engine())
    assert not ctl.enabled and "main" in msg
    assert _git(repo, "branch", "--show-current") == "main"


def test_a_posture_widened_after_enabling_stops_auto(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine(mode="default")
    ctl.enable(eng)
    eng.kit_permissions.mode = "acceptEdits"
    clock.advance(120)
    assert ctl.next_prompt(eng, work=[am.WorkItem(1, "task", "x", True)]) == ""
    assert not ctl.enabled and "posture" in ctl.stop_report


# ---------------------------------------------------------------------------
# idle and triggers
# ---------------------------------------------------------------------------

def test_nothing_starts_before_the_idle_time(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    ctl.enable(eng)
    work = [am.WorkItem(1, "task", "task #1 in progress: write the parser",
                        True)]
    ctl.note_turn_end()
    clock.advance(30)
    assert ctl.next_prompt(eng, work=work) == ""
    clock.advance(31)
    p = ctl.next_prompt(eng, work=work)
    assert isinstance(p, HarnessText) and "write the parser" in p
    assert "not the user" in p


def test_typed_text_and_a_pause_hold_it(repo, monkeypatch):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    ctl.enable(eng)
    work = [am.WorkItem(1, "task", "x", True)]
    clock.advance(120)
    assert ctl.next_prompt(eng, typed="half a sent", work=work) == ""
    from delfin.agent import session_pause
    monkeypatch.setattr(session_pause, "wake_blocked", lambda key: True)
    assert ctl.next_prompt(eng, work=work, pause_key="k") == ""


def test_off_means_off(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    clock.advance(1000)
    assert ctl.next_prompt(eng, work=[am.WorkItem(1, "task", "x", True)]) == ""


def test_an_every_trigger_waits_for_its_interval(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    d = repo / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "poll.md").write_text(
        "---\ngoal: poll\ntrigger: every 30m\n---\n1. Look.\n", encoding="utf-8")
    spec = loops.load("poll", repo)
    loops.confirm(spec)
    assert "started on branch" in ctl.start_loop(loops.load("poll", repo), eng)
    clock.advance(10)
    assert isinstance(ctl.next_prompt(eng), HarnessText)
    ctl.finish_iteration(eng, "LOOP: progress looked")
    clock.advance(600)
    assert ctl.next_prompt(eng) == ""
    clock.advance(1200)
    assert isinstance(ctl.next_prompt(eng), HarnessText)


def test_a_job_trigger_fires_on_a_wake_only(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    d = repo / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "jobs.md").write_text(
        "---\ngoal: evaluate jobs\ntrigger: job\n---\n1. Evaluate it.\n",
        encoding="utf-8")
    loops.confirm(loops.load("jobs", repo))
    ctl.start_loop(loops.load("jobs", repo), eng)
    clock.advance(3600)
    assert ctl.next_prompt(eng) == ""
    wake = "[watch — a system event, not the user] A job ... finished:\n- slurm 4711 [COMPLETED]"
    p = ctl.next_prompt(eng, wake=wake)
    assert isinstance(p, HarnessText) and "4711" in p and "Evaluate it" in p


def test_an_unconfirmed_loop_does_not_start(repo):
    d = repo / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "x.md").write_text("---\ngoal: g\n---\n1. a\n", encoding="utf-8")
    ctl = _ctl(repo)
    msg = ctl.start_loop(loops.load("x", repo), Engine())
    assert "confirm" in msg and ctl.loop is None


# ---------------------------------------------------------------------------
# an iteration: commit, check, keep or revert, log
# ---------------------------------------------------------------------------

def test_a_kept_iteration_is_committed_on_the_loop_branch_and_logged(repo):
    clock = Clock()
    checks = Checks((1, "missing=2"))
    ctl = _ctl(repo, clock, checks)
    eng = Engine()
    spec = _loop(repo)
    ctl.start_loop(spec, eng)
    assert _git(repo, "branch", "--show-current") == "loop/table"
    digest = _iterate(ctl, eng, clock,
                      write={"results/table.csv": "job,energy\nj1,-1.0\n"})
    assert "kept" in digest and "2" in digest
    log = _git(repo, "log", "--format=%s", "-1")
    assert "table" in log and "iteration 1" in log
    assert _git(repo, "show", "--stat", "--format=", "HEAD").count("|") == 1
    rows = [json.loads(l) for l in ctl.log_path.read_text().splitlines()]
    assert rows[-1]["kept"] is True and rows[-1]["metric"] == 2.0
    assert rows[-1]["changed"] == ["results/table.csv"]
    # The baseline before the first iteration, then the iteration's own.
    assert checks.calls == ["python scripts/check_table.py"] * 2


def test_a_worse_measurement_is_reverted(repo):
    clock = Clock()
    checks = Checks((1, "missing=3"), (1, "missing=2"), (1, "missing=4"))
    ctl = _ctl(repo, clock, checks)
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    _iterate(ctl, eng, clock, write={"results/table.csv": "a\n"})
    kept_head = _git(repo, "rev-parse", "HEAD")
    digest = _iterate(ctl, eng, clock, write={"results/table.csv": "b\n"})
    assert "revert" in digest.lower()
    assert (repo / "results" / "table.csv").read_text() == "a\n"
    assert _git(repo, "rev-parse", "HEAD~2") == kept_head
    rows = [json.loads(l) for l in ctl.log_path.read_text().splitlines()]
    assert rows[-1]["kept"] is False and "worse" in rows[-1]["why"]


def test_an_edit_of_the_measure_is_reverted(repo):
    clock = Clock()
    checks = Checks((1, "missing=3"), (0, "missing=0"))
    ctl = _ctl(repo, clock, checks)
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    digest = _iterate(ctl, eng, clock, write={
        "scripts/check_table.py": "print('missing=0')\n",
        "results/table.csv": "x\n"})
    assert "revert" in digest.lower()
    assert "missing=3" in (repo / "scripts" / "check_table.py").read_text()
    rows = [json.loads(l) for l in ctl.log_path.read_text().splitlines()]
    assert rows[-1]["kept"] is False and "measure" in rows[-1]["why"]
    # And it is not "goal met": a pass bought by editing the check is none.
    assert ctl.loop is not None or "goal" not in ctl.stop_report


def test_a_test_edit_outside_the_declared_measure_is_kept_and_flagged(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=2")))
    eng = Engine()
    d = repo / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "t.md").write_text(LOOP.format(its=5, plateau=3).replace(
        "touch: [results/]", "touch: []"), encoding="utf-8")
    loops.confirm(loops.load("t", repo))
    ctl.start_loop(loops.load("t", repo), eng)
    _iterate(ctl, eng, clock, write={
        "tests/test_table.py": "def test_rows():\n    pass\n"})
    rows = [json.loads(l) for l in ctl.log_path.read_text().splitlines()]
    assert rows[-1]["kept"] is True and rows[-1]["flagged"]
    assert any("tests/test_table.py" in n["text"] for n in ctl.needs_you)


def test_a_goal_met_at_the_start_does_not_start(repo):
    ctl = _ctl(repo, Clock(), Checks((0, "missing=0")))
    msg = ctl.start_loop(_loop(repo), Engine())
    assert ctl.loop is None and "already" in msg


def test_goal_met_stops_with_a_report(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3"), (0, "missing=0")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    _iterate(ctl, eng, clock, write={"results/table.csv": "done\n"})
    assert ctl.loop is None and "goal met" in ctl.stop_report
    assert any("goal met" in n for n in ctl._notes)


# ---------------------------------------------------------------------------
# stops: budget, plateau, repeated failure, the agent's own word
# ---------------------------------------------------------------------------

def test_the_iteration_budget_stops_it(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks(*[(1, f"missing={9 - i}") for i in range(9)]))
    eng = Engine()
    ctl.start_loop(_loop(repo, its=2), eng)
    for i in range(2):
        _iterate(ctl, eng, clock, write={"results/table.csv": f"{i}\n"})
    clock.advance(10)
    assert ctl.next_prompt(eng) == ""
    assert ctl.loop is None and "2 of 2 iterations" in ctl.stop_report


def test_the_token_budget_stops_it(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3"), (1, "missing=2"),
                                   (1, "missing=1")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    clock.advance(10)
    ctl.next_prompt(eng)
    eng.spend(60000)
    ctl.finish_iteration(eng, "LOOP: progress")
    clock.advance(10)
    assert ctl.next_prompt(eng) == ""
    assert "token" in ctl.stop_report


def test_the_wall_budget_stops_it(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=2"), (1, "missing=1")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    _iterate(ctl, eng, clock, write={"results/table.csv": "1\n"})
    clock.advance(3 * 3600)
    assert ctl.next_prompt(eng) == "" and "time" in ctl.stop_report


def test_a_plateau_stops_it(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3")))
    eng = Engine()
    ctl.start_loop(_loop(repo, plateau=2), eng)
    # Changes that leave the measurement where it was: no progress.
    _iterate(ctl, eng, clock, write={"results/table.csv": "a\n"})
    _iterate(ctl, eng, clock, write={"results/table.csv": "b\n"})
    assert ctl.loop is None and "plateau" in ctl.stop_report


def test_the_same_failure_twice_stops_it(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((2, "Traceback\nKeyError: 'energy'")))
    eng = Engine()
    ctl.start_loop(_loop(repo, plateau=9), eng)
    _iterate(ctl, eng, clock)
    assert ctl.loop is not None
    _iterate(ctl, eng, clock)
    assert ctl.loop is None and "same failure" in ctl.stop_report


def test_needs_human_pauses_and_queues_the_question(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    _iterate(ctl, eng, clock,
             answer="Two folders.\nLOOP: needs-human which folder is the run?")
    assert ctl.paused and any("which folder" in n["text"]
                              for n in ctl.needs_you)
    clock.advance(100)
    assert ctl.next_prompt(eng) == ""


# ---------------------------------------------------------------------------
# the user comes first
# ---------------------------------------------------------------------------

def test_user_input_preempts_after_the_current_tool_call(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=2")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    clock.advance(10)
    ctl.next_prompt(eng)
    assert eng.kit_permissions.auto_veto is not None
    # A tool call before the user typed: allowed.
    assert ctl.veto("read_file", {"path": "x"}) is None
    assert ctl.note_user_input(eng) is True
    # The next call is refused and the turn asked to stop -- between calls,
    # so nothing that already runs is killed.
    refusal = ctl.veto("read_file", {"path": "x"})
    assert refusal and "user" in refusal
    assert eng.stops == 1
    ctl.finish_iteration(eng, "stopped")
    assert eng.kit_permissions.auto_veto is None
    assert ctl.paused
    clock.advance(1000)
    assert ctl.next_prompt(eng) == ""
    ctl.resume()
    clock.advance(10)
    assert isinstance(ctl.next_prompt(eng), HarnessText)


def test_a_loop_that_says_resume_resumes_after_the_users_turn(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=2"), (1, "missing=1")))
    eng = Engine()
    d = repo / ".delfin" / "loops"
    d.mkdir(parents=True)
    (d / "r.md").write_text(LOOP.format(its=5, plateau=3).replace(
        "trigger: manual", "trigger: manual\nafter_user: resume"),
        encoding="utf-8")
    loops.confirm(loops.load("r", repo))
    ctl.start_loop(loops.load("r", repo), eng)
    ctl.note_user_input(eng)          # the user's own turn
    ctl.note_turn_end()
    clock.advance(10)
    assert isinstance(ctl.next_prompt(eng), HarnessText)


# ---------------------------------------------------------------------------
# what only the user may do is queued, not done
# ---------------------------------------------------------------------------

def test_a_push_in_an_iteration_is_queued_not_run(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    clock.advance(10)
    ctl.next_prompt(eng)
    out = ctl.veto("bash", {"command": "git push -u origin loop/table"})
    assert out and "queued" in out
    item = ctl.needs_you[-1]
    assert item["kind"] == "outward" and "git push" in item["text"]
    assert ctl._notes          # told out-of-band (attention inbox)
    # A measure write is refused, not queued.
    out = ctl.veto("write_file", {"path": "scripts/check_table.py"})
    assert out and "queued" not in out
    assert ctl.needs_you[-1]["kind"] == "outward"


def test_the_veto_is_disarmed_outside_an_iteration(repo):
    ctl = _ctl(repo)
    assert ctl.veto("bash", {"command": "git push"}) is None


# ---------------------------------------------------------------------------
# self-found work
# ---------------------------------------------------------------------------

def test_self_found_work_is_ranked_and_scoped(repo):
    from delfin.agent import agent_tasks
    store = agent_tasks.get_store(repo)
    store.create("pending thing", session_id="s1")
    t = store.create("started thing", session_id="s1")
    store.update(t["id"], status="in_progress")
    store.create("someone else's", session_id="other")
    (repo / "TODO.md").write_text("- [ ] tidy the docs\n- [x] done\n",
                                  encoding="utf-8")
    eng = Engine()
    eng.messages = [{"role": "user", "content": "gate tests/ -> 2 failed, 5 passed"}]
    ctl = _ctl(repo)
    work = ctl.find_work(eng)
    texts = [w.text for w in work]
    inscope = [w for w in work if w.in_scope]
    assert "started thing" in inscope[0].text
    assert any("pending thing" in w.text for w in inscope)
    assert any("2 failed" in w.text for w in inscope)
    out = [w for w in work if not w.in_scope]
    assert any("tidy the docs" in w.text for w in out)
    assert any("someone else's" in w.text for w in out)
    assert not any("done" == t for t in texts)


def test_out_of_scope_work_is_asked_about_not_started(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    ctl.enable(eng)
    clock.advance(120)
    p = ctl.next_prompt(eng, work=[am.WorkItem(5, "register", "tidy the docs",
                                               False)])
    assert p == ""
    assert any("tidy the docs" in n["text"] for n in ctl.needs_you)
    notes = ctl.pop_notices()
    assert any("nothing in scope" in n for n in notes)
    # Said once, not every poll.
    clock.advance(120)
    ctl.next_prompt(eng, work=[am.WorkItem(5, "register", "tidy the docs",
                                           False)])
    assert ctl.pop_notices() == []
    assert len([n for n in ctl.needs_you if "tidy" in n["text"]]) == 1


def test_the_ranked_list_is_shown_before_the_iteration(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    ctl.enable(eng)
    clock.advance(120)
    work = [am.WorkItem(2, "task", "task #2 in progress: parser", True),
            am.WorkItem(4, "task", "task #3 pending: docs", True)]
    p = ctl.next_prompt(eng, work=work)
    assert "parser" in p and "docs" not in p
    shown = "\n".join(ctl.pop_notices())
    assert "parser" in shown and "docs" in shown


def test_status_line_and_state(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=2")))
    eng = Engine()
    assert ctl.status_line() == ""
    ctl.start_loop(_loop(repo), eng)
    _iterate(ctl, eng, clock, write={"results/table.csv": "1\n"})
    line = ctl.status_line()
    assert "table" in line and "1/5" in line
    st = ctl.state()
    assert st["loop"] == "table" and st["iteration"] == 1
    assert st["max_iterations"] == 5 and st["last"].startswith("kept")
    assert "tokens_left" in st and "needs_you" in st


def test_an_iteration_whose_turn_never_started_is_taken_back(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    clock.advance(10)
    assert isinstance(ctl.next_prompt(eng), HarnessText)
    ctl.cancel_iteration("the surface held the turn")
    assert ctl.iteration is None and eng.kit_permissions.auto_veto is None
    assert ctl.state()["iteration"] == 0
    clock.advance(10)
    assert isinstance(ctl.next_prompt(eng), HarnessText)
    assert ctl.state()["iteration"] == 1


def test_the_harness_never_runs_a_publishing_check(repo):
    # A spec that slipped past validation (built in code) is still held:
    # the check is classified before the harness runs it.
    checks = Checks((1, "missing=3"))
    ctl = _ctl(repo, Clock(), checks)
    code, out = ctl._check(Engine(), "git push origin HEAD")
    assert code != 0 and "not run" in out
    assert checks.calls == []


def test_a_claim_the_check_contradicts_is_queued_for_review(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, Checks((1, "missing=3"), (1, "missing=2")))
    eng = Engine()
    ctl.start_loop(_loop(repo), eng)
    # The metric improved (kept), but the answer claims more than the
    # failing check shows.
    _iterate(ctl, eng, clock, write={"results/table.csv": "a\n"},
             answer="All tests pass now.\nLOOP: progress")
    rows = [json.loads(l) for l in ctl.log_path.read_text().splitlines()]
    assert rows[-1]["claims"] and rows[-1]["kept"] is True
    assert any("claims" in n["text"] for n in ctl.needs_you)


def test_a_queued_call_without_a_command_says_what_it_carried(repo):
    clock = Clock()
    ctl = _ctl(repo, clock)
    eng = Engine()
    ctl.enable(eng)
    clock.advance(120)
    ctl.next_prompt(eng, work=[am.WorkItem(1, "task", "x", True)])
    ctl.veto("session_message", {"to": "partner", "text": "the table is done"})
    assert "the table is done" in ctl.needs_you[-1]["text"]


def test_the_same_item_without_a_change_is_named_as_such(repo):
    clock = Clock()
    ctl = _ctl(repo, clock, repeat_failures=2)
    eng = Engine()
    ctl.enable(eng)
    work = [am.WorkItem(1, "task", "task 1 in progress: parser", True)]
    for _ in range(2):
        clock.advance(120)
        assert ctl.next_prompt(eng, work=work)
        ctl.finish_iteration(eng, "LOOP: progress none")
    assert not ctl.enabled
    assert "twice in a row without a change" in ctl.stop_report
