"""Every benchmark round is kept, and the trend can be read (package m1).

`bench_nightly.run` produced a RunRecord and nothing stored it, so no one
could say whether DELFIN got better last week. These tests drive the
product's own entry points — `bench_nightly.main` for a round,
`delfin.agent.cli.main(["bench", ...])` for the read side — never the
helpers, so a verb that exists but is not wired into the parser fails
here.

Red on origin/main for the reason each name gives: there is no
`delfin.agent.bench_history` module and no `bench history` / `bench trend`
/ `bench schedule` verb there, so every case that names one fails with
ImportError or argparse's "invalid choice". The two `test_CONTROL_` cases
pass on main and after — they are the never-worse cases, the seed gate and
the budget gate, which must survive the store being wired into run().

No model call anywhere: the runner and grader are fakes injected through
the glue, exactly as tests/test_hb1_cli.py does it.
"""

import json

import pytest

from delfin.agent import bench_nightly as bn
from delfin.agent import bench_holdout_glue as glue


# ---------------------------------------------------------------------------
# The faked pipeline (no engine, no key, no model)
# ---------------------------------------------------------------------------

class _FakeAccounting:
    def __init__(self, estimate=0.1):
        self._estimate = estimate
        self.charged = 0.0

    def estimate(self, task):
        return self._estimate

    def charge(self, core_hours):
        self.charged += core_hours


class _Sampler:
    """Draws `n` tasks of one capability, ids unique per draw."""

    def __init__(self, capability="redox_consistency", n=2):
        self.capability = capability
        self.n = n

    def sample(self, seed, n):
        from delfin.agent.benchmark import Task
        return [Task(id=f"task{i}", task_class=self.capability, mode="solo",
                     prompt=f"prompt {self.capability} {i}")
                for i in range(min(n, self.n))]


class _Grader:
    """Passes everything at a fixed quality, or fails the tasks named."""

    def __init__(self, fail=(), quality=80):
        self.fail = set(fail)
        self.quality = quality

    def grade(self, task, result):
        return (task.id not in self.fail, self.quality)


def _runner(*, unmeasured=(), tokens=(100, 20)):
    """A runner that returns a BenchmarkResult; names which tasks come
    back unmeasured (the run never happened)."""

    def runner(task, *, model="", **kw):
        from delfin.agent.benchmark import BenchmarkResult
        res = BenchmarkResult(task_id=task.id, task_class=task.task_class,
                              model=model)
        res.input_tokens, res.output_tokens = tokens
        if task.id in set(unmeasured):
            res.unmeasured = True
            res.error = "engine init failed: no API key"
        else:
            res.success = True
        return res

    return runner


@pytest.fixture
def fake_glue(monkeypatch):
    """Install a config builder that fakes sampler/runner/grader and keeps
    everything else of run() real: the budget gate, the seed gate, the
    record, and the store."""
    state = {"cfg": None}

    def fake_build(*, model, instances, mode, holdout_seed_path="",
                   provider="", public_seed=7, **kw):
        cfg = {
            "model": model, "provider": provider,
            "sampler": state.get("sampler") or _Sampler(),
            "grader": state.get("grader") or _Grader(),
            "accounting": _FakeAccounting(),
            "runner": state.get("runner") or _runner(),
            "instances": instances, "mode": mode,
            "holdout_seed_path": holdout_seed_path,
            "public_seed": public_seed,
        }
        cfg.update(state.get("extra", {}))
        state["cfg"] = cfg
        return cfg

    monkeypatch.setattr(glue, "build_default_config", fake_build)
    return state


@pytest.fixture
def history_dir(tmp_path, monkeypatch):
    """Point the store's DEFAULT root at this test's directory.

    Monkeypatching the resolver is the product's own seam: the name is
    registered in state_paths.USER_STATE_RESOLVERS precisely so the root
    can be moved wholesale. Without it these tests would append fixture
    rounds to the operator's real history — and an append-only store is
    one nothing later cleans.

    The import degrades to a no-op where the module does not exist, so the
    two CONTROL cases below can use this fixture for isolation and still
    pass on the baseline (which has no store to redirect).
    """
    try:
        from delfin.agent import bench_history as bh
    except ImportError:
        return tmp_path / "kept-rounds"

    root = tmp_path / "kept-rounds"
    monkeypatch.setattr(bh, "default_root", lambda: root)
    return root


def _round(fake_glue, tmp_path, *, capability="redox_consistency",
           practice=True, model="M", fail=(), unmeasured=(), seed=None,
           tokens=(100, 20), extra=None):
    """One round through `python -m delfin.agent.bench_nightly`'s main()."""
    fake_glue["sampler"] = _Sampler(capability=capability)
    fake_glue["grader"] = _Grader(fail=fail)
    fake_glue["runner"] = _runner(unmeasured=unmeasured, tokens=tokens)
    fake_glue["extra"] = dict(extra or {})
    argv = ["--model", model, "--instances", "2"]
    if practice:
        argv.append("--practice")
    else:
        seed = seed or (tmp_path / "seed.txt")
        seed.write_text("424242\n")
        argv += ["--seed-path", str(seed)]
    return bn.main(argv)


def _lines(path):
    if not path.is_file():
        return []
    return [json.loads(ln) for ln in
            path.read_text(encoding="utf-8").splitlines() if ln.strip()]


# ---------------------------------------------------------------------------
# Item 1 — the round is kept
# ---------------------------------------------------------------------------

def test_a_round_that_ran_is_kept_on_disk(fake_glue, history_dir, tmp_path,
                                          capsys):
    """RED on main: a round leaves no file at all — bench_nightly.main()
    printed its record and that was the end of it."""
    from delfin.agent import bench_history as bh

    assert _round(fake_glue, tmp_path, practice=True) == 0
    rounds = _lines(history_dir / bh.ROUNDS_FILE)
    instances = _lines(history_dir / bh.INSTANCES_FILE)
    assert len(rounds) == 1
    assert len(instances) == 2
    # Every field the brief names, on the round line.
    round_line = rounds[0]
    for field_name in ("round_id", "commit", "model", "provider", "seed_mode",
                       "task_set_version"):
        assert field_name in round_line, field_name
    # ... and per-instance: task, capability, pass/fail/score, unmeasured,
    # duration, tokens. The capability tag and the score are the record's
    # own fields, task_class and quality (bench_report scores per tool
    # family from exactly task_class) — see bench_history.inst_capability.
    inst = instances[0]
    for field_name in ("task_id", "task_class", "success", "quality",
                       "unmeasured", "duration_s", "input_tokens",
                       "output_tokens", "round_id"):
        assert field_name in inst, field_name
    assert inst["task_class"] == "redox_consistency"
    assert inst["success"] is True and inst["quality"] == 80
    assert (inst["input_tokens"], inst["output_tokens"]) == (100, 20)
    assert round_line["n_unmeasured"] == 0
    assert round_line["input_tokens"] == 200


def test_the_commit_and_round_id_written_are_the_ones_that_ran(fake_glue,
                                                              history_dir,
                                                              tmp_path):
    """The record says WHICH code was measured: without the commit, a flat
    trend cannot be attributed to anything."""
    from delfin.agent import bench_history as bh

    assert _round(fake_glue, tmp_path) == 0
    line = _lines(history_dir / bh.ROUNDS_FILE)[0]
    assert line["commit"], "no commit recorded"
    assert line["round_id"], "no round id recorded"


def test_two_rounds_are_both_kept_and_read_back(fake_glue, history_dir,
                                                tmp_path, capsys):
    """Two rounds run, then `delfin-agent bench history` shows both, each
    with its commit — driven through the CLI, not the helper."""
    from delfin.agent import cli

    assert _round(fake_glue, tmp_path) == 0
    assert _round(fake_glue, tmp_path) == 0
    assert capsys.readouterr().out          # drop the two round dumps
    rc = cli.main(["bench", "history"])
    out = capsys.readouterr().out
    assert rc == 0, out
    rounds = _lines(history_dir / "rounds.jsonl")
    assert len(rounds) == 2
    for line in rounds:
        assert line["round_id"] in out, "a kept round is missing from the list"
    assert "100%" in out


def test_history_names_the_per_task_pass_fail_next_to_the_pass_rate(
        fake_glue, history_dir, tmp_path, capsys):
    """Per capability AND per task: a rate alone hides which task moved."""
    from delfin.agent import cli

    assert _round(fake_glue, tmp_path, fail=("task1",)) == 0
    capsys.readouterr()
    assert cli.main(["bench", "history"]) == 0
    out = capsys.readouterr().out
    assert "task0" in out and "task1" in out
    assert "pass" in out and "fail" in out
    assert "50%" in out                      # 1 of 2, over the graded ones


# ---------------------------------------------------------------------------
# Item 2 — the trend, and only comparable rounds
# ---------------------------------------------------------------------------

def test_bench_trend_shows_the_pass_rate_per_round_for_a_capability(
        fake_glue, history_dir, tmp_path, capsys):
    from delfin.agent import cli

    assert _round(fake_glue, tmp_path) == 0
    assert _round(fake_glue, tmp_path, fail=("task0",)) == 0
    capsys.readouterr()
    rc = cli.main(["bench", "trend", "--capability", "redox_consistency"])
    out = capsys.readouterr().out
    assert rc == 0, out
    assert "100%" in out and "50%" in out


def test_bench_trend_without_a_capability_is_refused(fake_glue, history_dir,
                                                    capsys):
    """Guessing a capability would print a rate over a mixed population.

    Two ways to get there, both refused with exit 2: no flag at all
    (argparse's own required-argument exit, which is how the product
    answers a bare `bench trend`), and an empty value, which argparse
    accepts and only the handler can turn away.
    """
    from delfin.agent import cli

    with pytest.raises(SystemExit) as exc:
        cli.main(["bench", "trend"])
    assert exc.value.code == 2
    assert "capability" in capsys.readouterr().err.lower()

    assert cli.main(["bench", "trend", "--capability", ""]) == 2
    assert "capability" in capsys.readouterr().err.lower()


def test_a_round_on_another_seed_mode_is_excluded_and_says_why(
        fake_glue, history_dir, tmp_path, capsys):
    """A practice round and a holdout round are not the same measurement —
    the holdout is the one that must not be tuned on. Excluding it is
    only honest if the reader is told."""
    from delfin.agent import cli

    assert _round(fake_glue, tmp_path, practice=True) == 0
    assert _round(fake_glue, tmp_path, practice=False) == 0
    capsys.readouterr()
    assert cli.main(["bench", "trend", "--capability",
                     "redox_consistency"]) == 0
    out = capsys.readouterr().out
    assert "excluded" in out
    assert "seed_mode" in out, out


def test_a_round_of_another_model_is_excluded_and_says_why(
        fake_glue, history_dir, tmp_path, capsys):
    from delfin.agent import cli

    assert _round(fake_glue, tmp_path, model="model-A") == 0
    assert _round(fake_glue, tmp_path, model="model-B") == 0
    capsys.readouterr()
    assert cli.main(["bench", "trend", "--capability",
                     "redox_consistency"]) == 0
    out = capsys.readouterr().out
    assert "excluded" in out and "model" in out


def test_a_round_of_another_task_set_is_excluded_and_says_why(
        fake_glue, history_dir, tmp_path, capsys):
    """Same capability label, different prompt pack: the version is hashed
    from the drawn tasks, so a reworded pack cannot masquerade as the
    same measurement."""
    from delfin.agent import bench_history as bh, cli

    assert _round(fake_glue, tmp_path) == 0

    class _Reworded(_Sampler):
        def sample(self, seed, n):
            from delfin.agent.benchmark import Task
            # Same ids as _Sampler, different prompts: the pack moved.
            return [Task(id=f"task{i}", task_class="redox_consistency",
                         mode="solo", prompt=f"REWORDED {i}")
                    for i in range(2)]

    assert bn.run({"model": "M", "sampler": _Reworded(),
                   "grader": _Grader(), "accounting": _FakeAccounting(),
                   "runner": _runner(), "instances": 2,
                   "mode": "practice", "public_seed": 7}) is not None
    capsys.readouterr()
    assert cli.main(["bench", "trend", "--capability",
                     "redox_consistency"]) == 0
    out = capsys.readouterr().out
    assert "excluded" in out and "task_set_version" in out


def test_the_task_set_version_is_stable_and_moves_with_the_pack():
    """The axis itself: order-independent, and it changes when a prompt
    changes. A version that never moved would exclude nothing."""
    from delfin.agent.benchmark import Task

    a = [Task(id="t1", task_class="c", mode="solo", prompt="one"),
         Task(id="t2", task_class="c", mode="solo", prompt="two")]
    b = list(reversed(a))
    c = [Task(id="t1", task_class="c", mode="solo", prompt="one"),
         Task(id="t2", task_class="c", mode="solo", prompt="TWO")]
    assert bn.task_set_version(a) == bn.task_set_version(b)
    assert bn.task_set_version(a) != bn.task_set_version(c)


def test_two_rounds_are_comparable_when_only_their_results_differ(
        fake_glue, history_dir, tmp_path, capsys):
    """The other side of the gate: the same pack, mode and model must NOT
    be excluded, or the read side would drop every real trend."""
    from delfin.agent import cli

    assert _round(fake_glue, tmp_path) == 0
    assert _round(fake_glue, tmp_path, fail=("task1",)) == 0
    capsys.readouterr()
    assert cli.main(["bench", "trend", "--capability",
                     "redox_consistency"]) == 0
    out = capsys.readouterr().out
    assert "excluded" not in out


def test_the_change_between_two_named_rounds_is_named_per_task(
        fake_glue, history_dir, tmp_path, capsys):
    from delfin.agent import bench_history as bh, cli

    r1 = fake_glue
    assert _round(r1, tmp_path, fail=("task1",)) == 0
    first = _lines(history_dir / bh.ROUNDS_FILE)[0]["round_id"]
    assert _round(r1, tmp_path) == 0
    second = _lines(history_dir / bh.ROUNDS_FILE)[1]["round_id"]
    capsys.readouterr()
    assert cli.main(["bench", "diff", first, second]) == 0
    out = capsys.readouterr().out
    assert "fixed" in out, out           # task1 failed, then passed
    assert "50%" in out and "100%" in out


def test_a_diff_across_two_incomparable_rounds_is_refused(
        fake_glue, history_dir, tmp_path, capsys):
    """A delta between two models prints a number whose meaning no one can
    state. It must be refused, not shown."""
    from delfin.agent import bench_history as bh, cli

    assert _round(fake_glue, tmp_path, model="model-A") == 0
    a = _lines(history_dir / bh.ROUNDS_FILE)[0]["round_id"]
    assert _round(fake_glue, tmp_path, model="model-B") == 0
    b = _lines(history_dir / bh.ROUNDS_FILE)[1]["round_id"]
    capsys.readouterr()
    assert cli.main(["bench", "diff", a, b]) != 0
    assert "not comparable" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Item 1, hard cases — a lost record must be loud; an old byte untouched
# ---------------------------------------------------------------------------

def test_a_record_that_cannot_be_written_is_named_on_the_round(
        fake_glue, tmp_path, monkeypatch, capsys):
    """The root is a FILE, so nothing can be appended under it.

    The round RAN and was paid for; what it must not do is report itself kept.
    ``run()`` returns the record with the refusal on it, and ``main()`` turns
    that into exit 4 (the CLI test below). Raising out of ``run()`` instead
    would throw away the instances, the spend and the reason -- and leave a
    caller unable to tell a lost store from a crash in the round.
    """
    from delfin.agent import bench_history as bh

    blocked = tmp_path / "not-a-directory"
    blocked.write_text("stop")
    monkeypatch.setattr(bh, "default_root", lambda: blocked)

    record = bn.run({"model": "M", "sampler": _Sampler(), "grader": _Grader(),
                     "accounting": _FakeAccounting(), "runner": _runner(),
                     "instances": 2, "mode": "practice", "public_seed": 7})

    assert record.history_error, (
        "the round could not be kept and its record does not say so")
    assert record.instances, "the round ran and returned no instances"


def test_the_cli_says_the_round_was_not_kept(fake_glue, tmp_path,
                                             monkeypatch, capsys):
    """Same failure through `python -m delfin.agent.bench_nightly`: a
    non-zero exit and the words that say the round was lost.

    The rounds file is made a DIRECTORY: appending to it then fails on
    every open mode (an existing regular file would happily accept
    O_APPEND, which is why a file-in-the-way is not a blocked store).
    """
    from delfin.agent import bench_history as bh

    blocked = tmp_path / "kept"
    (blocked / bh.ROUNDS_FILE).mkdir(parents=True)
    monkeypatch.setattr(bh, "default_root", lambda: blocked)

    rc = _round(fake_glue, tmp_path)
    err = capsys.readouterr().err
    assert rc != 0, "a round whose record could not be written exited 0"
    assert "not kept" in err.lower()


def test_an_existing_line_is_never_rewritten(fake_glue, history_dir,
                                             tmp_path):
    """Append-only: the bytes already written stay byte-identical, and the
    second round only adds lines."""
    from delfin.agent import bench_history as bh

    assert _round(fake_glue, tmp_path) == 0
    rounds_path = history_dir / bh.ROUNDS_FILE
    before = rounds_path.read_bytes()
    assert _round(fake_glue, tmp_path) == 0
    after = rounds_path.read_bytes()
    assert after.startswith(before), "an existing byte was rewritten"
    assert len(after) > len(before)
    assert len(_lines(rounds_path)) == 2


def test_a_re_run_of_a_named_round_is_kept_beside_the_first_one(
        fake_glue, history_dir, tmp_path):
    """Two rounds with the SAME id are two measurements: each keeps its
    own instances, and neither pass rate is computed over the union."""
    from delfin.agent import bench_history as bh

    assert bn.run({"model": "M", "sampler": _Sampler(), "grader": _Grader(),
                   "accounting": _FakeAccounting(), "runner": _runner(),
                   "instances": 2, "mode": "practice", "public_seed": 7,
                   "round_id": "named"}) is not None
    assert bn.run({"model": "M", "sampler": _Sampler(),
                   "grader": _Grader(fail=("task0",)),
                   "accounting": _FakeAccounting(), "runner": _runner(),
                   "instances": 2, "mode": "practice", "public_seed": 7,
                   "round_id": "named"}) is not None
    rounds = bh.load(history_dir)
    assert len(rounds) == 2
    assert [len(r.instances) for r in rounds] == [2, 2]
    rates = sorted(r.pass_rate() for r in rounds)
    assert rates == [0.5, 1.0], "the two same-id rounds were merged"


def test_an_unmeasured_round_is_not_counted_as_a_failure(fake_glue,
                                                        history_dir,
                                                        tmp_path, capsys):
    """Every run died in engine init: that measured nothing. Reading it as
    a 0% round would push a permanent false regression into the trend."""
    from delfin.agent import bench_history as bh, cli

    assert _round(fake_glue, tmp_path,
                  unmeasured=("task0", "task1")) == 3   # HB2 exit code
    rounds = bh.load(history_dir)
    assert len(rounds) == 1, "a dead round is still a round: it must be kept"
    assert rounds[0].pass_rate() is None
    capsys.readouterr()
    assert cli.main(["bench", "trend", "--capability",
                     "redox_consistency"]) == 0
    out = capsys.readouterr().out
    assert "UNMEASURED" in out or "measured nothing" in out


# ---------------------------------------------------------------------------
# Item 3 — scheduled, off by default, budget gate kept
# ---------------------------------------------------------------------------

@pytest.fixture
def fake_systemd(tmp_path, monkeypatch):
    """Shadow systemctl and the unit directory.

    The host's systemd is never asked what it has: that would make the
    test report the machine, and on a host with a user manager it would
    install units into the operator's real ~/.config.
    """
    from delfin.agent import bench_history as bh

    calls = []

    class _Done:
        def __init__(self, rc=0, stdout="active", stderr=""):
            self.returncode, self.stdout, self.stderr = rc, stdout, stderr

    def fake_systemctl(*args, ok=True):
        calls.append(list(args))
        return _Done(0 if ok else 1)

    monkeypatch.setattr(bh, "_systemctl",
                        lambda *a: fake_systemctl(*a))
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    return calls


def test_no_round_is_scheduled_until_asked(fake_glue, fake_systemd, capsys):
    """Off by default: neither importing the module nor running the
    commands registers anything."""
    from delfin.agent import bench_history as bh, cli

    assert bh.schedule_status()["installed"] is False
    assert cli.main(["bench", "schedule"]) == 0
    out = capsys.readouterr().out
    assert "off" in out.lower()
    assert bh.schedule_status()["installed"] is False
    assert not any(c[0] == "enable" for c in fake_systemd)


def test_the_scheduled_command_keeps_the_budget_gate(fake_glue, fake_systemd,
                                                     capsys):
    """A fired round must refuse instances over the limit exactly like a
    hand-run one: the gate is part of the command it runs."""
    from delfin.agent import bench_history as bh, cli

    assert cli.main(["bench", "schedule", "--enable", "--model", "M",
                     "--budget", "7.5", "--instances", "3"]) == 0
    unit = (bh.unit_dir() / bh.UNIT).read_text()
    assert "--budget" in unit and "7.5" in unit
    assert "--instances" in unit and "3" in unit
    assert "-m delfin.agent.bench_nightly" in unit.replace('"', "")


def test_enabling_installs_a_timer_and_disabling_removes_it(fake_glue,
                                                            fake_systemd,
                                                            capsys):
    from delfin.agent import bench_history as bh, cli

    assert cli.main(["bench", "schedule", "--enable", "--model", "M"]) == 0
    assert bh.schedule_status()["installed"] is True
    assert (bh.unit_dir() / bh.TIMER).exists()
    assert ["enable", "--now", bh.TIMER] in fake_systemd
    capsys.readouterr()
    assert cli.main(["bench", "schedule", "--disable"]) == 0
    assert bh.schedule_status()["installed"] is False


def test_scheduling_without_a_model_is_refused(fake_glue, fake_systemd,
                                               capsys):
    """An ExecStart that cannot run must not be enabled."""
    from delfin.agent import bench_history as bh, cli

    assert cli.main(["bench", "schedule", "--enable"]) == 2
    assert bh.schedule_status()["installed"] is False


def test_scheduling_refuses_where_there_is_no_user_manager(fake_glue,
                                                           tmp_path,
                                                           monkeypatch,
                                                           capsys):
    """Writing units no supervisor will read is a round that silently never
    runs — the same failure this package is about."""
    from delfin.agent import bench_history as bh, cli

    class _Done:
        returncode, stdout, stderr = 1, "", "Failed to connect to bus"

    monkeypatch.setattr(bh, "_systemctl", lambda *a: _Done())
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))
    assert cli.main(["bench", "schedule", "--enable", "--model", "M"]) == 1
    assert "NOT scheduled" in capsys.readouterr().err
    assert bh.schedule_status()["installed"] is False


# ---------------------------------------------------------------------------
# The store never lives in a repository
# ---------------------------------------------------------------------------

def test_the_history_root_refuses_to_live_inside_a_repository(tmp_path,
                                                             monkeypatch):
    """A record committed beside the code would show up in every
    `git status` a session runs and get committed by someone."""
    from delfin.agent import bench_history as bh

    (tmp_path / ".git").mkdir()
    (tmp_path / ".git" / "HEAD").write_text("ref: refs/heads/main\n")
    monkeypatch.setattr(bh, "default_root",
                        lambda: tmp_path / "checkout" / "rounds")
    with pytest.raises(bh.HistoryWriteError) as exc:
        bh.root()
    assert "repository" in str(exc.value)


def test_an_explicit_root_in_a_repository_is_refused_too(tmp_path):
    """The same write arriving by a different argument obeys the same rule.
    """
    from delfin.agent import bench_history as bh

    (tmp_path / ".git").mkdir()
    (tmp_path / ".git" / "HEAD").write_text("ref: refs/heads/main\n")
    with pytest.raises(bh.HistoryWriteError):
        bh.write_round({"round_id": "x"},
                       history_root=tmp_path / "checkout" / "r")


def test_an_empty_dot_git_directory_is_not_a_repository(tmp_path):
    """An empty ``.git`` left in a shared parent is not a work tree; a
    worktree's ``gitdir:`` file is one."""
    from delfin.agent import bench_history as bh

    (tmp_path / ".git").mkdir()
    assert bh._git_toplevel(tmp_path / "a" / "r") == ""
    wt = tmp_path / "wt"
    wt.mkdir()
    (wt / ".git").write_text("gitdir: /nowhere\n")
    assert bh._git_toplevel(wt / "r") == str(wt)


def test_the_history_store_is_registered_as_a_redirectable_user_state():
    """The store must be in the product's user-state table, with the leaf
    the resolver actually uses.

    This is what keeps a test run's fixture rounds out of the history the
    operator reads trends from: state_paths.USER_STATE_RESOLVERS is the
    one table the suite and RedirectedUserState both follow, and an
    append-only store is one nothing later cleans. If the module or the
    leaf name drifts from the entry, the redirect sends a run's writes to
    the real store while the code writes somewhere else — so the drift
    fails here, not silently in the operator's history.
    """
    from delfin.agent import bench_history as bh, state_paths

    entry = [e for e in state_paths.USER_STATE_RESOLVERS
             if e[0] == "delfin.agent.bench_history"]
    assert entry, ("bench_history.default_root is not in "
                   "USER_STATE_RESOLVERS: a test would append fixture "
                   "rounds to the operator's real history")
    assert entry[0][1] == "default_root"
    assert entry[0][2] == bh.HISTORY_LEAF, (
        f"the table redirects {entry[0][2]!r} but the store lives in "
        f"{bh.HISTORY_LEAF!r}")
    assert callable(bh.default_root), (
        "the redirect replaces a CALLABLE resolver; a constant would be "
        "skipped and the writes would go to the real home")


# ---------------------------------------------------------------------------
# Controls: these pass on origin/main AND after the change. They are the
# never-worse cases — the seed gate and the budget gate must survive the
# store being wired into run().
# ---------------------------------------------------------------------------

def test_CONTROL_a_holdout_round_without_a_seed_file_is_still_refused(
        fake_glue, history_dir, tmp_path, capsys):
    rc = bn.run({"model": "M", "sampler": _Sampler(), "grader": _Grader(),
                 "accounting": _FakeAccounting(), "runner": _runner(),
                 "instances": 2, "mode": "holdout",
                 "holdout_seed_path": str(tmp_path / "absent.txt")})
    assert rc.holdout_refusal_reason
    assert rc.instances == []


def test_CONTROL_the_budget_gate_still_refuses_over_the_daily_limit(
        fake_glue, history_dir):
    class _Expensive(_FakeAccounting):
        def estimate(self, task):
            return 30.0

    rec = bn.run({"model": "M", "sampler": _Sampler(), "grader": _Grader(),
                  "accounting": _Expensive(), "runner": _runner(),
                  "instances": 2, "mode": "practice", "public_seed": 7,
                  "daily_budget_core_hours": 50.0})
    started = [i for i in rec.instances if i.ran]
    refused = [i for i in rec.instances if not i.ran]
    assert len(started) == 1 and len(refused) == 1
    assert "over daily budget" in refused[0].skip_reason


def test_the_seed_value_never_reaches_the_kept_history(fake_glue, history_dir,
                                                      tmp_path):
    """Not a new contract — the seed rule already binds the record, held
    there by test_pl4_budget_and_seed.py::test_seed_never_in_outputs and
    ::test_seed_path_not_in_record_or_instances. This one is red on main
    only because main keeps no history to leak into: it checks the NEW
    writer did not become the leak. The round line must carry the seed
    MODE, never the seed."""
    assert _round(fake_glue, tmp_path, practice=False) == 0
    blob = "".join(p.read_text(encoding="utf-8")
                   for p in history_dir.glob("*.jsonl"))
    assert "424242" not in blob
    assert '"seed_mode"' in blob


def test_the_no_systemd_message_names_a_command_that_works(fake_glue,
                                                           monkeypatch,
                                                           tmp_path,
                                                           capsys):
    """The refusal must not send the reader to a flag nobody implements.

    Measured on this branch BEFORE the fix: with no systemd user manager,
    `bench schedule --enable` said "the command is in `bench schedule
    --show`" — and `--show` is not an argument of `bench schedule` at all.
    The parser registers only --enable/--disable/--model/--provider/
    --instances/--budget/--practice/--seed-path/--on-calendar
    (cli.py:6267-6285), so following the message produced "unrecognized
    arguments: --show" and the reader was worse off than before.

    The message now names the bare command. This test holds it there AND
    checks that the bare command really is accepted, so the help text
    cannot drift away from the parser again.
    """
    from delfin.agent import bench_history as bh, cli

    class _NoBus:
        returncode = 1
        stdout = ""
        stderr = "Failed to connect to bus: Permission denied"

    # The host's systemd is never asked: the probe is shadowed.
    monkeypatch.setattr(bh, "_systemctl", lambda *a: _NoBus())
    monkeypatch.setenv("XDG_CONFIG_HOME", str(tmp_path / "config"))

    with pytest.raises(bh.HistoryWriteError) as exc:
        bh.install_schedule(model="M")
    message = str(exc.value)
    assert "no systemd user manager" in message

    # What the message points at must be a command the parser accepts.
    assert "--show" not in message
    capsys.readouterr()
    assert cli.main(["bench", "schedule"]) == 0
    out = capsys.readouterr().out
    assert "benchmark round schedule" in out


def test_the_history_header_names_no_round_it_does_not_have(fake_glue,
                                                            history_dir,
                                                            tmp_path, capsys):
    """A header that prints `round ?` is a header with nothing to say.

    RED before the fix: two rounds whose every instance is UNMEASURED are
    both excluded, `trend()["reference"]` is None, and the header read
    `benchmark history — comparable with round ?`. The `?` is not a round
    id: no round is called that, so the line named a comparison that did
    not exist. A reader cannot tell `?` from a round id they have not
    learned to recognise yet.

    The assertion is on the SET of round ids the store actually holds, so
    it fails for any placeholder, and it also fails if the header stops
    naming a round when there IS one comparable.
    """
    from delfin.agent import bench_history as bh
    from delfin.agent import cli

    # Every instance unmeasured -> both rounds excluded -> no reference.
    for _ in range(2):
        assert _round(fake_glue, tmp_path, practice=True,
                      unmeasured=("task0", "task1")) == 3
    capsys.readouterr()

    assert cli.main(["bench", "history"]) == 0
    out = capsys.readouterr().out
    header = out.splitlines()[0]
    stored_ids = {r["round_id"] for r in _lines(history_dir / bh.ROUNDS_FILE)}
    assert stored_ids, "the fixture wrote no round at all"
    assert "?" not in header, header

    # CONTROL: with a comparable round present the header still names it.
    assert _round(fake_glue, tmp_path, practice=True) == 0
    capsys.readouterr()
    assert cli.main(["bench", "history"]) == 0
    header = capsys.readouterr().out.splitlines()[0]
    stored_ids = {r["round_id"] for r in _lines(history_dir / bh.ROUNDS_FILE)}
    named = [rid for rid in stored_ids if rid in header]
    assert len(named) == 1, (header, sorted(stored_ids))
