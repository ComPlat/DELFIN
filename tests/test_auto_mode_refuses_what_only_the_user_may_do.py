"""An auto iteration publishes nothing, widens nothing, and leaves the measure.

Input: one tool call (bare or namespaced name, arguments) made while an
auto iteration runs, the loop's measure and touch lists. Output of
``auto_guard.classify_call``: None (the call goes on to the ordinary gate)
or a Verdict whose ``kind`` is ``outward`` (push, pull request, merge,
branch deletion, message or upload leaving the machine, a future turn
outside the budget), ``widen`` (a permission or a writable root added),
``measure`` (a write to a file the iteration is judged by) or ``scope`` (a
write outside the loop's ``touch`` list). Outward and widen verdicts are
queued for the user (``queue``), never run.
"""

from __future__ import annotations

import pytest

from delfin.agent import auto_guard as ag


def _c(name, args=None, **kw):
    kw.setdefault("workspace", "/w")
    return ag.classify_call(name, args or {}, **kw)


@pytest.mark.parametrize("cmd", [
    "git push",
    "git push origin loop/x",
    "git -C /w push --force-with-lease",
    "cd /w && git push -u origin HEAD",
    "bash -c 'git push origin main'",
    "gh pr create --fill",
    "gh pr merge 12 --squash",
    "gh pr comment 12 -b done",
    "gh issue create -t x -b y",
    "gh release create v1",
    "gh api -X POST repos/o/r/issues -f title=x",
    "git branch -D loop/old",
    "git branch --delete feature",
    "git push origin --delete feature",
    "git merge main",
    "git switch main",
    "git checkout master",
    "curl -X POST https://hooks.example/x -d @results.csv",
    "curl --data-binary @a https://x",
    "curl -F file=@a https://x",
    "wget --post-data=a=b https://x",
    "mail -s hi someone@example.org < note.txt",
    "sendmail someone@example.org",
    "scp results.csv me@cluster:/tmp/",
    "rsync -a out/ me@host:backup/",
    "ssh login1 'sbatch job.sh'",
])
def test_every_way_out_is_queued(cmd):
    v = _c("bash", {"command": cmd})
    assert v is not None and v.kind == "outward" and v.queue, cmd
    # Same command through a background shell and an MCP shell.
    assert _c("bash_background", {"command": cmd}).kind == "outward"
    assert _c("mcp__kit-coding__bash", {"command": cmd}).kind == "outward"


@pytest.mark.parametrize("cmd", [
    "git status",
    "git diff HEAD~1",
    "git log --oneline -5",
    "git commit -m 'loop: row added'",
    "git add results/energies.csv",
    "git switch -c loop/x",
    "git switch loop/x",
    "pytest -q tests/test_table.py",
    "python scripts/parse_energy.py out/job1",
    "curl -s https://example.org/data.json",
    "sbatch job.sh",
    "squeue -u $USER",
])
def test_local_work_passes(cmd):
    assert _c("bash", {"command": cmd}) is None, cmd


@pytest.mark.parametrize("name", [
    "session_message", "push_notification", "remote_trigger",
    "cron_create", "schedule_wakeup", "worktree_merge",
    "mcp__github__create_pull_request", "mcp__slack__send_message",
    "mcp__github__merge_pull_request",
])
def test_outward_tools_are_queued(name):
    v = _c(name, {"text": "x"})
    assert v is not None and v.kind == "outward" and v.queue, name


@pytest.mark.parametrize("name", [
    "remember_permission", "remember_permission_bundle",
    "exit_plan_mode", "enter_worktree",
])
def test_widening_tools_are_queued(name):
    v = _c(name, {})
    assert v is not None and v.kind == "widen" and v.queue


@pytest.mark.parametrize("name", [
    "read_file", "grep_file", "list_files", "task_update", "web_search",
    "mcp__github__get_issue", "subagent", "run_tests",
])
def test_reads_and_bookkeeping_pass(name):
    assert _c(name, {"path": "x"}) is None


MEASURE = ["tests/test_table.py", "scripts/check_table.py",
           ".delfin/loops/cluster-table.md"]


@pytest.mark.parametrize("name, args", [
    ("write_file", {"path": "tests/test_table.py", "content": "x"}),
    ("edit_file", {"path": "/w/scripts/check_table.py", "old_string": "a",
                   "new_string": "b"}),
    ("multi_edit", {"path": "./tests/test_table.py", "edits": []}),
    ("write_file", {"path": ".delfin/loops/cluster-table.md",
                    "content": "budget: {iterations: 500}"}),
    ("apply_patch", {"patch": "--- a/tests/test_table.py\n"
                              "+++ b/tests/test_table.py\n@@\n-x\n+y\n"}),
    ("bash", {"command": "sed -i 's/== 3/== 4/' tests/test_table.py"}),
    ("bash", {"command": "echo 'exit 0' > scripts/check_table.py"}),
    ("bash", {"command": "rm tests/test_table.py"}),
])
def test_the_measure_is_not_written(name, args):
    v = _c(name, args, measure=MEASURE)
    assert v is not None and v.kind == "measure", (name, args)
    assert not v.queue


def test_reading_the_measure_passes():
    assert _c("read_file", {"path": "tests/test_table.py"},
              measure=MEASURE) is None
    assert _c("bash", {"command": "pytest -q tests/test_table.py"},
              measure=MEASURE) is None
    assert _c("bash", {"command": "cat scripts/check_table.py"},
              measure=MEASURE) is None


def test_touch_limits_writes():
    touch = ["results/", "scripts/parse_energy.py"]
    assert _c("write_file", {"path": "results/energies.csv"},
              touch=touch) is None
    assert _c("edit_file", {"path": "scripts/parse_energy.py"},
              touch=touch) is None
    v = _c("write_file", {"path": "delfin/core.py"}, touch=touch)
    assert v is not None and v.kind == "scope" and not v.queue
    # .git is never in scope.
    v = _c("write_file", {"path": ".git/config"})
    assert v is not None and v.kind == "scope"


def test_mode_rank_orders_postures():
    r = ag.mode_rank
    assert r("plan") < r("default") <= r("diff_approval") < r("acceptEdits") \
        < r("bypassPermissions")
    assert ag.widened("default", "acceptEdits")
    assert ag.widened("acceptEdits", "bypassPermissions")
    assert not ag.widened("acceptEdits", "plan")
    assert not ag.widened("default", "default")


def test_measure_findings_hard_and_soft(tmp_path):
    changes = {
        "scripts/check_table.py": ("exit 1\n", "exit 0\n"),
        "tests/test_other.py": ("def test_a():\n    assert f() == 3\n",
                                "import pytest\n@pytest.mark.skip\n"
                                "def test_a():\n    pass\n"),
        "results/energies.csv": ("a,1\n", "a,1\nb,2\n"),
    }
    found = ag.measure_findings(changes, measure=MEASURE)
    by_path = {f.path: f for f in found}
    assert by_path["scripts/check_table.py"].hard
    assert "tests/test_other.py" in by_path
    assert not by_path["tests/test_other.py"].hard
    assert "results/energies.csv" not in by_path


def test_measure_findings_uses_edit_signals_when_present(monkeypatch):
    import sys
    import types

    calls = []
    fake = types.ModuleType("delfin.agent.edit_signals")

    def classify(path, old, new):
        calls.append(path)
        return {"measure": [("threshold", "tolerance loosened")]}

    fake.classify = classify
    monkeypatch.setitem(sys.modules, "delfin.agent.edit_signals", fake)
    found = ag.measure_findings({"conf/limits.yaml": ("a: 1\n", "a: 9\n")},
                                measure=[])
    assert calls == ["conf/limits.yaml"]
    assert found and found[0].detail == "tolerance loosened"
    assert ag.measure_source() == "edit_signals"
