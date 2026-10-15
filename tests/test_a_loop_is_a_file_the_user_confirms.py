"""A loop is a small declarative file, validated, and confirmed by its hash.

Input: ``.delfin/loops/<name>.md`` in a workspace or ``~/.delfin/loops/``
(YAML front matter, numbered steps in the body). Output: a ``LoopSpec`` or
the list of what is wrong with the file. A loop runs only after the user
confirmed exactly this content: an edit after the confirmation (by anyone,
the agent included) withdraws it.
"""

from __future__ import annotations

from pathlib import Path

import pytest

from delfin.agent import loops


GOOD = """---
name: cluster-table
goal: Every finished cluster job has its final energy in results/energies.csv
trigger: job
check: python scripts/check_table.py
metric: 'missing=(\\d+)'
better: lower
budget:
  iterations: 12
  wall: 6h
  tokens: 900000
stop:
  plateau: 2
  repeat_failures: 2
touch: [results/, scripts/parse_energy.py]
measure: [tests/test_table.py]
---
# Steps
1. Read the output folder of the job the wake names.
2. Extract the final single-point energy.
   run: python scripts/parse_energy.py <folder>
3. Add or update its row in results/energies.csv.
   check: python scripts/check_table.py
"""


@pytest.fixture()
def ws(tmp_path, monkeypatch):
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: home))
    w = tmp_path / "ws"
    (w / ".delfin" / "loops").mkdir(parents=True)
    return w


def _write(ws: Path, name: str, text: str) -> Path:
    p = ws / ".delfin" / "loops" / f"{name}.md"
    p.write_text(text, encoding="utf-8")
    return p


def test_a_good_file_parses_into_every_field(ws):
    _write(ws, "cluster-table", GOOD)
    spec = loops.load("cluster-table", ws)
    assert spec.errors == []
    assert spec.name == "cluster-table"
    assert spec.trigger == "job" and spec.every_seconds == 0
    assert [s.text.split()[0] for s in spec.steps] == ["Read", "Extract", "Add"]
    assert spec.steps[1].run == "python scripts/parse_energy.py <folder>"
    assert spec.steps[2].check == "python scripts/check_table.py"
    assert spec.max_iterations == 12
    assert spec.wall_seconds == 6 * 3600
    assert spec.max_tokens == 900000
    assert spec.plateau == 2 and spec.repeat_failures == 2
    assert spec.branch == "loop/cluster-table"
    assert spec.touch == ["results/", "scripts/parse_energy.py"]
    # The measure: the declared files, the loop file itself, and the
    # file the check command runs.
    assert "tests/test_table.py" in spec.measure
    assert ".delfin/loops/cluster-table.md" in spec.measure
    assert "scripts/check_table.py" in spec.measure
    assert spec.metric_value("ok missing=3 rows") == 3.0


def test_defaults_are_bounded(ws):
    _write(ws, "tidy", "---\ngoal: tidy the notes\n---\n1. Tidy one file.\n")
    spec = loops.load("tidy", ws)
    assert spec.errors == []
    assert spec.trigger == "manual"
    assert spec.max_iterations == loops.DEFAULT_ITERATIONS
    assert spec.wall_seconds == loops.DEFAULT_WALL_S
    assert spec.max_tokens == loops.DEFAULT_TOKENS


@pytest.mark.parametrize("front, needle", [
    ("goal: x\nbranch: main", "branch"),
    ("goal: x\nbranch: master", "branch"),
    ("goal: x\ntrigger: sometimes", "trigger"),
    ("goal: x\ntrigger: every 5", "trigger"),
    ("goal: x\nbudget:\n  iterations: 0", "iterations"),
    ("goal: x\nbudget:\n  iterations: 100000", "iterations"),
    ("goal: x\ntouch: [/etc/passwd]", "touch"),
    ("goal: x\ntouch: [../outside]", "touch"),
    ("goal: x\nmetric: '(['", "metric"),
    ("goal: x\nbetter: sideways", "better"),
    ("trigger: idle", "goal"),
    ("goal: x\nname: Bad Name!", "name"),
])
def test_what_is_wrong_is_said(ws, front, needle):
    _write(ws, "bad", f"---\n{front}\n---\n1. Do it.\n")
    spec = loops.load("bad", ws)
    assert spec.errors, front
    assert any(needle in e for e in spec.errors), spec.errors


def test_a_file_without_steps_is_refused(ws):
    _write(ws, "empty", "---\ngoal: nothing\n---\nNo numbered steps here.\n")
    assert any("step" in e for e in loops.load("empty", ws).errors)


def test_every_interval_trigger(ws):
    _write(ws, "poll", "---\ngoal: x\ntrigger: every 30m\n---\n1. Look.\n")
    spec = loops.load("poll", ws)
    assert spec.errors == [] and spec.trigger == "every"
    assert spec.every_seconds == 1800


def test_workspace_wins_over_personal_and_both_are_listed(ws):
    personal = Path.home() / ".delfin" / "loops"
    personal.mkdir(parents=True)
    (personal / "tidy.md").write_text(
        "---\ngoal: personal\n---\n1. a\n", encoding="utf-8")
    (personal / "mine.md").write_text(
        "---\ngoal: only mine\n---\n1. a\n", encoding="utf-8")
    _write(ws, "tidy", "---\ngoal: project\n---\n1. a\n")
    names = {s.name: s for s in loops.discover(ws)}
    assert set(names) == {"tidy", "mine"}
    assert names["tidy"].goal == "project"
    assert names["mine"].scope == "personal"


def test_a_loop_is_confirmed_for_its_content_only(ws):
    p = _write(ws, "cluster-table", GOOD)
    spec = loops.load("cluster-table", ws)
    assert not loops.is_confirmed(spec)
    loops.confirm(spec)
    assert loops.is_confirmed(loops.load("cluster-table", ws))
    # Any edit afterwards -- the agent's included -- withdraws it.
    p.write_text(GOOD.replace("iterations: 12", "iterations: 400"),
                 encoding="utf-8")
    assert not loops.is_confirmed(loops.load("cluster-table", ws))


def test_an_invalid_file_cannot_be_confirmed(ws):
    _write(ws, "bad", "---\ngoal: x\nbranch: main\n---\n1. a\n")
    with pytest.raises(ValueError):
        loops.confirm(loops.load("bad", ws))


def test_the_ledger_is_outside_the_workspace(ws):
    _write(ws, "cluster-table", GOOD)
    spec = loops.load("cluster-table", ws)
    loops.confirm(spec)
    ledger = loops.ledger_path()
    assert ledger.is_file()
    assert ws not in ledger.parents


def test_a_template_parses_and_is_not_confirmed(ws):
    path = loops.write_template(ws, "weekly-report",
                                goal="The weekly table is current")
    assert path == ws / ".delfin" / "loops" / "weekly-report.md"
    spec = loops.load("weekly-report", ws)
    assert spec.errors == []
    assert not loops.is_confirmed(spec)
    with pytest.raises(FileExistsError):
        loops.write_template(ws, "weekly-report", goal="again")


def test_an_ad_hoc_loop_from_a_command_line(ws):
    spec = loops.ad_hoc("30m", "run the suite and fix the first failure", ws)
    assert spec.errors == [] and spec.trigger == "every"
    assert spec.every_seconds == 1800 and len(spec.steps) == 1
    assert loops.is_confirmed(spec)      # typed by the user: confirmed
