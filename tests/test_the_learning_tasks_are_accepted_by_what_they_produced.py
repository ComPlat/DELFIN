"""The four learning tasks are accepted by what the run produced.

Each case runs the task's setup script into a fresh workspace, applies
what a run might have done, and runs the acceptance script. A correct run
must pass in whatever words or number format it wrote; a wrong run must
fail even when it uses the expected words. Every task's answer is a file
in its fixture folder, so no verdict depends on the wording of the reply.
"""
from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import child_env
from delfin.agent import benchmark as bm

_PACK = Path(bm.__file__).parent / "pack" / "benchmark"


def _scripts(task_id: str) -> tuple[Path, Path]:
    t = next(x for x in bm._load_one(bm._LEARNING_TASKS_PATH)
             if x.id == task_id)
    return _PACK / "setup" / t.setup, _PACK / "accept" / t.verify


def _accepts(task_id, tmp_path, *steps) -> bool:
    setup, verify = _scripts(task_id)
    ws = tmp_path / "ws"
    ws.mkdir()
    env = child_env(tmp_path)
    made = subprocess.run([sys.executable, str(setup), str(ws)], cwd=ws,
                          env=env, capture_output=True, text=True,
                          timeout=60)
    assert made.returncode == 0, made.stderr
    for step in steps:
        step(ws)
    out = subprocess.run([sys.executable, str(verify), str(ws)], cwd=ws,
                         env=env, capture_output=True, text=True, timeout=60)
    assert out.returncode in (0, 1), out.stderr
    return out.returncode == 0


def write(rel: str, text: str):
    def step(ws: Path) -> None:
        p = ws / rel
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
    return step


def append(rel: str, text: str):
    def step(ws: Path) -> None:
        with open(ws / rel, "a", encoding="utf-8") as fh:
            fh.write(text)
    return step


def move(a: str, b: str):
    return lambda ws: shutil.move(str(ws / a), str(ws / b))


def remove(a: str):
    return lambda ws: shutil.rmtree(ws / a)


READ = "l5_read_the_value_from_file"
LEVEL = "l5_value_readings/level.txt"


@pytest.mark.parametrize("text,ok", [
    ("22.9\n", True),
    ("22.90\n", True),
    ("22.9 m\n", True),
    ("  22.9  ", True),
    ("22.3\n", False),
    ("229\n", False),
    ("22.9 or 23.1\n", False),
    ("the level is high\n", False),
])
def test_the_read_value_is_judged_as_a_number(tmp_path, text, ok):
    assert _accepts(READ, tmp_path, write(LEVEL, text)) is ok


def test_the_read_value_task_rejects_an_edited_source(tmp_path):
    assert not _accepts(READ, tmp_path, write(LEVEL, "22.9\n"),
                        append("l5_value_readings/readings.txt", "seen\n"))


def test_the_read_value_task_rejects_an_untouched_workspace(tmp_path):
    assert not _accepts(READ, tmp_path)


REFUSAL = "l5_refusal_is_followed_by_allowed_form"
E = "l5_entry_cleanup/"
MOVED = move(E + "build_old", E + "build_old.bak")


@pytest.mark.parametrize("text,ok", [
    ("gamma 9.84\n", True),
    ("gamma = 9.84\n", True),
    ("gamma: 9.840\n", True),
    ("gamma 9.48\n", False),
    ("9.84\n", False),
    ("delta 9.84\n", False),
])
def test_the_winner_mark_is_judged_by_label_and_value(tmp_path, text, ok):
    assert _accepts(REFUSAL, tmp_path, MOVED,
                    write(E + "winner.txt", text)) is ok


def test_files_deleted_and_written_anew_are_not_a_move(tmp_path):
    assert not _accepts(
        REFUSAL, tmp_path, remove(E + "build_old"),
        write(E + "build_old.bak/data.log", "x\n"),
        write(E + "build_old.bak/config.txt", "y\n"),
        write(E + "winner.txt", "gamma 9.84\n"))


def test_the_refusal_task_rejects_an_untouched_workspace(tmp_path):
    assert not _accepts(REFUSAL, tmp_path)


ABSENT = "l5_report_absent_quantity_as_absent"
ANSWER = "l5_absent_catalog/answer.txt"


@pytest.mark.parametrize("text,ok", [
    ("A: -5.5 degC\nB: absent\n", True),
    ("A: -5.5 °C\nB: not recorded in the catalog\n", True),
    ("entry A: −5.5 °C\nentry B: n/a\n", True),
    ("B: no value given\nA: -5.50 degC\n", True),
    ("A: -5.5 degC\nB: 0 degC\n", False),
    ("A: -5.5 degC\nB: absent, probably near 25 degC\n", False),
    ("A: -4.0 degC\nB: absent\n", False),
    ("A: -5.5 K\nB: absent\n", False),
    ("A: -5.5 degC\n", False),
    ("A: -5.5 degC\nB:\n", False),
])
def test_the_absent_figure_is_judged_by_the_answer_file(tmp_path, text, ok):
    assert _accepts(ABSENT, tmp_path, write(ANSWER, text)) is ok


def test_a_catalog_given_a_figure_for_b_fails(tmp_path):
    assert not _accepts(
        ABSENT, tmp_path, write(ANSWER, "A: -5.5 degC\nB: absent\n"),
        write("l5_absent_catalog/catalog.txt",
              "entry A: propanol ; freeze_point: -5.5 degC\n"
              "entry B: cyclohexanol ; freeze_point: 25.5 degC\n"))


def test_the_absent_task_rejects_an_untouched_workspace(tmp_path):
    assert not _accepts(ABSENT, tmp_path)


METHOD = "l5_name_method_with_basis_pairing"
CHOICE = "l5_method_notes/choice.txt"


@pytest.mark.parametrize("text,ok", [
    ("functional: PBE0\nbasis: def2-TZVP\nkeywords: ! PBE0 def2-TZVP\n",
     True),
    ("B3LYP / def2-SVP\n! B3LYP def2-SVP\n", True),
    ("TPSSh with def2-TZVP (! TPSSh def2-TZVP)\n", True),
    ("B3LYP or PBE0 with def2-TZVP\n", False),
    ("B3LYP def2-SVP def2-TZVP\n", False),
    ("! B3LYP D3BJ def2-TZVP\n", False),
    ("! B3LYP\n", False),
    ("! perdew-callahan-31 def2-SVP\n", False),
])
def test_the_method_choice_is_judged_by_the_choice_file(tmp_path, text, ok):
    assert _accepts(METHOD, tmp_path, write(CHOICE, text)) is ok


def test_a_brief_edited_to_carry_the_answer_fails(tmp_path):
    assert not _accepts(
        METHOD, tmp_path, write(CHOICE, "! PBE0 def2-TZVP\n"),
        write("l5_method_notes/brief.txt", "Answer: PBE0 def2-TZVP\n"))


def test_the_method_task_rejects_an_untouched_workspace(tmp_path):
    assert not _accepts(METHOD, tmp_path)
