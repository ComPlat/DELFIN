"""M1 phase-3: the yaml suite loads, is complete, and its scripts exist.

load_tasks() picks up every tasks_*.yaml sibling automatically, so these
six tasks appear in the packaged suite without registration. This test
pins what that pickup must hold: six ids, each with prompt_en, local
setup/verify script paths that resolve inside the pack, and no network
dependency (the fixtures are all git-seeded).
"""
from pathlib import Path

from delfin.agent import benchmark as bm

_PACK = Path(bm.__file__).resolve().parent / "pack" / "benchmark"

IDS = ("w17_method_code_red",
       "w17_method_code_feature",
       "w17_method_data_audit",
       "w17_method_data_convert",
       "w17_method_web_manual",
       "w17_method_web_changelog")


def _task(task_id: str):
    for t in bm.load_tasks():
        if t.id == task_id:
            return t
    raise AssertionError(f"{task_id} not loaded from the packaged suite")


def test_all_six_tasks_load_from_the_packaged_suite():
    ids = {t.id for t in bm.load_tasks()}
    missing = [i for i in IDS if i not in ids]
    assert not missing, f"missing from suite: {missing}"


def test_every_task_has_its_scripts_inside_the_pack():
    for tid in IDS:
        t = _task(tid)
        if t.setup:
            p = (_PACK / "setup" / t.setup).resolve()
            assert p.is_file(), f"{tid}: setup missing {p}"
            assert str(p).startswith(str((_PACK / "setup").resolve()))
        p = (_PACK / "accept" / t.verify).resolve()
        assert p.is_file(), f"{tid}: verify missing {p}"
        assert str(p).startswith(str((_PACK / "accept").resolve()))


def test_every_task_prompt_names_the_brief_and_offline():
    for tid in IDS:
        t = _task(tid)
        assert t.prompt_en, f"{tid}: prompt_en empty"
        assert "m1_brief.txt" in t.prompt_en, f"{tid}: prompt misses brief"
        assert "offline" in t.prompt_en.lower(), f"{tid}: prompt misses offline"


def test_every_task_grades_by_budgets():
    for tid in IDS:
        t = _task(tid)
        assert t.max_tool_calls >= 25, f"{tid}: too few tool calls allowed"
        assert t.max_duration_s >= 300, f"{tid}: too short"
        assert t.max_cost_usd >= 0.5, f"{tid}: budget too small"
