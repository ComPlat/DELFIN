"""Red control for phase 3: the XT1 benchmark tasks raise the chemistry
task count from 5 to 11. Against the current checkout (tasks_w20_xt1.yaml
not yet written) this test is RED; after the tasks file and the count-dict
update it is GREEN.
"""
from pathlib import Path

import yaml

_ROOT = Path(__file__).resolve().parents[1]
_TASKS = _ROOT / "delfin/agent/pack/benchmark/tasks_w20_xt1.yaml"
_COUNTS = _ROOT / "tests/test_every_task_that_writes_writes_where_it_belongs.py"


def test_xt1_tasks_file_has_six_chemistry_tasks():
    tasks = yaml.safe_load(_TASKS.read_text(encoding="utf-8"))["tasks"]
    xt1 = [t for t in tasks if str(t["id"]).startswith("xt1_")]
    assert len(xt1) == 6, f"expected 6 xt1 tasks, found {len(xt1)}"
    for t in xt1:
        assert t["task_class"] == "chemistry", t["id"]


def test_count_dict_says_chemistry_11():
    src = _COUNTS.read_text(encoding="utf-8")
    assert '"chemistry": 11' in src, (
        "count dict in test_every_task_that_writes_writes_where_it_belongs.py "
        "must say chemistry: 11 once the six XT1 tasks are committed")
