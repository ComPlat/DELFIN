"""Red control for CB0B phase 2: the workflow tasks file.

tasks_w20_workflows.yaml carries CONCRETE practice tasks derived from the
holdout sampler's public-seed draws (two per family, round-robin), in the
standard benchmark task format, classed ``chemistry`` so the runner
confines their writes to the chemistry workspace fixture.

Written BEFORE the file and its setup scripts existed: the task-existence
and count cases are red on the previous commit.
"""
from __future__ import annotations

import pathlib
import re

import pytest
import yaml

from delfin.agent import bench_holdout as bh

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"
_FILE = _PACK / "tasks_w20_workflows.yaml"

WORKFLOW_WORDS = re.compile(
    r"(?i)\b(classic|redox|esd|occupier|crest|conformer.screen"
    r"|screening.workfl)\b")

TEST_SEED = bh.PUBLIC_SEED


def _tasks():
    data = yaml.safe_load(_FILE.read_text(encoding="utf-8"))
    return data.get("tasks", [])


def test_the_file_exists_and_carries_tasks():
    assert len(_tasks()) > 0


def test_one_concrete_task_per_family_from_the_public_seed():
    """The practice set: the FIRST public-seed draw of each family.
    The holdout stays the generator; this file pins the concrete
    practice instances a runner can execute without a cluster."""
    by_family: dict[str, int] = {}
    for t in _tasks():
        parts = t["id"].split("_")
        assert parts[0] == "w20"
        key = "_".join(parts[1:-1])
        by_family[key] = by_family.get(key, 0) + 1
    assert by_family == {f.id: 1 for f in bh.FAMILIES}


def test_each_prompt_is_the_sampler_prompt_with_fixture_pointer():
    drawn = {i.family_id: i for i in bh.sample(TEST_SEED, len(bh.FAMILIES))}
    for t in _tasks():
        key = "_".join(t["id"].split("_")[1:-1])
        inst = drawn[key]
        assert inst.prompt in t["prompt"], t["id"]
        assert "w20_" in t["prompt"], t["id"]  # fixture folder named


def test_no_task_prompt_names_the_workflow():
    for t in _tasks():
        hit = WORKFLOW_WORDS.search(t["prompt"])
        assert not hit, f"{t['id']} leaked {hit.group(0)!r}"


def test_every_task_is_classed_chemistry_and_confined():
    from delfin.agent.benchmark_runner import workspace_for
    for t in _tasks():
        assert t.get("task_class") == "chemistry", t["id"]
        assert workspace_for(_ROOT, task_class="chemistry") is not None


def test_every_task_has_an_existing_setup_script():
    for t in _tasks():
        setup = t.get("setup", "")
        assert setup.startswith("w20_"), t["id"]
        assert (_PACK / "setup" / setup).is_file(), setup


def test_every_task_carries_grading_levels_and_cost_record():
    for t in _tasks():
        levels = t.get("grading_levels", [])
        assert isinstance(levels, list) and levels, t["id"]
        assert set(levels) <= {"consistency", "reproduction",
                               "reference", "trap"}, t["id"]
        assert "cost" in t, t["id"]


def test_the_count_dict_matches_the_new_task_count():
    """The pinned counts in the 'writes where it belongs' test must grow
    by exactly the chemistry tasks this file adds."""
    import tests.test_every_task_that_writes_writes_where_it_belongs as wc
    n_mine = sum(1 for t in _tasks()
                 if t.get("task_class") == "chemistry")
    assert n_mine > 0
    # The dict is a module constant in that test; read it through the
    # test's own pinned expectation rather than duplicating it here.
    src = pathlib.Path(wc.__file__).read_text(encoding="utf-8")
    # The pinned dict is the LAST "chemistry": N occurrence; earlier
    # ones belong to other tests in the same file.
    matches = re.findall(r'"chemistry": (\d+)', src)
    assert matches, "chemistry count not pinned in the writes test"
    assert int(matches[-1]) == 5 + n_mine
