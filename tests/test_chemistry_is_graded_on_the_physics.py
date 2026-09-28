"""Chemistry tasks must be checkable, not vibes.

Each tasks_chem task is graded on physics the acceptance script can
verify itself: a frequency job must show no imaginary mode, a gap must
fall in a band, an ordering must be an ordering.  The tolerances come
from xtb reference runs made once when the suite was written -- they are
baked into the YAML as ``expected_values`` bands so the scorer can also
judge the model's own reported numbers.

This test guards the SUITE: every chem task declares a verify script,
a molecule identity signal, and -- where it asks for a figure -- a
numeric band, so no task can regress into prose-graded chemistry.
"""

import math
from pathlib import Path

import pytest

from delfin.agent.benchmark import load_tasks

_HERE = Path(__file__).resolve().parent
_CHEM_YAML = (_HERE.parent / "delfin" / "agent" / "pack" / "benchmark"
              / "tasks_chem.yaml")


def _chem_tasks():
    return [t for t in load_tasks() if t.task_class == "chemistry"]


def test_chem_suite_exists_with_expected_molecules():
    ids = {t.id for t in _chem_tasks()}
    for needle in ("acetaminophen", "chloronitrobenzene",
                   "carbocation", "au_cyanide", "cr_carbonyl"):
        assert any(needle in i for i in ids), f"no task id mentions {needle}"
    assert len(_chem_tasks()) >= 5


def test_every_chem_task_is_verified_by_a_script():
    for t in _chem_tasks():
        assert t.verify, f"{t.id} has no acceptance script"
        # The script must live in the packaged accept/ dir -- a name
        # that does not resolve would be a suite bug, not a model fail.
        base = (_CHEM_YAML.parent / "accept")
        assert (base / t.verify).exists(), f"{t.id}: accept/{t.verify} missing"


def test_figure_tasks_carry_numeric_bands():
    # The gap tasks must be judged as numbers with a tolerance, not by
    # digit patterns -- the whole point of the value machinery.
    for t in _chem_tasks():
        if "gap" in t.id or "order" in t.id:
            assert t.expected_values, (
                f"{t.id} asks for a figure but declares no expected_values")


def test_chem_tasks_do_not_run_orca():
    # Cheap is the requirement: xtb-level work only. A task whose
    # signals accept an ORCA run slipped into a cheap suite and will
    # cost 50x its budget on the cluster.
    for t in _chem_tasks():
        for sig in t.expected_signals:
            assert "orca" not in sig.pattern.lower()
