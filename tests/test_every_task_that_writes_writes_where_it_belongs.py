"""Only office tasks were confined; every other writer had the checkout.

``_default_engine_factory`` narrowed the workspace to
``tests/fixtures/office_workspace`` when the MODE was ``office``, and
handed everything else ``os.getcwd()`` — the whole repository, writable.

The mode cannot do better: behaviour tasks and generic-project tasks both
run as ``solo``. So the classes that edit toy files just as concretely as
office does were pointed at the checkout, and the only thing keeping a
run from leaving edits in ``delfin/`` was that the tasks happened to
behave. #110 added DETECTION for that; this is the prevention it named.

Counted over the four packaged task files:

    office            11   already confined
    behavior_*        12   were not
    generic_project    8   were not
    everything else   48   reads only — deliberately left alone, see below

The 48 keep the repository as their root. Confining them needs more than
a different directory: 12 of them name real ``delfin/`` paths and read
them, and ``repo_dir`` is what relative paths resolve against, so moving
it turns a legitimate read into "file not found". A previous attempt did
exactly that and broke fixture discovery. Reads need a read-only reach
into the checkout, which is a separate mechanism (``read_only_workspace_
dirs``) and a separate change.
"""

from __future__ import annotations

import pathlib

import pytest
import yaml

from delfin.agent.benchmark_runner import workspace_for

_ROOT = pathlib.Path(__file__).resolve().parents[1]
_PACK = _ROOT / "delfin" / "agent" / "pack" / "benchmark"


def _tasks():
    out = []
    for f in sorted(_PACK.glob("*.yaml")):
        data = yaml.safe_load(f.read_text(encoding="utf-8"))
        out.extend(data if isinstance(data, list) else (data.get("tasks") or []))
    return out


# ---------------------------------------------------------------------------
# The classes that write get their own directory
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("task_class,expected", [
    ("office", "office_workspace"),
    ("behavior_ask", "behavior_workspace"),
    ("behavior_plan", "behavior_workspace"),
    ("behavior_scout", "behavior_workspace"),
    ("behavior_verify", "behavior_workspace"),
    ("generic_project", "user_project_workspace"),
])
def test_a_writing_class_is_confined_to_its_fixture(task_class, expected):
    ws = workspace_for(_ROOT, task_class=task_class)
    assert ws is not None, f"{task_class} was handed the whole checkout"
    assert ws == _ROOT / "tests" / "fixtures" / expected


def test_the_office_mode_still_decides_on_its_own():
    """Office was keyed on the mode before the class existed as a key.
    Both routes have to keep working — the runner passes both."""
    assert workspace_for(_ROOT, mode="office") == (
        _ROOT / "tests" / "fixtures" / "office_workspace")


# ---------------------------------------------------------------------------
# The classes that only read are left where they were
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("task_class", [
    "dashboard_nav", "solo_research", "fact_verify", "plan_refactor",
    "verify_enforcement", "workflow_pattern", "debug", "fact_verify_auto",
])
def test_a_reading_class_keeps_the_repository(task_class):
    assert workspace_for(_ROOT, task_class=task_class) is None


def test_an_unknown_class_is_not_redirected_somewhere():
    assert workspace_for(_ROOT, task_class="something_new") is None


def test_an_empty_class_is_not_redirected_either():
    assert workspace_for(_ROOT) is None


# ---------------------------------------------------------------------------
# A missing fixture never silently redirects a run
# ---------------------------------------------------------------------------

def test_a_missing_fixture_directory_answers_none(tmp_path):
    """Better the old behaviour than a workspace that is not there: a
    path that does not exist would send every write somewhere new."""
    assert workspace_for(tmp_path, task_class="behavior_ask") is None


def test_a_present_fixture_directory_answers_it(tmp_path):
    (tmp_path / "tests" / "fixtures" / "behavior_workspace").mkdir(parents=True)
    assert workspace_for(tmp_path, task_class="behavior_ask") == (
        tmp_path / "tests" / "fixtures" / "behavior_workspace")


# ---------------------------------------------------------------------------
# Against the packaged tasks, not against my idea of them
# ---------------------------------------------------------------------------

def test_every_packaged_writing_class_has_a_workspace():
    missing = sorted({
        t.get("task_class", "") for t in _tasks()
        if (str(t.get("task_class", "")).startswith("behavior")
            or t.get("task_class") in ("office", "generic_project"))
        and workspace_for(_ROOT, mode=t.get("mode", ""),
                          task_class=t.get("task_class", "")) is None
    })
    assert not missing, f"writing classes with no workspace: {missing}"


def test_the_counts_still_match_what_was_measured():
    """If a task file grows a new writing class, this says so rather than
    letting it inherit the checkout unnoticed."""
    counts = {"office": 0, "behavior": 0, "generic_project": 0, "other": 0}
    for t in _tasks():
        c = str(t.get("task_class", ""))
        if c == "office":
            counts["office"] += 1
        elif c.startswith("behavior"):
            counts["behavior"] += 1
        elif c == "generic_project":
            counts["generic_project"] += 1
        else:
            counts["other"] += 1
    assert counts == {"office": 11, "behavior": 12,
                      "generic_project": 8, "other": 48}, counts


# ---------------------------------------------------------------------------
# The runner asks for it
# ---------------------------------------------------------------------------

def test_the_runner_hands_the_task_class_to_the_factory():
    """The decision is worthless if the call site does not pass the key."""
    from delfin.agent import benchmark_runner as br
    from delfin.agent.benchmark import Task

    seen: dict = {}

    def _factory(model, backend, provider, mode, task_class=""):
        seen["mode"], seen["class"] = mode, task_class
        raise RuntimeError("stop here — the arguments are the point")

    br.run_task(
        Task(id="t", task_class="behavior_ask", mode="solo", prompt="p"),
        model="m", engine_factory=_factory,
        run_once=lambda *a, **k: {"text": ""},
    )
    assert seen == {"mode": "solo", "class": "behavior_ask"}
