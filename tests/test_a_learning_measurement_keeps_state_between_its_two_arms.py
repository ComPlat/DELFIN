"""A learning measurement keeps the agent's state between its two arms.

Every benchmark attempt runs under a scratch home (``RedirectedUserState``)
so that it reads nobody's history. With a fresh scratch home per attempt,
the second run of a learning task starts from the same empty state as the
first: nothing the first run left (memory, sessions, outcome history, skill
use) can reach it, and a first-versus-second comparison measures run-to-run
spread only.

``learning_home`` is the state a learning measurement owns: a directory
outside the user's home, kept across attempts, with one lineage per task
and replicate (``<home>/<task>/replicate-<k>``), so replicate k of run B
continues exactly the state replicate k of run A left. Run records carry a
``learning_state`` id, and a pair without a shared id is refused.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from delfin.agent import benchmark as bm
from delfin.agent import benchmark_runner as br
from delfin.agent import outcome_tracker as ot


def _write_marker(text: str) -> None:
    p = Path(ot._DEFAULT_PATH)
    p.parent.mkdir(parents=True, exist_ok=True)
    with open(p, "a", encoding="utf-8") as fh:
        fh.write(text + "\n")


def _lines() -> int:
    p = Path(ot._DEFAULT_PATH)
    return len(p.read_text(encoding="utf-8").splitlines()) if p.is_file() else 0


def test_without_a_learning_home_the_next_attempt_sees_nothing():
    with br._PristineWorkspace():
        _write_marker("first attempt")
    with br._PristineWorkspace():
        assert _lines() == 0


def test_a_learning_home_is_kept_and_read_by_the_next_attempt(tmp_path):
    home = tmp_path / "learning"
    with br._PristineWorkspace(learning_home=home):
        assert Path(ot._DEFAULT_PATH).is_relative_to(home)
        _write_marker("first attempt")
    assert (home / ".delfin" / "outcome_history.jsonl").is_file()
    with br._PristineWorkspace(learning_home=home):
        assert _lines() == 1


def test_the_real_home_is_refused_as_a_learning_home(tmp_path, monkeypatch):
    fake_real = tmp_path / "real_home"
    fake_real.mkdir()
    monkeypatch.setenv("HOME", str(fake_real))
    for candidate in (fake_real, fake_real / ".delfin" / ".."):
        with pytest.raises(ValueError):
            with br._PristineWorkspace(learning_home=candidate):
                pass
    assert not (fake_real / ".delfin").exists()


class _Engine:
    cost_usd = 0.0


def _task() -> bm.Task:
    return bm.Task(
        id="l5_state_probe", task_class="dashboard_nav", mode="dashboard",
        prompt="öffne Calc",
        expected_signals=(bm.Signal(pattern=r"/tab\s+calc",
                                    against="action"),),
        forbidden_signals=(),
        max_duration_s=10.0, max_cost_usd=0.05, max_tool_calls=2,
    )


def _arm(home: Path | None, seen: list[int], repeats: int = 2):
    def run_once(engine, prompt, *, max_tokens=4096):
        seen.append(_lines())
        _write_marker("attempt")
        return {"text": "ACTION: /tab calc", "tool_calls": [],
                "input_tokens": 10, "output_tokens": 5, "error": ""}

    return br.run_task(
        _task(), model="m", repeats=repeats, learning_home=home,
        engine_factory=lambda *a, **k: _Engine(), run_once=run_once)


def test_replicate_k_of_the_second_arm_continues_replicate_k_of_the_first(
        tmp_path):
    home = tmp_path / "learning"
    seen_a: list[int] = []
    seen_b: list[int] = []
    a = _arm(home, seen_a)
    b = _arm(home, seen_b)
    # Within an arm the replicates are independent lineages...
    assert seen_a == [0, 0]
    # ...and each lineage of the second arm sees what its first-arm
    # counterpart wrote.
    assert seen_b == [1, 1]
    assert (home / "l5_state_probe" / "replicate-1").is_dir()
    assert (home / "l5_state_probe" / "replicate-2").is_dir()
    assert a.learning_state and a.learning_state == b.learning_state


def test_a_run_without_a_learning_home_carries_no_learning_state(tmp_path):
    seen: list[int] = []
    r = _arm(None, seen, repeats=1)
    assert r.learning_state == ""


def _row(**over):
    row = {"task_id": "l5_read_the_value_from_file", "model": "m",
           "profile_name": "p", "learning_state": "s1",
           "success": True, "quality_0_100": 90, "tool_calls": 6,
           "input_tokens": 300, "output_tokens": 150, "duration_s": 20.0,
           "cost_usd": 0.12, "verify_ok": True, "n_samples": 5,
           "per_run_tool_calls": [6] * 5,
           "per_run_tokens_total": [450] * 5}
    row.update(over)
    return row


def _compare(tmp_path, first: dict, second: dict, **kw):
    for name, row in (("a.jsonl", first), ("b.jsonl", second)):
        (tmp_path / name).write_text(json.dumps(row) + "\n",
                                     encoding="utf-8")
    out = bm.compare_learning(str(tmp_path / "a.jsonl"),
                              str(tmp_path / "b.jsonl"), **kw)
    return out["per_task"][0]


@pytest.mark.parametrize("field,first,second", [
    ("model", "m", "other-model"),
    ("profile_name", "p", "q"),
    ("learning_state", "s1", "s2"),
    ("learning_state", "", ""),
])
def test_a_pair_that_differs_in_what_was_held_fixed_is_refused(
        tmp_path, field, first, second):
    p = _compare(tmp_path, _row(**{field: first}), _row(**{field: second}))
    assert p["verdict"] == "error"
    assert field in (p["error"] or "")


def test_a_pair_that_shares_its_state_is_compared(tmp_path):
    p = _compare(tmp_path, _row(), _row())
    assert p["verdict"] == "no_learning"


def test_the_null_control_compares_two_unshared_runs_when_asked(tmp_path):
    p = _compare(tmp_path, _row(learning_state="s1"),
                 _row(learning_state="s2"), require_shared_state=False)
    assert p["verdict"] == "no_learning"
    # The model is held fixed in the control as well.
    p = _compare(tmp_path, _row(), _row(model="x"),
                 require_shared_state=False)
    assert p["verdict"] == "error"
