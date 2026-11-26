"""A capability tag is not part of the question a task asks.

``bench ab`` lets each SIDE hash its own task, so that a base predating
``bench_ab.py`` can still describe the question it asks. That only works while
the digest covers fields both sides define: a field one side carries and the
other cannot is read as a changed question.

``Task.capability`` is such a field, and it is an annotation -- which product
capability the task is meant to measure. It changes nothing the model is asked
or judged on. Hashed whole, it made every task "a DIFFERENT question" between a
tagged tree and an untagged one, so the verb refused every comparison against a
historical base and blamed the task's text.

The assertions here are invariant / control / agreement, in that order.
"""

from __future__ import annotations

import dataclasses

from delfin.agent.benchmark import Task, load_tasks


def _ab():
    """The module under test, imported on use (it does not exist on main)."""
    from delfin.agent import bench_ab

    return bench_ab


def _task(**kw) -> Task:
    fields = dict(id="t-alpha", task_class="generalist", mode="solo",
                  prompt="do the thing")
    fields.update(kw)
    return Task(**fields)


def test_a_capability_tag_does_not_change_the_task_digest():
    """The defect: tagging a task made the two sides disagree about it.

    Red before the fix: the digest serialised every field, so adding a tag
    produced a different hash for an unchanged question, and ``_check_pair``
    refused the pair with "the same task id is a DIFFERENT question".
    """
    ab = _ab()
    plain = _task()
    tagged = dataclasses.replace(plain, capability=("tool:grep_file",))

    assert ab.task_digest(tagged) == ab.task_digest(plain), (
        "a capability tag changed the task digest, so a comparison whose base "
        "predates the tag is refused as 'a DIFFERENT question' although the "
        "task asks exactly the same thing"
    )


def test_a_changed_prompt_is_still_a_different_question():
    """The control: the guarantee above must not be bought by hashing nothing.

    Without this, dropping every field from the payload would satisfy the test
    above and the refusal it guards would never fire again.
    """
    ab = _ab()
    assert ab.task_digest(_task(prompt="do the thing")) != \
        ab.task_digest(_task(prompt="do something else"))


def test_a_reworded_task_in_the_shipped_pack_keeps_its_digest_under_tagging():
    """And it holds for a real task, not only for a constructed one."""
    ab = _ab()
    tagged = next(t for t in load_tasks() if getattr(t, "capability", ()))
    stripped = dataclasses.replace(tagged, capability=())
    assert ab.task_digest(stripped) == ab.task_digest(tagged)


def test_the_worker_copies_the_parents_digest_rule():
    """Two copies of the digest exist on purpose; they must agree.

    The worker is source text, not an import: the side it runs on may have no
    ``bench_ab.py`` at all. A rule that drifted between the copies would report
    a difference for every task, so the payload both of them build is the one
    thing that has to stay in step. Run the worker verbatim and compare.
    """
    import json
    import subprocess
    import sys
    import tempfile
    import textwrap
    from pathlib import Path

    ab = _ab()
    task = next(t for t in load_tasks() if getattr(t, "capability", ()))
    prelude = textwrap.dedent("""
        import sys, types
        _stub = types.ModuleType("delfin.agent.benchmark_runner")

        class _R:
            pass

        def run_task(task, **kw):
            return _R()

        _stub.run_task = run_task
        sys.modules["delfin.agent.benchmark_runner"] = _stub
    """)
    source = prelude + "\n" + ab._WORKER_SOURCE
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "attempt.json"
        proc = subprocess.run(
            [sys.executable, "-P", "-c", source,
             str(task.id), "m", "", "api", "1", str(out)],
            cwd=tmp, capture_output=True, text=True, timeout=120,
            env={"PATH": "/usr/bin:/bin",
                 "PYTHONPATH": str(Path(ab.__file__).resolve().parents[2])})
        assert out.exists(), (
            "the worker wrote no result file; stderr was:\n"
            + (proc.stderr or proc.stdout or "")[-1500:])
        payload = json.loads(out.read_text(encoding="utf-8"))

    assert payload["task_digest"] == ab.task_digest(task), (
        "the worker and the parent disagree about the same task, so every "
        "side would report a changed question for an unchanged task")
