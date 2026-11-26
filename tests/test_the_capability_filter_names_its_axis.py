"""One word, one axis: ``--capability`` filters on the field it names.

``bench_history.inst_capability`` reads ``task_class`` from a stored instance
line, and its docstring says naming one accessor is what stops the read side
drifting from the writer. The CLI then documented the filter with
``redox_consistency``, which is neither a ``task_class`` in the pack nor a
capability tag -- so the example an operator copies can only ever select
nothing, and the same word ("capability") means two different things in one
verb.

These tests state the axis and keep the help text honest. They do not change
which field is filtered on: ``task_class`` is the stored vocabulary, and the
one ``bench_report`` already scores per tool family.
"""

from __future__ import annotations


def _history():
    from delfin.agent import bench_history as _bh

    return _bh


def test_the_filter_reads_the_field_the_writer_stores():
    """The accessor's answer is the stored ``task_class``, verbatim."""
    history = _history()
    inst = {"task_id": "t", "task_class": "science_analysis",
            "quality": 3, "ran": True}
    assert history.inst_capability(inst) == "science_analysis"


def test_an_instance_without_a_class_answers_empty_not_none():
    """A missing field must not read as the string ``None``."""
    history = _history()
    assert history.inst_capability({}) == ""


def test_a_documented_example_is_a_value_the_pack_carries():
    """The example in the help text must be a task family that exists.

    The filter is on ``task_class``, so an example that is not one sends the
    operator to a command that prints an empty table forever.
    """
    import re
    from pathlib import Path

    from delfin.agent import cli
    from delfin.agent.benchmark import load_tasks

    source = Path(cli.__file__).read_text(encoding="utf-8")
    # Only the --capability help strings: the paragraphs of source that
    # declare that option. Scanning the whole file picked up unrelated
    # examples ("e.g. add a flag") and proved nothing.
    blocks = re.findall(
        r"\"--capability\"[^)]*?help\s*=\s*((?:\s*\"[^\"]*\")+)",
        source, re.S)
    assert blocks, "the --capability help strings were not found"
    examples: set[str] = set()
    for block in blocks:
        for piece in re.findall(r"\"([^\"]*)\"", block):
            examples |= set(re.findall(r"e\.g\. ([a-z][a-z0-9_]*)", piece))
    assert examples, (
        "the --capability help text carries no example at all; a filter whose "
        "axis is not stated is a filter nobody can use")
    classes = {t.task_class for t in load_tasks()}
    unknown = sorted(examples - classes)
    assert not unknown, (
        "the help text names a task family the pack does not carry, so an "
        f"operator who copies it selects nothing: {unknown}")


def test_the_history_filter_selects_only_the_named_class():
    """End to end on the rows, not on the help text."""
    from delfin.agent.bench_history import Round

    history = _history()
    rnd = Round(round_id="r1", instances=[
        {"task_id": "a", "task_class": "coding", "ran": True},
        {"task_id": "b", "task_class": "chemistry", "ran": True}])
    rows = history.history_rows(rnd, capability="chemistry")
    assert [r["task_id"] for r in rows] == ["b"]


def test_no_filter_keeps_every_row():
    """The control: the filter must not drop rows when it is not used."""
    from delfin.agent.bench_history import Round

    history = _history()
    rnd = Round(round_id="r1", instances=[
        {"task_id": "a", "task_class": "coding", "ran": True},
        {"task_id": "b", "task_class": "chemistry", "ran": True}])
    assert len(history.history_rows(rnd)) == 2
