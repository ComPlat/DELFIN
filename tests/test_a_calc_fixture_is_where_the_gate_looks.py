"""A calc fixture only scopes the calc tools if it is where the gate looks.

``benchmark_runner._fixture_calc_dirs`` points ``search_calcs`` and
``calc_summary`` at a task's fixture ONLY when the setup built
``<workspace>/calc_archive/{calc,archive}``. That condition is what keeps a
benchmark off the operator's own corpus -- the docstring beside it says the
real folders hold hundreds of calculations "whose contents nobody chose".

The m2 search_calcs task built ``sv1/calc`` instead, so the gate was a no-op:
the task ran against the real corpus, and the verifier demanded energies that
corpus does not hold. Both halves are wrong, but the FIRST one is the
dangerous half -- a benchmark that silently reads a person's data and reports a
score for it.

This file pins the layout, not the tool: the fixture must build the tree the
gate recognises, and the task prompt must name it.
"""

from __future__ import annotations

from pathlib import Path


def _setup_script() -> Path:
    from delfin.agent.benchmark_runner import setup_path

    return setup_path("m2_capability_fixtures.py")


def _task(task_id: str = "m2_a6_search_calcs_within_method"):
    from delfin.agent.benchmark import load_tasks

    return next(t for t in load_tasks() if t.id == task_id)


def _build(tmp_path: Path) -> Path:
    """Run the pack's own setup script into a workspace."""
    import subprocess
    import sys

    ws = tmp_path / "ws"
    ws.mkdir()
    proc = subprocess.run([sys.executable, str(_setup_script()), str(ws)],
                          capture_output=True, text=True, timeout=300,
                          cwd=str(ws.parent))
    assert proc.returncode == 0, (
        f"the setup script failed:\n{proc.stdout}\n{proc.stderr}")
    return ws


def test_the_calc_fixture_lands_where_the_scoping_gate_looks(tmp_path):
    """The defect: sv1/calc was invisible to _fixture_calc_dirs."""
    ws = _build(tmp_path)
    assert (ws / "calc_archive" / "calc").is_dir(), (
        "the fixture built no calc_archive/calc, so _fixture_calc_dirs yields "
        "without redirecting and the calc tools read the operator's real "
        "corpus: " + ", ".join(sorted(p.name for p in ws.iterdir())))


def test_the_redirected_tools_then_see_the_fixture(tmp_path):
    """The gate, driven for real: enter it and ask what the tools would read.

    This is the property the layout exists for, and it is stated on the gate
    rather than on the fixture, so a future change to either side is caught.
    """
    from delfin.agent.api_client import _doc_executor
    from delfin.agent.benchmark_runner import _fixture_calc_dirs

    ws = _build(tmp_path)
    with _fixture_calc_dirs(ws):
        dirs = getattr(_doc_executor, "_calc_dirs", {})
    assert dirs.get("calc"), (
        "the gate did not redirect the calc tools at this workspace's "
        "fixture, so the task measures the operator's corpus")


def test_the_answer_folder_is_a_sibling_of_the_fixture_not_the_corpus(tmp_path):
    """The answer must not be written INTO the indexed tree.

    A task that writes its answer inside ``calc_archive/calc`` would change the
    corpus it is being scored against, and the next task would see it.
    """
    import re

    ws = _build(tmp_path)
    prompt = _task().prompt
    written = set(re.findall(r"[\w./-]*answer\.md", prompt))
    assert written, "the prompt names no answer file: " + prompt
    for path in written:
        assert "/calc/" not in "/" + path, (
            "the task writes its answer inside the indexed calc tree: " + path)
    assert (ws / "calc_archive").is_dir()


def test_the_prompt_names_the_path_the_fixture_built(tmp_path):
    """A prompt naming a directory the setup never built asks about nothing."""
    ws = _build(tmp_path)
    prompt = _task().prompt
    named = [tok for tok in ("calc_archive", "sv1") if tok in prompt]
    assert named, "the prompt names no fixture path: " + prompt
    for token in named:
        assert (ws / token).exists(), (
            f"the prompt says {token!r} but the setup built no such directory; "
            "the agent would report on an empty tree: "
            + ", ".join(sorted(p.name for p in ws.iterdir())))


def test_the_verifier_reads_the_answer_where_the_prompt_put_it():
    """The prompt and the verifier must name the same file."""
    from delfin.agent.benchmark_runner import acceptance_path

    task = _task()
    source = acceptance_path(task.verify).read_text(encoding="utf-8")
    assert "answer.md" in source
    for token in ("calc_archive", "sv1"):
        if token in task.prompt:
            assert token in source, (
                f"the prompt names {token!r} and the verifier does not, so a "
                "correct answer written where the prompt said is scored wrong")
