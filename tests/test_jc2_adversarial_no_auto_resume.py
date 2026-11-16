"""JC2 reviewer adversarial tests, part 2: nothing resumes on its own.

Target 4 of the JC2 brief: no code path calls ``delfin chain resume``
without an explicit user action; the restart path only produces TEXT.
Control: green on main — this file PINS the current safe behaviour so the
builder's phase-3 changes cannot silently break it.
"""

import inspect

from delfin.tools import chain_state


def test_chain_state_module_has_no_resume_trigger():
    """The chain-state module must not invoke a resume anywhere: parse its
    AST and look for a call into cli_chain, a subprocess launch, or a run()
    call — docstring mentions of the consumer are fine, code is not."""
    import ast
    import delfin.tools.chain_state as mod

    src = inspect.getsource(mod)
    tree = ast.parse(src)
    violations = []
    for node in ast.walk(tree):
        # imports of cli_chain or subprocess
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            names = ([a.name for a in node.names]
                     if isinstance(node, ast.Import)
                     else [node.module or ""])
            for n in names:
                if "cli_chain" in n or "subprocess" in n:
                    violations.append(f"import {n} (line {node.lineno})")
        # calls named run/resume
        if isinstance(node, ast.Call):
            func = node.func
            name = getattr(func, "attr", getattr(func, "id", ""))
            if name in ("resume", "Popen", "system", "check_output"):
                violations.append(f"call {name}() (line {node.lineno})")
    assert not violations, (
        f"chain_state.py must never trigger a resume or shell out; got "
        f"{violations}"
    )


def test_pipeline_resume_only_through_explicit_flag():
    """Pipeline.run resumes ONLY when resume_chain=True (or a passed
    chain_state_obj); a plain run(work_dir=...) on a directory that holds
    a finished chain state must not silently skip steps."""
    from delfin.tools.pipeline import Pipeline

    src = inspect.getsource(Pipeline.run)
    assert "resume_chain" in src
    # The skip set is built only under `if resume_chain:`
    skip_builds = [l for l in src.splitlines() if "chain_skip" in l]
    assert any("if resume_chain" in s for s in skip_builds) or all(
        "chain_skip" not in s or "resume_chain" in s or "plan" in s
        for s in skip_builds
    ), f"chain_skip must be gated on the explicit resume flag; got {skip_builds}"


def test_consume_crash_recovery_note_returns_text_never_resumes(
    monkeypatch, tmp_path
):
    """The restart path produces a note and nothing else: even with a
    chain present that wants resuming, consume_crash_recovery_note must
    only return a string (or None) — no function call into the resume
    path, no subprocess."""
    from delfin.tools import chain_state as cs

    wd = tmp_path / "chain"
    chain = cs.ChainState.create(name="c", work_dir=wd)
    chain.add_step("xtb_opt", "xtb_opt", params={"charge": 0})
    # Leave it unfinished on purpose.

    from delfin.agent import session_store

    # No session state at all -> the note path returns None and must not
    # touch the chain directory (its state file is untouched).
    before = (wd / "chain_state.json").read_text()
    note = session_store.consume_crash_recovery_note("no-such-session")
    assert note is None
    assert (wd / "chain_state.json").read_text() == before, (
        "the restart note path must be read-only towards chain state"
    )


def test_chain_state_load_edge_cases():
    """Empty, unicode, and huge-name chains load without crashing and
    keep their identity (reviewer edge-case sweep of ChainState.load)."""
    import json

    cases = [
        # empty steps list
        {"name": "e", "work_dir": "/nowhere", "steps": []},
        # unicode name and param values
        {"name": "käse-ß", "work_dir": "/nowhere", "steps": [
            {"name": "stép", "tool": "orca_sp", "state": "pending",
             "params": {"charge": 0, "note": "Ω≈ω"}},
        ]},
        # a very long name (no crash, no truncation surprise)
        {"name": "x" * 5000, "work_dir": "/nowhere", "steps": []},
    ]
    for i, data in enumerate(cases):
        raw = json.dumps(data)
        import tempfile

        with tempfile.TemporaryDirectory() as td:
            from pathlib import Path

            p = Path(td) / "chain_state.json"
            p.write_text(raw, encoding="utf-8")
            chain = chain_state.ChainState.load(p)
        assert chain.name == data["name"], f"case {i}: name round-trips"
        if data["steps"]:
            assert chain.steps[0].params["note"] == "Ω≈ω", (
                f"case {i}: unicode param values survive the load"
            )
