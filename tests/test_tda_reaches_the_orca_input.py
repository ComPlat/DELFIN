"""A decision recorded and not executed is worse than one never asked for.

The excited-state ladder asks an agent to say whether the Tamm-Dancoff
approximation is used, and names ``tda`` in the step's kwargs as what
closes that decision. Measured 2026-11-23: nothing read the kwarg.
``OrcaTddftAdapter.execute`` injected only ``TDDFT NRoots N`` into the
bang line, ``_run`` never looked at ``tda``, and the base
``validate_params`` is a no-op -- so a spec that said ``tda=false``
validated, linted clean, graded as justified, and ORCA ran its default
(TDA on) anyway. The published number then carried a provenance that was
false, which is the one failure mode a provenance record exists to
prevent.

What is asserted here is the input FILE, not a return value: the input is
written before the binary is called, so the claim holds whether or not
ORCA is installed on the machine running the test.

There is no dependable bang-line spelling for TDA, so it goes in the
``%tddft`` block.
"""
from __future__ import annotations

import pytest

from delfin.tools.adapters.orca import OrcaTddftAdapter


@pytest.fixture()
def run(tmp_path):
    """Drive the real adapter and hand back the input it wrote."""
    geo = tmp_path / "h2.xyz"
    geo.write_text("2\nh2\nH 0.0 0.0 0.0\nH 0.0 0.0 0.74\n", encoding="utf-8")
    counter = {"n": 0}

    def _run(**kwargs) -> str:
        counter["n"] += 1
        work = tmp_path / f"job{counter['n']}"
        work.mkdir()
        try:
            OrcaTddftAdapter().execute(
                work, geometry=geo, cores=1, charge=0, mult=1,
                method="PBE0", basis="def2-SVP", base_name="calc", **kwargs)
        except Exception:
            # The ORCA binary may be absent or may refuse the job; the
            # input file is written before it is called, and the input is
            # what this module is about.
            pass
        inp = work / "calc.inp"
        assert inp.is_file(), "no input file was written at all"
        return inp.read_text(encoding="utf-8")

    return _run


def _tddft_lines(text: str) -> list[str]:
    return [ln.strip() for ln in text.splitlines()
            if "tda" in ln.lower() or "%tddft" in ln.lower()]


# --- the defect ---------------------------------------------------------

@pytest.mark.parametrize("value,expected", [(False, "false"), (True, "true")])
def test_the_decision_reaches_the_input(run, value, expected):
    """Either value, written where ORCA reads it."""
    text = run(tda=value)
    assert "%tddft" in text.lower()
    assert f"tda {expected}" in text.lower(), _tddft_lines(text)


def test_not_saying_leaves_the_engine_default_alone(run):
    """The rung asks for a decision; it does not impose one. A step that
    says nothing must produce the input it produced before this change,
    so the ladder cannot silently alter runs that never opted in."""
    text = run()
    assert not any("tda" in ln.lower() for ln in _tddft_lines(text)), \
        _tddft_lines(text)


def test_a_block_the_caller_wrote_is_not_touched(run):
    """Two `%tddft` blocks is an ORCA error, and the caller who spelled
    the block out has already decided. Theirs wins."""
    text = run(tda=False, extra_blocks="%tddft\n  NRoots 3\nend")
    assert text.lower().count("%tddft") == 1
    assert "nroots 3" in text.lower()


def test_the_block_is_added_beside_other_blocks(run):
    """An unrelated block must survive: the TDA block is appended, not
    substituted for whatever the caller passed."""
    text = run(tda=True, extra_blocks="%geom TolE 1e-5 end")
    assert "%geom" in text.lower()
    assert "tda true" in text.lower()


# --- the contract the ladder depends on --------------------------------

def test_tda_is_a_declared_parameter_of_the_step():
    """The ladder tells the model to set `tda` in the step's kwargs. A
    parameter that is not declared is not part of the step's contract,
    so `validate_spec` could not report a typo in it and the resolver
    could not carry a default."""
    names = {p.name for p in OrcaTddftAdapter.params}
    assert "tda" in names, sorted(names)
    spec = next(p for p in OrcaTddftAdapter.params if p.name == "tda")
    assert spec.type == "bool"
    assert spec.description and "default" in spec.description.lower()


def test_nroots_still_reaches_the_bang_line(run):
    """What must not move with the change: the root count was already
    wired through the bang line and stays there."""
    text = run(nroots=7)
    assert "nroots 7" in text.lower()
