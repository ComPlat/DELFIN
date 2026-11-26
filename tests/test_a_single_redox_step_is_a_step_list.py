"""A single redox step is a step list, as the guide and the paper calculations write it.

`docs/USER_MANUAL.md` writes `oxidation_steps=1`, and four of the calculations
shipped under `examples/DELFIN_Paper_Calcualtions/` do the same. The CONTROL
parser splits a value only on a comma, so `oxidation_steps=1` reaches the
validator as the bare int 1 -- and `_as_list` rejected it with
"must be a list or comma-separated string". Because `run_OCCUPIER` and the
`OCCUPIER_parser` call sites catch that error, the run did not crash: it logged
and returned, planning no redox stage at all. A valid CONTROL silently produced
nothing.

The test drives `OCCUPIER_parser`, the entry point those code paths use (the
non-literal parse), over the shipped paper CONTROLs themselves; `read_control_file`
is the control case that already worked on the base.
"""

import pytest

from delfin.config import OCCUPIER_parser, read_control_file
from delfin.workflows.engine.occupier import _parse_step_list

PAPER = "examples/DELFIN_Paper_Calcualtions"
SINGLE_STEP = ["FeTPP_classic", "FeTPP_OCCUPIER", "HS_FeTpy_Me", "LS_FeTpy"]


@pytest.mark.parametrize("name", SINGLE_STEP)
def test_a_single_redox_step_survives_the_occupier_parse(name, monkeypatch):
    """The shipped paper calculations write oxidation_steps=1; it must parse.

    Red on the base: `OCCUPIER_parser` raised
    `ValueError: Invalid value for oxidation_steps: must be a list or
    comma-separated string` for every one of these four folders.
    """
    monkeypatch.chdir(_repo_root())
    config = OCCUPIER_parser(f"{PAPER}/{name}/CONTROL.txt")

    # One value means the one-element list -- the same meaning as "1," would
    # have, and what the workflow builders below consume.
    assert _parse_step_list(config["oxidation_steps"]) == [1]
    assert _parse_step_list(config["reduction_steps"]) == [1, 2]


@pytest.mark.parametrize("name", SINGLE_STEP)
def test_the_literal_parse_still_reads_the_same_step(name):
    """Control case, green on the base as well: `read_control_file` kept the
    text `'1'` and never hit the rejecting branch. It must be unchanged."""
    config = read_control_file(f"{_repo_root()}/{PAPER}/{name}/CONTROL.txt")
    assert _parse_step_list(config["oxidation_steps"]) == [1]


def test_a_two_step_list_is_unchanged(monkeypatch):
    """Control case, green on the base: a comma list already worked."""
    monkeypatch.chdir(_repo_root())
    config = OCCUPIER_parser(f"{PAPER}/ZnTPP/CONTROL.txt")
    assert _parse_step_list(config["oxidation_steps"]) == [1, 2]


def test_a_single_reduction_step_is_accepted(monkeypatch, tmp_path):
    """The same parser reads `reduction_steps`; a single value must work there too.

    `reduction_steps` goes through the same `_as_list` FieldSpec as
    `oxidation_steps`, so the one-value form has to be accepted on both keys --
    fixing only one would leave half the finding open.
    """
    monkeypatch.chdir(_repo_root())
    control = tmp_path / "CONTROL.txt"
    control.write_text(
        _paper_control("FeTPP_classic").replace(
            "reduction_steps=1,2", "reduction_steps=1"
        ),
        encoding="utf-8",
    )
    config = OCCUPIER_parser(str(control))
    assert _parse_step_list(config["reduction_steps"]) == [1]


def test_a_value_that_is_not_a_step_list_is_still_refused(monkeypatch, tmp_path):
    """The negative case: the guard must still refuse what is not a list.

    Without this, accepting the int would be indistinguishable from accepting
    everything. A dict is not a step list and must keep raising.
    """
    from delfin.common.control_validator import _as_list

    with pytest.raises(ValueError, match="must be a list or comma-separated string"):
        _as_list({"1": 2})


def _repo_root():
    from pathlib import Path

    return Path(__file__).resolve().parents[1]


def _paper_control(name):
    return (_repo_root() / PAPER / name / "CONTROL.txt").read_text(encoding="utf-8")
