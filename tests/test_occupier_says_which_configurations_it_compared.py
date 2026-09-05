"""OCCUPIER decides which electron configuration a molecule adopts.

Kohn-Sham DFT returns the energy of *a* self-consistent solution, not of the
ground state. For open-shell and transition-metal systems the SCF has several
local minima — different multiplicities, and broken-symmetry solutions where
alpha density localises on one centre and beta on another — and which one is
reached depends on the starting orbitals. OCCUPIER therefore enumerates
candidates, converges each deliberately from a named predecessor, and compares.

That predecessor chain is chemistry, not scheduling: starting a BS candidate
from the converged orbitals of its pure counterpart is what makes the BS
solution findable instead of collapsing back. So when one configuration fails,
everything downstream of it goes with it — and the comparison then names a
preferred configuration out of whatever survived.

These tests pin what the run says about that, and what may enter the
comparison in the first place. Failing to converge is a result, not a
malfunction: not every spin state exists for a given system.
"""

from pathlib import Path

import pytest

from delfin.occupier_flat_extraction import _parse_energy, _fallback_propagate_geometry
from delfin.reporting.occupier_reports import (
    _configuration_outcomes,
    _configuration_survey,
    _thermally_populated,
)
from delfin.occupier_auto import (
    BOLTZMANN_PARENT_THRESHOLD,
    merge_sequences,
    populated_parents,
)


_CONVERGED = """
         *                GEOMETRY OPTIMIZATION CYCLE   4            *
FINAL SINGLE POINT ENERGY     -75.512600000000
                    ***********************HURRAY********************
                    ***        THE OPTIMIZATION HAS CONVERGED     ***
                             ****ORCA TERMINATED NORMALLY****
"""

_RAN_OUT_OF_CYCLES = """
         *                GEOMETRY OPTIMIZATION CYCLE   1            *
FINAL SINGLE POINT ENERGY     -75.400000000000
                    The optimization did not converge but reached the maximum
                    number of optimization cycles.
                             ****ORCA TERMINATED NORMALLY****
"""

_KILLED_MIDWAY = """
         *                GEOMETRY OPTIMIZATION CYCLE   3            *
FINAL SINGLE POINT ENERGY     -75.300000000000
         *                GEOMETRY OPTIMIZATION CYCLE   4            *
FINAL SINGLE POINT ENERGY     -75.310000000000
"""


def _out(tmp_path: Path, name: str, body: str) -> Path:
    p = tmp_path / name
    p.write_text(body, encoding="utf-8")
    return p


# --------------------------------------------------------------------------
# what may enter the comparison
# --------------------------------------------------------------------------

def test_a_converged_optimisation_contributes_its_energy(tmp_path):
    assert _parse_energy(_out(tmp_path, "output.out", _CONVERGED), False) == pytest.approx(-75.5126)


def test_an_optimisation_that_ran_out_of_cycles_does_not_compete(tmp_path):
    """Its last energy belongs to a geometry that is not a stationary point.
    Ranked against converged states it can win on nothing but being unfinished."""
    assert _parse_energy(_out(tmp_path, "output4.out", _RAN_OUT_OF_CYCLES), False) is None


def test_a_run_killed_partway_does_not_compete(tmp_path):
    """The file still holds one FINAL SINGLE POINT ENERGY per completed cycle,
    and the last of them used to be taken as the configuration's energy."""
    assert _parse_energy(_out(tmp_path, "output7.out", _KILLED_MIDWAY), False) is None


def test_gibbs_energies_were_already_safe(tmp_path):
    """A crashed run has no "Final Gibbs free energy" line, which is why this
    only ever bit the FSPE path — the default."""
    assert _parse_energy(_out(tmp_path, "output.out", _KILLED_MIDWAY), True) is None


# --------------------------------------------------------------------------
# what the report says about the tree
# --------------------------------------------------------------------------

_SEQUENCE = [
    {"index": 1, "m": 1, "BS": "", "from": 0},
    {"index": 4, "m": 3, "BS": "", "from": 1},
    {"index": 7, "m": 5, "BS": "", "from": 4},
]


def test_each_planned_configuration_gets_a_verdict(tmp_path):
    _out(tmp_path, "output.out", _CONVERGED)
    _out(tmp_path, "output4.out", _RAN_OUT_OF_CYCLES)
    # no output7.out: it starts from 4 and was never reached

    outcomes = _configuration_outcomes(_SEQUENCE, [-75.5126, None, None], folder=tmp_path)

    assert [state for _, state, _ in outcomes] == ["compared", "not converged", "not reached"]
    assert "FoB 4" in outcomes[2][2], "it has to say which predecessor was missing"


def test_the_survey_states_how_much_of_the_tree_was_compared(tmp_path):
    _out(tmp_path, "output.out", _CONVERGED)
    _out(tmp_path, "output4.out", _RAN_OUT_OF_CYCLES)

    survey = _configuration_survey(
        _configuration_outcomes(_SEQUENCE, [-75.5126, None, None], folder=tmp_path)
    )

    assert "1 of 3 planned configurations" in survey
    assert "not of all of them" in survey
    assert "m=3" in survey and "not converged" in survey
    assert "m=5" in survey and "not reached" in survey


def test_a_complete_tree_says_nothing_extra(tmp_path):
    """The block exists to explain a gap. With no gap it would be noise."""
    outcomes = _configuration_outcomes(_SEQUENCE, [-75.51, -75.40, -75.30], folder=tmp_path)

    assert _configuration_survey(outcomes) == ""


def test_a_configuration_that_failed_on_its_own_is_not_blamed_on_a_predecessor(tmp_path):
    """FoB 1 starts from the input geometry; there is nothing upstream to blame."""
    _out(tmp_path, "output.out", _RAN_OUT_OF_CYCLES)

    outcomes = _configuration_outcomes(_SEQUENCE[:1], [None], folder=tmp_path)

    assert outcomes[0][1] == "not converged"


# --------------------------------------------------------------------------
# a geometry that did not come from OCCUPIER
# --------------------------------------------------------------------------

def test_a_fallback_geometry_is_marked_as_such(tmp_path):
    """Without a summary the run continues from the geometry OCCUPIER was
    handed: not optimised, spin state unresolved. Downstream steps read the
    .xyz and cannot tell, so it is said next to the file."""
    stage = tmp_path / "initial_OCCUPIER"
    stage.mkdir()
    (stage / "input.xyz").write_text("1\n\nH 0.0 0.0 0.0\n", encoding="utf-8")

    _fallback_propagate_geometry("initial_OCCUPIER", stage)

    assert (tmp_path / "initial.xyz").exists(), "the run has to keep going"
    marker = tmp_path / "initial.xyz.not-from-occupier"
    assert marker.exists()
    text = marker.read_text(encoding="utf-8")
    assert "did not come from OCCUPIER" in text
    assert "not optimised" in text


# --------------------------------------------------------------------------
# thermal mixtures seed the next redox step
# --------------------------------------------------------------------------

def test_one_dominant_state_yields_one_parent():
    populated = _thermally_populated(_SEQUENCE, {1: 0.97, 4: 0.03})

    assert [p["m"] for p in populated] == [1]


def test_two_populated_states_both_become_parents():
    """A 50:50 mixture is two species. Oxidising or reducing it can start from
    either, and the daughters reachable from one are not those from the other."""
    populated = _thermally_populated(_SEQUENCE, {1: 0.50, 4: 0.50})

    assert [p["m"] for p in populated] == [1, 3]


def test_the_threshold_is_where_it_is_documented_to_be():
    just_under = BOLTZMANN_PARENT_THRESHOLD - 0.001
    assert [p["m"] for p in _thermally_populated(_SEQUENCE, {1: 1 - just_under, 4: just_under})] == [1]
    just_over = BOLTZMANN_PARENT_THRESHOLD + 0.001
    assert len(_thermally_populated(_SEQUENCE, {1: 1 - just_over, 4: just_over})) == 2


def test_a_state_file_without_populations_still_works():
    """Runs recorded before this existed carry only the winner."""
    assert [p["m"] for p in populated_parents({"index": 1, "m": 1, "BS": ""})] == [1]


def test_populations_are_ordered_most_populated_first():
    populated = _thermally_populated(_SEQUENCE, {1: 0.30, 4: 0.70})

    assert [p["m"] for p in populated] == [3, 1]


# --------------------------------------------------------------------------
# merging what several parents propose
# --------------------------------------------------------------------------

def test_a_daughter_both_parents_propose_is_computed_once():
    """Otherwise every branch point doubles the cost for nothing."""
    merged = merge_sequences([
        [{"index": 1, "m": 2, "BS": "", "from": 0}],
        [{"index": 1, "m": 2, "BS": "", "from": 0}, {"index": 2, "m": 4, "BS": "", "from": 1}],
    ])

    assert [(e["m"], e["BS"]) for e in merged] == [(2, ""), (4, "")]


def test_the_predecessor_chain_survives_the_merge():
    """`from` names the orbitals a candidate starts from. Point it at the wrong
    entry and a broken-symmetry solution collapses into the pure one, which is
    the entire reason the chain exists."""
    merged = merge_sequences([
        [{"index": 1, "m": 2, "BS": "", "from": 0}, {"index": 2, "m": 4, "BS": "", "from": 1}],
        [{"index": 1, "m": 4, "BS": "", "from": 0}, {"index": 2, "m": 4, "BS": "3,1", "from": 1}],
    ])

    by_key = {(e["m"], e["BS"]): e for e in merged}
    pure_quartet = by_key[(4, "")]
    bs_quartet = by_key[(4, "3,1")]
    assert bs_quartet["from"] == pure_quartet["index"], "BS must start from its pure counterpart"


def test_indices_are_contiguous_after_merging():
    merged = merge_sequences([
        [{"index": 1, "m": 2, "BS": "", "from": 0}],
        [{"index": 5, "m": 4, "BS": "", "from": 0}],
        [{"index": 9, "m": 6, "BS": "", "from": 0}],
    ])

    assert [e["index"] for e in merged] == [1, 2, 3]


def test_a_from_pointing_at_the_input_geometry_stays_there():
    merged = merge_sequences([[{"index": 3, "m": 2, "BS": "", "from": 0}]])

    assert merged[0]["from"] == 0


def test_merging_nothing_is_not_an_error():
    assert merge_sequences([]) == []
    assert merge_sequences([[]]) == []
