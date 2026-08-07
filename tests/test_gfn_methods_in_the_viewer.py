"""GFN next to UFF: what Optimise minimises with, and what it needs to know.

UFF and MMFF94 live in the browser, which is what lets a drag follow the mouse.
GFN runs xtb on the server -- too slow for a drag, right for a press -- so
choosing it changes what Optimise does, not what dragging does.  And xtb needs
a charge and a spin: the charge can be read off a SMILES, the spin never can.
"""

from __future__ import annotations

import shutil

import pytest

from delfin.dashboard import gfn_optimize as gfn
from delfin.dashboard.context import DashboardContext

_WATER = "3\nwater\nO 0.0 0.0 0.0\nH 0.96 0.0 0.0\nH -0.24 0.93 0.0\n"
_needs_xtb = pytest.mark.skipif(not shutil.which("xtb"), reason="xtb not installed")


# ---------------------------------------------------------------------------
# the runner
# ---------------------------------------------------------------------------
def test_the_methods_offered_are_the_ones_xtb_knows():
    assert set(gfn.GFN_METHODS) == {"gfnff", "gfn2", "gfn1"}
    assert gfn.is_gfn_method("gfnff") and gfn.is_gfn_method("GFN2")
    assert not gfn.is_gfn_method("uff") and not gfn.is_gfn_method("")


def test_a_method_it_does_not_know_is_refused_not_guessed():
    result = gfn.optimize_with_gfn(_WATER, "b3lyp")
    assert result["ok"] is False
    assert "not a GFN method" in result["status"]
    assert result["xyz"] == _WATER


def test_a_structure_too_large_says_so_instead_of_starting():
    big = "5000\nbig\n" + "\n".join(f"C {i} 0 0" for i in range(5000))
    result = gfn.optimize_with_gfn(big, "gfn2")
    assert result["ok"] is False
    assert "past the GFN2-xTB limit" in result["status"]
    assert "submit it as a job" in result["status"]


def test_nothing_to_optimise_is_not_an_error_worth_running():
    assert gfn.optimize_with_gfn("1\nx\nH 0 0 0\n", "gfnff")["ok"] is False


@_needs_xtb
def test_it_relaxes_and_says_what_it_cost():
    result = gfn.optimize_with_gfn(_WATER, "gfnff", charge=0, uhf=0)

    assert result["ok"] is True
    assert result["xyz"] != _WATER, "the geometry came back unchanged"
    assert result["xyz"].splitlines()[0].strip() == "3"
    assert result["energy"] is not None
    assert result["seconds"] > 0
    assert "converged" in result["status"]


@_needs_xtb
def test_the_charge_and_the_spin_reach_xtb():
    """A different charge has to give a different energy, or they were dropped."""
    neutral = gfn.optimize_with_gfn(_WATER, "gfn2", charge=0, uhf=0)
    cation = gfn.optimize_with_gfn(_WATER, "gfn2", charge=1, uhf=1)

    assert neutral["ok"] and cation["ok"]
    assert neutral["energy"] != cation["energy"]


# ---------------------------------------------------------------------------
# in the tab
# ---------------------------------------------------------------------------
@pytest.fixture
def editor(tmp_path):
    pytest.importorskip("ipywidgets")
    from delfin.dashboard import tab_submit

    for name in ("calc", "archive", "office"):
        (tmp_path / name).mkdir()
    sent: list[str] = []
    ctx = DashboardContext(
        calc_dir=tmp_path / "calc",
        archive_dir=tmp_path / "archive",
        office_dir=tmp_path / "office",
    )
    ctx.run_js = sent.append
    _widget, refs = tab_submit.create_tab(ctx)
    refs["coords_widget"].value = _WATER
    return refs


def test_the_methods_stand_next_to_uff(editor):
    values = [v for _label, v in editor["submit_ff_dd"].options]
    assert values == ["uff", "mmff94", "gfnff", "gfn2"]


def test_charge_and_spin_appear_only_for_a_gfn_method(editor):
    assert editor["submit_gfn_charge"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfnff"
    assert editor["submit_gfn_charge"].layout.display == ""
    assert editor["submit_gfn_mult"].layout.display == ""

    editor["submit_ff_dd"].value = "uff"
    assert editor["submit_gfn_charge"].layout.display == "none"


def test_the_charge_is_read_off_the_smiles_when_there_is_one(editor):
    editor["editor_state"]["converted_xyz_cache"] = {
        "smiles": "[Fe+2]([NH3])([NH3])([NH3])([NH3])([NH3])[NH3]",
        "xyz": _WATER,
    }

    editor["submit_ff_dd"].value = "gfnff"

    assert editor["submit_gfn_charge"].value == 2


def test_without_a_smiles_the_charge_is_the_users_to_set(editor):
    editor["editor_state"]["converted_xyz_cache"] = {"smiles": None, "xyz": _WATER}
    editor["submit_gfn_charge"].value = -1

    editor["submit_ff_dd"].value = "gfn2"

    assert editor["submit_gfn_charge"].value == -1, "a pasted XYZ says nothing"


def test_dragging_keeps_a_force_field_that_lives_in_the_browser(editor):
    from delfin.dashboard import tab_submit

    source = open(tab_submit.__file__, encoding="utf-8").read()
    live = source.split("def _live_ff_method")[1].split("\n    def ")[0]
    assert "'uff' if _gfn.is_gfn_method(chosen) else chosen" in live
    # and the live export uses that, not the dropdown
    assert "method=_live_ff_method()," in source


# ---------------------------------------------------------------------------
# it really is xtb, and it says which Hamiltonian it ran
# ---------------------------------------------------------------------------
@_needs_xtb
def test_the_result_names_the_program_that_produced_it():
    """Passing --gfn 2 and being given GFN2 are two different claims."""
    result = gfn.optimize_with_gfn(_WATER, "gfn2")

    assert result["engine"] == "xtb"
    assert result["version"], "xtb did not report a version"
    assert result["hamiltonian"] == "GFN2-xTB"
    assert "xtb" in result["status"]


@_needs_xtb
def test_each_method_reports_its_own_hamiltonian():
    assert gfn.optimize_with_gfn(_WATER, "gfnff")["hamiltonian"] == "GFN-FF"
    assert gfn.optimize_with_gfn(_WATER, "gfn1")["hamiltonian"] == "GFN1-xTB"


# ---------------------------------------------------------------------------
# autospin
# ---------------------------------------------------------------------------
def test_the_parity_decides_which_multiplicities_are_possible():
    # water: 10 electrons, even -> singlet, triplet, quintet
    assert gfn.electron_parity(_WATER, 0) == 0
    # take one electron away and it can only be a doublet, quartet ...
    assert gfn.electron_parity(_WATER, 1) == 1


@_needs_xtb
def test_autospin_keeps_the_multiplicity_that_came_out_lowest():
    result = gfn.optimize_autospin(_WATER, "gfn2", charge=0)

    assert result["ok"] is True
    assert result["multiplicity"] == 1, "water is a singlet"
    assert len(result["tried"]) == 3
    energies = [e for _m, e, ok in result["tried"] if ok and e is not None]
    assert result["energy"] == min(energies)
    assert "Lowest of 3 multiplicities" in result["status"]


@_needs_xtb
def test_autospin_scans_the_parity_the_charge_implies():
    result = gfn.optimize_autospin(_WATER, "gfnff", charge=1)

    assert [m for m, _e, _ok in result["tried"]] == [2, 4, 6], (
        "an odd electron count cannot be a singlet"
    )


def test_the_checkbox_appears_with_the_method_and_takes_over_the_box(editor):
    assert editor["submit_gfn_autospin"].layout.display == "none"

    editor["submit_ff_dd"].value = "gfn2"
    assert editor["submit_gfn_autospin"].layout.display == ""

    editor["submit_gfn_autospin"].value = True
    assert editor["submit_gfn_mult"].disabled is True, (
        "a scanned multiplicity is not one the user is setting"
    )
