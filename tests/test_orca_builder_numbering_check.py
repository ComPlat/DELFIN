"""ORCA Builder: superposition for the numbering check, and the atom-number size.

The superposition is what "Check Numbering" reports as its RMSD and what the
"Aligned reference" view draws, so a rotation applied the wrong way round shows
a structure as several Angstrom away from an exact copy of itself.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from delfin.dashboard import structure_editor
from delfin.dashboard.tab_orca_builder import _orca_kabsch_align


def _z_rotation(angle):
    c, s = np.cos(angle), np.sin(angle)
    return np.array([[c, -s, 0.0], [s, c, 0.0], [0.0, 0.0, 1.0]])


def _cloud(n=30, seed=0):
    return np.random.default_rng(seed).normal(size=(n, 3)) * 3.0


def test_alignment_recovers_an_exact_rotation():
    ref = _cloud()
    target = ref @ _z_rotation(0.7).T + np.array([1.0, 2.0, 3.0])

    aligned, rmsd = _orca_kabsch_align(ref, target)

    assert rmsd == pytest.approx(0.0, abs=1e-9)
    assert np.allclose(aligned, target, atol=1e-9)


def test_alignment_does_not_absorb_a_mirror_image():
    # Only proper rotations are allowed: an enantiomer must not superpose.
    ref = _cloud()
    mirrored = ref * np.array([1.0, 1.0, -1.0])

    _aligned, rmsd = _orca_kabsch_align(ref, mirrored)

    assert rmsd > 1.0


def test_alignment_applies_the_mapping():
    ref = _cloud()
    perm = np.random.default_rng(1).permutation(len(ref))
    shuffled = ref[perm]                      # target atom k holds reference perm[k]
    mapping = np.argsort(perm)                # target index for each reference atom

    _aligned, rmsd = _orca_kabsch_align(ref, shuffled, mapping=mapping)

    assert rmsd == pytest.approx(0.0, abs=1e-9)


def test_alignment_rejects_a_mismatched_mapping_length():
    ref = _cloud(n=5)
    with pytest.raises(ValueError):
        _orca_kabsch_align(ref, ref, mapping=np.arange(4))


# ---------------------------------------------------------------------------
# atom-number size control
# ---------------------------------------------------------------------------
def _build_tab(tmp_path):
    from delfin.dashboard import tab_orca_builder

    executed_js = []
    ctx = SimpleNamespace(
        backend=SimpleNamespace(),
        calc_dir=Path(tmp_path),
        run_js=executed_js.append,
        add_init_js=executed_js.append,
    )
    _tab, widgets_map = tab_orca_builder.create_tab(ctx)
    return widgets_map, executed_js


def test_atom_number_size_resizes_without_redrawing_the_molecule(tmp_path):
    pytest.importorskip("ipywidgets")
    widgets_map, executed_js = _build_tab(tmp_path)
    size = widgets_map["submit_label_size"]

    assert size.value == structure_editor.LABEL_PX_DEFAULT
    assert size.min == structure_editor.LABEL_PX_MIN
    assert size.max == structure_editor.LABEL_PX_MAX

    executed_js.clear()
    size.value = structure_editor.LABEL_PX_MAX

    # one call into the page that rescales the existing labels -- no new viewer
    assert len(executed_js) == 1
    assert "__delfinAtomNumbers.setScale(" in executed_js[0]
    assert "_submitMolViewerByScope" in executed_js[0]
    assert "addModel" not in executed_js[0]


def test_a_saddle_search_is_offered_and_is_given_what_it_needs(tmp_path):
    """OptTS is here because there is now something to give it.

    The editor's scan walks a reaction and its path finder hands back an
    estimated transition state; what that estimate wants next is a real
    optimisation to a first-order saddle, and that is ORCA's.

    It is not Opt with a different name.  A saddle search follows one mode
    uphill and every other down, so it has to know which mode -- which means
    a Hessian.  Without one it starts from the identity matrix and walks to a
    minimum or to nothing, which is a job that runs for hours and answers a
    different question.  Choosing it turns Calc_Hess on, and recalculates
    every five steps because the mode being followed changes character as the
    structure moves.
    """
    pytest.importorskip("ipywidgets")
    widgets_map, _js = _build_tab(tmp_path)
    kinds = list(widgets_map["orca_job_type"].options)
    assert "OPTTS" in kinds and "OPTTS FREQ" in kinds, kinds
    # And the ordinary ones are still there: this is a switch, not a swap.
    assert "OPT" in kinds and "OPT FREQ" in kinds, kinds

    widgets_map["orca_coords"].value = "3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"
    build = widgets_map["generate_orca_input"]

    widgets_map["orca_job_type"].value = "OPT"
    plain = build()
    assert "%geom" not in plain, plain
    assert " OPT " in plain

    widgets_map["orca_job_type"].value = "OPTTS"
    saddle = build()
    assert " OPTTS " in saddle
    assert "%geom\n  Calc_Hess true\n  Recalc_Hess 5\nend" in saddle, saddle


def test_a_saddle_search_and_a_held_value_share_one_geom_block(tmp_path):
    """ORCA reads one %geom.  Two in an input is one of them thrown away
    silently, and which one is not something a user can see."""
    pytest.importorskip("ipywidgets")
    widgets_map, _js = _build_tab(tmp_path)
    widgets_map["orca_coords"].value = "3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"
    widgets_map["editor_state"]["constraints"] = [
        {"kind": "distance", "atoms": [0, 1], "value": 1.2, "mode": "fix"}]
    build = widgets_map["generate_orca_input"]

    # The hold alone, as before.
    widgets_map["orca_job_type"].value = "OPT"
    held = build()
    assert held.count("%geom") == 1
    assert "Constraints" in held and "Calc_Hess" not in held

    # And with a saddle search, both, in one block.
    widgets_map["orca_job_type"].value = "OPTTS FREQ"
    both = build()
    assert both.count("%geom") == 1, both
    assert "Calc_Hess true" in both
    assert "  Constraints" in both
    assert "{ B 0 1 1.2000 C }" in both
    # And it closes once for the constraints and once for the block.
    inside = both.split("%geom", 1)[1].split("* xyz", 1)[0]
    assert inside.count("end") == 2, inside
