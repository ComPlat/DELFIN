"""ORCA Builder: superposition for the numbering check, and the atom-number size.

The superposition is what "Check Numbering" reports as its RMSD and what the
"Aligned reference" view draws, so a rotation applied the wrong way round shows
a structure as several Angstrom away from an exact copy of itself.
"""

from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from delfin.dashboard.tab_orca_builder import (
    _LABEL_SCALE_DEFAULT,
    _orca_kabsch_align,
)


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
    size = widgets_map["orca_mol_label_size"]

    assert size.value == _LABEL_SCALE_DEFAULT
    assert dict(size.options)["XXL"] > dict(size.options)["S"]

    executed_js.clear()
    size.value = dict(size.options)["XXL"]

    # one call into the page that rescales the existing labels -- no new viewer
    assert len(executed_js) == 1
    assert "__delfinAtomNumbers.setScale(window._orcaBuildViewer," in executed_js[0]
    assert "addModel" not in executed_js[0]
