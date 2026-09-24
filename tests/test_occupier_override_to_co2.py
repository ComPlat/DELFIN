"""OCCUPIER manual-override → CO2 coordinator spin handoff tests.

A dashboard ``--occupier-override`` retags the winning FSPE entry in
OCCUPIER.txt from ``<-- PREFERRED VALUE`` to ``<-- OVERRIDE``. Two
consumers must see the override:

1. ``extract_preferred_spin`` (delfin/copy_helpers.py) — the OCCUPIER.txt
   fallback used by the CO2 chain when no state file exists.
2. ``.delfin_occ_auto_state.json`` — refreshed via
   ``record_auto_preference`` so the primary source
   (``_spin_from_state_json``) carries the override too.
"""

from __future__ import annotations

from pathlib import Path

from delfin.copy_helpers import extract_preferred_spin
from delfin.occupier_auto import _load_state, record_auto_preference, infer_parity_from_m
from delfin.occupier_sequences import infer_species_delta


OCCUPIER_TXT_WITH_OVERRIDE = """\
Method: PBE0 def2-SVP
Charge: 0
-------------
FINAL SINGLE POINT ENERGY (1)       = -2609.340 (H)
multiplicity 1
----------------------------------------------------------------
FINAL SINGLE POINT ENERGY (2)       = -2609.352 (H) <-- OVERRIDE
multiplicity 1, BrokenSym 1,1
----------------------------------------------------------------
FINAL SINGLE POINT ENERGY (3)       = -2609.352 (H) <-- was PREFERRED (pre-override)
multiplicity 3
----------------------------------------------------------------
Preferred Index: 2
"""

OCCUPIER_TXT_WITHOUT_OVERRIDE = """\
Method: PBE0 def2-SVP
Charge: 0
-------------
FINAL SINGLE POINT ENERGY (2)       = -2609.352 (H) <-- PREFERRED VALUE
multiplicity 3
----------------------------------------------------------------
Preferred Index: 2
"""


def test_extract_preferred_spin_honors_override_marker(tmp_path: Path):
    """The OVERRIDE tag marks the winner after a manual override."""
    folder = tmp_path / "red_step_2_OCCUPIER"
    folder.mkdir()
    (folder / "OCCUPIER.txt").write_text(OCCUPIER_TXT_WITH_OVERRIDE, encoding="utf-8")
    mult, bs = extract_preferred_spin(folder)
    assert mult == 1
    assert bs == "1,1"


def test_extract_preferred_spin_still_reads_plain_preferred(tmp_path: Path):
    """Unchanged behavior for files without an override."""
    folder = tmp_path / "red_step_2_OCCUPIER"
    folder.mkdir()
    (folder / "OCCUPIER.txt").write_text(OCCUPIER_TXT_WITHOUT_OVERRIDE, encoding="utf-8")
    mult, bs = extract_preferred_spin(folder)
    assert mult == 3
    # No BrokenSym in this file → the helper returns None for bs
    assert not bs


def test_record_auto_preference_updates_state_json(tmp_path: Path):
    """The state file must carry the override winner (m + BS), so the
    CO2 chain's primary spin source sees it."""
    # Simulate the pre-override state: m=3 won for delta=-2
    record_auto_preference("odd", 3, -2, m_value=3, bs_value="", root=tmp_path)
    state = _load_state(tmp_path)
    assert state["-2"]["odd"]["m"] == 3

    # Now the manual override picks index 2 (m=1, BS 1,1)
    parity = infer_parity_from_m(1)
    record_auto_preference(parity, 2, -2, m_value=1, bs_value="1,1", root=tmp_path)

    state = _load_state(tmp_path)
    assert state["-2"][parity]["m"] == 1
    assert state["-2"][parity]["BS"] == "1,1"
    assert state["-2"][parity]["index"] == 2


def test_infer_species_delta_from_folder_name(tmp_path: Path):
    """Folder-name-based delta inference for the override call site."""
    assert infer_species_delta(tmp_path / "red_step_2_OCCUPIER") == -2
    assert infer_species_delta(tmp_path / "ox_step_1_OCCUPIER") == 1
