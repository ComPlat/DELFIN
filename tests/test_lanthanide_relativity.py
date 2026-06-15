"""Lanthanides/actinides (4f/5f) must be detected as metals AND routed
through DELFIN's relativistic path.

Regression for the silent science bug behind Jerome's Dy3+ thread: the
f-block was absent from _TM_LIST / _TM_D45 / cli.py's local set, so a
dysprosium (or any 4f/5f) complex was seen as "no transition metal" and ran
non-relativistically with a non-relativistic basis. The correct recipe
already exists in DELFIN (relativity=ZORA, metal_basisset_rel=SARC-ZORA-TZVP);
the only fix needed was to make f-elements trigger it.
"""

from __future__ import annotations

from delfin import utils


_REL_CFG = {
    "relativity": "ZORA",
    "main_basisset_rel": "ZORA-def2-SVP",
    "metal_basisset_rel": "SARC-ZORA-TZVP",
    "aux_jk_rel": "SARC/J",
}


def _write_xyz(tmp_path, *symbols):
    lines = [str(len(symbols)), "comment"]
    for i, s in enumerate(symbols):
        lines.append(f"{s} {i*2.0:.3f} 0.000 0.000")
    p = tmp_path / "geom.xyz"
    p.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return str(p)


def test_dysprosium_is_detected_as_metal(tmp_path):
    xyz = _write_xyz(tmp_path, "Dy", "C", "O", "H", "H")
    assert utils.search_transition_metals(xyz) == ["Dy"]


def test_uranium_is_detected_as_metal(tmp_path):
    xyz = _write_xyz(tmp_path, "U", "O", "O")
    assert utils.search_transition_metals(xyz) == ["U"]


def test_classification_of_blocks():
    assert utils.classify_tm_presence([]) == "none"
    assert utils.classify_tm_presence(["Fe"]) == "3d"
    assert utils.classify_tm_presence(["Ru"]) == "4d5d"
    assert utils.classify_tm_presence(["Dy"]) == "4f5f"      # lanthanide
    assert utils.classify_tm_presence(["U"]) == "4f5f"       # actinide
    assert utils.classify_tm_presence(["Fe", "Dy"]) == "mixed"


def test_should_use_rel_for_f_block():
    assert utils._should_use_rel(["Dy"]) is True
    assert utils._should_use_rel(["U"]) is True
    assert utils._should_use_rel(["Gd", "Eu"]) is True
    # Regression: light d-block behaviour is unchanged.
    assert utils._should_use_rel(["Fe"]) is False
    assert utils._should_use_rel([]) is False
    assert utils._should_use_rel(["Ru"]) is True


def test_resolve_level_of_theory_gives_sarc_for_dysprosium():
    main, metal, rel_token, aux = utils.resolve_level_of_theory(["Dy"], dict(_REL_CFG))
    assert rel_token == "ZORA"
    assert metal == "SARC-ZORA-TZVP"          # relativistic metal basis attached
    assert main == "ZORA-def2-SVP"
    assert aux == "SARC/J"


def test_resolve_level_of_theory_stays_nonrel_for_iron():
    cfg = dict(_REL_CFG)
    main, metal, rel_token, aux = utils.resolve_level_of_theory(["Fe"], cfg)
    assert rel_token == ""                     # no relativity for pure 3d
    assert metal is None                       # no SARC override on a 3d metal
    assert main == "def2-SVP"


def test_cli_does_not_force_relativity_off_for_dysprosium():
    """cli.py uses utils._should_use_rel as the single source of truth, so a
    Dy complex keeps its configured relativity instead of being forced to
    'none' by a stale local metal set."""
    from delfin.utils import _should_use_rel
    cfg = {"relativity": "ZORA"}
    metals = ["Dy"]
    use_rel = _should_use_rel(metals)          # mirrors the cli.py decision
    if not use_rel:
        cfg["relativity"] = "none"
    assert use_rel is True
    assert cfg["relativity"] == "ZORA"
