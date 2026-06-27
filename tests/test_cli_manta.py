"""Fast unit tests for the delfin-manta CLI helpers and argument contract.

These do NOT invoke the converter (no structure building) so they stay fast and
CI-cheap. They lock in two invariants that matter for the completeness claim and
byte-identity: the CLI defaults to ALL isomers, and ranking is OFF by default.
"""

from delfin.cli_manta import _safe_name, _atom_lines, _to_xyz, _build_parser


def test_safe_name_sanitizes_labels():
    assert _safe_name("fac-TBP-5-1", 3) == "003__fac-TBP-5-1.xyz"
    assert _safe_name("trans-SQ Cl0+Cl0|N1", 0).endswith(".xyz")
    assert _safe_name("", 7) == "007__isomer.xyz"


def test_atom_lines_bare_block():
    assert _atom_lines("C 0 0 0\nH 1 0 0\n") == ["C 0 0 0", "H 1 0 0"]


def test_atom_lines_strips_existing_xyz_header():
    assert _atom_lines("2\ncomment\nC 0 0 0\nH 1 0 0\n") == ["C 0 0 0", "H 1 0 0"]


def test_to_xyz_adds_standard_header():
    lines = _to_xyz("C 0 0 0\nO 1 0 0\n", "ethanol | CCO").splitlines()
    assert lines[0] == "2"
    assert lines[1] == "ethanol | CCO"
    assert lines[2] == "C 0 0 0"


def test_to_xyz_idempotent_on_standard_xyz():
    lines = _to_xyz("2\nc\nC 0 0 0\nO 1 0 0", "x").splitlines()
    assert lines[0] == "2" and lines[2] == "C 0 0 0"


def test_default_emits_all_isomers():
    # completeness claim: no isomer cap unless the user asks for one
    args = _build_parser().parse_args(["CCO"])
    assert args.max_isomers is None


def test_ranking_off_by_default():
    # byte-identity: GFN2/GFN-FF ranking must be opt-in
    args = _build_parser().parse_args(["CCO"])
    assert args.rank is False


def test_rank_and_method_flags():
    args = _build_parser().parse_args(["CCO", "--rank", "--method", "gfnff"])
    assert args.rank is True and args.method == "gfnff"
