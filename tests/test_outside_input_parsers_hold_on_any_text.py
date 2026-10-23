"""The parsers of input from outside hold on any text: no crash, bounded time.

Input: generated text (hypothesis, ``database=None`` so nothing is written
into the checkout, ``derandomize=True`` so a CI run is reproducible) fed to
the four parsers that read what a user, a workspace or a program wrote:

- CONTROL.txt: ``config.parse_control_text`` (returns a dict) and
  ``config.validate_control_text`` (returns a list of messages; its
  contract is to collect errors, not raise them);
- XYZ coordinates: ``xyz_io._parse_xyz_atoms`` and
  ``xyz_io.normalize_xyz_body`` (lists), and the Submit tab's
  classification ``input_processing.clean_input_data`` (``(text, kind)``,
  kind one of smiles/xyz/empty);
- ORCA output: ``energies.electronic_energy_in_text`` and
  ``energies.gibbs_energy_in_text`` (a float or None);
- the SMILES entry path: ``mol_prep.is_smiles_string`` (bool) and
  ``classify.classify_smiles`` (its features; documented as pure and
  deterministic).

Output: a failure for any exception, a result of another type, or a run
over the per-example deadline (2 s). Separately, each parser is given a
1 MB input built to defeat a backtracking or quadratic scan and must
answer within 5 s: ORCA outputs and CONTROL files arrive from clusters
and other people's workspaces, and a parser that stalls on one stalls the
dashboard or the agent's turn with it.

Measured 2026-10-16 (Python 3.11, this host): the property tests take
about 6 s together, inside the fast suite's budget, so they are not marked
slow. ClusterFuzzLite (coverage-guided, PR mode) is a possible follow-up;
these properties are what it would start from.
"""

from __future__ import annotations

import math
import time

import pytest

hypothesis = pytest.importorskip("hypothesis", reason="hypothesis is in the [dev] extra")
from hypothesis import HealthCheck, given, settings  # noqa: E402
from hypothesis import strategies as st  # noqa: E402

from delfin import classify, config, energies, xyz_io  # noqa: E402
from delfin.dashboard import input_processing  # noqa: E402
from delfin.manta import mol_prep  # noqa: E402

_SETTINGS = settings(max_examples=150, deadline=2000, database=None,
                     derandomize=True,
                     suppress_health_check=[HealthCheck.too_slow])


def _lines(fragments: list[str], *, max_lines: int = 30):
    """Text made of lines joined from ``fragments`` and arbitrary text."""
    token = st.one_of(st.sampled_from(fragments), st.text(max_size=12))
    line = st.lists(token, max_size=8).map("".join)
    return st.lists(line, max_size=max_lines).map("\n".join)


_CONTROL_FRAGMENTS = [
    "charge", "multiplicity", "SMILES", "input_file", "oxidation_steps",
    "reduction_steps", "functional", "PAL", "maxcore", "solvent", "ESD_modul",
    "Species delta", "smiles_converter", "=", "=", ":", " ", "[", "]", "{",
    "}", "(", ")", ",", "'", '"', "-1,0,+1=[", "even_seq", "odd_seq",
    "INFOS:", "ESD MODULE:", "#", "---", "***", "0", "1", "-2", "+3", "1e999",
    "nan", "None", "True", "[" * 60, "(" * 60, "9" * 5000, "yes", "no",
    "[CHARGE]", "%maxcore", "\t", "\r",
]

_XYZ_FRAGMENTS = [
    "C", "Fe", "Fe(1)", "H", "Xx", "*", "$", " ", "  ", "\t", "0", "0.0",
    "-1.25", "1e5", "nan", "inf", "-inf", "1e400", "3", "12", "comment",
    "1.0.0", "--", "+",
]

_ORCA_FRAGMENTS = [
    "FINAL SINGLE POINT ENERGY", "Total Energy", ":", "Electronic energy",
    "Final Gibbs free energy", "Gibbs free energy", "Zero point energy",
    "...", " ", "Eh", "a.u.", "-1234.5678", "1e-3", "E+05", "nan",
    "$$$$$$$$ JOB NUMBER 1 $$$$$$$$", "$$$$$$$$ JOB NUMBER 2 $$$$$$$$",
    "|  1> ", "$new_job", '%base "T1"', '%base "x"', "\n",
]

_SMILES_FRAGMENTS = [
    "C", "c1ccccc1", "[Fe]", "[Fe-2]", "[Cu+]", "(", ")", "=", "#", "[",
    "]", "->", "<-", ".", "@", "@@", "/", "\\", "%10", "1", "2", "[nH]",
    "[Pt](Cl)(Cl)N", "B", "[B-]", "*", "#", "Cl", "Br", "[Ru]", "[U]",
    " ", "\n",
]


@_SETTINGS
@given(_lines(_CONTROL_FRAGMENTS))
def test_control_text_parses_to_a_dict(text):
    assert isinstance(config.parse_control_text(text), dict)


@_SETTINGS
@given(_lines(_CONTROL_FRAGMENTS))
def test_control_text_validates_to_a_list(text):
    found = config.validate_control_text(text, converts_smiles=False)
    assert isinstance(found, list)
    assert all(isinstance(m, str) for m in found)


@_SETTINGS
@given(_lines(_XYZ_FRAGMENTS))
def test_xyz_lines_parse_to_atoms(text):
    lines = text.splitlines()
    atoms = xyz_io._parse_xyz_atoms(lines)
    assert isinstance(atoms, list)
    for atom in atoms:
        assert isinstance(atom["elem"], str)
        assert all(isinstance(atom[k], float) for k in "xyz")
    assert all(isinstance(line, str) for line in xyz_io.normalize_xyz_body(lines))


@_SETTINGS
@given(st.one_of(_lines(_XYZ_FRAGMENTS), _lines(_SMILES_FRAGMENTS, max_lines=3)))
def test_the_submit_box_classifies_any_text(text):
    cleaned, kind = input_processing.clean_input_data(text)
    assert isinstance(cleaned, str)
    assert kind in ("smiles", "xyz", "empty")


@_SETTINGS
@given(_lines(_ORCA_FRAGMENTS, max_lines=40),
       st.sampled_from(["T1.out", "x.out", "", "a/b/T1_TDDFT.out"]))
def test_orca_output_gives_a_float_or_none(text, path):
    for reader in (energies.electronic_energy_in_text, energies.gibbs_energy_in_text):
        value = reader(text, path)
        assert value is None or isinstance(value, float)


@_SETTINGS
@given(_lines(_SMILES_FRAGMENTS, max_lines=3))
def test_smiles_detection_is_a_bool(text):
    assert isinstance(mol_prep.is_smiles_string(text), bool)


@_SETTINGS
@given(_lines(_SMILES_FRAGMENTS, max_lines=1))
def test_smiles_classification_is_total_and_deterministic(text):
    first = classify.classify_smiles(text)
    assert classify.classify_smiles(text) == first


# ---------------------------------------------------------------------------
# bounded time on inputs built against the scan
# ---------------------------------------------------------------------------
_MB = 1 << 20
_LIMIT_S = 5.0


def _timed(call, *args):
    started = time.perf_counter()
    result = call(*args)
    return result, time.perf_counter() - started


@pytest.mark.parametrize("keyword", ["Electronic energy x\n", "Final Gibbs free energy x\n",
                                     "Gibbs free energy x\n"])
def test_an_energy_keyword_without_a_number_is_read_in_bounded_time(keyword):
    """Every occurrence of the keyword started a scan to the end of the text
    for a number, so a megabyte of them took minutes (measured before the
    fix: 400 lines 0.046 s, x4 per doubling)."""
    text = keyword * (_MB // len(keyword))
    reader = (energies.electronic_energy_in_text if keyword.startswith("Electronic")
              else energies.gibbs_energy_in_text)
    value, took = _timed(reader, text, "x.out")
    assert value is None
    assert took < _LIMIT_S, f"{took:.1f} s for 1 MB"


def test_a_number_far_from_its_keyword_is_not_taken():
    """The window that bounds the scan also keeps a number in an unrelated
    part of the file from being read as the energy."""
    text = "Electronic energy\n" + "x" * 5000 + "\n-1234.5 Eh\n"
    assert energies.electronic_energy_in_text(text, "x.out") is None
    near = "Electronic energy         ...    -1234.5678 Eh\n"
    assert energies.electronic_energy_in_text(near, "x.out") == -1234.5678


@pytest.mark.parametrize("text", [
    "[" * _MB, "charge=" + "[" * _MB, ("a=[\n" * (_MB // 4)), "x" * _MB,
    ("charge=1\n" * (_MB // 9)),
], ids=["brackets", "open-list", "open-lists", "one-line", "many-keys"])
def test_a_megabyte_of_control_text_is_parsed_in_bounded_time(text):
    result, took = _timed(config.parse_control_text, text)
    assert isinstance(result, dict)
    assert took < _LIMIT_S, f"{took:.1f} s for 1 MB"


@pytest.mark.parametrize("text", [
    "C 0 0 0\n" * (_MB // 8), ("Fe " + "1.0 " * (_MB // 4)), "*\n" * (_MB // 2),
], ids=["atoms", "one-long-line", "stars"])
def test_a_megabyte_of_coordinates_is_parsed_in_bounded_time(text):
    atoms, took = _timed(xyz_io._parse_xyz_atoms, text.splitlines())
    assert isinstance(atoms, list)
    assert took < _LIMIT_S, f"{took:.1f} s for 1 MB"
    kind, took = _timed(input_processing.clean_input_data, text)
    assert took < _LIMIT_S, f"{took:.1f} s for 1 MB"


@pytest.mark.parametrize("text", ["C" * _MB, "(" * _MB, "[Fe]" * (_MB // 4),
                                  "c1ccccc1" * (_MB // 8)],
                         ids=["chain", "parens", "metals", "rings"])
def test_a_megabyte_of_smiles_is_classified_in_bounded_time(text):
    features, took = _timed(classify.classify_smiles, text)
    assert features is not None
    assert took < _LIMIT_S, f"{took:.1f} s for 1 MB"
    flag, took = _timed(mol_prep.is_smiles_string, text)
    assert isinstance(flag, bool)
    assert took < _LIMIT_S, f"{took:.1f} s for 1 MB"


def test_the_float_regex_still_reads_ordinary_energies():
    """The parsers' ordinary answers, so a fix for time does not move a value."""
    text = ("FINAL SINGLE POINT ENERGY      -1234.567890123\n"
            "Final Gibbs free energy         ...   -1234.12345 Eh\n")
    assert math.isclose(energies.electronic_energy_in_text(text, "x.out"), -1234.567890123)
    assert math.isclose(energies.gibbs_energy_in_text(text, "x.out"), -1234.12345)
