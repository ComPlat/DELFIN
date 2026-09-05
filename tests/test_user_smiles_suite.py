"""User-shared SMILES regression anchors.

This file collects every SMILES the user has sent so far in the
conversation history.  Each entry defines a baseline invariant that
must hold regardless of later refactors:

- ``min_isomers`` — a lower bound on the number of distinct isomers.
  Completeness work on the topology pipeline should **only grow**
  this count; a regression drops below it.
- ``required_labels`` — a set of labels that must appear in the
  output.  These capture constitutional isomers the user has
  explicitly asked to see (e.g. N-N-ax for Cd-histidine).
- ``forbidden_labels`` — labels that must NOT appear (e.g. no
  ``alt-C-isomer`` for Ir(ppy)2(acac)).

Determinism is asserted globally: two consecutive runs must return
identical label sets and identical counts.

New SMILES from the user should be added to :data:`USER_SMILES` as
additional dicts.  The suite does not encode exact isomer counts in
every case — the regression contract is directional (count never
shrinks, required labels always present) so the completeness work
can keep improving without tripping the tests.
"""
from __future__ import annotations

import os
import pytest

from delfin.smiles_converter import smiles_to_xyz_isomers


# ---------------------------------------------------------------------------
# User-shared SMILES collected from the development conversation.
# ---------------------------------------------------------------------------
USER_SMILES = [
    dict(
        name="Cd-histidine CN=7",
        smiles=(
            "[OH2+][Cd-5]12([OH2+])(OC(C3=C(C(O)=O)NC=[N+]13)=O)"
            "(OC(C4=C(C(O)=O)NC=[N+]24)=O)[OH2+]"
        ),
        min_isomers=5,
        # At least one PBP axial variant with two chelate-O trans
        # and one with N-N axial (from sampling classify labels) must surface.
        required_label_fragments=["-ax"],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Cd MA2B2C2 (triazolothiadiazine N4O2Cl2)",
        smiles=(
            "C[OH+][Cd-4](Cl)(Cl)([OH+]C)"
            "([N+]1=C2SC(C)=NN2C(C3=CC=CC=C3)=N1)"
            "([N+]4=C5SC(C)=NN5C(C6=CC=CC=C6)=N4)"
        ),
        min_isomers=5,
        required_label_fragments=[],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Ir(ppy)2(acac) CN=6",
        smiles=(
            "CC1=CC(C)=[O+][Ir-3]2([N+]3=C4C=CC=C3)"
            "(C5=CC=CC=C54)(O1)"
            "[N+]6=CC=CC=C6C7=C2C=CC=C7"
        ),
        min_isomers=3,
        required_label_fragments=["C-trans", "N-trans", "all-cis"],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Fe(CO)3(NHC)2 CN=5 — neutral C",
        smiles=(
            "O#[C][Fe-2]([C]#O)(C1=[N+](C)C=CN1C)([C]#O)"
            "C2=[N+](C)C=CN2C"
        ),
        min_isomers=3,
        required_label_fragments=["C0-C0-ax", "C1-C1-ax"],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Fe(CO)3(NHC)2 CN=5 — [C+]/[Fe-3] variant",
        smiles=(
            "O#[C+][Fe-3]([C+]#O)(C1=[N+](C)C=CN1C)([C+]#O)"
            "C2=[N+](C)C=CN2C"
        ),
        min_isomers=3,
        required_label_fragments=["C0-C0-ax", "C1-C1-ax"],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Fe(CO)3(NHC)2 CN=5 — [C+]/[Fe-5] variant",
        smiles=(
            "O#[C+][Fe-5]([C+]#O)(C1=[N+](C)C=CN1C)([C+]#O)"
            "C2=[N+](C)C=CN2C"
        ),
        min_isomers=3,
        required_label_fragments=["C0-C0-ax", "C1-C1-ax"],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Co(en)2Cl2 bidentate",
        smiles="Cl[Co+3]12(Cl)(NCCN1)NCCN2",
        min_isomers=2,
        required_label_fragments=[],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Fe2 (mu-Cl)2 bimetallic",
        smiles="[Fe+2]1(N)(N)(Cl[Fe+2](N)(N)Cl1)Cl",
        min_isomers=1,
        required_label_fragments=[],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Fe(H2O)7 CN=7 homoleptic",
        smiles="[Fe+2](O)(O)(O)(O)(O)(O)O",
        min_isomers=1,
        required_label_fragments=[],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Zr(H2O)8 CN=8 homoleptic",
        smiles="[Zr+4](O)(O)(O)(O)(O)(O)(O)O",
        min_isomers=1,
        required_label_fragments=[],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Fe/Sc(OTf)4(OH)(mu-O) cyclam bimetal",
        smiles=(
            "O=S(O[Sc](OS(=O)(C(F)(F)F)=O)(OS(=O)(C(F)(F)F)=O)"
            "(OS(=O)(C(F)(F)F)=O)(O)"
            "O[Fe-3]123[N@@+]4(C)CCC[N@@+]1(CC[N@@+]2(CCC[N@+]3(C)CC4)C)C)"
            "(C(F)(F)F)=O"
        ),
        # The user reports an earlier version produced more options.
        # Keep an honest floor while completeness work is ongoing.
        min_isomers=1,
        required_label_fragments=[],
        forbidden_label_fragments=[],
    ),
    dict(
        name="Cd MA2B2C2 octahedral (five OH isomers)",
        smiles=(
            "C[OH+][Cd-4]([Cl])([Cl])([OH+]C)"
            "([N+]1=C2SC(C)=NN2C(C2=CC=CC=C2)=N1)"
            "[N+]1=C2SC(C)=NN2C(C2=CC=CC=C2)=N1"
        ),
        min_isomers=5,
        required_label_fragments=["all-cis", "all-trans", "N-trans", "O-trans", "Cl-trans"],
        forbidden_label_fragments=[],
    ),
]


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------
def _run(smi: str):
    res, err = smiles_to_xyz_isomers(
        smi,
        apply_uff=True,
        deterministic=True,
        collapse_label_variants=True,
    )
    assert err is None, f"smiles_to_xyz_isomers returned error: {err}"
    return res


# One build, shared by every test that only reads it.
#
# Four tests below ask four different questions of the same structures, and
# each of them used to build those structures again from scratch. For the
# Fe/Sc entry that is four times twelve minutes to answer four questions
# about one answer. Determinism is the one property that genuinely needs a
# second build, and it still takes one — so an entry now costs two builds
# instead of five, and the nightly run shrinks with the gate.
_BUILD_CACHE: dict = {}


def _built(smi: str):
    if smi not in _BUILD_CACHE:
        _BUILD_CACHE[smi] = _run(smi)
    return _BUILD_CACHE[smi]


# ---------------------------------------------------------------------------
# What runs on every push, and what waits for the nightly run
# ---------------------------------------------------------------------------
# Seconds for ONE build, measured 2026-08-14. The spread is the whole point:
# a single entry costs more than eleven minutes while three others cost about
# a second each.
#
#   702.1  Fe/Sc(OTf)4(OH)(mu-O) cyclam bimetal
#   184.9  Ir(ppy)2(acac) CN=6
#   157.9  Cd-histidine CN=7
#    89.2  Cd MA2B2C2 octahedral (five OH isomers)
#    44.2  Cd MA2B2C2 (triazolothiadiazine N4O2Cl2)
#    31.0  Fe2 (mu-Cl)2 bimetallic
#    13.7  Co(en)2Cl2 bidentate
#     7.6  Zr(H2O)8 CN=8 homoleptic
#     4.3  Fe(H2O)7 CN=7 homoleptic
#     1.2  Fe(CO)3(NHC)2 CN=5, all three variants
#
# This file was marked slow as a whole, so the one twelve-minute case exiled
# eleven others with it and the SMILES contract was checked nowhere except a
# nightly run. The mark belongs on the entry.
#
# The numbers are here rather than measured at collection time because a
# selection that re-decides itself on every machine is not a contract: the
# gate would quietly shrink on a loaded runner and nobody would be told.
# Re-measure deliberately when the converter's cost changes.
_BUILD_SECONDS = {
    "Cd-histidine CN=7": 157.9,
    "Cd MA2B2C2 (triazolothiadiazine N4O2Cl2)": 44.2,
    "Ir(ppy)2(acac) CN=6": 184.9,
    "Fe(CO)3(NHC)2 CN=5 — neutral C": 1.2,
    "Fe(CO)3(NHC)2 CN=5 — [C+]/[Fe-3] variant": 1.2,
    "Fe(CO)3(NHC)2 CN=5 — [C+]/[Fe-5] variant": 1.3,
    "Co(en)2Cl2 bidentate": 13.7,
    "Fe2 (mu-Cl)2 bimetallic": 31.0,
    "Fe(H2O)7 CN=7 homoleptic": 4.3,
    "Zr(H2O)8 CN=8 homoleptic": 7.6,
    "Fe/Sc(OTf)4(OH)(mu-O) cyclam bimetal": 702.1,
    "Cd MA2B2C2 octahedral (five OH isomers)": 89.2,
}

# Per entry, for one build.
_GATE_BUDGET_S = 10.0

# Cheap enough for the gate, and held back anyway because they are red today.
#
# Both fail test_topology_invariants_for_every_output: the structure the
# converter emits for a homoleptic aqua complex does not satisfy
# _verify_topology_from_graph. They are two of the thirteen the nightly run
# reports, they are in the converter rather than in this suite, and a merge
# gate that is red on arrival teaches people to ignore it.
#
# Not marked xfail, which would be a claim about how they end. They are held
# out until somebody fixes them; then deleting a line here is the whole
# change, and the gate grows from about 8 s to about 31 s.
_GATE_HELD_BACK = {
    "Fe(H2O)7 CN=7 homoleptic",
    "Zr(H2O)8 CN=8 homoleptic",
}


def _in_gate(entry) -> bool:
    """Unmeasured entries wait for the nightly run.

    A new SMILES is expensive until somebody has shown otherwise, which is
    the safe direction: the cost of being wrong here is a slower merge gate
    for everyone.
    """
    if entry["name"] in _GATE_HELD_BACK:
        return False
    return _BUILD_SECONDS.get(entry["name"], float("inf")) <= _GATE_BUDGET_S


def _params(entries=None):
    entries = USER_SMILES if entries is None else entries
    return [
        pytest.param(
            e,
            id=e["name"],
            marks=[] if _in_gate(e) else [pytest.mark.slow],
        )
        for e in entries
    ]


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize("entry", _params())
def test_min_isomer_floor(entry):
    """The output must contain at least ``min_isomers`` distinct entries."""
    res = _built(entry["smiles"])
    assert len(res) >= entry["min_isomers"], (
        f"{entry['name']!r}: only {len(res)} isomers returned, "
        f"expected >= {entry['min_isomers']}: {[l for _, l in res]}"
    )


@pytest.mark.parametrize(
    "entry",
    _params([e for e in USER_SMILES if e.get("required_label_fragments")]),
)
def test_required_labels_present(entry):
    """Each required label fragment must be a substring of at least one output label."""
    res = _built(entry["smiles"])
    labels = [l for _, l in res]
    for needle in entry["required_label_fragments"]:
        assert any(needle in l for l in labels), (
            f"{entry['name']!r}: required label fragment {needle!r} "
            f"missing from output {labels}"
        )


@pytest.mark.parametrize(
    "entry",
    _params([e for e in USER_SMILES if e.get("forbidden_label_fragments")]),
)
def test_forbidden_labels_absent(entry):
    """Forbidden label fragments must never appear in the output."""
    res = _built(entry["smiles"])
    labels = [l for _, l in res]
    for needle in entry["forbidden_label_fragments"]:
        offenders = [l for l in labels if needle in l]
        assert not offenders, (
            f"{entry['name']!r}: forbidden fragment {needle!r} "
            f"appeared in {offenders}"
        )


@pytest.mark.parametrize("entry", _params())
def test_determinism_across_runs(entry):
    """Two consecutive runs must return the same (sorted) label set and count."""
    r1 = _built(entry["smiles"])
    r2 = _run(entry["smiles"])
    s1 = sorted(l for _, l in r1)
    s2 = sorted(l for _, l in r2)
    assert s1 == s2, (
        f"{entry['name']!r} not deterministic:\n  run1: {s1}\n  run2: {s2}"
    )
    assert len(r1) == len(r2), (
        f"{entry['name']!r}: count drift {len(r1)} vs {len(r2)}"
    )


@pytest.mark.parametrize("entry", _params())
def test_topology_invariants_for_every_output(entry):
    """Every output XYZ must pass the graph-based topology gate.

    This is the universal guard rail: no matter which branch produced
    the structure (sampling, topo enumeration, multi-metal coupled
    build, linkage rewiring), ``_verify_topology_from_graph`` must
    succeed.  If it does not, the pipeline has leaked a broken bond
    through the gate.
    """
    from rdkit import Chem
    from delfin.smiles_converter import (
        _verify_topology_from_graph,
        _normalize_metal_smiles,
    )

    smi = entry["smiles"]
    norm = _normalize_metal_smiles(smi) or smi
    mol = Chem.MolFromSmiles(norm)
    if mol is None:
        pytest.skip(f"{entry['name']!r}: SMILES failed to parse for template")
    mol = Chem.AddHs(mol)

    res = _built(smi)
    for xyz, lbl in res:
        assert _verify_topology_from_graph(xyz, mol), (
            f"{entry['name']!r}: output isomer {lbl!r} fails graph gate"
        )


# ---------------------------------------------------------------------------
# The split itself
# ---------------------------------------------------------------------------
# The value of the arrangement above is not in any one line of it, so these
# guard the shape: that something cheap actually reaches the merge gate, that
# it stays cheap, and that the whole-file mark does not come back and quietly
# take the SMILES contract out of every push again.

def test_the_gate_actually_gets_some_of_this_suite():
    """A split that promotes nothing is the state this replaces, and it looks
    identical from the outside: green gate, contract unchecked."""
    in_gate = [e["name"] for e in USER_SMILES if _in_gate(e)]

    assert in_gate, (
        "no entry runs on push, so the SMILES contract is once again checked "
        "nowhere but the nightly run"
    )


def test_the_gate_set_stays_within_its_budget():
    """Two builds per entry: one shared by the read-only tests, one more for
    the determinism comparison. Measured at 8.3 s for the current set."""
    total = 2 * sum(_BUILD_SECONDS[e["name"]] for e in USER_SMILES if _in_gate(e))

    assert total <= 45.0, (
        f"the gate portion of this file would cost about {total:.0f} s. The "
        "fast suite is what people wait on before a merge; move an entry back "
        "or raise this deliberately."
    )


def test_no_entry_reaches_the_gate_without_a_measured_cost():
    """The budget above is only meaningful while every promoted entry has a
    number behind it. Unmeasured means slow, and that is the safe direction."""
    for e in USER_SMILES:
        if _in_gate(e):
            assert e["name"] in _BUILD_SECONDS, e["name"]


def test_the_cost_table_still_describes_this_suite():
    """A renamed entry silently loses its measurement and drops to the
    nightly run, which is safe but invisible. Say it instead."""
    names = {e["name"] for e in USER_SMILES}
    stale = set(_BUILD_SECONDS) - names

    assert not stale, f"_BUILD_SECONDS names no entry in USER_SMILES: {stale}"


def test_the_whole_file_is_not_marked_slow_again():
    """conftest marks by FILENAME. Putting this file back in that set would
    undo the split without touching anything here — which is exactly how the
    twelve ended up behind one mark in the first place."""
    import tests.conftest as _cf

    assert "test_user_smiles_suite.py" not in _cf._SLOW_TEST_FILES, (
        "this file marks itself per entry; a whole-file mark exiles the cheap "
        "ones along with the twelve-minute one"
    )


def test_the_expensive_entries_did_not_quietly_join_the_gate():
    """The point of the split is that a twelve-minute build stays out of the
    path people wait on."""
    for e in USER_SMILES:
        cost = _BUILD_SECONDS.get(e["name"], float("inf"))
        if cost > _GATE_BUDGET_S:
            assert not _in_gate(e), f"{e['name']} costs {cost} s"
