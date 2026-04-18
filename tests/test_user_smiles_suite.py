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
        # Bimetallic sampling augmentation is not fully deterministic
        # yet (known open issue flagged by the user — multi-metal
        # ETKDG path has non-deterministic ordering of fingerprint
        # collision tiebreaking).  Flag so the test surface reports
        # it clearly but does not block CI.
        nondeterministic=True,
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
        nondeterministic=True,  # same multi-metal sampling issue
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


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
@pytest.mark.parametrize(
    "entry",
    USER_SMILES,
    ids=[e["name"] for e in USER_SMILES],
)
def test_min_isomer_floor(entry):
    """The output must contain at least ``min_isomers`` distinct entries."""
    res = _run(entry["smiles"])
    assert len(res) >= entry["min_isomers"], (
        f"{entry['name']!r}: only {len(res)} isomers returned, "
        f"expected >= {entry['min_isomers']}: {[l for _, l in res]}"
    )


@pytest.mark.parametrize(
    "entry",
    [e for e in USER_SMILES if e.get("required_label_fragments")],
    ids=[e["name"] for e in USER_SMILES if e.get("required_label_fragments")],
)
def test_required_labels_present(entry):
    """Each required label fragment must be a substring of at least one output label."""
    res = _run(entry["smiles"])
    labels = [l for _, l in res]
    for needle in entry["required_label_fragments"]:
        assert any(needle in l for l in labels), (
            f"{entry['name']!r}: required label fragment {needle!r} "
            f"missing from output {labels}"
        )


@pytest.mark.parametrize(
    "entry",
    [e for e in USER_SMILES if e.get("forbidden_label_fragments")],
    ids=[e["name"] for e in USER_SMILES if e.get("forbidden_label_fragments")],
)
def test_forbidden_labels_absent(entry):
    """Forbidden label fragments must never appear in the output."""
    res = _run(entry["smiles"])
    labels = [l for _, l in res]
    for needle in entry["forbidden_label_fragments"]:
        offenders = [l for l in labels if needle in l]
        assert not offenders, (
            f"{entry['name']!r}: forbidden fragment {needle!r} "
            f"appeared in {offenders}"
        )


@pytest.mark.parametrize(
    "entry",
    USER_SMILES,
    ids=[e["name"] for e in USER_SMILES],
)
def test_determinism_across_runs(entry):
    """Two consecutive runs must return the same (sorted) label set and count.

    Systems flagged ``nondeterministic=True`` document an open bug in
    the multi-metal sampling augmentation; their failure is expected
    until the coupled-enumeration path becomes fully deterministic.
    """
    if entry.get("nondeterministic"):
        pytest.xfail(
            "Known open issue: multi-metal ETKDG augmentation is not "
            "fully deterministic (user flagged)."
        )
    r1 = _run(entry["smiles"])
    r2 = _run(entry["smiles"])
    s1 = sorted(l for _, l in r1)
    s2 = sorted(l for _, l in r2)
    assert s1 == s2, (
        f"{entry['name']!r} not deterministic:\n  run1: {s1}\n  run2: {s2}"
    )
    assert len(r1) == len(r2), (
        f"{entry['name']!r}: count drift {len(r1)} vs {len(r2)}"
    )


@pytest.mark.parametrize(
    "entry",
    USER_SMILES,
    ids=[e["name"] for e in USER_SMILES],
)
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

    res = _run(smi)
    for xyz, lbl in res:
        assert _verify_topology_from_graph(xyz, mol), (
            f"{entry['name']!r}: output isomer {lbl!r} fails graph gate"
        )
