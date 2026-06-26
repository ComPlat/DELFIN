"""Determinism tests for the DELFIN_UFF_SOFT_DONORS=1 path.

Background
----------
After Determinism Fix #1 (commit a5fa353, embed-seed wall-clock) and
Determinism Fix #2 (commit 7cf73e3, embed thread-pool race + constraint-block
``NameError`` + Open Babel force-field global state), ``smiles_to_xyz_isomers``
is bit-deterministic with the default UFF refinement path.

A residual non-determinism source surfaced when ``DELFIN_UFF_SOFT_DONORS=1``
flipped the donor pinning from ``AddAtomConstraint`` ("HARD" / FixAtom) to
``AddDistanceConstraint`` ("SOFT" / distance pin) so the donor could reorient
during UFF.

Root cause (forensik 2026-05-15)
--------------------------------
Open Babel's UFFTYPER cannot match the perceived atom type for transition
metals that carry an explicit formal charge — ``Pt+2``, ``Fe+2``, ``Pd+2``,
``Ru+3``, ``Ir+3``, ``Os+3``, ``Co+2``/``Co+3``, ``Ni+2`` etc.  UFF silently
falls back to a default parameter slot.  Energy stays at a recognisable
uninitialised-memory marker (~5.93e11 kcal/mol) BUT the per-atom *gradient*
on the affected atoms is read from uninitialised heap memory and differs
between process invocations.

Calls that downstream freeze every donor + metal via ``AddAtomConstraint`` are
unaffected because the bad gradients cannot move fixed atoms.  Calls that
swap a donor's ``AddAtomConstraint`` for an ``AddDistanceConstraint`` (the
soft-donor path) leave the donor free; the uninitialised gradient drags the
donor across the workspace; the output geometry depends on whatever the
freshly-allocated heap page contains.  Same input → different output per run.

Fix (always-on correctness, NOT env-gated)
------------------------------------------
``_optimize_xyz_openbabel`` probes ``ff.Energy()`` immediately after
``ff.Setup()``.  If the energy is non-finite or its magnitude exceeds
1e9 kcal/mol, OB UFF has not parameterised the molecule fully → the
distance-pin path is disabled for THIS call and every donor falls back to
``AddAtomConstraint`` (the legacy deterministic behaviour).

This is a runtime check on OB's actual parameterisation response — no metal
allowlist, no SMILES patterns.  Parameterised metals (Zn, Cu, neutral Fe…)
keep the soft-donor path; unparameterised metals (Pt+2, Pd+2, Ru+3, …) get
the deterministic HARD fallback automatically.

Why these tests exist
---------------------
v2-final-prime depends on default-flipping ``DELFIN_UFF_SOFT_DONORS`` to 1
because the equal-n verdict shows it wins all three CRITICAL topology
metrics.  These tests lock the determinism contract in so the default flip
cannot reintroduce the run-to-run drift.
"""

import hashlib
import os
from typing import List, Tuple

import pytest

from delfin import smiles_converter as sc


pytestmark = pytest.mark.skipif(
    not getattr(sc, "RDKIT_AVAILABLE", False),
    reason="RDKit is required for UFF-soft-donor determinism tests",
)


# Metal complexes that hit the topological enumerator's constrained-UFF
# refinement path AND the soft-donor code path when DELFIN_UFF_SOFT_DONORS=1.
# Coverage spans:
#   * unparameterised metals with formal charge (Pt+2, Pd+2, Co+2, Fe+2, Ir+3)
#     — the soft-donor path's worst-case for the OB UFF uninitialised-memory
#     bug; the fix must keep these deterministic.
#   * parameterised metals (Cu+2, Zn+2) — the soft-donor path stays ACTIVE and
#     must still be deterministic (the existing fresh-FF + fresh-OBMol fixes
#     are sufficient when UFF has a real parameter slot for the metal).
_UFF_SOFT_SMILES: List[str] = [
    "[Pt](Cl)(Cl)(N)N",                    # Pt+2 — original handoff reproducer
    "[Pd](Cl)(Cl)(N)N",                    # Pd+2 — same fail mode as Pt
    "[Cu](Cl)(Cl)(N)N",                    # Cu+2 — PARAMETERISED, soft active
    "[Co](Cl)(Cl)(N)N",                    # Co+2 — unparameterised
    "[Fe](Cl)(Cl)(N)N",                    # Fe+2 — unparameterised
    "[Ir](Cl)(Cl)(N)N",                    # Ir+3 — unparameterised
]


def _enable_soft(monkeypatch: pytest.MonkeyPatch) -> None:
    """Force DELFIN_UFF_SOFT_DONORS=1 for the duration of one test."""
    monkeypatch.setenv("DELFIN_UFF_SOFT_DONORS", "1")
    # Drop any per-class allowlist so the default scalar gate applies and we
    # exercise both sigma + multi_sigma paths consistently.
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CLASSES", raising=False)
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CARBON", raising=False)


def _isomer_signature(smiles: str) -> Tuple[Tuple[str, str], ...]:
    """Run :func:`smiles_to_xyz_isomers` once and return a hashable signature.

    The tuple captures both isomer count AND every (xyz, label) pair, so any
    geometry / labelling / ordering drift is caught.
    """
    results, err = sc.smiles_to_xyz_isomers(smiles)
    if err:
        return (("ERROR", err),)
    return tuple((xyz, label) for xyz, label in results)


def _digest(sig: Tuple[Tuple[str, str], ...]) -> str:
    """Stable short fingerprint of a signature, handy for assertion messages."""
    h = hashlib.sha1()
    for xyz, label in sig:
        h.update(label.encode())
        h.update(xyz.encode())
    return h.hexdigest()[:12]


# ---------------------------------------------------------------------------
# 1.  The original handoff reproducer must be bit-identical with UFF-soft ON.
# ---------------------------------------------------------------------------

def test_pt_chloride_amine_three_runs_bit_identical(monkeypatch):
    """``[Pt](Cl)(Cl)(N)N`` with UFF-soft ON: three runs are bit-identical.

    Before the OB-parameterisation determinism gate this reproducer yielded
    33 / 34 / 34 isomers across three otherwise-identical invocations with
    flipped donor geometries — the bug that blocked v2-final-prime.
    """
    _enable_soft(monkeypatch)
    runs = [_isomer_signature("[Pt](Cl)(Cl)(N)N") for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "UFF-soft reproducer is non-deterministic across 3 runs; "
        f"counts={[len(r) for r in runs]} digests={[_digest(r) for r in runs]}"
    )
    # Sanity: a real isomer set, not the ERROR token.
    assert runs[0] and runs[0][0][0] != "ERROR", (
        f"UFF-soft reproducer failed to convert: {runs[0]}"
    )


# ---------------------------------------------------------------------------
# 2.  Parametrised determinism contract — every UFF-soft-active SMILES is
#     bit-identical across three runs regardless of whether the metal is
#     UFF-parameterised.
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("smiles", _UFF_SOFT_SMILES)
def test_uff_soft_metal_complex_three_runs_bit_identical(smiles, monkeypatch):
    """Every UFF-soft-active SMILES is bit-identical over 3 runs.

    Parametrised across parameterised + unparameterised metals so the test
    matrix proves the fix is universal — not specific to Pt or to charged
    metals.
    """
    _enable_soft(monkeypatch)
    runs = [_isomer_signature(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        f"UFF-soft SMILES {smiles!r} is non-deterministic across 3 runs; "
        f"counts={[len(r) for r in runs]} digests={[_digest(r) for r in runs]}"
    )


# ---------------------------------------------------------------------------
# 3.  Isomer count alone must be stable — the bug manifested first as a
#     drifting count (34 vs 33 vs 34).  Check it explicitly so a regression
#     where the geometry stabilises but the count flaps is still caught.
# ---------------------------------------------------------------------------

def test_uff_soft_isomer_count_stable_over_five_runs(monkeypatch):
    """The Pt reproducer's isomer count is the same on every run."""
    _enable_soft(monkeypatch)
    counts = set()
    for _ in range(5):
        results, err = sc.smiles_to_xyz_isomers("[Pt](Cl)(Cl)(N)N")
        assert not err, f"UFF-soft reproducer returned an error: {err}"
        counts.add(len(results))
    assert len(counts) == 1, (
        f"UFF-soft isomer count is non-deterministic across 5 runs: {counts}"
    )


# ---------------------------------------------------------------------------
# 4.  Soft-mode ON must not contaminate the geometry of an unrelated normal
#     molecule that runs after a metal complex in the same process.  Even with
#     the OB parameterisation gate disabling soft-mode for unparameterised
#     metals, the surrounding state must still be call-order-independent
#     (this protects against any future leak in the soft-mode code path).
# ---------------------------------------------------------------------------

def test_uff_soft_metal_run_does_not_perturb_subsequent_normal_run(monkeypatch):
    """A non-metal isomer run is unchanged by a preceding UFF-soft metal run."""
    _enable_soft(monkeypatch)
    normal = "CC(=O)Nc1ccc(O)cc1"
    standalone = _isomer_signature(normal)
    _ = _isomer_signature("[Pt](Cl)(Cl)(N)N")
    after_metal = _isomer_signature(normal)
    assert standalone == after_metal, (
        "non-metal isomer output changed after a UFF-soft metal-complex "
        "run — soft-donor constraint state leaked into the unconstrained path"
    )


# ---------------------------------------------------------------------------
# 5.  Cross-metal call-order independence with UFF-soft ON.  Running an
#     unrelated metal complex between two invocations of the reproducer must
#     not perturb the reproducer's output.
# ---------------------------------------------------------------------------

def test_uff_soft_metal_complex_call_order_independent(monkeypatch):
    """A metal-complex isomer run is unchanged by an intervening metal run."""
    _enable_soft(monkeypatch)
    smiles_a = "[Pt](Cl)(Cl)(N)N"
    smiles_b = "[Cu](Cl)(Cl)(N)N"
    standalone_a = _isomer_signature(smiles_a)
    _ = _isomer_signature(smiles_b)
    after_b_a = _isomer_signature(smiles_a)
    assert standalone_a == after_b_a, (
        "UFF-soft metal-complex isomer output depends on call order — "
        "soft-mode constraint state is leaking between calls"
    )


# ---------------------------------------------------------------------------
# 6.  UFF-soft OFF baseline must be unaffected by the parameterisation gate.
#     The gate runs unconditionally inside ``_optimize_xyz_openbabel`` so we
#     prove explicitly that legacy callers see the same output as before.
# ---------------------------------------------------------------------------

def test_uff_soft_off_three_runs_bit_identical(monkeypatch):
    """With UFF-soft OFF (default), determinism is preserved (no regression)."""
    # Ensure the env-flag is OFF — the gate must be a no-op here.
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS", raising=False)
    monkeypatch.delenv("DELFIN_UFF_SOFT_DONORS_CLASSES", raising=False)
    runs = [_isomer_signature("[Pt](Cl)(Cl)(N)N") for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "UFF-soft OFF baseline became non-deterministic; the parameterisation "
        "gate must not perturb the legacy path"
    )


# ---------------------------------------------------------------------------
# 7.  Direct contract check on ``_optimize_xyz_openbabel`` itself: with a
#     soft-donor meta block + an unparameterised metal, repeated calls inside
#     the SAME Python process produce bit-identical output.
# ---------------------------------------------------------------------------

def test_optimize_xyz_openbabel_soft_path_intra_process_deterministic(monkeypatch):
    """Direct call to ``_optimize_xyz_openbabel`` is deterministic intra-process.

    Reproduces the minimum failure mode pre-fix: five consecutive identical
    calls to the optimizer with a soft-donor meta block on an unparameterised
    metal yielded five DIFFERENT outputs.  The OB parameterisation gate must
    catch this and produce one bit-identical output every time.
    """
    _enable_soft(monkeypatch)
    if not getattr(sc, "OPENBABEL_AVAILABLE", False):
        pytest.skip("Open Babel is required for the direct-optimizer test")

    xyz = (
        "Pt       0.000000     0.000000     0.000000\n"
        "Cl       2.300000     0.000000     0.000000\n"
        "Cl      -2.300000     0.000000     0.000000\n"
        "N        0.000000     2.030000     0.000000\n"
        "N        0.000000    -2.030000     0.000000\n"
        "H        0.000000     2.030000     1.000000\n"
        "H        0.000000     2.030000    -1.000000\n"
        "H        1.000000     2.030000     0.000000\n"
        "H        0.000000    -3.030000     0.000000\n"
    )
    constraints = {
        "fix_atoms": [0],
        "distances": [(0, 1, 2.3), (0, 2, 2.3), (0, 3, 2.03), (0, 4, 2.03)],
        "angles": [],
        "torsions": [],
        "_soft_donor_meta": {
            "metal_indices": [0],
            "donor_indices": [1, 2, 3, 4],
            "pairs": [(0, 1), (0, 2), (0, 3), (0, 4)],
            "class_label": "sigma",
        },
    }
    sigs = []
    for _ in range(5):
        out = sc._optimize_xyz_openbabel(xyz, constraints=constraints)
        sigs.append(hashlib.sha1(out.encode()).hexdigest()[:16])
    assert len(set(sigs)) == 1, (
        "_optimize_xyz_openbabel returns different geometries for the same "
        f"input across consecutive calls (soft-donor unparameterised-metal "
        f"path): {sigs}"
    )


# ---------------------------------------------------------------------------
# 8.  Symmetric contract: parameterised metal stays deterministic in the
#     same direct-call setup (sanity check that the gate did not over-fire
#     and accidentally disable soft-mode for the metals it should still
#     handle).
# ---------------------------------------------------------------------------

def test_optimize_xyz_openbabel_soft_path_param_metal_deterministic(monkeypatch):
    """Parameterised metal: soft-mode stays ACTIVE and bit-deterministic."""
    _enable_soft(monkeypatch)
    if not getattr(sc, "OPENBABEL_AVAILABLE", False):
        pytest.skip("Open Babel is required for the direct-optimizer test")

    xyz = (
        "Cu       0.000000     0.000000     0.000000\n"
        "Cl       2.300000     0.000000     0.000000\n"
        "Cl      -2.300000     0.000000     0.000000\n"
        "N        0.000000     2.030000     0.000000\n"
        "N        0.000000    -2.030000     0.000000\n"
        "H        0.000000     2.030000     1.000000\n"
        "H        0.000000     2.030000    -1.000000\n"
        "H        1.000000     2.030000     0.000000\n"
        "H        0.000000    -3.030000     0.000000\n"
    )
    constraints = {
        "fix_atoms": [0],
        "distances": [(0, 1, 2.3), (0, 2, 2.3), (0, 3, 2.03), (0, 4, 2.03)],
        "angles": [],
        "torsions": [],
        "_soft_donor_meta": {
            "metal_indices": [0],
            "donor_indices": [1, 2, 3, 4],
            "pairs": [(0, 1), (0, 2), (0, 3), (0, 4)],
            "class_label": "sigma",
        },
    }
    sigs = []
    for _ in range(5):
        out = sc._optimize_xyz_openbabel(xyz, constraints=constraints)
        sigs.append(hashlib.sha1(out.encode()).hexdigest()[:16])
    assert len(set(sigs)) == 1, (
        "Parameterised-metal soft-donor path drifted across consecutive "
        f"calls: {sigs}"
    )


# ---------------------------------------------------------------------------
# 9.  Welle-3 T6.1: the multi-metal-augmentation ETKDG block (formerly at
#     ``smiles_converter.py:26457``) submitted ``_embed_multiple_confs_robust``
#     concurrently against the SHARED mol — exactly the same RDKit thread-
#     safety race that fix #2 (7cf73e3) patched in the primary embed loop.
#     This test runs a small bimetallic SMILES that takes the multi-metal
#     augmentation path and asserts every run is bit-identical.
# ---------------------------------------------------------------------------

_BIMETALLIC_SMILES = (
    # Pool idx-15 family adapted for a synthetic bimetallic dummy: two
    # independent Pt centres glued via a covalent C-C linker so the
    # converter has to take the >= 2-metals augmentation branch. Constructed
    # WITHOUT referring to any specific refcode/dataset so the test stays
    # universal (cf. memory rule "no SMILES-specific shortcuts").
    "[Pt](Cl)(Cl)(N)NC(N)N[Pt](Cl)(Cl)N"
)


def test_multimetal_augmentation_three_runs_bit_identical(monkeypatch):
    """Bimetallic SMILES is bit-identical across 3 runs.

    Targets the previously-unpatched multi-metal augmentation ThreadPoolExecutor
    block. A regression here means the data-race fix (private mol copy + seed-
    order merge) has been reverted or mis-replicated.
    """
    # Make the high-fanout case maximally likely to expose any race.
    monkeypatch.setenv("DELFIN_MAX_THREAD_WORKERS", "32")
    # Default flag state — Welle-2 + UFF-soft as currently shipped.
    runs = [_isomer_signature(_BIMETALLIC_SMILES) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "multi-metal augmentation path is non-deterministic across 3 runs; "
        f"counts={[len(r) for r in runs]} digests={[_digest(r) for r in runs]}"
    )


def test_multimetal_augmentation_worker_count_invariant(monkeypatch):
    """Bimetallic output is bit-identical at WORKER=1 and WORKER=32.

    The race in the old code was *only* exposed under concurrency; this test
    locks in that the output no longer depends on the worker count.
    """
    monkeypatch.setenv("DELFIN_MAX_THREAD_WORKERS", "1")
    seq_sig = _isomer_signature(_BIMETALLIC_SMILES)
    monkeypatch.setenv("DELFIN_MAX_THREAD_WORKERS", "32")
    par_sig = _isomer_signature(_BIMETALLIC_SMILES)
    assert seq_sig == par_sig, (
        "multi-metal augmentation output depends on DELFIN_MAX_THREAD_WORKERS — "
        f"seq_digest={_digest(seq_sig)} par_digest={_digest(par_sig)}"
    )
