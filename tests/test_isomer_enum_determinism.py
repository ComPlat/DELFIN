"""Determinism tests for the topological isomer enumeration pipeline.

Background
----------
After the embed-seed determinism fix (commit a5fa353) three residual
non-determinism sources remained in :func:`smiles_to_xyz_isomers`, all on the
metal-complex / topological-enumeration path:

1. **Embed thread-pool data race** — ``smiles_to_xyz_isomers`` embeds several
   ETKDG seeds in parallel via ``ThreadPoolExecutor``, but every worker called
   ``AllChem.EmbedMultipleConfs`` against the *same shared* ``mol``.  RDKit
   conformer embedding is not thread-safe on a shared molecule, so the
   conformer pool (and therefore the isomer count / geometry) drifted between
   runs.  Fix: each worker embeds into its own ``Chem.Mol`` copy and the main
   thread merges the conformers back in seed order.

2. **Constraint-block ``NameError``** — the Wave 1-7 soft-donor patch added a
   reference to a ``mol`` variable inside ``_optimize_xyz_openbabel`` that does
   not exist in that function's scope.  Every constrained UFF call raised
   ``NameError``, which the broad ``except`` around the constraint block
   swallowed — so ``ff.SetConstraints`` was never reached and every
   "constrained" optimisation silently ran *unconstrained*.  Unconstrained UFF
   on a metal complex wanders to a different geometry each call.  Fix: pass
   ``None`` (the documented fail-safe value) instead of the undefined ``mol``.

3. **Open Babel force-field global state** — ``_optimize_xyz_openbabel`` reused
   the process-global ``pybel._forcefields["uff"]`` object, which retains the
   previous call's constraint set.  A constrained metal optimisation left its
   fixed atoms / distance pins behind, so the *next* call (even an unrelated
   molecule) inherited them, making output depend on call order.  Fix: take a
   fresh force-field instance per call and explicitly install an empty
   constraint set on unconstrained calls.

DELFIN must be deterministic: ``same SMILES + same code + same env`` has to
produce bit-identical isomer output.  These tests lock that contract in for
the topological-enumeration path.
"""

import pytest

from delfin import smiles_converter as sc


# The original reproducer from the determinism handoff: a pyrrolyl-anion +
# Na(+) ionic complex that drove the embed thread-pool race.  Before the fix
# it yielded a varying isomer count (1 or 2) and varying geometry across
# identical runs.
_REPRODUCER_SMILES = "c1cc[n-]c1.[Na+]"

# Metal complexes that exercise the topological enumeration path
# (`_generate_topological_isomers` and the constrained-UFF refinement inside
# it).  Each must produce a bit-identical isomer set on every run.
_TOPO_ENUM_SMILES = [
    "c1cc[n-]c1.[Na+]",                              # ionic, embed-pool race
    "[Pt](Cl)(Cl)(N)N",                              # CN4, constrained UFF
    "C(=O)([O-])[O-].[Ca+2]",                        # carboxylate-Ca ionic
    "O=C([O-])C[NH2+]CC(=O)[O-].[Cu+2]",             # glycine-Cu chelate
    "[Co+3].[NH3:1].[NH3:1].[NH3:1].[Cl-].[Cl-].[Cl-]",  # Co(III) ammine
]

# Normal (non-metal) SMILES.  These do not hit the topological enumerator at
# all; they are included to prove the fix does not perturb the common path and
# that an isomer run is not contaminated by a preceding metal-complex run
# (Open Babel force-field global-state leak).
_NORMAL_SMILES = [
    "CCO",
    "c1ccccc1",
    "CC(=O)Nc1ccc(O)cc1",
    "O=C(O)c1ccccc1C(=O)O",
    "C1CCCCC1",
]


pytestmark = pytest.mark.skipif(
    not getattr(sc, "RDKIT_AVAILABLE", False),
    reason="RDKit is required for isomer-enumeration determinism tests",
)


def _isomer_token(smiles):
    """Run ``smiles_to_xyz_isomers`` once and return a hashable, comparable token.

    The token captures both the isomer count and every (xyz, label) pair, so
    any drift in count *or* geometry *or* labelling is detected.
    """
    results, err = sc.smiles_to_xyz_isomers(smiles)
    if err:
        return ("ERROR", err)
    # tuple-of-tuples is hashable and compares element-wise (bit-identical XYZ).
    return tuple((xyz, label) for xyz, label in results)


def test_reproducer_three_runs_bit_identical():
    """The handoff reproducer must yield a bit-identical isomer set over 3 runs.

    Before the embed thread-pool fix this SMILES produced 1 or 2 isomers with
    varying geometry depending on thread interleaving.
    """
    runs = [_isomer_token(_REPRODUCER_SMILES) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "topological-enumeration reproducer is non-deterministic across 3 runs"
    )
    # Sanity: it must actually produce at least one isomer (not an error token).
    assert not (runs[0] and runs[0][0] == "ERROR"), (
        f"reproducer failed to convert: {runs[0]}"
    )


@pytest.mark.parametrize("smiles", _TOPO_ENUM_SMILES)
def test_topo_enum_smiles_three_runs_bit_identical(smiles):
    """Every topological-enumeration SMILES is reproducible across 3 runs.

    Covers the embed thread-pool race, the constrained-UFF ``NameError`` and
    the Open Babel force-field state leak — all three surfaced as a drifting
    isomer count / geometry on these complexes.
    """
    runs = [_isomer_token(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        f"topological-enumeration SMILES {smiles!r} is non-deterministic "
        f"across 3 runs (counts: {[len(r) for r in runs]})"
    )


@pytest.mark.parametrize("smiles", _NORMAL_SMILES)
def test_normal_smiles_isomers_three_runs_bit_identical(smiles):
    """Normal (non-metal) SMILES stay deterministic through the isomer entry point.

    The determinism fixes only touch the metal / topological-enumeration path;
    normal molecules must continue to enumerate reproducibly and succeed.
    """
    runs = [_isomer_token(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        f"normal SMILES {smiles!r} is non-deterministic across 3 runs"
    )
    assert runs[0] and runs[0][0] != "ERROR", (
        f"normal SMILES {smiles!r} unexpectedly failed to convert: {runs[0]}"
    )


def test_isomer_count_is_stable_for_reproducer():
    """The reproducer's isomer *count* is stable (not just the geometry).

    The embed thread-pool race manifested most visibly as a flapping isomer
    count (1 vs 2) for ``c1cc[n-]c1.[Na+]``; this test pins the count directly
    so a regression in the enumeration cardinality is caught even if the XYZ
    comparison were ever loosened.
    """
    counts = set()
    for _ in range(5):
        results, err = sc.smiles_to_xyz_isomers(_REPRODUCER_SMILES)
        assert not err, f"reproducer returned an error: {err}"
        counts.add(len(results))
    assert len(counts) == 1, (
        f"reproducer isomer count is non-deterministic across 5 runs: {counts}"
    )


def test_metal_complex_isomers_order_independent():
    """A metal complex run is not contaminated by a preceding metal-complex run.

    The Open Babel force-field reused a process-global object that kept the
    previous call's constraints alive.  This test runs one metal complex, then
    a *different* metal complex, then re-runs the first — the first complex's
    isomer set must be identical to its standalone result regardless of what
    ran in between.
    """
    smiles_a = "[Pt](Cl)(Cl)(N)N"
    smiles_b = "O=C([O-])C[NH2+]CC(=O)[O-].[Cu+2]"

    standalone_a = _isomer_token(smiles_a)
    _ = _isomer_token(smiles_b)            # unrelated constrained-UFF run
    after_b_a = _isomer_token(smiles_a)

    assert standalone_a == after_b_a, (
        "metal-complex isomer output depends on call order — Open Babel "
        "force-field global state is leaking between calls"
    )


def test_normal_smiles_not_contaminated_by_metal_run():
    """A non-metal isomer run is unaffected by a preceding metal-complex run.

    Direct guard for the force-field constraint-state leak: a constrained
    metal UFF optimisation must not change the geometry of a subsequent
    unconstrained, non-metal molecule.
    """
    normal = "CC(=O)Nc1ccc(O)cc1"

    standalone = _isomer_token(normal)
    _ = _isomer_token("[Pt](Cl)(Cl)(N)N")  # constrained metal UFF run
    after_metal = _isomer_token(normal)

    assert standalone == after_metal, (
        "non-metal isomer output changed after a metal-complex run — "
        "Open Babel constraint state leaked into the unconstrained path"
    )
