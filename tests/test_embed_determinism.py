"""Determinism tests for the SMILES -> XYZ embed pipeline.

Background
---------
The ``ITER-E_vs_A_anomaly`` forensik found that exception-fallback embed paths
in :mod:`delfin.smiles_converter` bypassed the seeded ETKDG schedule and used
an unseeded RDKit / Open Babel embed.  One tested SMILES produced six unique
geometries in six identical runs, which silently poisoned every pool metric.

DELFIN must be deterministic: ``same SMILES + same code + same env`` has to
produce bit-identical XYZ output.  These tests lock that contract in:

* fallback-path SMILES (KekulizeException / valence errors) are reproducible
  across repeated runs,
* normal (non-fallback) SMILES are unchanged by the determinism fix,
* the SMILES-derived seed helper is itself deterministic and well-formed.
"""

import pytest

from delfin import smiles_converter as sc


# A SMILES that drives the converter into the exception-fallback machinery
# (its anionic pyridyl / aromatic fragments make RDKit's primary embed throw
# KekulizeException, so the seeded fallback paths must take over).  Named by
# the forensik report as the original "6 geometries in 6 runs" reproducer.
_FALLBACK_SMILES = (
    "CC(=O)/C=C(\\[O-])C.[Ir+3]."
    "c1ccc(-c2cccc[n-]2)cc1.c1ccc(-c2cccc[n-]2)cc1"
)

# More fallback-path SMILES: anionic aromatic heterocycles and metallocenes
# whose kekulization / valence handling routes them through the same
# exception-fallback embed code rather than the common seeded path.
_FALLBACK_SMILES_EXTRA = [
    "c1cc[n-]c1.[Na+]",                       # pyrrolyl anion
    "[CH-]1C=CC=C1.[Fe+2].[CH-]1C=CC=C1",     # ferrocene (cyclopentadienyl)
    "c1c[n-]cn1.[K+]",                        # imidazolyl anion
]

# Normal, well-behaved SMILES that take the common seeded ETKDG path and must
# NOT change as a result of the fallback-seed fix.
_NORMAL_SMILES = [
    "CCO",
    "c1ccccc1",
    "CC(=O)O",
    "C1CCCCC1",
    "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
]


pytestmark = pytest.mark.skipif(
    not getattr(sc, "RDKIT_AVAILABLE", False),
    reason="RDKit is required for embed-determinism tests",
)


def _convert(smiles):
    """Run the public converter once and return a hashable result token.

    Returns the XYZ string on success or an ``("ERROR", message)`` tuple on
    failure; both are deterministic-comparable across runs.
    """
    xyz, err = sc.smiles_to_xyz(smiles, apply_uff=False)
    if xyz is not None:
        return xyz
    return ("ERROR", err)


def test_fallback_smiles_three_runs_bit_identical():
    """The forensik reproducer must yield bit-identical output over 3 runs.

    Before the fix this SMILES produced a different geometry on every run
    because its KekulizeException fallback embed was unseeded.
    """
    runs = [_convert(_FALLBACK_SMILES) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "fallback-path SMILES is non-deterministic across 3 runs"
    )


@pytest.mark.parametrize("smiles", _FALLBACK_SMILES_EXTRA)
def test_extra_fallback_smiles_three_runs_bit_identical(smiles):
    """Additional KekulizeException / valence-fallback SMILES are reproducible."""
    runs = [_convert(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        f"fallback-path SMILES {smiles!r} is non-deterministic across 3 runs"
    )


@pytest.mark.parametrize("smiles", _NORMAL_SMILES)
def test_normal_smiles_unchanged_and_deterministic(smiles):
    """Normal (non-fallback) SMILES stay deterministic and succeed.

    The determinism fix only touches fallback-path seeding; the common seeded
    ETKDG path already passed a fixed seed, so normal molecules must continue
    to embed reproducibly and without error.
    """
    runs = [_convert(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        f"normal SMILES {smiles!r} is non-deterministic across 3 runs"
    )
    assert not isinstance(runs[0], tuple), (
        f"normal SMILES {smiles!r} unexpectedly failed to convert: {runs[0]}"
    )


def test_seed_derivation_is_deterministic():
    """``_deterministic_embed_seed`` returns the same seed for the same SMILES.

    The fallback paths derive their RDKit ``randomSeed`` from the SMILES
    string; that derivation must be a pure, stable function of the input so
    repeated runs (and process restarts) reproduce the same geometry.
    """
    for smiles in [_FALLBACK_SMILES, "CCO", "c1ccccc1"]:
        seeds = {sc._deterministic_embed_seed(smiles) for _ in range(5)}
        assert len(seeds) == 1, (
            f"seed derivation is non-deterministic for {smiles!r}: {seeds}"
        )


def test_seed_derivation_is_universal_and_well_formed():
    """The derived seed is non-negative, 31-bit, and varies with the SMILES.

    The derivation must be universal (a function of the string only, no
    per-SMILES special-casing) yet still spread distinct molecules across
    distinct seeds.  ``None`` / empty input falls back to the fixed default.
    """
    seed_cco = sc._deterministic_embed_seed("CCO")
    seed_benzene = sc._deterministic_embed_seed("c1ccccc1")
    assert seed_cco != seed_benzene, "distinct SMILES collapsed to one seed"
    for seed in (seed_cco, seed_benzene):
        assert isinstance(seed, int)
        assert 0 <= seed <= 0x7FFFFFFF, "seed is outside RDKit's 31-bit range"
    assert sc._deterministic_embed_seed(None) == sc._DEFAULT_EMBED_SEED
    assert sc._deterministic_embed_seed("") == sc._DEFAULT_EMBED_SEED


def test_kekulize_exception_path_is_seeded():
    """A KekulizeException-triggering SMILES converts deterministically.

    ``c1cc[n-]c1`` (pyrrolyl anion) cannot be kekulized by RDKit's primary
    embed, so it exercises the unsanitized / dearomatized fallback embed.
    That path must be seeded -> identical result on every run.
    """
    smiles = "c1cc[n-]c1.[Na+]"
    runs = [_convert(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "KekulizeException fallback path is not seeded"
    )


def test_atom_valence_exception_path_is_seeded():
    """An explicit-valence-error SMILES converts deterministically.

    Hypervalent / unusual-valence bracket atoms route through the
    explicit-valence fallback in ``_smiles_to_xyz_unsanitized_fallback``,
    whose ``TypeError`` branch previously dropped ``randomSeed``.  The result
    must now be reproducible across runs.
    """
    # [PH5] has an explicit valence (5) above RDKit's permitted P valence and
    # forces the no-sanitize / valence-error fallback branch.
    smiles = "[PH5]"
    runs = [_convert(smiles) for _ in range(3)]
    assert runs[0] == runs[1] == runs[2], (
        "AtomValenceException fallback path is not seeded"
    )
