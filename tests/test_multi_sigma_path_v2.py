"""Unit tests for the multi_sigma V2 path (Iter-multi_sigma audit).

Patch summary
-------------

The ``multi_sigma`` class (>=2 metals, σ-only coordination) is the worst
class on the 2026-05-13 voll-pool: 23 out of 33 SMILES (70%) hit the
external 600 s timeout, lowering coverage to 30.3% (Mode A) / 36.4%
(Mode B).  Forensik shows the wall-clock is consumed by a fixed seed
schedule (20 top-level + 10 multi-metal augmentation) × 25 s per-seed
``_MULTIEMBED_TIMEOUT`` on 50-125-atom complexes — the per-conformer
ETKDG cost scales super-linearly with atom count.

The V2 path tightens three knobs class-conditionally:

  * top-level ETKDG seed count          (``_qprof["seeds"]``)
  * multi-metal augmentation seed count (``_extra_seed_count``)
  * per-seed ETKDG timeout              (``_MULTIEMBED_TIMEOUT_OVERRIDE``)

A wall-clock budget is also enforced on the multi-metal augmentation
block so a single hung embedding cannot starve the remaining seeds.

Defaults
--------

* ``DELFIN_MULTI_SIGMA_PATH_V2`` defaults to ``0`` (OFF).
* When OFF, every helper is bit-exact pre-patch.
* When ON, the budget applies *only* to SMILES classified as
  ``multi_sigma`` by ``_classify_complex_class``; other classes are
  untouched.

The tests below cover that contract.  They never call the full
``smiles_to_xyz_isomers`` pipeline (expensive, RDKit/ETKDG bound) —
they exercise the public helpers directly so the CI run stays fast.
"""
from __future__ import annotations

import os
from typing import Iterable

import pytest

pytest.importorskip("rdkit", reason="RDKit required for multi_sigma V2 tests")
from rdkit import Chem  # noqa: E402

from delfin.smiles_converter import (  # noqa: E402
    DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES,
    _MULTIEMBED_TIMEOUT_OVERRIDE,
    _classify_complex_class,
    _multi_sigma_v2_active,
    _multi_sigma_v2_budget,
)


# ----------------------------------------------------------------------
# Fixtures
# ----------------------------------------------------------------------

#: SMILES samples per chemistry class.  Picked to keep parsing cheap and
#: to match the labels emitted by ``_classify_complex_class``.
_FIXTURES = {
    "sigma":       "[Pt](Cl)(Cl)(N)N",
    "multi_sigma": "[Mn]([C]#O)([C]#O)[Mn]([C]#O)([C]#O)",
    "no_metal":    "CC(=O)O",
}


def _mol(smi: str):
    m = Chem.MolFromSmiles(smi, sanitize=False)
    assert m is not None, f"SMILES did not parse: {smi}"
    return m


def _scrub_env(*names: Iterable[str]) -> None:
    for n in names:
        os.environ.pop(n, None)


@pytest.fixture(autouse=True)
def _clean_v2_env():
    """Strip V2 env-vars + thread-local override before every test."""
    _scrub_env(
        "DELFIN_MULTI_SIGMA_PATH_V2",
        "DELFIN_MULTI_SIGMA_PATH_V2_CLASSES",
    )
    if hasattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value"):
        delattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value")
    yield
    _scrub_env(
        "DELFIN_MULTI_SIGMA_PATH_V2",
        "DELFIN_MULTI_SIGMA_PATH_V2_CLASSES",
    )
    if hasattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value"):
        delattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value")


# ----------------------------------------------------------------------
# Test 1: budget is monotone-decreasing in heavy-atom count.
# ----------------------------------------------------------------------
def test_budget_monotone_in_size() -> None:
    """Larger molecules MUST get smaller seed pools + tighter timeouts.

    Otherwise the V2 path could spend more wall-clock on a 125-atom system
    than on a 20-atom one, which is the failure mode the patch fixes.
    """
    sizes = [20, 50, 80, 125]
    budgets = [_multi_sigma_v2_budget(n) for n in sizes]

    seeds_top = [b["seeds_top"] for b in budgets]
    seeds_aug = [b["seeds_mm_aug"] for b in budgets]
    timeouts = [b["embed_timeout"] for b in budgets]
    walltimes = [b["mm_walltime"] for b in budgets]

    assert seeds_top == sorted(seeds_top, reverse=True), (
        f"seeds_top must be non-increasing with size, got {seeds_top}"
    )
    assert seeds_aug == sorted(seeds_aug, reverse=True), (
        f"seeds_mm_aug must be non-increasing with size, got {seeds_aug}"
    )
    assert timeouts == sorted(timeouts, reverse=True), (
        f"embed_timeout must be non-increasing with size, got {timeouts}"
    )
    assert walltimes == sorted(walltimes, reverse=True), (
        f"mm_walltime must be non-increasing with size, got {walltimes}"
    )

    # Large-system worst case must stay well below the external 600 s pool
    # cap so the converter can still complete topology + UFF + final gate.
    largest = budgets[-1]
    assert largest["seeds_top"] <= 5, largest
    assert largest["mm_walltime"] <= 60.0, largest
    assert largest["embed_timeout"] <= 15.0, largest


# ----------------------------------------------------------------------
# Test 2: class-conditional default-ON contract (Phase 1.5 wire-in).
# ----------------------------------------------------------------------
def test_v2_default_on_for_multi_sigma_only() -> None:
    """With no env-var set, ``_multi_sigma_v2_active`` returns True for
    SMILES classified as ``multi_sigma`` and False for every other class.

    This is the Phase 1.5 wire-in contract: the patch is auto-active by
    default for the multi_sigma class (to recover the 70% timeout-driven
    coverage loss documented in ITER-multi_sigma_audit.md), but every
    other class remains bit-exact pre-patch so sigma/hapto pools cannot
    silently regress.
    """
    for class_label, smi in _FIXTURES.items():
        mol = _mol(smi)
        assert _classify_complex_class(mol) in {class_label, "no_metal"}, (
            f"fixture mismatch for {class_label!r}: classifier says "
            f"{_classify_complex_class(mol)!r}"
        )
        observed = _multi_sigma_v2_active(mol)
        if class_label in DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES:
            assert observed is True, (
                f"V2 must be ON-by-default for class={class_label!r} "
                "(Phase 1.5 wire-in contract)"
            )
        else:
            assert observed is False, (
                f"V2 must stay OFF-by-default for class={class_label!r} "
                "(only multi_sigma is auto-activated)"
            )

    # mol=None must never crash, must return False.
    assert _multi_sigma_v2_active(None) is False


# ----------------------------------------------------------------------
# Test 3: explicit OFF is the global escape hatch.
# ----------------------------------------------------------------------
def test_v2_explicit_off_disables_for_all_classes() -> None:
    """``DELFIN_MULTI_SIGMA_PATH_V2=0`` is the operator escape hatch —
    must disable V2 for every class including multi_sigma, restoring the
    pre-wire-in pipeline bit-exact.
    """
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2"] = "0"

    for class_label, smi in _FIXTURES.items():
        mol = _mol(smi)
        assert _multi_sigma_v2_active(mol) is False, (
            f"V2 must be OFF when scalar=0 for class={class_label!r}"
        )


# ----------------------------------------------------------------------
# Test 4: explicit ON is equivalent to the default for multi_sigma.
# ----------------------------------------------------------------------
def test_v2_scalar_flag_one_keeps_class_conditional_behaviour() -> None:
    """``DELFIN_MULTI_SIGMA_PATH_V2=1`` keeps the class-conditional
    semantics (same as unset): activates only for the default class set.

    sigma / no_metal must stay OFF — the patch is class-conditional,
    not a global toggle, regardless of the scalar value.
    """
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2"] = "1"
    assert "multi_sigma" in set(DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES), (
        "DEFAULT_CLASSES must include multi_sigma — that is the whole point"
    )

    for class_label, smi in _FIXTURES.items():
        mol = _mol(smi)
        observed = _multi_sigma_v2_active(mol)
        if class_label in DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES:
            assert observed is True, (
                f"V2 should be ON for class={class_label!r} when "
                "scalar flag = 1"
            )
        else:
            assert observed is False, (
                f"V2 must stay OFF for class={class_label!r} (only "
                f"{set(DELFIN_MULTI_SIGMA_PATH_V2_DEFAULT_CLASSES)} "
                "should trip the scalar default-on path)"
            )


# ----------------------------------------------------------------------
# Test 4: ``_CLASSES`` env wins over the scalar fallback.
# ----------------------------------------------------------------------
def test_v2_classes_env_takes_precedence() -> None:
    """The whitelist env-var overrides the scalar fallback so an
    operator can re-target the V2 path to *additional* classes (e.g.
    multi_hapto if a future regression appears there)."""
    # Whitelist 'sigma' explicitly — scalar fallback would say no, but
    # _CLASSES MUST win.
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2_CLASSES"] = "sigma"
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2"] = "0"  # scalar OFF
    sigma_mol = _mol(_FIXTURES["sigma"])
    multi_mol = _mol(_FIXTURES["multi_sigma"])
    nm_mol = _mol(_FIXTURES["no_metal"])
    assert _multi_sigma_v2_active(sigma_mol) is True, (
        "DELFIN_MULTI_SIGMA_PATH_V2_CLASSES=sigma must enable V2 for sigma"
    )
    assert _multi_sigma_v2_active(multi_mol) is False, (
        "multi_sigma is NOT in the whitelist → must be OFF"
    )
    assert _multi_sigma_v2_active(nm_mol) is False, (
        "no_metal is NOT in the whitelist → must be OFF"
    )


# ----------------------------------------------------------------------
# Test 5: thread-local override is set + cleaned up.
# ----------------------------------------------------------------------
def test_thread_local_override_cleanup() -> None:
    """The thread-local timeout override is the mechanism that pushes a
    tighter per-seed ``EmbedMultipleConfs`` timeout into
    ``_embed_multiple_confs_with_timeout``.  This test exercises the
    same lifecycle the production call-site does:

      1. start with no override on the calling thread,
      2. set it to a small float (simulating multi_sigma V2 entry),
      3. confirm the global default still applies on threads that
         did not touch the local.

    A regression here would either leak the tight timeout into unrelated
    SMILES on the same thread (causing spurious truncation) or leave
    every embedding running at the global 25 s default after V2 enters
    (defeating the patch).
    """
    # Pre-condition: nothing set on this thread.
    assert not hasattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value")

    # Simulate the V2 entry: write override.
    _MULTIEMBED_TIMEOUT_OVERRIDE.value = 12.0
    assert getattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value", None) == 12.0

    # A second thread MUST NOT see our override (TLS contract).
    import threading

    observed: list = []

    def _peek() -> None:
        observed.append(getattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value", None))

    t = threading.Thread(target=_peek)
    t.start()
    t.join(timeout=2.0)
    assert observed == [None], (
        "thread-local override must NOT bleed into other threads, "
        f"got {observed!r}"
    )

    # Simulate the V2 exit: delete the attribute, confirm cleared.
    delattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value")
    assert not hasattr(_MULTIEMBED_TIMEOUT_OVERRIDE, "value")


# ----------------------------------------------------------------------
# Test 6: malformed env-var values are safe (no crash, V2 stays OFF).
# ----------------------------------------------------------------------
def test_v2_malformed_env_safe() -> None:
    """Whitespace, trailing commas, unknown class labels and empty
    strings must all be tolerated without raising — the V2 path must
    never break the converter on a bad env.

    Under the Phase 1.5 wire-in:
      * empty CLASSES → falls through to the scalar / default-on path
        (multi_sigma class gets V2, every other class does not).
      * non-empty CLASSES with unknown labels → no class matches, so V2
        cannot activate.
    """
    multi_mol = _mol(_FIXTURES["multi_sigma"])
    sigma_mol = _mol(_FIXTURES["sigma"])

    # Empty string (caught explicitly so the helper doesn't take the
    # whitelist branch with an empty set — falls through to the
    # default-on resolver).
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2_CLASSES"] = ""
    assert _multi_sigma_v2_active(multi_mol) is True, (
        "empty CLASSES must fall through to default-on (multi_sigma class)"
    )
    assert _multi_sigma_v2_active(sigma_mol) is False, (
        "empty CLASSES → sigma still OFF-by-default"
    )

    # Trailing comma + whitespace.
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2_CLASSES"] = "  multi_sigma , ,  "
    assert _multi_sigma_v2_active(multi_mol) is True
    assert _multi_sigma_v2_active(sigma_mol) is False

    # Unknown label only — whitelist branch is taken (set is non-empty),
    # but no class matches → V2 stays OFF for every class.
    os.environ["DELFIN_MULTI_SIGMA_PATH_V2_CLASSES"] = "not_a_real_class"
    assert _multi_sigma_v2_active(multi_mol) is False
    assert _multi_sigma_v2_active(sigma_mol) is False
