"""A non-chemistry task must never classify as chemistry.

Measured 2026-11-01 (delfin.agent.briefing.classify_task): "add a SIMD
kernel to the Rust encoder crate" scores "chemistry" — the plain
substring check fires on "md " sitting inside "simd kernel" (the trigger
is "md " + space, any word ending in md followed by a space). A
chemistry misfire injects the "CHEMISTRY TASK DETECTED" protocol into a
non-chemistry prompt. This test pins the negative cases.
"""
from __future__ import annotations

import pytest

from delfin.agent.briefing import classify_task

# Tasks a chemist is NOT asking for, phrased the way an engineer does.
NOT_CHEMISTRY = (
    "add a SIMD kernel to the Rust encoder crate",
    "refactor the charge-coupled device driver module",
    "fix the scan method of the barcode parser",
    "improve the caching method in the api layer",
    "audit the energy accounting in the billing service",
    "add regression tests for the import pipeline",
)


@pytest.mark.parametrize("text", NOT_CHEMISTRY,
                         ids=[t[:40] for t in NOT_CHEMISTRY])
def test_engineering_tasks_are_not_chemistry(text):
    assert classify_task(text) != "chemistry", (
        f"{text!r} classified as chemistry — the chemistry protocol "
        "would leak into a non-chemistry prompt")


@pytest.mark.parametrize("text", [
    "optimize the DFT basis set for iron",
    "compute the redox potential with ORCA",
    "which spin state is lowest for my Fe complex",
    # Weak cues (energy, scan, charge, method, md) count only together
    # with a molecular cue or a strong keyword -- these plain tasks use
    # only a weak cue plus a molecular context:
    "calculate the energy of benzene",
    "relaxed scan of the C-C-C-C dihedral of butane",
    "what is the partial charge on the nitrogen",
    "run an md of water in a box",
    "which method for a Cu(II) complex",
])
def test_real_chemistry_tasks_still_classify_as_chemistry(text):
    assert classify_task(text) == "chemistry"
