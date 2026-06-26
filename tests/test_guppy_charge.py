"""Charge-from-SMILES derivation tests for GUPPY.

Covers: cations, anions, neutral complexes, tokenizer fallback, and
enforcement that an explicit --charge override never overrules the derived
charge in run_sampling().
"""
import logging
from pathlib import Path

import pytest

from delfin import guppy_sampling


@pytest.mark.parametrize(
    "smiles, expected",
    [
        # [Fe(bpy)3]2+ style: Fe(II) + 3 neutral bipyridines (each bpy written as 2 py).
        ("[Fe+2].N1=CC=CC=C1.N1=CC=CC=C1.N1=CC=CC=C1.N1=CC=CC=C1.N1=CC=CC=C1.N1=CC=CC=C1", 2),
        # [Cu(NH3)4]2+
        ("[Cu+2].N.N.N.N", 2),
        # Cisplatin Pt(NH3)2Cl2
        ("[Pt+2].[Cl-].[Cl-].N.N", 0),
        # [Fe(CN)6]4-
        ("[Fe-4].[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N.[C-]#N", -10),
        # Methylammonium
        ("C[NH3+]", 1),
        # Acetate
        ("CC(=O)[O-]", -1),
        # Neutral water
        ("O", 0),
        # Double-charge token '++'
        ("[Cu++]", 2),
        # Double-charge token '--'
        ("[S--]", -2),
    ],
)
def test_derive_charge_from_smiles(smiles: str, expected: int) -> None:
    assert guppy_sampling._derive_charge_from_smiles(smiles) == expected


def test_derive_charge_logs_per_bracket_contributions(caplog) -> None:
    caplog.set_level(logging.INFO)
    total = guppy_sampling._derive_charge_from_smiles("[Fe+2].[Cl-].[Cl-]")
    assert total == 0
    messages = "\n".join(rec.getMessage() for rec in caplog.records)
    assert "SMILES formal charges" in messages
    assert "total +0" in messages or "total 0" in messages


def test_derive_charge_token_fallback_when_rdkit_missing(monkeypatch) -> None:
    """Regex fallback must handle +n, -n, ++, -- even without RDKit."""
    monkeypatch.setattr(guppy_sampling, "RDKIT_AVAILABLE", False)
    assert guppy_sampling._derive_charge_from_smiles("[Fe+2]") == 2
    assert guppy_sampling._derive_charge_from_smiles("[Cu++]") == 2
    assert guppy_sampling._derive_charge_from_smiles("[Fe-4]") == -4
    assert guppy_sampling._derive_charge_from_smiles("[S--]") == -2
    assert guppy_sampling._derive_charge_from_smiles("O") == 0


def test_run_sampling_warns_when_override_disagrees(monkeypatch, tmp_path, caplog) -> None:
    """Verify the derived charge wins over a disagreeing --charge override."""
    input_file = tmp_path / "input.txt"
    input_file.write_text("[Fe+2].N.N.N.N\n", encoding="utf-8")

    # Short-circuit the heavy sampling path: force `_collect_start_geometries`
    # to return an empty list so `run_sampling` bails with code 1 before doing
    # any XTB work. We only care that the derived-vs-override warning fires
    # during the preamble (l.774-781 in guppy_sampling.py).
    monkeypatch.setattr(
        guppy_sampling, "_collect_start_geometries", lambda *args, **kwargs: []
    )

    caplog.set_level(logging.WARNING)
    rc = guppy_sampling.run_sampling(
            input_file=input_file,
            runs=1,
            charge=-99,  # deliberately wrong
            pal=1,
            maxcore=100,
            parallel_jobs=1,
            method="XTB2",
            output_file=tmp_path / "GUPPY_try.xyz",
            workdir=tmp_path / "GUPPY",
            seed=31,
            allow_partial=True,
        )

    assert rc == 1  # empty starts -> non-zero return
    messages = "\n".join(rec.getMessage() for rec in caplog.records)
    assert "Ignoring provided charge override" in messages
    assert "-99" in messages
