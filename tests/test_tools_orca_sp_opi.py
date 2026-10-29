"""Tests for the OPI-backed single-point DFT building block.

The building block (``delfin.tools.adapters.orca_sp_opi``) is DELFIN's first
ORCA step built on OPI (orca-pi 2.0.0): every keyword and block option is
checked against DELFIN's ORCA knowledge store before the input is written
(``delfin.orca_knowledge``), the input is assembled with OPI's input classes
rather than string templates, the run is launched through DELFIN's existing
``delfin.orca.run_orca`` (never a new runner), and the result is read back
with OPI while honouring the rule that no value is reported until the
convergence flags say it is meaningful.

These fixtures were measured once against ORCA 6.1.1 runs kept under
``tests/fixtures/orca611/`` (H2O opt+freq, O2 triplet, SCF-not-converged).
A fake ORCA script copies those runs into the working directory, so the tests
never invoke the real ORCA binary.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters.orca_sp_opi import (  # noqa: E402
    validate_method,
)


@pytest.fixture(scope="module")
def valid_method() -> dict:
    """A method the knowledge store accepts (B3LYP/def2-SVP/D4/RIJCOSX+def2/J)."""
    return {
        "functional": "B3LYP",
        "basis": "def2-SVP",
        "dispersion": "D4",
        "ri": "RIJCOSX",
        "aux_basis": "def2/J",
    }


class TestValidateMethod:
    """Every keyword is checked against delfin.orca_knowledge before use.
    An unknown keyword or block option is an error that stops the step; the
    knowledge-store cards' findings come back as warnings."""

    def test_known_method_produces_no_errors(self, valid_method):
        errors, warnings = validate_method(**valid_method)
        assert errors == []
        assert warnings == []

    def test_unknown_functional_is_an_error(self):
        errors, warnings = validate_method(functional="NOTAFUNCTIONAL", basis="def2-SVP")
        assert any("NOTAFUNCTIONAL" in e for e in errors)
        assert warnings == []

    def test_unknown_basis_is_an_error(self):
        errors, _ = validate_method(functional="B3LYP", basis="not-a-basis")
        assert any("not-a-basis" in e for e in errors)

    def test_unknown_dispersion_is_an_error(self):
        errors, _ = validate_method(functional="B3LYP", basis="def2-SVP", dispersion="D999")
        assert any("D999" in e for e in errors)

    def test_unknown_ri_is_an_error(self):
        errors, _ = validate_method(functional="B3LYP", basis="def2-SVP", ri="RI-NOTREAL")
        assert any("RI" in e for e in errors)

    def test_unknown_scf_option_is_an_error(self):
        errors, _ = validate_method(
            functional="B3LYP",
            basis="def2-SVP",
            scf_options={"maxiter_wrong": 10},
        )
        assert any("scf" in e for e in errors)

    def test_knowledge_card_warning_is_reported_not_silenced(self):
        """A finding from a verified knowledge card comes back as a warning,
        never silently dropped."""
        card = {
            "id": "test.card-warns",
            "topic": "scf",
            "kind": "recommendation",
            "statement": "A TIGHTSCF calc must set %scf maxiter (test card, not shipped).",
            "orca": ">=6.0",
            "status": "verified",
            "severity": "warning",
            "applies": {"keyword": "TIGHTSCF"},
            "check": {"requires_keyword": ["MAXITER"]},
            "sources": [{"type": "manual", "version": "6.1", "page": "x"}],
        }
        errors, warnings = validate_method(
            functional="B3LYP",
            basis="def2-SVP",
            extra_keywords=["TIGHTSCF"],
            cards=[card],
        )
        assert errors == []
        assert any(w.get("id") == "test.card-warns" for w in warnings)

