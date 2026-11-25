"""The method a caller asks for is the method that runs.

``orca_opt_opi`` is reached by the ``opi_opt_freq_dlpno`` template and by the
agent's ``pipeline-build`` skill, both of which name the method with the
template vocabulary ``dft_method`` / ``dft_basis`` (the template documents
``dft_method="r2SCAN-3c", dft_basis="def2-SVP"`` and every step of the chain
is handed those two keys).  The adapter reads ``functional`` / ``basis`` and
silently falls back to its own defaults for anything else, so the request
"r2SCAN-3c opt+freq" ran BP86/def2-SVP while the plan, the cost estimate and
the user all said r2SCAN-3c.

Red control (measured on the baseline): the ``!`` line of the input the step
writes is ``!bp86`` for a caller asking PBE0.  A silently substituted
functional is the failure mode CONTRIBUTING.md opens with — a number that
looks reasonable and is wrong — so the test asserts the written input, not a
return value: the input is what ran.

The control case (an explicit ``functional=``/``basis=`` pair) passes on both
trees, so a green here is not just "some keyword appears".
"""

from __future__ import annotations

from pathlib import Path

import pytest

# Same convention as the rest of tests/test_opi_*: the main CI jobs have no
# OPI and skip, the OPI job sets DELFIN_OPI_REQUIRED=1 which turns the skip
# into a failure, so this test cannot pass by not being run.
pytest.importorskip("opi")

from delfin.tools import run_step  # noqa: E402

WATER_E = ["O", "H", "H"]
WATER_C = [
    [0.0, 0.0, 0.0],
    [0.0, 0.757, 0.587],
    [0.0, -0.757, 0.587],
]


def _bang_lines(monkeypatch, tmp_path, **kwargs) -> list[str]:
    """The `!` keywords of the input the step writes, with the launch stubbed.

    ``run_orca`` is replaced (the building block imports it inside its run
    path), so no binary is started; the step returns FAILED at the read-back
    because no output exists, which is irrelevant here -- the input is already
    on disk and is what the assertion is about.
    """
    import delfin.orca as orca_mod

    seen: dict[str, str] = {}

    def fake_run_orca(inp, out, **kw):  # noqa: ANN001, ANN003 - signature of run_orca
        seen["input"] = Path(inp).read_text()
        return False

    monkeypatch.setattr(orca_mod, "run_orca", fake_run_orca)
    run_step(
        "orca_opt_opi",
        elements=WATER_E,
        coordinates=WATER_C,
        charge=0,
        work_dir=tmp_path,
        cores=2,
        **kwargs,
    )
    assert "input" in seen, "the step wrote no input at all"
    return [ln.strip().lower() for ln in seen["input"].splitlines()
            if ln.strip().startswith("!")]


def test_the_template_vocabulary_selects_the_method(monkeypatch, tmp_path):
    """A caller using the template's own keys gets the method it named."""
    bang = _bang_lines(monkeypatch, tmp_path,
                       dft_method="PBE0", dft_basis="def2-TZVP")
    assert any("pbe0" in ln for ln in bang), (
        "the caller asked for PBE0 via the template's dft_method/dft_basis keys; "
        f"the step wrote {bang}")
    assert any("def2-tzvp" in ln for ln in bang), (
        f"the caller asked for def2-TZVP; the step wrote {bang}")


def test_an_explicit_functional_is_still_honoured(monkeypatch, tmp_path):
    """Control: the adapter's own parameter names keep working (green on both
    trees -- it is the alias path above that was broken)."""
    bang = _bang_lines(monkeypatch, tmp_path,
                       functional="B3LYP", basis="def2-SVP")
    assert any("b3lyp" in ln for ln in bang), bang
