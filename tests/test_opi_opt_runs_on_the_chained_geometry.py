"""The OPT step runs on the geometry the chain hands it.

``delfin.tools.templates._opi_chain`` composes

    smiles_to_xyz -> orca_opt_opi -> orca_freq_opi -> orca_dlpno_opi

and ``Pipeline.run`` passes each step the PREVIOUS step's geometry as the
``geometry=`` argument (``current_geom``), never as ``elements`` /
``coordinates`` kwargs.  ``OrcaOptOpiAdapter.execute`` read only the kwargs,
so the first step of DELFIN's own OPI chain -- and of every spec the agent's
workflow builder emits for it -- failed with "missing required 'elements'"
before it wrote an input.  The sibling ``orca_opt_freq_opi`` (same template,
same chain) always read the chained geometry, which is why the chain looked
half-working.

Reachability (rule 7): the product path is
``PipelineTemplate("opi_opt_freq_dlpno").build(...).run(geometry=...)`` ->
``Pipeline.run`` -> ``run_step`` -> this adapter.  ``run_step`` is the
registry dispatch the pipeline and the CLI both call, so these tests drive it
rather than calling ``execute`` directly.

No ORCA binary: DELFIN's launcher is replaced, and the assertion is on the
input file the step wrote before launching anything.
"""

from __future__ import annotations

import pytest

pytest.importorskip("opi")

from delfin.tools import run_step  # noqa: E402

H2O_XYZ = "3\nwater\nO 0.0 0.0 0.0\nH 0.0 0.757 0.587\nH 0.0 -0.757 0.587\n"


def _written_bang(monkeypatch, tmp_path, **kwargs) -> list[str]:
    """The ``!`` keywords of the input the step wrote, with ORCA stubbed out."""
    import delfin.orca as orca_mod

    seen: dict[str, str] = {}

    def fake_run_orca(inp, out, **kw):  # noqa: ANN001, ANN003
        seen["input"] = open(inp, encoding="utf-8").read()
        return False

    monkeypatch.setattr(orca_mod, "run_orca", fake_run_orca)
    run_step("orca_opt_opi", work_dir=tmp_path, cores=1, **kwargs)
    assert "input" in seen, "the step wrote no input at all"
    return [ln.strip().lower() for ln in seen["input"].splitlines()
            if ln.strip().startswith("!")]


def test_a_chained_geometry_is_enough_to_write_the_input(monkeypatch, tmp_path):
    """The shape the template uses: geometry= only, from the step before."""
    xyz = tmp_path / "input.xyz"
    xyz.write_text(H2O_XYZ)

    bang = _written_bang(
        monkeypatch, tmp_path / "run",
        geometry=xyz, charge=0, dft_method="PBE0", dft_basis="def2-TZVP")

    assert any("opt" in ln for ln in bang), (
        f"the step wrote no optimisation keyword: {bang}")


def test_the_chained_geometry_carries_the_atoms(monkeypatch, tmp_path):
    """The structure in the input is the one the chain handed over."""
    xyz = tmp_path / "input.xyz"
    xyz.write_text(H2O_XYZ)

    import delfin.orca as orca_mod

    seen: dict[str, str] = {}

    def fake_run_orca(inp, out, **kw):  # noqa: ANN001, ANN003
        seen["input"] = open(inp, encoding="utf-8").read()
        return False

    monkeypatch.setattr(orca_mod, "run_orca", fake_run_orca)
    run_step("orca_opt_opi", geometry=xyz, work_dir=tmp_path / "run",
             cores=1, charge=0, dft_method="PBE0", dft_basis="def2-TZVP")

    text = seen.get("input", "")
    assert "0.757" in text, (
        "the chained geometry's coordinates are not in the input:\n" + text)
    assert "O" in text and "H" in text, text


def test_explicit_elements_still_work(monkeypatch, tmp_path):
    """Control: a caller passing elements/coordinates keeps working.

    Green on both trees -- it is the chained-geometry path that was broken.
    """
    bang = _written_bang(
        monkeypatch, tmp_path / "run",
        elements=["O", "H", "H"],
        coordinates=[[0.0, 0.0, 0.0], [0.0, 0.757, 0.587], [0.0, -0.757, 0.587]],
        charge=0, functional="B3LYP", basis="def2-SVP")
    assert any("b3lyp" in ln for ln in bang), bang
