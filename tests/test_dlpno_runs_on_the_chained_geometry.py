"""The DLPNO step runs on the geometry the chain hands it.

``delfin.tools.templates._opi_chain`` composes

    smiles_to_xyz -> orca_opt_opi -> orca_freq_opi -> orca_dlpno_opi

and ``Pipeline.run`` passes each step the PREVIOUS step's geometry as the
``geometry=`` argument (``current_geom``), never as ``elements`` /
``coordinates`` kwargs.  ``OrcaDlpnoOpiAdapter`` declares
``consumes = ("geometry",)`` and accepts the ``geometry=`` parameter -- and
then never reads it: ``execute`` looked only at ``kwargs["elements"]`` and
``kwargs["coordinates"]``, so the LAST step of DELFIN's own OPI chain returned
FAILED with "missing required 'elements'" before writing an input.

This is the same defect in the last step that a real run found in the second
step (``tests/test_opi_opt_runs_on_the_chained_geometry.py``); the shape of
the chain means a chain can only run if every step consumes what the step
before it produces.

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


def _written_input(monkeypatch, tmp_path, step="orca_dlpno_opi", **kwargs) -> str:
    """The input text the step wrote, with ORCA stubbed out."""
    import delfin.orca as orca_mod

    seen: dict[str, str] = {}

    def fake_run_orca(inp, out, **kw):  # noqa: ANN001, ANN003
        seen["input"] = open(inp, encoding="utf-8").read()
        return False

    monkeypatch.setattr(orca_mod, "run_orca", fake_run_orca)
    run_step(step, work_dir=tmp_path, cores=1, **kwargs)
    assert "input" in seen, "the step wrote no input at all"
    return seen["input"]


def test_a_chained_geometry_is_enough_to_write_the_input(monkeypatch, tmp_path):
    """The shape the template uses: geometry= only, from the step before."""
    xyz = tmp_path / "input.xyz"
    xyz.write_text(H2O_XYZ)

    text = _written_input(
        monkeypatch, tmp_path / "run",
        geometry=xyz, charge=0, basis="def2-TZVPP")

    assert text.strip(), "the step wrote an empty input"
    assert "dlpno" in text.lower(), f"not a DLPNO input:\n{text}"


def test_the_chained_geometry_carries_the_atoms(monkeypatch, tmp_path):
    """The structure in the input is the one the chain handed over."""
    xyz = tmp_path / "input.xyz"
    xyz.write_text(H2O_XYZ)

    text = _written_input(
        monkeypatch, tmp_path / "run",
        geometry=xyz, charge=0, basis="def2-TZVPP")

    assert "0.757" in text, (
        "the chained geometry's coordinates are not in the input:\n" + text)
    assert "O" in text and "H" in text, text


def test_explicit_elements_still_work(monkeypatch, tmp_path):
    """Control: a caller passing elements/coordinates keeps working.

    Green on both trees -- it is the chained-geometry path that was broken.
    """
    text = _written_input(
        monkeypatch, tmp_path / "run",
        elements=["O", "H", "H"],
        coordinates=[[0.0, 0.0, 0.0], [0.0, 0.757, 0.587], [0.0, -0.757, 0.587]],
        charge=0, basis="def2-TZVPP")

    assert "0.757" in text, text
