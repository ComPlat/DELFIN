"""The DLPNO step's MOREAD path builds an input, it does not raise.

The ``opi_opt_freq_dlpno`` template promises the DLPNO single point is
"started from the chain's orbitals (MOREAD)", and the engine's ``gbw ->
moread`` wiring (``delfin.tools._wiring``) supplies the upstream ``.gbw`` as
the ``moread`` kwarg.  ``tests/test_dlpno_step_consumes_moread_wiring.py``
pins that the kwarg is *read* -- by inspecting signatures and source.  A
source inspection cannot see that the code it finds raises the moment it
runs: ``Input.add_arbitrary_string`` is keyword-only after the string
(``opi`` 2.0.0: ``(self, string, /, *, pos=None)``), so passing the position
positionally is a ``TypeError``, and the first step of the chain that really
carries orbitals forward fails instead of starting from them.

These tests exercise the BUILD path, so they need no ORCA binary: they
assert the file ``build_dlpno_input`` writes.  They run where ``opi`` is
installed, and skip (with the import named) everywhere else.

Their reachability is narrower than the ``test_opi_*`` files': the CI "OPI
tests" job installs ``opi`` and passes ``tests/test_opi_*.py
tests/test_orca_knowledge*.py``, and this filename is not in that list, so
that job does not run these two ``test_dlpno_*`` files.  ``DELFIN_OPI_REQUIRED``
therefore does not reach this skip, and no CI job does: they run for a
developer with ``opi`` installed, and that is all.  Widening the job's path
list would change CI, which was not verified from here (``opi`` is not
installed on the reviewing host, so the newly selected tests were never
executed), and a test that is newly selected by a job nobody ran is the same
claim in a new place.

Manual ground truth for the written form (ORCA 6.1 manual, section 2.20.5
"MORead"/Table 2.68): MORead reads the MOs of a previous calculation from the
file named by ``%moinp`` in the simple-input line's block section, and
``AutoStart`` only picks up files of the same base name -- which a chained
step's upstream ``.gbw`` is not.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("opi")

from delfin.tools.adapters.orca_dlpno_opi import build_dlpno_input  # noqa: E402

H2O_E = ["O", "H", "H"]
H2O_C = [
    [0.0, 0.0, 0.0],
    [0.0, 0.757, 0.587],
    [0.0, -0.757, 0.587],
]


def _build(tmp_path: Path, **kwargs) -> str:
    written = build_dlpno_input(
        tmp_path,
        elements=H2O_E,
        coordinates=H2O_C,
        charge=0,
        mult=1,
        basis="def2-TZVPP",
        ncores=2,
        maxcore=1000,
        **kwargs,
    )
    return Path(written).read_text()


def test_moread_writes_the_manual_form(tmp_path):
    """The wired .gbw is copied beside the input and named by %moinp."""
    upstream = tmp_path / "upstream"
    upstream.mkdir()
    gbw = upstream / "dft_opt.gbw"
    gbw.write_bytes(b"gbw-bytes")

    text = _build(tmp_path, moread=str(gbw))
    lowered = text.lower()

    assert "moread" in lowered, (
        "the step was given start orbitals but its input never asks ORCA to "
        f"read them:\n{text}")
    assert '%moinp "dlpno.gbw"' in lowered, (
        "MOREAD without %moinp names no file to read; the manual's form is "
        f"`! moread` plus `%moinp \"<file>\"`:\n{text}")
    assert (tmp_path / "dlpno.gbw").is_file(), (
        "the .gbw was not copied next to the input, so ORCA cannot find it "
        "by the name %moinp gives")


def test_no_moread_leaves_the_input_alone(tmp_path):
    """Control: without the wiring the input names no start orbitals.

    Green on both trees -- ``moread=None`` is the not-chained case, and it
    must stay free of a dangling MORead reference.
    """
    text = _build(tmp_path)
    assert "moinp" not in text.lower(), text
    assert not (tmp_path / "dlpno.gbw").exists()


def test_a_dangling_moread_is_refused_not_written(tmp_path):
    """A .gbw path that does not exist is an error, not a dead guess."""
    with pytest.raises(FileNotFoundError):
        _build(tmp_path, moread=str(tmp_path / "missing.gbw"))
    assert not (tmp_path / "dlpno.inp").exists(), (
        "an input naming a nonexistent .gbw was written anyway; ORCA would "
        "then abort at guess time instead of the step saying what is wrong")
