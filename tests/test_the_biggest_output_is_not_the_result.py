"""An abandoned run is often the biggest file in the folder.

``delfin.api`` reads a folder's result through one picker, and eight
public functions in that module are built on it: orbital energies,
excited states, thermochemistry and the rest. It took the LARGEST .out.

Size picks the wrong file whenever a job was abandoned and rerun, and
abandoned runs tend to be bigger: they ran for a long time before being
killed and wrote every SCF cycle while doing so.

Measured on a real workspace (bench17 t13, 2026-11-23):

    tddft_pbe0.out       98,374 bytes   abandoned mid-SCF, no spectrum
    tddft_pbe0_v2.out    86,460 bytes   the result, with the spectrum

So the picker chose the failure, and every reader on top of it reported
"no absorption block (run with !TDDFT or !CIS)" while the answer sat
beside it in the same directory. After the fix the same folder reads 10
transitions with a first bright state at 4.194 eV.

A FINISHED OUTPUT IS PREFERRED, and the largest of the finished ones is
taken. Where none finished, the largest is still taken: a folder of
unfinished work reads exactly as it did before, so this prefers a result
where one exists and never refuses to read.

The outputs here are written by this test; the repository may not hold
real engine output (LICENSE_GUARD.md).
"""
from __future__ import annotations

import pytest

from delfin import api

_SPECTRUM = """\
-----------------------------------------------------------------------
         ABSORPTION SPECTRUM VIA TRANSITION ELECTRIC DIPOLE MOMENTS
-----------------------------------------------------------------------
      Transition      Energy     Energy  Wavelength     fosc        D2
                       (eV)      (cm-1)      (nm)                (au**2)
-----------------------------------------------------------------------
  0-1A  ->  1-1A     4.194000   33827.0     295.5   0.369000000   3.69000
"""

_FINISHED = _SPECTRUM + "                   ****ORCA TERMINATED NORMALLY****\n"


def _abandoned(size: int) -> str:
    """A long unfinished run: SCF cycles and no final block."""
    return ("           Program Version 6.1.0  -  RELEASE  -\n"
            + "*** Restarting incremental Fock matrix formation ***\n"
              * size)


def test_a_finished_output_wins_over_a_bigger_abandoned_one(tmp_path):
    """The measured case, as small as it goes."""
    (tmp_path / "job.out").write_text(_abandoned(4000), encoding="utf-8")
    (tmp_path / "job_v2.out").write_text(_FINISHED, encoding="utf-8")
    abandoned = (tmp_path / "job.out").stat().st_size
    result = (tmp_path / "job_v2.out").stat().st_size
    assert abandoned > result, (
        "the fixture must have the abandoned file bigger, or it tests "
        "nothing")
    target, text, err = api._read_largest_out(str(tmp_path))
    assert err is None
    assert target.name == "job_v2.out"
    assert "ABSORPTION SPECTRUM" in text


def test_the_reader_on_top_of_it_now_finds_the_spectrum(tmp_path):
    """What the defect cost: eight public functions read the failure."""
    (tmp_path / "job.out").write_text(_abandoned(4000), encoding="utf-8")
    (tmp_path / "job_v2.out").write_text(_FINISHED, encoding="utf-8")
    res = api.extract_excited_states(str(tmp_path))
    assert res.error is None, res.error
    assert len(res.transitions) == 1
    assert res.first_bright is not None
    assert res.first_bright["energy_ev"] == pytest.approx(4.194, abs=1e-3)


def test_the_largest_finished_output_is_taken(tmp_path):
    """Among finished outputs, size still decides: a short restart log
    that happens to end normally must not beat the full run."""
    (tmp_path / "short.out").write_text(_FINISHED, encoding="utf-8")
    (tmp_path / "full.out").write_text(
        _FINISHED + "x" * 50000, encoding="utf-8")
    target, _text, _err = api._read_largest_out(str(tmp_path))
    assert target.name == "full.out"


def test_a_folder_where_nothing_finished_reads_as_before(tmp_path):
    """No regression: unfinished work is still read, largest first, so a
    caller inspecting a crashed run sees the longest log."""
    (tmp_path / "a.out").write_text(_abandoned(10), encoding="utf-8")
    (tmp_path / "b.out").write_text(_abandoned(500), encoding="utf-8")
    target, text, err = api._read_largest_out(str(tmp_path))
    assert err is None
    assert target.name == "b.out"
    assert "Restarting" in text


def test_a_single_output_is_read_whatever_its_state(tmp_path):
    (tmp_path / "only.out").write_text(_abandoned(3), encoding="utf-8")
    target, _text, err = api._read_largest_out(str(tmp_path))
    assert err is None and target.name == "only.out"


def test_a_file_path_is_read_directly(tmp_path):
    """The picker also takes a file rather than a folder."""
    f = tmp_path / "one.out"
    f.write_text(_FINISHED, encoding="utf-8")
    target, text, err = api._read_largest_out(str(f))
    assert err is None and target == f and "ABSORPTION" in text


def test_a_missing_folder_and_an_empty_one_still_report(tmp_path):
    _t, _x, err = api._read_largest_out(str(tmp_path / "nowhere"))
    assert err == "folder missing"
    empty = tmp_path / "empty"
    empty.mkdir()
    _t, _x, err = api._read_largest_out(str(empty))
    assert err == "no ORCA output"


def test_an_unreadable_output_is_reported_not_skipped(tmp_path, monkeypatch):
    """A read that fails must name the file rather than look like an
    empty result."""
    f = tmp_path / "x.out"
    f.write_text(_FINISHED, encoding="utf-8")

    real = type(f).read_text

    def boom(self, *a, **k):
        if self.name == "x.out" and k.get("errors") == "replace" \
                and a == ():
            raise OSError("denied")
        return real(self, *a, **k)

    monkeypatch.setattr(type(f), "read_text", boom)
    target, text, err = api._read_largest_out(str(tmp_path))
    assert target is not None and text == ""
    assert err and "read failed" in err
