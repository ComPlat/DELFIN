"""A thermochemistry block more than one tail-length from the end is still read.

``orca_read.read`` reads the last ``_TAIL_BYTES`` (1 MiB) for the Gibbs, ZPE,
enthalpy, inner-energy and entropy quantities.  The final single-point energy
is deliberately exempt -- it is answered by ``energies.find_electronic_energy``
on the whole file, because the bounded tail was measured to pick a foreign job
on multi-job outputs (``orca_read._final_energy_text``).

The thermochemistry path has the same exposure and no such exemption.  DELFIN
itself writes multi-job outputs: ``T1.inp`` optimises the triplet under
``%base "T1"`` and appends a closed-shell TD-DFT check job under
``%base "T1_TDDFT"`` (``delfin/energies.py``, :func:`own_job_text`).  When the
appended job prints more than a tail's worth of output, the thermochemistry
block of the job that owns the file falls outside the window, and the reader
answers ``None`` for a Gibbs energy the file contains.

These tests are red before the fix: the window-only read returns ``None``
where ``energies.find_gibbs_energy`` -- the reader this module routes to and
the one callers had before -- returns the value.  The control cases (a block
inside the window, and a file with no block at all) pass before and after, so
the regression they guard is not confused with the fix.
"""
from pathlib import Path

import pytest

from delfin import orca_read
from delfin.energies import find_gibbs_energy, find_ZPE
from delfin.orca_read import Quantity, read

HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
GIBBS_LINE = "Final Gibbs free energy         ...  -1000.12345678 Eh\n"
ZPE_LINE = "Zero point energy               ...      0.02120238 Eh\n"

#: One output line of a second job; the filler below is built from it.
FILLER_LINE = "TD-DFT check job output line that is not thermochemistry\n"


def _block_then_long_job(tmp_path: Path, name: str = "chain") -> Path:
    """A job-1 thermochemistry block, then > 1 MiB of appended job-2 output."""
    filler = FILLER_LINE * 22000
    text = (
        HEADER
        + "Program Version 6.1.1 - RELEASE\n"
        + 'JOB NUMBER 1\n%base "chain"\n'
        + ZPE_LINE
        + GIBBS_LINE
        + "****ORCA TERMINATED NORMALLY****\n"
        + 'JOB NUMBER 2\n%base "chain_TDDFT"\n'
        + filler
        + "****ORCA TERMINATED NORMALLY****\n"
    )
    out = tmp_path / f"{name}.out"
    out.write_text(text, encoding="utf-8")
    # The test only means something if the block really is outside the window.
    assert out.stat().st_size > orca_read._TAIL_BYTES
    head = out.read_text(encoding="utf-8")
    assert len(head.encode()) - head.encode().rfind(GIBBS_LINE.encode()) > orca_read._TAIL_BYTES
    return out


def _block_inside_the_window(tmp_path: Path, name: str = "short") -> Path:
    """Control: a thermochemistry block near the end of the file."""
    out = tmp_path / f"{name}.out"
    out.write_text(
        HEADER
        + "Program Version 6.1.1 - RELEASE\n"
        + ZPE_LINE
        + GIBBS_LINE
        + "****ORCA TERMINATED NORMALLY****\n",
        encoding="utf-8",
    )
    assert out.stat().st_size < orca_read._TAIL_BYTES
    return out


@pytest.mark.parametrize("quantity,expected_reader", [
    (Quantity.GIBBS, find_gibbs_energy),
    (Quantity.ZPE, find_ZPE),
])
def test_a_block_beyond_the_tail_is_read_whole_file(tmp_path, quantity,
                                                    expected_reader):
    """The value the file contains, not None, when it sits past the window."""
    out = _block_then_long_job(tmp_path)
    expected = expected_reader(str(out))
    assert expected is not None, "fixture must carry the quantity"

    reading = read(tmp_path, out.stem, quantity)

    assert reading.value == expected, (
        f"{quantity.name} sits more than {orca_read._TAIL_BYTES} bytes from "
        f"the end of the output; the tail-only read answered "
        f"{reading.value!r} where the file contains {expected!r}"
    )


def test_the_block_inside_the_window_is_unchanged(tmp_path):
    """Control that passes before and after: the ordinary case."""
    out = _block_inside_the_window(tmp_path)
    reading = read(tmp_path, out.stem, Quantity.GIBBS)
    assert reading.value == find_gibbs_energy(str(out))
    assert reading.source == "text"


def test_an_output_without_a_block_stays_none(tmp_path):
    """Control that passes before and after: absence is still answered None."""
    out = tmp_path / "nothing.out"
    out.write_text(HEADER + "Program Version 6.1.1 - RELEASE\n", encoding="utf-8")
    reading = read(tmp_path, out.stem, Quantity.GIBBS)
    assert reading.value is None
    assert reading.source == "none"


def test_the_wider_read_is_paid_only_when_the_window_is_empty(tmp_path):
    """The ordinary run keeps the bounded read; only an empty window widens it.

    The whole-file fallback exists for the case above. It must not become the
    normal path: a 12 MB output read whole cost ~200 MB peak, which is what
    the bounded read was introduced for.
    """
    import tracemalloc

    budget = 8 * 1024 * 1024

    short = _block_inside_the_window(tmp_path)
    tracemalloc.start()
    read(tmp_path, short.stem, Quantity.GIBBS)
    _, peak = tracemalloc.get_traced_memory()
    tracemalloc.stop()
    assert peak < budget, (
        f"read {peak / 1024 / 1024:.1f} MB for a thermochemistry block inside "
        "the window -- the whole file is being read where the tail suffices"
    )

    # And with the block outside the window, the value is returned rather
    # than the read being abandoned.
    far = _block_then_long_job(tmp_path, name="far")
    assert read(tmp_path, far.stem, Quantity.GIBBS).value == find_gibbs_energy(str(far))
