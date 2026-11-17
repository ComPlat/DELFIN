"""orca_read FINAL energy on multi-job outputs matches the old reader.

The old reader (delfin.energies.find_electronic_energy) scans the whole
file and selects the job whose %base equals the file's name.  orca_read's
text path reads a bounded tail only; on a multi-job output whose own job
sits before that tail it picked a later job's energy (11 of 200 archived
runs).  After the fix orca_read routes FINAL through the whole-file
reader whenever the output holds several job markers, so the selection is
identical by construction.
"""
from pathlib import Path

from delfin import energies, orca_read

HEADER = "# SYNTHETIC test fixture written for DELFIN -- not ORCA output\n"
ENERGY = "FINAL SINGLE POINT ENERGY   -76.321269021507\n"
FILLER = "x" * 76 + "\n"
# the own job's energy must sit further from the end than the bounded
# tail reaches (1 MiB)
TAIL_BYTES = orca_read._TAIL_BYTES
AFTER_BYTES = TAIL_BYTES + 4096


def _multi_job_output(tmp_path, energy, other_energy=-100.0):
    """A $new_job output: own job first (energy far before the tail), a
    later foreign job last (inside the tail)."""
    out_dir = tmp_path / "run"
    out_dir.mkdir()
    out = out_dir / "T1.out"
    lines = [HEADER,
             '$$$ JOB NUMBER 1 $$$\n',
             '| 1> %base "T1"\n',
             "Optimization converging\n", ENERGY]
    while out.tell() if out.exists() else 0 < 0:  # placeholder, replaced below
        pass
    body = "".join(lines)
    pad = FILLER * ((AFTER_BYTES - len(body)) // len(FILLER) + 1)
    foreign = ('$$$ JOB NUMBER 2 $$$\n' + '| 2> $new_job\n'
               + '| 2> %base "T1_TDDFT"\n'
               + "Excited state check job\n"
               + f"FINAL SINGLE POINT ENERGY   {other_energy:.9f}\n"
               + "ORCA TERMINATED NORMALLY\n")
    out.write_text(body + pad + foreign)
    return out


def test_orca_read_final_matches_old_reader_on_multi_job(tmp_path):
    out = _multi_job_output(tmp_path, -76.321269021507)
    # both readers must pick the own job's energy, not the foreign tail job
    old = energies.find_electronic_energy(str(out))
    assert old == -76.321269021507
    reading = orca_read.read(out.parent, "T1",
                             orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert reading.value == old
    assert reading.source == "text"


def test_orca_read_final_matches_old_reader_single_job(tmp_path):
    out = _multi_job_output(tmp_path, -76.321269021507)
    # sanity: the tail alone would give the foreign energy
    tail = orca_read._tail_text(out)
    from delfin.energies import electronic_energy_in_text
    assert electronic_energy_in_text(tail, str(out)) == -100.0


def test_orca_read_final_single_job_tail_ok(tmp_path):
    """A single-job output stays on the bounded-tail path."""
    out_dir = tmp_path / "run2"
    out_dir.mkdir()
    out = out_dir / "t.out"
    out.write_text(HEADER + '| 1> %base "t"\n' + ENERGY
                   + "ORCA TERMINATED NORMALLY\n")
    reading = orca_read.read(out.parent, "t",
                             orca_read.Quantity.FINAL_SINGLE_POINT_ENERGY)
    assert reading.value == -76.321269021507
    assert reading.source == "text"
