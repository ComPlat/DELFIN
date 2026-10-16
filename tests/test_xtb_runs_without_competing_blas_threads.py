"""The dashboard's xtb runs give BLAS one thread; OpenMP gets the cores.

Input: an xtb stand-in on PATH that records the OPENBLAS_NUM_THREADS it was
started with, then each dashboard path that starts xtb for a structure on
screen -- an optimisation (gfn_optimize.optimize_with_gfn), the climb's
gradient by the command line (climb._CommandLine) and the climb's Hessian
(climb.Climb.modes_from_xtb).
Output: every recorded value is "1".

The paths set OMP_NUM_THREADS and MKL_NUM_THREADS to the cores they want and
left OpenBLAS alone. The conda-forge xtb -- the build DELFIN's installer
fetches and the one CI runs -- is linked against OpenBLAS, which then starts
as many threads again as the CPUs it may use, on top of xtb's OpenMP team.
Measured 2026-10-16, conda-forge xtb 6.7.1 gfortran_ha02e54d_5, GFN2,
relax_steps with 20 cycles, median of three, wall time in ms:

    CPUs   molecule (cores asked)   OpenBLAS unset   OPENBLAS_NUM_THREADS=1
       4   glycine (4)                       323                       34
       4   ibuprofen (8)                     421                      216
       4   eicosane (16)                     670                      319
       8   ibuprofen (8)                    9277                      269
       8   eicosane (16)                    6991                      368
     384   glycine (4)                        66                       68
     384   eicosane (16)                     320                      270

With CPUs to spare the setting changes nothing; with as many CPUs as threads
asked -- a four-core laptop, a SLURM allocation, a CI runner -- unset was 2x
to 34x slower. The xtb release binary (static, MKL) was not affected.

Rejected: OPENBLAS_NUM_THREADS equal to the cores. Measured on four CPUs it
is what unset already does (270 ms for glycine); the two thread pools then
compete for the same CPUs.
"""
from __future__ import annotations

import os
import stat

import pytest


_XYZ = "3\nwater\nO 0 0 0\nH 0.96 0 0\nH -0.24 0.93 0\n"


@pytest.fixture
def recording_xtb(tmp_path, monkeypatch):
    """An xtb on PATH that answers --version and records its BLAS threads."""
    from delfin.dashboard import gfn_optimize as gfn
    from delfin.qm_runtime import clear_resolver_cache

    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    marks = tmp_path / "marks.txt"
    script = bin_dir / "xtb"
    script.write_text(
        "#!/bin/sh\n"
        'case "$*" in *--version*) echo "xtb version 6.7.1 (recorder)"; '
        "exit 0;; esac\n"
        f'echo "${{OPENBLAS_NUM_THREADS-unset}},${{OMP_WAIT_POLICY-unset}}"'
        f' >> "{marks}"\n'
        "exit 1\n", encoding="utf-8")
    script.chmod(script.stat().st_mode | stat.S_IXUSR)
    monkeypatch.setenv("PATH", f"{bin_dir}{os.pathsep}/usr/bin{os.pathsep}/bin")
    monkeypatch.setenv("DELFIN_QM_TOOLS_ROOT", str(tmp_path / "no-tools"))
    monkeypatch.delenv("XTBHOME", raising=False)
    monkeypatch.delenv("OPENBLAS_NUM_THREADS", raising=False)
    monkeypatch.delenv("OMP_WAIT_POLICY", raising=False)
    gfn._XTB_JUDGED.clear()
    clear_resolver_cache()
    assert gfn.find_xtb() == str(script)
    yield marks
    gfn._XTB_JUDGED.clear()
    clear_resolver_cache()


def _recorded(marks, field=0):
    """One field per xtb start: 0 is OPENBLAS_NUM_THREADS, 1 OMP_WAIT_POLICY."""
    if not marks.exists():
        return []
    return [line.split(",")[field]
            for line in marks.read_text(encoding="utf-8").split()]


def test_an_optimisation_gives_blas_one_thread(recording_xtb):
    from delfin.dashboard import gfn_optimize as gfn

    gfn.optimize_with_gfn(_XYZ, "gfn2", max_steps=1, timeout=30)
    assert _recorded(recording_xtb) and set(_recorded(recording_xtb)) == {"1"}


def test_the_climb_gives_blas_one_thread(recording_xtb):
    from delfin.dashboard import climb

    walk = climb.Climb(_XYZ, "gfn2")
    try:
        walk.modes_from_xtb()
    finally:
        walk.close()
    line = climb._CommandLine([8, 1, 1], walk.bohr, "gfn2", 0, 0, None, 4)
    try:
        with pytest.raises(Exception):
            line(walk.bohr)
    finally:
        line.close()
    got = _recorded(recording_xtb)
    assert len(got) >= 2 and set(got) == {"1"}, got



def test_xtbs_openmp_threads_sleep_rather_than_spin_while_waiting(
        recording_xtb):
    """OMP_WAIT_POLICY=PASSIVE on the same xtb runs.

    Input and output as above, for the second variable. xtb's OpenMP
    threads busy-wait at every barrier by default; on CPUs something else
    also wants, that turns a step into a queue. Measured 2026-10-16, the
    climb test's 16-atom estimate, median of eight Climb.step() calls, 4
    CPUs: idle 18 ms default, 16 ms PASSIVE; with four busy processes on
    the same CPUs 280 ms default, 22 ms PASSIVE. relax_steps (20 GFN2
    cycles) with CPUs to spare: glycine 68/68 ms, ibuprofen 220/219 ms,
    eicosane 270/270 ms. The one case that lost: eicosane asking for 16
    threads on 8 CPUs, 319 ms default against 368 ms PASSIVE.
    """
    from delfin.dashboard import climb
    from delfin.dashboard import gfn_optimize as gfn

    gfn.optimize_with_gfn(_XYZ, "gfn2", max_steps=1, timeout=30)
    walk = climb.Climb(_XYZ, "gfn2")
    try:
        walk.modes_from_xtb()
    finally:
        walk.close()
    line = climb._CommandLine([8, 1, 1], walk.bohr, "gfn2", 0, 0, None, 4)
    try:
        with pytest.raises(Exception):
            line(walk.bohr)
    finally:
        line.close()
    got = _recorded(recording_xtb, field=1)
    assert len(got) >= 3 and set(got) == {"PASSIVE"}, got

def test_a_walk_under_load_minimises_on_one_blas_thread():
    """The parent's own BLAS pools are held to one thread while it drives xtb.

    Input: relax_under_load with an engine that records, at every call, the
    largest thread count of any BLAS pool loaded in this process
    (threadpoolctl: numpy's and scipy's OpenBLAS).
    Output: every recorded count is 1.

    relax_under_load alternates scipy's L-BFGS-B with one xtb process per
    gradient. Left at their defaults, numpy's and scipy's OpenBLAS pools have
    one thread per available CPU, and after each call those threads go on
    spinning while the next xtb runs on the same CPUs. Measured 2026-10-16,
    test_pulling_on_a_structure.py::
    test_a_push_stops_at_a_conformational_product_not_past_it on 4 CPUs,
    conda-forge xtb with OPENBLAS_NUM_THREADS=1 already in its environment:
    141 s; with both parent pools held to one thread, 21 s (1292 gradients,
    117 ms and 14 ms each); with CPUs to spare, 23 s either way.
    """
    pytest.importorskip("threadpoolctl")
    import numpy as np
    from threadpoolctl import threadpool_info

    from delfin.dashboard import under_load

    seen = []

    class _Engine:
        def __call__(self, bohr):
            pools = [one["num_threads"] for one in threadpool_info()
                     if one.get("user_api") == "blas"]
            seen.append(max(pools) if pools else 1)
            flat = np.asarray(bohr, dtype=float).ravel()
            return 0.5 * float(flat @ flat), flat.reshape(-1, 3)

    class _Loaded:
        engine = _Engine()

    start = np.full((3, 3), 0.3)
    got = under_load.relax_under_load(_Loaded(), start, np.zeros((3, 3)),
                                      steps=5)
    assert got["ok"], got["status"]
    assert seen and set(seen) == {1}, seen
