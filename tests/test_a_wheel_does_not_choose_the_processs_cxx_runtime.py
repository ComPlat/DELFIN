"""PyMuPDF does not decide which C++ runtime the whole process gets.

A process maps ``libstdc++.so.6`` once; every later library that needs it
gets that copy. PyMuPDF is a wheel with ``RUNPATH $ORIGIN`` and no copy of
its own, so the dynamic linker searches ``LD_LIBRARY_PATH`` first. On a
host where that variable names an older runtime (an ORCA installation puts
``/opt/orca/lib`` there, holding GCC 12's 6.0.32), the first
``import pymupdf`` maps the old copy, and from then on the interpreter's
own libraries that need a newer one fail to load: ``import sqlite3`` (ICU
78 needs ``CXXABI_1.3.15``), and with it ``stk`` / ``atomlite`` and
``delfin.smiles_converter``. Measured in one pytest process:
``tests/test_doctor.py`` (the document-backend probe imports pymupdf)
followed by ``tests/test_the_manual_every_recipe_runs.py`` failed three
tests that pass alone.

Input: a fresh interpreter whose ``LD_LIBRARY_PATH`` names a directory
holding a COPY of the interpreter's own ``libstdc++.so.6`` -- a stand-in
for the foreign runtime that exists on every machine, so the test does not
depend on an ORCA install. Output: after ``office._pymupdf()`` the mapped
``libstdc++`` is the interpreter's own file, not the copy.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pytest

from conftest import child_env

_OWN = Path(sys.prefix) / "lib" / "libstdc++.so.6"

_PROBE = r"""
import sys
sys.path.insert(0, sys.argv[1])
from delfin.agent import office
office._pymupdf()
with open("/proc/self/maps") as fh:
    print(sorted({l.split()[-1] for l in fh if "libstdc++" in l}))
"""


def _repo() -> str:
    return str(Path(__file__).resolve().parents[1])


@pytest.mark.skipif(not sys.platform.startswith("linux"),
                    reason="/proc/self/maps and LD_LIBRARY_PATH are Linux")
@pytest.mark.skipif(not _OWN.exists(),
                    reason="the interpreter brings no C++ runtime of its own; "
                           "there is no newer copy a foreign one could displace")
def test_pymupdf_maps_the_interpreters_own_runtime(tmp_path):
    pytest.importorskip("pymupdf")
    foreign = tmp_path / "foreign_lib"
    foreign.mkdir()
    shutil.copy2(_OWN.resolve(), foreign / "libstdc++.so.6")
    env = child_env(tmp_path)
    env["LD_LIBRARY_PATH"] = str(foreign)
    out = subprocess.run([sys.executable, "-c", _PROBE, _repo()],
                         env=env, capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr[-2000:]
    mapped = out.stdout.strip().splitlines()[-1]
    assert str(foreign) not in mapped, (
        f"pymupdf mapped the runtime LD_LIBRARY_PATH named: {mapped}")
    assert str(_OWN.resolve()) in mapped, mapped


def test_preloading_never_raises_and_says_what_is_mapped():
    """Called in a process that has already mapped a runtime (pytest has,
    or has not -- either way) it returns a path or "" and raises nothing."""
    from delfin.agent import native_libs
    first = native_libs.prefer_own_libstdcxx()
    second = native_libs.prefer_own_libstdcxx()
    assert isinstance(first, str) and first == second
