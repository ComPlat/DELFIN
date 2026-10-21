"""The interpreter's own C++ runtime, mapped before a wheel picks another.

Input: none. Output of ``prefer_own_libstdcxx``: the path of the
``libstdc++.so.6`` this process has mapped afterwards, or "" (not Linux,
nothing mapped and nothing to map). Never raises; idempotent.

Semantics: a process maps ``libstdc++.so.6`` once, and every library
loaded later that needs it gets that copy. A wheel built with
``RUNPATH $ORIGIN`` and no runtime of its own (PyMuPDF) is resolved
through ``LD_LIBRARY_PATH`` first. Where that variable names an older
runtime -- an ORCA installation puts ``/opt/orca/lib`` there, holding
GCC 12's 6.0.32 -- the first such import maps the old copy, and from then
on the interpreter's own libraries that need a newer one fail to load
(``import sqlite3``: ICU 78 needs ``CXXABI_1.3.15``; with it ``stk``,
``atomlite`` and ``delfin.smiles_converter``). Mapping the interpreter's
copy (``<sys.prefix>/lib/libstdc++.so.6``) first, when one exists and no
copy is mapped yet, makes the order irrelevant: the runtime is backward
compatible, so the wheel works with the newer one.

Called right before every ``import pymupdf`` in DELFIN. A runtime that
is already mapped is left alone: a process cannot swap it.

Rejected: removing the foreign directory from ``LD_LIBRARY_PATH``. The
variable is read once at process start, so changing it in-process does
not affect this process, and children that run ORCA need it. Rejected:
importing ``sqlite3`` first as a side effect. It maps the right runtime
only where ICU links it, which is a build detail of one distribution.
Cost: one read of ``/proc/self/maps`` (well under a millisecond) per
call; one ``dlopen`` the first time.
"""

from __future__ import annotations

import sys
from pathlib import Path

_SONAME = "libstdc++.so.6"


def _mapped() -> str:
    try:
        with open("/proc/self/maps", encoding="utf-8",
                  errors="replace") as fh:
            for line in fh:
                if "libstdc++" in line:
                    return line.split()[-1]
    except OSError:
        pass
    return ""


def prefer_own_libstdcxx() -> str:
    """Map the interpreter's ``libstdc++.so.6`` unless one is mapped."""
    if not sys.platform.startswith("linux"):
        return ""
    try:
        current = _mapped()
        if current:
            return current
        own = Path(sys.prefix) / "lib" / _SONAME
        if not own.exists():
            return ""
        import ctypes
        ctypes.CDLL(str(own), mode=ctypes.RTLD_GLOBAL)
        return _mapped()
    except Exception:
        return ""
