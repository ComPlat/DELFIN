"""A reproducible random stream for Open Babel's force fields (#527).

Open Babel moves two atoms apart along a random unit vector whenever they come closer
than 0.1 A inside a distance term (bond, van der Waals, distance constraint), and nudges
an angle term along a random vector when the angle falls below 2.5 degrees.  The vector
comes from ``vector3::randomUnitVector``, which draws from a function-static
``OBRandom(true)``: its constructor seeds the C library's ``rand()`` from the microseconds
of the wall clock (``srand(gettimeofday().tv_usec % 70092)``) the first time any such
nudge happens in the process, and every later nudge reads the next ``rand()`` values.

So whenever a minimisation passes through such a contact, the result depends on the
clock at the first contact of the process and on how many random numbers earlier calls
consumed.  The same input then gives a different geometry in every process (and in every
repeat inside one process), with no thread, hash seed or time budget involved.

Fix: construct the static generator once (one warm-up nudge, so its clock seed happens
before ours), then set ``rand()`` to a fixed seed before every Open Babel minimisation.
The nudge keeps its job of breaking a degenerate contact, but its direction now depends
only on the input of that call.

Switch: ``DELFIN_OB_DETERMINISTIC_RNG=1`` (default off; off leaves every call untouched).
"""
from __future__ import annotations

import ctypes
import os
import threading

_SEED = 527
_LOCK = threading.Lock()
_STATE = {"libc": None, "warm": False, "failed": False}


def ob_rng_enabled() -> bool:
    return os.environ.get("DELFIN_OB_DETERMINISTIC_RNG", "0") == "1"


def reseed_ob_rng(pybel_module) -> bool:
    """Seed the C ``rand()`` that Open Babel's random nudges read.  Returns True when done.

    No-op (False) when the switch is off, when ``pybel_module`` is None, or when the C
    library cannot be reached; the caller then runs exactly as before.
    """
    if not ob_rng_enabled() or pybel_module is None or _STATE["failed"]:
        return False
    try:
        with _LOCK:
            if _STATE["libc"] is None:
                # The process-wide symbol table: the same rand()/srand() Open Babel links to.
                _STATE["libc"] = ctypes.CDLL(None)
            if not _STATE["warm"]:
                # Constructs the function-static OBRandom(true) inside randomUnitVector now,
                # so its wall-clock srand() cannot happen later, in the middle of a minimisation.
                _v = pybel_module.ob.vector3()
                _v.randomUnitVector()
                _STATE["warm"] = True
            _STATE["libc"].srand(ctypes.c_uint(_SEED))
        return True
    except Exception:
        _STATE["failed"] = True
        return False
