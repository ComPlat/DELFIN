"""Whether a rate constant ORCA's ESD module printed is a converged number.

ORCA computes ISC, IC, fluorescence and phosphorescence rates by integrating
a correlation function over time.  Two things make the printed number not a
result, and ORCA flags only the first:

* a negative rate -- ORCA adds "WARNING: negative rates are unphysical! It
  means something went wrong with the CorrFunc integration";
* a time window that ends before the correlation function has decayed.  ORCA
  chooses its own window from the linewidth so that it decays to the cutoff
  (2934 fs for LINEW 50 cm-1 on 6.1.1); a MAXTIME set shorter truncates it
  without a word.  Measured on formaldehyde with MAXTIME 12000 a.u. (290 fs,
  7 % of the amplitude left): ISC S1>T1 5.70e7 s-1 against the converged
  2.46e7 s-1, IC S1>S0 -8.1e-3 against +3.79e-3 s-1.

The damping is read from the output: a Lorentzian linewidth (homogeneous,
HWHM) decays as exp(-g t), a Gaussian one (inhomogeneous, standard
deviation) as exp(-(s t)^2 / 2), a Voigt profile as their product.
"""

from __future__ import annotations

import math
import re
from typing import List, Optional

__all__ = ["rate_problem", "window_remainders", "TRUNCATION_LIMIT"]

_CM1_PER_HARTREE = 219474.6313705
_FS_PER_AU = 0.02418884326585747

#: Amplitude the correlation function may still have at the end of the
#: window.  At 1.8e-5 (4x the old window) the formaldehyde ISC rate agreed with
#: ORCA's own window to 3e-4; at 6.5e-2 (the old window) it was 2.3x too fast.
TRUNCATION_LIMIT = 1e-4

_NUM = r"([0-9]+(?:\.[0-9]*)?(?:[eE][-+]?[0-9]+)?)"
_HOMOGENEOUS = re.compile(rf"Homogeneous\s+linewidth\s+is:\s*{_NUM}\s*cm-1", re.I)
_INHOMOGENEOUS = re.compile(rf"Inhomogeneous\s+linewidth\s+is:\s*{_NUM}\s*cm-1", re.I)
_MAXTIME = re.compile(rf"Maximum\s+time:\s*{_NUM}\s*fs", re.I)
_RATE = re.compile(r"rate\s+constant\s+is\s+(-?[0-9.]+(?:[eE][-+]?[0-9]+)?)", re.I)
_NEGATIVE_WARNING = re.compile(r"negative\s+rates\s+are\s+unphysical", re.I)
_BANNER = re.compile(r"\${4,}\s+JOB NUMBER\s+\d+\s+\$+")


def _sections(text: str) -> List[str]:
    parts = _BANNER.split(text)
    return [p for p in parts if p.strip()] or [text]


def window_remainders(text: str) -> List[float]:
    """Amplitude left at the end of each correlation-function window in the output (one per job)."""
    remainders: List[float] = []
    for section in _sections(text):
        window = _MAXTIME.search(section)
        if not window:
            continue
        t_au = float(window.group(1)) / _FS_PER_AU
        homogeneous = _HOMOGENEOUS.search(section)
        inhomogeneous = _INHOMOGENEOUS.search(section)
        if not homogeneous and not inhomogeneous:
            continue
        exponent = 0.0
        if homogeneous:
            exponent += float(homogeneous.group(1)) / _CM1_PER_HARTREE * t_au
        if inhomogeneous:
            sigma = float(inhomogeneous.group(1)) / _CM1_PER_HARTREE
            exponent += 0.5 * (sigma * t_au) ** 2
        remainders.append(math.exp(-exponent) if exponent < 700 else 0.0)
    return remainders


def rate_problem(text: Optional[str]) -> Optional[str]:
    """Why the rate(s) in an ESD output are not a result, or None when they are."""
    if not text:
        return None
    if _NEGATIVE_WARNING.search(text) or any(float(m) < 0 for m in _RATE.findall(text)):
        return ("unphysical: ORCA's correlation-function integration gave a negative rate "
                "(ORCA: 'negative rates are unphysical')")
    left = max(window_remainders(text), default=0.0)
    if left > TRUNCATION_LIMIT:
        return (f"not converged: the correlation function still had {left:.1%} of its amplitude "
                f"when the time window (MAXTIME) ended")
    return None
