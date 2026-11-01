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

__all__ = [
    "rate_problem",
    "window_remainders",
    "TRUNCATION_LIMIT",
    "ISC_REFERENCE",
    "REFERENCE_LEVEL",
    "REFERENCE_RATIO_BAND",
    "isc_level_matches",
    "isc_reference_comparison",
]

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


#: The in-tree converged ISC reference (this module's docstring above,
#: measured on formaldehyde with ORCA 6.1.1, TD-DFT, LINEW 50,
#: MAXTIME auto): ISC S1>T1 converged 2.46e7 s-1.  A comparison is
#: made ONLY when the run's level matches REFERENCE_LEVEL.
ISC_REFERENCE = {
    "rate": 2.46e7,
    "source": "delfin/common/esd_numerics.py docstring "
    "(ORCA 6.1.1 formaldehyde, TD-DFT, LINEW 50, MAXTIME auto)",
}

#: The level the reference was measured at.  MAXTIME auto cannot be
#: read from the output (ORCA prints the chosen window instead), so
#: it is not checked; the docstring records it.
REFERENCE_LEVEL = {
    "functional": "PBE0",
    "disp_corr": "D4",
    "main_basisset": "def2-SVP",
    "solvent": "",
    "linew": 50.0,
}

#: How far the linewidth may differ and still count as the same
#: damping.
_LINEW_TOLERANCE_CM1 = 0.5

#: Ratio band against the reference: outside it the comparison FAILS.
#: Wide (0.2-5.0) to catch a wrong sign, an unscaled unit or a
#: truncated window reported as converged -- not to be a precision
#: benchmark.
REFERENCE_RATIO_BAND = (0.2, 5.0)

_LINEWIDTH_RE = re.compile(r"Homogeneous\s+linewidth\s+is:\s*([0-9.]+)\s*cm-1", re.I)


def isc_level_matches(control: dict, esd_out_text: Optional[str]) -> "tuple[bool, str]":
    """Verify the run's level against REFERENCE_LEVEL before comparing.

    Checks functional, dispersion, basis and solvent from the CONTROL
    dict and the homogeneous linewidth from the ESD output text.
    Every mismatch is named in the reason; the reason also says when
    the linewidth could not be checked.
    """
    mism = []
    if str(control.get("functional", "")).upper() != REFERENCE_LEVEL["functional"]:
        mism.append(f"functional={control.get('functional')}")
    if str(control.get("disp_corr", "")).upper() != REFERENCE_LEVEL["disp_corr"]:
        mism.append(f"disp_corr={control.get('disp_corr')}")
    if str(control.get("main_basisset", "")).upper() != REFERENCE_LEVEL["main_basisset"].upper():
        mism.append(f"basis={control.get('main_basisset')}")
    if str(control.get("solvent", "")) != REFERENCE_LEVEL["solvent"]:
        mism.append(f"solvent={control.get('solvent')!r}")
    linewidth_checked = False
    if esd_out_text:
        m = _LINEWIDTH_RE.search(esd_out_text)
        if m:
            linewidth_checked = True
            if abs(float(m.group(1)) - REFERENCE_LEVEL["linew"]) > _LINEW_TOLERANCE_CM1:
                mism.append(f"LINEW={m.group(1)}")
    if mism:
        return False, "; ".join(mism)
    if not linewidth_checked:
        return True, ("level matches PBE0-D4/def2-SVP gas phase "
                      "(linewidth not stated in the output; not checked)")
    return True, "level matches PBE0-D4/def2-SVP gas phase, LINEW 50"


def isc_reference_comparison(
    rate: Optional[float],
    problem: Optional[str],
    level_ok: bool,
    level_reason: str = "",
) -> dict:
    """Verdict record for ISC rate vs ISC_REFERENCE, one consumer's contract.

    The comparison is level-gated: without a matching level, a rate, or
    with a rate that has a problem (esd_numerics.rate_problem), the
    verdict is NOT GRADED with the reason -- never a number.
    """
    if not level_ok:
        return {
            "criterion": "isc_vs_reference",
            "verdict": "NOT GRADED",
            "value": f"level differs: {level_reason}",
            "unit": "",
            "source": "CONTROL.txt + ESD output via isc_level_matches",
        }
    if rate is None:
        return {
            "criterion": "isc_vs_reference",
            "verdict": "NOT GRADED",
            "value": "no rate to compare",
            "unit": "",
            "source": "ESD ISC output",
        }
    if problem:
        return {
            "criterion": "isc_vs_reference",
            "verdict": "NOT GRADED",
            "value": f"rate has problem: {problem}",
            "unit": "",
            "source": "esd_numerics.rate_problem on the ISC output",
        }
    ratio = rate / ISC_REFERENCE["rate"]
    low, high = REFERENCE_RATIO_BAND
    return {
        "criterion": "isc_vs_reference",
        "verdict": "PASS" if low <= ratio <= high else "FAIL",
        "value": f"rate/ref = {ratio:.2f} (ref {ISC_REFERENCE['rate']:.3g} s^-1, "
        "delfin/common/esd_numerics.py docstring)",
        "unit": "ratio",
        "source": f"in-tree measurement, matched level verified: {ISC_REFERENCE['source']}",
    }


#: Experimental adiabatic singlet-triplet splitting of formaldehyde and
#: its provenance.  Band 0.25-0.55 eV around the operator-provided
#: origins (T1 ~25190 cm-1, S1 ~28190 cm-1 -> ~3000 cm-1 = 0.37 eV):
#: the review is Clouthier & Ramsay, Annu. Rev. Phys. Chem. 34, 31-58
#: (1983), DOI 10.1146/annurev.pc.34.100183.000335 (DOI verified on
#: annualreviews.org via DELFIN's web tools); the specific origins are
#: operator-provided and were NOT independently verifiable from that
#: session (publisher 403).  A criterion that grades against these
#: numbers must say so, not silently assume them.
EST_REFERENCE = {
    "band_ev": (0.25, 0.55),
    "centre_ev": 0.37,
    "source": (
        "Clouthier & Ramsay, Annu. Rev. Phys. Chem. 34, 31 (1983), "
        "DOI 10.1146/annurev.pc.34.100183.000335; origins "
        "(T1 ~25190 cm-1, S1 ~28190 cm-1) operator-provided, not "
        "independently verified (publisher 403)"
    ),
}


def est_reference_comparison(delta_est_ev: Optional[float]) -> dict:
    """Verdict record for the adiabatic Delta E_ST vs the experimental band.

    Sign first: a Delta E_ST of the wrong sign kills the TADF ranking
    claim regardless of magnitude.  Then the cited experimental band:
    a value outside 0.25-0.55 eV FAILS, never passes a loose
    |dE_ST|<1 eV check (which accepted 0.889 eV -- 2.4x the
    experimental adiabatic splitting -- on job 8181368, reviewer
    finding F6).  Without a value the verdict is NOT GRADED.
    """
    if delta_est_ev is None:
        return {
            "criterion": "delta_E_ST_vs_reference",
            "verdict": "NOT GRADED",
            "value": "no Delta E_ST to compare",
            "unit": "",
            "source": "ESD state energies via collect_esd_results",
        }
    low, high = EST_REFERENCE["band_ev"]
    # Sign convention: Delta E_ST = E(S1) - E(T1), so TADF-favourable
    # (T1 below S1) is POSITIVE here; a negative value means T1 above
    # S1 and kills the ranking claim.
    if delta_est_ev < 0:
        return {
            "criterion": "delta_E_ST_vs_reference",
            "verdict": "FAIL",
            "value": (
                f"Delta E_ST = {delta_est_ev:.3f} eV with the wrong sign "
                "(T1 not below S1): the TADF ranking claim fails"
            ),
            "unit": "eV",
            "source": f"experimental band source: {EST_REFERENCE['source']}",
        }
    inside = low <= delta_est_ev <= high
    return {
        "criterion": "delta_E_ST_vs_reference",
        "verdict": "PASS" if inside else "FAIL",
        "value": (
            f"Delta E_ST = {delta_est_ev:.3f} eV vs experimental adiabatic "
            f"{EST_REFERENCE['centre_ev']:.2f} eV (band {low:.2f}-{high:.2f} eV)"
        ),
        "unit": "eV",
        "source": f"experimental band source: {EST_REFERENCE['source']}",
    }
