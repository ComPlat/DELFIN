"""Result checker for chemistry answers.

Checks a finished ORCA run against the failure modes known from
autonomous chemistry agents: values guessed instead of computed,
overconfident analysis, and wrong job outcomes accepted as correct.

:func:`check_run` performs four families of checks on a run:

* **Traceability** — every number the agent reports must trace to an
  output file it read. A reported value with no ``source`` is flagged
  as *not computed*. This is the agent-facing half of DELFIN's claim
  check: the agent attaches the file it read to each number it reports,
  and this checker refuses numbers that carry no source.
* **Job-type consistency** — a minimum has 0 imaginary modes, a
  transition state exactly 1 (and its mode is named, i.e. its index and
  frequency are reported); SCF and geometry are converged; for an
  unrestricted (UKS) run the ``<S**2>`` deviation sits within a stated
  bound; a T1 diagnostic is reported for a DLPNO-CCSD(T) run and
  flagged when it exceeds the manual's threshold.
* **Plausibility** — when the task supplies reference values, each
  reported number is compared against them within a stated tolerance.

Values are taken through :mod:`delfin.orca_read` where the quantity
exists (energies, ``<S**2>`` deviation, termination, geometry
convergence, multiplicity); imaginary modes, SCF convergence and the T1
diagnostic have no ``Quantity`` member, so the small text readers below
answer those — all without OPI, so this module runs on the CI
interpreter.

The report is a dict: ``ok`` is the overall verdict; ``checks`` lists
each individual check with ``passed``, a ``check`` name and a
``detail`` string. A check marked ``untested`` is a warning (it could
not be verified) and does not by itself flip the verdict to failed —
see ``s2_bound`` on a UKS run with no bound given.
"""
from __future__ import annotations

import re
from pathlib import Path

from delfin.imag import collect_imaginary_modes
from delfin.orca_read import Quantity, read

#: The manual's rule of thumb for a T1 diagnostic on a DLPNO-CCSD(T)
#: run: a value above this indicates significant static (single-
#: reference-invalidating) correlation and is flagged. This threshold
#: is the one cited throughout the ORCA manual's discussion of scalar
#: static-correlation diagnostics (e.g. in the FOD section, which
#: cross-references the T1 diagnostic).
T1_THRESHOLD = 0.02

_SCF_CONVERGED_MARKERS = (
    ("SCF CONVERGED AFTER", True),
    ("SCF CONVERGED", True),
    ("SCF NOT CONVERGED", False),
    ("SCF iterations did not converge", False),
)

_GEOM_CONVERGED_MARKERS = (
    ("OPTIMIZATION RUN DONE", True),
    ("THE OPTIMIZATION HAS CONVERGED", True),
    ("OPTIMIZATION DID NOT CONVERGE", False),
    ("FAILED TO CONVERGE THE GEOMETRY OPTIMIZATION", False),
)

# Both spellings ORCA prints: the DLPNO module's "T1(DLPNO) : 0.0104"
# and the FOD-style "T1 diagnostic ... 0.006287935".
_T1_RE = re.compile(
    r"T1(?:\s+diagnostic|\(DLPNO\))\s*(?:\.\.\.)?\s*:?\s*([-+]?\d+\.\d+)"
)

_JOB_TYPES = ("minimum", "ts", "sp", "other")


def _scf_converged(text: str):
    """True/False when the run states its SCF outcome, else None."""
    for marker, value in _SCF_CONVERGED_MARKERS:
        if marker in text:
            return value
    return None


def _geometry_converged(reading) -> bool:
    """Geometry opt converged; None when the run did not optimise."""
    return reading.value if isinstance(reading.value, bool) else None


def _imag_count_and_modes(out_path: Path):
    """(count, modes) of the imaginary modes an ORCA output prints.

    ``modes`` is the list of ``(index, frequency)`` for the imaginary
    modes (most negative first), using :func:`delfin.imag.collect_imag-
    inary_modes`, which reads the ``***imaginary mode***`` lines.
    """
    modes = collect_imaginary_modes(str(out_path))
    return len(modes), modes


def _t1_diagnostic(text: str):
    """The T1 diagnostic value an output prints, or None."""
    m = _T1_RE.search(text)
    return float(m.group(1)) if m else None


def _check(name: str, passed: bool, detail: str,
           untested: bool = False, **extra) -> dict:
    entry = {"check": name, "passed": passed, "detail": detail}
    if untested:
        entry["untested"] = True
    entry.update(extra)
    return entry


def check_run(work_dir, base_name: str, job_type: str,
              reported=None, s2_bound=None, references=None,
              method: str | None = None) -> dict:
    """Validate a finished ORCA run and the numbers an agent reports.

    Args:
        work_dir: directory holding ``<base_name>.out`` (and ``.inp``).
        base_name: the run's base name (its ``.out`` without extension).
        job_type: what the agent claims the run is: ``"minimum"``,
            ``"ts"``, ``"sp"`` (single point) or ``"other"``.
        reported: ``{label: {"value": number, "source": path-or-None}}``
            — the numbers the agent reports about this run and the
            output file each number was read from. A number with no
            ``source`` is flagged *not computed*.
        s2_bound: for an unrestricted (UKS) run, the largest acceptable
            ``<S**2>`` deviation. When a UKS run has no bound, the
            check is reported *untested* (a warning), not failed.
        references: ``{label: {"expected": number, "tolerance": number}}``
            — plausibility targets supplied by the task.
        method: the electronic-structure method (e.g. ``"B3LYP"``,
            ``"DLPNO-CCSD(T)"``). Only a DLPNO-CCSD(T) run must report a
            T1 diagnostic.

    Returns:
        A dict report (see the module docstring). Never raises for a
        missing or partial run: unreadable sub-results are ``None`` and
        the checks they feed fail or are marked untested.
    """
    work_dir = Path(work_dir)
    out_path = work_dir / f"{base_name}.out"
    reported = reported or {}
    references = references or {}
    checks = []
    result = {"work_dir": str(work_dir), "base_name": base_name,
              "job_type": job_type}

    # --- read-only sub-results via delfin.orca_read -----------------------
    reading_s2 = read(work_dir, base_name, Quantity.S2_DEVIATION)
    reading_term = read(work_dir, base_name, Quantity.TERMINATED_NORMALLY)
    reading_geom = read(work_dir, base_name, Quantity.GEOMETRY_CONVERGED)
    reading_mult = read(work_dir, base_name, Quantity.MULTIPLICITY)

    try:
        text = out_path.read_text(encoding="utf-8", errors="replace")
    except FileNotFoundError:
        text = ""
    except OSError:
        text = ""

    img_count, img_modes = _imag_count_and_modes(out_path) if out_path.is_file() else (0, [])
    t1 = _t1_diagnostic(text) if text else None

    # --- traceability ----------------------------------------------------
    trace = {}
    reported_ok = True
    for label, item in reported.items():
        source = None if isinstance(item, dict) else None
        if isinstance(item, dict):
            source = item.get("source")
        not_computed = source is None or str(source).strip() == ""
        trace[label] = {
            "value": item.get("value") if isinstance(item, dict) else item,
            "source": source,
            "not_computed": not_computed,
        }
        if not_computed:
            reported_ok = False
            checks.append(_check(
                "traceability", False,
                f"value '{label}' has no source file it was read from; "
                "flagged as not computed"))
    result["traceability"] = trace
    if reported_ok:
        checks.append(_check("traceability", True,
                             f"all {len(reported)} reported value(s) carry a source"))

    # --- convergence (SCF, geometry, termination) ------------------------
    scf = _scf_converged(text) if text else None
    geom = _geometry_converged(reading_geom)
    terminated = reading_term.value if isinstance(reading_term.value, bool) else None
    result["convergence"] = {
        "scf_converged": scf,
        "geometry_converged": geom,
        "terminated_normally": terminated,
    }
    if scf is True:
        checks.append(_check("scf_converged", True, "SCF converged"))
    elif scf is False:
        checks.append(_check("scf_converged", False, "SCF did not converge"))
    elif terminated is True:
        # A run that ended normally but did not state its SCF outcome is
        # itself suspicious — every finished ORCA job prints the marker.
        checks.append(_check("scf_converged", False,
                             "SCF convergence not stated, though the run "
                             "terminated normally"))
    else:
        checks.append(_check("scf_converged", False,
                             "SCF convergence not stated in the output",
                             untested=True))
    if job_type == "minimum":
        if geom is True:
            checks.append(_check("geometry_converged", True,
                                 "geometry optimisation converged"))
        elif geom is False:
            checks.append(_check("geometry_converged", False,
                                 "geometry optimisation did not converge"))
        elif terminated is True:
            checks.append(_check("geometry_converged", False,
                                 "geometry optimisation outcome not stated, "
                                 "though the run terminated normally"))
        else:
            checks.append(_check("geometry_converged", False,
                                 "geometry optimisation outcome not stated",
                                 untested=True))
    else:
        # A "ts" verdict comes from a frequency run, which need not
        # optimise; a "sp"/"other" run never optimises. Geometry
        # convergence is not required for these.
        result["convergence"]["geometry_required"] = False

    # --- job-type consistency (imaginary modes) --------------------------
    consistent_job = None
    if job_type == "minimum":
        consistent_job = (img_count == 0)
        reason = (f"{img_count} imaginary mode(s); a minimum must have 0"
                  if not consistent_job
                  else "no imaginary modes, as a minimum requires")
    elif job_type == "ts":
        consistent_job = (img_count == 1)
        reason = ("exactly one imaginary mode, which is reported as "
                  f"mode {img_modes[0][0]} at {img_modes[0][1]} cm**-1"
                  if consistent_job
                  else f"{img_count} imaginary mode(s); a TS must have "
                       "exactly 1, and its mode must be named")
    elif job_type == "sp":
        consistent_job = True
        reason = "imaginary-mode count not applicable to a single point"
    else:  # "other"
        consistent_job = True
        reason = "job type not checked against imaginary modes"
    result["job_type"] = {
        "consistent": consistent_job,
        "imag_freq_count": img_count,
        "imag_modes": img_modes,
        "reason": reason,
    }
    checks.append(_check("job_type", bool(consistent_job), reason))

    # --- spin: <S**2> deviation for UKS ----------------------------------
    mult = reading_mult.value if isinstance(reading_mult.value, int) else None
    s2 = reading_s2.value if isinstance(reading_s2.value, float) else None
    is_uks = mult not in (None, 1)
    if is_uks:
        if s2_bound is not None and s2 is not None:
            within = s2 <= s2_bound
            checks.append(_check(
                "s2_bound", within,
                f"<S**2> deviation {s2:.6g} within bound {s2_bound:g}: "
                f"{'yes' if within else 'no'}"))
        elif s2 is not None:
            checks.append(_check(
                "s2_bound", False,
                f"UKS run: <S**2> deviation {s2:.6g} needs an explicit "
                "bound to be judged",
                untested=True))
        else:
            checks.append(_check("s2_bound", False,
                                 "UKS run: <S**2> deviation not found",
                                 untested=True))
    else:
        checks.append(_check("s2_bound", True,
                             "not an unrestricted run; no <S**2> bound needed"))
    result["spin"] = {
        "multiplicity": mult,
        "s2_deviation": s2,
        "s2_bound": s2_bound,
        "within_bound": (s2 <= s2_bound) if (is_uks and s2_bound is not None and s2 is not None) else None,
    }

    # --- T1 diagnostic for DLPNO-CCSD(T) ---------------------------------
    is_dlpno = bool(method) and "DLPNO-CCSD(T)" in method
    if is_dlpno:
        if t1 is not None:
            flagged = t1 > T1_THRESHOLD
            if flagged:
                detail = (f"T1 diagnostic {t1:.8g} reported and flagged: it "
                          f"exceeds the manual threshold {T1_THRESHOLD} "
                          "(static correlation)")
            else:
                detail = (f"T1 diagnostic {t1:.8g} reported, below the "
                          f"manual threshold {T1_THRESHOLD}: "
                          "single-reference is fine")
            checks.append(_check(
                "t1_diagnostic", True, detail, flagged=flagged))
            result["t1_diagnostic"] = {
                "required": True, "present": True, "value": t1,
                "threshold": T1_THRESHOLD, "flagged": flagged,
            }
        else:
            checks.append(_check(
                "t1_diagnostic", False,
                "DLPNO-CCSD(T) run must report a T1 diagnostic, but the "
                "output has none"))
            result["t1_diagnostic"] = {
                "required": True, "present": False, "value": None,
                "threshold": T1_THRESHOLD, "flagged": False,
            }
    else:
        result["t1_diagnostic"] = {"required": False}

    # --- plausibility against reference values ---------------------------
    plaus = {}
    plaus_ok = True
    for label, ref in references.items():
        expected = ref.get("expected")
        tolerance = ref.get("tolerance")
        reported_value = _reported_value(reported, label)
        if reported_value is None:
            plaus[label] = {"within_tolerance": False,
                            "detail": "no reported value to compare"}
            plaus_ok = False
            checks.append(_check(
                "plausibility", False,
                f"reference for '{label}' but no reported value"))
            continue
        within = abs(reported_value - expected) <= tolerance
        plaus[label] = {
            "reported": reported_value, "expected": expected,
            "tolerance": tolerance, "within_tolerance": within,
        }
        if not within:
            plaus_ok = False
            checks.append(_check(
                "plausibility", False,
                f"'{label}' = {reported_value:.8g} deviates from reference "
                f"{expected:.8g} by more than the stated tolerance {tolerance:g}"))
    result["plausibility"] = plaus
    if plaus_ok and references:
        checks.append(_check("plausibility", True,
                             f"all {len(references)} reported value(s) within "
                             "the stated reference tolerances"))

    # --- overall verdict: any hard (non-untested) failed check -----------
    failed = [c for c in checks if c["passed"] is False and not c.get("untested")]
    result["checks"] = checks
    result["ok"] = not failed
    return result


def _reported_value(reported, label):
    """The numeric value an agent reported for *label*, or None."""
    item = reported.get(label)
    if isinstance(item, dict):
        return item.get("value")
    return item
