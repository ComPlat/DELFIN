"""Scientific-correctness critic for finished ORCA calculations.

A *result* can be numerically "done" (ORCA TERMINATED NORMALLY) yet
scientifically wrong: an optimisation that landed on a saddle point
(imaginary frequencies), a broken-symmetry/radical run with heavy spin
contamination, an SCF that never converged, a transition-state search
with the wrong number of imaginary modes. Trusting those numbers is the
expensive failure mode for a chemistry agent.

This module is the agent's **correctness gate**: a read-only scan that
flags the red flags BEFORE the user trusts the result. It reuses
DELFIN's own parsers (``imag``, ``parser``, ``energies``) rather than
re-implementing output parsing — the agent learns from DELFIN, it does
not rebuild it. Model-agnostic, deterministic, never raises.

The critic NEVER edits anything and NEVER weakens convergence — it only
reports. Acting on a flag (re-running tighter, reviewing the method) is
a separate, human-confirmed decision.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional


# Spin-contamination deviation thresholds (|<S^2> - S(S+1)|).
_SPIN_WARN = 0.10
_SPIN_ERROR = 1.0
# Imaginary modes with magnitude below this are usually numerical noise
# (low-frequency rotations / soft modes), not a real saddle point.
_IMAG_NOISE_CM = 50.0


@dataclass
class Critique:
    level: str        # "error" | "warn" | "ok"
    code: str         # short machine label
    message: str      # human-readable, English

    def __str__(self) -> str:
        icon = {"error": "❌", "warn": "⚠️", "ok": "✅"}.get(self.level, "•")
        return f"{icon} {self.message}"


def _read(out_path: Path, limit: int = 2_000_000) -> str:
    try:
        return out_path.read_text(encoding="utf-8", errors="replace")[:limit]
    except Exception:
        return ""


def _job_kind(content: str) -> set[str]:
    """Detect what the calculation intended, from ORCA's own echo."""
    low = content.lower()
    kinds: set[str] = set()
    if "geometry optimization cycle" in low or "*** optimization run ***" in low:
        kinds.add("opt")
    if "vibrational frequencies" in low:
        kinds.add("freq")
    # Transition-state search markers (input echo / driver banner).
    if "optts" in low or "transition state" in low or "ts optimization" in low:
        kinds.add("ts")
    return kinds


def _check_termination(content: str) -> list[Critique]:
    out: list[Critique] = []
    if "ORCA TERMINATED NORMALLY" in content:
        out.append(Critique("ok", "termination", "ORCA terminated normally."))
    else:
        out.append(Critique(
            "error", "no-termination",
            "No 'ORCA TERMINATED NORMALLY' — the run did not finish cleanly; "
            "any numbers below may be from an incomplete calculation."))
    if "SCF NOT CONVERGED" in content or "SCF not converged" in content:
        out.append(Critique(
            "error", "scf-not-converged",
            "SCF did not converge — the wavefunction is not valid; do NOT "
            "trust energies/properties. Investigate the method, not the "
            "convergence thresholds."))
    elif _has_soft_scf_warning(content):
        out.append(Critique(
            "warn", "scf-not-fully-converged",
            "SCF did not fully converge (soft warning in the output) — the "
            "run may still report itself green; treat energies/properties "
            "with caution. Loosening convergence thresholds is not the fix; "
            "reconsider the method or SCF settings."))
    return out


def _has_soft_scf_warning(content: str) -> bool:
    """A soft SCF-convergence warning, in DELFIN's own wording.

    Reuses the phrasing delfin/api.py:516 already classifies as an
    SCF-convergence failure ("SCF iterations did not converge") — no new
    vocabulary of our own. Excludes the hard failure already handled
    above; ORCA never prints both for the same SCF block, and when it
    does the hard error wins.
    """
    return ("SCF iterations did not converge" in content
            or "The SCF procedure does not converge" in content)


def _check_geometry(content: str, kinds: set[str]) -> list[Critique]:
    if "opt" not in kinds:
        return []
    low = content.lower()
    converged = ("the optimization has converged" in low
                 or "optimization converged" in low or "hurray" in low)
    if converged:
        findings = [Critique("ok", "geom-converged",
                             "Geometry optimization converged.")]
        # Wave-10/s1: a converged optimization WITHOUT a frequency
        # calculation has not confirmed its structure to be a minimum.
        # The freq check below only runs when a freq block exists, so
        # without this the run scored a clean 'ok' on the geometry.
        if "freq" not in kinds:
            findings.append(Critique(
                "warn", "opt-no-freq",
                "Optimization converged but no frequencies were computed — "
                "the structure is not confirmed to be a stationary point "
                "(minimum). Run a Freq calculation on this geometry."))
        return findings
    return [Critique(
        "warn", "geom-not-converged",
        "Optimization run detected but no convergence banner found — the "
        "geometry may not be a true stationary point yet.")]


def _check_frequencies(out_path: Path, kinds: set[str]) -> list[Critique]:
    if "freq" not in kinds:
        return []
    try:
        from delfin.imag import collect_imaginary_modes
        modes = collect_imaginary_modes(str(out_path)) or []
    except Exception:
        return []
    n = len(modes)
    if n == 0:
        if "ts" in kinds:
            return [Critique(
                "warn", "ts-no-imag",
                "Transition-state search but NO imaginary frequency — this is "
                "not a first-order saddle point; the TS was not located.")]
        return [Critique("ok", "no-imag",
                         "No imaginary frequencies — a genuine minimum.")]
    worst = min(f for _i, f in modes)        # most negative
    small = abs(worst) < _IMAG_NOISE_CM
    if "ts" in kinds:
        if n == 1:
            return [Critique("ok", "ts-one-imag",
                             f"Exactly one imaginary mode ({worst:.0f} cm^-1) — "
                             f"correct for a transition state.")]
        return [Critique(
            "error", "ts-many-imag",
            f"{n} imaginary modes ({worst:.0f} cm^-1 worst) — a transition "
            f"state must have exactly one; this is a higher-order saddle.")]
    # Intended (or assumed) minimum with imaginary modes present.
    level = "warn" if small else "error"
    noise = (" (small magnitude — likely numerical: try a finer integration "
             "grid or tighter opt, not looser convergence)" if small else "")
    return [Critique(
        level, "min-has-imag",
        f"{n} imaginary frequency(ies) ({worst:.0f} cm^-1 worst) on a "
        f"structure expected to be a minimum — it is a saddle point, not a "
        f"true minimum{noise}.")]


def _check_spin(out_path: Path) -> list[Critique]:
    try:
        from delfin.parser import extract_last_uhf_deviation
        dev = extract_last_uhf_deviation(str(out_path))
    except Exception:
        return []
    if dev is None:
        return []
    dev = abs(dev)
    if dev >= _SPIN_ERROR:
        return [Critique(
            "error", "spin-contamination",
            f"Severe spin contamination (<S^2> deviation {dev:.2f}) — the "
            f"UHF/UKS state is heavily mixed; energies are unreliable. "
            f"Reconsider the spin state / broken-symmetry setup.")]
    if dev >= _SPIN_WARN:
        return [Critique(
            "warn", "spin-contamination",
            f"Notable spin contamination (<S^2> deviation {dev:.2f}) — check "
            f"the multiplicity and whether the radical character is physical.")]
    return [Critique("ok", "spin-clean",
                     f"Spin contamination negligible (<S^2> deviation {dev:.2f}).")]


def _check_energy(out_path: Path) -> list[Critique]:
    try:
        from delfin.energies import find_electronic_energy
        e = find_electronic_energy(str(out_path))
    except Exception:
        return []
    if e is None:
        return []
    if e >= 0:
        return [Critique(
            "warn", "energy-nonneg",
            f"Final electronic energy is not negative ({e:.6f} Eh) — "
            f"molecular electronic energies should be < 0; suspicious.")]
    return [Critique("ok", "energy-ok",
                     f"Final electronic energy {e:.6f} Eh.")]


def _is_largest_out(out_path: Path) -> bool:
    """True when the critiqued file IS the largest .out in its folder.

    The DELFIN spectra parsers (delfin.api.extract_vibrational_modes /
    extract_excited_states) take a FOLDER and read the LARGEST .out in
    it, while this critic scans ONE file. To never silently parse a
    different file than the one being flagged, the spectra checks run
    only under this condition; otherwise the caller reports a visible
    'spectra-skipped' warning instead.
    """
    folder = out_path.parent
    try:
        outs = list(folder.glob("*.out"))
    except Exception:
        return False
    if not outs:
        return False
    largest = max(outs, key=lambda p: p.stat().st_size)
    return largest == out_path


def _n_atoms(content: str) -> int | None:
    """Atom count from ORCA's echo line ("Number of atoms ... N")."""
    import re as _re
    m = _re.search(r"Number of atoms\s+\.{2,}\s+(\d+)", content)
    return int(m.group(1)) if m else None


def _check_mode_count(content: str, modes) -> list[Critique]:
    """Real-mode count vs 3N-6 (non-linear) / 3N-5 (linear molecule).

    Source: standard vibrational analysis; ORCA manual 6.1.1 section
    4.6 "Vibrational Frequencies" (3N-6 vibrational normal modes for a
    non-linear molecule). Translational/rotational zero modes (the 6
    or 5 near-0 cm^-1 rows the parser returns) are excluded before
    counting, because the parser returns ALL rows (water fixture:
    9 rows = 6 zeros + 3 real). A linear marker in the output switches
    the expectation to 3N-5.
    """
    n = _n_atoms(content)
    if n is None:
        return []
    real = [m for m in modes if abs(m.frequency_cm) > 1.0]
    if not real:
        return []
    linear = "linear molecule" in content.lower()
    expected = 3 * n - 5 if linear else 3 * n - 6
    if len(real) == expected:
        return []
    return [Critique(
        "error", "freq-count-mismatch",
        f"{len(real)} real vibrational modes for {n} atoms"
        f"{' (linear)' if linear else ''} — expected 3N-{'5' if linear else '6'}"
        f" = {expected}; the Hessian is likely incomplete or the run "
        f"finished partially. Do not trust this spectrum.")]


def _check_ir_intensities(modes) -> list[Critique]:
    """Negative IR intensity is unphysical for an absorption spectrum.

    Beer-Lambert law: the absorption cross-section (and hence the
    integrated intensity in km/mol) is non-negative for an absorption
    spectrum; a negative value means a broken table, a wrong column
    mapping or a corrupted output — not a real spectrum.
    A zero intensity is physical (symmetry-forbidden mode) and is
    reported as an ok-level note, never as a warning or error.
    """
    out: list[Critique] = []
    neg = [m for m in modes if m.ir_intensity is not None and m.ir_intensity < 0]
    for m in neg:
        out.append(Critique(
            "error", "ir-neg-intensity",
            f"Mode {m.mode_index} at {m.frequency_cm:.2f} cm^-1 reports a "
            f"negative IR intensity ({m.ir_intensity:g} km/mol) — "
            f"absorption intensities are physically non-negative "
            f"(Beer-Lambert); the spectrum table is unreliable."))
    if any(m.ir_intensity == 0.0 for m in modes):
        out.append(Critique(
            "ok", "ir-zero-intensity",
            "At least one mode has zero IR intensity — physical for a "
            "symmetry-forbidden mode."))
    return out


def _check_ir_imaginary_peaks(content: str, modes) -> list[Critique]:
    """An imaginary mode listed as a spectrum peak is not a physical line.

    ORCA lists imaginary modes with a NEGATIVE frequency in the IR
    table's frequency column (delfin/ir_spectrum.py:95 regex accepts
    them). extract_vibrational_modes fills ``frequency_cm`` from the
    VIBRATIONAL FREQUENCIES block only, so the IR table's frequency
    column is scanned minimally here (QS ruling: first parser for a
    quantity nothing reads). A peak at a negative frequency is not an
    absorption line; reporting it in a spectrum is wrong.
    """
    import re as _re
    ir_re = _re.search(
        r"IR SPECTRUM.*?\n[-=]+\n(.*?)(?:\nThe first frequency|\n\n\n)",
        content, _re.DOTALL,
    )
    if not ir_re:
        return []
    row_re = _re.compile(r"^\s*(\d+):\s+(-?\d+\.\d+)", _re.MULTILINE)
    neg = [(int(m.group(1)), float(m.group(2)))
           for m in row_re.finditer(ir_re.group(1))
           if float(m.group(2)) < -1.0]
    if not neg:
        return []
    worst = min(f for _i, f in neg)
    return [Critique(
        "error", "ir-imaginary-as-peak",
        f"{len(neg)} imaginary mode(s) (down to {worst:.0f} cm^-1) appear "
        f"as IR spectrum peaks — imaginary frequencies are not physical "
        f"absorption lines; exclude them or fix the structure first.")]


def _largest_out_name(out_path: Path) -> str:
    """Name of the largest .out in the critiqued file's folder."""
    try:
        outs = list(out_path.parent.glob("*.out"))
        return max(outs, key=lambda p: p.stat().st_size).name
    except Exception:
        return "(unknown)"


def _check_spectra(out_path: Path, content: str) -> list[Critique]:
    """Spectrum plausibility checks, reusing DELFIN's own parsers.

    Reuses delfin.api.extract_vibrational_modes (folder-level, reads
    the LARGEST .out) — never re-parses ORCA output. Runs only when the
    critiqued file IS that largest .out (see _is_largest_out); otherwise
    a visible 'spectra-skipped' warning keeps 'file flagged' and 'file
    parsed' from diverging silently.
    """
    if not _is_largest_out(out_path):
        return [Critique(
            "warn", "spectra-skipped",
            f"Spectrum checks skipped: the folder's largest .out is "
            f"{_largest_out_name(out_path)}, not {out_path.name}, so "
            f"DELFIN's folder-level parsers would read a different file "
            f"than the one critiqued.")]
    findings: list[Critique] = []
    try:
        from delfin.api import extract_vibrational_modes
        res = extract_vibrational_modes(str(out_path.parent))
        if res.error is None and res.modes:
            findings += _check_mode_count(content, res.modes)
            findings += _check_ir_intensities(res.modes)
            findings += _check_ir_imaginary_peaks(content, res.modes)
    except Exception:
        findings.append(Critique(
            "warn", "spectra-parser-error",
            "Vibrational-spectrum check could not run (parser error); "
            "the spectrum was NOT verified."))
    return findings


def critique_output(out_path) -> list[Critique]:
    """Scan ONE ORCA ``.out`` for scientific red flags. Never raises."""
    try:
        p = Path(out_path)
        content = _read(p)
        if not content:
            return [Critique("warn", "empty", f"{p.name}: unreadable/empty.")]
        kinds = _job_kind(content)
        findings: list[Critique] = []
        findings += _check_termination(content)
        findings += _check_geometry(content, kinds)
        findings += _check_frequencies(p, kinds)
        findings += _check_spin(p)
        findings += _check_energy(p)
        findings += _check_spectra(p, content)
        return findings
    except Exception as exc:
        return [Critique("warn", "critic-error",
                         f"correctness scan failed: {exc}")]


def critique_folder(folder) -> dict:
    """Critique every ``.out`` in ``folder``. Returns {filename: [Critique]}."""
    try:
        d = Path(folder)
        if not d.is_dir():
            return {}
        return {p.name: critique_output(p) for p in sorted(d.glob("*.out"))}
    except Exception:
        return {}


def worst_level(critiques: list[Critique]) -> str:
    levels = {c.level for c in critiques}
    if "error" in levels:
        return "error"
    if "warn" in levels:
        return "warn"
    return "ok"


def format_report(by_file: dict) -> str:
    """Human-readable correctness report across a folder's outputs."""
    if not by_file:
        return "No .out files to check."
    lines: list[str] = []
    for fname, crits in by_file.items():
        wl = worst_level(crits)
        head = {"error": "❌ PROBLEMS", "warn": "⚠️ REVIEW",
                "ok": "✅ looks sound"}.get(wl, "•")
        lines.append(f"\n{fname} — {head}")
        for c in crits:
            lines.append(f"  {c}")
    lines.append(
        "\nThe critic only reports — it never edits results or weakens "
        "convergence. Act on ❌/⚠️ by reviewing the method or re-running "
        "tighter (with your confirmation).")
    return "\n".join(lines)
