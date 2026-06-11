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
    return out


def _check_geometry(content: str, kinds: set[str]) -> list[Critique]:
    if "opt" not in kinds:
        return []
    low = content.lower()
    converged = ("the optimization has converged" in low
                 or "optimization converged" in low or "hurray" in low)
    if converged:
        return [Critique("ok", "geom-converged",
                         "Geometry optimization converged.")]
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
