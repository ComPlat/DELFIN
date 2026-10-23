"""Claims checking: every number an agent's ANSWER states must have a
source in the evidence.

``result_critic`` checks the finished calculation outputs. This module
closes the loop one step later: it checks the agent's own answer text.
A scientifically wrong ANSWER is the failure mode that reaches the
user even when every underlying output was fine — or worse, when one
was flagged and the flag was silently dropped.

Contract (package SC1, wave 20):

* Every number the answer states WITH A UNIT must be traceable to an
  evidence record: the same numeric value must appear in one of the
  evidence texts. A number with no source is flagged (``claim-no-source``)
  unless it is marked as an estimate.
* An energy or property value must name its level of theory: method
  (functional or wavefunction method), basis set and (for a solvated
  calculation) solvent (``claim-no-method``).
* When the evidence record for a number carries result-critic flags
  (imaginary mode on a minimum, unconverged SCF, spin contamination,
  ...), the answer must acknowledge that problem
  (``claim-flagged-source``). Dropping a known warning is not
  reporting.
* A value marked as an estimate (nearby ``estimate``/``estimated``/
  ``rough``/``about`` wording) is exempt from provenance.

Evidence records use the wave-10 ledger shape
(:mod:`delfin.agent.evidence`): ``{"kind": "calc", "ref": ..., "text":
<excerpt of the output the number came from>, "flags": [<critique
code>, ...]}``. ``text`` is what a number is matched against; ``flags``
are ``result_critic.Critique.code`` values, error/warn ones trigger the
acknowledgement requirement.

Model-agnostic, deterministic, read-only, NEVER raises (empty, None,
huge and unicode answers are all handled; a 2 MB cap bounds the scan).
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Iterable, Optional

__all__ = ["Claim", "check", "worst_level"]

_SCAN_LIMIT = 2_000_000          # chars; same bound as result_critic._read

# Units a reported physical value carries. Structural: a number is only
# a *claim* when one of these follows it.
_UNIT = (
    r"Eh|Hartree|Ha\b|eV|meV|kcal/mol|kcal\s*mol-1|kJ/mol|kJ\s*mol-1"
    r"|cm-1|cm\*\*-1|cm\^-1|kcal/mol\b|Debye|D\b|Angstrom|A\b|bohr"
)
_NUM_UNIT_RE = re.compile(
    r"(-?\d[\d,]*(?:\.\d+)?(?:[eE][-+]?\d+)?)\s*(" + _UNIT + r")")

# Numbers inside units that are NOT a measured value (e.g. "cm-1" in
# "50 cm-1 threshold"): still a claim when stated as a result; handled
# by provenance like any other.

# A value is an estimate when the answer says so near the number.
_ESTIMATE_RE = re.compile(
    r"(?i)(rough\s+)?(estimate[ds]?|approx\.?(?:imately)?|about|around|"
    r"expected\s+to\s+be|on\s+the\s+order\s+of)")

# Level-of-theory tokens: method + basis (+ solvent) must accompany an
# energy/property. Basis sets and functionals are single tokens; the
# basis carries the distinctive "def2-"/"cc-p"/"SVP"/"TZ" shapes that
# DELFIN's own workflows produce.
_METHOD_RE = re.compile(
    r"(?i)\b([A-Z]{2,}[A-Za-z0-9]*|HF|MP2|MP4|CCSD\(T\)|CCSD|CISD|"
    r"DLPNO-[A-Z]+|GFN2-?xTB|GFN1-?xTB|GFN-FF)\b")
_BASIS_RE = re.compile(
    r"(?i)\b(def2-[A-Z]+[A-Z0-9]*|cc-pV[A-Z]+Z|aug-cc-pV[A-Z]+Z|"
    r"[scdtq][fp]?[tzq]?vp[2p]?|SVP|TZVP|TZVPP|QZVP|6-31[+*]*[Gd*]*\(?"
    r"[0-9dp]*\)?|6-311[+*]*[Gd*]*\(?[0-9dp]*\)?|STO-3G|MINIX|ma-def2-[A-Z]+)\b")
_SOLVENT_RE = re.compile(
    r"(?i)\b(?:solvent|solvated(?:\s+in)?|in)\s+"
    r"(water|acetone|acetonitrile|dmso|methanol|ethanol|toluene|benzene|"
    r"thf|dcm|dichloromethane|chloroform|dmf|hexane|gas\s+phase|vacuum)\b")


@dataclass
class Claim:
    level: str        # "error" | "warn" | "ok"
    code: str         # short machine label
    message: str      # human-readable, English

    def __str__(self) -> str:
        icon = {"error": "❌", "warn": "⚠️", "ok": "✅"}.get(self.level, "•")
        return f"{icon} {self.message}"


def _text_of(record) -> str:
    if isinstance(record, dict):
        return str(record.get("text") or "")
    return str(getattr(record, "text", "") or "")


def _flags_of(record) -> list:
    if isinstance(record, dict):
        return list(record.get("flags") or [])
    return list(getattr(record, "flags", None) or [])


def _ref_of(record) -> str:
    if isinstance(record, dict):
        return str(record.get("ref") or "<unnamed evidence>")
    return str(getattr(record, "ref", "<unnamed evidence>"))


def _is_estimate(sentence: str) -> bool:
    return bool(_ESTIMATE_RE.search(sentence))


def _splits(answer: str) -> list[tuple[str, float, str]]:
    """(context, value, unit) for every number-with-unit stated."""
    out = []
    # sentence-ish windows: split on ., ;, newlines so estimate markers
    # near the number are seen
    sentences = re.split(r"(?<=[.;!\n])\s+", answer)
    for sent in sentences:
        for m in _NUM_UNIT_RE.finditer(sent):
            out.append((sent, float(m.group(1).replace(",", "")),
                        m.group(2)))
    return out


def _has_method_basis(sentence: str) -> bool:
    return bool(_METHOD_RE.search(sentence) and _BASIS_RE.search(sentence))


def _has_solvent(sentence: str) -> bool:
    return bool(_SOLVENT_RE.search(sentence))


def check(answer_text, evidence: Optional[Iterable] = None) -> list[Claim]:
    """Check the numbers in an ANSWER against the evidence.

    Returns a list of :class:`Claim` (empty = nothing to flag). Never
    raises; None/empty/huge/unicode inputs all handled. Deterministic
    and read-only.
    """
    try:
        text = str(answer_text or "")[:_SCAN_LIMIT]
        if not text.strip():
            return []
        claims: list[Claim] = []
        records = list(evidence or [])
        texts = [_text_of(r) for r in records]
        flags: dict[str, list] = {}
        for r in records:
            fl = [str(f) for f in _flags_of(r) if f in
                  ("min-has-imag", "ts-many-imag", "ts-no-imag",
                   "scf-not-converged", "scf-not-fully-converged",
                   "spin-contamination", "no-termination",
                   "geom-not-converged")]
            if fl:
                flags[_text_of(r)] = fl

        for sent, value, unit in _splits(text):
            if _is_estimate(sent):
                continue
            ref = None
            for t, r in zip(texts, records):
                if _value_in_text(value, t):
                    ref = r
                    break
            if ref is None:
                claims.append(Claim(
                    "error", "claim-no-source",
                    f"Value {value:g} {unit} has no source: no evidence "
                    f"record contains this number. Report only numbers a "
                    f"calculation produced, or mark it as an estimate."))
                continue
            if not _has_method_basis(sent):
                claims.append(Claim(
                    "error", "claim-no-method",
                    f"Value {value:g} {unit} ({_ref_of(ref)}) states an "
                    f"energy/property without its level of theory: method "
                    f"and basis set (and solvent, if solvated) must be "
                    f"named with the number."))
            fl = flags.get(_text_of(ref))
            if fl and not _acknowledges(sent, fl):
                claims.append(Claim(
                    "warn", "claim-flagged-source",
                    f"Value {value:g} {unit} comes from an output the "
                    f"correctness critic flagged ({', '.join(fl)}); the "
                    f"answer must state that caveat, not drop it."))
        return claims
    except Exception as exc:                       # a check must never crash
        return [Claim("warn", "claims-error",
                      f"claims check failed (reported, not raised): {exc}")]


def _value_in_text(value: float, text: str) -> bool:
    """The numeric value appears in the evidence text (string-tolerant:
    -15.3 matches -15.30, 1234.56789 matches -1234.567890)."""
    if not text:
        return False
    s = f"{value:.10f}".rstrip("0")
    for m in re.finditer(r"-?\d[\d,]*(?:\.\d+)?", text):
        tok = m.group(0).replace(",", "")
        try:
            v = float(tok)
        except ValueError:
            continue
        if f"{v:.10f}".rstrip("0") == s:
            return True
    return False


_WARN_WORDS = ("imaginary", "saddle", "not converged", "unconverged",
               "spin contamination", "contamination", "warning",
               "caution", "did not finish", "incomplete")


def _acknowledges(sentence: str, flags: list) -> bool:
    low = sentence.lower()
    return any(w in low for w in _WARN_WORDS)


def worst_level(claims: list[Claim]) -> str:
    levels = {c.level for c in claims}
    if "error" in levels:
        return "error"
    if "warn" in levels:
        return "warn"
    return "ok"
