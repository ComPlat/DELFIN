"""Claims that a file was changed, in a turn where no change landed.

Input: the answer text and how many file changes of the turn landed and
how many were refused or failed (``TurnVerification``). Output:
``find(text)`` lists the sentences that report a change as done;
``caveat(claims, landed=..., failed=...)`` is the visible ``[verify]``
line, or "".

Measured 2026-10-19 (session 78c36f8c, gpt-5-nano, a 3D game in one HTML
file): for fifty minutes every write was refused -- apply_patch for a
missing ``diff``, edit_file for a file server "using" the file -- and the
model answered "Ich habe die zentrale Mittellinie aus dem Rendering
entfernt", "Was umgesetzt wurde: ...", turn after turn. The user tested,
saw nothing change, and was told the change was there. No caveat named
it: the verification-claim scanner judges "tested" and "tests pass", not
"I changed it".

Judged only when nothing landed in the turn: a turn that changed files
may summarise them however it likes. Sentences that plan, ask, negate or
quote are waived (``verification_claims._WAIVE_RE``). German and English
on the recognising side; the code stays English. Pure; never raises.
"""

from __future__ import annotations

import re

from . import verification_claims as _vc

# First person, perfect: "ich habe ... entfernt", "I've added", "I have
# updated"; passive result: "wurde ... entfernt", "has been fixed"; and the
# summary header a model writes over a list of done changes.
_DONE_VERBS_DE = (
    r"(?:geändert|entfernt|hinzugefügt|umgesetzt|implementiert|eingebaut|"
    r"angepasst|ersetzt|behoben|gefixt|gepatcht|aktualisiert|erweitert|"
    r"verbreitert|vergrößert|gelöscht|eingefügt|korrigiert|überarbeitet|"
    r"umgebaut|ergänzt|gesetzt)")
_DONE_VERBS_EN = (
    r"(?:changed|removed|added|implemented|updated|patched|fixed|replaced|"
    r"modified|rewrote|rewritten|refactored|deleted|inserted|widened|"
    r"adjusted|applied)")
_CHANGE_CLAIM_RE = re.compile(
    r"(?i)\b(?:ich\s+habe|hab(?:e)?\s+ich)\b[^.\n]{0,120}?\b" + _DONE_VERBS_DE
    + r"\b"
    r"|\b(?:wurde|wurden|ist|sind)\b[^.\n]{0,80}?\b" + _DONE_VERBS_DE
    + r"\b(?:\s+worden)?"
    r"|\bI(?:'ve|\s+have)?\s+(?:now\s+|just\s+|also\s+)?" + _DONE_VERBS_EN
    + r"\b"
    r"|\b(?:has|have)\s+been\s+" + _DONE_VERBS_EN + r"\b"
    r"|\bwas\s+ich\s+(?:umgesetzt|geändert|gemacht)\s+habe\b"
    r"|\bwas\s+umgesetzt\s+wurde\b|\bumgesetzt\s*:")


def find(text: str) -> list[str]:
    """The sentences of *text* that report a file change as done."""
    out: list[str] = []
    try:
        for s in _vc._sentences(text):
            if _vc._WAIVE_RE.search(s):
                continue
            if _CHANGE_CLAIM_RE.search(s):
                out.append(s.strip()[:300])
    except Exception:
        return []
    return out


def caveat(claims: list[str], *, landed: int, failed: int) -> str:
    """The ``[verify]`` line for change claims no landed change backs."""
    if not claims or landed:
        return ""
    quote = re.sub(r"[*`_#>]+", "", claims[0]).strip()
    if len(quote) > 90:
        quote = quote[:87] + "..."
    why = (f"{failed} write{'s' if failed != 1 else ''} in this turn "
           "failed or were refused" if failed
           else "this turn made no file change")
    return ("\n\n[verify] Caveat: this answer reports a change as made "
            f"('{quote}'), but no file change landed -- {why}. The files "
            "are as they were before this turn.")


__all__ = ["find", "caveat"]
