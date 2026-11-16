"""Engine-side verification decision for a finished answer.

The dashboard used to scan a finished answer for citations of paths that
do not exist, and force one correction turn when it found one; the
terminal had no such rerun for a bare (line-less) path citation. The
decision now lives here, and the engine's claim-grounding guard consults
it, so every surface reruns the same way.

The two note helpers are moved verbatim out of delfin/dashboard/
tab_agent.py so the dashboard keeps rendering the same line; the old
names remain importable there as one-line delegates.
"""
from __future__ import annotations

def self_verification_note(verified: bool) -> str:
    """The line shown after a forced self-correction turn.

    Two outcomes, two lines. The green check is a VERDICT — it may only
    appear when a re-scan of the corrected answer finds the claims
    grounded. Spending the correction budget is not the same event: a
    correction that raised a new claim, or that only rewrote the wording
    into hedges, is still an answer nobody checked, and printing the
    check for it told the reader the opposite of what happened.
    """
    if verified:
        return "✓ Self-verification: statements checked, answer corrected."
    return (
        "⚠️ Self-verification: the answer was revised, but the claims in it "
        "are still not backed by anything looked up this turn — treat them "
        "as unconfirmed."
    )



def correction_verdict_note(scan, text, *, keep=None, **kwargs) -> str:
    """Re-scan a corrected answer and return the note that matches the
    verdict, never the fact that a retry ran.

    ``keep`` filters the flags that count as "still ungrounded" (the code
    citation scan also returns soft, non-forcing flags). A scanner that
    raises counts as NOT verified: the check must be earned, and an
    unreadable answer has not earned it.
    """
    try:
        flags = list(scan(text, **kwargs) or ())
    except Exception:
        flags = ["rescan-failed"]
    if keep is not None:
        try:
            flags = [f for f in flags if keep(f)]
        except Exception:
            pass
    return self_verification_note(not flags)



def nonexistent_citations(text, *, repo_root, observed, named_in=""):
    """The cited paths in ``text`` that the repo does not contain, as
    location claims the engine's guard already knows how to correct.

    The engine's location scanner only flags ``file_line`` claims; a bare
    path (no line number) is invisible to it. This wraps the code-claim
    scanner so the engine's single correction turn also demands proof for
    a bare nonexistent path -- the case the dashboard used to check on its
    own -- and returns them in ``LocationClaimFlag`` shape so the existing
    feedback, re-scan and caveat read them without a second code path.
    """
    from delfin.agent import verify_guard as _vg
    try:
        flags = _vg.scan_for_ungrounded_code_claims(
            text, repo_root=repo_root,
            observed_files=set(observed or ()),
            named_in=named_in or "")
    except Exception:
        return []
    out = []
    for f in flags:
        if getattr(f, "kind", None) != "nonexistent":
            continue
        ref = f"{f.path}:{f.line}" if f.line else f.path
        out.append(_vg.LocationClaimFlag(
            claim=ref, path=f.path, line=f.line, kind="file_line"))
    return out
