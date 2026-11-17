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



def nonexistent_citations(text, *, repo_root, observed, named_in="",
                          acted_on=()):
    """The cited paths in ``text`` that the repo does not contain, as
    location claims the engine's guard already knows how to correct.

    The engine's location scanner only flags ``file_line`` claims; a bare
    path (no line number) is invisible to it. This wraps the code-claim
    scanner so the engine's single correction turn also demands proof for
    a bare nonexistent path -- the case the dashboard used to check on its
    own -- and returns them in ``LocationClaimFlag`` shape so the existing
    re-scan and caveat read them without a second code path.

    ``acted_on`` is the turn's executed-command ledger (the same evidence
    the functional-claim guard trusts): a path the turn RAN is not an
    ungrounded citation -- "python solver.py" grounds a claim about
    solver.py even though the file sits outside the scanned tree. Commands
    that REMOVE a path never ground anything: for "rm -rf solver.py" the
    turn's own action destroyed the file, so the conservative outcome is
    one extra correction turn, never a silently accepted fabrication.
    """
    from delfin.agent import verify_guard as _vg
    acted = _acted_tokens(acted_on)
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
        if f.path in acted:
            continue
        ref = f"{f.path}:{f.line}" if f.line else f.path
        out.append(_vg.LocationClaimFlag(
            claim=ref, path=f.path, line=f.line, kind="file_line"))
    return out


#: Tokens that DESTROY a path; a command containing one beside the token
#: never grounds a citation of it.
_DESTROY_TOKEN_RE = __import__("re").compile(
    r"(?i)(?:^|\s|-)(?:rm|rmdir|unlink|shred|git\s+(?:rm|clean)|mv)\b")


def _acted_tokens(commands):
    """The path tokens (full path and basename, word-bounded) that appear
    in an executed command, minus those the command destroyed.

    The ledger entries are already shaped like "bash python solver.py"
    (verify_guard.extract_exec_command) and only contain commands whose
    result showed they ran -- the same evidence the functional guard
    trusts. Matching is a whole-token match, not a free substring: a path
    named in prose or in a longer identifier does not ground a citation,
    while a quoted or bare argument does. Removing verbs disqualify the
    command (see nonexistent_citations).
    """
    import re as _re
    out = set()
    for cmd in (commands or ()):
        text = str(cmd or "")
        if not text.strip():
            continue
        destroyed = bool(_DESTROY_TOKEN_RE.search(text))
        for tok in _re.findall(r"[A-Za-z0-9_./\\-]+", text):
            tok = tok.strip("./\\")
            if not tok or "." not in tok.rsplit("/", 1)[-1]:
                # The cited object is a file: only tokens that name a file
                # (an extension in the last segment) can ground one.
                continue
            if destroyed:
                continue
            out.add(tok)
            out.add(tok.rsplit("/", 1)[-1])
    return out
