"""Verification claims in an answer, and whether the turn's runs back them.

Input: the answer text, and the turn's runs as ``Run`` records (what
executed, how it ended, whether it ran after the last code edit and
touched the edited file, whether the pass followed an edit to the failing
test). Output: ``find_claims(text)`` lists the claims; ``judge(claims,
runs, edited=...)`` returns the unbacked ones with the reason.

Two kinds of claim:

* ``test``: a statement about a test result -- "tests pass", "all green",
  "5 passed", "alle Tests bestanden". Backed by an executed, passing test
  run; when the turn edited code, one that ran after the last edit and
  touched the edited file. A figure ("5 passed") must equal a passed
  count the runs reported, when any was reported.
* ``behaviour``: "verified", "tested", "works", "funktioniert",
  "verifiziert", "getestet", a check mark. Checked only in a turn that
  edited code, where it is a statement about that code: backed by an
  executed, passing test, check or run of it after the last edit. In a
  turn without an edit "verified" is usually a reading (file:line), which
  the location scanner judges; it is out of scope here.

A pass obtained by editing the failing test (``Run.gamed``) backs
nothing. Sentence level; a sentence that negates, conditions, plans, asks
or reports someone else's claim is not a claim (``_WAIVE_RE``). German
only on the recognising side: the matchers read what a model or user
writes (framework rule), the code stays English.

Measured on 1,054 captured benchmark answers (2026-10-13/14) before it
was wired anywhere; the numbers are in the commit that introduced it.
Pure; never raises.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Iterable, Optional

_TEST_CLAIM_RE = re.compile(
    r"(?i)(?:\btests?\b|\btestsuite\b|\bsuite\b|\bspecs?\b)[^.\n]{0,40}?"
    r"\b(?:pass(?:es|ed)?|passing|green|grün|bestanden|laufen durch|"
    r"läuft durch|succeed(?:s|ed)?|erfolgreich)\b"
    r"|\b[1-9]\d*\s+(?:passed|tests?\s+(?:pass(?:ed)?|bestanden|grün|green))\b"
    r"|\b[1-9]\d*/[1-9]\d*\s+(?:tests?\s+)?(?:pass(?:ed)?|grün|green|bestanden)\b"
    r"|\balle\b[^.\n]{0,30}\b(?:grün|bestanden)\b"
    r"|\ball\b[^.\n]{0,30}\b(?:pass(?:ed)?|green)\b"
    r"|\b(?:grün|green)\s+(?:sind|ist|are|is)?\s*(?:alle|all)?\s*tests?\b")
# "works" / "funktioniert" only as a statement about the result of the
# turn ("works now", "funktioniert jetzt"); bare, they describe code read
# ("der Fall funktioniert zufällig, weil ...").
_BEHAVIOUR_CLAIM_RE = re.compile(
    r"(?i)\b(?:verified|verifiziert|getestet|tested)\b"
    r"|\b(?:works|funktioniert|läuft)\s+(?:jetzt|nun|wieder|now|again|"
    r"korrekt|correctly|einwandfrei|fehlerfrei|as expected|wie erwartet)\b"
    r"|✅|✔")
_COUNT_RE = re.compile(
    r"(?i)\b(\d+)\s+(?:passed|tests?\s+(?:pass(?:ed)?|bestanden|grün|green))\b"
    r"|\b(\d+)/\1\s+(?:tests?\s+)?(?:pass(?:ed)?|grün|green|bestanden)\b"
    r"|\balle\s+(\d+)\s+tests?\b|\ball\s+(\d+)\s+tests?\b")
# A sentence that does not assert a result of this turn: a negation, a
# condition, a plan or a request, a question, a quoted claim (the commit
# message says ..., laut ...), the harness's own [verify] lines.
_WAIVE_RE = re.compile(
    r"(?i)\bnicht\b|\bnot\b|\bkein\w*|\bno\b|\bnever\b|\bnie\b|\bweder\b"
    r"|\bcannot\b|\bcan't\b|\bcould not\b|\bkonnte\b|\bunverified\b"
    r"|\buntested\b|\bungetestet\b|\bunverifiziert\b|\bpending\b"
    r"|\bblock(?:ed|iert)\b|\brefused\b|\babgelehnt\b"
    r"|\bif\b|\bwenn\b|\bfalls\b|\bsobald\b|\bonce\b|\bunless\b"
    r"|\bshould\b|\bsoll(?:te|en)?\b|\bmust\b|\bmüss\w*|\bwould\b|\bwürde\w*"
    r"|\bwill\b|\bwerde\w*|\bneed\w*\b|\bbrauch\w*|\bplease\b|\bbitte\b"
    r"|\bto verify\b|\bzu (?:prüfen|testen|verifizieren)\b|\bum\b[^.]{0,60}\bzu\b"
    r"|\bexpect\w*|\berwart\w*|\bclaims?\b|\bclaimed\b|\bbehaupt\w*"
    r"|\blaut\b|\baccording to\b|\bcommit(?:-| )?message\b|\bgemeldet\b"
    r"|\breported\b|\bbisher\b|\bpreviously\b|\bbefore\b|\bvorher\b|\bzuvor\b"
    r"|\?|\[verify\]|\bnext step\b|\bnächste\w* schritt"
    # A plan step or an acceptance criterion: what is to hold, not what held.
    r"|schritt|\bsteps?\b|akzeptanz|acceptance|\bAC\b|kriteri|criteri"
    r"|\bgate\b|\bhalten\b|\bkeep\b|\bziel\b|\bgoal\b|\bplan"
    r"|\bverifikation:|\bverification:"
    # Whether, a transition, a quoted claim, a historic or a red result.
    r"|\bob\b|\bwhether\b|grün auf rot|rot auf grün|green to red|red to green"
    r"|„|“|”|\bpre-?fix\b|\bvor dem fix\b"
    r"|\b[1-9]\d*\s+failed\b|\brot\b|\bred\b"
    # Nothing ran: "0 passed", "5 skipped", "none of the tests ran".
    r"|\b0\s+passed\b|\bskipped\b|übersprungen|\bnone\b"
    # Distributive: a procedure ("per module: ... tests green").
    r"|\bpro (?:modul|datei|commit)\b|\bper (?:module|file|commit)\b"
    r"|\bjede[nmrs]?\b|\beach\b|\bevery\b"
    # What a run WOULD show: "der Nachweis dafür, dass cargo test grün ist".
    r"|nachweis|\bproof that\b|\bevidence that\b"
    # Checked by reading: a file:line citation, "Fundstellen", "Zeilen";
    # the location scanner judges those against what was read.
    r"|\b[\w./-]+\.\w{1,5}:\d+|fundstelle|zeilenangabe|\bcitation|\bline \d")
# An invariant ("stay green", "bleiben grün"): in a turn that ran and edited
# nothing it describes a goal; in a turn that did, it asserts the state.
_INVARIANT_RE = re.compile(r"(?i)\bbleib\w*|\bremain\w*|\bstay\w*")


@dataclass
class Claim:
    kind: str                      # "test" | "behaviour"
    sentence: str
    count: Optional[int] = None    # the N of "N passed", when stated
    invariant: bool = False        # "stays green": judged only with evidence


@dataclass
class Run:
    """One executed run of the turn, as a claim needs it.

    ``fresh``: no later edit to a file that existed before the turn (or to
    a test file, or an ``__init__``) falls inside this run's scope -- the
    run still describes the code. ``covers_last_edit``: it ran after the
    last code edit of the turn and touched that file."""
    kind: str                      # "test" | "check" | "run"
    outcome: str                   # verify_commands.run_outcome
    fresh: bool = True
    covers_last_edit: bool = True
    gamed: bool = False
    passed_count: Optional[int] = None


@dataclass
class Verdict:
    claim: Claim
    reason: str
    extra: dict = field(default_factory=dict)


# Abbreviations whose dot does not end a sentence.
_ABBREV_RE = re.compile(
    r"(?i)\b(?:inkl|bzw|ca|vgl|ggf|usw|etc|evtl|zzgl|e\.g|i\.e|z\.b|u\.a|d\.h"
    r"|vs|approx|incl|resp)\.")


def _sentences(text: str) -> list[str]:
    body = re.sub(r"```.*?```", " ", str(text or ""), flags=re.S)
    lines = [ln for ln in body.splitlines()
             if not ln.lstrip().startswith((">", "|"))]
    out: list[str] = []
    for ln in lines:
        ln = _ABBREV_RE.sub(lambda m: m.group(0).replace(".", "\x00"), ln)
        out.extend(s.replace("\x00", ".")
                   for s in re.split(r"(?<=[.!])\s+", ln) if s.strip())
    return out


def find_claims(text: str) -> list[Claim]:
    """Every verification claim in ``text`` (module docstring)."""
    out: list[Claim] = []
    for s in _sentences(text):
        if _WAIVE_RE.search(s):
            continue
        if _TEST_CLAIM_RE.search(s):
            m = _COUNT_RE.search(s)
            n = None
            if m:
                n = int(next(g for g in m.groups() if g))
            out.append(Claim("test", s.strip()[:300], n,
                             bool(_INVARIANT_RE.search(s))))
        elif _BEHAVIOUR_CLAIM_RE.search(s):
            out.append(Claim("behaviour", s.strip()[:300]))
    return out


def judge(claims: Iterable[Claim], runs: Iterable[Run], *,
          edited: bool, tried_tests: bool = True) -> list[Verdict]:
    """The claims the runs do not back, each with its reason.

    ``edited``: the turn edited code (a behaviour claim is judged only
    then). ``tried_tests``: the turn ran or tried a test; a test claim in
    a turn that did neither and edited nothing restates an earlier turn's
    result ("tests are green (10/10)" in a summary) and is not judged. A test claim needs a passing, fresh, non-gamed test run; a
    behaviour claim needs, after the last code edit and touching it, a
    passing test or check, or a program run that executed (its exit code
    is the program's own answer -- a search that exits 1 on no match ran).
    Reasons: ``gamed`` (only a pass after an edit to the failing test),
    ``stale`` (the passing run predates an edit inside its scope),
    ``no_passing_run``, ``no_covering_run``, ``count_mismatch``."""
    runs = list(runs)
    counts = {r.passed_count for r in runs if r.passed_count}
    out: list[Verdict] = []
    for c in claims:
        if c.invariant and not runs and not edited:
            continue
        if c.kind == "test" and not edited and not tried_tests:
            continue
        if c.kind == "test":
            passing = [r for r in runs if r.kind == "test"
                       and r.outcome == "passed"]
            good = [r for r in passing if r.fresh and not r.gamed]
            if not good:
                if any(r.gamed for r in passing):
                    reason = "gamed"
                elif passing:
                    reason = "stale"
                else:
                    reason = "no_passing_run"
                out.append(Verdict(c, reason))
            elif c.count is not None and counts and c.count not in counts:
                out.append(Verdict(c, "count_mismatch",
                                   {"observed": sorted(counts)}))
        elif edited:
            good = [r for r in runs if r.covers_last_edit and not r.gamed
                    and ((r.kind in ("test", "check")
                          and r.outcome == "passed")
                         or (r.kind == "run"
                             and r.outcome in ("passed", "failed")))]
            if not good:
                gamed = any(r.gamed and r.covers_last_edit for r in runs)
                out.append(Verdict(c, "gamed" if gamed else "no_covering_run"))
    return out


_REASON_TEXT = {
    "no_covering_run": "no test, check or run after the last edit touched "
                       "the edited code",
    "no_passing_run": "no test run in this turn passed",
    "stale": "the passing run came before a later edit inside its scope",
    "gamed": "the only passing run came after an edit to the failing test "
             "itself",
    "count_mismatch": "the runs reported a different number of passed tests",
}


def caveat(verdicts: list) -> str:
    """The visible caveat for unbacked verification claims ("" for none).
    One line per claim, the claim quoted short, the reason in words."""
    if not verdicts:
        return ""
    items = []
    for v in verdicts[:4]:
        quote = re.sub(r"[*`_#>]+", "", v.claim.sentence).strip()
        if len(quote) > 90:
            quote = quote[:87] + "..."
        reason = _REASON_TEXT.get(v.reason, v.reason)
        if v.reason == "count_mismatch" and v.extra.get("observed"):
            reason += f" ({', '.join(str(n) for n in v.extra['observed'])})"
        items.append(f"'{quote}' -- {reason}")
    return ("\n\n[verify] Caveat: these verification claims are not "
            "backed by a run in this turn: " + "; ".join(items)
            + ". Please treat them as unconfirmed.")


# The answer already says the change was not checked.
_DISCLOSED_RE = re.compile(
    r"(?i)\b(?:untested|ungetestet|unverified|unverifiziert|not (?:yet )?"
    r"(?:tested|verified|run|checked)|nicht (?:getestet|verifiziert|geprüft|"
    r"ausgeführt|gelaufen)|could not (?:run|test|verify)|konnte [^.\n]{0,60}"
    r"nicht (?:ausführen|testen|prüfen|laufen))\b|\bUNTESTED\b")


def discloses_unverified(text: str) -> bool:
    """True when the answer names its change as not checked."""
    return bool(_DISCLOSED_RE.search(str(text or "")))


__all__ = ["Claim", "Run", "Verdict", "find_claims", "judge",
           "discloses_unverified", "caveat"]
