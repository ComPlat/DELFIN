"""An answer that claims an improvement without a paired comparison.

Input: the model's finished answer and the session's execution ledger
(the commands that ran and produced output). Output of :func:`scan`:
the sentences that state a comparative effect as a measured figure --
"improved by 12 %", "3x faster", "the pass rate rose from 6/10 to 9/10",
"um 30 % schneller" -- when no paired comparison ran in the session;
of :func:`caveat`: the ``[verify]`` line naming them, or "".

Semantics.

* A claim is a comparative word (better, faster, improved, reduced, rose,
  fewer, ...; German forms) within the same sentence as a relative
  figure (a percentage, a factor such as ``3x`` / ``3-fold``, percentage
  points) or a "from X to Y" pair of figures.
* Not a claim: a sentence with an uncertainty marker (the shared hedge
  list), a conditional or a plan (would, could, should, expected, aim,
  target; German forms), a question, a figure attributed to a source
  (according to, the paper, the documentation, laut), fenced code,
  quoted lines and table rows.
* Evidence that silences it: a paired comparison in the session's ledger
  that exited cleanly -- ``delfin-agent compare-arms`` (a refused or
  blocked comparison exits 3 and is not in the ledger), ``delfin-agent
  bench compare`` / ``bench gate``, ``hyperfine`` with two commands,
  ``pyperf compare_to``.

Measured on 4,465 recorded answers (benchmark run excerpts, archived
and local sessions; read only): 9 firings, all in the answers of one
benchmark task that times two scripts with jitter and reports "~5 %
schneller" from unpaired runs -- the reading the task exists to catch.
Before table rows were exempt: 21, the other 12 cost columns of a
basis-set comparison and planning estimates written from general
knowledge.

The caveat names the claim and what a comparison needs; it does not
force a correction turn (a retry could only restate the figure).

Rejected: treating any two runs as a comparison (two runs of different
populations, or under a changed instrument, are exactly the readings
this exists to mark); a claim list in the prompt (a rule as a sentence
does not bind).

Limits: lexical. A comparison stated without a figure ("much faster")
is not judged; a figure in a table row without a comparative word is not
judged.
"""

from __future__ import annotations

import re
from typing import Iterable, Optional

_COMP = (
    r"(?:improv\w*|better|faster|slower|quicker|speed-?ups?|cheaper|"
    r"fewer|higher|lower|reduc\w*|increas\w*|decreas\w*|dropp\w*|rose|"
    r"rises|risen|gain\w*|boost\w*|halv\w*|doubl\w*|cut\s+by|"
    r"verbesser\w*|besser\w*|schneller\w*|langsamer\w*|weniger|"
    r"höher\w*|niedriger\w*|reduzier\w*|gesenkt|gesteigert|gestiegen|"
    r"gesunken|beschleunig\w*|verkürz\w*|halbiert|verdoppelt)")
_REL_FIG = (
    r"(?:[-+]?\d+(?:[.,]\d+)?\s*(?:%|percent\b|prozent\b|pp\b|"
    r"percentage\s+points?|prozentpunkte\w*)|\b\d+(?:[.,]\d+)?\s*(?:x|×)"
    r"(?![\w])|\b\d+(?:[.,]\d+)?[- ]?(?:fold|fach)\b|\b\d+(?:[.,]\d+)?\s+"
    r"times\b)")
_FROM_TO = (r"\b(?:from|von)\s+[-+]?\d[\d.,/%]*\s*(?:[a-zA-Z%/]+\s+)?"
            r"(?:to|auf|zu)\s+[-+]?\d")
_CLAIM_RE = re.compile(
    rf"{_COMP}[^.!?\n]{{0,60}}?{_REL_FIG}|{_REL_FIG}[^.!?\n]{{0,40}}?{_COMP}"
    rf"|{_COMP}[^.!?\n]{{0,60}}?{_FROM_TO}|{_FROM_TO}[^.!?\n]{{0,40}}?{_COMP}",
    re.IGNORECASE)
_CONDITIONAL_RE = re.compile(
    r"\b(?:would|could|might|may|should|expect\w*|aim\w*|target\w*|goal|"
    r"hope\w*|plan\w*|predict\w*|estimat\w*|if\b|unless|würde\w*|könnte\w*|"
    r"sollte\w*|soll\b|erwart\w*|ziel\w*|geschätzt|schätz\w*|falls|wenn)",
    re.IGNORECASE)
_SOURCE_RE = re.compile(
    r"\b(?:according\s+to|the\s+paper|the\s+authors?|the\s+docs?|"
    r"documentation|the\s+manual|reported\s+(?:by|in)|literature|"
    r"benchmark\s+(?:paper|study)|laut|gemäß|der\s+literatur|"
    r"im\s+handbuch|die\s+autoren)\b", re.IGNORECASE)
_PAIRED_RE = re.compile(
    r"\bdelfin-agent\s+(?:compare-arms|bench\s+(?:compare|gate))\b|"
    r"\bdelfin\.agent\.cli\s+(?:compare-arms|bench\s+(?:compare|gate))\b|"
    r"\bpyperf\s+compare_to\b")
# hyperfine options that take a value; every other word after the program
# is a command it times.
_HYPERFINE_VALUED = frozenset({
    "-w", "--warmup", "-m", "--min-runs", "-M", "--max-runs", "-r", "--runs",
    "-p", "--prepare", "-c", "--cleanup", "-P", "--parameter-scan",
    "-D", "--parameter-step-size", "-L", "--parameter-list", "-s", "--style",
    "-S", "--shell", "--export-json", "--export-csv", "--export-markdown",
    "--export-asciidoc", "--export-orgmode", "-u", "--time-unit",
    "--output", "--input", "-n", "--command-name", "--setup", "--sort"})


def _hyperfine_pair(cmd: str) -> bool:
    """hyperfine timing two or more commands: one interleaved run under
    one machine state, with its spread."""
    import shlex
    m = re.search(r"\bhyperfine\b(.*)", cmd)
    if not m:
        return False
    try:
        words = shlex.split(m.group(1))
    except ValueError:
        return False
    commands, i = 0, 0
    while i < len(words):
        w = words[i]
        if w in _HYPERFINE_VALUED:
            i += 2
            continue
        if not w.startswith("-"):
            commands += 1
        i += 1
    return commands >= 2


_FENCE_RE = re.compile(r"```.*?```", re.DOTALL)
_SENT_RE = re.compile(r"[^.!?\n]+(?:[.!?]+|$)", re.MULTILINE)
_MAX_FLAGS = 3


def _sentences(text: str) -> Iterable[str]:
    body = _FENCE_RE.sub(" ", text or "")
    for line in body.splitlines():
        # A quoted line, and a table row: measured on 4,465 recorded
        # answers, every table-row firing was a comparison of methods from
        # general knowledge (basis-set cost columns, planning estimates),
        # not an effect the session measured.
        if line.lstrip().startswith((">", "|")):
            continue
        for m in _SENT_RE.finditer(line):
            s = m.group(0).strip()
            if s:
                yield s


def paired_comparison_ran(exec_commands: Optional[Iterable[str]]) -> bool:
    """Whether the session ran a paired comparison."""
    return any(_PAIRED_RE.search(str(c)) or _hyperfine_pair(str(c))
               for c in (exec_commands or ()))


def scan(text: str, exec_commands: Optional[Iterable[str]] = None) -> list[str]:
    """The comparative claims of *text* that no paired comparison backs."""
    if not text or paired_comparison_ran(exec_commands):
        return []
    try:
        from .german import HEDGE_RE
    except Exception:                                   # pragma: no cover
        HEDGE_RE = None
    out: list[str] = []
    for s in _sentences(text):
        if not _CLAIM_RE.search(s):
            continue
        if s.endswith("?") or _CONDITIONAL_RE.search(s) or _SOURCE_RE.search(s):
            continue
        if HEDGE_RE is not None and HEDGE_RE.search(s):
            continue
        claim = " ".join(s.split())[:120]
        if claim not in out:
            out.append(claim)
        if len(out) >= _MAX_FLAGS:
            break
    return out


def caveat(claims: list[str]) -> str:
    if not claims:
        return ""
    items = "; ".join(f"'{c}'" for c in claims)
    return ("\n\n[verify] Caveat: comparative figure without a paired "
            f"comparison in this session: {items}. A comparison needs the "
            "same items under one instrument with the axis declared "
            "(delfin-agent compare-arms); until then treat the figure as "
            "unconfirmed.")


__all__ = ["scan", "caveat", "paired_comparison_ran"]
