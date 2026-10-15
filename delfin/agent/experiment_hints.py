"""Whether a user message calls for the ``experiment-discipline`` skill.

Input: the text of the message that opens a turn. Output of
:func:`user_hint`: one line pointing at the skill, or "".

Semantics. The pointer is set for a message that asks for a measured
comparison or acts on one: an A/B run, a benchmark, an ablation, a
baseline or control run, a pre-registration, "is it better / faster
now", "did it help", comparing runs or arms or before/after, an
optimisation loop, landing a measured change, or making a test, CI or
benchmark pass. A chemistry optimisation ("optimise the geometry") and
ordinary questions get nothing. English and German, the two languages
users write in (matchers are the one place German appears in code).

Cost (len/4, measured): the skill's listing entry, about 22 tokens per
request in a code session; the pointer, 71 tokens in the tail of a turn
that matched; the skill itself, about 1,010 tokens, only when called.

Rejected: putting the rules in the stable prompt head (about 1,010
tokens in every request for a situation most turns never meet); matching
"optimi[sz]e" alone (a geometry optimisation is the commonest chemistry
request and is no optimisation loop).
"""

from __future__ import annotations

import re

SKILL_NAME = "experiment-discipline"

_PATTERNS = (
    r"\bA\s*/\s*B\b(?!\s*[-+*/=])",
    r"\bablation\w*",
    r"\bbenchmark\w*",
    r"\b(?:baseline|control)\s+(?:run|arm|vs\.?|versus|gegen)\b",
    r"\bkontroll(?:lauf|läufe|arm)\w*",
    r"\bnull[- ]?(?:test|run|lauf)\w*",
    r"\bpre-?regist\w*|\bpräregistr\w*|\bvorregistr\w*",
    r"\b(?:is|was|ist|war)\s+(?:it|this|that|es|das)\s+(?:now\s+|jetzt\s+|"
    r"really\s+|wirklich\s+)?(?:better|faster|improved|besser|schneller)\b",
    r"\b(?:better|faster|besser|schneller)\s+(?:now|jetzt|than before|"
    r"als vorher|geworden)\b",
    r"\b(?:did|does|has|hat|hilft|half)\b[^.?!\n]{0,40}"
    r"\b(?:help|helped|improve\w*|geholfen|verbessert)\b",
    r"\b(?:compare|vergleich\w*)\b[^.?!\n]{0,50}"
    r"\b(?:runs?|arms?|läufe|lauf|before|after|vorher|nachher|baseline|"
    r"branch(?:es)?|commits?)\b",
    r"\boptimi[sz]ation\s+loop\b|\boptimi[sz]ing\s+(?:the\s+)?(?:metric|"
    r"score|rate)\b|\boptimierungs(?:schleife|lauf)\w*",
    r"\b(?:land|landen)\b[^.?!\n]{0,30}\b(?:change|fix|branch|änderung|"
    r"fix|zweig)\w*",
    r"\b(?:make|get)\s+(?:the\s+)?(?:failing\s+)?(?:tests?|ci|checks?|"
    r"benchmark|suite)\s+(?:to\s+)?(?:pass|green)\b",
    r"\b(?:tests?|ci|suite|checks?)\s+(?:wieder\s+)?(?:grün|green)\s+"
    r"(?:machen|bekommen|kriegen)\b",
    r"\b(?:metric|score|rate|kennzahl|metrik)\s+(?:improv\w*|verbesser\w*|"
    r"steigern|raise|erhöhen)\b",
)
_RX = re.compile("|".join(f"(?:{p})" for p in _PATTERNS), re.IGNORECASE)


def matches(text: str) -> bool:
    """Whether *text* asks for, or acts on, a measured comparison."""
    return bool(_RX.search(str(text or "")))


def user_hint(text: str) -> str:
    """The one-line pointer for a message that calls for the skill."""
    if not matches(text):
        return ""
    return (f"[{SKILL_NAME}] This turn measures, compares or lands a "
            f"change. Before the first run, call skill(name='{SKILL_NAME}') "
            "and follow it: pre-register, one instrument state, the same "
            "population, never-worse cases listed, and only the user's own "
            "message approves a landing.")


__all__ = ["SKILL_NAME", "matches", "user_hint"]
