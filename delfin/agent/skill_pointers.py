"""Which skill a piece of text calls for, from the skills' own metadata.

Input: a user message or a tool result, and the skill registry. Output:
at most two pointers naming a skill and the word that called for it.

THE DEFECT THIS CLOSES. 67 skills are installed and the agent is told
about none of them. The ``skill`` tool's whole schema is 107 tokens --
"Skill name, no slash", no enum -- and 6 of the 67 names appear in the
51k-character system prompt only as coincidental matches of common words
(``loops``, ``merge``, ``pr``, ``push``). ``pipeline-build``, which says
in as many words "work only through the MCP server delfin-tools; never
hand-roll pipeline logic", is not in it.

Measured consequence, on a real 33-minute run: the agent called
``write_file`` 20 times and ``bash`` 46 times to assemble and start ORCA
by hand, 12 of its 14 calculations were abandoned, and the 26 tools of
the ``delfin-tools`` server -- offered on that request -- were called
zero times. The field's two failure modes are exactly these ends: full
documentation floods the context, and no documentation produces invented
calls. This repository sat at the second end.

WHY A POINTER AND NOT A DIRECTORY. Measured here: names plus the
descriptions the registry already resolves cost 1,461 tokens, 12% of a
prompt head whose budget has no slack; names alone cost about 370. A
pointer costs about 60 tokens WHEN IT FIRES and nothing otherwise, and
it answers the question the directory answers badly -- not "what exists"
but "what applies to this". ``repo_error_hints`` already does exactly
this for one skill, section by section, so this is that mechanism
generalised rather than a second one invented beside it. The two differ
in granularity on purpose: that module points at a SECTION of
git-errors, this one points at a skill.

THE TRIGGERS ARE THE METADATA. No trigger table: a second list for 67
skills is a list that drifts from the first, which is the defect this
repository keeps producing. A skill's name and description ARE its
trigger, which also means a skill nobody can find has a DESCRIPTION
problem and is fixed by writing a better one -- not by adding a pattern
here. Ten skills currently carry fewer than two distinctive terms; a
test names them so the gap is visible and reducible.

WHAT COUNTS AS DISTINCTIVE is computed, not listed. A term in many
skills' metadata ("run", "method", "report", "delfin") cannot select
between them, so the cut is on document frequency over the registry
itself. A hand-written stopword list would be the same drift in
miniature.
"""
from __future__ import annotations

import re
from collections import Counter
from typing import Any, Iterable, Optional

#: A term occurring in more than this many skills' metadata cannot
#: select between them. Three of 67 is about 4%: `casscf` (1) selects,
#: `report` (many) does not. Computed against the registry, so adding
#: skills re-tunes it without an edit here.
MAX_SKILLS_PER_TERM = 3

#: Terms shorter than this are noise in both directions.
MIN_TERM_LEN = 4

#: At most two pointers, as repo_error_hints allows at most two sections
#: per result: a turn that is handed five playbooks reads none.
MAX_POINTERS = 2

#: A match needs this many distinctive terms. One is a coincidence --
#: `pipeline` occurs in prose that has nothing to do with assembling
#: one -- and requiring two was measured to be the difference between
#: pointing at `pipeline-build` for "set up a TD-DFT pipeline" and
#: pointing at it for "the CI pipeline is red".
MIN_TERMS_FOR_MATCH = 2

#: Only the tail of a long text is scanned, as repo_error_hints does:
#: the relevant words are near the end of a message or a result.
SCAN_TAIL_CHARS = 4000

_WORD = re.compile(r"[a-z][a-z0-9]{2,}")


def _words(text: Optional[str]) -> set[str]:
    return set(_WORD.findall(str(text or "").lower()))


def _metadata_words(skill: Any) -> set[str]:
    """The terms a skill's own metadata offers as triggers."""
    name = str(getattr(skill, "name", "") or "").replace("-", " ")
    desc = str(getattr(skill, "description", "") or "")
    return {w for w in (_words(name) | _words(desc)) if len(w) >= MIN_TERM_LEN}


def distinctive_terms(skills: Iterable[Any]) -> dict[str, set[str]]:
    """Per skill, the terms that select IT rather than the registry.

    Document frequency over the registry decides: a term in more than
    MAX_SKILLS_PER_TERM skills describes the project, not the skill.
    """
    per: dict[str, set[str]] = {}
    for s in skills:
        name = str(getattr(s, "name", "") or "")
        if name:
            per[name] = _metadata_words(s)
    df: Counter = Counter()
    for words in per.values():
        for w in words:
            df[w] += 1
    return {n: {w for w in words if df[w] <= MAX_SKILLS_PER_TERM}
            for n, words in per.items()}


def unfindable(skills: Iterable[Any]) -> list[str]:
    """Skills whose metadata cannot select them, in name order.

    Reported rather than patched around: a skill here needs a better
    DESCRIPTION, and adding a trigger pattern for it would start the
    second list this module exists to avoid.
    """
    terms = distinctive_terms(skills)
    return sorted(n for n, w in terms.items() if len(w) < MIN_TERMS_FOR_MATCH)


def match(text: str, skills: Iterable[Any], *,
          limit: int = MAX_POINTERS) -> list[tuple[str, tuple[str, ...]]]:
    """``[(skill name, the terms that matched)]``, best first.

    Empty for text that calls for nothing, which is most text. Ranked by
    how many distinctive terms matched, then by name so the order is
    stable across sessions -- a pointer that moves between runs is a
    pointer nobody can test.
    """
    body = str(text or "")[-SCAN_TAIL_CHARS:]
    found = _words(body)
    if not found:
        return []
    # TWO TERMS, ALWAYS, and the plaster that was here is gone.
    #
    # A single term unique to one skill's metadata used to be treated as
    # decisive, because "TD-DFT vertical excitation energies" matched
    # `excitation` and nothing else. That was a workaround for broken
    # metadata, not a rule: three skills carried ">-" as their whole
    # description (a YAML folded block the line-based front-matter
    # reader took literally) and 28 more were described by their title
    # alone. With the block scalars folded and the `>` cue read, the
    # same sentence matches `energies`, `excitation` and `state` -- so
    # the rule is unnecessary, and it cost two false positives that no
    # frequency trick removes: "the CI pipeline is red" selected
    # `pipeline-build` and "thanks, that looks good" selected
    # `provenance`, both on one ambiguous English word.
    out: list[tuple[int, str, tuple[str, ...]]] = []
    for name, terms in distinctive_terms(skills).items():
        hit = tuple(sorted(terms & found))
        if len(hit) >= MIN_TERMS_FOR_MATCH:
            out.append((len(hit), name, hit))
    out.sort(key=lambda row: (-row[0], row[1]))
    return [(name, hit) for _n, name, hit in out[:max(0, int(limit))]]


def pointer(matches: list[tuple[str, tuple[str, ...]]]) -> str:
    """The text a turn carries, or "" for no match.

    Names the skill, the tool that loads it and the words that called
    for it. The words matter: a pointer that cannot be argued with is a
    pointer that gets ignored when it is wrong.
    """
    if not matches:
        return ""
    parts = [f"`{name}` ({', '.join(terms)})" for name, terms in matches]
    return ("A playbook covers this: " + "; ".join(parts)
            + ". `skill(name=...)` loads one; it is a procedure that has "
              "been measured here, not a suggestion.")


def hint(text: str, skills: Iterable[Any], *,
         shown: Optional[set] = None) -> str:
    """``pointer`` for *text*, each skill named at most once per session.

    ``shown`` is the session's set, passed in rather than kept here:
    this module holds no state, so a test cannot be decided by an
    earlier test's run.
    """
    found = match(text, skills)
    if shown is not None:
        found = [(n, t) for n, t in found if n not in shown]
    if not found:
        return ""
    if shown is not None:
        for name, _terms in found:
            shown.add(name)
    return pointer(found)
