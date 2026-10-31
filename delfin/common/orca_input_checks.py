"""ORCA simple-input checks DELFIN runs before an input reaches SLURM.

Consumed by ``delfin.api.validate_orca_input`` (delfin/api.py, the
function ``delfin.agent.calc_policy.check_orca_text`` calls before every
submit and the ops server exposes). Every finding is returned as a plain
``Issue`` here and mapped to ``ValidationIssue`` by the caller.

Sources of truth:
- keyword/block vocabulary: ``delfin.orca_knowledge`` (OPI 2.0.0 catalog
  + the words ORCA's own parser accepted),
- ``delfin.orca_knowledge.probe.probe_input`` when the static check
  cannot decide.
"""
from __future__ import annotations

import re
from dataclasses import dataclass


@dataclass
class Issue:  # noqa: D401 - mirrors api.ValidationIssue fields
    severity: str   # "error" | "warning" | "info"
    code: str
    message: str
    suggestion: str = ""


def _closest(words, word):
    """The known keyword closest to ``word`` (case-blind), else "".

    Both the word and its reverse are compared: invented keywords are
    often a real keyword's tokens in the wrong order (OptTight vs
    TightOpt), which plain distance misses.
    """
    from difflib import SequenceMatcher
    pool = sorted(words)
    w = word.lower()
    # Token-aware ranking: an invented keyword is often a real one's
    # tokens in the wrong order (OptTight vs TightOpt). Score each
    # catalog word by how well it covers BOTH the full word and the
    # word's best split into two parts, so a word made of the same
    # tokens ("tight" + "opt") beats a superficially closer string.
    def _split_score(k):
        best = 0.0
        for i in range(1, len(w)):
            h1, h2 = w[:i], w[i:]
            best = max(best, min(SequenceMatcher(None, h1, k).ratio(),
                                 SequenceMatcher(None, h2, k).ratio()))
        return best
    ranked = sorted(
        ((-(SequenceMatcher(None, w, k).ratio() + _split_score(k)), k)
         for k in pool))
    return ranked[0][1] if ranked and -ranked[0][0] >= 0.5 else ""


def check_keyword_line(text: str) -> list[Issue]:
    """Check every word on the ``!`` line(s) against ORCA's vocabulary."""
    from delfin.orca_knowledge import known_keyword
    issues: list[Issue] = []
    for line in re.findall(r"^\s*!(.*)$", text, re.MULTILINE):
        for word in line.split():
            w = word.strip()
            if not w:
                continue
            if w.startswith("%"):
                issues.append(Issue(
                    severity="error", code="block_on_keyword_line",
                    message=f"`{w}` is a %-block written on the `!` line.",
                    suggestion="Blocks belong on their own line: "
                               f"`%{w.lstrip('%')} ... end`.",
                ))
                continue
            if known_keyword(w):
                continue
            from delfin.orca_knowledge import _keywords
            kws, _solvents, aliases, _blocks = _keywords()
            close = _closest(kws | aliases, w)
            issues.append(Issue(
                severity="error", code="unknown_keyword",
                message=f"`{w}` is not an ORCA keyword.",
                suggestion=(f"Did you mean `{close}`?" if close
                            else "Check the ORCA manual chapter 6."),
            ))
    return issues


#: Blocks ORCA knows but the OPI catalog does not list. Both are real
#: ORCA blocks DELFIN itself writes; %pal is range-checked in
#: delfin/api.py (validate_orca_input %pal branch) and %maxcore in the
#: %maxcore branch right below it, so the vocabulary check must not
#: flag them.
_CATALOG_GAP_BLOCKS = {"maxcore", "pal"}


def check_blocks(text: str) -> list[Issue]:
    """Check ``%block ... end`` names and options against the catalog."""
    from delfin.orca_knowledge import known_block_option
    issues: list[Issue] = []
    for m in re.finditer(r"^\s*%\s*([A-Za-z0-9_]+)(.*)$", text, re.MULTILINE):
        block = m.group(1)
        if block.lower() in _CATALOG_GAP_BLOCKS:
            continue
        if not known_block_option(block):
            issues.append(Issue(
                severity="error", code="unknown_block",
                message=f"`%{block}` is not an ORCA block.",
                suggestion="Check the ORCA manual chapter 7 for block names.",
            ))
            continue
        # options until the block's end: next standalone end line
        tail = text[m.end():]
        stop = re.search(r"^\s*end\s*$", tail, re.MULTILINE | re.IGNORECASE)
        body = tail[:stop.start()] if stop else tail
        for opt_m in re.finditer(r"^\s*([A-Za-z0-9_]+)\b", body, re.MULTILINE):
            option = opt_m.group(1)
            if option.lower() == "end":
                break
            if not known_block_option(block, option):
                # ORCA warns and continues on unknown block options; the
                # parser probe (phase 3) decides inputs this cannot.
                issues.append(Issue(
                    severity="warning", code="unknown_block_option",
                    message=f"`{option}` is not a known option of `%{block}`.",
                    suggestion="Verify against the ORCA manual or use "
                               "the parser probe for this input.",
                ))
    return issues


def check_all(text: str) -> list[Issue]:
    """All static simple-input checks of this module."""
    return check_keyword_line(text) + check_blocks(text)
