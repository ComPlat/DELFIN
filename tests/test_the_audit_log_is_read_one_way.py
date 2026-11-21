"""The audit log is read one way, and the deleted spellings stay deleted.

Provenance ("cite-sources"): every line in this file cites the code or the
command it was measured against.

``read_last_n`` (audit_log.py:621 before this change) was a second copy of
"return the last records of one log", reached by 7 test files and by
nothing in ``delfin/``. ``_diff_line_counts`` (audit_log.py:642) was
reached by nothing at all. Both were deleted and the test read-backs moved
to the module's one record parser ``_parse_log`` — the parse step behind
every product read (``build_changes_report`` -> ``_read_all_records`` ->
``_parse_log``, all in delfin/agent/audit_log.py).

Input: the source tree (``delfin/`` and ``tests/``) and a throwaway audit
log written through the product's own ``append(make_record(...))``.
Output: three refuted failure modes for the read path:
  1. nothing re-spells the deleted helpers anywhere in source or text;
  2. the module really has one parser and not the deleted two;
  3. parsing stays parsing and verification stays ``verify``'s duty.

Units and limits, so nothing here is read as more than it measures: the
units are records and file-order positions (no time series); the tree
scan reads text, so a byte-compiled reference is not counted (it is not
a call site either); the parser requires dict records where the deleted
``read_last_n`` accepted any JSON value -- the one behavioural
difference, spelled out in scratch/PR.md.
"""

from __future__ import annotations

import json
from pathlib import Path

from delfin.agent import audit_log

ROOT = Path(__file__).resolve().parents[1]

#: The deleted spellings. A tree-wide scan is the guard because an import
#: grep does not catch a call site spelled through an alias
#: (``as _al`` / ``as _audit`` / ``as _al_undo`` all appear in delfin/).
DEAD_NAMES = ("read_last_n", "_diff_line_counts")

#: Not source text: reading these would spray identifiers into the scan.
_BINARY = frozenset({
    ".png", ".jpg", ".jpeg", ".gif", ".pdf", ".ico", ".woff", ".woff2",
    ".zip", ".gz", ".so", ".pyc", ".pyo", ".xlsx", ".docx", ".pptx",
    ".odt", ".ods",
})

#: Text artefacts allowed to name what they measure: the capability
#: ledger keeps its own rows (its rule is that a name leaves the list by
#: being connected or deleted, so the two rows go with this change; the
#: upstream capability scanner on origin/main enforces that with its
#: stale-row check). This file is excluded for the same reason that
#: scanner excludes its ledger file: a scan must not measure itself.
EXEMPT = frozenset({
    "tests/test_the_audit_log_is_read_one_way.py",
    "tests/fixtures/unmounted_capabilities.txt",
    "scratch/unmounted_capabilities.txt",
})


def _tree_text_hits() -> list[tuple[str, str, int]]:
    """Every (relative path, name, line) where a dead spelling is text.

    Scans ``delfin/`` and ``tests/`` like upstream's capability scan does:
    Python files as tokens would be finer, but words in a string literal
    ARE name-based dispatch (yaml/json/md mentions included here on
    purpose -- the wave question asks for them).
    """
    hits: list[tuple[str, str, int]] = []
    for base in (ROOT / "delfin", ROOT / "tests"):
        for path in sorted(base.rglob("*")):
            if not path.is_file() or path.suffix.lower() in _BINARY:
                continue
            rel = str(path.relative_to(ROOT))
            if rel in EXEMPT or "__pycache__" in path.parts:
                continue
            try:
                text = path.read_text(encoding="utf-8", errors="replace")
            except OSError:
                continue
            for lineno, line in enumerate(text.splitlines(), 1):
                for name in DEAD_NAMES:
                    if name in line:
                        hits.append((rel, name, lineno))
    return hits


def test_the_deleted_helpers_stay_deleted():
    """One parser, not three.

    Fails at the baseline (both helpers importable), passes after the
    deletion -- so this, not a passing-only pin, is one of the two tests
    that carry the change.
    """
    assert not hasattr(audit_log, "read_last_n"), (
        "read_last_n is the second copy of the read path; the tests read "
        "audit_log._parse_log instead")
    assert not hasattr(audit_log, "_diff_line_counts"), (
        "_diff_line_counts was reached by nothing; return it only with the "
        "caller that needs it")
    assert callable(audit_log._parse_log), "the one parser must exist"
    assert callable(audit_log.verify), "verification must stay reachable"


def test_a_second_reader_does_not_come_back():
    """Repo-wide scanner for the deleted spellings.

    Fails at the baseline (19 source mentions across 7 test files plus the
    definition lines) and must never fail again. An alias-spelled call
    site, a name in a yaml/json/md dispatch or a re-import in a new test
    would all land here.
    """
    hits = _tree_text_hits()
    assert hits == [], (
        f"deleted audit_log helpers mentioned again at: {hits[:10]}")


def test_reading_the_last_records_is_parsing_and_verification_stays_a_separate_duty(tmp_path):
    """Never-worse pin (passes in both arms; not proof of the change).

    Drives the product's real write entry (``append(make_record(...))``)
    and its one parser; then shows the duty split the brief asks about:
    the parser returns what it can and fingers nothing, ``verify`` names
    the defect through its real dict contract. Assertion provenance:
    verify()'s docstring in delfin/agent/audit_log.py (kinds
    ``unparseable``/``chain``/``mac``/"head" and problems shaped
    ``{"line", "kind", "detail"}`` with a 1-based ``line``).
    """
    log = tmp_path / "audit.log"
    for i in range(5):
        audit_log.append(
            audit_log.make_record(tool="write_file", decision="ok",
                                  path=f"f{i}.py"),
            log_path=log)

    records = audit_log._parse_log(log)
    assert [r["path"] for r in records[-3:]] == ["f2.py", "f3.py", "f4.py"], (
        "last 3 of 5, file order -- what read_last_n(3) returned")
    assert len(audit_log._parse_log(log)) == 5

    # A torn line in the middle: the parser's answer is definition-based
    # (dict records only), the verifier's is a named defect.
    text = log.read_text(encoding="utf-8")
    lines = text.splitlines()
    mangled = json.loads(lines[2])
    mangled["path"] = "edited-after-the-fact"
    lines[2] = json.dumps(mangled, sort_keys=True)
    log.write_text("\n".join(lines) + "\n", encoding="utf-8")

    still_five = audit_log._parse_log(log)
    assert len(still_five) == 5, "the edited line still parses as a dict"
    assert still_five[2]["path"] == "edited-after-the-fact"

    report = audit_log.verify(log)
    assert report["ok"] is False, (
        "an edited line is exactly what verify is for")
    kinds = {p["kind"] for p in report["problems"]}
    assert "chain" in kinds or "mac" in kinds or "head" in kinds, report
    assert all(isinstance(p["line"], int) and p["line"] > 0
               for p in report["problems"]), report