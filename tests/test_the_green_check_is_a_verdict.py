"""The check mark said "checked", and meant "a retry was paid for".

The engine's claim-grounding guard keeps two facts apart: that the single
correction turn for this user turn has been SPENT, and that a re-scan of
the correction found the claims grounded. Only the second is a verdict,
and only the second may put "✓ Self-verification: statements checked,
answer corrected." in front of the reader.

The dashboard runs three forced corrections of its own -- for a citation
of a file that does not exist, for an ORCA keyword the manual does not
know, and for a physical quantity nothing in the turn measured. All three
printed the green check the moment the retry produced any text at all. A
correction that raised a NEW unsupported claim, or that only wrapped the
old one in hedges, was reported to the user as verified.

The line is now produced in exactly one place, from a boolean that a
re-scan of the corrected answer decides.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

from delfin.dashboard import tab_agent as T

_SRC = pathlib.Path(inspect.getfile(T)).read_text(encoding="utf-8")


class _Flag:
    """Stand-in for a verify_guard flag: the scanners return objects with
    a ``kind``, and the citation scan mixes hard and soft ones."""

    def __init__(self, kind: str = "nonexistent") -> None:
        self.kind = kind


# ---------------------------------------------------------------------------
# the verdict itself
# ---------------------------------------------------------------------------

def test_a_clean_rescan_earns_the_check():
    note = T._correction_verdict_note(lambda text: [], "the corrected answer")
    assert note.startswith("✓")


def test_a_correction_that_left_the_claims_ungrounded_says_so():
    note = T._correction_verdict_note(
        lambda text: [_Flag()], "the corrected answer")
    assert "✓" not in note
    assert "unconfirmed" in note


def test_a_correction_that_raised_a_new_claim_is_not_a_verification():
    """The re-scan reads the CORRECTION, so a fresh ungrounded claim in it
    denies the check exactly like a surviving old one."""
    seen: list[str] = []

    def _scan(text):
        seen.append(text)
        return [_Flag()] if "src/nope.py" in text else []

    note = T._correction_verdict_note(_scan, "now see src/nope.py:12")
    assert seen == ["now see src/nope.py:12"]
    assert "✓" not in note


def test_a_scanner_that_raises_is_not_a_pass():
    """Fail closed: the check has to be earned, and an answer that could
    not be scanned has not earned it."""
    def _boom(text):
        raise RuntimeError("scanner died")

    assert "✓" not in T._correction_verdict_note(_boom, "text")


def test_the_soft_flags_do_not_deny_the_check():
    """The citation scan also returns 'unread' flags. Those never forced
    the correction, so they must not veto its verdict either."""
    note = T._correction_verdict_note(
        lambda text: [_Flag("unread")], "text",
        keep=lambda f: f.kind == "nonexistent")
    assert note.startswith("✓")


def test_the_scanner_keeps_its_keyword_arguments():
    """The sites pass repo_root / observed_files / evidence_tools_used --
    a re-scan run without them judges against a different world than the
    scan that forced the correction."""
    got: dict = {}

    def _scan(text, **kwargs):
        got.update(kwargs)
        return []

    T._correction_verdict_note(_scan, "text", observed_files={"a.py"},
                               evidence_tools_used={"bash"})
    assert got == {"observed_files": {"a.py"},
                   "evidence_tools_used": {"bash"}}


# ---------------------------------------------------------------------------
# ...and every site reads it
# ---------------------------------------------------------------------------

def test_the_line_exists_in_exactly_one_place():
    """Three of the four sites were literals. A literal cannot be wrong
    about the verdict because it never asks."""
    assert _SRC.count("✓ Self-verification") == 1
    assert "✓ Self-verification" in inspect.getsource(
        T._self_verification_note)


def test_no_chat_message_states_the_verdict_as_a_constant():
    """Read off the calls, not off a line: no _append_system_message may
    hand the user a hard-coded self-verification claim."""
    tree = ast.parse(_SRC)
    bad: list[str] = []
    for node in ast.walk(tree):
        if not isinstance(node, ast.Call):
            continue
        name = getattr(node.func, "id", "") or getattr(node.func, "attr", "")
        if name != "_append_system_message":
            continue
        for arg in node.args:
            for sub in ast.walk(arg):
                if (isinstance(sub, ast.Constant)
                        and isinstance(sub.value, str)
                        and "Self-verification" in sub.value):
                    bad.append(f"line {node.lineno}")
    assert not bad, (
        "a self-verification verdict printed as a constant at " + ", ".join(bad))


def test_each_forced_correction_asks_for_a_verdict():
    """One call per dashboard-driven correction turn: citation, keyword,
    quantity. Fewer means a retry that reports itself."""
    tree = ast.parse(_SRC)
    calls = [n for n in ast.walk(tree)
             if isinstance(n, ast.Call)
             and getattr(n.func, "id", "") == "_correction_verdict_note"]
    assert len(calls) >= 3
    scanners = {getattr(c.args[0], "attr", "") for c in calls}
    assert scanners == {"scan_for_ungrounded_code_claims",
                        "scan_for_unverified_keywords",
                        "scan_for_unsourced_quantities"}
