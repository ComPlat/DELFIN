"""tests/fixtures/unmounted_capabilities.txt must not advertise a closed
class as an open one.

Input: the fixture file the capability suite reads, and the tree itself
(``delfin/``), read with the same word-boundary rule
tests/test_every_capability_has_a_caller.py uses.

Why this test exists. Two rows of that fixture were measured on
2026-11-22 and have since been repaired or deleted, but the prose beside
them kept describing the defects in the present tense. A wave-45 package
read that prose as a standing order ("three redirect forms currently RUN
unasked; the repair is the two classes above, in api_client.py; the
self-modification guard refuses it to an unattended session") and was
about to add a second judgement at ``_run_permission_gate`` -- which would
have REPLACED the write-path refusal that already covers every one of
those forms at the tool entry with a weaker question.

Two things are pinned here:

* a name whose module has been deleted is not described as an open wiring
  contract, because the next reader will look for the contract file and
  find it gone -- as this package did;
* the redirect class is recorded as CLOSED at the tool entry, with the
  measurement that closed it, rather than as running unasked.

Both are asserted over the file's own text, and the tree check is what
makes the first one falsifiable: if ``delfin.agent.shared_tmp`` comes
back, the name is no longer deleted and the marker stops being required.
"""
from __future__ import annotations

import ast
import sys
from pathlib import Path

BASELINE = (Path(__file__).resolve().parent / "fixtures"
            / "unmounted_capabilities.txt")
REPO = Path(__file__).resolve().parent.parent

#: Capabilities deleted with delfin/agent/shared_tmp.py in a92ed6e36
#: ("shared_tmp: delete the unwired module and its tests", 2026-11-22).
DELETED_WITH_SHARED_TMP = (
    "bwrap_substitution",
    "map_to_host",
    "enforce_quota",
    "sweep_stale",
)


def _text() -> str:
    return BASELINE.read_text(encoding="utf-8")


#: This file names the identifiers it looks for, so it must not count as a
#: mention of them -- otherwise the test fails on its own source instead of
#: on the fixture, which is a red for the wrong reason.
_SELF = Path(__file__).resolve()


def _mentions_outside_the_fixture(name: str) -> list[str]:
    """Files under delfin/ and tests/ (not the fixture, not this test) that
    name *name* in CODE -- an import, an attribute or a bare identifier.

    Comments and docstrings are not counted: they are prose about history,
    not a surface the product reaches, and the test below pins exactly that
    ("named by running code again"). Measured 2026-11-25 while integrating
    two arms: a ledger test that names these four symbols in its docstring,
    to explain why they are gone, tripped a scan that read raw text -- so
    the honest test of one arm failed the honest fixture of the other.
    """
    hits: list[str] = []
    for root in (REPO / "delfin", REPO / "tests"):
        for path in root.rglob("*.py"):
            if path.resolve() == _SELF:
                continue
            try:
                tree = ast.parse(path.read_text(encoding="utf-8",
                                                errors="replace"))
            except (OSError, SyntaxError):
                continue
            for node in ast.walk(tree):
                named = (
                    isinstance(node, ast.Name) and node.id == name
                    or isinstance(node, ast.Attribute) and node.attr == name
                    or isinstance(node, (ast.Import, ast.ImportFrom))
                    and any((a.name == name or a.name.rsplit(".", 1)[-1] == name)
                            for a in node.names)
                )
                if named:
                    hits.append(str(path.relative_to(REPO)))
                    break
    return hits


def _block_for(name: str) -> str:
    """The contiguous run of comment lines that opens with ``# <name>``."""
    lines = _text().splitlines()
    for i, line in enumerate(lines):
        if line.startswith(f"# {name}"):
            block = [line]
            for nxt in lines[i + 1:]:
                if not nxt.startswith("#"):
                    break
                block.append(nxt)
            return "\n".join(block)
    return ""


def test_the_deleted_shared_tmp_names_are_marked_as_history():
    """A capability that no longer exists is not described as awaiting a
    wiring commit -- this package looked for
    tests/test_shared_tmp_wiring.py and found nothing."""
    assert not (REPO / "delfin" / "agent" / "shared_tmp.py").exists(), (
        "delfin/agent/shared_tmp.py exists again: the names below are no "
        "longer deleted, and this test's premise has to be re-measured")

    unmarked = []
    for name in DELETED_WITH_SHARED_TMP:
        assert _mentions_outside_the_fixture(name) == [], (
            f"{name} is named by running code again; it is not a deleted "
            "name any more")
        block = _block_for(name)
        # No prose is not a misleading claim: only a block that DOES
        # describe the name has to say it is gone. sweep_stale had a data
        # row and no prose before a92ed6e36 removed it.
        if block and "DELETED" not in block:
            unmarked.append(name)
    assert unmarked == [], (
        "these names document a module deleted in a92ed6e36 as if it were "
        f"still there to wire, with no 'DELETED' marker: {unmarked}")


def test_the_redirect_class_is_recorded_as_closed_at_the_tool_entry():
    """The gate's redirect class was measured at the PREDICATE on
    2026-11-22 and at the TOOL ENTRY the same day, where the write-path
    gate refuses all of it. The fixture must say so; the present-tense
    "RUN, no question" reading is what a later package acted on."""
    block = _block_for("is_safe")
    assert block, f"no is_safe block in {BASELINE.name}"

    assert "would write to" in block, (
        "the is_safe block must name the wording that closes the redirect "
        "class at the tool entry, so a reader can re-measure it")
    assert "tool entry" in block, (
        "the block must say WHERE the class is closed; naming the gate "
        "alone is what let it read as an open defect")
    for stale in (
        "RUN, no question; _changes_outside_the_workspace does not read",
    ):
        assert stale not in block, (
            f"the block still claims {stale!r}, which was measured FALSE at "
            "the tool entry (tests/test_a_flag_that_runs_a_program_is_not_a_"
            "reading_command.py::test_the_write_path_gate_already_refuses_"
            "an_outside_target)")


def test_the_mention_scan_can_still_fail_on_code(tmp_path, monkeypatch):
    """Control on the narrowed scan: a real CODE reference is still caught.

    Without this, `_mentions_outside_the_fixture` could be relaxed to
    "always empty" and every test above would stay green. Failing input:
    a module that imports one of the deleted names -- which is exactly the
    regression the scan exists to catch.
    """
    fake = tmp_path / "delfin"
    fake.mkdir()
    (fake / "revived.py").write_text(
        "from delfin.agent.shared_tmp import map_to_host\n", encoding="utf-8")
    (fake / "prose.py").write_text(
        '"""map_to_host is mentioned here only in prose."""\n', encoding="utf-8")
    monkeypatch.setattr(
        sys.modules[__name__], "REPO", tmp_path)
    assert _mentions_outside_the_fixture("map_to_host") == ["delfin/revived.py"], (
        "the scan must catch an import of a deleted name and must NOT count "
        "a docstring that merely mentions it")
