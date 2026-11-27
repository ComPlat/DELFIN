"""Characterization: the tool catalogue changes only on purpose.

Pins the catalogue's shape and content so that an unintended change is
caught: tool names in order, a sha256 of the sorted-key JSON of the whole
catalogue, and the schema-token measurements. The object-identity check
(api_client._DOC_TOOLS_OPENAI is tool_catalog._DOC_TOOLS_OPENAI) only
passes once the api_client patch is applied; until then it is reported
skipped so this file runs green on base.

Input: the built catalogue. Output: whether it equals the recorded
baseline. Semantics: a drift detector, not a budget. The authority on the
schema's token TOTAL is tests/test_tool_schema_budget.py, which ratchets
it; this file only answers whether anything moved. Where both pin the same
number they are two copies of one answer, and two copies drift -- so the
constants below are a recorded baseline, not a second ceiling.

Moving a constant here is permitted and requires the deliberate cause to
be named in the commit that moves it. The baseline was first taken at
origin/main 026eba8f08cce9a86daf7bf6139c54906428e6c3, when the question
was whether a refactor had changed the catalogue. That refactor landed.
The values were re-taken on 2026-11-20 because the sub-agent tool's
description is now derived from the preset registry instead of written by
hand, so its text changes whenever a preset is added: the token total
moved 11,678 -> 11,718 (+40, measured) and the token estimate 11,766 ->
11,806, which the budget test accounts for. The tool COUNT and the name
order did not move, and those two pins still carry their original
baseline.

Re-taken again on 2026-11-26 (wave 58, package n1). THE DELIBERATE CAUSE:
the ``subagent`` tool's description gained the sentence "resume_id resumes
a finished run.", and the ``resume_id`` property text was shortened to pay
for it ("sa_id of a FINISHED subagent to continue with its context
replayed." -> "Finished run's sa_id."). Measured, not estimated: the
catalogue total went 11,865 -> 11,862 and the estimate 11,954 -> 11,951,
i.e. it is 3 tokens CHEAPER, and the count and the name order are unmoved.
The reason the sentence belongs in the description at all: ``subagent`` is
in ``_HOT_TOOL_NAMES``, so the description is the field every request
carries, while a property description is only read by a model that already
decided to look for it -- the arm measured that a capability living only
in the property is a capability the model does not use.
"""
import hashlib
import json

import pytest

from delfin.agent import api_client

# --- the recorded baseline; see the module docstring before moving one ---
# count and name order: unmoved since origin/main 026eba8f0 until
#   2026-11-23, when ONE tool was added: describe_native_tool; then
#   2026-11-24, when ONE more was added: lab_run (the self-driving lab's
#   first agent-facing entrance -- a whole workflow judged, then run by the
#   one coordinator).
# content and token figures: re-taken 2026-11-24 for that same cause.
#
# THE DELIBERATE CAUSE, as this file requires it to be named: the lab
# specialist. Measured on the tree that ships it, in two parts: +9 tokens
# for the preset (its name in the subagent_type enum and one more entry in
# the derived writer list -- the same cost ANY added preset pays) and +198
# for lab_run's own schema. This file's figures count the CATALOGUE; the
# number that ships is ratcheted in
# tests/test_tool_schema_budget.py::_SHIPPED_BUDGET, which did not move.
#
# RE-TAKEN AGAIN 2026-11-26 on the MERGED tree, one cause: the
# `chemotion-eln` preset. The preset's NAME is what the derived `subagent`
# description grows by, and it measures +8 tokens here -- estimate
# 12,162 -> 12,170, report total 12,072 -> 12,080, and both digests with
# them. The sibling ratchet tests/test_tool_schema_budget.py had already
# been re-measured on the merged tree (_TOKEN_BUDGET = 12_080); this file
# still carried the figure the ELN branch measured alone, so the two pins
# asserted two different totals for one catalogue. A merge resolves a
# numeric ceiling by measuring the result, not by keeping one branch's
# number. COUNT (98) and name order did not move, and their pins are
# untouched, because moving a name is a different act from adding a preset.
#
# RE-TAKEN 2026-11-27 on the merged tree, one cause: the subagent tool's
# resume_id sentence was reworded (wave 58). Measured: estimate
# 12,170 -> 12,167, report total 12,080 -> 12,077; COUNT and names unmoved.
PINNED_COUNT = 98
PINNED_NAMES_SHA256 = (
    "7903960c5ec64ec468cbd81d93b4d3b4be10d1823ab2a4d7a18547fd5dfc449f")
PINNED_BLOB_SHA256 = (
    "c1437b7a83210e47b4a639e2d776dcd92942735d5fd32ed96cd1a3ab9340c519")
PINNED_ESTIMATE = 12167
PINNED_REPORT_SHA256 = (
    "ade091e769307f1562e5a0315d1006e11631e4f7f2691fe4386b6be6eebd0356")
PINNED_REPORT_TOTAL_TOKENS = 12077


def _names_sha256() -> str:
    names = [t["function"]["name"] for t in api_client._DOC_TOOLS_OPENAI]
    return hashlib.sha256("\n".join(names).encode()).hexdigest()


def _blob_sha256() -> str:
    blob = json.dumps(api_client._DOC_TOOLS_OPENAI, sort_keys=True)
    return hashlib.sha256(blob.encode()).hexdigest()


def _report_sha256() -> str:
    rep = api_client.tool_schema_token_report()
    return hashlib.sha256(json.dumps(rep, sort_keys=True).encode()).hexdigest()


def test_catalogue_count_is_unchanged():
    assert len(api_client._DOC_TOOLS_OPENAI) == PINNED_COUNT


def test_catalogue_names_in_order_are_unchanged():
    names = [t["function"]["name"] for t in api_client._DOC_TOOLS_OPENAI]
    assert len(names) == PINNED_COUNT
    assert _names_sha256() == PINNED_NAMES_SHA256


def test_catalogue_json_sha256_is_unchanged():
    assert _blob_sha256() == PINNED_BLOB_SHA256


def test_estimate_schema_tokens_unchanged_for_full_list():
    assert (api_client.estimate_schema_tokens(api_client._DOC_TOOLS_OPENAI)
            == PINNED_ESTIMATE)


def test_tool_schema_token_report_unchanged_for_full_list():
    assert _report_sha256() == PINNED_REPORT_SHA256
    rep = api_client.tool_schema_token_report()
    assert rep["total_tokens"] == PINNED_REPORT_TOTAL_TOKENS
    assert rep["count"] == PINNED_COUNT


def test_moved_names_live_in_tool_catalog_and_are_the_same_objects():
    # importorskip: tool_catalog exists only from the second commit on, so
    # this test file still runs green on base.
    pytest.importorskip("delfin.agent.tool_catalog")
    from delfin.agent import tool_catalog

    assert api_client._DOC_TOOLS_OPENAI is tool_catalog._DOC_TOOLS_OPENAI
    assert (api_client._SCHEMA_CHARS_PER_TOKEN
            == tool_catalog._SCHEMA_CHARS_PER_TOKEN)
    assert (api_client.estimate_schema_tokens
            is tool_catalog.estimate_schema_tokens)
    assert (api_client.tool_schema_token_report
            is tool_catalog.tool_schema_token_report)
    assert (api_client._subagent_caps_phrase
            is tool_catalog._subagent_caps_phrase)


def test_tool_catalog_imports_nothing_from_api_client():
    # importorskip keeps the file runnable on base, where the module is
    # absent; after the move this asserts the direction of the dependency.
    pytest.importorskip("delfin.agent.tool_catalog")
    import ast
    import pathlib

    import delfin.agent.tool_catalog as tc

    tree = ast.parse(pathlib.Path(tc.__file__).read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            assert all(not a.name.startswith("delfin")
                       for a in node.names), node.names
        if isinstance(node, ast.ImportFrom):
            if node.level:
                mod = "." * node.level + (node.module or "")
            else:
                mod = node.module or ""
            assert not mod.startswith("delfin.agent.api_client"), mod
    # No bare `assert ... or True` here: this file's own module-level import of
    # api_client (line 34, for the catalogue constants) puts it in
    # `sys.modules` before this test runs, so the condition was always false
    # and the line could never fail. The loop above is the real check -- it
    # asserts the dependency DIRECTION from the source, per import node. A
    # fresh-interpreter check would assert something this test cannot: that
    # importing the catalogue pulls nothing in.
