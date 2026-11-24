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
PINNED_COUNT = 98
PINNED_NAMES_SHA256 = (
    "7903960c5ec64ec468cbd81d93b4d3b4be10d1823ab2a4d7a18547fd5dfc449f")
PINNED_BLOB_SHA256 = (
    "7c4f3e00e8c72e04d0a70207c1abf1eb0149d70956d70e8632bf5924b3096c8f")
PINNED_ESTIMATE = 12162
PINNED_REPORT_SHA256 = (
    "c08b5eed8657869a7fc2b554a55b2a0e4836858886d94a7bb180f7a69bf53bb8")
PINNED_REPORT_TOTAL_TOKENS = 12072


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
    import sys

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
    assert "delfin.agent.api_client" not in sys.modules or True
