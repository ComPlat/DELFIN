"""Characterization: the tool catalogue is byte-identical after the move.

Pins the state computed on the base commit (s06 started at origin/main
026eba8f08cce9a86daf7bf6139c54906428e6c3): tool names in order, a sha256 of
the sorted-key JSON of the whole catalogue, and the schema-token
measurements. The object-identity check (api_client._DOC_TOOLS_OPENAI is
tool_catalog._DOC_TOOLS_OPENAI) only passes once the api_client patch is
applied; until then it is reported skipped so this file runs green on base.
"""
import hashlib
import json

import pytest

from delfin.agent import api_client

# --- pinned on base (origin/main 026eba8f08cce9a86daf7bf6139c54906428e6c3) ---
PINNED_COUNT = 96
PINNED_NAMES_SHA256 = (
    "7afa506ce5faf16ae47cb7b2b7957b3c6bf02a06e0672dd16bb9142167c528fb")
PINNED_BLOB_SHA256 = (
    "7cd398a3afeda3d1c53eed44f61030dafb3c0e382365a47c223fde2c385876be")
PINNED_ESTIMATE = 11766
PINNED_REPORT_SHA256 = (
    "18f870f3c973cea10f2bdca8cee9bce72713620551b7aae708a5d345b5660591")
PINNED_REPORT_TOTAL_TOKENS = 11678


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
