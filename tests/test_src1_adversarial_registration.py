"""Adversarial: SRC1 P2 registration patch (reviewer nacht-s4).

Two of these tests are RED on the unpatched tree BY DESIGN
(test_executor_never_raises_through, test_unknown_action_is_data_not_crash):
they drive the real _DocToolExecutor.execute with "verify_source" and
pin the executor behaviour the patch adds; the control run shows
'Unknown tool: verify_source.' until the operator builds the patch.
Two are green already (the patch-text audit and the policy-escape
finding pinned below).

The policy case (test_network_policy_does_not_bind_verify_source) pins
the FINDING: verify_source performs outbound HTTP but is invisible to
web_fetch-shaped network rules, so a user's or managed network ceiling
does not stop it. It is green as written (the escape exists); if the
transport is later fixed to carry the host through the rule layer, this
test turns red and must be rewritten along its intent.
"""
from __future__ import annotations

import json

import pytest


def test_executor_wraps_result_untrusted():
    """(3) result always _wrap_untrusted — reads the patch directly.

    The success path wraps; the error paths deliberately return plain
    json.dumps (they carry no API content — an error message is our own
    text, not untrusted data; the same pattern _execute_web_fetch uses).
    So the assertion is: the payload return is wrapped AND every other
    return of _execute_verify_source is either an error json.dumps or
    the wrapped one — nothing escapes raw.
    """
    import pathlib
    patch = pathlib.Path(
        "/pfs/data6/home/ka/ka_ibcs/ka_ew7404/software/delfin/"
        ".delfin/worktrees/runde3-s3/.gate/"
        "src1_verify_source_registration.patch")
    text = patch.read_text()
    body = text.split("def _execute_verify_source")[1].split(
        "def _execute_session_search")[0]
    returns = [ln.strip().lstrip("+").strip() for ln in body.splitlines()
               if ln.strip().lstrip("+").strip().startswith("return ")]
    wrapped = [r for r in returns if "_wrap_untrusted" in r]
    errors = [r for r in returns if r.startswith(
        "return json.dumps({\"error\"")]
    assert len(wrapped) == 1, wrapped
    assert len(returns) == len(wrapped) + len(errors), returns


def test_executor_never_raises_through():
    """(2b) forced failure inside verify_doi -> structured error string,
    driven through the REAL executor entry point.

    RED on the unpatched tree: 'Unknown tool' — the tool is not wired
    yet. Green at the operator build, where the executor routes the
    call into _execute_verify_source.
    """
    from unittest.mock import patch as mockpatch
    from delfin.agent.api_client import _DocToolExecutor

    with mockpatch("delfin.agent.source_apis.verify_doi",
                   side_effect=RuntimeError("HTTP 503 broke")):
        result = _DocToolExecutor().execute(
            "verify_source", {"action": "doi", "value": "10.1039/x"})
    assert isinstance(result, str)
    assert "verify_source failed" in result
    assert "HTTP 503" in result


def test_unknown_action_is_data_not_crash():
    """RED on the unpatched tree (unknown tool); green at build."""
    from delfin.agent.api_client import _DocToolExecutor

    result = _DocToolExecutor().execute(
        "verify_source", {"action": "nonsense", "value": "x"})
    assert "unknown action" in result


def test_network_policy_does_not_bind_verify_source():
    """THE FINDING: network rules are web_fetch-shaped; verify_source
    escapes them. Red until the transport carries the host through the
    rule layer (or the executor asks permission like web_fetch)."""
    from delfin.agent import permission_rules as pr

    subj = pr._subjects(  # the same builder web_fetch uses
        "web_fetch", {"url": "https://api.crossref.org/works/10.1/x"},
        "/tmp/ws")
    rule = pr.parse_rule("web_fetch(domain:api.crossref.org)")
    assert pr._hits(rule, "web_fetch", subj, "/tmp/ws"), (
        "sanity: the rule does bind web_fetch")

    # now the new tool, same host, through the same rule layer
    subj2 = pr._subjects(
        "verify_source", {"action": "doi", "value": "10.1039/x"},
        "/tmp/ws")
    assert not pr._hits(rule, "verify_source", subj2, "/tmp/ws"), (
        "unexpected: the rule binds verify_source (fixed?)")
    # and the managed ceiling does not see it either: _ceiling_rule's
    # gate is the NAME check (permission_rules.py:858) — it fires only
    # when bare_tool_name(name) == "web_fetch". A real RuleSet whose
    # policy is the managed Policy dataclass, empty, yields None for
    # the new tool through the same branch.
    from delfin.agent import managed_policy
    rs = pr.RuleSet(deny=(), ask=(), allow=(),
                    policy=managed_policy.Policy(sources=()))
    assert pr._ceiling_rule(rs, "verify_source",
                            {"url": "https://pubchem.ncbi.nlm.nih.gov/x"}
                            ) is None, (
        "the ceiling must not bind verify_source (fixed if this fails)")
