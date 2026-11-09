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
    # The registration patch lives only in the reviewer's tree (never
    # committed); resolve it relative to the repository root and skip
    # when it is absent — CI must not depend on a machine path.
    patch = pathlib.Path(".gate", "src1_verify_source_registration.patch")
    if not patch.is_file():
        pytest.skip("registration patch not present in this tree")
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

    Needs a REAL KitToolPermissions with a tmp workspace: the executor
    refuses every catalogue tool without a permissions object (the
    refusal without permissions is pinned separately below). Red on the
    unpatched tree ('Unknown tool' even WITH permissions); green at the
    operator build, where the executor routes the call into
    _execute_verify_source.
    """
    import tempfile

    from unittest.mock import patch as mockpatch
    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    with tempfile.TemporaryDirectory() as ws:
        perms = KitToolPermissions(mode="default", workspace=ws)
        with mockpatch("delfin.agent.source_apis.verify_doi",
                       side_effect=RuntimeError("HTTP 503 broke")):
            result = _DocToolExecutor().execute(
                "verify_source",
                {"action": "doi", "value": "10.1039/x"}, perms)
    assert isinstance(result, str)
    assert "verify_source failed" in result
    assert "HTTP 503" in result


def test_unknown_action_is_data_not_crash():
    """Unknown action -> structured 'unknown action' string, not a raise.

    With a REAL permissions object (tmp workspace) — the executor
    refuses without one before the tool body runs. Red on the unpatched
    tree ('Unknown tool'); green at the operator build.
    """
    import tempfile

    from delfin.agent.api_client import KitToolPermissions, _DocToolExecutor

    with tempfile.TemporaryDirectory() as ws:
        perms = KitToolPermissions(mode="default", workspace=ws)
        result = _DocToolExecutor().execute(
            "verify_source", {"action": "nonsense", "value": "x"}, perms)
    assert "unknown action" in result


def test_executor_refuses_without_permissions():
    """Pins the executor's fail-closed contract the build must keep: a
    catalogue tool called WITHOUT a permissions object is refused with
    the sandbox/permissions refusal — never executed. Green on every
    tree (the contract predates the hook; that is why the two tests
    above need a permissions object to reach the tool body at all).
    """
    from delfin.agent.api_client import _DocToolExecutor

    out = _DocToolExecutor().execute("web_fetch",
                                     {"url": "https://example.org/x"})
    data = json.loads(out)
    assert "error" in data, out
    assert "needs a workspace sandbox" in data["error"], out


def test_network_policy_binds_verify_source():
    """The operator's decision (2026-11-01): the policy BINDS the new
    tool. verify_source doi->Crossref / compound->PubChem are outbound
    HTTP like web_fetch, so a web_fetch(domain:) deny/ask rule and the
    managed allowed_domains_max ceiling must stop it the same way.

    RED on the pre-fix tree (the escape this package found); green once
    s3's permission_rules.py fix teaches the rule layer the
    verify_source host. If it goes red again after a refactor, the
    escape is back.
    """
    from delfin.agent import managed_policy, permission_rules as pr

    # sanity: web_fetch itself is bound (unchanged behaviour)
    subj = pr._subjects(
        "web_fetch", {"url": "https://api.crossref.org/works/10.1/x"},
        "/tmp/ws")
    rule = pr.parse_rule("web_fetch(domain:api.crossref.org)")
    assert pr._hits(rule, "web_fetch", subj, "/tmp/ws"), (
        "sanity: the rule binds web_fetch")

    # the new tool, same host, must meet the same rule
    subj2 = pr._subjects(
        "verify_source", {"action": "doi", "value": "10.1039/x"},
        "/tmp/ws")
    assert pr._hits(rule, "verify_source", subj2, "/tmp/ws"), (
        "verify_source escaped the web_fetch domain rule")

    # and the managed ceiling must see it. The host comes from the
    # action (host_for_action), never from a user-supplied url — the
    # decided design; a url argument contributes nothing.
    pol = managed_policy.Policy(sources=(),
                                allowed_domains_max=(("api.crossref.org",),))
    rs = pr.RuleSet(deny=(), ask=(), allow=(), policy=pol)
    assert pr._ceiling_rule(rs, "verify_source",
                            {"action": "compound", "value": "aspirin"}
                            ) is not None, (
        "the managed ceiling must bind verify_source outside it")
    assert pr._ceiling_rule(rs, "verify_source",
                            {"action": "doi", "value": "10.1/x"}
                            ) is None, (
        "a host inside the ceiling stays allowed")
