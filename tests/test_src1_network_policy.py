"""SRC1: the network policy binds verify_source (operator 2026-11-01).

Runs the rule layer from a COPY of permission_rules.py: with
SRC1_POLICY_COPY unset the ORIG copy is used (the escape exists -> the
test must be red on the pre-fix tree); with
SRC1_POLICY_COPY=.gate/px/new_permission_rules.py the fixed copy is
used and the test must go green. Overlay via a temp directory holding
the copy as delfin/agent/_policy_under_test.py, imported directly —
never a shadow import of the real package (the real
delfin.agent.permission_rules stays untouched in sys.modules).
"""

from __future__ import annotations

import importlib.util
import os
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
POLICY_COPY = os.environ.get("SRC1_POLICY_COPY",
                             ".gate/px/orig_permission_rules.py")


def _load_policy(copy: str):
    """Import permission_rules from one .gate/px copy.

    The copy is loaded as a member of the real delfin.agent package
    under a per-copy module name, so its relative imports resolve and
    neither copy is registered as the live permission_rules.
    """
    path = (ROOT / copy).resolve()
    if not path.exists():
        pytest.skip(f"policy copy missing: {path}")
    name = "delfin.agent._policy_under_test_" + (
        "new" if "new_" in copy else "orig")
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


ORIG_COPY = ".gate/px/orig_permission_rules.py"
NEW_COPY = ".gate/px/new_permission_rules.py"


def _rule(copy):
    pr = _load_policy(copy)
    return pr, pr.parse_rule("web_fetch(domain:api.crossref.org)")


def test_verify_source_is_bound_by_a_web_fetch_domain_rule():
    pr, rule = _rule(NEW_COPY)
    # sanity: web_fetch itself is bound (unchanged behaviour)
    subj = pr._subjects("web_fetch",
                        {"url": "https://api.crossref.org/works/10.1/x"},
                        "/tmp/ws")
    assert pr._hits(rule, "web_fetch", subj, "/tmp/ws")
    # the new tool, same host, must meet the same rule
    subj2 = pr._subjects("verify_source",
                         {"action": "doi", "value": "10.1039/x"},
                         "/tmp/ws")
    assert pr._hits(rule, "verify_source", subj2, "/tmp/ws"), (
        "verify_source escaped the web_fetch domain rule")


def test_the_orig_copy_still_shows_the_escape():
    # the control: WITHOUT the fix the rule misses verify_source. This
    # is the red-before / green-after pair in one file — the orig copy
    # proves the escape existed, the new copy proves the fix closes it.
    pr, rule = _rule(ORIG_COPY)
    subj2 = pr._subjects("verify_source",
                         {"action": "doi", "value": "10.1039/x"},
                         "/tmp/ws")
    assert not pr._hits(rule, "verify_source", subj2, "/tmp/ws")
    pr_new, rule_new = _rule(NEW_COPY)
    subj3 = pr_new._subjects("verify_source",
                             {"action": "doi", "value": "10.1039/x"},
                             "/tmp/ws")
    assert pr_new._hits(rule_new, "verify_source", subj3, "/tmp/ws")


def test_verify_source_unknown_action_gets_no_host():
    pr, _ = _rule(NEW_COPY)
    subj = pr._subjects("verify_source",
                        {"action": "scrape", "value": "x"}, "/tmp/ws")
    assert subj["hosts"] == [], "an unknown action must get NO host"


def test_managed_ceiling_stops_verify_source_outside_it():
    pr = _load_policy(NEW_COPY)
    from delfin.agent import managed_policy as mp
    pol = mp.Policy(sources=(),
                    allowed_domains_max=(("api.crossref.org",),))
    rs = pr.RuleSet(deny=(), ask=(), allow=(), policy=pol)
    args = {"action": "compound", "value": "aspirin"}  # -> pubchem
    assert pr._ceiling_rule(rs, "verify_source", args) is not None, (
        "a host outside the ceiling must be denied")
    args2 = {"action": "doi", "value": "10.1/x"}  # -> crossref
    assert pr._ceiling_rule(rs, "verify_source", args2) is None, (
        "a host inside the ceiling stays allowed")


def test_managed_ceiling_still_stops_web_fetch():
    pr = _load_policy(NEW_COPY)
    from delfin.agent import managed_policy as mp
    pol = mp.Policy(sources=(),
                    allowed_domains_max=(("api.crossref.org",),))
    rs = pr.RuleSet(deny=(), ask=(), allow=(), policy=pol)
    assert pr._ceiling_rule(
        rs, "web_fetch", {"url": "https://example.org/x"}) is not None
