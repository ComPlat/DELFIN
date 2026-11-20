"""SRC1: the network policy binds verify_source (operator 2026-11-01).

The green half runs against the LIVE ``delfin/agent/permission_rules.py``:
the SRC1 fix is merged, so the live module is the "new" copy — what the
original test loaded from ``.gate/px/new_permission_rules.py`` (an
untracked, operator-built copy that exists only transiently in a
protected-file session and is absent in every clean checkout, which is
why those tests skipped everywhere).

The red half is the escape control. Its pre-fix subject is no longer in
any working tree, so it is loaded from the pinned git-history blob of
``permission_rules.py`` at the commit before the fix — the exact code
the escape lived in — never a monkeypatched reconstruction. It is a
regression lock on the history (full-depth clones, which CI is; a local
shallow checkout names the gap in its own skip).
"""

from __future__ import annotations

import importlib.util
import os
import subprocess
import sys
import tempfile

import pytest

from delfin.agent import managed_policy as mp
from delfin.agent import permission_rules as pr

# The commit before the SRC1 fix (8a16670f3): its permission_rules.py is
# the one verify_source escaped.  Pin it — re-pin if history rewrites it.
PRE_FIX = "541471d175883b4536ea9e4b9a8aba356e11fd32"
_LIVE_RULE = "web_fetch(domain:api.crossref.org)"


def _load_pre_fix():
    """Load the pre-fix permission_rules as a separate module.

    Written to a temp file and imported as a member of the real
    delfin.agent package, so the blob's relative imports (``.source_apis``)
    resolve, and the live permission_rules stays registered under its own
    name.
    """
    try:
        src = subprocess.check_output(
            ["git", "show", f"{PRE_FIX}:delfin/agent/permission_rules.py"],
        ).decode()
    except Exception as exc:  # shallow/missing history, or outside the repo
        pytest.skip(
            "the pre-fix commit is not available here "
            f"({exc.__class__.__name__}) — run inside a full-history clone")
    fd, path = tempfile.mkstemp(suffix="_permission_rules.py")
    with os.fdopen(fd, "w") as f:
        f.write(src)
    try:
        name = "delfin.agent._policy_under_test_pre_fix"
        spec = importlib.util.spec_from_file_location(name, path)
        mod = importlib.util.module_from_spec(spec)
        sys.modules[name] = mod
        spec.loader.exec_module(mod)
        return mod
    finally:
        os.unlink(path)


def test_verify_source_is_bound_by_a_web_fetch_domain_rule():
    """The live rule layer binds verify_source to a web_fetch domain rule."""
    rule = pr.parse_rule(_LIVE_RULE)
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


def test_the_pre_fix_code_shows_the_escape():
    """The red half of the pair: the pre-fix code let verify_source escape.

    Loaded from git history, not a reconstruction — it proves the escape
    was real at the pinned commit and that the live module (tested
    alongside) is what closes it.
    """
    pf = _load_pre_fix()
    rule = pf.parse_rule(_LIVE_RULE)
    # the pre-fix code: web_fetch is bound...
    subj = pf._subjects("web_fetch",
                        {"url": "https://api.crossref.org/works/10.1/x"},
                        "/tmp/ws")
    assert pf._hits(rule, "web_fetch", subj, "/tmp/ws")
    # ...but verify_source produced no host at all, so the rule missed it
    subj2 = pf._subjects("verify_source",
                         {"action": "doi", "value": "10.1039/x"},
                         "/tmp/ws")
    assert subj2["hosts"] == [], (
        f"pre-fix code should derive no host, got {subj2['hosts']}")
    assert not pf._hits(rule, "verify_source", subj2, "/tmp/ws"), (
        "the escape is gone from the pre-fix commit — re-pin PRE_FIX "
        "to the commit before the fix")
    # and the live module is the one that closes it: same args, it
    # derives the host, so the same rule now hits
    subj_live = pr._subjects("verify_source",
                             {"action": "doi", "value": "10.1039/x"},
                             "/tmp/ws")
    assert subj_live["hosts"] != [], (
        f"live module derives no host for a known action: {subj_live['hosts']}")
    assert pr._hits(pr.parse_rule(_LIVE_RULE), "verify_source", subj_live,
                    "/tmp/ws")


def test_verify_source_unknown_action_gets_no_host():
    subj = pr._subjects("verify_source",
                        {"action": "scrape", "value": "x"}, "/tmp/ws")
    assert subj["hosts"] == [], "an unknown action must get NO host"


def test_managed_ceiling_stops_verify_source_outside_it():
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
    pol = mp.Policy(sources=(),
                    allowed_domains_max=(("api.crossref.org",),))
    rs = pr.RuleSet(deny=(), ask=(), allow=(), policy=pol)
    assert pr._ceiling_rule(
        rs, "web_fetch", {"url": "https://example.org/x"}) is not None
