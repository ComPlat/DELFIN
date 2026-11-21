"""The tool-name judges agree, and a name one judge accepts is real.

Three places answer 'is this a real tool name': unknown_tool_names (the
CLI typo check), subagents._known_preset_tool_names (the preset
validator's catalogue) and the T4 grounding surface, all now reading
api_client.known_tool_names().  These tests are the pair lattice for that
agreement and for the sets that stay separate on purpose:

* the regression battery is what must FAIL at the baseline -- it pins the
  measured 2026-11-20 defect (the doc-list-only judgement called the ops
  server's geometry_build_smiles / geometry_read_xyz /
  geometry_add_hydrogens misspelt);
* the pins pass at both ends by construction and prove consistency, not
  the change:  api_client._WRITE_TOOL_NAMES stays inside the classifier's
  write set (the enforce-vs-classify split is deliberate; a hand copy
  drifted over make_plot), the lead surface stays exactly lead.py's two
  sets, and every documented catalogue name stays reachable.

Injected deliberately, no skipif: the ops introspection is the thing under
test and must fail loudly if the import cannot resolve -- the pre-change
code swallowed exactly that failure (bare except at
subagents.py:_known_preset_tool_names) and silently shrank its catalogue.
"""

from __future__ import annotations

import pytest

from delfin.agent import api_client as ac
from delfin.agent import subagents as SA
from delfin.agent.swarm import lead


# Names whose verdicts split the catalogues before the unification:
# ops tools (unknown to the advertised doc list and to the typo check)
# and lead tools (accepted by the typo check, refused by the preset
# validator's catalogue).
_OPS_NAMES = (
    # geometry_fix_hydrogens, not _add_: measured against the server
    # module, which defines 102 tool_* functions and no _add_hydrogens.
    # The branch asserted the misspelling must be ACCEPTED, which
    # would have weakened the typo judge for the case it exists for.
    "geometry_build_smiles", "geometry_read_xyz",
    "geometry_fix_hydrogens",
    "parse_orca_output", "find_orca_errors",
)
_LEAD_NAMES = tuple(sorted(ac._SWARM_LEAD_TOOL_NAMES))
_TYPOS = ("reed_file", "bahs")


def test_the_typo_check_reports_no_real_ops_tool():
    """The regression of 2026-11-20: the doc-tool list alone knew none of
    the server's 102 tool_* functions, so real tools were reported as
    misspelt.  Fails at the baseline."""
    assert ac.unknown_tool_names(list(_OPS_NAMES)) == []


def test_the_preset_catalogue_and_the_typo_check_agree_on_one_name():
    """A name one judge accepts must be known to the other -- this is the
    pair that diverged both ways (ops names one way, swarm_* the other).
    Fails at the baseline on swarm_status (and on _OPS_NAMES)."""
    for name in _OPS_NAMES + _LEAD_NAMES:
        assert ac.unknown_tool_names([name]) == [], name
        assert name in SA._known_preset_tool_names(), (
            f"{name}: accepted by the typo check but refused by the preset "
            "catalogue -- two answers to one question")


def test_a_misspelling_is_still_reported():
    """Tightening one judge must not loosen the other: a genuine typo
    still surfaces, an mcp__ name is still exempt (deliberate asymmetry:
    remote MCP servers appear only inside the turn).  Passes at both ends
    -- a pin of the behaviour that existed before."""
    # sorted(), so a set comparison; the point is WHICH names are
    # reported, not the order the fixture lists them in.
    assert set(ac.unknown_tool_names(["read_file"] + list(_TYPOS))) \
        == set(_TYPOS)
    assert ac.unknown_tool_names(["mcp__nowhere__x"]) == []
    for typo in _TYPOS:
        assert typo not in SA._known_preset_tool_names()


def test_every_enforced_write_still_counts_as_a_write_in_classification():
    """The two _WRITE_TOOL_NAMES sets are different on purpose (enforce
    vs classify) but one-directional: a tool the role denylists enforce
    must never fall out of the classifier's notion of a write.  Was FALSE
    at the baseline over make_plot -- a delegate's plot write scored less
    than its file write in verification evidence.  Rejected alternative:
    merging the two symbols (classification would then label foreign
    backend spellings as this executor's writes)."""
    missing = set(ac._WRITE_TOOL_NAMES) - set(SA._WRITE_TOOL_NAMES)
    assert missing == set(), (
        f"enforced writes invisible to classification: {sorted(missing)}")


def test_the_lead_surface_is_exactly_the_lead_module():
    """lead.py READ_ONLY_TOOLS | STATE_CHANGING_TOOLS is the source of the
    nine swarm_* names; the freezer in api_client must equal it, and the
    test's own copy must be that read.  Renaming a lead tool without the
    freezer now fails here.  Passes at both ends today -- a consistency
    pin, not a claim about the change."""
    expected = set(lead.READ_ONLY_TOOLS) | set(lead.STATE_CHANGING_TOOLS)
    assert set(ac._SWARM_LEAD_TOOL_NAMES) == expected
    assert set(_LEAD_NAMES) == expected


def test_every_documented_tool_name_is_a_known_tool():
    """The advertised catalogue's name fields must all be reachable in
    the unified answer.  Different shapes (schema dicts vs names) can
    only be pinned, never aliased.  Passes at both ends -- states the
    relation the doc-list-based judgement relied on."""
    docs = {
        str(t.get("name") or (t.get("function") or {}).get("name") or "")
        for t in ac._DOC_TOOLS_OPENAI
    }
    docs.discard("")
    missing = sorted(docs - set(ac.known_tool_names()))
    assert missing == [], f"documented but unknown: {missing}"


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(pytest.main([__file__, "-q"]))