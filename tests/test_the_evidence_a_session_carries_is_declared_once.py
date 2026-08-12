"""Evidence was exported by a hand-written list, and the list was short.

``_export_evidence`` named six fields. The engine keeps eight. The two it
forgot are both guards:

``_observed_ledger_available`` (engine.py:517, read at 2447 as
``ledger_available``) came back False after every resume while
``_last_observed_files`` came back fully populated -- so the location
claim guard was switched off for the first answer of every resumed
session, holding a ledger it had been told did not exist.

``role_test_evidence`` (engine.py:428, consumed at tab_agent.py:2863) is
the ledger behind "a PASS with NO recorded test execution is downgraded
to UNTESTED". It came back empty while ``role_verdicts`` came back
restored, so a restored test_agent PASS faced no veto.

Neither omission is interesting on its own. The shape is: a list of
fields maintained by hand, next to a set of fields maintained by need,
and nothing that notices when they drift apart. It had already drifted
twice.

So the fields are declared once, and export, restore and reset all read
that declaration. The test that matters is the last one here: it walks
the engine's own evidence attributes and fails when one is not covered
-- which is the failure that happened twice, made impossible rather than
fixed twice.

``saved_at`` is stamped as part of the same change. Nothing acts on it
yet, deliberately: a resumed session that asserts "the tests pass now"
is currently believed on a ledger recorded against a tree that may have
moved since, and the fix for that is to say where the evidence came
from, not to throw it away -- discarding it is what produced the false
"unverified" caveats that made exporting necessary in the first place.
"""

from __future__ import annotations

import pytest

from delfin.agent.engine import AgentEngine


def _engine() -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    eng.role_verdicts = {}
    eng.role_test_evidence = {}
    eng._trimmed_chars_since_floor = 0
    eng._last_observed_files = set()
    eng._observed_ledger_available = False
    eng._exec_commands_session = []
    eng._session_tool_names = set()
    eng._delegation_satisfied = False
    return eng


def _populated() -> AgentEngine:
    eng = _engine()
    eng.role_verdicts = {"test_agent": {"status": "pass"}}
    eng.role_test_evidence = {"test_agent": ["pytest -q: 12 passed"]}
    eng._trimmed_chars_since_floor = 4096
    eng._last_observed_files = {"delfin/agent/engine.py"}
    eng._observed_ledger_available = True
    eng._exec_commands_session = ["pytest -q"]
    eng._session_tool_names = {"read_file", "bash"}
    eng._delegation_satisfied = True
    return eng


def _round_trip(eng: AgentEngine) -> AgentEngine:
    data = {"evidence": eng._export_evidence()}
    out = _engine()
    out._restore_evidence(data)
    return out


# ---------------------------------------------------------------------------
# The two that were dropped
# ---------------------------------------------------------------------------

def test_the_observed_ledger_flag_survives():
    """False here means the location guard runs disabled while holding a
    populated ledger."""
    assert _round_trip(_populated())._observed_ledger_available is True


def test_the_test_evidence_survives():
    """Empty here means a restored PASS has nothing that can veto it."""
    assert _round_trip(_populated()).role_test_evidence == {
        "test_agent": ["pytest -q: 12 passed"]}


def test_a_flag_that_was_false_stays_false():
    """Restoring must not invent a ledger either."""
    assert _round_trip(_engine())._observed_ledger_available is False


# ---------------------------------------------------------------------------
# ...and everything that already worked still does
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("attr,expected", [
    ("_last_observed_files", {"delfin/agent/engine.py"}),
    ("_exec_commands_session", ["pytest -q"]),
    ("_session_tool_names", {"read_file", "bash"}),
    ("_delegation_satisfied", True),
    ("role_verdicts", {"test_agent": {"status": "pass"}}),
    ("_trimmed_chars_since_floor", 4096),
])
def test_the_existing_fields_round_trip(attr, expected):
    assert getattr(_round_trip(_populated()), attr) == expected


def test_a_session_saved_before_this_existed_still_loads():
    eng = _engine()
    eng._restore_evidence({})
    assert eng._last_observed_files == set()
    assert eng.role_test_evidence == {}


def test_malformed_evidence_does_not_raise():
    eng = _engine()
    eng._restore_evidence({"evidence": {"observed_files": 7,
                                        "role_test_evidence": "nonsense"}})
    assert isinstance(eng.role_test_evidence, dict)


# ---------------------------------------------------------------------------
# The drift is made impossible, not fixed twice
# ---------------------------------------------------------------------------

def test_every_evidence_field_is_covered_by_the_declaration():
    """The failure that happened twice: a field added to the engine and
    not to the export. It is caught here rather than discovered in a
    resumed session months later."""
    declared = {spec.attr for spec in AgentEngine._EVIDENCE_FIELDS}
    kept = {
        "role_verdicts", "role_test_evidence", "_trimmed_chars_since_floor",
        "_last_observed_files", "_observed_ledger_available",
        "_exec_commands_session", "_session_tool_names",
        "_delegation_satisfied",
    }
    assert kept <= declared, f"not carried across a resume: {kept - declared}"


def test_the_export_is_built_from_the_declaration():
    keys = set(_populated()._export_evidence())
    for spec in AgentEngine._EVIDENCE_FIELDS:
        assert spec.key in keys, f"{spec.attr} is declared but not exported"


def test_the_evidence_is_stamped_with_when_it_was_taken():
    """Nothing reads it yet; a guard that wants to say 'this came from an
    earlier session' needs it to exist first."""
    assert _populated()._export_evidence().get("saved_at")
