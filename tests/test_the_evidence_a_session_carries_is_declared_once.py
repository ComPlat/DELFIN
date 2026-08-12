"""Evidence was exported by a hand-written list, and the list was short.

``_export_evidence`` named six fields. The engine kept eight. The two it
forgot were both guards:

``_observed_ledger_available`` (read as ``ledger_available``) came back
False after every resume while ``_last_observed_files`` came back fully
populated -- so the location claim guard was switched off for the first
answer of every resumed session, holding a ledger it had been told did
not exist.

``role_test_evidence`` is the ledger behind "a PASS with NO recorded test
execution is downgraded to UNTESTED". It came back empty while
``role_verdicts`` came back restored, so a restored test_agent PASS faced
no veto.

Neither omission is interesting on its own. The shape is: a list of
fields maintained by hand, next to a set of fields maintained by need,
and nothing that notices when they drift apart.

The declaration was introduced to end that. It did not, because it was
one of THREE hand-written lists:

1. the declaration itself;
2. ``reset_cycle``, which hand-cleared three of the eight fields -- so a
   brand-new conversation started holding the previous one's tool-name
   set and a delegation flag already satisfied, which suppressed the
   delegation reminder for work that had never delegated anything;
3. this test, whose coverage case compared the declaration to a literal
   set written beside it in this file. It asserted ``kept <= declared``
   where ``kept`` was a copy of ``declared``. A ninth ledger passed
   green, while the docstring above it claimed the test walked the
   engine's attributes.

So: one declaration that export, restore AND reset iterate, and a test
that enumerates the engine for real. The enumeration found the ninth and
tenth ledger on its first run -- ``_tasklist_satisfied`` and
``_session_errors``, both session-scoped, neither carried across a
resume nor cleared on a new cycle.

``saved_at`` is stamped as part of the same change. Nothing acts on it
yet, deliberately: a resumed session that asserts "the tests pass now"
is currently believed on a ledger recorded against a tree that may have
moved since, and the fix for that is to say where the evidence came
from, not to throw it away -- discarding it is what produced the false
"unverified" caveats that made exporting necessary in the first place.
"""

from __future__ import annotations

import ast
import inspect
import pathlib

import pytest

from delfin.agent.engine import AgentEngine


def _engine() -> AgentEngine:
    eng = AgentEngine.__new__(AgentEngine)
    for spec in AgentEngine._SESSION_FIELDS:
        setattr(eng, spec.attr, spec.reset())
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
    eng._tasklist_satisfied = True
    eng._session_errors = [{"key": "k", "tool": "bash", "error": "boom",
                            "count": "2"}]
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
# The ninth and tenth ledger, which the old coverage test could not see
# ---------------------------------------------------------------------------

def test_the_task_list_flag_survives():
    """Not carried, the reminder re-fired after a resume for a task list
    that already existed."""
    assert _round_trip(_populated())._tasklist_satisfied is True


def test_the_session_error_memory_survives():
    """"do NOT repeat these" is worth nothing when the list comes back
    empty on the turn most likely to repeat them."""
    assert _round_trip(_populated())._session_errors == [
        {"key": "k", "tool": "bash", "error": "boom", "count": "2"}]


# ---------------------------------------------------------------------------
# The drift is made impossible, not fixed three times
# ---------------------------------------------------------------------------

def _engine_source_tree() -> ast.ClassDef:
    source = pathlib.Path(inspect.getfile(AgentEngine)).read_text(
        encoding="utf-8")
    tree = ast.parse(source)
    return next(n for n in tree.body
                if isinstance(n, ast.ClassDef) and n.name == "AgentEngine")


def _assigned_attributes(node: ast.AST) -> set[str]:
    """Every ``self.<name>`` the given tree assigns to."""
    out: set[str] = set()
    for sub in ast.walk(node):
        targets: list[ast.expr] = []
        if isinstance(sub, ast.Assign):
            targets = list(sub.targets)
        elif isinstance(sub, (ast.AnnAssign, ast.AugAssign)):
            targets = [sub.target]
        for t in targets:
            if (isinstance(t, ast.Attribute)
                    and isinstance(t.value, ast.Name)
                    and t.value.id == "self"):
                out.add(t.attr)
    return out


# Every attribute AgentEngine keeps that is deliberately NOT carried
# across a resume, with the reason it is not. This list is the enumeration
# working: adding an attribute to the engine now forces the question
# "does a resumed session need this?" to be answered here, in one line,
# instead of being answered by silence three releases later.
_NOT_CARRIED = {
    # -- construction: rebuilt from the launch arguments, never restored --
    "repo_dir", "backend", "provider", "effort", "client", "loader",
    "_lock", "_agent_workspace_dir", "_is_delfin_workspace", "_distiller",
    "_context_tracker", "_progressive_disclosure", "_live_state",
    "_turn_gate", "_trace_id",
    # -- restored explicitly, outside the declaration, with side effects --
    "mode",             # _load_mode() reloads the role definitions
    "mode_description",  # a product of _load_mode
    "route",             # a product of _load_mode
    "messages",          # the conversation itself
    "session_id",        # also re-binds the task list
    # -- per-TURN scratch: meaningless the moment the turn ends -----------
    "_last_structured_verdict", "_last_test_evidence", "_last_turn_tools",
    "_saw_message_start", "_floor_captured_this_turn",
    "_claim_guard_active", "_claim_guard_corrected", "_exec_pending",
    "_trace_pending", "_turn_in_flight", "_turn_start_cost",
    "_cost_cap_hit", "_cost_cap_value", "_ambiguous_columns_turn",
    "_truncated_tools_turn", "_stop_requested", "_steering_delivered",
    "_steering_refreshes",
    # -- the message being answered right now -----------------------------
    #    The figure guard reads it so a number the USER typed grounds
    #    itself. It is the current turn's input; carrying it across a
    #    resume would let a figure from a previous conversation excuse
    #    one stated in this one.
    "_last_user_message",
    # -- the name of THIS turn's office figure ledger ----------------------
    #    The ledger it names lives in memory in this process, so the token
    #    means nothing to another one. Carrying it would name a ledger that
    #    does not exist — and the empty-ledger reading is the strict one,
    #    so a resumed turn would caveat figures its own tools produced.
    "_figure_ledger_token",
    # -- the correction budget, spent at most once per turn ---------------
    #    Kept apart from the VERDICT (_claim_guard_corrected) so the green
    #    "answer corrected" line comes from a re-scan rather than from
    #    having tried. Neither is worth carrying: a resumed session has
    #    not spent anything yet.
    "_claim_guard_spent",
    # -- the stop's owner, and the counter it is stamped from -------------
    #    A stop belongs to the turn that asked for it — that is what keeps
    #    the next message from erasing it. Nothing about that survives the
    #    process: on resume there is no turn in flight to own a stop, and
    #    a serial carried across would let a stale id match a fresh turn.
    "_turn_serial", "_turn_id", "_stop_owner_turn",
    # -- the last turn's diagnosis, shown once and cleared next turn ------
    "last_empty_turn",
    # -- the live system prompt: its TEXT carries the injected memory and
    #    is deliberately never written to disk. Its SIZE is carried, as
    #    _system_prompt_chars, which is declared.
    "last_system_prompt",
    # -- resolved per process from the active model / settings ------------
    "context_window_tokens", "auto_compact_pct", "_active_capabilities",
    # -- reports about the last operation, not state to continue from -----
    "last_compaction_info",   # persisted by the callers, not the exporter
    "last_restore_report",
    # -- re-armed once per process on purpose ------------------------------
    "_foreign_tasks_shown",   # a resumed session may re-show the notice
    # The "no USD ceiling is in force" notice. The COUNTS behind it are
    # declared session state, because a resumed run's spend figure means
    # nothing without them; this flag is only "has the agent been told",
    # and the agent of a resumed session has not been.
    "_unmeasured_budget_notice_shown",
    "_prompt_session_serial",
}


def test_every_attribute_the_engine_keeps_is_classified():
    """The test that matters, and the one that used to be a lie.

    It walks the engine's own attribute assignments -- not a copy of the
    declaration written beside it -- and fails when one is neither
    declared as session state nor listed above as deliberately dropped.
    """
    declared = {spec.attr for spec in AgentEngine._SESSION_FIELDS}
    found = _assigned_attributes(_engine_source_tree())
    unclassified = found - declared - _NOT_CARRIED
    assert not unclassified, (
        f"AgentEngine keeps {sorted(unclassified)}, and nothing says "
        f"whether a resumed session needs it. Add it to _SESSION_FIELDS "
        f"(carried across a resume, exported, restored and reset from one "
        f"declaration) or to _NOT_CARRIED in this file, with the reason."
    )


def test_the_classification_list_does_not_rot():
    """A name listed as deliberately dropped that the engine no longer
    keeps is a stale entry hiding the next real one."""
    found = _assigned_attributes(_engine_source_tree())
    stale = _NOT_CARRIED - found
    assert not stale, f"no longer engine attributes: {sorted(stale)}"


def test_a_declared_field_is_never_also_listed_as_dropped():
    declared = {spec.attr for spec in AgentEngine._SESSION_FIELDS}
    assert not (declared & _NOT_CARRIED)


def test_the_export_is_built_from_the_declaration():
    state = _populated().export_state()
    evidence = state.get("evidence") or {}
    for spec in AgentEngine._SESSION_FIELDS:
        where = evidence if spec.section == "evidence" else state
        assert spec.key in where, (
            f"{spec.attr} is declared but not exported into "
            f"{spec.section}")


def test_the_evidence_is_stamped_with_when_it_was_taken():
    """Nothing reads it yet; a guard that wants to say 'this came from an
    earlier session' needs it to exist first."""
    assert _populated()._export_evidence().get("saved_at")


# ---------------------------------------------------------------------------
# The third list: the reset
# ---------------------------------------------------------------------------

def test_the_reset_returns_every_declared_field_to_a_fresh_value():
    """reset_cycle hand-cleared three of eight. The five it forgot were
    carried into the next conversation."""
    eng = _populated()
    eng._reset_session_fields()
    for spec in AgentEngine._SESSION_FIELDS:
        got = getattr(eng, spec.attr)
        if spec.attr == "_run_started_at":
            continue          # a clock: fresh means "now", not a constant
        assert got == spec.reset(), (
            f"{spec.attr} survived a reset as {got!r}")


def test_a_new_conversation_does_not_inherit_the_previous_tool_names():
    """The concrete consequence: the delegation reminder was suppressed
    for a cycle that had never delegated anything."""
    eng = _populated()
    eng._reset_session_fields()
    assert eng._session_tool_names == set()
    assert eng._delegation_satisfied is False


def test_a_mode_switch_keeps_the_evidence_for_the_history_it_keeps():
    """The narrow reset: role-keyed values go, the ledgers for the
    conversation that is being kept stay."""
    eng = _populated()
    eng._reset_session_fields(role_scoped_only=True)
    assert eng.role_verdicts == {}
    assert eng.current_role_index == 0
    assert eng._exec_commands_session == ["pytest -q"]
    assert eng._last_observed_files == {"delfin/agent/engine.py"}


def test_the_reset_path_hand_clears_nothing_beside_the_declaration():
    """The third hand-written list, caught mechanically."""
    cls = _engine_source_tree()
    reset = next(n for n in cls.body
                 if isinstance(n, ast.FunctionDef) and n.name == "reset_cycle")
    touched = _assigned_attributes(reset)
    # Names reset_cycle is allowed to set itself: they are not session
    # LEDGERS, they are the identity of the new cycle.
    allowed = {"messages", "session_id", "_foreign_tasks_shown",
               "_prompt_session_serial", "_stop_requested"}
    declared = {spec.attr for spec in AgentEngine._SESSION_FIELDS}
    assert not (touched & declared), (
        f"reset_cycle hand-clears {sorted(touched & declared)} beside the "
        f"declaration that already resets it")
    assert touched <= allowed, f"unexpected: {sorted(touched - allowed)}"
