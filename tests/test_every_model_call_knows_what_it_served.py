"""Every model call knows which molecule, workflow and task it served.

The metric rows say *how* a call went and *who answered* it. They cannot say
**what it was for**. ``turn_metrics`` and ``agent_metrics`` record provider,
model, role, tokens and wall time, and ``budget_ledger`` charges a delegate's
spend to the calling session — so a session of 200 turns over three molecules
and two workflow runs is one undifferentiated block, and "what did this
molecule cost?" has no answer in the data.

This adds an attribution context — molecule id, workflow id, task id, role,
subagent name, parent session — set where the work is known and inherited by
every model call made under it, including a delegated subagent and the
request a tool triggers. Four properties matter and are tested here:

* **Inherited, not re-supplied.** A child row carries the parent's
  ``molecule_id`` and ``workflow_id`` and a *different* role. A delegation
  that loses them makes the child's spend unattributable, which is the exact
  hole this fills — so a child row that drops the ids is a red test, not a
  cosmetic miss.
* **Never guessed.** Where nothing is known the columns stay ``None``. An
  invented id is worse than a null: it silently bills one molecule for
  another's compute, and nothing downstream can tell the two apart.
* **Identifiers only.** The stores are kept for months. A path, a SMILES
  string or a user name must not reach them, and truncating them is not a fix
  — a cut path still leaks its first directories. Such a value is dropped.
* **Additive.** The columns are null-defaulted and appended, so a row written
  before them still reads (``test_an_old_row_still_reads``, the control that
  passes on both trees).

Red on the unpatched tree: the rows have no such columns, and
:mod:`delfin.agent.attribution` does not exist.
"""

from __future__ import annotations

import json
import textwrap
import threading
import types
from unittest.mock import MagicMock, patch

import pytest

from delfin.agent import agent_metrics
from delfin.agent import state_paths
from delfin.agent import turn_metrics as tm
from delfin.agent.api_client import StreamEvent

COLUMNS = ("molecule_id", "workflow_id", "task_id",
           "agent_role", "subagent_name", "parent_session")


@pytest.fixture
def home(tmp_path, monkeypatch):
    """Move both metric stores into a temp home (state_paths ratchet).

    The ambient run and member are cleared for the same reason ``HOME`` is
    moved: ``attribution.from_env`` reads them, so a session that starts
    this suite inside a swarm member would put that member's ids on every
    row and these tests would report the harness, not the code.

    Clearing the variables is not enough on its own. ``cli.main`` adopts them
    into a module global that outlives the call (it must: it has to survive
    the thread hop), so an earlier test file that called ``main()`` inside a
    member leaves that member's identity installed no matter what the
    environment says. Measured: ``test_agent_cli.py`` followed by this file,
    started from inside a member, gave 7 failed / 25 passed before this line
    and 32 passed after -- with no change to the code under test.
    """
    from delfin.agent import attribution

    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    for name in ("DELFIN_SWARM_RUN", "DELFIN_AGENT_LABEL"):
        monkeypatch.delenv(name, raising=False)
    attribution.reset_process()
    with state_paths.RedirectedUserState(h):
        monkeypatch.setattr(tm, "_DIR", h / ".delfin" / "turn_metrics")
        try:
            yield h
        finally:
            attribution.reset_process()


# ---------------------------------------------------------------------------
# A fake endpoint: the real client, no network and no key
# ---------------------------------------------------------------------------

def _fake_stream(text: str = "ok"):
    content = types.SimpleNamespace(
        usage=None,
        choices=[types.SimpleNamespace(
            delta=types.SimpleNamespace(content=text, tool_calls=None,
                                        reasoning_content=None),
            finish_reason="stop")])
    usage = types.SimpleNamespace(
        usage=types.SimpleNamespace(
            prompt_tokens=1000, completion_tokens=5,
            prompt_tokens_details=types.SimpleNamespace(cached_tokens=900)),
        choices=[])

    class _Stream:
        def __iter__(self):
            return iter((content, usage))

        def close(self):
            pass

    return _Stream()


class _Completions:
    """Answers both shapes a turn issues: a streamed round and a verdict call.

    The streamed request is the turn's own model call; the non-streamed one is
    the action classifier's, which the client issues through
    ``client.with_options(...)``. A stub that serves only the first silently
    tests only half the paths that need attributing.
    """

    @staticmethod
    def create(**kw):
        if kw.get("stream"):
            return _fake_stream()
        return types.SimpleNamespace(choices=[types.SimpleNamespace(
            message=types.SimpleNamespace(content="ALLOW: routine"),
            finish_reason="stop")])


class _Endpoint:
    chat = types.SimpleNamespace(completions=_Completions())

    def with_options(self, **kw):
        return self


def _steerable(client):
    """The per-turn state ``stream_message`` expects a real client to have."""
    client._steer_msgs = []
    client._steer_pending = []
    client._turn_rows_base = 0
    client._turn_model_switches = []
    client._tool_context = None
    return client


def _client(session: str, model: str = "kit.fake-model"):
    """A real ``OpenAIClient`` over a fake endpoint."""
    from delfin.agent import api_client as ac

    c = ac.OpenAIClient.__new__(ac.OpenAIClient)
    c.client = _Endpoint()
    c.model = model
    c._provider = "kit"
    c._base_url = "https://endpoint.invalid/api/v1"
    c._api_key = "x"  # secret-guard: allow (synthetic DELFIN test sample)
    c.effort = ""
    c._permissions = None
    c.on_model_switched = None
    c._steer_lock = threading.Lock()
    c._steer_queue = []
    c._run_notes = []
    c._stop_flag = False
    c.metrics_session = session
    return _steerable(c)


def _one_request(client, session):
    """Drive the real ``stream_message`` once and return the row it wrote."""
    for ev in client.stream_message(
            messages=[{"role": "user", "content": "hi"}],
            system="t", max_tokens=16):
        if getattr(ev, "type", "") == "message_start":
            break
    rows = tm.read_requests(session)
    assert rows, "the model call wrote no request row at all"
    return rows[-1]


# ---------------------------------------------------------------------------
# 1. The columns exist on the rows the product writes
# ---------------------------------------------------------------------------

def test_a_request_row_carries_the_attribution_columns(home):
    """A request row names what it served — here, honestly nothing.

    Driven through ``stream_message`` (the model call itself), not through the
    writer helper. The columns must be PRESENT and null: an absent key means
    "this recorder had no such field", which is a different statement from
    "nothing was known" — the distinction turn_metrics already draws for
    ``first_event_ms`` and ``retries``.
    """
    entry = _one_request(_client("plain"), "plain")
    for col in COLUMNS:
        assert col in entry, f"request row is missing the {col!r} column"
        assert entry[col] is None


def test_a_turn_row_carries_the_attribution_columns(home, agent_tree):
    """The engine's own turn row carries them too (turn_metrics.record)."""
    eng = _engine(agent_tree, _quiet_stream)
    eng.stream_response("hi")
    entry = _last_turn(eng)
    for col in COLUMNS:
        assert col in entry, f"turn row is missing the {col!r} column"


def test_an_agent_metrics_row_carries_the_attribution_columns(home):
    """``agent_metrics`` states its schema is forward-compatible; a new
    null-defaulted column must appear in the JSONL record it writes."""
    agent_metrics.record_turn(agent_metrics.TurnMetrics(session_id="s"))
    (rec,) = agent_metrics.read_turns()
    for col in COLUMNS:
        assert col in rec, f"agent_metrics row is missing {col!r}"
        assert rec[col] is None


# ---------------------------------------------------------------------------
# 2. Set where the work is known
# ---------------------------------------------------------------------------

def test_a_workflow_run_stamps_the_ids_it_knows(home, agent_tree):
    """A workflow run knows its molecule, workflow and task, and every call
    made under it inherits them without being told again."""
    from delfin.agent import attribution

    with attribution.attribution_scope(molecule_id="mol-7",
                                       workflow_id="esd-chain",
                                       task_id="t1"):
        req = _one_request(_client("wf"), "wf")
        eng = _engine(agent_tree, _quiet_stream)
        eng.stream_response("hi")
        turn = _last_turn(eng)

    assert req["molecule_id"] == "mol-7"
    assert req["workflow_id"] == "esd-chain"
    assert req["task_id"] == "t1"
    assert turn["molecule_id"] == "mol-7"
    assert turn["workflow_id"] == "esd-chain"

    # ...and leaving the scope clears it: the next call knows nothing again.
    assert _one_request(_client("after"), "after")["molecule_id"] is None


def test_the_role_that_ran_is_recorded_apart_from_the_who(home, agent_tree):
    """``agent_role`` says who ran, and a delegate is not the parent.

    The parent's own call and the delegate's differ in ``agent_role`` while
    agreeing on the ids — that pair is what makes a workflow run's cost
    readable: one molecule, two owners of the spend.
    """
    from delfin.agent import attribution

    with attribution.attribution_scope(molecule_id="mol-11",
                                       role="session_manager"):
        eng = _engine(agent_tree, _quiet_stream)
        eng.stream_response("hi")
        row = _last_turn(eng)
    assert row["molecule_id"] == "mol-11"
    assert row["agent_role"] == "session_manager"
    # A turn nobody wrapped keeps its ID columns null. The engine knows its own
    # role and session and no molecule, and inventing one would bill this turn
    # to somebody else's compute.
    plain = _engine(agent_tree, _quiet_stream)
    plain.stream_response("hi")
    bare = _last_turn(plain)
    assert bare["molecule_id"] is None
    assert bare["workflow_id"] is None
    assert bare["task_id"] is None


def test_a_tool_triggered_call_is_attributed_without_being_told(home):
    """A call made deep inside a tool is attributed.

    This is the part an argument cannot reach: the classifier request is
    issued inside the client and its writer never heard of attribution. Nothing
    here passes a parameter to it.
    """
    from delfin.agent import attribution

    client = _client("tool")
    with attribution.attribution_scope(molecule_id="mol-9"):
        assert client._classifier_complete("say hello", "kit.fake-model")
    rows = [r for r in tm.read_requests("tool")
            if r.get("purpose") == "classifier"]
    assert rows, "the classifier request was not recorded at all"
    assert rows[-1]["molecule_id"] == "mol-9"


# ---------------------------------------------------------------------------
# 3. Inheritance through a delegation
# ---------------------------------------------------------------------------

def test_a_delegate_keeps_the_parents_ids_and_changes_its_role(home):
    """The product-level claim: parent and child rows share molecule and
    workflow ids and carry different roles.

    Driven through ``run_subagent`` — the delegation entry point itself — with
    the fake endpoint standing in for the model.

    The delegation deliberately runs AFTER the wrapping scope has closed. An
    earlier version ran it inside the still-open ``attribution_scope``, where
    the ambient context supplies molecule/workflow/task whether or not
    ``run_subagent`` inherits anything -- so a ``run_subagent`` that dropped
    the parent's ids entirely would still have produced a green test. With the
    scope closed, the client binding is the only route left, which is the
    claim: the parent's client is copied into the delegate and the ids travel
    with it.
    """
    from delfin.agent import attribution, subagents

    parent = _client("deleg")
    with attribution.attribution_scope(molecule_id="mol-3",
                                       workflow_id="opt-freq",
                                       task_id="t7",
                                       role="session_manager"):
        attribution.bind(parent)          # what the turn wrappers do at turn start

    # Outside the scope: nothing ambient is left to inherit.
    assert attribution.current() is None or attribution.current().is_empty(), (
        "the scope is still open, so this test could not see the delegation")
    assert attribution.bound_to(parent) is not None, (
        "the parent client lost its binding, so the ids could not travel")

    # The run's own call, then the delegation it hands work to: two rows, one
    # owner each, for the same piece of work.
    _one_request(parent, "deleg")
    res = subagents.run_subagent(
        subagent_type="explore", description="map it", prompt="go",
        parent_client=parent, parent_perms=None)
    assert res.error == "", f"delegation failed: {res.error}"

    rows = tm.read_requests("deleg")
    assert rows, "no request rows were written for the delegation"
    parents = [r for r in rows if r.get("subagent_name") is None]
    kids = [r for r in rows if r.get("subagent_name") == "explore"]
    assert parents, f"no parent-attributed row in {rows}"
    assert kids, f"no child-attributed row in {rows}"

    child = kids[-1]
    assert child["molecule_id"] == "mol-3"       # ids survive the delegation
    assert child["workflow_id"] == "opt-freq"
    assert child["task_id"] == "t7"
    assert child["subagent_name"] == "explore"   # and the row says who ran
    assert child["agent_role"] != parents[-1]["agent_role"]
    assert child["parent_session"] == "deleg"


def test_a_delegate_inherits_across_a_thread_boundary(home):
    """A background delegate runs on a NEW thread, and a context var does not
    cross that boundary.

    DELFIN spawns a background sub-agent with ``threading.Thread`` and calls
    the runner *inside* that thread (``api_client.py``: ``_bg_run`` ->
    ``perms.subagent_runner``), so an attribution living only in the ambient
    context would be ``None`` in the child and its spend would land
    unattributed. The first assertion pins that hazard down; the second is the
    regression this test exists for.
    """
    from delfin.agent import attribution, subagents

    parent = _client("bg")
    box: dict = {}

    with attribution.attribution_scope(molecule_id="mol-5",
                                       workflow_id="scan"):
        attribution.bind(parent)                  # the turn wrappers do this per turn

        def body():
            # A fresh thread: the caller's context is gone by definition.
            box["ambient"] = attribution.current()
            sub = subagents.run_subagent(
                subagent_type="explore", description="d", prompt="p",
                parent_client=parent, parent_perms=None)
            box["err"] = sub.error
            box["rows"] = tm.read_requests("bg")

        t = threading.Thread(target=body)
        t.start()
        t.join(120)

    assert not t.is_alive(), "the delegate thread hung"
    assert box["ambient"] is None, (
        "a new thread must not inherit the caller's context var — if this "
        "ever passes, the client-bound path below is dead code")
    assert box["err"] == "", f"delegation failed: {box['err']}"
    kids = [r for r in box["rows"] if r.get("subagent_name") == "explore"]
    assert kids, f"no child row in {box['rows']}"
    assert kids[-1]["molecule_id"] == "mol-5", (
        "the child lost the parent's molecule id across the thread boundary")
    assert kids[-1]["workflow_id"] == "scan"


# ---------------------------------------------------------------------------
# 4. Never guess, and never let a name or path through
# ---------------------------------------------------------------------------

def test_nothing_known_stays_null_rather_than_an_invented_id(home):
    """A call with no context records nulls — not a default, not a guess."""
    entry = _one_request(_client("noguess"), "noguess")
    assert entry["molecule_id"] is None
    assert entry["workflow_id"] is None
    assert entry["task_id"] is None
    assert entry["parent_session"] is None


@pytest.mark.parametrize("bad", [
    "/home/user/molecule.xyz",
    r"M:\scratch\job1",
    "CC(=O)Oc1ccccc1C(=O)O",          # a SMILES string
    "alice@institute.example",         # an account
    "../../etc/passwd",
    "benzene dimer complex",           # a molecule NAME
    "",
])
def test_a_path_or_smiles_is_never_stored_as_an_identifier(bad):
    """Identifiers only: the stores are kept for months.

    A value that does not look like an opaque identifier is dropped to
    ``None`` rather than truncated — a cut path still leaks its leading
    directories.
    """
    from delfin.agent import attribution

    with attribution.attribution_scope(molecule_id=bad) as got:
        cur = attribution.current()
    assert cur is None or cur.molecule_id is None, (
        f"{bad!r} was accepted as an identifier")
    # Not ``got is not None``: ``__enter__`` returns an Attribution instance on
    # every path (``_merged`` always builds one), so that held even when the
    # scope dropped the id -- the state this test is about. What the scope
    # must never do is raise, and an exception here fails the test anyway; what
    # is worth asserting is that the rejected id did not survive anywhere.
    assert got.molecule_id is None, f"{bad!r} survived on the scope's record"


def test_an_identifier_that_is_one_is_kept():
    """The counterpart: a real opaque id must survive, or the feature is
    decoration. These are the shapes DELFIN actually uses."""
    from delfin.agent import attribution

    with attribution.attribution_scope(molecule_id="mol_2026-12",
                                       workflow_id="esd-chain",
                                       task_id="t3a"):
        cur = attribution.current()
    assert cur is not None
    assert (cur.molecule_id, cur.workflow_id, cur.task_id) == (
        "mol_2026-12", "esd-chain", "t3a")


def test_a_child_cannot_forget_the_ids_it_inherited():
    """The child record is built by copying the parent's ids, so a caller
    cannot silently drop them by naming only its own role.

    Input that would make this fail: a ``child_of`` that started from an empty
    record instead of the parent's.
    """
    from delfin.agent import attribution

    with attribution.attribution_scope(molecule_id="mol-1",
                                       workflow_id="wf-1",
                                       task_id="t1"):
        parent = attribution.current()
        kid = attribution.child_of(parent, role="subagent",
                                   subagent_name="explore",
                                   parent_session="sess-1")
    assert kid.molecule_id == "mol-1"
    assert kid.workflow_id == "wf-1"
    assert kid.task_id == "t1"
    assert kid.subagent_name == "explore"
    assert kid.parent_session == "sess-1"
    # ...and the parent record is untouched: it is immutable.
    assert parent.subagent_name is None


def test_the_turn_binds_its_attribution_to_the_client_it_hands_out(home,
                                                                   agent_tree):
    """The turn binds the attribution to its own client WHILE it runs.

    This is the caller that makes the client-bound path real rather than a
    helper nobody reaches (rule 7). It is needed because a delegate is reached
    two ways: a foreground one runs in the caller's context, but a BACKGROUND
    one runs on a thread that has no context at all and receives the parent's
    client by shallow copy -- so only the binding carries the ids to it.

    Driven through ``repl.run_turn``, the turn every interactive and
    service-worker session ends with. The binding deliberately lives at that
    entry point and not inside ``engine.stream_response``: the engine is on the
    protected list, and two writable entry points can do the same job. Red
    before either existed: without a binding, every background delegate
    records nulls.
    """
    from delfin.agent import attribution, repl

    seen: dict = {}

    def _stream(system, messages, **kw):
        # Mid-turn, where a delegate or a tool would read it.
        seen["bound"] = attribution.bound_to(eng.client)
        yield StreamEvent(type="message_start", input_tokens=900)

    with attribution.attribution_scope(molecule_id="mol-13",
                                       workflow_id="wf-13"):
        eng = _engine(agent_tree, _stream)
        repl.run_turn(eng, "hi", sink=lambda item: None)

    got = seen.get("bound")
    assert got is not None, (
        "the turn never bound an attribution to its client, so a background "
        "delegate -- which runs on a thread with no context -- inherits nothing")
    assert got.molecule_id == "mol-13"
    assert got.workflow_id == "wf-13"
    assert got.parent_session                      # the session that paid
    # ...and it does not outlive the turn: the client is reused next turn.
    assert attribution.bound_to(eng.client) is None


# ---------------------------------------------------------------------------
# 5. Control: an old row still reads (passes before and after)
# ---------------------------------------------------------------------------

def test_an_old_row_still_reads(home):
    """Forward compatibility, the other direction.

    A row written before the columns existed must still read and still
    aggregate. This passes on the unpatched tree and after — the control that
    says the new columns are additive, not a schema break.
    """
    legacy = {"ts": 1700000000.0, "provider": "kit", "model": "m",
              "role": "session_manager", "total_ms": 1200,
              "ttft_ms": 300, "first_event_ms": 250, "output_chars": 12,
              "tool_calls": 0, "stopped": False, "error": "",
              "retries": None, "input_tokens": 10, "output_tokens": 2,
              "cached_tokens": 0}
    p = tm.metrics_path("legacy")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(legacy) + "\n", encoding="utf-8")

    (entry,) = tm.read("legacy")
    assert entry["total_ms"] == 1200
    assert entry["model"] == "m"
    assert tm.aggregate_turn_stats(window_days=36_500,
                                   dir_path=tm._DIR)["turns"] >= 1


# ---------------------------------------------------------------------------
# Engine harness (same shape as the turn-metrics engine tests)
# ---------------------------------------------------------------------------

@pytest.fixture
def agent_tree(tmp_path):
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    (lite_dir / "manifest.yaml").write_text(textwrap.dedent("""\
        pack_name: DELFIN_AGENT_LITE
        version: 1
        modes:
          - id: solo
            file: modes/solo.md
            route:
              - session_manager
        """))
    return tmp_path


def _quiet_stream(system, messages, **kw):
    yield StreamEvent(type="message_start", input_tokens=900)


def _engine(agent_tree, stream):
    from delfin.agent.engine import AgentEngine
    client = _client("eng")
    client.stream_message = MagicMock(side_effect=stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        eng = AgentEngine(repo_dir=agent_tree, backend="api", provider="kit",
                          mode="quick", pack_dir=agent_tree)
    eng.client = client
    return eng


def _last_turn(eng) -> dict:
    entries = tm.read(eng.trace_session())
    assert entries, "the turn should have been recorded"
    return entries[-1]
