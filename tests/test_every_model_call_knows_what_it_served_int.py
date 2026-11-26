"""Every model call knows which molecule, workflow, task and agent it served.

The stores say *how* a call went (tokens, ttft, tool calls, wall time) and
*who answered* it (provider, model, role). They could not say **what the call
was for**, so "the cost of this molecule", "the cost of this workflow" and
"what did this delegate spend on the task it was given" had no answer in the
data.

This file covers the part that makes the feature a capability rather than a
module: the attribution is set at a REAL ENTRY POINT and reaches the rows.

The four properties, and where each is measured:

* **Set where the work is known, inherited below it.** ``repl.run_turn`` is the
  turn every interactive and service-worker session ends with; bound there, a
  workflow-run scope above it is what the turn row records. A delegate keeps
  the ids and changes only who ran — driven through ``subagents.run_subagent``,
  the real delegation path, not a helper.
* **A background delegate keeps them.** ``api_client`` starts one on a
  ``threading.Thread`` and calls the runner from inside it; a new thread has an
  empty context, so only the binding on the client reaches it. Measured on a
  real thread.
* **A process DELFIN started reads its own identity.** A swarm member gets
  ``DELFIN_SWARM_RUN``/``DELFIN_AGENT_LABEL`` and a fresh context; ``cli.main``
  adopts them, so its rows name the run and the member.
* **Never guessed, identifiers only.** Nothing bound -> six nulls. A path, a
  SMILES string or an e-mail is dropped, never truncated.

Each test states the input that would make it fail.
"""

from __future__ import annotations

import json
import threading
import types

import pytest

from delfin.agent import agent_metrics
from delfin.agent import state_paths
from delfin.agent import turn_metrics as tm

COLUMNS = ("molecule_id", "workflow_id", "task_id",
           "agent_role", "subagent_name", "parent_session")


@pytest.fixture(autouse=True)
def _clean_attribution():
    """Every test starts with no attribution in force.

    ``attribution.set_attribution`` MERGES onto the current context by design
    (that is how a nested scope adds ids), so a test that installs a process
    attribution would otherwise leak it into the next test's "nothing known"
    case and make that case pass for the wrong reason. Resetting between tests
    is what keeps the null case an actual measurement.

    Tolerates the module being ABSENT rather than skipping: on the baseline
    there is no ``delfin.agent.attribution`` at all, and each test must then
    fail on its own assertion — naming what is missing — instead of erroring
    in a fixture, which would say only "setup failed".
    """
    try:
        from delfin.agent import attribution
    except ImportError:
        yield                     # the baseline: each test fails by name
        return
    # The process-level identity is a module global (it has to outlive the
    # thread hop), so it must be cleared here as well: otherwise a test that
    # adopts a member's run leaves a later "nothing known" case attributed,
    # and that case would pass for the wrong reason.
    attribution.reset_process()
    try:
        with attribution.attribution_scope(inherit=False):
            yield
    finally:
        attribution.reset_process()


@pytest.fixture()
def home(tmp_path, monkeypatch):
    """Both metric stores and the agent's state under a temp home."""
    h = tmp_path / "home"
    h.mkdir()
    monkeypatch.setenv("HOME", str(h))
    with state_paths.RedirectedUserState(h):
        monkeypatch.setattr(tm, "_DIR", h / ".delfin" / "turn_metrics")
        yield h


# ---------------------------------------------------------------------------
# A fake endpoint: the real client and the real turn, no network and no key
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
    """Answers both shapes a turn issues: a streamed round and a verdict call."""

    def create(self, **kw):
        if kw.get("stream"):
            return _fake_stream()
        return types.SimpleNamespace(choices=[types.SimpleNamespace(
            message=types.SimpleNamespace(content="ALLOW: routine"),
            finish_reason="stop")])


class _Endpoint:
    chat = types.SimpleNamespace(completions=_Completions())

    def with_options(self, **kw):
        return self


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
    c._steer_msgs = []
    c._steer_pending = []
    c._turn_rows_base = 0
    c._turn_model_switches = []
    c._tool_context = None
    return c


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


class _FakeEngine:
    """The shape ``repl.run_turn`` drives, with a real client underneath.

    Not a stand-in for the engine's behaviour: only the attributes ``run_turn``
    reads. The model call itself goes through the real client, so the row this
    produces is written by the real writer from a real request.
    """

    def __init__(self, session: str, role: str = "session_manager"):
        self._session = session
        self._role = role
        self.client = _client(session)
        self.model = self.client.model
        self.mode = "solo"
        self.token_usage = {"input": 0, "output": 0, "cached": 0}
        self.cost_usd = 0.0
        self.route = [role]
        self.current_role_index = 0

    @property
    def current_role(self) -> str:
        return self._role

    def trace_session(self) -> str:
        return self._session

    def stream_response(self, **kw):
        on_token = kw.get("on_token")
        ev = self._turn()
        if on_token:
            on_token(ev)
        self.token_usage["input"] += 1000
        self.token_usage["output"] += 5
        return ev

    def _turn(self) -> str:
        for ev in self.client.stream_message(
                messages=[{"role": "user", "content": "hi"}],
                system="t", max_tokens=16):
            if getattr(ev, "type", "") == "message_start":
                break
        return "ok"

    def estimate_context_tokens(self) -> int:
        return 0

    def delegate_spend(self):
        return agent_metrics.DelegateLedger()

    def get_status(self) -> dict:
        return {"input_tokens": self.token_usage["input"],
                "output_tokens": self.token_usage["output"],
                "total_cost_usd": self.cost_usd}


def _chat_rows(home) -> list[dict]:
    p = home / ".delfin" / "agent_metrics.jsonl"
    if not p.exists():
        return []
    return [json.loads(line) for line in
            p.read_text(encoding="utf-8").splitlines() if line.strip()]


@pytest.fixture()
def agent_tree(tmp_path):
    """A minimal pack the real ``AgentEngine`` can be built against."""
    lite_dir = tmp_path / "pack_lite"
    modes = lite_dir / "modes"
    modes.mkdir(parents=True)
    (modes / "solo.md").write_text("# solo mode")
    (lite_dir / "manifest.yaml").write_text(
        "pack_name: DELFIN_AGENT_LITE\n"
        "version: 1\n"
        "modes:\n"
        "  - id: solo\n"
        "    file: modes/solo.md\n"
        "    route:\n"
        "      - session_manager\n")
    return tmp_path


# ---------------------------------------------------------------------------
# 1. The turn every session ends with carries the attribution
# ---------------------------------------------------------------------------

def test_the_entry_turn_binds_its_attribution_to_the_client(home):
    """``repl.run_turn`` binds what the turn serves onto its client.

    The input that would make this fail: a client whose binding is still
    ``None`` after one turn, which is what happens when the scope is only
    ambient — and then a background delegate, whose thread has no context,
    records nothing.
    """
    from delfin.agent import attribution, repl

    engine = _FakeEngine("entry")
    seen: dict = {}
    real_turn = engine._turn

    def _probe():
        seen["bound"] = attribution.bound_to(engine.client)
        return real_turn()

    engine._turn = _probe
    with attribution.attribution_scope(molecule_id="mol-entry",
                                       workflow_id="wf-entry"):
        repl.run_turn(engine, "hi", sink=lambda item: None)

    got = seen.get("bound")
    assert got is not None, (
        "run_turn never bound an attribution to the client it uses, so a "
        "background delegate inherits nothing (no attribution module on the "
        "base: this is the red-on-base assertion)")
    assert got.molecule_id == "mol-entry"
    assert got.workflow_id == "wf-entry"
    assert got.role == "session_manager", (
        "the engine knows the role it is running as; the row must say it")
    assert got.parent_session == "entry"
    # ...and it does not outlive the turn: the client is reused next turn.
    assert attribution.bound_to(engine.client) is None, (
        "a binding left on the client would bill the NEXT turn to this one's "
        "molecule")


def test_the_turn_row_and_the_agent_row_carry_the_same_ids(home, agent_tree):
    """The per-turn store and the agent_metrics line agree on what was served.

    Driven through ``repl.run_turn`` over a REAL ``AgentEngine``. The turn row
    is written by the engine (``turn_metrics.record``) and the agent_metrics
    line by ``record_chat_turn`` at the end of the same function — the two
    writers that used to know nothing about the work. A fake engine would only
    prove that ``run_turn`` calls something; and the REQUEST row is measured
    against the real client in
    ``test_the_entry_turn_binds_its_attribution_to_the_client``, because it is
    written inside ``stream_message`` itself.
    """
    from unittest.mock import MagicMock, patch

    from delfin.agent import attribution, repl
    from delfin.agent.api_client import StreamEvent
    from delfin.agent.engine import AgentEngine

    client = _client("rows")

    def _stream(system, messages, **kw):
        yield StreamEvent(type="message_start", input_tokens=900)

    client.stream_message = MagicMock(side_effect=_stream)
    with patch("delfin.agent.engine.create_client", return_value=client):
        engine = AgentEngine(repo_dir=agent_tree, backend="api",
                             provider="kit", mode="quick", pack_dir=agent_tree)
    engine.client = client
    client.metrics_session = engine.trace_session()

    with attribution.attribution_scope(molecule_id="mol-2", workflow_id="wf-2",
                                       task_id="t9"):
        repl.run_turn(engine, "hi", sink=lambda item: None)

    turns = tm.read(engine.trace_session())
    assert turns, "no turn row written"
    turn = turns[-1]
    assert turn["molecule_id"] == "mol-2", "the turn row lost the molecule id"
    assert turn["workflow_id"] == "wf-2"
    assert turn["task_id"] == "t9"
    # ...and the agent_metrics line, written inside the same scope. It is
    # written AFTER the model call, so a scope that closed with the call would
    # leave it null — which is why this assertion is here.
    chat = _chat_rows(home)
    assert chat, "repl.run_turn wrote no agent_metrics line"
    assert chat[-1]["molecule_id"] == "mol-2", (
        "the agent_metrics row was written outside the attribution scope")
    assert chat[-1]["workflow_id"] == "wf-2"


# ---------------------------------------------------------------------------
# 2. A delegation inherits the ids and changes only who ran
# ---------------------------------------------------------------------------

def _delegation_rows(client, parent_types=("explore",)):
    return [r for r in tm.read_requests(getattr(client, "metrics_session", ""))
            if r.get("subagent_name") in parent_types]


def test_a_delegate_keeps_the_parents_ids_and_changes_its_role(home):
    """Parent and child rows share molecule/workflow; the roles differ.

    The input that would make this fail: a child row whose ``molecule_id`` is
    null while the parent's is set — the delegate's spend would fall outside
    "what did this molecule cost", which is the hole this feature fills.
    """
    from delfin.agent import attribution, subagents

    parent = _client("deleg-int")
    with attribution.attribution_scope(molecule_id="mol-3",
                                       workflow_id="opt-freq",
                                       role="session_manager"):
        attribution.bind(parent)
        # The run's own call first, then the delegation it hands work to: two
        # rows for the same piece of work, one owner each.
        _one_request(parent, "deleg-int")
        res = subagents.run_subagent(
            subagent_type="explore", description="d", prompt="p",
            parent_client=parent, parent_perms=None)
    assert res.error == "", f"delegation failed: {res.error}"

    rows = tm.read_requests("deleg-int")
    parents = [r for r in rows if r.get("subagent_name") is None]
    kids = [r for r in rows if r.get("subagent_name") == "explore"]
    assert parents, f"no parent-attributed row in {rows}"
    assert kids, f"no child-attributed row in {rows}"

    child = kids[-1]
    assert child["molecule_id"] == "mol-3", "the child lost the molecule id"
    assert child["workflow_id"] == "opt-freq"
    assert child["agent_role"] != parents[-1]["agent_role"], (
        "parent and child must be distinguishable by role, or a report cannot "
        "say who spent the tokens")
    assert child["parent_session"] == "deleg-int"


def test_a_delegate_on_a_new_thread_keeps_the_ids(home):
    """A BACKGROUND delegate runs on a thread with no context.

    ``api_client`` spawns it with ``threading.Thread`` and calls the runner
    inside it. The input that would make this fail: a design that only sets a
    ``ContextVar`` — the child's ``ambient`` below is then ``None`` and its row
    carries no ids.
    """
    from delfin.agent import attribution, subagents

    parent = _client("bg-int")
    box: dict = {}

    with attribution.attribution_scope(molecule_id="mol-5", workflow_id="scan"):
        attribution.bind(parent)

        def body():
            box["ambient"] = attribution.current()
            sub = subagents.run_subagent(
                subagent_type="explore", description="d", prompt="p",
                parent_client=parent, parent_perms=None)
            box["err"] = sub.error
            box["rows"] = tm.read_requests("bg-int")

        t = threading.Thread(target=body)
        t.start()
        t.join(120)

    assert not t.is_alive(), "the delegate thread hung"
    assert box["ambient"] is None, (
        "a new thread must not inherit the caller's context var — if this ever "
        "passes, the client-bound path below is dead code")
    assert box["err"] == "", f"delegation failed: {box['err']}"
    kids = [r for r in box["rows"] if r.get("subagent_name") == "explore"]
    assert kids, f"no child row in {box['rows']}"
    assert kids[-1]["molecule_id"] == "mol-5", (
        "the child lost the parent's molecule id across the thread boundary")
    assert kids[-1]["workflow_id"] == "scan"


# ---------------------------------------------------------------------------
# 3. A process DELFIN started knows which run it works for
# ---------------------------------------------------------------------------

def test_a_swarm_member_process_adopts_its_run_and_member():
    """``cli.main`` reads the identity a member process is launched with.

    The input that would make this fail: ``from_env`` returning an empty
    record for a process that was given ``DELFIN_SWARM_RUN`` — the member's
    rows would then name no run at all, one process further out than the
    defect this feature closes.
    """
    from delfin.agent import attribution

    got = attribution.from_env({"DELFIN_SWARM_RUN": "wave50",
                                "DELFIN_AGENT_LABEL": "wave50/t1-a1"})
    assert got.workflow_id == "wave50", (
        "the run a member works for is the unit whose cost is asked about")
    assert got.subagent_name == "t1-a1"
    assert got.molecule_id is None, "a member is told no molecule; never guess"


def test_an_ordinary_session_gains_nothing_from_the_environment():
    """No variables -> the same six nulls as before. Fail input: a default id."""
    from delfin.agent import attribution

    got = attribution.from_env({})
    assert got.is_empty()
    assert attribution.columns(got) == {k: None for k in COLUMNS}


def test_the_process_attribution_reaches_the_rows(home, monkeypatch):
    """``cli.main``'s adoption is what a member's rows carry.

    Driven through the real helper ``cli.main`` calls, then a real turn on top:
    the input that would make this fail is a helper that reads the environment
    but never installs it, leaving every later turn unattributed.
    """
    from delfin.agent import attribution, cli, repl

    monkeypatch.setenv("DELFIN_SWARM_RUN", "wave50")
    monkeypatch.setenv("DELFIN_AGENT_LABEL", "wave50/t1-a1")
    cli._process_attribution()
    try:
        engine = _FakeEngine("member")
        repl.run_turn(engine, "hi", sink=lambda item: None)
    finally:
        attribution.set_attribution()

    rows = tm.read_requests("member")
    assert rows, "no request row written"
    row = rows[-1]
    assert row["workflow_id"] == "wave50", (
        "the member's request row names no run, so a wave's spend cannot be "
        "totalled")
    assert row["subagent_name"] == "t1-a1"


def test_a_member_turn_run_on_the_repls_worker_thread_still_names_its_run(
        home, monkeypatch):
    """The turn runs on a NEW thread, and a ContextVar does not cross that.

    ``repl.py`` starts every turn on a ``threading.Thread`` named
    ``agent-turn`` (repl.py:1690). ``cli.main`` adopts the member's identity in
    the MAIN thread's context, so a scope that only sets a ContextVar is
    invisible where the turn actually runs: measured on the real product, a
    member process wrote ``workflow_id=None``.

    The input that would make this fail: an adoption that stores the record
    only in the main thread's context, with nothing on the engine's client —
    the thread then starts from an empty context and the row is unattributed.
    """
    from delfin.agent import attribution, cli, repl

    monkeypatch.setenv("DELFIN_SWARM_RUN", "wave50")
    monkeypatch.setenv("DELFIN_AGENT_LABEL", "wave50/t1-a2")
    cli._process_attribution()
    engine = _FakeEngine("worker")
    seen: dict = {}

    def _in_a_fresh_thread() -> None:
        # A brand-new thread: exactly the context the REPL's worker has.
        try:
            repl.run_turn(engine, "hi", sink=lambda item: None)
        except BaseException as exc:          # noqa: BLE001 - reported below
            seen["error"] = exc

    try:
        t = threading.Thread(target=_in_a_fresh_thread, name="agent-turn",
                             daemon=True)
        t.start()
        t.join(timeout=60)
    finally:
        attribution.set_attribution()

    assert not t.is_alive(), "the worker turn did not finish"
    assert "error" not in seen, f"the worker turn raised: {seen.get('error')}"
    rows = tm.read_requests("worker")
    assert rows, "no request row written from the worker thread"
    row = rows[-1]
    assert row["workflow_id"] == "wave50", (
        "the row written from the REPL's worker thread lost the run: the "
        "member's identity was installed only in the main thread's context")
    assert row["subagent_name"] == "t1-a2"


# ---------------------------------------------------------------------------
# 4. Never guess; identifiers only
# ---------------------------------------------------------------------------

def test_nothing_known_stays_null_rather_than_an_invented_id(home):
    """An unwrapped turn invents no id.

    ``agent_role`` is NOT null and must not be: the engine knows the role it is
    running as, and recording it is the point of the column. What must stay
    null is every id nobody supplied — a guessed ``molecule_id`` would bill one
    molecule for another's compute and nothing downstream could tell.
    """
    from delfin.agent import repl

    engine = _FakeEngine("bare")
    repl.run_turn(engine, "hi", sink=lambda item: None)
    row = tm.read_requests("bare")[-1]
    for col in COLUMNS:
        assert col in row, f"the row is missing the {col!r} column"
    for col in ("molecule_id", "workflow_id", "task_id", "subagent_name"):
        assert row[col] is None, f"{col} was guessed as {row[col]!r}"
    assert row["agent_role"] == "session_manager", (
        "the role the turn ran as is known and must be recorded")


@pytest.mark.parametrize("bad", [
    "/home/user/molecule.xyz",
    "M:\\scratch\\job1",
    "CC(=O)Oc1ccccc1C(=O)O",
    "alice@institute.example",
    "../../etc/passwd",
    "benzene dimer complex",
    "",
])
def test_a_path_or_smiles_is_never_stored_as_an_identifier(home, bad):
    """A value that is not an identifier becomes null, never a truncation.

    The input that would make this fail: a sanitiser that strips the offending
    characters instead of refusing the value — a cut path still names its first
    directories.
    """
    from delfin.agent import attribution

    with attribution.attribution_scope(molecule_id=bad):
        assert attribution.current() is None, (
            f"{bad!r} was accepted as an identifier")


def test_an_identifier_that_is_one_is_kept():
    """The gate must not reject what DELFIN actually uses."""
    from delfin.agent import attribution

    for good in ("mol_2026-12", "esd-chain", "t3a", "calc.4", "run:12"):
        with attribution.attribution_scope(molecule_id=good):
            got = attribution.current()
            assert got is not None and got.molecule_id == good


# ---------------------------------------------------------------------------
# 5. Control: a row written before the columns existed still reads
# ---------------------------------------------------------------------------

def test_an_old_row_still_reads(home):
    """Forward compatibility. This passes on the base too — it is the control."""
    p = tm.metrics_path("old")
    p.parent.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps({"ts": 1.0, "provider": "kit", "total_ms": 5}) + "\n")
    entries = tm.read("old")
    assert entries and entries[0]["provider"] == "kit"
    # An aggregate over a store mixing old and new rows must not care.
    tm.record("old", provider="kit", model="m", role="r", total_ms=9,
              ttft_ms=None, first_event_ms=None, output_chars=1,
              tool_calls=0, stopped=False)
    assert len(tm.read("old")) == 2
    stats = tm.aggregate_turn_stats.__wrapped__ if hasattr(
        tm.aggregate_turn_stats, "__wrapped__") else None
    assert stats is None or callable(stats)      # the API still exists

# ---------------------------------------------------------------------------
# 6. The process entry point is the caller that adopts the identity
# ---------------------------------------------------------------------------

def test_the_entry_point_adopts_the_members_identity(home, monkeypatch):
    """``main()``, not the helper it calls, is what a member's process runs.

    ``_process_attribution`` had only direct callers in the tests: every one
    of them called the private helper, so deleting the line in ``main()``
    that calls it left the whole suite green (measured: 38 passed with the
    call removed). The wiring is the part that decides whether a member's
    rows name their run at all, so it is asserted through ``main()`` itself.

    ``--version`` is the cheapest verb that reaches ``main()``'s adoption:
    argparse answers it and raises ``SystemExit(0)`` before any work that
    would need a terminal, a model or a network. That the identity is
    already adopted when that exit is raised is itself part of the claim:
    ``main()`` reads the environment BEFORE it parses anything, so a turn
    this process runs later cannot miss it.
    """
    from delfin.agent import attribution, cli

    monkeypatch.setenv("DELFIN_SWARM_RUN", "wave50")
    monkeypatch.setenv("DELFIN_AGENT_LABEL", "wave50/t1-a1")

    with pytest.raises(SystemExit) as exited:
        cli.main(["--version"])
    assert exited.value.code == 0

    adopted = attribution.current()
    assert adopted is not None, (
        "main() never adopted the member's identity: a process started by a "
        "swarm records six nulls for every call it makes")
    assert adopted.workflow_id == "wave50", (
        "the run a member works for is the unit whose cost is asked about")
    assert adopted.subagent_name == "t1-a1"


def test_an_ordinary_process_adopts_nothing(home, monkeypatch):
    """The control for the test above.

    Without it, a ``main()`` that adopted a constant would pass the positive
    case and attribute every interactive session to a run that is not its own.
    """
    from delfin.agent import attribution, cli

    monkeypatch.delenv("DELFIN_SWARM_RUN", raising=False)
    monkeypatch.delenv("DELFIN_AGENT_LABEL", raising=False)

    with pytest.raises(SystemExit):
        cli.main(["--version"])
    assert attribution.current() is None, (
        "an interactive session was attributed to a run it does not serve")
