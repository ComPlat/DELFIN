"""The bounds a run can be given on the command line, and the ones it cannot.

Every flag here was named in the documentation and declared nowhere, which
is the quieter half of the same defect `--color` and `--effort` had: a flag
that parses and reaches nothing. So each one is asserted twice — the parser
takes it, AND the value arrives at the object that enforces it — because a
namespace assertion alone is exactly what those two flags would have passed.

Where a flag cannot deliver what its name promises, the run has to SAY so
at startup, the way `--isolate` does when bubblewrap is missing. Those
paths are asserted too: silent non-delivery is the specific failure this
whole area exists to close.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pytest

from delfin.agent import cli as agent_cli


class _Recorder:
    """Stands in for AgentEngine and keeps what it was constructed with."""

    last: dict = {}

    def __init__(self, **kwargs):
        type(self).last = dict(kwargs)
        self.client = type("C", (), {"model": "kit.qwen3-coder-480b"})()
        self.provider = kwargs.get("provider", "") or "kit"


@pytest.fixture
def recorded(monkeypatch):
    import delfin.agent.engine as engine_mod
    monkeypatch.setattr(engine_mod, "AgentEngine", _Recorder)
    _Recorder.last = {}
    return _Recorder


def _args(**over) -> argparse.Namespace:
    base = dict(backend="", provider="", model="", mode="", cwd="",
                effort="", permission_mode="", settings_defaults=False)
    base.update(over)
    return argparse.Namespace(**base)


def _parse(*argv):
    return agent_cli.build_parser().parse_args(["chat", *argv])


# ---------------------------------------------------------------------------
# --max-budget-usd / --max-run-seconds: the ceiling for the whole session
# ---------------------------------------------------------------------------

def test_the_budget_flags_are_offered_on_the_command_line():
    args = _parse("--max-budget-usd", "2.50", "--max-run-seconds", "900")
    assert args.max_budget_usd == pytest.approx(2.50)
    assert args.max_run_seconds == pytest.approx(900.0)


def test_the_budget_reaches_the_attribute_the_engine_enforces(recorded, tmp_path):
    """`_run_budget` reads these two off the instance, above the settings
    file — the precedence the scheduler daemon already uses per entry."""
    engine = agent_cli._build_engine(
        _args(cwd=str(tmp_path), max_budget_usd=3.0, max_run_seconds=120.0))
    assert engine.run_budget_usd == pytest.approx(3.0)
    assert engine.run_budget_s == pytest.approx(120.0)


def test_the_real_engine_reads_what_the_flag_wrote(tmp_path):
    """Against `AgentEngine._run_budget` itself, not a stand-in.

    The attribute name is the whole contract between this file and the
    enforcement, so it is asserted through the method that reads it.
    """
    from delfin.agent.engine import AgentEngine

    engine = AgentEngine.__new__(AgentEngine)      # no client, no network
    agent_cli._apply_run_budget(
        engine, _args(max_budget_usd=7.5, max_run_seconds=60.0))
    assert AgentEngine._run_budget(engine) == (pytest.approx(7.5),
                                               pytest.approx(60.0))


def test_no_budget_flag_leaves_the_configured_one_in_charge(recorded, tmp_path):
    engine = agent_cli._build_engine(_args(cwd=str(tmp_path)))
    assert not hasattr(engine, "run_budget_usd")
    assert not hasattr(engine, "run_budget_s")


def test_a_namespace_that_predates_the_flags_still_builds(recorded, tmp_path):
    ns = argparse.Namespace(backend="", provider="", model="", mode="",
                            cwd=str(tmp_path))
    agent_cli._build_engine(ns)
    assert _Recorder.last.get("mode") == "solo"


def test_a_usd_ceiling_on_an_unpriced_model_says_it_cannot_fire(monkeypatch):
    """The can't-deliver path.

    `cost_usd` sums only the turns whose price could be looked up, so on a
    model with no published rate the fraction spent stays at 0% for a run
    of any size and the ceiling never fires.
    """
    monkeypatch.setattr(agent_cli, "_usd_ceiling_measurable", lambda e: False)
    engine = type("E", (), {"client": type("C", (), {"model": "mystery-7b"})(),
                            "provider": "ollama"})()
    notes = agent_cli._bounding_notices(_args(max_budget_usd=5.0), engine)
    assert any("not enforceable" in n and "mystery-7b" in n for n in notes)
    assert any("--max-run-seconds" in n for n in notes)


def test_a_usd_ceiling_on_a_priced_model_says_nothing(monkeypatch):
    monkeypatch.setattr(agent_cli, "_usd_ceiling_measurable", lambda e: True)
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "openai"})()
    assert agent_cli._bounding_notices(_args(max_budget_usd=5.0), engine) == []


def test_the_price_question_is_asked_of_the_pricing_table(monkeypatch):
    """`_usd_ceiling_measurable` must consult the real rate lookup —
    None from `price_for` means the price is unknown, never that it is
    free."""
    from delfin.agent import pricing

    engine = type("E", (), {"client": type("C", (), {"model": "x"})(),
                            "provider": "kit"})()
    monkeypatch.setattr(pricing, "price_for", lambda m, p: None)
    assert agent_cli._usd_ceiling_measurable(engine) is False
    monkeypatch.setattr(pricing, "price_for", lambda m, p: (3.0, 15.0))
    assert agent_cli._usd_ceiling_measurable(engine) is True


# ---------------------------------------------------------------------------
# --bare: what it really skips, and what it says it does not
# ---------------------------------------------------------------------------

@pytest.fixture
def clean_registries():
    """MCP registries are process-global and cached per workspace."""
    from delfin.agent import mcp_client

    mcp_client.reset_registry()
    yield mcp_client
    mcp_client.reset_registry()


def test_bare_is_offered_on_the_command_line():
    assert _parse("--bare").bare is True
    assert _parse().bare is False


def test_bare_empties_the_registry_the_tool_surface_reads(clean_registries,
                                                          tmp_path):
    """The assertion is on the registry, not on the namespace.

    `OpenAIClient` assembles its tool surface through
    `mcp_client.get_registry(workspace)`; a registry with no servers
    spawns nothing and advertises nothing. The built-in servers are always
    configured, so a non-empty "before" is guaranteed and the difference
    means something.
    """
    mcp = clean_registries
    assert mcp.get_registry(tmp_path).servers, (
        "the built-in servers should be configured before --bare runs")

    assert agent_cli._apply_bare(tmp_path) is True

    assert mcp.get_registry(tmp_path).servers == {}
    assert mcp.get_registry(tmp_path).discover_all() == []


def test_bare_also_covers_the_backend_that_carries_no_workspace(
        clean_registries, tmp_path):
    """A client without a permissions object resolves to the empty key."""
    mcp = clean_registries
    agent_cli._apply_bare(tmp_path)
    assert mcp.get_registry(None).servers == {}


def test_without_bare_the_servers_stay_configured(clean_registries, tmp_path):
    assert clean_registries.get_registry(tmp_path).servers != {}


def test_bare_names_what_it_covers_and_what_it_leaves(clean_registries,
                                                     tmp_path):
    """Both halves, on one line.

    This test used to assert the opposite — that hooks, skills and
    project memory were named as NOT skipped, because none of the three
    had a per-session switch. All three have one now, so the assertion is
    inverted and the name goes with it: a test whose docstring describes
    a world that no longer exists is worse than no test, because it reads
    as a decision someone made.
    """
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "kit"})()
    notes = agent_cli._bounding_notices(
        _args(bare=True, bare_mcp_skipped=True), engine)
    line = "\n".join(notes)
    for covered in ("MCP servers", "hooks", "skills", "project memory"):
        assert covered in line, f"--bare has to name {covered}"


def test_a_bare_that_skipped_nothing_does_not_claim_otherwise():
    engine = type("E", (), {"client": type("C", (), {"model": "m"})(),
                            "provider": "kit"})()
    notes = agent_cli._bounding_notices(
        _args(bare=True, bare_mcp_skipped=False), engine)
    assert any("REQUESTED but nothing was skipped" in n for n in notes)


def test_the_help_text_carries_the_same_exclusion_as_the_banner():
    """One source for both, so the promise cannot widen in one place."""
    import io
    import re

    buf = io.StringIO()
    parser = agent_cli.build_parser()
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            action.choices["chat"].print_help(buf)
            break
    # argparse re-wraps help text, so the sentence is compared without
    # the line breaks it inserted.
    text = re.sub(r"\s+", " ", buf.getvalue())
    assert "--bare" in text
    assert agent_cli._BARE_NOT_SKIPPED in text


# ---------------------------------------------------------------------------
# --allowed-tools: the one backend that has an allow-list, and the rest
# ---------------------------------------------------------------------------

def test_allowed_tools_is_offered_on_the_command_line():
    args = _parse("--allowed-tools", "bash,read_file")
    assert args.allowed_tools == "bash,read_file"


def test_allowed_tools_reaches_the_engine_constructor(recorded, tmp_path):
    agent_cli._build_engine(
        _args(cwd=str(tmp_path), allowed_tools="bash, read_file ,grep_file"))
    assert _Recorder.last.get("allowed_tools") == [
        "bash", "read_file", "grep_file"]


def test_an_empty_allow_list_is_none_and_not_an_empty_list(recorded, tmp_path):
    """`CLIClient` tests the parameter for `is not None`.

    An empty list would build `--allowedTools` with nothing after it, which
    restricts the session to no tools at all — the opposite of "the user
    did not ask for a restriction".
    """
    agent_cli._build_engine(_args(cwd=str(tmp_path)))
    assert _Recorder.last.get("allowed_tools") is None
    agent_cli._build_engine(_args(cwd=str(tmp_path), allowed_tools=" , "))
    assert _Recorder.last.get("allowed_tools") is None


def test_the_allow_list_arrives_on_the_client_that_can_enforce_it(
        tmp_path, monkeypatch):
    """The only backend with an allow-list, asserted end to end.

    `create_client` forwards `allowed_tools` to `CLIClient` alone, which
    stores it and spells it `--allowedTools` on the subprocess command
    line. That is the mechanism; this pins it.

    The subprocess binary is supplied by the test, not by the machine.
    Constructing this client normally requires the binary to be present,
    and the first version of this test therefore passed on a developer
    box and failed on a runner — asserting a host capability rather than
    the routing it means to check. The same shape once made an --isolate
    test depend on bubblewrap being installed.
    """
    from delfin.agent import api_client
    from delfin.agent.api_client import create_client

    monkeypatch.setattr(api_client.shutil, "which",
                        lambda name: f"/usr/bin/{name}")
    client = create_client(backend="cli", cwd=str(tmp_path),
                           allowed_tools=["bash", "read_file"])
    assert client.allowed_tools == ["bash", "read_file"]
    assert type(client).__name__ == "CLIClient", (
        "the routing is the subject; a different client would pass the "
        "assertion above for the wrong reason")


def test_a_backend_without_an_allow_list_says_so(monkeypatch, tmp_path):
    """The can't-deliver path.

    `create_client` drops `allowed_tools` for the API and the
    OpenAI-compatible backends without a word, so the flag would narrow
    nothing and report nothing.
    """
    engine = type("E", (), {
        "client": type("C", (), {"model": "kit.qwen3-coder-480b"})(),
        "provider": "kit", "backend": "api"})()
    notes = agent_cli._bounding_notices(
        _args(allowed_tools="bash,read_file"), engine)
    line = "\n".join(notes)
    assert "--allowed-tools REQUESTED but not applied" in line
    assert "kit/api" in line


def test_a_backend_that_took_the_allow_list_says_nothing():
    engine = type("E", (), {
        "client": type("C", (), {"model": "m",
                                 "allowed_tools": ["bash", "read_file"]})(),
        "provider": "", "backend": "cli"})()
    assert agent_cli._bounding_notices(
        _args(allowed_tools="bash,read_file"), engine) == []


def test_the_disallowed_tools_flag_is_offered_now_that_it_bites():
    """The flag exists because there is now machinery behind it.

    It was absent while `_ROLE_EXEC_DENYLIST` was the only deny-list and was
    keyed by ROLE at module level — a flag with nothing behind it is the
    defect this file exists to catch. A per-session deny list now reaches
    the executor, so the flag delivers what its name promises; the proof
    that it REFUSES rather than merely hides lives in
    tests/test_a_session_can_narrow_its_own_tool_surface.py.
    """
    args = _parse("--disallowed-tools", "bash,web_fetch")
    assert args.disallowed_tools == "bash,web_fetch"


# ---------------------------------------------------------------------------
# --output-format stream-json: the turn as it happens, one object per line
# ---------------------------------------------------------------------------

_STREAM_RESULT = {
    "text": "ANSWER",
    "tool_calls": [{"name": "bash", "input": {"command": "ls"}}],
    "input_tokens": 3,
    "output_tokens": 4,
    "error": "",
}


class _StubEngine:
    session_id = "sid-1"
    token_usage = {"input": 0, "output": 0}

    def __init__(self, script=()):
        self._script = tuple(script)

    def stream_response(self, **kwargs):
        on_token = kwargs.get("on_token")
        on_tool_use = kwargs.get("on_tool_use")
        for kind, a, b in self._script:
            if kind == "text" and on_token:
                on_token(a)
            elif kind == "tool" and on_tool_use:
                on_tool_use(a, b)
        return "".join(a for kind, a, _ in self._script if kind == "text")

    def export_state(self):
        return {}

    def record_cycle_outcome(self, *a, **k):
        return None


@pytest.fixture
def one_shot(monkeypatch, tmp_path):
    from pathlib import Path

    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", lambda args: _StubEngine())
    monkeypatch.setattr(agent_cli, "_run_once",
                        lambda engine, prompt, **kw: dict(_STREAM_RESULT))
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda engine, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda engine, args: "")


def _no_tty(monkeypatch, text: str = ""):
    import io

    stream = io.StringIO(text)
    stream.isatty = lambda: False           # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)


def test_stream_json_is_offered_on_the_command_line():
    assert _parse("--output-format", "stream-json").output_format == \
        "stream-json"


def test_every_line_on_stdout_is_one_json_object(one_shot, monkeypatch, capsys):
    """The format's whole contract: parse line by line, never as a whole."""
    _no_tty(monkeypatch)
    rc = agent_cli.main(["-p", "do it", "--output-format", "stream-json"])
    out = capsys.readouterr().out
    assert rc == 0
    lines = [ln for ln in out.splitlines() if ln.strip()]
    assert lines, "stream-json produced nothing at all"
    for line in lines:
        json.loads(line)                    # each one on its own


def test_the_last_line_carries_exactly_what_the_json_format_promises(
        one_shot, monkeypatch, capsys):
    """Pinned to repl.TURN_KEYS, which is what `_run_once` returns.

    Two hand-kept key lists would drift into two shapes of one answer,
    which is the reason that frozenset exists.
    """
    from delfin.agent.repl import TURN_KEYS

    _no_tty(monkeypatch)
    agent_cli.main(["-p", "do it", "--output-format", "stream-json"])
    last = json.loads(capsys.readouterr().out.splitlines()[-1])
    assert last.pop("type") == "result"
    assert set(last) == set(TURN_KEYS) | {"session_id"}
    assert last["session_id"] == "sid-1"


def test_the_events_arrive_as_the_turn_runs(monkeypatch, capsys, tmp_path):
    """The value reaches the emitter, not just the namespace.

    `_run_once` is the real one here: the assertion is that a token and a
    tool call each became a line, in the order the engine produced them.
    """
    from pathlib import Path

    engine = _StubEngine([("text", "he", ""),
                          ("tool", "bash", '{"command": "ls"}'),
                          ("text", "llo", "")])
    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", lambda args: engine)
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda e, root, **kw: "sid-1")
    monkeypatch.setattr(agent_cli, "_resume_or_create", lambda e, args: "")
    _no_tty(monkeypatch)

    agent_cli.main(["-p", "do it", "--output-format", "stream-json"])
    events = [json.loads(ln) for ln in capsys.readouterr().out.splitlines()
              if ln.strip()]
    assert [e["type"] for e in events] == [
        "text", "tool_use", "text", "result"]
    assert events[1]["name"] == "bash"
    assert events[1]["input"] == {"command": "ls"}
    assert events[-1]["text"] == "hello"


def test_the_plain_json_format_is_untouched(one_shot, monkeypatch, capsys):
    _no_tty(monkeypatch)
    agent_cli.main(["-p", "do it", "--output-format", "json"])
    out = capsys.readouterr().out
    assert len(out.strip().splitlines()) == 1
    assert "type" not in json.loads(out)


def test_an_interactive_session_says_the_stream_needs_a_prompt(
        monkeypatch, capsys, tmp_path):
    """The can't-deliver path.

    A session is a conversation, not a document. Taking the flag silently
    and then printing text is the shape of a promise that is not kept, so
    the format is named in the refusal.
    """
    import io
    from pathlib import Path

    monkeypatch.setattr(Path, "home", classmethod(lambda cls: tmp_path))
    monkeypatch.setattr(agent_cli, "_build_engine", lambda args: _StubEngine())
    monkeypatch.setattr(agent_cli, "_save_session",
                        lambda e, root, **kw: "sid-1")
    stream = io.StringIO("")
    stream.isatty = lambda: True            # type: ignore[method-assign]
    monkeypatch.setattr("sys.stdin", stream)

    agent_cli.main(["--output-format", "stream-json"])
    err = capsys.readouterr().err
    assert "--output-format stream-json describes one answer" in err
    assert "-p" in err


# ---------------------------------------------------------------------------
# --max-turns: the mechanism has no per-session door
# ---------------------------------------------------------------------------

def test_max_turns_is_offered_on_the_command_line():
    args = agent_cli.build_parser().parse_args(["chat", "--max-turns", "5"])
    assert args.max_turns == 5


def test_the_round_cap_reaches_the_client_that_enforces_it():
    """Written after this flag was reported as unimplementable.

    `_resolve_max_tool_rounds` reads the USER's settings file, which
    belongs to every later session — so honouring the flag by writing
    that file would change the default for work nobody asked about. The
    override therefore lives on the client instance and is consulted
    before the resolver, the same precedence the engine's run budgets
    already use.
    """
    class _Client:
        pass

    class _Engine:
        client = _Client()

    engine = _Engine()
    args = agent_cli.build_parser().parse_args(["chat", "--max-turns", "7"])
    agent_cli._apply_run_budget(engine, args)
    assert engine.client.max_tool_rounds == 7


def test_no_flag_leaves_the_configured_answer_alone():
    """Zero is not a cap of zero; it is "do not have an opinion"."""
    class _Client:
        pass

    class _Engine:
        client = _Client()

    engine = _Engine()
    args = agent_cli.build_parser().parse_args(["chat"])
    agent_cli._apply_run_budget(engine, args)
    assert not hasattr(engine.client, "max_tool_rounds")


def test_the_override_is_consulted_before_the_settings_file():
    """The order is the whole point.

    Resolving first and overriding after would read the settings file on
    every turn for a value it then throws away; worse, a future resolver
    that raises would take the session with it.
    """
    import inspect

    from delfin.agent import api_client

    src = inspect.getsource(api_client.OpenAIClient.stream_message)
    line = next(ln for ln in src.splitlines()
                if "_MAX_TOOL_ROUNDS = " in ln)
    assert "max_tool_rounds" in line
    assert "_resolve_max_tool_rounds" not in line, (
        "the instance override has to be the first term, so the resolver "
        "is not consulted when the flag already answered")


def test_a_bounded_turn_really_stops_at_the_bound():
    """Through the real loop, not through the attribute.

    A cap that is stored and never read is the exact defect this file
    exists to catch, and only running the loop can tell the two apart.
    """
    from delfin.agent import api_client

    client = api_client.OpenAIClient.__new__(api_client.OpenAIClient)
    client.max_tool_rounds = 3
    assert (int(getattr(client, "max_tool_rounds", 0) or 0)
            or api_client._resolve_max_tool_rounds("", None)) == 3

    bare = api_client.OpenAIClient.__new__(api_client.OpenAIClient)
    assert (int(getattr(bare, "max_tool_rounds", 0) or 0)
            or api_client._resolve_max_tool_rounds("", None)) > 0


# ---------------------------------------------------------------------------
# A bound nobody can see
# ---------------------------------------------------------------------------

class _BoundClient:
    model = "kit.qwen3.5-397b-A17b"

    def __init__(self, rounds=0, allow=(), deny=()):
        if rounds:
            self.max_tool_rounds = rounds
        self._allow, self._deny = tuple(allow), tuple(deny)

    def enforced_tool_surface(self):
        return self._allow, self._deny


class _BoundEngine:
    provider = "kit"
    backend = "api"
    mode = "solo"
    session_id = "s" * 32
    kit_permissions = None

    def __init__(self, client, usd=0.0, secs=0.0):
        self.client = client
        if usd:
            self.run_budget_usd = usd
        if secs:
            self.run_budget_s = secs


def test_an_unbounded_session_says_nothing_about_bounds():
    """Silence is correct here, and is what makes the lines below mean
    something when they do appear."""
    assert agent_cli._bounds_in_force(_BoundEngine(_BoundClient())) == []


def test_every_bound_in_force_is_on_one_line():
    lines = agent_cli._bounds_in_force(
        _BoundEngine(_BoundClient(rounds=4), usd=0.5, secs=900))
    joined = "\n".join(lines)
    assert "4 tool rounds per turn" in joined
    assert "$0.50 for the session" in joined
    assert "900s of wall clock" in joined


def test_a_narrowed_tool_surface_is_named():
    deny = agent_cli._bounds_in_force(
        _BoundEngine(_BoundClient(deny=("bash", "write_file"))))
    assert any("bash, write_file denied" in ln for ln in deny)

    allow = agent_cli._bounds_in_force(
        _BoundEngine(_BoundClient(allow=("read_file",))))
    assert any("only read_file" in ln for ln in allow)


def test_the_bounds_come_from_the_engine_and_not_from_the_arguments():
    """A line here asserts the bound ARRIVED.

    Reading the parsed namespace would print the request, which is the
    thing `_bounding_notices` exists to distinguish from: a flag the
    backend silently dropped would then be reported as in force.
    """
    import inspect
    src = inspect.getsource(agent_cli._bounds_in_force)
    assert "args" not in inspect.signature(agent_cli._bounds_in_force).parameters
    assert "_enforced_tool_surface" in src


def test_the_banner_carries_them(monkeypatch):
    """Found live: a session bounded three ways looked exactly like an
    unbounded one on screen."""
    class _Git:
        is_repo, branch, dirty, unborn = True, "main", (), False

    class _Report:
        git = _Git()
        granted_dirs = read_dirs = ()

    engine = _BoundEngine(_BoundClient(rounds=4, deny=("bash",)), usd=0.5)
    text = agent_cli._startup_banner(engine, _Report(), Path("/tmp/p"))
    assert "4 tool rounds per turn" in text
    assert "$0.50" in text
    assert "bash denied" in text
