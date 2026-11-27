"""Every entry path for a tool call passes the one entry where the gate decides.

A tool call can arrive from the streaming turn, from a Code Mode program over
its RPC socket, from an external CLI program's PreToolUse hook, or from a
person pressing a button outside any turn. This file drives one call down each
of those routes and asks one question of each: did the route consult the gate
DELFIN already has — ``_DocToolExecutor.execute`` for a native name,
``_mcp_call_refusal`` + ``_gate_mcp_tool`` for a namespaced one — before
anything ran?

The measurement is a refusal only the gate can produce, not a substring of
prose some check happens to print: the call's own tool name is put in the
ADMINISTRATOR's managed policy (``permissions.deny``, injected at a path —
never /etc). Both entries consult the rules, and only they do, and the refusal
text comes from ``managed_policy.deny_text``. The managed layer is the right
sentinel because it is the one that holds in every mode, with nobody to ask,
and across the user/session layer: a route that skips the entry cannot print
it, and a route that reaches the entry cannot skip it. Nothing can execute
before it answers either — a rule is matched against the name and arguments
before dispatch.

Routes (``ROUTES``), each a product entry point, not a helper:

``stream``  the streaming turn: a scripted model emits one tool call and
            ``OpenAIClient.stream_message`` dispatches it — the route a
            model's call actually takes, with the fake stream the MCP
            meta-tool tests already use;
``code``    a Code Mode program: ``code_mode.RpcServer.handle_payload``, the
            payload a caged child sends over its socket;
``cli``     the external CLI backend's PreToolUse hook: ``cli_gate.main`` fed
            the event JSON on stdin, the decision read off its stdout — the
            text that program obeys (same harness as
            ``test_the_cli_hook_decides_every_tool``);
``merge``   the session-worktree merge the dashboard's session list and the
            service's ``session_merge`` drive (``session_worktree.merge``);
``orphan``  the orchestration's merge of the worktrees its delegates kept
            (``_merge_orchestration_worktrees``, the ``orchestrate`` tool's
            completion step);
``mcpdisp`` ``execute_mcp_call``, the one MCP dispatch, with a fake registry.

Tools: ``worktree_merge`` (native, and it really writes), the two MCP
resource/prompt meta-tools the loop adds, and ``mcp__`` names.

Two measured facts per route, both falsifiable:

* with the name denied by policy the route ENDS THE CALL BEFORE IT, and for
  the merge routes the target file is provably unchanged; and
* with the policy empty the route RUNS the call, or refuses for the one reason
  ``DOCUMENTED_REFUSAL`` names. A route that refused everything would pass the
  first half for free, so the second half is what makes it mean a gate.

RED on the base commit, for the reason each row names: ``merge`` and
``orphan`` call ``_execute_worktree_merge`` on a bare executor below both
entries, so no rule is consulted and the merge applies the worktree's change;
the stream's two MCP meta-tools run inline with only a role check, so the
rules, the session lists and the mode are never consulted for them.
``cli`` and ``mcpdisp`` are the controls — both have always gone through an
entry, and their rows pass on both trees.

No model, no network, no ORCA, no SLURM. The registry is a fake; the stream is
a scripted list of chunks; no shell command of any grid is run.
"""

from __future__ import annotations

import ast
import contextlib
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO))

from delfin import user_settings as US  # noqa: E402
from delfin.agent import api_client as A  # noqa: E402
from delfin.agent import cli_gate  # noqa: E402
from delfin.agent import code_mode as CM  # noqa: E402
from delfin.agent import kit_settings  # noqa: E402
from delfin.agent import managed_policy as MP  # noqa: E402
from delfin.agent import mcp_client as MC  # noqa: E402
from delfin.agent import model_capabilities as MCAP  # noqa: E402
from delfin.agent import security_events  # noqa: E402
from delfin.agent import session_worktree as SW  # noqa: E402
from delfin.agent import tool_catalog  # noqa: E402
from delfin.agent.api_client import (  # noqa: E402
    KitToolPermissions, _DocToolExecutor,
)

#: The MCP resource/prompt meta-tools: no schema in the catalogue, added to
#: the advertised surface by the loop when a server exposes resources or
#: prompts (api_client, ``is_mcp_meta``). They are the one family that is
#: neither a catalogue tool nor an MCP-server tool, and the loop dispatches
#: them in a third branch.
META_TOOLS = frozenset({"mcp_read_resource", "mcp_get_prompt"})

#: A namespaced call for a server DELFIN does not ship. It routes to the MCP
#: entry, which is deny-by-default for anything off its read-only list.
FOREIGN_MCP = "mcp__lab__write_note"

#: A namespaced call the Code Mode allow-list admits, so a program can reach
#: the MCP entry at all.
PIPELINE_MCP = "mcp__delfin-tools__get_guide"

#: The refusal only the rule layer produces (``managed_policy.deny_text``).
RULE_MARK = "refused by the managed policy rule"


def registered_tools() -> set[str]:
    """Every tool name DELFIN registers for a model to call."""
    names = {str(t.get("function", {}).get("name") or "")
             for t in tool_catalog._DOC_TOOLS_OPENAI}
    names.discard("")
    return names | META_TOOLS


def bare(name: str) -> str:
    """A namespaced ``mcp__server__tool`` reduced to its tool name."""
    return name.rsplit("__", 1)[-1] if name.startswith("mcp__") else name


# ---------------------------------------------------------------------------
# Routes, and which tool each one can actually reach
# ---------------------------------------------------------------------------

ROUTES = {
    "stream": "the streaming turn (OpenAIClient.stream_message)",
    "code": "a Code Mode program (code_mode.RpcServer payload)",
    "cli": "the external CLI backend's PreToolUse hook (cli_gate.main)",
    "merge": "the dashboard/service session-worktree merge",
    "orphan": "the orchestration's merge of its delegates' worktrees",
    "mcpdisp": "the ONE MCP dispatch (execute_mcp_call)",
}

#: One entry per measured (route, tool). A route is measured only for a call it
#: can reach: the merge routes exist for one verb; a program may name only the
#: RPC server's read-only allow-list, so it is measured with an MCP name that
#: list admits; the hook sees every name the CLI program asks about.
MEASURED = (
    ("stream", "worktree_merge"),
    ("stream", "mcp_read_resource"),
    ("stream", "mcp_get_prompt"),
    ("code", PIPELINE_MCP),
    ("cli", FOREIGN_MCP),
    ("merge", "worktree_merge"),
    ("orphan", "worktree_merge"),
    ("mcpdisp", FOREIGN_MCP),
)

#: A route that does not run an un-denied call must give the reason it is
#: allowed to (asserted present). Without this a route could pass the gate row
#: by refusing at random or by dying quietly, which is not a gate. Empty today:
#: every measured route runs its call when nothing refuses it, which the
#: control row asserts. The mechanism is here so a route added later — one
#: that legitimately refuses, like the RPC server's read-only allow-list over a
#: name off it — states its reason instead of being excused.
DOCUMENTED_REFUSAL: dict[tuple[str, str], str] = {}


# ---------------------------------------------------------------------------
# The world: a git workspace, an injected managed policy, no ambient settings
# ---------------------------------------------------------------------------


@pytest.fixture
def env(tmp_path, monkeypatch):
    """A git workspace + a private HOME with no settings + an empty policy.

    The worktree the merge routes merge lives under
    ``<repo>/.delfin/worktrees`` — the layout ``session_worktree.create`` makes
    — so it sits inside the session's own root: the routes' containment check
    passes, and the only thing that can stop the merge is the gate.
    """
    home = tmp_path / "home"
    home.mkdir()
    monkeypatch.setenv("HOME", str(home))
    monkeypatch.setattr(
        US, "get_settings_path",
        lambda base_path=None: (Path(base_path).expanduser() if base_path
                                else home / US.SETTINGS_FILE_NAME))
    monkeypatch.setattr(kit_settings, "USER_SETTINGS_PATH",
                        home / ".delfin" / "settings.json")
    monkeypatch.delenv(MP.ENV_VAR, raising=False)
    etc = tmp_path / "etc" / "delfin"
    etc.mkdir(parents=True)
    etc.chmod(0o755)
    etc.parent.chmod(0o755)
    policy = etc / "managed_policy.json"
    monkeypatch.setattr(MP, "_PATHS_OVERRIDE", (str(policy),))
    monkeypatch.setattr(MP, "_ADMIN_UIDS_OVERRIDE", frozenset({0, os.getuid()}))
    monkeypatch.setattr(MP, "_ISOLATION_PROBE", lambda backend="": True)

    def write_policy(deny):
        """Put NAMES in the administrator's deny list (or empty it)."""
        policy.write_text(json.dumps({"permissions": {"deny": list(deny)}}),
                          encoding="utf-8")
        policy.chmod(0o644)
        MP.clear_cache()

    write_policy(())
    monkeypatch.setenv("DELFIN_CLASSIFIER", "off")
    monkeypatch.setattr(cli_gate, "HOOK_WAIT_S", 30.0)

    ws = tmp_path / "repo"
    (ws / "src").mkdir(parents=True)
    (ws / "src" / "a.py").write_text("A\n", encoding="utf-8")
    _git(ws, "init", "-q")
    _git(ws, "-c", "user.email=t@t", "-c", "user.name=t", "add", "-A")
    _git(ws, "-c", "user.email=t@t", "-c", "user.name=t", "commit", "-qm",
         "base")
    security_events.clear()
    yield {"ws": ws.resolve(), "tmp": tmp_path, "home": home,
           "policy": policy, "write_policy": write_policy}
    MP.clear_cache()
    security_events.clear()


def _git(cwd, *args):
    subprocess.run(["git", *args], cwd=str(cwd), check=True,
                   capture_output=True, text=True)


def _worktree(env):
    """A worktree of env['ws'] with one edited file (the merge's payload)."""
    ws = Path(env["ws"])
    root = ws / ".delfin" / "worktrees"
    n = len(list(root.glob("*"))) if root.is_dir() else 0
    wt = root / f"w{n}"
    _git(ws, "worktree", "add", "-q", "-b", f"session/probe{n}", str(wt))
    (wt / "src" / "a.py").write_text("CHANGED BY THE SESSION\n",
                                     encoding="utf-8")
    return wt


def _perms(env, *, mode="bypassPermissions", extra=()):
    """Permissions with nobody to ask.

    ``bypassPermissions`` is deliberate: its promise is that it asks nothing,
    so any refusal a route still makes is the gate's, not a dialog nobody
    answered.
    """
    p = KitToolPermissions(workspace=Path(env["ws"]), mode=mode,
                           extra_workspace_dirs=tuple(Path(e) for e in extra))
    p.confirm_callback = None
    p.skip_hook_discovery = True
    p.skip_skill_discovery = True
    return p


def _args(tool, env, wt):
    if tool == "worktree_merge":
        return {"path": str(wt), "target_dir": str(env["ws"])}
    if tool == "mcp_read_resource":
        return {"server": "remote", "uri": "file:///notes.md"}
    if tool == "mcp_get_prompt":
        return {"name": "mcp__remote__summarize", "arguments": {}}
    return {"path": "note.txt", "text": "hello"}


def _refused_by_rule(result) -> bool:
    return RULE_MARK in str(result)


def _merged(env) -> bool:
    """Did a worktree's change reach the target file?"""
    return (Path(env["ws"]) / "src" / "a.py").read_text(
        encoding="utf-8") != "A\n"


# ---------------------------------------------------------------------------
# A fake MCP registry that answers, so "it ran" is observable
# ---------------------------------------------------------------------------


class _Registry:
    def __init__(self):
        self.calls: list = []

    def discover_all(self):
        return []

    def discover_resources(self):
        return [MC.MCPResource(server="remote", uri="file:///notes.md",
                               name="Notes")]

    def discover_prompts(self):
        return [MC.MCPPrompt(server="remote", name="summarize",
                             description="Summarize")]

    def call(self, name, arguments, agent_session_id=""):
        self.calls.append(("call", name, dict(arguments)))
        return json.dumps({"ok": True, "tool": name})

    def read_resource(self, server, uri):
        self.calls.append(("read_resource", server, uri))
        return "RESOURCE-BODY"

    def get_prompt(self, name, arguments=None):
        self.calls.append(("get_prompt", name, dict(arguments or {})))
        return "PROMPT-BODY"


# ---------------------------------------------------------------------------
# The scripted model stream (the shape the MCP meta-tool tests use)
# ---------------------------------------------------------------------------


class _Fn:
    def __init__(self, name, arguments):
        self.name, self.arguments = name, arguments


class _TC:
    def __init__(self, name, arguments):
        self.index, self.id = 0, "c1"
        self.function = _Fn(name, arguments)


class _Delta:
    def __init__(self, tool_calls=None):
        self.content, self.tool_calls = None, tool_calls


class _Choice:
    def __init__(self, delta, finish=None):
        self.delta, self.finish_reason = delta, finish


class _Usage:
    prompt_tokens, completion_tokens = 5, 3


class _Chunk:
    def __init__(self, choices, usage=None):
        self.choices, self.usage = choices, usage


class _Stream:
    def __init__(self, chunks):
        self._chunks = chunks

    def __iter__(self):
        return iter(self._chunks)

    def close(self):
        pass


def _tool_call_stream(name, args_json):
    return _Stream([
        _Chunk([_Choice(_Delta([_TC(name, args_json)]))]),
        _Chunk([_Choice(_Delta(), finish="tool_calls")], usage=_Usage()),
    ])


def _text_stream(text="done"):
    class _D:
        content, tool_calls = text, None

    class _C:
        def __init__(self):
            self.delta, self.finish_reason = _D(), "stop"

    class _Ch:
        def __init__(self):
            self.choices, self.usage = [_C()], _Usage()

    return _Stream([_Ch()])


# ---------------------------------------------------------------------------
# The six routes, each driving ONE call
# ---------------------------------------------------------------------------


def _route_stream(tool, args, env, perms, reg, monkeypatch):
    """One model tool call through the real streaming loop."""
    caps = MCAP.ModelCapabilities(model="m", provider="ollama",
                                 context_window=200_000, supports_tools=True)
    monkeypatch.setattr(MCAP, "resolve", lambda *a, **k: caps)
    monkeypatch.setattr(MC, "get_registry", lambda *a, **k: reg)
    client = A.create_client(backend="api", provider="ollama",
                             model="qwen2.5-coder:7b", cwd=str(env["ws"]))
    client._permissions = perms
    streams = [_tool_call_stream(tool, json.dumps(args)), _text_stream()]
    n = {"i": 0}

    def _create(**kwargs):
        s = streams[min(n["i"], len(streams) - 1)]
        n["i"] += 1
        return s

    client.client.chat.completions.create = _create
    got: list[str] = []
    for ev in client.stream_message("sys", [{"role": "user", "content": "go"}],
                                    max_tokens=50):
        if getattr(ev, "type", "") == "tool_result":
            got.append(str(getattr(ev, "tool_output", "") or ""))
    return "\n".join(got)


def _route_code(tool, args, env, perms, reg, monkeypatch):
    """A program's request, as the caged child sends it to the RPC server."""
    monkeypatch.setattr(MC, "get_registry", lambda *a, **k: reg)
    server = CM.RpcServer(_DocToolExecutor(), permissions=perms)
    return server.handle_payload(json.dumps({"tool": tool,
                                             "arguments": dict(args)}))


def _route_cli(tool, args, env, perms, reg, monkeypatch):
    """The external CLI program's hook: its event in, its decision out."""
    monkeypatch.setattr(MC, "get_registry", lambda *a, **k: reg)
    token = cli_gate.new_token()
    cli_gate.write_state(token, workspace=str(env["ws"]),
                         mode=getattr(perms, "mode", "default"),
                         can_ask=False,
                         extra_dirs=tuple(str(p) for p in
                                          getattr(perms,
                                                  "extra_workspace_dirs", ())))
    event = {"hook_event_name": "PreToolUse", "tool_name": tool,
             "tool_input": dict(args), "cwd": str(env["ws"]),
             "permission_mode": getattr(perms, "mode", "default")}
    out = io.StringIO()
    old = sys.stdin
    sys.stdin = io.StringIO(json.dumps(event))
    try:
        with contextlib.redirect_stdout(out):
            assert cli_gate.main(["cli_gate", token]) == 0
    finally:
        sys.stdin = old
    printed = out.getvalue().strip()
    if not printed:
        return json.dumps({"decision": None, "reason": ""})
    spec = json.loads(printed)["hookSpecificOutput"]
    return json.dumps({"decision": spec["permissionDecision"],
                       "reason": spec.get("permissionDecisionReason", "")})


def _route_merge(tool, args, env, perms, reg, monkeypatch):
    """The dashboard / service route: session_worktree.merge."""
    return json.dumps(SW.merge(str(args["path"])))


def _route_orphan(tool, args, env, perms, reg, monkeypatch):
    """The orchestration's merge step, reached as a tool's completion work."""
    out = {"stages": {"s1": [{"description": "d",
                              "worktree": {"final_path": str(args["path"]),
                                           "branch": "session/probe"}}]}}
    return json.dumps(_DocToolExecutor()._merge_orchestration_worktrees(
        out, perms))


def _route_mcpdisp(tool, args, env, perms, reg, monkeypatch):
    """The ONE MCP dispatch, with a fake registry."""
    monkeypatch.setattr(MC, "get_registry", lambda *a, **k: reg)
    return A.execute_mcp_call(tool, dict(args), perms, _DocToolExecutor())


ROUTE_FN = {
    "stream": _route_stream,
    "code": _route_code,
    "cli": _route_cli,
    "merge": _route_merge,
    "orphan": _route_orphan,
    "mcpdisp": _route_mcpdisp,
}


def _drive(route, tool, env, monkeypatch, *, denied=True,
           mode="bypassPermissions"):
    """Drive one call down one route -> (result text, registry).

    ``denied`` puts the call's own tool name in the administrator's deny list.
    A fresh worktree per merge call, so the refused run and the live run each
    have a change to merge (or not).
    """
    env["write_policy"]((bare(tool),) if denied else ())
    wt = _worktree(env) if tool == "worktree_merge" else None
    perms = _perms(env, mode=mode,
                   extra=[str(wt)] if wt is not None else ())
    reg = _Registry()
    monkeypatch.setattr(MC, "get_registry", lambda *a, **k: reg)
    out = ROUTE_FN[route](tool, _args(tool, env, wt), env, perms, reg,
                          monkeypatch)
    return out, reg


def _reached(out, reg, route, tool, env):
    """The observable proof that this call ran, per route.

    Each route leaves a different trace, and every branch here is one a
    mistake would show up in: the registry records the namespaced calls, the
    merge writes a file, and the hook's own product is a decision (it never
    runs the call — the external program does what the decision says), so for
    that route an ``allow`` IS the call having been let through. There is no
    default "yes" branch: an unmeasured combination fails rather than passes.
    """
    if route == "cli":
        try:
            return json.loads(str(out))["decision"] == "allow"
        except (TypeError, ValueError, KeyError):
            return False
    if tool == "worktree_merge":
        return _merged(env)
    if route == "stream" and tool in META_TOOLS:
        return any(c[0] in ("read_resource", "get_prompt") for c in reg.calls)
    if tool.startswith("mcp__"):
        return any(c[0] == "call" for c in reg.calls)
    raise AssertionError(
        f"no observable for {route}/{tool}: add one rather than letting the "
        "row pass without evidence")


# ---------------------------------------------------------------------------
# 1. A call the gate ends does not run, on any route
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("route,tool", MEASURED)
def test_a_call_the_policy_denies_does_not_run_on_this_route(route, tool, env,
                                                             monkeypatch):
    """A managed deny on the name ends the call on every route.

    ``mcpdisp`` and ``cli`` are the controls (green on both trees — both have
    always consulted an entry). ``merge`` and ``orphan`` are RED on the base
    commit for ``worktree_merge``: they reach the executor method below both
    entries, so the rule is never read and the merge lands. The two meta-tool
    rows on ``stream`` are RED for the same reason.
    """
    out, reg = _drive(route, tool, env, monkeypatch)
    assert _refused_by_rule(out), (
        f"{route} ({ROUTES[route]}) ran '{tool}' although the administrator "
        "denies it — this route does not pass the gate, so a policy the "
        f"installation set cannot stop it. It returned: {str(out)[:500]}")
    assert not _reached(out, reg, route, tool, env), (
        f"{route} refused '{tool}' and reached its runner anyway: {reg.calls}")
    if tool == "worktree_merge":
        assert not _merged(env), (
            f"{route} merged the worktree's change into the target despite "
            "the denial")


@pytest.mark.parametrize("route,tool", MEASURED)
def test_each_route_either_runs_the_call_or_refuses_for_its_stated_reason(
        route, tool, env, monkeypatch):
    """CONTROL, green on both trees: with the policy empty the route is live.

    Either the call reaches the executor / registry, or the route refuses for
    the one reason ``DOCUMENTED_REFUSAL`` names. Without this row, the gate row
    above would be satisfied by a route that never dispatches anything.
    """
    out, reg = _drive(route, tool, env, monkeypatch, denied=False)
    text = str(out)
    assert not _refused_by_rule(out), text[:500]
    expected = DOCUMENTED_REFUSAL.get((route, tool))
    if expected is not None:
        assert expected in text, (
            f"{route}/{tool}: expected its documented refusal {expected!r}, "
            f"got {text[:500]}")
        return
    assert _reached(out, reg, route, tool, env), (
        f"{route}/{tool} neither ran the call nor stated a reason: "
        f"{text[:500]}")
    if tool == "worktree_merge" and route in ("merge", "orphan"):
        assert _merged(env), text[:500]


# ---------------------------------------------------------------------------
# 2. Every registered tool has an executor the single entry can reach
# ---------------------------------------------------------------------------


def _dispatch_body() -> ast.FunctionDef:
    """``_DocToolExecutor._dispatch`` — the table the entry runs."""
    tree = ast.parse(Path(A.__file__).read_text(encoding="utf-8"))
    for cls in ast.walk(tree):
        if not isinstance(cls, ast.ClassDef):
            continue
        for m in cls.body:
            if isinstance(m, ast.FunctionDef) and m.name == "_dispatch":
                return m
    raise AssertionError(
        "_dispatch not found: the entry's dispatch table moved — re-read the "
        "entry and update this test rather than trusting the old shape")


def _executed_methods() -> set[str]:
    """Every ``_execute_*`` method the executor class defines."""
    tree = ast.parse(Path(A.__file__).read_text(encoding="utf-8"))
    for cls in ast.walk(tree):
        if isinstance(cls, ast.ClassDef) and cls.name == "_DocToolExecutor":
            return {m.name for m in cls.body
                    if isinstance(m, (ast.FunctionDef, ast.AsyncFunctionDef))
                    and m.name.startswith("_execute_")}
    raise AssertionError("_DocToolExecutor not found")


def test_every_registered_tool_reaches_an_executor_from_the_entry():
    """A name on the surface executes from the entry — not merely is listed.

    Two changes in this repository added tool names to an allow-list, passed
    their own tests, and every call answered "Unknown tool". The entry's real
    contract is ``getattr(self, f"_execute_{name}")`` plus inline branches, so
    a row here fails for a registered name with neither — a tool the model can
    call and DELFIN cannot run.
    """
    methods = {m[len("_execute_"):] for m in _executed_methods()}
    assert methods, "no executor methods found — the scan would be vacuous"
    src = Path(A.__file__).read_text(encoding="utf-8")
    dispatch_src = src.split("def _dispatch", 1)[1].split(
        "def _execute_calc", 1)[0]
    unrunnable = []
    for name in sorted(registered_tools()):
        if name in META_TOOLS or name.startswith("mcp__"):
            continue
        if name in methods or f'"{name}"' in dispatch_src:
            continue
        unrunnable.append(name)
    assert not unrunnable, (
        f"registered with no executor the entry can reach: {unrunnable}")
    assert "worktree_merge" in methods, (
        "the verb this package's routes execute is not an executor method")


def test_the_surface_covers_both_entries():
    """The enumeration is the real surface, and both entries exist to judge it.

    Guards the claim "every tool" against an empty or stale surface: a
    catalogue that stopped exporting would otherwise make the rows above
    vacuous.
    """
    names = registered_tools()
    assert len(names) > 50, len(names)
    assert {"bash", "write_file", "worktree_merge"} <= names
    assert META_TOOLS <= names
    assert callable(A.execute_mcp_call), "the MCP entry moved"
    assert callable(A._mcp_call_refusal), "the MCP gate moved"
    assert callable(_DocToolExecutor.execute), "the native entry moved"


# ---------------------------------------------------------------------------
# 3. No route reaches an executor method below the entry
# ---------------------------------------------------------------------------


def _class_at(lines, lineno):
    for i in range(lineno - 1, -1, -1):
        if lines[i].startswith("class "):
            return lines[i][len("class "):].split("(")[0].split(":")[0].strip()
    return ""


def test_no_route_calls_a_tool_method_below_the_entry():
    """No module outside the entry runs a tool by calling its method.

    This is the scan that finds the defect this package is about: a route
    reaches ``_execute_<tool>`` directly, below ``execute()`` — skipping role,
    permission rules, the session lists, the mode and the hooks, none of which
    can be re-added at the method without duplicating them. Call sites are read
    from the AST of every file under ``delfin/``, so a new route turns this red
    in the commit that adds it. Calls inside the executor class itself are the
    entry's own dispatch, which is allowed.

    This file is skipped by path: it drives the bypass on purpose, to show it
    is refused. That is the only exception, and it is named rather than hidden
    in a comprehension.
    """
    methods = _executed_methods()
    assert methods, "no executor methods found — the scan would be vacuous"
    here = Path(__file__).resolve()
    offenders = []
    for py in sorted((REPO / "delfin").rglob("*.py")):
        text = py.read_text(encoding="utf-8", errors="replace")
        if "_execute_" not in text:
            continue
        try:
            tree = ast.parse(text)
        except SyntaxError:                            # pragma: no cover
            continue
        lines = text.splitlines()
        for node in ast.walk(tree):
            if not isinstance(node, ast.Call):
                continue
            f = node.func
            if not (isinstance(f, ast.Attribute) and f.attr in methods):
                continue
            if py.name == "api_client.py" \
                    and _class_at(lines, node.lineno) == "_DocToolExecutor":
                continue
            offenders.append(f"{py.relative_to(REPO)}:{node.lineno} {f.attr}")
    assert not offenders, (
        "routes that run a tool below the gate — route them through "
        f"_DocToolExecutor.execute: {offenders}")
    assert here.name not in {o.split(":")[0].split("/")[-1] for o in offenders}


# ---------------------------------------------------------------------------
# 4. Each route is a live entry point (so row 1 measures a gate, not a stub)
# ---------------------------------------------------------------------------


def test_the_code_route_is_the_server_a_program_talks_to(env, monkeypatch):
    """The Code Mode route answers an allow-listed tool and refuses the rest.

    ``handle_payload`` is what the socket's serve loop calls, so the route
    measured is the one ``_code_mode_child_env`` starts for a caged child.
    """
    reg = _Registry()
    monkeypatch.setattr(MC, "get_registry", lambda *a, **k: reg)
    server = CM.RpcServer(_DocToolExecutor(), permissions=_perms(env))
    allowed = server.handle_payload(json.dumps(
        {"tool": PIPELINE_MCP, "arguments": {}}))
    assert "not one of them" not in allowed, allowed
    refused = server.handle_payload(json.dumps(
        {"tool": "worktree_merge", "arguments": {"path": str(env["ws"])}}))
    assert "not one of them" in refused, refused


def test_the_stream_route_really_dispatches(env, monkeypatch):
    """The streaming route runs a native call end to end.

    Without this, an empty result could be misread as a refusal and the gate
    row would measure nothing. A harmless read-only native tool is dispatched
    by the real loop and its result reaches the caller.
    """
    out, _reg = _drive("stream", "check_environment", env, monkeypatch,
                       denied=False)
    assert out.strip(), out[:300]
    assert not _refused_by_rule(out), out[:300]


def test_the_merge_route_applies_a_real_merge(env, monkeypatch):
    """The session-merge route merges when nothing refuses.

    The control case of the merge rows: it proves ``merge`` writes, so its
    refusal row measures a gate and not a path that never worked.
    """
    out, _reg = _drive("merge", "worktree_merge", env, monkeypatch,
                       denied=False)
    assert _merged(env), out


def test_the_orphan_route_applies_a_real_merge(env, monkeypatch):
    """The orchestration's merge step merges too (same control, other route)."""
    out, _reg = _drive("orphan", "worktree_merge", env, monkeypatch,
                       denied=False)
    assert _merged(env), out


def test_the_cli_route_prints_a_decision_for_a_native_call(env, monkeypatch):
    """The hook answers for a native call: the control's control.

    A call DELFIN refuses is denied; a call it checked and admits is allowed,
    so the external program does not ask a second time.
    """
    out, _reg = _drive("cli", FOREIGN_MCP, env, monkeypatch, denied=False)
    data = json.loads(out)
    assert data["decision"] in ("allow", "deny"), data


# ---------------------------------------------------------------------------
# 5. The repository's own rule, applied to this file
# ---------------------------------------------------------------------------


def test_no_personal_paths_or_accounts_in_this_file():
    """No home path, account name or host name in a committed file (rule 16).

    The needles are built from parts, so the checking line is not itself a
    match: a test that quotes the string it scans fails on its own source and
    teaches whoever reads the red run nothing.
    """
    text = Path(__file__).read_text(encoding="utf-8")
    assert str(Path.home()) not in text, "a home path leaked into a test"
    sep = os.sep
    needles = [sep + "h" + "ome" + sep,
               sep + "U" + "sers" + sep,
               "n" + "oreply",
               "@" + "users"]
    hit = [n for n in needles if n in text]
    assert not hit, f"personal-data-shaped literal in the file: {hit}"
