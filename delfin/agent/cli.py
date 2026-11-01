"""CLI entrypoint for the .delfin agent — installed as ``delfin-agent``.

Two front doors, one parser. ``chat`` is the agent in the directory you
are standing in and is what a bare invocation routes to; ``run`` and the
rest are the headless commands, unchanged, and are what CI hooks, the
scheduler daemon and the benchmark drive.

Examples::

    # the agent, here, in this directory
    delfin-agent

    # one prompt, answer on stdout, exit
    delfin-agent -p "summarise the changes since main"
    echo "list failing tests" | delfin-agent

    # headless: one turn, auto-saved so --session continues it
    delfin-agent run "summarise the changes since main"
    delfin-agent run --session latest "any unresolved TODOs?"
    delfin-agent run --json "list failing tests"

    # init a fresh project
    delfin-agent init /path/to/repo

``python -m delfin.agent.cli <subcommand>`` keeps working for every
subcommand; only the implicit-``chat`` routing is specific to the
installed console script.
"""

from __future__ import annotations

import argparse
import datetime
import json
import os
import sys
import threading
from pathlib import Path
from typing import Any


def _engine_defaults() -> dict[str, str]:
    """Backend/provider/model/effort the user last chose in the dashboard.

    Read from ``~/.delfin_settings.json`` the same way the dashboard reads
    it. Without this the terminal and the dashboard talk to two different
    models out of one settings file, and the terminal picks the constructor
    default (Anthropic) rather than the provider the user actually
    configured. Applied only where the CLI flag is empty.
    """
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent", {}) or {}
    except Exception:
        return {}
    return {k: str(agent.get(k, "") or "")
            for k in ("backend", "provider", "model", "effort")}


def guard_principles_or_exit(pack_dir: "Path | None" = None) -> None:
    """Refuse to start when the principles were changed or removed.

    The engine enforces the same check for every entry point; here the
    terminal gets a clean exit with the reason instead of a traceback.
    """
    from .principles_guard import PrinciplesTampered, enforce
    try:
        enforce(pack_dir)
    except PrinciplesTampered as exc:
        raise SystemExit(f"refusing to start: {exc}") from None


def _resolve_provider(args: argparse.Namespace,
                      fallback: "dict | None" = None) -> str:
    """The provider a session is built with.

    Input: the parsed flags and the settings-file defaults (``fallback``;
    read here when None and the command takes them). Output: ``--provider``
    when given, else the settings file's provider, else the provider the
    model id names (``pricing.infer_provider``: ``kit.``/``azure.`` -> kit,
    ``name:tag`` -> ollama, ``gpt-``/``o1``/``o3``/``o4`` -> openai,
    ``claude`` -> claude), else "". An inferred provider is recorded on
    ``args`` for ``_provider_notes``.

    Rejected: leaving it empty, which sent ``--model kit.<id>`` to the
    Anthropic backend and its key error; a second prefix table here (the
    pricing table is the one answerer for "which provider is this id").
    """
    if fallback is None:
        fallback = (_engine_defaults()
                    if getattr(args, "settings_defaults", False) else {})
    try:
        args.provider_inferred_note = ""
    except Exception:
        pass
    provider = (getattr(args, "provider", "") or ""
                or fallback.get("provider", "") or "")
    if provider:
        return provider
    model = getattr(args, "model", "") or fallback.get("model", "") or ""
    try:
        from .pricing import infer_provider
        inferred = infer_provider(model) if model else ""
    except Exception:
        inferred = ""
    if inferred:
        try:
            args.provider_inferred_note = (
                f"provider   {inferred}, read from the model name {model} "
                "(no --provider given; --provider names another)")
        except Exception:
            pass
    return inferred


def _provider_notes(args: argparse.Namespace) -> list[str]:
    """The start note of a provider read from the model name, or []."""
    note = getattr(args, "provider_inferred_note", "") or ""
    return [note] if note else []


def _build_engine(args: argparse.Namespace):
    """Construct an AgentEngine for the given CLI args.

    AgentEngine creates its own client internally via ``create_client``,
    so we just hand it the resolved (backend, provider, model, mode)
    tuple and let it own the lifecycle.

    ``settings_defaults`` is opt-in per command: ``chat`` passes them so a
    bare ``delfin-agent`` uses the configured provider, while ``run`` keeps
    its historical defaults so the scheduler and the benchmark are not
    moved under anyone's feet by a settings file.
    """
    from .engine import AgentEngine

    fallback = _engine_defaults() if getattr(args, "settings_defaults", False) else {}
    backend = args.backend or fallback.get("backend", "") or "api"
    model = args.model or fallback.get("model", "")
    provider = _resolve_provider(args, fallback)
    mode = getattr(args, "mode", "") or "solo"
    cwd = getattr(args, "cwd", "") or os.getcwd()
    engine = AgentEngine(
        repo_dir=Path(cwd).expanduser().resolve(),
        backend=backend,
        provider=provider,
        model=model,
        mode=mode,
        # Declared on the parser since the first version of this file and
        # never forwarded: --effort parsed fine and changed nothing.
        effort=getattr(args, "effort", "") or fallback.get("effort", ""),
        permission_mode=getattr(args, "permission_mode", "") or "",
        extra_dirs=list(getattr(args, "extra_dirs", None) or []),
        read_only_dirs=list(getattr(args, "read_only_dirs", None) or []),
        # None, never []: CLIClient tests the parameter for ``is not None``
        # and an empty list would build ``--allowedTools`` with nothing
        # after it, which restricts the session to no tools at all.
        allowed_tools=_allowed_tools(args) or None,
    )
    _apply_run_budget(engine, args)
    _apply_tool_surface(engine, args)
    # After the client exists, because two of the three switches live on
    # the permissions object it carries. The MCP half runs earlier, before
    # anything can read the registry — see ``_apply_bare``.
    if getattr(args, "bare", False):
        _apply_bare_discovery(engine)
    return engine


def _allowed_tools(args: argparse.Namespace) -> list[str]:
    """``--allowed-tools a,b,c`` as a list, or empty for "no restriction"."""
    return _tool_name_list(getattr(args, "allowed_tools", ""))


def _disallowed_tools(args: argparse.Namespace) -> list[str]:
    """``--disallowed-tools a,b,c`` as a list, or empty for "deny nothing"."""
    return _tool_name_list(getattr(args, "disallowed_tools", ""))


def _tool_name_list(raw) -> list[str]:
    """A comma-separated flag value as a list of stripped, non-empty names."""
    raw = raw or ""
    if isinstance(raw, (list, tuple)):
        parts = [str(p) for p in raw]
    else:
        parts = str(raw).split(",")
    return [p.strip() for p in parts if p.strip()]


def _apply_tool_surface(engine, args: argparse.Namespace) -> None:
    """Hand this run's tool lists to the client that can enforce them.

    Applied to the constructed client rather than passed down the engine
    constructor, because the lists live on the permissions object the
    executor is handed on every call — that is what makes an excluded tool
    REFUSED when called instead of merely hidden from the model. Backends
    with no tool loop of their own store nothing; ``_bounding_notices``
    reads back what really arrived and says so.

    Nothing is written to any settings file: a bound asked for on this
    command line bounds this session and no later one.
    """
    allowed = _allowed_tools(args)
    denied = _disallowed_tools(args)
    if not allowed and not denied:
        return
    client = getattr(engine, "client", None)
    narrow = getattr(client, "narrow_tool_surface", None)
    if callable(narrow):
        try:
            narrow(allowed=allowed, denied=denied)
        except Exception:
            pass


def _apply_run_budget(engine, args: argparse.Namespace) -> None:
    """Carry ``--max-budget-usd`` / ``--max-run-seconds`` onto the engine.

    ``AgentEngine._run_budget`` reads ``run_budget_usd`` / ``run_budget_s``
    off the instance and falls back to the settings file only when the
    attribute is absent or zero — the precedence the scheduler daemon and
    ``benchmark_runner`` already use to bound one entry without touching
    anybody's configuration. Set the same way here, so a ceiling asked for
    on the command line is a ceiling for THIS run and nothing else.

    The ceiling is cumulative over the session; the per-turn cost breaker
    is a different limit and neither replaces the other. Zero or absent
    leaves the settings value in charge, which is what every caller that
    predates these flags hands over.
    """
    try:
        usd = float(getattr(args, "max_budget_usd", 0.0) or 0.0)
    except (TypeError, ValueError):
        usd = 0.0
    try:
        secs = float(getattr(args, "max_run_seconds", 0.0) or 0.0)
    except (TypeError, ValueError):
        secs = 0.0
    if usd > 0:
        engine.run_budget_usd = usd
    if secs > 0:
        engine.run_budget_s = secs
    try:
        turns = int(getattr(args, "max_run_turns", 0) or 0)
    except (TypeError, ValueError):
        turns = 0
    if turns > 0:
        engine.run_budget_turns = turns

    # `--max-turns` bounds the tool ROUNDS inside one turn, which is a
    # different dimension from the two above: a session can be cheap and
    # still spin, and a single turn can spin without the cumulative
    # ceilings ever being approached. The client resolves it from the
    # user's settings file, so the override lives on the instance — the
    # settings file belongs to every later session, and bounding one run
    # must not change what every run does.
    try:
        rounds = int(getattr(args, "max_turns", 0) or 0)
    except (TypeError, ValueError):
        rounds = 0
    if rounds > 0:
        client = getattr(engine, "client", None)
        if client is not None:
            client.max_tool_rounds = rounds


def _resume_or_create(engine, args: argparse.Namespace) -> str:
    """Restore the engine from a saved session if requested."""
    from . import session_store as _ss
    sid = (args.session or "").strip()
    if not sid:
        return ""
    if sid == "latest":
        # Workspace-scoped, which is what latest_session exists for:
        # continuing a project means continuing THAT project's last
        # conversation, not whatever ran last on the machine. resume_latest
        # asks the whole machine, so `run --session latest` in project A
        # could pick up project B's history and answer out of it.
        _ws = getattr(args, "cwd", "") or os.getcwd()
        data = _ss.latest_session(workspace=_ws) or _ss.resume_latest()
    else:
        data = _ss.load_session(sid)
    if not data:
        print(f"WARN: session '{sid}' not found, starting fresh.", file=sys.stderr)
        return ""
    try:
        report = engine.restore_state({
            # Forward everything the store holds, then apply the defaults
            # below. Listing the fields by hand was the bug: restore_state
            # learned to read the compaction summaries, the project pin and
            # the context floor, and this path kept handing over the same
            # seven keys it always had, so headless --resume came back
            # thinner than the session that was saved.
            **data,
            "mode": data.get("mode", args.mode or "solo"),
            "session_id": data.get("session_id", sid),
        })
    except _ss.SessionSchemaError as exc:
        # A stated refusal. Continuing fresh is safe: the engine keeps the
        # id it minted at construction, so the newer file is not written
        # over by this run.
        print(f"REFUSED: {exc}", file=sys.stderr)
        return ""
    except Exception as exc:
        print(f"WARN: restore_state failed ({exc}); continuing fresh.", file=sys.stderr)
        return ""
    # Say what actually came back. An almost-empty restore used to look
    # exactly like a complete one, and unattended runs have no one watching
    # the screen to notice the session had forgotten everything.
    try:
        if report is not None and (
                not report.complete or report.migrations
                or getattr(args, "verbose", False)):
            print(f"resume: {report.summary()}", file=sys.stderr)
    except Exception:
        pass
    return data.get("session_id", sid)


def _json_line(obj: dict) -> None:
    """One JSON object, one line of stdout, flushed.

    ``--output-format stream-json`` is consumed by a reader that parses
    line by line as the turn runs, so a buffered write would deliver the
    whole stream at the end and make the format pointless. ``sys.stdout``
    is looked up per call rather than captured, because the one-shot path
    is also driven from tests that replace it.
    """
    sys.stdout.write(json.dumps(obj, ensure_ascii=False) + "\n")
    sys.stdout.flush()


def _run_once(engine, prompt: str, *, max_tokens: int = 4096,
                 memory_context: str = "",
                 emit: Any = None,
                 on_tool_result: Any = None,
                 on_permission_denied: Any = None) -> dict[str, Any]:
    """Stream a single turn and collect text + tool-calls + token-usage.

    AgentEngine's ``stream_response`` is callback-driven, not event-
    iterable: text arrives through ``on_token``, tool calls through
    ``on_tool_use``.  The function also returns the assembled text as
    a string.  Token usage is read from ``engine.token_usage`` after
    the call (cumulative for the engine; each benchmark task gets a
    fresh engine so the cumulative IS per-turn).

    ``emit`` receives one dict per event as it happens, for the callers
    that want the turn as a stream rather than as one answer at the end.
    The return value is unchanged either way, so the summary a streaming
    caller prints last is the same object the JSON caller prints alone.

    ``on_tool_result`` (optional) receives ``(name, output)`` per tool
    result -- the harness hints (import origin, refusal reasons, …)
    live there.  It is PASSED to the engine only when given: a stub
    engine with a fixed ``stream_response`` signature (the scheduler
    and benchmark test doubles) must keep working, and the return
    contract stays the keys of ``repl.TURN_KEYS`` -- the benchmark
    runner collects results through this callback instead of a key of
    its own.

    ``on_permission_denied`` (optional) receives ``(tool_name)`` each
    time the backend refused a call -- same forward-only pattern, so
    the benchmark can count interventions without a key of its own.

    ``memory_context`` (optional) goes to the engine's system prompt;
    passed only when non-empty, for the same fixed-signature doubles.
    ``cmd_run`` fills it from ``turn_context``; the scheduler and the
    benchmark call this directly and send none.

    The gate's own decisions are collected by a ``turn_outcome``
    observer for the length of the call: they become
    ``permission_request`` / ``permission_denied`` events on ``emit``
    and the ``schema_version`` / ``stop_reason`` / ``denials`` /
    ``files_changed`` keys of the result.
    """
    from . import turn_outcome

    chunks: list[str] = []
    tool_calls: list[dict] = []
    # Set True once a tool_use event arrives; reset False on each text
    # chunk. True at return means the turn's final text was followed by a
    # tool call, i.e. the turn ended ON a tool round, not on the text.
    # Carried on the ENGINE, not in the returned dict: the dict's keys are
    # pinned to repl.TURN_KEYS by an existing test, so a private
    # attribute is the contract-preserving channel for this internal flag.
    _nonlocal: dict[str, bool] = {"tools_after_last_text": False}
    error = ""

    if emit is not None:
        # Gate decisions arrive on whichever thread runs the tool; one
        # line must not be written into the middle of another.
        _emit_lock = threading.Lock()
        _raw_emit = emit

        def emit(event: dict) -> None:  # noqa: F811 - the locked wrapper
            with _emit_lock:
                _raw_emit(event)

    def _mode() -> str:
        perms = getattr(engine, "kit_permissions", None)
        return str(getattr(perms, "mode", "") or "")

    observer = turn_outcome.TurnObserver(
        emit=emit, mode=_mode,
        workspace=(getattr(getattr(engine, "kit_permissions", None),
                           "workspace", None)
                   or getattr(engine, "repo_dir", None)))

    def _forward_tool_result(name: str, output: str) -> None:
        if on_tool_result is not None:
            try:
                on_tool_result(name, str(output or ""))
            except Exception:
                pass

    def _on_token(text: str) -> None:
        if text:
            chunks.append(text)
            # The turn's LAST text chunk is the final assistant message; a
            # tool call arriving after it means that message was followed
            # by a tool round, so it did not close the turn. Reset here,
            # set in _on_tool_use, read by _continue_announced.
            _nonlocal["tools_after_last_text"] = False
            if emit is not None:
                emit({"type": "text", "text": text})

    def _on_notice(text: str) -> None:
        # Harness speech (retry banners, a stop, a cost ceiling).
        # Deliberately NOT appended to chunks — a notice is not the
        # model's answer, and _run_once's return contract says "text"
        # is the answer. Emitted with its own type so stream-json
        # readers and the CLI renderer can tell the two apart.
        if text and emit is not None:
            emit({"type": "notice", "text": text})

    def _on_wait(text: str) -> None:
        # The first-byte tick. Its own event type, because it differs
        # from a notice in what the reader must DO with it: each tick
        # replaces the previous one, where a notice is a line of the
        # record. A stream-json consumer can now overwrite its status
        # field instead of appending sixty entries that differ only in
        # the number of seconds.
        if text and emit is not None:
            emit({"type": "wait", "text": text})

    def _on_tool_use(name: str, input_json: str) -> None:
        try:
            inp = json.loads(input_json) if input_json else {}
        except (json.JSONDecodeError, TypeError):
            inp = {"raw": str(input_json)}
        tool_calls.append({"name": name, "input": inp})
        _nonlocal["tools_after_last_text"] = True
        observer.saw_tool_use(name)
        if emit is not None:
            emit({"type": "tool_use", "name": name, "input": inp})

    def _on_backend_denied(name: str) -> None:
        # The CLI backend's gate reports a refusal with a tool name and
        # no category, so it gets the one code that says exactly that.
        turn_outcome.denied(name, "backend_denied", "")
        if on_permission_denied is not None:
            on_permission_denied(name)

    in_before = int((getattr(engine, "token_usage", {}) or {}).get("input", 0))
    out_before = int((getattr(engine, "token_usage", {}) or {}).get("output", 0))

    kwargs: dict[str, Any] = dict(
        user_message=prompt,
        on_token=_on_token,
        on_tool_use=_on_tool_use,
        max_tokens=max_tokens,
    )
    # on_notice only when the engine takes it. Two branches of one wave
    # disagreed here: one made the waiting line reach this caller, the
    # other pinned a test that _run_once must keep working against an
    # engine with a FIXED signature -- a test double, or a backend out of
    # this tree. Both are right, and the engine itself already resolves
    # exactly this for `no_tools` the same way.
    #
    # Losing the waiting line on an older engine costs a progress
    # message; passing it costs a TypeError and the whole turn.
    try:
        import inspect as _inspect
        _sig = _inspect.signature(engine.stream_response)
        _takes = ("on_notice" in _sig.parameters or any(
            prm.kind is _inspect.Parameter.VAR_KEYWORD
            for prm in _sig.parameters.values()))
    except (TypeError, ValueError):
        _takes = True
    if _takes:
        kwargs["on_notice"] = _on_notice
    # Same probe, same reason: on_wait is newer than on_notice, so an
    # engine double or an out-of-tree backend may accept one and not the
    # other. Without it the ticks fall back to on_notice and the caller
    # keeps exactly the behaviour it had -- visible, just not in place.
    try:
        import inspect as _inspect2
        _sig2 = _inspect2.signature(engine.stream_response)
        _takes_wait = ("on_wait" in _sig2.parameters or any(
            prm.kind is _inspect2.Parameter.VAR_KEYWORD
            for prm in _sig2.parameters.values()))
    except (TypeError, ValueError):
        _takes_wait = True
    if _takes_wait:
        kwargs["on_wait"] = _on_wait
    if on_tool_result is not None:
        kwargs["on_tool_result"] = _forward_tool_result
    if memory_context:
        kwargs["memory_context"] = memory_context
    # The same probe once more: the backend's refusals are wanted for the
    # result even when no caller asked for them, but only an engine that
    # takes the keyword may be handed it.
    try:
        import inspect as _inspect3
        _sig3 = _inspect3.signature(engine.stream_response)
        _takes_denied = ("on_permission_denied" in _sig3.parameters or any(
            prm.kind is _inspect3.Parameter.VAR_KEYWORD
            for prm in _sig3.parameters.values()))
    except (TypeError, ValueError):
        _takes_denied = True
    if on_permission_denied is not None or _takes_denied:
        kwargs["on_permission_denied"] = _on_backend_denied
    with observer:
        try:
            full_text = engine.stream_response(**kwargs) or ""
        except Exception as exc:
            error = str(exc)
            full_text = ""

    in_after = int((getattr(engine, "token_usage", {}) or {}).get("input", 0))
    out_after = int((getattr(engine, "token_usage", {}) or {}).get("output", 0))
    # The internal "did a tool follow the last text" flag rides on the
    # engine, not the returned dict (whose keys are pinned to TURN_KEYS).
    try:
        setattr(engine, "_tools_after_last_text",
                bool(_nonlocal["tools_after_last_text"]))
    except (AttributeError, TypeError):
        pass

    # The model switches this turn made live on the engine's client
    # (api_client._turn_model_switches: {before, after, reason}). Read them
    # defensively -- a stub engine without a client (the scheduler and
    # benchmark test doubles) must keep working and report []. The read is
    # of the client the engine owns, so `engine.client` is the only route.
    _model_switches = getattr(
        getattr(engine, "client", None), "_turn_model_switches", None)
    model_switches = (list(_model_switches) if _model_switches else [])

    return {
        "text": (full_text or "".join(chunks)).strip(),
        "tool_calls": tool_calls,
        "input_tokens": max(0, in_after - in_before),
        "output_tokens": max(0, out_after - out_before),
        "error": error,
        "model_switches": model_switches,
        **observer.summary(
            error=error,
            engine_stop_reason=str(
                getattr(engine, "last_turn_stop_reason", "") or "")),
    }


#: Continuations for a headless run whose answer ends by announcing work
#: it did not do. Env over ``agent.announce_continuations`` over this
#: default; invalid keeps the default, 0 disables. Read per run, like the
#: other run knobs, never written into anyone's settings file.
_ANNOUNCE_CONTINUATIONS_DEFAULT = 5


def _resolve_announce_continuations() -> int:
    try:
        env = str(os.environ.get("DELFIN_ANNOUNCE_CONTINUATIONS") or "").strip()
        if env:
            return max(0, int(env))
    except (TypeError, ValueError):
        pass
    try:
        from delfin.user_settings import load_settings
        raw = (load_settings() or {}).get("agent", {}) or {}
        val = int(raw.get("announce_continuations",
                          _ANNOUNCE_CONTINUATIONS_DEFAULT))
    except Exception:
        val = _ANNOUNCE_CONTINUATIONS_DEFAULT
    return val if val >= 0 else _ANNOUNCE_CONTINUATIONS_DEFAULT


_ANNOUNCE_CONT_PROMPT = (
    "You ended this run by announcing an action and did not take it. "
    "Take it NOW with the tool call you described, and continue until "
    "the run's task is done. If the action is risky or you are unsure "
    "it is wanted, say clearly in your answer that it is not done "
    "instead of announcing it.")


def _continue_announced(engine, prompt_sent: str, out: dict[str, Any], *,
                        max_tokens: int, emit: Any = None) -> dict[str, Any]:
    """Keep a headless run going while its final reply announces undone
    work.

    Input: ``out`` -- the ``TURN_KEYS``-shaped dict ``_run_once`` returned
    for the one-shot run's turn (``chat -p`` / ``run``), plus the internal
    flag ``tools_after_last_text``; ``engine`` -- the run's engine, whose
    ceilings ``run_limits.ceiling_reason`` reads; ``max_tokens`` -- the
    per-turn token bound cmd_run was given; ``emit`` -- the JSON-line
    renderer or None.

    Output: the LAST turn's dict, same shape. Each continuation is a full
    ``_run_once`` turn: its tokens are billed, its rounds draw the same
    per-turn tool-round budget, and ``ceiling_reason`` counts it against
    the run's ceilings -- the first reached ceiling ends the loop the way
    it would end any further turn.

    The turn's internal ``tools_after_last_text`` flag is read off the
    ENGINE (set by ``_run_once``), not from ``out``: the returned dict's
    keys are pinned to ``repl.TURN_KEYS`` and cannot carry it.

    Semantics: the decision is about the turn's FINAL assistant message
    only, not the turn as a whole -- a run that made many tool calls and
    then closed on an announcement is exactly the failure this fixes, so
    ``out["tool_calls"]`` must not veto. The turn continues when the
    final text announces an action and that text was the turn's last
    output (no tool call followed it); a turn that ended ON a tool round
    is already mid-work and gets none. Bounded by
    ``agent.announce_continuations`` (env ``DELFIN_ANNOUNCE_CONTINUATIONS``
    over settings; default 5; 0 disables). Stops on: a real final answer,
    a ceiling, a non-"done" stop reason, an error, or a no-tool
    continuation that STILL announces -- a real-final no-tool continuation
    is the task finishing, not a stall. Attended sessions never come
    through here: the terminal and the dashboard drive ``repl.run_turn``,
    not ``cmd_run``.

    Rejected: keying the decision on ``out["tool_calls"]`` being empty --
    that is the turn-wide count, and it would have left precisely the
    measured failure (tools, then an announcing close) unfixed. Editing
    api_client's in-turn ``_announced`` guard instead -- that path serves
    attended sessions too, and the headless loop belongs to the CLI.
    """
    allowed = _resolve_announce_continuations()
    try:
        from .api_client import _announces_an_action
    except Exception:
        return out
    if allowed <= 0 or not isinstance(out, dict) or out.get("error"):
        return out
    if str(out.get("stop_reason") or "done") != "done":
        return out
    if not _announces_an_action(str(out.get("text") or "")):
        return out
    # The announcement'd final message must carry NO following tool call.
    # Read off the engine (set by _run_once): a turn that ended ON a tool
    # round is already mid-work and gets no nudge.
    if bool(getattr(engine, "_tools_after_last_text", False)):
        return out
    used = 0
    while used < allowed:
        try:
            from .run_limits import ceiling_reason
            if ceiling_reason(engine):
                break
        except Exception:
            pass
        used += 1
        if emit is not None:
            try:
                emit({"type": "notice", "text":
                      "\n↻ announced an action → continuing the run "
                      f"({used}/{allowed})\n"})
            except Exception:
                pass
        out = _run_once(engine, _ANNOUNCE_CONT_PROMPT,
                        max_tokens=max_tokens, emit=emit)
        tools_after = bool(getattr(engine, "_tools_after_last_text", False))
        if out.get("error") or str(out.get("stop_reason") or "done") != "done":
            break
        if not _announces_an_action(str(out.get("text") or "")):
            # A real answer -- with or without tool calls -- ends it.
            break
        if tools_after:
            continue
        # No tool call AND still announcing: nudging again would only
        # repeat the same sentence.
        break
    return out


def _display_messages(engine, limit: int = 200) -> list[dict[str, str]]:
    """A plain role/content view of the conversation, for the session file.

    Two reasons it is a projection rather than the wire history. The store
    derives a session TITLE from the first user message and searches on it,
    so a message whose content is a list of tool blocks would break the
    title derivation outright. And `chat_messages=[]` is what every
    headless save has passed since this file existed, which is why every
    session it ever wrote is called "Untitled" and invisible to
    `session search`.
    """
    from .harness_rows import history_mark
    out: list[dict[str, str]] = []
    for msg in list(getattr(engine, "messages", []) or [])[-limit:]:
        if not isinstance(msg, dict):
            continue
        role = str(msg.get("role", "") or "")
        if role not in ("user", "assistant"):
            continue
        content = msg.get("content", "")
        if isinstance(content, list):
            parts = [str(b.get("text", "")) for b in content
                     if isinstance(b, dict) and b.get("type") == "text"]
            content = " ".join(p for p in parts if p)
        text = str(content or "").strip()
        if text:
            # Skill learning reads this projection: a harness row keeps
            # its mark (harness_rows).
            out.append(history_mark({"role": role, "content": text}, msg))
    return out


def _save_session(engine, repo_root: Path, *, title: str = "") -> str:
    """Auto-save so the next ``--session`` resumes cleanly.

    The exported state is forwarded WHOLESALE. Listing the fields by hand
    here was the bug — and the mirror image of the one already fixed on the
    restore side. The exporter produces the evidence ledgers, the
    compaction summaries, the project pin, the context floor and the
    system-prompt size; this call named eight keys and none of those, so
    ``save_session`` defaulted every one of them away. The restored session
    then believed it had read no files, run no commands and used no tools,
    while the "a ledger exists" flag was re-derived as True on the first
    turn — the ENFORCING branch. It judged as if it had checked. Every
    scheduled unattended run saves through this same path.
    """
    from . import session_store as _ss
    sid = getattr(engine, "session_id", "") or ""
    if not sid:
        import uuid
        sid = str(uuid.uuid4())
        engine.session_id = sid
    try:
        estate = dict(engine.export_state())
        estate["session_id"] = sid
        if title:
            estate.setdefault("title", title)
        # The model is the client's, not the engine's, so export_state
        # never carried it and every session this command saved recorded
        # an empty one. The dashboard passes both explicitly; a listing
        # that offers sessions to resume showed "?" for anything the CLI
        # had written (2026-09-18).
        estate.setdefault(
            "model", str(getattr(getattr(engine, "client", None), "model", "")
                         or ""))
        estate.setdefault("provider", str(getattr(engine, "provider", "") or ""))
        # The save guard compares this record with the stored one: a new
        # record explains a shorter history. The dashboard passes it; this
        # saver (terminal, -p/run, scheduler) did not, and restore_state
        # reads it back so a resumed session compares against its own.
        estate.setdefault("last_compaction_info",
                          getattr(engine, "last_compaction_info", None))
        _ss.save_session(chat_messages=_display_messages(engine),
                         workspace=str(repo_root), **estate)
    except Exception as exc:
        print(f"WARN: session save failed: {exc}", file=sys.stderr)
    return sid


def _end_session(engine, sid: str, repo: Path, *, task: str, failed: bool,
                 t0: float | None, run_stage=None) -> None:
    """The session-end stages, one list for ``run``, ``chat -p`` and the
    interactive session.

    Input: the engine, the saved session id, the workspace, the task the
    user gave (the -p prompt; the first typed prompt of a terminal
    session, "" when none was typed), whether the session ended on an
    error, and its monotonic start. Output: none; each stage is
    best-effort and a failing one is said on stderr.

    Stages, in order: session report, skill learning, indexing, cycle
    outcome (outcome history + provider profile), episodic memory, eval
    loop. Without a task -- a terminal session left before its first
    prompt -- the cycle outcome and the episode are skipped: there is no
    turn to judge, and a PASS for nothing would enter the history the
    briefing quotes. The terminal used to run the first three only.

    ``run_stage(name, fn)``, when given, runs each stage instead of a
    plain call (the terminal's attended leave: ``repl._LeaveStages``, a
    budget per stage and Ctrl+C to skip). None runs them inline, as for
    ``run``, ``chat -p`` and the scheduler.
    """
    verdict = "FAIL" if failed else "PASS"

    def _session_end_stage(name: str, fn) -> None:
        """One best-effort session-end stage, said when it fails.

        Best-effort stays right — a report write must never break the
        answer's exit code — but a stage that fails on EVERY scheduled
        unattended run vanished silently here: the run printed its
        answer and the failure was indistinguishable from a stage that
        did its job. The warning is the same shape as the save-failure
        line above.
        """
        def _guarded() -> None:
            try:
                fn()
            except Exception as exc:
                print(f"WARN: session-end stage {name} failed: {exc}",
                      file=sys.stderr)

        if run_stage is None:
            _guarded()
        else:
            run_stage(name, _guarded)

    # Session report: best-effort Markdown write; never breaks exit.
    def _report_stage(sid_, engine_):
        from .session_report import write_session_report
        write_session_report(sid_ or getattr(engine_, "session_id", ""))

    def _learn_stage(engine_, sid_):
        from .session_end import learn_at_session_end
        learn_at_session_end(
            _display_messages(engine_),
            session_id=str(sid_ or getattr(engine_, "session_id", "") or ""))

    def _index_stage(sid_):
        from .session_end import index_at_session_end
        index_at_session_end(str(sid_ or ""))

    def _episode_stage(engine_, sid_, repo_):
        from .episodes import build_episode_from_state, save_episode
        fields = build_episode_from_state(engine_.export_state(), [])
        save_episode(
            sid_,
            repo_root=repo_,
            verdict=verdict,
            **fields,
        )

    def _eval_stage():
        from .eval_loop import maybe_run_scheduled
        maybe_run_scheduled()

    _session_end_stage("report", lambda: _report_stage(sid, engine))
    # Shared session-end stage: skill learning, same hook as the CLI
    # chat and the dashboard. Best-effort; never breaks the run.
    _session_end_stage("skill learning", lambda: _learn_stage(engine, sid))
    # Session-end stage (Paket 5): the ended session becomes searchable.
    _session_end_stage("indexing", lambda: _index_stage(sid))

    # Learning signal: record the outcome so provider profiles learn from
    # CLI/headless usage too — previously only dashboard cycles fed the
    # profile, so KIT/CLI sessions contributed nothing.
    if task:
        _session_end_stage(
            "cycle outcome",
            lambda: engine.record_cycle_outcome(
                verdict,
                task,
                error_type=("cli_error" if failed else None),
                start_time=t0,
            ))
    # Episodic memory: persist one compact per-session record so a future
    # session can recall similar past work (best-effort, LLM-free) —
    # without this the saved session state is write-only.
    if task:
        _session_end_stage(
            "episodic memory",
            lambda: _episode_stage(engine, sid, repo))
    # Close the eval loop opportunistically (opt-in, LLM-free, max 1/day).
    _session_end_stage(
        "eval loop",
        lambda: _eval_stage())


def cmd_run(args: argparse.Namespace, *, prepare: Any = None) -> int:
    """One prompt, one answer, then the session-end stages.

    *prepare*, when given, replaces the plain ``--session`` restore: it is
    called with the built engine and opens the session itself (``chat -p``
    passes the interactive session opener and the ``--supervised``
    binding). False from it is a refusal to start, exit 2. ``run``, the
    scheduler and the benchmark pass nothing and keep the plain restore.
    """
    if not args.prompt:
        print("ERROR: prompt is required.", file=sys.stderr)
        return 2
    prompt = " ".join(args.prompt) if isinstance(args.prompt, list) else str(args.prompt)
    repo = Path(args.cwd or os.getcwd()).expanduser().resolve()
    if prepare is None:
        # The `run` door: the bypass rule chat -p applies in
        # _cmd_chat_one_shot, which has already decided the posture.
        mode, note, refused = _headless_posture(args, repo)
        if refused:
            print(note, file=sys.stderr)
            return 2
        if note:
            print(note, file=sys.stderr)
        args.permission_mode = mode
    if args.cwd:
        os.chdir(repo)

    try:
        engine = _build_engine(args)
    except Exception as exc:
        print(f"ERROR: engine init failed: {exc}", file=sys.stderr)
        return 3
    # The headless half has no banner, so the can't-deliver lines go to
    # stderr — stdout carries the answer and nothing else.
    for _note in _provider_notes(args) + _bounding_notices(args, engine):
        print(_note, file=sys.stderr)
    if prepare is None:
        _resume_or_create(engine, args)
        _bind_supervised_approvals(engine, args)
    elif not prepare(engine):
        # chat -p: the shared set-up binds --supervised (_bind_supervised).
        return 2
    # The start check cmd_chat runs too (isolation_posture), once per run
    # for both doors (`run` and chat -p; neither prepare path runs it):
    # without a layer, bypass is confirmed when a terminal can answer,
    # else the run continues in default and says so here. After the
    # session is opened, as in cmd_chat, so it judges the mode the run
    # will actually have.
    from .isolation_posture import apply_at_start
    for _note in apply_at_start(engine, can_ask=sys.stdin.isatty()):
        print(_note, file=sys.stderr)
    # A run nobody watches (--unattended, or no terminal on stdin) gets
    # the ceilings its flags and settings did not set (run_limits).
    from .run_limits import apply_unattended
    for _note in apply_unattended(
            engine, unattended=bool(getattr(args, "unattended", False))
            or not sys.stdin.isatty()):
        print(_note, file=sys.stderr)
    import time as _time
    _t0 = _time.monotonic()
    stream_json = (getattr(args, "output_format", "text") == "stream-json")
    # Passed only when it is wanted: `emit` is a new keyword, and every
    # caller that stands in for `_run_once` was written against the
    # signature without it. A turn that streams nothing is called exactly
    # as it always was.
    # A terminal that is watching gets the turn rendered as it happens --
    # the tool being called, then its verdict -- instead of silence until
    # the answer. cli_stream does the rendering; this only decides when
    # it applies: a tty, and no machine-readable stream asked for.
    _human_stream = bool(not stream_json and sys.stdout.isatty())
    if _human_stream:
        from .cli_stream import StreamRenderer

        # One row for the wait tick. The renderer returns it as a
        # carriage return plus an erase and no newline, so the next tick
        # lands on the same row; printed with the default end="\n" every
        # tick would stand and a stalled round would scroll sixty copies
        # of one sentence past the user. The state lives here and not in
        # the renderer because the renderer is called per event and
        # therefore cannot remember what it printed last.
        _on_wait_row = [False]

        def _render_event(event: dict) -> None:
            is_wait = event.get("type") == "wait"
            if _on_wait_row[0] and not is_wait:
                # Leaving the wait row: erase it so the answer does not
                # start halfway along "waiting for the model … 240s".
                print("\r\x1b[K", end="", flush=True)
                _on_wait_row[0] = False
            for line in StreamRenderer([event], is_tty=True):
                print(line, end="" if is_wait else "\n", flush=True)
                if is_wait:
                    _on_wait_row[0] = True

        _emit = {"emit": _render_event}
    else:
        _emit = {"emit": _json_line} if stream_json else {}
    # What the terminal and the dashboard send with a turn besides the
    # text: the project instruction files and the per-turn hints. One
    # helper for the three surfaces (turn_context).
    from . import turn_context
    _memory = turn_context.memory_context(engine, prompt)
    if _memory:
        _emit = {**_emit, "memory_context": _memory}
    # Over the context budget, the turn starts from a small fresh block,
    # as it does in the terminal (repl.run_turn) and the dashboard: the
    # same helper, under the budget a no-op. A restored session (--session,
    # --continue, --resume) is the case it exists for. Applied after the
    # hint, in run_turn's order.
    from .repl import prepare_fresh_turn
    prompt_sent = prepare_fresh_turn(
        engine, turn_context.user_message(prompt), workspace=repo)
    # 0 is "the role decides" (the engine's max_tokens_for_role), the value
    # chat declares as its default and run --max-tokens 0 asks for. It used
    # to be replaced by 4096 here, so chat -p ran at an eighth of the
    # terminal's solo budget. run keeps 4096 through its parser default.
    out = _run_once(engine, prompt_sent,
                    max_tokens=int(getattr(args, "max_tokens", 0) or 0),
                    **_emit)
    # A one-shot run stops for good after this answer; one that ends by
    # announcing work it did not do must keep going here, bounded (see
    # _continue_announced). Attended sessions are unaffected: they never
    # pass through cmd_run.
    out = _continue_announced(engine, prompt_sent, out,
                              max_tokens=int(getattr(args, "max_tokens", 0)
                                             or 0),
                              emit=_emit.get("emit"))
    sid = _save_session(engine, repo)

    _end_session(engine, sid, repo, task=prompt, failed=bool(out["error"]),
                 t0=_t0)

    from .repl import TURN_KEYS
    if stream_json:
        # The last line closes the stream and carries the same payload the
        # one-object format prints, so a reader that only wants the answer
        # can ignore everything before it. The key set comes from repl's
        # TURN_KEYS rather than from a second list here: that frozenset is
        # what `_run_once` returns and therefore what both formats promise,
        # and two hand-kept lists would drift into two shapes of one answer.
        _json_line({"type": "result",
                    **{k: out.get(k) for k in sorted(TURN_KEYS)},
                    "session_id": sid})
    elif args.json:
        # The same key list as the stream's last line, so the two formats
        # are one object with and without its "type".
        payload = {**{k: out.get(k) for k in sorted(TURN_KEYS)},
                   "session_id": sid}
        print(json.dumps(payload, ensure_ascii=False))
    else:
        if out["text"]:
            print(out["text"])
        if out["error"]:
            print(f"\n[error] {out['error']}", file=sys.stderr)
        if args.verbose:
            print(
                f"\n--- session={sid}  tokens={out['input_tokens']}/"
                f"{out['output_tokens']}  tools={len(out['tool_calls'])}",
                file=sys.stderr,
            )
    # One exit code per stop reason (turn_outcome.EXIT_CODES): 0 and 1
    # keep their meaning for a clean turn and an error; the ceilings and a
    # turn that ended on a refusal used to exit 0 and now say which.
    # A result without a stop reason (an injected _run_once) keeps the
    # old two-way answer.
    from .turn_outcome import EXIT_CODES
    reason = out.get("stop_reason") or ("error" if out["error"] else "done")
    return EXIT_CODES.get(reason, 1)


def _bind_supervised_approvals(engine, args: argparse.Namespace) -> bool:
    """``--supervised`` on a one-shot run: questions go to the approvals
    directory. True when a broker was bound.

    Without the flag nothing is bound and a call that needs a human is
    refused, which is what the scheduler and the benchmark rely on. With
    it the run asks through ``file_confirm.FileConfirmBroker``, whose
    request id the stream-json ``permission_request`` event carries.
    Lives outside ``cmd_run`` so that function binds no dialog of its own.
    """
    if not getattr(args, "supervised", False):
        return False
    return _bind_supervised(engine, args)


def _persisted_default_mode(workspace: Path) -> tuple[str, str]:
    """``(mode, the file that declared it)``, or ``("", "")`` for none.

    Read from the files rather than from ``kit_settings.load()``, and the
    difference is the whole point: with no settings anywhere the loader
    still returns ``default_mode="default"``, its own fallback, which is
    indistinguishable from a user having written that word down. Asking
    the merged view therefore reports a decision nobody made — and it
    silently decided the start posture, because resolve_posture
    correctly honours anything a user configured.

    The path travels with the mode because the banner names it. It used
    to be a literal in resolve_posture's signature, so a checked-out
    repository's settings file could raise the posture and the user was
    pointed at a file in their home directory that need not even exist.
    """
    import json

    # The user's own settings file first (``agent.permission_mode``, what
    # the dashboard's control and a hand edit write): a kit file's
    # ``default_mode`` -- one a permissions tool call persisted, or a
    # checked-out repository's -- does not override what the user set
    # there. Read raw, like the kit files: load_settings() fills in a
    # value nobody chose. A value that is no mode (an old profile name)
    # is no setting.
    try:
        from delfin.user_settings import get_settings_path
        from .launch_guard import _VALID_MODES
        own = get_settings_path()
        if own.exists():
            agent_block = json.loads(own.read_text()).get("agent") or {}
            declared = agent_block.get("permission_mode")
            if isinstance(declared, str) and declared in _VALID_MODES:
                return declared, _tilde(own)
    except Exception:
        pass
    try:
        from . import kit_settings
        candidates = [kit_settings.USER_SETTINGS_PATH,
                      kit_settings.repo_settings_path(workspace)]
    except Exception:
        return "", ""
    for path in candidates:
        if path is None:
            continue
        try:
            if not Path(path).exists():
                continue
            block = json.loads(Path(path).read_text()).get("kit") or {}
        except Exception:
            continue
        declared = str(block.get("default_mode", "") or "")
        if declared:
            return declared, _tilde(path)
    return "", ""


def _claim_session(sid: str) -> bool:
    """Take the writer lock now, not at the first save.

    save_session raises SessionLockedError when a different live process
    holds it, and the headless saver swallows that into one WARN line —
    correct for a scheduled run, wrong for a person, who then works for an
    hour and finds out at the end that nothing was written. Refusing at
    second zero costs nothing and says what to do instead.
    """
    if not sid:
        return True
    from . import session_store as _ss
    try:
        _ss.acquire_session_lock(sid)
    except _ss.SessionLockedError as exc:
        print(
            f"Session {sid} is being written by another process "
            f"(pid {exc.holder_pid}).\n"
            "  It is most likely open in the dashboard.\n"
            f"  delfin-agent --resume {sid} --fork-session   work on a copy\n"
            "  delfin-agent --new                            start fresh here",
            file=sys.stderr)
        return False
    except Exception:
        # A lock we could not take is not a reason to refuse to start; the
        # save path still guards the file itself.
        return True
    import atexit
    atexit.register(_ss.release_session_lock, sid)
    return True


def _pick_session(workspace: Path) -> str:
    """Bare --resume: choose from what this directory has done before."""
    from . import session_store as _ss
    rows = _ss.list_sessions(limit=20, workspace=str(workspace))
    if not rows:
        print("No previous sessions in this directory.", file=sys.stderr)
        return ""
    # The listing comes from cli_resume so the terminal and the picker
    # show the same table -- id, age, model and the task line -- instead
    # of two renderings of the same rows.
    try:
        from .cli_resume import render_sessions
        print(render_sessions(rows), file=sys.stderr)
    except Exception:
        for i, row in enumerate(rows, 1):
            title = str(row.get("title", "") or "(untitled)")[:60]
            when = str(row.get("updated_at", "") or "")[:16]
            print(f"  {i:2}. {when}  {title}", file=sys.stderr)
    try:
        raw = input("resume which? [1] ").strip() or "1"
        idx = int(raw)
    except (ValueError, EOFError, KeyboardInterrupt):
        return ""
    if not 1 <= idx <= len(rows):
        return ""
    return str(rows[idx - 1].get("session_id", "") or "")


def _open_session(engine, args: argparse.Namespace, workspace: Path) -> bool:
    """Resolve -c / -r / --fork-session into a restored, claimed session.

    Returns False when the session could not be opened, which is a refusal
    to start rather than a warning to ignore.
    """
    from . import session_store as _ss

    sid = ""
    if getattr(args, "new_session", False):
        sid = ""
    elif str(getattr(args, "session", "") or "").strip():
        # `--session <id>` is offered on this subcommand and its help says
        # "Session ID to resume". Interactively it was only ever WRITTEN
        # here, never read, so the flag parsed, printed nothing, and
        # started a fresh conversation — the resume that was asked for
        # simply did not happen.
        sid = str(args.session).strip()
    elif getattr(args, "continue_session", False):
        row = _ss.latest_session(workspace=str(workspace))
        sid = str((row or {}).get("session_id", "") or "")
        if not sid:
            print("No previous session in this directory; starting fresh.",
                  file=sys.stderr)
    elif getattr(args, "resume", None) is not None:
        sid = str(args.resume or "").strip() or _pick_session(workspace)
        if not sid:
            return True                      # nothing chosen: a fresh session

    if sid and getattr(args, "fork_session", False):
        try:
            # The copy's id (session_store.fork_session returns it, not a
            # record: reading it as a dict raised, the warning below said
            # "fork failed", and every --fork-session resumed the original).
            # An id of its own: the store's default has one-second steps,
            # and two copies of one session in a second would share a file.
            import uuid
            forked = _ss.fork_session(sid, new_id=uuid.uuid4().hex)
            if isinstance(forked, dict):
                forked = forked.get("session_id")
            if not forked:
                raise ValueError(f"session {sid} could not be read")
            sid = str(forked)
            print(f"Forked to {sid}; the original is untouched.",
                  file=sys.stderr)
        except Exception as exc:
            print(f"WARN: fork failed ({exc}); resuming the original.",
                  file=sys.stderr)

    if sid:
        args.session = sid
        _resume_or_create(engine, args)
        note = None
        try:
            note = _ss.consume_crash_recovery_note(sid)
        except Exception:
            note = None
        if note:
            print(note, file=sys.stderr)

    # One chain of decisions, in one order, so what was decided last
    # is what everything downstream reads:
    #   1. the engine mints or keeps a session id;
    #   2. the session id decides the presence key (the -n name wins,
    #      and the id keeps a nameless session reachable);
    #   3. the presence key is written to the permissions object LAST,
    #      so the operator note (which reads it for the sender) and
    #      session_message agree on which session they speak about.
    # The order is load-bearing in both directions: reading the key
    # before the resume would name the session the id had before the
    # restore, and writing it before the id settles would leave the
    # broker's key and the note's sender disagreeing.
    claimed = _claim_session(getattr(engine, "session_id", "") or sid)
    # Everything downstream -- broker, operator note, session_message
    # -- reads the key through the permissions object, so it is set
    # here, once, after both halves of the name exist. The
    # _presence_key_for call further down (which writes the same
    # value for the broker) stays: it is the same value by the same
    # rule, and removing it would leave the broker's key resting on
    # this distant block.
    try:
        _kp = engine.kit_permissions
        if _kp is not None:
            _kp.presence_key = _presence_key_for(
                getattr(args, "session_name", "") or "",
                str(getattr(engine, "session_id", "") or sid or ""))
    except Exception:
        pass
    return claimed


def _tilde(path: Path | str) -> str:
    """``/home/someone/work/thing`` as ``~/work/thing``.

    The banner is read at a glance, and on this machine the absolute form
    is long enough to push the part that identifies the directory past
    where anyone looks.
    """
    text = str(path)
    try:
        home = str(Path.home())
    except Exception:
        return text
    if home and (text == home or text.startswith(home + os.sep)):
        return "~" + text[len(home):]
    return text


def _parked_work_line(engine, workspace: Path) -> str:
    """One line for the user about open tasks left by earlier sessions.

    The engine already surfaces these — into the SYSTEM PROMPT, where the
    only way they can reach the person at the keyboard is the model
    choosing to read them out. That is how a greeting once came back as a
    fifteen-item backlog. The count belongs on screen, addressed to the
    user, at the one moment it is context rather than an interruption.
    """
    try:
        from .agent_tasks import open_foreign_tasks
        perms = engine.kit_permissions
        if perms is None:
            return ""
        sid = str(getattr(perms, "task_session_id", "") or "")
        summary = open_foreign_tasks(getattr(perms, "workspace", workspace), sid)
    except Exception:
        return ""
    count = int((summary or {}).get("count", 0) or 0)
    if count <= 0:
        return ""
    return (f"parked     {count} open task"
            f"{'s' if count != 1 else ''} from earlier sessions  "
            "(/tasks to list)")


def _launch_questions_answered(report) -> bool:
    """Put every ASK-level finding to the user before the session opens.

    ``LaunchFinding`` has three levels and the middle one is documented as
    "start only if the user says so" — but `LaunchReport.questions` had no
    caller, so an ASK degraded into a paragraph that scrolled past. A
    level whose meaning is enforced nowhere is a comment.

    Bare Enter is not consent, and neither is a pipe: a question nobody
    can answer is answered no. That is why the trust finding is a NOTICE
    rather than an ASK — withholding is already the safe state, so there
    is nothing to decide — but the level now works for a finding where
    there is.
    """
    questions = tuple(getattr(report, "questions", ()) or ())
    if not questions:
        return True
    for finding in questions:
        print(finding.message, file=sys.stderr)
        if finding.detail:
            print(finding.detail, file=sys.stderr)
    if not sys.stdin.isatty():
        print("Not a terminal, so this cannot be answered; not starting.",
              file=sys.stderr)
        return False
    try:
        answer = input("Start anyway? [y/N] ").strip().lower()
    except (EOFError, KeyboardInterrupt):
        print(file=sys.stderr)
        return False
    return answer in ("y", "yes")


#: The name at the top of a session. A box around a word, not a wordmark
#: drawn out of half-blocks: block art is a different size in every
#: terminal font and a smear in the ones that lack the glyphs, and it
#: costs rows of the work that scrolls away. A name is read in a glance
#: either way.
_WORDMARK_TEXT = "DELFIN"

#: Below this the box is wider than the screen, and a wrapped box reads
#: as damage. The line under it says the same thing anyway.
_WORDMARK_MIN_WIDTH = 20


def _terminal_width(default: int = 80) -> int:
    """Never raises: a banner is not worth a traceback."""
    try:
        import shutil as _sh
        return int(_sh.get_terminal_size((default, 24)).columns)
    except Exception:
        return default


def _wordmark(width: int) -> tuple[str, ...]:
    """The name in a box, or nothing when it would not fit."""
    if int(width or 0) < _WORDMARK_MIN_WIDTH:
        return ()
    inner = len(_WORDMARK_TEXT) + 2
    return (
        "╭" + "─" * inner + "╮",
        "│ " + _WORDMARK_TEXT + " │",
        "╰" + "─" * inner + "╯",
        "",
    )


def _startup_banner(engine, report, workspace: Path,
                    why: str = "", isolation_note: str = "",
                    notes: tuple[str, ...] = ()) -> str:
    """What the user is looking at, in the lines that decide safety."""
    from .repl import permission_mode as _permission_mode

    # The model lives on the client; the engine never held one.
    model = str(getattr(getattr(engine, "client", None), "model", "") or "?")
    provider = str(getattr(engine, "provider", "") or "?")
    role_mode = str(getattr(engine, "mode", "") or "?")
    perms_mode = _permission_mode(engine)
    sid = str(getattr(engine, "session_id", "") or "")

    git = report.git
    if git.is_repo:
        if getattr(git, "unborn", False):
            # `git init` and straight in is an ordinary way to start, and
            # it is worth naming: there is a repository, so /rewind works,
            # but nothing to diff against yet.
            where = f"git {git.branch}, no commits yet" if git.branch \
                else "git, no commits yet"
        elif git.branch:
            where = f"git {git.branch}"
        else:
            where = "git, detached HEAD"
        if git.dirty:
            where += f" · {len(git.dirty)} uncommitted"
    else:
        # Said once, plainly. Elsewhere an empty branch used to be printed
        # as a field with nothing after it, which reads as a failed lookup
        # rather than as "this directory is not a repository".
        where = "not a git repository"

    lines = list(_wordmark(_terminal_width()))
    lines += [
        f"delfin-agent · {provider}/{model} · {role_mode}",
        f"workspace  {_tilde(workspace)}  ({where})",
    ]
    # A model too small to drive an agent, said where the model is named
    # (model_capabilities.agent_fitness): the dashboard warned in its
    # preflight, the terminal started on gpt-5-nano without a word.
    try:
        from . import model_capabilities as _mc
        # By name only: the banner must not wait on a model listing.
        _warn = _mc.agent_fitness(provider, model, _mc.ModelCapabilities(
            model=model, provider=provider, context_window=0,
            supports_tools=True))
        if _warn:
            lines.append(f"model      ⚠ {_warn}")
    except Exception:
        pass
    if perms_mode:
        lines.append(f"approval   {perms_mode}"
                     + (f"  [{why}]" if why else ""))
        if any(str(n).startswith("isolation  ") for n in notes or ()):
            pass        # the start check's line says it, with its outcome
        elif isolation_note:
            lines.append(isolation_note)
        else:
            # isolation_status reads the mode (flag, then setting) and the
            # same three probes _bash_isolation_argv picks from, so the
            # line cannot drift from the behaviour it describes -- a
            # setting that switched it off included.
            try:
                from .isolation_posture import isolation_status
                lines.append(f"isolation  {isolation_status(engine.kit_permissions).line()}")
            except Exception:
                lines.append("isolation  unknown — the status could not "
                             "be read")
        # The network of what a command starts, said only where it is not
        # (fully) held: the line exists to keep the banner from implying
        # a protection this host does not give.
        try:
            from .api_client import child_network_enforcement
            _net = child_network_enforcement()
            if not _net["enforced"] or _net["partial"]:
                lines.append(f"network    {_net['notice']}")
        except Exception:
            pass
        # The model second opinion behind the rules: whether it runs, on
        # which model and scope, and its last failure.
        try:
            from .api_client import watcher_status_line
            lines.append("watcher    " + watcher_status_line(
                getattr(engine, "kit_permissions", None)))
        except Exception:
            pass
    else:
        # Not decoration: without a permissions object the write and shell
        # tools refuse, and a user staring at a silent agent deserves the
        # reason on screen rather than in a tool result.
        lines.append(
            f"approval   none — the {provider} backend carries no permission "
            "gate, so file and shell tools will refuse")
    for extra in (_grant_line("writable", getattr(report, "granted_dirs", ())),
                  _grant_line("readable", getattr(report, "read_dirs", ())),
                  *_bounds_in_force(engine),
                  _mcp_note(workspace),
                  _parked_work_line(engine, workspace)):
        if extra:
            lines.append(extra)
    # Every bound the user asked for that this run cannot impose. Beside
    # the approval and grant lines because it belongs to the same
    # question: what is actually stopping this session.
    for note in notes or ():
        if note:
            lines.append(note)
    if sid:
        # The head is what a person types after `-r`; the rest is for the
        # store. Printing all 32 characters put the one useful field on
        # the widest line of the banner.
        lines.append(f"session    {sid[:8]}   (/status for the full id)")
    # What the memory store could shed, said once and only when there is
    # something. Deterministic and measured at 43 ms over 1000 notes, so
    # it costs nothing to ask; a command nobody is told about is a
    # command nobody runs.
    try:
        from .memory_store import _delfin_memory_dir
        from .memory_tidy import hint as _tidy_hint
        _line = _tidy_hint(_delfin_memory_dir(Path(workspace)),
                           repo_root=Path(workspace))
        if _line:
            lines.append(_line)
    except Exception:
        pass
    lines.append("esc interrupt · shift+tab approval mode · /help · ctrl+d exit")
    return "\n".join(lines)


def _mcp_note(workspace) -> str:
    """The configured servers nothing contains, or "".

    The isolation line above it is about the shell DELFIN starts. A tool
    reached over MCP is not that shell, and under bypass — where the
    shell IS wrapped — the difference is at its widest while the banner
    said nothing about it. Reads the config only; starts no server.
    """
    try:
        from . import mcp_client, mcp_isolation
        return mcp_isolation.uncontained_note(
            mcp_client.effective_servers(workspace))
    except Exception:
        return ""


def _grant_line(label: str, dirs) -> str:
    """The extra roots this session was given, or nothing.

    A grant that only exists inside the permissions object is a grant
    nobody can audit: `--add-dir` widens what the agent may write to, and
    the banner is where the user finds out it took.
    """
    paths = [str(_tilde(p)) for p in (dirs or ())]
    if not paths:
        return ""
    shown = ", ".join(paths[:3])
    if len(paths) > 3:
        shown += f", +{len(paths) - 3} more"
    return f"{label:<10} {shown}"


def _usd_ceiling_measurable(engine) -> bool:
    """Whether a USD ceiling can be enforced against THIS model at all.

    ``cost_usd`` is a sum over turns whose price could be looked up. On a
    model with no published rate every turn adds 0.00, so the fraction
    spent stays at 0% for a run of any size and the ceiling never fires.
    The engine says this itself once an unpriced turn has run
    (``_unmeasured_budget_block``); asked here it can be said BEFORE the
    money is spent, which is the only moment the user can still pick a
    dimension that is measurable.
    """
    model = str(getattr(getattr(engine, "client", None), "model", "") or "")
    provider = str(getattr(engine, "provider", "") or "")
    try:
        from .pricing import price_for
        return price_for(model, provider) is not None
    except Exception:
        # A pricing table that cannot be read is not proof of a rate.
        return False


# The discovery ``--bare`` switches off, in the order a session meets it.
# Named once so the help text and the startup line cannot drift from each
# other, and reported per name only when the object that took the switch
# confirms it — see ``_bare_coverage``.
_BARE_SKIPS = ("MCP servers", "tool hooks", "skills", "project memory")

# What it still cannot reach, named the same way and for the same reason.
# PreToolUse and PostToolUse run from the executor, which is handed the
# permissions object carrying the switch (``api_client._session_hooks``);
# the turn-level pair is loaded by the engine from the workspace alone,
# never seeing that object, so it keeps running. A ``--bare`` that implied
# it covered them would be the silent non-delivery the flag exists to
# avoid.
_BARE_NOT_SKIPPED = "UserPromptSubmit and Stop hooks"


def _and_list(names) -> str:
    """``a, b and c`` — for naming what a bound did and did not reach."""
    names = [str(n) for n in names]
    if len(names) <= 1:
        return names[0] if names else ""
    return ", ".join(names[:-1]) + " and " + names[-1]


def _apply_bare_discovery(engine) -> None:
    """Switch off, for THIS session, the discovery that happens in a turn.

    Three loaders read from disk once the turn is running: hook
    definitions at every tool call, skill playbooks while the tool surface
    is assembled, and the memory stores while the system prompt is built.
    Each is switched at the object that does the work — the permissions
    object the executor is handed on every call, and the prompt loader the
    engine owns — the way ``--max-budget-usd`` and ``--max-turns`` are.
    Nothing is written to any settings file: those belong to every later
    session, and this bounds one.

    Not a module global for either half. A global would switch discovery
    off for a dashboard session running in the same process, which nobody
    asked for. The permissions half IS inherited by a sub-agent, because
    ``dataclasses.replace`` copies field values — deliberately: a child
    that re-read what its parent was told to skip would be a route around
    the flag.

    Silent where an object is absent. A backend with no permissions object
    ran none of this in the first place, and the startup line reports what
    actually took rather than what was asked for.
    """
    perms = getattr(engine, "kit_permissions", None) or getattr(
        getattr(engine, "client", None), "_permissions", None)
    if perms is not None:
        try:
            perms.skip_hook_discovery = True
            perms.skip_skill_discovery = True
        except Exception:
            pass
    loader = getattr(engine, "loader", None)
    if loader is not None:
        try:
            loader.skip_external_memory = True
        except Exception:
            pass


def _bare_coverage(args: argparse.Namespace, engine) -> tuple[list[str], list[str]]:
    """``(skipped, unreached)`` for this run, as ``_BARE_SKIPS`` names them.

    Read back off the objects that took the switch, never off the parsed
    arguments: naming something as skipped is a statement that the switch
    arrived, and only the object holding it can say that. The MCP half is
    the exception in form only — its registry is emptied before the engine
    exists, so what took is recorded on the namespace at that moment.
    """
    perms = getattr(engine, "kit_permissions", None) or getattr(
        getattr(engine, "client", None), "_permissions", None)
    loader = getattr(engine, "loader", None)
    took = {
        "MCP servers": bool(getattr(args, "bare_mcp_skipped", False)),
        "tool hooks": bool(getattr(perms, "skip_hook_discovery", False)),
        "skills": bool(getattr(perms, "skip_skill_discovery", False)),
        "project memory": bool(getattr(loader, "skip_external_memory", False)),
    }
    return ([n for n in _BARE_SKIPS if took.get(n)],
            [n for n in _BARE_SKIPS if not took.get(n)])


def _apply_bare(workspace: Path) -> bool:
    """Skip MCP server discovery for this run. True when it took.

    An MCP server definition is executable configuration: it is spawned
    with the parent environment while the tool surface is being ASSEMBLED,
    before any model output, and then answers every call routed to it.
    This is the half that runs before the engine exists, because the
    registry is cached per workspace and loads its configuration once —
    emptying it after that load means no server is started and no MCP tool
    is advertised, for this process only. Nothing on disk is touched, so
    the next run without the flag is unchanged.

    The other three ride on objects the engine builds; they are switched
    in :func:`_apply_bare_discovery`.
    """
    try:
        from . import mcp_client as _mcp
    except Exception:
        return False
    took = False
    # Two keys: the resolved workspace a permissions object carries, and
    # the empty one a backend without a permissions object resolves to.
    for ws in (workspace, None):
        try:
            registry = _mcp.get_registry(ws)
            registry.servers.clear()
            registry.sources.clear()
            registry.trust_notice = ""
            registry.loaded = True
            took = True
        except Exception:
            continue
    return took


def _enforced_tool_surface(client) -> tuple[list[str], list[str]]:
    """What *client* really enforces as ``(allow, deny)``.

    Falls back to the ``allowed_tools`` attribute for anything that predates
    ``enforced_tool_surface`` — a stand-in built by a caller, or a backend
    object that only ever stored the list.
    """
    fn = getattr(client, "enforced_tool_surface", None)
    if callable(fn):
        try:
            allow, deny = fn()
            return ([str(t) for t in (allow or ())],
                    [str(t) for t in (deny or ())])
        except Exception:
            pass
    return ([str(t) for t in (getattr(client, "allowed_tools", None) or ())], [])


def _bounds_in_force(engine) -> list[str]:
    """One banner line per bound this session is actually running under.

    The counterpart of :func:`_bounding_notices`, which says what could
    NOT be honoured. Both are needed and neither substitutes for the
    other: a session where a round cap, a spend ceiling and a denied tool
    set are all in force looked identical on screen to one with no bounds
    at all, so the banner named the approval mode — which decides what may
    happen — and said nothing about the limits deciding what CAN.

    Read off the constructed engine and client, never off the parsed
    arguments. A line here is a statement that the bound arrived, and only
    the object that took it can say that.
    """
    lines: list[str] = []
    client = getattr(engine, "client", None)

    parts: list[str] = []
    rounds = int(getattr(client, "max_tool_rounds", 0) or 0)
    if rounds > 0:
        parts.append(f"{rounds} tool rounds per turn")
    try:
        usd = float(getattr(engine, "run_budget_usd", 0.0) or 0.0)
    except (TypeError, ValueError):
        usd = 0.0
    if usd > 0:
        parts.append(f"${usd:.2f} for the session")
    try:
        secs = float(getattr(engine, "run_budget_s", 0.0) or 0.0)
    except (TypeError, ValueError):
        secs = 0.0
    if secs > 0:
        parts.append(f"{int(secs)}s of wall clock")
    if parts:
        lines.append("bounds     " + " · ".join(parts))

    allow, deny = _enforced_tool_surface(client)
    if allow:
        lines.append(f"tools      only {', '.join(sorted(allow))}")
    if deny:
        lines.append(f"tools      {', '.join(sorted(deny))} denied")
    return lines


def _bounding_notices(args: argparse.Namespace, engine) -> list[str]:
    """One line per bounding flag this run cannot actually honour.

    Same contract as the ``--isolate`` note: a flag that names a limit and
    then does not impose it is worse than no flag, because the user stops
    watching. Every line here is produced from what the constructed engine
    and client really are, never from a second copy of the rules that
    decided it.
    """
    notes: list[str] = []

    # The administrator's managed policy: what it enforces, and a requested
    # mode it replaced. Read from the policy itself, as the gate reads it.
    try:
        from . import managed_policy as _mp
        requested = str(getattr(args, "permission_mode", "") or "")
        if requested:
            _m, _note, _refuse = _mp.launch_decision(
                requested, explicit=True, attended=False)
            if _note:
                notes.append(_note)
        if _mp.isolation_unmet():
            notes.append("managed policy requires an isolation layer and "
                         "none works here: this run is capped at default "
                         "and every shell command is refused")
        notes.extend(_mp.summary_lines())
    except Exception:
        pass

    try:
        usd = float(getattr(args, "max_budget_usd", 0.0) or 0.0)
    except (TypeError, ValueError):
        usd = 0.0
    try:
        secs = float(getattr(args, "max_run_seconds", 0.0) or 0.0)
    except (TypeError, ValueError):
        secs = 0.0
    if usd > 0 and not _usd_ceiling_measurable(engine):
        model = str(getattr(getattr(engine, "client", None), "model", "")
                    or "this model")
        line = (f"budget     ${usd:.2f} REQUESTED but not enforceable — "
                f"{model} has no published rate, so spend cannot be "
                f"measured and the ceiling can never fire")
        line += ("; the wall-clock ceiling IS in force"
                 if secs > 0 else " (--max-run-seconds is measurable)")
        notes.append(line)

    # Tool lists the constructed client did not take. Asked of the client
    # itself rather than by repeating create_client's branch here: each
    # backend enforces what it can and reports it, so comparing against what
    # actually arrived is the only form of this check that cannot go stale.
    client = getattr(engine, "client", None)
    got_allow, got_deny = _enforced_tool_surface(client)
    provider = str(getattr(engine, "provider", "") or "?")
    backend = str(getattr(engine, "backend", "") or "?")
    for flag, wanted, got in (("--allowed-tools", _allowed_tools(args), got_allow),
                              ("--disallowed-tools", _disallowed_tools(args), got_deny)):
        if wanted and set(got) != set(wanted):
            notes.append(
                f"tools      {flag} REQUESTED but not applied — the "
                f"{provider}/{backend} backend carries no such tool list, so "
                f"all {len(wanted)} named tools and every other one stay "
                f"available")

    # A name that matches no tool is a typo, and a typo is silent in both
    # directions: on the allow list it narrows the session further than
    # anyone meant, on the deny list it stops nothing. Said at startup, not
    # discovered from a turn that could not do its work.
    try:
        from .api_client import unknown_tool_names
        stray = unknown_tool_names(_allowed_tools(args) + _disallowed_tools(args))
    except Exception:
        stray = []
    if stray:
        notes.append(
            f"tools      {', '.join(stray)} — no tool of that name exists, so "
            f"{'these names' if len(stray) > 1 else 'this name'} narrows "
            f"nothing and is most likely a typo")

    if getattr(args, "bare", False):
        skipped, unreached = _bare_coverage(args, engine)
        if not skipped:
            notes.append(
                "bare       REQUESTED but nothing was skipped — none of "
                f"{_and_list(_BARE_SKIPS)} could be reached here, so this "
                "session discovers exactly what it would without the flag")
        else:
            line = f"bare       {_and_list(skipped)} skipped"
            if unreached:
                line += (f"; {_and_list(unreached)} REQUESTED but not "
                         "reachable on this backend, so they are unchanged")
            notes.append(f"{line}; {_BARE_NOT_SKIPPED} still run")

    # Permission rules that are not in force: a broken deny or ask is a
    # protection the user believes they have, a broken allow grants nothing.
    try:
        from . import permission_rules
        where = (getattr(engine, "workspace", None)
                 or getattr(args, "cwd", "") or Path.cwd())
        rules = permission_rules.load(where)
        for line in permission_rules.error_lines(rules):
            notes.append(f"permissions {line}")
        scope = permission_rules.scope_notice(rules)
        if scope:
            notes.append(f"permissions {scope}")
        from .api_client import unknown_tool_names
        named = [r.value for r in rules.deny + rules.ask + rules.allow
                 if r.kind == "tool"]
        for name in unknown_tool_names(named):
            notes.append(
                f"permissions rule {name!r} names no tool, so it decides "
                "nothing and is most likely a typo")
    except Exception:
        pass

    return notes


def _apply_launch_grants(args: argparse.Namespace, report) -> None:
    """Hand the engine only the grants launch_guard accepted.

    It has already refused a credential directory, a parent of the
    workspace and a forbidden root by name, rather than letting the
    permissions object drop them silently the way it does. Shared by the
    interactive and the one-shot half of ``chat``: before, a ``-p`` run
    checked ``--add-dir``/``--read-dir`` and then granted nothing.
    """
    args.extra_dirs = [str(p) for p in report.granted_dirs]
    args.read_only_dirs = [str(p) for p in report.read_dirs]


def _apply_isolation_flags(args: argparse.Namespace) -> str:
    """Apply ``--isolate`` / ``--no-isolate``; return the line that says so.

    Output: "" when nothing needs saying. The line goes to the banner of
    an interactive session and to stderr on a one-shot run, because a run
    without containment should not look like one that has it.
    """
    try:
        from . import managed_policy as _mp
        _isolation_required = _mp.load().require_isolation
    except Exception:
        _isolation_required = False
    if getattr(args, "no_isolate", False) and _isolation_required:
        # The managed policy forces isolation in the shell gate whatever
        # this flag says; say so rather than claim it is off.
        return ("isolation  --no-isolate IGNORED: the managed policy "
                "requires isolation for every command")
    if getattr(args, "no_isolate", False):
        from .api_client import set_bash_isolation_override
        set_bash_isolation_override("off")
        return ("isolation  OFF by request (--no-isolate) — a command the "
                "agent runs can write outside the workspace")
    if getattr(args, "isolate", False):
        import shutil as _shutil
        from .api_client import set_bash_isolation_override
        set_bash_isolation_override("bwrap")
        # Asking for isolation the host cannot give must not be silent.
        # _bash_isolation_argv falls back to plain /bin/bash when bwrap is
        # missing — correct, since there is nothing else it can do — but a
        # flag that promises containment and delivers none without saying
        # so is the same defect, one layer up.
        if not _shutil.which("bwrap"):
            return ("isolation  REQUESTED but unavailable — bubblewrap is "
                    "not installed here, so commands run unisolated")
    return ""


def _bind_supervised(engine, args: argparse.Namespace) -> bool:
    """Bind the approvals-directory broker; True when the gate took it.

    No terminal, but somebody is watching the approvals directory.
    Without this the Self-Modification Guard has no callback and refuses
    outright, so a headless run cannot touch the agent's own safety layer
    at all. Asked for explicitly with ``--supervised``: a run nobody is
    watching keeps the refusal. ``--confirm-timeout`` bounds the wait
    (0 is the broker's own default).
    """
    bind = getattr(engine, "set_kit_confirm_callback", None)
    if not callable(bind):
        return False
    from .file_confirm import FileConfirmBroker
    file_broker = FileConfirmBroker(
        session_id=str(getattr(engine, "session_id", "") or ""),
        timeout_s=float(getattr(args, "confirm_timeout", 0) or 0))
    if not bind(file_broker.callback):
        return False
    print("[supervised] confirmations go to the approvals directory; "
          "answer with: delfin-agent approvals", file=sys.stderr)
    return True


def _headless_posture(args: argparse.Namespace,
                      workspace: Path) -> tuple[str, str, bool]:
    """``(mode, note, refused)`` for a headless run, before any engine.

    Called by ``chat -p`` (and a piped prompt) and by the ``run``
    subcommand, so both headless doors hold one rule.

    The interactive rule, held for a run nobody watches:

    * ``--permission-mode bypassPermissions`` needs ``--unattended``.
      A terminal session falls back to ``launch_guard.DEFAULT_MODE``
      (``default``); a headless run is refused (exit 2), because a
      fallen-back run would exit 0 having done less than the script
      asked, and nobody reads the banner.
    * A settings default of ``bypassPermissions`` is never inherited.
      Without ``--unattended``: refused. With it: the posture
      ``launch_guard.resolve_posture`` gives the terminal session
      (``DEFAULT_MODE``), said on stderr. Taking the settings default
      there would give a headless run more than an interactive session
      with the same flags.
    * Otherwise the flag, or ``""`` when none is passed, so the engine
      resolves the settings default (or ``default``) exactly as before --
      the same start posture the terminal now has.
    """
    from . import launch_guard

    flag = str(getattr(args, "permission_mode", "") or "").strip()
    unattended = bool(getattr(args, "unattended", False))
    if flag == "bypassPermissions" and not unattended:
        return "", ("REFUSED: --permission-mode bypassPermissions needs "
                    "--unattended, on a headless run as in a terminal "
                    "session."), True
    if flag:
        return flag, "", False
    persisted, source = _persisted_default_mode(workspace)
    if persisted != "bypassPermissions":
        return "", "", False
    source = source or "~/.delfin/settings.json"
    if not unattended:
        return "", (
            f"REFUSED: {source} sets default_mode = \"bypassPermissions\", "
            "and a headless run does not take that default without "
            "--unattended. "
            "Pass --permission-mode plan, default or acceptEdits, or "
            "--permission-mode bypassPermissions --unattended if you mean "
            "it."), True
    mode, _why = launch_guard.resolve_posture(
        flag_mode="", persisted_mode=persisted, unattended_opt_in=True,
        settings_path=source)
    return mode, (
        f"{source} sets default_mode = \"bypassPermissions\"; a settings "
        f"default of bypass is not inherited, so this run starts in {mode}. "
        "Pass --permission-mode bypassPermissions --unattended if you mean "
        "it."), False


def _cmd_chat_one_shot(args: argparse.Namespace, report, workspace: Path,
                       prompt: str) -> int:
    """``chat -p`` / a piped prompt: the interactive set-up, then ``cmd_run``.

    Input: the parsed flags, the accepted launch report, the workspace and
    the prompt. Output: ``cmd_run``'s exit code, or 2 when a check refuses
    before any engine exists.

    The same set-up an interactive session gets, in the same order:
    grants, isolation flags, posture, then -- with the engine built --
    ``_open_session`` (``--session``, ``--continue``, ``--resume``,
    ``--fork-session``, the writer lock) and the ``--supervised`` binding.
    It used to hand the prompt to ``cmd_run`` before all of it, so a
    one-shot run got more than a terminal with the same flags.

    A bare ``--resume`` is refused: the picker needs a terminal, and a
    fresh session in its place is not what was asked for.
    """
    if (getattr(args, "resume", None) == ""
            and not str(getattr(args, "session", "") or "").strip()
            and not getattr(args, "continue_session", False)
            and not getattr(args, "new_session", False)):
        print("REFUSED: --resume without an id needs a terminal to pick "
              "from; on a -p run pass --resume ID or --continue.",
              file=sys.stderr)
        return 2
    mode, note, refused = _headless_posture(args, workspace)
    if refused:
        print(note, file=sys.stderr)
        return 2
    _apply_launch_grants(args, report)
    isolation_note = _apply_isolation_flags(args)
    for line in (isolation_note, note):
        if line:
            print(line, file=sys.stderr)
    args.permission_mode = mode
    args.prompt = [prompt]
    args.session = getattr(args, "session", "") or ""
    args.json = (getattr(args, "output_format", "text") == "json")

    def _prepare(engine) -> bool:
        if not _open_session(engine, args, workspace):
            return False
        if getattr(args, "supervised", False):
            _bind_supervised(engine, args)
        return True

    return cmd_run(args, prepare=_prepare)


class InteractiveStartRefused(Exception):
    """An interactive session that does not start, and why.

    ``message`` is the text for the person ("" when the step that refused
    has said it already), ``code`` the exit status ``cmd_chat`` returns
    (2 refused, 3 the engine could not be built), ``detail`` the
    underlying exception's text where there is one.
    """

    def __init__(self, message: str, code: int = 2, detail: str = ""):
        super().__init__(message)
        self.message = message
        self.code = code
        self.detail = detail


class InteractiveStart:
    """What ``_start_interactive_session`` built: the engine, the reason
    for its starting mode, the isolation flag line and the start notes
    (isolation_posture.apply_at_start)."""

    def __init__(self, engine, why: str, isolation_note: str,
                 start_notes: tuple):
        self.engine = engine
        self.why = why
        self.isolation_note = isolation_note
        self.start_notes = start_notes

    def banner(self, args: argparse.Namespace, report, workspace: Path) -> str:
        return _startup_banner(self.engine, report, workspace, self.why,
                               self.isolation_note,
                               tuple(_bounding_notices(args, self.engine))
                               + self.start_notes)


def _start_interactive_session(args: argparse.Namespace, report,
                               workspace: Path, *,
                               can_ask: bool) -> InteractiveStart:
    """The set-up of an attended, interactive session, after the launch
    directory was inspected and its questions answered.

    Callers: ``cmd_chat`` (``can_ask`` = stdin is a terminal) and the agent
    service's worker (``service.worker.open_real_session``, ``can_ask``
    False: nobody can answer a yes/no before the window attaches). One
    function, so a step added here reaches both entry paths.

    Steps, in order: the managed policy's start decision (a mode named on
    the command line that the policy disallows, or a required isolation
    layer the host lacks, refuses the start); chdir to ``workspace``; the
    launch grants; the isolation flags; the engine; the session open
    (resume / continue / fork / claim); the starting posture
    (launch_guard.resolve_posture, then managed_policy.coerce_mode); the
    isolation start check (isolation_posture.apply_at_start).

    Output: an ``InteractiveStart``. Raises ``InteractiveStartRefused``;
    the working directory is then the one before the call. On success the
    process stays in ``workspace``; the caller moves it back.
    """
    from . import launch_guard

    # The administrator's managed policy, before anything is built: a mode
    # the person named on the command line that the policy disallows, or a
    # required isolation layer this host lacks, ends an attended start here
    # with the reason. A mode from a settings file falls back instead.
    try:
        from . import managed_policy as _mp
        _flag_mode = str(getattr(args, "permission_mode", "") or "")
        _, _managed_note, _managed_refuse = _mp.launch_decision(
            _flag_mode or "default", explicit=bool(_flag_mode),
            attended=True)
    except Exception:
        _managed_note, _managed_refuse = "", False
    if _managed_refuse:
        raise InteractiveStartRefused(_managed_note, 2)

    _cwd_before = os.getcwd()
    os.chdir(workspace)
    try:
        _apply_launch_grants(args, report)
        isolation_note = _apply_isolation_flags(args)
        try:
            engine = _build_engine(args)
        except Exception as exc:
            raise InteractiveStartRefused(f"ERROR: engine init failed: {exc}",
                                          3, detail=str(exc)) from exc
        if not _open_session(engine, args, workspace):
            raise InteractiveStartRefused("", 2)
    except BaseException:
        try:
            os.chdir(_cwd_before)
        except Exception:
            pass
        raise

    persisted, persisted_from = _persisted_default_mode(workspace)
    posture, why = launch_guard.resolve_posture(
        flag_mode=getattr(args, "permission_mode", "") or "",
        persisted_mode=persisted,
        unattended_opt_in=bool(getattr(args, "unattended", False)),
        settings_path=persisted_from or "~/.delfin/settings.json",
    )
    try:
        from . import managed_policy as _mp
        _managed_posture = _mp.coerce_mode(posture)
        if _managed_posture != posture:
            why = (f"{why}; {posture} disallowed by the managed policy, "
                   f"so {_managed_posture}")
            posture = _managed_posture
    except Exception:
        pass
    try:
        # The starting rule, not a switch: no turn is running to be told.
        try:
            engine.set_kit_permission_mode(posture, announce=False)
        except TypeError:
            # An engine without the keyword still gets the posture; only
            # the note is lost. Swallowing the TypeError below would leave
            # the session in the mode it was built with.
            engine.set_kit_permission_mode(posture)
    except Exception:
        pass          # a backend with no gate has no posture to set
    # The same start check as cmd_run: bypass without isolation is
    # confirmed once here, else the session starts in default.
    from .isolation_posture import apply_at_start
    start_notes = tuple(_provider_notes(args))
    start_notes += tuple(apply_at_start(engine, can_ask=can_ask))
    # --unattended: the default run ceilings. Not keyed on can_ask: the
    # window's worker starts with can_ask=False and a person attaches.
    from .run_limits import apply_unattended
    start_notes += tuple(apply_unattended(
        engine, unattended=bool(getattr(args, "unattended", False))))
    return InteractiveStart(engine, why, isolation_note, start_notes)


def _finish_interactive_session(engine, agent, workspace: Path, *,
                                title: str, t0: float | None,
                                learn: bool = True, run_stage=None) -> str:
    """The end of an interactive session: save it, then ``_end_session``'s
    stages (report, skill learning, indexing, cycle outcome, episode, eval
    loop). The task is the first prompt typed in this process, the verdict
    the last turn's. ``learn`` False saves only (a scripted session of the
    agent service's tests). Callers: ``cmd_chat`` and the agent service's
    worker. Returns the saved session id ("" when the save failed)."""
    sid = _save_session(engine, workspace, title=title)
    if not learn:
        return sid or ""
    _last = getattr(agent, "_last_turn_result", None)
    _end_session(
        engine, sid or str(getattr(engine, "session_id", "") or ""),
        workspace,
        task=str(getattr(agent, "first_task", "") or ""),
        failed=bool(getattr(_last, "error", "")),
        t0=t0, run_stage=run_stage)
    return sid or ""


def cmd_chat(args: argparse.Namespace) -> int:
    """The ``delfin-agent`` front door: a session in the current directory.

    The non-interactive half (``-p/--print``, or a piped prompt) is routed
    through ``cmd_run`` rather than reimplemented, so the JSON payload and
    the session-save semantics stay one contract instead of two.
    """
    from . import launch_guard

    prompt = (getattr(args, "print_prompt", "") or "").strip()
    positional = getattr(args, "prompt", None) or []
    if isinstance(positional, list) and positional:
        positional = " ".join(str(p) for p in positional).strip()
    else:
        positional = ""

    piped = not sys.stdin.isatty()
    if not prompt and piped:
        try:
            prompt = sys.stdin.read().strip()
        except Exception:
            prompt = ""

    workspace = Path(getattr(args, "cwd", "") or os.getcwd()).expanduser().resolve()

    # Before anything is built: is this a directory to work in at all?
    report = launch_guard.inspect_launch_dir(
        workspace,
        add_dirs=tuple(getattr(args, "add_dirs", None) or ()),
        read_dirs=tuple(getattr(args, "read_dirs", None) or ()),
    )
    if report.refused:
        print(report.render(), file=sys.stderr)
        return 2
    if not _launch_questions_answered(report):
        return 2

    # Before the engine exists, because the tool surface is assembled on
    # the first turn and both halves of this command run one. Recorded on
    # the namespace so the startup line reports what actually happened
    # rather than what was asked for.
    if getattr(args, "bare", False):
        args.bare_mcp_skipped = _apply_bare(workspace)

    # One shot, and out. Identical to `run`, because it IS `run`.
    # A positional prompt SEEDS an interactive session; it only becomes a
    # one-shot when there is no terminal to be interactive on.
    one_shot = prompt or (piped and positional)
    if one_shot:
        return _cmd_chat_one_shot(args, report, workspace,
                                  prompt or positional)

    if piped:
        print("Nothing on stdin, and no -p. Nothing to do.", file=sys.stderr)
        return 2

    # An interactive session is a conversation, not a document. The flag
    # is accepted on this subcommand because the one-shot half lives here
    # too; taking it silently and then emitting text is the shape of a
    # promise that is not kept.
    _fmt = getattr(args, "output_format", "text")
    if _fmt in ("json", "stream-json"):
        print(f"--output-format {_fmt} describes one answer, so it needs -p "
              "or a piped prompt; this session will print text.",
              file=sys.stderr)

    # Restored on the way out. The process usually ends right after, so
    # this looks unnecessary — but a function that moves the process and
    # never moves it back makes every later relative path in the caller
    # read a different file, and that is only harmless until it is called
    # from somewhere else.
    _cwd_before = os.getcwd()
    try:
        start = _start_interactive_session(args, report, workspace,
                                           can_ask=sys.stdin.isatty())
    except InteractiveStartRefused as refused:
        if refused.message:
            print(refused.message, file=sys.stderr)
        return refused.code
    engine = start.engine

    from .repl import ReplOptions, TerminalAgent
    from .terminal_confirm import TerminalConfirmBroker

    # Bound ONLY here, and only for a terminal. cmd_run, the benchmark and
    # the scheduler must keep the behaviour they have: without a callback
    # a write inside the workspace is allowed silently, which is right for
    # an unattended run and exactly what this command exists to change.
    broker = None
    bind = getattr(engine, "set_kit_confirm_callback", None)
    if getattr(args, "supervised", False) and not sys.stdin.isatty():
        _bind_supervised(engine, args)    # no terminal UI to drive
    elif sys.stdin.isatty() and callable(bind):
        _sid = str(getattr(engine, "session_id", "") or "")
        _key = _presence_key_for(getattr(args, "session_name", ""), _sid)
        broker = TerminalConfirmBroker(
            # 0 waits for an answer as long as it takes, which is right for
            # someone at the keyboard. A session nobody watches (a wave run
            # in tmux) passes --confirm-timeout: an unanswered question then
            # expires as "not now" -- never as a refusal -- and the session
            # carries on instead of standing still for hours.
            timeout_s=float(getattr(args, "confirm_timeout", 0) or 0),
            persist=lambda pat: engine.persist_kit_pattern(pat, kind="allow"),
            persist_read=lambda d: engine.add_kit_read_dir(d, persist=True),
            set_mode=engine.set_kit_permission_mode,
            # So a question waiting in this pane can be read -- and
            # answered -- from outside it, whole.
            session_id=_sid,
            session_key=_key,
        )
        # False means this provider carries no permissions object at all,
        # and on that backend the file and shell tools refuse outright —
        # so there is nothing to ask about and a broker would only add a
        # prompt that never fires.
        if not bind(broker.callback):
            broker = None
        else:
            perms = engine.kit_permissions
            perms.ask_user_callback = broker.ask_user
            perms.plan_approval_callback = broker.approve_plan
            # session_message reads this. The dashboard set it and the
            # CLI never did, so a terminal session was told it does not
            # exist when it went to warn a parallel one about work it
            # had found in the tree (2026-09-20).
            perms.presence_key = _key


    notices = report.render()
    agent = TerminalAgent(engine, broker=broker, opts=ReplOptions(
        cwd=workspace,
        max_tokens=getattr(args, "max_tokens", 0) or 0,
        show_thinking=bool(getattr(args, "verbose", False)),
        color=getattr(args, "color", "auto"),
        session_name=getattr(args, "session_name", "") or "",
        banner=(start.banner(args, report, workspace)
                + (f"\n\n{notices}" if notices else "")),
    ))
    import time as _time
    from . import repl_keys as _rk
    _t0 = _time.monotonic()
    # Held from the loop's start until the session is saved: every way out
    # is bounded by the leave deadline, and a repeated SIGHUP, SIGTERM or
    # SIGINT does not end the process in the middle of the save.
    with agent.leave_held():
        try:
            return agent.run(positional)
        except _rk.TerminalLeft as left:
            # A signal before the loop's own try took it.
            agent.begin_leaving(left.code, terminal_gone=True)
            return left.code
        finally:
            _finish_interactive_session(
                engine, agent, workspace,
                title=getattr(args, "session_name", "") or "", t0=_t0,
                run_stage=agent.leave_stages())
            try:
                os.chdir(_cwd_before)
            except Exception:
                pass


def _run_prcheck_subcommand(args: argparse.Namespace) -> int:
    """The `prcheck` subcommand body (package PC1): delegate to the
    owned module so the check itself is tested outside cli.py."""
    from .cli_prcheck import main as prcheck_main
    argv = []
    if getattr(args, "base", "origin/main") != "origin/main":
        argv += ["--base", args.base]
    for partner in getattr(args, "partner", []):
        argv += ["--partner", partner]
    return prcheck_main(argv)


def cmd_init(args: argparse.Namespace) -> int:
    from .project_init import init_project
    target = Path(args.path or os.getcwd()).expanduser().resolve()
    try:
        result = init_project(target, overwrite=bool(args.force))
    except Exception as exc:
        print(f"ERROR: /init failed: {exc}", file=sys.stderr)
        return 1
    p = result["profile"]
    print(f"Detected {p.language} project '{p.name}'.")
    if result["files"]:
        print("Created:")
        for f in result["files"]:
            print(f"  + {f}")
    if result["skipped"]:
        print("Skipped (use --force to overwrite):")
        for f in result["skipped"]:
            print(f"  · {f}")
    return 0


def cmd_bench(args: argparse.Namespace) -> int:
    """Canned-task benchmark suite — run, list, or compare."""
    from . import benchmark as _bm
    from . import benchmark_runner as _br
    action = getattr(args, "bench_action", "") or "run"

    if action == "list":
        tasks = _bm.load_tasks()
        if not tasks:
            print("No tasks found.", file=sys.stderr)
            return 1
        for t in tasks:
            print(f"  {t.id:<28} {t.mode:<10} {t.task_class:<22} "
                  f"{t.prompt[:60].replace(chr(10), ' ')}")
        print(f"\n{len(tasks)} tasks")
        return 0

    if action in ("gate", "set-baseline"):
        return _cmd_bench_baseline(args, action)

    if action == "audit":
        run_path = Path(getattr(args, "run", "")).expanduser().resolve()
        if not run_path.exists():
            # Try looking it up in the default runs dir
            rd = _bm.runs_dir()
            candidate = rd / getattr(args, "run", "")
            if candidate.exists():
                run_path = candidate
            else:
                print(f"ERROR: run file not found: {args.run}",
                      file=sys.stderr)
                return 2
        records = _bm.read_run(run_path)
        if not records:
            print(f"Run file empty or unreadable: {run_path}",
                  file=sys.stderr)
            return 1
        entries = _bm.audit_run(records)
        print(_bm.format_audit_report(entries))
        _print_behavior_rates(_bm, records)
        return 0

    if action == "latest":
        import time as _time
        rd = _bm.runs_dir()
        if not rd.exists():
            print("No runs yet — use `bench run --model X` first.")
            return 0
        files = sorted(
            rd.glob("*.jsonl"),
            key=lambda p: p.stat().st_mtime, reverse=True,
        )
        limit = int(getattr(args, "limit", 10) or 10)
        if not files:
            print("No runs yet — use `bench run --model X` first.")
            return 0
        for f in files[:limit]:
            ts = _time.strftime(
                "%Y-%m-%d %H:%M", _time.localtime(f.stat().st_mtime),
            )
            try:
                n_records = sum(
                    1 for line in f.open(encoding="utf-8")
                    if line.strip() and not line.startswith('{"run_header"')
                )
            except OSError:
                n_records = 0
            print(f"  {ts}  {f.name}  ({n_records} tasks)")
        return 0

    if action == "nightly":
        return _cmd_bench_nightly(args)

    if action == "compare" and (getattr(args, "model", "") or "").strip():
        # History mode: classify one run against the rolling per-task
        # baseline (bench_watch). Pure file comparison — no API spend.
        return _cmd_bench_compare_history(args)

    if action == "compare":
        if not getattr(args, "baseline", "") or not getattr(args, "candidate", ""):
            print("ERROR: pass either two run files (baseline candidate) "
                  "or --model for a rolling-baseline compare.",
                  file=sys.stderr)
            return 2
        baseline = Path(args.baseline).expanduser().resolve()
        candidate = Path(args.candidate).expanduser().resolve()
        if not baseline.exists() or not candidate.exists():
            print(f"ERROR: run file missing: "
                  f"{baseline if not baseline.exists() else candidate}",
                  file=sys.stderr)
            return 2
        cmp = _bm.compare_runs(baseline, candidate)
        # Whether one instrument produced both files. A changed scorer or
        # task pack makes the score difference partly the measure's, so
        # the comparison is refused unless asked for explicitly; a file
        # without a header (written before 2026-10-15) is said, not
        # refused. The settings that differ are always printed: they
        # should be exactly the axis the A/B declares.
        _prov = cmp.get("provenance") or {}
        if not args.json:
            for _line in _bm.format_provenance_check(_prov):
                print(_line)
        if _prov.get("instrument_match") is False and not getattr(
                args, "allow_instrument_change", False):
            print("ERROR: comparison refused -- the instrument differs "
                  "between the arms. Re-run one arm, or pass "
                  "--allow-instrument-change to compare anyway.",
                  file=sys.stderr)
            return 4
        if args.json:
            print(json.dumps(cmp, indent=2, ensure_ascii=False))
            return 0
        if getattr(args, "markdown", False):
            print(_bm.format_compare_markdown(
                cmp, baseline_path=baseline, candidate_path=candidate,
            ))
            return 0
        s = cmp["summary"]
        print(f"Verdict: {cmp['verdict'].upper()}")
        print(f"  overlap: {s['n_overlap']} tasks   "
              f"better: {s['n_better']}   "
              f"worse: {s['n_worse']}   "
              f"neutral: {s['n_neutral']}")
        old, new = s["old"], s["new"]
        print(f"\n  pass-rate  {old['pass_rate']:.0%} → {new['pass_rate']:.0%}")
        print(f"  avg-quality {old['avg_quality']:.1f} → {new['avg_quality']:.1f}")
        print(f"  cost(total) ${old['total_cost_usd']:.4f} "
              f"→ ${new['total_cost_usd']:.4f}")
        print(f"  duration(total) {old['total_duration_s']:.1f}s "
              f"→ {new['total_duration_s']:.1f}s")
        # The cost side of the same two runs. The verdict above is built
        # from quality; these can get worse while it does not move, which
        # is exactly the regression a wave of guard work produces.
        _oc, _nc = old.get("avg_caveats"), new.get("avg_caveats")
        if _oc is not None and _nc is not None:
            _worse = _nc > _oc + 0.25
            print(f"  caveats/answer {_oc:.1f} → {_nc:.1f}"
                  + ("   ← MORE HEDGING, treat as a regression" if _worse
                     else ""))
        _ot, _nt = old.get("avg_output_tokens"), new.get("avg_output_tokens")
        if _ot is not None and _nt is not None:
            print(f"  out-tok/task {_ot:.0f} → {_nt:.0f}"
                  + ("   ← LONGER ANSWERS for the same score"
                     if _ot and _nt > _ot * 1.15 else ""))
        print()
        print(f"  {'task_id':<28} {'cls':<8} {'qual':>9} "
              f"{'Δcost':>9} {'Δdur':>8}")
        for row in cmp["per_task"]:
            mark = {"better": "+", "worse": "-", "neutral": "="}.get(row["class"], "?")
            print(f"  {row['task_id']:<28} {mark} {row['class']:<6} "
                  f"{row['old_quality']:>3}→{row['new_quality']:<3}   "
                  f"{row['d_cost_usd']:+.4f}   {row['d_duration_s']:+.2f}s")
        return 0

    # action == "run"
    tasks = _bm.load_tasks()
    if not tasks:
        print("ERROR: no benchmark tasks found.", file=sys.stderr)
        return 1
    _asked = getattr(args, "task", None) or []
    if isinstance(_asked, str):          # a caller that set the field itself
        _asked = [_asked]
    # Exact ids and shell globs (quote them: --task 'generalist_*').
    tasks, missing = _bm.select_tasks(tasks, _asked)
    if not tasks:
        print(f"ERROR: no task matched {sorted(missing)}", file=sys.stderr)
        return 1
    # A name that matched nothing is a typo or a renamed task, and
    # running the rest without saying so hands back a number for a
    # different question than the one that was asked.
    if missing:
        print(f"WARNING: {len(missing)} requested task(s) do not exist "
              f"and were skipped: {', '.join(sorted(missing))}",
              file=sys.stderr)
    if (getattr(args, "lang", "de") or "de") == "en":
        _en = [_bm.english_variant(t) for t in tasks]
        _no_en = [t.id for t, e in zip(tasks, _en) if e is None]
        tasks = [e for e in _en if e is not None]
        if _no_en:
            print(f"WARNING: {len(_no_en)} task(s) have no English prompt "
                  f"and were skipped: {', '.join(_no_en)}", file=sys.stderr)
        if not tasks:
            print("ERROR: none of the selected tasks has an English prompt",
                  file=sys.stderr)
            return 1

    model = (getattr(args, "model", "") or "").strip()
    if not model:
        print("ERROR: --model is required for `bench run`",
              file=sys.stderr)
        return 2
    # 0 = the role's own output limit, as in a session (engine
    # _ROLE_MAX_TOKENS; solo 32768). 1024 cut every write_file longer than
    # about four kilobytes into invalid JSON in the benchmark only.
    max_tokens = int(getattr(args, "max_tokens", 0) or 0)
    backend = getattr(args, "backend", "") or "api"
    provider = getattr(args, "provider", "") or ""
    repeats = max(1, int(getattr(args, "repeats", 1) or 1))
    profile_name = _br.resolve_profile_name(model)

    def _progress(task, result):
        mark = "PASS" if result.success else "FAIL"
        extra = ""
        if result.n_samples > 1:
            extra = (f"  N={result.n_samples}  "
                     f"σ={result.quality_stdev:.1f}  "
                     f"rate={result.success_rate:.0%}  "
                     f"pass@k={result.pass_at_k}  "
                     f"pass^k={result.pass_pow_k}")
        print(f"  [{mark}] {task.id:<28} q={result.quality_0_100:>3}  "
              f"{result.duration_s:>5.1f}s  ${result.cost_usd:.4f}  "
              f"tool={result.tool_calls}{extra}",
              flush=True)
        _proc = getattr(result, "process", None) or {}
        if _proc:
            # Where the run went, not only whether it arrived: the
            # generalist tasks' trace signals (means over replicates).
            _vq = getattr(result, "verify_quality", None)
            print("        process: "
                  + (f"accept-q={_vq} " if _vq is not None else "")
                  + " ".join(f"{k}={_proc[k]:g}" for k in (
                      "planned", "scouted", "tested_after_last_edit",
                      "passing_test_after_last_edit", "test_runs",
                      "refused_test_runs", "edits", "measure_edits",
                      "measure_changed", "tool_errors", "nonzero_exits",
                      "tool_budget_used") if k in _proc),
                  flush=True)
        if result.violated_signals or result.missing_signals or result.error:
            for v in result.violated_signals:
                print(f"        violated: {v}", flush=True)
            _flaky = set(getattr(result, "flaky_signals", None) or ())
            for m in result.missing_signals:
                if m.endswith(":optional"):
                    continue
                # A signal that matched in some replicates is a different
                # fact from one that never matched, and printing them the
                # same way sends the reader to the wrong defect.
                label = "flaky:   " if m in _flaky else "missing: "
                print(f"        {label} {m}", flush=True)
            if result.error:
                print(f"        error:    {result.error[:120]}", flush=True)

    rep_note = f" × {repeats} repeats" if repeats > 1 else ""
    print(f"Benchmarking {model} (profile={profile_name}) "
          f"on {len(tasks)} tasks{rep_note}…")

    def _on_rep(task, idx, result):
        if repeats > 1:
            mark = "P" if result.success else "F"
            print(f"    sample {idx + 1}/{repeats} [{mark}] "
                  f"{task.id:<28} q={result.quality_0_100:>3}  "
                  f"{result.duration_s:>5.1f}s",
                  flush=True)

    _learning_home = getattr(args, "learning_home", "") or None
    if _learning_home is not None:
        try:
            _learning_home = _br._learning_home_checked(_learning_home)
        except ValueError as exc:
            print(f"\nbenchmark not started: {exc}", file=sys.stderr)
            return 2
        print(f"  learning home: {_learning_home} (kept across runs)")

    try:
        results = _br.run_suite(
            tasks,
            model=model,
            backend=backend,
            provider=provider,
            profile_name=profile_name,
            max_tokens=max_tokens,
            repeats=repeats,
            progress=_progress,
            on_replicate=_on_rep if repeats > 1 else None,
            learning_home=_learning_home,
        )
    except _br.BenchmarkRunInProgress as exc:
        # The lock's message already names the holder and the two ways
        # out. It was reaching the caller as a twenty-line traceback with
        # the sentence at the bottom, which is where a queued script's
        # log buries it -- three arms of an A/B died this way and the
        # output read as a crash rather than as a queue.
        print(f"\nbenchmark not started: {exc}", file=sys.stderr)
        return 3
    path = _bm.write_run(results, model=model, run_args={
        "max_tokens": max_tokens, "repeats": repeats, "backend": backend,
        "provider": provider})
    s = _bm.summarise_run(results)
    _pm = s.get("process_means") or {}
    if _pm:
        # Rates per model, not only per task line.
        print("process (mean over tasks): " + "  ".join(
            f"{k}={v['mean']:g} (n={v['n']})" for k, v in _pm.items()
            if k in ("tested_after_last_edit",
                     "passing_test_after_last_edit", "refused_test_runs",
                     "measure_edits", "measure_changed", "edits")))
    if s.get("channels_total"):
        print("channels fired: " + "  ".join(
            f"{k}={v}" for k, v in s["channels_total"].items()))
    print(f"\n{s['n_pass']}/{s['n_tasks']} passed "
          f"({s['pass_rate']:.0%})   "
          f"avg-quality {s['avg_quality']:.1f}   "
          f"total ${s['total_cost_usd']:.4f}   "
          f"{s['total_duration_s']:.1f}s")
    # The stricter, alongside readings over the same repeats (empty only
    # in a no-op run, so guard the keys a summarise_run always now emits).
    if s.get("n_tasks"):
        print(f"   pass@k {s['n_pass_at_k']}/{s['n_tasks']} "
              f"({s['pass_at_k_rate']:.0%} >=1 repeat)   "
              f"pass^k  {s['n_pass_pow_k']}/{s['n_tasks']} "
              f"({s['pass_pow_k_rate']:.0%} ALL repeats)")
    # Said out loud, never folded into the rate: a run that could not
    # reach the model measured nothing, and a summary that hides that is
    # how an outage becomes a baseline.
    if s.get("n_unmeasured"):
        print(f"  ⚠️  {s['n_unmeasured']} task(s) NOT MEASURED — the request "
              f"never reached the model (endpoint unavailable, connection "
              f"or auth error). Excluded from the rate and the average:")
        for _tid in s.get("unmeasured_tasks") or []:
            print(f"       {_tid}")
        print("     Re-run them before comparing this to anything.")
    # The cost side, printed whether or not anything went wrong. The rate
    # above cannot move when an answer grows a third hedge or a guard
    # starts refusing honest work -- these are the only numbers that can,
    # and a number nobody prints is a number nobody compares.
    _denials = s.get("total_denials")
    print(f"  cost side: {s['avg_caveats']:.1f} caveats/answer "
          f"(max {s['max_caveats']})   "
          f"{s['avg_output_tokens']:.0f} out-tok/task   "
          f"{s['avg_answer_chars']:.0f} chars/answer   "
          f"denials {'not observed' if _denials is None else _denials}")
    print("     Rising here is a regression even at an unchanged pass "
          "rate: an answer nobody finishes reading says nothing.")
    # Figures, split by why they did not match. The two failures need
    # different repairs and used to look identical in this summary: a
    # wrong figure is arithmetic, an absent one is usually an answer that
    # never arrived, and reading them as one number files a harness fault
    # as a model fault.
    _v_ok = s.get("values_matched", 0)
    _v_bad = s.get("values_wrong", 0)
    _v_gone = s.get("values_absent", 0)
    if _v_ok or _v_bad or _v_gone:
        print(f"  figures: {_v_ok} matched   {_v_bad} wrong   "
              f"{_v_gone} absent from the answer")
        if _v_gone:
            print("     An absent figure is usually not the model's "
                  "arithmetic — check the answer arrived before blaming it.")
    _print_behavior_rates(_bm, results)
    print(f"\nWritten to: {path}")
    return 0



# ---------------------------------------------------------------------------
# Reference standard (committed) — gate a run, or deliberately move the bar
# ---------------------------------------------------------------------------

def _default_baseline_path() -> Path:
    """The committed reference, resolved relative to the installed package."""
    return (Path(__file__).resolve().parents[2]
            / "tests" / "fixtures" / "office_baseline.json")


def _resolve_run_file(args, _bm) -> Path | None:
    """The run to judge: an explicit file, else the newest for a model."""
    runs_dir = (getattr(args, "runs_dir", "") or "").strip()
    base = Path(runs_dir).expanduser() if runs_dir else _bm.runs_dir()
    raw = (getattr(args, "run", "") or "").strip()
    if raw:
        p = Path(raw).expanduser()
        if p.exists():
            return p.resolve()
        candidate = base / raw
        return candidate.resolve() if candidate.exists() else None
    model = (getattr(args, "model", "") or "").strip()
    if not model:
        return None
    try:
        from . import bench_watch as _bw
        files = _bw.list_model_runs(model, runs_dir=base)
    except Exception:
        files = []
    return Path(files[0]).resolve() if files else None


def _cmd_bench_baseline(args, action: str) -> int:
    """`bench gate` and `bench set-baseline`."""
    from . import benchmark as _bm
    from . import benchmark_baseline as _bb

    raw_path = (getattr(args, "baseline_path", "") or "").strip()
    baseline_path = (Path(raw_path).expanduser() if raw_path
                     else _default_baseline_path())

    run_path = _resolve_run_file(args, _bm)
    if run_path is None:
        print("ERROR: no run file. Give one, or --model to take the newest.",
              file=sys.stderr)
        return 2
    rows = _bm.read_run(run_path)
    if not rows:
        print(f"ERROR: no results in {run_path}", file=sys.stderr)
        return 2

    if action == "set-baseline":
        note = (getattr(args, "note", "") or "").strip()
        commit = ""
        try:
            from . import safe_git
            commit = safe_git.out(Path(__file__).resolve().parents[2],
                                  "rev-parse", "HEAD").strip()
        except Exception:
            pass
        measured = datetime.date.fromtimestamp(
            run_path.stat().st_mtime).isoformat()
        new = _bb.baseline_from_results(
            rows, measured_at=measured, commit=commit, note=note)
        old = None
        try:
            old = _bb.load_baseline(baseline_path)
        except ValueError as exc:
            print(f"note: existing reference unreadable ({exc})")
        print(f"Run:       {run_path}")
        print(f"Reference: {baseline_path}")
        if old is not None:
            print(f"  replacing {old.suite_pass_rate:.1%} over "
                  f"{old.total_samples} samples (measured {old.measured_at})")
        print(f"  with      {new.suite_pass_rate:.1%} over "
              f"{new.total_samples} samples (measured {measured})")
        if old is not None and new.suite_pass_rate < old.suite_pass_rate:
            print("  WARNING: this LOWERS the standard. A reference that "
                  "follows the code downwards guards nothing.")
        if not getattr(args, "yes", False):
            print("\nNothing written. Re-run with --yes to move the standard.")
            return 0
        _bb.save_baseline(new, baseline_path)
        print(f"\nWritten. Commit it so the change is reviewed: {baseline_path}")
        return 0

    # --- gate ---------------------------------------------------------------
    try:
        baseline = _bb.load_baseline(baseline_path)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    if baseline is None:
        print(f"No reference at {baseline_path}. Nothing to judge against; "
              "create one with `bench set-baseline`.", file=sys.stderr)
        return 0        # absence of a standard is not a failed run

    result = _bb.compare_to_baseline(rows, baseline)
    if getattr(args, "json", False):
        print(json.dumps(
            {**result, "run": str(run_path), "baseline": str(baseline_path)},
            indent=2, ensure_ascii=False))
    else:
        print(f"Run: {run_path}")
        print(_bb.format_baseline_report(result, baseline))
    return 1 if result.get("verdict") == "regressed" else 0


_BENCH_SCHEDULE_HINT = (
    "Scheduling is opt-in — nothing was registered automatically.\n"
    "To run this nightly via the scheduler daemon, run exactly:\n"
    "  delfin-agent scheduler add-bench --model {model}"
    "{provider}{backend} --every 24h\n"
    "  delfin-agent scheduler start\n"
    "Cost estimate: ~$8 per nightly run for the 48-task KIT-Qwen suite at "
    "repeats=1 (repeats multiply cost; recheck adds up to $3 more on "
    "suspect days)."
)


def _bench_schedule_hint(model: str, provider: str = "",
                         backend: str = "") -> str:
    return _BENCH_SCHEDULE_HINT.format(
        model=model or "<model>",
        provider=f" --provider {provider}" if provider else "",
        backend=f" --backend {backend}" if backend and backend != "api" else "",
    )


def _cmd_bench_nightly(args: argparse.Namespace) -> int:
    """Full unattended cycle: run suite → compare vs rolling history →
    recheck suspects → report → attention on confirmed regressions."""
    from . import bench_watch as _bw

    model = (getattr(args, "model", "") or "").strip()
    if not model:
        print("ERROR: --model is required for `bench nightly`",
              file=sys.stderr)
        return 2
    provider = getattr(args, "provider", "") or ""
    backend = getattr(args, "backend", "") or "api"
    repeats = max(1, int(getattr(args, "repeats", 1) or 1))
    recheck = not bool(getattr(args, "no_recheck", False))
    last_k = max(1, int(getattr(args, "last_k", _bw.DEFAULT_LAST_K)
                        or _bw.DEFAULT_LAST_K))
    runs_dir_arg = (getattr(args, "runs_dir", "") or "").strip() or None

    print(f"Nightly benchmark cycle for {model} "
          f"(repeats={repeats}, recheck={'on' if recheck else 'off'}, "
          f"baseline last-k={last_k})…", flush=True)
    summary = _bw.nightly(
        model, provider, backend,
        repeats=repeats, recheck=recheck, last_k=last_k,
        runs_dir=runs_dir_arg,
    )

    comparison = summary.get("comparison") or {}
    counts = comparison.get("counts") or {}
    if comparison:
        print(f"  compared {comparison.get('n_tasks', 0)} tasks: "
              f"{counts.get('stable', 0)} stable, "
              f"{counts.get('improved', 0)} improved, "
              f"{counts.get('suspect_regression', 0)} suspect, "
              f"{counts.get('new_task', 0)} new")
    rc = summary.get("recheck") or {}
    if rc:
        print(f"  recheck: {len(rc.get('confirmed') or [])} confirmed, "
              f"{len(rc.get('noise') or [])} noise, "
              f"{len(rc.get('skipped') or [])} skipped by cap "
              f"(spent ${float(rc.get('spent_usd') or 0):.2f})")
    for c in summary.get("confirmed") or []:
        print(f"  CONFIRMED regression: {c['task_id']} "
              f"(recheck q={float(c.get('recheck_quality') or 0):.0f})")
    if summary.get("attention_id"):
        print(f"  attention event raised: {summary['attention_id']}")
    for err in summary.get("errors") or []:
        print(f"  [error] {err}", file=sys.stderr)
    if summary.get("report_path"):
        print(f"\nReport: {summary['report_path']}")
    print()
    print(_bench_schedule_hint(model, provider, backend))
    return 0 if summary.get("ok") else 1


def _cmd_bench_compare_history(args: argparse.Namespace) -> int:
    """Classify one run vs the rolling baseline. No API spend."""
    from . import bench_watch as _bw
    from . import benchmark as _bm

    model = (getattr(args, "model", "") or "").strip()
    last_k = max(1, int(getattr(args, "last_k", _bw.DEFAULT_LAST_K)
                        or _bw.DEFAULT_LAST_K))
    runs_dir_arg = (getattr(args, "runs_dir", "") or "").strip() or None
    base_dir = (Path(runs_dir_arg).expanduser() if runs_dir_arg
                else _bm.runs_dir())

    run_arg = (getattr(args, "run", "") or "").strip()
    if run_arg:
        run_path = Path(run_arg).expanduser()
        if not run_path.exists():
            candidate = base_dir / run_arg
            if candidate.exists():
                run_path = candidate
            else:
                print(f"ERROR: run file not found: {run_arg}",
                      file=sys.stderr)
                return 2
    else:
        files = _bw.list_model_runs(model, runs_dir=runs_dir_arg)
        if not files:
            print(f"ERROR: no run files for model {model!r} in "
                  f"{base_dir}", file=sys.stderr)
            return 2
        run_path = files[0]

    # The compared run must never sit inside its own baseline window.
    history = _bw.load_history(
        model, runs_dir=runs_dir_arg, last_k=last_k, exclude=run_path)
    cmp = _bw.compare_run(run_path, history=history)
    if getattr(args, "json", False):
        print(json.dumps(cmp, indent=2, ensure_ascii=False))
        return 0
    counts = cmp.get("counts") or {}
    print(f"Run {Path(cmp['run_path']).name} vs rolling baseline "
          f"({len(history.get('files') or [])} file(s), last-k={last_k}):")
    print(f"  {cmp.get('n_tasks', 0)} tasks — "
          f"{counts.get('stable', 0)} stable, "
          f"{counts.get('improved', 0)} improved, "
          f"{counts.get('suspect_regression', 0)} suspect, "
          f"{counts.get('new_task', 0)} new")
    for e in cmp.get("suspects") or []:
        base = e.get("baseline") or {}
        print(f"  SUSPECT {e['task_id']:<28} ({e.get('reason', '?')}) "
              f"q={e.get('quality', 0):.0f} vs median "
              f"{float(base.get('median_quality') or 0):.0f} "
              f"(threshold {float(e.get('threshold') or 0):.1f})")
    if not cmp.get("suspects"):
        print("  no suspects — within the noise band of the baseline")
    print("\nCompare-only: nothing was re-run and no API cost was "
          "incurred. Use `bench nightly --model ...` for the full cycle.")
    return 0


def _print_behavior_rates(_bm, records) -> None:
    """Print per-behaviour parity rates for a run, if any task carried a
    ``behavior:`` tag.  No-op for suites without behaviour tasks."""
    rates = _bm.behavior_rates(records)
    if not rates:
        return
    print("\nBehaviour rates (planning / scouting / verifying / asking):")
    for flag in ("planned", "scouted", "verified", "asked"):
        r = rates.get(flag)
        if r:
            print(f"  {flag:<9} {r['rate']:>4.0%}  (n={r['n']})")


def cmd_knowledge(args: argparse.Namespace) -> int:
    """``delfin-agent knowledge ...``: the /knowledge command. A write is
    confirmed with a code shown and read on /dev/tty, never on stdout."""
    from .knowledge_commands import run_cli
    words = list(getattr(args, "words", None) or ["list"])
    return run_cli(words)


def _cmd_bug_report(args: argparse.Namespace) -> int:
    """``bug report``: file a report for a session that already ended.

    Input: --session (``latest`` = newest session of --workspace, an id,
    or a unique id prefix), -m description, --yes. Output: the report
    directory, printed; surface ``cli``. Exit 0 written, 1 declined or
    failed, 2 no such session or no way to confirm.

    Without --yes the description and a y/N confirmation are read from
    the terminal after a preview of what goes in. With no terminal and no
    --yes nothing is written: a pipe cannot see the preview it would be
    agreeing to. The copy to the remote archive follows
    ``agent.bug_upload`` (``bug_filing.upload``): asked separately without
    --yes, never made under --yes unless the setting is ``always``.
    """
    from . import bug_filing
    from .bug_report import resolve_archive_dir
    workspace = Path(getattr(args, "workspace", "") or os.getcwd()
                     ).expanduser().resolve()
    ref = str(getattr(args, "session", "") or "latest")
    data = bug_filing.resolve_session(ref, workspace)
    if data is None:
        where = (f"in {workspace}" if ref == "latest"
                 else f"matching {ref!r} (unknown or ambiguous)")
        print(f"ERROR: no saved session {where}", file=sys.stderr)
        return 2
    assume_yes = bool(getattr(args, "yes", False))
    interactive = False
    try:
        interactive = sys.stdin.isatty()
    except Exception:
        interactive = False
    if not assume_yes and not interactive:
        print("ERROR: confirmation needs a terminal; pass --yes to file "
              "without one", file=sys.stderr)
        return 2
    description = str(getattr(args, "message", "") or "").strip()
    if not description and not assume_yes:
        try:
            description = input("Describe the problem in one line: ").strip()
        except (EOFError, KeyboardInterrupt):
            description = ""
        if not description:
            print("no description given; nothing was written")
            return 1
    collected = bug_filing.collect_saved(data)
    print(bug_filing.preview(collected, resolve_archive_dir()))
    if not assume_yes:
        try:
            answer = input("File this report? [y/N] ")
        except (EOFError, KeyboardInterrupt):
            answer = ""
        if answer.strip().lower() not in ("y", "yes"):
            print("nothing was written")
            return 1
    try:
        report_dir = bug_filing.file_report(
            collected, description=description, surface="cli")
    except Exception as exc:
        print(f"ERROR: bug report failed: {exc}", file=sys.stderr)
        return 1
    print(f"bug report written: {report_dir}")
    # --yes answers the questions about writing, not the one about
    # sending: with it nobody is asked, so the copy to the remote archive
    # happens only under agent.bug_upload: always.
    remote = bug_filing.upload(report_dir,
                               ask=None if assume_yes else input)
    if remote:
        print(remote)
    return 0


def cmd_bug(args: argparse.Namespace) -> int:
    """Triage / watch the user bug-report archive (auto-analyse → Solved)."""
    from . import bug_watcher as _bw
    from . import bug_report as _br
    action = getattr(args, "bug_action", "") or "ls"
    archive_arg = getattr(args, "archive", "") or ""
    archive = (Path(archive_arg).expanduser() if archive_arg
               else _bw.resolve_watch_archive())

    if action == "report":
        return _cmd_bug_report(args)

    if action == "ls":
        pending = _bw.find_unsolved(archive)
        if not pending:
            print(f"No un-triaged reports in {archive}")
            return 0
        for d in pending:
            try:
                rep = _br.load_report(d)
            except Exception:
                rep = {}
            desc = str(rep.get("description") or "").replace(chr(10), " ")[:48]
            print(f"  {d.name:<44} {_br.report_surface(rep):<9} "
                  f"{str(rep.get('model', '?')):<22} {desc}")
        print(f"\n{len(pending)} pending in {archive}")
        return 0

    if action == "triage":
        name = getattr(args, "name", "") or ""
        if name:
            target = archive / name
            if not (target / "report.json").exists():
                print(f"ERROR: no report.json in {target}", file=sys.stderr)
                return 2
            targets = [target]
        else:
            targets = _bw.find_unsolved(archive)
        if not targets:
            print(f"No reports to triage in {archive}")
            return 0
        print(_bw.send_line())
        consent = _bw.ask_consent(len(targets))
        for d in targets:
            print(f"→ analysing {d.name} …", flush=True)
            t = _bw.triage_report(d, consent=consent)
            if t.withheld:
                print(f"  [WITHHELD] ({t.withheld_reason})")
            elif t.error:
                print(f"  [ERROR] {t.error}")
            else:
                print(f"  [OK] {t.summary[:200].strip()}")
                print(f"       triage → {d / 'triage.md'}"
                      + ("  · fix → fix_proposal.patch" if t.has_fix else ""))
        print("\nReports stay put (not solved yet). After you apply + verify "
              "a fix, run `bug solve <name>` to move it to Solved.")
        return 0

    if action == "solve":
        name = getattr(args, "name", "") or ""
        target = archive / name
        if not target.is_dir():
            print(f"ERROR: no such report dir: {target}", file=sys.stderr)
            return 2
        dest = _bw.mark_solved(target, archive=archive)
        if dest:
            print(f"Moved → {dest}")
            return 0
        print(f"ERROR: could not move {target} (permissions?)", file=sys.stderr)
        return 1

    if action == "watch":
        # Build an optional settings override for --archive / --interval.
        settings = None
        interval = int(getattr(args, "interval", 0) or 0)
        if archive_arg or interval:
            try:
                from delfin.user_settings import load_settings
                settings = load_settings() or {}
            except Exception:
                settings = {}
            bw_cfg = settings.setdefault("agent", {}).setdefault("bug_watcher", {})
            if archive_arg:
                bw_cfg["archive_dir"] = str(archive)
            if interval:
                bw_cfg["interval_s"] = interval

        if getattr(args, "once", False):
            print(f"Triaging {archive} once …", flush=True)
            print(_bw.send_line(settings))
            consent = _bw.ask_consent(
                len(_bw.find_unsolved(archive, settings=settings)), settings)
            results = _bw.run_once(archive=archive, settings=settings,
                                   consent=consent)
            # One line per report from fixed fields only; the analysis
            # itself is in triage.md (bug triage prints its first line).
            for t in results:
                if t.withheld:
                    tag, detail = "WITHHELD", f"({t.withheld_reason})"
                elif t.error:
                    tag, detail = "ERROR", t.error[:150]
                else:
                    tag = "OK"
                    detail = f"triage → {Path(t.report_path) / 'triage.md'}"
                print(f"  [{tag}] {t.report_name}: {detail}")
            print(f"Done — {len(results)} report(s) triaged.")
            return 0

        cfg = _bw.watcher_settings(settings)
        if not cfg["enabled"]:
            print("bug_watcher is disabled (agent.bug_watcher.enabled=false). "
                  "Enable it in Settings, or use `bug watch --once` / "
                  "`bug triage`.")
            return 2
        if not _bw.acquire_pid_lock():
            print("bug_watcher is already running (PID lock).")
            return 3
        from . import lifeline as _lifeline
        _lifeline.guard_daemon()
        try:
            print(f"bug_watcher started (interval {cfg['interval_s']}s, "
                  f"archive={archive}).")
            print(_bw.send_line(settings))
            return _bw.run_loop(settings=settings)
        finally:
            _bw.release_pid_lock()

    return 0


def cmd_credentials(args: argparse.Namespace) -> int:
    """Manage API keys: the OS keyring or ~/.delfin/credentials.json (0600).

    Where `set` writes is agent.credentials.store (auto|keyring|file).
    No subcommand here echoes a stored value — `list` masks every key and
    names its store, so the user can verify WHICH credentials are stored
    without exposing them. `migrate` moves file entries into the keyring
    and removes each from the file only after an equal read-back. Input is
    read via getpass so the value is never visible on screen and never
    lands in shell history.
    """
    import getpass
    from . import credentials as _cred
    action = getattr(args, "cred_action", "list") or "list"

    if action == "list":
        items = _cred.list_credentials()
        if not items:
            print("No credentials configured.")
            print()
            print("To store one securely (input is hidden, never echoed):")
            print("  delfin-agent credentials set "
                  "KIT_TOOLBOX_API_KEY")
            print("Other well-known names: OPENAI_API_KEY, ANTHROPIC_API_KEY")
            return 0
        mode = _cred.configured_store()
        status = _cred.keyring_status() if mode != "file" else None
        if status is None:
            kr = "not used (agent.credentials.store = file)"
        elif status["usable"]:
            kr = status["backend"]
        else:
            kr = f"unavailable: {status['reason']}"
        print(f"Credentials (store setting: {mode})")
        print(f"  file:    {_cred.credentials_path()}")
        print(f"  keyring: {kr}")
        print()
        for name in sorted(items):
            info = items[name]
            src = info.get("source", "?")
            tag = f"[{src}]"
            extra = [s for s in info.get("stores") or [] if s != src]
            also = f"  (also in: {', '.join(extra)})" if extra else ""
            print(f"  {name:<28}  {info.get('value', ''):<14}  {tag}{also}")
        return 0

    if action == "migrate":
        names = list(getattr(args, "names", None) or [])
        try:
            rows = _cred.migrate_to_keyring(names)
        except _cred.CredentialStoreUnavailable as exc:
            print(f"Nothing moved: {exc}", file=sys.stderr)
            return 1
        if not rows:
            print("No file entries to move.")
            return 0
        bad = 0
        for row in rows:
            reason = f" ({row['reason']})" if row.get("reason") else ""
            print(f"  {row['name']:<28}  {row['action']}{reason}")
            bad += row["action"] in ("failed", "differs")
        return 1 if bad else 0

    name = (getattr(args, "name", "") or "").strip()
    if not name:
        print("ERROR: credential name required.", file=sys.stderr)
        return 2

    if action == "delete":
        try:
            removed = _cred.delete_credential(name)
        except _cred.CredentialStoreUnavailable as exc:
            print(f"{name} is still in the keyring: {exc}", file=sys.stderr)
            return 1
        if removed:
            print(f"Removed {name}.")
            return 0
        print(f"No credential named {name}.", file=sys.stderr)
        return 1

    if action == "set":
        try:
            value = getpass.getpass(
                f"Enter value for {name} (input hidden, no echo): "
            )
        except (KeyboardInterrupt, EOFError):
            print()
            return 130
        if not value:
            print("No value entered, aborting.", file=sys.stderr)
            return 1
        try:
            ok = _cred.set_credential(name, value)
        except (_cred.CredentialStoreNotPrivate,
                _cred.CredentialStoreUnavailable) as exc:
            print(f"Not stored: {exc}", file=sys.stderr)
            return 1
        if ok:
            # No part of the value, not even a masked one: the first and
            # last four characters of a key stay in the terminal's
            # scrollback (alert 58). `credentials list` shows which key it is.
            print(f"Stored {name} in {_cred.describe_store(name)}")
            if _cred.stored_in(name) == ["keyring", "file"]:
                print(f"An older copy of {name} is still in "
                      f"{_cred.credentials_path()}; `delfin-agent credentials "
                      "migrate` or `delete` removes it.")
        else:
            print(f"{name} already stored with this value (no change).")
        return 0

    print(f"Unknown credentials action: {action}", file=sys.stderr)
    return 2


def cmd_skills(args: argparse.Namespace) -> int:
    """Skill proposals: what is waiting, its full preview, the decision.

    The other side of ``skill_proposals``: `show` prints text, evidence
    and safety findings so the accept/reject decision is made with a
    real preview, `accept` activates a pending proposal, `reject
    --reason` turns it away (into ``rejected/``, never deleted).
    """
    from . import skill_proposals as sp

    action = getattr(args, "skills_action", "") or "proposals"

    if action == "proposals":
        rows = sp.list_proposals()
        if not rows:
            print("no skill proposals waiting — a proposal is created at "
                  "the end of a qualifying session and accepted by a "
                  "person")
            return 0
        for p in rows:
            print(f"{p.name:<24} [{p.status}] {len(p.evidence)} evidence")
        print("\nshow one with: delfin-agent skills show <name>")
        return 0

    name = getattr(args, "name", "")
    if action == "show":
        p = sp.get_proposal(name)
        if p is None:
            print(f"ERROR: no proposal named {name!r}", file=sys.stderr)
            return 2
        print(f"proposal {p.name} — status: {p.status}")
        print(f"source: {p.source}   created: {p.created}")
        if p.base_version:
            print(f"base version: {p.base_version}")
        print()
        print(p.text)
        print()
        print(f"evidence ({len(p.evidence)}):")
        for ev in p.evidence:
            line = f"  [{ev.kind}] {ev.ref}"
            if ev.detail:
                line += f" — {ev.detail}"
            print(line)
        if p.findings:
            print("safety findings:")
            for f in p.findings:
                print(f"  ! {f}")
        else:
            print("safety findings: none")
        return 0

    # Recorded as who accepted: the account, not $USER (os_account).
    from .os_account import name as _account
    who = _account()

    if action == "accept":
        try:
            target = sp.accept(name, by=who)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        print(f"accepted: {name} -> {target}")
        print("the skill is live for new sessions now")
        return 0

    if action == "reject":
        reason = (getattr(args, "reason", "") or "").strip()
        if not reason:
            print("ERROR: a rejection needs --reason: the next session "
                  "reads why, instead of guessing the same thing again",
                  file=sys.stderr)
            return 2
        try:
            target = sp.reject(name, reason=reason, by=who)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        print(f"rejected: {name} -> {target}")
        return 0

    print(f"ERROR: unknown skills action {action!r}", file=sys.stderr)
    return 2


def _code_file_line(request_id: str) -> str:
    """"code ABC-DEF ..." for a run with no terminal whose code file this
    person can read, else "". Never under the agent-session marker: the
    agent's shell is refused the file, and this is the same refusal."""
    from . import approval_code as _code
    from .approval_guard import inside_agent_session
    if inside_agent_session():
        return ""
    code = _code.read_code_file(request_id)
    if not code:
        return ""
    return (f"{_code.display(code)} (the session has no terminal; "
            "approve/deny read it, no --code needed)")


def _approval_code_for(args: argparse.Namespace):
    """The one-time code for answering ``args.request_id`` from here.

    Output: the normalised code from ``--code`` or typed at this
    terminal; None when no code is asked for (nothing waits under the id,
    the request predates codes or is exempt, or this very process asked
    it); False, with the refusal printed and recorded, when one is needed
    and none was given. The code is checked by the waiting process, not
    here: this process does not know it.
    """
    from . import approval_code as _code
    from . import file_confirm as _fc
    from . import terminal_confirm as _tc
    rid = str(getattr(args, "request_id", "") or "")
    given = _code.normalize(getattr(args, "code", "") or "")
    if given:
        return given
    record = None
    for row in list(_fc.pending()) + list(_tc.pending_for_supervisors()):
        if str(row.get("id") or "") == rid:
            record = row
            break
    if record is None or not record.get("needs_code") or _code.code_for(rid):
        return None
    # A run with no terminal left its code in the deny-read code file.
    # Only the person reaches here: cmd_approvals refused the marker first.
    from .approval_guard import inside_agent_session
    if not inside_agent_session():
        from_file = _code.read_code_file(rid)
        if from_file:
            return from_file
    try:
        if sys.stdin.isatty():
            typed = _code.normalize(input(
                f"code shown by the session waiting on {rid}: "))
            if typed:
                return typed
    except (EOFError, KeyboardInterrupt, OSError):
        pass
    try:
        from .api_client import _record_security_event
        _record_security_event(
            "approval_code_missing", "approvals",
            f"{rid}: an answer without the session's code was not written",
            blocked=True)
    except Exception:
        pass
    print(f"ERROR: answering {rid} needs the one-time code the waiting "
          "session shows (its terminal, or its row in the dashboard); "
          "pass it with --code", file=sys.stderr)
    return False


def cmd_approvals(args: argparse.Namespace) -> int:
    """The questions a headless session waits on, and the answers.

    A session with no terminal has nobody to ask, so the Self-Modification
    Guard refuses outright and the most valuable work is the work it
    cannot do. This is the other side of ``file_confirm``: what is
    waiting, what it wants to change, and the two words that let it
    through or turn it away. ``watch`` prints one line per new question,
    which is what a supervising process reads.
    """
    import time as _time

    from . import file_confirm as _fc

    action = getattr(args, "approvals_action", "") or "ls"

    if action in ("approve", "deny", "answer"):
        # Run by an agent's shell (the marker every agent child carries):
        # the person answers, not the session that waits nor another one.
        # The gates refuse the command before it runs; this holds where a
        # spelling got past them. ls / show / watch stay readable.
        from . import approval_guard as _guard
        if _guard.inside_agent_session():
            try:
                from .api_client import _record_security_event
                _record_security_event(
                    "approval_self_answer", "approvals",
                    f"{action} {args.request_id} from an agent session",
                    blocked=True)
            except Exception:
                pass
            print(_guard.refusal(
                "`delfin-agent approvals "
                f"{action}` was run by an agent session"), file=sys.stderr)
            return 2
        code = _approval_code_for(args)
        if code is False:
            return 2
    else:
        code = None

    if action == "answer":
        # A choice question (ask_user_question) published by a terminal
        # session. approve/deny is no answer to one: the operator had to
        # deny and hide the pick in the refusal reason. The choice itself
        # resolves against the question's own options; anything that is
        # not a choice question is refused here, not guessed at.
        from . import approval_answers as _ans
        from .os_account import name as _account
        who = _account()
        try:
            picks = _ans.answer(args.request_id, args.choice, by=who,
                                code=code)
        except _ans.ChoiceError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        print(f"answer: {args.request_id}")
        print(f"choice: {', '.join(picks)}")
        return 0

    if action in ("approve", "deny"):
        from .os_account import name as _account
        who = _account()
        # A refusal may say what to do instead, and the model reads it in
        # the same turn -- instead of guessing the same thing spelled
        # differently, which is one more dialog for the person here.
        reason = (getattr(args, "reason", "") or "") if action == "deny" else ""
        ok = _fc.answer(args.request_id, action == "approve", by=who,
                        reason=reason, code=code)
        if not ok:
            # The same id may belong to a session waiting at a terminal.
            # Answering it here is the point of publishing it whole: the
            # pane shows 24 cut lines, this shows the question.
            from . import terminal_confirm as _tc
            ok = _tc.answer_waiting(args.request_id,
                                    action == "approve", by=who,
                                    reason=reason, code=code)
        if not ok:
            print(f"ERROR: nothing waiting with id {args.request_id!r}",
                  file=sys.stderr)
            return 2
        print(f"{action}: {args.request_id}")
        if reason:
            print(f"reason: {reason}")
        return 0

    if action == "show":
        from . import approval_answers as _ansq
        for row in _fc.pending():
            if row.get("id") == args.request_id:
                print(f"id      {row['id']}")
                print(f"session {row.get('session_id', '')}")
                print(f"tool    {row.get('tool', '')}")
                print(f"path    {row.get('path', '')}")
                if row.get("needs_code"):
                    print("code    " + (
                        _code_file_line(row["id"]) or
                        "shown on that session's terminal; pass it with "
                        "--code"))
                print(f"waited  {int(_time.time() - float(row.get('asked_at') or 0))}s"
                      f" of {int(row.get('timeout_s') or 0)}s")
                print()
                body = _ansq.render_question(row)
                print(body or row.get("preview", ""))
                return 0
        from . import terminal_confirm as _tc
        # Every surface, not only the terminals: a dashboard window that
        # was reloaded cannot answer its own request any more, and this
        # is where it is answered instead.
        for row in _tc.pending_for_supervisors():
            if row.get("id") == args.request_id:
                where = ("that session's terminal"
                         if _tc.surface_of(row) == _tc.TERMINAL
                         else "that session's dashboard window")
                print(f"id      {row['id']}")
                print(f"session {row.get('session_key') or row.get('session_id', '')}")
                print(f"tool    {row.get('tool', '')}")
                print(f"command {row.get('command', '')}")
                print(f"waited  {int(_time.time() - float(row.get('asked_at') or 0))}s")
                if row.get("needs_code"):
                    print(f"code    shown at {where}; pass it with --code")
                print(f"answer  here or at {where} — "
                      "whichever comes first")
                print()
                body = _ansq.render_question(row)
                print(body or row.get("preview", ""))
                return 0
        print(f"ERROR: nothing waiting with id {args.request_id!r}",
              file=sys.stderr)
        return 2

    if action == "watch":
        seen: set = set()
        limit = float(getattr(args, "seconds", 0) or 0)
        deadline = _time.monotonic() + (limit if limit > 0 else 1e9)
        while _time.monotonic() < deadline:
            for row in _fc.pending():
                rid = row.get("id")
                if rid in seen:
                    continue
                seen.add(rid)
                print(f"{rid}  {row.get('tool', '')}  "
                      f"{row.get('path') or '-'}  "
                      f"session={row.get('session_id', '')}", flush=True)
            _time.sleep(1.0)
        return 0

    rows = _fc.pending()
    # A question a supervisor cannot read is one they cannot answer: `ls`
    # showed id/tool/path and seconds, so every question needed a second
    # `show`. One line per entry now says what it is about -- command,
    # path, or the question of a choice request -- clipped to the pane
    # width so a long command cannot push the rest off the screen.
    from . import terminal_confirm as _tc
    from . import approval_answers as _ans
    at_terminals = _tc.pending_at_terminals()
    in_dashboards = _tc.pending_in_dashboards()
    if not rows and not at_terminals and not in_dashboards:
        print("(nothing waiting)")
        return 0
    width = _ls_width()
    for row in rows:
        waited = int(_time.time() - float(row.get("asked_at") or 0))
        head = f"{row.get('id', ''):<24} {row.get('tool', ''):<12} "
        tail = f" {waited}s"
        subject = _ans.preview_line(
            row, max(10, width - len(head) - len(tail)))
        print((head + subject + tail)[:width])
        code_line = _code_file_line(str(row.get("id") or ""))
        if code_line:
            print(f"{'':<25}code {code_line}")
    def _group(heading: str, group_rows: list) -> None:
        if not group_rows:
            return
        print(heading)
        for row in group_rows:
            waited = int(_time.time() - float(row.get("asked_at") or 0))
            mark = " PROTECTED" if row.get("protected") else ""
            head = (f"  {row.get('id', ''):<24} "
                    f"{(row.get('session_key') or '?'):<14} "
                    f"{row.get('tool', ''):<12} ")
            tail = f" {waited}s{mark}"
            subject = _ans.preview_line(
                row, max(10, width - len(head) - len(tail)))
            print((head + subject + tail)[:width])

    _group("\nWaiting at a terminal — answer here or in that session:",
           at_terminals)
    # Its own heading because it is its own fact: a dashboard window can
    # be gone while its request stands, and then here is the only place
    # left to answer it.
    _group("\nWaiting in a dashboard — answer here or in that window:",
           in_dashboards)
    return 0


def _ls_width() -> int:
    """The pane width for `approvals ls`, 120 when there is none.

    COLUMNS is respected first -- inside a pipe it carries the width of
    the terminal the operator is actually looking at, which is the one
    the clipping is for.
    """
    try:
        cols = int(os.environ.get("COLUMNS", "") or 0)
        if cols > 0:
            return cols
    except ValueError:
        pass
    try:
        import shutil
        return shutil.get_terminal_size().columns or 120
    except Exception:
        return 120


def cmd_memory(args: argparse.Namespace) -> int:
    """Show what tidying the memory store would do, and do it on request.

    Merging near-duplicates happens today only as a memory is WRITTEN, so
    look-alikes already in a store stay there; and the decay prune
    deletes on its own schedule without showing anyone what it took. This
    reports first and changes nothing until --apply.
    """
    from pathlib import Path as _P

    from . import memory_tidy as _tidy
    from .memory_store import _delfin_memory_dir

    root = _P(getattr(args, "workspace", "") or os.getcwd())
    store = _delfin_memory_dir(root)
    proposal = _tidy.propose(store, repo_root=root)
    print(f"store: {store}")
    print(proposal.render())
    if not getattr(args, "apply", False):
        return 0
    if not (proposal.merges or proposal.retire):
        return 0
    done = _tidy.apply(proposal)
    print(f"merged {done['merged']}, retired {done['retired']}"
          + (f", failed {done['failed']}" if done["failed"] else ""))
    print(f"retired files are in {store / 'retired'} — nothing was deleted")
    return 0


def cmd_messages(args: argparse.Namespace) -> int:
    """Delivery receipts for one sender, from the terminal.

    Both subcommands are scoped to the caller's own ``--from`` key
    (``status(message_id, from_key)`` / ``ls(from_key)``), the security
    boundary the operator set: a caller that names no key, or a key that is
    not the sender, sees nothing -- ``unknown`` covers both 'no such message'
    and 'not yours', so a foreign message's existence is never revealed.
    Rows and receipts carry no message body, so the CLI never prints one.
    """
    from . import session_messages as _msgs
    from_key = getattr(args, "from_key", "") or ""
    if args.messages_action == "ls":
        rows = _msgs.ls(from_key)
        if not rows:
            print("(nothing to see)" if not from_key else "(no messages)")
            return 0
        for r in rows:
            print(f"{r['id'][:12]:<14} "
                  f"{str(r.get('to') or ''):<12} {r['status']}")
        return 0
    if args.messages_action == "status":
        mid = getattr(args, "message_id", "") or ""
        if not mid:
            print("ERROR: message_id is required", file=sys.stderr)
            return 2
        print(_msgs.status(mid, from_key))
        return 0
    print("ERROR: unknown messages action", file=sys.stderr)
    return 2


def cmd_session(args: argparse.Namespace) -> int:
    from . import session_store as _ss
    if args.session_action == "ls":
        rows = _ss.list_sessions(limit=args.limit or 20)
        if not rows:
            print("(no saved sessions)")
            return 0
        import time
        for r in rows:
            when = time.strftime("%Y-%m-%d %H:%M",
                                  time.localtime(r.get("updated_at", 0)))
            print(f"{when}  {r['session_id'][:16]:<18}  "
                  f"{r.get('title','')[:60]}")
        return 0
    if args.session_action == "search":
        if not args.query:
            print("ERROR: query is required", file=sys.stderr)
            return 2
        q = (args.query or "").lower()
        hits = []
        for r in _ss.list_sessions(limit=100):
            data = _ss.load_session(r["session_id"]) or {}
            for i, m in enumerate(data.get("chat_messages") or []):
                if q in str(m.get("content", "")).lower():
                    hits.append((r["session_id"], i, m.get("role", "?"),
                                 str(m.get("content", ""))[:120]))
        for sid, i, role, snippet in hits[:30]:
            print(f"  {sid[:12]}  msg#{i:<3}  {role:<10}  "
                  f"{snippet.replace(chr(10), ' ')}")
        if not hits:
            print(f"(no matches for {args.query!r})")
        return 0
    if args.session_action == "export":
        from . import session_export as _se
        sid = args.session_id
        if sid == "latest":
            row = _ss.latest_session()
            if not row:
                print("ERROR: no saved sessions to export", file=sys.stderr)
                return 1
            sid = row["session_id"]
        data = _ss.load_session(sid)
        if data is None:
            print(f"ERROR: session {sid!r} not found", file=sys.stderr)
            return 1
        out = args.notebook or f"{sid}.ipynb"
        path = _se.export_session_to_file(data, out)
        print(f"exported {sid} -> {path}")
        return 0
    print(f"ERROR: unknown session action {args.session_action!r}",
          file=sys.stderr)
    return 2


def _print_stop_all_check(report: dict) -> int:
    """Print ``stop_all.status()``. Returns 0 when nothing is running, 1
    when something is."""
    import time as _time

    def _clock(epoch) -> str:
        try:
            return _time.strftime("%m-%d %H:%M", _time.localtime(float(epoch)))
        except (TypeError, ValueError):
            return "?"

    from .stop_all import HEARTBEAT_FRESH_S

    host = report.get("host") or "?"
    print(f"DELFIN agents -- check only, nothing is changed "
          f"({host}, {_time.strftime('%H:%M:%S', _time.localtime(report['at']))})")
    print("Every login node (read from the shared home directory):")
    sessions = report.get("open_sessions") or []
    fresh = [r for r in sessions
             if r["seconds_since_heartbeat"] <= HEARTBEAT_FRESH_S]
    print(f"  Open dashboard sessions: "
          + (f"{len(fresh)} with a heartbeat in the last "
             f"{HEARTBEAT_FRESH_S // 60} min" if fresh else "none")
          + (f", {len(sessions) - len(fresh)} silent (ended)"
             if len(sessions) > len(fresh) else ""))
    for row in sessions:
        age = row["seconds_since_heartbeat"]
        print(f"    {row['host']}: {row['title']!r} "
              + (f"(heartbeat {age} s ago)" if age <= HEARTBEAT_FRESH_S
                 else f"(last heartbeat {age // 60} min ago -- ended)"))
    kept = report.get("kept_sessions") or []
    print(f"  Kept sessions:           {len(kept) or 'none'}")
    for row in kept:
        where = ("running on this machine" if row.get("alive_here")
                 else "its process is gone -- left over" if row.get("here")
                 else f"run this check on {row['host']} to be certain")
        print(f"    {row['session']} on {row['host']} since "
              f"{_clock(row['since'])} ({where})")
    schedules = report.get("active_schedules") or []
    print(f"  Active schedules:        {len(schedules) or 'none'}"
          + (f" (could not read them: {report['schedules_error']})"
             if report.get("schedules_error") else ""))
    for row in schedules:
        print(f"    {row['id']} next {_clock(row['next_fire_at'])}: {row['prompt']}")
    daemons = report.get("daemon_pid_files") or {}
    print(f"  Daemons (pid files):     {len(daemons) or 'none'}")
    for name, row in daemons.items():
        print(f"    {name}: pid {row['pid']} ("
              + ("running on this machine" if row.get("alive_here") else
                 "not running on this machine; left over unless it runs on "
                 "another login node") + ")")
    stop = report.get("last_stop") or {}
    print("  Last emergency stop:     "
          + (f"{_clock(stop.get('at'))} on {stop.get('host')}" if stop else "never"))
    print(f"This machine only ({host}) -- run this on every login node you used:")
    here = report.get("agent_processes_here") or []
    print(f"  Agent processes:         {len(here) or 'none'}")
    for row in here:
        print(f"    pid {row['pid']}: {row['command']}")
    outlived = report.get("outlived_their_start_here") or []
    print(f"  Outlived their start:    {len(outlived) or 'none'}")
    for row in outlived:
        print(f"    pid {row['pid']} ({row['why']}): {row['command']}")
    running = bool(
        fresh or schedules or here or outlived
        or any(r.get("alive_here") for r in kept)
        or any(r.get("alive_here") for r in daemons.values()))
    if running:
        print("Result: something is running or scheduled -- see above. "
              "To end all of it: delfin-agent stop-all")
        return 1
    leftovers = len(sessions) + len(kept) + len(daemons)
    if leftovers:
        print("Result: nothing runs on this machine, nothing is scheduled, "
              f"and no session has sent a heartbeat in the last "
              f"{HEARTBEAT_FRESH_S // 60} min. The records above are left "
              "over; where one names another login node, run this check "
              "there too.")
        return 0
    print("Result: no agent is running or scheduled.")
    return 0


def _refuse_pause_inside_agent(verb: str) -> bool:
    """True (and says so) when an agent session started this process: a
    pause or a resume decides over another session, which is the
    person's, as for ``ops pause`` (ops_interface._refuse_inside_agent)."""
    from . import approval_guard
    if not approval_guard.inside_agent_session():
        return False
    print(f"ERROR: {verb} is run by the person, not from inside an agent "
          "session", file=sys.stderr)
    return True


def cmd_pause(args: argparse.Namespace) -> int:
    """Pause one session (see ``session_pause``). It takes effect at the
    session's next tool call; nothing else is touched. Refused inside an
    agent session."""
    if _refuse_pause_inside_agent("delfin-agent pause"):
        return 2
    from . import session_pause
    session_pause.pause(args.session, args.reason or "delfin-agent pause")
    print(f"paused {args.session}: it stops at its next tool call; "
          f"resume with: delfin-agent resume {args.session}")
    return 0


def cmd_mcp_pins(args: argparse.Namespace) -> int:
    """The MCP tool descriptions pinned on first use, and the changes
    waiting for the person (mcp_client). ``accept`` is the person's:
    refused inside an agent session."""
    from . import mcp_client as _mcp
    if getattr(args, "pins_action", "") == "accept":
        if _refuse_pause_inside_agent("delfin-agent mcp-pins accept"):
            return 2
        done = _mcp.accept_pending(args.server, args.tool or "")
        if not done:
            print(f"nothing waits for {args.server}"
                  + (f"/{args.tool}" if args.tool else ""))
            return 1
        for name in done:
            print(f"accepted {name}: the next session shows it to the model")
        return 0
    from .session_messages import _plain
    rows = _mcp.pins_overview()
    if not rows:
        print("no MCP server pinned yet")
        return 0
    for row in rows:
        print(f"{_plain(row['name'])}  ({len(row['tools'])} tools pinned)")
        for tool, info in sorted(row["pending"].items()):
            state = "new tool" if info.get("state") == "new" else "changed"
            print(f"  {state}: {_plain(tool)} -- "
                  f"{_plain(info.get('description') or '')[:200]}")
            print(f"    accept with: delfin-agent mcp-pins accept "
                  f"{_plain(row['name'])} {_plain(tool)}")
    return 0


def cmd_resume(args: argparse.Namespace) -> int:
    """Lift a pause set by ``delfin-agent pause``. Refused inside an agent
    session."""
    if _refuse_pause_inside_agent("delfin-agent resume"):
        return 2
    from . import session_pause
    was = session_pause.resume(args.session)
    print(f"resumed {args.session}" if was else f"{args.session} was not paused")
    return 0


def cmd_stop_all(args: argparse.Namespace) -> int:
    """The emergency stop: end every agent of this user on every login node.

    See ``stop_all``. Asks first unless ``--yes``: it ends dashboards'
    kernels, not only the agent in them.
    """
    import time as _time

    from . import session_presence, stop_all

    if args.check:
        return _print_stop_all_check(stop_all.status())
    open_before = session_presence.open_sessions()
    if not args.yes:
        hosts = sorted({str(r.get("host") or "?").split(".")[0]
                        for r in open_before})
        print("This ends every DELFIN agent of yours on every login node: "
              "dashboard kernels (with their agents, shells and MCP "
              "servers) and the daemons. Cluster jobs keep running.")
        if hosts:
            print(f"Open sessions right now: {len(open_before)} "
                  f"on {', '.join(hosts)}.")
        try:
            answer = input("Stop everything? [y/N] ").strip().lower()
        except EOFError:
            answer = ""
        if answer not in ("y", "yes"):
            print("Nothing stopped.")
            return 1
    summary = stop_all.request(args.reason or "delfin-agent stop-all")
    when = _time.strftime("%H:%M:%S", _time.localtime(summary["at"]))
    print(f"Emergency stop given at {when} on {summary['host']}.")
    print("  Every kernel and daemon of yours that started before it ends "
          "within a few seconds, on every login node that shares this home.")
    found = summary.get("found_here") or []
    ended = summary.get("ended_here") or []
    print(f"  On this machine: {len(found)} agent process(es) found; "
          f"{len(found) - len(ended)} ended on their own, {len(ended)} had "
          "to be ended" + (f" ({', '.join(map(str, ended))})" if ended else "")
          + ".")
    print(f"  Schedules disabled: {summary.get('schedules_disabled', 0)}"
          + (f" (could not read them: {summary['schedules_error']})"
             if summary.get("schedules_error") else "") + ".")
    others: dict[str, int] = {}
    for record in open_before:
        host = str(record.get("host") or "?").split(".")[0]
        if host != summary["host"]:
            others[host] = others.get(host, 0) + 1
    if others:
        print("  Sessions on other machines, which end there on their own: "
              + ", ".join(f"{h}: {n}" for h, n in sorted(others.items())) + ".")
    print("  Nothing starts on its own again until you send a session "
          "something yourself. Cluster jobs keep running (squeue --me).")
    return 0


def cmd_scheduler(args: argparse.Namespace) -> int:
    """Control the headless scheduler daemon: start / status / stop.

    The daemon executes ``schedule_wakeup`` / ``cron_create`` entries
    (``~/.delfin/cron.json``) without an open dashboard. It only runs
    entries the user explicitly created — see
    :mod:`delfin.agent.scheduler_daemon` for the full contract.
    """
    from . import scheduler_daemon as _sd

    action = getattr(args, "scheduler_action", "") or "status"

    if action == "start":
        st = _sd.daemon_status()
        if st["running"]:
            print(f"Scheduler daemon already running (PID {st['pid']}).")
            return 0
        import subprocess as _sp
        log = Path.home() / ".delfin" / "scheduler_daemon.log"
        log.parent.mkdir(parents=True, exist_ok=True)
        from . import lifeline as _lifeline
        with log.open("a") as lf:
            daemon = _sp.Popen(
                [sys.executable, "-m", "delfin.agent.scheduler_daemon"],
                stdout=lf, stderr=lf,
                # Detached from Ctrl+C in the shell, but not from the
                # terminal: it follows this terminal (or the dashboard it
                # was started from) through the lifeline, and ends with it.
                start_new_session=True,
                env=_lifeline.child_env(),
            )
        _lifeline.record_child(int(getattr(daemon, "pid", 0) or 0),
                               "scheduler_daemon")
        print(f"Scheduler daemon started (detached); log: {log}")
        print("It ends when this terminal closes.")
        print("It executes ONLY entries you explicitly scheduled "
              "(schedule_wakeup / cron_create) — each fire is one paid "
              "agent turn.")
        return 0

    if action == "add-bench":
        # Explicit opt-in creation of a recurring benchmark entry — the
        # only way a nightly bench lands in the schedule. The entry is
        # executed by the daemon's bench hook (scheduler_daemon), which
        # calls bench_watch.nightly directly instead of an LLM turn.
        from . import scheduler as _sched
        from . import scheduler_daemon as _sdm

        model = (getattr(args, "model", "") or "").strip()
        if not model:
            print("ERROR: --model is required for `scheduler add-bench`",
                  file=sys.stderr)
            return 2
        every_token = getattr(args, "every", "") or "24h"
        every_s = _sched.parse_interval_seconds(every_token)
        if every_s is None:
            print(f"ERROR: bad --every interval {every_token!r} "
                  "(use e.g. 24h, 12h, 1d)", file=sys.stderr)
            return 2
        prompt = _sdm.format_bench_entry_prompt(
            model=model,
            provider=getattr(args, "provider", "") or "",
            backend=getattr(args, "backend", "") or "",
            repeats=max(1, int(getattr(args, "repeats", 1) or 1)),
        )
        ent = _sched.Scheduler().schedule_interval(
            every_seconds=every_s,
            prompt=prompt,
            reason=f"nightly benchmark watch ({model})",
            workspace=os.getcwd(),
        )
        print(f"Scheduled bench entry {ent.id}: every {every_token}, "
              f"workspace {ent.workspace}")
        print("Cost estimate: ~$8 per run for the 48-task KIT-Qwen suite "
              "at repeats=1 (repeats multiply cost; the suspect recheck "
              "adds up to $3 more on suspect days).")
        st = _sd.daemon_status()
        if not st["running"]:
            print("The scheduler daemon is NOT running — the entry fires "
                  "only after you start it:\n"
                  "  delfin-agent scheduler start")
        print(f"Remove anytime: delete entry {ent.id} via the dashboard "
              "scheduler tools, or edit ~/.delfin/cron.json.")
        return 0

    if action == "stop":
        st = _sd.daemon_status()
        if not st["running"]:
            print("No scheduler daemon running.")
            return 0
        import signal as _sig
        try:
            os.kill(st["pid"], _sig.SIGTERM)
        except OSError as exc:
            print(f"ERROR: stop failed: {exc}", file=sys.stderr)
            return 1
        print(f"SIGTERM sent to scheduler daemon (PID {st['pid']}); "
              "it exits after the current turn.")
        return 0

    # status (default)
    import time as _time
    from .scheduler import Scheduler
    st = _sd.daemon_status()
    state = f"running (PID {st['pid']})" if st["running"] else "not running"
    print(f"Scheduler daemon: {state}")
    entries = sorted(Scheduler().list_entries(), key=lambda e: e.next_fire_at)
    if not entries:
        print("No scheduled entries (~/.delfin/cron.json is empty).")
        return 0
    print(f"{len(entries)} entries ({st['disabled']} disabled):")
    now = _time.time()
    for e in entries:
        if e.disabled:
            due = f"DISABLED — {e.disabled_reason or 'no reason recorded'}"
        else:
            dt = int(e.next_fire_at - now)
            due = f"due in {dt}s" if dt > 0 else f"overdue by {-dt}s"
        prompt = (e.prompt or "").strip().replace("\n", " ")
        if len(prompt) > 60:
            prompt = prompt[:57] + "..."
        every = f", every {e.every_seconds}s" if e.kind == "interval" else ""
        ws = f" @ {e.workspace}" if e.workspace else ""
        print(f"  {e.id}  [{e.kind}{every}]  {due}{ws}  — {prompt}")
    return 0


def _doctor_repair(results: list) -> int:
    """The `--repair` ask-loop: one approval per fixable step, in order.

    Refuses to repair at all unless stdin is a terminal -- a piped or
    scripted call must never click through an approval the user did not
    give. Every step carries what it changes and how to undo it; the user
    approves each before it runs, and a declined step changes nothing.
    Returns the doctor's non-zero-exit convention for repair failures.
    """
    from . import repair as _repair

    steps = _repair.plan(results)
    if not steps:
        print("Nothing machine-actionable to repair "
              "(the failing checks need a human fix).")
        return 0
    if not sys.stdin.isatty():
        print("Repair needs a terminal to ask about each step; refusing. "
              "Run `doctor --repair` in a terminal to approve repairs one "
              "at a time.", file=sys.stderr)
        return 1
    failed = False
    for step in steps:
        print(f"[{step['id']}] {step['what']}")
        print(f"    undo: {step['undo']}")
        try:
            answer = input("Apply this repair step? [y/N] ").strip().lower()
        except (EOFError, KeyboardInterrupt):
            print(file=sys.stderr)
            print("Aborted; no further steps were applied.", file=sys.stderr)
            return 1
        if answer not in ("y", "yes"):
            print("Skipped.")
            continue
        outcome = _repair.apply(step, approved=True)
        if outcome.get("ok"):
            print(f"Applied ({outcome.get('status', 'ok')}); "
                  f"undo: {outcome.get('undo', '')}")
        else:
            failed = True
            status = outcome.get("status", "failed")
            print(f"Failed to apply: {status}", file=sys.stderr)
    return 1 if failed else 0


def cmd_audit(args: argparse.Namespace) -> int:
    """``delfin-agent audit verify``: 0 when every log checked is intact,
    1 when one is not, 2 without a subcommand. Reads only."""
    from . import audit_log as _al
    if getattr(args, "audit_action", "") != "verify":
        print("usage: delfin-agent audit verify [--log PATH] [--json]",
              file=sys.stderr)
        return 2
    log = str(getattr(args, "log", "") or "")
    reports = [_al.verify(Path(log).expanduser())] if log else _al.verify_all()
    if getattr(args, "json", False):
        print(json.dumps(reports, indent=2))
    else:
        print(_al.format_verify(reports))
    return 0 if all(r["ok"] for r in reports) else 1


def cmd_doctor(args: argparse.Namespace) -> int:
    """One-surface prerequisite report: docs index, credentials,
    binaries, Python deps, MCP servers, scheduler, attention inbox,
    benchmark ground truth, memory store, disk space.

    Exit code 1 when any check FAILs (warnings alone exit 0), so the
    command is scriptable as a pre-flight gate. With `--repair`, prompts
    in a terminal about each fixable step and applies those approved; the
    exit code is the repair's own (0 = nothing to repair or all steps
    applied/skipped, 1 = refusal or a failed step).
    """
    from . import doctor as _doc

    workspace = getattr(args, "workspace", "") or None
    results = _doc.run_doctor(workspace)
    if not getattr(args, "repair", False):
        print(_doc.format_doctor(results))
        return 1 if any(r.get("status") == "FAIL" for r in results) else 0
    print(_doc.format_doctor(results))
    return _doctor_repair(results)


def _load_limits() -> dict:
    """The user's load alarm thresholds, from `agent.load_limits` settings.

    ``~/.delfin/settings.json`` always, plus the trusted project files —
    the same sources hooks read. Anything unreadable or malformed falls
    back to session_load's defaults: thresholds are configuration, and
    broken configuration must not break the reporting.
    """
    try:
        from . import hooks as _hooks
        # The raw settings dict is needed for agent.*, which load_hooks'
        # merged HookConfig does not carry; read the user file directly.
        # The trusted-project override is a refinement a settings file of
        # this kind does not need yet.
        raw = _hooks._read_json_safe(_hooks._user_settings_path()) or {}
        limits = raw.get("agent", {}).get("load_limits")
        if isinstance(limits, dict):
            return dict(limits)
    except Exception:
        pass
    from . import session_load as _sl_default
    return dict(_sl_default.DEFAULT_LIMITS)


def format_load_row(key: str, load: dict, alarms: list[str]) -> str:
    """One watch line: session, procs, rss, cpu — alarms up front.

    cpu is a difference between two measurements; before the second one
    there is no number, and a dash says so instead of a made-up 0.0.
    """
    rss = float(load.get("rss_bytes") or 0)
    for unit, factor in (("G", 1 << 30), ("M", 1 << 20), ("k", 1 << 10)):
        if rss >= factor or unit == "k":
            rss_txt = f"{rss / factor:.1f}{unit}"
            break
    cores = load.get("cpu_cores")
    cores_txt = "-" if cores is None else f"{cores:.1f}"
    line = (f"  {key[:16]:<16} procs {int(load.get('procs') or 0):>4}  "
            f"rss {rss_txt:>6}  cpu {cores_txt:>4}")
    for alarm in alarms:
        line += f"  ! {alarm}"
    return line


def _print_session_load(watch_seconds: float) -> None:
    """The --load view: one line per open session on THIS host, repeated.

    Only this host's sessions have readable pids; others are listed
    without numbers so the absence is explained rather than silent.
    cpu needs two measurements a window apart, so the first pass shows
    a dash and the second onward shows cores. Never raises: a monitor
    that crashes on a dead session is worse than no monitor.
    """
    import socket as _socket
    import time as _time

    from . import session_load as _sl

    host = _socket.gethostname()
    limits = _load_limits()
    previous: dict[str, dict] = {}   # key -> last load_of, for cpu cores
    while True:
        try:
            from . import session_presence as _pres
            records = _pres.open_sessions()
        except Exception:
            records = []
        for rec in records:
            pid = int(rec.get("pid") or 0)
            key = str(rec.get("key") or "?")
            if str(rec.get("host") or host) != host or pid <= 0:
                print(f"  {key[:16]:<16} (on another host — not counted here)")
                continue
            load = _sl.load_of(pid, previous.get(key))
            print(format_load_row(key, load, _sl.alarms(load, limits)))
            previous[key] = load
        if not watch_seconds:
            return
        _time.sleep(watch_seconds)


def cmd_sessions(args: argparse.Namespace) -> int:
    """What ran here before, as a table that says how to come back."""
    from . import cli_resume as _cr

    watch = getattr(args, "watch", None)
    if getattr(args, "load", False) or watch:
        try:
            _print_session_load(float(watch) if watch else 0)
        except Exception as exc:  # a broken /proc must not hide history
            print(f"(load unavailable: {exc})")
        return 0

    # What is open NOW comes first. A table of history read as "nothing
    # is running" while five sessions were, and the supervisor had to
    # assemble the live picture from tmux panes and git logs instead.
    live_rows = _cr.open_now()
    live = _cr.render_open(live_rows)
    if live:
        print(live)
        # What the open ones cost this host, in the same breath: the
        # listing is where an operator looks first, so the load belongs
        # here too — one quiet line each, alarms and all.
        try:
            import socket as _socket
            from . import session_load as _sl
            limits = _load_limits()
            host = _socket.gethostname()
            for row in live_rows:
                pid = int(row.get("pid") or 0)
                key = str(row.get("key") or "?")
                if str(row.get("host") or host) != host or pid <= 0:
                    continue
                load = _sl.load_of(pid)
                print(format_load_row(key, load, _sl.alarms(load, limits)))
        except Exception:  # history answers even when /proc cannot
            pass
        print()
    rows = _cr.list_sessions(limit=max(1, int(getattr(args, "limit", 20) or 20)))
    if not rows:
        if not live:
            print("No sessions recorded yet.")
        return 0
    print(_cr.render_sessions(rows))
    return 0


def cmd_trace(args: argparse.Namespace) -> int:
    """One run read as a trace: turns, model calls, tool calls (OT1).

    Text timeline from the existing stores only (tool_trace, turn_metrics,
    session_store) — no new store, nothing sent anywhere.
    """
    from .otel_export import spans, timeline, write_otlp_json

    session_id = args.session_id
    if args.otlp:
        path = write_otlp_json(session_id, args.otlp)
        if path is None:
            print(f"(could not write {args.otlp})")
            return 1
        print(f"OTLP/JSON written: {path}")
    print(timeline(session_id))
    # A trace with no recorded tool calls and no turns is likely a wrong
    # session id; the note says so instead of leaving a blank screen.
    if not spans(session_id):
        print(f"(no spans for session {session_id!r} — check the id with "
              "'delfin-agent sessions')")
    return 0


def cmd_where(args: argparse.Namespace) -> int:
    """Say where the dashboard is and how to walk back into it."""
    from . import where as _where

    flags = [flag for flag in ("--token", "--attach", "--ssh-config", "--write")
             if getattr(args, flag[2:].replace("-", "_"), False)]
    if flags:
        return _where.main(flags)
    print(_where.format_text(_where.dashboard(), _where.sessions()))
    return 0


def _boundary_options(args: argparse.Namespace) -> dict:
    """The session options ``capabilities`` and ``gate explain`` share,
    as keyword arguments of ``boundaries``."""
    from . import boundaries as _b
    return dict(
        cwd=getattr(args, "cwd", "") or os.getcwd(),
        permission_mode=getattr(args, "permission_mode", "") or "",
        add_dirs=tuple(getattr(args, "add_dirs", None) or ()),
        read_dirs=tuple(getattr(args, "read_dirs", None) or ()),
        allowed_tools=_allowed_tools(args),
        disallowed_tools=_disallowed_tools(args),
        asker=(getattr(args, "asker", "") or _b.default_asker(
            bool(getattr(args, "supervised", False)))),
        provider=getattr(args, "provider", "") or "",
    )


def cmd_capabilities(args: argparse.Namespace) -> int:
    """What a session started with these options may do; changes nothing.
    Exit 0, or 2 when the options name no gated session."""
    from . import boundaries as _b
    try:
        caps = _b.capabilities(**_boundary_options(args))
    except ValueError as exc:
        print(f"capabilities: {exc}", file=sys.stderr)
        return 2
    if getattr(args, "json", False):
        print(json.dumps(caps, indent=2, ensure_ascii=False))
    else:
        print(_b.format_capabilities(caps))
    return 0


def cmd_gate_explain(args: argparse.Namespace) -> int:
    """What the gate would do with one call, without making it.
    Exit 0 whatever the verdict (it is in the output), 2 on a usage
    error."""
    from . import boundaries as _b
    tool = getattr(args, "tool", "") or "bash"
    words = [w for w in (getattr(args, "command", None) or [])]
    if words and words[0] == "--":
        words = words[1:]
    if tool == "bash":
        if not words:
            print("gate explain: name the command after --", file=sys.stderr)
            return 2
        arguments = {"command": " ".join(words)}
    elif tool == "write_file":
        if not getattr(args, "path", ""):
            print("gate explain: --tool write_file needs --path",
                  file=sys.stderr)
            return 2
        arguments = {"path": args.path,
                     "content": getattr(args, "content", "") or ""}
    else:
        arguments = {}
    try:
        said = _b.explain(tool, arguments, **_boundary_options(args))
    except ValueError as exc:
        print(f"gate explain: {exc}", file=sys.stderr)
        return 2
    if getattr(args, "json", False):
        print(json.dumps(said, indent=2, ensure_ascii=False))
    else:
        print(_b.format_explain(said))
    return 0


def _add_boundary_flags(p: argparse.ArgumentParser) -> None:
    """The session options of ``capabilities`` and ``gate explain``: the
    same spellings ``chat`` takes, plus who would be asked."""
    p.add_argument("--json", action="store_true",
                   help="Print JSON instead of text")
    p.add_argument("--permission-mode", default="", dest="permission_mode",
                   choices=["", "plan", "default", "acceptEdits",
                            "bypassPermissions"],
                   help="The mode to report on (default: what chat would "
                        "start in)")
    p.add_argument("--cwd", default="",
                   help="Workspace (default: the current directory)")
    p.add_argument("--add-dir", action="append", default=[],
                   dest="add_dirs", metavar="PATH",
                   help="As chat --add-dir (repeatable)")
    p.add_argument("--read-dir", action="append", default=[],
                   dest="read_dirs", metavar="PATH",
                   help="As chat --read-dir (repeatable)")
    p.add_argument("--allowed-tools", default="", dest="allowed_tools",
                   metavar="a,b,c", help="As chat --allowed-tools")
    p.add_argument("--disallowed-tools", default="", dest="disallowed_tools",
                   metavar="a,b,c", help="As chat --disallowed-tools")
    p.add_argument("--asker", default="",
                   choices=["", "terminal", "approvals", "nobody"],
                   help="Who could answer a question (default: a terminal "
                        "on stdin, else approvals with --supervised, else "
                        "nobody)")
    p.add_argument("--supervised", action="store_true",
                   help="As chat --supervised")
    p.add_argument("--provider", default="",
                   help="kit, openai or ollama (default: the configured "
                        "provider)")


def cmd_jobs(args: argparse.Namespace) -> int:
    """The cluster's jobs from a terminal: read-only, answer when asked.

    `jobs watch --on|--off|--status` manages the monitor daemon's
    setting; the listing itself never starts or signals anything. A host
    without squeue gets one line and exit 0 — missing is not failing.
    """
    import json as _json

    from . import cli_jobs as _cj

    action = getattr(args, "jobs_action", None)
    if action == "watch":
        if getattr(args, "on", False):
            print(_cj.watch_set("on"))
            return 0
        if getattr(args, "off", False):
            print(_cj.watch_set("off"))
            return 0
        print(_cj.watch_report())
        return 0

    try:
        rows = _cj.collect_job_rows()
    except _cj.SchedulerUnavailable as exc:
        print(_cj.scheduler_note(0, exc))
        return 0
    if getattr(args, "json", False):
        print(_json.dumps(rows, indent=2))
    else:
        print(_cj.render_jobs(rows))
    return 0


def _presence_key_for(name: str, session_id: str) -> str:
    """The address this session is known by, to a person and to a peer.

    The name from ``-n`` when there is one, because that is what a person
    types; the head of the session id otherwise, so no session is
    unreachable for want of a name. One function, because the operator's
    inbox, the published approval and ``session_message`` must all mean
    the same session by the same word.
    """
    name = str(name or "").strip()
    if name:
        return name
    return str(session_id or "")[:8]


def _session_id_for(name: str) -> str:
    """The audit-log session id behind a -n name, or the name itself.

    The log records the engine's id; a supervisor knows the session by
    the name it was started with. Presence carries both, so the lookup
    is one read and costs nothing when the name IS an id.
    """
    name = str(name or "").strip()
    if not name:
        return ""
    try:
        from . import session_presence as _pres
        for record in _pres.open_sessions():
            if str(record.get("key") or "") == name:
                return str(record.get("session_id") or "") or name
    except Exception:
        pass
    return name


def cmd_watch(args: argparse.Namespace) -> int:
    """Follow every step, and above all every refusal.

    Supervising five sessions meant grepping the audit log by hand, once
    per question. What was wanted twice over: what a session is doing
    now, and what it TRIED to do and was refused -- because that is the
    line to stop on and the moment to step in.
    """
    import time as _time

    from . import audit_log as _al

    session = _session_id_for(getattr(args, "session", ""))
    denied = bool(getattr(args, "denied", False))
    # Start a little way back, so the first screen has context rather
    # than an empty prompt that looks like nothing is happening.
    back = max(0, int(getattr(args, "since", 40) or 0))
    rows, offset = _al.steps_since(0, session=session, denied_only=denied)
    for row in rows[-back:] if back else []:
        print(_al.render_step(row), flush=True)
    limit = float(getattr(args, "seconds", 0) or 0)
    deadline = _time.monotonic() + (limit if limit > 0 else 1e9)
    try:
        while _time.monotonic() < deadline:
            rows, offset = _al.steps_since(offset, session=session,
                                           denied_only=denied)
            for row in rows:
                print(_al.render_step(row), flush=True)
            _time.sleep(1.0)
    except KeyboardInterrupt:
        pass
    return 0


def cmd_provenance(args: argparse.Namespace) -> int:
    """Print the provenance stamp of one arm as JSON (paired_compare).
    Reads code, git and an environment; writes nothing."""
    from . import paired_compare as _pc
    pid = None
    if args.pid:
        try:
            pid = int(args.pid)
        except ValueError:
            print(f"--pid must be a number, got {args.pid!r}", file=sys.stderr)
            return 2
    stamp = _pc.provenance(code=list(args.code or ()), command=args.command,
                           pid=pid, prefixes=list(args.prefix or ()),
                           cwd=args.cwd)
    if args.save:
        try:
            path = _pc.save_stamp(args.save, stamp)
        except ValueError as exc:
            print(str(exc), file=sys.stderr)
            return 2
        print(f"saved stamp {args.save!r}: {path}", file=sys.stderr)
    print(json.dumps(stamp, indent=1, sort_keys=True))
    return 0


def cmd_compare_arms(args: argparse.Namespace) -> int:
    """Judge two arm result files as a paired comparison
    (paired_compare.compare). Exit 0 when a comparison stands and a
    person should look, 3 when it was refused or blocked, 2 on bad
    input. Reads only."""
    from . import paired_compare as _pc

    _stamp = _pc.load_stamp

    defaults = {}
    for item in args.default or ():
        name, sep, value = item.partition("=")
        if not sep or not name:
            print(f"--default takes VAR=VALUE, got {item!r}", file=sys.stderr)
            return 2
        defaults[name] = value
    try:
        a = _pc.load_records(args.arm_a)
        b = _pc.load_records(args.arm_b)
        stamp_a, stamp_b = _stamp(args.stamp_a), _stamp(args.stamp_b)
        rerun = _pc.load_records(args.rerun_a) if args.rerun_a else None
    except (OSError, ValueError) as exc:
        print(f"cannot read an input: {exc}", file=sys.stderr)
        return 2
    verdict = _pc.compare(
        a, b, stamp_a=stamp_a, stamp_b=stamp_b, axis=list(args.axis or ()),
        key=args.key, ok=args.ok, defects=list(args.defect or ()),
        quality=args.quality, lower_is_better=args.lower_is_better,
        min_n=args.min_n, probe=args.probe, sample=args.sample,
        defaults=defaults,
        allow_shared_switches=True if args.allow_shared_switches else None,
        rerun_a=rerun, prereg=args.prereg)
    if args.json:
        print(json.dumps(verdict, indent=1, default=str))
    else:
        print(_pc.format_report(verdict))
    return _pc.exit_code(verdict)


def cmd_report(args: argparse.Namespace) -> int:
    """What one agent session actually did: tools, files, commands,
    tests, denials, cost. `--json` prints the SessionReport itself,
    otherwise the terminal rendering.

    `--since <time>` switches to the round report: ONE report over every
    session with activity since that point (relative like `12h`/`3d`, or
    a date `2026-09-18`), optionally filtered by `--name PREFIX`. The
    hand tally this replaces covered tool calls/failures, dialogues,
    denials, >5min calls, commits, tokens, ttft and endpoint errors.
    """
    import dataclasses
    import json as _json

    from . import session_report as _sr

    since = (getattr(args, "since", "") or "").strip()
    sid = (getattr(args, "session", "") or "").strip()
    if since and not sid:
        from . import round_report as _rpt
        try:
            since_s = _rpt.parse_since(since)
        except ValueError as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        data = _rpt.collect(
            since_s=since_s,
            name=(getattr(args, "name", "") or "").strip(),
            # The real compaction archive, so the per-session compaction
            # count in the round report is this machine's own number.
            archive_root=None)
        print(_rpt.render_text(data))
        return 0

    if not sid:
        # No id: the most recently updated session, via the same source
        # `session ls` reads (sorted by updated_at, missing dir -> []).
        from . import session_store as _ss
        sid = (_ss.latest_session() or {}).get("session_id", "")
        if not sid:
            print("ERROR: no sessions found", file=sys.stderr)
            return 1

    report = _sr.collect_session_report(sid)
    if getattr(args, "json", False):
        print(_json.dumps(dataclasses.asdict(report), indent=2))
    else:
        print(_sr.render_terminal(report))
    return 0


def _subcommand_names(parser: argparse.ArgumentParser) -> frozenset[str]:
    """The subcommands the parser really registers.

    Derived from the parser instead of duplicated in a constant. A hand-
    kept list is one edit away from routing a real subcommand into the
    chat prompt, which would turn a typo'd command into a billed turn.
    """
    for action in parser._actions:
        if isinstance(action, argparse._SubParsersAction):
            return frozenset(action.choices)
    return frozenset()


def _route_argv(argv: list[str], known: frozenset[str]) -> list[str]:
    """Prepend the implicit ``chat`` subcommand.

    A bare invocation and anything starting with a flag is the session in
    the current directory; a first token that names a real subcommand, or
    asks for help, is left exactly as typed.
    """
    if not argv:
        return ["chat"]
    if argv[0] in known or argv[0] in ("-h", "--help", "--version"):
        return list(argv)
    return ["chat", *argv]


def _add_agent_flags(p: argparse.ArgumentParser, *,
                     mode_default: str = "solo",
                     max_tokens_default: int = 4096) -> None:
    """Engine-selection flags shared by ``chat`` and ``run``.

    One definition, so the two front doors cannot drift into offering
    different spellings of the same choice.
    """
    # The old help advertised `quick`, one of the pipeline modes retired
    # from the dashboard months ago — so the CLI kept offering a mode the
    # product had dropped, and the flag took it because there is no
    # `choices=`. Named modes only, and an unknown one falls back to solo
    # rather than failing, which is what the engine already does.
    p.add_argument("--mode", default=mode_default,
                   help="Agent mode: solo / dashboard / office / research "
                        "(plan is a permission profile, not a mode)")
    p.add_argument("--backend", default="", choices=["", "api", "cli"],
                   help="api (direct Anthropic) or cli (subprocess)")
    p.add_argument("--provider", default="",
                   help="claude / openai / kit / ollama (not given: the "
                        "configured one for chat, else the one the model "
                        "id names)")
    p.add_argument("--model", default="",
                   help="Model name (provider-specific)")
    p.add_argument("--effort", default="",
                   help="low/medium/high/xhigh")
    # 0 means "the role decides" (engine: max_tokens or role budget).
    # The old 4096 default silently cut a solo session to an eighth of
    # its 32768 role budget -- on a model that thinks invisibly, a
    # large edit then ended at finish_reason length with no tool call
    # and the turn stood at the prompt (operator's point d, 2026-09-22).
    # Only chat passes 0: run keeps its historical default, and the
    # scheduler and the benchmark run through run.
    p.add_argument("--max-tokens", type=int, default=max_tokens_default,
                   dest="max_tokens",
                   help="Max output tokens; 0 = the role's budget")
    p.add_argument("--cwd", default="", help="Run in this directory")


def cmd_watchpost(args: argparse.Namespace) -> int:
    """The read-only look-out over the user's own account."""
    from delfin.watchpost import scan as wp_scan
    from delfin.watchpost import scheduler_entry
    home = Path(args.home) if args.home else Path.home()
    report = Path(args.report) if args.report else \
        home / ".delfin" / "watchpost"
    if args.watchpost_action == "enable":
        info = scheduler_entry.enable(
            every_minutes=args.every_minutes or 30)
        print(f"watchpost scheduled every {info['every_seconds'] // 60}m "
              f"(entry {info['id']})")
        return 0
    if args.watchpost_action == "disable":
        removed = scheduler_entry.disable()
        print("watchpost schedule removed" if removed
              else "no watchpost schedule was set")
        return 0
    if args.watchpost_action == "baseline":
        findings = wp_scan.run_scan(home, report, mode="baseline")
        print(f"baseline reset; {len(findings)} current findings recorded")
    elif args.watchpost_action == "diff":
        findings = wp_scan.run_scan(home, report, mode="diff")
        print(f"{len(findings)} new finding(s) since baseline")
    else:
        findings = wp_scan.run_scan(home, report, mode="report")
        print(f"{len(findings)} finding(s)")
    for f in findings:
        loc = f"{f.path}:{f.line}" if f.line else f.path
        print(f"[{f.severity}] {f.check} {loc} -- {f.what}")
    if args.json:
        import json as _json
        print(_json.dumps([f.as_dict() for f in findings], indent=2))
    return 0


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="delfin-agent",
        description="The .delfin agent, in a terminal",
    )
    sub = p.add_subparsers(dest="cmd", required=True)

    # chat — the implicit subcommand a bare `delfin-agent` routes to
    chat = sub.add_parser("chat", help="Agent session in the current directory")
    chat.add_argument("prompt", nargs="*",
                      help="Opening prompt (optional)")
    chat.add_argument("-p", "--print", default="", dest="print_prompt",
                      metavar="PROMPT",
                      help="Answer one prompt and exit (non-interactive)")
    chat.add_argument("--output-format", default="text",
                      choices=["text", "json", "stream-json"],
                      dest="output_format",
                      help="text (default), json (one object at the end) or "
                           "stream-json (one object per line as the turn "
                           "runs; needs -p or a piped prompt)")
    chat.add_argument("--session", default="",
                      help="Session ID to resume, or 'latest'")
    chat.add_argument("-c", "--continue", action="store_true",
                      dest="continue_session",
                      help="Continue this directory's last conversation")
    chat.add_argument("-r", "--resume", nargs="?", const="", default=None,
                      metavar="ID",
                      help="Resume a session by id, or pick from a list")
    chat.add_argument("--new", action="store_true", dest="new_session",
                      help="Start a fresh session (the default)")
    chat.add_argument("--fork-session", action="store_true",
                      dest="fork_session",
                      help="Work on a copy, leaving the original untouched")
    chat.add_argument("-n", "--name", default="", dest="session_name",
                      help="Name this session, so --resume can find it")
    # The one way a run with no terminal can ask anybody anything. Opt-in
    # on purpose: without it a headless run keeps the behaviour it has,
    # where the Self-Modification Guard refuses outright because there is
    # no callback to ask. See delfin/agent/file_confirm.py.
    chat.add_argument("--supervised", action="store_true",
                      help="Send confirmations to the approvals directory "
                           "instead of refusing them (for a run with no "
                           "terminal; answer with 'delfin-agent approvals')")
    chat.add_argument("--permission-mode", default="", dest="permission_mode",
                      choices=["", "plan", "default", "acceptEdits",
                               "bypassPermissions"],
                      help="Start in this approval posture (default: the "
                           "settings default_mode, else default)")
    chat.add_argument("--confirm-timeout", type=float, default=0.0,
                      dest="confirm_timeout", metavar="SECONDS",
                      help="Let an unanswered question in the terminal expire "
                           "after SECONDS as 'not now' (never as a refusal); "
                           "0 waits for an answer (default)")
    chat.add_argument("--unattended", action="store_true",
                      help="Required alongside --permission-mode "
                           "bypassPermissions; nothing inside the working "
                           "folders will be asked, and reading outside them "
                           "is refused unless an approval can reach you")
    chat.add_argument("--add-dir", action="append", default=[],
                      dest="add_dirs", metavar="PATH",
                      help="Also writable this session (repeatable, never "
                           "persisted)")
    chat.add_argument("--read-dir", action="append", default=[],
                      dest="read_dirs", metavar="PATH",
                      help="Readable this session — the right answer to "
                           "'let it read my other repo'")
    chat.add_argument("--isolate", action="store_true",
                      help="Force filesystem isolation, and refuse to run "
                           "commands where this host cannot provide it")
    # The other half of a default that is ON. A protection nobody can
    # switch off when it blocks real work gets switched off badly -- by
    # editing settings, or by not using the tool at all.
    chat.add_argument("--no-isolate", action="store_true", dest="no_isolate",
                      help="Run shell commands unisolated — they may then "
                           "write outside the workspace")
    # Named for what it does, not for what the word suggests: what it
    # covers and the one thing it still cannot reach are both in the help
    # itself rather than left to be discovered — see _BARE_SKIPS and
    # _BARE_NOT_SKIPPED.
    # Two ways to narrow the tool surface of ONE session. Both land on the
    # permissions object the executor is handed, so a named tool is refused
    # when called and not merely hidden from the model — hiding alone would
    # leave every other route to the executor (an MCP backend, a sub-agent)
    # wide open. The subprocess CLI backend takes the allow list only, as
    # its own command line, and has no deny list; a run that asks for one
    # where there is none is told so at startup rather than left believing
    # the surface was narrowed.
    #
    # Neither list can WIDEN anything: the per-role gates are evaluated
    # separately, so a tool the role forbids stays forbidden however it is
    # named here. An empty allow list means no restriction. Deny wins.
    chat.add_argument("--allowed-tools", default="", dest="allowed_tools",
                      metavar="a,b,c",
                      help="Restrict this session to these tools (empty: no "
                           "restriction; cannot grant what the role forbids)")
    chat.add_argument("--disallowed-tools", default="", dest="disallowed_tools",
                      metavar="a,b,c",
                      help="Refuse these tools for this session — wins over "
                           "--allowed-tools when a name is on both")
    chat.add_argument("--bare", action="store_true",
                      help="Skip discovery for this run: "
                           f"{_and_list(_BARE_SKIPS)} "
                           f"({_BARE_NOT_SKIPPED} still run)")
    # Cumulative over the SESSION, not per turn — turns compose into an
    # unbounded run cost without this, which is what makes an unattended
    # run undeployable. Both land on the engine attributes _run_budget
    # reads, so nothing is written to the user's settings file.
    chat.add_argument("--max-turns", type=int, default=0, dest="max_turns",
                      metavar="N",
                      help="Stop a turn after N tool rounds (0: the "
                           "per-model default). Bounds one turn, not the "
                           "session — see --max-budget-usd for that")
    chat.add_argument("--max-budget-usd", type=float, default=0.0,
                      dest="max_budget_usd", metavar="USD",
                      help="Stop this session once measured spend reaches "
                           "this much (0: use the configured budget; an "
                           "unattended run without one: "
                           "agent.unattended_budget, $25)")
    chat.add_argument("--max-run-seconds", type=float, default=0.0,
                      dest="max_run_seconds", metavar="SECONDS",
                      help="Wall-clock ceiling for this session — the "
                           "dimension that still holds on a model with no "
                           "published rate (an unattended run without one: "
                           "agent.unattended_budget, 8 h)")
    chat.add_argument("--max-run-turns", type=int, default=0,
                      dest="max_run_turns", metavar="N",
                      help="Turn ceiling for this session (0: the "
                           "configured one; an unattended run without one: "
                           "agent.unattended_budget, 200)")
    # The consumer existed without the producer: ReplOptions.color has been
    # read as getattr(args, "color", "auto") since the first version, and no
    # parser ever declared the flag — so the one setting that decides
    # whether the session emits escape codes at all could not be set.
    chat.add_argument("--color", default="auto",
                      choices=["auto", "always", "never"],
                      help="Colour output (auto: only on a terminal, and "
                           "never when NO_COLOR is set)")
    _add_agent_flags(chat, max_tokens_default=0)
    chat.add_argument("-v", "--verbose", action="store_true")
    # Only this front door inherits the dashboard's saved provider/model.
    chat.set_defaults(func=cmd_chat, settings_defaults=True)

    # run
    run = sub.add_parser("run", help="Run one agent turn")
    run.add_argument("prompt", nargs="+", help="The user prompt")
    run.add_argument("--session", default="",
                     help="Session ID to resume, or 'latest'")
    _add_agent_flags(run)
    # `_build_engine` has always read this; only `chat` could set it, so a
    # headless turn took whatever the settings file said. A workspace or
    # profile carrying `default_mode: plan` therefore ran every headless
    # turn as a plan — and plan turns are exempt from claim grounding, so
    # the answer-side guards were off with nothing on screen to say so.
    # That is not a setting a benchmark or a probe should inherit by
    # accident; it is a decision the caller should be able to state.
    run.add_argument("--permission-mode", default="", dest="permission_mode",
                     choices=["", "plan", "default", "acceptEdits",
                              "bypassPermissions"],
                     help="Approval mode for this turn (default: settings)")
    run.add_argument("--unattended", action="store_true",
                     help="Required alongside --permission-mode "
                          "bypassPermissions, as on chat")
    run.add_argument("--json", action="store_true",
                     help="Emit JSON instead of plain text")
    run.add_argument("-v", "--verbose", action="store_true")
    run.set_defaults(func=cmd_run)

    # init
    init = sub.add_parser("init", help="Scaffold AGENTS.md + .delfin/")
    init.add_argument("path", nargs="?", default=".")
    init.add_argument("--force", action="store_true",
                      help="Overwrite existing files")
    init.set_defaults(func=cmd_init)

    # new-adapter (WB3): the hook is exactly these two lines; everything
    # else lives in delfin/agent/cli_new_adapter.py.
    from delfin.agent import cli_new_adapter
    _na = cli_new_adapter.build_parser(sub)
    _na.set_defaults(func=cli_new_adapter.handler)

    # new-workflow: the hook is exactly these two lines; everything else
    # lives in delfin/agent/cli_new_workflow.py.
    from delfin.agent import cli_new_workflow
    _nw = cli_new_workflow.build_parser(sub)
    _nw.set_defaults(func=cli_new_workflow.handler)

    # prcheck (package PC1)
    prcheck = sub.add_parser(
        "prcheck", help="The pre-PR self-check: run, verify, then hand over")
    prcheck.add_argument(
        "--base", default="origin/main",
        help="Branch to diff against (default: origin/main)")
    prcheck.add_argument(
        "--partner", action="append", default=[],
        help="Reviewer branch that must be merged in; repeatable")
    prcheck.set_defaults(func=_run_prcheck_subcommand)

    # bench — canned-task benchmark suite
    bench = sub.add_parser(
        "bench",
        help="Run / list / compare the canned-task benchmark suite",
    )
    bench_sub = bench.add_subparsers(dest="bench_action", required=False)

    bench_run = bench_sub.add_parser("run", help="Run the suite vs a model")
    bench_run.add_argument("--model", required=True,
                           help="Model name (e.g. kit.deepseek-v4-flash, opus)")
    bench_run.add_argument("--backend", default="", choices=["", "api", "cli"])
    bench_run.add_argument("--provider", default="",
                           help="claude / openai / kit")
    # Repeatable AND comma-separated. It used to be a plain single-value
    # option, so `--task a --task b` silently kept only the last one:
    # asking for six tasks ran one, and the run announced "on 1 tasks"
    # as if that had been the request.
    bench_run.add_argument("--task", default=[], action="append",
                           help=("Task IDs (default: all). Comma-separated, "
                                 "and the flag may be repeated — both forms "
                                 "accumulate. Shell globs select a family: "
                                 "--task 'generalist_*'"))
    bench_run.add_argument("--lang", default="de", choices=["de", "en"],
                           help=("Prompt language. en runs each selected "
                                 "task's prompt_en under the id <id>.en and "
                                 "skips tasks without one (default: de)"))
    bench_run.add_argument("--max-tokens", type=int, default=0,
                           dest="max_tokens")
    bench_run.add_argument(
        "--repeats", type=int, default=1,
        help=("N retries per task to defeat single-sample noise; "
              "result is median quality + majority success + quality_stdev "
              "(default: 1)"),
    )
    bench_run.add_argument(
        "--learning-home", default="", dest="learning_home",
        help=("Keep the agent's state under DIR across runs, one lineage "
              "per task and replicate: the two arms of a learning "
              "measurement (PREREG_W15.md) pass the same DIR. Never the "
              "user's home. Default: a fresh scratch home per attempt"),
    )

    bench_ls = bench_sub.add_parser("list", help="List packaged tasks")
    bench_ls.set_defaults(bench_action="list")

    bench_audit = bench_sub.add_parser(
        "audit",
        help=("Diagnose failed tasks in a run: prints model output excerpt + "
              "missing/violated patterns; flags likely pattern-bugs"),
    )
    bench_audit.add_argument(
        "run",
        help="Run JSONL file (absolute path or name in ~/.delfin/benchmark_runs/)",
    )
    bench_audit.set_defaults(bench_action="audit")

    bench_latest = bench_sub.add_parser(
        "latest", help="List recent run files in ~/.delfin/benchmark_runs/",
    )
    bench_latest.add_argument("--limit", type=int, default=10)
    bench_latest.set_defaults(bench_action="latest")

    bench_cmp = bench_sub.add_parser(
        "compare",
        help=("Diff two run files (baseline candidate), or with --model "
              "classify one run against the rolling per-task baseline — "
              "both are pure file comparisons, no API spend"),
    )
    bench_cmp.add_argument("baseline", nargs="?", default="",
                           help="Baseline JSONL run file (two-file mode)")
    bench_cmp.add_argument("candidate", nargs="?", default="",
                           help="Candidate JSONL run file (two-file mode)")
    bench_cmp.add_argument(
        "--model", default="",
        help=("Rolling-baseline mode: classify a run of this model as "
              "stable/improved/suspect/new vs the last-k history "
              "(noise-aware; no re-runs, no cost)"),
    )
    bench_cmp.add_argument(
        "--run", default="",
        help=("Run file to classify in --model mode (absolute path or name "
              "in the runs dir; default: the newest run of the model)"),
    )
    bench_cmp.add_argument(
        "--last-k", type=int, default=5, dest="last_k",
        help="Baseline window: per task, samples from the newest K run "
             "files containing it (default: 5)",
    )
    bench_cmp.add_argument("--runs-dir", default="", dest="runs_dir",
                           help="Override the runs directory "
                                "(default: ~/.delfin/benchmark_runs)")
    bench_cmp.add_argument("--json", action="store_true",
                           help="Emit raw JSON")
    bench_cmp.add_argument("--markdown", action="store_true",
                           help="Emit a markdown report (PR-body ready, "
                                "annotates profile commits between runs)")
    bench_cmp.add_argument(
        "--allow-instrument-change", action="store_true",
        dest="allow_instrument_change",
        help=("Compare two run files whose headers show a different "
              "scorer or task pack (refused otherwise)"))

    bench_gate = bench_sub.add_parser(
        "gate",
        help=("Judge a run against the committed reference standard "
              "(tests/fixtures/office_baseline.json) — a file comparison, "
              "no API spend"),
        description=(
            "Unlike `compare --model`, which measures against the last few "
            "runs on this machine, the reference is committed to the repo: "
            "it survives a clone, it is reviewed when it changes, and it "
            "does not drift downwards with the code it is meant to guard."
        ),
    )
    bench_gate.add_argument(
        "run", nargs="?", default="",
        help="Run JSONL (path or name in the runs dir; default: newest)")
    bench_gate.add_argument("--model", default="",
                            help="Pick the newest run of this model")
    bench_gate.add_argument("--baseline", default="", dest="baseline_path",
                            help="Reference file (default: the committed one)")
    bench_gate.add_argument("--runs-dir", default="", dest="runs_dir")
    bench_gate.add_argument("--json", action="store_true")
    bench_gate.set_defaults(bench_action="gate")

    bench_rebase = bench_sub.add_parser(
        "set-baseline",
        help=("Write a run as the new reference standard. Deliberate and "
              "explicit: nothing else ever writes one"),
        description=(
            "A reference that appears by itself records whatever state the "
            "code happened to be in, which is how a benchmark ends up "
            "certifying its own regression. This command exists so that "
            "moving the standard is an act somebody performs and signs."
        ),
    )
    bench_rebase.add_argument(
        "run", nargs="?", default="",
        help="Run JSONL to promote (path or name in the runs dir)")
    bench_rebase.add_argument("--model", default="",
                              help="Pick the newest run of this model")
    bench_rebase.add_argument("--baseline", default="", dest="baseline_path",
                              help="Where to write (default: the committed one)")
    bench_rebase.add_argument("--runs-dir", default="", dest="runs_dir")
    bench_rebase.add_argument("--note", default="",
                              help="Why this run is the new standard")
    bench_rebase.add_argument(
        "--yes", action="store_true",
        help="Required. Without it the command reports what it would do.")
    bench_rebase.set_defaults(bench_action="set-baseline")

    bench_nightly = bench_sub.add_parser(
        "nightly",
        help=("Unattended cycle: run suite → compare vs rolling history → "
              "re-check suspects (capped) → report + attention on "
              "CONFIRMED regressions. Scheduling stays opt-in: use "
              "`scheduler add-bench` explicitly (~$8/run for the 48-task "
              "KIT-Qwen suite at repeats=1)"),
        description=(
            "Full unattended benchmark cycle: runs the suite against "
            "--model, classifies every task against the per-task rolling "
            "baseline (majority success + median quality + noise-aware "
            "threshold), re-runs only suspect tasks with repeats under "
            "hard task/cost caps, writes a markdown report to "
            "<runs-dir>/reports/, and raises ONE attention event only "
            "for regressions the recheck CONFIRMED.\n\n"
            "Scheduling is strictly opt-in — this command never registers "
            "a scheduler entry. To run it nightly, run exactly:\n"
            "  delfin-agent scheduler add-bench "
            "--model <model> [--provider P] [--backend B] --every 24h\n"
            "  delfin-agent scheduler start\n"
            "Cost estimate: ~$8 per nightly run for the 48-task KIT-Qwen "
            "suite at repeats=1 (repeats multiply cost; recheck adds up "
            "to $3 more on suspect days)."),
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    bench_nightly.add_argument("--model", required=True,
                               help="Model name (e.g. kit.deepseek-v4-flash)")
    bench_nightly.add_argument("--provider", default="",
                               help="claude / openai / kit")
    bench_nightly.add_argument("--backend", default="",
                               choices=["", "api", "cli"])
    bench_nightly.add_argument(
        "--repeats", type=int, default=1,
        help="Samples per task for the main suite run (default: 1; "
             "cost scales linearly)")
    bench_nightly.add_argument(
        "--no-recheck", action="store_true", dest="no_recheck",
        help="Skip the suspect recheck stage (suspects stay unconfirmed; "
             "no attention event is raised)")
    bench_nightly.add_argument(
        "--last-k", type=int, default=5, dest="last_k",
        help="Baseline window: per task, samples from the newest K run "
             "files containing it (default: 5)")
    bench_nightly.add_argument("--runs-dir", default="", dest="runs_dir",
                               help="Override the runs directory "
                                    "(default: ~/.delfin/benchmark_runs)")
    bench_nightly.set_defaults(bench_action="nightly")

    bench.set_defaults(func=cmd_bench, bench_action="run")

    # credentials — secure key management
    cred = sub.add_parser(
        "credentials",
        help=("Manage API keys in the OS keyring or ~/.delfin/credentials.json "
              "(0600); no stored value is ever echoed back"),
    )
    cred_sub = cred.add_subparsers(dest="cred_action", required=False)

    cred_ls = cred_sub.add_parser(
        "list", help="List stored credentials (masked) and their store")
    cred_ls.set_defaults(cred_action="list")

    cred_mig = cred_sub.add_parser(
        "migrate",
        help=("Move file entries into the OS keyring; each leaves the file "
              "only after an equal read-back"))
    cred_mig.add_argument("names", nargs="*",
                          help="Credential names (default: every file entry)")

    cred_set = cred_sub.add_parser(
        "set", help="Store a credential (value read via getpass, no echo)",
    )
    cred_set.add_argument("name",
                          help="Credential name, e.g. KIT_TOOLBOX_API_KEY")

    cred_del = cred_sub.add_parser("delete", help="Remove a stored credential")
    cred_del.add_argument("name")

    cred.set_defaults(func=cmd_credentials, cred_action="list")

    # session
    sess = sub.add_parser("session", help="Session inspection")
    sess_sub = sess.add_subparsers(dest="session_action", required=True)
    ls = sess_sub.add_parser("ls", help="List recent sessions")
    ls.add_argument("--limit", type=int, default=20)
    srch = sess_sub.add_parser("search", help="Grep across session chats")
    srch.add_argument("query")
    exp = sess_sub.add_parser(
        "export", help="Export a session as a replayable notebook")
    exp.add_argument("session_id",
                     help="Session id, or 'latest' for the most recent")
    exp.add_argument("--notebook", default="",
                     help="Write the notebook here (default: <session_id>.ipynb)")
    sess.set_defaults(func=cmd_session)

    # messages -- delivery receipts for one sender, scoped to --from
    msgs_p = sub.add_parser("messages",
                            help="Session-message delivery receipts")
    msgs_sub = msgs_p.add_subparsers(dest="messages_action", required=True)
    msgs_ls = msgs_sub.add_parser("ls",
                                  help="List the mail the sender knows about")
    msgs_ls.add_argument("--from", dest="from_key", default="",
                         help="Sender key; without it nothing is shown")
    msgs_status = msgs_sub.add_parser(
        "status", help="Receipt for one message (queued/delivered/read/"
                        "unknown)")
    msgs_status.add_argument(
        "message_id", help="The message id to ask about")
    msgs_status.add_argument("--from", dest="from_key", default="",
                             help="Sender key; without it status is unknown")
    msgs_p.set_defaults(func=cmd_messages)

    # watchpost — read-only look-out over the user's own account
    wp = sub.add_parser(
        "watchpost",
        help="Read-only intrusion-trace scan of your own account "
             "(never changes anything)")
    wp_sub = wp.add_subparsers(dest="watchpost_action", required=True)
    wp_run = wp_sub.add_parser(
        "run", help="Scan once and report findings")
    wp_run.add_argument("--home", default="",
                        help="Home directory to look at (default: yours)")
    wp_run.add_argument("--report", default="",
                        help="Report directory (default: ~/.delfin/watchpost)")
    wp_run.add_argument("--json", action="store_true",
                        help="Also print findings as JSON")
    wp_sub.add_parser(
        "diff", help="Report only what is new since the baseline")
    wp_sub.add_parser(
        "baseline", help="Reset the baseline (first run, or on purpose)")
    wp_en = wp_sub.add_parser(
        "enable", help="Schedule a recurring watchpost run")
    wp_en.add_argument("--every", dest="every_minutes", type=int, default=30,
                       metavar="30m", help="Interval in minutes (default 30)")
    wp_sub.add_parser("disable", help="Remove the watchpost schedule")
    wp.set_defaults(func=cmd_watchpost, home="", report="", json=False)

    # memory — tidy the typed-memory store
    mem = sub.add_parser(
        "memory",
        help="Inspect and tidy the per-project memory store")
    mem_sub = mem.add_subparsers(dest="memory_action", required=False)
    mem_tidy = mem_sub.add_parser(
        "tidy",
        help="Report near-duplicates and unrecalled model notes (default: "
             "changes nothing)")
    mem_tidy.add_argument("--apply", action="store_true",
                          help="Carry out the proposal; retired notes are "
                               "moved to <store>/retired, never deleted")
    mem_tidy.add_argument("--workspace", default="",
                          help="Which project's store (default: this "
                               "directory)")
    mem.set_defaults(func=cmd_memory, memory_action="tidy")

    # knowledge -- the global knowledge folder: views for everyone,
    # curation for the knowledge admin (code dialog on /dev/tty only)
    kno = sub.add_parser(
        "knowledge",
        help="Global knowledge: lint, index, dupes, stats; admin curation "
             "(merge, move, retire, rename, accept, reject) as local branches")
    kno.add_argument("words", nargs=argparse.REMAINDER,
                     help="lint | index | dupes | stats | list | search <w> | "
                          "queue | accept <q> | reject <q> <reason> | merge "
                          "<keep> <drop> | move <entry> <topic> | retire "
                          "<entry> <reason> | rename <entry> <new> | topic "
                          "add <slug> <text>")
    kno.set_defaults(func=cmd_knowledge)

    # approvals — the other side of a supervised headless session
    appr = sub.add_parser(
        "approvals",
        help="Answer the confirmations a supervised headless session is "
             "waiting on (see 'chat --supervised')")
    appr_sub = appr.add_subparsers(dest="approvals_action", required=False)
    appr_sub.add_parser("ls", help="What is waiting (the default)")
    appr_show = appr_sub.add_parser(
        "show", help="The full preview of one request, diff included")
    appr_show.add_argument("request_id")
    appr_ok = appr_sub.add_parser("approve", help="Let one request through")
    appr_ok.add_argument("request_id")
    appr_no = appr_sub.add_parser("deny", help="Turn one request away")
    appr_no.add_argument("request_id")
    appr_no.add_argument(
        "--reason", default="",
        help="Tell the agent why, and what to do instead — it reads this "
             "in the same turn instead of guessing")
    appr_ans = appr_sub.add_parser(
        "answer", help="Answer a choice question (ask_user_question): "
                       "pick by number or by option text")
    appr_ans.add_argument("request_id")
    appr_ans.add_argument(
        "choice", help="The pick: an option's number (2) or its text, "
                       "comma-separated when the question allows more")
    for _p in (appr_ok, appr_no, appr_ans):
        _p.add_argument(
            "--code", default="",
            help="The one-time code the waiting session shows (its "
                 "terminal, its dashboard row); asked for on a terminal "
                 "when left out")
    appr_watch = appr_sub.add_parser(
        "watch", help="One line per new request, as they arrive")
    appr_watch.add_argument("--seconds", type=float, default=0.0,
                            help="Stop after N seconds (0: keep watching)")
    appr.set_defaults(func=cmd_approvals, approvals_action="ls")

    # skills — the human side of skill proposals (accept/reject)
    sk = sub.add_parser(
        "skills",
        help="Skill proposals: what sessions learned, waiting for a "
             "person to accept or reject")
    sk_sub = sk.add_subparsers(dest="skills_action", required=False)
    sk_sub.add_parser("proposals",
                      help="What is waiting (the default)")
    sk_show = sk_sub.add_parser(
        "show", help="Full preview of one proposal: text, evidence, "
                     "safety findings")
    sk_show.add_argument("name")
    sk_ok = sk_sub.add_parser("accept", help="Make a pending proposal "
                                             "a live skill")
    sk_ok.add_argument("name")
    sk_no = sk_sub.add_parser("reject", help="Turn a proposal away "
                                             "(moved to rejected/, kept)")
    sk_no.add_argument("name")
    sk_no.add_argument(
        "--reason", default="",
        help="Why it is rejected — the next session reads this instead "
             "of guessing the same thing again")
    sk.set_defaults(func=cmd_skills, skills_action="proposals")

    # bug — triage / watch the user bug-report archive (maintainer tool)
    bug = sub.add_parser(
        "bug",
        help="File a bug report (bug report) or triage the archive: "
             "auto-analyse + propose a fix (read-only; never moves to "
             "Solved automatically)",
    )
    bug_sub = bug.add_subparsers(dest="bug_action", required=False)

    bug_rep = bug_sub.add_parser(
        "report",
        help="File a bug report for a saved session (the dashboard's "
             "report, from a shell)")
    bug_rep.add_argument("--session", default="latest",
                         help="Session id, unique id prefix, or 'latest' "
                              "(newest session of --workspace)")
    bug_rep.add_argument("-m", "--message", default="",
                         help="One-line description (asked for when absent)")
    bug_rep.add_argument("--yes", action="store_true",
                         help="Write without the preview confirmation")
    bug_rep.add_argument("--workspace", default="",
                         help="Workspace whose latest session is meant "
                              "(default: current directory)")
    bug_rep.set_defaults(bug_action="report")

    bug_ls = bug_sub.add_parser("ls", help="List un-triaged reports")
    bug_ls.add_argument("--archive", default="",
                        help="Archive dir (default: settings / ~/.delfin/agent_bugs)")
    bug_ls.set_defaults(bug_action="ls")

    bug_tri = bug_sub.add_parser(
        "triage",
        help="Analyse un-triaged reports + propose a fix (read-only; writes "
             "triage.md + fix_proposal.patch; does NOT move to Solved)",
    )
    bug_tri.add_argument("name", nargs="?",
                         help="Only this report dir (default: all pending)")
    bug_tri.add_argument("--archive", default="")
    bug_tri.set_defaults(bug_action="triage")

    bug_solve = bug_sub.add_parser(
        "solve",
        help="Move a report to Solved/ (do this only AFTER the fix is "
             "applied + verified)",
    )
    bug_solve.add_argument("name", help="Report dir name to mark solved")
    bug_solve.add_argument("--archive", default="")
    bug_solve.set_defaults(bug_action="solve")

    bug_watch = bug_sub.add_parser(
        "watch", help="Watcher daemon: poll archive + triage new reports",
    )
    bug_watch.add_argument("--interval", type=int, default=0,
                           help="Poll seconds (default: settings interval_s)")
    bug_watch.add_argument("--once", action="store_true",
                           help="One pass then exit (no daemon loop)")
    bug_watch.add_argument("--archive", default="")
    bug_watch.set_defaults(bug_action="watch")

    bug.set_defaults(func=cmd_bug, bug_action="ls")

    # schedule -- the user's own schedules, run by the agent service
    user_sched = sub.add_parser(
        "schedule",
        help="Your schedules (run by the agent service): list | add <when> <task> | "
             "show | ok | pause | resume | delete | run | log | autostart")
    user_sched.add_argument("words", nargs=argparse.REMAINDER,
                            help="e.g. add täglich 07:30 --cwd ~/proj Fasse die Jobs zusammen")

    def _cmd_user_schedule(args):
        from .schedule_commands import cli_main
        return cli_main(list(args.words or []))
    user_sched.set_defaults(func=_cmd_user_schedule)

    # scheduler — headless executor for schedule_wakeup / cron entries
    sched = sub.add_parser(
        "scheduler",
        help="Headless daemon that executes scheduled entries "
             "(~/.delfin/cron.json) without the dashboard",
    )
    sched_sub = sched.add_subparsers(dest="scheduler_action", required=False)

    sched_start = sched_sub.add_parser(
        "start", help="Start the daemon (detached; survives shell close)")
    sched_start.set_defaults(scheduler_action="start")

    sched_status = sched_sub.add_parser(
        "status", help="Daemon state + entries sorted by next fire time")
    sched_status.set_defaults(scheduler_action="status")

    sched_stop = sched_sub.add_parser(
        "stop", help="Stop the daemon (SIGTERM; finishes the current turn)")
    sched_stop.set_defaults(scheduler_action="stop")

    sched_bench = sched_sub.add_parser(
        "add-bench",
        help=("OPT-IN: register a recurring `bench nightly` entry "
              "(~$8 per run for the 48-task KIT-Qwen suite at repeats=1). "
              "Never registered automatically — running this command IS "
              "the consent"),
    )
    sched_bench.add_argument("--model", required=True,
                             help="Model to benchmark nightly")
    sched_bench.add_argument("--provider", default="",
                             help="claude / openai / kit")
    sched_bench.add_argument("--backend", default="",
                             choices=["", "api", "cli"])
    sched_bench.add_argument("--repeats", type=int, default=1,
                             help="Samples per task per run (default: 1; "
                                  "cost scales linearly)")
    sched_bench.add_argument("--every", default="24h",
                             help="Interval token: 24h / 12h / 1d "
                                  "(default: 24h)")
    sched_bench.set_defaults(scheduler_action="add-bench")

    sched.set_defaults(func=cmd_scheduler, scheduler_action="status")

    stop_all_p = sub.add_parser(
        "stop-all",
        help="Emergency stop: end every agent of yours on every login node "
             "(dashboard kernels, shells, MCP servers, daemons); schedules "
             "are disabled and nothing starts on its own afterwards",
    )
    stop_all_p.add_argument("--yes", action="store_true",
                            help="Do not ask first")
    stop_all_p.add_argument("--check", action="store_true",
                            help="Only show what is open, kept, scheduled or "
                                 "running on its own; change nothing. Exit "
                                 "code 1 when something is")
    stop_all_p.add_argument("--reason", default="",
                            help="Recorded with the stop")
    stop_all_p.set_defaults(func=cmd_stop_all)

    pause_p = sub.add_parser(
        "pause",
        help="Pause one session by its -n name: it stops at the next tool "
             "call, starts no new turn and is not woken until resumed")
    pause_p.add_argument("session", help="The session's -n name")
    pause_p.add_argument("--reason", default="", help="Recorded with the pause")
    pause_p.set_defaults(func=cmd_pause)
    resume_p = sub.add_parser("resume", help="Resume a paused session")
    resume_p.add_argument("session", help="The session's -n name")
    resume_p.set_defaults(func=cmd_resume)

    # swarm -- build and run agent swarms (package SW1). The body lives
    # in cli_swarm.py, owned and tested there; this is the only hunk in
    # cli.py (same pattern as the prcheck subcommand).
    from .cli_swarm import main as _swarm_main
    _swarm_p = sub.add_parser(
        "swarm",
        help="Build and run agent swarms (init/plan/start/status)")
    _swarm_p.set_defaults(func=_swarm_main)

    pins_p = sub.add_parser(
        "mcp-pins",
        help="MCP tool descriptions pinned on first use, and changes that "
             "wait for you (accept SERVER [TOOL])")
    pins_sub = pins_p.add_subparsers(dest="pins_action")
    pins_sub.add_parser("list", help="Pinned servers and waiting changes")
    accept_p = pins_sub.add_parser(
        "accept", help="Accept a changed or new tool of a server")
    accept_p.add_argument("server", help="The server's name")
    accept_p.add_argument("tool", nargs="?", default="",
                          help="One tool (default: every waiting one)")
    pins_p.set_defaults(func=cmd_mcp_pins)

    # doctor — aggregate prerequisite health report
    doctor = sub.add_parser(
        "doctor",
        help="Check all prerequisites in one report (docs index, "
             "credentials, binaries, deps, MCP, scheduler, inbox, "
             "benchmark truth, memory store, disk); exit 1 on any FAIL",
    )
    doctor.add_argument("--workspace", default="",
                        help="Workspace directory (default: current dir)")
    doctor.add_argument("--repair", action="store_true",
                        help="Ask about, and apply with approval, every "
                             "step that fixes a WARN/FAIL; refuses when stdin "
                             "is not a terminal")
    doctor.set_defaults(func=cmd_doctor)

    # audit verify -- the audit log's hash chain (audit_log.verify)
    audit_p = sub.add_parser(
        "audit", help="Check the audit log for edited, removed or cut-off "
                      "records")
    audit_sub = audit_p.add_subparsers(dest="audit_action")
    audit_verify = audit_sub.add_parser(
        "verify", help="Verify the hash chain of the audit log and its "
                       "rotated weeks; exit 1 when one is not intact")
    audit_verify.add_argument("--log", default="",
                              help="One log file (default: ~/.delfin/"
                                   "audit.log and every rotated week)")
    audit_verify.add_argument("--json", action="store_true",
                              help="Emit the reports as JSON")
    audit_p.set_defaults(func=cmd_audit)

    # capabilities / gate explain -- the boundary, read before a session
    # starts. Both read only; see delfin/agent/boundaries.py.
    caps = sub.add_parser(
        "capabilities",
        help="What a session here may do: mode, writable and readable "
             "roots, deny paths, isolation on this host, who is asked, git "
             "role, tool lists (changes nothing)",
    )
    _add_boundary_flags(caps)
    caps.set_defaults(func=cmd_capabilities)

    gate = sub.add_parser(
        "gate", help="Ask the permission gate about a call without making it")
    gate_sub = gate.add_subparsers(dest="gate_action", required=True)
    gate_ex = gate_sub.add_parser(
        "explain",
        help="allow / ask / refuse for one command or write, with the "
             "reason code; nothing is run, written or recorded",
    )
    _add_boundary_flags(gate_ex)
    gate_ex.add_argument("--tool", default="bash",
                         choices=["bash", "write_file"],
                         help="The tool to ask about (default: bash)")
    gate_ex.add_argument("--path", default="",
                         help="With --tool write_file: the target")
    gate_ex.add_argument("--content", default="",
                         help="With --tool write_file: the text (default: "
                              "empty)")
    gate_ex.add_argument("command", nargs="*",
                         help="With bash: the command, after --")
    gate_ex.set_defaults(func=cmd_gate_explain)

    # report — what one session actually did, from the recorded sources
    report = sub.add_parser(
        "report",
        help="What one agent session did: tools, files, commands, tests, "
             "denials, cost (default: the most recent session)",
    )
    report.add_argument("--session", default="",
                        help="Session ID (default: the most recent session)")
    report.add_argument(
        "--since", default="",
        help="Round report: one report over every session with activity "
             "since this point (12h / 3d / 2026-09-18) instead of one")
    report.add_argument(
        "--name", default="",
        help="With --since: only sessions whose id starts with this prefix")
    report.add_argument("--json", action="store_true",
                        help="Print the SessionReport as JSON instead of "
                             "the terminal rendering")
    report.set_defaults(func=cmd_report)

    # provenance / compare-arms -- a paired comparison of two runs that
    # refuses skewed pairs and null tests (paired_compare.py). Read only.
    prov = sub.add_parser(
        "provenance",
        help="Print the stamp of one experiment arm: code hash, git head "
             "and diff, platform, behaviour environment (reads only)")
    prov.add_argument("--code", action="append", default=[],
                      help="A file or directory of the instrument (repeat)")
    prov.add_argument("--command", default="",
                      help="The arm's command line; its VAR=value, env and "
                           "export assignments are applied to this shell's "
                           "environment")
    prov.add_argument("--pid", default="",
                      help="A runner already running: read its environment, "
                           "start time and directory")
    prov.add_argument("--prefix", action="append", default=[],
                      help="A further behaviour-variable prefix (repeat); "
                           "the checkout's own NAME_ is always taken")
    prov.add_argument("--cwd", default="",
                      help="Directory whose git checkout is stamped "
                           "(default: here, or the runner's)")
    prov.add_argument("--save", default="",
                      help="Also keep the stamp, once, under this name in "
                           "~/.delfin/provenance (compare-arms takes the name)")
    prov.set_defaults(func=cmd_provenance)

    cmp_p = sub.add_parser(
        "compare-arms",
        help="Judge two arm result files as one paired comparison: same "
             "population, one instrument, declared axis, never-worse cases "
             "(reads only; exit 3 when refused or blocked)")
    cmp_p.add_argument("arm_a", help="Arm A results (JSONL or JSON list)")
    cmp_p.add_argument("arm_b", help="Arm B results (JSONL or JSON list)")
    cmp_p.add_argument("--stamp-a", default="",
                       help="Arm A stamp: a saved name or a JSON file")
    cmp_p.add_argument("--stamp-b", default="",
                       help="Arm B stamp: a saved name or a JSON file")
    cmp_p.add_argument("--axis", action="append", default=[],
                       help="What differs on purpose: an env variable, "
                            "'code' or 'platform' (repeat)")
    cmp_p.add_argument("--default", action="append", default=[],
                       help="VAR=VALUE: the code's default of an axis "
                            "variable, so an unset arm resolves")
    cmp_p.add_argument("--key", default="", help="Item id field")
    cmp_p.add_argument("--ok", default="", help="Per-item success field")
    cmp_p.add_argument("--defect", action="append", default=[],
                       help="Per-item defect count/flag field (repeat)")
    cmp_p.add_argument("--quality", default="",
                       help="Per-item quality value (higher is better)")
    cmp_p.add_argument("--lower-is-better", action="store_true",
                       help="With --quality: lower values are better")
    cmp_p.add_argument("--min-n", type=int, default=100,
                       help="Paired items below which nothing is concluded "
                            "automatically (default 100)")
    cmp_p.add_argument("--probe", type=int, default=24,
                       help="Leading items whose identity flags a null test "
                            "(default 24)")
    cmp_p.add_argument("--sample", type=int, default=100,
                       help="Size of the stratified sample to look at "
                            "(default 100)")
    cmp_p.add_argument("--rerun-a", default="",
                       help="A second run of arm A: differing items are "
                            "nondeterminism and block")
    cmp_p.add_argument("--prereg", default="",
                       help="The pre-registration file; reported whether it "
                            "predates the first stamp")
    cmp_p.add_argument("--allow-shared-switches", action="store_true",
                       help="Accept behaviour variables set identically in "
                            "both arms")
    cmp_p.add_argument("--json", action="store_true",
                       help="Print the verdict as JSON")
    cmp_p.set_defaults(func=cmd_compare_arms)

    # watch — every step every session takes, as it takes it, and above
    # all the ones the gate refused.
    watch_p = sub.add_parser(
        "watch",
        help="Follow what the open sessions do, step by step; --denied "
             "keeps only what the gate refused",
    )
    watch_p.add_argument("--session", default="",
                         help="One session: its -n name or the head of its id")
    watch_p.add_argument("--denied", action="store_true",
                         help="Only the refusals — what a session tried and "
                              "was not allowed to do")
    watch_p.add_argument("--seconds", type=float, default=0.0,
                         help="Stop after this long (default: until Ctrl+C)")
    watch_p.add_argument("--since", type=int, default=40,
                         help="How many recent steps to print before "
                              "following (default 40)")
    watch_p.set_defaults(func=cmd_watch)

    # where — the dashboard runs on one node, in one multiplexer, behind
    # one forwarded port, and the next login lands anywhere.
    where_p = sub.add_parser(
        "where",
        help="Where the dashboard runs and how to get back into it: node, "
             "tmux session, port, and a link per kept session",
    )
    where_p.add_argument("--attach", action="store_true",
                         help="Walk into the tmux session the note names "
                              "(ssh + port forward when it runs elsewhere)")
    where_p.add_argument("--ssh-config", action="store_true",
                         help="Print a ready ssh config entry for it")
    where_p.add_argument("--token", action="store_true",
                         help="Print the dashboard's login token and nothing "
                              "else (for the login page; links carry none)")
    where_p.add_argument("--write", action="store_true",
                         help="With --ssh-config: append the entry instead "
                              "of printing (never overwrites)")
    where_p.set_defaults(func=cmd_where)

    # sessions — what ran here before, without starting anything
    sessions_p = sub.add_parser(
        "sessions",
        help="List earlier agent sessions: id, age, model and the task "
             "they were given (resume one with --resume <id>)",
    )
    sessions_p.add_argument("--limit", type=int, default=20,
                            help="How many to list (default 20)")
    sessions_p.add_argument("--load", action="store_true",
                            help="One line per open session on this host: "
                                 "processes, RSS, CPU cores, alarms "
                                 "(reporting only — nothing is signalled)")
    sessions_p.add_argument("--watch", type=int, default=None, metavar="SECONDS",
                            help="With --load: repeat every SECONDS "
                                 "(Ctrl-C to stop)")
    sessions_p.set_defaults(func=cmd_sessions)

    # trace — one run read as a trace (OT1): text timeline from the
    # existing tool_trace / turn_metrics / session_store records, plus
    # an OTLP/JSON file with --otlp. Read-only; nothing is sent anywhere.
    trace_p = sub.add_parser(
        "trace",
        help="One run as a trace: turns, model calls, tool calls "
             "(text timeline; --otlp writes an OTLP/JSON file)",
    )
    trace_p.add_argument("session_id",
                         help="Session id to render (see 'delfin-agent "
                              "sessions')")
    trace_p.add_argument("--otlp", default="", metavar="PATH",
                         help="Also write the span tree as an OTLP/JSON "
                              "file a collector can import (no network)")
    trace_p.set_defaults(func=cmd_trace)

    # jobs — the cluster's jobs, read-only, without starting a session.
    # Same query and watch files the monitor daemon uses (cli_jobs is
    # the surface, job_monitor the mechanism).
    jobs_p = sub.add_parser(
        "jobs",
        help="Your watched cluster jobs: id, state, elapsed, workspace "
             "(read-only); `jobs watch` manages the monitor daemon",
    )
    jobs_p.add_argument("--json", action="store_true",
                        help="Print the rows as JSON instead of the table")
    jobs_p.set_defaults(func=cmd_jobs)

    jobs_sub = jobs_p.add_subparsers(dest="jobs_action")
    jobs_watch = jobs_sub.add_parser(
        "watch",
        help="Start, stop or report the job-monitor daemon (its setting "
             "agent.job_monitor.enabled; the watch loop is LLM-free)",
    )
    jobs_watch.add_argument("--on", action="store_true",
                            help="Enable job monitoring (persists the "
                                 "setting; starts nothing)")
    jobs_watch.add_argument("--off", action="store_true",
                            help="Disable job monitoring (a running daemon "
                                 "exits on its next pass)")
    jobs_watch.add_argument("--status", action="store_true",
                            help="Report daemon, setting and watched count "
                                 "(default)")

    # ops status -- checkouts, dashboard servers, sessions and leftovers in
    # one read-only view (delfin/agent/ops_status.py). Every flag reads:
    # the bash auto-allow entry for this command depends on that, and
    # tests/test_ops_status_reads_the_machinery.py pins the flag set.
    def _cmd_ops(args):
        from delfin.agent import ops_status
        return ops_status.main(json_out=bool(args.json),
                               repo=str(args.repo or ""))

    ops_p = sub.add_parser(
        "ops",
        help="DELFIN's own machinery: `ops status` lists checkouts and "
             "worktrees, dashboard servers, sessions and leftovers with a "
             "safe next step per problem (read-only)",
    )
    ops_sub = ops_p.add_subparsers(dest="ops_action", required=True)
    ops_status_p = ops_sub.add_parser(
        "status",
        help="Read-only overview: checkouts (branch / detached / held by), "
             "dashboard servers of yours, session locks, stale locks and "
             "worktrees; one 'problem → next step' line each; changes nothing",
    )
    ops_status_p.add_argument("--json", action="store_true",
                              help="Print the status as JSON "
                                   "(schema delfin.ops_status/1)")
    ops_status_p.add_argument("--repo", default="",
                              help="A checkout of the repository to inspect "
                                   "(default: the current directory, else "
                                   "the checkout DELFIN runs from)")
    ops_status_p.set_defaults(func=_cmd_ops)

    # ops health -- the person sees at a glance which session is stuck and
    # why (UX2, wave 17). READ-ONLY: classify + collect read, nothing
    # writes; --wait is the loop the person or operator starts that emits
    # the attention events (emit_attention only, interval >= 60 s).
    def _cmd_ops_health(args):
        import time as _t

        from delfin.agent import ops_interface as _oi
        from delfin.agent import session_health as _sh
        from delfin.agent import stop_all as _sa
        try:
            from delfin import user_settings
            agent_cfg = (user_settings.load_settings().get("agent") or {})
        except Exception:
            agent_cfg = {}

        interval = max(60, int(getattr(args, "interval", 0) or 0))
        state = {"seen": {}}
        stop_seen = {}
        stop = _sa.last_stop()
        while True:
            try:
                records = _sh.collect()
            except KeyboardInterrupt:
                return 0
            if getattr(args, "wait", False):
                events = _sh.watch_events(records, _t.time(), state,
                                          settings=agent_cfg,
                                          stop_all=stop,
                                          stop_all_seen=stop_seen)
                for name, text in events:
                    print(f"attention: {text}", flush=True)
            rows = []
            for rec in records:
                verdict = _sh.classify(rec, _t.time(), agent_cfg)
                rows.append({
                    "name": rec.get("name"),
                    "session_id": rec.get("session_id"),
                    "verdict": verdict["verdict"],
                    "reason": verdict["reason"],
                    "next": verdict["next"],
                })
            if getattr(args, "json", False):
                print(_oi.dumps(rows), flush=True)
            else:
                for r in rows:
                    print(f"{r['name']}  {r['verdict']}: {r['reason']} "
                          f"[{r['next']}]", flush=True)
            if not getattr(args, "wait", False):
                return 0
            _sh.pause_watch(interval)
            stop = _sa.last_stop()

    ops_health_p = ops_sub.add_parser(
        "health",
        help="One line per open session: verdict (working / waiting on "
             "you / stalled / paused / gone / ...), why, and the next "
             "step; --wait loops and emits the attention events",
    )
    ops_health_p.add_argument("--json", action="store_true",
                              help="Print the rows as JSON")
    ops_health_p.add_argument("--wait", action="store_true",
                              help="Loop forever (interval >= 60 s) and "
                                   "emit ONE attention event per session "
                                   "entering a stuck state, deduped")
    ops_health_p.add_argument("--interval", type=int, default=60,
                              help="Seconds between --wait ticks (min 60)")
    ops_health_p.set_defaults(func=_cmd_ops_health)

    # The supervisor's side of running sessions, as data
    # (delfin/agent/ops_interface.py): events and state to read; a message,
    # an escalation to the person and a pause as the only actions. There is
    # no approve: a person answers approvals (#202, #207).
    def _cmd_ops_interface(args):
        from delfin.agent import ops_interface as _oi
        act = args.ops_action

        def _out(obj):
            print(_oi.dumps(obj), flush=True)

        try:
            if act == "sessions":
                for row in _oi.sessions():
                    _out(row)
                return 0
            if act == "session":
                row = _oi.session(args.name)
                if row is None:
                    print(f"ERROR: no open session {args.name!r}",
                          file=sys.stderr)
                    return 2
                _out(row)
                return 0
            if act == "approvals":
                for row in _oi.pending_approvals():
                    _out(row)
                return 0
            if act == "events":
                if not args.follow:
                    # Printed first, saved second: a reader that dies
                    # while printing reads the batch again.
                    state = _oi.load_cursor(args.cursor)
                    for ev in _oi.next_events(state, idle_s=args.idle_s):
                        _out(ev)
                    _oi.save_cursor(state, args.cursor)
                    return 0
                for ev in _oi.follow(args.cursor, idle_s=args.idle_s,
                                     interval=args.interval,
                                     seconds=args.seconds):
                    _out(ev)
                return 0
            if act == "wait":
                ev = _oi.wait(args.until, cursor=args.cursor,
                              timeout=args.timeout, idle_s=args.idle_s,
                              interval=args.interval)
                if ev is None:
                    return 3
                _out(ev)
                return 0
            if act == "send":
                text = args.text
                if args.file:
                    # Bounded: the inbox keeps 4000 characters, and a
                    # read of /dev/zero or a FIFO must not take the
                    # supervisor's memory or hang it.
                    with open(args.file, "rb") as fh:
                        text = fh.read(_oi.FILE_MAX).decode("utf-8",
                                                            "replace")
                _out(_oi.send(args.name, text or ""))
                return 0
            if act == "escalate":
                _out({"attention_id": _oi.escalate(args.request_id,
                                                   args.note or "")})
                return 0
            if act == "pause":
                _out(_oi.pause(args.name, args.reason))
                return 0
        except (_oi.OpsRefused, ValueError, OSError) as exc:
            print(f"ERROR: {exc}", file=sys.stderr)
            return 2
        except KeyboardInterrupt:
            return 0
        return 2

    for _name, _help in (
            ("sessions", "Open sessions on this host as JSON lines (read-only)"),
            ("approvals", "Approvals any session waits on, as JSON lines "
                          "(read-only; answering is the person's)"),
    ):
        ops_sub.add_parser(_name, help=_help).set_defaults(
            func=_cmd_ops_interface)
    _p = ops_sub.add_parser("session", help="One open session as JSON "
                                             "(read-only)")
    _p.add_argument("name", help="The -n name or the head of the id")
    _p.set_defaults(func=_cmd_ops_interface)
    for _name, _help in (
            ("events", "Events since the cursor as JSON lines: sessions "
                       "opened/exited/idle, turns ended, approvals "
                       "open/closed, messages, refused steps, commits"),
            ("wait", "Block until one event matches --until; exit 3 on "
                     "timeout"),
    ):
        _p = ops_sub.add_parser(_name, help=_help)
        _p.add_argument("--cursor", default="default",
                        help="Name of the reader's cursor (one per reader)")
        _p.add_argument("--idle-s", type=int, default=900,
                        help="A session with no activity for this long is "
                             "reported idle once")
        _p.add_argument("--interval", type=float, default=5.0)
        if _name == "events":
            _p.add_argument("--follow", action="store_true",
                            help="Keep printing new events")
            _p.add_argument("--seconds", type=float, default=0.0,
                            help="With --follow: stop after this long")
        else:
            _p.add_argument("--until", required=True,
                            help="kind=approval_open,session=s11 or "
                                 "idle_s>=1800; every term must hold")
            _p.add_argument("--timeout", type=float, default=0.0)
        _p.set_defaults(func=_cmd_ops_interface)
    _p = ops_sub.add_parser("send", help="Leave a message in a session's "
                                          "inbox, as the operator (never as "
                                          "the user)")
    _p.add_argument("name")
    _g = _p.add_mutually_exclusive_group(required=True)
    _g.add_argument("--text")
    _g.add_argument("--file")
    _p.set_defaults(func=_cmd_ops_interface)
    _p = ops_sub.add_parser("escalate", help="Put a waiting approval in "
                                              "front of the person; it stays "
                                              "open")
    _p.add_argument("request_id")
    _p.add_argument("--note", default="")
    _p.set_defaults(func=_cmd_ops_interface)
    _p = ops_sub.add_parser("pause", help="Pause a session at its next safe "
                                           "point")
    _p.add_argument("name")
    _p.add_argument("--reason", required=True)
    _p.set_defaults(func=_cmd_ops_interface)

    # T6: run a session's tests on a SLURM compute node from a node-local
    # copy. The handler lives in delfin/agent/slurm_tests.py; it is imported
    # on use so an import error there cannot break every other command.
    def _cmd_test_on_slurm(args):
        from delfin.agent.slurm_tests import cmd_test_on_slurm
        return cmd_test_on_slurm(args)

    t6 = sub.add_parser(
        "test-on-slurm",
        help="Render and submit a node-local SLURM job that runs the given "
             "test paths; prints the job id",
    )
    t6.add_argument("repo", help="Path to the repository")
    t6.add_argument("ref", help="Git ref to run at (copied node-local)")
    t6.add_argument("partition", help="SLURM partition, e.g. cpu")
    t6.add_argument("minutes", type=int, help="Job time limit in minutes")
    t6.add_argument("tests", nargs="+",
                    help="Test paths under tests/, e.g. tests/test_x.py")
    t6.add_argument("--run-dir", default=".",
                    help="Dir to write the job and collect logs")
    t6.set_defaults(func=_cmd_test_on_slurm)

    # The terminal UI and the service its sessions live in. The handlers
    # live in delfin/agent/service/cli_service.py, imported on use like
    # test-on-slurm's.
    def _cmd_tui(args):
        from delfin.agent.service.cli_service import cmd_tui
        return cmd_tui(args)

    def _cmd_service(args):
        from delfin.agent.service.cli_service import cmd_service
        return cmd_service(args)

    tui = sub.add_parser(
        "tui",
        help="The terminal UI: one window over all your sessions, which keep "
             "running when it closes (sessions live in `delfin-agent service`)",
    )
    tui.add_argument("--cwd", default="", help="Workspace of the first session")
    tui.add_argument("-n", "--name", dest="session_name", default="",
                     help="Name of the first session")
    tui.add_argument("--model", default="")
    tui.add_argument("--provider", default="")
    tui.add_argument("--effort", default="")
    tui.add_argument("--permission-mode", default="",
                     choices=["", "plan", "default", "acceptEdits"])
    tui.add_argument("--session", "--resume", dest="session", default="",
                     help="Resume this saved session")
    tui.add_argument("--continue", dest="continue_session", action="store_true",
                     help="Resume the latest session in the workspace")
    tui.add_argument("--add-dir", dest="add_dirs", action="append", default=[])
    tui.add_argument("--read-dir", dest="read_dirs", action="append", default=[])
    tui.add_argument("--isolate", action="store_true")
    tui.add_argument("--no-isolate", action="store_true")
    tui.add_argument("-v", "--verbose", action="store_true",
                     help="Show the model's thinking, as chat --verbose does")
    tui.add_argument("--inline", action="store_true",
                     help="Draw under the shell's last line, the conversation in the "
                          "terminal's own scrollback, instead of taking the whole terminal")
    tui.add_argument("--setup", action="store_true",
                     help="Set the window up (a pinned Bun and its packages, checked) and "
                          "exit; done by itself on the first start unless "
                          "DELFIN_TUI_AUTO_INSTALL=0 or agent.tui_auto_install=false")
    # Scripted sessions, no model: for working on the window and its tests.
    tui.add_argument("--fake", action="store_true", help=argparse.SUPPRESS)
    tui.set_defaults(func=_cmd_tui)

    svc = sub.add_parser(
        "service",
        help="The agent service behind `delfin-agent tui`: status (its sessions), "
             "end a session, stop, restart, logs",
    )
    svc.add_argument("action", nargs="?", default="status",
                     choices=["status", "end", "stop", "restart", "logs"])
    svc.add_argument("target", nargs="?", default="",
                     help="end: the session's id, or a unique prefix of it")
    svc.add_argument("--all-idle", dest="all_idle", action="store_true",
                     help="end: every idle session (none that works or waits for you)")
    svc.add_argument("--after-turn", dest="after_turn", action="store_true",
                     help="end: let a running turn finish first (without it, the turn "
                          "is stopped at its next safe point, as ctrl+x twice does)")
    svc.add_argument("--yes", "-y", action="store_true",
                     help="end: without asking (for scripts)")
    svc.add_argument("--json", action="store_true",
                     help="status: as JSON, for scripts")
    svc.add_argument("--lines", type=int, default=40,
                     help="Lines of the log to show (logs)")
    svc.set_defaults(func=_cmd_service)

    return p


def _show_security_notices() -> None:
    """Print, once, what a security update changed in the settings."""
    try:
        from delfin.user_settings import load_settings, take_security_notices
        load_settings()                     # applies a pending update
        for notice in take_security_notices():
            print(notice, file=sys.stderr)
    except Exception:
        pass


def _stop_own_background_shells() -> None:
    """End the background shells this process started. Never raises."""
    try:
        from . import bash_jobs as _bj
        _bj.get_registry().stop_running()
    except Exception:
        pass


def main(argv: list[str] | None = None) -> int:
    # A probe or an interview driven through here is nobody's session:
    # with DELFIN_SCRATCH_STATE set it keeps its history, memory, profile
    # learning and logs under that directory instead of the user's, and
    # reads none of theirs. Applied before anything imports a sink.
    from . import state_paths as _state_paths
    _state_paths.scratch_state_from_environment()
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(_route_argv(argv, _subcommand_names(parser)))
    if getattr(args, "func", None) in (cmd_chat, cmd_run):
        guard_principles_or_exit()
        # An agent in a terminal ends with an emergency stop given anywhere
        # (stop_all), and with the dashboard it was started from, if any.
        from . import lifeline as _lifeline
        # The background shells this terminal agent started end with it,
        # as a dashboard's end with its kernel: at a normal exit, at an
        # interrupt, and when the lifeline or an emergency stop ends it
        # (os._exit skips atexit, so that path stops them first).
        import atexit as _atexit
        _atexit.register(_stop_own_background_shells)
        # This process holds the model's key: no command it runs may read it.
        from . import process_guard as _process_guard
        _process_guard.protect("terminal agent")
        _show_security_notices()
        _exported = _process_guard.exported_provider_keys()
        # Unconditional, because the environment is not where the worst
        # case lives: an alias sets a key for one command and leaves the
        # line in a shell file for good, so nothing ever exports it and
        # the old `if _exported:` never looked.
        #
        # Do it rather than ask for it: an exported key goes into the
        # 0600 store and the export is commented out, with a copy of the
        # file kept. Only what could not be decided is left to say.
        _done = _process_guard.put_exported_keys_away(_exported)
        if _done:
            print(_done, file=sys.stderr)
        elif _exported:
            print("Warning: "
                  + _process_guard.exported_key_advice(_exported),
                  file=sys.stderr)

        def _end_with_everything() -> None:
            _stop_own_background_shells()
            _lifeline.exit_now()

        _lifeline.watch(_end_with_everything)
    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
