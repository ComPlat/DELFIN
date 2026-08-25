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
import subprocess
import sys
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
    provider = args.provider or fallback.get("provider", "")
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
    return engine


def _allowed_tools(args: argparse.Namespace) -> list[str]:
    """``--allowed-tools a,b,c`` as a list, or empty for "no restriction"."""
    raw = getattr(args, "allowed_tools", "") or ""
    if isinstance(raw, (list, tuple)):
        parts = [str(p) for p in raw]
    else:
        parts = str(raw).split(",")
    return [p.strip() for p in parts if p.strip()]


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


def _run_once(engine, prompt: str, *, max_tokens: int = 4096) -> dict[str, Any]:
    """Stream a single turn and collect text + tool-calls + token-usage.

    AgentEngine's ``stream_response`` is callback-driven, not event-
    iterable: text arrives through ``on_token``, tool calls through
    ``on_tool_use``.  The function also returns the assembled text as
    a string.  Token usage is read from ``engine.token_usage`` after
    the call (cumulative for the engine; each benchmark task gets a
    fresh engine so the cumulative IS per-turn).
    """
    chunks: list[str] = []
    tool_calls: list[dict] = []
    error = ""

    def _on_token(text: str) -> None:
        if text:
            chunks.append(text)

    def _on_tool_use(name: str, input_json: str) -> None:
        try:
            inp = json.loads(input_json) if input_json else {}
        except (json.JSONDecodeError, TypeError):
            inp = {"raw": str(input_json)}
        tool_calls.append({"name": name, "input": inp})

    in_before = int((getattr(engine, "token_usage", {}) or {}).get("input", 0))
    out_before = int((getattr(engine, "token_usage", {}) or {}).get("output", 0))

    try:
        full_text = engine.stream_response(
            user_message=prompt,
            on_token=_on_token,
            on_tool_use=_on_tool_use,
            max_tokens=max_tokens,
        ) or ""
    except Exception as exc:
        error = str(exc)
        full_text = ""

    in_after = int((getattr(engine, "token_usage", {}) or {}).get("input", 0))
    out_after = int((getattr(engine, "token_usage", {}) or {}).get("output", 0))

    return {
        "text": (full_text or "".join(chunks)).strip(),
        "tool_calls": tool_calls,
        "input_tokens": max(0, in_after - in_before),
        "output_tokens": max(0, out_after - out_before),
        "error": error,
    }


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
            out.append({"role": role, "content": text})
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
        _ss.save_session(chat_messages=_display_messages(engine),
                         workspace=str(repo_root), **estate)
    except Exception as exc:
        print(f"WARN: session save failed: {exc}", file=sys.stderr)
    return sid


def cmd_run(args: argparse.Namespace) -> int:
    if not args.prompt:
        print("ERROR: prompt is required.", file=sys.stderr)
        return 2
    prompt = " ".join(args.prompt) if isinstance(args.prompt, list) else str(args.prompt)
    repo = Path(args.cwd or os.getcwd()).expanduser().resolve()
    if args.cwd:
        os.chdir(repo)

    try:
        engine = _build_engine(args)
    except Exception as exc:
        print(f"ERROR: engine init failed: {exc}", file=sys.stderr)
        return 3
    # The headless half has no banner, so the can't-deliver lines go to
    # stderr — stdout carries the answer and nothing else.
    for _note in _bounding_notices(args, engine):
        print(_note, file=sys.stderr)
    _resume_or_create(engine, args)
    import time as _time
    _t0 = _time.monotonic()
    out = _run_once(engine, prompt, max_tokens=args.max_tokens or 4096)
    sid = _save_session(engine, repo)

    # Learning signal: record the outcome so provider profiles learn from
    # CLI/headless usage too — previously only dashboard cycles fed the
    # profile, so KIT/CLI sessions contributed nothing.
    try:
        engine.record_cycle_outcome(
            "FAIL" if out["error"] else "PASS",
            prompt,
            error_type=("cli_error" if out["error"] else None),
            start_time=_t0,
        )
    except Exception:
        pass
    # Episodic memory: persist one compact per-session record so a future
    # session can recall similar past work (best-effort, LLM-free) —
    # without this the saved session state is write-only.
    try:
        from .episodes import build_episode_from_state, save_episode
        fields = build_episode_from_state(engine.export_state(), [])
        save_episode(
            sid,
            repo_root=repo,
            verdict="FAIL" if out["error"] else "PASS",
            **fields,
        )
    except Exception:
        pass
    # Close the eval loop opportunistically (opt-in, LLM-free, max 1/day).
    try:
        from .eval_loop import maybe_run_scheduled
        maybe_run_scheduled()
    except Exception:
        pass

    if args.json:
        payload = {**out, "session_id": sid}
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
    return 0 if not out["error"] else 1


def _persisted_default_mode(workspace: Path) -> tuple[str, str]:
    """``(mode, the file that declared it)``, or ``("", "")`` for none.

    Read from the files rather than from ``kit_settings.load()``, and the
    difference is the whole point: with no settings anywhere the loader
    still returns ``default_mode="default"``, its own fallback, which is
    indistinguishable from a user having written that word down. Asking
    the merged view therefore reports a decision nobody made — and it
    silently defeated the plan-first default, because resolve_posture
    correctly honours anything a user configured.

    The path travels with the mode because the banner names it. It used
    to be a literal in resolve_posture's signature, so a checked-out
    repository's settings file could raise the posture and the user was
    pointed at a file in their home directory that need not even exist.
    """
    import json

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
            forked = _ss.fork_session(sid)
            sid = str((forked or {}).get("session_id", "") or sid)
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

    return _claim_session(getattr(engine, "session_id", "") or sid)


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

    lines = [
        f"delfin-agent · {provider}/{model} · {role_mode}",
        f"workspace  {_tilde(workspace)}  ({where})",
    ]
    if perms_mode:
        lines.append(f"approval   {perms_mode}"
                     + (f"  [{why}]" if why else ""))
        if isolation_note:
            lines.append(isolation_note)
        elif perms_mode in ("default", "acceptEdits"):
            # Nobody has been told this. _bash_isolation_argv engages bwrap
            # only under bypassPermissions or a locked scope; in the
            # attended modes with the shipped "auto" setting it returns a
            # plain /bin/bash -c. A user reading "workspace confinement"
            # reasonably assumes a sandbox, and what is actually there is
            # path checking plus a regex list.
            lines.append(
                "isolation  off — a command the agent runs can still write "
                "outside the workspace (--isolate)")
    else:
        # Not decoration: without a permissions object the write and shell
        # tools refuse, and a user staring at a silent agent deserves the
        # reason on screen rather than in a tool result.
        lines.append(
            f"approval   none — the {provider} backend carries no permission "
            "gate, so file and shell tools will refuse")
    for extra in (_grant_line("writable", getattr(report, "granted_dirs", ())),
                  _grant_line("readable", getattr(report, "read_dirs", ())),
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
    lines.append("esc interrupt · shift+tab approval mode · /help · ctrl+d exit")
    return "\n".join(lines)


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


# What ``--bare`` cannot reach from this file, named once so the help text
# and the startup line cannot drift from each other. All three are
# discovered inside the turn, from the permissions workspace, and none has
# a per-session switch: hooks through ``hooks.load_hooks`` at every hook
# point, skills through ``skills.discover_skills`` while the tool surface
# is assembled, project memory through ``PromptLoader`` while the system
# prompt is built. A ``--bare`` that implied it covered them would be the
# silent non-delivery the flag exists to avoid.
_BARE_NOT_SKIPPED = "hooks, skills and project memory"


def _apply_bare(workspace: Path) -> bool:
    """Skip MCP server discovery for this run. True when it took.

    An MCP server definition is executable configuration: it is spawned
    with the parent environment while the tool surface is being ASSEMBLED,
    before any model output, and then answers every call routed to it.
    That is the one piece of discovery reachable from here, because the
    registry is cached per workspace and loads its configuration once —
    emptying it after that load means no server is started and no MCP tool
    is advertised, for this process only. Nothing on disk is touched, so
    the next run without the flag is unchanged.
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


def _bounding_notices(args: argparse.Namespace, engine) -> list[str]:
    """One line per bounding flag this run cannot actually honour.

    Same contract as the ``--isolate`` note: a flag that names a limit and
    then does not impose it is worse than no flag, because the user stops
    watching. Every line here is produced from what the constructed engine
    and client really are, never from a second copy of the rules that
    decided it.
    """
    notes: list[str] = []

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

    # An allow-list that the constructed client did not take. Asked of the
    # client itself rather than by repeating create_client's branch here:
    # the parameter is forwarded only to the CLI backend, which stores it
    # and turns it into --allowedTools, and is dropped without a word for
    # every other one. Comparing against what actually arrived is the only
    # form of this check that cannot go stale.
    wanted = _allowed_tools(args)
    if wanted:
        client = getattr(engine, "client", None)
        got = getattr(client, "allowed_tools", None)
        if [str(t) for t in (got or ())] != wanted:
            provider = str(getattr(engine, "provider", "") or "?")
            backend = str(getattr(engine, "backend", "") or "?")
            notes.append(
                f"tools      --allowed-tools REQUESTED but not applied — the "
                f"{provider}/{backend} backend carries no tool allow-list, so "
                f"all {len(wanted)} named tools and every other one stay "
                f"available")

    if getattr(args, "bare", False):
        took = getattr(args, "bare_mcp_skipped", False)
        notes.append(
            (f"bare       MCP servers skipped — {_BARE_NOT_SKIPPED} are "
             "discovered inside the turn and cannot be skipped from here"
             ) if took else
            (f"bare       REQUESTED but nothing was skipped — the MCP "
             f"registry could not be reached, and {_BARE_NOT_SKIPPED} are "
             "discovered inside the turn"))

    return notes


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
        args.prompt = [prompt or positional]
        args.session = getattr(args, "session", "") or ""
        args.json = (getattr(args, "output_format", "text") == "json")
        return cmd_run(args)

    if piped:
        print("Nothing on stdin, and no -p. Nothing to do.", file=sys.stderr)
        return 2

    # An interactive session is a conversation, not a document. The flag
    # is accepted on this subcommand because the one-shot half lives here
    # too; taking it silently and then emitting text is the shape of a
    # promise that is not kept.
    if getattr(args, "output_format", "text") == "json":
        print("--output-format json describes one answer, so it needs -p "
              "or a piped prompt; this session will print text.",
              file=sys.stderr)

    # Restored on the way out. The process usually ends right after, so
    # this looks unnecessary — but a function that moves the process and
    # never moves it back makes every later relative path in the caller
    # read a different file, and that is only harmless until it is called
    # from somewhere else.
    _cwd_before = os.getcwd()
    os.chdir(workspace)
    # Only the grants launch_guard accepted. It has already refused a
    # credential directory, a parent of the workspace and a forbidden
    # root by name, rather than letting the permissions object drop them
    # silently the way it does.
    args.extra_dirs = [str(p) for p in report.granted_dirs]
    args.read_only_dirs = [str(p) for p in report.read_dirs]
    isolation_note = ""
    if getattr(args, "isolate", False):
        import shutil as _shutil
        from .api_client import set_bash_isolation_override
        set_bash_isolation_override("bwrap")
        # Asking for isolation the host cannot give must not be silent.
        # _bash_isolation_argv falls back to plain /bin/bash when bwrap is
        # missing — correct, since there is nothing else it can do — but a
        # flag that promises containment and delivers none without saying
        # so is the same defect this wave exists to close, one layer up.
        if not _shutil.which("bwrap"):
            isolation_note = (
                "isolation  REQUESTED but unavailable — bubblewrap is not "
                "installed here, so commands run unisolated")
    try:
        engine = _build_engine(args)
    except Exception as exc:
        print(f"ERROR: engine init failed: {exc}", file=sys.stderr)
        os.chdir(_cwd_before)
        return 3
    if not _open_session(engine, args, workspace):
        os.chdir(_cwd_before)
        return 2

    persisted, persisted_from = _persisted_default_mode(workspace)
    posture, why = launch_guard.resolve_posture(
        flag_mode=getattr(args, "permission_mode", "") or "",
        persisted_mode=persisted,
        unattended_opt_in=bool(getattr(args, "unattended", False)),
        settings_path=persisted_from or "~/.delfin/settings.json",
    )
    try:
        engine.set_kit_permission_mode(posture)
    except Exception:
        pass          # a backend with no gate has no posture to set

    from .repl import ReplOptions, TerminalAgent
    from .terminal_confirm import TerminalConfirmBroker

    # Bound ONLY here, and only for a terminal. cmd_run, the benchmark and
    # the scheduler must keep the behaviour they have: without a callback
    # a write inside the workspace is allowed silently, which is right for
    # an unattended run and exactly what this command exists to change.
    broker = None
    bind = getattr(engine, "set_kit_confirm_callback", None)
    if sys.stdin.isatty() and callable(bind):
        broker = TerminalConfirmBroker(
            persist=lambda pat: engine.persist_kit_pattern(pat, kind="allow"),
            set_mode=engine.set_kit_permission_mode,
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


    notices = report.render()
    agent = TerminalAgent(engine, broker=broker, opts=ReplOptions(
        cwd=workspace,
        max_tokens=getattr(args, "max_tokens", 0) or 0,
        show_thinking=bool(getattr(args, "verbose", False)),
        color=getattr(args, "color", "auto"),
        banner=(_startup_banner(engine, report, workspace, why,
                                isolation_note,
                                tuple(_bounding_notices(args, engine)))
                + (f"\n\n{notices}" if notices else "")),
    ))
    try:
        return agent.run(positional)
    finally:
        _save_session(engine, workspace,
                      title=getattr(args, "session_name", "") or "")
        try:
            os.chdir(_cwd_before)
        except Exception:
            pass


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
                    1 for line in f.open(encoding="utf-8") if line.strip()
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
    wanted = {t.strip() for chunk in _asked
              for t in str(chunk).split(",") if t.strip()}
    if wanted:
        tasks = [t for t in tasks if t.id in wanted]
        if not tasks:
            print(f"ERROR: no task matched {sorted(wanted)}", file=sys.stderr)
            return 1
        # A name that matched nothing is a typo or a renamed task, and
        # running the rest without saying so hands back a number for a
        # different question than the one that was asked.
        missing = sorted(wanted - {t.id for t in tasks})
        if missing:
            print(f"WARNING: {len(missing)} requested task(s) do not exist "
                  f"and were skipped: {', '.join(missing)}", file=sys.stderr)

    model = (getattr(args, "model", "") or "").strip()
    if not model:
        print("ERROR: --model is required for `bench run`",
              file=sys.stderr)
        return 2
    max_tokens = int(getattr(args, "max_tokens", 1024) or 1024)
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
                     f"rate={result.success_rate:.0%}")
        print(f"  [{mark}] {task.id:<28} q={result.quality_0_100:>3}  "
              f"{result.duration_s:>5.1f}s  ${result.cost_usd:.4f}  "
              f"tool={result.tool_calls}{extra}",
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
    )
    path = _bm.write_run(results, model=model)
    s = _bm.summarise_run(results)
    print(f"\n{s['n_pass']}/{s['n_tasks']} passed "
          f"({s['pass_rate']:.0%})   "
          f"avg-quality {s['avg_quality']:.1f}   "
          f"total ${s['total_cost_usd']:.4f}   "
          f"{s['total_duration_s']:.1f}s")
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
            commit = subprocess.run(
                ["git", "rev-parse", "HEAD"], capture_output=True, text=True,
                cwd=str(Path(__file__).resolve().parents[2]),
            ).stdout.strip()
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


def cmd_bug(args: argparse.Namespace) -> int:
    """Triage / watch the user bug-report archive (auto-analyse → Solved)."""
    from . import bug_watcher as _bw
    from . import bug_report as _br
    action = getattr(args, "bug_action", "") or "ls"
    archive_arg = getattr(args, "archive", "") or ""
    archive = (Path(archive_arg).expanduser() if archive_arg
               else _bw.resolve_watch_archive())

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
            desc = str(rep.get("description") or "").replace(chr(10), " ")[:56]
            print(f"  {d.name:<44} {str(rep.get('model', '?')):<22} {desc}")
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
        for d in targets:
            print(f"→ analysing {d.name} …", flush=True)
            t = _bw.triage_report(d)
            if t.error:
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
            results = _bw.run_once(archive=archive, settings=settings)
            for t in results:
                tag = "ERROR" if t.error else "OK"
                print(f"  [{tag}] {t.report_name}: "
                      f"{(t.error or t.summary)[:150].strip()}")
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
        try:
            print(f"bug_watcher started (interval {cfg['interval_s']}s, "
                  f"archive={archive}).")
            return _bw.run_loop(settings=settings)
        finally:
            _bw.release_pid_lock()

    return 0


def cmd_credentials(args: argparse.Namespace) -> int:
    """Manage credentials in ~/.delfin/credentials.json (chmod 0600).

    No subcommand here echoes a stored value — `list` masks every key
    so the user can verify WHICH credentials are stored without exposing
    them.  Input is read via getpass so the value is never visible on
    screen and never lands in shell history.
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
        print(f"Credentials (file: {_cred.credentials_path()})")
        print()
        for name in sorted(items):
            info = items[name]
            src = info.get("source", "?")
            tag = "[env]" if src == "env" else "[file]"
            print(f"  {name:<28}  {info.get('value', ''):<14}  {tag}")
        return 0

    name = (getattr(args, "name", "") or "").strip()
    if not name:
        print("ERROR: credential name required.", file=sys.stderr)
        return 2

    if action == "delete":
        if _cred.delete_credential(name):
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
        ok = _cred.set_credential(name, value)
        if ok:
            print(f"Stored {name} = {_cred.mask(value)} "
                  f"in {_cred.credentials_path()} (chmod 0600)")
        else:
            print(f"{name} already stored with this value (no change).")
        return 0

    print(f"Unknown credentials action: {action}", file=sys.stderr)
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
    print(f"ERROR: unknown session action {args.session_action!r}",
          file=sys.stderr)
    return 2


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
        with log.open("a") as lf:
            _sp.Popen(
                [sys.executable, "-m", "delfin.agent.scheduler_daemon"],
                stdout=lf, stderr=lf,
                start_new_session=True,  # survives shell/dashboard close
            )
        print(f"Scheduler daemon started (detached); log: {log}")
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


def cmd_doctor(args: argparse.Namespace) -> int:
    """One-surface prerequisite report: docs index, credentials,
    binaries, Python deps, MCP servers, scheduler, attention inbox,
    benchmark ground truth, memory store, disk space.

    Exit code 1 when any check FAILs (warnings alone exit 0), so the
    command is scriptable as a pre-flight gate.
    """
    from . import doctor as _doc

    workspace = getattr(args, "workspace", "") or None
    results = _doc.run_doctor(workspace)
    print(_doc.format_doctor(results))
    return 1 if any(r.get("status") == "FAIL" for r in results) else 0


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
                   help="claude / openai / kit / ollama")
    p.add_argument("--model", default="",
                   help="Model name (provider-specific)")
    p.add_argument("--effort", default="",
                   help="low/medium/high/xhigh")
    p.add_argument("--max-tokens", type=int, default=max_tokens_default,
                   dest="max_tokens")
    p.add_argument("--cwd", default="", help="Run in this directory")


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
                      choices=["text", "json"], dest="output_format",
                      help="text (default) or json")
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
    chat.add_argument("--permission-mode", default="", dest="permission_mode",
                      choices=["", "plan", "default", "acceptEdits",
                               "bypassPermissions"],
                      help="Start in this approval posture (default: plan)")
    chat.add_argument("--unattended", action="store_true",
                      help="Required alongside --permission-mode "
                           "bypassPermissions; nothing will be asked")
    chat.add_argument("--add-dir", action="append", default=[],
                      dest="add_dirs", metavar="PATH",
                      help="Also writable this session (repeatable, never "
                           "persisted)")
    chat.add_argument("--read-dir", action="append", default=[],
                      dest="read_dirs", metavar="PATH",
                      help="Readable this session — the right answer to "
                           "'let it read my other repo'")
    chat.add_argument("--isolate", action="store_true",
                      help="Run shell commands under filesystem isolation")
    # Named for what it does, not for what the word suggests. The three it
    # cannot reach are named in the help itself rather than left to be
    # discovered — see _BARE_NOT_SKIPPED.
    # Only the claude CLI backend has a tool allow-list (create_client
    # forwards this to CLIClient, which spells it --allowedTools); the API
    # and OpenAI-compatible backends drop it. A run that asks for one on a
    # backend that has none is told so at startup rather than left
    # believing the surface was narrowed. There is deliberately no
    # --disallowed-tools: no deny-list machinery exists on any backend,
    # and a flag with nothing behind it is the defect, not the fix.
    chat.add_argument("--allowed-tools", default="", dest="allowed_tools",
                      metavar="a,b,c",
                      help="Restrict the agent to these tools "
                           "(claude CLI backend only; other backends have "
                           "no allow-list and will say so)")
    chat.add_argument("--bare", action="store_true",
                      help="Skip MCP server discovery for this run "
                           f"({_BARE_NOT_SKIPPED} are discovered inside the "
                           "turn and are NOT skipped)")
    # Cumulative over the SESSION, not per turn — turns compose into an
    # unbounded run cost without this, which is what makes an unattended
    # run undeployable. Both land on the engine attributes _run_budget
    # reads, so nothing is written to the user's settings file.
    chat.add_argument("--max-budget-usd", type=float, default=0.0,
                      dest="max_budget_usd", metavar="USD",
                      help="Stop this session once measured spend reaches "
                           "this much (0: use the configured budget)")
    chat.add_argument("--max-run-seconds", type=float, default=0.0,
                      dest="max_run_seconds", metavar="SECONDS",
                      help="Wall-clock ceiling for this session — the "
                           "dimension that still holds on a model with no "
                           "published rate")
    # The consumer existed without the producer: ReplOptions.color has been
    # read as getattr(args, "color", "auto") since the first version, and no
    # parser ever declared the flag — so the one setting that decides
    # whether the session emits escape codes at all could not be set.
    chat.add_argument("--color", default="auto",
                      choices=["auto", "always", "never"],
                      help="Colour output (auto: only on a terminal, and "
                           "never when NO_COLOR is set)")
    _add_agent_flags(chat)
    chat.add_argument("-v", "--verbose", action="store_true")
    # Only this front door inherits the dashboard's saved provider/model.
    chat.set_defaults(func=cmd_chat, settings_defaults=True)

    # run
    run = sub.add_parser("run", help="Run one agent turn")
    run.add_argument("prompt", nargs="+", help="The user prompt")
    run.add_argument("--session", default="",
                     help="Session ID to resume, or 'latest'")
    _add_agent_flags(run)
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

    # bench — canned-task benchmark suite
    bench = sub.add_parser(
        "bench",
        help="Run / list / compare the canned-task benchmark suite",
    )
    bench_sub = bench.add_subparsers(dest="bench_action", required=False)

    bench_run = bench_sub.add_parser("run", help="Run the suite vs a model")
    bench_run.add_argument("--model", required=True,
                           help="Model name (e.g. kit.qwen3.5-397b-A17b, opus)")
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
                                 "accumulate"))
    bench_run.add_argument("--max-tokens", type=int, default=1024,
                           dest="max_tokens")
    bench_run.add_argument(
        "--repeats", type=int, default=1,
        help=("N retries per task to defeat single-sample noise; "
              "result is median quality + majority success + quality_stdev "
              "(default: 1)"),
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
                               help="Model name (e.g. kit.qwen3.5-397b-A17b)")
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
        help=("Manage API keys in ~/.delfin/credentials.json (chmod 0600); "
              "no stored value is ever echoed back"),
    )
    cred_sub = cred.add_subparsers(dest="cred_action", required=False)

    cred_ls = cred_sub.add_parser("list", help="List stored credentials (masked)")
    cred_ls.set_defaults(cred_action="list")

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
    sess.set_defaults(func=cmd_session)

    # bug — triage / watch the user bug-report archive (maintainer tool)
    bug = sub.add_parser(
        "bug",
        help="Triage the bug-report archive: auto-analyse + propose a fix "
             "(read-only; never moves to Solved automatically)",
    )
    bug_sub = bug.add_subparsers(dest="bug_action", required=False)

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

    # doctor — aggregate prerequisite health report
    doctor = sub.add_parser(
        "doctor",
        help="Check all prerequisites in one report (docs index, "
             "credentials, binaries, deps, MCP, scheduler, inbox, "
             "benchmark truth, memory store, disk); exit 1 on any FAIL",
    )
    doctor.add_argument("--workspace", default="",
                        help="Workspace directory (default: current dir)")
    doctor.set_defaults(func=cmd_doctor)

    return p


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    argv = list(sys.argv[1:] if argv is None else argv)
    args = parser.parse_args(_route_argv(argv, _subcommand_names(parser)))
    try:
        return int(args.func(args) or 0)
    except KeyboardInterrupt:
        print("\nInterrupted.", file=sys.stderr)
        return 130


if __name__ == "__main__":
    sys.exit(main())
