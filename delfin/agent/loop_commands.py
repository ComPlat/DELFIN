"""``/loop`` and ``/auto``: the words, for every surface.

Input: the argument text after the command, the session's
:class:`~delfin.agent.auto_mode.AutoController` (``None`` on a surface that
has none), the engine, the workspace, and optionally ``ask`` (read one
answer line) and ``say`` (print a line now). Output: ``(text, prompt)`` --
the text to show and, for ``/loop new <name>``, a prompt to send to the
model instead of the raw line.

The terminal reaches this through its command table (``repl_commands``),
the window through the same table (its worker runs the terminal loop), the
dashboard by calling :func:`loop_command` / :func:`auto_command` from its
own dispatcher. One handler, so a ``/loop run`` means the same thing on the
three surfaces.

``/loop run`` of an unconfirmed loop shows the file and asks; with nobody to
ask it refuses. Confirmation (:func:`loops.confirm`) happens only here, in a
command the user typed: the agent can write a loop file but cannot make it
run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Callable, Optional

from . import loops

__all__ = ["loop_command", "auto_command", "LOOP_USAGE", "AUTO_USAGE"]

LOOP_USAGE = ("/loop [list] · /loop show|confirm|run|log <name> · /loop stop · "
              "/loop resume · /loop new <name> [--blank] · "
              "/loop <5m|2h|1d> <prompt>")
AUTO_USAGE = "/auto [on|off] · /auto needs [clear]"

_NO_CONTROLLER = ("auto mode and loops are not available on this surface "
                  "(no session controller).")


def _status(ctl) -> str:
    if ctl is None:
        return _NO_CONTROLLER
    line = ctl.status_line()
    if line:
        st = ctl.state()
        extra = []
        if st.get("paused"):
            extra.append(f"paused: {st['paused']} (/auto on resumes)")
        if st.get("branch"):
            extra.append(f"branch {st['branch']}")
        extra.append(f"left: {st.get('tokens_left', 0)} tokens, "
                     f"{st.get('wall_left_s', 0) // 60} min")
        return "auto: " + line + "\n  " + " · ".join(extra)
    tail = f"\nlast stop: {ctl.stop_report}" if ctl.stop_report else ""
    return "auto mode is off. " + AUTO_USAGE + tail


def auto_command(args: str, ctl, engine) -> str:
    words = (args or "").split()
    sub = words[0].lower() if words else ""
    if ctl is None:
        return _NO_CONTROLLER
    if sub in ("", "status"):
        return _status(ctl)
    if sub == "on":
        return ctl.enable(engine)
    if sub == "off":
        return ctl.disable("turned off by the user")
    if sub == "needs":
        if len(words) > 1 and words[1].lower() == "clear":
            n = len(ctl.needs_you)
            ctl.needs_you.clear()
            return f"cleared {n} item(s)."
        if not ctl.needs_you:
            return "nothing needs you."
        return "needs you:\n" + "\n".join(
            f"  #{n['id']} [{n['kind']}] {n['text']}"
            + (f" — {n['why']}" if n.get("why") else "")
            for n in ctl.needs_you)
    return "usage: " + AUTO_USAGE


def _listing(ctl, workspace) -> str:
    specs = loops.discover(workspace)
    lines = []
    running = getattr(ctl, "loop", None) if ctl is not None else None
    if running is not None:
        lines.append("running: " + ctl.status_line())
    if not specs:
        lines.append("no loops defined. /loop new <name> drafts one from this "
                     "conversation; files live in .delfin/loops/ (project) or "
                     "~/.delfin/loops/ (personal).")
        return "\n".join(lines)
    lines.append("loops:")
    for s in specs:
        if s.errors:
            mark = "invalid: " + s.errors[0]
        elif loops.is_confirmed(s):
            mark = "confirmed"
        else:
            mark = "not confirmed"
        lines.append(f"  {s.name} [{s.scope}, {mark}] — {s.goal[:70]}")
    lines.append("usage: " + LOOP_USAGE)
    return "\n".join(lines)


def _show(spec) -> str:
    try:
        text = Path(spec.path).read_text(encoding="utf-8") if spec.path else ""
    except OSError as exc:
        text = f"(cannot read: {exc})"
    lines = [f"── {spec.path}", text.rstrip(), "──"]
    if spec.errors:
        lines.append("invalid: " + "; ".join(spec.errors))
    else:
        lines.append(f"{spec.summary()}; budget {spec.wall_seconds // 60} min, "
                     f"{spec.max_tokens} tokens; stops after {spec.plateau} "
                     f"iterations without progress or {spec.repeat_failures} "
                     "identical failures.")
        lines.append("measure (an iteration may not change): "
                     + (", ".join(spec.measure) or "—"))
        lines.append("touch (an iteration may change only): "
                     + (", ".join(spec.touch) or "the workspace"))
        lines.append(f"branch: {spec.branch}; "
                     + ("confirmed" if loops.is_confirmed(spec) else
                        f"not confirmed — /loop confirm {spec.name}"))
    return "\n".join(lines)


def _draft_prompt(name: str, rest: str) -> str:
    return (
        f"Write a loop definition for what this conversation has been "
        f"doing, as `.delfin/loops/{name}.md`, following the `loops` skill "
        f"(load it first). {('The user adds: ' + rest + ' ') if rest else ''}"
        "Use the conversation's own commands for `check` when one measures "
        "the goal; list in `measure` every test or script the result is "
        "judged by; keep the budget small. Do not run the loop and do not "
        "confirm it: when the file is written, show it and tell the user to "
        f"read it with `/loop show {name}` and start it with "
        f"`/loop run {name}` (which asks them to confirm).")


def loop_command(args: str, ctl, engine, workspace, *,
                 ask: Optional[Callable[[str], str]] = None,
                 say: Optional[Callable[[str], None]] = None
                 ) -> tuple[str, str]:
    """``(text, prompt)``; see the module docstring. Never raises."""
    try:
        return _loop(args, ctl, engine, Path(workspace), ask, say)
    except Exception as exc:
        return f"/loop failed: {exc}", ""


def _loop(args, ctl, engine, workspace, ask, say) -> tuple[str, str]:
    words = (args or "").split()
    sub = words[0].lower() if words else ""
    name = words[1] if len(words) > 1 else ""
    if sub in ("", "list", "ls"):
        return _listing(ctl, workspace), ""
    if sub == "new":
        if not name:
            return "usage: /loop new <name> [--blank] [what it should do]", ""
        rest = " ".join(w for w in words[2:] if w != "--blank")
        if "--blank" in words:
            try:
                path = loops.write_template(workspace, name, goal=rest)
            except FileExistsError:
                return f"{name} exists already: /loop show {name}", ""
            except ValueError as exc:
                return str(exc), ""
            return (f"wrote {path} (not confirmed). Edit it, then "
                    f"/loop show {name} and /loop run {name}."), ""
        if loops.load(name, workspace) is not None:
            return f"{name} exists already: /loop show {name}", ""
        return "", _draft_prompt(name, rest)
    if sub in ("show", "confirm", "log") and not name:
        return f"usage: /loop {sub} <name>", ""
    if sub == "show":
        spec = loops.load(name, workspace)
        return (_show(spec) if spec else f"no loop called {name}."), ""
    if sub == "confirm":
        spec = loops.load(name, workspace)
        if spec is None:
            return f"no loop called {name}.", ""
        if spec.errors:
            return f"{name} is not valid: " + "; ".join(spec.errors), ""
        loops.confirm(spec)
        return (f"{name} confirmed for its current content (any edit withdraws "
                f"it). /loop run {name} starts it."), ""
    if sub == "log":
        if ctl is None:
            return _NO_CONTROLLER, ""
        path = ctl._log_path(name)
        try:
            rows = path.read_text(encoding="utf-8").splitlines()[-10:]
        except OSError:
            return f"no log for {name} yet ({path}).", ""
        return f"{path}:\n" + "\n".join(rows), ""
    if ctl is None:
        return _NO_CONTROLLER, ""
    if sub == "stop":
        if ctl.loop is None and not ctl.enabled:
            return "no loop is running.", ""
        return ctl.stop("stopped by the user"), ""
    if sub == "resume":
        return ctl.resume(), ""
    if sub == "run":
        if not name:
            return "usage: /loop run <name>", ""
        spec = loops.load(name, workspace)
        if spec is None:
            return f"no loop called {name}. /loop lists them.", ""
        if spec.errors:
            return f"{name} is not valid: " + "; ".join(spec.errors), ""
        if not loops.is_confirmed(spec):
            if ask is None or say is None:
                return (f"{name} is not confirmed for its current content: "
                        f"/loop show {name}, then /loop confirm {name}."), ""
            say(_show(spec))
            answer = str(ask(f"Run loop {name} as shown? [y/N] ") or "")
            if answer.strip().lower() not in ("y", "yes", "j", "ja"):
                return f"loop {name} not started.", ""
            loops.confirm(spec)
            spec = loops.load(name, workspace)
        return ctl.start_loop(spec, engine), ""
    from .scheduler import parse_interval_seconds
    if parse_interval_seconds(sub):
        prompt = " ".join(words[1:])
        spec = loops.ad_hoc(sub, prompt, workspace)
        if spec.errors:
            return "; ".join(spec.errors), ""
        return ctl.start_loop(spec, engine), ""
    return "usage: " + LOOP_USAGE, ""


def controller_hint(ctl: Any) -> str:
    """The status line, or "" -- for a surface's mode line."""
    try:
        return ctl.status_line() if ctl is not None else ""
    except Exception:
        return ""
