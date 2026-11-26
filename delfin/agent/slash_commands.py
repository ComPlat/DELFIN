"""User-defined slash commands as markdown templates.

Drop a ``<name>.md`` file in any of the search locations below and the
chat will accept ``/<name>`` as a slash command that expands to the
file's body. Lighter weight than a full Skill — no Skill frontmatter,
no Skill-runtime, just text substitution.

Search order (later wins on collision):
  1. ``~/.delfin/commands/`` — user-global, survives projects
  2. ``<workspace>/.delfin/commands/`` — project-scoped

DELFIN reads its own namespace only; other agent frameworks' command
directories are not searched (see _search_dirs).

Frontmatter (optional, all keys optional)::

    ---
    description: Short one-line summary shown in /help and the palette.
    argument-hint: <pr-number>      # what the user types after /name
    ---

Inside the body, the following placeholders are substituted before the
text is injected into the chat:

  ``$ARGUMENTS``        full text the user typed after the command
  ``$1``, ``$2``, …    individual whitespace-separated tokens
  ``$@``                alias for ``$ARGUMENTS``

Unknown placeholders are left as-is so a literal ``$foo`` in a template
doesn't accidentally vanish.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable


_NAME_RE = re.compile(r"^[a-zA-Z][a-zA-Z0-9_\-]{0,40}$")
_FRONTMATTER_RE = re.compile(r"^---\s*\n(.*?)\n---\s*\n", re.DOTALL)


@dataclass(frozen=True)
class SlashTemplate:
    name: str
    description: str
    argument_hint: str
    body: str
    source: Path

    @property
    def slash(self) -> str:
        return f"/{self.name}"


def _search_dirs(workspace: Path | str | None) -> list[Path]:
    """Order: user-delfin → project-delfin. DELFIN only reads its own
    namespace — other agent frameworks' command dirs are not loaded."""
    out: list[Path] = [
        Path.home() / ".delfin" / "commands",
    ]
    if workspace is not None:
        out.append(Path(workspace) / ".delfin" / "commands")
    return out


def _parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
    m = _FRONTMATTER_RE.match(text or "")
    if not m:
        return {}, text or ""
    block, body = m.group(1), text[m.end():]
    meta: dict[str, str] = {}
    for line in block.splitlines():
        line = line.strip()
        if not line or line.startswith("#") or ":" not in line:
            continue
        k, _, v = line.partition(":")
        meta[k.strip()] = v.strip().strip("'").strip('"')
    return meta, body


def _load_one(path: Path) -> SlashTemplate | None:
    try:
        text = path.read_text(encoding="utf-8")
    except OSError:
        return None
    if not text.strip():
        return None
    meta, body = _parse_frontmatter(text)
    name = path.stem
    if not _NAME_RE.match(name):
        return None
    description = (meta.get("description") or "").strip()
    argument_hint = (meta.get("argument-hint") or meta.get("argument_hint") or "").strip()
    return SlashTemplate(
        name=name,
        description=description or _first_nonempty_line(body) or name,
        argument_hint=argument_hint,
        body=body.strip(),
        source=path,
    )


def _first_nonempty_line(text: str, *, max_len: int = 100) -> str:
    for line in (text or "").splitlines():
        s = line.strip().lstrip("#").strip()
        if s:
            return s[:max_len]
    return ""


def discover_commands(workspace: Path | str | None = None) -> list[SlashTemplate]:
    """Return every command markdown found across the search layers.

    Project-scoped files override user-global files with the same name
    (last-wins because we iterate in priority order, project after user).
    """
    by_name: dict[str, SlashTemplate] = {}
    for d in _search_dirs(Path(workspace) if workspace else None):
        if not d.is_dir():
            continue
        try:
            paths = sorted(d.glob("*.md"))
        except OSError:
            continue
        for p in paths:
            tpl = _load_one(p)
            if tpl is None:
                continue
            by_name[tpl.name] = tpl
    return sorted(by_name.values(), key=lambda t: t.name)


def get_command(
    name: str, workspace: Path | str | None = None,
) -> SlashTemplate | None:
    target = (name or "").strip().lstrip("/").lower()
    if not target:
        return None
    for tpl in discover_commands(workspace):
        if tpl.name.lower() == target:
            return tpl
    return None


def expand_template(template: str, args: str) -> str:
    """Substitute ``$ARGUMENTS`` / ``$@`` / ``$1`` … in ``template``.

    Unknown ``$<word>`` placeholders are left untouched so users can
    write ``"$0.05"`` inside a template without it being eaten. Numeric
    positional placeholders beyond available tokens become empty
    strings (matching POSIX shell semantics).
    """
    template = template or ""
    args = (args or "").strip()
    tokens = args.split() if args else []

    def _sub(m: re.Match[str]) -> str:
        ref = m.group(1)
        if ref in ("ARGUMENTS", "@"):
            return args
        if ref.isdigit():
            idx = int(ref)
            if idx == 0:
                return args
            return tokens[idx - 1] if idx <= len(tokens) else ""
        return m.group(0)

    pattern = re.compile(r"\$(\{[A-Za-z_@][A-Za-z0-9_@]*\}|[A-Za-z0-9_@]+)")

    def _norm(m: re.Match[str]) -> str:
        raw = m.group(1)
        ref = raw[1:-1] if raw.startswith("{") else raw
        if ref in ("ARGUMENTS", "@"):
            return args
        if ref.isdigit():
            idx = int(ref)
            if idx == 0:
                return args
            return tokens[idx - 1] if idx <= len(tokens) else ""
        return m.group(0)

    return pattern.sub(_norm, template)


def render_command(
    name_or_template: "str | SlashTemplate",
    args: str = "",
    workspace: Path | str | None = None,
) -> str:
    """Resolve a command name (or pre-loaded ``SlashTemplate``), expand
    ``$ARGUMENTS`` and positional placeholders, and return the final
    text ready to send to the agent as a user message."""
    tpl = (
        name_or_template
        if isinstance(name_or_template, SlashTemplate)
        else get_command(name_or_template, workspace)
    )
    if tpl is None:
        raise KeyError(f"unknown slash command: {name_or_template!r}")
    return expand_template(tpl.body, args)


# Built-in subagent slash commands: a direct lever for the USER to delegate
# to a sub-agent (instead of waiting for the agent to auto-delegate). They
# expand to a directive that steers the agent to call the subagent tool with
# the right preset. A user-defined commands/<name>.md OR a skill of the same
# name takes precedence (checked first by the dashboard), so these are
# overridable defaults.
_BUILTIN_SUBAGENT_COMMANDS: dict[str, tuple[str, str]] = {
    "explore": ("explore", "read-only investigation, reports findings"),
    "review": ("code-reviewer", "independent read-only code review"),
    "plan": ("plan", "step-by-step plan, makes no edits"),
    "delegate": ("general-purpose", "self-contained task, inherits permissions"),
}


#: A first token this many characters long is treated as the START OF THE
#: TASK, not as a mistyped preset name. Measured over the 18 registered
#: names: the shortest is 5 characters ("plan" is 4), and a task's first
#: word is usually shorter than a name ("run", "find", "such", "bitte").
#: A long token that is not a name is almost always a German compound or
#: an English noun the user believes is a specialist ("transition-state"),
#: and for those the right answer is still to name a specialist or say
#: none fits — never to refuse the line.
_MIN_NAMEISH = 6

#: How close a token must be to a registered name before it is treated as
#: that name mistyped. Measured over the 18 names with
#: ``difflib.SequenceMatcher``: real typos score 0.875-0.923
#: (``code-rewiewer``→``code-reviewer`` 0.923, ``plann``→``plan`` 0.889,
#: ``ts-find``→``ts-finder`` 0.875) and task words that merely sit near a
#: name score 0.640-0.783 (``extractor`` 0.783, ``reviewer`` 0.762,
#: ``verify`` 0.714, ``planning`` 0.667). 0.85 sits inside that gap. It is
#: HIGHER than the 0.6 a refusal may use, because here the suggestion
#: becomes the destination and the mistyped word is dropped from the
#: briefing: a wrong guess reroutes real work instead of only hinting.
#: ``nearest_specialist`` keeps 0.6 as its default for the refusal case.
_TYPO_CUTOFF = 0.85

#: Inflectional tails an English or German task word wears. Used to tell a
#: TYPO of a preset name from an INFLECTED ordinary word, which a raw
#: similarity ratio cannot do: measured over the 18 names, ``verified``
#: scores 0.875 against ``verifier`` and ``plans`` 0.889 against ``plan``
#: -- both above ``_TYPO_CUTOFF``, both ordinary verbs. `/delegate
#: verified that the SLURM path still works` then arrived as
#: `verifier` + "that the SLURM path still works": the user's first word
#: was eaten and the task was rerouted.
#:
#: THE RULE. Strip one tail and score again. An INFLECTED word is a stem
#: the name also contains, so stripping moves it CLOSER to the name
#: (``plans``->``plan`` 0.889->1.000, ``exploring``->``explore``
#: 0.750->0.923, ``verified`` 0.875->0.857 is level). A TYPO is not an
#: inflected name, so stripping moves it FURTHER or leaves no stem to
#: score (``code-rewiewer`` 0.923->0.833, ``ts-find``/``plann`` no tail).
#: ``_INFLECTION_SLACK`` is the band that counts as "level": the two
#: cases really are the same word. Only a token whose score does NOT rise
#: (beyond the slack) is treated as a mistyped name.
_INFLECTIONS = ("ing", "ed", "es", "s", "en", "er")
_INFLECTION_SLACK = 0.05


def _inflection_stem(token: str) -> str:
    """``token`` without one inflectional tail, or "" if it has none."""
    for tail in _INFLECTIONS:
        if token.endswith(tail) and len(token) - len(tail) >= 3:
            return token[:-len(tail)]
    return ""


def _delegate_fallback(first: str, rest: str, task: str
                       ) -> "tuple[str, str, str]":
    """Who takes ``/delegate <task>``, and what the directive should say.

    Input: the first token, the remainder, and the whole task text.
    Output: ``(subagent_type, task_text, gloss)``.

    Two cases the caller cannot tell apart by itself:

    * The token reads like a mistyped preset name ("code-rewiewer"). Then
      the nearest registered name is used IF one is close, and the gloss
      says it was a suggestion, so a wrong guess is visible in the line
      the model gets. The task stays the user's text.
    * The token is the task's own first word. Then
      ``subagents.suggest_specialist`` reads the whole request against the
      cues the model itself reads. A confident pick is offered WITH ITS
      REASON; anything else — a tie, no match — keeps ``general-purpose``
      and says so.

    The fallback never silently changes the destination for a name the
    user typed correctly, and never picks when it is unsure: a wrong
    specialist that announces itself is a bug someone can see, and a
    wrong specialist that announces nothing is how w41 measured 44 of 50
    dispatches landing on ``explore``.
    """
    from .subagents import nearest_specialist, suggest_specialist

    whole = (task or "").strip()
    near = (nearest_specialist(first, cutoff=_TYPO_CUTOFF)
            if len(first) >= _MIN_NAMEISH else "")
    if near:
        # A name this close can still be an ordinary INFLECTED word that
        # merely resembles one. Distinguish them by re-scoring the stem:
        # stripping a tail that a real typo does not wear moves an
        # inflected word closer to its name and a typo not at all. See
        # ``_INFLECTIONS``. Without this, "/delegate verified that the
        # SLURM path still works" ate "verified" and rerouted the task.
        import difflib
        from .subagents import SUBAGENT_PRESETS as _registry
        raw = difflib.SequenceMatcher(None, first, near).ratio()
        stem = _inflection_stem(first)
        stemmed = (difflib.SequenceMatcher(None, stem, near).ratio()
                   if stem else 0.0)
        if stemmed >= raw - _INFLECTION_SLACK:
            # The first word is the task's own: fall through to scoring.
            near = ""
    if near:
        # The mistyped token is DROPPED from the task: the user was naming
        # who should do the work, not writing part of it, and leaving
        # "code-rewiewer" in the briefing hands the delegate a word that
        # matches nothing. With no remainder the caller's own no-task path
        # fires, and it now names the corrected preset.
        return (near, (rest or "").strip(),
                f"suggested from the name {first!r} you typed, which is "
                f"not registered — the closest is {near!r}")
    name, reason, _ranked = suggest_specialist(whole)
    if name:
        return (name, whole,
                f"suggested by the request itself ({reason}); the user did "
                f"not name a specialist, so confirm this is the right one "
                f"before delegating")
    why = reason or "no specialist's description matches the request"
    return ("general-purpose", whole,
            f"the generalist preset, chosen because no specialist was a "
            f"clear match ({why})")


def builtin_subagent_command(name: str, args: str = "") -> str | None:
    """Expand a built-in subagent slash command to an agent directive.

    Returns the rewritten user message, or None if ``name`` isn't a built-in
    subagent command. ``/delegate`` accepts an optional leading subagent-type
    token (``/delegate explore find X``); otherwise general-purpose is used.
    """
    key = (name or "").strip().lower().lstrip("/")
    spec = _BUILTIN_SUBAGENT_COMMANDS.get(key)
    if spec is None:
        return None
    sa_type, gloss = spec
    task = (args or "").strip()

    # /delegate <type> <task> — let the user pick the preset explicitly.
    if key == "delegate" and task:
        first, _, rest = task.partition(" ")
        try:
            from .subagents import subagent_type_names
            known = set(subagent_type_names() or [])
        except Exception:
            known = set()
        if first in known:
            sa_type, task = first, rest.strip()
            gloss = f"the '{sa_type}' preset"
        elif key == "delegate":
            # The first token is not a preset. Two things can be true and
            # the old code could only say one of them ("general-purpose,
            # always"): the user meant a name and mistyped it, or they
            # wrote the task itself and left the choice to the agent.
            # Distinguish them, and name what was chosen either way.
            sa_type, task, gloss = _delegate_fallback(first, rest, task)

    if not task:
        return (
            f"The user invoked /{key} with no task. Ask them — in one short "
            f"line — what the `{sa_type}` subagent should work on."
        )
    return (
        f"The user asked you to delegate the following to a `{sa_type}` "
        f"subagent ({gloss}). Call subagent(subagent_type='{sa_type}', "
        f"description=<3-7 word label>, prompt=<self-contained briefing with "
        f"the goal, relevant files/paths, and the desired form of the answer>)."
        f" For long or independent work, add background=true and collect the "
        f"result later with subagent_result(sa_id). Then summarise the "
        f"subagent's findings for the user.\n\nTask: {task}"
    )


__all__ = [
    "SlashTemplate",
    "discover_commands",
    "get_command",
    "expand_template",
    "render_command",
    "builtin_subagent_command",
]
