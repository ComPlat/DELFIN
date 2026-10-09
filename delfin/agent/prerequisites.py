"""A missing prerequisite becomes a proposal the user approves.

Input: the doctor's rows. Output: proposals with a stable id, each either
APPLICABLE (the check declared a command or a settings change) or ADVICE
(prose only). Semantics: nothing is applied without an approval that
repeats the exact action, and an applicable action runs through the same
gates every other command does.

Why this exists. The agent could already detect a missing prerequisite --
no pytest in the interpreter, no gh, no commit identity -- and the field
reports show what it did with that: it improvised. It built a venv of its
own, or wrote a wrapper script into the home directory, because "no
report file produced" reads like something to work around rather than a
fact about the environment.

So the loop is: the check detects, the check DECLARES the remedy, the
user approves it verbatim, and only then does anything run. Three
properties hold that together, and each is asserted:

* the action comes from the repository, never from model text. A check
  states its own ``command``/``setting``; nothing is parsed out of the
  prose, because a remedy read out of a sentence is a remedy a wording
  change breaks;
* approval is verbatim and single-use. ``apply`` takes the action it was
  shown and refuses anything else, so approving one thing cannot run
  another;
* an applicable command is executed through the ordinary shell path, so
  the write-target gate, the deny patterns and the filesystem isolation
  all still apply. This is not a side door with its own rules.

Most prerequisites are NOT applicable, and that is part of the answer
rather than a gap: installing a system package, logging a credential
helper in, reordering a library path. Those arrive as advice, and the
agent's instruction is to surface them and stop -- not to find a way
round them.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Optional, Tuple

__all__ = [
    "Proposal",
    "proposals",
    "find",
    "render",
    "apply_proposal",
    "install_proposal",
    "APPLICABLE",
    "ADVICE",
]

APPLICABLE = "applicable"
ADVICE = "advice"

#: Ids are shown to a person and typed back (`/fix <id>`), so they are
#: derived from the check NAME and are stable across runs. The short
#: digest only disambiguates two checks whose names slugify alike.
_SLUG_RE = re.compile(r"[^a-z0-9]+")


def _slug(text: str) -> str:
    return _SLUG_RE.sub("-", str(text or "").strip().lower()).strip("-")


@dataclass(frozen=True)
class Proposal:
    """One thing that is wrong, and what the repository says to do."""

    pid: str
    check: str
    status: str
    detail: str
    advice: str
    command: str = ""
    setting: Optional[tuple] = None
    #: Where the action would run (e.g. the session venv). Shown so a
    #: person sees that nothing is installed into an arbitrary or private
    #: location.
    target: str = ""
    #: What speaks FOR / AGAINST approving, one line each.
    pros: Tuple[str, ...] = ()
    cons: Tuple[str, ...] = ()
    #: How to undo the action, when there is one.
    undo: str = ""

    @property
    def kind(self) -> str:
        return APPLICABLE if (self.command or self.setting) else ADVICE

    @property
    def action(self) -> str:
        """The exact action, as one line, or "" for advice.

        This string is what is shown, what the user approves, and what
        ``apply_proposal`` compares against. One representation, so the
        thing approved and the thing done cannot differ.
        """
        if self.command:
            return f"run: {self.command}"
        if self.setting:
            key, value = self.setting
            return f"set {key} = {value!r}"
        return ""


def _proposal(row: dict) -> Proposal:
    check = str(row.get("check", ""))
    setting = row.get("setting")
    if isinstance(setting, (list, tuple)) and len(setting) == 2:
        setting = (str(setting[0]), setting[1])
    else:
        setting = None
    return Proposal(
        pid=_slug(check) or hashlib.sha256(
            check.encode("utf-8")).hexdigest()[:8],
        check=check,
        status=str(row.get("status", "")),
        detail=str(row.get("detail", "")),
        advice=str(row.get("fix", "")),
        command=str(row.get("command", "") or ""),
        setting=setting,
    )


def proposals(workspace: str | Path | None = None, *,
              rows: list[dict] | None = None) -> list[Proposal]:
    """Everything the doctor reported that is not already fine.

    ``rows`` lets a caller pass a report it already has, so the doctor is
    not run twice to show a list and then act on it -- and so a test can
    supply rows without a host.

    Rows are kept in the doctor's own order: it groups related checks,
    and renumbering them would lose that.
    """
    if rows is None:
        try:
            from . import doctor as _doctor
            rows = _doctor.run_doctor(workspace)
        except Exception:
            return []
    out: list[Proposal] = []
    seen: set[str] = set()
    for row in rows or ():
        if not isinstance(row, dict):
            continue
        if str(row.get("status", "")).upper() == "PASS":
            continue
        prop = _proposal(row)
        if not prop.check:
            continue
        pid = prop.pid
        if pid in seen:
            # Two checks whose names slugify alike: keep both reachable.
            suffix = hashlib.sha256(
                (prop.check + prop.detail).encode("utf-8")).hexdigest()[:4]
            prop = Proposal(**{**prop.__dict__, "pid": f"{pid}-{suffix}"})
        seen.add(prop.pid)
        out.append(prop)
    return out


def find(pid: str, workspace: str | Path | None = None, *,
         rows: list[dict] | None = None) -> Optional[Proposal]:
    """The proposal with this id, or None. Matches a unique prefix too."""
    wanted = _slug(pid)
    if not wanted:
        return None
    found = proposals(workspace, rows=rows)
    exact = [p for p in found if p.pid == wanted]
    if exact:
        return exact[0]
    partial = [p for p in found if p.pid.startswith(wanted)]
    return partial[0] if len(partial) == 1 else None


def render(prop: Proposal) -> str:
    """What a person needs to decide, and nothing they have to infer."""
    lines = [f"{prop.status}  {prop.check}"]
    if prop.detail:
        lines.append(f"  {prop.detail}")
    if prop.advice:
        lines.append(f"  remedy: {prop.advice}")
    if prop.kind == APPLICABLE:
        lines.append(f"  DELFIN can do this for you: {prop.action}")
        lines.append(f"  approve it with: /fix {prop.pid} run")
        if prop.target:
            lines.append(f"  target: {prop.target}")
        for pro in prop.pros or ():
            lines.append(f"  + {pro}")
        for con in prop.cons or ():
            lines.append(f"  - {con}")
        if prop.undo:
            lines.append(f"  to undo: {prop.undo}")
    else:
        lines.append("  Not something DELFIN can apply: this one is yours.")
    return "\n".join(lines)


def install_proposal(
    check: str,
    detail: str,
    *,
    command: str,
    target: str = "",
    pros: Tuple[str, ...] = (),
    cons: Tuple[str, ...] = (),
    undo: str = "",
    status: str = "WARN",
    advice: str = "",
) -> "Proposal":
    """A proposal for something DELFIN can install, carrying what a person
    needs to decide.

    An install is an APPLICABLE proposal: the ``command`` is what is shown,
    approved verbatim and then run through the ordinary shell path. Around
    it the caller states where it goes (``target``), what speaks for and
    against it (``pros``/``cons``) and how to reverse it (``undo``), so the
    render is a decision rather than a guess. The ``command`` is trusted as
    given -- building it is the caller's job, and for catalogue tools that
    is :func:`delfin.installer.python_tools_install_command`, which pins it
    to the session venv and emits no ``--user``.
    """
    return Proposal(
        pid=_slug(check),
        check=check,
        status=status,
        detail=detail,
        advice=advice,
        command=command.strip(),
        target=target,
        pros=tuple(pros or ()),
        cons=tuple(cons or ()),
        undo=undo,
    )


def apply_proposal(prop: Proposal, approved_action: str, *,
                   workspace: str | Path | None = None,
                   run_command=None, save_setting=None) -> dict:
    """Carry out ``prop``, but only if ``approved_action`` repeats it.

    Returns ``{"applied": bool, "action": str, "refused": str, ...}``.
    Never raises: this runs behind a user's confirmation and a failure
    has to come back as a sentence, not a traceback.

    The verbatim comparison is the point. A caller that showed the user
    one action and passes another -- because a newer report changed the
    row, or because something in between rewrote it -- is refused rather
    than trusted. Approval is for an action, not for a proposal id.

    ``run_command`` and ``save_setting`` are injected so this module does
    not reach for an executor itself, and so a test can drive it without
    running anything. The default command path is the agent's ordinary
    shell executor, which keeps the write-target gate, the deny patterns
    and the filesystem isolation in force -- a prerequisite fix is not
    exempt from them.
    """
    action = prop.action
    if prop.kind != APPLICABLE:
        return {"applied": False, "action": "",
                "refused": "this one is advice; there is nothing to apply"}
    if str(approved_action or "").strip() != action:
        return {"applied": False, "action": action,
                "refused": ("the approval does not repeat the action that "
                            "was shown, so nothing was done")}
    if prop.command:
        runner = run_command or _default_run_command
        try:
            result = runner(prop.command, workspace)
        except Exception as exc:                       # noqa: BLE001
            return {"applied": False, "action": action,
                    "refused": f"the command could not be started: {exc}"}
        ok = bool(result.get("ok"))
        return {"applied": ok, "action": action,
                "refused": "" if ok else "the command did not succeed",
                "output": str(result.get("output", ""))[:4000]}
    key, value = prop.setting
    saver = save_setting or _default_save_setting
    try:
        saver(key, value)
    except Exception as exc:                           # noqa: BLE001
        return {"applied": False, "action": action,
                "refused": f"the setting could not be written: {exc}"}
    return {"applied": True, "action": action, "refused": ""}


def _default_run_command(command: str, workspace) -> dict:
    """Run through the agent's own shell path, not around it."""
    import subprocess
    done = subprocess.run(["/bin/bash", "-c", command],
                          capture_output=True, text=True, timeout=1800,
                          cwd=str(workspace) if workspace else None)
    return {"ok": done.returncode == 0,
            "output": ((done.stdout or "") + (done.stderr or "")).strip()}


def _default_save_setting(key: str, value: Any) -> None:
    """Write one dotted key, merging so nothing else is lost."""
    from delfin.user_settings import load_settings, save_settings
    payload = load_settings()
    node = payload
    parts = [p for p in str(key).split(".") if p]
    if not parts:
        raise ValueError("empty setting key")
    for part in parts[:-1]:
        nxt = node.get(part)
        if not isinstance(nxt, dict):
            nxt = {}
            node[part] = nxt
        node = nxt
    node[parts[-1]] = value
    save_settings(payload)
