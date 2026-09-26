"""Answering a CHOICE question (ask_user_question) from outside the pane.

A terminal session's ask_user_question is published to the terminal
confirmations directory. ``approvals`` used to know only approve/deny,
which is no answer at all for a numbered choice: the operator could only
deny and smuggle the pick into the refusal reason. This module resolves
what the operator typed (an option NUMBER or an option LABEL) against
the question's own options, and writes the answer beside the request
the same way file_confirm writes an approve/deny -- same directory,
same tmp-then-replace discipline, same owner-only mode -- so the
session's checked reader can pick it up.

Nothing here DECIDES anything: the choice comes from the operator, the
options come from the question, and the reading side stays in
file_confirm / terminal_confirm (security code).
"""

from __future__ import annotations

import json
import os
from pathlib import Path
from typing import Any, Optional


class ChoiceError(Exception):
    """The choice cannot be applied to this request. Message is the CLI
    text: it says what was wrong and what the valid picks are."""


def options_of(record: dict) -> list[str]:
    """The labels a choice question offers, in order. Empty when the
    record is not a choice question (or carries no usable options).

    A published ask question stores the ask_user_question payload under
    ``payload``; ``options`` is a list of {label, ...} dicts. Anything
    malformed is no options rather than a raise: the caller refuses the
    answer either way.
    """
    payload = record.get("payload")
    if not isinstance(payload, dict):
        return []
    raw = payload.get("options")
    if not isinstance(raw, list):
        return []
    labels = []
    for opt in raw:
        if isinstance(opt, dict) and isinstance(opt.get("label"), str) \
                and opt["label"]:
            labels.append(opt["label"])
    return labels


def is_choice_question(record: dict) -> bool:
    """True when this record is an ask_user_question (kind ``ask``) with
    at least two usable options. Anything else cannot be answered with
    a choice and is refused by the CLI."""
    if str(record.get("kind") or "") != "ask":
        return False
    return len(options_of(record)) >= 2


def is_multi_select(record: dict) -> bool:
    payload = record.get("payload")
    return bool(isinstance(payload, dict) and payload.get("multiSelect"))


def resolve_choice(record: dict, choice: str) -> list[str]:
    """One operator-typed choice -> the selected label list.

    Accepts the option's number (``2``) or its exact label. Case is kept
    for the label but matched case-insensitively, so the operator does
    not have to retype capitalisation the model chose. Multiple picks
    (comma-separated numbers or labels) are allowed only when the
    question allows multiSelect; duplicates collapse.

    Raises ChoiceError with the question's options listed when the pick
    matches nothing -- a bare "cannot answer" would send the operator
    back to `show` to count options by hand.
    """
    labels = options_of(record)
    multi = is_multi_select(record)
    parts = [p.strip() for p in str(choice or "").split(",")] if multi \
        else [str(choice or "").strip()]
    picks: list[str] = []
    for part in parts:
        if not part:
            continue
        # A number is the option's position, as the pane shows it.
        if part.isdigit():
            n = int(part)
            if not (1 <= n <= len(labels)):
                raise ChoiceError(
                    f"option number {n} is not one of the "
                    f"{len(labels)} options: "
                    + ", ".join(f"{i + 1}. {l}" for i, l in enumerate(labels)))
            picks.append(labels[n - 1])
            continue
        lower = part.lower()
        hit = [l for l in labels if l.lower() == lower]
        if not hit:
            raise ChoiceError(
                f"{part!r} is not one of the options: "
                + ", ".join(f"{i + 1}. {l}" for i, l in enumerate(labels)))
        picks.append(hit[0])
    if not picks:
        raise ChoiceError(
            "no option chosen; pick one of: "
            + ", ".join(f"{i + 1}. {l}" for i, l in enumerate(labels)))
    # Deduplicate, keep the chosen order.
    out: list[str] = []
    for p in picks:
        if p not in out:
            out.append(p)
    if not multi and len(out) > 1:
        raise ChoiceError(
            "the question allows one pick only (multiSelect is off)")
    return out


def write_choice_answer(record: dict, answers: list[str], *,
                        room: Optional[Path] = None,
                        by: str = "") -> Path:
    """Write the chosen labels beside the request, file_confirm-style.

    Same shape and same discipline as ``file_confirm.answer``: written
    under a temporary name, mode 0600, then moved into place, so the
    checked reader never sees half an answer. The decision key is
    ``"choose"`` -- distinct from APPROVE/DENY on purpose, so an old
    reader treats it as unrecognised (not an answer) rather than
    approving something.
    """
    from . import file_confirm as _fc
    room = Path(room) if room else _fc.requests_dir()
    request_id = str(record.get("id") or "")
    request_path = room / f"{request_id}.request.json"
    if not request_path.is_file():
        raise ChoiceError(f"no request file for id {request_id!r}")
    payload = {"id": request_id,
               "decision": "choose",
               "answers": [str(a) for a in answers],
               "by": str(by or "")}
    path = room / f"{request_id}.answer.json"
    tmp = room / f".{request_id}.answer.partial"
    _fc._own_dir(room)
    tmp.write_text(json.dumps(payload), encoding="utf-8")
    os.chmod(tmp, 0o600)
    tmp.replace(path)
    return path


def find_pending(request_id: str) -> Optional[dict]:
    """The pending record with this id, headless room first, then the
    terminal rooms. None when nothing waits under that id.

    Read-only; the rooms are the two places ``approvals ls`` already
    lists, so `answer` accepts exactly the ids `ls` shows.
    """
    from . import file_confirm as _fc
    for row in _fc.pending():
        if str(row.get("id") or "") == str(request_id):
            return row
    from . import terminal_confirm as _tc
    for row in _tc.pending_at_terminals():
        if str(row.get("id") or "") == str(request_id):
            return row
    return None


def answer(request_id: str, choice: str, *, by: str = "") -> list[str]:
    """End-to-end for the CLI: find the question, check it is one,
    resolve the pick, write the answer. Returns the chosen labels.

    Raises ChoiceError (never a bare return code) so the CLI can print
    the message and exit 2 without re-deriving what went wrong.
    """
    record = find_pending(request_id)
    if record is None:
        raise ChoiceError(f"nothing waiting with id {request_id!r}")
    if not is_choice_question(record):
        raise ChoiceError(
            f"request {request_id!r} is not a choice question -- "
            "use approve or deny for it")
    picks = resolve_choice(record, choice)
    write_choice_answer(record, picks, by=by)
    return picks


def preview_line(record: dict, width: int = 80) -> str:
    """One line saying what a pending request is about (Phase 3).

    The subject is the command, else the path, else the question of a
    choice request -- the three things an operator decides on. Clipped
    to *width* with an ellipsis so a long command cannot push the rest
    of the listing off the pane.
    """
    subject = str(record.get("command") or "").strip()
    if not subject:
        subject = str(record.get("path") or "").strip()
    if not subject:
        payload = record.get("payload")
        if isinstance(payload, dict):
            subject = str(payload.get("question") or "").strip()
    if not subject:
        subject = str(record.get("preview") or "").strip()
    subject = " ".join(subject.split())
    if len(subject) > width:
        subject = subject[:max(0, width - 1)] + "…"
    return subject
