"""Branch commands for the Pi's surfaces: the words, one spelling.

The REPL (``/branches [view]``, ``/branches test <n|name>``,
``/branches pull``, ``/branches cleanup``) and the
messenger gateway (``branches``, ``test <n>``, ``pull``, ``cleanup``,
``ok <n>``) answer the same questions with the same words. This module is
that shared layer: it formats ``branch_ops`` answers into a few lines and
resolves what a number means. The surfaces stay thin -- a handler calls
one function here and prints the text.

Numbering, one rule for every surface: a number names a row of the
``mine`` view as :func:`rows_text` printed it (1-based, ``main`` is 1).
A typed branch name always works too, whatever was listed last.

Deletion, the dashboard's rules: ``cleanup`` lists what
``branch_cleanup`` offers (merged branches and stale worktrees only, so
nothing with its own commits is shown by default); ``cleanup ok`` (TUI)
or ``ok <n>`` (phone) deletes through ``branch_cleanup.delete``, which
re-reads its rules at deletion time. A branch with commits not in the
base is not offered; naming it directly answers with the count and needs
the second confirmation (``confirm_unmerged=True``), which the phone
gets by answering ``ok <n>`` again after the count was shown.

Nothing here switches the shared checkout.
"""

from __future__ import annotations

from pathlib import Path

from delfin.agent import branch_ops


def _workspace(workspace) -> Path:
    return Path(workspace).expanduser()


def rows_text(workspace, view: str = "mine", text: str = "") -> str:
    """The branches of *view*, numbered: ``1. main``, ``2. feat (+2)``,
    one line each, plus the PR note a row carries."""
    answer = branch_ops.rows(_workspace(workspace), view, text)
    if not answer.get("ok"):
        return answer.get("text") or "not a git checkout"
    rows = answer.get("rows") or []
    if not rows:
        return "no branches in this view"
    lines = []
    for i, row in enumerate(rows, 1):
        parts = [f"{i}. {row['ref']}"]
        if row.get("ahead"):
            parts.append(f"+{row['ahead']}")
        if row.get("pr"):
            parts.append(f"PR #{row['pr']}" +
                         (f" {row['ci']}" if row.get("ci") else ""))
        if row.get("current"):
            parts.append("checked out here")
        if row.get("merged") and view == "merged":
            parts.append("merged")
        lines.append("  ".join(parts))
    return "\n".join(lines)


def resolve(workspace, token: str) -> str:
    """What *token* names: a number is a row of the ``mine`` view (the
    numbering :func:`rows_text` prints); anything else is the branch
    name as typed. Empty for a name that resolves to nothing."""
    token = (token or "").strip().lstrip("#")
    if not token:
        return ""
    if token.isdigit():
        n = int(token)
        answer = branch_ops.rows(_workspace(workspace), "mine")
        rows = (answer.get("rows") or []) if answer.get("ok") else []
        if not 1 <= n <= len(rows):
            return ""
        return str(rows[n - 1]["ref"])
    return token


def test_text(workspace, ref_or_number: str, *, keep: bool = False) -> str:
    """The answer to ``/branches test <n|name>``: the run's counts and
    the first
    failure, in a few lines."""
    ref = resolve(workspace, ref_or_number)
    if not ref:
        return f"no branch {ref_or_number!r} (a number names a row of " \
               "branches; see /branches)"
    result = branch_ops.test_branch(_workspace(workspace), ref, keep=keep)
    header = result.get("text") or ("tested" if result.get("ok")
                                    else "not tested")
    return header


def pull_text(workspace) -> str:
    """The answer to ``/branches pull``: main updated fast-forward-only, or why
    not."""
    result = branch_ops.pull_main(_workspace(workspace))
    return result.get("text") or ("updated" if result.get("ok")
                                  else "not pulled")


def cleanup_text(workspace) -> str:
    """The answer to ``/branches cleanup``: what may be deleted,
    numbered, and the
    word that confirms."""
    candidates, _kept = branch_ops.cleanup_candidates(_workspace(workspace))
    if candidates is None:
        return "not a git checkout"
    if not candidates:
        return "nothing to clean up: no merged branches, no stale worktrees."
    lines = []
    for i, c in enumerate(candidates, 1):
        lines.append(f"{i}. {c['name']}" +
                     (f" ({c['reason']})" if c.get("reason") else ""))
    lines.append("confirm with: /branches cleanup ok <number>")
    return "\n".join(lines)


def cleanup_names(workspace) -> list[str]:
    """The names :func:`cleanup_text` numbered, in its order."""
    candidates, _kept = branch_ops.cleanup_candidates(_workspace(workspace))
    if not candidates:
        return []
    return [str(c["name"]) for c in candidates]


def cleanup_confirm_text(workspace, token: str = "", *,
                         include_unmerged: bool = False) -> tuple[str, bool]:
    """The answer to ``cleanup ok <n|name>`` (empty token = everything
    offered). ``(text, deleted_something)``. Deletion goes through
    ``branch_cleanup.delete``: its rules are re-read at deletion, so a
    branch that moved or became the checkout's branch since the listing
    is kept with a reason."""
    names = cleanup_names(workspace)
    if token.strip():
        ref = resolve(workspace, token)
        if not ref:
            return (f"no branch {token!r} (a number names a row of the "
                    "cleanup listing)", False)
        chosen = [ref]
    else:
        chosen = names
        if not chosen:
            return ("nothing to clean up: no merged branches, no stale "
                    "worktrees.", False)
    result = branch_ops.cleanup_delete(
        _workspace(workspace), chosen, confirm=True,
        confirm_unmerged=True if include_unmerged else None)
    return (result.get("text") or ("deleted" if result.get("ok")
                                   else "not deleted"),
            bool(result.get("ok")))


def unmerged_note(workspace, name: str) -> str:
    """Why *name* is not in the cleanup listing, when it has commits of
    its own: the count, and what confirms its deletion."""
    from delfin.dashboard import branch_cleanup as BC
    repo = _workspace(workspace)
    if not BC.BL.ref_kind(repo, name):
        return ""
    n = BC.unmerged_count(repo, name) or 0
    if not n:
        return ""
    return (f"{name} has {n} commit(s) not in the base; it is not offered "
            "for cleanup. Delete it from a terminal or the dashboard "
            "(a second confirmation is asked there).")


__all__ = ["rows_text", "resolve", "test_text", "pull_text", "cleanup_text",
           "cleanup_names", "cleanup_confirm_text", "unmerged_note"]
