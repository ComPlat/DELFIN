"""Ground tool calls against the actual world before they run.

Wave 13, package T4. Models call tools that do not exist, edit paths that
do not exist, and retry a failing step many times — each costing an
operator decision or a turn. This module answers, before a call is
dispatched, "is this call grounded?" with a structured, always-optional
hint (``None`` when nothing is wrong).

It is deliberately PURE: read-only, no network, and free of any import
from api_client or the private sandbox, exactly as ``process_budget.py``
and ``principles_guard.py`` are, so the security core stays untouched. The
known tool names are INJECTED by the caller, because a checker that had to
import the protected {@link dispatch} would itself be protected code.

Two kinds of grounding are answered (phase 2):

* an unknown tool name — difflib near-miss over the injected known list,
  and
* a missing path in a file tool ('edit_file', 'write_file', 'read_file',
  ...) whose parent directory exists on disk — difflib over the sibling
  names that are actually there, mirroring the established
  ``office._near_miss`` behaviour (a case-miss alone once cost a real
  round).

The path check is bounded: it lists only the immediate parent directory of
the requested path (which must exist), never the whole workspace tree, so
a caller cannot turn it into a directory walk.

The third kind — a failure budget over repeated identical failing calls —
lives in :mod:`delfin.agent.action_grounding.FailureBudget` and is covered
by phase 3.
"""
from __future__ import annotations

import difflib
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional

# Tools whose primary argument is a file path. The path grounding only
# applies to these, so a harmless "bash" call with a command string is
# never misread as a missing file.
_FILE_TOOLS = frozenset({
    "read_file", "write_file", "edit_file", "multi_edit",
    "read_document", "grep_file", "list_files",
})

# Path-argument keys checked, in preference order. edit_file/write_file
# take "path"; some wrappers name it differently.
_PATH_KEYS = ("path", "file", "src", "folder", "directory")

# difflib cutoff for a "did you mean" path suggestion. Mirrors the
# established office._near_miss value; below this a suggestion is more
# noise than help.
_PATH_CUTOFF = 0.75
# difflib cutoff for a "closest tool" suggestion. The existing API-client
# hint (api_client._DOC_TOOLS dispatch) uses 0.5; keep the same bar so a
# weak miss still gets one useful suggestion rather than a wall of names.
_TOOL_CUTOFF = 0.5


@dataclass(frozen=True)
class GroundingHint:
    """A structured reason a tool call is not grounded.

    ``kind`` is one of ``"unknown_tool"``, ``"missing_path"``.
    ``closest`` is the read-only suggestion (the nearest known tool name,
    or the nearest sibling file name) — the caller decides whether to act
    on it. ``message`` is human-readable English.
    """

    kind: str
    message: str
    closest: Optional[str] = None


def _arg_path(args: object) -> Optional[str]:
    """The path a file tool would act on, if its arguments name one."""
    if not isinstance(args, dict):
        return None
    for key in _PATH_KEYS:
        val = args.get(key)
        if isinstance(val, str) and val.strip():
            return val.strip()
    return None


def _closest_name(want: str, names: Iterable[str],
                  cutoff: float, n: int) -> Optional[str]:
    """The nearest case-insensitive match in *names*, or None below cutoff."""
    lowered = want.lower()
    pool = {nm: nm.lower() for nm in names}
    exact = [nm for nm, low in pool.items() if low == lowered]
    if exact:
        return exact[0]
    close = difflib.get_close_matches(
        lowered, sorted(pool.values()), n=n, cutoff=cutoff)
    if not close:
        return None
    return next(nm for nm, low in pool.items() if low == close[0])


def _suggest_tool(tool: str, known: Iterable[str]) -> Optional[str]:
    """The nearest known tool name for an unknown *tool*, or None."""
    return _closest_name(tool, known, cutoff=_TOOL_CUTOFF, n=3)


def _suggest_path(path: Path, workspace: Path) -> Optional[str]:
    """The nearest existing sibling name for a missing *path*, or None.

    Only the immediate parent is listed and only when it exists on disk,
    so the check costs at most one directory listing and never walks the
    workspace. A missing parent means there is nothing to diff against.
    """
    parent = path.parent
    if not parent.is_absolute():
        parent = workspace / parent
    try:
        siblings = [c.name for c in parent.iterdir() if c.is_file()]
    except OSError:
        return None
    return _closest_name(path.name, siblings, cutoff=_PATH_CUTOFF, n=1)


def check(tool: object, args: dict, workspace: str | Path,
          known_tools: Optional[Iterable[str]] = None,
          ):
    """Ground one tool call; return a :class:`GroundingHint` or None.

    Read-only and safe to call on every dispatch. ``known_tools`` is the
    caller-injected list of tool names the agent may actually call — omit
    it (or pass ``None``) to disable the unknown-tool check, which cannot
    rule a name out without that list. ``workspace`` bounds where the path
    check looks.
    """
    wspace = Path(workspace).expanduser() if workspace else Path(".")
    if not isinstance(tool, str) or not tool.strip():
        return None
    tool = tool.strip()
    args = args if isinstance(args, dict) else {}

    # 1. Unknown tool name — only when the caller can name the known set.
    if known_tools is not None:
        known = [t for t in known_tools if isinstance(t, str) and t.strip()]
        if known and tool not in known:
            closest = _suggest_tool(tool, known)
            if closest:
                return GroundingHint(
                    kind="unknown_tool",
                    message=(
                        f"no such tool: {tool!r}; closest: {closest!r}. "
                        "Call that tool by its exact name, or pick one from "
                        f"the known tools: {', '.join(sorted(known[:8]))}"
                        + (" …" if len(known) > 8 else "") + "."
                    ),
                    closest=closest,
                )
            return None

    # 2. Missing path in a file tool, with a near-miss on disk.
    if tool in _FILE_TOOLS:
        path_str = _arg_path(args)
        if path_str:
            path = Path(path_str)
            # Resolve relative paths against the given workspace, not the
            # process cwd — the caller told us where the tree lives.
            full = path if path.is_absolute() else wspace / path
            if not full.exists():
                closest = _suggest_path(path, wspace)
                if closest:
                    return GroundingHint(
                        kind="missing_path",
                        message=(
                            f"no such file: {path_str!r}; did you mean "
                            f"{closest!r}? The folder has a name that close. "
                            "Check the exact spelling before acting on it."
                        ),
                        closest=closest,
                    )

    return None


# ---------------------------------------------------------------------------
# Phase 3: repeated-failure budget — a "stop and ask" when the SAME call
# (tool * normalized args) has already failed several times in one task.
#
# The existing per-call guard (api_client.py:11506-11543 -> failure_log.py)
# keys on a RAW [:80] command/error string, so two spellings of the same
# failing call that differ only in argument order or whitespace land in
# DIFFERENT buckets and reset each other's counter. This module is the pure,
# reusable counterpart: it keys on ``(tool, normalize_args(args))``, so a
# reordered / re-spaced / differently-cased repeat of a failing call is
# recognised as the SAME call. It owns no I/O and no import from api_client
# or the sandbox — it mirrors ``process_budget.py``'s shape (a frozen
# limits dataclass + a pure decision object), stays unit-testable, and the
# protected wiring patch feeds this module's *signature* into the existing
# failure_log hook instead of adding a parallel note.
# ---------------------------------------------------------------------------

def normalize_args(args: object) -> str:
    """A canonical key for a call's arguments.

    So reordered / re-spaced / differently-cased spellings of the same call
    collapse into ONE signature. Non-dicts fall back to ``repr`` (collapsed);
    dict keys are sorted, values lowered and whitespace-collapsed.
    """
    if isinstance(args, dict):
        parts = []
        for key in sorted(args):
            val = args[key]
            parts.append(f"{key}={_canonical(val)}")
        return " ".join(parts)
    return _canonical(args)


def _canonical(value: object) -> str:
    """Lowercase and collapse runs of whitespace — the only normalisation
    that is safe for paths AND command strings alike."""
    text = repr(value) if not isinstance(value, (str, int, float, bool)) \
        else str(value)
    return " ".join(text.lower().split())


@dataclass(frozen=True)
class FailureBudgetLimits:
    """How many identical failing calls before we stop and ask.

    Mirrors ``process_budget.BudgetLimits``: plain frozen dataclass of the
    knobs the pure decision function reads. ``identical_fail_limit`` is the
    count of consecutive identical ``(tool, normalized args)`` *failures*
    at which a subsequent identical call yields a "stop and ask" hint.
    """
    identical_fail_limit: int = 3


class FailureBudget:
    """Per-task gate: the same failing call a third time says "stop and ask".

    Pure and stateful but side-effect-free — all state is an in-memory
    counter keyed on ``(tool, normalize_args(args))``. ``retrospective``
    reports the failing signature after the fact; ``repeat_hint`` is called
    BEFORE a call runs and returns a :class:`GroundingHint` when the
    identical failing call has already fired ``identical_fail_limit`` times.
    A success for a signature resets its streak.
    """

    def __init__(self, limits: Optional[FailureBudgetLimits] = None):
        self._limits = limits or FailureBudgetLimits()
        self._fails: Counter = Counter()          # (tool, sig) -> fail count
        self._streaks: dict[tuple[str, str], int] = {}  # consecutive fails

    def _key(self, tool: object, args: object) -> tuple[str, str]:
        t = str(tool).strip()
        return t, normalize_args(args)

    def record(self, tool: object, args: object, ok: bool) -> None:
        """Register the outcome of one call. A success breaks the streak."""
        key = self._key(tool, args)
        if ok:
            self._fails[key] = 0
            self._streaks[key] = 0
            return
        self._fails[key] += 1
        self._streaks[key] = self._streaks.get(key, 0) + 1

    def repeat_hint(self, tool: object, args: object
                    ) -> Optional[GroundingHint]:
        """A "stop and ask" hint when this failing call has already failed
        ``identical_fail_limit`` times; None otherwise (run the call)."""
        key = self._key(tool, args)
        # Use consecutive-failure count, not total: a call that succeeded
        # in between is not a stuck loop.
        reached = self._streaks.get(key, 0) >= self._limits.identical_fail_limit
        if not reached:
            return None
        return GroundingHint(
            kind="repeated_failure",
            message=(
                f"this exact call ({key[0]}, normalized args) has failed "
                f"{self._limits.identical_fail_limit}+ times in a row in this "
                "task. Repeating it will keep failing. Change approach — pick "
                "a different argument, a different tool, or STOP and ask the "
                "user for help."
            ),
        )

    def failure_signature(self) -> Optional[str]:
        """The normalized ``(tool, args)`` string of the most-repeated failing
        call, for feeding the existing failure_log hook; None when none."""
        # A success sets a key's count to 0 — never report a reset (0-count)
        # key as "most repeated", it would mask a real repeating failure.
        positives = [(k, c) for k, c in self._fails.items() if c > 0]
        if not positives:
            return None
        (tool, sig) = max(positives, key=lambda kc: kc[1])[0]
        return f"{tool}|{sig}"
