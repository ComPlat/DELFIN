"""swarm role policy: role deny patterns as DATA (package SWP).

The role policies are data here, so a swarm's members get the same
enforced practice the operators' own waves run under (see
docs/SWARM_PRACTICE.md). A role policy is a list of deny RE patterns; a
member's effective policy is COMMON + role extras + NO_PUSH, composed at
runtime.

Nothing here hard-codes a home path, user name or host: the "other
people's homes" and the member's own ``~/.delfin`` reads are derived from
``Path.home()`` (or the caller-supplied home) at call time.

These are REPO-tier deny patterns (the member's own
``<worktree>/.delfin/settings.json`` is that worktree's repo settings
file, and ``kit_settings._merge`` makes deny_patterns tighten-only, so a
member cannot loosen its own policy).
"""

from __future__ import annotations

import re
from pathlib import Path


class PolicyError(ValueError):
    """A role policy cannot be built, with the reason in prose."""


# ---------------------------------------------------------------------------
# Member-wide COMMON patterns. Every role refuses all of these.
#
# The vocabulary mirrors the operators' proven per-worktree deny file
# (the write-tighten model) translated into pathless/runtime forms:
# nothing here is tied to a concrete machine, user or home.
# ---------------------------------------------------------------------------

COMMON_DENY_PATTERNS: tuple[str, ...] = (
    # delete nothing (mv to .gate/ instead)
    r"(?:^|[;|&`(]|\s)rm\s",
    r"\bunlink\b",
    r"rmtree",
    r"\bos\.remove\(",
    r"\.unlink\(",
    r"\bshred\b",
    # repository hygiene
    r"git\s+worktree\s+add\b",
    r"git\s+fetch\b",
    r"git\s+pull\b",
    # patch application writes files the write-scope check cannot see:
    # only the dry run (--check, taken as a whole flag) is allowed
    r"git\s+apply\b(?![^;|&]*--check(?:\s|$|[;&|]))",
    r"git\s+am\b",
    r"(?:^|[;&`(>\{\}]\s*)patch\s",
    # tests only through the project's test entry; a bare `pytest` is a
    # command only when it stands on its own (not inside `pytest.ini` or
    # `pytest-results.txt`)
    r"python\s+-m\s+pytest\b",
    r"python3\s+-m\s+pytest\b",
    r"(?:^|[;&`(>\s])pytest(?=\s|$|[;&|])",
    # no "||" chains, and no unicode / whitespace look-alikes that fake a
    # separator (FULLWIDTH SEMICOLON, NBSP, ZERO WIDTH SPACE, ...)
    r"\|\|",
    r"[\u00a0\u200b\u202f\uFF1B\uFD3E\uFD3F\u2212\u2011]",
    # cluster control is DELFIN's job list, never the raw commands
    r"(?:^|[;|&`(]\s*|\bxargs\s+|\bexec\s+)(?:squeue|sbatch|scancel|scontrol)\b",
    # secret env reads
    r"(?:^|[;|&`(]\s*)(?:env|printenv|set)\s*(?:$|[;|&>])",
    r"\bprintenv\b",
    r"\$\{?\w*(?:KEY|TOKEN|SECRET|PASSW)\w*",
    r"environ\b[^;|&]*(?:KEY|TOKEN|SECRET|PASSW)",
    r"getenv\([^)]*(?:KEY|TOKEN|SECRET|PASSW)",
    r"/proc/[^\s]*/environ",
    # secret stores (pathless so no home is hard-coded)
    r"\.ssh\b",
    r"\.env(?:\b|\.)",
    r"\.pem\b",
    r"(?:\.|/)(?:id_rsa|id_ed25519|id_dsa)\b",
    # own policy files and credential stores
    r"\.delfin/(?:settings|trusted_workspaces|mcp_servers)",
    # "credentials" is refused everywhere except a Python module named
    # credentials.py (code, not secrets)
    r"credentials(?!\.py\b)",
)

#: NO push for EVERY role: the push is done by the user/bot, never a member.
NO_PUSH_PATTERNS: tuple[str, ...] = (r"git\s+push\b",)

#: A reviewer (and QS) commits no product code -- only its adversarial
#: tests. `git add`/`git commit` that touch delfin/ (or stage everything)
#: are refused; staging ``tests/test_*_adversarial.py`` stays allowed.
REVIEWER_NO_CODE_PATTERNS: tuple[str, ...] = (
    r"git\s+commit\b[^;|&]*\s(?:--all\b|-[b-zA-Z]*a[a-zA-Z]*\b)",
    r"git\s+(?:add|commit\s+-a|commit\s+--all)\b[^;|&]*(?:\sdelfin/|\s-A\b|\s--all\b|\s\.(?:\s|$))",
)

#: A QS runs only the QS test entry, never a product test run.
QS_TEST_SCOPE_PATTERNS: tuple[str, ...] = (
    r"gate\s+tests/[^;\|&]*\.py",
)


#: Role -> extra deny patterns. The DATA a caller or a test reads. The
#: effective policy of a role = COMMON + home-derived + NO_PUSH + extras.
ROLE_DENY: dict[str, tuple[str, ...]] = {
    "builder": (),
    "reviewer": REVIEWER_NO_CODE_PATTERNS,
    # QS may commit neither product code nor (like a reviewer) touch it.
    # Its test running is restricted to the QS entry.
    "qs": REVIEWER_NO_CODE_PATTERNS + QS_TEST_SCOPE_PATTERNS,
}


def _home_owner(home: Path) -> str:
    """The first path component after a ``home`` mount, as the owner key.

    ``/home/<user>``                 -> ``<user>``
    ``/home/alice/team/one``          -> ``alice``   (nested home tree)
    ``/pfs/data6/home/alice/team/one``  -> ``alice``
    When no ``home`` component exists, the final part is the fallback.
    """
    parts = home.parts
    for i, part in enumerate(parts):
        if part == "home" and i + 1 < len(parts):
            return parts[i + 1]
    return home.name


def home_deny_patterns(home: Path | str | None = None) -> list[str]:
    """The home-derived deny patterns, built at runtime from ``home``.

    Returns the patterns that refuse other users' homes (absolute and via
    ``~user`` expansion) while protecting the member's own, plus reads of
    the member's own ``~/.delfin``. No concrete path is a literal here.
    """
    resolved = Path(home) if home is not None else Path.home()
    resolved = resolved.expanduser()
    owner = _home_owner(resolved)
    esc = re.escape(owner)
    return [
        # tilde / $HOME forms toward the member's own .delfin are refused,
        # like the absolute form below (the user's own state is off-limits)
        r"(?:~|\$HOME|\$\{HOME\})/\.delfin\b",
        # walking OUT of the member's own home via '..' lands in another
        # user's territory by definition -- refuse it.  This never blocks a
        # legitimate command: a member's own files live *inside* ~, and any
        # '../' past the home root is someone else's.
        r"(?:~|\$HOME|\$\{HOME\})/(?:\.\./)+",
        # the member's own home .delfin, derived (no literal user/path)
        re.escape(str(resolved)) + r"/\.delfin\b",
        # absolute other-user homes: /home/<OTHER>/... and
        # /pfs/dataN/home/<OTHER>/... ; the member's own is protected
        r"/home/(?!" + esc + r"/)[\w.-]+/",
        r"/pfs/data\d+/home/(?!" + esc + r"/)[\w.-]+/",
        # ~<OTHER>/ home expansion; the member's own ~<owner>/ is allowed
        r"~(?!" + esc + r"/)[\w.-]+/",
    ]


def deny_patterns_for(role: str, *, home: Path | str | None = None
                      ) -> list[str]:
    """The role's effective deny pattern list, deduped, in stable order.

    Raises :class:`PolicyError` for an unknown role -- an unknown role is
    a typo today and a permission widening tomorrow (same rule the spec's
    closed vocabulary enforces).
    """
    if role not in ROLE_DENY:
        raise PolicyError(
            f"unknown role {role!r} (valid: {sorted(ROLE_DENY)})")
    merged: list[str] = []
    for group in (COMMON_DENY_PATTERNS,
                  tuple(home_deny_patterns(home)),
                  NO_PUSH_PATTERNS,
                  ROLE_DENY[role]):
        for pattern in group:
            if pattern not in merged:
                merged.append(pattern)
    return merged


# The owning-user key for the current process's home, for callers that
# need to name it without repeating the derivation.
def home_owner(home: Path | str | None = None) -> str:
    return _home_owner(_resolve_home(home))


def _resolve_home(home: Path | str | None = None) -> Path:
    return (Path(home) if home is not None else Path.home()).expanduser()
