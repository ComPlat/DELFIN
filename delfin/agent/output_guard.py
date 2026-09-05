"""Output guardrail stage for the FINAL agent answer.

Every serious agent framework runs a guard stage over the finished answer
before it is stored/forwarded; until now DELFIN's verify_guard-style
scanners were dashboard-only and nothing inspected the final text in the
engine itself.  This module is that stage: a small ordered pipeline of
checks over the complete answer text.

Checks (each individually fault-isolated — a crashing check never breaks
the answer):

1. SECRET REDACTION (mutating).  High-precision patterns for credential
   material — AWS access keys, PEM private-key blocks, GitHub / Slack
   tokens, ``Bearer`` header tokens, and generic ``api_key/token/password``
   assignments whose value looks high-entropy.  The sensitive substring is
   replaced by ``[redacted:<kind>]`` and a finding is recorded.  Precision
   over recall: ordinary hex content hashes (git SHAs, sha256 sums), DOIs
   and chemical formula strings (SMILES) must NEVER be touched — they carry
   no assignment/prefix context, and pure-hex values of hash length are
   explicitly excluded even inside assignments.

2. ABSOLUTE-CERTAINTY SCAN (non-mutating, telemetry only).  Sentences
   claiming absolute certainty ("guaranteed", "definitely", "cannot fail",
   "100% sicher") are recorded as findings, capped at 3.  The integrity
   contract governs behaviour; this is observability.

Configuration lives under settings ``agent.output_guard``:
``{"enabled": true, "redact_secrets": true}`` (both default True).
Reading settings never raises.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field


@dataclass
class GuardResult:
    text: str
    findings: list = field(default_factory=list)
    changed: bool = False


# ---------------------------------------------------------------------------
# Secret patterns (precision over recall)
# ---------------------------------------------------------------------------

# Full PEM block (BEGIN..END) — non-greedy so it stops at the first END.
_PRIVATE_KEY_BLOCK = re.compile(
    r"-----BEGIN [A-Z ]*PRIVATE KEY-----"
    r".*?"
    r"-----END [A-Z ]*PRIVATE KEY-----",
    re.DOTALL,
)

# Dangling header (stream cut off before END): header plus following
# base64 runs of PEM line length (>= 40 chars) only — ordinary prose words
# after a stray header are never swallowed.
_PRIVATE_KEY_DANGLING = re.compile(
    r"-----BEGIN [A-Z ]*PRIVATE KEY-----(?:\s*[A-Za-z0-9+/=]{40,})*"
)

_AWS_ACCESS_KEY = re.compile(r"\bAKIA[0-9A-Z]{16}\b")

_GITHUB_TOKEN = re.compile(
    r"\b(?:gh[pousr]_[A-Za-z0-9]{20,}|github_pat_[A-Za-z0-9_]{20,})\b"
)

_SLACK_TOKEN = re.compile(r"\bxox[bpars]-[A-Za-z0-9\-]{10,}\b")

# Provider API keys in their own vendor-prefixed shape. The assignment
# rule below only fires when a KEY NAME precedes the value, so a key that
# simply appears in a traceback, a curl echo or an error body -- which is
# how one actually leaks -- matched nothing at all. The vendor prefixes
# are specific enough to need no entropy check: nothing else looks like
# them. Anthropic is first because it is the one this project ships with.
_PROVIDER_API_KEY = re.compile(
    r"\b(?:"
    r"sk-ant-[A-Za-z0-9\-_]{16,}"        # Anthropic
    r"|sk-proj-[A-Za-z0-9\-_]{16,}"      # OpenAI project
    r"|sk-[A-Za-z0-9]{32,}"              # OpenAI classic / compatible
    r"|AIza[0-9A-Za-z\-_]{30,}"          # Google
    r"|hf_[A-Za-z0-9]{20,}"              # Hugging Face
    r")\b"
)

# ``Bearer <token>`` in headers: only the token part is redacted, and only
# when the candidate passes the high-entropy validation below (so prose
# like "Bearer authentication requires ..." is never touched).
_BEARER_TOKEN = re.compile(r"(?i)\bbearer\s+([A-Za-z0-9._\-+/=]{16,})")

# Generic ``api_key/token/password = <value>`` assignments.  The key name
# anchors the match — bare high-entropy strings in prose never match, which
# is what keeps git SHAs, DOIs and SMILES safe.
_ASSIGNMENT = re.compile(
    # A lookbehind, not \b: in OPENAI_API_KEY the character before "API"
    # is an underscore, which IS a word character, so the boundary never
    # matched and every PROVIDER_API_KEY=... assignment -- the commonest
    # shape there is -- went unredacted. "not preceded by a letter or
    # digit" keeps the same protection against matching inside a longer
    # word while allowing the underscore-separated prefix.
    r"(?i)(?<![A-Za-z0-9])("
    r"api[_-]?key|apikey|access[_-]?key|secret[_-]?key|client[_-]?secret|"
    r"auth[_-]?token|token|password|passwd|secret"
    r")\b[\"']?\s*[:=]\s*[\"']?"
    r"([A-Za-z0-9+/=_.\-]{16,})"
)

# Pure-hex strings of standard digest lengths are content hashes
# (md5 / git sha1 / sha256), not credentials.
_HEX_DIGEST_LENGTHS = (32, 40, 64)


def _plausible_secret(value: str) -> bool:
    """True when ``value`` looks like real credential material: long enough,
    mixed character classes, not a content hash, not a placeholder."""
    if len(value) < 16:
        return False
    if value.startswith("[redacted:"):
        return False
    classes = sum((
        any(c.islower() for c in value),
        any(c.isupper() for c in value),
        any(c.isdigit() for c in value),
    ))
    if classes < 2:
        return False
    if (len(value) in _HEX_DIGEST_LENGTHS
            and re.fullmatch(r"[0-9a-fA-F]+", value)):
        return False
    if len(set(value)) < 6:            # aaaa..., 1212... placeholders
        return False
    return True


def _redact(text: str, pattern: re.Pattern, kind: str, findings: list,
            *, group: int = 0, validate: bool = False) -> str:
    """Replace each match (or match ``group``) with ``[redacted:<kind>]``,
    recording one finding per redaction."""

    def repl(m: re.Match) -> str:
        whole = m.group(0)
        if group == 0:
            findings.append({"check": "secret_redaction", "detail": kind})
            return f"[redacted:{kind}]"
        candidate = m.group(group)
        if validate and not _plausible_secret(candidate):
            return whole
        findings.append({"check": "secret_redaction", "detail": kind})
        s, e = m.span(group)
        base = m.start(0)
        return (whole[: s - base]
                + f"[redacted:{kind}]"
                + whole[e - base:])

    return pattern.sub(repl, text)


# Ordered pipeline: specific prefixes before the generic assignment check so
# each redaction carries the most precise kind (and is not re-matched).
_SECRET_CHECKS: list[tuple[re.Pattern, str, int, bool]] = [
    (_PRIVATE_KEY_BLOCK, "private-key", 0, False),
    (_PRIVATE_KEY_DANGLING, "private-key", 0, False),
    (_AWS_ACCESS_KEY, "aws-access-key", 0, False),
    (_GITHUB_TOKEN, "github-token", 0, False),
    (_SLACK_TOKEN, "slack-token", 0, False),
    (_PROVIDER_API_KEY, "provider-api-key", 0, False),
    (_BEARER_TOKEN, "bearer-token", 1, True),
    (_ASSIGNMENT, "credential-assignment", 2, True),
]


def _redact_secrets(text: str, findings: list) -> str:
    for pattern, kind, group, validate in _SECRET_CHECKS:
        try:
            text = _redact(text, pattern, kind, findings,
                           group=group, validate=validate)
        except Exception:
            continue
    return text


# ---------------------------------------------------------------------------
# Absolute-certainty scan (telemetry only)
# ---------------------------------------------------------------------------

_ABSOLUTE = re.compile(
    r"(?i)\b(?:"
    r"guaranteed|definitely|cannot\s+fail|can'?t\s+fail|never\s+fails|"
    r"always\s+works|garantiert|"
    r"100\s*%\s*(?:sicher|sure|certain)"
    r")\b"
)

_MAX_ABSOLUTE_FINDINGS = 3


def _scan_absolutes(text: str, findings: list) -> None:
    count = 0
    for sentence in re.split(r"(?<=[.!?])\s+|\n+", text):
        if count >= _MAX_ABSOLUTE_FINDINGS:
            break
        if sentence and _ABSOLUTE.search(sentence):
            findings.append({
                "check": "absolute_certainty",
                "detail": sentence.strip()[:120],
            })
            count += 1


# ---------------------------------------------------------------------------
# Config + pipeline entry point
# ---------------------------------------------------------------------------

def _load_config() -> dict:
    """Settings ``agent.output_guard`` as a dict; never raises."""
    try:
        from delfin.user_settings import load_settings
        agent_cfg = (load_settings() or {}).get("agent", {}) or {}
        cfg = agent_cfg.get("output_guard", {}) or {}
        if isinstance(cfg, dict):
            return cfg
    except Exception:
        pass
    return {}


def run_output_guards(text: str, *, config: dict | None = None) -> GuardResult:
    """Run the ordered guard pipeline over the final answer ``text``.

    Returns the (possibly redacted) text plus all findings.  ``changed`` is
    True only when the text was mutated (i.e. redactions happened) — the
    absolute-certainty scan never mutates.
    """
    if not text:
        return GuardResult(text=text or "")

    if config is None:
        config = _load_config()
    if not isinstance(config, dict):
        config = {}
    if not config.get("enabled", True):
        return GuardResult(text=text)

    findings: list[dict] = []
    guarded = text

    if config.get("redact_secrets", True):
        try:
            guarded = _redact_secrets(guarded, findings)
        except Exception:
            guarded = text  # never return a half-broken mutation

    try:
        _scan_absolutes(guarded, findings)
    except Exception:
        pass

    return GuardResult(text=guarded, findings=findings,
                       changed=guarded != text)
