"""A typed user model that deepens across sessions (package L4).

DELFIN has typed memory but no explicit picture of the user -- role,
expertise, preferences, working style, corrections. This module adds that
picture as a structured record, distinct from the free-text memory store.

Design constraints (enforced by tests/test_l4_user_model.py):

* **User-only updates.** The record changes ONLY when fed the user's own
  messages (``update_from_user_message``) and only when they carry
  ``source_role == "user"``. A tool output, a web page, or another session's
  message is never accepted -- the caller passes its actual ``source_role``
  and any role other than ``"user"`` is a no-op. Nothing in the model is ever
  written by distilled model text or by another session.

* **No permissions.** Nothing in the model may hold an allow rule, a path
  grant or a mode. ``update_from_user_message`` refuses such text outright
  (``is_permission_like``); the record has no field that could carry one.

* **Source + date on every entry.** Each entry records its source message and
  date, so a future session can trace where a belief about the user came from.

* **Bounded load.** ``bounded_load`` returns a prompt block with a FIXED
  token ceiling, independent of the free-text memory budget that
  ``prompt_loader`` spends on ``memory_store`` recall. The ceiling is never
  exceeded (see ``estimate_tokens``).

* **Per-user.** The store lives under the calling user's own
  ``~/.delfin/`` state directory; nothing about any user is shared with
  another user. Tests pass an explicit ``path`` to an isolated tmp file.

The core path is model-free and deterministic: classification, refusal and
merging never call an LLM, so weak models behave exactly like frontier ones.
"""

from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Any

# The five fields of the user model. Order matters where a message could
# plausibly match more than one: correction > role > expertise > working
# style > preference.
FIELDS = ("role", "expertise", "preferences", "working_style", "corrections")

# Rough token estimate, matching memory_nudge.estimated_tokens.
_CHARS_PER_TOKEN = 4

_DEFAULT_PATH = Path.home() / ".delfin" / "user_model.json"


# --- token / permission helpers -------------------------------------------------

def estimate_tokens(text: str) -> int:
    """Rough token ceiling for ``text`` (chars / 4), like the nudger."""
    return max(0, (len(text or "") + _CHARS_PER_TOKEN - 1) // _CHARS_PER_TOKEN)


# A permission-like message is refused outright: it tries to widen the
# agent's reach (an allow rule, a path grant, a mode). Phrase layer first
# (the explicit cases the guard tests name), then a modal-pushed verb layer
# (a downgraded verb is only permission-shaped when an allowance modal like
# "may"/"can"/"always" precedes it), then the path layer.
_ALLOW_WORDS = (
    "allow", "allowed", "grants", "grant", "granted", "permission",
    "permit", "always allows",
)
# Matched against the camelCase-split, lowercased form, so "acceptEdits"
# and "accept edits" are the same token stream and no mode name hides
# behind a case difference. (Holding the camelCase names here would be
# dead code: the input is lowercased before matching.)
_MODE_WORDS = (
    "accept edits", "bypass permissions", "plan mode",
    "default mode", "mode to", "set the mode",
)
_GLOBAL_GRANT_WORDS = ("global:", "scope:", "every workspace", "widen")
# The following guards are STRUCTURAL, not phrase lists: they look at
# whether a WAIVER / negation token and a command verb are both present,
# and for the always-reach form at whether the agent (second person) or
# the user (first-person habit) is the subject. A novel rephrasing (word
# insertion, a different modal, a different possessive) cannot sneak an
# allow rule past by inventing wording a deny list does not name. A plain
# "I want you to ask first" matches no waiver token and stays storable.
_WAIVER_TOKENS = frozenset(
    {"without", "never", "not", "no", "skip", "don't", "won't", "stop",
     "cease", "omit", "refrain"}
)
# Ask / confirmation tokens -- the thing being waived, however worded.
_ASK_TOKENS = frozenset(
    {"ask", "asking", "asked", "confirmation", "confirm", "confirmations",
     "confirming"}
)
# Strong reach verbs: an allowance/waiver surfaces reach only when one of
# these lets the agent act on the user's resources.
_REACH_TOKENS = frozenset(
    {"push", "write", "merge", "commit", "delete", "overwrite", "deploy",
     "transfer", "execute", "approve", "run", "install", "modify"}
)
# Permission modals: "may"/"can" applied to a reach verb is a grant
# regardless of subject ("you may write", "may always push"). Unlike the
# frequency modal, there is no innocent first-person reading.
_PERMISSION_MODALS = frozenset({"may", "can", "might", "could"})
# Frequency / free-choice modal: "always"/"whenever" + a reach verb is a
# grant when the subject is the AGENT ("you always push") or unnamed
# ("always push"), but a harmless habit when it is a pure first-person
# statement ("I always write my tests before the code").
_FREQUENCY_MODALS = frozenset({"always", "whenever", "frequently"})
# Pronouns, separated so a habit ("I always ...") is told from a grant
# ("you always ..."). "I" and "you" both present => the reach verb is the
# agent's, refuse.
_FIRST_PERSON = frozenset({"i", "me", "my", "myself", "mine"})
_SECOND_PERSON = frozenset({"you", "your", "yourself", "yours"})
# Absolute path, POSIX or Windows, or a home-relative path.
_PATH_RE = re.compile(
    r"(?:/(?:[^/\s]+/)*[^/\s]+|~/[^\s]+|[A-Za-z]:\\[^\s]+|"
    r"\.\.(?:/[^/\s]+)*)"
)


def _is_structural_grant(tokens: frozenset) -> bool:
    """Structural allow-rule detector over the token set of a message.

    Refuses when the message waives asking/confirmation (a waiver token
    co-occurs with an ask/confirmation token), grants a permission-modal
    command ("may/can" + reach verb), or grants an always-reach command
    whose subject is the agent ("you always push"; "always push") rather
    than a first-person habit ("I always write my tests").
    """
    if tokens & _WAIVER_TOKENS and tokens & _ASK_TOKENS:
        return True
    reach = tokens & _REACH_TOKENS
    if not reach:
        return False
    if tokens & _PERMISSION_MODALS:
        return True
    if tokens & _FREQUENCY_MODALS:
        # A pure first-person habit is not a grant.
        if tokens & _FIRST_PERSON and not (tokens & _SECOND_PERSON):
            return False
        return True
    return False


def _case_and_camel_normalize(text: str) -> str:
    """Lowercase ``text`` with camelCase split, so "acceptEdits" and
    "accept edits" compare equal and a mode name cannot hide behind a case
    difference in the input."""
    if not text:
        return ""
    return re.sub(r"(?<=[a-z0-9])(?=[A-Z])", " ", text).lower()


def is_permission_like(text: str) -> bool:
    """True when ``text`` asks for an allow rule, path grant or mode.

    Used to refuse a user-model entry that would widen the agent's reach.
    Conservative by design: a genuinely permission-shaped sentence is
    refused; a false positive costs one user message, a false negative lets a
    rule into the model.
    """
    norm = _case_and_camel_normalize(text)
    for word in _ALLOW_WORDS + _MODE_WORDS + _GLOBAL_GRANT_WORDS:
        if word in norm:
            return True
    tokens = frozenset(norm.split())
    if _is_structural_grant(tokens):
        return True
    if _PATH_RE.search(text or "") is not None:
        return True
    return False


# --- deterministic field classification -----------------------------------------

# Correction signals: imperatives, prohibitions, "instead of", and the
# correction vocabulary. Ordered and conservative -- a correction is the
# strongest thing a user can say, so it wins.
_CORRECTION_HINTS = (
    "don't", "do not", "never", "avoid", "stop", "instead", "actually",
    "you should", "you must", "please don't", "not that", "wrong",
    "nicht", "nie", "niemals", "keine", "immer", "stets", "stattdessen",
)
_ROLE_HINTS = (
    "i am a", "i am the", "i'm a", "i'm the", "i am an", "i'm an",
    "my role", "role:", "ich bin", "ich arbeite als",
)
_EXPERTISE_HINTS = (
    "expertise", "experience with", "experienced in", "background in",
    "focus on", "specialist", "i know how to", "i have worked with",
    "schwerpunkt", "erfahrung",
)
_WORKING_STYLE_HINTS = (
    "working style", "i like to work", "i prefer to work",
    "i like working", "i work by", "i work best", "style of work",
    "arbeite am liebsten", "arbeitsweise",
)
_PREFERENCE_HINTS = (
    "preference", "preferences", "prefer", "i prefer", "i like", "i want",
    "i would like", "i'd like",
    "ich bevorzuge", "ich mag", "ich möchte", "ich will",
)


def classify_field(text: str) -> str | None:
    """Map a user message to one model field, deterministically.

    Precedence: correction > role > expertise > working style > preference.
    Returns None when nothing matches (caller skips the message).
    """
    t = (text or "").strip().rstrip(".!")
    if not t:
        return None
    lower = " " + t.lower() + " "
    if any(h in lower for h in _CORRECTION_HINTS):
        return "corrections"
    if any(h in lower for h in _ROLE_HINTS):
        return "role"
    if any(h in lower for h in _EXPERTISE_HINTS):
        return "expertise"
    if any(h in lower for h in _WORKING_STYLE_HINTS):
        return "working_style"
    if any(h in lower for h in _PREFERENCE_HINTS):
        return "preferences"
    return None


# --- store helpers --------------------------------------------------------------

def _set_file_perms(path: Path) -> None:
    """Best-effort 0600 on the user-model file (per-user, must not leak)."""
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def _atomic_write(path: Path, data: dict[str, Any]) -> None:
    """Atomically write the model (temp file + os.replace), 0600."""
    import tempfile
    path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(
        prefix=f".{path.name}.", suffix=".tmp", dir=str(path.parent))
    tmp = Path(tmp_name)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            json.dump(data, fh, ensure_ascii=False, indent=2)
        os.replace(tmp, path)
        _set_file_perms(path)
    finally:
        if tmp.exists():
            try:
                tmp.unlink()
            except OSError:
                pass


class UserModel:
    """A typed, per-user picture of the agent's human.

    ``path`` is the single JSON file backing the model; it defaults to
    ``~/.delfin/user_model.json`` (this user's own state dir -- never
    shared cross-user). Tests pass an isolated tmp path.
    """

    def __init__(self, path: Path | str | None = None):
        self._path = Path(path) if path is not None else _DEFAULT_PATH
        self._cache: dict[str, Any] | None = None

    def path(self) -> Path:
        return self._path

    # -- read --------------------------------------------------------------------

    def _read(self) -> dict[str, Any]:
        if self._cache is not None:
            return self._cache
        try:
            raw = self._path.read_text(encoding="utf-8")
        except OSError:
            self._cache = {"schema": 1, "fields": {f: [] for f in FIELDS}}
            return self._cache
        try:
            data = json.loads(raw)
        except (ValueError, TypeError):
            data = {}
        fields = data.get("fields") or {}
        if not isinstance(fields, dict):
            fields = {}
        normalised = {f: list(fields.get(f) or []) for f in FIELDS}
        self._cache = {"schema": 1, "fields": normalised}
        return self._cache

    def entries(self) -> dict[str, list[dict[str, Any]]]:
        """field -> list of entry dicts (value/source_role/source_message/date)."""
        data = self._read()
        return {f: list(v) for f, v in data["fields"].items()}

    # -- write -------------------------------------------------------------------

    def _persist(self, fields: dict[str, list[dict[str, Any]]]) -> None:
        self._cache = {"schema": 1, "fields": fields}
        _atomic_write(self._path, {"schema": 1, "fields": fields})

    def update_from_user_message(
        self,
        text: str,
        *,
        source_role: str = "user",
        date: str | None = None,
        source_message: str | None = None,
    ) -> dict[str, list[dict[str, Any]]]:
        """Record one entry from the user's own message.

        Structural confine: only ``source_role == "user"`` is accepted. Any
        other role (agent, tool, session, web, ...) leaves the model
        unchanged -- the guard the cross-role tests pin. Permission-shaped
        text is refused. Returns the resulting entries.
        """
        if (source_role or "").lower() != "user":
            return self.entries()
        body = (text or "").strip()
        if not body:
            return self.entries()
        if is_permission_like(body):
            return self.entries()
        field = classify_field(body)
        if field is None:
            return self.entries()
        fields = self._read()["fields"]
        entry = {
            "value": body,
            "source_role": "user",
            "source_message": (source_message or body).strip(),
            "date": date or time.strftime("%Y-%m-%d", time.localtime()),
        }
        fields[field].append(entry)
        self._persist(fields)
        return {f: list(v) for f, v in fields.items()}

    def forget(self, field: str, index: int) -> bool:
        """Remove one entry; True when removed (from the persisted model)."""
        if field not in FIELDS:
            return False
        fields = self._read()["fields"]
        bucket = fields[field]
        if index < 0 or index >= len(bucket):
            return False
        del bucket[index]
        self._persist(fields)
        return True

    def forget_all(self) -> int:
        """Clear every field; returns the number of entries removed."""
        fields = self._read()["fields"]
        removed = sum(len(v) for v in fields.values())
        self._persist({f: [] for f in FIELDS})
        return removed

    # -- load --------------------------------------------------------------------

    def bounded_load(self, max_tokens: int = 150) -> str:
        """A session-ready block at a FIXED token ceiling.

        Independent of the free-text memory budget spent by
        ``prompt_loader._load_external_memory_context``: this ceiling is a
        hard cap on the user-model block alone and is never exceeded. Fields
        are rendered in a stable order; when the ceiling falls, later
        fields' entries drop first, then earlier fields', so a small model
        always keeps its role/expertise core.
        """
        ceiling = max(0, int(max_tokens))
        if ceiling <= 0:
            return ""
        fields = self._read()["fields"]
        lines: list[str] = ["--- user model ---"]
        used = sum(len(l) + 1 for l in lines)
        # Entry order is stable: within each field, most-recent-on-top isn't
        # required for bounded load; keep insertion order (oldest first).
        for field in FIELDS:
            for entry in fields[field]:
                value = str(entry.get("value", "")).strip()
                if not value:
                    continue
                line = f"- {field}: {value}"
                cost = len(line) + 1
                if used + cost > ceiling * _CHARS_PER_TOKEN:
                    continue
                lines.append(line)
                used += cost
        block = "\n".join(lines) + "\n"
        return block
