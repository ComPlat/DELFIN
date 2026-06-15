"""Persistent memory store for the DELFIN Agent.

Stores facts, preferences, and project context that persist across
agent cycles and dashboard sessions.

Files live exclusively in ``~/.delfin/`` (per-user, per-machine) and
are written with 0600 permissions — they MUST NOT be committed to
the repo.
"""

from __future__ import annotations

import json
import os
import re
import time
from pathlib import Path
from typing import Any


_DEFAULT_PATH = Path.home() / ".delfin" / "agent_memory.json"
_STOPWORDS = {
    "the", "and", "for", "with", "this", "that", "from", "into", "your",
    "about", "task", "user", "agent", "project", "repo", "code", "file",
}


def _tokenize(text: str) -> set[str]:
    return {
        token
        for token in re.split(r"[^a-z0-9]+", (text or "").lower())
        if len(token) >= 3 and token not in _STOPWORDS
    }


def _read(path: Path) -> dict[str, Any]:
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except (FileNotFoundError, json.JSONDecodeError):
        return {"facts": [], "updated_at": 0}


def _write(path: Path, data: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    data["updated_at"] = time.time()
    path.write_text(json.dumps(data, indent=2, ensure_ascii=False), encoding="utf-8")
    try:
        os.chmod(path, 0o600)
    except OSError:
        pass


def load_memories(path: Path | None = None) -> list[dict[str, Any]]:
    """Return list of stored memory facts."""
    return _read(path or _DEFAULT_PATH).get("facts", [])


def save_memory(
    text: str,
    source: str = "user",
    path: Path | None = None,
) -> int:
    """Append a memory fact. Returns the index of the new fact."""
    p = path or _DEFAULT_PATH
    data = _read(p)
    facts = data.get("facts", [])
    facts.append({
        "text": text.strip(),
        "source": source,
        "created_at": time.time(),
    })
    data["facts"] = facts
    _write(p, data)
    return len(facts) - 1


def delete_memory(index: int, path: Path | None = None) -> bool:
    """Remove a memory by index. Returns True if removed."""
    p = path or _DEFAULT_PATH
    data = _read(p)
    facts = data.get("facts", [])
    if 0 <= index < len(facts):
        facts.pop(index)
        data["facts"] = facts
        _write(p, data)
        return True
    return False


def format_memory_context(
    path: Path | None = None,
    *,
    task_text: str = "",
    max_entries: int = 4,
) -> str:
    """Format the most-relevant memories into a compact prompt block.

    Ranking is BM25-flavoured: per-fact score is the sum over the task's
    query tokens of  ``log((N - df + 0.5) / (df + 0.5))``  multiplied by
    a length-normalised term-frequency factor (k1=1.2, b=0.75 in the
    standard BM25 saturation curve). Bag-of-words is computed with the
    existing tokenizer so feedback-typed facts still slot in cleanly.

    Falls back to recency (last_n) when the task has no tokens after
    stopword removal, matching the previous behaviour.
    """
    facts = load_memories(path)
    if not facts:
        return ""

    task_tokens = _tokenize(task_text)
    if not task_tokens:
        selected: list[tuple[float, dict[str, Any]]] = [
            (float(i), fact)
            for i, fact in enumerate(facts[-max_entries:])
        ]
    else:
        import math as _math
        N = len(facts)
        # Pre-tokenize each fact + compute doc length / avg-doc length.
        docs = [_tokenize(str(f.get("text", ""))) for f in facts]
        doc_lens = [max(1, len(d)) for d in docs]
        avgdl = sum(doc_lens) / max(1, len(doc_lens))
        # Document frequency for each query token.
        df: dict[str, int] = {}
        for tok in task_tokens:
            df[tok] = sum(1 for d in docs if tok in d)
        k1, b = 1.2, 0.75
        scored: list[tuple[float, dict[str, Any]]] = []
        for i, (fact, doc, dlen) in enumerate(zip(facts, docs, doc_lens)):
            score = 0.0
            for tok in task_tokens:
                df_t = df.get(tok, 0)
                if df_t == 0 or tok not in doc:
                    continue
                idf = _math.log((N - df_t + 0.5) / (df_t + 0.5) + 1.0)
                # Token frequency is binary here (we don't track tf in a
                # set), so the saturation curve simplifies to:
                norm = (k1 + 1.0) / (1.0 + k1 * (1 - b + b * dlen / avgdl))
                score += idf * norm
            if score > 0:
                # Tiebreak on recency: later facts win on equal score.
                scored.append((score + 1e-6 * i, fact))
        scored.sort(key=lambda kv: kv[0], reverse=True)
        selected = scored[:max_entries]

    if not selected:
        return ""

    lines = []
    for _, fact in selected:
        src = fact.get("source", "?")
        lines.append(f"({src}) {fact['text']}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Typed project-memory writer (.delfin memory layout)
# ---------------------------------------------------------------------------
#
# The typed per-project memory lives under DELFIN's OWN namespace
# ``~/.delfin/projects/<slug>/memory/`` (MEMORY.md index + typed sidecar
# files) — NOT ``~/.claude`` (that is Claude Code's directory). An existing
# store at the legacy ~/.claude location is migrated on first access. This is
# the single source of truth the prompt loader recalls each turn; the flat
# ``~/.delfin/agent_memory.json`` store is kept only for backwards-compat.

_TYPE_LABELS = {
    "user": "User (background)",
    "feedback": "Feedback (how to work)",
    "project": "Project (current work)",
    "reference": "Reference (read only when relevant)",
}

_TYPE_PREFIX_RE = re.compile(
    r"^(user|feedback|project|reference)\s*[:\-]\s*", re.IGNORECASE
)


# Heuristic keyword sets for auto-classification when no explicit prefix
# is supplied. The first matching type wins. Order matters: feedback is
# checked first (most user-corrective utterances start with a prohibition),
# then project (time/deadline cues), then reference (external pointers),
# then user as the fallback. Patterns are intentionally conservative so a
# misclassification stays an editable text file the user can re-tag.
_FEEDBACK_HINTS = (
    "don't ", "do not ", "never ", "stop ", "always ",
    "nicht ", "keine ", "kein ", "immer ", "stets ",
    "use ... instead", "use x instead of",
    "prefer ", "should not ", "must not ", "muss nicht ",
    "muss immer", "darf nicht", "anti-pattern", "anti pattern",
)
_PROJECT_HINTS = (
    "deadline", "due ", "due:", "by friday", "by monday", "by tuesday",
    "by wednesday", "by thursday", "next week", "this week", "sprint",
    "milestone", "freeze", "release branch", "release on",
    "frist", "abgabe", "release", "cutoff",
    "is working on", "is doing", "ist dran", "arbeitet an",
)
_REFERENCE_HINTS = (
    "https://", "http://", "linear ", "jira ", "asana ",
    "grafana", "kibana", "datadog", "slack channel", "slack #",
    "github.com/", "gitlab.com/", "see issue ", "see pr ",
    "look at ", "siehe ", "documented at ", "dashboard at ",
)


def _classify_by_heuristic(text: str) -> str:
    """Pick a memory type when the user didn't supply a prefix.

    Returns one of ``"feedback" | "project" | "reference" | "user"``.
    """
    lowered = (text or "").lower()
    if any(h in lowered for h in _FEEDBACK_HINTS):
        return "feedback"
    if any(h in lowered for h in _PROJECT_HINTS):
        return "project"
    if any(h in lowered for h in _REFERENCE_HINTS):
        return "reference"
    return "user"


def parse_memory_type(text: str) -> tuple[str, str]:
    """Strip a leading ``<type>:`` prefix and return ``(memory_type, body)``.

    When no prefix is present we run a keyword heuristic across feedback /
    project / reference patterns and fall back to ``"user"``. The body is
    always the original text minus any consumed prefix.
    """
    m = _TYPE_PREFIX_RE.match(text or "")
    if m:
        return m.group(1).lower(), (text[m.end():]).strip()
    body = (text or "").strip()
    return _classify_by_heuristic(body), body


def _slugify(text: str, max_len: int = 60) -> str:
    s = re.sub(r"[^a-z0-9\s\-]", "", (text or "").lower())
    s = re.sub(r"\s+", "-", s).strip("-")
    return (s or "memory")[:max_len].rstrip("-") or "memory"


def _project_slug(repo_root: Path) -> str:
    return "-" + str(Path(repo_root).resolve()).replace("/", "-").lstrip("-")


def _migrate_legacy_dir(old: Path, new: Path) -> None:
    """One-time move of a per-project store out of the legacy ~/.claude path
    into DELFIN's own ~/.delfin namespace (same filesystem → cheap rename).
    Best-effort; never raises."""
    try:
        if old.is_dir() and not new.exists():
            new.parent.mkdir(parents=True, exist_ok=True)
            old.rename(new)
    except Exception:
        pass


def _delfin_memory_dir(repo_root: Path) -> Path:
    """DELFIN's OWN per-project memory store: ``~/.delfin/projects/<slug>/
    memory`` (not ``~/.claude`` — that is Claude Code's namespace). Migrates
    an existing store from the legacy ~/.claude location on first access."""
    slug = _project_slug(repo_root)
    new = Path.home() / ".delfin" / "projects" / slug / "memory"
    _migrate_legacy_dir(
        Path.home() / ".claude" / "projects" / slug / "memory", new)
    return new


def _delfin_plans_dir(repo_root: Path) -> Path:
    """Sibling of the memory dir — approved Plan-Mode plans land here so the
    user can re-open them in a later session. Under DELFIN's own ~/.delfin
    (migrated from the legacy ~/.claude location)."""
    slug = _project_slug(repo_root)
    new = Path.home() / ".delfin" / "projects" / slug / "plans"
    _migrate_legacy_dir(
        Path.home() / ".claude" / "projects" / slug / "plans", new)
    return new


# Back-compat aliases — the per-project store moved from ~/.claude (Claude
# Code's namespace) into DELFIN's own ~/.delfin. Old callers/tests keep working.
_claude_memory_dir = _delfin_memory_dir
_claude_plans_dir = _delfin_plans_dir


def list_plans(repo_root: Path | str) -> list[dict]:
    """Return one record per saved plan-file under the project's plans
    dir. Newest first (by created_at frontmatter; falls back to mtime)."""
    plans_dir = _delfin_plans_dir(Path(repo_root))
    if not plans_dir.is_dir():
        return []
    out: list[dict] = []
    for p in plans_dir.glob("*.md"):
        try:
            text = p.read_text(encoding="utf-8")
        except OSError:
            continue
        meta: dict[str, str] = {}
        body = text
        if text.startswith("---\n"):
            try:
                _, fm, body = text.split("---\n", 2)
                for line in fm.splitlines():
                    if ":" in line:
                        k, _, v = line.partition(":")
                        meta[k.strip()] = v.strip()
            except ValueError:
                body = text
        try:
            created_at = int(meta.get("created_at") or int(p.stat().st_mtime))
        except (TypeError, ValueError):
            created_at = int(p.stat().st_mtime)
        out.append({
            "file": p.name,
            "path": str(p),
            "name": meta.get("name") or p.stem,
            "description": meta.get("description") or "",
            "created_at": created_at,
            "body": body.strip(),
        })
    out.sort(key=lambda r: r["created_at"], reverse=True)
    return out


def get_plan(repo_root: Path | str, name_or_file: str) -> dict | None:
    """Return one plan by name or filename. None if not found."""
    target = (name_or_file or "").strip().lower()
    if not target:
        return None
    for rec in list_plans(repo_root):
        if (rec["file"].lower() == target
                or rec["file"].lower() == target + ".md"
                or rec["name"].lower() == target
                or rec["name"].lower() == target.removesuffix(".md")):
            return rec
    return None


def delete_plan(repo_root: Path | str, name_or_file: str) -> Path | None:
    """Delete a saved plan by name/filename. Returns the path that was
    removed, or None if not found."""
    rec = get_plan(repo_root, name_or_file)
    if rec is None:
        return None
    p = Path(rec["path"])
    try:
        p.unlink()
    except OSError:
        return None
    return p


def save_plan(
    plan_body: str,
    *,
    repo_root: Path | str,
    title: str | None = None,
) -> Path:
    """Persist an approved plan-mode plan to ``<plans-dir>/<slug>.md``.

    Frontmatter records the slug, a one-line description (derived from
    the first heading or first line of the body) and the creation
    timestamp. Returns the file path that was written.
    """
    body = (plan_body or "").strip()
    if not body:
        raise ValueError("plan body is empty")

    # Derive a title for the slug + description.
    first_meaningful = ""
    for line in body.splitlines():
        stripped = line.strip().lstrip("#").strip()
        if stripped:
            first_meaningful = stripped
            break
    display_title = (title or first_meaningful or "plan")[:80].strip() or "plan"
    slug = _slugify(display_title)

    plans_dir = _delfin_plans_dir(Path(repo_root))
    plans_dir.mkdir(parents=True, exist_ok=True)

    fname = f"{slug}.md"
    fpath = plans_dir / fname
    if fpath.exists():
        # Don't clobber — keep prior plans for diff / history.
        fname = f"{slug}_{int(time.time())}.md"
        fpath = plans_dir / fname

    description = first_meaningful[:160] or display_title
    front = (
        "---\n"
        f"name: {slug}\n"
        f"description: {description}\n"
        f"created_at: {int(time.time())}\n"
        "---\n\n"
        f"{body}\n"
    )
    fpath.write_text(front, encoding="utf-8")
    return fpath


def save_typed_memory(
    text: str,
    *,
    repo_root: Path | str,
    memory_type: str | None = None,
    title: str | None = None,
) -> tuple[Path, str, str]:
    """Persist a typed memory in the .delfin project-memory layout.

    Writes:
    - ``~/.claude/projects/<slug>/memory/<type>_<kebab-slug>.md`` with
      ``name:``, ``description:`` and ``metadata.type`` frontmatter
    - Prepends a one-line pointer to that file under the matching section
      of ``MEMORY.md``, creating the section if missing.

    Returns ``(file_path, slug, memory_type)``.
    """
    body = (text or "").strip()
    if memory_type is None:
        memory_type, body = parse_memory_type(body)
    memory_type = (memory_type or "user").lower()
    if memory_type not in _TYPE_LABELS:
        memory_type = "user"

    first_line = (body.splitlines() or [""])[0].strip()
    display_title = (title or first_line)[:80].strip() or "memory"
    slug = _slugify(display_title)

    memory_dir = _delfin_memory_dir(Path(repo_root))
    memory_dir.mkdir(parents=True, exist_ok=True)

    fname = f"{memory_type}_{slug}.md"
    fpath = memory_dir / fname
    if fpath.exists():
        # Don't clobber — disambiguate with a monotonic suffix.
        fname = f"{memory_type}_{slug}_{int(time.time())}.md"
        fpath = memory_dir / fname

    description = first_line[:160] or display_title
    front = (
        "---\n"
        f"name: {slug}\n"
        f"description: {description}\n"
        f"metadata:\n"
        f"  type: {memory_type}\n"
        "---\n\n"
        f"{body}\n"
    )
    fpath.write_text(front, encoding="utf-8")

    _update_memory_index(
        memory_dir,
        memory_type=memory_type,
        title=display_title,
        filename=fname,
        hook=description,
    )
    return fpath, slug, memory_type


def _update_memory_index(
    memory_dir: Path,
    *,
    memory_type: str,
    title: str,
    filename: str,
    hook: str,
) -> None:
    """Insert a one-line link under the matching ``## <Label>`` section of
    ``MEMORY.md``. Creates the index file or the section header if missing."""
    index = memory_dir / "MEMORY.md"
    label = _TYPE_LABELS[memory_type]
    section_header = f"## {label}"
    line = f"- [{title}]({filename}) — {hook}"

    if index.exists():
        content = index.read_text(encoding="utf-8")
    else:
        content = "# DELFIN Project Memory\n\n"

    if section_header in content:
        idx = content.index(section_header)
        try:
            eol = content.index("\n", idx) + 1
        except ValueError:
            eol = len(content)
        content = content[:eol] + line + "\n" + content[eol:]
    else:
        content = content.rstrip() + f"\n\n{section_header}\n{line}\n"

    index.write_text(content, encoding="utf-8")


def _type_from_filename(fname: str) -> str:
    for t in ("feedback", "project", "reference", "user"):
        if fname.startswith(t + "_"):
            return t
    return ""


def _parse_frontmatter(text: str) -> tuple[dict[str, str], str]:
    """Return ({name, description, type}, body) from a memory file's text."""
    meta: dict[str, str] = {}
    body = text
    if text.startswith("---\n"):
        try:
            _, fm, body = text.split("---\n", 2)
        except ValueError:
            return meta, text
        for line in fm.splitlines():
            s = line.strip()
            if s.startswith("type:"):
                meta["type"] = s.split(":", 1)[1].strip()
            elif ":" in s and not s.startswith(("metadata", "-")):
                k, _, v = s.partition(":")
                meta.setdefault(k.strip(), v.strip())
    return meta, body


def list_typed_memories(repo_root: Path | str) -> list[dict]:
    """One record per typed memory file under the project's memory dir.

    Sorted by type then name. Each record: file, path, name, description,
    type, body. This is the single source of truth for the agent's learned
    memories (the legacy flat JSON store is no longer used for recall)."""
    memory_dir = _delfin_memory_dir(Path(repo_root))
    if not memory_dir.is_dir():
        return []
    out: list[dict] = []
    for p in sorted(memory_dir.glob("*.md")):
        if p.name == "MEMORY.md":
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except OSError:
            continue
        meta, body = _parse_frontmatter(text)
        out.append({
            "file": p.name,
            "path": str(p),
            "name": meta.get("name") or p.stem,
            "description": meta.get("description") or "",
            "type": meta.get("type") or _type_from_filename(p.name) or "user",
            "body": body.strip(),
        })
    out.sort(key=lambda r: (r["type"], r["name"]))
    return out


def _remove_memory_index_line(memory_dir: Path, filename: str) -> None:
    """Drop the MEMORY.md pointer line that references ``filename``."""
    index = memory_dir / "MEMORY.md"
    if not index.exists():
        return
    try:
        lines = index.read_text(encoding="utf-8").splitlines(keepends=True)
    except OSError:
        return
    kept = [ln for ln in lines if f"({filename})" not in ln]
    if len(kept) != len(lines):
        try:
            index.write_text("".join(kept), encoding="utf-8")
        except OSError:
            pass


def delete_typed_memory(repo_root: Path | str, name_or_file: str) -> Path | None:
    """Delete a typed memory by name / slug / filename and remove its
    MEMORY.md pointer. Returns the deleted path, or None if not found."""
    target = (name_or_file or "").strip().lower()
    if not target:
        return None
    memory_dir = _delfin_memory_dir(Path(repo_root))
    match = None
    for rec in list_typed_memories(repo_root):
        if (rec["file"].lower() == target
                or rec["file"].lower() == target + ".md"
                or rec["name"].lower() == target
                or rec["name"].lower() == target.removesuffix(".md")):
            match = rec
            break
    if match is None:
        return None
    p = Path(match["path"])
    try:
        p.unlink()
    except OSError:
        return None
    _remove_memory_index_line(memory_dir, match["file"])
    return p


# ---------------------------------------------------------------------------
# Memory enrichment: cross-links, stale-ref verification, recall helpers
# ---------------------------------------------------------------------------

_WIKILINK_RE = re.compile(r"\[\[([a-zA-Z0-9][a-zA-Z0-9_\-]*)\]\]")
# Detect "path[:line]" style references inside memory text so the agent can
# verify they still exist before quoting them. Avoid grabbing trivial words
# by requiring at least one slash OR a dotted file-like suffix.
_PATH_REF_RE = re.compile(
    r"(?:(?<=\s)|(?<=^)|(?<=`))"           # left boundary
    r"((?:[\w./\-]+/[\w./\-]+|[\w./\-]+\.[a-zA-Z]{1,5}))"
    r"(?::(\d+))?"                          # optional :line
    r"(?=[\s)\.,;'\"`]|$)"                   # right boundary
)


def _find_memory_file_by_name(memory_dir: Path, name: str) -> Path | None:
    """Locate ``<any-type>_<name>.md`` under the memory dir; first hit wins.

    We don't know the type prefix at link-resolution time, so we glob the
    four legal prefixes in deterministic priority order. Missing dir → None.
    """
    if not memory_dir.is_dir():
        return None
    for kind in ("feedback", "project", "reference", "user"):
        cand = memory_dir / f"{kind}_{name}.md"
        if cand.is_file():
            return cand
    # Fallback: bare ``<name>.md`` so externally-curated MEMORY.md entries
    # without a type prefix still resolve.
    flat = memory_dir / f"{name}.md"
    return flat if flat.is_file() else None


def resolve_wikilinks(text: str, memory_dir: Path) -> str:
    """Expand each ``[[name]]`` in ``text`` to a resolved markdown link.

    - ``[[X]]`` with an existing memory file → ``[X](file.md) — description``
    - ``[[X]]`` without a target → ``[[X]] (not yet written)``

    Idempotent for resolved links (we only touch the ``[[…]]`` syntax).
    """
    if not text or "[[" not in text:
        return text

    def _replace(m: re.Match) -> str:
        name = m.group(1)
        target = _find_memory_file_by_name(memory_dir, name)
        if target is None:
            return f"[[{name}]] (not yet written)"
        # Pull the description line from frontmatter for the inline hook.
        try:
            head = target.read_text(encoding="utf-8")[:400]
        except OSError:
            return f"[[{name}]] (unreadable)"
        desc = ""
        for line in head.splitlines():
            if line.lower().startswith("description:"):
                desc = line.split(":", 1)[1].strip().strip("\"'")
                break
        if desc:
            return f"[{name}]({target.name}) — {desc}"
        return f"[{name}]({target.name})"

    return _WIKILINK_RE.sub(_replace, text)


def find_stale_references(text: str, repo_root: Path | str) -> list[str]:
    """Return a list of ``path[:line]`` references mentioned in the memory
    text that no longer exist on disk.

    Used by ``/memories verify`` to flag rotted recommendations. The check
    is conservative: only entries that *look* like paths (slash- or
    dot-shaped) are tested, so plain prose like "delete" or "method" is
    skipped.
    """
    root = Path(repo_root)
    stale: list[str] = []
    for match in _PATH_REF_RE.finditer(text or ""):
        ref = match.group(0)
        path_part = match.group(1)
        # Skip URL-like, version-y, or trivial fragments
        if "://" in path_part or path_part.startswith(("v", "V")) and path_part[1:2].isdigit():
            continue
        if "." not in path_part and "/" not in path_part:
            continue
        candidate = root / path_part if not Path(path_part).is_absolute() else Path(path_part)
        try:
            if not candidate.exists():
                stale.append(ref)
        except OSError:
            stale.append(ref)
    return stale


def verify_typed_memories(repo_root: Path | str) -> list[dict]:
    """Walk every memory file under the project's memory dir and return a
    list of ``{file, stale_refs}`` records for any file containing one or
    more dead references. Empty list = everything still resolves.
    """
    memory_dir = _delfin_memory_dir(Path(repo_root))
    if not memory_dir.is_dir():
        return []
    results: list[dict] = []
    for p in sorted(memory_dir.glob("*.md")):
        if p.name == "MEMORY.md":
            continue
        try:
            text = p.read_text(encoding="utf-8")
        except OSError:
            continue
        stale = find_stale_references(text, repo_root)
        if stale:
            results.append({"file": p.name, "stale_refs": stale})
    return results
