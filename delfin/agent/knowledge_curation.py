"""Curation of the global knowledge (``delfin/agent/pack/knowledge``).

Read-only tools, for anyone and for CI:

* :func:`lint` -- every entry against the template and the taxonomy
  (``knowledge.template_findings``), the public-repository content rules
  (``knowledge_check.check(public=True)``, without this machine's identity
  so the result is the same on every machine), exact duplicates, and the
  generated index in sync. ``tests/test_knowledge_lint.py`` runs it over
  the real folder, so a PR with a malformed or duplicate entry fails CI.
* :func:`dupes` -- near-duplicate pairs and contradictions
  (``knowledge.near_duplicate``/``contradicts``).
* :func:`stats` -- entries per topic, expired ``review_by``, entries never
  loaded on this machine (``~/.delfin/knowledge/usage.json``, local only),
  and broken references: a backticked path that is not in the repository,
  a ``name()`` that ``git grep`` finds no definition for, a
  ``delfin-agent <command>`` the CLI does not define.

Write actions, for the knowledge admin only (``knowledge_share.is_admin``):
merge, move, retire, rename, split, index regeneration, topic add, and the
review queue's accept, edit, merge-into and reject. Each builds a
:class:`Plan` -- the change set and the diff it produces against the base
-- and :func:`start` opens the admin's own confirmation
(``knowledge_share.start_admin``: refused for anyone else and inside an
agent session; two codes for the public DELFIN repository). A confirmed
plan becomes a local branch through ``knowledge_share.prepare_branch``
with INDEX.md regenerated; nothing is pushed, the working tree is not
touched. A rejection writes only the queue record, with the reason the
author sees.

Rejected: letting the curation tools edit the folder in place (an
unreviewed rewrite of other people's entries); counting usage centrally
(load events would leave the machine).
"""
from __future__ import annotations

import datetime as dt
import re
import subprocess
from collections import Counter
from dataclasses import dataclass, field, replace
from pathlib import Path
from typing import Callable, Optional

from . import knowledge as kn
from . import knowledge_share as ks

__all__ = ["lint", "dupes", "stats", "Plan", "start", "plan_merge",
           "plan_move", "plan_retire", "plan_rename", "plan_split",
           "plan_index", "plan_topic_add", "plan_accept", "plan_reject",
           "plan_queue_merge", "plan_queue_edit"]


# ---------------------------------------------------------------------------
# read-only
# ---------------------------------------------------------------------------

def lint(folder: Optional[Path] = None) -> list[str]:
    from . import knowledge_check as kc
    folder = Path(folder or kn.global_dir())
    topics = kn.read_topics(folder)
    entries = kn.read_entries(folder)
    out: list[str] = []
    if not topics:
        out.append(f"{kn.TOPICS_FILE}: missing or empty taxonomy")
    for e in entries:
        out += [f"{e.slug}: {f}" for f in kn.template_findings(
            e, topics=topics, layer=kn.GLOBAL)]
        out += [f"{e.slug}: {f.line()}" for f in kc.check(
            e.text, public=True, identity={}, unpublished_terms=())]
    seen: dict = {}
    for e in entries:
        keys = {("text", kn._norm(f"{e.situation} {e.lesson}")),
                ("title", kn._norm(e.title))}
        for key in keys:
            if key[1] and key in seen:
                out.append(f"duplicate: {e.slug} repeats {seen[key]} "
                           f"(same {key[0]})")
            seen.setdefault(key, e.slug)
    if not kn.index_in_sync(folder):
        out.append(f"{kn.INDEX_FILE}: out of sync with the entries; run "
                   "`delfin-agent knowledge index`")
    return out


@dataclass(frozen=True)
class Pair:
    kind: str           # near-duplicate | contradiction
    a: str
    b: str


def dupes(folder: Optional[Path] = None) -> list[Pair]:
    entries = kn.read_entries(Path(folder or kn.global_dir()))
    out = []
    for i, a in enumerate(entries):
        for b in entries[i + 1:]:
            if kn.near_duplicate(a, b):
                out.append(Pair("near-duplicate", a.slug, b.slug))
            elif kn.contradicts(a, b):
                out.append(Pair("contradiction", a.slug, b.slug))
    return out


@dataclass(frozen=True)
class BrokenRef:
    slug: str
    ref: str
    kind: str


_TICK_RE = re.compile(r"`([^`\n]{2,120})`")


def _git_grep(repo: Path, pattern: str) -> Optional[bool]:
    try:
        out = subprocess.run(["git", "-C", str(repo), "grep", "-q", "-E",
                              pattern], capture_output=True, timeout=30)
    except (OSError, subprocess.SubprocessError):
        return None
    if out.returncode in (0, 1):
        return out.returncode == 0
    return None


def _tracked(repo: Path, ref: str) -> bool:
    """Whether some tracked file ends with ``ref`` (a partial path such as
    ``tools.py`` counts when one file is called that)."""
    try:
        out = subprocess.run(["git", "-C", str(repo), "ls-files", "--",
                              f"*{ref}"], capture_output=True, text=True,
                             timeout=30)
        return bool(out.stdout.strip()) if out.returncode == 0 else True
    except (OSError, subprocess.SubprocessError):
        return True       # cannot tell: not reported as broken


def broken_refs(e: kn.Entry, repo_root: Path) -> list[BrokenRef]:
    out = []
    for m in _TICK_RE.finditer(f"{e.situation} {e.lesson} {e.why}"):
        ref = m.group(1).strip()
        if re.fullmatch(r"[A-Za-z_][\w.]*\(\)", ref):
            name = ref[:-2].rsplit(".", 1)[-1]
            hit = _git_grep(repo_root, rf"(def|class) {re.escape(name)}\b"
                            rf"|^{re.escape(name)} *=")
            if hit is False:
                out.append(BrokenRef(e.slug, ref, "function"))
        elif ref.startswith("delfin-agent "):
            word = ref.split()[1]
            cli = repo_root / "delfin" / "agent" / "cli.py"
            try:
                if f'"{word}"' not in cli.read_text(encoding="utf-8"):
                    out.append(BrokenRef(e.slug, ref, "command"))
            except OSError:
                out.append(BrokenRef(e.slug, ref, "command"))
        elif " " not in ref and ("/" in ref or re.search(
                r"\.(py|md|json|yaml|yml|toml|sh|txt)$", ref)):
            if not (repo_root / ref).exists() and not _tracked(
                    repo_root, ref):
                out.append(BrokenRef(e.slug, ref, "path"))
    return out


def stats(folder: Optional[Path] = None, *, usage_path: Optional[Path]
          = None, repo_root: Optional[Path] = None,
          today: Optional[dt.date] = None) -> dict:
    folder = Path(folder or kn.global_dir())
    repo_root = Path(repo_root or ks.delfin_root())
    entries = kn.read_entries(folder)
    usage = kn.read_usage(usage_path)
    broken = []
    for e in entries:
        broken += broken_refs(e, repo_root)
    return {
        "entries": len(entries),
        "by_topic": dict(Counter(e.topic for e in entries)),
        "expired": sorted(e.slug for e in entries if e.expired(today)),
        "never_loaded": sorted(e.slug for e in entries
                               if not (usage.get(e.slug) or {}).get("loads")),
        "broken": broken,
    }


# ---------------------------------------------------------------------------
# plans
# ---------------------------------------------------------------------------

@dataclass
class Plan:
    title: str
    changes: dict
    message: str
    branch_hint: str
    public: bool = True
    after: Optional[Callable[[str], None]] = None
    diff: str = ""
    notes: list = field(default_factory=list)

    def preview(self) -> str:
        lines = [self.title] + self.notes
        lines.append("--- the change, as a diff against the base branch ---")
        lines.append(self.diff.strip() or "(no file change)")
        return "\n".join(lines)


def _repo() -> Path:
    return ks.delfin_root()


def _rel(name: str) -> str:
    return f"{ks.GLOBAL_FOLDER}/{name}"


def _entry(slug: str) -> kn.Entry:
    for e in kn.read_entries(kn.global_dir()):
        if e.slug == slug:
            return e
    raise ValueError(f"no global entry {slug!r}")


def _finish(plan: Plan) -> Plan:
    if plan.changes or plan.branch_hint == "index":
        plan.diff = ks.preview_diff(_repo(), plan.changes,
                                    index_folder=ks.GLOBAL_FOLDER)
    plan.public = plan.public and (bool(plan.changes)
                                   or plan.branch_hint == "index")
    return plan


def plan_merge(keep: str, drop: str) -> Plan:
    a, b = _entry(keep), _entry(drop)
    merged = kn.merge_entries(a, b)
    return _finish(Plan(
        f"Merge {drop} into {keep} (both provenances kept)",
        {_rel(f"{keep}.md"): kn.render_entry(merged),
         _rel(f"{drop}.md"): None},
        f"Knowledge: merge {drop} into {keep}\n\nProvenance kept: "
        + "; ".join(p.render() for p in merged.provenance),
        f"merge-{drop}"))


def plan_move(slug: str, topic: str) -> Plan:
    if topic not in kn.read_topics(kn.global_dir()):
        raise ValueError(f"{topic!r} is not in the taxonomy; add it first "
                         "(knowledge topic add)")
    e = _entry(slug)
    return _finish(Plan(
        f"Move {slug} from {e.topic} to {topic}",
        {_rel(f"{slug}.md"): kn.render_entry(replace(e, topic=topic))},
        f"Knowledge: move {slug} to topic {topic}", f"move-{slug}"))


def plan_retire(slug: str, *, reason: str) -> Plan:
    if not (reason or "").strip():
        raise ValueError("a retirement needs a reason")
    _entry(slug)
    return _finish(Plan(
        f"Retire {slug}: {reason}", {_rel(f"{slug}.md"): None},
        f"Knowledge: retire {slug}\n\nReason: {reason.strip()}",
        f"retire-{slug}"))


def plan_rename(slug: str, new_slug: str) -> Plan:
    if not kn._SLUG_RE.match(new_slug or ""):
        raise ValueError("the new name is lower-case letters, digits, '-'")
    if (kn.global_dir() / f"{new_slug}.md").exists():
        raise ValueError(f"{new_slug!r} exists already")
    e = _entry(slug)
    return _finish(Plan(
        f"Rename {slug} to {new_slug}",
        {_rel(f"{slug}.md"): None,
         _rel(f"{new_slug}.md"): kn.render_entry(replace(e, slug=new_slug))},
        f"Knowledge: rename {slug} to {new_slug}", f"rename-{slug}"))


def plan_split(slug: str, parts: list) -> Plan:
    e = _entry(slug)
    topics = kn.read_topics(kn.global_dir())
    changes: dict = {_rel(f"{slug}.md"): None}
    for part in parts:
        new = replace(e, **{k: v for k, v in part.items()
                            if k in ("slug", "title", "situation", "lesson",
                                     "why", "topic")})
        problems = kn.template_findings(new, topics=topics, layer=kn.GLOBAL)
        if problems:
            raise ValueError(f"{new.slug}: " + "; ".join(problems))
        changes[_rel(f"{new.slug}.md")] = kn.render_entry(new)
    return _finish(Plan(
        f"Split {slug} into " + ", ".join(p["slug"] for p in parts),
        changes, f"Knowledge: split {slug}", f"split-{slug}"))


def plan_index() -> Plan:
    return _finish(Plan("Regenerate INDEX.md from the entries", {},
                        "Knowledge: regenerate the index", "index"))


def plan_topic_add(topic: str, description: str = "") -> Plan:
    if not kn._SLUG_RE.match(topic or ""):
        raise ValueError("a topic is a lower-case slug")
    topics = kn.read_topics(kn.global_dir())
    if topic in topics:
        raise ValueError(f"{topic!r} is in the taxonomy already")
    topics[topic] = (description or "").strip()
    import tempfile
    with tempfile.TemporaryDirectory() as tmp:
        text = kn.write_topics(Path(tmp), topics).read_text(encoding="utf-8")
    return _finish(Plan(
        f"Add the topic {topic} to the taxonomy",
        {_rel(kn.TOPICS_FILE): text},
        f"Knowledge: add topic {topic}", f"topic-{topic}"))


# -- the review queue --------------------------------------------------------

def _record(qid: str) -> dict:
    for rec in ks.read_queue():
        if rec.get("id") == qid:
            if rec.get("status") != "pending":
                raise ValueError(f"{qid} is {rec.get('status')} already")
            return rec
    raise ValueError(f"no queued submission {qid!r}")


def _queued_entry(rec: dict) -> kn.Entry:
    return kn.parse_entry(rec.get("entry", ""), slug=rec.get("slug", ""))


def _recheck(e: kn.Entry) -> None:
    from . import knowledge_check as kc
    problems = [f.line() for f in kc.check(e.text, public=True, identity={},
                                           unpublished_terms=())]
    problems += kn.template_findings(e, topics=kn.read_topics(
        kn.global_dir()), layer=kn.GLOBAL)
    if problems:
        raise ValueError("the entry does not pass the check: "
                         + "; ".join(problems))


def _decide(rec: dict, status: str, *, reason: str = "") -> Callable:
    def after(branch: str) -> None:
        import getpass
        rec.update(status=status, reason=reason, branch=branch,
                   decided_by=getpass.getuser())
        ks.write_queue_record(rec)
    return after


def plan_accept(qid: str) -> Plan:
    rec = _record(qid)
    e = _queued_entry(rec)
    _recheck(e)
    notes = []
    if rec.get("contradicts"):
        notes.append("contradicts: " + ", ".join(rec["contradicts"])
                     + " -- accepting keeps both; retire or edit the other")
    plan = Plan(f"Accept {qid}: {e.title} (by {rec.get('author')})",
                {_rel(f"{e.slug}.md"): kn.render_entry(e)},
                f"Knowledge: add '{e.title}'\n\nQueue {qid}; provenance "
                + "; ".join(p.render() for p in e.provenance),
                f"accept-{e.slug}", notes=notes)
    plan.after = _decide(rec, "accepted")
    return _finish(plan)


def plan_reject(qid: str, reason: str) -> Plan:
    if not (reason or "").strip():
        raise ValueError("a rejection needs a reason the author can read")
    rec = _record(qid)
    plan = Plan(f"Reject {qid} ({rec.get('author')}): {reason}", {},
                "", "reject", public=False)
    plan.after = _decide(rec, "rejected", reason=reason.strip())
    return plan


def plan_queue_merge(qid: str, into: str) -> Plan:
    rec = _record(qid)
    target = _entry(into)
    merged = kn.merge_entries(target, _queued_entry(rec))
    plan = Plan(f"Merge {qid} into {into} (both provenances kept)",
                {_rel(f"{into}.md"): kn.render_entry(merged)},
                f"Knowledge: extend {into} with queue {qid}",
                f"extend-{into}")
    plan.after = _decide(rec, "merged", reason=f"merged into {into}")
    return _finish(plan)


def plan_queue_edit(qid: str, **fields) -> Plan:
    rec = _record(qid)
    e = _queued_entry(rec)
    e = replace(e, **{k: " ".join(str(v).split()) for k, v in fields.items()
                      if k in ("title", "situation", "lesson", "why",
                               "topic")})
    _recheck(e)
    plan = Plan(f"Accept {qid} with the admin's wording",
                {_rel(f"{e.slug}.md"): kn.render_entry(e)},
                f"Knowledge: add '{e.title}' (edited in review)\n\nQueue "
                f"{qid}", f"accept-{e.slug}")
    plan.after = _decide(rec, "edited", reason="accepted with edits")
    return _finish(plan)


# ---------------------------------------------------------------------------
# run a plan
# ---------------------------------------------------------------------------

def start(plan: Plan, *, clock=None) -> "ks.Dialog":
    """Open the admin's confirmation for ``plan``."""
    def perform() -> tuple:
        if not plan.changes and plan.branch_hint != "index":
            if plan.after:
                plan.after("")
            return (f"done: {plan.title}", "")
        branch, sha, _diff = ks.prepare_branch(
            _repo(), f"delfin-knowledge/{plan.branch_hint}", plan.changes,
            plan.message, index_folder=ks.GLOBAL_FOLDER)
        if plan.after:
            plan.after(branch)
        return ("\n".join([
            f"prepared branch {branch} at {sha[:12]} (nothing pushed, the "
            "working tree untouched)",
            f"review: git show {branch}",
            f"publish when reviewed: git push -u origin {branch} && gh pr "
            f"create --head {branch} --fill"]), branch)
    kwargs = {"clock": clock} if clock is not None else {}
    return ks.start_admin(plan.title, preview_text=plan.preview(),
                          perform=perform, public=plan.public, **kwargs)
