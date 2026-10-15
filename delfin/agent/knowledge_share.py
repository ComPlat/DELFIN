"""Sharing a lesson: the person's act, confirmed by a code only they can read.

Input of :func:`start_share`: a private lesson id (``knowledge_private``),
a layer (``global`` or ``project``), the workspace, an optional scope,
topic and author handle, the snippets the person marks releasable. Output:
a :class:`Dialog` -- the exact entry text that would be written, the check
result with every finding and its span, the target (repository, folder,
scope), its visibility and who will see it, duplicates and contradictions
among the existing entries -- and, when nothing blocks, a one-time code.

The code (``approval_code.issue``) lives only in this process and is
shown only on the person's own channel: the REPL's chrome, the window's
frame, a dashboard notice, the CLI's ``/dev/tty``. :func:`confirm` takes
what the person typed. No agent tool reaches :func:`start_share`,
:func:`confirm` or the writers (pinned by a test over the tool surface),
inside an agent session (``approval_guard.ENV_MARKER`` set) nothing here
starts, and a code from another dialog, a file, a tool result or a web
page answers nothing: it was never there.

Public targets (``repo_visibility``: unknown counts as public) take two
confirmations. The second code is issued only after the first was typed,
together with the full final text, and an answer sooner than
:data:`PAUSE_S` after it is refused as typing, not reading. Private and
internal targets take one.

Where a confirmed share lands:

* GLOBAL (``delfin/agent/pack/knowledge``) -- the admin review queue
  (:func:`queue_dir`), never a branch. The admin accepts, edits, merges
  or rejects (``knowledge_curation``); the author's lesson shows the
  outcome and the reason.
* PROJECT (``<repo>/.delfin/knowledge``) -- a new local branch
  ``delfin-knowledge/<slug>`` committed with git plumbing (a temporary
  index, ``hash-object``/``commit-tree``/``update-ref``): the person's
  working tree and index are untouched, nothing is pushed, and the push
  and PR commands are printed for the person to run.
* Inside the DELFIN tree the project layer IS the global layer; a project
  share there becomes a global one with the scope choice
  (``always``/``delfin-development``).

Checks never rewrite: a finding blocks, the person edits the lesson
(``knowledge_private.edit``) and starts the share again, which checks the
new text. The swarm and personal stores have no path into the writers:
:func:`prepare_branch` is called with entries built here from a lesson
that passed the check, and the stores do not import this module.

Admin (:func:`is_admin`): ``knowledge_admins`` in the managed policy, else
``agent.knowledge_admins`` in the settings, matched against the OS
account and the gh-authenticated login -- never against the author handle
a person chooses. Admin actions use :func:`start_admin`, the same dialog.

Rejected: a y/N answer (a queued or pasted line answers it; a code that
exists only after the question was asked does not); writing into the
working tree and committing there (it would mix with the person's own
changes, and an agent's write could ride along); pushing and opening the
PR here (the person reviews the branch first).
"""
from __future__ import annotations

import datetime as dt
import hmac
import json
import os
import re
import shutil
import stat
import subprocess
import tempfile
import threading
import time
import uuid
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable, Optional

from . import knowledge as kn
from . import knowledge_private as kp
from . import repo_visibility as rv

PAUSE_S = 2.0
MAX_WRONG = 3
DIALOG_TTL_S = 15 * 60
PUBLIC_BANNER = ("PUBLIC: anyone can read this, and it cannot be taken back "
                 "after merging (git history keeps it).")
GLOBAL_FOLDER = "delfin/agent/pack/knowledge"
PROJECT_FOLDER = ".delfin/knowledge"

#: Tests point this at a temporary repository.
_DELFIN_ROOT: Optional[Path] = None
#: Tests replace the gh login probe.
_GH_LOGIN: Optional[Callable[[], str]] = None

_LOCK = threading.Lock()
_DIALOGS: dict = {}


class ShareRefused(ValueError):
    pass


# ---------------------------------------------------------------------------
# where
# ---------------------------------------------------------------------------

def delfin_root() -> Path:
    return Path(_DELFIN_ROOT) if _DELFIN_ROOT else kn.global_dir().parents[3]


def queue_dir() -> Path:
    """``DELFIN_KNOWLEDGE_QUEUE`` > ``agent.knowledge.queue_dir`` >
    ``~/.delfin/knowledge/queue``. A shared folder (group-readable) lets
    an admin on the same host review; the default is this account's."""
    env = os.environ.get("DELFIN_KNOWLEDGE_QUEUE", "").strip()
    if env:
        return Path(env)
    cfg = _knowledge_settings()
    if str(cfg.get("queue_dir") or "").strip():
        return Path(str(cfg["queue_dir"])).expanduser()
    return Path.home() / ".delfin" / "knowledge" / "queue"


def _knowledge_settings() -> dict:
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        cfg = agent.get("knowledge") or {}
        return cfg if isinstance(cfg, dict) else {}
    except Exception:
        return {}


@dataclass
class Target:
    layer: str
    repo_root: Path
    folder_rel: str
    scope: str
    visibility: str
    public: bool
    audience: str
    note: str = ""

    @property
    def folder(self) -> Path:
        return self.repo_root / self.folder_rel


def resolve_target(layer: str, *, workspace, scope: str = "",
                   probe: bool = True) -> Target:
    layer = (layer or "").strip().lower()
    if layer not in (kn.GLOBAL, kn.PROJECT):
        raise ShareRefused("layer must be 'project' or 'global'")
    note = ""
    if layer == kn.PROJECT and kn.is_delfin_repo(workspace):
        layer = kn.GLOBAL
        note = ("inside the DELFIN repository the project layer is the "
                "global layer; choose the scope: always (every DELFIN user) "
                "or delfin-development (only when working on DELFIN itself)")
    if layer == kn.GLOBAL:
        repo = delfin_root()
        folder = GLOBAL_FOLDER
        scope = scope or "always"
        if scope not in kn.SCOPES[kn.GLOBAL]:
            raise ShareRefused("global scope is 'always' or "
                               "'delfin-development'")
    else:
        repo = kn.repo_root_of(workspace)
        folder = PROJECT_FOLDER
        scope = "project"
        if not (repo / ".git").exists():
            raise ShareRefused(f"{repo.name} is not a git repository; the "
                               "project layer lives in one")
    vis = rv.visibility(repo, probe=probe)
    public = rv.treated_public(vis)
    if layer == kn.GLOBAL:
        audience = "every DELFIN user after git pull"
    else:
        audience = f"everyone with access to the repository {repo.name}"
    if public:
        audience += ", and anyone on the internet"
    return Target(layer, repo, folder, scope, vis, public, audience, note)


# ---------------------------------------------------------------------------
# identity and admin
# ---------------------------------------------------------------------------

def _refuse_in_agent_session() -> None:
    try:
        from .approval_guard import ENV_MARKER
    except Exception:
        ENV_MARKER = "DELFIN_AGENT_SESSION"
    if os.environ.get(ENV_MARKER):
        raise ShareRefused("this runs inside an agent session; sharing and "
                           "curating are the person's own acts")


def default_handle() -> str:
    return str(_knowledge_settings().get("handle") or "").strip()


def _policy_admins() -> tuple:
    try:
        from . import managed_policy
        return tuple(getattr(managed_policy.load(), "knowledge_admins", ())
                     or ())
    except Exception:
        return ()


def _settings_admins() -> tuple:
    try:
        from delfin.user_settings import load_settings
        agent = (load_settings() or {}).get("agent") or {}
        raw = agent.get("knowledge_admins") or []
        return tuple(str(x).strip() for x in raw if str(x).strip())
    except Exception:
        return ()


def admins() -> tuple:
    """The managed policy's list when it sets one, else the settings'."""
    return _policy_admins() or _settings_admins()


_LOGIN_CACHE: dict = {}


def _gh_login() -> str:
    if _GH_LOGIN is not None:
        return _GH_LOGIN()
    if "v" in _LOGIN_CACHE:
        return _LOGIN_CACHE["v"]
    try:
        out = subprocess.run(["gh", "api", "user", "-q", ".login"],
                             capture_output=True, text=True, timeout=15)
        login = out.stdout.strip() if out.returncode == 0 else ""
    except (OSError, subprocess.SubprocessError):
        login = ""
    _LOGIN_CACHE["v"] = login
    return login


def _identities() -> set:
    import getpass
    ids = set()
    try:
        ids.add(getpass.getuser())
    except Exception:
        pass
    login = _gh_login()
    if login:
        ids.add(login)
    return {i.lower() for i in ids if i}


def is_admin() -> bool:
    allowed = {a.lower() for a in admins()}
    return bool(allowed & _identities())


# ---------------------------------------------------------------------------
# the dialog
# ---------------------------------------------------------------------------

@dataclass
class Dialog:
    id: str
    kind: str                        # share | unshare | admin
    title: str
    target: Optional[Target]
    entry: Optional[kn.Entry] = None
    text: str = ""                   # the exact text that would be written
    mode: str = "new"                # new | extend | remove | admin
    findings: list = field(default_factory=list)
    released: list = field(default_factory=list)
    contradicts: list = field(default_factory=list)
    extra: list = field(default_factory=list)
    stage: int = 1
    code: str = ""
    wrong: int = 0
    created: float = 0.0
    stage2_at: float = 0.0
    lesson_id: int = 0
    handle: str = ""
    perform: Optional[Callable[[], tuple]] = None
    public: bool = True

    @property
    def blocked(self) -> bool:
        return bool(self.findings)

    @property
    def steps(self) -> int:
        return 2 if self.public else 1

    def _code_rid(self) -> str:
        return f"kn-{self.id}-{self.stage}"

    def preview(self) -> str:
        t = self.target
        lines = [f"{self.title}"]
        if t is not None:
            lines.append(f"target: {t.layer} layer, {t.folder_rel} in "
                         f"{t.repo_root.name}, scope {t.scope}")
            if t.note:
                lines.append(f"note: {t.note}")
        if self.public:
            lines.append(PUBLIC_BANNER)
        elif t is not None:
            lines.append(f"visibility: {t.visibility}")
        if t is not None:
            lines.append(f"who will see it: {t.audience}")
            if t.layer == kn.GLOBAL and self.kind == "share":
                lines.append("route: the knowledge admin's review queue "
                             "first; the admin prepares the branch and PR")
        if self.mode == "extend" and self.entry is not None:
            lines.append(f"duplicate: this would extend the existing entry "
                         f"'{self.entry.slug}' (both provenances kept) "
                         "instead of adding a new file")
        if self.contradicts:
            lines.append("contradicts existing entries: "
                         + ", ".join(self.contradicts)
                         + " -- shown to the admin, who decides")
        lines += self.extra
        if self.released:
            lines.append("released by you: " + ", ".join(
                repr(r) for r in self.released))
        if self.findings:
            lines.append("check: BLOCKED -- " + (
                "change the plan's wording:" if self.kind == "admin" else
                f"reword the lesson (/lessons edit {self.lesson_id} ...) "
                "and share again:"))
            lines += [f"  - {f}" for f in self.findings]
        elif self.kind == "admin":
            lines.append("check: passed (" + ("public" if self.public
                                              else "base") + " rules) for "
                         "every entry this change writes")
        else:
            lines.append("check: passed (" + ("public" if self.public
                                              else "base") + " rules)")
        if self.text:
            lines += ["--- exact text ---", self.text.rstrip(),
                      "--- end of text ---"]
        if self.code:
            from .approval_code import display
            step = f"step {self.stage} of {self.steps}"
            if self.stage == 2:
                lines.append("Second confirmation for a public target: read "
                             "the full final text above, then type the code "
                             f"{display(self.code)} ({step}); anything else "
                             "cancels.")
            else:
                lines.append(f"Type the code {display(self.code)} to confirm "
                             f"({step}); anything else cancels.")
        return "\n".join(lines)


@dataclass
class Outcome:
    status: str                      # next | done | refused | expired
    message: str = ""
    dialog: Optional[Dialog] = None
    ref: str = ""


def _issue(d: Dialog) -> None:
    from .approval_code import issue
    d.code = issue(d._code_rid())


def _forget(d: Dialog) -> None:
    try:
        from .approval_code import forget
        forget(d._code_rid())
    except Exception:
        pass


def _register(d: Dialog, clock) -> Dialog:
    d.created = clock()
    if not d.findings:
        _issue(d)
    with _LOCK:
        _DIALOGS[d.id] = d
    return d


def get_dialog(dialog_id: str) -> Optional[Dialog]:
    with _LOCK:
        return _DIALOGS.get(str(dialog_id))


def open_dialogs() -> list:
    with _LOCK:
        return list(_DIALOGS.values())


def cancel(dialog_id: str) -> bool:
    with _LOCK:
        d = _DIALOGS.pop(str(dialog_id), None)
    if d is not None:
        _forget(d)
    return d is not None


def confirm(dialog_id: str, typed: str, *,
            clock: Callable[[], float] = time.monotonic) -> Outcome:
    """Answer one dialog with what the person typed."""
    d = get_dialog(dialog_id)
    if d is None:
        return Outcome("expired", "no such dialog (answered, cancelled or "
                       "expired); start it again")
    if d.blocked:
        return Outcome("refused", "this dialog is blocked by the check; "
                       "reword and start again", d)
    if clock() - d.created > DIALOG_TTL_S:
        cancel(d.id)
        return Outcome("expired", "the dialog expired; start it again")
    if d.stage == 2 and clock() < d.stage2_at + PAUSE_S:
        return Outcome("refused", f"read the full text first; type the code "
                       f"again in {PAUSE_S:.0f} s", d)
    from .approval_code import normalize
    if not d.code or not hmac.compare_digest(normalize(typed),
                                             normalize(d.code)):
        d.wrong += 1
        if d.wrong >= MAX_WRONG:
            cancel(d.id)
            return Outcome("expired", "three wrong codes: the dialog is "
                           "closed and nothing was written")
        return Outcome("refused", "that is not the code; nothing was "
                       "written", d)
    _forget(d)
    if d.public and d.stage == 1:
        d.stage = 2
        d.wrong = 0
        _issue(d)
        d.stage2_at = clock()
        return Outcome("next", d.preview(), d)
    with _LOCK:
        _DIALOGS.pop(d.id, None)
    try:
        message, ref = d.perform() if d.perform else ("nothing to do", "")
    except ShareRefused as exc:
        return Outcome("refused", str(exc))
    except Exception as exc:
        return Outcome("refused", f"failed, nothing was written: {exc}")
    return Outcome("done", message, None, ref)


# ---------------------------------------------------------------------------
# share
# ---------------------------------------------------------------------------

def _project_names(les: kp.Lesson) -> list:
    ws = les.workspace
    if not ws or kn.is_delfin_repo(ws):
        return []
    root = kn.repo_root_of(ws)
    names = [root.name]
    try:
        url = subprocess.run(["git", "-C", str(root), "remote", "get-url",
                              "origin"], capture_output=True, text=True,
                             timeout=5).stdout.strip()
        if url:
            names.append(url.rstrip("/").rsplit("/", 1)[-1]
                         .removesuffix(".git"))
    except (OSError, subprocess.SubprocessError):
        pass
    return [n for n in dict.fromkeys(names) if n]


def build_entry(les: kp.Lesson, target: Target, handle: str,
                *, today: Optional[dt.date] = None) -> kn.Entry:
    today = today or dt.date.today()
    prov = [kn.Provenance(handle, today.isoformat(), ev)
            for ev in les.evidence[:3]]
    return kn.Entry(
        slug=kn.slugify(les.title), title=les.title, topic=les.topic,
        scope=target.scope, provenance=prov,
        review_by=(today + dt.timedelta(days=kn.DEFAULT_REVIEW_DAYS))
        .isoformat(),
        situation=les.situation, lesson=les.lesson, why=les.why,
        layer=target.layer)


_OPAQUE_POINTER_RE = re.compile(r"^(?:commit:[0-9a-f]{7,40}|pr:\d{1,7}"
                                r"|issue:\d{1,7})$")


def content_findings(entry: kn.Entry, *, public: bool, handle: str = "",
                     names: Iterable[str] = (), releasable=(),
                     identity: Optional[dict] = None,
                     unpublished_terms=None) -> list[str]:
    """What the content rules find in an entry: its text, its authors and
    evidence pointers, and -- for a public target -- every pointer that is
    not one a reader elsewhere can follow. The lint calls this with no
    machine identity and no personal terms, so it answers the same on
    every machine."""
    from . import knowledge_check as kc
    kw = dict(public=public, names=names, releasable=releasable,
              identity=identity, unpublished_terms=unpublished_terms)
    out = [f.line() for f in kc.check(entry.text, **kw)]
    authors = [handle] if handle else [p.author for p in entry.provenance]
    # A commit hash or a PR/issue number is a public identifier and its
    # digits are not a value or a node name; every other pointer is read.
    meta = " ".join(authors + [
        p.evidence for p in entry.provenance
        if not _OPAQUE_POINTER_RE.match(p.evidence)])
    out += [f"provenance: {f.line()}" for f in kc.check(meta, **kw)
            if f.rule != "user-specific"]
    if public:
        for p in entry.provenance:
            kind = p.evidence.split(":", 1)[0]
            if kind not in kn.PUBLIC_EVIDENCE_KINDS:
                out.append(
                    f"provenance: evidence {p.evidence!r} points at this "
                    "machine; a public entry needs a public pointer ("
                    + ", ".join(f"{k}:" for k in kn.PUBLIC_EVIDENCE_KINDS)
                    + "): /lessons edit <id> evidence=<pointer>")
    return out


def entry_findings(entry: kn.Entry, target: Target, *, handle: str,
                   names: Iterable[str] = (), releasable=(),
                   topics: Optional[dict] = None) -> list[str]:
    out = content_findings(entry, public=target.public, handle=handle,
                           names=names, releasable=releasable)
    if topics is None and target.layer == kn.GLOBAL:
        topics = kn.read_topics(target.folder)
    out += kn.template_findings(entry, topics=topics, layer=target.layer)
    return out


def start_share(lesson_id, layer: str, *, workspace, scope: str = "",
                topic: str = "", handle: str = "", releasable=(),
                clock: Callable[[], float] = time.monotonic) -> Dialog:
    _refuse_in_agent_session()
    les = kp.get(lesson_id)
    if les is None:
        raise ShareRefused(f"no lesson {lesson_id}; /memory lists them")
    if les.status in (kp.QUEUED, kp.PREPARED):
        raise ShareRefused(f"lesson {les.id} is already {les.status} "
                           f"({les.shared_ref})")
    if les.status == kp.WITHDRAWN:
        raise ShareRefused(f"lesson {les.id} was withdrawn")
    if topic:
        les = kp.edit(les.id, topic=topic)
    handle = (handle or default_handle()).strip()
    if not handle:
        raise ShareRefused("choose an author handle for the provenance: "
                           "/share <id> <layer> as <handle>, or set "
                           "agent.knowledge.handle")
    target = resolve_target(layer, workspace=workspace, scope=scope)
    entry = build_entry(les, target, handle)
    names = _project_names(les) if target.public else []
    findings = entry_findings(entry, target, handle=handle, names=names,
                              releasable=releasable)
    if target.layer == kn.GLOBAL:
        topics = kn.read_topics(target.folder)
        if les.topic and les.topic not in topics:
            findings.append(f"topic {les.topic!r} is not in the taxonomy; "
                            "only an admin adds a topic")
    existing = kn.read_entries(target.folder, layer=target.layer)
    mode, written = "new", entry
    dup = next((e for e in existing if kn.near_duplicate(entry, e)), None)
    if dup is not None:
        mode = "extend"
        written = kn.merge_entries(dup, entry)
        written.layer = target.layer
    contradictions = [e.slug for e in existing if kn.contradicts(entry, e)]
    d = Dialog(
        id=uuid.uuid4().hex[:12], kind="share",
        title=f"Share lesson {les.id} -> {target.layer}",
        target=target, entry=written, text=kn.render_entry(written),
        mode=mode, findings=findings, contradicts=contradictions,
        released=[str(r) for r in releasable or ()],
        lesson_id=les.id, handle=handle, public=target.public)

    def perform() -> tuple:
        if target.layer == kn.GLOBAL:
            qid = enqueue(d)
            kp.set_status(les.id, kp.QUEUED, shared_ref=qid)
            return (f"queued for the knowledge admin's review as {qid}; "
                    "/memory shows the outcome", qid)
        branch, sha, _diff = prepare_branch(
            target.repo_root, f"delfin-knowledge/{written.slug}",
            {f"{target.folder_rel}/{written.slug}.md": d.text},
            _commit_message(d))
        kp.set_status(les.id, kp.PREPARED, shared_ref=branch)
        return (_offer_pr(target, branch, sha, d), branch)
    d.perform = perform
    return _register(d, clock)


def _commit_message(d: Dialog) -> str:
    e = d.entry
    verb = {"extend": "Extend", "remove": "Withdraw"}.get(d.mode, "Add")
    lines = [f"Knowledge: {verb} '{e.title}'",
             "", f"Layer {d.target.layer}, topic {e.topic}, scope "
             f"{e.scope}.", "Provenance: " + "; ".join(
                 p.render() for p in e.provenance)]
    if d.contradicts:
        lines.append("Contradicts: " + ", ".join(d.contradicts))
    return "\n".join(lines)


def _offer_pr(target: Target, branch: str, sha: str, d: Dialog) -> str:
    folder = target.folder_rel
    return "\n".join([
        f"prepared branch {branch} at {sha[:12]} in {target.repo_root.name} "
        "(working tree and index untouched; nothing pushed)",
        f"review it: git -C {target.repo_root.name} show {branch}",
        f"then, if you want the pull request: git push -u origin {branch} "
        f"&& gh pr create --head {branch} --fill",
        f"suggested for the PR: guard {folder}/ in CODEOWNERS with the "
        "knowledge maintainers",
    ])


# ---------------------------------------------------------------------------
# the admin review queue
# ---------------------------------------------------------------------------

def enqueue(d: Dialog) -> str:
    q = queue_dir()
    shared = _queue_is_shared()
    q.mkdir(parents=True, exist_ok=True)
    if not shared:
        os.chmod(q, 0o700)
    qid = time.strftime("%Y%m%d-%H%M%S-") + uuid.uuid4().hex[:6]
    record = {
        "id": qid, "created": time.time(), "author": d.handle,
        "lesson_id": d.lesson_id, "slug": d.entry.slug,
        "topic": d.entry.topic, "scope": d.entry.scope, "mode": d.mode,
        "entry": d.text, "contradicts": list(d.contradicts),
        "released": list(d.released), "visibility": d.target.visibility,
        "status": "pending", "reason": "", "decided_by": "",
    }
    _write_record(q / f"{qid}.json", record, shared=shared)
    return qid


def _queue_is_shared() -> bool:
    return bool(os.environ.get("DELFIN_KNOWLEDGE_QUEUE")
                or _knowledge_settings().get("queue_dir"))


def _record_mode(shared: bool) -> int:
    """The mode of a queue record: 0o600, plus group read in a shared
    queue (the admin of the same host reads it through the group)."""
    return 0o640 if shared else 0o600


def _write_record(path: Path, data: dict, *, shared: bool,
                  cap: Optional[int] = None) -> None:
    """Write one queue record with an exact mode, atomically.

    Input: the final path, the record, whether the queue is shared, and
    ``cap``, the mode bits of the record being replaced (a rewrite never
    grants a bit the old record did not have). Output: ``path`` holds the
    record, mode ``_record_mode(shared) & cap``.

    Semantics: the bytes go to a fresh temporary file in the same folder
    (``O_CREAT | O_EXCL | O_NOFOLLOW``, so no existing file or planted
    link is opened), whose mode is set with ``fchmod`` before the first
    byte -- independent of the umask --, then fsync, then ``os.replace``
    onto ``path``, then a best-effort fsync of the folder. A reader sees
    the old record or the new one, never a partial one, and never a
    moment with a wider mode. On any failure the temporary file is
    removed and the old record is untouched.

    Rejected: ``Path.write_text`` + ``os.replace`` (the mode was 0o666
    minus the rewriter's umask: 0o644 under the common 022, so an admin
    decision made a record world-readable); ``chmod`` after the write (a
    window in which the wider mode exists); the umask-masked mode
    argument of ``os.open`` alone (under umask 077 a shared record was
    unreadable to the admin).
    """
    mode = _record_mode(shared)
    if cap is not None:
        mode &= (cap & 0o777) | 0o600
    payload = json.dumps(data, ensure_ascii=False, indent=1).encode("utf-8")
    tmp = path.with_name(f".{path.name}.{uuid.uuid4().hex[:8]}.partial")
    flags = (os.O_WRONLY | os.O_CREAT | os.O_EXCL
             | getattr(os, "O_NOFOLLOW", 0) | getattr(os, "O_CLOEXEC", 0))
    fd = os.open(str(tmp), flags, 0o600)
    try:
        if hasattr(os, "fchmod"):
            os.fchmod(fd, mode)
        view = memoryview(payload)
        while view:
            view = view[os.write(fd, view):]
        os.fsync(fd)
    except BaseException:
        os.close(fd)
        tmp.unlink(missing_ok=True)
        raise
    os.close(fd)
    try:
        os.replace(tmp, path)
    except BaseException:
        tmp.unlink(missing_ok=True)
        raise
    try:
        dfd = os.open(str(path.parent), os.O_RDONLY)
    except OSError:
        return
    try:
        os.fsync(dfd)
    except OSError:
        pass
    finally:
        os.close(dfd)


def read_queue() -> list:
    q = queue_dir()
    out = []
    for p in sorted(q.glob("*.json")) if q.is_dir() else []:
        if p.is_symlink():
            continue
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if isinstance(rec, dict) and rec.get("id") == p.stem:
            rec["_path"] = str(p)
            out.append(rec)
    return out


def write_queue_record(rec: dict) -> None:
    """Rewrite one record (an admin decision) through
    :func:`_write_record`, never wider than the record it replaces."""
    p = Path(rec["_path"])
    data = {k: v for k, v in rec.items() if not k.startswith("_")}
    try:
        st = os.lstat(p)
    except FileNotFoundError:
        st = None
    if st is not None and not stat.S_ISREG(st.st_mode):
        raise ShareRefused(f"{p.name} is not a regular file; not rewritten")
    _write_record(p, data, shared=_queue_is_shared(),
                  cap=None if st is None else stat.S_IMODE(st.st_mode))


def sync_own_lessons() -> list:
    """Carry the admin's decisions back into this account's lessons.
    Returns one line per lesson whose status changed."""
    by_id = {r.get("id"): r for r in read_queue()}
    lines = []
    for les in kp.list_lessons():
        if les.status != kp.QUEUED or les.shared_ref not in by_id:
            continue
        rec = by_id[les.shared_ref]
        st = rec.get("status")
        if st == "rejected":
            kp.set_status(les.id, kp.REJECTED, note=str(rec.get("reason")))
            lines.append(f"lesson {les.id} was rejected: {rec.get('reason')}")
        elif st in ("accepted", "merged", "edited"):
            kp.set_status(les.id, kp.PREPARED,
                          note=str(rec.get("reason") or st),
                          shared_ref=str(rec.get("branch") or les.shared_ref))
            lines.append(f"lesson {les.id} was {st} by the admin")
    return lines


# ---------------------------------------------------------------------------
# unshare
# ---------------------------------------------------------------------------

def _find_entry(slug: str, workspace) -> tuple:
    candidates = []
    if workspace and not kn.is_delfin_repo(workspace):
        root = kn.repo_root_of(workspace)
        candidates.append((kn.PROJECT, root, PROJECT_FOLDER))
    candidates.append((kn.GLOBAL, delfin_root(), GLOBAL_FOLDER))
    for layer, root, folder in candidates:
        for e in kn.read_entries(root / folder, layer=layer):
            if e.slug == slug:
                return layer, root, folder, e
    raise ShareRefused(f"no shared entry named {slug!r}")


def start_unshare(slug: str, *, workspace, handle: str = "",
                  clock: Callable[[], float] = time.monotonic) -> Dialog:
    """Withdraw one's own entry as a branch that removes it."""
    _refuse_in_agent_session()
    handle = (handle or default_handle()).strip()
    layer, root, folder, entry = _find_entry(slug, workspace)
    authors = {p.author for p in entry.provenance}
    if handle not in authors:
        raise ShareRefused(f"{slug!r} is not yours (authors: "
                           f"{', '.join(sorted(authors))}); only the knowledge "
                           "admin retires another person's entry")
    target = resolve_target(layer if layer == kn.GLOBAL else kn.PROJECT,
                            workspace=root if layer == kn.PROJECT
                            else workspace, scope=entry.scope
                            if layer == kn.GLOBAL else "")
    rel = f"{folder}/{slug}.md"
    rest = [p for p in entry.provenance if p.author != handle]
    if rest:
        kept = kn.Entry(**{**entry.__dict__, "provenance": rest})
        change = {rel: kn.render_entry(kept)}
        what = "removes your provenance; the co-authors' entry stays"
    else:
        change = {rel: None}
        what = "removes the entry file"
    d = Dialog(id=uuid.uuid4().hex[:12], kind="unshare",
               title=f"Withdraw {slug} from the {layer} layer ({what})",
               target=target, entry=entry, text=kn.render_entry(entry),
               mode="remove", handle=handle, public=target.public,
               extra=["a removal does not erase the text from the git history"
                      if target.public else ""])
    d.extra = [x for x in d.extra if x]

    def perform() -> tuple:
        branch, sha, _ = prepare_branch(
            root, f"delfin-knowledge/withdraw-{slug}", change,
            f"Knowledge: withdraw '{entry.title}'\n\nRequested by "
            f"{handle}.", index_folder=folder if layer == kn.GLOBAL
            else None)
        return (_offer_pr(target, branch, sha, d), branch)
    d.perform = perform
    return _register(d, clock)


# ---------------------------------------------------------------------------
# admin actions
# ---------------------------------------------------------------------------

def start_admin(title: str, *, preview_text: str,
                perform: Callable[[], tuple], public: bool = True,
                target: Optional[Target] = None, findings=(),
                clock: Callable[[], float] = time.monotonic) -> Dialog:
    """The dialog behind every curation action. Refused for a non-admin,
    and inside an agent session."""
    _refuse_in_agent_session()
    if not is_admin():
        raise ShareRefused("curation is for the knowledge admin "
                           "(knowledge_admins in the managed policy or "
                           "agent.knowledge_admins); you can propose, share "
                           "and withdraw your own entries")
    d = Dialog(id=uuid.uuid4().hex[:12], kind="admin", title=title,
               target=target, text=preview_text, mode="admin",
               public=public, perform=perform, findings=list(findings))
    return _register(d, clock)


# ---------------------------------------------------------------------------
# the branch writer: git plumbing, no working tree, no push
# ---------------------------------------------------------------------------

def _g(repo: Path, *args, env=None, input_text=None) -> str:
    out = subprocess.run(["git", "-C", str(repo), *args],
                         capture_output=True, text=True, env=env,
                         input=input_text, timeout=60)
    if out.returncode != 0:
        raise ShareRefused(f"git {args[0]} failed: "
                           f"{(out.stderr or out.stdout).strip()[:300]}")
    return out.stdout.strip()


def base_ref(repo: Path) -> str:
    for ref in ("origin/main", "origin/master", "HEAD"):
        try:
            return _g(repo, "rev-parse", "--verify", "-q",
                      f"{ref}^{{commit}}")
        except ShareRefused:
            continue
    raise ShareRefused("the repository has no commit to branch from")


def _materialise(repo: Path, base: str, folder: str, dst: Path) -> None:
    names = _g(repo, "ls-tree", "-r", "--name-only", base, "--", folder)
    for name in names.splitlines():
        if not name.endswith(".md"):
            continue
        p = dst / Path(name).relative_to(folder)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(_g(repo, "show", f"{base}:{name}") + "\n",
                     encoding="utf-8")


def prepare_branch(repo: Path | str, branch_hint: str, changes: dict,
                   message: str, *, index_folder: Optional[str] = None,
                   base: str = "") -> tuple:
    """Commit ``changes`` ({path: text or None to delete}) on a new branch.

    Output: (branch, commit sha, diff against the base). A temporary index
    is used, so the working tree and the person's index stay as they were;
    nothing is pushed. With ``index_folder`` the folder's INDEX.md is
    regenerated from the entries as they are on the new branch."""
    repo = Path(repo)
    base = base or base_ref(repo)
    changes = dict(changes)
    tmp = Path(tempfile.mkdtemp(prefix="delfin-kn-"))
    try:
        if index_folder:
            mat = tmp / "folder"
            mat.mkdir()
            _materialise(repo, base, index_folder, mat)
            for rel, text in changes.items():
                if not rel.startswith(index_folder + "/"):
                    continue
                p = mat / Path(rel).relative_to(index_folder)
                if text is None:
                    p.unlink(missing_ok=True)
                else:
                    p.parent.mkdir(parents=True, exist_ok=True)
                    p.write_text(text, encoding="utf-8")
            changes[f"{index_folder}/{kn.INDEX_FILE}"] = kn.render_index(mat)
        env = dict(os.environ, GIT_INDEX_FILE=str(tmp / "index"))
        _g(repo, "read-tree", base, env=env)
        for rel, text in changes.items():
            if text is None:
                _g(repo, "update-index", "--force-remove", "--", rel, env=env)
                continue
            blob = _g(repo, "hash-object", "-w", "--stdin", env=env,
                      input_text=text)
            _g(repo, "update-index", "--add", "--cacheinfo",
               f"100644,{blob},{rel}", env=env)
        tree = _g(repo, "write-tree", env=env)
        sha = _g(repo, "commit-tree", tree, "-p", base, "-m", message,
                 env=env)
        name = branch_hint
        i = 2
        while True:
            try:
                _g(repo, "rev-parse", "--verify", "-q",
                   f"refs/heads/{name}")
            except ShareRefused:
                break
            name = f"{branch_hint}-{i}"
            i += 1
        _g(repo, "update-ref", f"refs/heads/{name}", sha, "0" * 40)
        diff = _g(repo, "diff", base, sha)
        return name, sha, diff
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def preview_diff(repo: Path | str, changes: dict, *,
                 index_folder: Optional[str] = None) -> str:
    """The diff ``prepare_branch`` would commit, without keeping a ref."""
    repo = Path(repo)
    branch, sha, diff = prepare_branch(
        repo, f"delfin-knowledge/preview-{uuid.uuid4().hex[:6]}", changes,
        "preview", index_folder=index_folder)
    try:
        _g(repo, "update-ref", "-d", f"refs/heads/{branch}", sha)
    except ShareRefused:
        pass
    return diff
