"""Test a branch in a second dashboard, without moving the shared checkout.

Input: the shared checkout (the one the running dashboard imports DELFIN
from), a ref (a local branch, ``origin/<branch>``, a tag or a commit), the
interpreter of the running dashboard and its environment. Output: a
:class:`Preview` record -- port, token, address, worktree, the verified
``delfin.__file__`` -- or the reason it was refused.

What runs where:

* **Worktree.** ``<parent of the main worktree>/delfin-preview/<slug>``,
  outside every checkout. A free branch is checked out by name there; a
  branch another worktree holds (the shared checkout holds ``main``) is
  checked out detached at its head (``held_branch.held_elsewhere``). An
  existing preview worktree is reused: moved to the ref's current commit
  when it is clean, used as it is (and said so) when it is not. No
  ``git stash``, ``git switch`` or ``git checkout`` runs in the shared
  checkout.
* **Interpreter.** The running dashboard's own (``sys.executable``), with
  ``PYTHONPATH=<worktree>`` and nothing else on it, started in
  ``<worktree>/delfin_voila_runtime`` where its kernels start. Before the
  server starts, ``import delfin; print(delfin.__file__)`` runs there with
  exactly that environment; a file outside the worktree refuses the
  preview, because a worktree isolates files, not imports.
* **Server.** ``python -m delfin.cli_voila --no-browser --port <free>
  --strict-port`` from the worktree: the branch's own launcher, bound to
  127.0.0.1 with its own generated token (read back from the 0600 token
  file the launcher names), in a process group of its own, put on the
  running dashboard's lifeline ledger so it ends with it.
* **Environment.** The dashboard's, minus the server's token
  (``JUPYTER_TOKEN``), kernel and lifeline identifiers, kept-session and
  terminal switches, browser hand-offs (``BROWSER``, ``VSCODE_*``) and every
  provider key or variable named ``*_KEY``/``*_TOKEN``/``*SECRET*``/
  ``*PASSWORD*``; plus ``DELFIN_NO_BROWSER=1``, the preview markers and
  ``DELFIN_LAUNCH_CWD=<worktree>`` (the agent of a preview works in the
  preview's checkout). ``HOME`` and the settings are the real ones, as for
  the normal dashboard.
* **Browser.** Nothing is opened on the server; the dashboard shows a link.

:func:`stop` sends SIGTERM to the process group (SIGKILL after 20 s), ends
what the preview's agents put on its ledger, and removes the worktree when
it has no uncommitted changes and no commits on no branch; otherwise it
asks (``needs_confirm``), and a confirmed discard logs the head commit to
the clean-up log with the command that brings it back.

Rejected: a per-preview ``pip install -e`` venv (minutes per preview,
needs the network, and an editable install is one more resolver between
the worktree and the import that the check above already pins); switching
the shared checkout (moves every person and agent on it, and a kept agent
session keeps the old modules in memory).
"""

from __future__ import annotations

import json
import os
import re
import signal
import socket
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

from delfin.dashboard import branch_list as BL

PREVIEW_ENV_BRANCH = "DELFIN_PREVIEW_BRANCH"
PREVIEW_ENV_SHA = "DELFIN_PREVIEW_SHA"
PREVIEW_ENV_WORKTREE = "DELFIN_PREVIEW_WORKTREE"
SHARED_CHECKOUT_ENV = "DELFIN_SHARED_CHECKOUT"
PREVIEW_DIR_NAME = "delfin-preview"
RUNTIME_DIR = "delfin_voila_runtime"
NOTEBOOK = "delfin_dashboard.ipynb"

#: Removed from the environment handed to a preview server.
_DROP_EXACT = {
    "JUPYTER_TOKEN", "DELFIN_VOILA_PORT", "DELFIN_VOILA_ROOT_DIR",
    "DELFIN_LIFELINE_PID", "DELFIN_LIFELINE_TICKS", "DELFIN_RESUME_SESSION",
    "DELFIN_RESUME_URL_PATH", "DELFIN_CLAUDE_TERMINAL",
    "DELFIN_CLAUDE_TERMINAL_ARGV", "DELFIN_CLAUDE_TERMINAL_CWD",
    "DELFIN_KETCHER_URL", "DELFIN_KEEP_SESSIONS", "BROWSER", "TERM_PROGRAM",
    "TERM_PROGRAM_VERSION", "MPLBACKEND", "PYTHONPATH", "PYTHONHOME",
    "PYTHONSTARTUP", "JUPYTER_CONFIG_DIR", "JUPYTER_RUNTIME_DIR",
}
_DROP_PREFIX = ("JPY_", "VSCODE_", "KERNEL_", "DELFIN_PREVIEW_")
_SECRET_RE = re.compile(r"(_KEY|_TOKEN|SECRET|PASSWORD|PASSWD)", re.I)


@dataclass
class Target:
    ref: str
    sha: str
    branch: str     # checked out by name when free; "" for a detached ref
    kind: str


@dataclass
class Preview:
    slug: str
    ref: str
    sha: str
    branch: str
    repo: str
    worktree: str
    port: int
    token: str
    url_path: str
    pid: int
    proc_start: str
    ticks: Optional[int]
    host: str
    started_at: float
    delfin_file: str
    log: str
    attached: bool = False

    def url(self, host: str = "localhost", scheme: str = "http") -> str:
        return f"{scheme}://{host}:{self.port}{self.url_path}?token={self.token}"


# -- places ------------------------------------------------------------------

def _registry_dir() -> Path:
    """Resolved per call, so a redirected home moves it."""
    return Path.home() / ".delfin" / "branch_previews"


def slug(ref: str) -> str:
    s = re.sub(r"[^A-Za-z0-9._-]+", "-", ref or "").strip("-.")
    return (s or "preview")[:80]


def main_worktree(repo) -> Path:
    from delfin.dashboard.branch_cleanup import worktrees
    rows = worktrees(repo)
    if rows:
        return Path(rows[0]["worktree"])
    return Path(BL._out(repo, "rev-parse", "--show-toplevel") or repo)


def preview_root(repo) -> Path:
    return main_worktree(repo).parent / PREVIEW_DIR_NAME


def worktree_path(repo, ref: str) -> Path:
    return preview_root(repo) / slug(ref)


def is_preview_worktree(repo, path) -> bool:
    try:
        return Path(path).resolve().parent == preview_root(repo).resolve()
    except OSError:
        return False


# -- the ref and the worktree ------------------------------------------------

def resolve_target(repo, ref: str) -> tuple[Optional[Target], str]:
    ref = (ref or "").strip()
    if not BL.valid_ref(ref):
        return None, f"Not a branch, tag or commit name: {ref!r}"
    kind = BL.ref_kind(repo, ref)
    if not kind and ref.startswith("origin/"):
        if BL.fetch_remote_branch(repo, ref[len("origin/"):]):
            kind = "remote"
    if not kind:
        return None, f"{ref}: not found."
    if kind == "remote" and ref.startswith("origin/"):
        BL.fetch_remote_branch(repo, ref[len("origin/"):])
    sha = BL._out(repo, "rev-parse", "--verify", f"{ref}^{{commit}}")
    if not sha:
        return None, f"{ref}: does not name a commit."
    return Target(ref=ref, sha=sha, branch=ref if kind in ("local", "main") else "",
                  kind=kind), ""


def _worktree_row(repo, path: Path) -> Optional[dict]:
    from delfin.dashboard.branch_cleanup import worktrees, _same
    for row in worktrees(repo):
        if _same(row["worktree"], path):
            return row
    return None


def changes(path) -> str:
    """Uncommitted changes plus commits on no branch, as text ("" = clean)."""
    from delfin.dashboard.branch_cleanup import _changes, _commits_on_no_branch
    if not Path(path).is_dir():
        return ""
    out = _changes(str(path))
    lonely = 0
    if not BL._out(path, "branch", "--show-current"):
        lonely = _commits_on_no_branch(str(path))
    if lonely:
        out = (out + "\n" if out else "") + f"{lonely} commit(s) on no branch"
    return out


def ensure_worktree(repo, target: Target) -> tuple[bool, str, Path, bool]:
    """``(ok, message, path, attached)`` -- see the module docstring."""
    from delfin.agent import held_branch as HB
    path = worktree_path(repo, target.ref)
    holder = HB.held_elsewhere(repo, target.branch) if target.branch else None
    # held_elsewhere answers "another worktree than *repo*"; the shared
    # checkout itself holding the branch counts as held too.
    if target.branch and not holder and BL.current_branch(repo) == target.branch:
        holder = str(main_worktree(repo))
    row = _worktree_row(repo, path)
    if row is not None and holder and Path(holder).resolve() == path.resolve():
        holder = None                                 # held by this preview
    if path.exists() and row is None:
        return False, (f"{path} exists and is not a worktree of this "
                       "repository; nothing was changed."), path, False
    if row is not None:
        dirty = changes(path)
        if dirty:
            return True, (f"Reusing {path} as it is: it has uncommitted "
                          f"changes ({dirty.splitlines()[0]} ...)."), path, \
                bool(row.get("branch"))
        if target.branch and not holder:
            r = BL._git(path, "checkout", "--quiet", target.branch)
            attached = True
        else:
            r = BL._git(path, "checkout", "--quiet", "--detach", target.sha)
            attached = False
        if r.returncode != 0:
            return False, (r.stderr or r.stdout).strip(), path, False
        return True, f"Reusing {path} at {target.ref}.", path, attached
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return False, f"Cannot create {path.parent}: {exc}", path, False
    if target.branch and not holder:
        r = BL._git(repo, "worktree", "add", "--quiet", str(path),
                    target.branch, timeout=300.0)
        attached = True
    else:
        r = BL._git(repo, "worktree", "add", "--quiet", "--detach", str(path),
                    target.sha, timeout=300.0)
        attached = False
    if r.returncode != 0:
        return False, (r.stderr or r.stdout).strip() or "git worktree add failed", \
            path, False
    note = (f" ({target.branch} is checked out in {holder}, so this one is "
            "detached at its head)") if holder else ""
    return True, f"Created {path} at {target.ref}{note}.", path, attached


# -- environment and import check --------------------------------------------

def preview_env(base_env: dict, worktree, target: Target,
                shared_checkout: str = "") -> dict:
    env = {}
    for k, v in (base_env or {}).items():
        if k in _DROP_EXACT or k.startswith(_DROP_PREFIX):
            continue
        if _SECRET_RE.search(k):
            continue
        env[k] = v
    try:
        from delfin.agent.process_guard import PROVIDER_KEYS
        for k in PROVIDER_KEYS:
            env.pop(k, None)
    except Exception:
        pass
    wt = str(Path(worktree).resolve())
    env.update({
        "PYTHONPATH": wt,
        "PYTHONUNBUFFERED": "1",
        "DELFIN_NO_BROWSER": "1",
        "DELFIN_LAUNCH_CWD": wt,
        "DELFIN_VOILA_ROOT_DIR": wt,
        PREVIEW_ENV_BRANCH: target.ref,
        PREVIEW_ENV_SHA: target.sha,
        PREVIEW_ENV_WORKTREE: wt,
    })
    if shared_checkout:
        env[SHARED_CHECKOUT_ENV] = str(shared_checkout)
    return env


def verify_import(python: str, env: dict, worktree) -> tuple[bool, str, str]:
    """``(ok, delfin.__file__, message)`` as a preview kernel would import it."""
    wt = Path(worktree).resolve()
    cwd = wt / RUNTIME_DIR
    try:
        cwd.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return False, "", f"Cannot create {cwd}: {exc}"
    try:
        r = subprocess.run(
            [python, "-c", "import delfin; print(delfin.__file__)"],
            cwd=str(cwd), env=env, capture_output=True, text=True, timeout=120)
    except (OSError, subprocess.TimeoutExpired) as exc:
        return False, "", f"The import check could not run: {exc}"
    file = (r.stdout or "").strip().splitlines()[-1:] or [""]
    file = file[0]
    if r.returncode != 0 or not file:
        return False, "", ("The import check failed: "
                           + (r.stderr or r.stdout or "").strip()[-400:])
    try:
        inside = wt in Path(file).resolve().parents
    except OSError:
        inside = False
    if not inside:
        return False, file, (f"Refused: delfin would be imported from {file}, "
                             f"not from the preview worktree {wt}.")
    return True, file, f"delfin.__file__ = {file}"


# -- the registry ------------------------------------------------------------

def _record_path(s: str) -> Path:
    return _registry_dir() / f"{s}.json"


def _save(p: Preview) -> None:
    from delfin.agent.state_paths import write_text_atomic
    write_text_atomic(_record_path(p.slug), json.dumps(asdict(p), indent=1))


def _alive(p: Preview) -> bool:
    _reap(p.pid)
    try:
        from delfin.agent import proc_identity
        return bool(proc_identity.alive(p.pid, p.proc_start, p.host))
    except Exception:
        return False


def load(s: str) -> Optional[Preview]:
    try:
        data = json.loads(_record_path(s).read_text(encoding="utf-8"))
        return Preview(**data)
    except (OSError, ValueError, TypeError):
        return None


def list_previews(repo=None) -> list[Preview]:
    """Running previews (of *repo* when given); dead records are dropped."""
    out: list[Preview] = []
    try:
        files = sorted(_registry_dir().glob("*.json"))
    except OSError:
        return out
    key = BL.repo_key(repo) if repo else ""
    for f in files:
        p = load(f.stem)
        if p is None:
            continue
        if p.host != socket.gethostname():
            continue
        if not _alive(p):
            try:
                f.unlink()
            except OSError:
                pass
            continue
        if key and BL.repo_key(p.repo) != key:
            continue
        out.append(p)
    return out


def running_paths(repo=None) -> list[str]:
    return [p.worktree for p in list_previews(repo)]


# -- start / stop ------------------------------------------------------------

def free_port() -> int:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(("127.0.0.1", 0))
        return int(s.getsockname()[1])


def _port_open(port: int) -> bool:
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.settimeout(0.5)
        return s.connect_ex(("127.0.0.1", port)) == 0


_TOKEN_LINE = re.compile(r"cat (\S+voila-token-\S+)")


def _read_token(log: Path, pid: int, env: dict) -> str:
    try:
        text = log.read_text(encoding="utf-8", errors="replace")
    except OSError:
        text = ""
    m = _TOKEN_LINE.search(text)
    candidates = [Path(m.group(1))] if m else []
    dirs = [env.get("XDG_RUNTIME_DIR") or "",
            str(Path(env.get("XDG_CACHE_HOME") or (Path.home() / ".cache"))
                / "delfin" / "tokens")]
    for d in dirs:
        if d:
            candidates += sorted(Path(d).glob(f"voila-token-{pid}-*"))
    for c in candidates:
        try:
            st = c.stat()
            if hasattr(os, "getuid") and st.st_uid != os.getuid():
                continue
            tok = c.read_text(encoding="utf-8").strip()
            if tok:
                return tok
        except OSError:
            continue
    return ""


def _where_note() -> tuple[Optional[Path], Optional[bytes]]:
    try:
        from delfin.agent import where
        path = where.record_path()
    except Exception:
        return None, None
    try:
        return path, path.read_bytes()
    except OSError:
        return path, None


def _restore_where_note(path: Optional[Path], before: Optional[bytes]) -> None:
    """A launcher from before preview awareness writes the note that says
    where *the* dashboard runs; put the running dashboard's back."""
    if path is None:
        return
    try:
        now = path.read_bytes()
    except OSError:
        now = None
    if now == before:
        return
    try:
        if before is None:
            path.unlink()
        else:
            from delfin.agent.state_paths import write_text_atomic
            write_text_atomic(path, before.decode("utf-8"))
    except Exception:
        pass


def _url_path(worktree: Path) -> str:
    try:
        from delfin.cli_voila import _dashboard_default_url
        return _dashboard_default_url(str(worktree / RUNTIME_DIR / NOTEBOOK),
                                      str(worktree))
    except Exception:
        return f"/voila/render/{RUNTIME_DIR}/{NOTEBOOK}"


def _served_url_path(port: int, token: str) -> str:
    """The address ``/`` redirects to (the server's default_url), or "".

    Asked rather than computed: where the launcher puts the notebook
    depends on what it finds (a notebook on the search path is staged
    into ``delfin_voila_runtime/``; the packaged one inside the root is
    served in place), and that differs between branches."""
    import http.client
    from urllib.parse import quote, urlsplit
    try:
        conn = http.client.HTTPConnection("127.0.0.1", int(port), timeout=10)
        conn.request("GET", "/?token=" + quote(token, safe=""))
        resp = conn.getresponse()
        location = resp.getheader("Location") or ""
        conn.close()
    except (OSError, http.client.HTTPException):
        return ""
    if resp.status not in (301, 302, 303, 307, 308) or not location:
        return ""
    path = urlsplit(location).path
    if not path.startswith("/") or path.startswith("/login"):
        return ""
    return path


def start(repo, ref: str, *, python: Optional[str] = None,
          base_env: Optional[dict] = None, shared_checkout: str = "",
          wait_s: float = 180.0, progress=None) -> tuple[bool, str, Optional[Preview]]:
    """Open a preview of *ref*. ``(ok, message, preview)``."""
    say = progress or (lambda _m: None)
    python = python or sys.executable
    base_env = dict(os.environ if base_env is None else base_env)
    target, why = resolve_target(repo, ref)
    if target is None:
        return False, why, None
    s = slug(target.ref)
    running = load(s)
    if running is not None and _alive(running):
        return True, f"Already running: {target.ref} on port {running.port}.", running
    say(f"Preparing the worktree for {target.ref} ...")
    ok, msg, path, attached = ensure_worktree(repo, target)
    if not ok:
        return False, msg, None
    say(msg)
    head = BL._out(path, "rev-parse", "HEAD") or target.sha
    target = Target(ref=target.ref, sha=head, branch=target.branch, kind=target.kind)
    env = preview_env(base_env, path, target, shared_checkout or str(main_worktree(repo)))
    ok, delfin_file, msg = verify_import(python, env, path)
    if not ok:
        return False, msg, None
    say(msg)
    from delfin.agent.state_paths import ensure_dir, write_text
    ensure_dir(_registry_dir())
    log = _registry_dir() / f"{s}.log"
    write_text(log, "")
    note_path, note_before = _where_note()
    proc = None
    port = 0
    for _attempt in range(3):
        port = free_port()
        with open(log, "a", encoding="utf-8") as fh:
            proc = subprocess.Popen(
                [python, "-m", "delfin.cli_voila", "--no-browser", "--port",
                 str(port), "--strict-port"],
                cwd=str(path), env=env, stdin=subprocess.DEVNULL, stdout=fh,
                stderr=subprocess.STDOUT, start_new_session=True)
        _children[proc.pid] = proc
        try:
            from delfin.agent import lifeline
            lifeline.record_child(proc.pid, "branch-preview")
        except Exception:
            pass
        say(f"Starting the preview server on port {port} ...")
        deadline = time.monotonic() + wait_s
        while time.monotonic() < deadline:
            if proc.poll() is not None or _port_open(port):
                break
            time.sleep(0.25)
        if proc.poll() is None and _port_open(port):
            break
        text = log.read_text(encoding="utf-8", errors="replace")
        if proc.poll() is not None and "occupied" in text:
            continue
        _kill_group(proc.pid)
        _restore_where_note(note_path, note_before)
        return False, ("The preview server did not start:\n"
                       + "\n".join(text.strip().splitlines()[-15:])), None
    else:
        _restore_where_note(note_path, note_before)
        return False, "No free port could be held for the preview.", None
    token = ""
    for _ in range(40):
        token = _read_token(log, proc.pid, env)
        if token:
            break
        time.sleep(0.25)
    _restore_where_note(note_path, note_before)
    if not token:
        _kill_group(proc.pid)
        return False, f"The preview started but its token was not found (log {log}).", None
    from delfin.agent import proc_identity
    preview = Preview(
        slug=s, ref=target.ref, sha=target.sha, branch=target.branch,
        repo=str(main_worktree(repo)), worktree=str(path), port=port,
        token=token,
        url_path=_served_url_path(port, token) or _url_path(Path(path)),
        pid=proc.pid,
        proc_start=proc_identity.process_start(proc.pid),
        ticks=proc_identity.start_ticks(proc.pid),
        host=socket.gethostname(), started_at=time.time(),
        delfin_file=delfin_file, log=str(log), attached=attached)
    _save(preview)
    return True, (f"Preview of {target.ref} @ {target.sha[:8]} on port {port}; "
                  f"{msg}"), preview


def _group_alive(pgid: int) -> bool:
    try:
        os.killpg(pgid, 0)
        return True
    except ProcessLookupError:
        return False
    except PermissionError:
        return True
    except OSError:
        return False


#: The launchers this process started, so their exit is reaped (a zombie
#: still answers a signal and would read as running).
_children: dict = {}


def _reap(pid: int) -> None:
    proc = _children.get(pid)
    if proc is not None:
        try:
            proc.poll()
        except Exception:
            pass


def _kill_group(pgid: int, grace_s: float = 20.0) -> None:
    """SIGTERM the process group *pgid*, SIGKILL after *grace_s*.

    Called only with the group of a launcher this module started with
    ``start_new_session=True`` (its pid is the group id), after its
    identity was checked by pid and start time, or for a group whose
    leader is gone but whose members (the server) are left."""
    if pgid <= 1 or pgid == os.getpgid(0):
        return
    try:
        os.killpg(pgid, signal.SIGTERM)
    except OSError:
        return
    deadline = time.monotonic() + grace_s
    while time.monotonic() < deadline:
        _reap(pgid)
        if not _group_alive(pgid):
            break
        time.sleep(0.2)
    if _group_alive(pgid):
        try:
            os.killpg(pgid, signal.SIGKILL)
        except OSError:
            pass
    for _ in range(25):
        _reap(pgid)
        if not _group_alive(pgid):
            break
        time.sleep(0.2)


#: Upper bound for one step of a close (checking the server, the worktree
#: use, ``git status``, ``git worktree remove``). Measured need on a local
#: disk is well under a second; a home on a parallel filesystem can take
#: seconds, and a hung mount forever.
STEP_TIMEOUT_S = 30.0


@dataclass
class CloseResult:
    """The final state of a close. ``text`` is for the person (no port,
    token, pid or path); ``details`` holds those and, when a step failed,
    the commands that finish the job by hand."""
    ok: bool
    state: str          # closed | kept | ask | failed
    text: str
    details: str = ""
    ask: bool = False


def _bounded(fn, timeout: float, what: str):
    """``fn()`` in a daemon thread; TimeoutError after *timeout* seconds.

    The thread is left behind on a timeout (a hung syscall cannot be
    interrupted from here); the caller reports and moves on."""
    import threading
    box: dict = {}

    def _run():
        try:
            box["value"] = fn()
        except BaseException as exc:      # handed to the caller below
            box["error"] = exc
    t = threading.Thread(target=_run, daemon=True, name=f"preview-close:{what}")
    t.start()
    t.join(timeout)
    if t.is_alive():
        raise TimeoutError(f"{what} did not finish within {int(timeout)} s")
    if "error" in box:
        raise box["error"]
    return box.get("value")


def manual_commands(p: Preview, *, server: bool = True,
                    worktree: bool = True) -> list[str]:
    """The shell commands that do what a failed close did not."""
    import shlex
    q = shlex.quote
    cmds = []
    if server and p.pid > 1:
        cmds.append(f"kill -TERM -- -{int(p.pid)}    # the preview server's process group")
    if worktree:
        cmds.append(f"git -C {q(str(p.repo))} worktree remove {q(str(p.worktree))}")
    cmds.append(f"rm -f {q(str(_record_path(p.slug)))} {q(str(_pending_path(p.slug)))}")
    return cmds


def _details(p: Preview, *lines: str, server=False, worktree=True) -> str:
    body = [x for x in lines if x]
    body.append(f"worktree: {p.worktree}; port {p.port}; process group {p.pid}")
    body.append("To finish by hand:")
    body.extend("  " + c for c in manual_commands(p, server=server, worktree=worktree))
    return "\n".join(body)


def close(s: str, *, worktree: str = "if-clean", progress=None,
          step_timeout: float = STEP_TIMEOUT_S) -> CloseResult:
    """Close preview *s*: stop its server, then release its worktree.

    ``worktree``: ``if-clean`` (remove when clean, else ask), ``keep``, or
    ``discard`` (remove after the person confirmed; the head commit is
    logged). Every step is bounded by *step_timeout* and every failure ends
    in ``state="failed"`` (or ``kept``) with the commands that finish the
    job; nothing raises. A record whose server is gone skips the stop."""
    say = progress or (lambda _m: None)
    try:
        p = load(s)
    except Exception as exc:
        return CloseResult(False, "failed", f"Couldn't read the preview's record: {exc}.",
                           f"record: {_record_path(s)}")
    if p is None:
        return CloseResult(False, "failed", "This preview is no longer listed.",
                           f"record: {_record_path(s)}")
    if p.host != socket.gethostname():
        return CloseResult(False, "failed", (
            f"This preview runs on {p.host}; close it from the dashboard there."),
            _details(p, server=True, worktree=True))
    say("Stopping the preview server …")
    try:
        up = _bounded(lambda: _alive(p) or _group_alive(p.pid), step_timeout,
                      "checking the server")
        if up:
            grace = max(1.0, min(20.0, step_timeout / 2))
            _bounded(lambda: _kill_group(p.pid, grace_s=grace),
                     step_timeout + grace + 6, "stopping the server")
            if _group_alive(p.pid):
                raise RuntimeError("the server did not stop")
    except Exception as exc:
        return CloseResult(False, "failed", f"Couldn't stop the preview: {exc}.",
                           _details(p, str(exc), server=True))
    try:
        from delfin.agent import lifeline
        _bounded(lambda: lifeline.end_children((p.pid, p.ticks)), step_timeout,
                 "ending the preview's agents")
    except Exception:
        pass
    try:
        _record_path(s).unlink()
    except FileNotFoundError:
        pass
    except OSError as exc:
        return CloseResult(False, "failed", f"Stopped, but the record could not be removed: {exc}.",
                           _details(p, str(exc), worktree=False))
    return _release(p, worktree, say, step_timeout)


def _release(p: Preview, mode: str, say, step_timeout: float) -> CloseResult:
    path = Path(p.worktree)
    repo = p.repo
    try:
        exists = _bounded(path.exists, step_timeout, "looking at the worktree")
    except Exception as exc:
        return CloseResult(False, "failed", f"Stopped. Couldn't check its worktree: {exc}.",
                           _details(p, str(exc)))
    if mode == "keep" or not exists:
        _drop_pending(p.slug)
        return CloseResult(True, "closed", "Closed." + (
            " Its worktree is kept." if exists else ""),
            f"worktree: {path}" if exists else "")
    if not is_preview_worktree(repo, path):
        return CloseResult(True, "kept", "Closed. Its folder is not a preview worktree, so it stays.",
                           f"worktree: {path}")
    say("Checking the worktree …")
    try:
        from delfin.dashboard.worktree_use import in_use
        busy = _bounded(lambda: in_use([str(path)]).get(str(path)), step_timeout,
                        "checking which programs use the worktree")
    except Exception as exc:
        return CloseResult(True, "kept", (
            f"Closed. Its worktree is kept: checking it failed ({exc})."),
            _details(p, str(exc)))
    if busy:
        return CloseResult(True, "kept", (
            "Closed. Its worktree is kept because a program is still using it; "
            "close the preview again once that has ended."), _details(p, busy))
    try:
        dirty = _bounded(lambda: changes(path), step_timeout, "git status")
    except Exception as exc:
        return CloseResult(True, "kept", (
            f"Closed. Its worktree is kept: checking it for changes failed ({exc})."),
            _details(p, str(exc)))
    if dirty and mode != "discard":
        try:
            _save_pending(p)
        except Exception:
            pass
        return CloseResult(True, "ask", (
            "Closed. Its worktree has changes — remove it anyway, or keep it?"),
            f"{path}:\n{dirty}", ask=True)
    say("Removing the worktree …")
    try:
        head = _bounded(lambda: BL._out(path, "rev-parse", "HEAD"), step_timeout,
                        "git rev-parse")
        branch = _bounded(lambda: BL._out(path, "branch", "--show-current"),
                          step_timeout, "git branch")
        r = _bounded(lambda: BL._git(repo, "worktree", "remove",
                                     *(["--force"] if dirty else []), str(path),
                                     timeout=step_timeout),
                     step_timeout + 5, "git worktree remove")
    except Exception as exc:
        return CloseResult(False, "failed", f"Stopped, but couldn't remove its worktree: {exc}.",
                           _details(p, str(exc)))
    if r.returncode != 0:
        why = (r.stderr or r.stdout).strip() or "git worktree remove failed"
        return CloseResult(False, "failed", "Stopped, but couldn't remove its worktree.",
                           _details(p, why))
    _drop_pending(p.slug)
    from delfin.dashboard import branch_cleanup as BC
    import shlex
    undo = (f"git -C {shlex.quote(repo)} worktree add {shlex.quote(str(path))} "
            + (shlex.quote(branch) if branch else f"--detach {head}"))
    BC._log({"at": time.time(), "repo": repo, "kind": "worktree",
             "name": str(path), "path": str(path), "branch": branch,
             "sha": head, "reason": "preview closed"
             + (" (changes discarded)" if dirty else ""), "undo": undo})
    return CloseResult(True, "closed", "Closed." + (
        " Its changes were discarded." if dirty else ""),
        f"Removed {path}. Undo: {undo}")


def _release_worktree(p: Preview, msg: str, mode: str) -> tuple[bool, str, bool]:
    """:func:`_release` as ``(ok, message, needs_confirm)``."""
    r = _release(p, mode, lambda _m: None, STEP_TIMEOUT_S)
    text = (msg + " " if msg else "") + r.text
    return r.ok, text + ("\n" + r.details if r.details else ""), r.ask


def stop(s: str, *, worktree: str = "if-clean") -> tuple[bool, str, bool]:
    """:func:`close` as ``(ok, message, needs_confirm)``; the message
    carries the details too (the command line prints it)."""
    r = close(s, worktree=worktree)
    msg = r.text + ("\n" + r.details if r.details else "")
    return r.ok, msg, r.ask


def _pending_path(s: str) -> Path:
    return _registry_dir() / f"{s}.closing.json"


def _save_pending(p: Preview) -> None:
    from delfin.agent.state_paths import write_text_atomic
    write_text_atomic(_pending_path(p.slug), json.dumps(asdict(p)))


def _drop_pending(s: str) -> None:
    try:
        _pending_path(s).unlink()
    except OSError:
        pass


def finish_close(s: str, *, discard: bool,
                 step_timeout: float = STEP_TIMEOUT_S) -> tuple[bool, str]:
    """The answer to the question :func:`close` asked."""
    r = finish_close_result(s, discard=discard, step_timeout=step_timeout)
    return r.ok, r.text + ("\n" + r.details if r.details else "")


def finish_close_result(s: str, *, discard: bool,
                        step_timeout: float = STEP_TIMEOUT_S) -> CloseResult:
    try:
        p = Preview(**json.loads(_pending_path(s).read_text(encoding="utf-8")))
    except (OSError, ValueError, TypeError):
        return CloseResult(False, "failed", "Nothing is waiting for that answer.")
    if not discard:
        _drop_pending(s)
        return CloseResult(True, "closed", "Closed. Its worktree is kept.",
                           f"worktree: {p.worktree}")
    return _release(p, "discard", lambda _m: None, step_timeout)


def pending_closes(repo=None) -> list[Preview]:
    """Closed previews whose worktree question is still open."""
    out: list[Preview] = []
    try:
        files = sorted(_registry_dir().glob("*.closing.json"))
    except OSError:
        return out
    key = BL.repo_key(repo) if repo else ""
    for f in files:
        try:
            p = Preview(**json.loads(f.read_text(encoding="utf-8")))
        except (OSError, ValueError, TypeError):
            continue
        if key and BL.repo_key(p.repo) != key:
            continue
        out.append(p)
    return out


def in_preview(env=None) -> bool:
    env = os.environ if env is None else env
    return bool(env.get(PREVIEW_ENV_BRANCH) and env.get(PREVIEW_ENV_WORKTREE))


def banner_info(env=None) -> dict:
    """What the preview banner says, computed in the running kernel."""
    env = os.environ if env is None else env
    import delfin
    file = str(Path(delfin.__file__).resolve())
    code = Path(file).parent.parent
    wt = Path(env.get(PREVIEW_ENV_WORKTREE) or "")
    try:
        inside = wt.resolve() in Path(file).parents
    except OSError:
        inside = False
    sha = BL._out(code, "rev-parse", "--short", "HEAD") or (
        env.get(PREVIEW_ENV_SHA) or "")[:8]
    return {"branch": env.get(PREVIEW_ENV_BRANCH) or "", "sha": sha,
            "code": str(code), "delfin_file": file, "worktree": str(wt),
            "inside": inside}


def banner_html(env=None) -> str:
    import html as _html
    i = banner_info(env)
    ok = i["inside"]
    bg = "#e65100" if ok else "#b71c1c"
    warn = "" if ok else (" &mdash; WARNING: this code is NOT from the preview "
                          f"worktree {_html.escape(i['worktree'])}")
    return (
        f'<div class="delfin-preview-banner" data-delfin-file="{_html.escape(i["delfin_file"])}" '
        f'style="background:{bg}; color:#fff; padding:6px 12px; '
        'border-radius:4px; font-family:monospace; font-weight:600; margin:4px 0;">'
        f'PREVIEW: {_html.escape(i["branch"])} @ {_html.escape(i["sha"])} &mdash; '
        f'code from {_html.escape(i["code"])}{warn}'
        f'<div style="font-weight:400; font-size:11px;">delfin.__file__ = '
        f'{_html.escape(i["delfin_file"])} &middot; the shared checkout is not '
        'touched; branch tools are in the main dashboard</div></div>')


def _main(argv=None) -> int:
    """``python -m delfin.dashboard.branch_preview list|start <ref>|stop <slug>``."""
    import argparse
    ap = argparse.ArgumentParser(prog="python -m delfin.dashboard.branch_preview")
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    st = sub.add_parser("start")
    st.add_argument("ref")
    st.add_argument("--repo", default="")
    sp = sub.add_parser("stop")
    sp.add_argument("slug")
    sp.add_argument("--keep-worktree", action="store_true")
    args = ap.parse_args(argv)
    if args.cmd == "list":
        for p in list_previews():
            print(f"{p.slug}\t{p.ref}@{p.sha[:8]}\tport {p.port}\t{p.worktree}")
        return 0
    if args.cmd == "start":
        import delfin
        repo = args.repo or os.environ.get(SHARED_CHECKOUT_ENV) or str(
            Path(delfin.__file__).resolve().parent.parent)
        ok, msg, p = start(repo, args.ref, progress=print)
        print(msg)
        if ok and p:
            print(f"Open: {p.url()}")
        return 0 if ok else 1
    ok, msg, _ask = stop(args.slug, worktree="keep" if args.keep_worktree else "if-clean")
    print(msg)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(_main())
