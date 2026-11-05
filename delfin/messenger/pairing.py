"""Pairing a phone without admin rights: the person's side and the file channel.

A person starts pairing in their own session (``/pair``) or terminal
(``delfin-agent messenger pair``). The request is a file in the person's
``~/.delfin/messenger_pairing`` (0700 directory, 0600 files, masked from
agent file tools and shells like the activation records). The bridge, the
agent-service user's door for the gateway, collects the requests of every
person (:func:`collect`) and writes the gateway's answers back
(:func:`answer`); the gateway never reads a person's home and the person
never reaches the gateway.

Flow (each arrow one poll of the gateway, about 5 s):

1. ``pair``: the request carries ``sha256(nonce)``; the nonce stays in the
   memory of the process that asked. The gateway answers ``offered`` with a
   one-time SimpleX invitation (128-bit+ random link, 10 minutes, single
   use, bound to the person and this request).
2. The phone connects: the gateway binds the contact to the person,
   UNVERIFIED, sends the phone a 6-digit code and answers ``connected``.
3. ``verify``: the person types the code where they typed ``/pair``; the
   request carries the nonce, so only the asking process can confirm. The
   bridge checks ``sha256(nonce)`` against the pair request; the gateway
   checks the code (3 attempts, then the contact is removed).

Until step 3 the contact is told only to verify first. The link and the
SimpleX security code are never written to a file. The asking process
opens a pairing door (:class:`Door`): a 0600 Unix socket in its own private
runtime directory (0700; the agent service's socket lives there too), named
in the pair request. The bridge hands those two fields to that door only --
after checking that it lies in the person's own runtime directory and that
the listening process runs as this uid -- and the door accepts them only
from a peer of this uid that is not an agent process, only for its own
request id. They stay in the asking process's memory (:func:`read`); the
answer file carries the state alone, a log line the 12-character tag.

Why a socket and not the agent service's protocol: the asking process is
not always a service session (the terminal REPL, the dashboard kernel,
``delfin-agent messenger pair``), and the bridge must reach it directly
rather than through a durable session event log.
"""
from __future__ import annotations

import hashlib
import hmac
import json
import os
import re
import secrets
import socket
import stat
import threading
import time
from pathlib import Path

#: How long an invitation link is valid.
INVITATION_TTL_S = 600
#: How long a connected phone has to be confirmed.
VERIFY_TTL_S = 600
#: How long a request file is kept unanswered or after its final answer.
REQUEST_TTL_S = 900
#: Paired phones per person.
MAX_PHONES = 3
#: Invitations per person per window.
INVITES_PER_WINDOW = 3
INVITE_WINDOW_S = 600
#: The confirmation code the phone receives.
CODE_DIGITS = 6
MAX_ATTEMPTS = 3

KINDS = frozenset({"pair", "verify", "phones", "revoke", "cancel"})
FINAL_STATES = frozenset({"verified", "refused", "expired", "failed", "done", "cancelled"})
MAX_FILE_BYTES = 16384

_DIR = Path.home() / ".delfin" / "messenger_pairing"
_ID = re.compile(r"[0-9a-f]{32}")
_HEX = re.compile(r"[0-9a-f]{64}")

#: requester -> (pair request id, nonce, door): kept only in this process.
_OPEN: dict[str, tuple] = {}

#: The answer fields handed to the asking process's door, never written.
DOOR_FIELDS = ("link", "securityCode")
_OPEN_LOCK = threading.Lock()


def directory() -> Path:
    return _DIR


def _sha(text: str) -> str:
    return hashlib.sha256(text.encode()).hexdigest()


def digits(text: object) -> str:
    return re.sub(r"\D", "", str(text or ""))


def new_code() -> str:
    """The code sent to a connected phone: CODE_DIGITS random digits."""
    return "".join(secrets.choice("0123456789") for _ in range(CODE_DIGITS))


def mask(text: object) -> str:
    """``text`` without SimpleX links and long digit runs: for log lines."""
    text = str(text or "")
    text = re.sub(r"(simplex:/|https?://(?:www\.)?simplex\.chat/)\S*", r"\1…", text)
    text = re.sub(r"\d(?:[\d ]{4,})\d", "…", text)
    return text


def _door_dir(directory: Path | None = None) -> Path:
    """Where this process's pairing door goes: its private runtime
    directory (the one its agent service socket is in)."""
    from delfin.agent.service import rundir
    return rundir.ensure_runtime_dir()


def _peer_is_agent(sock) -> bool:
    """Whether the peer of ``sock`` may be an agent's process (the agent
    service's own check: a descendant of this process's agents, the agent
    marker in its environment, or an environment that cannot be read).
    Never this process itself, which is no agent."""
    from types import SimpleNamespace
    from delfin.agent.service.supervisor import Supervisor, _peer_pid
    pid = _peer_pid(sock)
    if pid == os.getpid():
        return False
    identity = SimpleNamespace(sessions={}, starting=[SimpleNamespace(proc=SimpleNamespace(pid=os.getpid()))])
    return Supervisor._from_an_agent(identity, pid)


def _start_time(pid) -> str | None:
    """A process's start time (``/proc/<pid>/stat`` field 22, clock ticks
    since boot) as text, or None where it cannot be read. With the pid it
    names one process: a pid reused later has another start time."""
    try:
        with open(f"/proc/{int(pid)}/stat", "rb") as stream:
            data = stream.read(4096).decode("ascii", "replace")
        return data[data.rindex(")") + 2:].split()[19]
    except (OSError, ValueError, IndexError):
        return None


def _peer_pid(sock) -> int:
    from delfin.agent.service.supervisor import _peer_pid as peer
    try:
        return int(peer(sock))
    except Exception:  # noqa: BLE001 - unknown peer: no pid
        return 0


def _peer_uid(sock):
    from delfin.agent.service import rundir
    return rundir.peer_uid(sock)


class Door:
    """The asking process's end of the hand-over: a 0600 Unix socket at
    ``<runtime dir>/pair-<random>.sock`` that takes one JSON line per
    connection -- ``{"id": <request id>, "link"?, "securityCode"?}`` -- from
    a peer of this uid that is not an agent, for this request only, and
    keeps the fields in memory. Answers ``{"ok": true|false}``."""

    def __init__(self, rid: str, where: Path):
        self.rid = rid
        self.values: dict = {}
        self.path = str(Path(where) / f"pair-{secrets.token_hex(8)}.sock")
        self._sock = socket.socket(socket.AF_UNIX, socket.SOCK_STREAM)
        old = os.umask(0o177)
        try:
            self._sock.bind(self.path)
        finally:
            os.umask(old)
        os.chmod(self.path, 0o600)
        self._sock.listen(4)
        self._sock.settimeout(0.5)
        self._closed = False
        self._thread = threading.Thread(target=self._serve, name="messenger-pair-door", daemon=True)
        self._thread.start()

    def _serve(self) -> None:
        while not self._closed:
            try:
                conn, _ = self._sock.accept()
            except socket.timeout:
                continue
            except OSError:
                return
            with conn:
                ok = False
                try:
                    conn.settimeout(5)
                    if _peer_uid(conn) == os.getuid() and not _peer_is_agent(conn):
                        line = conn.makefile("rb").readline(MAX_FILE_BYTES + 1)
                        data = json.loads(line) if len(line) <= MAX_FILE_BYTES else None
                        if isinstance(data, dict) and hmac.compare_digest(str(data.get("id") or ""), self.rid):
                            self.values.update({k: str(v) for k, v in data.items() if k in DOOR_FIELDS})
                            ok = True
                    conn.sendall(b'{"ok": true}\n' if ok else b'{"ok": false}\n')
                except (OSError, ValueError):
                    pass

    def close(self) -> None:
        self._closed = True
        try:
            self._sock.close()
        finally:
            try:
                os.unlink(self.path)
            except OSError:
                pass


def deliver(door: str, root: Path | None, rid: str, fields: dict, *,
            pid=None, start=None) -> bool:
    """The bridge's side: hand ``fields`` to the door ``door`` of request
    ``rid``. Only a socket inside ``root`` (the person's own runtime
    directory), owned by this uid, whose listening process (SO_PEERCRED)
    runs as this uid, is the process the request names (``pid`` and, where
    /proc says, its ``start`` time), and is not an agent's process. True
    when the door took them."""
    if root is None or not door or not pid:
        return False
    try:
        real = Path(door).resolve()
        real.relative_to(Path(root).resolve())
        info = os.lstat(door)
    except (OSError, ValueError):
        return False
    if not stat.S_ISSOCK(info.st_mode) or info.st_uid != os.getuid() or real.parent != Path(root).resolve():
        return False
    try:
        with socket.socket(socket.AF_UNIX, socket.SOCK_STREAM) as sock:
            sock.settimeout(5)
            sock.connect(str(real))
            if _peer_uid(sock) != os.getuid() or _peer_is_agent(sock):
                return False
            listener = _peer_pid(sock)
            if listener != int(pid):
                return False
            if start is not None and _start_time(listener) != str(start):
                return False
            sock.sendall((json.dumps({"id": rid, **fields}) + "\n").encode())
            reply = sock.makefile("rb").readline(100)
    except OSError:
        return False
    try:
        return json.loads(reply).get("ok") is True
    except (ValueError, AttributeError):
        return False


def opened(res: dict | None, requester: str) -> dict | None:
    """An answer with what the door received, when ``requester`` is the
    process that asked for this pairing; otherwise the answer alone."""
    if res is None:
        return None
    out = dict(res)
    with _OPEN_LOCK:
        entry = _OPEN.get(str(requester))
    if entry is not None and entry[2] is not None:
        out.update(entry[2].values)
    return out


def read(rid: str, requester: str, *, directory: Path | None = None) -> dict | None:
    """The answer to request ``rid`` as ``requester`` sees it."""
    res = result(rid, directory=directory)
    with _OPEN_LOCK:
        entry = _OPEN.get(str(requester))
    if entry is None or entry[0] != rid:
        return res
    return opened(res, requester)


def link_tag(link: str) -> str:
    """The 12-character hash an audit line may carry instead of the link."""
    return _sha(str(link))[:12]


# -- files ---------------------------------------------------------------------------

def _ensure(path: Path) -> Path:
    path.mkdir(parents=True, exist_ok=True, mode=0o700)
    try:
        if stat.S_IMODE(path.stat().st_mode) & 0o077:
            os.chmod(path, 0o700)
    except OSError:
        pass
    return path


def _write(path: Path, data: dict) -> None:
    _ensure(path.parent)
    tmp = path.with_name("." + path.name + "." + secrets.token_hex(4))
    fd = os.open(tmp, os.O_WRONLY | os.O_CREAT | os.O_EXCL | getattr(os, "O_NOFOLLOW", 0), 0o600)
    try:
        with os.fdopen(fd, "w") as stream:
            json.dump(data, stream)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(tmp, path)
    finally:
        if os.path.lexists(tmp):
            os.unlink(tmp)


def _read(path: Path) -> dict | None:
    """A private file of this uid, or None: a link, another owner, a mode
    others can read, a size over MAX_FILE_BYTES or no JSON object refuse."""
    try:
        fd = os.open(path, os.O_RDONLY | getattr(os, "O_NOFOLLOW", 0) | os.O_NONBLOCK)
    except OSError:
        return None
    try:
        info = os.fstat(fd)
        if (not stat.S_ISREG(info.st_mode) or info.st_uid != os.getuid()
                or stat.S_IMODE(info.st_mode) & 0o077 or info.st_size > MAX_FILE_BYTES):
            return None
        data = os.read(fd, MAX_FILE_BYTES + 1)
    finally:
        os.close(fd)
    try:
        value = json.loads(data)
    except ValueError:
        return None
    return value if isinstance(value, dict) else None


def submit(kind: str, fields: dict | None = None, *, directory: Path | None = None,
           rid: str | None = None) -> str:
    if kind not in KINDS:
        raise ValueError("unknown pairing request")
    rid = rid or secrets.token_hex(16)
    _write((directory or _DIR) / f"req-{rid}.json",
           {"id": rid, "kind": kind, "created": time.time(), **(fields or {})})
    return rid


def result(rid: str, *, directory: Path | None = None) -> dict | None:
    if not _ID.fullmatch(str(rid or "")):
        return None
    return _read((directory or _DIR) / f"res-{rid}.json")


def wait(rid: str, done, *, timeout: float, directory: Path | None = None,
         interval: float = 0.25, sleep=time.sleep, clock=time.monotonic) -> dict | None:
    """The first result for which ``done(result)`` holds, or None at timeout."""
    deadline = clock() + timeout
    while True:
        value = result(rid, directory=directory)
        if value is not None and done(value):
            return value
        if clock() >= deadline:
            return None
        sleep(interval)


# -- the person's side ---------------------------------------------------------------

def begin(requester: str, *, directory: Path | None = None) -> str:
    """A pair request bound to ``requester`` (a session id or a process tag)."""
    nonce = secrets.token_hex(16)
    rid = secrets.token_hex(16)
    forget(requester)
    door = Door(rid, _door_dir(directory))
    submit("pair", {"requester": str(requester)[:80], "nonce_sha": _sha(nonce),
                    "door": door.path, "pid": os.getpid(), "start": _start_time(os.getpid())},
           directory=directory, rid=rid)
    with _OPEN_LOCK:
        _OPEN[str(requester)] = (rid, nonce, door)
    return rid


def open_request(requester: str) -> str | None:
    with _OPEN_LOCK:
        entry = _OPEN.get(str(requester))
    return entry[0] if entry else None


def confirm(requester: str, code: str, *, directory: Path | None = None) -> str | None:
    """The verify request for the pairing this requester started, or None
    when it started none (another session's pairing cannot be confirmed)."""
    with _OPEN_LOCK:
        entry = _OPEN.get(str(requester))
    if entry is None:
        return None
    rid, nonce = entry[0], entry[1]
    return submit("verify", {"pair": rid, "nonce": nonce, "code": digits(code)[:16]},
                  directory=directory)


def forget(requester: str) -> None:
    """Drop this process's pairing for ``requester`` and close its door."""
    with _OPEN_LOCK:
        entry = _OPEN.pop(str(requester), None)
    if entry is not None and entry[2] is not None:
        entry[2].close()


def qr_text(link: str) -> str | None:
    """The link as a text QR code (qrencode -t utf8, an argument, never a
    shell), or None when qrencode is not installed."""
    import shutil
    import subprocess
    tool = shutil.which("qrencode")
    if tool is None:
        return None
    try:
        done = subprocess.run([tool, "-t", "utf8", "-m", "1", link], capture_output=True,
                              text=True, timeout=10, check=False)
    except (OSError, subprocess.SubprocessError):
        return None
    return done.stdout if done.returncode == 0 and done.stdout.strip() else None


STEPS = ("Connect your phone in three steps:\n"
         "  1. /pair in your session (done)\n"
         "  2. scan the QR code (or paste the link) in SimpleX: + > Scan / paste link\n"
         "  3. type here the 6-digit code DELFIN sends to your phone: /pair 123456\n"
         "Then: hol / get NAME brings a running session to the phone, neu / new NAME "
         "starts one, /sync connects this one.")


def offer_text(res: dict) -> str:
    link = str(res.get("link") or "")
    minutes = max(1, int((float(res.get("expires") or time.time()) - time.time()) // 60))
    qr = qr_text(link)
    parts = [STEPS, ""]
    if qr:
        parts.append(qr.rstrip("\n"))
    else:
        parts.append("(a QR code needs qrencode: apt install qrencode)")
    parts += [f"Link (one-time, valid {minutes} min, only for you):", link]
    return "\n".join(parts)


def state_text(res: dict | None) -> str:
    """One answer of the gateway, in words for the person."""
    if res is None:
        return ("No answer from the messenger gateway. Is it running? "
                "(systemctl status delfin-messenger-gateway)")
    state, kind = res.get("state"), res.get("kind")
    if kind == "phones":
        rows = res.get("phones") or []
        lines = [f"  {r.get('label') or '?'} · {r.get('state')} · since "
                 f"{time.strftime('%Y-%m-%d %H:%M', time.localtime(float(r.get('since') or 0)))}"
                 for r in rows]
        return "Your paired phones:\n" + "\n".join(lines) if lines else "No phone paired."
    if kind == "revoke":
        gone = res.get("revoked") or []
        return ("Unpaired: " + ", ".join(gone) + ". They can no longer reach your sessions."
                if gone else "No phone matched; nothing unpaired.")
    if state == "connected":
        code = str(res.get("securityCode") or "")
        text = (f"Phone connected ({res.get('label') or 'SimpleX contact'}). It is NOT yet allowed "
                "to do anything.\nDELFIN sent a 6-digit code to the phone. Type it here: /pair CODE "
                "(terminal: at the prompt).")
        if code:
            text += ("\nOptionally compare the security code on the phone (the contact > Verify "
                     f"security code): {code}")
        return text
    if state == "verified":
        return (f"Phone verified ✓ ({res.get('label') or 'phone'}). It reaches your sessions only.\n"
                "On the phone: sessions, hol / get NAME, neu / new NAME, help.")
    if state == "wrong":
        return f"Wrong code. {int(res.get('left') or 0)} attempt(s) left: /pair CODE"
    if state == "failed":
        return "Too many wrong codes: the phone was removed. Start again with /pair."
    if state == "expired":
        return "The invitation or the confirmation time expired. Start again with /pair."
    if state == "cancelled":
        return "Pairing cancelled."
    if state == "refused":
        reason = {"limit": f"at most {MAX_PHONES} phones per person: /pair revoke LABEL first",
                  "rate": "too many invitations; wait ten minutes",
                  "store": "the gateway's pairing store is refused; ask the admin",
                  }.get(str(res.get("reason") or ""), str(res.get("reason") or "refused"))
        return "Pairing refused: " + reason + "."
    if state == "none":
        return "No connected phone waits for a code. Start with /pair."
    return "Pairing: " + str(state or "?")


def _spawn(fn) -> None:
    threading.Thread(target=fn, name="messenger-pair", daemon=True).start()


def session_command(args: str, *, requester: str, show, spawn=_spawn,
                    directory: Path | None = None, timeout: float = 30.0) -> str:
    """``/pair [CODE | status | phones | revoke LABEL|all | cancel]`` in a session.

    Returns the immediate answer; ``show(text)`` is called later from a
    background thread with what the gateway answered (the invitation, the
    connection, the verification). The invitation text is meant for a
    display that is not kept (a live notice), never for a transcript."""
    words = str(args or "").strip()
    head, _, rest = words.partition(" ")
    head = head.lower()
    if not requester:
        return "Start a session first, then /pair."
    if re.fullmatch(r"[\d\s-]+", words) and digits(words):
        rid = confirm(requester, words, directory=directory)
        if rid is None:
            return "This session started no pairing. Start with /pair."

        def after_verify():
            res = wait(rid, lambda r: r.get("state") in FINAL_STATES | {"wrong", "none"},
                       timeout=timeout, directory=directory)
            show(state_text(res))
            if res and res.get("state") in FINAL_STATES:
                forget(requester)
        spawn(after_verify)
        return "Code sent to the gateway; the answer follows in a few seconds."
    if head in {"phones", "status", "list"}:
        rid = submit("phones", {}, directory=directory)
        spawn(lambda: show(state_text(wait(rid, lambda r: True, timeout=timeout,
                                            directory=directory))))
        return "Asking the gateway for your phones ..."
    if head in {"revoke", "unpair"}:
        label = rest.strip()
        if not label:
            return "/pair revoke LABEL (from /pair phones) or /pair revoke all"
        rid = submit("revoke", {"label": label[:64]}, directory=directory)
        spawn(lambda: show(state_text(wait(rid, lambda r: True, timeout=timeout,
                                            directory=directory))))
        return "Unpairing ..."
    if head == "cancel":
        old = open_request(requester)
        if old is None:
            return "No pairing open in this session."
        submit("cancel", {"pair": old}, directory=directory)
        forget(requester)
        return "Pairing cancelled."
    if head and head not in {"new", "neu", "start"}:
        return ("/pair · /pair CODE · /pair phones · /pair revoke LABEL|all · /pair cancel")
    rid = begin(requester, directory=directory)

    def after_pair():
        offered = opened(wait(rid, lambda r: r.get("state") in FINAL_STATES | {"offered"},
                              timeout=timeout, directory=directory), requester)
        if offered is None or offered.get("state") != "offered" or not offered.get("link"):
            show(state_text(offered))
            forget(requester)
            return
        show(offer_text(offered))
        connected = opened(wait(rid, lambda r: r.get("state") != "offered",
                                timeout=INVITATION_TTL_S + 30, directory=directory), requester)
        show(state_text(connected))
        if connected is None or connected.get("state") in FINAL_STATES:
            forget(requester)
    spawn(after_pair)
    return "Asking the messenger gateway for a one-time invitation ..."


CLI_USAGE = """delfin-agent messenger pair      connect your phone (QR code, then the code it receives)
delfin-agent messenger phones    your paired phones
delfin-agent messenger unpair LABEL|all
                                 unpair one of your phones, or all
delfin-agent messenger people    (admin) every person's phones
delfin-agent messenger revoke NAME
                                 (admin) unpair person NAME's phones (owner: the account owner's)"""


def _from_an_agent() -> bool:
    """True inside an agent's process tree (approval_guard's marker):
    pairing is the person's act, never an agent's."""
    try:
        from delfin.agent.approval_guard import ENV_MARKER
    except ImportError:
        return False
    return bool(os.environ.get(ENV_MARKER))


def cli_main(argv: list[str], *, ask=input, out=print, timeout: float = 30.0) -> int:
    """``delfin-agent messenger ...``: the person's pairing in a terminal,
    and the admin's listing and revocation (which need the gateway's store:
    run as the gateway user, e.g. ``sudo bash pi_messenger.sh --people``)."""
    words = list(argv or [])
    head = words[0] if words else ""
    if head in {"people", "revoke"}:
        from . import __main__ as admin
        return admin.main(words)
    if head not in {"pair", "phones", "unpair"}:
        out(CLI_USAGE)
        return 2
    if _from_an_agent():
        out("Refused: an agent session cannot pair, list or unpair phones.")
        return 2
    if head == "phones":
        rid = submit("phones", {})
        res = wait(rid, lambda r: True, timeout=timeout)
        out(state_text(res))
        for item in notes()[-5:]:
            out(time.strftime("  %Y-%m-%d %H:%M ", time.localtime(float(item.get("ts") or 0)))
                + str(item.get("text") or ""))
        return 0 if res else 1
    if head == "unpair":
        label = " ".join(words[1:]).strip()
        if not label:
            out("delfin-agent messenger unpair LABEL|all (labels: delfin-agent messenger phones)")
            return 2
        res = wait(submit("revoke", {"label": label[:64]}), lambda r: True, timeout=timeout)
        out(state_text(res))
        return 0 if res and res.get("revoked") else 1
    requester = f"cli:{os.getpid()}:{secrets.token_hex(4)}"
    try:
        return _cli_pair(requester, ask=ask, out=out, timeout=timeout)
    finally:
        forget(requester)


def _cli_pair(requester: str, *, ask, out, timeout: float) -> int:
    rid = begin(requester)
    out("Asking the messenger gateway for a one-time invitation ...")
    res = opened(wait(rid, lambda r: r.get("state") in FINAL_STATES | {"offered"},
                      timeout=timeout), requester)
    if res is None or res.get("state") != "offered" or not res.get("link"):
        out(state_text(res))
        return 1
    out(offer_text(res))
    out("waiting for the phone to connect ...")
    res = opened(wait(rid, lambda r: r.get("state") != "offered",
                      timeout=INVITATION_TTL_S + 30), requester)
    out(state_text(res) if res is None or res.get("state") != "connected"
        else state_text(res).replace("Type it here: /pair CODE (terminal: at the prompt).",
                                     "Type it below."))
    if res is None or res.get("state") != "connected":
        return 1
    for _ in range(MAX_ATTEMPTS):
        try:
            code = ask("Code from the phone: ")
        except (EOFError, KeyboardInterrupt, OSError):
            out("Not confirmed; the phone stays without rights and is removed after 10 minutes.")
            return 1
        vid = confirm(requester, code)
        answer_ = wait(vid, lambda r: r.get("state") in FINAL_STATES | {"wrong", "none"},
                       timeout=timeout) if vid else None
        out(state_text(answer_))
        if answer_ and answer_.get("state") == "verified":
            return 0
        if not answer_ or answer_.get("state") != "wrong":
            return 1
    return 1


# -- the bridge's side ---------------------------------------------------------------

def collect(path: Path, *, now: float | None = None) -> list[dict]:
    """The open requests in one person's pairing directory, checked.

    A request older than REQUEST_TTL_S, a file that is not private to this
    uid, a malformed one and a verify whose nonce does not match its pair
    request are removed (a verify that failed this check is answered
    ``refused`` so the asking process is told). The answer files of
    finished requests are removed after REQUEST_TTL_S too."""
    now = time.time() if now is None else now
    try:
        names = sorted(os.listdir(path))
    except OSError:
        return []
    requests: dict[str, dict] = {}
    for name in names:
        full = path / name
        if name.startswith("res-") and name.endswith(".json"):
            try:
                if now - os.lstat(full).st_mtime > REQUEST_TTL_S:
                    full.unlink()
            except OSError:
                pass
            continue
        if not (name.startswith("req-") and name.endswith(".json")):
            continue
        rid = name[4:-5]
        data = _read(full)
        ok = (data is not None and _ID.fullmatch(rid) and data.get("id") == rid
              and data.get("kind") in KINDS
              and isinstance(data.get("created"), (int, float))
              and 0 <= now - float(data["created"]) <= REQUEST_TTL_S)
        if not ok:
            try:
                full.unlink()
            except OSError:
                pass
            continue
        requests[rid] = data
    out = []
    for rid, data in requests.items():
        kind = data["kind"]
        item = {"id": rid, "kind": kind}
        if kind == "pair":
            if not _HEX.fullmatch(str(data.get("nonce_sha") or "")) or not data.get("door"):
                continue
            if result(rid, directory=path) and result(rid, directory=path).get("state") in FINAL_STATES:
                continue
        elif kind == "verify":
            pair = requests.get(str(data.get("pair") or ""))
            nonce = str(data.get("nonce") or "")
            if pair is None or pair.get("kind") != "pair" or not hmac.compare_digest(
                    _sha(nonce), str(pair.get("nonce_sha") or "")):
                answer(path, rid, {"kind": "verify", "state": "refused",
                                   "reason": "not the session that started this pairing"})
                continue
            item.update(pair=pair["id"], code=digits(data.get("code"))[:16])
        elif kind == "revoke":
            item["label"] = str(data.get("label") or "")[:64]
        elif kind == "cancel":
            item["pair"] = str(data.get("pair") or "")
        out.append(item)
    return out


def answer(path: Path, rid: str, value: dict, *, door_root: Path | None = None) -> dict:
    """Write the gateway's answer for request ``rid``; a final answer of a
    pair request, and any answer of the other kinds, removes the request.
    The link and the security code are handed to the request's door (inside
    ``door_root``, the person's runtime directory) and never written; an
    invitation whose door does not take it is answered ``refused``.
    Returns {"ok": written, "delivered": True | False | None (nothing to hand)}."""
    if not _ID.fullmatch(str(rid or "")):
        return {"ok": False, "delivered": None}
    req = path / f"req-{rid}.json"
    data = _read(req)
    if data is None:
        return {"ok": False, "delivered": None}
    value = {k: v for k, v in dict(value).items() if k in {
        "kind", "state", "reason", "link", "expires", "label", "securityCode", "left",
        "phones", "revoked"}}
    value["kind"] = data["kind"]
    # The link and the security code go to the asking process's door only.
    private = {k: str(value.pop(k)) for k in DOOR_FIELDS if k in value}
    delivered = None
    if private:
        delivered = data["kind"] == "pair" and deliver(
            str(data.get("door") or ""), door_root, rid, private,
            pid=data.get("pid") if isinstance(data.get("pid"), int) else None,
            start=data.get("start"))
        if not delivered and "link" in private:
            value.update(state="refused", reason="the asking window is gone")
    _write(path / f"res-{rid}.json", value)
    final = data["kind"] != "pair" or value.get("state") in FINAL_STATES
    if final:
        try:
            req.unlink()
        except OSError:
            pass
        # A finished pairing drops its confirmations still waiting.
        if data["kind"] == "pair":
            for other in os.listdir(path):
                if other.startswith("req-") and other.endswith(".json"):
                    item = _read(path / other)
                    if item and item.get("kind") == "verify" and item.get("pair") == rid:
                        try:
                            (path / other).unlink()
                        except OSError:
                            pass
    return {"ok": True, "delivered": delivered}


def note(path: Path, text: str) -> None:
    """A line for the person (a new phone, a revocation), kept 0600, the
    newest 50: shown by ``/pair phones`` and ``delfin-agent messenger phones``."""
    _ensure(path)
    file = path / "notes.json"
    old = (_read(file) or {}).get("notes") or []
    _write(file, {"notes": [*old, {"ts": time.time(), "text": mask(text)[:300]}][-50:]})


def notes(path: Path | None = None) -> list[dict]:
    return list((_read((path or _DIR) / "notes.json") or {}).get("notes") or [])
