"""Filesystem containment for stdio MCP servers, decided at the launch.

A tool reached over MCP is the one hole in the agent's containment, and
not by oversight. ``_gate_side_effecting_mcp_tool`` says so about itself:
the arguments belong to the server's own schema, so a path check at the
call site is reading a field whose meaning it does not know. And the
process is one DELFIN did not build a command line for, so the bwrap wrap
that the shell tool gets (``_bash_isolation_argv``) never covers it. Under
bypass this was measured: a read that native bash could not perform came
back through an MCP shell tool, from a directory nothing in the session
had named.

The one place where the question IS answerable is the launch. A server
started inside a namespace holding only the roots it was declared with
cannot reach outside them whatever it is later asked, by the model or by
anything that has learned to ask it. That moves the decision from the call
(where the data is opaque) to the process boundary (where it is a list of
paths the user wrote).

Opt-in, per server, in the MCP config::

    "docs-index": {
      "command": "…", "args": [...],
      "roots":      ["/data/project"],     read-write inside the namespace
      "read_roots": ["/opt/corpus"],       read-only
      "isolation":  "off"                  optional: run uncontained anyway
    }

Declaring nothing keeps today's behaviour exactly — and says so, in the
banner and in ``/mcp``, because "outside every containment" is a fact
worth reading rather than discovering.

Two deliberate refusals rather than a best effort:

* a declared root that does not exist stops the launch. bwrap would fail
  on it anyway, and a typo that silently yields *less* access than the
  line says is the failure this module exists to prevent.
* declared roots with no usable bwrap stop the launch too. The shell tool
  degrades to plain bash there and records why, which is right for a
  fallback nobody asked for; here the user wrote the containment down,
  and running without it is not a lesser version of what they asked for.
"""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

# Trees a server needs in order to RUN, as opposed to trees it might be
# asked to read. All read-only: an isolated server may execute the
# interpreter it was installed with and may not rewrite it. /opt is here
# because scientific software is installed there as a matter of course.
_SYSTEM_ROOTS = (
    "/usr", "/bin", "/sbin", "/lib", "/lib32", "/lib64", "/libx32",
    "/etc", "/opt",
)

# Directory names that mean "this is an install prefix, one level up".
_BIN_DIRS = ("bin", "sbin", "libexec")


@dataclass(frozen=True)
class Isolation:
    """What one server was declared to be allowed to see."""

    write_roots: tuple[str, ...] = ()
    read_roots: tuple[str, ...] = ()
    missing: tuple[str, ...] = ()

    @property
    def roots(self) -> tuple[str, ...]:
        return self.write_roots + self.read_roots

    def describe(self) -> str:
        """The roots, for a listing line. Never raises, never empty."""
        parts = [f"{r} (rw)" for r in self.write_roots]
        parts += [f"{r} (ro)" for r in self.read_roots]
        return ", ".join(parts) if parts else "no roots"


def parse_isolation(cfg: dict) -> Isolation | None:
    """The declared containment for one config entry, or ``None``.

    ``None`` means the server runs as it always has. That is the default
    and it is not a silent one — every caller that shows a server also
    shows whether this returned anything.
    """
    if not isinstance(cfg, dict):
        return None
    if str(cfg.get("isolation", "") or "").strip().lower() == "off":
        # An explicit escape hatch, kept next to the roots rather than in
        # a global setting: the server that needs it is one server.
        return None
    write_roots, read_roots, missing = [], [], []
    for key, sink in (("roots", write_roots), ("read_roots", read_roots)):
        raw = cfg.get(key)
        if isinstance(raw, str):
            raw = [raw]
        for entry in raw or ():
            text = str(entry or "").strip()
            if not text:
                continue
            try:
                path = Path(os.path.expandvars(text)).expanduser().resolve()
            except (OSError, ValueError):
                missing.append(text)
                continue
            if path.exists():
                sink.append(str(path))
            else:
                missing.append(str(path))
    if not (write_roots or read_roots or missing):
        return None
    return Isolation(tuple(write_roots), tuple(read_roots), tuple(missing))


# DELFIN's own servers are configured by the same settings file the
# dashboard writes, so their roots are READ rather than invented. Each
# entry is (settings path, writable) and resolves to nothing when the
# directory is absent — a user who does not keep an office folder should
# not have a server refuse to start over it.
def delfin_roots(settings: dict | None = None,
                 workspace: Path | str | None = None,
                 home: Path | str | None = None,
                 interpreter: str = "") -> Isolation | None:
    """The roots the built-in servers work in, from the user's settings.

    What goes in is what DELFIN is configured to use: the calculations
    folder, the office folder, the agent workspace, the two state
    directories the tools server persists into, and the launch workspace.
    Read-only: the archive, and the runtime trees that hold xtb and ORCA —
    a server may run the tools it was installed with and may not rewrite
    them.

    What stays OUT is as much of the decision as what goes in:

    * ``~/.delfin`` as a whole, because ``credentials.json`` lives there.
      The two sub-directories the servers persist into are bound by name.
    * ``transfer.remote_path``. It is another account's tree, reached over
      SSH by design — the network is untouched by the namespace, so the
      configured route keeps working, and the local shortcut into someone
      else's files stops being available. That shortcut is what the field
      report caught being taken.
    * the archive as a WRITE root. Every other route in the product
      refuses to write there, and a namespace that allowed it would grant
      what the path policy denies. ``move_to_archive`` is refused under
      isolation; a user who needs it says so with an explicit ``roots``
      entry of their own, which overrides this derivation.
    """
    try:
        if settings is None:
            from delfin.user_settings import load_settings
            settings = load_settings() or {}
    except Exception:
        settings = {}
    home_path = Path(home) if home is not None else Path.home()
    paths = (settings.get("paths") or {}) if isinstance(settings, dict) else {}
    runtime = (settings.get("runtime") or {}) if isinstance(settings, dict) else {}

    def configured(key: str, default: str) -> str:
        return str(paths.get(key) or "").strip() or str(home_path / default)

    write = [
        configured("calculations_dir", "calc"),
        configured("office_dir", "office"),
        str(home_path / "agent_workspace"),
        str(home_path / ".delfin" / "applications"),
        str(home_path / ".delfin" / "adapters"),
    ]
    if workspace:
        write.append(str(workspace))
    read = [configured("archive_dir", "archive")]
    qm_root = str(runtime.get("qm_tools_root") or "").strip()
    if qm_root:
        read.append(qm_root)
    orca = str(((runtime.get("local") or {}).get("orca_base") or "")).strip()
    if orca:
        read.append(orca)
    for binary in (runtime.get("tool_binaries") or {}).values():
        read.extend(_runtime_binds(str(binary or ""), home_path))

    def usable(entries) -> list[str]:
        out: list[str] = []
        for entry in entries:
            try:
                resolved = Path(entry).expanduser().resolve()
            except (OSError, ValueError):
                continue
            # An absent directory is dropped, not refused. Unlike a root
            # the user typed, this one was inferred — and an inference
            # that stops the server would be the derivation making policy.
            if resolved.is_dir() and str(resolved) not in out:
                out.append(str(resolved))
        return out

    write_roots = usable(write)
    read_roots = [r for r in usable(read) if r not in write_roots]
    if not write_roots and not read_roots:
        return None
    partial = Isolation(tuple(write_roots), tuple(read_roots))
    # DELFIN's own source, read-only, and ASKED FOR rather than assumed.
    # Measured: with only the data roots bound, the tools server died with
    # "No module named 'delfin'". Taking the answer from this process was
    # not enough either — on a checkout install the package resolves
    # through a .pth into a directory that has nothing to do with where
    # the launcher happens to be standing. The server's own interpreter,
    # starting where the server will start, is the only thing that knows.
    package = _package_root(interpreter or sys.executable,
                            _start_dir(partial, None))
    if package and package not in write_roots and package not in read_roots:
        read_roots.append(package)
    return Isolation(tuple(write_roots), tuple(read_roots))


_PACKAGE_PROBE = (
    "import os,delfin;print(os.path.dirname(os.path.dirname(delfin.__file__)))"
)
_PACKAGE_ROOTS: dict[tuple[str, str], str] = {}


def _package_root(interpreter: str, start_dir: str) -> str:
    """Where *interpreter*, started in *start_dir*, imports delfin from.

    One short subprocess per contained launch, and only when the derived
    containment is switched on. It fails to "" — an unanswerable question
    must not stop a server from starting; the missing bind then shows up
    as the server's own import error, which names the path.
    """
    if not interpreter:
        return ""
    key = (interpreter, start_dir)
    if key in _PACKAGE_ROOTS:
        # The listing, the banner and the doctor all build the same rows,
        # and each one used to pay for its own interpreter start.
        return _PACKAGE_ROOTS[key]
    try:
        done = subprocess.run([interpreter, "-c", _PACKAGE_PROBE],
                              capture_output=True, text=True, timeout=20,
                              cwd=start_dir or None)
        path = Path(done.stdout.strip()) if done.returncode == 0 else None
        found = str(path.resolve()) if path and path.is_dir() else ""
    except Exception:
        found = ""
    _PACKAGE_ROOTS[key] = found
    return found


def isolation_disabled(cfg: dict) -> bool:
    """True when the entry turned it off in so many words.

    ``parse_isolation`` returns None both for "nothing declared" and for
    "declared off", and those must not be treated alike: the first may be
    filled in from the settings, the second is the user saying no.
    """
    if not isinstance(cfg, dict):
        return False
    return str(cfg.get("isolation", "") or "").strip().lower() == "off"


def builtin_isolation_enabled(settings: dict | None = None) -> bool:
    """Is the derived containment switched on for DELFIN's own servers?

    ``agent.mcp_isolation``: "off" (the default) or "builtin". Opt-in on
    purpose — the roots are inferred, and an inference that is wrong takes
    out the dashboard's own engine. It graduates to a default once it has
    been run against a real session, not before.
    """
    try:
        if settings is None:
            from delfin.user_settings import load_settings
            settings = load_settings() or {}
        value = ((settings.get("agent") or {}).get("mcp_isolation", "off"))
    except Exception:
        return False
    return str(value or "off").strip().lower() == "builtin"


def _under(path: Path, root: str) -> bool:
    try:
        path.relative_to(root)
        return True
    except ValueError:
        return False


def _runtime_binds(command: str, home: Path) -> list[str]:
    """Read-only paths a server needs merely to start.

    The interpreter is routinely somewhere ``$HOME`` covers — a venv, a
    micromamba prefix, ``~/.local/bin`` — and this module empties ``$HOME``
    inside the namespace. Without this the containment would not be strict,
    it would simply mean nothing starts.

    An executable directly under ``$HOME`` (or above it) yields the file
    itself rather than its directory: binding the prefix there would undo
    the tmpfs it is meant to sit inside.
    """
    if not command:
        return []
    exe = command if os.path.isabs(command) else (shutil.which(command) or "")
    if not exe:
        return []
    try:
        resolved = Path(exe).resolve()
    except OSError:
        return []
    if not resolved.exists():
        return []
    if any(_under(resolved, root) for root in _SYSTEM_ROOTS):
        return []       # already covered by the system binds
    parent = resolved.parent
    prefix = parent.parent if parent.name in _BIN_DIRS else parent
    if prefix == home or _under(home, str(prefix)):
        return [str(resolved)]
    return [str(prefix)]


def _start_dir(iso: Isolation, cwd: Path | str | None) -> str:
    """Where the server starts, keeping the launcher's cwd when it can.

    Without isolation a server inherits the cwd of whatever started it,
    and plenty of them resolve work against it — `npx` finds its packages
    in `./node_modules`, an interpreter finds a package in `.`. Silently
    starting somewhere else would be a behaviour change smuggled in with
    a containment change, so the inherited directory is kept whenever it
    is inside the roots, and only otherwise does a root stand in for it.
    """
    if cwd:
        return str(cwd)
    try:
        here = Path(os.getcwd()).resolve()
        for root in iso.roots:
            if here == Path(root) or _under(here, root):
                return str(here)
    except OSError:                         # pragma: no cover - cwd unlinked
        pass
    return iso.roots[0] if iso.roots else "/"


def bwrap_argv(
    command: str,
    args,
    iso: Isolation,
    *,
    home: Path | str | None = None,
    cwd: Path | str | None = None,
) -> list[str]:
    """The argv that launches *command* seeing only what *iso* declares.

    Order is load-bearing: the system trees first, then an empty ``$HOME``
    over whatever the first step brought in, then the runtime paths and
    the declared roots — each of which may legitimately live under the
    home that was just blanked.
    """
    home_path = Path(home) if home is not None else Path.home()
    try:
        home_path = home_path.resolve()
    except OSError:            # pragma: no cover - resolve on a live $HOME
        pass
    argv: list[str] = ["bwrap", "--dev", "/dev", "--proc", "/proc",
                       "--tmpfs", "/tmp"]
    for root in _SYSTEM_ROOTS:
        if Path(root).exists():
            argv += ["--ro-bind", root, root]
    # A writable, empty home. Servers cache into it (npm, pip, uv) and a
    # read-only one turns a containment question into a startup failure
    # that reads like a broken install.
    argv += ["--tmpfs", str(home_path)]
    for path in _runtime_binds(command, home_path):
        argv += ["--ro-bind", path, path]
    for root in iso.read_roots:
        argv += ["--ro-bind", root, root]
    for root in iso.write_roots:
        argv += ["--bind", root, root]
    argv += ["--chdir", _start_dir(iso, cwd), "--die-with-parent", command]
    argv += [str(a) for a in (args or ())]
    return argv


_BWRAP_OK: bool | None = None


def bwrap_functional(probe: bool = True) -> bool:
    """True only if bwrap is installed and works in THIS namespace shape.

    Probed with the argv this module actually builds, not a minimal one:
    a stripped probe passes on hosts where the real wrap fails, and the
    consequence there would be a server refused at every start with the
    reason pointing at the wrong thing.
    """
    global _BWRAP_OK
    if _BWRAP_OK is not None:
        return _BWRAP_OK
    ok = False
    try:
        if shutil.which("bwrap"):
            true_bin = shutil.which("true") or "/bin/true"
            argv = bwrap_argv(true_bin, [], Isolation())
            if probe:
                ok = subprocess.run(
                    argv, capture_output=True, timeout=5).returncode == 0
            else:                                  # pragma: no cover
                ok = True
    except Exception:
        ok = False
    _BWRAP_OK = ok
    return ok


def reset_probe_cache() -> None:
    """Forget both probes. For tests, and for a doctor re-run."""
    global _BWRAP_OK
    _BWRAP_OK = None
    _PACKAGE_ROOTS.clear()


def refusal_reason(name: str, iso: Isolation) -> str:
    """Why a declared server was not started, or "" if it may start."""
    if iso.missing:
        return (f"isolation for '{name}' names a path that does not exist: "
                + ", ".join(iso.missing)
                + " — fix the root or set \"isolation\": \"off\"")
    if not bwrap_functional():
        return (f"'{name}' declares isolation ({iso.describe()}) but "
                "bubblewrap is not usable here, so the server was not "
                "started. Set \"isolation\": \"off\" to run it uncontained.")
    return ""


def uncontained_note(rows) -> str:
    """One line naming the servers no containment reaches, or "".

    Takes the ``effective_servers`` rows so the caller needs no registry
    and starts nothing. An HTTP server is counted too: it runs on a host
    DELFIN has no say over at all, which is further outside, not less.
    """
    names = []
    for row in rows or ():
        if not isinstance(row, dict) or not row.get("enabled", True):
            continue
        if row.get("isolation"):
            continue
        names.append(str(row.get("name", "?")))
    if not names:
        return ""
    shown = ", ".join(sorted(names)[:4])
    if len(names) > 4:
        shown += f", +{len(names) - 4} more"
    return (f"mcp        {shown} — tools reached through MCP run outside "
            "the shell's isolation (/mcp)")


__all__ = [
    "Isolation",
    "builtin_isolation_enabled",
    "bwrap_argv",
    "bwrap_functional",
    "delfin_roots",
    "isolation_disabled",
    "parse_isolation",
    "refusal_reason",
    "reset_probe_cache",
    "uncontained_note",
]
