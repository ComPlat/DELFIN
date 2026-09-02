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
    start_dir = str(cwd) if cwd else (iso.roots[0] if iso.roots else "/")
    argv += ["--chdir", start_dir, "--die-with-parent", command]
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
    """Forget the bwrap probe. For tests, and for a doctor re-run."""
    global _BWRAP_OK
    _BWRAP_OK = None


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
    "bwrap_argv",
    "bwrap_functional",
    "parse_isolation",
    "refusal_reason",
    "reset_probe_cache",
    "uncontained_note",
]
