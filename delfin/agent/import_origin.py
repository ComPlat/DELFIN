"""Which code a Python command in this workspace will REALLY import.

The false red of 2026-09-26: a session in a worktree ran
``python3 .gate/probe_login_shell_full.py`` and got a failure that was
not in its code. Python puts the SCRIPT's directory on ``sys.path[0]`` --
``.gate/``, never the worktree root -- so ``import delfin`` resolved
through the editable install (``__editable__.delfin_complat-*.pth`` in
the venv) into the MAIN checkout: the probe measured the unfixed code.
Nothing in the output says so; the fix was in the worktree all along.

``python -c`` and ``python -m`` run from the root do not have this hole
(they put the CWD on the path), and neither does the test gate (it sets
``PYTHONPATH`` to the worktree). A plain script in a subdirectory does.

This module is pure: it detects the shadow and words the note. The note
is attached to bash results in api_client, next to ``_pipe_exit_note``.
"""
import re
import subprocess
import sys
from pathlib import Path

__all__ = ["shadowed_packages", "note_for_command"]

#: How the probe resolves each package: it prints ``name<TAB>origin`` for
#: every workspace source package it can find a spec for. The probe runs
#: with ``-I`` (isolated: no PYTHONPATH, no script dir, no cwd injection)
#: but WITH site-packages, and appends any *extra_paths* -- staging an
#: "installed" twin in tests, and matching what a subdirectory script
#: sees once its own directory is dropped and site-packages' editable
#: finder takes over.
_PROBE = """
import importlib.util, json, sys
for p in json.loads(sys.argv[1]):
    if p not in sys.path:
        sys.path.append(p)
out = {}
for name in json.loads(sys.argv[2]):
    try:
        spec = importlib.util.find_spec(name)
    except (ImportError, AttributeError, ValueError):
        continue
    if spec is not None and spec.origin:
        out[name] = spec.origin
print(json.dumps(out))
"""

#: A Python entry that runs a FILE. `-c` and `-m` don't: `-c` puts the
#: cwd on the path, `-m` resolves through sys.path[0] as the caller set
#: it. Everything after a known script path is that script's own
#: arguments.
_PYTHON_BIN = re.compile(r"(?:^|[/\s=])python(?:\d(?:\.\d+)?)?(?:\.exe)?\s")
#: Flags before the script name. Unknown tokens stop the scan: safer to
#: miss a note than to call an argument a script.
_PYTHON_PRE_FLAGS = {"-B", "-c", "-E", "-I", "-i", "-m", "-O", "-OO",
                     "-P", "-s", "-S", "-u", "-v", "-W", "-X", "-q"}
#: Running the gate or pytest is covered: the gate exports PYTHONPATH to
#: the worktree, pytest imports through its own rootdir machinery.
_QUIET_TOOLS = {"pytest", "gate", "pip", "ruff", "mypy", "black"}


def _workspace_packages(workspace: str) -> list[str]:
    """Importable top-level package names that live in the workspace.

    Only directories -- a workspace whose package is a single module
    (``delfin.py``) is not how DELFIN-like checkouts are laid out, and
    scoping to directories keeps stray files out of the probe.
    """
    root = Path(workspace)
    if not root.is_dir():
        return []
    names = []
    for child in sorted(root.iterdir()):
        if ((child / "__init__.py").is_file()
                and child.name.isidentifier()
                and child.name not in {"tests", "docs", "build"}):
            names.append(child.name)
    return names


def shadowed_packages(workspace: str, extra_paths=None,
                      interpreter: str = "") -> dict:
    """Which workspace source packages a subdirectory script would import
    from somewhere else, as ``{name: installed_location}``.

    Resolved by a SHORT SUBPROCESS PROBE, not by mutating this process's
    ``sys.path`` and calling ``find_spec`` here: a pip editable install
    (21.3+) installs a ``__editable__`` meta-path FINDER that intercepts
    the import before ``sys.path`` is even consulted -- swapping
    ``sys.path`` in-process would silently miss it, and importing the
    workspace package here would run its code. The probe interpreter runs
    isolated (``-I``: no PYTHONPATH, no script directory) but with site
    enabled, which is exactly the view a script in ``.gate/`` has once
    its own directory -- and only that directory -- is dropped. Any
    failure returns ``{}``: an unanswerable question must not stop the
    command, it just means no note.
    """
    names = _workspace_packages(workspace)
    if not names:
        return {}
    argv = [interpreter or sys.executable, "-I", "-c", _PROBE,
            __import__("json").dumps([str(p) for p in (extra_paths or [])]),
            __import__("json").dumps(names)]
    try:
        done = subprocess.run(argv, capture_output=True, text=True,
                              timeout=20)
    except (OSError, subprocess.SubprocessError):
        return {}
    if done.returncode != 0:
        return {}
    try:
        found = __import__("json").loads(done.stdout.strip().splitlines()[-1])
    except (ValueError, IndexError):
        return {}
    root = str(Path(workspace).resolve())
    shadow = {}
    for name, origin in found.items():
        origin = str(origin)
        try:
            inside = str(Path(origin).resolve()).startswith(
                root + os_sep_str())
        except (OSError, ValueError):
            inside = False
        if not inside:
            shadow[name] = origin
    return shadow


def os_sep_str() -> str:
    """One os.sep, split out for the startswith check above (readability
    and a single place to mock in tests on exotic filesystems)."""
    import os
    return os.sep


def _script_path(cmd: str) -> str:
    """The script path a python command runs, or ``""`` when it runs none
    (``-c``, ``-m``, a tool like pytest) or when the shape is unknown."""
    hit = _PYTHON_BIN.search(cmd)
    if not hit:
        return ""
    rest = cmd[hit.end():].strip()
    script = ""
    while rest:
        token = rest.split(None, 1)[0]
        rest = rest[len(token):].strip()
        if token in {"-c", "-m"} or token.startswith(("-c", "-m")) and \
                len(token) > 2:
            return ""
        if token.startswith("--"):
            if "=" in token:
                continue
            return ""  # unknown long flag: give up rather than guess
        if token.startswith("-"):
            flag = token[:2]
            if flag in {"-W", "-X"}:
                continue  # these take a value
            if token in _PYTHON_PRE_FLAGS:
                continue
            return ""  # unknown flag
        script = token
        break
    return script


def note_for_command(cmd: str, cwd: str, workspace: str,
                     _shadow=None) -> str:
    """A short English note when *cmd* runs a Python SCRIPT from a
    subdirectory of the worktree while a workspace package is shadowed;
    ``""`` otherwise.

    Quiet for ``python -c`` / ``python -m`` (their path[0] is the cwd or
    the caller's, not the script's), for pytest and the gate (the gate
    exports PYTHONPATH to the worktree), and whenever the command already
    sets PYTHONPATH to cover the workspace. ``_shadow`` lets a caller
    reuse a cached ``shadowed_packages`` result -- the probe is not free,
    and api_client may run it once per session instead of per command.
    """
    cmd = (cmd or "").strip()
    if not cmd or "PYTHONPATH" in cmd:
        return ""
    script = _script_path(cmd)
    if not script:
        return ""
    try:
        script_abs = Path(script) if Path(script).is_absolute() else \
            Path(cwd or ".") / script
        script_abs = script_abs.resolve()
        root = Path(workspace).resolve()
    except (OSError, ValueError):
        return ""
    if script_abs.parent == root:
        return ""  # script in the root: its sys.path[0] IS the workspace
    try:
        script_abs.relative_to(root)
    except ValueError:
        return ""  # script outside the workspace: none of our packages
    shadow = _shadow if _shadow is not None else \
        shadowed_packages(workspace)
    if not shadow:
        return ""
    names = ", ".join(sorted(shadow))
    where = "; ".join(f"{n} <- {shadow[n]}" for n in sorted(shadow))
    return (
        f"[Shell] this script runs from {script_abs.parent.name}/, so "
        "Python puts that directory -- not the worktree -- first on "
        f"sys.path. Its `import` of {names} resolves to the installed "
        f"copy, not this worktree's code ({where}). The result measures "
        "code that is not yours; a fix here cannot show up in it. Run it "
        "as `PYTHONPATH=. python3 <script>` from the worktree root, or "
        "as `python3 -m` with the script on the path, or -- for a probe "
        "like this -- as a test through the gate, which sets PYTHONPATH "
        "itself.")
