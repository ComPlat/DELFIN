"""Set the terminal window up where it is first used: a pinned Bun and the
window's npm packages, into directories DELFIN owns.

Input: the window's sources (``launcher.tui_source_root()``: ``tui/`` of a
checkout, or the copy a wheel carries as package data), this machine
(system, architecture, AVX2), the environment (proxy variables, the
opt-out). Output: a window ``launcher.tui_command()`` finds -- Bun at
``~/.delfin/tools/bun-<version>/bun`` and ``node_modules`` in the folder
the window runs in (``launcher.tui_root()``).

What is downloaded, from where, checked how:

- Bun ``BUN_VERSION``: ``https://github.com/oven-sh/bun/releases/download/
  bun-v<version>/<asset>.zip``. Its size and SHA-256 are pinned here (the
  values of ``tui/scripts/install-bun.sh``, taken from the release's
  ``SHASUMS256.txt`` and checked against GitHub's per-asset digest; a
  test holds the two tables equal). The archive is read at most to the
  pinned size, hashed while it is written, and opened only when size and
  hash match; of it only the member ``<asset>/bun`` is extracted, and it
  is kept only when ``bun --version`` names the pinned version.
- The packages: ``bun install --frozen-lockfile --ignore-scripts
  --production`` in the run folder. ``--frozen-lockfile``: exactly the
  versions of ``bun.lock``, each tarball checked against its integrity
  hash there, and a lockfile that does not match ``package.json`` is an
  error, not an update. ``--ignore-scripts``: no lifecycle script runs --
  measured on the current lockfile, no dependency has an install script
  (``bun pm untrusted``: none; the five ``prepare`` scripts are not run
  for registry packages), so the flag changes nothing in what is
  installed and refuses one a later lockfile could bring. The patches of
  ``patchedDependencies`` are applied by ``bun install`` itself, scripts
  or not. ``--production``: the window needs no devDependency
  (typescript, @types/bun, @xterm/headless are for tests); 103 packages
  instead of 108.

Refused: as root, inside an agent session (``DELFIN_AGENT_SESSION``, which
every agent child process carries), in an agent's protected process
(``process_guard.is_protected``), over anything but https (a redirect to
plain http included). The setup runs in the person's own ``delfin-agent
tui`` process, before the window starts.

Not covered: the window then runs the prebuilt OpenTUI native library of
``@opentui/core``; the patched one (``scripts/build-opentui-native.sh``,
opentui#1630: rows can leave the scrollback of an xterm.js terminal when
the footer grows) needs Zig and a build of minutes and stays a manual
step.

Measured on linux-x64 (2026-10-15, a fresh HOME, no Bun): Bun 42.0 MB
download, 89 MB on disk; the packages 103 tarballs, 17.7 MB download,
about 90 MB in node_modules; the whole setup 1.6 s wall on this host's
link. Bun's package cache (``~/.delfin/tools/bun-cache``) is removed after
a setup: node_modules holds hard links or copies of it, and on a home
with a quota a second copy is what costs.

Rejected alternatives: Bun's own installer (``curl | sh``: an unpinned
script); node_modules or binaries in git or the wheel (89 MB per platform,
and the native library differs per platform); any Bun from the network
without a pinned hash; installing from inside an agent's turn (its sandbox
does not reach ``~/.delfin/tools``, and an agent must not change what the
person's window runs).
"""

from __future__ import annotations

import contextlib
import hashlib
import json
import os
import platform
import shutil
import subprocess
import sys
import tempfile
import time
import urllib.error
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path

BUN_VERSION = "1.4.3"

#: asset -> (bytes, SHA-256) of ``<asset>.zip`` of Bun ``BUN_VERSION``.
BUN_ASSETS: dict[str, tuple[int, str]] = {
    "bun-linux-x64": (42015139, "3e730528d89775f03f87f05c06787c4df4099d0a37f3f252e7b829cb23816e18"),
    "bun-linux-x64-baseline": (42015103, "1fc2edac843102909e3a1be1d8d9802cc6071cf074e67e88231f4ffe0f8b397b"),
    "bun-linux-aarch64": (41786424, "efa9813da5ed72423bf847f916e8d2c47c0d776add972354026a75e10da9aa21"),
    "bun-darwin-aarch64": (30530705, "80afd4c069b46a6fa8f3ec520bebc750a7843e4a8d5b9583d29ad436c00a3403"),
    "bun-darwin-x64": (33229351, "62660e6438eb872278877acca77e1152d5d031324df368d4a81b428211489d2d"),
}

BUN_URL = "https://github.com/oven-sh/bun/releases/download/bun-v{version}/{asset}.zip"

#: The packages' download, measured: the 103 tarballs a clean
#: ``bun install --production`` of the current lockfile fetched on
#: linux-x64, 17.7 MB.
PACKAGES_DOWNLOAD_MB = 18

#: What a run folder must hold for the window to start from source.
REQUIRED_PACKAGES = ("@opentui/core", "@opentui/solid", "solid-js")

INSTALL_ARGS = ("install", "--frozen-lockfile", "--ignore-scripts", "--production")

#: The marker of a setup made here, in node_modules: the lockfile it was
#: made for. A folder whose packages someone installed by hand has none
#: and is left alone.
MARKER = ".delfin-setup.json"

ENV_OPT_OUT = "DELFIN_TUI_AUTO_INSTALL"

DOWNLOAD_TIMEOUT_S = 60.0
INSTALL_TIMEOUT_S = 900.0
#: An extracted ``bun`` larger than this is not Bun (89 MB measured).
MAX_BUN_BYTES = 400 * 1024 * 1024


class SetupFailed(RuntimeError):
    """The window could not be set up; the message says what failed and
    how to do it by hand."""


# --- where ---------------------------------------------------------------------

def tools_dir() -> Path:
    return Path.home() / ".delfin" / "tools"


def pinned_bun_path() -> Path:
    return tools_dir() / f"bun-{BUN_VERSION}" / "bun"


def _executable(path: Path) -> bool:
    return path.is_file() and os.access(path, os.X_OK)


def find_bun() -> Path | None:
    """A Bun to run the window with: one on PATH, Bun's own place
    (``~/.bun/bin``), or the one this setup installs -- in that order."""
    on_path = shutil.which("bun")
    for candidate in (on_path, Path.home() / ".bun" / "bin" / "bun", pinned_bun_path()):
        if candidate and _executable(Path(candidate)):
            return Path(candidate)
    return None


def _lock_digest(root: Path) -> str:
    try:
        return hashlib.sha256((root / "bun.lock").read_bytes()).hexdigest()
    except OSError:
        return ""


def window_ready(root: Path) -> bool:
    """The window can start from source in ``root``: its entry point and
    its packages are there, and a setup made here was made for this
    lockfile."""
    if not (root / "src" / "main.tsx").is_file():
        return False
    modules = root / "node_modules"
    if not all((modules / pkg / "package.json").is_file() for pkg in REQUIRED_PACKAGES):
        return False
    try:
        marker = json.loads((modules / MARKER).read_text())
    except FileNotFoundError:
        return True
    except (OSError, ValueError):
        return False
    return isinstance(marker, dict) and marker.get("lock") == _lock_digest(root)


def _write_marker(root: Path, *, bun: Path) -> None:
    (root / "node_modules" / MARKER).write_text(json.dumps({
        "lock": _lock_digest(root), "bun": str(bun), "at": time.time()}))


#: What a copy of the sources holds: what the window runs, and what the
#: manual steps and the patched native build need.
_SKIP = {"node_modules", "dist", "test", "spike", ".git"}


def _source_files(source: Path) -> list[Path]:
    files = []
    for path in sorted(source.rglob("*")):
        rel = path.relative_to(source)
        if rel.parts and rel.parts[0] in _SKIP or not path.is_file():
            continue
        files.append(rel)
    return files


def content_digest(source: Path) -> str:
    """SHA-256 over the relative paths and contents of the sources: the
    name of the folder a copy of them runs in."""
    h = hashlib.sha256()
    for rel in _source_files(source):
        h.update(rel.as_posix().encode() + b"\0")
        h.update((source / rel).read_bytes())
        h.update(b"\0")
    return h.hexdigest()


def cache_root(source: Path) -> Path:
    """Where a copy of ``source`` runs: one folder per content, so an
    upgrade with other sources never runs the old packages."""
    return tools_dir() / f"tui-{content_digest(source)[:16]}"


# --- whether -------------------------------------------------------------------

def _is_root() -> bool:
    return hasattr(os, "geteuid") and os.geteuid() == 0


def _protected_process() -> bool:
    try:
        from .. import process_guard
        return process_guard.is_protected()
    except Exception:
        return False


def refusal(env=None) -> str:
    """Why the setup must not run in this process, or ""."""
    from .. import approval_guard

    env = os.environ if env is None else env
    if approval_guard.inside_agent_session(env):
        return ("this runs inside an agent session (DELFIN_AGENT_SESSION is set): the window "
                "is set up by you, in your own shell, never by an agent")
    if _protected_process():
        return "this is an agent's protected process: the window is set up from your own shell"
    if _is_root():
        return "this runs as root: set the window up as the user who runs it, not as root"
    return ""


def _settings_auto_install():
    """``agent.tui_auto_install`` of the user's settings: True, False, or
    None when unset or unreadable."""
    try:
        from delfin.user_settings import load_settings
        value = ((load_settings() or {}).get("agent") or {}).get("tui_auto_install")
    except Exception:
        return None
    if isinstance(value, bool):
        return value
    if isinstance(value, str) and value.strip().lower() in ("0", "false", "no", "off"):
        return False
    return None


def auto_install_enabled(env=None, settings="read") -> bool:
    """Whether ``delfin-agent tui`` sets the window up by itself. The
    environment wins over the settings: ``DELFIN_TUI_AUTO_INSTALL=0`` (or
    false/no/off) switches it off, ``=1`` on; else
    ``agent.tui_auto_install=false`` switches it off. On by default."""
    env = os.environ if env is None else env
    raw = str(env.get(ENV_OPT_OUT, "") or "").strip().lower()
    if raw in ("0", "false", "no", "off"):
        return False
    if raw in ("1", "true", "yes", "on"):
        return True
    value = _settings_auto_install() if settings == "read" else settings
    return value is not False


# --- what ----------------------------------------------------------------------

def _cpu_has_avx2() -> bool:
    try:
        with open("/proc/cpuinfo", encoding="ascii", errors="replace") as fh:
            for line in fh:
                if line.startswith("flags"):
                    return " avx2" in line
    except OSError:
        pass
    return True


def bun_asset(system: str | None = None, machine: str | None = None, *, avx2=None) -> str | None:
    """The Bun release asset for this machine, or None. A Linux x86_64
    CPU without AVX2 gets the baseline build (the other one stops there
    with "Illegal instruction")."""
    system = platform.system() if system is None else system
    machine = (platform.machine() if machine is None else machine).lower()
    arch = {"x86_64": "x64", "amd64": "x64", "aarch64": "aarch64", "arm64": "aarch64"}.get(machine)
    if system == "Linux" and arch == "x64":
        if avx2 is None:
            avx2 = _cpu_has_avx2()
        return "bun-linux-x64" if avx2 else "bun-linux-x64-baseline"
    if system == "Linux" and arch == "aarch64":
        return "bun-linux-aarch64"
    if system == "Darwin" and arch in ("x64", "aarch64"):
        return f"bun-darwin-{arch}"
    return None


@dataclass
class Plan:
    """What a setup does: the sources, the folder the window runs in
    (``copy``: the sources are copied there first), the Bun it uses (None:
    ``asset`` is downloaded), and whether the packages are installed."""
    source: Path
    root: Path
    copy: bool
    bun: Path | None
    asset: str | None
    packages: bool


def plan() -> Plan:
    from . import launcher

    source = launcher.tui_source_root()
    if source is None:
        raise SetupFailed(
            "the window's sources are not part of this installation (neither tui/ beside the "
            "delfin package nor delfin/agent/service/" + launcher.PACKAGED_TUI + "/). "
            "Reinstall DELFIN, or set DELFIN_TUI_BINARY to a delfin-tui-" + launcher.platform_name()
            + " binary.")
    root = launcher.tui_root()
    bun = find_bun()
    return Plan(source=source, root=root, copy=root != source, bun=bun,
                asset=None if bun else bun_asset(), packages=not window_ready(root))


def _tilde(path) -> str:
    text, home = str(path), str(Path.home())
    return "~" + text[len(home):] if text == home or text.startswith(home + os.sep) else text


def announce(p: Plan, *, explicit: bool = False) -> str:
    """The one line said before anything is downloaded; the opt-out is
    named unless the setup was asked for (``--setup``)."""
    parts = []
    if p.bun is None and p.asset:
        size = BUN_ASSETS[p.asset][0] / 1e6
        parts.append(f"Bun {BUN_VERSION} ({size:.0f} MB from github.com/oven-sh/bun, SHA-256 "
                     f"pinned) into {_tilde(pinned_bun_path().parent)}")
    if p.packages:
        parts.append(f"its npm packages (about {PACKAGES_DOWNLOAD_MB} MB, each checked against "
                     f"bun.lock, no install scripts) into {_tilde(p.root / 'node_modules')}")
    tail = "." if explicit else f". Opt out: {ENV_OPT_OUT}=0."
    return "Setting up the terminal window (once): " + "; then ".join(parts) + tail


def manual_commands(p: Plan | None = None) -> str:
    """The setup by hand, as shell lines for the person's own shell."""
    from . import launcher

    if p is None:
        try:
            p = plan()
        except SetupFailed:
            p = None
    name = f"delfin-tui-{launcher.platform_name()}"
    lines = ["  delfin-agent tui --setup                 # all of the below, checked"]
    if p is not None:
        bun = _tilde(p.bun) if p.bun else "~/.bun/bin/bun"
        if p.bun is None:
            lines.append(f"  bash {_tilde(p.source / 'scripts' / 'install-bun.sh')} {BUN_VERSION}"
                         "   # Bun, SHA-256 checked, into ~/.bun/bin")
        if p.copy:
            lines.append(f"  mkdir -p {_tilde(p.root)} && cp -R {_tilde(p.source)}/. {_tilde(p.root)}/")
        lines.append(f"  cd {_tilde(p.root)} && {bun} {' '.join(INSTALL_ARGS)}")
    lines.append(f"Or, without network: DELFIN_TUI_BINARY=/path/to/{name} (built with "
                 f"`bun scripts/build-binary.ts` in tui/ on a machine with network).")
    return "\n".join(lines)


def not_set_up(reason: str, p: Plan | None = None) -> str:
    return (f"The terminal window is not set up here, and {reason}. To set it up:\n"
            + manual_commands(p))


# --- how -----------------------------------------------------------------------

class _HttpsOnlyRedirect(urllib.request.HTTPRedirectHandler):
    """A redirect is followed only to another https address."""

    def redirect_request(self, req, fp, code, msg, headers, newurl):
        if urllib.parse.urlsplit(newurl).scheme != "https":
            raise urllib.error.URLError(f"a redirect to a non-https address was refused: {newurl[:200]}")
        return super().redirect_request(req, fp, code, msg, headers, newurl)


def _open(url: str, timeout: float):
    """urllib with the environment's proxies (HTTPS_PROXY, NO_PROXY: its
    default ProxyHandler) and https-only redirects."""
    from ... import __version__
    opener = urllib.request.build_opener(_HttpsOnlyRedirect)
    req = urllib.request.Request(url, headers={"User-Agent": f"delfin-tui-setup/{__version__}"})
    return opener.open(req, timeout=timeout)


def _offline_hint(host: str) -> str:
    return (f"could not reach {host}: this host may be offline or need a proxy "
            "(HTTPS_PROXY=http://proxy:port; NO_PROXY for exceptions). On a host without "
            "network, set DELFIN_TUI_BINARY to a binary built elsewhere")


def _fetch(url: str, dest: Path, size: int, sha256: str, opener) -> None:
    if urllib.parse.urlsplit(url).scheme != "https":
        raise SetupFailed(f"refused to download {url}: only https is used")
    host = urllib.parse.urlsplit(url).hostname or url
    h = hashlib.sha256()
    n = 0
    try:
        with opener(url, DOWNLOAD_TIMEOUT_S) as resp, open(dest, "wb") as out:
            length = (resp.headers or {}).get("Content-Length")
            try:
                if length is not None and int(length) >= 0 and int(length) != size:
                    raise SetupFailed(f"{url} announces {length} bytes; the pinned size is {size}: "
                                      "refused, nothing was extracted")
            except ValueError:
                pass
            while True:
                chunk = resp.read(1 << 20)
                if not chunk:
                    break
                n += len(chunk)
                if n > size:
                    raise SetupFailed(f"{url} sent more than its pinned size of {size} bytes: "
                                      "refused, nothing was extracted")
                h.update(chunk)
                out.write(chunk)
    except SetupFailed:
        raise
    except (urllib.error.URLError, OSError, TimeoutError) as exc:
        reason = getattr(exc, "reason", exc)
        raise SetupFailed(f"{_offline_hint(host)} ({reason})") from exc
    got = h.hexdigest()
    if n != size or got != sha256:
        raise SetupFailed(f"{url}: {n} bytes with SHA-256 {got}; pinned are {size} bytes with "
                          f"SHA-256 {sha256}. Refused: nothing was extracted.")


def download_bun(asset: str, *, opener=None, run=None) -> Path:
    """Bun ``BUN_VERSION`` for ``asset`` into ``pinned_bun_path()``,
    checked (module docstring). Returns its path."""
    import zipfile

    opener = opener or _open
    run = run or subprocess.run
    if asset not in BUN_ASSETS:
        raise SetupFailed(f"no pinned Bun {BUN_VERSION} for {asset}")
    size, sha = BUN_ASSETS[asset]
    url = BUN_URL.format(version=BUN_VERSION, asset=asset)
    final = pinned_bun_path()
    tools = _ensure_private_dir(tools_dir())
    work = Path(tempfile.mkdtemp(prefix=".download-", dir=tools))
    try:
        archive = work / f"{asset}.zip"
        _fetch(url, archive, size, sha, opener)
        member = f"{asset}/bun"
        try:
            with zipfile.ZipFile(archive) as z:
                info = z.getinfo(member)
                if info.is_dir() or info.file_size > MAX_BUN_BYTES:
                    raise SetupFailed(f"{url}: {member} is not a plausible Bun executable")
                staged = work / "bun"
                with z.open(info) as src, open(staged, "wb") as out:
                    shutil.copyfileobj(src, out, 1 << 20)
        except KeyError:
            raise SetupFailed(f"{url} holds no {member}; nothing was installed") from None
        except zipfile.BadZipFile as exc:
            raise SetupFailed(f"{url} is not a zip archive ({exc})") from None
        os.chmod(staged, 0o755)
        try:
            proc = run([str(staged), "--version"], capture_output=True, text=True,
                       timeout=60, stdin=subprocess.DEVNULL)
            said = ((proc.stdout or "") + (proc.stderr or "")).strip()
            ok = proc.returncode == 0 and said.splitlines()[:1] == [BUN_VERSION]
        except (OSError, subprocess.SubprocessError) as exc:
            ok, said = False, str(exc)
        if not ok:
            raise SetupFailed(f"the downloaded Bun ({asset}) does not run here: {said[:300] or 'no output'}")
        _ensure_private_dir(final.parent)
        os.replace(staged, final)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return final


def _ensure_private_dir(path: Path) -> Path:
    """Created owner-only (0700), as every DELFIN state folder is."""
    from ..state_paths import ensure_dir
    ensure_dir(path)
    if not path.is_dir():
        raise SetupFailed(f"cannot create {_tilde(path)}")
    return path


_SECRET_MARKERS = ("_API_KEY", "_SECRET", "_PASSWORD", "_TOKEN")


def install_env() -> dict[str, str]:
    """``bun install``'s environment: the person's, without credentials
    (provider keys, tokens), with Bun's telemetry off and its package
    cache in DELFIN's tools folder. Proxy variables are kept."""
    env = {k: v for k, v in os.environ.items()
           if not any(m in k.upper() for m in _SECRET_MARKERS)}
    env["DO_NOT_TRACK"] = "1"
    env["BUN_INSTALL_CACHE_DIR"] = str(tools_dir() / "bun-cache")
    return env


def install_packages(p: Plan, bun: Path, *, run=None, root: Path | None = None) -> None:
    run = run or subprocess.run
    where = root or p.root
    argv = [str(bun), *INSTALL_ARGS]
    try:
        proc = run(argv, cwd=str(where), env=install_env(), capture_output=True, text=True,
                   timeout=INSTALL_TIMEOUT_S, stdin=subprocess.DEVNULL)
        output = ((proc.stdout or "") + (proc.stderr or "")).strip()
        code = proc.returncode
    except (OSError, subprocess.SubprocessError) as exc:
        output, code = str(exc), -1
    if code != 0 or not all((where / "node_modules" / pkg / "package.json").is_file()
                            for pkg in REQUIRED_PACKAGES):
        tail = "\n".join(output.splitlines()[-12:])
        raise SetupFailed(f"`{' '.join(argv)}` in {_tilde(where)} failed (exit {code}):\n{tail}\n"
                          "By hand:\n" + manual_commands(p))
    _write_marker(where, bun=bun)
    shutil.rmtree(tools_dir() / "bun-cache", ignore_errors=True)


@contextlib.contextmanager
def _setup_lock():
    """One setup at a time per home (two windows started at once); a
    filesystem without locks (some NFS) goes without."""
    _ensure_private_dir(tools_dir())
    fh = open(tools_dir() / ".setup.lock", "a")
    try:
        try:
            import fcntl
            fcntl.flock(fh, fcntl.LOCK_EX)
        except Exception:
            pass
        yield
    finally:
        fh.close()


def set_up(*, out=None, opener=None, run=None, explicit: bool = False) -> Plan:
    """Set the window up (module docstring); says in one line what it
    does. Raises SetupFailed with what failed and the manual commands."""
    out = out if out is not None else sys.stderr
    why = refusal()
    if why:
        raise SetupFailed(f"not set up: {why}.\nBy hand, in your own shell:\n" + manual_commands())
    with _setup_lock():
        p = plan()
        if not p.packages and p.bun is not None:
            return p                               # another window set it up meanwhile
        if p.bun is None and p.asset is None:
            from . import launcher
            raise SetupFailed(f"there is no Bun release for this machine ({launcher.platform_name()}); "
                              "set DELFIN_TUI_BINARY to a binary built for it")
        print(announce(p, explicit=explicit), file=out, flush=True)
        started = time.monotonic()
        try:
            bun = p.bun or download_bun(p.asset, opener=opener, run=run)
            if p.packages:
                if p.copy:
                    _install_copy(p, bun, run=run)
                else:
                    install_packages(p, bun, run=run)
        except SetupFailed as exc:
            text = str(exc)
            if "By hand" not in text:
                text += "\nBy hand, in your own shell:\n" + manual_commands(p)
            raise SetupFailed(text) from None
        print(f"the terminal window is set up ({time.monotonic() - started:.0f} s)", file=out, flush=True)
        return p


def _install_copy(p: Plan, bun: Path, *, run=None) -> None:
    """The sources copied into a staging folder, the packages installed
    there, and the folder renamed into place: a failed setup leaves no
    half-made run folder."""
    tools = _ensure_private_dir(tools_dir())
    stage = Path(tempfile.mkdtemp(prefix=".tui-", dir=tools))
    try:
        work = stage / "tui"
        for rel in _source_files(p.source):
            target = work / rel
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(p.source / rel, target)
        install_packages(p, bun, run=run, root=work)
        if p.root.exists():
            shutil.rmtree(p.root)
        os.replace(work, p.root)
    finally:
        shutil.rmtree(stage, ignore_errors=True)


def native_note(root: Path) -> str:
    """On Linux, when the window's OpenTUI native library is the prebuilt
    one: how to build the patched one; else ""."""
    from . import launcher

    plat = launcher.platform_name()
    if not plat.startswith("linux-"):
        return ""
    if (root / "node_modules" / "@opentui" / f"core-{plat}" / ".delfin-patched").is_file():
        return ""
    return ("note: the window runs OpenTUI's prebuilt native library; in xterm.js terminals "
            "(VS Code, JupyterLab) rows can leave the scrollback when the footer grows. "
            f"`cd {_tilde(root)} && scripts/build-opentui-native.sh` builds the patched one "
            "(Zig 0.16 and the OpenTUI source, SHA-256 pinned; minutes).")


def describe() -> tuple[str, str]:
    """For the doctor: ("binary"|"source"|"missing", detail)."""
    from . import launcher

    binary = launcher.tui_binary()
    if binary is not None:
        return "binary", f"the window's binary {_tilde(binary)}"
    command = launcher.tui_command()
    if command is not None:
        return "source", (f"from source: Bun {_tilde(command[0])}, packages in "
                          f"{_tilde(launcher.tui_root() / 'node_modules')}")
    return "missing", "not set up here (no binary, no Bun with the window's packages)"
