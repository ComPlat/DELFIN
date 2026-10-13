"""The git role, resolved for the repository a command acts on.

Input: a repository directory, and optionally the remote word or URL a
command names. Output: ``maintainer`` or ``contributor``.

Settings (``~/.delfin_settings.json``, written only by the user):

* ``agent.git_roles`` -- ``{"<key>": "maintainer"|"contributor"}``. A key
  is a normalised remote URL (``github.com/owner/repo``: scheme, user,
  port and ``.git`` dropped, host lower-cased, path case kept) or, for a
  repository with no remote, the resolved absolute path of its top level.
* ``agent.git_role`` -- the role of every repository without an entry.

Resolution for one key: its entry, else the global role, else
contributor. An entry always wins, in both directions: ``contributor``
for one repository holds under a global ``maintainer``. A command that
reaches several keys (several push URLs, several remotes) is maintainer
only when every key is. Anything that cannot be resolved -- no
repository, an unparsable URL, a ``git_roles`` value that is not a
mapping, an entry that is neither role -- is contributor.

Rejected: keying by the checkout path alone. Two clones of one repository
would hold different roles, and a clone of a contributor repository
placed at a path the user once made maintainer would inherit it. The
remote is what a push reaches; the path is used only when there is none.

Nothing here writes except ``set_role``, which the terminal's ``/git
role`` command calls on the user's typed request.
"""

from __future__ import annotations

import os
import re
import subprocess
from pathlib import Path
from typing import Any, Callable, Iterable, Optional

ROLES = ("maintainer", "contributor")

#: Variables that make git act on another repository or with other
#: configuration than the one read here. A command run under any of them
#: cannot be resolved; the reader itself runs without them, but keeps the
#: two that only choose which global configuration file applies.
_REDIRECTING_GIT_ENV = ("GIT_DIR", "GIT_WORK_TREE", "GIT_CONFIG_PARAMETERS",
                        "GIT_CONFIG_COUNT", "GIT_NAMESPACE", "GIT_COMMON_DIR")
_KEPT_GIT_ENV = ("GIT_CONFIG_GLOBAL", "GIT_CONFIG_NOSYSTEM")

_SCP_RE = re.compile(r"^(?:[^@/\s]+@)?([^:/\s]+):(?!//)(.+)$")
_URL_RE = re.compile(
    r"^(?:ssh|git|https?|git\+ssh|ssh\+git)://(?:[^@/\s]+@)?"
    r"(\[[^\]]+\]|[^:/\s]+)(?::\d*)?(/.*)$", re.IGNORECASE)


def _git(where: Any, *args: str) -> Optional[str]:
    """stdout of a read-only git command, or None. Never raises."""
    env = {k: v for k, v in os.environ.items()
           if not k.startswith("GIT_") or k in _KEPT_GIT_ENV}
    try:
        r = subprocess.run(["git", *args], cwd=str(where), capture_output=True,
                           text=True, timeout=10, env=env)
    except Exception:
        return None
    return r.stdout.strip() if r.returncode == 0 else None


def env_redirects_git(env: Optional[dict] = None) -> bool:
    """True when the environment points git somewhere this module does not
    look."""
    env = os.environ if env is None else env
    return any(env.get(k) for k in _REDIRECTING_GIT_ENV) or any(
        k.startswith(("GIT_CONFIG_KEY_", "GIT_CONFIG_VALUE_")) for k in env)


def normalise_remote(url: Any, base: Any = None) -> Optional[str]:
    """The settings key for a remote URL; None when it is not one.

    ``git@GitHub.com:Me/Repo.git``, ``ssh://git@github.com:22/Me/Repo``
    and ``https://user@github.com/Me/Repo/`` are all
    ``github.com/Me/Repo``. A local path or ``file://`` URL is the
    resolved absolute path (relative to ``base``).

    A host with a character outside ASCII has no key. ``str.lower`` maps
    some of them onto ASCII (KELVIN SIGN to ``k``, ``İ`` to ``i`` plus a
    combining dot), so a host git or DNS reads as another name would
    match an entry written for this one. Rejected: IDNA encoding, whose
    mapping differs between IDNA 2003 and 2008 and between resolvers."""
    text = str(url or "").strip()
    if not text or any(c.isspace() for c in text):
        return None
    if text.lower().startswith("file://"):
        text = text[len("file://"):]
        return _path_key(text, base)
    m = _URL_RE.match(text)
    if m:
        host, path = m.group(1), m.group(2)
    elif "://" in text:
        return None
    else:
        m = _SCP_RE.match(text)
        if m and "/" not in m.group(1) and not text.startswith(("/", ".")):
            host, path = m.group(1), m.group(2)
        else:
            return _path_key(text, base)
    path = path.strip("/")
    if path.endswith(".git"):
        path = path[:-4]
    path = path.strip("/")
    if not host or not path or ".." in path.split("/") or not host.isascii():
        return None
    return f"{host.lower()}/{path}"


def _path_key(text: str, base: Any) -> Optional[str]:
    try:
        p = Path(text).expanduser()
        if not p.is_absolute():
            p = Path(base or ".") / p
        return str(p.resolve())
    except Exception:
        return None


def _rewrite_url(url: str, where: Any) -> str:
    """``url`` after the repository's ``pushInsteadOf`` / ``insteadOf``
    rules, the way git rewrites a URL given on the command line: the
    longest matching pushInsteadOf, else the longest insteadOf."""
    out = _git(where, "config", "--get-regexp",
               r"^url\..*\.(pushinsteadof|insteadof)$") or ""
    best: dict[str, tuple[int, str]] = {}
    for line in out.splitlines():
        name, _, prefix = line.partition(" ")
        kind = "push" if name.lower().endswith(".pushinsteadof") else "fetch"
        base = name[len("url."):name.rfind(".")]
        if prefix and url.startswith(prefix):
            if len(prefix) > best.get(kind, (-1, ""))[0]:
                best[kind] = (len(prefix), base + url[len(prefix):])
    hit = best.get("push") or best.get("fetch")
    return hit[1] if hit else url


def entry_key(key: Any) -> Optional[str]:
    """A ``git_roles`` key as the user may have written it, normalised:
    ``github.com/Me/Repo`` as it stands (host lower-cased), a URL through
    ``normalise_remote``, an absolute path as it stands."""
    text = str(key or "").strip()
    if text.startswith("/"):
        return text.rstrip("/") or None
    if "://" in text or ":" in text or "@" in text:
        return normalise_remote(text)
    host, _, path = text.partition("/")
    path = path.strip("/")
    if path.endswith(".git"):
        path = path[:-4]
    if (not host or not path or any(c.isspace() for c in text)
            or not host.isascii()):
        return None
    return f"{host.lower()}/{path}"


def toplevel(where: Any) -> Optional[str]:
    top = _git(where, "rev-parse", "--show-toplevel")
    return str(Path(top).resolve()) if top else None


def remote_keys(where: Any, remote: Optional[str] = None) -> Optional[list]:
    """The keys a push from ``where`` to ``remote`` reaches; None when it
    cannot be told.

    ``remote`` None: every push URL of every remote (git and gh pick one
    of them by rules this does not reproduce, so all of them count). No
    remote at all: the top-level path. A remote name: its push URLs. A
    word that names no remote: a URL or path, rewritten as git would."""
    top = toplevel(where)
    if top is None:
        return None
    if remote is None:
        names = (_git(where, "remote") or "").split()
        # An argument-less push takes its remote from these keys, and each
        # may hold a URL instead of a remote name. A value that names no
        # remote fails ``get-url`` below and leaves the push unresolved.
        branch = _git(where, "symbolic-ref", "--quiet", "--short", "HEAD")
        for cfg in ["remote.pushDefault"] + (
                [f"branch.{branch}.pushRemote", f"branch.{branch}.remote"]
                if branch else []):
            value = _git(where, "config", "--get", cfg)
            if value and value not in names:
                names.append(value)
        if not names:
            return [top]
    else:
        names = [remote]
    keys: list = []
    for name in names:
        urls = _git(where, "remote", "get-url", "--push", "--all", name)
        if urls is None:
            if remote is None:
                return None
            urls = _rewrite_url(name, where)
        for u in urls.splitlines():
            key = normalise_remote(u, base=where)
            if key is None:
                return None
            keys.append(key)
    return keys or None


def repo_key(where: Any) -> Optional[str]:
    """The key ``/git role`` writes for the repository at ``where``:
    origin's push URL, else the only remote's, else the top-level path."""
    top = toplevel(where)
    if top is None:
        return None
    names = (_git(where, "remote") or "").split()
    if not names:
        return top
    name = "origin" if "origin" in names else (names[0] if len(names) == 1
                                               else None)
    if name is None:
        return None
    urls = (_git(where, "remote", "get-url", "--push", name) or "").splitlines()
    return normalise_remote(urls[0], base=where) if urls else None


def _settings_agent(load: Optional[Callable] = None) -> dict:
    try:
        if load is None:
            from delfin import user_settings
            load = user_settings.load_settings
        return ((load() or {}).get("agent") or {})
    except Exception:
        return {}


def _role(value: Any) -> Optional[str]:
    role = str(value or "").strip().lower()
    return role if role in ROLES else None


def global_role(agent: Optional[dict] = None) -> str:
    agent = _settings_agent() if agent is None else agent
    return _role(agent.get("git_role", "contributor")) or "contributor"


def role_for_keys(keys: Optional[Iterable[str]], *,
                  fallback: Callable[[], str],
                  agent: Optional[dict] = None) -> tuple:
    """(role, source) for a command reaching ``keys``.

    ``fallback`` returns the installation-wide role for keys without an
    entry. Source is ``"repository"``, ``"installation"`` or
    ``"unresolved"``."""
    keys = list(keys or ())
    if not keys:
        return "contributor", "unresolved"
    agent = _settings_agent() if agent is None else agent
    table = agent.get("git_roles", {})
    if table is None:
        table = {}
    if not isinstance(table, dict):
        return "contributor", "unresolved"
    entries: dict = {}
    for k, v in table.items():
        nk = entry_key(k)
        if nk:
            entries[nk] = v
    roles, sources = [], set()
    for key in keys:
        if key in entries:
            r = _role(entries[key])
            if r is None:
                return "contributor", "unresolved"
            roles.append(r)
            sources.add("repository")
        else:
            try:
                r = _role(fallback()) or "contributor"
            except Exception:
                r = "contributor"
            roles.append(r)
            sources.add("installation")
    role = "maintainer" if all(r == "maintainer" for r in roles) else "contributor"
    return role, ("repository" if "repository" in sources else "installation")


def set_role(key: str, role: str, settings_path: Any = None) -> dict:
    """Write ``agent.git_roles[key] = role``; the saved settings.

    Called by the user's ``/git role`` command only."""
    if role not in ROLES:
        raise ValueError(f"role must be one of {ROLES}")
    if not key:
        raise ValueError("no repository key")
    from delfin import user_settings
    settings = user_settings.load_settings(settings_path)
    agent = dict(settings.get("agent") or {})
    table = agent.get("git_roles")
    table = dict(table) if isinstance(table, dict) else {}
    table[key] = role
    agent["git_roles"] = table
    settings["agent"] = agent
    return user_settings.save_settings(settings, settings_path)
