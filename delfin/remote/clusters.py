"""The clusters DELFIN sends calculations to: ctun's entries plus DELFIN's own settings.

Input: the ``remote`` block of the user's settings (``~/.delfin_settings.json``,
written only by the user) and ctun's config (its own default location, or
``ctun_config``)::

    "remote": {"clusters": {"kcist": {"work_dir": "~/delfin/jobs",
                                      "install_root": "~/delfin/installs",
                                      "software_dir": "~/delfin/software",
                                      "allowed_partitions": ["small", "batch"],
                                      "default_partition": "small"}},
               "install_sources": ["https://github.com/ComPlat/DELFIN"]}

A cluster is one of DELFIN's when the settings name it; ctun's config says how
it is reached. ``install_sources`` are the repositories a commit must be
published on before DELFIN installs it on a cluster unasked
(``install_sources()``; ComPlat/DELFIN unless the settings say otherwise). Output: a ``Cluster`` per name, with the problems that keep
DELFIN from using it in ``problems`` -- a bad entry is reported, never raised.

Semantics:

* ctun must have the cluster's host key pinned (``strict_host_key_checking:
  'yes'`` and a ``known_hosts_file``, as ``ctun -t <name> hostkey --pin``
  writes them). The logged-in master is the user's whole identity on the
  cluster, and a key accepted on first contact could be anyone's.
* A budget configured in ctun uses ``guard_scope: every-command``, so ctun's
  own guard also finds a submission inside ``bash -c '...'``.
* Every file ctun's ssh reads for the cluster -- ctun's config, its
  ``ssh_config``, the pinned ``known_hosts_file``, the budget script -- lies
  where no agent route reaches: ctun's own directories (``state_dirs``) or
  ``~/.ssh``. A ProxyCommand an agent wrote into one of them would run in
  DELFIN's own process, outside every cage. Not judged where ctun runs ssh
  inside WSL (native Windows): those are paths of another filesystem.
* ``work_dir``, ``install_root`` and ``software_dir`` (OpenMPI and the tools
  every install shares) are directories on the cluster, under its home (``~/...``, or a relative path, which means the same) or absolute;
  never the home itself, ``/`` or a path with ``..``. The work dir bounds what
  the file tools reach without asking, so it must not be everything.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path, PurePosixPath
from typing import Any, Optional, Union

from cluster_tunnel import api as ctun_api
from cluster_tunnel import config as ctun_config_mod

from . import state_dirs

PathLike = Union[str, Path]

DEFAULT_WORK_DIR = "~/delfin/jobs"
DEFAULT_INSTALL_ROOT = "~/delfin/installs"
DEFAULT_SOFTWARE_DIR = "~/delfin/software"
DEFAULT_INSTALL_SOURCES = ("https://github.com/ComPlat/DELFIN",)

_KEYS = ("work_dir", "install_root", "software_dir", "allowed_partitions",
         "default_partition")
_PARTITION = re.compile(r"[A-Za-z0-9][A-Za-z0-9_.-]*")


@dataclass(frozen=True)
class Cluster:
    """One of DELFIN's clusters, as the settings and ctun's config describe it."""

    name: str
    target: str                       # the ssh destination ctun uses; "" if ctun has no entry
    work_dir: str                     # on the cluster, "~/..." or absolute
    install_root: str                 # on the cluster, "~/..." or absolute
    software_dir: str                 # on the cluster: what every install shares
    allowed_partitions: tuple[str, ...]  # empty: any partition
    default_partition: Optional[str]  # None: the cluster's own default
    pinned: bool                      # ctun accepts only the pinned host key
    requires_password: bool           # what a login asks for (ctun's config)
    requires_otp: bool
    problems: tuple[str, ...]         # why DELFIN will not use it yet; empty when usable
    ctun_config: Optional[Path]       # the ctun config file that describes it

    @property
    def usable(self) -> bool:
        return not self.problems


def load(*, settings: Optional[dict] = None,
         ctun_config: Optional[PathLike] = None) -> dict[str, Cluster]:
    """DELFIN's clusters by name, in the order the settings list them."""
    if settings is None:
        from delfin.user_settings import load_settings
        settings = load_settings()
    remote = settings.get("remote") if isinstance(settings, dict) else None
    block = remote.get("clusters") if isinstance(remote, dict) else None
    if not isinstance(block, dict):
        return {}
    config_file = ctun_config_mod.resolve_config_path(
        str(ctun_config) if ctun_config is not None else None)
    try:
        config = ctun_api.load_config(config_file)
        config_error = ""
    except ctun_api.ConfigError as exc:
        config, config_error = None, str(exc)
    return {str(name): _cluster(str(name), entry, config, config_error, config_file)
            for name, entry in block.items()}


def _cluster(name: str, entry: Any, config, config_error: str,
             config_file: Path) -> Cluster:
    problems: list[str] = []
    if entry is None:
        entry = {}
    if not isinstance(entry, dict):
        problems.append("its entry in remote.clusters must be an object of settings")
        entry = {}
    for key in entry:
        if key not in _KEYS:
            problems.append(f"unknown setting {key!r} (known: {', '.join(_KEYS)})")
    work_dir = _remote_dir(entry.get("work_dir", DEFAULT_WORK_DIR), "work_dir", problems)
    install_root = _remote_dir(entry.get("install_root", DEFAULT_INSTALL_ROOT),
                               "install_root", problems)
    software_dir = _remote_dir(entry.get("software_dir", DEFAULT_SOFTWARE_DIR),
                               "software_dir", problems)
    partitions = _partitions(entry.get("allowed_partitions", []), problems)
    default_partition = _default_partition(entry.get("default_partition"), partitions,
                                           problems)

    target, pinned, requires_password, requires_otp = "", False, True, False
    if config is None:
        problems.append(f"ctun's config cannot be read ({config_error})")
    elif name not in config.clusters:
        problems.append(f"ctun's config ({config_file}) has no cluster {name!r}; "
                        f"add it there (`ctun config`)")
    else:
        spec = ctun_api.connection(name, config=config)
        target = spec.target
        pinned = spec.strict_host_key_checking == "yes" and bool(spec.known_hosts_file)
        if not pinned:
            problems.append(f"ctun has not pinned its host key; check the fingerprint "
                            f"and run `ctun -t {name} hostkey --pin`")
        requires_password = config.clusters[name].requires_password
        requires_otp = config.clusters[name].requires_otp
        budget = config.clusters[name].budget
        if budget is not None and budget.guard_scope != "every-command":
            problems.append(
                "ctun's budget guard looks only at the first word of a command; set "
                f"`guard_scope: every-command` under clusters.{name}.budget in ctun's "
                "config so it also finds sbatch inside `bash -c`")
        files = [("ctun's config", config_file), ("ssh_config", spec.ssh_config),
                 ("known_hosts_file", spec.known_hosts_file)]
        if budget is not None:
            files.append(("the budget script", ctun_config_mod.budget_script_path(
                config_file, name, budget.script)))
        problems += _files_out_of_reach(files)
    return Cluster(name=name, target=target, work_dir=work_dir,
                   install_root=install_root, software_dir=software_dir,
                   allowed_partitions=partitions, default_partition=default_partition,
                   pinned=pinned, requires_password=requires_password,
                   requires_otp=requires_otp, problems=tuple(problems),
                   ctun_config=config_file if config is not None else None)


def _remote_dir(value: Any, key: str, problems: list[str]) -> str:
    """``value`` as ``~/...`` or an absolute path, or "" with a problem noted."""
    if not isinstance(value, str) or not value.strip():
        problems.append(f"{key} must be a directory on the cluster, as text")
        return ""
    text = value.strip()
    if any(c in text for c in "\n\r\0"):
        problems.append(f"{key} must be a single line")
        return ""
    if text.startswith("~") and not (text == "~" or text.startswith("~/")):
        problems.append(f"{key} {text!r}: ~user is not supported, write ~/... or a full path")
        return ""
    absolute = text.startswith("/")
    rest = text[1:] if absolute else (text[2:] if text.startswith("~/") else
                                      "" if text == "~" else text)
    path = PurePosixPath(rest.lstrip("/"))  # "~//a" is ~/a, not /a under ~
    if ".." in path.parts:
        problems.append(f"{key} {text!r} must not contain '..'")
        return ""
    if str(path) in ("", "."):
        problems.append(f"{key} {text!r} must be a directory below "
                        f"{'/' if absolute else 'the home directory'}, not all of it")
        return ""
    return ("/" if absolute else "~/") + str(path)


def _partitions(value: Any, problems: list[str]) -> tuple[str, ...]:
    if value is None:
        return ()
    if not isinstance(value, list) or not all(isinstance(p, str) for p in value):
        problems.append("allowed_partitions must be a list of partition names")
        return ()
    bad = [p for p in value if not _PARTITION.fullmatch(p)]
    if bad:
        problems.append(f"allowed_partitions has names Slurm would not accept: {bad}")
        return ()
    return tuple(value)


def _default_partition(value: Any, allowed: tuple[str, ...],
                       problems: list[str]) -> Optional[str]:
    if value is None:
        return None
    if not isinstance(value, str) or not _PARTITION.fullmatch(value):
        problems.append(f"default_partition must be a partition name, not {value!r}")
        return None
    if allowed and value not in allowed:
        problems.append(f"default_partition {value!r} is not in allowed_partitions {list(allowed)}")
        return None
    return value


def install_sources(settings: Optional[dict] = None) -> tuple[str, ...]:
    """The repositories (https URLs, no trailing ``.git``) a commit must be
    published on for DELFIN to install it on a cluster unasked. An entry
    that is not an https URL is left out, never guessed at."""
    if settings is None:
        from delfin.user_settings import load_settings
        settings = load_settings()
    remote = settings.get("remote") if isinstance(settings, dict) else None
    value = remote.get("install_sources") if isinstance(remote, dict) else None
    if value is None:
        return DEFAULT_INSTALL_SOURCES
    if not isinstance(value, list):
        return ()
    out = []
    for url in value:
        if not isinstance(url, str) or not re.fullmatch(r"https://[^\s'\"\\]+", url.strip()):
            continue
        url = url.strip().rstrip("/")
        url = url[:-4] if url.endswith(".git") else url
        if url not in out:
            out.append(url)
    return tuple(out)


def _files_out_of_reach(files: list[tuple[str, Any]]) -> list[str]:
    """A problem for each file that lies where an agent route could write it."""
    try:
        from cluster_tunnel import wsl as ctun_wsl
        if ctun_wsl.enabled():
            return []
    except Exception:
        pass
    roots = [Path(d) for d in state_dirs()]
    try:
        roots.append((Path.home() / ".ssh").resolve())
    except OSError:
        pass
    problems = []
    for label, path in files:
        if not path:
            continue
        try:
            real = Path(path).expanduser().resolve()
        except OSError:
            problems.append(f"{label} {path} cannot be resolved")
            continue
        if not any(real == r or r in real.parents for r in roots):
            problems.append(
                f"{label} {path} lies outside ctun's directories and ~/.ssh, where an "
                f"agent's command could rewrite it; move it there")
    return problems
