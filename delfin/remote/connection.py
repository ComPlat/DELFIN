"""DELFIN's one way to act on a cluster: ctun's Python API, behind DELFIN's checks.

Input: a cluster name from the registry (``clusters.load``). Output: a
``Connection`` whose calls return ctun's results (``Result``: exit code,
stdout, stderr; ``Decision`` for the budget; ``Job`` rows) and raise only
``delfin.remote.errors``.

Semantics:

* ``connect`` refuses inside an agent session (``DELFIN_AGENT_SESSION``): a
  command an agent runs does not drive the user's ssh master, even where no
  cage hides the socket. The agent reaches clusters through the tools
  DELFIN's own process offers.
* ``connect`` refuses a name the settings do not list (``UnknownCluster``)
  and a listed cluster with problems (``ClusterNotUsable``): no pinned host
  key, a budget guard that reads only first words, a file ctun reads where
  an agent could write it.
* Commands are argv lists: each word reaches the cluster's shell as one
  argument. Output is always captured (``run`` has a timeout; nothing
  streams into DELFIN's own terminal).
* ctun's errors arrive as DELFIN's: LoginRequired (10), BudgetExhausted (11),
  BudgetGuardError (12), HostKeyMismatch, TransportError, and
  ClusterNotUsable for a config ctun itself refuses.
"""

from __future__ import annotations

from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional, Sequence, Union

from cluster_tunnel import api as ctun_api

from . import _inside_agent_session
from . import clusters as _clusters
from .errors import (BudgetExhausted, BudgetGuardError, ClusterNotUsable,
                     HostKeyMismatch, LoginFailed, LoginRequired, RemoteError,
                     RemoteRefused, TransportError, UnknownCluster)

PathLike = Union[str, Path]

#: Seconds a command may take before DELFIN stops waiting (the remote
#: process may keep running). Generous: a login node under load is slow.
DEFAULT_TIMEOUT = 300


def connect(name: str, *, settings: Optional[dict] = None,
            ctun_config: Optional[PathLike] = None) -> "Connection":
    """The connection to DELFIN's cluster ``name``, once every check passes."""
    if _inside_agent_session():
        raise RemoteRefused(
            "the cluster connection is DELFIN's own; a command an agent runs "
            "reaches a cluster only through the cluster tools")
    registry = _clusters.load(settings=settings, ctun_config=ctun_config)
    if name not in registry:
        known = ", ".join(registry) or "none"
        raise UnknownCluster(
            f"{name!r} is not one of DELFIN's clusters (remote.clusters in "
            f"~/.delfin_settings.json lists: {known})")
    cluster = registry[name]
    if not cluster.usable:
        raise ClusterNotUsable(name, cluster.problems)
    return Connection(cluster)


class Connection:
    """One usable cluster, driven through ctun's master."""

    def __init__(self, cluster: "_clusters.Cluster"):
        self.cluster = cluster

    @property
    def name(self) -> str:
        return self.cluster.name

    def _where(self) -> dict:
        return {"config_path": self.cluster.ctun_config}

    def login(self, password: str = "", otp: str = "", *, timeout: int = 120) -> bool:
        """Open the tunnel with a password and OTP a person typed into DELFIN's
        own masked field; they go to ctun (and through it to ssh's prompts)
        and nowhere else. True when the tunnel was opened now, False when it
        was up already. Raises LoginFailed with a message that holds neither."""
        with _ctun_errors(self.name):
            result = ctun_api.login(self.name,
                                    credentials=ctun_api.Credentials(password, otp),
                                    timeout=timeout, **self._where())
        return not result.already_live

    def logout(self) -> bool:
        """Close the tunnel; True when one was open."""
        with _ctun_errors(self.name):
            (result,) = ctun_api.logout(self.name, **self._where())
        return bool(getattr(result, "was_live", True))

    def is_live(self) -> bool:
        """Whether a person is logged in, i.e. ctun's master is alive."""
        with _ctun_errors(self.name):
            return ctun_api.is_live(self.name, **self._where())

    def run(self, argv: Sequence[str], *, timeout: Optional[float] = DEFAULT_TIMEOUT):
        """Run ``argv`` on the login node; the captured ``Result``."""
        with _ctun_errors(self.name):
            return ctun_api.run(self.name, list(argv), capture=True,
                                timeout=timeout, **self._where())

    def feed_script(self, script: PathLike, args: Sequence[object] = (), *,
                    timeout: Optional[float] = DEFAULT_TIMEOUT):
        """Run a local bash script on the login node, fed over the tunnel."""
        with _ctun_errors(self.name):
            return ctun_api.feed_script(self.name, script, args, timeout=timeout,
                                        **self._where())

    def upload(self, src: PathLike, dest: str, *, extra: Sequence[str] = ()):
        """Copy local ``src`` to ``dest`` on the cluster (rsync -r; ``extra``
        are further rsync flags, such as --files-from)."""
        with _ctun_errors(self.name):
            return ctun_api.upload(self.name, src, dest, extra=list(extra), capture=True,
                                   **self._where())

    def download(self, src: str, dest: PathLike, *, extra: Sequence[str] = ()):
        """Copy ``src`` on the cluster to local ``dest`` (rsync -r; ``extra``
        are further rsync flags, such as --files-from)."""
        with _ctun_errors(self.name):
            return ctun_api.download(self.name, src, dest, extra=list(extra), capture=True,
                                     **self._where())

    def jobs(self, *, since_minutes: Optional[int] = None):
        """The user's Slurm jobs on the cluster: running, pending, and those
        that ended in the last ``since_minutes`` (ctun's day by default)."""
        window = {} if since_minutes is None else {"since_minutes": int(since_minutes)}
        with _ctun_errors(self.name):
            (entry,) = ctun_api.jobs(self.name, **window, **self._where())
        if entry.error:
            if not entry.live:
                raise LoginRequired(f"no live tunnel to {self.name!r}: {entry.error}")
            raise TransportError(f"listing the jobs on {self.name!r} failed: {entry.error}")
        return entry.jobs

    def budget(self):
        """The session budget as it stands, probed over the tunnel now."""
        with _ctun_errors(self.name):
            return ctun_api.budget(self.name, **self._where())


@contextmanager
def _ctun_errors(name: str) -> Iterator[None]:
    """ctun's exceptions as DELFIN's, with ctun's message kept."""
    try:
        yield
    except RemoteError:
        raise
    except ctun_api.LoginRequired as exc:
        raise LoginRequired(str(exc)) from exc
    except ctun_api.BudgetExhausted as exc:
        raise BudgetExhausted(str(exc)) from exc
    except ctun_api.BudgetGuardError as exc:
        raise BudgetGuardError(str(exc)) from exc
    except ctun_api.HostKeyError as exc:
        raise HostKeyMismatch(str(exc)) from exc
    except ctun_api.LoginError as exc:
        raise LoginFailed(str(exc)) from exc
    except (ctun_api.ConfigError, ctun_api.UnknownClusterError) as exc:
        raise ClusterNotUsable(name, (str(exc),)) from exc
    except ctun_api.ClusterTunnelError as exc:
        raise TransportError(str(exc)) from exc
    except OSError as exc:
        raise TransportError(f"ssh could not be started: {exc}") from exc
