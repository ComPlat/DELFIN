"""What can stop DELFIN from acting on a cluster, as DELFIN's own errors.

ctun's exceptions are mapped onto these in ``connection`` so that the rest of
DELFIN (the agent's cluster tools, the job watcher, ``delfin doctor``) never
imports ctun's classes. ``code`` is the word a tool result carries
(``login_required`` and so on), stable across ctun releases.
"""

from __future__ import annotations


class RemoteError(Exception):
    """Base of every error of the remote layer."""

    code = "remote_error"


class UnknownCluster(RemoteError):
    """The name is not one of DELFIN's clusters (``remote.clusters``)."""

    code = "unknown_cluster"


class ClusterNotUsable(RemoteError):
    """A configured cluster that DELFIN refuses to use until it is fixed;
    ``problems`` says what to fix, one line each."""

    code = "cluster_not_usable"

    def __init__(self, name: str, problems: tuple[str, ...]):
        self.name = name
        self.problems = tuple(problems)
        super().__init__(f"cluster {name!r} cannot be used yet: " + "; ".join(self.problems))


class RemoteRefused(RemoteError):
    """Asked for from a command an agent runs. The cluster is reached only
    through DELFIN's own process and the cluster tools it offers."""

    code = "remote_refused"


class LoginRequired(RemoteError):
    """No live tunnel: a person has to log in (ctun's exit code 10)."""

    code = "login_required"


class BudgetExhausted(RemoteError):
    """The session's compute budget is used up (ctun's exit code 11)."""

    code = "budget_exhausted"


class BudgetGuardError(RemoteError):
    """The budget could not be checked and the cluster fails closed (12)."""

    code = "budget_guard_error"


class LoginFailed(RemoteError):
    """The login did not open the tunnel: a wrong password or OTP, a
    refused host key, a cluster that did not answer. The message never
    holds a secret."""

    code = "login_failed"


class HostKeyMismatch(RemoteError):
    """The cluster presents a host key other than the pinned one."""

    code = "host_key_mismatch"


class TransportError(RemoteError):
    """ssh or rsync failed for a reason of their own."""

    code = "transport_error"


class CommitNotPublished(RemoteError):
    """The commit the PC runs is on no branch of the install sources, or the
    PC does not run DELFIN from a git checkout at all; DELFIN installs only
    published commits on a cluster unasked."""

    code = "commit_not_published"


class InstallFailed(RemoteError):
    """install_delfin.sh failed on the cluster (``log`` holds its tail), or
    did not finish in time."""

    code = "install_failed"

    def __init__(self, message: str, log: str = ""):
        self.log = log
        super().__init__(message + (f"\n{log}" if log else ""))
