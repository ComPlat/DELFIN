"""A cluster that needs a login, or one an agent proposes, as an item in the
user's attention list.

``ask_for_login`` records a ``cluster_login`` item when a cluster tool finds
no live tunnel: the window lists it under NEEDS YOU and opens its masked
login field for it (never a text answer, which would go to the model). One
pending item per cluster: asking again while one is open adds none.
``logged_in`` resolves a cluster's items once a login opened the tunnel,
without an answer and acknowledged, so nothing of it is handed to a model.

``propose_cluster`` records a ``cluster_add`` item (cluster_propose): the
window opens its add form filled in with the proposal, and the person
checks and sends it -- nothing else writes a cluster. ``added`` resolves a
cluster's proposals once it was added, the same way.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union


def ask_for_login(cluster: str, *, reason: str = "", session_id: str = "",
                  workspace: Optional[Union[str, Path]] = None) -> Optional[str]:
    """The id of the cluster's pending login item, recorded now if there was
    none; None when the inbox cannot be written."""
    from delfin.agent import attention
    try:
        for event in attention.list_pending("cluster_login"):
            if event.get("cluster") == cluster:
                return str(event.get("id"))
        return attention.emit_attention(
            "cluster_login", cluster=cluster, session_id=session_id, workspace=workspace,
            title=f"Log in to {cluster}",
            detail=reason or f"A cluster tool needs a login to {cluster}.")
    except Exception:
        return None


def logged_in(cluster: str) -> int:
    """Resolve the cluster's pending login items; how many there were."""
    from delfin.agent import attention
    done = 0
    try:
        for event in attention.list_pending("cluster_login"):
            if event.get("cluster") == cluster and attention.resolve(
                    str(event.get("id")), answer=None, acknowledged=True):
                done += 1
    except Exception:
        pass
    return done


def propose_cluster(proposal: dict, *, session_id: str = "",
                    workspace: Optional[Union[str, Path]] = None) -> Optional[str]:
    """The id of the pending proposal for ``proposal["name"]`` (a new one
    replaces an earlier one's values); None when the inbox cannot be written."""
    from delfin.agent import attention
    name = str(proposal.get("name") or "")
    try:
        for event in attention.list_pending("cluster_add"):
            if event.get("cluster") == name:
                attention.resolve(str(event.get("id")), answer=None, acknowledged=True)
        what = proposal.get("preset") or proposal.get("host") or ""
        return attention.emit_attention(
            "cluster_add", cluster=name, session_id=session_id, workspace=workspace,
            proposal=proposal, title=f"Add the cluster {name}",
            detail=f"An agent proposes {name}" + (f" ({what})" if what else "")
                   + " for DELFIN: check it in the form and send it, or dismiss it.")
    except Exception:
        return None


def added(cluster: str) -> int:
    """Resolve the cluster's pending proposals; how many there were."""
    from delfin.agent import attention
    done = 0
    try:
        for event in attention.list_pending("cluster_add"):
            if event.get("cluster") == cluster and attention.resolve(
                    str(event.get("id")), answer=None, acknowledged=True):
                done += 1
    except Exception:
        pass
    return done
