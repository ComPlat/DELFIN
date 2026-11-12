"""Adding and removing DELFIN's clusters, as a person asks for it in the window.

Input: what the window's add form sends (``cluster.add``: a name, a ctun
preset or a host, the user name, DELFIN's partition and directories), a
fingerprint the person confirmed (``cluster.pin``), a name to remove
(``cluster.remove``). Output: plain dicts for the service's results; the
cluster's entry in ctun's config (through ctun's own API, which keeps the
file's comments and a ``.bak``) and in ``remote.clusters`` of DELFIN's
settings.

Semantics:

* Only DELFIN's own process writes these entries, and only for a request
  from a window (the supervisor's ACTING gate): they decide where DELFIN's
  ssh goes and with which options. An agent can at most propose a cluster
  (``propose``): a NEEDS YOU item that opens this form filled in, for the
  person to check and send.
* ``add`` writes ctun's entry first (``api.add_cluster``; a preset brings
  its login prompts and its budget guard, with ``guard_scope:
  every-command``; a cluster given by host has no budget guard), then
  DELFIN's. When DELFIN's cannot be written, ctun's is taken out again, so a
  failed add leaves neither. A name either config has already is refused.
* The host key: a preset with a fingerprint its operators publish is pinned
  at once when the server presents exactly that key (``pin_host_key``
  refuses any other). Otherwise -- no published key, or a server that
  presents a different one -- the key it presents is returned unpinned, and
  only ``pin`` with that fingerprint, which the person confirmed, pins it.
  Until then DELFIN refuses the cluster (``clusters.load``).
* ``remove`` takes the cluster out of DELFIN's settings and out of ctun's
  config (which closes a live tunnel and unpins its key). It refuses while
  DELFIN still follows a job of the cluster's: its results could not come
  back.
* Every value the person typed is checked here as ``clusters.load`` would
  check it, before anything is written; an error names the field.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from . import clusters as registry

#: The fields of DELFIN's own entry the form sets.
_DIR_KEYS = ("work_dir", "install_root")


def _ctun():
    from cluster_tunnel import api as ctun_api
    return ctun_api


def _text(value: Any) -> str:
    return value.strip() if isinstance(value, str) else ""


def presets() -> dict:
    """``{"presets": [...]}``: the clusters ctun knows, for the form's picker."""
    from delfin.agent.service.clusters_feed import clean
    try:
        found = _ctun().presets()
    except Exception as exc:  # noqa: BLE001 - an older ctun: the form offers "other" only
        return {"presets": [], "error": clean(f"ctun lists no presets: {exc}", 200)}
    out = []
    for p in found:
        budget = getattr(p, "budget", None)
        out.append({
            "name": clean(p.name, 64), "title": clean(p.title, 120),
            "host": clean(p.host or "", 120) or None,
            "requiresOtp": bool(p.requires_otp), "requiresPassword": bool(p.requires_password),
            "userHint": clean(p.user_hint or "", 160),
            "budget": None if budget is None else {
                "limit": budget.session_limit, "unit": clean(budget.unit, 16)},
            "fingerprint": clean(p.fingerprint or "", 120) or None,
            "fingerprintSource": clean(p.fingerprint_source or "", 200) or None,
            "partition": clean(PRESET_PARTITIONS.get(p.name, ""), 64) or None,
        })
    return {"presets": out, "error": None}


#: The partition DELFIN's jobs default to on a preset's cluster, where one
#: suits DELFIN's default of a few cores for hours (the form's pre-filled
#: value; the person may change it).
PRESET_PARTITIONS = {"bwuni3": "cpu", "kcist": "small", "horeka": "cpuonly"}


def _checked_entry(params: dict) -> tuple[dict, list[str]]:
    """DELFIN's entry from the form, and what is wrong with it."""
    problems: list[str] = []
    entry: dict[str, Any] = {}
    partition = _text(params.get("partition"))
    if partition:
        if not registry._PARTITION.fullmatch(partition):
            problems.append(f"partition {partition!r} is not a name Slurm would accept")
        else:
            entry["default_partition"] = partition
    for key, field in (("work_dir", "workDir"), ("install_root", "installRoot")):
        value = _text(params.get(field))
        if not value:
            continue
        before = len(problems)
        checked = registry._remote_dir(value, key, problems)
        if len(problems) == before and checked:
            entry[key] = checked
    return entry, problems


def add(params: dict, *, settings_path: Optional[str] = None,
        ctun_config: Optional[str] = None) -> dict:
    """``cluster.add``: ctun's entry and DELFIN's, then the host key.

    ``{"ok", "code", "error", "name", "hostKey"}``; ``hostKey`` is
    ``{"pinned", "fingerprint", "source", "offered", "expected", "error"}``.
    """
    from delfin.agent.service.clusters_feed import clean, forget
    from delfin.user_settings import load_settings, save_settings

    ctun = _ctun()
    name = _text(params.get("name"))
    user = _text(params.get("user"))
    preset = _text(params.get("preset")) or None
    host = _text(params.get("host")) or None
    # Whether the login asks for a one-time code after the password: the
    # person says it for another cluster; a preset (IntNano included) knows.
    otp = params.get("requiresOtp")
    requires_otp = otp if isinstance(otp, bool) and not preset else None

    def refused(code: str, why: str) -> dict:
        return {"ok": False, "code": code, "error": clean(why, 300), "name": name, "hostKey": None}

    if not name:
        return refused("invalid", "the cluster needs a name")
    if not user:
        return refused("invalid", "the user name on the cluster is empty")
    if not preset and not host:
        return refused("invalid", "pick a known cluster, or give the login host of another")
    entry, problems = _checked_entry(params)
    if problems:
        return refused("invalid", "; ".join(problems))
    try:
        settings = load_settings(settings_path)
    except Exception as exc:  # noqa: BLE001
        return refused("settings_unreadable", f"DELFIN's settings cannot be read: {exc}")
    remote = settings.get("remote") if isinstance(settings.get("remote"), dict) else {}
    known = remote.get("clusters") if isinstance(remote.get("clusters"), dict) else {}
    if name in known:
        return refused("exists", f"DELFIN has a cluster {name!r} already")
    try:
        added = ctun.add_cluster(name, user=user, preset=preset, host=host,
                                 requires_otp=requires_otp, config_path=ctun_config)
    except ctun.ClusterExistsError:
        return refused("exists", f"ctun has a cluster {name!r} already; pick another name, "
                                 f"or remove it there (`ctun remove {name}`)")
    except ctun.ConfigError as exc:
        return refused("invalid", str(exc))
    except Exception as exc:  # noqa: BLE001
        return refused("ctun_failed", f"ctun did not add it: {type(exc).__name__}: {exc}")
    try:
        new_remote = dict(remote)
        new_clusters = dict(known)
        new_clusters[name] = entry
        new_remote["clusters"] = new_clusters
        save_settings({**settings, "remote": new_remote}, settings_path)
    except Exception as exc:  # noqa: BLE001 - ctun's entry goes again: neither stays
        try:
            ctun.remove_cluster(name, config_path=ctun_config)
        except Exception:  # noqa: BLE001
            pass
        return refused("settings_unwritable", f"DELFIN's settings could not be written: {exc}")
    finally:
        forget(name)
    from . import attention as remote_attention
    remote_attention.added(name)
    key = _first_pin(ctun, name, added.fingerprint, added.fingerprint_source, ctun_config)
    return {"ok": True, "code": None, "error": None, "name": name, "hostKey": key}


def _first_pin(ctun, name: str, fingerprint: Optional[str], source: Optional[str],
               ctun_config: Optional[str]) -> dict:
    """The published key pinned, if the server presents it; else what it presents."""
    from delfin.agent.service.clusters_feed import clean
    out = {"pinned": False, "fingerprint": None, "source": clean(source or "", 200) or None,
           "offered": None, "expected": clean(fingerprint or "", 120) or None, "error": None}
    if fingerprint:
        try:
            ctun.pin_host_key(name, fingerprint, config_path=ctun_config)
            out.update(pinned=True, fingerprint=clean(fingerprint, 120))
            return out
        except ctun.HostKeyMismatch as exc:
            out["offered"] = clean(getattr(exc, "offered", "") or "", 120) or None
            return out
        except Exception as exc:  # noqa: BLE001 - not reached: the window may ask again
            out["error"] = clean(f"the server could not be asked for its key: {exc}", 300)
            return out
    try:
        status = ctun.host_key(name, config_path=ctun_config)
        out["offered"] = clean(status.offered.fingerprint, 120)
    except Exception as exc:  # noqa: BLE001
        out["error"] = clean(f"the server could not be asked for its key: {exc}", 300)
    return out


def host_key(name: str, *, ctun_config: Optional[str] = None) -> dict:
    """``cluster.hostkey``: the key the server presents now, for the person
    to compare (``{"offered", "pinned", "error"}``)."""
    from delfin.agent.service.clusters_feed import clean
    ctun = _ctun()
    try:
        status = ctun.host_key(name, config_path=ctun_config)
    except Exception as exc:  # noqa: BLE001
        return {"offered": None, "pinned": False,
                "error": clean(f"the server could not be asked for its key: {exc}", 300)}
    return {"offered": clean(status.offered.fingerprint, 120),
            "pinned": bool(status.strict and status.matches),
            "error": None}


def pin(name: str, fingerprint: str, *, ctun_config: Optional[str] = None) -> dict:
    """``cluster.pin``: the key the person confirmed, pinned -- only if the
    server still presents exactly that one."""
    from delfin.agent.service.clusters_feed import clean, forget
    ctun = _ctun()
    fingerprint = _text(fingerprint)
    if not fingerprint:
        return {"ok": False, "code": "invalid", "error": "no fingerprint was confirmed"}
    try:
        ctun.pin_host_key(name, fingerprint, config_path=ctun_config)
    except ctun.HostKeyMismatch as exc:
        offered = clean(getattr(exc, "offered", "") or "", 120)
        return {"ok": False, "code": "host_key_mismatch",
                "error": f"the server now presents {offered or 'another key'}, not the key you "
                         f"confirmed: nothing was pinned"}
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "code": "not_pinned", "error": clean(str(exc), 300)}
    finally:
        forget(name)
    return {"ok": True, "code": None, "error": None}


def _followed(name: str) -> int:
    """How many jobs of the cluster's DELFIN still follows."""
    try:
        from . import jobs, watch
        return sum(1 for record in jobs.records(cluster=name)
                   if watch.followed(name, record.get("job_id")))
    except Exception:  # noqa: BLE001 - no record to read: nothing is followed
        return 0


def remove(name: str, *, settings_path: Optional[str] = None,
           ctun_config: Optional[str] = None) -> dict:
    """``cluster.remove``: out of DELFIN's settings and ctun's config."""
    from delfin.agent.service.clusters_feed import clean, forget
    from delfin.user_settings import load_settings, save_settings

    ctun = _ctun()
    name = _text(name)
    following = _followed(name)
    if following:
        return {"ok": False, "code": "jobs_followed",
                "error": f"DELFIN still follows {following} job{'s' if following != 1 else ''} "
                         f"on {name}: its results could not come back. Remove it once they have."}
    try:
        settings = load_settings(settings_path)
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "code": "settings_unreadable", "error": clean(str(exc), 300)}
    remote = settings.get("remote") if isinstance(settings.get("remote"), dict) else {}
    known = remote.get("clusters") if isinstance(remote.get("clusters"), dict) else {}
    in_delfin = name in known
    try:
        in_ctun = name in ctun.clusters(config_path=ctun_config)
    except Exception:  # noqa: BLE001 - an unreadable config: DELFIN's entry still goes
        in_ctun = False
    if not in_delfin and not in_ctun:
        return {"ok": False, "code": "unknown", "error": f"no cluster {name!r} to remove"}
    try:
        if in_delfin:
            new_clusters = {k: v for k, v in known.items() if k != name}
            save_settings({**settings, "remote": {**remote, "clusters": new_clusters}},
                          settings_path)
        removed = ctun.remove_cluster(name, config_path=ctun_config) if in_ctun else None
    except Exception as exc:  # noqa: BLE001
        return {"ok": False, "code": "not_removed", "error": clean(str(exc), 300)}
    finally:
        forget(name)
    return {"ok": True, "code": None, "error": None,
            "loggedOut": bool(getattr(removed, "logged_out", False)),
            "keysRemoved": int(getattr(removed, "keys_removed", 0) or 0)}


def propose(*, name: Any, preset: Any = None, user: Any = None, host: Any = None,
            partition: Any = None, session_id: str = "", workspace: Any = None) -> dict:
    """An agent's proposal: a NEEDS YOU item that opens the add form filled
    in. Nothing is written but the item; the person checks and sends it."""
    from delfin.agent.service.clusters_feed import clean
    proposal = {"name": clean(_text(name), 64), "preset": clean(_text(preset), 64),
                "user": clean(_text(user), 64), "host": clean(_text(host), 120),
                "partition": clean(_text(partition), 64)}
    if not proposal["name"]:
        return {"ok": False, "error": "a proposal needs the cluster's name"}
    if not proposal["preset"] and not proposal["host"]:
        return {"ok": False, "error": "name a preset (cluster_propose lists none: ask the user) "
                                      "or the login host"}
    try:
        known = {p.name: p for p in _ctun().presets()}
    except Exception:  # noqa: BLE001
        known = {}
    if proposal["preset"] and known and proposal["preset"] not in known:
        return {"ok": False, "error": f"unknown preset {proposal['preset']!r}; known: "
                                      f"{', '.join(sorted(known))}"}
    from . import attention as remote_attention
    item = remote_attention.propose_cluster(proposal, session_id=session_id,
                                            workspace=workspace)
    if item is None:
        return {"ok": False, "error": "the attention inbox could not be written"}
    return {"ok": True, "item": item, "presets": sorted(known)}


#: The script that says where the installer would find ORCA (read-only).
ORCA_SCRIPT = Path(__file__).with_name("scripts") / "orca.sh"


def orca_status(conn) -> dict:
    """Where DELFIN's installer would find ORCA on ``conn``'s cluster:
    ``{"found", "path", "how", "hint"}``. ORCA is the user's own (its
    licence): DELFIN looks for it and says what to do, never copies it."""
    from delfin.agent.service.clusters_feed import clean
    result = conn.feed_script(ORCA_SCRIPT, [conn.cluster.software_dir], timeout=60)
    words = (result.stdout or "").strip().split(None, 1) if result.ok else []
    how = words[0] if words else ""
    path = clean(words[1], 200) if len(words) > 1 else ""
    if how == "dir":
        return {"found": True, "path": path, "how": "dir", "hint": ""}
    if how == "tarball":
        return {"found": True, "path": path, "how": "tarball",
                "hint": "the installer unpacks it at the first submit"}
    return {"found": False, "path": "", "how": "none",
            "hint": "put your ORCA into the cluster's home (an unpacked orca_* directory, or "
                    "its .tar.xz from orcaforum.kofo.mpg.de): DELFIN never copies ORCA"}
