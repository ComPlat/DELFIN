"""Which dashboard tabs show, and in which order, from the settings.

Input: the tab specs (``id``, ``default_order``, ``default_visible``,
``fixed``) and the person's ``ui.tabs`` settings -- ``order``, the
arrangement Settings > Dashboard Tabs saves (every tab of its day), and
``hidden``, the tabs switched off there.

Output: :func:`resolve_tab_preferences` -> ``(order, hidden)`` over the
movable tabs; the tab bar shows ``order`` minus ``hidden``, then the fixed
ones (Settings).

The first start decides once (:func:`first_start`, marker
``ui.first_start``): a new installation -- no saved arrangement and no
trace of earlier use -- saves the five :data:`FIRST_TABS` as its
arrangement; anyone who used DELFIN before saves the view they had (the
defaults of before 2026-11-13, kept on each spec as ``legacy_visible`` and
``legacy_order``), so nothing changes for them. From then on the saved
arrangement decides, and the defaults (:func:`apply_first_tabs`) only
decide a tab added after it was saved.

Traces of earlier use, chosen because only use writes them: a saved agent
session (the session store), a job in the local job list
(``~/.delfin_jobs.json``; SLURM jobs live in the calculations folder), an
entry in the calculations or archive folder. Rejected as a trace: the
settings file (the installer writes it before the dashboard ever starts),
the dashboard's own note ``dashboard_here.json`` (the launcher writes it
and removes it at exit), an empty job list (written at the first start).
Rejected: deciding by the defaults alone without a marker (an existing
person who never arranged the tabs would have lost eleven of them).

"""
from __future__ import annotations

import json
import time
from pathlib import Path

#: The tabs a person sees before arranging any, in this order (2026-11-13).
#: Every other tab stays available: Settings > Dashboard Tabs switches it on.
FIRST_TABS = ("submit", "job_status", "agent", "home", "calculations")


def apply_first_tabs(specs: list) -> list:
    """Set ``default_visible`` and the leading ``default_order`` of every
    movable spec from FIRST_TABS (in place; returns ``specs``). Fixed tabs
    (Settings) are left as they are: always shown."""
    for spec in specs:
        if spec.get("fixed"):
            continue
        tab_id = spec.get("id")
        spec.setdefault("legacy_visible", bool(spec.get("default_visible", True)))
        spec.setdefault("legacy_order", spec.get("default_order", 10_000))
        spec["default_visible"] = tab_id in FIRST_TABS
        if tab_id in FIRST_TABS:
            spec["default_order"] = FIRST_TABS.index(tab_id) + 1
    return specs


def resolve_tab_preferences(specs: list, prefs: dict | None) -> tuple[list, set]:
    """``(order, hidden)`` of the movable tabs for ``prefs`` (``ui.tabs``).

    The saved order first (unknown and fixed ids dropped), then every other
    tab by ``default_order``. A tab the saved order names is hidden when the
    saved hidden list says so; any other tab when its default says so."""
    prefs = prefs or {}
    spec_map = {spec.get("id"): spec for spec in specs or []}
    raw_order = [str(item) for item in (prefs.get("order", []) or []) if str(item).strip()]
    raw_hidden = {str(item) for item in (prefs.get("hidden", []) or []) if str(item).strip()}
    order = []
    seen = set()
    for tab_id in raw_order:
        spec = spec_map.get(tab_id)
        if not spec or spec.get("fixed") or tab_id in seen:
            continue
        seen.add(tab_id)
        order.append(tab_id)
    for spec in sorted(specs or [], key=lambda item: item.get("default_order", 10_000)):
        tab_id = spec.get("id")
        if not tab_id or spec.get("fixed") or tab_id in seen:
            continue
        seen.add(tab_id)
        order.append(tab_id)
    hidden = set()
    saved_ids = set(raw_order)
    for tab_id in order:
        spec = spec_map.get(tab_id) or {}
        if tab_id in saved_ids:
            if tab_id in raw_hidden:
                hidden.add(tab_id)
        elif not bool(spec.get("default_visible", True)):
            hidden.add(tab_id)
    return order, hidden


#: The settings key of the first-start marker (under ``ui``).
FIRST_START_KEY = "first_start"


def _has_entry(folder) -> bool:
    """A visible entry (not a dot file) in ``folder``."""
    try:
        return any(not p.name.startswith(".") for p in Path(folder).iterdir())
    except (OSError, TypeError):
        return False


def _has_a_job(jobs_file) -> bool:
    try:
        data = json.loads(Path(jobs_file).read_text(encoding="utf-8"))
    except (OSError, ValueError, TypeError):
        return False
    return bool(isinstance(data, dict) and data.get("jobs"))


def traces_of_use(*, calc_dir=None, archive_dir=None, sessions_dir=None,
                  jobs_file=None) -> list:
    """The traces of earlier use found, as names (no paths)."""
    found = []
    if sessions_dir is not None and _has_entry(sessions_dir):
        found.append("agent session")
    if jobs_file is not None and _has_a_job(jobs_file):
        found.append("job in the local job list")
    if calc_dir is not None and _has_entry(calc_dir):
        found.append("calc folder entry")
    if archive_dir is not None and _has_entry(archive_dir):
        found.append("archive folder entry")
    return found


def _arrangement(specs, *, visible, order_key):
    movable = sorted((s for s in specs if not s.get("fixed") and s.get("id")),
                     key=lambda s: s.get(order_key, 10_000))
    return {"order": [s["id"] for s in movable],
            "hidden": [s["id"] for s in movable if not visible(s)]}


def first_start(specs, *, settings_path=None, calc_dir=None, archive_dir=None,
                sessions_dir=None, jobs_file=None):
    """Decide the tab arrangement once, at the dashboard's first start.

    Returns what was decided -- ``"new"``, ``"existing"``, ``"arranged"`` --
    or None when the marker was already there (or the settings cannot be
    read). Writes ``ui.tabs`` (new, existing) and ``ui.first_start``
    through user_settings.save_settings; a failed write leaves the
    defaults in force for this start and tries again at the next.
    """
    from delfin import user_settings
    try:
        settings = user_settings.load_settings(settings_path)
    except Exception:
        return None
    ui = dict(settings.get("ui") or {})
    if ui.get(FIRST_START_KEY):
        return None
    tabs = ui.get("tabs") or {}
    marker = {"at": time.strftime("%Y-%m-%dT%H:%M:%S")}
    if tabs.get("order"):
        kind = "arranged"
    else:
        traces = traces_of_use(calc_dir=calc_dir, archive_dir=archive_dir,
                               sessions_dir=sessions_dir, jobs_file=jobs_file)
        if traces:
            kind = "existing"
            marker["traces"] = traces
            ui["tabs"] = _arrangement(
                specs, visible=lambda s: s.get("legacy_visible", True),
                order_key="legacy_order")
        else:
            kind = "new"
            ui["tabs"] = _arrangement(
                specs, visible=lambda s: s.get("id") in FIRST_TABS,
                order_key="default_order")
    marker["tabs"] = kind
    ui[FIRST_START_KEY] = marker
    settings["ui"] = ui
    try:
        user_settings.save_settings(settings, settings_path)
    except Exception:
        return None
    return kind
