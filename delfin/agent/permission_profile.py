"""The permission profile a person chose, and the mode it starts.

Input: ``agent.permission_profile`` as the dashboard's Perms control
stores it in the user settings (``plan``, ``ask_all``, ``repo_free``,
``all_free``, or a name an older version stored). Output: the backend
``permission_mode`` of that profile (:data:`PROFILE_TO_MODE`), and for a
session start the stored profile with the settings file it came from
(:func:`stored_profile`).

One table for every door. The dashboard builds its engine from it
(``tab_agent.PROFILE_TO_CLI_PERM`` is this object), and the terminal doors
(``cli._headless_posture`` for ``run`` / ``chat -p``, and
``cli._start_interactive_session`` for ``chat``, the TUI, the agent
service's worker and the sessions the messenger opens through it) take it
as the start mode when neither ``--permission-mode`` nor a mode setting
(``agent.permission_mode``, ``kit.default_mode``) names one. Before, only
the dashboard read the profile: a host whose installer set ``all_free``
for every person ran its terminal and headless sessions in ``default``,
and every command off the auto-allow list was refused there.

What a profile does not change: the administrator's managed policy
coerces the mode wherever it is held (``KitToolPermissions.mode``), the
isolation start check still lowers bypass on a host without a layer
(``isolation_posture.apply_at_start``), and bypass drops the asking, not
the roots, the deny list or the write path gate.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Optional

#: The settings key the dashboard writes.
SETTING = "permission_profile"

#: Each rung of the dashboard's ladder -> the backend permission_mode.
PROFILE_TO_MODE: dict[str, str] = {
    "plan":      "plan",                # read-only
    "ask_all":   "default",             # asks before every write and shell
    "repo_free": "acceptEdits",         # writes go through, shell asks
    "all_free":  "bypassPermissions",   # asks nothing inside the folders
}

#: Names older versions stored under the same key (or as a mode).
LEGACY_NAMES: dict[str, str] = {
    "default": "ask_all", "erlaubt": "repo_free", "full": "all_free",
    "acceptEdits": "repo_free", "auto": "all_free",
    "bypassPermissions": "all_free",
}


def normalise(value) -> str:
    """The current profile name for ``value``, or "" when it names none."""
    name = str(value or "").strip()
    name = LEGACY_NAMES.get(name, name)
    return name if name in PROFILE_TO_MODE else ""


def mode_for(value) -> str:
    """The permission mode of profile ``value``, or "" when it names none."""
    return PROFILE_TO_MODE.get(normalise(value), "")


def stored_profile(layers: Optional[list] = None) -> tuple[str, str]:
    """``(profile, settings file)`` of the first settings layer that names
    a profile, or ``("", "")``.

    Read from the raw files in precedence order (personal first, then the
    account's defaults: ``user_settings.settings_layers``), as the other
    start-mode readers do: the merged view fills in values nobody chose.
    A value that is no profile is no setting. Never raises."""
    try:
        if layers is None:
            from delfin.user_settings import settings_layers
            layers = settings_layers()
    except Exception:
        return "", ""
    for path in layers or ():
        try:
            path = Path(path)
            if not path.exists():
                continue
            block = json.loads(path.read_text()).get("agent") or {}
        except Exception:
            continue
        name = normalise(block.get(SETTING)) if isinstance(block, dict) else ""
        if name:
            return name, str(path)
    return "", ""
