"""The Ketcher tab: the editor at the size it was built to be used at.

Registered additively, so nothing built in changes and a tab that fails to
build cannot take the dashboard down with it -- the registry catches that and
marks it unavailable.  Named in :data:`tab_registry._BUILTIN_DYNAMIC_TABS`, or
the module is never imported and the registration never runs.

The panel itself is :mod:`delfin.dashboard.ketcher_panel`; this is only where
it is placed and how the rest of the dashboard finds it.  The Calculations
browser needs the second part: a drawing double-clicked over there is opened
here, and the registry throws a factory's refs away, so the panel publishes
itself onto the context rather than waiting to be handed on.
"""

from __future__ import annotations

from typing import Any, Dict

import ipywidgets as widgets

from . import ketcher_panel as _panel

__all__ = ['create_tab']


def create_tab(ctx) -> Any:
    """Build the Ketcher tab.  Returns ``(tab_widget, refs_dict)``."""
    panel = _panel.build(ctx, height='72vh', scope='delfin-ketcher-tab')
    refs: Dict[str, Any] = {
        'ketcher_panel': panel,
        'open_drawing': panel.open_text,
        'refresh_files': panel._refresh_files,
    }
    # The registry keeps only the widget (``result[0]``), so anything the rest
    # of the dashboard has to reach has to be put somewhere it can be reached.
    #
    # Filled in rather than replaced. The Archive and Office tabs are the
    # Calculations browser over a different root, made with
    # ``dataclasses.replace`` -- a second context that shares this dict -- and
    # they are built before the registered tabs are. Rebinding the attribute
    # here would leave both of them holding the empty dict they were made
    # with, and a drawing opened from Archive would say there was no Ketcher
    # tab while it was sitting two tabs along.
    existing = getattr(ctx, 'ketcher_refs', None)
    if isinstance(existing, dict):
        existing.clear()
        existing.update(refs)
        refs = existing
    else:
        ctx.ketcher_refs = refs
    box = widgets.VBox([panel.widget],
                       layout=widgets.Layout(width='100%', padding='10px'))
    return box, refs


try:                                                 # pragma: no cover
    from delfin.dashboard.tab_registry import register_tab
    register_tab('ketcher', 'Ketcher', create_tab, order=45)
except Exception:                                    # pragma: no cover
    pass
