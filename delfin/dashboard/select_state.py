"""New options for a selection widget, sent to the browser as one state.

ipywidgets' ``SelectMultiple`` keeps its selection as ``index``, a tuple of
positions into ``_options_labels``. Setting ``options`` in Python sends two
updates, one after the other: first the new ``_options_labels``, then
``index = ()``. The browser applies the first while the old ``index`` is
still set; when the old selection points past the end of the new list, the
view's ``updateSelection`` reads an option that is not there and the page
logs::

    Error setting state: Cannot set properties of undefined (setting 'selected')

(a file explorer whose listing shrinks under a selected file: entering a
folder by double click, Up, a filter that hides the selected file).

:func:`replace_options` holds the widget's sync while options and value
change, so the browser receives labels and selection in one message and
never holds a selection that points outside the list. Python sees exactly
the sequence it saw before -- options set (value reset to ``()``), then the
value -- so every ``observe`` handler runs as it did.
"""

from __future__ import annotations

from collections.abc import Iterable
from typing import Any


def replace_options(widget: Any, options: Iterable[Any], value: Any = ()) -> None:
    """Set ``widget.options`` and then ``widget.value``; one browser update.

    *value* defaults to no selection (``()``, what ``SelectMultiple`` falls
    back to whenever its options change). Pass ``None`` to leave the value
    to the widget's own reset (a single ``Select`` picks its first option).
    """
    with widget.hold_sync():
        widget.options = options
        if value is not None:
            widget.value = value
