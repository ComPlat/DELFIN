"""Where the structure editor's source is, now that it is a part of its own.

The editor -- the molecule preview, its toolbar and everything that acts on the
structure -- used to be written out inside ``tab_submit.create_tab``. It is
``delfin/dashboard/structure_editor.py`` now, so that a second tab can hold the
same editor over its own coordinates.

Tests that read the source to pin down how the editor behaves want both files:
what the toolbar does is a fact about the part, and what the tab around it does
is a fact about the tab. Which of the two a given line landed in is not what
those tests are about, so they read the pair.
"""

from delfin.dashboard import structure_editor as _structure_editor
from delfin.dashboard import tab_submit as _tab_submit


def _read(module):
    return open(module.__file__, encoding='utf-8').read()


#: The structure editor on its own.
EDITOR_SOURCE = _read(_structure_editor)

#: The Submit tab on its own.
TAB_SOURCE = _read(_tab_submit)

#: The Submit tab and the editor it holds, as one text.
SUBMIT_SOURCE = TAB_SOURCE + '\n' + EDITOR_SOURCE
