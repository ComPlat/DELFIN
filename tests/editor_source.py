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

from delfin.dashboard import molecule_viewer as _molecule_viewer
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

#: The one fullscreen: the sheet every molecule overlay is laid out by, and
#: the script that moves the members into it and back. The Submit tab used to
#: describe an overlay of its own here, beside a second implementation for the
#: other three tabs, so a test about what fullscreen does had to know which tab
#: it was asking about. It does not any more -- these are the answer for all of
#: them, and a test that reads them is testing every tab at once.
FULLSCREEN_CSS = _molecule_viewer.STRUCTURE_VIEWER_FULLSCREEN_CSS
FULLSCREEN_JS = _molecule_viewer.STRUCTURE_VIEWER_FULLSCREEN_BOOTSTRAP_JS
