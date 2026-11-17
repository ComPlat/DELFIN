"""The browser tab carries the DELFIN logo, not the blank icon Voila ships.

Voila's own page (``voila/templates/lab/index.html.j2``, line 60 in 0.5.13)
contains::

    <link rel="icon" href="data:;base64,=">

a blank icon it ships so the browser asks for none. It sits AFTER
``</head>``, which is later in the document than anything appended to the
head, and with two icon links a browser takes the last one. The dashboard's
own favicon was therefore added and then ignored, and the tab showed no
logo. The pin for ``voila`` allows anything below 0.6, so this can come
back with an upgrade; the fix does not depend on where Voila puts its link.

Decided in Python, as a browser fix must be: the emitted script is parsed
for the two decisions that matter -- every foreign icon link is dropped,
and the dashboard's own link is the one that gets the logo. The DOM
behaviour itself was checked once by driving the emitted script against a
page that carries Voila's link, and the old script left both links in
place while this one leaves only the logo.
"""

from __future__ import annotations

from delfin.dashboard import helpers as H

LOGO = "data:image/png;base64,LOGOBYTES"


def _script(logo: str = LOGO) -> str:
    sent: list[str] = []
    original = H._append_js
    H._append_js = lambda ctx, script: sent.append(script)
    try:
        H.apply_branding(object(), title="DELFIN Dashboard",
                         favicon_data_uri=logo)
    finally:
        H._append_js = original
    assert sent, "apply_branding sent nothing"
    return sent[0]


def test_every_foreign_icon_link_is_dropped():
    script = _script()
    assert "link[rel~='icon']" in script, (
        "the script must find EVERY icon link, not only its own: Voila's "
        "blank one carries no marker of ours")
    assert "removeChild" in script, (
        "a foreign icon link left in the page wins, because it sits after "
        "</head> and the browser takes the last one")


def test_the_dashboards_own_link_carries_the_logo():
    script = _script()
    assert "data-delfin-favicon" in script
    assert LOGO in script, "the logo never reaches the link"


def test_without_a_logo_nothing_is_touched():
    """No logo in the installation: the page keeps whatever it had, rather
    than losing its icon to an empty href."""
    script = _script(logo="")
    assert "if (!'') return;" in script.replace('"', "'"), script[:200]


def test_the_title_is_set_and_held():
    script = _script()
    assert "document.title = 'DELFIN Dashboard'" in script
    assert "MutationObserver" in script, (
        "the lab front-end rewrites the title after the page is up")
