"""Connecting a phone takes a button, a scan and a six-digit code.

Before: the QR code needed qrencode on the host (missing: "apt install
qrencode"), the dashboard showed the invitation as preformatted text, the
code was typed as "/pair 123456" into the chat, and a wrong code dropped
the panel.

Output, each with its control:

* delfin/messenger/qr.py draws the code with no tool installed. The matrix
  of two inputs is pinned (both decoded by an independent reader, OpenCV's
  QRCodeDetectorAruco, when the pin was taken); every code has its three
  finder patterns, its timing rows and a valid format word;
* the invitation text carries that code with qrencode absent;
* the dashboard panel shows the steps, the code as a picture (SVG) and the
  link; after the phone connects it asks for the code (stage "code"); a
  wrong code keeps the panel with the attempts left; a verified phone
  clears it and says so in the chat; the link never reaches the chat;
* the first message on the phone names what to try first.
"""
from __future__ import annotations

import hashlib

import pytest

from delfin.messenger import qr

PINNED = {
    "simplex:/invitation#/?v=2-7&smp=smp%3A%2F%2Fabc%40example.org%2Fq%23key":
        (33, "c51f6e35252237f994f531d5bb2369b5840a2fb5b2a05c73dba876af89956b00"),
    "DELFIN": (21, "0d52e9ecfb1fcf028734ef27bc3c165023a00c6f02c699e3d61d2864010f8076"),
}
LINK = "simplex:/invitation#/?v=2-7&smp=smp%3A%2F%2Fabc%40example.org%2Fq%23key"


@pytest.mark.parametrize("text", sorted(PINNED))
def test_the_code_of_a_known_input_is_the_pinned_one(text):
    m = qr.encode(text)
    bits = "".join("1" if v else "0" for row in m for v in row)
    assert (len(m), hashlib.sha256(bits.encode()).hexdigest()) == PINNED[text]


def _finder_at(m, r0, c0):
    for r in range(7):
        for c in range(7):
            ring = max(abs(r - 3), abs(c - 3))
            if m[r0 + r][c0 + c] != (ring != 2):
                return False
    return True


@pytest.mark.parametrize("n", [1, 40, 200, 600, 1200])
def test_every_code_has_finders_timing_and_a_valid_format_word(n):
    text = ("simplex:/x" + "ab9" * n)[:n + 10]
    m = qr.encode(text)
    size = len(m)
    assert (size - 17) % 4 == 0
    assert _finder_at(m, 0, 0) and _finder_at(m, 0, size - 7) and _finder_at(m, size - 7, 0)
    assert [m[6][c] for c in range(8, size - 8)] == [c % 2 == 0 for c in range(8, size - 8)]
    assert [m[r][6] for r in range(8, size - 8)] == [r % 2 == 0 for r in range(8, size - 8)]
    # Format word around the top-left finder (bits 14..0), XOR mask 0x5412:
    # a valid BCH(15,5) codeword.
    seq = [m[8][c] for c in (0, 1, 2, 3, 4, 5, 7, 8)] + [m[r][8] for r in (7, 5, 4, 3, 2, 1, 0)]
    word = int("".join("1" if b else "0" for b in seq), 2) ^ 0x5412
    rem = word
    for i in range(14, 9, -1):
        if rem >> i & 1:
            rem ^= 0x537 << (i - 10)
    assert rem == 0 and m[size - 8][8] is True        # and the dark module


def test_the_invitation_carries_a_code_without_qrencode(monkeypatch):
    import shutil
    from delfin.messenger import pairing
    monkeypatch.setattr(shutil, "which", lambda name: None)
    text = pairing.offer_text({"link": LINK, "expires": 9e12})
    assert "qrencode" not in text
    assert "█" in text or "▀" in text
    assert LINK in text


class _Panel:
    def __init__(self):
        from delfin.messenger.dashboard import PairPanel
        self.html, self.said, self.stages = [], [], []
        self.panel = PairPanel(render=self.html.append, say=self.said.append,
                               on_stage=self.stages.append, ttl_s=0)


def test_the_panel_shows_steps_a_picture_and_the_link():
    from delfin.messenger import pairing
    p = _Panel()
    p.panel.show(pairing.offer_text({"link": LINK, "expires": 9e12}))
    html = p.html[-1]
    assert "<svg" in html and "Scan" in html
    assert LINK.replace("&", "&amp;") in html
    assert "█" not in html and "▀" not in html, "the text code is not shown beside the picture"
    assert p.stages[-1] == "offer" and not p.said


def test_after_the_phone_connects_the_panel_asks_for_the_code():
    from delfin.messenger import pairing
    p = _Panel()
    p.panel.show(pairing.offer_text({"link": LINK, "expires": 9e12}))
    p.panel.show(pairing.state_text({"state": "connected", "label": "Pixel"}))
    assert p.stages[-1] == "code"
    assert "6-digit" in p.html[-1] and LINK not in p.html[-1]


def test_a_wrong_code_keeps_the_panel_with_the_attempts_left():
    from delfin.messenger import pairing
    p = _Panel()
    p.panel.show(pairing.state_text({"state": "connected", "label": "Pixel"}))
    p.panel.show(pairing.state_text({"state": "wrong", "left": 2}))
    assert p.stages[-1] == "code" and "2 attempt" in p.html[-1]
    assert not p.said


def test_a_verified_phone_clears_the_panel_and_says_so_in_the_chat():
    from delfin.messenger import pairing
    p = _Panel()
    p.panel.show(pairing.offer_text({"link": LINK, "expires": 9e12}))
    p.panel.show(pairing.state_text({"state": "verified", "label": "Pixel"}))
    assert p.html[-1] == "" and p.stages[-1] is None
    assert p.said and "verified" in p.said[-1]
    assert all("simplex:/" not in s for s in p.said)


def test_the_first_message_on_the_phone_names_what_to_try(tmp_path, monkeypatch):
    import asyncio
    from _messenger_people_rig import Account
    from delfin.messenger import store
    monkeypatch.delenv("DELFIN_MESSENGER_LANG", raising=False)
    acct = Account(tmp_path, monkeypatch)
    store.save("21", person="anna", label="anna phone")
    asyncio.run(acct.rooms._welcome("21", "anna", {"label": "anna phone"}))
    first = acct.daemon.to("@21")[0]
    for word in ("sessions", "get", "new", "help"):
        assert word in first, first
