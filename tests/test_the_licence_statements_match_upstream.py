"""What DELFIN says about other programs' licences matches their upstream.

Input: ``README.md`` and ``delfin/cli_banner.py``. Output: assertions.

The upstream licence files, read on 2026-10-15 through the GitHub API:
grimme-lab/xtb, crest-lab/crest and grimme-lab/xtb4stda carry ``COPYING``
(GPL-3.0) plus ``COPYING.LESSER`` (LGPL-3.0) and state "version 3 of the
License, or (at your option) any later version" in their READMEs, i.e.
LGPL-3.0-or-later. The README called xTB and CREST "free for academic use
under the GNU General Public License (GPLv3)" and the banner called CREST
GPL v3; the LGPL restricts neither academic nor commercial use.
"""

from __future__ import annotations

import re
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
README = (REPO / "README.md").read_text(encoding="utf-8")
BANNER = (REPO / "delfin" / "cli_banner.py").read_text(encoding="utf-8")


def _section(text: str, heading: str) -> str:
    start = text.index(heading)
    nxt = re.search(r"^#{1,6} ", text[start + len(heading):], re.M)
    return text[start:start + len(heading) + (nxt.start() if nxt else len(text))]


def test_xtb_and_crest_are_named_lgpl_in_the_readme():
    for heading in ("## xTB Requirements", "## CREST Requirements"):
        body = _section(README, heading)
        assert "LGPL-3.0-or-later" in body, heading
        assert "General Public License (GPLv3)" not in body, heading
        assert "free for academic use" not in body, heading


def test_the_banner_names_crest_and_xtb_lgpl():
    assert "GNU General Public License (GPL v3)" not in BANNER
    for opening in ("xTB 6.7.1+ is distributed", "CREST 3.0.2 is released"):
        paragraph = BANNER[BANNER.index(opening):][:200]
        assert "LGPL-3.0-or-later" in paragraph, opening


def test_the_readme_says_where_the_lgpl_text_is():
    body = _section(README, "## License")
    assert "COPYING.LESSER" in body and "COPYING" in body
    assert "in the file [LICENSE](LICENSE)" not in body


def test_the_readme_points_at_the_notices_file():
    assert "THIRD_PARTY_NOTICES.md" in _section(
        README, "## Other programs the installers download")
