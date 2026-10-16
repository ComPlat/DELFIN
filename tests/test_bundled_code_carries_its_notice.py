"""Code by other authors in the tree ships with the notice its licence requires.

Input: ``THIRD_PARTY_NOTICES.md``, ``delfin/dashboard/static/3Dmol-min.js.LICENSE.txt``
and the tracked file list. Output: assertions.

BSD-3-Clause clause 2 and the MIT licence both require the copyright notice
and the licence text in every redistributed copy. Until this file existed
the bundle's companion LICENSE.txt held only the 120-byte webpack banner,
and the OpenTUI patches carried no notice at all.
"""

from __future__ import annotations

import subprocess
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NOTICES = REPO / "THIRD_PARTY_NOTICES.md"
BUNDLE_LICENSE = REPO / "delfin" / "dashboard" / "static" / "3Dmol-min.js.LICENSE.txt"


def _tracked() -> list[str]:
    try:
        out = subprocess.run(["git", "ls-files"], cwd=REPO, capture_output=True,
                             text=True, check=True).stdout
    except (OSError, subprocess.CalledProcessError):
        # An unpacked sdist has no git; its file set is the tree itself.
        return [str(p.relative_to(REPO)) for p in REPO.rglob("*") if p.is_file()]
    return out.splitlines()


def test_the_bundle_licence_file_carries_the_full_3dmol_licence():
    text = BUNDLE_LICENSE.read_text(encoding="utf-8")
    assert text.startswith("/*!\n * 3dmol v2.5.3"), "the webpack banner stays first"
    for clause in (
        "Copyright (c) 2014, University of Pittsburgh and contributors",
        "Redistributions in binary form must reproduce the above copyright",
        "Neither the name of the copyright holder nor the names of its",
        'THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS',
    ):
        assert clause in text, clause


def test_every_component_inside_the_bundle_has_its_text():
    text = BUNDLE_LICENSE.read_text(encoding="utf-8")
    for holder in (
        "Vitaly Puzrin and Andrei Tuputcyn",      # pako, MIT
        "Jean-loup Gailly and Mark Adler",        # pako's zlib port, Zlib
        "Copyright (c) 2017 Photopea",            # UPNG.js
        "Michaël Zasso",                          # iobuffer
        "Copyright (c) 2016 cheminfo",            # netcdfjs
        "Copyright (c) 2017 Simon Rodriguez",     # Rendu FXAA shader
        "Mol* contributors",                      # Mol* port of it
    ):
        assert holder in text, holder
    assert text.count("Permission is hereby granted, free of charge") >= 6


def test_the_notices_file_and_the_bundle_file_quote_one_3dmol_text():
    """Two copies of one text drift; this pins them together."""
    bundle = BUNDLE_LICENSE.read_text(encoding="utf-8")
    start = bundle.index("3Dmol.js incorporates code from GLmol")
    end = bundle.index("=" * 72, start)
    quoted = bundle[start:end].rstrip("\n")
    assert quoted in NOTICES.read_text(encoding="utf-8")


def test_the_opentui_patches_carry_the_mit_notice():
    text = NOTICES.read_text(encoding="utf-8")
    assert "Copyright (c) 2025 opentui" in text
    assert "Permission is hereby granted, free of charge" in text
    for patch in sorted((REPO / "tui" / "patches").glob("*.patch")):
        assert patch.name in text, f"{patch.name} is not named in THIRD_PARTY_NOTICES.md"


def test_every_vendored_or_minified_file_is_named_in_the_notices():
    """A new vendored file without a notice is the defect this guards."""
    text = NOTICES.read_text(encoding="utf-8")
    vendored = [p for p in _tracked()
                if p.endswith((".min.js", "-min.js", ".min.css"))
                or "/vendor/" in p or "/third_party/" in p]
    assert vendored, "the scan found nothing; 3Dmol-min.js should be there"
    for path in vendored:
        assert path in text, f"{path} is not named in THIRD_PARTY_NOTICES.md"


def test_the_notices_say_the_downloaded_programs_are_not_distributed():
    text = NOTICES.read_text(encoding="utf-8")
    for program, licence in (("xtb4stda", "LGPL-3.0-or-later"),
                             ("CREST", "LGPL-3.0-or-later"),
                             ("g-xTB", "GPL-3.0"),
                             ("Bun", "LGPL-2.0")):
        row = next((ln for ln in text.splitlines()
                    if ln.startswith(f"| {program}")), "")
        assert licence in row, (program, row)
