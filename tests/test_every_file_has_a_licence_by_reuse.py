"""REUSE.toml gives every file an SPDX licence and copyright line.

Input: ``REUSE.toml``, ``LICENSES/``, ``LICENSE`` and ``THIRD_PARTY_NOTICES.md``.
Output: assertions. The full ``reuse lint`` runs in packaging.yml (reuse
6.2.0); this test checks the parts that drift without it: the files the
overrides name, the licence texts the expressions need, and the holder line
against LICENSE.
"""

from __future__ import annotations

import re
from pathlib import Path

import tomllib

REPO = Path(__file__).resolve().parents[1]
CONFIG = tomllib.loads((REPO / "REUSE.toml").read_text(encoding="utf-8"))
ANNOTATIONS = CONFIG["annotations"]


def _paths(annotation: dict) -> list[str]:
    path = annotation["path"]
    return [path] if isinstance(path, str) else list(path)


def test_the_first_annotation_covers_every_file_with_the_project_licence():
    first = ANNOTATIONS[0]
    assert first["path"] == "**"
    assert first["SPDX-License-Identifier"] == "LGPL-3.0-or-later"
    assert first["precedence"] == "closest"


def test_the_holder_line_is_the_one_license_states():
    stated = re.search(r"Copyright \(c\) (.+)", (REPO / "LICENSE").read_text())
    assert stated, "LICENSE has no copyright line"
    assert stated.group(1).strip() in ANNOTATIONS[0]["SPDX-FileCopyrightText"]


def test_every_override_names_files_that_exist_and_are_in_the_notices():
    notices = (REPO / "THIRD_PARTY_NOTICES.md").read_text(encoding="utf-8")
    for annotation in ANNOTATIONS[1:]:
        assert annotation["precedence"] == "override"
        for pattern in _paths(annotation):
            matches = list(REPO.glob(pattern))
            assert matches, f"REUSE.toml names {pattern}, which matches nothing"
            for match in matches:
                rel = match.relative_to(REPO).as_posix()
                assert rel.split("/")[-1] in notices, f"{rel} has no notice"


def test_every_licence_in_an_expression_has_its_text():
    used = set()
    for annotation in ANNOTATIONS:
        expression = annotation["SPDX-License-Identifier"]
        used |= {t for t in re.split(r"\s+|\(|\)", expression)
                 if t and t not in {"AND", "OR", "WITH"}}
    present = {p.stem for p in (REPO / "LICENSES").glob("*.txt")}
    assert used == present, f"used {sorted(used)}, texts {sorted(present)}"


def test_the_lgpl_text_in_licenses_is_the_shipped_one():
    assert ((REPO / "LICENSES" / "LGPL-3.0-or-later.txt").read_bytes()
            == (REPO / "COPYING.LESSER").read_bytes())
